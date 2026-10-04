/**
 * Per-cycle history of the snapshot stream, stored column-wise in display units
 * (bar, cm³, K, J/°CA, mm, g/s) for the crank-angle charts.
 *
 * The UI receives the playback snapshot plus a `recent` array every frame. Whatever
 * window `recent` covers, the store ingests each snapshot exactly once (strictly
 * increasing t) and never beyond the playback cursor, so the charts always show
 * exactly what has "happened" at the current playback time.
 */
import type { EngineSnapshot } from '../physics/core/snapshot';
import { upperBound } from './ring-buffer';
import { cyclePeriod, kgPerSToGPerS, m3ToCm3, mToMm, paToBar, wattsToJPerDeg } from './units';

/** Column indices. */
export const CH = {
  t: 0,
  theta: 1,
  /** bar */
  p: 2,
  /** cm³ */
  V: 3,
  /** K (NaN when the zone does not exist) */
  Tu: 4,
  Tb: 5,
  Tm: 6,
  /** burned mass fraction, – */
  xb: 7,
  /** J/°CA */
  hrr: 8,
  /** Livengood–Wu knock integral, – */
  lw: 9,
  /** mm */
  liftIn: 10,
  liftEx: 11,
  /** g/s (intake + into cylinder, exhaust + out of cylinder) */
  mIn: 12,
  mEx: 13,
  /** synthesised knock oscillation, bar */
  knock: 14,
  rpm: 15,
} as const;
export type ChannelKey = keyof typeof CH;
export const N_CH = 16;

/** Convert one snapshot into the display-unit sample vector `out` (length N_CH). */
export function sampleChannels(s: EngineSnapshot, out: Float64Array): Float64Array {
  out[CH.t] = s.t;
  out[CH.theta] = s.thetaDeg;
  out[CH.p] = paToBar(s.pressure);
  out[CH.V] = m3ToCm3(s.volume);
  out[CH.Tu] = s.temperatureUnburned > 0 ? s.temperatureUnburned : NaN;
  out[CH.Tb] = s.temperatureBurned > 0 ? s.temperatureBurned : NaN;
  out[CH.Tm] = s.temperatureMean > 0 ? s.temperatureMean : NaN;
  out[CH.xb] = s.massFractionBurned;
  out[CH.hrr] = wattsToJPerDeg(s.heatReleaseRate, s.rpm);
  out[CH.lw] = s.knock.integral;
  out[CH.liftIn] = mToMm(s.intakeLift);
  out[CH.liftEx] = mToMm(s.exhaustLift);
  out[CH.mIn] = kgPerSToGPerS(s.intakeMassFlow);
  out[CH.mEx] = kgPerSToGPerS(s.exhaustMassFlow);
  out[CH.knock] = paToBar(s.knock.oscillation);
  out[CH.rpm] = s.rpm;
  return out;
}

/**
 * Maps a snapshot to one sample row of a store. Channels 0 and 1 must be t and θ (CH.t, CH.theta): the
 * store splits traces on them and the θ plots use column 1 as x.
 */
export interface TraceSampler {
  readonly channels: number;
  sample(s: EngineSnapshot, out: Float64Array): Float64Array;
}

/** The default sampler: the 16 single-cylinder channels of CH. */
export const CYLINDER_SAMPLER: TraceSampler = { channels: N_CH, sample: sampleChannels };

/** Column of cylinder i's pressure in a cylinderPressureSampler store. */
export const overlayPressureChannel = (i: number): number => 2 + i;

/**
 * Sampler for the "all cylinders" pressure overlay: t, ENGINE θ (cylinder 1's angle), then every
 * cylinder's pressure in bar (NaN when the snapshot has no such cylinder).
 */
export function cylinderPressureSampler(cylinders: number): TraceSampler {
  return {
    channels: 2 + cylinders,
    sample(s: EngineSnapshot, out: Float64Array): Float64Array {
      out[CH.t] = s.t;
      out[CH.theta] = s.thetaDeg;
      const cyl = s.cylinders;
      for (let i = 0; i < cylinders; i++) {
        const c = cyl ? cyl[i] : i === 0 ? s : undefined;
        out[2 + i] = c ? paToBar(c.pressure) : NaN;
      }
      return out;
    },
  };
}

/** One engine cycle's samples, column-major, growable without per-sample allocation. */
export class CycleTrace {
  cycle = -1;
  n = 0;
  cols: Float64Array[];
  /** End-gas autoignition happened during this cycle. */
  knocked = false;
  /** Crank angle of the first sample flagged autoignited (NaN if none). */
  knockDeg = NaN;

  constructor(
    capacity = 2048,
    readonly channels = N_CH,
  ) {
    this.cols = [];
    for (let i = 0; i < channels; i++) this.cols.push(new Float64Array(Math.max(16, capacity)));
  }

  get capacity(): number {
    return this.cols[0].length;
  }

  reset(cycle: number): void {
    this.cycle = cycle;
    this.n = 0;
    this.knocked = false;
    this.knockDeg = NaN;
  }

  push(sample: Float64Array): void {
    const nc = this.channels;
    if (this.n === this.capacity) {
      const cap = this.capacity * 2;
      for (let i = 0; i < nc; i++) {
        const next = new Float64Array(cap);
        next.set(this.cols[i]);
        this.cols[i] = next;
      }
    }
    const k = this.n++;
    for (let i = 0; i < nc; i++) this.cols[i][k] = sample[i];
  }

  get tStart(): number {
    return this.n > 0 ? this.cols[CH.t][0] : NaN;
  }

  get tEnd(): number {
    return this.n > 0 ? this.cols[CH.t][this.n - 1] : NaN;
  }

  get lastTheta(): number {
    return this.n > 0 ? this.cols[CH.theta][this.n - 1] : NaN;
  }

  /** Column view of the first `n` samples (shares memory). */
  col(ch: number, n: number = this.n): Float64Array {
    return this.cols[ch].subarray(0, Math.min(n, this.n));
  }
}

export interface TraceSelection {
  /** Trace containing the playback cursor (null before any data). */
  current: CycleTrace | null;
  /** Number of samples of `current` at or before the playback time. */
  currentCount: number;
  /** Older cycles, newest first. */
  ghosts: CycleTrace[];
}

export function emptySelection(): TraceSelection {
  return { current: null, currentCount: 0, ghosts: [] };
}

/** Tolerance on simulated-time comparisons, s. */
const T_EPS = 1e-12;

export class CycleTraceStore {
  /** Oldest → newest. */
  readonly traces: CycleTrace[] = [];
  private readonly pool: CycleTrace[] = [];
  private readonly sample: Float64Array;
  /** Time of the last ingested snapshot. */
  lastT = -Infinity;
  /** Incremented whenever the stored data changes (cheap dirty check for charts). */
  version = 0;

  constructor(
    readonly maxCycles = 6,
    private readonly initialCapacity = 2048,
    private readonly sampler: TraceSampler = CYLINDER_SAMPLER,
  ) {
    this.sample = new Float64Array(sampler.channels);
  }

  clear(): void {
    while (this.traces.length) this.pool.push(this.traces.pop()!);
    this.lastT = -Infinity;
    this.version++;
  }

  get newest(): CycleTrace | null {
    return this.traces.length ? this.traces[this.traces.length - 1] : null;
  }

  get earliestT(): number {
    return this.traces.length ? this.traces[0].tStart : NaN;
  }

  /**
   * Clear the history if the playback cursor jumped backwards past what we hold
   * (simulator reset, or a rewind longer than two cycles). Returns true if cleared.
   */
  checkDiscontinuity(playheadT: number, playheadCycle: number, rpm: number): boolean {
    if (!this.traces.length) {
      return false;
    }
    const oldest = this.traces[0];
    const period = cyclePeriod(rpm);
    const rewind = this.lastT - playheadT;
    if (
      playheadT + T_EPS < oldest.tStart ||
      playheadCycle < oldest.cycle ||
      (Number.isFinite(period) && rewind > 2 * period)
    ) {
      this.clear();
      return true;
    }
    return false;
  }

  /**
   * Ingest snapshots with lastT < t ≤ playheadT. Returns the number of samples added.
   * Snapshots are expected in time order; out-of-order or duplicate ones are skipped.
   */
  ingest(snaps: readonly EngineSnapshot[], playheadT: number): number {
    let added = 0;
    for (let i = 0; i < snaps.length; i++) {
      const s = snaps[i];
      if (!(s.t > this.lastT + T_EPS) || s.t > playheadT + T_EPS) continue;
      this.lastT = s.t;
      let cur = this.newest;
      if (!cur || s.cycle !== cur.cycle || s.thetaDeg < cur.lastTheta - 180) {
        cur = this.startTrace(s.cycle);
      } else if (s.thetaDeg < cur.lastTheta) {
        continue; // keep θ monotonic within a trace (x must be sorted for the charts)
      }
      cur.push(this.sampler.sample(s, this.sample));
      if (s.knock.autoignited && !cur.knocked) {
        cur.knocked = true;
        cur.knockDeg = s.thetaDeg;
      }
      added++;
    }
    if (added) this.version++;
    return added;
  }

  private startTrace(cycle: number): CycleTrace {
    let tr: CycleTrace;
    if (this.traces.length >= this.maxCycles) tr = this.traces.shift()!;
    else tr = this.pool.pop() ?? new CycleTrace(this.initialCapacity, this.sampler.channels);
    tr.reset(cycle);
    this.traces.push(tr);
    return tr;
  }

  /** The trace under the playback cursor, its visible sample count, and older cycles. */
  select(playheadT: number, maxGhosts: number, out: TraceSelection = emptySelection()): TraceSelection {
    out.current = null;
    out.currentCount = 0;
    out.ghosts.length = 0;
    let ci = -1;
    for (let i = this.traces.length - 1; i >= 0; i--) {
      if (this.traces[i].n > 0 && this.traces[i].tStart <= playheadT + T_EPS) {
        ci = i;
        break;
      }
    }
    if (ci < 0) return out;
    const cur = this.traces[ci];
    out.current = cur;
    out.currentCount = upperBound(cur.cols[CH.t], cur.n, playheadT + T_EPS);
    for (let i = ci - 1; i >= 0 && out.ghosts.length < maxGhosts; i--) out.ghosts.push(this.traces[i]);
    return out;
  }
}

/**
 * Copy `src[0..n)` × `scale` into the reusable plain array `out`, mapping NaN/±Inf
 * to null (uPlot draws gaps at nulls; typed arrays cannot hold them).
 */
export function fillSeries(src: ArrayLike<number>, n: number, out: (number | null)[], scale = 1): (number | null)[] {
  out.length = n;
  for (let i = 0; i < n; i++) {
    const v = src[i] * scale;
    out[i] = Number.isFinite(v) ? v : null;
  }
  return out;
}

/**
 * Polytropic index n from a least-squares fit of ln p = c − n ln V over the samples
 * with θ ∈ [theta1, theta2]. NaN if fewer than 3 usable samples.
 */
export function polytropicIndex(trace: CycleTrace, count: number, theta1: number, theta2: number): number {
  const th = trace.cols[CH.theta];
  const p = trace.cols[CH.p];
  const v = trace.cols[CH.V];
  const n = Math.min(count, trace.n);
  let sx = 0;
  let sy = 0;
  let sxx = 0;
  let sxy = 0;
  let m = 0;
  for (let i = 0; i < n; i++) {
    const a = th[i];
    if (a < theta1 || a > theta2) continue;
    if (!(p[i] > 0) || !(v[i] > 0)) continue;
    const x = Math.log(v[i]);
    const y = Math.log(p[i]);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
    m++;
  }
  if (m < 3) return NaN;
  const den = m * sxx - sx * sx;
  if (Math.abs(den) < 1e-300) return NaN;
  return -(m * sxy - sx * sy) / den;
}
