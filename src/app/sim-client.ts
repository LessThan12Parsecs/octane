/**
 * Main-thread side of the simulation pipeline.
 *
 * `SimClient` spawns the simulation worker (src/worker/sim.worker.ts), buffers
 * the snapshots it streams in a bounded ring buffer, plays them back on a
 * simulated-time clock (usually in slow motion) and interpolates between the
 * two snapshots bracketing the playback time. It keeps the worker's look-ahead
 * small in wall-clock terms so operating-point changes show up quickly even in
 * slow motion (see `maybeDemand`).
 *
 * The buffer / interpolation logic is pure (`SnapshotBuffer`,
 * `interpolateSnapshot`) and unit-tested without a Worker; `SimClient` accepts
 * an injectable worker factory for the same reason.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import type { FromWorker, SimulatorOptions, ToWorker } from '../worker/protocol';

// ---------------------------------------------------------------------------
// Snapshot helpers
// ---------------------------------------------------------------------------

/** A zero-filled snapshot (used as a reusable interpolation target). */
export function createEmptySnapshot(): EngineSnapshot {
  return {
    t: 0,
    cycle: 0,
    thetaDeg: -360,
    rpm: 0,
    pistonDisplacement: 0,
    clearanceHeight: 0,
    rodAngle: 0,
    intakeLift: 0,
    exhaustLift: 0,
    phase: 'gas-exchange',
    volume: 0,
    pressure: 0,
    temperatureMean: 0,
    temperatureUnburned: 0,
    temperatureBurned: 0,
    massFractionBurned: 0,
    mass: 0,
    heatReleaseRate: 0,
    heatLossRate: 0,
    flame: {
      stage: 'none',
      radius: 0,
      center: [0, 0, 0],
      area: 0,
      laminarSpeed: 0,
      turbulentSpeed: 0,
      turbulenceIntensity: 0,
    },
    spark: {
      phase: 'off',
      primaryCurrent: 0,
      secondaryVoltage: 0,
      secondaryCurrent: 0,
      energyDelivered: 0,
      breakdownVoltage: 0,
    },
    intakeMassFlow: 0,
    exhaustMassFlow: 0,
    intakeManifoldPressure: 0,
    exhaustManifoldPressure: 0,
    knock: { integral: 0, autoignited: false, oscillation: 0 },
    burnedComposition: { CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0 },
    gasTorque: 0,
    netTorque: 0,
  };
}

/** Monotonic crank angle, deg: cycle·720 + θ + 360 (0 at t = 0 for a run starting at θ = −360). */
export function cumulativeCrankDeg(s: Pick<EngineSnapshot, 'cycle' | 'thetaDeg'>): number {
  return s.cycle * 720 + s.thetaDeg + 360;
}

/**
 * Interpolate a crank angle in the [−360, 360) convention, taking the short way
 * forward across the ±360 wrap. Returns the angle wrapped back into [−360, 360).
 */
export function lerpCrankAngleDeg(a: number, b: number, alpha: number): number {
  let d = b - a;
  if (d < -360) d += 720;
  else if (d > 360) d -= 720;
  let th = a + d * alpha;
  if (th >= 360) th -= 720;
  else if (th < -360) th += 720;
  return th;
}

const lerp = (x: number, y: number, u: number): number => x + (y - x) * u;

/**
 * Interpolate between two snapshots (a earlier, b later) into `out`:
 * linear for scalars, wrap-aware for thetaDeg; discrete fields (phase, flame
 * stage, spark phase, autoignited) and `cycle` come from the earlier sample —
 * or from the later one when the interpolated angle has crossed the cycle wrap,
 * so the result is always consistent with its own (cycle, thetaDeg).
 */
export function interpolateSnapshot(
  a: EngineSnapshot,
  b: EngineSnapshot,
  alpha: number,
  out: EngineSnapshot,
): EngineSnapshot {
  const u = alpha < 0 ? 0 : alpha > 1 ? 1 : alpha;
  const th = lerpCrankAngleDeg(a.thetaDeg, b.thetaDeg, u);
  const wrapped = b.cycle !== a.cycle && th < a.thetaDeg;
  const d = wrapped ? b : a;
  out.t = lerp(a.t, b.t, u);
  out.cycle = d.cycle;
  out.thetaDeg = th;
  out.rpm = lerp(a.rpm, b.rpm, u);
  out.pistonDisplacement = lerp(a.pistonDisplacement, b.pistonDisplacement, u);
  out.clearanceHeight = lerp(a.clearanceHeight, b.clearanceHeight, u);
  out.rodAngle = lerp(a.rodAngle, b.rodAngle, u);
  out.intakeLift = lerp(a.intakeLift, b.intakeLift, u);
  out.exhaustLift = lerp(a.exhaustLift, b.exhaustLift, u);
  out.phase = d.phase;
  out.volume = lerp(a.volume, b.volume, u);
  out.pressure = lerp(a.pressure, b.pressure, u);
  out.temperatureMean = lerp(a.temperatureMean, b.temperatureMean, u);
  out.temperatureUnburned = lerp(a.temperatureUnburned, b.temperatureUnburned, u);
  out.temperatureBurned = lerp(a.temperatureBurned, b.temperatureBurned, u);
  out.massFractionBurned = lerp(a.massFractionBurned, b.massFractionBurned, u);
  out.mass = lerp(a.mass, b.mass, u);
  out.heatReleaseRate = lerp(a.heatReleaseRate, b.heatReleaseRate, u);
  out.heatLossRate = lerp(a.heatLossRate, b.heatLossRate, u);

  const fa = a.flame;
  const fb = b.flame;
  const fo = out.flame;
  fo.stage = d.flame.stage;
  fo.radius = lerp(fa.radius, fb.radius, u);
  fo.center[0] = lerp(fa.center[0], fb.center[0], u);
  fo.center[1] = lerp(fa.center[1], fb.center[1], u);
  fo.center[2] = lerp(fa.center[2], fb.center[2], u);
  fo.area = lerp(fa.area, fb.area, u);
  fo.laminarSpeed = lerp(fa.laminarSpeed, fb.laminarSpeed, u);
  fo.turbulentSpeed = lerp(fa.turbulentSpeed, fb.turbulentSpeed, u);
  fo.turbulenceIntensity = lerp(fa.turbulenceIntensity, fb.turbulenceIntensity, u);

  const sa = a.spark;
  const sb = b.spark;
  const so = out.spark;
  so.phase = d.spark.phase;
  so.primaryCurrent = lerp(sa.primaryCurrent, sb.primaryCurrent, u);
  so.secondaryVoltage = lerp(sa.secondaryVoltage, sb.secondaryVoltage, u);
  so.secondaryCurrent = lerp(sa.secondaryCurrent, sb.secondaryCurrent, u);
  so.energyDelivered = lerp(sa.energyDelivered, sb.energyDelivered, u);
  so.breakdownVoltage = lerp(sa.breakdownVoltage, sb.breakdownVoltage, u);

  out.intakeMassFlow = lerp(a.intakeMassFlow, b.intakeMassFlow, u);
  out.exhaustMassFlow = lerp(a.exhaustMassFlow, b.exhaustMassFlow, u);
  out.intakeManifoldPressure = lerp(a.intakeManifoldPressure, b.intakeManifoldPressure, u);
  out.exhaustManifoldPressure = lerp(a.exhaustManifoldPressure, b.exhaustManifoldPressure, u);

  out.knock.integral = lerp(a.knock.integral, b.knock.integral, u);
  out.knock.autoignited = d.knock.autoignited;
  out.knock.oscillation = lerp(a.knock.oscillation, b.knock.oscillation, u);

  const ca = a.burnedComposition;
  const cb = b.burnedComposition;
  const co = out.burnedComposition;
  co.CO2 = lerp(ca.CO2, cb.CO2, u);
  co.H2O = lerp(ca.H2O, cb.H2O, u);
  co.CO = lerp(ca.CO, cb.CO, u);
  co.O2 = lerp(ca.O2, cb.O2, u);
  co.H2 = lerp(ca.H2, cb.H2, u);
  co.OH = lerp(ca.OH, cb.OH, u);
  co.H = lerp(ca.H, cb.H, u);
  co.O = lerp(ca.O, cb.O, u);
  co.NO = lerp(ca.NO, cb.NO, u);
  co.N2 = lerp(ca.N2, cb.N2, u);

  out.gasTorque = lerp(a.gasTorque, b.gasTorque, u);
  out.netTorque = lerp(a.netTorque, b.netTorque, u);
  return out;
}

/** Deep-copy `src` into `out` (no allocation). */
export function copySnapshot(src: EngineSnapshot, out: EngineSnapshot): EngineSnapshot {
  return interpolateSnapshot(src, src, 0, out);
}

// ---------------------------------------------------------------------------
// Ring buffer
// ---------------------------------------------------------------------------

/**
 * Time-ordered ring buffer of snapshots (strictly increasing `t`), growable
 * (power-of-two capacity) and trimmed from the front by the owner.
 */
export class SnapshotBuffer {
  private buf: (EngineSnapshot | undefined)[];
  private mask: number;
  private head = 0;
  private len = 0;

  constructor(initialCapacity = 4096) {
    let cap = 16;
    while (cap < initialCapacity) cap *= 2;
    this.buf = new Array<EngineSnapshot | undefined>(cap);
    this.mask = cap - 1;
  }

  get length(): number {
    return this.len;
  }
  get capacity(): number {
    return this.buf.length;
  }
  /** i-th sample from the oldest (0) to the newest (length − 1). */
  at(i: number): EngineSnapshot {
    return this.buf[(this.head + i) & this.mask]!;
  }
  get first(): EngineSnapshot | undefined {
    return this.len > 0 ? this.at(0) : undefined;
  }
  get last(): EngineSnapshot | undefined {
    return this.len > 0 ? this.at(this.len - 1) : undefined;
  }
  /** Time of the oldest / newest sample (NaN when empty). */
  get earliestTime(): number {
    return this.len > 0 ? this.at(0).t : NaN;
  }
  get latestTime(): number {
    return this.len > 0 ? this.at(this.len - 1).t : NaN;
  }

  clear(): void {
    this.buf.fill(undefined);
    this.head = 0;
    this.len = 0;
  }

  /** Append; samples not strictly later than the newest one are dropped (returns false). */
  push(s: EngineSnapshot): boolean {
    if (this.len > 0 && !(s.t > this.at(this.len - 1).t)) return false;
    if (this.len === this.buf.length) this.grow();
    this.buf[(this.head + this.len) & this.mask] = s;
    this.len++;
    return true;
  }

  pushBatch(batch: readonly EngineSnapshot[]): number {
    let n = 0;
    for (let i = 0; i < batch.length; i++) if (this.push(batch[i])) n++;
    return n;
  }

  private grow(): void {
    const old = this.buf;
    const cap = old.length * 2;
    const next = new Array<EngineSnapshot | undefined>(cap);
    for (let i = 0; i < this.len; i++) next[i] = old[(this.head + i) & this.mask];
    this.buf = next;
    this.mask = cap - 1;
    this.head = 0;
  }

  /** Index of the last sample with t ≤ time (−1 if none). */
  indexAtOrBefore(time: number): number {
    let lo = 0;
    let hi = this.len - 1;
    if (hi < 0 || time < this.at(0).t) return -1;
    if (time >= this.at(hi).t) return hi;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (this.at(mid).t <= time) lo = mid;
      else hi = mid - 1;
    }
    return lo;
  }

  /** Index of the last sample with cumulative crank angle ≤ phi (−1 if none). */
  indexAtOrBeforeAngle(phi: number): number {
    let lo = 0;
    let hi = this.len - 1;
    if (hi < 0 || phi < cumulativeCrankDeg(this.at(0))) return -1;
    if (phi >= cumulativeCrankDeg(this.at(hi))) return hi;
    while (lo < hi) {
      const mid = (lo + hi + 1) >> 1;
      if (cumulativeCrankDeg(this.at(mid)) <= phi) lo = mid;
      else hi = mid - 1;
    }
    return lo;
  }

  /** Interpolated state at `time` (clamped to the buffered range) into `out`; null if empty. */
  sample(time: number, out: EngineSnapshot): EngineSnapshot | null {
    if (this.len === 0) return null;
    const i = this.indexAtOrBefore(time);
    if (i < 0) return copySnapshot(this.at(0), out);
    if (i >= this.len - 1) return copySnapshot(this.at(this.len - 1), out);
    const a = this.at(i);
    const b = this.at(i + 1);
    interpolateSnapshot(a, b, (time - a.t) / (b.t - a.t), out);
    out.t = time;
    return out;
  }

  /** Simulated time at which the cumulative crank angle equals `phi` (clamped to the buffer). */
  timeAtCumulativeAngle(phi: number): number {
    if (this.len === 0) return NaN;
    const i = this.indexAtOrBeforeAngle(phi);
    if (i < 0) return this.at(0).t;
    if (i >= this.len - 1) return this.at(this.len - 1).t;
    const a = this.at(i);
    const b = this.at(i + 1);
    const pa = cumulativeCrankDeg(a);
    const pb = cumulativeCrankDeg(b);
    return pb > pa ? a.t + ((phi - pa) / (pb - pa)) * (b.t - a.t) : a.t;
  }

  /**
   * Drop samples older than `time`, but always keep the last sample at or
   * before `time` (the interpolation bracket). Returns the number dropped.
   */
  trimBefore(time: number): number {
    const i = this.indexAtOrBefore(time);
    if (i <= 0) return 0;
    for (let k = 0; k < i; k++) this.buf[(this.head + k) & this.mask] = undefined;
    this.head = (this.head + i) & this.mask;
    this.len -= i;
    return i;
  }

  /** Samples with t0 ≤ t ≤ t1 (oldest first) into `out` (cleared first). */
  collect(t0: number, t1: number, out: EngineSnapshot[]): EngineSnapshot[] {
    out.length = 0;
    if (this.len === 0 || t1 < t0) return out;
    let i = this.indexAtOrBefore(t0);
    if (i < 0) i = 0;
    else if (this.at(i).t < t0) i++;
    for (; i < this.len; i++) {
      const s = this.at(i);
      if (s.t > t1) break;
      out.push(s);
    }
    return out;
  }
}

// ---------------------------------------------------------------------------
// Client
// ---------------------------------------------------------------------------

/** The subset of `Worker` the client uses (injectable for tests). */
export interface WorkerLike {
  postMessage(msg: ToWorker): void;
  terminate(): void;
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null;
  onerror: ((ev: ErrorEvent) => void) | null;
}

export interface SimClientConfig {
  /** Worker factory (default: the Vite-bundled module worker). */
  createWorker?: () => WorkerLike;
  /** Cycles of history kept behind the playback cursor for charts (default 4). */
  historyCycles?: number;
  /** Look-ahead requested from the worker in wall-clock seconds at the current time scale (default 0.5). */
  aheadWallSeconds?: number;
  /** Minimum look-ahead, simulated seconds (default 0.02 — lets paused single-stepping proceed). */
  minAheadSeconds?: number;
}

function defaultCreateWorker(): WorkerLike {
  return new Worker(new URL('../worker/sim.worker.ts', import.meta.url), { type: 'module' });
}

/** Max cycle summaries held back waiting for playback to reach them. */
const MAX_PENDING_CYCLES = 256;

export class SimClient {
  /** Resolves when the worker has initialised the simulator. */
  readonly ready: Promise<void>;

  private readonly worker: WorkerLike;
  private readonly options: SimulatorOptions;
  private readonly historyCycles: number;
  private readonly aheadWall: number;
  private readonly minAhead: number;
  private readonly buffer = new SnapshotBuffer();
  private readonly cycleListeners: ((c: CycleSummary) => void)[] = [];
  private pendingCycles: CycleSummary[] = [];
  private lastEmittedCycle = -1;
  private playback = NaN;
  /** Data is ignored until the worker's next 'ready' (after init / reset). */
  private awaitingReady = true;
  private lastDemand = NaN;
  private timeScale = 0.1;
  private readonly outs: [EngineSnapshot, EngineSnapshot] = [createEmptySnapshot(), createEmptySnapshot()];
  private outIndex = 0;
  private readonly scratch: EngineSnapshot = createEmptySnapshot();
  private readonly recent: EngineSnapshot[] = [];
  private resolveReady!: () => void;
  private rejectReady!: (e: unknown) => void;
  private isReady = false;
  private disposed = false;

  constructor(spec: EngineSpec, op: OperatingPoint, options: SimulatorOptions, config: SimClientConfig = {}) {
    this.options = { ...options };
    this.historyCycles = Math.max(1, config.historyCycles ?? 4);
    this.aheadWall = Math.max(0, config.aheadWallSeconds ?? 0.5);
    this.minAhead = Math.max(0, config.minAheadSeconds ?? 0.02);
    this.ready = new Promise<void>((resolve, reject) => {
      this.resolveReady = resolve;
      this.rejectReady = reject;
    });
    // Avoid unhandled-rejection noise if nobody awaits `ready`.
    this.ready.catch(() => undefined);
    this.worker = (config.createWorker ?? defaultCreateWorker)();
    this.worker.onmessage = (ev) => this.handleMessage(ev.data);
    this.worker.onerror = (ev) => {
      console.error('[sim worker] error', ev.message ?? ev);
      if (!this.isReady) this.rejectReady(new Error(`simulation worker failed: ${ev.message}`));
    };
    this.post({ type: 'init', spec, operatingPoint: op, options: this.options });
    this.maybeDemand(true);
  }

  /** Current playback position, simulated seconds (NaN until data arrives). */
  get playbackTime(): number {
    return this.playback;
  }
  /** Newest buffered simulated time (NaN when empty). */
  get bufferedUntil(): number {
    return this.buffer.latestTime;
  }
  /** Number of snapshots currently held. */
  get bufferedCount(): number {
    return this.buffer.length;
  }

  onCycle(cb: (c: CycleSummary) => void): void {
    this.cycleListeners.push(cb);
  }

  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.post({ type: 'set-operating-point', patch });
  }

  /** Restart the simulation at t = 0; buffered data and pending summaries are discarded. */
  reset(): void {
    this.post({ type: 'reset' });
    this.buffer.clear();
    this.pendingCycles = [];
    this.lastEmittedCycle = -1;
    this.playback = NaN;
    this.awaitingReady = true;
    this.lastDemand = NaN;
    this.maybeDemand(true);
  }

  /**
   * Move the playback clock by dtWall × timeScale simulated seconds and return
   * the interpolated snapshot (null until data is available). Never runs ahead
   * of buffered data: if the worker lags, playback stalls at the newest sample.
   * The returned object is owned by the client and stays valid until the
   * call after next (two alternating scratch objects) — copy it to keep it.
   */
  advance(dtWall: number, timeScale: number): EngineSnapshot | null {
    if (this.disposed) return null;
    this.timeScale = Math.max(0, timeScale);
    if (this.buffer.length === 0) {
      this.maybeDemand(false);
      return null;
    }
    if (Number.isNaN(this.playback)) this.playback = this.buffer.earliestTime;
    const dt = Math.max(0, dtWall) * this.timeScale;
    this.playback = Math.min(this.playback + dt, this.buffer.latestTime);
    return this.finishMove();
  }

  /**
   * Single-step the (paused) playback by `deg` crank degrees (negative steps
   * back through the retained history). Clamped to the buffered range.
   */
  stepDegrees(deg: number): EngineSnapshot | null {
    if (this.disposed || this.buffer.length === 0) {
      this.maybeDemand(false);
      return null;
    }
    if (Number.isNaN(this.playback)) this.playback = this.buffer.earliestTime;
    const cur = this.buffer.sample(this.playback, this.scratch)!;
    this.playback = this.buffer.timeAtCumulativeAngle(cumulativeCrankDeg(cur) + deg);
    return this.finishMove();
  }

  /**
   * Buffered raw snapshots with sinceT ≤ t ≤ playbackTime, oldest first (for
   * charts). The array is reused by the next call; don't keep it.
   */
  recentSnapshots(sinceT: number): EngineSnapshot[] {
    if (Number.isNaN(this.playback)) {
      this.recent.length = 0;
      return this.recent;
    }
    return this.buffer.collect(sinceT, this.playback, this.recent);
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.worker.onmessage = null;
    this.worker.onerror = null;
    this.worker.terminate();
    this.buffer.clear();
    this.pendingCycles = [];
    this.cycleListeners.length = 0;
  }

  /** Handle one message from the worker (public for tests / custom transports). */
  handleMessage(msg: FromWorker): void {
    if (this.disposed) return;
    switch (msg.type) {
      case 'ready':
        this.awaitingReady = false;
        if (!this.isReady) {
          this.isReady = true;
          this.resolveReady();
        }
        break;
      case 'snapshots':
        if (!this.awaitingReady) this.buffer.pushBatch(msg.batch);
        break;
      case 'cycle':
        if (!this.awaitingReady && msg.summary.cycle > this.lastEmittedCycle) {
          this.pendingCycles.push(msg.summary);
          if (this.pendingCycles.length > MAX_PENDING_CYCLES) this.pendingCycles.shift();
        }
        break;
      case 'log': {
        const fn = msg.level === 'error' ? console.error : msg.level === 'warn' ? console.warn : console.info;
        fn(`[sim worker] ${msg.message}`);
        break;
      }
    }
  }

  // -------------------------------------------------------------------------

  private post(msg: ToWorker): void {
    if (!this.disposed) this.worker.postMessage(msg);
  }

  private nextOut(): EngineSnapshot {
    this.outIndex ^= 1;
    return this.outs[this.outIndex];
  }

  private finishMove(): EngineSnapshot {
    const out = this.buffer.sample(this.playback, this.nextOut())!;
    this.flushCycles(out.cycle);
    const rpm = out.rpm > 1 ? out.rpm : 600;
    const keep = (this.historyCycles * 120) / rpm;
    this.buffer.trimBefore(this.playback - keep);
    this.maybeDemand(false);
    return out;
  }

  /** Deliver summaries of cycles that playback has completed. */
  private flushCycles(currentCycle: number): void {
    while (this.pendingCycles.length > 0 && this.pendingCycles[0].cycle < currentCycle) {
      const c = this.pendingCycles.shift()!;
      this.lastEmittedCycle = c.cycle;
      for (const cb of this.cycleListeners) cb(c);
    }
  }

  /**
   * Keep the worker buffering `ahead` simulated seconds past playback, where
   * `ahead` ≈ timeScale × aheadWallSeconds (bounded by [minAhead, bufferAheadSeconds]).
   * The worker buffers to `demand.t + bufferAheadSeconds`, so we send
   * `demand.t = playback + ahead − bufferAheadSeconds`.
   */
  private maybeDemand(force: boolean): void {
    const bufferAhead = Math.max(0, this.options.bufferAheadSeconds);
    const ahead = Math.min(bufferAhead, Math.max(this.minAhead, this.timeScale * this.aheadWall));
    const base = Number.isNaN(this.playback) ? 0 : this.playback;
    const t = base + ahead - bufferAhead;
    if (force || Number.isNaN(this.lastDemand) || Math.abs(t - this.lastDemand) > 0.1 * Math.max(ahead, 1e-4)) {
      this.lastDemand = t;
      this.post({ type: 'demand', t });
    }
  }
}
