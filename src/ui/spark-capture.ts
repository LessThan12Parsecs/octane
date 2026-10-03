/**
 * Captures the ignition event (dwell → breakdown → arc → glow → done) from the
 * snapshot stream as a time trace for the µs–ms spark chart. Keeps the last
 * complete event on display until the next one reaches breakdown, then shows the
 * new one live as it unfolds.
 */
import type { EngineSnapshot, SparkPhase } from '../physics/core/snapshot';
import { Float64Column } from './ring-buffer';
import { ampsToMA, voltsToKV } from './units';

export const SPARK_PHASE_CODE: Readonly<Record<SparkPhase, number>> = {
  off: 0,
  charging: 1,
  breakdown: 2,
  arc: 3,
  glow: 4,
  done: 5,
};

export const SPARK_PHASE_BY_CODE: readonly SparkPhase[] = ['off', 'charging', 'breakdown', 'arc', 'glow', 'done'];

/** True for phases in which the gap conducts (breakdown, arc, glow). */
export const isDischargePhase = (p: SparkPhase): boolean => p === 'breakdown' || p === 'arc' || p === 'glow';

export class SparkEvent {
  cycle = -1;
  /** Dwell (charging) start, s. */
  tStart = NaN;
  /** Breakdown time (first discharge-phase sample), s; NaN until it happens. */
  tBreakdown = NaN;
  /** Time the discharge ended ('done'/'off' after breakdown), s. */
  tEnd = NaN;
  complete = false;
  /** Gas-side electrical energy delivered (J) at the last sample. */
  energy = 0;
  /** Crank angle at breakdown, deg. */
  breakdownDeg = NaN;
  readonly t = new Float64Column(1024);
  /** kV */
  readonly vSec = new Float64Column(1024);
  /** mA */
  readonly iSec = new Float64Column(1024);
  /** A */
  readonly iPrim = new Float64Column(1024);
  /** kV */
  readonly vBd = new Float64Column(1024);
  readonly phase = new Float64Column(1024);

  get length(): number {
    return this.t.length;
  }

  reset(cycle: number, tStart: number): void {
    this.cycle = cycle;
    this.tStart = tStart;
    this.tBreakdown = NaN;
    this.tEnd = NaN;
    this.complete = false;
    this.energy = 0;
    this.breakdownDeg = NaN;
    this.t.clear();
    this.vSec.clear();
    this.iSec.clear();
    this.iPrim.clear();
    this.vBd.clear();
    this.phase.clear();
  }

  push(s: EngineSnapshot): void {
    const sp = s.spark;
    this.t.push(s.t);
    this.vSec.push(voltsToKV(Math.abs(sp.secondaryVoltage)));
    this.iSec.push(ampsToMA(sp.secondaryCurrent));
    this.iPrim.push(sp.primaryCurrent);
    this.vBd.push(voltsToKV(sp.breakdownVoltage));
    this.phase.push(SPARK_PHASE_CODE[sp.phase] ?? 0);
    this.energy = sp.energyDelivered;
  }

  /** Reference time for the chart x axis (breakdown if known, else dwell start). */
  get tRef(): number {
    return Number.isFinite(this.tBreakdown) ? this.tBreakdown : this.tStart;
  }
}

export interface SparkPhaseDurations {
  /** s */
  dwell: number;
  breakdown: number;
  arc: number;
  glow: number;
  /** Total discharge duration (breakdown + arc + glow), s. */
  discharge: number;
}

/** Time spent in each phase, from sample timestamps (each sample owns the interval to the next). */
export function sparkPhaseDurations(ev: SparkEvent): SparkPhaseDurations {
  const d: SparkPhaseDurations = { dwell: 0, breakdown: 0, arc: 0, glow: 0, discharge: 0 };
  const t = ev.t.data;
  const ph = ev.phase.data;
  for (let i = 0; i + 1 < ev.length; i++) {
    const dt = t[i + 1] - t[i];
    switch (ph[i]) {
      case SPARK_PHASE_CODE.charging:
        d.dwell += dt;
        break;
      case SPARK_PHASE_CODE.breakdown:
        d.breakdown += dt;
        break;
      case SPARK_PHASE_CODE.arc:
        d.arc += dt;
        break;
      case SPARK_PHASE_CODE.glow:
        d.glow += dt;
        break;
    }
  }
  d.discharge = d.breakdown + d.arc + d.glow;
  return d;
}

/** Peak of a column, NaN if empty. */
export function columnMax(c: Float64Column): number {
  let m = -Infinity;
  for (let i = 0; i < c.length; i++) if (c.data[i] > m) m = c.data[i];
  return c.length ? m : NaN;
}

export interface SparkCaptureOptions {
  /** Keep recording this long after the discharge ends, s. */
  tail: number;
  /** Close an event this long after breakdown even if it never reports 'done', s. */
  maxAfterBreakdown: number;
}

const DEFAULTS: SparkCaptureOptions = { tail: 0.5e-3, maxAfterBreakdown: 8e-3 };

export class SparkCapture {
  private building = new SparkEvent();
  private shown = new SparkEvent();
  private active = false;
  private hasShown = false;
  private lastT = -Infinity;
  private lastEventCycle = NaN;
  private readonly opts: SparkCaptureOptions;
  /** Incremented when the displayed event changes. */
  version = 0;

  constructor(opts: Partial<SparkCaptureOptions> = {}) {
    this.opts = { ...DEFAULTS, ...opts };
  }

  clear(): void {
    this.active = false;
    this.hasShown = false;
    this.lastT = -Infinity;
    this.lastEventCycle = NaN;
    this.building.reset(-1, NaN);
    this.shown.reset(-1, NaN);
    this.version++;
  }

  /** The event to display: the live one once it has broken down, else the last complete one. */
  get display(): SparkEvent | null {
    if (this.active && Number.isFinite(this.building.tBreakdown)) return this.building;
    return this.hasShown ? this.shown : null;
  }

  /** The event being recorded (may still be charging), or null. */
  get live(): SparkEvent | null {
    return this.active ? this.building : null;
  }

  ingest(snaps: readonly EngineSnapshot[], playheadT: number): boolean {
    let changed = false;
    for (let i = 0; i < snaps.length; i++) {
      const s = snaps[i];
      if (!(s.t > this.lastT) || s.t > playheadT + 1e-12) continue;
      this.lastT = s.t;
      if (this.step(s)) changed = true;
    }
    if (changed) this.version++;
    return changed;
  }

  /** Returns true if the displayed event changed. */
  private step(s: EngineSnapshot): boolean {
    const phase = s.spark.phase;
    let restarted = false;

    if (!this.active) {
      // A new event starts with a dwell; a discharge only opens one if this cycle has
      // none yet (stream joined mid-spark), not after a forced close of a long discharge.
      if (phase === 'charging' || (isDischargePhase(phase) && s.cycle !== this.lastEventCycle)) {
        this.building.reset(s.cycle, s.t);
        this.active = true;
        this.lastEventCycle = s.cycle;
      } else {
        return false;
      }
    } else if (phase === 'charging' && Number.isFinite(this.building.tBreakdown)) {
      // A new dwell started before the previous event was closed: close it and restart.
      this.finish();
      this.building.reset(s.cycle, s.t);
      this.active = true;
      this.lastEventCycle = s.cycle;
      restarted = true;
    }

    const ev2 = this.building;
    ev2.push(s);
    let changed = restarted || Number.isFinite(ev2.tBreakdown);

    if (!Number.isFinite(ev2.tBreakdown) && isDischargePhase(phase)) {
      ev2.tBreakdown = s.t;
      ev2.breakdownDeg = s.thetaDeg;
      changed = true;
    }
    if (Number.isFinite(ev2.tBreakdown)) {
      if (!Number.isFinite(ev2.tEnd) && (phase === 'done' || phase === 'off')) ev2.tEnd = s.t;
      const tailDone = Number.isFinite(ev2.tEnd) && s.t - ev2.tEnd >= this.opts.tail;
      const tooLong = s.t - ev2.tBreakdown >= this.opts.maxAfterBreakdown;
      if (tailDone || tooLong) {
        if (!Number.isFinite(ev2.tEnd)) ev2.tEnd = s.t;
        this.finish();
        changed = true;
      }
    } else if (phase === 'done' || phase === 'off') {
      // Dwell ended without breakdown (e.g. insufficient voltage): keep it as a (failed) event.
      ev2.tEnd = s.t;
      this.finish();
      changed = true;
    }
    return changed;
  }

  private finish(): void {
    this.building.complete = true;
    const tmp = this.shown;
    this.shown = this.building;
    this.building = tmp;
    this.hasShown = true;
    this.active = false;
  }
}
