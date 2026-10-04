/**
 * Captures the ignition event (dwell → breakdown → arc → glow → done) from the
 * snapshot stream as a time trace for the µs–ms spark chart. Keeps the last
 * complete event on display until the next one reaches breakdown, then shows the
 * new one live as it unfolds.
 *
 * Trembler (vibrator) coils fire a SHOWER of sparks for as long as the ignition
 * timer grounds the coil: charging → breakdown → arc/glow → charging → … The
 * snapshot then carries `spark.timerClosed` (and `breakdownCount`,
 * `firstSparkDeg`); one event is one timer contact, from the contact closing
 * until it opens and the last discharge has died away, however many sparks it
 * contains. Feed one SparkCapture per cylinder (cylinder views) in a
 * multi-cylinder engine.
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
  /**
   * Local crank angle of the FIRST breakdown of the event, deg: the simulator's `spark.firstSparkDeg`
   * when it reports one, else the angle of the first discharge sample (NaN before breakdown).
   */
  firstSparkDeg = NaN;
  /** Gap breakdowns in the event (a trembler shower has many; an inductive spark one). */
  sparkCount = 0;
  /** The event is a trembler shower (the stream reports the timer contact). */
  shower = false;
  private inDischarge = false;
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
    this.firstSparkDeg = NaN;
    this.sparkCount = 0;
    this.shower = false;
    this.inDischarge = false;
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
    const discharging = isDischargePhase(sp.phase);
    if (discharging && !this.inDischarge) this.sparkCount++;
    this.inDischarge = discharging;
    if (sp.breakdownCount !== undefined && sp.breakdownCount > this.sparkCount) this.sparkCount = sp.breakdownCount;
    if (sp.firstSparkDeg !== undefined && Number.isFinite(sp.firstSparkDeg)) this.firstSparkDeg = sp.firstSparkDeg;
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
  /**
   * Trembler showers: close an event this long after its first breakdown even if the timer never
   * reports the contact open, s (an 87° contact lasts ≈ 0.1 s at 150 rpm).
   */
  maxShower: number;
}

const DEFAULTS: SparkCaptureOptions = { tail: 0.5e-3, maxAfterBreakdown: 8e-3, maxShower: 0.25 };

export class SparkCapture {
  private building = new SparkEvent();
  private shown = new SparkEvent();
  private active = false;
  private hasShown = false;
  private lastT = -Infinity;
  private lastEventCycle = NaN;
  /** The timer contact has been seen open since the last trembler event started (a new contact may begin). */
  private contactReopened = true;
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
    this.contactReopened = true;
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
    const sp = s.spark;
    const phase = sp.phase;
    /** Trembler stream: the snapshot reports the timer contact. */
    const timer = sp.timerClosed;
    let restarted = false;
    if (timer === false) this.contactReopened = true;

    if (!this.active) {
      // A new event starts with a dwell, or with the timer contact closing (trembler: only a NEW contact,
      // not the rest of one whose event was force-closed); a discharge only opens one if this cycle has
      // none yet (stream joined mid-spark), not after a forced close of a long discharge.
      const starts =
        timer === undefined
          ? phase === 'charging' || (isDischargePhase(phase) && s.cycle !== this.lastEventCycle)
          : timer && (this.contactReopened || s.cycle !== this.lastEventCycle);
      if (starts) {
        this.building.reset(s.cycle, s.t);
        this.building.shower = timer !== undefined;
        this.active = true;
        this.lastEventCycle = s.cycle;
        this.contactReopened = false;
      } else {
        return false;
      }
    } else if (phase === 'charging' && Number.isFinite(this.building.tBreakdown) && timer !== true) {
      // A new dwell started before the previous event was closed: close it and restart. (A trembler
      // re-charging while the timer contact is still closed is the next spark of the same shower.)
      this.finish();
      this.building.reset(s.cycle, s.t);
      this.building.shower = timer !== undefined;
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
      if (!Number.isFinite(ev2.firstSparkDeg)) ev2.firstSparkDeg = s.thetaDeg;
      changed = true;
    }
    // The discharge (shower) has ended: no conduction, and for a trembler the timer contact is open.
    const quiet = timer === true ? false : timer === false ? !isDischargePhase(phase) : phase === 'done' || phase === 'off';
    if (Number.isFinite(ev2.tBreakdown)) {
      if (!Number.isFinite(ev2.tEnd) && quiet) ev2.tEnd = s.t;
      else if (Number.isFinite(ev2.tEnd) && !quiet && ev2.shower) ev2.tEnd = NaN; // contact bounced closed again
      const tailDone = Number.isFinite(ev2.tEnd) && s.t - ev2.tEnd >= this.opts.tail;
      const tooLong = s.t - ev2.tBreakdown >= (ev2.shower ? this.opts.maxShower : this.opts.maxAfterBreakdown);
      if (tailDone || tooLong) {
        if (!Number.isFinite(ev2.tEnd)) ev2.tEnd = s.t;
        this.finish();
        changed = true;
      }
    } else if (quiet && (timer !== undefined || phase === 'done' || phase === 'off')) {
      // Dwell (contact) ended without breakdown (e.g. insufficient voltage): keep it as a (failed) event.
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
