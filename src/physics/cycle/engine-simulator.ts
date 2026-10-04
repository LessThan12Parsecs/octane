/**
 * EngineSimulator: the physics CycleModel behind the worker's SimulatorLike interface
 * (src/worker/simulator-like.ts) — snapshot scheduling and the EngineSnapshot / CycleSummary
 * render contract (src/physics/core/snapshot.ts).
 *
 * Snapshot cadence:
 *  - every options.snapshotEveryDeg crank degrees (engine angle) on a grid anchored at θ = −360;
 *  - an event snapshot at every cylinder's spark command (primary switch-off; trembler-magneto: the
 *    timer make, each vibrator trip and the timer break) and knock onset;
 *  - dense time-based samples (every cylinder's): inductive — 5 µs from switch-off until 100 µs after
 *    the first breakdown (or 0.5 ms without breakdown), then 50 µs until the discharge has ended (≤ 5 ms);
 *    trembler-magneto spark train — 5 µs for DENSE_TRIP_WINDOW (150 µs) after each points opening and
 *    after the timer break (ring-up, breakdown, early discharge), DENSE_TRAIN_DT (100 µs) otherwise while
 *    the train is active (timer contact closed or the gap conducting); knock ringing: 10 µs for the
 *    first 3 ms (multi-cylinder engines: 20 µs for 2 ms, DENSE_KNOCK_DT_MULTI).
 * Multi-cylinder engines (spec.cylinders > 1): EngineSnapshot.cylinders holds every cylinder at its
 * local angle (CylinderSnapshot); the top-level per-cylinder fields are cylinder 0's; gasTorque is the
 * engine's (sum), and frictionTorque / loadTorque / vehicleSpeed are added (plus magnetoEmf and
 * firingCylinder for a trembler-magneto ignition, whose spark fields carry the train: breakdownCount,
 * pointsOpen, timerClosed, firstSparkDeg, primaryVoltage). A single-cylinder spec (the CFR) gets exactly
 * the former snapshot (no cylinders array, no extra fields).
 * Conventions kept from the mock simulator (the front-end relies on them):
 * temperatureUnburned = temperatureMean and temperatureBurned = 0 without a burned zone;
 * flame.radius at 'done' covers the chamber; laminar/turbulent speeds 0 after burn-out;
 * rodAngle has the sign of sin θ (no pin offset); intake/exhaust mass flow positive into /
 * out of the cylinder; spark fields 'off' until this cycle's dwell starts; pressure includes the
 * synthesised knock oscillation (thermodynamics never sees it); knockEndGasFraction 0 without
 * knock; isfc NaN when the net work ≤ 0; pmep positive = pumping loss.
 */
import type { EngineSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import type { CycleSummary, CylinderPhase, CylinderSnapshot, EngineSnapshot } from '../core/snapshot';
import { DEG } from '../core/constants';
import { NS, SP } from '../core/species';
import { newCrankTrainState } from '../mechanics';
import { newChamberResult } from '../combustion';
import { magnetoEmf } from '../ignition';
import { CycleModel, I_OM, MODE_BURNED, MODE_OPEN, MODE_TWO, swapNO, type CycleTrace } from './cycle-model';
import { TIME_TARGET_REL_TOL } from './cycle-common';
import type { Cylinder } from './cylinder';
import type { CycleModelOptions } from './options';

/** Snapshot cadence options (structurally the worker's SimulatorOptions). */
export interface SnapshotOptions {
  /** Snapshot cadence in crank degrees (e.g. 0.5). */
  snapshotEveryDeg: number;
  /** Max simulated seconds to buffer ahead of playback (used by the worker host only). */
  bufferAheadSeconds?: number;
}

/** Dense sampling intervals, s. */
export const DENSE_SPARK_DT = 5e-6;
export const DENSE_SPARK_AFTER_BREAKDOWN = 100e-6;
export const DENSE_DISCHARGE_DT = 50e-6;
export const DENSE_DISCHARGE_MAX = 5e-3;
export const DENSE_KNOCK_DT = 10e-6;
export const DENSE_KNOCK_DURATION = 3e-3;
/**
 * Knock ringing of a MULTI-cylinder engine: 20 µs for 2 ms (every cylinder may knock in every cycle — four
 * 10 µs × 3 ms windows doubled the snapshot stream of the Model T at WOT; 20 µs still gives ≈ 13 samples per
 * period of the L-head's ≈ 3.8 kHz fundamental). MAPO itself is sampling-independent (cycle model).
 */
export const DENSE_KNOCK_DT_MULTI = 20e-6;
export const DENSE_KNOCK_DURATION_MULTI = 2e-3;
/** Trembler spark train: 5 µs (DENSE_SPARK_DT) samples this long after each points opening / the timer break, s. */
export const DENSE_TRIP_WINDOW = 150e-6;
/** Trembler spark train: sample interval otherwise while the train is active, s. */
export const DENSE_TRAIN_DT = 100e-6;
/**
 * Dense-sampling targets lie more than DENSE_TIME_REL_TOL·max(1, |t|) ahead of the model time t (numerical):
 * CycleModel.stepUntil treats a target within TIME_TARGET_REL_TOL·max(1, |t|) of t as reached and returns
 * without stepping, so a closer window end or grid point came back on every call — the snapshot stream stalled
 * once t ≳ 100 s (code review: reproduced at t ≈ 3e4 s). Twice stepUntil's tolerance covers the rounding of the
 * candidate times. Below ≈ 50 s (window ends) / 250 s (5 µs grid) the former absolute margins (1e-12 s, 1e-6
 * of the grid interval) are the larger ones, so ordinary runs sample at exactly the former times.
 */
export const DENSE_TIME_REL_TOL = 2 * TIME_TARGET_REL_TOL;

type FlameStage = EngineSnapshot['flame']['stage'];

/** Per-cylinder part of a snapshot with this cylinder's gas torque (built by cylinderSnapshot). */
type CylinderPart = Omit<CylinderSnapshot, 'index' | 'thetaDeg' | 'cycle'>;

export class EngineSimulator {
  readonly model: CycleModel;
  readonly snapshotEveryDeg: number;
  private initialPending = true;
  private lastTime = 0;
  private readonly cts = newCrankTrainState();
  private readonly fg = newChamberResult();
  private readonly Nb = new Float64Array(NS);

  /**
   * @param spec engine
   * @param op initial operating point
   * @param options snapshot cadence + any CycleModelOptions overrides
   */
  constructor(spec: EngineSpec, op: OperatingPoint, options: SnapshotOptions & Partial<CycleModelOptions>) {
    const d = options.snapshotEveryDeg;
    this.snapshotEveryDeg = d > 0 ? Math.min(30, Math.max(0.01, d)) : 0.5;
    this.model = new CycleModel(spec, op, options);
  }

  /** Simulated time of the most recently returned snapshot, s. */
  get time(): number {
    return this.lastTime;
  }

  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.model.setOperatingPoint(patch);
  }

  reset(): void {
    this.model.reset();
    this.initialPending = true;
    this.lastTime = 0;
  }

  drainCycleSummaries(): CycleSummary[] {
    return this.model.drainSummaries();
  }

  /** Run n complete engine cycles without snapshots and return their summaries (validation hook). */
  runCycles(n: number): CycleSummary[] {
    const s = this.model.runCycles(n);
    this.lastTime = this.model.t;
    return s;
  }

  /** Record a full-resolution trace of the given cycle (default: the next complete one). */
  recordTrace(cycle?: number): void {
    this.model.recordTrace(cycle);
  }

  /** The recorded trace (null before recordTrace). */
  get trace(): CycleTrace | null {
    return this.model.trace;
  }

  advanceToNextSnapshot(): EngineSnapshot {
    const m = this.model;
    if (this.initialPending) {
      this.initialPending = false;
    } else {
      const th = m.theta;
      const d = this.snapshotEveryDeg;
      let thTarget = -360 + (Math.floor((th + 360) / d + 1e-7) + 1) * d;
      if (thTarget > 360) thTarget = 360;
      // snapshot exactly at every cylinder's spark command (trembler: timer make; engine angle of its local angle)
      for (const c of m.cylinders) {
        if (!c.ign || !Number.isNaN(c.tSparkCmd)) continue;
        const sp = c.engineAngle(c.sparkEventDeg);
        if (sp > th + 1e-9 && sp < thTarget) thTarget = sp;
      }
      m.stepUntil(this.nextDenseTime(), thTarget);
    }
    this.lastTime = m.t;
    return this.makeSnapshot();
  }

  /**
   * Next dense-sampling time after the current time over all cylinders (Infinity if none). Every candidate lies
   * more than DENSE_TIME_REL_TOL·max(1, |t|) ahead of t (and at least the former absolute margins), so stepUntil
   * always takes a step towards it.
   */
  private nextDenseTime(): number {
    const m = this.model;
    const t = m.t;
    const tol = DENSE_TIME_REL_TOL * Math.max(1, Math.abs(t));
    let next = Infinity;
    for (const c of m.cylinders) next = Math.min(next, c.trembler ? denseTimeTrembler(c, t, tol) : denseTimeOf(c, t, tol));
    return next;
  }

  private makeSnapshot(): EngineSnapshot {
    const m = this.model;
    const cyls = m.cylinders;
    if (cyls.length === 1) return this.singleCylinderSnapshot();
    const om = m.y[I_OM];
    const alpha = m.op.speedMode === 'free' ? m.alpha : 0;
    const pcc = m.op.ambientPressure + m.opts.crankcaseGaugePressure;
    const list: CylinderSnapshot[] = [];
    let gas = 0;
    let other = 0;
    let firing = -1;
    for (const c of cyls) {
      const part = this.cylinderSnapshot(c);
      const th = c.theta;
      // (gas torque from the THERMODYNAMIC pressure, see singleCylinderSnapshot)
      const cts = c.dyn.evaluate(th * DEG, om, alpha, c.p, pcc, this.cts);
      part.gasTorque = cts.gasTorque;
      gas += cts.gasTorque;
      other += cts.inertiaTorque + cts.gravityTorque;
      if (firing < 0 && c.trembler && c.ign && c.ign.state.timerClosed) firing = c.index;
      list.push({ ...part, index: c.index, thetaDeg: th >= 360 ? th - 720 : th, cycle: c.cycle });
    }
    // whole-engine friction (negative = opposing) with every cylinder at its throw phase
    const fric = m.friction.torqueCylinders(m.crank.update(m.theta * DEG), om);
    let net = gas + other + fric;
    const free = m.op.speedMode === 'free';
    const ld = m.load;
    const tLoad = free ? ld.torque(om) : 0;
    if (free) net -= tLoad;
    const c0 = list[0];
    const th = m.theta;
    const snap: EngineSnapshot = {
      t: m.t,
      cycle: m.cycle,
      thetaDeg: th >= 360 ? th - 720 : th,
      rpm: m.rpm,
      pistonDisplacement: c0.pistonDisplacement,
      clearanceHeight: c0.clearanceHeight,
      rodAngle: c0.rodAngle,
      intakeLift: c0.intakeLift,
      exhaustLift: c0.exhaustLift,
      phase: c0.phase,
      volume: c0.volume,
      pressure: c0.pressure,
      temperatureMean: c0.temperatureMean,
      temperatureUnburned: c0.temperatureUnburned,
      temperatureBurned: c0.temperatureBurned,
      massFractionBurned: c0.massFractionBurned,
      mass: c0.mass,
      heatReleaseRate: c0.heatReleaseRate,
      heatLossRate: c0.heatLossRate,
      flame: c0.flame,
      spark: c0.spark,
      intakeMassFlow: c0.intakeMassFlow,
      exhaustMassFlow: c0.exhaustMassFlow,
      intakeManifoldPressure: m.pInt,
      exhaustManifoldPressure: m.pExh,
      knock: c0.knock,
      burnedComposition: c0.burnedComposition,
      gasTorque: gas,
      netTorque: net,
      frictionTorque: -fric,
      // load: the load model in free mode; in fixed mode the dynamometer absorbs the net crank torque
      // (it holds α = 0)
      loadTorque: free ? tLoad : net,
      cylinders: list,
    };
    // (the car's road speed: k ω in gear, its coasting speed in neutral)
    if (free && (ld.kind === 'vehicle' || ld.kind === 'neutral')) snap.vehicleSpeed = m.vehicleSpeed();
    const ig = m.spec.ignition;
    if (ig.type === 'trembler-magneto') {
      snap.magnetoEmf = magnetoEmf(ig.magneto, th, om);
      snap.firingCylinder = firing;
    }
    return snap;
  }

  /** The single-cylinder snapshot (the CFR contract, unchanged). */
  private singleCylinderSnapshot(): EngineSnapshot {
    const m = this.model;
    const c = m.c0;
    const part = this.cylinderSnapshot(c);
    const th = m.theta;
    const thDeg = th >= 360 ? th - 720 : th;
    // mechanics
    const om = m.y[I_OM];
    const alpha = m.op.speedMode === 'free' ? m.alpha : 0;
    const pcc = m.op.ambientPressure + m.opts.crankcaseGaugePressure;
    // gas torque from the THERMODYNAMIC pressure: the synthesised knock field is a sum of rigid-wall
    // modes with zero mean over the piston face (no net force; the free-speed dynamics use m.p too) —
    // the pickup pressure made the reported torque ring by up to ±9 N m (validation round 2)
    const cts = c.dyn.evaluate(th * DEG, om, alpha, c.p, pcc, this.cts);
    let net = cts.gasTorque + cts.inertiaTorque + cts.gravityTorque + m.friction.torque(c.ks.dxdTheta, om);
    if (m.op.speedMode === 'free') net -= m.load.torque(om);
    return {
      t: m.t,
      cycle: m.cycle,
      thetaDeg: thDeg,
      rpm: m.rpm,
      pistonDisplacement: part.pistonDisplacement,
      clearanceHeight: part.clearanceHeight,
      rodAngle: part.rodAngle,
      intakeLift: part.intakeLift,
      exhaustLift: part.exhaustLift,
      phase: part.phase,
      volume: part.volume,
      pressure: part.pressure,
      temperatureMean: part.temperatureMean,
      temperatureUnburned: part.temperatureUnburned,
      temperatureBurned: part.temperatureBurned,
      massFractionBurned: part.massFractionBurned,
      mass: part.mass,
      heatReleaseRate: part.heatReleaseRate,
      heatLossRate: part.heatLossRate,
      flame: part.flame,
      spark: part.spark,
      intakeMassFlow: part.intakeMassFlow,
      exhaustMassFlow: part.exhaustMassFlow,
      intakeManifoldPressure: m.pInt,
      exhaustManifoldPressure: m.pExh,
      knock: part.knock,
      burnedComposition: part.burnedComposition,
      gasTorque: cts.gasTorque,
      netTorque: net,
    };
  }

  private flameStage(c: Cylinder): FlameStage {
    if (c.mode === MODE_OPEN) return c.burnDone ? 'done' : 'none';
    if (c.burnDone) return 'done';
    if (c.flameActive) return c.frontAtWalls ? 'burnout' : 'turbulent';
    if (c.ign && c.dwellSeen && c.ign.state.kernel.stage === 'kernel') return 'kernel';
    if (c.mode === MODE_TWO && (this.model.opts.combustionModel === 'wiebe' || c.knockOnset)) return 'turbulent';
    return 'none';
  }

  /** Kinematics, gas state, flame, spark, gas exchange, knock and composition of one cylinder (local angle). */
  private cylinderSnapshot(c: Cylinder): CylinderPart {
    const spec = this.model.spec;
    const th = c.theta;
    const ks = c.kin.evaluate(th * DEG, c.ks);
    const osc = c.knockOscillation();
    const pressure = c.p + osc;
    const closed = c.mode !== MODE_OPEN;
    const hasBurned = c.mode === MODE_TWO || c.mode === MODE_BURNED;
    const stage = this.flameStage(c);
    // temperatures
    const Tmean = c.T;
    const Tu = c.mode === MODE_BURNED ? 0 : c.mode === MODE_TWO ? c.Tu : Tmean;
    const Tb = hasBurned ? c.Tb : 0;
    // flame
    const fc = spec.sparkPlug.gapCenter;
    let radius = 0;
    let area = 0;
    let SL = 0;
    let ST = 0;
    const kst = c.ign && c.dwellSeen ? c.ign.state.kernel : null;
    if (stage === 'kernel' && kst) {
      radius = c.flameRadius();
      c.chamber.evaluate(radius, c.h, this.fg);
      area = Math.max(0, this.fg.frontArea);
      SL = c.SL;
      ST = kst.turbulentSpeed;
    } else if (stage === 'turbulent' || stage === 'burnout') {
      radius = c.flameActive ? c.rf : c.rb;
      area = c.flameActive ? c.Af : 0;
      if (!c.flameActive && radius > 0) {
        c.chamber.evaluate(radius, c.h, this.fg);
        area = Math.max(0, this.fg.frontArea);
      }
      SL = c.SL;
      ST = c.flameActive ? c.burnSpeed : area > 0 && c.rhoU > 0 ? c.burnRateStep / (c.rhoU * area) : 0;
    } else if (stage === 'done') {
      radius = c.chamber.maxRadius(c.h);
    }
    // phase
    let phase: CylinderPhase;
    if (!closed) phase = 'gas-exchange';
    else if (stage === 'kernel' || stage === 'turbulent' || stage === 'burnout') phase = 'combustion';
    else if (c.burnDone || th > 0) phase = 'expansion';
    else phase = 'compression';
    // composition
    const comp: EngineSnapshot['burnedComposition'] = { CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0 };
    let X: Float64Array;
    if (hasBurned && this.model.y[c.ix.MB] > 0) {
      const r = c.closure.eq.result;
      const Nb = this.Nb;
      for (let k = 0; k < NS; k++) Nb[k] = r.N[k];
      const nb = r.nTotal;
      swapNO(Nb, c.burnedNOFraction() * nb);
      let n = 0;
      for (let k = 0; k < NS; k++) n += Nb[k];
      if (n > 0) {
        for (let k = 0; k < NS; k++) Nb[k] /= n;
        X = Nb;
      } else X = c.closure.Xu; // (no equilibrium result yet: never 0/0)
    } else if (closed) {
      X = c.closure.Xu;
    } else {
      X = c.cyl.state.X;
    }
    comp.CO2 = X[SP.CO2];
    comp.H2O = X[SP.H2O];
    comp.CO = X[SP.CO];
    comp.O2 = X[SP.O2];
    comp.H2 = X[SP.H2];
    comp.OH = X[SP.OH];
    comp.H = X[SP.H];
    comp.O = X[SP.O];
    comp.NO = X[SP.NO];
    comp.N2 = X[SP.N2];
    // spark
    const spark: EngineSnapshot['spark'] = {
      phase: 'off',
      primaryCurrent: 0,
      secondaryVoltage: 0,
      secondaryCurrent: 0,
      energyDelivered: 0,
      breakdownVoltage: 0,
    };
    if (c.ign) {
      const s = c.ign.state;
      spark.breakdownVoltage = s.breakdownVoltage;
      if (c.dwellSeen) {
        spark.phase = s.phase;
        spark.primaryCurrent = s.primaryCurrent;
        spark.secondaryVoltage = s.secondaryVoltage;
        spark.secondaryCurrent = s.secondaryCurrent;
        spark.energyDelivered = s.energyDelivered;
      }
      if (c.trembler) {
        // the spark train of this cycle's timer contact (counters 0 / NaN until the make)
        spark.breakdownCount = c.dwellSeen ? s.breakdownCount : 0;
        spark.pointsOpen = s.pointsOpen;
        spark.timerClosed = s.timerClosed;
        spark.firstSparkDeg = c.dwellSeen ? s.firstSparkDeg : Number.NaN;
        spark.primaryVoltage = s.primaryVoltage;
      }
    }
    return {
      pistonDisplacement: ks.x,
      clearanceHeight: ks.clearanceHeight,
      rodAngle: ks.beta,
      intakeLift: c.ivLift.lift(th),
      exhaustLift: c.evLift.lift(th),
      phase,
      volume: ks.volume,
      pressure,
      temperatureMean: Tmean,
      temperatureUnburned: Tu,
      temperatureBurned: Tb,
      massFractionBurned: c.xb,
      mass: c.mCyl,
      heatReleaseRate: c.hrr,
      heatLossRate: c.Qwall,
      flame: {
        stage,
        radius,
        center: [fc[0], fc[1], fc[2]],
        area,
        laminarSpeed: stage === 'done' ? 0 : SL,
        turbulentSpeed: stage === 'done' ? 0 : ST,
        turbulenceIntensity: c.uPrime,
      },
      spark,
      intakeMassFlow: closed ? 0 : c.mdotIv,
      exhaustMassFlow: closed ? 0 : c.mdotEv,
      knock: {
        // (0 when the delay model does not cover the fuel — the integral is then disarmed)
        integral: c.knockProgress,
        autoignited: c.knockOnset,
        oscillation: osc,
      },
      burnedComposition: comp,
      gasTorque: 0,
    };
  }
}

/**
 * Next dense-sampling time of one cylinder after t (Infinity if none): 5 µs from its spark command
 * (switch-off) until 100 µs after the first breakdown (0.5 ms without one), then 50 µs while the
 * discharge lasts (≤ 5 ms after switch-off); knock ringing after its onset (denseKnockTime). A trembler
 * spark train is sampled by denseTimeTrembler instead. `tol`: DENSE_TIME_REL_TOL·max(1, |t|) (see there).
 */
function denseTimeOf(c: Cylinder, t: number, tol: number): number {
  let next = Infinity;
  const ign = c.ign;
  if (ign && c.dwellSeen && !Number.isNaN(c.tSparkCmd)) {
    const ts = c.tSparkCmd;
    const s = ign.state;
    const bd = s.breakdownDelay;
    // 5 µs from switch-off until 100 µs after the first breakdown (0.5 ms without one) ...
    const tFineEnd = Number.isNaN(bd) ? ts + 0.5e-3 : ts + bd + DENSE_SPARK_AFTER_BREAKDOWN;
    if (t < tFineEnd - windowMargin(tol)) {
      next = Math.min(nextGridTime(ts, DENSE_SPARK_DT, t, tol), tFineEnd);
    } else {
      // ... then 50 µs while the discharge lasts (≤ 5 ms after switch-off)
      const active = s.phase === 'breakdown' || s.phase === 'arc' || s.phase === 'glow';
      if (active && t < ts + DENSE_DISCHARGE_MAX) next = nextGridTime(tFineEnd, DENSE_DISCHARGE_DT, t, tol);
    }
  }
  return Math.min(next, denseKnockTime(c, t, tol));
}

/** Next knock-ringing sample of a cylinder after t (10 µs for 3 ms after its onset; multi-cylinder 20 µs for 2 ms). */
function denseKnockTime(c: Cylinder, t: number, tol: number): number {
  if (!c.knockOnset || c.mode === MODE_OPEN) return Infinity;
  const tk = c.tKnockOnset;
  const multi = c.e.cylinders.length > 1;
  const dt = multi ? DENSE_KNOCK_DT_MULTI : DENSE_KNOCK_DT;
  if (t < tk + (multi ? DENSE_KNOCK_DURATION_MULTI : DENSE_KNOCK_DURATION)) return nextGridTime(tk, dt, t, tol);
  return Infinity;
}

/**
 * How far a dense window's end must lie ahead of t to be returned as a target, s: the former absolute 1e-12 s,
 * or the relative tolerance `tol` once it is larger (t > 50 s).
 */
function windowMargin(tol: number): number {
  return tol > 1e-12 ? tol : 1e-12;
}

/**
 * First point of the grid t0 + k·dt more than max(1e-6·dt, tol) after t (exact arithmetic), s. For
 * tol ≤ 1e-6·dt this is the former expression, operation for operation (bit-identical sample times).
 */
function nextGridTime(t0: number, dt: number, t: number, tol: number): number {
  const f = tol / dt;
  return t0 + (Math.floor((t - t0) / dt + (f > 1e-6 ? f : 1e-6)) + 1) * dt;
}

/**
 * Next dense-sampling time of a trembler-magneto cylinder after t (Infinity if none): 5 µs for
 * DENSE_TRIP_WINDOW after each vibrator trip (points opening: ring-up, breakdown and the early discharge)
 * and after the timer break, DENSE_TRAIN_DT otherwise while the spark train is active (timer contact closed
 * or the gap conducting), anchored at the timer make; knock ringing after its onset (denseKnockTime). The
 * trip instants come from the ignition clock (TremblerCoil.lastTripTime), converted to model time with the
 * clock offset at the last step (the cylinder splits advance the ignition with every step). `tol`:
 * DENSE_TIME_REL_TOL·max(1, |t|) (see there).
 */
function denseTimeTrembler(c: Cylinder, t: number, tol: number): number {
  let next = Infinity;
  const ign = c.ign;
  const tr = ign ? ign.trembler : null;
  if (ign && tr && c.dwellSeen && !Number.isNaN(c.tSparkCmd)) {
    const off = c.e.t - ign.time; // model time − ignition-clock time
    const wm = windowMargin(tol);
    // 5 µs windows after the latest trip of this event and after the timer break
    for (let k = 0; k < 2; k++) {
      const t0 = k === 0 ? (tr.tripCount > 0 ? tr.lastTripTime + off : Number.NaN) : c.tTimerBreak;
      if (Number.isNaN(t0)) continue;
      const tEnd = t0 + DENSE_TRIP_WINDOW;
      if (t >= t0 - wm && t < tEnd - wm) next = Math.min(next, nextGridTime(t0, DENSE_SPARK_DT, t, tol), tEnd);
    }
    if (ign.state.trainActive) next = Math.min(next, nextGridTime(c.tSparkCmd, DENSE_TRAIN_DT, t, tol));
  }
  return Math.min(next, denseKnockTime(c, t, tol));
}
