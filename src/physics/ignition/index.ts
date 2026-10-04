/**
 * Spark-ignition system (public API of src/physics/ignition). SI units.
 *
 * IgnitionSystem composes
 *   coil.ts          — inductive primary/secondary circuit ODEs (dwell, switch-off, ring-up, discharge),
 *   trembler-coil.ts — trembler (vibrator) coil circuit for 'trembler-magneto' ignition: supply
 *                      (source.ts: battery or AC magneto), timer, condenser, vibrator (vibrator.ts),
 *   breakdown.ts     — density-dependent breakdown voltage of the gap,
 *   discharge.ts     — breakdown / capacitive arc / arc / glow phases, extinction, restrike and
 *                      the electrical-to-gas energy transfer,
 *   kernel.ts        — Herweg–Maly spark-kernel growth to hand-off, quench / misfire,
 * and exposes a state object whose top-level fields match EngineSnapshot['spark'].
 *
 * Timing: `step(dt, thetaDeg, command, gas)` takes the crank angle at the END of the step;
 * the previous call's angle is the start. Dwell starts when the crank crosses
 * `command.dwellStartDeg` and the spark fires (primary switch opens) when it crosses
 * `command.sparkDeg` (firing-TDC convention, cycle [−360, 360)). Event times inside a step
 * are found by linear interpolation of θ(t) and the circuit step is split there.
 * `stepTimed(dt, dwellOn, gas)` drives the switch directly (tests, bench rigs).
 * A new dwell starts a new ignition event (per-cycle ledgers and kernel reset).
 *
 * Make spark: closing the switch induces ≈ k·√(L₂/L₁)·V_s (≈ 1.3 kV for the CFR-like coil,
 * ringing to ≈ 1.8 kV) on the plug. If the gap's breakdown voltage at dwell start is lower
 * (only at very low gap density: p*·d ≲ 0.3 bar·mm with p* the pressure reduced to 20 °C,
 * or a hot gap), the circuit produces a "make"
 * spark during the dwell, and its kernel counts as this cycle's kernel (reported phase
 * 'charging'). Coil-on-plug coils usually block it with a high-voltage diode in the
 * secondary: IgnitionSystemOptions.makeSparkDiode (off by default here; on in the cycle model).
 *
 * State carry-over: the state (kernel, energy ledgers, phase) keeps describing the PREVIOUS
 * ignition event until the next dwell starts (a new event begins at switch-on); callers that
 * read it between cycles must track the dwell start (the cycle model does: dwellSeen).
 *
 * ── Trembler-magneto ignition (Ford Model T; spec.type 'trembler-magneto') ─────────────────
 * One IgnitionSystem per cylinder (createIgnitionSystems): its own trembler coil
 * (TremblerMagnetoIgnitionSpec.coil with coilOverrides[cylinder]), gap and kernel, fed through its
 * timer segment by the SHARED supply (battery or magneto; the supply's impedance is lumped into the
 * connected loop, exact because the timer grounds one coil at a time).
 *  - Command: the IgnitionCommand window is the TIMER CONTACT: dwellStartDeg = timer make (the spark
 *    lever, −sparkAdvanceDeg), sparkDeg = timer break (make + contactArcDeg) — tremblerTimerCommand.
 *    The spark itself is an OUTPUT: the vibrator trips when the primary current has built up, the
 *    secondary rings up and the gap breaks down (state.firstSparkDeg); the points re-close and the
 *    coil buzzes, giving a spark TRAIN for as long as the timer contact lasts.
 *  - Drive: `step(dt, thetaDeg, command, gas, omega, thetaEngineDeg)` — omega the crank speed
 *    (rad/s) and thetaEngineDeg the ENGINE angle (cylinder 1, for the magneto phase) at the end of the
 *    step. Defaults: omega from the angle progression, thetaEngineDeg = thetaDeg + firingOffsetDeg
 *    (options). The supply is selected with `ignitionSource` (applied at the next timer make).
 *  - One ignition EVENT = one timer contact: the event (gap counters, kernel, ledgers) starts at the
 *    timer make; points openings/closings never start a new event. While the timer is closed the
 *    kernel's misfire checks are suspended (trainActive) and a quenched kernel is RE-SEEDED by the
 *    next breakdown of the train; misfire is judged after the timer break.
 *  - makeSparkDiode is ignored (a trembler coil has no diode; make sparks are part of its physics).
 *  - state: switchState = the POINTS ('closed' / 'open'), timerClosed, pointsOpen, breakdownCount
 *    (whole train), pointsBreakCount, firstSparkDeg / firstSparkDelay (timer make → first breakdown),
 *    firstTripDelay / firstTripCurrent (timer make → first points opening), conductingTime (summed;
 *    also sparkDuration), sourceEmf, magnetoEmf, trainActive; breakdownDelay = first points opening →
 *    first breakdown.
 */

import type { EngineSpec, IgnitionSystemSpec, SparkPlugSpec, TremblerCoilSpec, TremblerMagnetoIgnitionSpec } from '../core/engine-spec';
import { firingOffsetsDeg } from '../core/engine-spec';
import type { SparkPhase } from '../core/snapshot';
import { breakdownVoltage, type BreakdownOptions } from './breakdown';
import { IgnitionCoil, type CoilOptions, type PrimarySwitchState } from './coil';
import { SparkGap, type DischargeMode, type SparkGapOptions } from './discharge';
import { SparkKernel, type KernelGasState, type KernelStage, type SparkKernelOptions } from './kernel';
import { PrimarySupply, type SupplyKind } from './source';
import { TremblerCoil, type TremblerCoilOptions } from './trembler-coil';

export * from './breakdown';
export * from './coil';
export * from './discharge';
export * from './kernel';
export * from './source';
export * from './trembler-coil';
export * from './vibrator';

/** Crank-angle commands of the ignition driver (firing-TDC convention, degrees). */
export interface IgnitionCommand {
  /**
   * Crank angle at which the primary switch closes (dwell start), deg.
   * Trembler-magneto: the timer contact MAKE (the spark lever).
   */
  dwellStartDeg: number;
  /**
   * Crank angle at which the primary switch opens (spark), deg (= −spark advance).
   * Trembler-magneto: the timer contact BREAK (make + contactArcDeg).
   */
  sparkDeg: number;
}

/**
 * Gas state at the spark plug (inputs of the kernel + gap). SI units.
 * S_L, Markstein length, flame thickness, u′, l_I and σ are injected by the caller
 * (combustion module) — the ignition module does not compute flame speeds.
 */
export interface IgnitionGasState extends KernelGasState {
  /** Mean gas velocity across the gap (channel stretching), m/s (default 0). */
  flowVelocity?: number;
}

/** Construction options. */
export interface IgnitionSystemOptions {
  coil?: Partial<CoilOptions>;
  gap?: Partial<SparkGapOptions>;
  kernel?: Partial<SparkKernelOptions>;
  breakdown?: BreakdownOptions;
  /** Never break down (open-circuit coil test). */
  suppressBreakdown?: boolean;
  /**
   * High-voltage "make-spark" suppression diode in the secondary (German EFU diode of coil-on-plug
   * coils): no breakdown while the primary switch is closed (the switch-on transient cannot fire
   * the plug). Default false (bare coil). The cycle model enables it by default: the real CFR
   * ignition is capacitive discharge (CFR_CDI_IGNITION), which has no dwell and no make spark, and
   * the bare stand-in coil fired a make spark at dwell start above ≈ 4500 rpm (≈ 100° before the
   * commanded spark; validation round 1). Ignored for 'trembler-magneto' ignition.
   */
  makeSparkDiode?: boolean;
  /** Trembler-magneto: 0-based cylinder index (selects coilOverrides[cylinder]; default 0). */
  cylinder?: number;
  /**
   * Trembler-magneto: this cylinder's firing-TDC offset, deg (engine angle = local angle + offset when
   * `step` is given no engine angle; default 0).
   */
  firingOffsetDeg?: number;
  /** Trembler-magneto: initial supply, the dash switch MAG / BAT (default 'magneto'). */
  ignitionSource?: SupplyKind;
  /** Trembler-magneto: circuit / integrator options. */
  trembler?: Partial<TremblerCoilOptions>;
}

/** Kernel part of the ignition state. */
export interface IgnitionKernelState {
  stage: KernelStage;
  /** Kernel radius, m. */
  radius: number;
  /** Kernel volume, m³. */
  volume: number;
  /** Burned (kernel) mass, kg. */
  burnedMass: number;
  /** Kernel gas temperature, K. */
  temperature: number;
  /** Time since the kernel was created (the first breakdown; trembler: the last re-seed), s. */
  age: number;
  /** Herweg–Maly stretch factor I₀. */
  stretchFactor: number;
  /** Kernel burning velocity S_T,k, m/s. */
  turbulentSpeed: number;
  /** Plasma velocity, m/s. */
  plasmaSpeed: number;
  /** dr_k/dt, m/s. */
  growthRate: number;
  /** Abdel-Gayed–Bradley Karlovitz number. */
  karlovitz: number;
  /** Hand-off radius, m. */
  handoffRadius: number;
  /** Hand-off time since the kernel's creation, s (NaN if none). */
  handoffTime: number;
  /** Energy deposited in the kernel by the discharge, J; electrode conduction loss, J. */
  energyDeposited: number;
  energyElectrodeLoss: number;
  /** Reason for a quench ('' if none). */
  quenchReason: string;
}

/** Full ignition state. The first six fields match EngineSnapshot['spark']. */
export interface IgnitionState {
  phase: SparkPhase;
  /** Primary current, A. */
  primaryCurrent: number;
  /** Secondary (plug) voltage magnitude, V. */
  secondaryVoltage: number;
  /** Secondary current magnitude, A. */
  secondaryCurrent: number;
  /** Electrical energy delivered to the gap this cycle, J. */
  energyDelivered: number;
  /** Breakdown voltage required at the current gap-gas density, V. */
  breakdownVoltage: number;

  /** Primary switch (collector) voltage, V. Trembler: the condenser (points) voltage. */
  primaryVoltage: number;
  /** Primary switch state. Trembler: the vibrator points ('closed' / 'open'). */
  switchState: PrimarySwitchState;
  /** Energy stored in the coil (magnetic + capacitive), J. */
  coilEnergy: number;
  /** ½L₁I₁² at the last switch-off (trembler: the last points opening), J. */
  energyAtSwitchOff: number;
  /** Discharge mode and gap voltage/current magnitudes. */
  dischargeMode: DischargeMode;
  gapVoltage: number;
  gapCurrent: number;
  /** Mean electrical power into the gap / into the gas over the last step, W. */
  electricalPower: number;
  gasPower: number;
  /** Cycle energies: to gas, to electrodes (sheaths), radiated; breakdown & capacitive-arc parts, J. */
  energyToGas: number;
  energyToElectrodes: number;
  energyRadiated: number;
  energyBreakdown: number;
  energyCapacitiveArc: number;
  /** Number of breakdowns this cycle (1 + restrikes; trembler: the whole spark train). */
  breakdownCount: number;
  /** Voltage at the last breakdown, V. */
  lastBreakdownVoltage: number;
  /** Time from switch-off (trembler: the first points opening) to the first breakdown, s (NaN if none). */
  breakdownDelay: number;
  /**
   * Discharge duration (first breakdown → extinction, or → now while active), s. Trembler: the summed
   * conducting time of the train (= conductingTime).
   */
  sparkDuration: number;
  /** Peak |V₂| since switch-off (trembler: since the last points opening), V. */
  peakSecondaryVoltage: number;
  kernel: IgnitionKernelState;
  /**
   * True if this cycle's ignition failed: the kernel was quenched, or no breakdown occurred
   * by max(0.5 ms, 4× the ring-up peak time) after switch-off. Trembler: judged only after the timer
   * break (kernel quenched, or no breakdown within 0.5 ms of the break); false while the train is active.
   */
  misfire: boolean;

  // ---- spark train / trembler-magneto (neutral values for 'inductive' where noted) ----
  /** Local crank angle of the first breakdown of the current event, deg (NaN before it / stepTimed). */
  firstSparkDeg: number;
  /** Ignition-clock time at which the current kernel was created, s (NaN: none) — hand-off instant = kernelStartTime + kernel.handoffTime. */
  kernelStartTime: number;
  /** Summed conducting time of the event's discharges, s. */
  conductingTime: number;
  /** Timer contact closed (trembler; inductive: the switch is closed). */
  timerClosed: boolean;
  /** Vibrator points open (trembler; inductive: false). */
  pointsOpen: boolean;
  /** Points openings (vibrator trips) in the event (trembler; inductive: 0). */
  pointsBreakCount: number;
  /** Timer make → first points opening, s (NaN: none; inductive NaN). */
  firstTripDelay: number;
  /** |I₁| at the first points opening (the firing current), A (NaN: none). */
  firstTripCurrent: number;
  /** Timer make → first breakdown, s (NaN: none; inductive: dwell start → first breakdown). */
  firstSparkDelay: number;
  /** EMF of the supply in use at the end of the step, V (inductive: the supply voltage). */
  sourceEmf: number;
  /** Open-circuit magneto EMF at the end of the step, V (0 without a magneto). */
  magnetoEmf: number;
  /** Spark train active: timer closed or gap discharging (inductive: switch closed or discharging). */
  trainActive: boolean;
}

/** Map θ to (θ0, θ0 + 720] for crossing tests. */
function unwrapAfter(theta0: number, theta: number): number {
  let x = theta;
  while (x <= theta0) x += 720;
  while (x > theta0 + 720) x -= 720;
  return x;
}

/** Fraction f ∈ (0, 1] of a step θ0 → θ1 at which the angle θe is crossed, or −1. */
function crossing(theta0: number, theta1: number, thetaE: number): number {
  const t1 = unwrapAfter(theta0, theta1);
  if (t1 - theta0 >= 720) return -1;
  const te = unwrapAfter(theta0, thetaE);
  return te <= t1 ? (te - theta0) / (t1 - theta0) : -1;
}

/** True if θ lies in the (cyclic) half-open window [a, b). */
function inWindow(theta: number, a: number, b: number): boolean {
  const bb = unwrapAfter(a, b);
  const tt = unwrapAfter(a, theta);
  return tt === a + 720 ? true : tt < bb;
}

/** Wrap an angle to the cycle [−360, 360), deg. */
function wrapCycle(theta: number): number {
  return ((((theta + 360) % 720) + 720) % 720) - 360;
}

/** Forward angle progress θ0 → θ1 over one step, deg, in (−360, 360] (θ cyclic mod 720). */
function progress(theta0: number, theta1: number): number {
  let d = (((theta1 - theta0) % 720) + 720) % 720;
  if (d > 360) d -= 720;
  return d;
}

const DEG = Math.PI / 180;

/**
 * The coil of cylinder i (0-based) of a trembler-magneto ignition: `coil` with `coilOverrides[i]`
 * applied (the vibrator merged field by field). Returns `spec.coil` itself without an override.
 */
export function tremblerCoilOf(spec: TremblerMagnetoIgnitionSpec, i: number): TremblerCoilSpec {
  const o = spec.coilOverrides?.[i];
  if (!o) return spec.coil;
  return { ...spec.coil, ...o, vibrator: { ...spec.coil.vibrator, ...o.vibrator } };
}

/**
 * Timer contact window of a trembler-magneto ignition as an IgnitionCommand (local crank angles,
 * wrapped to [−360, 360)): make = −sparkAdvanceDeg with the advance clamped to the lever range
 * timer.advanceRangeDeg, break = make + timer.contactArcDeg.
 */
export function tremblerTimerCommand(spec: TremblerMagnetoIgnitionSpec, sparkAdvanceDeg: number, out?: IgnitionCommand): IgnitionCommand {
  const [lo, hi] = spec.timer.advanceRangeDeg;
  const adv = sparkAdvanceDeg < lo ? lo : sparkAdvanceDeg > hi ? hi : sparkAdvanceDeg;
  const cmd = out ?? { dwellStartDeg: 0, sparkDeg: 0 };
  cmd.dwellStartDeg = wrapCycle(-adv);
  cmd.sparkDeg = wrapCycle(-adv + spec.timer.contactArcDeg);
  return cmd;
}

/**
 * One IgnitionSystem per cylinder of `spec` (cylinder index and firing offset from spec.layout set;
 * trembler-magneto: per-cylinder coils, shared supply data). `opts` apply to every cylinder.
 */
export function createIgnitionSystems(spec: EngineSpec, opts: IgnitionSystemOptions = {}): IgnitionSystem[] {
  const offsets = firingOffsetsDeg(spec);
  const out: IgnitionSystem[] = [];
  for (let i = 0; i < spec.cylinders; i++) {
    out.push(new IgnitionSystem(spec.ignition, spec.sparkPlug, { ...opts, cylinder: i, firingOffsetDeg: offsets[i] ?? 0 }));
  }
  return out;
}

/**
 * Ignition system + spark gap + spark kernel for one cylinder: inductive, or trembler coil fed by a
 * battery / magneto through its timer contact (see the module doc).
 */
export class IgnitionSystem {
  readonly spec: IgnitionSystemSpec;
  readonly plug: SparkPlugSpec;
  /** The coil circuit: IgnitionCoil ('inductive') or TremblerCoil ('trembler-magneto'). */
  readonly coil: IgnitionCoil | TremblerCoil;
  readonly gap: SparkGap;
  readonly kernel: SparkKernel;
  readonly breakdownOptions: BreakdownOptions;
  /** State after the last step (reused object). */
  readonly state: IgnitionState;
  /** Trembler-magneto: the coil circuit (null for 'inductive'). */
  readonly trembler: TremblerCoil | null;
  /** Trembler-magneto: the supply (battery / magneto) feeding this cylinder's coil (null for 'inductive'). */
  readonly supply: PrimarySupply | null;
  /** Trembler-magneto: 0-based cylinder index. */
  readonly cylinder: number;
  /** Trembler-magneto: firing-TDC offset of this cylinder, deg. */
  readonly firingOffsetDeg: number;

  private readonly inductive: IgnitionCoil | null;
  private lastTheta = NaN;
  private dwellOn = false;
  private readonly suppressBase: boolean;
  private readonly makeSparkDiode: boolean;
  private switchOffTime = NaN;
  private firedThisCycle = false;
  private eventStartTime = NaN;
  private firstSparkDeg = NaN;
  private kernelStartTime = NaN;
  // trembler drive
  private lastOmega = NaN;
  private lastThetaEngine = NaN;

  constructor(spec: IgnitionSystemSpec, plug: SparkPlugSpec, opts: IgnitionSystemOptions = {}) {
    this.spec = spec;
    this.plug = plug;
    this.cylinder = opts.cylinder ?? 0;
    this.firingOffsetDeg = opts.firingOffsetDeg ?? 0;
    if (spec.type === 'inductive') {
      this.inductive = new IgnitionCoil(spec, opts.coil);
      this.coil = this.inductive;
      this.trembler = null;
      this.supply = null;
      this.makeSparkDiode = opts.makeSparkDiode ?? false;
    } else {
      this.supply = new PrimarySupply(spec.battery, spec.magneto, opts.ignitionSource ?? 'magneto');
      this.trembler = new TremblerCoil(tremblerCoilOf(spec, this.cylinder), this.supply, spec.timer.contactResistance, opts.trembler);
      this.coil = this.trembler;
      this.inductive = null;
      this.makeSparkDiode = false;
    }
    this.gap = new SparkGap(plug.gap, opts.gap);
    this.suppressBase = opts.suppressBreakdown ?? false;
    this.gap.suppressBreakdown = this.suppressBase;
    this.kernel = new SparkKernel(plug, opts.kernel);
    this.breakdownOptions = opts.breakdown ?? {};
    this.state = {
      phase: 'off',
      primaryCurrent: 0,
      secondaryVoltage: 0,
      secondaryCurrent: 0,
      energyDelivered: 0,
      breakdownVoltage: 0,
      primaryVoltage: this.coil.V1,
      switchState: this.inductive ? 'open' : 'closed',
      coilEnergy: 0,
      energyAtSwitchOff: 0,
      dischargeMode: 'open',
      gapVoltage: 0,
      gapCurrent: 0,
      electricalPower: 0,
      gasPower: 0,
      energyToGas: 0,
      energyToElectrodes: 0,
      energyRadiated: 0,
      energyBreakdown: 0,
      energyCapacitiveArc: 0,
      breakdownCount: 0,
      lastBreakdownVoltage: 0,
      breakdownDelay: NaN,
      sparkDuration: 0,
      peakSecondaryVoltage: 0,
      kernel: {
        stage: 'none',
        radius: 0,
        volume: 0,
        burnedMass: 0,
        temperature: 0,
        age: 0,
        stretchFactor: 0,
        turbulentSpeed: 0,
        plasmaSpeed: 0,
        growthRate: 0,
        karlovitz: 0,
        handoffRadius: 0,
        handoffTime: NaN,
        energyDeposited: 0,
        energyElectrodeLoss: 0,
        quenchReason: '',
      },
      misfire: false,
      firstSparkDeg: NaN,
      kernelStartTime: NaN,
      conductingTime: 0,
      timerClosed: false,
      pointsOpen: false,
      pointsBreakCount: 0,
      firstTripDelay: NaN,
      firstTripCurrent: NaN,
      firstSparkDelay: NaN,
      sourceEmf: this.inductive ? this.inductive.spec.supplyVoltage : 0,
      magnetoEmf: 0,
      trainActive: false,
    };
  }

  /** Simulated time of the ignition clock, s. */
  get time(): number {
    return this.coil.time;
  }

  /** Trembler-magneto: the selected supply (MAG / BAT); a change takes effect at the next timer make. */
  get ignitionSource(): SupplyKind | null {
    return this.supply ? this.supply.requested : null;
  }

  set ignitionSource(kind: SupplyKind | null) {
    if (this.supply && kind) this.supply.request(kind);
  }

  /**
   * Trembler-magneto bench drive for `stepTimed`: engine angle (deg) and constant crank speed (rad/s)
   * the magneto turns at from now on.
   */
  setDrive(thetaEngineDeg: number, omega: number): void {
    this.lastThetaEngine = thetaEngineDeg;
    this.lastOmega = omega;
  }

  /** Start a new ignition event: per-cycle ledgers, gap and kernel are reset. */
  private newCycle(): void {
    this.gap.reset();
    this.kernel.reset();
    this.coil.resetLedger();
    this.switchOffTime = NaN;
    this.firedThisCycle = false;
    this.eventStartTime = this.coil.time;
    this.firstSparkDeg = NaN;
    this.kernelStartTime = NaN;
    if (this.trembler) this.trembler.resetEvent();
  }

  /** Update the gap's gas-dependent parameters. */
  private updateGap(gas: IgnitionGasState): void {
    const k = this.kernel;
    const hot = (k.stage === 'kernel' || k.stage === 'handoff') && k.radius > 0.5 * this.plug.gap;
    const Tgap = hot ? k.temperature : gas.Tu;
    this.gap.pressure = gas.p;
    this.gap.flowVelocity = gas.flowVelocity ?? 0;
    this.gap.breakdownVoltage = breakdownVoltage(this.plug.gap, gas.p, Tgap, gas.X, this.breakdownOptions);
  }

  /** Circuit sub-interval with a fixed switch (trembler: timer) state. */
  private advanceCircuit(h: number, dwellOn: boolean): void {
    if (dwellOn && !this.dwellOn) this.newCycle();
    if (!dwellOn && this.dwellOn) {
      this.switchOffTime = this.coil.time;
      this.firedThisCycle = true;
    }
    this.dwellOn = dwellOn;
    this.gap.suppressBreakdown = this.suppressBase || (this.makeSparkDiode && dwellOn);
    if (this.inductive) {
      if (h > 0) this.inductive.step(h, dwellOn, this.gap);
    } else {
      // zero-length calls apply the timer edge at once
      this.trembler!.step(h > 0 ? h : 0, dwellOn, this.gap);
    }
  }

  /**
   * Trembler: angle/speed interpolation of the supply over the step [t0, t0 + dt] (see `step`).
   */
  private driveSupply(t0: number, dt: number, thetaDeg: number, omega: number | undefined, thetaEngineDeg: number | undefined): void {
    const sup = this.supply!;
    let w1: number;
    if (omega !== undefined) w1 = omega;
    else if (!Number.isNaN(this.lastTheta) && dt > 0) w1 = (progress(this.lastTheta, thetaDeg) * DEG) / dt;
    else w1 = Number.isNaN(this.lastOmega) ? 0 : this.lastOmega;
    const w0 = Number.isNaN(this.lastOmega) ? w1 : this.lastOmega;
    const thE1 = thetaEngineDeg ?? thetaDeg + this.firingOffsetDeg;
    const thE0 = Number.isNaN(this.lastThetaEngine)
      ? thE1 - (0.5 * (w0 + w1) * dt) / DEG
      : thE1 - progress(this.lastThetaEngine, thE1);
    sup.setStep(t0, dt, thE0, thE1, w0, w1);
    this.trembler!.invalidateEmf();
    this.lastOmega = w1;
    this.lastThetaEngine = thE1;
  }

  /**
   * Advance by dt (s) with crank-angle commands; thetaDeg is the crank angle at the END of
   * the step. Returns the (reused) state object.
   * Trembler-magneto: `omega` = crank speed (rad/s) and `thetaEngineDeg` = engine angle (cylinder 1,
   * deg) at the end of the step (defaults: from the angle progression / thetaDeg + firingOffsetDeg).
   */
  step(dt: number, thetaDeg: number, command: IgnitionCommand, gas: IgnitionGasState, omega?: number, thetaEngineDeg?: number): IgnitionState {
    const t0 = this.coil.time;
    const thStart = this.lastTheta;
    this.gap.beginStep();
    this.updateGap(gas);
    if (this.trembler) this.driveSupply(t0, dt, thetaDeg, omega, thetaEngineDeg);
    if (Number.isNaN(this.lastTheta)) {
      // First call: no crossing information; take the switch state from the window.
      const on = inWindow(thetaDeg, command.dwellStartDeg, command.sparkDeg);
      this.advanceCircuit(dt, on);
    } else {
      let fOn = crossing(this.lastTheta, thetaDeg, command.dwellStartDeg);
      let fOff = crossing(this.lastTheta, thetaDeg, command.sparkDeg);
      let done = 0;
      let on = this.dwellOn;
      // process up to two events in time order
      for (let e = 0; e < 2; e++) {
        let f = -1;
        let next = on;
        if (fOn >= 0 && (fOff < 0 || fOn <= fOff)) {
          f = fOn;
          next = true;
          fOn = -1;
        } else if (fOff >= 0) {
          f = fOff;
          next = false;
          fOff = -1;
        }
        if (f < 0) break;
        this.advanceCircuit((f - done) * dt, on);
        done = f;
        on = next;
        this.advanceCircuit(0, on);
      }
      this.advanceCircuit((1 - done) * dt, on);
    }
    this.lastTheta = thetaDeg;
    this.finishStep(t0, dt, gas);
    // local angle of the first breakdown of the event (θ linear over the step)
    const tb = this.gap.firstBreakdownTime;
    if (Number.isNaN(this.firstSparkDeg) && tb >= t0 - 1e-15) {
      if (Number.isNaN(thStart) || !(dt > 0)) this.firstSparkDeg = thetaDeg;
      else {
        const f = Math.min(1, Math.max(0, (tb - t0) / dt));
        this.firstSparkDeg = wrapCycle(thStart + f * progress(thStart, thetaDeg));
      }
      this.state.firstSparkDeg = this.firstSparkDeg;
    }
    return this.state;
  }

  /**
   * Advance by dt (s) with the primary switch (trembler: the timer contact) state given directly.
   * Trembler-magneto: the magneto turns at the speed of the last `step` / `setDrive`.
   */
  stepTimed(dt: number, dwellOn: boolean, gas: IgnitionGasState): IgnitionState {
    const t0 = this.coil.time;
    this.gap.beginStep();
    this.updateGap(gas);
    if (this.trembler) {
      const w = Number.isNaN(this.lastOmega) ? 0 : this.lastOmega;
      const th0 = Number.isNaN(this.lastThetaEngine) ? 0 : this.lastThetaEngine;
      const th1 = th0 + (w * dt) / DEG;
      this.supply!.setStep(t0, dt, th0, th1, w, w);
      this.trembler.invalidateEmf();
      this.lastOmega = w;
      this.lastThetaEngine = wrapCycle(th1);
    }
    this.advanceCircuit(dt, dwellOn);
    this.finishStep(t0, dt, gas);
    return this.state;
  }

  /** Kernel update and state assembly after the circuit has advanced over [t0, t0+dt]. */
  private finishStep(t0: number, dt: number, gas: IgnitionGasState): void {
    const g = this.gap;
    const k = this.kernel;
    const tEnd = this.coil.time;
    const sparkActive = g.conducting || g.awaitingReentry;
    const imp = g.stepImpulseGasEnergy;
    const cont = g.stepContinuousGasEnergy;
    if (this.trembler) {
      // spark train: no misfire decision while the timer contact lasts; a later breakdown re-seeds a
      // quenched kernel
      const train = sparkActive || this.trembler.timerClosed;
      if (imp > 0) {
        if (k.stage === 'quenched' && this.trembler.timerClosed) k.reseed();
        const tb = Number.isNaN(g.stepBreakdownTime) ? t0 : g.stepBreakdownTime;
        const wasNone = k.stage === 'none';
        k.deposit(imp, gas);
        if (wasNone && k.stage === 'kernel') this.kernelStartTime = tb;
        const span = wasNone ? tEnd - tb : dt;
        if (span > 0) k.step(span, cont / span, train, gas);
      } else if (k.stage === 'kernel') {
        k.step(dt, cont / dt, train, gas);
      }
      this.assembleTrembler(dt);
      return;
    }
    if (imp > 0) {
      const tb = Number.isNaN(g.stepBreakdownTime) ? t0 : g.stepBreakdownTime;
      const wasNone = k.stage === 'none';
      k.deposit(imp, gas);
      if (wasNone && k.stage === 'kernel') this.kernelStartTime = tb;
      const span = wasNone ? tEnd - tb : dt;
      if (span > 0) k.step(span, cont / span, sparkActive, gas);
    } else if (k.stage === 'kernel') {
      k.step(dt, cont / dt, sparkActive, gas);
    }
    this.assemble(dt, gas);
  }

  /** Fields common to both ignition types. */
  private assembleCommon(dt: number): void {
    const s = this.state;
    const c = this.coil;
    const g = this.gap;
    const k = this.kernel;
    s.primaryCurrent = c.I1;
    s.secondaryVoltage = Math.abs(c.V2);
    s.secondaryCurrent = Math.abs(c.I2);
    s.energyDelivered = g.energyTotal;
    s.breakdownVoltage = g.breakdownVoltage;
    s.primaryVoltage = c.V1;
    s.coilEnergy = c.storedEnergy();
    s.dischargeMode = g.conducting ? g.mode : 'open';
    s.gapVoltage = g.conducting ? g.voltage : 0;
    s.gapCurrent = g.conducting ? g.current : 0;
    s.electricalPower = dt > 0 ? g.stepElectricalEnergy / dt : 0;
    s.gasPower = dt > 0 ? (g.stepContinuousGasEnergy + g.stepImpulseGasEnergy) / dt : 0;
    s.energyToGas = g.energyToGas;
    s.energyToElectrodes = g.energyToElectrodes;
    s.energyRadiated = g.energyRadiated;
    s.energyBreakdown = g.energyBreakdown;
    s.energyCapacitiveArc = g.energyCapacitiveArc;
    s.breakdownCount = g.breakdownCount;
    s.lastBreakdownVoltage = g.lastBreakdownVoltage;
    // kernel
    const ks = s.kernel;
    ks.stage = k.stage;
    ks.radius = k.radius;
    ks.volume = k.volume;
    ks.burnedMass = k.burnedMass;
    ks.temperature = k.temperature;
    ks.age = k.age;
    ks.stretchFactor = k.stretchFactor;
    ks.turbulentSpeed = k.turbulentSpeed;
    ks.plasmaSpeed = k.plasmaSpeed;
    ks.growthRate = k.growthRate;
    ks.karlovitz = k.karlovitzAGB;
    ks.handoffRadius = k.handoffRadius;
    ks.handoffTime = k.handoffTime;
    ks.energyDeposited = k.energyDeposited;
    ks.energyElectrodeLoss = k.energyElectrodeLoss;
    ks.quenchReason = k.quenchReason;
    // spark-train outputs
    s.firstSparkDeg = this.firstSparkDeg;
    s.kernelStartTime = this.kernelStartTime;
    s.conductingTime = g.conductingTime;
    s.firstSparkDelay = g.firstBreakdownTime - this.eventStartTime;
  }

  private assemble(dt: number, _gas: IgnitionGasState): void {
    const s = this.state;
    const c = this.inductive!;
    const g = this.gap;
    const k = this.kernel;
    this.assembleCommon(dt);
    s.switchState = c.switchState;
    s.energyAtSwitchOff = c.energyAtSwitchOff;
    s.breakdownDelay = g.firstBreakdownTime - this.switchOffTime;
    s.sparkDuration = Number.isNaN(g.firstBreakdownTime)
      ? 0
      : (g.conducting || g.awaitingReentry || Number.isNaN(g.extinctionTime) ? c.time : g.extinctionTime) -
        g.firstBreakdownTime;
    s.peakSecondaryVoltage = c.peakSecondaryVoltage;
    s.timerClosed = c.switchState === 'closed';
    s.trainActive = s.timerClosed || g.conducting || g.awaitingReentry;
    // phase
    if (c.switchState === 'closed') s.phase = 'charging';
    else if (!Number.isNaN(g.stepBreakdownTime)) s.phase = 'breakdown';
    else if (g.conducting) s.phase = g.mode === 'arc' ? 'arc' : 'glow';
    else if (g.awaitingReentry) s.phase = 'glow';
    else if (g.breakdownCount > 0) s.phase = 'done';
    else if (this.firedThisCycle && c.storedEnergy() > 1e-6) s.phase = 'breakdown';
    else s.phase = this.firedThisCycle ? 'done' : 'off';
    // No breakdown: the first ring-up peak is the largest the coil can produce, so once it
    // is well past without a breakdown the available voltage was below the required one.
    const sinceOff = c.time - this.switchOffTime;
    const tPeak = c.peakSecondaryVoltageTime - this.switchOffTime;
    const noSpark = this.firedThisCycle && g.breakdownCount === 0 && sinceOff > Math.max(0.5e-3, 4 * tPeak);
    s.misfire = k.misfire || noSpark;
  }

  private assembleTrembler(dt: number): void {
    const s = this.state;
    const c = this.trembler!;
    const g = this.gap;
    const k = this.kernel;
    this.assembleCommon(dt);
    s.switchState = c.pointsOpen ? 'open' : 'closed';
    s.energyAtSwitchOff = c.energyAtTrip;
    s.breakdownDelay = g.firstBreakdownTime - c.firstTripTime;
    s.sparkDuration = g.conductingTime;
    s.peakSecondaryVoltage = c.peakSecondaryVoltage;
    s.timerClosed = c.timerClosed;
    s.pointsOpen = c.pointsOpen;
    s.pointsBreakCount = c.tripCount;
    s.firstTripDelay = c.firstTripTime - this.eventStartTime;
    s.firstTripCurrent = c.firstTripCurrent;
    s.sourceEmf = c.emfNow();
    s.magnetoEmf = this.supply!.magnetoEmfAt(c.time);
    const discharging = g.conducting || g.awaitingReentry;
    s.trainActive = c.timerClosed || discharging;
    // phase: 'breakdown' in the step of a breakdown and during a ring-up (points open, no breakdown
    // since they opened); 'charging' while the timer contact lasts otherwise
    if (!Number.isNaN(g.stepBreakdownTime)) s.phase = 'breakdown';
    else if (g.conducting) s.phase = g.mode === 'arc' ? 'arc' : 'glow';
    else if (g.awaitingReentry) s.phase = 'glow';
    else if (c.timerClosed)
      s.phase = c.pointsOpen && !(g.lastBreakdownTime >= c.lastTripTime) && c.storedEnergy() > 1e-6 ? 'breakdown' : 'charging';
    else s.phase = g.breakdownCount > 0 || this.firedThisCycle ? 'done' : 'off';
    // misfire: only once the train is over (timer open): kernel quenched, or no breakdown at all
    // within 0.5 ms of the timer break (UNVERIFIED margin; a timer-break spark follows within µs)
    const noSpark = this.firedThisCycle && g.breakdownCount === 0 && c.time - this.switchOffTime > 0.5e-3;
    s.misfire = !c.timerClosed && (k.misfire || noSpark);
  }
}
