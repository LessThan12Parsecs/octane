/**
 * Spark-ignition system (public API of src/physics/ignition). SI units.
 *
 * IgnitionSystem composes
 *   coil.ts      — primary/secondary circuit ODEs (dwell, switch-off, ring-up, discharge),
 *   breakdown.ts — density-dependent breakdown voltage of the gap,
 *   discharge.ts — breakdown / capacitive arc / arc / glow phases, extinction, restrike and
 *                  the electrical-to-gas energy transfer,
 *   kernel.ts    — Herweg–Maly spark-kernel growth to hand-off, quench / misfire,
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
 */

import type { IgnitionSystemSpec, SparkPlugSpec } from '../core/engine-spec';
import type { SparkPhase } from '../core/snapshot';
import { breakdownVoltage, type BreakdownOptions } from './breakdown';
import { IgnitionCoil, type CoilOptions, type PrimarySwitchState } from './coil';
import { SparkGap, type DischargeMode, type SparkGapOptions } from './discharge';
import { SparkKernel, type KernelGasState, type KernelStage, type SparkKernelOptions } from './kernel';

export * from './breakdown';
export * from './coil';
export * from './discharge';
export * from './kernel';

/** Crank-angle commands of the ignition driver (firing-TDC convention, degrees). */
export interface IgnitionCommand {
  /** Crank angle at which the primary switch closes (dwell start), deg. */
  dwellStartDeg: number;
  /** Crank angle at which the primary switch opens (spark), deg (= −spark advance). */
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
   * commanded spark; validation round 1).
   */
  makeSparkDiode?: boolean;
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
  /** Time since the first breakdown, s. */
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
  /** Hand-off time since breakdown, s (NaN if none). */
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

  /** Primary switch (collector) voltage, V. */
  primaryVoltage: number;
  switchState: PrimarySwitchState;
  /** Energy stored in the coil (magnetic + capacitive), J. */
  coilEnergy: number;
  /** ½L₁I₁² at the last switch-off, J. */
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
  /** Number of breakdowns this cycle (1 + restrikes). */
  breakdownCount: number;
  /** Voltage at the last breakdown, V. */
  lastBreakdownVoltage: number;
  /** Time from switch-off to the first breakdown, s (NaN if none). */
  breakdownDelay: number;
  /** Discharge duration (first breakdown → extinction, or → now while active), s. */
  sparkDuration: number;
  /** Peak |V₂| since switch-off, V. */
  peakSecondaryVoltage: number;
  kernel: IgnitionKernelState;
  /**
   * True if this cycle's ignition failed: the kernel was quenched, or no breakdown occurred
   * by max(0.5 ms, 4× the ring-up peak time) after switch-off.
   */
  misfire: boolean;
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

/**
 * Inductive ignition system + spark gap + spark kernel for one cylinder.
 */
export class IgnitionSystem {
  readonly spec: IgnitionSystemSpec;
  readonly plug: SparkPlugSpec;
  readonly coil: IgnitionCoil;
  readonly gap: SparkGap;
  readonly kernel: SparkKernel;
  readonly breakdownOptions: BreakdownOptions;
  /** State after the last step (reused object). */
  readonly state: IgnitionState;

  private lastTheta = NaN;
  private dwellOn = false;
  private readonly suppressBase: boolean;
  private readonly makeSparkDiode: boolean;
  private switchOffTime = NaN;
  private firedThisCycle = false;

  constructor(spec: IgnitionSystemSpec, plug: SparkPlugSpec, opts: IgnitionSystemOptions = {}) {
    this.spec = spec;
    this.plug = plug;
    if (spec.type !== 'inductive') throw new Error(`IgnitionSystem: ignition type '${spec.type}' is not implemented yet`);
    this.coil = new IgnitionCoil(spec, opts.coil);
    this.gap = new SparkGap(plug.gap, opts.gap);
    this.suppressBase = opts.suppressBreakdown ?? false;
    this.makeSparkDiode = opts.makeSparkDiode ?? false;
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
      switchState: 'open',
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
    };
  }

  /** Simulated time of the ignition clock, s. */
  get time(): number {
    return this.coil.time;
  }

  /** Start a new ignition event: per-cycle ledgers, gap and kernel are reset. */
  private newCycle(): void {
    this.gap.reset();
    this.kernel.reset();
    this.coil.resetLedger();
    this.switchOffTime = NaN;
    this.firedThisCycle = false;
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

  /** Circuit sub-interval with a fixed switch state. */
  private advanceCircuit(h: number, dwellOn: boolean): void {
    if (dwellOn && !this.dwellOn) this.newCycle();
    if (!dwellOn && this.dwellOn) {
      this.switchOffTime = this.coil.time;
      this.firedThisCycle = true;
    }
    this.dwellOn = dwellOn;
    this.gap.suppressBreakdown = this.suppressBase || (this.makeSparkDiode && dwellOn);
    if (h > 0) this.coil.step(h, dwellOn, this.gap);
  }

  /**
   * Advance by dt (s) with crank-angle commands; thetaDeg is the crank angle at the END of
   * the step. Returns the (reused) state object.
   */
  step(dt: number, thetaDeg: number, command: IgnitionCommand, gas: IgnitionGasState): IgnitionState {
    const t0 = this.coil.time;
    this.gap.beginStep();
    this.updateGap(gas);
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
    return this.state;
  }

  /** Advance by dt (s) with the primary switch state given directly. */
  stepTimed(dt: number, dwellOn: boolean, gas: IgnitionGasState): IgnitionState {
    const t0 = this.coil.time;
    this.gap.beginStep();
    this.updateGap(gas);
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
    if (imp > 0) {
      const tb = Number.isNaN(g.stepBreakdownTime) ? t0 : g.stepBreakdownTime;
      const wasNone = k.stage === 'none';
      k.deposit(imp, gas);
      const span = wasNone ? tEnd - tb : dt;
      if (span > 0) k.step(span, cont / span, sparkActive, gas);
    } else if (k.stage === 'kernel') {
      k.step(dt, cont / dt, sparkActive, gas);
    }
    this.assemble(dt, gas);
  }

  private assemble(dt: number, _gas: IgnitionGasState): void {
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
    s.switchState = c.switchState;
    s.coilEnergy = c.storedEnergy();
    s.energyAtSwitchOff = c.energyAtSwitchOff;
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
    s.breakdownDelay = g.firstBreakdownTime - this.switchOffTime;
    s.sparkDuration = Number.isNaN(g.firstBreakdownTime)
      ? 0
      : (g.conducting || g.awaitingReentry || Number.isNaN(g.extinctionTime) ? c.time : g.extinctionTime) -
        g.firstBreakdownTime;
    s.peakSecondaryVoltage = c.peakSecondaryVoltage;
    // phase
    if (c.switchState === 'closed') s.phase = 'charging';
    else if (!Number.isNaN(g.stepBreakdownTime)) s.phase = 'breakdown';
    else if (g.conducting) s.phase = g.mode === 'arc' ? 'arc' : 'glow';
    else if (g.awaitingReentry) s.phase = 'glow';
    else if (g.breakdownCount > 0) s.phase = 'done';
    else if (this.firedThisCycle && c.storedEnergy() > 1e-6) s.phase = 'breakdown';
    else s.phase = this.firedThisCycle ? 'done' : 'off';
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
    // No breakdown: the first ring-up peak is the largest the coil can produce, so once it
    // is well past without a breakdown the available voltage was below the required one.
    const sinceOff = c.time - this.switchOffTime;
    const tPeak = c.peakSecondaryVoltageTime - this.switchOffTime;
    const noSpark = this.firedThisCycle && g.breakdownCount === 0 && sinceOff > Math.max(0.5e-3, 4 * tPeak);
    s.misfire = k.misfire || noSpark;
  }
}
