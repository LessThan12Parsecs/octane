/**
 * Cycle simulator (public API of src/physics/cycle): ONE physically consistent, energy-conserving
 * quasi-dimensional model of the CFR engine that wires every physics module together.
 *
 * ── Formulation (details in cycle-model.ts and closure.ts) ──────────────────────────────────
 *  - Time is the independent variable; classical RK4 over all continuous states; steps never
 *    cross a discrete event (valve events, ±180°, TDC for instantaneous combustion, dwell/spark,
 *    Wiebe start, cycle wrap, snapshot instants; the Livengood–Wu crossing is located by
 *    redoing the step). Step limits: maxStepDeg (0.25°), fineStepDeg (0.05°) from the spark
 *    to 2 % burned, a Gershgorin bound on the orifice-flow stiffness (h ≤ 2/λ, where RK4's real
 *    stability interval is 2.78 and its amplification stays positive — no chatter),
 *    h ≤ τ_b (burn-up), h ≤ 0.2 m_u/ṁ_b, h ≤ ¼ τ_ab during the end-gas burn-up.
 *  - Open phase: single well-mixed cylinder zone (N, U, burned-gas scalar) between an intake
 *    plenum (fed through the carburettor venturi with the fresh charge) and an exhaust plenum
 *    (outlet orifice to ambient); lash-corrected valve lift × direction-dependent C_D ×
 *    compressible orifice flow; Woschni gas-exchange heat transfer.
 *  - Closed phase: two zones at common p — frozen unburned (U_tot, S_u, m_u conserved-form
 *    states, isentropic except wall heat loss) and burned gas in chemical equilibrium at
 *    (T_b, p), closed by a 2-D Newton on (p, T_b) with the analytic equilibrium Jacobian.
 *  - Combustion: IgnitionSystem (coil, breakdown, arc/glow, Herweg–Maly kernel) → entrainment +
 *    burn-up (Keck 1982) with the K–k turbulence, Taylor-microscale burn-up, S_L tables and the
 *    exact sphere ∩ disc flame geometry; alternatives 'instantaneous-at-tdc', 'wiebe', 'none'.
 *  - Knock: Livengood–Wu from IVC on the end gas ahead of the front (Douaud–Eyzat τ by default,
 *    CFR_KNOCK_DELAY_MODEL; LLNL detailed-chemistry tables, single- or two-stage, as options) →
 *    end-gas burn-up (τ_ab from the history-integrated ignition-time spread of a ΔT-stratified end
 *    gas) and the KnockOscillator acoustic modes (sub-stepped ≤ 2 µs with the exact exponential
 *    release; reported pressure only). Fuels the delay model does not cover (CH4, C3H8, C2H5OH with
 *    the PRF models) have knock disabled (CycleModel.knockAvailable, integral 0).
 *  - NO: rate-controlled extended Zeldovich on the burned equilibrium state (exact step);
 *    frozen at EVO with elements conserved.
 *  - Mechanics: exact slider crank; fixed speed or free crank dynamics with PNH friction.
 *  - Operator splits (first order in the step): spark energy / kernel mass, end-gas burn-up,
 *    NO, burn-out merge, LW integral.
 * Energy: dU_tot = −p dV − Q̇_wall + P_spark exactly (U_tot is a state; every mass transfer between
 * zones is at constant U_tot); the network of volumes conserves mass, species (motored) /
 * elements (fired) and enthalpy to round-off.
 *
 * ── API ─────────────────────────────────────────────────────────────────────────────────────
 *  - {@link EngineSimulator} (SimulatorLike for the worker): constructor(spec, op, options),
 *    setOperatingPoint (rpm/throttle/load/speed mode immediately, the rest — incl. CR — at the
 *    next cycle start θ = −360), reset(), advanceToNextSnapshot(), drainCycleSummaries(), time,
 *    runCycles(n), recordTrace(cycle?), trace, model.
 *  - {@link CycleModel}: the engine model itself (stepUntil, runCycles, recordTrace, conservation
 *    inventories systemMass / elementInventory / speciesInventory / systemEnergy, ledgers in y).
 *  - {@link runClosedCycle}: closed-cycle-only run from a prescribed IVC state (validation).
 *  - {@link CycleModelOptions}: heatTransfer, combustionModel, wiebe, knock, ignitionDelayModel,
 *    knockIntegral, turbulentFlameClosure, calibration parameters (burnRateMultiplier,
 *    taylorScaleMultiplier, turbulenceLengthScaleFactor, turbulenceProduction, woschniMultiplier,
 *    intakePortHeatTransferMultiplier, venturiDischargeCoefficient, intakeRestrictionArea,
 *    dischargeCoefficientMultiplier, knockDecayTime, knockStratificationDT, knockExcitationTime),
 *    numerical controls. Defaults = the ONE global CFR F-1 calibration ({@link CFR_CALIBRATION},
 *    calibration.ts: values, literature ranges, evidence); neutral values give the uncalibrated
 *    model.
 *
 * ── Summary definitions (CycleSummary) ──────────────────────────────────────────────────────
 * IMEP from the thermodynamic pressure (the synthesised knock oscillation is excluded);
 * gross = ∫p dV over −180…180°, net = over the full 720° cycle −360…360°; PMEP = (W_g − W_n)/V_d
 * (positive = pumping loss). Peak pressure (and its angle) = maximum of the cubic Hermite
 * interpolant of the thermodynamic pressure between step ends; maxPressureRiseRate = max over θ
 * of [p(θ + 0.1°) − p(θ)]/0.1° (window position on a 0.01° grid anchored at −360°, plus windows
 * starting at slope discontinuities) — both independent of the integration steps and of the
 * snapshot cadence; the knock ringing is summarised by MAPO. CA10/50/90 from
 * x_b = m_b/(m_u + m_b). Indicated efficiency = W_net/(m_fuel·LHV) with m_fuel = fuel trapped at
 * IVC and the LOWER heating value of the GASEOUS (vaporised) trapped fuel at 298.15 K
 * (thermo/fuels.lowerHeatingValue); ISFC = m_fuel/W_net (NaN if W_net ≤ 0). Trapped mass and
 * residual (burned-gas scalar) fraction at IVC; volumetric efficiency = DRY air inducted through
 * the carburettor over the cycle (net venturi flow × dry-air fraction of the fresh stream) /
 * (ambient dry-air partial density × V_d); fuelMass = m_fuel × x_b(EVO); NO ppm (wet, mole) and
 * CO mole fraction of the whole charge at EVO after the kinetic-NO swap; knock onset angle,
 * end-gas mass fraction at onset = unburned mass AHEAD of the flame front (m − m_e)/m (0 without
 * knock), MAPO = max |p_osc| at the pickup including the maxima between samples (cubic Hermite of
 * the modal sum, sub-steps ≤ 2 µs while the free ringing can still exceed it); misfire = kernel
 * quenched / no breakdown, or x_b(EVO) < 0.1 (a quenched kernel followed by autoignition of the
 * charge still reports misfire: the contract is "the spark kernel failed to develop"; a partial
 * burn 0.1 ≤ x_b(EVO) < 0.5 reports misfire false with CA50 NaN); heatLoss = cylinder wall heat
 * IVC → EVO (J); indicatedWorkGross (J).
 */
import type { EngineSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import { CycleModel, type ClosedCycleInit, type CycleTrace } from './cycle-model';
import type { CycleModelOptions } from './options';

export { EngineSimulator, type SnapshotOptions } from './engine-simulator';
export {
  CycleModel,
  newCycleTrace,
  swapNO,
  resolveDelayModel,
  VENTURI_DISCHARGE_COEFFICIENT,
  OUTLET_DISCHARGE_COEFFICIENT,
  MODE_OPEN,
  MODE_SINGLE,
  MODE_TWO,
  MODE_BURNED,
  type ClosedCycleInit,
  type CycleTrace,
  type CycleProfile,
} from './cycle-model';
export { ZoneClosure } from './closure';
export { CFR_CALIBRATION, CFR_KNOCK_DELAY_MODEL, type CalibratedParameter } from './calibration';
export {
  DEFAULT_CYCLE_MODEL_OPTIONS,
  resolveCycleOptions,
  sanitizeOperatingPoint,
  type CombustionModel,
  type CycleModelOptions,
  type IgnitionDelayModelId,
  type KnockIntegralMode,
  type TurbulentFlameClosure,
  type WiebeOptions,
} from './options';

/** Result of {@link runClosedCycle}. */
export interface ClosedCycleResult {
  /** Full-resolution trace (every integration step, plus the post-event states). */
  trace: CycleTrace;
  /** The model after the run (read ledgers / closure state). */
  model: CycleModel;
  /** Trapped mass, kg; fuel mass, kg. */
  mass: number;
  fuelMass: number;
  /** ∫p dV over the run, J; wall heat loss, J; net spark energy added to the gas, J. */
  work: number;
  heatLoss: number;
  sparkEnergy: number;
  /** Total internal energy at the start and at the end (before the EVO merge), J. */
  U0: number;
  U1: number;
}

/**
 * Closed-cycle-only run (IVC → EVO or any start/end angles) from a prescribed state, with the
 * full closed-cycle physics selected by `options` (heat transfer, combustion model, knock…).
 * No gas exchange; the trace records every step.
 */
export function runClosedCycle(
  spec: EngineSpec,
  op: OperatingPoint,
  init: ClosedCycleInit,
  options: Partial<CycleModelOptions> = {},
): ClosedCycleResult {
  const model = new CycleModel(spec, op, options, init);
  model.recordTrace(0);
  const U0 = model.trace!.U[0];
  const W0 = model.trace!.work[0];
  const Q0 = model.trace!.heatLoss[0];
  const E0 = model.sparkEnergy;
  while (!model.finished) model.stepUntil(Infinity, init.endDeg + 1);
  const tr = model.trace!;
  const n = tr.theta.length;
  // (the EVO merge records a post-event sample with the same U)
  const U1 = tr.U[n - 1];
  return {
    trace: tr,
    model,
    mass: model.mIvc,
    fuelMass: model.fuelMassIvc,
    work: tr.work[n - 1] - W0,
    heatLoss: tr.heatLoss[n - 1] - Q0,
    sparkEnergy: model.sparkEnergy - E0,
    U0,
    U1,
  };
}
