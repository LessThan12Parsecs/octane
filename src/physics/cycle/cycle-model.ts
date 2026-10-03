/**
 * The cycle model: one quasi-dimensional, energy-conserving model of the whole engine (cylinder,
 * intake and exhaust plenums, spark system, flame, knock, NO, heat transfer, mechanics) advanced
 * in TIME with a classical explicit Runge–Kutta (RK4) scheme that never steps across a discrete
 * event. EngineSimulator (engine-simulator.ts) wraps it for the worker; runClosedCycle and
 * CycleModel.runCycles are the validation hooks (index.ts documents the public API).
 *
 * ── Continuous state y (Float64Array) ────────────────────────────────────────────────────────
 *   θ (crank deg), ω (rad/s)                                   mechanics (dθ/dt = ω, dω/dt = α or 0)
 *   cylinder N[NS] (mol), U (J), m_bg (kg)                      open phase (single well-mixed zone)
 *   U_tot (J), S_u (J/K), m_u, m_b, m_e (kg)                    closed phase (two zones at common p)
 *   intake plenum N, U, m_bg;  exhaust plenum N, U, m_bg        always
 *   K, k (J), swirl (kg m²/s)                                   K–k turbulence (Poulos & Heywood 1983)
 *   external-EGR mass scalars m_egr (cylinder, intake plenum) ⊂ m_bg
 *   ledgers: ∫p dV, ∫Q̇_wall (net: cylinder loss − intake-port gain), valve mass flows,
 *   venturi/outlet enthalpy and species flows, gross plenum inflows (warm-up relaxation)
 * m_bg is a passive burned-gas (products: residual + external EGR) mass scalar carried with every
 * stream; m_egr marks the external-EGR part of it (reset at EVO, when everything becomes residual).
 *
 * ── Open phase (EVO → IVC) ───────────────────────────────────────────────────────────────────
 * Fresh charge (humid air + fuel vapour at φ + EGR, thermo/fuels.freshCharge) at the intake
 * mixture temperature and ambient pressure → carburettor venturi (orifice, C_D·A_t scaled by the
 * butterfly throttleArea ratio) → intake plenum → intake valve (ValveFlowModel effective area of
 * the lash-corrected ValveLiftProfile × orificeFlow, bidirectional, direction-dependent C_D) →
 * cylinder → exhaust valve → exhaust plenum → outlet orifice → ambient. Plenums and cylinder are
 * gas-exchange/Plenum volumes (conserved N, U; the stream carries the upstream stagnation
 * enthalpy and composition), so the network conserves mass, species and energy to round-off.
 * Heat: Woschni (gas-exchange constants) over head/valves/piston/liner (wallHeatLoss); intake-port
 * heat transfer into the intake plenum while gas flows through the port (Dittus–Boelter with the
 * port Reynolds number, × intakePortHeatTransferMultiplier; see evaluate).
 *
 * ── Closed phase (IVC → EVO) ─────────────────────────────────────────────────────────────────
 * Before the first burned gas: one zone, T from U_tot (frozen composition). With burned gas: two
 * zones at common p, conserved states U_tot (dU = −p dV − Q̇ + P_spark), S_u (dS_u = s_u dm_u −
 * Q̇_u/T_u), m_u, m_b; closure.ts recovers (p, T_u, T_b) by Newton (burned gas in chemical
 * equilibrium at (T_b, p)). After burn-out: burned zone only (UV equilibrium).
 * Combustion (default 'entrainment'): IgnitionSystem (coil circuit, breakdown, arc/glow, Herweg–
 * Maly kernel) → at hand-off the Keck/Tabaczynski entrainment + burn-up ODEs with u′ from the
 * K–k model (u_T = C_T u′; option: Keck 1982 empirical u_T, ℓ_T), λ = Taylor microscale, S_L from
 * the flame-speed tables (at the TRAPPED fuel blend and φ, see trappedMixture) and A_f from the exact
 * sphere ∩ disc FlameGeometry of the entrained volume V_e = V_b + (m_e − m_b)/ρ_u, mapped onto the
 * disc in proportion A_p h/V (the lumped crevice volume holds charge distributed over the chamber).
 * Heat: Woschni (C₂ term with the motored pressure of the unburned-zone isentrope) split over the
 * zones with the burned wetted areas of the flame geometry (wallHeatLossTwoZone).
 *
 * ── Operator splits (applied after each accepted RK step, first order in the step) ───────────
 *  1. IgnitionSystem.step over the step (it sub-steps the coil internally): the kernel's
 *     burned-mass increment is moved unburned → burned at constant U_tot (S_u −= s_u Δm); then the
 *     electrical energy delivered to the gas minus the kernel's electrode conduction loss is added
 *     to U (closed: U_tot; open: cylinder U) — in the closed phase it pulls the burned zone toward
 *     the ignition model's kernel temperature (never beyond) and the rest goes to the unburned zone
 *     (S_u += ΔE_u/T_u). The hand-off happens at a known instant inside the step (the kernel
 *     integrator locates it exactly): entrainment and burn-up are caught up over the rest of the
 *     step (O(Δt²)). The make-spark diode (IgnitionSystemOptions) is on.
 *  2. End-gas burn-up after autoignition is NOT split: dm_b/dt += (m − m_e)/τ_ab — the end gas
 *     AHEAD of the front; the brush keeps burning up on τ_b — is part of the RK right-hand side
 *     with steps ≤ τ_ab/4 (dm_e/dt += the same, so the brush is unchanged); ∫p dV sees the pressure
 *     rise. The mass it burns per step drives the KnockOscillator, sub-stepped at ≤ knockBurnStep
 *     (2 µs) with the exponential release profile propagated exactly (acoustic pressure is
 *     reported, never fed back to the thermodynamics); mode frequencies follow the sound speed
 *     linearly in time over the step; MAPO is the oscillator's between-sample (Hermite) peak.
 *  3. Burned-zone NO: newly burned gas enters with the unburned NO, then ZeldovichKinetics.
 *     advanceRateControlled (exact) on the equilibrium state at the step end.
 *  4. Burn-out: once m_u/m < burnoutFraction the remainder is merged into the burned zone.
 *  5. Livengood–Wu: the integral is advanced with τ(T_u, p) of the trapped mixture at the end of
 *     the RK step (before the splits; log-mean quadrature of livengood-wu.ts); a step that would
 *     cross I = 1 is redone up to the (interpolated) crossing, so the onset is an exact step
 *     boundary. Armed from IVC while end gas remains AHEAD of the front and the delay model covers
 *     the trapped fuel; the history integral J = ∫(1/τ)(∂lnτ/∂T)T dt gives the burn-up time τ_ab.
 * Everything else (flows, heat, entrainment, end-gas burn-up, turbulence, mechanics) is inside
 * the RK right-hand side. The misfire decision comes from the ignition split and is located to
 * within one (fine, ≤ 0.05°) step. At hand-off the entrained volume covers the chamber-clipped
 * kernel sphere and the brush mass makes Keck's burning speed equal the kernel's S_T,k (onHandoff).
 * Summary metrics independent of the step placement: peak pressure = max of the cubic Hermite
 * interpolant of p (analytic ṗ at step ends); max dp/dθ = max 0.1°-window secant on a 0.01° grid
 * plus windows starting at slope discontinuities; MAPO as above.
 *
 * Numerical guards (documented where applied): pressure-only Newton pre-iterations in the
 * closure; zone heat loss limited to m c_v |T − T_w|/zoneHeatLossMinTime for vanishing zones;
 * burn-out merge at m_u < burnoutFraction·m; knock burn-up ended below 1e-9 of the charge;
 * 'free' speed held ≥ FREE_MODE_MIN_RPM. Warm-up: after each warm-up cycle both plenums are
 * extrapolated to their periodic state (relaxPlenum).
 *
 * Hot path: step() allocates no objects or arrays (the optional trace recorder and the snapshots
 * of engine-simulator.ts do). V8 still boxes doubles passed to / returned from non-inlined calls
 * (≈ 17 MB per 600 rpm cycle, GC ≈ 1 % of the CPU; validation round 1).
 */
import { DEG, RAD2DEG, R_UNIVERSAL } from '../core/constants';
import type { EngineSpec, WallSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import type { CycleSummary } from '../core/snapshot';
import { EL, ELEMENT_COUNTS, FUEL_SPECIES, NE, NS, SP } from '../core/species';
import { EquilibriumSolver } from '../equilibrium';
import {
  autoignitionBurnTime,
  douaudEyzat,
  douaudEyzatLLNL,
  endGasBurnRate,
  KnockOscillator,
  LivengoodWuIntegrator,
  prfLLNLGasoline2011,
  prfLLNLv2,
  residualMoleFraction,
  TabulatedIgnitionDelay,
  twoZoneSoundSpeed,
  ZeldovichKinetics,
  type IgnitionDelayModel,
} from '../chemistry';
import {
  angularMomentumLengthScale,
  DEFAULT_LENGTH_SCALE_FRACTION,
  DEFAULT_TURBULENCE_PARAMS,
  DISSIPATION_COEFFICIENT,
  entrainmentRates,
  FlameGeometry,
  integralLengthScale,
  intakeJetVelocity,
  keckCharacteristicLength,
  keckCharacteristicSpeed,
  keckMeanInletSpeed,
  laminarFlameSpeed,
  lensArea,
  lewisNumbers,
  marksteinLengths,
  meanFlowVelocity,
  newEntrainmentInputs,
  newEntrainmentRates,
  newFlameGeometryResult,
  newTurbulenceInputs,
  newTurbulenceRates,
  taylorMicroscale,
  turbulenceDerivatives,
  turbulenceIntensity,
  type LewisNumbers,
  type MarksteinResult,
  type TurbulenceParams,
  type TurbulenceState,
} from '../combustion';
import {
  gasStateFromTPX,
  newGasState,
  newOrificeFlow,
  orificeFlow,
  Plenum,
  throttleArea,
  ValveFlowModel,
  ValveLiftProfile,
  newValveJet,
  type GasState,
} from '../gas-exchange';
import {
  flatChamberAreas,
  hohenbergCoefficient,
  newChamberAreas,
  WOSCHNI_CONSTANTS,
  newWallHeatResult,
  wallHeatLoss,
  wallHeatLossTwoZone,
  woschniCoefficient,
  type WoschniInputs,
} from '../heat-transfer';
import { IgnitionSystem, type IgnitionCommand, type IgnitionGasState } from '../ignition';
import { CrankTrainDynamics, FrictionTorqueModel, newKinematicState, pnhFmep, SliderCrank } from '../mechanics';
import { CFR_COOLANT_TEMPERATURE } from '../engines/cfr';
import {
  completeCombustionProducts,
  freshCharge,
  fuelFromSelection,
  humidAir,
  lowerHeatingValue,
  type FuelBlend,
} from '../thermo/fuels';
import {
  mixCpMass,
  mixCvMass,
  mixCvMolar,
  mixMolarMass,
  mixSMass,
  mixSMolar,
  mixThermalConductivity,
  mixHMolar,
  mixUMolar,
  mixViscosity,
  temperatureFromUMolar,
} from '../thermo';
import { MOLAR_MASS } from '../thermo/thermo';
import { noSwapDeltas, ZoneClosure } from './closure';
import { resolveCycleOptions, sanitizeOperatingPoint, type CycleModelOptions } from './options';

// =============================================================================================
// State layout
// =============================================================================================

export const I_TH = 0;
export const I_OM = 1;
export const I_CN = 2;
export const I_CU = I_CN + NS;
export const I_CBG = I_CU + 1;
export const I_UT = I_CBG + 1;
export const I_SU = I_UT + 1;
export const I_MU = I_SU + 1;
export const I_MB = I_MU + 1;
export const I_ME = I_MB + 1;
export const I_IN = I_ME + 1;
export const I_IU = I_IN + NS;
export const I_IBG = I_IU + 1;
export const I_EN = I_IBG + 1;
export const I_EU = I_EN + NS;
export const I_EBG = I_EU + 1;
export const I_TK = I_EBG + 1; // mean-flow kinetic energy K
export const I_TKE = I_TK + 1; // turbulent kinetic energy k
export const I_SW = I_TKE + 1;
export const I_W = I_SW + 1; // ∫ p dV
export const I_Q = I_W + 1; // ∫ net wall heat loss dt: cylinder walls − intake-port gain (W·s)
export const I_MIVI = I_Q + 1; // ∫ intake-valve flow into the cylinder
export const I_MIVO = I_MIVI + 1; // ∫ intake-valve backflow out of the cylinder
export const I_MEVO = I_MIVO + 1; // ∫ exhaust-valve flow out of the cylinder
export const I_MEVI = I_MEVO + 1; // ∫ exhaust-valve backflow into the cylinder
export const I_HV = I_MEVI + 1; // ∫ enthalpy flow through the venturi (+ into the intake plenum)
export const I_HO = I_HV + 1; // ∫ enthalpy flow through the outlet (+ out to ambient)
export const I_LV = I_HO + 1; // ∫ species moles through the venturi (+ into the plenum), NS
export const I_LO = I_LV + NS; // ∫ species moles through the outlet (+ out), NS
export const I_MK = I_LO + NS; // ∫ end-gas (autoignition) burn rate dt, kg
// gross INFLOW ledgers of the two plenums (warm-up acceleration, relaxPlenum): enthalpy (J),
// burned-gas scalar (kg), species moles (NS) of everything that entered the intake plenum (venturi
// inflow + intake-valve backflow) and the exhaust plenum (exhaust-valve outflow + outlet backflow)
export const I_GIH = I_MK + 1;
export const I_GIBG = I_GIH + 1;
export const I_GIN = I_GIBG + 1;
export const I_GEH = I_GIN + NS;
export const I_GEBG = I_GEH + 1;
export const I_GEN = I_GEBG + 1;
export const I_GIK = I_GEN + NS; // ∫ intake-port thermal conductance h·A dt (J/K), warm-up relaxation
// external-EGR mass scalars (kg) of the cylinder and the intake plenum: the part of the burned-gas
// scalar m_bg that entered as EGR through the carburettor (the products dilution of the TRAPPED
// charge is m_bg/m; the internal residual is (m_bg − m_egr)/m) — validation round 2
export const I_CEG = I_GIK + 1;
export const I_IEG = I_CEG + 1;
// ∫ −∂ṁ_venturi/∂p_plenum dt (kg/Pa): the intake plenum's pressure-mode relaxation rate (warm-up)
export const I_GIV = I_IEG + 1;
// ∫ gas-to-surface heat flow dt per surface (J): head, piston, liner, intake valves, exhaust valves
// (lumped wall-temperature model, updateWalls)
export const I_QS = I_GIV + 1;
// crevice zone (closed phase): unburned and burned (complete-combustion products) mass at the wall
// temperature and the cylinder pressure, kg
export const I_CRU = I_QS + 5;
export const I_CRB = I_CRU + 1;
export const NY = I_CRB + 1;

/** Cylinder modes. */
export const MODE_OPEN = 0;
export const MODE_SINGLE = 1;
export const MODE_TWO = 2;
export const MODE_BURNED = 3;

// Event kinds
const EV_IVC = 1;
const EV_EVO = 2;
const EV_VALVE = 3; // IVO / EVC (step boundary only)
const EV_BDC_START = 4; // −180
const EV_BDC_END = 5; // +180
const EV_TDC = 6; // 0 (instantaneous combustion)
const EV_IGN = 7; // dwell start / spark (step boundary only)
const EV_WIEBE = 8;
const EV_WRAP = 9; // 360
const EV_END = 10; // closed-cycle-only end
const MAX_EVENTS = 16;

/** ISO 5167-4:2022 classical Venturi tube with an "as cast" convergent section, C = 0.984 (fetched
 * summary of ISO 5167-4 §5.5.4, 2026-09-30) — the round-1 default and the upper end of the
 * calibration range of options.venturiDischargeCoefficient (calibrated CFR value 0.65,
 * calibration.ts): the 9/16 in CFR venturi (throat Re ≈ 1e4, fuel-nozzle bridge) is far below the
 * standard's 2e5 ≤ Re range and not an ISO tube. Modelled as an orifice (no diffuser recovery; the
 * throttle module estimates the unrecovered loss as negligible at CFR flows). */
export const VENTURI_DISCHARGE_COEFFICIENT = 0.984;
/** Exhaust outlet (pipe discharging to the ambient) discharge coefficient. UNVERIFIED: 1 (a plain
 * pipe exit has no vena contracta; its pressure drop at the CFR's ~3 g/s is ≈ 10 Pa anyway). */
export const OUTLET_DISCHARGE_COEFFICIENT = 1.0;
/** Burned-zone temperature below which the kinetic NO is frozen (no rate evaluation), K. The
 * extended-Zeldovich rates are negligible there (round-1 value kept). */
const NO_FREEZE_T = 1000;
/** Relaxation time of the crevice-zone constraint keeper, s (numerical; ≫ the steps, ≪ the cycle). */
const CREVICE_RELAX_TIME = 0.5e-3;
/** Warm-up: extra cycles (at most) while the intake-plenum mass changes by more than WARMUP_MASS_TOL per cycle. */
const WARMUP_EXTRA_MAX = 20;
const WARMUP_MASS_TOL = 1e-3;
/** Burned-mass fraction seeded at the Wiebe start (numerical, see handleEvents EV_WIEBE). */
export const WIEBE_SEED = 1e-5;
/** Lowest crank speed in 'free' mode, rev/min (numerical stall guard; the model needs ω > 0). */
export const FREE_MODE_MIN_RPM = 60;
const FREE_MODE_MIN_OMEGA = (FREE_MODE_MIN_RPM * 2 * Math.PI) / 60;
/** Initial cylinder / exhaust temperatures of a cold start, K — initial conditions only (the
 * warm-up cycles erase them; not physical constants). */
const INITIAL_CYLINDER_T = 900;
const INITIAL_EXHAUST_T = 800;

/** One sample per integration step of a recorded cycle (validation trace). */
export interface CycleTrace {
  theta: number[];
  t: number[];
  volume: number[];
  /** Thermodynamic pressure, Pa. */
  pressure: number[];
  /** Pressure + knock oscillation at the pickup, Pa. */
  pressureReported: number[];
  Tu: number[];
  Tb: number[];
  Tmean: number[];
  mu: number[];
  mb: number[];
  xb: number[];
  me: number[];
  /** Total internal energy of the cylinder charge, J. */
  U: number[];
  /** Cumulative ∫p dV, J (whole run). */
  work: number[];
  /** Cumulative wall heat loss, J (whole run). */
  heatLoss: number[];
  /** Cumulative electrical energy to the gas minus kernel electrode loss, J (whole run). */
  sparkEnergy: number[];
  SL: number[];
  uPrime: number[];
  flameRadius: number[];
  frontArea: number[];
  /** Step-mean burn rate dm_b/dt, kg/s. */
  burnRate: number[];
  heatReleaseRate: number[];
  heatLossRate: number[];
  lwIntegral: number[];
  /** Burned-zone kinetic NO mole fraction. */
  xNO: number[];
  /**
   * Kinetic NO moles carried in the burned-zone thermodynamics by the closure solve of this sample
   * (closure.noKinetic; 0 when options.burnedNOThermo resolves to 'equilibrium').
   */
  nNOClosure: number[];
  mdotIntake: number[];
  mdotExhaust: number[];
  pIntake: number[];
  pExhaust: number[];
  mass: number[];
  mode: number[];
  /** Integral length scale and Taylor microscale, m. */
  L: number[];
  lambda: number[];
}

const TRACE_KEYS: readonly (keyof CycleTrace)[] = [
  'theta', 't', 'volume', 'pressure', 'pressureReported', 'Tu', 'Tb', 'Tmean', 'mu', 'mb', 'xb', 'me', 'U',
  'work', 'heatLoss', 'sparkEnergy', 'SL', 'uPrime', 'flameRadius', 'frontArea', 'burnRate', 'heatReleaseRate',
  'heatLossRate', 'lwIntegral', 'xNO', 'nNOClosure', 'mdotIntake', 'mdotExhaust', 'pIntake', 'pExhaust', 'mass', 'mode', 'L', 'lambda',
];

export function newCycleTrace(): CycleTrace {
  const t = {} as Record<keyof CycleTrace, number[]>;
  for (const k of TRACE_KEYS) t[k] = [];
  return t as CycleTrace;
}

/** Prescribed initial (IVC-like) state of a closed-cycle-only run. */
export interface ClosedCycleInit {
  /** Start angle (the "IVC"), crank deg. */
  startDeg: number;
  /** End angle (the "EVO"), crank deg. */
  endDeg: number;
  /** Charge temperature (K) and pressure (Pa) at startDeg. */
  T: number;
  p: number;
  /** Charge mole fractions (Float64Array(NS), or a species-name map); normalised on use. */
  X: Float64Array | Readonly<Record<string, number>>;
  /** Burned-gas (residual) mass fraction contained in X, for the S_L dilution (default 0). */
  residualMassFraction?: number;
  /**
   * Initial turbulence intensity u′ at startDeg, m/s. Default 0.5 × mean piston speed: the
   * motored-TDC value of Bopp, Vafidis & Whitelaw 1986 (0.45–0.6 S̄p, turbulence.ts); UNVERIFIED as
   * an IVC value (the full-cycle model computes it from the intake jet instead).
   */
  uPrime?: number;
}

/** Per-sub-model CPU time accumulators (ms) when options.profile is set. */
export interface CycleProfile {
  rhsOpen: number;
  rhsClosed: number;
  closure: number;
  flame: number;
  heat: number;
  ignition: number;
  knock: number;
  no: number;
  steps: number;
  rhsEvals: number;
  eqSolves: number;
}

const now = (): number => performance.now();

/**
 * The engine cycle model (see file header). Construct, then drive with {@link stepUntil} /
 * {@link runCycles}; read the evaluated fields (p, Tu, …) for output.
 */
export class CycleModel {
  readonly spec: EngineSpec;
  readonly opts: CycleModelOptions;
  /** Operating point in effect for the current cycle, and the latest requested one. */
  op: OperatingPoint;
  pendingOp: OperatingPoint;
  readonly closedInit: ClosedCycleInit | null;

  // ---- sub-models ----
  readonly kin: SliderCrank;
  readonly ks = newKinematicState();
  readonly dyn: CrankTrainDynamics;
  readonly friction: FrictionTorqueModel;
  readonly ivLift: ValveLiftProfile;
  readonly evLift: ValveLiftProfile;
  readonly ivFlow: ValveFlowModel;
  readonly evFlow: ValveFlowModel;
  readonly flameGeom: FlameGeometry;
  readonly closure = new ZoneClosure();
  readonly eqAux = new EquilibriumSolver();
  ign: IgnitionSystem | null = null;
  readonly knockOsc: KnockOscillator;
  readonly lw = new LivengoodWuIntegrator();
  readonly zeld = new ZeldovichKinetics();
  delayModel: IgnitionDelayModel;
  readonly cyl: Plenum;
  readonly intake: Plenum;
  readonly exhaust: Plenum;
  readonly fresh: GasState = newGasState();
  readonly ambient: GasState = newGasState();
  fuel: FuelBlend;
  private readonly turbParams: TurbulenceParams;
  private readonly lengthFraction: number;

  // ---- time / state ----
  readonly y = new Float64Array(NY);
  t = 0;
  cycle = 0;
  mode = MODE_OPEN;
  /** True once a closed-cycle-only run reached its end angle. */
  finished = false;
  /** Set when a step wrapped the cycle (θ = 360 → −360). */
  wrapped = false;

  // RK scratch
  private readonly dy = new Float64Array(NY); // derivative at (t, y) (valid when derivValid)
  private readonly k2 = new Float64Array(NY);
  private readonly k3 = new Float64Array(NY);
  private readonly k4 = new Float64Array(NY);
  private readonly yt = new Float64Array(NY);
  private readonly y0 = new Float64Array(NY);
  private readonly dy0 = new Float64Array(NY);
  private readonly scratch = new Float64Array(NY);
  /** Kahan compensation of the RK state update (and its copy at the step start). */
  private readonly yComp = new Float64Array(NY);
  private readonly yComp0 = new Float64Array(NY);
  private derivValid = false;

  // ---- events of the current cycle ----
  private readonly evAngle = new Float64Array(MAX_EVENTS);
  private readonly evKind = new Int32Array(MAX_EVENTS);
  private nEv = 0;
  private evNext = 0;
  readonly ignCmd: IgnitionCommand = { dwellStartDeg: -360, sparkDeg: 0 };
  ivcDeg = -152;
  evoDeg = 141;

  // ---- per-cycle composition data ----
  readonly Xfresh = new Float64Array(NS);
  readonly airX = new Float64Array(NS);
  /** Residual mass fraction of the trapped charge (burned-gas scalar) and products dilution for S_L. */
  yRes = 0;
  /** External-EGR mass fraction of the trapped charge (part of yRes). */
  yEgr = 0;
  xDil = 0;
  xResMole = 0;
  private cdaVenturi = 0;
  /** Fresh-charge mass per unit mass of humid air through the venturi, 1/(1 − Y_fuel − Y_EGR). */
  private freshPerAir = 1;
  /** Intake-port heat transfer: mult·(k/D)·πDL (W/K per unit Nu), 4/(πDμ), Pr^0.4, wall T. */
  private portHA = 0;
  private portReCoef = 0;
  private portPr04 = 1;
  private portTw = 0;
  /** Intake-port heat into the intake plenum at the last evaluation, W. */
  Qport = 0;
  /** Keck-1982 closure: selected; inlet density ρ_i (intake plenum at IVC), mean inlet speed ū_i. */
  private readonly keck: boolean;
  rhoInlet = 0;
  keckInletSpeed = 0;
  /** Two-stage knock integral: selected (and supported by the delay model); current stage (1, 2). */
  private twoStage = false;
  lwStage = 1;
  private cdaOutlet = 0;
  private throttleFor = NaN;
  walls: WallSpec;
  private wallTavg = 450;
  /** Crevice zone: active, volume (m³), temperature (K) and per-kg properties at it (set at IVC). */
  readonly creviceOn: boolean;
  readonly Vcr: number;
  Tcr = 450;
  private uuCr = 0;
  private huCr = 0;
  private uccCr = 0;
  private hccCr = 0;
  private Rcc = 0;
  /** Crevice net inflow (kg/s, + into the crevice) and its heat to the walls (W) at the last evaluation. */
  mdotCr = 0;
  Qcr = 0;
  /** Lumped wall model active (spec.walls.thermalResistance and options.wallTemperatureModel 'lumped'). */
  private readonly wallsLumped: boolean;
  /** Wall excess over the coolant, K: head, piston, liner, intake valve, exhaust valve. */
  readonly wallExcess = new Float64Array(5);
  /** Surface heat ledgers and time at the start of the current cycle. */
  private readonly qsCycleStart = new Float64Array(5);
  private tCycleStart = 0;
  private pCrankcase = 101325;
  // viscosity table of the frozen unburned mixture (T grid)
  private readonly muTab = new Float64Array(100);
  private readonly muT0 = 200;
  private readonly muDT = 25;

  // ---- closed-cycle bookkeeping ----
  mIvc = 0;
  pIvc = 0;
  TIvc = 0;
  VIvc = 0;
  fuelMassIvc = 0;
  /** Burned-zone NO moles (kinetic). */
  nNO = 0;
  /** Flame/ignition status. */
  sparkFired = false;
  kernelMassPrev = 0;
  handedOff = false;
  flameActive = false;
  burnDone = false;
  misfire = false;
  kernelQuenched = false;
  knockOnset = false;
  knockBurning = false;
  /** End-gas autoignition burn-up time τ_ab of this cycle, s (NaN before onset). */
  tauAB = NaN;
  private lwTauPrev = NaN;
  /** Livengood–Wu history integral J = ∫(1/τ)(∂lnτ/∂T)T dt since IVC (see onAutoignition). */
  lwJ = 0;
  private lwRPrev = NaN;
  private lwGPrev = 0;
  /**
   * True when the end-gas ignition-delay model covers the trapped fuel (checked at IVC). The PRF
   * delay models have no data for CH4 / C3H8 / C2H5OH (τ = NaN, Douaud–Eyzat throws): the knock
   * integral is then disarmed for the cycle (knock cannot be predicted for that fuel; reported
   * integral 0) instead of propagating NaN into the snapshots.
   */
  knockAvailable = true;
  /**
   * Mixture descriptors of the TRAPPED charge (set at IVC from its element and fuel-species
   * content): fuel blend and fresh-charge equivalence ratio. During operating-point transitions the
   * intake plenum still delivers the previous mixture for several cycles; S_L, Markstein/Lewis
   * numbers, the ignition delay, the residual mole fraction and the LHV use these, not the
   * requested op.fuel / op.equivalenceRatio (validation round 1).
   */
  fuelTrapped: FuelBlend;
  phiTrapped = 1;
  /** Kinetic NO carried in the burned-zone thermodynamics (options.burnedNOThermo; closure.noCoupled). */
  readonly noCoupled: boolean;
  /** Failed closure solves (after the robust retry) since construction — diagnostic. */
  closureFailures = 0;
  knockOnsetDeg = NaN;
  knockEndGasFraction = 0;
  mapo = 0;
  knockQdot = 0;
  tKnockOnset = NaN;
  /** End-gas mass at the knock onset and the autoignition-burn ledger I_MK there, kg (source shells). */
  private knockEgOnset = 0;
  private knockMkOnset = 0;
  /** Knock-mode sound speed, γ and V at the end of the previous knock step (0: none yet). */
  private knockC0 = 0;
  private knockG0 = 0;
  private knockV0 = 0;
  private knockGamma = 1.3;
  private Lref = 0;
  private rhoRef = 0;
  private wiebeStart = 0;
  private sparkEnergyPrev = 0;
  private electrodeLossPrev = 0;
  /** Cumulative net spark energy added to the gas (to gas − kernel electrode loss), J. */
  sparkEnergy = 0;
  /** Cumulative kernel electrode loss, J. */
  electrodeLoss = 0;
  private suLast = 0;
  private TmotGuess = 400;

  // ---- evaluated quantities (last RHS evaluation; the accepted step end after each step) ----
  V = 0;
  Vdot = 0;
  h = 0;
  p = 0;
  T = 0;
  Tu = 0;
  Tb = 0;
  mCyl = 0;
  Qwall = 0;
  Qu = 0;
  Qb = 0;
  mdotIv = 0;
  mdotEv = 0;
  mdotV = 0;
  mdotO = 0;
  pInt = 0;
  pExh = 0;
  SL = 0;
  uPrime = 0;
  L = 0;
  lambda = 0;
  rf = 0;
  Af = 0;
  rb = 0;
  /** Disc share of the cylinder volume A_p h/V (the rest is the lumped crevice volume). */
  discScale = 1;
  mdotB = 0;
  mdotE = 0;
  burnSpeed = 0;
  tauB = Infinity;
  rhoU = 0;
  nuU = 1.6e-5;
  dpdt = 0;
  alpha = 0;
  stiffness = 0;
  pMot = 0;
  hcoef = 0;
  /** Step-mean burn rate (kg/s) and chemical heat-release rate (W) of the last step. */
  burnRateStep = 0;
  hrr = 0;
  /** Guesses. */
  private pGuess = 1e5;
  private TbGuess = 2400;
  private TuGuess = 400;
  private tEval = 0;
  private dpEval = 0;
  private dTbEval = 0;
  private rfGuess = 0;
  private rbGuess = 0;

  // ---- per-cycle summary accumulators ----
  private readonly summaries: CycleSummary[] = [];
  private wCycleStart = 0;
  private wBdcStart = 0;
  private wBdcEnd = 0;
  private qIvc = 0;
  private qEvo = 0;
  peakP = 0;
  peakPDeg = 0;
  private maxDp = 0;
  /** Pressure interpolation state at the start of the current step (post-split). */
  private hTheta = NaN;
  private hP = NaN;
  private hDp = 0;
  private hMode = MODE_OPEN;
  /** Fine dp/dθ grid: index of the last point (from −360°), filled count, ring buffer of p. */
  private dpGridK = -1;
  private dpGridN = 0;
  private readonly dpBuf = new Float64Array(DP_SUBDIV);
  /** Pending dp/dθ windows starting at slope discontinuities (angle, pressure). */
  private readonly dpKinkTh = new Float64Array(8);
  private readonly dpKinkP = new Float64Array(8);
  private dpKinkN = 0;
  ca10 = NaN;
  ca50 = NaN;
  ca90 = NaN;
  private volEff = 0;
  private mVentCycleStart = 0;
  private noPpm = 0;
  private coFrac = 0;
  private xbEvo = 0;
  private closedHappened = false;

  // ---- trace ----
  traceCycle = -1;
  trace: CycleTrace | null = null;
  // ---- profiling ----
  readonly prof: CycleProfile = {
    rhsOpen: 0, rhsClosed: 0, closure: 0, flame: 0, heat: 0, ignition: 0, knock: 0, no: 0, steps: 0, rhsEvals: 0, eqSolves: 0,
  };

  // ---- scratch objects ----
  private readonly ofIv = newOrificeFlow();
  private readonly ofEv = newOrificeFlow();
  private readonly ofV = newOrificeFlow();
  private readonly ofO = newOrificeFlow();
  private readonly jet = newValveJet();
  private readonly areas = newChamberAreas();
  private readonly heat = newWallHeatResult();
  private readonly fg = newFlameGeometryResult();
  private readonly fgB = newFlameGeometryResult();
  private readonly fgR = newFlameGeometryResult();
  private readonly woschni: WoschniInputs = {
    bore: 0, pressure: 0, temperature: 0, meanPistonSpeed: 0, phase: 'gas-exchange',
    motoredPressure: 0, displacedVolume: 0, refPressure: 0, refTemperature: 0, refVolume: 0, variant: 'woschni1967',
  };
  private readonly turbIn = newTurbulenceInputs();
  private readonly turbOut = newTurbulenceRates();
  private readonly turbState: TurbulenceState = { K: 0, k: 0, swirl: 0 };
  private readonly entIn = newEntrainmentInputs();
  private readonly entOut = newEntrainmentRates();
  private readonly gas: IgnitionGasState = {
    p: 1e5, Tu: 300, rhoU: 1, X: new Float64Array(NS), SL: 0, marksteinLength: 0, flameThickness: 0,
    uPrime: 0, integralScale: 1e-3, dissipationLength: 1e-3, expansionRatio: 1, lewisNumber: 1, kinematicViscosity: 1.6e-5, flowVelocity: 0,
  };
  private readonly mk: MarksteinResult = { burned: 0, unburned: 0, numberBurned: 0, lf: 0, sigma: 0, beta: 0, leEff: 0, SL: 0 };
  private readonly le: LewisNumbers = { fuel: 0, oxygen: 0, deficient: 0, alphaU: 0 };
  private readonly Nscr = new Float64Array(NS);
  private readonly Xscr = new Float64Array(NS);
  private readonly Nb = new Float64Array(NS);
  private readonly bScr = new Float64Array(NE);

  constructor(spec: EngineSpec, op: OperatingPoint, options: Partial<CycleModelOptions> = {}, closedInit: ClosedCycleInit | null = null) {
    this.spec = spec;
    this.opts = resolveCycleOptions(options);
    this.closedInit = closedInit;
    this.pendingOp = sanitizeOperatingPoint(spec, op);
    this.op = cloneOp(this.pendingOp);
    this.fuel = fuelFromSelection(this.op.fuel);
    this.fuelTrapped = this.fuel;
    this.phiTrapped = this.op.equivalenceRatio;
    this.delayModel = resolveDelayModel(this.opts.ignitionDelayModel);
    const nt = this.opts.burnedNOThermo;
    this.noCoupled = nt === 'kinetic' || (nt === 'auto' && this.opts.combustionModel === 'entrainment');
    this.keck = this.opts.turbulentFlameClosure === 'keck1982';
    this.turbParams = { ...DEFAULT_TURBULENCE_PARAMS, cBeta: this.opts.turbulenceProduction };
    this.lengthFraction = DEFAULT_LENGTH_SCALE_FRACTION * this.opts.turbulenceLengthScaleFactor;
    const g = spec.geometry;
    this.kin = new SliderCrank(g, this.op.compressionRatio);
    this.dyn = new CrankTrainDynamics(this.kin, spec.masses);
    this.friction = new FrictionTorqueModel(this.kin);
    const lash = this.opts.valveLash;
    this.ivLift = new ValveLiftProfile({ ...spec.intakeValve, lash });
    this.evLift = new ValveLiftProfile({ ...spec.exhaustValve, lash });
    this.ivFlow = new ValveFlowModel(spec.intakeValve, 'intake', { portDiameter: Math.max(spec.manifolds.intakePortDiameter, spec.intakeValve.stemDiameter * 1.01) });
    this.evFlow = new ValveFlowModel(spec.exhaustValve, 'exhaust', { portDiameter: Math.max(spec.manifolds.exhaustPortDiameter, spec.exhaustValve.stemDiameter * 1.01) });
    // Flame-geometry table must cover the tallest chamber: stroke + h_TDC at the lowest CR.
    const hMax = this.kin.pistonTravel + this.kin.clearanceHeightTDCForCR(g.compressionRatioRange[0]);
    this.flameGeom = new FlameGeometry(g.bore, spec.sparkPlug.gapCenter, { maxHeight: hMax * 1.02 });
    this.knockOsc = new KnockOscillator(g.bore, { decayTime: this.opts.knockDecayTime, sensor: this.opts.knockSensor, band: this.opts.mapoBand ?? undefined });
    const X0 = new Float64Array(NS);
    X0[SP.N2] = 1;
    this.cyl = new Plenum(this.kin.volume(-2 * Math.PI), 300, 1e5, X0);
    this.intake = new Plenum(spec.manifolds.intakeVolume, 300, 1e5, X0);
    this.exhaust = new Plenum(spec.manifolds.exhaustVolume, 300, 1e5, X0);
    this.walls = { ...spec.walls };
    this.wallsLumped = !!spec.walls.thermalResistance && this.opts.wallTemperatureModel === 'lumped';
    this.Vcr = spec.geometry.quenchCreviceVolume ?? 0;
    // (the crevice zone follows the pressure continuously; the idealised prescribed burns — Wiebe, the
    // instantaneous fuel–air-cycle burn — keep the round-1 single-volume definition)
    const cmod = this.opts.combustionModel;
    this.creviceOn = this.opts.creviceModel && this.opts.heatTransfer && this.Vcr > 0 && (cmod === 'entrainment' || cmod === 'none');
    this.initState();
  }

  // ===========================================================================================
  // Public control
  // ===========================================================================================

  /** Change operating conditions: rpm / throttle / load / speed mode now, the rest at the next cycle start. */
  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    // undefined keys of the patch keep the pending value (a spread would overwrite it with undefined);
    // non-finite values are repaired from the pending point by the sanitiser (validation round 2)
    const next = { ...this.pendingOp };
    for (const k of Object.keys(patch) as (keyof OperatingPoint)[]) {
      const v = patch[k];
      if (v !== undefined) (next as Record<string, unknown>)[k] = v;
    }
    if (patch.fuel) next.fuel = { ...patch.fuel };
    this.pendingOp = sanitizeOperatingPoint(this.spec, next, this.pendingOp);
    const o = this.op;
    o.rpm = this.pendingOp.rpm;
    o.throttle = this.pendingOp.throttle;
    o.loadTorque = this.pendingOp.loadTorque;
    o.speedMode = this.pendingOp.speedMode;
    if (o.speedMode === 'fixed') this.y[I_OM] = (o.rpm * 2 * Math.PI) / 60;
    this.updateFriction();
    this.derivValid = false;
  }

  /** Restart from t = 0 with the latest operating point (runs the warm-up cycles again). */
  reset(): void {
    this.op = cloneOp(this.pendingOp);
    this.initState();
  }

  /** Completed-cycle summaries since the last drain (oldest first). */
  drainSummaries(): CycleSummary[] {
    return this.summaries.splice(0, this.summaries.length);
  }

  /**
   * Knock-integral progress for output (autoignition at 1): the Livengood–Wu integral, or for the
   * two-stage integral (I₁)/2 during the first stage and (1 + I₂)/2 after it; 0 when the knock
   * model is unavailable for the fuel.
   */
  get knockProgress(): number {
    const I = this.lw.integral;
    if (!Number.isFinite(I)) return this.knockOnset ? 1 : 0;
    if (!this.twoStage) return I;
    return this.lwStage === 1 ? 0.5 * I : 0.5 + 0.5 * Math.min(I, 1);
  }

  /** Current crank angle, deg. */
  get theta(): number {
    return this.y[I_TH];
  }

  /** Current speed, rev/min. */
  get rpm(): number {
    return (this.y[I_OM] * 60) / (2 * Math.PI);
  }

  /** Burned-mass fraction m_b/m of the closed charge; after EVO (same cycle) the value at EVO. */
  get xb(): number {
    if (this.mode === MODE_OPEN) return this.xbEvo;
    const m = this.y[I_MU] + this.y[I_MB];
    return m > 0 ? this.y[I_MB] / m : 0;
  }

  /**
   * Run n complete cycles (to n more cycle boundaries) without producing snapshots and return
   * their summaries (they are removed from the drain queue).
   */
  runCycles(n: number): CycleSummary[] {
    const before = this.summaries.length;
    while (this.summaries.length - before < n && !this.finished) this.stepUntil(Infinity, 360);
    return this.summaries.splice(before, n);
  }

  /** Record a full-resolution trace of cycle `cycle` (default: the next complete cycle). */
  recordTrace(cycle?: number): void {
    this.traceCycle = cycle ?? (this.y[I_TH] <= -360 + 1e-9 ? this.cycle : this.cycle + 1);
    this.trace = newCycleTrace();
    if (this.traceCycle === this.cycle) {
      if (!this.derivValid) this.refresh();
      this.recordSample(this.p + this.knockOscillation());
    }
  }

  /**
   * Advance with steps until t ≥ tTarget or θ reaches thetaTarget (≤ 360, current cycle), or the
   * cycle wraps, or a closed-cycle-only run finishes.
   */
  stepUntil(tTarget: number, thetaTarget: number): void {
    this.wrapped = false;
    this.interrupt = false;
    let guard = 0;
    while (!this.finished && guard++ < 10_000_000) {
      if (this.t >= tTarget - 1e-14 * Math.max(1, Math.abs(tTarget))) return;
      if (this.y[I_TH] >= thetaTarget - 1e-9) return;
      this.step(tTarget, thetaTarget);
      if (this.wrapped || this.interrupt) return;
    }
  }

  // ===========================================================================================
  // Initialisation
  // ===========================================================================================

  private initState(): void {
    const y = this.y;
    y.fill(0);
    this.yComp.fill(0);
    this.t = 0;
    this.cycle = 0;
    this.finished = false;
    this.summaries.length = 0;
    this.closure.eq.reset();
    this.eqAux.reset();
    this.knockOsc.reset();
    this.lw.reset();
    this.ign = null;
    this.sparkEnergy = 0;
    this.electrodeLoss = 0;
    this.pGuess = 1e5;
    this.TbGuess = 2400;
    this.TuGuess = 400;
    this.tEval = 0;
    this.dpEval = 0;
    this.dTbEval = 0;
    this.rfGuess = 0;
    this.rbGuess = 0;
    this.TmotGuess = 400;
    this.stiffness = 0;
    this.alpha = 0;
    this.T = 0;
    this.Tu = 0;
    this.Tb = 0;
    this.pInt = 0;
    this.nuU = 1.6e-5;
    this.tauB = Infinity;
    this.mdotB = 0;
    this.lwForce = false;
    this.derivValid = false;
    this.muTab.fill(0);
    this.fuel = fuelFromSelection(this.op.fuel);
    this.kin.setCompressionRatio(this.op.compressionRatio);
    y[I_OM] = (this.op.rpm * 2 * Math.PI) / 60;
    {
      const w = this.spec.walls;
      const Tc = w.referenceCoolantTemperature ?? CFR_COOLANT_TEMPERATURE;
      this.wallExcess[0] = w.headTemperature - Tc;
      this.wallExcess[1] = w.pistonTemperature - Tc;
      this.wallExcess[2] = w.linerTemperature - Tc;
      this.wallExcess[3] = w.intakeValveTemperature - Tc;
      this.wallExcess[4] = w.exhaustValveTemperature - Tc;
      this.qsCycleStart.fill(0);
      this.tCycleStart = 0;
    }
    this.prepareCycleData();
    const op = this.op;
    if (this.closedInit) {
      this.initClosedOnly(this.closedInit);
      return;
    }
    // Cold start at θ = −360 (gas-exchange TDC): cylinder full of complete-combustion products at
    // ambient pressure, plenums at ambient pressure (fresh charge / products).
    y[I_TH] = -360;
    // (a motored engine never contains products: start it from fresh charge everywhere)
    const prod = this.opts.combustionModel === 'none' ? this.Xfresh : completeCombustionProducts(this.Xfresh);
    this.cyl.volume = this.kin.volume(-2 * Math.PI);
    this.cyl.setTPX(INITIAL_CYLINDER_T, op.ambientPressure, prod);
    this.intake.setTPX(op.intakeMixtureTemperature, op.ambientPressure, this.Xfresh);
    this.exhaust.setTPX(INITIAL_EXHAUST_T, op.ambientPressure, prod);
    for (let k = 0; k < NS; k++) {
      y[I_CN + k] = this.cyl.N[k];
      y[I_IN + k] = this.intake.N[k];
      y[I_EN + k] = this.exhaust.N[k];
    }
    y[I_CU] = this.cyl.U;
    y[I_IU] = this.intake.U;
    y[I_EU] = this.exhaust.U;
    const motored = this.opts.combustionModel === 'none';
    y[I_CBG] = motored ? 0 : this.cyl.mass();
    // the fresh charge in the intake plenum carries its external EGR (products) from the start
    y[I_IBG] = this.intake.mass() * this.op.egrFraction;
    y[I_IEG] = y[I_IBG];
    y[I_CEG] = 0;
    y[I_EBG] = motored ? 0 : this.exhaust.mass();
    this.mode = MODE_OPEN;
    this.resetCycleFlags();
    this.startCycle();
    this.refresh();
    // warm-up: silent cycles; after each, both plenums (the slow modes of the network) are
    // extrapolated to their periodic state (relaxPlenum).
    // (the warm-up is extended — up to WARMUP_EXTRA_MAX cycles — while the intake plenum's mass still
    // changes by more than WARMUP_MASS_TOL per cycle: throttled operation, validation round 2)
    let target = this.opts.warmupCycles;
    const L0 = new Float64Array(NY);
    L0.set(y);
    this.intakeMassChange = 0;
    this.warmingUp = true;
    while (this.cycle < target) {
      this.stepUntil(Infinity, 360); // (returns early on interrupts: spark command, knock onset)
      if (this.wrapped) {
        this.relaxPlenum(this.intake, L0, I_GIN, I_IN, I_IU, I_IBG, y[I_GIK] - L0[I_GIK], I_IEG);
        this.relaxIntakePressure(L0);
        this.relaxPlenum(this.exhaust, L0, I_GEN, I_EN, I_EU, I_EBG, 0, -1);
        L0.set(y);
        // (at most WARMUP_EXTRA_MAX extra cycles, and not more than twice the requested number)
        const extraMax = Math.min(WARMUP_EXTRA_MAX, 2 * this.opts.warmupCycles);
        if (this.cycle >= target && target > 0 && this.intakeMassChange > WARMUP_MASS_TOL && target < this.opts.warmupCycles + extraMax) target++;
      }
    }
    this.warmingUp = false;
    this.t = 0;
    this.cycle = 0;
    this.summaries.length = 0;
    this.tCycleStart = 0;
    this.derivValid = false;
    this.refresh();
  }

  /**
   * Warm-up acceleration (initialisation only — never while cycles are emitted): extrapolate a
   * plenum to its periodic state after one warm-up cycle. A well-mixed plenum of mass m_p through
   * which the mass m_in flows per cycle approaches its periodic state as a dilution mode,
   * x_{k+1} − x* = λ (x_k − x*) with λ = exp(−m_in/m_p) (residence-time estimate: 0.81 for the 10 L
   * exhaust plenum, 0.5 for the 1 L intake plenum at RON — as observed cycle by cycle in round 1),
   * so from its temperature and mass fractions at the start and end of the cycle
   *   x* = x_end + (x_end − x_start)·λ/(1 − λ)
   * (Aitken extrapolation with the physical λ; for the temperature the wall heat exchange of the
   * intake port adds ∫hA dt / c_p to the exchanged mass); set at the plenum's current pressure. At the
   * periodic state x_end = x_start and nothing changes. Round 1 relied on 3 plain warm-up cycles:
   * the adiabatic exhaust plenum went 800 K → 1273 K over ≈ 20 cycles and the first emitted RON
   * cycle was 0.5 % low in IMEP, 3.7 % in peak pressure and 29 % in MAPO.
   */
  private relaxPlenum(pl: Plenum, L0: Float64Array, iGin: number, oN: number, oU: number, oBg: number, wallConductance: number, oEg: number): void {
    const y = this.y;
    // mass that entered the plenum during the cycle (gross inflow ledger)
    let mIn = 0;
    for (let k = 0; k < NS; k++) mIn += (y[iGin + k] - L0[iGin + k]) * MOLAR_MASS[k];
    const N0 = this.Nscr;
    let m0 = 0;
    for (let k = 0; k < NS; k++) {
      N0[k] = L0[oN + k];
      m0 += N0[k] * MOLAR_MASS[k];
    }
    const T0 = temperatureFromUMolar(N0, L0[oU], pl.state.T);
    let m1 = 0;
    for (let k = 0; k < NS; k++) m1 += y[oN + k] * MOLAR_MASS[k];
    const T1 = temperatureFromUMolar(y.subarray(oN, oN + NS), y[oU], pl.state.T);
    if (!(mIn > 0 && m0 > 0 && m1 > 0 && T0 > 0 && T1 > 0)) return;
    const lam = Math.exp(-mIn / m1);
    const f = Math.min(lam / (1 - lam), 20);
    // temperature: dilution plus the heat exchange with the port walls (∫hA dt, J/K, over the cycle)
    const lamT = Math.exp(-(mIn + wallConductance / pl.state.cp) / m1);
    const fT = Math.min(lamT / (1 - lamT), 20);
    const T = T1 + (T1 - T0) * fT;
    if (!(T > 200 && T < 3000)) return;
    // mass fractions (as moles per kg) extrapolated the same way
    const X = this.Xscr;
    let n = 0;
    for (let k = 0; k < NS; k++) {
      const a = N0[k] / m0;
      const b = y[oN + k] / m1;
      const v = b + (b - a) * f;
      X[k] = v > 0 ? v : 0;
      n += X[k];
    }
    if (!(n > 0)) return;
    for (let k = 0; k < NS; k++) X[k] /= n;
    const yb0 = L0[oBg] / m0;
    const yb1 = y[oBg] / m1;
    const yb = Math.min(1, Math.max(0, yb1 + (yb1 - yb0) * f));
    const ye0 = oEg >= 0 ? L0[oEg] / m0 : 0;
    const ye1 = oEg >= 0 ? y[oEg] / m1 : 0;
    const ye = Math.min(yb, Math.max(0, ye1 + (ye1 - ye0) * f));
    pl.setTPX(T, pl.state.p, X);
    for (let k = 0; k < NS; k++) y[oN + k] = pl.N[k];
    y[oU] = pl.U;
    y[oBg] = yb * pl.mass();
    if (oEg >= 0) y[oEg] = ye * pl.mass();
    this.derivValid = false;
    this.refresh();
  }

  /**
   * Warm-up acceleration of the intake plenum's MASS (pressure) mode, after relaxPlenum: the plenum
   * mass obeys dm/dt = ṁ_venturi(p) − ṁ_engine(p) with p = mRT/V, a relaxation of rate
   *   k = (RT/V)(|∂ṁ_venturi/∂p| + ṁ_engine/p),
   * the engine flow ∝ p (trapped mass ∝ intake density). Over one warm-up cycle of duration T_c the
   * mode decays by λ = e^{−k T_c}, so the periodic mass is m* = m₁ + (m₁ − m₀) λ/(1 − λ) (m₀, m₁ at
   * the start and end of the cycle; |∂ṁ_v/∂p| and ṁ_engine are cycle means from the ledgers). Wide
   * open (the CFR): λ ≈ 0, nothing changes; throttled the mode is slow (λ ≈ 0.9 at throttle 0.1),
   * and round 1's first emitted cycle was 26–30 % off in IMEP at throttle ≤ 0.1 (validation round 2).
   */
  /** Relative change of the intake-plenum mass over the last warm-up cycle (before extrapolation). */
  private intakeMassChange = 0;
  /** True during the silent warm-up cycles of initState. */
  private warmingUp = false;

  private relaxIntakePressure(L0: Float64Array): void {
    const y = this.y;
    const pl = this.intake;
    let m0 = 0;
    let m1 = 0;
    for (let k = 0; k < NS; k++) {
      m0 += L0[I_IN + k] * MOLAR_MASS[k];
      m1 += y[I_IN + k] * MOLAR_MASS[k];
    }
    this.intakeMassChange = m0 > 0 ? Math.abs(m1 / m0 - 1) : 0;
    const Tc = (720 / RAD2DEG) / Math.max(y[I_OM], 1e-9); // (fixed-speed warm-up: one cycle)
    const dmdp = (y[I_GIV] - L0[I_GIV]) / Tc;
    const mEng = (y[I_MIVI] - L0[I_MIVI] - (y[I_MIVO] - L0[I_MIVO])) / Tc;
    pl.updateState();
    const st = pl.state;
    if (!(m0 > 0 && m1 > 0 && st.p > 0 && dmdp >= 0 && mEng >= 0)) return;
    const k = ((st.R * st.T) / pl.volume) * (dmdp + mEng / st.p);
    const lam = Math.exp(-k * Tc);
    const f = Math.min(lam / (1 - lam), 50);
    const mStar = m1 + (m1 - m0) * f;
    if (!(mStar > 0.2 * m1 && mStar < 5 * m1) || Math.abs(mStar / m1 - 1) < 1e-9) return;
    const X = this.Xscr;
    const n = st.p * pl.volume / (R_UNIVERSAL * st.T);
    for (let k2 = 0; k2 < NS; k2++) X[k2] = y[I_IN + k2];
    let nt = 0;
    for (let k2 = 0; k2 < NS; k2++) nt += X[k2];
    for (let k2 = 0; k2 < NS; k2++) X[k2] /= nt;
    const scale = mStar / m1;
    const yb = y[I_IBG] / m1;
    const ye = y[I_IEG] / m1;
    pl.setTPX(st.T, st.p * scale, X);
    for (let k2 = 0; k2 < NS; k2++) y[I_IN + k2] = pl.N[k2];
    y[I_IU] = pl.U;
    y[I_IBG] = yb * pl.mass();
    y[I_IEG] = ye * pl.mass();
    void n;
    this.derivValid = false;
    this.refresh();
  }

  private initClosedOnly(ci: ClosedCycleInit): void {
    const y = this.y;
    y[I_TH] = ci.startDeg;
    const X = this.Xscr;
    X.fill(0);
    if (ci.X instanceof Float64Array) {
      for (let k = 0; k < NS; k++) X[k] = ci.X[k];
    } else {
      for (const [name, v] of Object.entries(ci.X)) {
        const k = (SP as Record<string, number>)[name];
        if (k === undefined) throw new RangeError(`ClosedCycleInit: unknown species ${name}`);
        X[k] = v;
      }
    }
    const V = this.kin.volume(ci.startDeg * DEG);
    // plenums: irrelevant, keep at ambient
    this.intake.setTPX(this.op.intakeMixtureTemperature, this.op.ambientPressure, this.Xfresh);
    this.exhaust.setTPX(INITIAL_EXHAUST_T, this.op.ambientPressure, this.Xfresh);
    for (let k = 0; k < NS; k++) {
      y[I_IN + k] = this.intake.N[k];
      y[I_EN + k] = this.exhaust.N[k];
    }
    y[I_IU] = this.intake.U;
    y[I_EU] = this.exhaust.U;
    this.cyl.volume = V;
    this.cyl.setTPX(ci.T, ci.p, X);
    for (let k = 0; k < NS; k++) y[I_CN + k] = this.cyl.N[k];
    y[I_CU] = this.cyl.U;
    y[I_CBG] = this.cyl.mass() * (ci.residualMassFraction ?? 0);
    y[I_CEG] = 0;
    const Sp = (2 * this.spec.geometry.stroke * this.op.rpm) / 60;
    const up = ci.uPrime ?? 0.5 * Sp;
    this.mode = MODE_OPEN;
    this.resetCycleFlags();
    this.startCycle();
    this.onIvc();
    y[I_TKE] = 1.5 * this.mIvc * up * up;
    y[I_TK] = 0;
    this.derivValid = false;
    this.refresh();
  }

  private resetCycleFlags(): void {
    this.dwellSeen = false;
    this.tSparkCmd = NaN;
    this.sparkFired = false;
    this.kernelMassPrev = 0;
    this.handedOff = false;
    this.flameActive = false;
    this.burnDone = false;
    this.misfire = false;
    this.kernelQuenched = false;
    this.knockOnset = false;
    this.knockBurning = false;
    this.knockOnsetDeg = NaN;
    this.knockEndGasFraction = 0;
    this.mapo = 0;
    this.knockQdot = 0;
    this.tKnockOnset = NaN;
    this.knockC0 = 0;
    this.tauAB = NaN;
    this.sparkEnergyPrev = 0;
    this.electrodeLossPrev = 0;
    this.nNO = 0;
    this.ca10 = NaN;
    this.ca50 = NaN;
    this.ca90 = NaN;
    this.peakP = 0;
    this.peakPDeg = -360;
    this.maxDp = 0;
    this.dpGridK = -1;
    this.dpGridN = 0;
    this.dpKinkN = 0;
    this.noPpm = 0;
    this.coFrac = 0;
    this.xbEvo = 0;
    this.closedHappened = false;
    this.volEff = 0;
    this.hrr = 0;
    this.burnRateStep = 0;
  }

  /** Composition / boundary data that depend on the operating point (cycle start). */
  private prepareCycleData(): void {
    const op = this.op;
    this.fuel = fuelFromSelection(op.fuel);
    humidAir(op.ambientTemperature, op.ambientPressure, op.relativeHumidity, this.airX);
    freshCharge({ fuel: this.fuel, phi: op.equivalenceRatio, airX: this.airX, egrFraction: op.egrFraction }, this.Xfresh);
    gasStateFromTPX(op.intakeMixtureTemperature, op.ambientPressure, this.Xfresh, this.fresh);
    gasStateFromTPX(op.ambientTemperature, op.ambientPressure, this.airX, this.ambient);
    {
      let mt = 0;
      let mf = 0;
      for (let k = 0; k < NS; k++) mt += this.Xfresh[k] * MOLAR_MASS[k];
      for (const sName of FUEL_SPECIES) mf += this.Xfresh[SP[sName]] * MOLAR_MASS[SP[sName]];
      const yAir = 1 - mf / mt - op.egrFraction;
      this.freshPerAir = yAir > 1e-6 ? 1 / yAir : 1;
    }
    const m = this.spec.manifolds;
    this.cdaOutlet = OUTLET_DISCHARGE_COEFFICIENT * 0.25 * Math.PI * m.exhaustOutletDiameter * m.exhaustOutletDiameter;
    this.throttleFor = NaN;
    this.updateVenturi();
    // wall temperatures: lumped model (coolant + R_i Q̄_i, updateWalls), else WallSpec shifted with the
    // coolant (UNVERIFIED: wall-to-coolant ΔT held at its value for the CFR's 100 °C boiling jacket).
    const w = this.spec.walls;
    const Tref = w.referenceCoolantTemperature ?? CFR_COOLANT_TEMPERATURE;
    const dT = op.coolantTemperature - Tref;
    let Tp = w.intakePortTemperature;
    if (this.wallsLumped) {
      const Tc = op.coolantTemperature;
      const e = this.wallExcess;
      this.walls = {
        headTemperature: Tc + e[0],
        pistonTemperature: Tc + e[1],
        linerTemperature: Tc + e[2],
        intakeValveTemperature: Tc + e[3],
        exhaustValveTemperature: Tc + e[4],
      };
      // intake port (in the water-cooled head): its excess over the coolant scales with the head's
      // (UNVERIFIED assumption; the port gets its heat through the head casting)
      if (Tp !== undefined) Tp = Tc + (Tp - Tref) * (e[0] / Math.max(1, w.headTemperature - Tref)) - dT;
    } else {
      this.walls = {
        headTemperature: w.headTemperature + dT,
        pistonTemperature: w.pistonTemperature + dT,
        linerTemperature: w.linerTemperature + dT,
        intakeValveTemperature: w.intakeValveTemperature + dT,
        exhaustValveTemperature: w.exhaustValveTemperature + dT,
      };
    }
    this.wallTavg = (this.walls.headTemperature + this.walls.pistonTemperature + this.walls.linerTemperature) / 3;
    // intake-port heat transfer (Dittus–Boelter over the heated port; see evaluate)
    const Lp = m.intakePortLength ?? 0;
    const Dp = m.intakePortDiameter;
    if (this.opts.intakePortHeatTransferMultiplier > 0 && Tp !== undefined && Lp > 0 && Dp > 0) {
      this.portTw = Tp + dT;
      const Tf = 0.5 * (op.intakeMixtureTemperature + this.portTw); // film temperature
      const mu = mixViscosity(this.Xfresh, Tf);
      const k = mixThermalConductivity(this.Xfresh, Tf);
      const cp = mixCpMass(this.Xfresh, Tf);
      this.portReCoef = 4 / (Math.PI * Dp * mu);
      // h = mult · Nu k/D over the wetted area πDL, Nu = max(3.66, 0.023 Re^0.8 Pr^0.4)
      this.portHA = this.opts.intakePortHeatTransferMultiplier * (k / Dp) * Math.PI * Dp * Lp;
      this.portPr04 = Math.pow((mu * cp) / k, 0.4);
    } else {
      this.portHA = 0;
    }
    this.pCrankcase = op.ambientPressure + this.opts.crankcaseGaugePressure;
    this.delayModel = resolveDelayModel(this.opts.ignitionDelayModel);
    this.updateFriction();
    if (this.muTab[0] === 0) for (let i = 0; i < this.muTab.length; i++) this.muTab[i] = mixViscosity(this.Xfresh, this.muT0 + i * this.muDT);
  }

  private updateVenturi(): void {
    const op = this.op;
    if (op.throttle === this.throttleFor) return;
    const D = this.spec.manifolds.throttleDiameter;
    const At = 0.25 * Math.PI * D * D;
    // The CFR has no throttle plate: the venturi throat is the restriction; a partial opening
    // scales it by the butterfly open-area ratio throttleArea(u)/throttleArea(1).
    const cdaV = this.opts.venturiDischargeCoefficient * At * (throttleArea(D, op.throttle) / throttleArea(D, 1));
    const cdaR = this.opts.intakeRestrictionArea;
    // series restrictions (see options.intakeRestrictionArea)
    this.cdaVenturi = cdaR > 0 ? 1 / Math.sqrt(1 / (cdaV * cdaV) + 1 / (cdaR * cdaR)) : cdaV;
    this.throttleFor = op.throttle;
  }

  private updateFriction(): void {
    const fr = this.opts.friction;
    if (!fr) {
      this.friction.setFmep(0, 0, 0, 1);
      return;
    }
    const rpm = Math.max(60, this.rpm || this.op.rpm);
    const pI = this.pInt > 0 ? this.pInt : this.op.ambientPressure;
    const b = pnhFmep({ ...fr, compressionRatio: this.op.compressionRatio }, rpm, pI, this.op.ambientPressure);
    this.friction.setFromPnh(b, (rpm * 2 * Math.PI) / 60);
  }

  /** Cycle start (θ = −360): events and ignition commands of this cycle. */
  private startCycle(): void {
    const op = this.op;
    const omegaDeg = this.y[I_OM] * RAD2DEG;
    const sparkDeg = -op.sparkAdvanceDeg;
    this.ignCmd.sparkDeg = sparkDeg;
    let dwell = sparkDeg - op.dwellTime * omegaDeg;
    if (dwell < -359.9) dwell = -359.9;
    this.ignCmd.dwellStartDeg = dwell;
    this.ivcDeg = this.ivLift.seatCloseDeg;
    this.evoDeg = this.evLift.seatOpenDeg;
    const ci = this.closedInit;
    const cm = this.opts.combustionModel;
    this.nEv = 0;
    const add = (a: number, k: number): void => {
      this.evAngle[this.nEv] = a;
      this.evKind[this.nEv] = k;
      this.nEv++;
    };
    if (ci) {
      add(ci.endDeg, EV_END);
    } else {
      add(this.ivcDeg, EV_IVC);
      add(this.evoDeg, EV_EVO);
      add(this.ivLift.seatOpenDeg, EV_VALVE);
      add(this.evLift.seatCloseDeg, EV_VALVE);
      add(360, EV_WRAP);
    }
    add(-180, EV_BDC_START);
    add(180, EV_BDC_END);
    if (cm === 'instantaneous-at-tdc') add(0, EV_TDC);
    if (cm === 'entrainment') {
      add(dwell, EV_IGN);
      add(sparkDeg, EV_IGN);
    }
    if (cm === 'wiebe') {
      this.wiebeStart = this.opts.wiebe.startDeg ?? sparkDeg;
      add(this.wiebeStart, EV_WIEBE);
    }
    // sort (insertion)
    for (let i = 1; i < this.nEv; i++) {
      const a = this.evAngle[i];
      const k = this.evKind[i];
      let j = i - 1;
      while (j >= 0 && this.evAngle[j] > a) {
        this.evAngle[j + 1] = this.evAngle[j];
        this.evKind[j + 1] = this.evKind[j];
        j--;
      }
      this.evAngle[j + 1] = a;
      this.evKind[j + 1] = k;
    }
    this.evNext = 0;
    const th = this.y[I_TH];
    while (this.evNext < this.nEv && this.evAngle[this.evNext] <= th + 1e-9) {
      if (this.evKind[this.evNext] === EV_BDC_START) this.wBdcStart = this.y[I_W];
      this.evNext++;
    }
    this.wCycleStart = this.y[I_W];
    this.mVentCycleStart = this.ventMass();
    if (cm === 'entrainment' && !this.ign) {
      const io = this.opts.ignition;
      this.ign = new IgnitionSystem(this.spec.ignition, this.spec.sparkPlug, {
        makeSparkDiode: true,
        ...io,
        kernel: { handoffIntegralScaleMultiple: this.opts.kernelHandoffMultiple, ...(io.kernel ?? {}) },
      });
    }
    if (cm !== 'entrainment') this.ign = null;
  }

  // ===========================================================================================
  // Stepping
  // ===========================================================================================

  /** Evaluate the RHS at the current state (end-of-step fields and the cached derivative). */
  refresh(): void {
    this.evaluate(this.t, this.y, this.dy);
    this.derivValid = true;
  }

  private step(tTarget: number, thetaTarget: number): void {
    const y = this.y;
    if (!this.derivValid) this.refresh();
    const P = this.opts.profile;
    const th = y[I_TH];
    this.cacheNOStartState();
    // ---- step limit ----
    // (the angle limits are converted with the same kinematics as the landing time below — in 'free'
    // mode a limit maxStepDeg/ω stopped a hair short of the next grid angle while decelerating and a
    // sliver step followed: 36 % of all steps, validation round 2)
    let hMax = this.timeToAngle(this.opts.maxStepDeg);
    if (this.finePhase()) hMax = Math.min(hMax, this.timeToAngle(this.opts.fineStepDeg));
    if (this.knockBurning) hMax = Math.min(hMax, 0.25 * this.tauAB);
    if (this.stiffness > 0) hMax = Math.min(hMax, 2.0 / this.stiffness);
    if (this.flameActive && this.tauB < Infinity) hMax = Math.min(hMax, Math.max(this.tauB, 1e-7));
    if (this.mode === MODE_TWO) {
      hMax = Math.min(hMax, 2 * this.opts.zoneHeatLossMinTime);
      if (this.mdotB > 0) hMax = Math.min(hMax, Math.max((0.2 * y[I_MU]) / this.mdotB, 1e-9));
    }
    // ---- crank-angle landing: grid, events, target ----
    const grid = this.opts.maxStepDeg;
    let thLand = -360 + (Math.floor((th + 360) / grid + 1e-7) + 1) * grid;
    if (this.evNext < this.nEv && this.evAngle[this.evNext] < thLand + 1e-9) thLand = this.evAngle[this.evNext];
    if (thetaTarget < thLand) thLand = thetaTarget;
    const hLand = this.timeToAngle(thLand - th);
    let h = hMax;
    let landing = false;
    if (hLand <= h * (1 + 1e-9)) {
      h = hLand;
      landing = true;
    }
    const hT = tTarget - this.t;
    if (hT < h) {
      h = hT;
      landing = false;
    }
    if (!(h > 0)) h = 1e-12;
    // ---- tentative RK step (redo for LW crossing / free-mode landing) ----
    this.y0.set(y);
    this.dy0.set(this.dy);
    this.yComp0.set(this.yComp);
    const tStart = this.t;
    for (let attempt = 0; attempt < 6; attempt++) {
      this.rk4(tStart, h);
      if (landing && this.op.speedMode === 'free' && Math.abs(y[I_TH] - thLand) > 1e-10) {
        // secant on the step length to land exactly on the crank event
        const dth = y[I_TH] - th;
        if (dth > 0 && attempt < 5) {
          h *= (thLand - th) / dth;
          this.restore(tStart);
          continue;
        }
      }
      if (landing) y[I_TH] = thLand;
      // Livengood–Wu crossing inside the step → shorten to the crossing
      if (this.lwArmed()) {
        this.tauStepEnd = this.tauNow();
        const f = this.lwCrossingFraction(h, this.tauStepEnd);
        if (f >= 0 && f < 1 - 1e-6 && h * f > 1e-10 && attempt < 5) {
          h *= f;
          landing = false;
          this.restore(tStart);
          this.lwForce = true;
          continue;
        }
      }
      break;
    }
    this.t = tStart + h;
    if (P) this.prof.steps++;
    // ---- operator splits and bookkeeping ----
    this.afterStep(h, tStart);
    // ---- events at the landing angle ----
    if (landing) {
      const mode0 = this.mode;
      const mb0 = y[I_MB];
      this.handleEvents(y[I_TH]);
      // record the post-event state too when an event changed it (IVC/EVO/TDC burn)
      if (this.trace && this.cycle === this.traceCycle && (this.mode !== mode0 || y[I_MB] !== mb0)) {
        if (!this.derivValid) this.refresh();
        this.recordSample(this.p + this.knockOscillation());
      }
      // the pressure interpolation of the next step starts from the post-event state
      if (!this.derivValid) this.refresh();
      this.hTheta = y[I_TH];
      this.hP = this.p;
      this.hDp = this.mode === MODE_OPEN ? 0 : this.dpdt / (y[I_OM] * RAD2DEG);
      this.hMode = this.mode;
    }
  }

  private lwForce = false;
  private tauStepEnd = NaN;
  /** Burned-zone state at the start of the current step for the NO split (Strang): valid, T_b, p, total moles, X. */
  private noStartValid = false;
  private noTb0 = 0;
  private noP0 = 0;
  private noNb0 = 0;
  private readonly noX0 = new Float64Array(NS);
  /** Burned mass merged (with its NO) by mergeToBurned during the current step's splits, kg. */
  private noMergedMass = 0;

  /** Cache the burned-zone equilibrium state at the step start (the last evaluation is at it). */
  private cacheNOStartState(): void {
    const r = this.closure.eq.result;
    this.noStartValid = (this.mode === MODE_TWO || this.mode === MODE_BURNED) && this.y[I_MB] > 0 && r.nTotal > 0 && this.Tb > NO_FREEZE_T;
    if (!this.noStartValid) return;
    this.noTb0 = this.Tb;
    this.noP0 = this.p;
    this.noNb0 = r.nTotal;
    this.noX0.set(r.X);
  }
  /** True once this cycle's dwell has started (the ignition state then belongs to this cycle). */
  dwellSeen = false;
  /** Set by the model when an output-worthy event happened (spark command, knock onset); stepUntil returns. */
  interrupt = false;
  /** Time of this cycle's spark command (primary switch-off), s (NaN before). */
  tSparkCmd = NaN;

  private restore(tStart: number): void {
    this.y.set(this.y0);
    this.dy.set(this.dy0);
    this.yComp.set(this.yComp0);
    this.t = tStart;
  }

  private rk4(t: number, h: number): void {
    const y = this.y;
    const k1 = this.dy;
    const k2 = this.k2;
    const k3 = this.k3;
    const k4 = this.k4;
    const yt = this.yt;
    const hh = 0.5 * h;
    for (let i = 0; i < NY; i++) yt[i] = y[i] + hh * k1[i];
    this.evaluate(t + hh, yt, k2);
    for (let i = 0; i < NY; i++) yt[i] = y[i] + hh * k2[i];
    this.evaluate(t + hh, yt, k3);
    for (let i = 0; i < NY; i++) yt[i] = y[i] + h * k3[i];
    this.evaluate(t + h, yt, k4);
    const h6 = h / 6;
    // compensated (Kahan) summation of the step increments: the conserved states and ledgers exchange
    // many small increments per cycle (mass between the zones and the crevice, the plenums and the
    // ledgers), whose rounding otherwise accumulates to ≈ 1e-12 of the network mass over 12 cycles
    const cmp = this.yComp;
    for (let i = 0; i < NY; i++) {
      const inc = h6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) - cmp[i];
      const t2 = y[i] + inc;
      cmp[i] = t2 - y[i] - inc;
      y[i] = t2;
    }
    // end-of-step evaluation (becomes k1 of the next step unless a split changes the state)
    this.evaluate(t + h, y, this.dy);
    this.derivValid = true;
  }

  /** Time to advance the crank by dθ (deg): exact for fixed speed, quadratic in free mode. */
  private timeToAngle(dth: number): number {
    const om = this.y[I_OM] * RAD2DEG;
    if (this.op.speedMode === 'fixed' || this.alpha === 0) return dth / om;
    const a = 0.5 * this.alpha * RAD2DEG;
    const disc = om * om + 4 * a * dth;
    if (disc <= 0) return dth / om;
    return (2 * dth) / (om + Math.sqrt(disc));
  }

  private finePhase(): boolean {
    if (!this.ign) return false;
    const th = this.y[I_TH];
    if (th < this.ignCmd.sparkDeg - 1e-9 && !this.sparkFired) return false;
    if (this.misfire || this.kernelQuenched || this.burnDone || this.mode === MODE_OPEN) return false;
    if (!this.handedOff) return this.sparkFired || th >= this.ignCmd.sparkDeg - 1e-9;
    return this.xb < 0.02;
  }

  // ===========================================================================================
  // Right-hand side
  // ===========================================================================================

  private evaluate(t: number, y: Float64Array, dy: Float64Array): void {
    const P = this.opts.profile;
    const tp = P ? now() : 0;
    dy.fill(0);
    const th = y[I_TH];
    const om = y[I_OM];
    const ks = this.kin.evaluate(th * DEG, this.ks);
    this.V = ks.volume;
    this.Vdot = ks.dVdTheta * om;
    this.h = ks.clearanceHeight;
    this.updateVenturi();
    // ---- plenums ----
    const ip = this.intake;
    const ep = this.exhaust;
    for (let k = 0; k < NS; k++) {
      ip.N[k] = y[I_IN + k];
      ep.N[k] = y[I_EN + k];
    }
    ip.U = y[I_IU];
    ep.U = y[I_EU];
    ip.updateState();
    ep.updateState();
    this.pInt = ip.state.p;
    this.pExh = ep.state.p;
    // venturi: a = ambient AIR (humid, at the intake-air temperature), b = intake plenum. The CFR
    // carburettor meters air; the fuel (and EGR) join it in/after the throat and the MON mixture
    // heater is downstream, so the throat flow is air at T_amb and the fresh-charge stream is
    // ṁ_air/Y_air,fresh (round 1 flowed the heated, vaporised mixture through the throat: 18 % less
    // mass at a given depression in MON — validation round 2). Backflow: plenum gas.
    let mV = orificeFlow(this.cdaVenturi, this.ambient.p, this.ambient.T, this.ambient.R, this.ambient.gamma, ip.state.p, ip.state.T, ip.state.R, ip.state.gamma, this.ofV);
    if (mV > 0) {
      const f = this.freshPerAir;
      mV *= f;
      this.ofV.dmdotdpa *= f;
      this.ofV.dmdotdpb *= f;
    }
    // outlet: a = exhaust plenum, b = ambient
    const mO = orificeFlow(this.cdaOutlet, ep.state.p, ep.state.T, ep.state.R, ep.state.gamma, this.ambient.p, this.ambient.T, this.ambient.R, this.ambient.gamma, this.ofO);
    this.mdotV = mV;
    this.mdotO = mO;
    dy[I_GIV] = -this.ofV.dmdotdpb;
    if (this.mode === MODE_OPEN) this.evalOpen(y, dy, mV, mO);
    else this.evalClosed(t, y, dy, mV, mO);
    // ---- intake-port heat transfer into the intake plenum (the lumped runner + port volume):
    // Dittus–Boelter Nu = 0.023 Re^0.8 Pr^0.4 for turbulent pipe flow (Re > 1e4, 0.7 ≤ Pr ≤ 160),
    // laminar floor Nu = 3.66 (fully developed, constant wall temperature) — both as given by
    // Ferrantelli, Méloïs & Viljanen 2013, arXiv:1308.2784, eqs. 29 and 34 (fetched) — with the port
    // Reynolds number of the intake-valve flow, times intakePortHeatTransferMultiplier (Choi et al.
    // 2018 GT-Power TPA: 4, range 2–6) ----
    let Qp = 0;
    if (this.portHA > 0 && this.mdotIv !== 0) {
      // (only while gas flows through the port: the stagnant port gas of the lumped plenum is not
      // exchanged with the cylinder, and the closed-phase heat ledger stays the cylinder's)
      const Re = this.portReCoef * Math.abs(this.mdotIv);
      let Nu = 0.023 * Math.pow(Re, 0.8) * this.portPr04;
      if (Nu < 3.66) Nu = 3.66;
      Qp = this.portHA * Nu * (this.portTw - ip.state.T);
      ip.dUdt += Qp;
      dy[I_GIK] = this.portHA * Nu;
    }
    this.Qport = Qp;
    // ---- plenum boundary flows (both phases) ----
    ip.addFlow(mV, this.fresh);
    ep.addFlow(-mO, this.ambient);
    for (let k = 0; k < NS; k++) {
      dy[I_IN + k] = ip.dNdt[k];
      dy[I_EN + k] = ep.dNdt[k];
    }
    dy[I_IU] = ip.dUdt;
    dy[I_EU] = ep.dUdt;
    const mI = massOf(ip.N);
    const mE = massOf(ep.N);
    const ybgI = mI > 0 ? y[I_IBG] / mI : 0;
    const ybgE = mE > 0 ? y[I_EBG] / mE : 0;
    const yegI = mI > 0 ? y[I_IEG] / mI : 0;
    // the fresh charge carries its external EGR (complete-combustion products, freshCharge) as
    // burned gas (round 1 marked it 0: the trapped dilution came from op.egrFraction instead)
    const yEgrFresh = this.op.egrFraction;
    dy[I_IBG] += mV > 0 ? mV * yEgrFresh : mV * ybgI;
    dy[I_IEG] += mV > 0 ? mV * yEgrFresh : mV * yegI;
    dy[I_EBG] -= mO > 0 ? mO * ybgE : 0;
    // ledgers of the boundaries
    const sV = mV > 0 ? this.fresh : ip.state;
    const sO = mO > 0 ? ep.state : this.ambient;
    dy[I_HV] = mV * sV.h;
    dy[I_HO] = mO * sO.h;
    // gross plenum inflows through the boundaries (fresh charge into the intake, ambient backflow
    // into the exhaust)
    if (mV > 0) {
      dy[I_GIH] += mV * sV.h;
      for (let k = 0; k < NS; k++) dy[I_GIN + k] += (mV * sV.Y[k]) / MOLAR_MASS[k];
    }
    if (mO < 0) {
      dy[I_GEH] -= mO * sO.h;
      for (let k = 0; k < NS; k++) dy[I_GEN + k] -= (mO * sO.Y[k]) / MOLAR_MASS[k];
    }
    for (let k = 0; k < NS; k++) {
      dy[I_LV + k] = (mV * sV.Y[k]) / MOLAR_MASS[k];
      dy[I_LO + k] = (mO * sO.Y[k]) / MOLAR_MASS[k];
    }
    // ---- turbulence ----
    const ts = this.turbState;
    ts.K = y[I_TK];
    ts.k = y[I_TKE];
    ts.swirl = y[I_SW];
    turbulenceDerivatives(ts, this.turbIn, this.turbOut, this.turbParams);
    dy[I_TK] = this.turbOut.dK;
    dy[I_TKE] = this.turbOut.dk;
    dy[I_SW] = this.turbOut.dSwirl;
    // ---- mechanics ----
    dy[I_TH] = om * RAD2DEG;
    if (this.op.speedMode === 'free') {
      const ext = -this.op.loadTorque + this.friction.torque(ks.dxdTheta, om);
      let a = this.dyn.angularAcceleration(th * DEG, om, this.p, this.pCrankcase, ext);
      // Stall guard (numerical): a time-stepped cycle cannot represent a stopped crank, so the speed
      // is held at FREE_MODE_MIN_RPM when the torques would decelerate it further.
      if (om <= FREE_MODE_MIN_OMEGA && a < 0) a = 0;
      this.alpha = a;
      dy[I_OM] = a;
    } else {
      this.alpha = 0;
    }
    dy[I_W] = this.p * this.Vdot;
    // net wall heat LOSS of the gas network: cylinder walls minus the intake-port gain
    dy[I_Q] = this.Qwall - Qp;
    if (P) {
      this.prof.rhsEvals++;
      const dt = now() - tp;
      if (this.mode === MODE_OPEN) this.prof.rhsOpen += dt;
      else this.prof.rhsClosed += dt;
    }
    void t;
  }

  private evalOpen(y: Float64Array, dy: Float64Array, mV: number, mO: number): void {
    const P = this.opts.profile;
    const c = this.cyl;
    const ip = this.intake;
    const ep = this.exhaust;
    for (let k = 0; k < NS; k++) c.N[k] = y[I_CN + k];
    c.U = y[I_CU];
    c.volume = this.V;
    c.updateState();
    const sc = c.state;
    const si = ip.state;
    const se = ep.state;
    const m = massOf(c.N);
    this.mCyl = m;
    this.p = sc.p;
    this.T = sc.T;
    this.Tu = sc.T;
    this.Tb = 0;
    const th = y[I_TH];
    const cdm = this.opts.dischargeCoefficientMultiplier;
    // intake valve: a = intake plenum, b = cylinder (+ = into the cylinder)
    const Li = this.ivLift.lift(th);
    let mIv = 0;
    let cdaI = 0;
    if (Li > 0) {
      cdaI = cdm * this.ivFlow.effectiveArea(Li, sc.p > si.p);
      mIv = orificeFlow(cdaI, si.p, si.T, si.R, si.gamma, sc.p, sc.T, sc.R, sc.gamma, this.ofIv);
    } else {
      this.ofIv.dmdotdpa = 0;
      this.ofIv.dmdotdpb = 0;
    }
    // exhaust valve: a = cylinder, b = exhaust plenum (+ = out of the cylinder)
    const Le = this.evLift.lift(th);
    let mEv = 0;
    let cdaE = 0;
    if (Le > 0) {
      cdaE = cdm * this.evFlow.effectiveArea(Le, se.p > sc.p);
      mEv = orificeFlow(cdaE, sc.p, sc.T, sc.R, sc.gamma, se.p, se.T, se.R, se.gamma, this.ofEv);
    } else {
      this.ofEv.dmdotdpa = 0;
      this.ofEv.dmdotdpb = 0;
    }
    this.mdotIv = mIv;
    this.mdotEv = mEv;
    // stiffness bound (Gershgorin) of the linearised pressure dynamics
    const kc = (sc.gamma * sc.R * sc.T) / this.V;
    const ki = (si.gamma * si.R * si.T) / ip.volume;
    const ke = (se.gamma * se.R * se.T) / ep.volume;
    const aIv = Math.abs(this.ofIv.dmdotdpa) + Math.abs(this.ofIv.dmdotdpb);
    const aEv = Math.abs(this.ofEv.dmdotdpa) + Math.abs(this.ofEv.dmdotdpb);
    const aV = Math.abs(this.ofV.dmdotdpa) + Math.abs(this.ofV.dmdotdpb);
    const aO = Math.abs(this.ofO.dmdotdpa) + Math.abs(this.ofO.dmdotdpb);
    this.stiffness = Math.max(kc * (aIv + aEv), ki * (aV + aIv), ke * (aEv + aO));
    // heat transfer
    const tq = P ? now() : 0;
    let Q = 0;
    if (this.opts.heatTransfer) {
      const w = this.woschni;
      w.bore = this.spec.geometry.bore;
      w.pressure = sc.p;
      w.temperature = sc.T;
      w.meanPistonSpeed = this.meanPistonSpeed();
      w.phase = 'gas-exchange';
      w.variant = this.opts.woschniVariant;
      const hc = this.opts.woschniMultiplier * (this.opts.heatTransferCorrelation === 'hohenberg' ? hohenbergCoefficient(this.V, sc.p, sc.T, w.meanPistonSpeed) : woschniCoefficient(w));
      this.hcoef = hc;
      flatChamberAreas(this.spec.geometry.bore, this.h, this.spec.intakeValve, this.spec.exhaustValve, this.areas);
      Q = wallHeatLoss(hc, sc.T, this.areas, this.walls, this.heat);
      this.surfaceHeatRates(dy, 1);
    }
    if (P) this.prof.heat += now() - tq;
    this.Qwall = Q;
    this.Qu = Q;
    this.Qb = 0;
    // cylinder balances
    c.beginRates(-Q, this.Vdot);
    c.addFlow(mIv, si);
    c.addFlow(-mEv, se);
    for (let k = 0; k < NS; k++) dy[I_CN + k] = c.dNdt[k];
    dy[I_CU] = c.dUdt;
    const mI = massOf(ip.N);
    const mE = massOf(ep.N);
    const ybgC = m > 0 ? y[I_CBG] / m : 0;
    const ybgI = mI > 0 ? y[I_IBG] / mI : 0;
    const ybgE = mE > 0 ? y[I_EBG] / mE : 0;
    const fIv = mIv > 0 ? mIv * ybgI : mIv * ybgC; // burned gas into the cylinder via the intake valve
    const fEv = mEv > 0 ? mEv * ybgC : mEv * ybgE; // burned gas out of the cylinder via the exhaust valve
    dy[I_CBG] = fIv - fEv;
    dy[I_IBG] = -fIv;
    dy[I_EBG] = fEv;
    // external-EGR marker: in/out through the intake valve, out through the exhaust valve (exhaust
    // backflow is residual, unmarked)
    const yegC = m > 0 ? y[I_CEG] / m : 0;
    const yegI = mI > 0 ? y[I_IEG] / mI : 0;
    const eIv = mIv > 0 ? mIv * yegI : mIv * yegC;
    dy[I_CEG] = eIv - (mEv > 0 ? mEv * yegC : 0);
    dy[I_IEG] = -eIv;
    // gross plenum inflows through the valves (backflow into the intake, outflow into the exhaust)
    if (mIv < 0) {
      dy[I_GIH] -= mIv * sc.h;
      dy[I_GIBG] -= mIv * ybgC;
      for (let k = 0; k < NS; k++) dy[I_GIN + k] -= (mIv * sc.Y[k]) / MOLAR_MASS[k];
    }
    if (mEv > 0) {
      dy[I_GEH] += mEv * sc.h;
      dy[I_GEBG] += mEv * ybgC;
      for (let k = 0; k < NS; k++) dy[I_GEN + k] += (mEv * sc.Y[k]) / MOLAR_MASS[k];
    }
    // plenum valve-side flows (boundary flows are added by the caller)
    ip.beginRates(0, 0);
    ep.beginRates(0, 0);
    ip.addFlow(-mIv, sc);
    ep.addFlow(mEv, sc);
    dy[I_MIVI] = mIv > 0 ? mIv : 0;
    dy[I_MIVO] = mIv < 0 ? -mIv : 0;
    dy[I_MEVO] = mEv > 0 ? mEv : 0;
    dy[I_MEVI] = mEv < 0 ? -mEv : 0;
    // turbulence inputs: inflow jets (intake inflow, exhaust backflow)
    const ti = this.turbIn;
    let mIn = 0;
    let e2 = 0;
    let mOut = 0;
    let swirlIn = 0;
    if (mIv > 0) {
      const v = intakeJetVelocity(mIv, si.rho, cdaI);
      mIn += mIv;
      e2 += mIv * v * v;
      if (this.opts.swirlMomentumEfficiency > 0) swirlIn += this.opts.swirlMomentumEfficiency * this.ivFlow.inflowJet(mIv, v, this.jet).angularMomentumFlux;
    } else mOut -= mIv;
    if (mEv < 0) {
      const v = intakeJetVelocity(-mEv, se.rho, cdaE);
      mIn -= mEv;
      e2 += -mEv * v * v;
    } else mOut += mEv;
    ti.m = m;
    ti.mDotIn = mIn;
    ti.mDotOut = mOut;
    ti.vIn = mIn > 0 ? Math.sqrt(e2 / mIn) : 0;
    ti.L = integralLengthScale(this.h, this.spec.geometry.bore, this.lengthFraction);
    ti.dlnRhoDt = (mIn - mOut) / m - this.Vdot / this.V;
    ti.swirlTorqueIn = swirlIn;
    ti.rho = sc.rho;
    ti.mu = this.viscosityU(sc.T); // swirl wall friction only: viscosity of the last trapped charge
    ti.bore = this.spec.geometry.bore;
    ti.h = this.h;
    this.L = ti.L;
    this.uPrime = turbulenceIntensity(y[I_TKE], m);
    this.SL = 0;
    this.mdotB = 0;
    this.mdotE = 0;
    this.dpdt = 0;
    void mV;
    void mO;
  }

  private evalClosed(t: number, y: Float64Array, dy: Float64Array, mV: number, mO: number): void {
    const P = this.opts.profile;
    const cl = this.closure;
    // volume of the two zones: the cylinder minus the crevice zone
    const V = this.creviceOn ? this.V - this.Vcr : this.V;
    const mu = y[I_MU];
    const mb = y[I_MB];
    const m = mu + mb;
    const U = y[I_UT];
    this.mCyl = m;
    this.mdotIv = 0;
    this.mdotEv = 0;
    this.stiffness = 0;
    // plenum valve-side rates are zero
    this.intake.beginRates(0, 0);
    this.exhaust.beginRates(0, 0);
    // ---- closure ----
    const tc = P ? now() : 0;
    const e0 = cl.eqSolves;
    cl.noCoupled = this.noCoupled;
    cl.noKinetic = this.nNO;
    const dtp = t - this.tEval;
    const mode = this.mode;
    if (mode === MODE_SINGLE) {
      const T = cl.solveSingle(U, mu, V, this.TuGuess);
      this.TuGuess = T;
      this.Tu = T;
      this.Tb = 0;
      this.T = T;
    } else if (mode === MODE_TWO) {
      const pG = this.pGuess + this.dpEval * dtp;
      const TbG = this.TbGuess + this.dTbEval * dtp;
      if (!cl.solveTwoZone(U, y[I_SU], mu, mb, V, pG > 0 ? pG : this.pGuess, TbG > 300 ? TbG : this.TbGuess, this.TuGuess)) {
        // robust retry from the last converged pressure and the HP flame temperature
        const TbR = cl.flameTemperatureHP(this.eqAux, this.TuGuess, this.pGuess);
        if (!cl.solveTwoZone(U, y[I_SU], mu, mb, V, this.pGuess, TbR, this.TuGuess)) this.closureFailures++;
      }
      this.Tu = cl.Tu;
      this.Tb = cl.Tb;
      this.TuGuess = cl.Tu;
      this.T = (mu * cl.Tu + mb * cl.Tb) / m;
    } else {
      const pG = this.pGuess + this.dpEval * dtp;
      const TbG = this.TbGuess + this.dTbEval * dtp;
      if (!cl.solveBurnedOnly(U, mb, V, pG > 0 ? pG : this.pGuess, TbG > 300 ? TbG : this.TbGuess)) {
        if (!cl.solveBurnedOnly(U, mb, V, this.pGuess, this.TbGuess)) this.closureFailures++;
      }
      this.Tu = 0;
      this.Tb = cl.Tb;
      this.T = cl.Tb;
    }
    const p = cl.p;
    this.p = p;
    if (P) {
      this.prof.closure += now() - tc;
      this.prof.eqSolves += cl.eqSolves - e0;
    }
    const vu = mode === MODE_BURNED ? 0 : cl.vu;
    this.rhoU = vu > 0 ? 1 / vu : 0;
    if (mu > 0 && mode !== MODE_BURNED) this.nuU = this.viscosityU(this.Tu) * vu;
    // ---- flame (entrainment / Wiebe) ----
    const tf = P ? now() : 0;
    let mdB = 0;
    let mdE = 0;
    const Vb = mode === MODE_SINGLE ? 0 : mb * cl.vb;
    // turbulence length scale and intensity (unburned zone)
    const bore = this.spec.geometry.bore;
    let L = integralLengthScale(this.h, bore, this.lengthFraction);
    if (mode === MODE_TWO && this.rhoRef > 0 && this.rhoU > 0) {
      const La = angularMomentumLengthScale(this.Lref, this.rhoRef, this.rhoU);
      if (La < L) L = La;
    }
    this.L = L;
    const mTurb = mode === MODE_TWO ? mu : m;
    const up = turbulenceIntensity(y[I_TKE], mTurb);
    this.uPrime = up;
    this.SL = 0;
    if (mode !== MODE_BURNED && mu > 0) this.SL = laminarFlameSpeed(this.fuelTrapped, this.phiTrapped, this.Tu, p, this.xDil);
    this.lambda = this.opts.taylorScaleMultiplier * taylorMicroscale(L, up, this.nuU);
    this.Af = 0;
    this.tauB = Infinity;
    this.burnSpeed = 0;
    // Zone volumes are mapped onto the disc chamber of the flame geometry in proportion to its share
    // of the cylinder volume, discScale = A_p h / V = 1 − V_crevice/V: the lumped crevice volume
    // (part of V, SliderCrank) holds unburned charge that is distributed over the chamber, so the
    // front reaches the far corner exactly when the whole charge is entrained (round 1 mapped V_e
    // onto the disc alone: the last V_crevice·ρ_u ≈ 0.8 % of the charge was never entrained and
    // non-knocking cycles never burned out). Crevice gas thus burns with the charge (no crevice
    // storage / blow-by model).
    const discScale = (this.kin.boreArea * this.h) / V;
    this.discScale = discScale;
    void 0;
    if (mode === MODE_TWO) {
      // equivalent radius of the burned gas in the chamber (flame-centred sphere ∩ disc)
      this.rb = this.flameGeom.radiusForVolume(Vb * discScale, this.h, this.clampRadiusGuess(this.rbGuess, this.h));
      this.rbGuess = this.rb;
    }
    if (mode === MODE_TWO && this.flameActive) {
      const me = y[I_ME];
      const Ve = Vb + (me > mb ? (me - mb) * vu : 0);
      this.rf = this.flameGeom.radiusForVolume(Ve * discScale, this.h, this.clampRadiusGuess(this.rfGuess, this.h));
      this.rfGuess = this.rf;
      this.flameGeom.evaluate(this.rf, this.h, this.fg);
      this.Af = this.fg.frontArea > 0 ? this.fg.frontArea : 0;
      let uT = this.opts.burnRateMultiplier * up;
      if (this.keck && this.rhoInlet > 0) {
        // Keck 1982 eq. 4.10 / Fig. 15 empirical closures (entrainment.ts)
        uT = this.opts.burnRateMultiplier * keckCharacteristicSpeed(this.keckInletSpeed, this.rhoU, this.rhoInlet);
        this.lambda = this.opts.taylorScaleMultiplier * keckCharacteristicLength(this.spec.intakeValve.maxLift, this.rhoU, this.rhoInlet);
      }
      const ei = this.entIn;
      ei.rhoU = this.rhoU;
      ei.frontArea = this.Af;
      ei.uPrime = uT;
      ei.SL = this.SL;
      ei.me = me;
      ei.mb = mb;
      ei.lambda = this.lambda;
      ei.mTotal = m;
      entrainmentRates(ei, this.entOut);
      mdE = this.entOut.dme;
      mdB = this.entOut.dmb;
      this.tauB = this.entOut.tauB;
      this.burnSpeed = this.entOut.burningSpeed;
      if (mdB > mu / 1e-6) mdB = mu / 1e-6;
    } else if (mode === MODE_TWO && this.opts.combustionModel === 'wiebe') {
      const w = this.opts.wiebe;
      const u = (y[I_TH] - this.wiebeStart) / w.durationDeg;
      if (u > 0) {
        const um = Math.pow(u, w.m);
        const dxdth = ((w.a * (w.m + 1) * um) / w.durationDeg) * Math.exp(-w.a * um * u);
        mdB = m * dxdth * y[I_OM] * RAD2DEG;
      }
      if (mu <= 0) mdB = 0;
    }
    // end-gas autoignition burn-up: first order in the end gas AHEAD of the flame front, m − m_e
    // (τ_ab); the entrained-but-unburned brush behind the front keeps burning up on τ_b. The
    // autoignited gas counts as entrained (dm_e/dt += ṁ_K) so the brush mass m_e − m_b is unchanged.
    let mdK = 0;
    if (this.knockBurning && mode === MODE_TWO && mu > 0) {
      // ALL unburned gas autoignites: the end gas ahead of the front and the entrained-but-unburned
      // pockets of the flame brush share the unburned zone's (T_u, p) history (validation round 2:
      // round 1 let only the gas ahead of the front autoignite — 1–3 % of the charge at the CFR's
      // standard knock while the brush held ≈ 10 %). The end gas ahead of the front becomes
      // entrained (dm_e += (m − m_e)/τ_ab), the brush burns (its m_e − m_b decays on τ_ab too).
      const eg = m - y[I_ME];
      const egc = eg > 0 ? (eg < mu ? eg : mu) : 0;
      mdK = endGasBurnRate(this.opts.knockBrushAutoignition ? mu : egc, this.tauAB);
      mdB += mdK;
      mdE += endGasBurnRate(egc, this.tauAB);
      dy[I_MK] = mdK;
    }
    this.mdotB = mdB;
    this.mdotE = mdE;
    if (P) this.prof.flame += now() - tf;
    // ---- heat transfer ----
    const th0 = P ? now() : 0;
    let Qu = 0;
    let Qb = 0;
    if (this.opts.heatTransfer) {
      const w = this.woschni;
      w.bore = bore;
      w.pressure = p;
      w.temperature = this.T;
      w.meanPistonSpeed = this.meanPistonSpeed();
      w.variant = this.opts.woschniVariant;
      if (mode === MODE_SINGLE) {
        w.phase = 'compression';
        this.pMot = p;
      } else {
        w.phase = 'combustion';
        this.pMot = this.motoredPressure(mode === MODE_TWO ? y[I_SU] / mu : this.suLast, m, V);
        w.motoredPressure = this.pMot;
        w.c2 = this.opts.woschniCombustionTermMultiplier * WOSCHNI_CONSTANTS.c2Combustion;
        w.displacedVolume = this.kin.displacedVolume;
        w.refPressure = this.pIvc;
        w.refTemperature = this.TIvc;
        w.refVolume = this.VIvc;
      }
      const hc = this.opts.woschniMultiplier * (this.opts.heatTransferCorrelation === 'hohenberg' ? hohenbergCoefficient(this.V, p, this.T, w.meanPistonSpeed) : woschniCoefficient(w));
      this.hcoef = hc;
      flatChamberAreas(bore, this.h, this.spec.intakeValve, this.spec.exhaustValve, this.areas);
      if (mode === MODE_TWO) {
        this.flameGeom.evaluate(this.rb, this.h, this.fgB);
        wallHeatLossTwoZone(hc, this.Tu, this.Tb, this.areas, this.fgB, this.walls, this.heat);
        Qu = this.heat.unburned;
        Qb = this.heat.burned;
        const q0 = Qu + Qb;
        // numerical guard for vanishing zones (see options.zoneHeatLossMinTime)
        const tg = this.opts.zoneHeatLossMinTime;
        const qmu = (mu * cl.cvu * Math.abs(this.Tu - this.wallTavg)) / tg;
        if (Math.abs(Qu) > qmu) Qu = Math.sign(Qu) * qmu;
        const qmb = (mb * cl.cvb * Math.abs(this.Tb - this.wallTavg)) / tg;
        if (Math.abs(Qb) > qmb) Qb = Math.sign(Qb) * qmb;
        this.surfaceHeatRates(dy, q0 !== 0 ? (Qu + Qb) / q0 : 1);
      } else {
        const q = wallHeatLoss(hc, this.T, this.areas, this.walls, this.heat);
        if (mode === MODE_SINGLE) Qu = q;
        else Qb = q;
        this.surfaceHeatRates(dy, 1);
      }
    }
    if (P) this.prof.heat += now() - th0;
    this.Qu = Qu;
    this.Qb = Qb;
    // ---- balances ----
    let dU = -p * this.Vdot - Qu - Qb;
    let dmu = -mdB;
    let dmb = mdB;
    let dSnet = mode === MODE_TWO ? -Qu / this.Tu : 0;
    let dme = mdE;
    let Qcr = 0;
    let crIn = 0; // unburned mass flow bulk → crevice (for the turbulence of the unburned zone)
    let crRet = 0; // crevice → unburned zone
    this.mdotCr = 0;
    if (this.creviceOn) {
      // ---- crevice zone (Namazian & Heywood 1982): gas at T_cr and the cylinder pressure, m_cr R_cr =
      // p V_cr/T_cr. Its exchange with the zones follows ṗ: ṁ = (V_cr/(R_in T_cr)) ṗ (inflow) or
      // (V_cr/(R_cr T_cr)) |ṗ| (outflow), with ṗ itself depending on ṁ (linear closure rates):
      // ṁ = c ṗ₀/(1 − c ∂ṗ/∂ṁ). Inflow comes from the zone at the crevice mouth (the piston top land:
      // burned fraction f_b = burned-wetted share of the liner); outflow carries the crevice mix, its
      // unburned part back to the unburned zone (share 1 − f_b) or into the burned gas (share f_b,
      // where it burns at equilibrium), its burned part (complete-combustion products) to the burned
      // zone. Energy: the bulk loses/gains the stream enthalpy; the crevice stores u(T_cr); the rest is
      // heat to the piston/liner (Q_cr = Σ ṁ_in (h_in − u_cr) − Σ ṁ_out R T_cr). The returning unburned
      // gas enters the unburned zone with dS_u = ṁ (s_u + (h_cr − h_u)/T_u) — the energy-consistent
      // entropy of adding colder gas at constant p. ----
      const mcu = y[I_CRU];
      const mcb = y[I_CRB];
      const mcr = mcu + mcb;
      let fb = mode === MODE_BURNED ? 1 : 0;
      if (mode === MODE_TWO && this.areas.liner > 0) {
        fb = this.fgB.wettedLiner / this.areas.liner;
        fb = fb > 0 ? (fb < 1 ? fb : 1) : 0;
      }
      const dp0 = this.closedDp(mode, dU, dSnet, dmu, dmb, this.Vdot, mu, V, p);
      const hu = cl.hu;
      const hb = cl.hb;
      const Ru = cl.Ru;
      const Rcc = this.Rcc;
      const Tcr = this.Tcr;
      // constraint keeper (numerical): the RK integration of ṁ = c ṗ drifts from m_cr R_cr = p V_cr/T_cr
      // where the splits move p discontinuously and by the small bias of the linearised ṗ (≈ 1 % over
      // the expansion without it); the residual is relaxed on CREVICE_RELAX_TIME (J/K per s)
      const nRerr = (p * this.Vcr) / Tcr - (mcu * Ru + mcb * Rcc);
      const corr = nRerr / CREVICE_RELAX_TIME;
      if (dp0 * this.Vcr / Tcr + corr >= 0 || !(mcr > 0)) {
        // inflow (per unit mass): from unburned (1 − f_b) and burned (f_b)
        const eU = -((1 - fb) * hu + fb * hb);
        const dp1 = this.closedDp(mode, eU, 0, -(1 - fb), -fb, 0, mu, V, p);
        const Rin = (1 - fb) * Ru + fb * Rcc;
        const c = this.Vcr / (Rin * Tcr);
        let md = (c * dp0 + corr / Rin) / (1 - c * dp1);
        if (!(md > 0) || !Number.isFinite(md)) md = 0;
        dU += md * eU;
        dmu -= (1 - fb) * md;
        dmb -= fb * md;
        dme -= fb * md;
        dy[I_CRU] = (1 - fb) * md;
        dy[I_CRB] = fb * md;
        Qcr = md * ((1 - fb) * (hu - this.uuCr) + fb * (hb - this.uccCr));
        crIn = (1 - fb) * md;
        this.mdotCr = md;
      } else {
        // outflow (per unit mass) with the crevice composition
        const yu = mcu / mcr;
        const toU = mode === MODE_BURNED ? 0 : yu * (1 - fb);
        const toB = 1 - toU;
        const eU = yu * this.huCr + (1 - yu) * this.hccCr;
        const eS = mode === MODE_BURNED ? 0 : (toU * (this.huCr - hu)) / this.Tu;
        const dp1 = this.closedDp(mode, eU, eS, toU, toB, 0, mu, V, p);
        const Rout = yu * Ru + (1 - yu) * Rcc;
        const c = this.Vcr / (Rout * Tcr);
        let md = (-c * dp0 - corr / Rout) / (1 + c * dp1);
        if (!(md > 0) || !Number.isFinite(md)) md = 0;
        dU += md * eU;
        dSnet += md * eS;
        dmu += toU * md;
        dmb += toB * md;
        dme += toB * md;
        dy[I_CRU] = -yu * md;
        dy[I_CRB] = -(1 - yu) * md;
        Qcr = -md * (yu * Ru + (1 - yu) * Rcc) * Tcr;
        crRet = toU * md;
        this.mdotCr = -md;
      }
      dy[I_QS + 1] += 0.5 * Qcr;
      dy[I_QS + 2] += 0.5 * Qcr;
    }
    this.Qcr = Qcr;
    this.Qwall = Qu + Qb + Qcr;
    dy[I_UT] = dU;
    dy[I_MU] = dmu;
    dy[I_MB] = dmb;
    dy[I_ME] = dme;
    if (mode === MODE_TWO) {
      dy[I_SU] = cl.su * dmu + dSnet;
      cl.rates(dU, dSnet, dmu, dmb, this.Vdot, mu);
      this.dpdt = cl.dp;
    } else if (mode === MODE_BURNED) {
      cl.rates(dU, 0, 0, dmb, this.Vdot, 0);
      this.dpdt = cl.dp;
    } else {
      this.dpdt = this.closedDp(mode, dU, 0, dmu, 0, this.Vdot, mu, V, p);
    }
    // guesses for the next closure solve (linear predictor in time)
    this.tEval = t;
    this.pGuess = p;
    if (mode !== MODE_SINGLE) {
      this.TbGuess = cl.Tb;
      this.dpEval = this.dpdt;
      this.dTbEval = mode === MODE_SINGLE ? 0 : cl.dTb;
    } else {
      this.dpEval = this.dpdt;
      this.dTbEval = 0;
    }
    // ---- turbulence inputs (unburned zone during combustion, whole charge otherwise) ----
    const ti = this.turbIn;
    ti.m = mTurb;
    ti.mDotIn = crRet; // (crevice gas returning at rest)
    ti.mDotOut = (mode === MODE_TWO ? mdB : 0) + crIn;
    ti.vIn = 0;
    ti.L = L;
    if (mode === MODE_TWO) ti.dlnRhoDt = this.dpdt / (cl.gu * p) + Qu / (this.Tu * mu * cl.cpu);
    else ti.dlnRhoDt = -this.Vdot / V;
    ti.swirlTorqueIn = 0;
    ti.rho = m / V;
    ti.mu = this.viscosityU(this.Tu > 0 ? this.Tu : this.T); // swirl wall friction only
    ti.bore = bore;
    ti.h = this.h;
    void mV;
    void mO;
  }

  /**
   * dp/dt of the closed charge for the given state rates (linear in them): two zones / burned only via
   * the closure Jacobian (ZoneClosure.rates, overwrites its dp/dTb), single zone from U = m u(T),
   * p = m R T/V. dSnet = dS_u/dt − s_u dm_u/dt.
   */
  private closedDp(mode: number, dU: number, dSnet: number, dmu: number, dmb: number, dV: number, mu: number, V: number, p: number): number {
    const cl = this.closure;
    if (mode === MODE_TWO) {
      cl.rates(dU, dSnet, dmu, dmb, dV, mu);
      return cl.dp;
    }
    if (mode === MODE_BURNED) {
      cl.rates(dU, 0, 0, dmb, dV, 0);
      return cl.dp;
    }
    const dT = (dU - cl.uu * dmu) / (mu * cl.cvu);
    return p * (dmu / mu + dT / this.Tu - dV / V);
  }

  /** Per-surface heat-flow ledgers (lumped wall model), scaled by the zone guard factor f. */
  private surfaceHeatRates(dy: Float64Array, f: number): void {
    const q = this.heat;
    dy[I_QS] = f * q.head;
    dy[I_QS + 1] = f * q.piston;
    dy[I_QS + 2] = f * q.liner;
    dy[I_QS + 3] = f * q.intakeValves;
    dy[I_QS + 4] = f * q.exhaustValves;
  }

  /**
   * Lumped wall model (spec.walls.thermalResistance; validation round 2 — the walls had been the fired
   * standard-knock values of Pal 2018 Table 3 at every operating point, including the MOTORED ASTM
   * compression check): at the end of each cycle each surface's excess over the coolant becomes
   * R_i·Q̄_i, Q̄_i = the cycle-mean gas-to-surface heat flow (steady conduction through the wall to the
   * coolant; R_i fitted to Pal 2018 Table 3 at the Choi/Pal standard-knock state, cfr.ts). During warm-up
   * the walls jump to that value (fixed-point iteration, contraction ≈ R·hA ≈ 0.1); afterwards they
   * relax with options.wallTimeConstant (0 = quasi-steady: each cycle sees the walls of the previous
   * cycle's heat load).
   */
  private updateWalls(dtCycle: number, instant: boolean): void {
    const R = this.spec.walls.thermalResistance;
    if (!R || !this.wallsLumped || !(dtCycle > 0)) return;
    const y = this.y;
    const tau = this.opts.wallTimeConstant;
    const a = instant || !(tau > 0) ? 1 : -Math.expm1(-dtCycle / tau);
    const rs = [R.head, R.piston, R.liner, R.intakeValve, R.exhaustValve];
    for (let i = 0; i < 5; i++) {
      const qBar = (y[I_QS + i] - this.qsCycleStart[i]) / dtCycle;
      let target = rs[i] * qBar;
      if (target < -50) target = -50;
      if (target > 600) target = 600;
      if (Number.isFinite(target)) this.wallExcess[i] += a * (target - this.wallExcess[i]);
    }
  }

  /** Mean piston speed at the current ω, m/s. */
  private meanPistonSpeed(): number {
    return (2 * this.spec.geometry.stroke * Math.abs(this.y[I_OM])) / (2 * Math.PI);
  }

  /**
   * Motored pressure for Woschni's combustion term: the pressure the whole charge (mass m) would
   * have at volume V on the unburned-zone isentrope s_u (which carries the heat-loss history, so
   * p_mot = p up to the first burned gas and needs no polytropic-exponent constant).
   * Newton on T: s(T, mRT/V) = s_u, ds/dT|_v = c_v/T.
   */
  private motoredPressure(su: number, m: number, V: number): number {
    const cl = this.closure;
    const Xu = cl.Xu;
    let T = this.TmotGuess;
    for (let it = 0; it < 30; it++) {
      const pT = (m * cl.Ru * T) / V;
      const f = mixSMass(Xu, T, pT) - su;
      const d = mixCvMass(Xu, T) / T;
      const dT = -f / d;
      T += dT;
      if (Math.abs(dT) < 1e-9 * T) break;
    }
    this.TmotGuess = T;
    return (m * cl.Ru * T) / V;
  }

  /** Radius guess kept strictly inside the bracket of FlameGeometry.radiusForVolume (deterministic warm start). */
  private clampRadiusGuess(r: number, h: number): number {
    const fg = this.flameGeom;
    const b = fg.headDistance;
    const lo = h > b ? Math.min(b, h - b, fg.radius - fg.offset) : 0;
    const hi = fg.maxRadius(h);
    const span = hi - lo;
    const a = lo + 1e-6 * span;
    const z = hi - 1e-6 * span;
    return r > a ? (r < z ? r : z) : a;
  }

  /** Dynamic viscosity of the frozen unburned mixture from the per-IVC table, Pa s. */
  private viscosityU(T: number): number {
    const x = (T - this.muT0) / this.muDT;
    const n = this.muTab.length;
    let i = Math.floor(x);
    if (i < 0) i = 0;
    if (i > n - 2) i = n - 2;
    const f = x - i;
    return this.muTab[i] + f * (this.muTab[i + 1] - this.muTab[i]);
  }

  // ===========================================================================================
  // Operator splits, bookkeeping
  // ===========================================================================================

  private afterStep(h: number, tStart: number): void {
    const y = this.y;
    const mbStart = this.y0[I_MB];
    const muStart = this.y0[I_MU];
    const P = this.opts.profile;
    let changed = false;
    this.noMergedMass = 0;
    const mbRK = y[I_MB]; // burned mass at the end of the RK step, before the splits
    // End-of-RK-step state BEFORE the splits: the end-gas burn of this step happened in the
    // unburned zone at this temperature (a burn-out merge below sets T_u = 0 in MODE_BURNED —
    // evaluating the knock source there gave heatOfReaction(0) = NaN; validation round 1).
    const TuStep = this.Tu > 0 ? this.Tu : this.TIvc;
    const pStepEnd = this.p;
    const dpStepEnd = this.dpdt;
    // autoigniting gas at the step start / end (its decay shapes the knock source): all unburned gas
    // (knockBrushAutoignition) or the end gas ahead of the front
    const brush = this.opts.knockBrushAutoignition;
    const eg0 = brush ? this.y0[I_MU] : this.y0[I_MU] + this.y0[I_MB] - this.y0[I_ME];
    const eg1 = brush ? y[I_MU] : y[I_MU] + y[I_MB] - y[I_ME];
    // ---- 1. ignition ----
    if (this.ign) {
      const ti = P ? now() : 0;
      changed = this.ignitionSplit(h) || changed;
      if (P) this.prof.ignition += now() - ti;
    }
    // ---- 2. end-gas burn-up: inside the RK right-hand side (dm_b/dt += (m − m_e)/τ_ab with steps ≤
    // τ_ab/4); the mass it burned over this step drives the acoustic modes below ----
    const dmKnock = this.knockBurning ? y[I_MK] - this.y0[I_MK] : 0;
    if (this.knockBurning && this.mode === MODE_TWO) {
      // end of the autoignition burn-up: the unburned gas is consumed; the e^{−t/τ_ab} tail below 1e-9
      // of the charge is moved at constant U_tot (the burn-out merge then follows).
      const m = y[I_MU] + y[I_MB];
      const left = this.opts.knockBrushAutoignition ? y[I_MU] : Math.min(y[I_MU], m - y[I_ME]);
      if (left <= 1e-9 * m) {
        if (left > 0) {
          this.transferToBurned(left);
          changed = true;
        }
        if (this.mode === MODE_TWO && y[I_ME] < y[I_MU] + y[I_MB]) y[I_ME] = y[I_MU] + y[I_MB];
        this.knockBurning = false;
      }
    }
    // ---- 3. burn-out ----
    if (this.mode === MODE_TWO) {
      const mu = y[I_MU];
      const m = mu + y[I_MB];
      if (mu <= this.opts.burnoutFraction * m) {
        this.mergeToBurned();
        changed = true;
      }
    }
    if (changed) this.refresh();
    // ---- 4. NO (burned zone): Strang splitting — half a step of the rate-controlled kinetics on the
    // step-START burned-zone state, the gas burned during the step enters at mid-step with the unburned
    // NO, half a step on the step-END state (second order in h; round 1 advanced the whole step on the
    // end state with the new gas added at the start: first order, 4.6 % in-cycle error at 0.25° and a
    // summary NO that moved 0.5 % with the snapshot cadence — validation round 2) ----
    if ((this.mode === MODE_TWO || this.mode === MODE_BURNED) && y[I_MB] > 0) {
      const tn = P ? now() : 0;
      // burned gas that entered the crevice this step takes its NO along (frozen there, returned as
      // complete-combustion products — the crevice's NO is not tracked)
      const dcrb = y[I_CRB] - this.y0[I_CRB];
      if (dcrb > 0 && y[I_MB] > 0) this.nNO *= y[I_MB] / (y[I_MB] + dcrb);
      // (mass merged by mergeToBurned in this step already brought its NO)
      const dmb = y[I_MB] - mbStart - this.noMergedMass;
      const hh = 0.5 * h;
      if (this.noStartValid) {
        const x0 = this.nNO / this.noNb0;
        this.nNO = this.zeld.advanceRateControlled(this.noTb0, this.noP0, this.noX0, x0, hh) * this.noNb0;
      }
      if (dmb > 0) this.nNO += dmb * this.closure.noPerKg;
      const r = this.closure.eq.result;
      const nb = r.nTotal;
      if (nb > 0 && this.Tb > NO_FREEZE_T) {
        const x = this.nNO / nb;
        const x1 = this.zeld.advanceRateControlled(this.Tb, this.p, r.X, x, hh);
        this.nNO = x1 * nb;
      }
      if (P) this.prof.no += now() - tn;
    }
    // ---- 5. knock oscillator ----
    if (this.knockOnset && this.mode !== MODE_OPEN) {
      const tk = P ? now() : 0;
      const cl = this.closure;
      // chemical heat of the end gas burned in this step, at the unburned temperature of the RK step
      // end (captured before the splits: T_u = 0 after a burn-out merge)
      const qc = dmKnock > 0 ? cl.heatOfReaction(TuStep) : 0;
      this.knockQdot = dmKnock > 0 ? (dmKnock * qc) / h : 0;
      // mode frequencies (sound speed), γ and V vary linearly in time over the step, from the
      // previous step's end to this one's (a piecewise-constant frequency made the ringing
      // amplitude — and MAPO — depend on where the steps fell, ±2–4 %; validation round 1)
      const c1 = this.knockSoundSpeed();
      const g1 = this.knockGamma;
      const V1 = this.V;
      const c0 = this.knockC0 > 0 ? this.knockC0 : c1;
      const g0 = this.knockC0 > 0 ? this.knockG0 : g1;
      const V0 = this.knockC0 > 0 ? this.knockV0 : V1;
      const osc = this.knockOsc;
      const dsMax = this.opts.knockBurnStep;
      const burn = dmKnock > 0 && qc > 0;
      if (burn || osc.sensorEnvelope(c1) > this.mapo) {
        // Sub-step the acoustic modes at ≤ knockBurnStep: during the burn-up the end-gas mass burned
        // over this RK step is released with the exponential profile e^{−t/τ_ab} of its first-order
        // burn-up, which KnockOscillator.step propagates EXACTLY (no zero-order hold); and while the
        // free ringing can still exceed the current MAPO, so the tracked peak (cubic Hermite between
        // sub-steps, KnockOscillator.peakSensorPressure) does not depend on the integrator's steps.
        const n = Math.ceil(h / dsMax - 1e-9);
        const ds = h / n;
        // Q̇(t) = Q̇₀ e^{−t/τ} with ∫₀^h Q̇ dt = q_c Δm_K (propagated exactly by the oscillator); τ is
        // the decay time of the end gas over this step — τ_ab, shortened by the flame still entraining
        // it — taken from its values at both ends (exact for an exponential decay within the step)
        let tau = this.tauAB;
        if (eg0 > 0 && eg1 > 0 && eg1 < eg0) tau = h / Math.log(eg0 / eg1);
        const q0 = burn ? (qc * dmKnock) / (tau * -Math.expm1(-h / tau)) : 0;
        // sequential source: the autoignited fraction F of the onset end gas selects the shell that is
        // releasing its heat (F follows the same exponential profile within the step)
        const seq = burn && osc.shellCount > 0 && this.knockEgOnset > 0;
        const F0 = seq ? (this.y0[I_MK] - this.knockMkOnset) / this.knockEgOnset : 0;
        const dF = seq ? dmKnock / this.knockEgOnset : 0;
        const den = -Math.expm1(-h / tau);
        for (let i = 1; i <= n; i++) {
          const q = burn ? q0 * Math.exp((-(i - 1) * ds) / tau) : 0;
          const f = (i - 0.5) / n;
          if (seq) osc.setSequentialPosition(F0 + (dF * -Math.expm1((-(i - 0.5) * ds) / tau)) / den);
          osc.step(ds, c0 + f * (c1 - c0), g0 + f * (g1 - g0), V0 + f * (V1 - V0), q, burn ? tau : Infinity);
        }
      } else {
        osc.step(h, 0.5 * (c0 + c1), 0.5 * (g0 + g1), 0.5 * (V0 + V1), 0);
      }
      this.knockC0 = c1;
      this.knockG0 = g1;
      this.knockV0 = V1;
      // MAPO: max |p_osc| at the pickup including the between-sample maxima (sampling-independent)
      if (osc.peakSensorPressure > this.mapo) this.mapo = osc.peakSensorPressure;
      if (P) this.prof.knock += now() - tk;
    }
    // ---- 6. Livengood–Wu (official integral) and onset ----
    if (this.lwArmed()) {
      // τ at the end of the RK step, before the splits (the same value the crossing test used)
      const tau = this.tauStepEnd;
      this.lw.advance(h, tau);
      this.lwTauPrev = tau;
      // history integral J = ∫ (1/τ)(∂lnτ/∂T) T dt for the ignition-time spread of a stratified end
      // gas: 1/τ exponential and g = (∂lnτ/∂T)T linear over the step (exact for that interpolant,
      // consistent with the log-mean quadrature of the LW integral itself)
      const tauT = this.twoStage ? this.tauTotal() : tau;
      const r1 = tauT > 0 && Number.isFinite(tauT) ? 1 / tauT : 0;
      const g1 = this.lwSensitivity(TuStep, pStepEnd);
      if (!Number.isNaN(this.lwRPrev)) this.lwJ += h * expLinearIntegral(this.lwRPrev, r1, this.lwGPrev, g1);
      this.lwRPrev = r1;
      this.lwGPrev = g1;
      if (this.lw.integral >= 1 || this.lwForce) {
        // (trace: the pre-onset state at this instant too — the rates jump at the onset, e.g. the
        // crevice inflow that follows ṗ, and the trace must carry both sides)
        if (this.trace && this.cycle === this.traceCycle && !(this.twoStage && this.lwStage === 1)) this.recordSample(this.p + this.knockOscillation());
        if (this.twoStage && this.lwStage === 1) {
          // first-stage (cool-flame) crossing: start the second-stage integral here
          this.lwStage = 2;
          const t2 = this.tauNow();
          this.lw.reset(t2);
          this.lwTauPrev = t2;
          if (!(t2 > 0)) this.onAutoignition(); // single-stage chemistry: hot ignition now
        } else this.onAutoignition();
      }
    }
    this.lwForce = false;
    // ---- 7. rates of the step, cycle accumulators, trace ----
    if (!this.derivValid) this.refresh();
    const dmbStep = y[I_MB] - mbStart;
    this.burnRateStep = this.mode !== MODE_OPEN && dmbStep > 0 ? dmbStep / h : 0;
    const Tq = this.Tu > 0 ? this.Tu : TuStep;
    // heat-release rate for output: the INSTANTANEOUS burn rate of the end-of-step right-hand side
    // plus the step mean of the split transfers (kernel mass, burn-out merge) where they act — the
    // round-1 step mean depended on the step and sampling pattern by 6–11 % (validation round 2)
    const splitRate = this.mode !== MODE_OPEN && y[I_MB] > mbRK ? (y[I_MB] - mbRK) / h : 0;
    const rate = (this.mode === MODE_TWO ? this.mdotB : 0) + splitRate;
    this.hrr = rate > 0 ? rate * this.closure.heatOfReaction(Tq) : 0;
    this.bookkeeping(h, muStart, pStepEnd, dpStepEnd);
    void tStart;
  }

  /** Effective sound speed of the charge for the knock modes (sets knockGamma), m/s. */
  private knockSoundSpeed(): number {
    const cl = this.closure;
    const y = this.y;
    if (this.mode === MODE_TWO) {
      const Vb = y[I_MB] * cl.vb;
      const Vu = y[I_MU] * cl.vu;
      this.knockGamma = cl.gu;
      return twoZoneSoundSpeed(this.p, Vb, 1 / cl.vb, cl.gammaSb, Vu, 1 / cl.vu, cl.gu);
    }
    if (this.mode === MODE_BURNED) {
      this.knockGamma = cl.gammab;
      return Math.sqrt(cl.gammaSb * this.p * cl.vb);
    }
    this.knockGamma = cl.gu;
    return Math.sqrt(cl.gu * this.p * cl.vu);
  }

  /**
   * (∂lnτ/∂T)·T at the end-gas state (T, p) — with 1/τ the integrand of the history integral J
   * whose value at the onset gives the ignition-time spread of a thermally stratified end gas (see
   * onAutoignition). ∂lnτ/∂T by central differences (±1 K) of the delay model.
   */
  private lwSensitivity(T: number, p: number): number {
    if (!(T > 1)) return 0;
    const f = this.fuelTrapped;
    const phi = this.phiTrapped;
    const x = this.xResMole;
    const s = (Math.log(this.delayModel.tau(T + 1, p, phi, f, x)) - Math.log(this.delayModel.tau(T - 1, p, phi, f, x))) / 2;
    const g = s * T;
    return Number.isFinite(g) ? g : 0;
  }

  /**
   * The Livengood–Wu integral runs on the end gas AHEAD of the flame front: from IVC while the knock
   * model is available for the trapped fuel, until onset — or until the front has swept the whole
   * charge (no end gas left: the remaining brush behind the front burns up on τ_b and is not "end
   * gas"; validation round 1 found brush/crevice remainders "autoigniting" late in expansion).
   */
  private lwArmed(): boolean {
    return (
      this.opts.knock &&
      this.knockAvailable &&
      !this.knockOnset &&
      (this.mode === MODE_SINGLE || this.mode === MODE_TWO) &&
      // (with only the end gas autoigniting, the integral stops once the front has swept the chamber;
      // with the brush pockets autoigniting too it runs until burn-out, so the onset moves
      // continuously with the operating point instead of appearing with a finite brush mass)
      (this.opts.knockBrushAutoignition || !(this.flameActive && this.frontAtWalls)) &&
      !this.closedInitEndReached()
    );
  }

  private closedInitEndReached(): boolean {
    return this.finished;
  }

  /** End-gas ignition delay (total, hot ignition) at the current evaluation (unburned zone, trapped mixture), s. */
  private tauTotal(): number {
    return this.delayModel.tau(this.Tu, this.p, this.phiTrapped, this.fuelTrapped, this.xResMole);
  }

  /**
   * Delay governing the ACTIVE Livengood–Wu integral, s: the total τ ('single'), or for the
   * two-stage integral τ₁ (first stage) then τ − τ₁ (second stage, ≥ 0: 0 where the chemistry is
   * single-stage, τ₁ = τ, i.e. hot ignition follows the first-stage crossing immediately).
   */
  private tauNow(): number {
    if (!this.twoStage) return this.tauTotal();
    const m = this.delayModel as TabulatedIgnitionDelay;
    const t1 = m.tauFirstStage(this.Tu, this.p, this.phiTrapped, this.fuelTrapped, this.xResMole);
    if (this.lwStage === 1) return t1;
    const d = this.tauTotal() - t1;
    return d > 0 ? d : 0;
  }

  /** Fraction of the tentative step at which the LW integral reaches 1 (−1 if not). */
  private lwCrossingFraction(h: number, tau1: number): number {
    const r1 = tau1 > 0 ? 1 / tau1 : tau1 <= 0 ? Infinity : NaN;
    const tau0 = this.lwTauPrev;
    const r0 = Number.isNaN(tau0) ? r1 : tau0 > 0 ? 1 / tau0 : Infinity;
    const need = 1 - this.lw.integral;
    if (!(need > 0) || !Number.isFinite(r0) || !Number.isFinite(r1)) return -1;
    // Log-mean quadrature of livengood-wu.ts (ln(1/τ) linear over the step):
    const b = r0 > 0 && r1 > 0 ? Math.log(r1 / r0) : 0;
    const I = Math.abs(b) > 1e-8 ? (h * (r1 - r0)) / b : 0.5 * h * (r0 + r1);
    if (!(I >= need)) return -1;
    let f: number;
    if (Math.abs(b) > 1e-8) f = Math.log1p((need * b) / (r0 * h)) / b;
    else f = need / (0.5 * h * (r0 + r1));
    return f > 0 ? (f < 1 ? f : 1) : 0;
  }

  /**
   * End-gas autoignition (Livengood–Wu integral reached 1): start the end-gas burn-up and the
   * acoustic modes.
   *  - Autoigniting ("end-gas") fraction x_eg = m_u/m: ALL unburned gas — ahead of the front and the
   *    entrained-but-unburned brush pockets, which share the unburned zone's T_u(t), p(t) history
   *    (fixer round 2; round 1: only (m − m_e)/m, the gas ahead of the front).
   *  - Burn-up time τ_ab = τ_e + Δt, with Δt the ignition-time spread of an end gas whose
   *    temperature is stratified by ±ΔT (knockStratificationDT) at the onset, taken from the
   *    Livengood–Wu HISTORY: a parcel δ(t) = ΔT·T(t)/T_on hotter all along its (isentropic)
   *    compression reaches I = 1 earlier by Δt = (ΔT/T_on)·τ_on·|J|, J = ∫(1/τ)(∂lnτ/∂T)T dt
   *    (= τ|∂lnτ/∂T|ΔT of knock.ts autoignitionBurnTime for a constant state). Unlike the local
   *    derivative at the onset state, J does not vanish when the onset sits at the NTC turning
   *    point of τ(T) (validation round 1: τ_ab collapsed to τ_e ≈ 1 µs and MAPO jumped to 60–130 bar).
   *  - Acoustic source region: the planform outside a circle about the plug whose area fraction
   *    equals the autoigniting VOLUME fraction (the head-plane section of the flame sphere covers the
   *    whole head late in the burn although end gas remains near the piston: zero source, MAPO 0),
   *    released sequentially from the periphery inward (options.knockSourceShells).
   */
  private onAutoignition(): void {
    if (this.knockOnset) return;
    const y = this.y;
    const mu = y[I_MU];
    const m = mu + y[I_MB];
    if (!(mu > 0)) return;
    // autoigniting gas: all unburned gas (end gas ahead of the front + brush pockets; see evalClosed),
    // or only the end gas ahead of the front (options.knockBrushAutoignition false; round 1)
    const eg = this.opts.knockBrushAutoignition || this.mode === MODE_SINGLE ? mu : Math.min(mu, m - y[I_ME]);
    if (!(eg > 0)) return;
    this.knockOnset = true;
    this.knockBurning = true;
    this.interrupt = true;
    this.knockOnsetDeg = y[I_TH];
    this.tKnockOnset = this.t;
    this.knockEndGasFraction = eg / m;
    const tau = this.tauTotal(); // total (hot-ignition) delay: J is built on it (two-stage too)
    const Ton = this.Tu;
    this.tauAB = autoignitionBurnTime(tau, Ton > 0 ? this.lwJ / Ton : 0, this.opts.knockStratificationDT, this.opts.knockExcitationTime);
    if (!(this.tauAB > 0)) this.tauAB = this.opts.knockExcitationTime;
    // end-gas planform region with the end-gas volume fraction
    const Vg = this.creviceOn ? this.V - this.Vcr : this.V;
    const vEg = this.mode === MODE_SINGLE ? Vg : eg * this.closure.vu;
    const c = this.spec.sparkPlug.gapCenter;
    const rIn = this.endGasCircleRadius(Math.min(1, vEg / Vg));
    if (this.opts.knockSourceShells > 0) this.knockOsc.setSequentialEndGasRegion(c[0], c[2], rIn, this.opts.knockSourceShells);
    else this.knockOsc.setEndGasRegion(c[0], c[2], rIn);
    this.knockEgOnset = eg;
    this.knockMkOnset = y[I_MK];
    // autoignition before any burned gas (e.g. before the spark): seed the burned zone
    if (this.mode === MODE_SINGLE) this.transferToBurned(1e-9 * mu);
    this.derivValid = false;
    this.refresh();
    this.knockC0 = this.knockSoundSpeed();
    this.knockG0 = this.knockGamma;
    this.knockV0 = this.V;
    this.knockOsc.setBandReference(this.knockC0);
  }

  /**
   * Radius of the planform circle centred at the spark plug whose complement in the bore has the
   * area fraction f (0..1): πR² − lens(r) = f πR² (closed-form circle–circle lens, bisection).
   */
  private endGasCircleRadius(f: number): number {
    const R = 0.5 * this.spec.geometry.bore;
    const g = this.spec.sparkPlug.gapCenter;
    const d = Math.max(Math.hypot(g[0], g[2]), 1e-9);
    const target = (1 - f) * Math.PI * R * R; // lens area inside the circle
    let lo = 0;
    let hi = R + d;
    for (let i = 0; i < 60; i++) {
      const mid = 0.5 * (lo + hi);
      if (lensArea(mid, R, d) < target) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /** Ignition-system split; returns true if the thermodynamic state changed. */
  private ignitionSplit(h: number): boolean {
    const ign = this.ign!;
    const y = this.y;
    const th = y[I_TH];
    const gas = this.gas;
    const closed = this.mode !== MODE_OPEN;
    // gas state at the plug (unburned gas around the gap)
    gas.p = this.p;
    gas.Tu = this.Tu > 0 ? this.Tu : this.T;
    const X = closed ? this.closure.Xu : this.cyl.state.X;
    for (let k = 0; k < NS; k++) gas.X[k] = X[k];
    gas.rhoU = closed ? (this.rhoU > 0 ? this.rhoU : this.mCyl / this.V) : this.cyl.state.rho;
    gas.uPrime = this.uPrime;
    // the kernel's eddy size (HM size/time factors, hand-off radius r ≥ l_I, AGB R_L) is the INTEGRAL
    // scale C_ε·L; its strain rate uses the dissipation length L itself (round 1 passed L as the
    // integral scale: hand-off at 5.3 mm instead of ≈ 2.7 mm; validation round 2)
    gas.integralScale = DISSIPATION_COEFFICIENT * this.L;
    gas.dissipationLength = this.L;
    gas.flowVelocity = meanFlowVelocity(y[I_TK], this.mode === MODE_TWO ? y[I_MU] : this.mCyl);
    const kst = ign.state.kernel.stage;
    const kernelWindow = closed && this.mode !== MODE_BURNED && (th >= this.ignCmd.sparkDeg - 1e-9 || this.sparkFired) && (kst === 'none' || kst === 'kernel');
    if (kernelWindow) {
      const phi = this.phiTrapped;
      const fuel = this.fuelTrapped;
      const Tu = gas.Tu;
      const p = this.p;
      const Tad = this.closure.flameTemperatureHP(this.eqAux, Tu, p);
      const nbHP = this.eqAux.result.nTotal; // mol per kg (element vector per kg)
      let nu = 0;
      for (let k = 0; k < NS; k++) nu += this.closure.nu[k];
      gas.expansionRatio = nbHP > 0 ? (nbHP * Tad) / (nu * Tu) : 1;
      marksteinLengths(fuel, phi, Tu, p, this.xDil, Tad, this.mk);
      gas.SL = this.mk.SL;
      gas.marksteinLength = this.opts.marksteinMultiplier * this.mk.unburned;
      gas.kinematicViscosity = this.nuU;
      gas.flameThickness = gas.SL > 0 ? this.nuU / gas.SL : 0;
      gas.lewisNumber = lewisNumbers(fuel, phi, Tu, p, this.xDil, this.le).deficient;
    } else {
      gas.SL = this.SL;
      gas.expansionRatio = 1;
      gas.marksteinLength = 0;
      gas.flameThickness = 0;
      gas.kinematicViscosity = this.nuU;
    }
    const st = ign.step(h, th, this.ignCmd, gas);
    // The IgnitionSystem starts a new event (kernel, ledgers) at the dwell start; until then its
    // state still describes the previous cycle's spark.
    if (st.switchState === 'closed') this.dwellSeen = true;
    if (!this.dwellSeen) return false;
    if (!this.sparkFired && st.switchState !== 'closed' && th >= this.ignCmd.sparkDeg - 1e-9 && th < this.ignCmd.sparkDeg + 180) this.sparkFired = true;
    let changed = false;
    // electrical energy to the gas minus the kernel's electrode conduction loss over the step
    const dE = st.energyToGas - this.sparkEnergyPrev;
    const dL = st.kernel.energyElectrodeLoss - this.electrodeLossPrev;
    this.sparkEnergyPrev = st.energyToGas;
    this.electrodeLossPrev = st.kernel.energyElectrodeLoss;
    const net = (dE > 0 ? dE : 0) - (dL > 0 ? dL : 0);
    // (1) kernel mass → burned zone FIRST: a burned zone created now must take S_u from the charge
    // before this step's electrical energy is added, otherwise that energy would be booked to the
    // unburned zone (its entropy) instead of the kernel gas that received it.
    let handoff = false;
    if (closed && this.mode !== MODE_BURNED && !this.burnDone) {
      const km = st.kernel.burnedMass;
      if ((st.kernel.stage === 'kernel' || st.kernel.stage === 'handoff') && !this.handedOff) {
        const dm = km - this.kernelMassPrev;
        if (dm > 0) {
          const mu = y[I_MU];
          const take = dm < 0.5 * mu ? dm : 0.5 * mu;
          this.transferToBurned(take);
          changed = true;
        }
        if (km > this.kernelMassPrev) this.kernelMassPrev = km;
      }
      handoff = st.kernel.stage === 'handoff' && !this.handedOff;
    }
    // (2) then the energy. With a burned zone, the net spark energy pulls the burned zone toward the
    // kernel temperature T_k of the ignition model (≈ T_ad: the Herweg–Maly kernel converts spark
    // energy into kernel MASS at T_k, which is already moved to the burned zone and burned) but
    // never beyond it; the remainder heats (or, for a net electrode loss, cools) the unburned
    // charge (S_u += ΔE_u/T_u). Round 1 put it all into the burned zone: a 1e-11 kg burned zone from
    // an early "make" spark received tens of mJ, T_b was pinned at the closure's 5900 K limit and
    // thousands of closure solves failed. U_tot always receives the whole net energy (exact).
    if (net !== 0) {
      if (closed) {
        if (this.mode === MODE_TWO && y[I_MU] > 0) {
          if (changed) this.refresh(); // T_b after the kernel-mass transfer
          const Tk = st.kernel.temperature > 0 ? st.kernel.temperature : this.Tb;
          const cap = y[I_MB] * this.closure.cvb * (Tk - this.Tb); // J to bring T_b to T_k
          // a net LOSS (the kernel's electrode conduction after the discharge) cools the kernel gas =
          // the burned zone, at most down to T_u (round 1 took it from the unburned zone; validation
          // round 2 — magnitude < 0.01 K)
          const capLoss = -y[I_MB] * this.closure.cvb * Math.max(0, this.Tb - this.Tu);
          const toB = net > 0 ? (cap > 0 ? Math.min(net, cap) : 0) : Math.max(net, capLoss);
          const toU = net - toB;
          if (toU !== 0 && this.Tu > 0) y[I_SU] += toU / this.Tu;
        }
        y[I_UT] += net;
      } else y[I_CU] += net;
      this.sparkEnergy += net;
      if (dL > 0) this.electrodeLoss += dL;
      changed = true;
    }
    // (3) hand-off to the entrainment model
    if (handoff) {
      if (changed) this.refresh();
      this.onHandoff(st.kernel.radius, st.kernel.turbulentSpeed);
      // The kernel reached the hand-off radius at a known instant INSIDE this step (the kernel
      // integrator locates it exactly; its growth stopped there): catch up the entrainment and
      // burn-up over the remaining part of the step with the rates at the hand-off state (error
      // O(Δt²) instead of the round-1 O(Δt) stairs of CA50 / knock onset vs spark advance).
      const tHo = ign.gap.firstBreakdownTime + st.kernel.handoffTime;
      let late = ign.time - tHo;
      if (!(late > 0)) late = 0;
      if (late > h) late = h;
      if (late > 0 && this.mode === MODE_TWO) {
        this.refresh();
        const m = y[I_MU] + y[I_MB];
        const dme = this.mdotE * late;
        const dmb = Math.min(this.mdotB * late, 0.5 * y[I_MU]);
        this.transferToBurned(dmb);
        y[I_ME] = Math.min(m, Math.max(y[I_ME] + dme, y[I_MB]));
      }
      changed = true;
    }
    if (st.misfire && !this.handedOff) {
      this.misfire = true;
      if (st.kernel.stage === 'quenched') this.kernelQuenched = true;
    }
    return changed;
  }

  /**
   * Initialise the entrainment state at the kernel hand-off, consistent with the kernel:
   *  - the enflamed volume is at least the chamber-clipped volume of the kernel sphere (the
   *    Herweg–Maly kernel is wall-free), V_e ≥ V_geo(r_k);
   *  - the entrained-but-unburned brush mass makes the burning speed continuous: Keck's
   *    s_b = S_L + μ/(ρ_u A_f τ_b) equals the kernel's S_T,k, i.e. μ₀ = ρ_u A_f (S_T,k − S_L) τ_b
   *    (otherwise s_b drops to S_L at the hand-off — Keck 1982 eq. 4.8 starts from an empty brush).
   * @param rk kernel radius at hand-off, m; @param STk kernel burning velocity, m/s
   */
  private onHandoff(rk: number, STk: number): void {
    this.handedOff = true;
    if (this.mode !== MODE_TWO) return;
    const y = this.y;
    y[I_ME] = this.handoffEntrainedMass(rk, STk);
    this.rfGuess = rk;
    this.flameActive = true;
    this.derivValid = false;
  }

  /**
   * Entrained mass m_e that the entrainment model starts from when the kernel of radius rk and
   * burning velocity STk is handed off now (see onHandoff); written to this.fgR's radius as a side
   * product (kernel-stage flame radius, flameRadius()).
   */
  private handoffEntrainedMass(rk: number, STk: number): number {
    const y = this.y;
    const cl = this.closure;
    const mb = y[I_MB];
    const Vb = mb * cl.vb;
    const fg = this.fgR;
    this.flameGeom.evaluate(rk, this.h, fg);
    const sc = this.discScale > 0 ? this.discScale : 1;
    const Vk = fg.volume / sc; // cylinder volume represented by the clipped kernel sphere
    let me = Vk > Vb ? mb + (Vk - Vb) / cl.vu : mb;
    const tauB = this.SL > 0 && this.lambda > 0 && this.lambda < Infinity ? this.lambda / this.SL : Infinity;
    if (STk > this.SL && tauB < Infinity) {
      const rfb = this.flameGeom.radiusForVolume((Vb + (me - mb) * cl.vu) * sc, this.h, this.clampRadiusGuess(rk, this.h));
      this.flameGeom.evaluate(rfb, this.h, fg);
      const mu0 = (fg.frontArea / cl.vu) * (STk - this.SL) * tauB;
      if (mb + mu0 > me) me = mb + mu0;
    }
    const mt = y[I_MU] + mb;
    return me < mt ? me : mt;
  }

  /** Move Δm from the unburned to the burned zone at constant U_tot (creates the burned zone if needed). */
  private transferToBurned(dm: number): void {
    const y = this.y;
    if (!(dm > 0)) return;
    if (this.mode === MODE_SINGLE) this.createBurnedZone();
    if (this.mode !== MODE_TWO) return;
    const mu = y[I_MU];
    if (dm >= mu) {
      this.mergeToBurned();
      return;
    }
    const su = y[I_SU] / mu;
    y[I_SU] -= su * dm;
    y[I_MU] = mu - dm;
    y[I_MB] += dm;
    if (y[I_ME] < y[I_MB]) y[I_ME] = y[I_MB];
    this.derivValid = false;
  }

  /** Single zone → two zones: S_u from the current single-zone state; T_b seeded at the HP flame state. */
  private createBurnedZone(): void {
    const y = this.y;
    const cl = this.closure;
    const mu = y[I_MU];
    const T = cl.solveSingle(y[I_UT], mu, this.creviceOn ? this.V - this.Vcr : this.V, this.TuGuess);
    const p = cl.p;
    y[I_SU] = mu * cl.entropyAt(T, p);
    y[I_MB] = 0;
    this.mode = MODE_TWO;
    this.pGuess = p;
    this.TuGuess = T;
    this.TbGuess = cl.flameTemperatureHP(this.eqAux, T, p);
    this.dpEval = 0;
    this.dTbEval = 0;
    this.tEval = this.t;
    this.Lref = integralLengthScale(this.h, this.spec.geometry.bore, this.lengthFraction);
    this.rhoRef = p / (cl.Ru * T);
    this.rbGuess = 0;
    this.derivValid = false;
  }

  /** Two zones → burned only (the unburned remainder is burned at constant U_tot). */
  private mergeToBurned(): void {
    const y = this.y;
    if (this.mode === MODE_SINGLE) this.createBurnedZone();
    const mu = y[I_MU];
    const m = mu + y[I_MB];
    if (mu > 0) {
      this.suLast = y[I_SU] / mu;
      // turbulence: continue on the whole charge with the unburned specific values
      y[I_TKE] *= m / mu;
      y[I_TK] *= m / mu;
      this.nNO += mu * this.closure.noPerKg;
      this.noMergedMass += mu;
    }
    y[I_MB] = m;
    y[I_MU] = 0;
    y[I_SU] = 0;
    y[I_ME] = m;
    this.mode = MODE_BURNED;
    this.burnDone = true;
    this.flameActive = false;
    this.knockBurning = false;
    this.derivValid = false;
  }

  /** Current flame radius for output (kernel / entrained front / burned sphere / chamber). */
  flameRadius(): number {
    if (this.burnDone) return this.flameGeom.maxRadius(this.h);
    if (this.flameActive) return this.rf;
    // kernel stage: the front radius the entrainment model would start from if the kernel were
    // handed off now (the chamber-clipped kernel sphere plus the brush that makes the burning speed
    // continuous, onHandoff), so the reported radius is continuous through the hand-off (validation
    // round 2: reporting the burned-gas radius jumped 0.5 mm = 7 % at the hand-off)
    if (this.mode === MODE_TWO) {
      const k = this.ign && this.dwellSeen ? this.ign.state.kernel : null;
      if (k && k.stage === 'kernel' && k.radius > 0 && this.rhoU > 0) {
        const me = this.handoffEntrainedMass(k.radius, k.turbulentSpeed);
        const cl = this.closure;
        const sc = this.discScale > 0 ? this.discScale : 1;
        const Ve = this.y[I_MB] * cl.vb + (me - this.y[I_MB]) * cl.vu;
        return this.flameGeom.radiusForVolume(Ve * sc, this.h, this.clampRadiusGuess(k.radius, this.h));
      }
      return this.rb;
    }
    if (this.ign && this.dwellSeen && this.ign.state.kernel.stage === 'kernel') return this.ign.state.kernel.radius;
    return 0;
  }

  /** True once the entrained front has swept the chamber (only the brush is still burning up). */
  get frontAtWalls(): boolean {
    if (!this.flameActive) return false;
    const m = this.y[I_MU] + this.y[I_MB];
    return this.y[I_ME] >= m * (1 - 1e-6) || this.rf >= this.flameGeom.maxRadius(this.h) * (1 - 1e-6) || this.Af <= 1e-10;
  }

  /**
   * Per-step cycle accumulators. Peak pressure and max dp/dθ are properties of the continuous
   * thermodynamic pressure, independent of where the integrator's steps (or the snapshot
   * instants) fall: over each step p(θ) is the cubic Hermite interpolant of p and the analytic
   * dp/dθ (closure rates) at both ends (linear during gas exchange, where the rate is not
   * evaluated). The peak is the interpolant's maximum; max dp/dθ is the largest 0.1°-window
   * secant of the interpolant (DP_WINDOW_DEG, window position on a 0.01° grid).
   * (Round 1 sampled step ends / windows anchored at step ends: max dp/dθ moved 9–20 % with the
   * step size and differed between runCycles and the snapshot-driven worker.)
   */
  private bookkeeping(h: number, muStart: number, pEnd: number, dpdtEnd: number): void {
    const y = this.y;
    const th = y[I_TH];
    const pRep = this.p + this.knockOscillation();
    const omDeg = y[I_OM] * RAD2DEG;
    const th0 = this.y0[I_TH];
    const closedStep = this.mode !== MODE_OPEN && this.hMode !== MODE_OPEN;
    // interpolant on [θ0, θ1] (u ∈ [0, 1]): p0, p1 and slopes a = dp/du at both ends
    const p0 = this.hP;
    const dth = th - th0;
    const valid = !Number.isNaN(p0) && dth > 0 && Math.abs(this.hTheta - th0) < 1e-9;
    let a0 = 0;
    let a1 = 0;
    if (valid) {
      if (closedStep) {
        a0 = this.hDp * dth;
        a1 = (dpdtEnd / omDeg) * dth;
      } else {
        a0 = pEnd - p0;
        a1 = a0;
      }
    }
    // ---- peak of the thermodynamic pressure ----
    if (valid && closedStep) {
      const pk = hermiteMax(p0, a0, pEnd, a1);
      if (pk > this.peakP) {
        this.peakP = pk;
        this.peakPDeg = th0 + HSCR[0] * dth;
      }
    }
    if (pEnd > this.peakP) {
      this.peakP = pEnd;
      this.peakPDeg = th;
    }
    if (this.p > this.peakP) {
      this.peakP = this.p;
      this.peakPDeg = th;
    }
    // ---- max dp/dθ: max over the window position of the secant over DP_WINDOW_DEG, the window
    // sliding on a fixed fine grid (DP_WINDOW_DEG/DP_SUBDIV) anchored at −360° ----
    if (valid) {
      const G = DP_WINDOW_DEG / DP_SUBDIV;
      const buf = this.dpBuf;
      let k = Math.floor((th0 + 360) / G + 1e-9) + 1;
      for (;;) {
        const tk = -360 + k * G;
        if (tk > th + 1e-9) break;
        const u = (tk - th0) / dth;
        const pk = hermiteAt(p0, a0, pEnd, a1, u < 1 ? u : 1);
        if (this.dpGridK !== k - 1) this.dpGridN = 0; // (gap: restart the window)
        const slot = k % DP_SUBDIV;
        if (this.dpGridN >= DP_SUBDIV) {
          const r = (pk - buf[slot]) / DP_WINDOW_DEG;
          if (r > this.maxDp) this.maxDp = r;
        } else this.dpGridN++;
        buf[slot] = pk;
        this.dpGridK = k;
        k++;
      }
      // windows that START at a slope discontinuity (a split that raised dp/dθ: knock onset,
      // hand-off) — the max secant then begins exactly there, which the fixed grid would miss by up
      // to one grid spacing times the slope jump
      for (let i = 0; i < this.dpKinkN; i++) {
        const te = this.dpKinkTh[i] + DP_WINDOW_DEG;
        if (te > th0 && te <= th + 1e-12) {
          const r = (hermiteAt(p0, a0, pEnd, a1, (te - th0) / dth) - this.dpKinkP[i]) / DP_WINDOW_DEG;
          if (r > this.maxDp) this.maxDp = r;
          this.dpKinkTh[i] = this.dpKinkTh[--this.dpKinkN];
          this.dpKinkP[i] = this.dpKinkP[this.dpKinkN];
          i--;
        } else if (te <= th0) {
          this.dpKinkTh[i] = this.dpKinkTh[--this.dpKinkN];
          this.dpKinkP[i] = this.dpKinkP[this.dpKinkN];
          i--;
        }
      }
    }
    if (this.mode !== MODE_OPEN && this.dpdt > dpdtEnd + 1e-3 * Math.abs(dpdtEnd) && this.dpKinkN < this.dpKinkTh.length) {
      this.dpKinkTh[this.dpKinkN] = th;
      this.dpKinkP[this.dpKinkN] = this.p;
      this.dpKinkN++;
    }
    // start of the next step: the post-split state
    this.hTheta = th;
    this.hP = this.p;
    this.hDp = this.mode === MODE_OPEN ? 0 : this.dpdt / omDeg;
    this.hMode = this.mode;
    // CA10/50/90 (linear interpolation of x_b inside the step)
    if (this.mode !== MODE_OPEN) {
      const m = y[I_MU] + y[I_MB];
      const x1 = y[I_MB] / m;
      const x0 = 1 - muStart / m;
      if (Number.isNaN(this.ca10)) this.ca10 = crossing(x0, x1, th0, th, 0.1);
      if (Number.isNaN(this.ca50)) this.ca50 = crossing(x0, x1, th0, th, 0.5);
      if (Number.isNaN(this.ca90)) this.ca90 = crossing(x0, x1, th0, th, 0.9);
    }
    if (this.trace && this.cycle === this.traceCycle) this.recordSample(pRep);
    void h;
  }

  private recordSample(pRep: number): void {
    const tr = this.trace!;
    const y = this.y;
    const closed = this.mode !== MODE_OPEN;
    const mu = closed ? y[I_MU] : this.mCyl;
    const mb = closed ? y[I_MB] : 0;
    tr.theta.push(y[I_TH]);
    tr.t.push(this.t);
    tr.volume.push(this.V);
    tr.pressure.push(this.p);
    tr.pressureReported.push(pRep);
    tr.Tu.push(this.Tu);
    tr.Tb.push(this.Tb);
    tr.Tmean.push(this.T);
    tr.mu.push(mu);
    tr.mb.push(mb);
    tr.xb.push(closed ? mb / (mu + mb) : 0);
    tr.me.push(closed ? y[I_ME] : 0);
    tr.U.push(closed ? y[I_UT] + this.creviceEnergy() : y[I_CU]);
    tr.work.push(y[I_W]);
    tr.heatLoss.push(y[I_Q]);
    tr.sparkEnergy.push(this.sparkEnergy);
    tr.SL.push(this.SL);
    tr.uPrime.push(this.uPrime);
    tr.flameRadius.push(this.flameRadius());
    tr.frontArea.push(this.Af);
    tr.burnRate.push(this.burnRateStep);
    tr.heatReleaseRate.push(this.hrr);
    tr.heatLossRate.push(this.Qwall);
    tr.lwIntegral.push(this.lw.integral);
    tr.xNO.push(this.burnedNOFraction());
    tr.nNOClosure.push(this.noCoupled && closed ? this.closure.noKinetic : 0);
    tr.mdotIntake.push(this.mdotIv);
    tr.mdotExhaust.push(this.mdotEv);
    tr.pIntake.push(this.pInt);
    tr.pExhaust.push(this.pExh);
    tr.mass.push(this.mCyl);
    tr.mode.push(this.mode);
    tr.L.push(this.L);
    tr.lambda.push(this.lambda);
  }

  /** Synthesised knock pressure oscillation at the pickup, Pa (0 before onset / with open valves). */
  knockOscillation(): number {
    return this.knockOnset && this.mode !== MODE_OPEN ? this.knockOsc.sensorPressure() : 0;
  }

  /** Kinetic NO mole fraction of the burned zone (0 without one). */
  burnedNOFraction(): number {
    const nb = this.closure.eq.result.nTotal;
    return (this.mode === MODE_TWO || this.mode === MODE_BURNED) && nb > 0 ? this.nNO / nb : 0;
  }

  // ===========================================================================================
  // Events
  // ===========================================================================================

  private handleEvents(th: number): void {
    while (this.evNext < this.nEv && this.evAngle[this.evNext] <= th + 1e-9) {
      const k = this.evKind[this.evNext];
      this.evNext++;
      switch (k) {
        case EV_IVC:
          if (this.mode === MODE_OPEN) this.onIvc();
          break;
        case EV_EVO:
          if (this.mode !== MODE_OPEN) this.onEvo();
          break;
        case EV_BDC_START:
          this.wBdcStart = this.y[I_W];
          break;
        case EV_BDC_END:
          this.wBdcEnd = this.y[I_W];
          break;
        case EV_TDC:
          if (this.mode === MODE_SINGLE || this.mode === MODE_TWO) {
            const x0 = this.xb;
            this.mergeToBurned();
            this.refresh();
            this.markBurnJump(x0, this.xb, th);
          }
          break;
        case EV_IGN:
          if (Math.abs(this.evAngle[this.evNext - 1] - this.ignCmd.sparkDeg) < 1e-9 && Number.isNaN(this.tSparkCmd)) {
            this.tSparkCmd = this.t;
            this.interrupt = true;
          }
          break;
        case EV_WIEBE:
          if (this.mode === MODE_SINGLE) {
            // Seed of the burned zone: WIEBE_SEED of the charge (numerical; x_b starts at it). The
            // two-zone closure recovers T_b from U_tot − m_u u_u, so dT_b/dp ≈ (V_u/γ)/(m_b c_v): with
            // the round-1 seed 1e-8 (7e-12 kg) p had to be exact to 1e-4 Pa and 2–3 solves failed per
            // cycle (T_b pinned at a bound; validation round 2). At 1e-5 the requirement is 1e-9 of p.
            const x0 = this.xb;
            this.transferToBurned(WIEBE_SEED * this.y[I_MU]);
            this.refresh();
            this.markBurnJump(x0, this.xb, th);
          }
          break;
        case EV_END:
          this.onEvo();
          this.finished = true;
          break;
        case EV_WRAP:
          this.onWrap();
          return;
        default:
          break;
      }
    }
  }

  /**
   * A discrete event moved x_b from x0 to x1 at angle th (instantaneous burn at TDC): CA levels it
   * crossed are located at the event (bookkeeping only finds crossings INSIDE steps; validation round 2
   * found CA10/50/90 = NaN for 'instantaneous-at-tdc').
   */
  private markBurnJump(x0: number, x1: number, th: number): void {
    if (Number.isNaN(this.ca10) && x0 < 0.1 && x1 >= 0.1) this.ca10 = th;
    if (Number.isNaN(this.ca50) && x0 < 0.5 && x1 >= 0.5) this.ca50 = th;
    if (Number.isNaN(this.ca90) && x0 < 0.9 && x1 >= 0.9) this.ca90 = th;
  }

  /** IVC: the cylinder zone becomes the closed charge (frozen unburned composition). */
  private onIvc(): void {
    const y = this.y;
    const N = this.Nscr;
    let n = 0;
    let m = 0;
    let mf = 0;
    for (let k = 0; k < NS; k++) {
      N[k] = y[I_CN + k];
      n += N[k];
      m += N[k] * MOLAR_MASS[k];
    }
    for (const s of FUEL_SPECIES) mf += N[SP[s]] * MOLAR_MASS[SP[s]];
    for (let k = 0; k < NS; k++) this.Xscr[k] = N[k] / n;
    const cl = this.closure;
    cl.setUnburned(this.Xscr);
    const U = y[I_CU];
    y[I_UT] = U;
    y[I_MU] = m;
    y[I_MB] = 0;
    y[I_ME] = 0;
    y[I_SU] = 0;
    y[I_CRU] = 0;
    y[I_CRB] = 0;
    this.mode = MODE_SINGLE;
    const V = this.kin.volume(y[I_TH] * DEG);
    let T = cl.solveSingle(U, m, V, this.T > 0 ? this.T : 400);
    let Vg = V;
    if (this.creviceOn) {
      // crevice zone at the mean piston/liner temperature (the top land lies between them); per-kg
      // properties of the unburned charge and of its complete-combustion products at T_cr
      this.Tcr = 0.5 * (this.walls.pistonTemperature + this.walls.linerTemperature);
      this.uuCr = mixUMolar(cl.nu, this.Tcr);
      this.huCr = mixHMolar(cl.nu, this.Tcr);
      this.uccCr = mixUMolar(cl.ncc, this.Tcr);
      this.hccCr = mixHMolar(cl.ncc, this.Tcr);
      let ncc = 0;
      for (let k = 0; k < NS; k++) ncc += cl.ncc[k];
      this.Rcc = R_UNIVERSAL * ncc;
      // split the trapped charge: m_cr = p V_cr/(R_u T_cr) at the IVC pressure (fixed point; the gas
      // moved into the crevice takes u(T_cr) with it — an adiabatic re-partition, energy exact)
      Vg = V - this.Vcr;
      let mcr = 0;
      for (let it = 0; it < 6; it++) {
        T = cl.solveSingle(U - mcr * this.uuCr, m - mcr, Vg, T);
        mcr = (cl.p * this.Vcr) / (cl.Ru * this.Tcr);
      }
      T = cl.solveSingle(U - mcr * this.uuCr, m - mcr, Vg, T);
      y[I_UT] = U - mcr * this.uuCr;
      y[I_MU] = m - mcr;
      y[I_CRU] = mcr;
    }
    this.TuGuess = T;
    this.pIvc = cl.p;
    this.TIvc = T;
    this.VIvc = Vg;
    this.mIvc = m;
    this.fuelMassIvc = mf;
    this.TmotGuess = T;
    // mixture descriptors of the trapped charge (fuel blend, φ) — see fuelTrapped
    this.trappedMixture(N, mf);
    // Keck (1982) inlet state: ρ_i = intake-plenum density now, ū_i = ε_v (A_p/A_IV) 2NS with
    // ε_v = trapped fresh mass/(ρ_i V_d) and A_IV the flow area of the inlet valve at maximum lift
    this.rhoInlet = this.intake.state.rho;
    {
      const yb = m > 0 ? Math.min(1, Math.max(0, (y[I_CBG] - y[I_CEG]) / m)) : 0;
      const ev = (m * (1 - yb)) / (this.rhoInlet * this.kin.displacedVolume);
      const Aiv = this.ivFlow.flowArea(this.spec.intakeValve.maxLift);
      this.keckInletSpeed = keckMeanInletSpeed(ev, this.kin.boreArea, Aiv, (y[I_OM] * 60) / (2 * Math.PI), this.spec.geometry.stroke);
    }
    // residual / dilution
    // burned-gas (products) fraction of the TRAPPED charge — residual + external EGR — is the dilution
    // of S_L and the ignition delay (round 1 added the REQUESTED op.egrFraction: 2.2× too much in the
    // first cycle after an EGR step, which then misfired; validation round 2)
    this.yRes = m > 0 ? Math.min(1, Math.max(0, y[I_CBG] / m)) : 0;
    this.yEgr = m > 0 ? Math.min(this.yRes, Math.max(0, y[I_CEG] / m)) : 0;
    this.xDil = this.yRes;
    const fa = freshCharge({ fuel: this.fuelTrapped, phi: this.phiTrapped > 0 ? this.phiTrapped : 1, airX: this.airX });
    const Mfa = mixMolarMass(fa);
    const Mres = mixMolarMass(completeCombustionProducts(fa));
    this.xResMole = residualMoleFraction(this.xDil, Mres, Mfa);
    // unburned viscosity table
    for (let i = 0; i < this.muTab.length; i++) this.muTab[i] = mixViscosity(cl.Xu, this.muT0 + i * this.muDT);
    this.nuU = (this.viscosityU(T) * V) / m;
    // knock / NO / flame state
    this.lw.reset();
    this.lwStage = 1;
    this.twoStage = this.opts.knockIntegral === 'two-stage' && this.delayModel instanceof TabulatedIgnitionDelay;
    this.lwTauPrev = NaN;
    this.lwJ = 0;
    this.lwRPrev = NaN;
    this.lwGPrev = 0;
    this.knockOsc.reset();
    this.nNO = 0;
    this.qIvc = y[I_Q];
    this.closedHappened = true;
    this.rfGuess = 0;
    this.rbGuess = 0;
    this.derivValid = false;
    this.refresh();
    // is the knock model defined for this fuel? (PRF delay models: NaN / RangeError otherwise)
    let tau0 = NaN;
    try {
      tau0 = this.tauNow();
    } catch {
      tau0 = NaN;
    }
    this.knockAvailable = this.fuelMassIvc > 0 && (Number.isFinite(tau0) || tau0 === Infinity);
    this.lwTauPrev = this.knockAvailable ? tau0 : NaN;
  }

  /**
   * Fuel blend and fresh-charge equivalence ratio of the trapped charge (moles N, fuel mass mf):
   *  - fuel: the normalised fuel-species vapour composition (the selected FuelBlend object when it
   *    has the same composition, so its tables/caches are reused);
   *  - φ from the element balance, exact for any mix of fresh charge, residual and EGR of any φ
   *    (combustion conserves elements; H2O and CO2 of the air are neutral in D):
   *      D = 2 b_C + b_H/2 − b_O = 2 n_O2,air (φ − 1),  n_O2,air = (b_N/2)·(x_O2/x_N2)_air.
   */
  private trappedMixture(N: Float64Array, mf: number): void {
    const sel = this.fuel;
    let nf = 0;
    for (const sName of FUEL_SPECIES) nf += N[SP[sName]];
    if (!(mf > 0) || !(nf > 0)) {
      this.fuelTrapped = sel;
      this.phiTrapped = 0;
      return;
    }
    let same = true;
    for (const sName of FUEL_SPECIES) {
      const k = SP[sName];
      if (Math.abs(N[k] / nf - sel.X[k]) > 1e-9) same = false;
    }
    if (same) this.fuelTrapped = sel;
    else {
      const X = new Float64Array(NS);
      for (const sName of FUEL_SPECIES) X[SP[sName]] = N[SP[sName]] / nf;
      this.fuelTrapped = { label: `${sel.label} (trapped mix)`, X };
    }
    const b = elementMolesOf(N, this.bScr);
    const rO2N2 = this.airX[SP.O2] / this.airX[SP.N2];
    const nO2 = 0.5 * b[EL.N] * rO2N2;
    this.phiTrapped = nO2 > 0 ? Math.max(0, 1 + (2 * b[EL.C] + 0.5 * b[EL.H] - b[EL.O]) / (2 * nO2)) : 0;
  }

  /**
   * Mass fraction of DRY air in the fresh charge (fuel vapour + humid air + EGR, freshCharge):
   * 1 − Y_fuel − Y_EGR − Y_water,air, for the volumetric efficiency.
   */
  private freshMassFractions(): { dryAir: number } {
    const X = this.Xfresh;
    let m = 0;
    let mf = 0;
    for (let k = 0; k < NS; k++) m += X[k] * MOLAR_MASS[k];
    for (const s of FUEL_SPECIES) mf += X[SP[s]] * MOLAR_MASS[SP[s]];
    const humidAirFraction = 1 - mf / m - this.op.egrFraction;
    const xw = this.airX[SP.H2O];
    const yWaterInAir = (xw * MOLAR_MASS[SP.H2O]) / mixMolarMass(this.airX);
    return { dryAir: Math.max(0, humidAirFraction * (1 - yWaterInAir)) };
  }

  /** EVO: merge the zones (burned composition frozen, NO → kinetic value) into the gas-exchange zone. */
  private onEvo(): void {
    const y = this.y;
    const cl = this.closure;
    const mu = y[I_MU];
    const mb = y[I_MB];
    const m = mu + mb;
    const N = this.Nscr;
    for (let k = 0; k < NS; k++) N[k] = mu * cl.nu[k];
    // crevice gas back into the gas-exchange zone (unburned: frozen charge; burned: complete-combustion
    // products), with its energy at T_cr
    const mcu = y[I_CRU];
    const mcb = y[I_CRB];
    if (mcu > 0 || mcb > 0) for (let k = 0; k < NS; k++) N[k] += mcu * cl.nu[k] + mcb * cl.ncc[k];
    if (mb > 0 && this.mode !== MODE_SINGLE) {
      // burned zone at its equilibrium state (T_b, p), NO replaced by the kinetic value
      const bb = this.bScr;
      for (let e = 0; e < NE; e++) bb[e] = cl.bu[e] * mb;
      const r = cl.eq.solveTP(bb, this.Tb, this.p);
      const Nb = this.Nb;
      for (let k = 0; k < NS; k++) Nb[k] = r.N[k];
      swapNO(Nb, this.nNO);
      for (let k = 0; k < NS; k++) N[k] += Nb[k];
    }
    let n = 0;
    for (let k = 0; k < NS; k++) {
      y[I_CN + k] = N[k];
      n += N[k];
    }
    y[I_CU] = y[I_UT] + this.creviceEnergy();
    y[I_CBG] = mb + mcb + (mu + mcu) * this.yRes;
    y[I_CEG] = 0; // after the cycle everything in the cylinder is residual
    y[I_CRU] = 0;
    y[I_CRB] = 0;
    // turbulence on the whole charge
    if (this.mode === MODE_TWO && mu > 0) {
      y[I_TKE] *= m / mu;
      y[I_TK] *= m / mu;
    }
    this.noPpm = (N[SP.NO] / n) * 1e6;
    this.coFrac = N[SP.CO] / n;
    this.xbEvo = m > 0 ? mb / m : 0;
    if (mb > 0) this.burnDone = this.burnDone || this.xbEvo > 0.999;
    this.qEvo = y[I_Q];
    this.knockOsc.reset();
    this.lw.reset();
    this.mode = MODE_OPEN;
    this.flameActive = false;
    this.knockBurning = false;
    this.stiffness = 0;
    this.derivValid = false;
    this.refresh();
  }

  /** θ = 360: close the cycle summary, apply the pending operating point, θ → −360. */
  private onWrap(): void {
    const y = this.y;
    this.summaries.push(this.makeSummary());
    this.updateWalls(this.t - this.tCycleStart, this.warmingUp);
    for (let i = 0; i < 5; i++) this.qsCycleStart[i] = y[I_QS + i];
    this.tCycleStart = this.t;
    // new cycle
    const crOld = this.op.compressionRatio;
    const rpmNow = this.op.rpm;
    const mode = this.op.speedMode;
    this.op = cloneOp(this.pendingOp);
    if (mode === 'free' && this.op.speedMode === 'free') this.op.rpm = rpmNow; // keep ω (state)
    if (this.op.compressionRatio !== crOld) this.changeCompressionRatio(this.op.compressionRatio);
    this.prepareCycleData();
    y[I_TH] = -360;
    if (this.op.speedMode === 'fixed') y[I_OM] = (this.op.rpm * 2 * Math.PI) / 60;
    this.cycle++;
    this.resetCycleFlags();
    this.startCycle();
    this.wrapped = true;
    this.derivValid = false;
    this.refresh();
  }

  /**
   * Compression-ratio change at the cycle wrap (θ = −360, gas-exchange TDC). The CFR raises the
   * cylinder over seconds while running; here the clearance volume changes at once, so the cylinder
   * gas is compressed/expanded ISENTROPICALLY to the new volume and the work is booked in the ∫p dV
   * ledger (between cycles: in neither summary). Round 1 kept U at the new volume (no work): p jumped
   * 1.013 → 1.306 bar for CR 6.43 → 8 while the entropy DROPPED — a second-law violation (validation
   * round 2). Newton on T: s(T, nRT/V₂) = s₁, (∂s/∂T)_v = c_v/T (molar, frozen composition).
   */
  private changeCompressionRatio(cr: number): void {
    const y = this.y;
    const N = this.Nscr;
    const X = this.Xscr;
    let n = 0;
    for (let k = 0; k < NS; k++) {
      N[k] = y[I_CN + k];
      n += N[k];
    }
    const V1 = this.kin.volume(y[I_TH] * DEG);
    this.kin.setCompressionRatio(cr);
    const V2 = this.kin.volume(y[I_TH] * DEG);
    if (!(n > 0) || this.mode !== MODE_OPEN || !(V1 > 0) || Math.abs(V2 / V1 - 1) < 1e-15) return;
    for (let k = 0; k < NS; k++) X[k] = N[k] / n;
    const U1 = y[I_CU];
    const T1 = temperatureFromUMolar(N, U1, this.T > 0 ? this.T : 800);
    const s1 = mixSMolar(X, T1, (n * R_UNIVERSAL * T1) / V1);
    let T = T1 * Math.pow(V1 / V2, 0.3);
    for (let it = 0; it < 50; it++) {
      const f = mixSMolar(X, T, (n * R_UNIVERSAL * T) / V2) - s1;
      const dT = -f / (mixCvMolar(X, T) / T);
      T += dT;
      if (Math.abs(dT) < 1e-12 * T) break;
    }
    const U2 = n * mixUMolar(X, T);
    y[I_CU] = U2;
    y[I_W] -= U2 - U1; // ∫p dV of the isentropic change (−W_on_gas)
    this.cyl.volume = V2;
  }

  /** Net mass through the carburettor venturi since the start of the run (ledger I_LV), kg. */
  private ventMass(): number {
    let m = 0;
    for (let k = 0; k < NS; k++) m += this.y[I_LV + k] * MOLAR_MASS[k];
    return m;
  }

  /**
   * Volumetric efficiency of the cycle just completed (Heywood's definition, as measured on the
   * CFR): DRY air inducted through the carburettor over the cycle / (ambient dry-air partial
   * density × V_d). Net venturi mass over the cycle × dry-air mass fraction of the fresh stream.
   * (Round 1 used the trapped mass × (1 − burned-residual fraction), which counted retained
   * UNBURNED gas — motored / misfiring cycles — as inducted: η_v ≈ 1.09.)
   */
  private volumetricEfficiency(): number {
    const mVent = this.ventMass() - this.mVentCycleStart;
    const Yf = this.freshMassFractions();
    const xw = this.airX[SP.H2O];
    // dry-air partial density of the ambient (humid) air, kg/m³ (0 if the "air" is saturated steam)
    const rhoDry = xw < 1 ? ((mixMolarMass(this.airX) - xw * MOLAR_MASS[SP.H2O]) * this.op.ambientPressure) / (R_UNIVERSAL * this.op.ambientTemperature) : 0;
    return rhoDry > 0 ? (mVent * Yf.dryAir) / (rhoDry * this.kin.displacedVolume) : 0;
  }

  private makeSummary(): CycleSummary {
    const y = this.y;
    const Vd = this.kin.displacedVolume;
    const wNet = y[I_W] - this.wCycleStart;
    const wGross = this.wBdcEnd - this.wBdcStart;
    const lhv = lowerHeatingValue(this.fuelTrapped);
    this.volEff = this.closedInit ? 0 : this.volumetricEfficiency();
    const mf = this.fuelMassIvc;
    const burned = this.xbEvo;
    const misfire =
      this.opts.combustionModel === 'entrainment' ? this.misfire || burned < 0.1 : this.opts.combustionModel === 'none' ? true : burned < 0.1;
    return {
      cycle: this.cycle,
      imepGross: wGross / Vd,
      imepNet: wNet / Vd,
      pmep: (wGross - wNet) / Vd,
      peakPressure: this.peakP,
      peakPressureDeg: this.peakPDeg,
      maxPressureRiseRate: this.maxDp,
      ca10: this.ca10,
      ca50: this.ca50,
      ca90: this.ca90,
      indicatedEfficiency: mf > 0 ? wNet / (mf * lhv) : 0,
      isfc: wNet > 0 && mf > 0 ? mf / wNet : NaN,
      trappedMass: this.mIvc,
      residualFraction: this.yRes - this.yEgr,
      volumetricEfficiency: this.volEff,
      fuelMass: mf * burned,
      noPpm: this.noPpm,
      coFraction: this.coFrac,
      knockOnsetDeg: this.knockOnsetDeg,
      knockEndGasFraction: this.knockOnset ? this.knockEndGasFraction : 0,
      mapo: this.mapo,
      misfire,
      heatLoss: this.qEvo - this.qIvc,
      indicatedWorkGross: wGross,
    };
  }

  // ===========================================================================================
  // Inventories (conservation checks)
  // ===========================================================================================

  /** Internal energy of the crevice zone (closed phase), J. */
  creviceEnergy(): number {
    const y = this.y;
    return this.mode === MODE_OPEN ? 0 : y[I_CRU] * this.uuCr + y[I_CRB] * this.uccCr;
  }

  /** Moles of each species in cylinder + plenums (closed phase: unburned + burned equilibrium + crevice), mol. */
  speciesInventory(out: Float64Array = new Float64Array(NS)): Float64Array {
    const y = this.y;
    for (let k = 0; k < NS; k++) out[k] = y[I_IN + k] + y[I_EN + k];
    if (this.mode === MODE_OPEN) {
      for (let k = 0; k < NS; k++) out[k] += y[I_CN + k];
    } else {
      const cl = this.closure;
      for (let k = 0; k < NS; k++) out[k] += (y[I_MU] + y[I_CRU]) * cl.nu[k] + y[I_CRB] * cl.ncc[k];
      if (y[I_MB] > 0) {
        const r = cl.eq.result;
        const s = r.mass > 0 ? y[I_MB] / r.mass : 0;
        for (let k = 0; k < NS; k++) out[k] += r.N[k] * s;
      }
    }
    return out;
  }

  /** Element moles in cylinder + plenums (exact in every mode), mol. */
  elementInventory(out: Float64Array = new Float64Array(NE)): Float64Array {
    const y = this.y;
    const N = this.Nscr;
    for (let k = 0; k < NS; k++) N[k] = y[I_IN + k] + y[I_EN + k] + (this.mode === MODE_OPEN ? y[I_CN + k] : 0);
    const b = elementMolesOf(N, out);
    if (this.mode !== MODE_OPEN) {
      const m = y[I_MU] + y[I_MB] + y[I_CRU] + y[I_CRB];
      for (let e = 0; e < NE; e++) b[e] += m * this.closure.bu[e];
    }
    return b;
  }

  /** Internal energy of cylinder + plenums, J. */
  systemEnergy(): number {
    const y = this.y;
    return y[I_IU] + y[I_EU] + (this.mode === MODE_OPEN ? y[I_CU] : y[I_UT] + this.creviceEnergy());
  }

  /** Mass in cylinder + plenums, kg. */
  systemMass(): number {
    const y = this.y;
    let m = 0;
    for (let k = 0; k < NS; k++) m += (y[I_IN + k] + y[I_EN + k]) * MOLAR_MASS[k];
    if (this.mode === MODE_OPEN) for (let k = 0; k < NS; k++) m += y[I_CN + k] * MOLAR_MASS[k];
    else m += y[I_MU] + y[I_MB] + y[I_CRU] + y[I_CRB];
    return m;
  }
}

// =============================================================================================
// Helpers
// =============================================================================================

/**
 * CycleSummary.maxPressureRiseRate = max over θ of [p(θ + Δ) − p(θ)]/Δ with Δ = DP_WINDOW_DEG
 * (the usual finite-difference definition for 0.1°-sampled indicator data), the window position
 * sampled every Δ/DP_SUBDIV so the result does not depend on where a steep (knock) rise falls
 * relative to a coarse grid.
 */
export const DP_WINDOW_DEG = 0.1;
export const DP_SUBDIV = 10;

/**
 * ∫₀¹ r(u) g(u) du for r exponential (r0 → r1, both ≥ 0) and g linear (g0 → g1): with b = ln(r1/r0),
 * r0 [g0 (e^b − 1)/b + (g1 − g0)(e^b (b − 1) + 1)/b²] (series for small |b|); trapezoid if r0 or r1 is 0.
 */
export function expLinearIntegral(r0: number, r1: number, g0: number, g1: number): number {
  if (!(r0 > 0 && r1 > 0)) return 0.5 * (r0 * g0 + r1 * g1);
  const b = Math.log(r1 / r0);
  const dg = g1 - g0;
  if (Math.abs(b) < 1e-4) return r0 * (g0 * (1 + b / 2 + (b * b) / 6) + dg * (0.5 + b / 3 + (b * b) / 8));
  const eb = Math.exp(b);
  return r0 * (g0 * ((eb - 1) / b) + (dg * (eb * (b - 1) + 1)) / (b * b));
}

/** Scratch: argument u of the last hermiteMax. */
const HSCR = new Float64Array(1);

/** Cubic Hermite P(u) on [0, 1] with P(0) = p0, P′(0) = a0, P(1) = p1, P′(1) = a1. */
function hermiteAt(p0: number, a0: number, p1: number, a1: number, u: number): number {
  const d = p1 - p0;
  const c2 = 3 * d - 2 * a0 - a1;
  const c3 = a0 + a1 - 2 * d;
  return p0 + u * (a0 + u * (c2 + u * c3));
}

/** Maximum of the cubic Hermite interpolant over u ∈ [0, 1] (argument in HSCR[0]). */
function hermiteMax(p0: number, a0: number, p1: number, a1: number): number {
  let best = p0;
  HSCR[0] = 0;
  if (p1 > best) {
    best = p1;
    HSCR[0] = 1;
  }
  const d = p1 - p0;
  const c2 = 3 * d - 2 * a0 - a1;
  const c3 = a0 + a1 - 2 * d;
  // P′(u) = a0 + 2 c2 u + 3 c3 u²
  const A = 3 * c3;
  const B = 2 * c2;
  if (Math.abs(A) <= 1e-14 * (Math.abs(B) + Math.abs(a0))) {
    if (B !== 0) best = hermiteTry(p0, a0, c2, c3, -a0 / B, best);
  } else {
    const disc = B * B - 4 * A * a0;
    if (disc >= 0) {
      const sq = Math.sqrt(disc);
      const q = -0.5 * (B + (B >= 0 ? sq : -sq));
      if (q !== 0) {
        best = hermiteTry(p0, a0, c2, c3, q / A, best);
        best = hermiteTry(p0, a0, c2, c3, a0 / q, best);
      }
    }
  }
  return best;
}

function hermiteTry(p0: number, a0: number, c2: number, c3: number, u: number, best: number): number {
  if (!(u > 0 && u < 1)) return best;
  const v = p0 + u * (a0 + u * (c2 + u * c3));
  if (v > best) {
    HSCR[0] = u;
    return v;
  }
  return best;
}

/** Angle at which x crosses `lev` inside a step (linear in x), NaN if it does not. */
function crossing(x0: number, x1: number, th0: number, th1: number, lev: number): number {
  return x1 >= lev && x0 < lev ? th0 + ((lev - x0) / (x1 - x0)) * (th1 - th0) : NaN;
}

function cloneOp(op: OperatingPoint): OperatingPoint {
  return { ...op, fuel: { ...op.fuel } };
}

function massOf(N: Float64Array): number {
  let m = 0;
  for (let k = 0; k < NS; k++) m += N[k] * MOLAR_MASS[k];
  return m;
}

function elementMolesOf(N: Float64Array, b: Float64Array): Float64Array {
  b.fill(0);
  // ELEMENT_COUNTS via the thermo helper would allocate nothing either; inline for clarity
  for (let k = 0; k < NS; k++) {
    const n = N[k];
    if (n === 0) continue;
    b[EL.C] += n * EC[k * NE + EL.C];
    b[EL.H] += n * EC[k * NE + EL.H];
    b[EL.O] += n * EC[k * NE + EL.O];
    b[EL.N] += n * EC[k * NE + EL.N];
    b[EL.AR] += n * EC[k * NE + EL.AR];
  }
  return b;
}

const EC = new Float64Array(NS * NE);
for (let k = 0; k < NS; k++) for (let e = 0; e < NE; e++) EC[k * NE + e] = ELEMENT_COUNTS[k][e];

/**
 * Replace the NO moles of a burned-gas mole vector by nNO conserving elements (closure.ts
 * noSwapDeltas: N from N2; O from O2, then from CO2 → CO and H2O → H2 in proportion; a decrease
 * returns N2 + O2). The same swap is the kinetic-NO energy correction of the closure.
 */
export function swapNO(N: Float64Array, nNO: number): void {
  noSwapDeltas(N, nNO, SWAP_DN);
  N[SP.NO] += SWAP_DN[SP.NO];
  N[SP.N2] += SWAP_DN[SP.N2];
  N[SP.O2] += SWAP_DN[SP.O2];
  N[SP.CO2] += SWAP_DN[SP.CO2];
  N[SP.CO] += SWAP_DN[SP.CO];
  N[SP.H2O] += SWAP_DN[SP.H2O];
  N[SP.H2] += SWAP_DN[SP.H2];
}
const SWAP_DN = new Float64Array(NS);

export function resolveDelayModel(m: CycleModelOptions['ignitionDelayModel']): IgnitionDelayModel {
  if (typeof m === 'object') return m;
  if (m === 'douaud-eyzat') return douaudEyzat;
  if (m === 'douaud-eyzat-llnl') return douaudEyzatLLNL;
  if (m === 'llnl-prf-v2') return prfLLNLv2;
  return prfLLNLGasoline2011;
}
