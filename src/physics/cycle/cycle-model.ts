/**
 * The cycle model: one quasi-dimensional, energy-conserving model of the whole engine (cylinders,
 * intake and exhaust plenums, spark system, flame, knock, NO, heat transfer, mechanics) advanced
 * in TIME with a classical explicit Runge–Kutta (RK4) scheme that never steps across a discrete
 * event. EngineSimulator (engine-simulator.ts) wraps it for the worker; runClosedCycle and
 * CycleModel.runCycles are the validation hooks (index.ts documents the public API).
 *
 * ── Engine and cylinders (Model T integration) ──────────────────────────────────────────────
 * CycleModel is the ENGINE: time, crank angle θ (cylinder 1's angle) and speed ω, the shared intake
 * and exhaust plenums, carburettor venturi and exhaust outlet, the RK4/Kahan stepping with its event
 * landing, the warm-up, the conservation inventories and the engine summary. Each of the
 * spec.cylinders cylinders (cylinder.ts Cylinder) carries its own state block, sub-models, right-hand
 * side, operator splits, events and summaries at its LOCAL angle θ_i = θ − layout.firingOffsetDeg[i]
 * (wrapped into [−360, 360)). Cylinder 0 keeps the legacy single-cylinder state indices I_* and its
 * fields are re-exposed here (m.p, m.closure, m.knockOsc, …), so a one-cylinder spec (the CFR F-1)
 * runs exactly the former single-cylinder model, bit for bit. Coupling: every cylinder's valve flows
 * enter the shared plenum rates (begun once per evaluation); the step limit is the minimum over the
 * cylinders; splits run per cylinder; each cylinder's local wrap (θ_i = 360) emits its CycleSummary
 * (summary.cylinder = i for N > 1); the ENGINE wrap (θ = 360, cylinder 0's wrap) applies the pending
 * operating point and attaches the EngineCycleSummary (summary.engine). Speed: 'fixed' for any N;
 * 'free' through the rigid multi-cylinder crank train (mechanics MultiCylinderCrankTrain: every
 * cylinder's gas, inertia and gravity torque at its throw phase, Σ J_m), the whole-engine friction
 * (FrictionTorqueModel normalised by N·V_d, torqueCylinders) and the load model (mechanics LoadModel from
 * op.load: constant / brake / vehicle road load with the car's reflected inertia / neutral); for one
 * cylinder and a constant load this is bit for bit the former single-cylinder right-hand side. The car of a
 * 'vehicle' load has its own road speed: k ω in gear, coasting on its road load when declutched (integrated
 * after every step), and a gear change in free speed is an inelastic clutch engagement conserving the angular
 * momentum through the gear train (setOperatingPoint; clutchLoss ledger). After a
 * split in a closed cylinder of a fixed-speed multi-cylinder engine only that cylinder is re-evaluated
 * (refreshCylinder).
 *
 * ── Continuous state y (Float64Array, layout in cycle-common.ts) ─────────────────────────────
 *   θ (crank deg), ω (rad/s)                                   mechanics (dθ/dt = ω, dω/dt = α or 0)
 *   per cylinder: N[NS] (mol), U (J), m_bg (kg)                 open phase (single well-mixed zone)
 *   per cylinder: U_tot (J), S_u (J/K), m_u, m_b, m_e (kg)      closed phase (two zones at common p)
 *   intake plenum N, U, m_bg;  exhaust plenum N, U, m_bg        always
 *   per cylinder: K, k (J), swirl (kg m²/s)                     K–k turbulence (Poulos & Heywood 1983)
 *   external-EGR mass scalars m_egr (each cylinder, intake plenum) ⊂ m_bg
 *   ledgers: ∫p dV and ∫Q̇_wall (net: cylinder loss − intake-port gain) and valve mass flows per
 *   cylinder; venturi/outlet enthalpy and species flows, gross plenum inflows (warm-up relaxation)
 * m_bg is a passive burned-gas (products: residual + external EGR) mass scalar carried with every
 * stream; m_egr marks the external-EGR part of it (reset at EVO, when everything becomes residual).
 *
 * ── Open phase (EVO → IVC) ───────────────────────────────────────────────────────────────────
 * Fresh charge (humid air + fuel vapour at φ + EGR, thermo/fuels.freshCharge) at the intake
 * mixture temperature and ambient pressure → carburettor venturi (orifice, C_D·A_t scaled by the
 * butterfly throttleArea ratio; with manifolds.venturiDiameter the venturi and a butterfly throttle as two
 * compressible restrictions in series, gas-exchange CarburettorFlowModel) → intake plenum → intake valve
 * (ValveFlowModel effective area of the lash-corrected LiftProfile — polydyne or the cam's own lobe — ×
 * orificeFlow, bidirectional, direction-dependent C_D; side valves add the pocket → bore transfer) →
 * cylinder → exhaust valve → exhaust plenum → outlet orifice → ambient. Plenums and cylinder are
 * gas-exchange/Plenum volumes (conserved N, U; the stream carries the upstream stagnation
 * enthalpy and composition), so the network conserves mass, species and energy to round-off.
 * Heat: Woschni (gas-exchange constants) over head/valves/piston/liner (wallHeatLoss); intake-port
 * heat transfer into the intake plenum while gas flows through the port (Dittus–Boelter with the
 * port Reynolds number, × intakePortHeatTransferMultiplier; see Cylinder.evaluate).
 *
 * ── Closed phase (IVC → EVO) ─────────────────────────────────────────────────────────────────
 * Before the first burned gas: one zone, T from U_tot (frozen composition). With burned gas: two
 * zones at common p, conserved states U_tot (dU = −p dV − Q̇ + P_spark), S_u (dS_u = s_u dm_u −
 * Q̇_u/T_u), m_u, m_b; closure.ts recovers (p, T_u, T_b) by Newton (burned gas in chemical
 * equilibrium at (T_b, p)). After burn-out: burned zone only (UV equilibrium).
 * Combustion (default 'entrainment'): IgnitionSystem (coil circuit, breakdown, arc/glow, Herweg–
 * Maly kernel; trembler-magneto: a spark TRAIN per timer contact) → at hand-off the Keck/Tabaczynski
 * entrainment + burn-up ODEs with u′ from the K–k model (u_T = C_T u′; option: Keck 1982 empirical u_T,
 * ℓ_T), λ = Taylor microscale, S_L from the flame-speed tables (at the TRAPPED fuel blend and φ, see
 * trappedMixture) and A_f from the exact sphere ∩ chamber geometry (combustion CombustionChamber: disc or
 * L-head) of the entrained volume V_e = V_b + (m_e − m_b)/ρ_u, mapped onto the chamber in proportion
 * V_chamber/V (the lumped crevice volume holds charge distributed over the chamber).
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
 *
 * Sub-model factories: createChamber (combustion createChamber), createValveLift / createValveFlow
 * (createLiftProfile / createValveFlowModel), createKnockOscillator (chemistry createKnockOscillator),
 * createIgnitionSystem (inductive or one trembler-magneto system per cylinder); updateVenturi (the CFR
 * venturi-as-throttle or the series venturi + butterfly carburettor); the crank train, friction and load
 * model of the free-speed right-hand side.
 */
import { DEG, RAD2DEG, R_UNIVERSAL } from '../core/constants';
import type { EngineSpec, WallSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import type { CycleSummary, EngineCycleSummary } from '../core/snapshot';
import { FUEL_SPECIES, NE, NS, SP } from '../core/species';
import type { EquilibriumSolver } from '../equilibrium';
import { createKnockOscillator, type IgnitionDelayModel, type KnockOscillator, type LivengoodWuIntegrator, type ZeldovichKinetics } from '../chemistry';
import {
  createChamber,
  DEFAULT_LENGTH_SCALE_FRACTION,
  DEFAULT_TURBULENCE_PARAMS,
  type CombustionChamber,
  type FlameGeometry,
  type TurbulenceParams,
} from '../combustion';
import {
  CarburettorFlowModel,
  createLiftProfile,
  createValveFlowModel,
  gasStateFromTPX,
  hasSeparateThrottle,
  newGasState,
  newOrificeFlow,
  newSeriesOrificeFlow,
  orificeFlow,
  Plenum,
  resolveValveLash,
  throttleArea,
  type GasState,
  type LiftProfile,
  type OrificeFlow,
  type SeriesOrificeFlow,
  type ValveFlowModel,
} from '../gas-exchange';
import { IgnitionSystem, type IgnitionCommand } from '../ignition';
import {
  FrictionTorqueModel,
  LoadModel,
  MultiCylinderCrankTrain,
  pnhFmep,
  type CrankTrainDynamics,
  type KinematicState,
  type LoadModelKind,
  type SliderCrank,
} from '../mechanics';
import { completeCombustionProducts, freshCharge, fuelFromSelection, humidAir, lowerHeatingValue, type FuelBlend } from '../thermo/fuels';
import { mixMolarMass, temperatureFromUMolar } from '../thermo';
import { MOLAR_MASS } from '../thermo/thermo';
import type { ZoneClosure } from './closure';
import { Cylinder } from './cylinder';
import {
  blockLedgerIndex,
  cloneOp,
  cylinderStateIndex,
  elementMolesOf,
  FREE_MODE_MIN_OMEGA,
  I_EBG,
  I_EN,
  I_EU,
  I_GEH,
  I_GEN,
  I_GIH,
  I_GIK,
  I_GIN,
  I_GIV,
  I_HO,
  I_HV,
  I_IBG,
  I_IEG,
  I_IN,
  I_IU,
  I_LO,
  I_LV,
  I_OM,
  I_TH,
  INITIAL_EXHAUST_T,
  massOf,
  MODE_OPEN,
  now,
  OUTLET_DISCHARGE_COEFFICIENT,
  resolveDelayModel,
  stateLength,
  TIME_TARGET_REL_TOL,
  WARMUP_EXTRA_MAX,
  WARMUP_MASS_TOL,
  newCycleTrace,
  type ClosedCycleInit,
  type CycleProfile,
  type CycleTrace,
} from './cycle-common';
import { resolveCycleOptions, sanitizeOperatingPoint, type CycleModelOptions } from './options';

// historical exports of this module (the layout and helpers now live in cycle-common.ts)
export {
  I_TH, I_OM, I_CN, I_CU, I_CBG, I_UT, I_SU, I_MU, I_MB, I_ME, I_IN, I_IU, I_IBG, I_EN, I_EU, I_EBG, I_TK, I_TKE, I_SW,
  I_W, I_Q, I_MIVI, I_MIVO, I_MEVO, I_MEVI, I_HV, I_HO, I_LV, I_LO, I_MK, I_GIH, I_GIBG, I_GIN, I_GEH, I_GEBG, I_GEN,
  I_GIK, I_CEG, I_IEG, I_GIV, I_QS, I_CRU, I_CRB, NY,
  CYLINDER_BLOCK_SIZE, cylinderStateIndex, stateLength, blockLedgerIndex,
  MODE_OPEN, MODE_SINGLE, MODE_TWO, MODE_BURNED,
  OUTLET_DISCHARGE_COEFFICIENT, WIEBE_SEED, FREE_MODE_MIN_RPM, DP_WINDOW_DEG, DP_SUBDIV,
  newCycleTrace, expLinearIntegral, swapNO, resolveDelayModel,
  type CycleTrace, type ClosedCycleInit, type CycleProfile, type CylinderStateIndex,
} from './cycle-common';
export { VENTURI_DISCHARGE_COEFFICIENT } from './options';
export { Cylinder } from './cylinder';

/**
 * The engine cycle model (see file header). Construct, then drive with {@link stepUntil} /
 * {@link runCycles}; read the evaluated fields (p, Tu, … of cylinder 0; `cylinders[i]` for the others)
 * for output.
 */
export class CycleModel {
  readonly spec: EngineSpec;
  readonly opts: CycleModelOptions;
  /** Operating point in effect for the current cycle, and the latest requested one. */
  op: OperatingPoint;
  pendingOp: OperatingPoint;
  readonly closedInit: ClosedCycleInit | null;
  /** The cylinders (index = cylinder number − 1); cylinders[0] is the reference (legacy) cylinder. */
  readonly cylinders: Cylinder[];
  /** cylinders[0]. */
  readonly c0: Cylinder;

  // ---- sub-models (engine) ----
  /** Whole-engine friction (FMEP over N·V_d; per-cylinder piston terms summed by torqueCylinders). */
  readonly friction: FrictionTorqueModel;
  /** Rigid crank with every cylinder at its throw phase (free-speed dynamics, inertia). */
  readonly crank: MultiCylinderCrankTrain;
  /** Load on the crank in free-speed mode (op.load / loadTorque), configured from the operating point in effect. */
  readonly load: LoadModel;
  /** Scratch load model: setOperatingPoint validates a new load on it before committing. */
  private readonly loadCheck: LoadModel;
  /** Venturi + butterfly carburettor (spec.manifolds.venturiDiameter), or null (the CFR venturi-as-throttle). */
  readonly carb: CarburettorFlowModel | null;
  delayModel: IgnitionDelayModel;
  readonly intake: Plenum;
  readonly exhaust: Plenum;
  readonly fresh: GasState = newGasState();
  readonly ambient: GasState = newGasState();
  fuel: FuelBlend;
  /** K–k model constants and the dissipation-length fraction f_L (options). */
  readonly turbParams: TurbulenceParams;
  readonly lengthFraction: number;
  /** Keck-1982 turbulent-flame closure selected. */
  readonly keck: boolean;
  /** Kinetic NO carried in the burned-zone thermodynamics (options.burnedNOThermo; closure.noCoupled). */
  readonly noCoupled: boolean;

  // ---- time / state ----
  /** Length of y (NY + 42 per additional cylinder). */
  readonly ny: number;
  readonly y: Float64Array;
  t = 0;
  /** Completed engine cycles (θ wraps at +360°). */
  cycle = 0;
  /** True once a closed-cycle-only run reached its end angle. */
  finished = false;
  /** Set when a step wrapped the engine cycle (θ = 360 → −360). */
  wrapped = false;
  /** Set by the model when an output-worthy event happened (spark command, knock onset); stepUntil returns. */
  interrupt = false;

  // RK scratch
  private readonly dy: Float64Array; // derivative at (t, y) (valid when derivValid)
  private readonly k2: Float64Array;
  private readonly k3: Float64Array;
  private readonly k4: Float64Array;
  private readonly yt: Float64Array;
  /** State at the start of the current step (read by the cylinders' splits). */
  readonly y0: Float64Array;
  private readonly dy0: Float64Array;
  /** Kahan compensation of the RK state update (and its copy at the step start). */
  private readonly yComp: Float64Array;
  private readonly yComp0: Float64Array;
  /** True when dy is the derivative at (t, y). */
  derivValid = false;
  /** The only cylinder whose state changed since dy was last valid (null: none known / several / engine). */
  private dirtyCyl: Cylinder | null = null;

  // ---- per-cycle composition and boundary data ----
  readonly Xfresh = new Float64Array(NS);
  readonly airX = new Float64Array(NS);
  private cdaVenturi = 0;
  /** Fresh-charge mass per unit mass of humid air through the venturi, 1/(1 − Y_fuel − Y_EGR). */
  private freshPerAir = 1;
  private cdaOutlet = 0;
  private throttleFor = NaN;
  private pCrankcase = 101325;
  /** Per-cylinder thermodynamic pressures for the crank train (last evaluation), Pa. */
  private readonly pCyl: Float64Array;

  // ---- evaluated engine quantities (last RHS evaluation) ----
  mdotV = 0;
  mdotO = 0;
  pInt = 0;
  pExh = 0;
  alpha = 0;
  stiffness = 0;

  // ---- engine-cycle accumulators ----
  /** @internal Completed-cycle summaries (pushed by the cylinders). */
  readonly summaries: CycleSummary[] = [];
  private tEngineStart = 0;
  private ventEngineStart = 0;
  private omegaEngineStart = 0;
  /** ∫ T_load ω dt over the current engine cycle (free speed), J. */
  private loadWork = 0;
  /**
   * Kinetic-energy change of the coupled load's reflected inertia, ½J_L(ω_b² − ω_a²), summed over the pieces of
   * the current engine cycle between load (gear) changes, J — the piece in progress starts at loadPieceOmega.
   * Clutch-engagement jumps fall between the pieces (excluded).
   */
  private loadKEPieces = 0;
  private loadPieceOmega = 0;
  /** Engine-side kinetic-energy jumps ½J_e(ω′² − ω²) of the clutch engagements in the current engine cycle, J. */
  private engineKEJumps = 0;
  /** Energy dissipated by clutch engagements since the start of the run (free speed, 'vehicle' load), J. */
  clutchLoss = 0;
  private clutchLossEngineStart = 0;
  /** Road speed of the declutched car (load 'neutral'), m/s (in gear the road speed is k ω: vehicleSpeed()). */
  private coastSpeed = 0;
  /** Road distance of the car since the start of the run, m (engine summary: mean road speed). */
  private roadDistance = 0;
  private roadDistanceEngineStart = 0;

  // ---- trace ----
  traceCycle = -1;
  trace: CycleTrace | null = null;
  // ---- profiling ----
  readonly prof: CycleProfile = {
    rhsOpen: 0, rhsClosed: 0, closure: 0, flame: 0, heat: 0, ignition: 0, knock: 0, no: 0, steps: 0, rhsEvals: 0, eqSolves: 0,
  };

  // ---- scratch objects ----
  /** Venturi flow (a SeriesOrificeFlow, carrying the series solve's warm start, with a carburettor). */
  private readonly ofV: OrificeFlow;
  private readonly ofO = newOrificeFlow();
  private readonly Nscr = new Float64Array(NS);
  private readonly Xscr = new Float64Array(NS);

  constructor(spec: EngineSpec, op: OperatingPoint, options: Partial<CycleModelOptions> = {}, closedInit: ClosedCycleInit | null = null) {
    this.spec = spec;
    this.opts = resolveCycleOptions(options, spec);
    this.closedInit = closedInit;
    const n = cylinderCount(spec);
    if (closedInit && n > 1) throw new Error('CycleModel: closed-cycle-only runs (runClosedCycle) are single-cylinder');
    this.pendingOp = sanitizeOperatingPoint(spec, op);
    this.op = cloneOp(this.pendingOp);
    this.fuel = fuelFromSelection(this.op.fuel);
    this.delayModel = resolveDelayModel(this.opts.ignitionDelayModel);
    const nt = this.opts.burnedNOThermo;
    this.noCoupled = nt === 'kinetic' || (nt === 'auto' && this.opts.combustionModel === 'entrainment');
    this.keck = this.opts.turbulentFlameClosure === 'keck1982';
    this.turbParams = { ...DEFAULT_TURBULENCE_PARAMS, cBeta: this.opts.turbulenceProduction };
    this.lengthFraction = DEFAULT_LENGTH_SCALE_FRACTION * this.opts.turbulenceLengthScaleFactor;
    // (chambers with a block-deck surface — the L-head — append one heat ledger per cylinder)
    const blockLedgers = spec.geometry.chamber === 'l-head';
    const ny = stateLength(n) + (blockLedgers ? n : 0);
    this.ny = ny;
    this.y = new Float64Array(ny);
    this.dy = new Float64Array(ny);
    this.k2 = new Float64Array(ny);
    this.k3 = new Float64Array(ny);
    this.k4 = new Float64Array(ny);
    this.yt = new Float64Array(ny);
    this.y0 = new Float64Array(ny);
    this.dy0 = new Float64Array(ny);
    this.yComp = new Float64Array(ny);
    this.yComp0 = new Float64Array(ny);
    const cyls: Cylinder[] = [];
    // (identical cylinders share the chamber tables: createChamber caches them, each cylinder gets its own evaluator)
    for (let i = 0; i < n; i++) cyls.push(new Cylinder(this, i, cylinderStateIndex(i), blockLedgers ? blockLedgerIndex(n, i) : -1));
    const c0 = cyls[0];
    this.cylinders = cyls;
    this.c0 = c0;
    this.friction = new FrictionTorqueModel(c0.kin, 0.5, n);
    this.crank = MultiCylinderCrankTrain.fromSpec(spec, c0.kin);
    this.load = new LoadModel(spec.vehicle);
    this.loadCheck = new LoadModel(spec.vehicle);
    this.pCyl = new Float64Array(n);
    this.carb = hasSeparateThrottle(spec.manifolds)
      ? new CarburettorFlowModel(spec.manifolds, { venturiDischargeCoefficient: this.opts.venturiDischargeCoefficient, restrictionArea: this.opts.intakeRestrictionArea })
      : null;
    this.ofV = this.carb ? newSeriesOrificeFlow() : newOrificeFlow();
    const X0 = new Float64Array(NS);
    X0[SP.N2] = 1;
    this.intake = new Plenum(spec.manifolds.intakeVolume, 300, 1e5, X0);
    this.exhaust = new Plenum(spec.manifolds.exhaustVolume, 300, 1e5, X0);
    this.initState();
  }

  // ===========================================================================================
  // Sub-model factories (extension points)
  // ===========================================================================================

  /**
   * Lift profile of a cylinder's valve: the ValveSpec's cam (default polydyne; three-arc flat-follower lobe;
   * measured table) at the running clearance ValveSpec.lash, else options.valveLash (gas-exchange
   * createLiftProfile / resolveValveLash). Local crank angles.
   */
  createValveLift(c: Cylinder, kind: 'intake' | 'exhaust'): LiftProfile {
    const v = kind === 'intake' ? this.spec.intakeValve : this.spec.exhaustValve;
    void c;
    return createLiftProfile(v, resolveValveLash(v, this.opts.valveLash));
  }

  /**
   * Flow model (effective area, discharge coefficients, inflow jet) of a cylinder's valve (gas-exchange
   * createValveFlowModel: overhead valves as before; side valves with the pocket-roof masking stage and the
   * pocket → bore transfer section in series).
   */
  createValveFlow(c: Cylinder, kind: 'intake' | 'exhaust'): ValveFlowModel {
    const s = this.spec;
    void c;
    return createValveFlowModel(kind === 'intake' ? s.intakeValve : s.exhaustValve, kind, s);
  }

  /**
   * Chamber geometry of a cylinder (combustion createChamber: flat disc or L-head). The fast tables cover
   * the tallest bore column: stroke + h_TDC at the lowest compression ratio (+ 2 %); identical cylinders
   * share the immutable tables (createChamber's cache), each gets its own evaluator.
   */
  createChamber(c: Cylinder): CombustionChamber {
    const hMax = c.kin.pistonTravel + c.kin.clearanceHeightTDCForCR(this.spec.geometry.compressionRatioRange[0]);
    return createChamber(this.spec, { maxHeight: hMax * 1.02 });
  }

  /**
   * Acoustic knock modes of a cylinder (chemistry createKnockOscillator: Bessel modes of the bore for a
   * flat disc, the depth-averaged planform modes of an L-head — shared by the cylinders — with
   * options.knockSensor / mapoBand; null band = unfiltered).
   */
  createKnockOscillator(c: Cylinder): KnockOscillator {
    void c;
    return createKnockOscillator(this.spec, { decayTime: this.opts.knockDecayTime, sensor: this.opts.knockSensor, band: this.opts.mapoBand });
  }

  /**
   * Ignition system of a cylinder (built at its first cycle start; 'entrainment' combustion only). A
   * trembler-magneto ignition gets one system per cylinder — its own coil (coilOverrides[i]), gap and kernel
   * on its timer segment, the supply data shared (what ignition createIgnitionSystems builds per cylinder) —
   * fed from the operating point's ignition source.
   */
  createIgnitionSystem(c: Cylinder): IgnitionSystem {
    const io = this.opts.ignition;
    const s = this.spec;
    if (s.ignition.type === 'trembler-magneto') {
      return new IgnitionSystem(s.ignition, s.sparkPlug, {
        ...io,
        kernel: { handoffIntegralScaleMultiple: this.opts.kernelHandoffMultiple, ...(io.kernel ?? {}) },
        cylinder: c.index,
        firingOffsetDeg: c.offsetDeg,
        ignitionSource: this.op.ignitionSource ?? 'magneto',
      });
    }
    return new IgnitionSystem(s.ignition, s.sparkPlug, {
      makeSparkDiode: true,
      ...io,
      kernel: { handoffIntegralScaleMultiple: this.opts.kernelHandoffMultiple, ...(io.kernel ?? {}) },
    });
  }

  // ===========================================================================================
  // Public control
  // ===========================================================================================

  /**
   * Change operating conditions: rpm / throttle / load / speed mode / ignition source now, the rest at the next
   * cycle start. The load is validated on the candidate operating point BEFORE anything is committed: a load
   * the model cannot configure throws here and leaves the model unchanged (code review: an inherited gear name
   * used to poison pendingOp, and every later call threw).
   *
   * Gear changes of a 'vehicle' load (engaging, changing or leaving a gear, neutral included) act on the car's
   * road speed v, a state of its own (see {@link vehicleSpeed}):
   *  - leaving gear (to neutral or to a dynamometer load) keeps ω and the car's speed v = k_old ω — the car
   *    then coasts (neutral) or leaves the model (a constant / brake load);
   *  - free speed, engaging a gear while the car has a speed of its own (from a gear or from neutral): an
   *    INELASTIC clutch engagement through the rigid gear train. The impulse conserves the generalised
   *    momentum along the crank, J_e ω + J_L v/k = (J_e + J_L) ω′ (J_e = J_rot + Σ_i J_m(θ_i) at the current
   *    angle, J_L and k = roadSpeedPerOmega of the new gear, signed: reverse while rolling forward gives a
   *    negative v/k), and dissipates ½ J_e J_L/(J_e + J_L)(ω − v/k)² in the clutch (clutchLoss ledger). The
   *    stall guard holds ω′ ≥ FREE_MODE_MIN_RPM (as in the right-hand side); the ledger books the exact
   *    kinetic-energy change, so it then also holds the energy the guard adds. Round 1 kept ω and made the
   *    car's speed jump (a downshift destroyed ≈ 48 kJ of the car's kinetic energy, neutral → high created it);
   *  - engaging a gear from a dynamometer load (constant / brake: no car before) starts with the car rolling
   *    at the gear's speed k ω — a new test set-up, not a clutch event (no ledger entry);
   *  - fixed speed: the dynamometer is an ideal speed source; it holds ω through any gear change and the car
   *    follows at k ω in gear (the dynamometer supplies or absorbs the impulse: no clutch ledger, nothing in
   *    the engine summary, which in fixed speed has no kinetic-energy term).
   * The engine summary keeps the impulsive exchange out of the brake torque (see makeEngineSummary).
   */
  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    // undefined keys of the patch keep the pending value (a spread would overwrite it with undefined);
    // non-finite values are repaired from the pending point by the sanitiser (validation round 2)
    const next = { ...this.pendingOp };
    for (const k of Object.keys(patch) as (keyof OperatingPoint)[]) {
      const v = patch[k];
      if (v !== undefined) (next as Record<string, unknown>)[k] = v;
    }
    if (patch.fuel) next.fuel = { ...patch.fuel };
    const pending = sanitizeOperatingPoint(this.spec, next, this.pendingOp);
    // (validate on a scratch load model with exactly what `this.load.configure(o)` will see below; throws)
    this.loadCheck.configure({ loadTorque: pending.loadTorque, load: pending.load, ambientPressure: this.op.ambientPressure, ambientTemperature: this.op.ambientTemperature });
    this.pendingOp = pending;
    const o = this.op;
    o.rpm = this.pendingOp.rpm;
    o.throttle = this.pendingOp.throttle;
    o.loadTorque = this.pendingOp.loadTorque;
    o.speedMode = this.pendingOp.speedMode;
    if (this.pendingOp.load) o.load = { ...this.pendingOp.load };
    else if (o.load !== undefined) delete o.load;
    if (this.pendingOp.ignitionSource !== undefined) o.ignitionSource = this.pendingOp.ignitionSource;
    else if (o.ignitionSource !== undefined) delete o.ignitionSource;
    if (o.speedMode === 'fixed') this.y[I_OM] = (o.rpm * 2 * Math.PI) / 60;
    // the load model follows the new load at once; a gear change acts through the clutch (changeGear); the
    // ignition source switches at each coil's next timer make
    const ld = this.load;
    const prevKind = ld.kind;
    const prevGear = ld.gear;
    const prevInertia = ld.inertia();
    const vCar = this.vehicleSpeed();
    ld.configure(o);
    if (ld.kind !== prevKind || ld.gear !== prevGear) this.changeGear(prevKind, prevInertia, vCar);
    this.applyIgnitionSource();
    this.updateFriction();
    this.invalidateAll();
  }

  /**
   * Engine-side inertia about the crank at the current angle, J_e = J_rot + Σ_i J_m(θ_i), kg m² (crank +
   * flywheel and every piston/rod mechanism; not the load's reflected inertia).
   */
  engineInertia(): number {
    let J = this.c0.dyn.rotatingInertia;
    for (const c of this.cylinders) J += c.dyn.mechanismInertia(c.theta * DEG);
    return J;
  }

  /**
   * Road speed of the vehicle, m/s: k ω while a gear is engaged ('vehicle' load; negative in reverse), the
   * coasting speed when declutched (gear 'neutral': m_c dv/dt = −F_road(v), integrated with every step,
   * LoadModel.coastSpeedAfter; the car starts at rest at a cold start in neutral), 0 for loads without a car.
   */
  vehicleSpeed(): number {
    const ld = this.load;
    if (ld.kind === 'vehicle') return ld.vehicleSpeed(this.y[I_OM]);
    return ld.kind === 'neutral' ? this.coastSpeed : 0;
  }

  /**
   * The load's engagement changed (setOperatingPoint; `this.load` already configured for the new load): car
   * speed, clutch engagement and the engine-summary bookkeeping (see setOperatingPoint for the physics).
   * @param prevKind load kind before the change
   * @param prevInertia reflected load inertia before the change, kg m²
   * @param vCar road speed of the car before the change, m/s (0 without a car)
   */
  private changeGear(prevKind: LoadModelKind, prevInertia: number, vCar: number): void {
    const ld = this.load;
    const y = this.y;
    const w = y[I_OM];
    // the piece of the coupled load's kinetic energy that ends here (engine summary: loadInertiaTorque)
    this.loadKEPieces += 0.5 * prevInertia * (w * w - this.loadPieceOmega * this.loadPieceOmega);
    this.loadPieceOmega = w;
    if (ld.kind === 'neutral') {
      // declutched: the car keeps its speed and coasts (a car new to the model starts at rest)
      this.coastSpeed = prevKind === 'vehicle' || prevKind === 'neutral' ? vCar : 0;
      return;
    }
    if (ld.kind !== 'vehicle') return; // a dynamometer load: the car leaves the model
    const hadCar = prevKind === 'vehicle' || prevKind === 'neutral';
    if (!hadCar || this.op.speedMode !== 'free') return; // new set-up / dynamometer-held speed: v = k ω
    // inelastic clutch engagement: J_e ω + J_L v/k conserved
    const k = ld.roadSpeedPerOmega;
    const JL = ld.inertia();
    const Je = this.engineInertia();
    const wCar = vCar / k;
    let w1 = (Je * w + JL * wCar) / (Je + JL);
    if (w1 < FREE_MODE_MIN_OMEGA) w1 = FREE_MODE_MIN_OMEGA; // stall guard (numerical, as in evaluate)
    const keBefore = 0.5 * Je * w * w + 0.5 * JL * wCar * wCar;
    const keAfter = 0.5 * (Je + JL) * w1 * w1;
    this.clutchLoss += keBefore - keAfter;
    this.engineKEJumps += 0.5 * Je * (w1 * w1 - w * w);
    this.loadPieceOmega = w1;
    y[I_OM] = w1;
    this.yComp[I_OM] = 0;
  }

  /** Trembler-magneto: request the operating point's supply (MAG / BAT) on every cylinder's ignition. */
  private applyIgnitionSource(): void {
    const kind = this.op.ignitionSource ?? 'magneto';
    for (const c of this.cylinders) if (c.ign && c.ign.supply) c.ign.ignitionSource = kind;
  }

  /** Restart from t = 0 with the latest operating point (runs the warm-up cycles again). */
  reset(): void {
    this.op = cloneOp(this.pendingOp);
    this.initState();
  }

  /** Completed-cycle summaries since the last drain (oldest first; every cylinder's). */
  drainSummaries(): CycleSummary[] {
    return this.summaries.splice(0, this.summaries.length);
  }

  /** Current engine crank angle (cylinder 0's), deg. */
  get theta(): number {
    return this.y[I_TH];
  }

  /** Current speed, rev/min. */
  get rpm(): number {
    return (this.y[I_OM] * 60) / (2 * Math.PI);
  }

  /**
   * Run n complete engine cycles (to n more engine-cycle boundaries) without producing snapshots and
   * return the summaries completed meanwhile (one per cylinder and cycle; removed from the drain queue).
   */
  runCycles(n: number): CycleSummary[] {
    const before = this.summaries.length;
    const c0 = this.cycle;
    while (this.cycle - c0 < n && !this.finished) this.stepUntil(Infinity, 360);
    return this.summaries.splice(before, this.summaries.length - before);
  }

  /** Record a full-resolution trace of cylinder 0 over engine cycle `cycle` (default: the next complete cycle). */
  recordTrace(cycle?: number): void {
    this.traceCycle = cycle ?? (this.y[I_TH] <= -360 + 1e-9 ? this.cycle : this.cycle + 1);
    this.trace = newCycleTrace();
    if (this.traceCycle === this.cycle) {
      if (!this.derivValid) this.refresh();
      this.c0.recordSample(this.c0.p + this.c0.knockOscillation());
    }
  }

  /**
   * Advance with steps until t ≥ tTarget or θ reaches thetaTarget (≤ 360, current engine cycle), or the
   * engine cycle wraps, or a closed-cycle-only run finishes.
   */
  stepUntil(tTarget: number, thetaTarget: number): void {
    this.wrapped = false;
    this.interrupt = false;
    let guard = 0;
    while (!this.finished && guard++ < 10_000_000) {
      if (this.t >= tTarget - TIME_TARGET_REL_TOL * Math.max(1, Math.abs(tTarget))) return;
      // (θ = 360 is never "reached" without the engine wrap: a step that was limited short of a landing —
      // fine/stiffness steps that are not grid-aligned, e.g. a cylinder's fine phase across the wrap —
      // may end within 1e-9° below it, and the next step then lands on it)
      if (thetaTarget < 360 && this.y[I_TH] >= thetaTarget - 1e-9) return;
      this.step(tTarget, thetaTarget);
      if (this.wrapped || this.interrupt) return;
    }
  }

  // ===========================================================================================
  // Initialisation
  // ===========================================================================================

  private initState(): void {
    const y = this.y;
    const cyls = this.cylinders;
    y.fill(0);
    this.yComp.fill(0);
    this.t = 0;
    this.cycle = 0;
    this.finished = false;
    this.summaries.length = 0;
    for (const c of cyls) c.resetForInit();
    this.stiffness = 0;
    this.alpha = 0;
    this.pInt = 0;
    // the car (a 'vehicle' load): at rest when declutched, k ω in gear; no clutch losses yet
    this.coastSpeed = 0;
    this.roadDistance = 0;
    this.clutchLoss = 0;
    this.invalidateAll();
    this.fuel = fuelFromSelection(this.op.fuel);
    y[I_OM] = (this.op.rpm * 2 * Math.PI) / 60;
    this.prepareCycleData();
    for (const c of cyls) c.prepareCycleData();
    const op = this.op;
    if (this.closedInit) {
      this.c0.initClosedOnly(this.closedInit);
      return;
    }
    // Cold start at θ = −360 (cylinder 0 at gas-exchange TDC; cylinder i at its local angle, its first
    // local cycle partial): cylinders full of complete-combustion products at ambient pressure, plenums
    // at ambient pressure (fresh charge / products).
    y[I_TH] = -360;
    // (a motored engine never contains products: start it from fresh charge everywhere)
    const motored = this.opts.combustionModel === 'none';
    const prod = motored ? this.Xfresh : completeCombustionProducts(this.Xfresh);
    for (const c of cyls) c.initGas(prod, this.Xfresh, motored);
    this.intake.setTPX(op.intakeMixtureTemperature, op.ambientPressure, this.Xfresh);
    this.exhaust.setTPX(INITIAL_EXHAUST_T, op.ambientPressure, prod);
    for (let k = 0; k < NS; k++) {
      y[I_IN + k] = this.intake.N[k];
      y[I_EN + k] = this.exhaust.N[k];
    }
    y[I_IU] = this.intake.U;
    y[I_EU] = this.exhaust.U;
    // the fresh charge in the intake plenum carries its external EGR (products) from the start
    y[I_IBG] = this.intake.mass() * this.op.egrFraction;
    y[I_IEG] = y[I_IBG];
    y[I_EBG] = motored ? 0 : this.exhaust.mass();
    for (const c of cyls) {
      c.mode = MODE_OPEN;
      c.resetCycleFlags();
      c.startCycle();
    }
    this.engineCycleStart();
    this.refresh();
    // warm-up: silent cycles; after each, both plenums (the slow modes of the network) are
    // extrapolated to their periodic state (relaxPlenum).
    // (the warm-up is extended — up to WARMUP_EXTRA_MAX cycles — while the intake plenum's mass still
    // changes by more than WARMUP_MASS_TOL per cycle: throttled operation, validation round 2)
    let target = this.opts.warmupCycles;
    const L0 = new Float64Array(this.ny);
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
    const tEnd = this.t;
    this.t = 0;
    this.cycle = 0;
    this.summaries.length = 0;
    // (cylinder 0 is at its cycle start; cylinders i ≥ 1 are inside the local cycle that began during the
    // warm-up — it becomes their cycle 0 — and keep their model times relative to the new t = 0)
    this.c0.tCycleStart = 0;
    this.c0.cycle = 0;
    for (let i = 1; i < cyls.length; i++) {
      cyls[i].shiftTimes(-tEnd);
      cyls[i].cycle = 0;
    }
    this.engineCycleStart();
    this.invalidateAll();
    this.refresh();
  }

  /** Plenums of a closed-cycle-only run: irrelevant, kept at ambient. */
  initClosedOnlyPlenums(): void {
    const y = this.y;
    this.intake.setTPX(this.op.intakeMixtureTemperature, this.op.ambientPressure, this.Xfresh);
    this.exhaust.setTPX(INITIAL_EXHAUST_T, this.op.ambientPressure, this.Xfresh);
    for (let k = 0; k < NS; k++) {
      y[I_IN + k] = this.intake.N[k];
      y[I_EN + k] = this.exhaust.N[k];
    }
    y[I_IU] = this.intake.U;
    y[I_EU] = this.exhaust.U;
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
    this.invalidateAll();
    this.refresh();
  }

  /**
   * Warm-up acceleration of the intake plenum's MASS (pressure) mode, after relaxPlenum: the plenum
   * mass obeys dm/dt = ṁ_venturi(p) − ṁ_engine(p) with p = mRT/V, a relaxation of rate
   *   k = (RT/V)(|∂ṁ_venturi/∂p| + ṁ_engine/p),
   * the engine flow ∝ p (trapped mass ∝ intake density). Over one warm-up cycle of duration T_c the
   * mode decays by λ = e^{−k T_c}, so the periodic mass is m* = m₁ + (m₁ − m₀) λ/(1 − λ) (m₀, m₁ at
   * the start and end of the cycle; |∂ṁ_v/∂p| and ṁ_engine are cycle means from the ledgers; ṁ_engine
   * summed over the cylinders). Wide open (the CFR): λ ≈ 0, nothing changes; throttled the mode is slow
   * (λ ≈ 0.9 at throttle 0.1), and round 1's first emitted cycle was 26–30 % off in IMEP at throttle ≤ 0.1
   * (validation round 2).
   */
  /** Relative change of the intake-plenum mass over the last warm-up cycle (before extrapolation). */
  private intakeMassChange = 0;
  /** True during the silent warm-up cycles of initState. */
  warmingUp = false;

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
    const cyls = this.cylinders;
    const ix0 = cyls[0].ix;
    let mFlow = y[ix0.MIVI] - L0[ix0.MIVI] - (y[ix0.MIVO] - L0[ix0.MIVO]);
    for (let i = 1; i < cyls.length; i++) {
      const ix = cyls[i].ix;
      mFlow += y[ix.MIVI] - L0[ix.MIVI] - (y[ix.MIVO] - L0[ix.MIVO]);
    }
    const mEng = mFlow / Tc;
    pl.updateState();
    const st = pl.state;
    if (!(m0 > 0 && m1 > 0 && st.p > 0 && dmdp >= 0 && mEng >= 0)) return;
    const k = ((st.R * st.T) / pl.volume) * (dmdp + mEng / st.p);
    const lam = Math.exp(-k * Tc);
    const f = Math.min(lam / (1 - lam), 50);
    const mStar = m1 + (m1 - m0) * f;
    if (!(mStar > 0.2 * m1 && mStar < 5 * m1) || Math.abs(mStar / m1 - 1) < 1e-9) return;
    const X = this.Xscr;
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
    this.invalidateAll();
    this.refresh();
  }

  /** Engine-level data that depend on the operating point (engine-cycle start; then each cylinder's part). */
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
    this.pCrankcase = op.ambientPressure + this.opts.crankcaseGaugePressure;
    this.delayModel = resolveDelayModel(this.opts.ignitionDelayModel);
    // load (air density from the ambient state) and ignition supply of the operating point in effect
    this.load.configure(op);
    this.applyIgnitionSource();
    this.updateFriction();
  }

  /**
   * Effective flow area C_D·A of the carburettor (cached per throttle opening). A spec with a separate
   * butterfly (manifolds.venturiDiameter + throttle) sets the plate opening of its series venturi +
   * butterfly model (CarburettorFlowModel, solved per evaluation); otherwise throttleDiameter is the
   * venturi scaled by the butterfly open-area ratio (the CFR form).
   */
  private updateVenturi(): void {
    const op = this.op;
    if (op.throttle === this.throttleFor) return;
    if (this.carb) {
      this.carb.setOpening(op.throttle);
      this.throttleFor = op.throttle;
      return;
    }
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

  // ===========================================================================================
  // Engine wrap and engine summary
  // ===========================================================================================

  /** θ = 360 (cylinder 0's wrap): close the engine cycle, apply the pending operating point, θ → −360. */
  engineWrap(): void {
    const y = this.y;
    const c0 = this.c0;
    const cyls = this.cylinders;
    const s = c0.endCycle();
    if (s) s.engine = this.makeEngineSummary();
    // new cycle
    const rpmNow = this.op.rpm;
    const mode = this.op.speedMode;
    this.op = cloneOp(this.pendingOp);
    if (mode === 'free' && this.op.speedMode === 'free') this.op.rpm = rpmNow; // keep ω (state)
    c0.applyCompressionRatio();
    this.prepareCycleData();
    y[I_TH] = -360;
    for (let i = 1; i < cyls.length; i++) cyls[i].angleShift += 720;
    if (this.op.speedMode === 'fixed') y[I_OM] = (this.op.rpm * 2 * Math.PI) / 60;
    this.cycle++;
    c0.beginCycle();
    this.engineCycleStart();
    this.wrapped = true;
    this.invalidateAll();
    this.refresh();
  }

  /** Engine-cycle start bookkeeping (engine summary accumulators). */
  private engineCycleStart(): void {
    this.tEngineStart = this.t;
    this.ventEngineStart = this.ventMass();
    this.omegaEngineStart = this.y[I_OM];
    this.loadWork = 0;
    this.loadKEPieces = 0;
    this.loadPieceOmega = this.y[I_OM];
    this.engineKEJumps = 0;
    this.clutchLossEngineStart = this.clutchLoss;
    this.roadDistanceEngineStart = this.roadDistance;
    for (const c of this.cylinders) c.wEngineStart = this.y[c.ix.W];
  }

  /**
   * Engine-level results of the engine cycle just completed (θ −360 → 360): mean speed, indicated torque
   * (∫p dV of all cylinders / 4π), friction torque (cycle mean of the FrictionTorqueModel at the mean
   * speed), brake torque = indicated − friction − ΔE_kin/4π — the ENGINE's output at the crankshaft (free
   * speed: ΔE_kin is the kinetic-energy change of the engine's own rotating and reciprocating parts,
   * ½(J_rot + ΣJ_m)ω², over the cycle; θ wraps to the same phase, so ΣJ_m is the same at both ends; the
   * impulsive exchange of a clutch engagement, engineKEJumps, is excluded), mean effective pressures over the
   * total displacement, carburettor air and fuel flow from the venturi ledger, η_v, BSFC (0 when the brake
   * power is not positive) and brake efficiency (LHV of the selected fuel), and the load: in free mode the
   * constant loadTorque, else the cycle mean ∫T_L ω dt/4π of the load model; the dynamometer's absorbed
   * (= brake) torque in fixed mode. Free speed with a vehicle (EngineSpec.vehicle): the mean road speed (road
   * distance / cycle time, coasting included), the kinetic-energy change of the car while coupled
   * (loadInertiaTorque, ½J_L Δω² per piece between gear changes / 4π) — the energy balance closes as
   * brake ≈ load + loadInertia — and the clutch losses of the cycle.
   * Round 1 put the car's J_L into ΔE_kin: brake torque then collapsed to the road-load torque (≈ 6 N m while
   * the engine delivered ≈ 90 N m accelerating in low gear; BSFC ≈ 10 kg/kWh) and a gear change booked the
   * car's kinetic-energy jump as brake work (brake 3844 N m, η_b 7.6: code review).
   */
  private makeEngineSummary(): EngineCycleSummary {
    const y = this.y;
    const cyls = this.cylinders;
    const n = cyls.length;
    const dt = this.t - this.tEngineStart;
    const om = dt > 0 ? (4 * Math.PI) / dt : y[I_OM];
    let W = 0;
    for (const c of cyls) W += y[c.ix.W] - c.wEngineStart;
    const VdTot = n * this.c0.kin.displacedVolume;
    const indicated = W / (4 * Math.PI);
    const friction = this.meanFrictionTorque(om);
    let dKE = 0;
    let dKELoad = 0;
    const free = this.op.speedMode === 'free';
    if (free) {
      const J = this.engineInertia();
      const w1 = y[I_OM];
      const w0 = this.omegaEngineStart;
      dKE = 0.5 * J * (w1 * w1 - w0 * w0) - this.engineKEJumps;
      const wp = this.loadPieceOmega;
      dKELoad = this.loadKEPieces + 0.5 * this.load.inertia() * (w1 * w1 - wp * wp);
    }
    const brake = indicated - friction - dKE / (4 * Math.PI);
    const mV = this.ventMass() - this.ventEngineStart;
    let mt = 0;
    let mf = 0;
    for (let k = 0; k < NS; k++) mt += this.Xfresh[k] * MOLAR_MASS[k];
    for (const sName of FUEL_SPECIES) mf += this.Xfresh[SP[sName]] * MOLAR_MASS[SP[sName]];
    const yFuel = mt > 0 ? mf / mt : 0;
    const fuelFlow = dt > 0 ? (mV * yFuel) / dt : 0;
    const airFlow = dt > 0 ? (mV * (1 - yFuel - this.op.egrFraction)) / dt : 0;
    const power = brake * om;
    const lhv = lowerHeatingValue(this.fuel);
    const ld = this.load;
    const loadTorque = !free ? brake : ld.kind === 'constant' ? this.op.loadTorque : this.loadWork / (4 * Math.PI);
    const summary: EngineCycleSummary = {
      rpmMean: (om * 60) / (2 * Math.PI),
      indicatedTorque: indicated,
      frictionTorque: friction,
      brakeTorque: brake,
      brakePower: power,
      bmep: (4 * Math.PI * brake) / VdTot,
      imepNet: W / VdTot,
      fmep: (4 * Math.PI * friction) / VdTot,
      airMassFlow: airFlow,
      fuelMassFlow: fuelFlow,
      volumetricEfficiency: this.volumetricEfficiencyOf(mV),
      // (undefined without positive brake power: reported as 0 — summaries stay finite)
      bsfc: power > 0 && fuelFlow > 0 ? fuelFlow / power : 0,
      brakeEfficiency: fuelFlow > 0 && lhv > 0 ? power / (fuelFlow * lhv) : 0,
      loadTorque,
    };
    // mean road speed of the car: distance / time (in gear ∫k ω dt = k Δθ; coasting in neutral included)
    if (free && (ld.kind === 'vehicle' || ld.kind === 'neutral')) {
      summary.vehicleSpeed = dt > 0 ? (this.roadDistance - this.roadDistanceEngineStart) / dt : this.vehicleSpeed();
    }
    if (free && this.spec.vehicle) {
      summary.loadInertiaTorque = dKELoad / (4 * Math.PI);
      summary.clutchLoss = this.clutchLoss - this.clutchLossEngineStart;
    }
    return summary;
  }

  /**
   * Cycle-mean friction torque of the engine at constant ω (positive = loss): the mean of
   * FrictionTorqueModel.torqueCylinders over a revolution, (T_c + N F_c⟨|x′|⟩)·ω/√(ω² + ε²) + N c_v⟨x′²⟩ω
   * (T_c the whole engine's constant part, F_c and c_v per cylinder) — exact for the model's distribution
   * (⟨|x′|⟩ = travel/π, ⟨x′²⟩ tabulated).
   */
  meanFrictionTorque(om: number): number {
    const f = this.friction;
    const e = f.smoothingOmega;
    const s = om / Math.sqrt(om * om + e * e);
    const n = f.cylinders;
    // (one cylinder: the single-cylinder expression, unchanged)
    if (n === 1) return (f.constantTorque + f.coulombForce * f.meanAbsDxdTheta) * s + f.viscousCoefficient * f.meanSqDxdTheta * om;
    return (f.constantTorque + n * f.coulombForce * f.meanAbsDxdTheta) * s + n * f.viscousCoefficient * f.meanSqDxdTheta * om;
  }

  // ===========================================================================================
  // Stepping
  // ===========================================================================================

  /** Evaluate the RHS at the current state (end-of-step fields and the cached derivative). */
  refresh(): void {
    this.evaluate(this.t, this.y, this.dy);
    this.derivValid = true;
    this.dirtyCyl = null;
  }

  /** The cached derivative is invalid because of an engine-level change (every cylinder affected). */
  invalidateAll(): void {
    this.derivValid = false;
    this.dirtyCyl = null;
  }

  /** The cached derivative is invalid because cylinder c's state changed. */
  invalidate(c: Cylinder): void {
    if (this.derivValid) {
      this.derivValid = false;
      this.dirtyCyl = c;
    } else if (this.dirtyCyl !== c) this.dirtyCyl = null;
  }

  /**
   * Re-evaluate after a split / event in cylinder c. A closed cylinder of a fixed-speed multi-cylinder
   * engine is decoupled from the rest (no valve flow, no crank coupling), so when c alone changed only c
   * is re-evaluated (its block of dy recomputed in place); otherwise — and always for one cylinder — the
   * whole engine (refresh).
   */
  refreshCylinder(c: Cylinder): void {
    if (
      this.cylinders.length > 1 &&
      this.op.speedMode === 'fixed' &&
      c.mode !== MODE_OPEN &&
      c.evalMode !== MODE_OPEN &&
      (this.derivValid || this.dirtyCyl === c)
    ) {
      c.clearRates(this.dy);
      c.evaluate(this.t, this.y, this.dy, true);
      if (this.opts.profile) this.prof.rhsEvals++;
      this.derivValid = true;
      this.dirtyCyl = null;
      return;
    }
    this.refresh();
  }

  private step(tTarget: number, thetaTarget: number): void {
    const y = this.y;
    if (!this.derivValid) this.refresh();
    const P = this.opts.profile;
    const th = y[I_TH];
    const cyls = this.cylinders;
    const n = cyls.length;
    for (let i = 0; i < n; i++) cyls[i].cacheNOStartState();
    // ---- step limit (minimum over the cylinders) ----
    // (the angle limits are converted with the same kinematics as the landing time below — in 'free'
    // mode a limit maxStepDeg/ω stopped a hair short of the next grid angle while decelerating and a
    // sliver step followed: 36 % of all steps, validation round 2)
    let hMax = this.timeToAngle(this.opts.maxStepDeg);
    let fine = false;
    for (let i = 0; i < n && !fine; i++) fine = cyls[i].finePhase();
    if (fine) hMax = Math.min(hMax, this.timeToAngle(this.opts.fineStepDeg));
    for (let i = 0; i < n; i++) hMax = cyls[i].limitStep(hMax);
    if (this.stiffness > 0) hMax = Math.min(hMax, 2.0 / this.stiffness);
    // ---- crank-angle landing: grid, events (every cylinder's, in engine angles), target ----
    const grid = this.opts.maxStepDeg;
    let thLand = -360 + (Math.floor((th + 360) / grid + 1e-7) + 1) * grid;
    for (let i = 0; i < n; i++) {
      const a = cyls[i].nextEventEngineAngle();
      if (a < thLand + 1e-9) thLand = a;
    }
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
      // Livengood–Wu crossing inside the step (earliest over the cylinders) → shorten to the crossing;
      // the onset is forced at the end of the shortened step for every cylinder crossing there
      let fMin = 1;
      let hit = false;
      for (let i = 0; i < n; i++) {
        const c = cyls[i];
        c.lwFraction = -1;
        if (!c.lwArmed()) continue;
        c.tauStepEnd = c.tauNow();
        const f = c.lwCrossingFraction(h, c.tauStepEnd);
        if (f >= 0 && f < 1 - 1e-6 && h * f > 1e-10 && attempt < 5) {
          c.lwFraction = f;
          if (f < fMin) fMin = f;
          hit = true;
        }
      }
      if (hit) {
        h *= fMin;
        landing = false;
        this.restore(tStart);
        for (let i = 0; i < n; i++) {
          const f = cyls[i].lwFraction;
          cyls[i].lwForce = f >= 0 && f <= fMin * (1 + 1e-9);
        }
        continue;
      }
      break;
    }
    this.t = tStart + h;
    if (P) this.prof.steps++;
    // load energy ∫T_L ω dt of a speed-dependent load over the engine cycle (trapezoid; engine summary)
    if (this.op.speedMode === 'free' && this.load.kind !== 'constant') {
      const w0 = this.y0[I_OM];
      const w1 = y[I_OM];
      this.loadWork += 0.5 * h * (this.load.torque(w0) * w0 + this.load.torque(w1) * w1);
    }
    // the car's road distance (in gear k Δθ) and, declutched, its coasting speed (independent of the engine)
    if (this.load.kind === 'vehicle') {
      this.roadDistance += this.load.roadSpeedPerOmega * (y[I_TH] - this.y0[I_TH]) * DEG;
    } else if (this.load.kind === 'neutral') {
      const v0 = this.coastSpeed;
      const v1 = this.load.coastSpeedAfter(v0, h);
      this.roadDistance += 0.5 * h * (v0 + v1);
      this.coastSpeed = v1;
    }
    // ---- operator splits and bookkeeping ----
    for (let i = 0; i < n; i++) cyls[i].afterStep(h);
    // ---- events at the landing angle ----
    if (landing) {
      const c0 = this.c0;
      const mode0 = c0.mode;
      const mb0 = y[c0.ix.MB];
      for (let i = 0; i < n; i++) {
        const c = cyls[i];
        c.handleEvents(c.localAngle(y[I_TH]));
      }
      // record the post-event state too when an event changed it (IVC/EVO/TDC burn)
      if (this.trace && this.cycle === this.traceCycle && (c0.mode !== mode0 || y[c0.ix.MB] !== mb0)) {
        if (!this.derivValid) this.refresh();
        c0.recordSample(c0.p + c0.knockOscillation());
      }
      // the pressure interpolation of the next step starts from the post-event state
      if (!this.derivValid) this.refresh();
      for (let i = 0; i < n; i++) cyls[i].setInterpolationStart();
    }
  }

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
    const ny = this.ny;
    const hh = 0.5 * h;
    for (let i = 0; i < ny; i++) yt[i] = y[i] + hh * k1[i];
    this.evaluate(t + hh, yt, k2);
    for (let i = 0; i < ny; i++) yt[i] = y[i] + hh * k2[i];
    this.evaluate(t + hh, yt, k3);
    for (let i = 0; i < ny; i++) yt[i] = y[i] + h * k3[i];
    this.evaluate(t + h, yt, k4);
    const h6 = h / 6;
    // compensated (Kahan) summation of the step increments: the conserved states and ledgers exchange
    // many small increments per cycle (mass between the zones and the crevice, the plenums and the
    // ledgers), whose rounding otherwise accumulates to ≈ 1e-12 of the network mass over 12 cycles
    const cmp = this.yComp;
    for (let i = 0; i < ny; i++) {
      const inc = h6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) - cmp[i];
      const t2 = y[i] + inc;
      cmp[i] = t2 - y[i] - inc;
      y[i] = t2;
    }
    // end-of-step evaluation (becomes k1 of the next step unless a split changes the state)
    this.evaluate(t + h, y, this.dy);
    this.derivValid = true;
    this.dirtyCyl = null;
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

  // ===========================================================================================
  // Right-hand side
  // ===========================================================================================

  private evaluate(t: number, y: Float64Array, dy: Float64Array): void {
    const P = this.opts.profile;
    const tp = P ? now() : 0;
    dy.fill(0);
    const th = y[I_TH];
    const om = y[I_OM];
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
    // (with a separate butterfly: venturi and throttle plate as compressible restrictions in series)
    const amb = this.ambient;
    let mV = this.carb
      ? this.carb.flow(amb.p, amb.T, amb.R, amb.gamma, ip.state.p, ip.state.T, ip.state.R, ip.state.gamma, this.ofV as SeriesOrificeFlow)
      : orificeFlow(this.cdaVenturi, amb.p, amb.T, amb.R, amb.gamma, ip.state.p, ip.state.T, ip.state.R, ip.state.gamma, this.ofV);
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
    // ---- cylinders: each adds its valve flows (and intake-port heat) to the plenum rates begun here ----
    ip.beginRates(0, 0);
    ep.beginRates(0, 0);
    const cyls = this.cylinders;
    const n = cyls.length;
    let open = false;
    let sCyl = 0;
    let sIv = 0;
    let sEv = 0;
    for (let i = 0; i < n; i++) {
      const c = cyls[i];
      c.evaluate(t, y, dy, i > 0);
      if (c.evalMode === MODE_OPEN) {
        if (!open) {
          open = true;
          sCyl = c.stiffCyl;
          sIv = c.aIv;
          sEv = c.aEv;
        } else {
          if (c.stiffCyl > sCyl) sCyl = c.stiffCyl;
          sIv += c.aIv;
          sEv += c.aEv;
        }
      }
    }
    // stiffness bound (Gershgorin) of the linearised pressure dynamics: cylinder rows (max) and the
    // plenum rows with every open cylinder's valve conductance (0 while all valves are shut)
    if (open) {
      const si = ip.state;
      const se = ep.state;
      const ki = (si.gamma * si.R * si.T) / ip.volume;
      const ke = (se.gamma * se.R * se.T) / ep.volume;
      const aV = Math.abs(this.ofV.dmdotdpa) + Math.abs(this.ofV.dmdotdpb);
      const aO = Math.abs(this.ofO.dmdotdpa) + Math.abs(this.ofO.dmdotdpb);
      this.stiffness = Math.max(sCyl, ki * (aV + sIv), ke * (sEv + aO));
    } else this.stiffness = 0;
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
    // ---- mechanics ----
    dy[I_TH] = om * RAD2DEG;
    if (this.op.speedMode === 'free') {
      // rigid crank: Σ_i gas + gravity − ½J′_m ω² torques of the cylinders at their throw phases, the
      // whole-engine friction and −T_load, over J_rot + J_load + ΣJ_m (mechanics MultiCylinderCrankTrain,
      // FrictionTorqueModel.torqueCylinders, LoadModel); the THERMODYNAMIC cylinder pressures (the
      // synthesised knock field exerts no net force on the piston)
      const thr = th * DEG;
      const dx = this.crank.update(thr);
      const pc = this.pCyl;
      for (let i = 0; i < n; i++) pc[i] = cyls[i].p;
      const ld = this.load;
      const ext = -ld.torque(om) + this.friction.torqueCylinders(dx, om);
      let a = this.crank.angularAcceleration(thr, om, pc, this.pCrankcase, ext, ld.inertia());
      // Stall guard (numerical): a time-stepped cycle cannot represent a stopped crank, so the speed
      // is held at FREE_MODE_MIN_RPM when the torques would decelerate it further.
      if (om <= FREE_MODE_MIN_OMEGA && a < 0) a = 0;
      this.alpha = a;
      dy[I_OM] = a;
    } else {
      this.alpha = 0;
    }
    if (P) {
      this.prof.rhsEvals++;
      const dt = now() - tp;
      if (this.c0.mode === MODE_OPEN) this.prof.rhsOpen += dt;
      else this.prof.rhsClosed += dt;
    }
  }

  // ===========================================================================================
  // Ledgers and inventories (conservation checks)
  // ===========================================================================================

  /** Net mass through the carburettor venturi since the start of the run (ledger I_LV), kg. */
  ventMass(): number {
    let m = 0;
    for (let k = 0; k < NS; k++) m += this.y[I_LV + k] * MOLAR_MASS[k];
    return m;
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

  /**
   * Volumetric efficiency of a net venturi mass mVent inducted over one cycle (Heywood's definition, as
   * measured on the CFR): DRY air through the carburettor / (ambient dry-air partial density × total
   * displacement), the dry air = mVent × the dry-air mass fraction of the fresh stream.
   */
  volumetricEfficiencyOf(mVent: number): number {
    const Yf = this.freshMassFractions();
    const xw = this.airX[SP.H2O];
    // dry-air partial density of the ambient (humid) air, kg/m³ (0 if the "air" is saturated steam)
    const rhoDry = xw < 1 ? ((mixMolarMass(this.airX) - xw * MOLAR_MASS[SP.H2O]) * this.op.ambientPressure) / (R_UNIVERSAL * this.op.ambientTemperature) : 0;
    const n = this.cylinders.length;
    const Vd = n === 1 ? this.c0.kin.displacedVolume : n * this.c0.kin.displacedVolume;
    return rhoDry > 0 ? (mVent * Yf.dryAir) / (rhoDry * Vd) : 0;
  }

  /** Σ over the cylinders of ∫p dV since the start of the run (ledgers W), J. */
  workLedger(): number {
    let w = 0;
    for (const c of this.cylinders) w += this.y[c.ix.W];
    return w;
  }

  /** Σ over the cylinders of the net wall-heat ledgers Q (cylinder walls − intake-port gain), J. */
  heatLedger(): number {
    let q = 0;
    for (const c of this.cylinders) q += this.y[c.ix.Q];
    return q;
  }

  /** Moles of each species in the cylinders + plenums (closed phase: unburned + burned equilibrium + crevice), mol. */
  speciesInventory(out: Float64Array = new Float64Array(NS)): Float64Array {
    const y = this.y;
    for (let k = 0; k < NS; k++) out[k] = y[I_IN + k] + y[I_EN + k];
    for (const c of this.cylinders) c.addSpecies(out);
    return out;
  }

  /** Element moles in the cylinders + plenums (exact in every mode), mol. */
  elementInventory(out: Float64Array = new Float64Array(NE)): Float64Array {
    const y = this.y;
    const N = this.Nscr;
    for (let k = 0; k < NS; k++) N[k] = y[I_IN + k] + y[I_EN + k];
    for (const c of this.cylinders) c.addOpenMoles(N);
    const b = elementMolesOf(N, out);
    for (const c of this.cylinders) c.addClosedElements(b);
    return b;
  }

  /** Internal energy of the cylinders + plenums, J. */
  systemEnergy(): number {
    const y = this.y;
    let E = y[I_IU] + y[I_EU];
    for (const c of this.cylinders) E += c.energy();
    return E;
  }

  /** Mass in the cylinders + plenums, kg. */
  systemMass(): number {
    const y = this.y;
    let m = 0;
    for (let k = 0; k < NS; k++) m += (y[I_IN + k] + y[I_EN + k]) * MOLAR_MASS[k];
    for (const c of this.cylinders) m = c.addMass(m);
    return m;
  }

  /** Cumulative net spark energy added to the gas, all cylinders, J. */
  get sparkEnergy(): number {
    const cs = this.cylinders;
    let s = cs[0].sparkEnergy;
    for (let i = 1; i < cs.length; i++) s += cs[i].sparkEnergy;
    return s;
  }

  /** Cumulative kernel electrode loss, all cylinders, J. */
  get electrodeLoss(): number {
    const cs = this.cylinders;
    let s = cs[0].electrodeLoss;
    for (let i = 1; i < cs.length; i++) s += cs[i].electrodeLoss;
    return s;
  }

  /** Failed closure solves since construction, all cylinders — diagnostic. */
  get closureFailures(): number {
    let s = 0;
    for (const c of this.cylinders) s += c.closureFailures;
    return s;
  }

  /** Intake-port heat into the intake plenum at the last evaluation, all cylinders, W. */
  get Qport(): number {
    const cs = this.cylinders;
    let s = cs[0].Qport;
    for (let i = 1; i < cs.length; i++) s += cs[i].Qport;
    return s;
  }

  // ===========================================================================================
  // Cylinder 0 (the single cylinder of the legacy API): sub-models, status and evaluated fields
  // ===========================================================================================

  get kin(): SliderCrank { return this.c0.kin; }
  get ks(): KinematicState { return this.c0.ks; }
  get dyn(): CrankTrainDynamics { return this.c0.dyn; }
  get ivLift(): LiftProfile { return this.c0.ivLift; }
  get evLift(): LiftProfile { return this.c0.evLift; }
  get ivFlow(): ValveFlowModel { return this.c0.ivFlow; }
  get evFlow(): ValveFlowModel { return this.c0.evFlow; }
  get flameGeom(): FlameGeometry { return this.c0.flameGeom; }
  get closure(): ZoneClosure { return this.c0.closure; }
  get eqAux(): EquilibriumSolver { return this.c0.eqAux; }
  get ign(): IgnitionSystem | null { return this.c0.ign; }
  get knockOsc(): KnockOscillator { return this.c0.knockOsc; }
  set knockOsc(k: KnockOscillator) { this.c0.knockOsc = k; }
  get lw(): LivengoodWuIntegrator { return this.c0.lw; }
  get zeld(): ZeldovichKinetics { return this.c0.zeld; }
  get cyl(): Plenum { return this.c0.cyl; }
  get mode(): number { return this.c0.mode; }
  get ignCmd(): IgnitionCommand { return this.c0.ignCmd; }
  get ivcDeg(): number { return this.c0.ivcDeg; }
  get evoDeg(): number { return this.c0.evoDeg; }
  get yRes(): number { return this.c0.yRes; }
  get yEgr(): number { return this.c0.yEgr; }
  get xDil(): number { return this.c0.xDil; }
  get xResMole(): number { return this.c0.xResMole; }
  get rhoInlet(): number { return this.c0.rhoInlet; }
  get keckInletSpeed(): number { return this.c0.keckInletSpeed; }
  get lwStage(): number { return this.c0.lwStage; }
  get walls(): WallSpec { return this.c0.walls; }
  get creviceOn(): boolean { return this.c0.creviceOn; }
  get Vcr(): number { return this.c0.Vcr; }
  get Tcr(): number { return this.c0.Tcr; }
  get mdotCr(): number { return this.c0.mdotCr; }
  get Qcr(): number { return this.c0.Qcr; }
  get wallExcess(): Float64Array { return this.c0.wallExcess; }
  get mIvc(): number { return this.c0.mIvc; }
  get pIvc(): number { return this.c0.pIvc; }
  get TIvc(): number { return this.c0.TIvc; }
  get VIvc(): number { return this.c0.VIvc; }
  get fuelMassIvc(): number { return this.c0.fuelMassIvc; }
  get nNO(): number { return this.c0.nNO; }
  get sparkFired(): boolean { return this.c0.sparkFired; }
  get kernelMassPrev(): number { return this.c0.kernelMassPrev; }
  get handedOff(): boolean { return this.c0.handedOff; }
  get flameActive(): boolean { return this.c0.flameActive; }
  get burnDone(): boolean { return this.c0.burnDone; }
  get misfire(): boolean { return this.c0.misfire; }
  get kernelQuenched(): boolean { return this.c0.kernelQuenched; }
  get knockOnset(): boolean { return this.c0.knockOnset; }
  get knockBurning(): boolean { return this.c0.knockBurning; }
  get tauAB(): number { return this.c0.tauAB; }
  get lwJ(): number { return this.c0.lwJ; }
  get knockAvailable(): boolean { return this.c0.knockAvailable; }
  get fuelTrapped(): FuelBlend { return this.c0.fuelTrapped; }
  get phiTrapped(): number { return this.c0.phiTrapped; }
  get knockOnsetDeg(): number { return this.c0.knockOnsetDeg; }
  get knockEndGasFraction(): number { return this.c0.knockEndGasFraction; }
  get mapo(): number { return this.c0.mapo; }
  get knockQdot(): number { return this.c0.knockQdot; }
  get tKnockOnset(): number { return this.c0.tKnockOnset; }
  get dwellSeen(): boolean { return this.c0.dwellSeen; }
  get tSparkCmd(): number { return this.c0.tSparkCmd; }
  get V(): number { return this.c0.V; }
  get Vdot(): number { return this.c0.Vdot; }
  get h(): number { return this.c0.h; }
  get p(): number { return this.c0.p; }
  get T(): number { return this.c0.T; }
  get Tu(): number { return this.c0.Tu; }
  get Tb(): number { return this.c0.Tb; }
  get mCyl(): number { return this.c0.mCyl; }
  get Qwall(): number { return this.c0.Qwall; }
  get Qu(): number { return this.c0.Qu; }
  get Qb(): number { return this.c0.Qb; }
  get mdotIv(): number { return this.c0.mdotIv; }
  get mdotEv(): number { return this.c0.mdotEv; }
  get SL(): number { return this.c0.SL; }
  get uPrime(): number { return this.c0.uPrime; }
  get L(): number { return this.c0.L; }
  get lambda(): number { return this.c0.lambda; }
  get rf(): number { return this.c0.rf; }
  get Af(): number { return this.c0.Af; }
  get rb(): number { return this.c0.rb; }
  get discScale(): number { return this.c0.discScale; }
  get mdotB(): number { return this.c0.mdotB; }
  get mdotE(): number { return this.c0.mdotE; }
  get burnSpeed(): number { return this.c0.burnSpeed; }
  get tauB(): number { return this.c0.tauB; }
  get rhoU(): number { return this.c0.rhoU; }
  get nuU(): number { return this.c0.nuU; }
  get dpdt(): number { return this.c0.dpdt; }
  get pMot(): number { return this.c0.pMot; }
  get hcoef(): number { return this.c0.hcoef; }
  get burnRateStep(): number { return this.c0.burnRateStep; }
  get hrr(): number { return this.c0.hrr; }
  get peakP(): number { return this.c0.peakP; }
  get peakPDeg(): number { return this.c0.peakPDeg; }
  get ca10(): number { return this.c0.ca10; }
  get ca50(): number { return this.c0.ca50; }
  get ca90(): number { return this.c0.ca90; }
  /** Burned-mass fraction of cylinder 0 (Cylinder.xb). */
  get xb(): number { return this.c0.xb; }
  /** Knock-integral progress of cylinder 0 (Cylinder.knockProgress). */
  get knockProgress(): number { return this.c0.knockProgress; }
  /** True once cylinder 0's entrained front has swept the chamber. */
  get frontAtWalls(): boolean { return this.c0.frontAtWalls; }
  /** Cylinder 0's flame radius for output. */
  flameRadius(): number { return this.c0.flameRadius(); }
  /** Cylinder 0's synthesised knock oscillation at the pickup, Pa. */
  knockOscillation(): number { return this.c0.knockOscillation(); }
  /** Cylinder 0's burned-zone kinetic NO mole fraction. */
  burnedNOFraction(): number { return this.c0.burnedNOFraction(); }
  /** Cylinder 0's crevice-zone internal energy, J. */
  creviceEnergy(): number { return this.c0.creviceEnergy(); }
}

// =============================================================================================
// Helpers
// =============================================================================================

/** Number of cylinders of a spec, with its layout validated (firing offsets in [0, 720), cylinder 1 at 0). */
export function cylinderCount(spec: EngineSpec): number {
  const n = spec.cylinders;
  if (!(Number.isInteger(n) && n >= 1)) throw new RangeError(`${spec.name}: cylinders must be a positive integer`);
  const l = spec.layout;
  if (!l) {
    if (n > 1) throw new RangeError(`${spec.name}: a ${n}-cylinder spec needs layout.firingOffsetDeg`);
    return 1;
  }
  const off = l.firingOffsetDeg;
  if (off.length !== n) throw new RangeError(`${spec.name}: layout.firingOffsetDeg has ${off.length} entries for ${n} cylinders`);
  if (off[0] !== 0) throw new RangeError(`${spec.name}: layout.firingOffsetDeg[0] must be 0 (θ is cylinder 1's angle)`);
  for (const o of off) if (!(o >= 0 && o < 720)) throw new RangeError(`${spec.name}: firing offsets must lie in [0, 720) deg`);
  return n;
}
