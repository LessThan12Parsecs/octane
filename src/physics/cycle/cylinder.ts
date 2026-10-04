/**
 * One cylinder of the cycle model: its state block, sub-models (slider crank, valves, flame geometry,
 * zone closure, ignition, knock, NO), the right-hand side of its gas (evalOpen / evalClosed, port heat
 * transfer, turbulence), its operator splits, events, per-cycle accumulators and summary. The engine
 * (cycle-model.ts CycleModel) owns time, crank angle and speed, the shared intake/exhaust plenums, the
 * venturi and outlet, the RK4 stepping and the warm-up; see the header of cycle-model.ts for the
 * physics. Every formula here is the single-cylinder model of validation rounds 1–2, moved verbatim:
 * cylinder 0 of an engine is bit-identical to the former single-cylinder CycleModel.
 *
 * Crank angles: a cylinder works in its LOCAL angle θ_i = θ + angleShift (deg, firing TDC = 0, cycle
 * [−360, 360)), θ the engine angle (cylinder 0's). angleShift = −firingOffsetDeg (mod 720) and changes by
 * ±720 at the engine wrap (θ: 360 → −360) and at the cylinder's own wrap (θ_i: 360 → −360); it is 0
 * for cylinder 0. Every *Deg field of a cylinder (events, ignition command, CA10…, knock onset, peak
 * angle) is local.
 *
 * Extension points (Model T integration):
 *  - chamber: every FlameGeometry / flatChamberAreas / discScale / integralLengthScale call site is a
 *    method of this class — evalOpen and evalClosed (areas, length scale, discScale, flame radius and
 *    front area, two-zone heat split, crevice burned fraction), clampRadiusGuess, endGasCircleRadius,
 *    handoffEntrainedMass, flameRadius, frontAtWalls, createBurnedZone (Lref) — plus
 *    engine-simulator.ts cylinderSnapshot (kernel/flame radius → area);
 *  - ignition: CycleModel.createIgnitionSystem(cylinder) builds this.ign; prepareIgnitionCommand sets
 *    the per-cycle command (inductive: dwell start and switch-off); ignitionSplit is the per-step call;
 *    finePhase decides the fine steps; EngineSimulator.nextDenseTime the dense snapshot sampling;
 *  - valve lift / flow: CycleModel.createValveLift / createValveFlow (factories);
 *  - knock: CycleModel.createKnockOscillator (factory).
 */
import { DEG, RAD2DEG, R_UNIVERSAL } from '../core/constants';
import type { WallSpec } from '../core/engine-spec';
import type { CycleSummary } from '../core/snapshot';
import { EL, FUEL_SPECIES, NE, NS, SP } from '../core/species';
import { EquilibriumSolver } from '../equilibrium';
import {
  autoignitionBurnTime,
  endGasBurnRate,
  KnockOscillator,
  LivengoodWuIntegrator,
  residualMoleFraction,
  TabulatedIgnitionDelay,
  twoZoneSoundSpeed,
  ZeldovichKinetics,
} from '../chemistry';
import {
  angularMomentumLengthScale,
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
  type TurbulenceState,
} from '../combustion';
import { newOrificeFlow, orificeFlow, Plenum, ValveFlowModel, ValveLiftProfile, newValveJet } from '../gas-exchange';
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
import { CrankTrainDynamics, newKinematicState, SliderCrank } from '../mechanics';
import { completeCombustionProducts, freshCharge, lowerHeatingValue, type FuelBlend } from '../thermo/fuels';
import { mixCpMass, mixCvMass, mixCvMolar, mixHMolar, mixMolarMass, mixSMass, mixSMolar, mixThermalConductivity, mixUMolar, mixViscosity, temperatureFromUMolar } from '../thermo';
import { MOLAR_MASS } from '../thermo/thermo';
import { ZoneClosure } from './closure';
import {
  CREVICE_RELAX_TIME,
  crossing,
  DP_SUBDIV,
  DP_WINDOW_DEG,
  elementMolesOf,
  EV_BDC_END,
  EV_BDC_START,
  EV_END,
  EV_EVO,
  EV_IGN,
  EV_IVC,
  EV_TDC,
  EV_VALVE,
  EV_WIEBE,
  EV_WRAP,
  expLinearIntegral,
  hermiteAt,
  hermiteMax,
  HSCR,
  I_EBG,
  I_GEBG,
  I_GEH,
  I_GEN,
  I_GIBG,
  I_GIH,
  I_GIK,
  I_GIN,
  I_IBG,
  I_IEG,
  I_OM,
  I_TH,
  INITIAL_CYLINDER_T,
  massOf,
  MAX_EVENTS,
  MODE_BURNED,
  MODE_OPEN,
  MODE_SINGLE,
  MODE_TWO,
  NO_FREEZE_T,
  now,
  swapNO,
  WIEBE_SEED,
  type CylinderStateIndex,
  type ClosedCycleInit,
} from './cycle-common';
import type { CycleModel } from './cycle-model';

export class Cylinder {
  /** The engine this cylinder belongs to. */
  readonly e: CycleModel;
  /** 0-based cylinder index (cylinder number − 1). */
  readonly index: number;
  /** Firing-TDC offset after cylinder 0's, deg (spec.layout.firingOffsetDeg). */
  readonly offsetDeg: number;
  /** θ_i = θ + angleShift, deg (see file header). */
  angleShift = 0;
  /** Absolute indices of this cylinder's state block. */
  readonly ix: Readonly<CylinderStateIndex>;
  private readonly iCN: number;
  private readonly iCU: number;
  private readonly iCBG: number;
  private readonly iUT: number;
  private readonly iSU: number;
  private readonly iMU: number;
  private readonly iMB: number;
  private readonly iME: number;
  private readonly iTK: number;
  private readonly iTKE: number;
  private readonly iSW: number;
  private readonly iW: number;
  private readonly iQ: number;
  private readonly iMIVI: number;
  private readonly iMIVO: number;
  private readonly iMEVO: number;
  private readonly iMEVI: number;
  private readonly iMK: number;
  private readonly iCEG: number;
  private readonly iQS: number;
  private readonly iCRU: number;
  private readonly iCRB: number;

  // ---- sub-models ----
  readonly kin: SliderCrank;
  readonly ks = newKinematicState();
  readonly dyn: CrankTrainDynamics;
  readonly ivLift: ValveLiftProfile;
  readonly evLift: ValveLiftProfile;
  readonly ivFlow: ValveFlowModel;
  readonly evFlow: ValveFlowModel;
  readonly flameGeom: FlameGeometry;
  readonly closure = new ZoneClosure();
  readonly eqAux = new EquilibriumSolver();
  ign: IgnitionSystem | null = null;
  knockOsc: KnockOscillator;
  readonly lw = new LivengoodWuIntegrator();
  readonly zeld = new ZeldovichKinetics();
  /** Open-phase gas zone (a Plenum of volume V(θ_i)). */
  readonly cyl: Plenum;

  // ---- cycle ----
  mode = MODE_OPEN;
  /**
   * Local cycle number: completed local cycles since t = 0 (the cycle in progress at t = 0 is cycle 0 for
   * every cylinder; for cylinders i ≥ 1 it began before t = 0). CycleSummary.cycle of this cylinder.
   */
  cycle = 0;
  /** True while the current local cycle began mid-cycle at a cold start (its summary is not emitted). */
  partialCycle = false;
  /** Mode at this cylinder's last right-hand-side evaluation (partial-refresh eligibility). */
  evalMode = MODE_OPEN;

  // ---- events of the current local cycle ----
  private readonly evAngle = new Float64Array(MAX_EVENTS);
  private readonly evKind = new Int32Array(MAX_EVENTS);
  private nEv = 0;
  private evNext = 0;
  readonly ignCmd: IgnitionCommand = { dwellStartDeg: -360, sparkDeg: 0 };
  ivcDeg = -152;
  evoDeg = 141;

  // ---- per-cycle composition data ----
  /** Residual mass fraction of the trapped charge (burned-gas scalar) and products dilution for S_L. */
  yRes = 0;
  /** External-EGR mass fraction of the trapped charge (part of yRes). */
  yEgr = 0;
  xDil = 0;
  xResMole = 0;
  /** Intake-port heat transfer: mult·(k/D)·πDL (W/K per unit Nu), 4/(πDμ), Pr^0.4, wall T. */
  private portHA = 0;
  private portReCoef = 0;
  private portPr04 = 1;
  private portTw = 0;
  /** Intake-port heat into the intake plenum at the last evaluation, W. */
  Qport = 0;
  /** Keck-1982 closure: inlet density ρ_i (intake plenum at IVC), mean inlet speed ū_i. */
  rhoInlet = 0;
  keckInletSpeed = 0;
  /** Two-stage knock integral: selected (and supported by the delay model); current stage (1, 2). */
  private twoStage = false;
  lwStage = 1;
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
  /** Surface heat ledgers and time at the start of the current local cycle. */
  private readonly qsCycleStart = new Float64Array(5);
  tCycleStart = 0;
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
  /** True once this cycle's dwell has started (the ignition state then belongs to this cycle). */
  dwellSeen = false;
  /** Time of this cycle's spark command (primary switch-off), s (NaN before). */
  tSparkCmd = NaN;
  /** Livengood–Wu crossing located by the engine's step redo: force the onset at the end of this step. */
  lwForce = false;
  /** Delay τ at the end of the current RK step (set by the engine's crossing test), s. */
  tauStepEnd = NaN;
  /** Fraction of the tentative step at which this cylinder's LW integral reaches 1 (−1: none; engine scratch). */
  lwFraction = -1;
  /** Burned-zone state at the start of the current step for the NO split (Strang): valid, T_b, p, total moles, X. */
  private noStartValid = false;
  private noTb0 = 0;
  private noP0 = 0;
  private noNb0 = 0;
  private readonly noX0 = new Float64Array(NS);
  /** Burned mass merged (with its NO) by mergeToBurned during the current step's splits, kg. */
  private noMergedMass = 0;

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
  pMot = 0;
  hcoef = 0;
  /** Open phase: Gershgorin terms of the valve flows (engine stiffness bound): kc·(a_iv + a_ev), a_iv, a_ev. */
  stiffCyl = 0;
  aIv = 0;
  aEv = 0;
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
  /** Venturi ledger (net mass, kg) at the start of the current local cycle. */
  private mVentCycleStart = 0;
  private noPpm = 0;
  private coFrac = 0;
  private xbEvo = 0;
  private closedHappened = false;
  /** ∫p dV of this cylinder at the start of the current ENGINE cycle (engine summary), J. */
  wEngineStart = 0;

  // ---- scratch objects ----
  private readonly ofIv = newOrificeFlow();
  private readonly ofEv = newOrificeFlow();
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

  /**
   * @param engine the engine (its spec, options and operating point must be set)
   * @param index 0-based cylinder index; @param ix its state indices
   * @param flameGeom flame geometry (shared by identical cylinders), or null to build one
   */
  constructor(engine: CycleModel, index: number, ix: Readonly<CylinderStateIndex>, flameGeom: FlameGeometry | null) {
    const spec = engine.spec;
    const opts = engine.opts;
    this.e = engine;
    this.index = index;
    this.offsetDeg = index === 0 ? 0 : (spec.layout?.firingOffsetDeg[index] ?? 0);
    this.ix = ix;
    this.iCN = ix.CN;
    this.iCU = ix.CU;
    this.iCBG = ix.CBG;
    this.iUT = ix.UT;
    this.iSU = ix.SU;
    this.iMU = ix.MU;
    this.iMB = ix.MB;
    this.iME = ix.ME;
    this.iTK = ix.TK;
    this.iTKE = ix.TKE;
    this.iSW = ix.SW;
    this.iW = ix.W;
    this.iQ = ix.Q;
    this.iMIVI = ix.MIVI;
    this.iMIVO = ix.MIVO;
    this.iMEVO = ix.MEVO;
    this.iMEVI = ix.MEVI;
    this.iMK = ix.MK;
    this.iCEG = ix.CEG;
    this.iQS = ix.QS;
    this.iCRU = ix.CRU;
    this.iCRB = ix.CRB;
    this.fuelTrapped = engine.fuel;
    this.phiTrapped = engine.op.equivalenceRatio;
    const g = spec.geometry;
    this.kin = new SliderCrank(g, engine.op.compressionRatio);
    this.dyn = new CrankTrainDynamics(this.kin, spec.masses);
    this.ivLift = engine.createValveLift(this, 'intake');
    this.evLift = engine.createValveLift(this, 'exhaust');
    this.ivFlow = engine.createValveFlow(this, 'intake');
    this.evFlow = engine.createValveFlow(this, 'exhaust');
    // Flame-geometry table must cover the tallest chamber: stroke + h_TDC at the lowest CR.
    if (flameGeom) this.flameGeom = flameGeom;
    else {
      const hMax = this.kin.pistonTravel + this.kin.clearanceHeightTDCForCR(g.compressionRatioRange[0]);
      this.flameGeom = new FlameGeometry(g.bore, spec.sparkPlug.gapCenter, { maxHeight: hMax * 1.02 });
    }
    this.knockOsc = engine.createKnockOscillator(this);
    const X0 = new Float64Array(NS);
    X0[SP.N2] = 1;
    this.cyl = new Plenum(this.kin.volume(-2 * Math.PI), 300, 1e5, X0);
    this.walls = { ...spec.walls };
    this.wallsLumped = !!spec.walls.thermalResistance && opts.wallTemperatureModel === 'lumped';
    this.Vcr = spec.geometry.quenchCreviceVolume ?? 0;
    // (the crevice zone follows the pressure continuously; the idealised prescribed burns — Wiebe, the
    // instantaneous fuel–air-cycle burn — keep the round-1 single-volume definition)
    const cmod = opts.combustionModel;
    this.creviceOn = opts.creviceModel && opts.heatTransfer && this.Vcr > 0 && (cmod === 'entrainment' || cmod === 'none');
  }

  // ===========================================================================================
  // Angles
  // ===========================================================================================

  /** Local crank angle of this cylinder at engine angle th, deg. */
  localAngle(th: number): number {
    return this.angleShift === 0 ? th : th + this.angleShift;
  }

  /** Current local crank angle, deg. */
  get theta(): number {
    return this.localAngle(this.e.y[I_TH]);
  }

  /** Engine angle of a local angle a of the current local cycle, deg. */
  engineAngle(a: number): number {
    return this.angleShift === 0 ? a : a - this.angleShift;
  }

  /** Engine angle of this cylinder's next event (Infinity if none left this local cycle), deg. */
  nextEventEngineAngle(): number {
    return this.evNext < this.nEv ? this.engineAngle(this.evAngle[this.evNext]) : Infinity;
  }

  /** Burned-mass fraction m_b/m of the closed charge; after EVO (same cycle) the value at EVO. */
  get xb(): number {
    if (this.mode === MODE_OPEN) return this.xbEvo;
    const y = this.e.y;
    const m = y[this.iMU] + y[this.iMB];
    return m > 0 ? y[this.iMB] / m : 0;
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

  // ===========================================================================================
  // Initialisation
  // ===========================================================================================

  /** Reset the per-cylinder solver state and guesses (CycleModel.initState, before prepareCycleData). */
  resetForInit(): void {
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
    this.T = 0;
    this.Tu = 0;
    this.Tb = 0;
    this.nuU = 1.6e-5;
    this.tauB = Infinity;
    this.mdotB = 0;
    this.lwForce = false;
    this.muTab.fill(0);
    this.kin.setCompressionRatio(this.e.op.compressionRatio);
    // cold start: cylinder 0 at the start of its cycle (θ = −360); cylinder i at its local angle there
    const off = this.offsetDeg;
    this.angleShift = off === 0 ? 0 : 720 - off;
    this.cycle = 0;
    this.partialCycle = off !== 0;
    {
      const w = this.e.spec.walls;
      const Tc = w.referenceCoolantTemperature ?? this.e.op.coolantTemperature;
      this.wallExcess[0] = w.headTemperature - Tc;
      this.wallExcess[1] = w.pistonTemperature - Tc;
      this.wallExcess[2] = w.linerTemperature - Tc;
      this.wallExcess[3] = w.intakeValveTemperature - Tc;
      this.wallExcess[4] = w.exhaustValveTemperature - Tc;
      this.qsCycleStart.fill(0);
      this.tCycleStart = 0;
    }
  }

  /**
   * Cold-start gas at ambient pressure in the cylinder's volume at its current local angle: residual
   * `prod` (complete-combustion products; the fresh charge when motored) at INITIAL_CYLINDER_T — the
   * single-cylinder start at θ = −360 — except for a cylinder starting mid-cycle BEFORE its firing TDC,
   * which holds `fresh` charge at the intake-mixture temperature (as if just inducted; 900 K residual
   * trapped there would autoignite during the compression of the partial first cycle).
   */
  initGas(prod: Float64Array, fresh: Float64Array, motored: boolean): void {
    const e = this.e;
    const y = e.y;
    const th = this.theta;
    const inducted = this.partialCycle && th < 0;
    this.cyl.volume = this.kin.volume(th * DEG);
    this.cyl.setTPX(inducted ? e.op.intakeMixtureTemperature : INITIAL_CYLINDER_T, e.op.ambientPressure, inducted ? fresh : prod);
    for (let k = 0; k < NS; k++) y[this.iCN + k] = this.cyl.N[k];
    y[this.iCU] = this.cyl.U;
    if (inducted) {
      // the fresh charge carries its external EGR (products) as marked burned gas
      y[this.iCBG] = this.cyl.mass() * e.op.egrFraction;
      y[this.iCEG] = y[this.iCBG];
    } else {
      y[this.iCBG] = motored ? 0 : this.cyl.mass();
      y[this.iCEG] = 0;
    }
  }

  /** Closed-cycle-only start from a prescribed state (single-cylinder validation runs). */
  initClosedOnly(ci: ClosedCycleInit): void {
    const e = this.e;
    const y = e.y;
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
    e.initClosedOnlyPlenums();
    this.cyl.volume = V;
    this.cyl.setTPX(ci.T, ci.p, X);
    for (let k = 0; k < NS; k++) y[this.iCN + k] = this.cyl.N[k];
    y[this.iCU] = this.cyl.U;
    y[this.iCBG] = this.cyl.mass() * (ci.residualMassFraction ?? 0);
    y[this.iCEG] = 0;
    const Sp = (2 * e.spec.geometry.stroke * e.op.rpm) / 60;
    const up = ci.uPrime ?? 0.5 * Sp;
    this.mode = MODE_OPEN;
    this.resetCycleFlags();
    this.startCycle();
    this.onIvc();
    y[this.iTKE] = 1.5 * this.mIvc * up * up;
    y[this.iTK] = 0;
    e.derivValid = false;
    e.refresh();
  }

  resetCycleFlags(): void {
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

  /**
   * Per-cylinder data that depend on the operating point (cycle start; after the engine's part of
   * CycleModel.prepareCycleData): wall temperatures, intake-port heat-transfer constants, viscosity table.
   */
  prepareCycleData(): void {
    const e = this.e;
    const op = e.op;
    const m = e.spec.manifolds;
    // wall temperatures: lumped model (coolant + R_i Q̄_i, updateWalls), else WallSpec shifted with the
    // coolant (UNVERIFIED: wall-to-coolant ΔT held at its value for the CFR's 100 °C boiling jacket).
    // (a spec without a reference coolant temperature: its wall temperatures are taken as the values at
    // the operating point's coolant temperature)
    const w = e.spec.walls;
    const Tref = w.referenceCoolantTemperature ?? op.coolantTemperature;
    const dT = op.coolantTemperature - Tref;
    let Tp = w.intakePortTemperature;
    if (this.wallsLumped) {
      const Tc = op.coolantTemperature;
      const ex = this.wallExcess;
      this.walls = {
        headTemperature: Tc + ex[0],
        pistonTemperature: Tc + ex[1],
        linerTemperature: Tc + ex[2],
        intakeValveTemperature: Tc + ex[3],
        exhaustValveTemperature: Tc + ex[4],
      };
      // intake port (in the water-cooled head): its excess over the coolant scales with the head's
      // (UNVERIFIED assumption; the port gets its heat through the head casting)
      if (Tp !== undefined) Tp = Tc + (Tp - Tref) * (ex[0] / Math.max(1, w.headTemperature - Tref)) - dT;
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
    if (e.opts.intakePortHeatTransferMultiplier > 0 && Tp !== undefined && Lp > 0 && Dp > 0) {
      this.portTw = Tp + dT;
      const Tf = 0.5 * (op.intakeMixtureTemperature + this.portTw); // film temperature
      const mu = mixViscosity(e.Xfresh, Tf);
      const k = mixThermalConductivity(e.Xfresh, Tf);
      const cp = mixCpMass(e.Xfresh, Tf);
      this.portReCoef = 4 / (Math.PI * Dp * mu);
      // h = mult · Nu k/D over the wetted area πDL, Nu = max(3.66, 0.023 Re^0.8 Pr^0.4)
      this.portHA = e.opts.intakePortHeatTransferMultiplier * (k / Dp) * Math.PI * Dp * Lp;
      this.portPr04 = Math.pow((mu * cp) / k, 0.4);
    } else {
      this.portHA = 0;
    }
    if (this.muTab[0] === 0) for (let i = 0; i < this.muTab.length; i++) this.muTab[i] = mixViscosity(e.Xfresh, this.muT0 + i * this.muDT);
  }

  /**
   * Ignition command of the local cycle starting now (inductive: dwell start = switch-off − dwell time
   * at the current speed, switch-off at −sparkAdvanceDeg). Extension point for the trembler-magneto timer
   * (command = the timer contact window of this cylinder).
   */
  prepareIgnitionCommand(): number {
    const op = this.e.op;
    const omegaDeg = this.e.y[I_OM] * RAD2DEG;
    const sparkDeg = -op.sparkAdvanceDeg;
    this.ignCmd.sparkDeg = sparkDeg;
    let dwell = sparkDeg - op.dwellTime * omegaDeg;
    if (dwell < -359.9) dwell = -359.9;
    this.ignCmd.dwellStartDeg = dwell;
    return dwell;
  }

  /** Local cycle start (θ_i = −360, or the current angle of a cold start): events and ignition command. */
  startCycle(): void {
    const e = this.e;
    const sparkDeg = -e.op.sparkAdvanceDeg;
    const dwell = this.prepareIgnitionCommand();
    this.ivcDeg = this.ivLift.seatCloseDeg;
    this.evoDeg = this.evLift.seatOpenDeg;
    const ci = e.closedInit;
    const cm = e.opts.combustionModel;
    this.nEv = 0;
    if (ci) {
      this.addEvent(ci.endDeg, EV_END);
    } else {
      this.addEvent(this.ivcDeg, EV_IVC);
      this.addEvent(this.evoDeg, EV_EVO);
      this.addEvent(this.ivLift.seatOpenDeg, EV_VALVE);
      this.addEvent(this.evLift.seatCloseDeg, EV_VALVE);
      this.addEvent(360, EV_WRAP);
    }
    this.addEvent(-180, EV_BDC_START);
    this.addEvent(180, EV_BDC_END);
    if (cm === 'instantaneous-at-tdc') this.addEvent(0, EV_TDC);
    if (cm === 'entrainment') {
      this.addEvent(dwell, EV_IGN);
      this.addEvent(sparkDeg, EV_IGN);
    }
    if (cm === 'wiebe') {
      this.wiebeStart = e.opts.wiebe.startDeg ?? sparkDeg;
      this.addEvent(this.wiebeStart, EV_WIEBE);
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
    const y = e.y;
    const th = this.theta;
    while (this.evNext < this.nEv && this.evAngle[this.evNext] <= th + 1e-9) {
      if (this.evKind[this.evNext] === EV_BDC_START) this.wBdcStart = y[this.iW];
      this.evNext++;
    }
    this.wCycleStart = y[this.iW];
    this.mVentCycleStart = e.ventMass();
    if (cm === 'entrainment' && !this.ign) this.ign = e.createIgnitionSystem(this);
    if (cm !== 'entrainment') this.ign = null;
  }

  private addEvent(a: number, k: number): void {
    this.evAngle[this.nEv] = a;
    this.evKind[this.nEv] = k;
    this.nEv++;
  }

  // ===========================================================================================
  // Step control
  // ===========================================================================================

  finePhase(): boolean {
    if (!this.ign) return false;
    const th = this.theta;
    if (th < this.ignCmd.sparkDeg - 1e-9 && !this.sparkFired) return false;
    if (this.misfire || this.kernelQuenched || this.burnDone || this.mode === MODE_OPEN) return false;
    if (!this.handedOff) return this.sparkFired || th >= this.ignCmd.sparkDeg - 1e-9;
    return this.xb < 0.02;
  }

  /** This cylinder's step limits (besides the fine phase and the engine's stiffness bound), s. */
  limitStep(hMax: number): number {
    const y = this.e.y;
    if (this.knockBurning) hMax = Math.min(hMax, 0.25 * this.tauAB);
    if (this.flameActive && this.tauB < Infinity) hMax = Math.min(hMax, Math.max(this.tauB, 1e-7));
    if (this.mode === MODE_TWO) {
      hMax = Math.min(hMax, 2 * this.e.opts.zoneHeatLossMinTime);
      if (this.mdotB > 0) hMax = Math.min(hMax, Math.max((0.2 * y[this.iMU]) / this.mdotB, 1e-9));
    }
    return hMax;
  }

  /** Cache the burned-zone equilibrium state at the step start (the last evaluation is at it). */
  cacheNOStartState(): void {
    const r = this.closure.eq.result;
    this.noStartValid = (this.mode === MODE_TWO || this.mode === MODE_BURNED) && this.e.y[this.iMB] > 0 && r.nTotal > 0 && this.Tb > NO_FREEZE_T;
    if (!this.noStartValid) return;
    this.noTb0 = this.Tb;
    this.noP0 = this.p;
    this.noNb0 = r.nTotal;
    this.noX0.set(r.X);
  }

  /** Start of the next step after events at the landing angle: the post-event pressure state. */
  setInterpolationStart(): void {
    this.hTheta = this.theta;
    this.hP = this.p;
    this.hDp = this.mode === MODE_OPEN ? 0 : this.dpdt / (this.e.y[I_OM] * RAD2DEG);
    this.hMode = this.mode;
  }

  // ===========================================================================================
  // Right-hand side
  // ===========================================================================================

  /**
   * This cylinder's right-hand side at the stage state (t, y): kinematics at its local angle, gas
   * (open or closed phase), its valve flows into the shared plenums (whose rates the engine has begun),
   * intake-port heat, turbulence and its ledgers. `acc`: accumulate into the shared-plenum entries of
   * dy (false for the first cylinder, which assigns them — the single-cylinder operation order).
   */
  evaluate(t: number, y: Float64Array, dy: Float64Array, acc: boolean): void {
    const e = this.e;
    const th = this.localAngle(y[I_TH]);
    const om = y[I_OM];
    const ks = this.kin.evaluate(th * DEG, this.ks);
    this.V = ks.volume;
    this.Vdot = ks.dVdTheta * om;
    this.h = ks.clearanceHeight;
    this.evalMode = this.mode;
    if (this.mode === MODE_OPEN) this.evalOpen(y, dy, th, acc);
    else this.evalClosed(t, y, dy, th);
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
      const ip = e.intake;
      const Re = this.portReCoef * Math.abs(this.mdotIv);
      let Nu = 0.023 * Math.pow(Re, 0.8) * this.portPr04;
      if (Nu < 3.66) Nu = 3.66;
      Qp = this.portHA * Nu * (this.portTw - ip.state.T);
      ip.dUdt += Qp;
      dy[I_GIK] = acc ? dy[I_GIK] + this.portHA * Nu : this.portHA * Nu;
    }
    this.Qport = Qp;
    // ---- turbulence ----
    const ts = this.turbState;
    ts.K = y[this.iTK];
    ts.k = y[this.iTKE];
    ts.swirl = y[this.iSW];
    turbulenceDerivatives(ts, this.turbIn, this.turbOut, e.turbParams);
    dy[this.iTK] = this.turbOut.dK;
    dy[this.iTKE] = this.turbOut.dk;
    dy[this.iSW] = this.turbOut.dSwirl;
    dy[this.iW] = this.p * this.Vdot;
    // net wall heat LOSS of the gas network: cylinder walls minus the intake-port gain
    dy[this.iQ] = this.Qwall - Qp;
  }

  /** Zero this cylinder's block of dy (before a partial re-evaluation). */
  clearRates(dy: Float64Array): void {
    for (let k = 0; k < NS; k++) dy[this.iCN + k] = 0;
    dy[this.iCU] = 0;
    dy[this.iCBG] = 0;
    dy[this.iUT] = 0;
    dy[this.iSU] = 0;
    dy[this.iMU] = 0;
    dy[this.iMB] = 0;
    dy[this.iME] = 0;
    dy[this.iTK] = 0;
    dy[this.iTKE] = 0;
    dy[this.iSW] = 0;
    dy[this.iW] = 0;
    dy[this.iQ] = 0;
    dy[this.iMIVI] = 0;
    dy[this.iMIVO] = 0;
    dy[this.iMEVO] = 0;
    dy[this.iMEVI] = 0;
    dy[this.iMK] = 0;
    dy[this.iCEG] = 0;
    for (let i = 0; i < 5; i++) dy[this.iQS + i] = 0;
    dy[this.iCRU] = 0;
    dy[this.iCRB] = 0;
  }

  private evalOpen(y: Float64Array, dy: Float64Array, th: number, acc: boolean): void {
    const e = this.e;
    const opts = e.opts;
    const P = opts.profile;
    const c = this.cyl;
    const ip = e.intake;
    const ep = e.exhaust;
    for (let k = 0; k < NS; k++) c.N[k] = y[this.iCN + k];
    c.U = y[this.iCU];
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
    const cdm = opts.dischargeCoefficientMultiplier;
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
    // stiffness bound (Gershgorin) of the linearised pressure dynamics: the cylinder row here, the
    // plenum rows (summed over the cylinders) in CycleModel.evaluate
    const kc = (sc.gamma * sc.R * sc.T) / this.V;
    const aIv = Math.abs(this.ofIv.dmdotdpa) + Math.abs(this.ofIv.dmdotdpb);
    const aEv = Math.abs(this.ofEv.dmdotdpa) + Math.abs(this.ofEv.dmdotdpb);
    this.aIv = aIv;
    this.aEv = aEv;
    this.stiffCyl = kc * (aIv + aEv);
    // heat transfer
    const tq = P ? now() : 0;
    let Q = 0;
    if (opts.heatTransfer) {
      const w = this.woschni;
      w.bore = e.spec.geometry.bore;
      w.pressure = sc.p;
      w.temperature = sc.T;
      w.meanPistonSpeed = this.meanPistonSpeed();
      w.phase = 'gas-exchange';
      w.variant = opts.woschniVariant;
      const hc = opts.woschniMultiplier * (opts.heatTransferCorrelation === 'hohenberg' ? hohenbergCoefficient(this.V, sc.p, sc.T, w.meanPistonSpeed) : woschniCoefficient(w));
      this.hcoef = hc;
      flatChamberAreas(e.spec.geometry.bore, this.h, e.spec.intakeValve, e.spec.exhaustValve, this.areas);
      Q = wallHeatLoss(hc, sc.T, this.areas, this.walls, this.heat);
      this.surfaceHeatRates(dy, 1);
    }
    if (P) e.prof.heat += now() - tq;
    this.Qwall = Q;
    this.Qu = Q;
    this.Qb = 0;
    // cylinder balances
    c.beginRates(-Q, this.Vdot);
    c.addFlow(mIv, si);
    c.addFlow(-mEv, se);
    for (let k = 0; k < NS; k++) dy[this.iCN + k] = c.dNdt[k];
    dy[this.iCU] = c.dUdt;
    const mI = massOf(ip.N);
    const mE = massOf(ep.N);
    const ybgC = m > 0 ? y[this.iCBG] / m : 0;
    const ybgI = mI > 0 ? y[I_IBG] / mI : 0;
    const ybgE = mE > 0 ? y[I_EBG] / mE : 0;
    const fIv = mIv > 0 ? mIv * ybgI : mIv * ybgC; // burned gas into the cylinder via the intake valve
    const fEv = mEv > 0 ? mEv * ybgC : mEv * ybgE; // burned gas out of the cylinder via the exhaust valve
    dy[this.iCBG] = fIv - fEv;
    dy[I_IBG] = acc ? dy[I_IBG] - fIv : -fIv;
    dy[I_EBG] = acc ? dy[I_EBG] + fEv : fEv;
    // external-EGR marker: in/out through the intake valve, out through the exhaust valve (exhaust
    // backflow is residual, unmarked)
    const yegC = m > 0 ? y[this.iCEG] / m : 0;
    const yegI = mI > 0 ? y[I_IEG] / mI : 0;
    const eIv = mIv > 0 ? mIv * yegI : mIv * yegC;
    dy[this.iCEG] = eIv - (mEv > 0 ? mEv * yegC : 0);
    dy[I_IEG] = acc ? dy[I_IEG] - eIv : -eIv;
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
    // plenum valve-side flows (the engine began the plenum rates; boundary flows are added by it)
    ip.addFlow(-mIv, sc);
    ep.addFlow(mEv, sc);
    dy[this.iMIVI] = mIv > 0 ? mIv : 0;
    dy[this.iMIVO] = mIv < 0 ? -mIv : 0;
    dy[this.iMEVO] = mEv > 0 ? mEv : 0;
    dy[this.iMEVI] = mEv < 0 ? -mEv : 0;
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
      if (opts.swirlMomentumEfficiency > 0) swirlIn += opts.swirlMomentumEfficiency * this.ivFlow.inflowJet(mIv, v, this.jet).angularMomentumFlux;
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
    ti.L = integralLengthScale(this.h, e.spec.geometry.bore, e.lengthFraction);
    ti.dlnRhoDt = (mIn - mOut) / m - this.Vdot / this.V;
    ti.swirlTorqueIn = swirlIn;
    ti.rho = sc.rho;
    ti.mu = this.viscosityU(sc.T); // swirl wall friction only: viscosity of the last trapped charge
    ti.bore = e.spec.geometry.bore;
    ti.h = this.h;
    this.L = ti.L;
    this.uPrime = turbulenceIntensity(y[this.iTKE], m);
    this.SL = 0;
    this.mdotB = 0;
    this.mdotE = 0;
    this.dpdt = 0;
  }

  private evalClosed(t: number, y: Float64Array, dy: Float64Array, th: number): void {
    const e = this.e;
    const opts = e.opts;
    const P = opts.profile;
    const cl = this.closure;
    // volume of the two zones: the cylinder minus the crevice zone
    const V = this.creviceOn ? this.V - this.Vcr : this.V;
    const mu = y[this.iMU];
    const mb = y[this.iMB];
    const m = mu + mb;
    const U = y[this.iUT];
    this.mCyl = m;
    this.mdotIv = 0;
    this.mdotEv = 0;
    this.stiffCyl = 0;
    this.aIv = 0;
    this.aEv = 0;
    // ---- closure ----
    const tc = P ? now() : 0;
    const e0 = cl.eqSolves;
    cl.noCoupled = e.noCoupled;
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
      if (!cl.solveTwoZone(U, y[this.iSU], mu, mb, V, pG > 0 ? pG : this.pGuess, TbG > 300 ? TbG : this.TbGuess, this.TuGuess)) {
        // robust retry from the last converged pressure and the HP flame temperature
        const TbR = cl.flameTemperatureHP(this.eqAux, this.TuGuess, this.pGuess);
        if (!cl.solveTwoZone(U, y[this.iSU], mu, mb, V, this.pGuess, TbR, this.TuGuess)) this.closureFailures++;
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
      e.prof.closure += now() - tc;
      e.prof.eqSolves += cl.eqSolves - e0;
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
    const bore = e.spec.geometry.bore;
    let L = integralLengthScale(this.h, bore, e.lengthFraction);
    if (mode === MODE_TWO && this.rhoRef > 0 && this.rhoU > 0) {
      const La = angularMomentumLengthScale(this.Lref, this.rhoRef, this.rhoU);
      if (La < L) L = La;
    }
    this.L = L;
    const mTurb = mode === MODE_TWO ? mu : m;
    const up = turbulenceIntensity(y[this.iTKE], mTurb);
    this.uPrime = up;
    this.SL = 0;
    if (mode !== MODE_BURNED && mu > 0) this.SL = laminarFlameSpeed(this.fuelTrapped, this.phiTrapped, this.Tu, p, this.xDil);
    this.lambda = opts.taylorScaleMultiplier * taylorMicroscale(L, up, this.nuU);
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
    if (mode === MODE_TWO) {
      // equivalent radius of the burned gas in the chamber (flame-centred sphere ∩ disc)
      this.rb = this.flameGeom.radiusForVolume(Vb * discScale, this.h, this.clampRadiusGuess(this.rbGuess, this.h));
      this.rbGuess = this.rb;
    }
    if (mode === MODE_TWO && this.flameActive) {
      const me = y[this.iME];
      const Ve = Vb + (me > mb ? (me - mb) * vu : 0);
      this.rf = this.flameGeom.radiusForVolume(Ve * discScale, this.h, this.clampRadiusGuess(this.rfGuess, this.h));
      this.rfGuess = this.rf;
      this.flameGeom.evaluate(this.rf, this.h, this.fg);
      this.Af = this.fg.frontArea > 0 ? this.fg.frontArea : 0;
      let uT = opts.burnRateMultiplier * up;
      if (e.keck && this.rhoInlet > 0) {
        // Keck 1982 eq. 4.10 / Fig. 15 empirical closures (entrainment.ts)
        uT = opts.burnRateMultiplier * keckCharacteristicSpeed(this.keckInletSpeed, this.rhoU, this.rhoInlet);
        this.lambda = opts.taylorScaleMultiplier * keckCharacteristicLength(e.spec.intakeValve.maxLift, this.rhoU, this.rhoInlet);
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
    } else if (mode === MODE_TWO && opts.combustionModel === 'wiebe') {
      const w = opts.wiebe;
      const u = (th - this.wiebeStart) / w.durationDeg;
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
      const eg = m - y[this.iME];
      const egc = eg > 0 ? (eg < mu ? eg : mu) : 0;
      mdK = endGasBurnRate(opts.knockBrushAutoignition ? mu : egc, this.tauAB);
      mdB += mdK;
      mdE += endGasBurnRate(egc, this.tauAB);
      dy[this.iMK] = mdK;
    }
    this.mdotB = mdB;
    this.mdotE = mdE;
    if (P) e.prof.flame += now() - tf;
    // ---- heat transfer ----
    const th0 = P ? now() : 0;
    let Qu = 0;
    let Qb = 0;
    if (opts.heatTransfer) {
      const w = this.woschni;
      w.bore = bore;
      w.pressure = p;
      w.temperature = this.T;
      w.meanPistonSpeed = this.meanPistonSpeed();
      w.variant = opts.woschniVariant;
      if (mode === MODE_SINGLE) {
        w.phase = 'compression';
        this.pMot = p;
      } else {
        w.phase = 'combustion';
        this.pMot = this.motoredPressure(mode === MODE_TWO ? y[this.iSU] / mu : this.suLast, m, V);
        w.motoredPressure = this.pMot;
        w.c2 = opts.woschniCombustionTermMultiplier * WOSCHNI_CONSTANTS.c2Combustion;
        w.displacedVolume = this.kin.displacedVolume;
        w.refPressure = this.pIvc;
        w.refTemperature = this.TIvc;
        w.refVolume = this.VIvc;
      }
      const hc = opts.woschniMultiplier * (opts.heatTransferCorrelation === 'hohenberg' ? hohenbergCoefficient(this.V, p, this.T, w.meanPistonSpeed) : woschniCoefficient(w));
      this.hcoef = hc;
      flatChamberAreas(bore, this.h, e.spec.intakeValve, e.spec.exhaustValve, this.areas);
      if (mode === MODE_TWO) {
        this.flameGeom.evaluate(this.rb, this.h, this.fgB);
        wallHeatLossTwoZone(hc, this.Tu, this.Tb, this.areas, this.fgB, this.walls, this.heat);
        Qu = this.heat.unburned;
        Qb = this.heat.burned;
        const q0 = Qu + Qb;
        // numerical guard for vanishing zones (see options.zoneHeatLossMinTime)
        const tg = opts.zoneHeatLossMinTime;
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
    if (P) e.prof.heat += now() - th0;
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
      const mcu = y[this.iCRU];
      const mcb = y[this.iCRB];
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
        dy[this.iCRU] = (1 - fb) * md;
        dy[this.iCRB] = fb * md;
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
        dy[this.iCRU] = -yu * md;
        dy[this.iCRB] = -(1 - yu) * md;
        Qcr = -md * (yu * Ru + (1 - yu) * Rcc) * Tcr;
        crRet = toU * md;
        this.mdotCr = -md;
      }
      dy[this.iQS + 1] += 0.5 * Qcr;
      dy[this.iQS + 2] += 0.5 * Qcr;
    }
    this.Qcr = Qcr;
    this.Qwall = Qu + Qb + Qcr;
    dy[this.iUT] = dU;
    dy[this.iMU] = dmu;
    dy[this.iMB] = dmb;
    dy[this.iME] = dme;
    if (mode === MODE_TWO) {
      dy[this.iSU] = cl.su * dmu + dSnet;
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
    const i = this.iQS;
    dy[i] = f * q.head;
    dy[i + 1] = f * q.piston;
    dy[i + 2] = f * q.liner;
    dy[i + 3] = f * q.intakeValves;
    dy[i + 4] = f * q.exhaustValves;
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
    const R = this.e.spec.walls.thermalResistance;
    if (!R || !this.wallsLumped || !(dtCycle > 0)) return;
    const y = this.e.y;
    const tau = this.e.opts.wallTimeConstant;
    const a = instant || !(tau > 0) ? 1 : -Math.expm1(-dtCycle / tau);
    const rs = [R.head, R.piston, R.liner, R.intakeValve, R.exhaustValve];
    for (let i = 0; i < 5; i++) {
      const qBar = (y[this.iQS + i] - this.qsCycleStart[i]) / dtCycle;
      let target = rs[i] * qBar;
      if (target < -50) target = -50;
      if (target > 600) target = 600;
      if (Number.isFinite(target)) this.wallExcess[i] += a * (target - this.wallExcess[i]);
    }
  }

  /** Mean piston speed at the current ω (accepted state), m/s. */
  private meanPistonSpeed(): number {
    return (2 * this.e.spec.geometry.stroke * Math.abs(this.e.y[I_OM])) / (2 * Math.PI);
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

  /** Operator splits and per-step accumulators after an accepted RK step of length h (engine state y). */
  afterStep(h: number): void {
    const e = this.e;
    const y = e.y;
    const y0 = e.y0;
    const mbStart = y0[this.iMB];
    const muStart = y0[this.iMU];
    const P = e.opts.profile;
    let changed = false;
    this.noMergedMass = 0;
    const mbRK = y[this.iMB]; // burned mass at the end of the RK step, before the splits
    // End-of-RK-step state BEFORE the splits: the end-gas burn of this step happened in the
    // unburned zone at this temperature (a burn-out merge below sets T_u = 0 in MODE_BURNED —
    // evaluating the knock source there gave heatOfReaction(0) = NaN; validation round 1).
    const TuStep = this.Tu > 0 ? this.Tu : this.TIvc;
    const pStepEnd = this.p;
    const dpStepEnd = this.dpdt;
    // autoigniting gas at the step start / end (its decay shapes the knock source): all unburned gas
    // (knockBrushAutoignition) or the end gas ahead of the front
    const brush = e.opts.knockBrushAutoignition;
    const eg0 = brush ? y0[this.iMU] : y0[this.iMU] + y0[this.iMB] - y0[this.iME];
    const eg1 = brush ? y[this.iMU] : y[this.iMU] + y[this.iMB] - y[this.iME];
    // ---- 1. ignition ----
    if (this.ign) {
      const ti = P ? now() : 0;
      changed = this.ignitionSplit(h) || changed;
      if (P) e.prof.ignition += now() - ti;
    }
    // ---- 2. end-gas burn-up: inside the RK right-hand side (dm_b/dt += (m − m_e)/τ_ab with steps ≤
    // τ_ab/4); the mass it burned over this step drives the acoustic modes below ----
    const dmKnock = this.knockBurning ? y[this.iMK] - y0[this.iMK] : 0;
    if (this.knockBurning && this.mode === MODE_TWO) {
      // end of the autoignition burn-up: the unburned gas is consumed; the e^{−t/τ_ab} tail below 1e-9
      // of the charge is moved at constant U_tot (the burn-out merge then follows).
      const m = y[this.iMU] + y[this.iMB];
      const left = e.opts.knockBrushAutoignition ? y[this.iMU] : Math.min(y[this.iMU], m - y[this.iME]);
      if (left <= 1e-9 * m) {
        if (left > 0) {
          this.transferToBurned(left);
          changed = true;
        }
        if (this.mode === MODE_TWO && y[this.iME] < y[this.iMU] + y[this.iMB]) y[this.iME] = y[this.iMU] + y[this.iMB];
        this.knockBurning = false;
      }
    }
    // ---- 3. burn-out ----
    if (this.mode === MODE_TWO) {
      const mu = y[this.iMU];
      const m = mu + y[this.iMB];
      if (mu <= e.opts.burnoutFraction * m) {
        this.mergeToBurned();
        changed = true;
      }
    }
    if (changed) e.refreshCylinder(this);
    // ---- 4. NO (burned zone): Strang splitting — half a step of the rate-controlled kinetics on the
    // step-START burned-zone state, the gas burned during the step enters at mid-step with the unburned
    // NO, half a step on the step-END state (second order in h; round 1 advanced the whole step on the
    // end state with the new gas added at the start: first order, 4.6 % in-cycle error at 0.25° and a
    // summary NO that moved 0.5 % with the snapshot cadence — validation round 2) ----
    if ((this.mode === MODE_TWO || this.mode === MODE_BURNED) && y[this.iMB] > 0) {
      const tn = P ? now() : 0;
      // burned gas that entered the crevice this step takes its NO along (frozen there, returned as
      // complete-combustion products — the crevice's NO is not tracked)
      const dcrb = y[this.iCRB] - y0[this.iCRB];
      if (dcrb > 0 && y[this.iMB] > 0) this.nNO *= y[this.iMB] / (y[this.iMB] + dcrb);
      // (mass merged by mergeToBurned in this step already brought its NO)
      const dmb = y[this.iMB] - mbStart - this.noMergedMass;
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
      if (P) e.prof.no += now() - tn;
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
      const dsMax = e.opts.knockBurnStep;
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
        const F0 = seq ? (y0[this.iMK] - this.knockMkOnset) / this.knockEgOnset : 0;
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
      if (P) e.prof.knock += now() - tk;
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
        if (this.tracing() && !(this.twoStage && this.lwStage === 1)) this.recordSample(this.p + this.knockOscillation());
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
    if (!e.derivValid) e.refreshCylinder(this);
    const dmbStep = y[this.iMB] - mbStart;
    this.burnRateStep = this.mode !== MODE_OPEN && dmbStep > 0 ? dmbStep / h : 0;
    const Tq = this.Tu > 0 ? this.Tu : TuStep;
    // heat-release rate for output: the INSTANTANEOUS burn rate of the end-of-step right-hand side
    // plus the step mean of the split transfers (kernel mass, burn-out merge) where they act — the
    // round-1 step mean depended on the step and sampling pattern by 6–11 % (validation round 2)
    const splitRate = this.mode !== MODE_OPEN && y[this.iMB] > mbRK ? (y[this.iMB] - mbRK) / h : 0;
    const rate = (this.mode === MODE_TWO ? this.mdotB : 0) + splitRate;
    this.hrr = rate > 0 ? rate * this.closure.heatOfReaction(Tq) : 0;
    this.bookkeeping(muStart, pStepEnd, dpStepEnd);
  }

  /** True when this cylinder records the engine's validation trace now (cylinder 0 only). */
  private tracing(): boolean {
    const e = this.e;
    return this.index === 0 && e.trace !== null && e.cycle === e.traceCycle;
  }

  /** Effective sound speed of the charge for the knock modes (sets knockGamma), m/s. */
  private knockSoundSpeed(): number {
    const cl = this.closure;
    const y = this.e.y;
    if (this.mode === MODE_TWO) {
      const Vb = y[this.iMB] * cl.vb;
      const Vu = y[this.iMU] * cl.vu;
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
    const dm = this.e.delayModel;
    const s = (Math.log(dm.tau(T + 1, p, phi, f, x)) - Math.log(dm.tau(T - 1, p, phi, f, x))) / 2;
    const g = s * T;
    return Number.isFinite(g) ? g : 0;
  }

  /**
   * The Livengood–Wu integral runs on the end gas AHEAD of the flame front: from IVC while the knock
   * model is available for the trapped fuel, until onset — or until the front has swept the whole
   * charge (no end gas left: the remaining brush behind the front burns up on τ_b and is not "end
   * gas"; validation round 1 found brush/crevice remainders "autoigniting" late in expansion).
   */
  lwArmed(): boolean {
    const opts = this.e.opts;
    return (
      opts.knock &&
      this.knockAvailable &&
      !this.knockOnset &&
      (this.mode === MODE_SINGLE || this.mode === MODE_TWO) &&
      // (with only the end gas autoigniting, the integral stops once the front has swept the chamber;
      // with the brush pockets autoigniting too it runs until burn-out, so the onset moves
      // continuously with the operating point instead of appearing with a finite brush mass)
      (opts.knockBrushAutoignition || !(this.flameActive && this.frontAtWalls)) &&
      !this.e.finished
    );
  }

  /** End-gas ignition delay (total, hot ignition) at the current evaluation (unburned zone, trapped mixture), s. */
  private tauTotal(): number {
    return this.e.delayModel.tau(this.Tu, this.p, this.phiTrapped, this.fuelTrapped, this.xResMole);
  }

  /**
   * Delay governing the ACTIVE Livengood–Wu integral, s: the total τ ('single'), or for the
   * two-stage integral τ₁ (first stage) then τ − τ₁ (second stage, ≥ 0: 0 where the chemistry is
   * single-stage, τ₁ = τ, i.e. hot ignition follows the first-stage crossing immediately).
   */
  tauNow(): number {
    if (!this.twoStage) return this.tauTotal();
    const m = this.e.delayModel as TabulatedIgnitionDelay;
    const t1 = m.tauFirstStage(this.Tu, this.p, this.phiTrapped, this.fuelTrapped, this.xResMole);
    if (this.lwStage === 1) return t1;
    const d = this.tauTotal() - t1;
    return d > 0 ? d : 0;
  }

  /** Fraction of the tentative step at which the LW integral reaches 1 (−1 if not). */
  lwCrossingFraction(h: number, tau1: number): number {
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
    const e = this.e;
    const y = e.y;
    const mu = y[this.iMU];
    const m = mu + y[this.iMB];
    if (!(mu > 0)) return;
    // autoigniting gas: all unburned gas (end gas ahead of the front + brush pockets; see evalClosed),
    // or only the end gas ahead of the front (options.knockBrushAutoignition false; round 1)
    const eg = e.opts.knockBrushAutoignition || this.mode === MODE_SINGLE ? mu : Math.min(mu, m - y[this.iME]);
    if (!(eg > 0)) return;
    this.knockOnset = true;
    this.knockBurning = true;
    e.interrupt = true;
    this.knockOnsetDeg = this.theta;
    this.tKnockOnset = e.t;
    this.knockEndGasFraction = eg / m;
    const tau = this.tauTotal(); // total (hot-ignition) delay: J is built on it (two-stage too)
    const Ton = this.Tu;
    this.tauAB = autoignitionBurnTime(tau, Ton > 0 ? this.lwJ / Ton : 0, e.opts.knockStratificationDT, e.opts.knockExcitationTime);
    if (!(this.tauAB > 0)) this.tauAB = e.opts.knockExcitationTime;
    // end-gas planform region with the end-gas volume fraction
    const Vg = this.creviceOn ? this.V - this.Vcr : this.V;
    const vEg = this.mode === MODE_SINGLE ? Vg : eg * this.closure.vu;
    const c = e.spec.sparkPlug.gapCenter;
    const rIn = this.endGasCircleRadius(Math.min(1, vEg / Vg));
    if (e.opts.knockSourceShells > 0) this.knockOsc.setSequentialEndGasRegion(c[0], c[2], rIn, e.opts.knockSourceShells);
    else this.knockOsc.setEndGasRegion(c[0], c[2], rIn);
    this.knockEgOnset = eg;
    this.knockMkOnset = y[this.iMK];
    // autoignition before any burned gas (e.g. before the spark): seed the burned zone
    if (this.mode === MODE_SINGLE) this.transferToBurned(1e-9 * mu);
    e.invalidate(this);
    e.refreshCylinder(this);
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
    const R = 0.5 * this.e.spec.geometry.bore;
    const g = this.e.spec.sparkPlug.gapCenter;
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
    const e = this.e;
    const ign = this.ign!;
    const y = e.y;
    const th = this.theta;
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
    gas.flowVelocity = meanFlowVelocity(y[this.iTK], this.mode === MODE_TWO ? y[this.iMU] : this.mCyl);
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
      gas.marksteinLength = e.opts.marksteinMultiplier * this.mk.unburned;
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
          const mu = y[this.iMU];
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
        if (this.mode === MODE_TWO && y[this.iMU] > 0) {
          if (changed) e.refreshCylinder(this); // T_b after the kernel-mass transfer
          const Tk = st.kernel.temperature > 0 ? st.kernel.temperature : this.Tb;
          const cap = y[this.iMB] * this.closure.cvb * (Tk - this.Tb); // J to bring T_b to T_k
          // a net LOSS (the kernel's electrode conduction after the discharge) cools the kernel gas =
          // the burned zone, at most down to T_u (round 1 took it from the unburned zone; validation
          // round 2 — magnitude < 0.01 K)
          const capLoss = -y[this.iMB] * this.closure.cvb * Math.max(0, this.Tb - this.Tu);
          const toB = net > 0 ? (cap > 0 ? Math.min(net, cap) : 0) : Math.max(net, capLoss);
          const toU = net - toB;
          if (toU !== 0 && this.Tu > 0) y[this.iSU] += toU / this.Tu;
        }
        y[this.iUT] += net;
      } else y[this.iCU] += net;
      this.sparkEnergy += net;
      if (dL > 0) this.electrodeLoss += dL;
      changed = true;
    }
    // (3) hand-off to the entrainment model
    if (handoff) {
      if (changed) e.refreshCylinder(this);
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
        e.refreshCylinder(this);
        const m = y[this.iMU] + y[this.iMB];
        const dme = this.mdotE * late;
        const dmb = Math.min(this.mdotB * late, 0.5 * y[this.iMU]);
        this.transferToBurned(dmb);
        y[this.iME] = Math.min(m, Math.max(y[this.iME] + dme, y[this.iMB]));
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
    const y = this.e.y;
    y[this.iME] = this.handoffEntrainedMass(rk, STk);
    this.rfGuess = rk;
    this.flameActive = true;
    this.e.invalidate(this);
  }

  /**
   * Entrained mass m_e that the entrainment model starts from when the kernel of radius rk and
   * burning velocity STk is handed off now (see onHandoff); written to this.fgR's radius as a side
   * product (kernel-stage flame radius, flameRadius()).
   */
  private handoffEntrainedMass(rk: number, STk: number): number {
    const y = this.e.y;
    const cl = this.closure;
    const mb = y[this.iMB];
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
    const mt = y[this.iMU] + mb;
    return me < mt ? me : mt;
  }

  /** Move Δm from the unburned to the burned zone at constant U_tot (creates the burned zone if needed). */
  private transferToBurned(dm: number): void {
    const y = this.e.y;
    if (!(dm > 0)) return;
    if (this.mode === MODE_SINGLE) this.createBurnedZone();
    if (this.mode !== MODE_TWO) return;
    const mu = y[this.iMU];
    if (dm >= mu) {
      this.mergeToBurned();
      return;
    }
    const su = y[this.iSU] / mu;
    y[this.iSU] -= su * dm;
    y[this.iMU] = mu - dm;
    y[this.iMB] += dm;
    if (y[this.iME] < y[this.iMB]) y[this.iME] = y[this.iMB];
    this.e.invalidate(this);
  }

  /** Single zone → two zones: S_u from the current single-zone state; T_b seeded at the HP flame state. */
  private createBurnedZone(): void {
    const e = this.e;
    const y = e.y;
    const cl = this.closure;
    const mu = y[this.iMU];
    const T = cl.solveSingle(y[this.iUT], mu, this.creviceOn ? this.V - this.Vcr : this.V, this.TuGuess);
    const p = cl.p;
    y[this.iSU] = mu * cl.entropyAt(T, p);
    y[this.iMB] = 0;
    this.mode = MODE_TWO;
    this.pGuess = p;
    this.TuGuess = T;
    this.TbGuess = cl.flameTemperatureHP(this.eqAux, T, p);
    this.dpEval = 0;
    this.dTbEval = 0;
    this.tEval = e.t;
    this.Lref = integralLengthScale(this.h, e.spec.geometry.bore, e.lengthFraction);
    this.rhoRef = p / (cl.Ru * T);
    this.rbGuess = 0;
    e.invalidate(this);
  }

  /** Two zones → burned only (the unburned remainder is burned at constant U_tot). */
  private mergeToBurned(): void {
    const y = this.e.y;
    if (this.mode === MODE_SINGLE) this.createBurnedZone();
    const mu = y[this.iMU];
    const m = mu + y[this.iMB];
    if (mu > 0) {
      this.suLast = y[this.iSU] / mu;
      // turbulence: continue on the whole charge with the unburned specific values
      y[this.iTKE] *= m / mu;
      y[this.iTK] *= m / mu;
      this.nNO += mu * this.closure.noPerKg;
      this.noMergedMass += mu;
    }
    y[this.iMB] = m;
    y[this.iMU] = 0;
    y[this.iSU] = 0;
    y[this.iME] = m;
    this.mode = MODE_BURNED;
    this.burnDone = true;
    this.flameActive = false;
    this.knockBurning = false;
    this.e.invalidate(this);
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
        const y = this.e.y;
        const Ve = y[this.iMB] * cl.vb + (me - y[this.iMB]) * cl.vu;
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
    const y = this.e.y;
    const m = y[this.iMU] + y[this.iMB];
    return y[this.iME] >= m * (1 - 1e-6) || this.rf >= this.flameGeom.maxRadius(this.h) * (1 - 1e-6) || this.Af <= 1e-10;
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
  private bookkeeping(muStart: number, pEnd: number, dpdtEnd: number): void {
    const e = this.e;
    const y = e.y;
    const th = this.theta;
    const pRep = this.p + this.knockOscillation();
    const omDeg = y[I_OM] * RAD2DEG;
    const th0 = this.localAngle(e.y0[I_TH]);
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
      const m = y[this.iMU] + y[this.iMB];
      const x1 = y[this.iMB] / m;
      const x0 = 1 - muStart / m;
      if (Number.isNaN(this.ca10)) this.ca10 = crossing(x0, x1, th0, th, 0.1);
      if (Number.isNaN(this.ca50)) this.ca50 = crossing(x0, x1, th0, th, 0.5);
      if (Number.isNaN(this.ca90)) this.ca90 = crossing(x0, x1, th0, th, 0.9);
    }
    if (this.tracing()) this.recordSample(pRep);
  }

  /** Append the current state to the engine's validation trace (cylinder 0). */
  recordSample(pRep: number): void {
    const e = this.e;
    const tr = e.trace!;
    const y = e.y;
    const closed = this.mode !== MODE_OPEN;
    const mu = closed ? y[this.iMU] : this.mCyl;
    const mb = closed ? y[this.iMB] : 0;
    tr.theta.push(this.theta);
    tr.t.push(e.t);
    tr.volume.push(this.V);
    tr.pressure.push(this.p);
    tr.pressureReported.push(pRep);
    tr.Tu.push(this.Tu);
    tr.Tb.push(this.Tb);
    tr.Tmean.push(this.T);
    tr.mu.push(mu);
    tr.mb.push(mb);
    tr.xb.push(closed ? mb / (mu + mb) : 0);
    tr.me.push(closed ? y[this.iME] : 0);
    tr.U.push(closed ? y[this.iUT] + this.creviceEnergy() : y[this.iCU]);
    tr.work.push(y[this.iW]);
    tr.heatLoss.push(y[this.iQ]);
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
    tr.nNOClosure.push(e.noCoupled && closed ? this.closure.noKinetic : 0);
    tr.mdotIntake.push(this.mdotIv);
    tr.mdotExhaust.push(this.mdotEv);
    tr.pIntake.push(e.pInt);
    tr.pExhaust.push(e.pExh);
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

  /**
   * Events of the current local cycle up to the local angle th (the engine landed on one). Returns true
   * when the ENGINE wrapped (cylinder 0's EV_WRAP at θ = 360).
   */
  handleEvents(th: number): boolean {
    const e = this.e;
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
          this.wBdcStart = e.y[this.iW];
          break;
        case EV_BDC_END:
          this.wBdcEnd = e.y[this.iW];
          break;
        case EV_TDC:
          if (this.mode === MODE_SINGLE || this.mode === MODE_TWO) {
            const x0 = this.xb;
            this.mergeToBurned();
            e.refreshCylinder(this);
            this.markBurnJump(x0, this.xb, th);
          }
          break;
        case EV_IGN:
          if (Math.abs(this.evAngle[this.evNext - 1] - this.ignCmd.sparkDeg) < 1e-9 && Number.isNaN(this.tSparkCmd)) {
            this.tSparkCmd = e.t;
            e.interrupt = true;
          }
          break;
        case EV_WIEBE:
          if (this.mode === MODE_SINGLE) {
            // Seed of the burned zone: WIEBE_SEED of the charge (numerical; x_b starts at it). The
            // two-zone closure recovers T_b from U_tot − m_u u_u, so dT_b/dp ≈ (V_u/γ)/(m_b c_v): with
            // the round-1 seed 1e-8 (7e-12 kg) p had to be exact to 1e-4 Pa and 2–3 solves failed per
            // cycle (T_b pinned at a bound; validation round 2). At 1e-5 the requirement is 1e-9 of p.
            const x0 = this.xb;
            this.transferToBurned(WIEBE_SEED * e.y[this.iMU]);
            e.refreshCylinder(this);
            this.markBurnJump(x0, this.xb, th);
          }
          break;
        case EV_END:
          this.onEvo();
          e.finished = true;
          break;
        case EV_WRAP:
          if (this.index === 0) {
            e.engineWrap();
            return true;
          }
          this.onLocalWrap();
          return false;
        default:
          break;
      }
    }
    return false;
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
    const e = this.e;
    const y = e.y;
    const N = this.Nscr;
    let n = 0;
    let m = 0;
    let mf = 0;
    for (let k = 0; k < NS; k++) {
      N[k] = y[this.iCN + k];
      n += N[k];
      m += N[k] * MOLAR_MASS[k];
    }
    for (const s of FUEL_SPECIES) mf += N[SP[s]] * MOLAR_MASS[SP[s]];
    for (let k = 0; k < NS; k++) this.Xscr[k] = N[k] / n;
    const cl = this.closure;
    cl.setUnburned(this.Xscr);
    const U = y[this.iCU];
    y[this.iUT] = U;
    y[this.iMU] = m;
    y[this.iMB] = 0;
    y[this.iME] = 0;
    y[this.iSU] = 0;
    y[this.iCRU] = 0;
    y[this.iCRB] = 0;
    this.mode = MODE_SINGLE;
    const V = this.kin.volume(this.theta * DEG);
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
      y[this.iUT] = U - mcr * this.uuCr;
      y[this.iMU] = m - mcr;
      y[this.iCRU] = mcr;
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
    this.rhoInlet = e.intake.state.rho;
    {
      const yb = m > 0 ? Math.min(1, Math.max(0, (y[this.iCBG] - y[this.iCEG]) / m)) : 0;
      const ev = (m * (1 - yb)) / (this.rhoInlet * this.kin.displacedVolume);
      const Aiv = this.ivFlow.flowArea(e.spec.intakeValve.maxLift);
      this.keckInletSpeed = keckMeanInletSpeed(ev, this.kin.boreArea, Aiv, (y[I_OM] * 60) / (2 * Math.PI), e.spec.geometry.stroke);
    }
    // residual / dilution
    // burned-gas (products) fraction of the TRAPPED charge — residual + external EGR — is the dilution
    // of S_L and the ignition delay (round 1 added the REQUESTED op.egrFraction: 2.2× too much in the
    // first cycle after an EGR step, which then misfired; validation round 2)
    this.yRes = m > 0 ? Math.min(1, Math.max(0, y[this.iCBG] / m)) : 0;
    this.yEgr = m > 0 ? Math.min(this.yRes, Math.max(0, y[this.iCEG] / m)) : 0;
    this.xDil = this.yRes;
    const fa = freshCharge({ fuel: this.fuelTrapped, phi: this.phiTrapped > 0 ? this.phiTrapped : 1, airX: e.airX });
    const Mfa = mixMolarMass(fa);
    const Mres = mixMolarMass(completeCombustionProducts(fa));
    this.xResMole = residualMoleFraction(this.xDil, Mres, Mfa);
    // unburned viscosity table
    for (let i = 0; i < this.muTab.length; i++) this.muTab[i] = mixViscosity(cl.Xu, this.muT0 + i * this.muDT);
    this.nuU = (this.viscosityU(T) * V) / m;
    // knock / NO / flame state
    this.lw.reset();
    this.lwStage = 1;
    this.twoStage = e.opts.knockIntegral === 'two-stage' && e.delayModel instanceof TabulatedIgnitionDelay;
    this.lwTauPrev = NaN;
    this.lwJ = 0;
    this.lwRPrev = NaN;
    this.lwGPrev = 0;
    this.knockOsc.reset();
    this.nNO = 0;
    this.qIvc = y[this.iQ];
    this.closedHappened = true;
    this.rfGuess = 0;
    this.rbGuess = 0;
    e.derivValid = false;
    e.refresh();
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
    const sel = this.e.fuel;
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
    const airX = this.e.airX;
    const rO2N2 = airX[SP.O2] / airX[SP.N2];
    const nO2 = 0.5 * b[EL.N] * rO2N2;
    this.phiTrapped = nO2 > 0 ? Math.max(0, 1 + (2 * b[EL.C] + 0.5 * b[EL.H] - b[EL.O]) / (2 * nO2)) : 0;
  }

  /** EVO: merge the zones (burned composition frozen, NO → kinetic value) into the gas-exchange zone. */
  private onEvo(): void {
    const e = this.e;
    const y = e.y;
    const cl = this.closure;
    const mu = y[this.iMU];
    const mb = y[this.iMB];
    const m = mu + mb;
    const N = this.Nscr;
    for (let k = 0; k < NS; k++) N[k] = mu * cl.nu[k];
    // crevice gas back into the gas-exchange zone (unburned: frozen charge; burned: complete-combustion
    // products), with its energy at T_cr
    const mcu = y[this.iCRU];
    const mcb = y[this.iCRB];
    if (mcu > 0 || mcb > 0) for (let k = 0; k < NS; k++) N[k] += mcu * cl.nu[k] + mcb * cl.ncc[k];
    if (mb > 0 && this.mode !== MODE_SINGLE) {
      // burned zone at its equilibrium state (T_b, p), NO replaced by the kinetic value
      const bb = this.bScr;
      for (let el = 0; el < NE; el++) bb[el] = cl.bu[el] * mb;
      const r = cl.eq.solveTP(bb, this.Tb, this.p);
      const Nb = this.Nb;
      for (let k = 0; k < NS; k++) Nb[k] = r.N[k];
      swapNO(Nb, this.nNO);
      for (let k = 0; k < NS; k++) N[k] += Nb[k];
    }
    let n = 0;
    for (let k = 0; k < NS; k++) {
      y[this.iCN + k] = N[k];
      n += N[k];
    }
    y[this.iCU] = y[this.iUT] + this.creviceEnergy();
    y[this.iCBG] = mb + mcb + (mu + mcu) * this.yRes;
    y[this.iCEG] = 0; // after the cycle everything in the cylinder is residual
    y[this.iCRU] = 0;
    y[this.iCRB] = 0;
    // turbulence on the whole charge
    if (this.mode === MODE_TWO && mu > 0) {
      y[this.iTKE] *= m / mu;
      y[this.iTK] *= m / mu;
    }
    this.noPpm = (N[SP.NO] / n) * 1e6;
    this.coFrac = N[SP.CO] / n;
    this.xbEvo = m > 0 ? mb / m : 0;
    if (mb > 0) this.burnDone = this.burnDone || this.xbEvo > 0.999;
    this.qEvo = y[this.iQ];
    this.knockOsc.reset();
    this.lw.reset();
    this.mode = MODE_OPEN;
    this.flameActive = false;
    this.knockBurning = false;
    e.stiffness = 0;
    e.derivValid = false;
    e.refresh();
  }

  /** End of this cylinder's local cycle (θ_i = 360): summary, lumped walls, surface ledgers. */
  endCycle(): CycleSummary | null {
    const e = this.e;
    const y = e.y;
    const s = this.partialCycle ? null : this.makeSummary();
    if (s) e.summaries.push(s);
    this.updateWalls(e.t - this.tCycleStart, e.warmingUp);
    for (let i = 0; i < 5; i++) this.qsCycleStart[i] = y[this.iQS + i];
    this.tCycleStart = e.t;
    this.partialCycle = false;
    return s;
  }

  /** Apply the engine's compression ratio when it differs from this cylinder's (at θ_i = 360, before the wrap). */
  applyCompressionRatio(): void {
    const cr = this.e.op.compressionRatio;
    if (cr !== this.kin.compressionRatio) this.changeCompressionRatio(cr);
  }

  /** New local cycle at θ_i = −360 (after the angle wrap): per-cylinder data, flags, events, ignition command. */
  beginCycle(): void {
    this.prepareCycleData();
    this.cycle++;
    this.resetCycleFlags();
    this.startCycle();
  }

  /** Local wrap of a cylinder i ≥ 1 (θ_i = 360 → −360; the engine angle runs on). */
  private onLocalWrap(): void {
    const e = this.e;
    this.endCycle();
    this.applyCompressionRatio();
    this.angleShift -= 720;
    this.beginCycle();
    e.derivValid = false;
    e.refresh();
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
    const y = this.e.y;
    const N = this.Nscr;
    const X = this.Xscr;
    let n = 0;
    for (let k = 0; k < NS; k++) {
      N[k] = y[this.iCN + k];
      n += N[k];
    }
    const th = this.theta;
    const V1 = this.kin.volume(th * DEG);
    this.kin.setCompressionRatio(cr);
    const V2 = this.kin.volume(th * DEG);
    if (!(n > 0) || this.mode !== MODE_OPEN || !(V1 > 0) || Math.abs(V2 / V1 - 1) < 1e-15) return;
    for (let k = 0; k < NS; k++) X[k] = N[k] / n;
    const U1 = y[this.iCU];
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
    y[this.iCU] = U2;
    y[this.iW] -= U2 - U1; // ∫p dV of the isentropic change (−W_on_gas)
    this.cyl.volume = V2;
  }

  /**
   * Volumetric efficiency of the cycle just completed (Heywood's definition, as measured on the
   * CFR): DRY air inducted through the carburettor over the cycle / (ambient dry-air partial
   * density × V_d). Net venturi mass over the cycle × dry-air mass fraction of the fresh stream.
   * (Round 1 used the trapped mass × (1 − burned-residual fraction), which counted retained
   * UNBURNED gas — motored / misfiring cycles — as inducted: η_v ≈ 1.09.) Multi-cylinder: the venturi
   * feeds all cylinders, so this is the engine value over this cylinder's local cycle (total
   * displacement N·V_d).
   */
  private volumetricEfficiency(): number {
    const e = this.e;
    const mVent = e.ventMass() - this.mVentCycleStart;
    return e.volumetricEfficiencyOf(mVent);
  }

  private makeSummary(): CycleSummary {
    const e = this.e;
    const y = e.y;
    const Vd = this.kin.displacedVolume;
    const wNet = y[this.iW] - this.wCycleStart;
    const wGross = this.wBdcEnd - this.wBdcStart;
    const lhv = lowerHeatingValue(this.fuelTrapped);
    this.volEff = e.closedInit ? 0 : this.volumetricEfficiency();
    const mf = this.fuelMassIvc;
    const burned = this.xbEvo;
    const cm = e.opts.combustionModel;
    const misfire = cm === 'entrainment' ? this.misfire || burned < 0.1 : cm === 'none' ? true : burned < 0.1;
    const s: CycleSummary = {
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
    if (e.cylinders.length > 1) s.cylinder = this.index;
    return s;
  }

  // ===========================================================================================
  // Inventories (conservation checks)
  // ===========================================================================================

  /** Internal energy of the crevice zone (closed phase), J. */
  creviceEnergy(): number {
    const y = this.e.y;
    return this.mode === MODE_OPEN ? 0 : y[this.iCRU] * this.uuCr + y[this.iCRB] * this.uccCr;
  }

  /** Add this cylinder's species moles (closed phase: unburned + burned equilibrium + crevice) to out, mol. */
  addSpecies(out: Float64Array): void {
    const y = this.e.y;
    if (this.mode === MODE_OPEN) {
      for (let k = 0; k < NS; k++) out[k] += y[this.iCN + k];
    } else {
      const cl = this.closure;
      for (let k = 0; k < NS; k++) out[k] += (y[this.iMU] + y[this.iCRU]) * cl.nu[k] + y[this.iCRB] * cl.ncc[k];
      if (y[this.iMB] > 0) {
        const r = cl.eq.result;
        const s = r.mass > 0 ? y[this.iMB] / r.mass : 0;
        for (let k = 0; k < NS; k++) out[k] += r.N[k] * s;
      }
    }
  }

  /** Add this cylinder's open-zone moles to N (0 in the closed phase; elementInventory). */
  addOpenMoles(N: Float64Array): void {
    const y = this.e.y;
    const open = this.mode === MODE_OPEN;
    for (let k = 0; k < NS; k++) N[k] += open ? y[this.iCN + k] : 0;
  }

  /** Add this cylinder's closed-zone element moles to b (elementInventory). */
  addClosedElements(b: Float64Array): void {
    if (this.mode === MODE_OPEN) return;
    const y = this.e.y;
    const m = y[this.iMU] + y[this.iMB] + y[this.iCRU] + y[this.iCRB];
    for (let el = 0; el < NE; el++) b[el] += m * this.closure.bu[el];
  }

  /** Internal energy of the cylinder charge (+ crevice), J. */
  energy(): number {
    const y = this.e.y;
    return this.mode === MODE_OPEN ? y[this.iCU] : y[this.iUT] + this.creviceEnergy();
  }

  /** Add the cylinder mass to m, kg (same summation order as the single-cylinder inventory). */
  addMass(m: number): number {
    const y = this.e.y;
    if (this.mode === MODE_OPEN) for (let k = 0; k < NS; k++) m += y[this.iCN + k] * MOLAR_MASS[k];
    else m += y[this.iMU] + y[this.iMB] + y[this.iCRU] + y[this.iCRB];
    return m;
  }

  /** Shift the model times held by this cylinder by dt (the engine's warm-up end resets t to 0). */
  shiftTimes(dt: number): void {
    this.tCycleStart += dt;
    this.tEval += dt;
    if (!Number.isNaN(this.tSparkCmd)) this.tSparkCmd += dt;
    if (!Number.isNaN(this.tKnockOnset)) this.tKnockOnset += dt;
  }
}
