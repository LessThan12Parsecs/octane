/**
 * Options of the cycle model (src/physics/cycle): sub-model switches, calibration multipliers and
 * numerical controls. See index.ts for the public API.
 *
 * Calibration policy (DESIGN.md): only genuinely uncertain closures are exposed as multipliers
 * — burn rate (entrainment speed u_T = C·u′ and burn-up length), turbulence length scale,
 * Woschni coefficient, intake-port heat transfer, venturi / valve discharge coefficients, knock
 * options — and ONE parameter set per ENGINE serves all of its operating points (calibration.ts).
 *
 * Defaults are layered (Model T integration): {@link NEUTRAL_CYCLE_MODEL_OPTIONS} (engine-independent
 * switches, numerics and universal closures; the engine-fitted closures at their uncalibrated literature
 * values; no engine hardware) ← the per-engine defaults keyed by EngineSpec.id
 * ({@link ENGINE_CYCLE_OPTION_DEFAULTS}: calibration set, friction inputs, valve lash, crankcase pressure,
 * knock pickup and band) ← the caller's options ({@link resolveCycleOptions}). A spec without an id is
 * treated as a CFR derivative (the pre-registry behaviour); an unregistered id gets the neutral set.
 * {@link DEFAULT_CYCLE_MODEL_OPTIONS} is the resolved CFR F-1 set.
 */
import type { OperatingPoint, LoadSpec } from '../core/operating-point';
import type { EngineSpec } from '../core/engine-spec';
import type { IgnitionSystemOptions } from '../ignition';
import type { WoschniVariant } from '../heat-transfer';
import type { PnhFrictionInputs } from '../mechanics/friction';
import type { IgnitionDelayModel } from '../chemistry/ignition-delay';
import { CFR_CRANKCASE_GAUGE_PRESSURE, CFR_FRICTION, CFR_KNOCK_PICKUP, CFR_RON_CONDITIONS, CFR_VALVE_LASH } from '../engines/cfr';
import { MODEL_T, MODEL_T_FRICTION, MODEL_T_VALVE_LASH } from '../engines/model-t';
import { engineOfSpec } from '../engines';
import { END_GAS_STRATIFICATION_DT, EXCITATION_TIME_PRF, KNOCK_DECAY_TIME, L_HEAD_MAPO_BAND, virtualKnockSensor } from '../chemistry/knock';
import { CFR_CALIBRATION, CFR_KNOCK_DELAY_MODEL, MODEL_T_CALIBRATION, type CalibrationSet } from './calibration';
import { moistAirEnhancementFactor, waterSaturationPressure } from '../thermo/fuels';

/** ISO 5167-4:2022 classical Venturi tube with an "as cast" convergent section, C = 0.984 (fetched
 * summary of ISO 5167-4 §5.5.4, 2026-09-30) — the round-1 default and the upper end of the
 * calibration range of options.venturiDischargeCoefficient (calibrated CFR value 0.6,
 * calibration.ts): the 9/16 in CFR venturi (throat Re ≈ 1e4, fuel-nozzle bridge) is far below the
 * standard's 2e5 ≤ Re range and not an ISO tube. Modelled as an orifice (no diffuser recovery; the
 * throttle module estimates the unrecovered loss as negligible at CFR flows). The neutral
 * (uncalibrated) venturi coefficient. */
export const VENTURI_DISCHARGE_COEFFICIENT = 0.984;

/**
 * Burn-rate model:
 *  - 'entrainment': spark (IgnitionSystem: coil, breakdown, discharge, Herweg–Maly kernel) →
 *    hand-off → Keck/Tabaczynski entrainment + burn-up with the exact flame geometry (default);
 *  - 'instantaneous-at-tdc': the whole charge is burned at θ = 0 at constant volume (UV
 *    equilibrium) — the fuel–air-cycle limit (no spark system);
 *  - 'wiebe': prescribed Wiebe mass-fraction-burned curve (no spark system), see WiebeOptions;
 *  - 'none': motored (no spark, no combustion).
 */
export type CombustionModel = 'entrainment' | 'instantaneous-at-tdc' | 'wiebe' | 'none';

/** End-gas ignition-delay models (chemistry/ignition-delay*.ts). */
export type IgnitionDelayModelId = 'llnl-gasoline-2011' | 'llnl-prf-v2' | 'douaud-eyzat' | 'douaud-eyzat-llnl';

/**
 * Turbulent-flame closure of the entrainment model (Keck 1982 eqs. 4.1A/B):
 *  - 'kk-taylor': u_T = C_T·u′ from the K–k turbulence model, burn-up length = C_λ·Taylor
 *    microscale (Tabaczynski et al. 1977);
 *  - 'keck1982': Keck's own empirical closures fitted to several SI engines, u_T = C_T·0.08 ū_i
 *    (ρ_u/ρ_i)^½ (eq. 4.10, ±10 %) and ℓ_T = C_λ·0.8 L_IV (ρ_i/ρ_u)^¾ (Fig. 15, ±25 %), with the
 *    mean inlet-gas speed ū_i = ε_v (A_p/A_IV) 2NS of the current cycle.
 * C_T = burnRateMultiplier, C_λ = taylorScaleMultiplier.
 */
export type TurbulentFlameClosure = 'kk-taylor' | 'keck1982';

/**
 * End-gas autoignition integral: 'single' = Livengood–Wu on the total delay τ; 'two-stage' =
 * first-stage (cool-flame) integral on τ₁, then the second stage on τ − τ₁ (chemistry/
 * livengood-wu.ts TwoStageLivengoodWu; needs a tabulated model with τ₁ — others fall back to
 * 'single').
 */
export type KnockIntegralMode = 'single' | 'two-stage';

/**
 * Wiebe function x_b = 1 − exp(−a ((θ − θ₀)/Δθ)^{m+1}) (Wiebe 1956/1970, as used by Heywood 1988
 * ch. 9 — UNVERIFIED equation number).
 */
export interface WiebeOptions {
  /** Start of combustion θ₀, crank deg (default: the spark angle −sparkAdvanceDeg). */
  startDeg?: number;
  /** Burn duration Δθ, crank deg (default 50 — UNVERIFIED representative CFR value). */
  durationDeg: number;
  /** Efficiency parameter a (default −ln(0.001) = 6.908: x_b = 0.999 at θ₀ + Δθ; a definition). */
  a: number;
  /** Form factor m (default 2 — UNVERIFIED: the commonly quoted Heywood 1988 fit value). */
  m: number;
}

/** Cycle-model options (all optional in the constructors; defaults in DEFAULT_CYCLE_MODEL_OPTIONS). */
export interface CycleModelOptions {
  /** Wall heat transfer (Woschni) on/off. */
  heatTransfer: boolean;
  /** Burn-rate model, see {@link CombustionModel}. */
  combustionModel: CombustionModel;
  /** Wiebe parameters for combustionModel 'wiebe'. */
  wiebe: WiebeOptions;
  /** End-gas autoignition (Livengood–Wu + end-gas burn-up + acoustic modes) on/off. */
  knock: boolean;
  /** End-gas ignition-delay model, or any object implementing IgnitionDelayModel. */
  ignitionDelayModel: IgnitionDelayModelId | IgnitionDelayModel;
  /**
   * Multiplier C_T on the entrainment speed, u_T = C_T·u′ (Keck 1982 eq. 4.1B with u_T = u′
   * when 1) — the prime calibration knob. Default: CFR_CALIBRATION (5.25; 1 = uncalibrated).
   */
  burnRateMultiplier: number;
  /** Multiplier on the burn-up length (Taylor microscale λ, τ_b = λ/S_L). Default: CFR_CALIBRATION (2; 1 = uncalibrated). */
  taylorScaleMultiplier: number;
  /**
   * Spark-kernel hand-off radius as a multiple C of the integral scale, r_ho = max(1 mm, C·l_I)
   * (ignition/kernel.ts handoffIntegralScaleMultiple; Fluent C = 1, Forte C = 2). Default: CFR_CALIBRATION.
   * options.ignition.kernel.handoffIntegralScaleMultiple, when given, overrides it.
   */
  kernelHandoffMultiple: number;
  /**
   * Multiplier on the (Bechtold–Matalon asymptotic) unburned-gas Markstein length passed to the spark
   * kernel's stretch factor. Default: CFR_CALIBRATION (1 = the asymptotic theory).
   */
  marksteinMultiplier: number;
  /** Turbulent-flame closure (u_T, ℓ_T), see {@link TurbulentFlameClosure}. */
  turbulentFlameClosure: TurbulentFlameClosure;
  /** End-gas autoignition integral, see {@link KnockIntegralMode}. */
  knockIntegral: KnockIntegralMode;
  /**
   * Carburettor-venturi discharge coefficient (throat C_D; the CFR has no throttle plate). The
   * ISO 5167-4 classical-venturi value 0.984 applies at Re ≥ 2e5; the CFR throat runs at Re ≈ 1e4
   * with the fuel-nozzle bridge in the throat.
   */
  venturiDischargeCoefficient: number;
  /**
   * Effective flow area C_D·A (m²) of a fixed, venturi-independent restriction in series with the
   * carburettor venturi (air horn, MON mixture-heater housing, runner bends), combined with the
   * venturi as 1/(C_D A)²_eff = 1/(C_D A)²_venturi + 1/(C_D A)²_restriction (series orifices,
   * incompressible limit). 0 = none.
   */
  intakeRestrictionArea: number;
  /**
   * Intake-port heat-transfer multiplier on the Dittus–Boelter pipe correlation over the heated
   * intake port (spec.manifolds.intakePortDiameter × intakePortLength, wall at
   * spec.walls.intakePortTemperature). 0 = adiabatic intake. Choi et al. 2018 (GT-Power TPA of
   * the ANL CFR) needed an intake-port multiplier of 4 (DoE range 2–6) to match the measured
   * volumetric efficiency.
   */
  intakePortHeatTransferMultiplier: number;
  /** Multiplier on the turbulence dissipation-length fraction f_L (turbulence.DEFAULT_LENGTH_SCALE_FRACTION). Default 1. */
  turbulenceLengthScaleFactor: number;
  /** K–k production multiplier c_β (1 = standard C_μ = 0.09 eddy viscosity). Default 1. */
  turbulenceProduction: number;
  /**
   * Fraction of the ideal shrouded-valve jet angular momentum that becomes swirl (the ideal jet
   * gives swirl ratio ≈ 23 at IVC for the CFR — too strong, DESIGN.md). Swirl only feeds the
   * diagnostic swirl state (turbulence.ts swirlFrictionToTurbulence = 0), so it is off by default.
   */
  swirlMomentumEfficiency: number;
  /** Multiplier on the Woschni coefficient (Choi et al. 2018 needed 1–3 for the CFR). Default: CFR_CALIBRATION (1.3). */
  woschniMultiplier: number;
  /**
   * In-cylinder heat-transfer correlation: 'woschni' (1967, with the combustion-velocity term; default)
   * or 'hohenberg' (1979, heat-transfer/correlations.ts). woschniMultiplier scales either.
   */
  heatTransferCorrelation: 'woschni' | 'hohenberg';
  /**
   * Multiplier on Woschni's combustion-induced gas-velocity constant C₂ (3.24e-3 m/(s K)); the overall
   * coefficient is scaled by woschniMultiplier. Default: CFR_CALIBRATION.
   */
  woschniCombustionTermMultiplier: number;
  /** Woschni variant (heat-transfer/correlations.ts). Default the original 1967 form. */
  woschniVariant: WoschniVariant;
  /** Multiplier on both valves' effective flow areas (C_D uncertainty). Default 1. */
  dischargeCoefficientMultiplier: number;
  /** Knock mode (1,0) amplitude decay time, s (default chemistry/knock KNOCK_DECAY_TIME, UNVERIFIED 1 ms). */
  knockDecayTime: number;
  /** End-gas temperature stratification ΔT for the autoignition burn time, K (CFR_CALIBRATION, 60 K). */
  knockStratificationDT: number;
  /**
   * Sequential-autoignition acoustic source: the end gas is split into this many equal-area shells
   * released one after the other, farthest from the flame first (KnockOscillator.setSequentialEndGasRegion);
   * 0 = the round-1 uniform release over the whole end-gas region. Default 16 (numerical resolution of
   * the sweep; results change < 1 % from 16 to 32 shells — fixer round 2 check).
   */
  knockSourceShells: number;
  /**
   * At knock onset ALL unburned gas autoignites — the end gas ahead of the front and the entrained-
   * but-unburned pockets of the flame brush, which share the unburned zone's T_u, p history (default
   * true, fixer round 2) — or only the end gas ahead of the front (false: round-1 model).
   */
  knockBrushAutoignition: boolean;
  /** Autoignition excitation time τ_e, s (knock.ts EXCITATION_TIME_PRF). */
  knockExcitationTime: number;
  /**
   * Measurement band of CycleSummary.mapo, Hz: [f_lo, f_hi] — an ideal band-pass over the acoustic
   * modes (a mode counts when its frequency at the knock-onset sound speed lies in the band;
   * KnockOscillator.bandGain), or null for the unfiltered modal sum. The (4,0)/(1,1) modes sit at
   * ≈ 18.5–19 kHz, just above 18 kHz: a real filter of unknown order (plus the 0.1°-sampling Nyquist
   * limit of the measurements) would pass part of them — ≈ ±15 % MAPO uncertainty (validation round 2). Default [4000, 18000]: Hoth & Kolodziej 2021, SON 2023 and the ANL
   * critical-CR data (test/fixtures/cfr_validation_README.md); Rockstroh 2018 used 6–20 kHz. The
   * knock-induced bulk pressure rise is not band-passed into MAPO (≈ 5–9 % of MAPO in band for the
   * round-1 burn-up, offline scipy analysis, validation round 2).
   */
  mapoBand: readonly [number, number] | null;
  /** Knock-pickup position (x, z) on the head face, cylinder frame, m (default CFR D-1 pickup). */
  knockSensor: readonly [number, number];
  /** Base integration step, crank degrees (default 0.25). */
  maxStepDeg: number;
  /** Step limit during dwell/spark/kernel/early flame, crank degrees (default 0.05). */
  fineStepDeg: number;
  /** Step limit during the end-gas burn-up after autoignition, s (default 2 µs, knock.ts note). */
  knockBurnStep: number;
  /**
   * Silent cycles run after construction/reset (default 3). After each, both plenums are
   * extrapolated to their periodic state (CycleModel.relaxPlenum), so the first emitted cycle is
   * converged (RON: IMEP within 2e-5 of cycle 30).
   */
  warmupCycles: number;
  /** Valve running clearance for the lash-corrected lift profiles, m (default CFR 0.008 in). */
  valveLash: number;
  /** Crankcase gauge pressure (relative to ambient) for the gas torque, Pa (default CFR mid-band). */
  crankcaseGaugePressure: number;
  /** Friction-model inputs (Patton–Nitschke–Heywood); default CFR_FRICTION. null = no friction. */
  friction: PnhFrictionInputs | null;
  /** Burned-mass fraction at which the remaining unburned gas is merged into the burned zone (numerical, 1e-5). */
  burnoutFraction: number;
  /**
   * Numerical guard: a zone's wall heat loss is limited so that it cannot cool the zone faster
   * than this time constant, s (only acts on vanishing zones; default 0.1 ms).
   */
  zoneHeatLossMinTime: number;
  /** Options passed to the IgnitionSystem (coil / gap / kernel). */
  ignition: IgnitionSystemOptions;
  /** Collect per-sub-model CPU times (performance.now) — diagnostic, default false. */
  profile: boolean;
  /**
   * Burned-zone NO in the thermodynamics: 'kinetic' = the burned zone is the TP equilibrium with its
   * NO replaced by the rate-controlled (extended-Zeldovich) amount, energy and volume included
   * (closure.ts noCoupled); 'equilibrium' = full equilibrium, the kinetic NO a passive tracer written
   * into the gas at EVO (round-1 model; the fuel–air-cycle and two-zone Cantera oracles use it);
   * 'auto' (default) = 'kinetic' for the 'entrainment' combustion model, else 'equilibrium'.
   */
  burnedNOThermo: 'kinetic' | 'equilibrium' | 'auto';
  /**
   * Wall temperatures: 'lumped' (default when the spec provides walls.thermalResistance) = coolant +
   * R_i·Q̄_i per surface from the cycle-mean heat flows (CycleModel.updateWalls), 'fixed' = the spec
   * values shifted with the coolant temperature (round-1 model).
   */
  wallTemperatureModel: 'lumped' | 'fixed';
  /**
   * Thermal time constant of the lumped walls after the warm-up, s. 0 (default) = quasi-steady (each
   * cycle uses the walls of the previous cycle's heat load). The real CFR head/liner take minutes;
   * UNVERIFIED — set it for transient realism, keep 0 for steady-state validation.
   */
  wallTimeConstant: number;
  /**
   * Crevice-flow zone (Namazian & Heywood 1982, SAE 820088: gas in the narrow piston/ring crevices at
   * the wall temperature and the cylinder pressure, filled by the unburned or burned gas at the
   * crevice entrance and emptied back as the pressure falls; cycle-model.ts). Needs heat transfer and
   * spec.geometry.quenchCreviceVolume > 0. Default true.
   */
  creviceModel: boolean;
}

/**
 * Engine-independent defaults: sub-model switches, numerical controls and the UNIVERSAL closures
 * (default knock-delay model CFR_KNOCK_DELAY_MODEL, Markstein multiplier = the Bradley 1998
 * measured/theory ratio, kernel hand-off multiple = Forte's C_m1 — calibration.ts: not fitted to an
 * engine); the ENGINE-FITTED closures at their uncalibrated literature values (multipliers 1, venturi
 * C_D 0.984 = ISO 5167-4, adiabatic intake port, end-gas stratification 15 K = knock.ts
 * END_GAS_STRATIFICATION_DT); no engine hardware (no friction, zero valve lash, crankcase at ambient,
 * unfiltered MAPO). The knock pickup is a geometric placeholder replaced per spec
 * ({@link engineCycleOptionDefaults}). Pass these values to get the uncalibrated model of any engine.
 */
export const NEUTRAL_CYCLE_MODEL_OPTIONS: Readonly<CycleModelOptions> = Object.freeze({
  heatTransfer: true,
  combustionModel: 'entrainment' as CombustionModel,
  wiebe: Object.freeze({ durationDeg: 50, a: -Math.log(0.001), m: 2 }) as WiebeOptions,
  knock: true,
  ignitionDelayModel: CFR_KNOCK_DELAY_MODEL as IgnitionDelayModelId,
  burnRateMultiplier: 1,
  taylorScaleMultiplier: 1,
  kernelHandoffMultiple: CFR_CALIBRATION.kernelHandoffMultiple.value,
  marksteinMultiplier: CFR_CALIBRATION.marksteinMultiplier.value,
  turbulentFlameClosure: 'kk-taylor' as TurbulentFlameClosure,
  knockIntegral: 'single' as KnockIntegralMode,
  venturiDischargeCoefficient: VENTURI_DISCHARGE_COEFFICIENT,
  intakeRestrictionArea: 0,
  intakePortHeatTransferMultiplier: 0,
  turbulenceLengthScaleFactor: 1,
  turbulenceProduction: 1,
  swirlMomentumEfficiency: 0,
  woschniMultiplier: 1,
  heatTransferCorrelation: 'woschni' as const,
  woschniCombustionTermMultiplier: 1,
  woschniVariant: 'woschni1967' as WoschniVariant,
  dischargeCoefficientMultiplier: 1,
  knockDecayTime: KNOCK_DECAY_TIME,
  knockStratificationDT: END_GAS_STRATIFICATION_DT,
  knockExcitationTime: EXCITATION_TIME_PRF,
  knockSourceShells: 16,
  knockBrushAutoignition: true,
  knockSensor: Object.freeze([0, 0]) as readonly [number, number],
  mapoBand: null,
  maxStepDeg: 0.25,
  fineStepDeg: 0.05,
  knockBurnStep: 2e-6,
  warmupCycles: 3,
  valveLash: 0,
  crankcaseGaugePressure: 0,
  friction: null,
  burnoutFraction: 1e-5,
  zoneHeatLossMinTime: 1e-4,
  ignition: Object.freeze({}) as IgnitionSystemOptions,
  profile: false,
  burnedNOThermo: 'auto' as const,
  wallTemperatureModel: 'lumped' as const,
  wallTimeConstant: 0,
  creviceModel: true,
});

/** The calibrated closures of a calibration set as cycle-model options. */
export function calibrationOptions(cal: CalibrationSet): Partial<CycleModelOptions> {
  return {
    burnRateMultiplier: cal.burnRateMultiplier.value,
    taylorScaleMultiplier: cal.taylorScaleMultiplier.value,
    kernelHandoffMultiple: cal.kernelHandoffMultiple.value,
    marksteinMultiplier: cal.marksteinMultiplier.value,
    venturiDischargeCoefficient: cal.venturiDischargeCoefficient.value,
    intakePortHeatTransferMultiplier: cal.intakePortHeatTransferMultiplier.value,
    woschniMultiplier: cal.woschniMultiplier.value,
    knockStratificationDT: cal.knockStratificationDT.value,
  };
}

/**
 * Model T friction inputs of the cycle model: the PNH inputs of the engine (engines/model-t.ts
 * MODEL_T_FRICTION — geometry, oil, valvetrain; pure PNH) with the two factors of MODEL_T_CALIBRATION that
 * account for what PNH (fitted on 1980s engines) lacks: the 1/4 in cast-iron ring pack (ringTensionFactor)
 * and Ford's transmission-output measurement basis — fan, generator, magneto and the high-gear churning of
 * the flywheel/planetary in the oil pit (auxiliaryFactor). WOT FMEP ≈ 1.0–1.2 bar at 500–1800 rpm
 * (η_m 0.78 at 1500 rpm; ALAM/SAE assumption 0.75).
 */
export const MODEL_T_CALIBRATED_FRICTION: Readonly<PnhFrictionInputs> = Object.freeze({
  ...MODEL_T_FRICTION,
  ringTensionFactor: MODEL_T_CALIBRATION.ringTensionFactor.value,
  auxiliaryFactor: MODEL_T_CALIBRATION.auxiliaryFactor.value,
});

/**
 * Per-engine option defaults, keyed by EngineSpec.id (engines/index.ts registry ids): the engine's
 * calibration set plus its hardware — friction inputs, valve lash, crankcase pressure, knock pickup
 * and MAPO band.
 */
export const ENGINE_CYCLE_OPTION_DEFAULTS: Readonly<Record<string, Readonly<Partial<CycleModelOptions>>>> = Object.freeze({
  'cfr-f1': Object.freeze({
    ...calibrationOptions(CFR_CALIBRATION),
    ignitionDelayModel: CFR_KNOCK_DELAY_MODEL as IgnitionDelayModelId,
    // CFR D-1 detonation pickup in the head face (cfr.ts CFR_KNOCK_PICKUP); measurement band of the ANL data
    knockSensor: Object.freeze([CFR_KNOCK_PICKUP.position[0], CFR_KNOCK_PICKUP.position[2]]) as readonly [number, number],
    mapoBand: Object.freeze([4000, 18000]) as readonly [number, number],
    valveLash: CFR_VALVE_LASH,
    crankcaseGaugePressure: CFR_CRANKCASE_GAUGE_PRESSURE,
    friction: CFR_FRICTION,
  }),
  'ford-model-t': Object.freeze({
    ...calibrationOptions(MODEL_T_CALIBRATION),
    // (calibration phase) the Model T-only closure of MODEL_T_CALIBRATION: the turbulence dissipation length
    turbulenceLengthScaleFactor: MODEL_T_CALIBRATION.turbulenceLengthScaleFactor.value,
    ignitionDelayModel: CFR_KNOCK_DELAY_MODEL as IgnitionDelayModelId,
    // UNVERIFIED (no knock instrumentation exists for the Model T): a virtual plug-mounted transducer at the
    // spark plug's planform position over the valve pocket (chemistry/knock.ts virtualKnockSensor — the usual
    // retrofit for in-cylinder pressure on old heads); the L-head's own depth-averaged planform modes
    // (createKnockOscillator) are defined over the whole chamber, pocket included
    knockSensor: virtualKnockSensor(MODEL_T),
    // L_HEAD_MAPO_BAND [2, 18] kHz (UNVERIFIED convention): the CFR's 4–18 kHz ANL band would cut the
    // L-head's 3.4–4.0 kHz bore-to-pocket fundamental (knock.ts: 3812 Hz at c = 950 m/s)
    mapoBand: L_HEAD_MAPO_BAND,
    valveLash: MODEL_T_VALVE_LASH,
    // UNVERIFIED: the crankcase breathes to the atmosphere through the oil-filler breather (no PCV, no pump)
    crankcaseGaugePressure: 0,
    friction: MODEL_T_CALIBRATED_FRICTION,
  }),
});

/**
 * Defaults: the uncalibrated literature sub-models with the calibrated closures and hardware of the
 * CFR F-1 (calibration.ts CFR_CALIBRATION — values, literature ranges and evidence), i.e.
 * NEUTRAL_CYCLE_MODEL_OPTIONS ← ENGINE_CYCLE_OPTION_DEFAULTS['cfr-f1'].
 */
export const DEFAULT_CYCLE_MODEL_OPTIONS: Readonly<CycleModelOptions> = Object.freeze({
  ...NEUTRAL_CYCLE_MODEL_OPTIONS,
  ...ENGINE_CYCLE_OPTION_DEFAULTS['cfr-f1'],
} as CycleModelOptions);

const engineDefaultsCache = new Map<string, Readonly<CycleModelOptions>>();

/**
 * Resolved default options of an engine: the CFR set for 'cfr-f1' and for specs without an id (ad-hoc
 * CFR derivatives, the pre-registry behaviour); NEUTRAL ← ENGINE_CYCLE_OPTION_DEFAULTS[id] for a
 * registered id; the neutral set with a generic head-face knock pickup near the liner for any other id.
 */
export function engineCycleOptionDefaults(spec: EngineSpec): Readonly<CycleModelOptions> {
  const id = spec.id;
  if (id === undefined || id === 'cfr-f1') return DEFAULT_CYCLE_MODEL_OPTIONS;
  const eng = ENGINE_CYCLE_OPTION_DEFAULTS[id];
  if (eng) {
    let r = engineDefaultsCache.get(id);
    if (!r) {
      r = Object.freeze({ ...NEUTRAL_CYCLE_MODEL_OPTIONS, ...eng } as CycleModelOptions);
      engineDefaultsCache.set(id, r);
    }
    return r;
  }
  // unregistered engine: neutral, with a generic pickup in the head face 3 mm inside the liner (like the CFR D-1)
  const knockSensor = Object.freeze([0, -(0.5 * spec.geometry.bore - 0.003)]) as readonly [number, number];
  return Object.freeze({ ...NEUTRAL_CYCLE_MODEL_OPTIONS, knockSensor } as CycleModelOptions);
}

/**
 * Merge user options over the defaults (shallow; `wiebe` merged one level). With a spec the defaults are
 * that engine's ({@link engineCycleOptionDefaults}); without one, the CFR set (DEFAULT_CYCLE_MODEL_OPTIONS).
 */
export function resolveCycleOptions(o: Partial<CycleModelOptions> = {}, spec?: EngineSpec): CycleModelOptions {
  const d = spec ? engineCycleOptionDefaults(spec) : DEFAULT_CYCLE_MODEL_OPTIONS;
  const r: CycleModelOptions = { ...d, ...stripUndefined(o) } as CycleModelOptions;
  r.wiebe = { ...d.wiebe, ...(o.wiebe ?? {}) };
  r.maxStepDeg = clamp(r.maxStepDeg, 1e-3, 2);
  r.fineStepDeg = clamp(Math.min(r.fineStepDeg, r.maxStepDeg), 1e-4, 2);
  r.knockBurnStep = clamp(r.knockBurnStep, 1e-8, 1e-4);
  r.warmupCycles = Math.max(0, Math.floor(r.warmupCycles));
  return r;
}

function stripUndefined<T extends object>(o: T): Partial<T> {
  const out: Partial<T> = {};
  for (const k of Object.keys(o) as (keyof T)[]) if (o[k] !== undefined) out[k] = o[k];
  return out;
}

const clamp = (x: number, lo: number, hi: number): number => (x < lo ? lo : x > hi ? hi : x);

/** Numeric OperatingPoint fields repaired by sanitizeOperatingPoint when non-finite or missing. */
const NUMERIC_FIELDS = [
  'rpm', 'loadTorque', 'throttle', 'ambientPressure', 'ambientTemperature', 'relativeHumidity', 'intakeMixtureTemperature',
  'equivalenceRatio', 'sparkAdvanceDeg', 'dwellTime', 'compressionRatio', 'egrFraction', 'coolantTemperature',
] as const satisfies readonly (keyof OperatingPoint)[];
const CFR_RON_CONDITIONS_ON = 90;

/**
 * Largest water-vapour mole fraction of the intake air accepted by the sanitiser (x_H2O =
 * RH·f·p_sat(T_amb)/p_amb ≤ this). A validity limit, not a physical constant: 0.2 covers any
 * weather (saturated air at 60 °C and 1 atm has x ≈ 0.2); above it the "air" is mostly steam and
 * RH·p_sat can even exceed p (round 1: RH 1 at 400 K was accepted).
 */
export const MAX_INTAKE_WATER_FRACTION = 0.2;

/**
 * Sanitise an operating point into the model's valid range (a copy). φ is limited to [0, 2.5]:
 * φ = 0 is air only (the drained-carburettor ASTM compression-pressure check); above 2.5 the
 * burned-gas equilibrium (no condensed carbon) approaches infeasibility near φ ≈ 3.1 (DESIGN.md
 * integration notes). Relative humidity is limited so that x_H2O ≤ MAX_INTAKE_WATER_FRACTION.
 * Non-finite or missing numeric fields (and an invalid speed mode / fuel / load / ignition source) are
 * replaced by the corresponding field of `fallback` (default: the engine's default operating point,
 * {@link defaultOperatingPointOf}) before clamping. The optional `load` (copied, validated) and
 * `ignitionSource` are kept.
 */
export function sanitizeOperatingPoint(spec: EngineSpec, op: OperatingPoint, fallback?: OperatingPoint): OperatingPoint {
  const [crMin, crMax] = spec.geometry.compressionRatioRange;
  const o: OperatingPoint = { ...op, fuel: { ...op.fuel } };
  // Non-finite / missing numeric fields (NaN, ±Infinity, undefined) are repaired BEFORE clamping —
  // clamp() passes NaN through, and NaN rpm / φ threw inside the model (the worker then stopped) while
  // an undefined spark advance silently disabled the spark (validation round 2). The repair value is
  // the caller's previous operating point, else the engine's default operating point (CFR: the ASTM
  // D2699 Research-method preset).
  const fb = fallback ?? defaultOperatingPointOf(spec);
  if (o.load !== undefined) {
    const l = sanitizeLoad(spec, o.load, fb.load);
    if (l) o.load = l;
    else delete o.load;
  }
  if (o.ignitionSource !== undefined && o.ignitionSource !== 'magneto' && o.ignitionSource !== 'battery') {
    if (fb.ignitionSource === 'magneto' || fb.ignitionSource === 'battery') o.ignitionSource = fb.ignitionSource;
    else delete o.ignitionSource;
  }
  for (const k of NUMERIC_FIELDS) {
    const v = o[k] as unknown;
    if (typeof v !== 'number' || !Number.isFinite(v)) (o[k] as number) = fb[k] as number;
  }
  if (o.speedMode !== 'fixed' && o.speedMode !== 'free') o.speedMode = fb.speedMode;
  if (!o.fuel || (o.fuel.kind !== 'PRF' && o.fuel.kind !== 'pure')) o.fuel = { ...fb.fuel };
  if (o.fuel.kind === 'PRF' && !Number.isFinite(o.fuel.octaneNumber)) {
    o.fuel = { kind: 'PRF', octaneNumber: fb.fuel.kind === 'PRF' ? fb.fuel.octaneNumber : CFR_RON_CONDITIONS_ON };
  }
  o.rpm = clamp(o.rpm, 60, 6000);
  o.throttle = clamp(o.throttle, 0, 1);
  o.equivalenceRatio = clamp(o.equivalenceRatio, 0, 2.5);
  o.compressionRatio = clamp(o.compressionRatio, Math.max(1.5, crMin), crMax);
  o.sparkAdvanceDeg = clamp(o.sparkAdvanceDeg, -60, 120);
  o.dwellTime = clamp(o.dwellTime, 1e-4, 20e-3);
  o.egrFraction = clamp(o.egrFraction, 0, 0.6);
  o.ambientPressure = clamp(o.ambientPressure, 2e4, 3e5);
  o.ambientTemperature = clamp(o.ambientTemperature, 200, 400);
  o.intakeMixtureTemperature = clamp(o.intakeMixtureTemperature, 200, 600);
  o.relativeHumidity = clamp(o.relativeHumidity, 0, 1);
  const pvMax = (MAX_INTAKE_WATER_FRACTION * o.ambientPressure) / (moistAirEnhancementFactor(o.ambientTemperature, o.ambientPressure) * waterSaturationPressure(o.ambientTemperature));
  if (o.relativeHumidity > pvMax) o.relativeHumidity = pvMax;
  o.coolantTemperature = clamp(o.coolantTemperature, 250, 450);
  if (o.fuel.kind === 'PRF') o.fuel = { kind: 'PRF', octaneNumber: clamp(o.fuel.octaneNumber, 0, 100) };
  return o;
}

/**
 * Default operating point of an engine: the registry's (engines/index.ts) for a registered spec.id,
 * else the CFR Research-method preset (the pre-registry fallback).
 */
export function defaultOperatingPointOf(spec: EngineSpec): OperatingPoint {
  return engineOfSpec(spec)?.defaultOperatingPoint ?? CFR_RON_CONDITIONS;
}

/** Largest |road grade| accepted for the 'vehicle' load (rise/run). A validity limit: 0.5 ≈ 27°. */
const MAX_GRADE = 0.5;

/**
 * A validated copy of a load spec, or the fallback's (copied) when it is malformed, or undefined
 * (= constant loadTorque) when neither is valid. 'vehicle' needs spec.vehicle and a gear it lists (or
 * 'neutral'); 'brake' needs refRpm > 0 and a finite exponent (clamped to [0, 4]).
 */
function sanitizeLoad(spec: EngineSpec, l: LoadSpec, fb: LoadSpec | undefined): LoadSpec | undefined {
  const ok = (x: LoadSpec | undefined): LoadSpec | undefined => {
    if (!x || typeof x !== 'object') return undefined;
    if (x.kind === 'constant') return { kind: 'constant' };
    if (x.kind === 'brake') {
      if (!(Number.isFinite(x.refRpm) && x.refRpm > 0 && Number.isFinite(x.exponent))) return undefined;
      return { kind: 'brake', refRpm: x.refRpm, exponent: clamp(x.exponent, 0, 4) };
    }
    if (x.kind === 'vehicle') {
      const v = spec.vehicle;
      if (!v || typeof x.gear !== 'string' || !(x.gear === 'neutral' || x.gear in v.gears)) return undefined;
      return { kind: 'vehicle', gear: x.gear, grade: Number.isFinite(x.grade) ? clamp(x.grade, -MAX_GRADE, MAX_GRADE) : 0 };
    }
    return undefined;
  };
  return ok(l) ?? ok(fb);
}
