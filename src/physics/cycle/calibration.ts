/**
 * The ONE global calibration of the cycle model for the Waukesha CFR F-1 (DESIGN.md calibration
 * policy): every genuinely uncertain closure that is tuned, with its value, the literature range it
 * must stay in, the source of that range, and the validation evidence that set it. The same set
 * serves every operating point (RON, MON, all octane numbers, spark timings, fuels); nothing is
 * tuned per operating point. DEFAULT_CYCLE_MODEL_OPTIONS (options.ts) takes its defaults from here.
 * Multi-engine (Model T integration): one set PER ENGINE, selected by EngineSpec.id through the
 * engine-keyed option defaults of options.ts — CFR_CALIBRATION ('cfr-f1') and MODEL_T_CALIBRATION
 * ('ford-model-t', calibrated 2026-10 against Ford's WOT table, Upton's MBT and suction data; see below).
 *
 * Validation round 2 (fixer pass, 2026-09-30). The model changed under the calibration: crevice zone,
 * lumped wall temperatures, rate-controlled NO in the burned-zone energy, venturi metering AIR at the
 * intake-air temperature, integral-scale input of the spark kernel, soft kernel extinction, Douaud–
 * Eyzat with LLNL relative φ/residual/low-ON sensitivities, band-limited MAPO (4–18 kHz), all unburned
 * gas autoigniting at the knock onset with a sequential (front-like) acoustic source, 0.25 L intake
 * plenum (cfr.ts). The burn / heat-transfer set was re-chosen by Nelder–Mead on the CFR data metrics
 * (tmp/fix2/calib.ts, 60 + 45 evaluations; targets and weights: Choi 2018 Fig. 9 spark +5 apparent
 * CA10/50/90 ±1/1/1.5°, spark −13 (knock off) CA10/50/90 ±1/1/1.5°, KW17 and Choi 2018 Fig. 2 PRF98
 * CA50 ±1–1.5°, gIMEP(+5) ±2 %, gross work per unit fuel (Choi Table 6, KW17) ±2 %, gIMEP(−13)/gIMEP(+5)
 * +3.1 ± 1.5 %, apparent-heat-release peak ±3 % and post-peak decline ×(1 ± 0.3), trace RMS ±2 % of
 * peak). The cost valley is flat (≈ 22–25 over C_T 5.0–5.7, C_λ 2.2–3.3, C_m1 1.3–2.0, Woschni
 * 1.18–1.24); round values near its bottom were taken. ΔT was then set on the MAPO of the two ANL
 * standard-knock states and the knock-limited CR scale; it stops at its upper bound (60 K). "Apparent" burn angles = single-zone net heat release, γ = 1.30, the same routine on measured
 * and simulated pressure (test/validation/measured_cfr_data.test.ts).
 *
 * Before (round-1 set on the round-1 model) → after (round-2 set on the round-2 model):
 *   Choi PRF100 CR 7.55 λ 0.89, spark +5: apparent CA10/50/90 +1.6/+0.9/+1.6° → +1.3/+1.4/+3.4°, gIMEP
 *     +1.7 → +2.0 %, apparent-HR peak +5.9 → +5.0 %, post-peak decline ×2.09 → ×1.33;
 *   spark −13 (knock off): CA10/50 +0.4/−0.7° → −0.5/−1.3°, apparent-HR peak +4.6 → 0.0 %, decline ×1.45
 *     → ×1.00, gIMEP(−13)/gIMEP(+5) +1.5 → +1.7 % (measured +3.1 %);
 *   standard knock (spark −12.72°): onset 12.4° → 12.2° (knock point 11.0 ± 1.0°), MAPO 1.50 bar (all
 *     modes) → 1.16 bar (4–18 kHz; Hoth 1.43 bar), max dp/dθ 3.36 → 2.92 bar/° (average trace 3.88),
 *     peak 47.9 → 47.3 bar (47.9);
 *   KW17 PRF98 CA50 7.0 → 6.6° (8.8), knock onset 12.5 → 12.3° (12.0), MAPO 0.71 → 0.95 bar (0.63–0.69);
 *     fuel flow Choi/KW +1.6/−0.5 → +2.5/+0.6 %, work per fuel −0.3/0.0 → −1.2/−0.8 %; intake-port RMS
 *     0.025 → 0.014 bar;
 *   ASTM compression pressure (motored walls now) RON 930 +4.5 → +1.9 %, MON 930 +5.4 → +4.9 %, MON/RON
 *     ratio +0.9 → +2.9 %, venturi-size effect 9.0 → 7.2 % (8.0 %);
 *   knock-limited CR (MAPO 0.67 bar) − guide: RON 60/70/80/90/100 −0.84/−0.38/−0.02/+0.17/−0.17 →
 *     −0.94/−0.77/−0.45/−0.14/−0.37; MON 70/80/90/100 +0.36/+0.61/+0.55/+0.18 → −0.16/−0.22/−0.14/−0.38
 *     (mean |error| ON 70–100 0.29 → 0.27 CR); PRF 90 rated RON/MON 92.5/95.2 → 87.6/88.4 (sensitivity
 *     −2.7 → −0.9), PRF 100 98.9/101.5 → 97.6/97.0; RON PRF 90 CR_KL at spark 9/13/17°: 6.73/6.60/6.60
 *     (U-shaped) → 6.59/6.29/6.13 (−0.075/−0.041 CR/°); Rockstroh PRF90 λ 1: CR at 0.67 bar (33 °C) −0.24
 *     → −0.42, 33 → 150 °C −0.54 → −0.02 CR (measured −0.1…+0.12); MAPO(λ 1)/MAPO(λ 0.89) at fixed CR
 *     1.29 → 0.82 (measured 0.42).
 */

/** One calibrated parameter. */
export interface CalibratedParameter {
  /** Value used by default. */
  readonly value: number;
  /** Literature-supported range [lo, hi] the value must stay within. */
  readonly range: readonly [number, number];
  /** Where the range comes from. */
  readonly source: string;
  /** The validation evidence that set the value. */
  readonly evidence: string;
}

const BURN_EVIDENCE =
  'Nelder–Mead on the CFR burn metrics (file header). Choi et al. 2018 (SAE 2018-01-0848) Fig. 9 PRF100 CR 7.55 λ 0.89: ' +
  'spark +5 apparent CA10/50/90 +1.3/+1.4/+3.4°, spark −13 (knock off) −0.5/−1.3/+0.2°; KW17 PRF98 CA50 −2.2°, Choi 2018 ' +
  'Fig. 2 PRF98 CA50 −2.5°. Model form: no single set fits both spark timings (the +5 case burns ≈ 2.7° too slowly ' +
  'relative to −13) nor PRF98 at CR ≈ 7.2 and PRF100 at 7.55 (≈ 1.2° apart) — C_λ and C_m1 sit near their upper bounds.';

export const CFR_CALIBRATION = Object.freeze({
  /**
   * C_T: entrainment speed u_T = C_T·u′ (K–k u′), Keck 1982 eq. 4.1B.
   */
  burnRateMultiplier: Object.freeze({
    value: 5.3,
    range: [2.6, 13] as const,
    source:
      "Keck 1982 (19th Symp. Combust., p. 1451) eq. 4.10, u_T = 0.08 ū_i (ρ_u/ρ_i)^½ (±10 %, his Fig. 14), " +
      'evaluated on this model’s own burn (ū_i from ε_v, A_p and the inlet-valve area at maximum lift): ' +
      'u_T/u′ = 6.4–13.5 with the shrouded valve’s flow area, 2.9–6.1 with the unshrouded curtain area ' +
      '(scratch evaluation, validation round 1); DESIGN.md quoted ≈ 2.6 at CFR conditions.',
    evidence: BURN_EVIDENCE,
  }) as CalibratedParameter,
  /** C_λ: burn-up length ℓ_T = C_λ·λ (Taylor microscale), τ_b = ℓ_T/S_L. */
  taylorScaleMultiplier: Object.freeze({
    value: 3.0,
    range: [1.3, 3.3] as const,
    source:
      'Keck 1982 Fig. 15, ℓ_T = 0.8 L_IV (ρ_i/ρ_u)^¾ (±25 %), evaluated on the model’s burn: ' +
      'ℓ_T/λ = 1.7–2.6 (range × 0.75…1.25).',
    evidence: BURN_EVIDENCE + ' Round 1: 2.0 (with the round-1 kernel hand-off at 2.2× the integral scale).',
  }) as CalibratedParameter,
  /** Spark-kernel hand-off radius as a multiple of the integral scale, r_ho = max(1 mm, C·l_I). */
  kernelHandoffMultiple: Object.freeze({
    value: 2.0,
    range: [1, 2] as const,
    source:
      'ANSYS Fluent Theory Guide 2025 R1 §10.1.2 (Herweg–Maly spark model: kernel grown "until the length scale ' +
      'is reached", C = 1) and ANSYS Forte / [HM-R] SAE 18PFL-1037 eq. 13 (C_m1 = 2.0), both cited in ignition/kernel.ts.',
    evidence:
      BURN_EVIDENCE +
      ' With the kernel now given the true integral scale C_ε L (round 1 passed the dissipation length L ≈ 2 l_I), ' +
      'C = 2 hands off at ≈ 4.2 mm at RON (round 1: 5.3 mm; C = 1: 2.1 mm and the spark-timing split grows).',
  }) as CalibratedParameter,
  /** Multiplier on the asymptotic Markstein length used by the spark kernel (combustion/laminar-flame-speed.ts). */
  marksteinMultiplier: Object.freeze({
    value: 0.5,
    range: [0, 1] as const,
    source:
      'Bechtold–Matalon asymptotic theory (1 = as computed; Matalon 2011: measured Markstein numbers are usually ' +
      'smaller). Measured burned-gas Markstein length of stoichiometric iso-octane at 10 bar: Bradley et al. 1998 ' +
      '(358 K) ≈ 0.19 mm, Jerzembeck et al. 2009 (373 K) ≈ 0 (raster Fig. 4 of Jerzembeck 2009, as read in ' +
      'combustion/laminar-flame-speed.test.ts) vs 0.37 mm from the theory here → ratio 0.52 and ≈ 0.',
    evidence:
      'Set to the Bradley 1998 ratio (0.5), not tuned on engine data. With the theory value the Herweg–Maly kernel ' +
      'quenched (I₀ ≤ 0 at r ≈ 1.4–1.7 mm) at ordinary CFR points — φ 0.7/0.9/1.5, EGR 0.2, 1200–2400 rpm, spark 50°, ' +
      'throttle 0.1 — mostly as fire/misfire period-2 cycles (numerics / code review r2); with 0.5 all of them fire.',
  }) as CalibratedParameter,
  /** Woschni (1967) heat-transfer-coefficient multiplier (all phases). */
  woschniMultiplier: Object.freeze({
    value: 1.2,
    range: [1, 3] as const,
    source:
      'Choi et al. 2018 §Heat transfer model: the GT-Power TPA of the same ANL CFR needed a ' +
      'combustion-chamber convection multiplier of 1–3 (1.5 default in their DoE, Table 5).',
    evidence:
      BURN_EVIDENCE +
      ' Gross work per unit fuel −1.3/−0.9 % (Choi Table 6 / KW17), apparent-HR decline after the peak ×1.33 (spark +5) ' +
      'and ×1.00 (−13) of measured (round 1, no crevice zone: ×2.09 / ×1.45 at 1.3). Remaining: model MBT ≈ 5° BTDC vs ' +
      '9–10° measured (Hoth 2021b) — a lower C₂ (woschniCombustionTermMultiplier) or Hohenberg did not move it at ' +
      'matched work per fuel (fixer round 2 scans).',
  }) as CalibratedParameter,
  /** Intake-port heat-transfer multiplier on Dittus–Boelter over the heated port (0 = adiabatic). */
  intakePortHeatTransferMultiplier: Object.freeze({
    value: 6,
    range: [2, 6] as const,
    source:
      'Choi et al. 2018: intake-port heat-transfer multiplier 4 (DoE range 2–6, Table 5) to match the ' +
      'measured volumetric efficiency within ±1 %. Only multiplier × port length matters here and the ' +
      'port length (0.1 m, cfr.ts) is UNVERIFIED.',
    evidence:
      'Round 1 (unchanged): fuel flow Choi Table 6 PRF98 +2.5 % and KW17 +0.6 % (round 2), T_IVC 384 K ' +
      '(Choi/Pal TPA model: 406–410 K), trapped mass 668 mg (TPA 628 mg).',
  }) as CalibratedParameter,
  /** Carburettor-venturi (9/16 in throat) discharge coefficient. */
  venturiDischargeCoefficient: Object.freeze({
    value: 0.6,
    range: [0.6, 0.984] as const,
    source:
      'UNVERIFIED range: sharp-edged orifice ≈ 0.6 to the ISO 5167-4:2022 classical venturi 0.984 ' +
      '(Re ≥ 2e5); the CFR throat runs at Re ≈ 1e4 with the fuel-nozzle bridge across it.',
    evidence:
      'With the venturi now metering AIR at the intake-air temperature (code review r2): ASTM D2700 Fig. 2 (grade A) ' +
      '9/16 → 3/4 in venturi effect on the MON compression pressure 7.2 % at the bound (measured 8.0 %; 0.62: 6.5 %); ' +
      'MON/RON compression-pressure ratio +2.9 % (residual: a MON-only restriction downstream of the heater would be ' +
      'needed); fuel flow +2.5/+0.6 %.',
  }) as CalibratedParameter,
  /**
   * End-gas temperature stratification ΔT for the autoignition burn-up time (knock intensity).
   */
  knockStratificationDT: Object.freeze({
    value: 60,
    range: [10, 60] as const,
    source:
      '10–30 K: natural thermal stratification of the bulk charge (Sjöberg & Dec SAE 2005-01-0113, numbers not ' +
      'checked; Kokjohn, Musculus & Reitz, OSTI 1184578: bulk 825–845 K, ±5 K std, fetched); up to ≈ 60 K: end-gas ' +
      'temperature fluctuations "exceeding 20 K" measured by LIF in an SI end gas (Schießl & Maas, Combust. Flame ' +
      '133 (2003) 19, abstract) — the autoigniting gas now includes the flame-brush pockets and near-wall end gas ' +
      '(UNVERIFIED upper bound: ±20 K amplitude taken as a 40–60 K spread).',
    evidence:
      'At the upper bound. MAPO (4–18 kHz) of the two ANL standard-knock states, which differ by 2.2× in measured ' +
      'MAPO but only 1.2× in the model: Choi/Pal PRF100 CR 7.55 spark −12.72° 1.16 bar (Hoth & Kolodziej 2021b Fig. 1: ' +
      '1.43 bar) and KW17 PRF98 0.95 bar (0.63–0.69; with Rockstroh 0.67 and Hoth Table 6 0.69 bar the usual ' +
      'standard-knock equivalence, cfr_validation_README.md): geometric-mean ratio 1.08 (50 K: 1.39 / 1.14 bar, 1.29). ' +
      'Knock-limited CR (MAPO 0.67 bar) − guide RON 80/90/100 −0.45/−0.14/−0.37 (50 K: −0.55/−0.24/−0.48); ' +
      'the remaining offset needs ΔT beyond the range (model form: CFR_KNOCK_DELAY_MODEL note).',
  }) as CalibratedParameter,
});

/**
 * End-gas ignition-delay model of the default knock prediction (fixer round 2): Douaud & Eyzat (1978,
 * SAE 780080), the PRF correlation fitted to CFR knock, with the RELATIVE φ, residual and low-octane
 * (ON < 80) dependence of the LLNL-2011 detailed-chemistry table (chemistry/ignition-delay-llnl.ts
 * douaudEyzatLLNL; at φ 1.1, 6 % residual and ON ≥ 80 it is Douaud–Eyzat).
 * Evidence: Douaud–Eyzat's (ON/100)^3.402 → 0 made PRF 0–10 autoignite at IVC and PRF 60–70 knock far
 * too early; it has no φ term (model MAPO(λ 1)/MAPO(λ 0.89) 1.29 vs 0.42 measured). With the round-2
 * knock model (all unburned gas autoigniting, sequential source, ΔT 60 K): knock-limited CR (MAPO
 * 0.67 bar) − guide RON 60/70/80/90/100 −0.94/−0.77/−0.45/−0.14/−0.37, MON 70/80/90/100
 * −0.16/−0.22/−0.14/−0.38 (PRF 90 rated RON 87.6 / MON 88.4; PRF 100 97.6 / 97.0); spark 9 → 17° lowers
 * the RON PRF 90 knock-limited CR by 0.46; Rockstroh 33 → 150 °C −0.02 CR; λ 1/λ 0.89 MAPO ratio 0.82.
 * Known limitations (model form): an absolute offset of −0.14…−0.45 CR at ON 80–100 (within the ±0.4 CR
 * counter→CR spread between ANL campaigns and the 0.66–1.43 bar spread of measured standard-knock MAPO,
 * cfr_validation_README.md) because the model's MAPO varies too little between standard-knock states;
 * the low-ON RON scale (PRF 60–80) knocks 0.45–0.94 CR too early (PRF 80 rated RON 63 / MON 76); the λ
 * trend is too weak. The plain Douaud–Eyzat and the LLNL tables (single-
 * or two-stage) remain available as options.
 */
export const CFR_KNOCK_DELAY_MODEL = 'douaud-eyzat-llnl' as const;

/** A per-engine calibration set: the same parameters as {@link CFR_CALIBRATION}. */
export type CalibrationSet = { readonly [K in keyof typeof CFR_CALIBRATION]: CalibratedParameter };

/**
 * Transfer one CFR parameter to another engine: same value, range and source; the evidence states that
 * the value is NOT fitted to that engine and names the data that would set it.
 */
function transferred(p: CalibratedParameter, wouldBeSetBy: string): CalibratedParameter {
  return Object.freeze({
    value: p.value,
    range: p.range,
    source: p.source,
    evidence: `UNVERIFIED for this engine: the CFR F-1 value transferred unchanged (nothing fitted). Would be set by: ${wouldBeSetBy}`,
  });
}

/**
 * The Model T set: the CFR closures plus the ones only the Model T calibration uses — the turbulence
 * dissipation-length multiplier (options.turbulenceLengthScaleFactor; neutral 1 for the CFR) and the two
 * friction factors added to PNH (options.ts MODEL_T_CALIBRATED_FRICTION = MODEL_T_FRICTION with these).
 */
export type ModelTCalibrationSet = CalibrationSet & {
  readonly turbulenceLengthScaleFactor: CalibratedParameter;
  readonly ringTensionFactor: CalibratedParameter;
  readonly auxiliaryFactor: CalibratedParameter;
};

/**
 * Ford Model T calibration (calibration phase, 2026-10). No Model T cylinder-pressure trace exists, so the
 * set is fitted to GLOBAL period data, ONE set for every operating point of the engine (DESIGN.md policy):
 *
 *   1. Upton's WOT MBT spark advance a₀ = 0.108R/(1 + 0.001R) on a Cornell Ford engine (J. SAE Aug. 1923,
 *      "reduced to zero intake suction"; data to ≈ 1500 rpm) — sets the burn (burnRateMultiplier,
 *      taylorScaleMultiplier, turbulenceLengthScaleFactor). Simulated as Upton ran it: an ideal timed spark
 *      (his Atwater Kent battery-distributor ignition — the CFR inductive coil on the Model T spec), intake
 *      air and discharge water at 140 °F (333 K; mixture 318 K after fuel evaporation), barometer 29.35 inHg,
 *      wide open (engines/model-t.ts MODEL_T_UPTON_CONDITIONS); the model's MBT is the parabola vertex of its
 *      brake-torque-vs-spark scan (2.5° steps).
 *   2. Upton's Table 3 (same paper) intake suction at the lowest-suction (≈ wide-open) point of each speed:
 *      1.10 / 2.60 / 4.40 / 5.55 inHg at 600 / 800 / 1200 / 1400 rpm — sets the carburettor restriction
 *      (venturiDischargeCoefficient on the 23/32 in proxy throat). Cross-check: Upton's Fig. 4 peak BMEP
 *      ≈ 70 psi at 4 inHg (≈ 1100 rpm by Table 3) = Ford's 69 psi at 1000–1100 rpm.
 *   3. Ford's WOT brake-torque table (FSB Fig. 84, transmission output, "representative of the motors in
 *      general use") — the model's MAXIMUM brake torque over the spark-lever range (magneto, lever scanned
 *      in 2.5° steps: the lever→first-spark map is a magneto staircase) at Ford's dyno condition
 *      (engines/model-t.ts MODEL_T_FORD_DYNO) — sets the friction added to PNH (ringTensionFactor,
 *      auxiliaryFactor; MODEL_T_CALIBRATED_FRICTION in options.ts) and checks the shape.
 *   The remaining entries are kept at physically argued values and checked, not fitted.
 * Scripts: tmp/calib/ of the calibration worktree (sweep-lib.ts runner, gen.py / mk_grid2.py grids, rank.py),
 * ≈ 5000 fixed-speed runs; grid over C_T 2–6, C_λ 1.6–3, f_L × 0.35–1, venturi C_D 0.6–0.95, friction
 * variants. Warm-up 2 + 2 cycles; checked against 5 + 4 cycles (brake torque within 0.01 lb-ft).
 *
 * Before (transferred CFR set, venturi 0.95, pure PNH friction) → after:
 *   Upton MBT (model − Upton) 400/600/800/1000/1200 rpm: −16.6/−20.1/−23.7/−26.6/−30.0° → +2.5/+2.0/−0.3/−2.2/−2.7°
 *     (1400 rpm: flat optimum 63–73°, earlier inductive sparks misfire); Table 3 measured optima at the WOT points
 *     600/800/1200/1400 rpm: 42/53/55/60° vs model 42.5/47.7/56.2/63–73°. BMEP at the MBT spark, Upton's condition:
 *     67–68 psi at 600–800 rpm (his Fig. 4: ≈ 70 psi).
 *   WOT suction 600/800/1200/1400 rpm: 0.85/1.44/2.67/3.25 → 1.61/2.60/4.60/5.50 inHg (Upton 1.10/2.60/4.40/5.55).
 *   Ford WOT max brake torque, model/Ford − 1, 500/600/…/1900 rpm: before +47/…/+36 (700)/…/+27 (900)/…/+33
 *     (1200)/…/+37 (1500)/…/+52 (1800) % → +29.4/+25.5/+19.0/+14.2/+9.0/+9.2/+7.6/+5.5/+1.8/−2.3/−5.6/−8.2/
 *     −4.9/−4.5/−5.9 %: 900–1900 rpm within ±10 %; the 500–800 rpm excess is a documented model-form residual
 *     (test/fixtures/modelt_validation_README.md); peak 92.8 lb ft at 700 rpm (Ford 83 at 900), 19.3 hp at
 *     1400 rpm (Ford 20 hp at 1500–1600).
 *   Top speed (fixed-speed WOT curve, lever 64.5, against the 'vehicle' road load): 49.5 → 43.6 mph at 1784 rpm.
 *   Motored compression at 150 rpm: 67.9 psig before and after (Ford 60, Motor Age 55): a residual, see
 *     modelt_validation_README.md (no blow-by in the model; CR 3.78 would give 62.8 psig).
 *
 * Physics of the fitted burn: the L-head's turbulence is weak ("the air entering through the inlet valves
 * had to turn two right angles before it entered the cylinder and … lost much of its initial velocity",
 * Ricardo, The High-Speed IC Engine, ch. 6, on the pre-1919 side-valve slab chamber with the plug over the
 * inlet valve — the Model T's) and decays fast in the shallow pocket: with the fitted dissipation length the
 * K–k u′ at TDC is ≈ 0.25 S̄p (CFR ≈ 0.6) and the entrainment speed C_T u′ is ≈ 1/3 of Keck's (1982) eq. 4.10
 * correlation evaluated on the Model T's inlet flow, so S_L carries a large share of the burn and the
 * burn angle grows with speed as Upton measured (a₀ ∝ R/(1 + 0.001R)). Because u′ decays more in the
 * longer time of a low-speed cycle, u′ rises faster than rpm, which keeps the 1500–1900 rpm burn short
 * enough for the advance the magneto/trembler ignition can reach (≈ 51–64° first spark).
 */
export const MODEL_T_CALIBRATION: ModelTCalibrationSet = Object.freeze({
  burnRateMultiplier: Object.freeze({
    value: 4.0,
    range: [1.2, 16] as const,
    source:
      'Both ends evaluated on the Model T model’s own burn (K–k u′ at the calibrated length scale f_L = 0.4 × 0.25, ' +
      'S_L at CA10–CA50, 400–1800 rpm WOT): lower end Upton 1923 (J. SAE XIII(2) Fig. 10, the Ford engine’s ' +
      '"turbulence factor" S_T/S_L = 1 + 0.001R, i.e. u_T = 0.001R·S_L) → u_T/u′ ≈ 1.2–1.5; upper end Keck 1982 ' +
      '(19th Symp. Combust. p. 1451) eq. 4.10, u_T = 0.08 ū_i (ρ_u/ρ_i)^½ with ū_i from the unshrouded inlet valve at ' +
      'maximum lift → u_T/u′ = 10.7–15.5 (×1.1 for his ±10 %). Qualitatively (Ricardo, High-Speed IC Engine ch. 6) ' +
      'the pre-turbulent-head side-valve chamber "lacked turbulence".',
    evidence:
      'Upton MBT (ideal spark, MODEL_T_UPTON_CONDITIONS) model − Upton: 400/600/800/1000/1200 rpm +2.5/+2.0/−0.3/−2.2/−2.7° ' +
      '(rms 2.1°; 1400 rpm flat optimum 63–73°). Grid (C_T × C_λ × f_L, rank.py): C_T 3.5/4.5 with f_L 0.45/0.35 ' +
      'fit equally (rms 1.8/2.3°); C_T 2 with f_L 1 (the CFR length scale) fits the MBT (rms 1.2°) but its early-spark ' +
      'kernel waits for the hand-off radius (∝ L) to shrink and the WOT torque collapses above 1400 rpm ' +
      '(1800 rpm −24 % vs Ford). CFR transfer 5.3 (f_L 1): MBT 17–30° too early.',
  }) as CalibratedParameter,
  taylorScaleMultiplier: Object.freeze({
    value: 2.0,
    range: [1.15, 3.8] as const,
    source:
      'Keck 1982 Fig. 15, ℓ_T = 0.8 L_IV (ρ_i/ρ_u)^¾ (±25 %), evaluated on the Model T burn (L_IV = 5.70 mm net lift, ' +
      'raw Taylor microscale at CA50, 400–1800 rpm WOT, calibrated f_L): ℓ_T/λ = 1.5–3.0 (range × 0.75…1.25).',
    evidence:
      'Upton MBT with C_T/f_L (burnRateMultiplier): C_λ 1.6/2.0/2.5 at C_T 4, f_L × 0.4 give MBT rms 3.2/2.1/2.4° ' +
      '(calibration grid g3, venturi 0.62) and Ford 1500/1800 rpm +8.3/+15.7 %, +6.0/+11.9 %, +2.8/+7.1 % with pure ' +
      'PNH friction: 2.0 is the MBT optimum; C_λ moves the high-speed tail (burn-up time ∝ λ/S_L) more than the MBT.',
  }) as CalibratedParameter,
  turbulenceLengthScaleFactor: Object.freeze({
    value: 0.4,
    range: [0.4, 1.5] as const,
    source:
      'Multiplier on f_L = 0.25 in L = f_L·min(h, B) (combustion/turbulence.ts): L_int/h ≈ 0.10–0.15 (Aleiferis & ' +
      'Behringer 2017, Fuel 189:238) divided by C_ε = εL_int/u′³ ≈ 0.4–1.0 (shear flow ≈ 0.5, grid turbulence 0.7–1; ' +
      'UNVERIFIED, search summaries — turbulence.ts) → f_L = 0.10–0.375, i.e. × 0.4–1.5.',
    evidence:
      'At the lower bound (L_int/h ≈ 0.10 with grid-like C_ε ≈ 1: the bore column and the 13 mm deep valve pocket). ' +
      'With C_T 4, C_λ 2: Upton MBT rms 2.1° and Ford WOT 1500/1800 rpm −5.6/−4.5 %. f_L × 0.5/0.65/1.0, each with C_T ' +
      'refitted to the MBT, give 1800 rpm −5/−17/−24 % (calibration grids g2/w1, 0.6–1 bar more or less friction ' +
      'aside): the smaller L shortens the hand-off radius r_ho = 2·0.5·L of the early (50–64° BTDC) magneto sparks — ' +
      'at f_L 1 the kernel stalls at ≈ 2 mm for ≈ 30° at 1800 rpm — and makes u′ grow faster than rpm. ' +
      'u′(TDC)/S̄p ≈ 0.25 (Ricardo: weak side-valve turbulence).',
  }) as CalibratedParameter,
  kernelHandoffMultiple: transferred(
    CFR_CALIBRATION.kernelHandoffMultiple,
    'universal closure (Forte C_m1 = 2, not engine-fitted); C = 1 (Fluent) changed the calibrated WOT torque by ≤ 5 % ' +
      '(1800 rpm) and the MBT by < 1° in the calibration grid (with C_T 2, f_L 1), so it was left universal.',
  ),
  marksteinMultiplier: transferred(
    CFR_CALIBRATION.marksteinMultiplier,
    'universal closure (Bradley 1998 measured/theory ratio, not engine-fitted); unchanged.',
  ),
  woschniMultiplier: Object.freeze({
    value: CFR_CALIBRATION.woschniMultiplier.value,
    range: CFR_CALIBRATION.woschniMultiplier.range,
    source: CFR_CALIBRATION.woschniMultiplier.source,
    evidence:
      'Kept at the CFR value (per-area coefficient; the L-head pocket and block-deck areas are in the chamber geometry). ' +
      'Check against Ricardo’s heat balance (High-Speed IC Engine ch. 5: a well-designed 5:1 engine loses ≈ 13 % of ' +
      'the fuel heat to the walls during combustion + expansion): the calibrated Model T loses 20–25/16/14.5 % (closed ' +
      'cycle, LHV basis, φ 1.15) at 900/1500/1800 rpm WOT — already above it, as a slow-burning L-head should; ×2.0 ' +
      'would make it 28 % at 900 rpm and lowers the WOT torque ≈ 15 % at every speed (no shape change), so the torque ' +
      'level was assigned to friction (transmission-output basis) instead.',
  }) as CalibratedParameter,
  intakePortHeatTransferMultiplier: Object.freeze({
    value: CFR_CALIBRATION.intakePortHeatTransferMultiplier.value,
    range: CFR_CALIBRATION.intakePortHeatTransferMultiplier.range,
    source: CFR_CALIBRATION.intakePortHeatTransferMultiplier.source,
    evidence:
      'Kept at the CFR value (the siamesed ports are cast in the block next to the jacket; spec.manifolds.intakePortLength ' +
      '0.1 m UNVERIFIED). Sensitivity: 0 (adiabatic) raises the WOT torque by 3.5 % at 500–700 rpm, 2–3 % at 900–1200, ' +
      '≈ 0 at 1500 — too weak to carry any calibration target, so not refitted.',
  }) as CalibratedParameter,
  venturiDischargeCoefficient: Object.freeze({
    value: 0.62,
    range: CFR_CALIBRATION.venturiDischargeCoefficient.range,
    source: CFR_CALIBRATION.venturiDischargeCoefficient.source,
    evidence:
      'Fitted to Upton 1923 Table 3 (engines/model-t.ts MODEL_T_UPTON_TABLE3: the lowest intake suction at each speed, ' +
      'Cornell Ford engine with a "regular Holley" carburettor) at MODEL_T_UPTON_CONDITIONS and Upton’s MBT spark: ' +
      'cycle-mean manifold suction 1.61/2.60/4.60/5.50 inHg at 600/800/1200/1400 rpm vs 1.10/2.60/4.40/5.55 (least ' +
      'squares over 800–1400 rpm: 0.625; the 600 rpm point is ≥ 0.4 inHg high at any C_D ≤ 0.7). C_D 0.6/0.64/0.66/0.7/0.95: ' +
      'mean model/Upton ratio 1.057/0.971/0.933/0.863/0.58. Sensitivity to the UNVERIFIED 318 K mixture temperature of ' +
      'Upton’s condition: a 300 K mixture moves the fit to 0.64. The value multiplies the area of the ' +
      'UNVERIFIED 23/32 in Holley G proxy throat (manifolds.venturiDiameter), so it is the effective C_D·A of the Holley NH ' +
      'venturi + fuel-nozzle bridge + swayback passage + 1 in butterfly that the data measure, at Re ≈ 0.6–1.1e5.',
  }) as CalibratedParameter,
  knockStratificationDT: transferred(
    CFR_CALIBRATION.knockStratificationDT,
    'knock intensity — no Model T measurement exists. Ricardo (High-Speed IC Engine ch. 6): side-valve engines at ' +
      '≈ 4:1 on the 45–50 ON petrol of 1919–21 "would detonate heavily unless the ignition timing were constantly ' +
      'adjusted"; the calibrated model knocks at WOT up to 1400 rpm with the lever fully advanced (MAPO 0.01–0.42 bar at ' +
      'the virtual plug transducer; a light-knock lever, MAPO ≤ 0.05 bar, gives 2–15 % less torque at 500–900 rpm) ' +
      'and not at the cruise settings — onset is the delay model’s, the MAPO scale unvalidated.',
  ),
  ringTensionFactor: Object.freeze({
    value: 2.5,
    range: [1, 3.2] as const,
    source:
      'PNH ring-tension term ∝ ring tension (Patton et al. 1989 / Sandoval 2002 F_t/F_t0). Model T: three cast-iron ' +
      'rings 1/4 in wide (groove 1/4 × 13/64 in, eccentric 0.180→0.085 in [Dyke 1924 Instr. 70]) = 19.05 mm of ring face ' +
      'vs ≈ 6–7.5 mm (two ≈ 1.5 mm compression rings + a 3–4 mm oil ring) in the 1980s engines PNH was fitted to ' +
      '(UNVERIFIED typical widths) → × 2.5–3.2 at equal wall pressure; 1 = PNH as fitted.',
    evidence:
      'Physical estimate (face-width ratio 19.05/7.5), not fitted: + 0.20/0.14/0.11 bar FMEP at 500/900/1500 rpm. Its ' +
      '(1 + 1000/N) boundary-friction form is the only low-speed-heavy loss available; a factor at the bound (3.2) moves ' +
      'the 500 rpm torque by −1.8 % only (the 500–800 rpm excess is not a friction-sized residual).',
  }) as CalibratedParameter,
  auxiliaryFactor: Object.freeze({
    value: 3.0,
    range: [0.3, 3.5] as const,
    source:
      'Multiplier on PNH’s auxiliary term (oil + water pump + alternator of a 1980s engine: 0.49 kW at 1500 rpm for 2.9 L). ' +
      'Model T: no pumps — belt fan, generator (1919+) and magneto load ≈ 0.15 kW → 0.3 (MODEL_T_FRICTION, UNVERIFIED); ' +
      'Ford measured at the transmission output [Tulsa], so the high-gear churning of the magneto flywheel, clutch drum ' +
      'and planetary in the oil pit belongs here too: UNVERIFIED estimate ≤ 1.5 kW at 1500 rpm (≈ 5–7 % of the 20 hp; ' +
      'Tulsa assumes 20–25 % for transmission + axle on chassis dynos) → upper end 3.5 (1.7 kW).',
    evidence:
      'Fitted to the Ford WOT level at 900–1900 rpm with ringTensionFactor 2.5 (minimax error): max-over-lever torque ' +
      '+9.0/+9.2/+5.5/−5.6/−8.2/−4.5 % at 900/1000/1200/1500/1600/1800 rpm, every 900–1900 rpm point within ±10 % (aux 2.5: ' +
      '1000 rpm +10.4 %; aux 3.5: 1600 rpm −10.0 %; offline PNH re-evaluation of the same cycles). Total FMEP ' +
      '1.00/1.11/1.18 bar at 900/1500/1800 rpm, η_m 0.78 at 1500 rpm — the period ALAM/SAE rating assumption is η_m 0.75 ' +
      'at 1000 ft/min (90 psi IMEP, D²N/2.5 [Good 1922 pp. 37–38]); pure PNH gave 0.63 bar, η_m 0.88, and +5…+18 % at ' +
      '900–1900 rpm. Aux 3.0 ≈ 1.5 kW at 1500 rpm.',
  }) as CalibratedParameter,
});
