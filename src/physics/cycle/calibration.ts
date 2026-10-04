/**
 * The ONE global calibration of the cycle model for the Waukesha CFR F-1 (DESIGN.md calibration
 * policy): every genuinely uncertain closure that is tuned, with its value, the literature range it
 * must stay in, the source of that range, and the validation evidence that set it. The same set
 * serves every operating point (RON, MON, all octane numbers, spark timings, fuels); nothing is
 * tuned per operating point. DEFAULT_CYCLE_MODEL_OPTIONS (options.ts) takes its defaults from here.
 * Multi-engine (Model T integration): one set PER ENGINE, selected by EngineSpec.id through the
 * engine-keyed option defaults of options.ts — CFR_CALIBRATION ('cfr-f1') and MODEL_T_CALIBRATION
 * ('ford-model-t', transferred CFR values, nothing fitted yet; see below).
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
 * Ford Model T calibration (Model T integration, 2026-10). No measured Model T cylinder-pressure traces
 * exist, so NOTHING here is fitted: every value is the CFR F-1 value transferred (DESIGN.md calibration
 * policy, one set per engine; the universal closures marksteinMultiplier and kernelHandoffMultiple are the
 * same physics on both engines). Each entry names the Model T evidence that would set it. The global
 * targets are Ford's WOT torque/power table (engines/model-t.ts MODEL_T_FORD_WOT_TABLE, ±1.5 % rounding,
 * transmission-output basis) and Upton's MBT spark advance vs speed (modelTUptonMbtDeg). Reasons the CFR
 * values should NOT be expected to carry over: the CFR burn constants were fitted with a shrouded overhead
 * intake valve and a disc chamber, while the Model T has an unshrouded side valve feeding an L-head pocket
 * (weaker, differently structured turbulence, a longer flame path), a much larger surface/volume ratio, a
 * thermosyphon jacket, and a Holley/Kingston carburettor venturi of another size and Reynolds number.
 */
export const MODEL_T_CALIBRATION: CalibrationSet = Object.freeze({
  burnRateMultiplier: transferred(
    CFR_CALIBRATION.burnRateMultiplier,
    "Upton's MBT spark advance vs speed on a Ford engine at WOT (J. SAE 1923, 0.108R/(1 + 0.001R) deg: MBT puts CA50 " +
      'near 8–10° ATDC, so the MBT curve measures the burn duration and its growth with rpm) together with the shape of ' +
      "Ford's WOT torque curve (FSB Fig. 84) at a fixed volumetric efficiency.",
  ),
  taylorScaleMultiplier: transferred(
    CFR_CALIBRATION.taylorScaleMultiplier,
    'the same MBT-vs-rpm data as C_T (C_λ sets the burn-up tail CA50→CA90, which moves the MBT angle at low speed).',
  ),
  kernelHandoffMultiple: transferred(
    CFR_CALIBRATION.kernelHandoffMultiple,
    'universal closure (Forte C_m1 = 2, not engine-fitted); a trembler spark shower (several breakdowns per event) may ' +
      'need its own hand-off criterion: check the misfire limit at idle (≈ 400 rpm, throttle nearly shut, spark retarded).',
  ),
  marksteinMultiplier: transferred(
    CFR_CALIBRATION.marksteinMultiplier,
    'universal closure (Bradley 1998 measured/theory ratio, not engine-fitted); unchanged.',
  ),
  woschniMultiplier: transferred(
    CFR_CALIBRATION.woschniMultiplier,
    "the WOT brake efficiency (Ford's dyno table with a fuel-flow measurement, or period road-test fuel economy at a " +
      'known speed) and the jacket heat rejection of the thermosyphon system; the L-head pocket adds surface that the ' +
      'per-area Woschni coefficient does not know about.',
  ),
  intakePortHeatTransferMultiplier: transferred(
    CFR_CALIBRATION.intakePortHeatTransferMultiplier,
    'the WOT torque level at 500–900 rpm (FSB Fig. 84: 69–83 lb ft), which fixes the volumetric efficiency once the burn ' +
      'is set; the siamesed intake ports run through the water jacket of the block (spec.manifolds.intakePortLength UNVERIFIED).',
  ),
  venturiDischargeCoefficient: transferred(
    CFR_CALIBRATION.venturiDischargeCoefficient,
    'the WOT torque fall-off above 1000 rpm (FSB Fig. 84: 82 → 47 lb ft at 1000 → 1900 rpm), mostly carburettor and ' +
      'manifold restriction; a measured air flow of a Holley NH / Kingston L-4 venturi would set it directly.',
  ),
  knockStratificationDT: transferred(
    CFR_CALIBRATION.knockStratificationDT,
    'knock intensity — no measurements exist; the period practice of retarding the spark lever on hills (audible spark ' +
      'knock at WOT and low speed on ≈ 40–55 ON gasoline) only bounds the knock ONSET (the delay model), not the intensity.',
  ),
});
