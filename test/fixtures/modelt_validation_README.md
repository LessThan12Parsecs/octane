# Ford Model T validation data

Validation data for the Ford Model T engine model (`src/physics/engines/model-t.ts`, a 1924–25 high-head engine).
Machine-readable: `modelt_period_data.json`, written by `tools/reference/modelt_data_period.py`. Every number is a
transcription of a printed value or a closed-form relation printed by its author, so nothing is digitised from
a figure. Unlike the CFR data (`cfr_validation_README.md`), no in-cylinder pressure trace of a Model T is known
in open sources. Validation is therefore global: brake torque and power, spark timing, compression pressure,
ignition-circuit behaviour and vehicle performance.

## How the data were gathered

Five research agents swept period Ford service literature, period engineering journals, club technical pages
and restorer measurements in October 2026. Five independent skeptic agents then re-checked every physics-relevant
number against a different source, and corrected or flagged it where needed. The facts table for each area, with
the verdicts, is kept with the session notes. The dispositions that matter for modelling are below.

| quantity | value used | evidence | disposition |
|---|---|---|---|
| bore × stroke, rod | 3.750 × 4.000 in, 7.000 in | Ford 1923 Data Book; Good 1922; Dyke's 1924; Pagé 1929; Fahnestock; Tulsa piston-position table back-solves to 7.000 in | confirmed |
| compression ratio | 3.98 | McCalley (Ford records) 3.98; Gunnell 4.0; Allen 1987 4.0. Dyke's 1924 gives 3.6; Tulsa measured 3.8 on a 1917–18 head. The post-9/1918 chamber build-up gives 3.985 | **disputed** (3.6–4.0). 3.98 is used, with a 3.78–3.98 sensitivity range |
| crown height at TDC | +5/16 in above the deck | Ford Manual 1919 A22; Good 1922; Pagé 1929 | confirmed |
| crown-to-head distance at TDC | ≈ 1 in | Ford Service Bulletin (full-text snippet) | single source |
| gasket opening | 18.3 in² per cylinder | Tulsa 'Head Design' | single source |
| valve timing | IVO 12.7° ATDC, IVC 50.8° ABDC, EVO 37.9° BBDC, EVC at TDC | Ford piston positions (FM19, DB23, Pagé 1929) converted with L = 7, r = 2 | confirmed |
| cam lobe | three-arc, base 0.406 in, flank 1.260 in, nose 0.031 in, rise 0.250 in; lash 0.0256 in → 218°, 0.225 in | valvetrain research, verified | confirmed |
| coil | 3.3 mH / 0.295 Ω / 22 H / 3300 Ω / 0.40–0.45 µF; 212/16,600 turns | Ford 1916 data via Boggess & Patterson; ECCT; oscillograms | confirmed (k derived) |
| magneto | 8 cycles/rev, ≥ 7 V at 400 rpm, ≈ 0.0235 V_rms/rpm | Dyke's 1924; Ford service spec; restorer measurements | confirmed (impedance estimated) |
| timer | 87° contact, make 15.5° ATDC (full retard) to 64.5° BTDC, 80° lever | Patterson (Ford gauge 2-1/2 in, 1919–27) | single measurement |
| period fuel | PRF ≈ 45 surrogate (≈ 40–55 MON) | Bureau of Mines TP 328 distillation; Ricardo 1921 HUCRs (interpolated) | estimate |
| plug location | over the valve pocket | MTFCA forum (L. Young 2016); Green Engineering head is the unusual over-piston design | lower trust |
| bore spacing | 4-1/8, 5-1/4, 4-1/8 in | MTFCA forum snippet | **unverified** |

## Datasets (`modelt_period_data.json`)

| key | content | grade |
|---|---|---|
| `ford_wot_1918` | Ford's WOT torque/hp table, 300–1900 rpm, at the transmission output, described as "representative of the motors in general use" (best engines up to 22.5 hp). The 1400 rpm torque is back-computed from the hp column. | B |
| `upton_mbt_1923` | Upton's MBT relation a₀ = 0.108R/(1 + 0.001R): 31° at 400 rpm, 48° at 800, 59° at 1200 (data up to ≈ 1500 rpm). | C |
| `compression_pressure` | 60 psig (Ford), 55 psig (Motor Age 1922), 64 psig at CR 3.6 (Dyke's). | C |
| `trembler_coil` | Coil constants; firing time and current at 6/9/12 V; reclose time; HCCT currents; peak voltages; stored energy. | B |
| `magneto` | Service minimum; restorer open-circuit voltages; the 1976 V/A/Hz table; Patterson's 600 rpm spark ladder. | C |
| `timer_and_lever` | Contact arc, lever range, notches. | C |
| `vehicle_performance` | 40 rpm/mph in high, top speed 42–45 mph, grade consistency checks, fuel economy. | C/D |

## Ranked validation targets

1. **Ford WOT brake torque curve.** Its shape (peak 83 lb-ft at 900 rpm, roughly linear fall to 65 lb-ft at 1600)
   and its level (20 hp at 1500–1600) should be matched within ≈ ±10–15 %. That band covers the in-service vs
   best-engine spread (22.5 hp) and the unknown spark and fuel of the test.
2. **Upton MBT** up to 1500 rpm: the model's best-torque spark at WOT. Note that the stock ignition cannot reach
   MBT above ≈ 1000 rpm, because the timer make plus the coil firing time limits the advance.
3. **Trembler-coil firing time and current** vs supply voltage (3.5 ms / 2.5 ms / 2.0 ms at 6 / 9 / 12 V), and
   the magneto spark ladder at 600 rpm (22.5° steps).
4. **Compression pressure** 55–60 psig, motored at cranking speed and wide-open throttle.
5. **Vehicle**: top speed 42–45 mph in high on a level road, and 40 rpm per mph.

## Caveats

- The Ford table is a test of an "average" engine c. 1918, run with period gasoline and an unknown spark setting.
- No measured friction, volumetric efficiency, heat-transfer or pressure-trace data exist. Friction comes from
  the PNH model extrapolated to a 1920s babbitt-bearing engine, so the brake numbers inherit its uncertainty.
- Knock: the end-gas model is weakest for PRF ≤ 70 (DESIGN.md), and that is where period fuel sits. At CR ≈ 4,
  the Model T was not knock-limited on 1920s petrol (Ricardo's HUCRs were 4.3–6.0). Treat model knock at normal
  lever settings as a model-form artefact to investigate, not a validated prediction.

## Calibration results (calibration phase, October 2026)

The model was calibrated against the period data with ONE parameter set for every operating point
(`src/physics/cycle/calibration.ts` `MODEL_T_CALIBRATION`; each entry has its value, literature range, source
and evidence). The checks are `test/validation/measured_modelt_data.test.ts` (gated by `OCTANE_VALIDATION=1`)
and `test/validation/modelt_calibration.test.ts` (fast). The CFR F-1 calibration is unchanged: a 22-case IEEE-bit
golden of CFR runs is identical before and after.

One point corrects the caveat above. Ricardo (The High-Speed IC Engine, ch. 6) writes that pre-1919 side-valve
engines at ≈ 4:1 on the 45–50 ON petrol of the day "would detonate heavily unless the ignition timing were
constantly adjusted". Knock at WOT with the lever fully advanced and low speed is therefore realistic: it is
spark knock, which period drivers removed by retarding the lever (Ford Manual). Ricardo's HUCRs were measured at
optimum spark in his own engine.

### What was fitted, and to what

| parameter | before | after | range | set by |
|---|---|---|---|---|
| burnRateMultiplier C_T | 5.3 (CFR) | 4.0 | 1.2–16 (Upton 1923 turbulence factor … Keck 1982 eq. 4.10, both on the Model T's own burn) | Upton's MBT, 400–1200 rpm |
| taylorScaleMultiplier C_λ | 3.0 (CFR) | 2.0 | 1.15–3.8 (Keck 1982 Fig. 15 on the Model T burn) | Upton's MBT |
| turbulenceLengthScaleFactor | 1 | 0.4 | 0.4–1.5 (L_int/h 0.10–0.15 over C_ε 0.4–1) | Upton's MBT + the 1500–1900 rpm torque (early-spark kernel hand-off) |
| venturiDischargeCoefficient | 0.95 (estimate) | 0.62 | 0.6–0.984 | Upton's Table 3 WOT intake suction, 800–1400 rpm |
| ringTensionFactor (PNH) | 1 | 2.5 | 1–3.2 (ring-face width 19.05 mm vs ≈ 6–7.5 mm) | physical estimate, not fitted |
| auxiliaryFactor (PNH) | 0.3 | 3.0 | 0.3–3.5 (fan, generator, magneto, transmission churning) | Ford's WOT level at 900–1900 rpm; η_m 0.78 at 1500 rpm vs the ALAM 0.75 |
| woschniMultiplier, intakePortHeatTransferMultiplier, knockStratificationDT, kernelHandoffMultiple, marksteinMultiplier | CFR values | unchanged | — | checked, not refitted (see the entries) |

Data changes that came with the calibration (these are data, not calibration):
- Presets: the cruise lever 25 → 50; the full-throttle lever 40 → 64.5 (Pagé: "maximum speed: spark and
  throttle fully advanced"); the idle throttle 0.06 → 0.12 and lever 0 → 5.
- New sourced data in `engines/model-t.ts`: Upton's Table 3 (`MODEL_T_UPTON_TABLE3`) and his test condition
  (`MODEL_T_UPTON_CONDITIONS`).
- The compression ratio stays 3.98, the chamber build-up of the 1924–25 head.

### 1. Ford WOT brake torque (max over the lever range, magneto, `MODEL_T_FORD_DYNO`)

The lever was scanned from 15 to 64.5 in 2.5° steps (2 warm-up + 2 cycles; 5 + 4 cycles agree to 0.01 lb-ft).
The last column is a diagnostic, not the target: the best torque with at most light knock (MAPO ≤ 0.05 bar at the
virtual plug transducer), i.e. the lever a period tester would use to stop audible detonation.

| rpm | Ford lb-ft | model lb-ft | error % | lever | first spark °BTDC | η_v | IMEP bar | FMEP bar | PMEP bar | CA50 °ATDC | knock cyl/MAPO bar | light-knock max (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 69 | 89.3 | +29.4 | 42.0 | 32.8 | 0.790 | 6.27 | 1.01 | 0.11 | 17.7 | 4/0.15 | 80.2 (+16.3) |
| 600 | 73 | 91.6 | +25.5 | 47.0 | 38.4 | 0.793 | 6.39 | 0.99 | 0.15 | 13.9 | 4/0.19 | 77.5 (+6.2) |
| 700 | 78 | 92.8 | +19.0 | 54.5 | 40.8 | 0.794 | 6.45 | 0.99 | 0.19 | 13.3 | 4/0.19 | 88.7 (+13.7) |
| 800 | 81 | 92.5 | +14.2 | 54.5 | 40.5 | 0.793 | 6.43 | 0.99 | 0.24 | 15.6 | 4/0.15 | 86.3 (+6.5) |
| 900 | 83 | 90.5 | +9.0 | 54.5 | 40.2 | 0.787 | 6.32 | 1.00 | 0.28 | 17.7 | 4/0.10 | 88.5 (+6.6) |
| 1000 | 82 | 89.5 | +9.2 | 64.5 | 53.5 | 0.772 | 6.28 | 1.01 | 0.32 | 2.7 | 4/0.34 | 86.7 (+5.7) |
| 1100 | 81 | 87.2 | +7.6 | 64.5 | 53.1 | 0.754 | 6.15 | 1.02 | 0.37 | 5.2 | 4/0.27 | 81.8 (+1.0) |
| 1200 | 79 | 83.4 | +5.5 | 64.5 | 52.7 | 0.733 | 5.95 | 1.04 | 0.41 | 7.9 | 4/0.18 | 76.5 (−3.2) |
| 1300 | 77 | 78.4 | +1.8 | 64.5 | 52.3 | 0.709 | 5.67 | 1.06 | 0.45 | 10.8 | 4/0.09 | 70.7 (−8.2) |
| 1400 | 74 | 72.3 | −2.3 | 64.5 | 51.9 | 0.685 | 5.34 | 1.08 | 0.50 | 14.5 | 4/0.01 | 72.3 (−2.3) |
| 1500 | 70 | 66.1 | −5.6 | 64.5 | 51.5 | 0.660 | 4.99 | 1.11 | 0.54 | 18.2 | 0 | 66.1 (−5.6) |
| 1600 | 65 | 59.7 | −8.2 | 64.5 | 51.1 | 0.635 | 4.64 | 1.13 | 0.57 | 22.0 | 0 | 59.7 (−8.2) |
| 1700 | 60 | 57.1 | −4.9 | 64.5 | 63.9 | 0.610 | 4.51 | 1.16 | 0.60 | 17.4 | 0 | 57.1 (−4.9) |
| 1800 | 53 | 50.6 | −4.5 | 64.5 | 63.9 | 0.586 | 4.16 | 1.18 | 0.63 | 21.8 | 0 | 50.6 (−4.5) |
| 1900 | 47 | 44.2 | −5.9 | 62.0 | 61.3 | 0.564 | 3.81 | 1.21 | 0.66 | 26.3 | 0 | 44.2 (−5.9) |

Before calibration (the same scan): +47 / +36 / +27 / +33 / +37 / +52 % at 500 / 700 / 900 / 1200 / 1500 /
1800 rpm. After:
- 900–1900 rpm are within ±10 %; the worst points are −8.2 % at 1600 rpm and +9.2 % at 1000 rpm.
- Peak power is 19.3 hp at 1400 rpm (Ford: 20 hp at 1500–1600).
- Peak torque is 92.8 lb-ft at 700 rpm (Ford: 83 at 900).
- The dip at 1600 rpm is the magneto staircase: the first spark stays at 51° up to 1600 rpm, and the timer-make
  spark at 64° only appears from ≈ 1650 rpm.

**What makes the curve fall above 1000 rpm.** From 900 to 1800 rpm the brake torque falls 44 % (BMEP −2.3 bar):
- **Breathing, about 70 %.** η_v falls from 0.79 to 0.59. The carburettor (effective C_D·A 1.63 cm²) and the
  side-valve ports cannot pass the flow: the WOT manifold suction reaches 4.6 inHg at 1200 rpm and 5.5 at 1400,
  as Upton measured.
- **Net indicated efficiency, about 23 %.** Pumping MEP rises from 0.28 to 0.63 bar. On top of that the burn is
  late:
  - the burn takes more crank degrees as speed rises (the model's own MBT is 52° at 1000 rpm and 63–73° at
    1400 rpm, as Upton's MBT rises with speed);
  - the stock ignition stops at a first spark of 51–53° (timer make 64.5° less the coil firing time), or 64° once
    the magneto's make spark appears above ≈ 1650 rpm;
  - CA50 therefore sits at 18–22° ATDC from 1500 rpm up.
- **Friction, about 8 %.** FMEP rises from 1.00 to 1.18 bar.

### 2. Upton 1923 (Cornell Ford engine; `MODEL_T_UPTON_CONDITIONS`, ideal timed spark)

| rpm | Upton a₀ ° | model MBT ° (before → after) | Table 3 optimum ° | WOT suction inHg: model (Upton Table 3) | BMEP at a₀ psi |
|---|---|---|---|---|---|
| 400 | 30.9 | 14.2 → 33.4 | — | 0.69 | 63.1 |
| 600 | 40.5 | 20.4 → 42.5 | 42 | 1.61 (1.10) | 66.9 |
| 800 | 48.0 | 24.3 → 47.7 | 53 | 2.60 (2.60) | 67.5 |
| 1000 | 54.0 | 27.4 → 51.8 | — | 3.46 | 64.8 |
| 1200 | 58.9 | 28.9 → 56.2 | 55 | 4.60 (4.40) | 58.6 |
| 1400 | 63.0 | — → flat 63–73 | 60 | 5.50 (5.55) | 51.8 |

- The "before" MBT values are the uncalibrated model at 101.3 kPa.
- At 1400 rpm the torque curve is flat (within 0.1 lb-ft) from 63 to 73°, and inductive sparks earlier than 76°
  misfire: the early-spark kernel hands off late.
- Upton's Fig. 4 peak BMEP is ≈ 70 psi at 4 inHg (read by eye); the model gives 67–68 psi at 600–800 rpm.
- At low speed the model's peak pressure at MBT comes later than Upton's theoretical a₀/3 ATDC (34° vs 10° at
  400 rpm). That is his theory, not a measurement.

### 3–6. Compression, vehicle, knock, idle

- **Motored compression at 150 rpm, WOT:** 67.9 psig (66.9 at 100 rpm, 68.5 at 200 rpm), against Ford's 60 and
  Motor Age's 55. This is **a residual, not fitted**:
  - The model has no blow-by. Real engines read 35–45 psig when hand-cranked and 55–60 at starter speed, so
    leakage depends strongly on speed; the model barely changes with speed.
  - A gauge's dead volume (≈ 5–10 cm³ on 243 cm³) takes off another 1–3 psi.
  - CR 3.78 (the low end of the evidence) would give 62.8 psig, and a Woschni multiplier of 2 would give 61.6.
  - The CR was kept at 3.98, the chamber build-up of the 1924–25 head. Tulsa's 3.8 was measured on a 1917–18 head,
    before Ford reduced the chamber by 1/16 in.
- **Top speed** (fixed-speed WOT torque at lever 64.5 against the road load): 43.6 mph at 1784 rpm (target 42–45;
  49.5 mph before).
- **Cruise** (throttle 0.35, lever 50): 30.5 mph at 1250 rpm, first spark 39° BTDC, no knock. The old lever 25 sits
  on the magneto's 13–16° BTDC step and gives 25.9 mph with CA50 ≈ 60° ATDC.
- **Knock** (PRF 45; Ricardo gives 45–50 ON for 1919–21 petrol):
  - At WOT with the lever fully advanced, all four cylinders knock from 500 to 1400 rpm (MAPO ≈ 0.4 bar at
    500–900 rpm, 0.01 at 1400).
  - Nothing knocks from 1500 rpm up, or at the cruise settings at 1100 rpm and above. At 900 rpm there is a trace
    (MAPO 0.002 bar).
  - The onset comes from the delay model (which knocks early for PRF ≤ 70); the MAPO scale is not validated.
- **Idle** (preset throttle 0.12, lever 5, first spark ≈ 4° ATDC): 405 rpm declutched, no misfire. Throttle 0.07
  gives 230 rpm, 0.09 gives 300, 0.11 gives 370 and 0.13 gives 440; at 0.05 the engine stalls.

### Documented disagreements (`it.fails` in the gated suite)

1. **WOT torque at 500–800 rpm is +14…+29 % above Ford.** This holds for every burn set in the calibration grid
   (C_T 2–6, f_L × 0.35–1). Friction cannot close it: the ring factor at its bound moves the 500 rpm torque by
   −1.8 %, and reaching +10 % would take ≈ 0.8 bar more FMEP at 500 rpm (η_m ≈ 0.7 at WOT). Likely causes, none of them in the model:
   - **Mixture preparation.** At low air speed the heavy 1920 gasoline (T90 ≈ 190 °C) forms a film and is badly
     distributed in the long, cold manifold.
   - **Blow-by.** Leakage grows with the time per cycle, so it costs most at low speed (the same omission behind
     the compression residual).
   - **Period practice.** Ford's tester probably retarded the lever to stop the detonation Ricardo describes. The
     light-knock column above is within +16 % at 500 rpm and +7 % at 800 rpm.

   The club engineers' simulation of the stock engine (Sigworth 1999) showed the same pattern: +40 % at 500 rpm and
   within 1 % at 1500–2000 rpm, which they attributed to low-speed transmission and flywheel losses. As a result the
   model's torque peaks at 700 rpm, not 900.
2. **Upton's 600 rpm suction** is 1.61 inHg in the model against 1.10. That point stays ≥ 0.4 inHg high at any C_D
   up to 0.7; the difference is small in absolute terms.
3. **Motored compression is +13 % above Ford's 60 psig** (see above).
