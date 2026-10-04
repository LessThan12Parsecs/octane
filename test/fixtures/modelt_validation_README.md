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
