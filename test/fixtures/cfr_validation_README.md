# CFR F-1 experimental validation dataset

Experimental (and standard-method) data for validating the Waukesha CFR F-1 simulator against the real engine.
Every number is traceable to a document, table/figure and page; every digitised curve states its method and
uncertainty. Machine-readable entry point: `cfr_validation_index.json` (one entry per dataset with file, key,
kind, source, figure, extraction class, quality grade, operating conditions, number of points, plus a ranked
target list).

## Conventions (all `cfr_*.json`)

| quantity | convention |
|---|---|
| units | SI (Pa, K, m, s, kg; `*Deg` = crank degrees). Source units are converted in the generator; the source unit is stated in each file's `conventions`. |
| crank angle `thetaDeg` | 0 = **firing TDC**, cycle [-360, 360). Every source used "CAD aTDC" (firing) except the Neste motored trace (0-720, TDC = 360; converted with θ = CAD − 360). |
| pressure `pPa` | **absolute**. All in-cylinder sources were absolute already (intake-stroke level 0.9-1.1 bar, Kulite/AVL absolute or pegged by the authors); no offset was applied. ASTM compression pressures are **gauge** (psig) in the source and are given both as gauge and as absolute (= gauge + barometer). |
| spark | `sparkAdvanceDegBTDC` > 0 before TDC (= −spark timing in "CAD aTDC"), as `OperatingPoint.sparkAdvanceDeg`. |
| conditions | `conditions` objects use the `cfr_data_common.conditions()` keys (`method, fuel, compressionRatio, compressionRatioSource, digitalCounter, rpm, sparkAdvanceDegBTDC, lambda, phi, intakeAirTemperatureK, intakeMixtureTemperatureK, intakePressurePa, barometricPressurePa, coolantTemperatureK, oilTemperatureK, knockIntensityKU, notes`); `null` = not reported. |
| counter → CR | `compressionRatioRigidRaise` = `cfrCompressionRatioAtCounter` (src/physics/engines/cfr.ts: rigid cylinder raise, 0.0007 in/digit, h₉₃₀ fitted to Choi 2018 Fig. 4); `compressionRatioChoiPolynomial` = `cfrCompressionRatioFromCounter` (Choi 2018 cubic, 400-1400). Replicated in `tools/reference/cfr_data_common.py`; cross-checked against the TypeScript functions with a one-off vitest script (not committed): max \|ΔCR\| = 5.0e-6 over all 3 200+ guide-table rows (fixture rounding), MON spark table vs `cfrMonSparkAdvanceDeg` ≤ 0.023°. |
| MAPO | max of the band-passed, rectified pressure per cycle, 300-cycle mean. **Definitions differ**: Hoth/SON/critical-CR 4-18 kHz; Rockstroh 6-20 kHz (text) / 3-40 kHz (Fig. 3 caption); Pal 2018 = AVL IndiCom default. |

Quality grades: **A** exact (table transcription, or vector data with a reproduced published cross-check);
**B** vector-exact data but an operating condition (usually CR) not stated / inferred; **C** raster, calibrated,
complete (≤ 1-2 px); **D** raster with gaps or possible series mis-assignment (trends only).

## Inventory (11 files, 37 datasets, ≈ 81 000 samples)

| file | source | content | datasets / samples | grade |
|---|---|---|---|---|
| `cfr_astm_guide_tables.json` | ASTM D2699-15a, D2700-14 (archive.org public copies) | full guide tables at standard KI: RON 40.0-120.3 → counter (804 rows) + inverse (counter 450-1249) + dial (804); MON 9/16, 19/32, 3/4 in venturi (counter, inverse, dial; 3×804); barometric compensation RON 21.0-30.9 inHg incl. required IAT (100 rows), MON 22.0-30.9 (90 rows); MON spark-timing table (13 rows); CR by both relations; standard operating conditions with tolerances | 3 / 3 407 | A |
| `cfr_astm_compression_pressure.json` | D2699 Fig. 2 + Table A2.2; D2700 Fig. 2 + Table A2.2 | motored peak (gauge) pressure defining the basic cylinder height (counter 930) vs barometer, RON and MON (3 venturis) — straight lines, e.g. RON 202.2 psig, MON(9/16) 176.0 psig at 29.92 inHg; check points RON 93.4 ON (778 counts) 169 ± 2 psig, 105 ON (1061) 241 ± 4; MON 81.1 (578) 120 ± 2, 105 (1008) 194 ± 4 | 2 / 220 | A (table), C (figures, ±0.2 psi) |
| `cfr_astm_knock_intensity.json` | D2699/D2700 §3.1, 10.3.20-22, 11.3, A2.4-A2.5; Rockstroh 2018; Hoth 2021a | definition of standard knock intensity (analog meter 50 ± 2 divisions after adjusting METER READING on the guide-table PRF at max-knock F/A; digital ≈0.15 V p-p RON / 0.25 V MON), D-1 magnetostrictive pickup (∝ dp/dt), 501-C chain (input filter → threshold → stretch → multi-cycle integration), input low-pass ≈ 6.5 kHz (Swarts via Rockstroh), time-constant/spread defaults, spread characteristic Fig. A2.6 (RON & MON, 80.6-105 ON) | 2 / 488 | A / C |
| `cfr_choi2018_traces.json` | Choi et al. 2018, SAE 2018-01-0848 (OSTI 1501884) | **vector** 300-cycle-average traces at 0.1°, full cycle: PRF98 RON test cylinder + intake-port + exhaust-port pressure; PRF100 CR 7.55 λ 0.89 at spark −13 (standard knock) and +5 (non-knocking), each with spark-plug AND flush-mount transducer; Tables 2-6 | 6 / 43 205 | A (B for PRF98: CR not stated) |
| `cfr_pal2018_knock.json` | Pal et al. 2018, SAE 2018-01-0187 (OSTI 1572720) | KU and MAPO vs spark timing (21 points, −14 … +6), knocking cycle 169, resonance peaks of the average spectrum (6.46, 10.83, 14.64 kHz), knock points of 219 of 300 cycles (mean 11.03°, SD 0.97°), Table 3 (TPA wall T, T_IVC, trapped mass 0.628 g, RGF 6.04 %), Table 5 | 5 / 644 | C (tables A) |
| `cfr_kw2017_ron98.json` | Kolodziej & Wallner 2017 (OSTI 1394801) | 4 RON-98 fuels at standard RON: λ, fuel rate, mixture T and carburettor ΔT (→ upstream IAT), gIMEP, ITE, exhaust T, CA10/CA50/knock point, KU vs PPRR and knock-pressure peak; KU vs λ markers | 2 / 89 | C / D |
| `cfr_hoth_knock_metrics.json` | Hoth & Kolodziej 2021 Part 1 (OSTI 1880351) and Part 2 (OSTI 2423200); ANL 2025, SAE 2025-01-8451 (OSTI 2561394) | **vector**: PRF100 knocking cycle + **501-C input signal** (V) of the same cycle, full cycle 0.1°; KU and MAPO vs λ for PRF93/95/97 at the RON-95 CR; PRF93 near-stoich λ sensitivity; PRF100 MAPO vs spark; KU/MAPO/gIMEP vs spark for PRF87-100 at peak-knock λ and stoichiometric (regression lines + markers); representative PRF89/PRF97 stoich cycles (full cycle 0.1°) + RoHR; PRF98 representative knocking cycle + band-passed signal; critical CR at 0.6 bar MAPO for 12 fuels (Table 7) | 8 / 30 208 | B (critical CR A) |
| `cfr_son2023_critical_cr.json` | SON, SAE 2023-01-0251 (OSTI 1969816) | **vector**: critical CR (MAPO 0.6 bar, λ 1, spark 13° bTDC, 52 °C) vs intake pressure 1.013-1.5 bar for PRF100/96/93 and 7 other fuels | 1 / 54 | A |
| `cfr_rockstroh2018_knock_vs_cr.json` | Rockstroh et al. 2018, SAE 2018-01-0210 (DTU Orbit author version) | PRF90 CR sweeps at 33/90/150 °C × 1.0/1.1/1.2/1.28 bar (λ 1): KU, p at knock point, dp/dθ after KP, peak pressure, p at knock onset, MAPO; ASTM standard-knock point (987 mbar, IAT 45.3 °C, CR 6.51, λ 0.89, 50 KU, MAPO ≈ 0.67 bar) | 1 / 229 of 432 | D (standard point C) |
| `cfr_hcci_autoignition.json` | OSTI 1583125 (ANL HCCI in the standard CFR); Kalvakala et al. (OSTI 1962089) | PRF90 HCCI (λ ≈ 3, 600 rpm): CR for CA50 = 3° aTDC vs intake T (36-180 °C) at 0.9/1.0/1.15/1.3 bar (21 pts); PRF70-97 transfer functions (16 pts); LTHR heat-release rates at 52 °C (3 CRs); 300-cycle pressure at BRON (1.3 bar, 52 °C, CR 12.84) and BMON (1.0 bar, 149 °C, CR 14.05) | 6 / 2 452 | C (BMON D) |
| `cfr_neste_motored.json` | Bhattacharya et al. (Aalto/ANL/Neste, OSTI 2536666) Fig. 5 | motored (hot, spark off) 250-cycle cylinder pressure at RON-relevant intake (600 rpm, 1 atm, 35 °C), dotted-line samples + zoomed peak (10.4 bar at ≈ −1°); **CR not stated** | 1 / 102 | D |

## Ranked validation targets

1. **ASTM guide tables** (`astm_guide_tables`, A) — knock-limited CR (counter) at standard KI vs PRF ON, RON and
   MON. This *is* the octane scale. Compare the simulated counter/CR at standard KI for PRF 60-100; RON 90 → 726
   (CR 6.43), RON 100 → 919 (7.49), MON 90 → 749 (6.54). Absolute CR inherits the counter→CR calibration (see
   caveats) — compare the **slope** dCR/dON first (e.g. RON 80→100: counter 609→919).
2. **Choi 2018 Fig. 9, PRF100 spark +5, flush-mount** (A) — full-cycle, 0.1°, CR 7.55, λ 0.89, iso-octane,
   non-knocking. The compression from IVC to +5° is effectively a **motored compression at a known CR** with
   fired-engine walls/residual (p(IVC ≈ −152°) = 1.01 bar, p(TDC) = 13.87 bar); then flame propagation (peak
   26.5 bar at 37.3°), expansion and gas exchange. Recomputed gIMEP 7.902 bar vs published 7.892 (+0.13 %).
3. **Choi 2018 Fig. 9, PRF100 spark −13 (standard knock)** (A) — same state, knocking: burn rate + knock onset
   (peak 47.9 bar at 13.0°, gIMEP 8.146 vs 8.137 published). Pair with Pal Fig. 1/9 (knock point 11.0 ± 0.97°,
   KU 47.5, MAPO 2.27 bar at the logged −12.72°) and Hoth 2021b Fig. 1 (same sweep, 4-18 kHz MAPO 1.43 bar).
   The spark-plug-transducer twins bound the sensor uncertainty (−1.2 % / −2.2 % gIMEP).
4. **ASTM motored compression pressure** (A) — peak motored gauge pressure at known cylinder heights, RON and
   MON: 169 ± 2 psig at counter 778 (RON), 241 ± 4 at 1061; MON 120 ± 2 at 578, 194 ± 4 at 1008; basic
   setting 202 psig (RON) / 176 psig (MON 9/16) at counter 930 and 29.92 inHg. Tests heat loss, blow-by,
   valve timing and the counter→CR relation simultaneously. NB the gauge is a check-valve peak-hold device in
   the pickup hole, carburettor drained (air only).
5. **Choi 2018 Fig. 2, PRF98 RON** (B) — cylinder + **intake-port + exhaust-port** pressure over the full cycle:
   gas-exchange boundary conditions (mean port pressures 0.974 / 0.989 bar), trapped mass/residual via the
   model; CR not stated (guide-table 867 counts → 7.18 before barometric compensation).
6. **Pal 2018 Fig. 1** (C) — KU and MAPO vs spark at fixed CR: the knockmeter falls to 0 by −8° while MAPO
   decays smoothly to 0.05 bar at +2°; a direct test of the knock-onset model and of any knockmeter model.
7. **SON 2023 Fig. 6 + critical-CR Table 7** (A) — knock-limited CR vs PRF ON (93/95/96/97/100) and intake
   pressure at a pressure-based KI threshold (0.6 bar MAPO, λ 1).
8. **Hoth 2021a Fig. 2** (B) — the 501-C input voltage and the cylinder pressure of the same knocking cycle:
   identify the knockmeter input transfer function (the input is a smooth, pressure-like low-pass signal
   peaking at 13.0° vs pressure 13.2°, with no trace of the 6.5 kHz ringing).
9. **Hoth 2021b Fig. 3 / Hoth 2021a Fig. 6** (B) — KU, MAPO and gIMEP vs spark and vs λ for PRF87-100.
10. **OSTI 1583125 Fig. 6 / Kalvakala BRON** (C) — flame-free autoignition (HCCI λ 3): end-gas-chemistry +
    compression-heating check at CR 11-15.7, intake 0.9-1.3 bar, 36-180 °C.

Secondary: Pal Fig. 8 resonance frequencies (6.46/10.83/14.64 kHz, cf. Rockstroh 6/10/14 kHz) for the
acoustic-mode model; KW17 fuel rate/λ/mixture T for the carburettor model (PRF98: λ 0.886, 0.747 kg/h,
T_mix 28.3 °C, 13.1 K carburettor drop → IAT 41.3 °C); Rockstroh sweeps (trends of knock-point pressure,
dp/dθ and MAPO with CR, T and p).

## Extraction methods and uncertainties

* **Table transcription** (ASTM, text tables): PyMuPDF word boxes grouped into rows, values assigned to the
  printed column by x position; completeness and anchor values asserted in the generator (e.g. RON 90 → 726,
  RON 100 → 919, MON 90 → 749; D2700 Table 3 rows checked against the text layer). Exact.
* **Vector-path digitisation** (Choi, Hoth ×3, SON): the Excel charts are stored as PDF drawing operators; data =
  path end points, axes = least-squares fits through the chart's own tick segments or tick-label centres
  (residual < 0.05 pt). Coordinate rounding 0.001 pt ⇒ < 0.01° and < 0.002 bar. Several traces are stored for the
  full cycle although the printed axes clip them (Choi Fig. 9, Hoth Figs. 2 and 5, 2025 Fig. 2) — the stored data
  are exported. Choi Fig. 4 (CR vs counter) was *not* used: its markers are a dense rendering of the cubic, not
  the oil measurements.
* **Raster digitisation** (ASTM Figs. 2/A2.6, Pal, KW17, Rockstroh, HCCI, Neste): colour masks on the embedded
  images at native resolution, axes from detected grid lines/tick marks; markers by morphological opening +
  centroid, bars by top edge, curves by column-wise run centres with continuity tracking. Pixel scales and
  uncertainties are stored per dataset (typically 1 px: 0.04° / 0.01 bar for Pal Fig. 1, 0.15 psi for ASTM Fig. 2).
  Known losses: overlapping markers (Rockstroh 229/432 points recovered, KW17 Fig. 1 3 of 6 PRF98 markers,
  HCCI Fig. 4 16/18), Pal Fig. 7 ±SD band not extractable (hidden under four cycle lines).

## Caveats (read before writing tests)

* **CR calibration differs between ANL papers.** Choi/Pal quote CR 7.55 for PRF100 at standard knock (rigid-raise
  relation: counter 919 → 7.49 at 29.92 inHg; 7.55 corresponds to ≈ +9 counts, i.e. the D2699 Table A4.4 compensation for a ≈ 29.6 inHg barometer ✓). The 2025 paper
  quotes CR **7.26** for the standard RON-95 cylinder height (relation: 805 counts → 6.82), and SON 2023 finds a
  critical CR 7.82 for PRF100 at stoichiometry / MAPO 0.6. The Neste engine uses CR = 6345/(1850 − DCR) + 1
  (0.35 higher at 930). Treat reported CRs from different campaigns as ±0.4 systematic; prefer datasets with a
  stated CR and the same campaign (Choi/Pal 2017-18).
* **Standard knock intensity is a relative, meter-defined state**: 50 KU is set by adjusting the meter on the
  guide-table PRF. A physics model needs a knockmeter model (input filter < 6.5 kHz, threshold, integration) or a
  calibrated equivalence (e.g. ASTM standard knock ≈ MAPO 0.6-0.7 bar at 4-18 kHz for PRF90-98, Rockstroh
  0.67 bar, Hoth Table 6 0.69 bar; Pal's IndiCom MAPO is ~1.6× larger).
* The published ASTM tables contain a few typos, transcribed as printed and listed in `sourceErrata` of `cfr_astm_guide_tables.json` (D2700 MON dial entries at 42.6, 77.5, 117.7 ON (9/16) and 61.1 ON (19/32); MON 19/32 counter 1052 at 115.7 ON after 1053). The RON and MON 9/16 forward and inverse tables agree to ≤ 0.1 ON.
* Spark-plug transducers (AVL GU13Z-24) read 1-2 % low during combustion vs flush-mount (Choi Fig. 9).
* Cycle-averaged traces smear knock ringing; the single representative cycles (Pal 169, Hoth 2021a Fig. 2, Hoth
  2025 Fig. 2, Hoth 2021b Fig. 5) keep it but are sampled at 0.1° (36 kHz at 600 rpm).
* Pal Table 3 trapped mass/RGF/T_IVC and Choi TPA temperatures are GT-Power **model** outputs, not measurements.
* Logged spark timing differs from the nominal step by +0.2…0.35° (Pal/Hoth Fig. 1): use the digitised value.
* The Pal/Hoth 2021b Fig. 1 spark sweeps are almost certainly the same PRF100 data (identical logged timings).

## Gaps (not found in open-access sources)

* **MON in-cylinder data**: no open trace/metric at D2700 conditions (only the ASTM MON tables, spark table and
  compression checks). MON validation is limited to the guide table and motored peak pressure.
* **Motored traces at a stated CR**: only the ASTM peak values and the pre-spark part of Choi ST +5; the Neste
  full motored trace has no CR.
* Measured residual fraction, trapped mass, volumetric efficiency and wall heat flux (only TPA-model values);
  Choi Fig. 7 measured fuel rate/gIMEP pairs exist but without per-point conditions (not extracted).
* COV of IMEP (not reported numerically), NOₓ (none for PRFs on the CFR), knock onset vs ON *at standard KI*
  other than via the guide table.
* Paywalled/unreachable: Swarts & Yates (SAE 2005-01-2081, 2007-01-0008; bouncing-pin model and 501-C circuit),
  Swarts & Kalaskar 2020, Hoth et al. SAE 2019-01-0627 and 2018-01-1672, Keskinen et al. SAE 2022-01-1082,
  Foong 2013 PhD thesis (Univ. Melbourne; repository host not reachable, SciSpace copy 403), Waqas et al. (KAUST).

## Reproduction

Sources (≤ 5 MB each) are in `tools/reference/cfr_sources/`: `astm_D2699-15a.pdf`
(https://archive.org/details/gov.law.astm.D2699.15A), `astm_D2700-14.pdf`
(https://archive.org/details/gov.law.astm.D2700.14), `osti_<id>_*.pdf` (https://www.osti.gov/servlets/purl/<id>
for 1501884, 1572720, 1394801, 1880351, 2423200, 2561394, 1969816, 1583125, 1962089, 2536666) and
`rockstroh2018_sae2018-01-0210_dtu.pdf` (DTU Orbit, private-study copy — git-ignored, not for redistribution;
https://backend.orbit.dtu.dk/ws/files/236274325/KM_Investigation_CFR_Engine_SAE_Journal_Final.pdf).

```
uv pip install --python .venv/bin/python pymupdf     # + numpy, scipy (already in the venv)
.venv/bin/python tools/reference/cfr_data_astm.py      # and cfr_data_{choi2018,hcci,hoth,kw2017,neste_motored,
                                                       #   pal2018,rockstroh2018,son2023}.py, then cfr_data_zz_index.py
```

`generate_all.sh` runs them in the right (alphabetical) order; `cfr_data_pal2018.py` reads the Choi fixture and
`cfr_data_zz_index.py` reads all of them. All generators are deterministic.
