"""Index of the CFR F-1 experimental validation dataset: test/fixtures/cfr_validation_index.json.

Named 'zz_' so that tools/reference/generate_all.sh (alphabetical) runs it after every cfr_data_*.py generator; it only
reads the committed fixtures. For every dataset it records the file, key, kind, source, figure/table, extraction
class, a quality grade, the operating conditions (SI) and the number of data points, so that tests can iterate:

    const index = JSON.parse(readFileSync('test/fixtures/cfr_validation_index.json', 'utf8'));
    for (const e of index.datasets.filter((e) => e.kind === 'pressureTraceCycleAveraged' && e.quality <= 'B')) { ... }

Quality grades (see test/fixtures/cfr_validation_README.md):
  A  exact: table transcription, or vector-path data with a published cross-check reproduced (e.g. gIMEP)
  B  vector-exact data but an operating condition (usually CR) is not stated / inferred
  C  raster digitisation with calibrated axes, complete series, error <= 1-2 px
  D  raster digitisation with gaps / possible series mis-assignment; use for trends
"""
from __future__ import annotations

import json

from cfr_data_common import FIXTURES, write_fixture

QUALITY = {
    "astm_guide_tables": "A", "astm_barometric_compensation": "A", "astm_mon_spark_timing": "A",
    "astm_compression_pressure_checks": "A", "astm_basic_compression_pressure": "C", "astm_spread_characteristic": "C",
    "astm_knock_intensity_definition": "A",
    "choi2018_fig2_prf98_ron": "B",
    "choi2018_fig9_prf100_st_p5_sparkPlug": "A", "choi2018_fig9_prf100_st_p5_flushMount": "A",
    "choi2018_fig9_prf100_st_m13_sparkPlug": "A", "choi2018_fig9_prf100_st_m13_flushMount": "A",
    "choi2018_tables": "A",
    "pal2018_fig1_ku_mapo_vs_spark": "C", "pal2018_fig7_knocking_st_m13": "C", "pal2018_fig8_spectrum_peaks": "C",
    "pal2018_fig9_knock_points": "C", "pal2018_tables": "A",
    "kw2017_ron98_operating_metrics": "C", "kw2017_fig1_ku_vs_lambda": "D",
    "hoth2021a_fig2_prf100_detonation_meter_input": "B", "hoth2021a_fig6_ku_mapo_vs_lambda": "B",
    "hoth2021a_fig7_prf93_near_stoich": "B", "hoth2021b_fig1_prf100_mapo_vs_spark": "B",
    "hoth2021b_fig3_prf_spark_sweeps": "B", "hoth2021b_fig5_prf89_prf97_stoich_cycles": "B",
    "hoth2025_fig2_prf98_pkl_cycle": "B", "hoth2025_critical_cr": "A",
    "son2023_fig6_critical_cr_vs_intake_pressure": "A",
    "rockstroh2018_prf90_knock_metrics_vs_cr": "D",
    "ltr_fig6_prf90_hcci_cr_for_ca50": "C", "ltr_fig4_prf_hcci_transfer_function": "C", "ltr_fig10_prf90_hcci_lthr": "C",
    "kalvakala_fig4_prf90_bron": "C", "kalvakala_fig5_prf90_bmon": "D", "kalvakala_text": "A",
    "neste_motored_ron87": "D",
}

COND_KEYS = ["method", "compressionRatio", "digitalCounter", "rpm", "sparkAdvanceDegBTDC", "lambda",
             "intakeAirTemperatureK", "intakeMixtureTemperatureK", "intakePressurePa", "barometricPressurePa",
             "knockIntensityKU"]


def count_points(obj) -> int:  # noqa: C901
    """Number of data samples: in every dict the parallel numeric arrays count once (their maximum length), nested
    dicts are summed, numeric scalars count 1, lists of records count their length. Regression lines
    ('trendline'), the digitised average spectrum and '_'-prefixed bookkeeping entries are not counted."""
    if isinstance(obj, list):
        if obj and all(isinstance(t, dict) for t in obj):
            return len(obj)
        return len(obj) if obj and all(isinstance(t, (int, float)) or t is None for t in obj) else 0
    if not isinstance(obj, dict):
        return 1 if isinstance(obj, (int, float)) and not isinstance(obj, bool) else 0
    arrays, n = [], 0
    for k, v in obj.items():
        if k in ("trendline", "linearTrendline", "averageSpectrum", "fuel", "unit") or k.startswith("_"):
            continue
        if isinstance(v, list) and v and all(isinstance(t, (int, float)) or t is None for t in v):
            arrays.append(len(v))
        elif isinstance(v, (dict, list)):
            n += count_points(v)
    return n + (max(arrays) if arrays else 0)


def extraction_class(ds) -> str:
    e = ds.get("extraction", "")
    m = (e.get("method", "") if isinstance(e, dict) else str(e)).lower().lstrip()
    if m.startswith("vector"):
        return "vector"
    if m.startswith("raster") or "colour mask" in m or "raster" in m[:40]:
        return "raster"
    if "vector" in m[:60]:
        return "vector"
    return "table/text"


def summarise_conditions(c) -> dict:
    if not isinstance(c, dict):
        return {}
    if "method" not in c and all(isinstance(v, dict) for v in c.values()):  # e.g. {PKL: {...}, stoich: {...}}
        return {k: summarise_conditions(v) for k, v in c.items()}
    out = {k: c.get(k) for k in COND_KEYS if c.get(k) is not None}
    f = c.get("fuel")
    if isinstance(f, dict):
        out["fuel"] = f.get("label") or f.get("kind")
    return out


def main() -> None:
    entries = []
    files = sorted(p for p in FIXTURES.glob("cfr_*.json") if p.name != "cfr_validation_index.json")
    for p in files:
        d = json.loads(p.read_text())
        if "datasets" in d:
            src_default = (d.get("source") or {}).get("key")
            for ds in d["datasets"]:
                key = ds.get("key")
                if key == "kw2017_ron98_operating_metrics":  # bar-chart values per fuel (raster)
                    ds = dict(ds, extraction={"method": "raster bar tops"},
                              data=[v for f in ds["figures"] for q in f["values"].values()
                                    for k, v in q.items() if k in ("PRF98", "RON98Alk", "PRF71E30", "RON98E30")])
                entries.append({
                    "file": p.name, "key": key, "kind": ds.get("kind"),
                    "source": ds.get("source") or src_default,
                    "figure": ds.get("figure") or ds.get("tables"),
                    "extraction": extraction_class(ds),
                    "quality": QUALITY.get(key, "?"),
                    "conditions": summarise_conditions(ds.get("conditions")),
                    "nPoints": count_points(ds.get("data") if "data" in ds else
                                            {k: v for k, v in ds.items() if k not in ("conditions", "extraction", "fuels")}),
                    "path": f"datasets[key={key}]",
                })
        else:  # ASTM files: fixed structure
            if p.name == "cfr_astm_guide_tables.json":
                for k, sub, kind, q in [("astm_guide_tables", "RON_9_16 / MON", "guideTable", "A"),
                                        ("astm_barometric_compensation", "barometricCompensation", "guideTableCompensation", "A"),
                                        ("astm_mon_spark_timing", "monSparkTiming", "sparkTimingTable", "A")]:
                    n = 0
                    if k == "astm_guide_tables":
                        n = len(d["RON_9_16"]["onToCounter"]["octaneNumber"]) + sum(len(v["onToCounter"]["octaneNumber"]) for v in d["MON"].values())
                    elif k == "astm_barometric_compensation":
                        n = sum(len(d[sub][m][s]) for m in ("RON", "MON") for s in ("below2992", "above2992"))
                    else:
                        n = len(d[sub]["counter"])
                    entries.append({"file": p.name, "key": k, "kind": kind, "source": "ASTM", "figure": sub,
                                    "extraction": "table/text", "quality": q, "conditions": {"method": "RON/MON standard"},
                                    "nPoints": n, "path": sub})
            elif p.name == "cfr_astm_compression_pressure.json":
                entries.append({"file": p.name, "key": "astm_basic_compression_pressure", "kind": "motoredPeakPressure",
                                "source": "ASTM", "figure": "D2699 Fig. 2, D2700 Fig. 2", "extraction": "raster",
                                "quality": "C", "conditions": {"method": "motored, counter 930", "rpm": "600 (RON) / 900 (MON)"},
                                "nPoints": sum(len(e["barometerInHg"]) for e in d["basicCylinderHeightVsBarometer"]),
                                "path": "basicCylinderHeightVsBarometer"})
                entries.append({"file": p.name, "key": "astm_compression_pressure_checks", "kind": "motoredPeakPressure",
                                "source": "ASTM", "figure": "D2699/D2700 Table A2.2", "extraction": "table/text",
                                "quality": "A", "conditions": {"method": "motored"}, "nPoints": len(d["checkPoints"]),
                                "path": "checkPoints"})
            elif p.name == "cfr_astm_knock_intensity.json":
                entries.append({"file": p.name, "key": "astm_knock_intensity_definition", "kind": "knockMeterDefinition",
                                "source": "ASTM", "figure": "D2699 3.1.24-25, A2.4-A2.5", "extraction": "table/text",
                                "quality": "A", "conditions": {}, "nPoints": 0, "path": "standardKnockIntensity"})
                entries.append({"file": p.name, "key": "astm_spread_characteristic", "kind": "knockMeterSpread",
                                "source": "ASTM", "figure": "Fig. A2.6", "extraction": "raster", "quality": "C",
                                "conditions": {}, "nPoints": sum(len(s["octaneNumber"]) for s in d["spreadCharacteristic"]),
                                "path": "spreadCharacteristic"})
    ranking = [
        {"rank": 1, "key": "astm_guide_tables", "why": "the defining end-to-end target: CR (counter) at standard knock vs PRF ON, "
                                                       "RON and MON, exact; the simulator's knock-limited CR for PRF 60-100 must "
                                                       "track it (absolute CR subject to the counter->CR calibration)"},
        {"rank": 2, "key": "choi2018_fig9_prf100_st_p5_flushMount", "why": "full-cycle 0.1 deg cycle-averaged pressure, stated CR "
                                                                          "7.55, lambda 0.89, iso-octane, NON-knocking (spark +5): "
                                                                          "compression, flame propagation, expansion, gas exchange; "
                                                                          "gIMEP cross-check reproduced to 0.1 %"},
        {"rank": 3, "key": "choi2018_fig9_prf100_st_m13_flushMount", "why": "same engine state at standard knock (spark -13); "
                                                                           "burn rate + knock onset; spark-plug-sensor twin gives "
                                                                           "the transducer uncertainty"},
        {"rank": 4, "key": "astm_compression_pressure_checks", "why": "motored peak pressure at known cylinder heights (RON "
                                                                     "and MON): heat loss / blow-by / CR calibration check"},
        {"rank": 5, "key": "choi2018_fig2_prf98_ron", "why": "full-cycle cylinder + intake-port + exhaust-port pressure: "
                                                            "gas-exchange boundary conditions and residual"},
        {"rank": 6, "key": "pal2018_fig1_ku_mapo_vs_spark", "why": "knockmeter and MAPO vs spark at fixed CR: knock onset "
                                                                  "threshold and KI growth"},
        {"rank": 7, "key": "son2023_fig6_critical_cr_vs_intake_pressure", "why": "knock-limited CR vs PRF ON and boost at a "
                                                                                "pressure-based KI threshold"},
        {"rank": 8, "key": "hoth2021a_fig2_prf100_detonation_meter_input", "why": "the 501-C input signal vs the cylinder "
                                                                                 "pressure of the same cycle: identify the "
                                                                                 "knockmeter transfer function"},
        {"rank": 9, "key": "hoth2021b_fig3_prf_spark_sweeps", "why": "KU, MAPO, gIMEP vs spark for PRF87-100 at PKL and "
                                                                    "stoichiometry"},
        {"rank": 10, "key": "ltr_fig6_prf90_hcci_cr_for_ca50", "why": "autoignition chemistry + compression heating without a "
                                                                      "flame (HCCI, lambda 3)"},
    ]
    out = {
        "description": "Index of the CFR F-1 experimental validation fixtures (test/fixtures/cfr_*.json). Units SI; crank "
                       "angle deg with 0 = firing TDC in [-360, 360); pressures Pa absolute. See cfr_validation_README.md.",
        "generator": "tools/reference/cfr_data_zz_index.py",
        "qualityGrades": {"A": "exact (table or vector with reproduced published cross-check)",
                          "B": "vector-exact, but a condition (usually CR) not stated/inferred",
                          "C": "raster, calibrated, complete", "D": "raster with gaps / possible series mis-assignment"},
        "datasets": entries,
        "ranking": ranking,
        "totals": {"files": len(files), "datasets": len(entries), "points": sum(e["nPoints"] for e in entries)},
    }
    write_fixture("cfr_validation_index.json", out)
    for e in entries:
        print(f'{e["quality"]} {e["extraction"]:10s} {e["nPoints"]:6d} {e["file"]:38s} {e["key"]}')
    print(out["totals"])


if __name__ == "__main__":
    main()
