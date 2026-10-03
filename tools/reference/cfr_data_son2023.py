"""Supercharged Octane Number (SAE 2023-01-0251, ANL; OSTI 1969816): knock-limited ('critical') compression ratio of
PRF100/96/93, a TSF and six gasolines vs intake pressure 1.0-1.5 bar at a fixed MAPO threshold on the CFR engine.

Source: tools/reference/cfr_sources/osti_1969816_son2023.pdf, Fig. 6 (p. 6) - born-digital Excel chart ->
VECTOR-PATH DIGITISATION: marker = centre of each circle/asterisk marker path; axes from the tick-label centres
(y: CR 5.5-8.0, x: intake pressure 1.0-1.5 bar); series identified by fill colour and legend order; the two black-
circle series (PRF100, PRF96) are separated by drawing order and checked against the solid (PRF100) and dashed (PRF96)
linear trendlines.

Test conditions (Sec. 'Supercharge Octane Number Test Method', Table 4): 600 rpm, intake air 52 C, spark 13 deg bTDC,
stoichiometric (not peak-knock) mixture, compressed dry air at 1.0-1.5 bar abs, the compression ratio raised until the
300-cycle mean MAPO (4-18 kHz band-pass, rectified, spark-plug transducer) reaches 0.6 bar; all cycles knock.
Pressure: bar ABSOLUTE ('all intake pressures ... in bar absolute').

Writes test/fixtures/cfr_son2023_critical_cr.json.
"""
from __future__ import annotations

import numpy as np

from cfr_data_common import (
    BAR, T0, axis_from_labels, colour_eq, conditions, drawing_points, page_drawings, page_words, prf, source_pdf,
    write_fixture,
)

PDF = source_pdf("osti_1969816_son2023.pdf")
CITATION = ("'Development of a Supercharged Octane Number and a Supercharged Octane Index', SAE Technical Paper 2023-01-0251 "
            "(Argonne National Laboratory), accepted manuscript OSTI 1969816")
URL = "https://www.osti.gov/servlets/purl/1969816"
FUELS = {  # Table 5 (p. 5)
    "TSF96.9": {"RON": 96.9, "MON": 85.2, "composition": "isooctane 5 / n-heptane 21 / toluene 74 %v"},
    "RON98Aro": {"RON": 97.9, "MON": 87.3, "composition": "n-par 8, iso-par 38, arom 40, naph 8, olef 5 %v"},
    "RON98E30": {"RON": 97.6, "MON": 87.6, "composition": "n-par 13, iso-par 28, arom 14, naph 7, olef 6, ethanol 30 %v"},
    "RON98Ole": {"RON": 98.5, "MON": 87.8, "composition": "n-par 12, iso-par 44, arom 13, naph 3, olef 27 %v"},
    "RON98CA": {"RON": 97.6, "MON": 86.6, "composition": "n-par 8, iso-par 32, arom 33, naph 24, olef 2 %v"},
    "RON98Alk": {"RON": 97.9, "MON": 96.7, "composition": "n-par 3, iso-par 96, arom 1 %v"},
    "Tier3 EEE": {"RON": 92.1, "MON": 84.0, "composition": "saturates 63, aromatics 23, olefins 6, ethanol 10 %v"},
    "PRF100": {"RON": 100, "MON": 100, "composition": "iso-octane"},
    "PRF96": {"RON": 96, "MON": 96, "composition": "96 %v iso-octane / 4 %v n-heptane"},
    "PRF93": {"RON": 93, "MON": 93, "composition": "93 %v iso-octane / 7 %v n-heptane"},
}
COLOURS = {"TSF96.9": (0.969, 0.216, 0.753), "RON98Aro": (0.439, 0.188, 0.627), "RON98E30": (0.498, 0.498, 0.498),
           "RON98Ole": (0.18, 0.459, 0.714), "RON98CA": (0.929, 0.49, 0.192), "RON98Alk": (0.0, 0.69, 0.314),
           "Tier3 EEE": (0.753, 0.0, 0.0)}


def main() -> None:
    drs = page_drawings(PDF, 6)
    words = page_words(PDF, 6)
    xax = axis_from_labels(words, (355, 295, 575, 307), "x")
    yax = axis_from_labels(words, (340, 105, 362, 297), "y")
    plot = (366, 114, 567, 292)

    def in_plot(d):
        r = d["rect"]
        return plot[0] <= r.x0 and r.x1 <= plot[2] + 3 and plot[1] <= r.y0 and r.y1 <= plot[3] and r.width < 6

    legend_markers = lambda d: 495 < d["rect"].x0 < 510 and d["rect"].y1 < 195  # noqa: E731
    marks = [(i, d) for i, d in enumerate(drs) if in_plot(d) and len(d["items"]) in (3, 4) and not legend_markers(d)]
    series: dict[str, list] = {k: [] for k in FUELS}
    black = []
    for i, d in marks:
        r = d["rect"]
        c = ((r.x0 + r.x1) / 2, (r.y0 + r.y1) / 2)
        if d["type"] == "s" and len(d["items"]) == 3 and colour_eq(d.get("color"), (0, 0, 0)):
            series["PRF93"].append(c)  # asterisk markers
            continue
        if d["type"] != "fs":
            continue
        for name, col in COLOURS.items():
            if colour_eq(d.get("fill"), col, 0.01):
                series[name].append(c)
                break
        else:
            if colour_eq(d.get("fill"), (0, 0, 0)):
                black.append((i, c))
    # two black-circle series: first drawn group = PRF100, second = PRF96 (verified against the trendlines below)
    black.sort()
    series["PRF100"] = [c for _, c in black[:6]]
    series["PRF96"] = [c for _, c in black[6:12]]
    # trendlines: black solid (PRF100), black dashed [4.95 6.61] (PRF96), black dotted [0 2.48] (PRF93)
    trend = {}
    for d in drs:
        if len(d["items"]) == 1 and d["type"] == "s" and d["rect"].width > 150:
            p = drawing_points(d)
            x, y = xax(p[:, 0]), yax(p[:, 1])
            dash = d.get("dashes") or ""
            col = d.get("color")
            if colour_eq(col, (0, 0, 0)):
                key = "PRF100" if dash.startswith("[] ") else ("PRF96" if "4.95" in dash else "PRF93")
            else:
                key = next((n for n, c in COLOURS.items() if colour_eq(col, c, 0.02)), None)
                if key is None:  # two trendlines use colours that match no marker series (ambiguous): not exported
                    continue
            trend[key] = {"intakePressurePa": (x * BAR).tolist(), "compressionRatio": y.tolist()}
    out_series = {}
    for name, pts in series.items():
        pts = sorted(pts)
        x = xax(np.array([p[0] for p in pts]))
        y = yax(np.array([p[1] for p in pts]))
        e = {"intakePressurePa": x * BAR, "compressionRatio": y, "fuel": dict(FUELS[name], label=name)}
        if name in trend:
            t = trend[name]
            tx = np.array(t["intakePressurePa"]) / BAR
            ty = np.array(t["compressionRatio"])
            e["linearTrendline"] = t
            e["maxDeviationFromTrendline"] = float(np.max(np.abs(np.interp(x, tx, ty) - y))) if len(x) else None
        out_series[name] = e
    assert len(out_series["PRF100"]["compressionRatio"]) == 6 and out_series["PRF100"]["compressionRatio"][0] > \
        out_series["PRF96"]["compressionRatio"][0]
    ds = {
        "key": "son2023_fig6_critical_cr_vs_intake_pressure",
        "kind": "knockLimitedCompressionRatio",
        "figure": "Fig. 6", "page": 6,
        "caption": "Critical compression ratio (0.6 bar MAPO) vs intake pressure, stoichiometric, spark 13 deg bTDC",
        "conditions": conditions(
            engine="ANL CFR F1 (SON campaign)", method="SON (critical CR at MAPO 0.6 bar)", rpm=600, sparkAdvanceDegBTDC=13,
            **{"lambda": 1.0}, intakeAirTemperatureK=52 + T0, coolantTemperatureK=100 + T0,
            compressionRatioSource="measured quantity (y axis); ANL CR calibration of this era (see README caveat)",
            notes="Intake pressure (x) is absolute and varies per point (first column at 1.013 bar = ambient-level supply). "
                  "Knock criterion: 300-cycle mean MAPO = 0.6 bar (4-18 kHz), NOT the ASTM knockmeter."),
        "extraction": {"method": "vector: marker path centres (circles 4.1 pt, asterisks), axes from tick-label centres",
                       "axisResidual": {"barPerPt": xax.per_unit(), "crPerPt": yax.per_unit(), "x": xax.resid, "y": yax.resid},
                       "uncertainty": {"compressionRatio": 0.002, "intakePressureBar": 0.0005},
                       "note": "Tier3 EEE has 3 points and PRF93 3 points (as plotted)."},
        "data": out_series,
    }
    out = {
        "description": "Critical (MAPO 0.6 bar) compression ratio vs intake pressure for PRF100/96/93 and seven other fuels "
                       "on the ANL CFR engine (SON 2023). Pa absolute.",
        "generator": "tools/reference/cfr_data_son2023.py",
        "source": {"key": "SON2023", "citation": CITATION, "url": URL, "osti": "1969816",
                   "localPdf": "tools/reference/cfr_sources/osti_1969816_son2023.pdf"},
        "datasets": [ds],
    }
    write_fixture("cfr_son2023_critical_cr.json", out)
    for k, v in out_series.items():
        print(k, [round(p / BAR, 3) for p in v["intakePressurePa"]], [round(c, 3) for c in v["compressionRatio"]],
              round(v.get("maxDeviationFromTrendline") or -1, 3))


if __name__ == "__main__":
    main()
