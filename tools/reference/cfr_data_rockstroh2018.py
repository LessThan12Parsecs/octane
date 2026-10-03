"""Rockstroh, Kolodziej, Jespersen, Goldsborough & Wallner 2018, 'Insights into Engine Knock: Comparison of Knock
Metrics across Ranges of Intake Temperature and Pressure in the CFR Engine', SAE Int. J. Fuels Lubr. 11(4):545-561,
SAE 2018-01-0210 (peer-reviewed author version, DTU Orbit).

PRF90 on the ANL CFR engine, 600 rpm, spark 13 deg bTDC, lambda = 1 (sweeps) / 0.89 (ASTM standard knock), compression
ratio sweeps 5.1-7.1 at intake-mixture temperatures 33 / 90 / 150 C (MON-style mixture heater) and intake pressures
1.0 / 1.1 / 1.2 / 1.28 bar abs; knockmeter (D-1 + 501-C) and Kistler 6045AU20 (200 kHz) cylinder pressure.

Source file: tools/reference/cfr_sources/rockstroh2018_sae2018-01-0210_dtu.pdf (DTU Orbit terms: private study copy -
do not redistribute; re-download from the URL below). Figures are embedded RASTER images (matplotlib):
  * Fig. 5 (xref 117, 770x1599 px): (a) knockmeter KU, (b) pressure at knock point, (c) dp/dtheta after KP,
    (d) peak (3 kHz low-pass-filtered) pressure, each vs CR, 12 series (colour = T_in, marker = P_in) + the ASTM
    standard-knock point (black diamond).
  * Fig. 6 (xref 125, 525x1090 px): (a) pressure at knock onset, (b) MAPO (6-20 kHz band-pass... per Eq. 1 / Fig. 3 3-40
    kHz), vs CR.
Digitisation: saturated colour masks (blue/orange/red), binary opening with a 5-px disc (removes lines and the pale
error bars), connected components; marker shape (circle / left-triangle / down-triangle / square) from the fill ratio
and the width/height profile slopes of the unopened blob; overlapping markers are NOT separated (merged blobs are
dropped). Fig. 6 markers are too small for shape classification: a Fig. 6 marker is kept only if exactly one same-colour
series has a Fig. 5 point at that compression ratio (+-0.025). Axes: Fig. 5 from the detected grid lines; Fig. 6 from the panel
frames (x 5.0-7.5; MAPO 0-4 bar; P_KO 15-55 bar). Uncertainty ~1 px: CR 0.004, KU 0.35, 0.12 bar (P_KP), 0.005
bar/CAD, 0.1 bar (peak p), 0.17 bar (P_KO), 0.017 bar (MAPO).

Writes test/fixtures/cfr_rockstroh2018_knock_vs_cr.json.
"""
from __future__ import annotations

import collections

import numpy as np
from scipy import ndimage

from cfr_data_common import BAR, T0, Axis, blobs, conditions, find_lines, load_image, prf, source_pdf, write_fixture

URL = "https://backend.orbit.dtu.dk/ws/files/236274325/KM_Investigation_CFR_Engine_SAE_Journal_Final.pdf"
CITATION = ("Rockstroh T., Kolodziej C.P., Jespersen M.C., Goldsborough S.S., Wallner T., 'Insights into Engine Knock: "
            "Comparison of Knock Metrics across Ranges of Intake Temperature and Pressure in the CFR Engine', SAE Int. J. "
            "Fuels Lubr. 11(4):545-561, 2018, doi:10.4271/2018-01-0210 (peer-reviewed author version, DTU Orbit)")
TIN = {"blue": 33.0, "orange": 90.0, "red": 150.0}
PIN = {"ci": 1.0, "lt": 1.1, "dn": 1.2, "sq": 1.28}
SE = (lambda r: (np.add.outer(np.arange(-r, r + 1) ** 2, np.arange(-r, r + 1) ** 2) <= r * r))(2)


def colour_masks(img):
    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    return {"blue": (B > 190) & (R < 80) & (G < 80), "red": (R > 190) & (G < 80) & (B < 80),
            "orange": (R > 190) & (G > 120) & (G < 200) & (B < 80)}


def shape_features(m, b):
    y0, x0 = max(b["y0"] - 2, 0), max(b["x0"] - 2, 0)
    sub = m[y0:b["y0"] + b["h"] + 2, x0:b["x0"] + b["w"] + 2]
    lab, _ = ndimage.label(sub, structure=np.ones((3, 3)))
    cy, cx = int(round(b["cy"] - y0)), int(round(b["cx"] - x0))
    lbl = lab[cy, cx] if lab[cy, cx] else int(np.bincount(lab.ravel())[1:].argmax()) + 1
    s = lab == lbl
    ys, xs = np.nonzero(s)
    s = s[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    h, w = s.shape
    wr, hc = s.sum(axis=1).astype(float), s.sum(axis=0).astype(float)
    sr = np.polyfit(np.linspace(-0.5, 0.5, h), wr, 1)[0] / w if h > 2 else 0.0
    sc = np.polyfit(np.linspace(-0.5, 0.5, w), hc, 1)[0] / h if w > 2 else 0.0
    return h, w, float(s.sum() / (h * w)), float(sr), float(sc)


def classify(f):
    h, w, fill, sr, sc = f
    if max(h, w) > 13:
        return None  # merged markers
    if fill > 0.9:
        return "sq"
    if sr < -0.45:
        return "dn"
    if sc > 0.45:
        return "lt"
    if fill > 0.6:
        return "ci"
    return None


def panel_markers(img, rows, cols, classify_shapes=True, min_area=12):
    out = []
    for colour, m in colour_masks(img).items():
        m = m.copy()
        m[: rows[0], :] = False
        m[rows[1]:, :] = False
        m[:, : cols[0]] = False
        m[:, cols[1]:] = False
        op = ndimage.binary_opening(m, structure=SE)
        for b in blobs(op, min_area=min_area):
            shp = classify(shape_features(m, b)) if classify_shapes else "?"
            if classify_shapes and shp is None:
                continue
            if not classify_shapes and max(b["w"], b["h"]) > 11:
                continue
            out.append({"colour": colour, "shape": shp, "x": b["cx"], "y": b["cy"]})
    return out


def standard_point(img, rows, cols):
    k = img.max(axis=2) < 50
    k[: rows[0], :] = False
    k[rows[1]:, :] = False
    k[:, : cols[0] + 3] = False
    k[:, cols[1] - 3:] = False
    op = ndimage.binary_opening(k, structure=SE)
    bl = [b for b in blobs(op, min_area=20) if b["w"] <= 16 and b["h"] <= 18]
    return (bl[0]["cx"], bl[0]["cy"]) if bl else None


def snap(detected, approx):
    return [min(detected, key=lambda d: abs(d - a)) for a in approx]


def main() -> None:
    try:
        pdf = source_pdf("rockstroh2018_sae2018-01-0210_dtu.pdf")
    except FileNotFoundError as e:  # not redistributable: keep the committed fixture
        print(f"skip: {e}")
        return
    img5 = load_image(pdf, 117).astype(int)
    g = (img5.max(axis=2) - img5.min(axis=2) <= 6) & (img5.min(axis=2) >= 200) & (img5.max(axis=2) <= 240)
    hl, vl = find_lines(g, 0, 0.3), find_lines(g, 1, 0.2)
    xax = Axis.fit(snap(vl, [98, 228, 358, 489, 619, 750]), [5.0, 5.5, 6.0, 6.5, 7.0, 7.5])
    panels5 = {
        "knockUnits": ((68, 412), Axis.fit(snap(hl, [67, 124, 182, 240, 298, 355, 413]), [120, 100, 80, 60, 40, 20, 0]), 1.0),
        "pressureAtKnockPointPa": ((416, 788), Axis.fit(snap(hl, [487, 530, 574, 617, 660, 704, 747, 790]),
                                                        [50, 45, 40, 35, 30, 25, 20, 15]), BAR),
        "pressureRiseRateAfterKPPaPerDeg": ((793, 1165), Axis.fit(snap(hl, [850, 890, 929, 969, 1048, 1088, 1167]),
                                                                  [1.6, 1.4, 1.2, 1.0, 0.6, 0.4, 0.0]), BAR),
        "peakFilteredPressurePa": ((1170, 1542), Axis.fit(snap(hl, [1296, 1346, 1395, 1445, 1495, 1544]),
                                                          [45, 40, 35, 30, 25, 20]), BAR),
    }
    cols5 = (101, 748)
    fig5: dict = {}
    allpts = []
    for q, (rows, yax, scale) in panels5.items():
        pts = panel_markers(img5, rows, cols5)
        for p in pts:
            p.update(cr=float(xax(p["x"])), value=float(yax(p["y"])) * scale, quantity=q)
        allpts += pts
        sp = standard_point(img5, rows, cols5)
        fig5[q] = {"points": pts, "standard": {"cr": float(xax(sp[0])), "value": float(yax(sp[1])) * scale} if sp else None}
    # series from each point's own marker shape; per (series, quantity) points closer than 0.02 CR are conflicting
    # classifications (a series has one point per CR) -> all of them dropped
    series: dict = collections.defaultdict(lambda: collections.defaultdict(list))
    conflicts = 0
    lookup: dict = collections.defaultdict(list)  # (colour, shape) -> CR values seen in Fig. 5
    for p in allpts:
        name = f"Tin{TIN[p['colour']]:.0f}C_Pin{PIN[p['shape']]:.2f}bar"
        series[name][p["quantity"]].append((p["cr"], p["value"]))
        lookup[(p["colour"], p["shape"])].append(p["cr"])
    for name in series:
        for q in list(series[name]):
            v = sorted(series[name][q])
            keep = [v[i] for i in range(len(v)) if (i == 0 or v[i][0] - v[i - 1][0] > 0.02)
                    and (i == len(v) - 1 or v[i + 1][0] - v[i][0] > 0.02)]
            conflicts += len(v) - len(keep)
            series[name][q] = keep
    # Fig. 6: small markers; series = the unique (colour, shape) whose Fig. 5 points include this CR (+-0.025)
    img6 = load_image(pdf, 125).astype(int)
    x6 = Axis.fit([65, 512], [5.0, 7.5])
    panels6 = {"pressureAtKnockOnsetPa": ((48, 279), Axis.fit([281, 46], [15.0, 55.0]), BAR),
               "mapoPa": ((305, 536), Axis.fit([538, 302], [0.0, 4.0]), BAR)}
    unmatched = 0
    std6 = {}
    for q, (rows, yax, scale) in panels6.items():
        for p in panel_markers(img6, rows, (68, 509), classify_shapes=False, min_area=8):
            cr = float(x6(p["x"]))
            cands = [k for k, crs in lookup.items() if k[0] == p["colour"] and min(abs(c - cr) for c in crs) < 0.025]
            if len(cands) != 1:
                unmatched += 1
                continue
            crs = lookup[cands[0]]
            name = f"Tin{TIN[p['colour']]:.0f}C_Pin{PIN[cands[0][1]]:.2f}bar"
            series[name][q].append((float(min(crs, key=lambda c: abs(c - cr))), float(yax(p["y"])) * scale))
        sp = standard_point(img6, rows, (68, 509))
        std6[q] = {"cr": float(x6(sp[0])), "value": float(yax(sp[1])) * scale} if sp else None
    out_series = {}
    for name, qs in sorted(series.items()):
        tin = float(name.split("_")[0][3:-1])
        pin = float(name.split("_")[1][3:-3])
        e = {"intakeMixtureTemperatureK": tin + T0, "intakePressurePa": pin * BAR}
        for q, v in qs.items():
            v = sorted(v)
            e[q] = {"compressionRatio": [a for a, _ in v], "value": [b for _, b in v]}
        out_series[name] = e
    standard = {q: v["standard"] for q, v in fig5.items()}
    standard.update(std6)
    ds = {
        "key": "rockstroh2018_prf90_knock_metrics_vs_cr",
        "kind": "knockOnsetAndIntensityVsCompressionRatio",
        "figure": "Figs. 5-6", "pages": [8, 9],
        "caption": "Knock metrics vs compression ratio for PRF90 at 12 intake conditions (3 temperatures x 4 pressures), "
                   "lambda = 1, spark 13 deg bTDC, 600 rpm; ASTM standard-knock point (lambda 0.89, CR 6.5) separately",
        "conditions": conditions(
            engine="ANL CFR F1/F2 (RON set-up + MON mixture heater, boosted air)", method="CR sweeps (RON-like)",
            fuel=prf(90), rpm=600, sparkAdvanceDegBTDC=13, **{"lambda": 1.0}, coolantTemperatureK=100 + T0,
            compressionRatioSource="x axis; engine CR verified by SAE30-oil hydraulic calibration (paper p. 4)",
            notes="intakeMixtureTemperatureK / intakePressurePa per series below. Intake mixture temperature is set with "
                  "the MON mixture heater. Pressure at knock point: KP from 3 kHz low-pass + second derivative (Kistler, "
                  "200 kHz); knock onset: 6-20 kHz band-pass exceeding 2 SD of the pre-KP noise; MAPO per Eq. 1 "
                  "(6-20 kHz band-pass in the text, 3-40 kHz in the Fig. 3 caption); 300 cycles per point."),
        "standardKnockPoint": {
            "description": "ASTM D2699 standard knock for PRF90 on that day: barometer 987 mbar -> IAT 45.3 C, CR 6.51, "
                           "50 KU, peak-knock lambda 0.89 (text p. 4 and Fig. 1 caption); knock intensity ~0.6 bar (p. 7)",
            "barometricPressurePa": 98700.0, "intakeAirTemperatureK": 45.3 + T0, "compressionRatio": 6.51,
            "lambda": 0.89, "knockUnits": 50, "digitisedFromFigures": standard,
            "knockPointThreshold": "dp/dtheta after KP 3.11 bar/CAD (Swarts, Kistler) / 2.5 bar/CAD (AVL spark plug)",
            "knockPointAt": "~60 % of total apparent heat release (Fig. 2)"},
        "extraction": {
            "method": __doc__.split("Digitisation:")[1].split("Writes")[0].strip().replace("\n", " "),
            "recovered": {q: sum(len(v.get(q, {}).get("value", [])) for v in out_series.values())
                          for q in list(panels5) + list(panels6)},
            "expectedPerQuantity": 72, "fig5ConflictingShapeAssignmentsDropped": conflicts,
            "fig6AmbiguousOrUnmatchedDropped": unmatched,
            "note": "Only unambiguous markers are kept: overlapping markers (low CR, low KU) and Fig. 6 markers whose CR "
                    "is shared by several same-colour series are dropped. Points are published 300-cycle means.",
        },
        "data": out_series,
        "reported": {"dominantFrequenciesHz": [6000, 10000, 14000], "kistlerNaturalFrequencyHz": 80000,
                     "knockmeterLowPassHz": 6500,
                     "findings": "KU correlates poorly with MAPO/K400/KSE (Pearson ~0.45 with peak filtered pressure); "
                                 "KO ~ KP; boost raises P_KP, heating lowers it; KI up to 3 bar mean (8 bar single cycles) "
                                 "at 90 C NA before 50 KU is reached"},
    }
    out = {
        "description": "Rockstroh et al. 2018: PRF90 knock point / onset pressure, pressure-rise rate after KP, peak "
                       "pressure, knockmeter reading and MAPO vs compression ratio at 12 intake T/p conditions. SI units.",
        "generator": "tools/reference/cfr_data_rockstroh2018.py",
        "source": {"key": "Rockstroh2018", "citation": CITATION, "url": URL,
                   "localPdf": "tools/reference/cfr_sources/rockstroh2018_sae2018-01-0210_dtu.pdf",
                   "licenceNote": "DTU Orbit: one private copy; do not redistribute (URL only)"},
        "datasets": [ds],
    }
    write_fixture("cfr_rockstroh2018_knock_vs_cr.json", out)
    print("recovered", ds["extraction"]["recovered"], "unmatched fig6", unmatched)
    print("standard", {k: (round(v["cr"], 3), round(v["value"] / (BAR if "Pa" in k else 1), 3)) if v else None for k, v in standard.items()})
    for name, e in out_series.items():
        print(name, {q: len(v["value"]) for q, v in e.items() if isinstance(v, dict)})


if __name__ == "__main__":
    main()
