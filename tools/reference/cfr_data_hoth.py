"""ANL (Hoth, Kolodziej et al.) CFR RON-rating studies with simultaneous knockmeter and cylinder-pressure knock metrics:

  [Hoth21a] Hoth A., Kolodziej C.P., 'Effects of knock intensity measurement technique and fuel chemical composition on
            the research octane number (RON) of FACE gasolines: Part 1 - Lambda and knock characterization', Fuel 2021
            (accepted manuscript OSTI 1880351).
  [Hoth21b] Hoth A., Kolodziej C.P., '... Part 2 - Effects of spark timing' (accepted manuscript OSTI 2423200).
  [Hoth25]  'Effects of Critical Compression Ratio on Rating Gasoline Knock Propensity', SAE 2025-01-8451 (ANL,
            accepted manuscript OSTI 2561394).

All three are born-digital (Excel charts exported as PDF drawing operators) -> VECTOR-PATH DIGITISATION. Axes are
calibrated on the centres of the numeric tick labels (text layer) or on the grid-line segments; markers are the
centres of their marker paths (or of the clusters of short segments that form 'x' markers); regression lines
('trendline' paths, third-order polynomials per the authors, R^2 > 0.99) are extracted as published.

MAPO in these papers: 4-18 kHz band-pass (Hoth21a/b, Hoth25) of the AVL GU13Z-24 spark-plug transducer signal sampled
at 0.1 CAD (36 kHz at 600 rpm), rectified, max per cycle in -60..+60 CAD, averaged over 300 cycles.
Crank angle: source 'deg aTDC' = firing-TDC convention (identical). Pressure: bar absolute -> Pa.

Writes test/fixtures/cfr_hoth_knock_metrics.json.
"""
from __future__ import annotations

import numpy as np

from cfr_data_common import (
    BAR, T0, Axis, axis_from_labels, colour_eq, conditions, cr_rigid_raise, drawing_points, legend_names,
    page_drawings, page_words, prf, short_segment_markers, source_pdf, style_key, write_fixture,
)

P1 = source_pdf("osti_1880351_hoth2021_part1.pdf")
P2 = source_pdf("osti_2423200_hoth_part2.pdf")
P3 = source_pdf("osti_2561394_critical_cr.pdf")
SRC = {
    "Hoth2021a": {"citation": "Hoth A., Kolodziej C.P., 'Effects of knock intensity measurement technique and fuel chemical "
                              "composition on the research octane number (RON) of FACE gasolines: Part 1 - Lambda and knock "
                              "characterization', Fuel (2021), doi:10.1016/j.fuel.2021.120??? (accepted manuscript OSTI 1880351)",
                  "url": "https://www.osti.gov/servlets/purl/1880351", "localPdf": "tools/reference/cfr_sources/osti_1880351_hoth2021_part1.pdf"},
    "Hoth2021b": {"citation": "Hoth A., Kolodziej C.P., '... Part 2 - Effects of spark timing' (accepted manuscript OSTI 2423200)",
                  "url": "https://www.osti.gov/servlets/purl/2423200", "localPdf": "tools/reference/cfr_sources/osti_2423200_hoth_part2.pdf"},
    "Hoth2025": {"citation": "'Effects of Critical Compression Ratio on Rating Gasoline Knock Propensity', SAE Technical Paper "
                             "2025-01-8451 (ANL; accepted manuscript OSTI 2561394)",
                 "url": "https://www.osti.gov/servlets/purl/2561394", "localPdf": "tools/reference/cfr_sources/osti_2561394_critical_cr.pdf"},
}
SRC["Hoth2021a"]["citation"] = SRC["Hoth2021a"]["citation"].replace("doi:10.1016/j.fuel.2021.120??? ", "")
RON95_COUNTER = 805  # D2699-15a Table A4.1, RON 95.0 at 29.92 in Hg (cfr_astm_guide_tables.json)
BLUE = (0.267, 0.447, 0.769)
GREEN = (0.439, 0.678, 0.278)
ORANGE = (0.929, 0.49, 0.192)


def series_by_style(drawings, min_items=100):
    return [d for d in drawings if d["type"] == "s" and len(d["items"]) >= min_items]


def in_box(d, box):
    r = d["rect"]
    cx, cy = 0.5 * (r.x0 + r.x1), 0.5 * (r.y0 + r.y1)
    return box[0] <= cx <= box[2] and box[1] <= cy <= box[3]


def curve(d, xax, yax, lo=None, hi=None):
    pts = drawing_points(d)
    x, y = xax(pts[:, 0]), yax(pts[:, 1])
    o = np.argsort(x)
    x, y = x[o], y[o]
    if lo is not None:
        s = (x >= lo) & (x <= hi)
        x, y = x[s], y[s]
    return x, y


def assign_to_curves(markers_xy: np.ndarray, curves: dict[str, tuple[np.ndarray, np.ndarray]]):
    """Assign each marker (data coords) to the nearest regression curve (vertical distance at the marker x)."""
    out = {k: [] for k in curves}
    worst = {k: 0.0 for k in curves}
    for mx, my in markers_xy:
        best, bd = None, np.inf
        for k, (x, y) in curves.items():
            if mx < x.min() - 0.3 or mx > x.max() + 0.3:
                continue
            dd = abs(np.interp(mx, x, y) - my)
            if dd < bd:
                best, bd = k, dd
        if best is not None:
            out[best].append((float(mx), float(my)))
            worst[best] = max(worst[best], float(bd))
    return {k: sorted(v) for k, v in out.items()}, worst


# ───────────────────────────── Part 1 ─────────────────────────────
def p1_fig2() -> dict:
    drs = page_drawings(P1, 13)
    xg = [it[1].x for it in drs[3]["items"]]
    yg = [it[1].y for it in drs[5]["items"]]
    assert len(xg) == 11 and len(yg) == 8
    xax = Axis.fit(xg, [5 + 2.5 * k for k in range(11)])
    yp = Axis.fit(yg, [20, 25, 30, 35, 40, 45, 50, 55])
    yv = Axis.fit(yg, [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
    assert colour_eq(drs[8]["color"], BLUE) and colour_eq(drs[9]["color"], GREEN)
    th, p = curve(drs[8], xax, yp)
    tv, v = curve(drs[9], xax, yv)
    return {
        "key": "hoth2021a_fig2_prf100_detonation_meter_input",
        "kind": "pressureTraceSingleCycleWithKnockmeterSignal",
        "source": "Hoth2021a", "figure": "Fig. 2", "page": 13,
        "caption": "Crank-angle resolved 501-C detonation-meter INPUT signal (D-1 pickup after the meter's input filter) and "
                   "cylinder pressure, PRF100, one knocking cycle (ringing visible), full cycle stored in the PDF",
        "conditions": conditions(
            engine="ANL CFR F1 (Hoth 2021)", method="RON", fuel=prf(100), rpm=600, sparkAdvanceDegBTDC=13,
            compressionRatio=None,
            compressionRatioSource="not stated (standard RON test on PRF100 -> guide table RON 100.0 = counter 919 at "
                                   f"29.92 inHg, CR {cr_rigid_raise(919):.3f} rigid raise, before barometric compensation)",
            notes="Standard RON conditions (peak-knock lambda); single representative cycle; cylinder pressure from the AVL "
                  "GU13Z-24 indicating spark plug; meter input signal logged simultaneously (volts)."),
        "extraction": {"method": "vector (page 13 paths #8 pressure (blue) and #9 meter input (green), 7199 line segments "
                                 "each = 0.1 CAD, full cycle; x axis from 11 grid lines 5..30, y from 8 tick segments "
                                 "(left 20..55 bar, right -0.6..0.8 V))",
                       "axisResidual": {"deg": xax.resid, "bar": yp.resid, "V": yv.resid},
                       "uncertainty": "digitisation < 0.01 deg / 0.002 bar / 0.0001 V"},
        "data": {"thetaDeg": th, "pPa": p * BAR, "meterInputThetaDeg": tv, "meterInputV": v},
        "derived": {"peakPressurePa": float(p.max() * BAR), "peakThetaDeg": float(th[np.argmax(p)]),
                    "meterInputPeakV": float(v.max()), "meterInputPeakThetaDeg": float(tv[np.argmax(v)])},
    }


def p1_fig6() -> dict:
    drs = page_drawings(P1, 19)
    words = page_words(P1, 19)
    panels = {
        "knockUnits": {"box": (100, 230, 292, 400), "ylab": (80, 228, 97, 404), "xlab": (90, 400, 300, 412)},
        "mapoBar": {"box": (335, 230, 527, 401), "ylab": (315, 228, 333, 405), "xlab": (325, 401, 535, 413)},
    }
    names = legend_names(drs, words, (330, 355, 400, 395))
    assert len(names) == 3, names
    out = {}
    for q, pnl in panels.items():
        xax = axis_from_labels(words, pnl["xlab"], "x")
        yax = axis_from_labels(words, pnl["ylab"], "y")
        trend = [d for d in series_by_style(drs) if in_box(d, pnl["box"])]
        curves = {names[style_key(d)]: curve(d, xax, yax) for d in trend}
        mk = short_segment_markers(drs, (0, 0, 0), max_len=8.0, exclude=[(340, 356, 400, 395)])
        mk = np.array([m for m in mk if pnl["box"][0] <= m[0] <= pnl["box"][2] and pnl["box"][1] <= m[1] <= pnl["box"][3]])
        mxy = np.column_stack([xax(mk[:, 0]), yax(mk[:, 1])])
        pts, worst = assign_to_curves(mxy, curves)
        out[q] = {k: {"lambda": [a for a, _ in v], q: [b for _, b in v], "maxDistanceToTrendline": worst[k],
                      "trendline": {"lambda": curves[k][0], q: curves[k][1]}} for k, v in pts.items()}
        out[q]["axisResidual"] = {"x": xax.resid, "y": yax.resid}
    for k in list(out["mapoBar"]):
        if k.startswith("PRF"):
            e = out["mapoBar"][k]
            e["mapoPa"] = [v * BAR for v in e.pop("mapoBar")]
            e["trendline"]["mapoPa"] = [v * BAR for v in e["trendline"].pop("mapoBar")]
    return {
        "key": "hoth2021a_fig6_ku_mapo_vs_lambda",
        "kind": "knockIntensityVsLambda",
        "source": "Hoth2021a", "figure": "Fig. 6", "page": 19,
        "caption": "Knockmeter (A) and MAPO (B) vs lambda for PRF93/95/97 at the ASTM RON-95 compression ratio, 13 deg bTDC",
        "conditions": conditions(
            engine="ANL CFR F1 (Hoth 2021)", method="RON (lambda sweep)", rpm=600, sparkAdvanceDegBTDC=13,
            digitalCounter=RON95_COUNTER,
            compressionRatio=None,
            compressionRatioSource=f"ASTM D2699 guide-table cylinder height for RON 95 (counter {RON95_COUNTER} at 29.92 inHg -> "
                                   f"CR {cr_rigid_raise(RON95_COUNTER):.3f} rigid raise), compensated for barometer (value not given); "
                                   "the later ANL paper OSTI 2561394 quotes CR 7.26 for 'standard RON 95 at 1.0 bar'",
            notes="Intake air temperature per D2699 for the barometer (no temperature tuning needed); lambda from LSU 4.9. "
                  "Knockmeter calibrated 12-15 KU/ON. MAPO: 4-18 kHz band-pass, 300-cycle mean."),
        "extraction": {"method": "vector: 'x' markers = clusters of 2-segment paths (centre), trendlines = dashed paths "
                                 "(series identified from the legend dash patterns), axes from tick-label centres; each "
                                 "marker assigned to the nearest trendline",
                       "uncertainty": {"lambda": 0.001, "knockUnits": 0.2, "mapoBar": 0.003}},
        "data": out,
        "reported": {"table7": {"PRF97": {"peakKnockLambda": 0.88, "peakMapoLambda": 0.88, "peakKU": 23.0, "peakMapoBar": 0.65},
                                "PRF95": {"peakKnockLambda": 0.88, "peakMapoLambda": 0.88, "peakKU": 52.8, "peakMapoBar": 0.90},
                                "PRF93": {"peakKnockLambda": 0.88, "peakMapoLambda": 0.88, "peakKU": 76.1, "peakMapoBar": 1.12},
                                "TSF96.9": {"RON": 97.2, "peakKnockLambda": 0.93, "peakKU": 20.2, "peakMapoBar": 0.22},
                                "TSF93.4": {"RON": 93.7, "peakKnockLambda": 0.93, "peakKU": 67.3, "peakMapoBar": 0.43}},
                     "table6_PRF98_standardRON": {
                         "sparkPlugTransducer_18kHz": {"mapoBar": 0.69, "cycleSD": 0.21},
                         "sparkPlugTransducer_40kHz_100kHzSampling": {"mapoBar": 0.79, "cycleSD": 0.23},
                         "kistler6044A_in_knockmeter_port_18kHz": {"mapoBar": 0.63, "cycleSD": 0.17},
                         "kistler6044A_40kHz": {"mapoBar": 0.67, "cycleSD": 0.18}},
                     "knockThreshold": "MAPO 0.1 bar separates knocking from non-knocking cycles; at standard RON (RON-95 "
                                       "CR, PKL) every cycle of PRF95/97 knocks; at stoichiometry PRF97 still knocks in 93 % "
                                       "of cycles",
                     "lambdaSensitivity": "near stoichiometry d(lambda) = 0.02 changes MAPO by 0.1 bar and KU by 10 (PRF93, "
                                          "Fig. 7); lambda SD at steady state +-0.015"},
    }


def p1_fig7() -> dict:
    drs = page_drawings(P1, 20)
    words = page_words(P1, 20)
    # bottom panel labels 0.985..1.010 are exact; the top panel prints them rounded to 2 decimals -> use the same values
    xb = axis_from_labels(words, (90, 310, 310, 322), "x")
    top_lab = [w for w in words if 221 < 0.5 * (w[1] + w[3]) < 233 and w[4].startswith(("0.9", "1.0"))]
    xt = Axis.fit(sorted(0.5 * (w[0] + w[2]) for w in top_lab), [0.985, 0.990, 0.995, 1.000, 1.005, 1.010])
    yt = axis_from_labels(words, (85, 130, 104, 225), "y")
    yb = axis_from_labels(words, (88, 220, 102, 314), "y")
    mk = [d for d in drs if d["type"] == "fs" and all(it[0] == "c" for it in d["items"]) and colour_eq(d.get("fill"), (0, 0, 0))]
    cen = np.array([((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2) for d in mk])
    top = cen[cen[:, 1] < 215]
    bot = cen[cen[:, 1] > 235]
    top = top[np.argsort(top[:, 0])]
    bot = bot[np.argsort(bot[:, 0])]
    return {
        "key": "hoth2021a_fig7_prf93_near_stoich",
        "kind": "knockIntensityVsLambda",
        "source": "Hoth2021a", "figure": "Fig. 7", "page": 20,
        "caption": "Knock intensity response to small lambda changes near stoichiometry, PRF93 at an increased compression ratio",
        "conditions": conditions(engine="ANL CFR F1 (Hoth 2021)", method="RON-like (stoichiometric)", fuel=prf(93), rpm=600,
                                 sparkAdvanceDegBTDC=13, compressionRatio=None,
                                 compressionRatioSource="'increased compression ratio' (value not given)"),
        "extraction": {"method": "vector: filled-circle marker centres; x from tick-label centres (top-panel labels are "
                                 "printed rounded, the bottom-panel values are used)",
                       "uncertainty": {"lambda": 0.0002, "mapoBar": 0.002, "knockUnits": 0.1}},
        "data": {"mapo": {"lambda": xt(top[:, 0]), "mapoPa": yt(top[:, 1]) * BAR},
                 "knockmeter": {"lambda": xb(bot[:, 0]), "knockUnits": yb(bot[:, 1])}},
    }


# ───────────────────────────── Part 2 ─────────────────────────────
def p2_fig1() -> dict:
    drs = page_drawings(P2, 10)
    words = page_words(P2, 10)
    xax = axis_from_labels(words, (100, 305, 300, 320), "x")
    yax = axis_from_labels(words, (90, 200, 108, 310), "y")
    mk = [d for d in drs if d["type"] == "fs" and colour_eq(d.get("fill"), (0, 0, 0)) and len(d["items"]) == 4]
    cen = np.array(sorted(((d["rect"].x0 + d["rect"].x1) / 2, (d["rect"].y0 + d["rect"].y1) / 2) for d in mk))
    tr = [d for d in drs if len(d["items"]) > 100][0]
    tx, ty = curve(tr, xax, yax)
    return {
        "key": "hoth2021b_fig1_prf100_mapo_vs_spark",
        "kind": "knockIntensityVsSparkTiming",
        "source": "Hoth2021b", "figure": "Fig. 1", "page": 10,
        "caption": "MAPO vs spark timing for iso-octane (PRF100) at otherwise standard RON conditions, with the RON (-13) and "
                   "modern-SI KLSA (0.6 bar MAPO) markers",
        "conditions": conditions(engine="ANL CFR F1 (Hoth 2021)", method="RON (spark sweep)", fuel=prf(100), rpm=600,
                                 compressionRatio=None,
                                 compressionRatioSource="standard RON conditions for PRF100 per the caption; the campaign "
                                                        f"otherwise used the RON-95 guide-table CR (counter {RON95_COUNTER}); not stated",
                                 notes="peak-knock lambda; MAPO 4-18 kHz, 300-cycle mean. The digitised spark timings "
                                       "coincide with Pal et al. 2018 Fig. 1 within 0.05 deg (both +0.2..+0.35 deg off the "
                                       "nominal integers), so this is very probably the same ANL PRF100 sweep (CR 7.55, "
                                       "lambda 0.89) re-evaluated with the 4-18 kHz MAPO; Pal's IndiCom MAPO is ~1.6x larger "
                                       "at -13 deg (2.27 vs 1.43 bar)"),
        "extraction": {"method": "vector: 20 diamond marker centres + dashed third-order trendline; axes from tick labels",
                       "uncertainty": {"sparkDeg": 0.02, "mapoBar": 0.003}},
        "data": {"sparkTimingCADATDC": xax(cen[:, 0]), "sparkAdvanceDegBTDC": -xax(cen[:, 0]), "mapoPa": yax(cen[:, 1]) * BAR,
                 "trendline": {"sparkTimingCADATDC": tx, "mapoPa": ty * BAR}},
        "reported": {"knockThresholds": "MAPO 0.6 bar (= 1 bar/1000 rpm at 600 rpm) and 40 KU"},
    }


def p2_fig3() -> dict:
    drs = page_drawings(P2, 14)
    words = page_words(P2, 14)
    # panels: (y-label region, x-label region, plot box, legend box)
    P = {
        ("PKL", "mapoBar"): ((80, 115, 97, 245), (90, 243, 290, 253), (100, 115, 290, 241), (215, 120, 290, 175)),
        ("stoich", "mapoBar"): ((315, 115, 332, 245), (325, 243, 520, 252), (333, 115, 515, 241), (445, 120, 515, 180)),
        ("PKL", "knockUnits"): ((80, 245, 97, 373), (90, 372, 290, 381), (100, 245, 290, 370), (215, 250, 290, 303)),
        ("stoich", "knockUnits"): ((315, 245, 332, 373), (325, 371, 520, 380), (333, 245, 515, 369), (460, 248, 515, 302)),
        ("PKL", "gIMEPbar"): ((78, 375, 97, 503), (90, 500, 290, 510), (100, 375, 290, 498), (215, 380, 290, 420)),
        ("stoich", "gIMEPbar"): ((315, 374, 332, 503), (325, 500, 520, 509), (333, 374, 515, 498), (340, 445, 400, 495)),
    }
    out: dict = {}
    for (lam, q), (ylab, xlab, box, leg) in P.items():
        xax = axis_from_labels(words, xlab, "x")
        yax = axis_from_labels(words, ylab, "y")
        names = legend_names(drs, words, leg)
        trend = [d for d in series_by_style(drs) if in_box(d, box)]
        curves = {}
        for d in trend:
            k = style_key(d)
            if k in names:
                curves[names[k]] = curve(d, xax, yax)
        mk = short_segment_markers(drs, (0, 0, 0), max_len=7.0, exclude=[leg])
        mk = np.array([m for m in mk if box[0] <= m[0] <= box[2] and box[1] <= m[1] <= box[3]])
        pts, worst = assign_to_curves(np.column_stack([xax(mk[:, 0]), yax(mk[:, 1])]) if len(mk) else np.zeros((0, 2)), curves)
        out.setdefault(lam, {})[q] = {
            k: {"sparkTimingCADATDC": [a for a, _ in v], q: [b for _, b in v], "maxDistanceToTrendline": worst[k],
                "trendline": {"sparkTimingCADATDC": curves[k][0], q: curves[k][1]}} for k, v in pts.items()}
        out[lam][q]["_legend"] = sorted(names.values())
        out[lam][q]["_axisResidual"] = {"x": xax.resid, "y": yax.resid}
        out[lam][q]["_nMarkers"] = int(len(mk))
    # convert to SI
    for lam in out:
        for q, unit in (("mapoBar", "mapoPa"), ("gIMEPbar", "gIMEPPa")):
            for k, e in out[lam][q].items():
                if k.startswith("PRF"):
                    e[unit] = [v * BAR for v in e.pop(q)]
                    e["trendline"][unit] = [v * BAR for v in e["trendline"].pop(q)]
    return {
        "key": "hoth2021b_fig3_prf_spark_sweeps",
        "kind": "knockIntensityVsSparkTiming",
        "source": "Hoth2021b", "figure": "Fig. 3", "page": 14,
        "caption": "Spark-timing sweeps of PRF87-PRF100 at peak-knock lambda (PKL, RON-95 guide-table CR) and at stoichiometry "
                   "(higher CR, fixed): MAPO, knockmeter reading and gIMEP",
        "conditions": {
            "PKL": conditions(engine="ANL CFR F1 (Hoth 2021)", method="RON-like spark sweep", rpm=600, digitalCounter=RON95_COUNTER,
                              compressionRatioSource=f"ASTM RON-95 cylinder height (counter {RON95_COUNTER} at 29.92 inHg, CR "
                                                     f"{cr_rigid_raise(RON95_COUNTER):.3f} rigid raise), barometer-compensated",
                              notes="each PRF at its peak-knock lambda (~0.88); intake air per D2699"),
            "stoich": conditions(engine="ANL CFR F1 (Hoth 2021)", method="RON-like spark sweep", rpm=600, **{"lambda": 1.0},
                                 compressionRatioSource="increased CR chosen so that FACE-G reaches 0.6 bar MAPO at max-gIMEP "
                                                        "phasing (Fig. 2); value not given",
                                 notes="lambda = 1"),
        },
        "extraction": {"method": "vector: markers = clusters of short black segments (x/+ glyphs), trendlines = dashed paths "
                                 "(series from legend dash patterns), axes from tick-label centres, marker -> nearest "
                                 "trendline",
                       "uncertainty": {"sparkDeg": 0.03, "mapoBar": 0.005, "knockUnits": 0.3, "gIMEPbar": 0.003,
                                       "note": "markers of different PRFs that overlap may be assigned to the wrong series; "
                                               "maxDistanceToTrendline flags it"}},
        "data": out,
    }


def p2_fig5() -> dict:
    drs = page_drawings(P2, 17)
    words = page_words(P2, 17)
    # panel A (full combustion event): x labels -15..35, left y 0..45 bar, right y 0..360 J/deg
    xa = axis_from_labels(words, (90, 258, 340, 268), "x")
    yha = axis_from_labels(words, (334, 110, 350, 260), "y")
    # panel B (zoom -15..3 CAD): pressure and RoHR share the 0..18 axis; its paths are stored at 0.1 CAD for the FULL
    # cycle (7199 segments) -> used for the pressure trace
    xb = axis_from_labels(words, (395, 250, 546, 260), "x")
    yb = axis_from_labels(words, (386, 108, 398, 252), "y")
    na = legend_names(drs, words, (105, 125, 180, 200), pattern=r"^PRF\d+")
    nb = legend_names(drs, words, (430, 150, 530, 200), pattern=r"^PRF\d+")
    out: dict = {}
    for d in drs:
        n = len(d["items"])
        k = style_key(d)
        if n >= 7000 and k in nb:  # panel B pressure, full cycle
            x, y = curve(d, xb, yb)
            s = (x >= -360) & (x < 360)
            out.setdefault(nb[k], {}).update({"thetaDeg": x[s], "pPa": y[s] * BAR})
        elif 2000 <= n < 2100 and k in na:  # panel A RoHR
            x, y = curve(d, xa, yha)
            out.setdefault(na[k], {}).update({"rohrThetaDeg": x, "rohrJPerDeg": y})
        elif 2400 <= n < 2500 and k in nb:  # panel B RoHR (same data, finer axis)
            x, y = curve(d, xb, yb)
            s = (x >= -40) & (x <= 5)
            out.setdefault(nb[k], {}).update({"rohrZoomThetaDeg": x[s], "rohrZoomJPerDeg": y[s]})
    return {
        "key": "hoth2021b_fig5_prf89_prf97_stoich_cycles",
        "kind": "pressureTraceRepresentativeCycle",
        "source": "Hoth2021b", "figure": "Fig. 5", "page": 17,
        "caption": "Representative cycles (selected on CA50, knock point, MAPO, gIMEP, dp/dtheta_max, p at spark) of PRF89 and "
                   "PRF97 at stoichiometric conditions and identical spark timing (the green-circled points of Fig. 3): "
                   "cylinder pressure (full cycle, from the panel-B paths) and apparent rate of heat release (panel A; "
                   "panel-B zoom for the pre-spark LTHR)",
        "conditions": conditions(engine="ANL CFR F1 (Hoth 2021)", method="RON-like (stoichiometric)", rpm=600, **{"lambda": 1.0},
                                 compressionRatioSource="stoichiometric-campaign CR (not given; see Fig. 3 notes)",
                                 notes="spark timing identical for both fuels (value only as a marker in the figure)"),
        "extraction": {"method": "vector: panel-B pressure paths (7199 segments = 0.1 CAD over the full cycle; the published "
                                 "view is clipped to -15..3 CAD / 0..18 bar) mapped with panel-B axes; RoHR paths of panel "
                                 "A (right axis, J/deg); series by legend line style", "uncertainty": "< 0.01 deg, < 0.01 bar"},
        "data": out,
        "reported": {"text": "CA50 of PRF89 ~3.5 CAD earlier than PRF97; PRF89 knock point almost at CA50; PRF89 shows "
                             "pre-spark low-temperature heat release (LTHR) from ~-15 CAD, PRF97 none"},
    }


# ───────────────────────────── 2025 ─────────────────────────────
def p3_fig2() -> dict:
    drs = page_drawings(P3, 3)
    words = page_words(P3, 3)
    xax = axis_from_labels(words, (55, 270, 265, 280), "x")
    yp = axis_from_labels(words, (48, 140, 61, 272), "y")
    yf = axis_from_labels(words, (259, 140, 268, 272), "y")
    p = [d for d in drs if len(d["items"]) > 7000][0]
    f = [d for d in drs if 1000 <= len(d["items"]) < 2000][0]
    x, y = curve(p, xax, yp)
    s = (x >= -360) & (x < 360)
    xf, yf_ = curve(f, xax, yf)
    return {
        "key": "hoth2025_fig2_prf98_pkl_cycle",
        "kind": "pressureTraceRepresentativeCycle",
        "source": "Hoth2025", "figure": "Fig. 2", "page": 3,
        "caption": "Representative cylinder pressure trace for PRF98 at peak-knock lambda (standard RON test) with knock onset "
                   "and knock point, and the band-pass-filtered, rectified pressure (MAPO definition)",
        "conditions": conditions(engine="ANL CFR F1 (2025)", method="RON", fuel=prf(98), rpm=600, sparkAdvanceDegBTDC=13,
                                 intakePressurePa=1.0e5, intakeAirTemperatureK=52 + T0,
                                 compressionRatioSource="standard RON-98 setting (not stated; compressed-air intake at 1.0 bara)",
                                 notes="compressed building air at 1.0 bara replaces ambient air in this campaign"),
        "extraction": {"method": "vector: 7198-segment pressure path (full cycle, 0.1 CAD) and 1200-segment filtered path; "
                                 "axes from tick labels", "uncertainty": "< 0.01 deg, < 0.02 bar"},
        "data": {"thetaDeg": x[s], "pPa": y[s] * BAR, "filtered": {"thetaDeg": xf, "rectifiedBandPassPa": yf_ * BAR}},
    }


def p3_tables() -> dict:
    t7 = {  # Table 7: CR at 0.6 bar MAPO (stoich, CA50 = 12 aTDC, 600 rpm, 1.0 bara, 52 C), max PRR, spark timing
        "FACE-B": ("Paraffinic", 7.43, 5.9, -12.5), "FACE-D": ("Aromatic", 7.63, 9.9, -10.7),
        "FACE-F": ("Paraffinic", 7.26, 7.1, -11.7), "FACE-G": ("Aromatic", 7.88, 9.9, -10.6),
        "FACE-A+E15": ("Paraffinic+E15", 7.31, 6.6, -12.3), "FACE-C+E15": ("Paraffinic+E15", 7.40, 7.3, -12.0),
        "FACE-H+E15": ("Aromatic+E15", 7.54, 9.2, -10.1), "PRF93": ("All paraffinic", 7.16, 6.2, -13.5),
        "PRF95": ("All paraffinic", 7.34, 5.2, -13.2), "PRF97": ("All paraffinic", 7.50, 3.9, -13.4),
        "TSF93.4": ("Aromatic", 7.65, 9.9, -9.1), "TSF96.9": ("Aromatic", 8.17, 9.8, -9.5),
    }
    rows = []
    for fuel, (cls, cr, prr, st) in t7.items():
        rows.append({"fuel": fuel, "class": cls, "criticalCR": cr, "maxPressureRiseRatePaPerDeg": prr * BAR,
                     "sparkTimingCADATDC": st, "sparkAdvanceDegBTDC": -st,
                     "octaneNumber": float(fuel[3:]) if fuel.startswith("PRF") else None})
    return {
        "key": "hoth2025_critical_cr",
        "kind": "knockLimitedCompressionRatio",
        "source": "Hoth2025", "tables": "Tables 3, 4, 7 (pp. 3, 6)",
        "conditions": conditions(engine="ANL CFR F1 (2025)", method="critical-CR (MAPO 0.6 bar)", rpm=600, **{"lambda": 1.0},
                                 intakePressurePa=1.0e5, intakeAirTemperatureK=52 + T0, coolantTemperatureK=100 + T0,
                                 notes="CR raised until the 300-cycle mean MAPO = 0.6 bar at stoichiometry with the spark "
                                       "adjusted for CA50 = 12 deg aTDC and IMEP held at 7.6 bar (intake pressure trimmed "
                                       "+-0.005 bar)"),
        "data": {"table7": rows,
                 "table3_FACE-G_example": {
                     "columns": ["Std. RON", "Stoich. RON", "Knockmeter threshold", "MAPO threshold"],
                     "compressionRatio": [7.26, 7.26, 7.31, 7.88], "sparkTimingCADATDC": [-13, -13, -12.4, -10.8],
                     "CA50CADATDC": [10.3, 11.7, 12.1, 11.8], "lambda": [0.95, 1.00, 1.00, 1.01],
                     "IMEPPa": [7.7e5, 7.6e5, 7.6e5, 7.6e5], "knockUnits": [27, 22, 26, 86], "mapoPa": [0.33e5, 0.31e5, 0.33e5, 0.6e5]},
                 "standardRON95CR_at_1bar": 7.26,
                 "note": "The ASTM RON-95 guide-table cylinder height corresponds to CR 7.26 on this engine per the authors "
                         f"(vs {cr_rigid_raise(RON95_COUNTER):.2f} from the project's rigid-raise counter relation fitted to Choi "
                         "2018) - an unresolved CR-calibration discrepancy between ANL papers; use as a +-0.4 CR systematic."},
    }


def main() -> None:
    ds = [p1_fig2(), p1_fig6(), p1_fig7(), p2_fig1(), p2_fig3(), p2_fig5(), p3_fig2(), p3_tables()]
    out = {
        "description": "ANL CFR RON-rating studies with simultaneous knockmeter and cylinder-pressure knock metrics: 501-C "
                       "input signal vs pressure, KU and MAPO vs lambda and spark timing for PRF87-100, representative "
                       "knocking cycles, critical CR at a MAPO threshold. SI units, Pa absolute, deg with 0 = firing TDC.",
        "generator": "tools/reference/cfr_data_hoth.py",
        "sources": SRC,
        "conventions": {"thetaDeg": "deg, 0 = firing TDC (source deg aTDC)", "pressure": "Pa absolute",
                        "sparkAdvanceDegBTDC": "positive = before TDC"},
        "datasets": ds,
    }
    write_fixture("cfr_hoth_knock_metrics.json", out, sig=6)  # 6 digits: 10 Pa on 1e6 Pa, far below the 100 Pa digitisation step
    f6 = ds[1]["data"]
    for q in ("knockUnits", "mapoBar"):
        print(q, {k: (len(v["lambda"]), round(v["maxDistanceToTrendline"], 3)) for k, v in f6[q].items() if k.startswith("PRF")})
    f3 = ds[4]["data"]
    for lam in f3:
        for q in f3[lam]:
            print(lam, q, f3[lam][q]["_legend"], f3[lam][q]["_nMarkers"],
                  {k: (len(v["sparkTimingCADATDC"]), round(v["maxDistanceToTrendline"], 3)) for k, v in f3[lam][q].items() if k.startswith("PRF")})
    print("fig5", {k: (len(v.get("thetaDeg", [])), round(max(v.get("pPa", [0])) / 1e5, 2), sorted(v)) for k, v in ds[5]["data"].items()})
    print("fig1", len(ds[3]["data"]["mapoPa"]))
    print("p1fig7", {k: len(v[list(v)[0]]) for k, v in ds[2]["data"].items()})
    print("fig2 derived", ds[0]["derived"])


if __name__ == "__main__":
    main()
