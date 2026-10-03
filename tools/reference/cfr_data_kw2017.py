"""Kolodziej & Wallner 2017 (Combustion Engines 171(4):164-169, doi:10.19206/CE-2017-427; OSTI 1394801): four ~RON 98
fuels (PRF98, RON98Alk, PRF71E30, RON98E30) rated at standard RON conditions on the instrumented ANL CFR engine.

Source: tools/reference/cfr_sources/osti_1394801_kolodziej_wallner2017.pdf. All figures are embedded RASTER images
(Excel charts, 750-757 px wide) -> raster digitisation:
  * bar charts (Figs. 2-6, 9, 11): bar top edge (first bar-colour row from the top in the central 60 % of the bar),
    axes calibrated on the light-grey major grid lines (value of each detected line assigned from the printed labels,
    lines hidden behind bars are skipped). Uncertainty 1 px.
  * Fig. 1 knockmeter vs lambda: markers by colour; circles/diamonds/triangles/squares after morphological opening
    (removes the fitted trend lines of the same colour), centroid.
  * NOT extracted: Fig. 7 (P-T trajectory IVC -> spark; a derived ideal-gas bulk temperature), Figs. 8 and 10 (PRF98
    pressure around knock): the PRF98 navy line is hidden under the RON98Alk dashed line, the CA50 cross markers and
    the knock-point annotation over most of its length, and attempted traces were not reliable (errors > 1 bar).
    The vector-exact 300-cycle PRF98 RON trace of Choi et al. 2018 Fig. 2 (same ANL engine, same method, same era)
    covers this need (cfr_choi2018_traces.json).

Writes test/fixtures/cfr_kw2017_ron98.json.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from cfr_data_common import BAR, T0, Axis, blobs, conditions, find_lines, load_image, prf, source_pdf, write_fixture

PDF = source_pdf("osti_1394801_kolodziej_wallner2017.pdf")
CITATION = ("Kolodziej C.P., Wallner T., 'Combustion characteristics of various fuels during research octane number testing "
            "on an instrumented CFR F1/F2 engine', Combustion Engines 171(4):164-169, 2017, doi:10.19206/CE-2017-427; OSTI 1394801")
URL = "https://www.osti.gov/servlets/purl/1394801"
FUELS = ["PRF98", "RON98Alk", "PRF71E30", "RON98E30"]
FUEL_INFO = {  # Table 5
    "PRF98": {"kind": "PRF", "octaneNumber": 98, "label": "PRF98", "RON": 98.0, "MON": 98.0, "HoV_kJkg": 308, "LHV_MJkg": 43.3, "AFst": 15.1},
    "RON98Alk": {"kind": "gasoline", "label": "RON98Alk (alkylate: ~75 %v iso-octane, 20 %v other iso-paraffins, 3 %v n-butane)",
                 "RON": 97.8, "MON": 96.6, "HoV_kJkg": 309, "LHV_MJkg": 44.5, "AFst": 15.1},
    "PRF71E30": {"kind": "blend", "label": "PRF71E30 (70 %v PRF71 + 30 %v ethanol)", "RON": 97.8, "MON": None,
                 "HoV_kJkg": 519, "LHV_MJkg": 39.3, "AFst": 13.1},
    "RON98E30": {"kind": "gasoline", "label": "RON98E30 (full-boiling-range gasoline + 30 %v ethanol)", "RON": 97.4, "MON": 86.6,
                 "HoV_kJkg": 536, "LHV_MJkg": 38.2, "AFst": 12.9},
}
BLUE = (91, 155, 213)
ORANGE = (237, 125, 49)
GREY_BAR = (217, 217, 217)
NAVY = (50, 62, 82)


def cmask(img, rgb, tol=30.0):
    return np.sqrt(np.sum((img.astype(float) - np.array(rgb, float)) ** 2, axis=2)) <= tol


def grid_axis(img, assign: dict[int, float], lo=200, hi=235) -> Axis:
    """Fit a y axis on detected horizontal light-grey grid lines; `assign` maps approximate row -> value."""
    i = img.astype(int)
    g = (i.max(axis=2) - i.min(axis=2) <= 8) & (i.min(axis=2) >= lo) & (i.max(axis=2) <= hi)
    rows = find_lines(g, 0, 0.3)
    pairs = []
    for approx, val in assign.items():
        near = [r for r in rows if abs(r - approx) <= 4]
        if near:
            pairs.append((near[0], val))
    assert len(pairs) >= 2, (rows, assign)
    return Axis.fit([p[0] for p in pairs], [p[1] for p in pairs])


def bar_tops(img, rgb, yax: Axis, n: int, min_area=400, min_cy=0.0) -> list[float]:
    """Top edge of each bar of colour `rgb`. A bar may be split into several connected components by overlaid
    arrows/text boxes: components are grouped by x centre; the top is the highest row of the group's own components
    (so anti-aliased legend text of the same grey above the bar is ignored)."""
    m = cmask(img, rgb, 30.0)
    lab, _ = ndimage.label(m, structure=np.ones((3, 3), int))
    bl = [b for b in blobs(m, min_area=min_area) if b["h"] >= 1 and b["cy"] > min_cy and b["w"] > 25]
    groups: list[list[dict]] = []
    for b in sorted(bl, key=lambda b: b["cx"]):
        if groups and abs(b["cx"] - groups[-1][0]["cx"]) < 15:
            groups[-1].append(b)
        else:
            groups.append([b])
    assert len(groups) == n, (rgb, [[(round(b["cx"]), b["area"]) for b in g] for g in groups])
    out = []
    for g in groups:
        labels = set()
        for b in g:
            sub = lab[b["y0"]:b["y0"] + b["h"], b["x0"]:b["x0"] + b["w"]]
            vals, cnt = np.unique(sub[sub > 0], return_counts=True)
            labels.add(int(vals[np.argmax(cnt)]))
        gm = np.isin(lab, list(labels))
        x0 = min(b["x0"] for b in g)
        w = max(b["x0"] + b["w"] for b in g) - x0
        tops = [np.nonzero(gm[:, c])[0].min() for c in range(x0 + int(0.2 * w), x0 + int(0.8 * w)) if gm[:, c].any()]
        out.append(float(yax(np.median(tops) - 0.5)))  # top edge of the bar = boundary between rows
    return out


def bars(xref: int, assign: dict, series: dict[str, tuple], figure: str, caption: str, min_area=400, min_cy=0.0) -> dict:
    img = load_image(PDF, xref)
    yax = grid_axis(img, assign)
    data = {}
    for name, (rgb, unit_scale, unit) in series.items():
        v = bar_tops(img, rgb, yax, 4, min_area, min_cy)
        data[name] = {f: val * unit_scale for f, val in zip(FUELS, v)}
        data[name]["unit"] = unit
    return {"figure": figure, "caption": caption, "values": data,
            "extraction": {"method": f"raster bar tops (xref {xref}, {img.shape[1]}x{img.shape[0]} px), y from grid lines",
                           "valuePerPx": yax.per_unit(), "uncertaintyPx": 1.0}}


def knockmeter_bars(xref=79) -> dict:
    """Fig. 11: wide light-grey knockmeter bars behind blue (PPRR) and orange (KP_PK) bars; left axis 0-70 KU,
    right axis 0-8 bar/CAD or bar."""
    img = load_image(PDF, xref)
    yax = grid_axis(img, {18: 70, 72: 60, 126: 50, 180: 40, 234: 30, 288: 20, 342: 10, 396: 0})
    grey = cmask(img, GREY_BAR, 6.0)
    blue = cmask(img, BLUE, 30.0)
    orange = cmask(img, ORANGE, 30.0)
    front = blue | orange
    # knockmeter bar = grey or front pixels stacked from the baseline; find its top per bar column range
    bl_b = sorted([b for b in blobs(blue, min_area=400)], key=lambda b: b["cx"])
    bl_o = sorted([b for b in blobs(orange, min_area=400)], key=lambda b: b["cx"])
    assert len(bl_b) == 4 and len(bl_o) == 4
    ku = []
    for bb, bo in zip(bl_b, bl_o):
        # grey is visible to the left of the blue bar and between blue and orange bars
        x_lo, x_hi = bb["x0"] - 15, bo["x0"] + bo["w"] + 15
        tops = []
        for c in range(x_lo, x_hi):
            col = grey[:, c] | front[:, c]
            ys = np.nonzero(col[: int(yax.b) if False else 400])[0]
            if ys.size and grey[:, c].any():
                # contiguous from the baseline upwards
                run_top = 396
                for y in range(395, 0, -1):
                    if col[y]:
                        run_top = y
                    else:
                        break
                if run_top < 390:
                    tops.append(run_top)
        ku.append(float(yax(np.median(sorted(tops)[: max(3, len(tops) // 4)]) - 0.5)))
    right = lambda v: v * 8.0 / 70.0  # right axis 0..8 aligned with left 0..70
    pprr = [right(v) for v in bar_tops(img, BLUE, yax, 4)]
    kppk = [right(v) for v in bar_tops(img, ORANGE, yax, 4)]
    return {"figure": "Fig. 11", "caption": "ASTM knockmeter reading vs peak pressure rise rate (PPRR) and knocking pressure "
                                            "peak (KP_PK), spark-plug transducer, RON rating conditions",
            "values": {"knockmeterKU": dict(zip(FUELS, ku)),
                       "PPRR_barPerCAD": dict(zip(FUELS, pprr)), "PPRR_PaPerDeg": {f: v * BAR for f, v in zip(FUELS, pprr)},
                       "KP_PK_bar": dict(zip(FUELS, kppk)), "KP_PK_Pa": {f: v * BAR for f, v in zip(FUELS, kppk)}},
            "extraction": {"method": "raster (xref 79): blue/orange bar tops; knockmeter = top of the light-grey column "
                                     "stack behind them; right axis = left x 8/70",
                           "uncertainty": {"KU": 0.3, "barPerCAD": 0.03, "bar": 0.03}}}


def fig1_markers() -> dict:
    img = load_image(PDF, 43)
    i = img.astype(int)
    yax = grid_axis(img, {28: 70, 128: 60, 226: 50, 326: 40, 424: 30, 524: 20, 622: 10})
    g = (i.max(axis=2) - i.min(axis=2) <= 8) & (i.min(axis=2) >= 200) & (i.max(axis=2) <= 235)
    cols = find_lines(g, 1, 0.3)
    xax = Axis.fit([c for c in cols if 150 < c < 650][:5], [0.80, 0.85, 0.90, 0.95, 1.00])
    specs = {
        "PRF98": ((0, 0, 0), 60.0, 3),  # small black circles (~8 px) on a 4 px black line -> disc r=3
        "RON98Alk": (BLUE, 35.0, 5),  # diamonds
        "PRF71E30": ((165, 165, 165), 20.0, 5),  # grey triangles
        "RON98E30": (ORANGE, 35.0, 5),  # squares
    }
    out = {}
    masks = [(0, 0, 260, 190),  # legend (x0, y0, x1, y1) px
             (350, 435, 710, 500), (135, 535, 350, 620),  # the two annotation text boxes
             (400, 378, 700, 408), (340, 585, 485, 615)]  # the two black annotation arrows
    for fuel, (rgb, tol, r) in specs.items():
        m = cmask(img, rgb, tol)
        for x0b, y0b, x1b, y1b in masks:
            m[y0b:y1b, x0b:x1b] = False
        m[:, :105] = False  # y axis labels
        m[700:, :] = False  # x axis labels
        yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
        op = ndimage.binary_opening(m, structure=(xx ** 2 + yy ** 2) <= r * r)
        bl = [b for b in blobs(op, min_area=int(0.6 * np.pi * r * r))]
        # markers are compact; arrows / lines crossing are elongated
        bl = [b for b in bl if 2 * r <= b["w"] <= 24 and 2 * r <= b["h"] <= 24 and max(b["w"], b["h"]) < 1.8 * min(b["w"], b["h"])]
        bl.sort(key=lambda b: b["cx"])
        out[fuel] = {"lambda": [float(xax(b["cx"])) for b in bl], "knockUnits": [float(yax(b["cy"])) for b in bl]}
    return {"figure": "Fig. 1", "caption": "ASTM knockmeter readings for lambda sweep from 0.8 to 1.0 for each fuel "
                                           "(RON 98 compression ratio and intake conditions, spark 13 deg bTDC)",
            "values": out,
            "extraction": {"method": "raster (xref 43, 750x810 px): markers by colour + binary opening with a disc "
                                     "(removes the same-colour polynomial trend lines), centroid; axes from grid lines "
                                     "(y: 10-KU lines, x: 0.80-1.00 vertical lines)",
                           "uncertainty": {"lambda": 0.002, "KU": 0.5,
                                           "note": "markers hidden under other fuels' markers/lines are NOT recovered "
                                                   "(PRF98: only the 3 visible markers; others partially). The fitted "
                                                   "trend lines are not digitised."}},
            "reported": {"peakKnockLambda": {"PRF98": 0.89, "RON98Alk": 0.89, "PRF71E30": 0.93, "RON98E30": 0.93}}}


def main() -> None:
    cond = conditions(
        engine="ANL CFR F1/F2", method="RON", fuel=None, compressionRatio=None,
        compressionRatioSource="not stated; set per ASTM D2699 for a RON 98 fuel compensated for the day's barometer "
                               "(guide table RON 98.0 -> counter 867 at 29.92 in Hg)",
        rpm=600, sparkAdvanceDegBTDC=13, coolantTemperatureK=100 + T0, knockIntensityKU=50,
        notes="Standard ASTM D2699 knockmeter calibration on PRF98 (50 KU, spread 12-15 KU/ON); each fuel at its peak-knock "
              "lambda (Fig. 1/2). Intake air temperature set per D2699 for the barometer (mixture temperature measured "
              "after the carburettor, Fig. 3). Spark-plug pressure transducer + D-1/501-C knockmeter simultaneously.")
    fig2 = bars(46, {18: 1.00, 66: 0.95, 112: 0.90, 160: 0.85, 208: 0.80, 254: 0.75, 396: 0.60},
                {"lambda": (BLUE, 1.0, "-"), "fuelRate": (ORANGE, 1.0 / 3600.0, "kg/s")},
                "Fig. 2", "Lambda and fuel rate (FR) at RON rating conditions")
    fig3 = bars(44, {18: 30, 82: 25, 144: 20, 208: 15, 396: 0},
                {"mixtureTemperatureC": (BLUE, 1.0, "degC"), "carburettorTemperatureDropK": (ORANGE, 1.0, "K")},
                "Fig. 3", "Mixture temperature after the carburettor and temperature drop across the carburettor")
    fig3["values"]["mixtureTemperatureK"] = {f: v + T0 for f, v in fig3["values"]["mixtureTemperatureC"].items() if f != "unit"}
    fig3["values"]["intakeAirTemperatureK_upstream"] = {
        f: fig3["values"]["mixtureTemperatureC"][f] + fig3["values"]["carburettorTemperatureDropK"][f] + T0 for f in FUELS}
    fig4 = bars(56, {18: 8.1, 144: 8.0, 396: 7.8}, {"gIMEP": (BLUE, BAR, "Pa")}, "Fig. 4",
                "Gross IMEP at RON test conditions (peak-knock lambda)")
    fig5 = bars(60, {18: 28.0, 66: 27.5, 114: 27.0, 160: 26.5, 254: 25.5, 302: 25.0, 396: 24.0},
                {"indicatedThermalEfficiency": (ORANGE, 0.01, "-")}, "Fig. 5", "Indicated thermal efficiency")
    fig6 = bars(58, {20: 390, 60: 380, 102: 370, 144: 360, 186: 350, 228: 340, 270: 330, 396: 300},
                {"exhaustTemperatureC": (BLUE, 1.0, "degC")}, "Fig. 6", "Measured exhaust gas temperature")
    fig9 = bars(71, {18: 14, 72: 12, 126: 10, 180: 8, 342: 2, 396: 0},
                {"CA10": ((91, 155, 213), 1.0, "deg aTDC"), "CA50": (ORANGE, 1.0, "deg aTDC"),
                 "knockPoint": ((165, 165, 165), 1.0, "deg aTDC")},
                "Fig. 9", "CA10, CA50 and knock point (Swarts definition)", min_area=20, min_cy=45)
    fig11 = knockmeter_bars()
    fig1 = fig1_markers()

    datasets = [
        {"key": "kw2017_ron98_operating_metrics", "kind": "reportedScalarsPerFuel", "fuels": FUEL_INFO,
         "conditions": cond, "figures": [fig2, fig3, fig4, fig5, fig6, fig9, fig11],
         "reportedText": {"peakKnockLambda": "~0.89 iso-paraffinic, ~0.93 E30 fuels (Fig. 1)",
                          "RONrating": "PRF98 98.0 (reference), RON98Alk 97.8, PRF71E30 97.8, RON98E30 97.4",
                          "knockmeterLinearity": "501-C loses linearity < 20 KU and > 80 KU"}},
        {"key": "kw2017_fig1_ku_vs_lambda", "kind": "knockIntensityVsLambda", "conditions": cond, **fig1},
    ]
    out = {
        "description": "Kolodziej & Wallner 2017: ANL CFR engine at standard RON conditions for four RON-98 fuels: lambda, "
                       "fuel rate, mixture temperature, gIMEP, ITE, exhaust T, CA10/CA50/knock point, knockmeter vs PPRR "
                       "and knock pressure peak, knockmeter vs lambda sweeps.",
        "generator": "tools/reference/cfr_data_kw2017.py",
        "source": {"key": "KW2017", "citation": CITATION, "url": URL, "osti": "1394801",
                   "localPdf": "tools/reference/cfr_sources/osti_1394801_kolodziej_wallner2017.pdf"},
        "conventions": {"thetaDeg": "deg, 0 = firing TDC (source dATDC)", "pressure": "Pa absolute (source bar)",
                        "fuelRate": "kg/s (source kg/h)"},
        "datasets": datasets,
    }
    write_fixture("cfr_kw2017_ron98.json", out)
    for f in (fig2, fig3, fig4, fig5, fig6, fig9, fig11):
        print(f["figure"], {k: ({kk: round(vv, 4) for kk, vv in v.items() if kk != "unit"}) for k, v in f["values"].items()})
    print("fig1", {k: len(v["lambda"]) for k, v in fig1["values"].items()})


if __name__ == "__main__":
    main()
