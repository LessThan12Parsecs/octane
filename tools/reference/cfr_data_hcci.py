"""HCCI (lean autoignition) data from the standard, carburetted ANL CFR engine - autoignition-chemistry validation
targets at known compression ratio, intake temperature and pressure:

  [LTR]   'Characterization of low temperature reactions in the standard Cooperative Fuel Research (CFR) Engine' (ANL,
          accepted manuscript OSTI 1583125): PRF90 at lambda ~3, 600 rpm, compression ratio required for CA50 = 3 deg aTDC
          vs intake mixture temperature at 0.9/1.0/1.15/1.3 bar (Fig. 6); PRF 'transfer functions' HCCI fuel number vs CR
          (Fig. 4); PRF90 rate of heat release with LTHR at 52 C (Fig. 10).
  [Kal]   Kalvakala K., Pal P., et al., 'Numerical analysis of fuel effects on advanced compression ignition using a
          cooperative fuel research engine CFD model', ASME J. Energy Resour. Technol. (accepted manuscript OSTI 1962089):
          PRF90 at lambda 3, 600 rpm, BRON (1.3 bar, 52 C, CR 12.84) and BMON (1.0 bar, 149 C, CR 14.05), 300-cycle mean
          pressure and HRR (Figs. 4-5, experiment = solid lines).

All figures are RASTER images (MATLAB / Excel exports) -> raster digitisation:
  * markers: colour mask, holes filled (hollow markers), binary opening (removes lines and error bars), centroid;
    square/circle by bounding-box fill ratio. Axes from the light-grey grid lines (values assigned from the labels).
  * curves: colour mask; for Kal Figs. 4-5 the dashed CFD curve of the same colour is removed by discarding connected
    components shorter than 40 px; column-wise centre of the remaining run nearest to the previous point; the LTHR inset
    of Kal Fig. 4 is masked. Axes from the frame (x -40..40 CAD) and the tick marks (Kal: 1 MPa = 5 minor ticks).
Pressure: MPa/bar absolute (IVC ~ intake pressure) -> Pa. Crank angle: 'CAD ATDC' firing (identical convention).

Writes test/fixtures/cfr_hcci_autoignition.json.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from cfr_data_common import BAR, T0, Axis, blobs, conditions, find_lines, load_image, prf, source_pdf, write_fixture

LTR = source_pdf("osti_1583125_cfr_hcci_ltr.pdf")
KAL = source_pdf("osti_1962089_kalvakala.pdf")


def cmask(img, rgb, tol):
    return np.sqrt(np.sum((img.astype(float) - np.array(rgb, float)) ** 2, axis=2)) <= tol


def hue_mask(img, name):
    """Anti-aliasing-tolerant masks for saturated primary colours and black."""
    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    if name == "black":
        return (img.max(axis=2) < 170) & (img.max(axis=2) - img.min(axis=2) < 30)
    if name == "red":
        return (R > 150) & (R - G > 90) & (R - B > 90)
    if name == "blue":
        return (B > 150) & (B - R > 80) & (B - G > 40)
    if name == "green":
        return (G > 150) & (G - R > 90) & (G - B > 90)
    if name == "darkgreen":
        return (G > 80) & (G - R > 50) & (G - B > 50)
    if name == "navy":
        return (B > 80) & (B - R > 50) & (B - G > 50) & (R < 90)
    raise KeyError(name)


def grid_lines(img):
    g = (img.max(axis=2) - img.min(axis=2) <= 8) & (img.min(axis=2) >= 200) & (img.max(axis=2) <= 240)
    return find_lines(g, 0, 0.3), find_lines(g, 1, 0.3)


def snap(detected, approx, tol=4):
    out = []
    for a in approx:
        d = min(detected, key=lambda v: abs(v - a))
        out.append(d if abs(d - a) <= tol else a)
    return out


def hollow_markers(img, rgb, tol, r=4, min_area=60, box=None):
    m = hue_mask(img, rgb) if isinstance(rgb, str) else cmask(img, rgb, tol)
    if box is not None:
        keep = np.zeros_like(m)
        keep[box[1]:box[3], box[0]:box[2]] = True
        m &= keep
    f = ndimage.binary_fill_holes(m)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    op = ndimage.binary_opening(f, structure=(xx ** 2 + yy ** 2) <= r * r)
    return [b for b in blobs(op, min_area=min_area) if b["w"] < 40 and b["h"] < 40]


def ltr_fig6() -> dict:
    img = load_image(LTR, 55).astype(int)
    hl, vl = grid_lines(img)
    xax = Axis.fit(snap(vl, [115, 241, 496, 623, 750, 876]), [10, 11, 13, 14, 15, 16])
    yax = Axis.fit(snap(hl, [94, 165, 238, 383, 455, 527, 599]), [180, 160, 140, 100, 80, 60, 40])
    series = {"1.3 bar": ("black", 0, 1.3), "1.15 bar": ("red", 0, 1.15), "1.0 bar": ("blue", 0, 1.0),
              "0.9 bar": ("green", 0, 0.9)}
    legend = (120, 480, 330, 630)
    out = {}
    for name, (rgb, tol, pin) in series.items():
        bl = hollow_markers(img, rgb, tol, box=(118, 60, 875, 632))
        bl = [b for b in bl if not (legend[0] <= b["cx"] <= legend[2] and legend[1] <= b["cy"] <= legend[3])]
        bl.sort(key=lambda b: -b["cy"])
        out[name] = {"intakePressurePa": pin * BAR,
                     "intakeMixtureTemperatureK": [float(yax(b["cy"])) + T0 for b in bl],
                     "compressionRatio": [float(xax(b["cx"])) for b in bl]}
    return {
        "key": "ltr_fig6_prf90_hcci_cr_for_ca50",
        "kind": "autoignitionCompressionRatio",
        "source": "LTR", "figure": "Fig. 6", "page": 9,
        "caption": "Combinations of intake mixture temperature and compression ratio giving CA50 = 3 deg aTDC for PRF90 in "
                   "HCCI mode at four intake pressures",
        "conditions": conditions(engine="ANL CFR F1 (standard carburettor, MON mixture heater, boosted air)", method="HCCI",
                                 fuel=prf(90), rpm=600, **{"lambda": 3.0}, coolantTemperatureK=100 + T0,
                                 notes="excess-air ratio ~3 (Lund-Chevron HCCI protocol; lambda not re-stated per point); "
                                       "spark off; CA50 held at 3 +- 1 deg aTDC by the CR; intake temperature +-5 C; "
                                       "daily PRF90 repeat: CR uncertainty ~2 %"),
        "extraction": {"method": "raster (xref 55, 954x722 px): hollow markers filled, opened, centroid; axes from grid "
                                 "lines (CR 10-16, T 40-180 C)",
                       "uncertainty": {"compressionRatio": 0.01, "temperatureK": 0.3}},
        "data": out,
    }


def ltr_fig4() -> dict:
    img = load_image(LTR, 47).astype(int)
    hl, vl = grid_lines(img)
    xax = Axis.fit(snap(vl, [199, 284, 368, 454, 538, 623]), [8, 9, 10, 11, 12, 13])
    yax = Axis.fit(snap(hl, [130, 274, 346, 419, 491]), [100, 90, 85, 80, 75])
    cols = {"HCCI#1 (52 C)": "darkgreen", "HCCI#2 (180 C)": "navy"}
    legend = (130, 70, 390, 205)
    out = {}
    for name, rgb in cols.items():
        bl = hollow_markers(img, rgb, 70, r=5, min_area=80, box=(116, 60, 878, 634))
        bl = [b for b in bl if not (legend[0] <= b["cx"] <= legend[2] and legend[1] <= b["cy"] <= legend[3])]
        for b in bl:
            on = float(yax(b["cy"]))
            # PRF 88/93 were only run at 1 bar (squares), PRF 97 only at 1.3 bar (circles); otherwise shape decides
            if min(abs(on - 88), abs(on - 93)) < 1.0:
                shape = "square"
            elif abs(on - 97) < 1.0:
                shape = "circle"
            else:
                shape = "square" if b["fill"] > 0.88 else "circle"
            pin = 1.0 if shape == "square" else 1.3
            key = f"{name}, {pin} bar"
            e = out.setdefault(key, {"intakePressurePa": pin * BAR, "compressionRatio": [], "hcciFuelNumber": []})
            e["compressionRatio"].append(float(xax(b["cx"])))
            e["hcciFuelNumber"].append(float(yax(b["cy"])))
    # a series has one point per PRF: if a colour has two same-PRF points in one shape class, the lower-CR one is the
    # 1.3 bar point and the higher-CR one the 1.0 bar point (boost lowers the CR needed for CA50 = 3; see the figure)
    for name in cols:
        k1, k13 = f"{name}, 1.0 bar", f"{name}, 1.3 bar"
        for ka, kb in ((k1, k13), (k13, k1)):
            if ka not in out:
                continue
            e = out[ka]
            ons = [round(v) for v in e["hcciFuelNumber"]]
            for on in set(ons):
                idx = [i for i, v in enumerate(ons) if v == on]
                if len(idx) == 2:
                    lo, hi = sorted(idx, key=lambda i: e["compressionRatio"][i])
                    move = lo if ka == k1 else hi
                    f = out.setdefault(kb, {"intakePressurePa": (1.3 if kb == k13 else 1.0) * BAR, "compressionRatio": [],
                                            "hcciFuelNumber": []})
                    f["compressionRatio"].append(e["compressionRatio"].pop(move))
                    f["hcciFuelNumber"].append(e["hcciFuelNumber"].pop(move))
                    ons.pop(move)
    for e in out.values():
        o = np.argsort(e["compressionRatio"])
        e["compressionRatio"] = [e["compressionRatio"][i] for i in o]
        e["hcciFuelNumber"] = [round(e["hcciFuelNumber"][i]) for i in o]  # PRF octane numbers (integers)
    return {
        "key": "ltr_fig4_prf_hcci_transfer_function",
        "kind": "autoignitionCompressionRatio",
        "source": "LTR", "figure": "Fig. 4", "page": 7,
        "caption": "HCCI 'transfer function': PRF (iso-octane vol %) vs compression ratio giving CA50 = 3 deg aTDC; HCCI#1 = "
                   "600 rpm, 52 C; HCCI#2 = 600 rpm, 180 C; at 1.0 and 1.3 bar intake",
        "conditions": conditions(engine="ANL CFR F1", method="HCCI", rpm=600, **{"lambda": 3.0},
                                 notes="intake mixture temperature 52 C (#1) or 180 C (#2); PRF 70-97; y values are the PRF "
                                       "octane numbers tested (rounded to integers after digitisation)"),
        "extraction": {"method": "raster (xref 47, 954x722 px): hollow markers (square = 1 bar, circle = 1.3 bar) by colour",
                       "uncertainty": {"compressionRatio": 0.015}},
        "data": out,
    }


def trace(mask, xax, yax, cols, rows, min_comp=0, start_row=None, max_jump=25):
    m = mask.copy()
    m[: rows[0], :] = False
    m[rows[1]:, :] = False
    m[:, : cols[0]] = False
    m[:, cols[1]:] = False
    if min_comp:
        lab, n = ndimage.label(m, structure=np.ones((3, 3)))
        sl = ndimage.find_objects(lab)
        for i, s in enumerate(sl, start=1):
            if s is None:
                continue
            ext = max(s[0].stop - s[0].start, s[1].stop - s[1].start)
            if ext < min_comp:
                m[lab == i] = False
    xs, ys = [], []
    prev = start_row
    for x in range(cols[0], cols[1]):
        yy = np.nonzero(m[:, x])[0]
        if yy.size == 0:
            continue
        runs = np.split(yy, np.nonzero(np.diff(yy) > 1)[0] + 1)
        cs = [0.5 * (r[0] + r[-1]) for r in runs]
        y = min(cs, key=lambda c: abs(c - prev)) if prev is not None else max(cs)
        if prev is not None and abs(y - prev) > max_jump:
            continue
        prev = y
        xs.append(float(xax(x)))
        ys.append(float(yax(y)))
    return np.array(xs), np.array(ys)


def ltr_fig10() -> dict:
    img = load_image(LTR, 69).astype(int)
    hl, vl = grid_lines(img)
    xax = Axis.fit(snap(vl, [115, 241, 496, 623, 750, 876]), [-30, -25, -15, -10, -5, 0])
    yax = Axis.fit(snap(hl, [140, 223, 305, 387, 470, 552, 634]), [25, 20, 15, 10, 5, 0, -5])
    legend = (125, 70, 435, 185)
    out = {}
    zero_row = float((0.0 - yax.b) / yax.a)
    for name, col, cr in (("1.3 bar", "black", 12.54), ("1.15 bar", "red", 13.55), ("1.0 bar", "blue", 15.00)):
        m = hue_mask(img, col)
        m[legend[1]:legend[3], legend[0]:legend[2]] = False
        # symbols (circles/triangles/squares) sit on the curves; they are small -> the trace passes through them.
        # Start at the zero line (-30 CAD, before any heat release); skip the tick marks next to the left frame.
        x, y = trace(m, xax, yax, (126, 874), (60, 632), start_row=zero_row, max_jump=12)
        g = np.round(np.arange(-30.0, -2.95, 0.1), 2)  # the steep main-ignition rise near TDC is not tracked
        out[name] = {"compressionRatio": cr, "intakePressurePa": float(name.split()[0]) * BAR,
                     "thetaDeg": g, "rohrJPerDeg": np.interp(g, x, y)}
    return {
        "key": "ltr_fig10_prf90_hcci_lthr",
        "kind": "heatReleaseRate",
        "source": "LTR", "figure": "Fig. 10", "page": 12,
        "caption": "Rate of heat release (300-cycle average of per-cycle analyses) showing LTHR for PRF90 in HCCI mode at "
                   "52 C intake and three intake pressures (CR 12.54 / 13.55 / 15.00 for CA50 = 3 deg aTDC)",
        "conditions": conditions(engine="ANL CFR F1", method="HCCI", fuel=prf(90), rpm=600, **{"lambda": 3.0},
                                 intakeMixtureTemperatureK=52 + T0,
                                 notes="first-law apparent heat release from the Kistler 6045 pressure; LTHR start = 0.2 J/CAD, "
                                       "main HR threshold 5 J/CAD (markers in the figure)"),
        "extraction": {"method": "raster (xref 69, 954x722 px), colour masks (black solid, red dashed, blue dash-dot), "
                                 "continuity-tracked column centres, linear interpolation over dash gaps to 0.1 deg",
                       "uncertainty": {"thetaDeg": 0.05, "rohrJPerDeg": 0.1}},
        "data": out,
    }


def kal_fig(xref: int, fig: str, page: int, frame: tuple, pmax_mpa: float, hrr_max: float, zero_row: float,
            top_row: float, cond: dict, inset=None, topmost=False) -> dict:
    img = load_image(KAL, xref).astype(int)
    xax = Axis.fit([frame[0], frame[1]], [-40.0, 40.0])
    yp = Axis.fit([zero_row, top_row], [0.0, pmax_mpa])
    yh = Axis.fit([zero_row, top_row], [0.0, hrr_max])
    blue = hue_mask(img, "blue")
    if inset is not None:
        blue[inset[1]:inset[3], inset[0]:inset[2]] = False
    # pressure: upper part of the plot (rows above the HRR baseline region)
    if topmost:  # BMON: the experimental (solid) PRF90 pressure is the highest blue curve at every crank angle
        xs, ys = [], []
        prev = None
        for x in range(frame[0] + 2, frame[1] - 2):
            yy = np.nonzero(blue[int(top_row):int(zero_row) - 12, x])[0]
            if yy.size:
                run = np.split(yy, np.nonzero(np.diff(yy) > 1)[0] + 1)[0]
                yc = int(top_row) + 0.5 * (run[0] + run[-1])
                if prev is not None and abs(yc - prev) > 70:  # hidden under another curve in this column
                    continue
                prev = yc
                xs.append(float(xax(x)))
                ys.append(float(yp(yc)))
        xp, ypv = np.array(xs), np.array(ys)
    else:
        xp, ypv = trace(blue, xax, yp, (frame[0] + 2, frame[1] - 2), (int(top_row), int(zero_row) - 12), min_comp=40,
                        max_jump=40)
    g = np.round(np.arange(-40.0, 40.01, 0.1), 2)
    p = np.interp(g, xp, ypv)
    return {"figure": fig, "page": page, "conditions": cond,
            "data": {"thetaDeg": g, "pPa": p * 1e6},
            "extraction": {"method": f"raster (xref {xref}, {img.shape[1]}x{img.shape[0]} px), blue solid line (dashed CFD "
                                     "line removed as short components (BRON) / topmost blue curve with continuity (BMON), "
                                     "legend and inset masked), frame x -40..40, y from the tick marks",
                           "uncertainty": {"thetaDeg": 0.15, "pressureMPa": 0.03},
                           "note": "HRR curve not extracted (solid and dashed PRF90 HRR overlap at the peak)"}}


def main() -> None:
    kal_bron = conditions(engine="ANL CFR F1/F2 (Kistler 6045-AU20 in the knockmeter port)", method="HCCI (BRON)",
                          fuel=prf(90), rpm=600, **{"lambda": 3.0}, compressionRatio=12.84,
                          compressionRatioSource="stated (Kalvakala et al., Model Validation section)",
                          intakePressurePa=1.3 * BAR, intakeMixtureTemperatureK=52 + T0,
                          notes="CA50 = 3 +- 1 deg aTDC; experimental trace = 300-cycle average; CFD needed CR 11.5")
    kal_bmon = conditions(engine="ANL CFR F1/F2 (Kistler 6045-AU20 in the knockmeter port)", method="HCCI (BMON)",
                          fuel=prf(90), rpm=600, **{"lambda": 3.0}, compressionRatio=14.05,
                          compressionRatioSource="stated (Kalvakala et al., Model Validation section)",
                          intakePressurePa=1.0 * BAR, intakeMixtureTemperatureK=149 + T0,
                          notes="CA50 = 3 +- 1 deg aTDC; experimental trace = 300-cycle average; CFD needed CR 11.92")
    k4 = kal_fig(159, "Fig. 4", 35, (80, 713), 8.0, 400.0, 663.0, 15.0, kal_bron, inset=(95, 90, 400, 345))
    k5 = kal_fig(162, "Fig. 5", 36, (76, 593), 5.0, 300.0, 555.0, 10.0, kal_bmon, inset=(420, 40, 600, 125), topmost=True)
    k4.update(key="kalvakala_fig4_prf90_bron", kind="pressureTraceCycleAveraged", source="Kal",
              caption="PRF90, BRON-L3-600rpm: experimental (solid) cylinder pressure, 300-cycle average")
    k5.update(key="kalvakala_fig5_prf90_bmon", kind="pressureTraceCycleAveraged", source="Kal",
              caption="PRF90, BMON-L3-600rpm: experimental (solid) cylinder pressure, 300-cycle average")
    kal_text = {"key": "kalvakala_text", "kind": "reportedScalars", "source": "Kal",
                "experimentalCR_for_CA50_3": {"PRF90": {"BRON": 12.84, "BMON": 14.05}, "TH90": {"BRON": 13.62, "BMON": 13.33}},
                "cfdCR": {"PRF90": {"BRON": 11.5, "BMON": 11.92}, "TH90": {"BRON": 12.2, "BMON": 11.42}},
                "fuels": {"PRF90": "RON 90, S 0", "TH90": "toluene/n-heptane, RON 90, S 10.8"},
                "table1_valveTimingAsPrinted": {"IVO": "10 ATDC", "IVC": "34 ABDC", "EVO": "40 BBDC", "EVC": "15 ATDC",
                                                "note": "degree signs lost in the PDF text ('100 ATDC' = 10 deg)"},
                "pressureSampling": "Kistler 6045-AU20 replacing the knockmeter, 0.1 CAD encoder"}
    ds = [ltr_fig6(), ltr_fig4(), ltr_fig10(), k4, k5, kal_text]
    out = {
        "description": "HCCI autoignition on the standard ANL CFR engine: PRF90 CR for CA50 = 3 deg aTDC vs intake T and p, "
                       "PRF transfer functions, LTHR heat-release rates, and 300-cycle mean pressure traces at BRON/BMON "
                       "conditions (lambda 3). SI units.",
        "generator": "tools/reference/cfr_data_hcci.py",
        "sources": {
            "LTR": {"citation": "'Characterization of low temperature reactions in the standard Cooperative Fuel Research (CFR) "
                                "Engine' (Argonne), accepted manuscript OSTI 1583125",
                    "url": "https://www.osti.gov/servlets/purl/1583125", "localPdf": "tools/reference/cfr_sources/osti_1583125_cfr_hcci_ltr.pdf"},
            "Kal": {"citation": "Kalvakala K., Pal P., Wu Y., Kukkadapu G., et al., 'Numerical analysis of fuel effects on advanced "
                                "compression ignition using a cooperative fuel research engine CFD model', ASME J. Energy Resour. "
                                "Technol., accepted manuscript OSTI 1962089",
                    "url": "https://www.osti.gov/servlets/purl/1962089", "localPdf": "tools/reference/cfr_sources/osti_1962089_kalvakala.pdf"},
        },
        "datasets": ds,
    }
    write_fixture("cfr_hcci_autoignition.json", out)
    for k, v in ds[0]["data"].items():
        print(k, [round(t - T0, 1) for t in v["intakeMixtureTemperatureK"]], [round(c, 2) for c in v["compressionRatio"]])
    for k, v in ds[1]["data"].items():
        print(k, [round(c, 2) for c in v["compressionRatio"]], v["hcciFuelNumber"])
    for k, v in ds[2]["data"].items():
        r = np.array(v["rohrJPerDeg"]); t = np.array(v["thetaDeg"])
        s = (t > -26) & (t < -18)
        print(k, "LTHR peak", round(r[s].max(), 2), "at", t[s][np.argmax(r[s])])
    for d in (k4, k5):
        p = np.array(d["data"]["pPa"]); t = np.array(d["data"]["thetaDeg"])
        print(d["key"], "p(-40)", round(p[0] / 1e5, 2), "pmax", round(p.max() / 1e5, 2), "at", t[np.argmax(p)])


if __name__ == "__main__":
    main()
