"""Pal et al. 2018 (SAE 2018-01-0187, SAE Int. J. Engines 11(6); OSTI 1572720): ANL CFR engine, iso-octane
(PRF100) at RON conditions, CR 7.55, lambda 0.89, spark-timing sweep -14..+6 CAD ATDC.

Source: tools/reference/cfr_sources/osti_1572720_pal2018.pdf (accepted manuscript printed to PDF; all figures are
embedded RASTER images, so everything below is raster digitisation; the tables are transcribed from the text layer).

Datasets
  * Fig. 1  knockmeter reading (KU) and MAPO vs spark timing: diamond markers located by colour mask + morphological
            opening with a diamond structuring element (removes the connecting lines), centroid. Axes: light-grey grid
            lines detected per panel. The marker x positions are NOT integers: the logged spark timing (coil current
            clamp) differs from the nominal 1-CAD steps by +0.2..+0.35 CAD; both are given.
  * Fig. 7  knocking case (ST -13): 300-cycle average (white), +-1 SD band (dark grey) and all-cycle envelope (light
            grey) between -10 and 30 CAD; representative cycle 169 (black).
  * Fig. 8  pressure spectra: resonance peaks of the 300-cycle average spectrum (white curve).
  * Fig. 9  knock point of each of the 300 experimental cycles (grey dots), mean (dashed) and +-SD band.
  * Tables 3 and 5, text.
Crank angle: source 'CAD ATDC' = firing-TDC convention (identical). Pressure: MPa in the source (absolute; the traces
are the same data as Choi 2018 Fig. 9, which are absolute) -> Pa.

Writes test/fixtures/cfr_pal2018_knock.json.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from cfr_data_common import (
    BAR, Axis, T0, blobs, conditions, find_lines, load_image, prf, source_pdf, write_fixture,
)

PDF = source_pdf("osti_1572720_pal2018.pdf")
CITATION = ("Pal P., Kolodziej C., Choi S., Som S., Broatch A., Gomez-Soriano J., Wu Y., Lu T., See Y.C., 'Development "
            "of a Virtual CFR Engine Model for Knocking Combustion Analysis', SAE Int. J. Engines 11(6):1069-1082, 2018, "
            "SAE 2018-01-0187; accepted manuscript OSTI 1572720")
URL = "https://www.osti.gov/servlets/purl/1572720"

BASE_COND = dict(
    engine="ANL CFR F1/F2", fuel=prf(100), compressionRatio=7.55,
    compressionRatioSource="stated (Pal 2018 p. 2: 'fixed at 7.55, for which standard knock intensity (~50 KU) was achieved')",
    rpm=600, **{"lambda": 0.89},
    notes="RON test conditions (600 rpm, intake air ~52 C compensated for barometer, 13 deg bTDC standard); lambda fixed "
          "at 0.89 = peak-knock lambda of iso-octane (prior lambda sweep). Cylinder pressure: AVL spark-plug transducer. "
          "Knock units: D-1 pickup + 501-C meter; MAPO from AVL IndiCom standard algorithm.",
)


def grey_mask(img, lo, hi, sat=10):
    i = img.astype(int)
    return (i.max(axis=2) - i.min(axis=2) <= sat) & (i.min(axis=2) >= lo) & (i.max(axis=2) <= hi)


def merge_close(bl: list[dict], dist: float) -> list[tuple[float, float]]:
    """Merge blob fragments whose centroids are within `dist` px (area-weighted)."""
    groups: list[list[dict]] = []
    for b in sorted(bl, key=lambda b: b["cx"]):
        for g in groups:
            if any(np.hypot(b["cx"] - o["cx"], b["cy"] - o["cy"]) < dist for o in g):
                g.append(b)
                break
        else:
            groups.append([b])
    out = []
    for g in groups:
        a = np.array([o["area"] for o in g], float)
        out.append((float(np.sum(a * [o["cx"] for o in g]) / a.sum()), float(np.sum(a * [o["cy"] for o in g]) / a.sum())))
    return out


def fig1() -> dict:
    img = load_image(PDF, 42).astype(int)
    grid = grey_mask(img, 195, 235, 8)
    blue = (img[:, :, 2] - img[:, :, 0]) > 60
    r = 5
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    op = ndimage.binary_opening(blue, structure=(np.abs(yy) + np.abs(xx)) <= r)
    pts = merge_close(blobs(op, min_area=20), 24.0)
    hl = find_lines(grid[:, 200:930], 0, 0.3)
    upper_y = [y for y in hl if 40 < y < 440]  # 60..0 KU
    lower_y = [y for y in hl if 470 < y < 842]  # 3.0..0.0 bar
    assert len(upper_y) == 7 and len(lower_y) == 7, (upper_y, lower_y)
    vu = find_lines(grid[40:420, :], 1, 0.3)
    vl = find_lines(grid[480:830, :], 1, 0.3)
    # grid columns at -10, -5, 0, 10 (the +5 line is hidden under the blue ST=5 marker line; -15 merges with the axis)
    xu = Axis.fit([v for v in vu if 300 < v][:3] + [vu[-1]], [-10, -5, 0, 10])
    xl = Axis.fit([v for v in vl if 300 < v][:3] + [vl[-1]], [-10, -5, 0, 10])
    yu = Axis.fit(upper_y, [60, 50, 40, 30, 20, 10, 0])
    ylo = Axis.fit(lower_y, [3.0, 2.5, 2.0, 1.5, 1.0, 0.5, 0.0])
    up = sorted([p for p in pts if p[1] < 450], key=lambda p: p[0])
    lo = sorted([p for p in pts if p[1] > 450], key=lambda p: p[0])
    assert len(up) == 21 and len(lo) == 21, (len(up), len(lo))
    st_u = xu([p[0] for p in up])
    st_l = xl([p[0] for p in lo])
    assert np.max(np.abs(st_u - st_l)) < 0.15
    st = 0.5 * (st_u + st_l)
    ku = yu([p[1] for p in up])
    mapo = ylo([p[1] for p in lo])
    nominal = np.arange(-14, 7)
    return {
        "key": "pal2018_fig1_ku_mapo_vs_spark",
        "kind": "knockIntensityVsSparkTiming",
        "figure": "Fig. 1", "page": 3,
        "caption": "Effect of spark timing retard on CFR knockmeter reading and knock overpressure (MAPO)",
        "conditions": conditions(**BASE_COND),
        "extraction": {
            "method": "raster (embedded image xref 42, 988x959 px): blue diamond markers via colour mask + binary opening "
                      "(diamond SE r=5 px), fragments merged within 25 px, area-weighted centroid; axes from light-grey grid "
                      "lines (upper: 7 lines 60..0 KU; lower: 7 lines 3.0..0 bar; x: -10,-5,0,10)",
            "pixelScale": {"degPerPx": xu.per_unit(), "kuPerPx": yu.per_unit(), "barPerPx": ylo.per_unit()},
            "uncertainty": {"sparkTimingDeg": 0.05, "knockUnits": 0.3, "mapoBar": 0.01,
                            "note": "+-1 px marker centroid + calibration residual; values are the published 300-cycle "
                                    "(KU: time-averaged) means"},
        },
        "data": {
            "sparkTimingNominalCADATDC": nominal,
            "sparkTimingDigitisedCADATDC": st,
            "sparkAdvanceDegBTDC": -st,
            "knockUnits": ku,
            "mapoPa": mapo * BAR,
        },
        "reported": {"text": "standard RON spark timing (-13 CAD ATDC) = 50 KU; knock not detectable by the knockmeter for ST "
                             "later than -8 CAD ATDC although the pressure transducer still detected slight knock; no knock "
                             "at +5 CAD ATDC (p. 3)"},
    }


def choi_average_st_m13():
    """300-cycle average of the same data (Choi 2018 Fig. 9, ST -13, spark-plug transducer, vector-exact), used as the
    anchor for the band detection. Generated by cfr_data_choi2018.py (runs first alphabetically)."""
    import json

    from cfr_data_common import FIXTURES

    d = json.loads((FIXTURES / "cfr_choi2018_traces.json").read_text())
    e = next(x for x in d["datasets"] if x["key"] == "choi2018_fig9_prf100_st_m13_sparkPlug")
    return np.array(e["data"]["thetaDeg"]), np.array(e["data"]["pPa"]) / 1e6


def _run_around(rows: np.ndarray, anchor: float, gap: int = 8):
    """Contiguous run (gaps <= `gap` px bridged: lines crossing the band) of `rows` that contains/nearest to `anchor`."""
    if rows.size < 3:
        return None
    runs = np.split(rows, np.nonzero(np.diff(rows) > gap)[0] + 1)
    best = min(runs, key=lambda r: 0 if r[0] - 2 <= anchor <= r[-1] + 2 else min(abs(r[0] - anchor), abs(r[-1] - anchor)))
    if not (best[0] - 6 <= anchor <= best[-1] + 6):
        return None
    return best


def fig7() -> dict:
    img = load_image(PDF, 156).astype(int)
    white = img.min(axis=2) >= 245
    bgm = grey_mask(img, 226, 240, 6)
    rows = np.nonzero(bgm.mean(axis=1) > 0.3)[0]
    cols = np.nonzero(bgm.mean(axis=0) > 0.3)[0]
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    sub = white[r0:r1 + 1, c0:c1 + 1]
    hl = [v + r0 for v in find_lines(sub, 0, 0.5)]
    vl = [v + c0 for v in find_lines(sub, 1, 0.5)]
    assert len(vl) == 7 and len(hl) == 4, (vl, hl)
    xax = Axis.fit(vl, [-5, 0, 5, 10, 15, 20, 25])
    yax = Axis.fit(hl, [5.0, 4.0, 3.0, 2.0])
    inv_y = Axis.fit([5.0, 4.0, 3.0, 2.0], hl)  # MPa -> row
    grid_cols = set()
    for v in vl:
        grid_cols.update(range(int(v) - 3, int(v) + 4))
    all_band = grey_mask(img, 172, 205, 10)
    sd_band = grey_mask(img, 135, 168, 10)
    black = img.max(axis=2) < 60
    th_avg, p_avg = choi_average_st_m13()
    th, sdlo, sdhi, alllo, allhi, c169 = [], [], [], [], [], []
    # legend box (upper left, image px x 150-490, y < 340): the data there lie below it (p < 2.9 MPa)
    legend_x, legend_ybot = (150, 490), 340
    for x in range(c0 + 2, c1 - 1):
        if x in grid_cols:
            continue
        top = legend_ybot if legend_x[0] <= x <= legend_x[1] else r0 + 2
        t = float(xax(x))
        anchor = float(inv_y(np.interp(t, th_avg, p_avg)))
        sd = _run_around(np.nonzero(sd_band[top:r1, x])[0] + top, anchor)
        al = _run_around(np.nonzero(all_band[top:r1, x] | sd_band[top:r1, x] | white[top:r1, x])[0] + top, anchor)
        bk = np.nonzero(black[top:r1, x])[0] + top
        th.append(t)
        sdlo.append(float(yax(sd[-1])) if sd is not None else np.nan)
        sdhi.append(float(yax(sd[0])) if sd is not None else np.nan)
        alllo.append(float(yax(al[-1])) if al is not None else np.nan)
        allhi.append(float(yax(al[0])) if al is not None else np.nan)
        if bk.size:
            runs = np.split(bk, np.nonzero(np.diff(bk) > 1)[0] + 1)
            run = min(runs, key=lambda r: abs(0.5 * (r[0] + r[-1]) - anchor))  # cycle 169 stays near the mean
            c169.append(float(yax(0.5 * (run[0] + run[-1]))))
        else:
            c169.append(np.nan)
    th = np.array(th)
    g = np.round(np.arange(np.ceil(th.min() * 10) / 10, th.max(), 0.1), 3)

    def grid(v):
        v = np.array(v)
        ok = np.isfinite(v)
        out = np.interp(g, th[ok], v[ok])
        out[(g < th[ok].min()) | (g > th[ok].max())] = np.nan
        return out

    sl, sh, al_, ah, c = grid(sdlo), grid(sdhi), grid(alllo), grid(allhi), grid(c169)
    avg = np.interp(g, th_avg, p_avg)
    return {
        "key": "pal2018_fig7_knocking_st_m13",
        "kind": "pressureTraceStatistics",
        "figure": "Fig. 7", "page": 7,
        "caption": "Representative knocking cycle 169 (ST -13 CAD ATDC), with the 300-cycle average of the same data set "
                   "from the vector trace of Choi 2018 Fig. 9 for reference. The +-SD band and all-cycle envelope of Fig. 7 "
                   "could NOT be digitised reliably (their edges are covered by four coloured cycle lines) and are omitted.",
        "conditions": conditions(**dict(BASE_COND, sparkAdvanceDegBTDC=13, knockIntensityKU=50)),
        "extraction": {
            "method": "raster (embedded image xref 156, 1052x654 px), per pixel column: cycle 169 = centre of the black "
                      "pixel run nearest the 300-cycle mean (Choi 2018 Fig. 9 vector trace, identical data set); legend "
                      "region masked; resampled to 0.1 deg",
            "pixelScale": {"degPerPx": xax.per_unit(), "barPerPx": 10 * yax.per_unit()},
            "uncertainty": {"thetaDeg": 0.05, "cycle169Bar": 0.1,
                            "note": "line half-width 1-2 px (0.01-0.02 MPa); cycle 169 is covered by other cycles in places "
                                    "(then linearly interpolated). The published cycle is sampled at 0.1 CAD, so the 6 kHz "
                                    "ringing (0.56 CAD period at 600 rpm) is only coarsely resolved in the source itself."},
        },
        "data": {
            "thetaDeg": g,
            "averagePa_ChoiFig9": avg * 1e6,
            "cycle169Pa": c * 1e6,
        },
        "reported": {"cycle169": "most representative cycle (least deviation in (dP/dt)max, energy of resonance, MAPO + "
                                 "spectrum); Table 5 experiment: MAPO 2.3 bar, knock point 11.15 CAD ATDC, E_res 15.86 kPa^2 s"},
    }


def fig8() -> dict:
    img = load_image(PDF, 157).astype(int)
    bgm = grey_mask(img, 226, 240, 6)
    rows = np.nonzero(bgm.mean(axis=1) > 0.3)[0]
    cols = np.nonzero(bgm.mean(axis=0) > 0.3)[0]
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    white = img.min(axis=2) >= 245
    sub = white[r0:r1 + 1, c0:c1 + 1]
    hl = [v + r0 for v in find_lines(sub, 0, 0.5)]
    vl = [v + c0 for v in find_lines(sub, 1, 0.5)]
    assert len(vl) == 4 and len(hl) == 5, (vl, hl)
    # log axis: grid at 10^-0.5, 10^0, 10^0.5, 10^1 kHz (labelled 0.3, 1.0, 3.2, 10.0)
    xax = Axis.fit(vl, [10 ** -0.5, 1.0, 10 ** 0.5, 10.0], log=True)
    yax = Axis.fit(hl, [180, 160, 140, 120, 100])
    bad_rows = set()
    for h in hl:
        bad_rows.update(range(int(h) - 3, int(h) + 4))
    f, spl = [], []
    for x in range(c0 + 2, c1 - 1):
        if min(abs(x - v) for v in vl) <= 3:
            continue
        ys = [y for y in np.nonzero(white[r0 + 2:r1 - 1, x])[0] + r0 + 2 if y not in bad_rows]
        if not ys:
            continue
        # legend box (upper-left) excluded below by frequency window
        f.append(float(xax(x)))
        spl.append(float(yax(np.median(ys))))
    f, spl = np.array(f), np.array(spl)
    sel = f > 2.5
    fs, ss = f[sel], spl[sel]
    # local maxima of the average spectrum above 2.5 kHz: 5-px median filter (removes single-column spikes, e.g. the
    # legend's white sample line), prominence >= 4 dB, peaks closer than 8 % in frequency merged (keep the highest)
    from scipy.ndimage import median_filter
    from scipy.signal import find_peaks

    sm = median_filter(ss, size=5, mode="nearest")
    pk, prop = find_peaks(sm, prominence=4.0)
    cand = sorted(zip(pk, prop["prominences"]), key=lambda t: fs[t[0]])
    merged: list[list] = []
    for i, pr in cand:
        if merged and fs[i] / fs[merged[-1][0]] < 1.08:
            if sm[i] > sm[merged[-1][0]]:
                merged[-1] = [i, max(pr, merged[-1][1])]
        else:
            merged.append([i, pr])
    peaks = [{"frequencyHz": float(fs[i] * 1e3), "splDb": float(sm[i]), "prominenceDb": float(p)} for i, p in merged]
    return {
        "key": "pal2018_fig8_spectrum_peaks",
        "kind": "knockResonanceFrequencies",
        "figure": "Fig. 8", "page": 7,
        "caption": "In-cylinder pressure spectra (ST -13 CAD ATDC): resonance peaks of the 300-cycle average spectrum",
        "conditions": conditions(**dict(BASE_COND, sparkAdvanceDegBTDC=13, knockIntensityKU=50)),
        "extraction": {
            "method": "raster (xref 157, 1052x646 px): white average-spectrum pixels (grid rows/columns excluded), "
                      "median row per column; log-x calibration on the half-decade grid lines; scipy find_peaks "
                      "(prominence >= 6 dB) above 2.5 kHz",
            "uncertainty": {"frequencyRelative": 0.01, "splDb": 1.0,
                            "note": "1 px = 0.57% in frequency; peak position +-2 px"},
        },
        "data": {"peaks": peaks, "averageSpectrum": {"frequencyHz": f * 1e3, "splDb": spl}},
        "context": {"expectedModes": "cylinder acoustic modes f = alpha_mn c / (pi B), B = 82.55 mm (Draper); "
                                     "Rockstroh 2018 reports dominant 6, 10, 14 kHz for PRF90 at standard knock"},
    }


def fig9() -> dict:
    img = load_image(PDF, 160).astype(int)
    bgm = grey_mask(img, 226, 240, 6)
    rows = np.nonzero(bgm.mean(axis=1) > 0.2)[0]
    cols = np.nonzero(bgm.mean(axis=0) > 0.2)[0]
    r0, r1, c0, c1 = rows.min(), rows.max(), cols.min(), cols.max()
    white = img.min(axis=2) >= 245
    sub = white[r0:r1 + 1, c0:c1 + 1]
    hl = [v + r0 for v in find_lines(sub, 0, 0.4)]
    vl = [v + c0 for v in find_lines(sub, 1, 0.4)]
    # vertical grid 50..250 cycles; horizontal: 13, 12, (11, 10 hidden by the SD band), 9 ...
    xax = Axis.fit(vl, [50, 100, 150, 200, 250][: len(vl)])
    # y grid visible: 13 (72), 12 (130), 10?..: use the two outer pairs found; 1 CAD = 58 px
    hvals = {72.0: 13, 130.0: 12, 246.0: 10, 304.0: 9}
    ys = [h for h in hl if round(h) in hvals]
    yax = Axis.fit(ys, [hvals[round(h)] for h in ys])
    dots = grey_mask(img, 90, 145, 12)
    dots[:, : c0 + 3] = False
    bl = blobs(dots, min_area=15)
    single = [b for b in bl if b["area"] < 90]
    multi = [b for b in bl if b["area"] >= 90]
    cyc = [float(xax(b["cx"])) for b in single]
    kp = [float(yax(b["cy"])) for b in single]
    # dashed mean line: darkest long horizontal structure
    dash = img.max(axis=2) < 80
    dash[:, : c0 + 3] = False
    mean_rows = find_lines(dash[r0:r1, c0:c1], 0, 0.3)
    band = grey_mask(img, 176, 200, 8)
    band_rows = np.nonzero(band[r0:r1, c0 + 5:c1 - 5].mean(axis=1) > 0.3)[0] + r0
    return {
        "key": "pal2018_fig9_knock_points",
        "kind": "knockOnsetDistribution",
        "figure": "Fig. 9", "page": 8,
        "caption": "Knock points (|d2p/dtheta2| > 2 MPa/CAD^2 first exceeded) of the 300 experimental cycles, ST -13",
        "conditions": conditions(**dict(BASE_COND, sparkAdvanceDegBTDC=13, knockIntensityKU=50)),
        "extraction": {
            "method": "raster (xref 160, 702x433 px): grey dot blobs (area < 90 px = single dots; larger blobs = "
                      "overlapping dots, not resolved), centroid; dashed mean line and SD band rows",
            "nSingleDots": len(single), "nMergedBlobs": len(multi),
            "uncertainty": {"knockPointDeg": 0.03, "cycle": 0.5},
        },
        "data": {
            "cycle": cyc, "knockPointCADATDC": kp,
            "meanKnockPointCADATDC_dashedLine": [float(yax(r + r0)) for r in mean_rows],
            "sdBandCADATDC": [float(yax(band_rows.max())), float(yax(band_rows.min()))] if band_rows.size else None,
            "statsOfSingleDots": {"mean": float(np.mean(kp)), "sd": float(np.std(kp, ddof=1)), "n": len(kp)},
        },
        "reported": {"criterion": "knock point = crank angle where |d2p/dtheta2| first exceeds 2 MPa/CAD^2 (same as Foong "
                                  "et al.)", "cycle169KnockPointCADATDC": 11.15, "table5AverageKnockPoint": 11.15},
    }


def main() -> None:
    tables = {
        "key": "pal2018_tables",
        "kind": "reportedScalars",
        "table1": {"RON": {"rpm": 600, "intakeAirTemperatureC": 52, "mixtureTemperature": "not controlled",
                           "sparkBTDC": 13}, "MON": {"rpm": 900, "intakeAirTemperatureC": 38, "mixtureTemperatureC": 149,
                                                     "sparkBTDC": [19, 26]}},
        "table3_GTPower_TPA": {
            "note": "GT-Power TPA (three-pressure-analysis) model estimates from measured cylinder/intake/exhaust pressure, "
                    "not direct measurements",
            "sparkCADATDC": [-13, 5],
            "linerK": [430.1, 426.8], "headK": [525.2, 493.2], "pistonK": [485.2, 461.8], "intakePortK": [387.8, 382.9],
            "exhaustPortK": [436.4, 441.6], "intakeValveK": [464.9, 440.4], "exhaustValveK": [519.6, 501.9],
            "TIVC_K": [136.9 + T0, 132.9 + T0], "trappedMassKg": [0.6280e-3, 0.6329e-3], "residualGasFraction": [0.06044, 0.05781],
        },
        "table5_knockCharacteristics": {
            "experimentMostRepresentativeCycle": {"mapoPa": 2.3e5, "knockPointCADATDC": 11.15, "energyOfResonancePa2s": 15.86e6},
            "simulation": {"mapoPa": 2.57e5, "knockPointCADATDC": 11.30, "energyOfResonancePa2s": 17.2e6},
        },
        "text": {
            "nonKnocking_ST5": "CA10 and CA50 differences simulation vs experiment 5.0% and 3.0%; peak pressure magnitude and "
                               "location predicted accurately (Fig. 6)",
            "knockPointCriterion": "|d2p/dtheta2| > 2 MPa/CAD^2",
            "cycles": "300 consecutive cycles per operating point",
        },
    }
    out = {
        "description": "Pal et al. 2018 (ANL CFR, iso-octane, CR 7.55, lambda 0.89, RON conditions): knockmeter and MAPO vs "
                       "spark timing, knocking-cycle statistics, resonance frequencies, knock-point distribution. SI units.",
        "generator": "tools/reference/cfr_data_pal2018.py",
        "source": {"key": "Pal2018", "citation": CITATION, "url": URL, "osti": "1572720",
                   "localPdf": "tools/reference/cfr_sources/osti_1572720_pal2018.pdf"},
        "conventions": {"thetaDeg": "deg, 0 = firing TDC (source 'CAD ATDC', identical)",
                        "pressure": "Pa absolute (source MPa/bar absolute)",
                        "sparkAdvanceDegBTDC": "positive = before TDC (source spark timing in CAD ATDC = -advance)"},
        "datasets": [fig1(), fig7(), fig8(), fig9(), tables],
    }
    write_fixture("cfr_pal2018_knock.json", out)
    f1 = out["datasets"][0]["data"]
    for s, k, m in zip(f1["sparkTimingDigitisedCADATDC"], f1["knockUnits"], f1["mapoPa"]):
        print(f"ST {s:6.2f}  KU {k:5.1f}  MAPO {m / 1e5:4.2f} bar")
    print("peaks", [(round(p["frequencyHz"]), round(p["splDb"], 1)) for p in out["datasets"][2]["data"]["peaks"]])
    f9 = out["datasets"][3]
    print("fig9", f9["extraction"]["nSingleDots"], f9["extraction"]["nMergedBlobs"], f9["data"]["statsOfSingleDots"],
          f9["data"]["meanKnockPointCADATDC_dashedLine"], f9["data"]["sdBandCADATDC"])


if __name__ == "__main__":
    main()
