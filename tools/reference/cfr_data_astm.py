"""ASTM D2699-15a (RON) and D2700-14 (MON) guide tables, barometric compensation, motored compression-
pressure specifications and knock-intensity definitions for the CFR F-1/F-2 engine.

Sources (public copies incorporated by reference, archive.org; PDFs in tools/reference/cfr_sources/):
  [D2699] ASTM D2699-15a, archive.org item gov.law.astm.D2699.15A (astm.D2699.15A.pdf, 47 pp, born-digital XPP)
  [D2700] ASTM D2700-14,  archive.org item gov.law.astm.D2700.14  (astm.D2700.14.pdf,  58 pp, born-digital XPP)

Extraction
  * Tables: transcribed programmatically from the PDF text layer (PyMuPDF word boxes grouped into rows) ->
    exact (the numbers are the published digits). Row/column completeness is asserted below.
  * D2699 Fig. 2 / D2700 Fig. 2 (basic compression pressure vs barometric pressure) and Fig. A2.6 (spread
    characteristic): embedded raster images (no vector data), digitised by colour mask; axes calibrated on the
    detected major grid lines / tick marks. Uncertainty stated per dataset.
  * Counter -> CR: cfr.ts cfrCompressionRatioAtCounter (rigid raise, replicated in cfr_data_common.py) and, for
    comparison, the Choi et al. 2018 polynomial (cfrCompressionRatioFromCounter).

Writes test/fixtures/cfr_astm_guide_tables.json, cfr_astm_compression_pressure.json, cfr_astm_knock_intensity.json.
"""
from __future__ import annotations

import collections
import re

import numpy as np

from cfr_data_common import (
    BAR, IN_HG, MM_H2O, PSI, T0, Axis, blobs, colour_mask, counter_from_dial, cr_choi_polynomial, cr_rigid_raise,
    find_lines, load_image, source_pdf, trace_curve, write_fixture,
)

D2699 = source_pdf("astm_D2699-15a.pdf")
D2700 = source_pdf("astm_D2700-14.pdf")
URL_D2699 = "https://archive.org/details/gov.law.astm.D2699.15A"
URL_D2700 = "https://archive.org/details/gov.law.astm.D2700.14"
NUM = re.compile(r"^-?\d+(\.\d+)?$")


def page_rows(pdf, page_no: int, ytol: float = 1.5, with_x: bool = False) -> list[list]:
    """Words of a page grouped into text rows (sorted by x); with_x -> (x_centre, text) tuples."""
    import pymupdf

    words = pymupdf.open(pdf)[page_no - 1].get_text("words")
    rows: dict[float, list] = collections.OrderedDict()
    for w in sorted(words, key=lambda w: (round(w[3], 1), w[0])):
        key = next((k for k in rows if abs(k - w[3]) <= ytol), None)
        if key is None:
            key = w[3]
            rows[key] = []
        rows[key].append((0.5 * (w[0] + w[2]), w[4].replace("\u2013", "-").replace("\u2212", "-")))
    out = []
    for _, r in sorted(rows.items()):
        r = sorted(r)
        out.append(r if with_x else [t for _, t in r])
    return out


def grid_table(pdf, pages: list[int], label_is_float: bool, value_float: bool) -> dict[float, list]:
    """Rows 'label v0 .. v9' (v may be '...'); values are assigned to the 10 columns by their x position
    relative to the column-header row ('0.0 0.1 .. 0.9' or '0 1 .. 9'): returns {label: [10 values or None]}."""
    out: dict[float, list] = {}
    for pg in pages:
        rows = page_rows(pdf, pg, with_x=True)
        hdr = None
        for r in rows:
            toks = [t for _, t in r]
            if toks[-10:] in (["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9"],
                              [str(i) for i in range(10)]):
                hdr = [x for x, _ in r[-10:]]
                break
        assert hdr is not None, f"no column header on page {pg}"
        for r in rows:
            toks = [t for _, t in r]
            if len(toks) < 2 or not NUM.match(toks[0]):
                continue
            vals = r[1:]
            if not all(NUM.match(t) or t == "..." for _, t in vals) or len(vals) > 10 or len(vals) < 3:
                continue
            label = float(toks[0])
            conv = (lambda v: float(v)) if value_float else (lambda v: int(v))
            row: list = [None] * 10
            for x, t in vals:
                k = int(np.argmin([abs(x - h) for h in hdr]))
                assert row[k] is None and abs(x - hdr[k]) < 12, (pg, toks)
                row[k] = conv(t) if t != "..." else None
            out[label] = row
    return out


def on_to_counter_table(pdf, pages, what: str) -> dict:
    """Guide table 'octane number -> counter (or dial)', ON rows 40..120 x 10 tenths."""
    t = grid_table(pdf, pages, True, what == "dial")
    on, val = [], []
    for base in sorted(t):
        for k, v in enumerate(t[base]):
            if v is not None:
                on.append(round(base + 0.1 * k, 1))
                val.append(v)
    assert np.all(np.diff(on) > 0)
    return {"octaneNumber": on, what: val}


def counter_to_on_table(pdf, pages) -> dict:
    """Inverse guide table 'counter -> ON', rows of 10 counts."""
    t = grid_table(pdf, pages, False, True)
    c, on = [], []
    for base in sorted(t):
        for k, v in enumerate(t[base]):
            if v is not None:
                c.append(int(base) + k)
                on.append(v)
    assert np.all(np.diff(c) == 1), "counter rows must be contiguous"
    return {"counter": c, "octaneNumber": on}


def baro_table(pdf, page: int, with_iat: bool) -> list[dict]:
    """Barometric compensation rows. Returns list of {inHg, counterCorrection, dialCorrection[, iatC, iatF]}."""
    rows = page_rows(pdf, page)
    out = []
    i = 0
    base = None
    pending: dict[str, list] = {}

    def flush():
        if base is None or "dc" not in pending:
            return
        for k in range(10):
            e = {"barometerInHg": round(base + 0.1 * k, 1), "counterCorrection": pending["dc"][k],
                 "dialCorrectionIn": pending["dial"][k]}
            if with_iat:
                e["intakeAirTemperatureC"] = pending["iatC"][k]
                e["intakeAirTemperatureF"] = pending["iatF"][k]
            out.append(e)

    for r in rows:
        s = " ".join(r)
        m = re.match(r"^(\d\d\.0) \(([\d.]+)\) (Digital counter correction|Dial indicator correction) (.*)$", s)
        if m:
            if m.group(3).startswith("Digital"):
                flush()
                pending = {}
            base = float(m.group(1))
            key = "dc" if m.group(3).startswith("Digital") else "dial"
            vals = m.group(4).split()
        elif s.startswith("Digital counter correction"):
            flush()
            pending = {}
            key, vals = "dc", s.split()[3:]
        elif s.startswith("Dial indicator correction"):
            key, vals = "dial", s.split()[3:]
        elif s.startswith("IAT, °C"):
            key, vals = "iatC", s.split()[2:]
        elif s.startswith("IAT, °F"):
            key, vals = "iatF", s.split()[2:]
        else:
            continue
        vals = [v.replace("–", "-") for v in vals]
        assert len(vals) == 10, (s, vals)
        pending[key] = [int(v) if key in ("dc", "iatF") else float(v) for v in vals]
    flush()
    return out


# ───────────────────────── raster figures ─────────────────────────
def fig2_compression_pressure() -> list[dict]:
    """D2699 Fig. 2 (one line) and D2700 Fig. 2 (three lines): basic compression pressure (psig) at the basic
    cylinder height (counter 930 / dial 0.352 in) vs barometric pressure (in Hg)."""
    out = []
    specs = [
        (D2699, 50, "D2699-15a Fig. 2", 10, "RON (D2699), 9/16 in venturi", [(21, 31)]),
        (D2700, 53, "D2700-14 Fig. 2", 11, "MON (D2700)", [(21, 27), (26, 29), (28, 31)]),
    ]
    venturi_lbl = {0: "3/4 in venturi, 26.00 in Hg base", 1: "19/32 in venturi, 28.00 in Hg base",
                   2: "9/16 in venturi, 29.92 in Hg base"}
    for pdf, xref, fig, page, label, segs in specs:
        img = load_image(pdf, xref)
        g = img.mean(axis=2)
        # major grid lines (darker than minor): x = 21..31 in Hg, y = ymin..ymax psig in 10 psig steps
        vx = find_lines(g < 100, 1, 0.5)
        if len(vx) != 11:  # D2700: text boxes interrupt some major lines -> accept lighter pixels
            vx = find_lines(g < 200, 1, 0.5)
        hy = find_lines(g < 200, 0, 0.5)
        assert len(vx) == 11, vx
        xax = Axis.fit(vx, list(range(21, 32)))
        if pdf is D2699:
            ylines = [y for y in find_lines(g < 100, 0, 0.5)][:8]  # 210 .. 140
            yvals = list(range(210, 139, -10))
        else:
            ylines = [y for y in hy if y not in ()]
            ylines = [9.5, 76.5, 144.0, 213.5, 279.0, 346.5, 415.0, 481.5]
            ylines = [min(hy, key=lambda h: abs(h - y)) for y in ylines]  # snap to detected
            yvals = list(range(190, 119, -10))
        yax = Axis.fit(ylines, yvals)
        # curve pixels: black; remove the major grid lines (also near-black in D2699) and keep the points that lie
        # within 3 psig of a coarse visual guess of each straight segment (rejects label text), then least squares.
        m = g < 60
        for v in vx:
            m[:, max(0, int(round(v)) - 2): int(round(v)) + 3] = False
        for h in ylines:
            m[max(0, int(round(h)) - 2): int(round(h)) + 3, :] = False
        x0, x1 = int(vx[0]) + 3, int(vx[-1]) - 2
        y0, y1 = int(ylines[0]) + 3, int(ylines[-1]) - 2
        px, py = [], []
        for x in range(x0, x1):
            col = np.nonzero(m[y0:y1, x])[0]
            if col.size == 0:
                continue
            runs = np.split(col, np.nonzero(np.diff(col) > 1)[0] + 1)
            for r in runs:
                px.append(x)
                py.append(y0 + 0.5 * (r[0] + r[-1]))
        px, py = np.asarray(px, float), np.asarray(py, float)
        guesses = [(21, 147.7, 31, 209.0)] if pdf is D2699 else \
            [(21, 135.5, 27, 172.3), (26, 159.0, 29, 176.8), (28, 165.0, 31, 182.5)]
        lines = []
        for (xa, ya, xb, yb) in guesses:
            X, Y = xax(px), yax(py)
            yg = ya + (yb - ya) * (X - xa) / (xb - xa)
            sel = (X >= xa + 0.05) & (X <= xb - 0.05) & (np.abs(Y - yg) < 3.0)
            lines.append((px[sel], py[sel]))
        assert len(lines) == len(segs), (fig, len(lines))
        for k, (xs, ys) in enumerate(lines):
            inhg = xax(xs)
            psig = yax(ys)
            A = np.vstack([inhg, np.ones_like(inhg)]).T
            (slope, icpt), *_ = np.linalg.lstsq(A, psig, rcond=None)
            resid = psig - A @ np.array([slope, icpt])
            grid = np.round(np.arange(np.ceil(inhg.min() * 10) / 10, np.floor(inhg.max() * 10) / 10 + 1e-9, 0.1), 1)
            fit = slope * grid + icpt
            ent = {
                "figure": fig, "page": page, "method": label if pdf is D2699 else f"MON (D2700), {venturi_lbl[k]}",
                "counter": 930, "dialIndicatorIn": 0.352,
                "linearFit": {"psigPerInHg": slope, "psigAt0InHg": icpt, "rmsResidualPsi": float(np.sqrt(np.mean(resid ** 2))),
                              "maxResidualPsi": float(np.max(np.abs(resid)))},
                "barometerInHg": grid.tolist(),
                "barometerPa": (grid * IN_HG).tolist(),
                "compressionPressurePsig": fit.tolist(),
                "compressionPressureGaugePa": (fit * PSI).tolist(),
                "compressionPressureAbsPa": (fit * PSI + grid * IN_HG).tolist(),
                "extraction": {
                    "method": "raster digitisation (embedded image xref %d, %dx%d px, black-line mask, column-wise centre, "
                              "straight-line least-squares fit; axes from the 11 vertical / 8 horizontal major grid lines)"
                              % (xref, img.shape[1], img.shape[0]),
                    "pixelsPerInHg": 1 / xax.per_unit(), "pixelsPerPsi": 1 / yax.per_unit(),
                    "uncertaintyPsi": float(max(0.5 * yax.per_unit() * 2, np.max(np.abs(resid)))),
                    "note": "The published curve is a straight line; the fit residual is the digitisation noise. "
                            "Uncertainty = max(1 px, max residual).",
                },
            }
            out.append(ent)
    return out


def spread_figure(pdf, xref, fig, page, xticks, xvals, yticks, yvals, axis_x, axis_y) -> dict:
    """Fig. A2.6 'Typical detonation meter spread characteristic' (K.I. divisions per O.N. vs O.N.).
    Non-uniform printed x axis -> piecewise-linear calibration through the tick marks."""
    img = load_image(pdf, xref)
    g = img.mean(axis=2)
    m = g < 110
    m[:, : int(axis_x) + 10] = False  # remove y axis and its inside tick marks
    m[int(axis_y) - 10:, :] = False  # remove x axis and its inside tick marks
    xs, ys, hw = trace_curve(m, x_range=(int(axis_x) + 10, img.shape[1]))
    on = np.interp(xs, xticks, xvals)
    spread = Axis.fit(yticks, yvals)(ys)  # uniform y ticks: linear fit (extrapolates below the lowest tick)
    grid = np.round(np.arange(np.ceil(on.min() * 10) / 10, np.floor(on.max() * 10) / 10 + 1e-9, 0.1), 1)
    sp = np.interp(grid, on, spread)
    return {
        "figure": fig, "page": page,
        "octaneNumber": grid.tolist(), "spreadKIDivisionsPerON": sp.tolist(),
        "extraction": {
            "method": "raster digitisation (embedded image xref %d, %dx%d px, dark-pixel mask, column-wise centre of "
                      "the largest run; x piecewise-linear through tick marks %s px = %s O.N. (the printed axis is not "
                      "uniform), y linear through ticks %s px = %s div/ON)" % (xref, img.shape[1], img.shape[0], xticks,
                                                                                xvals, yticks, yvals),
            "uncertaintyKIDivPerON": 0.4,
            "note": "Schematic 'typical' curve (ASTM wording); line half-width ~1.5 px = 0.15 div/ON; the tick "
                    "spacing irregularity (up to 5 px = 0.3 ON) is absorbed by the piecewise calibration.",
        },
    }


def main() -> None:
    # ── guide tables ──
    ron_c = on_to_counter_table(D2699, [31, 32], "counter")
    ron_c2on = counter_to_on_table(D2699, [33, 34])
    ron_dial = on_to_counter_table(D2699, [35, 36], "dial")
    assert ron_c["octaneNumber"][0] == 40.0 and ron_c["octaneNumber"][-1] == 120.3 and len(ron_c["counter"]) == 804
    assert ron_c["counter"][ron_c["octaneNumber"].index(90.0)] == 726  # cfr.ts CFR_GUIDE_COUNTER_PRF90.RON
    assert ron_c["counter"][ron_c["octaneNumber"].index(100.0)] == 919
    mon = {}
    for key, pages_c, pages_inv, pages_dial in [
        ("9/16", [32, 33], [34, 35], [36, 37]),
        ("19/32", [38, 39], [40, 41], [42, 43]),
        ("3/4", [44, 45], None, [46, 47]),
    ]:
        mon[key] = {"onToCounter": on_to_counter_table(D2700, pages_c, "counter"),
                    "onToDial": on_to_counter_table(D2700, pages_dial, "dial")}
        if pages_inv:
            mon[key]["counterToOn"] = counter_to_on_table(D2700, pages_inv)
    assert mon["9/16"]["onToCounter"]["counter"][mon["9/16"]["onToCounter"]["octaneNumber"].index(90.0)] == 749

    def add_cr(tbl):
        c = tbl["counter"]
        tbl["compressionRatioRigidRaise"] = [cr_rigid_raise(v) for v in c]
        tbl["compressionRatioChoiPolynomial"] = [cr_choi_polynomial(v) if 400 <= v <= 1400 else None for v in c]
        return tbl

    def add_cr_dial(tbl):
        c = [counter_from_dial(v) for v in tbl["dial"]]
        tbl["equivalentCounter"] = c
        tbl["compressionRatioRigidRaise"] = [cr_rigid_raise(v) for v in c]
        return tbl

    add_cr(ron_c)
    add_cr(ron_c2on)
    add_cr_dial(ron_dial)
    for v in mon.values():
        add_cr(v["onToCounter"])
        add_cr_dial(v["onToDial"])
        if "counterToOn" in v:
            add_cr(v["counterToOn"])
    # consistency: dial tables vs counter tables (dial = 1.012 - c/1410)
    dmax = max(abs(a - b) for a, b in zip(ron_dial["equivalentCounter"], ron_c["counter"]) if a is not None)

    # Errata of the published tables (the transcription is faithful; these are inconsistencies IN the standard):
    # dial entries whose equivalent counter (1.012 - dial) x 1410 differs from the counter table by > 4 counts, and
    # non-monotone counter entries.
    errata = []

    def dial_errata(label, cnt_tbl, dial_tbl):
        for on, eq, dial, c in zip(dial_tbl["octaneNumber"], dial_tbl["equivalentCounter"], dial_tbl["dial"],
                                   cnt_tbl["counter"]):
            if abs(eq - c) > 4:
                errata.append({"table": label, "octaneNumber": on, "publishedDialIn": dial,
                               "dialConsistentWithCounterTableIn": round(1.012 - c / 1410.0, 3),
                               "counterTable": c, "note": "dial entry inconsistent with the counter table (typo in the standard)"})

    def mono_errata(label, tbl):
        c = tbl["counter"]
        for i in range(1, len(c)):
            if c[i] < c[i - 1]:
                errata.append({"table": label, "octaneNumber": tbl["octaneNumber"][i], "publishedCounter": c[i],
                               "previousCounter": c[i - 1], "note": "counter decreases with increasing ON (typo in the standard)"})

    dial_errata("D2699 Table A4.3 (RON dial)", ron_c, ron_dial)
    mono_errata("D2699 Table A4.1 (RON counter)", ron_c)
    for key, lab in (("9/16", "A4.1/A4.3"), ("19/32", "A4.4/A4.6"), ("3/4", "A4.7/A4.8")):
        if len(mon[key]["onToCounter"]["counter"]) == len(mon[key]["onToDial"]["dial"]):
            dial_errata(f"D2700 Table {lab} (MON {key} dial)", mon[key]["onToCounter"], mon[key]["onToDial"])
        mono_errata(f"D2700 Table {lab.split('/')[0]} (MON {key} counter)", mon[key]["onToCounter"])

    baro_ron_below = baro_table(D2699, 37, True)
    baro_ron_above = baro_table(D2699, 38, True)
    baro_mon_below = baro_table(D2700, 48, False)
    baro_mon_above = [e for e in baro_mon_below if e["barometerInHg"] >= 30.0]
    baro_mon_below = [e for e in baro_mon_below if e["barometerInHg"] < 30.0]
    assert len(baro_ron_below) == 90 and len(baro_ron_above) == 10 and len(baro_mon_below) == 80 and len(baro_mon_above) == 10
    for e in baro_ron_below + baro_ron_above + baro_mon_below + baro_mon_above:
        e["barometerPa"] = e["barometerInHg"] * IN_HG
        if "intakeAirTemperatureC" in e:
            e["intakeAirTemperatureK"] = e["intakeAirTemperatureC"] + T0

    mon_spark = [(264, 0.825, 26), (337, 0.773, 25), (410, 0.721, 24), (484, 0.669, 23), (556, 0.617, 22),
                 (630, 0.565, 21), (704, 0.513, 20), (777, 0.461, 19), (851, 0.408, 18), (925, 0.356, 17),
                 (998, 0.304, 16), (1072, 0.252, 15), (1145, 0.200, 14)]
    # verify against the text layer of D2700 page 9 (Table 3)
    txt = " ".join(" ".join(r) for r in page_rows(D2700, 10))
    for c, dial, adv in mon_spark:
        assert f"{c} {dial:.3f} {adv}" in txt, (c, dial, adv)

    guide = {
        "description": (
            "ASTM D2699-15a (Research) and D2700-14 (Motor) guide tables of cylinder height (uncompensated digital "
            "counter / dial indicator) vs octane number at STANDARD KNOCK INTENSITY, 9/16 in venturi (MON also 19/32 and "
            "3/4 in), at standard barometric pressure 29.92 in Hg (ASTM writes '101.0 kPa'; 29.92 in Hg = 101.32 kPa), "
            "barometric compensation tables, MON spark-timing table, converted to compression ratio with the project "
            "relation (src/physics/engines/cfr.ts cfrCompressionRatioAtCounter, rigid raise 0.0007 in/digit) and the "
            "Choi et al. 2018 polynomial (ANL engine, valid 400-1400). Octane numbers > 100 are iso-octane + TEL blends "
            "(D2699 3.1.19.2, Table A3.4) - not representable with PRF chemistry."
        ),
        "generator": "tools/reference/cfr_data_astm.py",
        "sources": {
            "D2699": {"citation": "ASTM D2699-15a, Standard Test Method for Research Octane Number of Spark-Ignition Engine Fuel",
                      "url": URL_D2699, "localPdf": "tools/reference/cfr_sources/astm_D2699-15a.pdf"},
            "D2700": {"citation": "ASTM D2700-14, Standard Test Method for Motor Octane Number of Spark-Ignition Engine Fuel",
                      "url": URL_D2700, "localPdf": "tools/reference/cfr_sources/astm_D2700-14.pdf"},
        },
        "extraction": "table transcription from the PDF text layer (exact published digits); counter->CR computed",
        "sourceErrata": errata,
        "counterToCompressionRatio": {
            "rigidRaise": "CR = 1 + S/(h930 - 0.0007 in (c - 930)), S = 4.5 in, h930 = 0.01742116238087605 m "
                          "(cfr.ts CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER, fitted to Choi 2018 Fig. 4 oil-measured CR)",
            "choiPolynomial": "CR = 1.126e-8 c^3 - 2.126e-5 c^2 + 1.694e-2 c + 1.024 (Choi et al. 2018 Fig. 4, 400 <= c <= 1400)",
            "dialToCounter": "equivalent counter = (1.012 - dial [in]) x 1410 (D2699 Table A4.3 footnote A; D2700 Table A4.3)",
            "caveat": "The guide table fixes the cylinder HEIGHT; the CR at a given counter is engine-specific (indexed "
                      "by the motored compression pressure of Fig. 2, cfr_astm_compression_pressure.json). Published ANL "
                      "CR values differ between campaigns (e.g. RON 95 standard CR 7.26 in OSTI 2561394 vs 6.82 here).",
            "maxDialVsCounterTableDifferenceCounts": dmax,
        },
        "standardConditions": {
            "barometricPressureInHg": 29.92, "barometricPressurePa": 29.92 * IN_HG, "astmStatedKPa": 101.0,
            "RON": {"rpm": 600, "rpmTolerance": 6, "sparkAdvanceDegBTDC": 13,
                    "intakeAirTemperatureK": 52 + T0, "intakeAirTemperatureToleranceK": 1,
                    "intakeAirTemperatureNote": "at 29.92 in Hg; other barometers: Table A4.4/A4.5 IAT column",
                    "intakeMixtureTemperature": "not controlled", "coolantTemperatureK": 100 + T0,
                    "oilTemperatureK": 57 + T0, "oilPressurePa": [172e3, 207e3],
                    "humidityRatioKgPerKg": [0.00356, 0.00712],
                    "crankcaseGaugePressurePa": [-150 * MM_H2O, -25 * MM_H2O],
                    "exhaustBackPressureMaxGaugePa": 255 * MM_H2O,
                    "valveClearanceHotM": 0.008 * 0.0254, "sparkPlugGapM": 0.51e-3,
                    "fuelAirRatio": "maximum knock intensity (fuel level 0.7-1.7 in above venturi centreline)",
                    "clause": "D2699-15a 10.2-10.3"},
            "MON": {"rpm": 900, "rpmTolerance": 9, "sparkAdvanceDegBTDC": "Table 3 (linked to cylinder height)",
                    "intakeAirTemperatureK": 38 + T0, "intakeAirTemperatureToleranceK": 2.8,
                    "intakeMixtureTemperatureK": 149 + T0, "intakeMixtureTemperatureToleranceK": 1,
                    "intakeMixtureTemperatureTuningRangeK": [141 + T0, 163 + T0],
                    "coolantTemperatureK": 100 + T0, "oilTemperatureK": 57 + T0,
                    "humidityRatioKgPerKg": [0.00356, 0.00712],
                    "clause": "D2700-14 10.2-10.3"},
        },
        "RON_9_16": {"onToCounter": ron_c, "counterToOn": ron_c2on, "onToDial": ron_dial,
                     "tables": "D2699-15a Tables A4.1 (pp. 31-32), A4.2 (pp. 33-34), A4.3 (pp. 35-36)"},
        "MON": {k: dict(v, tables={"9/16": "D2700-14 Tables A4.1-A4.3 (pp. 32-37)",
                                   "19/32": "D2700-14 Tables A4.4-A4.6 (pp. 38-43)",
                                   "3/4": "D2700-14 Tables A4.7-A4.8 (pp. 44-47)"}[k]) for k, v in mon.items()},
        "barometricCompensation": {
            "RON": {"table": "D2699-15a Tables A4.4 (p. 37) and A4.5 (p. 38)",
                    "rule": "below 29.92 in Hg ADD the counter correction to the guide-table counter (higher CR), above "
                            "29.92 SUBTRACT; IAT column = required intake-air temperature at that barometer",
                    "below2992": baro_ron_below, "above2992": baro_ron_above},
            "MON": {"table": "D2700-14 Tables A4.9 and A4.10 (p. 48)",
                    "rule": "as RON; no IAT change for MON (IAT fixed 38 C, mixture 149 C)",
                    "below2992": baro_mon_below, "above2992": baro_mon_above},
        },
        "monSparkTiming": {
            "table": "D2700-14 Table 3 (p. 10), uncompensated cylinder height; linear, 1 deg per ~73.4 counts",
            "counter": [c for c, _, _ in mon_spark], "dialIn": [d for _, d, _ in mon_spark],
            "sparkAdvanceDegBTDC": [a for _, _, a in mon_spark],
            "compressionRatioRigidRaise": [cr_rigid_raise(c) for c, _, _ in mon_spark],
        },
        "monVenturi": {"table": "D2700-14 Table 2 (p. 8)",
                       "rows": [{"altitudeM": [0, 500], "venturiIn": 9 / 16, "barometerInHg": [28.0, 31.0]},
                                {"altitudeM": [500, 1000], "venturiIn": 19 / 32, "barometerInHg": [26.0, 29.0]},
                                {"altitudeM": [1000, None], "venturiIn": 3 / 4, "barometerInHg": [None, 27.0]}],
                       "RON": "9/16 in regardless of barometer (D2699-15a 10.2.6)"},
    }
    write_fixture("cfr_astm_guide_tables.json", guide)

    # ── motored compression pressure ──
    fig2 = fig2_compression_pressure()
    checks = [
        {"method": "RON", "clause": "D2699-15a A2.3, Table A2.2 (p. 25)", "octaneNumber": 93.4, "compensatedCounter": 778,
         "compensatedDialIn": 0.460, "psig": 169, "psigTol": 2, "mpaGaugeStated": 1.16, "mpaGaugeTol": 0.01},
        {"method": "RON", "clause": "D2699-15a A2.3, Table A2.2 (p. 25)", "octaneNumber": 105, "compensatedCounter": 1061,
         "compensatedDialIn": 0.259, "psig": 241, "psigTol": 4, "mpaGaugeStated": 1.66, "mpaGaugeTol": 0.02},
        {"method": "MON", "clause": "D2700-14 A2.3, Table A2.2 (p. 26)", "octaneNumber": 81.1, "compensatedCounter": 578,
         "compensatedDialIn": 0.602, "venturiIn": 9 / 16, "psig": 120, "psigTol": 2},
        {"method": "MON", "clause": "D2700-14 A2.3, Table A2.2 (p. 26)", "octaneNumber": 105, "compensatedCounter": 1008,
         "compensatedDialIn": 0.297, "venturiIn": 9 / 16, "psig": 194, "psigTol": 4},
    ]
    for c in checks:
        c["compressionRatioRigidRaise"] = cr_rigid_raise(c["compensatedCounter"])
        c["compressionRatioChoiPolynomial"] = cr_choi_polynomial(c["compensatedCounter"])
        c["gaugePa"] = c["psig"] * PSI
        c["gaugeTolPa"] = c["psigTol"] * PSI
        c["absPaAtStandardBarometer"] = c["psig"] * PSI + 29.92 * IN_HG
    comp = {
        "description": (
            "Motored ('compression') pressure specifications of the CFR engine: the peak motored cylinder pressure read "
            "on the ASTM compression-pressure gauge (check-valve peak-hold gauge screwed into the detonation-pickup hole, "
            "D2699 Fig. A2.4) with the engine hot, motored (ignition off, carburettor drained) at the method speed. "
            "Fig. 2: pressure that defines the basic cylinder height (counter 930 / dial 0.352 in) vs barometric pressure. "
            "Table A2.2: check values at two compensated cylinder heights. psig = gauge; absolute = gauge + barometer."
        ),
        "generator": "tools/reference/cfr_data_astm.py",
        "sources": guide["sources"],
        "conditions": {
            "RON": {"rpm": 600, "intakeAirTemperature": "per Table A4.4/A4.5 for the barometer (A2.2.3.1); A2.3: 51.7 +- 1 C",
                    "coolantTemperatureK": 100 + T0, "state": "hot engine, motored immediately after firing on a typical fuel"},
            "MON": {"rpm": 900, "intakeMixtureTemperatureK": 149 + T0, "intakeAirTemperatureK": 38 + T0,
                    "coolantTemperatureK": 100 + T0, "state": "hot engine, motored"},
            "note": "The gauge includes a hose and check valve (dead volume, D2699 Fig. A2.4) and reads the PEAK pressure "
                    "in the pickup hole, not the in-cylinder trace; treat as peak motored pressure +- tolerance. The "
                    "carburettor is drained, so the charge is AIR (no fuel), at the method intake temperature.",
        },
        "basicCylinderHeightVsBarometer": fig2,
        "checkPoints": checks,
    }
    write_fixture("cfr_astm_compression_pressure.json", comp)

    # ── knock intensity definitions + spread characteristic ──
    spread = [
        spread_figure(D2699, 155, "D2699-15a Fig. A2.6", 26, [66.5, 154.5, 249.5, 337.5, 426.0, 514.5],
                      [80, 85, 90, 95, 100, 105], [10.0, 99.0, 187.5, 276.5], [40, 30, 20, 10], 66.5, 360.0),
        spread_figure(D2700, 161, "D2700-14 Fig. A2.6", 27, [64.5, 154.5, 241.5, 333.5, 423.5, 514.0],
                      [80, 85, 90, 95, 100, 105], [63.5, 146.5, 234.5], [30, 20, 10], 64.5, 324.0),
    ]
    spread[0]["method"] = "RON"
    spread[1]["method"] = "MON"
    ki = {
        "description": "Definition and instrumentation of CFR 'standard knock intensity' (ASTM D2699-15a / D2700-14), "
                       "plus the published characteristics of the D-1 pickup / 501-C detonation meter chain.",
        "generator": "tools/reference/cfr_data_astm.py",
        "sources": dict(guide["sources"], Rockstroh2018={
            "citation": "Rockstroh T., Kolodziej C.P., Jespersen M.C., Goldsborough S.S., Wallner T., 'Insights into Engine "
                        "Knock: Comparison of Knock Metrics across Ranges of Intake Temperature and Pressure in the CFR "
                        "Engine', SAE Int. J. Fuels Lubr. 11(4), 2018, SAE 2018-01-0210 (peer-reviewed author version, DTU Orbit)",
            "url": "https://backend.orbit.dtu.dk/ws/files/236274325/KM_Investigation_CFR_Engine_SAE_Journal_Final.pdf"},
            Hoth2021={"citation": "Hoth A., Kolodziej C.P., 'Effects of knock intensity measurement technique and fuel chemical "
                                  "composition on the RON of FACE gasolines: Part 1', Fuel (2021), OSTI 1880351",
                      "url": "https://www.osti.gov/servlets/purl/1880351"}),
        "standardKnockIntensity": {
            "definition": "the knock level produced by a PRF of the guide-table octane number, at the fuel-air ratio of "
                          "MAXIMUM knock intensity, with the cylinder height at the guide-table value (compensated for "
                          "barometer); the detonation meter is then ADJUSTED so that the reading is 50 (D2699 3.1.24-3.1.25, "
                          "11.3, A2.4; D2700 3.1.28-3.1.29). It is therefore a relative scale: 50 = standard KI by construction.",
            "analogKnockmeterDivisions": 50, "analogTolerance": 2,
            "analogLinearRange": [20, 80], "digitalRange": [0, 999],
            "digitalPeakToPeakVoltsTypical": {"RON": 0.15, "MON": 0.25, "RONNote6Range": [0.05, 0.20]},
            "digitalDefaults": {"spread": 0, "timeConstant": {"RON": 25, "MON": 35}, "units": "meter setting (dimensionless)"},
            "analogZeroCheckTimeConstantSwitch": 3,
            "spreadSetting": {"RON": "12-15 K.I. divisions per O.N. at 90 O.N. (D2699 A2.5.2, 14.3.4)",
                              "MON": "approximately 12-14 at 90 O.N. (D2700 A2.5.2)"},
            "origin": "guide tables generated by setting the cylinder height to the former bouncing-pin value at 85 O.N. "
                      "and using that knock intensity as reference for PRF 40-100 (D2699 footnote 15; D2700 footnote 16)",
            "ratingInterpolation": "ON_s = ON_LRF + (KI_LRF - KI_s)/(KI_LRF - KI_HRF) (ON_HRF - ON_LRF) (bracketing, D2699 "
                                   "Procedure A; Hoth 2021 Eq. 2)",
        },
        "measurementChain": {
            "pickup": "D-1 detonation pickup: magnetostrictive transducer threaded into the cylinder head (pickup hole, "
                      "opposite the spark plug), exposed to chamber pressure; output voltage proportional to dp/dt "
                      "(D2699 3.1.6). Torque 30 lbf ft (Table A2.1).",
            "detonationMeter": "501-C: input filter shaped to 'simulate the octane rating characteristics of the original "
                               "bouncing-pin instrumentation'; threshold subtraction of the knock-free part (METER READING "
                               "control); amplification + pulse stretching (SPREAD control); integration over multiple "
                               "cycles (TIME CONSTANT); dc output to the knockmeter (D2699 A2.5.1, Fig. A2.5).",
            "inputFilterBandwidth": {"value": "low-pass, 'tuned for frequencies lower than 6.5 kHz' (from the 501-C circuit "
                                              "diagram per Swarts & Yates 2007, SAE 2007-01-0008, as reported by Rockstroh 2018 "
                                              "p. 2 and p. 8; Hoth 2021 Sec. 3.1: R-C input filter removes most knock "
                                              "oscillations)",
                                     "cutoffHz": 6500, "confidence": "secondary citation; circuit values not published openly"},
            "consequence": "knockmeter responds to the low-frequency pressure development after the knock point (autoignition "
                           "'knock point' pressure rise), not to the 6/10/14 kHz acoustic ringing; KU correlates poorly with "
                           "MAPO (Rockstroh 2018 Figs. 7-8; Hoth 2021 Fig. 12).",
            "dominantKnockFrequenciesHz": {"values": [6000, 10000, 14000], "source": "Rockstroh 2018 Fig. 4 (PRF90, standard RON "
                                           "knock, Kistler 6045AU20; ~2 orders of magnitude more energy at 6 kHz)"},
        },
        "spreadCharacteristic": spread,
    }
    write_fixture("cfr_astm_knock_intensity.json", ki)

    for e in fig2:
        print(e["figure"], e["method"], {k: round(v, 3) for k, v in e["linearFit"].items()},
              "p(29.92)=", round(np.interp(29.92, e["barometerInHg"], e["compressionPressurePsig"]), 2) if 29.92 <= max(e["barometerInHg"]) and 29.92 >= min(e["barometerInHg"]) else "")
    print("RON 90 ->", ron_c["counter"][ron_c["octaneNumber"].index(90.0)], "CR", round(cr_rigid_raise(726), 4))
    for e in guide["sourceErrata"]:
        print("erratum", e)


if __name__ == "__main__":
    main()
