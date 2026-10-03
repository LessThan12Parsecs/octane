"""Digitise the measured CFR F1/F2 valve-lift profiles of Choi et al. (2018), Fig. 3.

Source: Choi, Kolodziej, Wallner & Hoth, "Development and Validation of a Three Pressure
Analysis (TPA) GT-Power Model of the CFR F1/F2 Engine for Estimating Cylinder Conditions",
SAE 2018-01-0848, author manuscript OSTI 1501884 (https://www.osti.gov/servlets/purl/1501884).
Fig. 3 = "Measured intake and exhaust valve lift profiles at valve lash '0'" (dial gauge on
the valve, ANL CFR engine whose cam timing satisfies ASTM D2699-15a A2.1.2: 0.054 in lifter
rise at 30 ± 2° — the digitised intake curve crosses 0.054 in at 30.3°).

Method (no manual clicking): the figure is VECTOR graphics in the PDF. The two curve paths
(solid = intake, dotted = exhaust; 44 cubic Béziers each) and the grid lines are read with
PyMuPDF; the axes are calibrated from the grid lines / tick labels (x: −360…360 CAD,
y: 0…0.3 in); each Bézier is sampled densely and the result is resampled on a 1° grid.
The figure's crank-angle axis has 0 at the GAS-EXCHANGE TDC (intake event around 0…240);
it is converted to the firing-TDC convention of DESIGN.md (intake: θ − 360, exhaust:
θ + 360, kept unwrapped so the exhaust table is continuous through 360°).

Cross-check printed by the script: crossings of the zero-lash profile at the 0.008 in hot
running clearance reproduce Choi et al. Table 2 (IVO/IVC 376/−152, EVO/EVC 141/373) to ≤ 2°.

Writes src/physics/engines/cfr-valve-lift.ts (GENERATED) and
test/fixtures/mechanics_cfr_valve_lift.json (crossing angles and peaks, for tests).

Requires PyMuPDF (`pip install pymupdf`), which is NOT part of .venv: if it (or the PDF)
is unavailable the script prints a notice and exits 0, leaving the committed outputs as they are.
Usage: python tools/reference/mechanics_cfr_valve_lift.py [path/to/osti_1501884.pdf]
"""
from __future__ import annotations

import json
import pathlib
import sys
import urllib.request

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT_TS = ROOT / "src" / "physics" / "engines" / "cfr-valve-lift.ts"
OUT_JSON = ROOT / "test" / "fixtures" / "mechanics_cfr_valve_lift.json"
URL = "https://www.osti.gov/servlets/purl/1501884"
IN = 0.0254


def load_pdf(argv: list[str]):
    try:
        import pymupdf  # noqa: F401
    except ImportError:
        print("mechanics_cfr_valve_lift.py: PyMuPDF not installed — skipped (committed outputs kept)")
        sys.exit(0)
    import pymupdf

    path = pathlib.Path(argv[1]) if len(argv) > 1 else pathlib.Path("/tmp/osti_1501884.pdf")
    if not path.exists():
        try:
            req = urllib.request.Request(URL, headers={"User-Agent": "Mozilla/5.0"})
            path.write_bytes(urllib.request.urlopen(req, timeout=60).read())
        except Exception as e:  # pragma: no cover - network
            print(f"mechanics_cfr_valve_lift.py: could not download {URL} ({e}) — skipped")
            sys.exit(0)
    return pymupdf.open(str(path)), pymupdf


def digitise(doc, pymupdf):
    page = doc[2]  # "Page 3 of 8"
    clip = pymupdf.Rect(30, 450, 250, 590)  # Fig. 3 on the left column
    drawings = [d for d in page.get_drawings() if clip.intersects(d["rect"])]
    # grid lines: horizontal (y = 0.05 … 0.3 in) and vertical (x = −270 … 360) thin strokes
    hy, vx = [], []
    for d in drawings:
        # thin solid strokes inside the plot area (excludes the outer figure frame and the legend box)
        inside = d["rect"].x0 > clip.x0 + 30 and d["rect"].y0 > clip.y0 + 5 and d["rect"].y1 < clip.y1 - 20
        if inside and d["type"] == "s" and d["width"] < 0.5 and not d.get("dashes", "").strip("[] 0"):
            for it in d["items"]:
                if it[0] != "l":
                    continue
                p, q = it[1], it[2]
                if abs(p.y - q.y) < 1e-3:
                    hy.append(p.y)
                elif abs(p.x - q.x) < 1e-3:
                    vx.append(p.x)
    hy = sorted(set(round(v, 3) for v in hy))
    vx = sorted(set(round(v, 3) for v in vx))
    # Least-squares axis calibration on all grid lines (Excel snaps each line to ±0.3 pt, so a
    # two-point calibration would carry ~1° of error).
    assert len(hy) == 7 and len(vx) == 8, (hy, vx)
    lift_ticks = np.arange(0.3, -0.01, -0.05)  # top grid line (smallest y) = 0.3 in … axis = 0
    ang_ticks = np.arange(-270.0, 361.0, 90.0)
    ky, cy = np.polyfit(lift_ticks, np.array(hy), 1)  # y_pt = ky·lift + cy
    kx, cx = np.polyfit(ang_ticks, np.array(vx), 1)  # x_pt = kx·θ + cx
    assert np.max(np.abs(np.polyval([ky, cy], lift_ticks) - hy)) < 0.5
    assert np.max(np.abs(np.polyval([kx, cx], ang_ticks) - vx)) < 0.5
    curves = [d for d in drawings if d["type"] == "s" and d["width"] > 1.0 and len(d["items"]) > 20]
    assert len(curves) == 2, "expected the two lift curves"

    def cal(x, y):
        return (x - cx) / kx, (y - cy) / ky

    out = {}
    for d in curves:
        dotted = bool(d.get("dashes", "").strip("[] 0"))
        pts = []
        for it in d["items"]:
            if it[0] == "l":
                pts += [cal(it[1].x, it[1].y), cal(it[2].x, it[2].y)]
            else:
                P = np.array([[q.x, q.y] for q in it[1:5]])
                for t in np.linspace(0, 1, 41):
                    b = (1 - t) ** 3 * P[0] + 3 * (1 - t) ** 2 * t * P[1] + 3 * (1 - t) * t**2 * P[2] + t**3 * P[3]
                    pts.append(cal(*b))
        pts = np.array(pts)
        pts = pts[np.argsort(pts[:, 0], kind="stable")]
        out["exhaust" if dotted else "intake"] = pts
    return out


def crossing(theta, lift, thr):
    idx = np.where(lift > thr)[0]
    a, b = idx[0], idx[-1]
    o = theta[a - 1] + (thr - lift[a - 1]) * (theta[a] - theta[a - 1]) / (lift[a] - lift[a - 1])
    c = theta[b] + (thr - lift[b]) * (theta[b + 1] - theta[b]) / (lift[b + 1] - lift[b])
    return float(o), float(c)


def main() -> None:
    doc, pymupdf = load_pdf(sys.argv)
    raw = digitise(doc, pymupdf)
    tables = {}
    fixture = dict(
        description="Choi et al. 2018 (SAE 2018-01-0848, OSTI 1501884) Fig. 3 zero-lash valve lift, digitised from the PDF vector paths.",
        generator="tools/reference/mechanics_cfr_valve_lift.py",
    )
    for name, shift in (("intake", -360.0), ("exhaust", 360.0)):
        pts = raw[name]
        th = pts[:, 0] + shift
        lift_in = np.clip(pts[:, 1], 0.0, None)
        # 1° grid over the support (first/last integer degree with the curve defined)
        g = np.arange(np.ceil(th[0]), np.floor(th[-1]) + 1.0)
        L = np.interp(g, th, lift_in)
        L[0] = 0.0
        L[-1] = 0.0
        tables[name] = (float(g[0]), L * IN)
        peak_i = int(np.argmax(L))
        rec = dict(startDeg=float(g[0]), stepDeg=1.0, peakLiftIn=float(L.max()), peakDeg=float(g[peak_i]))
        for thr in (0.003, 0.008, 0.054):
            o, c = crossing(g, L, thr)
            rec[f"cross_{thr}"] = [o, c]
        fixture[name] = rec
        print(name, rec)

    def fmt(a):
        return ",".join(f"{v:.4e}".replace("e-0", "e-").replace("e+0", "e") if v else "0" for v in a)

    ts = f"""/**
 * GENERATED by tools/reference/mechanics_cfr_valve_lift.py — do not edit by hand.
 *
 * Measured CFR F1/F2 valve-lift profiles at ZERO valve lash, Choi, Kolodziej, Wallner & Hoth,
 * SAE 2018-01-0848 (OSTI 1501884), Fig. 3, digitised from the PDF's vector paths (axes calibrated
 * on the grid lines). Tables: lift (m) on a 1° grid in the firing-TDC crank-angle convention
 * (intake around the gas-exchange TDC at −360°, exhaust unwrapped through +360°).
 * Peaks: intake {fixture['intake']['peakLiftIn']:.4f} in, exhaust {fixture['exhaust']['peakLiftIn']:.4f} in
 * (ASTM D2699-15a §10.2.4: resulting valve lift 0.238 ± 0.002 in).
 */

/** First tabulated crank angle of the intake table, deg (1° step). */
export const CHOI2018_INTAKE_LIFT_START_DEG = {tables['intake'][0]:.1f};
/** Intake valve lift at zero lash, m, at CHOI2018_INTAKE_LIFT_START_DEG + i (deg). */
export const CHOI2018_INTAKE_LIFT_ZERO_LASH: readonly number[] = [{fmt(tables['intake'][1])}];
/** First tabulated crank angle of the exhaust table, deg (1° step, unwrapped past 360°). */
export const CHOI2018_EXHAUST_LIFT_START_DEG = {tables['exhaust'][0]:.1f};
/** Exhaust valve lift at zero lash, m, at CHOI2018_EXHAUST_LIFT_START_DEG + i (deg). */
export const CHOI2018_EXHAUST_LIFT_ZERO_LASH: readonly number[] = [{fmt(tables['exhaust'][1])}];
"""
    OUT_TS.write_text(ts)
    OUT_JSON.write_text(json.dumps(fixture, separators=(",", ":")) + "\n")
    print(f"wrote {OUT_TS} and {OUT_JSON}")


if __name__ == "__main__":
    main()
