"""Published EXPERIMENTAL laminar burning velocities (iso-octane, n-heptane, methane, propane)
for src/physics/combustion/laminar-flame-speed.ts.

Source of the numbers: the CaltechMech validation reports of Blanquart et al. (Caltech, 2015),
https://www.theforce.caltech.edu/CaltechMech/reports/report_SL_{IC8,IC8_P,NC7,NC7_P,CH4,C3H8}.pdf,
which re-plot the original measurements (cited per point below) as VECTOR graphics. The points
are recovered exactly from the PDF drawing operators (marker centres mapped through the axis
calibration given by the tick marks and tick labels), so the only error is the plotting
precision of the report (≈ 0.05 pt ≈ 0.02-0.05 cm/s), not a manual read-off.

Mapping panel → conditions (from the report captions/annotations):
  SL_IC8, SL_NC7, SL_CH4, SL_C3H8 (S_L vs φ): each marker is assigned to the nearest of the
      mechanism curves drawn in the panel (IC8: T = 298 / 400 K; NC7: 298 / 353 / 400 K; CH4:
      298 K at 1 / 5 atm; C3H8: 298 K 1 bar, 343 K 1 bar, 298 K 5 bar); the dataset's own
      temperature is then used (Glaude 2012: 298 / 358 / 398 K; Kumar 2007: 400 K; Ji 2010 and
      Kelley 2011: 353 K; Veloo 2011: 343 K), all at 1 atm unless stated.
  SL_*_T08/T10/T13 (S_L vs T at φ = 0.8 / 1.0 / 1.3, 1 atm): T is the abscissa.
  SL_IC8_P*, SL_NC7_P* (S_L vs φ at 2-25 bar): Jerzembeck et al. 2009 (blue, 373 K, oxidiser
      X_O2 = 0.205) and Kelley et al. 2011 (red, 353 K). Pressures as labelled in the report
      ("bar"; Kelley et al. quote atm — a 1.3 % difference in p, ≈ 0.5 % in S_L).
  SL_CH4_p (S_L vs p at φ = 1, 298 K, log axis; pressure in atm).

Writes test/fixtures/flamespeed_experiments.json. Run:
  .venv/bin/python tools/reference/flamespeed_experiments.py [dir-with-cached-pdfs]
"""
from __future__ import annotations

import json
import math
import os
import re
import sys
import tempfile
import urllib.request
import zlib

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "test", "fixtures", "flamespeed_experiments.json")
BASE = "https://www.theforce.caltech.edu/CaltechMech/reports/report_SL_{}.pdf"
REPORTS = ["IC8", "IC8_P", "NC7", "NC7_P", "CH4", "C3H8"]
ATM = 101325.0

REFS = {
    "Davis 98": ("S.G. Davis, C.K. Law, Combust. Sci. Technol. 140 (1998) 427-449",
                 "counterflow twin flames, non-linear extrapolation"),
    "Huang 04": ("Y. Huang, C.J. Sung, J.A. Eng, Combust. Flame 139 (2004) 239-251",
                 "counterflow twin flames, linear extrapolation"),
    "Kumar 07": ("K. Kumar, J.E. Freeh, C.J. Sung, Y. Huang, J. Propul. Power 23 (2007) 428-436",
                 "counterflow twin flames, linear extrapolation (~2 cm/s above non-linear)"),
    "Ji 10": ("C. Ji, E. Dames, Y.L. Wang, H. Wang, F.N. Egolfopoulos, Combust. Flame 157 (2010) 277-287",
              "counterflow twin flames, non-linear extrapolation"),
    "Kelley 11": ("A.P. Kelley, A.J. Smallbone, D.L. Zhu, C.K. Law, Proc. Combust. Inst. 33 (2011) 963-970",
                  "outwardly propagating spherical flames"),
    "van Lipzig 11": ("J.P.J. van Lipzig, E.J.K. Nilsson, L.P.H. de Goey, A.A. Konnov, Fuel 90 (2011) 2773-2781",
                      "heat flux method"),
    "Glaude 12": ("P.A. Glaude, O. Herbinet, P. Dirrenberger, H. Le Gall, R. Bounaceur, F. Battin-Leclerc, "
                  "A. Pires da Cruz, A.A. Konnov, Proc. Combust. Inst. (2012) W2P028", "heat flux method"),
    "Jerzembeck 09": ("S. Jerzembeck, N. Peters, P. Pepiot-Desjardins, H. Pitsch, Combust. Flame 156 (2009) 292-301",
                      "outwardly propagating spherical flames, constant-volume bomb, air X_O2 = 0.205"),
    "Vagelopoulos 98": ("C.M. Vagelopoulos, F.N. Egolfopoulos, Proc. Combust. Inst. 27 (1998) 513-519",
                        "counterflow, direct (near-zero strain) determination"),
    "Hassan 98": ("M.I. Hassan, K.T. Aung, G.M. Faeth, Combust. Flame 115 (1998) 539-550",
                  "outwardly propagating spherical flames"),
    "Rozenchan 02": ("G. Rozenchan, D.L. Zhu, C.K. Law, S.D. Tse, Proc. Combust. Inst. 29 (2002) 1461-1469",
                     "outwardly propagating spherical flames"),
    "Bosschaart 04": ("K.J. Bosschaart, L.P.H. de Goey, Combust. Flame 136 (2004) 261-269", "heat flux method"),
    "Lowry 10": ("W. Lowry, J. de Vries, M. Krejci, E.L. Petersen, Z. Serinyel, W. Metcalfe, H. Curran, "
                 "G. Bourque, J. Eng. Gas Turb. Power 133 (2010) 019102", "spherical flames, constant-volume vessel"),
    "Vagelopoulos 94": ("C.M. Vagelopoulos, F.N. Egolfopoulos, C.K. Law, Proc. Combust. Inst. 25 (1994) 1341-1347",
                        "counterflow twin flames, extrapolation to zero strain"),
    "Jomaas 05": ("G. Jomaas, X.L. Zheng, D.L. Zhu, C.K. Law, Proc. Combust. Inst. 30 (2005) 193-200",
                  "outwardly propagating spherical flames"),
    "Veloo 11": ("P.S. Veloo, F.N. Egolfopoulos, Combust. Flame 158 (2011) 501-510", "counterflow, non-linear extrapolation"),
}
REFS["Lowri 10"] = REFS["Lowry 10"]  # label typo in the CH4 report

# ---------------------------------------------------------------------------------------
# minimal vector-PDF interpreter (Form XObjects written by pdfTeX from matplotlib/gnuplot)
# ---------------------------------------------------------------------------------------
NUM = r"-?\d*\.?\d+(?:[eE][-+]?\d+)?"


def form_streams(data: bytes) -> dict[str, str]:
    objs = {int(m.group(1)): m.group(2) for m in re.finditer(rb"(\d+) 0 obj(.*?)endobj", data, re.S)}
    out = {}
    for _, body in sorted(objs.items()):
        m = re.match(rb"\s*<<(.*?)>>\s*stream\r?\n", body, re.S)
        if not m or b"/Subtype /Form" not in m.group(1):
            continue
        name = re.search(rb"PTEX.FileName \((.*?)\)", m.group(1))
        raw = body[m.end():]
        raw = raw[: raw.rfind(b"endstream")]
        key = os.path.basename(name.group(1).decode()) if name else str(len(out))
        out[key] = zlib.decompress(raw).decode("latin1")
    return out


def mul(a, b):
    return [a[0] * b[0] + a[1] * b[2], a[0] * b[1] + a[1] * b[3], a[2] * b[0] + a[3] * b[2],
            a[2] * b[1] + a[3] * b[3], a[4] * b[0] + a[5] * b[2] + b[4], a[4] * b[1] + a[5] * b[3] + b[5]]


def interpret(s: str):
    toks = re.findall(r"\((?:[^()\\]|\\.)*\)|/[^\s/\[\]()<>]+|\[|\]|" + NUM + r"|[A-Za-z*'\"]+", s)
    st, ctm, stack = [], [1, 0, 0, 1, 0, 0], []
    fill = stroke = (0.0, 0.0, 0.0)
    path, cur, markers, lines, texts = [], [], [], [], []
    tmat = lm = [1, 0, 0, 1, 0, 0]
    seq = 0

    def tr(x, y):
        return (ctm[0] * x + ctm[2] * y + ctm[4], ctm[1] * x + ctm[3] * y + ctm[5])

    for t in toks:
        seq += 1
        if re.fullmatch(NUM, t):
            st.append(float(t))
            continue
        if t[0] in "(/":
            st.append(t)
            continue
        op = t
        if op == "q":
            stack.append((ctm[:], fill, stroke))
        elif op == "Q":
            ctm, fill, stroke = stack.pop()
        elif op == "cm":
            ctm = mul(st[-6:], ctm)
        elif op == "rg":
            fill = tuple(st[-3:])
        elif op == "g":
            fill = (st[-1],) * 3
        elif op == "RG":
            stroke = tuple(st[-3:])
        elif op == "G":
            stroke = (st[-1],) * 3
        elif op == "m":
            if cur:
                path.append(cur)
            cur = [tr(st[-2], st[-1])]
        elif op == "l":
            cur.append(tr(st[-2], st[-1]))
        elif op == "c":
            cur += [tr(st[-2], st[-1]), tr(st[-4], st[-3]), tr(st[-6], st[-5])]
        elif op == "re":
            x, y, w, h = st[-4:]
            if cur:
                path.append(cur)
            path.append([tr(x, y), tr(x + w, y), tr(x + w, y + h), tr(x, y + h)])
            cur = []
        elif op in ("f", "F", "f*", "S", "s", "B", "b", "n"):
            if cur:
                path.append(cur)
                cur = []
            for pth in path:
                xs = [p[0] for p in pth]
                ys = [p[1] for p in pth]
                w, h = max(xs) - min(xs), max(ys) - min(ys)
                filled = op in ("f", "F", "f*", "B", "b")
                rec = dict(i=seq, x=(max(xs) + min(xs)) / 2, y=(max(ys) + min(ys)) / 2, n=len(pth),
                           op=op, color=tuple(round(c, 2) for c in (fill if filled else stroke)))
                if w < 8 and h < 8 and len(pth) >= 3:
                    markers.append(rec)
                elif w < 8 and h < 8 and len(pth) == 2 and not filled:
                    last = markers[-1] if markers else None
                    if last and last.get("strokes") and abs(last["x"] - rec["x"]) < 0.3 and \
                            abs(last["y"] - rec["y"]) < 0.3 and last["color"] == rec["color"]:
                        last["strokes"] += 1
                    else:
                        rec["strokes"] = 1
                        rec["n"] = "x"
                        markers.append(rec)
                    lines.append(dict(pts=pth))
                elif not filled:
                    lines.append(dict(pts=pth))
            path = []
        elif op == "BT":
            tmat = lm = [1, 0, 0, 1, 0, 0]
        elif op == "Tm":
            tmat = lm = st[-6:]
        elif op == "Td":
            lm = mul([1, 0, 0, 1, st[-2], st[-1]], lm)
            tmat = lm[:]
        elif op == "Tj":
            full = mul(tmat, ctm)
            texts.append(dict(i=seq, t=st[-1][1:-1], x=full[4], y=full[5]))
        st = []
    markers = [m for m in markers if m.get("strokes", 2) >= 2]
    return markers, lines, texts


def isnum(s: str) -> bool:
    try:
        float(s.strip())
        return True
    except ValueError:
        return False


class Panel:
    """Axis calibration from the frame, the major tick marks and the numeric tick labels."""

    def __init__(self, s: str, logx: bool = False):
        self.markers, self.lines, self.texts = interpret(s)
        frames = [l["pts"] for l in self.lines if len(l["pts"]) == 4]
        fr = max(frames, key=lambda p: (max(q[0] for q in p) - min(q[0] for q in p)) *
                 (max(q[1] for q in p) - min(q[1] for q in p)))
        self.x0, self.x1 = min(q[0] for q in fr), max(q[0] for q in fr)
        self.y0, self.y1 = min(q[1] for q in fr), max(q[1] for q in fr)
        two = [l["pts"] for l in self.lines if len(l["pts"]) == 2]
        L = lambda p: math.hypot(p[1][0] - p[0][0], p[1][1] - p[0][1])  # noqa: E731
        xt_all = [p for p in two if abs(p[0][1] - self.y0) < 0.01 and abs(p[0][0] - p[1][0]) < 0.01]
        yt_all = [p for p in two if abs(p[0][0] - self.x0) < 0.01 and abs(p[0][1] - p[1][1]) < 0.01]
        mx, my = max(map(L, xt_all)), max(map(L, yt_all))
        xt = sorted({round(p[0][0], 3) for p in xt_all if L(p) > 0.8 * mx})
        yt = sorted({round(p[0][1], 3) for p in yt_all if L(p) > 0.8 * my})
        xv = [float(t["t"]) for t in sorted((t for t in self.texts if isnum(t["t"]) and t["y"] < self.y0 - 8),
                                             key=lambda t: t["x"])]
        yv = [float(t["t"]) for t in sorted((t for t in self.texts if isnum(t["t"]) and t["x"] < self.x0
                                              and t["y"] >= self.y0 - 8), key=lambda t: t["y"])]
        if len(xt) != len(xv) or len(yt) != len(yv):
            raise ValueError("tick/label mismatch")
        self.logx = logx
        self.ax, self.bx = self._lsq(xt, [math.log(v) for v in xv] if logx else xv)
        self.ay, self.by = self._lsq(yt, yv)

    @staticmethod
    def _lsq(t, v):
        n = len(t)
        mt, mv = sum(t) / n, sum(v) / n
        a = sum((ti - mt) * (vi - mv) for ti, vi in zip(t, v)) / sum((ti - mt) ** 2 for ti in t)
        b = mv - a * mt
        if max(abs(a * ti + b - vi) for ti, vi in zip(t, v)) > 0.003 * (max(v) - min(v)):
            raise ValueError("non-linear axis")
        return a, b

    def X(self, x):
        v = self.ax * x + self.bx
        return math.exp(v) if self.logx else v

    def Y(self, y):
        return self.ay * y + self.by

    def inside(self, m):
        return self.x0 - 1 <= m["x"] <= self.x1 + 1 and self.y0 - 1 <= m["y"] <= self.y1 + 1

    def curves(self):
        """Long stroked polylines (mechanism results) as lists of (x, y) in data units."""
        return [[(self.X(x), self.Y(y)) for x, y in l["pts"]] for l in self.lines if len(l["pts"]) > 15]


def style(m):
    return (m["color"], m["n"], m["op"])


def legend_map(panel: Panel, names) -> dict:
    out = {}
    for t in panel.texts:
        if t["t"] not in names:
            continue
        best = None
        for m in panel.markers:
            dx, dy = m["x"] - t["x"], m["y"] - (t["y"] + 2.5)
            if 0 < dx < 90 and abs(dy) < 2.5 and (best is None or abs(dy) < best[0]):
                best = (abs(dy), m)
        if best:
            out[style(best[1])] = (t["t"], best[1]["x"], best[1]["y"])
    return out


def dataset_points(panel: Panel, legend: dict):
    pts = []
    for m in panel.markers:
        if not panel.inside(m):
            continue
        entry = legend.get(style(m))
        if entry is None:
            continue
        name, lx, ly = entry
        if abs(m["x"] - lx) < 1 and abs(m["y"] - ly) < 1:
            continue  # the legend sample itself
        pts.append((name, panel.X(m["x"]), panel.Y(m["y"])))
    return pts


def nearest_curve(curves, x, y):
    """Index of the curve closest to (x, y) in log(S_L) at abscissa x (curves are extended
    with their end values beyond their abscissa range)."""
    best, k = None, -1
    for i, c in enumerate(curves):
        c = sorted(c)
        if x <= c[0][0]:
            yc = c[0][1]
        elif x >= c[-1][0]:
            yc = c[-1][1]
        else:
            yc = next(ya + (yb - ya) * (x - xa) / (xb - xa)
                      for (xa, ya), (xb, yb) in zip(c, c[1:]) if xa <= x <= xb and xb > xa)
        d = abs(math.log(max(y, 1e-3) / max(yc, 1e-3)))
        if best is None or d < best:
            best, k = d, i
    return k


# ---------------------------------------------------------------------------------------
# per-report mapping to (fuel, phi, T, p)
# ---------------------------------------------------------------------------------------
def collect(pdfs: dict[str, bytes]):
    rows = []

    def add(fuel, name, phi, T, p, SL, panel):
        ref, method = REFS[name]
        key = "Lowry 10" if name == "Lowri 10" else name
        rows.append(dict(fuel=fuel, phi=round(phi, 3), Tu=round(T, 1), p=round(p, 1), SL=round(SL / 100, 5),
                         dataset=key, ref=ref, method=method, digitizedFrom=panel))

    for fuel, rep, curveT, dsT in (
        ("IC8H18", "IC8", [298, 400], {"Kumar 07": {400: 400}, "Glaude 12": {400: 398}}),
        ("NC7H16", "NC7", [298, 353, 400], {"Kumar 07": {400: 400}, "Glaude 12": {353: 358, 400: 398},
                                            "Ji 10": {353: 353}, "Kelley 11": {353: 353}}),
    ):
        forms = form_streams(pdfs[rep])
        main = Panel(forms[f"SL_{rep}.pdf"])
        names = [t["t"] for t in main.texts if t["t"] in REFS]
        leg = legend_map(main, names)
        curves = main.curves()
        order = sorted(range(len(curves)), key=lambda i: sum(y for _, y in curves[i]) / len(curves[i]))
        # φ-sweep points (except φ = 0.8/1.0/1.3, which come from the T-sweep panels)
        for name, phi, SL in dataset_points(main, leg):
            if any(abs(phi - q) < 0.01 for q in (0.8, 1.0, 1.3)):
                continue
            if name not in dsT or name in ("Ji 10", "Kelley 11"):
                T = dsT[name][353] if name in dsT else 298  # single-temperature datasets
            else:
                k = nearest_curve(curves, phi, SL)
                Tc = curveT[order.index(k)]
                T = dsT[name].get(Tc, 298 if Tc == 298 else None)
                if T is None:
                    raise ValueError(f"{rep}: {name} point at φ={phi:.2f} near the {Tc} K curve")
            add(fuel, name, phi, T, ATM, SL, f"report_SL_{rep}.pdf / SL_{rep}")
        for tag, phi in (("T08", 0.8), ("T10", 1.0), ("T13", 1.3)):
            pan = Panel(forms[f"SL_{rep}_{tag}.pdf"])
            for name, T, SL in dataset_points(pan, leg):
                add(fuel, name, phi, round(T), ATM, SL, f"report_SL_{rep}.pdf / SL_{rep}_{tag}")

    for fuel, rep in (("IC8H18", "IC8_P"), ("NC7H16", "NC7_P")):
        for key, s in form_streams(pdfs[rep]).items():
            pan = Panel(s)
            pbar = float(re.search(r"P(\d+)", key).group(1))
            for m in pan.markers:
                if not pan.inside(m):
                    continue
                if m["color"] == (0.0, 0.0, 1.0):
                    name, T = "Jerzembeck 09", 373
                elif m["color"] == (1.0, 0.0, 0.0):
                    name, T = "Kelley 11", 353
                else:
                    continue
                add(fuel, name, pan.X(m["x"]), T, pbar * 1e5, pan.Y(m["y"]), f"report_SL_{rep}.pdf / {key}")

    # methane: φ sweep at 1 atm and 5 atm (298 K), pressure sweep at φ = 1
    forms = form_streams(pdfs["CH4"])
    main = Panel(forms["SL_CH4.pdf"])
    leg = legend_map(main, [t["t"] for t in main.texts if t["t"] in REFS])
    curves = main.curves()
    order = sorted(range(len(curves)), key=lambda i: sum(y for _, y in curves[i]) / len(curves[i]))
    for name, phi, SL in dataset_points(main, leg):
        p = ATM  # single-pressure datasets; Rozenchan 02 and Lowry 10 also measured at 5 atm
        if name in ("Rozenchan 02", "Lowri 10", "Lowry 10"):
            p = 5 * ATM if order.index(nearest_curve(curves, phi, SL)) == 0 else ATM
        add("CH4", name, phi, 298, p, SL, "report_SL_CH4.pdf / SL_CH4")
    pan = Panel(forms["SL_CH4_p.pdf"], logx=True)
    for name, p, SL in dataset_points(pan, leg):
        if abs(p - 1) < 0.01 or abs(p - 5) < 0.05:
            continue  # already in the φ sweep
        add("CH4", name, 1.0, 298, p * ATM, SL, "report_SL_CH4.pdf / SL_CH4_p")

    # propane: 298 K 1 bar, 343 K 1 bar, 298 K 5 bar
    forms = form_streams(pdfs["C3H8"])
    main = Panel(forms["SL_C3H8.pdf"])
    leg = legend_map(main, [t["t"] for t in main.texts if t["t"] in REFS])
    curves = main.curves()
    order = sorted(range(len(curves)), key=lambda i: sum(y for _, y in curves[i]) / len(curves[i]))
    for name, phi, SL in dataset_points(main, leg):
        T, p = (343, 1e5) if name == "Veloo 11" else (298, 1e5)
        if name in ("Jomaas 05", "Lowry 10"):  # measured at 1 and 5 bar
            k = order.index(nearest_curve(curves, phi, SL))
            T, p = [(298, 5e5), (298, 1e5), (298, 1e5)][k]
        add("C3H8", name, phi, T, p, SL, "report_SL_C3H8.pdf / SL_C3H8")
    return rows


def main() -> None:
    cache = sys.argv[1] if len(sys.argv) > 1 else None
    pdfs = {}
    try:
        with tempfile.TemporaryDirectory() as tmp:
            for rep in REPORTS:
                path = os.path.join(cache or tmp, f"report_SL_{rep}.pdf")
                if not os.path.exists(path):
                    with urllib.request.urlopen(BASE.format(rep), timeout=60) as r, open(path, "wb") as fh:
                        fh.write(r.read())
                with open(path, "rb") as fh:
                    pdfs[rep] = fh.read()
    except OSError as err:
        if os.path.exists(OUT):
            print(f"flamespeed_experiments: could not fetch the reports ({err}); keeping the committed fixture")
            return
        raise
    rows = collect(pdfs)
    rows.sort(key=lambda r: (r["fuel"], r["p"], r["Tu"], r["dataset"], r["phi"]))
    data = dict(
        source="Experimental laminar burning velocities digitised (exactly, from vector graphics) from the "
               "CaltechMech validation reports (Blanquart et al., Caltech, 2015); original references per point. "
               "See tools/reference/flamespeed_experiments.py.",
        units=dict(SL="m/s", Tu="K", p="Pa"),
        points=rows,
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=None, separators=(",", ":"))
        fh.write("\n")
    by = {}
    for r in rows:
        by.setdefault((r["fuel"], r["dataset"]), 0)
        by[(r["fuel"], r["dataset"])] += 1
    for k, v in sorted(by.items()):
        print(f"  {k[0]:7s} {k[1]:16s} {v}")
    print(f"wrote {os.path.relpath(OUT, ROOT)} ({len(rows)} points)")


if __name__ == "__main__":
    main()
