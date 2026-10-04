"""Oracle for the side-valve (L-head) flow stages of valve-flow.ts.

For each (valve geometry, pocket) case and lift:
 1. Heywood (1988) §6.3.1 seat/port stages (gasex_valve_flow.heywood, independent of the TS code),
    then the roof-masking stage A_roof = pi D_v (h_c - L) entering the minimum:
      A_m = min(A_heywood(L), A_roof(L)),  stage 4 when the roof is the minimum.
 2. The pocket/bore boundary length: the arc of the bore circle inside the rounded-rectangle
    pocket plan, found from the roots of the plan's SIGNED DISTANCE FUNCTION along the circle
    (dense sampling + brentq) — a different construction from the TS clamp/bisection test.
 3. The transfer area A_t = w_t * max(0, h_p - max(0, crown height above deck)) for a set of crown
    heights, and the series combination 1/(CdA)^2 = 1/(CdA_v)^2 + 1/(Cd_t A_t)^2 for given CdA_v.
Writes test/fixtures/gasex_lhead_flow.json.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gasex_valve_flow import heywood  # noqa: E402
from thermo_common import rnd, write_fixture  # noqa: E402

IN = 0.0254

# Model T-like inputs (src/physics/engines/model-t.ts at the time of writing; inputs only): valve
# 1.47 in head, 1-5/16 in throat, 45 deg, 0.311 in stem, 1-1/8 in manifold port; 3.75 in bore; pocket
# plan as MODEL_T; roof 12.92 mm above the deck.
_VD = math.pi / 4 * (3.75 * IN) ** 2 * 4 * IN
_VC = _VD / (3.98 - 1)
_POCKET_H = (_VC - 2.5e-6 - math.pi / 4 * (3.75 * IN) ** 2 * 1.0 * IN) / 46.03020157589403e-4

CASES = [
    dict(name="model-t-like", Dv=1.47 * IN, Di=1.3125 * IN, beta=math.radians(45), Ds=0.311 * IN, Dp=1.125 * IN,
         bore=3.75 * IN, pocket=dict(xMin=-0.0915, xMax=-0.03, zMin=-0.0475, zMax=0.0475, cornerRadius=0.02),
         roofClearance=_POCKET_H, pocketHeight=_POCKET_H, crownAtTDC=5 / 16 * IN),
    # low roof (binds at high lift), sharp-cornered pocket crossing the straight edge x = xMax
    dict(name="low-roof", Dv=0.036, Di=0.032, beta=math.radians(45), Ds=0.008, Dp=0.032,
         bore=0.080, pocket=dict(xMin=-0.075, xMax=-0.025, zMin=-0.030, zMax=0.030, cornerRadius=0.004),
         roofClearance=0.006, pocketHeight=0.007, crownAtTDC=0.003),
]


def sdf_rounded_rect(p, x, z):
    rc = min(p["cornerRadius"], 0.5 * (p["xMax"] - p["xMin"]), 0.5 * (p["zMax"] - p["zMin"]))
    cx = 0.5 * (p["xMin"] + p["xMax"])
    cz = 0.5 * (p["zMin"] + p["zMax"])
    hx = 0.5 * (p["xMax"] - p["xMin"]) - rc
    hz = 0.5 * (p["zMax"] - p["zMin"]) - rc
    qx = abs(x - cx) - hx
    qz = abs(z - cz) - hz
    outside = math.hypot(max(qx, 0.0), max(qz, 0.0))
    return outside + min(max(qx, qz), 0.0) - rc


def boundary_length(c):
    R = c["bore"] / 2
    f = lambda a: sdf_rounded_rect(c["pocket"], R * math.cos(a), R * math.sin(a))  # noqa: E731
    a = np.linspace(0, 2 * math.pi, 200001)
    v = np.array([f(t) for t in a])
    roots = []
    for i in range(len(a) - 1):
        if v[i] == 0.0:
            roots.append(a[i])
        elif v[i] * v[i + 1] < 0:
            roots.append(brentq(f, a[i], a[i + 1], xtol=1e-15))
    # integrate the inside (sdf < 0) arcs between consecutive roots
    total = 0.0
    pts = [0.0] + roots + [2 * math.pi]
    for lo, hi in zip(pts[:-1], pts[1:]):
        if hi > lo and f(0.5 * (lo + hi)) < 0:
            total += hi - lo
    return total * R


def main():
    out = []
    for c in CASES:
        g = dict(Dv=c["Dv"], Di=c["Di"], beta=c["beta"], Ds=c["Ds"], Dp=c["Dp"])
        rows = []
        for L in np.concatenate([np.linspace(0.0002, 0.9 * c["roofClearance"], 24), [c["roofClearance"] * 0.999]]):
            L = float(L)
            A, st = heywood(g, L)
            Ar = math.pi * c["Dv"] * (c["roofClearance"] - L)
            if Ar < A:
                A, st = Ar, 4
            rows.append({"lift": rnd(L), "area": rnd(A), "stage": st, "roofArea": rnd(Ar)})
        w = boundary_length(c)
        tr = []
        for crown in [-0.02, 0.0, 0.5 * c["crownAtTDC"], c["crownAtTDC"], c["pocketHeight"] + 0.001]:
            At = w * max(0.0, c["pocketHeight"] - max(0.0, crown))
            series = []
            for cdav, cdt in [(1e-4, 1.0), (3e-4, math.pi / (math.pi + 2)), (5e-5, 0.8)]:
                cdat = cdt * At
                series.append({"cdaValve": cdav, "cdT": rnd(cdt),
                               "cdaSeries": rnd(0.0 if cdat <= 0 else (cdav ** -2 + cdat ** -2) ** -0.5)})
            tr.append({"crownAboveDeck": rnd(crown), "transferArea": rnd(At), "series": series})
        out.append({"name": c["name"], "headDiameter": c["Dv"], "seatInnerDiameter": c["Di"], "seatAngle": c["beta"],
                    "stemDiameter": c["Ds"], "portDiameter": c["Dp"], "bore": c["bore"], "pocket": c["pocket"],
                    "roofClearance": rnd(c["roofClearance"]), "pocketHeight": rnd(c["pocketHeight"]),
                    "boundaryLength": rnd(w), "rows": rows, "transfer": tr})
    write_fixture("gasex_lhead_flow.json", {"source": "tools/reference/gasex_lhead_flow.py", "cases": out})


if __name__ == "__main__":
    main()
