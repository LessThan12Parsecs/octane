"""Oracle for the minimum geometric flow area of a poppet valve (valve-flow.ts).

Two independent evaluations for each (geometry, lift):
 1. Heywood (1988) §6.3.1 three-stage formulas, coded here independently of the TS code:
      stage 1  (w/(sinB cosB) > L > 0):  A = pi L cosB (Dv - 2w + L/2 sin 2B)
      stage 2:                            A = pi Dm sqrt((L - w tanB)^2 + w^2),  Dm = Dv - w
      stage 3  (port limited):            A = pi/4 (Dp^2 - Ds^2)
    (formulas and Fig. 6-12 geometry checked against the transcription at
    rgmracing.free.fr/luc/heywood1, fetched 2026-09-29), taking A = min(stage 1|2, stage 3).
 2. The exact minimum-area conical frustum between the two seat faces, found numerically:
    head seat from (r = Di/2, y = w tanB) to (Dv/2, 0); valve face = the same segment moved
    down by L. Frustum area pi (r1 + r2) |P1 - P2| minimised over both endpoints
    (scipy L-BFGS-B from a grid start). Heywood's stage-1 frustum (normal to the seat from the
    valve's inner edge) and stage-2 frustum (valve inner edge to seat outer edge) are
    candidates of this family, so exact <= Heywood; the gap measures the formulas' idealisation.
Writes test/fixtures/gasex_valve_flow.json.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402

IN = 0.0254

GEOMETRIES = [
    # CFR-like intake (sizes as in src/physics/engines/cfr.ts at the time of writing; inputs only)
    dict(name="cfr-intake-like", Dv=1.346 * IN + 0.125 * IN, Di=1.346 * IN, beta=math.radians(45),
         Ds=11 / 32 * IN, Dp=1.346 * IN),
    # automotive-like: 30 deg seat, narrow seat, small port (reaches stage 3)
    dict(name="auto-30deg", Dv=0.036, Di=0.0325, beta=math.radians(30), Ds=0.007, Dp=0.030),
    # wide 45 deg seat
    dict(name="wide-seat-45deg", Dv=0.040, Di=0.034, beta=math.radians(45), Ds=0.008, Dp=0.034),
]


def heywood(g, L):
    Dv, Di, b, Ds, Dp = g["Dv"], g["Di"], g["beta"], g["Ds"], g["Dp"]
    w = (Dv - Di) / 2
    Dm = Dv - w
    Aport = math.pi / 4 * (Dp * Dp - Ds * Ds)
    if L <= 0:
        return 0.0, 0
    if L < w / (math.sin(b) * math.cos(b)):
        A, st = math.pi * L * math.cos(b) * (Dv - 2 * w + L / 2 * math.sin(2 * b)), 1
    else:
        A, st = math.pi * Dm * math.sqrt((L - w * math.tan(b)) ** 2 + w * w), 2
    if A >= Aport:
        return Aport, 3
    return A, st


def exact_min_frustum(g, L):
    Dv, Di, b = g["Dv"], g["Di"], g["beta"]
    w = (Dv - Di) / 2
    h = w * math.tan(b)
    ri, ro = Di / 2, Dv / 2

    def pt(s, dy):  # s in [0,1] along seat from inner edge to outer edge
        return ri + s * (ro - ri), h * (1 - s) + dy

    def area(v):
        s1, s2 = v
        r1, y1 = pt(s1, 0.0)
        r2, y2 = pt(s2, -L)
        return math.pi * (r1 + r2) * math.hypot(r1 - r2, y1 - y2)

    best = None
    for s1 in np.linspace(0, 1, 11):
        for s2 in np.linspace(0, 1, 11):
            res = minimize(area, [s1, s2], bounds=[(0, 1), (0, 1)], method="L-BFGS-B",
                           options={"ftol": 1e-15, "gtol": 1e-12})
            if best is None or res.fun < best:
                best = res.fun
    return float(best)


def main():
    cases = []
    for g in GEOMETRIES:
        w = (g["Dv"] - g["Di"]) / 2
        rows = []
        for lod in [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4]:
            L = lod * g["Di"]
            A, st = heywood(g, L)
            rows.append({"lift": rnd(L), "areaHeywood": rnd(A), "stage": st,
                         "areaExactFrustum": rnd(exact_min_frustum(g, L)) if st < 3 else None})
        cases.append({"name": g["name"], "headDiameter": g["Dv"], "seatInnerDiameter": g["Di"],
                      "seatAngle": g["beta"], "stemDiameter": g["Ds"], "portDiameter": g["Dp"],
                      "seatWidth": rnd(w), "stage12Lift": rnd(w / (math.sin(g["beta"]) * math.cos(g["beta"]))),
                      "rows": rows})
    write_fixture("gasex_valve_flow.json", {"source": "tools/reference/gasex_valve_flow.py", "cases": cases})


if __name__ == "__main__":
    main()
