"""Oracle for the three-arc cam on a flat-faced translating follower (gas-exchange/cam-lift.ts).

Independent of the closed-form lift equations of the TS code:
 1. The arc layout is found numerically: the flank-circle centre C_f is the point (below the nose
    axis) with |C_f| = rho - r_b (internal tangency to the base circle) and |N - C_f| = rho - r_n
    (internal tangency to the nose circle), N = (r_b + rise - r_n, 0); scipy.optimize.fsolve.
    The tangent points (junctions) are the points where the circles touch.
 2. The lobe outline is POLYGONISED: base arc, two flank arcs and the nose arc are sampled densely
    between their tangent points (about 4e5 vertices).
 3. The follower lift of a flat face perpendicular to the follower axis at cam angle phi is the
    support function of the outline, max_i p_i . u(phi), minus the base radius; computed by brute
    force over the vertices. Seat-to-seat durations for a lash are found by brentq on that envelope.
 4. The contact-point offset (distance of the touching vertex from the follower axis; resolved to the
    vertex spacing, rho/60000 on the flanks) gives the smallest admissible follower-face radius.
Polygon error <= r (1 - cos(dpsi/2)) ~ 1e-10 m for the vertex spacing used.
Cases: the stock Model T lobe (MTFC Tulsa design_stock.htm Table 1, inch values), and a generic
modern-proportion three-arc lobe (inputs only).
Writes test/fixtures/mechanics_cam_three_arc.json.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.optimize import brentq, fsolve

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402

IN = 0.0254

CASES = [
    # MTFC Tulsa 'Stock 1912+ Model T Camshaft', Table 1 (Ford drawing via Post 1997 / Turkish 1946)
    dict(name="model-t-mtfc", baseRadius=0.4060 * IN, flankRadius=1.2601 * IN, noseRadius=0.0313 * IN,
         rise=0.2502 * IN, lashes=[0.0, 0.010 * IN, 0.015625 * IN, 0.025 * IN, 0.0256 * IN, 0.03125 * IN,
                                   0.050 * IN, 0.065 * IN]),
    # generic proportions (inputs only): 15 mm base, 60 mm flanks, 4 mm nose, 8 mm rise
    dict(name="generic", baseRadius=0.015, flankRadius=0.060, noseRadius=0.004, rise=0.008,
         lashes=[0.0, 0.0002, 0.0005, 0.002]),
]


def layout(c):
    rb, rf, rn, rise = c["baseRadius"], c["flankRadius"], c["noseRadius"], c["rise"]
    N = np.array([rb + rise - rn, 0.0])

    def eqs(v):
        x, y = v
        return [math.hypot(x, y) - (rf - rb), math.hypot(N[0] - x, N[1] - y) - (rf - rn)]

    # flank for phi > 0 (upper flank): centre below the axis, behind the cam axis
    Cf = np.array(fsolve(eqs, [-0.3 * (rf - rb), -0.9 * (rf - rb)], xtol=1e-13))
    assert Cf[1] < 0
    ub = -Cf / np.linalg.norm(Cf)  # base tangent direction (from O through T, away from C_f)
    un = (N - Cf) / np.linalg.norm(N - Cf)
    phi_b = math.atan2(ub[1], ub[0])
    phi_n = math.atan2(un[1], un[0])
    return N, Cf, phi_b, phi_n


def outline(c, N, Cf, phi_b, phi_n, n_per_rad=60000):
    rb, rf, rn = c["baseRadius"], c["flankRadius"], c["noseRadius"]
    pts = []

    def arc(center, r, a0, a1):
        m = max(16, int(abs(a1 - a0) * n_per_rad))
        a = np.linspace(a0, a1, m)
        pts.append(np.stack([center[0] + r * np.cos(a), center[1] + r * np.sin(a)], axis=1))

    arc(N, rn, -phi_n, phi_n)  # nose
    arc(Cf, rf, phi_n, phi_b)  # upper flank (normal angle phi_n..phi_b)
    Cf_low = np.array([Cf[0], -Cf[1]])
    arc(Cf_low, rf, -phi_b, -phi_n)  # lower flank
    arc(np.zeros(2), rb, phi_b, 2 * math.pi - phi_b)  # base circle
    return np.concatenate(pts)


def support(P, phi):
    u = np.array([math.cos(phi), math.sin(phi)])
    d = P @ u
    i = int(np.argmax(d))
    # contact offset = component of the touching point perpendicular to u (signed, toward +phi)
    off = float(P[i] @ np.array([-math.sin(phi), math.cos(phi)]))
    return float(d[i]), off


def main():
    out = []
    for c in CASES:
        N, Cf, phi_b, phi_n = layout(c)
        P = outline(c, N, Cf, phi_b, phi_n)
        rb = c["baseRadius"]
        rows = []
        max_off = 0.0
        for cam_deg in np.concatenate([np.arange(0.0, 90.0, 0.5), [math.degrees(phi_n), math.degrees(phi_b)]]):
            phi = math.radians(float(cam_deg))
            h, off = support(P, phi)
            rows.append({"camDeg": rnd(float(cam_deg)), "camLift": rnd(max(h - rb, 0.0)), "contactOffset": rnd(off)})
            max_off = max(max_off, abs(off))
        # finer search for the max contact offset near the junction
        for phi in np.linspace(phi_n - 0.01, phi_n + 0.01, 401):
            max_off = max(max_off, abs(support(P, float(phi))[1]))
        durs = []
        for lash in c["lashes"]:
            if lash == 0:
                phi_s = phi_b
            else:
                phi_s = brentq(lambda p: support(P, p)[0] - rb - lash, 0.0, phi_b + 0.05, xtol=1e-14)
            durs.append({"lash": rnd(lash), "durationCrankDeg": rnd(4 * math.degrees(phi_s))})
        out.append({
            "name": c["name"], "baseRadius": c["baseRadius"], "flankRadius": c["flankRadius"],
            "noseRadius": c["noseRadius"], "rise": c["rise"],
            "baseAngleDeg": rnd(math.degrees(phi_b)), "noseAngleDeg": rnd(math.degrees(phi_n)),
            "maxContactOffset": rnd(max_off), "vertices": int(P.shape[0]),
            "rows": rows, "durations": durs,
        })
    write_fixture("mechanics_cam_three_arc.json", {"source": "tools/reference/mechanics_cam_three_arc.py", "cases": out})


if __name__ == "__main__":
    main()
