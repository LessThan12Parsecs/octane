"""Flat-follower lift of the stock Model T three-arc cam lobe (oracle for render/engine-modelt/cam.ts).

Method (deliberately different from cam.ts, which uses the closed-form piecewise support function):
  * the flank-circle centre is found NUMERICALLY (scipy brentq) on the locus |C_f| = r_f - r_b
    (internal tangency with the base circle) such that |C_f - C_n| = r_f - r_n (internal tangency with
    the nose circle);
  * the lobe outline is then sampled densely as the boundary of the convex region (base arc, both flank
    arcs, nose arc), with each arc's extent decided by which circle is outermost in that direction;
  * the follower lift is the brute-force support function max_p (p . u) - r_b over the sampled outline
    (no piecewise formulas; 40000 samples per circle, error < 1e-9 m).
  * seat-to-seat durations at a given lash are found with brentq on lift(f) = lash.

Lobe radii: Ford drawing via MTFC Tulsa 'Stock 1912+ Model T Camshaft', Table 1
(tildentechnologies.com/mtfctulsa/Cams/design_stock.htm): base 0.4060 in, flank 1.2601 in,
nose 0.0313 in, rise 0.2502 in. The spec (physics/engines/model-t.ts) rounds them to 0.406 / 1.260 /
0.031 / 0.250 in; both sets are written.

Writes test/fixtures/render_modelt_cam.json.
"""
from __future__ import annotations

import json
import math
import pathlib

import numpy as np
from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "render_modelt_cam.json"
IN = 0.0254


def lobe(rb: float, rf: float, rn: float, rise: float):
    d = rb + rise - rn
    a = rf - rb

    def g(f0: float) -> float:
        cf = np.array([-a * math.cos(f0), -a * math.sin(f0)])
        cn = np.array([d, 0.0])
        return float(np.linalg.norm(cf - cn) - (rf - rn))

    f0 = brentq(g, 1e-6, math.pi - 1e-6, xtol=1e-15)
    cfp = np.array([-a * math.cos(f0), -a * math.sin(f0)])  # flank centre for the +f side
    cfm = np.array([cfp[0], -cfp[1]])
    cn = np.array([d, 0.0])
    circles = [(np.zeros(2), rb), (cn, rn), (cfp, rf), (cfm, rf)]
    # tangency chords: beyond the nose/flank tangency only the nose circle bounds the lobe, below the
    # base/flank tangency only the base circle does; elsewhere both flank discs do.
    tn = cn + rn * (cn - cfp) / np.linalg.norm(cn - cfp)
    xn = tn[0]
    xb = rb * math.cos(f0)

    def inside(p: np.ndarray) -> bool:
        tol = 1e-12
        if np.linalg.norm(p - cfp) > rf + tol or np.linalg.norm(p - cfm) > rf + tol:
            return False
        if p[0] > xn and np.linalg.norm(p - cn) > rn + tol:
            return False
        if p[0] < xb and np.linalg.norm(p) > rb + tol:
            return False
        return True

    pts = []
    for c, r in circles:
        for t in np.linspace(0, 2 * math.pi, 40000, endpoint=False):
            p = c + r * np.array([math.cos(t), math.sin(t)])
            if inside(p):
                pts.append(p)
    outline = np.array(pts)

    def lift(f: float) -> float:
        u = np.array([math.cos(f), math.sin(f)])
        return float(np.max(outline @ u)) - rb

    return f0, lift


def build(name: str, rb: float, rf: float, rn: float, rise: float) -> dict:
    f0, lift = lobe(rb, rf, rn, rise)
    crank = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 135, 140, 150, 180]
    lifts = [lift(math.radians(c / 2)) for c in crank]
    durations = []
    for lash_in in [0.010, 0.020, 0.0256, 0.030, 0.040, 0.050]:
        lash = lash_in * IN
        fz = brentq(lambda f: lift(f) - lash, 1e-4, f0 - 1e-6, xtol=1e-12)
        durations.append({"lashIn": lash_in, "durationDeg": 4 * math.degrees(fz)})
    return {
        "name": name,
        "baseRadius": rb,
        "flankRadius": rf,
        "noseRadius": rn,
        "rise": rise,
        "flankBaseAngleDeg": math.degrees(f0),
        "crankDegFromCentre": crank,
        "grossLift": lifts,
        "seatToSeat": durations,
    }



def main() -> None:
    cases = [
        build("spec", 0.406 * IN, 1.26 * IN, 0.031 * IN, 0.25 * IN),
        build("mtfc", 0.4060 * IN, 1.2601 * IN, 0.0313 * IN, 0.2502 * IN),
    ]
    OUT.write_text(json.dumps({"source": "tools/reference/render_modelt_cam.py", "cases": cases}, indent=1))
    print(f"wrote {OUT}")
    for c in cases:
        print(c["name"], "f0 =", round(c["flankBaseAngleDeg"], 4), "lift0 =", c["grossLift"][0] / IN, "in",
              [round(d["durationDeg"], 2) for d in c["seatToSeat"]])


if __name__ == "__main__":
    main()
