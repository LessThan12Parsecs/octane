"""Oracle for the butterfly-throttle open area (throttle.ts).

Independent of the closed form: the open area is integrated numerically over the bore
cross-section. Thin elliptical plate that exactly fills the bore at the closed angle psi0
(plate plane measured from the plane normal to the bore axis), rotated to psi about a shaft
of diameter d = a D lying along x. Projected onto the bore cross-section:
  bore      : x^2 + y^2 <= R^2
  plate     : x^2 + (y/k)^2 <= R^2,   k = cos(psi)/cos(psi0)
  shaft     : |y| <= a R
Open area = integral over y of the chord length of the bore minus the plate/shaft-covered
part (scipy.integrate.quad with the kinks as breakpoints).
Writes test/fixtures/gasex_throttle.json.
"""
from __future__ import annotations

import math
import os
import sys

from scipy.integrate import quad

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402


def open_area(D, a, psi0, psi):
    R = D / 2
    k = math.cos(psi) / math.cos(psi0)

    def chord_open(y):
        c = math.sqrt(max(R * R - y * y, 0.0))
        if abs(y) <= a * R:
            return 0.0
        e = math.sqrt(max(R * R - (y / k) ** 2, 0.0)) if k > 0 else 0.0
        return 2 * max(c - e, 0.0)

    pts = sorted({a * R, min(k * R, R)})
    val, err = quad(chord_open, 0, R, points=pts, limit=200, epsabs=1e-16, epsrel=1e-13)
    return 2 * val


def main():
    cases = []
    for D, a, psi0d in [(0.05, 0.12, 8.0), (0.0143, 0.2, 5.0), (0.06, 0.0, 10.0)]:
        psi0 = math.radians(psi0d)
        rows = []
        for psid in [psi0d, psi0d + 0.01, psi0d + 1, 15, 20, 30, 45, 60, 75, 80, 85, 88, 89.5, 90]:
            psi = math.radians(psid)
            rows.append({"plateAngleDeg": psid, "area": rnd(open_area(D, a, psi0, psi))})
        cases.append({"diameter": D, "shaftRatio": a, "closedAngleDeg": psi0d, "rows": rows})
    write_fixture("gasex_throttle.json", {"source": "tools/reference/gasex_throttle.py", "cases": cases})


if __name__ == "__main__":
    main()
