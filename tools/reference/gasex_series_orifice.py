"""Oracle for two compressible restrictions in series (series-orifice.ts, carburettor.ts).

The same physical problem, solved independently with scipy:
  a -- [A] -- m -- [B] -- b,  each restriction an isentropic nozzle (Heywood App. C form)
  mdot = CdA p0/sqrt(R T0) Psi(x),  x = 1 - p_down/p0,
  Psi exact for x >= delta (subsonic branch / choked Psi*), and inside the low-dp band x < delta the
  regularisation of orifice.ts: Psi = G(x) sqrt(delta) P(x/delta), P(t) = (45 t - 18 t^3 + 5 t^5)/32,
  G(x)^2 = 2g/(g-1) (1-x)^(2/g) (1 - (1-x)^k)/x, k = (g-1)/g  (coded here from those formulas).
The node is a stagnation state at the upstream (T, R, gamma). The node pressure is the root of
mdot_1 - mdot_2 (scipy.optimize.brentq on [p_down, p_up]); dmdot/dpa and dmdot/dpb are central
finite differences of the whole solve. Cases cover subsonic, part throttle, a choked throttle (idle),
a choked venturi, both choked, reverse flow with a different gas, and the regularised band.
Writes test/fixtures/gasex_series_orifice.json.
"""
from __future__ import annotations

import math
import os
import sys

from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402

DELTA = 1e-3


def psi(x, g, delta=DELTA):
    if x <= 0:
        return 0.0
    k = (g - 1) / g
    xc = 1 - (2 / (g + 1)) ** (g / (g - 1))
    if x >= xc:
        return math.sqrt(g) * (2 / (g + 1)) ** ((g + 1) / (2 * (g - 1)))
    r = 1 - x
    if x >= delta:
        return math.sqrt(2 * g / (g - 1) * (r ** (2 / g) - r ** ((g + 1) / g)))
    E = -math.expm1(k * math.log1p(-x)) / x
    G = math.sqrt(2 * g / (g - 1) * r ** (2 / g) * E)
    t = x / delta
    return G * math.sqrt(delta) * (45 * t - 18 * t ** 3 + 5 * t ** 5) / 32


def nozzle(cda, p0, T0, R, g, pd):
    if pd >= p0:
        return 0.0
    return cda * p0 * psi(1 - pd / p0, g) / math.sqrt(R * T0)


def solve(cA, cB, pa, Ta, Ra, ga, pb, Tb, Rb, gb):
    if pa >= pb:
        pu, pd, T, R, g, c1, c2, sgn = pa, pb, Ta, Ra, ga, cA, cB, 1.0
    else:
        pu, pd, T, R, g, c1, c2, sgn = pb, pa, Tb, Rb, gb, cB, cA, -1.0
    if pu == pd:
        return 0.0, pu
    f = lambda p: nozzle(c1, pu, T, R, g, p) - nozzle(c2, p, T, R, g, pd)  # noqa: E731
    pm = brentq(f, pd, pu, xtol=1e-13, rtol=1e-15, maxiter=500)
    m = 0.5 * (nozzle(c1, pu, T, R, g, pm) + nozzle(c2, pm, T, R, g, pd))
    return sgn * m, pm


AIR = (288.7, 287.0, 1.4)
CASES = [
    dict(name="wot-subsonic", cA=1.571e-4, cB=2.3e-4, pa=101325.0, a=AIR, pb=95000.0, b=(300.0, 287.0, 1.4)),
    dict(name="part-throttle", cA=1.571e-4, cB=3e-5, pa=101325.0, a=AIR, pb=60000.0, b=(300.0, 287.0, 1.4)),
    dict(name="idle-choked-throttle", cA=1.571e-4, cB=4e-6, pa=101325.0, a=AIR, pb=30000.0, b=(300.0, 287.0, 1.4)),
    dict(name="venturi-choked", cA=2e-5, cB=1e-3, pa=101325.0, a=AIR, pb=40000.0, b=(300.0, 287.0, 1.4)),
    dict(name="both-choked", cA=1e-5, cB=1e-4, pa=101325.0, a=AIR, pb=5000.0, b=(300.0, 287.0, 1.4)),
    dict(name="reverse-other-gas", cA=1.571e-4, cB=2.3e-4, pa=101325.0, a=AIR, pb=104000.0, b=(330.0, 290.0, 1.36)),
    dict(name="regularised-band", cA=1.571e-4, cB=2.3e-4, pa=101325.0, a=AIR, pb=101325.0 - 50.0, b=(300.0, 287.0, 1.4)),
    dict(name="reverse-choked-a-side", cA=3e-6, cB=2e-4, pa=20000.0, a=(300.0, 287.0, 1.4), pb=101325.0, b=(350.0, 287.0, 1.38)),
]


def main():
    rows = []
    for c in CASES:
        Ta, Ra, ga = c["a"]
        Tb, Rb, gb = c["b"]
        args = lambda pa, pb: (c["cA"], c["cB"], pa, Ta, Ra, ga, pb, Tb, Rb, gb)  # noqa: E731
        m, pm = solve(*args(c["pa"], c["pb"]))
        h = min(1.0, 1e-3 * abs(c["pa"] - c["pb"]))
        dpa = (solve(*args(c["pa"] + h, c["pb"]))[0] - solve(*args(c["pa"] - h, c["pb"]))[0]) / (2 * h)
        dpb = (solve(*args(c["pa"], c["pb"] + h))[0] - solve(*args(c["pa"], c["pb"] - h))[0]) / (2 * h)
        pu = max(c["pa"], c["pb"])
        pd = min(c["pa"], c["pb"])
        g = ga if c["pa"] >= c["pb"] else gb
        rc = (2 / (g + 1)) ** (g / (g - 1))
        first_choked = pm / pu <= rc
        second_choked = pd / pm <= rc
        rows.append({"name": c["name"], "cdaA": c["cA"], "cdaB": c["cB"], "pa": c["pa"], "Ta": Ta, "Ra": Ra, "gammaA": ga,
                     "pb": c["pb"], "Tb": Tb, "Rb": Rb, "gammaB": gb, "mdot": rnd(m), "intermediatePressure": rnd(pm),
                     "dmdotdpa": rnd(dpa), "dmdotdpb": rnd(dpb),
                     "chokedFirst": bool(first_choked), "chokedSecond": bool(second_choked)})
    write_fixture("gasex_series_orifice.json", {"source": "tools/reference/gasex_series_orifice.py", "delta": DELTA, "cases": rows})


if __name__ == "__main__":
    main()
