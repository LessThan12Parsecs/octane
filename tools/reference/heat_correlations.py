"""Oracle for src/physics/heat-transfer (independent re-implementation + Cantera transport).

Correlations, coded here in their PUBLISHED units (the TS code works in SI internally):
 * Woschni (1967), SAE 670931, as quoted by the review "Revising engine heat transfer",
   Annals of the Faculty of Engineering Hunedoara VI(3) 2008, eq. (4) (fetched):
       h [kW/m^2K] = 0.820 D^-0.2 p^0.8 w^0.8 T^-0.53,  D [m], p [MPa], T [K], w [m/s]
       w = C1 c_m + C2 (V_s T1/(p1 V1)) (p - p0)
       C1 = 6.18 (+0.417 c_u/c_m) gas exchange, 2.28 (+0.308 c_u/c_m) otherwise;
       C2 = 0 (gas exchange, compression), 3.24e-3 m/(s K) (combustion, expansion)
 * Heywood (1988) form (quoted e.g. by the Chalmers thesis research.chalmers.se/publication/508534,
   eq. 12, fetched): h [W/m^2K] = 3.26 B^-0.2 p^0.8 T^-0.55 w^0.8, p [kPa]
 * Hohenberg (1979), SAE 790825: h [W/m^2K] = 130 V^-0.06 p^0.8 T^-0.4 (c_m + 1.4)^0.8,
   V [m^3], p [bar] (xarin.com catoolRT documentation quotes 129.8)
 * Annand (1963): q/A = a (k/B) Re^b (T - Tw) + c (T^4 - Tw^4), Re = rho c_m B / mu, b = 0.7,
   k and mu of the gas at T from Cantera mixture-averaged transport (same LJ data as TS).
Writes test/fixtures/heat_correlations.json.
"""
from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import build_phase, rnd, write_fixture  # noqa: E402

SIGMA = 5.670374419e-8


def woschni_w(phase, cm, Vs=None, T1=None, p1=None, V1=None, p=None, p0=None, swirl=0.0):
    if phase == "gas-exchange":
        return (6.18 + 0.417 * swirl) * cm
    w = (2.28 + 0.308 * swirl) * cm
    if phase == "combustion":
        w += 3.24e-3 * Vs * T1 / (p1 * V1) * (p - p0)
    return w


def main():
    rows = []
    B = 0.08255
    Vs = 6.117e-4
    for phase, p, T, p0, cm, swirl in [
        ("gas-exchange", 0.98e5, 400.0, 0.98e5, 2.286, 0.0),
        ("gas-exchange", 1.2e5, 900.0, 1.2e5, 3.429, 1.5),
        ("compression", 8e5, 600.0, 8e5, 2.286, 0.0),
        ("compression", 12e5, 700.0, 12e5, 2.286, 2.0),
        ("combustion", 35e5, 2400.0, 14e5, 2.286, 0.0),
        ("combustion", 20e5, 1800.0, 6e5, 3.429, 0.5),
    ]:
        p1, T1, V1 = 0.95e5, 340.0, 6.6e-4
        w = woschni_w(phase, cm, Vs, T1, p1, V1, p, p0, swirl)
        h_w = 0.820e3 * B ** -0.2 * (p / 1e6) ** 0.8 * w ** 0.8 * T ** -0.53
        h_hey = 3.26 * B ** -0.2 * (p / 1e3) ** 0.8 * T ** -0.55 * w ** 0.8
        V = 2.0e-4
        h_hoh = 130.0 * V ** -0.06 * (p / 1e5) ** 0.8 * T ** -0.4 * (cm + 1.4) ** 0.8
        rows.append({"phase": phase, "bore": B, "p": p, "T": T, "pMotored": p0, "meanPistonSpeed": cm,
                     "swirlRatio": swirl, "Vs": Vs, "p1": p1, "T1": T1, "V1": V1, "V": V,
                     "w": rnd(w), "hWoschni1967": rnd(h_w), "hHeywood": rnd(h_hey), "hHohenberg": rnd(h_hoh)})

    gas = build_phase(transport=True, trange=(250.0, 3500.0))
    air = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}
    gas.set_equivalence_ratio(1.0, "IC8H18:1", air)
    gas.TP = 2400.0, 40e5
    gas.equilibrate("TP")
    burned = {k: float(x) for k, x in zip(gas.species_names, gas.X) if x > 1e-12}
    annand = []
    for label, X, T, p, Tw, cm, a, c in [
        ("air-compression", air, 700.0, 10e5, 450.0, 2.286, 0.49, 0.0),
        ("burned-expansion", burned, 2400.0, 40e5, 480.0, 2.286, 0.49, 0.075 * SIGMA),
        ("burned-hot-fast", burned, 2000.0, 25e5, 500.0, 3.429, 0.7, 0.075 * SIGMA),
    ]:
        gas.TPX = T, p, X
        k = gas.thermal_conductivity
        mu = gas.viscosity
        rho = gas.density
        Re = rho * cm * B / mu
        hconv = a * k / B * Re ** 0.7
        q = hconv * (T - Tw) + c * (T ** 4 - Tw ** 4)
        annand.append({"label": label, "X": {s: rnd(v) for s, v in zip(gas.species_names, gas.X) if v > 0},
                       "T": T, "p": p, "Tw": Tw, "meanPistonSpeed": cm, "bore": B, "a": a, "c": c,
                       "k": rnd(k), "mu": rnd(mu), "rho": rnd(rho), "Re": rnd(Re), "h": rnd(hconv), "flux": rnd(q)})
    write_fixture("heat_correlations.json", {"source": "tools/reference/heat_correlations.py",
                                             "woschniHohenberg": rows, "annand": annand})


if __name__ == "__main__":
    main()
