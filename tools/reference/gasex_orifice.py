"""Oracle for src/physics/gas-exchange/orifice.ts.

1. Constant-gamma isentropic restriction flow, evaluated independently of the TS code in
   50-digit decimal arithmetic:
      Psi(r) = sqrt(2g/(g-1) (r^(2/g) - r^((g+1)/g)))   r >= r*
      Psi*   = sqrt(g) (2/(g+1))^((g+1)/(2(g-1)))        r <  r*
   plus a numerical maximisation of Psi(r) (scipy, bounded Brent) to confirm that the
   maximum sits at r* = (2/(g+1))^(g/(g-1)) with value Psi*.
2. Variable-cp (real ideal-gas mixture, frozen composition) isentropic nozzle with Cantera:
   s(T_t, p_t) = s0, v_t = sqrt(2 (h0 - h_t)), mass flux G = rho_t v_t, throat pressure
   p_t = max(p_d, p*) where p* maximises G (the true sonic point). Gas: dry air at 300 K and
   frozen equilibrium burned gas (iso-octane/air, phi = 1) at 1600 K — quantifies the error
   of the constant-gamma (gamma at T0) formula used in the TS hot path.
   Species thermo = our species set from nasa_gas.yaml (same data as the TS thermo).

Writes test/fixtures/gasex_orifice.json.
"""
from __future__ import annotations

import os
import sys
from decimal import Decimal, getcontext

import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import build_phase, rnd, write_fixture  # noqa: E402

getcontext().prec = 50


def psi_dec(r: Decimal, g: Decimal) -> Decimal:
    one = Decimal(1)
    two = Decimal(2)
    rc = (two / (g + one)) ** (g / (g - one))
    if r >= one:
        return Decimal(0)
    if r <= rc:
        return g.sqrt() * (two / (g + one)) ** ((g + one) / (two * (g - one)))
    a = (two / g) * r.ln()
    b = ((g + one) / g) * r.ln()
    return ((two * g / (g - one)) * (a.exp() - b.exp())).sqrt()


def psi_float(r: float, g: float) -> float:
    return float(psi_dec(Decimal(r), Decimal(g)))  # exact binary values, as the TS code sees them


def analytic_cases():
    rows = []
    for g in (1.25, 1.3, 1.35, 1.4, 1.6667):
        gd = Decimal(g)
        rc = float((Decimal(2) / (gd + 1)) ** (gd / (gd - 1)))
        psic = float(gd.sqrt() * (Decimal(2) / (gd + 1)) ** ((gd + 1) / (2 * (gd - 1))))
        # numerical maximisation of the subsonic branch formula (continued below r*)
        def neg(r):
            return -np.sqrt(2 * g / (g - 1) * (r ** (2 / g) - r ** ((g + 1) / g)))
        opt = minimize_scalar(neg, bounds=(0.05, 0.99), method="bounded", options={"xatol": 1e-12})
        pts = []
        for r in (0.0, 0.1, 0.3, 0.5, rc * 0.999, rc, rc * 1.001, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99,
                  0.999, 0.9999, 1 - 1e-6, 1 - 1e-9):
            r = min(r, 1.0)
            pts.append({"r": rnd(r), "psi": rnd(psi_float(r, g))})
        rows.append({
            "gamma": g, "rCrit": rnd(rc), "psiCrit": rnd(psic),
            "rArgmaxNumeric": rnd(float(opt.x)), "psiMaxNumeric": rnd(float(-opt.fun)),
            "points": pts,
        })
    return rows


def cantera_nozzle(gas, T0, p0, X, ratios):
    gas.TPX = T0, p0, X
    h0 = gas.h
    s0 = gas.s
    Y = gas.Y.copy()

    def flux(pt):
        gas.SPY = s0, pt, Y
        v = np.sqrt(max(2 * (h0 - gas.h), 0.0))
        return gas.density * v, v, gas.T

    # true sonic point: maximise G over the throat pressure
    opt = minimize_scalar(lambda x: -flux(x * p0)[0], bounds=(0.3, 0.8), method="bounded",
                          options={"xatol": 1e-12})
    r_star = float(opt.x)
    out = []
    for r in ratios:
        rt = max(r, r_star)
        G, v, Tt = flux(rt * p0)
        out.append({"r": rnd(r), "massFlux": rnd(G), "jetVelocity": rnd(v), "throatT": rnd(Tt),
                    "choked": bool(r < r_star)})
    gas.TPX = T0, p0, X
    return {
        "T0": T0, "p0": p0,
        "X": {k: rnd(float(x)) for k, x in zip(gas.species_names, gas.X) if x > 0},
        "R": rnd(8.31446261815324 / gas.mean_molecular_weight * 1000.0),
        "gammaT0": rnd(gas.cp_mass / gas.cv_mass),
        "rStarVariableCp": rnd(r_star),
        "rows": out,
    }


def equilibrium_burned(gas, fuel, phi, air, T, p):
    """Frozen equilibrium products (mole fractions > 1e-14) of fuel/air at (phi, T, p)."""
    gas.set_equivalence_ratio(phi, fuel, air)
    gas.TP = T, p
    gas.equilibrate("TP")
    return {k: float(x) for k, x in zip(gas.species_names, gas.X) if x > 1e-14}


def fresh_charge(gas, fuel, phi, air, residual=None, residual_fraction=0.0):
    """Unburned fuel-vapour/air mixture (optionally diluted by a mole fraction of `residual`)."""
    gas.set_equivalence_ratio(phi, fuel, air)
    X = dict(zip(gas.species_names, gas.X))
    if residual:
        X = {k: (1 - residual_fraction) * X.get(k, 0.0) + residual_fraction * residual.get(k, 0.0)
             for k in set(X) | set(residual)}
    return {k: float(v) for k, v in X.items() if v > 0}


def main():
    gas = build_phase()
    ratios = [0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99, 0.999]
    air = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}
    # Burned gas: equilibrium products of iso-octane/air at phi = 1, 1600 K, 4 bar (then frozen).
    burned = equilibrium_burned(gas, "IC8H18:1", 1.0, air, 1600.0, 4e5)
    # Harder conditions (review extension): hot early-EVO blowdown of dissociated lean products,
    # rich (CO/H2-laden) products, the real intake fluids (fuel vapour, residual dilution, the
    # 149 C MON mixture temperature, ethanol), and near-unity pressure ratios.
    burned_hot_lean = equilibrium_burned(gas, "IC8H18:1", 0.7, air, 2200.0, 8e5)
    burned_rich = equilibrium_burned(gas, "IC8H18:1", 1.3, air, 1200.0, 3e5)
    fresh_ron = fresh_charge(gas, "IC8H18:0.9, NC7H16:0.1", 1.1, air, burned, 0.08)
    fresh_ethanol = fresh_charge(gas, "C2H5OH:1", 1.2, air)
    ratios_fine = ratios + [0.9999]
    variable_cp = [
        {"label": "air-300K", **cantera_nozzle(gas, 300.0, 101325.0, air, ratios)},
        {"label": "burned-1600K", **cantera_nozzle(gas, 1600.0, 4e5, burned, ratios)},
        {"label": "burned-lean-2200K-8bar", **cantera_nozzle(gas, 2200.0, 8e5, burned_hot_lean, ratios_fine)},
        {"label": "burned-rich-1200K-3bar", **cantera_nozzle(gas, 1200.0, 3e5, burned_rich, ratios_fine)},
        {"label": "fresh-prf90-phi1.1-8pct-residual-302K", **cantera_nozzle(gas, 302.3, 0.95e5, fresh_ron, ratios_fine)},
        {"label": "fresh-prf90-MON-422K", **cantera_nozzle(gas, 422.15, 1.0e5, fresh_ron, ratios_fine)},
        {"label": "fresh-ethanol-phi1.2-330K", **cantera_nozzle(gas, 330.0, 0.9e5, fresh_ethanol, ratios_fine)},
    ]
    write_fixture("gasex_orifice.json", {
        "source": "tools/reference/gasex_orifice.py",
        "analytic": analytic_cases(),
        "variableCp": variable_cp,
    })


if __name__ == "__main__":
    main()
