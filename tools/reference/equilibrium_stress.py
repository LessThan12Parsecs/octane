"""Adversarial-review oracle fixtures for src/physics/equilibrium (beyond equilibrium_oracle.py).

Writes test/fixtures/equilibrium_stress.json with
  tp:       TP equilibria at the edges of (and beyond) the contract envelope: T 200–6000 K,
            p 1 kPa – 1000 bar, phi 0.1 … the O/C limit (iso-octane phi = 3.05, O/C = 1.025),
            50 % EGR, 10 % humidity, argon-diluted charge (Ar 'power cycle' style, 70 % of the
            oxidiser), pure air up to 6000 K (N2 dissociation, NO, N);
  hp / uv:  adiabatic flames from harsh start states (end-gas-like 900 K / 100 bar, cold lean
            300 K / 0.5 bar at phi 0.3, rich phi 2.5 at 800 K / 40 bar, Ar-diluted);
  textbook: stoichiometric fuel + (O2 + 3.76 N2) at 298.15 K / 1 atm for fuels of Turns,
            "An Introduction to Combustion", Table B.1 — including n-octane, which is NOT one of
            our species: the reactant enthalpy is computed here from the NASA Glenn n-octane fit
            (nasa_gas.yaml 'C8H18,n-octane', 1-bar basis) and handed to the TS solver as (b, H).

Method: identical to equilibrium_oracle.py (Cantera 3.2 equilibrate with solver fallbacks, then
an element-potential Newton polish against the exact b; HP/UV temperatures refined by Newton on
polished equilibria). Also records how much the adiabatic temperature and the major/NO mole
fractions move when 15 extra NASA species (HO2, H2O2, NO2, N2O, HNO, NH3, HCN, CH4, O3, HCO,
NH2, NH, CN, HNCO, NCO) are allowed — a bound on the truncation error of the 12-species set.

Run: .venv/bin/python tools/reference/equilibrium_stress.py
"""
from __future__ import annotations

import math
import os
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import equilibrium_oracle as eo  # noqa: E402
from thermo_common import P_REF_DATA, write_fixture  # noqa: E402

BAR = 1e5
P_ATM = 101325.0


def textbook_air_fuel(fuel: dict, phi: float = 1.0) -> np.ndarray:
    """fuel + (nu/phi) (O2 + 3.76 N2), mole fractions over our SPECIES."""
    f = eo.norm(eo.vec(fuel))
    nu = eo.nu_o2(f)
    x = f + (nu / phi) * eo.vec({"O2": 1.0, "N2": 3.76})
    return eo.norm(x)


def ar_diluted(fuel: dict, phi: float, ar_frac: float) -> np.ndarray:
    """Fuel in an O2/Ar/N2 oxidiser with Ar mole fraction ar_frac of the oxidiser, O2 21 %."""
    ox = eo.norm(eo.vec({"O2": 0.21, "AR": ar_frac, "N2": max(0.79 - ar_frac, 0.0)}))
    return eo.fuel_air(fuel, phi, ox)


ISO = {"IC8H18": 1}
STRESS_COMPS = [
    ("iC8H18-air phi=0.1", eo.fuel_air(ISO, 0.1, eo.air())),
    ("iC8H18-air phi=3.05 (O/C=1.025)", eo.fuel_air(ISO, 3.05, eo.air())),
    ("iC8H18-air phi=1 + 50% EGR", eo.with_egr(eo.fuel_air(ISO, 1.0, eo.air()), 0.50)),
    ("CH4-humid air(10% H2O) phi=0.8", eo.fuel_air({"CH4": 1}, 0.8, eo.air(0.10))),
    ("iC8H18 in O2/Ar(70%)/N2 phi=1", ar_diluted(ISO, 1.0, 0.70)),
    ("dry air", eo.air()),
    ("C2H5OH-air phi=2.8", eo.fuel_air({"C2H5OH": 1}, 2.8, eo.air())),
]
T_GRID = [200.0, 250.0, 300.0, 1000.0, 2500.0, 4000.0, 4500.0, 5000.0, 5500.0, 6000.0]
P_GRID = [1e3, 1e5, 1e7, 1e8]


def tp_cases():
    comps, cases = [], []
    for ci, (label, x) in enumerate(STRESS_COMPS):
        b = eo.b_per_kg(x)
        comps.append({"label": label, "b": [float(v) for v in b]})
        for T in T_GRID:
            for p in P_GRID:
                n = eo.eq_tp(b, T, p)
                row = {"c": ci, "T": T, "p": p, "X": eo.xrow(n)}
                if eo.LAST_ILL:
                    row["ill"] = list(eo.LAST_ILL)
                cases.append(row)
    return comps, cases


FLAMES = [
    ("iC8H18-air phi=1, end gas 900 K / 100 bar", eo.fuel_air(ISO, 1.0, eo.air()), 900.0, 100 * BAR),
    ("CH4-air phi=0.3, 300 K / 0.5 bar", eo.fuel_air({"CH4": 1}, 0.3, eo.air()), 300.0, 0.5 * BAR),
    ("iC8H18-air phi=2.5, 800 K / 40 bar", eo.fuel_air(ISO, 2.5, eo.air()), 800.0, 40 * BAR),
    ("iC8H18 O2/Ar(70%)/N2 phi=1, 700 K / 20 bar", ar_diluted(ISO, 1.0, 0.70), 700.0, 20 * BAR),
    ("PRF90-air phi=1 + 50% EGR, 750 K / 30 bar",
     eo.with_egr(eo.fuel_air({"IC8H18": eo.PRF90_XISO, "NC7H16": 1 - eo.PRF90_XISO}, 1.0, eo.air()), 0.5),
     750.0, 30 * BAR),
]


def flame(label, x, T0, p0, mode):
    b = eo.b_per_kg(x)
    eo.GAS_ALL.TPX = T0, p0, x
    h0, u0, v0 = eo.GAS_ALL.enthalpy_mass, eo.GAS_ALL.int_energy_mass, eo.GAS_ALL.volume_mass
    gas, idx, _, _ = eo.phase_for(b)
    maj = eo.majors_from_b(b)[idx]
    gas.TPX = 1500.0, p0, maj / maj.sum()
    if mode == "HP":
        gas.HP = h0, p0
    else:
        gas.UV = u0, v0
    eo.cantera_equilibrate(gas, mode)
    n_ct = np.zeros(eo.N_EQ)
    n_ct[idx] = gas.X / (gas.mean_molecular_weight / 1000.0)
    T, n = eo.refine_T(mode, b, h0 if mode == "HP" else u0, p0 if mode == "HP" else v0, gas.T, n_ct)
    row = {"label": label, "Xr": [float(v) for v in x], "T0": T0, "p0": p0, "T": T,
           "T_cantera": gas.T, "X": eo.xrow(n)}
    if mode == "UV":
        row["p"] = n.sum() * eo.R * T / v0
    return row


# ---- textbook (Turns Table B.1) fuels, incl. n-octane via its own NASA fit -----------------
_NASA = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}


def _sp(name):
    src = _NASA[name]
    th = src.thermo
    s = ct.Species(name, src.composition)
    s.thermo = ct.NasaPoly2(th.min_temp, th.max_temp, P_REF_DATA, th.coeffs)
    return s


PROD_CT = ["N2", "O2", "Ar", "CO2", "H2O", "CO", "H2", "OH", "H", "O", "NO", "N"]
EXTRA_CT = ["HO2", "H2O2", "NO2", "N2O", "HNO", "NH3", "HCN", "CH4", "O3", "HCO", "NH2", "NH", "CN",
            "HNCO", "NCO"]
# (label, nasa_gas.yaml fuel name, C, H, O atoms)
TEXTBOOK = [
    ("CH4", "CH4", 1, 4, 0),
    ("C3H8", "C3H8", 3, 8, 0),
    ("n-C7H16", "C7H16,n-heptane", 7, 16, 0),
    ("n-C8H18 (not in SPECIES)", "C8H18,n-octane", 8, 18, 0),
    ("iso-C8H18", "C8H18,isooctane", 8, 18, 0),
]


def textbook_cases():
    out = []
    for label, fname, c, h, o in TEXTBOOK:
        nu = c + h / 4 - o / 2
        mol = {fname: 1.0, "O2": nu, "N2": 3.76 * nu}
        g = ct.Solution(thermo="ideal-gas", species=[_sp(s) for s in dict.fromkeys(PROD_CT + [fname])])
        g.TPX = 298.15, P_ATM, mol
        ntot = 1.0 + 4.76 * nu
        H0 = g.enthalpy_mole / 1000.0 * ntot  # J per (1 mol fuel + air)
        b = [c, h, o + 2 * nu, 2 * 3.76 * nu, 0.0]  # C H O N AR, mol
        # 12-species (+ the fuel as a possible product, which stays < 1e-15) Cantera value, then the
        # polished/refined value on exactly our species set
        g.equilibrate("HP")
        T_ct = g.T
        bb = np.array(b)
        gas, idx, _, _ = eo.phase_for(bb)
        n_ct = np.zeros(eo.N_EQ)
        for k, s in enumerate(eo.PROD):
            n_ct[k] = g["Ar" if s == "AR" else s].X[0] * ntot  # start estimate only (polished below)
        T, n = eo.refine_T("HP", bb, H0, P_ATM, T_ct, n_ct)
        # extended species set (truncation bound)
        gx = ct.Solution(thermo="ideal-gas",
                         species=[_sp(s) for s in dict.fromkeys(PROD_CT + EXTRA_CT + [fname])])
        gx.TPX = 298.15, P_ATM, mol
        gx.equilibrate("HP")
        out.append({"label": label, "b": b, "H": H0, "p": P_ATM, "T": T, "T_cantera": T_ct,
                    "T_extended_species": gx.T, "X": eo.xrow(n)})
    return out


def truncation_bound():
    """Engine-like states: 12-species vs 12 + 15 extra NASA species (iso-octane-air phi=1)."""
    rows = []
    for (T, p) in [(2600.0, 60 * BAR), (2900.0, 100 * BAR), (2200.0, 20 * BAR), (1800.0, 5 * BAR)]:
        res = {}
        for tag, extra in (("12", []), ("27", EXTRA_CT)):
            g = ct.Solution(thermo="ideal-gas",
                            species=[_sp(s) for s in dict.fromkeys(PROD_CT + extra + ["C8H18,isooctane"])])
            g.TPX = T, p, {"C8H18,isooctane": 1.0, "O2": 12.5, "N2": 47.0}
            g.equilibrate("TP")
            res[tag] = {k: float(g[k].X[0]) for k in ["CO2", "H2O", "CO", "H2", "O2", "OH", "NO"]}
            if extra:
                res[tag].update({k: float(g[k].X[0]) for k in ["HO2", "NO2", "N2O", "H2O2"]})
        rows.append({"T": T, "p": p, "x12": res["12"], "x27": res["27"]})
    return rows


def main() -> None:
    comps, cases = tp_cases()
    hp = [flame(l, x, T0, p0, "HP") for (l, x, T0, p0) in FLAMES]
    uv = [flame(l, x, T0, p0, "UV") for (l, x, T0, p0) in FLAMES]
    tb = textbook_cases()
    tr = truncation_bound()
    write_fixture("equilibrium_stress.json", {
        "description": "Adversarial-review equilibrium fixtures (tools/reference/equilibrium_stress.py): "
                       "TP states T 200-6000 K, p 1 kPa-1000 bar at extreme compositions (b mol/kg, "
                       "ELEMENTS order); HP/UV flames from harsh start states (Xr over SPECIES); "
                       "textbook stoichiometric fuel+(O2+3.76N2) flames at 298.15 K/1 atm given as "
                       "(b mol, H J, p Pa) incl. n-octane from its NASA fit; 12- vs 27-species "
                       "truncation comparison. Cantera 3.2 + element-potential polish.",
        "species": eo.PROD,
        "compositions": comps,
        "tp": cases,
        "hp": hp,
        "uv": uv,
        "textbook": tb,
        "truncation": tr,
    })
    print("stats:", {k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in eo.STATS.items()})
    for r in tb:
        print(f"  {r['label']:28s} T_ad = {r['T']:.2f} K (Cantera {r['T_cantera']:.2f}, "
              f"27 species {r['T_extended_species']:.2f})")


if __name__ == "__main__":
    main()
