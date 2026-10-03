"""Fuel-air cycle oracle (Heywood 1988 §5.5, Ferguson & Kirkpatrick ch. 4).

Ideal limiting cycle computed with Cantera using EXACTLY the project's species set
and NASA thermo data (nasa_gas.yaml coefficients on their true 1-bar standard state,
via thermo_common.nasa_species):

  1 -> 2  isentropic compression, frozen composition (fresh charge + residual)
  2 -> 3  adiabatic constant-volume combustion to chemical equilibrium (UV)
  3 -> 4  isentropic expansion in shifting chemical equilibrium (SV)
          (+ a frozen-expansion variant for comparison)

The cycle integrator, run with heat transfer off, instantaneous combustion at TDC
and equilibrium burned gas, must reproduce these pressures/temperatures/work.

Writes test/fixtures/cycle_fuel_air.json.
"""
import json
import pathlib

import cantera as ct
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "cycle_fuel_air.json"

# project species order (src/physics/core/species.ts) -> nasa_gas.yaml names
SPECIES = [
    ("N2", "N2"), ("O2", "O2"), ("AR", "Ar"), ("CO2", "CO2"), ("H2O", "H2O"),
    ("CO", "CO"), ("H2", "H2"), ("OH", "OH"), ("H", "H"), ("O", "O"),
    ("NO", "NO"), ("N", "N"), ("IC8H18", "C8H18,isooctane"),
    ("NC7H16", "C7H16,n-heptane"), ("CH4", "CH4"), ("C3H8", "C3H8"),
    ("C2H5OH", "C2H5OH"),
]
N_EQ = 12

import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from thermo_common import nasa_species  # noqa: E402  (1-bar standard state, our species names)

_SP = nasa_species([p for p, _ in SPECIES])
full = ct.Solution(thermo="ideal-gas", species=_SP)
products = ct.Solution(thermo="ideal-gas", species=_SP[:N_EQ])
proj_names = [p for p, _ in SPECIES]

# Dry air used by this oracle (stated explicitly so the TS test can use the same X).
AIR = {"N2": 0.780848, "O2": 0.209476, "AR": 0.009365, "CO2": 0.000311}
FUEL_C_H = {"IC8H18": (8, 18), "CH4": (1, 4)}


def reactants_X(fuel, phi, residual_mass_fraction):
    """Fresh charge at phi mixed with stoichiometric-equilibrium-free complete products."""
    c, h = FUEL_C_H[fuel]
    o2_per_fuel = c + h / 4.0
    x = dict(AIR)
    n_air_per_fuel = o2_per_fuel / phi / AIR["O2"]
    fresh = {k: v * n_air_per_fuel for k, v in x.items()}
    fresh[fuel] = 1.0
    tot = sum(fresh.values())
    fresh = {k: v / tot for k, v in fresh.items()}
    if residual_mass_fraction <= 0:
        return fresh
    # residual = complete-combustion products of the same fresh charge (lean/stoich only)
    prod = dict(fresh)
    nf = prod.pop(fuel)
    prod["CO2"] = prod.get("CO2", 0) + c * nf
    prod["H2O"] = prod.get("H2O", 0) + h / 2 * nf
    prod["O2"] = prod["O2"] - o2_per_fuel * nf
    assert prod["O2"] >= -1e-12, "oracle residual model supports phi <= 1 only"
    full.TPX = 300, ct.one_atm, _cantera_X(fresh)
    m_fresh = full.mean_molecular_weight
    full.TPX = 300, ct.one_atm, _cantera_X(_norm(prod))
    m_prod = full.mean_molecular_weight
    # mixing by mass: moles_i = Y_mass * X_i / M
    wf = (1 - residual_mass_fraction) / m_fresh
    wr = residual_mass_fraction / m_prod
    pn = _norm(prod)
    mix = {}
    for k in set(fresh) | set(pn):
        mix[k] = wf * fresh.get(k, 0) + wr * pn.get(k, 0)
    return _norm(mix)


def _norm(d):
    t = sum(d.values())
    return {k: v / t for k, v in d.items() if v > 0}


def _cantera_X(proj_x):
    return dict(proj_x)


def _proj_X(gas):
    return {proj_names[i]: float(gas.X[i]) for i in range(gas.n_species) if gas.X[i] > 0}


def element_moles_per_kg(gas):
    """Element moles per kg of mixture in project order C,H,O,N,AR."""
    out = []
    for el in ["C", "H", "O", "N", "Ar"]:
        out.append(sum(
            gas.n_atoms(k, el) * gas.Y[k] / gas.molecular_weights[k] * 1000.0
            for k in range(gas.n_species)))
    return out


def run_case(case):
    bore, stroke, rod = case["bore"], case["stroke"], case["rod"]
    cr = case["cr"]
    vd = np.pi / 4 * bore**2 * stroke
    vc = vd / (cr - 1)
    v1 = vc + vd
    X0 = reactants_X(case["fuel"], case["phi"], case["residual"])
    full.TPX = case["T1"], case["p1"], _cantera_X(X0)
    m = full.density * v1  # kg trapped
    s1, u1 = full.s, full.int_energy_mass
    rho1 = full.density
    fuel_idx = proj_names.index(case["fuel"])
    fuel_mass = m * full.Y[fuel_idx]
    # LHV of the gaseous fuel at 298.15 K (complete combustion, H2O vapour)
    lhv = lhv_gaseous(case["fuel"])

    # slider-crank volume for pressure traces (zero pin offset)
    a = stroke / 2

    def volume(theta):
        s = a * np.cos(theta) + np.sqrt(rod**2 - (a * np.sin(theta)) ** 2)
        return vc + np.pi / 4 * bore**2 * (rod + a - s)

    # compression -180 -> 0 (frozen, isentropic)
    comp = []
    for deg in np.linspace(-180, 0, 37):
        v = volume(np.radians(deg))
        full.SV = s1, v / m
        comp.append({"deg": float(deg), "V": float(v), "p": float(full.P), "T": float(full.T)})
    full.SV = s1, vc / m
    u2 = full.int_energy_mass
    T2, p2 = full.T, full.P

    # constant-volume combustion to equilibrium (only product species can form):
    # equilibrate the element composition over the product set at fixed u, v.
    b = element_moles_per_kg(full)
    # product phase with the same element composition (complete-combustion surrogate)
    products.X = _elements_to_products_guess(b)
    products.UV = u2, vc / m
    products.equilibrate("UV")
    _check_elements(products, b)
    T3, p3, u3, s3 = products.T, products.P, products.int_energy_mass, products.s
    X3 = {proj_names[i]: float(products.X[i]) for i in range(products.n_species)}

    # equilibrium isentropic expansion 0 -> 180
    exp_eq = []
    for deg in np.linspace(0, 180, 37):
        v = volume(np.radians(deg))
        products.SV = s3, v / m
        products.equilibrate("SV")
        exp_eq.append({"deg": float(deg), "V": float(v), "p": float(products.P), "T": float(products.T)})
    products.SV = s3, v1 / m
    products.equilibrate("SV")
    u4 = products.int_energy_mass
    T4, p4 = products.T, products.P

    # frozen expansion variant
    products.UV = u3, vc / m
    products.equilibrate("UV")
    Xfrozen = products.X.copy()
    products.SV = s3, v1 / m  # frozen: composition held
    u4f = products.int_energy_mass

    w = m * ((u1 - u2) + (u3 - u4))
    wf = m * ((u1 - u2) + (u3 - u4f))
    return {
        "input": case,
        "X_reactants": X0,
        "element_moles_per_kg": b,
        "Vd": vd, "Vc": vc, "mass": m, "fuel_mass": fuel_mass, "lhv_gaseous": lhv,
        "state1": {"T": case["T1"], "p": case["p1"], "rho": rho1, "u": u1, "s": s1},
        "state2": {"T": T2, "p": p2, "u": u2},
        "state3": {"T": T3, "p": p3, "u": u3, "s": s3, "X": X3},
        "state4_equilibrium": {"T": T4, "p": p4, "u": u4},
        "work_equilibrium_expansion": w,
        "work_frozen_expansion": wf,
        "efficiency_equilibrium": w / (fuel_mass * lhv),
        "efficiency_frozen": wf / (fuel_mass * lhv),
        "imep_equilibrium": w / vd,
        "compression_trace": comp,
        "expansion_trace_equilibrium": exp_eq,
    }


def _elements_to_products_guess(b):
    """Element-conserving major-product guess: lean CO2/H2O/O2; rich CO2/CO/H2O; very rich (o < c + h/2)
    CO/H2O/H2. (Validation round 1: the very-rich branch used to put CO = c AND H2O = h/2, which needs
    more O than available — Cantera then equilibrated the GUESS's element composition; latent, all
    committed cases are phi <= 1. Same fix as validation_numerics_fuel_air.product_guess.)"""
    c, h, o, n, ar = b
    assert o >= c, "C/O > 1 cannot be represented without condensed carbon"
    x = np.zeros(products.n_species)
    names = proj_names[:N_EQ]
    x[names.index("N2")] = n / 2
    x[names.index("AR")] = ar
    if o >= 2 * c + h / 2:
        x[names.index("CO2")] = c
        x[names.index("H2O")] = h / 2
        x[names.index("O2")] = (o - 2 * c - h / 2) / 2
    elif o >= c + h / 2:  # rich: part of the carbon as CO
        co2 = o - c - h / 2
        x[names.index("CO2")] = co2
        x[names.index("CO")] = c - co2
        x[names.index("H2O")] = h / 2
    else:  # very rich: all carbon as CO, the remaining O as H2O, the rest of H as H2
        w = o - c
        x[names.index("CO")] = c
        x[names.index("H2O")] = w
        x[names.index("H2")] = (h - 2 * w) / 2
    return x / x.sum()


def _check_elements(phase, b):
    """Assert that `phase` carries the element ratios of b (C, H, O, N, Ar per kg) after equilibrate."""
    names = ["C", "H", "O", "N"]
    ref = b[3]
    for e, v in zip(names, b[:4]):
        if v > 0:
            r = phase.elemental_mole_fraction(e) / phase.elemental_mole_fraction("N")
            assert abs(r / (v / ref) - 1) < 1e-9, (e, r, v / ref)


def lhv_gaseous(fuel):
    c, h = FUEL_C_H[fuel]
    T = 298.15
    g = ct.Solution(thermo="ideal-gas", species=_SP)
    def hmol(name):
        g.TPX = T, ct.one_atm, {name: 1.0}
        return g.enthalpy_mole  # J/kmol
    dh = c * hmol("CO2") + h / 2 * hmol("H2O") - hmol(fuel) - (c + h / 4) * hmol("O2")
    g.TPX = T, ct.one_atm, {fuel: 1.0}
    mw = g.mean_molecular_weight  # kg/kmol
    return -dh / mw  # J/kg


IN = 0.0254
CFR = {"bore": 3.25 * IN, "stroke": 4.5 * IN, "rod": 10.0 * IN}
CASES = []
for cr in [5.0, 7.0, 10.0, 14.0]:
    CASES.append({**CFR, "cr": cr, "fuel": "IC8H18", "phi": 1.0, "residual": 0.0, "T1": 350.0, "p1": 1.0e5})
for phi in [0.6, 0.8]:
    CASES.append({**CFR, "cr": 8.0, "fuel": "IC8H18", "phi": phi, "residual": 0.0, "T1": 350.0, "p1": 1.0e5})
CASES.append({**CFR, "cr": 8.0, "fuel": "IC8H18", "phi": 1.0, "residual": 0.08, "T1": 380.0, "p1": 0.95e5})
CASES.append({**CFR, "cr": 10.0, "fuel": "CH4", "phi": 1.0, "residual": 0.0, "T1": 350.0, "p1": 1.0e5})

if __name__ == "__main__":
    results = [run_case(c) for c in CASES]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps({
        "description": "Fuel-air cycle oracle (frozen isentropic compression, UV equilibrium combustion, "
                       "equilibrium isentropic expansion) with Cantera + nasa_gas.yaml; SI units, "
                       "energies per J, crank deg with 0 = firing TDC.",
        "cantera_version": ct.__version__,
        "air_X": AIR,
        "cases": results,
    }, indent=1))
    for r in results:
        i = r["input"]
        print(f"{i['fuel']:7s} CR={i['cr']:4.1f} phi={i['phi']:.1f} xr={i['residual']:.2f}: "
              f"T3={r['state3']['T']:.0f} K p3={r['state3']['p']/1e5:.1f} bar "
              f"eta_eq={r['efficiency_equilibrium']:.4f} eta_frozen={r['efficiency_frozen']:.4f} "
              f"IMEP={r['imep_equilibrium']/1e5:.2f} bar")
