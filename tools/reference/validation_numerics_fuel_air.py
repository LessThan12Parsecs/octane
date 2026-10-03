"""Extended fuel-air-cycle oracle for the cycle-integrator validation (validation round 1).

Same ideal limiting cycle as tools/reference/cycle_fuel_air_oracle.py (the fuel-air cycle of Heywood 1988
ch. 5 / Ferguson & Kirkpatrick ch. 4 — UNVERIFIED section numbers, copied from that script), computed with Cantera on EXACTLY the project's species set
and NASA-7 thermo on its true 1-bar standard state (thermo_common.nasa_species):

  1 -> 2  isentropic compression, frozen composition (fresh charge + residual)
  2 -> 3  adiabatic constant-volume combustion to chemical equilibrium (UV) at TDC
  3 -> 4  isentropic expansion in shifting chemical equilibrium (SV)

This script widens the envelope of the original 8 cases (all phi <= 1, iso-octane/methane, dry
air, complete-combustion residual) to where the burned-gas equilibrium and the closure are
stressed: rich (phi 1.2 / 1.5 / 2.0) and very lean (0.4) mixtures, PRF blends, n-heptane,
ethanol, propane, humid air, EQUILIBRIUM residual (contains CO/H2 when rich), CR 4..18,
throttled (0.3 bar) and boosted (1.5 bar) intake, hot MON-like charge.

The reactant mole fractions are written to the fixture and used verbatim by the TypeScript test
(test/validation/numerics_fuel_air.test.ts), so the oracle does not depend on the project's
fuel/air/residual conversions.

Writes test/fixtures/validation_numerics_fuel_air.json. Run with .venv/bin/python.
"""
import json
import pathlib
import sys

import cantera as ct
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "test" / "fixtures" / "validation_numerics_fuel_air.json"
sys.path.insert(0, str(HERE))
from thermo_common import nasa_species, read_species  # noqa: E402

NAMES = read_species()
N_EQ = 12  # burned-gas equilibrium species = first 12 (core/species.ts)
_SP = nasa_species(NAMES)
full = ct.Solution(thermo="ideal-gas", species=_SP)
products = ct.Solution(thermo="ideal-gas", species=_SP[:N_EQ])

# Dry air (same as cycle_fuel_air_oracle.py; the TS test uses the X written here).
DRY_AIR = {"N2": 0.780848, "O2": 0.209476, "AR": 0.009365, "CO2": 0.000311}
FUEL_CHO = {"IC8H18": (8, 18, 0), "NC7H16": (7, 16, 0), "CH4": (1, 4, 0), "C3H8": (3, 8, 0), "C2H5OH": (2, 6, 1)}


def norm(d):
    t = sum(v for v in d.values() if v > 0)
    return {k: v / t for k, v in d.items() if v > 0}


def fresh_X(fuel_X, phi, x_h2o):
    """Fresh charge: fuel blend (mole fractions fuel_X) + humid air at equivalence ratio phi."""
    o2_st = sum(x * (c + h / 4 - o / 2) for s, x in fuel_X.items() for (c, h, o) in [FUEL_CHO[s]])
    air = {k: v * (1 - x_h2o) for k, v in DRY_AIR.items()}
    if x_h2o > 0:
        air["H2O"] = x_h2o
    n_air = o2_st / phi / air["O2"]  # mol air per mol fuel blend
    mix = {k: v * n_air for k, v in air.items()}
    for s, x in fuel_X.items():
        mix[s] = mix.get(s, 0) + x
    return norm(mix)


def product_guess(el):
    """Element-conserving major-product guess (mole fractions over the product species) for the
    element mole fractions el (C, H, O, N, Ar): lean CO2/H2O/O2; rich CO2/CO/H2O; very rich CO/H2O/H2."""
    c, h, o, n, ar = (el[e] for e in ["C", "H", "O", "N", "Ar"])
    assert o >= c, "C/O > 1 cannot be represented without condensed carbon"
    g = {"N2": n / 2, "AR": ar}
    if o >= 2 * c + h / 2:
        g.update(CO2=c, H2O=h / 2, O2=(o - 2 * c - h / 2) / 2)
    elif o >= c + h / 2:
        co2 = o - c - h / 2
        g.update(CO2=co2, CO=c - co2, H2O=h / 2)
    else:
        w = o - c
        g.update(CO=c, H2O=w, H2=(h - 2 * w) / 2)
    return norm(g)


def check_elements(phase, el):
    for e in ["C", "H", "O", "N"]:
        if el[e] > 0:
            # elemental_mole_fraction is per atom total; compare ratios to N (always present)
            r_ph = phase.elemental_mole_fraction(e) / phase.elemental_mole_fraction("N")
            r_ch = el[e] / el["N"]
            assert abs(r_ph / r_ch - 1) < 1e-9, (e, r_ph, r_ch)


def equilibrium_residual(X_fresh, T_res, p_res):
    """TP-equilibrium products (12 product species) of the fresh charge's elements at (T_res, p_res)."""
    full.TPX = 300, 1e5, X_fresh
    el = {e: full.elemental_mole_fraction(e) for e in ["C", "H", "O", "N", "Ar"]}
    products.TPX = T_res, p_res, product_guess(el)
    products.equilibrate("TP")
    check_elements(products, el)
    return {NAMES[i]: float(products.X[i]) for i in range(N_EQ) if products.X[i] > 0}


def mix_by_mass(Xa, Xb, yb):
    full.TPX = 300, 1e5, Xa
    Ma = full.mean_molecular_weight
    full.TPX = 300, 1e5, Xb
    Mb = full.mean_molecular_weight
    wa, wb = (1 - yb) / Ma, yb / Mb
    out = {}
    for k in set(Xa) | set(Xb):
        out[k] = wa * Xa.get(k, 0) + wb * Xb.get(k, 0)
    return norm(out)


def lhv_blend(fuel_X):
    """LHV of the gaseous fuel blend at 298.15 K (complete combustion, H2O vapour), J/kg of fuel."""
    T = 298.15

    def hmol(name):
        full.TPX = T, 1e5, {name: 1.0}
        return full.enthalpy_mole  # J/kmol

    dh = 0.0
    mw = 0.0
    for s, x in fuel_X.items():
        c, h, o = FUEL_CHO[s]
        dh += x * (c * hmol("CO2") + h / 2 * hmol("H2O") - hmol(s) - (c + h / 4 - o / 2) * hmol("O2"))
        full.TPX = T, 1e5, {s: 1.0}
        mw += x * full.mean_molecular_weight
    return -dh / mw


def run_case(case):
    bore, stroke, rod, cr = case["bore"], case["stroke"], case["rod"], case["cr"]
    vd = np.pi / 4 * bore**2 * stroke
    vc = vd / (cr - 1)
    v1 = vc + vd
    Xf = fresh_X(case["fuel_X"], case["phi"], case["x_h2o"])
    if case["residual"] > 0:
        Xr = equilibrium_residual(Xf, case["T_res"], 1e5)
        X0 = mix_by_mass(Xf, Xr, case["residual"])
    else:
        X0 = Xf
    lhv = lhv_blend(case["fuel_X"])  # (changes the state of `full`: call before setting the charge)
    full.TPX = case["T1"], case["p1"], X0
    m = full.density * v1
    s1, u1 = full.s, full.int_energy_mass
    fuel_mass = m * sum(full.Y[NAMES.index(s)] for s in FUEL_CHO if s in NAMES and full.Y[NAMES.index(s)] > 0)
    a = stroke / 2

    def volume(theta):
        s = a * np.cos(theta) + np.sqrt(rod**2 - (a * np.sin(theta)) ** 2)
        return vc + np.pi / 4 * bore**2 * (rod + a - s)

    comp = []
    for deg in np.linspace(-180, 0, 19):
        v = volume(np.radians(deg))
        full.SV = s1, v / m
        assert max(abs(full[k].X[0] - X0.get(k, 0)) for k in NAMES) < 1e-12  # frozen composition
        comp.append({"deg": float(deg), "p": float(full.P), "T": float(full.T)})
    full.SV = s1, vc / m
    u2 = full.int_energy_mass
    T2, p2 = full.T, full.P
    # UV equilibrium over the product species with the charge's element composition
    el = {e: full.elemental_mole_fraction(e) for e in ["C", "H", "O", "N", "Ar"]}
    products.TPX = 2500, p2, product_guess(el)
    products.UV = u2, vc / m
    products.equilibrate("UV")
    check_elements(products, el)
    T3, p3, u3, s3 = products.T, products.P, products.int_energy_mass, products.s
    X3 = {NAMES[i]: float(products.X[i]) for i in range(N_EQ)}
    exp_eq = []
    for deg in np.linspace(0, 180, 19):
        v = volume(np.radians(deg))
        products.SV = s3, v / m
        products.equilibrate("SV")
        exp_eq.append({"deg": float(deg), "p": float(products.P), "T": float(products.T)})
    products.SV = s3, v1 / m
    products.equilibrate("SV")
    u4, T4, p4 = products.int_energy_mass, products.T, products.P
    w = m * ((u1 - u2) + (u3 - u4))
    return {
        "input": case,
        "X_reactants": X0,
        "mass": m,
        "fuel_mass": fuel_mass,
        "lhv_gaseous": lhv,
        "state2": {"T": T2, "p": p2},
        "state3": {"T": T3, "p": p3, "X": X3},
        "state4": {"T": T4, "p": p4},
        "work": w,
        "efficiency": w / (fuel_mass * lhv),
        "compression_trace": comp,
        "expansion_trace_equilibrium": exp_eq,
    }


IN = 0.0254
CFR = {"bore": 3.25 * IN, "stroke": 4.5 * IN, "rod": 10.0 * IN}
ISO = {"IC8H18": 1.0}
# PRF 90 blend as MOLE fractions of the liquid-volume blend (iso-octane 0.9 vol).
# UNVERIFIED: liquid densities 691.9 / 683.8 kg/m3 (iso-octane / n-heptane, 20 °C, from memory) and molar
# masses 114.23 / 100.20 g/mol. The exact conversion is irrelevant here: the TS test takes X_reactants
# verbatim and the oracle LHV is computed for this very blend.
_v_iso, _v_hep = 0.9 * 691.9 / 114.23, 0.1 * 683.8 / 100.20
PRF90 = {"IC8H18": _v_iso / (_v_iso + _v_hep), "NC7H16": _v_hep / (_v_iso + _v_hep)}
BASE = {**CFR, "x_h2o": 0.0, "residual": 0.0, "T_res": 1000.0, "T1": 350.0, "p1": 1.0e5}
CASES = [
    {**BASE, "label": "iso-octane phi 1.2 CR 8", "fuel_X": ISO, "phi": 1.2, "cr": 8.0},
    {**BASE, "label": "iso-octane phi 1.5 CR 8", "fuel_X": ISO, "phi": 1.5, "cr": 8.0},
    {**BASE, "label": "iso-octane phi 2.0 CR 8", "fuel_X": ISO, "phi": 2.0, "cr": 8.0},
    {**BASE, "label": "iso-octane phi 0.4 CR 8", "fuel_X": ISO, "phi": 0.4, "cr": 8.0},
    {**BASE, "label": "PRF90 phi 1.1 CR 6.43 humid, 6 % eq. residual (RON-like)", "fuel_X": PRF90, "phi": 1.1,
     "cr": 6.43, "x_h2o": 0.0125, "residual": 0.06, "T1": 370.0, "p1": 1.0e5},
    {**BASE, "label": "PRF90 phi 1.05 CR 6.54 humid, 8 % eq. residual, 450 K (MON-like)", "fuel_X": PRF90, "phi": 1.05,
     "cr": 6.54, "x_h2o": 0.02, "residual": 0.08, "T1": 450.0, "p1": 0.97e5},
    {**BASE, "label": "n-heptane phi 1.0 CR 10", "fuel_X": {"NC7H16": 1.0}, "phi": 1.0, "cr": 10.0},
    {**BASE, "label": "ethanol phi 1.0 CR 12", "fuel_X": {"C2H5OH": 1.0}, "phi": 1.0, "cr": 12.0},
    {**BASE, "label": "propane phi 0.9 CR 10", "fuel_X": {"C3H8": 1.0}, "phi": 0.9, "cr": 10.0},
    {**BASE, "label": "iso-octane phi 1.0 CR 18 boosted 1.5 bar 400 K", "fuel_X": ISO, "phi": 1.0, "cr": 18.0,
     "T1": 400.0, "p1": 1.5e5},
    {**BASE, "label": "iso-octane phi 1.0 CR 4 throttled 0.3 bar", "fuel_X": ISO, "phi": 1.0, "cr": 4.0,
     "T1": 300.0, "p1": 0.3e5},
    {**BASE, "label": "methane phi 1.3 CR 14, 15 % eq. residual 420 K", "fuel_X": {"CH4": 1.0}, "phi": 1.3,
     "cr": 14.0, "residual": 0.15, "T1": 420.0},
]

if __name__ == "__main__":
    results = [run_case(c) for c in CASES]
    OUT.write_text(json.dumps({
        "description": "Extended fuel-air-cycle oracle (frozen isentropic compression, UV equilibrium at TDC, "
                       "SV equilibrium expansion), Cantera + nasa_gas.yaml at 1 bar standard state; SI; crank deg, "
                       "0 = firing TDC. Generated by tools/reference/validation_numerics_fuel_air.py.",
        "cantera_version": ct.__version__,
        "cases": results,
    }, indent=1))
    for r in results:
        print(f"{r['input']['label']:70s} T3={r['state3']['T']:.0f} K p3={r['state3']['p']/1e5:6.1f} bar "
              f"eta={r['efficiency']:.4f} W={r['work']:.2f} J")
