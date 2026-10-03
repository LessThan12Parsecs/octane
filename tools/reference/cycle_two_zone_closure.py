"""Two-zone closure oracle for src/physics/cycle/closure.ts.

For a set of two-zone states (frozen unburned charge at (T_u, p); burned gas of the SAME element
composition in chemical equilibrium over the 12 product species at (T_b, p)) Cantera computes the
conserved quantities the cycle integrator carries — U_tot = m_u u_u + m_b u_b, S_u = m_u s_u,
V = m_u v_u + m_b v_b — and the TypeScript closure must recover (p, T_u, T_b) and the burned
composition from (U_tot, S_u, m_u, m_b, V).

Thermo: the project's species and NASA-7 data on the true 1-bar standard state
(thermo_common.nasa_species), so agreement is expected to ~1e-9 (Newton tolerance).

Writes test/fixtures/cycle_two_zone_closure.json.
"""
import pathlib
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from thermo_common import nasa_species, read_species, rnd, write_fixture  # noqa: E402

NAMES = read_species()
N_EQ = 12
SP = nasa_species(NAMES)
full = ct.Solution(thermo="ideal-gas", species=SP)
products = ct.Solution(thermo="ideal-gas", species=SP[:N_EQ])

# Dry air (same reference as src/physics/thermo/fuels.ts DRY_AIR_REFERENCE, renormalised)
AIR = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}


def fresh_charge(fuel: dict, phi: float, residual_mass: float) -> dict:
    """Air + fuel vapour at phi, plus residual = complete-combustion products (mass fraction)."""
    c = sum(x * {"IC8H18": 8, "NC7H16": 7}[k] for k, x in fuel.items())
    h = sum(x * {"IC8H18": 18, "NC7H16": 16}[k] for k, x in fuel.items())
    nu = c + h / 4
    air_sum = sum(AIR.values())
    x_o2 = AIR["O2"] / air_sum
    n_fuel = phi * x_o2 / nu  # per mol air
    mix = {k: v / air_sum for k, v in AIR.items()}
    for k, x in fuel.items():
        mix[k] = mix.get(k, 0.0) + n_fuel * x
    if residual_mass > 0:
        full.TPX = 300, ct.one_atm, mix
        m_fresh = full.mean_molecular_weight
        # complete-combustion products of the same charge (lean or stoichiometric only)
        prod = dict(mix)
        for k in fuel:
            prod.pop(k)
        prod["CO2"] = prod.get("CO2", 0) + c * n_fuel
        prod["H2O"] = prod.get("H2O", 0) + h / 2 * n_fuel
        prod["O2"] = prod["O2"] - nu * n_fuel
        assert prod["O2"] > -1e-12
        full.TPX = 300, ct.one_atm, prod
        m_prod = full.mean_molecular_weight
        tot_f = sum(mix.values())
        tot_p = sum(prod.values())
        wf = (1 - residual_mass) / m_fresh
        wr = residual_mass / m_prod
        out = {}
        for k in set(mix) | set(prod):
            out[k] = wf * mix.get(k, 0) / tot_f + wr * prod.get(k, 0) / tot_p
        mix = out
    s = sum(mix.values())
    return {k: v / s for k, v in mix.items() if v > 0}


def products_guess(full_gas) -> np.ndarray:
    """Element-consistent product guess (complete combustion, CO for oxygen deficit)."""
    el = {e: 0.0 for e in ["C", "H", "O", "N", "Ar"]}
    for k in range(full_gas.n_species):
        for e in el:
            el[e] += full_gas.n_atoms(k, e) * full_gas.X[k]
    x = np.zeros(N_EQ)
    names = NAMES[:N_EQ]
    x[names.index("CO2")] = el["C"]
    x[names.index("H2O")] = el["H"] / 2
    x[names.index("N2")] = el["N"] / 2
    x[names.index("AR")] = el["Ar"]
    o_left = el["O"] - 2 * el["C"] - el["H"] / 2
    if o_left >= 0:
        x[names.index("O2")] = o_left / 2
    else:
        x[names.index("CO2")] = el["C"] + o_left
        x[names.index("CO")] = -o_left
    return x / x.sum()


def case(label, fuel, phi, residual, Tu, Tb, p, m, xb):
    Xu = fresh_charge(fuel, phi, residual)
    full.TPX = Tu, p, Xu
    uu, su, vu = full.int_energy_mass, full.entropy_mass, 1 / full.density
    products.TPX = Tb, p, products_guess(full)
    products.equilibrate("TP")
    ub, vb = products.int_energy_mass, 1 / products.density
    mb = xb * m
    mu = m - mb
    Xb = [float(products.X[k]) for k in range(N_EQ)] + [0.0] * (len(NAMES) - N_EQ)
    return {
        "label": label,
        "Xu": [rnd(Xu.get(n, 0.0)) for n in NAMES],
        "mu": rnd(mu), "mb": rnd(mb),
        "U": rnd(mu * uu + mb * ub),
        "S": rnd(mu * su),
        "V": rnd(mu * vu + mb * vb),
        "p": rnd(p), "Tu": rnd(Tu), "Tb": rnd(Tb),
        "Vu": rnd(mu * vu), "Vb": rnd(mb * vb),
        "Xb": [rnd(x) for x in Xb],
    }


PRF90 = {"IC8H18": 0.8962, "NC7H16": 0.1038}  # ≈ PRF 90 vapour (fuels.ts prfIsooctaneMoleFraction(90))
ISO = {"IC8H18": 1.0}
CASES = [
    case("spark-like: tiny burned zone", PRF90, 1.0, 0.06, 640.0, 2400.0, 10.5e5, 7.3e-4, 1e-6),
    case("early flame", PRF90, 1.0, 0.06, 660.0, 2450.0, 12.5e5, 7.3e-4, 0.02),
    case("mid burn, lean", ISO, 0.8, 0.05, 700.0, 2250.0, 25e5, 7.0e-4, 0.4),
    case("mid burn, rich", PRF90, 1.2, 0.0, 720.0, 2500.0, 30e5, 7.5e-4, 0.5),
    case("late burn, high p", ISO, 1.0, 0.10, 850.0, 2750.0, 70e5, 1.0e-3, 0.9),
    case("end-gas remnant", PRF90, 1.0, 0.06, 700.0, 2300.0, 18e5, 7.3e-4, 0.9995),
    case("burned only (expansion)", PRF90, 1.0, 0.06, 600.0, 1700.0, 5e5, 7.3e-4, 1.0),
]

if __name__ == "__main__":
    write_fixture("cycle_two_zone_closure.json", {
        "description": "Two-zone closure oracle: (U, S_u, m_u, m_b, V) -> (p, T_u, T_b, X_b); Cantera "
                       + ct.__version__ + ", project species on the 1-bar NASA standard state. SI units.",
        "species": NAMES,
        "cases": CASES,
    })
