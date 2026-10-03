"""Gas states at the spark for the ignition-kernel tests (Cantera oracle).

CFR-like conditions (Waukesha CFR F-1 geometry: bore 3.25 in, stroke 4.5 in, rod 254 mm;
Pal et al. SAE 2018-01-0187 Table 2), RON-type operation: IVC at -146 deg (ASTM 34 deg ABDC),
spark at -13 deg (13 deg BTDC), fresh iso-octane/dry-air charge at T_IVC = 340 K,
p_IVC = 1.0 bar, no residual, isentropic (frozen, adiabatic) compression — an idealised
"CFR-like" state, not a measured one.

For each case (compression ratio, equivalence ratio [, EGR, IVC state]) we compute with
Cantera and the project's NASA-7 data (thermo_common.build_phase, identical to
src/physics/thermo):
  - unburned state at the spark: p, T_u, rho_u, X (project species order);
  - adiabatic constant-pressure equilibrium burned state from (h_u, p) over the 12
    product species: T_ad, rho_b, expansion ratio sigma = rho_u / rho_b, M_b;
  - dh_heat = h_u(T_ad) - h_u(T_u) (frozen reactants; the Herweg-Maly plasma-velocity
    denominator);
  - unburned kinematic viscosity nu_u (Cantera mixture-averaged transport with the
    project's Lennard-Jones data).

Writes test/fixtures/ignition_kernel_gas.json.
"""
from __future__ import annotations

import math
import os
import sys

import cantera as ct

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import build_phase, read_species, rnd, write_fixture  # noqa: E402

IN = 0.0254
BORE, STROKE, ROD = 3.25 * IN, 4.5 * IN, 0.254
IVC_DEG, SPARK_DEG = -146.0, -13.0
T_IVC, P_IVC = 340.0, 1.0e5
# Dry air as used by src/physics/thermo/fuels.ts DRY_AIR_REFERENCE (renormalised).
AIR = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}
# (CR, phi, EGR mole fraction, T_IVC K, p_IVC Pa). The first four are the nominal CFR-like
# states used by the tests; the rest extend the oracle to harder conditions (reviewer):
# rich (water-gas-shift products), 20 % stoichiometric-product dilution, boosted/hot high-CR,
# very lean low-CR.
CASES = [
    (7.0, 1.0, 0.0, T_IVC, P_IVC), (7.0, 0.7, 0.0, T_IVC, P_IVC), (7.0, 0.5, 0.0, T_IVC, P_IVC),
    (10.0, 1.0, 0.0, T_IVC, P_IVC),
    (7.0, 1.3, 0.0, T_IVC, P_IVC), (7.0, 1.0, 0.2, T_IVC, P_IVC), (12.0, 0.9, 0.0, 380.0, 1.5e5),
    (4.5, 0.6, 0.0, 320.0, 0.8e5),
]


def volume(theta_deg: float, cr: float) -> float:
    a = STROKE / 2
    area = math.pi * BORE ** 2 / 4
    vd = area * STROKE
    vc = vd / (cr - 1)
    th = math.radians(theta_deg)
    x = a + ROD - a * math.cos(th) - math.sqrt(ROD ** 2 - (a * math.sin(th)) ** 2)
    return vc + area * x


def main():
    species = read_species()
    gas = build_phase(transport=True, trange=(250.0, 3500.0))
    products = ct.Solution(thermo="ideal-gas", species=[gas.species(k) for k in species[:12]])
    s_air = sum(AIR.values())
    out = []
    for cr, phi, egr, t_ivc, p_ivc in CASES:
        # fresh charge: phi * stoich; iso-octane needs 12.5 O2 per mol
        x = {k: v / s_air for k, v in AIR.items()}
        n_fuel = phi * x["O2"] / 12.5
        comp = dict(x)
        comp["IC8H18"] = n_fuel
        if egr > 0:
            # dilute with complete-combustion products of a stoichiometric charge (N2, CO2,
            # H2O, Ar), `egr` = product mole fraction of the unburned mixture
            n0 = x["O2"] / 12.5
            prod = {"N2": x["N2"], "AR": x["AR"], "CO2": x["CO2"] + 8 * n0, "H2O": 9 * n0}
            sp = sum(prod.values())
            sf = sum(comp.values())
            comp = {k: (1 - egr) * v / sf for k, v in comp.items()}
            for k, v in prod.items():
                comp[k] = comp.get(k, 0.0) + egr * v / sp
        gas.TPX = t_ivc, p_ivc, comp
        s0, v0 = gas.s, gas.v
        rv = volume(IVC_DEG, cr) / volume(SPARK_DEG, cr)
        gas.SV = s0, v0 / rv
        Tu, p, rho_u, h_u = gas.T, gas.P, gas.density, gas.enthalpy_mass
        X = [float(v) for v in gas.X]
        nu = gas.viscosity / rho_u
        mu_u = gas.mean_molecular_weight
        # Equilibrium HP products (fuel species excluded) with the reactants' element totals.
        n_el = {el: 0.0 for el in gas.element_names}
        for k, sp in enumerate(gas.species_names):
            for el, nat in gas.species(sp).composition.items():
                n_el[el] += X[k] * nat
        # element-conserving initial product guess: complete combustion when lean/stoich;
        # for rich charges CO2/CO/H2O/H2 with the oxygen shared out (then equilibrated)
        nC, nH, nO = n_el["C"], n_el["H"], n_el["O"]
        if nO >= 2 * nC + nH / 2:
            guess = {"CO2": nC, "H2O": nH / 2, "O2": max(1e-12, (nO - 2 * nC - nH / 2) / 2)}
        else:
            a_min = max(0.0, nO - nC - nH / 2)  # CO2 needed so that H2 >= 0
            a = a_min + 0.5 * (min(nC, nO - nC) - a_min)
            guess = {"CO2": a, "CO": nC - a, "H2O": nO - nC - a, "H2": nH / 2 - (nO - nC - a)}
        guess.update({"N2": n_el["N"] / 2, "AR": n_el.get("Ar", 0.0)})
        for k, v in guess.items():
            assert v >= 0, (k, v)
        products.TPX = Tu, p, guess
        # match the reactant enthalpy (per unit mass: element totals are mass-consistent)
        products.HP = h_u, p
        products.equilibrate("HP")
        T_ad, rho_b, M_b = products.T, products.density, products.mean_molecular_weight
        gas.TPX = T_ad, p, X
        dh = gas.enthalpy_mass - h_u
        gas.TPX = Tu, p, X
        out.append(dict(
            cr=cr, phi=phi, egr=egr, Tivc=t_ivc, pivc=p_ivc, rv=rnd(rv), p=rnd(p), Tu=rnd(Tu), rhoU=rnd(rho_u), X=[rnd(v) for v in X],
            nuU=rnd(nu), Mu=rnd(mu_u), Tad=rnd(T_ad), rhoB=rnd(rho_b), Mb=rnd(M_b),
            sigma=rnd(rho_u / rho_b), dhHeat=rnd(dh),
        ))
        print(f"CR {cr} phi {phi} egr {egr}: p {p/1e5:.2f} bar Tu {Tu:.1f} K rho_u {rho_u:.3f} "
              f"T_ad {T_ad:.1f} K sigma {rho_u/rho_b:.3f} dh {dh/1e6:.3f} MJ/kg nu {nu:.3e}")
    write_fixture("ignition_kernel_gas.json", dict(
        conditions=dict(bore=BORE, stroke=STROKE, rod=ROD, ivcDeg=IVC_DEG, sparkDeg=SPARK_DEG,
                        T_ivc=T_IVC, p_ivc=P_IVC, note="per-case Tivc/pivc/egr override"),
        species=species, cases=out))


if __name__ == "__main__":
    main()
