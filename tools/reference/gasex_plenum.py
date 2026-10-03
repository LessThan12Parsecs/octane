"""Oracle for src/physics/gas-exchange/plenum.ts (Cantera 3.2, same NASA-7 species data).

1. stateRecovery: given species moles N (mol), internal energy U (J) and volume V (m^3),
   Cantera's UV solve (gas.UVY) gives T and p.
2. filling: adiabatic rigid tank (T1, p1, X1) filled from an infinite reservoir (T0, p0, X0)
   until p = p0. Energy balance U2 = U1 + dm h0 with N2 = N1 + dm Y0/M; dm found with brentq so
   that p(N2, U2, V) = p0. Path independent: the TS ODE must converge to it.
3. mixing: two rigid volumes (A, B) that equalise pressure adiabatically through an orifice.
   The final equal-pressure state is NOT unique for a quasi-steady orifice process (it depends on
   how much of each gas crossed), so only the invariants are stored: total N_k and U, and the
   single-volume "fully mixed" state (N_A+N_B, U_A+U_B, V_A+V_B) for reference.
Writes test/fixtures/gasex_plenum.json.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import build_phase, rnd, write_fixture  # noqa: E402

RU = 8.31446261815324  # J/(mol K)
AIR = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}


def moles(gas, T, p, X, V):
    gas.TPX = T, p, X
    n = p * V / (RU * T)
    N = n * gas.X
    U = gas.int_energy_mole / 1000.0 * n  # J/kmol -> J/mol
    return N, U


def state_from_NU(gas, N, U, V):
    M = gas.molecular_weights / 1000.0  # kg/mol
    m = float(np.dot(N, M))
    Y = N * M / m
    gas.UVY = U / m, V / m, Y
    return gas.T, gas.P


def vec(gas, N):
    return {k: rnd(float(v)) for k, v in zip(gas.species_names, N) if v != 0}


def main():
    gas = build_phase()
    # burned-gas-like composition: equilibrium products at 1500 K, 1 bar, phi = 1.1 iso-octane
    gas.set_equivalence_ratio(1.1, "IC8H18:1", AIR)
    fresh = dict(zip(gas.species_names, gas.X))
    gas.TP = 1500.0, 1e5
    gas.equilibrate("TP")
    burned = dict(zip(gas.species_names, gas.X))

    # Harder states (review extension): strongly dissociated frozen products at 3000 K / 50 bar and
    # 2500 K / 0.5 bar, a rich ethanol charge, a PRF charge diluted with 15 % residual, cold thin air.
    gas.TP = 3000.0, 50e5
    gas.equilibrate("TP")
    burned_hot = {k: float(v) for k, v in zip(gas.species_names, gas.X) if v > 1e-14}
    gas.set_equivalence_ratio(0.8, "IC8H18:1", AIR)
    gas.TP = 2500.0, 0.5e5
    gas.equilibrate("TP")
    burned_lowp = {k: float(v) for k, v in zip(gas.species_names, gas.X) if v > 1e-14}
    gas.set_equivalence_ratio(1.5, "C2H5OH:1", AIR)
    ethanol_rich = dict(zip(gas.species_names, gas.X))
    gas.set_equivalence_ratio(1.1, "IC8H18:0.9, NC7H16:0.1", AIR)
    prf = dict(zip(gas.species_names, gas.X))
    prf_egr = {k: 0.85 * prf.get(k, 0.0) + 0.15 * burned.get(k, 0.0) for k in set(prf) | set(burned)}

    recovery = []
    for label, T, p, X, V in [("air", 300.0, 101325.0, AIR, 1e-3),
                               ("fresh-charge", 330.0, 0.95e5, fresh, 5e-4),
                               ("burned", 1200.0, 3e5, burned, 2e-4),
                               ("burned-dissociated-3000K-50bar", 3000.0, 50e5, burned_hot, 5e-5),
                               ("burned-lean-2500K-0.5bar", 2500.0, 0.5e5, burned_lowp, 6e-4),
                               ("ethanol-rich-phi1.5", 330.0, 0.9e5, ethanol_rich, 5e-4),
                               ("prf90-15pct-residual", 400.0, 1.1e5, prf_egr, 7e-4),
                               ("air-cold-thin", 230.0, 0.2e5, AIR, 1e-2)]:
        N, U = moles(gas, T, p, X, V)
        # perturb U to make the target temperature non-trivial
        U2 = U + 50.0
        T2, p2 = state_from_NU(gas, N, U2, V)
        recovery.append({"label": label, "V": V, "N": vec(gas, N), "U": rnd(U2), "T": rnd(T2), "p": rnd(p2)})

    filling = []
    for label, (T1, p1, X1), (T0, p0, X0), V in [
        ("air-into-air", (300.0, 0.5e5, AIR), (300.0, 1.0e5, AIR), 1e-3),
        ("fresh-into-burned", (900.0, 0.9e5, burned), (320.0, 1.2e5, fresh), 5e-4),
        # exhaust backflow into the intake plenum (hot products into a fresh charge)
        ("burned-into-fresh", (320.0, 0.95e5, fresh), (1000.0, 1.1e5, burned), 1e-3),
    ]:
        N1, U1 = moles(gas, T1, p1, X1, V)
        gas.TPX = T0, p0, X0
        h0 = gas.enthalpy_mass
        Y0 = gas.Y.copy()
        M = gas.molecular_weights / 1000.0

        def resid(dm):
            N2 = N1 + dm * Y0 / M
            U2 = U1 + dm * h0
            return state_from_NU(gas, N2, U2, V)[1] - p0

        dm = brentq(resid, 0.0, 1.0, xtol=1e-18, rtol=1e-15)
        N2 = N1 + dm * Y0 / M
        U2 = U1 + dm * h0
        T2, p2 = state_from_NU(gas, N2, U2, V)
        filling.append({"label": label, "V": V,
                        "tank": {"T": T1, "p": p1, "X": {k: rnd(float(v)) for k, v in X1.items() if v > 0}},
                        "reservoir": {"T": T0, "p": p0, "X": {k: rnd(float(v)) for k, v in X0.items() if v > 0},
                                      "h0": rnd(h0)},
                        "admittedMass": rnd(dm), "T2": rnd(T2), "p2": rnd(p2)})

    NA, UA = moles(gas, 1000.0, 3e5, burned, 2e-4)
    NB, UB = moles(gas, 320.0, 1e5, fresh, 1e-3)
    Tm, pm = state_from_NU(gas, NA + NB, UA + UB, 1.2e-3)
    mixing = {"A": {"T": 1000.0, "p": 3e5, "V": 2e-4}, "B": {"T": 320.0, "p": 1e5, "V": 1e-3},
              "Ntotal": vec(gas, NA + NB), "Utotal": rnd(UA + UB), "mixedT": rnd(Tm), "mixedP": rnd(pm)}
    write_fixture("gasex_plenum.json", {
        "source": "tools/reference/gasex_plenum.py",
        "compositions": {"air": AIR,
                         "fresh": {k: rnd(float(v)) for k, v in fresh.items() if v > 0},
                         "burned": {k: rnd(float(v)) for k, v in burned.items() if v > 1e-14}},
        "stateRecovery": recovery, "filling": filling, "mixing": mixing,
    })


if __name__ == "__main__":
    main()
