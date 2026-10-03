"""Oracle fixture for src/physics/thermo/transport.ts.

Two oracles are written:

1. `cantera` — Cantera 3.2 mixture-averaged transport of a phase with exactly our species,
   NASA-7 thermo (nasa_gas.yaml, 1-bar standard state) and the same Lennard-Jones
   parameters as species-data.ts (gri30.yaml + LLNL values for the liquid fuels; see
   thermo_common.py). Cantera evaluates Chapman–Enskog pure-species properties with the
   Monchick–Mason collision-integral tables and the Warnatz/Kee conductivity model, fits
   them in ln T (degree-4 polynomials) and mixes with Wilke / Mathur–Saxena. The phase's
   validity range is narrowed to 300–3000 K (same coefficients) so that Cantera's own fit
   error stays ≤ 0.5 % (over 200–6000 K it reaches ~6 % for the fuel molecules near 250 K);
   the reported fitting errors are stored in the fixture.

2. `model` — what transport.ts is meant to compute: the same Cantera pure-species values
   EXCEPT H2O, whose viscosity and conductivity come from the IAPWS dilute-gas correlations
   (IAPWS R12-08 / R15-11; the Chapman–Enskog + Warnatz values with the GRI-Mech H2O
   parameters over-predict the steam conductivity by 20–40 %; above 2500 K the kinetic-
   theory value scaled to IAPWS at 2500 K, as in transport.ts), combined here in numpy with
   Wilke (1950) and Mathur–Saxena (1967) exactly as Cantera's MixTransport does. For
   H2O-free mixtures the script asserts model == cantera (≤ 1e-9), which validates the
   numpy mixing rules against Cantera itself.

Writes test/fixtures/thermo_transport.json:
  pure:     Cantera μ_k, λ_k for every species on a T grid (+ the H2O values of `model`)
  mixtures: μ, λ, thermal diffusivity (cantera and model) for air, fresh charges, EGR and
            burned gas, 300–3000 K
Run: .venv/bin/python tools/reference/thermo_transport.py
"""
from __future__ import annotations

import os
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import (  # noqa: E402
    IAPWS_T_MAX, IAPWS_T_MIN, build_phase, read_species, rnd, water_conductivity_dilute,
    water_viscosity_dilute, write_fixture)

SPECIES = read_species()
IDX = {s: i for i, s in enumerate(SPECIES)}
NS = len(SPECIES)
DRY_AIR = {"N2": 0.780848, "O2": 0.209390, "AR": 0.009332, "CO2": 0.00040}


def vec(d: dict) -> np.ndarray:
    x = np.zeros(NS)
    for k, v in d.items():
        x[IDX[k]] += v
    return x / x.sum()


def main() -> None:
    gas = build_phase(transport=True, trange=(300.0, 3000.0))
    mw = gas.molecular_weights
    air = vec(DRY_AIR)

    def charge(fuel: dict, nu: float, phi: float, extra: dict | None = None) -> np.ndarray:
        x = air * 1.0
        nf = phi * air[IDX["O2"]] / nu
        for k, v in fuel.items():
            x[IDX[k]] += nf * v
        if extra:
            for k, v in extra.items():
                x[IDX[k]] += v
        return x / x.sum()

    stoich_iso = charge({"IC8H18": 1.0}, 12.5, 1.0)
    prf90_humid = charge({"IC8H18": 0.8887, "NC7H16": 0.1113}, 0.8887 * 12.5 + 0.1113 * 11, 0.9,
                         {"H2O": 0.015})
    ethanol = charge({"C2H5OH": 1.0}, 3.0, 1.0)
    eq = ct.Solution(thermo="ideal-gas", species=[gas.species(s) for s in SPECIES[:12]])

    def burned(phi: float, T: float, p: float) -> np.ndarray:
        """TP-equilibrium products of iso-octane/air at φ over the 12 burned-gas species."""
        r = charge({"IC8H18": 1.0}, 12.5, phi)
        nf = r[IDX["IC8H18"]]
        seed = {s: float(r[IDX[s]]) for s in ["N2", "AR", "CO2", "O2"]}
        seed["CO"] = 8 * nf
        seed["H2"] = 9 * nf
        eq.TPX = T, p, seed
        eq.equilibrate("TP")
        out = np.zeros(NS)
        for s in SPECIES[:12]:
            out[IDX[s]] = eq.X[eq.species_index(s)]
        return out

    burned_stoich = burned(1.0, 2500.0, 40e5)
    burned_rich = burned(1.4, 2400.0, 60e5)
    # 25 % (by mole) of cooled complete products mixed into a stoichiometric charge
    egr_charge = 0.75 * stoich_iso + 0.25 * burned(1.0, 1500.0, 1e5)
    mixtures = [
        ("dry air", air),
        ("stoichiometric iso-octane/air", stoich_iso),
        ("PRF90 phi=0.9 humid", prf90_humid),
        ("ethanol/air phi=1", ethanol),
        ("burned gas (eq. 2500K 40bar)", burned_stoich),
        ("rich burned gas phi=1.4 (eq. 2400K 60bar)", burned_rich),
        ("stoichiometric charge + 25% products", egr_charge / egr_charge.sum()),
        ("pure H2O", vec({"H2O": 1.0})),
    ]
    tgrid = [300.0, 400.0, 600.0, 800.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
    p = 101325.0

    def pure_props(T: float) -> tuple[np.ndarray, np.ndarray]:
        mu = np.zeros(NS)
        lam = np.zeros(NS)
        for k, s in enumerate(SPECIES):
            gas.TPX = T, p, {s: 1.0}
            mu[k] = gas.viscosity
            lam[k] = gas.thermal_conductivity
        return mu, lam

    def wilke(x: np.ndarray, mu: np.ndarray) -> float:
        # Cantera MixTransport clips mole fractions at Tiny = 1e-20 (negligible); use x > 0.
        idx = np.nonzero(x > 0)[0]
        out = 0.0
        for k in idx:
            phi = (1 + np.sqrt(mu[k] / mu[idx]) * (mw[idx] / mw[k]) ** 0.25) ** 2 \
                / np.sqrt(8 * (1 + mw[k] / mw[idx]))
            out += x[k] * mu[k] / np.dot(x[idx], phi)
        return float(out)

    def mathur_saxena(x: np.ndarray, lam: np.ndarray) -> float:
        idx = np.nonzero(x > 0)[0]
        return float(0.5 * (np.dot(x[idx], lam[idx]) + 1 / np.dot(x[idx], 1 / lam[idx])))

    pure_cache = {T: pure_props(T) for T in set(tgrid) | {IAPWS_T_MAX}}
    assert min(tgrid) >= IAPWS_T_MIN
    kw = IDX["H2O"]

    def water(T: float) -> tuple[float, float]:
        """H2O μ, λ as in transport.ts: IAPWS dilute gas up to IAPWS_T_MAX, above it the
        kinetic-theory (here: Cantera) value scaled by its IAPWS ratio at IAPWS_T_MAX."""
        if T <= IAPWS_T_MAX:
            return water_viscosity_dilute(T), water_conductivity_dilute(T)
        mu_t, lam_t = pure_cache[T] if T in pure_cache else pure_props(T)
        mu_l, lam_l = pure_cache[IAPWS_T_MAX]
        return (water_viscosity_dilute(IAPWS_T_MAX) * mu_t[kw] / mu_l[kw],
                water_conductivity_dilute(IAPWS_T_MAX) * lam_t[kw] / lam_l[kw])

    cases = []
    worst_self = 0.0
    for name, x in mixtures:
        rows = []
        for T in tgrid:
            gas.TPX = T, p, x
            rho_cp = gas.density * gas.cp_mass
            mu_ct, lam_ct = gas.viscosity, gas.thermal_conductivity
            mu_k, lam_k = (a.copy() for a in pure_cache[T])
            mu_np, lam_np = wilke(x, mu_k), mathur_saxena(x, lam_k)
            mu_k[kw], lam_k[kw] = water(T)
            mu_m, lam_m = wilke(x, mu_k), mathur_saxena(x, lam_k)
            if x[IDX["H2O"]] == 0:
                worst_self = max(worst_self, abs(mu_np / mu_ct - 1), abs(lam_np / lam_ct - 1))
            rows.append(dict(T=T,
                             viscosity=rnd(mu_ct), conductivity=rnd(lam_ct),
                             thermalDiffusivity=rnd(lam_ct / rho_cp),
                             modelViscosity=rnd(mu_m), modelConductivity=rnd(lam_m),
                             modelThermalDiffusivity=rnd(lam_m / rho_cp)))
        cases.append(dict(name=name, X=[rnd(v) for v in x], p=p, states=rows))
    # numpy Wilke / Mathur–Saxena must reproduce Cantera's MixTransport for H2O-free mixtures
    assert worst_self < 1e-9, worst_self
    print(f"numpy mixing rules vs Cantera MixTransport (H2O-free mixtures): {worst_self:.1e}")

    pure_T = [300.0, 400.0, 500.0, 700.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0]
    mu = [[0.0] * len(pure_T) for _ in range(NS)]
    lam = [[0.0] * len(pure_T) for _ in range(NS)]
    for j, T in enumerate(pure_T):
        m, l_ = pure_props(T)
        for k in range(NS):
            mu[k][j] = rnd(m[k])
            lam[k][j] = rnd(l_[k])
    write_fixture("thermo_transport.json", dict(
        cantera=ct.__version__, species=SPECIES,
        canteraFitRange=[gas.min_temp, gas.max_temp],
        canteraFittingErrors={k: rnd(v) for k, v in gas.transport_fitting_errors.items()},
        pureT=pure_T, pureViscosity=mu, pureConductivity=lam,
        waterModelViscosity=[rnd(water(T)[0]) for T in pure_T],
        waterModelConductivity=[rnd(water(T)[1]) for T in pure_T],
        mixtures=cases))


if __name__ == "__main__":
    main()
