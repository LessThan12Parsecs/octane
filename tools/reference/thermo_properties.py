"""Oracle fixtures for src/physics/thermo (species, mixtures, fuels) from Cantera 3.2.

Writes:
  test/fixtures/thermo_species.json   cp/R, h/RT, s°/R of every species on a T grid (200–6000 K)
                                      + Cantera molecular weights
  test/fixtures/thermo_mixtures.json  mixture cp, cv, h, u, s (mass and molar), M, γ for air,
                                      fresh charges, complete- and equilibrium-product mixtures
  test/fixtures/thermo_fuels.json     LHV (gas, 298.15 K) per fuel, water-gas-shift K(T),
                                      WGS-equilibrium rich products, stoichiometric AFR

The Cantera phase uses the same NASA-7 data (nasa_gas.yaml) as species-data.ts, rebuilt with
the data's true 1-bar standard state (thermo_common.nasa_species; Cantera's nasa_gas.yaml
labels it 1 atm). Mixture states cover 200–6000 K and 0.02–300 bar.
Run: .venv/bin/python tools/reference/thermo_properties.py
"""
from __future__ import annotations

import os
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import P_REF_DATA, build_phase, read_species, rnd, write_fixture  # noqa: E402

SPECIES = read_species()
IDX = {s: i for i, s in enumerate(SPECIES)}
NS = len(SPECIES)
P_ATM = 101325.0
T_REF = 298.15

# Picard et al. 2008 (CIPM-2007) Table 1, renormalised over the 4 species we model
# (must match fuels.ts DRY_AIR_REFERENCE).
DRY_AIR = {"N2": 0.780848, "O2": 0.209390, "AR": 0.009332, "CO2": 0.00040}


def vec(d: dict) -> np.ndarray:
    x = np.zeros(NS)
    for k, v in d.items():
        x[IDX[k]] += v
    return x


def normalise(x: np.ndarray) -> np.ndarray:
    return x / x.sum()


def dry_air() -> np.ndarray:
    return normalise(vec(DRY_AIR))


def fuel_air(fuel: dict, phi: float, air: np.ndarray, nu_o2: float) -> np.ndarray:
    """Mole fractions of air + fuel at φ (n_fuel/n_air = φ x_O2/ν)."""
    n_f = phi * air[IDX["O2"]] / nu_o2
    return normalise(air + n_f * vec(fuel))


def main() -> None:
    gas = build_phase(transport=False)
    assert list(gas.species_names) == SPECIES
    mw = gas.molecular_weights / 1000.0  # kg/mol

    # ---------------- species ----------------
    tgrid = [200.0, 250.0, 298.15, 300.0, 350.0, 432.1, 500.0, 600.0, 700.0, 800.0, 900.0,
             999.999, 1000.0, 1000.001, 1100.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0,
             2718.28, 3000.0, 3500.0, 4000.0, 4500.0, 5000.0, 5500.0, 5999.0, 6000.0]
    cpR = [[0.0] * len(tgrid) for _ in range(NS)]
    hRT = [[0.0] * len(tgrid) for _ in range(NS)]
    sR = [[0.0] * len(tgrid) for _ in range(NS)]
    # NB: Cantera's standard_* properties of an ideal gas are evaluated at the CURRENT
    # pressure (s°(T, p) = s°(T) − R ln(p/p_ref)), so set p = p_ref (1 bar) here.
    p_ref = gas.reference_pressure
    assert p_ref == P_REF_DATA, p_ref
    for j, T in enumerate(tgrid):
        gas.TP = T, p_ref
        c, h, s = gas.standard_cp_R, gas.standard_enthalpies_RT, gas.standard_entropies_R
        for k in range(NS):
            cpR[k][j], hRT[k][j], sR[k][j] = rnd(c[k]), rnd(h[k]), rnd(s[k])
    write_fixture("thermo_species.json", dict(
        cantera=ct.__version__, species=SPECIES, T=tgrid, molarMass=[rnd(m) for m in mw],
        pRef=rnd(gas.species(0).thermo.reference_pressure), cpR=cpR, hRT=hRT, s0R=sR))

    # ---------------- mixtures ----------------
    air = dry_air()
    iso = {"IC8H18": 1.0}
    stoich_iso = fuel_air(iso, 1.0, air, 12.5)
    ch4_lean = fuel_air({"CH4": 1.0}, 0.8, air, 2.0)
    # complete-combustion products of stoichiometric iso-octane/air (lean branch exact)
    b_C = stoich_iso[IDX["IC8H18"]] * 8
    b_H = stoich_iso[IDX["IC8H18"]] * 18
    prod = stoich_iso.copy()
    prod[IDX["IC8H18"]] = 0
    prod[IDX["O2"]] = 0
    prod[IDX["CO2"]] += b_C
    prod[IDX["H2O"]] += b_H / 2
    prod = normalise(prod)

    # equilibrium products over the 12 burned-gas species (Cantera TP equilibrium)
    prod_species = SPECIES[:12]
    eq_phase = ct.Solution(thermo="ideal-gas", species=[gas.species(s) for s in prod_species])

    def equil(x: np.ndarray, T: float, p: float) -> np.ndarray:
        eq_phase.TPX = T, p, {s: float(x[IDX[s]]) for s in prod_species if x[IDX[s]] > 0}
        eq_phase.equilibrate("TP")
        out = np.zeros(NS)
        for s in prod_species:
            out[IDX[s]] = eq_phase.X[eq_phase.species_index(s)]
        return out

    def rich_seed(phi: float) -> np.ndarray:
        """Element-equivalent product seed (C→CO, H→H2; equilibrate redistributes)."""
        seed = fuel_air(iso, phi, air, 12.5)
        n_f = seed[IDX["IC8H18"]]
        seed[IDX["IC8H18"]] = 0
        seed[IDX["CO"]] += 8 * n_f
        seed[IDX["H2"]] += 9 * n_f
        return normalise(seed)

    eq_2500 = equil(prod, 2500.0, 50e5)
    eq_1500 = equil(prod, 1500.0, 1e5)
    eq_rich_2800 = equil(rich_seed(1.25), 2800.0, 30e5)
    # harder conditions: very rich, strongly dissociated, heavily diluted, water-rich
    eq_rich2_2200 = equil(rich_seed(2.0), 2200.0, 100e5)
    eq_3500_low_p = equil(prod, 3500.0, 1e4)
    lean_diluted = normalise(0.7 * fuel_air(iso, 0.3, air, 12.5) + 0.3 * prod)
    humid = normalise(air * (1 - 0.07) + 0.07 * vec({"H2O": 1.0}))

    mixtures = [
        ("dry air", air),
        ("stoichiometric iso-octane/air", stoich_iso),
        ("methane/air phi=0.8", ch4_lean),
        ("complete products iso-octane phi=1", prod),
        ("equilibrium products 2500K 50bar", eq_2500),
        ("equilibrium products 1500K 1bar", eq_1500),
        ("equilibrium products rich phi=1.25 2800K 30bar", eq_rich_2800),
        ("PRF-like blend", normalise(vec({"IC8H18": 0.01, "NC7H16": 0.003, "N2": 0.77,
                                          "O2": 0.2, "H2O": 0.01, "CO2": 0.004, "AR": 0.003,
                                          "C2H5OH": 0.001, "CH4": 0.001, "C3H8": 0.001}))),
        ("equilibrium products rich phi=2 2200K 100bar", eq_rich2_2200),
        ("equilibrium products 3500K 0.1bar (dissociated)", eq_3500_low_p),
        ("lean phi=0.3 charge + 30% products (diluted)", lean_diluted),
        ("humid air x_H2O=0.07", humid),
        ("pure H2O", vec({"H2O": 1.0})),
    ]
    states_T = [200.0, 250.0, 300.0, 500.0, 800.0, 1000.0, 1200.0, 1500.0, 2000.0, 2500.0,
                3000.0, 4000.0, 5000.0, 6000.0]
    states_p = [2e3, 1e5, 5e6, 3e7]
    cases = []
    for name, x in mixtures:
        states = []
        for T in states_T:
            for p in states_p:
                gas.TPX = T, p, x
                states.append(dict(
                    T=T, p=p,
                    cpMass=rnd(gas.cp_mass), cvMass=rnd(gas.cv_mass),
                    hMass=rnd(gas.enthalpy_mass), uMass=rnd(gas.int_energy_mass),
                    sMass=rnd(gas.entropy_mass),
                    cpMole=rnd(gas.cp_mole / 1000), hMole=rnd(gas.enthalpy_mole / 1000),
                    sMole=rnd(gas.entropy_mole / 1000),
                    molarMass=rnd(gas.mean_molecular_weight / 1000),
                    gamma=rnd(gas.cp_mass / gas.cv_mass)))
        cases.append(dict(name=name, X=[rnd(v) for v in x], states=states))
    write_fixture("thermo_mixtures.json", dict(cantera=ct.__version__, species=SPECIES, cases=cases))

    # ---------------- fuels ----------------
    R = ct.gas_constant / 1000.0
    gas.TP = T_REF, p_ref
    h298 = gas.standard_enthalpies_RT * R * T_REF  # J/mol
    comp = {s: gas.species(s).composition for s in SPECIES}
    lhv = {}
    for f in ["IC8H18", "NC7H16", "CH4", "C3H8", "C2H5OH"]:
        c = comp[f].get("C", 0.0)
        h = comp[f].get("H", 0.0)
        o = comp[f].get("O", 0.0)
        nu = c + h / 4 - o / 2
        dH = (c * h298[IDX["CO2"]] + h / 2 * h298[IDX["H2O"]]) - (h298[IDX[f]] + nu * h298[IDX["O2"]])
        lhv[f] = dict(molar=rnd(-dH), mass=rnd(-dH / mw[IDX[f]]), nuO2=nu)

    # stoichiometric AFR of iso-octane with our dry air
    air_mw = float(np.dot(air, mw))
    afr_iso = 12.5 / air[IDX["O2"]] * air_mw / mw[IDX["IC8H18"]]

    # water-gas shift constant CO2 + H2 = CO + H2O
    wgs = []
    for T in [1000.0, 1500.0, 1740.0, 2000.0, 2500.0]:
        gas.TP = T, p_ref
        g = gas.standard_gibbs_RT
        K = np.exp(-(g[IDX["CO"]] + g[IDX["H2O"]] - g[IDX["CO2"]] - g[IDX["H2"]]))
        wgs.append(dict(T=T, K=rnd(K)))

    # rich complete-combustion products = WGS equilibrium at 1740 K over {N2, AR, CO2, H2O, CO, H2}
    wgs_species = ["N2", "AR", "CO2", "H2O", "CO", "H2"]
    wgs_phase = ct.Solution(thermo="ideal-gas", species=[gas.species(s) for s in wgs_species])
    rich_cases = []
    for phi in [1.05, 1.25, 1.5, 2.0, 2.5]:
        reac = fuel_air(iso, phi, air, 12.5)
        # seed with the right element totals (the solver only needs those)
        nf = reac[IDX["IC8H18"]]
        bC = 8 * nf + reac[IDX["CO2"]]
        bH = 18 * nf
        bO = 2 * reac[IDX["O2"]] + 2 * reac[IDX["CO2"]]
        seed = np.zeros(NS)
        seed[IDX["N2"]] = reac[IDX["N2"]]
        seed[IDX["AR"]] = reac[IDX["AR"]]
        if bO - bC <= bH / 2:
            seed[IDX["CO"]] = bC
            seed[IDX["H2O"]] = bO - bC
            seed[IDX["H2"]] = bH / 2 - (bO - bC)
        else:
            seed[IDX["H2O"]] = bH / 2
            seed[IDX["CO2"]] = bO - bC - bH / 2
            seed[IDX["CO"]] = bC - seed[IDX["CO2"]]
        assert seed.min() >= 0
        wgs_phase.TPX = 1740.0, P_ATM, {s: float(seed[IDX[s]]) for s in wgs_species}
        wgs_phase.equilibrate("TP")
        out = np.zeros(NS)
        for s in wgs_species:
            out[IDX[s]] = wgs_phase.X[wgs_phase.species_index(s)]
        rich_cases.append(dict(phi=phi, reactantsX=[rnd(v) for v in reac],
                               productsX=[rnd(v) for v in out]))

    write_fixture("thermo_fuels.json", dict(
        cantera=ct.__version__, species=SPECIES, lhvGas=lhv, dryAirX=[rnd(v) for v in air],
        stoichAfrIsooctaneDryAir=rnd(afr_iso), wgsK=wgs, richProducts1740=rich_cases))


if __name__ == "__main__":
    main()
