"""Oracle for src/physics/chemistry/zeldovich.ts (extended-Zeldovich thermal NO).

Writes test/fixtures/chemistry_zeldovich.json with three parts:

rates      Equilibrium burned-gas states (T 1800-2800 K, p 10-80 bar, phi 0.8-1.2, PRF 90/air,
           project species + nasa_gas.yaml thermo, Cantera equilibrate('TP')). For each state
           and alpha = [NO]/[NO]e in {0, 0.5}:
             - Cantera net NO and N production rates of the three Zeldovich reactions with
               GRI-Mech 3.0 rate parameters on the SAME (NASA) thermo as the TS code, at the
               composition X = Xeq except NO = alpha NO_e and N = its quasi-steady value
               (-> must equal Heywood eq. 11.8 exactly), and at N = N_e (tests the full
               forward/reverse form);
             - the same with the GRI-3.0 thermo (zeldovich_gri30.yaml; reported only);
             - the Heywood/Hanson-Salimian rate set evaluated the same way (reported only).
stress     Same construction (Heywood eq. 11.8 state, N quasi-steady) far outside the nominal
           range: T 1200-3400 K, p 1-250 bar, phi 0.4-2.0 (+ 25 % complete-combustion EGR at
           phi = 1), alpha in {0, 0.5, 2} (alpha = 2: decomposition, negative rate).
evolution  Constant-(T, p) NO(t) from NO-free equilibrium burned gas:
             - rate-controlled model (Heywood eq. 11.8, O/OH/H/O2/N2 at equilibrium,
               integrated exactly);
             - full GRI-Mech 3.0 (53 species, 325 reactions; all O/H radicals, N2O and NNH
               routes, NO2) in a constant-T,p Cantera reactor started from the GRI-3.0
               equilibrium composition with every N-containing species except N2 removed.
"""
from __future__ import annotations

import math

import cantera as ct
import numpy as np

import os

from chemistry_common import FIXTURES, ZELDOVICH_YAML, fresh_charge, write_json

PROJ = ["N2", "O2", "AR", "CO2", "H2O", "CO", "H2", "OH", "H", "O", "NO", "N",
        "IC8H18", "NC7H16", "CH4", "C3H8", "C2H5OH"]
NASA_NAME = {"AR": "Ar", "IC8H18": "C8H18,isooctane", "NC7H16": "C7H16,n-heptane"}
N_EQ = 12

# GRI-Mech 3.0 (reactions 178-180), exothermic direction, cm3/mol/s, cal/mol
GRI = [("N + NO <=> N2 + O", 2.7e13, 0.0, 355.0),
       ("N + O2 <=> NO + O", 9.0e9, 1.0, 6500.0),
       ("N + OH <=> NO + H", 3.36e13, 0.0, 385.0)]
# Heywood (1988) Table 11.1 / Hanson & Salimian (1984), same directions, cm3/mol/s, Ta in K
HEYWOOD = [("N + NO <=> N2 + O", 1.6e13, 0.0, 0.0),
           ("N + O2 <=> NO + O", 6.4e9, 1.0, 3150.0),
           ("N + OH <=> NO + H", 4.1e13, 0.0, 0.0)]


def nasa_phase(rset=None, ta_units=False):
    """Project species (nasa_gas thermo), optionally with the 3 Zeldovich reactions."""
    allsp = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    sps = []
    for n in PROJ:
        src = allsp[NASA_NAME.get(n, n)]
        s = ct.Species(n, src.composition)
        s.thermo = src.thermo
        sps.append(s)
    if rset is None:
        return ct.Solution(thermo="ideal-gas", species=sps)
    rx = []
    for eq, A, b, E in rset:
        Ea = E * ct.gas_constant if ta_units else E * 4184.0  # J/kmol
        rx.append(ct.Reaction(equation=eq, rate=ct.ArrheniusRate(A * 1e-3, b, Ea)))
    return ct.Solution(thermo="ideal-gas", kinetics="gas", species=sps, reactions=rx)


def equilibrium_state(eq_gas, T, p, phi):
    X0 = fresh_charge(phi, 90.0)
    eq_gas.TPX = T, p, X0
    eq_gas.equilibrate("TP")
    X = np.array(eq_gas.X)
    X[N_EQ:] = 0.0
    return X / X.sum()


def qss_concentration_state(gas, T, p, Xeq, alpha):
    """Concentrations (kmol/m3) with every species at its EQUILIBRIUM concentration
    c x_eq (c = p/RT) except NO = alpha [NO]e and N at its quasi-steady value -- exactly the
    assumptions of Heywood eq. 11.8. Set on `gas` via the concentrations setter (T kept)."""
    gas.TPX = T, p, Xeq
    c = p / (ct.gas_constant * T)
    C = np.array(Xeq) * c
    iNO, iN = PROJ.index("NO"), PROJ.index("N")
    C[iNO] = alpha * Xeq[iNO] * c
    kf = gas.forward_rate_constants
    kr = gas.reverse_rate_constants
    i = PROJ.index
    num = kr[0] * C[i("O")] * C[i("N2")] + kr[1] * C[iNO] * C[i("O")] + kr[2] * C[iNO] * C[i("H")]
    den = kf[0] * C[iNO] + kf[1] * C[i("O2")] + kf[2] * C[i("OH")]
    C[iN] = num / den
    gas.TP = T, p
    gas.concentrations = C
    assert abs(gas.T - T) < 1e-9
    return C


def rates_at(gas, T, p, X):
    gas.TPX = T, p, X
    w = gas.net_production_rates  # kmol/m3/s
    return w[gas.species_index("NO")] * 1e3, w[gas.species_index("N")] * 1e3  # mol/m3/s


def rates_now(gas):
    w = gas.net_production_rates
    return w[gas.species_index("NO")] * 1e3, w[gas.species_index("N")] * 1e3


def heywood_rate(T, p, Xeq, alpha, rset, ta_units):
    """Heywood eq. 11.8 in Python (reference for the report)."""
    c = p / (8.31446261815324 * T)  # mol/m3
    ks = []
    for _, A, b, E in rset:
        Ta = E if ta_units else E * 4184.0 / 8314.46261815324
        ks.append(A * 1e-6 * T**b * math.exp(-Ta / T))  # m3/mol/s
    N = Xeq[PROJ.index("N")] * c
    R1 = ks[0] * N * Xeq[PROJ.index("NO")] * c
    R2 = ks[1] * N * Xeq[PROJ.index("O2")] * c
    R3 = ks[2] * N * Xeq[PROJ.index("OH")] * c
    return 2 * R1 * (1 - alpha**2) / (1 + alpha * R1 / (R2 + R3)), R1, R2, R3


def alpha_exact(s, K):
    """alpha(t) of d alpha/dt = c (1 - a^2)/(1 + K a), alpha(0) = 0, at s = c t:
    F(a) = -(1+K)/2 ln(1-a) + (1-K)/2 ln(1+a) = s  (bisection on [0, 1))."""
    lo, hi = 0.0, 1.0 - 1e-16
    for _ in range(200):
        m = 0.5 * (lo + hi)
        F = -(1 + K) / 2 * math.log(1 - m) + (1 - K) / 2 * math.log(1 + m)
        if F < s:
            lo = m
        else:
            hi = m
    return 0.5 * (lo + hi)


def main():
    eq_gas = nasa_phase()
    zg = nasa_phase(GRI)
    zh = nasa_phase(HEYWOOD, ta_units=True)
    zgri = ct.Solution(ZELDOVICH_YAML)

    Ts = [1800.0, 2000.0, 2200.0, 2400.0, 2600.0, 2800.0]
    ps = [10e5, 20e5, 40e5, 80e5]
    phis = [0.8, 0.9, 1.0, 1.1, 1.2]
    iNO, iN = PROJ.index("NO"), PROJ.index("N")
    cases = []
    states = []
    worst_gri_thermo = 0.0
    ratio_heywood = []
    for T in Ts:
        for p in ps:
            for phi in phis:
                Xeq = equilibrium_state(eq_gas, T, p, phi)
                Xeq = np.array([float(f"{v:.12e}") for v in Xeq])  # exactly what TS reads
                states.append(dict(T=T, p=p, phi=phi, Xeq=[float(v) for v in Xeq[:N_EQ]]))
                for alpha in (0.0, 0.5):
                    # X_full: equilibrium with NO scaled by alpha (N = N_e), renormalised.
                    Xn = Xeq.copy()
                    Xn[iNO] = alpha * Xeq[iNO]
                    Xn /= Xn.sum()
                    wNO_f, wN_f = rates_at(zg, T, p, Xn)
                    # Heywood state: equilibrium concentrations, NO = alpha NO_e, N quasi-steady
                    Cq = qss_concentration_state(zg, T, p, Xeq, alpha)
                    wNO_q, wN_q = rates_now(zg)
                    xNq = Cq[iN] / (p / (ct.gas_constant * T))
                    # GRI thermo (reverse rates from GRI-3.0 NASA fits), same concentrations
                    gidx = [zgri.species_index(n) for n in PROJ[:N_EQ]]
                    Cg = np.zeros(zgri.n_species)
                    Cg[gidx] = Cq[:N_EQ]
                    zgri.TPX = T, p, {n: max(v, 1e-30) for n, v in zip(PROJ[:N_EQ], Xeq[:N_EQ])}
                    zgri.concentrations = Cg
                    wNO_gri = zgri.net_production_rates[zgri.species_index("NO")] * 1e3
                    worst_gri_thermo = max(worst_gri_thermo, abs(wNO_gri / wNO_q - 1))
                    # Heywood rate set, same construction
                    qss_concentration_state(zh, T, p, Xeq, alpha)
                    wNO_h, _ = rates_now(zh)
                    ratio_heywood.append(wNO_h / wNO_q)
                    cases.append(dict(
                        s=len(states) - 1, alpha=alpha,
                        wNOqss=float(f"{wNO_q:.12e}"), wNqss=float(f"{wN_q:.3e}"),
                        xNqss=float(f"{xNq:.12e}"),
                        wNOfull=float(f"{wNO_f:.12e}"), wNfull=float(f"{wN_f:.12e}"),
                        wNOgriThermo=float(f"{wNO_gri:.8e}"),
                        wNOheywood=float(f"{wNO_h:.8e}"),
                    ))
    print(f"rates: {len(cases)} cases; GRI-thermo vs NASA-thermo max rel diff {worst_gri_thermo:.3e}")

    # ------------------------------------------------------------ stress (extreme states)
    stress = []
    for T in (1200.0, 1500.0, 3000.0, 3400.0):
        for p in (1e5, 5e5, 150e5, 250e5):
            for phi, egr in ((0.4, 0.0), (0.6, 0.0), (1.0, 0.25), (1.5, 0.0), (2.0, 0.0)):
                X0 = fresh_charge(phi, 90.0)
                if egr > 0:
                    from chemistry_common import complete_combustion_products
                    res = complete_combustion_products(X0)
                    X0 = {k: (1 - egr) * X0.get(k, 0.0) + egr * res.get(k, 0.0) for k in set(X0) | set(res)}
                eq_gas.TPX = T, p, X0
                eq_gas.equilibrate("TP")
                Xeq = np.array(eq_gas.X)
                Xeq[N_EQ:] = 0.0
                Xeq = Xeq / Xeq.sum()
                Xeq = np.array([float(f"{v:.12e}") for v in Xeq])
                for alpha in (0.0, 0.5, 2.0):
                    qss_concentration_state(zg, T, p, Xeq, alpha)
                    wNO_q, wN_q = rates_now(zg)
                    stress.append(dict(T=T, p=p, phi=phi, egr=egr, alpha=alpha,
                                       Xeq=[float(v) for v in Xeq[:N_EQ]],
                                       wNOqss=float(f"{wNO_q:.12e}")))
    print(f"stress: {len(stress)} cases")
    rh = np.array(ratio_heywood)
    print(f"Heywood/Hanson-Salimian vs GRI-3.0 rate ratio: min {rh.min():.3f} max {rh.max():.3f}")

    # ------------------------------------------------------------ time evolution
    gri = ct.Solution("gri30.yaml")
    evo = []
    nsp_keep_zero = [s for s in gri.species_names if "N" in gri.species(s).composition and s != "N2"]
    for (T, p, phi) in [(2400.0, 40e5, 1.0), (2600.0, 60e5, 0.9), (2200.0, 20e5, 0.8),
                        (2800.0, 80e5, 1.1), (2000.0, 40e5, 1.0), (2500.0, 50e5, 1.2)]:
        # our model: NASA-thermo equilibrium of the PRF 90 products
        Xeq = equilibrium_state(eq_gas, T, p, phi)
        # GRI equilibrium of the same element composition, then remove N species except N2
        fc = fresh_charge(phi, 90.0)
        # GRI-3.0 has no PRF: use the element-equivalent CO2/H2O/O2/N2/AR products
        # (identical element vector) to equilibrate with GRI thermo.
        eq_gas.TPX = T, p, fc
        eq_gas.equilibrate("TP")
        Xg = {("AR" if n == "AR" else n): v for n, v in zip(PROJ, eq_gas.X) if v > 0 and n in gri.species_names}
        gri.TPX = T, p, Xg
        gri.equilibrate("TP")
        XeqG = gri.X.copy()
        xNOeG = XeqG[gri.species_index("NO")]
        X0 = XeqG.copy()
        for s in nsp_keep_zero:
            X0[gri.species_index(s)] = 0.0
        gri.TPX = T, p, X0
        r = ct.IdealGasConstPressureReactor(gri, energy="off", clone=False)
        net = ct.ReactorNet([r])
        net.rtol, net.atol = 1e-10, 1e-22
        # model characteristic time (Heywood eq. 11.10): tau_NO = [NO]e / (4 R1/(1+K))
        w0, R1, R2, R3 = heywood_rate(T, p, Xeq, 0.0, GRI, False)
        c = p / (8.31446261815324 * T)
        K = R1 / (R2 + R3)
        tauNO = Xeq[iNO] * c / (4 * R1 / (1 + K))
        times = [tauNO * f for f in (0.01, 0.03, 0.1, 0.3, 1.0, 3.0)]
        gri_no = []
        gri_no_frac = []
        for t in times:
            net.advance(t)
            gri_no.append(float(r.phase.X[gri.species_index("NO")]))
            # also report total NOx-N (NO + NO2 + N2O*2 ...) ~ NO here
            gri_no_frac.append(float(r.phase.X[gri.species_index("NO")] / xNOeG))
        model_alpha = [alpha_exact(t * (1 + K) / (2 * tauNO), K) for t in times]
        evo.append(dict(T=T, p=p, phi=phi, Xeq=[float(f"{v:.12e}") for v in Xeq[:N_EQ]],
                        tauNO=tauNO, times=times, xNOeGri=float(xNOeG),
                        griNO=gri_no, griAlpha=gri_no_frac, modelAlpha=model_alpha))
        print(f"evolution T={T} p={p/1e5}bar phi={phi}: tau_NO={tauNO*1e3:.3f} ms, "
              f"NOe(NASA)={Xeq[iNO]*1e6:.0f} ppm, NOe(GRI)={xNOeG*1e6:.0f} ppm, "
              f"GRI alpha(t)={[round(a, 4) for a in gri_no_frac]} model={[round(a, 4) for a in model_alpha]}")

    write_json(os.path.join(FIXTURES, "chemistry_zeldovich.json"), dict(
        species=PROJ,
        gri30=[dict(eq=e, A_cm3=A, b=b, Ea_cal=E) for e, A, b, E in GRI],
        heywood=[dict(eq=e, A_cm3=A, b=b, Ta=E) for e, A, b, E in HEYWOOD],
        states=states, rates=cases, evolution=evo, stress=stress,
        summary=dict(griThermoMaxRelDiff=worst_gri_thermo,
                     heywoodOverGriMin=float(rh.min()), heywoodOverGriMax=float(rh.max())),
    ))


if __name__ == "__main__":
    main()
