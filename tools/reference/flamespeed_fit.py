"""Build the detailed-chemistry laminar-burning-velocity tables used by
src/physics/combustion/laminar-flame-speed.ts from the Cantera FreeFlame oracles, fit the
overall activation energies, and generate transport reference data.

Tables: ln(S_L) on the grid φ ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6} ×
Tu ∈ {300, 500, 700, 900} K × p ∈ {1, 5, 20, 60} bar, undiluted fuel–air (O2 + 3.76 N2):
  CH4           GRI-Mech 3.0            (flamespeed_cantera.json)
  IC8H18, NC7H16 PRF_HT (Jerzembeck 2009) (flamespeed_prf_cantera.json)
Grid points whose flame did not converge are filled from a per-φ least-squares fit of
ln S_L = c0 + c1 θ + c2 π + c3 θπ + c4 θ² + c5 π² (θ = ln Tu/300 K, π = ln p/1 bar) to the
other points at that φ (listed in the fixture under "filled").

Activation energies of the residual-dilution model (also used for the Zeldovich number of the
Markstein model): E_a(φ, p, x) = −2R ln[S_L(x)/S_L(0)] / (1/T_ad(x) − 1/T_ad(0)) from the diluted
flames of flamespeed_dilution.py (fit set) and, for ethanol, from the residual correction of
Vancoillie et al. (2012).

Ethanol S_L table: Vancoillie et al. (2012) correlation (see vancoillie_sl), power-law continued
below its 5-bar range limit.

Transport reference: Cantera mixture-averaged diffusion coefficients / thermal diffusivity
for fresh charges using exactly our species data (thermo_common.build_phase).

Writes test/fixtures/flamespeed_fit.json and rewrites the generated blocks of
laminar-flame-speed.ts. Run: .venv/bin/python tools/reference/flamespeed_fit.py
"""
from __future__ import annotations

import json
import math
import os
import re
import sys

import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
FIX = os.path.join(ROOT, "test", "fixtures")
OUT = os.path.join(FIX, "flamespeed_fit.json")
TS = os.path.join(ROOT, "src", "physics", "combustion", "laminar-flame-speed.ts")
R = 8.31446261815324

PHI = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]
TU = [300.0, 500.0, 700.0, 900.0]
P = [1e5, 5e5, 20e5, 60e5]
FUELS = ["IC8H18", "NC7H16", "CH4", "C3H8", "C2H5OH"]  # FUEL_SPECIES order


def load(name):
    path = os.path.join(FIX, name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)["flames"]


def fuel_key(r):
    f = r["fuel"]
    if isinstance(f, dict):
        return next(iter(f)) if len(f) == 1 else None
    return f


def key(phi, Tu, p):
    return (round(phi, 4), round(Tu, 1), round(p / 1e5, 3))


def build_table(rows):
    have = {key(r["phi"], r["Tu"], r["p"]): math.log(r["SL"]) for r in rows}
    tab, filled = [], []
    for phi in PHI:
        pts = [(Tu, p, have[key(phi, Tu, p)]) for Tu in TU for p in P if key(phi, Tu, p) in have]
        coef = None
        if len(pts) < 16:
            if len(pts) < 8:
                raise RuntimeError(f"too few flames at phi={phi}: {len(pts)}")
            A = np.array([[1, math.log(T / 300), math.log(p / 1e5), math.log(T / 300) * math.log(p / 1e5),
                           math.log(T / 300) ** 2, math.log(p / 1e5) ** 2] for T, p, _ in pts])
            coef, *_ = np.linalg.lstsq(A, np.array([v for *_, v in pts]), rcond=None)
        for Tu in TU:
            for p in P:
                k = key(phi, Tu, p)
                if k in have:
                    tab.append(have[k])
                else:
                    th, pi = math.log(Tu / 300), math.log(p / 1e5)
                    tab.append(float(np.dot(coef, [1, th, pi, th * pi, th * th, pi * pi])))
                    filled.append(dict(phi=phi, Tu=Tu, p=p))
    return tab, filled


# ---------------------------------------------------------------------------------------
# Effective activation energy of the dilution model, E_a(φ, p, x)
# ---------------------------------------------------------------------------------------
EA_PHI = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4]  # = the S_L-table φ nodes in 0.6-1.4
EA_TP = [(300.0, 1e5), (500.0, 5e5), (700.0, 20e5), (900.0, 60e5)]  # (Tu, p) of the fitted series
EA_X = [0.15, 0.3]  # = X_DIL[1:], the diluted nodes of the TypeScript tables
# Ethanol (no mechanism oracle): Vancoillie et al. (2012) residual correction F, Eqs. 5-6, inside
# its fitted range (400-900 K, 5-85 bar); the 1-bar row repeats the 5-bar row.
EA_TP_ETHANOL = [(500.0, 5e5), (500.0, 5e5), (700.0, 20e5), (850.0, 60e5)]


def diluent_volume_fraction(fuel, phi, x):
    """Mole fraction of the complete-combustion-product diluent at mass fraction x (air O2 + 3.76 N2)."""
    N = {"O2": 1.0, "N2": 3.76, fuel: phi / NU[fuel]}
    nd = x / (1 - x) * sum(complete_products(N).values())
    return nd / (sum(N.values()) + nd)


def ea_tables():
    """E_a = −2R ln(S_L(x)/S_L(0)) / (1/T_ad(x) − 1/T_ad(0)) at every (x, φ, p) of the fitting set
    of flamespeed_dilution.py (T_ad: HP equilibrium with OUR thermo and the dry-air fresh charge
    of the tables, i.e. exactly the T_ad the TypeScript model interpolates). One value per diluted
    node x ∈ EA_X, so the model reproduces both diluted flames of every fitted series exactly and
    the T_ad dependence carries the thermal effect to other (Tu, p). Missing (failed) flames are
    filled from the nearest φ of the same (x, p)."""
    path = os.path.join(FIX, "flamespeed_dilution.json")
    with open(path, encoding="utf-8") as fh:
        flames = [r for r in json.load(fh)["flames"] if r["set"] == "fit"]
    idx = {}
    for r in flames:
        idx[(r["fuel"], round(r["phi"], 4), r["Tu"], round(r["p"]), r["ydil"])] = r["SL"]
    tables, diag = {}, {}
    for fu in FUELS:
        tab = np.full((len(EA_X), len(EA_PHI), len(EA_TP)), np.nan)
        for ix, x in enumerate(EA_X):
            for i, phi in enumerate(EA_PHI):
                for k in range(len(EA_TP)):
                    if fu == "C2H5OH":
                        Tu, p = EA_TP_ETHANOL[k]
                        ratio = vancoillie_F(phi, Tu, p / 1e5, diluent_volume_fraction(fu, phi, x))
                    else:
                        Tu, p = EA_TP[k]
                        s0 = idx.get((fu, phi, Tu, round(p), 0.0))
                        sx = idx.get((fu, phi, Tu, round(p), x))
                        if not (s0 and sx):
                            continue
                        ratio = sx / s0
                    T0 = adiabatic_T(fu, phi, Tu, p, 0.0)
                    Tx = adiabatic_T(fu, phi, Tu, p, x)
                    tab[ix, i, k] = -2 * R * math.log(ratio) / (1 / Tx - 1 / T0)
        missing = [(EA_X[a], EA_PHI[b], EA_TP[c][1]) for a, b, c in zip(*np.where(np.isnan(tab)))]
        if len(missing) > 6:
            raise RuntimeError(f"{fu}: too many missing diluted flames: {missing}")
        for a, b, c in zip(*np.where(np.isnan(tab))):
            cand = [(abs(j - b), tab[a, j, c]) for j in range(len(EA_PHI)) if not np.isnan(tab[a, j, c])]
            tab[a, b, c] = min(cand)[1]
        tables[fu] = [[[float(tab[a, b, c]) for c in range(len(EA_TP))] for b in range(len(EA_PHI))]
                      for a in range(len(EA_X))]
        diag[fu] = dict(missing=missing)
        print(f"E_a {fu} [kJ/mol], rows φ {EA_PHI}, columns p [bar] {[tp[1] / 1e5 for tp in EA_TP]}:")
        for a, x in enumerate(EA_X):
            print(f"  x = {x}: " + "; ".join(", ".join(f"{v / 1e3:.0f}" for v in tab[a, b]) for b in range(len(EA_PHI))))
    return tables, diag


def transport_reference():
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from thermo_common import build_phase, read_species  # noqa: E402
    species = read_species()
    gas = build_phase(transport=True, trange=(300.0, 3000.0))
    air = {"N2": 0.780848, "O2": 0.209390, "AR": 0.009332, "CO2": 0.00040}
    s = sum(air.values())
    air = {k: v / s for k, v in air.items()}
    nu = {"IC8H18": 12.5, "NC7H16": 11.0, "CH4": 2.0, "C3H8": 5.0, "C2H5OH": 3.0}
    out = []
    for fuel in FUELS:
        for phi in (0.7, 1.0, 1.3):
            for Tu, p in ((300.0, 1e5), (600.0, 20e5), (900.0, 60e5)):
                X = dict(air)
                X[fuel] = phi * air["O2"] / nu[fuel]
                tot = sum(X.values())
                x = np.zeros(len(species))
                for k, v in X.items():
                    x[species.index(k)] = v / tot
                gas.TPX = Tu, p, x
                D = gas.mix_diff_coeffs
                Dbin = gas.binary_diff_coeffs
                iF, iO, iN = species.index(fuel), species.index("O2"), species.index("N2")
                out.append(dict(
                    fuel=fuel, phi=phi, Tu=Tu, p=p,
                    alpha=float(gas.thermal_conductivity / (gas.density * gas.cp_mass)),
                    Dfuel=float(D[iF]), DO2=float(D[iO]),
                    DfuelN2=float(Dbin[iF, iN]), DO2N2=float(Dbin[iO, iN]),
                ))
    return out


# ---------------------------------------------------------------------------------------
# Adiabatic flame temperatures with OUR thermo (Cantera HP equilibrium over our species set)
# ---------------------------------------------------------------------------------------
X_DIL = [0.0, 0.15, 0.3]
DRY_AIR = {"N2": 0.780848, "O2": 0.209390, "AR": 0.009332, "CO2": 0.00040}
NU = {"IC8H18": 12.5, "NC7H16": 11.0, "CH4": 2.0, "C3H8": 5.0, "C2H5OH": 3.0}
ATOMS = {"N2": (0, 0, 0, 2, 0), "O2": (0, 0, 2, 0, 0), "AR": (0, 0, 0, 0, 1), "CO2": (1, 0, 2, 0, 0),
         "H2O": (0, 2, 1, 0, 0), "CO": (1, 0, 1, 0, 0), "H2": (0, 2, 0, 0, 0),
         "IC8H18": (8, 18, 0, 0, 0), "NC7H16": (7, 16, 0, 0, 0), "CH4": (1, 4, 0, 0, 0),
         "C3H8": (3, 8, 0, 0, 0), "C2H5OH": (2, 6, 1, 0, 0)}
_phase = None


def phase():
    global _phase
    if _phase is None:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from thermo_common import build_phase  # noqa: E402
        _phase = build_phase()
    return _phase


def complete_products(N):
    """Complete-combustion products (moles) of reactant moles N — mirrors thermo/fuels.ts
    completeCombustionMoles (lean: CO2, H2O, O2; rich: water-gas shift frozen at 1740 K)."""
    C = sum(v * ATOMS[k][0] for k, v in N.items())
    H = sum(v * ATOMS[k][1] for k, v in N.items())
    O = sum(v * ATOMS[k][2] for k, v in N.items())
    Nn = sum(v * ATOMS[k][3] for k, v in N.items())
    Ar = sum(v * ATOMS[k][4] for k, v in N.items())
    out = {"N2": Nn / 2, "AR": Ar}
    if O >= 2 * C + H / 2:
        out.update(CO2=C, H2O=H / 2, O2=(O - 2 * C - H / 2) / 2)
    else:
        g = phase()
        g.TP = 1740.0, 1e5
        gRT = dict(zip(g.species_names, g.standard_gibbs_RT))
        K = math.exp(-(gRT["CO"] + gRT["H2O"] - gRT["CO2"] - gRT["H2"]))
        A = H / 2 - O + C
        D = max(0.0, O - C)
        B = K * A + C + D
        CD = C * D
        a = 2 * CD / (B + math.sqrt(max(0.0, B * B + 4 * (K - 1) * CD))) if CD > 0 else 0.0
        a = min(max(a, 0.0), min(C, D))
        out.update(CO2=a, CO=C - a, H2O=max(0.0, D - a), H2=max(0.0, A + a))
    return out


def fresh_charge(fuel, phi, x, air=None, stoich_products=False):
    """Mole fractions: dry air (or `air`) + fuel at φ + diluent of MASS fraction x (complete
    products of the same mixture, or of the stoichiometric one)."""
    air = air or DRY_AIR
    s = sum(air.values())
    N = {k: v / s for k, v in air.items()}
    N[fuel] = N.get(fuel, 0.0) + phi * N["O2"] / NU[fuel]
    if x > 0:
        base = dict(N)
        if stoich_products:
            base = {k: v / s for k, v in air.items()}
            base[fuel] = base["O2"] / NU[fuel]
        prod = complete_products(base)
        g = phase()
        M = dict(zip(g.species_names, g.molecular_weights))
        m_fresh = sum(v * M[k] for k, v in N.items())
        m_prod = sum(v * M[k] for k, v in prod.items())
        scale = x / (1 - x) * m_fresh / m_prod
        for k, v in prod.items():
            N[k] = N.get(k, 0.0) + scale * v
    tot = sum(N.values())
    return {k: v / tot for k, v in N.items() if v > 0}


def adiabatic_T(fuel, phi, Tu, p, x, **kw):
    g = phase()
    g.TPX = Tu, p, fresh_charge(fuel, phi, x, **kw)
    g.equilibrate("HP")
    return float(g.T)


def tad_tables():
    out = {}
    for fu in FUELS:
        vals = []
        for phi in PHI:
            for Tu in TU:
                for p in P:
                    vals.append([math.log(adiabatic_T(fu, phi, Tu, p, x)) for x in X_DIL])
        out[fu] = vals  # [node][ix]
    return out


# ---------------------------------------------------------------------------------------
# Ethanol: Vancoillie, Demuynck, Galle, Verhelst & van Oijen, "A laminar burning velocity and
# flame thickness correlation for ethanol-air mixtures valid at spark-ignition engine
# conditions", Fuel 102 (2012) 460-469 (Ghent University preprint, 9 May 2012):
#   u_l = u_l0(φ,p) (Tu/T0)^α(φ,p) F(φ,Tu,p,f),  p0 = 1 bar, T0 = 300 K, u_l in cm/s  (Eq. 2)
#   α     = a1 + a2 φ + a3 P + a4 φ P + a5 φ² + a6 P² + a7 φ³ + a8 P³ + a9 φ² P + a10 φ² P²
#           + a11/φ + a12 P/φ                                                        (Eq. 3)
#   ln u_l0 = b1 + b2 φ + b3 P + b4 φ P + b5 φ² + b6 P² + b7 φ³ + b8 P³ + b9 φ² P
#           + b10 φ² P² + b11/φ                                                      (Eq. 4)
#   F1 = c1 + c2 φ + c3 t + c4 P + c5 f + c6 φ² + c7 t² + c8 f² + c9 φ t + c10 φ P + c11 φ f
#        + c12 t f + c13 φ³ + c14 f³ + c15 t φ² + c16 P φ² + c17 φ t² + c18 f t² + c19 φ f²
#        + c20 t f²,   F = min(1, F1)                                                (Eqs. 5-6)
# with P = p/p0, t = Tu/T0, f = VOLUME fraction of residuals. Coefficients: Tables 4 and 5 of
# the preprint (transcribed from the text layer; equation forms read from the typeset pages).
# Fitted to ~1500 flames computed with the Li et al. (2007) ethanol mechanism at 400-900 K,
# 5-85 bar, φ 0.5-2, f 0-0.5 (mean |residual| 7.4 %); the authors note that the cubic-in-p
# form under-predicts u_l below ~10 bar.
VANC_A = [3.717600e+00, -9.398400e+00, 3.980000e-02, -1.860000e-02, 8.413800e+00, -2.832200e-04,
          -2.055000e+00, 1.401100e-06, -8.349800e-04, 4.319800e-05, 1.332500e+00, -6.523800e-03]
VANC_B = [-3.123300e+00, 2.054070e+01, -5.880000e-02, 9.329600e-03, -1.617230e+01, 7.285600e-04,
          3.633400e+00, -3.097700e-06, 4.372200e-03, -1.179500e-04, -1.434600e+00]
VANC_C = [1.087600e+00, 1.088600e+00, -4.133000e-01, -2.787700e-03, -6.703000e+00, -7.413000e-01,
          1.250000e-01, 8.049200e+00, -2.486000e-01, 3.037800e-03, 5.054000e-01, 1.851200e+00,
          1.803000e-01, -2.108700e+00, 1.680000e-02, -5.028100e-04, 4.040000e-02, -2.813000e-01,
          -7.175000e-01, -1.323600e+00]


def vancoillie_alpha(phi, P):
    a = VANC_A
    return (a[0] + a[1] * phi + a[2] * P + a[3] * phi * P + a[4] * phi ** 2 + a[5] * P ** 2 + a[6] * phi ** 3
            + a[7] * P ** 3 + a[8] * phi ** 2 * P + a[9] * phi ** 2 * P ** 2 + a[10] / phi + a[11] * P / phi)


def vancoillie_ln_ul0(phi, P):
    b = VANC_B
    return (b[0] + b[1] * phi + b[2] * P + b[3] * phi * P + b[4] * phi ** 2 + b[5] * P ** 2 + b[6] * phi ** 3
            + b[7] * P ** 3 + b[8] * phi ** 2 * P + b[9] * phi ** 2 * P ** 2 + b[10] / phi)


def vancoillie_F(phi, Tu, P, f):
    c = VANC_C
    t = Tu / 300.0
    F1 = (c[0] + c[1] * phi + c[2] * t + c[3] * P + c[4] * f + c[5] * phi ** 2 + c[6] * t ** 2 + c[7] * f ** 2
          + c[8] * phi * t + c[9] * phi * P + c[10] * phi * f + c[11] * t * f + c[12] * phi ** 3 + c[13] * f ** 3
          + c[14] * t * phi ** 2 + c[15] * P * phi ** 2 + c[16] * phi * t ** 2 + c[17] * f * t ** 2
          + c[18] * phi * f ** 2 + c[19] * t * f ** 2)
    return min(1.0, F1)


def vancoillie_sl(phi, Tu, p, f=0.0):
    """Ethanol-air laminar burning velocity, m/s (Vancoillie et al. 2012, Eqs. 2-6). F ≡ 1 for
    f = 0: F is the residual correction (fitted as diluted/undiluted ratio); its polynomial F1
    scatters by −12/+9 % around 1 at f = 0, which the min(1, F1) clip would turn into a spurious
    reduction of undiluted flames."""
    P = p / 1e5
    F = vancoillie_F(phi, Tu, P, f) if f > 0 else 1.0
    return math.exp(vancoillie_ln_ul0(phi, P)) * (Tu / 300.0) ** vancoillie_alpha(phi, P) * F / 100.0


VANC_P_MIN = 5e5  # lower end of the correlation's fitted pressure range (5-85 bar), Pa


def vancoillie_sl_ext(phi, Tu, p):
    """Undiluted ethanol S_L, m/s: the correlation for p ≥ 5 bar; below, its local power law
    at 5 bar (C¹ continuation, β = ∂ln u_l/∂ln p at 5 bar) instead of the cubic-in-p polynomial,
    which the authors show to under-predict at low pressure (at 300 K / 1 bar, φ = 1:
    0.299 m/s from the polynomial, 0.348 m/s continued; measured ≈ 0.36-0.40 m/s)."""
    if p >= VANC_P_MIN:
        return vancoillie_sl(phi, Tu, p)
    h = 1e-4
    s0 = vancoillie_sl(phi, Tu, VANC_P_MIN)
    beta = (math.log(vancoillie_sl(phi, Tu, VANC_P_MIN * (1 + h))) - math.log(s0)) / math.log(1 + h)
    return s0 * (p / VANC_P_MIN) ** beta


def ethanol_table():
    return [math.log(vancoillie_sl_ext(phi, Tu, p)) for phi in PHI for Tu in TU for p in P]


def ethanol_reference_points():
    """Off-grid points of the correlation (inside its fitted range) to test the table."""
    out = []
    for phi in (0.7, 0.9, 1.1, 1.3):
        for Tu, p in ((400.0, 8e5), (550.0, 12e5), (650.0, 30e5), (800.0, 45e5)):
            out.append(dict(phi=phi, Tu=Tu, p=p, SL=vancoillie_sl(phi, Tu, p)))
    return out


def fmt(v):
    return f"{v:.6g}"


def rewrite_ts(tables, sources, ea, tad):
    src = open(TS, encoding="utf-8").read()
    L = ["// <generated:flame-speed-tables> — written by tools/reference/flamespeed_fit.py; do not edit.",
         "/** Equivalence-ratio nodes of the tables. */",
         f"export const TABLE_PHI: readonly number[] = [{', '.join(fmt(v) for v in PHI)}];",
         "/** Unburned-temperature nodes of the tables, K. */",
         f"export const TABLE_TU: readonly number[] = [{', '.join(fmt(v) for v in TU)}];",
         "/** Pressure nodes of the tables, Pa. */",
         f"export const TABLE_P: readonly number[] = [{', '.join(fmt(v) for v in P)}];",
         "/** ln(S_L / (m/s)) tables, index ((iφ·NT + iT)·NP + ip) (sources: TABLE_SOURCES, see header). */",
         "const LN_SL_TABLES: Readonly<Record<'CH4' | 'IC8H18' | 'NC7H16' | 'C2H5OH', readonly number[]>> = {"]
    for fu in ("CH4", "IC8H18", "NC7H16", "C2H5OH"):
        vals = tables.get(fu)
        if vals is None:
            L.append(f"  {fu}: [],")
            continue
        L.append(f"  {fu}: [")
        for i in range(0, len(vals), 8):
            L.append("    " + ", ".join(f"{v:.5f}" for v in vals[i:i + 8]) + ",")
        L.append("  ],")
    L.append("};")
    L.append("/** Oracle of each S_L table. */")
    L.append("export const TABLE_SOURCES: Readonly<Record<string, string>> = {")
    for fu, sname in sources.items():
        L.append(f"  {fu}: '{sname}',")
    L.append("};")
    L.append(f"/** Diluent mass fractions of the adiabatic-temperature tables. */")
    L.append(f"export const TABLE_XDIL: readonly number[] = [{', '.join(fmt(v) for v in X_DIL)}];")
    L.append("/**")
    L.append(" * ln(T_ad / K) of the fresh charge (dry air + fuel + complete-combustion products of the same")
    L.append(" * mixture at diluent mass fraction TABLE_XDIL), HP chemical equilibrium over our species set")
    L.append(" * with our NASA thermo (Cantera 3.2). Per fuel (FUEL_SPECIES order); per node 3 values.")
    L.append(" */")
    L.append("const LN_TAD_TABLES: readonly (readonly number[])[] = [")
    for fu in FUELS:
        flat = [v for node in tad[fu] for v in node]
        L.append(f"  // {fu}")
        L.append("  [")
        for i in range(0, len(flat), 9):
            L.append("    " + ", ".join(f"{v:.6f}" for v in flat[i:i + 9]) + ",")
        L.append("  ],")
    L.append("];")
    L.append("// </generated:flame-speed-tables>")
    block = "\n".join(L)
    src = re.sub(r"// <generated:flame-speed-tables>.*?// </generated:flame-speed-tables>", lambda m: block, src,
                 flags=re.S)
    rows = []
    for fu in FUELS:
        note = (" — Vancoillie et al. 2012 residual correction F (Eqs. 5-6); 1-bar row = 5-bar row"
                if fu == "C2H5OH" else " — Cantera diluted flames, flamespeed_dilution.json")
        rows.append(f"  // {fu}{note}")
        rows.append("  [")
        for a, x in enumerate(EA_X):
            for b, phi in enumerate(EA_PHI):
                vals = ", ".join(f"{v:.5g}" for v in ea[fu][a][b])
                rows.append(f"    {vals}, // x = {x}, φ = {phi}")
        rows.append("  ],")
    ea_block = ("// <generated:activation-energies>\n"
                "/** Equivalence-ratio nodes of the activation-energy tables. */\n"
                f"export const EA_PHI: readonly number[] = [{', '.join(fmt(v) for v in EA_PHI)}];\n"
                "/** Pressure nodes of the activation-energy tables, Pa. */\n"
                f"export const EA_P: readonly number[] = [{', '.join(fmt(tp[1]) for tp in EA_TP)}];\n"
                "/**\n"
                " * Effective activation energy E_a (J/mol) of the thermal-theory dilution factor, per fuel\n"
                " * (FUEL_SPECIES order), index (ix·EA_PHI.length + iφ)·EA_P.length + ip with ix = 0, 1 for the\n"
                " * diluted nodes TABLE_XDIL[1], TABLE_XDIL[2] (tools/reference/flamespeed_fit.py).\n"
                " */\n"
                "export const ACTIVATION_ENERGY_TABLE: readonly (readonly number[])[] = [\n"
                + "\n".join(rows) + "\n];\n// </generated:activation-energies>")
    src = re.sub(r"// <generated:activation-energies>.*?// </generated:activation-energies>", lambda m: ea_block, src,
                 flags=re.S)
    open(TS, "w", encoding="utf-8").write(src)


def main() -> None:
    gri = load("flamespeed_cantera.json")
    prf = load("flamespeed_prf_cantera.json")
    tables, sources, filled = {}, {}, {}
    ch4 = [r for r in gri if r["fuel"] == "CH4" and r["ydil"] == 0]
    tables["CH4"], filled["CH4"] = build_table(ch4)
    sources["CH4"] = "GRI-Mech 3.0, Cantera 3.2 FreeFlame, mixture-averaged"
    for fu in ("IC8H18", "NC7H16"):
        rows = [r for r in prf if r.get("set") == "grid" and fuel_key(r) == fu and r["ydil"] == 0]
        try:
            tables[fu], filled[fu] = build_table(rows)
            sources[fu] = "PRF_HT (Jerzembeck et al. 2009), Cantera 3.2 FreeFlame, mixture-averaged"
        except RuntimeError as err:  # oracle incomplete: keep the closed-form fallback
            print(f"{fu}: no table ({err})")
    tables["C2H5OH"] = ethanol_table()
    sources["C2H5OH"] = "Vancoillie et al. 2012 correlation (Fuel 102:460, Eqs. 2-4), power law < 5 bar"
    ea, ea_diag = ea_tables()
    for fu, f in filled.items():
        print(f"{fu}: table {len(tables[fu])} entries, {len(f)} filled by per-phi fit")
    tad = tad_tables()
    out = dict(
        source="tools/reference/flamespeed_fit.py",
        phi=PHI, Tu=TU, p=P, xDil=X_DIL, tables=tables, tableSources=sources, filled=filled,
        lnTad=tad, activationEnergy=dict(phi=EA_PHI, p=[tp[1] for tp in EA_TP], x=EA_X, table=ea,
                                         diagnostics=ea_diag),
        transport=transport_reference(),
        ethanolVancoillie2012=dict(
            source="Vancoillie et al., Fuel 102 (2012) 460, Eqs. 2-4 (F = 1), evaluated off the table grid",
            points=ethanol_reference_points()),
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, separators=(",", ":"))
        fh.write("\n")
    rewrite_ts(tables, sources, ea, tad)
    print(f"wrote {os.path.relpath(OUT, ROOT)} and the generated blocks of {os.path.relpath(TS, ROOT)}")


if __name__ == "__main__":
    main()
