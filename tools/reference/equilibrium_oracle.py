"""Oracle fixtures for src/physics/equilibrium from Cantera 3.2 (+ a numpy polish).

Writes:
  test/fixtures/equilibrium_tp.json      TP equilibria on a (composition × T × p) grid
  test/fixtures/equilibrium_energy.json  HP / UV adiabatic flames, SP expansions, TV states
  test/fixtures/equilibrium_props.json   equilibrium h, u, s, M and finite-difference cp_eq,
                                         (∂lnV/∂lnT)_p, (∂lnV/∂lnp)_T at selected states

Thermo data: Cantera ideal-gas phases built from nasa_gas.yaml with EXACTLY our species
(thermo_common.build_phase), so species data are bit-identical to species-data.ts. Each
product phase contains only the species of the elements present (e.g. no CO/CO2 without C),
as the TypeScript solver does.

Method per state:
  1. Cantera equilibrate (element_potential; fall back to vcs, then gibbs — ChemEquil fails
     below ~400 K), started from a complete-combustion composition with the target elements;
  2. polish: Newton on the element potentials π and ln n (n_j = n exp(−g°_j/RT − ln(p/p°) +
     Σ_e a_ej π_e)) to machine precision against the exact element amounts b (independent
     numpy implementation; reports how far Cantera's answer moved — ≤ 1e-8 relative for
     X > 1e-12 is asserted);
  3. HP/UV/SP: Cantera's T is refined by a Newton iteration on T with polished TP (or TV)
     equilibria and a central-difference equilibrium cp (cv) until |ΔT| < 1e-9 K.

Element amounts b are stored per kg of mixture (mol/kg) in ELEMENTS order [C, H, O, N, AR];
extensive targets (H, U, S, V) are then per kg. Deterministic.

Run: .venv/bin/python tools/reference/equilibrium_oracle.py
"""
from __future__ import annotations

import math
import os
import sys

import cantera as ct
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import build_phase, read_species, write_fixture  # noqa: E402

SPECIES = read_species()
NS = len(SPECIES)
N_EQ = 12
PROD = SPECIES[:N_EQ]
IDX = {s: i for i, s in enumerate(SPECIES)}
ELEMENTS = ["C", "H", "O", "N", "AR"]
CT_EL = {"C": "C", "H": "H", "O": "O", "N": "N", "AR": "Ar"}
R = ct.gas_constant / 1000.0  # J/(mol K)
P_ATM = 101325.0
BAR = 1e5

GAS_ALL = build_phase()
assert list(GAS_ALL.species_names) == SPECIES
P_REF = GAS_ALL.reference_pressure
MW = GAS_ALL.molecular_weights / 1000.0  # kg/mol
# atoms: A[k, e] in our element order
A_ALL = np.array([[GAS_ALL.n_atoms(k, CT_EL[e]) for e in ELEMENTS] for k in range(NS)])

# Dry air: Picard et al. 2008 (CIPM-2007) Table 1 over our 4 species (= fuels.ts DRY_AIR_REFERENCE).
DRY_AIR = {"N2": 0.780848, "O2": 0.209390, "AR": 0.009332, "CO2": 0.00040}


# ---------------------------------------------------------------------------------------------
# Compositions
# ---------------------------------------------------------------------------------------------

def vec(d: dict) -> np.ndarray:
    x = np.zeros(NS)
    for k, v in d.items():
        x[IDX[k]] += v
    return x


def norm(x: np.ndarray) -> np.ndarray:
    return x / x.sum()


def air(h2o: float = 0.0) -> np.ndarray:
    a = norm(vec(DRY_AIR))
    return norm(a * (1 - h2o) + h2o * vec({"H2O": 1.0}))


def nu_o2(fuel: np.ndarray) -> float:
    b = A_ALL.T @ fuel
    return (b[0] + b[1] / 4 - b[2] / 2) / fuel.sum()


def fuel_air(fuel: dict, phi: float, airx: np.ndarray) -> np.ndarray:
    f = norm(vec(fuel))
    n_f = phi * airx[IDX["O2"]] / nu_o2(f)
    return norm(airx + n_f * f)


def complete_products(x: np.ndarray) -> np.ndarray:
    """Element-conserving major-species allocation (mol) for moles x (any species)."""
    b = A_ALL.T @ x
    return majors_from_b(b)


def majors_from_b(b: np.ndarray, floor: float = 0.0) -> np.ndarray:
    C, H, O, N, Ar = b
    n = np.zeros(NS)
    n[IDX["N2"]] = N / 2
    n[IDX["AR"]] = Ar
    if O >= 2 * C + H / 2:
        n[IDX["CO2"]] = C
        n[IDX["H2O"]] = H / 2
        n[IDX["O2"]] = (O - 2 * C - H / 2) / 2
    else:
        assert O > C, "infeasible (O <= C)"
        co = C
        o_left = O - C
        h2o = min(H / 2, o_left)
        o_left -= h2o
        co2 = min(co, o_left)
        co -= co2
        n[IDX["CO"]] = co
        n[IDX["CO2"]] = co2
        n[IDX["H2O"]] = h2o
        n[IDX["H2"]] = H / 2 - h2o
    return n


def with_egr(fa: np.ndarray, egr_mass: float) -> np.ndarray:
    """Fresh charge fa (mole fractions) diluted with its complete-combustion products to an EGR
    mass fraction egr_mass of the total."""
    prod = complete_products(fa)
    prod = prod / prod.sum()
    m_fa = fa @ MW
    m_p = prod @ MW
    n_egr = egr_mass / (1 - egr_mass) * m_fa / m_p
    return norm(fa + n_egr * prod)


PRF90_XISO = 0.8958  # illustrative PRF blend (mole fraction iso-octane); b is what is stored
PRF80_XISO = 0.7910

COMPOSITIONS = [
    ("iC8H18-air phi=1", fuel_air({"IC8H18": 1}, 1.0, air())),
    ("iC8H18-air phi=0.2", fuel_air({"IC8H18": 1}, 0.2, air())),
    ("iC8H18-air phi=0.6", fuel_air({"IC8H18": 1}, 0.6, air())),
    ("iC8H18-air phi=1.3", fuel_air({"IC8H18": 1}, 1.3, air())),
    ("iC8H18-air phi=3", fuel_air({"IC8H18": 1}, 3.0, air())),
    ("CH4-air phi=1", fuel_air({"CH4": 1}, 1.0, air())),
    ("CH4-air phi=0.5", fuel_air({"CH4": 1}, 0.5, air())),
    ("CH4-air phi=2", fuel_air({"CH4": 1}, 2.0, air())),
    ("C2H5OH-air phi=1", fuel_air({"C2H5OH": 1}, 1.0, air())),
    ("C2H5OH-air phi=1.5", fuel_air({"C2H5OH": 1}, 1.5, air())),
    ("PRF90-air phi=1", fuel_air({"IC8H18": PRF90_XISO, "NC7H16": 1 - PRF90_XISO}, 1.0, air())),
    ("PRF80-humid air(2% H2O) phi=0.9",
     fuel_air({"IC8H18": PRF80_XISO, "NC7H16": 1 - PRF80_XISO}, 0.9, air(0.02))),
    ("iC8H18-air phi=1 + 20% EGR", with_egr(fuel_air({"IC8H18": 1}, 1.0, air()), 0.20)),
    ("PRF90-humid air(3% H2O) phi=1.1 + 10% EGR",
     with_egr(fuel_air({"IC8H18": PRF90_XISO, "NC7H16": 1 - PRF90_XISO}, 1.1, air(0.03)), 0.10)),
    ("H2-air(N2/O2 only) phi=1", norm(vec({"H2": 2 * 0.21, "O2": 0.21, "N2": 0.79}))),
    ("dry air", air()),
    ("H2-O2 phi=0.8", norm(vec({"H2": 1.6, "O2": 1.0}))),
    ("CO-O2-N2 (no H) phi=0.9", norm(vec({"CO": 1.8, "O2": 1.0, "N2": 3.76}))),
]
# Note: H2 and CO are allowed as "reactants" here only to set the element amounts b.


def b_per_kg(x: np.ndarray) -> np.ndarray:
    """Element amounts, mol per kg of mixture, of the moles/fractions x."""
    return (A_ALL.T @ x) / (x @ MW)


# ---------------------------------------------------------------------------------------------
# Equilibrium machinery (per kg basis, mol/kg)
# ---------------------------------------------------------------------------------------------

_PHASES: dict[int, tuple] = {}


def phase_for(b: np.ndarray):
    """Product phase with the species of the present elements; returns (gas, prod_idx, A)."""
    mask = 0
    tot = b.sum()
    for e in range(5):
        if b[e] > 1e-20 * tot:
            mask |= 1 << e
    if mask not in _PHASES:
        idx = [k for k in range(N_EQ)
               if all(A_ALL[k, e] == 0 or (mask >> e) & 1 for e in range(5))]
        sp = [GAS_ALL.species(SPECIES[k]) for k in idx]
        gas = ct.Solution(thermo="ideal-gas", species=sp)
        _PHASES[mask] = (gas, np.array(idx), A_ALL[idx][:, [e for e in range(5) if (mask >> e) & 1]],
                         [e for e in range(5) if (mask >> e) & 1])
    return _PHASES[mask]


def g_rt(gas, T: float) -> np.ndarray:
    gas.TP = T, P_REF
    return gas.standard_gibbs_RT.copy()


def h_rt(gas, T: float) -> np.ndarray:
    gas.TP = T, P_REF
    return gas.standard_enthalpies_RT.copy()


def s_r(gas, T: float) -> np.ndarray:
    gas.TP = T, P_REF
    return gas.standard_entropies_R.copy()


def _newton_ls(fun, z0, max_it=200):
    """Damped Newton with backtracking on the 2-norm of the (scaled) residual.
    fun(z) -> (F, J). Solves J dz = −F by least squares (tolerates near-singular J)."""
    z = z0.copy()
    F, J = fun(z)
    nF = np.linalg.norm(F)
    for _ in range(max_it):
        if not np.isfinite(nF) or nF < 1e-15:
            break
        step = np.linalg.lstsq(J, -F, rcond=None)[0]
        alpha = 1.0
        while True:
            z2 = z + alpha * step
            F2, J2 = fun(z2)
            nF2 = np.linalg.norm(F2)
            if np.isfinite(nF2) and nF2 <= (1 - 1e-4 * alpha) * nF:
                break
            alpha *= 0.5
            if alpha < 1e-12:
                return z
        z, F, J, nF = z2, F2, J2, nF2
        if np.max(np.abs(alpha * step)) < 1e-15:
            break
    return z


def _pi_start(Aa, mu, x):
    """Initial element potentials: weighted least squares of μ_j/RT = Σ_e a_je π_e."""
    ok = x > 1e-280
    w = np.sqrt(np.where(ok, x, 0.0)) + 1e-3 * ok
    return np.linalg.lstsq(Aa[ok] * w[ok, None], mu[ok] * w[ok], rcond=None)[0]


def polish_tp(Aa, b, g, lnP, n0):
    """Element-potential Newton at fixed (T, p) on z = (π, ln n), n_j = n exp(−g_j − lnP + a_j·π):
    residuals (Σ_j a_ej n_j − b_e)/b_e and Σ_j n_j/n − 1. Returns (n, rel. element residual,
    per-species conditioning Σ_e |∂ln n_j/∂ln b_e|)."""
    ne = Aa.shape[1]
    nt = n0.sum()
    x = n0 / nt
    mu = g + np.log(np.where(x > 1e-280, x, 1.0)) + lnP
    z0 = np.concatenate([_pi_start(Aa, mu, x), [math.log(nt)]])

    def fun(z):
        pi, lnn = z[:ne], z[ne]
        nj = np.exp(lnn - g - lnP + Aa @ pi)
        n = math.exp(lnn)
        F = np.concatenate([(Aa.T @ nj - b) / b, [nj.sum() / n - 1.0]])
        J = np.zeros((ne + 1, ne + 1))
        J[:ne, :ne] = ((Aa.T * nj) @ Aa) / b[:, None]
        J[:ne, ne] = (Aa.T @ nj) / b
        J[ne, :ne] = (Aa.T @ nj) / n
        return F, J

    z = _newton_ls(fun, z0)
    pi, lnn = z[:ne], z[ne]
    nj = np.exp(lnn - g - lnP + Aa @ pi)
    rel_bal = np.max(np.abs(Aa.T @ nj - b) / b)
    # conditioning: unscaled Jacobian [[G, bcur], [bcur^T, 0]] [dπ; dln n] = [db; 0]
    G = np.zeros((ne + 1, ne + 1))
    G[:ne, :ne] = (Aa.T * nj) @ Aa
    G[:ne, ne] = Aa.T @ nj
    G[ne, :ne] = Aa.T @ nj
    rhs = np.zeros((ne + 1, ne))
    rhs[:ne, :ne] = np.diag(b)
    sol = np.linalg.lstsq(G, rhs, rcond=None)[0]  # columns: response to d ln b_e
    S = Aa @ sol[:ne, :] + sol[ne, :][None, :]
    cond = np.sum(np.abs(S), axis=1)
    return nj, rel_bal, cond


def polish_tv(Aa, b, g, lnC, n0):
    """Element-potential Newton at fixed (T, V): n_j = exp(lnC − g_j + a_j·π), lnC = ln(p°V/RT)."""
    x = n0 / n0.sum()
    mu = g + np.log(np.where(x > 1e-280, n0, 1.0)) - lnC
    z0 = _pi_start(Aa, mu, x)

    def fun(pi):
        nj = np.exp(lnC - g + Aa @ pi)
        return (Aa.T @ nj - b) / b, ((Aa.T * nj) @ Aa) / b[:, None]

    pi = _newton_ls(fun, z0)
    nj = np.exp(lnC - g + Aa @ pi)
    rel_bal = np.max(np.abs(Aa.T @ nj - b) / b)
    return nj, rel_bal


STATS = {"max_polish_shift": 0.0, "max_rel_bal": 0.0, "fallbacks": 0, "max_T_refine": 0.0,
         "ill_conditioned_species": 0}
LAST_ILL: list = []


def cantera_equilibrate(gas, mode, fallback=True):
    for solver in (["element_potential", "vcs", "gibbs"] if fallback else ["element_potential"]):
        try:
            gas.equilibrate(mode, solver=solver, rtol=1e-12, max_iter=20000)
            if solver != "element_potential":
                STATS["fallbacks"] += 1
            return solver
        except ct.CanteraError:
            continue
    raise RuntimeError("all Cantera solvers failed")


def eq_tp(b: np.ndarray, T: float, p: float, n_start=None):
    """Polished TP equilibrium. Returns full-length (N_EQ) moles per kg."""
    gas, idx, Aa, els = phase_for(b)
    ba = b[els]
    if n_start is None:
        maj = majors_from_b(b)[idx]
        gas.TPX = T, p, maj / maj.sum()
        cantera_equilibrate(gas, "TP")
        n_ct = gas.X / (gas.mean_molecular_weight / 1000.0)  # mol/kg
    else:
        n_ct = n_start[idx]
    g = g_rt(gas, T)
    nj, rel, cond = polish_tp(Aa, ba, g, math.log(p / P_REF), n_ct)
    xp = nj / nj.sum()
    well = cond * 1e-15 <= 1e-5  # species not set by the last bits of b
    if n_start is None:
        xs = n_ct / n_ct.sum()
        big = (xp > 1e-12) & well
        shift = np.max(np.abs(xs[big] / xp[big] - 1)) if big.any() else 0.0
        STATS["max_polish_shift"] = max(STATS["max_polish_shift"], shift)
        assert shift < 1e-5, f"Cantera vs polished differ by {shift} at T={T} p={p}"
    assert rel < 1e-12, f"polish did not converge: rel. element residual {rel} at T={T} p={p}"
    STATS["max_rel_bal"] = max(STATS["max_rel_bal"], rel)
    out = np.zeros(N_EQ)
    out[idx] = nj
    ill = [int(idx[k]) for k in range(len(idx)) if not well[k] and xp[k] > 1e-12]
    if ill:
        STATS["ill_conditioned_species"] += len(ill)
    LAST_ILL[:] = ill
    return out


def eq_tv(b: np.ndarray, T: float, v: float, n_start):
    gas, idx, Aa, els = phase_for(b)
    g = g_rt(gas, T)
    lnC = math.log(P_REF * v / (R * T))
    nj, rel = polish_tv(Aa, b[els], g, lnC, n_start[idx])
    assert rel < 1e-12, f"TV polish did not converge: {rel}"
    STATS["max_rel_bal"] = max(STATS["max_rel_bal"], rel)
    out = np.zeros(N_EQ)
    out[idx] = nj
    return out


GAS_P = ct.Solution(thermo="ideal-gas", species=[GAS_ALL.species(s) for s in PROD])


def H_of(n, T):
    return R * T * float(n @ h_rt(GAS_P, T))


def U_of(n, T):
    return R * T * float(n @ (h_rt(GAS_P, T) - 1.0))


def S_of(n, T, p):
    nt = n.sum()
    s0 = s_r(GAS_P, T)
    pos = n > 0
    return R * float(np.sum(n[pos] * (s0[pos] - np.log(n[pos] / nt) - math.log(p / P_REF))))


def refine_T(kind: str, b, target, pv, T0, n0):
    """Newton on T with polished TP/TV equilibria and central-difference slopes.
    kind: 'HP' (target h J/kg, pv = p), 'SP' (s J/kg/K, p), 'UV' (u J/kg, pv = v m³/kg)."""
    T = T0
    n = n0
    for _ in range(60):
        def f_at(TT, nstart):
            if kind == "UV":
                nn = eq_tv(b, TT, pv, nstart)
                return U_of(nn, TT) - target, nn
            nn = eq_tp(b, TT, pv, n_start=nstart)
            if kind == "HP":
                return H_of(nn, TT) - target, nn
            return S_of(nn, TT, pv) - target, nn
        f, n = f_at(T, n)
        dT = 1e-4 * T
        fp_, _ = f_at(T + dT, n)
        fm_, _ = f_at(T - dT, n)
        slope = (fp_ - fm_) / (2 * dT)
        step = -f / slope
        T += step
        if abs(step) < 1e-9:
            break
    f, n = (f_at(T, n))
    return T, n


# ---------------------------------------------------------------------------------------------
# Fixture generation
# ---------------------------------------------------------------------------------------------

def fmt(x: float) -> float:
    return float(f"{x:.12g}")


def xrow(n: np.ndarray) -> list:
    x = n / n.sum()
    return [fmt(v) for v in x]


def tp_fixture():
    T_GRID = [300, 600, 1000, 1500, 2000, 2500, 3000, 3500, 4000]
    P_GRID = [0.1, 1, 10, 50, 150]
    comps = []
    cases = []
    for ci, (label, x) in enumerate(COMPOSITIONS):
        b = b_per_kg(x)
        comps.append({"label": label, "b": [float(v) for v in b]})
        for T in T_GRID:
            plist = P_GRID + ([300] if ci in (0, 5, 13) else [])
            for pbar in plist:
                n = eq_tp(b, T, pbar * BAR)
                row = {"c": ci, "T": T, "p": pbar * BAR, "X": xrow(n)}
                if LAST_ILL:
                    row["ill"] = list(LAST_ILL)
                cases.append(row)
    write_fixture("equilibrium_tp.json", {
        "description": "TP equilibria over the 12 product species (Cantera 3.2 + element-potential "
                       "polish). b = element amounts in mol/kg, ELEMENTS order [C,H,O,N,AR]; X = mole "
                       "fractions in SPECIES[0:12] order; ill = indices of species with X > 1e-12 whose amount is set by the last bits of b (sum_e |dln n_j/dln b_e| > 1e10: exactly stoichiometric, cold) - not comparable to 1e-4. Generated by tools/reference/equilibrium_oracle.py.",
        "species": PROD,
        "compositions": comps,
        "cases": cases,
    })


def reactant_mixtures():
    fuels = {
        "iC8H18": {"IC8H18": 1},
        "CH4": {"CH4": 1},
        "C2H5OH": {"C2H5OH": 1},
        "nC7H16": {"NC7H16": 1},
        "C3H8": {"C3H8": 1},
        "PRF90": {"IC8H18": PRF90_XISO, "NC7H16": 1 - PRF90_XISO},
    }
    out = []
    for fname, f in fuels.items():
        for phi in (0.6, 0.8, 1.0, 1.2, 1.5):
            out.append((f"{fname}-air phi={phi}", fuel_air(f, phi, air())))
    out.append(("PRF90-humid air phi=1 + 15% EGR",
                with_egr(fuel_air(fuels["PRF90"], 1.0, air(0.02)), 0.15)))
    out.append(("iC8H18-air phi=2.5", fuel_air(fuels["iC8H18"], 2.5, air())))
    out.append(("CH4-air phi=0.3", fuel_air(fuels["CH4"], 0.3, air())))
    return out


def energy_fixture():
    hp, uv, sp, tv = [], [], [], []
    mixes = reactant_mixtures()
    for label, x in mixes:
        b = b_per_kg(x)
        for (T0, p0) in [(298.15, P_ATM), (700.0, 20 * BAR)]:
            GAS_ALL.TPX = T0, p0, x
            h0 = GAS_ALL.enthalpy_mass
            u0 = GAS_ALL.int_energy_mass
            v0 = GAS_ALL.volume_mass
            # --- HP ---
            gas, idx, Aa, els = phase_for(b)
            maj = majors_from_b(b)[idx]
            gas.TPX = 1500.0, p0, maj / maj.sum()
            gas.HP = h0, p0
            cantera_equilibrate(gas, "HP")
            T_ct = gas.T
            n_ct = np.zeros(N_EQ)
            n_ct[idx] = gas.X / (gas.mean_molecular_weight / 1000.0)
            T, n = refine_T("HP", b, h0, p0, T_ct, n_ct)
            STATS["max_T_refine"] = max(STATS["max_T_refine"], abs(T - T_ct))
            hp.append({"label": label, "Xr": [float(v) for v in x], "T0": T0, "p0": p0,
                       "T": T, "T_cantera": T_ct, "X": xrow(n)})
            # --- UV (only from the engine-like state and 298 K for a subset) ---
            if T0 == 700.0 or label.endswith("phi=1.0"):
                gas.TPX = 1500.0, p0, maj / maj.sum()
                gas.UV = u0, v0
                cantera_equilibrate(gas, "UV")
                T_ct = gas.T
                n_ct = np.zeros(N_EQ)
                n_ct[idx] = gas.X / (gas.mean_molecular_weight / 1000.0)
                T, n = refine_T("UV", b, u0, v0, T_ct, n_ct)
                STATS["max_T_refine"] = max(STATS["max_T_refine"], abs(T - T_ct))
                p = n.sum() * R * T / v0
                uv.append({"label": label, "Xr": [float(v) for v in x], "T0": T0, "p0": p0,
                           "T": T, "T_cantera": T_ct, "p": p, "X": xrow(n)})
    # --- SP: isentropic equilibrium expansion of burned gas (expansion-stroke states) ---
    for label, x in [mixes[2], mixes[12], mixes[-3], mixes[-2]]:
        b = b_per_kg(x)
        for (T1, p1, p2) in [(2700.0, 60 * BAR, 3 * BAR), (2400.0, 40 * BAR, 1 * BAR),
                             (3000.0, 100 * BAR, 20 * BAR)]:
            n1 = eq_tp(b, T1, p1)
            s1 = S_of(n1, T1, p1)
            gas, idx, Aa, els = phase_for(b)
            gas.TPX = T1, p1, n1[idx] / n1[idx].sum()
            gas.SP = s1, p2
            cantera_equilibrate(gas, "SP")
            T_ct = gas.T
            n_ct = np.zeros(N_EQ)
            n_ct[idx] = gas.X / (gas.mean_molecular_weight / 1000.0)
            T, n = refine_T("SP", b, s1, p2, T_ct, n_ct)
            STATS["max_T_refine"] = max(STATS["max_T_refine"], abs(T - T_ct))
            sp.append({"label": label, "b": [float(v) for v in b], "S": s1, "p": p2,
                       "T1": T1, "p1": p1, "T": T, "T_cantera": T_ct, "X": xrow(n)})
    # --- TV: fixed (T, V) equilibria ---
    for label, x in [mixes[2], mixes[7], COMPOSITIONS[14], COMPOSITIONS[17]]:
        b = b_per_kg(x)
        for (T, rho) in [(2800.0, 20.0), (2200.0, 5.0), (1800.0, 0.5), (3500.0, 50.0)]:
            v = 1.0 / rho
            gas, idx, Aa, els = phase_for(b)
            maj = majors_from_b(b)[idx]
            gas.TDX = T, rho, maj / maj.sum()
            cantera_equilibrate(gas, "TV")
            n_ct = np.zeros(N_EQ)
            n_ct[idx] = gas.X / (gas.mean_molecular_weight / 1000.0)
            n = eq_tv(b, T, v, n_ct)
            p = n.sum() * R * T / v
            tv.append({"label": label, "b": [float(val) for val in b], "T": T, "V": v, "p": p,
                       "p_cantera": gas.P, "X": xrow(n)})
    write_fixture("equilibrium_energy.json", {
        "description": "HP/UV adiabatic flames (Xr = reactant mole fractions over all SPECIES, "
                       "T0 K, p0 Pa), SP isentropic equilibrium expansions (b mol/kg, S J/kg/K, "
                       "p Pa) and TV states (b mol/kg, V m3/kg). T = Cantera value refined by "
                       "Newton on polished equilibria (T_cantera = raw Cantera). "
                       "Generated by tools/reference/equilibrium_oracle.py.",
        "species": PROD,
        "hp": hp,
        "uv": uv,
        "sp": sp,
        "tv": tv,
    })


def props_fixture():
    cases = []
    picks = [0, 3, 4, 5, 9, 12, 14, 16]
    for ci in picks:
        label, x = COMPOSITIONS[ci]
        b = b_per_kg(x)
        for T in (1100.0, 1800.0, 2400.0, 3000.0, 3600.0):  # not 1000 K: FD across the NASA-7 Tmid break
            for pbar in (1.0, 40.0):
                p = pbar * BAR
                n = eq_tp(b, T, p)
                dT = 1e-4 * T
                npl = eq_tp(b, T + dT, p, n_start=n)
                nmi = eq_tp(b, T - dT, p, n_start=n)
                cp = (H_of(npl, T + dT) - H_of(nmi, T - dT)) / (2 * dT)
                dlnvdlnt = (math.log(npl.sum() * (T + dT)) - math.log(nmi.sum() * (T - dT))) / (
                    math.log((T + dT) / (T - dT)))
                ep = 1e-4
                npp = eq_tp(b, T, p * (1 + ep), n_start=n)
                npm = eq_tp(b, T, p * (1 - ep), n_start=n)
                dlnvdlnp = (math.log(npp.sum() / (p * (1 + ep))) - math.log(npm.sum() / (p * (1 - ep)))) / (
                    math.log((1 + ep) / (1 - ep)))
                m = 1.0  # per kg
                M = m / n.sum()
                cases.append({
                    "c": ci, "label": label, "b": [float(v) for v in b], "T": T, "p": p,
                    "M": M, "h": H_of(n, T), "u": U_of(n, T), "s": S_of(n, T, p),
                    "cp_fd": cp, "dlnV_dlnT_fd": dlnvdlnt, "dlnV_dlnp_fd": dlnvdlnp,
                })
    write_fixture("equilibrium_props.json", {
        "description": "Equilibrium properties per kg at (T, p): M kg/mol, h, u J/kg, s J/kg/K from "
                       "the polished equilibrium; cp_fd, dlnV_dlnT_fd, dlnV_dlnp_fd by central "
                       "differences of polished equilibria (dT = 1e-4 T, dp = 1e-4 p; truncation "
                       "error ~1e-7 relative). Generated by tools/reference/equilibrium_oracle.py.",
        "cases": cases,
    })


def main() -> None:
    tp_fixture()
    energy_fixture()
    props_fixture()
    print("stats:", {k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in STATS.items()})


if __name__ == "__main__":
    main()
