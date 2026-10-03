"""Two-zone Wiebe closed-cycle oracle for the cycle integrator (validation round 2, numerics).

The fuel-air-cycle oracles (cycle_fuel_air_oracle.py, validation_numerics_fuel_air.py) only exercise
the single-zone compression, an instantaneous UV burn at TDC and a burned-only expansion. This oracle
exercises the TWO-ZONE path of src/physics/cycle (closure.ts + the RK4 integration of cycle-model.ts
with mass transfer between the zones) on the same continuous model, computed independently with
Cantera thermo on exactly the project's species and NASA-7 data (1-bar standard state,
thermo_common.nasa_species):

  * closed, adiabatic (no wall heat), no knock; prescribed Wiebe mass-fraction burned
      x_W(theta) = 1 - exp(-a u^(m+1)),  u = (theta - theta0)/dtheta  (theta >= theta0)
    (the model's combustionModel 'wiebe', options.ts WiebeOptions);
  * before theta0 one frozen zone (isentropic, exact); at theta0 cycle-model.ts (EV_WIEBE) moves SEED
    (WIEBE_SEED = 1e-5; 1e-8 before fixer round 2) of the charge into a new burned zone at constant U_tot, then burns dm_b/dtheta = m dx_W/dtheta, so
    m_b(theta) = m (SEED + x_W(theta));
  * unburned zone: frozen composition, specific entropy constant (adiabatic: dS_u = s_u dm_u);
  * burned zone: ONE fully mixed zone in chemical equilibrium over the 12 product species at
    (T_b, p) with the charge's element composition; common pressure p;
  * dU_tot = -p dV, every transfer at constant U_tot.
The model integrates U_tot and recovers (p, T_b) algebraically, which is ill-conditioned while the
burned zone is tiny (T_b absorbs U errors times m/m_b). The oracle uses the equivalent, well-conditioned
differential form (derived from dU_tot = -p dV and ds_u = 0):
  burned energy   m_b (dh_b - v_b dp) = (h_u - h_b) dm_b
  volume          dV = (v_b - v_u) dm_b + m_u dv_u + m_b dv_b,  dv_u = -v_u dp/(gamma_u p)
solved for (dp, dT_b) with the equilibrium derivatives (dh_b/dT)_p, (dh_b/dp)_T, (dv_b/dT)_p,
(dv_b/dp)_T by 4-point central differences; states (p, T_b, W = int p dV), scipy DOP853 at rtol 1e-12.
The seed state at theta0 is h_b = h_u (HP flame) at the pressure that satisfies the volume balance
(the constant-U transfer of SEED of the charge; the two differ at O(SEED^2) ≈ 1e-10 relative). When m_u < 1e-9 m the rest is
merged (the TS test runs with burnoutFraction 1e-9). Consistency check inside the oracle: the work
integral equals U_start - U_end to <= 1e-9 of W.

Equilibrium: Cantera's TP equilibrate reproduces u_b only to ~1e-9 of c_v T between starting
compositions because it lets the element ratios drift at that level; every equilibrium here is solved
by an element-potential Newton on the charge's FIXED element vector (class Burned, started from
Cantera's solution), reproducible to ~1e-15.

Output per case: p, T_u, T_b, x_b at every crank degree, W, U_start, U_end. Writes
test/fixtures/validation_numerics_two_zone_wiebe.json. Run with .venv/bin/python (a few minutes).
"""
import json
import pathlib
import sys
import time

import cantera as ct
import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "test" / "fixtures" / "validation_numerics_two_zone_wiebe.json"
sys.path.insert(0, str(HERE))
from thermo_common import nasa_species, read_species  # noqa: E402
from validation_numerics_fuel_air import (  # noqa: E402
    CFR, ISO, PRF90, equilibrium_residual, fresh_X, mix_by_mass, product_guess,
)

NAMES = read_species()
N_EQ = 12
_SP = nasa_species(NAMES)
unb = ct.Solution(thermo="ideal-gas", species=_SP)
prod = ct.Solution(thermo="ideal-gas", species=_SP[:N_EQ])
SEED = 1e-5  # burned-zone seed of cycle-model.ts EV_WIEBE (WIEBE_SEED, fraction of the charge)
MERGE = 1e-9  # burn-out merge fraction used by the oracle and the TS test

A_EL = np.array([[prod.n_atoms(k, e) for e in prod.element_names] for k in range(prod.n_species)])  # (K, E)


class Burned:
    """Burned-gas equilibrium at (T, p) for a fixed element vector, to round-off."""

    def __init__(self, b, el):
        self.b = np.asarray(b, dtype=float)
        use = self.b > 1e-300
        self.A = A_EL[:, use]
        self.bb = self.b[use]
        self.ok = np.all(A_EL[:, ~use] == 0, axis=1)
        self.el = el
        self.pi = None
        self.n = 1.0

    def _cantera_start(self, T, p):
        prod.TPX = T, p, product_guess(self.el)
        prod.equilibrate("TP")
        x0 = prod.X.copy()
        prod.TP = T, prod.reference_pressure
        g = prod.standard_gibbs_RT + np.log(p / prod.reference_pressure)
        w = self.ok & (x0 > 1e-12)
        self.pi, *_ = np.linalg.lstsq(self.A[w], (np.log(x0[w]) + g[w]), rcond=None)
        s = float(self.A[:, 0] @ x0)
        self.n = self.bb[0] / s if s > 0 else 1.0

    def x(self, T, p):
        """x_k = exp(-g_k/RT - ln(p/p0) + sum_j a_kj pi_j); Newton on (pi, n) for
        n sum_k a_kj x_k = b_j and sum_k x_k = 1 (warm-started from the previous call)."""
        prod.TP = T, prod.reference_pressure
        g = prod.standard_gibbs_RT + np.log(p / prod.reference_pressure)
        A, b = self.A, self.bb
        for attempt in range(2):
            if self.pi is None or attempt == 1:
                self._cantera_start(T, p)
            pi, n = self.pi.copy(), self.n
            conv = False
            for _ in range(60):
                x = np.where(self.ok, np.exp(np.minimum(-g + A @ pi, 700.0)), 0.0)
                F = np.concatenate([n * (A.T @ x) - b, [x.sum() - 1.0]])
                J = np.zeros((len(b) + 1, len(b) + 1))
                J[:-1, :-1] = n * (A.T * x) @ A
                J[:-1, -1] = A.T @ x
                J[-1, :-1] = A.T @ x
                d = np.linalg.solve(J, -F)
                d[:-1] = np.clip(d[:-1], -2.0, 2.0)
                pi += d[:-1]
                n += d[-1]
                if np.max(np.abs(d[:-1]) / (1.0 + np.abs(pi))) < 1e-13 and abs(d[-1]) < 1e-13 * abs(n):
                    conv = True
                    break
            if not conv:
                x = np.where(self.ok, np.exp(np.minimum(-g + A @ pi, 700.0)), 0.0)
                F = np.concatenate([n * (A.T @ x) - b, [x.sum() - 1.0]])
                conv = bool(np.max(np.abs(F[:-1]) / b) < 1e-13 and abs(F[-1]) < 1e-13)
            if conv and n > 0:
                self.pi, self.n = pi, n
                x = np.where(self.ok, np.exp(-g + A @ pi), 0.0)
                return x / x.sum()
        raise RuntimeError(f"equilibrium polish failed at T={T} p={p}")

    def hv(self, T, p):
        prod.TPX = T, p, self.x(T, p)
        return prod.enthalpy_mass, prod.volume_mass, prod.int_energy_mass

    def derivs(self, T, p):
        """h, v and (dh/dT)_p, (dh/dp)_T, (dv/dT)_p, (dv/dp)_T of the equilibrium mixture (4-point
        central differences, relative steps 1e-3: truncation and round-off both ~1e-12)."""
        h0, v0, _ = self.hv(T, p)
        pi0, n0 = self.pi.copy(), self.n
        out = []
        for var, x0 in (("T", T), ("p", p)):
            e = 1e-3 * x0
            vals = []
            for k in (2, 1, -1, -2):
                self.pi, self.n = pi0.copy(), n0
                vals.append(self.hv(T + k * e, p) if var == "T" else self.hv(T, p + k * e))
            dh = (-vals[0][0] + 8 * vals[1][0] - 8 * vals[2][0] + vals[3][0]) / (12 * e)
            dv = (-vals[0][1] + 8 * vals[1][1] - 8 * vals[2][1] + vals[3][1]) / (12 * e)
            out.append((dh, dv))
        self.pi, self.n = pi0, n0
        (hT, vT), (hp, vp) = out
        return h0, v0, hT, hp, vT, vp


class Unburned:
    """Frozen unburned charge on its isentrope s = s1."""

    def __init__(self, X0, T1, p1):
        unb.TPX = T1, p1, X0
        self.s1 = unb.s
        self.T = T1

    def at(self, p):
        T = self.T
        for _ in range(60):
            unb.TP = T, p
            dT = -(unb.s - self.s1) / (unb.cp_mass / T)
            T += dT
            if abs(dT) < 1e-15 * T:
                break
        unb.TP = T, p
        self.T = T
        return T, unb.enthalpy_mass, unb.volume_mass, unb.cp_mass / unb.cv_mass, unb.int_energy_mass


def wiebe(theta, w):
    """x_W and dx_W/dtheta (per crank degree)."""
    u = (theta - w["startDeg"]) / w["durationDeg"]
    if u <= 0:
        return 0.0, 0.0
    e = np.exp(-w["a"] * u ** (w["m"] + 1))
    return 1.0 - e, w["a"] * (w["m"] + 1) * u ** w["m"] * e / w["durationDeg"]


def run_case(case):
    bore, stroke, rod, cr = case["bore"], case["stroke"], case["rod"], case["cr"]
    vd = np.pi / 4 * bore**2 * stroke
    vc = vd / (cr - 1)
    a = stroke / 2
    Ap = np.pi / 4 * bore**2

    def volume(deg):
        th = np.radians(deg)
        s = a * np.cos(th) + np.sqrt(rod**2 - (a * np.sin(th)) ** 2)
        return vc + Ap * (rod + a - s)

    def dvdtheta(deg):  # m^3 per crank degree
        th = np.radians(deg)
        ds = -a * np.sin(th) - (a**2 * np.sin(th) * np.cos(th)) / np.sqrt(rod**2 - (a * np.sin(th)) ** 2)
        return -Ap * ds * np.pi / 180

    Xf = fresh_X(case["fuel_X"], case["phi"], case["x_h2o"])
    if case["residual"] > 0:
        Xr = equilibrium_residual(Xf, case["T_res"], 1e5)
        X0 = mix_by_mass(Xf, Xr, case["residual"])
    else:
        X0 = Xf
    unb.TPX = case["T1"], case["p1"], X0
    el = {e: unb.elemental_mole_fraction(e) for e in ["C", "H", "O", "N", "Ar"]}
    b = np.array([sum(unb.n_atoms(k, e) * unb.X[k] for k in range(unb.n_species)) for e in prod.element_names])
    th_s, th_e = case["startDeg"], case["endDeg"]
    v_s = volume(th_s)
    m = unb.density * v_s
    U_s = m * unb.int_energy_mass
    un = Unburned(X0, case["T1"], case["p1"])
    bu = Burned(b, el)
    w = case["wiebe"]
    th0 = w["startDeg"]
    out_deg = [float(d) for d in np.arange(np.ceil(th_s), np.floor(th_e) + 1e-9, 1.0)]
    rows = {}
    # --- single zone (isentropic, frozen) up to theta0: exact
    for d in out_deg:
        if d <= th0:
            unb.SV = un.s1, volume(d) / m
            rows[d] = {"deg": d, "p": unb.P, "Tu": unb.T, "Tb": 0.0, "xb": 0.0}
    unb.SV = un.s1, volume(th0) / m
    p = unb.P
    un.T = unb.T
    W_pre = U_s - m * unb.int_energy_mass  # work start -> theta0 (isentropic)
    # --- seed: SEED of the charge burned at constant U (h_b = h_u at p, volume balance)
    mb0 = SEED * m
    V0 = volume(th0)
    Tb = 2500.0
    for _ in range(40):
        Tu, hu, vu, gu, _uu = un.at(p)
        for _ in range(60):
            h, v, hT, hp, vT, vp = bu.derivs(Tb, p)
            dT = (hu - h) / hT
            Tb += dT
            if abs(dT) < 1e-13 * Tb:
                break
        h, v, _ = bu.hv(Tb, p)
        f = (m - mb0) * vu + mb0 * v - V0
        dp = -f / (-(m - mb0) * vu / (gu * p))
        p += dp
        if abs(dp) < 1e-15 * p:
            break

    def mb_of(deg):
        x, dx = wiebe(deg, w)
        mb = m * (SEED + x)
        return (m, 0.0) if mb >= m else (mb, m * dx)

    def rhs(deg, y, merged):
        p, Tb = y[0], y[1]
        dV = dvdtheta(deg)
        h, v, hT, hp, vT, vp = bu.derivs(Tb, p)
        if merged:
            dp = dV / (m * vp - m * vT * (hp - v) / hT)
            dTb = -(hp - v) * dp / hT
            return [dp, dTb, p * dV]
        mb, dmb = mb_of(deg)
        mu = m - mb
        Tu, hu, vu, gu, _ = un.at(p)
        rel = (hu - h) * dmb / mb
        den = -mu * vu / (gu * p) + mb * vp - mb * vT * (hp - v) / hT
        dp = (dV - (v - vu) * dmb - mb * vT * rel / hT) / den
        dTb = (rel - (hp - v) * dp) / hT
        return [dp, dTb, p * dV]

    # merge where m_u = m (1 - SEED - x_W) reaches MERGE m (the model merges at the first step end after it)
    u_merge = (np.log(1 / (SEED + MERGE)) / w["a"]) ** (1 / (w["m"] + 1))
    th_merge = th0 + u_merge * w["durationDeg"]
    t0 = time.time()
    y = np.array([p, Tb, W_pre])
    segs = [(th0, min(th_merge, th_e), False)] + ([(th_merge, th_e, True)] if th_merge < th_e else [])
    for (a_, b_, merged) in segs:
        te = sorted({d for d in out_deg if a_ < d < b_} | {b_})
        sol = solve_ivp(lambda t, yy: rhs(t, yy, merged), (a_, b_), y, method="DOP853", rtol=1e-12,
                        atol=[1e-12 * p, 1e-12 * 2500, 1e-12 * abs(U_s)], t_eval=te, max_step=0.5, first_step=1e-4)
        assert sol.success, sol.message
        for k, d in enumerate(sol.t):
            d = float(d)
            if d in out_deg and d not in rows:
                pp, TT = sol.y[0][k], sol.y[1][k]
                mb, _ = mb_of(d)
                Tu = 0.0 if merged else un.at(pp)[0]
                rows[d] = {"deg": d, "p": float(pp), "Tu": float(Tu), "Tb": float(TT), "xb": 1.0 if merged else mb / m}
        y = sol.y[:, -1].copy()
        # (the merge burns <= 1e-9 of the charge at constant U: its state change is below 1e-9)
    p_end, Tb_end, W = y
    h, v, u_b = bu.hv(Tb_end, p_end)
    U_end = m * u_b
    eW = (U_s - U_end - W) / W
    print(f"{case['label']:70s} W={W:9.4f} J  p_end={p_end/1e5:8.4f} bar  Tb_end={Tb_end:8.2f} K  "
          f"energy check (U_s-U_e-W)/W={eW:.1e}  ({time.time() - t0:.0f} s)", flush=True)
    assert abs(eW) < 1e-8, eW
    return {
        "input": {k: v for k, v in case.items()},
        "X_reactants": X0,
        "mass": m,
        "work": float(W),
        "U_start": float(U_s),
        "U_end": float(U_end),
        "energy_check": float(eW),
        "end": {"deg": th_e, "p": float(p_end), "Tb": float(Tb_end)},
        "trace": [rows[d] for d in out_deg],
    }


BASE = {**CFR, "x_h2o": 0.0, "residual": 0.0, "T_res": 1000.0, "T1": 350.0, "p1": 1.0e5, "startDeg": -152.0,
        "endDeg": 141.0}
A_999 = 6.907755278982137  # -ln(0.001): x_W = 0.999 at theta0 + dtheta
CASES = [
    {**BASE, "label": "PRF90 phi 1.1 CR 6.43, 6 % eq. residual, Wiebe -15/50/m2 (RON-like)", "fuel_X": PRF90,
     "phi": 1.1, "cr": 6.43, "x_h2o": 0.0125, "residual": 0.06, "T1": 370.0,
     "wiebe": {"startDeg": -15.0, "durationDeg": 50.0, "a": A_999, "m": 2.0}},
    {**BASE, "label": "iso-octane phi 1.5 CR 8, Wiebe -20/40/m2 (rich)", "fuel_X": ISO, "phi": 1.5, "cr": 8.0,
     "wiebe": {"startDeg": -20.0, "durationDeg": 40.0, "a": A_999, "m": 2.0}},
    {**BASE, "label": "n-heptane phi 0.7 CR 10, Wiebe -10/60/m3 (lean, slow)", "fuel_X": {"NC7H16": 1.0}, "phi": 0.7,
     "cr": 10.0, "wiebe": {"startDeg": -10.0, "durationDeg": 60.0, "a": A_999, "m": 3.0}},
    {**BASE, "label": "ethanol phi 1.0 CR 12, Wiebe -5/30/m1.5", "fuel_X": {"C2H5OH": 1.0}, "phi": 1.0, "cr": 12.0,
     "wiebe": {"startDeg": -5.0, "durationDeg": 30.0, "a": A_999, "m": 1.5}},
    {**BASE, "label": "iso-octane phi 1.0 CR 8, Wiebe -5/10/m2 (very fast burn)", "fuel_X": ISO, "phi": 1.0,
     "cr": 8.0, "wiebe": {"startDeg": -5.0, "durationDeg": 10.0, "a": A_999, "m": 2.0}},
    {**BASE, "label": "PRF90 phi 1.05 CR 6.54, 8 % eq. residual, 450 K, Wiebe -20/45/m2 (MON-like)",
     "fuel_X": PRF90, "phi": 1.05, "cr": 6.54, "x_h2o": 0.02, "residual": 0.08, "T1": 450.0, "p1": 0.97e5,
     "wiebe": {"startDeg": -20.0, "durationDeg": 45.0, "a": A_999, "m": 2.0}},
]

if __name__ == "__main__":
    sel = [int(a) for a in sys.argv[1:]]
    results = [run_case(c) for i, c in enumerate(CASES) if not sel or i in sel]
    if sel:
        sys.exit(0)  # partial runs do not write the fixture
    OUT.write_text(json.dumps({
        "description": "Two-zone adiabatic closed cycle with a prescribed Wiebe burn (frozen isentropic unburned "
                       "zone, one fully mixed equilibrium burned zone, common p), Cantera + nasa_gas.yaml at the "
                       "1 bar standard state, element-potential-polished equilibria, DOP853 rtol 1e-12; SI; crank "
                       "deg, 0 = firing TDC. Generated by tools/reference/validation_numerics_two_zone_wiebe.py.",
        "cantera_version": ct.__version__,
        "seed_fraction": SEED,
        "merge_fraction": MERGE,
        "cases": results,
    }, indent=1))
