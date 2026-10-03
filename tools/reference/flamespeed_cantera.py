"""Oracle fixture for src/physics/combustion/laminar-flame-speed.ts: Cantera FreeFlame
(1-D freely propagating, unstretched, adiabatic premixed flame) with GRI-Mech 3.0.

Mechanism: gri30.yaml (Smith et al., GRI-Mech 3.0, 1999) as bundled with Cantera 3.2.
GRI-3.0 is optimised for natural gas (CH4) at 0.013-10 atm, 1000-2500 K; it contains C3H8
chemistry but is NOT optimised for propane — propane points are a weaker oracle (GRI-3.0
over-predicts propane S_L by ~10 % at 1 atm; see the test tolerances).

Transport: mixture-averaged (no Soret). Oxidiser: "air" = O2 + 3.76 N2 (the usual
experimental convention of the correlations being tested).

Grid convergence: each flame is solved with refine criteria (slope, curve) = (0.06, 0.12),
then re-solved from that solution with (0.03, 0.06), (0.015, 0.03) and (0.0075, 0.015),
stopping early once S_L changes by < 0.2 % between successive levels. S_L is the last value;
`dS` is |S_last - S_previous| (an upper estimate of the remaining discretisation error).
Run time per flame is logged.

Sets:
  grid:     fuel ∈ {CH4, C3H8} × φ ∈ {0.7, 0.8, 1.0, 1.2, 1.4} × Tu ∈ {300, 500, 700} K
            × p ∈ {1, 5, 20} bar, undiluted
  dilution: fuel ∈ {CH4, C3H8}, φ = 1, (Tu, p) ∈ {(300 K, 1 bar), (500 K, 5 bar)},
            burned-gas diluent mass fraction ∈ {0.1, 0.2, 0.3}; diluent = complete-combustion
            products of the stoichiometric mixture (CO2 + H2O + N2).
  engine:   fuel ∈ {CH4, C3H8}, φ ∈ {0.8, 1.0, 1.2} at Tu = 800 K, p = 40 bar (beyond the
            GRI-3.0 validation range; used only to check extrapolation behaviour).
  CH4 ext:  φ ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6} × Tu ∈ {300, 500, 700, 900} K
            × p ∈ {1, 5, 20, 60} bar — data for the fitted methane engine correlation.

For every flame the fixture also stores the burned temperature T_b and density ratio
ρ_u/ρ_b at the domain end, the thermal thickness (T_b − T_u)/max|dT/dx|, the unburned
thermal diffusivity α_u = λ/(ρ c_p), and the mixture-averaged mass diffusivity of the
deficient reactant (fuel if φ ≤ 1 else O2) for the effective Lewis number.

Writes test/fixtures/flamespeed_cantera.json. Run: .venv/bin/python tools/reference/flamespeed_cantera.py
(≈ 30 min on 8 cores from scratch; incremental re-runs only compute new points).
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import time

import cantera as ct
import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "test", "fixtures", "flamespeed_cantera.json")

PHI = [0.7, 0.8, 1.0, 1.2, 1.4]
TU = [300.0, 500.0, 700.0]
P_BAR = [1.0, 5.0, 20.0]
FUELS = ["CH4", "C3H8"]
DIL_CASES = [(300.0, 1.0), (500.0, 5.0)]
DIL_Y = [0.1, 0.2, 0.3]
ENGINE_PHI = [0.8, 1.0, 1.2]
ENGINE_TP = [(800.0, 40.0)]
EXT_PHI = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]
EXT_TU = [300.0, 500.0, 700.0, 900.0]
EXT_P = [1.0, 5.0, 20.0, 60.0]
AIR = "O2:1.0, N2:3.76"
LEVELS = [(3.0, 0.06, 0.12, 0.03), (2.0, 0.03, 0.06, 0.015), (2.0, 0.015, 0.03, 0.0075),
          (2.0, 0.0075, 0.015, 0.004)]
CONV_TOL = 0.002  # stop refining once S_L changes by less than 0.2 % between levels
# FLAMESPEED_MAXLEVELS caps the refinement levels (default all 4); stored per flame ("levels").
MAX_LEVELS = int(os.environ.get("FLAMESPEED_MAXLEVELS", str(len(LEVELS))))


def sig(x: float, n: int = 6) -> float:
    if not math.isfinite(x):
        raise ValueError(x)
    return float(f"{x:.{n}g}")


def make_gas() -> ct.Solution:
    gas = ct.Solution("gri30.yaml")
    gas.transport_model = "mixture-averaged"
    return gas


def set_mixture(gas: ct.Solution, fuel: str, phi: float, Tu: float, p: float, ydil: float) -> None:
    gas.set_equivalence_ratio(phi, f"{fuel}:1", AIR)
    gas.TP = Tu, p
    if ydil > 0:
        y_fresh = gas.Y.copy()
        prod = make_gas()
        prod.set_equivalence_ratio(1.0, f"{fuel}:1", AIR)
        # Complete-combustion products of the stoichiometric mixture: C -> CO2, H -> H2O, N2.
        nC = prod.elemental_mole_fraction("C")
        nH = prod.elemental_mole_fraction("H")
        nN = prod.elemental_mole_fraction("N")
        prod.TPX = Tu, p, {"CO2": nC, "H2O": nH / 2, "N2": nN / 2}
        gas.TPY = Tu, p, (1 - ydil) * y_fresh + ydil * prod.Y


def solve_one(args):
    fuel, phi, Tu, p_bar, ydil = args
    p = p_bar * 1e5
    gas = make_gas()
    set_mixture(gas, fuel, phi, Tu, p, ydil)
    rho_u = gas.density
    alpha_u = gas.thermal_conductivity / (gas.density * gas.cp_mass)
    dmix = gas.mix_diff_coeffs_mass
    kdef = gas.species_index(fuel) if phi <= 1.0 else gas.species_index("O2")
    d_def = float(dmix[kdef])
    # HP equilibrium (adiabatic flame temperature) for reference
    eq = make_gas()
    eq.TPY = gas.T, gas.P, gas.Y
    eq.equilibrate("HP")
    T_ad, rho_ad = eq.T, eq.density

    width = max(0.003, 0.03 / p_bar ** 0.7)
    t0 = time.time()
    f = ct.FreeFlame(gas, width=width)
    f.transport_model = "mixture-averaged"
    f.set_max_grid_points(f.flame, 4000)
    speeds, npts = [], []
    for i, (ratio, slope, curve, prune) in enumerate(LEVELS[:MAX_LEVELS]):
        if i >= 2 and abs(speeds[-1] - speeds[-2]) <= CONV_TOL * speeds[-1]:
            break
        f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
        f.solve(loglevel=0, auto=(i == 0))
        speeds.append(float(f.velocity[0]))
        npts.append(len(f.grid))
    dt = time.time() - t0
    T = f.T
    x = f.grid
    dTdx = np.max(np.abs(np.diff(T) / np.diff(x)))
    Tb = float(T[-1])
    return dict(
        fuel=fuel, phi=phi, Tu=Tu, p=p, ydil=ydil,
        SL=sig(speeds[-1]), dS=sig(abs(speeds[-1] - speeds[-2]), 3),
        Tb=sig(Tb), sigma=sig(rho_u / float(f.density[-1])),
        Tad=sig(T_ad), sigmaEq=sig(rho_u / rho_ad),
        deltaT=sig((Tb - Tu) / dTdx), alphaU=sig(alpha_u), Ddef=sig(d_def),
        points=npts[-1], seconds=round(dt, 1), levels=len(speeds),
    )


def all_cases():
    cases = [(f, phi, Tu, p, 0.0) for f in FUELS for p in P_BAR for Tu in TU for phi in PHI]
    cases += [(f, 1.0, Tu, p, y) for f in FUELS for (Tu, p) in DIL_CASES for y in DIL_Y]
    cases += [(f, phi, Tu, p, 0.0) for f in FUELS for (Tu, p) in ENGINE_TP for phi in ENGINE_PHI]
    # methane extension for the fitted engine correlation (GRI-3.0 is a good CH4 oracle):
    cases += [("CH4", phi, Tu, p, 0.0) for phi in EXT_PHI for Tu in EXT_TU for p in EXT_P
              if not (phi == 0.5 and Tu == 300.0 and p >= 20.0)]  # sub-limit; filled by the fit
    seen, out = set(), []
    for c in cases:
        k = (c[0], round(c[1], 4), c[2], c[3], c[4])
        if k not in seen:
            seen.add(k)
            out.append(c)
    return out


def main() -> None:
    """Incremental: flames already in the fixture (same fuel, φ, Tu, p, y_dil) are reused, so
    extending a grid only computes the new points. Set FLAMESPEED_FULL=1 to recompute all."""
    cases = all_cases()
    old = {}
    if os.path.exists(OUT) and not os.environ.get("FLAMESPEED_FULL"):
        with open(OUT, encoding="utf-8") as fh:
            for r in json.load(fh)["flames"]:
                old[(r["fuel"], round(r["phi"], 4), r["Tu"], r["p"] / 1e5, r["ydil"])] = r
    todo = [c for c in cases if (c[0], round(c[1], 4), c[2], c[3], c[4]) not in old]
    results = [old[k] for k in ((c[0], round(c[1], 4), c[2], c[3], c[4]) for c in cases) if k in old]
    print(f"{len(results)} flames reused, {len(todo)} to compute", flush=True)
    todo.sort(key=lambda c: (-c[3], -abs(c[1] - 1.05)))  # slow (high-p) first for load balance
    t0 = time.time()
    nproc = int(os.environ.get("FLAMESPEED_NPROC", max(1, min(8, (os.cpu_count() or 2) - 2))))
    failed = []
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(solve_one_safe, todo):
            if "error" in r:
                failed.append(r)
                print(f"FAILED {r}", flush=True)
                continue
            results.append(r)
            print(f"{r['fuel']:5s} phi={r['phi']:.1f} Tu={r['Tu']:.0f} p={r['p']/1e5:4.0f} bar "
                  f"ydil={r['ydil']:.1f}: SL={r['SL']:.4f} m/s dS={r['dS']:.1e} "
                  f"({r['points']} pts, {r['seconds']} s)", flush=True)
            if len(results) % 4 == 0:
                write_fixture(results, failed)  # checkpoint (fixture always valid, incremental)
    write_fixture(results, failed)
    print(f"wrote {os.path.relpath(OUT, ROOT)} ({len(results)} flames, {len(failed)} failed, "
          f"{time.time() - t0:.0f} s wall)")


def write_fixture(results, failed):
    key = lambda r: (r["fuel"], r["ydil"], r["p"], r["Tu"], r["phi"])  # noqa: E731
    data = dict(
        source="Cantera %s FreeFlame, gri30.yaml (GRI-Mech 3.0), mixture-averaged transport, "
               "air = O2 + 3.76 N2; see tools/reference/flamespeed_cantera.py" % ct.__version__,
        units=dict(SL="m/s", Tu="K", p="Pa", Tb="K", deltaT="m", alphaU="m^2/s", Ddef="m^2/s"),
        cpuSeconds=round(sum(r.get("seconds", 0) for r in results), 1),
        failed=[{k: r[k] for k in ("fuel", "phi", "Tu", "p", "ydil", "error")} for r in failed],
        flames=sorted(results, key=key),
    )
    tmp = OUT + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, separators=(",", ":"), allow_nan=False)
        fh.write("\n")
    os.replace(tmp, OUT)


def solve_one_safe(case):
    try:
        return solve_one(case)
    except Exception as err:  # noqa: BLE001 — record and drop flames that do not converge
        fuel, phi, Tu, p_bar, ydil = case
        return dict(fuel=fuel, phi=phi, Tu=Tu, p=p_bar * 1e5, ydil=ydil, error=str(err)[:200])


if __name__ == "__main__":
    sys.exit(main())
