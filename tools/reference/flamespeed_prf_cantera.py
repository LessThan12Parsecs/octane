"""Oracle fixture for the iso-octane / n-heptane / PRF laminar burning velocity in
src/physics/combustion/laminar-flame-speed.ts: Cantera FreeFlame (1-D freely propagating,
unstretched, adiabatic) with the PRF_HT reduced mechanism of Jerzembeck, Peters,
Pepiot-Desjardins & Pitsch, Combust. Flame 156 (2009) 292-301 (99 species / 601 reactions,
reduced from the LLNL n-heptane and iso-octane mechanisms; validated by the authors against
PRF burning velocities at 1 atm and their own 10-25 bar, 373 K bomb data; loaded from
test/fixtures/flamespeed_mech_prf_ht.json, see flamespeed_prf_mechanism.py).

It is a HIGH-TEMPERATURE mechanism (no low-temperature/NTC chemistry), so the fresh mixture
does not autoignite upstream of the flame even at Tu = 900 K — i.e. these are laminar
burning velocities in the classical sense, which is what flame-propagation models need.

Transport: mixture-averaged. Oxidiser: O2 + 3.76 N2. Grid convergence: refine criteria
(slope, curve) = (0.06, 0.12) → (0.03, 0.06) → (0.015, 0.03), stopping early if S_L changes
by < 0.3 %; `dS` = |S_last − S_previous|.

Sets:
  grid:     fuel ∈ {IC8H18, NC7H16} × φ ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6}
            × Tu ∈ {300, 500, 700, 900} K × p ∈ {1, 5, 20, 60} bar (except φ = 0.5 at 300 K and
            ≥ 20 bar: below the lean limit, left to the per-φ fill of flamespeed_fit.py)
  blend:    iso-octane mole fraction in the fuel x ∈ {0.5, 0.9}, φ ∈ {0.8, 1.0, 1.2},
            (Tu, p) ∈ {(500 K, 5 bar), (700 K, 20 bar)}  — for the PRF mixing rule
  dilution: IC8H18 and NC7H16, φ = 1, (Tu, p) ∈ {(500 K, 5 bar), (700 K, 20 bar)}, burned-gas
            diluent mass fraction ∈ {0.1, 0.2, 0.3} (complete-combustion products CO2+H2O+N2
            of the stoichiometric mixture of the same fuel)
  propane:  C3H8, φ ∈ {0.8, 1.0, 1.2}, (Tu, p) ∈ {(300 K, 1 bar), (500 K, 5 bar), (700 K, 20 bar)}
            — cross-check of the GRI-3.0 propane oracle
  jerzembeck: IC8H18, NC7H16 at 373 K, 10 and 25 bar, φ ∈ {0.7, 0.8, 1.0, 1.2} — off-grid
            interpolation checks at the conditions of the Jerzembeck et al. experiments
  jerzembeck-o2: the same flames with the experiments' oxidiser (X_O2 = 0.205, rest N2): the
            ratio to "jerzembeck" converts model values to the experimental oxidiser (−5…−10 %)

Writes test/fixtures/flamespeed_prf_cantera.json (incremental: existing flames are reused).
Run (≈ 1.5-3 h on 8 cores from scratch):
  .venv/bin/python tools/reference/flamespeed_prf_cantera.py
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os
import sys
import time
import warnings

import cantera as ct
import numpy as np

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
MECH = os.path.join(ROOT, "test", "fixtures", "flamespeed_mech_prf_ht.json")
OUT = os.path.join(ROOT, "test", "fixtures", "flamespeed_prf_cantera.json")

AIR = "O2:1.0, N2:3.76"
NAME = {"IC8H18": "IXC8H18", "NC7H16": "NXC7H16", "C3H8": "C3H8"}
LEVELS = [(3.0, 0.06, 0.12, 0.03), (2.0, 0.03, 0.06, 0.015), (2.0, 0.015, 0.03, 0.0075)]
# FLAMESPEED_MAXLEVELS=2 trades ~0.5 % accuracy (level 2 → 3 changes S_L by 0.3-1 %) for ~2x speed;
# the number of levels actually used is stored per flame ("levels").
MAX_LEVELS = int(os.environ.get("FLAMESPEED_MAXLEVELS", "3"))
CONV_TOL = 0.003

PHI = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.6]
TU = [300.0, 500.0, 700.0, 900.0]
P_BAR = [1.0, 5.0, 20.0, 60.0]

_yaml = None


def sig(x: float, n: int = 6) -> float:
    if not math.isfinite(x):
        raise ValueError(x)
    return float(f"{x:.{n}g}")


def make_gas() -> ct.Solution:
    global _yaml
    if _yaml is None:
        with open(MECH, encoding="utf-8") as fh:
            _yaml = json.load(fh)["yaml"]
    warnings.filterwarnings("ignore")
    gas = ct.Solution(yaml=_yaml, name="gas")
    gas.transport_model = "mixture-averaged"
    return gas


def fuel_string(fuel: dict) -> str:
    return ", ".join(f"{NAME[k]}:{v}" for k, v in fuel.items())


# Oxidiser of the Jerzembeck et al. (2009) bomb: "X_O2^air = 0.205" (figure captions).
AIR_JERZEMBECK = "O2:0.205, N2:0.795"


def set_mixture(gas, fuel: dict, phi, Tu, p, ydil, air: str = AIR) -> None:
    gas.set_equivalence_ratio(phi, fuel_string(fuel), air)
    gas.TP = Tu, p
    if ydil > 0:
        y_fresh = gas.Y.copy()
        prod = make_gas()
        prod.set_equivalence_ratio(1.0, fuel_string(fuel), AIR)
        nC = prod.elemental_mole_fraction("C")
        nH = prod.elemental_mole_fraction("H")
        nN = prod.elemental_mole_fraction("N")
        prod.TPX = Tu, p, {"CO2": nC, "H2O": nH / 2, "N2": nN / 2}
        gas.TPY = Tu, p, (1 - ydil) * y_fresh + ydil * prod.Y


def solve_one(case):
    tag, fuel, phi, Tu, p_bar, ydil = case
    p = p_bar * 1e5
    t0 = time.time()
    try:
        gas = make_gas()
        set_mixture(gas, fuel, phi, Tu, p, ydil, AIR_JERZEMBECK if tag == "jerzembeck-o2" else AIR)
        rho_u = gas.density
        alpha_u = gas.thermal_conductivity / (gas.density * gas.cp_mass)
        dmix = gas.mix_diff_coeffs_mass
        if phi <= 1.0:
            # deficient reactant = fuel; for blends use the mole-fraction-weighted mean D
            xs = {NAME[k]: v for k, v in fuel.items()}
            d_def = sum(v * dmix[gas.species_index(k)] for k, v in xs.items()) / sum(xs.values())
        else:
            d_def = dmix[gas.species_index("O2")]
        eq = make_gas()
        eq.TPY = gas.T, gas.P, gas.Y
        eq.equilibrate("HP")
        T_ad, rho_ad = eq.T, eq.density

        width = max(0.002, 0.03 / p_bar ** 0.7)
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
        T = f.T
        dTdx = np.max(np.abs(np.diff(T) / np.diff(f.grid)))
        Tb = float(T[-1])
        return dict(
            set=tag, fuel=fuel, phi=phi, Tu=Tu, p=p, ydil=ydil,
            SL=sig(speeds[-1]), dS=sig(abs(speeds[-1] - speeds[-2]), 3),
            Tb=sig(Tb), sigma=sig(rho_u / float(f.density[-1])),
            Tad=sig(T_ad), sigmaEq=sig(rho_u / rho_ad),
            deltaT=sig((Tb - Tu) / dTdx), alphaU=sig(alpha_u), Ddef=sig(float(d_def)), levels=len(speeds),
            points=npts[-1], seconds=round(time.time() - t0, 1),
        )
    except Exception as err:  # noqa: BLE001 — record and drop failed flames
        return dict(set=tag, fuel=fuel, phi=phi, Tu=Tu, p=p, ydil=ydil, error=str(err)[:200],
                    seconds=round(time.time() - t0, 1))


def cases():
    out = []
    for fu in ("IC8H18", "NC7H16"):
        for p in P_BAR:
            for Tu in TU:
                for phi in PHI:
                    if phi == 0.5 and Tu == 300.0 and p >= 20.0:
                        continue  # below the lean limit and irrelevant; filled by the per-φ fit
                    out.append(("grid", {fu: 1.0}, phi, Tu, p, 0.0))
    for x in (0.5, 0.9):
        for (Tu, p) in ((500.0, 5.0), (700.0, 20.0)):
            for phi in (0.8, 1.0, 1.2):
                out.append(("blend", {"IC8H18": x, "NC7H16": round(1 - x, 6)}, phi, Tu, p, 0.0))
    for fu in ("IC8H18", "NC7H16"):
        for (Tu, p) in ((500.0, 5.0), (700.0, 20.0)):
            for y in (0.1, 0.2, 0.3):
                out.append(("dilution", {fu: 1.0}, 1.0, Tu, p, y))
    for (Tu, p) in ((300.0, 1.0), (500.0, 5.0), (700.0, 20.0)):
        for phi in (0.8, 1.0, 1.2):
            out.append(("propane", {"C3H8": 1.0}, phi, Tu, p, 0.0))
    # off-grid checks at the conditions of the Jerzembeck et al. (2009) bomb experiments
    for fu in ("IC8H18", "NC7H16"):
        for p in (10.0, 25.0):
            for phi in (0.7, 0.8, 1.0, 1.2):
                out.append(("jerzembeck", {fu: 1.0}, phi, 373.0, p, 0.0))
                out.append(("jerzembeck-o2", {fu: 1.0}, phi, 373.0, p, 0.0))
    # hardest (slowest) first for better load balance: high pressure, lean/rich
    out.sort(key=lambda c: (-c[4], -abs(c[2] - 1.05)))
    return out


def case_key(tag, fuel, phi, Tu, p_bar, ydil):
    return (tag, json.dumps(fuel, sort_keys=True), round(phi, 4), float(Tu), round(p_bar, 4), float(ydil))


def main() -> None:
    """Incremental: flames already in the fixture are reused (FLAMESPEED_FULL=1 recomputes)."""
    allc = cases()
    old = {}
    if os.path.exists(OUT) and not os.environ.get("FLAMESPEED_FULL"):
        with open(OUT, encoding="utf-8") as fh:
            for r in json.load(fh)["flames"]:
                old[case_key(r["set"], r["fuel"], r["phi"], r["Tu"], r["p"] / 1e5, r["ydil"])] = r
    todo = [c for c in allc if case_key(*c) not in old]
    reused = [old[case_key(*c)] for c in allc if case_key(*c) in old]
    print(f"{len(reused)} flames reused, {len(todo)} to compute", flush=True)
    t0 = time.time()
    nproc = int(os.environ.get("FLAMESPEED_NPROC", max(1, min(8, (os.cpu_count() or 2) - 2))))
    results = list(reused)
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(solve_one, todo):
            results.append(r)
            msg = f"SL={r['SL']:.4f} dS={r['dS']:.1e} ({r['points']} pts)" if "SL" in r else f"FAILED {r['error']}"
            print(f"[{len(results)}/{len(todo) + len(reused)}] {r['set']:8s} {r['fuel']} phi={r['phi']:.2f} Tu={r['Tu']:.0f} "
                  f"p={r['p']/1e5:4.0f} bar ydil={r['ydil']:.1f}: {msg} {r['seconds']} s", flush=True)
            if len(results) % 4 == 0:
                write_fixture(results)  # checkpoint: the fixture is always valid and incremental
    ok, failed = write_fixture(results)
    print(f"wrote {os.path.relpath(OUT, ROOT)} ({len(ok)} flames, {len(failed)} failed, {time.time() - t0:.0f} s wall)")


def write_fixture(results):
    results = sorted(results, key=lambda r: (r["set"], json.dumps(r["fuel"], sort_keys=True), r["ydil"], r["p"],
                                             r["Tu"], r["phi"]))
    ok = [r for r in results if "SL" in r]
    failed = [r for r in results if "SL" not in r]
    data = dict(
        source="Cantera %s FreeFlame, PRF_HT mechanism (Jerzembeck et al., Combust. Flame 156 "
               "(2009) 292), mixture-averaged transport, air = O2 + 3.76 N2; see "
               "tools/reference/flamespeed_prf_cantera.py" % ct.__version__,
        units=dict(SL="m/s", Tu="K", p="Pa", Tb="K", deltaT="m", alphaU="m^2/s", Ddef="m^2/s"),
        cpuSeconds=round(sum(r.get("seconds", 0) for r in ok), 1),
        failed=[{k: r[k] for k in ("set", "fuel", "phi", "Tu", "p", "ydil", "error")} for r in failed],
        flames=ok,
    )
    tmp = OUT + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, separators=(",", ":"), allow_nan=False)
        fh.write("\n")
    os.replace(tmp, OUT)
    return ok, failed


if __name__ == "__main__":
    sys.exit(main())
