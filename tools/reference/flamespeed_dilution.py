"""Residual/EGR-diluted laminar flames for the dilution model of
src/physics/combustion/laminar-flame-speed.ts.

Why: the model's thermal-theory factor F = exp[−E_a/(2R) (1/T_ad(x) − 1/T_ad(0))] was first
calibrated on STOICHIOMETRIC diluted flames only (flamespeed_cantera.py / flamespeed_prf_cantera.py
"dilution" sets). Diluted flames at φ = 0.8 / 1.2 showed that a single E_a(p) per fuel
under-predicts S_L(x)/S_L(0) of off-stoichiometric mixtures by up to 0.1 (35 % relative at
x = 0.3), so the effective activation energy must depend on φ (and on the thermodynamic state).

Diluent definition (identical to the model: laminar-flame-speed.ts `unburnedComposition`,
thermo/fuels.ts `completeCombustionMoles`): the complete-combustion products of the SAME
fuel–air mixture (lean: CO2, H2O, O2, N2; rich: CO2, CO, H2O, H2, N2 with the water-gas shift
frozen at 1740 K, flamespeed_fit.complete_products), at MASS fraction x of the charge.

Sets (every diluted flame is paired with an undiluted flame at the same (φ, Tu, p) computed with
the same refinement sequence, so ratios carry little discretisation bias):
  fit:   CH4 (GRI-Mech 3.0), IC8H18, NC7H16, C3H8 (PRF_HT, Jerzembeck et al. 2009) ×
         φ ∈ {0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4} (the φ nodes of the S_L tables in 0.6-1.4:
         the effective E_a has a narrow maximum near φ = 1, so it must not be interpolated in φ
         between 0.8 and 1.2) × (Tu, p) ∈ {(300 K, 1 bar), (500 K, 5 bar), (700 K, 20 bar),
         (900 K, 60 bar)}
         × x ∈ {0, 0.15, 0.3}   → E_a(φ, p) tables (flamespeed_fit.py)
  check: conditions OFF the fitting set (other φ, Tu, p, x) for independent validation in the
         TypeScript tests.
Numerics as in the parent oracles: mixture-averaged transport, air = O2 + 3.76 N2, refine
criteria (slope, curve) (0.06, 0.12) → (0.03, 0.06) → (0.015, 0.03), stop when S_L changes by
< 0.3 % between levels. Flames that fail to converge are recorded under "failed" and skipped.

Writes test/fixtures/flamespeed_dilution.json (incremental: existing flames are reused).
Run: .venv/bin/python tools/reference/flamespeed_dilution.py   (≈ 3 h on 7 cores under load)
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
import time

import cantera as ct

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import flamespeed_prf_cantera as prf  # noqa: E402
from flamespeed_fit import NU, complete_products  # noqa: E402

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "test", "fixtures", "flamespeed_dilution.json")
LEVELS = [(3.0, 0.06, 0.12, 0.03), (2.0, 0.03, 0.06, 0.015), (2.0, 0.015, 0.03, 0.0075)]
CONV_TOL = 0.003
MECH_NAME = {"CH4": "CH4", "IC8H18": "IXC8H18", "NC7H16": "NXC7H16", "C3H8": "C3H8"}
MECH = {"CH4": "GRI-Mech 3.0", "IC8H18": "PRF_HT", "NC7H16": "PRF_HT", "C3H8": "PRF_HT"}

FIT_PHI = [0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4]
FIT_TP = [(300.0, 1.0), (500.0, 5.0), (700.0, 20.0), (900.0, 60.0)]
FIT_X = [0.0, 0.15, 0.3]


def cases():
    out = []
    for fuel in ("CH4", "IC8H18", "NC7H16", "C3H8"):
        for phi in FIT_PHI:
            for Tu, p in FIT_TP:
                for x in FIT_X:
                    out.append(("fit", fuel, phi, Tu, p, x))
    chk = []
    chk += [("CH4", phi, 600.0, 40.0, x) for phi in (0.7, 1.1) for x in (0.0, 0.1, 0.25)]
    chk += [("CH4", 0.9, 400.0, 2.0, x) for x in (0.0, 0.2)]
    chk += [("IC8H18", 1.0, 373.0, 10.0, x) for x in (0.0, 0.2)]
    chk += [("IC8H18", phi, 800.0, 40.0, x) for phi in (0.7, 1.3) for x in (0.0, 0.2)]
    chk += [("NC7H16", 0.9, 600.0, 10.0, x) for x in (0.0, 0.1, 0.25)]
    chk += [("C3H8", 1.1, 450.0, 3.0, x) for x in (0.0, 0.2)]
    out += [("check",) + c for c in chk]
    return out


def make_gas(fuel: str) -> ct.Solution:
    if fuel == "CH4":
        gas = ct.Solution("gri30.yaml")
        gas.transport_model = "mixture-averaged"
        return gas
    return prf.make_gas()


def fresh_mole_fractions(fuel: str, phi: float, x: float) -> dict:
    """Air (O2 + 3.76 N2) + fuel at φ + complete-combustion products of that same mixture at
    diluent MASS fraction x (products of N have the mass of N, so moles scale by x/(1−x))."""
    N = {"O2": 1.0, "N2": 3.76, fuel: phi / NU[fuel]}
    if x > 0:
        prod = complete_products(N)
        scale = x / (1 - x)
        for k, v in prod.items():
            N[k] = N.get(k, 0.0) + scale * v
    tot = sum(N.values())
    return {MECH_NAME.get(k, k): v / tot for k, v in N.items() if v > 0}


def solve(case):
    tag, fuel, phi, Tu, p_bar, x = case
    t0 = time.time()
    base = dict(set=tag, fuel=fuel, phi=phi, Tu=Tu, p=p_bar * 1e5, ydil=x, mech=MECH[fuel])
    try:
        gas = make_gas(fuel)
        gas.TPX = Tu, p_bar * 1e5, fresh_mole_fractions(fuel, phi, x)
        width = max(0.002, 0.03 / p_bar ** 0.7)
        f = ct.FreeFlame(gas, width=width)
        f.transport_model = "mixture-averaged"
        f.set_max_grid_points(f.flame, 4000)
        speeds = []
        for i, (ratio, slope, curve, prune) in enumerate(LEVELS):
            if i >= 2 and abs(speeds[-1] - speeds[-2]) <= CONV_TOL * speeds[-1]:
                break
            f.set_refine_criteria(ratio=ratio, slope=slope, curve=curve, prune=prune)
            f.solve(loglevel=0, auto=(i == 0))
            speeds.append(float(f.velocity[0]))
        return dict(base, SL=float(f"{speeds[-1]:.6g}"), dS=float(f"{abs(speeds[-1] - speeds[-2]):.3g}"),
                    levels=len(speeds), points=len(f.grid), seconds=round(time.time() - t0, 1))
    except Exception as err:  # noqa: BLE001 — record and skip
        return dict(base, error=str(err)[:200], seconds=round(time.time() - t0, 1))


def key(tag, fuel, phi, Tu, p_bar, x):
    return (tag, fuel, round(phi, 4), float(Tu), round(p_bar, 4), float(x))


def main() -> None:
    allc = cases()
    old = {}
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as fh:
            for r in json.load(fh)["flames"]:
                old[key(r["set"], r["fuel"], r["phi"], r["Tu"], r["p"] / 1e5, r["ydil"])] = r
    todo = [c for c in allc if key(*c) not in old]
    results = [old[key(*c)] for c in allc if key(*c) in old]
    failed = []
    print(f"{len(results)} reused, {len(todo)} to compute", flush=True)
    todo.sort(key=lambda c: (c[1] == "CH4", -c[4], c[5] == 0))  # slow PRF_HT / high-p first
    nproc = int(os.environ.get("FLAMESPEED_NPROC", "6"))
    with mp.Pool(nproc) as pool:
        for r in pool.imap_unordered(solve, todo):
            print(r, flush=True)
            (results if "SL" in r else failed).append(r)
            if "SL" in r and len(results) % 3 == 0:
                write(results, failed)
    write(results, failed)


def write(results, failed) -> None:
    results = sorted(results, key=lambda r: (r["set"], r["fuel"], r["phi"], r["Tu"], r["p"], r["ydil"]))
    data = dict(
        source="Cantera %s FreeFlame (GRI-Mech 3.0 for CH4; PRF_HT of Jerzembeck et al. 2009 for "
               "IC8H18, NC7H16, C3H8), mixture-averaged, air = O2 + 3.76 N2, diluent = "
               "complete-combustion products of the same mixture at mass fraction ydil; see "
               "tools/reference/flamespeed_dilution.py" % ct.__version__,
        units=dict(SL="m/s", Tu="K", p="Pa", ydil="mass fraction"),
        cpuSeconds=round(sum(r.get("seconds", 0) for r in results), 1),
        failed=[{k: r[k] for k in ("set", "fuel", "phi", "Tu", "p", "ydil", "error")} for r in failed],
        flames=results,
    )
    tmp = OUT + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(data, fh, separators=(",", ":"), allow_nan=False)
        fh.write("\n")
    os.replace(tmp, OUT)


if __name__ == "__main__":
    main()
