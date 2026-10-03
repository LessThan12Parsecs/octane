"""PRF ignition-delay tables from detailed chemistry (LLNL mechanisms).

Writes src/physics/chemistry/data/prf-ignition-<mechanism>.json for the mechanisms in MECHANISMS
(LLNL gasoline surrogate 2011 = default model, LLNL PRF v2), consumed by
src/physics/chemistry/ignition-delay-llnl.ts (tabulated IgnitionDelayModel).

Each entry: adiabatic CONSTANT-VOLUME 0-D ignition (Cantera 3.2 IdealGasMoleReactor with
AdaptivePreconditioner, rtol 1e-6, atol 1e-15; tighter on CVODES failure) from (T, p) of the
unburned mixture
    X = (1 - x_res) * [dry air + PRF vapour at phi] + x_res * [complete-combustion products]
(mole fractions; compositions exactly as src/physics/thermo/fuels.ts, see chemistry_common.py).
  tau  = time of max dT/dt (main/hot ignition) -- the definition used for reflected-shock data
         (max dp/dt, e.g. Fieweger et al. 1997);
  tau1 = first-stage (cool-flame) dT/dt maximum when a distinct first stage exists, else tau.
tau is capped at the grid's t_max; capped entries are listed in "capped".

Grids (multilinear interpolation in ln tau in the TS model):
  'extended' (llnl-gasoline-2011, the default model), t_max = 100 s:
    1000/T : the 17 'legacy' nodes (uniform, 1100 -> 650 K) + 600 K and 550 K. The end gas spends
             most of the compression stroke below 650 K / 10 bar, where extrapolating the legacy
             table made tau 2-30x too short (verified against direct Cantera runs: 550-700 K,
             3-6 bar, PRF 0/90/100) and added 0.05-0.4 of spurious Livengood-Wu integral before
             650 K was even reached.
    ln p   : 7 nodes uniform (step ln(8)/4) from 10 bar * 8^(-1/2) = 3.54 bar to 80 bar
             (the legacy 10-80 bar nodes + 5.95 and 3.54 bar; MON-like compressions pass
             550 K at ~3 bar).
  'legacy' (llnl-prf-v2, comparison only), t_max = 2 s: 1000/T 17 nodes uniform 1100 -> 650 K,
    ln p 5 nodes uniform 10 -> 80 bar.
  Both: phi 0.5 0.75 1.0 1.25 1.5; ON 0 20 40 60 70 80 90 100 (PRF liquid-volume % iso-octane);
  x_res 0 0.15 (residual MOLE fraction of the unburned gas; ln tau is nearly linear in dilution
  -- checked against off-grid points).
Axes are written as {"nodes": [...]} (extended) or {"start", "step", "n"} (legacy, uniform).
Flat index = ((((iT * nP + iP) * nPhi + iPhi) * nOn + iOn) * nRes + iRes).

Runtime: 1-17 s per reactor run (1382 species); the extended grid has 10640 points. Results are
cached per point in a JSONL file in the system temp dir (keyed by mechanism hash and the point),
so extending a grid only computes the new points; a table whose inputHash matches is skipped.
Usage: .venv/bin/python tools/reference/chemistry_ignition_table.py [--force] [--workers N]
       [--mech llnl-gasoline-2011|llnl-prf-v2]
"""
from __future__ import annotations

import hashlib
import json
import math
import multiprocessing as mp
import os
import sys
import tempfile
import time
import zlib

from chemistry_common import (CHEM_DATA_DIR, LLNL_GS_YAML, LLNL_PRF_YAML, file_sha256,
                              ignition_delay_robust, to_mech_names, unburned_mixture, write_json)

MECHANISMS = {
    # default model in the TS code: updated PRF chemistry, better vs Fieweger et al. (1997)
    "llnl-gasoline-2011": dict(
        yaml=LLNL_GS_YAML, lower=False, out="prf-ignition-llnl-gasoline-2011.json", grid="extended",
        cite=("LLNL gasoline-surrogate detailed mechanism (Mehl, Pitz, Westbrook & Curran, Proc. Combust. "
              "Inst. 33 (2011) 193; LLNL-MI-536371), PRF sub-mechanisms; 7 pentene LT species without "
              "thermo removed")),
    "llnl-prf-v2": dict(
        yaml=LLNL_PRF_YAML, lower=True, out="prf-ignition-llnl-prf-v2.json", grid="legacy",
        cite=("LLNL PRF v2 detailed mechanism (Curran, Gaffuri, Pitz & Westbrook, Combust. Flame 129 "
              "(2002) 253; LLNL prf_2d, UCRL-WEB-208393)")),
}

# legacy grid (VERSION 1; its input hash is unchanged so existing legacy tables stay valid)
VERSION = 1
T_MAX = 2.0
N_T = 17
INV_T0 = 1000.0 / 1100.0
INV_T1 = 1000.0 / 650.0
N_P = 5
LNP0 = math.log(10e5)
LNP1 = math.log(80e5)
PHI = [0.5, 0.75, 1.0, 1.25, 1.5]
ON = [0.0, 20.0, 40.0, 60.0, 70.0, 80.0, 90.0, 100.0]
XRES = [0.0, 0.15]
D_INV = (INV_T1 - INV_T0) / (N_T - 1)
D_LNP = (LNP1 - LNP0) / (N_P - 1)

# extended grid (see module docstring)
EXT_VERSION = 2
EXT_T_MAX = 100.0
EXT_INV_T = [INV_T0 + i * D_INV for i in range(N_T)] + [1000.0 / 600.0, 1000.0 / 550.0]
EXT_LNP = [LNP0 + i * D_LNP for i in range(-2, N_P)]


def grid_spec(kind: str):
    """(invT nodes, lnP nodes, t_max, axes JSON for invT, axes JSON for lnP)."""
    if kind == "legacy":
        inv = [INV_T0 + i * D_INV for i in range(N_T)]
        lnp = [LNP0 + i * D_LNP for i in range(N_P)]
        return (inv, lnp, T_MAX,
                dict(start=INV_T0, step=D_INV, n=N_T, unit="1000/T, 1/K"),
                dict(start=LNP0, step=D_LNP, n=N_P, unit="ln(p/Pa)"))
    return (EXT_INV_T, EXT_LNP, EXT_T_MAX,
            dict(nodes=EXT_INV_T, unit="1000/T, 1/K"), dict(nodes=EXT_LNP, unit="ln(p/Pa)"))


def grid(kind: str):
    inv, lnp, _, _, _ = grid_spec(kind)
    pts = []
    for it in inv:
        T = 1000.0 / it
        for lp in lnp:
            p = math.exp(lp)
            for phi in PHI:
                for on in ON:
                    for xr in XRES:
                        pts.append((T, p, phi, on, xr))
    return pts


def input_hash(yaml: str, kind: str) -> str:
    if kind == "legacy":
        spec = json.dumps(dict(v=VERSION, tmax=T_MAX, nT=N_T, T=[INV_T0, INV_T1], nP=N_P,
                               P=[LNP0, LNP1], phi=PHI, on=ON, xres=XRES,
                               mech=file_sha256(yaml)), sort_keys=True)
    else:
        spec = json.dumps(dict(v=EXT_VERSION, tmax=EXT_T_MAX, invT=EXT_INV_T, lnP=EXT_LNP, phi=PHI,
                               on=ON, xres=XRES, mech=file_sha256(yaml)), sort_keys=True)
    return hashlib.sha256(spec.encode()).hexdigest()[:16]


def point_key(pt) -> str:
    T, p, phi, on, xr = pt
    return f"{T:.9g}|{p:.9g}|{phi:g}|{on:g}|{xr:g}"


def cache_path(yaml: str) -> str:
    return os.path.join(tempfile.gettempdir(), f"octane_prf_ignition_pts_{file_sha256(yaml)[:16]}.jsonl")


def load_cache(path: str) -> dict:
    """{point_key: (tau, tau1, ignited, t_max)} from the per-point JSONL cache."""
    done = {}
    if os.path.exists(path):
        for line in open(path):
            try:
                k, tau, tau1, ign, tmax = json.loads(line)
                done[k] = (tau, tau1, ign, tmax)
            except Exception:
                pass
    return done


_gas = None
_lower = True
_tmax = T_MAX


def _init(yaml: str, lower: bool, tmax: float):
    global _gas, _lower, _tmax
    import cantera as ct

    ct.suppress_thermo_warnings()
    import warnings

    warnings.filterwarnings("ignore")
    _gas = ct.Solution(yaml)
    _lower = lower
    _tmax = tmax


def _run(args):
    key, (T, p, phi, on, xr) = args
    X = to_mech_names(unburned_mixture(phi, on, xr), lower=_lower)
    r = ignition_delay_robust(_gas, T, p, X, t_max=_tmax)
    return key, r["tau"], r["tau1"], r["ignited"], _tmax


def main():
    force = "--force" in sys.argv
    workers = 6
    if "--workers" in sys.argv:
        workers = int(sys.argv[sys.argv.index("--workers") + 1])
    names = list(MECHANISMS)
    if "--mech" in sys.argv:
        names = [sys.argv[sys.argv.index("--mech") + 1]]
    for name in names:
        build(name, MECHANISMS[name], force, workers)


def build(name: str, mech: dict, force: bool, workers: int) -> None:
    OUT = os.path.join(CHEM_DATA_DIR, mech["out"])
    kind = mech["grid"]
    h = input_hash(mech["yaml"], kind)
    if not force and os.path.exists(OUT):
        try:
            if json.load(open(OUT))["inputHash"] == h:
                print("up to date:", os.path.relpath(OUT))
                return
        except Exception:
            pass
    _, _, tmax, axT, axP = grid_spec(kind)
    pts = grid(kind)
    keys = [point_key(pt) for pt in pts]
    ckpt = cache_path(mech["yaml"])
    done = load_cache(ckpt)
    # a cached result is reusable if it ignited, or was run at least as long as t_max
    ok = {k: v for k, v in done.items() if v[2] or v[3] >= tmax}
    todo = [(keys[i], pts[i]) for i in range(len(pts)) if keys[i] not in ok]
    # interleave so slow low-T cases are spread over the run
    todo.sort(key=lambda a: zlib.crc32(a[0].encode()))
    print(f"[{name}] {len(pts)} points, {len(pts) - len(todo)} cached, {len(todo)} to run on {workers} "
          f"workers; cache {ckpt}")
    t0 = time.time()
    if todo:
        with mp.Pool(workers, initializer=_init, initargs=(mech["yaml"], mech["lower"], tmax)) as pool, \
                open(ckpt, "a") as f:
            for n, (k, tau, tau1, ign, tm) in enumerate(pool.imap_unordered(_run, todo, chunksize=1)):
                ok[k] = (tau, tau1, ign, tm)
                f.write(json.dumps([k, tau, tau1, ign, tm]) + "\n")
                f.flush()
                if n % 100 == 0:
                    el = time.time() - t0
                    print(f"  {n + 1}/{len(todo)}  {el:.0f}s elapsed, eta {el / (n + 1) * (len(todo) - n - 1):.0f}s",
                          flush=True)
    lnTau = [round(math.log(ok[k][0]), 4) for k in keys]
    lnTau1 = [round(math.log(ok[k][1]), 4) for k in keys]
    capped = [i for i, k in enumerate(keys) if not ok[k][2]]
    data = dict(
        description=("ln(ignition delay / s) of PRF/air/residual mixtures, adiabatic constant-volume "
                     "0-D reactor, " + mech["cite"] + ", Cantera 3.2. "
                     "tau = max dT/dt; tau1 = first-stage dT/dt max (= tau if single-stage). "
                     "Generated by tools/reference/chemistry_ignition_table.py."),
        mechanism=name,
        inputHash=h,
        reactor="constant-volume adiabatic",
        tauMax=tmax,
        axes=dict(invT=axT, lnP=axP, phi=PHI, on=ON, xres=XRES),
        order=["invT", "lnP", "phi", "on", "xres"],
        lnTau=lnTau,
        lnTau1=lnTau1,
        capped=capped,
    )
    write_json(OUT, data)
    print(f"done in {time.time() - t0:.0f}s; {len(capped)} capped at {tmax}s")


if __name__ == "__main__":
    main()
