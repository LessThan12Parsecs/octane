"""Validation data for the tabulated PRF ignition-delay models (LLNL mechanisms).

Writes test/fixtures/chemistry_ignition_validation.json with, per mechanism
(llnl-gasoline-2011 = Mehl et al. 2011, llnl-prf-v2 = Curran et al. 2002):

fieweger   Shock-tube n-heptane/air ignition delays at 40 bar, phi = 1 (Fieweger, Blumenthal &
           Adomeit, Combust. Flame 109 (1997) 599; tabular values from the ChemKED database
           file n-heptane/Fieweger 1997/st_fieweger_1997.yaml, converted from ReSpecTh
           st_fieweger_1997.xml; copy in tools/reference/mechanisms/experimental/), with the
           mechanism's constant-volume prediction for the exact experimental mixture
           (1.874 % nC7H16, O2 20.615 %, N2 77.511 %; criterion max dT/dt ~ max dp/dt).
           The Fieweger iso-octane and PRF 60/80/90 data could not be obtained in tabular
           form (ReSpecTh requires a login; ChemKED only holds n-heptane).
offgrid    150 random points inside the table domain (numpy default_rng(20260929)) with
           direct Cantera ignition delays -> interpolation error of the table.
sweep      T sweeps 650-1100 K at 40 bar, phi = 1, x_res = 0, for PRF 0/80/90/100 (direct
           Cantera; tau and tau1) -> Douaud-Eyzat vs detailed chemistry (NTC).
offgridLow (llnl-gasoline-2011 only, the extended-grid default table) 60 random points in the
           compression-stroke region T 550-700 K, p 3.54-10 bar (numpy default_rng(20260930),
           phi 0.5-1.5, PRF 0-100, x_res 0-0.15), plus the 24-point probe T = 550/600/650/700 K,
           p = 3/6 bar, PRF 0/90/100, phi = 1, x_res = 0 that exposed the extrapolation error
           of the legacy (650 K / 10 bar-bounded) table; t_max = 100 s.
trajectory (llnl-gasoline-2011 only) direct Cantera tau (t_max 100 s) at crank angles
           -60/-40/-30/-20/-10/0 deg on the 7 frozen engine compressions of
           test/fixtures/chemistry_lw_engine.json (unburned T, p, phi, PRF; x_res = 0): checks the
           table exactly where the Livengood-Wu integral is accumulated.
A part whose mechanism YAML hash matches the existing fixture entry is not recomputed.
"""
from __future__ import annotations

import json
import math
import multiprocessing as mp
import os

import numpy as np

from chemistry_common import (FIXTURES, LLNL_GS_YAML, LLNL_PRF_YAML, MECH_DIR, file_sha256,
                              ignition_delay_robust, to_mech_names, unburned_mixture, write_json)

FIEWEGER = os.path.join(MECH_DIR, "experimental", "fieweger1997_nc7h16_chemked.yaml")
OUT = os.path.join(FIXTURES, "chemistry_ignition_validation.json")
MECHS = {"llnl-gasoline-2011": (LLNL_GS_YAML, False), "llnl-prf-v2": (LLNL_PRF_YAML, True)}

_gas = None


def _init(yaml):
    global _gas
    import warnings

    import cantera as ct

    ct.suppress_thermo_warnings()
    warnings.filterwarnings("ignore")
    _gas = ct.Solution(yaml)


def _run(args):
    kind, T, p, X = args[:4]
    t_max = args[4] if len(args) > 4 else 2.0
    r = ignition_delay_robust(_gas, T, p, X, t_max=t_max)
    return kind, T, p, r


def read_fieweger():
    from ruamel.yaml import YAML

    d = YAML(typ="safe").load(open(FIEWEGER))
    comp = {s["species-name"]: float(s["amount"][0]) for s in d["common-properties"]["composition"]["species"]}
    pts = []
    for dp in d["datapoints"]:
        tau_us = float(str(dp["ignition-delay"][0]).split()[0])
        T = float(str(dp["temperature"][0]).split()[0])
        p = float(str(dp["pressure"][0]).split()[0]) * 1e5
        pts.append((T, p, tau_us * 1e-6))
    return comp, pts


def offgrid_points():
    rng = np.random.default_rng(20260929)
    off = []
    for _ in range(150):
        invT = rng.uniform(1000 / 1100, 1000 / 650)
        lnp = rng.uniform(math.log(10e5), math.log(80e5))
        off.append(dict(T=1000 / invT, p=math.exp(lnp), phi=rng.uniform(0.5, 1.5),
                        on=rng.uniform(0, 100), xres=rng.uniform(0, 0.15)))
    return off


def offgrid_low_points():
    rng = np.random.default_rng(20260930)
    off = []
    for _ in range(60):
        invT = rng.uniform(1000 / 700, 1000 / 550)
        lnp = rng.uniform(math.log(10e5 / math.sqrt(8)), math.log(10e5))
        off.append(dict(T=1000 / invT, p=math.exp(lnp), phi=rng.uniform(0.5, 1.5),
                        on=rng.uniform(0, 100), xres=rng.uniform(0, 0.15)))
    for T in (550.0, 600.0, 650.0, 700.0):
        for p in (3e5, 6e5):
            for on in (0.0, 90.0, 100.0):
                off.append(dict(T=T, p=p, phi=1.0, on=on, xres=0.0))
    return off


def run_offgrid_low(yaml, lower, workers):
    names = lambda X: to_mech_names(X, lower=lower)  # noqa: E731
    off = offgrid_low_points()
    jobs = [(("l", i), o["T"], o["p"], names(unburned_mixture(o["phi"], o["on"], o["xres"])), 100.0)
            for i, o in enumerate(off)]
    res = {}
    with mp.Pool(workers, initializer=_init, initargs=(yaml,)) as pool:
        for kind, T, p, r in pool.imap_unordered(_run, jobs, chunksize=1):
            res[kind] = r
            print(f"  low {kind[1]}: T={T:.0f} p={p/1e5:.2f} bar tau={r['tau']:.4g}", flush=True)
    return [dict(o, tau=res[("l", i)]["tau"], tau1=res[("l", i)]["tau1"], ignited=res[("l", i)]["ignited"])
            for i, o in enumerate(off)]


TRAJ_THETA = [-60.0, -40.0, -30.0, -20.0, -10.0, 0.0]


def run_trajectory(yaml, lower, workers):
    eng = json.load(open(os.path.join(FIXTURES, "chemistry_lw_engine.json")))
    names = lambda X: to_mech_names(X, lower=lower)  # noqa: E731
    pts = []
    for ci, c in enumerate(eng["cases"]):
        for th in TRAJ_THETA:
            i = int(round((th - c["theta0"]) / c["dTheta"]))
            pts.append(dict(case=ci, theta=th, T=c["T"][i], p=c["p"][i], phi=c["phi"], on=c["on"], xres=0.0))
    jobs = [(("t", i), o["T"], o["p"], names(unburned_mixture(o["phi"], o["on"], 0.0)), 100.0)
            for i, o in enumerate(pts)]
    res = {}
    with mp.Pool(workers, initializer=_init, initargs=(yaml,)) as pool:
        for kind, T, p, r in pool.imap_unordered(_run, jobs, chunksize=1):
            res[kind] = r
    return [dict(o, tau=res[("t", i)]["tau"], tau1=res[("t", i)]["tau1"], ignited=res[("t", i)]["ignited"])
            for i, o in enumerate(pts)]


SWEEP_T = [1000 / x for x in np.linspace(1000 / 1100, 1000 / 650, 25)]
SWEEP_ON = [0.0, 80.0, 90.0, 100.0]


def run_mechanism(yaml, lower, comp, fpts, off, workers):
    names = lambda X: to_mech_names(X, lower=lower)  # noqa: E731
    Xf = names({"NC7H16": comp["nC7H16"], "O2": comp["O2"], "N2": comp["N2"]})
    jobs = []
    for i, (T, p, _) in enumerate(fpts):
        jobs.append((("f", i), T, p, Xf))
    for i, o in enumerate(off):
        jobs.append((("o", i), o["T"], o["p"], names(unburned_mixture(o["phi"], o["on"], o["xres"]))))
    for j, on in enumerate(SWEEP_ON):
        for i, T in enumerate(SWEEP_T):
            jobs.append((("s", j * 1000 + i), T, 40e5, names(unburned_mixture(1.0, on, 0.0))))
    res = {}
    with mp.Pool(workers, initializer=_init, initargs=(yaml,)) as pool:
        for kind, T, p, r in pool.imap_unordered(_run, jobs, chunksize=1):
            res[kind] = r
    fie = []
    for i, (T, p, tau) in enumerate(fpts):
        r = res[("f", i)]
        fie.append(dict(T=T, p=p, tauExp=tau, tau=r["tau"], tau1=r["tau1"]))
        print(f"  Fieweger nC7 T={T:.0f} K: exp {tau*1e6:.0f} us, model {r['tau']*1e6:.0f} us "
              f"(ratio {r['tau']/tau:.2f})")
    offgrid = [dict(o, tau=res[("o", i)]["tau"], tau1=res[("o", i)]["tau1"], ignited=res[("o", i)]["ignited"])
               for i, o in enumerate(off)]
    sweep = []
    for j, on in enumerate(SWEEP_ON):
        sweep.append(dict(on=on, p=40e5, phi=1.0, xres=0.0, T=SWEEP_T,
                          tau=[res[("s", j * 1000 + i)]["tau"] for i in range(len(SWEEP_T))],
                          tau1=[res[("s", j * 1000 + i)]["tau1"] for i in range(len(SWEEP_T))]))
    return dict(yamlSha256=file_sha256(yaml), fieweger=fie, offgrid=offgrid, sweep=sweep)


def main():
    comp, fpts = read_fieweger()
    off = offgrid_points()
    old = {}
    if os.path.exists(OUT):
        try:
            old = json.load(open(OUT)).get("mechanisms", {})
        except Exception:
            old = {}
    workers = int(os.environ.get("OCTANE_WORKERS", "3"))
    mechs = {}
    for name, (yaml, lower) in MECHS.items():
        if name in old and old[name].get("yamlSha256") == file_sha256(yaml):
            print(f"[{name}] up to date")
            mechs[name] = old[name]
        else:
            print(f"[{name}] running")
            mechs[name] = run_mechanism(yaml, lower, comp, fpts, off, workers)
        if name == "llnl-gasoline-2011" and "offgridLow" not in mechs[name]:
            print(f"[{name}] running offgridLow")
            mechs[name]["offgridLow"] = run_offgrid_low(yaml, lower, workers)
        if name == "llnl-gasoline-2011" and "trajectory" not in mechs[name]:
            print(f"[{name}] running trajectory")
            mechs[name]["trajectory"] = run_trajectory(yaml, lower, workers)
    write_json(OUT, dict(
        fieweger=dict(source="Fieweger, Blumenthal & Adomeit, Combust. Flame 109 (1997) 599; "
                             "ChemKED n-heptane/Fieweger 1997 (ReSpecTh st_fieweger_1997.xml)",
                      mixture=comp, T=[p[0] for p in fpts], p=[p[1] for p in fpts],
                      tauExp=[p[2] for p in fpts]),
        mechanisms=mechs))


if __name__ == "__main__":
    main()
