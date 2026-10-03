"""Livengood-Wu + tabulated tau vs detailed chemistry in an engine-like compression.

For CFR geometry (bore 82.55 mm, stroke 114.3 mm, rod 254 mm) at 600 rpm, an adiabatic
closed cylinder is compressed from BDC (theta = -180 deg) by the exact slider-crank volume V(t)
(Cantera IdealGasReactor + moving Wall). For each case:
  - REACTING run with the detailed mechanism -> autoignition crank angle (max dT/dt), and the
    first-stage (cool-flame) angle if present;
  - FROZEN run (all rate multipliers 0) -> the unreacted T(theta), p(theta) along which the TS
    test integrates the Livengood-Wu integral of the tabulated tau (and of tau1 for the
    two-stage variant) and predicts the autoignition angle.
This tests the LW hypothesis + table under engine-like time-varying T, p (the classic RCM /
engine validation; cf. Livengood & Wu 1955).
Writes test/fixtures/chemistry_lw_engine.json.  Mechanism: LLNL gasoline surrogate 2011 (the
default table mechanism).
"""
from __future__ import annotations

import math
import multiprocessing as mp
import os

import numpy as np

from chemistry_common import FIXTURES, LLNL_GS_YAML, to_mech_names, unburned_mixture, write_json

BORE = 3.25 * 0.0254
STROKE = 4.5 * 0.0254
ROD = 10.0 * 0.0254
RPM = 600.0
OMEGA = RPM * 2 * math.pi / 60
A_PISTON = math.pi * BORE**2 / 4
VD = A_PISTON * STROKE
TH0 = -180.0
TH1 = 40.0
DTH_OUT = 0.5

# (ON, phi, CR, T_BDC K, p_BDC Pa): HCCI-like cases chosen to autoignite near TDC
CASES = [
    (0.0, 0.5, 8.0, 360.0, 1.0e5),
    (60.0, 0.5, 10.0, 405.0, 1.0e5),
    (0.0, 1.0, 7.5, 390.0, 1.0e5),
    (90.0, 0.5, 12.0, 400.0, 1.2e5),
    (100.0, 0.5, 13.0, 420.0, 1.2e5),
    (90.0, 1.0, 10.0, 420.0, 1.5e5),
    (100.0, 1.0, 11.0, 450.0, 1.5e5),
]


def volume(theta_deg, cr):
    a = STROKE / 2
    th = math.radians(theta_deg)
    vc = VD / (cr - 1)
    return vc + A_PISTON * (ROD + a - a * math.cos(th) - math.sqrt(ROD**2 - (a * math.sin(th))**2))


def dvdt(theta_deg):
    a = STROKE / 2
    th = math.radians(theta_deg)
    s = math.sin(th)
    return A_PISTON * a * s * (1 + a * math.cos(th) / math.sqrt(ROD**2 - (a * s)**2)) * OMEGA


def run_once(case, frozen, skip=True, rtol=1e-8):
    import warnings

    import cantera as ct

    warnings.filterwarnings("ignore")
    ct.suppress_thermo_warnings()
    on, phi, cr, T0, p0 = case
    X = to_mech_names(unburned_mixture(phi, on, 0.0), lower=False)
    full = ct.Solution(LLNL_GS_YAML)
    if frozen:
        # non-reacting: only the species present, same thermo data
        gas = ct.Solution(thermo="ideal-gas", species=[full.species(k) for k in X])
    else:
        gas = full
    gas.TPX = T0, p0, X
    r = ct.IdealGasMoleReactor(gas, clone=False)
    r.volume = volume(TH0, cr)
    env = ct.Reservoir(ct.Solution("air.yaml"))
    # velocity of the wall into the environment = dV/dt / A (A = 1 m^2): left volume grows
    ct.Wall(r, env, A=1.0, velocity=lambda t: dvdt(TH0 + math.degrees(OMEGA * t)))
    net = ct.ReactorNet([r])
    net.rtol, net.atol = rtol, 1e-15
    if not frozen:
        net.preconditioner = ct.AdaptivePreconditioner()
        if skip:
            net.derivative_settings = {"skip-third-bodies": True, "skip-falloff": True}
    t_end = math.radians(TH1 - TH0) / OMEGA
    ths = np.arange(TH0, TH1 + 1e-9, DTH_OUT)
    out_T, out_p = [], []
    ts, Ts = [0.0], [T0]
    for th in ths[1:] if frozen else []:
        net.advance(math.radians(th - TH0) / OMEGA)
        out_T.append(r.T)
        out_p.append(r.phase.P)
    if frozen:
        return dict(T=[T0] + out_T, p=[p0] + out_p)
    while net.time < t_end:
        net.step()
        ts.append(net.time)
        Ts.append(r.T)
        if r.T > 2200:
            break
    ts, Ts = np.array(ts), np.array(Ts)
    rate = np.diff(Ts) / np.diff(ts)
    im = int(np.argmax(rate))
    ignited = Ts.max() > 1500
    th_ign = TH0 + math.degrees(OMEGA * 0.5 * (ts[im] + ts[im + 1])) if ignited else None
    # first stage: dT/dt local max with a clear dip afterwards (as in chemistry_common)
    th1 = None
    run_max, i1 = 0.0, -1
    for i in range(im):
        if rate[i] > run_max:
            run_max, i1 = rate[i], i
        elif i1 >= 0 and rate[i] < 0.3 * run_max and run_max > 5e3:
            th1 = TH0 + math.degrees(OMEGA * 0.5 * (ts[i1] + ts[i1 + 1]))
            break
    return dict(thetaIgn=th_ign, thetaFirstStage=th1, Tmax=float(Ts.max()))


def run(case, frozen):
    """run_once with fallbacks for CVODES corrector failures at the ignition front."""
    import cantera as ct

    last = None
    for skip, rtol in ((True, 1e-8), (False, 1e-8), (False, 1e-7), (False, 1e-9)):
        try:
            return run_once(case, frozen, skip, rtol)
        except ct.CanteraError as e:  # pragma: no cover - solver robustness
            last = e
    raise last


def _job(args):
    case, frozen = args
    return case, frozen, run(case, frozen)


def main():
    jobs = [(c, f) for c in CASES for f in (True, False)]
    res = {}
    with mp.Pool(int(os.environ.get("OCTANE_WORKERS", "4"))) as pool:
        for case, frozen, r in pool.imap_unordered(_job, jobs):
            res[(case, frozen)] = r
    out = []
    for c in CASES:
        fr, rx = res[(c, True)], res[(c, False)]
        on, phi, cr, T0, p0 = c
        print(f"ON={on} phi={phi} CR={cr} T0={T0}: ignition {rx['thetaIgn']} deg, first stage "
              f"{rx['thetaFirstStage']}, T_TDC(frozen) {fr['T'][int((0 - TH0) / DTH_OUT)]:.0f} K")
        out.append(dict(on=on, phi=phi, cr=cr, T0=T0, p0=p0, rpm=RPM, theta0=TH0, dTheta=DTH_OUT,
                        T=[round(x, 4) for x in fr["T"]], p=[round(x, 2) for x in fr["p"]],
                        thetaIgn=rx["thetaIgn"], thetaFirstStage=rx["thetaFirstStage"]))
    write_json(os.path.join(FIXTURES, "chemistry_lw_engine.json"),
               dict(geometry=dict(bore=BORE, stroke=STROKE, rod=ROD), mechanism="llnl-gasoline-2011",
                    cases=out))


if __name__ == "__main__":
    main()
