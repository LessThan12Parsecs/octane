"""Ignition-coil circuit oracle (independent of the TypeScript trapezoidal integrator).

Integrates the SAME circuit equations as src/physics/ignition/coil.ts with scipy's
implicit Radau IIA (5th order, adaptive, rtol 1e-10) and exact event location, as a
numerical oracle for the TS integrator:

  L1 dI1/dt + M dI2/dt = Vs - R1 I1 - v_sw          (v_sw = 0 closed, V1 open, Vclamp clamped)
  M dI1/dt + L2 dI2/dt = V2 - R2 I2
  C1 dV1/dt = I1                                     (switch open, not clamped)
  C2 dV2/dt = -I2                                    (gap open)
  V2 = -sgn(I2) * Vg(|I2|)                           (gap conducting, glow)
  Vg(i) = Vsheath + Ccol * l[mm] * p[bar]^0.51 * i^-0.32

Cases (CFR-like inductive coil, see PARAMS):
  dwell      - switch closed from rest for t_dwell (secondary open): I1(t), V2(t) samples,
               plus the analytic RL limit Vs/R1 (1 - exp(-R1 t / L1)).
  open       - switch opened at the end of the dwell, gap never breaks down: V2 ring-up
               with the collector clamp (+Vclamp / -Vrev): peak |V2|, its time, energies.
  discharge  - as 'open' but the gap breaks down at |V2| = V_bd, then a glow discharge
               until |I2| = I_ext: breakdown delay, I2 at breakdown, duration, gap and R2
               energies, I2 samples.
  lc         - k = 0, R2 = 0, gap open, V2(0) = V0: free L2-C2 ring (frequency check).

Writes test/fixtures/ignition_circuit.json.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402

PARAMS = dict(
    Vs=13.5, L1=4e-3, R1=1.5, L2=40.0, R2=8e3, C2=60e-12, k=0.98, Ilim=7.0,
    C1=1e-9, Vclamp=400.0, Vrev=24.0, t_dwell=3e-3,
    # gap (glow): normal cathode fall 3(B/A) ln(1+1/gamma), A = 14.6/(cm Torr), B = 365 V/(cm Torr)
    Vsheath=3 * 365.0 / 14.6 * math.log(1 + 1 / 0.02), Ccol=40.46, nI=-0.32, npres=0.51,
    p=10.4e5, l=0.508e-3, Vbd=9172.0, Iext=2e-3,
)
P = PARAMS
M = P["k"] * math.sqrt(P["L1"] * P["L2"])
LMAT = np.array([[P["L1"], M], [M, P["L2"]]])
LINV = np.linalg.inv(LMAT)
RTOL, ATOL = 1e-11, 1e-14


def vg(i: float) -> float:
    i = max(i, 1e-4)
    return P["Vsheath"] + P["Ccol"] * (P["l"] * 1e3) * (P["p"] * 1e-5) ** P["npres"] * i ** P["nI"]


def stored(I1, I2, V1, V2):
    return (0.5 * P["L1"] * I1 ** 2 + M * I1 * I2 + 0.5 * P["L2"] * I2 ** 2
            + 0.5 * P["C1"] * V1 ** 2 + 0.5 * P["C2"] * V2 ** 2)


# state y = [I1, I2, V1, V2, E_R1, E_R2, E_clamp, E_gap]
def rhs_factory(pmode: str, smode: str, clamp_level: float = 0.0, sgn: float = 1.0):
    def f(t, y):
        I1, I2, V1, V2 = y[0], y[1], y[2], y[3]
        if smode == "cond":
            V2 = -sgn * vg(abs(I2))
        vsw = 0.0 if pmode == "closed" else (V1 if pmode == "open" else clamp_level)
        a = P["Vs"] - P["R1"] * I1 - vsw
        b = V2 - P["R2"] * I2
        dI1, dI2 = LINV @ np.array([a, b])
        dV1 = I1 / P["C1"] if pmode == "open" else 0.0
        dV2 = -I2 / P["C2"] if smode == "open" else 0.0
        eclamp = clamp_level * I1 if pmode == "clamp" else 0.0
        egap = -V2 * I2 if smode == "cond" else 0.0
        return [dI1, dI2, dV1, dV2, P["R1"] * I1 ** 2, P["R2"] * I2 ** 2, eclamp, egap]
    return f


def run_dwell():
    y0 = [0, 0, 0, 0, 0, 0, 0, 0]
    ts = np.linspace(0, P["t_dwell"], 31)
    sol = solve_ivp(rhs_factory("closed", "open"), (0, P["t_dwell"]), y0, method="Radau",
                    rtol=RTOL, atol=ATOL, t_eval=ts, max_step=2e-6)
    assert sol.success
    analytic = [P["Vs"] / P["R1"] * (1 - math.exp(-P["R1"] * t / P["L1"])) for t in ts]
    return sol, dict(
        t=[rnd(t) for t in ts], I1=[rnd(v) for v in sol.y[0]], V2=[rnd(v) for v in sol.y[3]],
        I1_analytic=[rnd(v) for v in analytic],
        end=[rnd(v) for v in sol.y[:4, -1]],
    )


def run_after_switch_off(y_end, breakdown: bool, t_max: float):
    """Open-switch dynamics with clamp events; optional breakdown + glow until extinction."""
    y = list(y_end[:4]) + [0.0, 0.0, 0.0, 0.0]
    y[2] = 0.0  # V1 continuous from the closed switch
    t = 0.0
    pmode, smode, clamp = "open", "open", 0.0
    sgn = 1.0
    W0 = stored(*y[:4])
    samples_t, samples = [], []
    peak, tpeak = 0.0, 0.0
    t_bd, i_bd, e_dump, t_ext = None, None, 0.0, None
    guard = 0
    while t < t_max and guard < 10000:
        guard += 1
        f = rhs_factory(pmode, smode, clamp, sgn)
        events = []

        def ev_clamp_up(tt, yy):
            return yy[2] - P["Vclamp"]
        ev_clamp_up.terminal = True
        ev_clamp_up.direction = 1

        def ev_clamp_dn(tt, yy):
            return yy[2] + P["Vrev"]
        ev_clamp_dn.terminal = True
        ev_clamp_dn.direction = -1

        def ev_release(tt, yy):
            return yy[0]
        ev_release.terminal = True
        ev_release.direction = -1 if clamp > 0 else 1

        def ev_bd(tt, yy):
            return abs(yy[3]) - P["Vbd"]
        ev_bd.terminal = True
        ev_bd.direction = 1

        def ev_ext(tt, yy):
            return abs(yy[1]) - P["Iext"]
        ev_ext.terminal = True
        ev_ext.direction = -1

        if pmode == "open":
            events += [ev_clamp_up, ev_clamp_dn]
        elif pmode == "clamp":
            events += [ev_release]
        if smode == "open" and breakdown and t_bd is None:
            events += [ev_bd]
        if smode == "cond":
            events += [ev_ext]
        max_step = 0.05e-6 if smode == "open" else 1e-6
        sol = solve_ivp(f, (t, t_max), y, method="Radau", rtol=RTOL, atol=ATOL,
                        events=events, max_step=max_step, dense_output=True)
        assert sol.status >= 0, sol.message
        # peak |V2| while the gap is open
        if smode == "open":
            k = int(np.argmax(np.abs(sol.y[3])))
            if abs(sol.y[3][k]) > peak:
                peak, tpeak = abs(sol.y[3][k]), sol.t[k]
        samples_t.append(sol.t)
        samples.append(sol.y)
        t = sol.t[-1]
        y = list(sol.y[:, -1])
        if sol.status == 0:
            break
        # which event fired?
        fired = [i for i, te in enumerate(sol.t_events) if len(te) > 0][0]
        name = events[fired].__name__
        if name == "ev_clamp_up":
            pmode, clamp = "clamp", P["Vclamp"]
            y[2] = P["Vclamp"]
        elif name == "ev_clamp_dn":
            pmode, clamp = "clamp", -P["Vrev"]
            y[2] = -P["Vrev"]
        elif name == "ev_release":
            pmode = "open"
        elif name == "ev_bd":
            t_bd, i_bd = t, abs(y[1])
            sgn = 1.0 if y[1] > 0 else -1.0
            v_new = -sgn * vg(abs(y[1]))
            e_dump = 0.5 * P["C2"] * (y[3] ** 2 - v_new ** 2)
            y[3] = v_new
            smode = "cond"
        elif name == "ev_ext":
            t_ext = t
            # store V2 at extinction (algebraic value) and stop
            y[3] = -sgn * vg(abs(y[1]))
            break
    tt = np.concatenate(samples_t)
    yy = np.concatenate(samples, axis=1)
    W1 = stored(*y[:4])
    out = dict(
        peak_V2=rnd(peak), t_peak=rnd(tpeak), W0=rnd(W0), W_end=rnd(W1),
        E_R1=rnd(y[4]), E_R2=rnd(y[5]), E_clamp=rnd(y[6]),
    )
    if breakdown:
        # gap energy = dump + integral(-V2 I2) - change of the algebraic C2 energy
        v_start = -sgn * vg(i_bd)
        e_cond = y[7] - 0.5 * P["C2"] * (y[3] ** 2 - v_start ** 2)
        out.update(t_bd=rnd(t_bd), I2_bd=rnd(i_bd), E_dump=rnd(e_dump), E_gap=rnd(e_dump + e_cond),
                   t_ext=rnd(t_ext), duration=rnd(t_ext - t_bd))
        # I2 samples at fixed times after breakdown
        ts = [t_bd + d for d in (0.1e-3, 0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3, 2.5e-3)]
        out["I2_t"] = [rnd(v) for v in ts]
        out["I2_s"] = [rnd(float(abs(np.interp(v, tt, yy[1])))) for v in ts]
    else:
        ts = [2e-6, 5e-6, 10e-6, 20e-6, 40e-6, 60e-6, 80e-6]
        out["V2_t"] = [rnd(v) for v in ts]
        out["V2_s"] = [rnd(float(np.interp(v, tt, yy[3]))) for v in ts]
        out["I1_s"] = [rnd(float(np.interp(v, tt, yy[0]))) for v in ts]
    return out


def run_lc():
    L2, C2, V0 = P["L2"], P["C2"], 1000.0
    f0 = 1 / (2 * math.pi * math.sqrt(L2 * C2))
    return dict(L2=L2, C2=C2, V0=V0, f0=rnd(f0))


def main():
    sol, dwell = run_dwell()
    y_end = sol.y[:, -1]
    open_case = run_after_switch_off(y_end, breakdown=False, t_max=200e-6)
    disc = run_after_switch_off(y_end, breakdown=True, t_max=10e-3)
    write_fixture("ignition_circuit.json", dict(
        params={k: rnd(v) for k, v in P.items()},
        dwell=dwell, open=open_case, discharge=disc, lc=run_lc(),
    ))
    print("dwell I1_end", dwell["end"][0], "open peak", open_case["peak_V2"], "at", open_case["t_peak"])
    print("discharge", {k: disc[k] for k in ("t_bd", "I2_bd", "E_dump", "E_gap", "duration", "E_R2")})


if __name__ == "__main__":
    main()
