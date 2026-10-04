"""Trembler-coil (Ford Model T) ignition oracle, independent of the TypeScript trapezoidal integrator.

Integrates the SAME equations as src/physics/ignition/trembler-coil.ts + vibrator.ts with scipy's
implicit Radau IIA (adaptive, rtol 1e-9) and exact event location:

  Lp dI1/dt + M dI2/dt = e(t) - R I1 - v_pts        (timer closed; v_pts = V1 points open, 0 closed)
  I1 = 0                                            (timer open; at the break I2 += M I1 / L2)
  M dI1/dt + L2 dI2/dt = V2 - R2 I2
  C1 dV1/dt = I1                                    (points open; V1 = 0 closed, discharged on closing)
  C2 dV2/dt = -I2                                   (gap open)
  V2 = -sgn(I2) Vg(|I2|)                            (gap conducting, glow)
  x'' = w^2 [ xp ((I1/Ip)^2 (g0/(g0 - x))^2 - 1) - x ] - 2 zeta w x',   xp = g0/2,   0 <= x <= xmax
  points open when x rises through xb, close when it falls through xb (inelastic stops at 0 and xmax)
  e(t) = V (battery)  or  k w sin(N (theta(t) - phi) pi/180) (magneto, theta = theta0 + 6 rpm t)

Cases: Ford/K-W coil (Ford 1916 data via Boggess & Patterson 1999, engines/model-t.ts) with the
vibrator parameters fitted here (`--fit`, scipy least_squares) to the DC firing times of Cool386 / ECCT
(3.5 ms at 6 V, 2.5 ms at 9 V, 2.0 ms at 12 V, points re-closing 1.8 ms after the fire at 6 V):
  bench6/9/12  DC bench (battery + 0.05 ohm, no timer), glow gap at a fixed 10 kV breakdown voltage
               (bench air gap): trips, firing currents, re-closes, breakdowns, extinctions, energies;
  open6/12     6 / 12 V, gap never breaks down: open-circuit ring-up after the first trip (peak V1, V2);
  magneto1000  magneto at 1000 rpm through the timer (0.1 ohm), make at -25 deg, break after 87 deg
               (ideal break with the secondary flux linkage conserved), engine-like gap (8 kV);
  slow150      magneto at 150 rpm (5 V peak, 20 Hz), make at a rising EMF zero: firing current on slow
               pulses (cf. 3.0-4.4 A on the hand-cranked coil tester, Kossor 2017);
  ladder600    first-spark angle vs timer make at 600 rpm on the magneto (Patterson & Coniff 2003).

Writes test/fixtures/ignition_trembler.json. `--fit` re-runs the vibrator fit (slow, ~5 min) and
prints the result; the fitted values are the constants below.
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import rnd, write_fixture  # noqa: E402

IN = 0.0254
COIL = dict(L1=3.3e-3, R1=0.295, L2=22.0, R2=3300.0, C2=40e-12, k=0.9, C1=0.43e-6)
# fitted vibrator (see --fit); breakTravel / maxTravel sourced (cushion 0.005 in, point gap 1/32 in)
VIB = dict(Ip=3.53, fn=133.0, zeta=0.34, g0=1.16e-3, xb=0.005 * IN, xmax=0.005 * IN + IN / 32)
VSHEATH = 3 * 365.0 / 14.6 * math.log(1 + 1 / 0.02)
BENCH_GAP = dict(Vbd=10e3, p=1e5, l=3e-3)
ENGINE_GAP = dict(Vbd=8e3, p=5e5, l=IN / 32)
GLOW = dict(Vsheath=VSHEATH, Ccol=40.46, nI=-0.32, npres=0.51, Iext=2e-3)
MAGNETO = dict(k=0.317, N=8, phi=-8.5, Rs=0.3, Ls=3.05e-3)
TIMER_R = 0.1
BATTERY_R = 0.05
RTOL, ATOL = 1e-9, 1e-12


class Trembler:
    def __init__(self, P, source):
        """P: coil + vibrator + gap + supply (Ls, Rs, Rt); source ('dc', V) or ('magneto', k, N, phi, rpm, theta0)."""
        self.P = P
        self.source = source
        self.M = P["k"] * math.sqrt(P["L1"] * P["L2"])
        self.Lp = P["L1"] + P["Ls"]
        self.R = P["R1"] + P["Rs"] + P["Rt"]
        det = self.Lp * P["L2"] - self.M ** 2
        self.i11, self.i12, self.i22 = P["L2"] / det, -self.M / det, self.Lp / det
        self.w = 2 * math.pi * P["fn"]
        self.xp = P["g0"] / 2

    def emf(self, t):
        s = self.source
        if s[0] == "dc":
            return s[1]
        _, k, N, phi, rpm, th0 = s
        w = 2 * math.pi * rpm / 60
        return k * w * math.sin(N * (th0 + 6 * rpm * t - phi) * math.pi / 180)

    def vg(self, i):
        P = self.P
        i = max(i, 1e-4)
        return P["Vsheath"] + P["Ccol"] * (P["l"] * 1e3) * (P["p"] * 1e-5) ** P["npres"] * i ** P["nI"]

    def accel(self, x, v, I):
        P, w = self.P, self.w
        g = P["g0"] / (P["g0"] - x)
        return w * w * (self.xp * ((I / P["Ip"]) ** 2 * g * g - 1) - x) - 2 * P["zeta"] * w * v

    def stored(self, y):
        I1, I2, V1, V2 = y[0], y[1], y[2], y[3]
        P = self.P
        return (0.5 * self.Lp * I1 ** 2 + self.M * I1 * I2 + 0.5 * P["L2"] * I2 ** 2
                + 0.5 * P["C1"] * V1 ** 2 + 0.5 * P["C2"] * V2 ** 2)

    def run(self, t_end, t_break=math.inf, gap=True, max_events=100000, stop_after_bd=False, samples=False):
        """Timer closed on [0, t_break). Returns event times and energies."""
        P = self.P
        y = [0.0] * 9  # I1 I2 V1 V2 x v Esup ER ER2
        t, timer, pts, smode, vib, sgn = 0.0, "closed", "closed", "open", "stop0", 1.0
        out = dict(trips=[], tripI=[], closes=[], bds=[], bdI2=[], exts=[], peakV1_first=0.0, peakV2_first=0.0,
                   tpeakV2_first=0.0)
        E_pts = E_timer = 0.0
        W0 = 0.0
        ts_all, ys_all = [], []
        n = 0
        while t < t_end and n < max_events:
            n += 1
            t_stop = min(t_end, t_break) if timer == "closed" else t_end

            def f(tt, yy, timer=timer, pts=pts, smode=smode, vib=vib, sgn=sgn):
                I1, I2, V1, V2, x, v = yy[0], yy[1], yy[2], yy[3], yy[4], yy[5]
                if smode == "cond":
                    V2 = -sgn * self.vg(abs(I2))
                e = self.emf(tt)
                if timer == "closed":
                    a = e - self.R * I1 - (V1 if pts == "open" else 0.0)
                    b = V2 - P["R2"] * I2
                    dI1, dI2 = self.i11 * a + self.i12 * b, self.i12 * a + self.i22 * b
                else:
                    dI1, dI2 = 0.0, (V2 - P["R2"] * I2) / P["L2"]
                dV1 = I1 / P["C1"] if pts == "open" else 0.0
                dV2 = -I2 / P["C2"] if smode == "open" else 0.0
                dx, dv = (v, self.accel(x, v, I1)) if vib == "free" else (0.0, 0.0)
                return [dI1, dI2, dV1, dV2, dx, dv, e * I1, self.R * I1 * I1, P["R2"] * I2 * I2]

            events = []

            def ev_pts(tt, yy):
                return yy[4] - P["xb"]
            ev_pts.terminal, ev_pts.direction = True, (1 if pts == "closed" else -1)
            events.append(ev_pts)
            if vib == "free":
                def ev_stop0(tt, yy):
                    return yy[4]
                ev_stop0.terminal, ev_stop0.direction = True, -1

                def ev_stopm(tt, yy):
                    return yy[4] - P["xmax"]
                ev_stopm.terminal, ev_stopm.direction = True, 1
                events += [ev_stop0, ev_stopm]
            elif vib == "stop0":
                def ev_rel0(tt, yy):
                    return abs(yy[0]) - P["Ip"]
                ev_rel0.terminal, ev_rel0.direction = True, 1
                events.append(ev_rel0)
            else:
                def ev_relm(tt, yy):
                    return self.accel(P["xmax"], 0.0, yy[0])
                ev_relm.terminal, ev_relm.direction = True, -1
                events.append(ev_relm)
            if gap and smode == "open":
                def ev_bd(tt, yy):
                    return abs(yy[3]) - P["Vbd"]
                ev_bd.terminal, ev_bd.direction = True, 1
                events.append(ev_bd)
            if smode == "cond":
                def ev_ext(tt, yy):
                    return abs(yy[1]) - P["Iext"]
                ev_ext.terminal, ev_ext.direction = True, -1
                events.append(ev_ext)
            max_step = 1e-6 if (pts == "open" or smode == "cond") else 10e-6
            sol = solve_ivp(f, (t, t_stop), y, method="Radau", rtol=RTOL, atol=ATOL, events=events,
                            max_step=max_step)
            assert sol.status >= 0, sol.message
            if len(out["trips"]) == 1 and not out["bds"]:
                out["peakV1_first"] = max(out["peakV1_first"], float(np.max(np.abs(sol.y[2]))))
                k = int(np.argmax(np.abs(sol.y[3])))
                if abs(sol.y[3][k]) > out["peakV2_first"]:
                    out["peakV2_first"] = float(abs(sol.y[3][k]))
                    out["tpeakV2_first"] = float(sol.t[k] - out["trips"][0])
            if samples:
                ts_all.append(sol.t)
                ys_all.append(sol.y)
            t = sol.t[-1]
            y = list(sol.y[:, -1])
            if sol.status == 0:
                if t >= t_end - 1e-15:
                    break
                if timer == "closed" and t >= t_break - 1e-15:
                    w0 = self.stored(y)
                    y[1] += self.M * y[0] / P["L2"]
                    y[0] = 0.0
                    E_timer += w0 - self.stored(y)
                    timer = "open"
                    if smode == "cond" and y[1] * sgn <= P["Iext"]:
                        smode = "open"
                        out["exts"].append(t)
                continue
            name = events[[i for i, te in enumerate(sol.t_events) if len(te) > 0][0]].__name__
            if name == "ev_pts":
                if pts == "closed":
                    pts = "open"
                    if timer == "closed":
                        out["trips"].append(t)
                        out["tripI"].append(abs(y[0]))
                else:
                    pts = "closed"
                    E_pts += 0.5 * P["C1"] * y[2] ** 2
                    y[2] = 0.0
                    out["closes"].append(t)
                    # ledger at the re-close (gap energy by balance; the gap is not conducting here)
                    W = self.stored(y)
                    out.setdefault("atClose", []).append(dict(
                        E_sup=y[6], E_R=y[7], E_R2=y[8], E_pts=E_pts, E_gap=y[6] - y[7] - y[8] - E_pts - E_timer - W, W=W))
                y[4] = P["xb"]
            elif name == "ev_stop0":
                vib, y[4], y[5] = "stop0", 0.0, 0.0
            elif name == "ev_stopm":
                vib, y[4], y[5] = "stopm", P["xmax"], 0.0
            elif name in ("ev_rel0", "ev_relm"):
                vib = "free"
            elif name == "ev_bd":
                sgn = 1.0 if y[1] > 0 else -1.0
                y[3] = -sgn * self.vg(abs(y[1]))
                smode = "cond"
                out["bds"].append(t)
                out["bdI2"].append(abs(y[1]))
                if stop_after_bd:
                    break
            elif name == "ev_ext":
                y[3] = -sgn * self.vg(abs(y[1]))
                smode = "open"
                out["exts"].append(t)
        W1 = self.stored(y)
        out.update(E_sup=y[6], E_R=y[7], E_R2=y[8], E_pts=E_pts, E_timer=E_timer, W_end=W1, t_end=t,
                   y_end=list(y[:6]))
        # gap energy from the energy balance (the dumps and the algebraic C2 energy are implied)
        out["E_gap"] = y[6] - y[7] - y[8] - E_pts - E_timer - (W1 - W0)
        if samples:
            out["t"] = np.concatenate(ts_all)
            out["Y"] = np.concatenate(ys_all, axis=1)
        return out


def params(gap, Ls=0.0, Rs=BATTERY_R, Rt=0.0, **over):
    P = dict(COIL, **VIB, **GLOW, **gap, Ls=Ls, Rs=Rs, Rt=Rt)
    P.update(over)
    return P


def ms_list(v):
    return [rnd(float(x)) for x in v]


def bench(V, t_end):
    o = Trembler(params(BENCH_GAP), ("dc", V)).run(t_end)
    first_close = {k: rnd(v) for k, v in o["atClose"][0].items()}
    return dict(V=V, t_end=rnd(t_end), trips=ms_list(o["trips"]), tripI=ms_list(o["tripI"]), closes=ms_list(o["closes"]),
                bds=ms_list(o["bds"]), bdI2=ms_list(o["bdI2"]), exts=ms_list(o["exts"]), atFirstClose=first_close,
                E_sup=rnd(o["E_sup"]), E_R=rnd(o["E_R"]), E_R2=rnd(o["E_R2"]), E_pts=rnd(o["E_pts"]),
                E_gap=rnd(o["E_gap"]), W_end=rnd(o["W_end"]))


def open_circuit(V):
    o = Trembler(params(BENCH_GAP), ("dc", V)).run(4.5e-3, gap=False)
    return dict(V=V, trip=rnd(o["trips"][0]), tripI=rnd(o["tripI"][0]), peakV1=rnd(o["peakV1_first"]),
                peakV2=rnd(o["peakV2_first"]), tpeakV2=rnd(o["tpeakV2_first"]))


def magneto_case(rpm, make, t_after_break=1e-3, gap=ENGINE_GAP, t_end=None, **over):
    m = MAGNETO
    t_break = 87 / (6 * rpm)
    t_end = t_end if t_end is not None else t_break + t_after_break
    P = params(gap, Ls=m["Ls"], Rs=m["Rs"], Rt=TIMER_R, **over)
    o = Trembler(P, ("magneto", m["k"], m["N"], m["phi"], rpm, make)).run(t_end, t_break=t_break)
    deg = lambda ts: [rnd(make + 6 * rpm * float(x)) for x in ts]
    return dict(rpm=rpm, make=make, t_break=rnd(t_break), t_end=rnd(t_end),
                trips=ms_list(o["trips"]), tripI=ms_list(o["tripI"]), closes=ms_list(o["closes"]),
                bds=ms_list(o["bds"]), exts=ms_list(o["exts"]), firstSparkDeg=(deg(o["bds"][:1]) or [None])[0],
                E_sup=rnd(o["E_sup"]), E_R=rnd(o["E_R"]), E_R2=rnd(o["E_R2"]), E_pts=rnd(o["E_pts"]),
                E_timer=rnd(o["E_timer"]), E_gap=rnd(o["E_gap"]), W_end=rnd(o["W_end"]))


def ladder(rpm=600, makes=None):
    m = MAGNETO
    makes = makes if makes is not None else [15.5 - 2.5 * i for i in range(33)]
    P = params(ENGINE_GAP, Ls=m["Ls"], Rs=m["Rs"], Rt=TIMER_R)
    out = []
    for mk in makes:
        o = Trembler(P, ("magneto", m["k"], m["N"], m["phi"], rpm, mk)).run(87 / (6 * rpm), stop_after_bd=True)
        out.append(rnd(mk + 6 * rpm * o["bds"][0]) if o["bds"] else None)
    return dict(rpm=rpm, makes=[rnd(x) for x in makes], firstSparkDeg=out)


# ---------------------------------------------------------------- vibrator fit
TARGETS = dict(t6=3.5e-3, t9=2.5e-3, t12=2.0e-3, reclose6=1.8e-3)


def fit_residuals(q, verbose=False):
    Ip, fn, zeta, g0mm = q
    P = params(BENCH_GAP, Ip=Ip, fn=fn, zeta=zeta, g0=g0mm * 1e-3)
    o = {V: Trembler(P, ("dc", V)).run(te, max_events=60) for V, te in ((6, 7.5e-3), (9, 4e-3), (12, 3.5e-3))}
    t = {V: (o[V]["trips"][0] if o[V]["trips"] else 2 * te) for V, te in ((6, 4e-3), (9, 2.5e-3), (12, 2e-3))}
    rc = (o[6]["closes"][0] - t[6]) if o[6]["closes"] else 5e-3
    r = [(t[6] - TARGETS["t6"]) / 1e-4, (t[9] - TARGETS["t9"]) / 1e-4, (t[12] - TARGETS["t12"]) / 1e-4,
         (rc - TARGETS["reclose6"]) / 1e-4]
    if verbose:
        print([round(float(a), 4) for a in q], [round(float(x), 4) for x in r], flush=True)
    return r


def run_fit():
    sol = least_squares(lambda q: fit_residuals(q, True), [3.0, 150.0, 0.1, 1.2],
                        bounds=([1.0, 30, 0.0, 0.95], [6.0, 1000, 1.0, 10.0]), diff_step=1e-3, max_nfev=60)
    print("fit:", sol.x, "residuals (0.1 ms units):", sol.fun)


def main():
    if "--fit" in sys.argv:
        run_fit()
        return
    b6, b9, b12 = bench(6, 12e-3), bench(9, 8e-3), bench(12, 8e-3)
    oc6, oc12 = open_circuit(6), open_circuit(12)
    mag = magneto_case(1000, -25.0)
    slow = magneto_case(150, MAGNETO["phi"], t_end=20e-3)
    lad = ladder()
    fit = dict(targets={k: rnd(v) for k, v in TARGETS.items()},
               achieved=dict(t6=b6["trips"][0], t9=b9["trips"][0], t12=b12["trips"][0],
                             reclose6=rnd(b6["closes"][0] - b6["trips"][0]),
                             I6=b6["tripI"][0], I9=b9["tripI"][0], I12=b12["tripI"][0]))
    write_fixture("ignition_trembler.json", dict(
        params=dict(coil={k: rnd(v) for k, v in COIL.items()}, vibrator={k: rnd(v) for k, v in VIB.items()},
                    glow={k: rnd(v) for k, v in GLOW.items()}, benchGap={k: rnd(v) for k, v in BENCH_GAP.items()},
                    engineGap={k: rnd(v) for k, v in ENGINE_GAP.items()}, magneto={k: rnd(v) for k, v in MAGNETO.items()},
                    timerResistance=TIMER_R, batteryResistance=BATTERY_R),
        bench6=b6, bench9=b9, bench12=b12, open6=oc6, open12=oc12, magneto1000=mag, slow150=slow, ladder600=lad, fit=fit,
    ))
    for name, c in (("bench6", b6), ("bench9", b9), ("bench12", b12)):
        print(name, "trips", [round(x * 1e3, 4) for x in c["trips"]], "I", [round(x, 3) for x in c["tripI"]],
              "closes", [round(x * 1e3, 4) for x in c["closes"]], "bds", len(c["bds"]), "E_gap", round(c["E_gap"] * 1e3, 3), "mJ")
    print("open6", oc6, "\nopen12", oc12)
    print("magneto1000 trips", [round(x * 1e3, 4) for x in mag["trips"]], "I", [round(x, 3) for x in mag["tripI"]],
          "bds", len(mag["bds"]), "first spark", mag["firstSparkDeg"], "E_timer", mag["E_timer"])
    print("slow150 trips", [round(x * 1e3, 4) for x in slow["trips"]], "I", [round(x, 3) for x in slow["tripI"]],
          "first spark", slow["firstSparkDeg"])
    print("ladder600", list(zip(lad["makes"], lad["firstSparkDeg"])))


if __name__ == "__main__":
    main()
