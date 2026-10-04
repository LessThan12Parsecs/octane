"""Crank-train dynamics oracle (independent numerical Lagrangian).

Method (deliberately different from src/physics/mechanics/dynamics.ts):
  * positions of the wrist pin and the rod CG come from the NUMERICAL loop-closure
    solution (brentq on the rod angle), not closed forms;
  * the reflected mechanism inertia J_m(θ) = 2·KE/ω² is obtained from
    finite-difference velocities of those positions (piston translation, rod CG
    translation and rod rotation);
  * the inertia torque follows from Lagrange's equation,
        T_inertia = −(J_m α + ½ dJ_m/dθ ω²),
    with dJ_m/dθ by finite differences of J_m;
  * the gravity torque is −dU/dθ with U the potential energy of piston + rod
    (finite differences).

Conventions identical to tools/reference/mechanics_kinematics.py
(x_w = +pinOffset, θ = 0 at the exact TDC found numerically).

Multi-cylinder cases ("multiCases"; MultiCylinderCrankTrain): N identical mechanisms on one rigid
crank, cylinder i at θ − radians(firingOffsetDeg[i]) (the raw 720° offset: the 2π-periodicity that
the TypeScript reduces to a throw phase is NOT assumed here). Totals J = Σ J_i, dJ/dθ by finite
differences of the TOTAL J, gravity −dU_total/dθ, per-cylinder gas torque (p_i − p_cc)·A·dx/dθ with
dx/dθ = −dy_p/dθ from finite differences of the numerical wrist-pin height, and the free
acceleration α = (Σ T_gas + T_grav − ½ dJ/dθ ω² + T_ext)/(J_rot + J_load + J).

Writes test/fixtures/mechanics_dynamics.json.
"""
from __future__ import annotations

import json
import math
import pathlib

from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "mechanics_dynamics.json"

IN = 0.0254
G = 9.80665

CASES = [
    dict(
        name="cfr_like",
        geometry=dict(bore=3.25 * IN, stroke=4.5 * IN, conRodLength=0.254, pinOffset=0.0, creviceVolume=1.0e-6),
        masses=dict(piston=1.6, conRod=1.9, conRodCgFromBigEnd=0.075, conRodInertiaCg=0.012, rotatingInertia=2.0),
        states=[dict(omega=2 * math.pi * 10, alpha=0.0), dict(omega=2 * math.pi * 15, alpha=150.0)],
    ),
    dict(
        name="offset",
        geometry=dict(bore=0.086, stroke=0.086, conRodLength=0.143, pinOffset=0.004, creviceVolume=0.5e-6),
        masses=dict(piston=0.45, conRod=0.55, conRodCgFromBigEnd=0.04, conRodInertiaCg=0.0018, rotatingInertia=0.25),
        states=[dict(omega=2 * math.pi * 50, alpha=-800.0)],
    ),
    # Stress case (reviewer): offset at 90 % of l − a, very short rod (l/a = 1.6), heavy rod with its CG
    # near the small end and a rod inertia well above the two-mass equivalent; reversed rotation (ω < 0).
    dict(
        name="extreme_offset_shortrod",
        geometry=dict(bore=0.070, stroke=0.100, conRodLength=0.080, pinOffset=0.027, creviceVolume=0.2e-6),
        masses=dict(piston=0.9, conRod=1.2, conRodCgFromBigEnd=0.06, conRodInertiaCg=0.004, rotatingInertia=0.1),
        states=[dict(omega=2 * math.pi * 40, alpha=2000.0), dict(omega=-2 * math.pi * 25, alpha=-300.0)],
    ),
]

THETA_DEG = [-360.0 + 20.0 * i for i in range(36)] + [-181.0, -3.3, 7.7, 91.0, 263.4]


def fd1(f, x, h):
    return (-f(x - 3 * h) + 9 * f(x - 2 * h) - 45 * f(x - h) + 45 * f(x + h) - 9 * f(x + 2 * h) + f(x + 3 * h)) / (60 * h)


class Mechanism:
    def __init__(self, geometry: dict, masses: dict):
        self.a = geometry["stroke"] / 2
        self.l = geometry["conRodLength"]
        self.xw = geometry["pinOffset"]
        self.m = masses
        yp = lambda p: self.positions(p)[1]
        dy = lambda p: fd1(yp, p, 1e-3)
        grid = [-math.pi / 2 + i * (2 * math.pi / 720) for i in range(721)]
        i_max = max(range(len(grid)), key=lambda i: yp(grid[i]))
        dg = grid[1] - grid[0]
        self.phi_tdc = brentq(dy, grid[i_max] - 2 * dg, grid[i_max] + 2 * dg, xtol=1e-16, rtol=1e-15)

    def rod_angle(self, phi: float) -> float:
        f = lambda b: self.a * math.sin(phi) - self.l * math.sin(b) - self.xw
        return brentq(f, -math.pi / 2 + 1e-12, math.pi / 2 - 1e-12, xtol=1e-16, rtol=1e-15, maxiter=500)

    def positions(self, phi: float):
        """(x_w, y_p, x_g, y_g, beta) at axis angle phi."""
        b = self.rod_angle(phi)
        cx, cy = self.a * math.sin(phi), self.a * math.cos(phi)
        wx, wy = cx - self.l * math.sin(b), cy + self.l * math.cos(b)
        lam = self.m["conRodCgFromBigEnd"] / self.l
        gx, gy = (1 - lam) * cx + lam * wx, (1 - lam) * cy + lam * wy
        return wx, wy, gx, gy, b

    def J(self, theta: float) -> float:
        phi = theta + self.phi_tdc
        h = 1e-3
        comp = lambda i: (lambda p: self.positions(p)[i])
        ypd = fd1(comp(1), phi, h)
        gxd = fd1(comp(2), phi, h)
        gyd = fd1(comp(3), phi, h)
        bd = fd1(comp(4), phi, h)
        m = self.m
        return m["piston"] * ypd ** 2 + m["conRod"] * (gxd ** 2 + gyd ** 2) + m["conRodInertiaCg"] * bd ** 2

    def U(self, theta: float) -> float:
        _, wy, _, gy, _ = self.positions(theta + self.phi_tdc)
        return G * (self.m["piston"] * wy + self.m["conRod"] * gy)


def build(case: dict) -> dict:
    mech = Mechanism(case["geometry"], case["masses"])
    rows = []
    for tdeg in THETA_DEG:
        th = math.radians(tdeg)
        J = mech.J(th)
        dJ = fd1(mech.J, th, 2e-3)
        Tg = -fd1(mech.U, th, 1e-3)
        row = dict(thetaDeg=tdeg, J=J, dJ=dJ, gravityTorque=Tg, inertiaTorque=[])
        for st in case["states"]:
            row["inertiaTorque"].append(-(J * st["alpha"] + 0.5 * dJ * st["omega"] ** 2))
        rows.append(row)
    return dict(name=case["name"], geometry=case["geometry"], masses=case["masses"], gravity=G, states=case["states"], rows=rows)


MULTI_CASES = [
    # Ford Model T-like inline four, firing order 1-2-4-3 (offsets 0/180/540/360 -> throws 0/180/180/0),
    # masses as engines/model-t.ts; one state with a vehicle-like load inertia.
    dict(
        name="inline4_1243_modelt_like",
        geometry=dict(bore=3.75 * IN, stroke=4.0 * IN, conRodLength=7.0 * IN, pinOffset=0.0, creviceVolume=2.5e-6),
        masses=dict(piston=1.04, conRod=0.6, conRodCgFromBigEnd=7.0 * IN / 3, conRodInertiaCg=0.0042, rotatingInertia=0.65),
        firingOffsetDeg=[0.0, 180.0, 540.0, 360.0],
        pCrankcase=1.0e5,
        states=[
            dict(omega=2 * math.pi * 1000 / 60, alpha=0.0, externalTorque=-20.0, loadInertia=0.0,
                 pressures=[25.0e5, 1.1e5, 0.6e5, 4.0e5]),
            dict(omega=2 * math.pi * 1800 / 60, alpha=-35.0, externalTorque=-60.0, loadInertia=9.8,
                 pressures=[8.0e5, 0.9e5, 1.3e5, 12.0e5]),
        ],
    ),
    # Generic three-cylinder (offsets 0/240/480 -> throws 0/240/120) with a pin offset and a rod that is
    # not dynamically equivalent to two masses; reversed rotation in the second state.
    dict(
        name="three_cyl_offset_pin",
        geometry=dict(bore=0.080, stroke=0.090, conRodLength=0.150, pinOffset=0.004, creviceVolume=0.5e-6),
        masses=dict(piston=0.42, conRod=0.5, conRodCgFromBigEnd=0.045, conRodInertiaCg=0.0016, rotatingInertia=0.12),
        firingOffsetDeg=[0.0, 240.0, 480.0],
        pCrankcase=0.98e5,
        states=[
            dict(omega=2 * math.pi * 3000 / 60, alpha=400.0, externalTorque=-8.0, loadInertia=0.3,
                 pressures=[40.0e5, 2.0e5, 0.5e5]),
            dict(omega=-2 * math.pi * 1200 / 60, alpha=-150.0, externalTorque=5.0, loadInertia=0.0,
                 pressures=[1.2e5, 15.0e5, 0.8e5]),
        ],
    ),
]


def build_multi(case: dict) -> dict:
    mech = Mechanism(case["geometry"], case["masses"])
    g = case["geometry"]
    area = math.pi * g["bore"] ** 2 / 4
    offs = [math.radians(d) for d in case["firingOffsetDeg"]]
    m = case["masses"]

    def j_total(th):
        return sum(mech.J(th - o) for o in offs)

    def u_total(th):
        return sum(mech.U(th - o) for o in offs)

    def dxdtheta(th):
        yp = lambda p: mech.positions(p)[1]
        return -fd1(yp, th + mech.phi_tdc, 1e-3)

    rows = []
    for tdeg in THETA_DEG:
        th = math.radians(tdeg)
        J = j_total(th)
        dJ = fd1(j_total, th, 2e-3)
        Tg = -fd1(u_total, th, 1e-3)
        x1 = [dxdtheta(th - o) for o in offs]
        row = dict(thetaDeg=tdeg, J=J, dJ=dJ, gravityTorque=Tg, dxdTheta=x1, gasTorque=[], inertiaTorque=[], alpha=[])
        for st in case["states"]:
            gas = [(p - case["pCrankcase"]) * area * x for p, x in zip(st["pressures"], x1)]
            row["gasTorque"].append(gas)
            row["inertiaTorque"].append(-(J * st["alpha"] + 0.5 * dJ * st["omega"] ** 2))
            num = sum(gas) + Tg - 0.5 * dJ * st["omega"] ** 2 + st["externalTorque"]
            row["alpha"].append(num / (m["rotatingInertia"] + st["loadInertia"] + J))
        rows.append(row)
    return dict(
        name=case["name"], geometry=g, masses=m, gravity=G, firingOffsetDeg=case["firingOffsetDeg"],
        pCrankcase=case["pCrankcase"], states=case["states"], rows=rows,
    )


def round_floats(o):
    if isinstance(o, float):
        return float(f"{o:.15g}")
    if isinstance(o, dict):
        return {k: round_floats(v) for k, v in o.items()}
    if isinstance(o, list):
        return [round_floats(v) for v in o]
    return o


def main() -> None:
    out = dict(
        description=(
            "Crank-train dynamics oracle: reflected mechanism inertia J_m(theta) from finite-difference velocities of the "
            "numerical loop-closure solution, inertia torque -(J_m*alpha + 0.5*dJ_m/dtheta*omega^2) (Lagrange), gravity "
            "torque -dU/dtheta. SI units; theta in crank degrees, 0 = firing TDC; torques positive in the rotation direction."
        ),
        generator="tools/reference/mechanics_dynamics.py",
        cases=[build(c) for c in CASES],
        multiDescription=(
            "Multi-cylinder crank train on one rigid crank: cylinder i at theta - firingOffsetDeg[i] (raw 720-deg "
            "offsets). Totals J, dJ/dtheta, gravity torque; per-cylinder dx/dtheta and gas torque "
            "(p_i - pCrankcase)*A*dx/dtheta per state; total inertia torque -(J*alpha + 0.5*dJ*omega^2); free "
            "acceleration (sum gas + gravity - 0.5*dJ*omega^2 + externalTorque)/(rotatingInertia + loadInertia + J)."
        ),
        multiCases=[build_multi(c) for c in MULTI_CASES],
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(round_floats(out), separators=(",", ":")) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
