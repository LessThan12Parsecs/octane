"""Slider-crank kinematics oracle (independent of the TypeScript closed forms).

Method (deliberately different from src/physics/mechanics/kinematics.ts):
  * the mechanism position is found by NUMERICALLY solving the loop-closure
    equation for the connecting-rod angle (scipy brentq to ~1e-15), not by the
    closed-form sqrt expression;
  * TDC / BDC are found by numerical maximisation / minimisation of the pin
    height (golden-section to ~1e-13 rad), not from asin(e/(l±a));
  * all derivatives (dV/dθ, d²x/dθ², dβ/dθ, d²β/dθ²) are 6th-order central
    finite differences of that numerical solution;
  * for zero pin offset the volume is also evaluated with Heywood (1988)
    eq. 2.6 as an additional textbook reference.

Frame/sign conventions (must match kinematics.ts):
  crank axis at origin, cylinder axis = +y, crank-pin angle φ measured from +y
  toward +x in the direction of rotation; crank pin at (a sinφ, a cosφ);
  the wrist pin moves on the line x = x_w with x_w = +pinOffset
  (the MAJOR-thrust wall is -x because the crank pin is at +x during the
  expansion stroke, so an offset toward the major-thrust side is negative).
  θ = φ - φ_TDC (θ = 0 at firing TDC). β = rod angle from the cylinder axis,
  positive when the crank pin is on the +x side of the wrist pin.

Writes test/fixtures/mechanics_kinematics.json.
"""
from __future__ import annotations

import json
import math
import pathlib

import numpy as np
from scipy.optimize import brentq

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "mechanics_kinematics.json"

IN = 0.0254

CASES = [
    # CFR F-1: 3.25 x 4.5 in, 254 mm rod (Pal et al. 2018 Table 2), no offset.
    dict(name="cfr_cr7", bore=3.25 * IN, stroke=4.5 * IN, rod=0.254, pinOffset=0.0, crevice=1.0e-6, cr=7.0),
    dict(name="cfr_cr18", bore=3.25 * IN, stroke=4.5 * IN, rod=0.254, pinOffset=0.0, crevice=1.0e-6, cr=18.0),
    # Synthetic automotive-like geometry with a large pin offset (exercises offset terms).
    dict(name="offset_pos", bore=0.086, stroke=0.086, rod=0.143, pinOffset=0.004, crevice=0.5e-6, cr=10.5),
    dict(name="offset_neg_shortrod", bore=0.080, stroke=0.090, rod=0.135, pinOffset=-0.006, crevice=0.0, cr=9.0),
    # Stress cases (reviewer): offset at 90 % of the geometric limit l − a with a very short rod
    # (l/a = 1.6), a tiny model-engine geometry and a large marine-diesel-like geometry at high CR.
    dict(name="extreme_offset_shortrod", bore=0.070, stroke=0.100, rod=0.080, pinOffset=0.027, crevice=0.2e-6, cr=6.0),
    dict(name="tiny_engine", bore=0.012, stroke=0.010, rod=0.020, pinOffset=0.0005, crevice=0.0, cr=15.0),
    dict(name="huge_engine", bore=0.98, stroke=2.4, rod=3.6, pinOffset=0.0, crevice=2.0e-4, cr=22.0),
]

THETA_DEG = np.arange(-360.0, 360.0, 15.0).tolist() + [-359.3, -181.1, -0.4, 0.4, 13.7, 179.9, 181.2, 359.6]


def rod_angle_numeric(phi: float, a: float, l: float, xw: float) -> float:
    """Solve crank_pin_x + l*(-sin β) = xw for β in (-π/2, π/2)."""
    f = lambda b: a * math.sin(phi) - l * math.sin(b) - xw
    return brentq(f, -math.pi / 2 + 1e-12, math.pi / 2 - 1e-12, xtol=1e-16, rtol=1e-15, maxiter=500)


def pin_height(phi: float, a: float, l: float, xw: float) -> float:
    b = rod_angle_numeric(phi, a, l, xw)
    return a * math.cos(phi) + l * math.cos(b)


def fd1(f, x, h):
    # 6th-order central first derivative
    return (-f(x - 3 * h) + 9 * f(x - 2 * h) - 45 * f(x - h) + 45 * f(x + h) - 9 * f(x + 2 * h) + f(x + 3 * h)) / (60 * h)


def fd2(f, x, h):
    # 6th-order central second derivative
    return (2 * f(x - 3 * h) - 27 * f(x - 2 * h) + 270 * f(x - h) - 490 * f(x) + 270 * f(x + h) - 27 * f(x + 2 * h) + 2 * f(x + 3 * h)) / (180 * h * h)


def build(case: dict) -> dict:
    B, S, l, e, Vcr, cr = case["bore"], case["stroke"], case["rod"], case["pinOffset"], case["crevice"], case["cr"]
    a = S / 2
    xw = e
    A = math.pi * B * B / 4
    yp = lambda phi: pin_height(phi, a, l, xw)
    # Extrema of the pin height: roots of the (finite-difference) derivative. A root of
    # dy/dφ is located to ~1e-15 rad, whereas maximising y itself is limited to ~sqrt(eps).
    dy = lambda p: fd1(yp, p, 1e-3)
    # bracket the extrema from a coarse grid search (no closed-form knowledge), then refine the root of dy/dφ
    grid = np.linspace(-math.pi / 2, 3 * math.pi / 2, 721)
    ys = np.array([yp(p) for p in grid])
    i_max, i_min = int(np.argmax(ys)), int(np.argmin(ys))
    dg = grid[1] - grid[0]
    phi_tdc = brentq(dy, grid[i_max] - 2 * dg, grid[i_max] + 2 * dg, xtol=1e-16, rtol=1e-15)
    phi_bdc = brentq(dy, grid[i_min] - 2 * dg, grid[i_min] + 2 * dg, xtol=1e-16, rtol=1e-15)
    y_tdc = yp(phi_tdc)
    y_bdc = yp(phi_bdc)
    travel = y_tdc - y_bdc
    Vd = A * travel
    Vc = Vd / (cr - 1)
    h_tdc = (Vc - Vcr) / A

    x_of = lambda th: y_tdc - yp(th + phi_tdc)
    beta_of = lambda th: rod_angle_numeric(th + phi_tdc, a, l, xw)
    V_of = lambda th: Vc + A * x_of(th)
    hstep = 2e-3

    rows = []
    for tdeg in THETA_DEG:
        th = math.radians(tdeg)
        x = x_of(th)
        V = V_of(th)
        row = dict(
            thetaDeg=tdeg,
            x=x,
            h=h_tdc + x,
            V=V,
            dVdTheta=fd1(V_of, th, hstep),
            d2xdTheta2=fd2(x_of, th, hstep),
            dxdTheta=fd1(x_of, th, hstep),
            beta=beta_of(th),
            dbetadTheta=fd1(beta_of, th, hstep),
            d2betadTheta2=fd2(beta_of, th, hstep),
            linerArea=math.pi * B * (h_tdc + x),
        )
        if e == 0.0:
            R = l / a
            # Heywood (1988) eq. 2.6: V/Vc = 1 + ½(rc − 1)[R + 1 − cosθ − (R² − sin²θ)^½]
            row["V_heywood"] = Vc * (1 + 0.5 * (cr - 1) * (R + 1 - math.cos(th) - math.sqrt(R * R - math.sin(th) ** 2)))
        rows.append(row)

    return dict(
        name=case["name"],
        input=dict(bore=B, stroke=S, conRodLength=l, pinOffset=e, creviceVolume=Vcr, compressionRatio=cr),
        tdcAxisAngle=phi_tdc,
        bdcAxisAngle=phi_bdc,
        bdcTheta=phi_bdc - phi_tdc,
        pistonTravel=travel,
        displacedVolume=Vd,
        clearanceVolume=Vc,
        clearanceHeightTDC=h_tdc,
        rows=rows,
    )


def round_floats(o):
    """Round to 15 significant digits (compact, deterministic fixture)."""
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
            "Slider-crank kinematics oracle: numerical loop-closure solution (brentq), numerical TDC/BDC search, "
            "6th-order finite-difference derivatives; Heywood (1988) eq. 2.6 volume for zero offset. SI units, "
            "theta in crank degrees with 0 = firing TDC; x = crown displacement below TDC; beta = rod angle."
        ),
        generator="tools/reference/mechanics_kinematics.py",
        cases=[build(c) for c in CASES],
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(round_floats(out), separators=(",", ":")) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
