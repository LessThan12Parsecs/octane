"""Oracle for src/physics/chemistry/knock.ts (cylinder acoustic modes).

Writes test/fixtures/chemistry_knock_modes.json:
  jnpZeros[m] : first 6 positive zeros of J'_m, m = 0..8 (scipy.special.jnp_zeros; for m = 0
                scipy omits x = 0, matching our convention)
  besselJ     : J_n(x) samples (scipy.special.jv) for n = 0..8, x in (0, 60]
  projection  : modal source projections <psi>_eg / <psi^2>_disc for an end-gas region outside
                a flame circle, by adaptive 2-D quadrature (scipy.integrate.dblquad in polar
                coordinates) -- checks the TS midpoint quadrature.
"""
from __future__ import annotations

import math
import os

import numpy as np
from scipy import integrate, special

from chemistry_common import FIXTURES, write_json


def main():
    zeros = {str(m): [float(z) for z in special.jnp_zeros(m, 6)] for m in range(9)}
    xs = [0.01, 0.5, 1.0, 1.8412, 3.0, 5.0, 7.5, 10.0, 15.0, 22.0, 35.0, 60.0]
    bj = [[float(special.jv(n, x)) for x in xs] for n in range(9)]

    # projection check: R = 1, flame circle centre (x0, z0) radius rf
    R = 1.0
    x0, z0, rf = -0.55, 0.1, 0.9
    modes = [(1, 0, 1), (1, 0, 0), (2, 0, 1), (2, 0, 0), (0, 1, 1), (3, 0, 1), (1, 1, 1)]
    proj = []
    for (m, n, cosine) in modes:
        a = float(special.jnp_zeros(m, n + 1)[n]) if m > 0 else float(special.jnp_zeros(0, n)[n - 1])

        def psi(r, th):
            trig = math.cos(m * th) if cosine else math.sin(m * th)
            return special.jv(m, a * r / R) * trig

        def eg_intervals(th):
            """r-intervals of [0, R] outside the flame circle along the ray at angle th."""
            b = x0 * math.cos(th) + z0 * math.sin(th)
            c = x0 * x0 + z0 * z0 - rf * rf
            disc = b * b - c
            if disc <= 0:
                return [(0.0, R)]
            r1, r2 = b - math.sqrt(disc), b + math.sqrt(disc)
            out = []
            if r1 > 0:
                out.append((0.0, min(r1, R)))
            if r2 < R:
                out.append((max(r2, 0.0), R))
            return [(u, v) for u, v in out if v > u]

        def ray(th, weight):
            tot = 0.0
            for u, v in eg_intervals(th):
                tot += integrate.quad(lambda r: weight(r, th) * r, u, v, epsabs=1e-13, epsrel=1e-12)[0]
            return tot

        # angular breakpoints where the ray becomes tangent to / enters the flame circle
        d0 = math.hypot(x0, z0)
        phi0 = math.atan2(z0, x0)
        pts = [phi0 % (2 * math.pi)]
        if d0 > rf:
            half = math.asin(rf / d0)
            pts += [(phi0 + half) % (2 * math.pi), (phi0 - half) % (2 * math.pi)]
        opts = dict(epsabs=1e-12, epsrel=1e-11, limit=400, points=sorted(pts))
        num = integrate.quad(lambda th: ray(th, psi), 0, 2 * math.pi, **opts)[0]
        area = integrate.quad(lambda th: ray(th, lambda r, t: 1.0), 0, 2 * math.pi, **opts)[0]
        ms = (1 - m * m / (a * a)) * special.jv(m, a) ** 2 * (1.0 if m == 0 else 0.5)
        proj.append(dict(m=m, n=n, cosine=cosine, alpha=a, shape=(num / area) / ms))
        print(m, n, cosine, a, (num / area) / ms)
    write_json(os.path.join(FIXTURES, "chemistry_knock_modes.json"), dict(
        jnpZeros=zeros, besselX=xs, besselJ=bj,
        projection=dict(R=R, x0=x0, z0=z0, rf=rf, endGasAreaFraction=area / (math.pi * R * R),
                        modes=proj)))


if __name__ == "__main__":
    main()
