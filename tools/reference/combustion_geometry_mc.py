"""Oracle for src/physics/combustion/flame-geometry.ts.

Sphere (centre c, radius r) ∩ flat-disc chamber {x² + z² ≤ R², −h ≤ y ≤ 0}:
volume, flame-front area (sphere surface inside the chamber), head / piston / liner areas inside the ball.

Two independent references:
  1. Monte-Carlo (numpy, fixed seed) with ≥ 1e7 samples per quantity and a 1-σ binomial error bar.
     Volume: uniform samples in the bounding box of ball ∩ cylinder. Front area: uniform directions on the
     sphere. Head / piston: uniform samples in the square bounding the cut circle on that plane.
     Liner: uniform samples on the wall patch (angle × height) that can lie inside the ball.
  2. scipy.integrate.quad of the slice integrals (disc slices of the ball are circle–circle lenses with the
     bore; written here independently with the classical acos formulas), relative error ~1e-12.
The script cross-checks quad against MC (|Δ| ≤ 5σ) and writes
  test/fixtures/combustion_geometry_mc.json    — MC cases (value, σ) + quad values
  test/fixtures/combustion_geometry_quad.json  — dense deterministic (r, h) samples, quad values only
Run: .venv/bin/python tools/reference/combustion_geometry_mc.py
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy import integrate

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = 0.0254

# (name, bore, centre) — cylinder frame, metres
GEOMETRIES = [
    ("cfr", 3.25 * IN, (0.9 * (3.25 * IN) / 2, -0.002, 0.0)),
    ("offset", 3.25 * IN, (0.3 * (3.25 * IN) / 2, -0.006, 0.4 * (3.25 * IN) / 2)),
    ("central", 0.100, (0.0, -0.005, 0.0)),
    # added by review: the actual CFR_F1 spec gap centre (engines/cfr.ts: [0, −2 mm, R − 1 mm]), a spark
    # 0.4 mm from the liner, and a near-axis spark (d = 1e-6 R) that exercises the small-d numerics
    ("cfr_spec", 3.25 * IN, (0.0, -0.002, (3.25 * IN) / 2 - 0.001)),
    ("wall099", 3.25 * IN, (0.99 * (3.25 * IN) / 2, -0.002, 0.0)),
    ("near_axis", 0.100, (0.05e-6, -0.005, 0.0)),
]


# ------------------------------------------------------------------------------------------------
# quadrature reference
# ------------------------------------------------------------------------------------------------
def _cos_ab(rho: float, R: float, d: float):
    """cos of the chord half-angles at the slice centre (a) and at the axis (b). Law of cosines written
    with P1 = ρ² − (R − d)² so that d² + ρ² − R² = P1 − 2d(R − d) and d² + R² − ρ² = 2dR − P1 carry no
    R² − ρ² cancellation for a near-axis spark (d ≪ R)."""
    P1 = (rho - (R - d)) * (rho + (R - d))
    ca = max(-1.0, min(1.0, (P1 - 2 * d * (R - d)) / (2 * d * rho)))
    cb = max(-1.0, min(1.0, (2 * d * R - P1) / (2 * d * R)))
    return ca, cb


def lens(rho: float, R: float, d: float) -> float:
    if rho <= 0:
        return 0.0
    if rho <= R - d:
        return math.pi * rho * rho
    if rho >= R + d:
        return math.pi * R * R
    ca, cb = _cos_ab(rho, R, d)
    k = (-d + rho + R) * (d + rho - R) * (d - rho + R) * (d + rho + R)
    return rho * rho * math.acos(ca) + R * R * math.acos(cb) - 0.5 * math.sqrt(max(k, 0.0))


def phi(rho: float, R: float, d: float) -> float:
    """Angle of the slice circle inside the bore."""
    if rho <= R - d:
        return 2 * math.pi
    if rho >= R + d:
        return 0.0
    return 2 * math.acos(_cos_ab(rho, R, d)[0])


def psi(rho: float, R: float, d: float) -> float:
    """Angle of the bore circle inside the slice disc."""
    if rho <= R - d:
        return 0.0
    if rho >= R + d:
        return 2 * math.pi
    return 2 * math.acos(_cos_ab(rho, R, d)[1])


def quad_ref(bore: float, c, r: float, h: float) -> dict:
    R = bore / 2
    cx, cy, cz = c
    d = math.hypot(cx, cz)
    lo = max(-h - cy, -r)
    hi = min(-cy, r)
    out = {"V": 0.0, "Af": 0.0, "Wl": 0.0}
    b = -cy
    ta = abs(h + cy)
    out["Wh"] = lens(math.sqrt(r * r - b * b), R, d) if r > b else 0.0
    out["Wp"] = lens(math.sqrt(r * r - ta * ta), R, d) if r > ta else 0.0
    if hi <= lo:
        return out
    pts = []
    for R0 in (R - d, R + d):
        if r > R0:
            t = math.sqrt(r * r - R0 * R0)
            pts += [t, -t]
    pts += [0.0]
    pts = sorted(p for p in set(pts) if lo < p < hi)
    rho = lambda t: math.sqrt(max(r * r - t * t, 0.0))
    kw = dict(points=pts or None, limit=500, epsabs=0.0, epsrel=1e-12)
    out["V"] = integrate.quad(lambda t: lens(rho(t), R, d), lo, hi, **kw)[0]
    out["Af"] = integrate.quad(lambda t: r * phi(rho(t), R, d), lo, hi, **kw)[0]
    out["Wl"] = integrate.quad(lambda t: R * psi(rho(t), R, d), lo, hi, **kw)[0]
    return out


# ------------------------------------------------------------------------------------------------
# Monte-Carlo reference
# ------------------------------------------------------------------------------------------------
N_MC = 10_000_000
CHUNK = 2_000_000


def mc_ref(bore: float, c, r: float, h: float, rng: np.random.Generator) -> dict:
    R = bore / 2
    cx, cy, cz = map(float, c)
    d = math.hypot(cx, cz)
    res = {}

    def run(sampler, measure):
        hits = 0
        for _ in range(N_MC // CHUNK):
            hits += int(sampler(CHUNK))
        p = hits / N_MC
        return measure * p, measure * math.sqrt(max(p * (1 - p), 1.0 / N_MC) / N_MC)

    # volume: box bounding ball ∩ cylinder
    x0, x1 = max(cx - r, -R), min(cx + r, R)
    z0, z1 = max(cz - r, -R), min(cz + r, R)
    y0, y1 = max(cy - r, -h), min(cy + r, 0.0)
    if x1 > x0 and z1 > z0 and y1 > y0:
        box = (x1 - x0) * (y1 - y0) * (z1 - z0)

        def vol(n):
            x = rng.uniform(x0, x1, n)
            y = rng.uniform(y0, y1, n)
            z = rng.uniform(z0, z1, n)
            inside = ((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r) & (x * x + z * z <= R * R)
            return np.count_nonzero(inside)

        res["V"] = run(vol, box)
    else:
        res["V"] = (0.0, 0.0)

    # front area: uniform directions
    def front(n):
        v = rng.normal(size=(3, n))
        v /= np.sqrt((v * v).sum(axis=0))
        x = cx + r * v[0]
        y = cy + r * v[1]
        z = cz + r * v[2]
        return np.count_nonzero((x * x + z * z <= R * R) & (y <= 0.0) & (y >= -h))

    res["Af"] = run(front, 4 * math.pi * r * r)

    # planar wetted areas
    def plane(yp):
        dy = yp - cy
        if abs(dy) >= r:
            return (0.0, 0.0)
        rc = math.sqrt(r * r - dy * dy)
        a0, a1 = max(cx - rc, -R), min(cx + rc, R)
        b0, b1 = max(cz - rc, -R), min(cz + rc, R)
        if a1 <= a0 or b1 <= b0:
            return (0.0, 0.0)

        def f(n):
            x = rng.uniform(a0, a1, n)
            z = rng.uniform(b0, b1, n)
            return np.count_nonzero(((x - cx) ** 2 + (z - cz) ** 2 <= rc * rc) & (x * x + z * z <= R * R))

        return run(f, (a1 - a0) * (b1 - b0))

    res["Wh"] = plane(0.0)
    res["Wp"] = plane(-h)

    # liner: angle measured from the direction of the centre
    yl0, yl1 = max(-h, cy - r), min(0.0, cy + r)
    if d > 0:
        cosm = (R * R + d * d - r * r) / (2 * R * d)
        dpsi = math.pi if cosm <= -1 else (0.0 if cosm >= 1 else math.acos(cosm))
        psic = math.atan2(cz, cx)
    else:
        dpsi = math.pi if r > R else 0.0
        psic = 0.0
    if yl1 > yl0 and dpsi > 0:

        def lin(n):
            ps = rng.uniform(psic - dpsi, psic + dpsi, n)
            y = rng.uniform(yl0, yl1, n)
            x = R * np.cos(ps)
            z = R * np.sin(ps)
            return np.count_nonzero((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r)

        res["Wl"] = run(lin, R * 2 * dpsi * (yl1 - yl0))
    else:
        res["Wl"] = (0.0, 0.0)
    return res


def max_radius(bore, c, h):
    R = bore / 2
    d = math.hypot(c[0], c[2])
    b = -c[1]
    return math.hypot(R + d, max(b, abs(b - h)))


def mc_cases(bore, c):
    """Hand-picked (r, h) covering every regime."""
    R = bore / 2
    d = math.hypot(c[0], c[2])
    b = -c[1]
    c1, c2 = R - d, R + d
    cases = []
    for h in (0.0003, 0.007, 0.019, 0.13):
        rm = max_radius(bore, c, h)
        for r in (
            0.5 * min(b, h - b, c1) if h > b and min(b, h - b, c1) > 0 else abs(b - h) + 0.5 * h,  # inside / thin
            1.5 * b if b > 0 else 2e-3,  # touching the head
            c1 + 0.3 * (c2 - c1) if d > 0 else 0.6 * R,  # touching the liner
            c1 + 0.4 * min(c1, c2 - c1) if d > 0 else 0.3 * R,  # just after first liner contact
            c1 + 2.0 * min(c1, c2 - c1) if d > 0 else 0.4 * R,  # early flame along the liner
            0.5 * rm,
            c2 + 0.5 * (rm - c2) if rm > c2 else 0.9 * rm,
            0.97 * rm,  # near the far corner
        ):
            if 0 < r < rm:
                cases.append((r, h))
    return cases


def sig(x: float, n: int = 10) -> float:
    return float(f"{x:.{n}g}")


def build_geometry(gi: int):
    """MC cases + dense quad rows for GEOMETRIES[gi]. Each geometry has its own fixed-seed generators, so
    the output is deterministic and independent of the process pool."""
    name, bore, c = GEOMETRIES[gi]
    rng = np.random.default_rng(20260929 + 1000 * gi)
    worst_sigma = 0.0
    cases = []
    lines = []
    for r, h in mc_cases(bore, c):
        r, h = sig(r, 12), sig(h, 12)  # the fixture stores exactly these values
        q = quad_ref(bore, c, r, h)
        m = mc_ref(bore, c, r, h, rng)
        row = {"r": r, "h": h}
        for k in ("V", "Af", "Wh", "Wp", "Wl"):
            val, se = m[k]
            row[k] = [sig(val, 8), sig(se, 3), sig(q[k], 12)]
            if se > 0:
                worst_sigma = max(worst_sigma, abs(val - q[k]) / se)
                assert abs(val - q[k]) <= 5 * se + 1e-18, (name, r, h, k, val, se, q[k])
        cases.append(row)
        lines.append(f"{name:9s} r={r*1e3:8.3f}mm h={h*1e3:7.3f}mm  V={q['V']:.6e} Af={q['Af']:.6e} Wl={q['Wl']:.6e}")
    mc_geom = {"name": name, "bore": bore, "center": list(c), "cases": cases}

    # dense deterministic samples (quad only) for the table-accuracy test
    g = np.random.default_rng(12345 + gi)
    pts = []
    for _ in range(400):
        h = sig(0.0003 * math.exp(g.uniform() * math.log(0.16 / 0.0003)), 12)
        rm = max_radius(bore, c, h)
        r = sig(g.uniform() * rm, 12)
        q = quad_ref(bore, c, r, h)
        pts.append([r, h] + [sig(q[k], 12) for k in ("V", "Af", "Wh", "Wp", "Wl")])
    quad_geom = {"name": name, "bore": bore, "center": list(c), "columns": ["r", "h", "V", "Af", "Wh", "Wp", "Wl"], "rows": pts}
    return mc_geom, quad_geom, worst_sigma, lines


def main():
    from multiprocessing import Pool

    mc_out = {"description": __doc__.splitlines()[0], "nSamples": N_MC, "geometries": []}
    quad_out = {"description": "quad reference of sphere ∩ disc (see combustion_geometry_mc.py)", "geometries": []}
    worst_sigma = 0.0
    with Pool(len(GEOMETRIES)) as pool:
        results = pool.map(build_geometry, range(len(GEOMETRIES)))
    for mc_geom, quad_geom, ws, lines in results:
        print("\n".join(lines))
        mc_out["geometries"].append(mc_geom)
        quad_out["geometries"].append(quad_geom)
        worst_sigma = max(worst_sigma, ws)
    print(f"max |MC − quad| / σ = {worst_sigma:.2f}")
    mc_out["maxDeviationSigma"] = round(worst_sigma, 3)
    fx = os.path.join(ROOT, "test", "fixtures")
    with open(os.path.join(fx, "combustion_geometry_mc.json"), "w") as f:
        json.dump(mc_out, f, separators=(",", ":"))
    with open(os.path.join(fx, "combustion_geometry_quad.json"), "w") as f:
        json.dump(quad_out, f, separators=(",", ":"))


if __name__ == "__main__":
    main()
