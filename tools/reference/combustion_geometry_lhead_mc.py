"""Oracle for src/physics/combustion/chamber.ts (LHeadChamber) and the d > R branch of flame-geometry.ts.

Side-valve ('l-head') chamber, cylinder frame (core/engine-spec.ts LHeadChamberSpec):
    bore column  {x² + z² ≤ R², −h ≤ y ≤ 0}
    valve pocket {(x, z) ∈ rounded rectangle, x² + z² > R², deckY ≤ y ≤ roofY}
Sphere (centre c = spark gap, radius r) ∩ chamber: volume, flame-front area (sphere surface inside the
chamber), and the burned-wetted areas of the wall surfaces (definitions as in chamber.ts):
    head   = roof over the bore (y = 0) + pocket roof (y = roofY) + bore-circle wall above max(deckY, −h)
             minus its part open to the pocket (bore arc inside the rectangle, y ≤ roofY) + pocket side walls
    block  = pocket floor (y = deckY) minus the valve faces
    piston = crown (y = −h) + its side facing the pocket while the crown is above the deck
    liner  = bore wall below the deck
    intakeValve / exhaustValve = valve-head discs in the pocket floor
plus the bore-column parts alone (sphere ∩ disc: Vcol, Afcol, Wh, Wp, Wl) for the FlameGeometry d > R branch.

Two independent references, written WITHOUT the Green's-theorem / table machinery of chamber.ts:
  1. Dense quadrature (scipy.integrate.quad) of the horizontal slice integrals (Archimedes' hat-box
     theorem). Bore-column slices are circle–circle lenses (classical acos formulas). Pocket slice areas
     are 1-D integrals over x of the chord length {z : in slice circle, in rectangle, outside the bore};
     slice-circle arcs and wall lengths inside the slice circle come from dense angular / arc-length
     sampling of the inside test refined by bisection. Relative accuracy ~1e-8.
  2. Monte Carlo (numpy, fixed seed, ≥ 2e6 samples per quantity) on a subset of cases, with 1-σ binomial
     error bars; the script asserts |MC − quad| ≤ 5σ.
It also checks the oracle itself: A_f = dV/dr by central differences of the quadrature volume.
Writes test/fixtures/combustion_geometry_lhead_mc.json.
Run: .venv/bin/python tools/reference/combustion_geometry_lhead_mc.py
"""
from __future__ import annotations

import json
import math
import os

import numpy as np
from scipy import integrate

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
IN = 0.0254
PI = math.pi

# ------------------------------------------------------------------------------------------------
# geometry
# ------------------------------------------------------------------------------------------------


def rr_zrange(x, pk):
    """z-range of the rounded rectangle at x (or None)."""
    xa, xb, za, zb, c = pk["xMin"], pk["xMax"], pk["zMin"], pk["zMax"], pk["cornerRadius"]
    c = min(c, 0.5 * (xb - xa), 0.5 * (zb - za))
    if x < xa or x > xb:
        return None
    dz = 0.0
    if x < xa + c:
        e = xa + c - x
        dz = c - math.sqrt(max(c * c - e * e, 0.0))
    elif x > xb - c:
        e = x - (xb - c)
        dz = c - math.sqrt(max(c * c - e * e, 0.0))
    return za + dz, zb - dz


def in_rr(x, z, pk):
    """Vectorised rounded-rectangle inside test."""
    xa, xb, za, zb, c = pk["xMin"], pk["xMax"], pk["zMin"], pk["zMax"], pk["cornerRadius"]
    hx, hz = 0.5 * (xb - xa), 0.5 * (zb - za)
    c = min(c, hx, hz)
    ax = np.abs(x - 0.5 * (xa + xb))
    az = np.abs(z - 0.5 * (za + zb))
    ex = ax - (hx - c)
    ez = az - (hz - c)
    corner = (ex > 0) & (ez > 0)
    return (ax <= hx) & (az <= hz) & (~corner | (ex * ex + ez * ez <= c * c))


def rr_perimeter_point(s, pk):
    """Point at arc length s ∈ [0, L) along the rounded-rectangle boundary (vectorised), and L."""
    xa, xb, za, zb, c = pk["xMin"], pk["xMax"], pk["zMin"], pk["zMax"], pk["cornerRadius"]
    c = min(c, 0.5 * (xb - xa), 0.5 * (zb - za))
    ly = (zb - za) - 2 * c
    lx = (xb - xa) - 2 * c
    q = 0.5 * PI * c
    segs = [
        ("l", (xb, za + c), (xb, zb - c), ly),
        ("a", (xb - c, zb - c), 0.0, q),
        ("l", (xb - c, zb), (xa + c, zb), lx),
        ("a", (xa + c, zb - c), 0.5 * PI, q),
        ("l", (xa, zb - c), (xa, za + c), ly),
        ("a", (xa + c, za + c), PI, q),
        ("l", (xa + c, za), (xb - c, za), lx),
        ("a", (xb - c, za + c), 1.5 * PI, q),
    ]
    L = sum(sg[3] for sg in segs)
    s = np.mod(np.asarray(s, dtype=float), L)
    x = np.empty_like(s)
    z = np.empty_like(s)
    s0 = 0.0
    for kind, a, b, ln in segs:
        m = (s >= s0) & (s < s0 + ln) if ln > 0 else np.zeros_like(s, dtype=bool)
        if np.any(m):
            u = (s[m] - s0) / ln
            if kind == "l":
                x[m] = a[0] + u * (b[0] - a[0])
                z[m] = a[1] + u * (b[1] - a[1])
            else:
                ang = b + u * 0.5 * PI
                x[m] = a[0] + c * np.cos(ang)
                z[m] = a[1] + c * np.sin(ang)
        s0 += ln
    return x, z, L


class Geom:
    def __init__(self, name, bore, deck_y, roof_y, crown, pocket, spark, valves):
        self.name = name
        self.bore = bore
        self.R = bore / 2
        self.deck = deck_y
        self.roof = roof_y
        self.crown = crown
        self.pk = pocket
        self.c = spark
        self.valves = valves  # [(x, z, radius, count)] intake, exhaust
        self.d = math.hypot(spark[0], spark[2])
        self.pocket_area = pocket_plan_area(self)
        self.crit = pocket_critical_radii(self)

    def in_pocket(self, x, z):
        return in_rr(x, z, self.pk) & (x * x + z * z > self.R * self.R)


def pocket_plan_area(g: Geom) -> float:
    """Exact area of rectangle ∖ bore disc by 1-D quadrature over x (independent of chamber.ts)."""

    def length(x):
        zr = rr_zrange(x, g.pk)
        if zr is None:
            return 0.0
        lo, hi = zr
        L = hi - lo
        if abs(x) < g.R:
            cz = math.sqrt(g.R * g.R - x * x)
            L -= max(0.0, min(hi, cz) - max(lo, -cz))
        return L

    pk = g.pk
    pts = [pk["xMin"] + pk["cornerRadius"], pk["xMax"] - pk["cornerRadius"], -g.R, g.R]
    pts = sorted(p for p in pts if pk["xMin"] < p < pk["xMax"])
    return integrate.quad(length, pk["xMin"], pk["xMax"], points=pts or None, limit=500, epsabs=0, epsrel=1e-13)[0]


def boundary_samples(g: Geom, n=40000):
    """Dense samples of the pocket's plan boundary: (x, z, kind) with kind 0 = side wall, 1 = open arc."""
    s = np.linspace(0, 1, n, endpoint=False)
    x, z, L = rr_perimeter_point(s * 1e9, g.pk)  # placeholder to get L
    s = np.linspace(0, L, n, endpoint=False)
    x, z, _ = rr_perimeter_point(s, g.pk)
    keep = x * x + z * z > g.R * g.R
    th = np.linspace(0, 2 * PI, n, endpoint=False)
    bx, bz = g.R * np.cos(th), g.R * np.sin(th)
    keep2 = in_rr(bx, bz, g.pk)
    return np.concatenate([x[keep], bx[keep2]]), np.concatenate([z[keep], bz[keep2]])


def pocket_critical_radii(g: Geom):
    """Approximate critical radii (local distance extrema along the boundary) — quadrature breakpoints only."""
    x, z = boundary_samples(g, 20000)
    dist = np.hypot(x - g.c[0], z - g.c[2])
    out = {float(dist.min()), float(dist.max())}
    for i in range(1, len(dist) - 1):
        if (dist[i] - dist[i - 1]) * (dist[i + 1] - dist[i]) <= 0:
            out.add(float(dist[i]))
    return sorted(out)


# ------------------------------------------------------------------------------------------------
# closed forms for the bore column (classical acos lens formulas, any d ≥ 0)
# ------------------------------------------------------------------------------------------------


def lens(rho, R, d):
    if rho <= 0:
        return 0.0
    if rho <= R - d:
        return PI * rho * rho
    if rho >= R + d:
        return PI * R * R
    if rho <= d - R:
        return 0.0
    ca = max(-1.0, min(1.0, (d * d + rho * rho - R * R) / (2 * d * rho)))
    cb = max(-1.0, min(1.0, (d * d + R * R - rho * rho) / (2 * d * R)))
    k = (-d + rho + R) * (d + rho - R) * (d - rho + R) * (d + rho + R)
    return rho * rho * math.acos(ca) + R * R * math.acos(cb) - 0.5 * math.sqrt(max(k, 0.0))


def phi_bore(rho, R, d):
    """Angle of the slice circle inside the bore disc."""
    if rho <= R - d:
        return 2 * PI
    if rho >= R + d or rho <= d - R:
        return 0.0
    return 2 * math.acos(max(-1.0, min(1.0, (d * d + rho * rho - R * R) / (2 * d * rho))))


def psi_bore(rho, R, d):
    """Angle of the bore circle inside the slice disc."""
    if rho <= R - d or rho <= d - R:
        return 0.0
    if rho >= R + d:
        return 2 * PI
    return 2 * math.acos(max(-1.0, min(1.0, (d * d + R * R - rho * rho) / (2 * d * R))))


# ------------------------------------------------------------------------------------------------
# pocket slice quantities (independent numerics)
# ------------------------------------------------------------------------------------------------


def pocket_slice_area(g: Geom, rho: float) -> float:
    if rho <= 0:
        return 0.0
    px, pz = g.c[0], g.c[2]
    pk = g.pk
    xa, xb = max(px - rho, pk["xMin"]), min(px + rho, pk["xMax"])
    if xb <= xa:
        return 0.0
    R = g.R

    def length(x):
        w2 = rho * rho - (x - px) ** 2
        if w2 <= 0:
            return 0.0
        w = math.sqrt(w2)
        zr = rr_zrange(x, pk)
        if zr is None:
            return 0.0
        lo, hi = max(pz - w, zr[0]), min(pz + w, zr[1])
        if hi <= lo:
            return 0.0
        L = hi - lo
        if abs(x) < R:
            cz = math.sqrt(R * R - x * x)
            L -= max(0.0, min(hi, cz) - max(lo, -cz))
        return L

    pts = [pk["xMin"] + pk["cornerRadius"], pk["xMax"] - pk["cornerRadius"], -R, R]
    pts = sorted(p for p in pts if xa < p < xb)
    return integrate.quad(length, xa, xb, points=pts or None, limit=800, epsabs=0, epsrel=1e-11)[0]


def bisect_bool(f, a, b, fa, iters=55):
    for _ in range(iters):
        m = 0.5 * (a + b)
        if f(m) == fa:
            a = m
        else:
            b = m
    return 0.5 * (a + b)


def measure_set(param_to_inside, lo, hi, n=6000):
    """Measure of {s ∈ [lo, hi) : inside(s)} for a periodic-or-not parametrisation: dense sampling + bisection."""
    s = np.linspace(lo, hi, n + 1)
    ins = param_to_inside(s)
    total = 0.0
    start = lo if ins[0] else None
    for i in range(n):
        if ins[i] != ins[i + 1]:
            a, b = s[i], s[i + 1]
            fa = bool(ins[i])
            t = bisect_bool(lambda q: bool(param_to_inside(np.array([q]))[0]), a, b, fa)
            if fa:
                total += t - start
                start = None
            else:
                start = t
    if start is not None:
        total += hi - start
    return total


def pocket_arc_angle(g: Geom, rho: float) -> float:
    if rho <= 0:
        return 0.0
    px, pz = g.c[0], g.c[2]
    return measure_set(lambda th: g.in_pocket(px + rho * np.cos(th), pz + rho * np.sin(th)), 0.0, 2 * PI)


def open_arc_length(g: Geom, rho: float) -> float:
    """Length of the bore arc inside the rectangle (open boundary) within the slice circle of radius ρ."""
    if rho <= 0:
        return 0.0
    px, pz = g.c[0], g.c[2]
    R = g.R

    def ins(th):
        x, z = R * np.cos(th), R * np.sin(th)
        return in_rr(x, z, g.pk) & ((x - px) ** 2 + (z - pz) ** 2 < rho * rho)

    return R * measure_set(ins, 0.0, 2 * PI)


def side_wall_length(g: Geom, rho: float) -> float:
    """Length of the rectangle outline outside the bore disc within the slice circle of radius ρ."""
    if rho <= 0:
        return 0.0
    px, pz = g.c[0], g.c[2]
    _, _, L = rr_perimeter_point(np.array([0.0]), g.pk)

    def ins(s):
        x, z, _ = rr_perimeter_point(s, g.pk)
        return (x * x + z * z > g.R * g.R) & ((x - px) ** 2 + (z - pz) ** 2 < rho * rho)

    return measure_set(ins, 0.0, L)


# ------------------------------------------------------------------------------------------------
# quadrature reference
# ------------------------------------------------------------------------------------------------


def ybreaks(g: Geom, r: float, y0: float, y1: float, radii):
    cy = g.c[1]
    pts = []
    for rho in radii:
        if r > rho:
            t = math.sqrt(r * r - rho * rho)
            pts += [cy + t, cy - t]
    pts += [cy]
    return sorted(p for p in set(pts) if y0 < p < y1)


def yquad(g: Geom, r: float, y0: float, y1: float, f, radii, rel=1e-10):
    cy = g.c[1]
    lo, hi = max(y0, cy - r), min(y1, cy + r)
    if hi <= lo:
        return 0.0
    pts = ybreaks(g, r, lo, hi, radii)
    rho = lambda y: math.sqrt(max(r * r - (y - cy) ** 2, 0.0))
    return integrate.quad(lambda y: f(rho(y)), lo, hi, points=pts or None, limit=800, epsabs=0, epsrel=rel)[0]


def quad_ref(g: Geom, r: float, h: float) -> dict:
    R, d, cy = g.R, g.d, g.c[1]
    col_radii = [abs(R - d), R + d]
    out = {}
    rho_at = lambda y: math.sqrt(r * r - (y - cy) ** 2) if r > abs(y - cy) else 0.0
    # bore column (sphere ∩ disc)
    out["Vcol"] = yquad(g, r, -h, 0.0, lambda p: lens(p, R, d), col_radii)
    out["Afcol"] = yquad(g, r, -h, 0.0, lambda p: r * phi_bore(p, R, d), col_radii)
    out["Wl"] = yquad(g, r, -h, 0.0, lambda p: R * psi_bore(p, R, d), col_radii)
    out["Wh"] = lens(rho_at(0.0), R, d)
    out["Wp"] = lens(rho_at(-h), R, d)
    # pocket prism
    radii = g.crit
    out["Vp"] = yquad(g, r, g.deck, g.roof, lambda p: pocket_slice_area(g, p), radii, rel=1e-9)
    out["Afp"] = yquad(g, r, g.deck, g.roof, lambda p: r * pocket_arc_angle(g, p), radii, rel=1e-9)
    out["V"] = out["Vcol"] + out["Vp"]
    out["Af"] = out["Afcol"] + out["Afp"]
    # walls
    ylow = max(g.deck, -h)
    bore_above = yquad(g, r, ylow, 0.0, lambda p: R * psi_bore(p, R, d), col_radii)
    open_above = yquad(g, r, ylow, g.roof, lambda p: open_arc_length(g, p), radii, rel=1e-9)
    side = yquad(g, r, g.deck, g.roof, lambda p: side_wall_length(g, p), radii, rel=1e-9)
    roof = pocket_slice_area(g, rho_at(g.roof))
    floor = pocket_slice_area(g, rho_at(g.deck))
    rf = rho_at(g.deck)
    vf = [cnt * lens(rf, rv, math.hypot(vx - g.c[0], vz - g.c[2])) for (vx, vz, rv, cnt) in g.valves]
    out["head"] = out["Wh"] + roof + max(0.0, bore_above - open_above) + side
    out["block"] = max(0.0, floor - vf[0] - vf[1])
    out["intakeValve"] = vf[0]
    out["exhaustValve"] = vf[1]
    pside = 0.0
    if -h > g.deck:
        pside = yquad(g, r, g.deck, min(-h, g.roof), lambda p: open_arc_length(g, p), radii, rel=1e-9)
    out["piston"] = out["Wp"] + pside
    out["liner"] = yquad(g, r, -h, g.deck, lambda p: R * psi_bore(p, R, d), col_radii) if h > -g.deck else 0.0
    return out


# ------------------------------------------------------------------------------------------------
# Monte-Carlo reference
# ------------------------------------------------------------------------------------------------
N_MC = 4_000_000
CHUNK = 1_000_000


def mc_ref(g: Geom, r: float, h: float, rng: np.random.Generator) -> dict:
    R, cx, cy, cz = g.R, g.c[0], g.c[1], g.c[2]
    pk = g.pk
    res = {}

    def run(sampler, measure):
        if measure <= 0:
            return (0.0, 0.0)
        hits = 0
        for _ in range(N_MC // CHUNK):
            hits += int(sampler(CHUNK))
        p = hits / N_MC
        return measure * p, measure * math.sqrt(max(p * (1 - p), 1.0 / N_MC) / N_MC)

    def in_chamber(x, y, z):
        col = (x * x + z * z <= R * R) & (y >= -h) & (y <= 0)
        poc = g.in_pocket(x, z) & (y >= g.deck) & (y <= g.roof)
        return col | poc

    # volume: box bounding ball ∩ chamber
    x0, x1 = max(cx - r, min(-R, pk["xMin"])), min(cx + r, max(R, pk["xMax"]))
    z0, z1 = max(cz - r, min(-R, pk["zMin"])), min(cz + r, max(R, pk["zMax"]))
    y0, y1 = max(cy - r, min(-h, g.deck)), min(cy + r, 0.0)

    def vol(n):
        x = rng.uniform(x0, x1, n)
        y = rng.uniform(y0, y1, n)
        z = rng.uniform(z0, z1, n)
        return np.count_nonzero(((x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r) & in_chamber(x, y, z))

    res["V"] = run(vol, (x1 - x0) * (y1 - y0) * (z1 - z0))

    def front(n):
        v = rng.normal(size=(3, n))
        v /= np.sqrt((v * v).sum(axis=0))
        return np.count_nonzero(in_chamber(cx + r * v[0], cy + r * v[1], cz + r * v[2]))

    res["Af"] = run(front, 4 * PI * r * r)

    def ball2(x, y, z):
        return (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r * r

    # planes: uniform in the plane's bounding box, indicator of the surface and the ball
    def plane(yp, surf, box):
        (a0, a1), (b0, b1) = box

        def f(n):
            x = rng.uniform(a0, a1, n)
            z = rng.uniform(b0, b1, n)
            return np.count_nonzero(surf(x, z) & ball2(x, yp, z))

        return run(f, (a1 - a0) * (b1 - b0))

    bore_box = ((-R, R), (-R, R))
    pk_box = ((pk["xMin"], pk["xMax"]), (pk["zMin"], pk["zMax"]))
    in_disc = lambda x, z: x * x + z * z <= R * R
    roof = plane(g.roof, g.in_pocket, pk_box)
    wh = plane(0.0, in_disc, bore_box)
    wp = plane(-h, in_disc, bore_box)

    def not_valve(x, z):
        m = np.ones_like(x, dtype=bool)
        for vx, vz, rv, _ in g.valves:
            m &= (x - vx) ** 2 + (z - vz) ** 2 > rv * rv
        return m

    block = plane(g.deck, lambda x, z: g.in_pocket(x, z) & not_valve(x, z), pk_box)
    vals = []
    for vx, vz, rv, cnt in g.valves:
        a = plane(g.deck, lambda x, z, vx=vx, vz=vz, rv=rv: (x - vx) ** 2 + (z - vz) ** 2 <= rv * rv, ((vx - rv, vx + rv), (vz - rv, vz + rv)))
        vals.append((cnt * a[0], cnt * a[1]))

    # walls
    def bore_wall(ya, yb, keep):
        if yb <= ya:
            return (0.0, 0.0)

        def f(n):
            th = rng.uniform(0, 2 * PI, n)
            y = rng.uniform(ya, yb, n)
            x, z = R * np.cos(th), R * np.sin(th)
            return np.count_nonzero(keep(x, y, z) & ball2(x, y, z))

        return run(f, 2 * PI * R * (yb - ya))

    ylow = max(g.deck, -h)
    head_wall = bore_wall(ylow, 0.0, lambda x, y, z: ~(in_rr(x, z, pk) & (y <= g.roof)))
    _, _, Lrr = rr_perimeter_point(np.array([0.0]), pk)

    def side_f(n):
        s = rng.uniform(0, Lrr, n)
        y = rng.uniform(g.deck, g.roof, n)
        x, z, _ = rr_perimeter_point(s, pk)
        return np.count_nonzero((x * x + z * z > R * R) & ball2(x, y, z))

    side = run(side_f, Lrr * (g.roof - g.deck))
    pside = bore_wall(g.deck, min(-h, g.roof), lambda x, y, z: in_rr(x, z, pk)) if -h > g.deck else (0.0, 0.0)
    liner = bore_wall(-h, g.deck, lambda x, y, z: np.ones_like(x, dtype=bool)) if h > -g.deck else (0.0, 0.0)

    def add(*terms):
        return (sum(t[0] for t in terms), math.sqrt(sum(t[1] ** 2 for t in terms)))

    res["head"] = add(wh, roof, head_wall, side)
    res["block"] = block
    res["intakeValve"] = vals[0]
    res["exhaustValve"] = vals[1]
    res["piston"] = add(wp, pside)
    res["liner"] = liner
    return res


# ------------------------------------------------------------------------------------------------
# cases
# ------------------------------------------------------------------------------------------------


def model_t_geometry():
    """The Model T spec chamber (engines/model-t.ts), rebuilt here from its sourced inputs."""
    bore = 3.75 * IN
    stroke = 4.0 * IN
    ap = PI * bore * bore / 4
    vd = ap * stroke
    vc = vd / (3.98 - 1)
    crevice = 2.0e-6 + 0.5e-6
    crown = (5 / 16) * IN
    h_tdc = 1.0 * IN
    deck = -(h_tdc + crown)
    pocket = {"xMin": -0.0915, "xMax": -0.03, "zMin": -0.0475, "zMax": 0.0475, "cornerRadius": 0.02}
    tmp = Geom("tmp", bore, deck, deck + 0.01, crown, pocket, (-0.06, deck + 0.005, 0.0), [])
    area = tmp.pocket_area
    height = (vc - crevice - ap * h_tdc) / area
    roof = deck + height
    rv = 0.5 * 1.47 * IN
    valves = [(-0.068, -0.024, rv, 1), (-0.068, 0.024, rv, 1)]
    return dict(bore=bore, deck=deck, roof=roof, crown=crown, pocket=pocket, valves=valves, h_tdc=h_tdc, stroke=stroke)


def build_geometries():
    mt = model_t_geometry()
    gs = []
    gs.append(Geom("modelt", mt["bore"], mt["deck"], mt["roof"], mt["crown"], mt["pocket"], (-0.06, mt["roof"] - 0.004, 0.0), mt["valves"]))
    # spark inside the bore planform (plug over the piston), same chamber
    gs.append(Geom("modelt_bore_spark", mt["bore"], mt["deck"], mt["roof"], mt["crown"], mt["pocket"], (-0.02, -0.005, 0.01), mt["valves"]))
    # spark near the far end of the pocket, off the centre line
    gs.append(Geom("modelt_far", mt["bore"], mt["deck"], mt["roof"], mt["crown"], mt["pocket"], (-0.08, mt["roof"] - 0.003, 0.02), mt["valves"]))
    # a different L-head: smaller bore, sharper pocket corners, crown below the deck at TDC
    pk2 = {"xMin": -0.075, "xMax": -0.025, "zMin": -0.03, "zMax": 0.035, "cornerRadius": 0.008}
    rv2 = 0.014
    gs.append(Geom("generic", 0.08, -0.02, -0.008, -0.003, pk2, (-0.05, -0.012, 0.004), [(-0.057, -0.013, rv2, 1), (-0.057, 0.017, rv2, 1)]))
    return gs, mt


def heights(g: Geom):
    htdc = -g.deck - g.crown
    return [htdc, 0.5 * (htdc - g.deck), -g.deck, -g.deck + 0.012, 0.08, 0.125]


def radii_for(g: Geom, h: float):
    R, d, cy = g.R, g.d, g.c[1]
    rmax_col = math.hypot(R + d, max(-cy, abs(-cy - h)))
    x, z = boundary_samples(g, 8000)
    rho_max = float(np.hypot(x - g.c[0], z - g.c[2]).max())
    rmax_p = math.hypot(rho_max, max(abs(g.deck - cy), abs(g.roof - cy)))
    rm = max(rmax_col, rmax_p)
    cand = [0.003, 0.006, abs(R - d) + 0.002, abs(R - d) + 0.01, 0.03, 0.045, 0.06, 0.075, 0.09, 0.5 * rm, 0.8 * rm, 0.95 * rm, 0.99 * rm]
    return sorted({sig(r, 12) for r in cand if 0 < r < rm}), rm


def sig(x: float, n: int = 10) -> float:
    return float(f"{x:.{n}g}")


KEYS = ["V", "Af", "head", "piston", "liner", "block", "intakeValve", "exhaustValve"]
COL_KEYS = ["Vcol", "Afcol", "Wh", "Wp", "Wl", "Vp", "Afp"]


def build_geometry(gi: int):
    gs, _ = build_geometries()
    g = gs[gi]
    rng = np.random.default_rng(20261003 + 1000 * gi)
    cases = []
    lines = []
    worst_sigma = 0.0
    worst_fd = 0.0
    for hi, h in enumerate(heights(g)):
        h = sig(h, 12)
        rs, rm = radii_for(g, h)
        for ri, r in enumerate(rs):
            q = quad_ref(g, r, h)
            row = {"r": r, "h": h}
            do_mc = (ri % 3 == 1) and hi in (0, 2, 4)
            m = mc_ref(g, r, h, rng) if do_mc else None
            for k in KEYS:
                if m is not None:
                    val, se = m[k]
                    row[k] = [sig(val, 8), sig(se, 3), sig(q[k], 12)]
                    if se > 0:
                        worst_sigma = max(worst_sigma, abs(val - q[k]) / se)
                        assert abs(val - q[k]) <= 5 * se + 1e-15, (g.name, r, h, k, val, se, q[k])
                else:
                    row[k] = [None, None, sig(q[k], 12)]
            for k in COL_KEYS:
                row[k] = sig(q[k], 12)
            # oracle self-check: A_f = dV/dr
            if ri % 4 == 2:
                dr = 1e-5 * r
                vp = quad_ref_volume(g, r + dr, h)
                vm = quad_ref_volume(g, r - dr, h)
                fd = (vp - vm) / (2 * dr)
                e = abs(fd - q["Af"]) / max(q["Af"], 1e-9)
                worst_fd = max(worst_fd, e)
                assert e < 2e-5, (g.name, r, h, fd, q["Af"])
            cases.append(row)
            lines.append(f"{g.name:18s} h={h*1e3:7.3f}mm r={r*1e3:8.3f}mm V={q['V']:.6e} Af={q['Af']:.6e} head={q['head']:.4e}")
    geom = {
        "name": g.name,
        "bore": g.bore,
        "deckY": g.deck,
        "roofY": g.roof,
        "crownAboveDeckAtTDC": g.crown,
        "pocket": g.pk,
        "spark": list(g.c),
        "valves": [{"position": [vx, vz], "headDiameter": 2 * rv, "count": cnt} for (vx, vz, rv, cnt) in g.valves],
        "pocketArea": g.pocket_area,
        "cases": cases,
    }
    return geom, worst_sigma, worst_fd, lines


def quad_ref_volume(g: Geom, r: float, h: float) -> float:
    R, d = g.R, g.d
    v = yquad(g, r, -h, 0.0, lambda p: lens(p, R, d), [abs(R - d), R + d], rel=1e-12)
    v += yquad(g, r, g.deck, g.roof, lambda p: pocket_slice_area(g, p), g.crit, rel=1e-11)
    return v


def main():
    from multiprocessing import Pool

    gs, mt = build_geometries()
    with Pool(len(gs)) as pool:
        results = pool.map(build_geometry, range(len(gs)))
    out = {
        "description": __doc__.splitlines()[0],
        "nSamples": N_MC,
        "modelT": {"pocketArea": gs[0].pocket_area, "roofY": mt["roof"], "deckY": mt["deck"]},
        "geometries": [],
    }
    worst_sigma = 0.0
    worst_fd = 0.0
    for geom, ws, wf, lines in results:
        print("\n".join(lines))
        out["geometries"].append(geom)
        worst_sigma = max(worst_sigma, ws)
        worst_fd = max(worst_fd, wf)
    print(f"max |MC − quad| / σ = {worst_sigma:.2f};  max |FD dV/dr − A_f| / A_f = {worst_fd:.2e}")
    print(f"Model T pocket plan area = {gs[0].pocket_area!r} m²")
    out["maxDeviationSigma"] = round(worst_sigma, 3)
    out["maxOracleFdError"] = float(f"{worst_fd:.3g}")
    with open(os.path.join(ROOT, "test", "fixtures", "combustion_geometry_lhead_mc.json"), "w") as f:
        json.dump(out, f, separators=(",", ":"))


if __name__ == "__main__":
    main()
