"""Oracle for the L-head knock acoustics of src/physics/chemistry/knock.ts (GridModeSet, GridEndGasProjector).

The Ford Model T chamber (src/physics/engines/model-t.ts): bore column (disc of the bore, roof y = 0 down to
the crown at TDC) + valve pocket (rounded rectangle minus the bore disc, deck to pocket roof), joined across the
bore circle only through the window max(deck, -depth) <= y <= roof (at TDC the crown stands 7.9 mm above the deck:
a 5 mm opening). Transverse modes of the depth-weighted Helmholtz problem div(H grad psi) + k^2 H psi = 0,
H dpsi/dn = 0 on the outline, with the window's junction resistance r on the bore <-> pocket faces
([psi] = r H dpsi/dn: knock.ts windowResistance / lHeadPlanform.interfaceResistance).

Writes test/fixtures/chemistry_knock_lhead_modes.json:
  geometry     : the Model T numbers, recomputed here with the expressions of model-t.ts
  window       : junction resistance r of two flat layers joined through a window in a zero-thickness wall: the
                 knock.ts Ritz method (Thomson's principle, Chebyshev basis with the 1/sqrt edge weight, image
                 Green's functions) re-implemented here in matrix form; the exact conformal-map value of a
                 mirror-symmetric slit (equipotential window: r = -(2/pi) ln[(cos(pi y0/H) - cos(pi y1/H))/2]);
                 and an independent 2-D finite-volume solve of the vertical section (graded mesh, potential flow)
  sameGrid     : the TS discretisation (cut-cell finite volumes, cell = B/64, window resistance on the faces
                 between bore-side and pocket-side cells) rebuilt independently here and
                 solved with scipy.sparse.linalg.eigsh (ARPACK shift-invert, not the TS block Lanczos):
                 alpha_j = k_j B/2, cell / sub-column counts, chamber volume, the end-gas volume outside
                 a few ball radii and the sign-invariant modal weights psi_j(sensor) * S_j (S_j = <psi_j>_eg /
                 <psi_j^2>) of the gas outside the ball about the plug that leaves 20 % end gas
  sameGridNoWindow : alpha_j of the same grid with every region boundary an ordinary depth step (the model before
                 the window resistance: a planform without interfaceResistance)
  converged    : the same problem on finer grids (B/128, B/192): grid convergence of alpha_j and the weights
  depth15      : alpha_j at the bore-column depth of 15 deg ATDC (B/96)
  threeD       : the 3-D Helmholtz problem of the same rigid flat-roofed columns (no depth averaging: finite
                 volumes on voxels = planform staircase x y-layers refined toward the window edges, z-symmetric
                 halves, scipy eigsh) at TDC and 15 deg ATDC on B/64 and B/96 staircases, with the depth-averaged
                 model on the same staircase with and without the window resistance
  disc         : same-grid alpha of a flat disc of the CFR bore (82.55 mm), scipy jnp_zeros for reference
  channel      : two-depth channel (L1, H1 | L2, H2): analytic k from H1 tan(k L1) + H2 tan(k L2) = 0 (brentq)
  endGasVolume : V_out(r) of the EXACT chamber geometry (ball about the plug), midpoint quadrature on a
                 0.025 mm planform grid with the analytic chord -- independent of the cut-cell construction
Run: .venv/bin/python tools/reference/chemistry_knock_lhead_modes.py   (about 8-10 min, most of it the 3-D solves)
"""
from __future__ import annotations

import math
import os
import time

import numpy as np
from scipy import optimize, sparse, special
from scipy.sparse.linalg import eigsh, spsolve

from chemistry_common import FIXTURES, write_json

NS = 8  # sub-samples per cell side (knock.ts GRID_SUBSAMPLES)
UNIFORM = 0.75  # knock.ts GRID_UNIFORM_DISTANCE

# --------------------------------------------------------------------------------------------
# Model T chamber (engines/model-t.ts)
# --------------------------------------------------------------------------------------------
IN = 0.0254
BORE = 3.75 * IN
STROKE = 4.0 * IN
A_PISTON = (math.pi * BORE * BORE) / 4
VD = A_PISTON * STROKE
CR = 3.98
VC = VD / (CR - 1)
CREVICE = 2.0e-6 + 0.5e-6
CROWN = (5 / 16) * IN
H_TDC = 1.0 * IN
DECK_Y = -(H_TDC + CROWN)
POCKET = dict(xMin=-0.0915, xMax=-0.03, zMin=-0.0475, zMax=0.0475, cornerRadius=0.02)
POCKET_HEIGHT = (VC - CREVICE - A_PISTON * H_TDC) / 46.03020157589403e-4
ROOF_Y = DECK_Y + POCKET_HEIGHT
GAP = (-0.06, DECK_Y + POCKET_HEIGHT - 0.004, 0.0)
ROD = 7.0 * IN


def depth_at(deg):
    """Bore-column depth (roof -> crown) at deg ATDC: slider crank, crank radius stroke/2, rod 7 in."""
    th = math.radians(deg)
    rc = STROKE / 2
    return H_TDC + rc * (1 - math.cos(th)) + ROD - math.sqrt(ROD**2 - (rc * math.sin(th)) ** 2)


# --------------------------------------------------------------------------------------------
# Window junction resistance (knock.ts windowResistance)
# --------------------------------------------------------------------------------------------
RITZ_N, RITZ_Q = 12, 256  # knock.ts WINDOW_RITZ_BASIS, WINDOW_RITZ_NODES


def window_r(lo1, hi1, lo2, hi2, w=None, N=RITZ_N, Q=RITZ_Q):
    """Ritz value of r = min E(u) over int u = 1 (Thomson's principle): the knock.ts method in matrix form."""
    w0, w1 = w if w else (max(lo1, lo2), min(hi1, hi2))
    if not w1 > w0:
        return math.inf
    if lo1 == lo2 and hi1 == hi2 and w0 == lo1 and w1 == hi1:
        return 0.0
    c, h = 0.5 * (w0 + w1), 0.5 * (w1 - w0)
    th = (np.arange(Q) + 0.5) * math.pi / Q
    y = c + h * np.cos(th)
    T = np.cos(np.outer(np.arange(N), th))
    D = y[:, None] - y[None, :]
    S = np.zeros((Q, Q))
    for lo, H in ((lo1, hi1 - lo1), (lo2, hi2 - lo2)):
        near = np.abs(D) <= 1e-12 * H
        with np.errstate(divide="ignore"):  # the D = 0 entries take the limit ln(pi/H)
            S -= np.where(near, math.log(math.pi / H), np.log(np.abs(2 * np.sin(math.pi * D / (2 * H))) / np.where(near, 1.0, np.abs(D))))
        S -= np.log(np.abs(2 * np.sin(math.pi * (y[:, None] + y[None, :] - 2 * lo) / (2 * H))))
    S /= math.pi
    M = (h * math.pi / Q) ** 2 * (T @ S @ T.T)
    M[0, 0] += -2 * math.pi * h * h * math.log(h / 2)
    j = np.arange(1, N)
    M[j, j] += math.pi * h * h / j
    s = np.zeros(N)
    s[0] = math.pi * h
    return float(1 / (s @ np.linalg.solve(M, s)))


def slit_exact(H, y0, y1):
    """Mirror-symmetric slit [y0, y1] in a zero-thickness diaphragm across a layer [0, H] (conformal map)."""
    return -(2 / math.pi) * math.log((math.cos(math.pi * y0 / H) - math.cos(math.pi * y1 / H)) / 2)


def _graded(a, b, pts, h0, hmax, growth):
    nodes = [a]
    x = a
    while x < b - 1e-15 * max(1.0, abs(b)):
        d = min(abs(x - p) for p in pts)
        nxt = x + min(hmax, h0 + growth * d)
        for p in pts:
            if x + 1e-15 < p < nxt - 1e-15:
                nxt = p
        nodes.append(min(nxt, b))
        x = nodes[-1]
    return np.array(nodes)


def section_fv(lo1, hi1, lo2, hi2, w=None, h0=1.25e-5, growth=0.08, hfrac=40, Lfac=6):
    """Potential flow through the vertical section by finite volumes (mesh refined toward the window edges and
    the wall): unit flux in at x = -L (uniform), phi = 0 at x = +L; r = phi(-L) - L/H1 - L/H2."""
    w0, w1 = w if w else (max(lo1, lo2), min(hi1, hi2))
    H1, H2 = hi1 - lo1, hi2 - lo2
    hmax = min(H1, H2) / hfrac
    L = Lfac * max(H1, H2)
    ye = _graded(min(lo1, lo2), max(hi1, hi2), sorted({lo1, hi1, lo2, hi2, w0, w1}), h0, hmax, growth)
    xr = _graded(0.0, L, [0.0], h0, hmax, growth)
    xe = np.concatenate([-xr[::-1], xr[1:]])
    i0 = len(xr) - 1  # xe[i0] = 0: the wall, open only over the window
    xc, yc = 0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])
    dx, dy = np.diff(xe), np.diff(ye)
    nx, ny = len(xc), len(yc)
    X, Y = np.meshgrid(xc, yc, indexing="ij")
    act = ((X < 0) & (Y > lo1) & (Y < hi1)) | ((X > 0) & (Y > lo2) & (Y < hi2))
    idx = -np.ones((nx, ny), dtype=int)
    n = int(act.sum())
    idx[act] = np.arange(n)
    m = act[:-1] & act[1:]
    m[i0 - 1, :] &= (yc > w0) & (yc < w1)
    gx = dy[None, :] / (xc[1:] - xc[:-1])[:, None] * np.ones((nx - 1, ny))
    m2 = act[:, :-1] & act[:, 1:]
    gy = dx[:, None] / (yc[1:] - yc[:-1])[None, :] * np.ones((nx, ny - 1))
    r = np.concatenate([idx[:-1][m], idx[:, :-1][m2]])
    c = np.concatenate([idx[1:][m], idx[:, 1:][m2]])
    v = np.concatenate([gx[m], gy[m2]])
    diag = np.bincount(r, v, n) + np.bincount(c, v, n)
    mr, ml = act[-1], act[0]
    diag[idx[-1][mr]] += dy[mr] / (0.5 * dx[-1])
    rhs = np.zeros(n)
    rhs[idx[0][ml]] += dy[ml] / H1
    ii = np.arange(n)
    A = sparse.coo_matrix((np.concatenate([-v, -v, diag]), (np.concatenate([r, c, ii]), np.concatenate([c, r, ii]))), shape=(n, n)).tocsc()
    phi = spsolve(A, rhs)
    phiL = (phi[idx[0][ml]] * dy[ml]).sum() / dy[ml].sum() + 0.5 * dx[0] / H1
    return float(phiL - (L / H1 + L / H2))


def rr_dist(x, z, p):
    hx = 0.5 * (p["xMax"] - p["xMin"]) - p["cornerRadius"]
    hz = 0.5 * (p["zMax"] - p["zMin"]) - p["cornerRadius"]
    qx = abs(x - 0.5 * (p["xMin"] + p["xMax"])) - hx
    qz = abs(z - 0.5 * (p["zMin"] + p["zMax"])) - hz
    return math.hypot(max(qx, 0.0), max(qz, 0.0)) + min(max(qx, qz), 0.0) - p["cornerRadius"]


class LHead:
    def __init__(self, bore, depth, window=True):
        self.R = bore / 2
        R = self.R
        p = POCKET
        self.xMin, self.xMax = min(-R, p["xMin"]), max(R, p["xMax"])
        self.zMin, self.zMax = min(-R, p["zMin"]), max(R, p["zMax"])
        self.yLo = [-depth, DECK_Y]
        self.yHi = [0.0, ROOF_Y]
        self.rwin = window_r(-depth, 0.0, DECK_Y, ROOF_Y) if window else 0.0

    def iface(self, r0, r1, x, z, ex, ez):
        """Junction resistance of a bore <-> pocket face at (x, z), normal (ex, ez): r/|n.e|, n radial."""
        if r0 == r1 or r0 < 0 or r1 < 0 or not self.rwin > 0:
            return 0.0
        rho = math.hypot(x, z)
        ne = abs(x * ex + z * ez) / rho if rho > 0 else 1.0
        return self.rwin / max(ne, 1e-12)

    def region(self, x, z):
        if x * x + z * z <= self.R * self.R:
            return 0
        return 1 if rr_dist(x, z, POCKET) <= 0 else -1

    def bdist(self, x, z):
        return min(abs(math.hypot(x, z) - self.R), abs(rr_dist(x, z, POCKET)))


class Disc:
    def __init__(self, radius, depth):
        self.R = radius
        self.xMin, self.xMax, self.zMin, self.zMax = -radius, radius, -radius, radius
        self.yLo = [-depth]
        self.yHi = [0.0]

    def region(self, x, z):
        return 0 if x * x + z * z <= self.R * self.R else -1

    def bdist(self, x, z):
        return abs(math.hypot(x, z) - self.R)


class Channel:
    """Strip [0, L1 + L2] x [0, W]: depth H1 for x < L1, H2 beyond (depth step across the strip)."""

    def __init__(self, L1, L2, H1, H2, W):
        self.L1, self.L, self.W = L1, L1 + L2, W
        self.xMin, self.xMax, self.zMin, self.zMax = 0.0, L1 + L2, 0.0, W
        self.yLo = [-H1, -H2]
        self.yHi = [0.0, 0.0]

    def region(self, x, z):
        if x < 0 or x > self.L or z < 0 or z > self.W:
            return -1
        return 0 if x < self.L1 else 1

    def bdist(self, x, z):
        return min(abs(x), abs(x - self.L1), abs(x - self.L), abs(z), abs(z - self.W))


# --------------------------------------------------------------------------------------------
# The knock.ts GridModeSet discretisation (independent re-implementation)
# --------------------------------------------------------------------------------------------
def build(pf, dx):
    nx = max(1, math.ceil((pf.xMax - pf.xMin) / dx - 1e-9))
    nz = max(1, math.ceil((pf.zMax - pf.zMin) / dx - 1e-9))
    gx0 = 0.5 * (pf.xMin + pf.xMax) - 0.5 * nx * dx
    gz0 = 0.5 * (pf.zMin + pf.zMax) - 0.5 * nz * dx
    H = [hi - lo for lo, hi in zip(pf.yLo, pf.yHi)]
    hof = lambda r: H[r] if r >= 0 else 0.0
    dUni = UNIFORM * dx
    off = [((a + 0.5) / NS - 0.5) * dx for a in range(NS)]
    sub_area = dx * dx / (NS * NS)
    cell_index = -np.ones((nx, nz), dtype=int)
    caps, cxs, czs, sides = [], [], [], []
    cols = []  # (x, z, area, yLo, yHi, cell)
    for k in range(nz):
        for i in range(nx):
            x = gx0 + (i + 0.5) * dx
            z = gz0 + (k + 0.5) * dx
            mine = []
            side = -1  # the region holding most of the cell's plan area (ties: the lower index)
            if pf.bdist(x, z) > dUni:
                r = pf.region(x, z)
                side = r
                if r >= 0:
                    mine.append((x, z, dx * dx, pf.yLo[r], pf.yHi[r]))
            else:
                acc = {}
                for b in range(NS):
                    zs = z + off[b]
                    for a in range(NS):
                        xs = x + off[a]
                        r = pf.region(xs, zs)
                        if r < 0:
                            continue
                        c = acc.setdefault(r, [0, 0.0, 0.0])
                        c[0] += 1
                        c[1] += xs
                        c[2] += zs
                for r in sorted(acc):
                    n, sx, sz = acc[r]
                    mine.append((sx / n, sz / n, n * sub_area, pf.yLo[r], pf.yHi[r]))
                if acc:
                    side = max(sorted(acc), key=lambda r: acc[r][0])  # first maximum = lowest index on ties
            cap = sum(a * (hi - lo) for (_, _, a, lo, hi) in mine)
            if cap > 0:
                c = len(caps)
                cell_index[i, k] = c
                caps.append(cap)
                cxs.append(x)
                czs.append(z)
                sides.append(side)
                for m in mine:
                    cols.append(m + (c,))

    def face(xf, zf, ex, ez, rj):
        # rj: junction resistance of the face (cells on two sides of the window), in series on every line
        if pf.bdist(xf, zf) > dUni:
            h = hof(pf.region(xf, zf))
            return dx / (dx / h + rj) if rj > 0 and h > 0 else h
        tot = 0.0
        for u in off:
            px, pz = xf + u * ez, zf + u * ex
            hf = hof(pf.region(px, pz))
            if not hf > 0:
                continue
            inv, n_in = 0.0, 0
            for v in off:
                h = hof(pf.region(px + v * ex, pz + v * ez))
                if h > 0:
                    inv += 1 / h
                    n_in += 1
            if n_in == 0:
                tot += dx / (dx / hf + rj) if rj > 0 else hf
            else:
                tot += dx / (dx * inv / n_in + rj) if rj > 0 else n_in / inv
        return tot / NS

    n = len(caps)
    rows, colsI, vals = [], [], []
    diag = np.zeros(n)
    for k in range(nz):
        for i in range(nx):
            c = cell_index[i, k]
            if c < 0:
                continue
            for (di, dk) in ((1, 0), (0, 1)):
                ii, kk = i + di, k + dk
                if ii >= nx or kk >= nz or cell_index[ii, kk] < 0:
                    continue
                d = cell_index[ii, kk]
                xf, zf = cxs[c] + 0.5 * dx * di, czs[c] + 0.5 * dx * dk
                rj = pf.iface(sides[c], sides[d], xf, zf, di, dk) if hasattr(pf, "iface") and sides[c] != sides[d] else 0.0
                G = face(xf, zf, di, dk, rj)
                if not G > 0:
                    continue
                rows += [c, d]
                colsI += [d, c]
                vals += [-G, -G]
                diag[c] += G
                diag[d] += G
    K = sparse.coo_matrix((vals, (rows, colsI)), shape=(n, n)).tocsr() + sparse.diags(diag)
    C = np.array(caps)
    cols = np.array(cols)
    return dict(K=K, C=C, nx=nx, nz=nz, gx0=gx0, gz0=gz0, dx=dx, cell_index=cell_index, cols=cols, n=n, pf=pf)


def modes(g, bore, max_alpha=7.1, kreq=48):
    K, C = g["K"], g["C"]
    ci = 1 / np.sqrt(C)
    A = sparse.diags(ci) @ K @ sparse.diags(ci)
    k = min(kreq, g["n"] - 2)
    lam, v = eigsh(A, k=k, sigma=-1 / bore**2, which="LM", tol=1e-13)
    o = np.argsort(lam)
    lam, v = lam[o], v[:, o]
    lam_cut = (2 * max_alpha / bore) ** 2
    assert lam[-1] > lam_cut, "request more eigenpairs"
    keep = (lam <= lam_cut)
    keep[0] = False  # the uniform mode
    lam, v = lam[keep], v[:, keep]
    psi = ci[:, None] * v
    # normalise: max |psi_j| = 1, largest entry positive (knock.ts)
    for j in range(psi.shape[1]):
        i = int(np.argmax(np.abs(psi[:, j])))
        psi[:, j] /= psi[i, j]
    alpha = np.sqrt(lam) * bore / 2
    ms = (C[:, None] * psi**2).sum(0) / C.sum()
    return alpha, psi, ms


def mode_shape(g, psi, x, z):
    if g["pf"].region(x, z) < 0:
        return np.zeros(psi.shape[1])
    dx = g["dx"]
    fx = (x - g["gx0"]) / dx - 0.5
    fz = (z - g["gz0"]) / dx - 0.5
    i0, k0 = math.floor(fx), math.floor(fz)
    tx, tz = fx - i0, fz - k0
    s = np.zeros(psi.shape[1])
    w = 0.0
    for dk in (0, 1):
        for di in (0, 1):
            i, k = i0 + di, k0 + dk
            if i < 0 or k < 0 or i >= g["nx"] or k >= g["nz"]:
                continue
            c = g["cell_index"][i, k]
            if c < 0:
                continue
            wt = (tx if di else 1 - tx) * (tz if dk else 1 - tz)
            s += wt * psi[c]
            w += wt
    return s / w if w > 0 else s


def outside(g, psi, spark, r):
    """End-gas volume outside the ball of radius r about the spark and its modal integrals (sub-column chords)."""
    cols = g["cols"]
    x, z, a, lo, hi, cell = cols.T
    rho2 = (x - spark[0]) ** 2 + (z - spark[2]) ** 2
    dlo, dhi = lo - spark[1], hi - spark[1]
    w = np.sqrt(np.maximum(0.0, r * r - rho2))
    chord = np.clip(np.minimum(dhi, w) - np.maximum(dlo, -w), 0.0, None)
    chord[r * r <= rho2] = 0.0
    veg = a * (dhi - dlo - chord)
    integ = (veg[:, None] * psi[cell.astype(int)]).sum(0)
    return veg.sum(), integ


def projector_checks(g, alpha, psi, ms, spark, sensor, radii):
    vtot = g["C"].sum()
    vout = [outside(g, psi, spark, r)[0] for r in radii]
    r20 = optimize.brentq(lambda r: outside(g, psi, spark, r)[0] - 0.2 * vtot, 0.0, 0.2, xtol=1e-15)
    v, integ = outside(g, psi, spark, r20)
    S = integ / v / ms
    ps = mode_shape(g, psi, sensor[0], sensor[1])
    return dict(vOut=vout, r20=r20, weights=(ps * S).tolist(), psiSensorAbs=np.abs(ps).tolist())


def exact_end_gas_volume(depth, radii, h=2.5e-5):
    """V_out(r) of the exact chamber: midpoint rule on an h-grid over the planform, analytic chord per point."""
    R = BORE / 2
    p = POCKET
    xs = np.arange(min(-R, p["xMin"]) + 0.5 * h, R, h)
    zs = np.arange(-R + 0.5 * h, R, h)
    out = np.zeros(len(radii))
    vtot = 0.0
    hx = 0.5 * (p["xMax"] - p["xMin"]) - p["cornerRadius"]
    hz = 0.5 * (p["zMax"] - p["zMin"]) - p["cornerRadius"]
    for x in xs:
        z = zs
        inb = x * x + z * z <= R * R
        qx = abs(x - 0.5 * (p["xMin"] + p["xMax"])) - hx
        qz = np.abs(z - 0.5 * (p["zMin"] + p["zMax"])) - hz
        d = np.hypot(np.maximum(qx, 0), np.maximum(qz, 0)) + np.minimum(np.maximum(qx, qz), 0) - p["cornerRadius"]
        inp = (~inb) & (d <= 0)
        lo = np.where(inb, -depth, DECK_Y)
        hi = np.where(inb, 0.0, ROOF_Y)
        m = inb | inp
        lo, hi, zz = lo[m], hi[m], z[m]
        vtot += ((hi - lo).sum()) * h * h
        rho2 = (x - GAP[0]) ** 2 + (zz - GAP[2]) ** 2
        dlo, dhi = lo - GAP[1], hi - GAP[1]
        for i, r in enumerate(radii):
            w = np.sqrt(np.maximum(0.0, r * r - rho2))
            chord = np.clip(np.minimum(dhi, w) - np.maximum(dlo, -w), 0.0, None)
            out[i] += (hi - lo - chord).sum() * h * h
    return vtot, out.tolist()


def channel_roots(L1, L2, H1, H2, n=5):
    f = lambda k: H1 * math.sin(k * L1) * math.cos(k * L2) + H2 * math.cos(k * L1) * math.sin(k * L2)
    ks, k, step = [], 1e-6, 0.002 / (L1 + L2)
    prev = f(k)
    while len(ks) < n:
        k2 = k + step
        cur = f(k2)
        if prev * cur < 0:
            ks.append(optimize.brentq(f, k, k2, xtol=1e-14))
        k, prev = k2, cur
    return ks


# --------------------------------------------------------------------------------------------
# 3-D Helmholtz of the same chamber (no depth averaging) and the depth-averaged model on its staircase
# --------------------------------------------------------------------------------------------
def staircase(nper):
    """Planform staircase (region of each cell centre), z-symmetric grid of cell B/nper over the chamber's box."""
    h = BORE / nper
    R = BORE / 2
    x0, x1, z0, z1 = POCKET["xMin"], R, -R, R
    nx, nz = math.ceil((x1 - x0) / h), math.ceil((z1 - z0) / h)
    if nz % 2:
        nz += 1
    gx0, gz0 = 0.5 * (x0 + x1) - 0.5 * nx * h, -0.5 * nz * h
    X, Z = np.meshgrid(gx0 + (np.arange(nx) + 0.5) * h, gz0 + (np.arange(nz) + 0.5) * h, indexing="ij")
    p = POCKET
    hx = 0.5 * (p["xMax"] - p["xMin"]) - p["cornerRadius"]
    hz = 0.5 * (p["zMax"] - p["zMin"]) - p["cornerRadius"]
    qx = np.abs(X - 0.5 * (p["xMin"] + p["xMax"])) - hx
    qz = np.abs(Z - 0.5 * (p["zMin"] + p["zMax"])) - hz
    rr = np.hypot(np.maximum(qx, 0), np.maximum(qz, 0)) + np.minimum(np.maximum(qx, qz), 0) - p["cornerRadius"]
    reg = np.full(X.shape, -1, dtype=int)
    inb = X**2 + Z**2 <= R * R
    reg[(~inb) & (rr <= 0)] = 1
    reg[inb] = 0
    return reg, X, Z, h


def _laplacian(n, pairs, extra_diag=None):
    r = np.concatenate([a for a, _, _ in pairs])
    c = np.concatenate([b for _, b, _ in pairs])
    g = np.concatenate([w for _, _, w in pairs])
    diag = np.bincount(r, g, n) + np.bincount(c, g, n)
    if extra_diag is not None:
        diag += extra_diag
    ii = np.arange(n)
    return sparse.coo_matrix((np.concatenate([-g, -g, diag]), (np.concatenate([r, c, ii]), np.concatenate([c, r, ii]))), shape=(n, n)).tocsr()


def _alphas(K, M, nev):
    lam = eigsh(K, k=nev, M=M, sigma=-1.0, which="LM")[0]
    return np.sqrt(np.maximum(np.sort(lam), 0)) * BORE / 2


def staircase_longwave(reg, X, Z, h, depth, rwin):
    """Depth-averaged finite volumes on the staircase: harmonic-mean faces, r/|n.e| on bore <-> pocket faces."""
    act = reg >= 0
    idx = -np.ones(reg.shape, dtype=int)
    n = int(act.sum())
    idx[act] = np.arange(n)
    H = np.where(reg == 0, depth, np.where(reg == 1, ROOF_Y - DECK_Y, 0.0))
    pairs = []
    for (di, dk) in ((1, 0), (0, 1)):
        a = (slice(0, reg.shape[0] - di), slice(0, reg.shape[1] - dk))
        b = (slice(di, None), slice(dk, None))
        m = act[a] & act[b]
        res = h * (0.5 / H[a][m] + 0.5 / H[b][m])
        xf, zf = 0.5 * (X[a][m] + X[b][m]), 0.5 * (Z[a][m] + Z[b][m])
        ne = np.abs((xf if di else zf) / np.hypot(xf, zf))
        res = res + (reg[a][m] != reg[b][m]) * rwin / np.maximum(ne, 1e-12)
        pairs.append((idx[a][m], idx[b][m], h / res))
    K = _laplacian(n, pairs)
    return _alphas(K, sparse.diags((H * h * h)[act]), 7)[1:]


def _layers(segs, dymax, dymin, sing):
    """y-layer edges over the column intervals, spacing dymin + 0.3 d (<= dymax) at distance d from the window edges."""
    br = sorted({v for s in segs for v in s})
    edges = [br[0]]
    for a, b in zip(br[:-1], br[1:]):
        y = a
        while y < b - 1e-12:
            step = min(dymax, dymin + 0.3 * min(abs(y - s) for s in sing))
            nxt = min(b, y + step)
            if b - nxt < 0.3 * step:
                nxt = b
            edges.append(nxt)
            y = nxt
    return np.array(edges)


def voxel_helmholtz(reg, h, depth, dymax=1.0e-3, dymin=0.25e-3, nev=6):
    """Lowest alpha_j of the 3-D Neumann Helmholtz problem in the rigid flat-roofed columns (bore column over
    [-depth, 0], pocket over [deck, roof]: they meet only across the window). 7-point finite volumes on
    staircase cells x y-layers; z >= 0 half solved twice (even modes: no flux on z = 0; odd: psi = 0 there)."""
    reg = reg[:, reg.shape[1] // 2:]
    edges = _layers([(-depth, 0.0), (DECK_Y, ROOF_Y)], dymax, dymin, [-depth, ROOF_Y, DECK_Y])
    yc, dy = 0.5 * (edges[1:] + edges[:-1]), np.diff(edges)
    nx, nz, ny = reg.shape[0], reg.shape[1], len(yc)
    act = ((reg == 0)[:, :, None] & (yc > -depth) & (yc < 0.0)) | ((reg == 1)[:, :, None] & (yc > DECK_Y) & (yc < ROOF_Y))
    idx = -np.ones(act.shape, dtype=int)
    n = int(act.sum())
    idx[act] = np.arange(n)
    DY = np.broadcast_to(dy, (nx, nz, ny))
    pairs = []
    m = act[:-1] & act[1:]
    pairs.append((idx[:-1][m], idx[1:][m], DY[:-1][m]))
    m = act[:, :-1] & act[:, 1:]
    pairs.append((idx[:, :-1][m], idx[:, 1:][m], DY[:, :-1][m]))
    m = act[:, :, :-1] & act[:, :, 1:]
    pairs.append((idx[:, :, :-1][m], idx[:, :, 1:][m], np.broadcast_to(h * h / (0.5 * (dy[:-1] + dy[1:])), (nx, nz, ny - 1))[m]))
    M = sparse.diags((h * h * DY)[act])
    even = _alphas(_laplacian(n, pairs), M, nev + 1)[1:]
    mirror = np.zeros(n)
    m0 = act[:, 0, :]
    mirror[idx[:, 0, :][m0]] = 2 * DY[:, 0, :][m0]  # ghost psi = -psi across z = 0
    odd = _alphas(_laplacian(n, pairs, mirror), M, nev)
    return np.sort(np.concatenate([even, odd]))[:nev], n, ny


def main():
    depth = H_TDC
    sensor = (GAP[0], GAP[2])
    radii = [0.0, 0.02, 0.04, 0.06, 0.08, 0.09, 0.1, 0.105, 0.11]
    res = dict(geometry=dict(bore=BORE, deckY=DECK_Y, crownAboveDeckAtTDC=CROWN, pocket=dict(POCKET, roofY=ROOF_Y),
                             gapCenter=list(GAP), depthTDC=depth, chamberVolumeExact=VC - CREVICE))
    # window junction resistance: Ritz (knock.ts method), exact slits, 2-D finite-volume section solves
    t0 = time.time()
    H, a = H_TDC, 0.00498
    lay = lambda lo1, hi1, lo2, hi2: dict(lo1=lo1, hi1=hi1, lo2=lo2, hi2=hi2, wLo=max(lo1, lo2), wHi=min(hi1, hi2))
    cases = [dict(name="mirror slit at one wall", deg=None, **dict(lay(0.0, H, 0.0, H), wLo=H - a), exact=slit_exact(H, H - a, H)),
             dict(name="mirror slit inside", deg=None, **dict(lay(0.0, H, 0.0, H), wLo=0.4 * H, wHi=0.4 * H + a),
                  exact=slit_exact(H, 0.4 * H, 0.4 * H + a)),
             dict(name="ordinary step at the roof", deg=None, **lay(-H_TDC, 0.0, -POCKET_HEIGHT, 0.0), exact=None)]
    for deg in (0, 15, 25, 30, 40):
        d = depth_at(deg)
        cases.append(dict(name=f"Model T {deg} deg ATDC", deg=deg, **lay(-d, 0.0, DECK_Y, ROOF_Y), exact=None))
    for c in cases:
        w = (c["wLo"], c["wHi"])
        c["ritz"] = window_r(c["lo1"], c["hi1"], c["lo2"], c["hi2"], w)
        c["ritzFine"] = window_r(c["lo1"], c["hi1"], c["lo2"], c["hi2"], w, N=24, Q=1024)
        c["fv"] = section_fv(c["lo1"], c["hi1"], c["lo2"], c["hi2"], w)
        print(f"window {c['name']:28s}: Ritz {c['ritz']:.6f} (N 24: {c['ritzFine']:.6f})  FV {c['fv']:.6f}" +
              (f"  exact {c['exact']:.6f}" if c["exact"] is not None else ""))
        assert abs(c["fv"] / c["ritz"] - 1) < 3e-3
    res["window"] = dict(ritzBasis=RITZ_N, ritzNodes=RITZ_Q, fvFinest=1.25e-5, cases=cases)
    print(f"  ({time.time() - t0:.0f} s)")
    # same grid as the TS default
    g = build(LHead(BORE, depth), BORE / 64)
    alpha, psi, ms = modes(g, BORE)
    pc = projector_checks(g, alpha, psi, ms, GAP, sensor, radii)
    res["sameGrid"] = dict(cellSize=BORE / 64, nCells=g["n"], nColumns=len(g["cols"]), volume=g["C"].sum(),
                           windowResistance=g["pf"].rwin, alpha=alpha.tolist(), meanSquare=ms.tolist(), radii=radii, **pc)
    print("same grid:", g["n"], "cells; alpha", np.round(alpha[:8], 5), "V", g["C"].sum())
    g0 = build(LHead(BORE, depth, window=False), BORE / 64)
    a0, _, _ = modes(g0, BORE)
    res["sameGridNoWindow"] = dict(cellSize=BORE / 64, alpha=a0.tolist())
    print("same grid, no window: alpha", np.round(a0[:8], 5))
    conv = []
    for nper in (128, 192):
        gf = build(LHead(BORE, depth), BORE / nper)
        af, pf_, mf = modes(gf, BORE)
        pcf = projector_checks(gf, af, pf_, mf, GAP, sensor, radii)
        conv.append(dict(cellsPerBore=nper, nCells=gf["n"], volume=gf["C"].sum(), alpha=af.tolist(), **pcf))
        print("B/%d:" % nper, gf["n"], "cells; alpha", np.round(af[:8], 5))
    res["converged"] = conv
    # 15 deg ATDC bore-column depth (slider-crank, L = 7 in, r = 2 in)
    d15 = depth_at(15)
    g15 = build(LHead(BORE, d15), BORE / 96)
    a15, _, _ = modes(g15, BORE)
    res["depth15"] = dict(depth=d15, cellsPerBore=96, alpha=a15.tolist())
    print("15 deg ATDC B/96: alpha", np.round(a15[:8], 5))
    # 3-D Helmholtz (no depth averaging) vs the depth-averaged model on the same staircase
    three = []
    for nper in (64, 96):
        reg, X, Z, h = staircase(nper)
        for deg in (0, 15):
            t0 = time.time()
            d = depth_at(deg)
            rw = window_r(-d, 0.0, DECK_Y, ROOF_Y)
            a3, n3, ny = voxel_helmholtz(reg, h, d)
            lw = staircase_longwave(reg, X, Z, h, d, 0.0)[:6]
            lwr = staircase_longwave(reg, X, Z, h, d, rw)[:6]
            three.append(dict(deg=deg, depth=d, window=ROOF_Y - max(DECK_Y, -d), windowResistance=rw, cellsPerBore=nper,
                              dyMax=1.0e-3, dyMin=0.25e-3, layers=ny, voxelsHalf=n3, alpha3d=a3.tolist(),
                              alphaLongWave=lw.tolist(), alphaLongWaveWindow=lwr.tolist()))
            print(f"3-D B/{nper} {deg:2d} deg ATDC ({n3} voxels per half, {ny} layers, {time.time() - t0:.0f} s): alpha {np.round(a3, 4)}\n"
                  f"   long wave + window {np.round(lwr, 4)} ({np.round(100 * (lwr / a3 - 1), 2)} %); no window {np.round(lw, 4)} ({np.round(100 * (lw / a3 - 1), 2)} %)")
            assert np.all(np.abs(lwr[:5] / a3[:5] - 1) < 0.01)
    res["threeD"] = three
    # disc (CFR bore), same grid as the TS default
    Bc = 0.08255
    gd = build(Disc(Bc / 2, 0.01), Bc / 64)
    ad, _, _ = modes(gd, Bc)
    ref = sorted([z for m in range(0, 9) for z in special.jnp_zeros(m, 4) if z <= 7.1 for _ in range(1 if m == 0 else 2)])
    res["disc"] = dict(bore=Bc, depth=0.01, cellSize=Bc / 64, alpha=ad.tolist(), besselReference=ref)
    print("disc err %", np.round(100 * (ad / np.array(ref[: len(ad)]) - 1), 3))
    # two-depth channel
    L1, L2, H1, H2, W = 0.0631, 0.0457, 0.0254, 0.0129, 0.006
    gc = build(Channel(L1, L2, H1, H2, W), 0.0015)
    ac, _, _ = modes(gc, 0.1, max_alpha=8.0, kreq=12)
    res["channel"] = dict(L1=L1, L2=L2, H1=H1, H2=H2, width=W, cellSize=0.0015, kExact=channel_roots(L1, L2, H1, H2),
                          kSameGrid=(2 * ac / 0.1).tolist())
    print("channel k exact", np.round(res["channel"]["kExact"], 4), "grid", np.round(res["channel"]["kSameGrid"][:5], 4))
    # exact end-gas volume
    vt, vo = exact_end_gas_volume(depth, radii)
    res["endGasVolume"] = dict(radii=radii, volume=vt, vOut=vo)
    print("exact V", vt, "vs Vc - crevice", VC - CREVICE, "; V_out", np.round(np.array(vo) / vt, 5))
    print("grid  V_out", np.round(np.array(pc["vOut"]) / g["C"].sum(), 5))
    write_json(os.path.join(FIXTURES, "chemistry_knock_lhead_modes.json"), res)


if __name__ == "__main__":
    main()
