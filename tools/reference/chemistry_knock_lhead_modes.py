"""Oracle for the L-head knock acoustics of src/physics/chemistry/knock.ts (GridModeSet, GridEndGasProjector).

The Ford Model T chamber (src/physics/engines/model-t.ts): bore column (disc of the bore, roof y = 0 down to
the crown at TDC) + valve pocket (rounded rectangle minus the bore disc, deck to pocket roof). Transverse
modes of the depth-weighted Helmholtz problem div(H grad psi) + k^2 H psi = 0, H dpsi/dn = 0 on the outline.

Writes test/fixtures/chemistry_knock_lhead_modes.json:
  geometry     : the Model T numbers, recomputed here with the expressions of model-t.ts
  sameGrid     : the TS discretisation (cut-cell finite volumes, cell = B/64) rebuilt independently here and
                 solved with scipy.sparse.linalg.eigsh (ARPACK shift-invert, not the TS block Lanczos):
                 alpha_j = k_j B/2, cell / sub-column counts, chamber volume, the end-gas volume outside
                 a few ball radii and the sign-invariant modal weights psi_j(sensor) * S_j (S_j = <psi_j>_eg /
                 <psi_j^2>) of the gas outside the ball about the plug that leaves 20 % end gas
  converged    : the same problem on finer grids (B/128, B/192): grid convergence of alpha_j and the weights
  depth15      : alpha_j at the bore-column depth of 15 deg ATDC (B/96)
  disc         : same-grid alpha of a flat disc of the CFR bore (82.55 mm), scipy jnp_zeros for reference
  channel      : two-depth channel (L1, H1 | L2, H2): analytic k from H1 tan(k L1) + H2 tan(k L2) = 0 (brentq)
  endGasVolume : V_out(r) of the EXACT chamber geometry (ball about the plug), midpoint quadrature on a
                 0.025 mm planform grid with the analytic chord -- independent of the cut-cell construction
Run: .venv/bin/python tools/reference/chemistry_knock_lhead_modes.py   (about 1-2 min)
"""
from __future__ import annotations

import math
import os

import numpy as np
from scipy import optimize, sparse, special
from scipy.sparse.linalg import eigsh

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
POCKET_HEIGHT = (VC - CREVICE - A_PISTON * H_TDC) / 46.016e-4
ROOF_Y = DECK_Y + POCKET_HEIGHT
GAP = (-0.06, DECK_Y + POCKET_HEIGHT - 0.004, 0.0)
ROD = 7.0 * IN


def rr_dist(x, z, p):
    hx = 0.5 * (p["xMax"] - p["xMin"]) - p["cornerRadius"]
    hz = 0.5 * (p["zMax"] - p["zMin"]) - p["cornerRadius"]
    qx = abs(x - 0.5 * (p["xMin"] + p["xMax"])) - hx
    qz = abs(z - 0.5 * (p["zMin"] + p["zMax"])) - hz
    return math.hypot(max(qx, 0.0), max(qz, 0.0)) + min(max(qx, qz), 0.0) - p["cornerRadius"]


class LHead:
    def __init__(self, bore, depth):
        self.R = bore / 2
        R = self.R
        p = POCKET
        self.xMin, self.xMax = min(-R, p["xMin"]), max(R, p["xMax"])
        self.zMin, self.zMax = min(-R, p["zMin"]), max(R, p["zMax"])
        self.yLo = [-depth, DECK_Y]
        self.yHi = [0.0, ROOF_Y]

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
    caps, cxs, czs = [], [], []
    cols = []  # (x, z, area, yLo, yHi, cell)
    for k in range(nz):
        for i in range(nx):
            x = gx0 + (i + 0.5) * dx
            z = gz0 + (k + 0.5) * dx
            mine = []
            if pf.bdist(x, z) > dUni:
                r = pf.region(x, z)
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
            cap = sum(a * (hi - lo) for (_, _, a, lo, hi) in mine)
            if cap > 0:
                c = len(caps)
                cell_index[i, k] = c
                caps.append(cap)
                cxs.append(x)
                czs.append(z)
                for m in mine:
                    cols.append(m + (c,))

    def face(xf, zf, ex, ez):
        if pf.bdist(xf, zf) > dUni:
            return hof(pf.region(xf, zf))
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
            tot += n_in / inv if n_in > 0 else hf
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
                G = face(cxs[c] + 0.5 * dx * di, czs[c] + 0.5 * dx * dk, di, dk)
                if not G > 0:
                    continue
                d = cell_index[ii, kk]
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


def main():
    depth = H_TDC
    sensor = (GAP[0], GAP[2])
    radii = [0.0, 0.02, 0.04, 0.06, 0.08, 0.09, 0.1, 0.105, 0.11]
    res = dict(geometry=dict(bore=BORE, deckY=DECK_Y, crownAboveDeckAtTDC=CROWN, pocket=dict(POCKET, roofY=ROOF_Y),
                             gapCenter=list(GAP), depthTDC=depth, chamberVolumeExact=VC - CREVICE))
    # same grid as the TS default
    g = build(LHead(BORE, depth), BORE / 64)
    alpha, psi, ms = modes(g, BORE)
    pc = projector_checks(g, alpha, psi, ms, GAP, sensor, radii)
    res["sameGrid"] = dict(cellSize=BORE / 64, nCells=g["n"], nColumns=len(g["cols"]), volume=g["C"].sum(),
                           alpha=alpha.tolist(), meanSquare=ms.tolist(), radii=radii, **pc)
    print("same grid:", g["n"], "cells; alpha", np.round(alpha[:8], 5), "V", g["C"].sum())
    conv = []
    for nper in (128, 192):
        gf = build(LHead(BORE, depth), BORE / nper)
        af, pf_, mf = modes(gf, BORE)
        pcf = projector_checks(gf, af, pf_, mf, GAP, sensor, radii)
        conv.append(dict(cellsPerBore=nper, nCells=gf["n"], volume=gf["C"].sum(), alpha=af.tolist(), **pcf))
        print("B/%d:" % nper, gf["n"], "cells; alpha", np.round(af[:8], 5))
    res["converged"] = conv
    # 15 deg ATDC bore-column depth (slider-crank, L = 7 in, r = 2 in)
    th = math.radians(15)
    rc = STROKE / 2
    x15 = rc * (1 - math.cos(th)) + ROD - math.sqrt(ROD**2 - (rc * math.sin(th)) ** 2)
    g15 = build(LHead(BORE, depth + x15), BORE / 96)
    a15, _, _ = modes(g15, BORE)
    res["depth15"] = dict(depth=depth + x15, cellsPerBore=96, alpha=a15.tolist())
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
