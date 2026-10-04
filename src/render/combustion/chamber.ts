/**
 * Chamber SHAPE of the in-cylinder visuals (pure geometry, no three.js), cylinder frame (engine-spec.ts):
 *
 *  - 'flat-disc' (CFR pancake): x²+z² ≤ R², −h ≤ y ≤ 0;
 *  - 'l-head' (side valves, LHeadChamberSpec): the union of
 *      the BORE COLUMN  x²+z² ≤ R², −h ≤ y ≤ 0 (crown → roof of the head cavity over the bore), and
 *      the VALVE POCKET (rounded-rectangle plan Ω ∖ bore disc) × [deckY, roofY] (valves seat in its floor).
 *
 * h is the depth of the piston crown below y = 0 (chamberDepth). The flame drawn by the gas volume is
 * sphere(flame.center, flame.radius) ∩ this shape — the same set the physics' entrainment model burns
 * (DESIGN.md: "a sphere centred at the spark gap intersected with the chamber").
 *
 * The GLSL in volume-shader.ts mirrors `rayChamberIntervals` (ray ∩ shape = at most two disjoint
 * intervals, each tagged with the proxy that owns its exit: EXIT_*); chamber.test.ts checks the TS
 * versions against brute-force sampling. The disc is the special case without a pocket and keeps the
 * original code paths everywhere (geometry.ts rayChamberInterval / insideChamber).
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { V3 } from './geometry';

/** Exit of a gas interval through the bore column's walls (liner, head roof over the bore, crown). */
export const EXIT_BORE = 0;
/** Exit through a valve-pocket wall, its roof or the deck/seat plane (outside the bore disc). */
export const EXIT_POCKET = 1;
/** Exit from the pocket INTO the bore circle below the crown (the piston's top land above the deck). */
export const EXIT_TRANSFER = 2;

/** Tracer/light wall margin: fraction of R kept clear of the liner (as the disc's 0.995 R clamp). */
export const WALL_MARGIN_FRACTION = 0.005;

export interface ChamberShape {
  readonly kind: 'flat-disc' | 'l-head';
  /** Bore radius, m. */
  readonly R: number;
  /** L-head pocket plan: inner rectangle (corner-arc centres) [ix0, ix1] × [iz0, iz1] and corner radius rc, m. */
  readonly ix0: number;
  readonly ix1: number;
  readonly iz0: number;
  readonly iz1: number;
  readonly rc: number;
  /** Deck (valve-seat) plane and pocket roof, cylinder-frame y, m. */
  readonly deckY: number;
  readonly roofY: number;
  /** L-head: crown depth below y = 0 at TDC = −deckY − crownAboveDeckAtTDC, m (disc: NaN, CR-dependent). */
  readonly depthTDC: number;
  /** Plan area of pocket ∖ bore disc, m², its volume (× pocket height), m³, and its plan centroid (x, z). */
  readonly pocketPlanArea: number;
  readonly pocketVolume: number;
  readonly pocketCentroid: readonly [number, number];
  /**
   * Footprint (bore disc ∪ pocket plan) axis: unit vector (ex, ez) from the bore axis toward the pocket
   * centroid, and the footprint's extent [s0, s1] along it (s = x·ex + z·ez), m. Disc: (1, 0), [−R, R].
   */
  readonly axis: readonly [number, number];
  readonly s0: number;
  readonly s1: number;
  /** Footprint width across the axis (pocket width; disc: 2R), m. */
  readonly width: number;
  /** Radius of the bounding circle of the footprint about the bore axis, m. */
  readonly planRadius: number;
  /**
   * Plan angles (atan2(z, x), rad) of the bore-circle arc that lies inside the pocket plan — the
   * pocket/bore transfer boundary [a0, a1] (a1 ≥ a0; empty: a1 < a0).
   */
  readonly transferArc: readonly [number, number];
}

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

/** The chamber shape of a spec ('l-head' requires geometry.lHead). */
export function chamberShapeOf(spec: EngineSpec): ChamberShape {
  const R = spec.geometry.bore / 2;
  const lh = spec.geometry.lHead;
  if (spec.geometry.chamber !== 'l-head' || !lh) {
    return {
      kind: 'flat-disc', R, ix0: 0, ix1: 0, iz0: 0, iz1: 0, rc: 0, deckY: 0, roofY: 0, depthTDC: NaN,
      pocketPlanArea: 0, pocketVolume: 0, pocketCentroid: [0, 0], axis: [1, 0], s0: -R, s1: R, width: 2 * R,
      planRadius: R, transferArc: [0, -1],
    };
  }
  const p = lh.pocket;
  const rc = Math.max(0, Math.min(p.cornerRadius, (p.xMax - p.xMin) / 2, (p.zMax - p.zMin) / 2));
  const base = {
    kind: 'l-head' as const, R,
    ix0: p.xMin + rc, ix1: p.xMax - rc, iz0: p.zMin + rc, iz1: p.zMax - rc, rc,
    deckY: lh.deckY, roofY: p.roofY, depthTDC: -lh.deckY - lh.crownAboveDeckAtTDC,
  };
  const plan = pocketPlanIntegrals(base);
  const height = Math.max(p.roofY - lh.deckY, 0);
  // axis: bore axis → pocket centroid (falls back to the rectangle centre, then +x)
  let ex = plan.cx, ez = plan.cz;
  if (!(Math.hypot(ex, ez) > 1e-9)) { ex = (p.xMin + p.xMax) / 2; ez = (p.zMin + p.zMax) / 2; }
  const n = Math.hypot(ex, ez);
  if (n > 1e-12) { ex /= n; ez /= n; } else { ex = 1; ez = 0; }
  const s1 = Math.max(R, roundRectSupport(base, ex, ez));
  const s0 = -Math.max(R, roundRectSupport(base, -ex, -ez));
  const width = roundRectSupport(base, -ez, ex) + roundRectSupport(base, ez, -ex);
  let planRadius = R;
  for (const [qx, qz] of [[base.ix0, base.iz0], [base.ix1, base.iz0], [base.ix0, base.iz1], [base.ix1, base.iz1]]) {
    planRadius = Math.max(planRadius, Math.hypot(qx, qz) + rc);
  }
  const shape: ChamberShape = {
    ...base,
    pocketPlanArea: plan.area,
    pocketVolume: plan.area * height,
    pocketCentroid: [plan.cx, plan.cz],
    axis: [ex, ez], s0, s1, width, planRadius,
    transferArc: [0, -1],
  };
  return { ...shape, transferArc: transferArcOf(shape) };
}

type RoundRect = Pick<ChamberShape, 'ix0' | 'ix1' | 'iz0' | 'iz1' | 'rc'>;

/** Support function of the rounded rectangle: max over its points of x·dx + z·dz (|d| = 1). */
function roundRectSupport(s: RoundRect, dx: number, dz: number): number {
  return Math.max(s.ix0 * dx, s.ix1 * dx) + Math.max(s.iz0 * dz, s.iz1 * dz) + s.rc;
}

/** z-interval of the rounded rectangle at abscissa x ([1, −1] if none). */
function roundRectSliceZ(s: RoundRect, x: number, out: [number, number]): [number, number] {
  const dx = x < s.ix0 ? s.ix0 - x : x > s.ix1 ? x - s.ix1 : 0;
  if (dx > s.rc) { out[0] = 1; out[1] = -1; return out; }
  const w = dx > 0 ? Math.sqrt(Math.max(s.rc * s.rc - dx * dx, 0)) : s.rc;
  out[0] = s.iz0 - w;
  out[1] = s.iz1 + w;
  return out;
}

/** 8-point Gauss–Legendre nodes/weights on [−1, 1] (Abramowitz & Stegun 1964, Table 25.4). */
const GL_X = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
  0.1834346424956498, 0.5255324099163290, 0.7966664774136267, 0.9602898564975363];
const GL_W = [0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
  0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763];

/**
 * Area and centroid of (rounded rectangle ∖ bore disc) by composite Gauss–Legendre quadrature over x of
 * the exact z-slice lengths, with panel breaks at every abscissa where a slice function has a kink
 * (corner-arc starts, x = ±R). Construction-time only (allocates).
 */
function pocketPlanIntegrals(s: RoundRect & { R: number }): { area: number; cx: number; cz: number } {
  const xa = s.ix0 - s.rc, xb = s.ix1 + s.rc;
  const brk = [xa, s.ix0, s.ix1, xb, -s.R, s.R].filter((x) => x >= xa && x <= xb).sort((a, b) => a - b);
  const sl: [number, number] = [0, 0];
  let A = 0, Mx = 0, Mz = 0;
  const PANELS = 600;
  for (let k = 0; k + 1 < brk.length; k++) {
    const a = brk[k], b = brk[k + 1];
    if (!(b > a)) continue;
    const hw = (b - a) / PANELS;
    for (let p = 0; p < PANELS; p++) {
      const m = a + (p + 0.5) * hw;
      for (let g = 0; g < 8; g++) {
        const x = m + 0.5 * hw * GL_X[g];
        const w = 0.5 * hw * GL_W[g];
        roundRectSliceZ(s, x, sl);
        if (!(sl[1] > sl[0])) continue;
        let len = sl[1] - sl[0];
        let mz = 0.5 * (sl[1] * sl[1] - sl[0] * sl[0]);
        if (Math.abs(x) < s.R) {
          const q = Math.sqrt(s.R * s.R - x * x);
          const lo = Math.max(sl[0], -q), hi = Math.min(sl[1], q);
          if (hi > lo) { len -= hi - lo; mz -= 0.5 * (hi * hi - lo * lo); }
        }
        A += w * len;
        Mx += w * x * len;
        Mz += w * mz;
      }
    }
  }
  return { area: A, cx: A > 0 ? Mx / A : 0, cz: A > 0 ? Mz / A : 0 };
}

/** Plan-angle range of the bore circle inside the pocket plan (bisection on the boundary crossings). */
function transferArcOf(s: ChamberShape): [number, number] {
  const inside = (a: number): boolean => insideRoundRect(s.R * Math.cos(a), s.R * Math.sin(a), s);
  const ac = Math.atan2(s.axis[1], s.axis[0]);
  if (!inside(ac)) return [0, -1];
  const edge = (dir: number): number => {
    let lo = 0, hi = Math.PI;
    if (inside(ac + dir * hi)) return dir * Math.PI;
    for (let i = 0; i < 60; i++) {
      const m = 0.5 * (lo + hi);
      if (inside(ac + dir * m)) lo = m; else hi = m;
    }
    return dir * lo;
  };
  return [ac + edge(-1), ac + edge(1)];
}

// ---------------------------------------------------------------------------
// Crown position and volumes
// ---------------------------------------------------------------------------

/**
 * Depth h of the piston crown below the cylinder-frame y = 0 (crown at y = −h), m.
 *  - flat disc: the snapshot's clearanceHeight (= h_TDC(CR) + x(θ), kinematics.ts; the CR is an
 *    operating-point variable, so h_TDC is not a spec constant);
 *  - L-head: from the spec and the piston displacement, h = −deckY − crownAboveDeckAtTDC + x(θ)
 *    (Model T: crown 7.94 mm above the deck, 25.4 mm below the roof at TDC). Independent of how the
 *    physics defines clearanceHeight; with the engine-spec frame convention the two agree (the bore-
 *    column depth from the roof).
 */
export function chamberDepth(s: ChamberShape, clearanceHeight: number, pistonDisplacement: number): number {
  return s.kind === 'l-head' ? s.depthTDC + pistonDisplacement : clearanceHeight;
}

/** Gas volume of the chamber at crown depth h (bore column + pocket; crevices excluded), m³. */
export function chamberVolume(s: ChamberShape, h: number): number {
  return Math.PI * s.R * s.R * Math.max(h, 0) + s.pocketVolume;
}

/**
 * Mean chamber depth V/A_plan, m: h for the disc; for the L-head the chamber volume over the footprint
 * area (the Poulos–Heywood "4V/(πB²)"-type length generalised to the planform).
 */
export function effectiveDepth(s: ChamberShape, h: number): number {
  if (s.kind !== 'l-head') return h;
  return chamberVolume(s, h) / (Math.PI * s.R * s.R + s.pocketPlanArea);
}

/** A distance beyond which a ray origin is surely outside the chamber (ortho back-off), m. */
export function boundingSize(s: ChamberShape, h: number): number {
  return s.planRadius + Math.max(h, s.kind === 'l-head' ? -s.deckY : 0);
}

// ---------------------------------------------------------------------------
// Inside tests
// ---------------------------------------------------------------------------

/** Rounded-rectangle plan test (Minkowski sum of the inner rectangle and a disc of radius rc). */
export function insideRoundRect(x: number, z: number, s: RoundRect, inset = 0): boolean {
  const dx = Math.max(s.ix0 - x, 0, x - s.ix1);
  const dz = Math.max(s.iz0 - z, 0, z - s.iz1);
  const r = s.rc - inset;
  if (r >= 0) return dx * dx + dz * dz <= r * r;
  // inset deeper than the corner radius: the inner rectangle shrunk by (inset − rc)
  return dx === 0 && dz === 0 && Math.min(x - s.ix0, s.ix1 - x, z - s.iz0, s.iz1 - z) >= -r;
}

/** True if (x, y, z) is in the chamber gas at crown depth h (boundaries inclusive). */
export function insideChamberShape(x: number, y: number, z: number, s: ChamberShape, h: number): boolean {
  if (x * x + z * z <= s.R * s.R) return y <= 0 && y >= -h;
  return s.kind === 'l-head' && y >= s.deckY && y <= s.roofY && insideRoundRect(x, z, s);
}

/** Burned region = flame sphere ∩ chamber shape. */
export function insideBurnedShape(
  x: number, y: number, z: number, s: ChamberShape, h: number, c: ArrayLike<number>, r: number,
): boolean {
  if (!(r > 0)) return false;
  const dx = x - c[0], dy = y - c[1], dz = z - c[2];
  return dx * dx + dy * dy + dz * dz <= r * r && insideChamberShape(x, y, z, s, h);
}

// ---------------------------------------------------------------------------
// Ray intersections (mirrored in GLSL)
// ---------------------------------------------------------------------------

const BIG = 1e30;

/** Ray (o + t d) ∩ y-slab [y0, y1] → [t0, t1] in out (empty: t0 > t1). */
function raySlabY(oy: number, dy: number, y0: number, y1: number, out: number[]): void {
  if (Math.abs(dy) < 1e-18) {
    if (oy < y0 || oy > y1) { out[0] = BIG; out[1] = -BIG; } else { out[0] = -BIG; out[1] = BIG; }
    return;
  }
  const a = (y0 - oy) / dy, b = (y1 - oy) / dy;
  out[0] = Math.min(a, b);
  out[1] = Math.max(a, b);
}

/** 2-D ray ∩ circle (centre (cx, cz), radius r) in the xz plane → [t0, t1] (empty: t0 > t1). */
function rayCircle2(ox: number, oz: number, dx: number, dz: number, cx: number, cz: number, r: number, out: number[]): void {
  const qx = ox - cx, qz = oz - cz;
  const a = dx * dx + dz * dz;
  const c = qx * qx + qz * qz - r * r;
  if (a < 1e-18) {
    if (c > 0) { out[0] = BIG; out[1] = -BIG; } else { out[0] = -BIG; out[1] = BIG; }
    return;
  }
  const b = 2 * (qx * dx + qz * dz);
  const disc = b * b - 4 * a * c;
  if (disc < 0) { out[0] = BIG; out[1] = -BIG; return; }
  const sq = Math.sqrt(disc);
  out[0] = (-b - sq) / (2 * a);
  out[1] = (-b + sq) / (2 * a);
}

/** 2-D ray ∩ axis-aligned box [x0, x1] × [z0, z1] → [t0, t1] (empty: t0 > t1). */
function rayBox2(ox: number, oz: number, dx: number, dz: number, x0: number, x1: number, z0: number, z1: number, out: number[]): void {
  let t0 = -BIG, t1 = BIG;
  if (Math.abs(dx) < 1e-18) {
    if (ox < x0 || ox > x1) { out[0] = BIG; out[1] = -BIG; return; }
  } else {
    const a = (x0 - ox) / dx, b = (x1 - ox) / dx;
    t0 = Math.max(t0, Math.min(a, b));
    t1 = Math.min(t1, Math.max(a, b));
  }
  if (Math.abs(dz) < 1e-18) {
    if (oz < z0 || oz > z1) { out[0] = BIG; out[1] = -BIG; return; }
  } else {
    const a = (z0 - oz) / dz, b = (z1 - oz) / dz;
    t0 = Math.max(t0, Math.min(a, b));
    t1 = Math.min(t1, Math.max(a, b));
  }
  out[0] = t0;
  out[1] = t1;
}

const sTmp = [0, 0];

/**
 * 2-D ray ∩ rounded rectangle (convex) → [t0, t1] in out (empty: t0 > t1). The rounded rectangle is the
 * union of two crossed boxes and four corner discs; being convex, its chord is the hull of the pieces'.
 */
export function rayRoundRect(ox: number, oz: number, dx: number, dz: number, s: RoundRect, out: number[]): void {
  let t0 = BIG, t1 = -BIG;
  const r = s.rc;
  rayBox2(ox, oz, dx, dz, s.ix0 - r, s.ix1 + r, s.iz0, s.iz1, sTmp);
  if (sTmp[0] <= sTmp[1]) { t0 = Math.min(t0, sTmp[0]); t1 = Math.max(t1, sTmp[1]); }
  rayBox2(ox, oz, dx, dz, s.ix0, s.ix1, s.iz0 - r, s.iz1 + r, sTmp);
  if (sTmp[0] <= sTmp[1]) { t0 = Math.min(t0, sTmp[0]); t1 = Math.max(t1, sTmp[1]); }
  if (r > 0) {
    for (let k = 0; k < 4; k++) {
      rayCircle2(ox, oz, dx, dz, k & 1 ? s.ix1 : s.ix0, k & 2 ? s.iz1 : s.iz0, r, sTmp);
      if (sTmp[0] <= sTmp[1]) { t0 = Math.min(t0, sTmp[0]); t1 = Math.max(t1, sTmp[1]); }
    }
  }
  out[0] = t0;
  out[1] = t1;
}

const rC = [0, 0], rS = [0, 0], rP = [0, 0];

/** Append [t0, t1, kind] to the sorted interval list (merging a touching/overlapping piece); returns the count. */
function pushInterval(out: number[] | Float64Array, n: number, t0: number, t1: number, kind: number): number {
  if (!(t1 > t0)) return n;
  if (n > 0 && t0 <= out[3 * n - 2]) {
    if (t1 > out[3 * n - 2]) { out[3 * n - 2] = t1; out[3 * n - 1] = kind; }
    return n;
  }
  if (n >= 3) return n;
  out[3 * n] = t0; out[3 * n + 1] = t1; out[3 * n + 2] = kind;
  return n + 1;
}

/**
 * Ray ∩ chamber shape at crown depth h, t ≥ 0. Writes up to two disjoint, sorted intervals as
 * out = [t0, t1, exit, t0', t1', exit'] and returns their count. `exit` (EXIT_*) names the proxy face
 * the interval leaves through. `rd` need not be normalised (t in units of |rd|). Allocation-free.
 *
 * Construction: C = the bore's infinite cylinder, A = C ∩ [−h, 0] (bore column), P = pocket plan prism
 * ∩ [deckY, roofY]. The pocket piece P ∖ C splits into L = [p0, min(p1, c0)] before the bore and
 * Rt = [max(p0, c1), p1] after it, so the pieces arrive already sorted: L, A, Rt; touching pieces merge.
 * (At most two survive: L and Rt both non-empty need the whole chord of C inside the pocket slab, and
 * then A — the part of that chord above the crown — touches c0 or c1, y being monotone along the ray.)
 */
export function rayChamberIntervals(
  ro: ArrayLike<number>, rd: ArrayLike<number>, s: ChamberShape, h: number, out: number[] | Float64Array,
): number {
  let n = 0;
  rayCircle2(ro[0], ro[2], rd[0], rd[2], 0, 0, s.R, rC);
  const cHit = rC[0] <= rC[1];
  const c0 = rC[0], c1 = rC[1];
  let a0 = BIG, a1 = -BIG;
  if (cHit) {
    raySlabY(ro[1], rd[1], -h, 0, rS);
    a0 = Math.max(c0, rS[0]);
    a1 = Math.min(c1, rS[1]);
  }
  if (s.kind === 'l-head') {
    rayRoundRect(ro[0], ro[2], rd[0], rd[2], s, rP);
    raySlabY(ro[1], rd[1], s.deckY, s.roofY, rS);
    const p0 = Math.max(rP[0], rS[0]), p1 = Math.min(rP[1], rS[1]);
    if (p1 > p0) {
      if (!cHit) n = pushInterval(out, n, p0, p1, EXIT_POCKET);
      else n = pushInterval(out, n, p0, Math.min(p1, c0), c0 < p1 ? EXIT_TRANSFER : EXIT_POCKET);
    }
    n = pushInterval(out, n, a0, a1, EXIT_BORE);
    if (p1 > p0 && cHit) n = pushInterval(out, n, Math.max(p0, c1), p1, EXIT_POCKET);
  } else {
    n = pushInterval(out, n, a0, a1, EXIT_BORE);
  }
  // keep t ≥ 0
  let m = 0;
  for (let i = 0; i < n; i++) {
    const t1 = out[3 * i + 1];
    if (!(t1 > 0)) continue;
    out[3 * m] = Math.max(out[3 * i], 0);
    out[3 * m + 1] = t1;
    out[3 * m + 2] = out[3 * i + 2];
    m++;
  }
  return m;
}

// ---------------------------------------------------------------------------
// Nearest interior point (tracer walls, light clamp) and footprint extremes
// ---------------------------------------------------------------------------

const qTmp: V3 = [0, 0, 0];

/**
 * Nearest point of the chamber kept `margin`·R clear of the walls in plan (liner r ≤ (1 − m)R, pocket
 * plan inset by m·R, piston top land r ≥ (1 + m)R; the planes are not inset), written to out. Where the
 * bore column and the pocket are open to each other (pocket slab above the crown) the two pieces meet
 * seamlessly at r = (1 − m)R. Returns the squared distance moved (0: the point was inside and is copied
 * unchanged). Allocation-free.
 */
export function projectIntoChamber(
  x: number, y: number, z: number, s: ChamberShape, h: number, out: V3, margin = WALL_MARGIN_FRACTION,
): number {
  const Rb = (1 - margin) * s.R;
  const Rp = (1 + margin) * s.R;
  const r2 = x * x + z * z;
  const inBore = r2 <= Rb * Rb && y <= 0 && y >= -h;
  const lh = s.kind === 'l-head';
  const inset = margin * s.R;
  if (inBore || (lh && y >= s.deckY && y <= s.roofY && (r2 >= Rp * Rp || (r2 >= Rb * Rb && y >= -h))
    && insideRoundRect(x, z, s, inset))) {
    out[0] = x; out[1] = y; out[2] = z;
    return 0;
  }
  // candidate 1: bore column
  const r = Math.sqrt(r2);
  const f = r > Rb ? Rb / r : 1;
  out[0] = x * f;
  out[1] = Math.min(Math.max(y, -Math.max(h, 0)), 0);
  out[2] = z * f;
  let best = (out[0] - x) ** 2 + (out[1] - y) ** 2 + (out[2] - z) ** 2;
  if (!lh) return best;
  // candidate 2: pocket (plan point pulled into the inset rounded rectangle, pushed out of the bore disc)
  const py = Math.min(Math.max(y, s.deckY), s.roofY);
  const rMin = py >= -h ? Rb : Rp;
  const rin = Math.max(s.rc - inset, 0);
  const shrink = Math.max(inset - s.rc, 0);
  let px = Math.min(Math.max(x, s.ix0 + shrink), s.ix1 - shrink);
  let pz = Math.min(Math.max(z, s.iz0 + shrink), s.iz1 - shrink);
  const ddx = x - px, ddz = z - pz;
  const dd = Math.hypot(ddx, ddz);
  if (dd > rin) { px += (ddx * rin) / dd; pz += (ddz * rin) / dd; } else { px = x; pz = z; }
  const pr = Math.hypot(px, pz);
  if (pr < rMin) {
    if (pr < 1e-12) return best;
    const g = (rMin * (1 + 1e-9)) / pr;
    px *= g; pz *= g;
    if (!insideRoundRect(px, pz, s, inset * (1 - 1e-6))) return best;
  }
  qTmp[0] = px;
  qTmp[1] = py;
  qTmp[2] = pz;
  const d2 = (qTmp[0] - x) ** 2 + (qTmp[1] - y) ** 2 + (qTmp[2] - z) ** 2;
  if (d2 < best) {
    best = d2;
    out[0] = qTmp[0]; out[1] = qTmp[1]; out[2] = qTmp[2];
  }
  return best;
}

/**
 * Footprint (disc ∪ pocket plan) point farthest from (cx, cz) — the end-gas site opposite the flame —
 * written to out = [x, z]; returns the distance. Disc: the liner point opposite the centre.
 */
export function farthestFootprintPoint(cx: number, cz: number, s: ChamberShape, out: [number, number]): number {
  const n = Math.hypot(cx, cz);
  const ux = n > 1e-6 ? -cx / n : -1, uz = n > 1e-6 ? -cz / n : 0;
  out[0] = ux * s.R;
  out[1] = uz * s.R;
  let best = n + s.R;
  if (s.kind !== 'l-head') return best;
  for (let k = 0; k < 4; k++) {
    const qx = k & 1 ? s.ix1 : s.ix0, qz = k & 2 ? s.iz1 : s.iz0;
    const dx = qx - cx, dz = qz - cz;
    const d = Math.hypot(dx, dz);
    if (d + s.rc > best) {
      best = d + s.rc;
      const k2 = d > 1e-12 ? s.rc / d : 0;
      out[0] = qx + dx * k2;
      out[1] = qz + dz * k2;
    }
  }
  return best;
}
