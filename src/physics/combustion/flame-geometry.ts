/**
 * Exact geometry of a spherical flame intersected with a flat-disc combustion chamber.
 *
 * Chamber at any instant (cylinder frame, DESIGN.md §Frames):
 *   D(h) = { x² + z² ≤ R², −h ≤ y ≤ 0 }  (flat head at y = 0, flat piston crown at y = −h).
 * Flame: ball of radius r centred at the spark gap c = (cx, cy, cz), cy ≤ 0.
 *
 * For ball ∩ D we return (all SI):
 *   volume        m³  enflamed volume |B(c,r) ∩ D|
 *   frontArea     m²  area of the sphere surface lying inside D (the flame front) — equals ∂volume/∂r
 *   wettedHead    m²  area of the head disc (y = 0) inside the ball
 *   wettedPiston  m²  area of the piston disc (y = −h) inside the ball
 *   wettedLiner   m²  area of the liner (x² + z² = R², −h ≤ y ≤ 0) inside the ball
 *
 * ── Exact formulation ────────────────────────────────────────────────────────────────────────────
 * Slice the ball by planes y = cy + t. The slice is a disc of radius ρ(t) = √(r² − t²) whose centre is
 * a horizontal distance d = √(cx² + cz²) from the cylinder axis; its intersection with the bore disc
 * (radius R) is a circle–circle lens. With c₁ = R − d and c₂ = R + d:
 *   L(ρ) = lens area = πρ² (ρ ≤ c₁), πR² (ρ ≥ c₂), ρ² seg(α) + R² seg(β) otherwise, seg(x) = x − sin x cos x
 *   Φ(ρ) = angle of the slice circle inside the bore = 2α   (2π for ρ ≤ c₁, 0 for ρ ≥ c₂)
 *   Ψ(ρ) = angle of the bore circle inside the slice = 2β   (0 for ρ ≤ c₁, 2π for ρ ≥ c₂)
 * where α, β are the half-angles subtended by the common chord at the two centres.
 * Archimedes' hat-box theorem (dA = r dφ dt on a sphere) gives the 1-D integrals over t ∈ [a, b],
 * a = max(−h − cy, −r), b = min(−cy, r):
 *   V = ∫ L(ρ(t)) dt,   A_f = ∫ r Φ(ρ(t)) dt,   A_liner = ∫ R Ψ(ρ(t)) dt,
 * and the planar wetted areas are single lenses: A_head = L(√(r² − cy²)), A_piston = L(√(r² − (h+cy)²)).
 * ∂V/∂r = A_f exactly because dL/dρ = ρΦ(ρ). (Keck 1982, 19th Symp. (Int.) Combust. pp. 1451–1466,
 * eq. 3.9, uses this identity to *define* the spherical burning area A_b = ∂V_b/∂r_b.)
 *
 * ── Fast path ────────────────────────────────────────────────────────────────────────────────────
 * The integrands depend on t only through ρ(t) (even in t), so V, A_f, A_liner are sums of one-side
 * integrals S(r, τ) = ∫₀^min(τ,r) (·) dt with τ = −cy (head side) and τ = |h + cy| (piston side, added
 * or subtracted by which side of the centre the piston plane lies).
 *  • r ≤ c₁ (flame not touching the liner): every slice lies inside the bore — all closed form.
 *  • r > c₁: L = Ref(ρ) − Dev(ρ) with a smooth reference Ref = πR²w + πρ²(1 − w),
 *    w = 1 − ((c₂² − ρ²)/(c₂² − c₁²))² on [c₁, c₂] (Ref = L outside the lens range). Ref is a polynomial in
 *    ρ² = r² − t², so its t-integral and r-derivative are closed form. The deviation integrals
 *    F = ∫ Dev dt and Y = ∫ RΨ dt are stored as bicubic-Hermite tables in normalised coordinates (p, w):
 *    p maps r ∈ [c₁, c₂] (quadratic clustering at both ends) and r ∈ [c₂, r_max] (quadratic clustering at
 *    c₂); w maps τ ∈ [t₂, t₁] (t₁ = √(r² − c₁²), t₂ = √(r² − c₂²)⁺; for τ ≤ t₂ the slices cover the bore
 *    and F = 0) with quadratic clustering at both ends. Every kink of the problem (liner first touch
 *    r = c₁, equator passing the far wall r = c₂, cut plane at the ρ = c₁ and ρ = c₂ circles) lies on a
 *    table edge where the stored functions are C¹–C⁵. Dev vanishes at both ends of the lens range and its
 *    r-derivative content is ∝ that of rΦ near ρ = c₂, so interpolation errors in A_f scale with A_f.
 *  • Node values AND exact analytic derivatives (∂/∂p, ∂/∂w, ∂²/∂p∂w) come from adaptive Gauss–Kronrod
 *    quadrature; frontArea is the exact r-derivative of the interpolated volume, so dV/dr = A_f holds to
 *    rounding in the fast path too (the entrainment model advances V_f with A_f and inverts V_f → r).
 *  • Two thin radius bands where the front has √-type structure finer than any fixed table resolves —
 *    |r − c₂| < 0.01 R and the first two table cells after r = c₁ — are evaluated by (looser-tolerance)
 *    exact quadrature instead. The flame crosses them in well under 1 % of its growth.
 *  • Spark near the liner (c₁ = R − d ≪ R, e.g. the CFR side plug 1 mm from the wall): just after first
 *    contact the stored integrals vary on the scale c₁, so the inner-table size grows as √(2d/c₁)
 *    ({@link defaultInnerCells}); fast/exact front-area error ≤ 1e-3 over the whole (r, h) range and
 *    ≤ 3e-4 in the early flame for the CFR_F1 spec geometry (review, 2026-09).
 *  • Near-axis spark (d ≪ R): all lens formulas are written in P1 = ρ² − c₁² (no ρ² − R² cancellation),
 *    so d down to the central threshold (1e-9 R) is accurate and the table build stays fast.
 *
 * ── Spark outside the bore planform (d > R; opt-in `allowSparkOutsideBore`) ──────────────────────
 * A side-valve ('l-head') chamber has its plug over the valve pocket beside the bore; this class then
 * describes the BORE COLUMN part of that chamber (combustion/chamber.ts adds the pocket). With
 * e₁ = d − R = |c₁| the lens range is ρ ∈ [e₁, c₂]: L = 0, Φ = Ψ = 0 for ρ ≤ e₁ (the slice circle has not
 * reached the bore), the lens formulas above hold unchanged with the SIGNED c₁ = R − d (they only use
 * c₁ and c₁² = e₁²), and nothing is "inside" (no πρ² regime). The reference becomes Ref = πR²·w with the
 * same w (w(e₁) = 0, w(c₂) = 1, w′(c₂) = 0), so Dev = w·A_out − (1 − w)·L again vanishes at both ends of
 * the lens range; the table maps r ∈ [e₁, c₂] (Δ₁ = c₂ − e₁ = 2R) exactly as [c₁, c₂] above, and the
 * first-contact band uses e₁ as its geometric scale ({@link defaultInnerCellsOutside}). For d < R every
 * expression evaluates exactly as before (the CFR path is bit-identical).
 */

const PI = Math.PI;
const TWO_PI = 2 * Math.PI;

/** Result of {@link FlameGeometry.evaluate}. All SI. */
export interface FlameGeometryResult {
  /** Enflamed volume (ball ∩ chamber), m³. */
  volume: number;
  /** Flame-front area: sphere surface inside the chamber, m². Equals ∂volume/∂r. */
  frontArea: number;
  /** Head (fire-deck, y = 0) area inside the ball, m². */
  wettedHead: number;
  /** Piston-crown (y = −h) area inside the ball, m². */
  wettedPiston: number;
  /** Liner area (between −h and 0) inside the ball, m². */
  wettedLiner: number;
}

/** Allocate a zeroed {@link FlameGeometryResult} (once, outside hot loops). */
export const newFlameGeometryResult = (): FlameGeometryResult => ({
  volume: 0,
  frontArea: 0,
  wettedHead: 0,
  wettedPiston: 0,
  wettedLiner: 0,
});

export interface FlameGeometryOptions {
  /**
   * Accept a spark centre OUTSIDE the bore planform (horizontal offset d > R, e.g. a side-valve chamber's
   * plug over the valve pocket): the ball then meets the bore column only for r > d − R (see the module
   * comment). Default false: such a centre throws, as does a centre ON the bore circle (|d − R| < 1e-9 R).
   */
  allowSparkOutsideBore?: boolean;
  /**
   * Reuse the immutable fast tables of an identical FlameGeometry (same bore, centre and table options)
   * instead of building them again — e.g. the identical cylinders of a multi-cylinder engine, each with
   * its own scratch state ({@link FlameGeometry.clone}). Throws if the geometry differs.
   */
  shareTablesWith?: FlameGeometry;
  /**
   * Largest clearance height h (m) for which the fast table must be valid (it covers flame radii up to
   * the farthest chamber corner at this height). Larger h still works via the slow exact path.
   * Default 2.5 × bore. For an engine pass stroke + clearance height at the lowest compression ratio.
   */
  maxHeight?: number;
  /**
   * Table cells along p for r ∈ [R−d, R+d] (even). Default {@link defaultInnerCells}: 64, raised to
   * 16·√(2d/(R−d)) (≤ 512) when the spark is close to the liner, so the cells just after the flame first
   * touches the liner stay small compared with the local geometric scale R − d.
   */
  cellsInner?: number;
  /** Table cells along p for r ∈ [R+d, r_max] (default 64). */
  cellsOuter?: number;
  /** Table cells along the normalised τ coordinate w (even; default 64). */
  cellsTau?: number;
  /**
   * Half-width (m) of the band |r − (R + d)| in which the fast path falls back to exact quadrature
   * (default 0.01 R). The front there is a thin sliver whose area has √-type structure the table
   * cannot resolve to 0.5 %.
   */
  cornerBand?: number;
}

/**
 * Default number of inner-table cells for bore radius R and spark offset d (even, 64…512).
 *
 * Right after the flame first touches the liner (r = c1 = R − d) the stored deviation integrals vary on
 * the length scale c1 in r, while the quadratic end clustering of the inner coordinate gives cells of
 * width Δr ≈ 2.8 √(2d·c1)/n there — i.e. Δr/c1 ≈ 2.8 √(2d/c1)/n. With the historical fixed n = 64 the
 * front-area error just after liner contact grew as (2d/c1)^{3/2}: 2.5e-3 for the CFR spec plug
 * (d = R − 1 mm), 8.6e-3 at d = 0.99R (vs exact quadrature). n = 16√(2d/c1) keeps Δr/c1 ≲ 0.18 and the
 * error ≲ 3e-4 (tested); capped at 512 (d ≥ 0.9976 R then degrades gracefully).
 */
export function defaultInnerCells(R: number, d: number): number {
  const c1 = R - d;
  if (!(c1 > 0)) return 512;
  const n = 2 * Math.ceil(8 * Math.sqrt((2 * d) / c1));
  return n < 64 ? 64 : n > 512 ? 512 : n;
}

/**
 * Default number of inner-table cells for a spark OUTSIDE the bore planform (d > R): the table spans
 * r ∈ [e₁, c₂] (e₁ = d − R, width 2R) and just after first contact the stored integrals vary on the scale
 * e₁ (two circles touching externally), so the same rule as {@link defaultInnerCells} with d → R, c₁ → e₁:
 * n = 16·√(2R/e₁), even, 64…512.
 */
export function defaultInnerCellsOutside(R: number, d: number): number {
  const e1 = d - R;
  if (!(e1 > 0)) return 512;
  const n = 2 * Math.ceil(8 * Math.sqrt((2 * R) / e1));
  return n < 64 ? 64 : n > 512 ? 512 : n;
}

// ---------------------------------------------------------------------------------------------
// Circle–circle lens (bore disc radius R centred on the axis; slice disc radius ρ at offset d)
// ---------------------------------------------------------------------------------------------

/** seg(x) = x − sin x cos x (twice the area of a unit-circle segment of half-angle x), accurate for small x. */
function seg(x: number): number {
  if (x < 0.05) {
    const x2 = x * x;
    // (2x)³/12 − (2x)⁵/240 + (2x)⁷/10080 − (2x)⁹/725760
    return x * x2 * (2 / 3 - x2 * (2 / 15 - x2 * (4 / 315 - x2 * (2 / 2835))));
  }
  return x - Math.sin(x) * Math.cos(x);
}

/**
 * Area (m²) of the intersection of a disc of radius `rho` with a disc of radius `R` whose centres are
 * `d` apart (circle–circle lens). Exact, closed form; any d ≥ 0 (d > R: 0 until ρ > d − R).
 */
export function lensArea(rho: number, R: number, d: number): number {
  if (rho <= 0) return 0;
  const c1 = R - d;
  const c2 = R + d;
  if (rho <= c1) return PI * rho * rho;
  if (rho >= c2) return PI * R * R;
  if (rho <= -c1) return 0; // d > R: the discs do not meet yet (never true for d < R)
  // common-chord half-length ℓ = √((ρ² − c1²)(c2² − ρ²)) / (2d); half-angles α (at the slice centre),
  // β (at the axis): ρ cos α = (ρ² − c1c2)/(2d), R cos β = (d² + R² − ρ²)/(2d); L = ρ²α + R²β − dℓ.
  // With P1 = ρ² − c1² these are P1/(2d) − c1 and R − P1/(2d) (no R² − ρ² cancellation for small d).
  const P1 = (rho - c1) * (rho + c1);
  const q = P1 / (2 * d);
  const l = Math.sqrt(P1 * (c2 - rho) * (c2 + rho)) / (2 * d);
  const al = Math.atan2(l, q - c1);
  const be = Math.atan2(l, R - q);
  return rho * rho * al + R * R * be - d * l;
}

/** Angle (rad, 0..2π) of the circle of radius `rho` (centre offset d from the axis) lying inside the bore disc R. */
export function arcInsideBore(rho: number, R: number, d: number): number {
  const c1 = R - d;
  const c2 = R + d;
  if (rho <= c1) return TWO_PI;
  if (rho >= c2) return 0;
  if (rho <= -c1) return 0; // d > R: the circle lies outside the bore disc
  const P1 = (rho - c1) * (rho + c1);
  const l = Math.sqrt(P1 * (c2 - rho) * (c2 + rho)) / (2 * d);
  return 2 * Math.atan2(l, P1 / (2 * d) - c1);
}

/** Angle (rad, 0..2π) of the bore circle R lying inside the disc of radius `rho` centred at offset d. */
export function arcOfBoreInside(rho: number, R: number, d: number): number {
  const c1 = R - d;
  const c2 = R + d;
  if (rho <= c1) return 0;
  if (rho >= c2) return TWO_PI;
  if (rho <= -c1) return 0; // d > R: the disc has not reached the bore circle
  const P1 = (rho - c1) * (rho + c1);
  const l = Math.sqrt(P1 * (c2 - rho) * (c2 + rho)) / (2 * d);
  return 2 * Math.atan2(l, R - P1 / (2 * d));
}

// ---------------------------------------------------------------------------------------------
// Adaptive Gauss–Kronrod (7/15) for small vector integrands. Build-time / exact path only.
// Abscissae and weights: Piessens et al. 1983, QUADPACK, routine QK15.
// ---------------------------------------------------------------------------------------------

const XGK = [
  0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
  0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
  0.586087235467691130294144845693013, 0.405845151377397166906606412076961,
  0.207784955007898467600689403773245, 0.0,
];
const WGK = [
  0.02293532201052922496373200805897, 0.063092092629978553290700663189204,
  0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
  0.16900472663926790282658342659855, 0.190350578064785409913256402421014,
  0.204432940075298892414161999234649, 0.209482141084727828012999174891714,
];
/** Gauss 7-point weights for XGK[1], XGK[3], XGK[5], XGK[7]. */
const WG = [
  0.129484966168869693270611432679082, 0.27970539148927666790146777142378,
  0.381830050505118944950369775488975, 0.417959183673469387755102040816327,
];

/**
 * Vector integrand at abscissa t. `dLo` = t − lo and `dHi` = hi − t are supplied to full relative
 * precision (so integrands with square-root endpoint behaviour can be evaluated without cancellation).
 */
type VecIntegrand = (t: number, dLo: number, dHi: number, out: Float64Array) => void;

const NC_MAX = 4;

/** 4-point Gauss–Legendre on [−1, 1] (exact for polynomials of degree ≤ 7). */
const GL4X = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
const GL4W = [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538];

/**
 * ∫_lo^hi f dt for a vector integrand that may behave like √ or 1/√ at BOTH ends: the substitution
 * t = lo + (hi − lo)(1 − cos πs)/2 makes such integrands smooth in s; then adaptive GK15 bisection.
 * Accumulates into acc[0..nc).
 */
function integrateCosSub(
  f: VecIntegrand,
  nc: number,
  lo: number,
  hi: number,
  acc: Float64Array,
  absTol: Float64Array,
  relTol: number,
): void {
  if (!(hi > lo)) return;
  const half = 0.5 * (hi - lo);
  const fv = new Float64Array(NC_MAX);
  const resK = new Float64Array(NC_MAX);
  const resG = new Float64Array(NC_MAX);
  const stack: number[] = [0, 1, 0]; // [s0, s1, depth]
  let evals = 0;
  while (stack.length > 0) {
    const depth = stack.pop()!;
    const s1 = stack.pop()!;
    const s0 = stack.pop()!;
    const c = 0.5 * (s0 + s1);
    const hl = 0.5 * (s1 - s0);
    resK.fill(0);
    resG.fill(0);
    for (let i = 0; i < 15; i++) {
      const k = i < 8 ? i : 14 - i;
      const s = c + hl * (i < 8 ? -XGK[k] : XGK[k]);
      const sh = Math.sin(0.5 * PI * s);
      const ch = Math.cos(0.5 * PI * s);
      const dLo = 2 * half * sh * sh;
      const dHi = 2 * half * ch * ch;
      const jac = half * PI * 2 * sh * ch;
      f(lo + dLo, dLo, dHi, fv);
      evals++;
      const wk = WGK[k];
      const isGauss = (k & 1) === 1;
      const wg = isGauss ? WG[(k - 1) >> 1] : 0;
      for (let j = 0; j < nc; j++) {
        const v = fv[j] * jac;
        resK[j] += wk * v;
        if (isGauss) resG[j] += wg * v;
      }
    }
    let ok = true;
    for (let j = 0; j < nc; j++) {
      resK[j] *= hl;
      resG[j] *= hl;
      if (Math.abs(resK[j] - resG[j]) > Math.max(absTol[j] * hl, relTol * Math.abs(resK[j]))) ok = false;
    }
    if (ok || depth >= 40 || evals > 400000) {
      for (let j = 0; j < nc; j++) acc[j] += resK[j];
    } else {
      stack.push(s0, c, depth + 1, c, s1, depth + 1);
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Coordinate maps (piecewise quadratic, C¹, cheap inverse)
// ---------------------------------------------------------------------------------------------

/** Clustering at both ends: σ(p) = 2p² (p ≤ ½), 1 − 2(1 − p)² (p > ½). */
const mapBoth = (p: number): number => (p <= 0.5 ? 2 * p * p : 1 - 2 * (1 - p) * (1 - p));
const mapBothD = (p: number): number => (p <= 0.5 ? 4 * p : 4 * (1 - p));
const mapBothInv = (q: number): number => (q <= 0.5 ? Math.sqrt(0.5 * q) : 1 - Math.sqrt(0.5 * (1 - q)));
/** 1 − σ(p), accurate near p = 1. */
const mapBothC = (p: number): number => (p <= 0.5 ? 1 - 2 * p * p : 2 * (1 - p) * (1 - p));

/** Nudge of node data off w = 0, 1 (where ∂τ/∂w = 0 makes some derivative data 0·∞ limits). */
const EPS_NODE = 1e-7;
/** Lookups are clamped this far from the map singular points (σ′ = 0). */
const EPS_LOOK = 1e-6;
/** Relative tolerance of the exact quadrature used inside the fallback bands. */
const BAND_TOL = 1e-10;
/** Default half-width of the exact-quadrature band around r = R + d, relative to R. */
const CORNER_BAND = 0.01;
/** Below this horizontal offset (relative to R) the spark is treated as central (d = 0, closed form). */
const D_CENTRAL = 1e-9;
/** A spark centre within this distance (relative to R) of the bore circle is rejected (degenerate lens range). */
const D_ON_BORE = 1e-9;

/** Maximum bisection depth of the allocation-free band quadrature. */
const STACK_DEPTH = 48;

/** Node layout: F̂, F̂_p·hp, F̂_w·hw, F̂_pw·hp·hw, Ŷ, Ŷ_p·hp, Ŷ_w·hw, Ŷ_pw·hp·hw. */
const NODE = 8;

/**
 * Ball (centre at the spark gap, radius r) ∩ flat-disc chamber of height h: volume, flame-front area and
 * wetted wall areas — exact ({@link evaluateExact}, quadrature) and fast ({@link evaluate}, table).
 */
export class FlameGeometry {
  /** Cylinder bore, m. */
  readonly bore: number;
  /** Bore radius R, m. */
  readonly radius: number;
  /** Flame centre (spark gap) in the cylinder frame, m. */
  readonly center: readonly [number, number, number];
  /** Horizontal distance d of the flame centre from the cylinder axis, m. */
  readonly offset: number;
  /** Distance from the flame centre up to the head plane, m (= −center[1] ≥ 0). */
  readonly headDistance: number;
  /** Clearance height up to which the fast table is valid, m. */
  readonly maxHeight: number;
  /** Largest flame radius covered by the fast table, m. */
  readonly maxTableRadius: number;

  /** True when the spark centre lies outside the bore planform (d > R; opt-in). */
  readonly outside: boolean;

  /** Signed R − d (the lens algebra uses c1 and c1² = (d − R)² for either sign). */
  private readonly c1: number;
  private readonly c2: number;
  /** Lower end of the lens range |R − d| (= c1 bit for bit when d < R; e₁ = d − R when outside). */
  private readonly lo: number;
  private readonly central: boolean;
  private readonly options: FlameGeometryOptions;
  /** c2² − c1² = 4Rd (blend normalisation), m². */
  private readonly dc2: number;
  /** Half-width of the band around r = R + d evaluated by exact quadrature, m. */
  private readonly cornerBand: number;
  /** Radii in (R − d, touchBand) (first table cells after the flame touches the liner) use exact quadrature, m. */
  private readonly touchBand: number;
  private readonly n1: number;
  private readonly n2: number;
  private readonly nw: number;
  private readonly delta1: number;
  private readonly delta2: number;
  private readonly tab1: Float64Array;
  private readonly tab2: Float64Array;

  // lookup scratch (no allocation in hot paths)
  private cV = 0;
  private cA = 0;
  private cW = 0;
  private sV = 0;
  private sA = 0;
  private sW = 0;
  private sF = 0;
  private sFr = 0;
  private sY = 0;
  private readonly ang = new Float64Array(4);
  private lastRadius = 0;
  private readonly qAcc = new Float64Array(NC_MAX);
  /** Depth-first bisection stack of {@link lensSide}: ≤ 2 open panels per level × 3 numbers. */
  private readonly stk = new Float64Array(6 * STACK_DEPTH + 3);

  /**
   * @param bore cylinder bore, m
   * @param sparkCenter flame (spark-gap) centre in the cylinder frame, m; requires y ≤ 0 and a horizontal
   *        offset d < bore/2 (accuracy of the fast table degrades only as d → bore/2), or d > bore/2 with
   *        `allowSparkOutsideBore`
   * @param options table range/resolution
   */
  constructor(bore: number, sparkCenter: readonly [number, number, number], options: FlameGeometryOptions = {}) {
    if (!(bore > 0)) throw new Error('FlameGeometry: bore must be > 0');
    const R = 0.5 * bore;
    const [cx, cy, cz] = sparkCenter;
    if (cy > 0) throw new Error('FlameGeometry: spark centre must satisfy y ≤ 0 (below the head plane)');
    let d = Math.hypot(cx, cz);
    let outside = false;
    if (!(d < R)) {
      if (!options.allowSparkOutsideBore) throw new Error('FlameGeometry: spark centre must lie inside the bore (d < R)');
      if (!(d > R * (1 + D_ON_BORE))) throw new Error('FlameGeometry: spark centre on the bore circle (|d − R| < 1e-9 R)');
      outside = true;
    }
    this.outside = outside;
    this.options = options;
    this.central = d < D_CENTRAL * R;
    if (this.central) d = 0;
    this.bore = bore;
    this.radius = R;
    this.center = [cx, cy, cz];
    this.offset = d;
    this.headDistance = -cy;
    this.c1 = R - d;
    this.c2 = R + d;
    this.lo = outside ? d - R : R - d;
    this.maxHeight = options.maxHeight ?? 2.5 * bore;
    const b = this.headDistance;
    const vmax = Math.max(b, Math.abs(this.maxHeight - b));
    this.maxTableRadius = Math.sqrt(this.c2 * this.c2 + vmax * vmax) * (1 + 1e-9);
    this.n1 = options.cellsInner ?? (outside ? defaultInnerCellsOutside(R, d) : defaultInnerCells(R, d));
    this.n2 = options.cellsOuter ?? 64;
    this.nw = options.cellsTau ?? 64;
    if (this.n1 % 2 !== 0 || this.nw % 2 !== 0 || this.n1 < 2 || this.n2 < 1 || this.nw < 2) {
      throw new Error('FlameGeometry: cellsInner and cellsTau must be even and ≥ 2, cellsOuter ≥ 1');
    }
    this.dc2 = 4 * R * d; // = c2² − c1² exactly (no cancellation for small d)
    this.cornerBand = options.cornerBand ?? CORNER_BAND * R;
    this.delta1 = outside ? 2 * R : 2 * d; // = c2 − lo exactly
    // first two p-cells of the inner table: r − lo < 2·D1·(2/n1)²
    this.touchBand = this.lo + 2 * this.delta1 * (2 / this.n1) * (2 / this.n1);
    this.delta2 = Math.max(this.maxTableRadius - this.c2, 1e-9 * R);
    const src = options.shareTablesWith;
    if (src) {
      if (
        src.bore !== bore ||
        src.center[0] !== cx ||
        src.center[1] !== cy ||
        src.center[2] !== cz ||
        src.maxHeight !== this.maxHeight ||
        src.n1 !== this.n1 ||
        src.n2 !== this.n2 ||
        src.nw !== this.nw ||
        src.cornerBand !== this.cornerBand
      ) {
        throw new Error('FlameGeometry: shareTablesWith needs an identical geometry and table layout');
      }
      this.tab1 = src.tab1;
      this.tab2 = src.tab2;
    } else if (this.central) {
      this.tab1 = new Float64Array(0);
      this.tab2 = new Float64Array(0);
    } else {
      this.tab1 = new Float64Array((this.n1 + 1) * (this.nw + 1) * NODE);
      this.tab2 = new Float64Array((this.n2 + 1) * (this.nw + 1) * NODE);
      this.buildTable(1);
      this.buildTable(2);
    }
  }

  /**
   * A FlameGeometry with the same geometry that SHARES this instance's immutable tables (no rebuild) and has
   * its own scratch / warm-start state — one per cylinder of a multi-cylinder engine.
   */
  clone(): FlameGeometry {
    return new FlameGeometry(this.bore, this.center, { ...this.options, shareTablesWith: this });
  }

  /** Distance from the flame centre to the farthest chamber corner at clearance height h, m. */
  maxRadius(h: number): number {
    const b = this.headDistance;
    const v = Math.max(b, Math.abs(b - h));
    return Math.sqrt(this.c2 * this.c2 + v * v);
  }

  /**
   * Fast evaluation (table + closed forms; allocation-free). ~0.2 µs typical; ~5–15 µs for radii in the
   * two narrow exact-quadrature bands (see the module comment). Writes into `out`, returns it.
   * @param r flame radius, m (≥ 0)
   * @param h instantaneous clearance height, m (> 0)
   */
  evaluate(r: number, h: number, out: FlameGeometryResult): FlameGeometryResult {
    const R = this.radius;
    const b = this.headDistance;
    const ta = Math.abs(b - h);
    this.core(r, h);
    out.volume = this.cV;
    out.frontArea = this.cA;
    out.wettedLiner = this.cW;
    if (!(r > 0) || !(h > 0)) {
      out.wettedHead = 0;
      out.wettedPiston = 0;
    } else if (this.cA === 0 && this.cV > 0) {
      out.wettedHead = PI * R * R;
      out.wettedPiston = PI * R * R;
    } else {
      out.wettedHead = r > b ? lensArea(Math.sqrt((r - b) * (r + b)), R, this.offset) : 0;
      out.wettedPiston = r > ta ? lensArea(Math.sqrt((r - ta) * (r + ta)), R, this.offset) : 0;
    }
    return out;
  }

  /**
   * Fast volume, front area and liner area only (no planar lenses; wettedHead / wettedPiston of `out` are
   * left untouched). Same numbers as {@link evaluate}. Allocation-free; returns `out`.
   */
  evaluateCore(r: number, h: number, out: FlameGeometryResult): FlameGeometryResult {
    this.core(r, h);
    out.volume = this.cV;
    out.frontArea = this.cA;
    out.wettedLiner = this.cW;
    return out;
  }

  /** Fast volume, front area and liner area into cV, cA, cW (no planar wetted areas). */
  private core(r: number, h: number): void {
    const R = this.radius;
    const b = this.headDistance;
    // Piston plane at signed axial offset a = b − h from the centre (a < 0: below the centre).
    const a = b - h;
    if (!(r > 0) || !(h > 0)) {
      // no flame, or a degenerate chamber (h ≤ 0 would otherwise give negative volumes)
      this.cV = 0;
      this.cA = 0;
      this.cW = 0;
      return;
    }
    const ta = a < 0 ? -a : a;
    const vmax = b > ta ? b : ta;
    if (r * r >= this.c2 * this.c2 + vmax * vmax) {
      // beyond the farthest corner: the ball contains the whole chamber
      this.cV = PI * R * R * h;
      this.cA = 0;
      this.cW = TWO_PI * R * h;
      return;
    }
    // head side: t ∈ [0, b]; piston side: add [a, 0] if a < 0, else subtract [0, a]
    this.side(r, b);
    const V = this.sV;
    const Af = this.sA;
    const Wl = this.sW;
    const sgn = a < 0 ? 1 : -1;
    this.side(r, ta);
    this.cV = V + sgn * this.sV;
    this.cA = Af + sgn * this.sA;
    this.cW = Wl + sgn * this.sW;
  }

  /** Enflamed volume only (fast path), m³. */
  volume(r: number, h: number): number {
    this.core(r, h);
    return this.cV;
  }

  /** Flame-front area only (fast path), m². */
  frontArea(r: number, h: number): number {
    this.core(r, h);
    return this.cA;
  }

  /**
   * Inverse of the fast volume: the flame radius r (m) with volume(r, h) = V. Safeguarded Newton using
   * dV/dr = frontArea (consistent with the table), warm-started from the previous call or `rGuess`.
   * Returns 0 for V ≤ 0 and {@link maxRadius}(h) for V ≥ πR²h.
   * @param V enflamed volume, m³
   * @param h clearance height, m
   * @param rGuess optional initial guess, m
   */
  radiusForVolume(V: number, h: number, rGuess?: number): number {
    if (!(V > 0) || !(h > 0)) return 0;
    const R = this.radius;
    const Vch = PI * R * R * h;
    const rMax = this.maxRadius(h);
    if (V >= Vch) return rMax;
    // closed form while the ball is strictly inside the chamber (never for a spark outside the bore,
    // whose ball first meets the column at r = d − R)
    const b = this.headDistance;
    const rIn = this.outside ? 0 : h > b ? Math.min(b, h - b, this.c1) : 0;
    const rs = Math.cbrt((3 * V) / (4 * PI));
    if (rs <= rIn) return rs;
    let lo = this.outside ? this.lo : rIn;
    let hi = rMax;
    let r = rGuess !== undefined && rGuess > lo && rGuess < hi ? rGuess : this.lastRadius;
    if (!(r > lo && r < hi)) r = Math.min(Math.max(rs, lo), 0.5 * (lo + hi));
    for (let it = 0; it < 200; it++) {
      this.core(r, h);
      const f = this.cV - V;
      if (f > 0) hi = r;
      else lo = r;
      const tolR = 2e-15 * r;
      if (f === 0 || hi - lo <= tolR) break;
      let rn = this.cA > 0 ? r - f / this.cA : NaN;
      if (!(rn > lo && rn < hi)) rn = 0.5 * (lo + hi);
      if (Math.abs(rn - r) <= tolR) {
        r = rn;
        break;
      }
      r = rn;
    }
    this.lastRadius = r;
    return r;
  }

  /**
   * Exact evaluation by adaptive Gauss–Kronrod quadrature of the slice integrals (no table; relative
   * accuracy ~1e-13; ~10–50 µs). Reference path for validation and for h beyond `maxHeight`.
   */
  evaluateExact(r: number, h: number, out: FlameGeometryResult): FlameGeometryResult {
    const R = this.radius;
    const d = this.offset;
    const b = this.headDistance;
    const a = b - h;
    const ta = a < 0 ? -a : a;
    out.wettedHead = r > b ? lensArea(Math.sqrt((r - b) * (r + b)), R, d) : 0;
    out.wettedPiston = r > ta ? lensArea(Math.sqrt((r - ta) * (r + ta)), R, d) : 0;
    if (!(r > 0) || !(h > 0)) {
      out.wettedHead = 0;
      out.wettedPiston = 0;
      out.volume = 0;
      out.frontArea = 0;
      out.wettedLiner = 0;
      return out;
    }
    this.sideExact(r, b);
    let V = this.sV;
    let Af = this.sA;
    let Wl = this.sW;
    const sgn = a < 0 ? 1 : -1;
    this.sideExact(r, ta);
    out.volume = V + sgn * this.sV;
    out.frontArea = Af + sgn * this.sA;
    out.wettedLiner = Wl + sgn * this.sW;
    return out;
  }

  // -------------------------------------------------------------------------------------------
  // internals
  //
  // Per side of the flame centre we need, for 0 ≤ τ (distance from the centre plane to the cut plane):
  //   S_V(r, τ) = ∫₀^m L(ρ(t)) dt,  S_A = ∂S_V/∂r|τ = ∫₀^m rΦ dt,  S_W = ∫₀^m RΨ dt,  m = min(τ, r).
  // For r > c1 we split L = Ref(ρ) − Dev(ρ) with a smooth reference whose t-integral is a polynomial:
  //   Ref(ρ) = πR² w + πρ² (1 − w),   w(ρ) = 1 − ((c2² − ρ²)/(c2² − c1²))²   on c1 ≤ ρ ≤ c2,
  //   Ref = πR² (ρ ≥ c2), πρ² (ρ ≤ c1)  (continuous; = L outside the lens range),
  //   Dev = w·A_out + (1 − w)·E,  A_out = πR² − L (bore area outside the slice), E = πρ² − L.
  // Dev vanishes at both ends of the lens range and its r-derivative content ∝ that of rΦ near ρ = c2,
  // so the tabulated F = ∫ Dev dt carries interpolation errors that scale with the front area itself.
  // -------------------------------------------------------------------------------------------

  /**
   * Lens-regime slice angles from P1 = ρ² − c1² ≥ 0 and P2 = c2² − ρ² ≥ 0 (given to full relative
   * precision). Writes o = [α, π − α, β, π − β] (half-angles of the common chord seen from the slice
   * centre and from the axis) and returns the chord half-length ℓ.
   */
  private lensAngles(P1: number, P2: number, o: Float64Array): number {
    const R = this.radius;
    const d = this.offset;
    const c1 = this.c1;
    const l = Math.sqrt(P1 * P2) / (2 * d);
    // ρ cos α = (ρ² − c1c2)/(2d) = P1/(2d) − c1,  R cos β = (d² + R² − ρ²)/(2d) = R − P1/(2d):
    // written without the ρ² − R² cancellation, which loses ~log10(R/d) digits for a near-axis spark
    // (and stalled the adaptive table build for d ≲ 1e-4 R).
    const q = P1 / (2 * d);
    const xa = q - c1; // ρ cos α
    const xb = R - q; // R cos β
    o[0] = Math.atan2(l, xa);
    o[1] = Math.atan2(l, -xa);
    o[2] = Math.atan2(l, xb);
    o[3] = Math.atan2(l, -xb);
    return l;
  }

  /**
   * Exact one-side integrals (adaptive quadrature): sV = ∫₀^m L dt, sA = ∫₀^m rΦ dt, sW = ∫₀^m RΨ dt,
   * m = min(τ, r).
   */
  private sideExact(r: number, tau: number, relTol = 1e-14): void {
    this.sV = 0;
    this.sA = 0;
    this.sW = 0;
    if (!(tau > 0) || !(r > 0)) return;
    if (this.outside && r <= this.lo) return; // the ball has not reached the bore column
    const R = this.radius;
    const c1 = this.c1;
    const c2 = this.c2;
    const m = tau < r ? tau : r;
    const t1 = r > c1 ? Math.sqrt((r - c1) * (r + c1)) : 0;
    const t2 = r > c2 ? Math.sqrt((r - c2) * (r + c2)) : 0;
    // ρ ≥ c2 (t ≤ t2): slice covers the bore — L = πR², Φ = 0, Ψ = 2π
    const aEnd = m < t2 ? m : t2;
    this.sV += PI * R * R * aEnd;
    this.sW += TWO_PI * R * aEnd;
    // ρ ≤ c1 (t ≥ t1): slice inside the bore — L = πρ², Φ = 2π, Ψ = 0 (outside spark: ρ ≤ e1, L = 0)
    if (m > t1 && !this.outside) {
      this.sV += PI * (r * r * (m - t1) - (m * m * m - t1 * t1 * t1) / 3);
      this.sA += TWO_PI * r * (m - t1);
    }
    // lens range t ∈ [t2, min(m, t1)]
    const hiT = m < t1 ? m : t1;
    if (this.central || !(hiT > t2)) return;
    this.lensSide(r, t1, t2, hiT, relTol);
  }

  /**
   * Adds ∫_{t2}^{hi} [L, rΦ, RΨ] dt (lens range) to sV, sA, sW: cos substitution + adaptive GK15 with
   * the integrand inlined and a preallocated stack — allocation-free (used in the fast path's bands).
   */
  private lensSide(r: number, t1: number, t2: number, hi: number, relTol: number): void {
    const R = this.radius;
    const c1 = this.c1;
    const c2 = this.c2;
    const ang = this.ang;
    const half = 0.5 * (hi - t2);
    const t1mHi = t1 - hi;
    const c2r = (c2 - r) * (c2 + r);
    const above = r >= c2;
    const tolV = 0.1 * relTol * PI * R * R * t1;
    const tolA = 0.1 * relTol * TWO_PI * r * t1;
    const tolW = 0.1 * relTol * TWO_PI * R * t1;
    const stk = this.stk;
    let sp = 0;
    stk[sp++] = 0;
    stk[sp++] = 1;
    stk[sp++] = 0;
    let evals = 0;
    while (sp > 0) {
      const depth = stk[--sp];
      const s1 = stk[--sp];
      const s0 = stk[--sp];
      const c = 0.5 * (s0 + s1);
      const hl = 0.5 * (s1 - s0);
      let kV = 0;
      let kA = 0;
      let kW = 0;
      let gV = 0;
      let gA = 0;
      let gW = 0;
      for (let i = 0; i < 15; i++) {
        const k = i < 8 ? i : 14 - i;
        const sv = c + hl * (i < 8 ? -XGK[k] : XGK[k]);
        const sh = Math.sin(0.5 * PI * sv);
        const ch = Math.cos(0.5 * PI * sv);
        const dLo = 2 * half * sh * sh;
        const dHi = 2 * half * ch * ch;
        const jac = half * PI * 2 * sh * ch;
        const t = t2 + dLo;
        const P1 = (t1mHi + dHi) * (t1 + t); // ρ² − c1²
        const P2 = above ? dLo * (t + t2) : t * t + c2r; // c2² − ρ²
        this.lensAngles(P1, P2, ang);
        const rho2 = c1 * c1 + P1;
        const fV = (rho2 * seg(ang[0]) + R * R * seg(ang[2])) * jac; // L
        const fA = r * 2 * ang[0] * jac; // rΦ
        const fW = R * 2 * ang[2] * jac; // RΨ
        const wk = WGK[k];
        kV += wk * fV;
        kA += wk * fA;
        kW += wk * fW;
        if ((k & 1) === 1) {
          const wg = WG[(k - 1) >> 1];
          gV += wg * fV;
          gA += wg * fA;
          gW += wg * fW;
        }
      }
      evals += 15;
      kV *= hl;
      kA *= hl;
      kW *= hl;
      const ok =
        Math.abs(kV - gV * hl) <= Math.max(tolV * hl, relTol * Math.abs(kV)) &&
        Math.abs(kA - gA * hl) <= Math.max(tolA * hl, relTol * Math.abs(kA)) &&
        Math.abs(kW - gW * hl) <= Math.max(tolW * hl, relTol * Math.abs(kW));
      if (ok || depth >= STACK_DEPTH - 1 || evals > 200000) {
        this.sV += kV;
        this.sA += kA;
        this.sW += kW;
      } else {
        stk[sp++] = s0;
        stk[sp++] = c;
        stk[sp++] = depth + 1;
        stk[sp++] = c;
        stk[sp++] = s1;
        stk[sp++] = depth + 1;
      }
    }
  }

  /** Sets sV, sA, sW for one side (see the block comment above). */
  private side(r: number, tau: number): void {
    if (!(tau > 0)) {
      this.sV = 0;
      this.sA = 0;
      this.sW = 0;
      return;
    }
    const R = this.radius;
    const c1 = this.c1;
    const m = tau < r ? tau : r;
    if (r <= c1) {
      // every slice lies inside the bore: L = πρ², Φ = 2π, Ψ = 0
      this.sV = PI * (r * r * m - (m * m * m) / 3);
      this.sA = TWO_PI * r * m;
      this.sW = 0;
      return;
    }
    if (this.outside && r <= this.lo) {
      // spark outside the bore: no slice has reached the bore column yet
      this.sV = 0;
      this.sA = 0;
      this.sW = 0;
      return;
    }
    if (this.central) {
      // d = 0: L = π min(ρ, R)², closed form
      const t1c = Math.sqrt((r - R) * (r + R));
      const mc = tau < t1c ? tau : t1c;
      this.sV = PI * (r * r * m - (m * m * m) / 3) - PI * ((r * r - R * R) * mc - (mc * mc * mc) / 3);
      this.sA = TWO_PI * r * (m - mc);
      this.sW = TWO_PI * R * mc;
      return;
    }
    const c2 = this.c2;
    if ((r > c2 - this.cornerBand && r < c2 + this.cornerBand) || r < this.touchBand) {
      // |r − (R + d)| small: the table cannot resolve the √-type structure of the front there
      this.sideExact(r, tau, BAND_TOL);
      return;
    }
    // ---- reference part (closed form) ----
    const t1 = Math.sqrt((r - c1) * (r + c1));
    const t2 = r > c2 ? Math.sqrt((r - c2) * (r + c2)) : 0;
    let V = 0;
    let A = 0;
    const la = m < t2 ? m : t2;
    V += PI * R * R * la; // ρ ≥ c2: Ref = πR² (independent of r)
    const tb = m < t1 ? m : t1;
    if (tb > t2 && !this.outside) {
      // Ref = πR² − (π/Δ²)(R² − ρ²)(c2² − ρ²)², a degree-6 polynomial in t: 4-point Gauss is exact
      const R2mr2 = R * R - r * r;
      const c2mr2 = (c2 - r) * (c2 + r);
      const hc = 0.5 * (tb - t2);
      const mid = 0.5 * (tb + t2);
      let sv = 0;
      let sa = 0;
      for (let k = 0; k < 4; k++) {
        const t = mid + hc * GL4X[k];
        const P2 = r > c2 ? (t - t2) * (t + t2) : t * t + c2mr2; // c2² − ρ²
        const Q = R2mr2 + t * t; // R² − ρ²
        sv += GL4W[k] * Q * P2 * P2;
        sa += GL4W[k] * P2 * (P2 + 2 * Q);
      }
      const inv = PI / (this.dc2 * this.dc2);
      V += PI * R * R * (tb - t2) - inv * hc * sv;
      A += 2 * r * inv * hc * sa; // ∂/∂r of −(π/Δ²)(R² − ρ²)P2² = (2πr/Δ²) P2 (P2 + 2(R² − ρ²))
    } else if (tb > t2) {
      // outside spark: Ref = πR² w = πR² − (πR²/Δ²)(c2² − ρ²)², degree 4 in t: 4-point Gauss is exact
      const c2mr2 = (c2 - r) * (c2 + r);
      const hc = 0.5 * (tb - t2);
      const mid = 0.5 * (tb + t2);
      let sv = 0;
      let sa = 0;
      for (let k = 0; k < 4; k++) {
        const t = mid + hc * GL4X[k];
        const P2 = r > c2 ? (t - t2) * (t + t2) : t * t + c2mr2; // c2² − ρ²
        sv += GL4W[k] * P2 * P2;
        sa += GL4W[k] * P2;
      }
      const inv = (PI * R * R) / (this.dc2 * this.dc2);
      V += PI * R * R * (tb - t2) - inv * hc * sv;
      A += 4 * r * inv * hc * sa; // ∂/∂r of −(πR²/Δ²)P2² = (4πR² r/Δ²) P2
    }
    if (m > t1 && !this.outside) {
      // ρ ≤ c1: Ref = πρ² = π(r² − t²)
      V += PI * (r * r * (m - t1) - (m * m * m - t1 * t1 * t1) / 3);
      A += TWO_PI * r * (m - t1);
    }
    // ---- tabulated deviation ----
    this.lookup(r, tau);
    this.sV = V - this.sF;
    this.sA = A - this.sFr;
    this.sW = this.sY;
  }

  /**
   * Sets sF = F(r, τ) = ∫₀^min(τ,t₁) Dev dt, sFr = ∂F/∂r|τ, sY = Y(r, τ) = ∫₀^min(τ,t₁) RΨ dt
   * (fast path; exact fallback beyond the table). Requires r > c1, τ > 0, non-central.
   */
  private lookup(r: number, tau: number): void {
    const c1 = this.c1;
    const c2 = this.c2;
    // ---- table coordinate p with the consistent r_eff, dr/dp and t1(r_eff) ----
    let tab: Float64Array;
    let n: number;
    let p: number;
    let rp: number;
    let re: number;
    let t1: number;
    let t2 = 0;
    if (r < c2) {
      // (lo = c1 for a spark inside the bore, e1 = d − R outside)
      const lo = this.lo;
      tab = this.tab1;
      n = this.n1;
      p = mapBothInv((r - lo) / this.delta1);
      if (p < EPS_LOOK) p = EPS_LOOK;
      else if (p > 1 - EPS_LOOK) p = 1 - EPS_LOOK;
      const dr = this.delta1 * mapBoth(p);
      re = lo + dr;
      rp = this.delta1 * mapBothD(p);
      t1 = Math.sqrt(dr * (re + lo));
    } else {
      tab = this.tab2;
      n = this.n2;
      p = Math.sqrt((r - c2) / this.delta2);
      if (p > 1) {
        this.lookupExact(r, tau);
        return;
      }
      if (p < EPS_LOOK) p = EPS_LOOK;
      const dr = this.delta2 * p * p;
      re = c2 + dr;
      rp = 2 * this.delta2 * p;
      t1 = Math.sqrt((re - c1) * (re + c1));
      t2 = Math.sqrt(dr * (re + c2));
      if (tau <= t2) {
        // every slice up to τ covers the bore: Dev = 0, Ψ = 2π
        this.sF = 0;
        this.sFr = 0;
        this.sY = TWO_PI * this.radius * tau;
        return;
      }
    }
    const nw = this.nw;
    const saturated = tau >= t1;
    let w = 1;
    if (!saturated) {
      w = mapBothInv((tau - t2) / (t1 - t2));
      if (w < EPS_LOOK) w = EPS_LOOK;
      else if (w > 1 - EPS_LOOK) w = 1 - EPS_LOOK;
    }
    // ---- cell and local coordinates ----
    let i = Math.floor(p * n);
    if (i >= n) i = n - 1;
    let j = Math.floor(w * nw);
    if (j >= nw) j = nw - 1;
    const u = p * n - i;
    const v = w * nw - j;
    const stride = nw + 1;
    const i00 = (i * stride + j) * NODE;
    const i01 = i00 + NODE;
    const i10 = i00 + stride * NODE;
    const i11 = i10 + NODE;
    // cubic Hermite basis H0 = 2u³−3u²+1, H1 = u³−2u²+u, H2 = −2u³+3u², H3 = u³−u², and u-derivatives
    const u2 = u * u;
    const u3 = u2 * u;
    const v2 = v * v;
    const v3 = v2 * v;
    const hu0 = 2 * u3 - 3 * u2 + 1;
    const hu1 = u3 - 2 * u2 + u;
    const hu2 = -2 * u3 + 3 * u2;
    const hu3 = u3 - u2;
    const du0 = 6 * u2 - 6 * u;
    const du1 = 3 * u2 - 4 * u + 1;
    const du2 = -du0;
    const du3 = 3 * u2 - 2 * u;
    const hv0 = 2 * v3 - 3 * v2 + 1;
    const hv1 = v3 - 2 * v2 + v;
    const hv2 = -2 * v3 + 3 * v2;
    const hv3 = v3 - v2;
    // F: edge functions along v at u = 0 (nodes 00, 01) and u = 1 (nodes 10, 11)
    const g0 = tab[i00] * hv0 + tab[i00 + 2] * hv1 + tab[i01] * hv2 + tab[i01 + 2] * hv3;
    const gp0 = tab[i00 + 1] * hv0 + tab[i00 + 3] * hv1 + tab[i01 + 1] * hv2 + tab[i01 + 3] * hv3;
    const g1 = tab[i10] * hv0 + tab[i10 + 2] * hv1 + tab[i11] * hv2 + tab[i11 + 2] * hv3;
    const gp1 = tab[i10 + 1] * hv0 + tab[i10 + 3] * hv1 + tab[i11 + 1] * hv2 + tab[i11 + 3] * hv3;
    const F = g0 * hu0 + gp0 * hu1 + g1 * hu2 + gp1 * hu3;
    const Fu = g0 * du0 + gp0 * du1 + g1 * du2 + gp1 * du3;
    // ∂F/∂r|τ = F̂_p / r_p − F̂_w (∂τ/∂r|w) / (∂τ/∂w),  τ = t2(r) + (t1(r) − t2(r))·m(w)
    let Fr = (Fu * n) / rp;
    if (!saturated) {
      const dv0 = 6 * v2 - 6 * v;
      const dv1 = 3 * v2 - 4 * v + 1;
      const dv2 = -dv0;
      const dv3 = 3 * v2 - 2 * v;
      const g0v = tab[i00] * dv0 + tab[i00 + 2] * dv1 + tab[i01] * dv2 + tab[i01 + 2] * dv3;
      const gp0v = tab[i00 + 1] * dv0 + tab[i00 + 3] * dv1 + tab[i01 + 1] * dv2 + tab[i01 + 3] * dv3;
      const g1v = tab[i10] * dv0 + tab[i10 + 2] * dv1 + tab[i11] * dv2 + tab[i11 + 2] * dv3;
      const gp1v = tab[i10 + 1] * dv0 + tab[i10 + 3] * dv1 + tab[i11 + 1] * dv2 + tab[i11 + 3] * dv3;
      const Fv = g0v * hu0 + gp0v * hu1 + g1v * hu2 + gp1v * hu3;
      const mw = mapBoth(w);
      const tauR = t2 > 0 ? (re / t2) * (1 - mw) + (re / t1) * mw : (re / t1) * mw;
      Fr -= ((Fv * nw) / ((t1 - t2) * mapBothD(w))) * tauR;
    }
    const Y =
      (tab[i00 + 4] * hv0 + tab[i00 + 6] * hv1 + tab[i01 + 4] * hv2 + tab[i01 + 6] * hv3) * hu0 +
      (tab[i00 + 5] * hv0 + tab[i00 + 7] * hv1 + tab[i01 + 5] * hv2 + tab[i01 + 7] * hv3) * hu1 +
      (tab[i10 + 4] * hv0 + tab[i10 + 6] * hv1 + tab[i11 + 4] * hv2 + tab[i11 + 6] * hv3) * hu2 +
      (tab[i10 + 5] * hv0 + tab[i10 + 7] * hv1 + tab[i11 + 5] * hv2 + tab[i11 + 7] * hv3) * hu3;
    this.sF = F;
    this.sFr = Fr;
    this.sY = Y;
  }

  /** Exact F, ∂F/∂r|τ, Y (slow; beyond-table fallback). */
  private lookupExact(r: number, tau: number): void {
    const out = this.qAcc;
    out.fill(0);
    this.integrateFY(r, 0, tau, out, false);
    this.sF = out[0];
    this.sFr = out[1];
    this.sY = out[2];
  }

  /**
   * Slice quantities in the lens range for P1 = ρ² − c1², P2 = c2² − ρ² (both > 0): writes
   * o = [Dev, q = Dev′(ρ)/ρ, Ψ, ℓ].
   */
  private devAt(P1: number, P2: number, o: Float64Array): void {
    const R = this.radius;
    const c1 = this.c1;
    const ang = this.ang;
    const l = this.lensAngles(P1, P2, ang);
    const rho2 = c1 * c1 + P1;
    const Aout = R * R * seg(ang[3]) - rho2 * seg(ang[0]); // πR² − L
    const s = P2 / this.dc2;
    const omw = s * s; // 1 − w
    if (this.outside) {
      // Dev = πR² w − L = w·A_out − (1 − w)·L;  Dev′/ρ = πR² w′/ρ − Φ = 4πR² P2/Δ² − Φ
      const L = rho2 * seg(ang[0]) + R * R * seg(ang[2]);
      o[0] = (1 - omw) * Aout - omw * L;
      o[1] = (4 * PI * R * R * P2) / (this.dc2 * this.dc2) - 2 * ang[0];
      o[2] = 2 * ang[2];
      o[3] = l;
      return;
    }
    const E = rho2 * seg(ang[1]) - R * R * seg(ang[2]); // πρ² − L
    o[0] = (1 - omw) * Aout + omw * E;
    // Dev′/ρ = w′π(R² − ρ²)/ρ + (1 − w)2π − Φ,  w′ = 4ρP2/Δ²
    // R² − ρ² = d(2R − d) − P1 (no R² − ρ² cancellation for small d)
    const d = this.offset;
    o[1] = (4 * PI * P2 * (d * (2 * R - d) - P1)) / (this.dc2 * this.dc2) + TWO_PI * omw - 2 * ang[0];
    o[2] = 2 * ang[2];
    o[3] = l;
  }

  /**
   * Accumulates into out[0..3] the integrals over t ∈ [tFrom, min(tTo, t₁)] (0 ≤ tFrom) of
   * [Dev, r·q, RΨ, RΨ′ r/ρ] — the integrands of F, ∂F/∂r|τ, Y, ∂Y/∂r|τ. `withYr` = false skips the
   * last one (log-singular at r = R + d).
   */
  private integrateFY(r: number, tFrom: number, tTo: number, out: Float64Array, withYr: boolean): void {
    const R = this.radius;
    const d = this.offset;
    const c1 = this.c1;
    const c2 = this.c2;
    if (r <= c1) return;
    if (this.outside && r <= this.lo) return;
    const t1 = Math.sqrt((r - c1) * (r + c1));
    const t2 = r > c2 ? Math.sqrt((r - c2) * (r + c2)) : 0;
    const hiT = Math.min(tTo, t1);
    if (!(hiT > tFrom)) return;
    // ρ ≥ c2 for t ≤ t2: Dev = 0, q = 0, Ψ = 2π, Ψ′ = 0
    const aEnd = Math.min(hiT, t2);
    if (aEnd > tFrom) out[2] += TWO_PI * R * (aEnd - tFrom);
    const lo = Math.max(tFrom, t2);
    if (!(hiT > lo)) return;
    const acc = new Float64Array(NC_MAX);
    const tol = new Float64Array(NC_MAX);
    // absolute tolerances per unit of the substituted variable, scaled by full-τ magnitudes
    tol[0] = 1e-15 * PI * R * R * t1;
    tol[1] = 1e-15 * TWO_PI * r * t1;
    tol[2] = 1e-15 * TWO_PI * R * t1;
    tol[3] = (1e-11 * TWO_PI * r * t1) / d;
    const t1mHi = t1 - hiT;
    const loMt2 = lo - t2;
    const f: VecIntegrand = (t, dLo, dHi, o) => {
      const P1 = (t1mHi + dHi) * (t1 + t); // ρ² − c1²
      const P2 = r >= c2 ? (loMt2 + dLo) * (t + t2) : t * t + (c2 - r) * (c2 + r); // c2² − ρ²
      this.devAt(P1, P2, o);
      const l = o[3];
      o[1] = r * o[1];
      o[2] = R * o[2];
      o[3] = withYr && l > 0 ? (2 * r * R) / (d * l) : 0; // RΨ′(ρ) r/ρ = 2r/(d sin β), sin β = ℓ/R
    };
    integrateCosSub(f, withYr ? 4 : 3, lo, hiT, acc, tol, 1e-12);
    out[0] += acc[0];
    out[1] += acc[1];
    out[2] += acc[2];
    out[3] += acc[3];
  }

  /**
   * Fill one table segment (1: r ∈ [c1, c2], clustered at both ends; 2: r ∈ [c2, r_max], clustered at
   * c2). Node (p, w) ↦ r = r(p), τ = t1(r)·m(w). Stored per node: F̂, F̂_p·hp, F̂_w·hw, F̂_pw·hp·hw and the
   * same for Ŷ, with exact analytic derivatives (chain rule through the maps).
   */
  private buildTable(seg: 1 | 2): void {
    const tab = seg === 1 ? this.tab1 : this.tab2;
    const n = seg === 1 ? this.n1 : this.n2;
    const nw = this.nw;
    const hp = 1 / n;
    const hw = 1 / nw;
    const R = this.radius;
    const d = this.offset;
    const c1 = this.c1;
    const c2 = this.c2;
    const D1 = this.delta1;
    const D2 = this.delta2;
    const I = new Float64Array(NC_MAX);
    const o = new Float64Array(NC_MAX);
    for (let i = 0; i <= n; i++) {
      const p = i * hp;
      let r: number;
      let rp: number;
      let t1: number;
      let rpT1r: number; // (dr/dp)·(∂t1/∂r), regularised
      let t2 = 0;
      let rpT2r = 0; // (dr/dp)·(∂t2/∂r), regularised
      if (seg === 1) {
        if (i === 0) continue; // r = lo: F ≡ Y ≡ 0 with all derivatives (F ∝ (r−lo)², Y ∝ (r−lo)); zero-filled
        const lo = this.lo; // = c1 (spark inside the bore), e1 = d − R (outside)
        const sig = mapBoth(p);
        r = lo + D1 * sig;
        rp = D1 * mapBothD(p);
        t1 = Math.sqrt(D1 * sig * (r + lo));
        // for p ≤ ½: σ′/√σ = 4p/(√2 p) = 2√2 exactly
        rpT1r = p <= 0.5 ? (r * D1 * 2 * Math.SQRT2) / Math.sqrt(D1 * (r + lo)) : (rp * r) / t1;
      } else {
        const dr = D2 * p * p;
        r = c2 + dr;
        rp = 2 * D2 * p;
        t1 = Math.sqrt((r - c1) * (r + c1));
        rpT1r = (rp * r) / t1;
        t2 = Math.sqrt(dr * (r + c2));
        rpT2r = (2 * r * D2) / Math.sqrt(D2 * (r + c2)); // 2D2 p · r / (p √(D2 (r + c2)))
      }
      const needYr = rp > 0;
      I.fill(0);
      // the part τ ≤ t2 is closed form (Dev = 0, Ψ = 2π): start the cumulative integrals there
      I[2] = TWO_PI * R * t2;
      let tPrev = t2;
      for (let j = 0; j <= nw; j++) {
        let w = j * hw;
        if (w < EPS_NODE) w = EPS_NODE;
        if (w > 1 - EPS_NODE) w = 1 - EPS_NODE;
        const mw = mapBoth(w);
        const mwD = mapBothD(w);
        const span = t1 - t2;
        const tau = t2 + span * mw;
        // integrate the last piece exactly to t1 so the 1/√ singularity of Ψ′ sits on an endpoint
        const tauInt = j === nw ? t1 : tau;
        this.integrateFY(r, tPrev, tauInt, I, needYr);
        tPrev = tauInt;
        const tauW = span * mwD;
        // slice quantities at τ (P1, P2 to full relative precision)
        const P1 = span * mapBothC(w) * (t1 + tau); // t1² − τ²
        const P2 = seg === 2 ? span * mw * (tau + t2) : tau * tau + (c2 - r) * (c2 + r); // c2² − ρ²
        const rho = Math.sqrt(c1 * c1 + P1);
        let Ftau = 0;
        let qTau = 0;
        let Ytau = TWO_PI * R;
        let psiD = 0; // dΨ/dρ
        if (P2 > 0) {
          this.devAt(P1, P2, o);
          Ftau = o[0];
          qTau = o[1];
          Ytau = R * o[2];
          const sb = Math.max(o[3] / R, 1e-12);
          psiD = (2 * rho) / (d * R * sb);
        }
        const k = (i * (nw + 1) + j) * NODE;
        const rpTauR = rpT2r * (1 - mw) + rpT1r * mw; // r_p ∂τ/∂r|w
        const rpTauRW = (rpT1r - rpT2r) * mwD; // r_p ∂²τ/∂r∂w
        tab[k] = I[0];
        tab[k + 1] = (I[1] * rp + Ftau * rpTauR) * hp;
        tab[k + 2] = Ftau * tauW * hw;
        tab[k + 3] = (qTau * tauW * (rp * r - tau * rpTauR) + Ftau * rpTauRW) * hp * hw;
        tab[k + 4] = I[2];
        tab[k + 5] = ((needYr ? I[3] * rp : 0) + Ytau * rpTauR) * hp;
        tab[k + 6] = Ytau * tauW * hw;
        tab[k + 7] = ((R * psiD * tauW * (rp * r - tau * rpTauR)) / rho + Ytau * rpTauRW) * hp * hw;
      }
    }
  }
}
