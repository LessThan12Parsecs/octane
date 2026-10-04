/**
 * Combustion-chamber geometry behind one interface ({@link CombustionChamber}) for every burn / heat
 * geometry call of the cycle model: the spherical flame of radius r centred at the spark gap intersected
 * with the chamber at the current piston position (volume, front area = ∂V/∂r, burned-wetted wall areas),
 * its inverse r(V), the chamber volume and wall areas, and the turbulence length-scale height.
 *
 *  • 'flat-disc' — {@link DiscChamber}: an adapter over FlameGeometry (sphere ∩ disc) and
 *    heat-transfer/wall-heat flatChamberAreas that calls the same code in the same float order, so the
 *    CFR F-1 is bit-identical (chamber.test.ts compares every quantity with toBe).
 *  • 'l-head'   — {@link LHeadChamber}: the side-valve chamber of LHeadChamberSpec (core/engine-spec.ts),
 *    cylinder frame with y = 0 at the roof of the head cavity over the bore:
 *      bore column  D × [−h(θ), 0],          D = { x² + z² ≤ R² }   (moving piston crown at y = −h)
 *      valve pocket Ω × [deckY, roofY],      Ω = rounded rectangle ∖ D (fixed prism beside the bore)
 *    The two plan regions are disjoint, so for the ball B(c, r) (c = spark gap, anywhere in the chamber —
 *    the Model T plug sits over the pocket, OUTSIDE the bore planform)
 *      V(r, h)  = |B ∩ column|(r, h) + |B ∩ prism|(r)
 *      A_f(r, h) = ∂V/∂r = the same sum of front areas.
 *    The column term is FlameGeometry (with its d > R branch when the spark is beside the bore). The
 *    prism term is a 1-D function of r only: slicing the ball by planes y = c_y + t (Archimedes' hat-box
 *    theorem, as in flame-geometry.ts) gives V_Ω(r) = ∫_{t_floor}^{t_roof} a_Ω(√(r² − t²)) dt with
 *    a_Ω(ρ) = |disc(c_xz, ρ) ∩ Ω|, and A_Ω = ∂V_Ω/∂r = ∫ 2r a_Ω′(u) dt (u = ρ², a_Ω′(u) = ℓ_Ω(ρ)/(2ρ),
 *    ℓ_Ω = length of the slice circle inside Ω: the co-area identity dA/dρ = arc length holds for any
 *    planform). a_Ω, ℓ_Ω and the wall lengths inside the slice circle are exact closed forms (Green's
 *    theorem ½∮(x dz − z dx) over the boundary of disc ∩ Ω: circle arcs, straight edges and arcs of the
 *    rounded rectangle, the bore arc).
 *
 * ── Tables (built once per engine, immutable, shared by identical cylinders; ≈ 50–150 ms) ───────────
 *  • plan profile g(u) = a_Ω(√u) on u = ρ² ∈ [0, ρ_max²]: cubic Hermite with the EXACT derivative
 *    ℓ/(2ρ) at the nodes; nodes at every critical radius of Ω about the spark (vertices, edge feet,
 *    nearest/farthest points of the arcs — the radii where the slice circle's crossings with ∂Ω change),
 *    refined by bisection until g and g′ match the closed forms at the ¼, ½, ¾ points (1e-9 A_Ω, 1e-5 π).
 *    The two wall lengths (pocket side walls, open bore arc) are interpolated linearly on the same nodes.
 *  • layer table V_Ω(r), A_Ω(r): cubic Hermite in r whose node data are the layer integrals of the
 *    profile, EXACT for the piecewise profile (on each profile cell the integrand is a polynomial of degree
 *    ≤ 6 in t, integrated by 4-point Gauss–Legendre), with nodes at every kink radius r = √(ρ_k² + τ²)
 *    (τ = the floor, roof and — when the spark lies inside the layer — the spark plane) and bisection until
 *    V and A match at mid-cell (1e-9 V_Ω, 1e-5 of the layer's circumference × height) and each cell is
 *    monotone (Fritsch & Carlson 1980, SIAM J. Numer. Anal. 17, 238–246, sufficient condition
 *    α² + β² ≤ 9). A_f is the r-derivative of the interpolated V, so
 *    dV/dr = A_f holds to rounding (Newton inversion, burn-out invariant), V(r ≥ r_max) is EXACTLY
 *    chamberVolume(h) = πR²h + V_fixed (the last node stores V_fixed = A_Ω·(roofY − deckY) itself), and
 *    V(r ≤ r_inscribed) = 4/3 π r³ to rounding (a cubic: reproduced exactly by the Hermite cell).
 *  • Accuracy (chamber.test.ts against tools/reference/combustion_geometry_lhead_mc.py — Monte Carlo
 *    and independent dense quadrature, four chambers / spark positions): V ≲ 3e-7, A_f ≲ 2e-4 (the bore
 *    column's FlameGeometry table), wetted areas ≲ 6e-4 relative.
 *
 * ── Wall surfaces (heat-transfer/wall-heat WALL_* order) of the L-head ─────────────────────────────
 *  head   = cavity roof over the bore (πR², y = 0) + pocket roof (A_Ω, y = roofY) + bore-circle wall of
 *           the head cavity above max(deckY, −h) except where it is open to the pocket (deckY…roofY,
 *           the bore arc inside the rounded rectangle) + pocket side walls (deckY…roofY)
 *  block  = pocket floor (A_Ω at y = deckY) minus the valve faces
 *  piston = crown (πR²) + its cylindrical side facing the pocket while the crown is above the deck
 *  liner  = bore wall below the deck, 2πR·max(0, h + deckY)
 *  valves = count · π D_head²/4, in the pocket floor
 *  Burned-wetted areas: planes are slice lenses (bore disc: closed form; Ω: the profile), the valve faces
 *  exact lenses at the floor plane, the bore walls from FlameGeometry's liner integral at h and at the
 *  deck height, the side walls / open arc from the layer table (GL-4 quadrature for the piston-side
 *  band while the crown is above the deck). The crevice-mouth burned fraction of the L-head is the
 *  share of the bore circle at the crown plane inside the ball.
 *
 * Limitations: the flame is a Euclidean sphere — in the non-convex L-head it "sees" the column through
 * the deck strip and the head metal above the pocket roof (the pocket burns slightly early); the true
 * fix is a geodesic flame (UNVERIFIED magnitude). Turbulence length-scale height: Poulos & Heywood (1983,
 * SAE 830334, via Xu & Filipi 2020 eq. 32; turbulence.ts) use the disc chamber's h = 4V/(πB²), its volume
 * over its plan area; meanDepth(h) = chamberVolume(h)/plan area generalises that (= h for the disc;
 * UNVERIFIED for side-valve chambers — at the Model T's TDC 20.5 mm against the 34 mm of a bore-sized disc
 * of the same volume).
 *
 * Hot paths (evaluate / volume / frontArea / radiusForVolume) allocate nothing; each instance owns its
 * scratch, {@link CombustionChamber.clone} gives another evaluator over the same immutable tables.
 */
import type { EngineGeometrySpec, EngineSpec, LHeadChamberSpec, ValveSpec } from '../core/engine-spec';
import type { SliderCrankGeometry } from '../mechanics/kinematics';
import {
  N_WALL_SURFACES,
  WALL_BLOCK,
  WALL_EXHAUST_VALVE,
  WALL_HEAD,
  WALL_INTAKE_VALVE,
  WALL_LINER,
  WALL_PISTON,
  flatChamberAreas,
  newChamberAreas,
  type ChamberAreas,
} from '../heat-transfer/wall-heat';
import {
  FlameGeometry,
  arcOfBoreInside,
  lensArea,
  newFlameGeometryResult,
  type FlameGeometryOptions,
  type FlameGeometryResult,
} from './flame-geometry';

const PI = Math.PI;
const TWO_PI = 2 * Math.PI;

/** 4-point Gauss–Legendre on [−1, 1]. */
const GL4X = [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526];
const GL4W = [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538];

// ---------------------------------------------------------------------------------------------
// Interface
// ---------------------------------------------------------------------------------------------

/** Chamber shapes (EngineGeometrySpec.chamber). */
export type ChamberKind = EngineGeometrySpec['chamber'];

/** Result of {@link CombustionChamber.evaluate}. All SI. */
export interface ChamberResult {
  /** Enflamed volume (ball ∩ chamber), m³. */
  volume: number;
  /** Flame-front area inside the chamber, m² (= ∂volume/∂r). */
  frontArea: number;
  /** Burned-wetted area of each wall surface inside the ball, m² (wall-heat WALL_* order). */
  wetted: Float64Array;
  /**
   * Burned fraction of each surface ∈ [0, 1] for wallHeatLossTwoZoneSurfaces, with respect to
   * surfaceAreas(h). Flat disc: the legacy split (both valve faces at the head-disc fraction).
   */
  burnedFraction: Float64Array;
  /** Burned fraction of the gas at the top-land crevice mouth ∈ [0, 1] (crevice-zone inflow). */
  creviceBurnedFraction: number;
}

/** Allocate a zeroed {@link ChamberResult} (once, outside hot loops). */
export function newChamberResult(): ChamberResult {
  return {
    volume: 0,
    frontArea: 0,
    wetted: new Float64Array(N_WALL_SURFACES),
    burnedFraction: new Float64Array(N_WALL_SURFACES),
    creviceBurnedFraction: 0,
  };
}

/**
 * Chamber geometry for the flame, turbulence and wall-heat models. `h` is the clearance height of
 * SliderCrank (flat disc: head-to-crown height; L-head: depth of the bore column below the head-cavity
 * roof, with SliderCrankGeometry.fixedChamberVolume = {@link fixedVolume}).
 */
export interface CombustionChamber {
  readonly kind: ChamberKind;
  /** Cylinder bore, m. */
  readonly bore: number;
  /** Spark-gap centre (flame centre), cylinder frame, m. */
  readonly sparkCenter: readonly [number, number, number];
  /** Chamber volume outside the bore column (L-head valve pocket; 0 for the disc), m³. */
  readonly fixedVolume: number;
  /** Plan area of the chamber (bore disc + pocket), m². */
  readonly planformArea: number;
  /** The bore-column sphere ∩ disc geometry (snapshot / legacy consumers). */
  readonly flame: FlameGeometry;
  /** Volume, front area, burned-wetted surfaces and fractions at flame radius r and height h; returns `out`. */
  evaluate(r: number, h: number, out: ChamberResult): ChamberResult;
  /** Enflamed volume only, m³. */
  volume(r: number, h: number): number;
  /** Flame-front area only, m². */
  frontArea(r: number, h: number): number;
  /** Flame radius r with volume(r, h) = V (safeguarded Newton on dV/dr = A_f), m; 0 for V ≤ 0, maxRadius(h) for V ≥ chamberVolume(h). */
  radiusForVolume(V: number, h: number, rGuess?: number): number;
  /** Distance from the spark to the farthest chamber point at height h, m. */
  maxRadius(h: number): number;
  /** Radius below which the ball lies entirely inside the gas (V = 4/3 π r³), m (a lower bound). */
  inscribedRadius(h: number): number;
  /** Flame-reachable chamber volume A_p·h + fixedVolume = V − V_crevice, m³. */
  chamberVolume(h: number): number;
  /** Turbulence length-scale height chamberVolume(h)/planformArea, m (= h for the disc). */
  meanDepth(h: number): number;
  /** Wall surface areas at height h, m² (wall-heat WALL_* order); returns `out`. */
  surfaceAreas(h: number, out: Float64Array): Float64Array;
  /** Another evaluator sharing this chamber's immutable tables (one per cylinder). */
  clone(): CombustionChamber;
}

/** Fast-table options of {@link createChamber} (FlameGeometryOptions of the bore column). */
export interface ChamberOptions {
  /**
   * Largest clearance height h the fast tables must cover, m. Cycle model: (pistonTravel + h_TDC(CR_min))
   * × 1.02. Default: 2.5 × bore (flat disc, as FlameGeometry); stroke + h_TDC with 2 % margin (L-head).
   */
  maxHeight?: number;
  cellsInner?: number;
  cellsOuter?: number;
  cellsTau?: number;
  cornerBand?: number;
  /** Reuse identical immutable tables built earlier in this process (default true; numerically invisible). */
  cache?: boolean;
}

/** The spec fields a chamber is built from. */
export type ChamberSpec = Pick<EngineSpec, 'geometry' | 'sparkPlug' | 'intakeValve' | 'exhaustValve'>;

const clamp01 = (f: number): number => (f > 0 ? (f < 1 ? f : 1) : 0);

// ---------------------------------------------------------------------------------------------
// Flat disc: adapter over FlameGeometry + flatChamberAreas (bit-identical to the legacy calls)
// ---------------------------------------------------------------------------------------------

/**
 * 'flat-disc' chamber: FlameGeometry(bore, gap) and the 5-surface flat split, with the cycle model's own
 * expressions (chamberVolume = (πB²/4)·h as SliderCrank.boreArea·h, the clampRadiusGuess inscribed bound,
 * the head-disc burned fraction for both valve faces, the liner fraction at the crevice mouth).
 */
export class DiscChamber implements CombustionChamber {
  readonly kind = 'flat-disc' as const;
  readonly bore: number;
  readonly sparkCenter: readonly [number, number, number];
  readonly fixedVolume = 0;
  readonly planformArea: number;
  readonly flame: FlameGeometry;
  /** πB²/4 written as SliderCrank.boreArea. */
  readonly boreArea: number;
  private readonly intake: Pick<ValveSpec, 'count' | 'headDiameter'>;
  private readonly exhaust: Pick<ValveSpec, 'count' | 'headDiameter'>;
  private readonly fr: FlameGeometryResult = newFlameGeometryResult();
  private readonly areas: ChamberAreas = newChamberAreas();

  /**
   * @param flame the bore's FlameGeometry (built with the gap centre and table options)
   * @param intake intake valve (count, head diameter) — head-disc split
   * @param exhaust exhaust valve
   */
  constructor(flame: FlameGeometry, intake: Pick<ValveSpec, 'count' | 'headDiameter'>, exhaust: Pick<ValveSpec, 'count' | 'headDiameter'>) {
    this.flame = flame;
    this.bore = flame.bore;
    this.sparkCenter = flame.center;
    this.boreArea = (Math.PI * flame.bore * flame.bore) / 4;
    this.planformArea = this.boreArea;
    this.intake = intake;
    this.exhaust = exhaust;
  }

  clone(): DiscChamber {
    return new DiscChamber(this.flame.clone(), this.intake, this.exhaust);
  }

  evaluate(r: number, h: number, out: ChamberResult): ChamberResult {
    const g = this.flame.evaluate(r, h, this.fr);
    out.volume = g.volume;
    out.frontArea = g.frontArea;
    // wallHeatLossTwoZone's fractions, expression by expression
    const a = flatChamberAreas(this.bore, h, this.intake, this.exhaust, this.areas);
    const headDisc = a.head + a.intakeValves + a.exhaustValves;
    let fh = headDisc > 0 ? g.wettedHead / headDisc : 0;
    fh = fh > 0 ? (fh < 1 ? fh : 1) : 0;
    let fp = a.piston > 0 ? g.wettedPiston / a.piston : 0;
    fp = fp > 0 ? (fp < 1 ? fp : 1) : 0;
    let fl = a.liner > 0 ? g.wettedLiner / a.liner : 0;
    fl = fl > 0 ? (fl < 1 ? fl : 1) : 0;
    const f = out.burnedFraction;
    f[WALL_HEAD] = fh;
    f[WALL_PISTON] = fp;
    f[WALL_LINER] = fl;
    f[WALL_INTAKE_VALVE] = fh;
    f[WALL_EXHAUST_VALVE] = fh;
    f[WALL_BLOCK] = 0;
    const w = out.wetted;
    w[WALL_HEAD] = fh * a.head;
    w[WALL_PISTON] = g.wettedPiston;
    w[WALL_LINER] = g.wettedLiner;
    w[WALL_INTAKE_VALVE] = fh * a.intakeValves;
    w[WALL_EXHAUST_VALVE] = fh * a.exhaustValves;
    w[WALL_BLOCK] = 0;
    // cycle model (crevice zone): f_b = wettedLiner / liner where the liner is exposed
    out.creviceBurnedFraction = fl;
    return out;
  }

  volume(r: number, h: number): number {
    return this.flame.volume(r, h);
  }

  frontArea(r: number, h: number): number {
    return this.flame.frontArea(r, h);
  }

  radiusForVolume(V: number, h: number, rGuess?: number): number {
    return this.flame.radiusForVolume(V, h, rGuess);
  }

  maxRadius(h: number): number {
    return this.flame.maxRadius(h);
  }

  inscribedRadius(h: number): number {
    const fg = this.flame;
    const b = fg.headDistance;
    return h > b ? Math.min(b, h - b, fg.radius - fg.offset) : 0;
  }

  chamberVolume(h: number): number {
    return this.boreArea * h;
  }

  meanDepth(h: number): number {
    return h;
  }

  surfaceAreas(h: number, out: Float64Array): Float64Array {
    const a = flatChamberAreas(this.bore, h, this.intake, this.exhaust, this.areas);
    out[WALL_HEAD] = a.head;
    out[WALL_PISTON] = a.piston;
    out[WALL_LINER] = a.liner;
    out[WALL_INTAKE_VALVE] = a.intakeValves;
    out[WALL_EXHAUST_VALVE] = a.exhaustValves;
    out[WALL_BLOCK] = 0;
    return out;
  }
}

// ---------------------------------------------------------------------------------------------
// Plan geometry (build time): regions bounded by straight edges and circular arcs
// ---------------------------------------------------------------------------------------------

/** Wall tags of plan-boundary pieces. */
const TAG_SIDE = 0; // pocket side wall (head casting)
const TAG_OPEN = 1; // bore circle inside the pocket outline: open boundary between pocket and bore column

/**
 * Boundary piece of a plan region in the (x, z) plane, oriented with the region on its LEFT for the
 * positive sense (angle increasing from +x toward +z), so ½∮(x dz − z dx) is the region's area.
 * Lines run from (x0, z0) to (x1, z1); arcs have centre (cx, cz), radius s, start angle a0 and signed
 * sweep sw (their end points are stored in x0…z1 too).
 */
interface PlanPrim {
  line: boolean;
  x0: number;
  z0: number;
  x1: number;
  z1: number;
  cx: number;
  cz: number;
  s: number;
  a0: number;
  sw: number;
  tag: number;
}

function linePrim(x0: number, z0: number, x1: number, z1: number, tag: number): PlanPrim {
  return { line: true, x0, z0, x1, z1, cx: 0, cz: 0, s: 0, a0: 0, sw: 0, tag };
}

function arcPrim(cx: number, cz: number, s: number, a0: number, sw: number, tag: number): PlanPrim {
  const a1 = a0 + sw;
  return {
    line: false,
    x0: cx + s * Math.cos(a0),
    z0: cz + s * Math.sin(a0),
    x1: cx + s * Math.cos(a1),
    z1: cz + s * Math.sin(a1),
    cx,
    cz,
    s,
    a0,
    sw,
    tag,
  };
}

function primTranslate(p: PlanPrim, dx: number, dz: number): PlanPrim {
  return { ...p, x0: p.x0 + dx, z0: p.z0 + dz, x1: p.x1 + dx, z1: p.z1 + dz, cx: p.line ? 0 : p.cx + dx, cz: p.line ? 0 : p.cz + dz };
}

/** Point of the piece at parameter t ∈ [0, 1] into out[0..1]. */
function primPoint(p: PlanPrim, t: number, out: Float64Array): void {
  if (p.line) {
    out[0] = p.x0 + t * (p.x1 - p.x0);
    out[1] = p.z0 + t * (p.z1 - p.z0);
  } else {
    const a = p.a0 + t * p.sw;
    out[0] = p.cx + p.s * Math.cos(a);
    out[1] = p.cz + p.s * Math.sin(a);
  }
}

function primSub(p: PlanPrim, t0: number, t1: number): PlanPrim {
  if (p.line) {
    const ex = p.x1 - p.x0;
    const ez = p.z1 - p.z0;
    return linePrim(p.x0 + t0 * ex, p.z0 + t0 * ez, p.x0 + t1 * ex, p.z0 + t1 * ez, p.tag);
  }
  return arcPrim(p.cx, p.cz, p.s, p.a0 + t0 * p.sw, (t1 - t0) * p.sw, p.tag);
}

/** ½∫(x dz − z dx) along the piece over t ∈ [t0, t1] (Green's theorem area contribution), m². */
function primGreen(p: PlanPrim, t0: number, t1: number): number {
  if (p.line) {
    const ex = p.x1 - p.x0;
    const ez = p.z1 - p.z0;
    const ax = p.x0 + t0 * ex;
    const az = p.z0 + t0 * ez;
    const bx = p.x0 + t1 * ex;
    const bz = p.z0 + t1 * ez;
    return 0.5 * (ax * bz - az * bx);
  }
  // x = cx + s cos α, z = cz + s sin α: x dz − z dx = (s(cx cos α + cz sin α) + s²) dα
  const a0 = p.a0 + t0 * p.sw;
  const a1 = p.a0 + t1 * p.sw;
  const s = p.s;
  return 0.5 * (s * p.cx * (Math.sin(a1) - Math.sin(a0)) - s * p.cz * (Math.cos(a1) - Math.cos(a0)) + s * s * (a1 - a0));
}

function primLength(p: PlanPrim, t0: number, t1: number): number {
  if (p.line) return Math.hypot(p.x1 - p.x0, p.z1 - p.z0) * (t1 - t0);
  return p.s * Math.abs(p.sw) * (t1 - t0);
}

/** Parameter tolerance of the circle–piece crossings (extra crossings are harmless, missed ones are not). */
const T_TOL = 1e-9;

/** Parameter t ∈ [0, 1] of angle α on an arc (a0, sw), or NaN if α is not on the arc. */
function arcParam(alpha: number, a0: number, sw: number): number {
  let x = sw > 0 ? alpha - a0 : a0 - alpha;
  x -= TWO_PI * Math.floor(x / TWO_PI); // [0, 2π)
  const aw = Math.abs(sw);
  const t = x / aw;
  if (t <= 1 + T_TOL) return t < 1 ? t : 1;
  if ((x - TWO_PI) / aw >= -T_TOL) return 0;
  return Number.NaN;
}

/** Crossing parameters of the piece with the circle |q| = ρ about the coordinate origin (appended to out). */
function primCircleParams(p: PlanPrim, rho: number, out: number[]): void {
  if (p.line) {
    const ex = p.x1 - p.x0;
    const ez = p.z1 - p.z0;
    const A = ex * ex + ez * ez;
    const B = p.x0 * ex + p.z0 * ez;
    const C = p.x0 * p.x0 + p.z0 * p.z0 - rho * rho;
    const disc = B * B - A * C;
    if (!(disc >= 0) || !(A > 0)) return;
    const sq = Math.sqrt(disc);
    const q = -(B + (B >= 0 ? sq : -sq));
    const ta = q / A;
    const tb = q !== 0 ? C / q : ta;
    for (const t of [ta, tb]) if (t >= -T_TOL && t <= 1 + T_TOL) out.push(t < 0 ? 0 : t > 1 ? 1 : t);
    return;
  }
  const cm = Math.hypot(p.cx, p.cz);
  if (!(cm > 1e-15)) return; // arc concentric with the circle: no isolated crossings
  const k = (rho * rho - cm * cm - p.s * p.s) / (2 * p.s * cm);
  if (k > 1 || k < -1) return;
  const gam = Math.atan2(p.cz, p.cx);
  const del = Math.acos(k);
  for (const a of [gam - del, gam + del]) {
    const t = arcParam(a, p.a0, p.sw);
    if (!Number.isNaN(t)) out.push(t);
  }
}

/** Distances from the coordinate origin to the piece's end points and interior distance extrema (appended). */
function primCriticalRadii(p: PlanPrim, out: number[]): void {
  out.push(Math.hypot(p.x0, p.z0), Math.hypot(p.x1, p.z1));
  if (p.line) {
    const ex = p.x1 - p.x0;
    const ez = p.z1 - p.z0;
    const A = ex * ex + ez * ez;
    if (A > 0) {
      const tf = -(p.x0 * ex + p.z0 * ez) / A;
      if (tf > 0 && tf < 1) out.push(Math.hypot(p.x0 + tf * ex, p.z0 + tf * ez));
    }
    return;
  }
  const cm = Math.hypot(p.cx, p.cz);
  if (!(cm > 1e-15)) {
    out.push(p.s);
    return;
  }
  const tn = arcParam(Math.atan2(-p.cz, -p.cx), p.a0, p.sw); // nearest point of the full circle
  if (!Number.isNaN(tn)) out.push(Math.abs(cm - p.s));
  const tf = arcParam(Math.atan2(p.cz, p.cx), p.a0, p.sw); // farthest point
  if (!Number.isNaN(tf)) out.push(cm + p.s);
}

/** Plan region: boundary pieces (absolute cylinder-frame x, z) + an inside test; exact area and wall lengths. */
interface PlanRegion {
  prims: PlanPrim[];
  inside: (x: number, z: number) => boolean;
  /** Exact area (Green's theorem over the boundary), m². */
  area: number;
  /** Total length of the TAG_SIDE pieces (pocket side walls), m. */
  sideLength: number;
  /** Total length of the TAG_OPEN pieces (bore arc open to the pocket), m. */
  openLength: number;
}

type PocketPlanSpec = Pick<LHeadChamberSpec['pocket'], 'xMin' | 'xMax' | 'zMin' | 'zMax' | 'cornerRadius'>;

/** Rounded-rectangle inside test (corner radius clamped to the half-widths). */
function inRoundedRect(x: number, z: number, pk: PocketPlanSpec): boolean {
  const hx = 0.5 * (pk.xMax - pk.xMin);
  const hz = 0.5 * (pk.zMax - pk.zMin);
  const c = Math.min(pk.cornerRadius, hx, hz);
  const ax = Math.abs(x - 0.5 * (pk.xMax + pk.xMin));
  const az = Math.abs(z - 0.5 * (pk.zMax + pk.zMin));
  if (ax > hx || az > hz) return false;
  const ex = ax - (hx - c);
  const ez = az - (hz - c);
  if (ex > 0 && ez > 0) return ex * ex + ez * ez <= c * c;
  return true;
}

/** Counter-clockwise boundary of the rounded rectangle (tag TAG_SIDE). */
function roundedRectPrims(pk: PocketPlanSpec): PlanPrim[] {
  const { xMin: xa, xMax: xb, zMin: za, zMax: zb } = pk;
  const c = Math.min(pk.cornerRadius, 0.5 * (xb - xa), 0.5 * (zb - za));
  const out: PlanPrim[] = [];
  const edge = (x0: number, z0: number, x1: number, z1: number): void => {
    if (Math.hypot(x1 - x0, z1 - z0) > 0) out.push(linePrim(x0, z0, x1, z1, TAG_SIDE));
  };
  const corner = (cx: number, cz: number, a0: number): void => {
    if (c > 0) out.push(arcPrim(cx, cz, c, a0, 0.5 * PI, TAG_SIDE));
  };
  edge(xb, za + c, xb, zb - c);
  corner(xb - c, zb - c, 0);
  edge(xb - c, zb, xa + c, zb);
  corner(xa + c, zb - c, 0.5 * PI);
  edge(xa, zb - c, xa, za + c);
  corner(xa + c, za + c, PI);
  edge(xa + c, za, xb - c, za);
  corner(xb - c, za + c, 1.5 * PI);
  return out;
}

/**
 * Plan region Ω = rounded rectangle ∖ bore disc (L-head valve pocket): the rectangle's boundary outside
 * the disc (side walls) and the bore arc inside the rectangle traversed clockwise (open boundary).
 */
function pocketPlanRegion(bore: number, pk: PocketPlanSpec): PlanRegion {
  if (!(pk.xMax > pk.xMin && pk.zMax > pk.zMin && pk.cornerRadius >= 0)) throw new Error('L-head pocket: invalid outline');
  const R = 0.5 * bore;
  const R2 = R * R;
  const prims: PlanPrim[] = [];
  const angles: number[] = [];
  const ts: number[] = [];
  const q = new Float64Array(2);
  for (const p of roundedRectPrims(pk)) {
    ts.length = 0;
    primCircleParams(p, R, ts);
    ts.sort((a, b) => a - b);
    let t0 = 0;
    for (let k = 0; k <= ts.length; k++) {
      const t1 = k < ts.length ? ts[k] : 1;
      if (t1 - t0 > 1e-14) {
        primPoint(p, 0.5 * (t0 + t1), q);
        if (q[0] * q[0] + q[1] * q[1] > R2) prims.push(primSub(p, t0, t1));
      }
      if (t1 > t0) t0 = t1;
    }
    for (const t of ts) {
      primPoint(p, t, q);
      angles.push(Math.atan2(q[1], q[0]));
    }
  }
  if (angles.length === 0) {
    if (inRoundedRect(R, 0, pk)) prims.push(arcPrim(0, 0, R, TWO_PI, -TWO_PI, TAG_OPEN)); // bore inside the outline
  } else {
    angles.sort((a, b) => a - b);
    for (let i = 0; i < angles.length; i++) {
      const a = angles[i];
      const b = i + 1 < angles.length ? angles[i + 1] : angles[0] + TWO_PI;
      if (b - a > 1e-14) {
        const m = 0.5 * (a + b);
        if (inRoundedRect(R * Math.cos(m), R * Math.sin(m), pk)) prims.push(arcPrim(0, 0, R, b, a - b, TAG_OPEN));
      }
    }
  }
  let area = 0;
  let side = 0;
  let open = 0;
  for (const p of prims) {
    area += primGreen(p, 0, 1);
    if (p.tag === TAG_SIDE) side += primLength(p, 0, 1);
    else open += primLength(p, 0, 1);
  }
  if (!(area > 0)) throw new Error('L-head pocket: the outline lies inside the bore (empty valve pocket)');
  return { prims, inside: (x, z) => inRoundedRect(x, z, pk) && x * x + z * z > R2, area, sideLength: side, openLength: open };
}

/** Exact cut of a plan region by the disc of radius ρ about the coordinate origin (pieces given relative to it). */
interface CircleCut {
  area: number;
  arc: number;
  side: number;
  open: number;
}

function circleCut(prims: PlanPrim[], insideRel: (x: number, z: number) => boolean, rho: number, res: CircleCut): void {
  res.area = 0;
  res.arc = 0;
  res.side = 0;
  res.open = 0;
  const rho2 = rho * rho;
  const angles: number[] = [];
  const ts: number[] = [];
  const q = new Float64Array(2);
  for (const p of prims) {
    ts.length = 0;
    primCircleParams(p, rho, ts);
    ts.sort((a, b) => a - b);
    let t0 = 0;
    for (let k = 0; k <= ts.length; k++) {
      const t1 = k < ts.length ? ts[k] : 1;
      if (t1 > t0) {
        primPoint(p, 0.5 * (t0 + t1), q);
        if (q[0] * q[0] + q[1] * q[1] < rho2) {
          res.area += primGreen(p, t0, t1);
          if (p.tag === TAG_SIDE) res.side += primLength(p, t0, t1);
          else res.open += primLength(p, t0, t1);
        }
        t0 = t1;
      }
    }
    for (const t of ts) {
      primPoint(p, t, q);
      angles.push(Math.atan2(q[1], q[0]));
    }
  }
  if (angles.length === 0) {
    if (insideRel(rho, 0)) {
      res.area += PI * rho2;
      res.arc += TWO_PI * rho;
    }
    return;
  }
  angles.sort((a, b) => a - b);
  for (let i = 0; i < angles.length; i++) {
    const a = angles[i];
    const b = i + 1 < angles.length ? angles[i + 1] : angles[0] + TWO_PI;
    if (b - a > 1e-15) {
      const m = 0.5 * (a + b);
      if (insideRel(rho * Math.cos(m), rho * Math.sin(m))) {
        res.area += 0.5 * rho2 * (b - a);
        res.arc += rho * (b - a);
      }
    }
  }
}

// ---------------------------------------------------------------------------------------------
// Plan profile g(u) = |disc(c, √u) ∩ Ω| about the spark (+ wall lengths), u = ρ²
// ---------------------------------------------------------------------------------------------

/** Profile tolerances: g (relative to A_Ω), g′ = ℓ/(2ρ) (relative to π), wall lengths (relative to their total). */
const PROFILE_TOL_G = 1e-9;
const PROFILE_TOL_GP = 1e-5;
const PROFILE_TOL_W = 1e-5;
const MAX_DEPTH = 48;

/** Immutable plan profile (see the module comment). */
class PlanProfile {
  /** Node u = ρ², m² (u[0] = 0, last = ρ_max²). */
  readonly u: Float64Array;
  /** a_Ω(√u), m². */
  readonly g: Float64Array;
  /** dg/du = ℓ_Ω/(2ρ), dimensionless. */
  readonly gp: Float64Array;
  /** Side-wall length inside the circle, m. */
  readonly ws: Float64Array;
  /** Open-arc length inside the circle, m. */
  readonly wo: Float64Array;
  /** Critical u (kinks), m², sorted, incl. 0 and u_max. */
  readonly critical: Float64Array;
  readonly area: number;
  readonly side: number;
  readonly open: number;
  readonly uMax: number;
  private readonly n: number;

  constructor(region: PlanRegion, px: number, pz: number) {
    const rel = region.prims.map((p) => primTranslate(p, -px, -pz));
    const insideRel = (x: number, z: number): boolean => region.inside(x + px, z + pz);
    const P0inside = insideRel(0, 0);
    this.area = region.area;
    this.side = region.sideLength;
    this.open = region.openLength;
    const radii: number[] = [];
    for (const p of rel) primCriticalRadii(p, radii);
    const rhoMax = Math.max(...radii);
    const uMax = rhoMax * rhoMax;
    this.uMax = uMax;
    const crit = [0, ...radii.map((r) => r * r)].sort((a, b) => a - b);
    const cu: number[] = [];
    for (const v of crit) if (cu.length === 0 || v - cu[cu.length - 1] > 1e-13 * uMax) cu.push(v);
    cu[cu.length - 1] = uMax;
    this.critical = Float64Array.from(cu);
    const cut: CircleCut = { area: 0, arc: 0, side: 0, open: 0 };
    const exact = (u: number): [number, number, number, number] => {
      if (u >= uMax) return [region.area, 0, region.sideLength, region.openLength];
      if (!(u > 0)) return [0, P0inside ? PI : 0, 0, 0];
      const rho = Math.sqrt(u);
      circleCut(rel, insideRel, rho, cut);
      return [cut.area, cut.arc / (2 * rho), cut.side, cut.open];
    };
    const tolG = PROFILE_TOL_G * region.area;
    const tolGp = PROFILE_TOL_GP * PI;
    const tolW = PROFILE_TOL_W * (region.sideLength + region.openLength);
    const U: number[] = [];
    const F: [number, number, number, number][] = [];
    const ok = (ua: number, fa: number[], ub: number, fb: number[], s: number, fm: number[]): boolean => {
      const D = ub - ua;
      const s2 = s * s;
      const s3 = s2 * s;
      const H = (2 * s3 - 3 * s2 + 1) * fa[0] + (s3 - 2 * s2 + s) * D * fa[1] + (-2 * s3 + 3 * s2) * fb[0] + (s3 - s2) * D * fb[1];
      const Hp = (6 * s - 6 * s2) * ((fb[0] - fa[0]) / D) + (3 * s2 - 4 * s + 1) * fa[1] + (3 * s2 - 2 * s) * fb[1];
      // (fb − fa)/D carries the rounding of the stored values, ≈ ε|g|/D: below that no cell can do better
      const noise = (16 * Number.EPSILON * (Math.abs(fa[0]) + Math.abs(fb[0]))) / D;
      return (
        Math.abs(H - fm[0]) <= tolG &&
        Math.abs(Hp - fm[1]) <= tolGp + noise &&
        Math.abs(fa[2] + s * (fb[2] - fa[2]) - fm[2]) <= tolW &&
        Math.abs(fa[3] + s * (fb[3] - fa[3]) - fm[3]) <= tolW
      );
    };
    const refine = (ua: number, fa: [number, number, number, number], ub: number, fb: [number, number, number, number], depth: number): void => {
      const um = 0.5 * (ua + ub);
      const fm = exact(um);
      let good = depth >= MAX_DEPTH || ub - ua <= 1e-15 * uMax;
      if (!good) good = ok(ua, fa, ub, fb, 0.5, fm) && ok(ua, fa, ub, fb, 0.25, exact(ua + 0.25 * (ub - ua))) && ok(ua, fa, ub, fb, 0.75, exact(ua + 0.75 * (ub - ua)));
      if (good) return;
      refine(ua, fa, um, fm, depth + 1);
      U.push(um);
      F.push(fm);
      refine(um, fm, ub, fb, depth + 1);
    };
    let fPrev = exact(cu[0]);
    U.push(cu[0]);
    F.push(fPrev);
    for (let k = 1; k < cu.length; k++) {
      const fk = exact(cu[k]);
      refine(cu[k - 1], fPrev, cu[k], fk, 0);
      U.push(cu[k]);
      F.push(fk);
      fPrev = fk;
    }
    const n = U.length;
    this.n = n;
    this.u = Float64Array.from(U);
    this.g = Float64Array.from(F, (f) => f[0]);
    this.gp = Float64Array.from(F, (f) => f[1]);
    this.ws = Float64Array.from(F, (f) => f[2]);
    this.wo = Float64Array.from(F, (f) => f[3]);
  }

  /** Index i of the cell [u_i, u_{i+1}] containing x (u_0 ≤ x < u_max). */
  private cell(x: number): number {
    const u = this.u;
    let lo = 0;
    let hi = this.n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (u[mid] <= x) lo = mid;
      else hi = mid;
    }
    return lo;
  }

  /** g̃(x) (plane area of Ω inside the circle of radius √x), m². */
  value(x: number): number {
    if (!(x > 0)) return 0;
    if (x >= this.uMax) return this.area;
    const i = this.cell(x);
    const u = this.u;
    const D = u[i + 1] - u[i];
    const s = (x - u[i]) / D;
    const s2 = s * s;
    const s3 = s2 * s;
    return (2 * s3 - 3 * s2 + 1) * this.g[i] + (s3 - 2 * s2 + s) * D * this.gp[i] + (-2 * s3 + 3 * s2) * this.g[i + 1] + (s3 - s2) * D * this.gp[i + 1];
  }

  /** Open-arc length inside the circle of radius √x (linear), m. */
  openLength(x: number): number {
    if (!(x > 0)) return 0;
    if (x >= this.uMax) return this.open;
    const i = this.cell(x);
    const s = (x - this.u[i]) / (this.u[i + 1] - this.u[i]);
    return this.wo[i] + s * (this.wo[i + 1] - this.wo[i]);
  }

  /** out = [g̃, g̃′, side length, open length] at x. */
  lookup(x: number, out: Float64Array): void {
    if (!(x > 0)) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
      return;
    }
    if (x >= this.uMax) {
      out[0] = this.area;
      out[1] = 0;
      out[2] = this.side;
      out[3] = this.open;
      return;
    }
    const i = this.cell(x);
    const u = this.u;
    const D = u[i + 1] - u[i];
    const s = (x - u[i]) / D;
    const s2 = s * s;
    const s3 = s2 * s;
    const g0 = this.g[i];
    const g1 = this.g[i + 1];
    const m0 = D * this.gp[i];
    const m1 = D * this.gp[i + 1];
    out[0] = (2 * s3 - 3 * s2 + 1) * g0 + (s3 - 2 * s2 + s) * m0 + (-2 * s3 + 3 * s2) * g1 + (s3 - s2) * m1;
    out[1] = ((6 * s2 - 6 * s) * g0 + (3 * s2 - 4 * s + 1) * m0 + (6 * s - 6 * s2) * g1 + (3 * s2 - 2 * s) * m1) / D;
    out[2] = this.ws[i] + s * (this.ws[i + 1] - this.ws[i]);
    out[3] = this.wo[i] + s * (this.wo[i + 1] - this.wo[i]);
  }

  /** Number of nodes. */
  get size(): number {
    return this.n;
  }
}

// ---------------------------------------------------------------------------------------------
// Layer table: ball ∩ (Ω × [floor, roof]) as functions of r
// ---------------------------------------------------------------------------------------------

const LAYER_TOL_V = 1e-9;
const LAYER_TOL_A = 1e-5;
const LAYER_TOL_W = 1e-4;

/** Immutable layer table (see the module comment). */
class LayerTable {
  readonly r: Float64Array;
  readonly V: Float64Array;
  readonly A: Float64Array;
  readonly ws: Float64Array;
  readonly wo: Float64Array;
  /** Radius beyond which the ball contains the whole layer, m. */
  readonly rMax: number;
  /** Layer volume A_Ω·(t_roof − t_floor), m³ (exact last-node value). */
  readonly fullV: number;
  readonly fullSide: number;
  readonly fullOpen: number;
  private readonly n: number;

  /**
   * @param prof plan profile about the spark
   * @param tA floor-plane offset from the spark plane (y_floor − c_y), m
   * @param tB roof-plane offset (y_roof − c_y) > tA, m
   * @param fullV the exact layer volume A_Ω·(y_roof − y_floor), m³
   */
  constructor(prof: PlanProfile, tA: number, tB: number, fullV: number) {
    const H = tB - tA;
    this.fullV = fullV;
    this.fullSide = prof.side * H;
    this.fullOpen = prof.open * H;
    const tm = Math.max(tA * tA, tB * tB);
    const rMax = Math.sqrt(prof.uMax + tm);
    this.rMax = rMax;
    const crit = prof.critical;
    // kink radii: √(u_k + τ²) for the floor, roof and (inside the layer) the spark plane
    const taus = [tA, tB];
    if (tA < 0 && tB > 0) taus.push(0);
    const kr: number[] = [0];
    for (const u of crit) for (const t of taus) kr.push(Math.sqrt(u + t * t));
    kr.sort((a, b) => a - b);
    const K: number[] = [];
    for (const r of kr) if (r <= rMax && (K.length === 0 || r - K[K.length - 1] > 1e-12 * rMax)) K.push(r);
    K[K.length - 1] = rMax;
    // Layer integrals of the profile, EXACT for the piecewise profile: on a profile cell g̃ is a cubic in
    // u = r² − t², i.e. a degree-6 polynomial in t (g̃′: degree 4, the linear walls: degree 2), so 4-point
    // Gauss–Legendre per cell is exact. One-side integrals S(r, τ) = ∫₀^min(τ,r) [g̃, 2r g̃′, w_s, w_o] dt.
    const rhoMax = Math.sqrt(prof.uMax);
    const Aref = TWO_PI * rhoMax * H;
    const wallRef = this.fullSide + this.fullOpen;
    const U = prof.u;
    const G = prof.g;
    const GP = prof.gp;
    const WS = prof.ws;
    const WO = prof.wo;
    const nu = U.length;
    const S = new Float64Array(4);
    const side = (r: number, tau: number, sgn: number, acc: Float64Array): void => {
      const m = tau < r ? tau : r;
      if (!(m > 0)) return;
      const r2 = r * r;
      const uLo = (r - m) * (r + m);
      S.fill(0);
      // u ≥ u_max (t ≤ √(r² − u_max)): the whole plan region
      if (r2 > prof.uMax) {
        const ts = Math.min(m, Math.sqrt(r2 - prof.uMax));
        S[0] += prof.area * ts;
        S[2] += prof.side * ts;
        S[3] += prof.open * ts;
      }
      for (let i = 0; i + 1 < nu; i++) {
        const ua = U[i];
        const ub = U[i + 1];
        if (ub <= uLo || ua >= r2) continue;
        const lo = ua > uLo ? ua : uLo;
        const hi = ub < r2 ? ub : r2;
        if (!(hi > lo)) continue;
        const ta = Math.sqrt(r2 - hi);
        const tb = Math.sqrt(r2 - lo);
        const hc = 0.5 * (tb - ta);
        const mid = 0.5 * (tb + ta);
        const D = ub - ua;
        const g0 = G[i];
        const g1 = G[i + 1];
        const m0 = D * GP[i];
        const m1 = D * GP[i + 1];
        let sg = 0;
        let sa = 0;
        let ss = 0;
        let so = 0;
        for (let k = 0; k < 4; k++) {
          const t = mid + hc * GL4X[k];
          const s = ((r - t) * (r + t) - ua) / D;
          const s2 = s * s;
          const s3 = s2 * s;
          const wk = GL4W[k];
          sg += wk * ((2 * s3 - 3 * s2 + 1) * g0 + (s3 - 2 * s2 + s) * m0 + (-2 * s3 + 3 * s2) * g1 + (s3 - s2) * m1);
          sa += wk * ((6 * s2 - 6 * s) * g0 + (3 * s2 - 4 * s + 1) * m0 + (6 * s - 6 * s2) * g1 + (3 * s2 - 2 * s) * m1);
          ss += wk * (WS[i] + s * (WS[i + 1] - WS[i]));
          so += wk * (WO[i] + s * (WO[i + 1] - WO[i]));
        }
        S[0] += hc * sg;
        S[1] += (hc * sa * 2 * r) / D;
        S[2] += hc * ss;
        S[3] += hc * so;
      }
      for (let j = 0; j < 4; j++) acc[j] += sgn * S[j];
    };
    const acc = new Float64Array(4);
    const exact = (r: number): [number, number, number, number] => {
      if (!(r > 0)) return [0, 0, 0, 0];
      if (r >= rMax) return [fullV, 0, this.fullSide, this.fullOpen];
      acc.fill(0);
      // ∫_{tA}^{tB} = sgn(tB) S(|tB|) − sgn(tA) S(|tA|)
      side(r, Math.abs(tB), tB >= 0 ? 1 : -1, acc);
      side(r, Math.abs(tA), tA >= 0 ? -1 : 1, acc);
      return [acc[0], acc[1], acc[2], acc[3]];
    };
    const tolV = LAYER_TOL_V * fullV;
    const tolA = LAYER_TOL_A * Aref;
    const tolW = LAYER_TOL_W * wallRef;
    const Rn: number[] = [];
    const Fn: [number, number, number, number][] = [];
    const refine = (ra: number, fa: number[], rb: number, fb: number[], depth: number): void => {
      const rm = 0.5 * (ra + rb);
      const fm = exact(rm);
      const D = rb - ra;
      // Hermite at mid-cell: (fa + fb)/2 + D(f′a − f′b)/8; derivative 1.5 (fb − fa)/D − (f′a + f′b)/4
      const Hm = 0.5 * (fa[0] + fb[0]) + (D * (fa[1] - fb[1])) / 8;
      const Hpm = (1.5 * (fb[0] - fa[0])) / D - 0.25 * (fa[1] + fb[1]);
      const noise = (16 * Number.EPSILON * (Math.abs(fa[0]) + Math.abs(fb[0]))) / D; // rounding of (fb − fa)/D
      let good =
        Math.abs(Hm - fm[0]) <= tolV &&
        Math.abs(Hpm - fm[1]) <= tolA + noise &&
        Math.abs(0.5 * (fa[2] + fb[2]) - fm[2]) <= tolW &&
        Math.abs(0.5 * (fa[3] + fb[3]) - fm[3]) <= tolW;
      if (good) {
        // monotone cubic (Fritsch & Carlson 1980: α² + β² ≤ 9 with α, β = end slopes / secant)
        const sec = (fb[0] - fa[0]) / D;
        if (sec > 0) {
          const al = fa[1] / sec;
          const be = fb[1] / sec;
          good = al >= 0 && be >= 0 && al * al + be * be <= 9;
        } else good = fa[1] === 0 && fb[1] === 0;
      }
      if (good || depth >= MAX_DEPTH || D <= 1e-13 * rMax) return;
      refine(ra, fa, rm, fm, depth + 1);
      Rn.push(rm);
      Fn.push(fm);
      refine(rm, fm, rb, fb, depth + 1);
    };
    let fPrev = exact(K[0]);
    Rn.push(K[0]);
    Fn.push(fPrev);
    for (let k = 1; k < K.length; k++) {
      // start long kink intervals from a few equal cells (guards the mid-cell test)
      const m = Math.max(1, Math.ceil((K[k] - K[k - 1]) / (rMax / 24)));
      for (let j = 1; j <= m; j++) {
        const rj = j === m ? K[k] : K[k - 1] + ((K[k] - K[k - 1]) * j) / m;
        const fj = exact(rj);
        refine(Rn[Rn.length - 1], fPrev, rj, fj, 0);
        Rn.push(rj);
        Fn.push(fj);
        fPrev = fj;
      }
    }
    this.n = Rn.length;
    this.r = Float64Array.from(Rn);
    this.V = Float64Array.from(Fn, (v) => v[0]);
    this.A = Float64Array.from(Fn, (v) => v[1]);
    this.ws = Float64Array.from(Fn, (v) => v[2]);
    this.wo = Float64Array.from(Fn, (v) => v[3]);
    this.V[0] = 0;
    this.A[0] = 0;
    this.V[this.n - 1] = fullV;
    this.A[this.n - 1] = 0;
  }

  /** out = [V, A (= dV/dr of the interpolant), side-wall wetted, open-arc wetted] at r. Allocation-free. */
  lookup(r: number, out: Float64Array): void {
    if (!(r > 0)) {
      out[0] = 0;
      out[1] = 0;
      out[2] = 0;
      out[3] = 0;
      return;
    }
    if (r >= this.rMax) {
      out[0] = this.fullV;
      out[1] = 0;
      out[2] = this.fullSide;
      out[3] = this.fullOpen;
      return;
    }
    const x = this.r;
    let lo = 0;
    let hi = this.n - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (x[mid] <= r) lo = mid;
      else hi = mid;
    }
    const D = x[lo + 1] - x[lo];
    const s = (r - x[lo]) / D;
    const s2 = s * s;
    const s3 = s2 * s;
    const v0 = this.V[lo];
    const v1 = this.V[lo + 1];
    const m0 = D * this.A[lo];
    const m1 = D * this.A[lo + 1];
    out[0] = (2 * s3 - 3 * s2 + 1) * v0 + (s3 - 2 * s2 + s) * m0 + (-2 * s3 + 3 * s2) * v1 + (s3 - s2) * m1;
    out[1] = ((6 * s2 - 6 * s) * v0 + (3 * s2 - 4 * s + 1) * m0 + (6 * s - 6 * s2) * v1 + (3 * s2 - 2 * s) * m1) / D;
    out[2] = this.ws[lo] + s * (this.ws[lo + 1] - this.ws[lo]);
    out[3] = this.wo[lo] + s * (this.wo[lo + 1] - this.wo[lo]);
  }

  /** Number of nodes. */
  get size(): number {
    return this.n;
  }
}

// ---------------------------------------------------------------------------------------------
// L-head chamber
// ---------------------------------------------------------------------------------------------

/** Plan metrics of an L-head chamber (exact, Green's theorem). */
export interface LHeadPlanMetrics {
  /** Plan area of the valve pocket outside the bore disc A_Ω, m². */
  pocketArea: number;
  /** Length of the pocket outline outside the bore disc (side walls), m. */
  sideWallLength: number;
  /** Length of the bore circle inside the pocket outline (open to the pocket), m. */
  openArcLength: number;
  /** Gasket-face opening: bore disc + pocket, m². */
  planformArea: number;
  /** Fixed (pocket) chamber volume A_Ω·(roofY − deckY), m³. */
  fixedVolume: number;
}

/** Exact plan metrics of an L-head chamber spec. */
export function lHeadPlanMetrics(bore: number, lHead: LHeadChamberSpec): LHeadPlanMetrics {
  const reg = pocketPlanRegion(bore, lHead.pocket);
  const R = 0.5 * bore;
  return {
    pocketArea: reg.area,
    sideWallLength: reg.sideLength,
    openArcLength: reg.openLength,
    planformArea: PI * R * R + reg.area,
    fixedVolume: reg.area * (lHead.pocket.roofY - lHead.deckY),
  };
}

interface LHeadValve {
  count: number;
  /** Valve-head radius, m. */
  radius: number;
  /** Plan distance spark → valve centre, m. */
  dist: number;
  /** Face area count·π r², m². */
  area: number;
}

/** Immutable data and tables of an L-head chamber (shared by clones). */
class LHeadTables {
  readonly bore: number;
  readonly R: number;
  readonly sparkCenter: readonly [number, number, number];
  readonly deckY: number;
  readonly roofY: number;
  readonly region: PlanRegion;
  readonly profile: PlanProfile;
  readonly layer: LayerTable;
  readonly flame: FlameGeometry;
  readonly fixedVolume: number;
  readonly planformArea: number;
  /** Spark-plane offsets of the pocket floor and roof (y − c_y), m. */
  readonly tauFloor: number;
  readonly tauRoof: number;
  readonly sparkInPocket: boolean;
  /** Inscribed radius of a spark in the pocket, m. */
  readonly pocketInscribed: number;
  /** Squared critical radii of the open bore arc about the spark (kinks of its length inside a slice circle), m². */
  readonly openCritical: Float64Array;
  readonly intake: LHeadValve;
  readonly exhaust: LHeadValve;

  constructor(spec: ChamberSpec, opts: ChamberOptions) {
    const g = spec.geometry;
    const lh = g.lHead;
    if (!lh) throw new Error("L-head chamber: geometry.lHead is required for chamber 'l-head'");
    const bore = g.bore;
    const R = 0.5 * bore;
    const deckY = lh.deckY;
    const roofY = lh.pocket.roofY;
    if (!(deckY < 0 && roofY > deckY && roofY <= 0)) throw new Error('L-head chamber: needs deckY < roofY ≤ 0');
    this.bore = bore;
    this.R = R;
    this.deckY = deckY;
    this.roofY = roofY;
    const [cx, cy, cz] = spec.sparkPlug.gapCenter;
    this.sparkCenter = [cx, cy, cz];
    const region = pocketPlanRegion(bore, lh.pocket);
    this.region = region;
    this.fixedVolume = region.area * (roofY - deckY);
    this.planformArea = PI * R * R + region.area;
    const d = Math.hypot(cx, cz);
    const inBore = d < R && cy <= 0;
    const inPocket = region.inside(cx, cz) && cy >= deckY && cy <= roofY;
    if (!inBore && !inPocket) throw new Error('L-head chamber: the spark gap must lie inside the chamber (bore column or valve pocket)');
    this.sparkInPocket = !inBore;
    this.tauFloor = deckY - cy;
    this.tauRoof = roofY - cy;
    // the bore column (FlameGeometry; spark beside the bore → its d > R branch)
    const hTdc = -deckY - lh.crownAboveDeckAtTDC;
    const maxHeight = opts.maxHeight ?? (g.stroke + Math.max(hTdc, 0)) * 1.02;
    const fgOpts: FlameGeometryOptions = { maxHeight, allowSparkOutsideBore: true };
    if (opts.cellsInner !== undefined) fgOpts.cellsInner = opts.cellsInner;
    if (opts.cellsOuter !== undefined) fgOpts.cellsOuter = opts.cellsOuter;
    if (opts.cellsTau !== undefined) fgOpts.cellsTau = opts.cellsTau;
    if (opts.cornerBand !== undefined) fgOpts.cornerBand = opts.cornerBand;
    this.flame = new FlameGeometry(bore, [cx, cy, cz], fgOpts);
    // the pocket prism
    this.profile = new PlanProfile(region, cx, cz);
    this.layer = new LayerTable(this.profile, this.tauFloor, this.tauRoof, this.fixedVolume);
    // inscribed radius of a pocket spark: roof, floor, side walls, and the bore circle (beyond which the
    // column may be piston metal near TDC)
    let side = Infinity;
    const radii: number[] = [];
    for (const p of region.prims) {
      if (p.tag !== TAG_SIDE) continue;
      radii.length = 0;
      primCriticalRadii(primTranslate(p, -cx, -cz), radii);
      for (const v of radii) side = Math.min(side, v);
    }
    this.pocketInscribed = Math.max(0, Math.min(roofY - cy, cy - deckY, side, d - R));
    const oc: number[] = [0];
    for (const p of region.prims) {
      if (p.tag !== TAG_OPEN) continue;
      radii.length = 0;
      primCriticalRadii(primTranslate(p, -cx, -cz), radii);
      for (const v of radii) oc.push(v * v);
    }
    this.openCritical = Float64Array.from(oc.sort((a, b) => a - b));
    const valve = (v: ValveSpec, name: string): LHeadValve => {
      const rv = 0.5 * v.headDiameter;
      const [vx, vz] = v.position;
      const dv = Math.hypot(vx, vz);
      if (!(dv - rv >= R * (1 - 1e-12)) || !region.inside(vx, vz)) throw new Error(`L-head chamber: ${name} valve must lie in the pocket outside the bore`);
      for (let k = 0; k < 16; k++) {
        const a = (TWO_PI * k) / 16;
        if (!inRoundedRect(vx + rv * Math.cos(a), vz + rv * Math.sin(a), lh.pocket)) throw new Error(`L-head chamber: ${name} valve head extends beyond the pocket outline`);
      }
      return { count: v.count, radius: rv, dist: Math.hypot(vx - cx, vz - cz), area: v.count * PI * rv * rv };
    };
    this.intake = valve(spec.intakeValve, 'intake');
    this.exhaust = valve(spec.exhaustValve, 'exhaust');
  }
}

/**
 * Side-valve 'l-head' chamber (module comment). Construct with {@link createChamber}; `clone()` shares
 * the tables.
 */
export class LHeadChamber implements CombustionChamber {
  readonly kind = 'l-head' as const;
  readonly bore: number;
  readonly sparkCenter: readonly [number, number, number];
  readonly fixedVolume: number;
  readonly planformArea: number;
  readonly flame: FlameGeometry;
  private readonly t: LHeadTables;
  private readonly fr: FlameGeometryResult = newFlameGeometryResult();
  private readonly fr2: FlameGeometryResult = newFlameGeometryResult();
  private readonly lv = new Float64Array(4);
  private readonly areas = new Float64Array(N_WALL_SURFACES);
  /** Split points of the piston-side band quadrature. */
  private readonly ys: Float64Array;
  private cV = 0;
  private cA = 0;
  private lastRadius = 0;

  /** @internal use {@link createChamber} */
  constructor(tables: LHeadTables, flame?: FlameGeometry) {
    this.t = tables;
    this.flame = flame ?? tables.flame;
    this.bore = tables.bore;
    this.sparkCenter = tables.sparkCenter;
    this.fixedVolume = tables.fixedVolume;
    this.planformArea = tables.planformArea;
    this.ys = new Float64Array(2 * tables.openCritical.length + 2);
  }

  clone(): LHeadChamber {
    return new LHeadChamber(this.t, this.t.flame.clone());
  }

  /** Plan area of the valve pocket outside the bore A_Ω, m². */
  get pocketArea(): number {
    return this.t.region.area;
  }

  /** Table sizes (diagnostics): profile nodes, layer nodes. */
  get tableSizes(): { profile: number; layer: number } {
    return { profile: this.t.profile.size, layer: this.t.layer.size };
  }

  /** Volume and front area into cV, cA. */
  private core(r: number, h: number): void {
    const fr = this.flame.evaluateCore(r, h, this.fr);
    this.t.layer.lookup(r, this.lv);
    this.cV = fr.volume + this.lv[0];
    this.cA = fr.frontArea + this.lv[1];
  }

  volume(r: number, h: number): number {
    this.core(r, h);
    return this.cV;
  }

  frontArea(r: number, h: number): number {
    this.core(r, h);
    return this.cA;
  }

  chamberVolume(h: number): number {
    const R = this.flame.radius;
    return PI * R * R * h + this.fixedVolume;
  }

  meanDepth(h: number): number {
    return this.chamberVolume(h) / this.planformArea;
  }

  maxRadius(h: number): number {
    const a = this.flame.maxRadius(h);
    const b = this.t.layer.rMax;
    return a > b ? a : b;
  }

  inscribedRadius(h: number): number {
    if (this.t.sparkInPocket) return this.t.pocketInscribed;
    const fg = this.flame;
    const b = fg.headDistance;
    return h > b ? Math.min(b, h - b, fg.radius - fg.offset) : 0;
  }

  radiusForVolume(V: number, h: number, rGuess?: number): number {
    if (!(V > 0) || !(h > 0)) return 0;
    const rMax = this.maxRadius(h);
    if (V >= this.chamberVolume(h)) return rMax;
    const rIn = this.inscribedRadius(h);
    const rs = Math.cbrt((3 * V) / (4 * PI));
    if (rs <= rIn) return rs;
    let lo = rIn;
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
      let rn = this.cA > 0 ? r - f / this.cA : Number.NaN;
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

  surfaceAreas(h: number, out: Float64Array): Float64Array {
    const t = this.t;
    const R = t.R;
    const A = t.region.area;
    const disc = PI * R * R;
    const yLow = t.deckY > -h ? t.deckY : -h;
    const openH = t.roofY - yLow > 0 ? t.roofY - yLow : 0;
    const boreWall = TWO_PI * R * -yLow - t.region.openLength * openH;
    out[WALL_HEAD] = disc + A + (boreWall > 0 ? boreWall : 0) + t.region.sideLength * (t.roofY - t.deckY);
    const top = -h < t.roofY ? -h : t.roofY;
    out[WALL_PISTON] = disc + t.region.openLength * (top - t.deckY > 0 ? top - t.deckY : 0);
    out[WALL_LINER] = TWO_PI * R * (h + t.deckY > 0 ? h + t.deckY : 0);
    out[WALL_INTAKE_VALVE] = t.intake.area;
    out[WALL_EXHAUST_VALVE] = t.exhaust.area;
    out[WALL_BLOCK] = A - t.intake.area - t.exhaust.area;
    return out;
  }

  /**
   * ∫_{y0}^{y1} (open-arc length inside the slice circle) dy, m²: split where the slice radius crosses a
   * critical radius of the open arc (√-type kinks of its length), 4-point Gauss–Legendre in the cosine
   * variable y = a + (b − a)(1 − cos πs)/2 on each piece (smooth across √ end behaviour). Allocation-free.
   */
  private openBand(r: number, y0: number, y1: number): number {
    if (!(y1 > y0)) return 0;
    const t = this.t;
    const cy = t.sparkCenter[1];
    const ys = this.ys;
    const oc = t.openCritical;
    const r2 = r * r;
    let n = 0;
    ys[n++] = y0;
    for (let k = 0; k < oc.length; k++) {
      if (oc[k] >= r2) break;
      const dt = Math.sqrt(r2 - oc[k]);
      if (cy + dt > y0 && cy + dt < y1) ys[n++] = cy + dt;
      if (cy - dt > y0 && cy - dt < y1) ys[n++] = cy - dt;
    }
    ys[n++] = y1;
    // insertion sort (few points)
    for (let i = 1; i < n; i++) {
      const v = ys[i];
      let j = i - 1;
      while (j >= 0 && ys[j] > v) {
        ys[j + 1] = ys[j];
        j--;
      }
      ys[j + 1] = v;
    }
    let total = 0;
    const ocMax = oc[oc.length - 1];
    for (let p = 0; p + 1 < n; p++) {
      const a = ys[p];
      const b = ys[p + 1];
      if (!(b > a)) continue;
      const da = a - cy;
      const db = b - cy;
      const tm = da * da > db * db ? da * da : db * db;
      if (r2 - tm >= ocMax) {
        // the slice circle covers the whole open arc across this piece
        total += t.profile.open * (b - a);
        continue;
      }
      const half = 0.5 * (b - a);
      let s = 0;
      for (let k = 0; k < 4; k++) {
        const sv = 0.5 * (1 + GL4X[k]);
        const c = Math.cos(PI * sv);
        const y = a + half * (1 - c);
        const jac = half * PI * Math.sin(PI * sv);
        const dy = Math.abs(y - cy);
        if (r > dy) s += GL4W[k] * jac * t.profile.openLength((r - dy) * (r + dy));
      }
      total += 0.5 * s;
    }
    return total;
  }

  evaluate(r: number, h: number, out: ChamberResult): ChamberResult {
    const t = this.t;
    const R = t.R;
    const cy = t.sparkCenter[1];
    const fr = this.flame.evaluate(r, h, this.fr);
    const lv = this.lv;
    t.layer.lookup(r, lv);
    out.volume = fr.volume + lv[0];
    out.frontArea = fr.frontArea + lv[1];
    const a = this.surfaceAreas(h, this.areas);
    const w = out.wetted;
    // bore-circle wall above max(deckY, −h): FlameGeometry's liner integral of the column cut at the deck
    const hDeck = -t.deckY;
    let cwAbove = fr.wettedLiner;
    let liner = 0;
    if (h > hDeck) {
      cwAbove = this.flame.evaluateCore(r, hDeck, this.fr2).wettedLiner;
      liner = fr.wettedLiner - cwAbove;
    }
    // piston side facing the pocket while the crown is above the deck
    const yc = -h;
    const pSide = yc > t.deckY ? this.openBand(r, t.deckY, yc < t.roofY ? yc : t.roofY) : 0;
    const openAbove = lv[3] - pSide; // open arc between max(deckY, −h) and roofY
    const headWall = cwAbove - openAbove;
    // planes of the pocket
    const tr = Math.abs(t.tauRoof);
    const tf = Math.abs(t.tauFloor);
    const roof = r > tr ? t.profile.value((r - tr) * (r + tr)) : 0;
    let floor = 0;
    let iv = 0;
    let ev = 0;
    if (r > tf) {
      const u = (r - tf) * (r + tf);
      floor = t.profile.value(u);
      const rho = Math.sqrt(u);
      iv = t.intake.count * lensArea(rho, t.intake.radius, t.intake.dist);
      ev = t.exhaust.count * lensArea(rho, t.exhaust.radius, t.exhaust.dist);
    }
    w[WALL_HEAD] = fr.wettedHead + roof + (headWall > 0 ? headWall : 0) + lv[2];
    w[WALL_PISTON] = fr.wettedPiston + pSide;
    w[WALL_LINER] = liner > 0 ? liner : 0;
    w[WALL_INTAKE_VALVE] = iv;
    w[WALL_EXHAUST_VALVE] = ev;
    const blk = floor - iv - ev;
    w[WALL_BLOCK] = blk > 0 ? blk : 0;
    const f = out.burnedFraction;
    for (let i = 0; i < N_WALL_SURFACES; i++) f[i] = a[i] > 0 ? clamp01(w[i] / a[i]) : 0;
    // crevice mouth: share of the bore circle at the crown plane inside the ball
    const tc = Math.abs(h + cy);
    out.creviceBurnedFraction = r > tc ? arcOfBoreInside(Math.sqrt((r - tc) * (r + tc)), R, this.flame.offset) / TWO_PI : 0;
    return out;
  }
}

// ---------------------------------------------------------------------------------------------
// Factory and helpers
// ---------------------------------------------------------------------------------------------

/** Process-wide cache of built chambers (immutable tables), keyed by every input that shapes them. */
const CACHE = new Map<string, CombustionChamber>();
const CACHE_MAX = 16;

function cacheKey(spec: ChamberSpec, opts: ChamberOptions): string {
  const g = spec.geometry;
  const v = (x: ValveSpec) => [x.count, x.headDiameter, x.position[0], x.position[1]];
  return JSON.stringify([
    g.chamber,
    g.bore,
    g.stroke,
    g.lHead ?? null,
    spec.sparkPlug.gapCenter,
    v(spec.intakeValve),
    v(spec.exhaustValve),
    opts.maxHeight ?? null,
    opts.cellsInner ?? null,
    opts.cellsOuter ?? null,
    opts.cellsTau ?? null,
    opts.cornerBand ?? null,
  ]);
}

/**
 * The chamber of an engine spec, keyed on geometry.chamber: 'flat-disc' → {@link DiscChamber} over
 * `new FlameGeometry(bore, gapCenter, { maxHeight, … })` (the cycle model's legacy construction),
 * 'l-head' → {@link LHeadChamber}. With `cache` (default) a second call for the same inputs returns a
 * clone sharing the first one's tables (bit-identical numbers, no rebuild).
 */
export function createChamber(spec: ChamberSpec, opts: ChamberOptions = {}): CombustionChamber {
  const useCache = opts.cache ?? true;
  const key = useCache ? cacheKey(spec, opts) : '';
  if (useCache) {
    const hit = CACHE.get(key);
    if (hit) return hit.clone();
  }
  let ch: CombustionChamber;
  if (spec.geometry.chamber === 'l-head') {
    ch = new LHeadChamber(new LHeadTables(spec, opts));
  } else {
    const fgOpts: FlameGeometryOptions = {};
    if (opts.maxHeight !== undefined) fgOpts.maxHeight = opts.maxHeight;
    if (opts.cellsInner !== undefined) fgOpts.cellsInner = opts.cellsInner;
    if (opts.cellsOuter !== undefined) fgOpts.cellsOuter = opts.cellsOuter;
    if (opts.cellsTau !== undefined) fgOpts.cellsTau = opts.cellsTau;
    if (opts.cornerBand !== undefined) fgOpts.cornerBand = opts.cornerBand;
    const fg = new FlameGeometry(spec.geometry.bore, spec.sparkPlug.gapCenter, fgOpts);
    ch = new DiscChamber(fg, spec.intakeValve, spec.exhaustValve);
  }
  if (useCache) {
    if (CACHE.size >= CACHE_MAX) CACHE.delete(CACHE.keys().next().value as string);
    CACHE.set(key, ch);
    return ch.clone();
  }
  return ch;
}

/**
 * Chamber volume outside the bore column, m³: 0 for 'flat-disc'; 'l-head': pocket plan area × (roofY −
 * deckY), the same number as {@link LHeadChamber.fixedVolume}.
 */
export function chamberFixedVolume(g: EngineGeometrySpec): number {
  if (g.chamber !== 'l-head') return 0;
  if (!g.lHead) throw new Error("chamberFixedVolume: geometry.lHead is required for chamber 'l-head'");
  return lHeadPlanMetrics(g.bore, g.lHead).fixedVolume;
}

/**
 * SliderCrank geometry of a spec with `fixedChamberVolume` filled in from the chamber shape, so that
 * clearanceHeight(θ) is the depth of the bore column (L-head: h_TDC = −deckY − crownAboveDeckAtTDC when
 * the spec's compression ratio closes the clearance volume). Flat disc: fixedChamberVolume 0 (identical
 * to passing the spec geometry).
 */
export function sliderCrankGeometry(g: EngineGeometrySpec): SliderCrankGeometry {
  return {
    bore: g.bore,
    stroke: g.stroke,
    conRodLength: g.conRodLength,
    pinOffset: g.pinOffset,
    creviceVolume: g.creviceVolume,
    fixedChamberVolume: chamberFixedVolume(g),
  };
}
