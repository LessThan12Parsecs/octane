/**
 * Exact slider-crank kinematics with wrist-pin offset.
 *
 * ── Frames and signs ────────────────────────────────────────────────────────
 * Mechanism plane = world x-y plane (crank axis = world z through the origin,
 * cylinder axis = world +y, pointing from the crank toward the head).
 *
 *  - φ ("axis angle") is the crank-pin angle measured from the +y axis toward +x
 *    in the direction of rotation: crank pin at (a sinφ, a cosφ). Rotation is
 *    therefore CLOCKWISE when viewed from +z (angular-velocity vector along −z).
 *  - θ is the crank angle of the physics convention: θ = 0 at the exact firing
 *    TDC (DESIGN.md), θ = φ − φ_TDC. The mechanism is 2π-periodic; the 4-stroke
 *    cycle spans 4π, so θ and θ ± 2π give the same geometry.
 *  - The wrist pin moves on the line x = x_w = pinOffset (same convention as
 *    src/render/engine/kinematics.ts). The MAJOR-thrust wall is −x (during the
 *    expansion stroke the crank pin is at +x and the gas-loaded, compressed rod
 *    presses the piston against the −x wall), so the usual anti-slap offset
 *    toward the major-thrust side is a NEGATIVE pinOffset here.
 *  - Rod angle β: angle of the rod (big end → small end) from the cylinder axis,
 *    sinβ = (a sinφ − x_w)/l, positive when the crank pin is on the +x side of
 *    the wrist-pin line; β is the rod's rotation about world +z (CCW positive),
 *    i.e. the direction vector big end → small end is (−sinβ, cosβ).
 *  - x(θ) = piston (crown) displacement below its TDC position, ≥ 0.
 *
 * ── Volume / crevice convention ─────────────────────────────────────────────
 * The geometric compression ratio is CR = (V_c + V_d)/V_c where V_c is the TOTAL
 * clearance volume at TDC — the quantity measured by oil filling (e.g. Choi et
 * al. 2018, SAE 2018-01-0848) — and therefore INCLUDES the top-land/ring
 * crevice, plug and pickup cavities (lumped in `creviceVolume`).
 *   volume(θ)          = V_c + A_p·x(θ)             total cylinder gas volume
 *   clearanceHeight(θ) = h_TDC + x(θ),  h_TDC = (V_c − V_crevice − V_fixed)/A_p
 * so volume(θ) = A_p·clearanceHeight(θ) + V_crevice + V_fixed exactly: the bore
 * column of height h (the flat-disc "pancake" chamber when V_fixed = 0, the
 * default) plus a constant crevice volume and a constant non-bore chamber volume
 * V_fixed = `fixedChamberVolume` (the valve pocket of a side-valve 'l-head'
 * chamber, combustion/chamber.ts chamberFixedVolume; then h is the depth of the
 * bore column below the head-cavity roof). Piston crown and the head face over
 * the bore are flat discs of bore diameter.
 *
 * Validation (kinematics.test.ts): numerical loop-closure oracle
 * (tools/reference/mechanics_kinematics.py), Heywood (1988) eqs. 2.4–2.6 for
 * zero offset, finite differences of every analytic derivative.
 *
 * Hot-path methods allocate nothing.
 */
import type { EngineGeometrySpec } from '../core/engine-spec';

/** The subset of EngineGeometrySpec the kinematics needs (all SI, m / m³). */
export type SliderCrankGeometry = Pick<
  EngineGeometrySpec,
  'bore' | 'stroke' | 'conRodLength' | 'pinOffset' | 'creviceVolume'
> & {
  /**
   * Chamber volume outside the bore column (side-valve pocket), part of the clearance volume, m³.
   * Absent / 0: the flat-disc convention (h_TDC = (V_c − V_crevice)/A_p). See the file header.
   */
  fixedChamberVolume?: number;
};

/** All kinematic quantities at one crank angle (SI; reused output record, see SliderCrank.evaluate). */
export interface KinematicState {
  /** Crank angle θ, rad (firing-TDC convention). */
  theta: number;
  /** sin φ, cos φ of the crank-pin axis angle φ = θ + φ_TDC. */
  sinPhi: number;
  cosPhi: number;
  /** Piston displacement below TDC x, m; dx/dθ, m/rad; d²x/dθ², m/rad². */
  x: number;
  dxdTheta: number;
  d2xdTheta2: number;
  /** Clearance height h (flat disc: head-to-crown height; L-head: bore-column depth below the cavity roof), m. */
  clearanceHeight: number;
  /** Cylinder volume V, m³; dV/dθ, m³/rad. */
  volume: number;
  dVdTheta: number;
  /** Rod angle β, rad; sin β, cos β; dβ/dθ; d²β/dθ². */
  beta: number;
  sinBeta: number;
  cosBeta: number;
  dBetadTheta: number;
  d2BetadTheta2: number;
}

/** Allocate a zeroed KinematicState (once, outside the hot loop). */
export const newKinematicState = (): KinematicState => ({
  theta: 0,
  sinPhi: 0,
  cosPhi: 1,
  x: 0,
  dxdTheta: 0,
  d2xdTheta2: 0,
  clearanceHeight: 0,
  volume: 0,
  dVdTheta: 0,
  beta: 0,
  sinBeta: 0,
  cosBeta: 1,
  dBetadTheta: 0,
  d2BetadTheta2: 0,
});

export class SliderCrank {
  /** Cylinder bore, m. */
  readonly bore: number;
  /** Crank throw a = stroke/2 (EngineSpec.stroke is defined as 2 × crank radius), m. */
  readonly crankRadius: number;
  /** Connecting-rod length l (centre to centre), m. */
  readonly rodLength: number;
  /** Wrist-pin offset = x-coordinate of the wrist-pin line (major-thrust wall is −x, see file header), m. */
  readonly pinOffset: number;
  /** x-coordinate of the wrist-pin line, x_w = pinOffset, m. */
  readonly pinLineX: number;
  /** Lumped crevice volume, included in the clearance volume, m³. */
  readonly creviceVolume: number;
  /** Chamber volume outside the bore column (L-head valve pocket), included in the clearance volume, m³ (0: flat disc). */
  readonly fixedChamberVolume: number;
  /** Bore cross-section πB²/4 (= flat piston-crown area = flat head-face area), m². */
  readonly boreArea: number;
  /** Flat head fire-deck area exposed to the gas (πB²/4, valve faces included), m². */
  readonly headArea: number;
  /** Flat piston-crown area (πB²/4), m². */
  readonly pistonCrownArea: number;
  /** Crank-pin axis angle φ at the exact TDC, rad (0 for zero offset). */
  readonly tdcAxisAngle: number;
  /** Crank-pin axis angle φ at the exact BDC, rad (π for zero offset). */
  readonly bdcAxisAngle: number;
  /** Crank angle θ (firing-TDC convention) of the BDC following θ = 0, rad (π for zero offset). */
  readonly bdcAngle: number;
  /** Exact TDC→BDC piston travel √((l+a)²−x_w²) − √((l−a)²−x_w²), m (= 2a without offset). */
  readonly pistonTravel: number;
  /** Swept (displaced) volume A_p × pistonTravel, m³. */
  readonly displacedVolume: number;
  /** Wrist-pin distance from the crank axis along +y at TDC, m. */
  readonly pinHeightTDC: number;

  private cr = 0;
  private vc = 0;
  private hTdc = 0;

  /**
   * @param geometry bore, stroke (= 2 × crank radius), rod length, pin offset, crevice volume (SI)
   * @param compressionRatio geometric compression ratio (V_c + V_d)/V_c, > 1
   */
  constructor(geometry: SliderCrankGeometry, compressionRatio: number) {
    const { bore, stroke, conRodLength, pinOffset, creviceVolume } = geometry;
    if (!(bore > 0 && stroke > 0 && conRodLength > 0)) throw new RangeError('SliderCrank: bore, stroke and rod length must be > 0');
    const a = stroke / 2;
    const l = conRodLength;
    const xw = pinOffset;
    if (!(l - a > Math.abs(xw))) throw new RangeError('SliderCrank: mechanism cannot pass BDC (l − a ≤ |offset|)');
    if (!(creviceVolume >= 0)) throw new RangeError('SliderCrank: creviceVolume must be ≥ 0');
    const fixedChamberVolume = geometry.fixedChamberVolume ?? 0;
    if (!(fixedChamberVolume >= 0)) throw new RangeError('SliderCrank: fixedChamberVolume must be ≥ 0');
    this.bore = bore;
    this.crankRadius = a;
    this.rodLength = l;
    this.pinOffset = pinOffset;
    this.pinLineX = xw;
    this.creviceVolume = creviceVolume;
    this.fixedChamberVolume = fixedChamberVolume;
    this.boreArea = (Math.PI * bore * bore) / 4;
    this.headArea = this.boreArea;
    this.pistonCrownArea = this.boreArea;
    // TDC: crank and rod collinear and extended → sinφ = x_w/(l + a).
    // BDC: collinear and folded (crank pin below the axis) → φ = π + asin(x_w/(l − a)).
    this.tdcAxisAngle = Math.asin(xw / (l + a));
    this.bdcAxisAngle = Math.PI + Math.asin(xw / (l - a));
    this.bdcAngle = this.bdcAxisAngle - this.tdcAxisAngle;
    this.pinHeightTDC = Math.sqrt((l + a) * (l + a) - xw * xw);
    this.pistonTravel = this.pinHeightTDC - Math.sqrt((l - a) * (l - a) - xw * xw);
    this.displacedVolume = this.boreArea * this.pistonTravel;
    this.setCompressionRatio(compressionRatio);
  }

  /** Geometric compression ratio (V_c + V_d)/V_c. */
  get compressionRatio(): number {
    return this.cr;
  }

  /** Total clearance volume at TDC (crevice included) V_d/(CR − 1), m³. */
  get clearanceVolume(): number {
    return this.vc;
  }

  /** Clearance height at TDC, (V_c − V_crevice − V_fixed)/A_p, m. */
  get clearanceHeightTDC(): number {
    return this.hTdc;
  }

  /**
   * Change the geometric compression ratio (variable-CR engines: the CFR raises
   * the cylinder+head; the mechanism below the crown is unchanged).
   * @throws RangeError if CR ≤ 1 or the crevice volume would exceed the clearance volume.
   */
  setCompressionRatio(compressionRatio: number): void {
    if (!(compressionRatio > 1)) throw new RangeError('SliderCrank: compression ratio must be > 1');
    const vc = this.displacedVolume / (compressionRatio - 1);
    const h = (vc - this.creviceVolume - this.fixedChamberVolume) / this.boreArea;
    if (!(h > 0)) throw new RangeError('SliderCrank: crevice (+ fixed chamber) volume ≥ clearance volume at this compression ratio');
    this.cr = compressionRatio;
    this.vc = vc;
    this.hTdc = h;
  }

  /** Clearance volume (crevice included) that gives compression ratio CR (> 1) with this mechanism, m³ (NaN if CR ≤ 1). */
  clearanceVolumeForCR(compressionRatio: number): number {
    return compressionRatio > 1 ? this.displacedVolume / (compressionRatio - 1) : Number.NaN;
  }

  /** Clearance height at TDC for compression ratio CR (> 1), m (see volume convention; NaN if CR ≤ 1). */
  clearanceHeightTDCForCR(compressionRatio: number): number {
    return (this.clearanceVolumeForCR(compressionRatio) - this.creviceVolume - this.fixedChamberVolume) / this.boreArea;
  }

  /** Wrist-pin distance from the crank axis along the cylinder axis at crank angle θ (rad), m. */
  pinHeight(theta: number): number {
    const phi = theta + this.tdcAxisAngle;
    const s = this.crankRadius * Math.sin(phi) - this.pinLineX;
    return this.crankRadius * Math.cos(phi) + Math.sqrt(this.rodLength * this.rodLength - s * s);
  }

  /** Piston (crown) displacement below its TDC position at crank angle θ (rad), m (0 at TDC). */
  pistonDisplacement(theta: number): number {
    return this.pinHeightTDC - this.pinHeight(theta);
  }

  /** Instantaneous clearance height (bore-column depth) h(θ) = h_TDC + x(θ), m. */
  clearanceHeight(theta: number): number {
    return this.hTdc + this.pistonDisplacement(theta);
  }

  /** Total cylinder gas volume V(θ) = V_c + A_p x(θ) (crevice included), m³. */
  volume(theta: number): number {
    return this.vc + this.boreArea * this.pistonDisplacement(theta);
  }

  /** dx/dθ, m/rad (positive while the piston moves away from the head). */
  dxdTheta(theta: number): number {
    const a = this.crankRadius;
    const phi = theta + this.tdcAxisAngle;
    const sinp = Math.sin(phi);
    const cosp = Math.cos(phi);
    const s = a * sinp - this.pinLineX;
    const c = Math.sqrt(this.rodLength * this.rodLength - s * s);
    return a * sinp + (s * a * cosp) / c;
  }

  /** d²x/dθ², m/rad². */
  d2xdTheta2(theta: number): number {
    const a = this.crankRadius;
    const phi = theta + this.tdcAxisAngle;
    const sinp = Math.sin(phi);
    const cosp = Math.cos(phi);
    const s = a * sinp - this.pinLineX;
    const c2 = this.rodLength * this.rodLength - s * s;
    const c = Math.sqrt(c2);
    const ds = a * cosp; // ds/dφ
    // x' = a sinφ + s s'/c ;  x'' = a cosφ + (s'² + s s'')/c + s² s'²/c³, s'' = −a sinφ
    return a * cosp + (ds * ds - s * a * sinp) / c + (s * s * ds * ds) / (c2 * c);
  }

  /** Analytic dV/dθ = A_p dx/dθ, m³/rad. */
  dVdTheta(theta: number): number {
    return this.boreArea * this.dxdTheta(theta);
  }

  /** Analytic d²V/dθ², m³/rad². */
  d2VdTheta2(theta: number): number {
    return this.boreArea * this.d2xdTheta2(theta);
  }

  /** Rod angle β from the cylinder axis, rad (sign: see file header). */
  rodAngle(theta: number): number {
    const phi = theta + this.tdcAxisAngle;
    return Math.asin((this.crankRadius * Math.sin(phi) - this.pinLineX) / this.rodLength);
  }

  /** dβ/dθ, rad/rad. */
  dRodAngledTheta(theta: number): number {
    const a = this.crankRadius;
    const phi = theta + this.tdcAxisAngle;
    const s = a * Math.sin(phi) - this.pinLineX;
    return (a * Math.cos(phi)) / Math.sqrt(this.rodLength * this.rodLength - s * s);
  }

  /** d²β/dθ², rad/rad². */
  d2RodAngledTheta2(theta: number): number {
    const a = this.crankRadius;
    const phi = theta + this.tdcAxisAngle;
    const sinp = Math.sin(phi);
    const cosp = Math.cos(phi);
    const s = a * sinp - this.pinLineX;
    const c2 = this.rodLength * this.rodLength - s * s;
    const c = Math.sqrt(c2);
    // β' = s'/c ; β'' = s''/c + s s'²/c³
    return (-a * sinp) / c + (s * a * a * cosp * cosp) / (c2 * c);
  }

  /** Rod angular velocity dβ/dt at crank speed ω (rad/s), rad/s. */
  rodAngularVelocity(theta: number, omega: number): number {
    return this.dRodAngledTheta(theta) * omega;
  }

  /** Rod angular acceleration d²β/dt² at crank speed ω and acceleration α (rad/s²), rad/s². */
  rodAngularAcceleration(theta: number, omega: number, alpha: number): number {
    return this.d2RodAngledTheta2(theta) * omega * omega + this.dRodAngledTheta(theta) * alpha;
  }

  /** Piston velocity dx/dt, m/s (positive away from the head), for crank speed ω, rad/s. */
  pistonVelocity(theta: number, omega: number): number {
    return this.dxdTheta(theta) * omega;
  }

  /** Piston acceleration d²x/dt², m/s² (positive away from the head), for ω (rad/s) and α (rad/s²). */
  pistonAcceleration(theta: number, omega: number, alpha: number): number {
    return this.d2xdTheta2(theta) * omega * omega + this.dxdTheta(theta) * alpha;
  }

  /** Mean piston speed 2·travel·N/60, m/s, at rpm (Heywood 1988 eq. 2.7 with exact travel). */
  meanPistonSpeed(rpm: number): number {
    return (2 * this.pistonTravel * rpm) / 60;
  }

  /** Exposed liner (cylinder-wall) area between head face and crown, πB·h(θ), m² (flat disc only). */
  linerArea(theta: number): number {
    return Math.PI * this.bore * this.clearanceHeight(theta);
  }

  /** Total disc-chamber surface area head + crown + liner (Heywood 1988 eq. 2.8), m² (flat disc only; see combustion/chamber.ts surfaceAreas). */
  chamberSurfaceArea(theta: number): number {
    return this.headArea + this.pistonCrownArea + this.linerArea(theta);
  }

  /**
   * Everything at once (one sin/cos/sqrt): the hot-path entry point for the cycle
   * integrator and the dynamics. Writes into and returns `out`.
   */
  evaluate(theta: number, out: KinematicState): KinematicState {
    const a = this.crankRadius;
    const l = this.rodLength;
    const phi = theta + this.tdcAxisAngle;
    const sinp = Math.sin(phi);
    const cosp = Math.cos(phi);
    const s = a * sinp - this.pinLineX;
    const c2 = l * l - s * s;
    const c = Math.sqrt(c2);
    const ds = a * cosp;
    const x = this.pinHeightTDC - (a * cosp + c);
    const x1 = a * sinp + (s * ds) / c;
    const x2 = a * cosp + (ds * ds - s * a * sinp) / c + (s * s * ds * ds) / (c2 * c);
    out.theta = theta;
    out.sinPhi = sinp;
    out.cosPhi = cosp;
    out.x = x;
    out.dxdTheta = x1;
    out.d2xdTheta2 = x2;
    out.clearanceHeight = this.hTdc + x;
    out.volume = this.vc + this.boreArea * x;
    out.dVdTheta = this.boreArea * x1;
    out.beta = Math.asin(s / l);
    out.sinBeta = s / l;
    out.cosBeta = c / l;
    out.dBetadTheta = ds / c;
    out.d2BetadTheta2 = (-a * sinp) / c + (s * ds * ds) / (c2 * c);
    return out;
  }

  /**
   * Crank-pin centre in the mechanism plane (world x, y; origin on the crank axis), m.
   * Writes into `out[0..1]` and returns it.
   */
  crankPinPosition(theta: number, out: Float64Array | number[]): Float64Array | number[] {
    const phi = theta + this.tdcAxisAngle;
    out[0] = this.crankRadius * Math.sin(phi);
    out[1] = this.crankRadius * Math.cos(phi);
    return out;
  }

  /** Wrist-pin centre (world x, y), m. Writes into `out[0..1]` and returns it. */
  wristPinPosition(theta: number, out: Float64Array | number[]): Float64Array | number[] {
    out[0] = this.pinLineX;
    out[1] = this.pinHeight(theta);
    return out;
  }
}
