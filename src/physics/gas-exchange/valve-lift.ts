/**
 * Cam-derived poppet-valve lift profile L(θ) in the firing-TDC crank-angle convention.
 *
 * ── Profile ─────────────────────────────────────────────────────────────────────────────
 * One symmetric lift event per 720° cycle, written as a single even "polydyne"-type
 * polynomial of the normalised cam angle x = (θ − θ_c)/H ∈ [−1, 1] (θ_c = event centre =
 * nose of the cam, H = cam half duration from the base circle, crank deg):
 *
 *   L_cam(θ) = (L_max + lash) · s(x),   s(x) = 1 + Σ_i c_i x^{p_i}   (even powers p_i),
 *   L_cam = 0 for |x| ≥ 1,
 *
 * with the n coefficients c_i fixed by s(1) = 0 and s^{(k)}(1) = 0 for k = 1 … n−1. With the
 * default powers (2, 10, 18, 26) the cam rise is C³ everywhere (velocity, acceleration and
 * jerk all vanish at the base circle; the valve lift is C³ at the seat only when lash = 0) and:
 *  - the nose acceleration is negative, s''(0) = 2c₁ = −3.047 (spring-controlled deceleration);
 *  - the peak positive (flank) acceleration is 9.98 = 3.3 × |nose| (typical of valve cams);
 *  - fullness ∫s dx / 2 = 0.553; s is monotone on each flank; the lift near the seat grows
 *    as ∝ (1 − |x|)⁴, i.e. a built-in gentle opening/closing ramp.
 * This is the polynomial "polydyne" family of Dudley (1948, Trans. ASME 70) and Stoddart
 * (1953, Machine Design 25, "Polydyne cam design"), as presented in Norton (2009), "Cam
 * Design and Manufacturing Handbook", ch. on polynomial cam functions.
 * UNVERIFIED: the specific power set 2-10-18-26 as the Stoddart example (from memory); the
 * shape metrics above are computed here, not quoted.
 *
 * ── Valve lash (running clearance) ──────────────────────────────────────────────────────
 * The polynomial describes the CAM rise transmitted to the valve at zero clearance,
 * L_cam = (maxLift + lash)·s(x). A mechanical (solid-lifter) valvetrain runs with a hot
 * clearance that the cam takes up before the valve moves (CFR: 0.008 in, ASTM D2699 §10.3.2.1,
 * taken up on the lobe's 0.008–0.010 in quieting ramps, D2699 A2.1), so the VALVE lift is
 *   L(θ) = max(0, L_cam(θ) − lash),
 * and the valve leaves / lands on its seat with the finite cam velocity at L_cam = lash — not
 * with the zero velocity, acceleration and jerk of the cam's base-circle departure. With
 * lash = 0 the seat events coincide with the cam's base-circle departure and the valve creeps
 * off the seat as (1 − |x|)⁴: against the measured CFR profiles of Choi et al. (2018, SAE
 * 2018-01-0848, Fig. 3; engines/cfr-valve-lift.ts) that under-predicts the lift 10° after IVO /
 * before IVC by 6.5× / 4.7×, the lift within 20° of the valve events by 70 % and the lift
 * integral by 11 % (intake) / 15 % (exhaust). With the 0.008 in lash the same cam reproduces the
 * measured intake profile to 0.095 mm rms (1.6 % of the peak), its lift integral to 0.7 % and
 * the lift near the events to 4 % (exhaust: 0.27 mm rms, −5 %, +25 %) (valve-lift.test.ts).
 * `lash` is an optional field of the timing spec (default 0; ValveSpec.lash).
 *
 * ── Timing ──────────────────────────────────────────────────────────────────────────────
 * The CFR (and most engines) quote timing at a small lift threshold of the VALVE (or at the
 * seat with the running clearance): s(x_th) = (timingLiftThreshold + lash)/(maxLift + lash) is
 * solved once, and the cam half duration is H = (quoted duration/2)/x_th, so
 * L(openDeg) = L(closeDeg) = timingLiftThreshold exactly. With timingLiftThreshold = 0 the
 * quoted angles are the valve seat-off / seat-on angles.
 *
 * ── Wrap ────────────────────────────────────────────────────────────────────────────────
 * Angles are taken modulo 720°; the quoted duration is (closeDeg − openDeg) mod 720 ∈ (0, 720),
 * so an event may straddle the ±360° seam (e.g. exhaust 140° → −345° ≡ 375°).
 *
 * ── Lift profiles from the spec ─────────────────────────────────────────────────────────
 * {@link LiftProfile} is the interface the cycle model and the render consume. Implementations:
 * this polydyne {@link ValveLiftProfile} (ValveSpec.cam absent or 'polydyne'), the exact three-arc
 * lobe on a flat follower and the measured table (cam-lift.ts). {@link createLiftProfile} builds the
 * one the ValveSpec asks for, with the running clearance resolved by {@link resolveValveLash}.
 *
 * Units: lift m, angles crank degrees, lift rate m/deg, lift acceleration m/deg².
 * Hot paths allocate nothing.
 */
import type { CamSpec, ValveSpec } from '../core/engine-spec';
import { TabulatedLiftProfile, ThreeArcFlatFollowerProfile } from './cam-lift';

/**
 * Valve lift L(θ) of one valve over the 720° cycle (θ = the cylinder's local crank angle, deg,
 * firing-TDC convention), all methods allocation-free.
 */
export interface LiftProfile {
  /** Maximum VALVE lift (after the running clearance), m. */
  readonly maxLift: number;
  /** Running clearance (lash) between follower and valve, m. */
  readonly lash: number;
  /** Maximum cam rise at the valve at zero clearance (maxLift + lash), m. */
  readonly camMaxLift: number;
  /** Crank angle of maximum lift (cam nose), deg in [−360, 360). */
  readonly centerDeg: number;
  /** Valve seat-off angle, deg in [−360, 360). */
  readonly seatOpenDeg: number;
  /** Valve seat-on angle, deg in [−360, 360). */
  readonly seatCloseDeg: number;
  /** Valve lift at crank angle θ (deg), m (0 when seated). */
  lift(thetaDeg: number): number;
  /** Cam rise transmitted to the valve at zero clearance (tappet motion; render), m. */
  camLift(thetaDeg: number): number;
  /** dL/dθ of the valve, m per crank degree (0 while seated). */
  liftRate(thetaDeg: number): number;
  /** d²L/dθ² of the valve, m per crank degree² (0 while seated). */
  liftAccel(thetaDeg: number): number;
  /** True while the valve is off its seat. */
  isOpen(thetaDeg: number): boolean;
}

/** Even-power lift polynomial s(x) = 1 + Σ coeffs[i]·x^powers[i] on |x| ≤ 1. */
export interface CamPolynomial {
  /** Even, strictly increasing exponents (first = 2). */
  readonly powers: readonly number[];
  /** Coefficients (same length as powers). */
  readonly coeffs: readonly number[];
}

/**
 * Coefficients of the even polynomial s(x) = 1 + Σ c_i x^{p_i} with s(1) = 0 and
 * s^{(k)}(1) = 0 for k = 1 … n−1 (n = powers.length): a C^{n−1} closure onto zero lift.
 * Solves the n×n linear system by Gaussian elimination with partial pivoting.
 */
export function polydyneCoefficients(powers: readonly number[]): number[] {
  const n = powers.length;
  if (n < 2) throw new RangeError('polydyneCoefficients: need ≥ 2 powers');
  for (let i = 0; i < n; i++) {
    const p = powers[i];
    if (!(p >= 2 && p % 2 === 0) || (i > 0 && !(p > powers[i - 1]))) {
      throw new RangeError('polydyneCoefficients: powers must be even, ≥ 2, strictly increasing');
    }
  }
  const A: number[][] = [];
  const b: number[] = [];
  for (let k = 0; k < n; k++) {
    const row: number[] = [];
    for (let i = 0; i < n; i++) {
      // k-th derivative of x^p at x = 1: p (p−1) … (p−k+1)
      let f = 1;
      for (let j = 0; j < k; j++) f *= powers[i] - j;
      row.push(f);
    }
    A.push(row);
    b.push(k === 0 ? -1 : 0);
  }
  for (let c = 0; c < n; c++) {
    let piv = c;
    for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[piv][c])) piv = r;
    [A[c], A[piv]] = [A[piv], A[c]];
    [b[c], b[piv]] = [b[piv], b[c]];
    for (let r = c + 1; r < n; r++) {
      const f = A[r][c] / A[c][c];
      for (let j = c; j < n; j++) A[r][j] -= f * A[c][j];
      b[r] -= f * b[c];
    }
  }
  const x = new Array<number>(n).fill(0);
  for (let r = n - 1; r >= 0; r--) {
    let s = b[r];
    for (let j = r + 1; j < n; j++) s -= A[r][j] * x[j];
    x[r] = s / A[r][r];
  }
  return x;
}

/** Build a {@link CamPolynomial} (C^{n−1} at the seat) from even powers. */
export function polydyneCam(powers: readonly number[]): CamPolynomial {
  return Object.freeze({ powers: Object.freeze([...powers]), coeffs: Object.freeze(polydyneCoefficients(powers)) });
}

/** Default valve cam: 2-10-18-26 polydyne (C³ at the seat; see file header). */
export const POLYDYNE_2_10_18_26: CamPolynomial = polydyneCam([2, 10, 18, 26]);

/** Wrap crank degrees into [−360, 360). */
export function wrapDeg720(deg: number): number {
  let d = (deg + 360) % 720;
  if (d < 0) d += 720;
  return d - 360;
}

/**
 * The ValveSpec fields the lift profile needs, plus the optional valve running clearance
 * `lash` (m, default 0; ValveSpec.lash — see file header).
 */
export type ValveTimingSpec = Pick<ValveSpec, 'maxLift' | 'openDeg' | 'closeDeg' | 'timingLiftThreshold'> & {
  /** Valve (tappet) running clearance taken up on the cam before the valve moves, m (default 0). */
  readonly lash?: number;
};

/**
 * Precomputed lift profile of one valve (recommended for the hot path: no per-call setup).
 * All methods are allocation-free.
 */
export class ValveLiftProfile implements LiftProfile {
  /** Maximum VALVE lift L_max, m. */
  readonly maxLift: number;
  /** Valve running clearance (lash), m. */
  readonly lash: number;
  /** Maximum cam rise at the valve at zero clearance, maxLift + lash, m. */
  readonly camMaxLift: number;
  /** Event centre (cam nose) θ_c, crank deg in [−360, 360). */
  readonly centerDeg: number;
  /** Cam half duration H (base circle to nose), crank deg. */
  readonly halfDurationDeg: number;
  /** Valve seat-off angle θ_c − H·x_seat, wrapped to [−360, 360). */
  readonly seatOpenDeg: number;
  /** Valve seat-on angle θ_c + H·x_seat, wrapped to [−360, 360). */
  readonly seatCloseDeg: number;
  /** Normalised cam position x_th at which the quoted timing threshold lift is reached. */
  readonly thresholdX: number;
  /** Normalised cam position x_seat at which the lash is taken up (1 when lash = 0). */
  readonly seatX: number;
  private readonly pw: Float64Array;
  private readonly cf: Float64Array;
  private readonly n: number;

  /**
   * @param spec max valve lift (m), open/close angles (crank deg) quoted at the valve lift
   *   timingLiftThreshold (m), optional running clearance `lash` (m)
   * @param cam lift polynomial (default 2-10-18-26 polydyne)
   */
  constructor(spec: ValveTimingSpec, cam: CamPolynomial = POLYDYNE_2_10_18_26) {
    const { maxLift, openDeg, closeDeg } = spec;
    const thr = spec.timingLiftThreshold > 0 ? spec.timingLiftThreshold : 0;
    const lash = spec.lash !== undefined && spec.lash > 0 ? spec.lash : 0;
    if (!(maxLift > 0)) throw new RangeError('ValveLiftProfile: maxLift must be > 0');
    if (!(thr < maxLift)) throw new RangeError('ValveLiftProfile: timingLiftThreshold must be < maxLift');
    if (!Number.isFinite(lash)) throw new RangeError('ValveLiftProfile: lash must be finite');
    const dur = (((closeDeg - openDeg) % 720) + 720) % 720;
    if (!(dur > 0)) throw new RangeError('ValveLiftProfile: zero event duration');
    this.maxLift = maxLift;
    this.lash = lash;
    const camMax = maxLift + lash;
    this.camMaxLift = camMax;
    this.n = cam.powers.length;
    this.pw = Float64Array.from(cam.powers);
    this.cf = Float64Array.from(cam.coeffs);
    const xth = this.solveShape((thr + lash) / camMax);
    this.thresholdX = xth;
    this.seatX = thr > 0 ? this.solveShape(lash / camMax) : xth;
    const half = dur / 2 / xth;
    if (!(half < 360)) throw new RangeError('ValveLiftProfile: cam duration exceeds 720°');
    this.halfDurationDeg = half;
    this.centerDeg = wrapDeg720(openDeg + dur / 2);
    this.seatOpenDeg = wrapDeg720(this.centerDeg - half * this.seatX);
    this.seatCloseDeg = wrapDeg720(this.centerDeg + half * this.seatX);
  }

  /**
   * Root x ∈ (0, 1] of s(x) = target (0 ≤ target < 1); s is decreasing on [0, 1] →
   * safeguarded Newton. target ≤ 0 → 1 (the base-circle departure).
   */
  private solveShape(target: number): number {
    if (!(target > 0)) return 1;
    let lo = 0;
    let hi = 1;
    let x = 0.9;
    for (let it = 0; it < 200; it++) {
      const f = this.shape(x) - target;
      if (f > 0) lo = x;
      else hi = x;
      const df = this.shapeD1(x);
      let xn = df < 0 ? x - f / df : 0.5 * (lo + hi);
      if (!(xn > lo && xn < hi)) xn = 0.5 * (lo + hi);
      if (Math.abs(xn - x) < 1e-15) return xn;
      x = xn;
    }
    return x;
  }

  /** s(x) for |x| ≤ 1 (dimensionless lift fraction). */
  shape(x: number): number {
    const u = x * x;
    let s = 1;
    let prevHalf = 0;
    let up = 1;
    for (let i = 0; i < this.n; i++) {
      const h = this.pw[i] / 2;
      for (let j = prevHalf; j < h; j++) up *= u;
      prevHalf = h;
      s += this.cf[i] * up;
    }
    return s;
  }

  /** ds/dx. */
  shapeD1(x: number): number {
    // Σ c p x^{p−1} = x Σ c p u^{p/2 − 1}
    const u = x * x;
    let s = 0;
    let prevHalf = 1;
    let up = 1;
    for (let i = 0; i < this.n; i++) {
      const h = this.pw[i] / 2;
      for (let j = prevHalf; j < h; j++) up *= u;
      prevHalf = h;
      s += this.cf[i] * this.pw[i] * up;
    }
    return x * s;
  }

  /** d²s/dx². */
  shapeD2(x: number): number {
    // Σ c p (p−1) u^{p/2 − 1}
    const u = x * x;
    let s = 0;
    let prevHalf = 1;
    let up = 1;
    for (let i = 0; i < this.n; i++) {
      const h = this.pw[i] / 2;
      for (let j = prevHalf; j < h; j++) up *= u;
      prevHalf = h;
      s += this.cf[i] * this.pw[i] * (this.pw[i] - 1) * up;
    }
    return s;
  }

  /** Normalised cam position x = (θ − θ_c)/H with the 720° wrap (|x| ≥ 1: valve seated). */
  camX(thetaDeg: number): number {
    return wrapDeg720(thetaDeg - this.centerDeg) / this.halfDurationDeg;
  }

  /** Valve lift at crank angle θ (deg), m (0 when seated): max(0, L_cam − lash). */
  lift(thetaDeg: number): number {
    const x = this.camX(thetaDeg);
    const xs = this.seatX;
    if (!(x > -xs && x < xs)) return 0;
    const l = this.camMaxLift * this.shape(x) - this.lash;
    return l > 0 ? l : 0;
  }

  /**
   * Cam rise transmitted to the valve at zero clearance, L_cam (m): the tappet/rocker motion
   * (for rendering the valvetrain). Equals lift + lash while the valve is open.
   */
  camLift(thetaDeg: number): number {
    const x = this.camX(thetaDeg);
    if (!(x > -1 && x < 1)) return 0;
    const s = this.shape(x);
    return s > 0 ? this.camMaxLift * s : 0;
  }

  /** dL/dθ of the valve, m per crank degree (0 while seated; jumps at the seat when lash > 0). */
  liftRate(thetaDeg: number): number {
    const x = this.camX(thetaDeg);
    const xs = this.seatX;
    if (!(x > -xs && x < xs)) return 0;
    return (this.camMaxLift * this.shapeD1(x)) / this.halfDurationDeg;
  }

  /** d²L/dθ² of the valve, m per crank degree² (0 while seated). */
  liftAccel(thetaDeg: number): number {
    const x = this.camX(thetaDeg);
    const xs = this.seatX;
    if (!(x > -xs && x < xs)) return 0;
    return (this.camMaxLift * this.shapeD2(x)) / (this.halfDurationDeg * this.halfDurationDeg);
  }

  /** True while the valve is off its seat. */
  isOpen(thetaDeg: number): boolean {
    const x = this.camX(thetaDeg);
    const xs = this.seatX;
    return x > -xs && x < xs;
  }
}

// ---- lift profile from a ValveSpec ------------------------------------------------------

/** The ValveSpec fields a lift profile needs (any ValveSpec qualifies). */
export type CamValveSpec = ValveTimingSpec & { readonly cam?: CamSpec };

/**
 * Running clearance of a valve: the spec's own ValveSpec.lash when present (contract: it overrides
 * the cycle option), else `optionLash` (CycleModelOptions.valveLash), else 0. m.
 */
export function resolveValveLash(valve: { readonly lash?: number }, optionLash?: number): number {
  return valve.lash ?? optionLash ?? 0;
}

/**
 * Lift profile of a valve as its ValveSpec.cam asks for:
 *  - absent: the default 2-10-18-26 polydyne ValveLiftProfile ({...valve, lash}: bit-identical to the
 *    cycle model's historical `new ValveLiftProfile({ ...spec.intakeValve, lash })`);
 *  - 'polydyne': ValveLiftProfile with the given powers (default 2-10-18-26);
 *  - 'three-arc-flat-follower': exact flat-follower lift of the three-arc lobe, lobe centre at the
 *    midpoint of openDeg/closeDeg (cam-lift.ts ThreeArcFlatFollowerProfile);
 *  - 'table': the measured zero-lash table (cam-lift.ts TabulatedLiftProfile).
 * @param lashOverride running clearance to use, m; absent → valve.lash ?? 0. The cycle model passes
 *   resolveValveLash(valve, options.valveLash) so a spec lash wins over the option.
 */
export function createLiftProfile(valve: CamValveSpec, lashOverride?: number): LiftProfile {
  const lash = lashOverride ?? valve.lash ?? 0;
  const cam = valve.cam;
  if (cam === undefined) return new ValveLiftProfile({ ...valve, lash });
  switch (cam.kind) {
    case 'polydyne':
      return new ValveLiftProfile({ ...valve, lash }, cam.powers === undefined ? POLYDYNE_2_10_18_26 : polydyneCam(cam.powers));
    case 'three-arc-flat-follower':
      return new ThreeArcFlatFollowerProfile({ ...valve, lash }, cam);
    case 'table':
      return new TabulatedLiftProfile(cam, lash);
    default:
      throw new RangeError(`createLiftProfile: unknown cam kind '${(cam as { kind: string }).kind}'`);
  }
}

// ---- functional API (contract) with a per-spec cache ------------------------------------

interface CacheEntry {
  maxLift: number;
  openDeg: number;
  closeDeg: number;
  thr: number;
  lash: number | undefined;
  cam: CamSpec | undefined;
  profile: LiftProfile;
}
const cache = new WeakMap<CamValveSpec, CacheEntry>();

/**
 * Cached {@link LiftProfile} for a spec object (rebuilt if its timing fields, lash or cam object
 * changed): {@link createLiftProfile} with the spec's own optional `lash` (default polydyne when the
 * spec has no `cam`).
 */
export function valveLiftProfile(spec: CamValveSpec): LiftProfile {
  const e = cache.get(spec);
  if (
    e !== undefined &&
    e.maxLift === spec.maxLift &&
    e.openDeg === spec.openDeg &&
    e.closeDeg === spec.closeDeg &&
    e.thr === spec.timingLiftThreshold &&
    e.lash === spec.lash &&
    e.cam === spec.cam
  ) {
    return e.profile;
  }
  const profile = createLiftProfile(spec);
  cache.set(spec, {
    maxLift: spec.maxLift,
    openDeg: spec.openDeg,
    closeDeg: spec.closeDeg,
    thr: spec.timingLiftThreshold,
    lash: spec.lash,
    cam: spec.cam,
    profile,
  });
  return profile;
}

/** Valve lift at crank angle θ (deg, firing-TDC convention), m. */
export function valveLift(spec: CamValveSpec, thetaDeg: number): number {
  return valveLiftProfile(spec).lift(thetaDeg);
}

/** Valve lift rate dL/dθ at crank angle θ (deg), m per crank degree. */
export function valveLiftRate(spec: CamValveSpec, thetaDeg: number): number {
  return valveLiftProfile(spec).liftRate(thetaDeg);
}
