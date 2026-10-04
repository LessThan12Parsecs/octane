/**
 * Valve lift from real cam geometry (CamSpec 'three-arc-flat-follower') and from a measured lift
 * table (CamSpec 'table'). Both implement {@link LiftProfile} (valve-lift.ts), like the default
 * polydyne {@link ValveLiftProfile}; `createLiftProfile` (valve-lift.ts) selects by ValveSpec.cam.
 *
 * ── Three-arc lobe on a flat-faced translating follower ─────────────────────────────────
 * The lobe outline (cam frame, nose along +x, cam axis at the origin) is the convex union of
 *  - the BASE circle, radius r_b, centred on the cam axis;
 *  - the NOSE arc, radius r_n, centre N = (d, 0) with d = r_b + rise − r_n (so the nose tip stands
 *    `rise` above the base circle);
 *  - two FLANK arcs, radius ρ, each tangent to the base circle and to the nose arc from the inside
 *    (the flank circle contains both). Internal tangency to the base circle at the cam angle φ_b puts
 *    the flank centre at C_f = −(ρ − r_b)·(cos φ_b, sin φ_b); tangency to the nose circle requires
 *    |N − C_f| = ρ − r_n, i.e.
 *      cos φ_b = ((ρ − r_n)² − d² − (ρ − r_b)²) / (2 d (ρ − r_b)),
 *    and the flank meets the nose at φ_n = atan2(−C_f,y, d − C_f,x) (the direction of N − C_f).
 * A flat follower face perpendicular to the follower axis touches the outline where the outward
 * normal points along the axis, so the follower displacement is the SUPPORT FUNCTION of the outline
 * in that direction (exact, no approximation): for a circle of centre c and radius r it is
 * c·u + r. With φ the cam angle between the nose axis and the follower axis (|φ| by symmetry):
 *   nose  (|φ| ≤ φ_n):        L_cam = d cos φ + r_n − r_b       = rise − d (1 − cos φ)
 *   flank (φ_n ≤ |φ| ≤ φ_b):  L_cam = (ρ − r_b)(1 − cos(φ_b − |φ|))
 *   base  (|φ| ≥ φ_b):        L_cam = 0
 * (the three-arc / "circular-arc" cam of the cam-design texts, e.g. Rothbart, Cam Design Handbook
 * 2004, and Norton, Cam Design and Manufacturing Handbook 2009 — UNVERIFIED chapter/equation numbers;
 * the formulas above are derived here and checked against a numerical envelope of the polygonised
 * outline, tools/reference/mechanics_cam_three_arc.py). They reproduce the MTFC Tulsa lift equations
 * of the stock Model T lobe (design_stock.htm: nose L = 0.2502 − 0.6250(1 − cos f), flank
 * L = 0.8541[1 − cos(68.57° − |f|)], junction 40.3°). The follower velocity dL/dφ is continuous; its
 * acceleration (ρ − r_b) cos(φ_b − |φ|) on the flank and −d cos φ on the nose steps at both
 * junctions (no quieting ramp). The contact point sits dL/dφ off the follower axis; its maximum,
 * d sin φ_n (at the flank/nose junction), is the smallest admissible follower-face radius.
 *
 * Valve lift = max(0, L_cam − lash): the running clearance is taken up on the flank, so the valve
 * leaves and lands on its seat with the finite flank velocity. The cam turns at half crank speed:
 * φ = (θ − θ_c)/2 with θ_c the lobe-centre crank angle. The lobe is symmetric, so θ_c is placed at
 * the midpoint of the quoted ValveSpec events (openDeg, closeDeg); the events themselves then follow
 * from the geometry and the lash, and `timingErrorDeg` reports how far the quoted half duration
 * differs from the geometric one (a geometric cam cannot be forced to an arbitrary duration).
 *
 * ── Tabulated lift ──────────────────────────────────────────────────────────────────────
 * A zero-lash valve-lift table uniform in crank angle (local cylinder angle), linearly interpolated,
 * periodic in 720°, valve lift = max(0, L_table − lash). With a 1° step this is bit-for-bit the
 * engines/cfr.ts cfrValveLiftChoi2018 evaluation.
 *
 * Units: lift m, angles crank degrees (cam angles inside in rad), lift rate m/deg, lift
 * acceleration m/deg². Hot paths allocate nothing.
 */
import type { CamSpec } from '../core/engine-spec';
import type { LiftProfile, ValveTimingSpec } from './valve-lift';

/** Crank degrees → cam radians (cam turns at half crank speed). */
const CAM_RAD_PER_CRANK_DEG = Math.PI / 360;

/** Wrap crank degrees into [−360, 360) (same arithmetic as valve-lift.ts wrapDeg720). */
function wrap720(deg: number): number {
  let d = (deg + 360) % 720;
  if (d < 0) d += 720;
  return d - 360;
}

/** Three-arc lobe dimensions (CamSpec 'three-arc-flat-follower' without the tag), m. */
export type ThreeArcLobe = Omit<Extract<CamSpec, { kind: 'three-arc-flat-follower' }>, 'kind'>;

/** Derived arc layout of a three-arc lobe (cam frame, nose along +x). */
export interface ThreeArcLobeGeometry {
  /** Nose-arc centre distance from the cam axis d = r_b + rise − r_n, m. */
  noseCenterDistance: number;
  /** Flank-arc centre distance from the cam axis ρ − r_b (on the far side of the axis), m. */
  flankCenterDistance: number;
  /** Cam angle from the nose axis where the flank leaves the base circle, φ_b, rad. */
  baseAngle: number;
  /** Cam angle from the nose axis where the flank meets the nose arc, φ_n, rad. */
  noseAngle: number;
  /** Cam lift at the flank/nose junction, m. */
  junctionLift: number;
  /** Largest follower-contact offset from the follower axis, d sin φ_n, m (minimum face radius). */
  maxContactOffset: number;
}

/**
 * Arc layout of a three-arc lobe from its radii and rise (see file header). Throws when the radii
 * admit no convex lobe tangent to all three arcs.
 */
export function threeArcLobeGeometry(lobe: ThreeArcLobe): ThreeArcLobeGeometry {
  const { baseRadius: rb, flankRadius: rf, noseRadius: rn, rise } = lobe;
  if (!(rb > 0 && rise > 0 && rn > 0)) throw new RangeError('threeArcLobeGeometry: radii and rise must be > 0');
  const d = rb + rise - rn;
  const k = rf - rb;
  if (!(d > 0 && k > 0 && rf > rn)) throw new RangeError('threeArcLobeGeometry: need flankRadius > baseRadius, noseRadius and a nose inside the lobe');
  const cb = ((rf - rn) * (rf - rn) - d * d - k * k) / (2 * d * k);
  if (!(cb > -1 && cb < 1)) throw new RangeError('threeArcLobeGeometry: flank radius cannot join the base and nose circles');
  const phiB = Math.acos(cb);
  // N − C_f with C_f = −k (cos φ_b, sin φ_b)
  const phiN = Math.atan2(k * Math.sin(phiB), d + k * cb);
  if (!(phiN > 0 && phiN < phiB)) throw new RangeError('threeArcLobeGeometry: degenerate flank (junction outside the flank)');
  return {
    noseCenterDistance: d,
    flankCenterDistance: k,
    baseAngle: phiB,
    noseAngle: phiN,
    junctionLift: k * (1 - Math.cos(phiB - phiN)),
    maxContactOffset: d * Math.sin(phiN),
  };
}

/**
 * Exact valve lift of a three-arc lobe on a flat-faced translating follower (mushroom tappet),
 * with running clearance. Lobe centre at the midpoint of the quoted events. Allocation-free.
 */
export class ThreeArcFlatFollowerProfile implements LiftProfile {
  /** Maximum VALVE lift rise − lash, m (geometric; ValveSpec.maxLift is not used). */
  readonly maxLift: number;
  /** Running clearance, m. */
  readonly lash: number;
  /** Cam rise at the follower (lift at zero clearance), m. */
  readonly camMaxLift: number;
  /** Lobe centre (nose) θ_c, crank deg in [−360, 360). */
  readonly centerDeg: number;
  /** Valve seat-off angle, crank deg in [−360, 360). */
  readonly seatOpenDeg: number;
  /** Valve seat-on angle, crank deg in [−360, 360). */
  readonly seatCloseDeg: number;
  /** Half of the seat-to-seat valve duration, crank deg (= 2 φ_seat in degrees). */
  readonly halfDurationDeg: number;
  /** Half of the base-circle-to-base-circle cam duration, crank deg (= 2 φ_b in degrees). */
  readonly camHalfDurationDeg: number;
  /**
   * Geometric minus quoted half duration at the timing threshold, crank deg: the cam reaches the
   * threshold lift this much before openDeg and after closeDeg (negative: after / before).
   */
  readonly timingErrorDeg: number;
  /** Arc layout. */
  readonly geometry: ThreeArcLobeGeometry;
  private readonly d: number;
  private readonly k: number;
  private readonly rnMinusRb: number;
  private readonly phiB: number;
  private readonly phiN: number;
  private readonly phiSeat: number;

  /**
   * @param spec quoted valve events (crank deg) at the valve lift timingLiftThreshold (m) and the
   *   running clearance `lash` (m); maxLift is ignored (geometry decides)
   * @param lobe base/flank/nose radii and rise, m
   */
  constructor(spec: ValveTimingSpec, lobe: ThreeArcLobe) {
    const g = threeArcLobeGeometry(lobe);
    const lash = spec.lash !== undefined && spec.lash > 0 ? spec.lash : 0;
    const thr = spec.timingLiftThreshold > 0 ? spec.timingLiftThreshold : 0;
    if (!(lash + thr < lobe.rise)) throw new RangeError('ThreeArcFlatFollowerProfile: lash + timingLiftThreshold must be < rise');
    const dur = (((spec.closeDeg - spec.openDeg) % 720) + 720) % 720;
    if (!(dur > 0)) throw new RangeError('ThreeArcFlatFollowerProfile: zero event duration');
    this.geometry = g;
    this.d = g.noseCenterDistance;
    this.k = g.flankCenterDistance;
    this.rnMinusRb = lobe.noseRadius - lobe.baseRadius;
    this.phiB = g.baseAngle;
    this.phiN = g.noseAngle;
    this.lash = lash;
    this.camMaxLift = lobe.rise;
    this.maxLift = lobe.rise - lash;
    this.phiSeat = this.camAngleAtLift(lash);
    const phiTh = this.camAngleAtLift(lash + thr);
    this.centerDeg = wrap720(spec.openDeg + dur / 2);
    this.halfDurationDeg = (2 * this.phiSeat) / (Math.PI / 180);
    this.camHalfDurationDeg = (2 * this.phiB) / (Math.PI / 180);
    this.timingErrorDeg = (2 * phiTh) / (Math.PI / 180) - dur / 2;
    this.seatOpenDeg = wrap720(this.centerDeg - this.halfDurationDeg);
    this.seatCloseDeg = wrap720(this.centerDeg + this.halfDurationDeg);
  }

  /** Cam angle |φ| (rad) at which the cam lift falls to `l` (0 ≤ l < rise; l ≤ 0 → φ_b). */
  camAngleAtLift(l: number): number {
    if (!(l > 0)) return this.phiB;
    if (l < this.geometry.junctionLift) return this.phiB - Math.acos(1 - l / this.k);
    return Math.acos((l - this.rnMinusRb) / this.d);
  }

  /** Cam lift at cam angle |φ| (rad), m (zero clearance). */
  camLiftAtAngle(phi: number): number {
    const a = phi < 0 ? -phi : phi;
    if (!(a < this.phiB)) return 0;
    if (a <= this.phiN) return this.d * Math.cos(a) + this.rnMinusRb;
    return this.k * (1 - Math.cos(this.phiB - a));
  }

  /** Signed cam angle φ (rad) from the nose at crank angle θ (deg), in [−π, π) (720° crank wrap). */
  camAngle(thetaDeg: number): number {
    return wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
  }

  /** Valve lift at crank angle θ (deg), m: max(0, L_cam − lash). */
  lift(thetaDeg: number): number {
    let a = wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
    if (a < 0) a = -a;
    if (!(a < this.phiSeat)) return 0;
    const l = (a <= this.phiN ? this.d * Math.cos(a) + this.rnMinusRb : this.k * (1 - Math.cos(this.phiB - a))) - this.lash;
    return l > 0 ? l : 0;
  }

  /** Cam lift at the follower (zero clearance), m: the tappet motion (render). */
  camLift(thetaDeg: number): number {
    return this.camLiftAtAngle(wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG);
  }

  /** dL/dθ of the valve, m per crank degree (0 while seated; jumps at the seat when lash > 0). */
  liftRate(thetaDeg: number): number {
    const x = wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
    const a = x < 0 ? -x : x;
    if (!(a < this.phiSeat)) return 0;
    // dL/d|φ| (≤ 0) × d|φ|/dθ
    const dl = a <= this.phiN ? -this.d * Math.sin(a) : -this.k * Math.sin(this.phiB - a);
    return (x < 0 ? -dl : dl) * CAM_RAD_PER_CRANK_DEG;
  }

  /** d²L/dθ² of the valve, m per crank degree² (0 while seated; steps at the arc junctions). */
  liftAccel(thetaDeg: number): number {
    let a = wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
    if (a < 0) a = -a;
    if (!(a < this.phiSeat)) return 0;
    const dd = a <= this.phiN ? -this.d * Math.cos(a) : this.k * Math.cos(this.phiB - a);
    return dd * CAM_RAD_PER_CRANK_DEG * CAM_RAD_PER_CRANK_DEG;
  }

  /**
   * Follower-face contact offset from the follower axis, m (signed, = dL_cam/dφ): where the lobe
   * touches the mushroom foot (render; |offset| ≤ geometry.maxContactOffset).
   */
  contactOffset(thetaDeg: number): number {
    const x = wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
    const a = x < 0 ? -x : x;
    if (!(a < this.phiB)) return 0;
    const dl = a <= this.phiN ? -this.d * Math.sin(a) : -this.k * Math.sin(this.phiB - a);
    return x < 0 ? -dl : dl;
  }

  /** True while the valve is off its seat. */
  isOpen(thetaDeg: number): boolean {
    let a = wrap720(thetaDeg - this.centerDeg) * CAM_RAD_PER_CRANK_DEG;
    if (a < 0) a = -a;
    return a < this.phiSeat;
  }
}

/** Measured zero-lash lift table (CamSpec 'table' without the tag). */
export type LiftTable = Omit<Extract<CamSpec, { kind: 'table' }>, 'kind'>;

/**
 * Valve lift from a zero-lash table uniform in crank angle, linear interpolation, periodic in 720°,
 * with running clearance. Allocation-free. Tables should start and end seated (the lift is 0
 * outside the table span).
 */
export class TabulatedLiftProfile implements LiftProfile {
  /** Maximum valve lift max(table) − lash, m. */
  readonly maxLift: number;
  /** Running clearance, m. */
  readonly lash: number;
  /** Maximum zero-lash lift, m. */
  readonly camMaxLift: number;
  /** Crank angle of the table maximum, deg in [−360, 360). */
  readonly centerDeg: number;
  /** First angle where the zero-lash lift exceeds the lash, deg in [−360, 360). */
  readonly seatOpenDeg: number;
  /** Last angle where it does, deg in [−360, 360). */
  readonly seatCloseDeg: number;
  /** First tabulated angle, deg. */
  readonly startDeg: number;
  /** Table step, deg. */
  readonly stepDeg: number;
  private readonly t: Float64Array;
  private readonly nSeg: number;

  /** @param table startDeg, stepDeg (crank deg) and zero-lash lift (m) @param lash running clearance, m */
  constructor(table: LiftTable, lash = 0) {
    const n = table.zeroLashLift.length;
    const step = table.stepDeg;
    if (!(n >= 2)) throw new RangeError('TabulatedLiftProfile: need ≥ 2 points');
    if (!(step > 0) || !((n - 1) * step <= 720)) throw new RangeError('TabulatedLiftProfile: table must span (0, 720] deg');
    if (!Number.isFinite(table.startDeg)) throw new RangeError('TabulatedLiftProfile: startDeg must be finite');
    const lsh = lash > 0 ? lash : 0;
    this.t = Float64Array.from(table.zeroLashLift);
    this.nSeg = n - 1;
    this.startDeg = table.startDeg;
    this.stepDeg = step;
    this.lash = lsh;
    const t = this.t;
    let iMax = 0;
    for (let i = 1; i < n; i++) if (t[i] > t[iMax]) iMax = i;
    if (!(t[iMax] > lsh)) throw new RangeError('TabulatedLiftProfile: lash exceeds the tabulated lift');
    this.camMaxLift = t[iMax];
    this.maxLift = t[iMax] - lsh;
    this.centerDeg = wrap720(table.startDeg + iMax * step);
    // seat events: first upward / last downward crossing of the lash (linear segments)
    let open = table.startDeg;
    if (!(t[0] > lsh)) {
      for (let i = 0; i < n - 1; i++) {
        if (t[i + 1] > lsh) {
          open = table.startDeg + step * (i + (lsh - t[i]) / (t[i + 1] - t[i]));
          break;
        }
      }
    }
    let close = table.startDeg + (n - 1) * step;
    if (!(t[n - 1] > lsh)) {
      for (let i = n - 1; i > 0; i--) {
        if (t[i - 1] > lsh) {
          close = table.startDeg + step * (i - (lsh - t[i]) / (t[i - 1] - t[i]));
          break;
        }
      }
    }
    this.seatOpenDeg = wrap720(open);
    this.seatCloseDeg = wrap720(close);
  }

  /** Zero-lash (cam) lift at crank angle θ (deg), m. */
  camLift(thetaDeg: number): number {
    let u = (thetaDeg - this.startDeg) % 720;
    if (u < 0) u += 720;
    const s = u / this.stepDeg;
    if (!(s < this.nSeg)) return 0;
    const i = Math.floor(s);
    const t = this.t;
    return t[i] + (s - i) * (t[i + 1] - t[i]);
  }

  /** Valve lift at crank angle θ (deg), m: max(0, L_table − lash). */
  lift(thetaDeg: number): number {
    // same arithmetic as engines/cfr.ts cfrValveLiftChoi2018 (bit-identical for stepDeg = 1)
    let u = (thetaDeg - this.startDeg) % 720;
    if (u < 0) u += 720;
    const s = u / this.stepDeg;
    if (!(s < this.nSeg)) return 0;
    const i = Math.floor(s);
    const t = this.t;
    const l0 = t[i] + (s - i) * (t[i + 1] - t[i]);
    return l0 > this.lash ? l0 - this.lash : 0;
  }

  /** dL/dθ, m per crank degree (segment slope while open, else 0). */
  liftRate(thetaDeg: number): number {
    let u = (thetaDeg - this.startDeg) % 720;
    if (u < 0) u += 720;
    const s = u / this.stepDeg;
    if (!(s < this.nSeg)) return 0;
    const i = Math.floor(s);
    const t = this.t;
    const l0 = t[i] + (s - i) * (t[i + 1] - t[i]);
    return l0 > this.lash ? (t[i + 1] - t[i]) / this.stepDeg : 0;
  }

  /** d²L/dθ²: 0 (piecewise-linear table; the curvature sits in the nodes). */
  liftAccel(_thetaDeg: number): number {
    return 0;
  }

  /** True while the valve is off its seat. */
  isOpen(thetaDeg: number): boolean {
    return this.lift(thetaDeg) > 0;
  }
}
