/**
 * Cam geometry for the Model T mechanism (render side, pure math).
 *
 * The Ford lobe is a convex three-arc ("circular-arc", "harmonic") cam acting on a flat-footed
 * mushroom tappet (engine-spec.ts CamSpec 'three-arc-flat-follower'): a base circle r_b, two flank
 * arcs of radius r_f internally tangent to the base circle, and a nose arc r_n internally tangent to
 * both flanks, with total rise = r_b + ... at the nose. The lift of a flat follower is the SUPPORT
 * FUNCTION h(f) of the lobe outline in the follower direction minus r_b (exact for any convex lobe):
 *
 *   nose   |f| ≤ f₁:   h = d·cos f + r_n,                 d = r_b + rise − r_n (nose-centre distance)
 *   flank f₁ ≤ |f| ≤ f₀: h = r_f − a·cos(f₀ − |f|),        a = r_f − r_b
 *   base   |f| ≥ f₀:   h = r_b
 *
 * f = cam angle from the lobe centreline. Tangency fixes the angles: the flank centre sits at
 * −a·e(f₀), so |a·e(f₀) + d·e(0)| = r_f − r_n  →  cos f₀ = ((r_f − r_n)² − a² − d²)/(2ad), and the
 * flank/nose junction is the direction of (C_n − C_f): f₁ = atan2(a sin f₀, d + a cos f₀).
 * These are the MTFC Tulsa lift equations for the stock 1913–27 lobe (design_stock.htm: nose
 * L = 0.2502 − 0.6250(1 − cos f) for |f| < 40.3°, flank L = 0.8541[1 − cos(68.57° − |f|)]), derived
 * here from the radii instead of the rounded constants (cam.test.ts checks both).
 *
 * The lobe centre is derived from the quoted valve timing (engine-spec.ts CamSpec: "its lobe centre is
 * derived from them"): the midpoint of the open → close interval, the lobe being symmetric
 * [Good22 p. 33: intake and exhaust lobes identical and symmetric].
 */
import type { CamSpec, ValveSpec } from '../../physics/core/engine-spec';
import { modelValveLift, wrapDeg720 } from '../engine/cam-profile';

export interface ThreeArcGeometry {
  baseRadius: number;
  flankRadius: number;
  noseRadius: number;
  rise: number;
  /** Nose-centre distance from the cam axis, m. */
  d: number;
  /** r_f − r_b, m. */
  a: number;
  /** Cam angle (rad) where the flank meets the base circle. */
  f0: number;
  /** Cam angle (rad) where the nose meets the flank. */
  f1: number;
}

export type ThreeArcCam = Extract<CamSpec, { kind: 'three-arc-flat-follower' }>;

/** Tangency construction of a three-arc lobe; throws for radii that cannot close a convex lobe. */
export function threeArcGeometry(c: Pick<ThreeArcCam, 'baseRadius' | 'flankRadius' | 'noseRadius' | 'rise'>): ThreeArcGeometry {
  const rb = c.baseRadius, rf = c.flankRadius, rn = c.noseRadius;
  const d = rb + c.rise - rn;
  const a = rf - rb;
  if (!(rb > 0 && rn > 0 && a > 0 && d > 0 && rf > rn)) throw new RangeError('threeArcGeometry: invalid radii');
  const cf0 = ((rf - rn) ** 2 - a * a - d * d) / (2 * a * d);
  if (!(cf0 > -1 && cf0 < 1)) throw new RangeError('threeArcGeometry: flank cannot be tangent to base and nose');
  const f0 = Math.acos(cf0);
  const f1 = Math.atan2(a * Math.sin(f0), d + a * Math.cos(f0));
  return { baseRadius: rb, flankRadius: rf, noseRadius: rn, rise: c.rise, d, a, f0, f1 };
}

/** Support function h(f) of the lobe (distance of the follower face from the cam axis), m. f in rad. */
export function threeArcSupport(g: ThreeArcGeometry, f: number): number {
  let x = Math.abs(f) % (2 * Math.PI);
  if (x > Math.PI) x = 2 * Math.PI - x;
  if (x <= g.f1) return g.d * Math.cos(x) + g.noseRadius;
  if (x <= g.f0) return g.flankRadius - g.a * Math.cos(g.f0 - x);
  return g.baseRadius;
}

/** Gross follower (tappet) lift h(f) − r_b, m. */
export function threeArcLift(g: ThreeArcGeometry, f: number): number {
  return threeArcSupport(g, f) - g.baseRadius;
}

/** Lobe-centre crank angle (local cylinder angle, deg): midpoint of the open → close interval. */
export function lobeCentreDeg(v: Pick<ValveSpec, 'openDeg' | 'closeDeg'>): number {
  const dur = (((v.closeDeg - v.openDeg) % 720) + 720) % 720 || 720;
  return wrapDeg720(v.openDeg + dur / 2);
}

/**
 * Tappet (follower) lift of a valve's cam at the local crank angle θ, m, and the valve lift that goes
 * with it. For a three-arc cam both are exact geometry (valve = max(0, tappet − lash)); for any other
 * cam kind the render falls back to the sin² model curve of the ValveSpec with zero lash.
 */
export class ValveCam {
  readonly geometry: ThreeArcGeometry | null;
  readonly centreDeg: number;
  readonly lash: number;
  readonly baseRadius: number;

  constructor(readonly spec: ValveSpec, fallbackBaseRadius: number) {
    const c = spec.cam;
    this.geometry = c && c.kind === 'three-arc-flat-follower' ? threeArcGeometry(c) : null;
    this.centreDeg = lobeCentreDeg(spec);
    this.lash = this.geometry ? Math.max(0, spec.lash ?? 0) : 0;
    this.baseRadius = this.geometry ? this.geometry.baseRadius : fallbackBaseRadius;
  }

  /** Tappet lift at crank angle offset `dThetaDeg` from the lobe centre (any value; wraps), m. */
  tappetLiftFromCentre(dThetaDeg: number): number {
    const d = wrapDeg720(dThetaDeg);
    if (this.geometry) return threeArcLift(this.geometry, (d * Math.PI) / 360);
    return modelValveLift(this.spec, d + this.centreDeg);
  }

  /** Tappet lift at the local crank angle θ (deg), m. */
  tappetLift(thetaLocalDeg: number): number {
    return this.tappetLiftFromCentre(thetaLocalDeg - this.centreDeg);
  }

  /** Valve lift at the local crank angle θ (deg), m. */
  valveLift(thetaLocalDeg: number): number {
    return Math.max(0, this.tappetLift(thetaLocalDeg) - this.lash);
  }
}
