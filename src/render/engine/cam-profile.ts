/**
 * Cam-lobe geometry consistent with the valve lifts the physics reports.
 *
 * The valves (and, through the rocker, the tappets) are driven ONLY by the
 * snapshot lifts. The camshaft angle is derived from θ (half crank speed).
 * For the rendered lobe to actually touch its flat-faced tappet, its outline
 * must be the flat-follower envelope of the lift curve:
 *
 *   support function p(φ) = R_b + s(θ(φ)),   θ = 180° − 2φ   (cam-local angle φ,
 *                                                             follower along world +y)
 *
 * The lift curve s(θ) is first a model from the ValveSpec (open/close angles,
 * max lift) and is then replaced by the curve LEARNED from the snapshot stream
 * (LiftProfileLearner), so the lobe matches whatever profile the physics uses.
 *
 * The lobe outline is built as the intersection of the half-planes
 * {x · n(φ_j) ≤ p(φ_j)} (dual convex hull), which is exact for convex
 * (realisable) profiles and degrades gracefully for noisy data.
 */
import type { ValveSpec } from '../../physics/core/engine-spec';

/** Wrap crank degrees into [−360, 360). */
export function wrapDeg720(deg: number): number {
  let d = (deg + 360) % 720;
  if (d < 0) d += 720;
  return d - 360;
}

/**
 * Model valve lift (m) from the ValveSpec: a sin² (cosine) profile between
 * the quoted open/close angles, extended so the lift equals
 * `timingLiftThreshold` at the quoted angles. Render fallback only.
 */
export function modelValveLift(v: ValveSpec, thetaDeg: number): number {
  const dur = ((v.closeDeg - v.openDeg) % 720 + 720) % 720 || 720;
  const thr = Math.min(Math.max(v.timingLiftThreshold, 0), 0.9 * v.maxLift);
  // Extend the duration so that lift(open) = thr: sin²(π δ / D') = thr/max.
  let ext = 0;
  if (thr > 0) {
    const a = Math.asin(Math.sqrt(thr / v.maxLift)) / Math.PI;
    // D' = dur + 2 δ, δ = a D'  ⇒  D' = dur / (1 − 2a)
    ext = (dur / (1 - 2 * a) - dur) / 2;
  }
  const open = v.openDeg - ext;
  const D = dur + 2 * ext;
  const t = ((((thetaDeg - open) % 720) + 720) % 720) / D;
  if (t >= 1) return 0;
  const s = Math.sin(Math.PI * t);
  return v.maxLift * s * s;
}

/**
 * Learns lift(θ) from a stream of (θ, lift) samples, one node per 1° bin
 * (the sample closest to the bin centre is kept). Linear interpolation
 * between nodes. Detects a changed profile and restarts.
 */
export class LiftProfileLearner {
  static readonly BINS = 720;
  private readonly theta = new Float64Array(LiftProfileLearner.BINS);
  private readonly lift = new Float64Array(LiftProfileLearner.BINS);
  private readonly filled = new Uint8Array(LiftProfileLearner.BINS);
  private nFilled = 0;
  private maxGapCache = 720;
  private gapDirty = true;
  private mismatch = 0;
  /** Incremented when the learned curve is reset (profile changed). */
  resets = 0;

  constructor(
    /** Max node spacing (crank deg) for the profile to count as usable. */
    readonly readyGapDeg = 8,
    /** Deviation from the learned curve (m) that counts as a mismatch. */
    readonly mismatchTol = 1.5e-4,
  ) {}

  reset(): void {
    this.filled.fill(0);
    this.nFilled = 0;
    this.gapDirty = true;
    this.maxGapCache = 720;
    this.mismatch = 0;
  }

  add(thetaDeg: number, lift: number): void {
    if (!Number.isFinite(thetaDeg) || !Number.isFinite(lift)) return;
    const th = wrapDeg720(thetaDeg);
    const l = Math.max(0, lift);
    if (this.isReady()) {
      // Once learned, samples that disagree are not stored (protects against
      // e.g. an interpolated sample straddling the ±360° wrap); a sustained
      // disagreement means the profile changed → start over.
      const pred = this.sample(th);
      if (Math.abs(pred - l) > this.mismatchTol) {
        if (++this.mismatch > 24) {
          this.reset();
          this.resets++;
        } else {
          return;
        }
      } else if (this.mismatch > 0) {
        this.mismatch--;
      }
    }
    let b = Math.floor(th + 360);
    if (b >= LiftProfileLearner.BINS) b = LiftProfileLearner.BINS - 1;
    const c = b - 360 + 0.5;
    if (!this.filled[b]) {
      this.filled[b] = 1;
      this.nFilled++;
      this.theta[b] = th;
      this.lift[b] = l;
      this.gapDirty = true;
    } else if (Math.abs(th - c) <= Math.abs(this.theta[b] - c)) {
      this.theta[b] = th;
      this.lift[b] = l;
    }
  }

  /** Largest spacing between consecutive nodes, crank degrees (cyclic). */
  maxGapDeg(): number {
    if (!this.gapDirty) return this.maxGapCache;
    this.gapDirty = false;
    if (this.nFilled < 2) return (this.maxGapCache = 720);
    let first = -1, prev = -1, maxGap = 0;
    for (let b = 0; b < LiftProfileLearner.BINS; b++) {
      if (!this.filled[b]) continue;
      if (first < 0) first = b;
      else maxGap = Math.max(maxGap, this.theta[b] - this.theta[prev]);
      prev = b;
    }
    maxGap = Math.max(maxGap, this.theta[first] + 720 - this.theta[prev]);
    return (this.maxGapCache = maxGap);
  }

  isReady(): boolean {
    return this.maxGapDeg() <= this.readyGapDeg;
  }

  /** 0 = not usable yet, 1 = usable (gaps ≤ readyGapDeg), 2 = every 1° bin filled. */
  stage(): 0 | 1 | 2 {
    if (!this.isReady()) return 0;
    return this.nFilled === LiftProfileLearner.BINS ? 2 : 1;
  }

  get filledBins(): number {
    return this.nFilled;
  }

  /** Learned lift at θ (linear between nodes). Only meaningful when ready. */
  sample(thetaDeg: number): number {
    const th = wrapDeg720(thetaDeg);
    const N = LiftProfileLearner.BINS;
    let b = Math.floor(th + 360);
    if (b >= N) b = N - 1;
    // previous node with θ_i ≤ th (cyclic), next node with θ_j > th
    let i = b, iθ = 0, found = false;
    for (let k = 0; k < N; k++) {
      const bb = (b - k + N) % N;
      if (this.filled[bb]) {
        let t = this.theta[bb];
        if (k === 0 && t > th) continue;
        if (bb > b || (k > 0 && t > th)) t -= 720;
        i = bb; iθ = t; found = true;
        break;
      }
    }
    if (!found) return 0;
    let j = b, jθ = 0;
    found = false;
    for (let k = 0; k < N; k++) {
      const bb = (b + k) % N;
      if (this.filled[bb]) {
        let t = this.theta[bb];
        if (k === 0 && t <= th) continue;
        if (bb < b || (k > 0 && t <= th)) t += 720;
        j = bb; jθ = t; found = true;
        break;
      }
    }
    if (!found) return this.lift[i];
    const span = jθ - iθ;
    if (span <= 1e-12) return this.lift[i];
    const w = (th - iθ) / span;
    return this.lift[i] * (1 - w) + this.lift[j] * w;
  }
}

export interface LobeProfile {
  /** Closed outline, cam-local x,y pairs (counter-clockwise). */
  outline: Float64Array;
  /** Base-circle radius actually used (may exceed the requested one for convexity). */
  baseRadius: number;
  /** Max lateral offset of the contact point on the flat tappet face, m. */
  maxContactOffset: number;
}

/**
 * Cam-local polar angle (rad) whose outward normal points at the follower
 * when the crank is at θ (deg), for a cam turning at +θ/2 about +z with the
 * follower along +y.
 */
export function followerPhi(thetaDeg: number): number {
  return Math.PI / 2 - (thetaDeg * Math.PI) / 360;
}

/**
 * Build a flat-follower cam lobe from a tappet-lift function s(θ) (m).
 * @param minCurvatureRadius minimum radius of curvature enforced by growing R_b.
 */
export function buildLobeProfile(
  tappetLift: (thetaDeg: number) => number,
  baseRadius: number,
  n = 720,
  minCurvatureRadius = 0.003,
): LobeProfile {
  const dphi = (2 * Math.PI) / n;
  const s = new Float64Array(n);
  for (let j = 0; j < n; j++) {
    const phi = j * dphi;
    const theta = wrapDeg720(180 - (2 * phi * 180) / Math.PI);
    s[j] = Math.max(0, tappetLift(theta));
  }
  // Convexity: ρ = p + p'' ≥ ρ_min, with p'' from a wide stencil (robust to node noise).
  const k = Math.max(1, Math.round(n / 240));
  const h = k * dphi;
  let minRho = Infinity;
  for (let j = 0; j < n; j++) {
    const spp = (s[(j + k) % n] - 2 * s[j] + s[(j - k + n) % n]) / (h * h);
    minRho = Math.min(minRho, s[j] + spp);
  }
  let Rb = baseRadius;
  if (Rb + minRho < minCurvatureRadius) Rb = minCurvatureRadius - minRho;

  // Dual points q_j = n_j / p_j; the convex hull picks the active half-planes.
  const qx = new Float64Array(n), qy = new Float64Array(n);
  for (let j = 0; j < n; j++) {
    const phi = j * dphi;
    const p = Rb + s[j];
    qx[j] = Math.cos(phi) / p;
    qy[j] = Math.sin(phi) / p;
  }
  const hull = convexHullIndices(qx, qy);
  // Vertices = intersections of consecutive active lines.
  const m = hull.length;
  const outline = new Float64Array(2 * m);
  for (let a = 0; a < m; a++) {
    const i = hull[a], j = hull[(a + 1) % m];
    const pi = Rb + s[i], pj = Rb + s[j];
    const ci = Math.cos(i * dphi), si = Math.sin(i * dphi);
    const cj = Math.cos(j * dphi), sj = Math.sin(j * dphi);
    const det = ci * sj - si * cj;
    let x: number, y: number;
    if (Math.abs(det) < 1e-12) {
      x = pi * ci; y = pi * si;
    } else {
      x = (pi * sj - pj * si) / det;
      y = (ci * pj - cj * pi) / det;
    }
    outline[2 * a] = x;
    outline[2 * a + 1] = y;
  }
  // Contact offset on the tappet face = p'(φ).
  let maxOff = 0;
  for (let j = 0; j < n; j++) {
    const dp = (s[(j + 1) % n] - s[(j - 1 + n) % n]) / (2 * dphi);
    maxOff = Math.max(maxOff, Math.abs(dp));
  }
  return { outline, baseRadius: Rb, maxContactOffset: maxOff };
}

/** Andrew's monotone chain; returns hull vertex indices in CCW order. */
function convexHullIndices(x: Float64Array, y: Float64Array): number[] {
  const n = x.length;
  const idx = Array.from({ length: n }, (_, i) => i);
  idx.sort((a, b) => (x[a] - x[b]) || (y[a] - y[b]));
  const cross = (o: number, a: number, b: number) =>
    (x[a] - x[o]) * (y[b] - y[o]) - (y[a] - y[o]) * (x[b] - x[o]);
  const lower: number[] = [];
  for (const i of idx) {
    while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], i) <= 0) lower.pop();
    lower.push(i);
  }
  const upper: number[] = [];
  for (let t = n - 1; t >= 0; t--) {
    const i = idx[t];
    while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], i) <= 0) upper.pop();
    upper.push(i);
  }
  upper.pop();
  lower.pop();
  const hull = lower.concat(upper);
  // Order by polar angle of the direction (so the outline goes CCW starting near φ=0).
  hull.sort((a, b) => a - b);
  return hull;
}

/**
 * Height of the lobe along +y after rotating by `camAngle` about +z, i.e.
 * the tappet-face height above the cam axis (for tests / diagnostics).
 */
export function lobeHeightAt(outline: Float64Array, camAngle: number): number {
  const c = Math.cos(camAngle), s = Math.sin(camAngle);
  let best = -Infinity;
  for (let i = 0; i < outline.length; i += 2) {
    const y = s * outline[i] + c * outline[i + 1];
    if (y > best) best = y;
  }
  return best;
}
