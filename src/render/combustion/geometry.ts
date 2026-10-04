/**
 * Pure geometry / physics helpers for the in-cylinder visuals (no three.js).
 * Everything in the CYLINDER frame: origin at the centre of the head
 * fire-deck face, +y toward the head, gas in x²+z² ≤ R², −h ≤ y ≤ 0.
 *
 * The GLSL in volume-shader.ts mirrors `rayChamberInterval`, `raySphere`,
 * `besselJ1`, `axialModeShape` and the brush model; tests exercise the
 * TypeScript versions. The disc helpers here are the flat-disc special case;
 * the general chamber shape (the L-head's bore column ∪ valve pocket) and its
 * ray intervals live in chamber.ts.
 */
import type { ValveSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import {
  ALPHA_10,
  INTEGRAL_SCALE_FRACTION,
  KNOCK_GAMMA,
  NU_AIR_REF,
} from './constants';

export type V3 = [number, number, number];

// ---------------------------------------------------------------------------
// Inside tests
// ---------------------------------------------------------------------------

/** True if (x, y, z) is inside the gas disc of radius R and height h. */
export function insideChamber(x: number, y: number, z: number, R: number, h: number): boolean {
  return x * x + z * z <= R * R && y <= 0 && y >= -h;
}

/** True if (x, y, z) is inside the flame sphere (centre c, radius r). */
export function insideFlameSphere(x: number, y: number, z: number, c: ArrayLike<number>, r: number): boolean {
  const dx = x - c[0], dy = y - c[1], dz = z - c[2];
  return r > 0 && dx * dx + dy * dy + dz * dz <= r * r;
}

/** Burned region = flame sphere ∩ chamber (what the entrainment model burns). */
export function insideBurnedRegion(
  x: number, y: number, z: number, R: number, h: number, c: ArrayLike<number>, r: number,
): boolean {
  return insideChamber(x, y, z, R, h) && insideFlameSphere(x, y, z, c, r);
}

/** Signed distance to the (smooth) flame sphere surface; negative inside. */
export function flameSignedDistance(x: number, y: number, z: number, c: ArrayLike<number>, r: number): number {
  const dx = x - c[0], dy = y - c[1], dz = z - c[2];
  return Math.sqrt(dx * dx + dy * dy + dz * dz) - r;
}

// ---------------------------------------------------------------------------
// Ray intersections (mirrored in GLSL)
// ---------------------------------------------------------------------------

/**
 * Parametric interval [t0, t1] (t ≥ 0) of a ray inside the chamber disc,
 * or null if it misses. `rd` need not be normalised (t is in units of |rd|).
 */
export function rayChamberInterval(
  ro: ArrayLike<number>, rd: ArrayLike<number>, R: number, h: number,
): [number, number] | null {
  let t0 = -Infinity, t1 = Infinity;
  const a = rd[0] * rd[0] + rd[2] * rd[2];
  const b = 2 * (ro[0] * rd[0] + ro[2] * rd[2]);
  const c = ro[0] * ro[0] + ro[2] * ro[2] - R * R;
  if (a < 1e-18) {
    if (c > 0) return null;
  } else {
    const disc = b * b - 4 * a * c;
    if (disc < 0) return null;
    const s = Math.sqrt(disc);
    t0 = (-b - s) / (2 * a);
    t1 = (-b + s) / (2 * a);
  }
  if (Math.abs(rd[1]) < 1e-18) {
    if (ro[1] > 0 || ro[1] < -h) return null;
  } else {
    let y0 = (-h - ro[1]) / rd[1];
    let y1 = (0 - ro[1]) / rd[1];
    if (y0 > y1) { const tmp = y0; y0 = y1; y1 = tmp; }
    t0 = Math.max(t0, y0);
    t1 = Math.min(t1, y1);
  }
  t0 = Math.max(t0, 0);
  return t0 < t1 ? [t0, t1] : null;
}

/** Ray-sphere entry/exit parameters (rd normalised), or null. */
export function raySphere(ro: ArrayLike<number>, rd: ArrayLike<number>, c: ArrayLike<number>, r: number): [number, number] | null {
  const ox = ro[0] - c[0], oy = ro[1] - c[1], oz = ro[2] - c[2];
  const b = ox * rd[0] + oy * rd[1] + oz * rd[2];
  const cc = ox * ox + oy * oy + oz * oz - r * r;
  const hh = b * b - cc;
  if (hh < 0) return null;
  const s = Math.sqrt(hh);
  return [-b - s, -b + s];
}

// ---------------------------------------------------------------------------
// Turbulent flame brush
// ---------------------------------------------------------------------------

export interface BrushInputs {
  /** Turbulence intensity u', m/s. */
  uPrime: number;
  /** Laminar burning velocity S_L, m/s. */
  laminarSpeed: number;
  /** Clearance height h, m. */
  clearanceHeight: number;
  /** Flame radius, m. */
  radius: number;
  /** Cylinder pressure, Pa. */
  pressure: number;
  /** Unburned temperature, K. */
  unburnedTemperature: number;
}

export interface BrushParams {
  /** Integral length scale L_I, m. */
  integralScale: number;
  /** Taylor microscale λ, m. */
  taylorScale: number;
  /** Laminar flame (thermal) thickness δ_L, m. */
  laminarThickness: number;
  /** Total brush thickness δ (reaction zone from leading edge back), m. */
  thickness: number;
  /** Max radial excursion of the wrinkled leading edge, m. */
  wrinkleAmplitude: number;
  /** Wavelength of the dominant wrinkles (≈ L_I), m. */
  wrinkleWavelength: number;
  /** Kernel development factor 1 − exp(−r/L_I) (0 = laminar kernel, 1 = fully turbulent). */
  development: number;
}

export function emptyBrush(): BrushParams {
  return {
    integralScale: 0, taylorScale: 0, laminarThickness: 0, thickness: 0,
    wrinkleAmplitude: 0, wrinkleWavelength: 1, development: 0,
  };
}

/**
 * Brush model consistent with the Blizard–Keck / Tabaczynski entrainment
 * model the physics uses:
 *  - L_I = C·h (INTEGRAL_SCALE_FRACTION),
 *  - ν_u ≈ ν_air(300 K, 1 bar)·(T_u/300)^1.7·(1 bar/p),
 *  - λ = L_I·sqrt(15 / Re_L), Re_L = u' L_I / ν_u (Tabaczynski et al. 1977),
 *  - laminar thickness δ_L ≈ 2 α/S_L with α ≈ 1.4 ν (Pr ≈ 0.7),
 *  - burn-up zone: eddies of size λ burn in τ_b = λ/S_L while the front moves
 *    on at (u' + S_L), so δ_b = (u' + S_L) τ_b = λ (1 + u'/S_L),
 *  - kernel development (Herweg & Maly 1992): turbulence only wrinkles the
 *    kernel once it is comparable with L_I: factor 1 − exp(−r/L_I),
 *  - wrinkle amplitude = δ_b/2 (the brush *is* the envelope of excursions).
 * δ is capped at 2 L_I and at the radius.
 */
export function brushParams(inp: BrushInputs, out: BrushParams = emptyBrush()): BrushParams {
  const h = Math.max(inp.clearanceHeight, 1e-4);
  const L = INTEGRAL_SCALE_FRACTION * h;
  const SL = Math.max(inp.laminarSpeed, 1e-3);
  const up = Math.max(inp.uPrime, 0);
  const Tu = inp.unburnedTemperature > 0 ? inp.unburnedTemperature : 300;
  const p = inp.pressure > 0 ? inp.pressure : 1e5;
  const nu = NU_AIR_REF * Math.pow(Tu / 300, 1.7) * (1e5 / p);
  const reL = (up * L) / nu;
  const lambda = reL > 1e-6 ? L * Math.sqrt(15 / reL) : L;
  const dL = Math.min((2 * 1.4 * nu) / SL, 2e-3);
  const dev = inp.radius > 0 ? 1 - Math.exp(-inp.radius / L) : 0;
  const db = up > 0 ? Math.min(lambda, L) * (1 + up / SL) * dev : 0;
  const r = Math.max(inp.radius, 0);
  const thick = Math.min(dL + db, 2 * L, Math.max(r, dL));
  out.integralScale = L;
  out.taylorScale = lambda;
  out.laminarThickness = dL;
  out.thickness = thick;
  out.wrinkleAmplitude = Math.min(0.5 * db, 0.5 * r);
  out.wrinkleWavelength = L;
  out.development = dev;
  return out;
}

// ---------------------------------------------------------------------------
// Acoustics (knock)
// ---------------------------------------------------------------------------

/** Bessel J1 by its power series to x¹¹ (error < 1e-7 for |x| ≤ 2; used on [0, α_10]). */
export function besselJ1(x: number): number {
  const x2 = x * x;
  return x * (0.5 - x2 * (1 / 16 - x2 * (1 / 384 - x2 * (1 / 18432 - x2 * (1 / 1474560 - x2 / 176947200)))));
}

/** Draper (1938) first circumferential mode frequency f_10 = α_10 c / (π B), Hz. */
export function knockModeFrequency(soundSpeed: number, bore: number): number {
  return (ALPHA_10 * soundSpeed) / (Math.PI * bore);
}

const SPECIES_M: Record<string, number> = {
  CO2: 44.0095e-3, H2O: 18.01528e-3, CO: 28.0101e-3, O2: 31.9988e-3, H2: 2.01588e-3,
  OH: 17.00734e-3, H: 1.00794e-3, O: 15.9994e-3, NO: 30.0061e-3, N2: 28.0134e-3,
};

/** Mean molar mass (kg/mol) of the snapshot's composition (normalised over the listed species). */
export function molarMassOf(comp: EngineSnapshot['burnedComposition']): number {
  let sx = 0, sm = 0;
  for (const k in SPECIES_M) {
    const x = (comp as unknown as Record<string, number>)[k] ?? 0;
    if (x > 0) { sx += x; sm += x * SPECIES_M[k]; }
  }
  return sx > 0.5 ? sm / sx : 0.0286;
}

/** Ideal-gas sound speed, m/s. */
export function soundSpeed(T: number, molarMass: number, gamma = KNOCK_GAMMA): number {
  return Math.sqrt((gamma * 8.31446261815324 * Math.max(T, 1)) / molarMass);
}

/** Normalised (1,0) mode shape at (x, z): J1(α r/R)/J1(α) · cos(φ − φ0). */
export function knockModeShape(x: number, z: number, R: number, axisAngle: number): number {
  const r = Math.min(Math.sqrt(x * x + z * z) / R, 1);
  return (besselJ1(ALPHA_10 * r) / besselJ1(ALPHA_10)) * Math.cos(Math.atan2(z, x) - axisAngle);
}

/**
 * Lowest non-trivial acoustic mode of a rigid-walled duct of length L (Neumann ends): p ∝ cos(π s/L),
 * f = c/(2L) (Rayleigh, The Theory of Sound, 1896, §255; Kinsler et al., Fundamentals of Acoustics,
 * 4th ed. 2000, §9.2), Hz. Used for the L-head footprint (bore + valve pocket, length L along the
 * bore → pocket axis): exact for a rectangular footprint of uniform depth; UNVERIFIED for the real
 * chamber, whose depth differs between the bore column and the pocket (the physics' own L-head modes,
 * when they exist, should replace it).
 */
export function axialModeFrequency(soundSpeed: number, length: number): number {
  return soundSpeed / (2 * length);
}

/** Normalised axial mode cos(π (s − s0)/L) along the unit axis (ex, ez), clamped to [s0, s0 + L] (mirrored in GLSL). */
export function axialModeShape(x: number, z: number, ex: number, ez: number, s0: number, L: number): number {
  const u = Math.min(Math.max((x * ex + z * ez - s0) / L, 0), 1);
  return Math.cos(Math.PI * u);
}

// ---------------------------------------------------------------------------
// Valves / flow
// ---------------------------------------------------------------------------

/** The open (unshrouded) arc of a valve curtain: centre angle and half width, rad. */
export function openArc(v: ValveSpec): { center: number; halfWidth: number } {
  const shroud = Math.min(Math.max(v.shroudArcDeg, 0), 359);
  return { center: v.shroudDirection, halfWidth: (Math.PI * (360 - shroud)) / 360 };
}

/**
 * Mean swirl moment arm of the curtain jets about the cylinder axis, m:
 * angular momentum flux about +y = ṁ · v_h · arm, where v_h is the
 * horizontal jet speed. Jets leave radially from the open arc; the valve
 * offset gives (r_v × ê(ψ))_y = z_v cos ψ − x_v sin ψ, averaged over the arc.
 * An unshrouded (360°) valve gives 0.
 */
export function swirlArm(v: ValveSpec): number {
  const { center, halfWidth } = openArc(v);
  if (halfWidth >= Math.PI - 1e-9) return 0;
  const sinc = halfWidth > 1e-9 ? Math.sin(halfWidth) / halfWidth : 1;
  const [x, z] = v.position;
  return sinc * (z * Math.cos(center) - x * Math.sin(center));
}

/** Geometric curtain area of the open arc, m² (π D_v L × open fraction). */
export function openCurtainArea(v: ValveSpec, lift: number): number {
  const { halfWidth } = openArc(v);
  return Math.PI * v.seatInnerDiameter * Math.max(lift, 0) * (halfWidth / Math.PI);
}

/** Map a unit random number to an angle on the valve's open arc. */
export function sampleOpenArc(v: ValveSpec, u: number): number {
  const { center, halfWidth } = openArc(v);
  return center + (2 * u - 1) * halfWidth;
}
