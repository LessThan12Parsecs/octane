/**
 * Minimal colour science for turning emission SPECTRA into display colours.
 *
 * - CIE 1931 2° colour-matching functions: multi-lobe Gaussian fit of
 *   Wyman, Sloan & Shirley, "Simple Analytic Approximations to the CIE XYZ
 *   Color Matching Functions", JCGT 2(2), 2013 (eq. 2; ≲1 % of peak error).
 * - XYZ → linear sRGB (D65): IEC 61966-2-1:1999 matrix.
 * - Planck's law.
 *
 * Photometric convention: XYZ values returned here are in cd/m² (i.e. already
 * multiplied by K_m = 683 lm/W) when they describe a radiance.
 * Pure TypeScript, no three.js — unit tested.
 */
import { K_M } from '../constants';

const H_PLANCK = 6.62607015e-34; // J s
const C_LIGHT = 299792458; // m/s
const K_B = 1.380649e-23; // J/K

/** Piecewise Gaussian g(λ; μ, σ1, σ2) of Wyman et al. (2013). */
function g(lambda: number, mu: number, s1: number, s2: number): number {
  const t = (lambda - mu) / (lambda < mu ? s1 : s2);
  return Math.exp(-0.5 * t * t);
}

/** CIE 1931 x̄(λ), λ in nm. */
export function cmfX(l: number): number {
  return 1.056 * g(l, 599.8, 37.9, 31.0) + 0.362 * g(l, 442.0, 16.0, 26.7) - 0.065 * g(l, 501.1, 20.4, 26.2);
}
/** CIE 1931 ȳ(λ) = photopic V(λ), λ in nm. */
export function cmfY(l: number): number {
  return 0.821 * g(l, 568.8, 46.9, 40.5) + 0.286 * g(l, 530.9, 16.3, 31.1);
}
/** CIE 1931 z̄(λ), λ in nm. */
export function cmfZ(l: number): number {
  return 1.217 * g(l, 437.0, 11.8, 36.0) + 0.681 * g(l, 459.0, 26.0, 13.8);
}

export type Vec3 = [number, number, number];

/** XYZ → linear sRGB (D65 white). Negative (out-of-gamut) components are NOT clipped here. */
export function xyzToLinearSrgb(xyz: ArrayLike<number>, out: Vec3 = [0, 0, 0]): Vec3 {
  const X = xyz[0], Y = xyz[1], Z = xyz[2];
  out[0] = 3.2406 * X - 1.5372 * Y - 0.4986 * Z;
  out[1] = -0.9689 * X + 1.8758 * Y + 0.0415 * Z;
  out[2] = 0.0557 * X - 0.204 * Y + 1.057 * Z;
  return out;
}

/** Gamut-map by clipping negative components (keeps the hue of saturated violets/greens). */
export function clipToGamut(rgb: Vec3): Vec3 {
  if (rgb[0] < 0) rgb[0] = 0;
  if (rgb[1] < 0) rgb[1] = 0;
  if (rgb[2] < 0) rgb[2] = 0;
  return rgb;
}

/** CIE xy chromaticity of an XYZ triple. */
export function chromaticity(xyz: ArrayLike<number>): [number, number] {
  const s = xyz[0] + xyz[1] + xyz[2];
  return s > 0 ? [xyz[0] / s, xyz[1] / s] : [1 / 3, 1 / 3];
}

/** Planck spectral radiance B_λ(λ, T), W m⁻² sr⁻¹ m⁻¹ (λ in metres). */
export function planckRadiance(lambdaM: number, T: number): number {
  if (T <= 0) return 0;
  const x = (H_PLANCK * C_LIGHT) / (lambdaM * K_B * T);
  if (x > 700) return 0;
  return (2 * H_PLANCK * C_LIGHT * C_LIGHT) / (lambdaM ** 5 * Math.expm1(x));
}

// ---------------------------------------------------------------------------
// Spectra
// ---------------------------------------------------------------------------

/** One spectral component. `weight` = fraction of the emitter's radiant power (in 250–1000 nm). */
export type SpectrumComponent =
  | { kind: 'line'; nm: number; fwhm: number; weight: number }
  | { kind: 'gauss'; nm: number; sigma: number; weight: number }
  | { kind: 'blackbody'; T: number; weight: number };

const L_MIN = 250;
const L_MAX = 1000;
const L_STEP = 1;

function componentShape(c: SpectrumComponent, l: number): number {
  switch (c.kind) {
    case 'line': {
      const s = c.fwhm / 2.3548;
      const t = (l - c.nm) / s;
      return Math.exp(-0.5 * t * t);
    }
    case 'gauss': {
      const t = (l - c.nm) / c.sigma;
      return Math.exp(-0.5 * t * t);
    }
    case 'blackbody':
      return planckRadiance(l * 1e-9, c.T);
  }
}

/**
 * XYZ (cd/m²) produced by 1 W m⁻² sr⁻¹ of radiance with the given spectral
 * shape. Each component is normalised to unit power over 250–1000 nm and
 * weighted; so UV (e.g. OH* 309 nm) or IR power correctly yields ~no light.
 */
export function spectrumXYZPerWatt(components: readonly SpectrumComponent[]): Vec3 {
  const out: Vec3 = [0, 0, 0];
  let wsum = 0;
  for (const c of components) wsum += c.weight;
  if (wsum <= 0) return out;
  for (const c of components) {
    let p = 0, X = 0, Y = 0, Z = 0;
    for (let l = L_MIN; l <= L_MAX; l += L_STEP) {
      const s = componentShape(c, l);
      p += s;
      X += s * cmfX(l);
      Y += s * cmfY(l);
      Z += s * cmfZ(l);
    }
    if (p <= 0) continue;
    const k = (c.weight / wsum) * K_M / p;
    out[0] += X * k;
    out[1] += Y * k;
    out[2] += Z * k;
  }
  return out;
}

/**
 * Absolute XYZ (cd/m²) of a black body at T, i.e. ∫ B_λ(T) cmf(λ) dλ · K_m.
 * Multiply by a (small) emissivity for gray-body gas.
 */
export function blackbodyXYZ(T: number, out: Vec3 = [0, 0, 0]): Vec3 {
  let X = 0, Y = 0, Z = 0;
  const dl = 2; // nm
  for (let l = 360; l <= 830; l += dl) {
    const b = planckRadiance(l * 1e-9, T) * 1e-9 * dl; // W m⁻² sr⁻¹ in this bin
    X += b * cmfX(l);
    Y += b * cmfY(l);
    Z += b * cmfZ(l);
  }
  out[0] = X * K_M;
  out[1] = Y * K_M;
  out[2] = Z * K_M;
  return out;
}

/**
 * Cached black-body table for per-frame lookups, 300–4500 K in 10 K steps
 * (linear interpolation of log XYZ between entries).
 */
const BB_T0 = 300;
const BB_DT = 10;
const BB_N = 421;
let bbTable: Float64Array | null = null;

function ensureBBTable(): Float64Array {
  if (bbTable) return bbTable;
  const t = new Float64Array(BB_N * 3);
  const v: Vec3 = [0, 0, 0];
  for (let i = 0; i < BB_N; i++) {
    blackbodyXYZ(BB_T0 + i * BB_DT, v);
    t[3 * i] = Math.log(Math.max(v[0], 1e-300));
    t[3 * i + 1] = Math.log(Math.max(v[1], 1e-300));
    t[3 * i + 2] = Math.log(Math.max(v[2], 1e-300));
  }
  bbTable = t;
  return t;
}

/** Fast black-body XYZ (cd/m²) for 300 ≤ T ≤ 4500 K (exact routine outside). Allocation-free with `out`. */
export function blackbodyXYZFast(T: number, out: Vec3): Vec3 {
  if (!(T >= BB_T0) || T >= BB_T0 + (BB_N - 1) * BB_DT) {
    if (!(T > 0)) { out[0] = out[1] = out[2] = 0; return out; }
    return blackbodyXYZ(T, out);
  }
  const tab = ensureBBTable();
  const f = (T - BB_T0) / BB_DT;
  const i = Math.floor(f);
  const w = f - i;
  for (let k = 0; k < 3; k++) {
    out[k] = Math.exp(tab[3 * i + k] * (1 - w) + tab[3 * (i + 1) + k] * w);
  }
  return out;
}
