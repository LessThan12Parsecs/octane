/**
 * Emission spectra of the light sources inside the cylinder, and their colour
 * (linear sRGB, in cd/m² per W m⁻² sr⁻¹ of radiance).
 *
 * Band positions are standard spectroscopy (Pearse & Gaydon, "The
 * Identification of Molecular Spectra"; Gaydon, "The Spectroscopy of
 * Flames"). RELATIVE band strengths are rough and marked UNVERIFIED: they set
 * hue, not brightness.
 */
import { clipToGamut, spectrumXYZPerWatt, xyzToLinearSrgb, type SpectrumComponent, type Vec3 } from './cie';

/**
 * Premixed hydrocarbon flame-front chemiluminescence at equivalence ratio φ.
 * OH* A–X 306–310 nm (UV, invisible), CH* A–X 431.4 nm, C2* Swan Δv = +1/0/−1
 * at 473.7/516.5/563.5 nm, CO2* continuum (broad, ~350–600 nm).
 * The C2* : CH* ratio grows with φ (used as an optical φ sensor, e.g.
 * Kojima et al. 2005; Hardalupas & Orain 2004).
 * UNVERIFIED: C2* : CH* ≈ 0.6·exp(2.5(φ−1)), CO2*:OH*:CH* ≈ 1.5:1:0.5.
 */
export function flameFrontSpectrum(phi: number): SpectrumComponent[] {
  const p = Math.min(Math.max(phi, 0.5), 1.6);
  const ch = 0.5;
  const c2 = ch * 0.6 * Math.exp(2.5 * (p - 1));
  return [
    { kind: 'line', nm: 308.9, fwhm: 6, weight: 1.0 },
    { kind: 'line', nm: 431.4, fwhm: 6, weight: ch },
    { kind: 'line', nm: 516.5, fwhm: 6, weight: 0.6 * c2 },
    { kind: 'line', nm: 473.7, fwhm: 6, weight: 0.2 * c2 },
    { kind: 'line', nm: 563.5, fwhm: 6, weight: 0.2 * c2 },
    { kind: 'gauss', nm: 420, sigma: 90, weight: 1.5 },
  ];
}

/**
 * End-gas autoignition: CO2* continuum, formaldehyde HCHO* (Emeléus bands
 * 395–470 nm, cool-flame emission), CH*, OH*. UNVERIFIED relative weights.
 */
export const AUTOIGNITION_SPECTRUM: SpectrumComponent[] = [
  { kind: 'gauss', nm: 420, sigma: 90, weight: 1.5 },
  { kind: 'line', nm: 395.2, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 412.3, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 423.2, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 432.6, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 457.0, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 470.1, fwhm: 4, weight: 0.12 },
  { kind: 'line', nm: 431.4, fwhm: 6, weight: 0.3 },
  { kind: 'line', nm: 308.9, fwhm: 6, weight: 0.8 },
];

/** N2 second positive (C³Πu–B³Πg) and N2⁺ first negative (B–X) band heads — the violet of glow discharges in air. */
const N2_BANDS: SpectrumComponent[] = [
  { kind: 'line', nm: 315.9, fwhm: 3, weight: 0.25 },
  { kind: 'line', nm: 337.1, fwhm: 3, weight: 0.45 },
  { kind: 'line', nm: 357.7, fwhm: 3, weight: 0.3 },
  { kind: 'line', nm: 380.5, fwhm: 3, weight: 0.15 },
  { kind: 'line', nm: 399.8, fwhm: 3, weight: 0.06 },
  { kind: 'line', nm: 405.9, fwhm: 3, weight: 0.04 },
  { kind: 'line', nm: 391.4, fwhm: 3, weight: 0.12 },
  { kind: 'line', nm: 427.8, fwhm: 3, weight: 0.05 },
];

function scaled(cs: SpectrumComponent[], w: number): SpectrumComponent[] {
  const s = cs.reduce((a, c) => a + c.weight, 0);
  return cs.map((c) => ({ ...c, weight: (c.weight / s) * w }));
}

/** Glow discharge (low current, cathode fall ~300–500 V): N2 molecular bands only. UNVERIFIED weights. */
export const GLOW_SPECTRUM: SpectrumComponent[] = N2_BANDS;

/**
 * Arc (high current, thermal plasma ~6000 K): continuum + N II (500.5, 567.9,
 * 594.2 nm), Hα 656.3 nm (water in the charge), O I 777.2 nm, and N2 bands
 * from the cooler sheath. UNVERIFIED weights.
 */
export function arcSpectrum(T: number): SpectrumComponent[] {
  return [
    { kind: 'blackbody', T, weight: 0.75 },
    { kind: 'line', nm: 500.5, fwhm: 2, weight: 0.02 },
    { kind: 'line', nm: 567.9, fwhm: 2, weight: 0.02 },
    { kind: 'line', nm: 594.2, fwhm: 2, weight: 0.01 },
    { kind: 'line', nm: 656.3, fwhm: 2, weight: 0.02 },
    { kind: 'line', nm: 777.2, fwhm: 2, weight: 0.08 },
    ...scaled(N2_BANDS, 0.1),
  ];
}

/** Breakdown plasma: very hot continuum + N2 bands + N II lines → white-violet. UNVERIFIED weights. */
export function breakdownSpectrum(T: number): SpectrumComponent[] {
  return [
    { kind: 'blackbody', T, weight: 0.8 },
    { kind: 'line', nm: 399.5, fwhm: 2, weight: 0.02 },
    { kind: 'line', nm: 463.1, fwhm: 2, weight: 0.015 },
    { kind: 'line', nm: 500.5, fwhm: 2, weight: 0.015 },
    ...scaled(N2_BANDS, 0.15),
  ];
}

/** Linear sRGB (cd/m² per W m⁻² sr⁻¹) of a spectrum, gamut-clipped. */
export function spectrumRGBPerWatt(components: readonly SpectrumComponent[]): Vec3 {
  return clipToGamut(xyzToLinearSrgb(spectrumXYZPerWatt(components)));
}

/** Rec. 709 luminance of a linear sRGB triple. */
export function luminance(rgb: ArrayLike<number>): number {
  return 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2];
}

/** Cache of flame colours by φ (0.02 resolution). */
const flameCache = new Map<number, Vec3>();
export function flameRGBPerWatt(phi: number): Vec3 {
  const key = Math.round(Math.min(Math.max(phi, 0.5), 1.6) * 50);
  let v = flameCache.get(key);
  if (!v) {
    v = spectrumRGBPerWatt(flameFrontSpectrum(key / 50));
    flameCache.set(key, v);
  }
  return v;
}

let aiCache: Vec3 | null = null;
export function autoignitionRGBPerWatt(): Vec3 {
  return (aiCache ??= spectrumRGBPerWatt(AUTOIGNITION_SPECTRUM));
}
let glowCache: Vec3 | null = null;
export function glowRGBPerWatt(): Vec3 {
  return (glowCache ??= spectrumRGBPerWatt(GLOW_SPECTRUM));
}
const arcCache = new Map<number, Vec3>();
export function arcRGBPerWatt(T: number): Vec3 {
  let v = arcCache.get(T);
  if (!v) arcCache.set(T, (v = spectrumRGBPerWatt(arcSpectrum(T))));
  return v;
}
const bdCache = new Map<number, Vec3>();
export function breakdownRGBPerWatt(T: number): Vec3 {
  let v = bdCache.get(T);
  if (!v) bdCache.set(T, (v = spectrumRGBPerWatt(breakdownSpectrum(T))));
  return v;
}

/**
 * Equivalence ratio from burned-gas mole fractions by element balance:
 * φ = O atoms needed for complete oxidation (2C + H/2) / O atoms present.
 * Returns NaN if the composition has no products (e.g. before combustion).
 */
export function equivalenceRatioFromBurned(c: {
  CO2: number; H2O: number; CO: number; O2: number; H2: number; OH: number; H: number; O: number; NO: number;
}): number {
  const C = c.CO2 + c.CO;
  const H = 2 * c.H2O + 2 * c.H2 + c.OH + c.H;
  const O = 2 * c.CO2 + c.H2O + c.CO + 2 * c.O2 + c.OH + c.O + c.NO;
  if (!(C + H > 1e-6) || !(O > 0)) return NaN;
  return (2 * C + 0.5 * H) / O;
}
