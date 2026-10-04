/**
 * Knock: cylinder acoustic modes excited by end-gas autoignition, and the end-gas burn-up.
 *
 * ACOUSTIC MODES (Draper 1935, NACA Report 493, "The physical effects of detonation in a closed
 * cylindrical chamber"; Draper 1938, J. Aeronaut. Sci. 5:219): the transverse modes of a
 * rigid-walled cylinder of bore B are p' ∝ J_m(α_mn r/R) {cos mθ, sin mθ} with J'_m(α_mn) = 0,
 * f_mn = α_mn c / (π B)  (ω = 2α c / B). m = circumferential, n = radial order. Axial modes
 * (f ≥ c/(2h) ≳ 25 kHz near TDC) are not modelled.
 *
 * EXCITATION (first principles): linearising the energy and momentum equations of an inviscid
 * ideal gas with a volumetric heat-release rate q (W/m³),
 *   ∂p'/∂t + γ p ∇·u' = (γ − 1) q',   ρ ∂u'/∂t = −∇p'
 *   ⇒ ∂²p'/∂t² − c² ∇²p' = (γ − 1) ∂q'/∂t,
 * the classical inhomogeneous wave equation of thermoacoustics (e.g. Dowling & Stow 2003,
 * J. Propulsion Power 19:751, eq. 2; Culick 2006, RTO AG-AVT-039 §3). Expanding p' = Σ η_j ψ_j
 * over the Neumann eigenfunctions (∫ψ_j dV = 0 for j ≠ 0) and projecting (Galerkin):
 *   η̈_j + 2δ_j η̇_j + ω_j² η_j = (γ − 1) (⟨ψ_j⟩_eg / (V ⟨ψ_j²⟩_cyl)) dQ̇/dt,
 * where Q̇ (W) is the end-gas heat-release rate distributed uniformly over the end-gas region
 * (⟨·⟩_eg = mean over that region of the planform, ⟨·⟩_cyl = over the bore cross-section).
 * The j = 0 (uniform) mode reproduces the 0-D mean pressure rise dp/dt = (γ − 1)Q̇/V, which the
 * cycle model already contains, so only j ≥ 1 is synthesised here. For an instantaneous release
 * of Q the initial field is p'_0 = (γ−1)(Q/V)(1_eg/v_eg − 1) and the modal amplitudes are its
 * exact projections — the acoustic energy is ∫ p'² /(2ρc²) dV (Parseval, tested). A finite
 * release time reduces the excitation of each mode by the Fourier content of Q̇(t) at ω_j,
 * automatically (no tuning constant). The end gas is taken as the planform region outside the
 * flame circle (flame front assumed to span the clearance height — disc chamber).
 *
 * DAMPING: δ_j = (1/τ_d) (α_j/α_10)^½; τ_d = amplitude e-folding time of the (1,0) mode.
 * UNVERIFIED: τ_d = 1.0 ms default (knock oscillations are observed to decay within ~1–2 ms;
 * no fetched source with a number for the CFR; Di Gaeta et al., Fuel 104 (2013) 230 and Li &
 * Zhu, Int. J. Engine Res. (2019) doi:10.1177/1468087419869161 calibrate a damping term in the
 * wave equation but their values were not usable). The ∝ √ω scaling is that of wall
 * boundary-layer (Kirchhoff) losses — an assumption for the total damping. Our own estimate of
 * the purely laminar boundary-layer decay time is ≈ 5–10 ms, so the observed decay must be
 * dominated by other losses (turbulence, crevices, piston motion).
 * The Bessel-derivative zeros are checked against scipy (knock.test.ts) and against Li & Zhu
 * (2019) Table 1 (1.8412, 3.8317, 5.3314, 7.0156 …). UNVERIFIED (from memory, not fetched):
 * the volume/page of Draper (1938) and of Dowling & Stow (2003); the Culick (2006) section.
 * VERIFIED (review 2026-09-30, NTRS 19930091567 PDF): Draper (1935) NACA Report 493 eqs. (5)–(9):
 * ∂φ/∂r = 0 at r = a ⇒ J'_s(βa) = 0, β² = 4π²n²/c² − p² (p = 0 for transverse modes, so
 * n = (βa) c/(2πa) = α c/(π B)), with βa = 1.841 (s = 1), 3.054 (s = 2), 3.832 (s = 0), 5.332,
 * 7.016 …; its Fig. 2 applies them to the Waukesha CFR chamber.
 *
 * NON-CYLINDRICAL CHAMBERS (side-valve 'l-head'): per mode the integrator needs only the wavenumber, the
 * volume mean ⟨ψ_j²⟩, ψ_j at the sensor and the end-gas projections ⟨ψ_j⟩_eg, so KnockOscillator is a
 * shape-agnostic modal integrator over an AcousticModeSet (DiscModeSet: the Bessel modes above;
 * GridModeSet: numerical modes of a chamber made of flat-roofed vertical gas columns) and an
 * EndGasProjector (DiscEndGasProjector: the exact planform-circle quadrature the CFR is calibrated with;
 * GridEndGasProjector: the depth-weighted gas outside the flame ball about the plug). createKnockOscillator
 * builds the pair from an EngineSpec (geometry.chamber).
 * Long-wave acoustics of a chamber of local height H(x, z): with the pressure uniform over the height
 * (vertical equilibration; the first vertical resonance c/(2H) ≈ 19 kHz of the Model T bore column at TDC
 * lies above the modes of interest), integrating the linearised mass and momentum balances over the
 * height gives ∂p'/∂t = −(ρc²/H) ∇·(H ū), ρ ∂ū/∂t = −∇p' — the depth-weighted Helmholtz problem
 * ∇·(H∇ψ) + k²Hψ = 0 in the planform with H ∂ψ/∂n = 0 on its outline and ψ, H ∂ψ/∂n continuous across a
 * depth step (the long-wave equation of a basin of variable depth, Lamb 1932, Hydrodynamics ch. VIII;
 * Webster's horn equation in 1-D — UNVERIFIED article numbers, from memory). Its modes are orthogonal in
 * ∫Hψ_iψ_j dA = V⟨ψ_iψ_j⟩, so the Galerkin projection above holds with volume means (⟨ψ_j⟩_eg =
 * (1/V_eg)∫H_eg ψ_j dA, H_eg the end-gas length of each column). With a constant depth it is the disc
 * problem exactly (tested: α_mn within 0.2 %, GridModeSet). The approximation degrades where a depth step
 * meets a short wavelength: the Model T modes above ≈ 15 kHz are qualitative. Not modelled: the
 * inertance of the window between the bore column and the pocket (at TDC the crown stands 7.9 mm above
 * the deck, leaving a ≈ 5 mm opening between 25.4 and 12.9 mm tall columns) — an O(kH) end correction that
 * would lower the bore↔pocket modes by an estimated few per cent (UNVERIFIED estimate: slit end
 * correction ≈ 2–5 mm per side on a ≈ 140 mm half-wavelength).
 * Model T (engines/model-t.ts, at TDC): fundamental α = 1.20, f = α c/(πB) = 3.8 kHz at c = 950 m/s (the
 * bore's own (1,0) mode would ring at 5.8 kHz), a bore↔pocket sloshing mode largest at the plug; 26 modes
 * up to α 7.1. Grid modes are damped by the same δ ∝ √ω law with α_ref = the bore's α_10 (the decay rate
 * is the same function of frequency as the CFR's); the side-valve chamber's larger surface/volume ratio
 * would raise boundary-layer losses (UNVERIFIED, not modelled).
 *
 * END-GAS BURN-UP after the Livengood–Wu integral reaches 1: thermally stratified end gas
 * autoignites sequentially; a spread ΔT of temperatures maps to a spread of ignition times
 * Δt = τ |∂lnτ/∂T| ΔT (Zeldovich 1980 "spontaneous" propagation, u_a = (∂τ/∂x)⁻¹; Bradley &
 * Kalghatgi, Combust. Flame 156 (2009) 2307), plus the chemical excitation (heat-release)
 * time τ_e (Lutz, Kee, Miller, Dwyer & Oppenheim, Proc. Combust. Inst. 22 (1988) 1683),
 * which we computed for PRF/air with the LLNL gasoline-surrogate 2011 mechanism (see
 * EXCITATION_TIME_PRF).
 * In the NTC region ∂τ/∂T → 0 and the whole end gas ignites within ~τ_e (violent knock).
 * UNVERIFIED (from memory, not fetched): the journal/volume/page details of Zeldovich (1980),
 * Bradley & Kalghatgi (2009) and Lutz et al. (1988).
 */

import { lensArea } from '../combustion/flame-geometry';
import type { EngineSpec, LHeadChamberSpec } from '../core/engine-spec';

// ---------------------------------------------------------------------------------------
// Bessel functions
// ---------------------------------------------------------------------------------------

/**
 * Bessel function of the first kind J_n(x), integer n ≥ 0, x ≥ 0 — Miller's backward
 * recurrence normalised by 1 = J0 + 2ΣJ_2k (Abramowitz & Stegun 1964, 9.1.27 and 9.1.46;
 * Press et al., Numerical Recipes, §6.5). Relative accuracy ~1e-14 for x ≲ 200.
 * Not a hot-path routine (O(x + n) work).
 */
export function besselJ(n: number, x: number): number {
  if (x === 0) return n === 0 ? 1 : 0;
  const ax = Math.abs(x);
  const nMax = Math.max(n, ax);
  let top = Math.ceil(nMax + 30 + 3 * Math.sqrt(nMax));
  if (top % 2 === 1) top++;
  let jp = 0; // J_{k+1}
  let j = 1e-300; // J_k
  let sum = 0;
  let res = 0;
  const twoOverX = 2 / ax;
  for (let k = top; k > 0; k--) {
    const jm = k * twoOverX * j - jp; // J_{k-1}
    jp = j;
    j = jm;
    if (Math.abs(j) > 1e250) {
      j *= 1e-250;
      jp *= 1e-250;
      res *= 1e-250;
      sum *= 1e-250;
    }
    if (k - 1 === n) res = j;
    if (k - 1 > 0 && (k - 1) % 2 === 0) sum += 2 * j;
  }
  sum += j; // J0
  const v = res / sum;
  return x < 0 && n % 2 === 1 ? -v : v;
}

/** Derivative J'_n(x) = (J_{n−1}(x) − J_{n+1}(x))/2, J'_0 = −J_1. */
export function besselJPrime(n: number, x: number): number {
  if (n === 0) return -besselJ(1, x);
  return 0.5 * (besselJ(n - 1, x) - besselJ(n + 1, x));
}

const zeroCache = new Map<number, { zeros: number[]; scannedTo: number }>();

/**
 * Positive zeros of J'_m in ascending order up to xMax (x = 0 excluded, so for m = 0 the
 * first entry is 3.8317). Bracketing scan (step 0.05) + bisection to 1e-14. Cached per m.
 */
export function besselJPrimeZeros(m: number, xMax: number): number[] {
  let entry = zeroCache.get(m);
  if (!entry) {
    entry = { zeros: [], scannedTo: 0.25 };
    zeroCache.set(m, entry);
  }
  const dx = 0.05;
  let a = entry.scannedTo;
  let fa = besselJPrime(m, a);
  while (a < xMax) {
    const b = a + dx;
    const fb = besselJPrime(m, b);
    if (fa * fb < 0) {
      let lo = a;
      let hi = b;
      let flo = fa;
      for (let i = 0; i < 100 && hi - lo > 1e-14 * hi; i++) {
        const mid = 0.5 * (lo + hi);
        const fm = besselJPrime(m, mid);
        if (flo * fm <= 0) hi = mid;
        else {
          lo = mid;
          flo = fm;
        }
      }
      entry.zeros.push(0.5 * (lo + hi));
    }
    a = b;
    fa = fb;
  }
  entry.scannedTo = Math.max(entry.scannedTo, a);
  return entry.zeros.filter((v) => v <= xMax);
}

/**
 * α_mn: the (n+1)-th positive zero of J'_m for m ≥ 1 (α_10 = 1.8412, α_20 = 3.0542,
 * α_30 = 4.2012, α_11 = 5.3314), and the n-th positive zero of J'_0 for m = 0, n ≥ 1
 * (α_01 = 3.8317). Draper's (1935, 1938) notation (βa, s = m).
 */
export function modeZero(m: number, n: number): number {
  if (m === 0 && n < 1) throw new RangeError('mode (0,0) is the uniform (bulk) mode');
  const idx = m === 0 ? n - 1 : n;
  let xMax = 10;
  for (;;) {
    const z = besselJPrimeZeros(m, xMax);
    if (z.length > idx) return z[idx];
    xMax *= 2;
  }
}

/**
 * Frequency of the (m, n) transverse acoustic mode of a cylinder, Hz (Draper 1935 eq. 7–9, 1938):
 * f = α_mn c / (π B). c: sound speed, m/s; bore: m.
 */
export function cylinderModeFrequency(m: number, n: number, c: number, bore: number): number {
  return (modeZero(m, n) * c) / (Math.PI * bore);
}

// ---------------------------------------------------------------------------------------
// Knock oscillator
// ---------------------------------------------------------------------------------------

/**
 * Effective (long-wave) sound speed of the two-zone charge, m/s — Wood's (1930, A Textbook of
 * Sound) mixture rule: volume-averaged compressibility and density,
 *   1/(ρ̄ c̄²) = Σ φ_i/(γ_i p),  ρ̄ = Σ φ_i ρ_i  ⇒  c̄² = p / (ρ̄ Σ φ_i/γ_i),
 * with zone volume fractions φ_i (both zones at the common pressure p, ρ_i c_i² = γ_i p).
 * An approximation for modes whose wavelength (≈ πB/α) is comparable to the zones; late in
 * combustion (where knock occurs) the burned zone dominates and c̄ → c_burned.
 * Vb, Vu: zone volumes (m³); rhoB, rhoU: densities (kg/m³); gammaB, gammaU: cp/cv; p in Pa.
 */
export function twoZoneSoundSpeed(
  p: number,
  Vb: number,
  rhoB: number,
  gammaB: number,
  Vu: number,
  rhoU: number,
  gammaU: number,
): number {
  const V = Vb + Vu;
  const fb = Vb / V;
  const fu = Vu / V;
  const rho = fb * rhoB + fu * rhoU;
  return Math.sqrt(p / (rho * (fb / gammaB + fu / gammaU)));
}

/** Options for KnockOscillator. */
export interface KnockOscillatorOptions {
  /**
   * Disc modes: include every mode with α_mn ≤ maxAlpha (default 7.1: (1,0) (2,0) (0,1) (3,0) (4,0) (1,1)
   * (5,0) (2,1) (0,2)). Ignored when `modes` is given.
   */
  maxAlpha?: number;
  /** Amplitude e-folding time of the (1,0) mode, s (default 1e-3, UNVERIFIED — see file header). */
  decayTime?: number;
  /** Pressure-transducer position (x, z) on the head face, cylinder frame, m (default (0.9 R, 0)). */
  sensor?: readonly [number, number];
  /**
   * Measurement band [f_lo, f_hi] (Hz) of {@link KnockOscillator.peakSensorPressure} (MAPO): each mode
   * is weighted by the zero-phase Butterworth band-pass magnitude at its current frequency (see
   * bandGain). Default: none (all modes, gain 1). The reported sensorPressure() is always unfiltered.
   */
  band?: readonly [number, number];
  /**
   * Butterworth order of each band edge (default Infinity = ideal band-pass: a mode is in or out
   * according to its frequency at the reference sound speed of {@link KnockOscillator.setBandReference};
   * the measurement filters are not specified in the sources — UNVERIFIED).
   */
  bandOrder?: number;
  /**
   * Disc source-projection quadrature: [Gauss–Legendre nodes per radial sub-interval, unused] (default
   * [24, 0]; the angular integrals are exact — see DiscEndGasProjector). Round 1: midpoint cells
   * (radial × angular, 48 × 96).
   */
  quadrature?: readonly [number, number];
  /**
   * Acoustic modes of the chamber (default: the Bessel modes of the bore disc, DiscModeSet(bore,
   * maxAlpha)). Its reference length must equal `bore`. Shareable between cylinders (immutable).
   */
  modes?: AcousticModeSet;
  /**
   * End-gas source projector (default: the exact disc projector of DiscModeSet modes; required with
   * other mode sets). Shareable between cylinders.
   */
  projector?: EndGasProjector;
  /**
   * Spark-gap centre (x, y, z), cylinder frame, m: the flame centre of
   * {@link KnockOscillator.setEndGasVolumeFraction} for the disc projector (default the bore axis; the
   * grid projector has its own, fixed at construction).
   */
  spark?: readonly [number, number, number];
}

/** Default (1,0) amplitude decay time, s. UNVERIFIED (see file header). */
export const KNOCK_DECAY_TIME = 1.0e-3;

// ---------------------------------------------------------------------------------------
// Acoustic mode sets and end-gas source projectors
// ---------------------------------------------------------------------------------------

/**
 * Transverse acoustic modes ψ_j of a chamber (the uniform j = 0 mode excluded): everything the modal
 * integrator (KnockOscillator) needs to know about the chamber shape besides the source projection —
 * ω_j = 2 α_j c / L, the volume-mean ⟨ψ_j²⟩ and ψ_j at any planform point. Immutable: one instance can
 * serve every cylinder of an engine.
 */
export interface AcousticModeSet {
  readonly nModes: number;
  /** Reference length L, m (the bore): ω_j = 2 α_j c / L, f_j = α_j c / (π L). */
  readonly bore: number;
  /** Dimensionless eigen-wavenumber α_j = k_j L / 2 (k_j: wavenumber, 1/m), ascending. */
  readonly alpha: Float64Array;
  /** Volume-mean square ⟨ψ_j²⟩ of the (unnormalised) mode shape over the chamber. */
  readonly meanSquare: Float64Array;
  /** Relative damping (α_j/α_ref)^½, δ_j = dampScale_j / τ_d (file header). */
  readonly dampScale: Float64Array;
  /** Circumferential order m of a Bessel mode (−1: numerical mode). */
  readonly m: Int32Array;
  /** Radial order n of a Bessel mode (numerical mode: its index j). */
  readonly n: Int32Array;
  /** 1 = cos(mθ) orientation, 0 = sin(mθ) (numerical modes: 1). */
  readonly cosine: Uint8Array;
  /** Mode shape ψ_j at planform point (x, z) (cylinder frame, m); 0 outside the chamber. */
  modeShape(j: number, x: number, z: number): number;
}

/** Transverse Bessel modes of a rigid cylinder of bore B (Draper 1935, 1938; file header): the flat-disc chamber. */
export class DiscModeSet implements AcousticModeSet {
  readonly bore: number;
  readonly radius: number;
  readonly nModes: number;
  readonly m: Int32Array;
  readonly n: Int32Array;
  readonly alpha: Float64Array;
  readonly cosine: Uint8Array;
  /** Mean-square mode shape over the bore cross-section ⟨ψ_j²⟩, ψ_j = J_m(α_j r/R)·{cos, sin}(mθ) (unnormalised). */
  readonly meanSquare: Float64Array;
  /** (α_j/α_10)^½. */
  readonly dampScale: Float64Array;

  /** Every mode with α_mn ≤ maxAlpha (default 7.1), sorted by frequency, cos before sin. */
  constructor(bore: number, maxAlpha = 7.1) {
    this.bore = bore;
    this.radius = bore / 2;
    const ms: number[] = [];
    const ns: number[] = [];
    const as: number[] = [];
    const cs: number[] = [];
    // the first zero of J'_m exceeds m, so m ≤ maxAlpha bounds the search
    for (let m = 0; m <= maxAlpha; m++) {
      const z = besselJPrimeZeros(m, maxAlpha);
      z.forEach((a, i) => {
        const n = m === 0 ? i + 1 : i;
        ms.push(m);
        ns.push(n);
        as.push(a);
        cs.push(1);
        if (m > 0) {
          ms.push(m);
          ns.push(n);
          as.push(a);
          cs.push(0);
        }
      });
    }
    // sort by frequency (stable: cos before sin)
    const order = as.map((_, i) => i).sort((i, j) => as[i] - as[j] || cs[j] - cs[i]);
    const N = order.length;
    this.nModes = N;
    this.m = Int32Array.from(order.map((i) => ms[i]));
    this.n = Int32Array.from(order.map((i) => ns[i]));
    this.alpha = Float64Array.from(order.map((i) => as[i]));
    this.cosine = Uint8Array.from(order.map((i) => cs[i]));
    this.meanSquare = new Float64Array(N);
    this.dampScale = new Float64Array(N);
    const a10 = modeZero(1, 0);
    for (let j = 0; j < N; j++) {
      const a = this.alpha[j];
      const m = this.m[j];
      // ⟨ψ²⟩ over the disc: (1/πR²) ∫ J_m(αr/R)² r dr ∫ trig² dθ
      //   = (1 − m²/α²) J_m(α)² × (1 for m = 0, ½ for m ≥ 1)
      const jm = besselJ(m, a);
      this.meanSquare[j] = (1 - (m * m) / (a * a)) * jm * jm * (m === 0 ? 1 : 0.5);
      this.dampScale[j] = Math.sqrt(a / a10);
    }
  }

  /** Mode shape ψ_j at planform point (x, z) (cylinder frame, m); 0 outside the bore. */
  modeShape(j: number, x: number, z: number): number {
    const r = Math.hypot(x, z);
    if (r > this.radius * (1 + 1e-12)) return 0;
    const th = Math.atan2(z, x);
    const m = this.m[j];
    const radial = besselJ(m, (this.alpha[j] * r) / this.radius);
    return radial * (this.cosine[j] ? Math.cos(m * th) : Math.sin(m * th));
  }
}

/**
 * Projection of a distributed end-gas heat release on the modes of an AcousticModeSet (the source
 * projector). The end gas is the part of the source domain farther than a flame radius r from the flame
 * centre; its "measure" is an area (disc projector: the planform, flame front spanning the clearance
 * height) or a volume (grid projector: depth-weighted columns). The modal source shape of a region is
 * S_j = ⟨ψ_j⟩_eg / ⟨ψ_j²⟩ (file header) with ⟨ψ_j⟩_eg = outsideIntegrals / outsideMeasure.
 */
export interface EndGasProjector {
  /** Measure of the whole source domain (m² or m³). */
  readonly totalMeasure: number;
  /** Flame radius beyond which nothing is outside, m. */
  readonly maxRadius: number;
  /** Measure of the end-gas region outside flame radius r (m² or m³). */
  outsideMeasure(r: number): number;
  /** acc[j] = ∫ over the region outside r of ψ_j (measure-weighted, unnormalised); returns the region's measure. */
  outsideIntegrals(r: number, acc: Float64Array): number;
  /** Flame radius in [lo, hi] whose outside measure is `target` (outsideMeasure(lo) ≥ target ≥ outsideMeasure(hi)). */
  radiusForOutsideMeasure(target: number, lo: number, hi: number): number;
  /** Flame radius whose outside measure is the fraction f (0..1) of totalMeasure. */
  radiusForOutsideFraction(f: number): number;
}

/**
 * Exact end-gas projector of the flat-disc chamber: the end gas is the planform outside a circle of
 * radius r about the flame centre (x0, z0) projected on the head face (setCentre), the flame front
 * spanning the clearance height.
 *
 * Quadrature (validation round 2): the angular integrals over the part of each ring r outside the
 * circle are EXACT — with d = |c|, θ_c = atan2(z0, x0), the ring is outside the circle for
 * |θ − θ_c| > α(r), cos α = (r² + d² − r_f²)/(2 r d), so ∫cos mθ dθ = −2 cos(mθ_c) sin(mα)/m,
 * ∫sin mθ dθ = −2 sin(mθ_c) sin(mα)/m, ∫dθ = 2π − 2α — and the radial integrals use Gauss–Legendre
 * on the sub-intervals between the breakpoints |d − r_f| and d + r_f, with r = a + (b − a)(3u² − 2u³)
 * on the partial-ring intervals (the arc length grows like √(r − a) at a breakpoint; the substitution
 * makes the integrand smooth). Round 1 used a 48 × 96 midpoint grid: an end-gas crescent thinner
 * than the outermost cell row (≈ 0.43 mm) contained no midpoint, so a detected knock rang with
 * MAPO = 0 exactly, and MAPO was biased by 1–4 % elsewhere. `quad[0]` (default 24) is the
 * number of Gauss–Legendre nodes per radial sub-interval.
 */
export class DiscEndGasProjector implements EndGasProjector {
  readonly modes: DiscModeSet;
  readonly radius: number;
  /** Bore cross-section πR², m². */
  readonly totalMeasure: number;
  /** [Gauss–Legendre nodes per radial sub-interval, unused] (KnockOscillatorOptions.quadrature). */
  quad: readonly [number, number];
  private x0 = 0;
  private z0 = 0;
  private d = 0;
  private readonly cosMScratch: Float64Array;
  private readonly sinMScratch: Float64Array;
  private readonly bpScratch = new Float64Array(4);

  constructor(modes: DiscModeSet, quadrature: readonly [number, number] = [24, 0]) {
    this.modes = modes;
    this.radius = modes.radius;
    this.totalMeasure = Math.PI * this.radius * this.radius;
    this.quad = quadrature;
    this.cosMScratch = new Float64Array(modes.nModes);
    this.sinMScratch = new Float64Array(modes.nModes);
  }

  /** Flame centre projected on the head face (x0, z0), cylinder frame, m. */
  setCentre(x0: number, z0: number): void {
    this.x0 = x0;
    this.z0 = z0;
    this.d = Math.hypot(x0, z0);
  }

  get maxRadius(): number {
    return this.d + this.radius;
  }

  outsideMeasure(r: number): number {
    return this.totalMeasure - circleLensArea(r, this.radius, this.d);
  }

  /** Bisection (60 halvings) on the closed-form lens area. */
  radiusForOutsideMeasure(target: number, lo: number, hi: number): number {
    for (let it = 0; it < 60; it++) {
      const mid = 0.5 * (lo + hi);
      if (this.outsideMeasure(mid) > target) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /**
   * Radius of the planform circle about the centre whose complement in the bore has the area fraction
   * f: πR² − lens(r) = f πR² (closed-form circle–circle lens of combustion/flame-geometry, bisection) —
   * the same operations as cycle-model.ts endGasCircleRadius (the CFR acoustic source region).
   */
  radiusForOutsideFraction(f: number): number {
    const R = this.radius;
    const d = Math.max(this.d, 1e-9);
    const target = (1 - f) * Math.PI * R * R; // lens area inside the circle
    let lo = 0;
    let hi = R + d;
    for (let i = 0; i < 60; i++) {
      const mid = 0.5 * (lo + hi);
      if (lensArea(mid, R, d) < target) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /**
   * ∫ over the planform outside the circle (radius rf about (x0, z0)) of J_m(α r/R)·trig(mθ) per mode
   * into acc (area-weighted, unnormalised); returns the region's area, m². Quadrature: see the class.
   */
  outsideIntegrals(flameRadius: number, acc: Float64Array): number {
    const md = this.modes;
    const R = this.radius;
    const N = md.nModes;
    const x0 = this.x0;
    const z0 = this.z0;
    acc.fill(0);
    const gl = gaussLegendre(Math.max(4, Math.round(this.quad[0])));
    const rf = Math.max(0, flameRadius);
    const d = Math.hypot(x0, z0);
    const thc = Math.atan2(z0, x0);
    const cm = this.cosMScratch;
    const sm = this.sinMScratch;
    for (let j = 0; j < N; j++) {
      cm[j] = Math.cos(md.m[j] * thc);
      sm[j] = Math.sin(md.m[j] * thc);
    }
    // radial sub-intervals of [0, R] between the breakpoints
    const bp = this.bpScratch;
    let nb = 0;
    bp[nb++] = 0;
    for (const b of [Math.abs(d - rf), d + rf]) if (b > 0 && b < R) bp[nb++] = b;
    bp[nb++] = R;
    if (nb === 4 && bp[1] > bp[2]) {
      const t = bp[1];
      bp[1] = bp[2];
      bp[2] = t;
    }
    let aEg = 0;
    for (let iv = 0; iv + 1 < nb; iv++) {
      const a = bp[iv];
      const b = bp[iv + 1];
      if (!(b > a)) continue;
      // classify the interval at its midpoint: 1 = whole ring outside, 0 = inside, 0.5 = partial
      const cMid = ringCos(0.5 * (a + b), d, rf);
      const partial = cMid > -1 && cMid < 1;
      if (!partial && cMid <= -1) continue; // whole ring inside the flame circle
      for (let q = 0; q < gl.x.length; q++) {
        const u = 0.5 * (gl.x[q] + 1);
        let r: number;
        let jac: number;
        if (partial) {
          r = a + (b - a) * u * u * (3 - 2 * u);
          jac = 6 * u * (1 - u) * (b - a) * 0.5 * gl.w[q];
        } else {
          r = a + (b - a) * u;
          jac = (b - a) * 0.5 * gl.w[q];
        }
        let alpha = 0; // half-angle of the arc INSIDE the circle
        if (partial) {
          const c = ringCos(r, d, rf);
          alpha = c >= 1 ? 0 : c <= -1 ? Math.PI : Math.acos(c);
        }
        const wr = jac * r;
        aEg += wr * (2 * Math.PI - 2 * alpha);
        for (let j = 0; j < N; j++) {
          const m = md.m[j];
          let ang: number;
          if (m === 0) ang = 2 * Math.PI - 2 * alpha;
          else if (!partial) ang = 0;
          else ang = (-2 * (md.cosine[j] ? cm[j] : sm[j]) * Math.sin(m * alpha)) / m;
          if (ang === 0) continue;
          acc[j] += wr * ang * besselJ(m, (md.alpha[j] * r) / R);
        }
      }
    }
    return aEg > 0 ? aEg : 0;
  }
}

// ---------------------------------------------------------------------------------------
// Knock oscillator (shape-agnostic modal integrator)
// ---------------------------------------------------------------------------------------

/**
 * Damped transverse acoustic modes of the combustion chamber driven by the end-gas
 * heat-release rate (see file header). Allocation-free step(); one instance per cylinder.
 * The chamber enters only through its AcousticModeSet (options.modes; default the bore's Bessel
 * modes) and EndGasProjector (options.projector) — createKnockOscillator builds both from an EngineSpec.
 */
export class KnockOscillator {
  readonly bore: number;
  readonly radius: number;
  readonly nModes: number;
  /** Circumferential order m of mode j (−1: numerical mode). */
  readonly m: Int32Array;
  /** Radial order n of mode j (numerical mode: j). */
  readonly n: Int32Array;
  /** α_j (ω_j = 2 α_j c / B). */
  readonly alpha: Float64Array;
  /** 1 = cos(mθ) orientation, 0 = sin(mθ). */
  readonly cosine: Uint8Array;
  /** Volume-mean square mode shape ⟨ψ_j²⟩ (disc: over the bore cross-section; unnormalised ψ_j). */
  readonly meanSquare: Float64Array;
  /** Modal source shape ⟨ψ_j⟩_eg / ⟨ψ_j²⟩ for the current end-gas region (dimensionless). */
  readonly sourceShape: Float64Array;
  /** Mode shape at the sensor ψ_j(x_s, z_s). */
  readonly psiSensor: Float64Array;
  /** Modal amplitudes η_j (Pa) and rates η̇_j (Pa/s). */
  readonly eta: Float64Array;
  readonly etaDot: Float64Array;
  /** The chamber's modes and end-gas projector (immutable, shareable between cylinders). */
  readonly modes: AcousticModeSet;
  readonly projector: EndGasProjector;
  /** Spark-gap centre, cylinder frame, m (flame centre of setEndGasVolumeFraction with the disc projector). */
  readonly spark: readonly [number, number, number];
  /** Relative damping factor (α_j/α_ref)^½. */
  private readonly dampScale: Float64Array;
  private readonly decayTime: number;
  private readonly bandLo: number;
  private readonly bandHi: number;
  private readonly bandOrder: number;
  /** Band gains of the modes at the sound speed of the current step (see bandGain). */
  private readonly gainScratch: Float64Array;
  private qPrev = 0;
  private readonly projScratch: Float64Array;
  /** Modal forcing Q̇·S_j at the end of the previous step, and the step-start impulse scratch. */
  private readonly qsPrev: Float64Array;
  private readonly impScratch: Float64Array;
  /** Sequential-autoignition shells: projections (K × nModes, row k = shell k) and their count. */
  private shellShapes: Float64Array = new Float64Array(0);
  shellCount = 0;
  /** Shell whose projection is the current sourceShape (−1: the whole end-gas region). */
  currentShell = -1;
  /**
   * End-gas fraction of the source domain for the current source (0..1): of the planform area (disc
   * projector) or of the chamber volume (grid projector).
   */
  endGasAreaFraction = 0;

  constructor(bore: number, opts: KnockOscillatorOptions = {}) {
    this.bore = bore;
    this.radius = bore / 2;
    this.decayTime = opts.decayTime ?? KNOCK_DECAY_TIME;
    this.bandLo = opts.band ? opts.band[0] : 0;
    this.bandHi = opts.band ? opts.band[1] : Infinity;
    this.bandOrder = opts.bandOrder ?? Infinity;
    const modes = opts.modes ?? new DiscModeSet(bore, opts.maxAlpha);
    if (modes.bore !== bore) throw new RangeError(`KnockOscillator: mode-set reference length ${modes.bore} m ≠ bore ${bore} m`);
    let projector = opts.projector;
    if (!projector) {
      if (!(modes instanceof DiscModeSet)) throw new Error('KnockOscillator: a non-disc mode set needs an end-gas projector (options.projector)');
      projector = new DiscEndGasProjector(modes, opts.quadrature ?? [24, 0]);
    }
    this.modes = modes;
    this.projector = projector;
    this.spark = opts.spark ?? [0, 0, 0];
    const N = modes.nModes;
    this.nModes = N;
    this.m = modes.m;
    this.n = modes.n;
    this.alpha = modes.alpha;
    this.cosine = modes.cosine;
    this.meanSquare = modes.meanSquare;
    this.dampScale = modes.dampScale;
    this.sourceShape = new Float64Array(N);
    this.psiSensor = new Float64Array(N);
    this.eta = new Float64Array(N);
    this.etaDot = new Float64Array(N);
    this.projScratch = new Float64Array(N);
    this.qsPrev = new Float64Array(N);
    this.impScratch = new Float64Array(N);
    this.gainScratch = new Float64Array(N).fill(1);
    const s = opts.sensor ?? [0.9 * this.radius, 0];
    for (let j = 0; j < N; j++) this.psiSensor[j] = this.modeShape(j, s[0], s[1]);
  }

  /** Mode shape ψ_j at planform point (x, z) (cylinder frame, m); 0 outside the chamber. */
  modeShape(j: number, x: number, z: number): number {
    return this.modes.modeShape(j, x, z);
  }

  /**
   * Source-projection quadrature of the disc projector (KnockOscillatorOptions.quadrature; validation
   * tests refine it through this private accessor); undefined for other projectors.
   */
  private get quad(): readonly [number, number] | undefined {
    const p = this.projector;
    return p instanceof DiscEndGasProjector ? p.quad : undefined;
  }

  private set quad(q: readonly [number, number] | undefined) {
    const p = this.projector;
    if (q && p instanceof DiscEndGasProjector) p.quad = q;
  }

  /**
   * Measurement-band weight of mode j at sound speed c: zero-phase Butterworth magnitude
   * [1 + (f_lo/f)^{2n}]^{−½} [1 + (f/f_hi)^{2n}]^{−½} at f = f_j(c) (1 without a band; n = ∞: 0 or 1
   * by the frequency at the reference sound speed, setBandReference). Validation
   * round 2: MAPO had summed all 16 modes up to ≈ 25 kHz while every measured MAPO is band-passed
   * (4–18 kHz Hoth/SON/critical CR; 6–20 kHz Rockstroh) — the out-of-band modes added 30–55 %.
   */
  bandGain(j: number, c: number): number {
    if (this.bandLo <= 0 && this.bandHi === Infinity) return 1;
    if (this.bandOrder === Infinity) {
      // ideal band-pass: mode in/out at the reference sound speed (fixed per knock event, so the
      // tracked signal has no switching discontinuity when c drifts)
      const fr = this.modeFrequency(j, this.bandRefC > 0 ? this.bandRefC : c);
      return fr >= this.bandLo && fr <= this.bandHi ? 1 : 0;
    }
    const f = this.modeFrequency(j, c);
    const n2 = 2 * this.bandOrder;
    const lo = this.bandLo > 0 ? Math.pow(this.bandLo / f, n2) : 0;
    const hi = this.bandHi < Infinity ? Math.pow(f / this.bandHi, n2) : 0;
    return 1 / Math.sqrt((1 + lo) * (1 + hi));
  }

  /** Sound speed at which the ideal band-pass decides which modes are in band (0: the step's own c). */
  private bandRefC = 0;

  /** Set the reference sound speed of the ideal band-pass (call at each knock onset), m/s. */
  setBandReference(c: number): void {
    this.bandRefC = c > 0 ? c : 0;
  }

  /** Mode frequency f_j = α_j c/(π B), Hz. */
  modeFrequency(j: number, c: number): number {
    return (this.alpha[j] * c) / (Math.PI * this.bore);
  }

  /** Zero all amplitudes and the source memory (new cycle). */
  reset(): void {
    this.eta.fill(0);
    this.etaDot.fill(0);
    this.qPrev = 0;
    this.qsPrev.fill(0);
    this.peakSensorPressure = 0;
    this.shellCount = 0;
    this.currentShell = -1;
  }

  /** The disc projector (planform-circle API); throws for other chambers. */
  private discProjector(api: string): DiscEndGasProjector {
    const p = this.projector;
    if (!(p instanceof DiscEndGasProjector)) throw new Error(`KnockOscillator.${api}: planform-circle end-gas regions need the disc projector; use setEndGasVolumeFraction`);
    return p;
  }

  /**
   * Set the end-gas region (disc chamber): the planform outside the flame circle of radius flameRadius
   * (m) centred at (x0, z0) (flame centre projected on the head face). Computes the modal projections
   * sourceShape_j = ⟨ψ_j⟩_eg/⟨ψ_j²⟩ (not a hot path; quadrature: DiscEndGasProjector) and returns the
   * end-gas area fraction.
   */
  setEndGasRegion(x0: number, z0: number, flameRadius: number): number {
    this.discProjector('setEndGasRegion').setCentre(x0, z0);
    return this.setEndGasOutsideRadius(flameRadius);
  }

  /**
   * Sequential autoignition source (validation round 2), disc chamber: the end-gas region outside the
   * flame circle (radius rInner about (x0, z0)) split into K equal-area shells ordered by distance from
   * the flame centre, shell 0 the farthest. See setSequentialEndGasOutsideRadius.
   */
  setSequentialEndGasRegion(x0: number, z0: number, rInner: number, K = 16): number {
    this.discProjector('setSequentialEndGasRegion').setCentre(x0, z0);
    return this.setSequentialEndGasOutsideRadius(rInner, K);
  }

  /**
   * Chamber-agnostic knock source (the cycle model's call at autoignition): the end gas is the fraction
   * f (0..1) of the source domain farthest from the spark — the planform outside a circle about the
   * plug with area fraction f (disc chamber, exactly the CFR source of cycle-model.ts onAutoignition)
   * or the gas outside a flame ball about the plug with volume fraction f (grid projector) — released
   * in K sequential shells (K = 0: uniformly). Returns the end-gas fraction of the source domain.
   * flameCentre (cylinder frame, m; default options.spark): the flame centre — the disc projector uses its
   * (x, z); a grid projector's is fixed at construction and a different one throws.
   */
  setEndGasVolumeFraction(f: number, K = 16, flameCentre?: readonly [number, number, number]): number {
    const P = this.projector;
    if (P instanceof DiscEndGasProjector) {
      const c = flameCentre ?? this.spark;
      P.setCentre(c[0], c[2]);
    } else if (flameCentre && P instanceof GridEndGasProjector) {
      const s = P.spark;
      if (flameCentre[0] !== s[0] || flameCentre[1] !== s[1] || flameCentre[2] !== s[2]) {
        throw new Error('KnockOscillator.setEndGasVolumeFraction: the grid projector was built for another spark position');
      }
    }
    const r = P.radiusForOutsideFraction(f);
    return K > 0 ? this.setSequentialEndGasOutsideRadius(r, K) : this.setEndGasOutsideRadius(r);
  }

  /**
   * Set the end-gas region outside flame radius r about the projector's flame centre (disc: planform
   * circle; grid: ball about the plug) — uniform release over the whole region. Returns the end-gas
   * fraction of the source domain.
   */
  setEndGasOutsideRadius(flameRadius: number): number {
    const N = this.nModes;
    const acc = this.projScratch;
    const aEg = this.projector.outsideIntegrals(flameRadius, acc);
    const aCyl = this.projector.totalMeasure;
    this.endGasAreaFraction = aEg / aCyl;
    for (let j = 0; j < N; j++) this.sourceShape[j] = aEg > 0 ? acc[j] / aEg / this.meanSquare[j] : 0;
    this.shellCount = 0;
    this.currentShell = -1;
    return this.endGasAreaFraction;
  }

  /**
   * Sequential autoignition source (validation round 2): the end-gas region outside flame radius rInner
   * split into K shells of equal measure ordered by distance from the flame centre, shell 0 the farthest
   * (the end gas at the chamber periphery opposite the plug, which the flame reaches last). The burn-up
   * releases its heat shell after shell (setSequentialShell), i.e. an autoignition front sweeps the end
   * gas toward the flame over the burn-up time τ_ab — the Zeldovich (1980) spontaneous-ignition front of
   * a stratified end gas (file header). Round 1 released it uniformly over the whole region, whose
   * projection on every transverse mode vanishes as the region approaches the whole chamber: the
   * heaviest knock (near-homogeneous autoignition of low-ON fuels) rang with MAPO → 0. UNVERIFIED: the
   * sweep direction (outside-in); real end gases autoignite from several centres (König & Sheppard 1990,
   * SAE 902135, abstract), any non-uniform sequence excites the modes. Sets sourceShape to shell 0 and
   * returns the end-gas fraction of the source domain.
   */
  setSequentialEndGasOutsideRadius(rInner: number, K = 16): number {
    const N = this.nModes;
    const P = this.projector;
    const aCyl = P.totalMeasure;
    const aEg = P.outsideMeasure(rInner);
    this.endGasAreaFraction = aEg / aCyl;
    if (!(aEg > 0) || K < 1) {
      this.setEndGasOutsideRadius(rInner);
      return this.endGasAreaFraction;
    }
    if (this.shellShapes.length < K * N) this.shellShapes = new Float64Array(K * N);
    const prev = new Float64Array(N);
    const cur = new Float64Array(N);
    let aPrev = 0;
    let sPrev = P.maxRadius; // nothing outside
    for (let k = 0; k < K; k++) {
      // radius s_{k+1} with (k+1)/K of the end-gas measure outside it
      const s = k === K - 1 ? rInner : P.radiusForOutsideMeasure(((k + 1) / K) * aEg, rInner, sPrev);
      const a = P.outsideIntegrals(s, cur);
      const da = a - aPrev;
      for (let j = 0; j < N; j++) this.shellShapes[k * N + j] = da > 0 ? (cur[j] - prev[j]) / da / this.meanSquare[j] : 0;
      prev.set(cur);
      aPrev = a;
      sPrev = s;
    }
    this.shellCount = K;
    this.currentShell = -1; // force the copy of the new shell 0 (a previous event may have left shell 0 current)
    this.setSequentialShell(0);
    return this.endGasAreaFraction;
  }

  /**
   * Current source shape at sweep position F ∈ [0, 1] (the autoignited fraction of the end gas):
   * linear interpolation between the projections of adjacent shells taken at their mid-positions
   * (k + ½)/K — the front moves continuously (no shell-switch impulses whose timing would depend on
   * the sub-steps).
   */
  setSequentialPosition(F: number): void {
    const K = this.shellCount;
    if (K === 0) return;
    let x = F * K - 0.5;
    if (x < 0) x = 0;
    if (x > K - 1) x = K - 1;
    const k0 = Math.floor(x);
    const k1 = k0 < K - 1 ? k0 + 1 : k0;
    const f = x - k0;
    const N = this.nModes;
    for (let j = 0; j < N; j++) this.sourceShape[j] = (1 - f) * this.shellShapes[k0 * N + j] + f * this.shellShapes[k1 * N + j];
    this.currentShell = k0;
  }

  /** Make shell k (0 = first to autoignite) of the sequential source the current source shape. */
  setSequentialShell(k: number): void {
    if (this.shellCount === 0) return;
    const kk = k < 0 ? 0 : k >= this.shellCount ? this.shellCount - 1 : k;
    if (kk === this.currentShell) return;
    const N = this.nModes;
    for (let j = 0; j < N; j++) this.sourceShape[j] = this.shellShapes[kk * N + j];
    this.currentShell = kk;
  }
  /**
   * Advance all modes by dt (s). c: current sound speed (m/s) setting ω_j = 2α_j c/B;
   * gamma: ratio of specific heats; volume: cylinder volume (m³); qDot: end-gas heat-release
   * rate over this step (W, piecewise constant — its jump at the step start enters as the
   * exact impulse (γ−1)(⟨ψ⟩_eg/(V⟨ψ²⟩)) ΔQ̇ on η̇). Exact damped-oscillator propagation:
   * unconditionally stable for any dt. Allocation-free.
   * qDecayTime (s, optional): the source decays as Q̇(t) = qDot e^{−t/qDecayTime} over the step
   * (first-order end-gas burn-up) — propagated exactly (particular solution), so the response does
   * not depend on how the burn-up is split into steps; the next step's start value is qDot e^{−dt/τ}.
   */
  step(dt: number, c: number, gamma: number, volume: number, qDot: number, qDecayTime = Infinity): void {
    // a non-finite source would poison every modal amplitude until the next reset(): ignore it
    // (the caller must never produce one; this only keeps the reported pressure finite)
    if (!Number.isFinite(qDot)) qDot = this.qPrev;
    // exponentially decaying source Q̇(t) = qDot e^{−t/τ_q} over the step (τ_q = ∞: constant):
    // exact particular solution η_p = P_j e^{−t/τ_q}, P_j = −K_j qDot/τ_q / (ω_j² − 2δ_j/τ_q + 1/τ_q²)
    const expQ = qDecayTime < Infinity && qDecayTime > 0 && qDot !== 0;
    const eq = expQ ? Math.exp(-dt / qDecayTime) : 1;
    const invTq = expQ ? 1 / qDecayTime : 0;
    this.qPrev = qDot * eq;
    const gv = (gamma - 1) / volume;
    // impulse on η̇_j from the jump of the modal forcing Q̇·S_j at the step start (the source shape
    // S_j may change between steps: sequential autoignition, setSequentialShell)
    const imp = this.impScratch;
    const qs = this.qsPrev;
    for (let j = 0; j < this.nModes; j++) imp[j] = gv * (qDot * this.sourceShape[j] - qs[j]);
    const w0 = (2 * c) / this.bore;
    const d0 = 1 / this.decayTime;
    // sensor pressure and its rate at the start of the step (after the source impulse), for the
    // between-sample peak (cubic Hermite) — only when there is something to track
    const track = this.trackPeak;
    let p0 = 0;
    let r0 = 0;
    const gb = this.gainScratch;
    if (track) {
      for (let j = 0; j < this.nModes; j++) {
        gb[j] = this.bandGain(j, c);
        const ws = this.psiSensor[j] * gb[j];
        p0 += ws * this.eta[j];
        r0 += ws * (this.etaDot[j] + imp[j]);
      }
    }
    for (let j = 0; j < this.nModes; j++) {
      let v = this.etaDot[j] + imp[j];
      let e = this.eta[j];
      const w = w0 * this.alpha[j];
      const d = d0 * this.dampScale[j];
      let P = 0;
      if (expQ && this.sourceShape[j] !== 0) {
        const den = w * w - 2 * d * invTq + invTq * invTq;
        if (den > 1e-300) {
          P = (-gv * this.sourceShape[j] * qDot * invTq) / den;
          e -= P; // homogeneous part
          v += P * invTq;
        }
      }
      if (e === 0 && v === 0 && P === 0) continue;
      // η(t) = e^{−dt} [η0 C + (η̇0 + d η0) S],  η̇(t) = e^{−dt} [η̇0 C − (ω² η0 + d η̇0) S]
      // with C = cos(ω_d t), S = sin(ω_d t)/ω_d (underdamped, ω_d² = ω² − d² > 0),
      // C = cosh(s t), S = sinh(s t)/s (overdamped, s² = d² − ω² > 0), their common limit
      // C = 1, S = t when ω_d → 0 (overdamped modes occur for decayTime ≲ 1/ω or c → 0).
      const wd2 = w * w - d * d;
      const wd = Math.sqrt(Math.abs(wd2));
      const x = wd * dt;
      let eC: number; // e^{−d dt} C
      let eS: number; // e^{−d dt} S
      if (x < 1e-4) {
        const ex = Math.exp(-d * dt);
        const x2 = wd2 * dt * dt; // signed: > 0 underdamped
        eC = ex * (1 - 0.5 * x2);
        eS = ex * dt * (1 - x2 / 6);
      } else if (wd2 > 0) {
        const ex = Math.exp(-d * dt);
        eC = ex * Math.cos(x);
        eS = (ex * Math.sin(x)) / wd;
      } else {
        // overdamped: combine the exponentials so e^{−d dt}·cosh never forms ∞·0
        const ep = Math.exp((wd - d) * dt);
        const em = Math.exp(-(wd + d) * dt);
        eC = 0.5 * (ep + em);
        eS = (0.5 * (ep - em)) / wd;
      }
      this.eta[j] = e * eC + (v + d * e) * eS + P * eq;
      this.etaDot[j] = v * eC - (w * w * e + d * v) * eS - P * invTq * eq;
    }
    for (let j = 0; j < this.nModes; j++) this.qsPrev[j] = qDot * eq * this.sourceShape[j];
    if (track) {
      let p1 = 0;
      let r1 = 0;
      for (let j = 0; j < this.nModes; j++) {
        const ws = this.psiSensor[j] * gb[j];
        p1 += ws * this.eta[j];
        r1 += ws * this.etaDot[j];
      }
      const a = hermitePeakAbs(p0, r0, p1, r1, dt);
      if (a > this.peakSensorPressure) this.peakSensorPressure = a;
    }
  }

  /**
   * Maximum |sensor pressure| since the last reset(), Pa, including the maxima BETWEEN step ends:
   * on each step the sensor signal is interpolated by the cubic Hermite polynomial of its values
   * and exact rates (Σψη, Σψη̇) at both ends, whose local extrema are evaluated. With steps ≲ 1/8 of
   * the shortest mode period the result is independent of the step placement to ≲ 1e-4 (the
   * interpolation error of a sinusoid is (ωΔt)⁴/384). Tracked only while `trackPeak` is true.
   * With a measurement band (options.band) the tracked signal is Σ_j g_j ψ_j η_j (bandGain; the gains
   * are those of each step's sound speed), i.e. MAPO in the band; sensorEnvelope stays the all-mode bound.
   */
  peakSensorPressure = 0;
  /** Enable the between-sample peak tracking of {@link peakSensorPressure} (default true). */
  trackPeak = true;

  /**
   * Upper bound of |sensor pressure| at any later time without further excitation, Pa:
   * Σ_j |ψ_j(sensor)| A_j with the free-oscillation amplitude A_j = √(η_j² + ((η̇_j + δ_j η_j)/ω_d,j)²)
   * of each (underdamped) mode, which only decays (overdamped modes: |η| + |η̇|/δ bound).
   * c: current sound speed, m/s.
   */
  sensorEnvelope(c: number): number {
    const w0 = (2 * c) / this.bore;
    const d0 = 1 / this.decayTime;
    let s = 0;
    for (let j = 0; j < this.nModes; j++) {
      const e = this.eta[j];
      const v = this.etaDot[j];
      if (e === 0 && v === 0) continue;
      const w = w0 * this.alpha[j];
      const d = d0 * this.dampScale[j];
      const wd2 = w * w - d * d;
      let A: number;
      if (wd2 > 0) {
        const q = (v + d * e) / Math.sqrt(wd2);
        A = Math.sqrt(e * e + q * q);
      } else A = Math.abs(e) + Math.abs(v) / d;
      s += Math.abs(this.psiSensor[j]) * A;
    }
    return s;
  }

  /** Synthesised pressure oscillation at the configured sensor, Pa. */
  sensorPressure(): number {
    let p = 0;
    for (let j = 0; j < this.nModes; j++) p += this.eta[j] * this.psiSensor[j];
    return p;
  }

  /** Pressure oscillation at planform point (x, z), Pa (evaluates Bessel functions; not hot). */
  pressureAt(x: number, z: number): number {
    let p = 0;
    for (let j = 0; j < this.nModes; j++) if (this.eta[j] !== 0) p += this.eta[j] * this.modeShape(j, x, z);
    return p;
  }

  /**
   * Total acoustic energy in the modes, J: Σ V⟨ψ²⟩ (η² + (η̇/ω)²) / (2ρc²)
   * (potential + kinetic energy of standing waves). rho: kg/m³, c: m/s, volume: m³.
   */
  acousticEnergy(rho: number, c: number, volume: number): number {
    const w0 = (2 * c) / this.bore;
    let e = 0;
    for (let j = 0; j < this.nModes; j++) {
      const w = w0 * this.alpha[j];
      const ed = this.etaDot[j] / w;
      e += this.meanSquare[j] * (this.eta[j] * this.eta[j] + ed * ed);
    }
    return (e * volume) / (2 * rho * c * c);
  }
}

/**
 * Maximum of |P(s)| over s ∈ [0, h] for the cubic Hermite interpolant with P(0) = p0, P′(0) = r0,
 * P(h) = p1, P′(h) = r1 (end values and interior extrema, i.e. the real roots of P′ in (0, h)).
 */
export function hermitePeakAbs(p0: number, r0: number, p1: number, r1: number, h: number): number {
  const m = Math.max(Math.abs(p0), Math.abs(p1));
  if (!(h > 0)) return m;
  // P(u) = p0 + a u + c2 u² + c3 u³ on u ∈ [0, 1];  P′(u) = a + 2 c2 u + 3 c3 u²
  const a = r0 * h;
  const b = r1 * h;
  const d = p1 - p0;
  const c2 = 3 * d - 2 * a - b;
  const c3 = a + b - 2 * d;
  const A = 3 * c3;
  const B = 2 * c2;
  let u1 = -1;
  let u2 = -1;
  if (Math.abs(A) <= 1e-14 * (Math.abs(B) + Math.abs(a))) {
    if (B !== 0) u1 = -a / B;
  } else {
    const disc = B * B - 4 * A * a;
    if (disc >= 0) {
      const sq = Math.sqrt(disc);
      const q = -0.5 * (B + (B >= 0 ? sq : -sq));
      if (q !== 0) {
        u1 = q / A;
        u2 = a / q;
      }
    }
  }
  return Math.max(m, cubicAbsAt(p0, a, c2, c3, u1), cubicAbsAt(p0, a, c2, c3, u2));
}

/** |p0 + a u + c2 u² + c3 u³| for u ∈ (0, 1), else 0. */
function cubicAbsAt(p0: number, a: number, c2: number, c3: number, u: number): number {
  return u > 0 && u < 1 ? Math.abs(p0 + u * (a + u * (c2 + u * c3))) : 0;
}

// ---------------------------------------------------------------------------------------
// End-gas burn-up after autoignition
// ---------------------------------------------------------------------------------------

/**
 * Excitation time of PRF/air autoignition, s: time from 5 % of the peak heat-release rate to
 * the peak (Lutz et al. 1988) in an adiabatic constant-volume reactor, LLNL gasoline-surrogate
 * 2011 mechanism (Mehl et al.), φ = 1, 5 % residual, T0 850–1050 K, 30–60 bar, PRF 90/100:
 * geometric mean 1.01 µs of test/fixtures/chemistry_knock_excitation.json (range 0.69–1.47 µs,
 * decreasing with p and T; LLNL PRF v2 gives 0.72–1.47 µs; tools/reference/chemistry_knock_excitation.py).
 */
export const EXCITATION_TIME_PRF = 1.01e-6;

/**
 * Temperature spread of the end gas that autoignites sequentially, K.
 * UNVERIFIED: 15 K — order of the "natural thermal stratification" of the bulk charge
 * reported for HCCI engines (Sjöberg & Dec, SAE 2005-01-0113; Dec, Hwang & Sjöberg,
 * SAE 2006-01-1518); numbers not checked against the papers.
 */
export const END_GAS_STRATIFICATION_DT = 15;

/**
 * End-gas autoignition burn-up time, s: τ_ab = τ_e + τ |∂lnτ/∂T| ΔT — sequential autoignition
 * of end gas with a temperature spread ΔT (K) given the current ignition delay τ (s) and its
 * temperature sensitivity dLnTauDT (1/K), plus the excitation time τ_e (s).
 */
export function autoignitionBurnTime(
  tau: number,
  dLnTauDT: number,
  deltaT: number = END_GAS_STRATIFICATION_DT,
  excitationTime: number = EXCITATION_TIME_PRF,
): number {
  return excitationTime + tau * Math.abs(dLnTauDT) * deltaT;
}

/**
 * End-gas burn rate after autoignition, kg/s: first-order (volumetric) consumption
 * dm_b/dt = m_u / τ_ab of the remaining unburned mass mUnburned (kg).
 */
export function endGasBurnRate(mUnburned: number, burnTime: number): number {
  return mUnburned > 0 ? mUnburned / burnTime : 0;
}

/**
 * End-gas mass burned over a step dt (s), kg: the EXACT integral of dm_u/dt = −m_u/τ_ab with m_u
 * (kg) at the step start, Δm = m_u (1 − e^{−dt/τ_ab}). τ_ab (≈ 1–100 µs) is usually far shorter
 * than the cycle time step (0.1° = 28 µs at 600 rpm), where the explicit rate endGasBurnRate
 * is unstable (Euler overshoots to m_u < 0 once dt > τ_ab); use this (operator split) instead,
 * with the step-mean heat-release rate Q̇ = Δm·q/dt fed to KnockOscillator.step. burnTime ≤ 0
 * burns everything in the step; NaN propagates.
 */
export function endGasBurnedMass(mUnburned: number, burnTime: number, dt: number): number {
  if (!(mUnburned > 0) || !(dt > 0)) return 0;
  if (burnTime <= 0) return mUnburned;
  return -mUnburned * Math.expm1(-dt / burnTime); // NaN burnTime → NaN
}

/** cos of the inside-arc half-angle of a ring of radius r about a circle (centre distance d, radius rf). */
function ringCos(r: number, d: number, rf: number): number {
  if (d < 1e-15 * (r + rf) || r === 0) return r <= rf ? -2 : 2; // concentric: all inside / outside
  return (r * r + d * d - rf * rf) / (2 * r * d);
}

const glCache = new Map<number, { x: Float64Array; w: Float64Array }>();

/**
 * Gauss–Legendre nodes and weights on [−1, 1] (n points; Newton on P_n from the Tricomi initial
 * guesses, Press et al., Numerical Recipes §4.6 gauleg). Cached per n.
 */
export function gaussLegendre(n: number): { x: Float64Array; w: Float64Array } {
  let g = glCache.get(n);
  if (g) return g;
  const x = new Float64Array(n);
  const w = new Float64Array(n);
  const m = (n + 1) >> 1;
  for (let i = 0; i < m; i++) {
    let z = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5));
    let pp = 0;
    for (let it = 0; it < 100; it++) {
      let p1 = 1;
      let p2 = 0;
      for (let j = 0; j < n; j++) {
        const p3 = p2;
        p2 = p1;
        p1 = ((2 * j + 1) * z * p2 - j * p3) / (j + 1);
      }
      pp = (n * (z * p1 - p2)) / (z * z - 1);
      const dz = p1 / pp;
      z -= dz;
      if (Math.abs(dz) < 1e-15) break;
    }
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = 2 / ((1 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
  g = { x, w };
  glCache.set(n, g);
  return g;
}

/** Intersection area of a circle of radius r (centre at distance d) with the disc of radius R, m². */
export function circleLensArea(r: number, R: number, d: number): number {
  if (r <= 0) return 0;
  if (d + r <= R) return Math.PI * r * r;
  if (d + R <= r) return Math.PI * R * R;
  if (d >= r + R) return 0;
  const a = Math.acos(Math.max(-1, Math.min(1, (d * d + r * r - R * R) / (2 * d * r))));
  const b = Math.acos(Math.max(-1, Math.min(1, (d * d + R * R - r * r) / (2 * d * R))));
  return r * r * a + R * R * b - 0.5 * Math.sqrt(Math.max(0, (-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R)));
}

// ---------------------------------------------------------------------------------------
// Non-cylindrical chambers: depth-weighted planform modes on a cut-cell grid
// ---------------------------------------------------------------------------------------

/**
 * A chamber made of vertical gas columns over a planform: at planform point (x, z) of region k the gas
 * fills yLo[k] ≤ y ≤ yHi[k] (cylinder frame; flat roofs and floors, vertical walls). The flat-disc
 * chamber (one region) and the side-valve L-head (bore column + valve pocket) are of this kind; the
 * local chamber height is H(x, z) = yHi[k] − yLo[k].
 */
export interface ChamberPlanform {
  /** Planform bounding box, cylinder frame, m. */
  readonly xMin: number;
  readonly xMax: number;
  readonly zMin: number;
  readonly zMax: number;
  /** Gas interval of each region along the cylinder axis, m. */
  readonly yLo: readonly number[];
  readonly yHi: readonly number[];
  /** Region index at planform point (x, z); −1 outside the chamber. */
  region(x: number, z: number): number;
  /**
   * A lower bound of the distance from (x, z) to the nearest region boundary (the chamber outline or a
   * depth step), m: grid cells farther than this from every boundary are uniform.
   */
  boundaryDistance(x: number, z: number): number;
}

/** The flat-disc chamber of radius R and depth h (head face y = 0) as a ChamberPlanform. */
export function discPlanform(radius: number, depth: number): ChamberPlanform {
  if (!(radius > 0) || !(depth > 0)) throw new RangeError('discPlanform: radius and depth must be positive');
  const r2 = radius * radius;
  return {
    xMin: -radius,
    xMax: radius,
    zMin: -radius,
    zMax: radius,
    yLo: [-depth],
    yHi: [0],
    region: (x, z) => (x * x + z * z <= r2 ? 0 : -1),
    boundaryDistance: (x, z) => Math.abs(Math.hypot(x, z) - radius),
  };
}

/** Bore-column depth of an L-head chamber at TDC (roof of the head cavity over the bore → crown), m. */
export function lHeadDepthAtTDC(lHead: LHeadChamberSpec): number {
  return -(lHead.deckY + lHead.crownAboveDeckAtTDC);
}

/** Signed distance from (x, z) to a rounded rectangle (negative inside), m (exact). */
export function roundedRectDistance(x: number, z: number, xMin: number, xMax: number, zMin: number, zMax: number, cornerRadius: number): number {
  const hx = 0.5 * (xMax - xMin) - cornerRadius;
  const hz = 0.5 * (zMax - zMin) - cornerRadius;
  const qx = Math.abs(x - 0.5 * (xMin + xMax)) - hx;
  const qz = Math.abs(z - 0.5 * (zMin + zMax)) - hz;
  return Math.hypot(Math.max(qx, 0), Math.max(qz, 0)) + Math.min(Math.max(qx, qz), 0) - cornerRadius;
}

/**
 * The side-valve chamber of engine-spec.ts LHeadChamberSpec as a ChamberPlanform: region 0 = the bore
 * column x² + z² ≤ (B/2)², −h ≤ y ≤ 0 (h: roof-to-crown depth, default TDC); region 1 = the valve pocket
 * (rounded rectangle minus the bore disc), deckY ≤ y ≤ pocket.roofY. The two connect across the bore
 * circle; in the depth-averaged model only their heights matter (file header).
 */
export function lHeadPlanform(bore: number, lHead: LHeadChamberSpec, depth = lHeadDepthAtTDC(lHead)): ChamberPlanform {
  const R = bore / 2;
  const p = lHead.pocket;
  const rc = p.cornerRadius;
  if (!(depth > 0)) throw new RangeError(`lHeadPlanform: bore-column depth ${depth} m must be positive`);
  if (!(p.roofY > lHead.deckY)) throw new RangeError('lHeadPlanform: pocket roof must lie above the deck');
  if (!(rc >= 0) || 2 * rc > p.xMax - p.xMin || 2 * rc > p.zMax - p.zMin) throw new RangeError('lHeadPlanform: invalid pocket corner radius');
  const r2 = R * R;
  return {
    xMin: Math.min(-R, p.xMin),
    xMax: Math.max(R, p.xMax),
    zMin: Math.min(-R, p.zMin),
    zMax: Math.max(R, p.zMax),
    yLo: [-depth, lHead.deckY],
    yHi: [0, p.roofY],
    region: (x, z) => (x * x + z * z <= r2 ? 0 : roundedRectDistance(x, z, p.xMin, p.xMax, p.zMin, p.zMax, rc) <= 0 ? 1 : -1),
    boundaryDistance: (x, z) => Math.min(Math.abs(Math.hypot(x, z) - R), Math.abs(roundedRectDistance(x, z, p.xMin, p.xMax, p.zMin, p.zMax, rc))),
  };
}

/** Options of GridModeSet. */
export interface GridModeSetOptions {
  /** Grid cell size, m (default bore/64: the disc's α_mn within 0.2 %, knock-lhead.test.ts). */
  cellSize?: number;
  /** Keep every mode with α_j = k_j L/2 ≤ maxAlpha (default 7.1, as the disc). */
  maxAlpha?: number;
  /** α_ref of the damping law δ_j = (α_j/α_ref)^½/τ_d (default the bore's (1,0) zero α_10 = 1.8412). */
  dampingReferenceAlpha?: number;
  /** Krylov block size (default 2: resolves double eigenvalues of symmetric planforms). */
  blockSize?: number;
  /** Relative residual of converged Ritz pairs (default 1e-11). */
  tolerance?: number;
}

/** Sub-samples per cell side for cut cells and faces (8 × 8 per cell, 8 per face). */
const GRID_SUBSAMPLES = 8;
/** A cell or face whose centre is farther than this × Δ from every region boundary is uniform (> √2/2). */
const GRID_UNIFORM_DISTANCE = 0.75;

/**
 * Transverse acoustic modes of a ChamberPlanform: the depth-weighted Helmholtz problem
 * ∇·(H ∇ψ) + k² H ψ = 0 with H ∂ψ/∂n = 0 on the outline (file header, NON-CYLINDRICAL CHAMBERS),
 * discretised by cut-cell finite volumes on a Cartesian grid of cell size Δ:
 *  - capacity of a cell c_i = ∫_cell H dA (exact for uniform cells; 8 × 8 midpoint sub-samples where an
 *    outline or a depth step crosses it) — the gas volume of the cell column;
 *  - conductance of the face between neighbours G = (1/Δ)∫_face H_h ds, H_h = harmonic mean of H along the
 *    segment joining the two cell centres (series resistance across a depth step), face points outside
 *    the chamber contributing nothing (the rigid outline: Neumann condition with no extra term);
 *  - K ψ = k² C ψ (K: the conductance Laplacian, C = diag c_i) — symmetric, with the constant ψ_0 for
 *    k = 0, so every other mode satisfies Σ c_i ψ_j,i = 0 exactly (a uniform release excites nothing).
 * The lowest modes come from a shift-invert block Lanczos iteration with full reorthogonalisation on an
 * envelope Cholesky factor of C^{-½}(K − σC)C^{-½} (σ < 0), Rayleigh–Ritz by Householder + implicit QL.
 * Accuracy (knock-lhead.test.ts): the disc's α_mn within 0.2 % at Δ = B/64, the two-depth channel within
 * 0.2 % of its analytic roots, the Model T α_j within 0.3 % of a B/192 scipy solution.
 * Mode shapes are normalised to max |ψ_j| = 1 over the cells (largest entry positive). Immutable.
 */
export class GridModeSet implements AcousticModeSet {
  readonly bore: number;
  readonly nModes: number;
  readonly alpha: Float64Array;
  readonly meanSquare: Float64Array;
  readonly dampScale: Float64Array;
  readonly m: Int32Array;
  readonly n: Int32Array;
  readonly cosine: Uint8Array;
  readonly planform: ChamberPlanform;
  readonly cellSize: number;
  /** Grid of nx × nz cells with lower-left corner (gx0, gz0); cell (i, k) is centred at gx0 + (i + ½)Δ, gz0 + (k + ½)Δ. */
  readonly nx: number;
  readonly nz: number;
  readonly gx0: number;
  readonly gz0: number;
  /** Active-cell index of grid cell (i, k) at [i + nx·k]; −1 outside the chamber. */
  readonly cellIndex: Int32Array;
  readonly nCells: number;
  /** Cell capacities ∫_cell H dA (gas volume of each cell column), m³. */
  readonly capacity: Float64Array;
  /** Chamber volume Σ c_i, m³. */
  readonly volume: number;
  /** Mode shapes at the cell centres, cell-major: psi[c·nModes + j]. */
  readonly psi: Float64Array;
  /**
   * Gas sub-columns (one per cell and region present in it — the source projector's quadrature points):
   * planform centroid (x, z), plan area, gas interval [yLo, yHi], cell index.
   */
  readonly columns: {
    readonly count: number;
    readonly x: Float64Array;
    readonly z: Float64Array;
    readonly area: Float64Array;
    readonly yLo: Float64Array;
    readonly yHi: Float64Array;
    readonly cell: Int32Array;
  };
  /**
   * Solver diagnostics: Krylov dimension, largest relative Ritz residual of the kept modes, k² of the
   * discarded uniform mode (1/m², ≈ 0), build time (ms).
   */
  readonly stats: { krylovDimension: number; maxResidual: number; nullEigenvalue: number; buildMs: number };

  constructor(bore: number, planform: ChamberPlanform, opts: GridModeSetOptions = {}) {
    const t0 = typeof performance !== 'undefined' ? performance.now() : Date.now();
    this.bore = bore;
    this.planform = planform;
    const dx = opts.cellSize ?? bore / 64;
    if (!(dx > 0)) throw new RangeError('GridModeSet: cell size must be positive');
    this.cellSize = dx;
    const nx = Math.max(1, Math.ceil((planform.xMax - planform.xMin) / dx - 1e-9));
    const nz = Math.max(1, Math.ceil((planform.zMax - planform.zMin) / dx - 1e-9));
    this.nx = nx;
    this.nz = nz;
    this.gx0 = 0.5 * (planform.xMin + planform.xMax) - 0.5 * nx * dx;
    this.gz0 = 0.5 * (planform.zMin + planform.zMax) - 0.5 * nz * dx;
    const hOf = (r: number): number => (r >= 0 ? planform.yHi[r] - planform.yLo[r] : 0);
    const ns = GRID_SUBSAMPLES;
    const dUni = GRID_UNIFORM_DISTANCE * dx;
    // ---- cells: sub-columns and capacities ----
    const cellIndex = new Int32Array(nx * nz).fill(-1);
    const capList: number[] = [];
    const cxList: number[] = [];
    const czList: number[] = [];
    const col = { x: [] as number[], z: [] as number[], area: [] as number[], yLo: [] as number[], yHi: [] as number[], cell: [] as number[] };
    const nReg = planform.yLo.length;
    const cnt = new Float64Array(nReg);
    const sx = new Float64Array(nReg);
    const sz = new Float64Array(nReg);
    const subArea = (dx * dx) / (ns * ns);
    for (let k = 0; k < nz; k++) {
      for (let i = 0; i < nx; i++) {
        const x = this.gx0 + (i + 0.5) * dx;
        const z = this.gz0 + (k + 0.5) * dx;
        const first = col.x.length;
        if (planform.boundaryDistance(x, z) > dUni) {
          const r = planform.region(x, z);
          if (r >= 0) pushColumn(col, x, z, dx * dx, planform.yLo[r], planform.yHi[r]);
        } else {
          cnt.fill(0);
          sx.fill(0);
          sz.fill(0);
          for (let b = 0; b < ns; b++) {
            const zs = z + ((b + 0.5) / ns - 0.5) * dx;
            for (let a = 0; a < ns; a++) {
              const xs = x + ((a + 0.5) / ns - 0.5) * dx;
              const r = planform.region(xs, zs);
              if (r < 0) continue;
              cnt[r]++;
              sx[r] += xs;
              sz[r] += zs;
            }
          }
          for (let r = 0; r < nReg; r++) {
            if (cnt[r] > 0) pushColumn(col, sx[r] / cnt[r], sz[r] / cnt[r], cnt[r] * subArea, planform.yLo[r], planform.yHi[r]);
          }
        }
        let cap = 0;
        for (let q = first; q < col.x.length; q++) cap += col.area[q] * (col.yHi[q] - col.yLo[q]);
        if (cap > 0) {
          const c = capList.length;
          cellIndex[i + nx * k] = c;
          capList.push(cap);
          cxList.push(x);
          czList.push(z);
          for (let q = first; q < col.x.length; q++) col.cell[q] = c;
        } else {
          col.x.length = col.z.length = col.area.length = col.yLo.length = col.yHi.length = col.cell.length = first;
        }
      }
    }
    const nC = capList.length;
    if (nC < 2) throw new RangeError('GridModeSet: the planform covers fewer than 2 cells');
    this.cellIndex = cellIndex;
    this.nCells = nC;
    this.capacity = Float64Array.from(capList);
    let vol = 0;
    for (let c = 0; c < nC; c++) vol += this.capacity[c];
    this.volume = vol;
    this.columns = {
      count: col.x.length,
      x: Float64Array.from(col.x),
      z: Float64Array.from(col.z),
      area: Float64Array.from(col.area),
      yLo: Float64Array.from(col.yLo),
      yHi: Float64Array.from(col.yHi),
      cell: Int32Array.from(col.cell),
    };
    // ---- faces: conductances to the +x and +z neighbours ----
    const gxp = new Float64Array(nC); // to (i+1, k)
    const gzp = new Float64Array(nC); // to (i, k+1)
    const face = (xf: number, zf: number, ex: number, ez: number): number => {
      if (planform.boundaryDistance(xf, zf) > dUni) return hOf(planform.region(xf, zf));
      let tot = 0;
      for (let a = 0; a < ns; a++) {
        const u = ((a + 0.5) / ns - 0.5) * dx;
        const px = xf + u * ez; // tangential direction (ez, ex)
        const pz = zf + u * ex;
        const hf = hOf(planform.region(px, pz));
        if (!(hf > 0)) continue;
        let inv = 0;
        let nIn = 0;
        for (let b = 0; b < ns; b++) {
          const v = ((b + 0.5) / ns - 0.5) * dx;
          const h = hOf(planform.region(px + v * ex, pz + v * ez));
          if (h > 0) {
            inv += 1 / h;
            nIn++;
          }
        }
        tot += nIn > 0 ? nIn / inv : hf;
      }
      return tot / ns;
    };
    for (let k = 0; k < nz; k++) {
      for (let i = 0; i < nx; i++) {
        const c = cellIndex[i + nx * k];
        if (c < 0) continue;
        const x = cxList[c];
        const z = czList[c];
        if (i + 1 < nx && cellIndex[i + 1 + nx * k] >= 0) gxp[c] = face(x + 0.5 * dx, z, 1, 0);
        if (k + 1 < nz && cellIndex[i + nx * (k + 1)] >= 0) gzp[c] = face(x, z + 0.5 * dx, 0, 1);
      }
    }
    // ---- eigen-solve ----
    const maxAlpha = opts.maxAlpha ?? 7.1;
    const lamCut = ((2 * maxAlpha) / bore) ** 2;
    const sol = lowestGridModes(this, gxp, gzp, lamCut, opts.blockSize ?? 2, opts.tolerance ?? 1e-11);
    const N = sol.lambda.length;
    this.nModes = N;
    this.alpha = new Float64Array(N);
    this.meanSquare = new Float64Array(N);
    this.dampScale = new Float64Array(N);
    this.m = new Int32Array(N).fill(-1);
    this.n = Int32Array.from({ length: N }, (_, j) => j);
    this.cosine = new Uint8Array(N).fill(1);
    this.psi = new Float64Array(nC * N);
    const aRef = opts.dampingReferenceAlpha ?? modeZero(1, 0);
    for (let j = 0; j < N; j++) {
      const a = (Math.sqrt(Math.max(sol.lambda[j], 0)) * bore) / 2;
      this.alpha[j] = a;
      this.dampScale[j] = Math.sqrt(a / aRef);
      const v = sol.vectors[j];
      let big = 0;
      let iBig = 0;
      for (let c = 0; c < nC; c++) {
        if (Math.abs(v[c]) > big) {
          big = Math.abs(v[c]);
          iBig = c;
        }
      }
      const s = v[iBig] < 0 ? -1 / big : 1 / big;
      let ms = 0;
      for (let c = 0; c < nC; c++) {
        const p = v[c] * s;
        this.psi[c * N + j] = p;
        ms += this.capacity[c] * p * p;
      }
      this.meanSquare[j] = ms / vol;
    }
    const t1 = typeof performance !== 'undefined' ? performance.now() : Date.now();
    this.stats = { krylovDimension: sol.krylovDimension, maxResidual: sol.maxResidual, nullEigenvalue: sol.nullEigenvalue, buildMs: t1 - t0 };
  }

  /**
   * Mode shape ψ_j at planform point (x, z): bilinear interpolation of the cell-centre values over the
   * active cells among the four nearest (weights renormalised); 0 outside the chamber.
   */
  modeShape(j: number, x: number, z: number): number {
    if (this.planform.region(x, z) < 0) return 0;
    const dx = this.cellSize;
    const fx = (x - this.gx0) / dx - 0.5;
    const fz = (z - this.gz0) / dx - 0.5;
    const i0 = Math.floor(fx);
    const k0 = Math.floor(fz);
    const tx = fx - i0;
    const tz = fz - k0;
    let s = 0;
    let w = 0;
    for (let dk = 0; dk <= 1; dk++) {
      const k = k0 + dk;
      if (k < 0 || k >= this.nz) continue;
      for (let di = 0; di <= 1; di++) {
        const i = i0 + di;
        if (i < 0 || i >= this.nx) continue;
        const c = this.cellIndex[i + this.nx * k];
        if (c < 0) continue;
        const wt = (di ? tx : 1 - tx) * (dk ? tz : 1 - tz);
        s += wt * this.psi[c * this.nModes + j];
        w += wt;
      }
    }
    return w > 0 ? s / w : 0;
  }
}

function pushColumn(
  col: { x: number[]; z: number[]; area: number[]; yLo: number[]; yHi: number[]; cell: number[] },
  x: number,
  z: number,
  area: number,
  yLo: number,
  yHi: number,
): void {
  col.x.push(x);
  col.z.push(z);
  col.area.push(area);
  col.yLo.push(yLo);
  col.yHi.push(yHi);
  col.cell.push(-1);
}

/**
 * Lowest eigenpairs (k² ≤ lamCut, the k = 0 mode dropped) of K ψ = k² C ψ for the cut-cell grid: shift-invert
 * block Lanczos (block Arnoldi with two classical Gram–Schmidt passes against the whole basis, i.e. full
 * reorthogonalisation) on S = C^{-½}(K − σC)C^{-½}, σ = −(1/L)², with S factored by an envelope Cholesky in
 * an ordering that runs along the shorter grid side (envelope ≈ that many cells). Converged when every Ritz
 * pair below the cut has ‖S⁻¹y − θy‖ ≤ tol·θ and at least one block of converged pairs lies above it, or
 * when the basis is complete (n vectors: the exact discrete eigenpairs). A numerically invariant Krylov space
 * before that (the start block's components along the highest modes decay like (θ_min/θ_max)^k) is
 * continued from fresh random directions. Model T at B/64: n = 5418, Krylov dimension ≈ 100 for 26 modes,
 * ≈ 0.15–0.3 s. Build-time only (allocates).
 */
function lowestGridModes(
  g: GridModeSet,
  gxp: Float64Array,
  gzp: Float64Array,
  lamCut: number,
  blockSize: number,
  tol: number,
): { lambda: Float64Array; vectors: Float64Array[]; krylovDimension: number; maxResidual: number; nullEigenvalue: number } {
  const n = g.nCells;
  const { nx, nz, cellIndex, capacity } = g;
  // ordering: inner loop along the shorter grid side
  const pos = new Int32Array(n);
  let p = 0;
  if (nz <= nx) {
    for (let i = 0; i < nx; i++) for (let k = 0; k < nz; k++) if (cellIndex[i + nx * k] >= 0) pos[cellIndex[i + nx * k]] = p++;
  } else {
    for (let k = 0; k < nz; k++) for (let i = 0; i < nx; i++) if (cellIndex[i + nx * k] >= 0) pos[cellIndex[i + nx * k]] = p++;
  }
  // scaled, shifted matrix S (by position): diagonal and the (up to 4) neighbours
  const sigma = -1 / (g.bore * g.bore);
  const isq = new Float64Array(n); // 1/√c by position
  for (let c = 0; c < n; c++) isq[pos[c]] = 1 / Math.sqrt(capacity[c]);
  const diag = new Float64Array(n);
  const nbP: number[][] = Array.from({ length: n }, () => []);
  const nbV: number[][] = Array.from({ length: n }, () => []);
  const link = (ca: number, cb: number, G: number): void => {
    if (!(G > 0)) return;
    const a = pos[ca];
    const b = pos[cb];
    diag[a] += G;
    diag[b] += G;
    const v = -G * isq[a] * isq[b];
    nbP[a].push(b);
    nbV[a].push(v);
    nbP[b].push(a);
    nbV[b].push(v);
  };
  for (let k = 0; k < nz; k++) {
    for (let i = 0; i < nx; i++) {
      const c = cellIndex[i + nx * k];
      if (c < 0) continue;
      if (gxp[c] > 0) link(c, cellIndex[i + 1 + nx * k], gxp[c]);
      if (gzp[c] > 0) link(c, cellIndex[i + nx * (k + 1)], gzp[c]);
    }
  }
  for (let a = 0; a < n; a++) diag[a] = diag[a] * isq[a] * isq[a] - sigma;
  // envelope Cholesky S = L Lᵀ (row a holds columns first[a]..a)
  const first = new Int32Array(n);
  const off = new Int32Array(n + 1);
  for (let a = 0; a < n; a++) {
    let f = a;
    for (const b of nbP[a]) if (b < f) f = b;
    first[a] = f;
    off[a + 1] = off[a] + (a - f + 1);
  }
  const L = new Float64Array(off[n]);
  for (let a = 0; a < n; a++) {
    L[off[a] + a - first[a]] = diag[a];
    for (let t = 0; t < nbP[a].length; t++) {
      const b = nbP[a][t];
      if (b < a) L[off[a] + b - first[a]] = nbV[a][t];
    }
  }
  for (let a = 0; a < n; a++) {
    const fa = first[a];
    const oa = off[a] - fa;
    for (let b = fa; b < a; b++) {
      const fb = first[b];
      const ob = off[b] - fb;
      let s = L[oa + b];
      for (let t = fa > fb ? fa : fb; t < b; t++) s -= L[oa + t] * L[ob + t];
      L[oa + b] = s / L[ob + b];
    }
    let s = L[oa + a];
    for (let t = fa; t < a; t++) s -= L[oa + t] * L[oa + t];
    if (!(s > 0)) throw new Error('GridModeSet: shifted operator not positive definite');
    L[oa + a] = Math.sqrt(s);
  }
  const solve = (x: Float64Array): void => {
    for (let a = 0; a < n; a++) {
      const oa = off[a] - first[a];
      let s = x[a];
      for (let t = first[a]; t < a; t++) s -= L[oa + t] * x[t];
      x[a] = s / L[oa + a];
    }
    for (let a = n - 1; a >= 0; a--) {
      const oa = off[a] - first[a];
      const xa = x[a] / L[oa + a];
      x[a] = xa;
      for (let t = first[a]; t < a; t++) x[t] -= L[oa + t] * xa;
    }
  };
  // ---- block Krylov with full reorthogonalisation ----
  const bs = Math.max(1, Math.min(blockSize, n));
  const mMax = n;
  const Q: Float64Array[] = [];
  const Hc: Float64Array[] = []; // H column j: coefficients of S⁻¹q_j on q_0, q_1, … (Arnoldi)
  const hAt = (i: number, j: number): number => (i < Hc[j].length ? Hc[j][i] : 0);
  let seed = 0x2545f491;
  const rand = (): number => {
    seed = (Math.imul(seed, 1103515245) + 12345) >>> 0;
    return seed / 4294967296 - 0.5;
  };
  const dot = (u: Float64Array, v: Float64Array): number => {
    let s = 0;
    for (let a = 0; a < n; a++) s += u[a] * v[a];
    return s;
  };
  const axpy = (y: Float64Array, s: number, x: Float64Array): void => {
    for (let a = 0; a < n; a++) y[a] += s * x[a];
  };
  /** Orthonormalise w against the basis (two CGS passes), append it if independent; returns the coefficients. */
  const extend = (w: Float64Array, coef: number[] | null): boolean => {
    const w0 = Math.sqrt(dot(w, w));
    for (let pass = 0; pass < 2; pass++) {
      const h = new Float64Array(Q.length);
      for (let i = 0; i < Q.length; i++) h[i] = dot(Q[i], w);
      for (let i = 0; i < Q.length; i++) axpy(w, -h[i], Q[i]);
      if (coef) for (let i = 0; i < Q.length; i++) coef[i] = (coef[i] ?? 0) + h[i];
    }
    const nw = Math.sqrt(dot(w, w));
    if (!(nw > 1e-12 * w0) || Q.length >= mMax) return false;
    for (let a = 0; a < n; a++) w[a] /= nw;
    if (coef) coef[Q.length] = nw;
    Q.push(w);
    return true;
  };
  const addRandom = (count: number): void => {
    for (let b = 0; b < count && Q.length < mMax; b++) {
      for (let attempt = 0; attempt < 4; attempt++) {
        const v = new Float64Array(n);
        for (let a = 0; a < n; a++) v[a] = rand();
        if (extend(v, null)) break;
      }
    }
  };
  addRandom(bs);
  let done = 0; // basis vectors whose images S⁻¹q have been orthogonalised
  let result: { theta: Float64Array; S: Float64Array; P: number; res: Float64Array } | null = null;
  let lastCheck = 0;
  for (;;) {
    const blockEnd = Q.length;
    for (let j = done; j < blockEnd; j++) {
      const w = Float64Array.from(Q[j]);
      solve(w);
      const coef: number[] = [];
      extend(w, coef);
      Hc[j] = Float64Array.from(coef, (v) => v ?? 0);
    }
    done = blockEnd;
    const exhausted = Q.length === done;
    if ((done >= 2 * bs && done - lastCheck >= 4 * bs) || exhausted) {
      lastCheck = done;
      const P = done;
      const T = new Float64Array(P * P);
      for (let i = 0; i < P; i++) for (let j = 0; j < P; j++) T[i * P + j] = 0.5 * (hAt(i, j) + hAt(j, i));
      const theta = new Float64Array(P);
      symmetricEigen(T, P, theta); // T ← eigenvectors (columns), theta ascending
      const res = new Float64Array(P);
      for (let e = 0; e < P; e++) {
        let r2 = 0;
        for (let i = P; i < Q.length; i++) {
          let s = 0;
          for (let j = 0; j < P; j++) s += hAt(i, j) * T[j * P + e];
          r2 += s * s;
        }
        res[e] = Math.sqrt(r2);
      }
      result = { theta, S: T, P, res };
      // convergence: every Ritz value with λ ≤ λcut converged, and ≥ bs converged above the cut
      let ok = true;
      let above = 0;
      for (let e = P - 1; e >= 0; e--) {
        const th = theta[e];
        if (!(th > 0)) continue;
        const lam = sigma + 1 / th;
        const conv = res[e] <= tol * th;
        if (lam <= lamCut) {
          if (!conv) ok = false;
        } else if (conv) above++;
      }
      if (ok && above >= bs) break;
      if (exhausted) {
        if (Q.length >= mMax) break; // complete basis: the Ritz pairs are the exact discrete eigenpairs
        // numerically invariant Krylov space: the start block's components along the highest modes decay
        // like (θ_min/θ_max)^k and are lost to rounding — continue from fresh random directions
        addRandom(bs);
      }
    }
  }
  const { theta, S, P, res } = result!;
  // Ritz pairs, ascending λ (descending θ), drop the k = 0 mode
  const lam: number[] = [];
  const vecs: Float64Array[] = [];
  let maxRes = 0;
  let nullEig = NaN;
  for (let e = P - 1; e >= 0; e--) {
    const th = theta[e];
    if (!(th > 0)) continue;
    const l = sigma + 1 / th;
    if (l > lamCut) break;
    if (Number.isNaN(nullEig)) {
      nullEig = l; // the smallest: the uniform mode
      continue;
    }
    maxRes = Math.max(maxRes, res[e] / th);
    const y = new Float64Array(n);
    for (let j = 0; j < P; j++) {
      const s = S[j * P + e];
      if (s !== 0) axpy(y, s, Q[j]);
    }
    const v = new Float64Array(n);
    for (let c = 0; c < n; c++) v[c] = y[pos[c]] * isq[pos[c]];
    lam.push(l);
    vecs.push(v);
  }
  return { lambda: Float64Array.from(lam), vectors: vecs, krylovDimension: P, maxResidual: maxRes, nullEigenvalue: nullEig };
}

/**
 * Eigen-decomposition of the symmetric n × n matrix A (row-major, overwritten by the eigenvectors as
 * columns), eigenvalues into d in ascending order: Householder tridiagonalisation and the implicit QL
 * algorithm with eigenvector accumulation (EISPACK tred2/tql2, Wilkinson & Reinsch 1971, Handbook for
 * Automatic Computation II; in the form of the public-domain JAMA library). Build-time only.
 */
export function symmetricEigen(A: Float64Array, n: number, d: Float64Array): void {
  const V = A;
  const e = new Float64Array(n);
  const at = (i: number, j: number): number => i * n + j;
  // ---- tred2 ----
  for (let j = 0; j < n; j++) d[j] = V[at(n - 1, j)];
  for (let i = n - 1; i > 0; i--) {
    let scale = 0;
    let h = 0;
    for (let k = 0; k < i; k++) scale += Math.abs(d[k]);
    if (scale === 0) {
      e[i] = d[i - 1];
      for (let j = 0; j < i; j++) {
        d[j] = V[at(i - 1, j)];
        V[at(i, j)] = 0;
        V[at(j, i)] = 0;
      }
    } else {
      for (let k = 0; k < i; k++) {
        d[k] /= scale;
        h += d[k] * d[k];
      }
      let f = d[i - 1];
      let g = Math.sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (let j = 0; j < i; j++) e[j] = 0;
      for (let j = 0; j < i; j++) {
        f = d[j];
        V[at(j, i)] = f;
        g = e[j] + V[at(j, j)] * f;
        for (let k = j + 1; k <= i - 1; k++) {
          g += V[at(k, j)] * d[k];
          e[k] += V[at(k, j)] * f;
        }
        e[j] = g;
      }
      f = 0;
      for (let j = 0; j < i; j++) {
        e[j] /= h;
        f += e[j] * d[j];
      }
      const hh = f / (h + h);
      for (let j = 0; j < i; j++) e[j] -= hh * d[j];
      for (let j = 0; j < i; j++) {
        f = d[j];
        g = e[j];
        for (let k = j; k <= i - 1; k++) V[at(k, j)] -= f * e[k] + g * d[k];
        d[j] = V[at(i - 1, j)];
        V[at(i, j)] = 0;
      }
    }
    d[i] = h;
  }
  for (let i = 0; i < n - 1; i++) {
    V[at(n - 1, i)] = V[at(i, i)];
    V[at(i, i)] = 1;
    const h = d[i + 1];
    if (h !== 0) {
      for (let k = 0; k <= i; k++) d[k] = V[at(k, i + 1)] / h;
      for (let j = 0; j <= i; j++) {
        let g = 0;
        for (let k = 0; k <= i; k++) g += V[at(k, i + 1)] * V[at(k, j)];
        for (let k = 0; k <= i; k++) V[at(k, j)] -= g * d[k];
      }
    }
    for (let k = 0; k <= i; k++) V[at(k, i + 1)] = 0;
  }
  for (let j = 0; j < n; j++) {
    d[j] = V[at(n - 1, j)];
    V[at(n - 1, j)] = 0;
  }
  V[at(n - 1, n - 1)] = 1;
  e[0] = 0;
  // ---- tql2 ----
  for (let i = 1; i < n; i++) e[i - 1] = e[i];
  e[n - 1] = 0;
  let f = 0;
  let tst1 = 0;
  const eps = 2 ** -52;
  for (let l = 0; l < n; l++) {
    tst1 = Math.max(tst1, Math.abs(d[l]) + Math.abs(e[l]));
    let m = l;
    while (m < n) {
      if (Math.abs(e[m]) <= eps * tst1) break;
      m++;
    }
    if (m > l) {
      let iter = 0;
      do {
        if (++iter > 200) throw new Error('symmetricEigen: QL iteration did not converge');
        let g = d[l];
        let pp = (d[l + 1] - g) / (2 * e[l]);
        let r = Math.hypot(pp, 1);
        if (pp < 0) r = -r;
        d[l] = e[l] / (pp + r);
        d[l + 1] = e[l] * (pp + r);
        const dl1 = d[l + 1];
        let h = g - d[l];
        for (let i = l + 2; i < n; i++) d[i] -= h;
        f += h;
        pp = d[m];
        let c = 1;
        let c2 = c;
        let c3 = c;
        const el1 = e[l + 1];
        let s = 0;
        let s2 = 0;
        for (let i = m - 1; i >= l; i--) {
          c3 = c2;
          c2 = c;
          s2 = s;
          g = c * e[i];
          h = c * pp;
          r = Math.hypot(pp, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = pp / r;
          pp = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (let k = 0; k < n; k++) {
            h = V[at(k, i + 1)];
            V[at(k, i + 1)] = s * V[at(k, i)] + c * h;
            V[at(k, i)] = c * V[at(k, i)] - s * h;
          }
        }
        pp = (-s * s2 * c3 * el1 * e[l]) / dl1;
        e[l] = s * pp;
        d[l] = c * pp;
      } while (Math.abs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }
  for (let i = 0; i < n - 1; i++) {
    let k = i;
    let pp = d[i];
    for (let j = i + 1; j < n; j++) {
      if (d[j] < pp) {
        k = j;
        pp = d[j];
      }
    }
    if (k !== i) {
      d[k] = d[i];
      d[i] = pp;
      for (let j = 0; j < n; j++) {
        const t = V[at(j, i)];
        V[at(j, i)] = V[at(j, k)];
        V[at(j, k)] = t;
      }
    }
  }
}

/**
 * End-gas projector of a GridModeSet chamber: the end gas is the gas farther than the flame radius r from
 * the spark-gap centre (a Euclidean ball, as in the flame geometry), column by column — sub-column q (plan
 * area a_q, gas interval [yLo, yHi], height H_q) holds a_q (H_q − ℓ_q(r)) of end gas, ℓ_q(r) the length of
 * the ball's chord inside its gas interval. The end-gas volume and its modal integrals Σ_q V_eg,q ψ_j(cell q)
 * are therefore continuous in r (fractional coverage: no thin-crescent zero, validation round 2) and use
 * the same capacity-weighted inner product as the modes (r = 0 gives Σ c_i ψ_j,i = 0: a uniform release
 * excites nothing). Shells (KnockOscillator.setSequentialEndGasOutsideRadius) are spherical shells about
 * the plug. Sub-columns are sorted by their nearest distance r_on to the spark, with suffix sums of the
 * volume and the modal integrals: an evaluation visits only the columns the sphere cuts
 * (r − max(r_full − r_on) < r_on < r). The flame ball reaching through the deck metal between the bore and
 * the pocket is the flame-geometry model's own (non-convex chamber) approximation. Immutable.
 */
export class GridEndGasProjector implements EndGasProjector {
  readonly modes: GridModeSet;
  /** Spark-gap centre (x, y, z), cylinder frame, m. */
  readonly spark: readonly [number, number, number];
  /** Chamber volume (the modes' reference depth), m³. */
  readonly totalMeasure: number;
  readonly maxRadius: number;
  private readonly count: number;
  private readonly rOn: Float64Array;
  private readonly rFull: Float64Array;
  private readonly rho2: Float64Array;
  private readonly dLo: Float64Array;
  private readonly dHi: Float64Array;
  private readonly area: Float64Array;
  private readonly cell: Int32Array;
  /** max(r_full − r_on): the radial width of the columns the sphere can cut. */
  private readonly band: number;
  private readonly suffixV: Float64Array;
  private readonly suffixPsi: Float64Array;

  constructor(modes: GridModeSet, spark: readonly [number, number, number]) {
    this.modes = modes;
    this.spark = [spark[0], spark[1], spark[2]];
    const col = modes.columns;
    const nq = col.count;
    const N = modes.nModes;
    const rOn = new Float64Array(nq);
    const rFull = new Float64Array(nq);
    const rho2 = new Float64Array(nq);
    const dLo = new Float64Array(nq);
    const dHi = new Float64Array(nq);
    for (let q = 0; q < nq; q++) {
      const ex = col.x[q] - spark[0];
      const ez = col.z[q] - spark[2];
      rho2[q] = ex * ex + ez * ez;
      dLo[q] = col.yLo[q] - spark[1];
      dHi[q] = col.yHi[q] - spark[1];
      const dMin = dLo[q] > 0 ? dLo[q] : dHi[q] < 0 ? -dHi[q] : 0;
      const dMax = Math.max(Math.abs(dLo[q]), Math.abs(dHi[q]));
      rOn[q] = Math.sqrt(rho2[q] + dMin * dMin);
      rFull[q] = Math.sqrt(rho2[q] + dMax * dMax);
    }
    const order = Array.from({ length: nq }, (_, q) => q).sort((a, b) => rOn[a] - rOn[b] || a - b);
    this.count = nq;
    this.rOn = Float64Array.from(order, (q) => rOn[q]);
    this.rFull = Float64Array.from(order, (q) => rFull[q]);
    this.rho2 = Float64Array.from(order, (q) => rho2[q]);
    this.dLo = Float64Array.from(order, (q) => dLo[q]);
    this.dHi = Float64Array.from(order, (q) => dHi[q]);
    this.area = Float64Array.from(order, (q) => col.area[q]);
    this.cell = Int32Array.from(order, (q) => col.cell[q]);
    let band = 0;
    let rMax = 0;
    for (let q = 0; q < nq; q++) {
      band = Math.max(band, this.rFull[q] - this.rOn[q]);
      rMax = Math.max(rMax, this.rFull[q]);
    }
    this.band = band;
    this.maxRadius = rMax;
    this.suffixV = new Float64Array(nq + 1);
    this.suffixPsi = new Float64Array((nq + 1) * N);
    const psi = modes.psi;
    for (let q = nq - 1; q >= 0; q--) {
      const w = this.area[q] * (this.dHi[q] - this.dLo[q]);
      this.suffixV[q] = this.suffixV[q + 1] + w;
      const c = this.cell[q] * N;
      for (let j = 0; j < N; j++) this.suffixPsi[q * N + j] = this.suffixPsi[(q + 1) * N + j] + w * psi[c + j];
    }
    this.totalMeasure = this.suffixV[0];
  }

  /** First sorted sub-column with r_on ≥ r (strict: r_on > r). */
  private bound(r: number, strict: boolean): number {
    let lo = 0;
    let hi = this.count;
    const a = this.rOn;
    while (lo < hi) {
      const mid = (lo + hi) >> 1;
      if (strict ? a[mid] <= r : a[mid] < r) lo = mid + 1;
      else hi = mid;
    }
    return lo;
  }

  /** End-gas volume of sorted sub-column q (cut by the sphere: r_on < r < r_full). */
  private cutVolume(q: number, r: number): number {
    const w = Math.sqrt(Math.max(0, r * r - this.rho2[q]));
    const lo = this.dLo[q];
    const hi = this.dHi[q];
    const top = hi < w ? hi : w;
    const bot = lo > -w ? lo : -w;
    const chord = top > bot ? top - bot : 0;
    return this.area[q] * (hi - lo - chord);
  }

  outsideMeasure(r: number): number {
    if (!(r > 0)) return this.totalMeasure;
    const k0 = this.bound(r, false);
    let v = this.suffixV[k0];
    for (let q = this.bound(r - this.band, true); q < k0; q++) if (this.rFull[q] > r) v += this.cutVolume(q, r);
    return v;
  }

  outsideIntegrals(r: number, acc: Float64Array): number {
    const N = this.modes.nModes;
    const psi = this.modes.psi;
    const rr = r > 0 ? r : 0;
    const k0 = rr > 0 ? this.bound(rr, false) : 0;
    let v = this.suffixV[k0];
    for (let j = 0; j < N; j++) acc[j] = this.suffixPsi[k0 * N + j];
    if (rr > 0) {
      for (let q = this.bound(rr - this.band, true); q < k0; q++) {
        if (!(this.rFull[q] > rr)) continue;
        const vq = this.cutVolume(q, rr);
        if (vq === 0) continue;
        v += vq;
        const c = this.cell[q] * N;
        for (let j = 0; j < N; j++) acc[j] += vq * psi[c + j];
      }
    }
    return v;
  }

  /**
   * Root of outsideMeasure(r) = target in [lo, hi] (non-increasing in r): Illinois (modified regula
   * falsi) to |Δr| ≤ 1e-12·maxRadius or a residual ≤ 1e-14·totalMeasure.
   */
  radiusForOutsideMeasure(target: number, lo: number, hi: number): number {
    let a = lo;
    let b = hi;
    let fa = this.outsideMeasure(a) - target;
    let fb = this.outsideMeasure(b) - target;
    if (!(fa > 0)) return a;
    if (!(fb < 0)) return b;
    const tolR = 1e-12 * this.maxRadius;
    const tolV = 1e-14 * this.totalMeasure;
    let side = 0;
    for (let it = 0; it < 200 && b - a > tolR; it++) {
      let c = (a * fb - b * fa) / (fb - fa);
      if (!(c > a && c < b)) c = 0.5 * (a + b);
      const fc = this.outsideMeasure(c) - target;
      if (Math.abs(fc) <= tolV) return c;
      if (fc > 0) {
        a = c;
        fa = fc;
        if (side === 1) fb *= 0.5;
        side = 1;
      } else {
        b = c;
        fb = fc;
        if (side === -1) fa *= 0.5;
        side = -1;
      }
    }
    return 0.5 * (a + b);
  }

  /** Radius of the ball about the plug with the end-gas volume fraction f (0..1) outside it. */
  radiusForOutsideFraction(f: number): number {
    if (!(f < 1)) return 0;
    if (!(f > 0)) return this.maxRadius;
    return this.radiusForOutsideMeasure(f * this.totalMeasure, 0, this.maxRadius);
  }
}

// ---------------------------------------------------------------------------------------
// Engine factory
// ---------------------------------------------------------------------------------------

/** Knock acoustics of an L-head engine: modes + end-gas projector, built once and shared by every cylinder. */
export interface LHeadKnockAcoustics {
  readonly modes: GridModeSet;
  readonly projector: GridEndGasProjector;
  /** Bore-column depth (roof → crown) the modes were computed at, m. */
  readonly depth: number;
}

/** Options of lHeadKnockAcoustics. */
export interface LHeadKnockAcousticsOptions extends GridModeSetOptions {
  /**
   * Bore-column depth (roof of the head cavity → piston crown) at which the modes and the end-gas
   * geometry are evaluated, m (default TDC, lHeadDepthAtTDC). Knock occurs near TDC; the Model T
   * fundamental moves +1 % from TDC to 15° ATDC (knock-lhead.test.ts).
   */
  depth?: number;
}

const L_HEAD_ACOUSTICS_CACHE = new Map<string, LHeadKnockAcoustics>();
const L_HEAD_ACOUSTICS_CACHE_SIZE = 8;

/**
 * The L-head knock acoustics of `spec` (geometry.chamber 'l-head'): GridModeSet of the bore + valve-pocket
 * planform (lHeadPlanform) at the given bore-column depth and the GridEndGasProjector about the spark gap.
 * Memoised by geometry (≈ 0.1–0.3 s to build; the instances are immutable and shared).
 */
export function lHeadKnockAcoustics(spec: EngineSpec, opts: LHeadKnockAcousticsOptions = {}): LHeadKnockAcoustics {
  const g = spec.geometry;
  if (g.chamber !== 'l-head' || !g.lHead) throw new Error(`${spec.name}: lHeadKnockAcoustics needs geometry.chamber 'l-head' with geometry.lHead`);
  const depth = opts.depth ?? lHeadDepthAtTDC(g.lHead);
  const key = JSON.stringify([
    g.bore,
    g.lHead,
    spec.sparkPlug.gapCenter,
    depth,
    opts.cellSize ?? null,
    opts.maxAlpha ?? null,
    opts.dampingReferenceAlpha ?? null,
    opts.blockSize ?? null,
    opts.tolerance ?? null,
  ]);
  const hit = L_HEAD_ACOUSTICS_CACHE.get(key);
  if (hit) return hit;
  const modes = new GridModeSet(g.bore, lHeadPlanform(g.bore, g.lHead, depth), opts);
  const projector = new GridEndGasProjector(modes, spec.sparkPlug.gapCenter);
  const ac: LHeadKnockAcoustics = { modes, projector, depth };
  if (L_HEAD_ACOUSTICS_CACHE.size >= L_HEAD_ACOUSTICS_CACHE_SIZE) {
    const oldest = L_HEAD_ACOUSTICS_CACHE.keys().next().value;
    if (oldest !== undefined) L_HEAD_ACOUSTICS_CACHE.delete(oldest);
  }
  L_HEAD_ACOUSTICS_CACHE.set(key, ac);
  return ac;
}

/**
 * MAPO measurement band for L-head engines, Hz. UNVERIFIED convention: no Model T knock measurement or
 * pickup exists, so the band only has to contain the chamber's modes up to the CFR convention's upper edge
 * (18 kHz, Hoth 2021 / SON 2023 — which is also about the first axial resonance c/(2h) ≈ 18.7 kHz of the
 * Model T bore column at TDC, where the depth-averaged modes lose validity). The lower edge sits below the
 * Model T fundamental (α ≈ 1.20: 3.2–4.0 kHz for c = 800–1000 m/s; 2 kHz needs c < 500 m/s) — the
 * CFR's 4 kHz would cut it off.
 */
export const L_HEAD_MAPO_BAND: readonly [number, number] = Object.freeze([2000, 18000]) as readonly [number, number];

/**
 * Virtual knock sensor of an engine without a pickup (Model T): the spark plug's planform position (x, z),
 * cylinder frame, m — a plug-mounted transducer, the usual retrofit for in-cylinder pressure on old heads.
 */
export function virtualKnockSensor(spec: EngineSpec): readonly [number, number] {
  const g = spec.sparkPlug.gapCenter;
  return [g[0], g[2]];
}

/** Options of createKnockOscillator: KnockOscillatorOptions (band null = unfiltered) + the L-head acoustics options. */
export type KnockFactoryOptions = Omit<KnockOscillatorOptions, 'band'> &
  LHeadKnockAcousticsOptions & {
    /** MAPO band, Hz; null = unfiltered; absent: none (flat disc) or L_HEAD_MAPO_BAND (L-head). */
    band?: readonly [number, number] | null;
  };

/**
 * The knock oscillator of one cylinder of `spec`, dispatching on geometry.chamber:
 *  - 'flat-disc': the Bessel modes of the bore and the exact disc projector — identical to
 *    `new KnockOscillator(bore, opts)` (the CFR), plus the spark gap for setEndGasVolumeFraction;
 *  - 'l-head': the shared lHeadKnockAcoustics(spec, opts) (numerical depth-weighted modes of the
 *    bore + valve-pocket planform, end gas outside the flame ball about the plug), with the per-engine
 *    defaults sensor = virtualKnockSensor(spec) (at the plug) and band = L_HEAD_MAPO_BAND.
 * The cycle model calls setEndGasVolumeFraction(f, K) at the knock onset for either chamber.
 */
export function createKnockOscillator(spec: EngineSpec, opts: KnockFactoryOptions = {}): KnockOscillator {
  const g = spec.geometry;
  const spark = spec.sparkPlug.gapCenter;
  const { band, ...rest } = opts;
  if (g.chamber === 'l-head') {
    const ac = rest.modes && rest.projector ? null : lHeadKnockAcoustics(spec, opts);
    return new KnockOscillator(g.bore, {
      ...rest,
      sensor: rest.sensor ?? virtualKnockSensor(spec),
      band: band === null ? undefined : band ?? L_HEAD_MAPO_BAND,
      modes: rest.modes ?? ac!.modes,
      projector: rest.projector ?? ac!.projector,
      spark: rest.spark ?? spark,
    });
  }
  return new KnockOscillator(g.bore, { ...rest, band: band ?? undefined, spark: rest.spark ?? spark });
}
