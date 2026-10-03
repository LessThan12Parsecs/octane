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
  /** Include every mode with α_mn ≤ maxAlpha (default 7.1: (1,0) (2,0) (0,1) (3,0) (4,0) (1,1) (5,0) (2,1) (0,2)). */
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
   * Source-projection quadrature: [Gauss–Legendre nodes per radial sub-interval, unused] (default
   * [24, 0]; the angular integrals are exact — see setEndGasRegion). Round 1: midpoint cells
   * (radial × angular, 48 × 96).
   */
  quadrature?: readonly [number, number];
}

/** Default (1,0) amplitude decay time, s. UNVERIFIED (see file header). */
export const KNOCK_DECAY_TIME = 1.0e-3;

/**
 * Damped transverse acoustic modes of the combustion chamber driven by the end-gas
 * heat-release rate (see file header). Allocation-free step(); one instance per cylinder.
 */
export class KnockOscillator {
  readonly bore: number;
  readonly radius: number;
  readonly nModes: number;
  /** Circumferential order m of mode j. */
  readonly m: Int32Array;
  /** Radial order n of mode j. */
  readonly n: Int32Array;
  /** α_mn of mode j. */
  readonly alpha: Float64Array;
  /** 1 = cos(mθ) orientation, 0 = sin(mθ). */
  readonly cosine: Uint8Array;
  /** Mean-square mode shape over the bore cross-section ⟨ψ_j²⟩, ψ_j = J_m(α_j r/R)·{cos, sin}(mθ) (unnormalised). */
  readonly meanSquare: Float64Array;
  /** Modal source shape ⟨ψ_j⟩_eg / ⟨ψ_j²⟩ for the current end-gas region (dimensionless). */
  readonly sourceShape: Float64Array;
  /** Mode shape at the sensor ψ_j(x_s, z_s). */
  readonly psiSensor: Float64Array;
  /** Modal amplitudes η_j (Pa) and rates η̇_j (Pa/s). */
  readonly eta: Float64Array;
  readonly etaDot: Float64Array;
  /** Relative damping factor (α_j/α_10)^½. */
  private readonly dampScale: Float64Array;
  private readonly decayTime: number;
  private readonly quad: readonly [number, number];
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
  private readonly cosMScratch: Float64Array;
  private readonly sinMScratch: Float64Array;
  private readonly bpScratch = new Float64Array(4);
  /** End-gas area fraction of the planform for the current source (0..1). */
  endGasAreaFraction = 0;

  constructor(bore: number, opts: KnockOscillatorOptions = {}) {
    this.bore = bore;
    this.radius = bore / 2;
    this.decayTime = opts.decayTime ?? KNOCK_DECAY_TIME;
    this.quad = opts.quadrature ?? [24, 0];
    this.bandLo = opts.band ? opts.band[0] : 0;
    this.bandHi = opts.band ? opts.band[1] : Infinity;
    this.bandOrder = opts.bandOrder ?? Infinity;
    const maxAlpha = opts.maxAlpha ?? 7.1;
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
    this.sourceShape = new Float64Array(N);
    this.psiSensor = new Float64Array(N);
    this.eta = new Float64Array(N);
    this.etaDot = new Float64Array(N);
    this.dampScale = new Float64Array(N);
    this.projScratch = new Float64Array(N);
    this.qsPrev = new Float64Array(N);
    this.impScratch = new Float64Array(N);
    this.gainScratch = new Float64Array(N).fill(1);
    this.cosMScratch = new Float64Array(N);
    this.sinMScratch = new Float64Array(N);
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
    const s = opts.sensor ?? [0.9 * this.radius, 0];
    for (let j = 0; j < N; j++) this.psiSensor[j] = this.modeShape(j, s[0], s[1]);
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

  /**
   * Set the end-gas region: the planform outside the flame circle of radius flameRadius (m)
   * centred at (x0, z0) (flame centre projected on the head face). Computes the modal
   * projections sourceShape_j = ⟨ψ_j⟩_eg/⟨ψ_j²⟩ (not a hot path) and returns the end-gas area
   * fraction.
   *
   * Quadrature (validation round 2): the angular integrals over the part of each ring r outside the
   * circle are EXACT — with d = |c|, θ_c = atan2(z0, x0), the ring is outside the circle for
   * |θ − θ_c| > α(r), cos α = (r² + d² − r_f²)/(2 r d), so ∫cos mθ dθ = −2 cos(mθ_c) sin(mα)/m,
   * ∫sin mθ dθ = −2 sin(mθ_c) sin(mα)/m, ∫dθ = 2π − 2α — and the radial integrals use Gauss–Legendre
   * on the sub-intervals between the breakpoints |d − r_f| and d + r_f, with r = a + (b − a)(3u² − 2u³)
   * on the partial-ring intervals (the arc length grows like √(r − a) at a breakpoint; the substitution
   * makes the integrand smooth). Round 1 used a 48 × 96 midpoint grid: an end-gas crescent thinner
   * than the outermost cell row (≈ 0.43 mm) contained no midpoint, so a detected knock rang with
   * MAPO = 0 exactly, and MAPO was biased by 1–4 % elsewhere. `quadrature[0]` (default 24) is the
   * number of Gauss–Legendre nodes per radial sub-interval.
   */
  setEndGasRegion(x0: number, z0: number, flameRadius: number): number {
    const N = this.nModes;
    const acc = this.projScratch;
    const aEg = this.outsideIntegrals(x0, z0, flameRadius, acc);
    const aCyl = Math.PI * this.radius * this.radius;
    this.endGasAreaFraction = aEg / aCyl;
    for (let j = 0; j < N; j++) this.sourceShape[j] = aEg > 0 ? acc[j] / aEg / this.meanSquare[j] : 0;
    this.shellCount = 0;
    this.currentShell = -1;
    return this.endGasAreaFraction;
  }

  /**
   * Sequential autoignition source (validation round 2): the end-gas region outside the flame circle
   * (radius rInner about (x0, z0)) split into K equal-area shells ordered by distance from the flame
   * centre, shell 0 the farthest (the end gas at the chamber periphery opposite the plug, which the
   * flame reaches last). The burn-up releases its heat shell after shell (setSequentialShell), i.e. an
   * autoignition front sweeps the end gas toward the flame over the burn-up time τ_ab — the Zeldovich
   * (1980) spontaneous-ignition front of a stratified end gas (file header). Round 1 released it
   * uniformly over the whole region, whose projection on every transverse mode vanishes as the region
   * approaches the whole bore: the heaviest knock (near-homogeneous autoignition of low-ON fuels) rang
   * with MAPO → 0. UNVERIFIED: the sweep direction (outside-in); real end gases autoignite from
   * several centres (König & Sheppard 1990, SAE 902135, abstract), any non-uniform sequence excites
   * the modes. Sets sourceShape to shell 0 and returns the end-gas area fraction.
   */
  setSequentialEndGasRegion(x0: number, z0: number, rInner: number, K = 16): number {
    const N = this.nModes;
    const R = this.radius;
    const d = Math.hypot(x0, z0);
    const aCyl = Math.PI * R * R;
    const outsideArea = (r: number): number => aCyl - circleLensArea(r, R, d);
    const aEg = outsideArea(rInner);
    this.endGasAreaFraction = aEg / aCyl;
    if (!(aEg > 0) || K < 1) {
      this.setEndGasRegion(x0, z0, rInner);
      return this.endGasAreaFraction;
    }
    if (this.shellShapes.length < K * N) this.shellShapes = new Float64Array(K * N);
    const prev = new Float64Array(N);
    const cur = new Float64Array(N);
    let aPrev = 0;
    let sPrev = d + R; // nothing outside
    for (let k = 0; k < K; k++) {
      // radius s_{k+1} with (k+1)/K of the end-gas area outside it (bisection on the lens area)
      let s: number;
      if (k === K - 1) s = rInner;
      else {
        const target = ((k + 1) / K) * aEg;
        let lo = rInner;
        let hi = sPrev;
        for (let it = 0; it < 60; it++) {
          const mid = 0.5 * (lo + hi);
          if (outsideArea(mid) > target) lo = mid;
          else hi = mid;
        }
        s = 0.5 * (lo + hi);
      }
      const a = this.outsideIntegrals(x0, z0, s, cur);
      const da = a - aPrev;
      for (let j = 0; j < N; j++) this.shellShapes[k * N + j] = da > 0 ? (cur[j] - prev[j]) / da / this.meanSquare[j] : 0;
      prev.set(cur);
      aPrev = a;
      sPrev = s;
    }
    this.shellCount = K;
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

  /** Make shell k (0 = first to autoignite) of setSequentialEndGasRegion the current source shape. */
  setSequentialShell(k: number): void {
    if (this.shellCount === 0) return;
    const kk = k < 0 ? 0 : k >= this.shellCount ? this.shellCount - 1 : k;
    if (kk === this.currentShell) return;
    const N = this.nModes;
    for (let j = 0; j < N; j++) this.sourceShape[j] = this.shellShapes[kk * N + j];
    this.currentShell = kk;
  }

  /**
   * ∫ over the planform outside the circle (radius rf about (x0, z0)) of J_m(α r/R)·trig(mθ) per mode
   * into acc (area-weighted, unnormalised); returns the region's area, m². Quadrature: see setEndGasRegion.
   */
  private outsideIntegrals(x0: number, z0: number, flameRadius: number, acc: Float64Array): number {
    const R = this.radius;
    const N = this.nModes;
    acc.fill(0);
    const gl = gaussLegendre(Math.max(4, Math.round(this.quad[0])));
    const rf = Math.max(0, flameRadius);
    const d = Math.hypot(x0, z0);
    const thc = Math.atan2(z0, x0);
    const cm = this.cosMScratch;
    const sm = this.sinMScratch;
    for (let j = 0; j < N; j++) {
      cm[j] = Math.cos(this.m[j] * thc);
      sm[j] = Math.sin(this.m[j] * thc);
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
          const m = this.m[j];
          let ang: number;
          if (m === 0) ang = 2 * Math.PI - 2 * alpha;
          else if (!partial) ang = 0;
          else ang = (-2 * (this.cosine[j] ? cm[j] : sm[j]) * Math.sin(m * alpha)) / m;
          if (ang === 0) continue;
          acc[j] += wr * ang * besselJ(m, (this.alpha[j] * r) / R);
        }
      }
    }
    return aEg > 0 ? aEg : 0;
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
