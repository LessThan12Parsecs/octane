/**
 * Quasi-steady, one-dimensional, isentropic compressible flow through a restriction
 * (valve, throttle, orifice) with an effective area C_D·A_R.
 *
 * ── Model ───────────────────────────────────────────────────────────────────────────────
 * Upstream stagnation state (p₀, T₀, R₀, γ₀), downstream static pressure p_d, pressure
 * ratio r = p_d/p₀. Frozen composition, constant γ (evaluated by the caller at T₀),
 * isentropic expansion to the throat, throat static pressure = max(p_d, p*) — the
 * standard restriction-flow model of Heywood (1988), App. C, eqs. (C.8)–(C.9):
 *
 *   ṁ = C_D A_R p₀/√(R₀T₀) · Ψ(r)
 *   Ψ(r) = √( 2γ/(γ−1) · (r^{2/γ} − r^{(γ+1)/γ}) )            r ≥ r*   (subsonic)
 *   Ψ*   = √γ · (2/(γ+1))^{(γ+1)/(2(γ−1))}                     r < r*   (choked)
 *   r*   = (2/(γ+1))^{γ/(γ−1)}   (0.5283 for γ = 1.4)
 *
 * The choked form and r* are verified against MIT OCW 2.61 (Spring 2017) lecture 8,
 * "Choking effect" (fetched); the subsonic form is the textbook isentropic nozzle and is
 * checked numerically (tools/reference/gasex_orifice.py): Ψ is maximal exactly at r*,
 * where it equals Ψ*, so the two branches join with a continuous first derivative.
 * UNVERIFIED: the Heywood appendix equation numbers (C.8)/(C.9) (quoted from memory).
 * Frozen-γ(T₀) error vs a variable-cp isentropic nozzle (Cantera, same species data;
 * tools/reference/gasex_orifice.py, r = 0.2 … 0.9999): air 300 K 0.017 %; PRF/ethanol fuel-vapour
 * charges 302–422 K 0.11–0.15 %; burned gas 1200–2200 K 0.08–0.15 % (jet speed ≤ 0.25 %).
 *
 * ── Low-Δp regularisation ───────────────────────────────────────────────────────────────
 * With x = 1 − r = Δp/p₀ (fractional pressure drop) the exact law is Ψ = G(x)·√x with
 *   G(x)² = 2γ/(γ−1) · (1−x)^{2/γ} · E(x),   E(x) = (1 − (1−x)^k)/x,   k = (γ−1)/γ,
 * a smooth (analytic) function with G(0) = √2 (incompressible limit ṁ = C_D A √(2ρΔp)).
 * √x has an infinite slope at x = 0, which makes the cylinder/plenum ODE infinitely stiff
 * as |Δp| → 0 (e.g. near every flow reversal). Inside the band x < δ we replace √x by
 *   √δ · P(x/δ),   P(t) = (45 t − 18 t³ + 5 t⁵)/32   (= 1.40625 t − 0.5625 t³ + 0.15625 t⁵),
 * the odd quintic that matches √t in value, slope and curvature at t = 1. Consequences:
 *  - EXACT (bit-for-bit the unregularised law) for x ≥ δ; C² at x = δ; C^∞ (odd) through
 *    x = 0 for a given upstream state; finite slope ∂ṁ/∂Δp = C_D A √(2/(R₀T₀)) · 1.40625/√δ
 *    at Δp = 0 (use it to bound an explicit integrator's step: the pressure-equalisation
 *    time constant of a volume V is τ ≈ ρV/(γ p ∂ṁ/∂Δp)).
 *  - 0 ≤ √t − P(t) ≤ 0.17889 on [0, 1] (max at t ≈ 0.127), so inside the band the
 *    regularised flow UNDER-estimates |ṁ| by at most 0.1789 · ṁ_exact(x = δ)·(1 + O(δ)).
 *  - The bidirectional form below is continuous through Δp = 0 (ṁ = 0 on both sides) but
 *    its slope jumps when the two sides have different √(RT) (Lipschitz, harmless for RK).
 * ORIFICE_REG_DEFAULT = 1e-3 (band |Δp| < 0.1 % of p₀ ≈ 100 Pa at 1 atm).
 *
 * Hot paths allocate nothing.
 */

/** Default regularisation half-width δ: fraction of the upstream stagnation pressure. */
export const ORIFICE_REG_DEFAULT = 1e-3;

/** Largest δ accepted (the band must stay far below the critical drop x* ≈ 0.47). */
const REG_MAX = 0.05;

// P(t) = A1 t + A3 t³ + A5 t⁵: P(1) = 1, P'(1) = 1/2, P''(1) = −1/4.
const A1 = 45 / 32;
const A3 = -18 / 32;
const A5 = 5 / 32;

/** max over t∈[0,1] of √t − P(t): bound of the in-band under-estimate (see header). */
export const ORIFICE_REG_MAX_DEFICIT = 0.178886;

/** Critical (sonic-throat) pressure ratio r* = (2/(γ+1))^{γ/(γ−1)}. */
export function criticalPressureRatio(gamma: number): number {
  return Math.pow(2 / (gamma + 1), gamma / (gamma - 1));
}

/** Choked mass-flux function Ψ* = √γ (2/(γ+1))^{(γ+1)/(2(γ−1))} (dimensionless). */
export function chokedFluxFunction(gamma: number): number {
  return Math.sqrt(gamma) * Math.pow(2 / (gamma + 1), (gamma + 1) / (2 * (gamma - 1)));
}

/**
 * Exact (unregularised) mass-flux function Ψ(r) for r = p_d/p₀ (dimensionless):
 * ṁ = C_D A p₀/√(R₀T₀) · Ψ. Returns 0 for r ≥ 1 and Ψ* for r ≤ r*.
 */
export function massFluxFunction(r: number, gamma: number): number {
  if (!(r < 1)) return 0;
  const rc = criticalPressureRatio(gamma);
  if (r <= rc) return chokedFluxFunction(gamma);
  const x = 1 - r;
  return gFactor(x, gamma) * Math.sqrt(x);
}

// E(x) = (1 − (1−x)^k)/x, accurate for all 0 < x < 1 (expm1/log1p), E(0) = k.
function eFactor(x: number, k: number): number {
  if (x < 1e-12) return k * (1 + 0.5 * (1 - k) * x);
  return -Math.expm1(k * Math.log1p(-x)) / x;
}

// dE/dx; Taylor series for small x (cancellation), direct formula otherwise.
function eFactorDeriv(x: number, k: number): number {
  if (x < 1e-3) {
    const c1 = 0.5 * k * (1 - k);
    return c1 + (2 / 3) * c1 * (2 - k) * x + (1 / 8) * k * (1 - k) * (2 - k) * (3 - k) * x * x;
  }
  return (k * Math.pow(1 - x, k - 1) - eFactor(x, k)) / x;
}

/** G(x) = Ψ(1 − x)/√x (smooth, G(0) = √2). */
function gFactor(x: number, gamma: number): number {
  const k = (gamma - 1) / gamma;
  return Math.sqrt(((2 * gamma) / (gamma - 1)) * Math.pow(1 - x, 2 / gamma) * eFactor(x, k));
}

/**
 * Scratch output of {@link fluxFunctionReg}: [0] Ψ, [1] dΨ/dx, [2] choked (1/0),
 * [3] 1 − r_t^{(γ−1)/γ} at the throat (for the jet speed).
 * Module-level (not re-entrant; single worker thread).
 */
const FLUX = new Float64Array(4);

// Two-slot memo of the γ-only constants (the two sides of a restriction alternate), in a typed
// array: module-level `let` doubles are boxed as heap numbers on every write (≈ 40 B per call;
// validation round 1 — hot paths must not allocate).
// [0] γ_A, [1] x_c,A, [2] Ψ*_A, [3] γ_B, [4] x_c,B, [5] Ψ*_B, [6] current x_c, [7] current Ψ*
const GM = Float64Array.of(NaN, 0, 0, NaN, 0, 0, 0, 0);
function gammaConstants(gamma: number): void {
  const m = GM;
  if (gamma === m[0]) {
    m[6] = m[1];
    m[7] = m[2];
    return;
  }
  if (gamma === m[3]) {
    m[6] = m[4];
    m[7] = m[5];
    return;
  }
  m[3] = m[0];
  m[4] = m[1];
  m[5] = m[2];
  m[0] = gamma;
  m[1] = 1 - criticalPressureRatio(gamma);
  m[2] = chokedFluxFunction(gamma);
  m[6] = m[1];
  m[7] = m[2];
}

/**
 * Regularised mass-flux function of the fractional pressure drop x = 1 − p_d/p₀ ≥ 0,
 * and its derivative dΨ/dx, written to FLUX. x ≤ 0 → 0.
 */
function fluxFunctionReg(x: number, gamma: number, delta: number): void {
  if (!(x > 0)) {
    // x = 0: Ψ = 0 with the finite regularised slope; x < 0 is not this direction.
    FLUX[0] = 0;
    FLUX[1] = x === 0 ? (Math.SQRT2 * A1) / Math.sqrt(delta) : 0;
    FLUX[2] = 0;
    FLUX[3] = 0;
    return;
  }
  gammaConstants(gamma);
  if (x >= GM[6]) {
    FLUX[0] = GM[7];
    FLUX[1] = 0;
    FLUX[2] = 1;
    FLUX[3] = (gamma - 1) / (gamma + 1); // 1 − r*^k = 1 − 2/(γ+1)
    return;
  }
  const k = (gamma - 1) / gamma;
  const r = 1 - x;
  const pre = (2 * gamma) / (gamma - 1);
  const lr = Math.log1p(-x); // ln r, accurate for small x
  const omk = x < 1e-12 ? k * x * (1 + 0.5 * (1 - k) * x) : -Math.expm1(k * lr); // 1 − r^k
  const E = omk / x;
  const r2g = Math.exp((2 / gamma) * lr); // r^{2/γ}
  const g2 = pre * r2g * E;
  const G = Math.sqrt(g2);
  FLUX[2] = 0;
  FLUX[3] = omk;
  if (x >= delta) {
    const psi = Math.sqrt(pre * r2g * omk); // = G √x
    FLUX[0] = psi;
    // dΨ/dx = −dΨ/dr = −(2 r^{2/γ−1} − (γ+1) r^{1/γ}) / ((γ−1) Ψ)
    FLUX[1] = -(2 * (r2g / r) - (gamma + 1) * Math.sqrt(r2g)) / ((gamma - 1) * psi);
    return;
  }
  // In the band: Ψ = G(x)·S(x), S = √δ P(x/δ).
  const sd = Math.sqrt(delta);
  const t = x / delta;
  const t2 = t * t;
  const S = sd * t * (A1 + t2 * (A3 + t2 * A5));
  const dS = (A1 + t2 * (3 * A3 + 5 * A5 * t2)) / sd;
  // d(G²)/dx = pre [−(2/γ) r^{2/γ−1} E + r^{2/γ} E']
  const dg2 = pre * (-(2 / gamma) * (r2g / r) * E + r2g * eFactorDeriv(x, k));
  const dG = dg2 / (2 * G);
  FLUX[0] = G * S;
  FLUX[1] = dG * S + G * dS;
}

/**
 * Regularised mass-flux function Ψ_reg(x) of the fractional pressure drop x = Δp/p₀
 * (dimensionless; 0 for x ≤ 0). Exact for x ≥ δ. Exposed for tests / analysis.
 */
export function regularizedFluxFunction(x: number, gamma: number, delta = ORIFICE_REG_DEFAULT): number {
  fluxFunctionReg(x, gamma, clampDelta(delta));
  return FLUX[0];
}

function clampDelta(delta: number): number {
  return delta > REG_MAX ? REG_MAX : delta > 1e-12 ? delta : 1e-12;
}

/**
 * One-directional restriction mass flow, kg/s (≥ 0): isentropic nozzle from the upstream
 * stagnation state to the downstream static pressure, subsonic or choked, with the
 * low-Δp regularisation of width `delta` (fraction of p₀; see file header).
 * Returns 0 when pDown ≥ p0, CdA ≤ 0 or p0 ≤ 0.
 *
 * @param CdA effective area C_D·A_R, m²
 * @param p0 upstream stagnation pressure, Pa
 * @param T0 upstream stagnation temperature, K
 * @param R0 upstream specific gas constant, J/(kg K)
 * @param gamma0 upstream ratio of specific heats (frozen, at T0)
 * @param pDown downstream static pressure, Pa
 */
export function orificeMassFlow(
  CdA: number,
  p0: number,
  T0: number,
  R0: number,
  gamma0: number,
  pDown: number,
  delta = ORIFICE_REG_DEFAULT,
): number {
  if (!(CdA > 0 && p0 > 0) || !(pDown < p0)) return 0;
  const x = pDown > 0 ? 1 - pDown / p0 : 1;
  fluxFunctionReg(x, gamma0, clampDelta(delta));
  return (CdA * p0 * FLUX[0]) / Math.sqrt(R0 * T0);
}

/**
 * Isentropic throat (vena-contracta) jet speed, m/s, for an upstream stagnation state and a
 * pressure ratio r = p_d/p₀ (throat pressure = max(p_d, p*)):
 *   v = √( 2 γ R T₀/(γ−1) · (1 − r_t^{(γ−1)/γ}) ).
 * Not regularised (∝ √Δp near Δp = 0); finite everywhere, 0 for r ≥ 1.
 */
export function isentropicJetVelocity(T0: number, R0: number, gamma0: number, r: number): number {
  if (!(r < 1)) return 0;
  const rc = criticalPressureRatio(gamma0);
  const rt = r > rc ? r : rc;
  const k = (gamma0 - 1) / gamma0;
  return Math.sqrt(((2 * R0 * T0) / k) * -Math.expm1(k * Math.log(rt)));
}

/** Result of {@link orificeFlow} (reuse one instance; allocation-free). */
export interface OrificeFlow {
  /** Signed mass flow, kg/s: positive from side a to side b. */
  mdot: number;
  /** Which side's stagnation state feeds the flow: 0 = a (mdot ≥ 0), 1 = b (mdot < 0). */
  upstream: 0 | 1;
  /** True when the throat is sonic. */
  choked: boolean;
  /** Fractional pressure drop x = |Δp|/p_up (dimensionless). */
  pressureDrop: number;
  /** Isentropic throat jet speed (unsigned), m/s. */
  jetVelocity: number;
  /** ∂mdot/∂p_a at fixed upstream T, R, γ, kg/(s Pa). */
  dmdotdpa: number;
  /** ∂mdot/∂p_b at fixed upstream T, R, γ, kg/(s Pa). */
  dmdotdpb: number;
}

/** Allocate an {@link OrificeFlow} result holder. */
export function newOrificeFlow(): OrificeFlow {
  return { mdot: 0, upstream: 0, choked: false, pressureDrop: 0, jetVelocity: 0, dmdotdpa: 0, dmdotdpb: 0 };
}

/**
 * Bidirectional restriction flow between two gas volumes a and b (each at rest, so static =
 * stagnation). The side with the higher pressure is upstream; its (T, R, γ) set the flow.
 * The caller carries the enthalpy/composition flux with the upstream side's state:
 * Ḣ = mdot·h₀(upstream). Writes into `out` and returns out.mdot (kg/s, + = a → b).
 *
 * @param CdA effective area, m² (≤ 0 → no flow)
 * @param pa, Ta, Ra, gammaA side-a pressure (Pa), temperature (K), gas constant (J/(kg K)), γ
 * @param pb, Tb, Rb, gammaB side-b state
 * @param delta regularisation half-width (fraction of the upstream pressure)
 */
export function orificeFlow(
  CdA: number,
  pa: number,
  Ta: number,
  Ra: number,
  gammaA: number,
  pb: number,
  Tb: number,
  Rb: number,
  gammaB: number,
  out: OrificeFlow,
  delta = ORIFICE_REG_DEFAULT,
): number {
  const d = clampDelta(delta);
  const aUp = pa >= pb;
  const pu = aUp ? pa : pb;
  const pd = aUp ? pb : pa;
  const Tu = aUp ? Ta : Tb;
  const Ru = aUp ? Ra : Rb;
  const gu = aUp ? gammaA : gammaB;
  out.upstream = aUp ? 0 : 1;
  if (!(CdA > 0 && pu > 0)) {
    out.mdot = 0;
    out.choked = false;
    out.pressureDrop = 0;
    out.jetVelocity = 0;
    out.dmdotdpa = 0;
    out.dmdotdpb = 0;
    return 0;
  }
  const x = pd > 0 ? 1 - pd / pu : 1;
  fluxFunctionReg(x, gu, d);
  const psi = FLUX[0];
  const dpsi = FLUX[1];
  const c = CdA / Math.sqrt(Ru * Tu);
  const m = c * pu * psi; // magnitude
  // ∂m/∂p_up = c(Ψ + Ψ' p_d/p_u), ∂m/∂p_down = −c Ψ'
  const dUp = c * (psi + (dpsi * pd) / pu);
  const dDown = -c * dpsi;
  out.choked = FLUX[2] === 1;
  out.pressureDrop = x;
  // v = √(2 R T₀/k · (1 − r_t^k)), k = (γ−1)/γ (same as isentropicJetVelocity, no extra pow)
  out.jetVelocity = Math.sqrt(((2 * gu) / (gu - 1)) * Ru * Tu * FLUX[3]);
  if (aUp) {
    out.mdot = m;
    out.dmdotdpa = dUp;
    out.dmdotdpb = dDown;
  } else {
    out.mdot = -m;
    out.dmdotdpa = -dDown;
    out.dmdotdpb = -dUp;
  }
  return out.mdot;
}
