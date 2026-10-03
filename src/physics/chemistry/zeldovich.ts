/**
 * Thermal NO: the extended Zeldovich mechanism (Zeldovich 1946; Lavoie, Heywood & Keck 1970).
 *
 * Reactions, all written in the EXOTHERMIC direction (forward = as written):
 *   R1  N + NO  ⇌ N2 + O
 *   R2  N + O2  ⇌ NO + O
 *   R3  N + OH  ⇌ NO + H
 * Reverse rate constants follow from detailed balance with the project's NASA-7 thermo
 * (Kc = Kp because Δn = 0 for all three), so the kinetics is consistent with the burned-gas
 * equilibrium computed from the same data.
 *
 * Two forms are provided:
 *  - `fullRates` — the elementary forward/reverse rates for ARBITRARY concentrations of
 *    N2, O2, O, N, NO, OH, H (d[NO]/dt and d[N]/dt);
 *  - the rate-controlled form of Heywood (1988), Internal Combustion Engine Fundamentals,
 *    eq. (11.8): with O, OH, H, O2, N2 at their equilibrium values and N in quasi-steady state
 *      d[NO]/dt = 2 R1 (1 − α²) / (1 + α R1/(R2 + R3)),   α = [NO]/[NO]e,
 *    R1 = k1[N]e[NO]e, R2 = k2[N]e[O2]e, R3 = k3[N]e[OH]e (one-way equilibrium rates; with the
 *    exothermic-direction constants used here these need no equilibrium constants at all).
 *    This is exact for the three-reaction system under those two assumptions (verified against
 *    Cantera to ~1e-10, see zeldovich.test.ts). The formula is confirmed in a fetched secondary
 *    source quoting Heywood (Chindaprasert, Dissertation Univ. Rostock 2008, eq. 3.57);
 *    UNVERIFIED: the equation number 11.8 in Heywood (1988) itself.
 *
 * Units: T in K, p in Pa, mole fractions, concentrations [·] in mol/m³, rates in mol/(m³ s)
 * (the chemical source per unit volume: dN_NO/dt = ω_NO · V_burned).
 *
 * Scope: thermal NO only. Compared with full GRI-Mech 3.0 from NO-free equilibrium burned gas
 * at constant T, p (2000–2800 K, 20–80 bar, φ 0.8–1.2), this model forms NO 1.2–1.8× more
 * slowly in the first ~τ_NO (the N2O-intermediate and NNH routes are missing; largest at
 * φ = 1.2 and at 2000 K) and reaches α = 0.93–0.95 at 3 τ_NO where GRI reaches 0.90–0.99
 * (test/fixtures/chemistry_zeldovich.json → evolution). Prompt (Fenimore) NO in the flame
 * front and NO2 are not modelled.
 * UNVERIFIED (from memory, not fetched): the bibliographic details of Lavoie, Heywood & Keck
 * (1970, Combust. Sci. Technol. 1:313) and the book/chapter of Hanson & Salimian (1984); the
 * title "Survey of rate constants in the N/H/O system" was confirmed (Semantic Scholar).
 */
import { R_UNIVERSAL } from '../core/constants';
import { SP } from '../core/species';
import { speciesHRT, speciesS0R } from '../thermo/thermo';

/**
 * Rate constants k_i = A_i T^b_i exp(−Ta_i / T) of R1..R3 in the exothermic direction.
 * A in m³/(mol s) (with T in K), Ta = E_a/R in K.
 */
export interface ZeldovichRateSet {
  readonly label: string;
  readonly A: readonly [number, number, number];
  readonly b: readonly [number, number, number];
  readonly Ta: readonly [number, number, number];
}

/** cal → J (thermochemical calorie, as Cantera). */
const CAL = 4.184;
/** cm³/(mol s) → m³/(mol s). */
const CM3 = 1e-6;

/**
 * DEFAULT set — GRI-Mech 3.0 (Smith, Golden, Frenklach, Moriarty, Eiteneer, Goldenberg,
 * Bowman, Hanson, Song, Gardiner, Lissianski & Qin 1999, http://combustion.berkeley.edu/gri-mech/),
 * reactions 178–180 as bundled in Cantera 3.2 gri30.yaml (values asserted by
 * tools/reference/chemistry_mechanisms.py):
 *   R178  N + NO ⇌ N2 + O   A = 2.70e13 cm³/(mol s), b = 0,   Ea = 355 cal/mol
 *   R179  N + O2 ⇌ NO + O   A = 9.00e9,              b = 1.0, Ea = 6500 cal/mol
 *   R180  N + OH ⇌ NO + H   A = 3.36e13,             b = 0,   Ea = 385 cal/mol
 * Chosen because (i) it is a critically evaluated, optimised set, (ii) it is what the
 * Cantera oracle (and full GRI-3.0 comparisons) use.
 */
export const ZELDOVICH_GRI30: ZeldovichRateSet = Object.freeze({
  label: 'GRI-Mech 3.0 (R178–R180)',
  A: [2.7e13 * CM3, 9.0e9 * CM3, 3.36e13 * CM3] as const,
  b: [0, 1, 0] as const,
  Ta: [(355 * CAL) / R_UNIVERSAL, (6500 * CAL) / R_UNIVERSAL, (385 * CAL) / R_UNIVERSAL] as const,
});

/**
 * ALTERNATIVE set — Heywood (1988) Table 11.1, from Hanson & Salimian (1984), "Survey of rate
 * constants in the N/H/O system", in Gardiner (ed.), Combustion Chemistry, Springer, ch. 6.
 * Heywood lists both directions (cm³/(mol s), T in K):
 *   k1+ (O + N2 → NO + N) = 7.6e13 exp(−38000/T)      k1− (N + NO → N2 + O) = 1.6e13
 *   k2+ (N + O2 → NO + O) = 6.4e9 T exp(−3150/T)       k2− (O + NO → O2 + N) = 1.5e9 T exp(−19500/T)
 *   k3+ (N + OH → NO + H) = 4.1e13                     k3− (H + NO → OH + N) = 2.0e14 exp(−23650/T)
 * k1+ = 7.6e13 exp(−38000/T) confirmed in W. Cheng, MIT 2.61 (2017) lectures 11–12 ("See
 * table 11.1 for rates"); k1+, k2−, k3− pre-factors/activation temperatures confirmed in
 * Chindaprasert, Dissertation Univ. Rostock (2008), eq. (3.57) ff.
 * UNVERIFIED: k1− = 1.6e13, k2+ = 6.4e9 T exp(−3150/T), k3+ = 4.1e13 and the T factor of k2−
 * (from memory; the web-search budget of the 2026-09-30 review was exhausted, so still not
 * fetched). Detailed-balance check with our NASA thermo (k+/k− vs Kc, 1500–3000 K): pair 1 agrees
 * to 4–5 %, pair 2 to 4–15 %, but pair 3 is inconsistent by a factor 2.1–2.3 (4.1e13 /
 * (2.0e14 e^{−23650/T}) = 0.44–0.48 Kc3) — the two directions of R3 were fitted independently
 * over different T ranges (300–2500 K vs 2200–4500 K). Only the exothermic-direction constants
 * are used here; reverses follow from our thermo, so this set's k3− is ≈ 2.2× Heywood's.
 * Relative to GRI-3.0 this set gives 0.63–0.82 × the NO formation rate over 1800–2800 K
 * (test/fixtures/chemistry_zeldovich.json summary).
 */
export const ZELDOVICH_HEYWOOD: ZeldovichRateSet = Object.freeze({
  label: 'Heywood (1988) Table 11.1 / Hanson & Salimian (1984)',
  A: [1.6e13 * CM3, 6.4e9 * CM3, 4.1e13 * CM3] as const,
  b: [0, 1, 0] as const,
  Ta: [0, 3150, 0] as const,
});

const iN2 = SP.N2;
const iO2 = SP.O2;
const iO = SP.O;
const iN = SP.N;
const iNO = SP.NO;
const iOH = SP.OH;
const iH = SP.H;

/** g°/(RT) of species k at T (dimensionless, standard state P_REF_THERMO). */
const g0RT = (k: number, T: number): number => speciesHRT(k, T) - speciesS0R(k, T);

/**
 * Extended-Zeldovich kinetics with cached rate constants. One instance per thread of use
 * (holds scratch state; methods allocate nothing).
 */
export class ZeldovichKinetics {
  readonly rates: ZeldovichRateSet;
  /** Forward (exothermic-direction) rate constants at the cached T, m³/(mol s). */
  readonly kf = new Float64Array(3);
  /** Reverse rate constants at the cached T, m³/(mol s) (valid after setTemperatureWithReverse). */
  readonly kr = new Float64Array(3);
  /** One-way equilibrium rates R1, R2, R3 of the last rateControlled/equilibriumRates call, mol/(m³ s). */
  readonly R = new Float64Array(3);
  private lastT = NaN;
  private lastTr = NaN;

  constructor(rates: ZeldovichRateSet = ZELDOVICH_GRI30) {
    this.rates = rates;
  }

  /** Evaluate the forward constants kf at T (K); cached, so repeated calls at the same T are free. */
  setTemperature(T: number): void {
    if (T === this.lastT) return;
    this.lastT = T;
    this.lastTr = NaN;
    const r = this.rates;
    const lnT = Math.log(T);
    for (let i = 0; i < 3; i++) this.kf[i] = r.A[i] * Math.exp(r.b[i] * lnT - r.Ta[i] / T);
  }

  /** Evaluate kf and the reverse constants kr (detailed balance with our thermo) at T (K); cached. */
  setTemperatureWithReverse(T: number): void {
    this.setTemperature(T);
    if (T === this.lastTr) return;
    this.lastTr = T;
    const gN2 = g0RT(iN2, T);
    const gO2 = g0RT(iO2, T);
    const gO = g0RT(iO, T);
    const gN = g0RT(iN, T);
    const gNO = g0RT(iNO, T);
    const gOH = g0RT(iOH, T);
    const gH = g0RT(iH, T);
    // Kc = exp(−ΔG°/RT) (Δn = 0): R1 N+NO→N2+O, R2 N+O2→NO+O, R3 N+OH→NO+H
    this.kr[0] = this.kf[0] * Math.exp(gN2 + gO - gN - gNO);
    this.kr[1] = this.kf[1] * Math.exp(gNO + gO - gN - gO2);
    this.kr[2] = this.kf[2] * Math.exp(gNO + gH - gN - gOH);
  }

  /**
   * One-way equilibrium rates R1 = k1[N]e[NO]e, R2 = k2[N]e[O2]e, R3 = k3[N]e[OH]e, mol/(m³ s),
   * written to this.R (and returned) for equilibrium mole fractions Xeq (SPECIES order) at T, p.
   */
  equilibriumRates(T: number, p: number, Xeq: ArrayLike<number>): Float64Array {
    this.setTemperature(T);
    const c = p / (R_UNIVERSAL * T);
    const cN = Xeq[iN] * c;
    this.R[0] = this.kf[0] * cN * Xeq[iNO] * c;
    this.R[1] = this.kf[1] * cN * Xeq[iO2] * c;
    this.R[2] = this.kf[2] * cN * Xeq[iOH] * c;
    return this.R;
  }

  /**
   * Rate-controlled NO production (Heywood 1988 eq. 11.8), mol/(m³ s):
   * d[NO]/dt = 2R1(1 − α²)/(1 + αR1/(R2 + R3)), α = xNO / Xeq[NO].
   * Xeq: burned-gas EQUILIBRIUM mole fractions at (T, p) (NO and N included); xNO: actual NO
   * mole fraction. Returns 0 if the equilibrium NO or N mole fraction underflows (≤ 1e-30).
   * A negative xNO (explicit-integrator undershoot) is treated as 0: for α < −1/K the formula's
   * denominator changes sign and would drive NO further negative.
   */
  rateControlled(T: number, p: number, Xeq: ArrayLike<number>, xNO: number): number {
    const xNOe = Xeq[iNO];
    if (!(xNOe > 1e-30) || !(Xeq[iN] > 1e-30)) return 0;
    const R = this.equilibriumRates(T, p, Xeq);
    const alpha = Math.max(xNO, 0) / xNOe;
    const K = R[0] / (R[1] + R[2]);
    return (2 * R[0] * (1 - alpha * alpha)) / (1 + alpha * K);
  }

  /**
   * Characteristic relaxation time of NO toward equilibrium, s: τ_NO = [NO]e (1 + K)/(4 R1),
   * K = R1/(R2 + R3) (linearisation of eq. 11.8 about α = 1; the initial formation rate from
   * α = 0 is 2R1, i.e. α ≈ (1 + K) t / (2 τ_NO) for t ≪ τ_NO). Infinity if R1 underflows.
   */
  relaxationTime(T: number, p: number, Xeq: ArrayLike<number>): number {
    const R = this.equilibriumRates(T, p, Xeq);
    if (!(R[0] > 0)) return Infinity;
    const cNOe = (Xeq[iNO] * p) / (R_UNIVERSAL * T);
    return (cNOe * (1 + R[0] / (R[1] + R[2]))) / (4 * R[0]);
  }

  /**
   * Advance the rate-controlled NO mole fraction EXACTLY over dt (s) with T, p and Xeq frozen
   * (analytic integral of eq. 11.8; unconditionally stable, for operator-split sub-cycling):
   *   ∫ (1 + Kα)/(1 − α²) dα = (2R1/[NO]e) t,
   *   F(α) = −(1+K)/2 ln|1 − α| + (1−K)/2 ln(1 + α).
   * Works for formation (α < 1) and decomposition (α > 1); α never crosses 1. A negative xNO
   * is treated as 0.
   * Returns the new NO mole fraction (total moles assumed constant over the step).
   */
  advanceRateControlled(T: number, p: number, Xeq: ArrayLike<number>, xNO: number, dt: number): number {
    const xNOe = Xeq[iNO];
    if (!(xNOe > 1e-30) || !(Xeq[iN] > 1e-30) || !(dt > 0)) return xNO;
    const R = this.equilibriumRates(T, p, Xeq);
    const cNOe = (xNOe * p) / (R_UNIVERSAL * T);
    const K = R[0] / (R[1] + R[2]);
    const s = ((2 * R[0]) / cNOe) * dt;
    const a0 = Math.max(xNO, 0) / xNOe;
    if (a0 === 1 || !(s > 0)) return xNO;
    return xNOe * advanceAlpha(a0, s, K);
  }

  /**
   * Elementary forward/reverse rates for actual concentrations (no equilibrium or steady-state
   * assumption). X: mole fractions (SPECIES order; N2, O2, O, N, NO, OH, H used). Writes
   * out[0] = d[NO]/dt, out[1] = d[N]/dt, out[2] = d[O]/dt (mol/(m³ s)); returns out.
   */
  fullRates(T: number, p: number, X: ArrayLike<number>, out: Float64Array): Float64Array {
    this.setTemperatureWithReverse(T);
    const c = p / (R_UNIVERSAL * T);
    const N2 = X[iN2] * c;
    const O2 = X[iO2] * c;
    const O = X[iO] * c;
    const N = X[iN] * c;
    const NO = X[iNO] * c;
    const OH = X[iOH] * c;
    const H = X[iH] * c;
    const q1 = this.kf[0] * N * NO - this.kr[0] * N2 * O; // N + NO → N2 + O
    const q2 = this.kf[1] * N * O2 - this.kr[1] * NO * O; // N + O2 → NO + O
    const q3 = this.kf[2] * N * OH - this.kr[2] * NO * H; // N + OH → NO + H
    out[0] = -q1 + q2 + q3;
    out[1] = -q1 - q2 - q3;
    out[2] = q1 + q2;
    return out;
  }

  /**
   * Quasi-steady N mole fraction (d[N]/dt = 0 in fullRates) for the other species in X:
   * [N] = (k1r[N2][O] + k2r[NO][O] + k3r[NO][H]) / (k1f[NO] + k2f[O2] + k3f[OH]).
   */
  quasiSteadyN(T: number, p: number, X: ArrayLike<number>): number {
    this.setTemperatureWithReverse(T);
    const num = this.kr[0] * X[iN2] * X[iO] + this.kr[1] * X[iNO] * X[iO] + this.kr[2] * X[iNO] * X[iH];
    const den = this.kf[0] * X[iNO] + this.kf[1] * X[iO2] + this.kf[2] * X[iOH];
    // [N] = c·num/den with num, den in mole fractions ⇒ x_N = num/den (p cancels)
    return den > 0 ? num / den : 0;
  }
}

/**
 * α(s) solving dα/ds = (1 − α²)/(1 + Kα) from α0 over s ≥ 0 (see advanceRateControlled).
 * Safeguarded Newton in y = −ln|1 − α| (F is strictly monotone in y; slope ≥ min(1, (1+K)/2)
 * for α < 1 and > K for α > 1), converged to 1e-14.
 */
export function advanceAlpha(a0: number, s: number, K: number): number {
  const h1 = 0.5 * (1 + K);
  const h2 = 0.5 * (1 - K);
  if (a0 < 1) {
    // F(y) = h1 y + h2 ln(2 − e^−y), α = 1 − e^−y, increasing in y
    const y0 = -Math.log1p(-a0);
    const F1 = h1 * y0 + h2 * Math.log(2 - Math.exp(-y0)) + s;
    let y = y0 + s / Math.max(h1, 1e-300);
    let lo = y0;
    let hi = y0 + s / Math.min(1, h1);
    if (!(y <= hi)) y = hi;
    for (let it = 0; it < 60; it++) {
      const e = Math.exp(-y);
      const f = h1 * y + h2 * Math.log(2 - e) - F1;
      if (f > 0) hi = y;
      else lo = y;
      const df = h1 + (h2 * e) / (2 - e);
      let yn = y - f / df;
      if (!(yn > lo && yn < hi)) yn = 0.5 * (lo + hi);
      if (Math.abs(yn - y) <= 1e-14 * (1 + Math.abs(y))) {
        y = yn;
        break;
      }
      y = yn;
    }
    return -Math.expm1(-y);
  }
  // α > 1: G(y) = h1 y + h2 ln(2 + e^−y), α = 1 + e^−y, increasing in y (α decreasing to 1)
  const y0 = -Math.log(a0 - 1);
  const G1 = h1 * y0 + h2 * Math.log(2 + Math.exp(-y0)) + s;
  let lo = y0;
  let hi = y0 + s / Math.max(Math.min(h1, K), 1e-300);
  if (!Number.isFinite(hi)) hi = y0 + 1e3;
  let y = y0 + s / h1;
  if (!(y < hi)) y = 0.5 * (lo + hi);
  for (let it = 0; it < 100; it++) {
    const e = Math.exp(-y);
    const g = h1 * y + h2 * Math.log(2 + e) - G1;
    if (g > 0) hi = y;
    else lo = y;
    const dg = h1 - (h2 * e) / (2 + e);
    let yn = dg > 0 ? y - g / dg : 0.5 * (lo + hi);
    if (!(yn > lo && yn < hi)) yn = 0.5 * (lo + hi);
    if (Math.abs(yn - y) <= 1e-14 * (1 + Math.abs(y))) {
      y = yn;
      break;
    }
    y = yn;
  }
  return 1 + Math.exp(-y);
}

const DEFAULT_KINETICS = new ZeldovichKinetics(ZELDOVICH_GRI30);

/**
 * Module contract (DESIGN.md): rate-controlled extended-Zeldovich NO production rate,
 * mol/(m³ s), at T (K), p (Pa) for burned-gas equilibrium mole fractions Xeq (SPECIES order)
 * and actual NO mole fraction xNO. GRI-Mech 3.0 rate constants (ZELDOVICH_GRI30).
 */
export function zeldovichNORate(T: number, p: number, Xeq: ArrayLike<number>, xNO: number): number {
  return DEFAULT_KINETICS.rateControlled(T, p, Xeq, xNO);
}
