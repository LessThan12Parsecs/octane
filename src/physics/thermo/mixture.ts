/**
 * Ideal-gas mixture properties and composition conversions.
 *
 * Composition vectors are Float64Array(NS) in SPECIES order: N = moles (mol), X = mole
 * fractions, Y = mass fractions. Element vectors are Float64Array(NE) in ELEMENTS order.
 *
 * Linearity: the molar functions mixCpMolar/mixCvMolar/mixHMolar/mixUMolar and the
 * temperatureFrom{H,U}Molar inverses are LINEAR in the composition vector (they never
 * normalise it), so passing moles N instead of X gives extensive values (J/K, J) and
 * lets an inverse solve U(T) = U_total directly. Entropy is not linear (ideal-mixing
 * term), so mixSMolar/mixSMass require true mole fractions.
 *
 * Species with X_k = 0 are skipped (and contribute exactly 0 to the mixing entropy).
 * Hot-path functions allocate nothing.
 */
import { R_UNIVERSAL } from '../core/constants';
import { ELEMENT_COUNTS, NE, NS } from '../core/species';
import { MIX_SUMS, MOLAR_MASS, mixSums, P_REF_THERMO } from './thermo';

// ---------------------------------------------------------------------------------------
// Composition conversions
// ---------------------------------------------------------------------------------------

/** Total moles Σ N_k, mol. */
export function totalMoles(N: Float64Array): number {
  let s = 0;
  for (let k = 0; k < NS; k++) s += N[k];
  return s;
}

/** Total mass Σ N_k M_k, kg (N in mol). */
export function totalMass(N: Float64Array): number {
  let s = 0;
  for (let k = 0; k < NS; k++) s += N[k] * MOLAR_MASS[k];
  return s;
}

/** Mole fractions from moles (mol). Writes into X if given (may alias N). Returns X. */
export function molesToX(N: Float64Array, X: Float64Array = new Float64Array(NS)): Float64Array {
  const inv = 1 / totalMoles(N);
  for (let k = 0; k < NS; k++) X[k] = N[k] * inv;
  return X;
}

/** Mass fractions from mole fractions. Writes into Y if given (may alias X). Returns Y. */
export function xToY(X: Float64Array, Y: Float64Array = new Float64Array(NS)): Float64Array {
  let m = 0;
  for (let k = 0; k < NS; k++) m += X[k] * MOLAR_MASS[k];
  const inv = 1 / m;
  for (let k = 0; k < NS; k++) Y[k] = X[k] * MOLAR_MASS[k] * inv;
  return Y;
}

/** Mole fractions from mass fractions. Writes into X if given (may alias Y). Returns X. */
export function yToX(Y: Float64Array, X: Float64Array = new Float64Array(NS)): Float64Array {
  let s = 0;
  for (let k = 0; k < NS; k++) s += Y[k] / MOLAR_MASS[k];
  const inv = 1 / s;
  for (let k = 0; k < NS; k++) X[k] = Y[k] / MOLAR_MASS[k] * inv;
  return X;
}

/**
 * Element amounts b_e = Σ_k a_ke N_k (mol of atoms, ELEMENTS order). Also valid with X
 * (atoms per mole of mixture). Writes into b if given. Returns b.
 */
export function elementMoles(N: Float64Array, b: Float64Array = new Float64Array(NE)): Float64Array {
  b.fill(0);
  for (let k = 0; k < NS; k++) {
    const n = N[k];
    if (n === 0) continue;
    const a = ELEMENT_COUNTS[k];
    for (let e = 0; e < NE; e++) b[e] += a[e] * n;
  }
  return b;
}

// ---------------------------------------------------------------------------------------
// Mixture molar mass and gas constant
// ---------------------------------------------------------------------------------------

/** Mean molar mass Σ X_k M_k, kg/mol (X normalised). */
export function mixMolarMass(X: Float64Array): number {
  let m = 0;
  for (let k = 0; k < NS; k++) m += X[k] * MOLAR_MASS[k];
  return m;
}

/** Specific gas constant R_u / M_mix, J/(kg K). */
export const mixGasConstant = (X: Float64Array): number => R_UNIVERSAL / mixMolarMass(X);

// ---------------------------------------------------------------------------------------
// Molar properties (linear in X — see header)
// ---------------------------------------------------------------------------------------

// All property functions below make ONE pass over the non-zero species (thermo.mixSums).
const S = MIX_SUMS;

/** Mixture molar cp at T (K), J/(mol K). Linear in X (N → J/K). */
export function mixCpMolar(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return R_UNIVERSAL * S[0];
}

/** Mixture molar cv = cp − R at T (K), J/(mol K). Linear in X (N → J/K). */
export function mixCvMolar(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return R_UNIVERSAL * (S[0] - S[3]);
}

/** Mixture molar enthalpy at T (K), J/mol (absolute). Linear in X (N → J). */
export function mixHMolar(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return R_UNIVERSAL * T * S[1];
}

/** Mixture molar internal energy u = h − R T at T (K), J/mol. Linear in X (N → J). */
export function mixUMolar(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return R_UNIVERSAL * T * (S[1] - S[3]);
}

/**
 * Mixture molar entropy of an ideal-gas mixture at T (K), p (Pa), J/(mol K):
 *   s = Σ X_k [s°_k(T) − R ln(X_k p / p_ref)],  species with X_k = 0 contribute 0.
 * X must be normalised mole fractions.
 */
export function mixSMolar(X: Float64Array, T: number, p: number): number {
  mixSums(X, T, true);
  return R_UNIVERSAL * (S[2] - S[3] * Math.log(p / P_REF_THERMO));
}

// ---------------------------------------------------------------------------------------
// Mass-specific properties (X normalised)
// ---------------------------------------------------------------------------------------

/** Mixture cp at T (K), J/(kg K). */
export function mixCpMass(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return (R_UNIVERSAL * S[0]) / S[4];
}
/** Mixture cv at T (K), J/(kg K). */
export function mixCvMass(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return (R_UNIVERSAL * (S[0] - S[3])) / S[4];
}
/** Mixture enthalpy at T (K), J/kg (absolute). */
export function mixHMass(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return (R_UNIVERSAL * T * S[1]) / S[4];
}
/** Mixture internal energy at T (K), J/kg (absolute). */
export function mixUMass(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return (R_UNIVERSAL * T * (S[1] - S[3])) / S[4];
}
/** Mixture entropy at T (K), p (Pa), J/(kg K) (ideal mixing; X normalised). */
export function mixSMass(X: Float64Array, T: number, p: number): number {
  mixSums(X, T, true);
  return (R_UNIVERSAL * (S[2] - S[3] * Math.log(p / P_REF_THERMO))) / S[4];
}

/** Ratio of specific heats γ = cp/cv at T (K) (frozen composition). */
export function mixGamma(X: Float64Array, T: number): number {
  mixSums(X, T, false);
  return S[0] / (S[0] - S[3]);
}

/** Mass-specific mixture state (output of mixStateMass). */
export interface MixState {
  /** J/(kg K) */
  cp: number;
  /** J/(kg K) */
  cv: number;
  /** J/kg (absolute) */
  h: number;
  /** J/kg (absolute) */
  u: number;
  /** cp/cv */
  gamma: number;
  /** Mean molar mass, kg/mol. */
  molarMass: number;
}

/**
 * Mass-specific cp, cv, h, u and γ at T (K) in one pass over the species (X normalised).
 * Writes into out (allocation-free when supplied).
 */
export function mixStateMass(
  X: Float64Array,
  T: number,
  out: MixState = { cp: 0, cv: 0, h: 0, u: 0, gamma: 0, molarMass: 0 },
): MixState {
  mixSums(X, T, false);
  const inv = 1 / S[4];
  out.molarMass = S[4] / S[3];
  out.cp = R_UNIVERSAL * S[0] * inv;
  out.cv = R_UNIVERSAL * (S[0] - S[3]) * inv;
  out.h = R_UNIVERSAL * T * S[1] * inv;
  out.u = R_UNIVERSAL * T * (S[1] - S[3]) * inv;
  out.gamma = S[0] / (S[0] - S[3]);
  return out;
}

/** Frozen speed of sound sqrt(γ R T / M), m/s. */
export function mixSoundSpeed(X: Float64Array, T: number): number {
  return Math.sqrt((mixGamma(X, T) * R_UNIVERSAL * T) / mixMolarMass(X));
}

/** Ideal-gas mass density p M / (R T), kg/m³. */
export function mixDensity(X: Float64Array, T: number, p: number): number {
  return (p * mixMolarMass(X)) / (R_UNIVERSAL * T);
}

// ---------------------------------------------------------------------------------------
// Inverses: safeguarded Newton (bisection fallback) — allocation-free
// ---------------------------------------------------------------------------------------

/** Lowest temperature the inverses will return, K. */
export const T_INVERSE_MIN = 1;
/** Highest temperature the inverses will return, K. */
export const T_INVERSE_MAX = 1e5;
/** Relative convergence tolerance of the inverses (|ΔT|/T). */
export const T_INVERSE_RTOL = 1e-10;
const MAX_ITER = 100;

/** Diagnostics of the most recent inverse call (for tests / debugging; no allocation). */
export const lastInverse = { iterations: 0, converged: false };

// kind: 0 = h (molar, linear), 1 = u (molar, linear), 2 = s (molar, X normalised).
// Evaluates f(T) and f'(T) in one mixSums pass; f' is left in DFDT[0] (a typed-array slot: a
// module-level `let` double is boxed as a heap number on every write — hot paths must not allocate).
const DFDT = new Float64Array(1);
function residual(kind: number, X: Float64Array, T: number, lnPr: number): number {
  if (kind === 2) {
    mixSums(X, T, true);
    DFDT[0] = (R_UNIVERSAL * S[0]) / T;
    return R_UNIVERSAL * (S[2] - S[3] * lnPr);
  }
  mixSums(X, T, false);
  if (kind === 0) {
    DFDT[0] = R_UNIVERSAL * S[0];
    return R_UNIVERSAL * T * S[1];
  }
  DFDT[0] = R_UNIVERSAL * (S[0] - S[3]);
  return R_UNIVERSAL * T * (S[1] - S[3]);
}

/**
 * Solve f(T) = target for a strictly increasing f (h, u or s) with Newton steps kept inside
 * a bracket that is tightened every iteration; a step leaving the bracket is replaced by
 * bisection. Converges when the Newton correction |ΔT| ≤ rtol·T (the returned T is then
 * accurate to ~machine precision, quadratic convergence) or, failing that, when the bracket
 * is narrower than rtol·T. Note: at a polynomial Tmid the (tiny) jump in h/s means the
 * returned T can sit on the other side of Tmid by |Δh|/cp ≲ 1e-5 K.
 * If the target lies outside [f(T_INVERSE_MIN), f(T_INVERSE_MAX)] the bound is returned
 * and lastInverse.converged = false.
 *
 * The convergence test on the Newton correction comes BEFORE the bracket test: once the
 * correction drops below one ulp of T, T + ΔT == T, which equals the bracket end that was
 * just set, so a bracket-first test would reject the converged step and fall into ~35
 * bisection steps (seen in ~1 % of calls, even from warm starts, before this ordering).
 */
function invert(kind: number, X: Float64Array, target: number, p: number, Tguess: number): number {
  const lnPr = kind === 2 ? Math.log(p / P_REF_THERMO) : 0;
  let T = Tguess > T_INVERSE_MIN && Tguess < T_INVERSE_MAX ? Tguess : 1000;
  let lo = T_INVERSE_MIN;
  let hi = T_INVERSE_MAX;
  let fLoKnown = false;
  let fHiKnown = false;
  lastInverse.converged = false;
  for (let it = 1; it <= MAX_ITER; it++) {
    const f = residual(kind, X, T, lnPr) - target;
    if (f === 0) {
      lastInverse.iterations = it;
      lastInverse.converged = true;
      return T;
    }
    if (f < 0) {
      lo = T;
      fLoKnown = true;
    } else {
      hi = T;
      fHiKnown = true;
    }
    const dT = -f / DFDT[0];
    if (Math.abs(dT) <= T_INVERSE_RTOL * T) {
      // Newton correction below tolerance: |f| ≤ rtol·T·f′, i.e. at the root to within rtol.
      lastInverse.iterations = it;
      lastInverse.converged = true;
      return T + dT;
    }
    const Tn = T + dT;
    if (Tn > lo && Tn < hi) {
      T = Tn;
    } else {
      // Newton left the bracket: bisect (geometric mean while it still spans decades).
      T = hi / lo > 4 ? Math.sqrt(lo * hi) : 0.5 * (lo + hi);
    }
    if (hi - lo <= T_INVERSE_RTOL * lo) {
      lastInverse.iterations = it;
      lastInverse.converged = fLoKnown && fHiKnown;
      return 0.5 * (lo + hi);
    }
  }
  lastInverse.iterations = MAX_ITER;
  return T;
}

/**
 * Temperature (K) at which the mixture molar enthalpy equals h (J/mol). Linear in X, so
 * (N, H_total [J]) works too. Tguess (K) is a warm start.
 */
export const temperatureFromHMolar = (X: Float64Array, h: number, Tguess = 1000): number =>
  invert(0, X, h, 0, Tguess);

/**
 * Temperature (K) at which the mixture molar internal energy equals u (J/mol). Linear in X,
 * so (N, U_total [J]) works too. Tguess (K) is a warm start.
 */
export const temperatureFromUMolar = (X: Float64Array, u: number, Tguess = 1000): number =>
  invert(1, X, u, 0, Tguess);

/** Temperature (K) at which the mixture mass-specific enthalpy equals h (J/kg). */
export const temperatureFromH = (X: Float64Array, h: number, Tguess = 1000): number =>
  invert(0, X, h * mixMolarMass(X), 0, Tguess);

/** Temperature (K) at which the mixture mass-specific internal energy equals u (J/kg). */
export const temperatureFromU = (X: Float64Array, u: number, Tguess = 1000): number =>
  invert(1, X, u * mixMolarMass(X), 0, Tguess);

/**
 * Temperature (K) at which the mixture mass-specific entropy at pressure p (Pa) equals
 * s (J/(kg K)). X must be normalised.
 */
export const temperatureFromS = (X: Float64Array, s: number, p: number, Tguess = 1000): number =>
  invert(2, X, s * mixMolarMass(X), p, Tguess);
