/**
 * Ideal-gas species thermodynamics from NASA-7 polynomials (NASA Glenn database as
 * bundled with Cantera, McBride, Gordon & Reno, NASA TM-4513, 1993). Molar SI units:
 * J/mol, J/(mol K). Absolute enthalpies (formation enthalpy included, 298.15 K reference);
 * standard state = ideal gas at P_REF_THERMO.
 *
 * Polynomial ranges: every species has [Tmin, Tmid] and (Tmid, Tmax] regions (Cantera
 * convention: T ≤ Tmid uses the low set). Outside [Tmin, Tmax] the heat capacity is held
 * constant at its boundary value and h, s are continued consistently:
 *   h(T) = h(Tb) + cp(Tb)(T − Tb),  s°(T) = s°(Tb) + cp(Tb) ln(T/Tb)
 * so dh/dT = cp and ds/dT = cp/T hold everywhere and nothing blows up.
 *
 * Tmid discontinuity: the two NASA-7 sets of a species are fitted jointly but are not
 * forced to match exactly at Tmid, so cp, h and s can jump there. For the species in this
 * set the measured jumps at 1000 K are |Δcp/cp| ≤ 1.4e-8, |Δh| ≤ 1e-3 J/mol,
 * |Δs| ≤ 3e-6 J/(mol K) (thermo.test.ts measures and bounds them). Cantera has exactly the
 * same jumps; the Newton inverses in mixture.ts are bracketed so they cannot stall on them.
 */
import { ATOMIC_WEIGHT, R_UNIVERSAL } from '../core/constants';
import { ELEMENT_COUNTS, ELEMENTS, NS, SPECIES } from '../core/species';
import { P_REF_NASA, SPECIES_NASA7 } from './species-data';

/**
 * Standard-state (reference) pressure of the polynomial data, Pa: 1 bar = 1e5 Pa.
 * NASA TM-4513 (McBride et al. 1993) defines the ideal-gas standard state at 1 bar, and the
 * fits' s°(298.15 K) reproduce the CODATA 1-bar key values (thermo.test.ts). NB: Cantera's
 * nasa_gas.yaml labels the same coefficients 1 atm (a ck2yaml default, CANTERA_NASA_GAS_P_REF
 * in species-data.ts); using 1 atm would bias every entropy by R ln(1.01325) =
 * 0.109 J/(mol K), dissociation equilibria by 0.4–0.9 % and a UV-equilibrium temperature at
 * ~3000 K by ~0.9 K. Cantera oracles must rebuild the species at 1 bar
 * (tools/reference/thermo_common.py: nasa_species / build_phase).
 */
export const P_REF_THERMO: number = P_REF_NASA;

/** Species molar masses, kg/mol, from core ATOMIC_WEIGHT × ELEMENT_COUNTS (SPECIES order). */
export const MOLAR_MASS: Float64Array = (() => {
  const w = new Float64Array(NS);
  for (let k = 0; k < NS; k++) {
    let m = 0;
    for (let e = 0; e < ELEMENTS.length; e++) m += ELEMENT_COUNTS[k][e] * ATOMIC_WEIGHT[ELEMENTS[e]];
    w[k] = m;
  }
  return w;
})();

/** Reciprocal molar masses, mol/kg. */
export const INV_MOLAR_MASS: Float64Array = MOLAR_MASS.map((m) => 1 / m);

// ---------------------------------------------------------------------------------------
// Flat coefficient table (Float64Array, stride STRIDE per species):
//   0 tmin, 1 tmid, 2 tmax,
//   LO = 3 .. 16  low-T block, HI = 17 .. 30 high-T block; each block (RB = 14 values):
//     a0 a1 a2 a3 a4 a5 a6 | a1/2 a2/3 a3/4 a4/5 | a2/2 a3/3 a4/4
//   31 cp/R(tmin), 32 h/R(tmin) [K], 33 s/R(tmin), 34 cp/R(tmax), 35 h/R(tmax) [K],
//   36 s/R(tmax), 37 ln(tmin), 38 ln(tmax)
// Pre-scaled coefficients let every property be one Horner polynomial (no divisions).
// ---------------------------------------------------------------------------------------
const STRIDE = 39;
const LO = 3;
const HI = 17;
const C = new Float64Array(NS * STRIDE);

/** cp/R from the block at offset o. */
function polyCpR(o: number, T: number): number {
  return C[o] + T * (C[o + 1] + T * (C[o + 2] + T * (C[o + 3] + T * C[o + 4])));
}
/** h/(RT) from the block at offset o. */
function polyHRT(o: number, T: number): number {
  return C[o] + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] / T;
}
/** s°/R from the block at offset o. */
function polySR(o: number, T: number, lnT: number): number {
  return C[o] * lnT + T * (C[o + 1] + T * (C[o + 11] + T * (C[o + 12] + T * C[o + 13]))) + C[o + 6];
}

for (let k = 0; k < NS; k++) {
  const r = SPECIES_NASA7[k];
  if (r.name !== SPECIES[k]) {
    throw new Error(`species-data.ts order mismatch at ${k}: ${r.name} vs ${SPECIES[k]} (regenerate)`);
  }
  const b = k * STRIDE;
  C[b] = r.tmin;
  C[b + 1] = r.tmid;
  C[b + 2] = r.tmax;
  for (const [off, a] of [[LO, r.low], [HI, r.high]] as const) {
    const o = b + off;
    for (let i = 0; i < 7; i++) C[o + i] = a[i];
    C[o + 7] = a[1] / 2;
    C[o + 8] = a[2] / 3;
    C[o + 9] = a[3] / 4;
    C[o + 10] = a[4] / 5;
    C[o + 11] = a[2] / 2;
    C[o + 12] = a[3] / 3;
    C[o + 13] = a[4] / 4;
  }
  const oLo = b + LO;
  const oHi = r.tmax <= r.tmid ? b + LO : b + HI; // Tmax evaluated with the set that covers it
  C[b + 31] = polyCpR(oLo, r.tmin);
  C[b + 32] = polyHRT(oLo, r.tmin) * r.tmin;
  C[b + 33] = polySR(oLo, r.tmin, Math.log(r.tmin));
  C[b + 34] = polyCpR(oHi, r.tmax);
  C[b + 35] = polyHRT(oHi, r.tmax) * r.tmax;
  C[b + 36] = polySR(oHi, r.tmax, Math.log(r.tmax));
  C[b + 37] = Math.log(r.tmin);
  C[b + 38] = Math.log(r.tmax);
}

/** Lower validity limit of species k's polynomial fit, K (constant-cp extrapolation below). */
export const speciesTmin = (k: number): number => C[k * STRIDE];
/** Upper validity limit of species k's polynomial fit, K (constant-cp extrapolation above). */
export const speciesTmax = (k: number): number => C[k * STRIDE + 2];
/** Common-range break temperature of species k's fit, K. */
export const speciesTmid = (k: number): number => C[k * STRIDE + 1];

// ---- single-species nondimensional kernels (with constant-cp extrapolation) -----------

/** cp°/R of species k at T (K). Dimensionless. */
export function speciesCpR(k: number, T: number): number {
  const b = k * STRIDE;
  if (T < C[b]) return C[b + 31];
  if (T > C[b + 2]) return C[b + 34];
  return polyCpR(T <= C[b + 1] ? b + LO : b + HI, T);
}

/** h/(R T) of species k at T (K). Dimensionless. */
export function speciesHRT(k: number, T: number): number {
  const b = k * STRIDE;
  if (T < C[b]) return (C[b + 32] + C[b + 31] * (T - C[b])) / T;
  if (T > C[b + 2]) return (C[b + 35] + C[b + 34] * (T - C[b + 2])) / T;
  return polyHRT(T <= C[b + 1] ? b + LO : b + HI, T);
}

/** s°/R of species k at T (K) and P_REF_THERMO. Dimensionless. */
export function speciesS0R(k: number, T: number): number {
  const b = k * STRIDE;
  const lnT = Math.log(T);
  if (T < C[b]) return C[b + 33] + C[b + 31] * (lnT - C[b + 37]);
  if (T > C[b + 2]) return C[b + 36] + C[b + 34] * (lnT - C[b + 38]);
  return polySR(T <= C[b + 1] ? b + LO : b + HI, T, lnT);
}

/** Molar heat capacity cp° of species k at T (K), J/(mol K). */
export const speciesCp = (k: number, T: number): number => R_UNIVERSAL * speciesCpR(k, T);
/** Molar enthalpy h of species k at T (K), J/mol (absolute, formation included). */
export const speciesH = (k: number, T: number): number => R_UNIVERSAL * T * speciesHRT(k, T);
/** Standard-state molar entropy s° of species k at T (K) and P_REF_THERMO, J/(mol K). */
export const speciesS0 = (k: number, T: number): number => R_UNIVERSAL * speciesS0R(k, T);
/** Standard-state molar Gibbs energy g° = h − T s° of species k at T (K), J/mol. */
export const speciesG0 = (k: number, T: number): number =>
  R_UNIVERSAL * T * (speciesHRT(k, T) - speciesS0R(k, T));

// ---- all-species hot paths --------------------------------------------------------------

/**
 * Evaluate cp°/R, h/(RT) and s°/R of ALL species at T (K) in one pass (hot path; allocates
 * nothing). Outputs are Float64Array(NS) in SPECIES order, dimensionless; s° at P_REF_THERMO.
 */
export function evalAllNondim(T: number, cpR: Float64Array, hRT: Float64Array, s0R: Float64Array): void {
  const lnT = Math.log(T);
  const invT = 1 / T;
  for (let k = 0, b = 0; k < NS; k++, b += STRIDE) {
    if (T < C[b]) {
      const c = C[b + 31];
      cpR[k] = c;
      hRT[k] = (C[b + 32] + c * (T - C[b])) * invT;
      s0R[k] = C[b + 33] + c * (lnT - C[b + 37]);
    } else if (T > C[b + 2]) {
      const c = C[b + 34];
      cpR[k] = c;
      hRT[k] = (C[b + 35] + c * (T - C[b + 2])) * invT;
      s0R[k] = C[b + 36] + c * (lnT - C[b + 38]);
    } else {
      const o = T <= C[b + 1] ? b + LO : b + HI;
      const a0 = C[o];
      cpR[k] = a0 + T * (C[o + 1] + T * (C[o + 2] + T * (C[o + 3] + T * C[o + 4])));
      hRT[k] = a0 + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] * invT;
      s0R[k] = a0 * lnT + T * (C[o + 1] + T * (C[o + 11] + T * (C[o + 12] + T * C[o + 13]))) + C[o + 6];
    }
  }
}

/**
 * cp°/R and h/(RT) of all species at T (K) — as evalAllNondim without the entropy
 * (no logarithm). Outputs Float64Array(NS), dimensionless.
 */
export function evalAllCpH(T: number, cpR: Float64Array, hRT: Float64Array): void {
  const invT = 1 / T;
  for (let k = 0, b = 0; k < NS; k++, b += STRIDE) {
    if (T < C[b]) {
      const c = C[b + 31];
      cpR[k] = c;
      hRT[k] = (C[b + 32] + c * (T - C[b])) * invT;
    } else if (T > C[b + 2]) {
      const c = C[b + 34];
      cpR[k] = c;
      hRT[k] = (C[b + 35] + c * (T - C[b + 2])) * invT;
    } else {
      const o = T <= C[b + 1] ? b + LO : b + HI;
      const a0 = C[o];
      cpR[k] = a0 + T * (C[o + 1] + T * (C[o + 2] + T * (C[o + 3] + T * C[o + 4])));
      hRT[k] = a0 + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] * invT;
    }
  }
}

/**
 * Standard-state g°/(RT) = h/(RT) − s°/R of all species at T (K) (for equilibrium /
 * equilibrium constants). Output Float64Array(NS), dimensionless, at P_REF_THERMO.
 */
export function evalAllG0RT(T: number, g0RT: Float64Array): void {
  const lnT = Math.log(T);
  const invT = 1 / T;
  for (let k = 0, b = 0; k < NS; k++, b += STRIDE) {
    if (T < C[b]) {
      const c = C[b + 31];
      g0RT[k] = (C[b + 32] + c * (T - C[b])) * invT - (C[b + 33] + c * (lnT - C[b + 37]));
    } else if (T > C[b + 2]) {
      const c = C[b + 34];
      g0RT[k] = (C[b + 35] + c * (T - C[b + 2])) * invT - (C[b + 36] + c * (lnT - C[b + 38]));
    } else {
      const o = T <= C[b + 1] ? b + LO : b + HI;
      const a0 = C[o];
      const h = a0 + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] * invT;
      const s = a0 * lnT + T * (C[o + 1] + T * (C[o + 11] + T * (C[o + 12] + T * C[o + 13]))) + C[o + 6];
      g0RT[k] = h - s;
    }
  }
}

// ---- fused mixture sums (used by mixture.ts) --------------------------------------------

/**
 * Results of the last mixSums call:
 *   [0] Σ x_k cp_k/R, [1] Σ x_k h_k/(RT), [2] Σ x_k s°_k/R − Σ_{x_k>0} x_k ln x_k (only if
 *   withS; the mixing term is skipped for x_k ≤ 0, so d[2]/dT = [0]/T holds for any x),
 *   [3] Σ x_k, [4] Σ x_k M_k (kg/mol).
 * Module-level scratch (not re-entrant). Prefer the mixture.ts functions.
 */
export const MIX_SUMS = new Float64Array(5);

/**
 * One pass over the species with x_k ≠ 0 accumulating the MIX_SUMS above at T (K). The
 * entropy sum (with the ideal-mixing −ln x_k term) is only meaningful for normalised mole
 * fractions and is skipped unless withS (it costs a logarithm per species).
 */
export function mixSums(X: Float64Array, T: number, withS: boolean): void {
  if (withS) {
    mixSumsS(X, T);
    return;
  }
  const invT = 1 / T;
  let sc = 0;
  let sh = 0;
  let sx = 0;
  let sm = 0;
  for (let k = 0, b = 0; k < NS; k++, b += STRIDE) {
    const x = X[k];
    if (x === 0) continue;
    let cp: number;
    let h: number;
    if (T < C[b]) {
      cp = C[b + 31];
      h = (C[b + 32] + cp * (T - C[b])) * invT;
    } else if (T > C[b + 2]) {
      cp = C[b + 34];
      h = (C[b + 35] + cp * (T - C[b + 2])) * invT;
    } else {
      const o = T <= C[b + 1] ? b + LO : b + HI;
      const a0 = C[o];
      cp = a0 + T * (C[o + 1] + T * (C[o + 2] + T * (C[o + 3] + T * C[o + 4])));
      h = a0 + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] * invT;
    }
    sc += x * cp;
    sh += x * h;
    sx += x;
    sm += x * MOLAR_MASS[k];
  }
  MIX_SUMS[0] = sc;
  MIX_SUMS[1] = sh;
  MIX_SUMS[2] = 0;
  MIX_SUMS[3] = sx;
  MIX_SUMS[4] = sm;
}

function mixSumsS(X: Float64Array, T: number): void {
  const invT = 1 / T;
  const lnT = Math.log(T);
  let sc = 0;
  let sh = 0;
  let ss = 0;
  let sx = 0;
  let sm = 0;
  for (let k = 0, b = 0; k < NS; k++, b += STRIDE) {
    const x = X[k];
    if (x === 0) continue;
    let cp: number;
    let h: number;
    let s: number;
    if (T < C[b]) {
      cp = C[b + 31];
      h = (C[b + 32] + cp * (T - C[b])) * invT;
      s = C[b + 33] + cp * (lnT - C[b + 37]);
    } else if (T > C[b + 2]) {
      cp = C[b + 34];
      h = (C[b + 35] + cp * (T - C[b + 2])) * invT;
      s = C[b + 36] + cp * (lnT - C[b + 38]);
    } else {
      const o = T <= C[b + 1] ? b + LO : b + HI;
      const a0 = C[o];
      cp = a0 + T * (C[o + 1] + T * (C[o + 2] + T * (C[o + 3] + T * C[o + 4])));
      h = a0 + T * (C[o + 7] + T * (C[o + 8] + T * (C[o + 9] + T * C[o + 10]))) + C[o + 5] * invT;
      s = a0 * lnT + T * (C[o + 1] + T * (C[o + 11] + T * (C[o + 12] + T * C[o + 13]))) + C[o + 6];
    }
    sc += x * cp;
    sh += x * h;
    ss += x * s;
    if (x > 0) ss -= x * Math.log(x);
    sx += x;
    sm += x * MOLAR_MASS[k];
  }
  MIX_SUMS[0] = sc;
  MIX_SUMS[1] = sh;
  MIX_SUMS[2] = ss;
  MIX_SUMS[3] = sx;
  MIX_SUMS[4] = sm;
}

/** Index lookup helper for tests/UI: species name → index (−1 if absent). */
export const speciesIndex = (name: string): number => (SPECIES as readonly string[]).indexOf(name);
