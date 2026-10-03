/**
 * Diagnostics for equilibrium results (tests, debugging, assertions in development builds).
 * Not hot-path code: they allocate small scratch arrays.
 */
import { ELEMENT_COUNTS, N_EQ, NE, NS } from '../core/species';
import { elementMoles } from '../thermo/mixture';
import { evalAllNondim, P_REF_THERMO } from '../thermo/thermo';
import { EQ_B_REL_MIN, type EqResult } from './solver';

/**
 * Largest relative element-balance error max_e |Σ_k a_ke N_k − b_e| / b_e of a result versus the
 * requested element amounts b (mol, ELEMENTS order), over the elements the solver treats as
 * present (b_e/Σb > EQ_B_REL_MIN). Dimensionless.
 */
export function elementBalanceError(res: EqResult, b: Float64Array): number {
  const bb = elementMoles(res.N, new Float64Array(NE));
  let tot = 0;
  for (let e = 0; e < NE; e++) tot += Math.max(b[e], 0);
  let worst = 0;
  for (let e = 0; e < NE; e++) {
    if (!(b[e] > EQ_B_REL_MIN * tot)) continue;
    worst = Math.max(worst, Math.abs(bb[e] - b[e]) / b[e]);
  }
  return worst;
}

/**
 * Largest violation of the equilibrium condition μ_j/RT = Σ_e a_je π_e (RP-1311 eq. 2.13 with
 * CEA's modified multipliers) over the species present in `res` (X_j > 0, which includes trace
 * species down to ~1e-300), using the element potentials res.pi. Dimensionless (≈ relative
 * error of each X_j).
 */
export function stationarityError(res: EqResult): number {
  const cpR = new Float64Array(NS);
  const hRT = new Float64Array(NS);
  const s0R = new Float64Array(NS);
  evalAllNondim(res.T, cpR, hRT, s0R);
  const lnP = Math.log(res.p / P_REF_THERMO);
  let worst = 0;
  for (let k = 0; k < N_EQ; k++) {
    if (!(res.X[k] > 0)) continue;
    const mu = hRT[k] - s0R[k] + Math.log(res.X[k]) + lnP;
    let s = 0;
    for (let e = 0; e < NE; e++) s += ELEMENT_COUNTS[k][e] * res.pi[e];
    worst = Math.max(worst, Math.abs(mu - s));
  }
  return worst;
}
