/**
 * Equilibrium thermodynamic properties of the burned gas at (T, p): mass-specific h, u, s,
 * the EQUILIBRIUM (composition-shifting) cp and cv, and the derivatives (∂lnV/∂lnT)_p and
 * (∂lnV/∂lnp)_T needed by the two-zone cylinder model — all analytic (NASA RP-1311 §2.5).
 */
import type { EqProperties, EquilibriumSolver } from './solver';
import { newEqProperties } from './solver';

/**
 * Solve the TP equilibrium for element amounts b (mol, ELEMENTS order) at T (K), p (Pa) with
 * `solver` (warm-started) and return its equilibrium properties (see EqProperties; SI,
 * mass-specific). Pass `out` to stay allocation-free. out.valid is false if the solve failed.
 */
export function equilibriumProperties(
  solver: EquilibriumSolver,
  b: Float64Array,
  T: number,
  p: number,
  out: EqProperties = newEqProperties(),
): EqProperties {
  const res = solver.solveTP(b, T, p);
  solver.properties(out);
  if (!res.converged) out.valid = false;
  return out;
}
