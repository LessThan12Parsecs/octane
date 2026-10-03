/**
 * Adiabatic (equilibrium) flame temperature of a reactant mixture burned at constant pressure
 * (HP) or constant volume (UV).
 */
import { R_UNIVERSAL } from '../core/constants';
import { NE, NS } from '../core/species';
import { elementMoles, mixHMolar, mixUMolar, totalMoles } from '../thermo/mixture';
import { EquilibriumSolver } from './solver';

/** Result of adiabaticFlameTemperature. */
export interface FlameResult {
  /** Adiabatic flame temperature, K. */
  T: number;
  /** Final pressure, Pa (= p0 for HP; n_b R T / V0 for UV). */
  p: number;
  /** Equilibrium product mole fractions, Float64Array(NS). */
  X: Float64Array;
  /** Product moles per mole of reactants, mol/mol. */
  molesPerMoleReactant: number;
  converged: boolean;
}

let defaultSolver: EquilibriumSolver | null = null;

/**
 * Adiabatic flame temperature of reactants with mole fractions Xreactants (Float64Array(NS),
 * fuel species allowed; need not be normalised) initially at T0 (K), p0 (Pa):
 *  - 'HP': constant pressure, H_products(T, p0) = H_reactants(T0);
 *  - 'UV': constant volume, U_products(T, V0) = U_reactants(T0), V0 = n R T0/p0.
 * Products are in chemical equilibrium over the N_EQ burned-gas species. Not a hot path
 * (allocates the result); pass a solver to reuse its warm start.
 */
export function adiabaticFlameTemperature(
  Xreactants: Float64Array,
  T0: number,
  p0: number,
  mode: 'HP' | 'UV',
  solver: EquilibriumSolver = (defaultSolver ??= new EquilibriumSolver()),
): FlameResult {
  const nR = totalMoles(Xreactants);
  const b = elementMoles(Xreactants, new Float64Array(NE));
  let res;
  if (mode === 'HP') {
    const H0 = mixHMolar(Xreactants, T0); // J for nR mol (linear in X)
    res = solver.solveHP(b, H0, p0, 2000);
  } else {
    const U0 = mixUMolar(Xreactants, T0);
    const V0 = (nR * R_UNIVERSAL * T0) / p0;
    res = solver.solveUV(b, U0, V0, 2500);
  }
  const X = new Float64Array(NS);
  X.set(res.X);
  return { T: res.T, p: res.p, X, molesPerMoleReactant: res.nTotal / nR, converged: res.converged };
}
