/**
 * Chemical equilibrium of the burned gas (public API of src/physics/equilibrium).
 * See DESIGN.md §Module contracts → equilibrium. SI units throughout (mol, not kmol).
 */
export { EquilibriumSolver, newEqProperties, EQ_T_MIN, EQ_T_MAX, EQ_B_REL_MIN, EQ_B_NEG_TOL } from './solver';
export type { EqResult, EqProperties } from './solver';
export { equilibriumProperties } from './properties';
export { adiabaticFlameTemperature } from './flame';
export type { FlameResult } from './flame';
export { elementBalanceError, stationarityError } from './diagnostics';
