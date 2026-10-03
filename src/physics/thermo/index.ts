/**
 * Thermodynamic and transport properties (public API of src/physics/thermo).
 * See DESIGN.md §Module contracts → thermo. SI units throughout (mol, not kmol).
 */
export * from './thermo';
export * from './mixture';
export * from './fuels';
export * from './transport';
export { CANTERA_NASA_GAS_P_REF, P_REF_NASA, SPECIES_NASA7, SPECIES_TRANSPORT } from './species-data';
export type { Nasa7Record, TransportRecord } from './species-data';
