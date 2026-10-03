/**
 * Gas exchange (public API of src/physics/gas-exchange): valve lift, valve and throttle flow
 * areas / discharge coefficients, compressible restriction flow, 0D plenums.
 * See DESIGN.md §Module contracts → gas-exchange. SI units throughout (mol, not kmol).
 */
export * from './valve-lift';
export * from './valve-flow';
export * from './orifice';
export * from './throttle';
export * from './plenum';
export { Pchip } from './pchip';
