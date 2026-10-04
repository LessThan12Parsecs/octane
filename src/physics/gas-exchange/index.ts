/**
 * Gas exchange (public API of src/physics/gas-exchange): valve lift (polydyne, three-arc
 * flat-follower cam, measured table; createLiftProfile), valve and throttle flow areas / discharge
 * coefficients (incl. the side-valve roof and pocket-transfer stages; createValveFlowModel),
 * compressible restriction flow (single and two-in-series), the venturi + butterfly carburettor,
 * 0D plenums. See DESIGN.md §Module contracts → gas-exchange. SI units throughout (mol, not kmol).
 */
export * from './valve-lift';
export * from './cam-lift';
export * from './valve-flow';
export * from './orifice';
export * from './series-orifice';
export * from './throttle';
export * from './carburettor';
export * from './plenum';
export { Pchip } from './pchip';
