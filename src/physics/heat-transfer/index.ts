/**
 * Wall heat transfer (public API of src/physics/heat-transfer): Woschni, Hohenberg and Annand
 * correlations and the head/valve/piston/liner (and burned/unburned) split of the heat loss.
 * See DESIGN.md §Module contracts → heat-transfer. SI units throughout.
 */
export * from './correlations';
export * from './wall-heat';
