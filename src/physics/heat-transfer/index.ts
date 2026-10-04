/**
 * Wall heat transfer (public API of src/physics/heat-transfer): Woschni, Hohenberg and Annand
 * correlations and the head/valve/piston/liner (and burned/unburned) split of the heat loss — the
 * 5-surface flat-disc functions and the chamber-generic 6-surface path (WALL_* indices incl. the
 * L-head block deck, per-surface burned fractions from combustion/chamber.ts).
 * See DESIGN.md §Module contracts → heat-transfer. SI units throughout.
 */
export * from './correlations';
export * from './wall-heat';
