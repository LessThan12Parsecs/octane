/**
 * Contract entry point (DESIGN.md §Module contracts → heat-transfer): the correlations live in
 * correlations.ts, the surface / zone split in wall-heat.ts; this file re-exports them.
 */
export {
  ANNAND_DEFAULTS,
  ANNAND_RADIATION,
  annandCoefficient,
  annandFlux,
  hohenbergCoefficient,
  WOSCHNI_CONSTANTS,
  woschniCoefficient,
  woschniMotoredPressure,
  woschniVelocity,
  type WoschniInputs,
  type WoschniPhase,
  type WoschniVariant,
} from './correlations';
export {
  flatChamberAreas,
  newChamberAreas,
  newWallHeatResult,
  wallHeatLoss,
  wallHeatLossTwoZone,
  type BurnedWettedAreas,
  type ChamberAreas,
  type WallHeatResult,
} from './wall-heat';
