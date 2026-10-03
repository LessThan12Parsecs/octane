/**
 * User-controllable operating conditions. SI units unless the field name says otherwise.
 */

export type FuelSelection =
  /** Primary reference fuel blend. octaneNumber = liquid-volume % iso-octane in n-heptane (0..100). */
  | { kind: 'PRF'; octaneNumber: number }
  | { kind: 'pure'; species: 'IC8H18' | 'NC7H16' | 'CH4' | 'C3H8' | 'C2H5OH' };

export interface OperatingPoint {
  /** 'fixed' = dynamometer/synchronous motor holds speed (CFR); 'free' = crank dynamics with load torque. */
  speedMode: 'fixed' | 'free';
  /** Target (fixed) or initial (free) speed, rev/min. */
  rpm: number;
  /** External load torque on the crank in 'free' mode, N m (positive opposes rotation). */
  loadTorque: number;
  /** Throttle opening 0 (closed, leakage only) .. 1 (wide open). */
  throttle: number;
  /** Ambient (and exhaust back-) pressure, Pa. */
  ambientPressure: number;
  /** Ambient air temperature, K. */
  ambientTemperature: number;
  /** Relative humidity of intake air, 0..1. */
  relativeHumidity: number;
  /** Fresh-charge temperature at the intake port (after heater/carburettor), K. */
  intakeMixtureTemperature: number;
  fuel: FuelSelection;
  /** Fuel/air equivalence ratio φ of the fresh charge. */
  equivalenceRatio: number;
  /** Spark timing, crank degrees BEFORE firing TDC (positive = advanced). */
  sparkAdvanceDeg: number;
  /** Coil dwell time, s. */
  dwellTime: number;
  /** Geometric compression ratio (variable-CR engines). */
  compressionRatio: number;
  /** External EGR mass fraction of the fresh intake charge, 0..1. */
  egrFraction: number;
  /** Coolant temperature, K (sets wall temperatures together with the WallSpec). */
  coolantTemperature: number;
}
