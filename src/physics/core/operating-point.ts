/**
 * User-controllable operating conditions. SI units unless the field name says otherwise.
 */

export type FuelSelection =
  /** Primary reference fuel blend. octaneNumber = liquid-volume % iso-octane in n-heptane (0..100). */
  | { kind: 'PRF'; octaneNumber: number }
  | { kind: 'pure'; species: 'IC8H18' | 'NC7H16' | 'CH4' | 'C3H8' | 'C2H5OH' };

/**
 * Load on the crank in 'free' speed mode (positive torque opposes rotation):
 *  - 'constant': T = loadTorque;
 *  - 'brake': fan / water brake, T = loadTorque·(n/refRpm)^exponent (exponent 2 ≈ fan or hydraulic brake);
 *  - 'vehicle': the engine drives EngineSpec.vehicle through `gear` (clutch engaged, rigid driveline):
 *    road load (rolling + aerodynamic + grade) reflected to the crank, and the vehicle's translating mass
 *    reflected as crank inertia; gear 'neutral' = declutched (no load). loadTorque is ignored.
 */
export type LoadSpec =
  | { kind: 'constant' }
  | { kind: 'brake'; refRpm: number; exponent: number }
  | { kind: 'vehicle'; gear: string; /** Road grade, rise / run. */ grade: number };

export interface OperatingPoint {
  /** 'fixed' = dynamometer/synchronous motor holds speed (CFR); 'free' = crank dynamics with load torque. */
  speedMode: 'fixed' | 'free';
  /** Target (fixed) or initial (free) speed, rev/min. */
  rpm: number;
  /** External load torque on the crank in 'free' mode, N m (positive opposes rotation). */
  loadTorque: number;
  /** Load model in 'free' mode; absent = { kind: 'constant' } (loadTorque). */
  load?: LoadSpec;
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
  /**
   * Spark timing, crank degrees BEFORE firing TDC (positive = advanced). Inductive ignition: the
   * switch-off (spark) angle. Trembler-magneto ignition: the spark LEVER — the timer-contact MAKE
   * angle; the first spark follows after the coil's firing time (an output, not an input).
   */
  sparkAdvanceDeg: number;
  /** Coil dwell time, s (inductive ignition only). */
  dwellTime: number;
  /** Ignition supply of a trembler-magneto system: the dash switch 'MAG' / 'BAT'. Absent: 'magneto'. */
  ignitionSource?: 'magneto' | 'battery';
  /** Geometric compression ratio (variable-CR engines). */
  compressionRatio: number;
  /** External EGR mass fraction of the fresh intake charge, 0..1. */
  egrFraction: number;
  /** Coolant temperature, K (sets wall temperatures together with the WallSpec). */
  coolantTemperature: number;
}
