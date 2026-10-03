/**
 * Engine description consumed by physics AND rendering. Pure data, SI units.
 *
 * Frames (see DESIGN.md §Frames):
 *  - Crank angle θ: 0 = firing TDC, cycle spans [-360°, 360°):
 *      intake -360→-180, compression -180→0, expansion 0→180, exhaust 180→360.
 *    Angles in *Deg fields are crank degrees in this convention.
 *  - Cylinder frame (CYL): origin at the centre of the cylinder-head fire-deck
 *    face, +y along the cylinder axis pointing away from the piston (toward the
 *    head). Gas occupies x²+z² ≤ (B/2)², −h(θ) ≤ y ≤ 0, where h is the
 *    instantaneous head-to-crown clearance height.
 *  - World/render frame: +y up, crankshaft axis = world z axis through origin,
 *    cylinder axis = world y axis.
 */

export interface ValveSpec {
  /** Number of valves of this kind per cylinder. */
  count: number;
  /** Outer diameter of the valve head, m. */
  headDiameter: number;
  /** Inner seat diameter D_v (reference for L/D and curtain area), m. */
  seatInnerDiameter: number;
  /** Seat angle measured from the valve face plane, rad (typ. 45°). */
  seatAngle: number;
  /** Valve stem diameter, m. */
  stemDiameter: number;
  /** Maximum lift, m. */
  maxLift: number;
  /** Opening angle (crank deg, firing-TDC convention) at `timingLiftThreshold`. */
  openDeg: number;
  /** Closing angle (crank deg, firing-TDC convention) at `timingLiftThreshold`. */
  closeDeg: number;
  /** Lift at which open/close timing is quoted (0 = true seat contact), m. */
  timingLiftThreshold: number;
  /** Valve centre (x, z) on the head face, cylinder frame, m. */
  position: [number, number];
  /** Shroud (masked arc) on the valve head in degrees of arc, 0 = none. CFR intake valve is shrouded. */
  shroudArcDeg: number;
  /** Direction (rad, in x-z plane from +x toward +z) the shroud opening faces; used for swirl. */
  shroudDirection: number;
}

export interface SparkPlugSpec {
  /** Centre of the spark gap in the cylinder frame, m (y ≤ 0: inside the chamber). */
  gapCenter: [number, number, number];
  /** Electrode gap, m. */
  gap: number;
  /** Centre-electrode diameter, m. */
  centerElectrodeDiameter: number;
  /** Ground-electrode width, m. */
  groundElectrodeWidth: number;
  /** Unit vector of the spark-plug axis in cylinder frame (pointing into the chamber). */
  axis: [number, number, number];
  /** Thread (shell) diameter, m (render + head bore). */
  threadDiameter: number;
}

/** Inductive (Kettering / transistorised) ignition circuit. */
export interface IgnitionSystemSpec {
  type: 'inductive';
  /** Supply (battery/alternator) voltage, V. */
  supplyVoltage: number;
  /** Primary winding inductance, H. */
  primaryInductance: number;
  /** Total primary-circuit resistance (winding + ballast + switch), Ω. */
  primaryResistance: number;
  /** Secondary winding inductance, H. */
  secondaryInductance: number;
  /** Total secondary resistance (winding + suppressor), Ω. */
  secondaryResistance: number;
  /** Lumped secondary capacitance (coil + HT lead + plug), F. */
  secondaryCapacitance: number;
  /** Magnetic coupling coefficient k (0..1). */
  couplingCoefficient: number;
  /** Primary current limit of the driver (Inf if none), A. */
  primaryCurrentLimit: number;
  /** Default dwell time, s. */
  dwellTime: number;
}

export interface EngineGeometrySpec {
  /** Cylinder bore, m. */
  bore: number;
  /** Stroke, m (= 2 × crank radius). */
  stroke: number;
  /** Connecting-rod length, centre to centre, m. */
  conRodLength: number;
  /** Nominal geometric compression ratio (can be overridden by the operating point for variable-CR engines). */
  compressionRatio: number;
  /** Allowed compression-ratio range (variable-CR engines), [min, max]. */
  compressionRatioRange: [number, number];
  /** Wrist-pin offset from cylinder axis, m (0 = none). */
  pinOffset: number;
  /** Wrist-pin axis to piston crown distance, m (render + piston top position). */
  compressionHeight: number;
  /** Combustion-chamber shape. Only the flat head + flat crown disc is modelled for now. */
  chamber: 'flat-disc';
  /** Top-land crevice volume (piston/ring/liner), m³. Included in clearance volume accounting. */
  creviceVolume: number;
  /**
   * Part of creviceVolume made of NARROW crevices the flame cannot enter (piston top land, ring
   * grooves), m³ — the crevice-flow zone of the cycle model (gas at wall temperature, filled and emptied
   * with the cylinder pressure). Optional (0 / absent: no crevice zone); the rest of creviceVolume
   * (plug and pickup cavities) stays part of the chamber.
   */
  quenchCreviceVolume?: number;
}

export interface MassSpec {
  /** Piston + rings + wrist pin, kg. */
  piston: number;
  /** Connecting rod total mass, kg. */
  conRod: number;
  /** Distance of rod centre of gravity from the big-end (crank pin) centre, m. */
  conRodCgFromBigEnd: number;
  /** Rod moment of inertia about its CG, kg m². */
  conRodInertiaCg: number;
  /** Crankshaft + flywheel + everything rigidly rotating with the crank (incl. dyno coupling), kg m². */
  rotatingInertia: number;
}

export interface WallSpec {
  /** Cylinder-head fire-deck temperature, K. */
  headTemperature: number;
  /** Piston crown temperature, K. */
  pistonTemperature: number;
  /** Liner temperature, K. */
  linerTemperature: number;
  /** Valve-face temperatures, K. */
  intakeValveTemperature: number;
  exhaustValveTemperature: number;
  /** Intake-port wall temperature, K (intake-port heat transfer; optional). */
  intakePortTemperature?: number;
  /**
   * Optional lumped wall model: thermal resistance from each surface to the coolant, K/W. With it the
   * cycle model sets each surface to T_coolant + R_i·Q̄_i (Q̄_i = cycle-mean gas-to-surface heat flow),
   * so the wall temperatures follow the operating point (fired/motored, load, speed); the temperatures
   * above are then the reference state at which R_i were fitted (and the initial values).
   */
  thermalResistance?: WallThermalResistance;
  /** Coolant temperature of the reference state of the temperatures above, K (lumped model). */
  referenceCoolantTemperature?: number;
}

/** Surface-to-coolant thermal resistances of the lumped wall model, K/W. */
export interface WallThermalResistance {
  head: number;
  piston: number;
  liner: number;
  intakeValve: number;
  exhaustValve: number;
}

export interface ManifoldSpec {
  /** Intake plenum/runner volume downstream of the throttle, m³. */
  intakeVolume: number;
  /** Exhaust plenum/runner volume upstream of the ambient/back-pressure orifice, m³. */
  exhaustVolume: number;
  /** Throttle bore diameter, m. */
  throttleDiameter: number;
  /** Exhaust outlet (to ambient) effective diameter, m. */
  exhaustOutletDiameter: number;
  /** Intake runner diameter at the port, m (for inflow velocity / turbulence production). */
  intakePortDiameter: number;
  exhaustPortDiameter: number;
  /** Heated length of the intake port (runner in the head), m — intake-port heat transfer (optional). */
  intakePortLength?: number;
}

export interface EngineSpec {
  /** Display name. */
  name: string;
  /** Literature sources for the numbers in this spec (free text list). */
  sources: string[];
  cycle: 'four-stroke-si';
  cylinders: number;
  geometry: EngineGeometrySpec;
  sparkPlug: SparkPlugSpec;
  intakeValve: ValveSpec;
  exhaustValve: ValveSpec;
  manifolds: ManifoldSpec;
  masses: MassSpec;
  walls: WallSpec;
  ignition: IgnitionSystemSpec;
}
