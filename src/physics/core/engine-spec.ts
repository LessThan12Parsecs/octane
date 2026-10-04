/**
 * Engine description consumed by physics AND rendering. Pure data, SI units.
 *
 * Frames (see DESIGN.md §Frames):
 *  - Crank angle θ: 0 = firing TDC, cycle spans [-360°, 360°):
 *      intake -360→-180, compression -180→0, expansion 0→180, exhaust 180→360.
 *    Angles in *Deg fields are crank degrees in this convention. In a multi-cylinder engine θ is the
 *    angle of cylinder 1 (the ENGINE angle); cylinder i runs at its local angle
 *    θ_i = wrap(θ − layout.firingOffsetDeg[i]) and every per-cylinder *Deg field of this spec
 *    (valve timing, …) is in that cylinder's local angle.
 *  - Cylinder frame (CYL): origin at the centre of the cylinder-head fire-deck
 *    face, +y along the cylinder axis pointing away from the piston (toward the
 *    head). Gas occupies x²+z² ≤ (B/2)², −h(θ) ≤ y ≤ 0, where h is the
 *    instantaneous head-to-crown clearance height ('flat-disc' chamber). For an
 *    'l-head' chamber the origin is on the bore axis at the roof of the head cavity
 *    over the bore, h(θ) is the roof-to-crown depth of the bore column, and the
 *    valve pocket beside the bore is a fixed prism (see LHeadChamberSpec).
 *  - World/render frame: +y up, crankshaft axis = world z axis through origin,
 *    +z = FRONT of the engine, cylinder axis parallel to world y. Cylinder i's
 *    frame is the CYL frame translated to (0, ·, layout.axisZ[i]) and, where
 *    layout.mirrorZ[i], mirrored z → −z.
 *  - Crank rotation: clockwise seen from the front (+z), for every engine
 *    (mechanics/kinematics.ts).
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
  /**
   * Valve centre (x, z), cylinder frame, m. Overhead valves: on the head face inside the bore.
   * Side valves ('l-head'): on the block deck beside the bore, inside the valve pocket.
   */
  position: [number, number];
  /** Shroud (masked arc) on the valve head in degrees of arc, 0 = none. CFR intake valve is shrouded. */
  shroudArcDeg: number;
  /** Direction (rad, in x-z plane from +x toward +z) the shroud opening faces; used for swirl. */
  shroudDirection: number;
  /**
   * Cylinder-frame y of the valve-seat plane, m. Absent: 0 (overhead valve seated in the head face).
   * Side valves: the block deck / pocket floor (LHeadChamberSpec.deckY).
   */
  seatY?: number;
  /**
   * Direction the valve head moves when it opens along the cylinder-frame y axis: −1 = down into the
   * cylinder (overhead valve, default), +1 = up out of the block into the head pocket (side valve).
   */
  liftDirection?: -1 | 1;
  /**
   * Running (hot) valve clearance between cam follower and valve stem, m. When present it is the
   * clearance the lift profile subtracts from the cam lift (overrides CycleModelOptions.valveLash).
   */
  lash?: number;
  /** Cam-lobe shape. Absent: the default polydyne (gas-exchange/valve-lift.ts). */
  cam?: CamSpec;
}

/**
 * Cam-lobe (and follower) description. maxLift/openDeg/closeDeg of the ValveSpec stay the quoted
 * valve-lift timing (at `lash`); a geometric cam must reproduce them (tests), and its lobe centre is
 * derived from them.
 */
export type CamSpec =
  /** Even-power polydyne family x(φ) = Σ c_k φ^k (valve-lift.ts polydyneCam); default powers [2, 10, 18, 26]. */
  | { kind: 'polydyne'; powers?: number[] }
  /**
   * Convex three-arc ("circular-arc") lobe acting on a FLAT-faced translating follower (mushroom tappet),
   * no quieting ramp: base circle r_b, flank arcs of radius r_f tangent to the base circle, nose arc r_n;
   * total cam rise = `rise` (the follower lift at the nose, before lash). Lift of a flat follower is the
   * support function of the lobe outline in the follower direction (exact; Norton, Cam Design and
   * Manufacturing Handbook 2009, ch. 6 — UNVERIFIED chapter).
   */
  | { kind: 'three-arc-flat-follower'; baseRadius: number; flankRadius: number; noseRadius: number; rise: number }
  /** Measured zero-lash valve lift table, uniform in crank angle (local cylinder angle), m. */
  | { kind: 'table'; startDeg: number; stepDeg: number; zeroLashLift: readonly number[] };

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
export interface InductiveIgnitionSpec {
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

/**
 * One trembler (vibrator, "buzz") coil: an open-core induction coil whose primary is interrupted by
 * its own magnetic vibrator. While the timer grounds the primary, current builds until the core's pull
 * on the spring-steel armature beats the spring; the armature moves, the points separate after the
 * cushion-spring clearance is taken up, the primary current rings into the condenser and the secondary
 * fires the plug; the armature springs back, the points re-close and the cycle repeats (a buzz of
 * sparks for as long as the timer contact lasts).
 */
export interface TremblerCoilSpec {
  /** Primary inductance with the secondary open, H. */
  primaryInductance: number;
  /** Primary winding resistance, Ω. */
  primaryResistance: number;
  /** Secondary inductance with the primary open, H. */
  secondaryInductance: number;
  /** Secondary winding resistance, Ω. */
  secondaryResistance: number;
  /** Lumped secondary capacitance (winding + HT lead + plug), F. */
  secondaryCapacitance: number;
  /** Magnetic coupling coefficient k (0..1); an open straight core couples loosely. */
  couplingCoefficient: number;
  /** Condenser across the vibrator points, F. */
  condenserCapacitance: number;
  /** Vibrator (armature + points + cushion spring). */
  vibrator: TremblerVibratorSpec;
}

/**
 * 1-DOF vibrator: armature (effective mass m, spring k_s = m(2πf_n)², damping ζ) pulled toward the
 * core by F = K·I₁²/(g₀ − x)² against the spring preload, with K set so that F equals the preload at
 * x = 0 and I₁ = `pullCurrent` (the static pull-in current). The points separate once the armature has
 * moved `breakTravel` (the cushion-spring clearance) and re-close when it returns below it; the
 * armature stops at `maxTravel`.
 */
export interface TremblerVibratorSpec {
  /** Static primary current at which the magnetic pull equals the spring preload, A. */
  pullCurrent: number;
  /** Natural frequency of the armature on its spring, Hz. */
  naturalFrequency: number;
  /** Damping ratio ζ of the armature motion. */
  dampingRatio: number;
  /** Armature travel at which the points separate (cushion-spring clearance), m. */
  breakTravel: number;
  /** Armature travel to the stop, m. */
  maxTravel: number;
  /** Core-to-armature air gap at rest, m. */
  airGap: number;
}

/** Low-tension AC flywheel magneto (permanent magnets on the flywheel sweeping a stationary coil ring). */
export interface MagnetoSpec {
  /** Electrical cycles per crank revolution (Model T: 16 magnets, 16 coils → 8). */
  cyclesPerRevolution: number;
  /** Open-circuit EMF amplitude per crank speed, V_peak per (rad/s): E = k·ω·sin(N·θ + φ). */
  emfConstant: number;
  /**
   * Phase φ of the EMF zero crossing, crank degrees: E = k ω sin(N (θ_crank − phaseDeg)·π/180), with
   * θ_crank the ENGINE angle modulo 360°.
   */
  phaseDeg: number;
  /** Source (coil-ring) resistance, Ω. */
  internalResistance: number;
  /** Source (coil-ring) inductance, H. */
  internalInductance: number;
}

/** Ignition timer (commutator): a roller grounding one coil's primary at a time, driven at cam speed. */
export interface IgnitionTimerSpec {
  /** Crank-angle length of each cylinder's contact, deg. */
  contactArcDeg: number;
  /** Spark-lever range: timer MAKE advance (deg before the cylinder's firing TDC) [full retard, full advance]. */
  advanceRangeDeg: [number, number];
  /** Contact (roller + wiring) resistance, Ω. */
  contactResistance: number;
}

/** Magneto + timer + four trembler coils (Ford Model T), optionally battery-fed. */
export interface TremblerMagnetoIgnitionSpec {
  type: 'trembler-magneto';
  /** The coil fitted to every cylinder. */
  coil: TremblerCoilSpec;
  /**
   * Optional per-cylinder deviations (index = cylinder − 1) — a mis-adjusted vibrator, a weak
   * condenser; absent entries use `coil`.
   */
  coilOverrides?: Partial<Omit<TremblerCoilSpec, 'vibrator'> & { vibrator: Partial<TremblerVibratorSpec> }>[];
  magneto: MagnetoSpec;
  /** Battery for 'battery' operation (OperatingPoint.ignitionSource). */
  battery: { voltage: number; internalResistance: number };
  timer: IgnitionTimerSpec;
}

/** Ignition system variants. */
export type IgnitionSystemSpec = InductiveIgnitionSpec | TremblerMagnetoIgnitionSpec;

/**
 * Side-valve ("L-head", "flathead") chamber, cylinder frame. The chamber is the union of
 *  - the BORE COLUMN: disc x²+z² ≤ (B/2)² from the piston crown (y = −h(θ)) up to the roof of the head
 *    cavity over the bore (y = 0); above the deck plane this is the head cavity, assumed to have the
 *    bore's plan outline;
 *  - the VALVE POCKET: a vertical prism of plan `pocket` (rounded rectangle) MINUS the bore disc, from
 *    the deck plane / valve-seat plane (y = deckY) up to the pocket roof (y = pocket.roofY). The
 *    valves seat in its floor and open upward.
 * Gas flows between them across the shared vertical boundary (the bore circle) for deckY ≤ y ≤ roofY.
 * Fixed chamber volume (pocket) + bore column + crevice must equal the clearance volume of
 * geometry.compressionRatio (engines/*.test.ts checks it).
 */
export interface LHeadChamberSpec {
  /** Cylinder-frame y of the block deck (gasket face of the block) = valve-seat plane, m (< 0). */
  deckY: number;
  /** Height of the piston crown above the deck plane at TDC, m (> 0: the crown rises into the head). */
  crownAboveDeckAtTDC: number;
  /** Valve pocket plan (rounded rectangle, cylinder frame) and roof. */
  pocket: {
    xMin: number;
    xMax: number;
    zMin: number;
    zMax: number;
    /** Plan corner radius, m. */
    cornerRadius: number;
    /** Cylinder-frame y of the pocket roof (underside of the head over the valves), m (deckY < roofY ≤ 0). */
    roofY: number;
  };
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
  /** Allowed compression-ratio range (variable-CR engines), [min, max]. A fixed-CR engine has [CR, CR]. */
  compressionRatioRange: [number, number];
  /** Wrist-pin offset from cylinder axis, m (0 = none). */
  pinOffset: number;
  /** Wrist-pin axis to piston crown distance, m (render + piston top position). */
  compressionHeight: number;
  /**
   * Combustion-chamber shape: 'flat-disc' (flat head + flat crown, the CFR pancake) or 'l-head'
   * (side valves; requires `lHead`).
   */
  chamber: 'flat-disc' | 'l-head';
  /** L-head chamber description (chamber 'l-head'). */
  lHead?: LHeadChamberSpec;
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
   * Block deck / valve-pocket floor temperature, K ('l-head' chambers: the water-jacketed block surface
   * between the bore and the valve seats). Absent: the liner temperature.
   */
  blockTemperature?: number;
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
  /** Block deck / pocket floor ('l-head'). */
  block?: number;
}

/** Butterfly throttle plate (gas-exchange/throttle.ts throttleArea options). */
export interface ThrottlePlateSpec {
  /** Throttle shaft diameter / bore. */
  shaftRatio: number;
  /** Plate angle from the bore cross-section when closed, deg. */
  closedAngleDeg: number;
  /** Leakage area with the plate closed (clearance + idle bypass), m². */
  leakageArea: number;
}

export interface ManifoldSpec {
  /** Intake plenum/runner volume downstream of the throttle, m³. */
  intakeVolume: number;
  /** Exhaust plenum/runner volume upstream of the ambient/back-pressure orifice, m³. */
  exhaustVolume: number;
  /**
   * Throttle bore diameter, m. Without `venturiDiameter` this is the ONLY intake restriction (the CFR:
   * its 9/16 in carburettor venturi, scaled by the throttle opening). With `venturiDiameter` it is the
   * bore of a butterfly throttle (`throttle`) in series downstream of the venturi.
   */
  throttleDiameter: number;
  /** Carburettor venturi throat diameter when the carburettor has a separate butterfly throttle, m. */
  venturiDiameter?: number;
  /** Butterfly plate details (with venturiDiameter). */
  throttle?: ThrottlePlateSpec;
  /** Exhaust outlet (to ambient) effective diameter, m. */
  exhaustOutletDiameter: number;
  /** Intake runner diameter at the port, m (for inflow velocity / turbulence production). */
  intakePortDiameter: number;
  exhaustPortDiameter: number;
  /** Heated length of the intake port (runner in the head), m — intake-port heat transfer (optional). */
  intakePortLength?: number;
}

/** Arrangement of the cylinders of a multi-cylinder engine. */
export interface CylinderLayoutSpec {
  /** Firing order, 1-based cylinder numbers (cylinder 1 = front). */
  firingOrder: number[];
  /**
   * Firing-TDC offset of each cylinder (index = cylinder − 1) after cylinder 1's firing TDC, crank deg,
   * in [0, 720): θ_i = wrap(θ − firingOffsetDeg[i]). Inline four 1-2-4-3: [0, 180, 540, 360].
   */
  firingOffsetDeg: number[];
  /** Cylinder-axis position along the crank axis (world z, +z = front), m. */
  axisZ: number[];
  /**
   * Render-only: true where the cylinder's valve/port arrangement is the mirror image (z → −z) of the
   * spec's (e.g. side-valve order E-I-I-E-E-I-I-E). Physics is mirror-invariant.
   */
  mirrorZ: boolean[];
  /** Main bearings: crank-axis positions (world z), m. */
  mainBearingZ: number[];
}

/** Vehicle driven by the engine ('vehicle' load model). */
export interface VehicleSpec {
  /** Laden vehicle mass, kg. */
  mass: number;
  /** Driven-wheel rolling radius, m. */
  wheelRadius: number;
  /** Final-drive (axle) ratio. */
  finalDrive: number;
  /** Transmission ratios (engine/driveshaft) by gear name. */
  gears: Record<string, number>;
  /** Driveline efficiency (engine → wheel), 0..1. */
  drivelineEfficiency: number;
  /** Rolling-resistance coefficient. */
  rollingResistance: number;
  /** Drag area C_d·A, m². */
  dragArea: number;
}

export interface EngineSpec {
  /** Registry key (engines/index.ts), e.g. 'cfr-f1', 'ford-model-t'. Absent: an ad-hoc spec. */
  id?: string;
  /** Display name. */
  name: string;
  /** Literature sources for the numbers in this spec (free text list). */
  sources: string[];
  cycle: 'four-stroke-si';
  cylinders: number;
  /** Multi-cylinder arrangement. Absent: one cylinder at the origin (firing offset 0). */
  layout?: CylinderLayoutSpec;
  geometry: EngineGeometrySpec;
  sparkPlug: SparkPlugSpec;
  intakeValve: ValveSpec;
  exhaustValve: ValveSpec;
  manifolds: ManifoldSpec;
  masses: MassSpec;
  walls: WallSpec;
  ignition: IgnitionSystemSpec;
  /** The vehicle the engine drives (enables the 'vehicle' load model). */
  vehicle?: VehicleSpec;
}

/** Firing-TDC offsets of all cylinders, deg ([0] for a single-cylinder spec). */
export function firingOffsetsDeg(spec: EngineSpec): readonly number[] {
  return spec.layout?.firingOffsetDeg ?? SINGLE_OFFSET;
}
const SINGLE_OFFSET: readonly number[] = Object.freeze([0]);

/** Local crank angle of cylinder `i` (0-based) at engine angle θ, wrapped to [-360, 360). */
export function cylinderAngleDeg(spec: EngineSpec, i: number, thetaDeg: number): number {
  const off = spec.layout?.firingOffsetDeg[i] ?? 0;
  let t = thetaDeg - off;
  t = ((((t + 360) % 720) + 720) % 720) - 360;
  return t;
}

/** True when the compression ratio can be changed at run time (range not degenerate). */
export function hasVariableCompressionRatio(spec: EngineSpec): boolean {
  const [lo, hi] = spec.geometry.compressionRatioRange;
  return hi > lo;
}

/** Lumped secondary (HT side) capacitance of the spec's ignition coil, F. */
export function ignitionSecondaryCapacitance(ig: IgnitionSystemSpec): number {
  return ig.type === 'inductive' ? ig.secondaryCapacitance : ig.coil.secondaryCapacitance;
}

/** The spec's ignition as an inductive system; throws for other ignition types (inductive-only consumers). */
export function requireInductiveIgnition(spec: EngineSpec): InductiveIgnitionSpec {
  const ig = spec.ignition;
  if (ig.type !== 'inductive') throw new Error(`${spec.name}: '${ig.type}' ignition where an inductive coil is required`);
  return ig;
}
