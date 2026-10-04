/**
 * Render/UI contract. The simulator (in a Web Worker) emits a stream of
 * EngineSnapshot at a fixed simulated-time or crank-angle cadence plus one
 * CycleSummary per completed cycle. The renderer only ever reads these —
 * it never calls physics code directly (except pure geometry helpers).
 *
 * Plain JSON-able data (no typed arrays) so it survives structured cloning cheaply.
 * All SI units; *Deg fields are crank degrees in the firing-TDC convention.
 */

export type SparkPhase = 'off' | 'charging' | 'breakdown' | 'arc' | 'glow' | 'done';

export type CylinderPhase = 'gas-exchange' | 'compression' | 'combustion' | 'expansion';

/**
 * Multi-cylinder engines: the top-level per-cylinder fields of EngineSnapshot (kinematics, in-cylinder
 * state, flame, spark, gas exchange, knock, composition) describe CYLINDER 1, and `cylinders[i]` holds
 * every cylinder (including cylinder 1, identical to the top level) at its own local crank angle.
 * Engine-level fields (t, cycle, thetaDeg = engine angle = cylinder 1's angle, rpm, manifold
 * pressures, gasTorque, netTorque, frictionTorque, loadTorque, magnetoEmf) exist only at the top level.
 */
export type CylinderSnapshot = Pick<
  EngineSnapshot,
  | 'pistonDisplacement' | 'clearanceHeight' | 'rodAngle' | 'intakeLift' | 'exhaustLift'
  | 'phase' | 'volume' | 'pressure' | 'temperatureMean' | 'temperatureUnburned' | 'temperatureBurned'
  | 'massFractionBurned' | 'mass' | 'heatReleaseRate' | 'heatLossRate'
  | 'flame' | 'spark' | 'intakeMassFlow' | 'exhaustMassFlow' | 'knock' | 'burnedComposition'
> & {
  /** 0-based cylinder index (cylinder number − 1). */
  index: number;
  /** Local crank angle of this cylinder, deg, [-360, 360), 0 = its firing TDC. */
  thetaDeg: number;
  /** Completed cycles of this cylinder (its local angle wraps at +360°). */
  cycle: number;
  /** This cylinder's gas-pressure torque on the crank, N m. */
  gasTorque: number;
};

export interface EngineSnapshot {
  /** Simulated time since start, s. */
  t: number;
  /** Completed cycles since start. */
  cycle: number;
  /** Crank angle, degrees, [-360, 360), 0 = firing TDC. */
  thetaDeg: number;
  /** Instantaneous speed, rev/min. */
  rpm: number;

  // ---- kinematics (for rendering the mechanism) ----
  /** Piston crown displacement below its TDC position, m. */
  pistonDisplacement: number;
  /** Instantaneous head-to-crown clearance height h(θ), m. */
  clearanceHeight: number;
  /** Connecting-rod angle from the cylinder axis, rad. */
  rodAngle: number;
  /** Intake / exhaust valve lift, m. */
  intakeLift: number;
  exhaustLift: number;

  // ---- in-cylinder thermodynamic state ----
  phase: CylinderPhase;
  /** Cylinder volume, m³. */
  volume: number;
  /** Cylinder pressure, Pa (including synthesised knock oscillation). */
  pressure: number;
  /** Mass-averaged gas temperature, K. */
  temperatureMean: number;
  /** Unburned-zone temperature, K. */
  temperatureUnburned: number;
  /** Burned-zone temperature, K (0 when no burned zone). */
  temperatureBurned: number;
  /** Burned mass fraction (0..1). */
  massFractionBurned: number;
  /** Total trapped mass, kg. */
  mass: number;
  /** Instantaneous chemical heat release rate, W. */
  heatReleaseRate: number;
  /** Total wall heat-loss rate, W (positive = gas → wall). */
  heatLossRate: number;

  // ---- flame ----
  flame: {
    /** 'none' | 'kernel' (spark-kernel growth) | 'turbulent' (fully developed) | 'burnout' | 'done' */
    stage: 'none' | 'kernel' | 'turbulent' | 'burnout' | 'done';
    /** Flame (entrainment) front radius from the gap centre, m. */
    radius: number;
    /** Centre of the flame sphere, cylinder frame, m (spark gap centre, possibly convected). */
    center: [number, number, number];
    /** Flame front area inside the chamber, m². */
    area: number;
    /** Laminar and turbulent burning velocities at the front, m/s. */
    laminarSpeed: number;
    turbulentSpeed: number;
    /** Turbulence intensity u', m/s. */
    turbulenceIntensity: number;
  };

  // ---- spark ----
  spark: {
    phase: SparkPhase;
    /** Primary current, A. */
    primaryCurrent: number;
    /** Secondary (gap) voltage magnitude, V. */
    secondaryVoltage: number;
    /** Secondary (gap) current, A. */
    secondaryCurrent: number;
    /** Cumulative electrical energy delivered to the gap this cycle, J. */
    energyDelivered: number;
    /** Breakdown voltage required at current gas density, V. */
    breakdownVoltage: number;
    /**
     * Trembler / multi-spark ignition (optional): gap breakdowns so far in the current ignition event
     * (one timer contact), cumulative; points open now; timer contact closed now; local crank angle of
     * the FIRST breakdown of the event (NaN before it).
     */
    breakdownCount?: number;
    /** Re-ignitions of a spark at the current zeros of its condenser ring (not counted as breakdowns). Optional. */
    reignitionCount?: number;
    pointsOpen?: boolean;
    timerClosed?: boolean;
    firstSparkDeg?: number;
    /** Primary (points) voltage, V (optional). */
    primaryVoltage?: number;
  };

  // ---- gas exchange ----
  /** Mass flow through intake / exhaust valves, kg/s (positive = into cylinder for intake, out of cylinder for exhaust). */
  intakeMassFlow: number;
  exhaustMassFlow: number;
  /** Intake / exhaust manifold pressure, Pa. */
  intakeManifoldPressure: number;
  exhaustManifoldPressure: number;

  // ---- knock ----
  knock: {
    /** Livengood–Wu integral ∫dt/τ in the end gas (autoignition at 1). */
    integral: number;
    /** True once end-gas autoignition occurred this cycle. */
    autoignited: boolean;
    /** Synthesised acoustic pressure oscillation added to `pressure`, Pa. */
    oscillation: number;
  };

  // ---- composition (burned zone, mole fractions; unburned when no burned zone) ----
  burnedComposition: {
    CO2: number; H2O: number; CO: number; O2: number; H2: number; OH: number;
    H: number; O: number; NO: number; N2: number;
  };

  // ---- mechanics ----
  /** Indicated (gas-pressure) torque on the crank, N m — the whole engine (sum over cylinders). */
  gasTorque: number;
  /** Net crank torque including inertia and friction, N m. */
  netTorque: number;
  /** Friction torque on the crank, N m (positive = opposing rotation). Optional. */
  frictionTorque?: number;
  /**
   * Torque absorbed by the load, N m (positive = opposing rotation): the load model in 'free' mode,
   * the dynamometer holding torque in 'fixed' mode. Optional.
   */
  loadTorque?: number;
  /** Vehicle road speed, m/s ('vehicle' load). Optional. */
  vehicleSpeed?: number;

  // ---- multi-cylinder / magneto (optional) ----
  /** Every cylinder at its local angle (multi-cylinder engines; absent for single-cylinder specs). */
  cylinders?: CylinderSnapshot[];
  /** Index of the cylinder whose ignition timer contact is closed (-1 = none). Optional. */
  firingCylinder?: number;
  /** Magneto open-circuit EMF, V (trembler-magneto ignition). Optional. */
  magnetoEmf?: number;
}

/** Engine-level results over one engine cycle (720° of the engine angle). */
export interface EngineCycleSummary {
  /** Mean speed, rev/min. */
  rpmMean: number;
  /** Mean indicated (gas) torque, N m (net, all cylinders). */
  indicatedTorque: number;
  /** Mean friction torque, N m (positive = loss). */
  frictionTorque: number;
  /**
   * Mean brake torque = indicated − friction − (inertia, ≈ 0 at steady speed), N m: the ENGINE's output at the
   * crankshaft. The inertia term is the kinetic-energy change of the engine's own rotating and reciprocating
   * parts (crank + flywheel and the piston/rod mechanisms) over the cycle / 4π; a driven car's inertia is
   * downstream of the output (see loadInertiaTorque) and the impulse of a clutch engagement is excluded.
   */
  brakeTorque: number;
  /** Brake power, W. */
  brakePower: number;
  /** Brake / indicated (net, mean over cylinders) / friction mean effective pressure, Pa. */
  bmep: number;
  imepNet: number;
  fmep: number;
  /** Air and fuel mass flow through the carburettor, kg/s. */
  airMassFlow: number;
  fuelMassFlow: number;
  /** Volumetric efficiency over the total displacement (ambient-referenced). */
  volumetricEfficiency: number;
  /** Brake specific fuel consumption, kg/J, and brake thermal efficiency (LHV). */
  bsfc: number;
  brakeEfficiency: number;
  /**
   * Mean load torque, N m: free speed — the load model's ∫T_L ω dt / 4π (a car: road load + driveline loss
   * while in gear, 0 in neutral); fixed speed — the dynamometer's absorbed torque (= brake). Mean vehicle speed,
   * m/s ('vehicle' load in free speed: road distance / cycle time, coasting in neutral included).
   */
  loadTorque: number;
  vehicleSpeed?: number;
  /**
   * Free speed, engines with a vehicle (EngineSpec.vehicle): the kinetic-energy change of the car while coupled
   * to the crank (its inertia reflected through the engaged gear, ½J_L Δω² between gear changes) / 4π, N m;
   * the energy balance closes as brakeTorque ≈ loadTorque + loadInertiaTorque. Optional.
   */
  loadInertiaTorque?: number;
  /** Energy dissipated in the clutch by gear engagements during the cycle (free speed), J. Optional (as loadInertiaTorque). */
  clutchLoss?: number;
}

export interface CycleSummary {
  cycle: number;
  /** Gross IMEP (compression + expansion strokes), Pa. */
  imepGross: number;
  /** Net IMEP (full 720°), Pa. */
  imepNet: number;
  /** Pumping MEP, Pa. */
  pmep: number;
  /** Peak pressure, Pa, and its location, crank deg. */
  peakPressure: number;
  peakPressureDeg: number;
  /** Max rate of pressure rise, Pa/deg. */
  maxPressureRiseRate: number;
  /** 10 / 50 / 90 % mass-fraction-burned angles, crank deg (NaN if not reached). */
  ca10: number;
  ca50: number;
  ca90: number;
  /** Net indicated thermal efficiency (fuel LHV basis). */
  indicatedEfficiency: number;
  /** Indicated specific fuel consumption, kg/J. */
  isfc: number;
  /** Mass trapped at IVC, kg; residual mass fraction; volumetric efficiency (ambient-referenced). */
  trappedMass: number;
  residualFraction: number;
  volumetricEfficiency: number;
  /** Fuel mass burned this cycle, kg. */
  fuelMass: number;
  /** Exhaust NO at EVO, ppm (wet, mole basis). */
  noPpm: number;
  /** Exhaust CO at EVO, mole fraction. */
  coFraction: number;
  /**
   * Knock: onset angle (NaN if none), autoigniting ("end-gas") mass fraction at onset — all unburned gas
   * by default, the gas ahead of the flame front with CycleModelOptions.knockBrushAutoignition false —
   * and max amplitude of pressure oscillation (MAPO), Pa.
   */
  knockOnsetDeg: number;
  knockEndGasFraction: number;
  mapo: number;
  /** True if the spark kernel failed to develop into a flame. */
  misfire: boolean;
  /** Wall heat loss over the closed cycle, J; gross indicated work, J. */
  heatLoss: number;
  indicatedWorkGross: number;
  /** 0-based cylinder this summary belongs to (multi-cylinder engines; absent = 0). */
  cylinder?: number;
  /** Local crank angle of the first spark (gap breakdown) of the cycle, deg (NaN: none). Optional. */
  sparkDeg?: number;
  /** Gap breakdowns (distinct sparks) in the cycle's ignition event. Optional. */
  sparkCount?: number;
  /** Re-ignitions at current zeros in the cycle's ignition event (trembler; not sparks). Optional. */
  reignitionCount?: number;
  /** Engine-level results, attached to cylinder 1's summary at the end of each engine cycle. Optional. */
  engine?: EngineCycleSummary;
}
