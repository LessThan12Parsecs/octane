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
  /** Indicated (gas-pressure) torque on the crank, N m. */
  gasTorque: number;
  /** Net crank torque including inertia and friction, N m. */
  netTorque: number;
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
}
