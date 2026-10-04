/**
 * Zero-allocation per-cylinder view of an EngineSnapshot (pure TypeScript).
 *
 * Multi-cylinder snapshots carry every cylinder in `cylinders[i]` (snapshot.ts CylinderSnapshot: its
 * local crank angle and cycle count, kinematics, in-cylinder state, flame, spark, knock, composition);
 * the top-level per-cylinder fields describe cylinder 1. `select(s, i)` returns an EngineSnapshot whose
 * per-cylinder fields are cylinder i's and whose engine-level fields (t, rpm, manifold pressures,
 * torques, magneto) are the engine's — by copying numbers and object REFERENCES into one reused
 * object, so the in-cylinder visuals read it exactly like a single-cylinder snapshot. Without
 * `cylinders[i]` (single-cylinder specs, or a worker that only reports cylinder 1) it returns `s`.
 */
import type { EngineSnapshot } from '../../physics/core/snapshot';

export class CylinderSnapshotView {
  /** The reused view object (valid until the next select). */
  readonly view: EngineSnapshot;

  constructor() {
    this.view = {
      t: 0, cycle: 0, thetaDeg: 0, rpm: 0,
      pistonDisplacement: 0, clearanceHeight: 0, rodAngle: 0, intakeLift: 0, exhaustLift: 0,
      phase: 'compression', volume: 0, pressure: 0, temperatureMean: 0, temperatureUnburned: 0,
      temperatureBurned: 0, massFractionBurned: 0, mass: 0, heatReleaseRate: 0, heatLossRate: 0,
      flame: { stage: 'none', radius: 0, center: [0, 0, 0], area: 0, laminarSpeed: 0, turbulentSpeed: 0, turbulenceIntensity: 0 },
      spark: { phase: 'off', primaryCurrent: 0, secondaryVoltage: 0, secondaryCurrent: 0, energyDelivered: 0, breakdownVoltage: 0 },
      intakeMassFlow: 0, exhaustMassFlow: 0, intakeManifoldPressure: 0, exhaustManifoldPressure: 0,
      knock: { integral: 0, autoignited: false, oscillation: 0 },
      burnedComposition: { CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0 },
      gasTorque: 0, netTorque: 0,
    };
  }

  /** Cylinder `i` (0-based) of `s`, or `s` itself when it has no per-cylinder record for i. */
  select(s: EngineSnapshot, i: number): EngineSnapshot {
    const c = s.cylinders?.[i];
    if (!c) return s;
    const v = this.view;
    // engine level
    v.t = s.t;
    v.rpm = s.rpm;
    v.intakeManifoldPressure = s.intakeManifoldPressure;
    v.exhaustManifoldPressure = s.exhaustManifoldPressure;
    v.netTorque = s.netTorque;
    v.frictionTorque = s.frictionTorque;
    v.loadTorque = s.loadTorque;
    v.vehicleSpeed = s.vehicleSpeed;
    v.firingCylinder = s.firingCylinder;
    v.magnetoEmf = s.magnetoEmf;
    v.cylinders = undefined;
    // this cylinder (local crank angle and cycle count)
    v.cycle = c.cycle;
    v.thetaDeg = c.thetaDeg;
    v.gasTorque = c.gasTorque;
    v.pistonDisplacement = c.pistonDisplacement;
    v.clearanceHeight = c.clearanceHeight;
    v.rodAngle = c.rodAngle;
    v.intakeLift = c.intakeLift;
    v.exhaustLift = c.exhaustLift;
    v.phase = c.phase;
    v.volume = c.volume;
    v.pressure = c.pressure;
    v.temperatureMean = c.temperatureMean;
    v.temperatureUnburned = c.temperatureUnburned;
    v.temperatureBurned = c.temperatureBurned;
    v.massFractionBurned = c.massFractionBurned;
    v.mass = c.mass;
    v.heatReleaseRate = c.heatReleaseRate;
    v.heatLossRate = c.heatLossRate;
    v.flame = c.flame;
    v.spark = c.spark;
    v.intakeMassFlow = c.intakeMassFlow;
    v.exhaustMassFlow = c.exhaustMassFlow;
    v.knock = c.knock;
    v.burnedComposition = c.burnedComposition;
    return v;
  }
}
