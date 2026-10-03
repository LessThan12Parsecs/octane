/** Test-only factories for EngineSnapshot / CycleSummary (used by src/ui/*.test.ts). */
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';

export function makeSnapshot(over: Partial<EngineSnapshot> = {}): EngineSnapshot {
  const base: EngineSnapshot = {
    t: 0,
    cycle: 0,
    thetaDeg: -360,
    rpm: 600,
    pistonDisplacement: 0,
    clearanceHeight: 0.01,
    rodAngle: 0,
    intakeLift: 0,
    exhaustLift: 0,
    phase: 'gas-exchange',
    volume: 100e-6,
    pressure: 1e5,
    temperatureMean: 350,
    temperatureUnburned: 350,
    temperatureBurned: 0,
    massFractionBurned: 0,
    mass: 5e-4,
    heatReleaseRate: 0,
    heatLossRate: 0,
    flame: { stage: 'none', radius: 0, center: [0, 0, 0], area: 0, laminarSpeed: 0, turbulentSpeed: 0, turbulenceIntensity: 0 },
    spark: { phase: 'off', primaryCurrent: 0, secondaryVoltage: 0, secondaryCurrent: 0, energyDelivered: 0, breakdownVoltage: 10e3 },
    intakeMassFlow: 0,
    exhaustMassFlow: 0,
    intakeManifoldPressure: 1e5,
    exhaustManifoldPressure: 1e5,
    knock: { integral: 0, autoignited: false, oscillation: 0 },
    burnedComposition: { CO2: 0, H2O: 0, CO: 0, O2: 0.21, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0.79 },
    gasTorque: 0,
    netTorque: 0,
  };
  return { ...base, ...over };
}

/**
 * Snapshots at a fixed crank-angle cadence starting at (cycle, θ0), advancing time
 * consistently with rpm. `fn` may customise each snapshot.
 */
export function makeRun(opts: {
  n: number;
  dDeg: number;
  theta0?: number;
  t0?: number;
  cycle0?: number;
  rpm?: number;
  fn?: (s: EngineSnapshot, i: number) => void;
}): EngineSnapshot[] {
  const rpm = opts.rpm ?? 600;
  let theta = opts.theta0 ?? -360;
  let cycle = opts.cycle0 ?? 0;
  let t = opts.t0 ?? 0;
  const out: EngineSnapshot[] = [];
  for (let i = 0; i < opts.n; i++) {
    const s = makeSnapshot({ t, cycle, thetaDeg: theta, rpm });
    opts.fn?.(s, i);
    out.push(s);
    theta += opts.dDeg;
    t += opts.dDeg / (6 * rpm);
    if (theta >= 360) {
      theta -= 720;
      cycle++;
    }
  }
  return out;
}

export function makeCycle(over: Partial<CycleSummary> = {}): CycleSummary {
  return {
    cycle: 0,
    imepGross: 10e5,
    imepNet: 9.5e5,
    pmep: -0.5e5,
    peakPressure: 35e5,
    peakPressureDeg: 15,
    maxPressureRiseRate: 2e5,
    ca10: 2,
    ca50: 10,
    ca90: 25,
    indicatedEfficiency: 0.32,
    isfc: 7e-8,
    trappedMass: 6e-4,
    residualFraction: 0.05,
    volumetricEfficiency: 0.85,
    fuelMass: 4e-5,
    noPpm: 1500,
    coFraction: 0.005,
    knockOnsetDeg: NaN,
    knockEndGasFraction: NaN,
    mapo: 0,
    misfire: false,
    heatLoss: 80,
    indicatedWorkGross: 600,
    ...over,
  };
}
