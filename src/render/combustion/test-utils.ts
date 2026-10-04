/**
 * Test fixtures for the combustion visuals (a fixed CFR-like spec independent
 * of the concurrently refined CFR_F1 numbers, and a snapshot factory).
 * Only imported by *.test.ts files.
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { CylinderSnapshot, EngineSnapshot } from '../../physics/core/snapshot';
import { MODEL_T } from '../../physics/engines/model-t';

const IN = 0.0254;
export const TEST_BORE = 3.25 * IN;
export const TEST_R = TEST_BORE / 2;

export function testSpec(): EngineSpec {
  return {
    name: 'test CFR-like',
    sources: [],
    cycle: 'four-stroke-si',
    cylinders: 1,
    geometry: {
      bore: TEST_BORE, stroke: 4.5 * IN, conRodLength: 10 * IN, compressionRatio: 7,
      compressionRatioRange: [4, 18], pinOffset: 0, compressionHeight: 2 * IN, chamber: 'flat-disc', creviceVolume: 0,
    },
    sparkPlug: {
      gapCenter: [0, -0.002, TEST_R - 0.004], gap: 0.5e-3, centerElectrodeDiameter: 2.5e-3,
      groundElectrodeWidth: 2.5e-3, axis: [0, 0, -1], threadDiameter: 18e-3,
    },
    intakeValve: {
      count: 1, headDiameter: 1.346 * IN, seatInnerDiameter: 1.2 * IN, seatAngle: Math.PI / 4, stemDiameter: 0.34 * IN,
      maxLift: 0.25 * IN, openDeg: -345, closeDeg: -152, timingLiftThreshold: 0, position: [-0.02, 0],
      shroudArcDeg: 180, shroudDirection: Math.PI / 2,
    },
    exhaustValve: {
      count: 1, headDiameter: 1.356 * IN, seatInnerDiameter: 1.2 * IN, seatAngle: Math.PI / 4, stemDiameter: 0.34 * IN,
      maxLift: 0.25 * IN, openDeg: 141, closeDeg: -347, timingLiftThreshold: 0, position: [0.02, 0],
      shroudArcDeg: 0, shroudDirection: 0,
    },
    manifolds: {
      intakeVolume: 2e-3, exhaustVolume: 3e-3, throttleDiameter: IN, exhaustOutletDiameter: 1.5 * IN,
      intakePortDiameter: 1.2 * IN, exhaustPortDiameter: 1.2 * IN,
    },
    masses: { piston: 1.2, conRod: 1.5, conRodCgFromBigEnd: 0.08, conRodInertiaCg: 0.01, rotatingInertia: 2 },
    walls: { headTemperature: 450, pistonTemperature: 500, linerTemperature: 400, intakeValveTemperature: 500, exhaustValveTemperature: 900 },
    ignition: {
      type: 'inductive', supplyVoltage: 13.5, primaryInductance: 4e-3, primaryResistance: 1.5, secondaryInductance: 20,
      secondaryResistance: 5e3, secondaryCapacitance: 60e-12, couplingCoefficient: 0.98, primaryCurrentLimit: 8, dwellTime: 3e-3,
    },
  };
}

type DeepPartial<T> = { [K in keyof T]?: T[K] extends object ? (T[K] extends unknown[] ? T[K] : DeepPartial<T[K]>) : T[K] };

/** A mid-compression, no-flame snapshot; override any field (nested objects merge). */
export function snap(o: DeepPartial<EngineSnapshot> = {}): EngineSnapshot {
  const base: EngineSnapshot = {
    t: 0, cycle: 1, thetaDeg: -20, rpm: 600,
    pistonDisplacement: 0.005, clearanceHeight: 0.02, rodAngle: 0, intakeLift: 0, exhaustLift: 0,
    phase: 'compression', volume: Math.PI * TEST_R * TEST_R * 0.02, pressure: 12e5,
    temperatureMean: 700, temperatureUnburned: 700, temperatureBurned: 0, massFractionBurned: 0,
    mass: 6e-4, heatReleaseRate: 0, heatLossRate: 0,
    flame: { stage: 'none', radius: 0, center: [0, -0.002, TEST_R - 0.004], area: 0, laminarSpeed: 0, turbulentSpeed: 0, turbulenceIntensity: 1.2 },
    spark: { phase: 'off', primaryCurrent: 0, secondaryVoltage: 0, secondaryCurrent: 0, energyDelivered: 0, breakdownVoltage: 10e3 },
    intakeMassFlow: 0, exhaustMassFlow: 0, intakeManifoldPressure: 1e5, exhaustManifoldPressure: 1.05e5,
    knock: { integral: 0, autoignited: false, oscillation: 0 },
    burnedComposition: { CO2: 0.12, H2O: 0.13, CO: 0.005, O2: 0.005, H2: 0.002, OH: 0.002, H: 0.0005, O: 0.0003, NO: 0.003, N2: 0.73 },
    gasTorque: 0, netTorque: 0,
  };
  const out = { ...base, ...o } as EngineSnapshot;
  out.flame = { ...base.flame, ...(o.flame ?? {}) } as EngineSnapshot['flame'];
  out.spark = { ...base.spark, ...(o.spark ?? {}) } as EngineSnapshot['spark'];
  out.knock = { ...base.knock, ...(o.knock ?? {}) } as EngineSnapshot['knock'];
  out.burnedComposition = { ...base.burnedComposition, ...(o.burnedComposition ?? {}) } as EngineSnapshot['burnedComposition'];
  return out;
}

/** A developed turbulent flame snapshot. */
export function flameSnap(o: DeepPartial<EngineSnapshot> = {}): EngineSnapshot {
  return snap({
    phase: 'combustion', thetaDeg: 5, pressure: 20e5, temperatureUnburned: 750, temperatureBurned: 2400,
    massFractionBurned: 0.1, heatReleaseRate: 1.2e5,
    ...o,
    flame: { stage: 'turbulent', radius: 0.025, area: 2e-3, laminarSpeed: 0.8, turbulentSpeed: 3, turbulenceIntensity: 1.2, ...(o.flame ?? {}) },
  });
}

// ---------------------------------------------------------------------------
// Model T (L-head, four cylinders, trembler-magneto ignition)
// ---------------------------------------------------------------------------

/** MODEL_T with a sharp-cornered valve pocket (exact rectangle cross-sections for flux checks). */
export function sharpLHeadSpec(): EngineSpec {
  const g = MODEL_T.geometry;
  const lh = g.lHead!;
  return { ...MODEL_T, geometry: { ...g, lHead: { ...lh, pocket: { ...lh.pocket, cornerRadius: 0 } } } };
}

const MT_R = MODEL_T.geometry.bore / 2;
const MT_LH = MODEL_T.geometry.lHead!;
/** Model T bore-column depth at TDC (crown below the roof over the bore), m. */
export const MT_DEPTH_TDC = -MT_LH.deckY - MT_LH.crownAboveDeckAtTDC;
/** Model T valve-pocket volume (spec plan area × height), m³. */
const MT_POCKET_V = 46.016e-4 * (MT_LH.pocket.roofY - MT_LH.deckY);

/** A Model T cylinder-1 snapshot at piston displacement x (mid-compression, no flame); override any field. */
export function modelTSnap(o: DeepPartial<EngineSnapshot> = {}, x = 0.005): EngineSnapshot {
  const h = MT_DEPTH_TDC + x;
  return snap({
    rpm: 1000, pistonDisplacement: x, clearanceHeight: h, volume: Math.PI * MT_R * MT_R * h + MT_POCKET_V,
    mass: 4e-4, pressure: 6e5, temperatureMean: 600, temperatureUnburned: 600,
    ...o,
    flame: { center: [...MODEL_T.sparkPlug.gapCenter] as [number, number, number], ...(o.flame ?? {}) },
  });
}

/** The per-cylinder record of a snapshot (as the worker builds cylinders[i]). */
export function cylinderOf(s: EngineSnapshot, index: number, thetaDeg = s.thetaDeg, cycle = s.cycle): CylinderSnapshot {
  return {
    index, thetaDeg, cycle, gasTorque: s.gasTorque,
    pistonDisplacement: s.pistonDisplacement, clearanceHeight: s.clearanceHeight, rodAngle: s.rodAngle,
    intakeLift: s.intakeLift, exhaustLift: s.exhaustLift, phase: s.phase, volume: s.volume, pressure: s.pressure,
    temperatureMean: s.temperatureMean, temperatureUnburned: s.temperatureUnburned, temperatureBurned: s.temperatureBurned,
    massFractionBurned: s.massFractionBurned, mass: s.mass, heatReleaseRate: s.heatReleaseRate, heatLossRate: s.heatLossRate,
    flame: s.flame, spark: s.spark, intakeMassFlow: s.intakeMassFlow, exhaustMassFlow: s.exhaustMassFlow,
    knock: s.knock, burnedComposition: s.burnedComposition,
  };
}

/** An engine snapshot whose top level is cylinder 1 (= cylinders[0]) plus the given per-cylinder states. */
export function engineSnap(cyl: EngineSnapshot[], thetas: number[] = cyl.map((c) => c.thetaDeg)): EngineSnapshot {
  const top = { ...cyl[0] };
  top.cylinders = cyl.map((c, i) => cylinderOf(c, i, thetas[i], c.cycle));
  return top;
}
