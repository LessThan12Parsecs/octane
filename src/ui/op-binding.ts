/**
 * Two-way mapping between the SI OperatingPoint and the user-facing control values
 * (°C, kPa, %, ms …). Pure; the lil-gui panel binds to an OpParams object.
 */
import type { FuelSelection, LoadSpec, OperatingPoint } from '../physics/core/operating-point';
import {
  celsiusToKelvin,
  fractionToPct,
  kelvinToCelsius,
  kPaToPa,
  msToS,
  paToKPa,
  pctToFraction,
  sToMs,
} from './units';

export type PureFuel = Extract<FuelSelection, { kind: 'pure' }>['species'];
export type FuelChoice = 'PRF' | PureFuel;
export type LoadModel = LoadSpec['kind'];
export type IgnitionSource = NonNullable<OperatingPoint['ignitionSource']>;

/** Dropdown label → value. */
export const FUEL_OPTIONS: Readonly<Record<string, FuelChoice>> = {
  'PRF blend (iso-octane + n-heptane)': 'PRF',
  'iso-octane (C₈H₁₈)': 'IC8H18',
  'n-heptane (C₇H₁₆)': 'NC7H16',
  'methane (CH₄)': 'CH4',
  'propane (C₃H₈)': 'C3H8',
  'ethanol (C₂H₅OH)': 'C2H5OH',
};

export interface OpParams {
  speedMode: 'fixed' | 'free';
  rpm: number;
  loadTorque: number;
  throttlePct: number;
  sparkAdvanceDeg: number;
  equivalenceRatio: number;
  dwellMs: number;
  compressionRatio: number;
  egrPct: number;
  fuel: FuelChoice;
  /** PRF octane number (vol-% iso-octane), used when fuel = 'PRF'. */
  octaneNumber: number;
  intakeMixtureC: number;
  ambientKPa: number;
  ambientC: number;
  humidityPct: number;
  coolantC: number;
  /** Load model in 'free' speed mode (OperatingPoint.load.kind; absent = 'constant'). */
  loadModel: LoadModel;
  /** 'vehicle' load: gear name (EngineSpec.vehicle.gears key or 'neutral'). */
  gear: string;
  /** 'vehicle' load: road grade, % (rise / run × 100). */
  gradePct: number;
  /** 'brake' load: reference speed at which the brake absorbs loadTorque, rpm. */
  brakeRefRpm: number;
  /** 'brake' load: speed exponent (2 ≈ fan or hydraulic brake). */
  brakeExponent: number;
  /** Trembler-magneto ignition supply ('MAG' / 'BAT' switch). */
  ignitionSource: IgnitionSource;
}

export type OpParamKey = keyof OpParams;

/** Round for display so the controls don't show 324.99999999 °C artefacts. */
const r = (v: number, decimals: number): number => {
  const f = Math.pow(10, decimals);
  return Math.round(v * f) / f;
};

export function paramsFromOp(op: OperatingPoint, prev?: OpParams): OpParams {
  const fuel: FuelChoice = op.fuel.kind === 'PRF' ? 'PRF' : op.fuel.species;
  const octane = op.fuel.kind === 'PRF' ? op.fuel.octaneNumber : (prev?.octaneNumber ?? 90);
  const load: LoadSpec = op.load ?? { kind: 'constant' };
  return {
    speedMode: op.speedMode,
    rpm: op.rpm,
    loadTorque: op.loadTorque,
    throttlePct: r(fractionToPct(op.throttle), 3),
    sparkAdvanceDeg: op.sparkAdvanceDeg,
    equivalenceRatio: op.equivalenceRatio,
    dwellMs: r(sToMs(op.dwellTime), 4),
    compressionRatio: op.compressionRatio,
    egrPct: r(fractionToPct(op.egrFraction), 3),
    fuel,
    octaneNumber: octane,
    intakeMixtureC: r(kelvinToCelsius(op.intakeMixtureTemperature), 3),
    ambientKPa: r(paToKPa(op.ambientPressure), 3),
    ambientC: r(kelvinToCelsius(op.ambientTemperature), 3),
    humidityPct: r(fractionToPct(op.relativeHumidity), 2),
    coolantC: r(kelvinToCelsius(op.coolantTemperature), 3),
    // Sub-settings of the load models not in use are remembered from `prev`.
    loadModel: load.kind,
    gear: load.kind === 'vehicle' ? load.gear : (prev?.gear ?? 'high'),
    gradePct: load.kind === 'vehicle' ? r(fractionToPct(load.grade), 3) : (prev?.gradePct ?? 0),
    brakeRefRpm: load.kind === 'brake' ? load.refRpm : (prev?.brakeRefRpm ?? op.rpm),
    brakeExponent: load.kind === 'brake' ? load.exponent : (prev?.brakeExponent ?? 2),
    ignitionSource: op.ignitionSource ?? 'magneto',
  };
}

/** The LoadSpec selected by the load-model controls. */
export function loadFromParams(p: Pick<OpParams, 'loadModel' | 'gear' | 'gradePct' | 'brakeRefRpm' | 'brakeExponent'>): LoadSpec {
  switch (p.loadModel) {
    case 'vehicle':
      return { kind: 'vehicle', gear: p.gear, grade: pctToFraction(p.gradePct) };
    case 'brake':
      return { kind: 'brake', refRpm: Math.max(1, p.brakeRefRpm), exponent: p.brakeExponent };
    case 'constant':
      return { kind: 'constant' };
  }
}

export function fuelSelectionFromParams(p: Pick<OpParams, 'fuel' | 'octaneNumber'>): FuelSelection {
  if (p.fuel === 'PRF') return { kind: 'PRF', octaneNumber: Math.min(100, Math.max(0, p.octaneNumber)) };
  return { kind: 'pure', species: p.fuel };
}

/** The OperatingPoint patch produced by a change of one control. */
export function patchForParam(p: OpParams, key: OpParamKey): Partial<OperatingPoint> {
  switch (key) {
    case 'speedMode':
      return { speedMode: p.speedMode };
    case 'rpm':
      return { rpm: p.rpm };
    case 'loadTorque':
      return { loadTorque: p.loadTorque };
    case 'throttlePct':
      return { throttle: pctToFraction(p.throttlePct) };
    case 'sparkAdvanceDeg':
      return { sparkAdvanceDeg: p.sparkAdvanceDeg };
    case 'equivalenceRatio':
      return { equivalenceRatio: p.equivalenceRatio };
    case 'dwellMs':
      return { dwellTime: msToS(p.dwellMs) };
    case 'compressionRatio':
      return { compressionRatio: p.compressionRatio };
    case 'egrPct':
      return { egrFraction: pctToFraction(p.egrPct) };
    case 'fuel':
    case 'octaneNumber':
      return { fuel: fuelSelectionFromParams(p) };
    case 'intakeMixtureC':
      return { intakeMixtureTemperature: celsiusToKelvin(p.intakeMixtureC) };
    case 'ambientKPa':
      return { ambientPressure: kPaToPa(p.ambientKPa) };
    case 'ambientC':
      return { ambientTemperature: celsiusToKelvin(p.ambientC) };
    case 'humidityPct':
      return { relativeHumidity: pctToFraction(p.humidityPct) };
    case 'coolantC':
      return { coolantTemperature: celsiusToKelvin(p.coolantC) };
    case 'loadModel':
    case 'gear':
    case 'gradePct':
    case 'brakeRefRpm':
    case 'brakeExponent':
      return { load: loadFromParams(p) };
    case 'ignitionSource':
      return { ignitionSource: p.ignitionSource };
  }
}

/** Full OperatingPoint from the controls (inverse of paramsFromOp). */
export function opFromParams(p: OpParams): OperatingPoint {
  return {
    speedMode: p.speedMode,
    rpm: p.rpm,
    loadTorque: p.loadTorque,
    throttle: pctToFraction(p.throttlePct),
    ambientPressure: kPaToPa(p.ambientKPa),
    ambientTemperature: celsiusToKelvin(p.ambientC),
    relativeHumidity: pctToFraction(p.humidityPct),
    intakeMixtureTemperature: celsiusToKelvin(p.intakeMixtureC),
    fuel: fuelSelectionFromParams(p),
    equivalenceRatio: p.equivalenceRatio,
    sparkAdvanceDeg: p.sparkAdvanceDeg,
    dwellTime: msToS(p.dwellMs),
    compressionRatio: p.compressionRatio,
    egrFraction: pctToFraction(p.egrPct),
    coolantTemperature: celsiusToKelvin(p.coolantC),
    load: loadFromParams(p),
    ignitionSource: p.ignitionSource,
  };
}

/** Merge a patch into an operating point (fuel and load replaced as a whole). */
export function mergeOp(op: OperatingPoint, patch: Partial<OperatingPoint>): OperatingPoint {
  const out: OperatingPoint = { ...op, ...patch, fuel: patch.fuel ? { ...patch.fuel } : op.fuel };
  if (patch.load) out.load = { ...patch.load };
  return out;
}

/** Short fuel description, e.g. "PRF 90" or "iso-octane". */
export function fuelLabel(f: FuelSelection): string {
  if (f.kind === 'PRF') return `PRF ${Number(f.octaneNumber.toFixed(1))}`;
  switch (f.species) {
    case 'IC8H18':
      return 'iso-octane';
    case 'NC7H16':
      return 'n-heptane';
    case 'CH4':
      return 'methane';
    case 'C3H8':
      return 'propane';
    case 'C2H5OH':
      return 'ethanol';
  }
}
