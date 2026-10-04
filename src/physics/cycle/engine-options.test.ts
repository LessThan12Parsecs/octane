/**
 * Engine-keyed option defaults and operating-point sanitising (Model T integration): neutral physics
 * defaults ← per-engine defaults (EngineSpec.id) ← caller options; the CFR set is unchanged; the Model T
 * gets its own friction, valve lash, crankcase pressure, knock pickup/band (the virtual plug transducer and
 * the L-head band) and MODEL_T_CALIBRATION (calibration phase: fitted to Upton's MBT and suction data and
 * Ford's WOT table, plus the Model T-only length-scale and friction factors); sanitizeOperatingPoint falls
 * back to the engine's registry default operating point and keeps the load model and the ignition source.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../core/engine-spec';
import { CFR_CRANKCASE_GAUGE_PRESSURE, CFR_F1, CFR_FRICTION, CFR_KNOCK_PICKUP, CFR_RON_CONDITIONS, CFR_VALVE_LASH } from '../engines/cfr';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_FRICTION, MODEL_T_VALVE_LASH } from '../engines/model-t';
import { L_HEAD_MAPO_BAND, virtualKnockSensor } from '../chemistry/knock';
import { CFR_CALIBRATION, MODEL_T_CALIBRATION } from './calibration';
import { CycleModel } from './cycle-model';
import {
  DEFAULT_CYCLE_MODEL_OPTIONS,
  defaultOperatingPointOf,
  ENGINE_CYCLE_OPTION_DEFAULTS,
  MODEL_T_CALIBRATED_FRICTION,
  NEUTRAL_CYCLE_MODEL_OPTIONS,
  resolveCycleOptions,
  sanitizeOperatingPoint,
  VENTURI_DISCHARGE_COEFFICIENT,
} from './options';

const CAL_KEYS = Object.keys(CFR_CALIBRATION) as (keyof typeof CFR_CALIBRATION)[];

describe('engine-keyed cycle-model option defaults', () => {
  it('the CFR set is the calibrated CFR F-1 set (with or without the spec)', () => {
    expect(resolveCycleOptions({}, CFR_F1)).toEqual(resolveCycleOptions({}));
    expect(resolveCycleOptions({}, CFR_F1)).toEqual({ ...DEFAULT_CYCLE_MODEL_OPTIONS, wiebe: { ...DEFAULT_CYCLE_MODEL_OPTIONS.wiebe } });
    const d = DEFAULT_CYCLE_MODEL_OPTIONS;
    for (const k of CAL_KEYS) expect(d[k], k).toBe(CFR_CALIBRATION[k].value);
    expect(d.friction).toBe(CFR_FRICTION);
    expect(d.valveLash).toBe(CFR_VALVE_LASH);
    expect(d.crankcaseGaugePressure).toBe(CFR_CRANKCASE_GAUGE_PRESSURE);
    expect(d.knockSensor).toEqual([CFR_KNOCK_PICKUP.position[0], CFR_KNOCK_PICKUP.position[2]]);
    expect(d.mapoBand).toEqual([4000, 18000]);
    // a spec without an id is a CFR derivative (pre-registry behaviour)
    expect(resolveCycleOptions({}, { ...CFR_F1, id: undefined })).toEqual(resolveCycleOptions({}));
  });

  it('the Model T gets its own hardware and its calibration', () => {
    const o = resolveCycleOptions({}, MODEL_T);
    // (calibration phase) the PNH inputs of MODEL_T_FRICTION with the two calibrated friction factors
    expect(o.friction).toBe(MODEL_T_CALIBRATED_FRICTION);
    expect(o.friction).toEqual({
      ...MODEL_T_FRICTION,
      ringTensionFactor: MODEL_T_CALIBRATION.ringTensionFactor.value,
      auxiliaryFactor: MODEL_T_CALIBRATION.auxiliaryFactor.value,
    });
    expect(o.turbulenceLengthScaleFactor).toBe(MODEL_T_CALIBRATION.turbulenceLengthScaleFactor.value);
    expect(o.friction!.cylinders).toBe(4);
    expect(o.valveLash).toBe(MODEL_T_VALVE_LASH);
    expect(o.crankcaseGaugePressure).toBe(0);
    // (Model T integration) the L-head acoustics: the virtual plug-mounted pickup over the valve pocket and
    // the 2–18 kHz band that keeps the bore-to-pocket fundamental
    expect(o.mapoBand).toEqual(L_HEAD_MAPO_BAND);
    for (const k of CAL_KEYS) expect(o[k], k).toBe(MODEL_T_CALIBRATION[k].value);
    expect(o.knockSensor).toEqual(virtualKnockSensor(MODEL_T));
    expect(Math.hypot(o.knockSensor[0], o.knockSensor[1])).toBeGreaterThan(MODEL_T.geometry.bore / 2);
    // caller options still win
    expect(resolveCycleOptions({ burnRateMultiplier: 4, friction: null }, MODEL_T)).toMatchObject({ burnRateMultiplier: 4, friction: null, valveLash: MODEL_T_VALVE_LASH });
    expect(ENGINE_CYCLE_OPTION_DEFAULTS['ford-model-t'].friction).toBe(MODEL_T_CALIBRATED_FRICTION);
    // the CFR keeps the neutral length scale
    expect(DEFAULT_CYCLE_MODEL_OPTIONS.turbulenceLengthScaleFactor).toBe(1);
  });

  it('MODEL_T_CALIBRATION (calibration phase): every entry inside its literature range with source and evidence; the universal closures are the CFR values', () => {
    const keys = Object.keys(MODEL_T_CALIBRATION) as (keyof typeof MODEL_T_CALIBRATION)[];
    expect(keys).toEqual(expect.arrayContaining([...CAL_KEYS, 'turbulenceLengthScaleFactor', 'ringTensionFactor', 'auxiliaryFactor']));
    for (const k of keys) {
      const p = MODEL_T_CALIBRATION[k];
      expect(p.value, k).toBeGreaterThanOrEqual(p.range[0]);
      expect(p.value, k).toBeLessThanOrEqual(p.range[1]);
      expect(p.source.length, k).toBeGreaterThan(20);
      expect(p.evidence.length, k).toBeGreaterThan(20);
    }
    // fitted to the Model T data (Upton MBT and suction, Ford's WOT table): no longer the transferred CFR values
    for (const k of ['burnRateMultiplier', 'taylorScaleMultiplier', 'venturiDischargeCoefficient'] as const) {
      expect(MODEL_T_CALIBRATION[k].value, k).not.toBe(CFR_CALIBRATION[k].value);
      expect(MODEL_T_CALIBRATION[k].evidence.startsWith('UNVERIFIED'), k).toBe(false);
    }
    // universal closures (not engine-fitted) stay the CFR values, marked as transferred
    for (const k of ['kernelHandoffMultiple', 'marksteinMultiplier'] as const) {
      expect(MODEL_T_CALIBRATION[k].value, k).toBe(CFR_CALIBRATION[k].value);
      expect(MODEL_T_CALIBRATION[k].evidence).toMatch(/Would be set by: \S/);
    }
  });

  it('an unregistered engine gets the neutral (uncalibrated) set with a pickup inside its bore', () => {
    const spec: EngineSpec = { ...CFR_F1, id: 'some-test-engine' };
    const o = resolveCycleOptions({}, spec);
    expect(o).toMatchObject({
      burnRateMultiplier: 1,
      taylorScaleMultiplier: 1,
      woschniMultiplier: 1,
      intakePortHeatTransferMultiplier: 0,
      venturiDischargeCoefficient: VENTURI_DISCHARGE_COEFFICIENT,
      friction: null,
      valveLash: 0,
      crankcaseGaugePressure: 0,
      mapoBand: null,
    });
    expect(o.marksteinMultiplier).toBe(NEUTRAL_CYCLE_MODEL_OPTIONS.marksteinMultiplier);
    expect(Math.hypot(o.knockSensor[0], o.knockSensor[1])).toBeLessThan(spec.geometry.bore / 2);
  });
});

describe('sanitizeOperatingPoint with the engine registry', () => {
  it('repairs non-finite fields from the engine default operating point', () => {
    expect(defaultOperatingPointOf(MODEL_T)).toBe(MODEL_T_CRUISE);
    expect(defaultOperatingPointOf(CFR_F1)).toBe(CFR_RON_CONDITIONS);
    const o = sanitizeOperatingPoint(MODEL_T, { ...MODEL_T_CRUISE, rpm: Number.NaN, sparkAdvanceDeg: undefined as unknown as number });
    expect(o.rpm).toBe(MODEL_T_CRUISE.rpm);
    expect(o.sparkAdvanceDeg).toBe(MODEL_T_CRUISE.sparkAdvanceDeg);
    // fixed compression ratio: forced to the spec's
    expect(sanitizeOperatingPoint(MODEL_T, { ...MODEL_T_CRUISE, compressionRatio: 8 }).compressionRatio).toBe(MODEL_T.geometry.compressionRatio);
  });

  it('keeps (copies) the load model and the ignition source, repairing invalid ones', () => {
    const op = { ...MODEL_T_CRUISE, ignitionSource: 'battery' as const, load: { kind: 'vehicle' as const, gear: 'low', grade: 0.05 } };
    const o = sanitizeOperatingPoint(MODEL_T, op);
    expect(o.load).toEqual(op.load);
    expect(o.load).not.toBe(op.load);
    expect(o.ignitionSource).toBe('battery');
    const bad = sanitizeOperatingPoint(MODEL_T, { ...op, ignitionSource: 'MAG' as unknown as 'magneto', load: { kind: 'vehicle', gear: 'overdrive', grade: 0 } });
    expect(bad.ignitionSource).toBe('magneto');
    expect(bad.load).toEqual(MODEL_T_CRUISE.load);
    expect(sanitizeOperatingPoint(MODEL_T, { ...op, load: { kind: 'brake', refRpm: 1000, exponent: 9 } }).load).toEqual({ kind: 'brake', refRpm: 1000, exponent: 4 });
    expect(sanitizeOperatingPoint(MODEL_T, { ...op, load: { kind: 'vehicle', gear: 'high', grade: Number.NaN } }).load).toEqual({ kind: 'vehicle', gear: 'high', grade: 0 });
    // the CFR has neither: no keys are added
    const c = sanitizeOperatingPoint(CFR_F1, CFR_RON_CONDITIONS);
    expect('load' in c).toBe(false);
    expect('ignitionSource' in c).toBe(false);
  });
});

describe('wall temperatures without a reference coolant temperature', () => {
  it('take the spec values at the operating point coolant temperature (no CFR fallback)', () => {
    const walls = { ...CFR_F1.walls, thermalResistance: undefined, referenceCoolantTemperature: undefined };
    const spec: EngineSpec = { ...CFR_F1, walls };
    const m = new CycleModel(spec, { ...CFR_RON_CONDITIONS, coolantTemperature: 350 }, { combustionModel: 'none', warmupCycles: 0 });
    expect(m.walls.headTemperature).toBe(walls.headTemperature);
    expect(m.walls.linerTemperature).toBe(walls.linerTemperature);
  });
});
