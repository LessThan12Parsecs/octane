import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import { fuelLabel, fuelSelectionFromParams, mergeOp, opFromParams, paramsFromOp, patchForParam, type OpParamKey } from './op-binding';

const OP: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 0.8,
  ambientPressure: 101325,
  ambientTemperature: 298.15,
  relativeHumidity: 0.3,
  intakeMixtureTemperature: 325.15,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: 1.05,
  sparkAdvanceDeg: 13,
  dwellTime: 3e-3,
  compressionRatio: 7.2,
  egrFraction: 0.05,
  coolantTemperature: 373.15,
};

describe('operating point ↔ control values', () => {
  it('converts to UI units', () => {
    const p = paramsFromOp(OP);
    expect(p.throttlePct).toBe(80);
    expect(p.dwellMs).toBe(3);
    expect(p.intakeMixtureC).toBe(52);
    expect(p.ambientKPa).toBe(101.325);
    expect(p.coolantC).toBe(100);
    expect(p.humidityPct).toBe(30);
    expect(p.egrPct).toBe(5);
    expect(p.fuel).toBe('PRF');
    expect(p.octaneNumber).toBe(90);
  });

  it('round-trips', () => {
    const back = opFromParams(paramsFromOp(OP));
    for (const k of Object.keys(OP) as (keyof OperatingPoint)[]) {
      if (typeof OP[k] !== 'number') expect(back[k]).toEqual(OP[k]);
      else expect(back[k] as number).toBeCloseTo(OP[k] as number, 9);
    }
  });

  it('each control produces an SI patch for its own field', () => {
    const p = paramsFromOp(OP);
    p.intakeMixtureC = 149;
    expect(patchForParam(p, 'intakeMixtureC')).toEqual({ intakeMixtureTemperature: 422.15 });
    p.dwellMs = 2.5;
    expect(patchForParam(p, 'dwellMs').dwellTime).toBeCloseTo(2.5e-3, 15);
    p.ambientKPa = 95;
    expect(patchForParam(p, 'ambientKPa')).toEqual({ ambientPressure: 95000 });
    p.throttlePct = 25;
    expect(patchForParam(p, 'throttlePct')).toEqual({ throttle: 0.25 });
    const keys: OpParamKey[] = ['speedMode', 'rpm', 'loadTorque', 'sparkAdvanceDeg', 'equivalenceRatio', 'compressionRatio', 'egrPct', 'ambientC', 'humidityPct', 'coolantC'];
    for (const k of keys) expect(Object.keys(patchForParam(p, k))).toHaveLength(1);
  });

  it('fuel selection', () => {
    expect(fuelSelectionFromParams({ fuel: 'PRF', octaneNumber: 120 })).toEqual({ kind: 'PRF', octaneNumber: 100 });
    expect(fuelSelectionFromParams({ fuel: 'C2H5OH', octaneNumber: 90 })).toEqual({ kind: 'pure', species: 'C2H5OH' });
    const p = paramsFromOp({ ...OP, fuel: { kind: 'pure', species: 'CH4' } }, paramsFromOp(OP));
    expect(p.fuel).toBe('CH4');
    expect(p.octaneNumber).toBe(90); // remembered from before
    expect(fuelLabel({ kind: 'PRF', octaneNumber: 92.5 })).toBe('PRF 92.5');
    expect(fuelLabel({ kind: 'pure', species: 'IC8H18' })).toBe('iso-octane');
  });

  it('mergeOp copies the fuel', () => {
    const f = { kind: 'PRF' as const, octaneNumber: 80 };
    const m = mergeOp(OP, { fuel: f, rpm: 900 });
    expect(m.rpm).toBe(900);
    expect(m.fuel).toEqual(f);
    expect(m.fuel).not.toBe(f);
    expect(OP.rpm).toBe(600);
  });
});
