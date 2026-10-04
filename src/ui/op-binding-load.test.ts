/** Load-model and ignition-source controls (free-speed engines, the Model T's MAG/BAT switch). */
import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import { CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { ENGINES } from '../physics/engines/index';
import { MODEL_T_CRUISE } from '../physics/engines/model-t';
import { loadFromParams, mergeOp, opFromParams, paramsFromOp, patchForParam } from './op-binding';

describe('load model ↔ controls', () => {
  it('reads the vehicle load and the ignition source', () => {
    const p = paramsFromOp(MODEL_T_CRUISE);
    expect(p.loadModel).toBe('vehicle');
    expect(p.gear).toBe('high');
    expect(p.gradePct).toBe(0);
    expect(p.ignitionSource).toBe('magneto');
    // absent load = constant; absent source = magneto
    const c = paramsFromOp(CFR_RON_CONDITIONS);
    expect(c.loadModel).toBe('constant');
    expect(c.ignitionSource).toBe('magneto');
  });

  it('every load control produces one `load` patch with the whole LoadSpec', () => {
    const p = paramsFromOp(MODEL_T_CRUISE);
    p.gradePct = 10;
    expect(patchForParam(p, 'gradePct')).toEqual({ load: { kind: 'vehicle', gear: 'high', grade: 0.1 } });
    p.gear = 'low';
    expect(patchForParam(p, 'gear')).toEqual({ load: { kind: 'vehicle', gear: 'low', grade: 0.1 } });
    p.loadModel = 'brake';
    p.brakeRefRpm = 1200;
    p.brakeExponent = 2;
    expect(patchForParam(p, 'loadModel')).toEqual({ load: { kind: 'brake', refRpm: 1200, exponent: 2 } });
    expect(patchForParam(p, 'brakeExponent')).toEqual(patchForParam(p, 'brakeRefRpm'));
    p.loadModel = 'constant';
    expect(patchForParam(p, 'loadModel')).toEqual({ load: { kind: 'constant' } });
    p.ignitionSource = 'battery';
    expect(patchForParam(p, 'ignitionSource')).toEqual({ ignitionSource: 'battery' });
  });

  it('remembers the settings of the load models not in use', () => {
    const p = paramsFromOp(MODEL_T_CRUISE);
    p.gradePct = 5;
    const braked = paramsFromOp({ ...MODEL_T_CRUISE, load: { kind: 'brake', refRpm: 1500, exponent: 1.5 } }, p);
    expect(braked.loadModel).toBe('brake');
    expect(braked.gradePct).toBe(5);
    expect(braked.gear).toBe('high');
    const back = paramsFromOp(MODEL_T_CRUISE, braked);
    expect(back.brakeRefRpm).toBe(1500);
    expect(back.brakeExponent).toBe(1.5);
  });

  it('brake reference speed never reaches zero', () => {
    expect(loadFromParams({ loadModel: 'brake', gear: 'high', gradePct: 0, brakeRefRpm: 0, brakeExponent: 2 })).toEqual({ kind: 'brake', refRpm: 1, exponent: 2 });
  });

  it('round-trips every Model T preset through the controls', () => {
    for (const pr of ENGINES['ford-model-t'].presets) {
      const op: OperatingPoint = mergeOp(MODEL_T_CRUISE, pr.op);
      const back = opFromParams(paramsFromOp(op));
      for (const k of Object.keys(op) as (keyof OperatingPoint)[]) {
        const want = op[k] ?? (k === 'load' ? { kind: 'constant' } : undefined);
        if (typeof want === 'number') expect(back[k] as number, `${pr.id}.${k}`).toBeCloseTo(want, 9);
        else expect(back[k], `${pr.id}.${k}`).toEqual(want);
      }
    }
  });

  it('mergeOp copies the load spec', () => {
    const load = { kind: 'vehicle' as const, gear: 'low', grade: 0.05 };
    const m = mergeOp(MODEL_T_CRUISE, { load });
    expect(m.load).toEqual(load);
    expect(m.load).not.toBe(load);
    expect(mergeOp(CFR_RON_CONDITIONS, {})).not.toHaveProperty('load');
    expect(mergeOp(MODEL_T_CRUISE, { load: undefined }).load).toBeUndefined();
  });
});
