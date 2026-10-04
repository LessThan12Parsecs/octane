import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../physics/core/engine-spec';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { ENGINES } from '../physics/engines/index';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import { controlVisibility, definitionForSpec, engineTitle, gearOptions, loadModelOptions, usesCfrRatingPresets } from './engine-ui';
import { paramsFromOp, type OpParams } from './op-binding';

const CFR = ENGINES['cfr-f1'];
const T = ENGINES['ford-model-t'];
const params = (over: Partial<OpParams>, op = MODEL_T_CRUISE): OpParams => ({ ...paramsFromOp(op), ...over });

describe('controls per engine profile', () => {
  it('CFR: CR slider and dwell, no ignition switch; constant/brake loads in free mode', () => {
    const fixed = controlVisibility(CFR, paramsFromOp(CFR_RON_CONDITIONS));
    expect(fixed).toMatchObject({ compressionRatio: true, dwell: true, ignitionSource: false, octaneNumber: true });
    expect(fixed).toMatchObject({ loadModel: false, loadTorque: false, gear: false, grade: false, brake: false });
    const free = controlVisibility(CFR, params({ speedMode: 'free', loadModel: 'constant' }, CFR_RON_CONDITIONS));
    expect(free).toMatchObject({ loadModel: true, loadTorque: true, gear: false, brake: false });
    // A vehicle load is not offered for the CFR: falls back to the constant-torque controls.
    expect(controlVisibility(CFR, params({ speedMode: 'free', loadModel: 'vehicle' }, CFR_RON_CONDITIONS))).toMatchObject({ gear: false, loadTorque: true });
    expect(loadModelOptions(CFR.ui)).toEqual({ 'Constant torque': 'constant', 'Brake (T ∝ nᵏ)': 'brake' });
  });

  it('Model T: fixed head (no CR), no dwell, MAG/BAT switch; gear and grade for the vehicle load', () => {
    const v = controlVisibility(T, paramsFromOp(MODEL_T_CRUISE));
    expect(v).toMatchObject({ compressionRatio: false, dwell: false, ignitionSource: true });
    expect(v).toMatchObject({ loadModel: true, loadTorque: false, gear: true, grade: true, brake: false });
    expect(controlVisibility(T, params({ gear: 'neutral' }))).toMatchObject({ gear: true, grade: false });
    expect(controlVisibility(T, params({ loadModel: 'brake' }))).toMatchObject({ gear: false, brake: true, loadTorque: true });
    expect(controlVisibility(T, params({ loadModel: 'constant' }))).toMatchObject({ gear: false, brake: false, loadTorque: true });
    expect(controlVisibility(T, params({ speedMode: 'fixed' }))).toMatchObject({ loadModel: false, loadTorque: false, gear: false });
    expect(controlVisibility(T, params({ fuel: 'IC8H18' })).octaneNumber).toBe(false);
  });

  it('gears come from the vehicle, plus neutral', () => {
    const g = gearOptions(MODEL_T);
    expect(Object.values(g).sort()).toEqual(['high', 'low', 'neutral', 'reverse']);
    expect(Object.keys(g)).toContain('Low (2.75:1)');
    expect(gearOptions(CFR_F1)).toEqual({ 'Neutral (declutched)': 'neutral' });
  });

  it('presets: the CFR keeps its rating presets, other engines use the registry list', () => {
    expect(usesCfrRatingPresets(CFR)).toBe(true);
    expect(usesCfrRatingPresets(T)).toBe(false);
    expect(T.presets.length).toBeGreaterThan(0);
    expect(engineTitle(CFR)).toBe('Octane · CFR F-1');
    expect(engineTitle(T)).toBe('Octane · Ford Model T');
  });
});

describe('definitionForSpec', () => {
  it('registry specs map to their entry (modified copies keep the profile with their own spec)', () => {
    expect(definitionForSpec(CFR_F1)).toBe(CFR);
    expect(definitionForSpec(MODEL_T)).toBe(T);
    const mod: EngineSpec = { ...CFR_F1, geometry: { ...CFR_F1.geometry, pinOffset: 0.004 } };
    const d = definitionForSpec(mod);
    expect(d.spec).toBe(mod);
    expect(d.ui).toBe(CFR.ui);
  });

  it('ad-hoc specs get a profile matching what they can do', () => {
    const { id: _id, ...noId } = MODEL_T;
    void _id;
    const d = definitionForSpec(noId as EngineSpec);
    expect(d.ui.showCompressionRatio).toBe(false);
    expect(d.ui.showDwell).toBe(false);
    expect(d.ui.showIgnitionSource).toBe(true);
    expect(d.ui.loadModels).toContain('vehicle');
    expect(d.label).toBe(MODEL_T.name);
    const { id: _c, ...cfrNoId } = CFR_F1;
    void _c;
    const c = definitionForSpec(cfrNoId as EngineSpec);
    expect(c.ui.showCompressionRatio).toBe(true);
    expect(c.ui.showDwell).toBe(true);
    expect(c.presets).toEqual([]);
  });
});
