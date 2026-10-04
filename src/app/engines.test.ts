/** Engine selection: URL knobs, start-up resolution, initial operating points, registry sanity. */
import { describe, expect, it } from 'vitest';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { DEFAULT_ENGINE_ID, ENGINE_IDS, ENGINES } from '../physics/engines/index';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import { engineChoices, initialOperatingPoint, isEngineId, resolveEngine } from './engines';
import { clampOperatingPoint, errorSummary, OP_KNOBS, parseUrlOptions, searchWithEngine } from './sync';

describe('URL engine knobs', () => {
  it('reads ?engine, ?cyl and ?sim outside the operating-point overrides', () => {
    const o = parseUrlOptions('?engine=ford-model-t&cyl=3&sim=mock&rpm=900');
    expect(o.engine).toBe('ford-model-t');
    expect(o.focusCylinder).toBe(2);
    expect(o.simulator).toBe('mock');
    expect(o.op).toEqual({ rpm: 900 }); // the engine never leaks into op
  });

  it('ignores malformed values', () => {
    const o = parseUrlOptions('?engine=../../x&cyl=0&sim=fast');
    expect(o.engine).toBeUndefined();
    expect(o.focusCylinder).toBeUndefined();
    expect(o.simulator).toBeUndefined();
    expect(parseUrlOptions('?cyl=1.5').focusCylinder).toBeUndefined();
    expect(parseUrlOptions('?engine=').engine).toBeUndefined();
  });

  it('writes the engine back, dropping the previous engine’s operating-point knobs and focus', () => {
    expect(searchWithEngine('?cr=7.2&on=95&rpm=900&spark=20&phi=1.0&ts=0.01&view=engine&cyl=2', 'ford-model-t')).toBe(
      '?ts=0.01&view=engine&engine=ford-model-t',
    );
    expect(searchWithEngine('', 'cfr-f1')).toBe('?engine=cfr-f1');
    expect(searchWithEngine('?engine=cfr-f1&paused=1', 'ford-model-t')).toBe('?engine=ford-model-t&paused=1');
    expect(searchWithEngine('?rpm=900', 'cfr-f1', false)).toBe('?rpm=900&engine=cfr-f1');
    // every knob parseUrlOptions maps into op is dropped
    const all = OP_KNOBS.map((k) => `${k}=1`).join('&');
    expect(parseUrlOptions(searchWithEngine(`?${all}`, 'cfr-f1')).op).toEqual({});
  });
});

describe('errorSummary', () => {
  it('drops stack frames and the repeated message, and bounds the length', () => {
    const msg =
      'simulator failed: FlameGeometry: spark centre must lie inside the bore (d < R)\n' +
      'Error: FlameGeometry: spark centre must lie inside the bore (d < R)\n' +
      '    at new FlameGeometry (http://x/flame-geometry.ts:323:23)\n' +
      '    at new CycleModel (http://x/cycle-model.ts:642:20)';
    expect(errorSummary(msg)).toBe('simulator failed: FlameGeometry: spark centre must lie inside the bore (d < R)');
    expect(errorSummary('a\nb\nc\nd')).toBe('a\nb\nc');
    expect(errorSummary('x'.repeat(500)).length).toBe(360);
    expect(errorSummary('')).toBe('');
  });
});

describe('resolveEngine', () => {
  it('takes the first registry id among explicit option, URL and remembered choice', () => {
    expect(resolveEngine('ford-model-t', 'cfr-f1').id).toBe('ford-model-t');
    expect(resolveEngine(undefined, 'ford-model-t', 'cfr-f1').id).toBe('ford-model-t');
    expect(resolveEngine(undefined, 'nope', 'ford-model-t').id).toBe('ford-model-t');
    expect(resolveEngine(undefined, null, null).id).toBe(DEFAULT_ENGINE_ID);
    expect(resolveEngine('toString', '__proto__').id).toBe(DEFAULT_ENGINE_ID); // no prototype lookups
  });

  it('isEngineId only accepts registry ids', () => {
    for (const id of ENGINE_IDS) expect(isEngineId(id)).toBe(true);
    expect(isEngineId('constructor')).toBe(false);
    expect(isEngineId(undefined)).toBe(false);
  });

  it('offers every registry engine in the picker', () => {
    expect(engineChoices().map((c) => c.id)).toEqual(ENGINE_IDS);
    for (const c of engineChoices()) expect(c.label).toBe(ENGINES[c.id as keyof typeof ENGINES].label);
  });
});

describe('initialOperatingPoint', () => {
  it('CFR: identical to the pre-registry app start (RON conditions + URL knobs, CR clamped)', () => {
    const url = parseUrlOptions('?cr=40&on=95&rpm=900').op;
    expect(initialOperatingPoint(ENGINES['cfr-f1'], url)).toEqual(clampOperatingPoint(CFR_F1, { ...CFR_RON_CONDITIONS, ...url }));
    expect(initialOperatingPoint(ENGINES['cfr-f1'])).toEqual(CFR_RON_CONDITIONS);
  });

  it('Model T: starts from its default; a URL CR cannot move the fixed head', () => {
    const op = initialOperatingPoint(ENGINES['ford-model-t'], parseUrlOptions('?cr=7&spark=10').op);
    expect(op.compressionRatio).toBe(MODEL_T.geometry.compressionRatio);
    expect(op.sparkAdvanceDeg).toBe(10);
    expect(op.load).toEqual(MODEL_T_CRUISE.load);
    expect(op.ignitionSource).toBe('magneto');
  });

  it('an explicit base operating point replaces the default', () => {
    const base = { ...CFR_RON_CONDITIONS, rpm: 777 };
    expect(initialOperatingPoint(ENGINES['cfr-f1'], {}, base).rpm).toBe(777);
  });
});

describe('engine registry (as the app uses it)', () => {
  for (const id of ENGINE_IDS) {
    const def = ENGINES[id];
    it(`${id}: default operating point is inside its own limits and the spec survives structured cloning`, () => {
      expect(clampOperatingPoint(def.spec, def.defaultOperatingPoint)).toEqual(def.defaultOperatingPoint);
      expect(def.defaultOperatingPoint.rpm).toBeGreaterThanOrEqual(def.ui.rpmRange[0]);
      expect(def.defaultOperatingPoint.rpm).toBeLessThanOrEqual(def.ui.rpmRange[1]);
      expect(def.defaultOperatingPoint.sparkAdvanceDeg).toBeGreaterThanOrEqual(def.ui.sparkRange[0]);
      expect(def.defaultOperatingPoint.sparkAdvanceDeg).toBeLessThanOrEqual(def.ui.sparkRange[1]);
      expect(structuredClone(def.spec)).toEqual(def.spec);
      expect(def.spec.id).toBe(id);
      for (const p of def.presets) expect(structuredClone(p.op)).toEqual(p.op);
    });
  }
});
