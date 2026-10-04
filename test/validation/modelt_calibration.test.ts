/**
 * Fast sanity checks of the Model T calibration (calibration phase, 2026-10; the full comparison with the period
 * data is test/validation/measured_modelt_data.test.ts, gated by OCTANE_VALIDATION):
 *  - the CFR calibration set is untouched (its values pinned literally: the CFR must stay bit-identical);
 *  - every Model T calibrated parameter lies inside its documented range and the engine defaults use it;
 *  - one wide-open point (1500 rpm, lever fully advanced — the torque maximum there) is within ±10 % of Ford's
 *    table and the friction is in the calibrated band (≈ 5 s).
 */
import { describe, expect, it } from 'vitest';
import { CFR_CALIBRATION, MODEL_T_CALIBRATION } from '../../src/physics/cycle/calibration';
import { CycleModel } from '../../src/physics/cycle/index';
import { MODEL_T_CALIBRATED_FRICTION, resolveCycleOptions } from '../../src/physics/cycle/options';
import { CFR_FRICTION } from '../../src/physics/engines/cfr';
import { MODEL_T, MODEL_T_FORD_DYNO, MODEL_T_FORD_WOT_TABLE, MODEL_T_FRICTION } from '../../src/physics/engines/model-t';

const LBFT = 1.3558179483314004;

describe('Model T calibration: parameters', () => {
  it('the CFR F-1 calibration is unchanged by the Model T calibration', () => {
    expect(Object.fromEntries(Object.entries(CFR_CALIBRATION).map(([k, p]) => [k, p.value]))).toEqual({
      burnRateMultiplier: 5.3,
      taylorScaleMultiplier: 3.0,
      kernelHandoffMultiple: 2.0,
      marksteinMultiplier: 0.5,
      woschniMultiplier: 1.2,
      intakePortHeatTransferMultiplier: 6,
      venturiDischargeCoefficient: 0.6,
      knockStratificationDT: 60,
    });
    expect(resolveCycleOptions({}).turbulenceLengthScaleFactor).toBe(1);
    expect(resolveCycleOptions({}).friction).toBe(CFR_FRICTION);
  });

  it('every Model T parameter is inside its range; the defaults of the Model T are those values', () => {
    const o = resolveCycleOptions({}, MODEL_T);
    for (const [k, p] of Object.entries(MODEL_T_CALIBRATION)) {
      expect(p.value, k).toBeGreaterThanOrEqual(p.range[0]);
      expect(p.value, k).toBeLessThanOrEqual(p.range[1]);
      if (k in o && typeof (o as unknown as Record<string, unknown>)[k] === 'number') expect((o as unknown as Record<string, number>)[k], k).toBe(p.value);
    }
    expect(o.friction).toBe(MODEL_T_CALIBRATED_FRICTION);
    expect(o.friction!.ringTensionFactor).toBe(MODEL_T_CALIBRATION.ringTensionFactor.value);
    expect(o.friction!.auxiliaryFactor).toBe(MODEL_T_CALIBRATION.auxiliaryFactor.value);
    // the PNH inputs themselves (geometry, oil, valvetrain) are the spec's
    expect({ ...o.friction, ringTensionFactor: MODEL_T_FRICTION.ringTensionFactor, auxiliaryFactor: MODEL_T_FRICTION.auxiliaryFactor }).toEqual(MODEL_T_FRICTION);
  });
});

describe('Model T calibration: one wide-open point', () => {
  it("1500 rpm, lever 64.5: brake torque within ±10 % of Ford's 70 lb-ft, FMEP 0.9–1.3 bar, no knock, no misfire", () => {
    const m = new CycleModel(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm: 1500, sparkAdvanceDeg: 64.5 }, { warmupCycles: 2 });
    m.runCycles(1);
    const s = m.runCycles(1);
    const e = s.find((x) => x.engine)!.engine!;
    const ford = MODEL_T_FORD_WOT_TABLE.find((r) => r[0] === 1500)![1];
    const tb = e.brakeTorque / LBFT;
    console.log(`[Model T] 1500 rpm WOT: ${tb.toFixed(1)} lb-ft (Ford ${ford}), FMEP ${(e.fmep / 1e5).toFixed(2)} bar, η_v ${e.volumetricEfficiency.toFixed(3)}`);
    expect(Math.abs(tb / ford - 1)).toBeLessThan(0.1);
    expect(e.fmep / 1e5).toBeGreaterThan(0.9);
    expect(e.fmep / 1e5).toBeLessThan(1.3);
    for (const x of s) {
      expect(x.misfire).toBe(false);
      expect(Number.isFinite(x.knockOnsetDeg)).toBe(false);
    }
  });
});
