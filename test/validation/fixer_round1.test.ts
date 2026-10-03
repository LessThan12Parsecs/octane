/**
 * Validation round 1 — fixer pass: regression tests for the fixes and additions that the validators'
 * test files do not already cover (each fails on the round-1 code). Runtime ≈ 15–30 s.
 *
 *  - closure robustness: no failed closure solve at 4500/6000 rpm (round 1: the coil's switch-on
 *    transient fired a "make" spark ≈ 100° early, a 1e-11 kg burned zone received tens of mJ, T_b was
 *    pinned at 5900 K and 4436 of 4477 solves failed at 6000 rpm);
 *  - spark-energy placement: the kernel-phase burned zone stays at the kernel temperature (round 1:
 *    3279 K on the first kernel step vs T_ad ≈ 2420 K);
 *  - sanitiser: φ = 0 (air only, ASTM drained carburettor) is accepted; RH is limited so x_H2O ≤ 0.2;
 *  - calibration table consistency (defaults = CFR_CALIBRATION, each value inside its range);
 *  - new options: intake-port heat transfer, Keck-1982 closure, two-stage Livengood–Wu;
 *  - KnockOscillator: exact exponential source and between-sample peak.
 */
import { describe, expect, it } from 'vitest';
import { SP } from '../../src/physics/core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { humidAir } from '../../src/physics/thermo/fuels';
import { KnockOscillator } from '../../src/physics/chemistry/knock';
import {
  CFR_CALIBRATION,
  CycleModel,
  DEFAULT_CYCLE_MODEL_OPTIONS,
  MODE_SINGLE,
  MODE_TWO,
  sanitizeOperatingPoint,
} from '../../src/physics/cycle/index';

describe('closure robustness and spark-energy placement', () => {
  for (const rpm of [4500, 6000]) {
    it(`${rpm} rpm: no failed closure solve and no burned zone before the commanded spark`, () => {
      const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, rpm }, { warmupCycles: 1 });
      const f0 = m.closureFailures;
      let early = false;
      while (m.cycle < 2) {
        m.stepUntil(Infinity, Math.min(360, m.theta + 1));
        if (m.mode === MODE_TWO && m.theta < m.ignCmd.sparkDeg - 1e-6 && m.theta > -150) early = true;
      }
      expect(m.closureFailures - f0).toBe(0);
      expect(early).toBe(false);
    });
  }

  it('RON: the burned zone never exceeds the kernel temperature during the kernel phase', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1 });
    let maxExcess = -Infinity;
    const c0 = m.cycle;
    while (m.cycle === c0) {
      m.stepUntil(Infinity, Math.min(360, m.theta + 0.05));
      const k = m.ign!.state.kernel;
      if (m.mode === MODE_TWO && m.dwellSeen && k.stage === 'kernel' && k.temperature > 0) maxExcess = Math.max(maxExcess, m.Tb - k.temperature);
    }
    expect(Number.isFinite(maxExcess)).toBe(true);
    expect(maxExcess).toBeLessThan(25); // (the burned zone is then compressed with the charge)
  });
});

describe('operating-point sanitiser', () => {
  it('accepts φ = 0 (air only): no fuel trapped, knock disarmed, a finite motored-like cycle', () => {
    const op = sanitizeOperatingPoint(CFR_F1, { ...CFR_RON_CONDITIONS, equivalenceRatio: 0 });
    expect(op.equivalenceRatio).toBe(0);
    const m = new CycleModel(CFR_F1, op, { warmupCycles: 1 });
    const s = m.runCycles(1)[0];
    expect(m.fuelMassIvc).toBe(0);
    expect(m.knockAvailable).toBe(false);
    expect(Number.isFinite(s.peakPressure)).toBe(true);
  });

  it('limits the relative humidity so the intake air holds at most 20 % water vapour', () => {
    const op = sanitizeOperatingPoint(CFR_F1, { ...CFR_RON_CONDITIONS, ambientTemperature: 400, relativeHumidity: 1 });
    const x = humidAir(op.ambientTemperature, op.ambientPressure, op.relativeHumidity)[SP.H2O];
    expect(x).toBeLessThanOrEqual(0.2 + 1e-12);
    expect(x).toBeGreaterThan(0.19);
  });
});

describe('calibration (calibration.ts)', () => {
  it('the defaults are the CFR calibration and every value lies in its literature range', () => {
    const d = DEFAULT_CYCLE_MODEL_OPTIONS as unknown as Record<string, number>;
    for (const [k, p] of Object.entries(CFR_CALIBRATION)) {
      expect(d[k], k).toBe(p.value);
      expect(p.value, k).toBeGreaterThanOrEqual(p.range[0]);
      expect(p.value, k).toBeLessThanOrEqual(p.range[1]);
      expect(p.source.length, k).toBeGreaterThan(20);
      expect(p.evidence.length, k).toBeGreaterThan(20);
    }
  });
});

describe('new sub-model options', () => {
  it('intake-port heat transfer warms the charge (T_IVC up, trapped mass down); 0 = adiabatic', () => {
    const a = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 2, intakePortHeatTransferMultiplier: 0 });
    const b = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 2 });
    const sa = a.runCycles(1)[0];
    const sb = b.runCycles(1)[0];
    expect(b.TIvc - a.TIvc).toBeGreaterThan(10);
    expect(sb.trappedMass).toBeLessThan(sa.trappedMass * 0.98);
  });

  it("turbulentFlameClosure 'keck1982' burns with Keck's empirical u_T and ℓ_T", () => {
    const k = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 2, turbulentFlameClosure: 'keck1982', burnRateMultiplier: 1, taylorScaleMultiplier: 1 });
    const s = k.runCycles(1)[0];
    expect(k.keckInletSpeed).toBeGreaterThan(5);
    expect(Number.isFinite(s.ca50)).toBe(true);
    expect(s.misfire).toBe(false);
  });

  it("knockIntegral 'two-stage' (LLNL tables) passes through the first stage and stays finite", () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 60 }, compressionRatio: 7 }, {
      warmupCycles: 1,
      ignitionDelayModel: 'llnl-gasoline-2011',
      knockIntegral: 'two-stage',
    });
    let sawStage2 = false;
    const c0 = m.cycle;
    while (m.cycle === c0) {
      m.stepUntil(Infinity, Math.min(360, m.theta + 1));
      const p = m.knockProgress;
      expect(p >= 0 && p <= 1).toBe(true);
      if (m.lwStage === 2) sawStage2 = true;
    }
    expect(sawStage2).toBe(true);
    const s = m.drainSummaries().pop()!;
    expect(Number.isFinite(s.imepNet)).toBe(true);
  });
});

describe('KnockOscillator (chemistry/knock.ts)', () => {
  it('an exponentially decaying source is propagated exactly: one step equals many sub-steps', () => {
    const run = (n: number) => {
      const o = new KnockOscillator(0.08255, { quadrature: [16, 32] });
      o.setEndGasRegion(0.02, 0, 0.03);
      const T = 200e-6;
      const tau = 60e-6;
      const q0 = 2e6;
      const dt = T / n;
      for (let i = 0; i < n; i++) o.step(dt, 950, 1.3, 1e-4, q0 * Math.exp((-i * dt) / tau), tau);
      o.step(1e-4, 950, 1.3, 1e-4, 0);
      return { p: o.sensorPressure(), peak: o.peakSensorPressure };
    };
    const a = run(1);
    const b = run(50);
    expect(Math.abs(a.p - b.p)).toBeLessThan(1e-9 * Math.abs(b.p) + 1e-9);
  });

  it('peakSensorPressure catches the maximum between coarse samples (cubic Hermite)', () => {
    const coarse = new KnockOscillator(0.08255, { quadrature: [16, 32] });
    const fine = new KnockOscillator(0.08255, { quadrature: [16, 32] });
    for (const o of [coarse, fine]) o.setEndGasRegion(0.02, 0, 0.03);
    coarse.step(2e-6, 950, 1.3, 1e-4, 5e7);
    fine.step(2e-6, 950, 1.3, 1e-4, 5e7);
    let dense = 0;
    for (let i = 0; i < 1000; i++) {
      fine.step(1e-6, 950, 1.3, 1e-4, 0);
      dense = Math.max(dense, Math.abs(fine.sensorPressure()));
    }
    for (let i = 0; i < 500; i++) coarse.step(2e-6, 950, 1.3, 1e-4, 0);
    expect(Math.abs(coarse.peakSensorPressure / dense - 1)).toBeLessThan(2e-3);
  });

  it('a non-finite source is ignored (never poisons the modal amplitudes)', () => {
    const o = new KnockOscillator(0.08255, { quadrature: [16, 32] });
    o.setEndGasRegion(0.02, 0, 0.03);
    o.step(2e-6, 950, 1.3, 1e-4, 1e6);
    o.step(2e-6, 950, 1.3, 1e-4, NaN);
    expect(Number.isFinite(o.sensorPressure())).toBe(true);
  });
});

void MODE_SINGLE;
