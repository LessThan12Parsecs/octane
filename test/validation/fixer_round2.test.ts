/**
 * Validation round 2 — fixer pass: regression tests for fixes and model additions that the round-2
 * validators' test files (numerics_round2, code-review-r2, knock_octane_round2, measured_cfr_round2)
 * do not already cover. Each fails on the round-1 code (the reason is given per test). Runtime ≈ 20–40 s.
 */
import { describe, expect, it } from 'vitest';
import { NS } from '../../src/physics/core/species';
import { R_UNIVERSAL } from '../../src/physics/core/constants';
import {
  CFR_ANL2020S_CLEARANCE_OFFSET,
  CFR_F1,
  CFR_RON_CONDITIONS,
  cfrCompressionRatioFromAnl2020s,
} from '../../src/physics/engines/cfr';
import { CycleModel, MODE_OPEN, MODE_TWO } from '../../src/physics/cycle/index';
import { I_CRB, I_CRU } from '../../src/physics/cycle/cycle-model';

describe('flame radius is continuous through the kernel hand-off (code review r2 #14)', () => {
  // round 1 reported the burned-gas radius during the kernel stage and the entrained front after the
  // hand-off: +0.48 mm (6.7 %) in one step at RON
  it('RON: the radius increment of the hand-off step is continuous with the kernel-stage increments', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1 });
    m.recordTrace();
    m.runCycles(1);
    const tr = m.trace!;
    let ratio = NaN;
    for (let i = 2; i < tr.theta.length; i++) {
      // the hand-off: first sample with a front area after kernel-stage samples (equal 0.05° steps)
      if (tr.frontArea[i - 1] === 0 && tr.frontArea[i] > 0 && tr.flameRadius[i - 2] > 0) {
        const dPrev = tr.flameRadius[i - 1] - tr.flameRadius[i - 2];
        const dHo = tr.flameRadius[i] - tr.flameRadius[i - 1];
        ratio = dHo / dPrev;
        break;
      }
    }
    // round 1: the hand-off step added ≈ 16 steps' worth of growth (0.48 mm vs 0.03 mm per step)
    expect(ratio).toBeGreaterThan(0.5);
    expect(ratio).toBeLessThan(2);
  });
});

describe('lumped wall temperatures (measured-data r2 #4)', () => {
  // round 1 used the fired standard-knock walls of Pal 2018 Table 3 at EVERY operating point
  it('motored walls sit a few K above the coolant; fired walls well above; the reference state reproduces Pal 2018 Table 3', () => {
    const mot = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, equivalenceRatio: 0 }, { combustionModel: 'none', warmupCycles: 6 });
    mot.runCycles(1);
    const Tc = CFR_RON_CONDITIONS.coolantTemperature;
    expect(mot.walls.headTemperature - Tc).toBeLessThan(20);
    expect(mot.walls.pistonTemperature - Tc).toBeLessThan(20);
    const ref = new CycleModel(
      CFR_F1,
      { ...CFR_RON_CONDITIONS, compressionRatio: 7.55, equivalenceRatio: 1 / 0.89, fuel: { kind: 'PRF', octaneNumber: 100 }, ambientPressure: 100240, ambientTemperature: 322, intakeMixtureTemperature: 302, sparkAdvanceDeg: 12.72 },
      { warmupCycles: 6 },
    );
    ref.runCycles(1);
    const w = CFR_F1.walls;
    for (const [a, b] of [
      [ref.walls.headTemperature, w.headTemperature],
      [ref.walls.pistonTemperature, w.pistonTemperature],
      [ref.walls.linerTemperature, w.linerTemperature],
      [ref.walls.exhaustValveTemperature, w.exhaustValveTemperature],
    ]) expect(Math.abs(a - b)).toBeLessThan(6); // (resistances fitted on the calibrated model)
    expect(ref.walls.headTemperature - mot.walls.headTemperature).toBeGreaterThan(80);
  });
});

describe('crevice zone (measured-data r2 #1)', () => {
  // round 1 had no crevice storage: the gas in the piston top-land/ring crevices burned with the charge
  it('RON: the crevice holds p·V_cr/(R_cr T_cr) through the closed phase, never negative, ≈ 2–6 % of the charge at peak pressure', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1 });
    expect(m.creviceOn).toBe(true);
    let worst = 0;
    let atPeak = 0;
    let pPeak = 0;
    let neg = false;
    const c0 = m.cycle;
    while (m.cycle === c0) {
      m.stepUntil(Infinity, Math.min(360, m.theta + 0.5));
      if (m.mode === MODE_OPEN) continue;
      const cl = m.closure;
      const mcu = m.y[I_CRU];
      const mcb = m.y[I_CRB];
      if (mcu < -1e-15 || mcb < -1e-15) neg = true;
      let ncc = 0;
      for (let k = 0; k < NS; k++) ncc += cl.ncc[k];
      const nR = mcu * cl.Ru + mcb * R_UNIVERSAL * ncc;
      const ref = (m.p * m.Vcr) / m.Tcr;
      worst = Math.max(worst, Math.abs(nR / ref - 1));
      if (m.p > pPeak) {
        pPeak = m.p;
        atPeak = (mcu + mcb) / m.mIvc;
      }
    }
    expect(neg).toBe(false);
    // (the constraint is kept by the ODE ṁ = c ṗ; the splits that move p discontinuously — kernel mass,
    // hand-off catch-up, burn-out merge, the NO energy — leave sub-percent offsets)
    expect(worst).toBeLessThan(5e-3);
    expect(atPeak).toBeGreaterThan(0.02);
    expect(atPeak).toBeLessThan(0.06);
  });
});

describe('knock intensity grows with the autoignited charge (knock-octane r2 #3, code review r2 #3)', () => {
  // round 1: the uniform source's projection vanished as the end gas approached the whole bore; PRF 40
  // (whole-charge autoignition) rang with MAPO 0, PRF 0 "autoignited" at IVC
  it('RON CR 6.43: MAPO does not fall as the octane number falls from 80 to 0', () => {
    let prev = 0;
    for (const on of [80, 60, 40, 20, 0]) {
      const s = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: on } }, { warmupCycles: 2 }).runCycles(1)[0];
      expect(s.knockOnsetDeg).toBeGreaterThan(-30);
      expect(s.mapo).toBeGreaterThanOrEqual(0.9 * prev);
      prev = s.mapo;
    }
    expect(prev).toBeGreaterThan(1e5);
  });
});

describe('snapshot heat-release rate (numerics r2 #11)', () => {
  // round 1 reported the step-mean burn rate: 6–11 % dependent on the step pattern
  it('RON: the heat-release rate at 5° ATDC agrees within 2 % between 0.25° and 0.0625° steps', () => {
    const at = (maxStepDeg: number): number => {
      const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1, knock: false, maxStepDeg, fineStepDeg: Math.min(0.05, maxStepDeg) });
      const c0 = m.cycle;
      m.stepUntil(Infinity, 360);
      while (m.cycle === c0 || m.theta < 5 - 1e-9) {
        m.stepUntil(Infinity, m.cycle === c0 ? 360 : 5);
        if (m.cycle > c0 && Math.abs(m.theta - 5) < 1e-9) break;
      }
      expect(m.mode).toBe(MODE_TWO);
      return m.hrr;
    };
    const a = at(0.25);
    const b = at(0.0625);
    expect(Math.abs(a / b - 1)).toBeLessThan(0.02);
  });
});

describe('ANL 2020s compression-ratio scale (measured-data r2 #5)', () => {
  it('clearance offset 7.38 cm³ (Hoth 2025 RON-95 height = CR 7.26 vs counter 805 here); SON PRF100 7.82 → ≈ 7.3', () => {
    expect(CFR_ANL2020S_CLEARANCE_OFFSET * 1e6).toBeCloseTo(7.38, 1);
    expect(cfrCompressionRatioFromAnl2020s(7.26)).toBeCloseTo(6.819, 2);
    expect(cfrCompressionRatioFromAnl2020s(7.82)).toBeGreaterThan(7.2);
    expect(cfrCompressionRatioFromAnl2020s(7.82)).toBeLessThan(7.4);
  });
});

describe('external EGR is tracked through the intake (code review r2 #5)', () => {
  it('10 % EGR: the trapped products dilution = residual + EGR and the reported residual excludes the EGR', () => {
    const a = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 3 });
    const b = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, egrFraction: 0.1 }, { warmupCycles: 3 });
    const sa = a.runCycles(2)[1];
    const sb = b.runCycles(2)[1];
    expect(b.yEgr).toBeGreaterThan(0.08);
    expect(b.yEgr).toBeLessThan(0.1);
    expect(b.xDil).toBeCloseTo(b.yRes, 12);
    expect(Math.abs(sb.residualFraction - sa.residualFraction)).toBeLessThan(0.02);
  });
});
