import { describe, expect, it } from 'vitest';
import { GapTransition, SparkGap } from './discharge';

describe('SparkGap voltage–current characteristic', () => {
  it('glow: V = V_cf + 40.46·l[mm]·p[bar]^0.51·I^−0.32 (Kim & Anderson 1995)', () => {
    const g = new SparkGap(1e-3);
    g.pressure = 1e5;
    g.mode = 'glow';
    expect(g.gapVoltage(0.1)).toBeCloseTo(g.opts.glowCathodeFall + 40.46 * Math.pow(0.1, -0.32), 9);
    // falling characteristic, slope = n·V_col/I
    const i = 0.05;
    const num = (g.gapVoltage(i * 1.0001) - g.gapVoltage(i * 0.9999)) / (i * 0.0002);
    expect(g.gapVoltageSlope(i) / num).toBeCloseTo(1, 5);
    expect(g.gapVoltageSlope(i)).toBeLessThan(0);
  });

  it('reproduces Herweg & Maly TCI glow levels (> 400 V at 5 bar, 1 mm, ≤ 80 mA)', () => {
    const g = new SparkGap(1e-3);
    g.pressure = 5e5;
    g.mode = 'glow';
    for (const i of [0.08, 0.05, 0.02]) expect(g.gapVoltage(i)).toBeGreaterThan(400);
  });

  it('arc: low voltage (< 100 V, Heywood 1988 p. 428) with a 15 V sheath (Shaffer et al. 2023)', () => {
    const g = new SparkGap(0.5e-3);
    g.pressure = 10e5;
    g.mode = 'arc';
    expect(g.sheathVoltage('arc')).toBe(15);
    expect(g.gapVoltage(0.2)).toBeLessThan(100);
    expect(g.gapVoltage(1)).toBeLessThan(g.gapVoltage(0.2));
  });

  it('emergent glow energy transfer to the gas ≈ Maly & Vogel 30 % (25–50 % at 1–15 bar)', () => {
    const g = new SparkGap(0.8e-3);
    g.mode = 'glow';
    for (const p of [1e5, 5e5, 10e5, 15e5]) {
      g.pressure = p;
      const v = g.gapVoltage(0.05);
      const eta = 1 - g.sheathVoltage('glow') / v;
      expect(eta).toBeGreaterThan(0.2);
      expect(eta).toBeLessThan(0.55);
    }
  });
});

describe('SparkGap breakdown event', () => {
  it('splits the capacitor dump into breakdown (½C_plug V², 94 %) and capacitive arc (50 %)', () => {
    const g = new SparkGap(0.5e-3);
    g.pressure = 10e5;
    g.breakdownVoltage = 10e3;
    const c2 = 60e-12;
    const vAfter = g.startConduction(10e3, 0.05, c2, 1.0);
    expect(g.mode).toBe('glow');
    const eDump = 0.5 * c2 * (1e8 - vAfter * vAfter);
    const eBd = 0.5 * 10e-12 * 1e8; // 0.5 mJ
    expect(g.energyBreakdown).toBeCloseTo(eBd, 12);
    expect(g.energyCapacitiveArc).toBeCloseTo(eDump - eBd, 12);
    expect(g.stepImpulseGasEnergy).toBeCloseTo(0.94 * eBd + 0.5 * (eDump - eBd), 12);
    expect(g.energyToGas + g.energyToElectrodes + g.energyRadiated).toBeCloseTo(g.energyTotal, 14);
    expect(g.breakdownCount).toBe(1);
    expect(g.firstBreakdownTime).toBe(1.0);
  });

  it('starts in the arc mode above the arc–glow current (0.1 A) and hands over to glow below it', () => {
    const g = new SparkGap(0.5e-3);
    g.pressure = 10e5;
    g.breakdownVoltage = 10e3;
    g.startConduction(10e3, 0.3, 60e-12, 0);
    expect(g.mode).toBe('arc');
    expect(g.afterConductingStep(0.2, 0)).toBe(GapTransition.None);
    expect(g.afterConductingStep(0.09, 0)).toBe(GapTransition.Recharge);
    expect(g.mode).toBe('glow');
    expect(g.awaitingReentry).toBe(true);
    // re-entry at the glow voltage: no energy dump
    const e0 = g.energyTotal;
    const v = g.startConduction(g.threshold(0.09), 0.09, 60e-12, 0);
    expect(v).toBeCloseTo(g.modeVoltage('glow', 0.09), 9);
    expect(g.energyTotal).toBe(e0);
    expect(g.afterConductingStep(0.001, 0)).toBe(GapTransition.Extinguish);
    expect(g.mode).toBe('open');
  });

  it('restrikes when the flow-stretched channel needs more than a fresh breakdown', () => {
    const g = new SparkGap(0.5e-3);
    g.pressure = 10e5;
    g.breakdownVoltage = 3e3;
    g.flowVelocity = 10;
    g.startConduction(3e3, 0.05, 60e-12, 0);
    // stretch the channel for 1 ms at 10 m/s: l = d + 2 u t = 20.5 mm
    for (let i = 0; i < 1000; i++) g.accumulate(0, 0, 0.05, 1e-6);
    expect(g.channelLength).toBeCloseTo(Math.min(20 * 0.5e-3, 0.5e-3 + 2 * 10 * 1e-3), 9);
    expect(g.afterConductingStep(0.05, 0)).toBe(GapTransition.Restrike);
    expect(g.channelLength).toBe(0.5e-3);
  });
});

describe('SparkGap consistency constraints (reviewer fixes)', () => {
  it('a (re)breakdown needs ≥ (1+h)·V_glow(d, I_ext), above every sustainable gap voltage', () => {
    const g = new SparkGap(0.18e-3);
    g.pressure = 0.55e5;
    g.breakdownVoltage = 300; // hot, low-density gap gas: the breakdown law is at the cathode-fall scale
    const vMin = g.minimumBreakdownVoltage();
    g.mode = 'glow';
    const vGlowExt = g.modeVoltage('glow', g.opts.extinctionCurrent);
    expect(vMin).toBeCloseTo((1 + g.opts.transitionHysteresis) * vGlowExt, 9);
    expect(g.threshold(0.01)).toBe(vMin);
    for (const i of [0.003, 0.01, 0.05, 0.2]) expect(g.modeVoltage('glow', i)).toBeLessThan(g.threshold(i));
    g.mode = 'open';
    // at engine densities the breakdown voltage itself is the threshold
    g.breakdownVoltage = 9800;
    g.pressure = 10.85e5;
    expect(g.threshold(0.05)).toBe(9800);
  });

  it('never returns energy at breakdown (V_after ≤ V_bd) and books a capacitive spark only when I ≤ I_ext', () => {
    const g = new SparkGap(0.107e-3);
    g.pressure = 1.75e5;
    g.breakdownVoltage = 345;
    // regression: V_after = V_glow(1.1 mA) (+ regularised I^−0.32) exceeded the 345 V breakdown
    const vb = g.threshold(1.1e-3);
    const vAfter = g.startConduction(vb, 1.1e-3, 22e-12, 0);
    expect(vAfter).toBeLessThanOrEqual(vb);
    expect(g.energyTotal).toBeGreaterThanOrEqual(0);
    expect(g.conducting).toBe(false); // circuit current below the sustaining current
    expect(g.breakdownCount).toBe(1);
    expect(g.extinctionTime).toBe(0);
    expect(g.energyToGas + g.energyToElectrodes + g.energyRadiated).toBeCloseTo(g.energyTotal, 18);
    // with a sustainable current the discharge continues in glow at V_glow(I) < V_bd
    const h = new SparkGap(0.107e-3);
    h.pressure = 1.75e5;
    h.breakdownVoltage = 345;
    const v2 = h.startConduction(h.threshold(0.03), 0.03, 22e-12, 0);
    expect(h.conducting).toBe(true);
    expect(v2).toBeCloseTo(h.modeVoltage('glow', 0.03), 9);
  });

  it('a fresh (unstretched) channel never re-triggers a restrike, even in hot low-density gas', () => {
    const g = new SparkGap(0.18e-3);
    g.pressure = 0.55e5;
    g.breakdownVoltage = 300;
    g.startConduction(g.threshold(0.05), 0.05, 46e-12, 0);
    for (const i of [0.1, 0.05, 0.01, 0.003]) expect(g.afterConductingStep(i, 0)).not.toBe(GapTransition.Restrike);
  });
});
