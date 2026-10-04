import { describe, expect, it } from 'vitest';
import { MODEL_T_IGNITION } from '../engines/model-t';
import { magnetoEmf, PrimarySupply } from './source';

const M = MODEL_T_IGNITION.magneto;
const rad = (rpm: number): number => (rpm * 2 * Math.PI) / 60;

describe('magneto EMF (Ford Model T flywheel magneto)', () => {
  it('e = k ω sin(N(θ − φ)): 8 cycles per revolution, zero at φ, peak k ω a quarter period later', () => {
    const w = rad(1000);
    expect(magnetoEmf(M, M.phaseDeg, w)).toBeCloseTo(0, 12);
    expect(magnetoEmf(M, M.phaseDeg + 45 / 4, w)).toBeCloseTo(M.emfConstant * w, 12);
    // 360°/8 = 45° per electrical cycle: every 180° firing offset sees the same phase
    for (const th of [-300, -17.3, 0, 41, 222]) {
      expect(magnetoEmf(M, th + 45, w)).toBeCloseTo(magnetoEmf(M, th, w), 9);
      expect(magnetoEmf(M, th + 180, w)).toBeCloseTo(magnetoEmf(M, th, w), 9);
    }
  });

  it('meets the Ford service minimum (≥ 7 V at 400 rpm) and the healthy-magneto data (≈ 10 V rms at 407 rpm)', () => {
    const vrms = (rpm: number): number => (M.emfConstant * rad(rpm)) / Math.SQRT2;
    expect(vrms(400)).toBeGreaterThan(7);
    expect(vrms(407)).toBeGreaterThan(9);
    expect(vrms(407)).toBeLessThan(11);
    expect(vrms(1220)).toBeGreaterThan(24.4 * 0.9);
    expect(vrms(1220)).toBeLessThan(24.4 * 1.25);
  });

  it('source impedance R_s + jω_e L_s fits the 1976 table |Z| = V/I (Payne; L_s least squares) within 0.2 Ω', () => {
    // rpm, V (open circuit), A (taken as short-circuit current) — Gas Engine Magazine, 1 Sep 1976
    const table: [number, number, number][] = [[400, 9.8, 7.9], [600, 14.4, 8.5], [800, 18.8, 8.8], [1000, 22.8, 8.9], [1200, 26.2, 9.0]];
    for (const [rpm, V, I] of table) {
      const z = Math.hypot(M.internalResistance, M.cyclesPerRevolution * rad(rpm) * M.internalInductance);
      expect(Math.abs(z - V / I)).toBeLessThan(0.2);
    }
    // short-circuit current at speed → k/(N L_s) ("gives a constant current", Ford test-stand pamphlet)
    const iInf = M.emfConstant / (M.cyclesPerRevolution * M.internalInductance) / Math.SQRT2;
    expect(iInf).toBeGreaterThan(7.9);
    expect(iInf).toBeLessThan(10);
  });
});

describe('PrimarySupply', () => {
  it('interpolates angle and speed linearly over the step and switches BAT/MAG only when applied', () => {
    const s = new PrimarySupply(MODEL_T_IGNITION.battery, M, 'magneto');
    s.setStep(1, 0.01, 10, 70, rad(1000), rad(1200));
    expect(s.thetaAt(1.005)).toBeCloseTo(40, 9);
    expect(s.omegaAt(1.005)).toBeCloseTo(rad(1100), 9);
    expect(s.emf(1.005)).toBeCloseTo(magnetoEmf(M, 40, rad(1100)), 9);
    expect(s.inductance).toBe(M.internalInductance);
    expect(s.resistance).toBe(M.internalResistance);
    s.request('battery');
    expect(s.kind).toBe('magneto');
    s.applyRequested();
    expect(s.emf(1.005)).toBe(MODEL_T_IGNITION.battery.voltage);
    expect(s.inductance).toBe(0);
    expect(s.resistance).toBe(MODEL_T_IGNITION.battery.internalResistance);
    expect(s.magnetoEmfAt(1.005)).toBeCloseTo(magnetoEmf(M, 40, rad(1100)), 9);
    expect(() => PrimarySupply.dc(6).request('magneto')).toThrow();
  });
});
