import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_throttle.json';
import { DEG } from '../core/constants';
import {
  THROTTLE_DEFAULTS,
  throttleArea,
  throttleGeometricArea,
  throttlePlateAngle,
  venturiPressureLoss,
  venturiThroatDepression,
} from './throttle';

interface Case {
  diameter: number;
  shaftRatio: number;
  closedAngleDeg: number;
  rows: { plateAngleDeg: number; area: number }[];
}
const cases = (fixture as unknown as { cases: Case[] }).cases;

describe('butterfly throttle area', () => {
  for (const c of cases) {
    it(`D = ${c.diameter} m, a = ${c.shaftRatio}, ψ0 = ${c.closedAngleDeg}°: closed form = numerical integration`, () => {
      const A0 = 0.25 * Math.PI * c.diameter ** 2;
      for (const r of c.rows) {
        const A = throttleGeometricArea(c.diameter, r.plateAngleDeg * DEG, c.shaftRatio, c.closedAngleDeg * DEG);
        expect(Math.abs(A - r.area) / A0).toBeLessThan(1e-11);
      }
    });
  }

  it('closed → 0 (+ leakage), monotone in opening, wide open = bore − shaft', () => {
    const D = 0.05;
    const leak = 1e-6;
    expect(throttleArea(D, 0, { leakageArea: leak })).toBe(leak);
    expect(throttleArea(D, -1, { leakageArea: 0 })).toBe(0);
    let prev = -1;
    for (let u = 0; u <= 1.0000001; u += 0.001) {
      const A = throttleArea(D, u, { leakageArea: 0 });
      expect(A).toBeGreaterThanOrEqual(prev - 1e-18);
      prev = A;
    }
    const a = THROTTLE_DEFAULTS.shaftRatio;
    const wide = 0.25 * Math.PI * D * D - D * D * 0.5 * (a * Math.sqrt(1 - a * a) + Math.asin(a));
    expect(throttleArea(D, 1, { leakageArea: 0 })).toBeCloseTo(wide, 15);
    expect(throttleArea(D, 2, { leakageArea: 0 })).toBe(throttleArea(D, 1, { leakageArea: 0 }));
    // default leakage = π D c
    expect(throttleArea(D, 0)).toBeCloseTo(Math.PI * D * THROTTLE_DEFAULTS.leakageClearance, 18);
    expect(throttlePlateAngle(0.5, 10 * DEG)).toBeCloseTo(50 * DEG, 14);
    // no shaft: A = A0 (1 − cosψ/cosψ0)
    const psi = 40 * DEG;
    expect(throttleGeometricArea(D, psi, 0, 10 * DEG)).toBeCloseTo(0.25 * Math.PI * D * D * (1 - Math.cos(psi) / Math.cos(10 * DEG)), 16);
  });

  it('venturi: throat depression and loss (CFR-like magnitude check)', () => {
    const IN = 0.0254;
    const dt = (9 / 16) * IN;
    const rho = 1.09;
    const mdot = 2.7e-3; // mean CFR air+fuel flow at 600 rpm (≈ ρ V_d n/2 η_v)
    const dp = venturiThroatDepression(mdot, rho, dt);
    const v = mdot / (rho * 0.25 * Math.PI * dt * dt);
    expect(dp).toBeCloseTo(0.5 * rho * v * v, 10);
    expect(dp).toBeGreaterThan(80);
    expect(dp).toBeLessThan(200);
    expect(venturiPressureLoss(-mdot, rho, dt)).toBeCloseTo(-0.15 * dp, 10);
  });
});
