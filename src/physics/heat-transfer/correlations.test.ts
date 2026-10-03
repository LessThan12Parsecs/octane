import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/heat_correlations.json';
import { SIGMA_SB } from '../core/constants';
import { NS, SP, type SpeciesName } from '../core/species';
import { mixDensity, mixThermalConductivity, mixViscosity } from '../thermo';
import {
  ANNAND_RADIATION,
  annandCoefficient,
  annandFlux,
  hohenbergCoefficient,
  WOSCHNI_CONSTANTS,
  woschniCoefficient,
  woschniMotoredPressure,
  woschniVelocity,
  type WoschniInputs,
  type WoschniPhase,
} from './woschni';

interface WRow {
  phase: WoschniPhase;
  bore: number;
  p: number;
  T: number;
  pMotored: number;
  meanPistonSpeed: number;
  swirlRatio: number;
  Vs: number;
  p1: number;
  T1: number;
  V1: number;
  V: number;
  w: number;
  hWoschni1967: number;
  hHeywood: number;
  hHohenberg: number;
}
interface ARow {
  label: string;
  X: Record<string, number>;
  T: number;
  p: number;
  Tw: number;
  meanPistonSpeed: number;
  bore: number;
  a: number;
  c: number;
  k: number;
  mu: number;
  rho: number;
  h: number;
  flux: number;
}
const F = fixture as unknown as { woschniHohenberg: WRow[]; annand: ARow[] };
const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

describe('Woschni and Hohenberg vs independent re-implementation in published units', () => {
  const inp: WoschniInputs = { bore: 0, pressure: 0, temperature: 0, meanPistonSpeed: 0, phase: 'compression' };
  for (const r of F.woschniHohenberg) {
    it(`${r.phase}, p = ${r.p / 1e5} bar, T = ${r.T} K, swirl ${r.swirlRatio}`, () => {
      inp.bore = r.bore;
      inp.pressure = r.p;
      inp.temperature = r.T;
      inp.meanPistonSpeed = r.meanPistonSpeed;
      inp.phase = r.phase;
      inp.motoredPressure = r.pMotored;
      inp.displacedVolume = r.Vs;
      inp.refPressure = r.p1;
      inp.refTemperature = r.T1;
      inp.refVolume = r.V1;
      inp.swirlRatio = r.swirlRatio;
      inp.variant = undefined;
      expect(rel(woschniVelocity(inp), r.w)).toBeLessThan(1e-13);
      expect(rel(woschniCoefficient(inp), r.hWoschni1967)).toBeLessThan(1e-12);
      inp.variant = 'heywood1988';
      expect(rel(woschniCoefficient(inp), r.hHeywood)).toBeLessThan(1e-12);
      expect(rel(hohenbergCoefficient(r.V, r.p, r.T, r.meanPistonSpeed), r.hHohenberg)).toBeLessThan(1e-12);
    });
  }

  it('unit systems: 0.820 kW/(m²K)·MPa^−0.8 = 3.26 W/(m²K)·kPa^−0.8 (to 0.14 %); exponent difference only', () => {
    const K = WOSCHNI_CONSTANTS;
    expect(rel(K.cWoschni1967 / Math.pow(1000, 0.8), K.cHeywood)).toBeLessThan(1.5e-3);
    const i: WoschniInputs = { bore: 0.0826, pressure: 3e6, temperature: 2000, meanPistonSpeed: 2.3, phase: 'compression' };
    const hw = woschniCoefficient(i);
    const hh = woschniCoefficient({ ...i, variant: 'heywood1988' });
    const expected = (K.cHeywood / (K.cWoschni1967 / Math.pow(1000, 0.8))) * Math.pow(2000, -0.02);
    expect(rel(hh / hw, expected)).toBeLessThan(1e-12);
  });

  it('combustion term: zero below the motored pressure; motored-pressure helper', () => {
    const base: WoschniInputs = {
      bore: 0.0826,
      pressure: 1e6,
      temperature: 800,
      meanPistonSpeed: 2.3,
      phase: 'combustion',
      motoredPressure: 1.2e6,
      displacedVolume: 6.1e-4,
      refPressure: 1e5,
      refTemperature: 340,
      refVolume: 6.6e-4,
    };
    expect(woschniVelocity(base)).toBeCloseTo(2.28 * 2.3, 14);
    expect(woschniVelocity({ ...base, phase: 'compression' })).toBeCloseTo(2.28 * 2.3, 14);
    expect(woschniVelocity({ ...base, phase: 'gas-exchange' })).toBeCloseTo(6.18 * 2.3, 14);
    expect(woschniVelocity({ ...base, pressure: 2e6, c2: WOSCHNI_CONSTANTS.c2DividedChamber })).toBeCloseTo(
      2.28 * 2.3 + (6.22e-3 * 6.1e-4 * 340 * 0.8e6) / (1e5 * 6.6e-4),
      12,
    );
    expect(woschniMotoredPressure(1e5, 6.6e-4, 1e-4, 1.3)).toBeCloseTo(1e5 * Math.pow(6.6, 1.3), 6);
    expect(woschniCoefficient({ ...base, pressure: 0 })).toBe(0);
    // a missing reference/motored input must not silently drop the combustion term
    for (const key of ['motoredPressure', 'displacedVolume', 'refPressure', 'refTemperature', 'refVolume'] as const) {
      expect(() => woschniVelocity({ ...base, pressure: 2e6, [key]: undefined })).toThrow(RangeError);
      expect(() => woschniCoefficient({ ...base, pressure: 2e6, [key]: undefined })).toThrow(RangeError);
    }
    expect(woschniVelocity({ phase: 'compression', bore: 0.08, pressure: 2e6, temperature: 800, meanPistonSpeed: 2.3 })).toBeCloseTo(2.28 * 2.3, 14);
    // 'expansion' (core CylinderPhase) = 'combustion' for Woschni
    expect(woschniVelocity({ ...base, pressure: 2e6, phase: 'expansion' })).toBe(woschniVelocity({ ...base, pressure: 2e6 }));
    expect(() => woschniVelocity({ ...base, pressure: 2e6, phase: 'expansion', refVolume: undefined })).toThrow(RangeError);
  });
});

describe('Annand (1963) with our transport properties vs Cantera transport', () => {
  for (const r of F.annand) {
    it(`${r.label}`, () => {
      const X = new Float64Array(NS);
      let s = 0;
      for (const [k, v] of Object.entries(r.X)) {
        X[SP[k as SpeciesName]] = v;
        s += v;
      }
      for (let k = 0; k < NS; k++) X[k] /= s;
      const mu = mixViscosity(X, r.T);
      const kc = mixThermalConductivity(X, r.T);
      const rho = mixDensity(X, r.T, r.p);
      expect(rel(rho, r.rho)).toBeLessThan(1e-9);
      // exact correlation with the oracle's own properties
      expect(rel(annandCoefficient(r.bore, r.meanPistonSpeed, r.rho, r.mu, r.k, r.a), r.h)).toBeLessThan(1e-12);
      expect(rel(annandFlux(r.h, r.T, r.Tw, r.c), r.flux)).toBeLessThan(1e-12);
      // with our thermo/transport: μ within 0.6 %, but the burned-gas mixture conductivity at
      // 2000–2400 K is ≈ 2 % below Cantera's (thermo/transport.ts, not this module) → h within 3 %
      const h = annandCoefficient(r.bore, r.meanPistonSpeed, rho, mu, kc, r.a);
      expect(rel(mu, r.mu)).toBeLessThan(6e-3);
      expect(rel(kc, r.k)).toBeLessThan(2.5e-2);
      expect(rel(h, r.h)).toBeLessThan(3e-2);
    });
  }
  it('radiation constants and zero flux at T = T_w', () => {
    expect(ANNAND_RADIATION.sparkIgnition).toBeCloseTo(0.075 * SIGMA_SB, 20);
    expect(annandFlux(500, 480, 480, ANNAND_RADIATION.sparkIgnition)).toBe(0);
  });
});
