import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_valve_flow.json';
import type { ValveSpec } from '../core/engine-spec';
import {
  INFLOW_CURTAIN_CD,
  newValveJet,
  OUTFLOW_CURTAIN_CD,
  ValveFlowModel,
  valveDischargeCoefficient,
  valveEffectiveArea,
  valveFlowArea,
  valveFlowModel,
} from './valve-flow';

interface Case {
  name: string;
  headDiameter: number;
  seatInnerDiameter: number;
  seatAngle: number;
  stemDiameter: number;
  portDiameter: number;
  seatWidth: number;
  stage12Lift: number;
  rows: { lift: number; areaHeywood: number; stage: number; areaExactFrustum: number | null }[];
}
const cases = (fixture as unknown as { cases: Case[] }).cases;

const IN = 0.0254;
const mkSpec = (o: Partial<ValveSpec> = {}): ValveSpec => ({
  count: 1,
  headDiameter: 1.346 * IN + 0.125 * IN,
  seatInnerDiameter: 1.346 * IN,
  seatAngle: Math.PI / 4,
  stemDiameter: (11 / 32) * IN,
  maxLift: 0.246 * IN,
  openDeg: -350,
  closeDeg: -146,
  timingLiftThreshold: 0,
  position: [-0.02, 0],
  shroudArcDeg: 0,
  shroudDirection: Math.PI / 2,
  ...o,
});

const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

describe('valve minimum flow area (Heywood 1988 §6.3.1)', () => {
  for (const c of cases) {
    it(`${c.name}: stages 1–3 match the independent oracle; ≤ 0.3 % above the exact min frustum`, () => {
      const spec = mkSpec({
        headDiameter: c.headDiameter,
        seatInnerDiameter: c.seatInnerDiameter,
        seatAngle: c.seatAngle,
        stemDiameter: c.stemDiameter,
      });
      const m = new ValveFlowModel(spec, 'intake', { portDiameter: c.portDiameter });
      expect(rel(m.seatWidth, c.seatWidth)).toBeLessThan(1e-13);
      expect(rel(m.stage12Lift, c.stage12Lift)).toBeLessThan(1e-13);
      for (const r of c.rows) {
        expect(rel(m.flowArea(r.lift), r.areaHeywood)).toBeLessThan(1e-13);
        expect(m.stage(r.lift)).toBe(r.stage);
        if (r.areaExactFrustum !== null) {
          const gap = r.areaHeywood / r.areaExactFrustum - 1;
          expect(gap).toBeGreaterThan(-1e-9);
          expect(gap).toBeLessThan(3e-3);
        }
      }
    });
  }

  it('is continuous across the stage boundaries and zero when closed', () => {
    const m = new ValveFlowModel(mkSpec(), 'intake', { portDiameter: 0.03 });
    const L12 = m.stage12Lift;
    expect(rel(m.seatArea(L12 * (1 - 1e-12)), m.seatArea(L12))).toBeLessThan(1e-10);
    for (const L of [0, -1e-3, NaN]) expect(m.flowArea(L)).toBe(0);
    // port limit: flat and continuous
    let prev = 0;
    for (let L = 1e-6; L < 0.02; L += 1e-5) {
      const a = m.flowArea(L);
      expect(a).toBeGreaterThanOrEqual(prev - 1e-18);
      expect(a).toBeLessThanOrEqual(m.portArea + 1e-18);
      prev = a;
    }
    expect(m.flowArea(0.02)).toBe(m.portArea);
  });

  it('shroud masks the curtain fraction (CFR 180°) but not the port limit; count scales', () => {
    const open = new ValveFlowModel(mkSpec(), 'intake');
    const shr = new ValveFlowModel(mkSpec({ shroudArcDeg: 180 }), 'intake');
    expect(shr.openFraction).toBe(0.5);
    for (const L of [1e-4, 1e-3, 3e-3, 6e-3]) expect(rel(shr.flowArea(L), 0.5 * open.flowArea(L))).toBeLessThan(1e-14);
    expect(shr.flowArea(0.1)).toBe(open.flowArea(0.1));
    // finite lip: above it the masked arc opens
    const lip = new ValveFlowModel(mkSpec({ shroudArcDeg: 180 }), 'intake', { shroudHeight: 2e-3 });
    expect(lip.flowArea(1.5e-3)).toBe(shr.flowArea(1.5e-3));
    expect(lip.flowArea(4e-3)).toBeGreaterThan(shr.flowArea(4e-3));
    const two = new ValveFlowModel(mkSpec({ count: 2 }), 'exhaust');
    expect(two.flowArea(2e-3)).toBeCloseTo(2 * open.flowArea(2e-3), 18);
    expect(two.effectiveArea(2e-3, false)).toBeCloseTo(2 * new ValveFlowModel(mkSpec(), 'exhaust').effectiveArea(2e-3, false), 18);
  });
});

describe('discharge coefficients and effective area', () => {
  it('reproduces the tabulated curtain coefficients at the nodes (unshrouded)', () => {
    const m = new ValveFlowModel(mkSpec(), 'intake');
    const D = m.innerDiameter;
    INFLOW_CURTAIN_CD.lOverD.forEach((x, i) => {
      if (x === 0) return;
      const L = x * D;
      expect(m.curtainDischargeCoefficient(L, false)).toBeCloseTo(INFLOW_CURTAIN_CD.cd[i], 14);
      expect(rel(m.effectiveArea(L, false), INFLOW_CURTAIN_CD.cd[i] * Math.PI * D * L)).toBeLessThan(1e-13);
      // C_D referred to A_m times A_m gives the same effective area
      expect(rel(m.dischargeCoefficient(L, false) * m.flowArea(L), m.effectiveArea(L, false))).toBeLessThan(1e-13);
    });
    OUTFLOW_CURTAIN_CD.lOverD.forEach((x, i) => {
      if (x === 0) return;
      expect(m.curtainDischargeCoefficient(x * D, true)).toBeCloseTo(OUTFLOW_CURTAIN_CD.cd[i], 14);
    });
  });

  it('direction logic: intake reverse = exhaust forward = outflow table', () => {
    const s = mkSpec();
    for (const L of [5e-4, 3e-3, 6e-3]) {
      expect(valveDischargeCoefficient(s, L, true, 'intake')).toBe(valveDischargeCoefficient(s, L, false, 'exhaust'));
      expect(valveDischargeCoefficient(s, L, false, 'intake')).toBe(valveDischargeCoefficient(s, L, true, 'exhaust'));
    }
    expect(valveFlowModel(s, 'intake')).toBe(valveFlowModel(s, 'intake'));
    expect(valveFlowModel(s, 'intake')).not.toBe(valveFlowModel(s, 'exhaust'));
  });

  it('finite, non-NaN at zero lift, tiny lift and beyond the table; saturates the flow coefficient', () => {
    const s = mkSpec({ shroudArcDeg: 180 });
    const m = new ValveFlowModel(s, 'intake');
    for (const rev of [false, true]) {
      expect(m.dischargeCoefficient(0, rev)).toBeCloseTo(m.curtainDischargeCoefficient(0, rev) / Math.cos(s.seatAngle), 14);
      expect(rel(m.dischargeCoefficient(1e-12, rev), m.dischargeCoefficient(0, rev))).toBeLessThan(1e-9);
      expect(m.effectiveArea(0, rev)).toBe(0);
      expect(m.effectiveArea(5e-324, rev)).toBe(0);
      expect(Number.isFinite(m.dischargeCoefficient(5e-324, rev))).toBe(true);
      expect(valveEffectiveArea(s, 0, rev)).toBe(0);
      expect(Number.isFinite(m.dischargeCoefficient(0.05, rev))).toBe(true);
    }
    const D = m.innerDiameter;
    const cf = (L: number) => m.curtainDischargeCoefficient(L, false) * 4 * (L / D);
    expect(rel(cf(0.3 * D), cf(0.25 * D))).toBeLessThan(1e-14);
    expect(rel(cf(0.6 * D), cf(0.25 * D))).toBeLessThan(1e-14);
    expect(valveFlowArea(s, 3e-3)).toBe(m.flowArea(3e-3));
  });
});

describe('shrouded-valve inflow jet (swirl source)', () => {
  it('matches a direct quadrature of r × v over the open arc; zero when unshrouded or for outflow', () => {
    const spec = mkSpec({ shroudArcDeg: 180, shroudDirection: Math.PI / 2, position: [-0.02, 0.003] });
    const m = new ValveFlowModel(spec, 'intake');
    const jet = newValveJet();
    const mdot = 3e-3;
    const vj = 40;
    m.inflowJet(mdot, vj, jet);
    // quadrature: open arc centred on ψ with half-angle π/2; uniform mass flux; velocity radial
    const vh = vj * Math.cos(spec.seatAngle);
    const R = spec.headDiameter / 2;
    const [x0, z0] = spec.position;
    const N = 100000;
    let H = 0;
    for (let i = 0; i < N; i++) {
      const phi = spec.shroudDirection - Math.PI / 2 + ((i + 0.5) / N) * Math.PI;
      const x = x0 + R * Math.cos(phi);
      const z = z0 + R * Math.sin(phi);
      const vx = vh * Math.cos(phi);
      const vz = vh * Math.sin(phi);
      H += (mdot / N) * (z * vx - x * vz); // (r × v)_y
    }
    expect(rel(jet.angularMomentumFlux, H)).toBeLessThan(1e-8);
    expect(jet.angularMomentumFlux).toBeGreaterThan(0); // open side facing +z at x < 0 → swirl about +y
    expect(jet.kineticEnergyFlux).toBeCloseTo(0.5 * mdot * vj * vj, 14);
    const open = new ValveFlowModel(mkSpec(), 'intake');
    expect(open.inflowJet(mdot, vj, jet).angularMomentumFlux).toBe(0);
    expect(m.inflowJet(-mdot, vj, jet).angularMomentumFlux).toBe(0);
    expect(jet.kineticEnergyFlux).toBe(0);
  });
});
