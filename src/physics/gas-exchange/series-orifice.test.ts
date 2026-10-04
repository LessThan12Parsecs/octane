import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_series_orifice.json';
import { CFR_F1 } from '../engines/cfr';
import { MODEL_T } from '../engines/model-t';
import {
  CarburettorFlowModel,
  hasSeparateThrottle,
  throttleEffectiveArea,
  throttlePlateOptions,
} from './carburettor';
import { newOrificeFlow, orificeFlow } from './orifice';
import { newSeriesOrificeFlow, seriesOrificeFlow } from './series-orifice';
import { throttleArea } from './throttle';
import { seriesEffectiveArea, SLOT_CONTRACTION_COEFFICIENT } from './valve-flow';

interface Case {
  name: string;
  cdaA: number;
  cdaB: number;
  pa: number;
  Ta: number;
  Ra: number;
  gammaA: number;
  pb: number;
  Tb: number;
  Rb: number;
  gammaB: number;
  mdot: number;
  intermediatePressure: number;
  dmdotdpa: number;
  dmdotdpb: number;
  chokedFirst: boolean;
  chokedSecond: boolean;
}
const cases = (fixture as unknown as { cases: Case[] }).cases;
const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

describe('two compressible restrictions in series (oracle gasex_series_orifice.py)', () => {
  for (const c of cases) {
    it(`${c.name}: flow, node pressure, derivatives and choking = scipy two-nozzle solve`, () => {
      const out = newSeriesOrificeFlow();
      const m = seriesOrificeFlow(c.cdaA, c.cdaB, c.pa, c.Ta, c.Ra, c.gammaA, c.pb, c.Tb, c.Rb, c.gammaB, out);
      expect(m).toBe(out.mdot);
      expect(rel(m, c.mdot)).toBeLessThan(1e-10);
      expect(rel(out.intermediatePressure, c.intermediatePressure)).toBeLessThan(1e-11);
      const scale = Math.abs(c.dmdotdpa) + Math.abs(c.dmdotdpb);
      expect(Math.abs(out.dmdotdpa - c.dmdotdpa)).toBeLessThan(2e-6 * scale);
      expect(Math.abs(out.dmdotdpb - c.dmdotdpb)).toBeLessThan(2e-6 * scale);
      const fwd = c.pa >= c.pb;
      expect(fwd ? out.chokedA : out.chokedB).toBe(c.chokedFirst);
      expect(fwd ? out.chokedB : out.chokedA).toBe(c.chokedSecond);
      expect(out.upstream).toBe(fwd ? 0 : 1);
      // both restrictions carry the same flow at the node (upstream stagnation T, R, γ)
      const [Tu, Ru, gu] = fwd ? [c.Ta, c.Ra, c.gammaA] : [c.Tb, c.Rb, c.gammaB];
      const o = newOrificeFlow();
      const pu = fwd ? c.pa : c.pb;
      const pd = fwd ? c.pb : c.pa;
      const m1 = orificeFlow(fwd ? c.cdaA : c.cdaB, pu, Tu, Ru, gu, out.intermediatePressure, Tu, Ru, gu, o);
      const m2 = orificeFlow(fwd ? c.cdaB : c.cdaA, out.intermediatePressure, Tu, Ru, gu, pd, Tu, Ru, gu, o);
      expect(rel(m1, Math.abs(m))).toBeLessThan(1e-12);
      expect(rel(m2, Math.abs(m))).toBeLessThan(1e-12);
      // warm start: the second solve converges at once to the same state
      const it0 = out.iterations;
      seriesOrificeFlow(c.cdaA, c.cdaB, c.pa, c.Ta, c.Ra, c.gammaA, c.pb, c.Tb, c.Rb, c.gammaB, out);
      expect(out.iterations).toBeLessThanOrEqual(Math.min(it0, 2));
      expect(rel(out.mdot, c.mdot)).toBeLessThan(1e-10);
    });
  }

  it('limits: incompressible series area at small Δp, single orifice for an absent restriction, no flow', () => {
    const out = newSeriesOrificeFlow();
    const o = newOrificeFlow();
    const cA = 1.5e-4;
    const cB = 2.5e-4;
    // Δp = 0.4 % of p (outside the regularisation band): ≈ incompressible series combination
    const m = seriesOrificeFlow(cA, cB, 1e5, 300, 287, 1.4, 0.996e5, 300, 287, 1.4, out);
    const mi = orificeFlow(seriesEffectiveArea(cA, cB), 1e5, 300, 287, 1.4, 0.996e5, 300, 287, 1.4, o);
    expect(rel(m, mi)).toBeLessThan(3e-3);
    // ∞ restriction → exactly the single orifice
    const s = seriesOrificeFlow(cA, Infinity, 1e5, 300, 287, 1.4, 0.7e5, 310, 287, 1.4, out);
    expect(s).toBe(orificeFlow(cA, 1e5, 300, 287, 1.4, 0.7e5, 310, 287, 1.4, o));
    expect(out.dmdotdpb).toBe(o.dmdotdpb);
    expect(seriesOrificeFlow(Infinity, cB, 0.7e5, 300, 287, 1.4, 1e5, 310, 287, 1.4, out)).toBe(orificeFlow(cB, 0.7e5, 300, 287, 1.4, 1e5, 310, 287, 1.4, o));
    expect(out.chokedB).toBe(o.choked);
    expect(seriesOrificeFlow(0, cB, 1e5, 300, 287, 1.4, 0.5e5, 300, 287, 1.4, out)).toBe(0);
    expect(out.dmdotdpa).toBe(0);
    // equal pressures: no flow; slopes = series of the two regularised conductances
    const z = seriesOrificeFlow(cA, cB, 1e5, 300, 287, 1.4, 1e5, 300, 287, 1.4, out);
    expect(z).toBe(0);
    orificeFlow(cA, 1e5, 300, 287, 1.4, 1e5, 300, 287, 1.4, o);
    const g1 = o.dmdotdpa;
    orificeFlow(cB, 1e5, 300, 287, 1.4, 1e5, 300, 287, 1.4, o);
    const g2 = o.dmdotdpa;
    expect(rel(out.dmdotdpa, (g1 * g2) / (g1 + g2))).toBeLessThan(1e-12);
    expect(rel(-out.dmdotdpb, (g1 * g2) / (g1 + g2))).toBeLessThan(1e-12);
  });

  it('seeded stress sweep: finite, below either restriction alone, exact derivatives, bounded iterations', () => {
    let s = 7;
    const r = () => {
      s = (s * 1103515245 + 12345) % 2147483648;
      return s / 2147483648;
    };
    const lu = (a: number, b: number) => Math.exp(Math.log(a) + r() * Math.log(b / a));
    const out = newSeriesOrificeFlow();
    const tmp = newSeriesOrificeFlow();
    const o = newOrificeFlow();
    let worst = 0;
    let maxIt = 0;
    for (let i = 0; i < 4000; i++) {
      const pa = lu(1e3, 3e6);
      const pb = r() < 0.2 ? pa * (1 + (r() - 0.5) * 4e-3) : lu(1e3, 3e6);
      const Ta = lu(250, 2500);
      const Tb = lu(250, 2500);
      const ga = 1.2 + 0.25 * r();
      const gb = 1.2 + 0.25 * r();
      const cA = lu(1e-8, 1e-2);
      const cB = lu(1e-8, 1e-2);
      if (r() < 0.5) out.intermediateFraction = Number.NaN; // cold and warm starts
      const m = seriesOrificeFlow(cA, cB, pa, Ta, Ra(), ga, pb, Tb, Rb(), gb, out);
      expect(Number.isFinite(m) && Number.isFinite(out.dmdotdpa) && Number.isFinite(out.dmdotdpb)).toBe(true);
      expect(Math.sign(m)).toBe(pa > pb ? 1 : pa < pb ? -1 : 0);
      expect(out.dmdotdpa).toBeGreaterThanOrEqual(0);
      expect(out.dmdotdpb).toBeLessThanOrEqual(0);
      const single = Math.min(Math.abs(orificeFlow(cA, pa, Ta, 287, ga, pb, Tb, 290, gb, o)), Math.abs(orificeFlow(cB, pa, Ta, 287, ga, pb, Tb, 290, gb, o)));
      expect(Math.abs(m)).toBeLessThanOrEqual(single * (1 + 1e-12));
      maxIt = Math.max(maxIt, out.iterations);
      // derivative check by central differences (away from Δp = 0)
      if (Math.abs(pa - pb) > 0.01 * Math.max(pa, pb)) {
        const h = 1e-6 * Math.max(pa, pb);
        const fd = (seriesOrificeFlow(cA, cB, pa + h, Ta, 287, ga, pb, Tb, 290, gb, tmp) - seriesOrificeFlow(cA, cB, pa - h, Ta, 287, ga, pb, Tb, 290, gb, tmp)) / (2 * h);
        const sc = Math.abs(out.dmdotdpa) + Math.abs(out.dmdotdpb) + 1e-30;
        worst = Math.max(worst, Math.abs(fd - out.dmdotdpa) / sc);
      }
    }
    expect(worst).toBeLessThan(1e-4); // FD across the C¹ choking kink
    expect(maxIt).toBeLessThan(30);
    function Ra() {
      return 287;
    }
    function Rb() {
      return 290;
    }
  });

  it('is continuous and monotone in the downstream pressure through choking and reversal', () => {
    const out = newSeriesOrificeFlow();
    let prev = Infinity;
    for (let pb = 1000; pb <= 1.2e5; pb += 137) {
      const m = seriesOrificeFlow(1.6e-4, 5e-6, 1.01325e5, 288, 287, 1.4, pb, 300, 287, 1.4, out);
      expect(m).toBeLessThanOrEqual(prev + 1e-15);
      expect(out.dmdotdpb).toBeLessThanOrEqual(0);
      expect(out.iterations).toBeLessThan(25);
      prev = m;
    }
  });
});

describe('carburettor: venturi + butterfly in series', () => {
  it('throttle effective area: Kirchhoff slot contraction at small openings, Borda–Carnot recovery wide open', () => {
    const D = 0.0254;
    const opts = throttlePlateOptions(MODEL_T.manifolds.throttle);
    const Ab = 0.25 * Math.PI * D * D;
    const a0 = throttleArea(D, 0, opts);
    expect(rel(throttleEffectiveArea(D, 0, opts), SLOT_CONTRACTION_COEFFICIENT * a0 / (1 - (SLOT_CONTRACTION_COEFFICIENT * a0) / Ab))).toBeLessThan(1e-14);
    expect(throttleEffectiveArea(D, 0, opts) / a0).toBeCloseTo(SLOT_CONTRACTION_COEFFICIENT, 2); // A_o ≪ A_bore
    const aw = throttleArea(D, 1, opts);
    const cw = throttleEffectiveArea(D, 1, opts) / aw;
    expect(cw).toBeGreaterThan(1.05);
    expect(cw).toBeLessThan(1.2);
    // increasing until the plate hides behind the shaft (cos ψ ≤ a cos ψ₀), then constant
    let prev = 0;
    for (let u = 0; u <= 1; u += 0.01) {
      const a = throttleEffectiveArea(D, u, opts);
      expect(a).toBeGreaterThanOrEqual(prev);
      if (u < 0.85) expect(a).toBeGreaterThan(prev);
      prev = a;
    }
    expect(throttleEffectiveArea(D, 1, opts, 1)).toBeGreaterThan(throttleEffectiveArea(D, 1, opts));
  });

  it('Model T carburettor: chokes at idle (∂ṁ/∂p_manifold = 0), plausible WOT air flow, monotone in the opening', () => {
    expect(hasSeparateThrottle(MODEL_T.manifolds)).toBe(true);
    expect(hasSeparateThrottle(CFR_F1.manifolds)).toBe(false);
    expect(() => new CarburettorFlowModel(CFR_F1.manifolds, { venturiDischargeCoefficient: 0.6 })).toThrow();
    const carb = new CarburettorFlowModel(MODEL_T.manifolds, { venturiDischargeCoefficient: 0.6 });
    const out = newSeriesOrificeFlow();
    const amb = [101325, 288.7, 287.05, 1.4] as const;
    // idle: closed plate, manifold at ≈ 0.3 bar
    carb.setOpening(0);
    const mi = carb.flow(...amb, 30000, 300, 287.05, 1.4, out);
    expect(mi).toBeGreaterThan(0);
    expect(out.chokedB).toBe(true);
    expect(out.chokedA).toBe(false);
    expect(out.dmdotdpb).toBe(0);
    expect(out.intermediatePressure / amb[0]).toBeGreaterThan(0.999); // the venturi sees almost no flow
    // same call through the generic solver
    const ref = newSeriesOrificeFlow();
    expect(seriesOrificeFlow(carb.venturiCdA, carb.throttleCdA, ...amb, 30000, 300, 287.05, 1.4, ref)).toBe(mi);
    // WOT at 1600 rpm: Ford's 20 bhp needs ≈ 0.018–0.022 m³/s of air (induction research). With the CFR
    // venturi C_D (0.6) the 23/32 in venturi + 1 in plate pass 0.02 m³/s at ≈ 12 kPa manifold depression
    // (a period "strangling tube" sized for velocity; the Model T venturi C_D is a calibration item).
    carb.setOpening(1);
    const rhoAmb = amb[0] / (amb[2] * amb[1]);
    let lo = 1000;
    let hi = 50000;
    for (let k = 0; k < 60; k++) {
      const dp = 0.5 * (lo + hi);
      if (carb.flow(...amb, amb[0] - dp, 300, 287.05, 1.4, out) / rhoAmb > 0.02) hi = dp;
      else lo = dp;
    }
    expect(lo).toBeGreaterThan(5000);
    expect(lo).toBeLessThan(20000);
    const mw = carb.flow(...amb, 96000, 300, 287.05, 1.4, out);
    // the open plate costs a little relative to the venturi alone
    const o = newOrificeFlow();
    const mv = orificeFlow(carb.venturiCdA, ...amb, 96000, 300, 287.05, 1.4, o);
    expect(mw / mv).toBeGreaterThan(0.85);
    expect(mw / mv).toBeLessThan(1);
    // part throttle: monotone in the opening at a fixed manifold pressure
    let prev = 0;
    for (let u = 0; u <= 1.0001; u += 0.05) {
      carb.setOpening(u);
      const m = carb.flow(...amb, 60000, 300, 287.05, 1.4, out);
      expect(m).toBeGreaterThanOrEqual(prev);
      if (u < 0.85) expect(m).toBeGreaterThan(prev); // constant once the plate hides behind the shaft
      prev = m;
    }
    // backflow passes the plate first and carries the manifold gas
    carb.setOpening(0.5);
    expect(carb.flow(...amb, 103000, 330, 290, 1.36, out)).toBeLessThan(0);
    expect(out.upstream).toBe(1);
    // optional fixed restriction is folded into the venturi
    const carbR = new CarburettorFlowModel(MODEL_T.manifolds, { venturiDischargeCoefficient: 0.6, restrictionArea: 2e-4 });
    expect(carbR.venturiCdA).toBeCloseTo(seriesEffectiveArea(carb.venturiCdA, 2e-4), 18);
  });
});
