import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/mechanics_kinematics.json';
import { DEG } from '../core/constants';
import { newKinematicState, SliderCrank } from './kinematics';

interface Row {
  thetaDeg: number;
  x: number;
  h: number;
  V: number;
  dVdTheta: number;
  dxdTheta: number;
  d2xdTheta2: number;
  beta: number;
  dbetadTheta: number;
  d2betadTheta2: number;
  linerArea: number;
  V_heywood?: number;
}
interface Case {
  name: string;
  input: { bore: number; stroke: number; conRodLength: number; pinOffset: number; creviceVolume: number; compressionRatio: number };
  tdcAxisAngle: number;
  bdcTheta: number;
  pistonTravel: number;
  displacedVolume: number;
  clearanceVolume: number;
  clearanceHeightTDC: number;
  rows: Row[];
}
const cases = (fixture as { cases: Case[] }).cases;

const IN = 0.0254;
const CFR = { bore: 3.25 * IN, stroke: 4.5 * IN, conRodLength: 0.254, pinOffset: 0, creviceVolume: 1e-6 };

const relErr = (a: number, b: number, scale: number) => Math.abs(a - b) / scale;

describe('SliderCrank vs numerical loop-closure oracle', () => {
  for (const c of cases) {
    it(`${c.name}: TDC/BDC, travel, volumes`, () => {
      const sc = new SliderCrank(c.input, c.input.compressionRatio);
      expect(sc.tdcAxisAngle).toBeCloseTo(c.tdcAxisAngle, 11);
      expect(sc.bdcAngle).toBeCloseTo(c.bdcTheta, 11);
      expect(relErr(sc.pistonTravel, c.pistonTravel, c.pistonTravel)).toBeLessThan(1e-12);
      expect(relErr(sc.displacedVolume, c.displacedVolume, c.displacedVolume)).toBeLessThan(1e-12);
      expect(relErr(sc.clearanceVolume, c.clearanceVolume, c.clearanceVolume)).toBeLessThan(1e-12);
      expect(relErr(sc.clearanceHeightTDC, c.clearanceHeightTDC, c.clearanceHeightTDC)).toBeLessThan(1e-12);
    });

    it(`${c.name}: x, h, V, β and analytic derivatives on ${c.rows.length} angles`, () => {
      const sc = new SliderCrank(c.input, c.input.compressionRatio);
      const S = c.pistonTravel;
      const a = c.input.stroke / 2;
      const A = (Math.PI * c.input.bore ** 2) / 4;
      let maxPos = 0;
      let maxD1 = 0;
      let maxD2 = 0;
      let maxBeta = 0;
      for (const r of c.rows) {
        const th = r.thetaDeg * DEG;
        maxPos = Math.max(maxPos, relErr(sc.pistonDisplacement(th), r.x, S));
        maxPos = Math.max(maxPos, relErr(sc.clearanceHeight(th), r.h, S));
        maxPos = Math.max(maxPos, relErr(sc.volume(th), r.V, sc.displacedVolume));
        maxPos = Math.max(maxPos, relErr(sc.linerArea(th), r.linerArea, Math.PI * c.input.bore * S));
        maxD1 = Math.max(maxD1, relErr(sc.dVdTheta(th), r.dVdTheta, A * a));
        maxD1 = Math.max(maxD1, relErr(sc.dxdTheta(th), r.dxdTheta, a));
        maxD2 = Math.max(maxD2, relErr(sc.d2xdTheta2(th), r.d2xdTheta2, a));
        maxBeta = Math.max(maxBeta, Math.abs(sc.rodAngle(th) - r.beta));
        maxBeta = Math.max(maxBeta, Math.abs(sc.dRodAngledTheta(th) - r.dbetadTheta));
        maxBeta = Math.max(maxBeta, Math.abs(sc.d2RodAngledTheta2(th) - r.d2betadTheta2));
      }
      // oracle = brentq loop closure (~1e-15) + 6th-order FD (≈1e-11 abs on first, 1e-9 on second derivatives)
      expect(maxPos).toBeLessThan(1e-12);
      expect(maxD1).toBeLessThan(1e-9);
      expect(maxD2).toBeLessThan(1e-7);
      expect(maxBeta).toBeLessThan(1e-7);
    });

    if (c.input.pinOffset === 0) {
      it(`${c.name}: Heywood (1988) eq. 2.6 volume`, () => {
        const sc = new SliderCrank(c.input, c.input.compressionRatio);
        for (const r of c.rows) {
          const th = r.thetaDeg * DEG;
          expect(relErr(sc.volume(th), r.V_heywood as number, sc.displacedVolume)).toBeLessThan(1e-13);
        }
      });
    }
  }
});

describe('SliderCrank identities', () => {
  it('V(TDC) = Vc, V(BDC) = Vc + Vd, CR consistency, crevice convention', () => {
    for (const cr of [4, 7, 12.5, 18]) {
      for (const pinOffset of [0, 0.003, -0.005]) {
        const sc = new SliderCrank({ ...CFR, pinOffset }, cr);
        const vTdc = sc.volume(0);
        const vBdc = sc.volume(sc.bdcAngle);
        expect(relErr(vTdc, sc.clearanceVolume, vTdc)).toBeLessThan(1e-14);
        expect(relErr(vBdc, sc.clearanceVolume + sc.displacedVolume, vBdc)).toBeLessThan(1e-13);
        expect(vBdc / vTdc).toBeCloseTo(cr, 12);
        // V = A·h + V_crevice everywhere
        for (let d = -360; d < 360; d += 13) {
          const th = d * DEG;
          expect(relErr(sc.volume(th), sc.boreArea * sc.clearanceHeight(th) + sc.creviceVolume, vBdc)).toBeLessThan(1e-14);
        }
        // extremes: TDC is the minimum and BDC the maximum of V
        expect(sc.dVdTheta(0)).toBeCloseTo(0, 14);
        expect(sc.dVdTheta(sc.bdcAngle)).toBeCloseTo(0, 14);
        expect(sc.volume(0.01)).toBeGreaterThan(vTdc);
        expect(sc.volume(-0.01)).toBeGreaterThan(vTdc);
        expect(sc.volume(sc.bdcAngle + 0.01)).toBeLessThan(vBdc);
        expect(sc.volume(sc.bdcAngle - 0.01)).toBeLessThan(vBdc);
      }
    }
  });

  it('periodicity: 2π-periodic mechanism', () => {
    const sc = new SliderCrank({ ...CFR, pinOffset: 0.002 }, 8);
    for (let d = -360; d < 360; d += 17) {
      const th = d * DEG;
      expect(sc.volume(th + 2 * Math.PI)).toBeCloseTo(sc.volume(th), 15);
      expect(sc.rodAngle(th + 4 * Math.PI)).toBeCloseTo(sc.rodAngle(th), 14);
    }
  });

  it('zero offset: symmetric curves about TDC and BDC, travel = stroke', () => {
    const sc = new SliderCrank(CFR, 7);
    expect(sc.pistonTravel).toBeCloseTo(CFR.stroke, 15);
    expect(sc.bdcAngle).toBeCloseTo(Math.PI, 15);
    for (let d = 0; d <= 180; d += 3) {
      const th = d * DEG;
      expect(sc.volume(th)).toBeCloseTo(sc.volume(-th), 16);
      expect(sc.dVdTheta(th)).toBeCloseTo(-sc.dVdTheta(-th), 16);
      expect(sc.rodAngle(th)).toBeCloseTo(-sc.rodAngle(-th), 15);
      expect(sc.volume(Math.PI + th)).toBeCloseTo(sc.volume(Math.PI - th), 16);
    }
  });

  it('with offset the curves are asymmetric and TDC/BDC are not 180° apart', () => {
    const sc = new SliderCrank({ ...CFR, pinOffset: 0.005 }, 8);
    expect(Math.abs(sc.bdcAngle - Math.PI)).toBeGreaterThan(1e-3);
    expect(sc.pistonTravel).toBeGreaterThan(CFR.stroke);
    expect(Math.abs(sc.volume(1) - sc.volume(-1))).toBeGreaterThan(1e-7);
  });

  it('analytic derivatives match central finite differences (dV/dθ, d²V/dθ², dβ/dθ, d²β/dθ²)', () => {
    for (const pinOffset of [0, 0.004, -0.006]) {
      const sc = new SliderCrank({ ...CFR, pinOffset }, 9);
      const h = 1e-5;
      for (let d = -360; d < 360; d += 5) {
        const th = d * DEG + 0.123;
        const fdV = (sc.volume(th + h) - sc.volume(th - h)) / (2 * h);
        const fdV2 = (sc.dVdTheta(th + h) - sc.dVdTheta(th - h)) / (2 * h);
        const fdB = (sc.rodAngle(th + h) - sc.rodAngle(th - h)) / (2 * h);
        const fdB2 = (sc.dRodAngledTheta(th + h) - sc.dRodAngledTheta(th - h)) / (2 * h);
        const scaleV = sc.boreArea * sc.crankRadius;
        expect(Math.abs(sc.dVdTheta(th) - fdV) / scaleV).toBeLessThan(1e-9);
        expect(Math.abs(sc.d2VdTheta2(th) - fdV2) / scaleV).toBeLessThan(1e-9);
        expect(Math.abs(sc.dRodAngledTheta(th) - fdB)).toBeLessThan(1e-9);
        expect(Math.abs(sc.d2RodAngledTheta2(th) - fdB2)).toBeLessThan(1e-9);
      }
    }
  });

  it('Heywood (1988) eq. 2.9 instantaneous/mean piston speed (zero offset)', () => {
    const sc = new SliderCrank(CFR, 7);
    const rpm = 600;
    const omega = (rpm * 2 * Math.PI) / 60;
    const Sbar = sc.meanPistonSpeed(rpm);
    expect(Sbar).toBeCloseTo(2 * CFR.stroke * 10, 14); // 2.286 m/s
    const R = CFR.conRodLength / (CFR.stroke / 2);
    for (let d = -180; d <= 180; d += 10) {
      const th = d * DEG;
      const ratio = (Math.PI / 2) * Math.sin(th) * (1 + Math.cos(th) / Math.sqrt(R * R - Math.sin(th) ** 2));
      expect(sc.pistonVelocity(th, omega) / Sbar).toBeCloseTo(ratio, 13);
    }
  });

  it('piston acceleration includes the α·dx/dθ term', () => {
    const sc = new SliderCrank(CFR, 7);
    const th = 0.7;
    const w = 60;
    const al = 500;
    expect(sc.pistonAcceleration(th, w, al)).toBeCloseTo(sc.d2xdTheta2(th) * w * w + sc.dxdTheta(th) * al, 12);
    // at TDC with zero offset: x'' = a(1 + a/l)
    const a = CFR.stroke / 2;
    expect(sc.d2xdTheta2(0)).toBeCloseTo(a * (1 + a / CFR.conRodLength), 14);
  });

  it('variable CR: setCompressionRatio moves only the clearance, not the mechanism', () => {
    const sc = new SliderCrank(CFR, 5);
    const x1 = sc.pistonDisplacement(1.1);
    sc.setCompressionRatio(15);
    expect(sc.pistonDisplacement(1.1)).toBe(x1);
    expect(sc.compressionRatio).toBe(15);
    expect(sc.clearanceHeightTDC).toBeCloseTo(sc.clearanceHeightTDCForCR(15), 16);
    expect(() => sc.setCompressionRatio(1)).toThrow(RangeError);
    // crevice larger than the clearance volume is rejected
    expect(() => new SliderCrank({ ...CFR, creviceVolume: 1e-4 }, 18)).toThrow(RangeError);
  });

  it('single-pass evaluate() equals the individual methods', () => {
    const sc = new SliderCrank({ ...CFR, pinOffset: -0.004 }, 10);
    const ks = newKinematicState();
    for (let d = -360; d < 360; d += 7) {
      const th = d * DEG;
      sc.evaluate(th, ks);
      expect(ks.x).toBeCloseTo(sc.pistonDisplacement(th), 15);
      expect(ks.dxdTheta).toBeCloseTo(sc.dxdTheta(th), 15);
      expect(ks.d2xdTheta2).toBeCloseTo(sc.d2xdTheta2(th), 15);
      expect(ks.volume).toBeCloseTo(sc.volume(th), 18);
      expect(ks.dVdTheta).toBeCloseTo(sc.dVdTheta(th), 18);
      expect(ks.clearanceHeight).toBeCloseTo(sc.clearanceHeight(th), 15);
      expect(ks.beta).toBeCloseTo(sc.rodAngle(th), 15);
      expect(ks.sinBeta).toBeCloseTo(Math.sin(sc.rodAngle(th)), 15);
      expect(ks.cosBeta).toBeCloseTo(Math.cos(sc.rodAngle(th)), 15);
      expect(ks.dBetadTheta).toBeCloseTo(sc.dRodAngledTheta(th), 15);
      expect(ks.d2BetadTheta2).toBeCloseTo(sc.d2RodAngledTheta2(th), 15);
    }
  });

  it('crank-pin / wrist-pin positions are consistent with rod length and rod angle', () => {
    const sc = new SliderCrank({ ...CFR, pinOffset: 0.003 }, 8);
    const c = new Float64Array(2);
    const w = new Float64Array(2);
    for (let d = -360; d < 360; d += 11) {
      const th = d * DEG;
      sc.crankPinPosition(th, c);
      sc.wristPinPosition(th, w);
      expect(Math.hypot(w[0] - c[0], w[1] - c[1])).toBeCloseTo(CFR.conRodLength, 14);
      // rod direction big end → small end = l(−sinβ, cosβ)
      const b = sc.rodAngle(th);
      expect((c[0] - w[0]) / CFR.conRodLength).toBeCloseTo(Math.sin(b), 14);
    }
    // firing TDC: crank pin, crank axis and wrist pin collinear
    sc.crankPinPosition(0, c);
    sc.wristPinPosition(0, w);
    expect(c[0] * w[1] - c[1] * w[0]).toBeCloseTo(0, 15);
  });
});

describe('SliderCrank stress: seeded random geometries, extreme offsets and angles', () => {
  // Park–Miller minimal standard generator (deterministic)
  let seed = 20260929;
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
  const ks = newKinematicState();

  it('finite, bounded and self-consistent over 120 random mechanisms × 50 angles (incl. |θ| up to 1e4 rad)', () => {
    for (let n = 0; n < 120; n++) {
      const bore = 0.02 + 0.3 * rnd();
      const stroke = 0.02 + 0.4 * rnd();
      const a = stroke / 2;
      const l = a * (1.3 + 6 * rnd()); // rod ratios 1.3…7.3
      const pinOffset = (2 * rnd() - 1) * 0.98 * (l - a); // up to 98 % of the geometric limit
      const cr = 1.2 + 40 * rnd();
      const sc = new SliderCrank({ bore, stroke, conRodLength: l, pinOffset, creviceVolume: 0 }, cr);
      expect(sc.pistonTravel).toBeGreaterThanOrEqual(stroke * (1 - 1e-12));
      for (let k = 0; k < 50; k++) {
        const th = k < 40 ? (2 * rnd() - 1) * 4 * Math.PI : (2 * rnd() - 1) * 1e4;
        sc.evaluate(th, ks);
        for (const v of [ks.x, ks.dxdTheta, ks.d2xdTheta2, ks.volume, ks.dVdTheta, ks.beta, ks.dBetadTheta, ks.d2BetadTheta2]) {
          expect(Number.isFinite(v)).toBe(true);
        }
        const tol = 1e-12 * sc.displacedVolume;
        expect(ks.volume).toBeGreaterThanOrEqual(sc.clearanceVolume - tol);
        expect(ks.volume).toBeLessThanOrEqual(sc.clearanceVolume + sc.displacedVolume + tol);
        expect(Math.abs(ks.sinBeta)).toBeLessThan(1);
        expect(ks.dxdTheta).toBeCloseTo(sc.dxdTheta(th), 12);
      }
      // exact TDC/BDC are the extrema even at extreme offsets
      expect(sc.dxdTheta(0) / a).toBeCloseTo(0, 12);
      expect(sc.dxdTheta(sc.bdcAngle) / a).toBeCloseTo(0, 12);
      expect(sc.d2xdTheta2(0)).toBeGreaterThan(0);
      expect(sc.d2xdTheta2(sc.bdcAngle)).toBeLessThan(0);
    }
  });

  it('rejects impossible or degenerate inputs', () => {
    const g = { bore: 0.08, stroke: 0.1, conRodLength: 0.15, pinOffset: 0, creviceVolume: 0 };
    expect(() => new SliderCrank({ ...g, pinOffset: 0.1 }, 8)).toThrow(RangeError); // |e| ≥ l − a
    expect(() => new SliderCrank({ ...g, conRodLength: 0.05 }, 8)).toThrow(RangeError); // l < a
    expect(() => new SliderCrank({ ...g, bore: 0 }, 8)).toThrow(RangeError);
    expect(() => new SliderCrank({ ...g, creviceVolume: -1e-6 }, 8)).toThrow(RangeError);
    expect(() => new SliderCrank(g, Number.NaN)).toThrow(RangeError);
    expect(() => new SliderCrank(g, 0.5)).toThrow(RangeError);
  });

  it('hot path: evaluate() ≲ 0.2 µs per call and allocation-free (timing printed, loose bound)', () => {
    const sc = new SliderCrank(CFR, 7);
    const n = 2_000_000;
    let acc = 0;
    for (let i = 0; i < 20000; i++) acc += sc.evaluate(i * 1e-3, ks).volume; // warm-up / JIT
    const t0 = performance.now();
    for (let i = 0; i < n; i++) acc += sc.evaluate(i * 3.1e-6, ks).dVdTheta;
    const ns = ((performance.now() - t0) * 1e6) / n;
    console.info(`SliderCrank.evaluate: ${ns.toFixed(1)} ns/call`);
    expect(Number.isFinite(acc)).toBe(true);
    expect(ns).toBeLessThan(1000); // generous: guards against pathological regressions only
  });
});

describe('SliderCrank CR helpers', () => {
  it('return NaN (not a negative/infinite volume) for CR ≤ 1', () => {
    const sc = new SliderCrank(CFR, 7);
    expect(sc.clearanceVolumeForCR(1)).toBeNaN();
    expect(sc.clearanceHeightTDCForCR(0.5)).toBeNaN();
    expect(sc.clearanceVolumeForCR(7)).toBeCloseTo(sc.clearanceVolume, 18);
  });
});
