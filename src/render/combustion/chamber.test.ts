import { describe, expect, it } from 'vitest';
import { MODEL_T, MODEL_T_POCKET_PLAN_AREA } from '../../physics/engines/model-t';
import {
  boundingSize,
  chamberDepth,
  chamberShapeOf,
  chamberVolume,
  effectiveDepth,
  EXIT_BORE,
  EXIT_POCKET,
  EXIT_TRANSFER,
  farthestFootprintPoint,
  insideBurnedShape,
  insideChamberShape,
  insideRoundRect,
  projectIntoChamber,
  rayChamberIntervals,
  rayRoundRect,
  type ChamberShape,
} from './chamber';
import { insideChamber, rayChamberInterval, type V3 } from './geometry';
import { sharpLHeadSpec, testSpec, TEST_R } from './test-utils';

const IN = 0.0254;
const T = chamberShapeOf(MODEL_T);
const H_TDC = T.depthTDC;
const H_MID = H_TDC + 0.05;
const H_BDC = H_TDC + MODEL_T.geometry.stroke;

function lcg(seed: number): () => number {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296);
}

describe('Model T L-head shape from the spec', () => {
  it('crown 5/16 in above the deck at TDC; bore-column depth 1 in (the spec frame convention)', () => {
    expect(T.kind).toBe('l-head');
    expect(H_TDC).toBeCloseTo(1.0 * IN, 12);
    expect(-H_TDC - T.deckY).toBeCloseTo((5 / 16) * IN, 12); // crown y − deck y
    expect(chamberDepth(T, 123, 0)).toBeCloseTo(H_TDC, 12); // independent of clearanceHeight
    expect(chamberDepth(T, 123, 0.03)).toBeCloseTo(H_TDC + 0.03, 12);
    const disc = chamberShapeOf(testSpec());
    expect(disc.kind).toBe('flat-disc');
    expect(chamberDepth(disc, 0.0213, 0.004)).toBe(0.0213);
    expect(effectiveDepth(disc, 0.0213)).toBe(0.0213);
  });

  it('pocket ∖ bore plan area: scipy quad oracle; the spec’s 4000² grid value is 0.03 % low', () => {
    // scipy.integrate.quad of the exact z-slice length of (rounded rectangle ∖ disc) over x, split at
    // xMin, xMin + rc, xMax − rc, xMax, −R (epsrel 1e-12): 4.603020157589311e-3 m².
    expect(T.pocketPlanArea).toBeCloseTo(4.603020157589311e-3, 10); // 1e-8 relative
    expect(Math.abs(T.pocketPlanArea - MODEL_T_POCKET_PLAN_AREA) / MODEL_T_POCKET_PLAN_AREA).toBeLessThan(5e-4);
    // centroid on the valve side (−x), on the axis (pocket symmetric in z)
    expect(T.pocketCentroid[0]).toBeLessThan(-T.R);
    expect(Math.abs(T.pocketCentroid[1])).toBeLessThan(1e-9);
    expect(T.axis[0]).toBeCloseTo(-1, 12);
    expect(T.s0).toBeCloseTo(-T.R, 12);
    expect(T.s1).toBeCloseTo(0.0915, 12);
    expect(T.width).toBeCloseTo(0.095, 12);
  });

  it('conservation: bore column + pocket + crevice at TDC = the clearance volume of the spec CR', () => {
    const g = MODEL_T.geometry;
    const Vd = (Math.PI / 4) * g.bore ** 2 * g.stroke;
    const Vc = Vd / (g.compressionRatio - 1);
    expect(Math.abs(chamberVolume(T, H_TDC) + g.creviceVolume - Vc) / Vc).toBeLessThan(1e-4);
  });

  it('Monte Carlo volume of insideChamberShape = chamberVolume (TDC and mid-stroke)', () => {
    for (const h of [H_TDC, H_MID]) {
      const rnd = lcg(11);
      const x0 = -0.0915, x1 = T.R, z0 = -T.R, z1 = T.R, y0 = Math.min(-h, T.deckY), y1 = 0;
      const box = (x1 - x0) * (z1 - z0) * (y1 - y0);
      let hit = 0;
      const N = 400000;
      for (let i = 0; i < N; i++) {
        if (insideChamberShape(x0 + (x1 - x0) * rnd(), y0 + (y1 - y0) * rnd(), z0 + (z1 - z0) * rnd(), T, h)) hit++;
      }
      const V = (hit / N) * box;
      expect(Math.abs(V - chamberVolume(T, h)) / chamberVolume(T, h)).toBeLessThan(0.01);
    }
  });

  it('inside test: bore column, pocket (also below the crown), piston top land, head metal', () => {
    const yMid = 0.5 * (T.deckY + T.roofY);
    expect(insideChamberShape(-0.07, yMid, 0, T, H_TDC)).toBe(true); // pocket over the valves
    expect(insideChamberShape(-0.07, T.roofY + 1e-4, 0, T, H_TDC)).toBe(false); // head over the pocket
    expect(insideChamberShape(0, -0.005, 0, T, H_TDC)).toBe(true); // head cavity over the bore
    expect(insideChamberShape(0, -0.005, 0, T, H_TDC)).toBe(true);
    expect(insideChamberShape(-0.07, -0.005, 0, T, H_TDC)).toBe(false); // above the pocket roof, outside the bore
    // between the deck and the raised crown: pocket outside the bore = gas, inside the bore = piston
    const yLow = T.deckY + 0.003;
    expect(yLow).toBeLessThan(-H_TDC);
    expect(insideChamberShape(-0.05, yLow, 0, T, H_TDC)).toBe(true);
    expect(insideChamberShape(-0.045, yLow, 0, T, H_TDC)).toBe(false);
    // pocket corners are rounded
    expect(insideRoundRect(-0.0915 + 1e-3, 0.0475 - 1e-3, T)).toBe(false);
    expect(insideRoundRect(-0.0915 + 1e-3, 0, T)).toBe(true);
    // flame sphere ∩ chamber: the part of the ball inside the piston is not burned gas
    const c: V3 = [-0.06, T.roofY - 0.004, 0];
    expect(insideBurnedShape(-0.05, T.roofY - 0.004, 0, T, H_TDC, c, 0.02)).toBe(true);
    expect(insideBurnedShape(-0.045, yLow, 0, T, H_TDC, c, 0.03)).toBe(false);
  });

  it('bounding size covers every chamber point', () => {
    for (const h of [H_TDC, H_BDC]) {
      const B = boundingSize(T, h);
      expect(B).toBeGreaterThanOrEqual(Math.hypot(-0.0915, 0.0275) + 0);
      expect(B).toBeGreaterThanOrEqual(h);
    }
  });
});

describe('ray ∩ chamber shape (mirrored in the GLSL)', () => {
  /** Brute-force inside runs along the ray (t ≥ 0). */
  function brute(ro: number[], rd: number[], s: ChamberShape, h: number, tMax: number, dt: number): [number, number][] {
    const runs: [number, number][] = [];
    let start = -1, last = -1;
    for (let t = 0; t <= tMax; t += dt) {
      const ins = insideChamberShape(ro[0] + rd[0] * t, ro[1] + rd[1] * t, ro[2] + rd[2] * t, s, h);
      if (ins) {
        if (start < 0) start = t;
        last = t;
      } else if (start >= 0) {
        runs.push([start, last]);
        start = -1;
      }
    }
    if (start >= 0) runs.push([start, last]);
    return runs;
  }

  it('the disc case equals rayChamberInterval', () => {
    const disc = chamberShapeOf(testSpec());
    const rnd = lcg(5);
    const out = new Float64Array(9);
    for (let n = 0; n < 500; n++) {
      const ro = [(rnd() - 0.5) * 0.3, (rnd() - 0.5) * 0.3, (rnd() - 0.5) * 0.3];
      const rd = [rnd() - 0.5, rnd() - 0.5, rnd() - 0.5];
      const ref = rayChamberInterval(ro, rd, TEST_R, 0.018);
      const k = rayChamberIntervals(ro, rd, disc, 0.018, out);
      if (!ref) {
        expect(k).toBe(0);
      } else {
        expect(k).toBe(1);
        expect(out[0]).toBe(ref[0]);
        expect(out[1]).toBe(ref[1]);
        expect(out[2]).toBe(EXIT_BORE);
      }
    }
  });

  it('L-head intervals agree with brute-force sampling (TDC with the crown above the deck, mid-stroke, BDC)', () => {
    const out = new Float64Array(9);
    let two = 0, kinds = new Set<number>();
    for (const h of [H_TDC, H_MID, H_BDC]) {
      const rnd = lcg(Math.round(h * 1e5));
      for (let n = 0; n < 700; n++) {
        // rays aimed through a random chamber-ish point so most of them hit
        const tgt = [-0.09 + 0.14 * rnd(), -h * rnd(), -0.05 + 0.1 * rnd()];
        const ro = [(rnd() - 0.5) * 0.4, (rnd() - 0.5) * 0.4, (rnd() - 0.5) * 0.4];
        let rd = [tgt[0] - ro[0], tgt[1] - ro[1], tgt[2] - ro[2]];
        const L = Math.hypot(rd[0], rd[1], rd[2]);
        rd = rd.map((v) => v / L);
        const k = rayChamberIntervals(ro, rd, T, h, out);
        expect(k).toBeLessThanOrEqual(2);
        const runs = brute(ro, rd, T, h, 0.7, 5e-5);
        // every analytic interval longer than a sampling step matches one sampled run and vice versa
        // (grazing slivers below the sampling resolution are ignored on both sides)
        for (let i = 0; i < k; i++) {
          const a = out[3 * i], b = out[3 * i + 1];
          if (b - a < 4e-4) continue;
          const m = runs.find(([ra, rb]) => Math.abs(ra - a) < 1.5e-4 && Math.abs(rb - b) < 1.5e-4);
          expect(m).toBeDefined();
        }
        for (const [ra, rb] of runs) {
          if (rb - ra < 4e-4) continue;
          let found = false;
          for (let i = 0; i < k; i++) found ||= Math.abs(out[3 * i] - ra) < 1.5e-4 && Math.abs(out[3 * i + 1] - rb) < 1.5e-4;
          expect(found).toBe(true);
        }
        if (k === 2) two++;
        for (let i = 0; i < k; i++) kinds.add(out[3 * i + 2]);
      }
    }
    expect(two).toBeGreaterThan(5); // the union really produces disjoint intervals
    expect([...kinds].sort()).toEqual([EXIT_BORE, EXIT_POCKET, EXIT_TRANSFER]);
  });

  it('pocket → piston top land → pocket: two intervals, the first leaving into the piston (transfer exit)', () => {
    const out = new Float64Array(9);
    const y = T.deckY + 0.003; // below the raised crown at TDC
    const k = rayChamberIntervals([-0.04, y, -0.2], [0, 0, 1], T, H_TDC, out);
    expect(k).toBe(2);
    const q = Math.sqrt(T.R ** 2 - 0.04 ** 2);
    expect(out[1]).toBeCloseTo(0.2 - q, 12);
    expect(out[2]).toBe(EXIT_TRANSFER);
    expect(out[3]).toBeCloseTo(0.2 + q, 12);
    expect(out[5]).toBe(EXIT_POCKET);
    // the same ray above the crown: one interval through pocket, bore column and pocket
    const k2 = rayChamberIntervals([-0.04, -H_TDC + 0.002, -0.2], [0, 0, 1], T, H_TDC, out);
    expect(k2).toBe(1);
    expect(out[2]).toBe(EXIT_POCKET);
    // camera inside the chamber starts at t = 0
    const k3 = rayChamberIntervals([0, -0.01, 0], [0, -1, 0], T, H_TDC, out);
    expect(k3).toBe(1);
    expect(out[0]).toBe(0);
    expect(out[1]).toBeCloseTo(H_TDC - 0.01, 12);
  });

  it('2-D rounded-rectangle chord agrees with sampling', () => {
    const rnd = lcg(9);
    const o = [0, 0];
    for (let n = 0; n < 300; n++) {
      const ox = -0.2 + 0.25 * rnd(), oz = -0.15 + 0.3 * rnd();
      const a = 2 * Math.PI * rnd();
      const dx = Math.cos(a), dz = Math.sin(a);
      rayRoundRect(ox, oz, dx, dz, T, o);
      let first = NaN, last = NaN;
      for (let t = -0.4; t <= 0.4; t += 2e-5) {
        if (insideRoundRect(ox + dx * t, oz + dz * t, T)) { if (Number.isNaN(first)) first = t; last = t; }
      }
      if (Number.isNaN(first)) expect(o[1] - o[0] < 1e-4).toBe(true);
      else {
        expect(Math.abs(o[0] - first)).toBeLessThan(5e-5);
        expect(Math.abs(o[1] - last)).toBeLessThan(5e-5);
      }
    }
  });
});

describe('proxy ownership (TS mirror of the GLSL fragment rule)', () => {
  /**
   * The gas volume draws three proxies with back faces / front faces: the bore cylinder (0.985 R, y in
   * [−h+g, −g]), the pocket prism (plan inset by g, y in [deck+g, roof−g]) and the transfer arc strip
   * (radius R+g, front faces = inward crossings). Each produces ≤ 1 fragment per ray; a fragment
   * survives iff its proxy owns the exit of the last visible interval. Check that the owner's fragment
   * exists, sits in the gas just before the exit (so the depth test sees the walls behind the gas).
   */
  const g = 2e-4;
  const Rb = 0.985 * T.R;
  const inset: ChamberShape = { ...T, rc: T.rc - g };
  const tmp = [0, 0];

  function exitOfConvex(ro: number[], rd: number[], kind: number, h: number): number {
    const slab = (y0: number, y1: number): [number, number] => {
      if (Math.abs(rd[1]) < 1e-15) return ro[1] >= y0 && ro[1] <= y1 ? [-1e30, 1e30] : [1e30, -1e30];
      const a = (y0 - ro[1]) / rd[1], b = (y1 - ro[1]) / rd[1];
      return [Math.min(a, b), Math.max(a, b)];
    };
    if (kind === EXIT_BORE) {
      const a = rd[0] ** 2 + rd[2] ** 2, b = 2 * (ro[0] * rd[0] + ro[2] * rd[2]), c = ro[0] ** 2 + ro[2] ** 2 - Rb * Rb;
      const disc = b * b - 4 * a * c;
      if (disc < 0) return NaN;
      const s = slab(-h + g, -g);
      const t0 = Math.max((-b - Math.sqrt(disc)) / (2 * a), s[0]), t1 = Math.min((-b + Math.sqrt(disc)) / (2 * a), s[1]);
      return t1 > t0 && t1 > 0 ? t1 : NaN;
    }
    if (kind === EXIT_POCKET) {
      rayRoundRect(ro[0], ro[2], rd[0], rd[2], inset, tmp);
      const s = slab(T.deckY + g, T.roofY - g);
      const t0 = Math.max(tmp[0], s[0]), t1 = Math.min(tmp[1], s[1]);
      return t1 > t0 && t1 > 0 ? t1 : NaN;
    }
    // transfer strip: inward crossing of the circle R+g within the strip's height and arc (+2° margin)
    const Rs = T.R + g;
    const a = rd[0] ** 2 + rd[2] ** 2, b = 2 * (ro[0] * rd[0] + ro[2] * rd[2]), c = ro[0] ** 2 + ro[2] ** 2 - Rs * Rs;
    const disc = b * b - 4 * a * c;
    if (disc < 0 || c < 0) return NaN;
    const t = (-b - Math.sqrt(disc)) / (2 * a);
    const p = [ro[0] + rd[0] * t, ro[1] + rd[1] * t, ro[2] + rd[2] * t];
    const ang = Math.atan2(p[2], p[0]);
    const [a0, a1] = T.transferArc;
    let d = ang - 0.5 * (a0 + a1);
    d = Math.atan2(Math.sin(d), Math.cos(d));
    const half = 0.5 * (a1 - a0) + (2 * Math.PI) / 180;
    return t > 0 && p[1] >= T.deckY + g && p[1] <= T.roofY - g && Math.abs(d) <= half ? t : NaN;
  }

  it('exactly one proxy owns each ray’s last visible exit, and its fragment lies in the gas', () => {
    const out = new Float64Array(9);
    let rays = 0, ok = 0, inGas = 0;
    const byKind = [0, 0, 0], hitByKind = [0, 0, 0];
    for (const h of [H_TDC, H_MID]) {
      const rnd = lcg(Math.round(h * 1e5) + 1);
      for (let n = 0; n < 4000; n++) {
        const tgt = [-0.09 + 0.14 * rnd(), -h * rnd(), -0.05 + 0.1 * rnd()];
        const ro = [(rnd() - 0.5) * 0.5, (rnd() - 0.2) * 0.4, (rnd() - 0.5) * 0.5];
        let rd = [tgt[0] - ro[0], tgt[1] - ro[1], tgt[2] - ro[2]];
        const L = Math.hypot(rd[0], rd[1], rd[2]);
        rd = rd.map((v) => v / L);
        const k = rayChamberIntervals(ro, rd, T, h, out);
        if (k === 0 || out[1] - out[0] < 1e-3) continue; // grazing slivers: a missing pixel at worst
        rays++;
        // default rule (no cut region): the last visible interval is the first one
        const tExit = out[1], kind = out[2];
        byKind[kind]++;
        const tf = exitOfConvex(ro, rd, kind, h);
        if (Number.isNaN(tf)) continue;
        ok++;
        hitByKind[kind]++;
        expect(tf).toBeLessThanOrEqual(tExit + 1e-9);
        if (tExit - tf < 0.01) {
          const p = [ro[0] + rd[0] * tf, ro[1] + rd[1] * tf, ro[2] + rd[2] * tf];
          if (insideChamberShape(p[0], p[1], p[2], T, h)) inGas++;
        }
      }
    }
    expect(rays).toBeGreaterThan(3000);
    expect(byKind[EXIT_BORE]).toBeGreaterThan(100);
    expect(byKind[EXIT_POCKET]).toBeGreaterThan(100);
    expect(byKind[EXIT_TRANSFER]).toBeGreaterThan(5);
    // The bore proxy is the CFR's 0.985 R cylinder: rays that only graze the outer 0.7 mm ring near
    // the liner/crown miss it (a thin silhouette, the same for the flat disc). The new proxies:
    expect(hitByKind[EXIT_BORE] / byKind[EXIT_BORE]).toBeGreaterThan(0.98);
    expect(hitByKind[EXIT_POCKET] / byKind[EXIT_POCKET]).toBeGreaterThan(0.995);
    expect(hitByKind[EXIT_TRANSFER] / byKind[EXIT_TRANSFER]).toBeGreaterThan(0.95);
    expect(inGas / ok).toBeGreaterThan(0.995);
  });
});

describe('nearest interior point and footprint extremes', () => {
  it('disc: identical to the tracers’ radial + plane clamp', () => {
    const disc = chamberShapeOf(testSpec());
    const out: V3 = [0, 0, 0];
    const rnd = lcg(21);
    for (let n = 0; n < 200; n++) {
      const x = (rnd() - 0.5) * 0.15, y = (rnd() - 0.7) * 0.05, z = (rnd() - 0.5) * 0.15;
      projectIntoChamber(x, y, z, disc, 0.02, out);
      const lim = 0.995 * TEST_R, r = Math.hypot(x, z), f = r > lim ? lim / r : 1;
      expect(out[0]).toBeCloseTo(x * f, 14);
      expect(out[2]).toBeCloseTo(z * f, 14);
      expect(out[1]).toBeCloseTo(Math.min(Math.max(y, -0.02), 0), 14);
      expect(insideChamber(out[0], out[1], out[2], TEST_R, 0.02)).toBe(true);
    }
  });

  it('L-head: every projected point is in the chamber; interior points do not move', () => {
    const out: V3 = [0, 0, 0];
    const rnd = lcg(22);
    for (const h of [H_TDC, H_MID]) {
      for (let n = 0; n < 3000; n++) {
        const x = -0.11 + 0.17 * rnd(), y = -h - 0.01 + (h + 0.02) * rnd(), z = (rnd() - 0.5) * 0.12;
        const d2 = projectIntoChamber(x, y, z, T, h, out);
        expect(insideChamberShape(out[0], out[1], out[2], T, h)).toBe(true);
        if (d2 === 0) expect([out[0], out[1], out[2]]).toEqual([x, y, z]);
      }
    }
    // a point in the piston top land next to the pocket goes out to the pocket, not up onto the crown
    const yLow = T.deckY + 0.003;
    projectIntoChamber(-0.045, yLow, 0, T, H_TDC, out);
    expect(out[1]).toBe(yLow);
    expect(Math.hypot(out[0], out[2])).toBeGreaterThan(0.99 * T.R);
  });

  it('end-gas site: the footprint point farthest from the flame (brute force over the outline)', () => {
    const o: [number, number] = [0, 0];
    const disc = chamberShapeOf(testSpec());
    expect(farthestFootprintPoint(0, TEST_R - 0.004, disc, o)).toBeCloseTo(2 * TEST_R - 0.004, 12);
    expect(o[0]).toBeCloseTo(0, 12);
    expect(o[1]).toBeCloseTo(-TEST_R, 12);
    const rnd = lcg(3);
    for (let n = 0; n < 50; n++) {
      const cx = -0.09 + 0.13 * rnd(), cz = (rnd() - 0.5) * 0.08;
      const d = farthestFootprintPoint(cx, cz, T, o);
      let best = 0;
      for (let i = 0; i < 20000; i++) {
        const a = (2 * Math.PI * i) / 20000;
        best = Math.max(best, Math.hypot(T.R * Math.cos(a) - cx, T.R * Math.sin(a) - cz));
        for (let k = 0; k < 4; k++) {
          const qx = (k & 1 ? T.ix1 : T.ix0) + T.rc * Math.cos(a), qz = (k & 2 ? T.iz1 : T.iz0) + T.rc * Math.sin(a);
          best = Math.max(best, Math.hypot(qx - cx, qz - cz));
        }
      }
      expect(d).toBeCloseTo(best, 6);
      expect(Math.hypot(o[0] - cx, o[1] - cz)).toBeCloseTo(d, 9);
    }
    // the stock plug over the valves: the end gas is at the far side of the bore
    farthestFootprintPoint(MODEL_T.sparkPlug.gapCenter[0], MODEL_T.sparkPlug.gapCenter[2], T, o);
    expect(o[0]).toBeCloseTo(T.R, 9);
  });

  it('transfer arc: the bore-circle arc inside the pocket plan, ends on the pocket outline', () => {
    const [a0, a1] = T.transferArc;
    expect(a1).toBeGreaterThan(a0);
    expect(Math.cos(0.5 * (a0 + a1))).toBeCloseTo(-1, 9); // centred on the valve side
    for (const a of [a0, a1]) {
      const inn = 1e-6 * Math.sign(0.5 * (a0 + a1) - a);
      expect(insideRoundRect(T.R * Math.cos(a + inn), T.R * Math.sin(a + inn), T)).toBe(true);
      expect(insideRoundRect(T.R * Math.cos(a - inn), T.R * Math.sin(a - inn), T)).toBe(false);
    }
    // a sharp-cornered variant: the arc ends where the circle crosses the pocket's inner wall x = −0.03
    const S = chamberShapeOf(sharpLHeadSpec());
    expect(S.rc).toBe(0);
    expect(S.R * Math.cos(S.transferArc[0])).toBeCloseTo(-0.03, 9);
    expect(S.R * Math.cos(S.transferArc[1])).toBeCloseTo(-0.03, 9);
  });
});
