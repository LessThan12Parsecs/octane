import { describe, expect, it } from 'vitest';
import lheadFx from '../../../test/fixtures/combustion_geometry_lhead_mc.json';
import type { EngineSpec, ValveSpec } from '../core/engine-spec';
import { CFR_F1 } from '../engines/cfr';
import { MODEL_T } from '../engines/model-t';
import {
  N_WALL_SURFACES,
  WALL_BLOCK,
  WALL_EXHAUST_VALVE,
  WALL_HEAD,
  WALL_INTAKE_VALVE,
  WALL_LINER,
  WALL_PISTON,
  flatChamberAreas,
  newChamberAreas,
  newWallHeatResult,
  newWallSurfaceArray,
  newWallSurfaceHeatResult,
  wallHeatLoss,
  wallHeatLossSurfaces,
  wallHeatLossTwoZone,
  wallHeatLossTwoZoneSurfaces,
  wallSurfaceTemperatures,
} from '../heat-transfer/wall-heat';
import { SliderCrank } from '../mechanics/kinematics';
import {
  DiscChamber,
  LHeadChamber,
  chamberFixedVolume,
  createChamber,
  lHeadPlanMetrics,
  newChamberResult,
  sliderCrankGeometry,
  type ChamberResult,
  type CombustionChamber,
} from './chamber';
import { entrainmentRates, newEntrainmentInputs, newEntrainmentRates } from './entrainment';
import { rapidDistortionIntensity } from './turbulence';
import {
  FlameGeometry,
  arcInsideBore,
  arcOfBoreInside,
  defaultInnerCellsOutside,
  lensArea,
  newFlameGeometryResult,
} from './flame-geometry';

const PI = Math.PI;
type Spec = Pick<EngineSpec, 'geometry' | 'sparkPlug' | 'intakeValve' | 'exhaustValve'>;

function lcg(seed: number): () => number {
  let s = seed;
  return () => (s = (s * 16807) % 2147483647) / 2147483647;
}

// ---------------------------------------------------------------------------------------------
// fixture geometries → chamber specs
// ---------------------------------------------------------------------------------------------

type FxGeom = (typeof lheadFx.geometries)[number];
type FxCase = FxGeom['cases'][number];

function specOf(G: FxGeom): Spec {
  const valve = (base: ValveSpec, v: FxGeom['valves'][number]): ValveSpec => ({
    ...base,
    position: [v.position[0], v.position[1]],
    headDiameter: v.headDiameter,
    count: v.count,
  });
  return {
    geometry: {
      ...MODEL_T.geometry,
      bore: G.bore,
      lHead: { deckY: G.deckY, crownAboveDeckAtTDC: G.crownAboveDeckAtTDC, pocket: { ...G.pocket, roofY: G.roofY } },
    },
    sparkPlug: { ...MODEL_T.sparkPlug, gapCenter: [G.spark[0], G.spark[1], G.spark[2]] },
    intakeValve: valve(MODEL_T.intakeValve, G.valves[0]),
    exhaustValve: valve(MODEL_T.exhaustValve, G.valves[1]),
  };
}

const chambers = new Map<string, LHeadChamber>();
function chamberOf(G: FxGeom): LHeadChamber {
  let c = chambers.get(G.name);
  if (!c) {
    c = createChamber(specOf(G)) as LHeadChamber;
    chambers.set(G.name, c);
  }
  return c;
}

const SURF = [
  ['head', WALL_HEAD],
  ['piston', WALL_PISTON],
  ['liner', WALL_LINER],
  ['block', WALL_BLOCK],
  ['intakeValve', WALL_INTAKE_VALVE],
  ['exhaustValve', WALL_EXHAUST_VALVE],
] as const;

// ---------------------------------------------------------------------------------------------
// DiscChamber: bit-identical adapter for the CFR
// ---------------------------------------------------------------------------------------------

describe('DiscChamber (CFR F-1): bit-identical to FlameGeometry + flatChamberAreas + the cycle-model expressions', () => {
  const g = CFR_F1.geometry;
  const kin = new SliderCrank(g, g.compressionRatio);
  const hMax = (kin.pistonTravel + kin.clearanceHeightTDCForCR(g.compressionRatioRange[0])) * 1.02;
  const fg = new FlameGeometry(g.bore, CFR_F1.sparkPlug.gapCenter, { maxHeight: hMax });
  const ch = createChamber(CFR_F1, { maxHeight: hMax });
  const legacyAreas = newChamberAreas();
  const fr = newFlameGeometryResult();
  const res = newChamberResult();
  const areas = newWallSurfaceArray();

  it('is a DiscChamber; geometry queries equal the legacy expressions', () => {
    expect(ch).toBeInstanceOf(DiscChamber);
    expect(ch.kind).toBe('flat-disc');
    expect(ch.fixedVolume).toBe(0);
    expect(chamberFixedVolume(g)).toBe(0);
    const rnd = lcg(11);
    for (let i = 0; i < 500; i++) {
      const h = 0.0005 + rnd() * hMax;
      // cycle model: discScale = (kin.boreArea · h)/V, integralLengthScale(h, …), clampRadiusGuess bounds
      expect(ch.chamberVolume(h)).toBe(kin.boreArea * h);
      expect(ch.meanDepth(h)).toBe(h);
      expect(ch.maxRadius(h)).toBe(fg.maxRadius(h));
      const b = fg.headDistance;
      expect(ch.inscribedRadius(h)).toBe(h > b ? Math.min(b, h - b, fg.radius - fg.offset) : 0);
      flatChamberAreas(g.bore, h, CFR_F1.intakeValve, CFR_F1.exhaustValve, legacyAreas);
      ch.surfaceAreas(h, areas);
      expect(areas[WALL_HEAD]).toBe(legacyAreas.head);
      expect(areas[WALL_PISTON]).toBe(legacyAreas.piston);
      expect(areas[WALL_LINER]).toBe(legacyAreas.liner);
      expect(areas[WALL_INTAKE_VALVE]).toBe(legacyAreas.intakeValves);
      expect(areas[WALL_EXHAUST_VALVE]).toBe(legacyAreas.exhaustValves);
      expect(areas[WALL_BLOCK]).toBe(0);
    }
  });

  it('volume, front area, radiusForVolume, burned fractions and the crevice-mouth fraction (toBe, 3000 points)', () => {
    const rnd = lcg(5);
    for (let i = 0; i < 3000; i++) {
      const h = 0.0005 + rnd() * hMax;
      const r = rnd() * fg.maxRadius(h) * 1.02;
      fg.evaluate(r, h, fr);
      ch.evaluate(r, h, res);
      expect(res.volume).toBe(fr.volume);
      expect(res.frontArea).toBe(fr.frontArea);
      expect(ch.volume(r, h)).toBe(fr.volume);
      expect(ch.frontArea(r, h)).toBe(fr.frontArea);
      // wallHeatLossTwoZone's fractions and the cycle model's crevice f_b
      const a = flatChamberAreas(g.bore, h, CFR_F1.intakeValve, CFR_F1.exhaustValve, legacyAreas);
      const headDisc = a.head + a.intakeValves + a.exhaustValves;
      let fh = headDisc > 0 ? fr.wettedHead / headDisc : 0;
      fh = fh > 0 ? (fh < 1 ? fh : 1) : 0;
      let fp = a.piston > 0 ? fr.wettedPiston / a.piston : 0;
      fp = fp > 0 ? (fp < 1 ? fp : 1) : 0;
      let fb = 0;
      if (a.liner > 0) {
        fb = fr.wettedLiner / a.liner;
        fb = fb > 0 ? (fb < 1 ? fb : 1) : 0;
      }
      expect(res.burnedFraction[WALL_HEAD]).toBe(fh);
      expect(res.burnedFraction[WALL_INTAKE_VALVE]).toBe(fh);
      expect(res.burnedFraction[WALL_EXHAUST_VALVE]).toBe(fh);
      expect(res.burnedFraction[WALL_PISTON]).toBe(fp);
      expect(res.burnedFraction[WALL_LINER]).toBe(fb);
      expect(res.creviceBurnedFraction).toBe(fb);
      const V = rnd() * ch.chamberVolume(h);
      const guess = rnd() * fg.maxRadius(h);
      expect(ch.radiusForVolume(V, h, guess)).toBe(fg.radiusForVolume(V, h, guess));
    }
  });

  it('6-surface wall heat with the disc fractions = the 5-surface wallHeatLoss / wallHeatLossTwoZone (toBe)', () => {
    const walls = CFR_F1.walls;
    const Tw = wallSurfaceTemperatures(walls, newWallSurfaceArray());
    const q5 = newWallHeatResult();
    const q6 = newWallSurfaceHeatResult();
    const rnd = lcg(17);
    for (let i = 0; i < 2000; i++) {
      const h = 0.0005 + rnd() * 0.12;
      const r = rnd() * fg.maxRadius(h);
      const hc = 50 + 3000 * rnd();
      const Tu = 300 + 1200 * rnd();
      const Tb = 1200 + 1700 * rnd();
      const c = i % 3 === 0 ? 4.3e-9 : 0;
      fg.evaluate(r, h, fr);
      const a = flatChamberAreas(g.bore, h, CFR_F1.intakeValve, CFR_F1.exhaustValve, legacyAreas);
      ch.evaluate(r, h, res);
      ch.surfaceAreas(h, areas);
      wallHeatLossTwoZone(hc, Tu, Tb, a, fr, walls, q5, c);
      wallHeatLossTwoZoneSurfaces(hc, Tu, Tb, areas, res.burnedFraction, Tw, q6, c);
      expect(q6.total).toBe(q5.total);
      expect(q6.burned).toBe(q5.burned);
      expect(q6.unburned).toBe(q5.unburned);
      expect(q6.surface[WALL_HEAD]).toBe(q5.head);
      expect(q6.surface[WALL_PISTON]).toBe(q5.piston);
      expect(q6.surface[WALL_LINER]).toBe(q5.liner);
      expect(q6.surface[WALL_INTAKE_VALVE]).toBe(q5.intakeValves);
      expect(q6.surface[WALL_EXHAUST_VALVE]).toBe(q5.exhaustValves);
      wallHeatLoss(hc, Tu, a, walls, q5, c);
      wallHeatLossSurfaces(hc, Tu, areas, Tw, q6, c);
      expect(q6.total).toBe(q5.total);
      expect(q6.unburned).toBe(q5.unburned);
      expect(q6.burned).toBe(0);
    }
  });

  it('clone() and the factory cache share the tables and give identical numbers', () => {
    const t0 = performance.now();
    const c2 = createChamber(CFR_F1, { maxHeight: hMax }); // cache hit: no rebuild
    const c3 = c2.clone();
    expect(performance.now() - t0).toBeLessThan(100);
    const r1 = newChamberResult();
    const r2 = newChamberResult();
    const rnd = lcg(3);
    for (let i = 0; i < 300; i++) {
      const h = 0.001 + rnd() * 0.1;
      const r = rnd() * 0.12;
      ch.evaluate(r, h, r1);
      c3.evaluate(r, h, r2);
      expect(r2).toEqual(r1);
    }
    expect(() => new FlameGeometry(g.bore, [0, -0.003, 0.01], { shareTablesWith: fg })).toThrow();
  });
});

// ---------------------------------------------------------------------------------------------
// FlameGeometry with the spark outside the bore planform (d > R)
// ---------------------------------------------------------------------------------------------

describe('FlameGeometry, spark outside the bore (d > R)', () => {
  const B = MODEL_T.geometry.bore;
  const R = B / 2;

  it('lens helpers for d > R match the classical acos formulas and the chord identity', () => {
    const d = 0.06;
    for (const rho of [0.005, 0.0124, 0.013, 0.02, 0.05, 0.08, 0.1, 0.107, 0.11]) {
      const ca = (d * d + rho * rho - R * R) / (2 * d * rho);
      const cb = (d * d + R * R - rho * rho) / (2 * d * R);
      const k = (-d + rho + R) * (d + rho - R) * (d - rho + R) * (d + rho + R);
      const lens = rho <= d - R ? 0 : rho >= R + d ? PI * R * R : rho * rho * Math.acos(ca) + R * R * Math.acos(cb) - 0.5 * Math.sqrt(k);
      expect(lensArea(rho, R, d)).toBeCloseTo(lens, 14);
      const phi = rho <= d - R || rho >= R + d ? 0 : 2 * Math.acos(ca);
      const psi = rho <= d - R ? 0 : rho >= R + d ? 2 * PI : 2 * Math.acos(cb);
      expect(arcInsideBore(rho, R, d)).toBeCloseTo(phi, 12);
      expect(arcOfBoreInside(rho, R, d)).toBeCloseTo(psi, 12);
      if (rho > d - R && rho < d + R) {
        // dL/dρ = ρΦ (co-area) and ρ sin(Φ/2) = R sin(Ψ/2)
        const dr = 1e-7;
        const fd = (lensArea(rho + dr, R, d) - lensArea(rho - dr, R, d)) / (2 * dr);
        expect(Math.abs(fd / (rho * arcInsideBore(rho, R, d)) - 1)).toBeLessThan(1e-6);
        expect(Math.abs(rho * Math.sin(0.5 * arcInsideBore(rho, R, d)) - R * Math.sin(0.5 * arcOfBoreInside(rho, R, d)))).toBeLessThan(1e-14);
      }
    }
  });

  it('needs the opt-in; rejects a centre on the bore circle', () => {
    expect(() => new FlameGeometry(B, [-0.06, -0.02, 0])).toThrow();
    expect(() => new FlameGeometry(B, [-R, -0.02, 0], { allowSparkOutsideBore: true })).toThrow();
    expect(new FlameGeometry(B, [-0.06, -0.02, 0], { allowSparkOutsideBore: true, maxHeight: 0.13 }).outside).toBe(true);
    expect(defaultInnerCellsOutside(R, 0.06)).toBe(64);
    expect(defaultInnerCellsOutside(R, R + 0.001)).toBeGreaterThan(64);
  });

  for (const c of [
    [-0.06, -0.0244, 0],
    [-0.05, -0.01, 0.012],
    [-0.08, -0.02, 0.03],
  ] as [number, number, number][]) {
    it(`fast table vs exact quadrature and dV/dr = A_f (centre ${c.map((v) => v * 1e3).join(', ')} mm)`, () => {
      const g = new FlameGeometry(B, c, { maxHeight: 0.13, allowSparkOutsideBore: true });
      const f = newFlameGeometryResult();
      const e = newFlameGeometryResult();
      const rnd = lcg(29);
      let worstA = 0;
      let worstV = 0;
      let worstFD = 0;
      for (let i = 0; i < 1500; i++) {
        const h = 0.02 + rnd() * 0.11;
        const r = rnd() * g.maxRadius(h);
        g.evaluate(r, h, f);
        g.evaluateExact(r, h, e);
        const Vch = PI * R * R * h;
        if (e.volume > 1e-3 * Vch) worstV = Math.max(worstV, Math.abs(f.volume / e.volume - 1));
        if (e.frontArea > 1e-2 * 2 * R * h) worstA = Math.max(worstA, Math.abs(f.frontArea / e.frontArea - 1));
        for (const k of ['wettedHead', 'wettedPiston', 'wettedLiner'] as const) {
          if (e[k] > 1e-3 * PI * R * R) expect(Math.abs(f[k] / e[k] - 1)).toBeLessThan(1e-3);
        }
        const dr = 1e-7 * r;
        const fd = (g.volume(r + dr, h) - g.volume(r - dr, h)) / (2 * dr);
        const A = g.frontArea(r, h);
        if (A > 1e-3 * 2 * R * h) worstFD = Math.max(worstFD, Math.abs(fd - A) / A);
        // the ball first meets the column at r = d − R
        if (r <= Math.hypot(c[0], c[2]) - R) expect(f.volume).toBe(0);
      }
      console.log(
        `[chamber] FG outside spark: fast vs exact max rel V ${worstV.toExponential(2)}, A_f ${worstA.toExponential(2)}; FD dV/dr ${worstFD.toExponential(2)}`,
      );
      expect(worstV).toBeLessThan(1e-4);
      expect(worstA).toBeLessThan(1e-3);
      expect(worstFD).toBeLessThan(1e-5);
      // radiusForVolume round trip
      for (const h of [0.0254, 0.06, 0.12]) {
        for (const frac of [1e-6, 1e-3, 0.1, 0.5, 0.9, 0.9999]) {
          const V = frac * PI * R * R * h;
          const r = g.radiusForVolume(V, h);
          expect(Math.abs(g.volume(r, h) / V - 1)).toBeLessThan(1e-9);
        }
      }
    });
  }
});

// ---------------------------------------------------------------------------------------------
// L-head chamber vs the oracle
// ---------------------------------------------------------------------------------------------

/** Scale of quantity k at height h for relative tolerances: its maximum over r (surface area for walls). */
function scaleOf(ch: CombustionChamber, k: string, h: number): number {
  const a = ch.surfaceAreas(h, new Float64Array(N_WALL_SURFACES));
  const s = SURF.find(([n]) => n === k);
  if (s) return a[s[1]];
  if (k === 'V') return ch.chamberVolume(h);
  // front area: max over a scan
  let m = 0;
  const rm = ch.maxRadius(h);
  for (let i = 1; i <= 400; i++) m = Math.max(m, ch.frontArea((rm * i) / 400, h));
  return m;
}

function fastValue(res: ChamberResult, k: string): number {
  if (k === 'V') return res.volume;
  if (k === 'Af') return res.frontArea;
  const s = SURF.find(([n]) => n === k)!;
  return res.wetted[s[1]];
}

describe('L-head chamber vs the oracle (test/fixtures/combustion_geometry_lhead_mc.json)', () => {
  it(`Monte Carlo (${lheadFx.nSamples} samples) agrees with the quadrature oracle (max ${lheadFx.maxDeviationSigma}σ); oracle A_f = dV/dr`, () => {
    expect(lheadFx.maxDeviationSigma).toBeLessThan(5);
    expect(lheadFx.maxOracleFdError).toBeLessThan(2e-5);
  });

  const KEYS = ['V', 'Af', 'head', 'piston', 'liner', 'block', 'intakeValve', 'exhaustValve'] as const;
  // fast table vs dense quadrature: |Δ| ≤ tol·max(|ref|, 1 % of the quantity's scale)
  const TOL: Record<(typeof KEYS)[number], number> = {
    V: 2e-4,
    Af: 2e-3,
    head: 2e-3,
    piston: 2e-3,
    liner: 2e-3,
    block: 2e-3,
    intakeValve: 2e-3,
    exhaustValve: 2e-3,
  };

  for (const G of lheadFx.geometries) {
    it(`${G.name}: pocket plan area (Green) = independent quadrature; fast vs quad and MC; column part vs FlameGeometry`, () => {
      const t0 = performance.now();
      const ch = chamberOf(G);
      const ms = performance.now() - t0;
      const m = lHeadPlanMetrics(G.bore, specOf(G).geometry.lHead!);
      expect(Math.abs(m.pocketArea / G.pocketArea - 1)).toBeLessThan(1e-12);
      const res = newChamberResult();
      const fr = newFlameGeometryResult();
      const fe = newFlameGeometryResult();
      const worst: Record<string, number> = {};
      let worstMc = 0;
      const scales = new Map<string, number>();
      for (const c of G.cases as FxCase[]) {
        ch.evaluate(c.r, c.h, res);
        for (const k of KEYS) {
          const [mc, se, ref] = (c as unknown as Record<string, (number | null)[]>)[k] as [number | null, number | null, number];
          const sk = `${k}@${c.h}`;
          let sc = scales.get(sk);
          if (sc === undefined) {
            sc = scaleOf(ch, k, c.h);
            scales.set(sk, sc);
          }
          const v = fastValue(res, k);
          const e = Math.abs(v - ref) / Math.max(Math.abs(ref), 0.01 * sc, 1e-12);
          worst[k] = Math.max(worst[k] ?? 0, e);
          expect(e, `${G.name} ${k} r=${c.r} h=${c.h}: ${v} vs ${ref}`).toBeLessThan(TOL[k]);
          if (mc !== null && se !== null && mc > 0 && mc >= 0.01 * sc) {
            const tol = 5e-3 + (5 * se) / mc;
            worstMc = Math.max(worstMc, Math.abs(v - mc) / mc / tol);
            expect(Math.abs(v - mc) / mc, `${G.name} ${k} vs MC r=${c.r} h=${c.h}`).toBeLessThan(tol);
          }
        }
        // the bore column alone: FlameGeometry fast (table) and exact paths
        const fg = ch.flame;
        fg.evaluate(c.r, c.h, fr);
        fg.evaluateExact(c.r, c.h, fe);
        const R = G.bore / 2;
        const col: [number, number, number, number][] = [
          [fr.volume, fe.volume, c.Vcol, PI * R * R * c.h],
          [fr.frontArea, fe.frontArea, c.Afcol, 2 * PI * R * c.h],
          [fr.wettedHead, fe.wettedHead, c.Wh, PI * R * R],
          [fr.wettedPiston, fe.wettedPiston, c.Wp, PI * R * R],
          [fr.wettedLiner, fe.wettedLiner, c.Wl, 2 * PI * R * c.h],
        ];
        for (const [fast, ex, ref, sc] of col) {
          expect(Math.abs(ex - ref) / Math.max(Math.abs(ref), 1e-6 * sc)).toBeLessThan(1e-7);
          expect(Math.abs(fast - ref) / Math.max(Math.abs(ref), 0.01 * sc)).toBeLessThan(1e-3);
        }
      }
      console.log(
        `[chamber] ${G.name}: build ${ms.toFixed(0)} ms (tables ${JSON.stringify(ch.tableSizes)}), ${G.cases.length} cases — max rel err vs quad ` +
          KEYS.map((k) => `${k} ${(worst[k] ?? 0).toExponential(1)}`).join(', ') +
          `; vs MC max |Δ|/(0.5 % + 5σ) ${worstMc.toFixed(2)}`,
      );
    });
  }
});

// ---------------------------------------------------------------------------------------------
// L-head invariants
// ---------------------------------------------------------------------------------------------

describe('L-head chamber invariants (Model T spec)', () => {
  const ch = createChamber(MODEL_T) as LHeadChamber;
  const lh = MODEL_T.geometry.lHead!;
  const [, cy] = MODEL_T.sparkPlug.gapCenter;
  const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
  const HS = [hTdc, 0.03, -lh.deckY, 0.05, 0.09, 0.127];

  it('table build is fast enough to share per engine (reported) and clones share it', () => {
    const t0 = performance.now();
    const fresh = createChamber(MODEL_T, { cache: false }) as LHeadChamber;
    const ms = performance.now() - t0;
    console.log(`[chamber] Model T L-head table build ${ms.toFixed(0)} ms, sizes ${JSON.stringify(fresh.tableSizes)}`);
    expect(ms).toBeLessThan(3000);
    const t1 = performance.now();
    const c2 = ch.clone();
    expect(performance.now() - t1).toBeLessThan(20);
    expect(c2.flame).not.toBe(ch.flame);
    const a = newChamberResult();
    const b = newChamberResult();
    const rnd = lcg(9);
    for (let i = 0; i < 400; i++) {
      const h = hTdc + rnd() * 0.1;
      const r = rnd() * ch.maxRadius(h);
      fresh.evaluate(r, h, a);
      c2.evaluate(r, h, b);
      expect(b).toEqual(a);
    }
  });

  it('dV/dr = A_f on the fast path (central differences, 3000 random points)', () => {
    const rnd = lcg(7);
    let worst = 0;
    for (let n = 0; n < 3000; n++) {
      const h = hTdc + rnd() * (0.127 - hTdc);
      const r = (0.001 + 0.998 * rnd()) * ch.maxRadius(h);
      const dr = 1e-7 * r;
      const fd = (ch.volume(r + dr, h) - ch.volume(r - dr, h)) / (2 * dr);
      const A = ch.frontArea(r, h);
      const scale = Math.max(A, 1e-3 * 4 * PI * r * r);
      const e = Math.abs(fd - A) / scale;
      // points within dr of a FlameGeometry fallback-band edge see its tiny table/exact jump
      if (e < 1e-3) worst = Math.max(worst, e);
      else expect(e * scale * 2 * dr).toBeLessThan(1e-9 * ch.chamberVolume(h));
    }
    console.log(`[chamber] L-head |FD dV/dr − A_f| / A_f max ${worst.toExponential(2)}`);
    expect(worst).toBeLessThan(1e-5);
  });

  it('V monotone and continuous in r; V(r ≥ r_max) = chamberVolume(h) exactly; wetted → surface areas', () => {
    const res = newChamberResult();
    const areas = newWallSurfaceArray();
    for (const h of HS) {
      const rm = ch.maxRadius(h);
      let prev = 0;
      let prevA = 0;
      let jump = 0;
      let Amax = 0;
      const n = 6000;
      for (let i = 1; i <= n; i++) Amax = Math.max(Amax, ch.frontArea((rm * i) / n, h));
      for (let i = 1; i <= n; i++) {
        const r = (rm * i) / n;
        const v = ch.volume(r, h);
        // (FlameGeometry's table may leave A_f ≈ −2e-5·max just short of the far corner, cf. its 1e-3 spec)
        expect(v).toBeGreaterThanOrEqual(prev - 1e-12 * ch.chamberVolume(h));
        const A = ch.frontArea(r, h);
        expect(A).toBeGreaterThanOrEqual(-1e-4 * Amax);
        if (i > 1) jump = Math.max(jump, Math.abs(A - prevA));
        prev = v;
        prevA = A;
      }
      // A_f has no jumps beyond what the step dr·dA/dr allows (≈ 4π r dr)
      expect(jump).toBeLessThan(4 * PI * rm * (rm / n) * 3);
      for (const r of [rm, rm * 1.0001, 1]) {
        ch.evaluate(r, h, res);
        expect(res.volume).toBe(ch.chamberVolume(h));
        expect(Math.abs(res.frontArea)).toBeLessThan(1e-20);
        ch.surfaceAreas(h, areas);
        for (let s = 0; s < N_WALL_SURFACES; s++) {
          expect(Math.abs(res.wetted[s] - areas[s])).toBeLessThan(1e-9 * (areas[s] + 1e-6));
          expect(res.burnedFraction[s]).toBeCloseTo(areas[s] > 0 ? 1 : 0, 9);
        }
        expect(res.creviceBurnedFraction).toBeCloseTo(1, 6);
      }
      expect(ch.radiusForVolume(ch.chamberVolume(h), h)).toBe(rm);
    }
  });

  it('continuous in h (piston-crown slab ≈ crown lens × dh) and at the deck (crown passing the block face)', () => {
    const res1 = newChamberResult();
    const res2 = newChamberResult();
    for (const h of [hTdc, 0.03, -lh.deckY, 0.06]) {
      for (const r of [0.02, 0.05, 0.08]) {
        const dh = 1e-7;
        const v1 = ch.volume(r, h - dh);
        const v2 = ch.volume(r, h + dh);
        ch.evaluate(r, h, res1);
        const lens = ch.flame.evaluate(r, h, newFlameGeometryResult()).wettedPiston;
        expect(Math.abs((v2 - v1) / (2 * dh) - lens)).toBeLessThan(1e-4 * PI * (MODEL_T.geometry.bore / 2) ** 2);
        ch.evaluate(r, h - 1e-9, res1);
        ch.evaluate(r, h + 1e-9, res2);
        for (let s = 0; s < N_WALL_SURFACES; s++) expect(Math.abs(res2.wetted[s] - res1.wetted[s])).toBeLessThan(1e-6 * 0.03);
      }
    }
  });

  it('analytic limits: whole sphere below the pocket roof distance, sphere minus a cap until the next wall', () => {
    const tRoof = lh.pocket.roofY - cy; // 4 mm
    const tFloor = cy - lh.deckY;
    const dBore = Math.hypot(MODEL_T.sparkPlug.gapCenter[0], MODEL_T.sparkPlug.gapCenter[2]) - MODEL_T.geometry.bore / 2;
    const rNext = Math.min(tFloor, dBore);
    for (const h of [hTdc, 0.08]) {
      expect(ch.inscribedRadius(h)).toBeCloseTo(tRoof, 12);
      for (const r of [1e-4, 0.001, 0.002, 0.0039]) {
        const V = (4 / 3) * PI * r ** 3;
        expect(Math.abs(ch.volume(r, h) / V - 1)).toBeLessThan(1e-10);
        expect(Math.abs(ch.frontArea(r, h) / (4 * PI * r * r) - 1)).toBeLessThan(1e-9);
        expect(ch.radiusForVolume(V, h)).toBeCloseTo(r, 15);
      }
      for (const r of [0.0041, 0.005, 0.007, 0.9 * rNext]) {
        const a = r - tRoof; // cap cut off by the pocket roof
        const V = (4 / 3) * PI * r ** 3 - (PI * a * a * (3 * r - a)) / 3;
        const A = 4 * PI * r * r - 2 * PI * r * a;
        expect(Math.abs(ch.volume(r, h) / V - 1), `V r=${r}`).toBeLessThan(1e-7);
        expect(Math.abs(ch.frontArea(r, h) / A - 1), `A r=${r}`).toBeLessThan(1e-5);
        // wetted pocket roof = the cap's base disc πρ² (head surface)
        const res = ch.evaluate(r, h, newChamberResult());
        expect(Math.abs(res.wetted[WALL_HEAD] / (PI * (r * r - tRoof * tRoof)) - 1)).toBeLessThan(1e-7);
      }
    }
  });

  it('radiusForVolume round trip (1e-9) across regimes and heights', () => {
    for (const h of HS) {
      const Vch = ch.chamberVolume(h);
      for (const frac of [1e-9, 1e-6, 1e-3, 0.02, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99, 0.99999]) {
        const V = frac * Vch;
        const r = ch.radiusForVolume(V, h, 0.5 * ch.maxRadius(h));
        expect(Math.abs(ch.volume(r, h) / V - 1)).toBeLessThan(1e-9);
      }
      expect(ch.radiusForVolume(0, h)).toBe(0);
    }
  });

  it('chamber volume and depth: V − V_crevice consistent with SliderCrank; meanDepth(h) = V_ch/plan area', () => {
    const g = MODEL_T.geometry;
    const kin = new SliderCrank(sliderCrankGeometry(g), g.compressionRatio);
    for (const th of [0, 0.5, 1.5, Math.PI]) {
      const h = kin.clearanceHeight(th);
      expect(Math.abs(ch.chamberVolume(h) / (kin.volume(th) - g.creviceVolume) - 1)).toBeLessThan(1e-12);
      expect(ch.meanDepth(h)).toBeCloseTo(ch.chamberVolume(h) / ch.planformArea, 15);
    }
    // at TDC the mean depth is ≈ 20 mm (pancake-equivalent h_TDC would be ≈ 34 mm)
    expect(ch.meanDepth(hTdc)).toBeGreaterThan(0.019);
    expect(ch.meanDepth(hTdc)).toBeLessThan(0.022);
  });

  it('surface areas at TDC: ≈ 340 cm² (vs 245 cm² for the pancake equivalent); crown above the deck exposes the piston side', () => {
    const a = ch.surfaceAreas(hTdc, newWallSurfaceArray());
    let total = 0;
    for (let s = 0; s < N_WALL_SURFACES; s++) total += a[s];
    console.log(`[chamber] Model T surface areas at TDC, cm²: ${Array.from(a, (v) => (v * 1e4).toFixed(1)).join(', ')} (total ${(total * 1e4).toFixed(1)})`);
    expect(total * 1e4).toBeGreaterThan(300);
    expect(total * 1e4).toBeLessThan(380);
    expect(a[WALL_LINER]).toBe(0);
    const R = MODEL_T.geometry.bore / 2;
    expect(a[WALL_PISTON]).toBeGreaterThan(PI * R * R); // + side facing the pocket
    const b = ch.surfaceAreas(0.06, newWallSurfaceArray());
    expect(b[WALL_PISTON]).toBeCloseTo(PI * R * R, 15);
    expect(b[WALL_LINER]).toBeCloseTo(2 * PI * R * (0.06 + lh.deckY), 15);
  });

  it('spark in the bore column: the pocket is not reached before the bore circle — identical to the disc column', () => {
    const spec: Spec = { ...MODEL_T, sparkPlug: { ...MODEL_T.sparkPlug, gapCenter: [-0.02, -0.005, 0.01] } };
    const cb = createChamber(spec);
    const d = Math.hypot(-0.02, 0.01);
    const R = MODEL_T.geometry.bore / 2;
    for (const h of [hTdc, 0.06]) {
      for (const r of [0.002, 0.004, 0.9 * (R - d)]) {
        expect(cb.volume(r, h)).toBe(cb.flame.volume(r, h));
        expect(cb.frontArea(r, h)).toBe(cb.flame.frontArea(r, h));
      }
    }
  });

  it('rejects a spark outside the chamber and valves outside the pocket', () => {
    expect(() => createChamber({ ...MODEL_T, sparkPlug: { ...MODEL_T.sparkPlug, gapCenter: [-0.06, -0.01, 0] } }, { cache: false })).toThrow();
    expect(() =>
      createChamber({ ...MODEL_T, intakeValve: { ...MODEL_T.intakeValve, position: [-0.05, -0.02] } }, { cache: false }),
    ).toThrow();
  });

  it('allocation-free hot paths (heap growth over 2e5 calls stays small)', () => {
    const res = newChamberResult();
    const run = (n: number): number => {
      let s = 0;
      for (let i = 0; i < n; i++) {
        const h = hTdc + (i % 97) * 1e-3;
        const r = 0.001 + (i % 1013) * 1e-4;
        ch.evaluate(r, h, res);
        s += res.volume + ch.radiusForVolume(0.3 * ch.chamberVolume(h), h, r);
      }
      return s;
    };
    run(20000); // warm up (JIT)
    const g0 = (globalThis as { gc?: () => void }).gc;
    g0?.();
    const m0 = process.memoryUsage().heapUsed;
    const t0 = performance.now();
    const s = run(200000);
    const us = ((performance.now() - t0) / 200000) * 1e3;
    const grow = process.memoryUsage().heapUsed - m0;
    console.log(`[chamber] L-head evaluate + radiusForVolume ${us.toFixed(2)} µs per pair, heap growth ${(grow / 1e6).toFixed(2)} MB (${s > 0 ? 'ok' : ''})`);
    // V8 may box doubles across calls (DESIGN.md); real per-call allocation would be ≥ 200000 × 32 B = 6.4 MB
    expect(grow).toBeLessThan(64e6);
  });
});

// ---------------------------------------------------------------------------------------------
// Model T closed-vessel burn at TDC (two-zone ideal gas + the chamber geometry + Keck entrainment)
// ---------------------------------------------------------------------------------------------

describe('closed-vessel burn in the Model T chamber at TDC (entrainment model on the chamber interface)', () => {
  // as entrainment.test.ts (CFR) — constant volume, isentropic unburned gas, q = 2.75 MJ/kg of mixture
  const lh = MODEL_T.geometry.lHead!;
  const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
  const gamma = 1.3;
  const p0 = 8e5;
  const Tu0 = 600;
  const Rg = 287;
  const q = 2.75e6;
  const rhoU0 = p0 / (Rg * Tu0);

  function run(ch: CombustionChamber, h: number, u0: number, SL: number, lambda: number) {
    const V = ch.chamberVolume(h);
    const m = rhoU0 * V;
    const inp = newEntrainmentInputs();
    inp.mTotal = m;
    const rates = newEntrainmentRates();
    let rfPrev = 0;
    const f = (me: number, mb: number, o: number[]) => {
      const p = p0 + ((gamma - 1) * mb * q) / V;
      const rhoU = rhoU0 * Math.pow(p / p0, 1 / gamma);
      const Vf = Math.min(V - (m - me) / rhoU, V);
      const rf = ch.radiusForVolume(Math.max(Vf, 0), h, rfPrev > 0 ? rfPrev : undefined);
      rfPrev = rf;
      inp.rhoU = rhoU;
      inp.uPrime = rapidDistortionIntensity(u0, rhoU0, rhoU);
      inp.SL = SL;
      inp.lambda = lambda;
      inp.frontArea = ch.frontArea(rf, h);
      inp.me = me;
      inp.mb = mb;
      entrainmentRates(inp, rates);
      o[0] = rates.dme;
      o[1] = rates.dmb;
    };
    const rhoB0 = (rhoU0 * Tu0) / (Tu0 + q / (Rg / (gamma - 1)));
    let me = rhoB0 * (4 / 3) * PI * 1e-9; // 1 mm burned kernel
    let mb = me;
    const k1 = [0, 0];
    const k2 = [0, 0];
    const k3 = [0, 0];
    const k4 = [0, 0];
    const dt = 4e-6;
    let t10 = NaN;
    let t90 = NaN;
    let ok = true;
    let rfAtFull = 0;
    for (let i = 0; i < 15000 && mb < 0.9999 * m; i++) {
      f(me, mb, k1);
      f(me + (dt / 2) * k1[0], mb + (dt / 2) * k1[1], k2);
      f(me + (dt / 2) * k2[0], mb + (dt / 2) * k2[1], k3);
      f(me + dt * k3[0], mb + dt * k3[1], k4);
      me += (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
      mb += (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
      if (me > m) me = m;
      if (!(mb <= me * (1 + 1e-12) && mb >= 0)) ok = false;
      const t = (i + 1) * dt;
      if (Number.isNaN(t10) && mb >= 0.1 * m) t10 = t;
      if (Number.isNaN(t90) && mb >= 0.9 * m) t90 = t;
      if (me >= m && rfAtFull === 0) rfAtFull = ch.radiusForVolume(V, h);
    }
    return { m, me, mb, t10, t90, ok, rfAtFull };
  }

  it('burns out; the front reaches the far corner exactly when all charge is entrained; slower than a central plug in an equal pancake', () => {
    const ch = createChamber(MODEL_T);
    const lt = run(ch, hTdc, 1.0, 0.5, 1e-3);
    expect(lt.ok).toBe(true);
    expect(lt.mb / lt.m).toBeGreaterThan(0.999);
    expect(lt.rfAtFull).toBe(ch.maxRadius(hTdc));
    // the same clearance volume as a bore-diameter pancake (the 'flat-disc' surrogate) with a central plug
    const R = MODEL_T.geometry.bore / 2;
    const hEq = ch.chamberVolume(hTdc) / (PI * R * R);
    const disc = createChamber(
      {
        ...MODEL_T,
        geometry: { ...MODEL_T.geometry, chamber: 'flat-disc', lHead: undefined },
        sparkPlug: { ...MODEL_T.sparkPlug, gapCenter: [0, -hEq / 2, 0] },
      },
      { maxHeight: 0.05 },
    );
    const pc = run(disc, hEq, 1.0, 0.5, 1e-3);
    expect(Math.abs(disc.chamberVolume(hEq) / ch.chamberVolume(hTdc) - 1)).toBeLessThan(1e-12);
    expect(pc.mb / pc.m).toBeGreaterThan(0.999);
    expect(lt.t90 - lt.t10).toBeGreaterThan(1.3 * (pc.t90 - pc.t10)); // ≈ 110 mm flame travel vs ≈ 50 mm
    console.log(
      `[chamber] closed-vessel TDC burn (u′₀ 1 m/s, S_L 0.5 m/s, λ 1 mm): Model T pocket plug 0–10 % ${(lt.t10 * 1e3).toFixed(2)} ms, ` +
        `10–90 % ${((lt.t90 - lt.t10) * 1e3).toFixed(2)} ms; equal-volume pancake, central plug 0–10 % ${(pc.t10 * 1e3).toFixed(2)} ms, ` +
        `10–90 % ${((pc.t90 - pc.t10) * 1e3).toFixed(2)} ms`,
    );
  });
});
