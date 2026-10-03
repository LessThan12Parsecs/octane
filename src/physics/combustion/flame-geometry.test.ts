import { describe, expect, it } from 'vitest';
import mcFx from '../../../test/fixtures/combustion_geometry_mc.json';
import quadFx from '../../../test/fixtures/combustion_geometry_quad.json';
import {
  FlameGeometry,
  arcInsideBore,
  arcOfBoreInside,
  defaultInnerCells,
  lensArea,
  newFlameGeometryResult,
  type FlameGeometryResult,
} from './flame-geometry';

const KEYS = ['volume', 'frontArea', 'wettedHead', 'wettedPiston', 'wettedLiner'] as const;
const FX_KEYS = ['V', 'Af', 'Wh', 'Wp', 'Wl'] as const;
type Key = (typeof KEYS)[number];

const PI = Math.PI;
const IN = 0.0254;
const CFR_BORE = 3.25 * IN;
const CFR_CENTER: [number, number, number] = [0.9 * (CFR_BORE / 2), -0.002, 0];

const geoms = new Map<string, FlameGeometry>();
function geometry(name: string, bore: number, center: number[]): FlameGeometry {
  let g = geoms.get(name);
  if (!g) {
    g = new FlameGeometry(bore, center as [number, number, number], { maxHeight: 0.16 });
    geoms.set(name, g);
  }
  return g;
}

/** Max of each quantity over r ∈ (0, r_max(h)] at fixed h (fast path scan). */
function maxima(g: FlameGeometry, h: number): Record<Key, number> {
  const out = newFlameGeometryResult();
  const mx: Record<Key, number> = { volume: 0, frontArea: 0, wettedHead: 0, wettedPiston: 0, wettedLiner: 0 };
  const rm = g.maxRadius(h);
  for (let i = 1; i <= 800; i++) {
    g.evaluate((rm * i) / 800, h, out);
    for (const k of KEYS) mx[k] = Math.max(mx[k], out[k]);
  }
  return mx;
}

describe('lensArea', () => {
  it('matches the classical acos formula and the limits', () => {
    const R = 0.04;
    const d = 0.03;
    for (const rho of [0.005, 0.0101, 0.02, 0.05, 0.069, 0.0699]) {
      const ca = (d * d + rho * rho - R * R) / (2 * d * rho);
      const cb = (d * d + R * R - rho * rho) / (2 * d * R);
      const k = (-d + rho + R) * (d + rho - R) * (d - rho + R) * (d + rho + R);
      const ref =
        rho <= R - d
          ? PI * rho * rho
          : rho >= R + d
            ? PI * R * R
            : rho * rho * Math.acos(ca) + R * R * Math.acos(cb) - 0.5 * Math.sqrt(k);
      expect(lensArea(rho, R, d)).toBeCloseTo(ref, 14);
    }
    expect(lensArea(0.009, R, d)).toBeCloseTo(PI * 0.009 * 0.009, 15);
    expect(lensArea(0.08, R, d)).toBeCloseTo(PI * R * R, 15);
  });
});

describe('review regressions', () => {
  // CFR_F1 spec gap centre (engines/cfr.ts): [0, −2 mm, R − 1 mm] — 1 mm from the liner.
  const R = CFR_BORE / 2;
  const SPEC: [number, number, number] = [0, -0.002, R - 0.001];

  it('lens helpers are free of the ρ² − R² cancellation for a near-axis spark (d = 1e-6 R)', () => {
    // dL/dρ = ρΦ(ρ) and the arc identities, sampled across the 2d-wide lens range by finite differences.
    // The former (ρ² − c1c2)/(2d) form lost log10(R/d) digits: the FD check failed at the 1e-2 level.
    const d = 1e-6 * R;
    const c1 = R - d;
    let worst = 0;
    let prev = 0;
    for (let i = 1; i < 200; i++) {
      const rho = c1 + (2 * d * i) / 200;
      const L = lensArea(rho, R, d);
      expect(L).toBeGreaterThan(prev);
      prev = L;
      const dr = 1e-4 * d;
      const fd = (lensArea(rho + dr, R, d) - lensArea(rho - dr, R, d)) / (2 * dr);
      worst = Math.max(worst, Math.abs(fd / (rho * arcInsideBore(rho, R, d)) - 1));
      // chord seen from both centres: ρ sin(Φ/2) = R sin(Ψ/2)
      const s1 = rho * Math.sin(0.5 * arcInsideBore(rho, R, d));
      const s2 = R * Math.sin(0.5 * arcOfBoreInside(rho, R, d));
      expect(Math.abs(s1 - s2)).toBeLessThan(1e-12 * R);
    }
    expect(worst).toBeLessThan(1e-4);
  });

  it('near-axis spark (d = 1e-7 R) builds promptly and agrees with the central closed form', () => {
    // Before the fix the adaptive build could not reach its tolerance through the cancellation noise and
    // ran into its evaluation caps: > 10 min at d = 1e-6 R, 14 s at d = 3e-5 R.
    const t0 = performance.now();
    const gn = new FlameGeometry(0.1, [0.05e-7, -0.005, 0]);
    expect(performance.now() - t0).toBeLessThan(20000);
    const gc = new FlameGeometry(0.1, [0, -0.005, 0]);
    const a = newFlameGeometryResult();
    const b = newFlameGeometryResult();
    for (const r of [0.01, 0.03, 0.049, 0.0499, 0.0501, 0.051, 0.06, 0.07]) {
      for (const h of [0.004, 0.02, 0.08]) {
        gn.evaluate(r, h, a);
        gc.evaluate(r, h, b);
        for (const k of KEYS) {
          if (b[k] > 0) expect(Math.abs(a[k] - b[k]) / b[k], `${k} r=${r} h=${h}`).toBeLessThan(1e-6);
          expect(Number.isFinite(a[k])).toBe(true);
        }
      }
    }
  });

  it('inner-table resolution adapts to a spark near the liner (defaultInnerCells)', () => {
    expect(defaultInnerCells(R, 0)).toBe(64);
    expect(defaultInnerCells(R, 0.5 * R)).toBe(64);
    expect(defaultInnerCells(R, 0.9 * R)).toBe(68);
    expect(defaultInnerCells(R, R - 0.001)).toBe(144);
    expect(defaultInnerCells(R, 0.99 * R)).toBe(226);
    expect(defaultInnerCells(R, 0.9999 * R)).toBe(512);
  });

  it('CFR_F1 spec spark (1 mm from the liner): fast front area ≤ 5e-4 of exact in the early flame', () => {
    // The flame touches the liner at r = 1 mm, i.e. the whole turbulent flame lives in the lens regime.
    // With the former fixed 64 inner cells the error just after contact reached 2.5e-3 (8.6e-3 at 0.99R).
    const g = new FlameGeometry(CFR_BORE, SPEC, { maxHeight: 0.16 });
    const c1 = 0.001;
    const f = newFlameGeometryResult();
    const e = newFlameGeometryResult();
    let worst = 0;
    for (let i = 0; i < 600; i++) {
      const h = 0.0025 + 0.0375 * ((i * 0.7548776662) % 1);
      const r = c1 + 0.01 * ((i * 0.569840291) % 1);
      g.evaluate(r, h, f);
      g.evaluateExact(r, h, e);
      worst = Math.max(worst, Math.abs(f.frontArea / e.frontArea - 1));
      expect(Math.abs(f.volume / e.volume - 1)).toBeLessThan(5e-5);
    }
    console.log(`[flame-geometry] CFR spec early flame: max |A_f fast/exact − 1| = ${worst.toExponential(2)}`);
    expect(worst).toBeLessThan(5e-4);
  });

  it('degenerate clearance height h ≤ 0 gives an empty chamber, not negative volumes', () => {
    const g = new FlameGeometry(CFR_BORE, SPEC, { maxHeight: 0.16 });
    const out = newFlameGeometryResult();
    for (const h of [0, -0.001]) {
      for (const ev of ['evaluate', 'evaluateExact'] as const) {
        g[ev](0.01, h, out);
        for (const k of KEYS) expect(out[k]).toBe(0);
      }
      expect(g.radiusForVolume(1e-6, h)).toBe(0);
    }
  });
});

describe('FlameGeometry vs Monte-Carlo oracle (test/fixtures/combustion_geometry_mc.json)', () => {
  it(`MC (${mcFx.nSamples} samples/quantity) agrees with the quadrature oracle (max ${mcFx.maxDeviationSigma}σ)`, () => {
    expect(mcFx.maxDeviationSigma).toBeLessThan(5);
  });

  for (const G of mcFx.geometries) {
    it(`${G.name}: exact path = quad (1e-9) and within 5σ of MC; fast path ≤ 0.5 % of MC where ≥ 1 % of max`, () => {
      const g = geometry(G.name, G.bore, G.center);
      const ex = newFlameGeometryResult();
      const fast = newFlameGeometryResult();
      let worstFast = 0;
      let worstFastQuad = 0;
      let worstExact = 0;
      let worstSigma = 0;
      for (const c of G.cases) {
        g.evaluateExact(c.r, c.h, ex);
        g.evaluate(c.r, c.h, fast);
        const mx = maxima(g, c.h);
        KEYS.forEach((k, i) => {
          const [mc, se, quad] = (c as unknown as Record<string, number[]>)[FX_KEYS[i]];
          const scale = Math.max(mx[k], 1e-30);
          // exact vs quadrature oracle
          const eEx = Math.abs(ex[k] - quad) / Math.max(Math.abs(quad), 1e-6 * scale);
          worstExact = Math.max(worstExact, eEx);
          expect(eEx).toBeLessThan(1e-8);
          // exact vs MC
          if (se > 0) {
            const ns = Math.abs(ex[k] - mc) / se;
            worstSigma = Math.max(worstSigma, ns);
            expect(ns).toBeLessThan(5);
          }
          // fast vs MC where the quantity is ≥ 1 % of its max over r: |Δ| ≤ 0.5 % + 5σ
          if (mc >= 0.01 * scale) {
            const e = Math.abs(fast[k] - mc) / mc;
            const tol = 5e-3 + (5 * se) / mc;
            worstFast = Math.max(worstFast, e / tol);
            expect(e, `${k} r=${c.r} h=${c.h}`).toBeLessThan(tol);
            worstFastQuad = Math.max(worstFastQuad, Math.abs(fast[k] - quad) / quad);
          }
        });
      }
      console.log(
        `[flame-geometry] ${G.name}: ${G.cases.length} MC cases — fast: max |Δ|/(0.5 % + 5σ) vs MC ${worstFast.toFixed(3)}, ` +
          `max rel vs quad ${worstFastQuad.toExponential(2)}; exact vs quad max rel ${worstExact.toExponential(2)}, ` +
          `exact vs MC max ${worstSigma.toFixed(2)}σ`,
      );
    });
  }
});

describe('FlameGeometry fast table vs dense quadrature oracle (combustion_geometry_quad.json)', () => {
  for (const G of quadFx.geometries) {
    it(`${G.name}: 400 random (r, h), h ∈ [0.3, 160] mm — fast ≤ 0.5 % where ≥ 1 % of max; exact ≤ 1e-8`, () => {
      const g = geometry(G.name, G.bore, G.center);
      const ex = newFlameGeometryResult();
      const fast = newFlameGeometryResult();
      const worst: Record<Key, number> = { volume: 0, frontArea: 0, wettedHead: 0, wettedPiston: 0, wettedLiner: 0 };
      for (const row of G.rows) {
        const [r, h] = row;
        g.evaluate(r, h, fast);
        g.evaluateExact(r, h, ex);
        const mx = maxima(g, h);
        KEYS.forEach((k, i) => {
          const ref = row[2 + i];
          const scale = Math.max(mx[k], 1e-30);
          expect(Math.abs(ex[k] - ref) / Math.max(Math.abs(ref), 1e-6 * scale)).toBeLessThan(1e-8);
          if (ref >= 0.01 * scale) {
            const e = Math.abs(fast[k] - ref) / ref;
            worst[k] = Math.max(worst[k], e);
            expect(e, `${k} r=${r} h=${h}`).toBeLessThan(5e-3);
          }
        });
      }
      console.log(
        `[flame-geometry] ${G.name} dense: max rel err (fast) ` +
          KEYS.map((k) => `${k} ${worst[k].toExponential(2)}`).join(', '),
      );
    });
  }
});

describe('FlameGeometry analytic limits', () => {
  const g = new FlameGeometry(CFR_BORE, CFR_CENTER, { maxHeight: 0.16 });
  const R = CFR_BORE / 2;
  const d = CFR_CENTER[0];
  const b = -CFR_CENTER[1];
  const out = newFlameGeometryResult();

  it('ball strictly inside: V = 4πr³/3, A_f = 4πr², no wetted walls', () => {
    const h = 0.03;
    for (const r of [1e-5, 1e-4, 1e-3, 0.0019]) {
      for (const ev of ['evaluate', 'evaluateExact'] as const) {
        g[ev](r, h, out);
        expect(out.volume).toBeCloseTo((4 / 3) * PI * r ** 3, 15);
        expect(out.frontArea / (4 * PI * r * r)).toBeCloseTo(1, 13);
        expect(out.wettedHead + out.wettedPiston + out.wettedLiner).toBe(0);
      }
    }
  });

  it('ball cut only by the head plane: spherical cap formulas', () => {
    const r = 0.0035; // > b = 2 mm, < R − d = 4.13 mm
    const cap = (PI * (r - b) ** 2 * (3 * r - (r - b))) / 3;
    g.evaluate(r, 0.03, out);
    expect(out.volume / ((4 / 3) * PI * r ** 3 - cap)).toBeCloseTo(1, 12);
    expect(out.frontArea / (2 * PI * r * (r + b))).toBeCloseTo(1, 12);
    expect(out.wettedHead / (PI * (r * r - b * b))).toBeCloseTo(1, 12);
  });

  it('centre on the head: hemisphere', () => {
    const gh = new FlameGeometry(CFR_BORE, [0.2 * R, 0, 0.1 * R]);
    for (const r of [0.001, 0.01, 0.02]) {
      for (const ev of ['evaluate', 'evaluateExact'] as const) {
        gh[ev](r, 0.05, out);
        expect(out.volume / ((2 / 3) * PI * r ** 3)).toBeCloseTo(1, 12);
        expect(out.frontArea / (2 * PI * r * r)).toBeCloseTo(1, 12);
        expect(out.wettedHead / (PI * r * r)).toBeCloseTo(1, 12);
      }
    }
  });

  it('r ≥ farthest corner: V = πR²h, A_f = 0, all walls wetted', () => {
    for (const h of [0.0003, 0.019, 0.13]) {
      const rc = Math.hypot(R + d, Math.max(b, Math.abs(h - b)));
      expect(g.maxRadius(h)).toBeCloseTo(rc, 15);
      for (const r of [rc, 1.2 * rc]) {
        for (const ev of ['evaluate', 'evaluateExact'] as const) {
          g[ev](r, h, out);
          expect(out.volume / (PI * R * R * h)).toBeCloseTo(1, 10);
          expect(out.frontArea).toBeLessThan(1e-12);
          expect(out.wettedHead / (PI * R * R)).toBeCloseTo(1, 10);
          expect(out.wettedPiston / (PI * R * R)).toBeCloseTo(1, 10);
          expect(out.wettedLiner / (2 * PI * R * h)).toBeCloseTo(1, 10);
        }
      }
      // continuity approaching the corner
      g.evaluate(rc * (1 - 1e-7), h, out);
      expect(out.volume / (PI * R * R * h)).toBeCloseTo(1, 5);
    }
  });

  it('central spark (d = 0): closed form of ball ∩ cylinder slab', () => {
    const Rc = 0.05;
    const gc = new FlameGeometry(0.1, [0, -0.01, 0]);
    const h = 0.03;
    for (const r of [0.004, 0.012, 0.03, 0.05, 0.052, 0.0535]) {
      // V = ∫ π min(ρ, R)² dt over t ∈ [max(−0.02, −r), min(0.01, r)]
      const lo = Math.max(-0.02, -r);
      const hi = Math.min(0.01, r);
      const t1 = r > Rc ? Math.sqrt(r * r - Rc * Rc) : 0;
      const F = (t: number) => (Math.abs(t) < t1 ? PI * Rc * Rc * t : PI * (r * r * t - (t * t * t) / 3));
      // piecewise antiderivative (continuous at ±t1)
      const Fc = (t: number) => {
        if (t1 === 0) return F(t);
        const s = Math.sign(t);
        const at = Math.abs(t);
        const inner = PI * Rc * Rc * Math.min(at, t1);
        const outer = at > t1 ? PI * (r * r * (at - t1) - (at ** 3 - t1 ** 3) / 3) : 0;
        return s * (inner + outer);
      };
      const V = Fc(hi) - Fc(lo);
      // A_f = 2πr × length of {t ∈ [lo, hi] : |t| ≥ t1} (Φ = 2π inside the bore, 0 where ρ > R)
      const Af = 2 * PI * r * (hi - lo - Math.max(0, Math.min(hi, t1) - Math.max(lo, -t1)));
      gc.evaluate(r, h, out);
      expect(out.volume / V).toBeCloseTo(1, 12);
      expect(out.frontArea / Af).toBeCloseTo(1, 12);
    }
  });
});

describe('FlameGeometry consistency', () => {
  const g = new FlameGeometry(CFR_BORE, CFR_CENTER, { maxHeight: 0.16 });

  it('dV/dr = frontArea (fast path: central differences, 2000 random points)', () => {
    let seed = 7;
    const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
    let worst = 0;
    for (let n = 0; n < 2000; n++) {
      const h = 0.0003 * Math.exp(rnd() * Math.log(0.16 / 0.0003));
      const r = (0.001 + 0.998 * rnd()) * g.maxRadius(h);
      const dr = 1e-7 * r;
      const fd = (g.volume(r + dr, h) - g.volume(r - dr, h)) / (2 * dr);
      const A = g.frontArea(r, h);
      const scale = Math.max(A, 1e-3 * 4 * PI * r * r);
      // points within dr of a fallback-band edge see the tiny table/exact jump; exclude them
      const e = Math.abs(fd - A) / scale;
      if (e < 1e-3) worst = Math.max(worst, e);
      else expect(e * scale * 2 * dr).toBeLessThan(1e-9 * PI * (CFR_BORE / 2) ** 2 * h); // jump ≤ 1e-9 V_ch
    }
    console.log(`[flame-geometry] |FD dV/dr − A_f| / A_f max ${worst.toExponential(2)} (fast path)`);
    expect(worst).toBeLessThan(1e-5);
  });

  it('dV/dr = frontArea (exact path)', () => {
    const e = newFlameGeometryResult();
    for (const [r, h] of [
      [0.003, 0.02],
      [0.03, 0.0005],
      [0.05, 0.019],
      [0.0784, 0.012],
      [0.1, 0.1],
    ]) {
      const dr = 1e-6 * r;
      const vp = g.evaluateExact(r + dr, h, e).volume;
      const vm = g.evaluateExact(r - dr, h, e).volume;
      const A = g.evaluateExact(r, h, e).frontArea;
      expect(Math.abs((vp - vm) / (2 * dr) - A) / A).toBeLessThan(1e-7);
    }
  });

  it('radiusForVolume inverts volume (round trip to 1e-9) across regimes', () => {
    const R = CFR_BORE / 2;
    for (const h of [0.0003, 0.0067, 0.019, 0.06, 0.15]) {
      const Vch = PI * R * R * h;
      for (const frac of [1e-9, 1e-6, 1e-3, 0.05, 0.3, 0.7, 0.95, 0.9999]) {
        const V = frac * Vch;
        const r = g.radiusForVolume(V, h);
        // (1e-9: for h < 2 mm the spark centre lies below the piston and a tiny V is a difference of
        // two half-ball integrals — cancellation floor ~5e-10 at V = 1e-9 V_chamber)
        expect(Math.abs(g.volume(r, h) - V) / V).toBeLessThan(1e-9);
      }
      expect(g.radiusForVolume(Vch, h)).toBeCloseTo(g.maxRadius(h), 12);
      expect(g.radiusForVolume(0, h)).toBe(0);
    }
  });

  it('monotone in r, continuous across the liner-touch radius and the fallback bands', () => {
    const out = newFlameGeometryResult();
    const R = CFR_BORE / 2;
    const c1 = R - CFR_CENTER[0];
    const c2 = R + CFR_CENTER[0];
    const h = 0.019;
    for (const rc of [c1, c1 + 0.00015, c2 - 0.01 * R, c2, c2 + 0.01 * R]) {
      const a: FlameGeometryResult = { ...g.evaluate(rc * (1 - 1e-9), h, out) };
      const bb = g.evaluate(rc * (1 + 1e-9), h, out);
      for (const k of KEYS) {
        const scale = Math.max(Math.abs(a[k]), 1e-6 * PI * R * R);
        // table ↔ exact switch at the band edges: jump = table error there (~1e-4 of A_f)
        expect(Math.abs(bb[k] - a[k]) / scale, `${k} at ${rc}`).toBeLessThan(5e-4);
      }
    }
    let prev = 0;
    for (let i = 1; i <= 2000; i++) {
      const v = g.volume((i / 2000) * g.maxRadius(h), h);
      expect(v).toBeGreaterThanOrEqual(prev - 1e-18);
      prev = v;
    }
  });

  it('evaluation cost (reported)', () => {
    const out = newFlameGeometryResult();
    const N = 200000;
    let s = 0;
    for (let i = 0; i < 20000; i++) s += g.evaluate(0.002 + (i % 2000) * 4e-5, 0.008 + (i % 97) * 3e-4, out).volume;
    let t0 = performance.now();
    for (let i = 0; i < N; i++) {
      const x = i / N; // combustion-like sweep: r grows, h grows
      s += g.evaluate(0.002 + 0.09 * x, 0.0067 + 0.03 * x * x, out).volume;
    }
    const nsEval = ((performance.now() - t0) / N) * 1e6;
    t0 = performance.now();
    const R = CFR_BORE / 2;
    for (let i = 0; i < N / 4; i++) {
      const x = (4 * i) / N;
      const h = 0.0067 + 0.03 * x * x;
      s += g.radiusForVolume(x * PI * R * R * h, h);
    }
    const nsInv = ((performance.now() - t0) / (N / 4)) * 1e6;
    t0 = performance.now();
    const gNew = new FlameGeometry(CFR_BORE, CFR_CENTER, { maxHeight: 0.16 });
    const msBuild = performance.now() - t0;
    console.log(
      `[flame-geometry] evaluate ${nsEval.toFixed(0)} ns/call, radiusForVolume ${nsInv.toFixed(0)} ns/call, ` +
        `table build ${msBuild.toFixed(0)} ms (${s > 0 && gNew ? 'ok' : ''})`,
    );
    expect(nsEval).toBeLessThan(5000);
  });

  it('rejects invalid spark positions', () => {
    expect(() => new FlameGeometry(0.08, [0.05, -0.001, 0])).toThrow();
    expect(() => new FlameGeometry(0.08, [0, 0.001, 0])).toThrow();
  });
});
