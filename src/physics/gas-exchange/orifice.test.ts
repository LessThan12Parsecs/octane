import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_orifice.json';
import { NS, SP, type SpeciesName } from '../core/species';
import { mixGamma, mixGasConstant } from '../thermo';
import {
  chokedFluxFunction,
  criticalPressureRatio,
  isentropicJetVelocity,
  massFluxFunction,
  newOrificeFlow,
  ORIFICE_REG_DEFAULT,
  ORIFICE_REG_MAX_DEFICIT,
  orificeFlow,
  orificeMassFlow,
  regularizedFluxFunction,
} from './orifice';

interface Analytic {
  gamma: number;
  rCrit: number;
  psiCrit: number;
  rArgmaxNumeric: number;
  psiMaxNumeric: number;
  points: { r: number; psi: number }[];
}
interface VarCp {
  label: string;
  T0: number;
  p0: number;
  X: Record<string, number>;
  R: number;
  gammaT0: number;
  rStarVariableCp: number;
  rows: { r: number; massFlux: number; jetVelocity: number; throatT: number; choked: boolean }[];
}
const F = fixture as unknown as { analytic: Analytic[]; variableCp: VarCp[] };

const rel = (a: number, b: number) => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

describe('isentropic restriction flow vs 50-digit analytic oracle', () => {
  for (const c of F.analytic) {
    it(`γ = ${c.gamma}: r*, Ψ*, Ψ(r)`, () => {
      expect(rel(criticalPressureRatio(c.gamma), c.rCrit)).toBeLessThan(1e-14);
      expect(rel(chokedFluxFunction(c.gamma), c.psiCrit)).toBeLessThan(1e-14);
      // Ψ is maximal at r* (numerical maximisation of the subsonic formula)
      expect(Math.abs(c.rArgmaxNumeric - c.rCrit)).toBeLessThan(1e-7);
      expect(rel(c.psiMaxNumeric, c.psiCrit)).toBeLessThan(1e-13);
      for (const p of c.points) {
        const psi = massFluxFunction(p.r, c.gamma);
        if (p.psi === 0) expect(psi).toBe(0);
        else expect(rel(psi, p.psi)).toBeLessThan(2e-13);
      }
    });
  }

  it('regularised flux is exact outside the band, bounded inside, C² at the band edge', () => {
    const g = 1.4;
    const d = ORIFICE_REG_DEFAULT;
    for (const x of [d, 1.5 * d, 0.01, 0.1, 0.3, 0.47, 0.6, 1]) {
      expect(rel(regularizedFluxFunction(x, g, d), massFluxFunction(1 - x, g))).toBeLessThan(1e-13);
    }
    // in-band deficit ≤ ORIFICE_REG_MAX_DEFICIT · Ψ(δ)
    const psiD = massFluxFunction(1 - d, g);
    let maxDef = 0;
    for (let i = 1; i <= 2000; i++) {
      const x = (i / 2000) * d;
      const def = massFluxFunction(1 - x, g) - regularizedFluxFunction(x, g, d);
      expect(def).toBeGreaterThanOrEqual(-1e-15);
      maxDef = Math.max(maxDef, def / psiD);
    }
    expect(maxDef).toBeLessThan(ORIFICE_REG_MAX_DEFICIT * (1 + 2e-3));
    expect(maxDef).toBeGreaterThan(ORIFICE_REG_MAX_DEFICIT * 0.99);
    // value, slope and curvature continuity at x = δ (finite differences)
    const h = d * 1e-5; // one-sided 2nd differences: O(h·f‴) bias, f‴ jumps at δ (C² only)
    const f = (x: number) => regularizedFluxFunction(x, g, d);
    const exact = (x: number) => massFluxFunction(1 - x, g);
    expect(rel(f(d - h), exact(d - h))).toBeLessThan(1e-7); // O(h³) mismatch only
    const d2in = (f(d) - 2 * f(d - h) + f(d - 2 * h)) / (h * h);
    const d2out = (exact(d + 2 * h) - 2 * exact(d + h) + exact(d)) / (h * h);
    expect(Math.abs(d2in - d2out) / Math.abs(d2out)).toBeLessThan(2e-3);
    // finite slope at 0, odd-symmetric extension via the bidirectional API
    const s0 = (f(1e-9) - 0) / 1e-9;
    expect(s0).toBeCloseTo((Math.SQRT2 * 1.40625) / Math.sqrt(d), 3);
  });

  it('orificeMassFlow: incompressible limit, choked limit, no NaN at zero area / Δp / reverse', () => {
    const R = 287;
    const T = 300;
    const p0 = 1e5;
    const A = 1e-4;
    // small but out-of-band drop: ṁ → A √(2ρΔp) (1 + O(x))
    const dp = 50; // x = 5e-4
    const m = orificeMassFlow(A, p0, T, R, 1.4, p0 - dp, 1e-5);
    const rho = p0 / (R * T);
    expect(rel(m, A * Math.sqrt(2 * rho * dp))).toBeLessThan(5e-4);
    // choked
    const mc = orificeMassFlow(A, p0, T, R, 1.4, 0.3 * p0);
    expect(rel(mc, (A * p0 * chokedFluxFunction(1.4)) / Math.sqrt(R * T))).toBeLessThan(1e-14);
    expect(orificeMassFlow(A, p0, T, R, 1.4, 0)).toBe(mc);
    expect(orificeMassFlow(A, p0, T, R, 1.4, p0)).toBe(0);
    expect(orificeMassFlow(A, p0, T, R, 1.4, 2 * p0)).toBe(0);
    expect(orificeMassFlow(0, p0, T, R, 1.4, 0.5 * p0)).toBe(0);
    const out = newOrificeFlow();
    for (const [pa, pb, CdA] of [
      [1e5, 1e5, 1e-4],
      [1e5, 1e5 + 1e-9, 1e-4],
      [1e5, 2e5, 1e-4],
      [1e5, 0.2e5, 0],
      [0, 0, 1e-4],
    ]) {
      orificeFlow(CdA, pa, 300, 287, 1.4, pb, 900, 290, 1.3, out);
      for (const v of Object.values(out)) expect(Number.isFinite(Number(v))).toBe(true);
    }
  });

  it('orificeFlow is antisymmetric, uses the upstream state, and its pressure derivatives are exact', () => {
    const out = newOrificeFlow();
    const a = { p: 1.02e5, T: 320, R: 288, g: 1.39 };
    const b = { p: 0.99e5, T: 1100, R: 290, g: 1.31 };
    const A = 2e-4;
    const m1 = orificeFlow(A, a.p, a.T, a.R, a.g, b.p, b.T, b.R, b.g, out);
    expect(out.upstream).toBe(0);
    expect(m1).toBeCloseTo(orificeMassFlow(A, a.p, a.T, a.R, a.g, b.p), 15);
    const m2 = orificeFlow(A, b.p, b.T, b.R, b.g, a.p, a.T, a.R, a.g, out);
    expect(out.upstream).toBe(1);
    expect(m2).toBe(-m1);
    // jet speed reported by orificeFlow = isentropicJetVelocity of the upstream state
    for (const pd of [0.99e5, 0.5e5, 1.019e5, 1.0199999e5]) {
      orificeFlow(A, a.p, a.T, a.R, a.g, pd, b.T, b.R, b.g, out);
      expect(rel(out.jetVelocity, isentropicJetVelocity(a.T, a.R, a.g, pd / a.p))).toBeLessThan(1e-9);
    }
    // derivatives vs central differences, in and out of the band, subsonic and choked
    for (const [pa, pb] of [
      [1.02e5, 0.99e5],
      [1.00005e5, 1e5],
      [0.99997e5, 1e5],
      [3e5, 1e5],
      [1e5, 1.5e5],
    ]) {
      orificeFlow(A, pa, a.T, a.R, a.g, pb, b.T, b.R, b.g, out);
      const dA = out.dmdotdpa;
      const dB = out.dmdotdpb;
      const h = 1e-3;
      const f = (x: number, y: number) => orificeFlow(A, x, a.T, a.R, a.g, y, b.T, b.R, b.g, newOrificeFlow());
      const fdA = (f(pa + h, pb) - f(pa - h, pb)) / (2 * h);
      const fdB = (f(pa, pb + h) - f(pa, pb - h)) / (2 * h);
      expect(Math.abs(dA - fdA)).toBeLessThan(1e-6 * Math.abs(fdA) + 1e-14);
      expect(Math.abs(dB - fdB)).toBeLessThan(1e-6 * Math.abs(fdB) + 1e-14);
    }
    // continuity through Δp = 0
    const eps = 1e-6;
    expect(Math.abs(orificeFlow(A, 1e5 + eps, 300, 287, 1.4, 1e5, 900, 290, 1.3, out))).toBeLessThan(1e-9);
    expect(Math.abs(orificeFlow(A, 1e5 - eps, 300, 287, 1.4, 1e5, 900, 290, 1.3, out))).toBeLessThan(1e-9);
  });

  it('jet velocity: isentropic, sonic when choked', () => {
    const g = 1.4;
    const R = 287;
    const T0 = 300;
    const rc = criticalPressureRatio(g);
    const vStar = Math.sqrt(((g * R * T0) * 2) / (g + 1)); // c* = √(γ R T*) , T* = 2T0/(γ+1)
    expect(rel(isentropicJetVelocity(T0, R, g, 0.1), vStar)).toBeLessThan(1e-13);
    expect(rel(isentropicJetVelocity(T0, R, g, rc), vStar)).toBeLessThan(1e-13);
    expect(isentropicJetVelocity(T0, R, g, 1)).toBe(0);
    // small drop: v → √(2Δp/ρ)
    const x = 1e-6;
    expect(rel(isentropicJetVelocity(T0, R, g, 1 - x), Math.sqrt(2 * x * R * T0))).toBeLessThan(1e-6);
  });
});

describe('constant-γ restriction flow vs Cantera variable-cp isentropic nozzle', () => {
  for (const c of F.variableCp) {
    it(`${c.label}: frozen-γ(T0) error bounded`, () => {
      const X = new Float64Array(NS);
      for (const [k, v] of Object.entries(c.X)) X[SP[k as SpeciesName]] = v;
      let s = 0;
      for (let k = 0; k < NS; k++) s += X[k];
      for (let k = 0; k < NS; k++) X[k] /= s;
      const g = mixGamma(X, c.T0);
      const R = mixGasConstant(X);
      expect(rel(g, c.gammaT0)).toBeLessThan(1e-9);
      expect(rel(R, c.R)).toBeLessThan(1e-9);
      let maxErr = 0;
      let maxVErr = 0;
      for (const row of c.rows) {
        // band δ = 1e-6 < every tested drop (r ≤ 0.9999): this measures the frozen-γ error only
        const flux = orificeMassFlow(1, c.p0, c.T0, R, g, row.r * c.p0, 1e-6);
        maxErr = Math.max(maxErr, rel(flux, row.massFlux));
        maxVErr = Math.max(maxVErr, rel(isentropicJetVelocity(c.T0, R, g, row.r), row.jetVelocity));
      }
      // Frozen-γ error: air 300 K < 0.05 %; fuel-vapour charges < 0.2 % (the fuel's cp varies
      // strongly with T: measured 0.11–0.15 %, jet speed 0.18–0.25 %); burned gas 1200–2200 K
      // < 0.5 % (measured ≤ 0.15 %, jet speed ≤ 0.25 %).
      const tol = c.label.startsWith('air') ? 5e-4 : c.label.startsWith('fresh') ? 2e-3 : 5e-3;
      expect(maxErr).toBeLessThan(tol);
      expect(maxVErr).toBeLessThan(c.label.startsWith('fresh') ? 3e-3 : tol);
      console.info(`[orifice] ${c.label}: max |Δṁ|/ṁ = ${maxErr.toExponential(2)}, |Δv|/v = ${maxVErr.toExponential(2)}`);
    });
  }
});
