import { describe, expect, it } from 'vitest';
import { P_ATM, R_UNIVERSAL } from '../core/constants';
import { EL, N_EQ, NE, NS, SP } from '../core/species';
import { evalAllNondim, P_REF_THERMO } from '../thermo/thermo';
import tpFx from '../../../test/fixtures/equilibrium_tp.json';
import enFx from '../../../test/fixtures/equilibrium_energy.json';
import { EquilibriumSolver } from './solver';
import { elementBalanceError as elementError, stationarityError } from './diagnostics';
import { adiabaticFlameTemperature } from './flame';

/** Fixture targets (task spec): |ΔX| ≤ 1e-7 for X > 1e-3; |ΔX|/X ≤ 1e-4 for X ≥ 1e-12. */
const ABS_MAJOR = 1e-7;
const REL_MINOR = 1e-4;

interface XErr {
  absMajor: number;
  relMinor: number;
  nCompared: number;
}
const newErr = (): XErr => ({ absMajor: 0, relMinor: 0, nCompared: 0 });

function compareX(X: Float64Array, ref: readonly number[], err: XErr, skip: readonly number[] = []): void {
  for (let k = 0; k < N_EQ; k++) {
    if (skip.includes(k)) continue;
    const r = ref[k];
    if (r > 1e-3) err.absMajor = Math.max(err.absMajor, Math.abs(X[k] - r));
    if (r >= 1e-12) {
      err.relMinor = Math.max(err.relMinor, Math.abs(X[k] - r) / r);
      err.nCompared++;
    }
  }
}

const comps = tpFx.compositions.map((c) => ({ label: c.label, b: Float64Array.from(c.b) }));

describe('solveTP vs Cantera (tools/reference/equilibrium_oracle.py)', () => {
  it(`${tpFx.cases.length} states, 300–4000 K, 0.1–300 bar, 18 compositions: cold starts`, () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let worstBal = 0;
    let worstStat = 0;
    let maxIt = 0;
    for (const c of tpFx.cases) {
      s.reset();
      const b = comps[c.c].b;
      const r = s.solveTP(b, c.T, c.p);
      expect(r.converged, `${comps[c.c].label} T=${c.T} p=${c.p}`).toBe(true);
      compareX(r.X, c.X, err, (c as { ill?: number[] }).ill ?? []);
      worstBal = Math.max(worstBal, elementError(r, b));
      worstStat = Math.max(worstStat, stationarityError(r));
      maxIt = Math.max(maxIt, r.iterations);
      for (let k = N_EQ; k < NS; k++) expect(r.N[k]).toBe(0);
    }
    console.log(
      `TP cold: max|ΔX| (X>1e-3) = ${err.absMajor.toExponential(2)}, max rel (X≥1e-12) = ` +
        `${err.relMinor.toExponential(2)} over ${err.nCompared} values; element balance ` +
        `${worstBal.toExponential(2)}; stationarity ${worstStat.toExponential(2)}; max iterations ${maxIt}`,
    );
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
    expect(worstBal).toBeLessThan(1e-12);
    expect(worstStat).toBeLessThan(1e-8);
  });

  it('same states with warm starts carried across the whole grid', () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let its = 0;
    for (const c of tpFx.cases) {
      const r = s.solveTP(comps[c.c].b, c.T, c.p);
      expect(r.converged).toBe(true);
      compareX(r.X, c.X, err, (c as { ill?: number[] }).ill ?? []);
      its += r.iterations;
    }
    console.log(
      `TP warm grid: max|ΔX| = ${err.absMajor.toExponential(2)}, max rel = ${err.relMinor.toExponential(2)}, ` +
        `mean iterations ${(its / tpFx.cases.length).toFixed(2)}`,
    );
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
  });

  it('drops the species of absent elements (no C, no Ar, no N, no H)', () => {
    const s = new EquilibriumSolver();
    const byLabel = (l: string): Float64Array => comps.find((c) => c.label.startsWith(l))!.b;
    let r = s.solveTP(byLabel('H2-air'), 2500, 10e5);
    expect(r.converged).toBe(true);
    expect(r.N[SP.CO]).toBe(0);
    expect(r.N[SP.CO2]).toBe(0);
    expect(r.N[SP.AR]).toBe(0);
    expect(r.X[SP.NO]).toBeGreaterThan(1e-4);
    r = s.solveTP(byLabel('H2-O2'), 3000, 1e5);
    expect(r.converged).toBe(true);
    for (const k of [SP.N2, SP.NO, SP.N, SP.CO, SP.CO2, SP.AR]) expect(r.N[k]).toBe(0);
    r = s.solveTP(byLabel('CO-O2-N2'), 3000, 1e5);
    expect(r.converged).toBe(true);
    for (const k of [SP.H2O, SP.H2, SP.OH, SP.H, SP.AR]) expect(r.N[k]).toBe(0);
    expect(r.X[SP.CO]).toBeGreaterThan(0.01);
  });

  it('is invariant to the scale of b (extensive in N, intensive in X)', () => {
    const s = new EquilibriumSolver();
    const b = comps[0].b;
    const r1 = s.solveTP(b, 2600, 30e5);
    const X1 = Float64Array.from(r1.X);
    const N1 = Float64Array.from(r1.N);
    const b2 = b.map((v) => v * 3.7e-6);
    const r2 = s.solveTP(b2, 2600, 30e5);
    for (let k = 0; k < NS; k++) {
      expect(r2.X[k]).toBeCloseTo(X1[k], 14);
      expect(Math.abs(r2.N[k] - N1[k] * 3.7e-6)).toBeLessThanOrEqual(1e-11 * N1[k] * 3.7e-6 + 1e-300);
    }
    expect(r2.V).toBeCloseTo((r2.nTotal * R_UNIVERSAL * 2600) / 30e5, 18);
  });

  it('reports infeasible carbon-rich element sets (O ≤ C) without throwing', () => {
    const s = new EquilibriumSolver();
    const b = new Float64Array(NE);
    b[EL.C] = 1;
    b[EL.H] = 2;
    b[EL.O] = 0.9;
    b[EL.N] = 4;
    const r = s.solveTP(b, 2000, 1e5);
    expect(r.converged).toBe(false);
    // and the solver recovers on the next feasible problem
    b[EL.O] = 3;
    expect(s.solveTP(b, 2000, 1e5).converged).toBe(true);
  });

  it('result H, U, S, V, mass are consistent with N', () => {
    const s = new EquilibriumSolver();
    const r = s.solveTP(comps[12].b, 2300, 45e5);
    const cpR = new Float64Array(NS);
    const hRT = new Float64Array(NS);
    const s0R = new Float64Array(NS);
    evalAllNondim(r.T, cpR, hRT, s0R);
    let H = 0;
    let S = 0;
    for (let k = 0; k < N_EQ; k++) {
      H += r.N[k] * hRT[k] * R_UNIVERSAL * r.T;
      if (r.X[k] > 0) S += r.N[k] * R_UNIVERSAL * (s0R[k] - Math.log(r.X[k]) - Math.log(r.p / P_REF_THERMO));
    }
    expect(r.H).toBeCloseTo(H, 6);
    expect(r.U).toBeCloseTo(H - r.nTotal * R_UNIVERSAL * r.T, 6);
    expect(r.S).toBeCloseTo(S, 9);
    expect(r.mass).toBeCloseTo(1, 12); // b is per kg
  });
});

describe('solveTV / solveHP / solveUV / solveSP vs Cantera', () => {
  it('TV states (p within 1e-9, X within spec)', () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let worstP = 0;
    for (const c of enFx.tv) {
      const r = s.solveTV(Float64Array.from(c.b), c.T, c.V);
      expect(r.converged).toBe(true);
      compareX(r.X, c.X, err);
      worstP = Math.max(worstP, Math.abs(r.p / c.p - 1));
    }
    console.log(`TV: max|ΔX| ${err.absMajor.toExponential(2)}, rel ${err.relMinor.toExponential(2)}, p ${worstP.toExponential(2)}`);
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
    expect(worstP).toBeLessThan(1e-9);
  });

  it('HP adiabatic flames at 298 K/1 atm and 700 K/20 bar (T within 0.1 K)', () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let worstT = 0;
    let outer = 0;
    for (const c of enFx.hp) {
      const f = adiabaticFlameTemperature(Float64Array.from(c.Xr), c.T0, c.p0, 'HP', s);
      expect(f.converged, c.label).toBe(true);
      worstT = Math.max(worstT, Math.abs(f.T - c.T));
      compareX(f.X, c.X, err);
      outer = Math.max(outer, s.result.outerIterations);
    }
    console.log(`HP: max|ΔT| ${worstT.toExponential(2)} K, max|ΔX| ${err.absMajor.toExponential(2)}, rel ${err.relMinor.toExponential(2)}, max outer its ${outer}`);
    expect(worstT).toBeLessThan(0.1);
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
  });

  it('UV constant-volume combustion (T within 0.1 K, p within 1e-7)', () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let worstT = 0;
    let worstP = 0;
    for (const c of enFx.uv) {
      const f = adiabaticFlameTemperature(Float64Array.from(c.Xr), c.T0, c.p0, 'UV', s);
      expect(f.converged, c.label).toBe(true);
      worstT = Math.max(worstT, Math.abs(f.T - c.T));
      worstP = Math.max(worstP, Math.abs(f.p / c.p - 1));
      compareX(f.X, c.X, err);
    }
    console.log(`UV: max|ΔT| ${worstT.toExponential(2)} K, max rel Δp ${worstP.toExponential(2)}, max|ΔX| ${err.absMajor.toExponential(2)}, rel ${err.relMinor.toExponential(2)}`);
    expect(worstT).toBeLessThan(0.1);
    expect(worstP).toBeLessThan(1e-7);
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
  });

  it('SP isentropic equilibrium expansions (T within 0.1 K)', () => {
    const s = new EquilibriumSolver();
    const err = newErr();
    let worstT = 0;
    for (const c of enFx.sp) {
      const r = s.solveSP(Float64Array.from(c.b), c.S, c.p, 0.6 * c.T1);
      expect(r.converged, c.label).toBe(true);
      worstT = Math.max(worstT, Math.abs(r.T - c.T));
      compareX(r.X, c.X, err);
      expect(r.S).toBeCloseTo(c.S, 6);
    }
    console.log(`SP: max|ΔT| ${worstT.toExponential(2)} K, max|ΔX| ${err.absMajor.toExponential(2)}, rel ${err.relMinor.toExponential(2)}`);
    expect(worstT).toBeLessThan(0.1);
    expect(err.absMajor).toBeLessThanOrEqual(ABS_MAJOR);
    expect(err.relMinor).toBeLessThanOrEqual(REL_MINOR);
  });

  it('HP/UV/SP hit their targets and round-trip through TP', () => {
    const s = new EquilibriumSolver();
    const b = comps[11].b; // PRF80, humid air, φ = 0.9 (per kg)
    const r0 = s.solveTP(b, 2450, 55e5);
    const H = r0.H;
    const U = r0.U;
    const S = r0.S;
    const V = r0.V;
    for (const Tg of [600, 1500, 3900]) {
      s.reset();
      expect(s.solveHP(b, H, 55e5, Tg).T).toBeCloseTo(2450, 6);
      expect(Math.abs(s.result.H / H - 1)).toBeLessThan(1e-12); // review fix: was one-sided
      s.reset();
      const ru = s.solveUV(b, U, V, Tg);
      expect(ru.T).toBeCloseTo(2450, 6);
      expect(ru.p).toBeCloseTo(55e5, 0);
      s.reset();
      expect(s.solveSP(b, S, 55e5, Tg).T).toBeCloseTo(2450, 6);
    }
  });

  it('HP target outside the 200–6000 K bracket → converged = false at the bound', () => {
    const s = new EquilibriumSolver();
    const b = comps[0].b;
    const rHot = s.solveTP(b, 5990, P_ATM);
    const Hhot = rHot.H;
    const r = s.solveHP(b, Hhot * 1.5 + 5e6, P_ATM, 2000);
    expect(r.converged).toBe(false);
    expect(r.T).toBe(6000);
  });
});
