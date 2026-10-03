import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_plenum.json';
import { NS, SP, type SpeciesName } from '../core/species';
import { MOLAR_MASS } from '../thermo';
import { newOrificeFlow, orificeFlow } from './orifice';
import { gasStateFromNU, gasStateFromTPX, newGasState, Plenum, type GasState } from './plenum';

type Comp = Record<string, number>;
interface Fx {
  compositions: { air: Comp; fresh: Comp; burned: Comp };
  stateRecovery: { label: string; V: number; N: Comp; U: number; T: number; p: number }[];
  filling: {
    label: string;
    V: number;
    tank: { T: number; p: number; X: Comp };
    reservoir: { T: number; p: number; X: Comp; h0: number };
    admittedMass: number;
    T2: number;
    p2: number;
  }[];
  mixing: {
    A: { T: number; p: number; V: number };
    B: { T: number; p: number; V: number };
    Ntotal: Comp;
    Utotal: number;
    mixedT: number;
    mixedP: number;
  };
}
const F = fixture as unknown as Fx;

const vec = (c: Comp): Float64Array => {
  const v = new Float64Array(NS);
  for (const [k, x] of Object.entries(c)) v[SP[k as SpeciesName]] = x;
  return v;
};
const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

/** Classic RK4 on a packed state [N(NS), U] of several plenums. */
function rk4(y: Float64Array, dt: number, f: (y: Float64Array, dy: Float64Array) => void, s: Float64Array[]): void {
  const [k1, k2, k3, k4, t] = s;
  const n = y.length;
  f(y, k1);
  for (let i = 0; i < n; i++) t[i] = y[i] + 0.5 * dt * k1[i];
  f(t, k2);
  for (let i = 0; i < n; i++) t[i] = y[i] + 0.5 * dt * k2[i];
  f(t, k3);
  for (let i = 0; i < n; i++) t[i] = y[i] + dt * k3[i];
  f(t, k4);
  for (let i = 0; i < n; i++) y[i] += (dt / 6) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
}
const scratch = (n: number) => Array.from({ length: 5 }, () => new Float64Array(n));

describe('gas state recovery vs Cantera UV', () => {
  for (const c of F.stateRecovery) {
    it(`${c.label}: T, p from (N, U, V)`, () => {
      const s = gasStateFromNU(vec(c.N), c.U, c.V, 500, newGasState());
      expect(rel(s.T, c.T)).toBeLessThan(1e-9);
      expect(rel(s.p, c.p)).toBeLessThan(1e-9);
      let sy = 0;
      for (let k = 0; k < NS; k++) sy += s.Y[k];
      expect(Math.abs(sy - 1)).toBeLessThan(1e-14);
      expect(rel(s.rho, s.p / (s.R * s.T))).toBeLessThan(1e-15);
    });
  }

  it('Plenum.setTPX ↔ updateState round trip', () => {
    const pl = new Plenum(1e-3, 345.6, 0.87e5, vec(F.compositions.fresh));
    pl.state.T = 900; // bad warm start
    const s = pl.updateState();
    expect(rel(s.T, 345.6)).toBeLessThan(1e-11);
    expect(rel(s.p, 0.87e5)).toBeLessThan(1e-11);
    expect(rel(pl.mass(), s.rho * pl.volume)).toBeLessThan(1e-12);
  });
});

describe('plenum filling and emptying', () => {
  for (const c of F.filling) {
    it(`${c.label}: adiabatic filling from a reservoir reaches the Cantera end state; energy exact`, () => {
      const V = c.V;
      const tank = new Plenum(V, c.tank.T, c.tank.p, vec(c.tank.X));
      const res = gasStateFromTPX(c.reservoir.T, c.reservoir.p, vec(c.reservoir.X), newGasState());
      expect(rel(res.h, c.reservoir.h0)).toBeLessThan(1e-10);
      const CdA = 1e-5;
      const fl = newOrificeFlow();
      // packed state: [N(NS), U, admitted mass]
      const n = NS + 2;
      const y = new Float64Array(n);
      y.set(tank.N, 0);
      y[NS] = tank.U;
      const m0 = tank.mass();
      const U0 = tank.U;
      const rhs = (yy: Float64Array, dy: Float64Array) => {
        tank.N.set(yy.subarray(0, NS));
        tank.U = yy[NS];
        const s = tank.updateState();
        const mdot = orificeFlow(CdA, res.p, res.T, res.R, res.gamma, s.p, s.T, s.R, s.gamma, fl);
        tank.beginRates();
        tank.addFlow(mdot, res);
        dy.set(tank.dNdt, 0);
        dy[NS] = tank.dUdt;
        dy[NS + 1] = mdot;
      };
      const sc = scratch(n);
      let t = 0;
      const dt = 1e-4;
      for (; t < 5; t += dt) {
        rk4(y, dt, rhs, sc);
        tank.N.set(y.subarray(0, NS));
        tank.U = y[NS];
        const s = tank.updateState();
        if (Math.abs(s.p - res.p) / res.p < 1e-11) break;
        // no NaN anywhere
        expect(Number.isFinite(s.T)).toBe(true);
      }
      const s = tank.updateState();
      expect(rel(s.p, c.p2)).toBeLessThan(1e-10);
      expect(rel(s.T, c.T2)).toBeLessThan(1e-8);
      expect(rel(tank.mass() - m0, c.admittedMass)).toBeLessThan(1e-7);
      // energy: U − U0 = h0 · admitted mass (exact for any integrator)
      expect(Math.abs(tank.U - U0 - res.h * y[NS + 1]) / Math.abs(res.h * y[NS + 1])).toBeLessThan(1e-11);
      // species: every admitted kg carries the reservoir mass fractions
      for (let k = 0; k < NS; k++) {
        const expected = (y[NS + 1] * res.Y[k]) / MOLAR_MASS[k];
        const got = tank.N[k] - vec(c.tank.X)[k] * ((c.tank.p * V) / (8.31446261815324 * c.tank.T)) / sumC(c.tank.X);
        expect(Math.abs(got - expected)).toBeLessThan(1e-12 * (1 + Math.abs(expected)) + 1e-18);
      }
      if (c.label === 'air-into-air') {
        // constant-cp textbook result: T2 = p2 / (p1/T1 + (p2 − p1)/(γ T0))  (≈ 350.0 K here);
        // the variable-cp Cantera answer differs by ~0.03 %.
        const g = res.gamma;
        const T2cp = c.p2 / (c.tank.p / c.tank.T + (c.p2 - c.tank.p) / (g * c.reservoir.T));
        expect(rel(s.T, T2cp)).toBeLessThan(1e-3);
      }
    });
  }

  it('two plenums equalising through an orifice conserve species moles and energy to rounding', () => {
    const mx = F.mixing;
    const A = new Plenum(mx.A.V, mx.A.T, mx.A.p, vec(F.compositions.burned));
    const B = new Plenum(mx.B.V, mx.B.T, mx.B.p, vec(F.compositions.fresh));
    const Ntot = vec(mx.Ntotal);
    for (let k = 0; k < NS; k++) {
      expect(Math.abs(A.N[k] + B.N[k] - Ntot[k])).toBeLessThan(1e-13 * (Ntot[k] + 1e-6));
    }
    expect(rel(A.U + B.U, mx.Utotal)).toBeLessThan(1e-12);
    const CdA = 3e-6;
    const fl = newOrificeFlow();
    const n = 2 * (NS + 1);
    const y = new Float64Array(n);
    y.set(A.N, 0);
    y[NS] = A.U;
    y.set(B.N, NS + 1);
    y[2 * NS + 1] = B.U;
    const load = (yy: Float64Array) => {
      A.N.set(yy.subarray(0, NS));
      A.U = yy[NS];
      B.N.set(yy.subarray(NS + 1, 2 * NS + 1));
      B.U = yy[2 * NS + 1];
      A.updateState();
      B.updateState();
    };
    const rhs = (yy: Float64Array, dy: Float64Array) => {
      load(yy);
      const sa: GasState = A.state;
      const sb: GasState = B.state;
      const m = orificeFlow(CdA, sa.p, sa.T, sa.R, sa.gamma, sb.p, sb.T, sb.R, sb.gamma, fl); // A → B
      A.beginRates();
      B.beginRates();
      A.addFlow(-m, sb);
      B.addFlow(m, sa);
      dy.set(A.dNdt, 0);
      dy[NS] = A.dUdt;
      dy.set(B.dNdt, NS + 1);
      dy[2 * NS + 1] = B.dUdt;
    };
    const sc = scratch(n);
    for (let i = 0; i < 3000; i++) rk4(y, 1e-4, rhs, sc);
    load(y);
    expect(Math.abs(A.state.p - B.state.p) / B.state.p).toBeLessThan(1e-6);
    for (let k = 0; k < NS; k++) {
      expect(Math.abs(A.N[k] + B.N[k] - Ntot[k])).toBeLessThan(1e-13 * (Ntot[k] + 1e-3));
    }
    expect(Math.abs(A.U + B.U - mx.Utotal)).toBeLessThan(1e-12 * Math.abs(mx.Utotal) + 1e-10);
    // mass-weighted sanity: final pressure between the fully mixed single-volume pressure bounds
    expect(A.state.p).toBeGreaterThan(1e5);
    expect(A.state.p).toBeLessThan(3e5);
  });

  it('beginRates adds heat and boundary work; addFlow picks the upstream state by sign', () => {
    const pl = new Plenum(1e-3, 300, 1e5, vec(F.compositions.air));
    const other = gasStateFromTPX(800, 2e5, vec(F.compositions.burned), newGasState());
    pl.beginRates(10, 1e-4);
    expect(pl.dUdt).toBeCloseTo(10 - 1e5 * 1e-4, 9);
    pl.beginRates();
    pl.addFlow(0.01, other);
    expect(pl.dUdt).toBeCloseTo(0.01 * other.h, 9);
    pl.beginRates();
    pl.addFlow(-0.01, other);
    expect(pl.dUdt).toBeCloseTo(-0.01 * pl.state.h, 9);
    let dm = 0;
    for (let k = 0; k < NS; k++) dm += pl.dNdt[k] * MOLAR_MASS[k];
    expect(dm).toBeCloseTo(-0.01, 15);
    pl.beginRates();
    pl.addFlow(0, other);
    expect(pl.dUdt).toBe(0);
  });
});

function sumC(c: Comp): number {
  let s = 0;
  for (const v of Object.values(c)) s += v;
  return s;
}
