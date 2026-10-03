import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/cycle_two_zone_closure.json';
import { NS, SP } from '../core/species';
import { ZoneClosure } from './closure';

/**
 * Two-zone closure vs the Cantera oracle (tools/reference/cycle_two_zone_closure.py): from the
 * conserved (U, S_u, m_u, m_b, V) the Newton closure must recover (p, T_u, T_b), the zone
 * volumes and the burned equilibrium composition — started from deliberately poor guesses.
 */
describe('ZoneClosure vs Cantera two-zone states', () => {
  for (const c of fixture.cases) {
    it(c.label, () => {
      const cl = new ZoneClosure();
      cl.setUnburned(Float64Array.from(c.Xu));
      let ok: boolean;
      if (c.mu > 0) {
        ok = cl.solveTwoZone(c.U, c.S, c.mu, c.mb, c.V, 0.7 * c.p, 0.85 * c.Tb, 0.9 * c.Tu);
        expect(cl.Tu / c.Tu - 1).toBeCloseTo(0, 9);
        expect((c.mu * cl.vu) / c.Vu - 1).toBeCloseTo(0, 8);
      } else {
        ok = cl.solveBurnedOnly(c.U, c.mb, c.V, 0.7 * c.p, 0.85 * c.Tb);
      }
      expect(ok).toBe(true);
      expect(cl.p / c.p - 1).toBeCloseTo(0, 9);
      expect(cl.Tb / c.Tb - 1).toBeCloseTo(0, 8);
      expect((c.mb * cl.vb) / c.Vb - 1).toBeCloseTo(0, 7);
      const X = cl.eq.result.X;
      for (const s of ['CO2', 'H2O', 'CO', 'O2', 'OH', 'NO'] as const) {
        const ref = c.Xb[SP[s]];
        expect(Math.abs(X[SP[s]] - ref)).toBeLessThan(1e-8 + 1e-6 * ref);
      }
      expect(cl.iterations).toBeLessThan(15);
    });
  }

  it('rates(): (ṗ, Ṫ_b) equal finite differences of the closure along the state rates', () => {
    const c = fixture.cases[3];
    const cl = new ZoneClosure();
    cl.setUnburned(Float64Array.from(c.Xu));
    // arbitrary state rates: compression, heat loss from the unburned zone, burning
    const dU = -1500;
    const Qu = 300;
    const dmb = 0.02;
    const dV = -2e-3;
    const solveAt = (e: number): [number, number] => {
      const mu = c.mu - dmb * e;
      cl.solveTwoZone(c.U + dU * e, c.S * (mu / c.mu) - (Qu / c.Tu) * e, mu, c.mb + dmb * e, c.V + dV * e, c.p, c.Tb, c.Tu);
      return [cl.p, cl.Tb];
    };
    const h = 1e-6;
    const [pp, tp] = solveAt(h);
    const [pm, tm] = solveAt(-h);
    solveAt(0);
    cl.rates(dU, -Qu / cl.Tu, -dmb, dmb, dV, c.mu);
    // S_u(e) above keeps s_u fixed except for the heat term, i.e. dS = s_u dm_u − Q̇_u/T_u
    expect(cl.dp / ((pp - pm) / (2 * h)) - 1).toBeCloseTo(0, 5);
    expect(cl.dTb / ((tp - tm) / (2 * h)) - 1).toBeCloseTo(0, 5);
  });

  it('single zone: T from U, p from the ideal-gas law', () => {
    const c = fixture.cases[0];
    const cl = new ZoneClosure();
    cl.setUnburned(Float64Array.from(c.Xu));
    const m = c.mu + c.mb;
    cl.unburnedAtT(c.Tu, c.p);
    const U = m * cl.uu;
    const V = m * cl.vu;
    const T2 = cl.solveSingle(U, m, V, 500);
    expect(T2 / c.Tu - 1).toBeCloseTo(0, 10);
    expect(cl.p / c.p - 1).toBeCloseTo(0, 10);
    expect(NS).toBe(17);
  });
});
