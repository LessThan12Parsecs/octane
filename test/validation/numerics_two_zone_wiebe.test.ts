/**
 * Validation round 2 (numerics): the TWO-ZONE path of the cycle integrator against an independent
 * Cantera oracle (tools/reference/validation_numerics_two_zone_wiebe.py →
 * test/fixtures/validation_numerics_two_zone_wiebe.json).
 *
 * The fuel–air-cycle oracles only exercise one zone (frozen compression, UV burn at TDC, burned-only
 * expansion). Here the charge burns along a prescribed Wiebe curve (combustionModel 'wiebe', no wall
 * heat, no knock), so the whole run between the Wiebe start and burn-out is two-zone: frozen isentropic
 * unburned zone + one fully mixed equilibrium burned zone at common p, U_tot integrated with −p dV and
 * mass moved between the zones at constant U_tot (closure.ts + cycle-model.ts RK4). The oracle solves
 * the same continuous formulation with Cantera thermo, an element-potential-polished TP equilibrium
 * (Cantera's own equilibrate reproduces u_b only to ~1e-9 of c_v T because it lets the element ratios
 * drift at that level) and DOP853 at rtol 1e-12.
 *
 * Checks: p, T_u, T_b at every crank degree, the work integral, U at start/end (absolute NASA
 * reference), and the RK4 step-size order on the two-zone phase. Runtime ≈ 3–5 s.
 */
import { describe, expect, it } from 'vitest';
import fixture from '../fixtures/validation_numerics_two_zone_wiebe.json';
import type { EngineSpec } from '../../src/physics/core/engine-spec';
import type { FuelSelection } from '../../src/physics/core/operating-point';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { MODE_OPEN, MODE_SINGLE, runClosedCycle } from '../../src/physics/cycle/index';

type Case = (typeof fixture.cases)[number];

function run(c: Case, maxStepDeg = 0.25) {
  const inp = c.input;
  const spec: EngineSpec = {
    ...CFR_F1,
    geometry: { ...CFR_F1.geometry, bore: inp.bore, stroke: inp.stroke, conRodLength: inp.rod, pinOffset: 0, compressionRatioRange: [2, 20] },
  };
  const fuels = Object.keys(inp.fuel_X);
  const fuel: FuelSelection = fuels.length > 1 ? { kind: 'PRF', octaneNumber: 90 } : { kind: 'pure', species: fuels[0] as 'IC8H18' };
  const op = { ...CFR_RON_CONDITIONS, compressionRatio: inp.cr, fuel, equivalenceRatio: inp.phi };
  return runClosedCycle(
    spec,
    op,
    { startDeg: inp.startDeg, endDeg: inp.endDeg, T: inp.T1, p: inp.p1, X: c.X_reactants as unknown as Readonly<Record<string, number>> },
    {
      heatTransfer: false,
      knock: false,
      combustionModel: 'wiebe',
      wiebe: { startDeg: inp.wiebe.startDeg, durationDeg: inp.wiebe.durationDeg, a: inp.wiebe.a, m: inp.wiebe.m },
      burnoutFraction: fixture.merge_fraction,
      maxStepDeg,
      fineStepDeg: maxStepDeg,
    },
  );
}

interface Err {
  p: number;
  Tu: number;
  Tb: number;
  at: string;
}

function compare(c: Case, maxStepDeg = 0.25): { e: Err; work: number; U0: number; U1: number } {
  const r = run(c, maxStepDeg);
  const tr = r.trace;
  const e: Err = { p: 0, Tu: 0, Tb: 0, at: '' };
  for (const row of c.trace) {
    // the LAST closed-phase sample at this angle (post-event), except at the Wiebe start where the
    // oracle row is the pre-seed single-zone state; at the end angle the EVO merge adds an open-phase
    // sample (kinetic-NO swap) that is not part of the closed cycle
    let i = -1;
    for (let k = 0; k < tr.theta.length; k++) if (Math.abs(tr.theta[k] - row.deg) < 1e-7 && tr.mode[k] !== MODE_OPEN) i = k;
    if (i < 0) throw new Error(`no trace sample at ${row.deg}°`);
    if (row.Tb === 0 && row.Tu > 0) while (i > 0 && tr.mode[i] !== MODE_SINGLE && Math.abs(tr.theta[i - 1] - row.deg) < 1e-7) i--;
    const ep = Math.abs(tr.pressure[i] / row.p - 1);
    if (ep > e.p) {
      e.p = ep;
      e.at = `p @${row.deg}°`;
    }
    // T_u while the unburned zone holds ≥ 1e-6 of the charge (the model merges the rest at a step end)
    if (row.Tu > 0 && row.xb <= 1 - 1e-6) e.Tu = Math.max(e.Tu, Math.abs(tr.Tu[i] / row.Tu - 1));
    // T_b only where the burned zone holds ≥ 1e-4 of the charge (a 1e-8 seed has no meaningful T_b)
    if (row.Tb > 0 && row.xb >= 1e-4) e.Tb = Math.max(e.Tb, Math.abs(tr.Tb[i] / row.Tb - 1));
  }
  return { e, work: r.work / c.work - 1, U0: (r.U0 - c.U_start) / Math.abs(c.work), U1: (r.U1 - c.U_end) / Math.abs(c.work) };
}

describe('two-zone Wiebe closed cycle vs the Cantera two-zone oracle', () => {
  for (const c of fixture.cases) {
    it(c.input.label, () => {
      const { e, work, U0, U1 } = compare(c);
      console.log(
        `[validation] two-zone Wiebe ${c.input.label}: max rel. error p ${e.p.toExponential(2)} (${e.at}), T_u ${e.Tu.toExponential(2)}, T_b ${e.Tb.toExponential(2)}, ` +
          `work ${work.toExponential(2)}, U0 ${U0.toExponential(2)}, U_end ${U1.toExponential(2)} (of W)`,
      );
      expect(e.p).toBeLessThan(2e-6);
      expect(e.Tu).toBeLessThan(2e-6);
      expect(e.Tb).toBeLessThan(2e-6);
      expect(Math.abs(work)).toBeLessThan(2e-6);
      expect(Math.abs(U0)).toBeLessThan(1e-8);
      expect(Math.abs(U1)).toBeLessThan(2e-6);
    });
  }

  it('RK4 step-size order on the two-zone phase (very fast burn, 10° Wiebe)', () => {
    const c = fixture.cases.find((x) => x.input.wiebe.durationDeg === 10)!;
    const r = [1, 0.5, 0.25, 0.125].map((d) => compare(c, d));
    const ew = r.map((x) => Math.abs(x.work));
    const ep = r.map((x) => x.e.p);
    const ord = (a: number[]) => a.slice(1).map((x, i) => Math.log2(a[i] / x).toFixed(2));
    console.log(
      `[validation] two-zone Wiebe step order (maxStepDeg 1/0.5/0.25/0.125): work error ${ew.map((x) => x.toExponential(2)).join(' / ')} (orders ${ord(ew).join(', ')}); ` +
        `max p error ${ep.map((x) => x.toExponential(2)).join(' / ')} (orders ${ord(ep).join(', ')})`,
    );
    // at default steps (0.25°) the fast burn is resolved to ≤ 1e-6 in work and pressure
    expect(ew[2]).toBeLessThan(1e-6);
    expect(ep[2]).toBeLessThan(2e-6);
  });
});
