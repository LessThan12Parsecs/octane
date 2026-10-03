/**
 * Validation round 1 (numerics): extended fuel–air-cycle oracle
 * (tools/reference/validation_numerics_fuel_air.py → test/fixtures/validation_numerics_fuel_air.json).
 * Closed cycle −180 → 180°, no heat transfer, 'instantaneous-at-tdc' (UV equilibrium at TDC),
 * shifting-equilibrium expansion; rich/lean/very rich, PRF blend, n-heptane, ethanol, propane,
 * humid air, equilibrium residual, CR 4–18, 0.3–1.5 bar. p, T at 37 oracle points, work, efficiency
 * (oracle LHV basis). Also the RK4 step-size order on the same limit cycle (maxStepDeg 2 / 1 / 0.25).
 */
import { describe, expect, it } from 'vitest';
import fixture from '../fixtures/validation_numerics_fuel_air.json';
import type { EngineSpec } from '../../src/physics/core/engine-spec';
import type { FuelSelection } from '../../src/physics/core/operating-point';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { MODE_BURNED, runClosedCycle } from '../../src/physics/cycle/index';

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
    { startDeg: -180, endDeg: 180, T: inp.T1, p: inp.p1, X: c.X_reactants as unknown as Readonly<Record<string, number>> },
    { heatTransfer: false, combustionModel: 'instantaneous-at-tdc', knock: false, maxStepDeg },
  );
}

describe('closed cycle vs the extended Cantera fuel–air cycle', () => {
  const tol = 1e-6;
  let worst = 0;
  let worstLabel = '';
  for (const c of fixture.cases) {
    it(c.input.label, () => {
      const r = run(c);
      const tr = r.trace;
      expect(Math.abs(r.mass / c.mass - 1)).toBeLessThan(1e-9);
      const find = (deg: number, burned: boolean): number => {
        for (let i = 0; i < tr.theta.length; i++) if (Math.abs(tr.theta[i] - deg) < 1e-7 && (tr.mode[i] === MODE_BURNED) === burned) return i;
        throw new Error(`no trace sample at ${deg}°`);
      };
      const errs: number[] = [];
      for (const s of c.compression_trace) {
        const i = find(s.deg, false);
        errs.push(Math.abs(tr.pressure[i] / s.p - 1), Math.abs(tr.Tmean[i] / s.T - 1));
      }
      for (const s of c.expansion_trace_equilibrium) {
        const i = find(s.deg, true);
        errs.push(Math.abs(tr.pressure[i] / s.p - 1), Math.abs(tr.Tb[i] / s.T - 1));
      }
      errs.push(Math.abs(r.work / c.work - 1));
      errs.push(Math.abs(r.work / (r.fuelMass * c.lhv_gaseous) / c.efficiency - 1));
      const e = Math.max(...errs);
      if (e > worst) {
        worst = e;
        worstLabel = c.input.label;
      }
      console.log(`[validation] fuel–air ${c.input.label}: worst rel. error ${e.toExponential(2)} (W ${r.work.toFixed(3)} vs ${c.work.toFixed(3)} J)`);
      expect(e).toBeLessThan(tol);
      expect(Math.abs(r.U1 - (r.U0 - r.work)) / Math.abs(r.work)).toBeLessThan(1e-9);
    });
  }
  it('worst case', () => {
    console.log(`[validation] fuel–air extended oracle: worst ${worst.toExponential(2)} (${worstLabel})`);
    expect(worst).toBeLessThan(tol);
  });

  it('RK4 step-size order on the limit cycle (work, CR 18 boosted case)', () => {
    const c = fixture.cases.find((x) => x.input.cr === 18)!;
    const e = [2, 1, 0.5, 0.25].map((d) => Math.abs(run(c, d).work / c.work - 1));
    console.log(`[validation] fuel–air work error vs maxStepDeg 2/1/0.5/0.25°: ${e.map((x) => x.toExponential(2)).join(' / ')}`);
    expect(e[3]).toBeLessThan(1e-7);
  });
});
