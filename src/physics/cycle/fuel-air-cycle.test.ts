import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/cycle_fuel_air.json';
import type { EngineSpec } from '../core/engine-spec';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { fuelFromSelection, lowerHeatingValue } from '../thermo/fuels';
import { MODE_BURNED, runClosedCycle } from './index';

/**
 * Fuel–air-cycle oracle (tools/reference/cycle_fuel_air_oracle.py, Cantera): frozen isentropic
 * compression from BDC, adiabatic constant-volume combustion to equilibrium at TDC, isentropic
 * expansion in shifting equilibrium. The cycle model runs the same limit in closed-cycle mode
 * with heat transfer off and 'instantaneous-at-tdc' combustion; p and T at every oracle point,
 * the work and the efficiency must agree within 0.1 % (they agree to ~1e-8).
 */
describe('closed cycle vs the Cantera fuel–air cycle', () => {
  const tol = 1e-3;
  let worst = 0;
  for (const c of fixture.cases) {
    const inp = c.input;
    const label = `${inp.fuel} CR ${inp.cr} φ ${inp.phi} x_r ${inp.residual}`;
    it(label, () => {
      const spec: EngineSpec = {
        ...CFR_F1,
        geometry: { ...CFR_F1.geometry, bore: inp.bore, stroke: inp.stroke, conRodLength: inp.rod, pinOffset: 0, compressionRatioRange: [2, 20] },
      };
      const fuel = inp.fuel === 'CH4' ? ({ kind: 'pure', species: 'CH4' } as const) : ({ kind: 'pure', species: 'IC8H18' } as const);
      const op = { ...CFR_RON_CONDITIONS, compressionRatio: inp.cr, fuel, equivalenceRatio: inp.phi };
      const r = runClosedCycle(
        spec,
        op,
        { startDeg: -180, endDeg: 180, T: inp.T1, p: inp.p1, X: c.X_reactants as unknown as Readonly<Record<string, number>> },
        { heatTransfer: false, combustionModel: 'instantaneous-at-tdc', knock: false },
      );
      const tr = r.trace;
      expect(r.mass / c.mass - 1).toBeCloseTo(0, 9);
      const find = (deg: number, burned: boolean): number => {
        for (let i = 0; i < tr.theta.length; i++) {
          if (Math.abs(tr.theta[i] - deg) < 1e-7 && (tr.mode[i] === MODE_BURNED) === burned) return i;
        }
        throw new Error(`no trace sample at ${deg}° (burned ${burned})`);
      };
      for (const s of c.compression_trace) {
        const i = find(s.deg, false);
        const ep = Math.abs(tr.pressure[i] / s.p - 1);
        const eT = Math.abs(tr.Tmean[i] / s.T - 1);
        worst = Math.max(worst, ep, eT);
        expect(ep).toBeLessThan(tol);
        expect(eT).toBeLessThan(tol);
      }
      for (const s of c.expansion_trace_equilibrium) {
        const i = find(s.deg, true);
        const ep = Math.abs(tr.pressure[i] / s.p - 1);
        const eT = Math.abs(tr.Tb[i] / s.T - 1);
        worst = Math.max(worst, ep, eT);
        expect(ep).toBeLessThan(tol);
        expect(eT).toBeLessThan(tol);
      }
      // work and efficiency (fuel LHV of the gaseous fuel)
      const ew = Math.abs(r.work / c.work_equilibrium_expansion - 1);
      const fb = fuelFromSelection(fuel);
      const lhv = lowerHeatingValue(fb);
      expect(lhv / c.lhv_gaseous - 1).toBeCloseTo(0, 6);
      const eta = r.work / (r.fuelMass * lhv);
      const eeta = Math.abs(eta / c.efficiency_equilibrium - 1);
      worst = Math.max(worst, ew, eeta);
      expect(ew).toBeLessThan(tol);
      expect(eeta).toBeLessThan(tol);
      // adiabatic, no spark: U_end = U_0 − W exactly
      expect(Math.abs(r.U1 - (r.U0 - r.work)) / Math.abs(r.work)).toBeLessThan(1e-9);
    });
  }
  it('worst relative deviation', () => {
    console.log(`[cycle] fuel–air oracle: worst relative p/T/work/efficiency error ${worst.toExponential(2)}`);
    expect(worst).toBeLessThan(tol);
  });
});
