import { describe, expect, it } from 'vitest';
import { NS, SP } from '../core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { completeCombustionProducts, freshCharge, fuelFromSelection, humidAir } from '../thermo/fuels';
import { mixMolarMass, mixUMolar } from '../thermo/mixture';
import { ZoneClosure } from './closure';
import { MODE_BURNED, MODE_SINGLE, MODE_TWO, runClosedCycle, swapNO, type ClosedCycleResult } from './index';

/** RON-like trapped charge at IVC: PRF 90 φ 1.1 fresh charge + 6 % (mass) residual products. */
function ronCharge(): { X: Float64Array; yRes: number } {
  const op = CFR_RON_CONDITIONS;
  const fuel = fuelFromSelection(op.fuel);
  const air = humidAir(op.ambientTemperature, op.ambientPressure, op.relativeHumidity);
  const fresh = freshCharge({ fuel, phi: op.equivalenceRatio, airX: air });
  const prod = completeCombustionProducts(fresh);
  const yRes = 0.06;
  const Mf = mixMolarMass(fresh);
  const Mp = mixMolarMass(prod);
  const X = new Float64Array(NS);
  for (let k = 0; k < NS; k++) X[k] = ((1 - yRes) / Mf) * fresh[k] + (yRes / Mp) * prod[k];
  return { X, yRes };
}

function run(opts: { heatTransfer: boolean; knock: boolean; cr?: number }): ClosedCycleResult {
  const { X, yRes } = ronCharge();
  const op = { ...CFR_RON_CONDITIONS, compressionRatio: opts.cr ?? CFR_RON_CONDITIONS.compressionRatio };
  return runClosedCycle(
    CFR_F1,
    op,
    { startDeg: -152, endDeg: 141, T: 370, p: 1.0e5, X, residualMassFraction: yRes, uPrime: 5 },
    { heatTransfer: opts.heatTransfer, knock: opts.knock, combustionModel: 'entrainment' },
  );
}

/** Independent re-evaluation of the zone energies from the traced (p, T_u, T_b). */
function zoneEnergy(r: ClosedCycleResult, i: number, cl: ZoneClosure): number {
  const tr = r.trace;
  let U = 0;
  if (tr.mu[i] > 0) {
    cl.unburnedAtT(tr.Tu[i], tr.pressure[i]);
    U += tr.mu[i] * cl.uu;
  }
  if (tr.mb[i] > 0) {
    const b = new Float64Array(cl.bu.length);
    for (let e = 0; e < b.length; e++) b[e] = cl.bu[e] * tr.mb[i];
    const res = cl.eq.solveTP(b, tr.Tb[i], tr.pressure[i]);
    // (fixer round 2: the burned zone carries the rate-controlled NO in its thermodynamics — the
    // equilibrium composition with its NO swapped for the kinetic amount the closure used)
    const N = Float64Array.from(res.N);
    swapNO(N, tr.nNOClosure[i]);
    U += mixUMolar(N, tr.Tb[i]);
  }
  return U;
}

describe('closed-cycle energy conservation with spark, entrainment and knock (adiabatic)', () => {
  for (const knock of [false, true]) {
    it(`U + ∫p dV − E_spark is conserved (knock ${knock ? 'on, CR 9' : 'off'})`, () => {
      const r = run({ heatTransfer: false, knock, cr: knock ? 9 : undefined });
      const tr = r.trace;
      const n = tr.theta.length;
      const m = r.model;
      // it really burned (and knocked when asked); the last sample is the post-EVO merge
      expect(tr.xb[n - 2]).toBeGreaterThan(0.95);
      expect(r.sparkEnergy).toBeGreaterThan(1e-3);
      expect(r.sparkEnergy).toBeLessThan(0.1);
      if (knock) {
        expect(m.knockOnset).toBe(true);
        expect(m.mapo).toBeGreaterThan(0);
      }
      // (1) the ledger identity at every step
      const scale = r.fuelMass * 44e6; // chemical energy scale, J
      let worst = 0;
      for (let i = 0; i < n; i++) {
        const res = tr.U[i] + (tr.work[i] - tr.work[0]) - (tr.sparkEnergy[i] - tr.sparkEnergy[0]) - r.U0;
        worst = Math.max(worst, Math.abs(res) / Math.abs(r.U0));
      }
      expect(worst).toBeLessThan(1e-6);
      expect(Math.abs(r.U1 + r.work - r.sparkEnergy - r.U0) / scale).toBeLessThan(1e-9);
      // (2) closure consistency: zone energies recomputed from (p, T_u, T_b) equal U_tot
      const cl = new ZoneClosure();
      cl.setUnburned(m.closure.Xu);
      let worstZ = 0;
      for (let i = 1; i < n; i += 37) {
        if (tr.mode[i] !== MODE_SINGLE && tr.mode[i] !== MODE_TWO && tr.mode[i] !== MODE_BURNED) continue;
        worstZ = Math.max(worstZ, Math.abs(zoneEnergy(r, i, cl) - tr.U[i]) / scale);
      }
      expect(worstZ).toBeLessThan(1e-8);
      // (3) the RK work integral equals an independent trapezoidal ∫p dV over the step samples
      let wt = 0;
      for (let i = 1; i < n; i++) wt += 0.5 * (tr.pressure[i] + tr.pressure[i - 1]) * (tr.volume[i] - tr.volume[i - 1]);
      expect(Math.abs(wt - r.work) / Math.abs(r.work)).toBeLessThan(1e-4);
      console.log(
        `[cycle] adiabatic closed cycle (knock ${knock}): W = ${r.work.toFixed(2)} J, E_spark = ${(r.sparkEnergy * 1e3).toFixed(2)} mJ, ` +
          `ledger residual ${worst.toExponential(1)} (rel. |U0|), zone-energy residual ${worstZ.toExponential(1)} (rel. m_f LHV), ` +
          `trapezoid vs RK work ${((wt - r.work) / r.work).toExponential(1)}` +
          (knock ? `, knock onset ${m.knockOnsetDeg.toFixed(2)}°, end gas ${m.knockEndGasFraction.toFixed(3)}, MAPO ${(m.mapo / 1e5).toFixed(2)} bar` : ''),
      );
    });
  }

  it('with wall heat transfer: U + ∫p dV + Q_wall − E_spark is conserved', () => {
    const r = run({ heatTransfer: true, knock: true });
    const scale = r.fuelMass * 44e6;
    expect(r.heatLoss).toBeGreaterThan(0.05 * scale);
    expect(Math.abs(r.U1 + r.work + r.heatLoss - r.sparkEnergy - r.U0) / scale).toBeLessThan(1e-9);
    // burned-gas NO formed and frozen by EVO
    expect(r.model.closure.Xu[SP.NO]).toBe(0);
  });
});
