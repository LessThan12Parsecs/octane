import { describe, expect, it } from 'vitest';
import { NS } from '../core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { completeCombustionProducts, freshCharge, fuelFromSelection, humidAir } from '../thermo/fuels';
import { mixMolarMass } from '../thermo/mixture';
import { CycleModel, MODE_OPEN, runClosedCycle, type CycleModelOptions, type CycleTrace } from './index';

/** RON-like trapped charge (PRF 90, φ 1.1, 6 % residual) for closed-cycle runs. */
function charge(): Float64Array {
  const op = CFR_RON_CONDITIONS;
  const fresh = freshCharge({ fuel: fuelFromSelection(op.fuel), phi: op.equivalenceRatio, airX: humidAir(325, 101325, 0.06) });
  const prod = completeCombustionProducts(fresh);
  const X = new Float64Array(NS);
  for (let k = 0; k < NS; k++) X[k] = (0.94 / mixMolarMass(fresh)) * fresh[k] + (0.06 / mixMolarMass(prod)) * prod[k];
  return X;
}

function closed(opts: Partial<CycleModelOptions>) {
  return runClosedCycle(
    CFR_F1,
    CFR_RON_CONDITIONS,
    { startDeg: -152, endDeg: 141, T: 370, p: 1.0e5, X: charge(), residualMassFraction: 0.06, uPrime: 5 },
    { knock: false, ...opts },
  );
}

/** Angle at which the traced x_b first reaches `lev` (linear interpolation). */
function caAt(tr: CycleTrace, lev: number): number {
  for (let i = 1; i < tr.xb.length; i++) {
    if (tr.mode[i] === MODE_OPEN) break;
    if (tr.xb[i] >= lev && tr.xb[i - 1] < lev) return tr.theta[i - 1] + ((lev - tr.xb[i - 1]) / (tr.xb[i] - tr.xb[i - 1])) * (tr.theta[i] - tr.theta[i - 1]);
  }
  return NaN;
}

describe('validation hooks: combustion models and calibration multipliers', () => {
  it("'wiebe' reproduces the prescribed mass-fraction-burned curve", () => {
    const w = { startDeg: -10, durationDeg: 50, a: -Math.log(0.001), m: 2 };
    const r = closed({ combustionModel: 'wiebe', wiebe: w });
    for (const lev of [0.1, 0.5, 0.9]) {
      const expected = w.startDeg + w.durationDeg * Math.pow(-Math.log(1 - lev) / w.a, 1 / (w.m + 1));
      expect(caAt(r.trace, lev)).toBeCloseTo(expected, 1);
    }
  });

  it('burnRateMultiplier (u_T = C·u′) advances the burn; taylorScaleMultiplier retards it', () => {
    // (explicit neutral baseline: the defaults are the calibrated CFR set since validation round 1)
    const neutral = { heatTransfer: false, burnRateMultiplier: 1, taylorScaleMultiplier: 1 };
    const base = caAt(closed(neutral).trace, 0.5);
    const fast = caAt(closed({ ...neutral, burnRateMultiplier: 2.6 }).trace, 0.5);
    const slowBurnUp = caAt(closed({ ...neutral, taylorScaleMultiplier: 2 }).trace, 0.5);
    console.log(`[cycle] closed-cycle CA50: default ${base.toFixed(1)}°, burnRateMultiplier 2.6 → ${fast.toFixed(1)}°, taylorScaleMultiplier 2 → ${slowBurnUp.toFixed(1)}°`);
    expect(fast).toBeLessThan(base - 5);
    expect(slowBurnUp).toBeGreaterThan(base);
  });

  it('woschniMultiplier scales the wall heat loss; heatTransfer false removes it', () => {
    // (explicit ×1 baseline: the default is the calibrated CFR value since validation round 1)
    const a = closed({ combustionModel: 'wiebe', woschniMultiplier: 1 });
    const b = closed({ combustionModel: 'wiebe', woschniMultiplier: 2 });
    const c = closed({ combustionModel: 'wiebe', heatTransfer: false });
    expect(b.heatLoss / a.heatLoss).toBeGreaterThan(1.6);
    expect(b.heatLoss / a.heatLoss).toBeLessThan(2.05);
    expect(c.heatLoss).toBe(0);
    expect(c.work).toBeGreaterThan(a.work);
  });

  it('ignitionDelayModel selects the knock model (Douaud–Eyzat vs LLNL tables)', () => {
    const op = { ...CFR_RON_CONDITIONS, compressionRatio: 9 };
    const run = (id: 'douaud-eyzat' | 'llnl-gasoline-2011') =>
      runClosedCycle(CFR_F1, op, { startDeg: -152, endDeg: 141, T: 370, p: 1e5, X: charge(), residualMassFraction: 0.06, uPrime: 5 }, { ignitionDelayModel: id }).model;
    const de = run('douaud-eyzat');
    const ll = run('llnl-gasoline-2011');
    expect(de.knockOnset || ll.knockOnset).toBe(true);
    expect(de.knockOnsetDeg).not.toBe(ll.knockOnsetDeg);
  });
});

describe('gas exchange integration', () => {
  it('valve flows do not chatter near Δp ≈ 0 (at most two physical reversals per valve event)', () => {
    for (const rpm of [600, 1800]) {
      const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, rpm }, { warmupCycles: 1 });
      m.recordTrace();
      m.runCycles(1);
      const tr = m.trace!;
      for (const f of [tr.mdotIntake, tr.mdotExhaust]) {
        let flips = 0;
        let alternating = 0;
        for (let i = 2; i < f.length; i++) {
          if (f[i] * f[i - 1] < 0) flips++;
          if (f[i] * f[i - 1] < 0 && f[i - 1] * f[i - 2] < 0) alternating++;
        }
        expect(flips).toBeLessThanOrEqual(2);
        expect(alternating).toBe(0);
      }
    }
  });
});
