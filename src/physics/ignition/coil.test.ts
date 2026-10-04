import { describe, expect, it } from 'vitest';
import type { InductiveIgnitionSpec } from '../core/engine-spec';
import fx from '../../../test/fixtures/ignition_circuit.json';
import { IgnitionCoil } from './coil';
import { SparkGap } from './discharge';

const P = fx.params;

/** Circuit spec identical to the oracle's parameters. */
function specFromFixture(over: Partial<InductiveIgnitionSpec> = {}): InductiveIgnitionSpec {
  return {
    type: 'inductive',
    supplyVoltage: P.Vs,
    primaryInductance: P.L1,
    primaryResistance: P.R1,
    secondaryInductance: P.L2,
    secondaryResistance: P.R2,
    secondaryCapacitance: P.C2,
    couplingCoefficient: P.k,
    primaryCurrentLimit: P.Ilim,
    dwellTime: P.t_dwell,
    ...over,
  };
}

const coilOpts = { primaryCapacitance: P.C1, primaryClampVoltage: P.Vclamp, primaryReverseClampVoltage: P.Vrev };

/** Gap with the oracle's glow model (fixed breakdown voltage, glow only). */
function oracleGap(suppress = false): SparkGap {
  const g = new SparkGap(P.l, { glowCathodeFall: P.Vsheath, glowAnodeFall: 0, extinctionCurrent: P.Iext });
  g.pressure = P.p;
  g.breakdownVoltage = P.Vbd;
  g.suppressBreakdown = suppress;
  return g;
}

function charge(coil: IgnitionCoil, gap: SparkGap, t: number, dt = 20e-6): void {
  const n = Math.round(t / dt);
  for (let i = 0; i < n; i++) coil.step(dt, true, gap);
}

describe('IgnitionCoil — dwell (primary RL charging)', () => {
  it('matches the Radau oracle, and the analytic RL law (Heywood eq. 9.53) once the make-ring has decayed', () => {
    const coil = new IgnitionCoil(specFromFixture(), coilOpts);
    const gap = oracleGap(true);
    let maxErr = 0;
    for (let i = 1; i < fx.dwell.t.length; i++) {
      coil.step(fx.dwell.t[i] - fx.dwell.t[i - 1], true, gap);
      // trapezoidal integrator vs the Radau oracle of the same ODEs (incl. the "make" ring
      // of the secondary LC excited by the switch-on step)
      const err = Math.abs(coil.I1 - fx.dwell.I1[i]);
      maxErr = Math.max(maxErr, err);
      // second-order trapezoidal error of the (under-resolved) make-ring decays with it
      if (fx.dwell.t[i] >= 1.5e-3) expect(err).toBeLessThan(1e-5);
      // after ≈ 1 ms the ring has decayed and I₁ follows (V_s/R₁)(1 − e^{−R₁t/L₁})
      if (fx.dwell.t[i] >= 1e-3) expect(Math.abs(coil.I1 / fx.dwell.I1_analytic[i] - 1)).toBeLessThan(1.5e-3);
    }
    expect(maxErr).toBeLessThan(1e-2);
    expect(coil.I1 / fx.dwell.end[0]).toBeCloseTo(1, 6);
    if (process.env.IGNITION_REPORT) console.log('dwell: max |I1 - I1_oracle| =', maxErr, 'A');
  });

  it('holds the driver current limit and books the driver dissipation (energy closes)', () => {
    const coil = new IgnitionCoil(specFromFixture({ primaryCurrentLimit: 5 }), coilOpts);
    const gap = oracleGap(true);
    charge(coil, gap, 12e-3);
    expect(coil.I1).toBeCloseTo(5, 6);
    expect(coil.currentLimited).toBe(true);
    // in limited operation the driver drops V_s − R₁ I_lim on average (the lightly damped
    // make-ring exchanges ±M·I_lim·ΔI₂ ≈ ±3 mJ with the supply, hence a long window)
    const before = coil.energyDriver;
    charge(coil, gap, 20e-3);
    expect((coil.energyDriver - before) / 20e-3 / ((P.Vs - P.R1 * 5) * 5)).toBeCloseTo(1, 1);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });
});

describe('IgnitionCoil — free secondary LC ring', () => {
  it('rings at f = 1/(2π√(L₂C₂)) with k = 0, R₂ = 0 and conserves ½C₂V² exactly', () => {
    const coil = new IgnitionCoil(specFromFixture({ couplingCoefficient: 0, secondaryResistance: 0 }), {
      ...coilOpts,
      tailSubstep: 1e-6,
    });
    const gap = oracleGap(true);
    coil.V2 = fx.lc.V0;
    const W0 = coil.storedEnergy();
    const dt = 1e-6;
    const crossings: number[] = [];
    let prev = coil.V2;
    let t = 0;
    for (let i = 0; i < 2000; i++) {
      coil.step(dt, false, gap);
      t += dt;
      if (prev > 0 && coil.V2 <= 0) crossings.push(t - (dt * coil.V2) / (coil.V2 - prev));
      prev = coil.V2;
    }
    const periods = (crossings[crossings.length - 1] - crossings[0]) / (crossings.length - 1);
    expect(1 / periods / fx.lc.f0).toBeCloseTo(1, 4);
    expect(Math.abs(coil.storedEnergy() / W0 - 1)).toBeLessThan(1e-9);
  });
});

describe('IgnitionCoil — open-circuit secondary (breakdown suppressed)', () => {
  it('rings up to the oracle peak voltage (≈ 20–45 kV) with the collector clamp active', () => {
    const coil = new IgnitionCoil(specFromFixture(), coilOpts);
    const gap = oracleGap(true);
    charge(coil, gap, P.t_dwell);
    coil.resetLedger();
    const dt = 1e-6;
    const samples: number[] = [];
    let t = 0;
    let si = 0;
    for (let i = 0; i < 200; i++) {
      coil.step(dt, false, gap);
      t += dt;
      while (si < fx.open.V2_t.length && t >= fx.open.V2_t[si] - 1e-12) {
        samples.push(coil.V2);
        si++;
      }
    }
    expect(coil.peakSecondaryVoltage / fx.open.peak_V2).toBeCloseTo(1, 2);
    expect(Math.abs(coil.peakSecondaryVoltageTime - P.t_dwell - fx.open.t_peak)).toBeLessThan(1e-6);
    for (let k = 0; k < samples.length; k++) {
      expect(Math.abs(samples[k] - fx.open.V2_s[k])).toBeLessThan(0.01 * Math.abs(fx.open.peak_V2));
    }
    expect(coil.peakSecondaryVoltage).toBeGreaterThan(20e3);
    expect(coil.peakSecondaryVoltage).toBeLessThan(45e3);
    expect(coil.energyClamp / fx.open.E_clamp).toBeCloseTo(1, 1);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });
});

describe('IgnitionCoil — breakdown and glow discharge', () => {
  function fire(dt: number) {
    const coil = new IgnitionCoil(specFromFixture(), coilOpts);
    const gap = oracleGap(false);
    charge(coil, gap, P.t_dwell);
    const Wsw = coil.storedEnergy();
    const halfLI2 = 0.5 * P.L1 * coil.I1 * coil.I1;
    coil.resetLedger();
    const tOff = coil.time;
    const i2: number[] = [];
    const tEnd = tOff + 8e-3;
    let si = 0;
    while (coil.time < tEnd - 1e-12) {
      // step exactly onto the oracle's sample times once they are known (after breakdown)
      let h = Math.min(dt, tEnd - coil.time);
      if (!Number.isNaN(gap.firstBreakdownTime) && si < fx.discharge.I2_t.length) {
        const ts = tOff + fx.discharge.I2_t[si];
        if (coil.time + h >= ts - 1e-12) h = ts - coil.time;
      }
      coil.step(h, false, gap);
      if (si < fx.discharge.I2_t.length && Math.abs(coil.time - tOff - fx.discharge.I2_t[si]) < 1e-12) {
        i2.push(Math.abs(coil.I2));
        si++;
      }
    }
    return { coil, gap, Wsw, halfLI2, tOff, i2 };
  }

  it('matches the Radau oracle: breakdown delay, duration, gap and R₂ energies, I₂(t)', () => {
    const { coil, gap, tOff, i2 } = fire(20e-6);
    const d = fx.discharge;
    expect(Math.abs(gap.firstBreakdownTime - tOff - d.t_bd)).toBeLessThan(0.02 * d.t_bd);
    const dur = gap.extinctionTime - gap.firstBreakdownTime;
    expect(dur / d.duration).toBeCloseTo(1, 3);
    expect(gap.energyTotal / d.E_gap).toBeCloseTo(1, 3);
    expect(gap.energyBreakdown + gap.energyCapacitiveArc).toBeCloseTo(d.E_dump, 5);
    expect(i2.length).toBe(d.I2_s.length);
    // 0.1 ms after breakdown the (under-resolved, physical) 400 kHz primary leakage ring still
    // modulates I₂ by ≈ 1 %; once it has decayed the waveforms agree to < 0.5 %
    for (let k = 0; k < i2.length; k++) expect(Math.abs(i2[k] / d.I2_s[k] - 1)).toBeLessThan(k === 0 ? 0.02 : 0.005);
    expect(coil.energySecondaryResistance).toBeGreaterThan(d.E_R2); // + post-extinction ring
    expect(coil.energySecondaryResistance / d.E_R2).toBeLessThan(1.01);
  });

  it('closes the energy balance: ½L₁I₁²(switch-off) = gap + R₁ + R₂ + clamp + stored (≤ 1e-9)', () => {
    const { coil, Wsw, halfLI2 } = fire(20e-6);
    const out =
      coil.energyGap +
      coil.energyPrimaryResistance +
      coil.energySecondaryResistance +
      coil.energyClamp +
      coil.energyDriver +
      coil.storedEnergy() -
      coil.energySupplied;
    expect(Math.abs(out / Wsw - 1)).toBeLessThan(1e-9);
    // the primary magnetic energy is > 99 % of the stored energy at switch-off
    expect(Math.abs(halfLI2 / Wsw - 1)).toBeLessThan(0.01);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('is insensitive to the caller step (10 µs vs 50 µs)', () => {
    const a = fire(10e-6);
    const b = fire(50e-6);
    expect(a.gap.energyTotal / b.gap.energyTotal).toBeCloseTo(1, 4);
    expect((a.gap.extinctionTime - a.gap.firstBreakdownTime) / (b.gap.extinctionTime - b.gap.firstBreakdownTime)).toBeCloseTo(1, 4);
  });
});

describe('IgnitionCoil — post-spark ringing tail (reviewer change)', () => {
  it('snaps the ringing to rest only once it can no longer reach 25 V on the plug; ledger exact', () => {
    const coil = new IgnitionCoil(specFromFixture(), coilOpts);
    const gap = oracleGap(false);
    charge(coil, gap, P.t_dwell);
    coil.resetLedger();
    const dt = 20e-6;
    let tRest = NaN;
    let maxV2AfterExt = 0;
    for (let t = 0; t < 300e-3; t += dt) {
      coil.step(dt, false, gap);
      if (!Number.isNaN(gap.extinctionTime)) maxV2AfterExt = Math.max(maxV2AfterExt, Math.abs(coil.V2));
      if (Number.isNaN(tRest) && coil.deviationEnergy() === 0) tRest = coil.time;
    }
    expect(gap.breakdownCount).toBe(1);
    // at rest well within the cycle (the R₂-damped L₂C₂ ring alone would need ≈ 0.2 s)
    expect(tRest - gap.extinctionTime).toBeLessThan(0.1);
    expect(coil.I1).toBe(0);
    expect(coil.V2).toBe(0);
    expect(coil.V1).toBe(P.Vs);
    // what was dropped could never have reached a breakdown (≥ 267 V): ≤ ½C₂(25 V)²
    expect(Math.abs(coil.energyTailTruncated)).toBeLessThanOrEqual(0.5 * P.C2 * 25 * 25);
    expect(maxV2AfterExt).toBeLessThan(P.Vbd);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });
});
