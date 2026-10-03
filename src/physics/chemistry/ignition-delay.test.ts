import { describe, expect, it } from 'vitest';
import { P_ATM } from '../core/constants';
import { fuelFromSelection, prfIsooctaneMoleFraction } from '../thermo/fuels';
import {
  douaudEyzat,
  douaudEyzatTau,
  type IgnitionDelayTable,
  ignitionDelayTemperatureSensitivity,
  prfOctaneNumber,
  residualMoleFraction,
  TabulatedIgnitionDelay,
} from './ignition-delay';
import { NS, SP } from '../core/species';

describe('Douaud–Eyzat (1978)', () => {
  it('τ = 17.68 ms (ON/100)^3.402 (p/atm)^−1.7 exp(3800/T)', () => {
    // hand value: ON 90, 40 atm, 900 K
    const ref = 17.68e-3 * Math.pow(0.9, 3.402) * Math.pow(40, -1.7) * Math.exp(3800 / 900);
    expect(douaudEyzatTau(900, 40 * P_ATM, 90)).toBeCloseTo(ref, 15);
    expect(ref).toBeCloseTo(1.592e-3, 5); // 17.68e-3 · 0.6990 · 1.892e-3 · 68.18 ≈ 1.59 ms
    const fuel = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    expect(douaudEyzat.tau(900, 40 * P_ATM, 1, fuel, 0.1)).toBe(douaudEyzatTau(900, 40 * P_ATM, 90));
  });

  it('pressure and temperature sensitivities are the published exponents', () => {
    const fuel = fuelFromSelection({ kind: 'PRF', octaneNumber: 95 });
    const T = 850;
    const s = ignitionDelayTemperatureSensitivity(douaudEyzat, T, 3e6, 1, fuel, 0, 0.5);
    expect(s).toBeCloseTo(-3800 / (T * T), 8);
    const r = douaudEyzatTau(T, 6e6, 95) / douaudEyzatTau(T, 3e6, 95);
    expect(Math.log(r) / Math.log(2)).toBeCloseTo(-1.7, 12);
  });

  it('rejects non-PRF fuels', () => {
    const methane = fuelFromSelection({ kind: 'pure', species: 'CH4' });
    expect(() => douaudEyzat.tau(900, 4e6, 1, methane, 0)).toThrow(RangeError);
  });
});

describe('prfOctaneNumber / residualMoleFraction', () => {
  it('uses octaneNumber when present and inverts the liquid-volume blending otherwise', () => {
    expect(prfOctaneNumber(fuelFromSelection({ kind: 'PRF', octaneNumber: 87 }))).toBe(87);
    for (const on of [0, 12.5, 50, 90, 100]) {
      const X = new Float64Array(NS);
      const x = prfIsooctaneMoleFraction(on);
      X[SP.IC8H18] = x;
      X[SP.NC7H16] = 1 - x;
      expect(prfOctaneNumber({ label: 'x', X })).toBeCloseTo(on, 10);
    }
    const X = new Float64Array(NS);
    X[SP.C2H5OH] = 0.1;
    X[SP.IC8H18] = 0.9;
    expect(prfOctaneNumber({ label: 'E10-ish', X })).toBeNaN();
  });

  it('mass → mole fraction of residual', () => {
    expect(residualMoleFraction(0.1, 0.0288, 0.0288)).toBeCloseTo(0.1, 14);
    // heavier fresh charge (fuel vapour) → larger residual mole fraction
    const x = residualMoleFraction(0.1, 0.0285, 0.0303);
    expect(x).toBeCloseTo(0.1 / 0.0285 / (0.1 / 0.0285 + 0.9 / 0.0303), 14);
  });
});

/** Synthetic table of a function that is linear in each variable separately (multilinear). */
function syntheticTable(f: (invT: number, lnP: number, phi: number, on: number, x: number) => number): IgnitionDelayTable {
  const invT = { start: 1000 / 1100, step: (1000 / 650 - 1000 / 1100) / 6, n: 7 };
  const lnP = { start: Math.log(1e6), step: Math.log(8) / 3, n: 4 };
  const phi = [0.5, 1, 1.5];
  const on = [0, 50, 80, 100];
  const xres = [0, 0.1, 0.2];
  const lnTau: number[] = [];
  for (let i = 0; i < invT.n; i++)
    for (let j = 0; j < lnP.n; j++)
      for (const ph of phi)
        for (const o of on)
          for (const x of xres) lnTau.push(f(invT.start + i * invT.step, lnP.start + j * lnP.step, ph, o, x));
  return {
    description: 'synthetic', reactor: 'test', tauMax: 2,
    axes: { invT, lnP, phi, on, xres },
    order: ['invT', 'lnP', 'phi', 'on', 'xres'],
    lnTau, lnTau1: lnTau.map((v) => v - 0.5),
  };
}

describe('TabulatedIgnitionDelay (interpolation mechanics)', () => {
  const f = (it: number, lp: number, ph: number, on: number, x: number): number =>
    -12 + 9 * it - 1.3 * (lp - 15) - 0.4 * ph + 0.01 * on * (1 + 0.5 * it) + 2 * x + 0.3 * it * (lp - 15) * ph;
  const model = new TabulatedIgnitionDelay('synthetic', syntheticTable(f));

  it('is exact for multilinear functions inside the table', () => {
    let seed = 1;
    const rnd = (): number => ((seed = (seed * 16807) % 2147483647) / 2147483647);
    for (let k = 0; k < 500; k++) {
      const it = 1000 / 1100 + rnd() * (1000 / 650 - 1000 / 1100);
      const lp = Math.log(1e6) + rnd() * Math.log(8);
      const ph = 0.5 + rnd();
      const on = 100 * rnd();
      const x = 0.2 * rnd();
      const v = model.lnTauAt(1000 / it, Math.exp(lp), ph, on, x);
      expect(v).toBeCloseTo(f(it, lp, ph, on, x), 10);
      expect(model.lnTauFirstStageAt(1000 / it, Math.exp(lp), ph, on, x)).toBeCloseTo(f(it, lp, ph, on, x) - 0.5, 10);
    }
  });

  it('clamps φ, ON, x_res; extrapolates T and p linearly with limited slopes', () => {
    const T = 900;
    const p = 3e6;
    expect(model.lnTauAt(T, p, 3, 50, 0.1)).toBeCloseTo(model.lnTauAt(T, p, 1.5, 50, 0.1), 12);
    expect(model.lnTauAt(T, p, 1, 150, 0.1)).toBeCloseTo(model.lnTauAt(T, p, 1, 100, 0.1), 12);
    expect(model.lnTauAt(T, p, 1, 50, -1)).toBeCloseTo(model.lnTauAt(T, p, 1, 50, 0), 12);
    // x_res above the table (last node 0.2): linear extrapolation up to 0.4, clamped beyond
    expect(model.lnTauAt(T, p, 1, 50, 0.35)).toBeCloseTo(f(1000 / T, Math.log(p), 1, 50, 0.35), 9);
    expect(model.lnTauAt(T, p, 1, 50, 0.9)).toBeCloseTo(f(1000 / T, Math.log(p), 1, 50, 0.4), 9);
    const mNeg = new TabulatedIgnitionDelay('neg', syntheticTable((a, b, c, d, x) => f(a, b, c, d, 0) - 3 * x));
    expect(mNeg.lnTauAt(T, p, 1, 50, 0.35)).toBeCloseTo(mNeg.lnTauAt(T, p, 1, 50, 0.2) - 3 * 0.15, 9);
    // f has ∂/∂invT > 0 and ∂/∂lnP < 0 everywhere: plain linear extrapolation
    const it = 1000 / 550;
    expect(model.lnTauAt(550, p, 1, 50, 0.1)).toBeCloseTo(f(it, Math.log(p), 1, 50, 0.1), 9);
    expect(model.lnTauAt(900, 2e7, 1, 50, 0.1)).toBeCloseTo(f(1000 / 900, Math.log(2e7), 1, 50, 0.1), 9);
    // wrong-signed edge slopes are limited to zero
    const g = (it2: number, lp: number): number => -5 - 2 * (it2 - 1.2) * (it2 - 1.2) + 0.2 * (lp - 15) ** 2;
    const m2 = new TabulatedIgnitionDelay('g', syntheticTable((a, b) => g(a, b)));
    const edgeLow = m2.lnTauAt(650, 1e6, 1, 50, 0);
    expect(m2.lnTauAt(500, 1e6, 1, 50, 0)).toBeCloseTo(edgeLow, 12); // g decreasing in invT at the edge → flat
    const edgeP = m2.lnTauAt(1000, 8e6, 1, 50, 0);
    expect(m2.lnTauAt(1000, 3e7, 1, 50, 0)).toBeCloseTo(edgeP, 12); // increasing with p at the edge → flat
  });

  it('tau(fuel) goes through prfOctaneNumber; non-PRF → NaN', () => {
    const fuel = fuelFromSelection({ kind: 'PRF', octaneNumber: 80 });
    expect(model.tau(800, 4e6, 1, fuel, 0.05)).toBeCloseTo(Math.exp(model.lnTauAt(800, 4e6, 1, 80, 0.05)), 14);
    expect(model.tau(800, 4e6, 1, fuelFromSelection({ kind: 'pure', species: 'C3H8' }), 0)).toBeNaN();
  });

  it('non-uniform (node-list) 1000/T and ln p axes: exact for multilinear f, same extrapolation', () => {
    const u = syntheticTable(f);
    const inv = [1000 / 1100, 1.0, 1.2, 1000 / 650, 1000 / 600, 1000 / 550];
    const lnp = [Math.log(3e5), Math.log(1e6), Math.log(2.5e6), Math.log(8e6)];
    const t: IgnitionDelayTable = {
      ...u,
      axes: { ...u.axes, invT: { nodes: inv }, lnP: { nodes: lnp } },
      lnTau: [],
      lnTau1: [],
    };
    for (const it of inv)
      for (const lp of lnp)
        for (const ph of u.axes.phi)
          for (const o of u.axes.on) for (const x of u.axes.xres) t.lnTau.push(f(it, lp, ph, o, x));
    t.lnTau1 = t.lnTau.slice();
    const m = new TabulatedIgnitionDelay('nodes', t);
    let seed = 3;
    const rnd = (): number => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (let k = 0; k < 300; k++) {
      // inside, and up to 30 % outside both ends of each extrapolated axis
      const it = inv[0] + (rnd() * 1.6 - 0.3) * (inv[inv.length - 1] - inv[0]);
      const lp = lnp[0] + (rnd() * 1.6 - 0.3) * (lnp[lnp.length - 1] - lnp[0]);
      const ph = 0.5 + rnd();
      const on = 100 * rnd();
      const x = 0.2 * rnd();
      // f is increasing in 1000/T and decreasing in ln p over this box ⇒ no limiter active
      expect(m.lnTauAt(1000 / it, Math.exp(lp), ph, on, x)).toBeCloseTo(f(it, lp, ph, on, x), 9);
    }
    const bad = { ...t, axes: { ...t.axes, invT: { nodes: [1, 0.9] } } };
    expect(() => new TabulatedIgnitionDelay('bad', bad)).toThrow(RangeError);
  });

  it('rejects inconsistent tables', () => {
    const t = syntheticTable(f);
    t.lnTau.pop();
    expect(() => new TabulatedIgnitionDelay('bad', t)).toThrow(RangeError);
  });
});
