import { describe, expect, it } from 'vitest';
import { R_UNIVERSAL } from '../core/constants';
import { EL, NE, NS, SP } from '../core/species';
import mixFx from '../../../test/fixtures/thermo_mixtures.json';
import {
  elementMoles,
  lastInverse,
  mixCpMass,
  mixCpMolar,
  mixCvMass,
  mixCvMolar,
  mixDensity,
  mixGamma,
  mixGasConstant,
  mixHMass,
  mixHMolar,
  mixMolarMass,
  mixSMass,
  mixSMolar,
  mixSoundSpeed,
  mixUMass,
  mixUMolar,
  molesToX,
  temperatureFromH,
  temperatureFromHMolar,
  temperatureFromS,
  temperatureFromU,
  temperatureFromUMolar,
  totalMass,
  totalMoles,
  xToY,
  yToX,
} from './mixture';
import { MOLAR_MASS, P_REF_THERMO, speciesS0R } from './thermo';

const rel = (a: number, b: number, scale = Math.abs(b)): number => Math.abs(a - b) / Math.max(scale, 1e-300);

const cases = mixFx.cases.map((c) => ({ ...c, Xa: Float64Array.from(c.X) }));

describe('composition conversions', () => {
  const N = new Float64Array(NS);
  N[SP.N2] = 3.1;
  N[SP.O2] = 0.8;
  N[SP.IC8H18] = 0.05;
  N[SP.H2O] = 0.2;
  N[SP.C2H5OH] = 0.01;

  it('molesToX / xToY / yToX round-trip and normalise', () => {
    const X = molesToX(N);
    expect(X.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 15);
    const Y = xToY(X);
    expect(Y.reduce((a, b) => a + b, 0)).toBeCloseTo(1, 15);
    const X2 = yToX(Y);
    for (let k = 0; k < NS; k++) expect(X2[k]).toBeCloseTo(X[k], 15);
    // Y_k = N_k M_k / m
    const m = totalMass(N);
    for (let k = 0; k < NS; k++) expect(Y[k]).toBeCloseTo((N[k] * MOLAR_MASS[k]) / m, 15);
    expect(totalMoles(N)).toBeCloseTo(4.16, 14);
    // in-place (aliasing) conversions
    const Z = Float64Array.from(X);
    xToY(Z, Z);
    for (let k = 0; k < NS; k++) expect(Z[k]).toBeCloseTo(Y[k], 15);
  });

  it('elementMoles counts atoms', () => {
    const b = elementMoles(N);
    expect(b.length).toBe(NE);
    expect(b[EL.C]).toBeCloseTo(0.05 * 8 + 0.01 * 2, 14);
    expect(b[EL.H]).toBeCloseTo(0.05 * 18 + 0.2 * 2 + 0.01 * 6, 14);
    expect(b[EL.O]).toBeCloseTo(0.8 * 2 + 0.2 + 0.01, 14);
    expect(b[EL.N]).toBeCloseTo(6.2, 14);
    expect(b[EL.AR]).toBe(0);
  });
});

describe('mixture properties vs Cantera', () => {
  it('cp, cv, h, u, s, M, γ (mass and molar) to 1e-9 relative, 200–6000 K, 1 and 50 bar', () => {
    let worst = 0;
    for (const c of cases) {
      const X = c.Xa;
      for (const st of c.states) {
        const { T, p } = st;
        const cpScale = st.cpMass * T; // h/u scale (h crosses 0 near 298 K)
        const errs = [
          rel(mixMolarMass(X), st.molarMass),
          rel(mixCpMass(X, T), st.cpMass),
          rel(mixCvMass(X, T), st.cvMass),
          rel(mixHMass(X, T), st.hMass, Math.max(Math.abs(st.hMass), cpScale)),
          rel(mixUMass(X, T), st.uMass, Math.max(Math.abs(st.uMass), cpScale)),
          rel(mixSMass(X, T, p), st.sMass),
          rel(mixCpMolar(X, T), st.cpMole),
          rel(mixHMolar(X, T), st.hMole, Math.max(Math.abs(st.hMole), st.cpMole * T)),
          rel(mixSMolar(X, T, p), st.sMole),
          rel(mixGamma(X, T), st.gamma),
        ];
        const e = Math.max(...errs);
        worst = Math.max(worst, e);
        expect(e, `${c.name} T=${T} p=${p}: ${errs.map((x) => x.toExponential(1)).join(',')}`).toBeLessThan(1e-9);
      }
    }
    console.info(`[mixture] max rel. error vs Cantera over ${cases.length} mixtures: ${worst.toExponential(2)}`);
  });

  it('derived quantities are consistent', () => {
    const X = cases[1].Xa;
    for (const T of [300, 900, 2500]) {
      expect(mixCvMolar(X, T)).toBeCloseTo(mixCpMolar(X, T) - R_UNIVERSAL, 10);
      expect(rel(mixGasConstant(X), R_UNIVERSAL / mixMolarMass(X))).toBeLessThan(1e-15);
      expect(rel(mixUMolar(X, T), mixHMolar(X, T) - R_UNIVERSAL * T, mixCpMolar(X, T) * T)).toBeLessThan(1e-14);
      const a = mixSoundSpeed(X, T);
      expect(rel(a * a, (mixGamma(X, T) * R_UNIVERSAL * T) / mixMolarMass(X))).toBeLessThan(1e-14);
      expect(rel(mixDensity(X, T, 1e5), (1e5 * mixMolarMass(X)) / (R_UNIVERSAL * T))).toBeLessThan(1e-15);
    }
  });

  it('molar functions are linear in the composition vector (N → extensive)', () => {
    const X = cases[4].Xa; // equilibrium products
    const n = 0.0123;
    const N = X.map((x) => x * n);
    for (const T of [400, 2200, 7000]) {
      expect(rel(mixHMolar(N, T), n * mixHMolar(X, T), n * mixCpMolar(X, T) * T)).toBeLessThan(1e-14);
      expect(rel(mixUMolar(N, T), n * mixUMolar(X, T), n * mixCpMolar(X, T) * T)).toBeLessThan(1e-14);
      expect(rel(mixCpMolar(N, T), n * mixCpMolar(X, T))).toBeLessThan(1e-14);
      expect(rel(mixCvMolar(N, T), n * mixCvMolar(X, T))).toBeLessThan(1e-14);
    }
  });

  it('ideal-mixing entropy: X = 0 species contribute nothing; pressure term −R ln(p/p_ref)', () => {
    const X = cases[0].Xa; // dry air: most species zero
    const T = 700;
    let s = 0;
    for (let k = 0; k < NS; k++) if (X[k] > 0) s += X[k] * (speciesS0R(k, T) - Math.log(X[k]));
    expect(rel(mixSMolar(X, T, P_REF_THERMO), R_UNIVERSAL * s)).toBeLessThan(1e-14);
    const ds = mixSMolar(X, T, 10 * P_REF_THERMO) - mixSMolar(X, T, P_REF_THERMO);
    expect(ds).toBeCloseTo(-R_UNIVERSAL * Math.log(10), 10);
  });

  it('ds/dT = cp/T even with a slightly negative (integrator round-off) mole fraction', () => {
    const X = Float64Array.from(cases[3].Xa); // complete products
    X[SP.CO] = -1e-4; // a species that is otherwise absent
    for (const T of [700, 2500]) {
      const dT = 1e-3 * T;
      const dsdT = (mixSMolar(X, T + dT, 2e6) - mixSMolar(X, T - dT, 2e6)) / (2 * dT);
      expect(rel(dsdT, mixCpMolar(X, T) / T)).toBeLessThan(1e-6);
    }
  });
});

describe('inverse functions', () => {
  const Ts = [120, 200, 250, 298.15, 555.5, 999.9999, 1000, 1000.0001, 1600, 2700, 4400, 6000, 6100, 9000];
  const guesses = [300, 3000, 50000, -5, NaN];

  it('temperatureFromH/U/S round-trip to 1e-10 from poor guesses (incl. extrapolated range)', () => {
    let worst = 0;
    let worstTmid = 0;
    let maxIt = 0;
    for (const c of cases) {
      const X = c.Xa;
      for (const T of Ts) {
        const h = mixHMass(X, T);
        const u = mixUMass(X, T);
        const p = 3.3e6;
        const s = mixSMass(X, T, p);
        for (const g of guesses) {
          for (const [Tx, name] of [
            [temperatureFromH(X, h, g), 'h'],
            [temperatureFromU(X, u, g), 'u'],
            [temperatureFromS(X, s, p, g), 's'],
          ] as const) {
            expect(lastInverse.converged).toBe(true);
            const e = rel(Tx, T);
            if (T === 1000) worstTmid = Math.max(worstTmid, e);
            else worst = Math.max(worst, e);
            maxIt = Math.max(maxIt, lastInverse.iterations);
            // At exactly Tmid the NASA-7 jump in h/s (|Δh| ≤ 1e-3 J/mol per species, see
            // thermo.test.ts) moves the root by |Δh|/c_p ≤ 1e-3/29 K ≈ 3.4e-5 K, i.e. ≤ 3.4e-8
            // relative (pure H2O: 1.2e-8); everywhere else the 1e-10 tolerance applies.
            const tol = T === 1000 ? 3.5e-8 : 1e-10;
            expect(e, `${c.name} ${name} T=${T} guess=${g}`).toBeLessThan(tol);
          }
        }
      }
    }
    console.info(
      `[mixture] inverse round-trip max rel. error ${worst.toExponential(2)} ` +
        `(${worstTmid.toExponential(2)} exactly at Tmid), max iterations ${maxIt}`,
    );
  });

  it('seeded random sweep: any composition, 60–15000 K, 0.1 kPa–1 GPa, any guess → ≤ 10 iterations, 1e-12', () => {
    // Regression test for a stall: when the final Newton correction dropped below one ulp,
    // T + ΔT == T == the bracket end just set, the step was rejected and ~35 bisection
    // steps followed (≈ 1 % of calls, even from warm starts; results only 5e-11 accurate).
    let seed = 12345;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const X = new Float64Array(NS);
    let worst = 0;
    let maxIt = 0;
    let maxItWarm = 0;
    for (let trial = 0; trial < 4000; trial++) {
      X.fill(0);
      const nsp = 1 + Math.floor(rnd() * NS);
      for (let j = 0; j < nsp; j++) X[Math.floor(rnd() * NS)] = Math.pow(10, -12 * rnd());
      if (rnd() < 0.2) X[Math.floor(rnd() * NS)] = -1e-18; // integrator round-off
      let sx = 0;
      for (let k = 0; k < NS; k++) sx += X[k];
      for (let k = 0; k < NS; k++) X[k] /= sx;
      const T = Math.exp(Math.log(60) + rnd() * Math.log(15000 / 60));
      const p = Math.exp(Math.log(1e2) + rnd() * Math.log(1e9 / 1e2));
      const guesses = [NaN, -1, 0, Infinity, 1e6, 300, 3000];
      const g = guesses[Math.floor(rnd() * guesses.length)];
      const warm = T * (1 + 0.1 * (rnd() - 0.5));
      for (const guess of [g, warm]) {
        const res = [
          temperatureFromH(X, mixHMass(X, T), guess),
          temperatureFromU(X, mixUMass(X, T), guess),
          temperatureFromS(X, mixSMass(X, T, p), p, guess),
        ];
        for (const Tx of res) worst = Math.max(worst, rel(Tx, T));
        // (lastInverse reflects the S call; check all three)
        for (const f of [
          () => temperatureFromH(X, mixHMass(X, T), guess),
          () => temperatureFromU(X, mixUMass(X, T), guess),
          () => temperatureFromS(X, mixSMass(X, T, p), p, guess),
        ]) {
          f();
          expect(lastInverse.converged).toBe(true);
          maxIt = Math.max(maxIt, lastInverse.iterations);
          if (guess === warm) maxItWarm = Math.max(maxItWarm, lastInverse.iterations);
        }
      }
    }
    console.info(`[mixture] random inverse sweep: max rel. error ${worst.toExponential(2)}, max iterations ${maxIt} (warm start ${maxItWarm})`);
    expect(worst).toBeLessThan(1e-12);
    expect(maxIt).toBeLessThanOrEqual(10);
    expect(maxItWarm).toBeLessThanOrEqual(5);
  });

  it('reports convergence and typical iteration counts with a warm start', () => {
    const X = cases[4].Xa;
    const T = 2345.6;
    temperatureFromU(X, mixUMass(X, T), 2300);
    expect(lastInverse.converged).toBe(true);
    expect(lastInverse.iterations).toBeLessThanOrEqual(5);
  });

  it('molar (extensive) variants solve U_total(T) = U for a mole vector', () => {
    const X = cases[1].Xa;
    const N = X.map((x) => x * 0.004);
    const T = 812.3;
    expect(rel(temperatureFromUMolar(N, mixUMolar(N, T), 400), T)).toBeLessThan(1e-10);
    expect(rel(temperatureFromHMolar(N, mixHMolar(N, T), 400), T)).toBeLessThan(1e-10);
  });

  it('returns the bound (converged = false) for unreachable targets', () => {
    const X = cases[0].Xa;
    const T = temperatureFromH(X, -1e12, 300);
    expect(T).toBeLessThanOrEqual(1.0000001);
    expect(lastInverse.converged).toBe(false);
  });
});
