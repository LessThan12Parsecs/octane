import { describe, expect, it } from 'vitest';
import { N_EQ, NE, NS, SP } from '../core/species';
import { dryAir, freshCharge, fuelFromSelection, humidAir, type FuelBlend } from '../thermo/fuels';
import { elementMoles, mixHMolar, mixUMolar, totalMoles } from '../thermo/mixture';
import { R_UNIVERSAL } from '../core/constants';
import { EquilibriumSolver } from './solver';
import { elementBalanceError as elementError, stationarityError } from './diagnostics';

/** Mulberry32 PRNG (deterministic, seeded). */
function rng(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const logUniform = (u: number, lo: number, hi: number): number => lo * Math.pow(hi / lo, u);

const FUELS: FuelBlend[] = [
  fuelFromSelection({ kind: 'pure', species: 'IC8H18' }),
  fuelFromSelection({ kind: 'pure', species: 'NC7H16' }),
  fuelFromSelection({ kind: 'pure', species: 'CH4' }),
  fuelFromSelection({ kind: 'pure', species: 'C3H8' }),
  fuelFromSelection({ kind: 'pure', species: 'C2H5OH' }),
  fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }),
  fuelFromSelection({ kind: 'PRF', octaneNumber: 60 }),
];

/** Random reactant mixture (mole fractions) covering φ 0.2–3, humidity, EGR, absent elements. */
function randomReactants(r: () => number): { X: Float64Array; kind: string } {
  const u = r();
  const phi = logUniform(r(), 0.2, 3);
  if (u < 0.06) {
    // hydrogen in N2/O2 air: no carbon, no argon
    const X = new Float64Array(NS);
    const o2 = 0.21;
    X[SP.O2] = o2;
    X[SP.N2] = 0.79;
    X[SP.H2] = 2 * o2 * phi;
    return { X, kind: 'H2-air' };
  }
  if (u < 0.1) {
    // CO in dry air (no hydrogen)
    const X = dryAir();
    X[SP.CO] = 2 * X[SP.O2] * phi;
    return { X, kind: 'CO-air' };
  }
  if (u < 0.13) return { X: humidAir(300, 1e5, r()), kind: 'air' };
  if (u < 0.16) {
    // H2/O2 only (no N, C, Ar)
    const X = new Float64Array(NS);
    X[SP.O2] = 1;
    X[SP.H2] = 2 * phi;
    return { X, kind: 'H2-O2' };
  }
  const fuel = FUELS[Math.floor(r() * FUELS.length)];
  const airX = r() < 0.5 ? dryAir() : humidAir(280 + 40 * r(), 1e5, r());
  const egr = r() < 0.4 ? 0.3 * r() : 0;
  return { X: freshCharge({ fuel, phi, airX, egrFraction: egr }), kind: `${fuel.label} φ=${phi.toFixed(2)} egr=${egr.toFixed(2)}` };
}

describe('robustness: seeded random sweep', () => {
  it('2400 random (b, T, p) points, φ 0.2–3, 300–4000 K, 0.1–300 bar: zero non-convergence', () => {
    const r = rng(20260929);
    const warm = new EquilibriumSolver();
    const cold = new EquilibriumSolver();
    const b = new Float64Array(NE);
    const Xw = new Float64Array(NS);
    const n = 2400;
    let fails = 0;
    let worstBal = 0;
    let worstStat = 0;
    let worstDiff = 0;
    let itsWarm = 0;
    let itsCold = 0;
    let maxCold = 0;
    let maxWarm = 0;
    for (let i = 0; i < n; i++) {
      const { X, kind } = randomReactants(r);
      elementMoles(X, b);
      const scale = logUniform(r(), 1e-8, 1e2); // mol-scale of the burned zone
      for (let e = 0; e < NE; e++) b[e] *= scale;
      const T = logUniform(r(), 300, 4000);
      const p = logUniform(r(), 1e4, 3e7);
      const rw = warm.solveTP(b, T, p);
      if (!rw.converged) {
        fails++;
        console.log('warm fail', kind, T, p);
      }
      itsWarm += rw.iterations;
      maxWarm = Math.max(maxWarm, rw.iterations);
      worstBal = Math.max(worstBal, elementError(rw, b));
      worstStat = Math.max(worstStat, stationarityError(rw));
      Xw.set(rw.X);
      cold.reset();
      const rc = cold.solveTP(b, T, p);
      if (!rc.converged) {
        fails++;
        console.log('cold fail', kind, T, p);
      }
      itsCold += rc.iterations;
      maxCold = Math.max(maxCold, rc.iterations);
      worstBal = Math.max(worstBal, elementError(rc, b));
      for (let k = 0; k < N_EQ; k++) {
        const ref = rc.X[k];
        if (ref > 1e-3) worstDiff = Math.max(worstDiff, Math.abs(Xw[k] - ref) / 1e-3);
        else if (ref > 1e-12) worstDiff = Math.max(worstDiff, Math.abs(Xw[k] / ref - 1));
      }
    }
    console.log(
      `random sweep n=${n}: failures ${fails}; element balance ${worstBal.toExponential(2)}; ` +
        `stationarity ${worstStat.toExponential(2)}; warm-vs-cold ${worstDiff.toExponential(2)}; ` +
        `iterations warm mean ${(itsWarm / n).toFixed(2)} max ${maxWarm}, cold mean ${(itsCold / n).toFixed(2)} max ${maxCold}`,
    );
    expect(fails).toBe(0);
    expect(worstBal).toBeLessThan(1e-12);
    expect(worstStat).toBeLessThan(1e-8);
    expect(worstDiff).toBeLessThan(1e-6);
  });

  it('400 random HP and UV flames from 300–900 K, 0.5–60 bar reactants: all converge, energy conserved', () => {
    const r = rng(7);
    const s = new EquilibriumSolver();
    const b = new Float64Array(NE);
    let worstE = 0;
    let outerMax = 0;
    for (let i = 0; i < 400; i++) {
      const { X } = randomReactants(r);
      elementMoles(X, b);
      const T0 = 300 + 600 * r();
      const p0 = logUniform(r(), 0.5e5, 60e5);
      const H0 = mixHMolar(X, T0);
      const rh = s.solveHP(b, H0, p0, 300 + 3000 * r());
      expect(rh.converged).toBe(true);
      worstE = Math.max(worstE, Math.abs(rh.H - H0) / (rh.nTotal * R_UNIVERSAL * rh.T));
      outerMax = Math.max(outerMax, rh.outerIterations);
      const U0 = mixUMolar(X, T0);
      const V0 = (totalMoles(X) * R_UNIVERSAL * T0) / p0;
      const ru = s.solveUV(b, U0, V0, 300 + 3000 * r());
      expect(ru.converged).toBe(true);
      worstE = Math.max(worstE, Math.abs(ru.U - U0) / (ru.nTotal * R_UNIVERSAL * ru.T));
      expect(Math.abs(ru.V / V0 - 1)).toBeLessThan(1e-12);
      outerMax = Math.max(outerMax, ru.outerIterations);
    }
    console.log(`random HP/UV: max |ΔE|/(nRT) = ${worstE.toExponential(2)}, max outer iterations ${outerMax}`);
    expect(worstE).toBeLessThan(1e-9);
  });

  it('elements present only in traces (1e-18 … 1e-8 of the atoms) are conserved to 1e-12 relative', () => {
    const s = new EquilibriumSolver();
    const base = elementMoles(freshCharge({ fuel: FUELS[0], phi: 0.9, airX: dryAir() }), new Float64Array(NE));
    const b = new Float64Array(NE);
    let tot = 0;
    for (let e = 0; e < NE; e++) tot += base[e];
    for (const e of [0, 3, 4]) {
      // C, N, Ar in traces (C needs O, which is always plentiful here)
      for (const frac of [1e-18, 1e-14, 1e-10, 1e-8]) {
        b.set(base);
        b[e] = frac * tot;
        for (const T of [300, 1200, 2500, 3800]) {
          for (const p of [1e4, 1e6, 3e7]) {
            s.reset();
            const r = s.solveTP(b, T, p);
            expect(r.converged).toBe(true);
            expect(elementError(r, b)).toBeLessThan(1e-12);
            expect(stationarityError(r)).toBeLessThan(1e-8);
            expect(s.solveTP(b, T * 1.02, p * 0.97).converged).toBe(true);
          }
        }
      }
    }
    // below EQ_B_REL_MIN (1e-20) the element is dropped with its species
    b.set(base);
    b[0] = 1e-22 * tot;
    const r = s.solveTP(b, 2000, 1e5);
    expect(r.converged).toBe(true);
    expect(r.N[SP.CO2]).toBe(0);
    expect(r.N[SP.CO]).toBe(0);
  });

  it('exactly stoichiometric mixtures at 300–700 K (O2/CO/H2 at the round-off level) converge', () => {
    const s = new EquilibriumSolver();
    const b = new Float64Array(NE);
    for (const fuel of FUELS) {
      const X = freshCharge({ fuel, phi: 1, airX: dryAir() });
      elementMoles(X, b);
      for (const T of [300, 350, 400, 500, 700]) {
        for (const p of [1e4, 1e5, 1e7, 3e7]) {
          s.reset();
          const rc = s.solveTP(b, T, p);
          expect(rc.converged).toBe(true);
          expect(elementError(rc, b)).toBeLessThan(1e-12);
          expect(s.solveTP(b, T * 1.003, p * 1.01).converged).toBe(true);
        }
      }
    }
  });
});

// ---------------------------------------------------------------------------------------------
// Warm-start sequences along a synthetic expansion stroke
// ---------------------------------------------------------------------------------------------

/** Slider-crank volume ratio V/Vc for crank angle θ (rad after TDC), CR 8, rod/crank 4.44 (CFR-like). */
function volumeRatio(theta: number, cr = 8, lr = 4.44): number {
  const s = 1 - Math.cos(theta) + lr - Math.sqrt(lr * lr - Math.sin(theta) ** 2);
  return 1 + ((cr - 1) / 2) * s;
}

function burnedB(): Float64Array {
  const X = freshCharge({ fuel: fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }), phi: 1, airX: humidAir(298, 1e5, 0.5), egrFraction: 0.08 });
  const b = elementMoles(X, new Float64Array(NE));
  for (let e = 0; e < NE; e++) b[e] *= 0.012; // ≈ burned-zone moles of a CFR charge
  return b;
}

describe('warm starts along an expansion stroke', () => {
  const b = burnedB();
  // polytropic expansion from 10° ATDC to 150° ATDC in 0.25° steps
  const thetas: number[] = [];
  for (let d = 10; d <= 150; d += 0.25) thetas.push((d * Math.PI) / 180);
  const V0 = volumeRatio(thetas[0]);
  const T0 = 2750;
  const p0 = 65e5;
  const Tof = (th: number): number => T0 * Math.pow(V0 / volumeRatio(th), 0.28);
  const pof = (th: number): number => p0 * Math.pow(V0 / volumeRatio(th), 1.28);

  it('TP: warm solves match cold solves, few iterations per step', () => {
    const warm = new EquilibriumSolver();
    const cold = new EquilibriumSolver();
    let its = 0;
    let maxIts = 0;
    let worst = 0;
    warm.solveTP(b, Tof(thetas[0]) * 1.01, pof(thetas[0])); // prime the warm start
    for (const th of thetas) {
      const T = Tof(th);
      const p = pof(th);
      const rw = warm.solveTP(b, T, p);
      expect(rw.converged).toBe(true);
      its += rw.iterations;
      maxIts = Math.max(maxIts, rw.iterations);
      const Xw = Float64Array.from(rw.X);
      cold.reset();
      const rc = cold.solveTP(b, T, p);
      for (let k = 0; k < N_EQ; k++) if (rc.X[k] > 1e-12) worst = Math.max(worst, Math.abs(Xw[k] / rc.X[k] - 1));
    }
    console.log(`expansion TP (${thetas.length} steps of 0.25°): mean iterations ${(its / thetas.length).toFixed(2)}, max ${maxIts}, warm vs cold rel ${worst.toExponential(2)}`);
    expect(maxIts).toBeLessThanOrEqual(4);
    expect(worst).toBeLessThan(1e-8);
  });

  it('UV: burned zone expanding with p dV work, warm-started (energy bookkeeping exact)', () => {
    const s = new EquilibriumSolver();
    const Vc = 5e-5; // m³, clearance-volume scale
    let V = Vc * V0;
    let r = s.solveTP(b, T0, p0);
    // start from the equilibrium state at (T0, p0) with its volume
    V = r.V;
    const scaleV = V / V0;
    let U = r.U;
    let T = r.T;
    let p = r.p;
    let its = 0;
    let outer = 0;
    let maxIts = 0;
    for (let i = 1; i < thetas.length; i++) {
      const Vn = scaleV * volumeRatio(thetas[i]);
      U -= p * (Vn - V); // explicit work extraction (the cycle model does this with RK)
      V = Vn;
      r = s.solveUV(b, U, V, T);
      expect(r.converged).toBe(true);
      expect(Math.abs(r.U - U)).toBeLessThan(1e-9 * r.nTotal * R_UNIVERSAL * r.T);
      its += r.iterations;
      outer += r.outerIterations;
      maxIts = Math.max(maxIts, r.iterations);
      T = r.T;
      p = r.p;
    }
    console.log(`expansion UV: mean inner iterations ${(its / thetas.length).toFixed(2)} (max ${maxIts}), mean outer ${(outer / thetas.length).toFixed(2)}; end T ${T.toFixed(1)} K, p ${(p / 1e5).toFixed(2)} bar`);
    expect(T).toBeLessThan(T0);
    expect(maxIts).toBeLessThanOrEqual(12);
  });
});

// ---------------------------------------------------------------------------------------------
// Timing (reported; bounds are generous to stay robust on slow CI machines)
// ---------------------------------------------------------------------------------------------

describe('timing', () => {
  it('µs per warm-started TP / UV / HP solve and per cold TP solve', () => {
    const b = burnedB();
    const steps: { T: number; p: number; V: number; U: number }[] = [];
    const s = new EquilibriumSolver();
    for (let d = 10; d <= 150; d += 0.25) {
      const th = (d * Math.PI) / 180;
      const T = 2750 * Math.pow(volumeRatio(0.1745) / volumeRatio(th), 0.28);
      const p = 65e5 * Math.pow(volumeRatio(0.1745) / volumeRatio(th), 1.28);
      const r = s.solveTP(b, T, p);
      steps.push({ T, p, V: r.V, U: r.U });
    }
    const reps = 40;
    const time = (f: () => void): number => {
      f(); // warm-up (JIT)
      const t0 = performance.now();
      for (let k = 0; k < reps; k++) f();
      return ((performance.now() - t0) * 1000) / (reps * steps.length);
    };
    let sink = 0;
    const tp = time(() => {
      for (const st of steps) sink += s.solveTP(b, st.T, st.p).X[SP.NO];
    });
    const uv = time(() => {
      for (const st of steps) sink += s.solveUV(b, st.U, st.V, st.T * 1.001).T;
    });
    const hp = time(() => {
      for (const st of steps) sink += s.solveHP(b, st.U + st.p * st.V, st.p, st.T * 1.001).T;
    });
    const c = new EquilibriumSolver();
    const cold = time(() => {
      for (const st of steps) {
        c.reset();
        sink += c.solveTP(b, st.T, st.p).T;
      }
    });
    console.log(
      `timing: warm TP ${tp.toFixed(2)} µs, warm UV ${uv.toFixed(2)} µs, warm HP ${hp.toFixed(2)} µs, cold TP ${cold.toFixed(2)} µs per solve`,
    );
    expect(sink).not.toBeNaN();
    expect(tp).toBeLessThan(100);
    expect(uv).toBeLessThan(400);
  });
});
