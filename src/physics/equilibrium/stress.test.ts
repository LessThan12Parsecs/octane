/**
 * Adversarial-review tests: states at and beyond the contract envelope (Cantera oracle
 * tools/reference/equilibrium_stress.py), invalid inputs, the textbook flame temperatures,
 * smoothness for the ODE integrator and hot-path allocation.
 */
import { getHeapStatistics } from 'node:v8';
import { describe, expect, it } from 'vitest';
import { P_ATM, T_REF } from '../core/constants';
import { EL, N_EQ, NE, NS, SP } from '../core/species';
import { dryAir, freshCharge, fuelFromSelection } from '../thermo/fuels';
import { elementMoles } from '../thermo/mixture';
import stFx from '../../../test/fixtures/equilibrium_stress.json';
import { elementBalanceError, stationarityError } from './diagnostics';
import { adiabaticFlameTemperature } from './flame';
import { equilibriumProperties } from './properties';
import { EQ_B_NEG_TOL, EquilibriumSolver, newEqProperties } from './solver';

const ABS_MAJOR = 1e-7;
const REL_MINOR = 1e-4;

function xErr(X: Float64Array, ref: readonly number[], skip: readonly number[] = []): [number, number] {
  let a = 0;
  let r = 0;
  for (let k = 0; k < N_EQ; k++) {
    if (skip.includes(k)) continue;
    const v = ref[k];
    if (v > 1e-3) a = Math.max(a, Math.abs(X[k] - v));
    if (v >= 1e-12) r = Math.max(r, Math.abs(X[k] - v) / v);
  }
  return [a, r];
}

const bIso = (phi: number): Float64Array =>
  elementMoles(
    freshCharge({ fuel: fuelFromSelection({ kind: 'pure', species: 'IC8H18' }), phi, airX: dryAir() }),
    new Float64Array(NE),
  );

describe('stress fixtures vs Cantera (equilibrium_stress.json)', () => {
  const comps = stFx.compositions.map((c) => Float64Array.from(c.b));

  it(`TP: ${stFx.tp.length} states, 200–6000 K, 1 kPa–1000 bar, φ 0.1…O/C limit, 50 % EGR, Ar-diluted, air`, () => {
    const cold = new EquilibriumSolver();
    const warm = new EquilibriumSolver();
    let wa = 0;
    let wr = 0;
    let wBal = 0;
    let wStat = 0;
    for (const c of stFx.tp) {
      const b = comps[c.c];
      const skip = (c as { ill?: number[] }).ill ?? [];
      cold.reset();
      for (const s of [cold, warm]) {
        const r = s.solveTP(b, c.T, c.p);
        expect(r.converged, `${stFx.compositions[c.c].label} T=${c.T} p=${c.p}`).toBe(true);
        const [a, rel] = xErr(r.X, c.X, skip);
        wa = Math.max(wa, a);
        wr = Math.max(wr, rel);
        wBal = Math.max(wBal, elementBalanceError(r, b));
        wStat = Math.max(wStat, stationarityError(r));
      }
    }
    console.log(`stress TP: max|ΔX| ${wa.toExponential(2)}, rel ${wr.toExponential(2)}, balance ${wBal.toExponential(2)}, stationarity ${wStat.toExponential(2)}`);
    expect(wa).toBeLessThanOrEqual(ABS_MAJOR);
    expect(wr).toBeLessThanOrEqual(REL_MINOR);
    expect(wBal).toBeLessThan(1e-12);
    expect(wStat).toBeLessThan(1e-8);
  });

  it('HP and UV flames from harsh start states (end gas 900 K/100 bar, φ 0.3 cold, φ 2.5, Ar, 50 % EGR)', () => {
    const s = new EquilibriumSolver();
    let wT = 0;
    let wP = 0;
    let wa = 0;
    let wr = 0;
    for (const [mode, list] of [['HP', stFx.hp], ['UV', stFx.uv]] as const) {
      for (const c of list) {
        const f = adiabaticFlameTemperature(Float64Array.from(c.Xr), c.T0, c.p0, mode, s);
        expect(f.converged, `${mode} ${c.label}`).toBe(true);
        wT = Math.max(wT, Math.abs(f.T - c.T));
        if (mode === 'UV') wP = Math.max(wP, Math.abs(f.p / (c as { p: number }).p - 1));
        const [a, r] = xErr(f.X, c.X);
        wa = Math.max(wa, a);
        wr = Math.max(wr, r);
      }
    }
    console.log(`stress flames: max|ΔT| ${wT.toExponential(2)} K, UV p ${wP.toExponential(2)}, |ΔX| ${wa.toExponential(2)}, rel ${wr.toExponential(2)}`);
    expect(wT).toBeLessThan(0.1);
    expect(wP).toBeLessThan(1e-7);
    expect(wa).toBeLessThanOrEqual(ABS_MAJOR);
    expect(wr).toBeLessThanOrEqual(REL_MINOR);
  });

  it('12-species truncation vs 27 NASA species at engine states is small (oracle-side bound)', () => {
    // Not a test of our code: documents that HO2, NO2, N2O, H2O2, … stay ≲ 1e-5 and shift the
    // majors/NO by < 0.2 % at 1800–2900 K, 5–100 bar (iso-octane–air φ = 1).
    for (const r of stFx.truncation) {
      for (const k of Object.keys(r.x12) as (keyof typeof r.x12)[]) {
        expect(Math.abs(r.x27[k] / r.x12[k] - 1), `${k} at ${r.T} K`).toBeLessThan(2e-3);
      }
      expect(r.x27.HO2 + r.x27.NO2 + r.x27.N2O + r.x27.H2O2).toBeLessThan(5e-5);
    }
  });
});

describe('textbook adiabatic flame temperatures (fuel + O2 + 3.76 N2, 298.15 K, 1 atm)', () => {
  // Turns, "An Introduction to Combustion" (2nd ed., 2000), Appendix B, Table B.1 lists the
  // constant-pressure equilibrium T_ad of stoichiometric fuel–air: CH4 2226 K, C3H8 2267 K,
  // n-C7H16 2274 K, C8H18 2275 K (the C8H18 row is n-octane, h_f = −208 447 kJ/kmol).
  // UNVERIFIED: the table could not be fetched; values from memory. Corroboration fetched
  // 2026-09-29: NIST WebBook ΔfH°gas n-octane −208.4 ± 0.67 kJ/mol (Prosen & Rossini 1945),
  // 2,2,4-trimethylpentane −224.1 ± 1.3 kJ/mol (same source) — so iso-octane burns ≈ 3.7 K cooler
  // than n-octane; Cantera + GRI-Mech 3.0 gives 2224.25 K for CH4–air (arXiv:2503.11826).
  // n-octane is not one of our species: the oracle hands its (b, H) from the NASA fit.
  const turns: Record<string, number> = { CH4: 2226, C3H8: 2267, 'n-C7H16': 2274, 'n-C8H18 (not in SPECIES)': 2275 };

  it('matches the Cantera oracle to 0.1 K and Turns Table B.1 to ≤ 1.1 K', () => {
    const s = new EquilibriumSolver();
    for (const c of stFx.textbook) {
      const r = s.solveHP(Float64Array.from(c.b), c.H, c.p, 2000);
      expect(r.converged, c.label).toBe(true);
      expect(Math.abs(r.T - c.T), c.label).toBeLessThan(0.1);
      // extended (27-species) set moves T_ad by < 0.02 K here
      expect(Math.abs(c.T_extended_species - c.T)).toBeLessThan(0.02);
      const ref = turns[c.label];
      console.log(`${c.label}: T_ad = ${r.T.toFixed(2)} K (Cantera ${c.T.toFixed(2)}, Turns ${ref ?? '—'})`);
      if (ref !== undefined) expect(Math.abs(r.T - ref), c.label).toBeLessThan(1.1);
    }
  });

  it('iso-octane (our species) burns 3–4.5 K cooler than n-octane (ΔfH° 15.7 kJ/mol more negative)', () => {
    const X = new Float64Array(NS);
    X[SP.IC8H18] = 1;
    X[SP.O2] = 12.5;
    X[SP.N2] = 47;
    const iso = adiabaticFlameTemperature(X, T_REF, P_ATM, 'HP');
    const nOct = stFx.textbook.find((c) => c.label.startsWith('n-C8H18'))!.T;
    const isoRef = stFx.textbook.find((c) => c.label.startsWith('iso-C8H18'))!.T;
    expect(Math.abs(iso.T - isoRef)).toBeLessThan(0.1); // our fuel thermo + solver = oracle
    expect(nOct - iso.T).toBeGreaterThan(3);
    expect(nOct - iso.T).toBeLessThan(4.5);
  });
});

describe('invalid inputs fail loudly (converged = false, NaN X), never "converge"', () => {
  const b = bIso(1);
  const s = new EquilibriumSolver();
  const bad = (r: { converged: boolean; X: Float64Array }): void => {
    expect(r.converged).toBe(false);
    expect(Number.isNaN(r.X[SP.N2])).toBe(true);
  };

  it('NaN / ±Inf / non-physically negative element amounts (were silently dropped → converged)', () => {
    for (const [e, v] of [[EL.H, NaN], [EL.N, Infinity], [EL.O, -Infinity], [EL.C, -1e-3], [EL.H, -0.5]] as const) {
      const bb = Float64Array.from(b);
      bb[e] = v;
      bad(s.solveTP(bb, 2000, 1e5));
      bad(s.solveTV(bb, 2000, 1e-3));
      bad(s.solveHP(bb, -1e5, 1e5, 2000));
      bad(s.solveUV(bb, -1e5, 1e-3, 2000));
      expect(s.properties().valid).toBe(false);
      expect(s.speciesDerivatives(new Float64Array(NS), new Float64Array(NS))).toBe(false);
    }
    // round-off-level negatives are treated as absent (like b_e/Σb < EQ_B_REL_MIN)
    const bb = Float64Array.from(b);
    let tot = 0;
    for (let e = 0; e < NE; e++) tot += bb[e];
    bb[EL.AR] = -0.5 * EQ_B_NEG_TOL * tot;
    const r = s.solveTP(bb, 2000, 1e5);
    expect(r.converged).toBe(true);
    expect(r.N[SP.AR]).toBe(0);
  });

  it('T, p, V not finite and positive', () => {
    for (const T of [NaN, 0, -5, Infinity]) bad(s.solveTP(b, T, 1e5));
    for (const p of [NaN, 0, -1, Infinity]) bad(s.solveTP(b, 2000, p));
    for (const V of [NaN, 0, -1, Infinity]) {
      bad(s.solveTV(b, 2000, V));
      bad(s.solveUV(b, -1e5, V, 2000));
    }
    for (const p of [NaN, 0, -1]) bad(s.solveHP(b, -1e5, p, 2000));
    // and a valid solve afterwards still works (warm start kept)
    expect(s.solveTP(b, 2000, 1e5).converged).toBe(true);
  });

  it('NaN/∞ H, U, S targets (a NaN target used to "converge" at 6000 K, −∞ at 200 K)', () => {
    for (const t of [NaN, Infinity, -Infinity]) {
      bad(s.solveHP(b, t, 1e5, 2000));
      bad(s.solveUV(b, t, 1e-3, 2000));
      bad(s.solveSP(b, t, 1e5, 2000));
    }
  });
});

describe('robustness beyond the contract envelope', () => {
  it('φ 0.1–3.05 (O/C → 1.02), 200–6000 K every 50 K, 100 Pa – 1000 bar: zero failures, warm and cold', () => {
    const warm = new EquilibriumSolver();
    const cold = new EquilibriumSolver();
    let n = 0;
    let wBal = 0;
    let maxIt = 0;
    for (const phi of [0.1, 0.5, 1, 1.5, 2.5, 3.05]) {
      const b = bIso(phi);
      for (let T = 200; T <= 6000; T += 50) {
        for (const p of [100, 1e4, 1e6, 3e7, 1e8]) {
          const rw = warm.solveTP(b, T, p);
          expect(rw.converged).toBe(true);
          wBal = Math.max(wBal, elementBalanceError(rw, b));
          cold.reset();
          const rc = cold.solveTP(b, T, p);
          expect(rc.converged).toBe(true);
          maxIt = Math.max(maxIt, rw.iterations, rc.iterations);
          n++;
        }
      }
    }
    console.log(`extended sweep: ${n} states × (warm, cold), element balance ${wBal.toExponential(2)}, max iterations ${maxIt}`);
    expect(wBal).toBeLessThan(1e-12);
  });

  it('single-element and near-limit systems (pure Ar, N, O, H; O/C = 1 + 1e-9)', () => {
    const s = new EquilibriumSolver();
    for (const b of [[0, 0, 0, 0, 1], [0, 0, 0, 1, 0], [0, 0, 1, 0, 0], [0, 1, 0, 0, 0], [1, 0, 1 + 1e-9, 0, 0], [1, 2.25, 1.0001, 7.52, 0]]) {
      const bb = Float64Array.from(b);
      for (const T of [200, 1000, 3000, 6000]) {
        for (const p of [1e3, 1e5, 3e7]) {
          s.reset();
          const r = s.solveTP(bb, T, p);
          expect(r.converged, `${b} ${T} ${p}`).toBe(true);
          expect(elementBalanceError(r, bb)).toBeLessThan(1e-12);
        }
      }
    }
  });

  it('HP/UV from any guess in and outside [200, 6000] K, including targets near 200 K', () => {
    const s = new EquilibriumSolver();
    const b = bIso(1);
    for (const Tt of [210, 2800]) {
      const r0 = s.solveTP(b, Tt, 60e5);
      const { U, V, H } = r0;
      for (const Tg of [150, 200, 201, 300, 1000, 5999, 6000, 7000, NaN]) {
        s.reset();
        expect(Math.abs(s.solveUV(b, U, V, Tg).T - Tt)).toBeLessThan(1e-6);
        expect(s.result.converged).toBe(true);
        s.reset();
        expect(Math.abs(s.solveHP(b, H, 60e5, Tg).T - Tt)).toBeLessThan(1e-6);
        expect(s.result.converged).toBe(true);
      }
    }
  });
});

describe('smoothness for the cycle integrator', () => {
  it('ln X_k(T) and H(T) along a 0.5 K grid have no noise above the analytic third difference', () => {
    const s = new EquilibriumSolver();
    const b = bIso(1);
    const lnX: Float64Array[] = [];
    const H: number[] = [];
    let cpT = 0;
    const out = newEqProperties();
    for (let i = 0; i <= 600; i++) {
      const T = 2000 + 0.5 * i;
      const r = s.solveTP(b, T, 30e5);
      lnX.push(Float64Array.from(r.X.subarray(0, N_EQ), Math.log));
      H.push(r.H);
      if (i === 300) cpT = s.properties(out).cp * r.mass * T;
    }
    // third differences of smooth functions over h = 0.5 K are ~1e-9 (ln X) and ~1e-9·cpT (H);
    // Newton-tolerance noise would show as ≥ 1e-8 spikes
    let w = 0;
    let wH = 0;
    for (let i = 3; i < lnX.length; i++) {
      for (let k = 0; k < N_EQ; k++) {
        w = Math.max(w, Math.abs(lnX[i][k] - 3 * lnX[i - 1][k] + 3 * lnX[i - 2][k] - lnX[i - 3][k]));
      }
      wH = Math.max(wH, Math.abs(H[i] - 3 * H[i - 1] + 3 * H[i - 2] - H[i - 3]) / cpT);
    }
    console.log(`smoothness: max third difference ln X ${w.toExponential(2)}, H/(cp T) ${wH.toExponential(2)}`);
    expect(w).toBeLessThan(1e-8);
    expect(wH).toBeLessThan(1e-9);
  });

  it('T(U) of warm-started UV solves is smooth to ~1e-10 K and path-independent (warm vs cold)', () => {
    // The outer Newton stops at |ΔT| ≤ 1e-6 T and finishes with the first-order composition
    // predictor; the returned T must still be a smooth function of U (the integrator differences it).
    const s = new EquilibriumSolver();
    const b = bIso(1);
    const r0 = s.solveTP(b, 2600, 50e5);
    const U0 = r0.U;
    const V0 = r0.V;
    const T: number[] = [];
    let Tg = 2600;
    for (let i = 0; i < 300; i++) T.push((Tg = s.solveUV(b, U0 + 0.003 * i, V0, Tg).T)); // ΔT ≈ 6e-5 K/step
    let w = 0;
    for (let i = 3; i < T.length; i++) w = Math.max(w, Math.abs(T[i] - 3 * T[i - 1] + 3 * T[i - 2] - T[i - 3]));
    let wc = 0;
    const c = new EquilibriumSolver();
    for (let i = 0; i < T.length; i += 23) wc = Math.max(wc, Math.abs(c.solveUV(b, U0 + 0.003 * i, V0, 1500).T - T[i]));
    console.log(`UV smoothness: third difference of T ${w.toExponential(2)} K, warm vs cold ${wc.toExponential(2)} K`);
    expect(w).toBeLessThan(1e-9);
    expect(wc).toBeLessThan(1e-9);
  });
});

describe('hot-path allocation (V8 total_allocated_bytes)', () => {
  // The solver allocates no arrays/objects per call; what remains is V8 boxing doubles passed
  // across non-inlined calls (≈ 1 heap number per warm TP solve, ~10–15 per UV/HP solve). The
  // bounds catch a regression that allocates a typed array or object per call (≥ ~150 B each).
  it('warm TP < 100 B/call, warm UV/HP < 600 B/call, TP+properties < 150 B/call', () => {
    const b = bIso(1);
    for (let e = 0; e < NE; e++) b[e] *= 0.012;
    const s = new EquilibriumSolver();
    const out = newEqProperties();
    const r0 = s.solveTP(b, 2600, 50e5);
    const U0 = r0.U;
    const V0 = r0.V;
    const H0 = r0.H;
    const run = (n: number, kind: number): number => {
      let sink = 0;
      let Tg = 2600;
      for (let i = 0; i < n; i++) {
        const j = i % 400;
        if (kind === 0) sink += s.solveTP(b, 2600 - 2 * j, 50e5 - 1000 * j).T;
        else if (kind === 1) sink += Tg = s.solveUV(b, U0 - 0.2 * j, V0 * (1 + 0.004 * j), Tg).T;
        else if (kind === 2) sink += Tg = s.solveHP(b, H0 - 0.2 * j, 50e5, Tg).T;
        else sink += equilibriumProperties(s, b, 2600 - 2 * j, 50e5, out).cp;
      }
      return sink;
    };
    const perCall: number[] = [];
    for (let kind = 0; kind < 4; kind++) {
      for (let w = 0; w < 3; w++) run(10000, kind); // JIT warm-up (interpreter frames box doubles)
      const a0 = getHeapStatistics().total_allocated_bytes;
      const n = 20000;
      expect(run(n, kind)).toBeGreaterThan(0);
      perCall.push((getHeapStatistics().total_allocated_bytes - a0) / n);
    }
    console.log(`allocation B/call: warm TP ${perCall[0].toFixed(1)}, UV ${perCall[1].toFixed(1)}, HP ${perCall[2].toFixed(1)}, TP+props ${perCall[3].toFixed(1)}`);
    expect(perCall[0]).toBeLessThan(100);
    expect(perCall[1]).toBeLessThan(600);
    expect(perCall[2]).toBeLessThan(600);
    expect(perCall[3]).toBeLessThan(150);
  });
});

