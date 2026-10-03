import { describe, expect, it } from 'vitest';
import fx from '../../../test/fixtures/chemistry_zeldovich.json';
import { R_UNIVERSAL } from '../core/constants';
import { N_EQ, NS, SP } from '../core/species';
import {
  advanceAlpha,
  ZELDOVICH_GRI30,
  ZELDOVICH_HEYWOOD,
  ZeldovichKinetics,
  zeldovichNORate,
} from './zeldovich';

interface State { T: number; p: number; phi: number; Xeq: number[] }
interface Case {
  s: number; alpha: number; wNOqss: number; wNqss: number; xNqss: number;
  wNOfull: number; wNfull: number; wNOgriThermo: number; wNOheywood: number;
}
interface Evo {
  T: number; p: number; phi: number; Xeq: number[]; tauNO: number; times: number[];
  xNOeGri: number; griNO: number[]; griAlpha: number[]; modelAlpha: number[];
}
const states = fx.states as State[];
const cases = fx.rates as Case[];
const evolution = fx.evolution as Evo[];

const toX = (x: number[]): Float64Array => {
  const X = new Float64Array(NS);
  for (let k = 0; k < N_EQ; k++) X[k] = x[k];
  return X;
};
const rel = (a: number, b: number): number => Math.abs(a - b) / Math.max(Math.abs(b), 1e-300);

describe('zeldovich rate constants', () => {
  it('GRI-3.0 set reproduces the fixture Arrhenius parameters', () => {
    fx.gri30.forEach((r, i) => {
      expect(ZELDOVICH_GRI30.A[i]).toBeCloseTo(r.A_cm3 * 1e-6, 12);
      expect(ZELDOVICH_GRI30.b[i]).toBe(r.b);
      expect(ZELDOVICH_GRI30.Ta[i] * R_UNIVERSAL).toBeCloseTo(r.Ea_cal * 4.184, 9);
    });
  });
});

describe('zeldovich vs Cantera (GRI-3.0 rates on NASA thermo)', () => {
  const kin = new ZeldovichKinetics();
  const out = new Float64Array(3);

  it('Heywood eq. 11.8 equals the quasi-steady-N mechanism rate (240 states, NO = 0 and 0.5 NO_e)', () => {
    let worst = 0;
    for (const c of cases) {
      const st = states[c.s];
      const Xeq = toX(st.Xeq);
      const w = zeldovichNORate(st.T, st.p, Xeq, c.alpha * Xeq[SP.NO]);
      worst = Math.max(worst, rel(w, c.wNOqss));
      // the Cantera state is a true N steady state
      expect(Math.abs(c.wNqss)).toBeLessThan(1e-9 * Math.abs(c.wNOqss));
    }
    expect(worst).toBeLessThan(1e-9);
  });

  it('quasiSteadyN matches the Cantera steady-state N', () => {
    let worst = 0;
    for (const c of cases) {
      const st = states[c.s];
      const X = toX(st.Xeq);
      X[SP.NO] *= c.alpha;
      worst = Math.max(worst, rel(kin.quasiSteadyN(st.T, st.p, X), c.xNqss));
    }
    expect(worst).toBeLessThan(1e-9);
  });

  it('full forward/reverse form matches Cantera net rates (actual N = N_e, NO scaled)', () => {
    let worstNO = 0;
    let worstN = 0;
    for (const c of cases) {
      const st = states[c.s];
      const X = toX(st.Xeq);
      X[SP.NO] *= c.alpha;
      let s = 0;
      for (let k = 0; k < NS; k++) s += X[k];
      for (let k = 0; k < NS; k++) X[k] /= s;
      kin.fullRates(st.T, st.p, X, out);
      worstNO = Math.max(worstNO, rel(out[0], c.wNOfull));
      worstN = Math.max(worstN, rel(out[1], c.wNfull));
    }
    expect(worstNO).toBeLessThan(1e-9);
    expect(worstN).toBeLessThan(1e-9);
  });

  it('Heywood/Hanson–Salimian set matches Cantera with the same constants', () => {
    const hk = new ZeldovichKinetics(ZELDOVICH_HEYWOOD);
    let worst = 0;
    for (const c of cases) {
      const st = states[c.s];
      const Xeq = toX(st.Xeq);
      worst = Math.max(worst, rel(hk.rateControlled(st.T, st.p, Xeq, c.alpha * Xeq[SP.NO]), c.wNOheywood));
    }
    expect(worst).toBeLessThan(2e-8); // fixture stores 9 significant digits
  });

  it('GRI-3.0-thermo reverse rates differ from NASA-thermo by ≲ 1.2 % (reported, not a defect)', () => {
    expect(fx.summary.griThermoMaxRelDiff).toBeLessThan(0.012);
    expect(fx.summary.heywoodOverGriMin).toBeGreaterThan(0.6);
    expect(fx.summary.heywoodOverGriMax).toBeLessThan(0.85);
  });
});

describe('zeldovich time evolution at constant T, p', () => {
  it('exact integrator reproduces the analytic α(t) of eq. 11.8', () => {
    const kin = new ZeldovichKinetics();
    let worst = 0;
    for (const e of evolution) {
      const Xeq = toX(e.Xeq);
      let x = 0;
      let t = 0;
      e.times.forEach((ti, j) => {
        x = kin.advanceRateControlled(e.T, e.p, Xeq, x, ti - t);
        t = ti;
        worst = Math.max(worst, Math.abs(x / Xeq[SP.NO] - e.modelAlpha[j]));
      });
      expect(kin.relaxationTime(e.T, e.p, Xeq)).toBeCloseTo(e.tauNO, 12 - Math.ceil(Math.log10(e.tauNO + 1e-30)));
      expect(rel(kin.relaxationTime(e.T, e.p, Xeq), e.tauNO)).toBeLessThan(1e-10);
    }
    expect(worst).toBeLessThan(1e-9);
  });

  it('full GRI-Mech 3.0 forms NO faster early on (N2O/NNH routes) but within a factor 2', () => {
    // Honest comparison: rate-controlled Zeldovich vs full GRI-3.0 (53 sp.) from NO-free
    // equilibrium burned gas. Early (t ≤ 0.1 τ_NO) GRI/model = 1.2–1.8; at 3 τ_NO both ≥ 0.9.
    for (const e of evolution) {
      for (let j = 0; j < 3; j++) {
        const r = e.griAlpha[j] / e.modelAlpha[j];
        expect(r).toBeGreaterThan(1.0);
        expect(r).toBeLessThan(2.0);
      }
      expect(e.griAlpha[5]).toBeGreaterThan(0.89);
      expect(e.modelAlpha[5]).toBeGreaterThan(0.92);
    }
  });

  it('advanceAlpha: step-splitting invariance, formation and decomposition vs RK4', () => {
    const rhs = (a: number, K: number): number => (1 - a * a) / (1 + K * a);
    const rk4 = (a0: number, s: number, K: number, n = 20000): number => {
      let a = a0;
      const h = s / n;
      for (let i = 0; i < n; i++) {
        const k1 = rhs(a, K);
        const k2 = rhs(a + 0.5 * h * k1, K);
        const k3 = rhs(a + 0.5 * h * k2, K);
        const k4 = rhs(a + h * k3, K);
        a += (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      }
      return a;
    };
    for (const K of [0.05, 0.5, 1, 3, 20]) {
      for (const [a0, s] of [[0, 0.3], [0.2, 2], [0.9, 5], [3, 0.2], [1.5, 1], [10, 0.05]] as const) {
        const one = advanceAlpha(a0, s, K);
        const two = advanceAlpha(advanceAlpha(a0, s / 2, K), s / 2, K);
        expect(Math.abs(one - two)).toBeLessThan(1e-12);
        expect(Math.abs(one - rk4(a0, s, K))).toBeLessThan(1e-9);
        if (a0 < 1) expect(one).toBeLessThan(1);
        else expect(one).toBeGreaterThan(1);
      }
    }
    // stiff limit: huge step goes to equilibrium without overshoot
    expect(advanceAlpha(0, 1e6, 0.7)).toBeCloseTo(1, 12);
    expect(advanceAlpha(5, 1e6, 0.7)).toBeCloseTo(1, 12);
  });

  it('returns 0 / leaves NO unchanged when equilibrium NO or N underflows', () => {
    const X = new Float64Array(NS);
    X[SP.N2] = 0.79;
    X[SP.O2] = 0.21;
    expect(zeldovichNORate(1500, 1e6, X, 1e-4)).toBe(0);
    const kin = new ZeldovichKinetics();
    expect(kin.advanceRateControlled(1500, 1e6, X, 1e-4, 1e-3)).toBe(1e-4);
  });
});

describe('zeldovich robustness (review 2026-09-30)', () => {
  const st = states[40];
  const Xeq = toX(st.Xeq);

  it('negative NO input (explicit-integrator undershoot) is treated as NO = 0', () => {
    const kin = new ZeldovichKinetics();
    const w0 = kin.rateControlled(st.T, st.p, Xeq, 0);
    expect(w0).toBeGreaterThan(0);
    // regression: for α < −1/K the old formula flipped sign (denominator 1 + αK < 0)
    const K = kin.R[0] / (kin.R[1] + kin.R[2]);
    const aBad = -2 / K;
    expect(kin.rateControlled(st.T, st.p, Xeq, aBad * Xeq[SP.NO])).toBe(w0);
    expect(kin.advanceRateControlled(st.T, st.p, Xeq, -1.5 * Xeq[SP.NO], 1e-4)).toBe(
      kin.advanceRateControlled(st.T, st.p, Xeq, 0, 1e-4),
    );
  });

  it('advanceAlpha: seeded sweep K 1e-3…1e3, α0 0…50, s 1e-12…1e8 stays finite, monotone toward 1', () => {
    let seed = 11;
    const rnd = (): number => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (let i = 0; i < 20000; i++) {
      const K = Math.exp(Math.log(1e-3) + rnd() * Math.log(1e6));
      const a0 = rnd() < 0.5 ? rnd() : 1 + Math.exp(-8 + rnd() * 12);
      const s = Math.exp(Math.log(1e-12) + rnd() * Math.log(1e20));
      const a = advanceAlpha(a0, s, K);
      expect(Number.isFinite(a)).toBe(true);
      if (a0 < 1) {
        expect(a).toBeGreaterThanOrEqual(a0 - 1e-15);
        expect(a).toBeLessThanOrEqual(1);
      } else {
        expect(a).toBeLessThanOrEqual(a0 * (1 + 1e-15));
        expect(a).toBeGreaterThanOrEqual(1);
      }
    }
  });
});

describe('zeldovich vs Cantera far outside the nominal range (review 2026-09-30)', () => {
  interface Stress { T: number; p: number; phi: number; egr: number; alpha: number; Xeq: number[]; wNOqss: number }
  const stress = (fx as unknown as { stress: Stress[] }).stress;

  it('T 1200–3400 K, p 1–250 bar, φ 0.4–2.0, 25 % EGR, α ∈ {0, 0.5, 2}: eq. 11.8 = Cantera, finite, right sign', () => {
    expect(stress.length).toBe(240);
    const kin = new ZeldovichKinetics();
    let worst = 0;
    for (const c of stress) {
      const Xeq = toX(c.Xeq);
      const w = kin.rateControlled(c.T, c.p, Xeq, c.alpha * Xeq[SP.NO]);
      expect(Number.isFinite(w)).toBe(true);
      expect(Math.sign(w)).toBe(Math.sign(c.wNOqss)); // formation for α < 1, decomposition for α > 1
      worst = Math.max(worst, rel(w, c.wNOqss));
      const x = kin.advanceRateControlled(c.T, c.p, Xeq, c.alpha * Xeq[SP.NO], 1e-3);
      expect(Number.isFinite(x)).toBe(true);
      if (c.alpha < 1) expect(x).toBeLessThanOrEqual(Xeq[SP.NO]);
      else expect(x).toBeGreaterThanOrEqual(Xeq[SP.NO]);
    }
    expect(worst).toBeLessThan(1e-9);
  });
});
