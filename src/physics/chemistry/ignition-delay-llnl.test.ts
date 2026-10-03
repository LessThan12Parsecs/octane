import { describe, expect, it } from 'vitest';
import eng from '../../../test/fixtures/chemistry_lw_engine.json';
import val from '../../../test/fixtures/chemistry_ignition_validation.json';
import { axisNodes, douaudEyzatTau, type TabulatedIgnitionDelay } from './ignition-delay';
import { DE_LLNL_PHI_REF, DE_LLNL_XRES_REF, douaudEyzatLLNL, prfDetailedChemistry, prfLLNLGasoline2011, prfLLNLv2 } from './ignition-delay-llnl';
import { LivengoodWuIntegrator } from './livengood-wu';

interface Off { T: number; p: number; phi: number; on: number; xres: number; tau: number; tau1: number; ignited: boolean }
interface Sweep { on: number; p: number; phi: number; xres: number; T: number[]; tau: number[]; tau1: number[] }
interface Fie { T: number; p: number; tauExp: number; tau: number; tau1: number }
interface Traj { case: number; theta: number; T: number; p: number; phi: number; on: number; tau: number; ignited: boolean }
interface MechVal { fieweger: Fie[]; offgrid: Off[]; sweep: Sweep[]; offgridLow?: Off[]; trajectory?: Traj[] }
interface EngineCase {
  on: number; phi: number; cr: number; T0: number; p0: number; rpm: number; theta0: number; dTheta: number;
  T: number[]; p: number[]; thetaIgn: number; thetaFirstStage: number | null;
}
const mechs = val.mechanisms as unknown as Record<string, MechVal>;
const engineCases = eng.cases as unknown as EngineCase[];

const quantile = (a: number[], q: number): number => {
  const s = [...a].sort((x, y) => x - y);
  return s[Math.min(s.length - 1, Math.floor(q * s.length))];
};
const geoMean = (r: number[]): number => Math.exp(r.reduce((s, x) => s + Math.log(x), 0) / r.length);

// [name, model, max |Δln τ| thresholds: median, 90 %, max; own-mechanism check at Fieweger points]
const cases: [string, TabulatedIgnitionDelay, number, number, number, number][] = [
  // measured: median 0.026, p90 0.064, max 0.100 (2.6 % / 6.6 % / 10.5 % in τ)
  ['llnl-gasoline-2011', prfLLNLGasoline2011, 0.035, 0.08, 0.12, 0.05],
  // measured: median 0.043, p90 0.148, max 0.276 (worst: rich, low T, low residual — see ignition-delay-llnl.ts)
  ['llnl-prf-v2', prfLLNLv2, 0.05, 0.16, 0.3, 0.07],
];

describe('LLNL PRF ignition-delay tables', () => {
  it('default detailed model is the 2011 LLNL gasoline-surrogate PRF chemistry', () => {
    expect(prfDetailedChemistry).toBe(prfLLNLGasoline2011);
    expect(axisNodes(prfDetailedChemistry.table.axes.invT).length).toBeGreaterThanOrEqual(17);
  });

  for (const [name, model, med, p90, max, fie] of cases) {
    const v = mechs[name];
    describe(name, () => {
      it('multilinear interpolation vs direct Cantera at 150 random off-grid points', () => {
        const err = v.offgrid.map((o) => Math.abs(model.lnTauAt(o.T, o.p, o.phi, o.on, o.xres) - Math.log(o.tau)));
        expect(quantile(err, 0.5)).toBeLessThan(med);
        expect(quantile(err, 0.9)).toBeLessThan(p90);
        expect(Math.max(...err)).toBeLessThan(max);
      });

      it('reproduces its own mechanism at the Fieweger n-heptane points inside the table', () => {
        for (const f of v.fieweger) {
          if (f.T > 1100) continue;
          expect(Math.abs(Math.log(model.tauPRF(f.T, f.p, 1, 0, 0) / f.tau))).toBeLessThan(fie);
        }
      });

      it('τ decreases with pressure at (almost) every table node (T ≥ 600 K)', () => {
        // At the 550 K nodes of the extended table τ (1–5 s) is nearly p-independent (±10 %
        // between 3.5 and 80 bar, slightly increasing for n-heptane) — excluded here.
        const t = model.table;
        const a = t.axes;
        const sP = a.phi.length * a.on.length * a.xres.length;
        const inv = axisNodes(a.invT);
        const nP = axisNodes(a.lnP).length;
        let bad = 0;
        let tot = 0;
        for (let i = 0; i < t.lnTau.length; i++) {
          if (Math.floor(i / sP) % nP === nP - 1) continue;
          if (1000 / inv[Math.floor(i / (sP * nP))] < 599) continue;
          tot++;
          if (t.lnTau[i + sP] > t.lnTau[i] + 1e-3) bad++;
        }
        expect(bad / tot).toBeLessThan(0.01);
      });
    });
  }
});

describe('default table covers the compression stroke (review 2026-09-30)', () => {
  // The end gas spends most of the compression stroke at 350–650 K and 1–10 bar, and
  // Livengood–Wu integrates from IVC. The first (legacy) table stopped at 650 K / 10 bar and
  // its slope-limited extrapolation was 2–33× too reactive there (median |Δln τ| 0.32, p90 1.64,
  // max 3.50, mean −0.59 on the points below), adding 0.05–0.4 of spurious integral before the
  // end gas even reached 650 K in the engine cases. The table now has nodes at 600 and 550 K and
  // at 5.95 and 3.54 bar.
  const low = mechs['llnl-gasoline-2011'].offgridLow!;

  it('grid spans 550–1100 K and 3.54–80 bar', () => {
    const inv = axisNodes(prfDetailedChemistry.table.axes.invT);
    const lnp = axisNodes(prfDetailedChemistry.table.axes.lnP);
    expect(1000 / inv[inv.length - 1]).toBeCloseTo(550, 6);
    expect(1000 / inv[0]).toBeCloseTo(1100, 6);
    expect(Math.exp(lnp[0]) / 1e5).toBeCloseTo(10 / Math.sqrt(8), 6);
    expect(Math.exp(lnp[lnp.length - 1]) / 1e5).toBeCloseTo(80, 6);
    // 14 entries did not ignite within t_max = 100 s: all lean/diluted (x_res = 0.15) at
    // T ≤ 650 K and p ≤ 5.95 bar, i.e. τ ≥ 100 s — irrelevant on engine time scales
    const t = prfDetailedChemistry.table;
    const nIn = t.axes.phi.length * t.axes.on.length * t.axes.xres.length;
    const capped = t.capped ?? [];
    expect(capped.length).toBeLessThanOrEqual(14);
    for (const i of capped) {
      const iT = Math.floor(i / (nIn * lnp.length));
      const iP = Math.floor(i / nIn) % lnp.length;
      expect(1000 / inv[iT]).toBeLessThanOrEqual(650 + 1e-9);
      expect(Math.exp(lnp[iP])).toBeLessThanOrEqual(6e5);
      expect(i % t.axes.xres.length).toBe(t.axes.xres.length - 1);
      expect(Math.exp(t.lnTau[i])).toBeCloseTo(t.tauMax, 2);
    }
  });

  it('table vs direct Cantera at 84 points, 550–700 K, 3–10 bar (60 random + 24-point probe)', () => {
    expect(low.length).toBe(84);
    const err = low.map((o) => prfDetailedChemistry.lnTauAt(o.T, o.p, o.phi, o.on, o.xres) - Math.log(o.tau));
    const abs = err.map(Math.abs);
    // measured: median 0.049, p90 0.116, max 0.255, mean +0.043 (legacy: 0.32 / 1.64 / 3.50 / −0.59)
    expect(quantile(abs, 0.5)).toBeLessThan(0.07);
    expect(quantile(abs, 0.9)).toBeLessThan(0.15);
    expect(Math.max(...abs)).toBeLessThan(0.3);
    // no systematic bias toward reactivity (the legacy table's mean error was −0.59)
    expect(Math.abs(err.reduce((a, b) => a + b, 0) / err.length)).toBeLessThan(0.08);
  });

  it('table = direct Cantera along the 7 engine compressions (−60…0 CAD) wherever τ < 1 s', () => {
    const traj = mechs['llnl-gasoline-2011'].trajectory!;
    expect(traj.length).toBe(42);
    // measured: |Δln τ| ≤ 0.041 for T ≥ 640 K (34 points); 0.05–0.16 at 578–624 K, inside the
    // coarse 600–650 K / 550–600 K intervals (τ 0.2–1 s there: ≲ 0.02 of LW integral)
    let worstHot = 0;
    let worstCold = 0;
    let n = 0;
    for (const o of traj) {
      const t = prfDetailedChemistry.tauPRF(o.T, o.p, o.phi, o.on, 0);
      if (!o.ignited || o.tau > 1) {
        expect(t).toBeGreaterThan(0.5); // slow states stay slow (their LW contribution ≲ 1e-2)
        continue;
      }
      n++;
      const e = Math.abs(Math.log(t / o.tau));
      if (o.T >= 640) worstHot = Math.max(worstHot, e);
      else worstCold = Math.max(worstCold, e);
    }
    expect(n).toBe(40);
    expect(worstHot).toBeLessThan(0.05);
    expect(worstCold).toBeLessThan(0.2);
  });

  it('extrapolation below 550 K keeps the low-temperature activation (τ ≥ seconds at 500 K)', () => {
    for (const on of [0, 90, 100]) {
      for (const p of [3e5, 6e5]) {
        const t550 = prfDetailedChemistry.tauPRF(550, p, 1, on, 0);
        const t500 = prfDetailedChemistry.tauPRF(500, p, 1, on, 0);
        expect(t550).toBeGreaterThan(0.5);
        // apparent activation temperature of the 600→550 K edge interval ≥ 10 000 K
        expect(Math.log(t500 / t550) / (1000 / 500 - 1000 / 550)).toBeGreaterThan(10);
      }
    }
  });
});

describe('table robustness far outside the domain (review 2026-09-30)', () => {
  it('seeded sweep T 250–2500 K, p 0.1–400 bar, φ 0–3, ON −10…120, x_res −0.1…0.6: finite, continuous', () => {
    let seed = 5;
    const rnd = (): number => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (const model of [prfLLNLGasoline2011, prfLLNLv2]) {
      for (let k = 0; k < 20000; k++) {
        const T = 250 + 2250 * rnd();
        const p = Math.exp(Math.log(1e4) + rnd() * Math.log(4e3));
        const phi = 3 * rnd();
        const on = -10 + 130 * rnd();
        const x = -0.1 + 0.7 * rnd();
        const v = model.lnTauAt(T, p, phi, on, x);
        const v1 = model.lnTauFirstStageAt(T, p, phi, on, x);
        expect(Number.isFinite(v) && Number.isFinite(v1)).toBe(true);
        // C0: a 1e-7 relative perturbation of T and p moves ln τ by < 1e-3 (no jumps; far
        // outside the table the multilinear extrapolation amplifies slopes, up to ~1e3 in ln T)
        const vp = model.lnTauAt(T * (1 + 1e-7), p * (1 + 1e-7), phi, on, x);
        expect(Math.abs(vp - v)).toBeLessThan(1e-3);
      }
      // ... and exactly at every 1000/T and ln p node edge
      for (const it of axisNodes(model.table.axes.invT)) {
        for (const lp of axisNodes(model.table.axes.lnP)) {
          const T = 1000 / it;
          const p = Math.exp(lp);
          const a = model.lnTauAt(T * (1 - 1e-9), p * (1 - 1e-9), 1, 90, 0.05);
          const b = model.lnTauAt(T * (1 + 1e-9), p * (1 + 1e-9), 1, 90, 0.05);
          expect(Math.abs(a - b)).toBeLessThan(1e-5);
        }
      }
    }
  });
});

describe('detailed chemistry vs experiment and vs Douaud–Eyzat', () => {
  it('n-heptane, 40 bar, φ = 1 vs Fieweger et al. (1997): 2011 model 1.17–1.73× slow, PRF v2 1.07–2.44×', () => {
    const r2011 = mechs['llnl-gasoline-2011'].fieweger.map((f) => f.tau / f.tauExp);
    const rV2 = mechs['llnl-prf-v2'].fieweger.map((f) => f.tau / f.tauExp);
    expect(geoMean(r2011)).toBeCloseTo(1.33, 1);
    expect(geoMean(rV2)).toBeCloseTo(1.57, 1);
    expect(Math.max(...r2011)).toBeLessThan(1.8);
    expect(Math.min(...r2011)).toBeGreaterThan(1.1);
  });

  it('NTC: n-heptane and iso-octane show a negative temperature coefficient at 40 bar; Douaud–Eyzat cannot', () => {
    const sw = mechs['llnl-gasoline-2011'].sweep;
    const at = (on: number, T: number): number => {
      const s = sw.find((x) => x.on === on)!;
      let b = 0;
      s.T.forEach((t, i) => (Math.abs(t - T) < Math.abs(s.T[b] - T) ? (b = i) : 0));
      return s.tau[b];
    };
    expect(at(0, 961)).toBeGreaterThan(at(0, 854)); // n-heptane NTC ~850–960 K
    expect(at(100, 835)).toBeGreaterThan(at(100, 768)); // iso-octane NTC ~770–835 K
    expect(douaudEyzatTau(835, 40e5, 100)).toBeLessThan(douaudEyzatTau(768, 40e5, 100));
    // the table reproduces the swept curve (its own mechanism) at 40 bar, φ = 1
    for (const s of sw) {
      s.T.forEach((T, i) => {
        expect(Math.abs(Math.log(prfDetailedChemistry.tauPRF(T, s.p, 1, s.on, 0) / s.tau[i]))).toBeLessThan(0.06);
      });
    }
  });

  it('Douaud–Eyzat is within ~2× of detailed chemistry for PRF 80–100 at 750–1000 K, 40 bar, but 3–5× short at 650 K', () => {
    const sw = mechs['llnl-gasoline-2011'].sweep.filter((s) => s.on >= 80);
    for (const s of sw) {
      s.T.forEach((T, i) => {
        const r = douaudEyzatTau(T, s.p, s.on) / s.tau[i];
        if (T >= 750 && T <= 1000) {
          expect(r).toBeGreaterThan(0.35);
          expect(r).toBeLessThan(2.0);
        }
        if (T < 655) expect(r).toBeLessThan(0.3);
      });
    }
  });
});

describe('Livengood–Wu + table vs detailed chemistry in an engine-like compression (CFR, 600 rpm)', () => {
  /** Autoignition crank angle predicted by LW on the frozen (unreacted) trajectory. */
  const predict = (c: EngineCase, tau: (T: number, p: number) => number): number => {
    const lw = new LivengoodWuIntegrator();
    lw.reset(tau(c.T[0], c.p[0]));
    const dt = ((c.dTheta / 360) * 60) / c.rpm;
    for (let i = 1; i < c.T.length; i++) {
      lw.advance(dt, tau(c.T[i], c.p[i]));
      if (lw.autoignited) return c.theta0 + (lw.ignitionTime / dt) * c.dTheta;
    }
    return NaN;
  };

  it('single-stage (PRF 90/100) ignition: LW + 2011 table within 6 CAD (mean 2) of the reacting Cantera run', () => {
    // measured errors: +5.6 (PRF 90, φ 0.5), +1.2, −0.3, +1.0 CAD. With the legacy (650 K / 10 bar)
    // table they were +2.8, +0.8, −2.4, −1.8: its too-reactive extrapolation below 650 K
    // compensated the lateness of LW on a FROZEN trajectory (the reacting run's low-temperature
    // heat release, strongest for lean PRF 90, advances ignition). The table itself is now within
    // 3 % of direct Cantera along these trajectories (next describe block), so the residual error
    // is the Livengood–Wu hypothesis, not the ignition-delay data.
    const errs: number[] = [];
    for (const c of engineCases.filter((x) => x.thetaFirstStage === null)) {
      const th = predict(c, (T, p) => prfDetailedChemistry.tauPRF(T, p, c.phi, c.on, 0));
      errs.push(th - c.thetaIgn);
      expect(Math.abs(th - c.thetaIgn)).toBeLessThan(6);
    }
    expect(errs.reduce((a, b) => a + Math.abs(b), 0) / errs.length).toBeLessThan(2.5);
  });

  it('two-stage (PRF 0/60) ignition: single-integral LW fires EARLY (known limitation, 2.4–16.5 CAD)', () => {
    // measured: −16.5 (PRF 0, φ 0.5), −8.2 (PRF 60), −2.4 (PRF 0, φ 1); legacy table: up to −20
    for (const c of engineCases.filter((x) => x.thetaFirstStage !== null)) {
      const th = predict(c, (T, p) => prfDetailedChemistry.tauPRF(T, p, c.phi, c.on, 0));
      expect(th).toBeLessThan(c.thetaIgn);
      expect(c.thetaIgn - th).toBeLessThan(17);
    }
  });

  it('integral accumulated before the end gas reaches 650 K is real chemistry, not extrapolation', () => {
    // measured 0.008–0.094 (largest: n-heptane, whose first stage is active at 600–650 K);
    // the legacy table's extrapolation gave 0.049–0.407 on the same trajectories
    for (const c of engineCases) {
      const lw = new LivengoodWuIntegrator();
      lw.reset(prfDetailedChemistry.tauPRF(c.T[0], c.p[0], c.phi, c.on, 0));
      const dt = ((c.dTheta / 360) * 60) / c.rpm;
      for (let i = 1; i < c.T.length && c.T[i] < 650; i++) {
        lw.advance(dt, prfDetailedChemistry.tauPRF(c.T[i], c.p[i], c.phi, c.on, 0));
      }
      expect(lw.integral).toBeLessThan(0.1);
    }
  });

  it('Douaud–Eyzat + LW fires 2–10 CAD earlier than detailed chemistry in these HCCI-like cases', () => {
    for (const c of engineCases.filter((x) => x.on >= 90)) {
      const th = predict(c, (T, p) => douaudEyzatTau(T, p, c.on));
      expect(th).toBeLessThan(c.thetaIgn - 1.5);
      expect(th).toBeGreaterThan(c.thetaIgn - 13);
    }
  });
});

describe('douaudEyzatLLNL: Douaud–Eyzat with LLNL relative φ / residual / low-ON sensitivities (fixer round 2)', () => {
  const prf = (on: number) => ({ label: `PRF${on}`, X: new Float64Array(17), octaneNumber: on });
  it('is Douaud–Eyzat at the reference mixture (φ 1.1, x_res 0.06) for ON ≥ 80', () => {
    for (const [T, p] of [[700, 20e5], [880, 41e5], [1005, 36e5], [450, 2e5]]) {
      for (const on of [80, 90, 100]) {
        expect(douaudEyzatLLNL.tau(T, p, DE_LLNL_PHI_REF, prf(on), DE_LLNL_XRES_REF) / douaudEyzatTau(T, p, on) - 1).toBeCloseTo(0, 12);
      }
    }
  });
  it('carries the detailed-chemistry φ and residual trends (leaner / more dilute → longer delay)', () => {
    const t = (phi: number, x: number): number => douaudEyzatLLNL.tau(880, 41e5, phi, prf(95), x);
    expect(t(1.0, 0.06) / t(1.1, 0.06)).toBeGreaterThan(1.05);
    expect(t(0.9, 0.06)).toBeGreaterThan(t(1.0, 0.06));
    expect(t(1.1, 0.12)).toBeGreaterThan(t(1.1, 0.06));
  });
  it('has a finite, physical delay for n-heptane and low-ON PRFs (Douaud–Eyzat alone: τ = 0 at ON 0)', () => {
    expect(douaudEyzatTau(400, 1.2e5, 0)).toBe(0);
    for (const on of [0, 10, 40, 60]) {
      const cold = douaudEyzatLLNL.tau(400, 1.2e5, 1.1, prf(on), 0.06);
      const hot = douaudEyzatLLNL.tau(880, 41e5, 1.1, prf(on), 0.06);
      expect(cold).toBeGreaterThan(1); // no autoignition in milliseconds at 400 K
      expect(hot).toBeGreaterThan(0.2 * prfLLNLGasoline2011.tauPRF(880, 41e5, 1.1, on, 0.06));
      expect(hot).toBeLessThan(douaudEyzatLLNL.tau(880, 41e5, 1.1, prf(80), 0.06));
    }
    // monotone in ON
    let prev = 0;
    for (const on of [0, 20, 40, 60, 80, 90, 100]) {
      const v = douaudEyzatLLNL.tau(880, 41e5, 1.1, prf(on), 0.06);
      expect(v).toBeGreaterThan(prev);
      prev = v;
    }
  });
});
