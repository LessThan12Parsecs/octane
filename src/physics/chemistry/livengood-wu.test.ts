import { describe, expect, it } from 'vitest';
import { LivengoodWuIntegrator, TwoStageLivengoodWu } from './livengood-wu';

describe('LivengoodWuIntegrator', () => {
  it('constant τ: I = t/τ and ignition at t = τ (both quadratures)', () => {
    for (const q of ['logarithmic', 'trapezoidal'] as const) {
      const lw = new LivengoodWuIntegrator(q);
      lw.reset();
      const tau = 3.3e-3;
      for (let i = 0; i < 100; i++) lw.advance(1e-4, tau);
      expect(lw.integral).toBeCloseTo(1e-2 / tau, 12);
      expect(lw.ignitionTime).toBeCloseTo(tau, 14);
      expect(lw.autoignited).toBe(true);
    }
  });

  it('exponentially shrinking τ: logarithmic quadrature is exact for any step, trapezoid O(dt²)', () => {
    // τ(t) = τ0 e^{−t/θ} ⇒ I(t) = θ (e^{t/θ} − 1)/τ0, ignition at t* = θ ln(1 + τ0/θ)
    const tau0 = 50e-3;
    const theta = 2e-3;
    const tStar = theta * Math.log(1 + tau0 / theta);
    const run = (q: 'logarithmic' | 'trapezoidal', dt: number): [number, number] => {
      const lw = new LivengoodWuIntegrator(q);
      lw.reset(tau0);
      let t = 0;
      while (t < 8e-3 - 1e-15) {
        t += dt;
        lw.advance(dt, tau0 * Math.exp(-t / theta));
      }
      return [lw.integral, lw.ignitionTime];
    };
    const exact = (theta * (Math.exp(8e-3 / theta) - 1)) / tau0;
    for (const dt of [1e-3, 2e-4, 2e-5]) {
      const [I, ti] = run('logarithmic', dt);
      expect(I / exact - 1).toBeCloseTo(0, 11);
      expect(ti).toBeCloseTo(tStar, 12);
    }
    const e1 = Math.abs(run('trapezoidal', 4e-4)[0] / exact - 1);
    const e2 = Math.abs(run('trapezoidal', 2e-4)[0] / exact - 1);
    expect(e1 / e2).toBeCloseTo(4, 0); // second order
    expect(e1).toBeGreaterThan(1e-3);
  });

  it('Arrhenius τ along a compression ramp: converges, logarithmic beats trapezoid at coarse steps', () => {
    // τ = 17.68e-3 · 0.9^3.402 · (p/atm)^−1.7 · exp(3800/T), T and p rising in time
    const tau = (t: number): number => {
      const T = 700 + 2.5e4 * t; // K
      const p = 20 + 2e3 * t; // atm
      return 17.68e-3 * Math.pow(0.9, 3.402) * Math.pow(p, -1.7) * Math.exp(3800 / T);
    };
    const run = (q: 'logarithmic' | 'trapezoidal', n: number): number => {
      const lw = new LivengoodWuIntegrator(q);
      lw.reset(tau(0));
      const dt = 10e-3 / n;
      for (let i = 1; i <= n; i++) lw.advance(dt, tau(i * dt));
      return lw.integral;
    };
    const ref = run('logarithmic', 200000);
    const eLog = Math.abs(run('logarithmic', 20) / ref - 1);
    const eTrap = Math.abs(run('trapezoidal', 20) / ref - 1);
    // 20 steps over 10 ms (0.5 ms ≈ 1.8° at 600 rpm): log-mean error 2.9e-4, trapezoid 8.5e-4
    expect(eLog).toBeLessThan(5e-4);
    expect(eTrap).toBeGreaterThan(2 * eLog);
  });

  it('τ = ∞ contributes nothing; τ ≤ 0 is immediate autoignition', () => {
    const lw = new LivengoodWuIntegrator();
    lw.reset();
    lw.advance(1e-3, Infinity);
    lw.advance(1e-3, Infinity);
    expect(lw.integral).toBe(0);
    lw.advance(1e-4, 0);
    expect(lw.autoignited).toBe(true);
    expect(lw.ignitionTime).toBeGreaterThanOrEqual(2e-3);
    expect(lw.ignitionTime).toBeLessThanOrEqual(2.1e-3);
  });

  it('reset restarts the integral and clears the ignition time', () => {
    const lw = new LivengoodWuIntegrator();
    lw.reset();
    lw.advance(2e-3, 1e-3);
    expect(lw.autoignited).toBe(true);
    lw.reset();
    expect(lw.integral).toBe(0);
    expect(Number.isNaN(lw.ignitionTime)).toBe(true);
    expect(lw.time).toBe(0);
  });
});

describe('TwoStageLivengoodWu', () => {
  it('constant delays: first stage at τ1, hot ignition at τ', () => {
    const lw = new TwoStageLivengoodWu();
    lw.reset();
    const tau = 4e-3;
    const tau1 = 2.5e-3;
    let t = 0;
    let tFirst = NaN;
    let tHot = NaN;
    const dt = 1e-5;
    while (t < 6e-3) {
      lw.advance(dt, tau, tau1);
      t += dt;
      if (Number.isNaN(tFirst) && lw.first.autoignited) tFirst = lw.first.ignitionTime;
      if (Number.isNaN(tHot) && lw.autoignited) tHot = tFirst + lw.second.ignitionTime;
    }
    expect(tFirst).toBeCloseTo(tau1, 10);
    expect(tHot).toBeCloseTo(tau, 9);
    expect(lw.progress).toBe(2);
  });
});

describe('LivengoodWuIntegrator robustness (review 2026-09-30)', () => {
  it('NaN τ (e.g. PRF table asked about a non-PRF fuel) gives a NaN integral, never autoignition', () => {
    // regression: rate(NaN) used to fall into the "τ ≤ 0 ⇒ ∞" branch → knock at the first step
    const lw = new LivengoodWuIntegrator();
    lw.reset();
    lw.advance(1e-5, 1e-2);
    lw.advance(1e-5, NaN);
    expect(Number.isNaN(lw.integral)).toBe(true);
    expect(lw.autoignited).toBe(false);
    expect(Number.isNaN(lw.ignitionTime)).toBe(true);
    lw.advance(1e-5, 1e-6); // stays NaN until reset
    expect(lw.autoignited).toBe(false);
    lw.reset();
    lw.advance(1e-3, 1e-3);
    expect(lw.autoignited).toBe(true);
    const two = new TwoStageLivengoodWu();
    two.reset();
    two.advance(1e-5, NaN, NaN);
    expect(two.autoignited).toBe(false);
    expect(Number.isNaN(two.progress)).toBe(true);
  });

  it('seeded random sweep: integral finite, non-decreasing, crossing time inside the crossing step', () => {
    let seed = 7;
    const rnd = (): number => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (const q of ['logarithmic', 'trapezoidal'] as const) {
      for (let run = 0; run < 200; run++) {
        const lw = new LivengoodWuIntegrator(q);
        lw.reset(Math.exp(-12 + 14 * rnd()));
        let prev = 0;
        for (let i = 0; i < 400; i++) {
          const dt = Math.exp(-16 + 12 * rnd()); // 0.1 µs … 20 ms
          const tau = rnd() < 0.02 ? Infinity : Math.exp(-14 + 18 * rnd()); // 1e-6 s … 60 s, jumps
          const t0 = lw.time;
          const wasIgn = lw.autoignited;
          const I = lw.advance(dt, tau);
          expect(Number.isFinite(I)).toBe(true);
          expect(I).toBeGreaterThanOrEqual(prev);
          prev = I;
          if (!wasIgn && lw.autoignited) {
            expect(lw.ignitionTime).toBeGreaterThanOrEqual(t0);
            expect(lw.ignitionTime).toBeLessThanOrEqual(t0 + dt * (1 + 1e-12));
          }
        }
      }
    }
  });
});
