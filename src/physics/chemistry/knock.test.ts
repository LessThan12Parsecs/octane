import { describe, expect, it } from 'vitest';
import fx from '../../../test/fixtures/chemistry_knock_modes.json';
import exc from '../../../test/fixtures/chemistry_knock_excitation.json';
import { CFR_F1, CFR_KNOCK_PICKUP } from '../engines/cfr';
import {
  autoignitionBurnTime,
  besselJ,
  besselJPrimeZeros,
  cylinderModeFrequency,
  endGasBurnRate,
  endGasBurnedMass,
  EXCITATION_TIME_PRF,
  KnockOscillator,
  modeZero,
  twoZoneSoundSpeed,
} from './knock';

const zeros = fx.jnpZeros as Record<string, number[]>;

describe('Bessel functions and Draper mode constants', () => {
  it('J_n(x) matches scipy.special.jv for n = 0..8, x ≤ 60', () => {
    fx.besselJ.forEach((row, n) => {
      row.forEach((ref, i) => {
        const x = fx.besselX[i];
        expect(Math.abs(besselJ(n, x) - ref)).toBeLessThan(1e-13 + 1e-12 * Math.abs(ref));
      });
    });
  });

  it("zeros of J'_m match scipy.special.jnp_zeros (m = 0..8, first 6)", () => {
    for (const [ms, ref] of Object.entries(zeros)) {
      const z = besselJPrimeZeros(Number(ms), ref[ref.length - 1] + 0.01);
      expect(z.length).toBe(ref.length);
      z.forEach((v, i) => expect(Math.abs(v - ref[i])).toBeLessThan(1e-11));
    }
  });

  it('α_mn in Draper notation: α10 1.841, α20 3.054, α01 3.832, α30 4.201, α11 5.331', () => {
    expect(modeZero(1, 0)).toBeCloseTo(1.8412, 4);
    expect(modeZero(2, 0)).toBeCloseTo(3.0542, 4);
    expect(modeZero(0, 1)).toBeCloseTo(3.8317, 4);
    expect(modeZero(3, 0)).toBeCloseTo(4.2012, 4);
    expect(modeZero(1, 1)).toBeCloseTo(5.3314, 4);
    expect(() => modeZero(0, 0)).toThrow(RangeError);
  });

  it('CFR bore (82.55 mm) at c = 1000 m/s: first mode ≈ 7.1 kHz', () => {
    const B = CFR_F1.geometry.bore;
    expect(B).toBeCloseTo(0.08255, 6);
    const f10 = cylinderModeFrequency(1, 0, 1000, B);
    expect(f10).toBeCloseTo((1.8411837813406595 * 1000) / (Math.PI * B), 6);
    expect(f10).toBeGreaterThan(7000);
    expect(f10).toBeLessThan(7200);
    expect(cylinderModeFrequency(2, 0, 1000, B)).toBeCloseTo(11777, -1);
    expect(cylinderModeFrequency(0, 1, 1000, B)).toBeCloseTo(14775, -1);
  });
});

describe('KnockOscillator', () => {
  const P = fx.projection;

  it('default mode set is sorted by frequency with degenerate cos/sin pairs', () => {
    const ko = new KnockOscillator(0.08255);
    const names = Array.from(ko.alpha, (_, j) => `${ko.m[j]}${ko.n[j]}${ko.cosine[j] ? 'c' : 's'}`);
    expect(names.slice(0, 7)).toEqual(['10c', '10s', '20c', '20s', '01c', '30c', '30s']);
    for (let j = 1; j < ko.nModes; j++) expect(ko.alpha[j]).toBeGreaterThanOrEqual(ko.alpha[j - 1]);
  });

  it('end-gas modal projections match adaptive quadrature (scipy)', () => {
    const fine = new KnockOscillator(2 * P.R, { maxAlpha: 5.4, quadrature: [400, 800] });
    const coarse = new KnockOscillator(2 * P.R, { maxAlpha: 5.4 });
    const vf = fine.setEndGasRegion(P.x0, P.z0, P.rf);
    const vc = coarse.setEndGasRegion(P.x0, P.z0, P.rf);
    expect(Math.abs(vf - P.endGasAreaFraction)).toBeLessThan(2e-4);
    expect(Math.abs(vc - P.endGasAreaFraction)).toBeLessThan(5e-3);
    for (const md of P.modes) {
      let j = -1;
      for (let k = 0; k < fine.nModes; k++) {
        if (fine.m[k] === md.m && fine.n[k] === md.n && fine.cosine[k] === md.cosine) j = k;
      }
      expect(j).toBeGreaterThanOrEqual(0);
      expect(Math.abs(fine.sourceShape[j] - md.shape)).toBeLessThan(3e-3 * Math.max(1, Math.abs(md.shape)));
      expect(Math.abs(coarse.sourceShape[j] - md.shape)).toBeLessThan(4e-2 * Math.max(1, Math.abs(md.shape)));
    }
  });

  it('instantaneous release: modal energy → acoustic energy of the initial pressure pulse (Parseval)', () => {
    // p'_0 = (γ−1)(Q/V)(1_eg/v − 1) ⇒ E = V (γ−1)² (Q/V)² (1 − v)/v / (2ρc²)
    const B = 0.08255;
    const gamma = 1.3;
    const V = 1e-4;
    const Q = 50;
    const rho = 20;
    const c = 900;
    const ratios: number[] = [];
    for (const maxAlpha of [10, 20, 40]) {
      const ko = new KnockOscillator(B, { maxAlpha, decayTime: 1e9, quadrature: [240, 480] });
      const v = ko.setEndGasRegion(-0.012, 0.004, 0.036);
      const dt = 1e-13;
      ko.step(dt, c, gamma, V, Q / dt);
      ko.step(dt, c, gamma, V, 0);
      // every mode amplitude equals its projection of the pulse
      for (let j = 0; j < ko.nModes; j++) {
        expect(ko.eta[j]).toBeCloseTo(((gamma - 1) * Q * ko.sourceShape[j]) / V, 6);
      }
      const eExact = (V * (gamma - 1) ** 2 * (Q / V) ** 2 * (1 - v)) / v / (2 * rho * c * c);
      ratios.push(ko.acousticEnergy(rho, c, V) / eExact);
    }
    // Bessel's inequality: partial sums increase toward 1 from below
    expect(ratios[0]).toBeLessThan(ratios[1]);
    expect(ratios[1]).toBeLessThan(ratios[2]);
    expect(ratios[2]).toBeLessThan(1.0005);
    expect(ratios[2]).toBeGreaterThan(0.9);
    // the default low-order set carries most of the energy of a peripheral end-gas pulse
    expect(ratios[0]).toBeGreaterThan(0.6);
  });

  it('free oscillation: frequency α c/(πB), amplitude decays as exp(−t/τ_d)', () => {
    const B = 0.08255;
    const c = 1000;
    const tauD = 1e-3;
    const ko = new KnockOscillator(B, { maxAlpha: 1.9, decayTime: tauD, sensor: [B / 2, 0] });
    ko.setEndGasRegion(-0.03, 0, 0.03);
    ko.step(1e-12, c, 1.3, 1e-4, 1e14); // impulse
    ko.step(1e-12, c, 1.3, 1e-4, 0);
    const e0 = ko.acousticEnergy(1, c, 1e-4);
    // sample the sensor over 2 ms; interpolated zero crossings give the (damped) frequency
    const dt = 2e-6;
    let prev = ko.sensorPressure();
    const tc: number[] = [];
    for (let i = 0; i < 1000; i++) {
      ko.step(dt, c, 1.3, 1e-4, 0);
      const p = ko.sensorPressure();
      if (p * prev < 0) tc.push(i * dt + (dt * prev) / (prev - p));
      prev = p;
    }
    const fMeas = (tc.length - 1) / (2 * (tc[tc.length - 1] - tc[0]));
    const w = 2 * Math.PI * cylinderModeFrequency(1, 0, c, B);
    const fDamped = Math.sqrt(w * w - 1 / (tauD * tauD)) / (2 * Math.PI);
    expect(Math.abs(fMeas / fDamped - 1)).toBeLessThan(2e-4);
    const e2 = ko.acousticEnergy(1, c, 1e-4);
    expect(e2 / e0).toBeCloseTo(Math.exp(-2 * (2e-3 / tauD)), 2);
  });

  it('exact propagator: results independent of the step size for piecewise-constant heat release', () => {
    const run = (nPer50: number): number[] => {
      const ko = new KnockOscillator(0.08255);
      ko.setEndGasRegion(-0.01, 0.005, 0.035);
      const out: number[] = [];
      const dt = 50e-6 / nPer50;
      for (let i = 0; i < 8 * nPer50; i++) {
        const q = i < 2 * nPer50 ? 2e6 : 0; // 100 µs burn
        ko.step(dt, 950, 1.3, 1.2e-4, q);
        if ((i + 1) % nPer50 === 0) out.push(ko.sensorPressure());
      }
      return out;
    };
    const a = run(1);
    const b = run(50);
    expect(a.length).toBe(8);
    a.forEach((v, i) => expect(Math.abs(v - b[i])).toBeLessThan(1e-7 * Math.max(1, Math.abs(v))));
  });

  it('a finite release time filters the excitation: |1/(1 + iωτ_b)| for exponential release', () => {
    const B = 0.08255;
    const c = 1000;
    const V = 1e-4;
    const Q = 30;
    const amp = (tauB: number): number => {
      const ko = new KnockOscillator(B, { maxAlpha: 1.9, decayTime: 1e9 });
      ko.setEndGasRegion(-0.03, 0, 0.03);
      const dt = tauB > 0 ? tauB / 400 : 1e-13;
      if (tauB === 0) {
        ko.step(dt, c, 1.3, V, Q / dt);
        ko.step(dt, c, 1.3, V, 0);
      } else {
        for (let t = 0; t < 30 * tauB; t += dt) {
          const qMid = (Q / tauB) * Math.exp(-(t + 0.5 * dt) / tauB);
          ko.step(dt, c, 1.3, V, qMid);
        }
        ko.step(1e-9, c, 1.3, V, 0);
      }
      const w = (2 * ko.alpha[0] * c) / B;
      return Math.hypot(ko.eta[0], ko.etaDot[0] / w);
    };
    const a0 = amp(0);
    for (const tauB of [20e-6, 100e-6]) {
      const w = (2 * modeZero(1, 0) * c) / B;
      expect(amp(tauB) / a0).toBeCloseTo(1 / Math.hypot(1, w * tauB), 3);
    }
  });

  it('two-zone (Wood) sound speed reduces to the single-zone value and lies between the zones', () => {
    const p = 4e6;
    const cOf = (rho: number, g: number): number => Math.sqrt((g * p) / rho);
    expect(twoZoneSoundSpeed(p, 1e-4, 5.5, 1.25, 0, 17, 1.33)).toBeCloseTo(cOf(5.5, 1.25), 10);
    expect(twoZoneSoundSpeed(p, 0, 5.5, 1.25, 1e-4, 17, 1.33)).toBeCloseTo(cOf(17, 1.33), 10);
    const c = twoZoneSoundSpeed(p, 0.8e-4, 5.5, 1.25, 0.2e-4, 17, 1.33);
    expect(c).toBeLessThan(cOf(5.5, 1.25));
    expect(c).toBeGreaterThan(cOf(17, 1.33));
  });

  it('mode frequency tracks the bulk sound speed', () => {
    const ko = new KnockOscillator(0.08255);
    expect(ko.modeFrequency(0, 800) / ko.modeFrequency(0, 1000)).toBeCloseTo(0.8, 12);
  });
});

describe('CFR knock event (integration of the pieces)', () => {
  it('10 % end gas opposite the plug, 150 J: MAPO at the D-1 pickup falls with the burn-up time', () => {
    // Spark plug at +z, knock pickup at −z (cfr.ts); end gas = 10 % of the planform near the
    // pickup. Mean pressure rise (γ−1)Q/V = 4.1 bar. Results (bar): τ_ab = 1 µs → 32,
    // 20 µs → 13, 50 µs → 6.7, 100 µs → 3.7 — the range of light-to-heavy knock.
    const B = CFR_F1.geometry.bore;
    const sp = CFR_F1.sparkPlug.gapCenter;
    const ko = new KnockOscillator(B, { sensor: [CFR_KNOCK_PICKUP.position[0], CFR_KNOCK_PICKUP.position[2]] });
    let lo = 0;
    let hi = 2 * B;
    for (let i = 0; i < 50; i++) {
      const mid = 0.5 * (lo + hi);
      if (ko.setEndGasRegion(sp[0], sp[2], mid) > 0.1) lo = mid;
      else hi = mid;
    }
    expect(ko.setEndGasRegion(sp[0], sp[2], lo)).toBeCloseTo(0.1, 2); // midpoint-cell quantisation
    // the end gas lies on the z axis ⇒ the sin(θ) (1,0) mode dominates, cos(θ) ≈ 0
    expect(Math.abs(ko.sourceShape[1])).toBeGreaterThan(50 * Math.abs(ko.sourceShape[0]));
    const Q = 150;
    const V = 1.1e-4;
    const mapo = (tauAb: number): number => {
      ko.reset();
      let m = 0;
      const dt = 1e-6;
      for (let i = 0; i < 3000; i++) {
        ko.step(dt, 900, 1.3, V, (Q / tauAb) * Math.exp(-((i + 0.5) * dt) / tauAb));
        m = Math.max(m, Math.abs(ko.sensorPressure()));
      }
      return m;
    };
    const m = [1e-6, 20e-6, 50e-6, 100e-6].map(mapo);
    for (let i = 1; i < m.length; i++) expect(m[i]).toBeLessThan(m[i - 1]);
    expect(m[0] / 1e5).toBeGreaterThan(25);
    expect(m[0] / 1e5).toBeLessThan(40);
    expect(m[3] / 1e5).toBeGreaterThan(2);
    expect(m[3] / 1e5).toBeLessThan(6);
  });
});

describe('end-gas burn-up', () => {
  it('excitation time constant is the geometric mean of the Cantera fixture', () => {
    expect(EXCITATION_TIME_PRF).toBeCloseTo(exc.tauEGeoMean, 8);
    expect(exc.tauEMin).toBeGreaterThan(5e-7);
    expect(exc.tauEMax).toBeLessThan(2e-6);
  });

  it('burn time = τ_e + τ|∂lnτ/∂T|ΔT; rate = m_u/τ_ab', () => {
    expect(autoignitionBurnTime(1e-3, -0.005, 15, 1e-6)).toBeCloseTo(1e-6 + 75e-6, 12);
    // NTC: ∂τ/∂T = 0 → the whole end gas goes on the excitation time
    expect(autoignitionBurnTime(1e-3, 0)).toBe(EXCITATION_TIME_PRF);
    expect(endGasBurnRate(2e-5, 1e-4)).toBeCloseTo(0.2, 12);
    expect(endGasBurnRate(0, 1e-4)).toBe(0);
  });
});

describe('knock robustness (review 2026-09-30)', () => {
  /** Analytic free response of η̈ + 2dη̇ + ω²η = 0 (all damping regimes). */
  const exact = (e0: number, v0: number, w: number, d: number, t: number): number => {
    const q = w * w - d * d;
    if (q > 0) {
      const wd = Math.sqrt(q);
      return Math.exp(-d * t) * (e0 * Math.cos(wd * t) + ((v0 + d * e0) / wd) * Math.sin(wd * t));
    }
    if (q < 0) {
      const s = Math.sqrt(-q);
      return Math.exp(-d * t) * (e0 * Math.cosh(s * t) + ((v0 + d * e0) / s) * Math.sinh(s * t));
    }
    return Math.exp(-d * t) * (e0 + (v0 + d * e0) * t);
  };

  it('overdamped and critically damped modes follow the exact solution for any step', () => {
    // regression: ω_d² was floored at 1e-30, so overdamped modes (short decayTime, c → 0)
    // evolved as critically damped ones (wrong decay rate by up to ~ω/(d − √(d² − ω²)))
    const B = 0.08255;
    const c = 1000;
    const w = (2 * modeZero(1, 0) * c) / B;
    for (const ratio of [4, 1.5, 1 + 1e-12, 1]) {
      const d = ratio * w;
      for (const nSteps of [1, 7, 200]) {
        const ko = new KnockOscillator(B, { maxAlpha: 1.9, decayTime: 1 / d, sensor: [B / 2, 0] });
        ko.setEndGasRegion(-0.03, 0, 0.03);
        ko.eta[0] = 1e5;
        ko.etaDot[0] = -3e9;
        const T = 60e-6;
        for (let i = 0; i < nSteps; i++) ko.step(T / nSteps, c, 1.3, 1e-4, 0);
        const ref = exact(1e5, -3e9, w, d, T);
        expect(Math.abs(ko.eta[0] - ref)).toBeLessThan(1e-9 * 1e5);
      }
    }
    // c = 0 (no restoring force): pure decay η̇ → 0, η → η0 + η̇0/(2d)(1 − e^{−2dT})
    const ko = new KnockOscillator(B, { maxAlpha: 1.9, decayTime: 1e-3 });
    ko.eta[0] = 1;
    ko.etaDot[0] = 100;
    ko.step(1e-3, 0, 1.3, 1e-4, 0);
    expect(ko.eta[0]).toBeCloseTo(1 + (100 / 2e3) * (1 - Math.exp(-2)), 12);
    // a step far longer than every time scale decays to 0 without NaN
    ko.step(10, 0, 1.3, 1e-4, 0);
    ko.step(10, 1000, 1.3, 1e-4, 0);
    expect(Number.isFinite(ko.eta[0]) && Math.abs(ko.eta[0]) < 1e-12).toBe(true);
  });

  it('endGasBurnedMass integrates dm/dt = −m/τ exactly; stable where the explicit rate overshoots', () => {
    const m = 3e-5;
    const tau = 2e-6;
    const dt = 28e-6; // 0.1° at 600 rpm ≫ τ_ab
    expect(endGasBurnedMass(m, tau, dt)).toBeCloseTo(m * (1 - Math.exp(-dt / tau)), 18);
    expect(endGasBurnedMass(m, tau, dt)).toBeLessThanOrEqual(m);
    expect(endGasBurnRate(m, tau) * dt).toBeGreaterThan(10 * m); // explicit Euler would overshoot 14×
    // splitting invariance: two half steps = one step
    const h1 = endGasBurnedMass(m, tau, dt / 2);
    expect(h1 + endGasBurnedMass(m - h1, tau, dt / 2)).toBeCloseTo(endGasBurnedMass(m, tau, dt), 18);
    expect(endGasBurnedMass(m, 0, dt)).toBe(m);
    expect(endGasBurnedMass(0, tau, dt)).toBe(0);
    expect(endGasBurnedMass(m, tau, 0)).toBe(0);
    expect(Number.isNaN(endGasBurnedMass(m, NaN, dt))).toBe(true);
  });

  it('heat release held over one cycle step excites mode j by sinc(ω_j dt/2) (index.ts sub-step note)', () => {
    const B = 0.08255;
    const c = 1000;
    const V = 1e-4;
    const Q = 40;
    const amp = (dt: number): Float64Array => {
      const ko = new KnockOscillator(B, { decayTime: 1e9 });
      ko.setEndGasRegion(0, 0.035, 0.05);
      ko.step(dt, c, 1.3, V, Q / dt);
      ko.step(1e-12, c, 1.3, V, 0);
      const a = new Float64Array(ko.nModes);
      for (let j = 0; j < ko.nModes; j++) {
        const w = (2 * ko.alpha[j] * c) / B;
        a[j] = Math.hypot(ko.eta[j], ko.etaDot[j] / w);
      }
      return a;
    };
    const a0 = amp(1e-12);
    const a1 = amp(28e-6);
    const ko = new KnockOscillator(B);
    for (let j = 0; j < ko.nModes; j++) {
      if (a0[j] < 1e-6 * Math.max(...a0)) continue;
      const x = (ko.alpha[j] * c * 28e-6) / B; // ω dt / 2
      expect(a1[j] / a0[j]).toBeCloseTo(Math.sin(x) / x, 6);
    }
  });

  it('geometry extremes: tiny and huge bores give finite, scale-consistent projections', () => {
    for (const B of [1e-3, 0.08255, 2]) {
      const ko = new KnockOscillator(B);
      const v = ko.setEndGasRegion(0, 0.45 * B, 0.7 * B);
      expect(v).toBeGreaterThan(0);
      expect(v).toBeLessThan(1);
      for (let j = 0; j < ko.nModes; j++) expect(Number.isFinite(ko.sourceShape[j])).toBe(true);
      // projections are dimensionless: identical for geometrically similar cases
      const ref = new KnockOscillator(1);
      ref.setEndGasRegion(0, 0.45, 0.7);
      for (let j = 0; j < ko.nModes; j++) expect(ko.sourceShape[j]).toBeCloseTo(ref.sourceShape[j], 9);
    }
    // flame covering the whole bore / no flame
    const ko = new KnockOscillator(0.08255);
    expect(ko.setEndGasRegion(0, 0, 1)).toBe(0);
    for (let j = 0; j < ko.nModes; j++) expect(ko.sourceShape[j]).toBe(0);
    // uniform release excites no mode (exactly 0 analytically; midpoint quadrature leaves ≲ 1e-3)
    expect(ko.setEndGasRegion(0, 0, 0)).toBeCloseTo(1, 12);
    for (let j = 0; j < ko.nModes; j++) expect(Math.abs(ko.sourceShape[j])).toBeLessThan(1e-3);
  });
});
