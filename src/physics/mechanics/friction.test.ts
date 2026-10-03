import { describe, expect, it } from 'vitest';
import {
  chenFlynnFmep,
  FrictionTorqueModel,
  meanFrictionTorque,
  newPnhFmepBreakdown,
  oilViscosityCst,
  pnhFmep,
  type PnhFrictionInputs,
} from './friction';
import { SliderCrank } from './kinematics';

/** A round-number engine so every PNH term can be checked by hand (units of the correlation: kPa, mm, rpm). */
const ROUND: PnhFrictionInputs = {
  bore: 0.1, // B = 100 mm
  stroke: 0.1, // S = 100 mm
  cylinders: 1,
  compressionRatio: 8,
  mainBearings: { count: 2, diameter: 0.05, length: 0.02 },
  rodBearings: { count: 1, diameter: 0.04, length: 0.02 },
  camBearings: 2,
  valves: 2,
  maxValveLift: 0.01,
  valvetrain: 'OHV',
  follower: 'flat',
};

describe('Patton–Nitschke–Heywood (1989) FMEP', () => {
  it('reproduces every component term by hand at 1000 rpm, p_i = p_a (PNH/Sandoval 2002 eqs. 7a–11)', () => {
    const b = pnhFmep(ROUND, 1000, 1e5, 1e5);
    const B = 100;
    const S = 100;
    const N = 1000;
    const Sp = 2 * 0.1 * N / 60; // 3.333 m/s
    const B2S = B * B * S;
    const kPa = 1e3;
    expect(b.crankSeals).toBeCloseTo(((1.22e5 * 50) / B2S) * kPa, 8); // 6.1 kPa
    expect(b.mainBearings).toBeCloseTo(((3.03e-4 * N * 50 ** 3 * 20 * 2) / B2S) * kPa, 8);
    expect(b.turbulentDissipation).toBeCloseTo(1.35e-10 * 50 ** 2 * N * N * 2 * kPa, 8);
    expect(b.pistonSkirt).toBeCloseTo(((2.94e2 * Sp) / B) * kPa, 8);
    expect(b.rings).toBeCloseTo(4.06e4 * (1 + 1000 / N) / (B * B) * kPa, 8); // 8.12 kPa
    expect(b.rodBearings).toBeCloseTo(((3.03e-4 * N * 40 ** 3 * 20) / B2S) * kPa, 8);
    expect(b.ringGasLoading).toBeCloseTo(6.89 * (0.088 * 8 + 0.182 * Math.pow(8, 1.33 - 2.38e-2 * Sp)) * kPa, 8);
    expect(b.camBearings).toBeCloseTo(((244 * N * 2) / B2S + 4.12) * kPa, 8);
    expect(b.follower).toBeCloseTo(((400 * (1 + 1000 / N) * 2) / S) * kPa, 8); // 16 kPa
    expect(b.oscillatingHydrodynamic).toBeCloseTo(((0.5 * 10 ** 1.5 * Math.sqrt(N) * 2) / (B * S)) * kPa, 8);
    expect(b.oscillatingMixed).toBeCloseTo(((32.1 * (1 + 1000 / N) * 10 * 2) / S) * kPa, 8);
    expect(b.auxiliaries).toBeCloseTo((6.23 + 5.22e-3 * N - 1.79e-7 * N * N) * kPa, 8);
    const sum =
      b.crankSeals + b.mainBearings + b.turbulentDissipation + b.pistonSkirt + b.rings + b.rodBearings + b.ringGasLoading +
      b.camBearings + b.follower + b.oscillatingHydrodynamic + b.oscillatingMixed + b.auxiliaries;
    expect(b.total).toBeCloseTo(sum, 8);
    expect(b.constantPart + b.pistonPart).toBeCloseTo(b.total, 8);
    // sanity: an OHV engine of this size at 1000 rpm, warm: ~1 bar
    expect(b.total).toBeGreaterThan(0.6e5);
    expect(b.total).toBeLessThan(1.4e5);
  });

  it('roller followers, throttling, viscosity scaling and the 1 + 500/N option', () => {
    const base = pnhFmep(ROUND, 2000, 1e5, 1e5);
    const roller = pnhFmep({ ...ROUND, valvetrain: 'SOHC-rocker', follower: 'roller' }, 2000, 1e5, 1e5);
    expect(roller.follower).toBeCloseTo(((0.0151 * 2000 * 2) / 100) * 1e3, 8);
    // throttled: ring gas loading ∝ p_i/p_a
    const thr = pnhFmep(ROUND, 2000, 0.4e5, 1e5);
    expect(thr.ringGasLoading / base.ringGasLoading).toBeCloseTo(0.4, 12);
    // hydrodynamic terms scale with √(ν/ν₀) (Sandoval 2002 eq. 3, App. A-2), boundary terms do not
    const visc = pnhFmep({ ...ROUND, viscosityRatio: 4 }, 2000, 1e5, 1e5);
    expect(visc.mainBearings / base.mainBearings).toBeCloseTo(2, 12);
    expect(visc.pistonSkirt / base.pistonSkirt).toBeCloseTo(2, 12);
    expect(visc.rodBearings / base.rodBearings).toBeCloseTo(2, 12);
    expect(visc.oscillatingHydrodynamic / base.oscillatingHydrodynamic).toBeCloseTo(2, 12);
    expect((visc.camBearings - 4.12e3) / (base.camBearings - 4.12e3)).toBeCloseTo(2, 12);
    expect(visc.turbulentDissipation).toBe(base.turbulentDissipation);
    expect(visc.rings).toBe(base.rings);
    expect(visc.crankSeals).toBe(base.crankSeals);
    expect(visc.follower).toBe(base.follower);
    expect(visc.oscillatingMixed).toBe(base.oscillatingMixed);
    expect(visc.auxiliaries).toBe(base.auxiliaries);
    const s03 = pnhFmep({ ...ROUND, boundaryRpm: 500 }, 1000, 1e5, 1e5, newPnhFmepBreakdown());
    expect(s03.rings / pnhFmep(ROUND, 1000, 1e5, 1e5).rings).toBeCloseTo(1.5 / 2, 12);
  });

  it('Sandoval (2002) App. A-2 modified total fmep, term by term (√(μ/μ₀) scaling, 1 + 500/N)', () => {
    // Independent transcription of the rendered App. A-2 expression (kPa, mm, rpm, m/s), pumping terms omitted.
    const mu = 3.1; // μ/μ₀
    const N = 1500;
    const pi = 0.8e5;
    const pa = 1e5;
    const B = 100;
    const S = 100;
    const nc = 1;
    const Sp = (2 * 0.1 * N) / 60;
    const r = Math.sqrt(mu);
    const [Cff, , Coh, Com] = [400, 0, 0.5, 32.1]; // OHV, flat follower (Table 4.2)
    const Db = 50, Lb = 20, nb = 2, Dr = 40, Lr = 20, rc = 8, Lv = 10, nv = 2, ncam = 2;
    const K = 2.38e-2;
    const expected =
      1.22e5 * (Db / (B * B * S * nc)) +
      3.03e-4 * r * ((N * Db ** 3 * Lb * nb) / (B * B * S * nc)) +
      1.35e-10 * ((Db * Db * N * N * nb) / nc) +
      2.94e2 * r * (Sp / B) +
      4.06e4 * (1 + 500 / N) * (1 / (B * B)) +
      3.03e-4 * r * ((N * Dr ** 3 * Lr * 1) / (B * B * S * nc)) +
      6.89 * (pi / pa) * (0.088 * r * rc + 0.182 * Math.pow(rc, 1.33 - K * Sp)) + // PNH K (Sandoval's 2K not used)
      4.12 +
      244 * r * ((N * ncam) / (B * B * S * nc)) +
      Cff * (1 + 500 / N) * (nv / (S * nc)) +
      Coh * r * ((Lv ** 1.5 * Math.sqrt(N) * nv) / (B * S * nc)) +
      Com * (1 + 500 / N) * ((Lv * nv) / (S * nc)) +
      6.23 + 5.22e-3 * N - 1.79e-7 * N * N;
    const b = pnhFmep({ ...ROUND, viscosityRatio: mu, boundaryRpm: 500 }, N, pi, pa);
    expect(b.total / 1e3).toBeCloseTo(expected, 9);
  });

  it('rejects rpm ≤ 0 (boundary terms diverge) instead of returning Infinity/NaN', () => {
    expect(() => pnhFmep(ROUND, 0, 1e5, 1e5)).toThrow(RangeError);
    expect(() => pnhFmep(ROUND, -100, 1e5, 1e5)).toThrow(RangeError);
    expect(() => pnhFmep(ROUND, Number.NaN, 1e5, 1e5)).toThrow(RangeError);
  });

  it('oil viscosity (Vogel, Sandoval 2002 App. A.2): SAE 30 within the SAE J300 100 °C band', () => {
    const v100 = oilViscosityCst('SAE30', 373.15);
    expect(v100).toBeGreaterThan(9.3); // SAE J300: SAE 30 = 9.3–12.5 cSt at 100 °C
    expect(v100).toBeLessThan(12.5);
    expect(oilViscosityCst('SAE30', 330.15)).toBeGreaterThan(40); // 57 °C CFR sump: ≈ 46 cSt
    expect(oilViscosityCst('10W30', 363.15)).toBeGreaterThan(8); // ≈ PNH reference oil at 90 °C
    expect(oilViscosityCst('10W30', 363.15)).toBeLessThan(12);
  });
});

describe('Chen–Flynn FMEP', () => {
  it('evaluates A + B·p_max + C·S̄p + D·S̄p²', () => {
    const c = { A: 0.4e5, B: 0.005, C: 0.09e5, D: 0.0009e5 };
    expect(chenFlynnFmep(c, 60e5, 10)).toBeCloseTo(0.4e5 + 0.3e5 + 0.9e5 + 0.09e5, 6);
  });
});

describe('FrictionTorqueModel', () => {
  const IN = 0.0254;
  const kin = new SliderCrank({ bore: 3.25 * IN, stroke: 4.5 * IN, conRodLength: 0.254, pinOffset: 0, creviceVolume: 1e-6 }, 7);

  it('cycle-mean torque equals FMEP·V_d/(4π) for each component at the reference speed', () => {
    const f = new FrictionTorqueModel(kin);
    const w = (600 * 2 * Math.PI) / 60;
    const parts: [number, number, number][] = [
      [0.7e5, 0, 0],
      [0, 0.45e5, 0],
      [0, 0, 0.35e5],
      [0.7e5, 0.45e5, 0.35e5],
    ];
    for (const [c, co, v] of parts) {
      f.setFmep(c, co, v, w);
      const N = 4096;
      let sum = 0;
      for (let i = 0; i < N; i++) sum += f.torque(kin.dxdTheta((4 * Math.PI * i) / N), w);
      const mean = sum / N;
      const expected = -meanFrictionTorque(c + co + v, kin.displacedVolume);
      // the only deviation is the sign smoothing at ω = 62.8 rad/s vs 0.5 rad/s: 1 − ω/√(ω²+ε²) ≈ 3e-5
      expect(Math.abs(mean - expected) / Math.abs(expected)).toBeLessThan(4e-5);
    }
  });

  it('opposes rotation, viscous part ∝ ω, smooth through ω = 0', () => {
    const f = new FrictionTorqueModel(kin);
    f.setFmep(0.5e5, 0.3e5, 0.4e5, 60);
    const x1 = kin.dxdTheta(1.2);
    expect(f.torque(x1, 60)).toBeLessThan(0);
    expect(f.torque(x1, -60)).toBeGreaterThan(0);
    expect(f.torque(x1, 0) + 0).toBe(0); // (−0 → 0)
    const coul = f.constantTorque + f.coulombForce * Math.abs(x1);
    const visc = (w: number) => -f.torque(x1, w) - coul * (w / Math.sqrt(w * w + 0.25));
    expect(visc(120) / visc(60)).toBeCloseTo(2, 12);
  });

  it('setFromPnh maps skirt → viscous, rings + gas loading → Coulomb, rest → constant', () => {
    const f = new FrictionTorqueModel(kin);
    const b = pnhFmep(ROUND, 600, 1e5, 1e5);
    f.setFromPnh(b, 62.83);
    const k = kin.displacedVolume / (4 * Math.PI);
    expect(f.constantTorque).toBeCloseTo(b.constantPart * k, 12);
    expect(f.coulombForce).toBeCloseTo(((b.rings + b.ringGasLoading) * k) / (kin.pistonTravel / Math.PI), 10);
    // ⟨x'²⟩ for zero offset ≈ a²/2·(1 + a²/(4l²)+…): check against the exact series to 1e-6
    const a = kin.crankRadius;
    const l = kin.rodLength;
    expect(f.meanSqDxdTheta / (a * a / 2)).toBeCloseTo(1 + (a * a) / (4 * l * l), 3);
  });
});

describe('FrictionTorqueModel stress', () => {
  const IN = 0.0254;
  const kin = new SliderCrank({ bore: 3.25 * IN, stroke: 4.5 * IN, conRodLength: 0.254, pinOffset: 0, creviceVolume: 1e-6 }, 7);

  it('continuous in ω through 0 (no jump for an ODE integrator) and odd in ω', () => {
    const f = new FrictionTorqueModel(kin);
    f.setFmep(0.6e5, 0.4e5, 0.3e5, 62.8);
    const x1 = kin.dxdTheta(0.9);
    let prev = f.torque(x1, -1);
    for (let w = -1; w <= 1; w += 1e-3) {
      const t = f.torque(x1, w);
      expect(Math.abs(t - prev)).toBeLessThan(0.1); // N m per 1e-3 rad/s: Lipschitz, no sign jump
      expect(t).toBeCloseTo(-f.torque(x1, -w), 12);
      prev = t;
    }
  });

  it('hot path: torque() timing (printed)', () => {
    const f = new FrictionTorqueModel(kin);
    f.setFmep(0.6e5, 0.4e5, 0.3e5, 62.8);
    const n = 2_000_000;
    let acc = 0;
    const t0 = performance.now();
    for (let i = 0; i < n; i++) acc += f.torque(0.01 * Math.sin(i * 1e-3), 62.8);
    const ns = ((performance.now() - t0) * 1e6) / n;
    console.info(`FrictionTorqueModel.torque: ${ns.toFixed(1)} ns/call`);
    expect(Number.isFinite(acc)).toBe(true);
    expect(ns).toBeLessThan(1000);
  });
});
