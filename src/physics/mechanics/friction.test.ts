import { describe, expect, it } from 'vitest';
import { CFR_F1, CFR_FRICTION } from '../engines/cfr';
import { MODEL_T, MODEL_T_FORD_WOT_TABLE, MODEL_T_FRICTION } from '../engines/model-t';
import { MultiCylinderCrankTrain } from './dynamics';
import {
  chenFlynnFmep,
  FrictionTorqueModel,
  meanFrictionTorque,
  newPnhFmepBreakdown,
  oilViscosityCst,
  PNH_VALVETRAIN_CONSTANTS,
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

  it('constructor rejects a non-integer or < 1 cylinder count', () => {
    expect(() => new FrictionTorqueModel(kin, 0.5, 0)).toThrow(RangeError);
    expect(() => new FrictionTorqueModel(kin, 0.5, 2.5)).toThrow(RangeError);
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

describe('FrictionTorqueModel, multi-cylinder normalisation', () => {
  const kinT = new SliderCrank(MODEL_T.geometry, MODEL_T.geometry.compressionRatio);
  const OFF = MODEL_T.layout.firingOffsetDeg;
  const w1000 = (1000 * 2 * Math.PI) / 60;

  /** Cycle-mean (720°) friction torque of all cylinders, with the per-cylinder x′ from the crank train. */
  function cycleMean(f: FrictionTorqueModel, crank: MultiCylinderCrankTrain, omega: number): number {
    const N = 4096;
    let s = 0;
    for (let i = 0; i < N; i++) s += f.torqueCylinders(crank.update((4 * Math.PI * i) / N), omega);
    return s / N;
  }

  it('one cylinder: torqueCylinders is bit-identical to torque, setFmep to the single-cylinder formulas', () => {
    const k1 = new SliderCrank(CFR_F1.geometry, 7);
    const f = new FrictionTorqueModel(k1);
    const fN = new FrictionTorqueModel(k1, 0.5, 1);
    const b = pnhFmep(CFR_FRICTION, 600, 0.95e5, 1e5);
    const w = (600 * 2 * Math.PI) / 60;
    f.setFromPnh(b, w);
    fN.setFromPnh(b, w);
    const k = k1.displacedVolume / (4 * Math.PI);
    expect(f.constantTorque).toBe(b.constantPart * k);
    expect(f.coulombForce).toBe(((b.rings + b.ringGasLoading) * k) / f.meanAbsDxdTheta);
    expect(f.viscousCoefficient).toBe((b.pistonSkirt * k) / (Math.abs(w) * f.meanSqDxdTheta));
    expect(f.totalDisplacedVolume).toBe(k1.displacedVolume);
    const dx = new Float64Array(1);
    for (let d = -360; d < 360; d += 2.9) {
      dx[0] = k1.dxdTheta(d * Math.PI / 180);
      for (const om of [w, -w, 0.3, 0]) expect(fN.torqueCylinders(dx, om)).toBe(f.torque(dx[0], om));
    }
  });

  it('4 cylinders: cycle-mean torque = FMEP·N·V_d/(4π) for each component, for any firing offsets', () => {
    const b = pnhFmep(MODEL_T_FRICTION, 1000, 1e5, 1e5);
    const f = new FrictionTorqueModel(kinT, 0.5, 4);
    expect(f.totalDisplacedVolume).toBeCloseTo(4 * kinT.displacedVolume, 18);
    expect(f.totalDisplacedVolume).toBeCloseTo(2.8958e-3, 6); // 176.7 in³
    for (const offsets of [OFF, [0, 90, 270, 450]]) {
      const crank = new MultiCylinderCrankTrain(kinT, MODEL_T.masses, offsets);
      const parts: [number, number, number][] = [
        [b.constantPart, 0, 0],
        [0, b.rings + b.ringGasLoading, 0],
        [0, 0, b.pistonSkirt],
        [b.constantPart, b.rings + b.ringGasLoading, b.pistonSkirt],
      ];
      for (const [c, co, v] of parts) {
        f.setFmep(c, co, v, w1000);
        const expected = -meanFrictionTorque(c + co + v, f.totalDisplacedVolume);
        expect(Math.abs(cycleMean(f, crank, w1000) - expected) / Math.abs(expected)).toBeLessThan(4e-5); // sign smoothing only
      }
    }
  });

  it('the trap it fixes: a single-cylinder model fed the whole-engine FMEP gives exactly 1/4 of the 4-cylinder friction', () => {
    const b = pnhFmep(MODEL_T_FRICTION, 1000, 1e5, 1e5);
    const crank = new MultiCylinderCrankTrain(kinT, MODEL_T.masses, OFF);
    const f4 = new FrictionTorqueModel(kinT, 0.5, 4);
    const f1 = new FrictionTorqueModel(kinT); // the cycle model's construction before multi-cylinder support
    f4.setFromPnh(b, w1000);
    f1.setFromPnh(b, w1000);
    const N = 4096;
    let s1 = 0;
    for (let i = 0; i < N; i++) s1 += f1.torque(kinT.dxdTheta((4 * Math.PI * i) / N), w1000);
    const m1 = s1 / N;
    const m4 = cycleMean(f4, crank, w1000);
    console.info(`Model T friction at 1000 rpm: ${m4.toFixed(2)} N m (4-cylinder) vs ${m1.toFixed(2)} N m (single-cylinder normalisation)`);
    expect(m4 / m1).toBeCloseTo(4, 9);
    expect(-m4).toBeCloseTo(meanFrictionTorque(b.total, 4 * kinT.displacedVolume), 3);
  });

  it('torqueCylinders = −sgn_ε(ω)[T_c + F_c Σ|x′_i|] − c_v Σ x′_i² ω (explicit), smooth and odd in ω', () => {
    const f = new FrictionTorqueModel(kinT, 0.5, 4);
    f.setFmep(0.4e5, 0.15e5, 0.1e5, w1000);
    const crank = new MultiCylinderCrankTrain(kinT, MODEL_T.masses, OFF);
    for (let d = -360; d < 360; d += 13) {
      const dx = crank.update((d * Math.PI) / 180);
      let sa = 0;
      let s2 = 0;
      for (let i = 0; i < 4; i++) {
        sa += Math.abs(dx[i]);
        s2 += dx[i] * dx[i];
      }
      for (const om of [w1000, -37, 0.2]) {
        const ref = (-(f.constantTorque + f.coulombForce * sa) * om) / Math.sqrt(om * om + 0.25) - f.viscousCoefficient * s2 * om;
        expect(f.torqueCylinders(dx, om)).toBeCloseTo(ref, 10);
        expect(f.torqueCylinders(dx, -om)).toBeCloseTo(-f.torqueCylinders(dx, om), 12);
      }
    }
  });
});

describe("PNH 'L-head' valvetrain and the Model T FMEP", () => {
  it("'L-head' uses the direct-acting (SOHC-direct) constants: identical PNH breakdown", () => {
    expect(PNH_VALVETRAIN_CONSTANTS['L-head']).toEqual(PNH_VALVETRAIN_CONSTANTS['SOHC-direct']);
    for (const rpm of [400, 1000, 2000]) {
      const a = pnhFmep({ ...MODEL_T_FRICTION, valvetrain: 'L-head' }, rpm, 0.8e5, 1e5);
      const b = pnhFmep({ ...MODEL_T_FRICTION, valvetrain: 'SOHC-direct' }, rpm, 0.8e5, 1e5);
      expect(a).toEqual(b);
      // and well below the OHV train (pushrod + rocker) of the same engine
      expect(a.valvetrain).toBeLessThan(0.6 * pnhFmep({ ...MODEL_T_FRICTION, valvetrain: 'OHV' }, rpm, 0.8e5, 1e5).valvetrain);
    }
    expect(MODEL_T_FRICTION.valvetrain).toBe('L-head');
  });

  it('Model T FMEP at 400/1000/1600/2000 rpm (printed) vs the Ford WOT brake table', () => {
    // Ford WOT brake torque, lb-ft [FSB Fig. 84 via MODEL_T_FORD_WOT_TABLE]; 2000 rpm: 40 lb-ft, Ford curve as
    // tabulated by Sigworth 1999 ('simstock', row 'Stock T 1913 + .250 lift') — not in Ford's printed table.
    const LBFT = 1.3558179483314004;
    const fordTorque = (rpm: number) => (rpm === 2000 ? 40 : MODEL_T_FORD_WOT_TABLE.find((r) => r[0] === rpm)![1]) * LBFT;
    const vd = 4 * (Math.PI / 4) * MODEL_T.geometry.bore ** 2 * MODEL_T.geometry.stroke;
    // Ideal fuel-air-cycle IMEP of a full cylinder at CR 3.98, φ 1.15 (iso-octane), 330 K, 0.95 bar:
    // 12.0 bar, η 0.269 (tools/reference/cycle_fuel_air_oracle.py run_case, Cantera 3.2) — the
    // "indicated" scale: IMEP_net = (η_i/η_fa)·(m_trapped/m_ideal)·12.0 bar.
    const IMEP_FA = 12.0e5;
    const rows: string[] = [];
    for (const rpm of [400, 1000, 1600, 2000]) {
      const b = pnhFmep(MODEL_T_FRICTION, rpm, 1e5, 1e5);
      const bmep = (4 * Math.PI * fordTorque(rpm)) / vd;
      const etaM = bmep / (bmep + b.total);
      const needed = (bmep + b.total) / IMEP_FA; // (η_i/η_fa)·trapping ratio PNH implies
      rows.push(
        `${rpm} rpm: FMEP ${(b.total / 1e5).toFixed(3)} bar (piston ${(b.pistonPart / 1e5).toFixed(3)}, valvetrain ${(b.valvetrain / 1e5).toFixed(3)}), ` +
          `Ford BMEP ${(bmep / 1e5).toFixed(2)} bar → η_m ${etaM.toFixed(3)}, IMEP_net ${((bmep + b.total) / 1e5).toFixed(2)} bar = ${needed.toFixed(3)} × fuel-air`,
      );
      expect(b.total).toBeGreaterThan(0.5e5);
      expect(b.total).toBeLessThan(0.8e5);
      if (rpm <= 1600) {
        expect(etaM).toBeGreaterThan(0.8);
        expect(etaM).toBeLessThan(0.92);
      }
    }
    console.info(`Model T PNH friction (L-head, SAE 30 at 70 °C, WOT):\n  ${rows.join('\n  ')}`);
    // Period cross-check: the ALAM/SAE rating assumed η_m = 0.75 at 1000 ft/min (1500 rpm for a 4 in stroke)
    // [Good 1922 pp. 37–38]; with Ford's 20 hp at 1500 rpm that is FMEP ≈ 1.37 bar, about twice PNH.
    const bmep1500 = (4 * Math.PI * 70 * LBFT) / vd;
    const fmepAlam = bmep1500 * (1 / 0.75 - 1);
    const pnh1500 = pnhFmep(MODEL_T_FRICTION, 1500, 1e5, 1e5).total;
    console.info(`1500 rpm: ALAM η_m 0.75 → FMEP ${(fmepAlam / 1e5).toFixed(2)} bar vs PNH ${(pnh1500 / 1e5).toFixed(2)} bar`);
    expect(fmepAlam / pnh1500).toBeGreaterThan(1.5);
    expect(fmepAlam / pnh1500).toBeLessThan(2.6);
  });
});
