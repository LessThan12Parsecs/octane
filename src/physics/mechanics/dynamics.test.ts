import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/mechanics_dynamics.json';
import { DEG } from '../core/constants';
import { CFR_F1 } from '../engines/cfr';
import { MODEL_T } from '../engines/model-t';
import {
  CrankTrainDynamics,
  MultiCylinderCrankTrain,
  newCrankTrainState,
  newMultiCylinderCrankTrainState,
  type CrankTrainMasses,
} from './dynamics';
import { newKinematicState, SliderCrank, type SliderCrankGeometry } from './kinematics';

interface DynCase {
  name: string;
  geometry: SliderCrankGeometry;
  masses: CrankTrainMasses;
  gravity: number;
  states: { omega: number; alpha: number }[];
  rows: { thetaDeg: number; J: number; dJ: number; gravityTorque: number; inertiaTorque: number[] }[];
}
const cases = (fixture as { cases: DynCase[] }).cases;

interface MultiCase {
  name: string;
  geometry: SliderCrankGeometry;
  masses: CrankTrainMasses;
  gravity: number;
  firingOffsetDeg: number[];
  pCrankcase: number;
  states: { omega: number; alpha: number; externalTorque: number; loadInertia: number; pressures: number[] }[];
  rows: {
    thetaDeg: number;
    J: number;
    dJ: number;
    gravityTorque: number;
    dxdTheta: number[];
    gasTorque: number[][];
    inertiaTorque: number[];
    alpha: number[];
  }[];
}
const multiCases = (fixture as unknown as { multiCases: MultiCase[] }).multiCases;

const IN = 0.0254;
const CFR_GEOM: SliderCrankGeometry = { bore: 3.25 * IN, stroke: 4.5 * IN, conRodLength: 0.254, pinOffset: 0, creviceVolume: 1e-6 };
const MASSES: CrankTrainMasses = { piston: 1.6, conRod: 1.9, conRodCgFromBigEnd: 0.075, conRodInertiaCg: 0.012, rotatingInertia: 2.0 };

/** A smooth, strongly non-trivial synthetic 720° pressure trace (Pa). */
function syntheticPressure(sc: SliderCrank, theta: number): number {
  // wrap to [-π·2, 2π)
  let t = theta;
  while (t < -2 * Math.PI) t += 4 * Math.PI;
  while (t >= 2 * Math.PI) t -= 4 * Math.PI;
  const vBdc = sc.clearanceVolume + sc.displacedVolume;
  const poly = 0.95e5 * Math.pow(vBdc / sc.volume(t), 1.32);
  const burn = 1 + 2.2 * Math.exp(-(((t - 15 * DEG) / (25 * DEG)) ** 2));
  const exchange = t < -Math.PI || t > Math.PI ? 0.35 : 1; // lower pressure during gas exchange strokes
  return poly * burn * exchange + 0.2e5 * Math.sin(3 * t);
}

describe('CrankTrainDynamics vs numerical-Lagrangian oracle', () => {
  for (const c of cases) {
    it(`${c.name}: J_m, dJ_m/dθ, inertia and gravity torques`, () => {
      const kin = new SliderCrank(c.geometry, 9);
      const dyn = new CrankTrainDynamics(kin, c.masses, { gravity: c.gravity });
      const out = newCrankTrainState();
      let eJ = 0;
      let eT = 0;
      let eG = 0;
      let Jscale = 0;
      let Tscale = 0;
      for (const r of c.rows) Jscale = Math.max(Jscale, Math.abs(r.J));
      for (const r of c.rows) for (const t of r.inertiaTorque) Tscale = Math.max(Tscale, Math.abs(t));
      for (const r of c.rows) {
        const th = r.thetaDeg * DEG;
        c.states.forEach((st, i) => {
          dyn.evaluate(th, st.omega, st.alpha, 1e5, 1e5, out);
          eJ = Math.max(eJ, Math.abs(out.mechanismInertia - r.J) / Jscale, Math.abs(out.mechanismInertiaDerivative - r.dJ) / Jscale);
          eT = Math.max(eT, Math.abs(out.inertiaTorque - r.inertiaTorque[i]) / Tscale);
          eG = Math.max(eG, Math.abs(out.gravityTorque - r.gravityTorque));
        });
      }
      expect(eJ).toBeLessThan(1e-8);
      expect(eT).toBeLessThan(1e-8);
      expect(eG).toBeLessThan(1e-8); // N m
    });
  }
});

describe('CrankTrainDynamics identities', () => {
  it('Newton–Euler crank torque = gas + inertia + gravity torque (independent routes)', () => {
    for (const pinOffset of [0, 0.004, -0.006]) {
      const kin = new SliderCrank({ ...CFR_GEOM, pinOffset }, 8);
      const dyn = new CrankTrainDynamics(kin, MASSES);
      const out = newCrankTrainState();
      for (let d = -360; d < 360; d += 7) {
        const th = d * DEG;
        const p = syntheticPressure(kin, th);
        dyn.evaluate(th, 94.2, -120, p, 1.01e5, out);
        const sum = out.gasTorque + out.inertiaTorque + out.gravityTorque;
        expect(Math.abs(out.crankTorque - sum)).toBeLessThan(1e-9 * (1 + Math.abs(out.gasTorque)));
      }
    }
  });

  it('energy check: ∮ T_gas dθ (Newton–Euler force path) = ∮ (p − p_cc) dV over the 720° cycle', () => {
    const kin = new SliderCrank(CFR_GEOM, 7);
    const massless = new CrankTrainDynamics(kin, { piston: 0, conRod: 0, conRodCgFromBigEnd: 0.1, conRodInertiaCg: 0, rotatingInertia: 1 }, { gravity: 0 });
    const out = newCrankTrainState();
    const N = 28800; // 0.025° steps
    const h = (4 * Math.PI) / N;
    const pcc = 1.013e5;
    let workTorque = 0; // Simpson on torque
    let workPdV = 0; // trapezoid-in-V with midpoint pressure: Σ p(θ_{i+½}) (V_{i+1} − V_i)
    for (let i = 0; i <= N; i++) {
      const th = -2 * Math.PI + i * h;
      massless.evaluate(th, 60, 0, syntheticPressure(kin, th), pcc, out);
      const w = i === 0 || i === N ? 1 : i % 2 === 1 ? 4 : 2;
      workTorque += (w * out.crankTorque * h) / 3;
      if (i < N) {
        const tm = th + h / 2;
        workPdV += (syntheticPressure(kin, tm) - pcc) * (kin.volume(th + h) - kin.volume(th));
      }
    }
    expect(Math.abs(workTorque)).toBeGreaterThan(100); // non-trivial net work (J)
    expect(Math.abs(workTorque - workPdV) / Math.abs(workTorque)).toBeLessThan(1e-7);
  });

  it('inertia and gravity torques integrate to zero over a revolution at constant ω', () => {
    for (const pinOffset of [0, 0.005]) {
      const kin = new SliderCrank({ ...CFR_GEOM, pinOffset }, 7);
      const dyn = new CrankTrainDynamics(kin, MASSES);
      const out = newCrankTrainState();
      const N = 3600;
      const h = (2 * Math.PI) / N;
      let intIn = 0;
      let intIn2 = 0;
      let intG = 0;
      let peak = 0;
      for (let i = 0; i < N; i++) {
        dyn.evaluate(0.3 + i * h, 2 * Math.PI * 15, 0, 1e5, 1e5, out); // periodic → rectangle rule is spectrally accurate
        intIn += out.inertiaTorque * h;
        intIn2 += out.inertiaTorqueTwoMass * h;
        intG += out.gravityTorque * h;
        peak = Math.max(peak, Math.abs(out.inertiaTorque));
      }
      expect(peak).toBeGreaterThan(20); // ≈ 41 N m at 900 rpm for these masses
      expect(Math.abs(intIn)).toBeLessThan(1e-10 * peak);
      expect(Math.abs(intIn2)).toBeLessThan(1e-10 * peak);
      expect(Math.abs(intG)).toBeLessThan(1e-12);
    }
  });

  it('two-mass approximation: differs from the exact rigid rod only by the rod inertia couple', () => {
    const kin = new SliderCrank(CFR_GEOM, 7);
    const dyn = new CrankTrainDynamics(kin, MASSES);
    const out = newCrankTrainState();
    const dI = dyn.rodInertiaDefect;
    const l = kin.rodLength;
    expect(dI).toBeCloseTo(MASSES.conRodInertiaCg - MASSES.conRod * MASSES.conRodCgFromBigEnd * (l - MASSES.conRodCgFromBigEnd), 15);
    let maxRel = 0;
    let peak = 0;
    for (let d = -180; d < 180; d += 3) {
      const th = d * DEG;
      const w = 2 * Math.PI * 15;
      const al = 200;
      dyn.evaluate(th, w, al, 1e5, 1e5, out);
      const b1 = kin.dRodAngledTheta(th);
      const b2 = kin.d2RodAngledTheta2(th);
      const couple = -dI * (b1 * b1 * al + b1 * b2 * w * w);
      expect(out.inertiaTorque - out.inertiaTorqueTwoMass).toBeCloseTo(couple, 10);
      peak = Math.max(peak, Math.abs(out.inertiaTorque));
      maxRel = Math.max(maxRel, Math.abs(out.inertiaTorque - out.inertiaTorqueTwoMass));
    }
    // expected order: |ΔI|(a/l)²ω² relative to m_rec a² ω²  → a few % for this rod
    const a = kin.crankRadius;
    const order = (Math.abs(dI) * (a / l) ** 2) / ((MASSES.piston + dyn.rodReciprocatingMass) * a * a);
    expect(maxRel / peak).toBeLessThan(3 * order);
    expect(maxRel / peak).toBeGreaterThan(0.1 * order);
    // and they coincide when the rod is dynamically equivalent to the two masses
    const lg = MASSES.conRodCgFromBigEnd;
    const eq = new CrankTrainDynamics(kin, { ...MASSES, conRodInertiaCg: MASSES.conRod * lg * (l - lg) });
    for (let d = -180; d < 180; d += 9) {
      eq.evaluate(d * DEG, 150, 300, 1e5, 1e5, out);
      expect(out.inertiaTorque).toBeCloseTo(out.inertiaTorqueTwoMass, 10);
    }
  });

  it('massless mechanism: rod force = F/cosβ, side thrust = F·tanβ (quasi-static textbook result)', () => {
    const kin = new SliderCrank({ ...CFR_GEOM, pinOffset: 0.003 }, 8);
    const dyn = new CrankTrainDynamics(kin, { piston: 0, conRod: 0, conRodCgFromBigEnd: 0.1, conRodInertiaCg: 0, rotatingInertia: 1 }, { gravity: 0 });
    const out = newCrankTrainState();
    for (let d = -170; d < 180; d += 10) {
      const th = d * DEG;
      dyn.evaluate(th, 100, 0, 40e5, 1e5, out);
      const b = kin.rodAngle(th);
      expect(out.rodAxialForce).toBeCloseTo(out.gasForce / Math.cos(b), 8);
      expect(out.sideThrust).toBeCloseTo(out.gasForce * Math.tan(b), 8);
      // expansion stroke (crank pin at +x): piston bears on the major-thrust (−x) wall → sideThrust > 0
      if (d > 5 && d < 175) expect(out.sideThrust).toBeGreaterThan(0);
    }
  });

  it('free coasting (no gas, no friction) conserves kinetic + potential energy (RK4)', () => {
    const kin = new SliderCrank(CFR_GEOM, 7);
    const dyn = new CrankTrainDynamics(kin, { ...MASSES, rotatingInertia: 0.05 }); // light flywheel → strong speed ripple
    let th = 0.1;
    let w = 2 * Math.PI * 5;
    const p = 1e5;
    const E0 = dyn.kineticEnergy(th, w) + dyn.potentialEnergy(th);
    const dt = 2e-5;
    let wMin = Infinity;
    let wMax = -Infinity;
    for (let i = 0; i < 20000; i++) {
      const k1t = w;
      const k1w = dyn.angularAcceleration(th, w, p, p, 0);
      const k2t = w + 0.5 * dt * k1w;
      const k2w = dyn.angularAcceleration(th + 0.5 * dt * k1t, k2t, p, p, 0);
      const k3t = w + 0.5 * dt * k2w;
      const k3w = dyn.angularAcceleration(th + 0.5 * dt * k2t, k3t, p, p, 0);
      const k4t = w + dt * k3w;
      const k4w = dyn.angularAcceleration(th + dt * k3t, k4t, p, p, 0);
      th += (dt / 6) * (k1t + 2 * k2t + 2 * k3t + k4t);
      w += (dt / 6) * (k1w + 2 * k2w + 2 * k3w + k4w);
      wMin = Math.min(wMin, w);
      wMax = Math.max(wMax, w);
    }
    const E1 = dyn.kineticEnergy(th, w) + dyn.potentialEnergy(th);
    expect((wMax - wMin) / w).toBeGreaterThan(0.05); // the mechanism inertia really modulates the speed
    expect(Math.abs(E1 - E0) / E0).toBeLessThan(1e-9);
  });

  it('holdingTorque: applying it as an external torque yields the requested acceleration', () => {
    const kin = new SliderCrank(CFR_GEOM, 7);
    const dyn = new CrankTrainDynamics(kin, MASSES);
    for (let d = -360; d < 360; d += 30) {
      const th = d * DEG;
      const p = syntheticPressure(kin, th);
      const ext = -12; // e.g. friction
      const Th = dyn.holdingTorque(th, 62.8, 0, p, 1e5, ext);
      expect(dyn.angularAcceleration(th, 62.8, p, 1e5, ext + Th)).toBeCloseTo(0, 9);
      const Th2 = dyn.holdingTorque(th, 62.8, 37, p, 1e5, ext);
      expect(dyn.angularAcceleration(th, 62.8, p, 1e5, ext + Th2)).toBeCloseTo(37, 9);
    }
  });
});

describe('CrankTrainDynamics stress and performance', () => {
  let seed = 4242;
  const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;

  it('random mechanisms/states: finite outputs, Newton–Euler = Lagrange, J_m > 0', () => {
    const out = newCrankTrainState();
    for (let n = 0; n < 80; n++) {
      const stroke = 0.03 + 0.3 * rnd();
      const a = stroke / 2;
      const l = a * (1.4 + 5 * rnd());
      const geom = { bore: 0.03 + 0.2 * rnd(), stroke, conRodLength: l, pinOffset: (2 * rnd() - 1) * 0.9 * (l - a), creviceVolume: 0 };
      const kin = new SliderCrank(geom, 2 + 20 * rnd());
      const mr = 0.1 + 3 * rnd();
      const lg = l * (0.1 + 0.5 * rnd());
      const masses = { piston: 0.1 + 3 * rnd(), conRod: mr, conRodCgFromBigEnd: lg, conRodInertiaCg: mr * lg * (l - lg) * (0.5 + 0.7 * rnd()), rotatingInertia: 0.01 + 2 * rnd() };
      const dyn = new CrankTrainDynamics(kin, masses);
      for (let k = 0; k < 20; k++) {
        const th = (2 * rnd() - 1) * 4 * Math.PI;
        const w = (2 * rnd() - 1) * 1000; // up to ±9500 rpm, both directions
        const al = (2 * rnd() - 1) * 1e5;
        const p = 1e3 + 300e5 * rnd();
        dyn.evaluate(th, w, al, p, 1e5, out);
        for (const v of Object.values(out)) expect(Number.isFinite(v)).toBe(true);
        expect(out.mechanismInertia).toBeGreaterThan(0);
        const sum = out.gasTorque + out.inertiaTorque + out.gravityTorque;
        expect(Math.abs(out.crankTorque - sum)).toBeLessThan(1e-9 * (Math.abs(out.gasTorque) + Math.abs(out.inertiaTorque) + 1));
        const acc = dyn.angularAcceleration(th, w, p, 1e5, 0);
        expect(Number.isFinite(acc)).toBe(true);
        // at the free acceleration the holding torque vanishes
        expect(dyn.holdingTorque(th, w, acc, p, 1e5, 0)).toBeCloseTo(0, 6);
      }
    }
  });

  it('hot path timings (printed): evaluate, angularAcceleration', () => {
    const kin = new SliderCrank(CFR_GEOM, 7);
    const dyn = new CrankTrainDynamics(kin, MASSES);
    const out = newCrankTrainState();
    const n = 1_000_000;
    let acc = 0;
    for (let i = 0; i < 20000; i++) acc += dyn.evaluate(i * 1e-3, 62.8, 0, 2e6, 1e5, out).crankTorque;
    let t0 = performance.now();
    for (let i = 0; i < n; i++) acc += dyn.evaluate(i * 3.1e-6, 62.8, 3, 2e6, 1e5, out).inertiaTorque;
    const nsEval = ((performance.now() - t0) * 1e6) / n;
    t0 = performance.now();
    for (let i = 0; i < n; i++) acc += dyn.angularAcceleration(i * 3.1e-6, 62.8, 2e6, 1e5, -5);
    const nsAcc = ((performance.now() - t0) * 1e6) / n;
    console.info(`CrankTrainDynamics.evaluate: ${nsEval.toFixed(1)} ns/call, angularAcceleration: ${nsAcc.toFixed(1)} ns/call`);
    expect(Number.isFinite(acc)).toBe(true);
    expect(nsEval).toBeLessThan(2000);
    expect(nsAcc).toBeLessThan(2000);
  });
});

/** 720°-periodic synthetic pressure of cylinder i (its local angle θ − offset_i), Pa. */
function cylinderPressure(sc: SliderCrank, theta: number, offsetDeg: number): number {
  return syntheticPressure(sc, theta - offsetDeg * DEG);
}

describe('MultiCylinderCrankTrain vs multi-cylinder numerical-Lagrangian oracle', () => {
  for (const c of multiCases) {
    it(`${c.name}: Σ J_m, Σ dJ_m/dθ, gravity, inertia, per-cylinder x′ and gas torque, free α`, () => {
      const kin = new SliderCrank(c.geometry, 5);
      const crank = new MultiCylinderCrankTrain(kin, c.masses, c.firingOffsetDeg, { gravity: c.gravity });
      const out = newMultiCylinderCrankTrainState(crank.cylinders);
      const N = crank.cylinders;
      let Jscale = 0;
      let Tscale = 0;
      let Gscale = 0;
      let Ascale = 0;
      for (const r of c.rows) {
        Jscale = Math.max(Jscale, Math.abs(r.J));
        for (const t of r.inertiaTorque) Tscale = Math.max(Tscale, Math.abs(t));
        for (const g of r.gasTorque) for (const t of g) Gscale = Math.max(Gscale, Math.abs(t));
        for (const a of r.alpha) Ascale = Math.max(Ascale, Math.abs(a));
      }
      let eJ = 0;
      let eT = 0;
      let eG = 0;
      let eGas = 0;
      let eX = 0;
      let eA = 0;
      for (const r of c.rows) {
        const th = r.thetaDeg * DEG;
        const dx = crank.update(th);
        for (let i = 0; i < N; i++) eX = Math.max(eX, Math.abs(dx[i] - r.dxdTheta[i]) / (c.geometry.stroke / 2));
        eJ = Math.max(eJ, Math.abs(crank.mechanismInertia(th) - r.J) / Jscale);
        c.states.forEach((st, s) => {
          crank.evaluate(th, st.omega, st.alpha, st.pressures, c.pCrankcase, out);
          eJ = Math.max(eJ, Math.abs(out.mechanismInertia - r.J) / Jscale, Math.abs(out.mechanismInertiaDerivative - r.dJ) / Jscale);
          eT = Math.max(eT, Math.abs(out.inertiaTorque - r.inertiaTorque[s]) / Tscale);
          eG = Math.max(eG, Math.abs(out.gravityTorque - r.gravityTorque));
          for (let i = 0; i < N; i++) eGas = Math.max(eGas, Math.abs(out.cylinders[i].gasTorque - r.gasTorque[s][i]) / Gscale);
          const a = crank.angularAcceleration(th, st.omega, st.pressures, c.pCrankcase, st.externalTorque, st.loadInertia);
          eA = Math.max(eA, Math.abs(a - r.alpha[s]) / Ascale);
          for (let i = 0; i < N; i++) eGas = Math.max(eGas, Math.abs(crank.gasTorques[i] - r.gasTorque[s][i]) / Gscale);
        });
      }
      expect(eX).toBeLessThan(1e-9);
      expect(eJ).toBeLessThan(1e-8);
      expect(eT).toBeLessThan(1e-8);
      expect(eG).toBeLessThan(1e-8); // N m
      expect(eGas).toBeLessThan(1e-9);
      expect(eA).toBeLessThan(1e-8);
    });
  }
});

describe('MultiCylinderCrankTrain identities', () => {
  const MT_GEOM: SliderCrankGeometry = MODEL_T.geometry;
  const OFFSETS = MODEL_T.layout.firingOffsetDeg;

  it('one cylinder at offset 0 with J_load = 0 is bit-identical to CrankTrainDynamics', () => {
    for (const pinOffset of [0, 0.004, -0.006]) {
      const kin = new SliderCrank({ ...CFR_GEOM, pinOffset }, 7);
      const single = new CrankTrainDynamics(kin, MASSES);
      const multi = new MultiCylinderCrankTrain(kin, MASSES, [0]);
      const p = new Float64Array(1);
      for (let d = -360; d < 360; d += 3.7) {
        const th = d * DEG;
        p[0] = syntheticPressure(kin, th);
        for (const [w, ext] of [[62.8, -12.3], [-40, 3.1], [157, 0]]) {
          expect(multi.angularAcceleration(th, w, p, 1.013e5, ext)).toBe(single.angularAcceleration(th, w, p[0], 1.013e5, ext));
          const tg = single.gasTorque(th, p[0], 1.013e5);
          expect(Math.abs(multi.gasTorques[0] - tg)).toBeLessThan(1e-13 * (1 + Math.abs(tg)));
          const h1 = single.holdingTorque(th, w, 17, p[0], 1.013e5, ext);
          expect(Math.abs(multi.holdingTorque(th, w, 17, p, 1.013e5, ext) - h1)).toBeLessThan(1e-12 * (1 + Math.abs(h1)));
        }
        expect(multi.update(th)[0]).toBe(kin.evaluate(th, newKinematicState()).dxdTheta);
        expect(Math.abs(multi.mechanismInertia(th) - single.mechanismInertia(th))).toBeLessThan(1e-15);
      }
    }
  });

  it('throw phases: offsets reduce mod 360° (2π-periodic mechanism) and equal cylinders share an evaluation', () => {
    const kin = new SliderCrank(MT_GEOM, 3.98);
    const a = new MultiCylinderCrankTrain(kin, MODEL_T.masses, OFFSETS);
    const b = new MultiCylinderCrankTrain(kin, MODEL_T.masses, [0, 180, 180, 0]);
    const c = new MultiCylinderCrankTrain(kin, MODEL_T.masses, [720, -180, 900, -360]);
    expect(a.throwGroups).toBe(2);
    expect(Array.from(a.throwPhase)).toEqual([0, Math.PI, Math.PI, 0]);
    expect(Array.from(c.throwPhase)).toEqual(Array.from(a.throwPhase));
    const single = new CrankTrainDynamics(kin, MODEL_T.masses);
    const p = new Float64Array(4);
    for (let d = -360; d < 360; d += 11) {
      const th = d * DEG;
      for (let i = 0; i < 4; i++) p[i] = cylinderPressure(kin, th, OFFSETS[i]);
      const aa = a.angularAcceleration(th, 120, p, 1e5, -15, 2.5);
      expect(b.angularAcceleration(th, 120, p, 1e5, -15, 2.5)).toBe(aa);
      expect(c.angularAcceleration(th, 120, p, 1e5, -15, 2.5)).toBe(aa);
      // explicit sum of four single-cylinder mechanisms
      let num = -15;
      let J = MODEL_T.masses.rotatingInertia + 2.5;
      const sc = newCrankTrainState();
      for (let i = 0; i < 4; i++) {
        const thi = th - a.throwPhase[i];
        single.evaluate(thi, 120, 0, p[i], 1e5, sc);
        num += sc.gasTorque + sc.gravityTorque - 0.5 * sc.mechanismInertiaDerivative * 120 * 120;
        J += sc.mechanismInertia;
      }
      expect(Math.abs(aa - num / J)).toBeLessThan(1e-12 * (1 + Math.abs(aa)));
    }
  });

  it('evaluate: Newton–Euler Σ crank torque = Σ (gas + inertia + gravity); gas torques match angularAcceleration', () => {
    const kin = new SliderCrank({ ...MT_GEOM, pinOffset: 0.002 }, 3.98);
    const crank = new MultiCylinderCrankTrain(kin, MODEL_T.masses, OFFSETS);
    const out = newMultiCylinderCrankTrainState(4);
    const p = new Float64Array(4);
    for (let d = -360; d < 360; d += 7) {
      const th = d * DEG;
      for (let i = 0; i < 4; i++) p[i] = cylinderPressure(kin, th, OFFSETS[i]);
      const al = crank.angularAcceleration(th, 150, p, 1e5, -20, 9.8);
      crank.evaluate(th, 150, al, p, 1e5, out);
      const sum = out.gasTorque + out.inertiaTorque + out.gravityTorque;
      expect(Math.abs(out.crankTorque - sum)).toBeLessThan(1e-9 * (1 + Math.abs(out.gasTorque)));
      expect(Math.abs(out.gasTorque - crank.gasTorqueSum)).toBeLessThan(1e-12 * (1 + Math.abs(out.gasTorque)));
      for (let i = 0; i < 4; i++) expect(out.cylinders[i].gasTorque).toBeCloseTo(crank.gasTorques[i], 9);
      expect(out.mechanismInertia).toBeCloseTo(crank.mechanismInertia(th), 14);
      // Lagrange: Σ T_inertia = −(Σ J_m α + ½ Σ J_m′ ω²) at the free acceleration
      expect(out.inertiaTorque).toBeCloseTo(-(out.mechanismInertia * al + 0.5 * out.mechanismInertiaDerivative * 150 * 150), 9);
      // the holding torque vanishes at the free acceleration
      expect(crank.holdingTorque(th, 150, al, p, 1e5, -20, 9.8)).toBeCloseTo(0, 9);
      const Th = crank.holdingTorque(th, 150, 0, p, 1e5, -20, 9.8);
      expect(crank.angularAcceleration(th, 150, p, 1e5, -20 + Th, 9.8)).toBeCloseTo(0, 10);
    }
  });

  it('inline four 1-2-4-3: Σ J_m, inertia and gravity torques are π-periodic (odd orders cancel, order 2 = 4× one cylinder)', () => {
    const kin = new SliderCrank(MT_GEOM, 3.98);
    const crank = new MultiCylinderCrankTrain(kin, MODEL_T.masses, OFFSETS);
    const single = new CrankTrainDynamics(kin, MODEL_T.masses);
    const out = newMultiCylinderCrankTrainState(4);
    const out2 = newMultiCylinderCrankTrainState(4);
    const sc = newCrankTrainState();
    const p = new Float64Array(4).fill(1e5);
    let jMin = Infinity;
    let jMax = -Infinity;
    let peak4 = 0;
    const w = (1800 * 2 * Math.PI) / 60;
    // Fourier sine/cosine coefficients (per revolution) of the inertia torque, orders 1–3
    const M = 720;
    const c4 = [0, 0, 0, 0];
    const s4 = [0, 0, 0, 0];
    const c1 = [0, 0, 0, 0];
    const s1 = [0, 0, 0, 0];
    for (let j = 0; j < M; j++) {
      const th = (2 * Math.PI * j) / M;
      crank.evaluate(th, w, 0, p, 1e5, out);
      crank.evaluate(th + Math.PI, w, 0, p, 1e5, out2);
      expect(out2.mechanismInertia).toBeCloseTo(out.mechanismInertia, 14);
      expect(out2.inertiaTorque).toBeCloseTo(out.inertiaTorque, 9);
      expect(out2.gravityTorque).toBeCloseTo(out.gravityTorque, 11);
      jMin = Math.min(jMin, out.mechanismInertia);
      jMax = Math.max(jMax, out.mechanismInertia);
      peak4 = Math.max(peak4, Math.abs(out.inertiaTorque));
      single.evaluate(th, w, 0, 1e5, 1e5, sc);
      for (let k = 1; k <= 3; k++) {
        c4[k] += (2 / M) * out.inertiaTorque * Math.cos(k * th);
        s4[k] += (2 / M) * out.inertiaTorque * Math.sin(k * th);
        c1[k] += (2 / M) * sc.inertiaTorque * Math.cos(k * th);
        s1[k] += (2 / M) * sc.inertiaTorque * Math.sin(k * th);
      }
    }
    const amp = (c: number[], s: number[], k: number) => Math.hypot(c[k], s[k]);
    console.info(
      `Model T Σ J_m = ${jMin.toFixed(4)}–${jMax.toFixed(4)} kg m² (J_rot ${MODEL_T.masses.rotatingInertia}); ` +
        `peak Σ inertia torque at 1800 rpm ${peak4.toFixed(1)} N m; orders 1/2/3: ` +
        `${amp(c4, s4, 1).toExponential(1)}/${amp(c4, s4, 2).toFixed(1)}/${amp(c4, s4, 3).toExponential(1)} N m ` +
        `(one cylinder ${amp(c1, s1, 1).toFixed(1)}/${amp(c1, s1, 2).toFixed(1)}/${amp(c1, s1, 3).toFixed(1)})`,
    );
    // two throw pairs: Σ J_m ≈ 4 m_rec a²/2 (1 ± …); small next to the 0.65 kg m² flywheel
    expect(jMin).toBeGreaterThan(0.003);
    expect(jMax).toBeLessThan(0.03);
    // odd orders cancel between the 0° and 180° throws, even orders add: order 2 = 4 × one cylinder's
    expect(amp(c4, s4, 1)).toBeLessThan(1e-9 * amp(c1, s1, 1));
    expect(amp(c4, s4, 3)).toBeLessThan(1e-9 * amp(c1, s1, 3) + 1e-9);
    expect(amp(c4, s4, 2) / amp(c1, s1, 2)).toBeCloseTo(4, 9);
  });

  it('cycle work: ∮ Σ T_gas dθ over 720° = N × one cylinder\'s ∮ (p − p_cc) dV (phase-shifted traces)', () => {
    const kin = new SliderCrank(MT_GEOM, 3.98);
    const crank = new MultiCylinderCrankTrain(kin, { ...MODEL_T.masses }, OFFSETS);
    const single = new CrankTrainDynamics(kin, MODEL_T.masses);
    const p = new Float64Array(4);
    const N = 14400;
    const h = (4 * Math.PI) / N;
    let w4 = 0;
    let w1 = 0;
    for (let i = 0; i < N; i++) {
      const th = -2 * Math.PI + i * h;
      for (let c = 0; c < 4; c++) p[c] = cylinderPressure(kin, th, OFFSETS[c]);
      crank.angularAcceleration(th, 100, p, 1e5, 0);
      w4 += crank.gasTorqueSum * h; // periodic integrand → rectangle rule is spectrally accurate
      w1 += single.gasTorque(th, cylinderPressure(kin, th, 0), 1e5) * h;
    }
    expect(Math.abs(w1)).toBeGreaterThan(50);
    expect(Math.abs(w4 - 4 * w1) / Math.abs(4 * w1)).toBeLessThan(1e-9);
  });

  it('free coasting with a reflected load inertia conserves kinetic + potential energy (RK4)', () => {
    const kin = new SliderCrank(MT_GEOM, 3.98);
    const crank = new MultiCylinderCrankTrain(kin, { ...MODEL_T.masses, rotatingInertia: 0.02 }, OFFSETS); // light → visible ripple
    const p = new Float64Array(4).fill(1e5);
    const Jl = 0.01;
    let th = 0.2;
    let w = 2 * Math.PI * 8;
    const acc = (t: number, om: number) => crank.angularAcceleration(t, om, p, 1e5, 0, Jl);
    const E0 = crank.kineticEnergy(th, w, Jl) + crank.potentialEnergy(th);
    const dt = 2e-5;
    let wMin = Infinity;
    let wMax = -Infinity;
    for (let i = 0; i < 20000; i++) {
      const k1t = w;
      const k1w = acc(th, w);
      const k2t = w + 0.5 * dt * k1w;
      const k2w = acc(th + 0.5 * dt * k1t, k2t);
      const k3t = w + 0.5 * dt * k2w;
      const k3w = acc(th + 0.5 * dt * k2t, k3t);
      const k4t = w + dt * k3w;
      const k4w = acc(th + dt * k3t, k4t);
      th += (dt / 6) * (k1t + 2 * k2t + 2 * k3t + k4t);
      w += (dt / 6) * (k1w + 2 * k2w + 2 * k3w + k4w);
      wMin = Math.min(wMin, w);
      wMax = Math.max(wMax, w);
    }
    const E1 = crank.kineticEnergy(th, w, Jl) + crank.potentialEnergy(th);
    expect((wMax - wMin) / w).toBeGreaterThan(0.05);
    expect(Math.abs(E1 - E0) / E0).toBeLessThan(1e-9);
  });

  it('fromSpec: Model T (4 cylinders, 2 throw groups), CFR (1 cylinder); layout/cylinder mismatch throws', () => {
    const mt = MultiCylinderCrankTrain.fromSpec(MODEL_T, new SliderCrank(MT_GEOM, 3.98));
    expect(mt.cylinders).toBe(4);
    expect(mt.throwGroups).toBe(2);
    expect(mt.rotatingInertia).toBe(MODEL_T.masses.rotatingInertia);
    const cfr = MultiCylinderCrankTrain.fromSpec(CFR_F1, new SliderCrank(CFR_F1.geometry, 7));
    expect(cfr.cylinders).toBe(1);
    expect(() => MultiCylinderCrankTrain.fromSpec({ ...MODEL_T, cylinders: 6 }, new SliderCrank(MT_GEOM, 3.98))).toThrow(RangeError);
    expect(() => MultiCylinderCrankTrain.fromSpec({ ...CFR_F1, cylinders: 2 }, new SliderCrank(CFR_F1.geometry, 7))).toThrow(RangeError);
    expect(() => new MultiCylinderCrankTrain(new SliderCrank(MT_GEOM, 3.98), MODEL_T.masses, [])).toThrow(RangeError);
    expect(() => new MultiCylinderCrankTrain(new SliderCrank(MT_GEOM, 3.98), MODEL_T.masses, [0, Number.NaN])).toThrow(RangeError);
  });

  it('hot path timing (printed): 4-cylinder angularAcceleration', () => {
    const kin = new SliderCrank(MT_GEOM, 3.98);
    const crank = new MultiCylinderCrankTrain(kin, MODEL_T.masses, OFFSETS);
    const p = new Float64Array([20e5, 1e5, 0.6e5, 3e5]);
    const n = 1_000_000;
    let acc = 0;
    for (let i = 0; i < 20000; i++) acc += crank.angularAcceleration(i * 1e-3, 120, p, 1e5, -5, 9.8);
    const t0 = performance.now();
    for (let i = 0; i < n; i++) acc += crank.angularAcceleration(i * 3.1e-6, 120, p, 1e5, -5, 9.8);
    const ns = ((performance.now() - t0) * 1e6) / n;
    console.info(`MultiCylinderCrankTrain.angularAcceleration (4 cyl, 2 throws): ${ns.toFixed(1)} ns/call`);
    expect(Number.isFinite(acc)).toBe(true);
    expect(ns).toBeLessThan(3000);
  });
});
