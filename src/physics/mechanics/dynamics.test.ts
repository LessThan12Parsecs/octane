import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/mechanics_dynamics.json';
import { DEG } from '../core/constants';
import { CrankTrainDynamics, newCrankTrainState, type CrankTrainMasses } from './dynamics';
import { SliderCrank, type SliderCrankGeometry } from './kinematics';

interface DynCase {
  name: string;
  geometry: SliderCrankGeometry;
  masses: CrankTrainMasses;
  gravity: number;
  states: { omega: number; alpha: number }[];
  rows: { thetaDeg: number; J: number; dJ: number; gravityTorque: number; inertiaTorque: number[] }[];
}
const cases = (fixture as { cases: DynCase[] }).cases;

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
