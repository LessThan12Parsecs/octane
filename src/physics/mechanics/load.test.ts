import { describe, expect, it } from 'vitest';
import { P_ATM, R_UNIVERSAL } from '../core/constants';
import type { VehicleSpec } from '../core/engine-spec';
import type { LoadSpec } from '../core/operating-point';
import { CFR_F1, CFR_FRICTION } from '../engines/cfr';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_FORD_WOT_TABLE, MODEL_T_VEHICLE } from '../engines/model-t';
import { CrankTrainDynamics, G_STANDARD, MultiCylinderCrankTrain } from './dynamics';
import { FrictionTorqueModel, pnhFmep } from './friction';
import { newKinematicState, SliderCrank } from './kinematics';
import { DRY_AIR_MOLAR_MASS, equilibriumOmega, LoadModel, type LoadOperatingPoint } from './load';

const MPH = 0.44704;
const LBFT = 1.3558179483314004;
const RPM = (2 * Math.PI) / 60;

const op = (load: LoadSpec | undefined, loadTorque = 0, T = 288.7, p = P_ATM): LoadOperatingPoint => ({
  load,
  loadTorque,
  ambientPressure: p,
  ambientTemperature: T,
});

/** Ford's WOT brake torque at the transmission output (linear interpolation of the table), N m. */
function fordTorque(omega: number): number {
  const rpm = omega / RPM;
  const t = MODEL_T_FORD_WOT_TABLE;
  if (rpm <= t[0][0]) return t[0][1] * LBFT;
  for (let i = 1; i < t.length; i++) {
    if (rpm <= t[i][0]) {
      const f = (rpm - t[i - 1][0]) / (t[i][0] - t[i - 1][0]);
      return (t[i - 1][1] + f * (t[i][1] - t[i - 1][1])) * LBFT;
    }
  }
  return t[t.length - 1][1] * LBFT;
}

describe('LoadModel: constant and brake', () => {
  it("'constant' (and absent load) returns loadTorque at any speed, no inertia, no road speed", () => {
    for (const ls of [undefined, { kind: 'constant' } as const]) {
      const m = new LoadModel().configure(op(ls, 7.5));
      expect(m.kind).toBe('constant');
      for (const w of [0, 10, 300, -50]) expect(m.torque(w)).toBe(7.5);
      expect(m.inertia()).toBe(0);
      expect(m.vehicleSpeed(100)).toBe(0);
    }
  });

  it("'brake': T = loadTorque·(|n|/refRpm)^e·sgn(ω), opposing rotation in both directions", () => {
    const m = new LoadModel().configure(op({ kind: 'brake', refRpm: 1200, exponent: 2 }, 40));
    expect(m.torque(1200 * RPM)).toBeCloseTo(40, 12);
    expect(m.torque(600 * RPM)).toBeCloseTo(10, 12);
    expect(m.torque(-600 * RPM)).toBeCloseTo(-10, 12);
    expect(m.torque(0)).toBe(0);
    const m15 = new LoadModel().configure(op({ kind: 'brake', refRpm: 1000, exponent: 1.5 }, 20));
    expect(m15.torque(4000 * RPM)).toBeCloseTo(20 * 8, 10);
    const m0 = new LoadModel().configure(op({ kind: 'brake', refRpm: 1000, exponent: 0 }, 20));
    expect(m0.torque(3 * RPM)).toBe(20);
    expect(m0.torque(-3 * RPM)).toBe(-20);
    expect(m.inertia()).toBe(0);
    expect(() => new LoadModel().configure(op({ kind: 'brake', refRpm: 0, exponent: 2 }, 1))).toThrow(RangeError);
    expect(() => new LoadModel().configure(op({ kind: 'brake', refRpm: 100, exponent: -1 }, 1))).toThrow(RangeError);
  });
});

describe('LoadModel: vehicle road load through the gearing', () => {
  const V: VehicleSpec = {
    mass: 1000,
    wheelRadius: 0.35,
    finalDrive: 4,
    gears: { low: 3, high: 1, reverse: 4 },
    drivelineEfficiency: 0.8,
    rollingResistance: 0.015,
    dragArea: 2,
  };

  it('level road, by hand: v = ωr/(G G_f), T = (m g C_rr + ½ρC_dA v²)·r/(G G_f η), J = m r²/(G G_f)²', () => {
    const T = 293.15;
    const m = new LoadModel(V).configure(op({ kind: 'vehicle', gear: 'low', grade: 0 }, 99, T));
    const rho = (P_ATM * DRY_AIR_MOLAR_MASS) / (R_UNIVERSAL * T);
    expect(rho).toBeCloseTo(1.2041, 3); // dry air at 20 °C, 1 atm
    expect(m.airDensity).toBe(rho);
    const ratio = 3 * 4;
    for (const rpm of [500, 1500, 3000]) {
      const w = rpm * RPM;
      const v = (w * 0.35) / ratio;
      expect(m.vehicleSpeed(w)).toBeCloseTo(v, 12);
      const sg = v / Math.sqrt(v * v + 0.05 * 0.05);
      const F = 1000 * G_STANDARD * 0.015 * sg + 0.5 * rho * 2 * v * v;
      expect(m.roadLoadForce(v)).toBeCloseTo(F, 9);
      expect(m.torque(w)).toBeCloseTo((F * 0.35) / (ratio * 0.8), 10);
      // power balance: engine power × η = road-load power
      expect(m.power(w) * 0.8).toBeCloseTo(F * v, 8);
    }
    expect(m.inertia()).toBeCloseTo((1000 * 0.35 * 0.35) / (ratio * ratio), 14);
    // loadTorque is ignored by the vehicle load
    expect(new LoadModel(V).configure(op({ kind: 'vehicle', gear: 'low', grade: 0 }, 0, T)).torque(100)).toBe(m.torque(100));
  });

  it('grade: F = m g (C_rr cos a + sin a) with a = atan(grade); kinetic energy of the reflected inertia', () => {
    const wi = 4.0;
    const m = new LoadModel(V, { wheelInertia: wi, airDensity: 1.2 }).configure(op({ kind: 'vehicle', gear: 'high', grade: 0.08 }));
    const a = Math.atan(0.08);
    const v = 15;
    const F = 1000 * G_STANDARD * (0.015 * Math.cos(a) * (v / Math.hypot(v, 0.05)) + Math.sin(a)) + 0.5 * 1.2 * 2 * v * v;
    expect(m.roadLoadForce(v)).toBeCloseTo(F, 9);
    const w = (v * 4) / 0.35;
    expect(m.vehicleSpeed(w)).toBeCloseTo(v, 12);
    // ½ J_L ω² = ½ m v² + ½ J_w ω_w²
    expect(0.5 * m.inertia() * w * w).toBeCloseTo(0.5 * 1000 * v * v + 0.5 * wi * (v / 0.35) ** 2, 8);
  });

  it('driveline efficiency divides the torque while driving and multiplies it on overrun (losses always dissipate)', () => {
    const down = new LoadModel(V, { airDensity: 1.2 }).configure(op({ kind: 'vehicle', gear: 'high', grade: -0.15 }));
    const lossless = (mm: LoadModel, w: number) => mm.roadLoadForce(mm.vehicleSpeed(w)) * mm.roadSpeedPerOmega;
    for (const rpm of [800, 1500, 2500, 3500]) {
      const w = rpm * RPM;
      const tw = lossless(down, w);
      const t = down.torque(w);
      if (tw > 0) expect(t).toBeCloseTo(tw / 0.8, 10);
      else expect(t).toBeCloseTo(tw * 0.8, 10);
      // the engine never receives more power than the wheels deliver, nor delivers less than they absorb
      const pw = down.roadLoadForce(down.vehicleSpeed(w)) * down.vehicleSpeed(w);
      const pe = down.power(w);
      expect(pe - pw).toBeGreaterThanOrEqual(-1e-9 * Math.abs(pw)); // P_engine − P_wheel = driveline loss ≥ 0
    }
    // 15 % down: overrun at low speed (car drives the engine), driving at high speed (drag wins)
    expect(down.torque(800 * RPM)).toBeLessThan(0);
    expect(down.torque(5000 * RPM)).toBeGreaterThan(0);
  });

  it('reverse drives the car backwards (v < 0) and still loads the engine; neutral is declutched', () => {
    const rev = new LoadModel(V, { airDensity: 1.2 }).configure(op({ kind: 'vehicle', gear: 'reverse', grade: 0 }));
    const w = 1000 * RPM;
    expect(rev.vehicleSpeed(w)).toBeCloseTo((-w * 0.35) / 16, 12);
    expect(rev.torque(w)).toBeGreaterThan(0);
    expect(rev.power(w) * 0.8).toBeCloseTo(rev.roadLoadForce(rev.vehicleSpeed(w)) * rev.vehicleSpeed(w), 8);
    // in reverse on a 20 % (forward-uphill) grade the car backs DOWN the slope: gravity drives the engine (overrun)
    const revUp = new LoadModel(V, { airDensity: 1.2 }).configure(op({ kind: 'vehicle', gear: 'reverse', grade: 0.2 }));
    expect(revUp.torque(200 * RPM)).toBeLessThan(0);
    const n = new LoadModel(V).configure(op({ kind: 'vehicle', gear: 'neutral', grade: 0.1 }, 50));
    expect(n.kind).toBe('neutral');
    expect(n.torque(150)).toBe(0);
    expect(n.inertia()).toBe(0);
    expect(n.vehicleSpeed(150)).toBe(0);
  });

  it('continuous in ω through 0 on a level road (rolling-resistance sign smoothing)', () => {
    const m = new LoadModel(V, { airDensity: 1.2 }).configure(op({ kind: 'vehicle', gear: 'low', grade: 0 }));
    let prev = m.torque(-1);
    for (let w = -1; w <= 1; w += 1e-3) {
      const t = m.torque(w);
      expect(Math.abs(t - prev)).toBeLessThan(0.5); // N m per 1e-3 rad/s
      expect(t).toBeCloseTo(-m.torque(-w), 12);
      prev = t;
    }
  });

  it('validation: vehicle load without a vehicle, unknown gear, invalid data', () => {
    expect(() => new LoadModel().configure(op({ kind: 'vehicle', gear: 'high', grade: 0 }))).toThrow(RangeError);
    expect(() => new LoadModel(V).configure(op({ kind: 'vehicle', gear: 'overdrive', grade: 0 }))).toThrow(RangeError);
    expect(() => new LoadModel({ ...V, drivelineEfficiency: 0 }).configure(op({ kind: 'vehicle', gear: 'high', grade: 0 }))).toThrow(RangeError);
    expect(() => new LoadModel(V).configure(op({ kind: 'vehicle', gear: 'high', grade: Number.NaN }))).toThrow(RangeError);
    expect(() => new LoadModel(V, { rollingSmoothingSpeed: 0 })).toThrow(RangeError);
  });
});

describe('Model T vehicle (MODEL_T_VEHICLE) against period road data', () => {
  it('high gear: 40–41 engine rpm per mph [Motor Age 1917: 1230 rpm at 30 mph; Ford table: 40]; low ≈ 14 mph at 1600 rpm', () => {
    const high = new LoadModel(MODEL_T_VEHICLE).configure(MODEL_T_CRUISE);
    const rpmPerMph = MPH / high.roadSpeedPerOmega / RPM;
    console.info(`Model T high gear: ${rpmPerMph.toFixed(2)} rpm/mph, J_load ${high.inertia().toFixed(2)} kg m² (J_rot ${MODEL_T.masses.rotatingInertia})`);
    expect(rpmPerMph).toBeGreaterThan(40);
    expect(rpmPerMph).toBeLessThan(41.2);
    const low = new LoadModel(MODEL_T_VEHICLE).configure({ ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'low', grade: 0 } });
    expect(low.vehicleSpeed(1600 * RPM) / MPH).toBeCloseTo(14.2, 0);
    expect(low.inertia() / high.inertia()).toBeCloseTo(1 / 2.75 ** 2, 12);
  });

  it('top speed in high with Ford\'s WOT torque curve: 42–45 mph at ≈ 1700–1850 rpm', () => {
    const high = new LoadModel(MODEL_T_VEHICLE).configure(MODEL_T_CRUISE);
    const w = equilibriumOmega(fordTorque, high, 1000 * RPM, 1900 * RPM);
    const mph = high.vehicleSpeed(w) / MPH;
    console.info(`Model T top speed (Ford WOT curve, ${MODEL_T_VEHICLE.mass} kg, C_rr ${MODEL_T_VEHICLE.rollingResistance}, C_dA ${MODEL_T_VEHICLE.dragArea}): ${mph.toFixed(1)} mph at ${(w / RPM).toFixed(0)} rpm`);
    expect(mph).toBeGreaterThan(42);
    expect(mph).toBeLessThan(45);
    expect(w / RPM).toBeGreaterThan(1700);
    expect(w / RPM).toBeLessThan(1850);
    // stable equilibrium: the load torque rises faster than the engine torque through it
    const h = 5 * RPM;
    expect(high.torque(w + h) - high.torque(w - h)).toBeGreaterThan(fordTorque(w + h) - fordTorque(w - h));
  });

  it("reproduces Tulsa MTFC's Ford-curve road calculation (850 kg, C_rr 0.01, C_dA 1.86 m², 20 % driveline loss)", () => {
    // Tulsa 'Power and Torque' (tildentechnologies.com/mtfctulsa): with the Ford curve, 45 mph level, 35 mph on
    // 5 %, 30 mph on 6.7 %, 8.1 % maximum grade in high (1875 lb); 6.9 % at 2200 lb. Their air density and
    // tyre radius are not stated: MODEL_T_VEHICLE's 41 rpm/mph and 1.225 kg/m³ assumed here (with Ford's
    // 40 rpm/mph the grade results drop by ≈ 0.3 % / 0.7 mph).
    const tulsa: VehicleSpec = { ...MODEL_T_VEHICLE, mass: 850, rollingResistance: 0.01, dragArea: 1.86, drivelineEfficiency: 0.8 };
    const speedOn = (grade: number, veh = tulsa) => {
      const m = new LoadModel(veh, { airDensity: 1.225 }).configure(op({ kind: 'vehicle', gear: 'high', grade }));
      return m.vehicleSpeed(equilibriumOmega(fordTorque, m, 950 * RPM, 1900 * RPM)) / MPH;
    };
    const maxGrade = (veh: VehicleSpec) => {
      let lo = 0;
      let hi = 0.3;
      for (let k = 0; k < 40; k++) {
        const g = 0.5 * (lo + hi);
        const m = new LoadModel(veh, { airDensity: 1.225 }).configure(op({ kind: 'vehicle', gear: 'high', grade: g }));
        let best = -Infinity;
        for (let rpm = 300; rpm <= 1900; rpm += 5) best = Math.max(best, fordTorque(rpm * RPM) - m.torque(rpm * RPM));
        if (best >= 0) lo = g;
        else hi = g;
      }
      return lo;
    };
    const s0 = speedOn(0);
    const s5 = speedOn(0.05);
    const s67 = speedOn(0.067);
    const g850 = maxGrade(tulsa);
    const g998 = maxGrade({ ...tulsa, mass: 998 });
    console.info(
      `Tulsa check: level ${s0.toFixed(1)} mph (45), 5 % ${s5.toFixed(1)} (35), 6.7 % ${s67.toFixed(1)} (30), ` +
        `max grade ${(100 * g850).toFixed(1)} % (8.1), at 998 kg ${(100 * g998).toFixed(1)} % (6.9)`,
    );
    expect(Math.abs(s0 - 45)).toBeLessThan(1.5);
    expect(Math.abs(s5 - 35)).toBeLessThan(1.5);
    expect(Math.abs(s67 - 30)).toBeLessThan(1.5);
    expect(Math.abs(100 * g850 - 8.1)).toBeLessThan(0.4);
    expect(Math.abs(100 * g998 - 6.9)).toBeLessThan(0.4);
  });

  it('free running in high gear settles on the equilibrium speed (RK4 on the reflected-inertia crank equation)', () => {
    const high = new LoadModel(MODEL_T_VEHICLE).configure({ ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'high', grade: 0.02 } });
    const J = MODEL_T.masses.rotatingInertia + high.inertia();
    const f = (w: number) => (fordTorque(w) - high.torque(w)) / J;
    const wEq = equilibriumOmega(fordTorque, high, 600 * RPM, 1900 * RPM);
    let w = 900 * RPM;
    const dt = 0.05;
    for (let i = 0; i < 4000; i++) {
      const k1 = f(w);
      const k2 = f(w + 0.5 * dt * k1);
      const k3 = f(w + 0.5 * dt * k2);
      const k4 = f(w + dt * k3);
      w += (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    expect(Math.abs(w - wEq) / wEq).toBeLessThan(1e-6);
    expect(high.vehicleSpeed(w) / MPH).toBeGreaterThan(30);
  });

  it('hot path timing (printed): vehicle torque()', () => {
    const m = new LoadModel(MODEL_T_VEHICLE).configure(MODEL_T_CRUISE);
    const n = 2_000_000;
    let acc = 0;
    const t0 = performance.now();
    for (let i = 0; i < n; i++) acc += m.torque(100 + (i & 1023) * 0.1);
    const ns = ((performance.now() - t0) * 1e6) / n;
    console.info(`LoadModel.torque (vehicle): ${ns.toFixed(1)} ns/call`);
    expect(Number.isFinite(acc)).toBe(true);
    expect(ns).toBeLessThan(1000);
  });
});

describe('Free-speed right-hand side composed from the mechanics primitives', () => {
  it('CFR (one cylinder, constant load): crank + friction + load is bit-identical to the single-cylinder RHS', () => {
    const kin = new SliderCrank(CFR_F1.geometry, 7);
    const dyn = new CrankTrainDynamics(kin, CFR_F1.masses);
    const fr1 = new FrictionTorqueModel(kin);
    const crank = MultiCylinderCrankTrain.fromSpec(CFR_F1, kin);
    const frN = new FrictionTorqueModel(kin, 0.5, CFR_F1.cylinders);
    const load = new LoadModel(CFR_F1.vehicle).configure(op({ kind: 'constant' }, 2));
    const b = pnhFmep(CFR_FRICTION, 600, 0.95e5, 1e5);
    fr1.setFromPnh(b, 600 * RPM);
    frN.setFromPnh(b, 600 * RPM);
    const pCyl = new Float64Array(1);
    const pcc = 1e5 - 858;
    const ks = newKinematicState();
    for (let d = -360; d < 360; d += 1.3) {
      const th = (d * Math.PI) / 180;
      pCyl[0] = 1e5 + 30e5 * Math.exp(-(((d - 15) / 30) ** 2));
      for (const om of [600 * RPM, 900 * RPM, 6.3]) {
        // today's cycle-model.ts free-mode RHS: ks = kin.evaluate(θ), ext = −loadTorque + friction.torque(ks.dxdTheta, ω)
        kin.evaluate(th, ks);
        const legacy = dyn.angularAcceleration(th, om, pCyl[0], pcc, -2 + fr1.torque(ks.dxdTheta, om));
        const dx = crank.update(th);
        const now = crank.angularAcceleration(th, om, pCyl, pcc, -load.torque(om) + frN.torqueCylinders(dx, om), load.inertia());
        expect(now).toBe(legacy);
      }
    }
  });

  it('Model T in high gear: the reflected car dominates the crank inertia (≈ 15× the flywheel)', () => {
    const kin = new SliderCrank(MODEL_T.geometry, MODEL_T.geometry.compressionRatio);
    const crank = MultiCylinderCrankTrain.fromSpec(MODEL_T, kin);
    const load = new LoadModel(MODEL_T.vehicle).configure(MODEL_T_CRUISE);
    const Jm = crank.mechanismInertia(0.3);
    const ratio = load.inertia() / MODEL_T.masses.rotatingInertia;
    expect(ratio).toBeGreaterThan(10);
    expect(ratio).toBeLessThan(20);
    expect(Jm / (MODEL_T.masses.rotatingInertia + load.inertia())).toBeLessThan(0.003);
  });
});
