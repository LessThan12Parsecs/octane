/**
 * The car of a 'vehicle' load in free speed (Ford Model T; code review 'gasx-mech-load' / 'cycle'):
 *  - Road speed is a state: k ω in gear; declutched the car coasts on its own road load, m dv/dt = −F_road(v)
 *    (oracle, level road, sgn_ε ≈ 1: v(t) = √(a/b)·tan(atan(v₀√(b/a)) − √(ab)·t), a = g C_rr, b = ½ρC_dA/m).
 *  - A gear change is an inelastic clutch engagement through the rigid gear train: the generalised momentum
 *    J_e ω + J_L v/k is conserved, ω′ = (J_e ω + J_L v/k)/(J_e + J_L) (J_e = J_rot + ΣJ_m), and the clutch
 *    dissipates ½J_e J_L/(J_e + J_L)(ω − v/k)² (clutchLoss ledger); leaving gear keeps ω and the car's speed.
 *    Round 1 kept ω and made the car's speed jump: a downshift at cruise destroyed ≈ 48 kJ of the car's kinetic
 *    energy, and the engine summary of that cycle reported brake 3844 N m, η_b 7.6.
 *  - The engine summary's brake torque is the ENGINE's output (ΔE_kin of J_rot + ΣJ_m only): accelerating at
 *    wide-open throttle in low gear brake ≈ indicated − friction − J_e α (round 1: ≈ 6 N m = the road load, BSFC
 *    ≈ 10 kg/kWh), BSFC plausible, η_b < η_i; the energy balance closes with the car's inertia term,
 *    brake ≈ load + loadInertia, also over cycles containing gear changes.
 * The model is uncalibrated: the checks are mechanics and bookkeeping, not Ford's numbers. Runtime ≈ 30–60 s.
 */
import { describe, expect, it } from 'vitest';
import type { CycleSummary, EngineCycleSummary } from '../core/snapshot';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_VEHICLE } from '../engines/model-t';
import { G_STANDARD } from '../mechanics/dynamics';
import { FREE_MODE_MIN_OMEGA } from './cycle-common';
import { CycleModel, I_OM } from './cycle-model';

const rel = (a: number, b: number): number => Math.abs(a / b - 1);
const MPH = 0.44704;

/** Car kinetic energy ½ m v² (MODEL_T_VEHICLE, no wheel inertia: LoadModel default), J. */
const carKE = (v: number): number => 0.5 * MODEL_T_VEHICLE.mass * v * v;

/** Engine summary of the summaries completed by one runCycles call. */
const engineOf = (s: CycleSummary[]): EngineCycleSummary => s.find((x) => x.engine)!.engine!;

/** Step to engine angle th (deg) in the current engine cycle (stepUntil also returns at spark / knock interrupts). */
function stepTo(m: CycleModel, th: number): void {
  const c = m.cycle;
  while (m.cycle === c && m.theta < th - 1e-9) m.stepUntil(Infinity, th);
  expect(m.cycle).toBe(c);
}

/** The energy balance of a free-speed vehicle cycle: brake ≈ load + loadInertia (friction mean at the mean speed). */
function closure(e: EngineCycleSummary): number {
  return Math.abs(e.brakeTorque - (e.loadTorque + e.loadInertiaTorque!)) / e.frictionTorque;
}

describe('clutch engagement in free speed: angular momentum through the gear train, loss to the clutch ledger', () => {
  const m = new CycleModel(MODEL_T, MODEL_T_CRUISE, { warmupCycles: 1 });
  m.runCycles(1);

  it('high → low at cruise (mid-cycle): ω′ from momentum conservation, the car keeps its momentum share, ledger = ΔE_kin', () => {
    stepTo(m, 0);
    const w = m.y[I_OM];
    const v = m.vehicleSpeed();
    expect(rel(v, m.load.roadSpeedPerOmega * w)).toBeLessThan(1e-15);
    expect(v / MPH).toBeGreaterThan(10);
    const Je = m.engineInertia();
    expect(Je).toBeGreaterThan(MODEL_T.masses.rotatingInertia);
    const loss0 = m.clutchLoss;
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'low', grade: 0 } });
    const k = m.load.roadSpeedPerOmega;
    const JL = m.load.inertia();
    const w1 = m.y[I_OM];
    // momentum: (J_e + J_L) ω′ = J_e ω + J_L v/k
    expect(rel((Je + JL) * w1, Je * w + JL * (v / k))).toBeLessThan(1e-12);
    expect(rel(m.vehicleSpeed(), k * w1)).toBeLessThan(1e-15);
    // the engine is spun up, the car slowed (no free kinetic energy): v′ < v, ω′ > ω
    expect(w1).toBeGreaterThan(w);
    expect(m.vehicleSpeed()).toBeLessThan(v);
    // energy: ½J_e ω² + ½ m v² = ½(J_e + J_L) ω′² + loss, loss = ½ J_e J_L/(J_e + J_L) (ω − v/k)² > 0
    const loss = m.clutchLoss - loss0;
    const keBefore = 0.5 * Je * w * w + carKE(v);
    const keAfter = 0.5 * (Je + JL) * w1 * w1;
    expect(Math.abs(loss - (keBefore - keAfter))).toBeLessThan(1e-9 * keBefore);
    expect(rel(loss, ((0.5 * Je * JL) / (Je + JL)) * (w - v / k) ** 2)).toBeLessThan(1e-9);
    expect(loss).toBeGreaterThan(1000);
    // the summary of the cycle with the shift: no spike (round 1: brake 3844 N m, η_b 7.6), the balance closes
    const e = engineOf(m.runCycles(1));
    expect(e.clutchLoss!).toBeCloseTo(loss, 6);
    expect(e.brakeEfficiency).toBeGreaterThan(0);
    expect(e.brakeEfficiency).toBeLessThan(0.4);
    expect(e.brakeTorque).toBeLessThan(e.indicatedTorque);
    // (the residual is the summary's friction approximation — the cycle mean at the cycle-mean speed — which a
    // 2× speed jump within the cycle puts ≈ 8 % off; steady or accelerating cycles close to < 1e-3)
    expect(closure(e)).toBeLessThan(0.15);
    console.log(
      `[clutch] high → low at ${(v / MPH).toFixed(1)} mph: ${((w * 30) / Math.PI).toFixed(0)} → ${((w1 * 30) / Math.PI).toFixed(0)} rpm, ` +
        `car ${(v / MPH).toFixed(2)} → ${((k * w1) / MPH).toFixed(2)} mph, clutch loss ${loss.toFixed(0)} J; cycle: brake ${e.brakeTorque.toFixed(1)} N m, ` +
        `load ${e.loadTorque.toFixed(1)} + inertia ${e.loadInertiaTorque!.toFixed(1)} N m (residual ${closure(e).toExponential(1)} of friction), η_b ${e.brakeEfficiency.toFixed(3)}`,
    );
  });

  it('neutral: ω and the car speed continue, the car coasts on its road load (analytic), the summary reports the coasting speed', () => {
    stepTo(m, -200);
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'high', grade: 0 } });
    m.runCycles(1);
    stepTo(m, -200);
    const w = m.y[I_OM];
    const v0 = m.vehicleSpeed();
    const loss0 = m.clutchLoss;
    // declutch and open the throttle wide: the engine revs freely while the car coasts
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'neutral', grade: 0 }, throttle: 1 });
    expect(m.y[I_OM]).toBe(w);
    expect(m.vehicleSpeed()).toBe(v0);
    expect(m.clutchLoss).toBe(loss0);
    expect(m.load.inertia()).toBe(0);
    const t0 = m.t;
    stepTo(m, 300);
    const dt = m.t - t0;
    // level-road coasting: dv/dt = −(a + b v²) (sgn_ε(v) = 1 − O(ε²/v²) ≈ 1 − 1e-5 at 10 m/s)
    const a = G_STANDARD * MODEL_T_VEHICLE.rollingResistance;
    const b = (0.5 * m.load.airDensity * MODEL_T_VEHICLE.dragArea) / MODEL_T_VEHICLE.mass;
    const vExact = Math.sqrt(a / b) * Math.tan(Math.atan(v0 * Math.sqrt(b / a)) - Math.sqrt(a * b) * dt);
    expect(m.vehicleSpeed()).toBeLessThan(v0);
    expect(Math.abs(m.vehicleSpeed() - vExact)).toBeLessThan(1e-6 * v0);
    // re-engaging high: momentum conservation with the COASTING car speed
    const wN = m.y[I_OM];
    const vN = m.vehicleSpeed();
    const Je = m.engineInertia();
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'high', grade: 0 }, throttle: 0.35 });
    const k = m.load.roadSpeedPerOmega;
    const JL = m.load.inertia();
    expect(wN).toBeGreaterThan(1.05 * (vN / k)); // the engine has run away from the car
    expect(rel((Je + JL) * m.y[I_OM], Je * wN + JL * (vN / k))).toBeLessThan(1e-12);
    const loss = m.clutchLoss - loss0;
    expect(Math.abs(loss - (0.5 * Je * wN * wN + carKE(vN) - 0.5 * (Je + JL) * m.y[I_OM] ** 2))).toBeLessThan(1e-9 * carKE(vN));
    expect(rel(loss, ((0.5 * Je * JL) / (Je + JL)) * (wN - vN / k) ** 2)).toBeLessThan(1e-9);
    expect(loss).toBeGreaterThan(10);
    // this cycle: high (−360 → −200), neutral (−200 → 300), high (300 → 360): the balance still closes
    const e = engineOf(m.runCycles(1));
    expect(e.clutchLoss!).toBeCloseTo(loss, 6);
    expect(closure(e)).toBeLessThan(0.05);
    expect(e.brakeEfficiency).toBeLessThan(0.4);
    expect(e.vehicleSpeed!).toBeGreaterThan(0.9 * vN);
    console.log(`[clutch] high → neutral → high: clutch loss ${loss.toFixed(0)} J, brake ${e.brakeTorque.toFixed(1)} N m, residual ${closure(e).toExponential(1)} of friction`);
    // a full cycle in neutral: the summary reports the coasting car (round 1: 0)
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'neutral', grade: 0 } });
    const v1 = m.vehicleSpeed();
    const eN = engineOf(m.runCycles(1));
    expect(eN.vehicleSpeed!).toBeGreaterThan(0.95 * v1);
    expect(eN.vehicleSpeed!).toBeLessThan(v1);
    expect(eN.loadTorque).toBe(0);
    expect(eN.loadInertiaTorque!).toBe(0);
    expect(eN.clutchLoss!).toBe(0);
  });

  it('reverse while rolling forward: the stall guard holds ω′ at FREE_MODE_MIN_RPM, the ledger books the exact energy change', () => {
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'high', grade: 0 } });
    stepTo(m, 100);
    const w = m.y[I_OM];
    const v = m.vehicleSpeed();
    const Je = m.engineInertia();
    const loss0 = m.clutchLoss;
    expect(v).toBeGreaterThan(1);
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'reverse', grade: 0 } });
    expect(m.y[I_OM]).toBe(FREE_MODE_MIN_OMEGA);
    expect(m.vehicleSpeed()).toBeLessThan(0);
    const JL = m.load.inertia();
    const loss = m.clutchLoss - loss0;
    expect(Math.abs(loss - (0.5 * Je * w * w + carKE(v) - 0.5 * (Je + JL) * FREE_MODE_MIN_OMEGA ** 2))).toBeLessThan(1e-9 * carKE(v));
  });
});

describe('gear changes outside a free-running car', () => {
  it('fixed speed: the dynamometer holds ω, the car follows at k ω, no clutch loss', () => {
    const m = new CycleModel(MODEL_T, { ...MODEL_T_CRUISE, speedMode: 'fixed', rpm: 1200 }, { warmupCycles: 0 });
    stepTo(m, -300);
    const w = m.y[I_OM];
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'low', grade: 0 } });
    expect(m.y[I_OM]).toBe(w);
    expect(m.clutchLoss).toBe(0);
    expect(m.vehicleSpeed()).toBe(m.load.roadSpeedPerOmega * w);
    // a dynamometer load → a car in gear in free speed: a new set-up, the car rolls at k ω (no clutch event)
    m.setOperatingPoint({ load: { kind: 'brake', refRpm: 1000, exponent: 2 }, loadTorque: 40, speedMode: 'free' });
    const w2 = m.y[I_OM];
    m.setOperatingPoint({ load: { kind: 'vehicle', gear: 'high', grade: 0 } });
    expect(m.y[I_OM]).toBe(w2);
    expect(m.clutchLoss).toBe(0);
    expect(m.vehicleSpeed()).toBe(m.load.roadSpeedPerOmega * w2);
  });
});

describe('engine summary: brake torque is the engine output (wide-open throttle in low gear, the car accelerating)', () => {
  it('brake ≈ indicated − friction − J_e α; load + loadInertia closes the balance; BSFC plausible; η_b < η_i', () => {
    const m = new CycleModel(MODEL_T, { ...MODEL_T_CRUISE, rpm: 800, throttle: 1, load: { kind: 'vehicle', gear: 'low', grade: 0 } }, { warmupCycles: 1 });
    for (let k = 0; k < 2; k++) {
      const w0 = m.y[I_OM];
      const Je = m.engineInertia(); // at θ = −360 (the same phase as the cycle end)
      const sums = m.runCycles(1);
      const w1 = m.y[I_OM];
      const e = engineOf(sums);
      const om = (e.rpmMean * 2 * Math.PI) / 60;
      const dt = (4 * Math.PI) / om;
      const alpha = (w1 - w0) / dt;
      expect(alpha).toBeGreaterThan(10); // the car accelerates (rad/s² at the crank)
      const JeAlpha = Je * alpha;
      // the definition: indicated − friction − ½J_e(ω₁² − ω₀²)/4π (no shift in this cycle) ...
      expect(Math.abs(e.brakeTorque - (e.indicatedTorque - e.frictionTorque - (0.5 * Je * (w1 * w1 - w0 * w0)) / (4 * Math.PI)))).toBeLessThan(1e-9 * e.indicatedTorque);
      // ... ≈ indicated − friction − J_e α (½Δω²/4π = J_e ω̄_arith Δω/(ω̄_time dt))
      expect(Math.abs(e.brakeTorque - (e.indicatedTorque - e.frictionTorque - JeAlpha))).toBeLessThan(0.02 * JeAlpha + 0.1);
      // most of the output accelerates the car: brake ≫ road load, and brake = load + loadInertia
      expect(e.brakeTorque).toBeGreaterThan(5 * e.loadTorque);
      expect(e.loadInertiaTorque!).toBeGreaterThan(0.5 * e.brakeTorque);
      expect(closure(e)).toBeLessThan(5e-3);
      // BSFC 250–1000 g/kWh (period Ford ≈ 0.7–1 lb/hp·h ≈ 425–600 g/kWh; uncalibrated model); round 1: ≈ 10 000
      const bsfc = e.bsfc * 3.6e9;
      expect(bsfc).toBeGreaterThan(250);
      expect(bsfc).toBeLessThan(1000);
      // brake efficiency below the cylinders' indicated efficiency
      const etaI = sums.reduce((a, s) => a + s.indicatedEfficiency, 0) / sums.length;
      expect(e.brakeEfficiency).toBeGreaterThan(0.4 * etaI);
      expect(e.brakeEfficiency).toBeLessThan(etaI);
      console.log(
        `[WOT low gear] ${((w0 * 30) / Math.PI).toFixed(0)} → ${((w1 * 30) / Math.PI).toFixed(0)} rpm: indicated ${e.indicatedTorque.toFixed(1)}, ` +
          `friction ${e.frictionTorque.toFixed(1)}, J_e α ${JeAlpha.toFixed(1)}, brake ${e.brakeTorque.toFixed(1)} N m = load ${e.loadTorque.toFixed(1)} + car ` +
          `${e.loadInertiaTorque!.toFixed(1)} (residual ${closure(e).toExponential(1)} of friction); BSFC ${bsfc.toFixed(0)} g/kWh, η_b ${e.brakeEfficiency.toFixed(3)} (η_i ${etaI.toFixed(3)})`,
      );
    }
  });
});
