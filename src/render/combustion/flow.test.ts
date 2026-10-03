import { describe, expect, it } from 'vitest';
import { CURTAIN_DISCHARGE_COEFF, MAX_SWIRL_RATIO } from './constants';
import { curtainJetSpeed, GasFlowModel, MAX_JET_SPEED } from './flow';
import { snap, testSpec, TEST_R } from './test-utils';

describe('curtain jet speed', () => {
  it('v = ṁ / (ρ C_d A), clamped', () => {
    expect(curtainJetSpeed(0.01, 1, 1e-4)).toBeCloseTo(0.01 / (1 * CURTAIN_DISCHARGE_COEFF * 1e-4), 9);
    expect(curtainJetSpeed(-0.01, 1, 1e-4)).toBeCloseTo(0.01 / (1 * CURTAIN_DISCHARGE_COEFF * 1e-4), 9);
    expect(curtainJetSpeed(1, 1, 1e-6)).toBe(MAX_JET_SPEED);
    expect(curtainJetSpeed(0.01, 1, 0)).toBe(0);
  });
});

describe('GasFlowModel', () => {
  it('axial velocity: zero at the head, piston velocity at the crown', () => {
    const f = new GasFlowModel(testSpec());
    f.update(snap({ t: 0, clearanceHeight: 0.02 }), 0);
    f.update(snap({ t: 1e-4, clearanceHeight: 0.0199 }), 1e-4); // piston rising at 1 m/s
    expect(f.hdot).toBeCloseTo(-1, 9);
    const v: [number, number, number] = [0, 0, 0];
    f.velocityAt(0, 0, 0, v);
    expect(v[1]).toBeCloseTo(0, 12);
    f.velocityAt(0, -0.0199, 0, v);
    expect(v[1]).toBeCloseTo(1, 9); // crown moves up (+y) at 1 m/s
  });

  it('shrouded intake builds positive swirl; decays when the valve closes; capped', () => {
    const spec = testSpec();
    const f = new GasFlowModel(spec);
    let t = 0;
    const dt = 1e-4;
    f.update(snap({ t, intakeLift: 5e-3, intakeMassFlow: 0.03, mass: 3e-4, phase: 'gas-exchange' }), 0);
    for (let i = 0; i < 150; i++) {
      t += dt;
      f.update(snap({ t, intakeLift: 5e-3, intakeMassFlow: 0.03, mass: 3e-4 + 0.03 * t, phase: 'gas-exchange' }), dt);
    }
    expect(f.omega).toBeGreaterThan(0);
    const wMax = (MAX_SWIRL_RATIO * 2 * Math.PI * 600) / 60;
    expect(f.omega).toBeLessThanOrEqual(wMax + 1e-9);
    const w0 = f.omega;
    for (let i = 0; i < 100; i++) {
      t += dt;
      f.update(snap({ t, mass: 7.5e-3 }), dt);
    }
    expect(f.omega).toBeLessThan(w0);
    expect(f.omega).toBeGreaterThan(0);
    // swirl field: v = ω ŷ × r
    const v: [number, number, number] = [0, 0, 0];
    f.velocityAt(TEST_R, -0.01, 0, v);
    expect(v[2]).toBeCloseTo(-f.omega * TEST_R, 9);
    expect(v[0]).toBeCloseTo(0, 12);
  });

  it('a shroud facing the axis gives no swirl; a mirrored shroud gives the opposite sense', () => {
    const run = (dir: number): number => {
      const spec = testSpec();
      spec.intakeValve.shroudDirection = dir;
      const f = new GasFlowModel(spec);
      f.update(snap({ t: 0, intakeLift: 5e-3, intakeMassFlow: 0.02, mass: 5e-4 }), 0);
      for (let i = 1; i <= 50; i++) f.update(snap({ t: i * 1e-4, intakeLift: 5e-3, intakeMassFlow: 0.02, mass: 5e-4 }), 1e-4);
      return f.omega;
    };
    expect(run(0)).toBeCloseTo(0, 9);
    const a = run(Math.PI / 2), b = run(-Math.PI / 2);
    expect(a).toBeGreaterThan(0);
    expect(b).toBeCloseTo(-a, 9);
  });

  it("falls back to u' = ½ S̄p when the snapshot has no turbulence intensity", () => {
    const spec = testSpec();
    const f = new GasFlowModel(spec);
    f.update(snap({ rpm: 600, flame: { turbulenceIntensity: 0 } }), 0);
    expect(f.uPrime).toBeCloseTo(0.5 * (2 * spec.geometry.stroke * 600) / 60, 12);
    f.update(snap({ t: 1e-4, flame: { turbulenceIntensity: 3 } }), 1e-4);
    expect(f.uPrime).toBe(3);
  });
});
