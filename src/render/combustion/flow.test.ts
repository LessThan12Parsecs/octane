import { describe, expect, it } from 'vitest';
import { MODEL_T } from '../../physics/engines/model-t';
import { chamberShapeOf, effectiveDepth } from './chamber';
import { CURTAIN_DISCHARGE_COEFF, INTEGRAL_SCALE_FRACTION, MAX_SWIRL_RATIO } from './constants';
import { curtainJetSpeed, GasFlowModel, MAX_JET_SPEED } from './flow';
import { modelTSnap, sharpLHeadSpec, snap, testSpec, TEST_R } from './test-utils';

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

describe('GasFlowModel, L-head (Model T): uniform dilatation with bore → pocket transfer', () => {
  const T = chamberShapeOf(MODEL_T);
  const R = MODEL_T.geometry.bore / 2;
  const Ap = Math.PI * R * R;
  /** Flow at piston displacement x with the crown rising at 1 m/s. */
  function rising(x: number): GasFlowModel {
    const f = new GasFlowModel(MODEL_T);
    f.update(modelTSnap({ t: 0 }, x + 1e-4), 0);
    f.update(modelTSnap({ t: 1e-4 }, x), 1e-4);
    return f;
  }
  const v: [number, number, number] = [0, 0, 0];

  it('crown depth from the spec, D = A_p ḣ / V, Q = −D·V_pocket', () => {
    const f = rising(0.02);
    expect(f.h).toBeCloseTo(T.depthTDC + 0.02, 12);
    expect(f.hdot).toBeCloseTo(-1, 9);
    const V = Ap * f.h + T.pocketVolume;
    expect(f.dilatation).toBeCloseTo((Ap * -1) / V, 9);
    expect(f.transferFlux).toBeCloseTo(-f.dilatation * T.pocketVolume, 15);
    expect(f.transferFlux).toBeGreaterThan(0); // compression pushes gas into the pocket
    expect(f.integralScale).toBeCloseTo(INTEGRAL_SCALE_FRACTION * effectiveDepth(T, f.h), 12);
  });

  it('walls: crown moves with the piston, the roof over the bore and the pocket floor/roof are fixed', () => {
    for (const x of [0.0005, 0.02, 0.06]) {
      const f = rising(x);
      f.velocityAt(0.01, -f.h, -0.01, v);
      expect(v[1]).toBeCloseTo(1, 9); // crown rises at 1 m/s
      f.velocityAt(0.01, 0, 0.01, v);
      expect(v[1]).toBeCloseTo(0, 12); // head roof over the bore
      for (const y of [T.deckY, T.roofY]) {
        f.velocityAt(-0.07, y, 0.02, v);
        expect(v[1]).toBe(0);
      }
      // pocket end wall: no flow through it
      f.velocityAt(T.s1 * T.axis[0], 0.5 * (T.deckY + T.roofY), 0, v);
      expect(v[0] * T.axis[0] + v[2] * T.axis[1]).toBeCloseTo(0, 12);
    }
  });

  it('continuity: every piece changes volume at D; the window delivers Q, the pocket receives Q', () => {
    for (const x of [0.0005, 0.03]) {
      const f = rising(x);
      const D = f.dilatation, Q = f.transferFlux, h = f.h;
      const yc = -h, ya = Math.max(T.deckY, yc), yb = T.roofY;
      // bore column below the window: A_p (v(y) − v_crown) = D A_p (y − y_c)
      const y1 = yc + 0.5 * (ya - yc);
      f.velocityAt(0, y1, 0, v);
      expect(Ap * (v[1] - 1)).toBeCloseTo(D * Ap * (y1 - yc), 12);
      // window slab: net outflow A_p (v(y_b) − v(y_a)) + Q = D·A_p (y_b − y_a)
      f.velocityAt(0, yb, 0, v);
      const vb = v[1];
      f.velocityAt(0, ya, 0, v);
      expect(Ap * (vb - v[1]) + Q).toBeCloseTo(D * Ap * (yb - ya), 12);
      // flux through the window cross-section at the effective interface = Q
      const sI = T.s1 - T.pocketPlanArea / T.width;
      f.velocityAt(sI * T.axis[0], 0.5 * (ya + yb), 0, v);
      expect(Math.hypot(v[0], v[2])).toBeGreaterThan(0);
      if (x > 0.01) {
        // crown below the deck: the window is the whole pocket height, pocket side matches
        expect((v[0] * T.axis[0] + v[2] * T.axis[1]) * T.width * (yb - ya)).toBeCloseTo(Q, 12);
      }
      f.velocityAt((sI + 1e-9) * T.axis[0] - 1e-7, 0.5 * (T.deckY + T.roofY), 0.04, v);
      // pocket side (outside the bore disc): −D (s_end − s) along the axis
      const s = ((sI + 1e-9) * T.axis[0] - 1e-7) * T.axis[0] + 0.04 * T.axis[1];
      expect(v[0] * T.axis[0] + v[2] * T.axis[1]).toBeCloseTo(-D * (T.s1 - s), 12);
    }
  });

  it('pocket cross-section flux = −D × pocket volume beyond it (sharp-cornered pocket: exact)', () => {
    const spec = sharpLHeadSpec();
    const S = chamberShapeOf(spec);
    const f = new GasFlowModel(spec);
    f.update(modelTSnap({ t: 0 }, 0.0301), 0);
    f.update(modelTSnap({ t: 1e-4 }, 0.03), 1e-4);
    for (const sx of [0.055, 0.07, 0.085]) {
      // integrate the axial velocity over the plane s = sx (beyond the bore circle)
      const n = 40;
      let flux = 0;
      const lh = spec.geometry.lHead!;
      const dz = (lh.pocket.zMax - lh.pocket.zMin) / n, dy = (S.roofY - S.deckY) / n;
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          f.velocityAt(-sx, S.deckY + (j + 0.5) * dy, lh.pocket.zMin + (i + 0.5) * dz, v);
          flux += -v[0] * dz * dy;
        }
      }
      const beyond = (S.s1 - sx) * (lh.pocket.zMax - lh.pocket.zMin) * (S.roofY - S.deckY);
      expect(flux).toBeCloseTo(-f.dilatation * beyond, 12);
    }
  });

  it('disc limit: the same field reduces exactly to y·ḣ/h (no pocket)', () => {
    const f = new GasFlowModel(testSpec());
    f.update(snap({ t: 0, clearanceHeight: 0.02 }), 0);
    f.update(snap({ t: 1e-4, clearanceHeight: 0.0199 }), 1e-4);
    expect(f.dilatation).toBe(0);
    f.velocityAt(0.01, -0.01, 0, v);
    expect(v[1]).toBe((-0.01 * f.hdot) / f.h);
  });
});
