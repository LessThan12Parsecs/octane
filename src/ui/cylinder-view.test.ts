import { describe, expect, it } from 'vitest';
import { cylinderAngleDeg } from '../physics/core/engine-spec';
import { CFR_F1 } from '../physics/engines/cfr';
import { MODEL_T } from '../physics/engines/model-t';
import { createEmptySnapshot } from '../app/sim-client';
import {
  clampCylinder,
  createCylinderViewTarget,
  cylinderCount,
  cylinderView,
  CylinderViewPool,
  firingOrder,
  firingStrip,
} from './cylinder-view';
import { makeSnapshot } from './test-helpers';

function multi(theta = 10) {
  const s = createEmptySnapshot(4);
  s.t = 1.5;
  s.cycle = 3;
  s.thetaDeg = theta;
  s.rpm = 1000;
  s.intakeManifoldPressure = 0.7e5;
  s.gasTorque = 99;
  s.loadTorque = 40;
  s.vehicleSpeed = 11;
  s.firingCylinder = 1;
  s.magnetoEmf = 12;
  s.cylinders!.forEach((c, i) => {
    c.thetaDeg = cylinderAngleDeg(MODEL_T, i, theta);
    c.cycle = 3 - (i === 1 ? 0 : 1);
    c.pressure = (i + 1) * 1e5;
    c.gasTorque = i;
    c.massFractionBurned = i / 10;
  });
  return s;
}

describe('cylinderView', () => {
  it('is the snapshot itself for a single-cylinder stream', () => {
    const s = makeSnapshot({ pressure: 3e5 });
    expect(cylinderView(s, 0, createCylinderViewTarget())).toBe(s);
    expect(cylinderView(s, 2, createCylinderViewTarget())).toBe(s);
  });

  it("maps one cylinder's fields and keeps the engine-level ones, sharing nested objects", () => {
    const s = multi();
    const v = cylinderView(s, 2, createCylinderViewTarget());
    const c = s.cylinders![2];
    expect(v.thetaDeg).toBe(c.thetaDeg);
    expect(v.cycle).toBe(c.cycle);
    expect(v.pressure).toBe(3e5);
    expect(v.gasTorque).toBe(2);
    expect(v.massFractionBurned).toBeCloseTo(0.2, 12);
    expect(v.flame).toBe(c.flame);
    expect(v.spark).toBe(c.spark);
    expect(v.knock).toBe(c.knock);
    expect(v.burnedComposition).toBe(c.burnedComposition);
    // engine level
    expect(v.t).toBe(1.5);
    expect(v.rpm).toBe(1000);
    expect(v.intakeManifoldPressure).toBe(0.7e5);
    expect(v.loadTorque).toBe(40);
    expect(v.vehicleSpeed).toBe(11);
    expect(v.firingCylinder).toBe(1);
    expect(v.magnetoEmf).toBe(12);
    expect(v.cylinders).toBe(s.cylinders);
  });

  it('out-of-range indices fall back to the snapshot', () => {
    const s = multi();
    expect(cylinderView(s, 7, createCylinderViewTarget())).toBe(s);
    expect(cylinderView(s, -1, createCylinderViewTarget())).toBe(s);
  });

  it('the pool maps arrays without allocating new views', () => {
    const pool = new CylinderViewPool();
    const snaps = [multi(10), multi(11), multi(12)];
    const a = pool.map(snaps, 1);
    const first = a[0];
    expect(a.map((v) => v.pressure)).toEqual([2e5, 2e5, 2e5]);
    const b = pool.map(snaps, 3);
    expect(b).toBe(a);
    expect(b[0]).toBe(first);
    expect(b.map((v) => v.thetaDeg)).toEqual(snaps.map((s) => s.cylinders![3].thetaDeg));
    expect(pool.map([], 0)).toHaveLength(0);
  });
});

describe('firing order strip', () => {
  it('Model T fires 1-2-4-3; the CFR is one cylinder', () => {
    expect(firingOrder(MODEL_T)).toEqual([1, 2, 4, 3]);
    expect(firingOrder(CFR_F1)).toEqual([1]);
    expect(cylinderCount(MODEL_T)).toBe(4);
    expect(cylinderCount(CFR_F1)).toBe(1);
  });

  it('each cylinder in firing order is 180° behind the previous one (strokes from the firing offsets)', () => {
    // Engine at −360 (cylinder 1 starts its intake stroke): cylinder 4 is at its firing TDC, 3 is
    // compressing (fires next), 1 fires after it, 2 is exhausting.
    const strip = firingStrip(MODEL_T, makeSnapshot({ thetaDeg: -360 }));
    expect(strip.map((e) => e.number)).toEqual([1, 2, 4, 3]);
    expect(strip.map((e) => e.stroke)).toEqual(['intake', 'exhaust', 'power', 'compression']);
    // Local angles step by −180° (mod 720) along the firing order.
    for (let k = 1; k < strip.length; k++) {
      const d = (((strip[k - 1].thetaDeg - strip[k].thetaDeg) % 720) + 720) % 720;
      expect(d).toBe(180);
    }
  });

  it('reads local angles and spark state from cylinders[] when present', () => {
    const s = multi(100);
    s.cylinders![3].spark.phase = 'arc';
    s.cylinders![2].spark.timerClosed = true;
    s.firingCylinder = 3;
    const strip = firingStrip(MODEL_T, s);
    const byNum = new Map(strip.map((e) => [e.number, e]));
    expect(byNum.get(4)!.sparking).toBe(true);
    expect(byNum.get(4)!.timerClosed).toBe(true); // firingCylinder
    expect(byNum.get(3)!.timerClosed).toBe(true); // spark.timerClosed
    expect(byNum.get(1)!.sparking).toBe(false);
    expect(byNum.get(2)!.thetaDeg).toBe(s.cylinders![1].thetaDeg);
  });

  it('reuses its output entries', () => {
    const out = firingStrip(MODEL_T, makeSnapshot({ thetaDeg: 0 }));
    const e0 = out[0];
    expect(firingStrip(MODEL_T, makeSnapshot({ thetaDeg: 90 }), out)[0]).toBe(e0);
  });

  it('clampCylinder keeps indices in range', () => {
    expect(clampCylinder(5, 4)).toBe(3);
    expect(clampCylinder(-2, 4)).toBe(0);
    expect(clampCylinder(NaN, 4)).toBe(0);
    expect(clampCylinder(2.7, 4)).toBe(2);
    expect(clampCylinder(3, 1)).toBe(0);
  });
});
