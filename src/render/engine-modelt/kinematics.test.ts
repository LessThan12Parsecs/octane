import { describe, expect, it } from 'vitest';
import { MODEL_T } from '../../physics/engines/model-t';
import { cylinderAngleDeg } from '../../physics/core/engine-spec';
import { computeModelTLayout, pistonDisplacementAt } from './layout';
import { computeModelTPose, createModelTInput, createModelTPose } from './kinematics';

const DEG = Math.PI / 180;

describe('Model T pose kinematics', () => {
  const L = computeModelTLayout(MODEL_T);

  it('firing order 1-2-4-3: each cylinder is at firing TDC at its offset, crown at −h_TDC in its frame', () => {
    const inp = createModelTInput(L);
    const P = createModelTPose(L);
    for (let i = 0; i < 4; i++) {
      inp.thetaDeg = MODEL_T.layout.firingOffsetDeg[i];
      computeModelTPose(L, inp, P);
      expect(P.cylinders[i].thetaDeg).toBeCloseTo(0, 9);
      expect(P.cylinders[i].crownYCyl).toBeCloseTo(-L.hTdc, 12);
      // its throw points straight up
      expect(P.cylinders[i].crankPin[0]).toBeCloseTo(0, 12);
      expect(P.cylinders[i].crankPin[1]).toBeCloseTo(L.crankRadius, 12);
    }
  });

  it('derives exact slider-crank kinematics per cylinder (rod span = rod length, rod angle = physics β)', () => {
    const inp = createModelTInput(L);
    const P = createModelTPose(L);
    for (let th = -360; th < 360; th += 13) {
      inp.thetaDeg = th;
      computeModelTPose(L, inp, P);
      for (let i = 0; i < 4; i++) {
        const c = P.cylinders[i];
        const ti = cylinderAngleDeg(MODEL_T, i, th);
        expect(c.pistonDisplacement).toBeCloseTo(pistonDisplacementAt(L, ti), 12);
        expect(c.rodSpan).toBeCloseTo(L.rodLength, 9);
        expect(c.rodRotation).toBeCloseTo(Math.asin((L.crankRadius * Math.sin(ti * DEG)) / L.rodLength), 9);
        // crown in the cylinder frame = −(h_TDC + displacement)
        expect(c.crownYCyl).toBeCloseTo(-(L.hTdc + c.pistonDisplacement), 12);
      }
      // throws 1&4 coincide, 2&3 coincide
      expect(P.cylinders[0].crankPin[1]).toBeCloseTo(P.cylinders[3].crankPin[1], 12);
      expect(P.cylinders[1].crankPin[0]).toBeCloseTo(P.cylinders[2].crankPin[0], 12);
      expect(P.cylinders[0].crankPin[0]).toBeCloseTo(-P.cylinders[1].crankPin[0], 12);
    }
  });

  it('uses snapshot values when given; tappets follow the cam, lash closes while the valve is open', () => {
    const inp = createModelTInput(L);
    const P = createModelTPose(L);
    inp.pistonDisplacement.fill(0.02);
    inp.intakeLift.fill(0.003);
    inp.thetaDeg = 40;
    computeModelTPose(L, inp, P);
    for (const c of P.cylinders) expect(c.pistonDisplacement).toBe(0.02);
    for (const v of L.valves) if (v.kind === 'intake') expect(P.valves[v.index].lift).toBe(0.003);
    // derived (NaN) lifts: valve = max(0, cam − lash); the stem end meets the tappet when open
    const inp2 = createModelTInput(L);
    for (let th = -360; th < 360; th += 2.5) {
      inp2.thetaDeg = th;
      computeModelTPose(L, inp2, P);
      for (const v of L.valves) {
        const vp = P.valves[v.index];
        const stemEnd = L.seatLineY + L.valve.stemEndY + vp.lift;
        const tappetTop = vp.tappetFaceY + L.tappet.length;
        const gap = stemEnd - tappetTop;
        expect(gap).toBeGreaterThan(-1e-12);
        if (vp.lift > 0) expect(gap).toBeCloseTo(0, 12);
        else expect(gap).toBeLessThanOrEqual(L.lash + 1e-12);
        expect(vp.springScale).toBeCloseTo((L.spring.installed - vp.lift) / L.spring.installed, 12);
      }
    }
    expect(P.camAngle).toBeCloseTo(((360 - 2.5) * DEG) / 2, 12);
    expect(P.fanAngle).toBeCloseTo(-(360 - 2.5) * DEG * L.fan.ratio, 9);
  });

  it('is allocation-free in steady state and tolerates non-finite inputs', () => {
    const inp = createModelTInput(L);
    const P = createModelTPose(L);
    const pin = P.cylinders[2].crankPin;
    const v0 = P.valves[0];
    inp.thetaDeg = NaN;
    inp.intakeLift[1] = -1;
    computeModelTPose(L, inp, P);
    expect(P.cylinders[2].crankPin).toBe(pin);
    expect(P.valves[0]).toBe(v0);
    expect(Number.isFinite(P.crankAngle)).toBe(true);
    expect(P.valves[2].lift).toBe(0);
  });
});
