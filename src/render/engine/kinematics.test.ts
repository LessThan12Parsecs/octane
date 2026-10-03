import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { EngineSpec } from '../../physics/core/engine-spec';
import { computeLayout, ROCKER_RATIO } from './layout';
import { computePose, createPose, sliderCrank } from './kinematics';

function withOffset(e: number): EngineSpec {
  return { ...CFR_F1, geometry: { ...CFR_F1.geometry, pinOffset: e } };
}

describe('snapshot → part poses', () => {
  for (const e of [0, 0.0025, -0.0015]) {
    it(`slider-crank consistency (pin offset ${e} m)`, () => {
      const L = computeLayout(withOffset(e));
      const pose = createPose();
      const cr = 8.5;
      const headY = L.headY(cr);
      for (let th = -360; th < 360; th += 3.7) {
        const sc = sliderCrank(L, th);
        computePose(L, headY, { thetaDeg: th, pistonDisplacement: sc.pistonDisplacement, intakeLift: 0, exhaustLift: 0 }, [0.025, 0.025], pose);
        // crank pin on its circle
        expect(Math.hypot(pose.crankPin[0], pose.crankPin[1])).toBeCloseTo(L.crankRadius, 12);
        // piston on its axis (wrist pin at the pin offset)
        expect(pose.wristPin[0]).toBe(e);
        // rod length preserved
        expect(Math.abs(pose.rodSpan - L.rodLength)).toBeLessThan(1e-9);
        // rod rotation maps rod-local (0, L) onto W − P
        const tipX = pose.crankPin[0] - L.rodLength * Math.sin(pose.rodRotation);
        const tipY = pose.crankPin[1] + L.rodLength * Math.cos(pose.rodRotation);
        expect(Math.abs(tipX - pose.wristPin[0])).toBeLessThan(1e-9);
        expect(Math.abs(tipY - pose.wristPin[1])).toBeLessThan(1e-9);
        // rod rotation equals the physics rod angle β (sin β = (r sin φ − e)/L, CCW about +z)
        expect(Math.abs(pose.rodRotation - sc.rodAngle)).toBeLessThan(1e-9);
        // crown sits at −h in the cylinder frame, h = h_TDC + displacement
        expect(pose.crownYCyl).toBeCloseTo(-(L.clearanceAtTdc(cr) + sc.pistonDisplacement), 12);
        expect(pose.crownY).toBeCloseTo(pose.wristPin[1] + L.compressionHeight, 12);
        // cam at half speed, crank at −θ
        expect(pose.camAngle).toBeCloseTo((th * Math.PI) / 360, 12);
        expect(pose.crankAngle).toBeCloseTo((-th * Math.PI) / 180 - L.phiTdc, 12);
      }
    });
  }

  it('TDC and BDC displacements match the stroke', () => {
    const L = computeLayout(CFR_F1);
    expect(sliderCrank(L, 0).pistonDisplacement).toBeCloseTo(0, 12);
    expect(sliderCrank(L, 180).pistonDisplacement).toBeCloseTo(L.stroke, 12);
    expect(sliderCrank(L, -360).pistonDisplacement).toBeCloseTo(0, 12);
  });

  it('head face sits h_TDC = S/(CR−1) − V_crevice/A above the crown at TDC', () => {
    const L = computeLayout(CFR_F1);
    const A = (Math.PI * L.bore ** 2) / 4;
    for (const cr of [4, 5, 7, 8, 10, 14, 18]) {
      expect(L.headY(cr) - L.crownTdcY).toBeCloseTo(L.stroke / (cr - 1) - CFR_F1.geometry.creviceVolume / A, 12);
    }
    // with a pin offset the travel (not the stroke) sets the clearance
    const Lo = computeLayout(withOffset(0.004));
    expect(Lo.pistonTravel).toBeGreaterThan(Lo.stroke);
    expect(sliderCrank(Lo, 180).pistonDisplacement).toBeCloseTo(Lo.pistonTravel, 5);
  });

  it('TDC is exact with a pin offset (θ = 0 maximises the wrist-pin height)', () => {
    const L = computeLayout(withOffset(0.004));
    const x0 = sliderCrank(L, 0).pistonDisplacement;
    expect(x0).toBeCloseTo(0, 12);
    expect(sliderCrank(L, 0.5).pistonDisplacement).toBeGreaterThan(x0);
    expect(sliderCrank(L, -0.5).pistonDisplacement).toBeGreaterThan(x0);
  });

  it('valve train: rocker, pushrod and tappet follow the snapshot lift', () => {
    const L = computeLayout(CFR_F1);
    const pose = createPose();
    const headY = L.headY(7);
    const rb: [number, number] = [0.024, 0.026];
    computePose(L, headY, { thetaDeg: 0, pistonDisplacement: 0, intakeLift: 0, exhaustLift: 0 }, rb, pose);
    const len0 = pose.valves.map((v) => Math.hypot(v.pushrodTop[0] - v.pushrodBottom[0], v.pushrodTop[1] - v.pushrodBottom[1], v.pushrodTop[2] - v.pushrodBottom[2]));
    for (const lift of [0.0005, 0.002, 0.004, CFR_F1.intakeValve.maxLift]) {
      computePose(L, headY, { thetaDeg: 0, pistonDisplacement: 0, intakeLift: lift, exhaustLift: lift * 0.5 }, rb, pose);
      for (let i = 0; i < 2; i++) {
        const v = L.valves[i];
        const vp = pose.valves[i];
        const l = i === 0 ? lift : lift * 0.5;
        expect(vp.lift).toBe(l);
        // valve pad on the rocker descends exactly by the lift
        const padY = headY + v.rockerPivot[1] - v.armValve * Math.sin(vp.rockerAngle);
        expect(padY).toBeCloseTo(headY + v.tipY - l, 12);
        // tappet lift = lift / rocker ratio
        expect(vp.tappetLift).toBeCloseTo(l / ROCKER_RATIO, 12);
        expect(vp.tappetFaceY).toBeCloseTo(L.camY + rb[i] + l / ROCKER_RATIO, 12);
        // pushrod (tappet cup → rocker cup) keeps its length (only rocker-cup scrub, < 0.05 mm)
        const len = Math.hypot(vp.pushrodTop[0] - vp.pushrodBottom[0], vp.pushrodTop[1] - vp.pushrodBottom[1], vp.pushrodTop[2] - vp.pushrodBottom[2]);
        expect(Math.abs(len - len0[i])).toBeLessThan(5e-5);
      }
    }
  });

  it('is allocation-free in steady state (reuses the output object)', () => {
    const L = computeLayout(CFR_F1);
    const pose = createPose();
    const a = pose.valves[0].pushrodTop;
    computePose(L, L.headY(7), { thetaDeg: 12, pistonDisplacement: 0.001, intakeLift: 0.001, exhaustLift: 0 }, [0.025, 0.025], pose);
    expect(pose.valves[0].pushrodTop).toBe(a);
  });

  it('tolerates non-finite inputs', () => {
    const L = computeLayout(CFR_F1);
    const pose = createPose();
    computePose(L, L.headY(7), { thetaDeg: NaN, pistonDisplacement: NaN, intakeLift: NaN, exhaustLift: -1 }, [0.025, 0.025], pose);
    expect(Number.isFinite(pose.crankAngle)).toBe(true);
    expect(pose.valves[0].lift).toBe(0);
    expect(pose.valves[1].lift).toBe(0);
  });
});
