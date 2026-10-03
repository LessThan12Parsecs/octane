import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { EngineModel, recommendedCameraView } from './index';
import { sliderCrank } from './kinematics';
import { lobeHeightAt, modelValveLift } from './cam-profile';

function snap(L: EngineModel['layout'], thetaDeg: number, cr: number): EngineSnapshot {
  const sc = sliderCrank(L, thetaDeg);
  return {
    thetaDeg,
    pistonDisplacement: sc.pistonDisplacement,
    rodAngle: sc.rodAngle,
    clearanceHeight: L.clearanceAtTdc(cr) + sc.pistonDisplacement,
    intakeLift: modelValveLift(CFR_F1.intakeValve, thetaDeg),
    exhaustLift: modelValveLift(CFR_F1.exhaustValve, thetaDeg),
  } as unknown as EngineSnapshot;
}

function worldBox(o: THREE.Object3D, frame: THREE.Object3D): THREE.Box3 {
  let top = o;
  while (top.parent) top = top.parent;
  top.updateMatrixWorld(true);
  const inv = new THREE.Matrix4().copy(frame.matrixWorld).invert();
  const box = new THREE.Box3();
  const v = new THREE.Vector3();
  o.traverse((c) => {
    const m = c as THREE.Mesh;
    if (!m.isMesh) return;
    const P = m.geometry.getAttribute('position');
    for (let i = 0; i < P.count; i++) {
      v.fromBufferAttribute(P as THREE.BufferAttribute, i).applyMatrix4(m.matrixWorld).applyMatrix4(inv);
      box.expandByPoint(v);
    }
  });
  return box;
}

describe('EngineModel', () => {
  const cr = 8;
  const model = new EngineModel(CFR_F1, cr);
  const L = model.layout;

  it('places the cylinder frame at the fire deck for the CR', () => {
    const A = (Math.PI * CFR_F1.geometry.bore ** 2) / 4;
    const h = (c: number) => CFR_F1.geometry.stroke / (c - 1) - CFR_F1.geometry.creviceVolume / A;
    expect(model.cylinderFrame.position.y).toBeCloseTo(L.crownTdcY + h(cr), 12);
    model.setCompressionRatio(12);
    expect(model.cylinderFrame.position.y).toBeCloseTo(L.crownTdcY + h(12), 12);
    model.setCompressionRatio(99);
    expect(model.compressionRatio).toBe(L.crRange[1]);
    model.setCompressionRatio(cr);
  });

  it('piston crown coincides with −clearanceHeight in the cylinder frame', () => {
    for (const th of [-300, -180, -30, 0, 45, 170]) {
      const s = snap(L, th, cr);
      model.update(s);
      const piston = model.root.getObjectByName('piston')!;
      const body = piston.getObjectByName('piston-body')!;
      const b = worldBox(body, model.cylinderFrame);
      expect(b.max.y).toBeCloseTo(-s.clearanceHeight, 6);
      expect(model.pose.crownYCyl).toBeCloseTo(-s.clearanceHeight, 12);
    }
  });

  it('spark-plug electrodes bracket the gap at gapCenter', () => {
    const sp = CFR_F1.sparkPlug;
    const ax = sp.axis;
    const k = [0, 1, 2].reduce((a, b) => (Math.abs(ax[b]) > Math.abs(ax[a]) ? b : a), 0);
    expect(Math.abs(ax[k])).toBeGreaterThan(0.99); // test assumes an axis-aligned plug
    const comp = (v: THREE.Vector3) => (k === 0 ? v.x : k === 1 ? v.y : v.z);
    const sgn = Math.sign(ax[k]);
    const bce = worldBox(model.root.getObjectByName('center-electrode')!, model.cylinderFrame);
    const bge = worldBox(model.root.getObjectByName('ground-electrode')!, model.cylinderFrame);
    // centre-electrode tip (its most inboard point along the axis) at gapCenter − axis·gap/2
    const ceTip = sgn < 0 ? comp(bce.min) : comp(bce.max);
    expect(ceTip).toBeCloseTo(sp.gapCenter[k] - sgn * sp.gap / 2, 6);
    // ground-electrode bridge underside at gapCenter + axis·gap/2 (it is the most inboard part)
    const geFace = sgn < 0 ? comp(bge.min) : comp(bge.max);
    expect(geFace).toBeCloseTo(sp.gapCenter[k] + sgn * (sp.gap / 2 + Math.max(0.0008, 0.5 * sp.groundElectrodeWidth)), 6);
    // electrode centred on the gap in the other two coordinates
    for (const j of [0, 1, 2].filter((j) => j !== k)) {
      const c = (j === 0 ? bce.min.x + bce.max.x : j === 1 ? bce.min.y + bce.max.y : bce.min.z + bce.max.z) / 2;
      expect(c).toBeCloseTo(sp.gapCenter[j], 6);
    }
  });

  it('cam lobes touch the tappets (model curve, then the learned curve)', () => {
    const check = () => {
      let maxErr = 0;
      for (let th = -360; th < 360; th += 5) {
        model.update(snap(L, th, cr));
        for (let i = 0; i < 2; i++) {
          const lobe = model.root.getObjectByName(`${L.valves[i].key}-lobe`) as THREE.Mesh;
          const P = lobe.geometry.getAttribute('position');
          const pts = new Float64Array(P.count * 2);
          for (let k = 0; k < P.count; k++) { pts[2 * k] = P.getX(k); pts[2 * k + 1] = P.getY(k); }
          const h = lobeHeightAt(pts, model.pose.camAngle);
          const face = model.pose.valves[i].tappetFaceY - L.camY;
          maxErr = Math.max(maxErr, Math.abs(h - face));
        }
      }
      return maxErr;
    };
    expect(check()).toBeLessThan(1e-4);
    // feed a dense cycle so the learner completes, then re-check
    for (let th = -360; th < 360; th += 0.5) model.update(snap(L, th + 0.25, cr));
    expect(check()).toBeLessThan(1e-4);
  });

  it('toggles cutaway variants', () => {
    model.setCutaway(false);
    expect(model.root.getObjectByName('crankcase-full')!.visible).toBe(true);
    expect(model.root.getObjectByName('crankcase-cut')!.visible).toBe(false);
    model.setCutaway(true);
    expect(model.root.getObjectByName('cylinder-cut')!.visible).toBe(true);
    expect(model.isCutaway).toBe(true);
  });

  it('recommended views look at the engine', () => {
    const v = recommendedCameraView(CFR_F1, cr, 'chamber');
    expect(v.position.distanceTo(v.target)).toBeGreaterThan(0.1);
    const e = recommendedCameraView(CFR_F1, cr);
    expect(e.position.z).toBeGreaterThan(0);
  });

  it('disposes', () => {
    model.dispose();
    expect(model.root.parent).toBeNull();
  });
});
