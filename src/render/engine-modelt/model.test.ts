import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../../physics/engines/model-t';
import { cylinderAngleDeg, type EngineSpec } from '../../physics/core/engine-spec';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { CylinderSnapshot, EngineSnapshot } from '../../physics/core/snapshot';
import { throttlePlateAngle } from '../../physics/gas-exchange/throttle';
// through the render barrel: it registers the Model T builder with createEngineModel
import { CfrEngineRenderModel, createEngineModel, ModelTEngineModel, type EngineRenderModel } from '../index';
import { sliderCrank } from '../engine/kinematics';
import { lobeHeightAt, modelValveLift } from '../engine/cam-profile';
import { pistonDisplacementAt } from './layout';

/** Bounding box of a subtree's meshes, expressed in `frame`'s local coordinates. */
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

function outlineXY(m: THREE.Mesh): Float64Array {
  const P = m.geometry.getAttribute('position');
  const pts = new Float64Array(P.count * 2);
  for (let k = 0; k < P.count; k++) { pts[2 * k] = P.getX(k); pts[2 * k + 1] = P.getY(k); }
  return pts;
}

interface Case {
  name: string;
  spec: EngineSpec;
  op: OperatingPoint;
  build(): EngineRenderModel;
  /** Snapshot at engine angle θ (exact kinematics); `perCylinder` fills s.cylinders. */
  snap(model: EngineRenderModel, thetaDeg: number, perCylinder: boolean): EngineSnapshot;
  /** Expected clearance height h_i of cylinder i for that snapshot. */
  h(model: EngineRenderModel, s: EngineSnapshot, i: number): number;
  pistonBody(model: EngineRenderModel, i: number): THREE.Object3D;
  /** Max |lobe support − tappet face| over every lobe at the current pose, m. */
  lobeError(model: EngineRenderModel): number;
}

const cfr: Case = {
  name: 'CFR F-1',
  spec: CFR_F1,
  op: CFR_RON_CONDITIONS,
  build: () => createEngineModel(CFR_F1, CFR_RON_CONDITIONS),
  snap(model, th) {
    const L = (model as CfrEngineRenderModel).model.layout;
    const sc = sliderCrank(L, th);
    return {
      thetaDeg: th,
      pistonDisplacement: sc.pistonDisplacement,
      rodAngle: sc.rodAngle,
      clearanceHeight: L.clearanceAtTdc(model.compressionRatio) + sc.pistonDisplacement,
      intakeLift: modelValveLift(CFR_F1.intakeValve, th),
      exhaustLift: modelValveLift(CFR_F1.exhaustValve, th),
    } as unknown as EngineSnapshot;
  },
  h: (_m, s) => s.clearanceHeight,
  pistonBody: (m) => m.root.getObjectByName('piston')!.getObjectByName('piston-body')!,
  lobeError(model) {
    const em = (model as CfrEngineRenderModel).model;
    const L = em.layout;
    let e = 0;
    for (let i = 0; i < 2; i++) {
      const lobe = em.root.getObjectByName(`${L.valves[i].key}-lobe`) as THREE.Mesh;
      const h = lobeHeightAt(outlineXY(lobe), em.pose.camAngle);
      e = Math.max(e, Math.abs(h - (em.pose.valves[i].tappetFaceY - L.camY)));
    }
    return e;
  },
};

const modelT: Case = {
  name: 'Ford Model T',
  spec: MODEL_T,
  op: MODEL_T_CRUISE,
  build: () => createEngineModel(MODEL_T, MODEL_T_CRUISE),
  snap(model, th, perCylinder) {
    const L = (model as ModelTEngineModel).layout;
    const iv = L.valves[0].cam, ev = L.valves[1].cam;
    const cyl: CylinderSnapshot[] = L.axisZ.map((_, i) => {
      const ti = cylinderAngleDeg(MODEL_T, i, th);
      const d = pistonDisplacementAt(L, ti);
      return { index: i, thetaDeg: ti, pistonDisplacement: d, clearanceHeight: L.hTdc + d, intakeLift: iv.valveLift(ti), exhaustLift: ev.valveLift(ti) } as unknown as CylinderSnapshot;
    });
    return {
      thetaDeg: th,
      pistonDisplacement: cyl[0].pistonDisplacement,
      clearanceHeight: cyl[0].clearanceHeight,
      intakeLift: cyl[0].intakeLift,
      exhaustLift: cyl[0].exhaustLift,
      cylinders: perCylinder ? cyl : undefined,
    } as unknown as EngineSnapshot;
  },
  h(model, s, i) {
    const L = (model as ModelTEngineModel).layout;
    return L.hTdc + pistonDisplacementAt(L, cylinderAngleDeg(MODEL_T, i, s.thetaDeg));
  },
  pistonBody: (m, i) => m.root.getObjectByName(`cyl${i + 1}-piston`)!.getObjectByName('piston-body')!,
  lobeError(model) {
    const mt = model as ModelTEngineModel;
    const cam = mt.root.getObjectByName('camshaft')!;
    let e = 0;
    for (const v of mt.layout.valves) {
      const lobe = cam.getObjectByName(`lobe-${v.index}`) as THREE.Mesh;
      const h = lobeHeightAt(outlineXY(lobe), cam.rotation.z + lobe.rotation.z);
      const tappet = mt.root.getObjectByName(`tappet-${v.index}`)!;
      e = Math.max(e, Math.abs(h - (tappet.position.y - cam.position.y)));
    }
    return e;
  },
};

for (const C of [cfr, modelT]) {
  describe(`EngineRenderModel contract: ${C.name}`, () => {
    const model = C.build();

    it('has one cylinder frame per cylinder', () => {
      expect(model.cylinderFrames).toHaveLength(C.spec.cylinders);
      for (const f of model.cylinderFrames) expect(f.parent).not.toBeNull();
    });

    it('piston crowns sit at −h in every cylinder frame', () => {
      for (const per of C.spec.cylinders > 1 ? [true, false] : [false]) {
        for (const th of [-300, -180, -30, 0, 45, 170, 333]) {
          const s = C.snap(model, th, per);
          model.update(s);
          model.cylinderFrames.forEach((f, i) => {
            const b = worldBox(C.pistonBody(model, i), f);
            expect(b.max.y).toBeCloseTo(-C.h(model, s, i), 6);
          });
        }
      }
    });

    it('every cam lobe touches its tappet over a cycle', () => {
      let err = 0;
      for (let th = -360; th < 360; th += 5) {
        model.update(C.snap(model, th, true));
        err = Math.max(err, C.lobeError(model));
      }
      expect(err).toBeLessThan(1e-4);
    });

    it('spark-plug electrodes bracket the gap in every cylinder frame', () => {
      const sp = C.spec.sparkPlug;
      const k = [0, 1, 2].reduce((a, b) => (Math.abs(sp.axis[b]) > Math.abs(sp.axis[a]) ? b : a), 0);
      const sgn = Math.sign(sp.axis[k]);
      const comp = (v: THREE.Vector3) => (k === 0 ? v.x : k === 1 ? v.y : v.z);
      const ces: THREE.Object3D[] = [], ges: THREE.Object3D[] = [];
      model.root.traverse((o) => { if (o.name === 'center-electrode') ces.push(o); else if (o.name === 'ground-electrode') ges.push(o); });
      expect(ces).toHaveLength(C.spec.cylinders);
      model.cylinderFrames.forEach((f, i) => {
        const bce = worldBox(ces[i], f);
        const ceTip = sgn < 0 ? comp(bce.min) : comp(bce.max);
        expect(ceTip).toBeCloseTo(sp.gapCenter[k] - (sgn * sp.gap) / 2, 6);
        const bge = worldBox(ges[i], f);
        expect(sgn < 0 ? comp(bge.min) : comp(bge.max)).toBeCloseTo(sp.gapCenter[k] + sgn * (sp.gap / 2 + Math.max(0.0008, 0.5 * sp.groundElectrodeWidth)), 6);
      });
    });

    it('toggles the cutaway', () => {
      model.setCutaway(false);
      expect(model.isCutaway).toBe(false);
      model.setCutaway(true);
      expect(model.isCutaway).toBe(true);
    });

    it('camera views are finite and look at the engine', () => {
      const box = new THREE.Box3().setFromObject(model.root).expandByScalar(0.05);
      const views = [model.cameraView('engine'), ...model.cylinderFrames.map((_, i) => model.cameraView('chamber', i))];
      for (const v of views) {
        for (const c of [...v.position.toArray(), ...v.target.toArray(), v.fov, v.near, v.far]) expect(Number.isFinite(c)).toBe(true);
        expect(v.position.distanceTo(v.target)).toBeGreaterThan(0.1);
        expect(box.containsPoint(v.target)).toBe(true);
        expect(box.containsPoint(v.position)).toBe(false);
      }
    });

    it('chamberShift follows the compression ratio (0 for a fixed head)', () => {
      const cr0 = model.compressionRatio;
      const [lo, hi] = C.spec.geometry.compressionRatioRange;
      if (hi > lo) expect(model.chamberShift(lo + 1, lo)).toBeLessThan(0);
      else {
        model.setCompressionRatio(cr0 + 2);
        expect(model.compressionRatio).toBe(cr0);
        expect(model.chamberShift(cr0 + 2, cr0)).toBe(0);
      }
    });

    it('disposes every geometry and detaches', () => {
      const geoms = new Set<THREE.BufferGeometry>();
      model.root.traverse((o) => { const m = o as THREE.Mesh; if (m.isMesh) geoms.add(m.geometry); });
      let disposed = 0;
      for (const g of geoms) g.addEventListener('dispose', () => disposed++);
      new THREE.Group().add(model.root);
      model.dispose();
      expect(disposed).toBe(geoms.size);
      expect(model.root.parent).toBeNull();
    });
  });
}

describe('ModelTEngineModel specifics', () => {
  it('is what createEngineModel builds for the Model T spec (registered through the render barrel)', () => {
    const m = createEngineModel(MODEL_T, MODEL_T_CRUISE);
    expect(m).toBeInstanceOf(ModelTEngineModel);
    m.dispose();
  });

  it('mirrored cylinder frames carry the spec valve positions onto the world valves', () => {
    const m = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const L = m.layout;
    m.root.updateMatrixWorld(true);
    for (const v of L.valves) {
      const f = m.cylinderFrames[v.cyl];
      expect(f.matrixWorld.determinant() < 0).toBe(L.mirror[v.cyl]);
      const w = new THREE.Vector3(v.spec.position[0], 0, v.spec.position[1]).applyMatrix4(f.matrixWorld);
      const valve = m.root.getObjectByName(`valve-${v.index}`)!;
      expect(w.x).toBeCloseTo(valve.position.x, 12);
      expect(w.z).toBeCloseTo(valve.position.z, 12);
      // closed valve seat line sits on the spec seat plane (deck), within the countersink depth
      expect(Math.abs(valve.position.y - m.pose.valves[v.index].lift - (f.position.y + v.spec.seatY!))).toBeLessThan(0.001);
    }
    m.dispose();
  });

  it('repeated parts share one geometry (clones)', () => {
    const m = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const g = (name: string) => ((m.root.getObjectByName(name)!.getObjectByName('piston-body') ?? m.root.getObjectByName(name)!.children[0]) as THREE.Mesh).geometry;
    expect(g('cyl2-piston')).toBe(g('cyl1-piston'));
    expect(g('cyl4-piston')).toBe(g('cyl1-piston'));
    const lobes: THREE.Mesh[] = [];
    m.root.getObjectByName('camshaft')!.traverse((o) => { if (o.name.startsWith('lobe-')) lobes.push(o as THREE.Mesh); });
    expect(lobes).toHaveLength(8);
    expect(new Set(lobes.map((l) => l.geometry)).size).toBe(1); // identical intake/exhaust lobes
    m.dispose();
  });

  it('without s.cylinders: cylinders 2–4 use the lift curve learned from cylinder 1 (cam model before)', () => {
    const m = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const L = m.layout;
    const cam = L.valves[0].cam;
    const top = (th: number) => modelT.snap(m, th, false);
    m.update(top(100));
    // before learning: the exact cam model
    expect(m.pose.valves[2].lift).toBeCloseTo(cam.valveLift(cylinderAngleDeg(MODEL_T, 1, 100)), 12);
    // feed a cycle whose cylinder-1 intake lift is scaled: the learned curve takes over for the others
    for (let th = -360; th < 360; th += 0.5) {
      const s = top(th + 0.25);
      (s as { intakeLift: number }).intakeLift *= 0.9;
      m.update(s);
    }
    let maxDev = 0;
    for (let th = -360; th < 360; th += 3.7) {
      const s = top(th);
      (s as { intakeLift: number }).intakeLift *= 0.9;
      m.update(s);
      for (const i of [1, 2, 3]) {
        const ti = cylinderAngleDeg(MODEL_T, i, th);
        maxDev = Math.max(maxDev, Math.abs(m.pose.valves[2 * i].lift - 0.9 * cam.valveLift(ti)));
      }
    }
    expect(maxDev).toBeLessThan(2e-4);
    m.dispose();
  });

  it('setControls: timer case turns by half the spark advance (make at θ_i = −advance), butterfly follows the throttle', () => {
    const m = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const L = m.layout;
    const caseG = m.root.getObjectByName('timer-case')!;
    const camG = m.root.getObjectByName('camshaft')!;
    const DEG = Math.PI / 180;
    for (const adv of [-15.5, 25, 64.5]) {
      m.setControls({ sparkAdvanceDeg: adv });
      expect(caseG.rotation.z).toBeCloseTo((-adv * DEG) / 2, 12);
      L.firingOffsetDeg.forEach((off) => {
        m.update(modelT.snap(m, off - adv, true));
        // rotor arm (cam-local +y) reaches the leading edge of this cylinder's segment
        const rotor = Math.PI / 2 + camG.rotation.z;
        const lead = Math.PI / 2 + (off * DEG) / 2 + caseG.rotation.z;
        expect(Math.cos(rotor - lead)).toBeCloseTo(1, 9);
      });
    }
    const thr = m.root.getObjectByName('throttle')!;
    m.setControls({ throttle: 0 });
    expect(thr.rotation.z).toBeCloseTo(throttlePlateAngle(0, L.carb.plateClosedAngle), 12);
    m.setControls({ throttle: 1 });
    expect(thr.rotation.z).toBeCloseTo(Math.PI / 2, 12);
    m.setControls({ throttle: NaN });
    expect(thr.rotation.z).toBeCloseTo(Math.PI / 2, 12);
    m.dispose();
  });

  it('build cost and triangle budget', () => {
    const t0 = performance.now();
    const m = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const ms = performance.now() - t0;
    let tris = 0, meshes = 0, cutTris = 0;
    m.root.traverse((o) => {
      const mm = o as THREE.Mesh;
      if (!mm.isMesh) return;
      meshes++;
      const n = mm.geometry.getAttribute('position').count / 3;
      tris += n;
      let p: THREE.Object3D | null = o, hiddenInCut = false;
      while (p) { if (p.name.endsWith('-full')) hiddenInCut = true; p = p.parent; }
      if (!hiddenInCut) cutTris += n;
    });
    const t1 = performance.now();
    for (let k = 0; k < 2000; k++) m.update(modelT.snap(m, (k % 720) - 360, true));
    const usPerUpdate = ((performance.now() - t1) / 2000) * 1000;
    console.info(`Model T render model: build ${ms.toFixed(0)} ms, ${meshes} meshes, ${tris} triangles (${cutTris} drawn in cutaway), update+snapshot ${usPerUpdate.toFixed(1)} µs`);
    expect(tris).toBeLessThan(600_000);
    expect(ms).toBeLessThan(5000); // generous: the machine is shared
    m.dispose();
  });
});
