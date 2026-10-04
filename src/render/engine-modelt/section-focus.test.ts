/**
 * The Model T's cutaway section follows the focus cylinder (ModelTEngineModel.setSectionCylinder).
 *
 * Regression for the review finding "chamber close-up and the featured visuals (tracers + light) of
 * cylinders 2–4 sit inside closed metal: the cut exists only at cylinder 1". The reviewer measured, from
 * cameraView('chamber', k) to cylinder k's bore axis under the head roof, 0 / 6 / 4 / 6 hits on the cut
 * housings ('block-cut') for k = 1…4. With the section on cylinder k those sight lines must be clear, the
 * whole view must be cylinder 1's view moved to cylinder k (identical hits on every sight line, the plug
 * gap and the roof point included — both graze the intact pocket roof behind the section plane exactly
 * as they always did for cylinder 1), and the region handed to the in-cylinder visuals (cutRegion) must
 * match the rebuilt geometry.
 *
 * Oracles: ray casts against the BUILT cut housings (front faces of the clipped solids and the back-face
 * section caps), the built vertices against cutRegion (independent of the plane formula), and the scene
 * graph (cylinder frame matrices, mirrored frames included).
 */
import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { MODEL_T, MODEL_T_CRUISE } from '../../physics/engines/model-t';
import type { CutPlanes } from '../combustion/index';
import { ModelTEngineModel } from './index';

function inRegion(planes: CutPlanes, p: THREE.Vector3, eps = 0): boolean {
  return planes.every(([nx, ny, nz, d]) => nx * p.x + ny * p.y + nz * p.z > d + eps);
}

/** World positions of every 3rd vertex of the meshes under `o`. */
function vertices(o: THREE.Object3D): THREE.Vector3[] {
  const out: THREE.Vector3[] = [];
  o.traverse((c) => {
    const m = c as THREE.Mesh;
    if (!m.isMesh) return;
    const P = m.geometry.getAttribute('position');
    for (let i = 0; i < P.count; i += 3) out.push(new THREE.Vector3().fromBufferAttribute(P as THREE.BufferAttribute, i).applyMatrix4(m.matrixWorld));
  });
  return out;
}

/** Hits of the segment a → b on the meshes under `o` (material side respected: section caps are back faces). */
function segmentHits(o: THREE.Object3D, a: THREE.Vector3, b: THREE.Vector3): THREE.Intersection[] {
  const dir = b.clone().sub(a);
  const far = dir.length();
  const ray = new THREE.Raycaster(a, dir.normalize(), 0, far - 1e-4);
  return ray.intersectObject(o, true);
}

describe('ModelTEngineModel section follows the focus cylinder', () => {
  const model = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
  const L = model.layout;
  const n = L.nCyl;

  /** Gas points of cylinder k's chamber the close-up must see unobstructed (world). */
  function openPoints(k: number): THREE.Vector3[] {
    const front = L.valves.find((v) => v.cyl === k && v.z === L.sectionZ[k])!;
    return [
      new THREE.Vector3(0, -0.01, 0).applyMatrix4(model.cylinderFrames[k].matrixWorld), // bore axis 10 mm under the roof
      new THREE.Vector3(front.x, L.deckY + 0.004, front.z - 0.005), // pocket over the front valve, behind the section
      model.cameraView('chamber', k).target,
    ];
  }

  /** Every sight line of the close-up: the open points, the bore axis AT the roof and the plug gap (the light). */
  function sightPoints(k: number): THREE.Vector3[] {
    const f = model.cylinderFrames[k];
    const gap = MODEL_T.sparkPlug.gapCenter;
    return [
      ...openPoints(k),
      new THREE.Vector3(0, 0, 0).applyMatrix4(f.matrixWorld),
      new THREE.Vector3(gap[0], gap[1], gap[2]).applyMatrix4(f.matrixWorld),
    ];
  }

  /** "name@distance" of every cut-housing hit on the close-up's sight lines of cylinder k. */
  function sightHits(k: number): string[][] {
    const cut = model.root.getObjectByName('block-cut')!;
    const cam = model.cameraView('chamber', k).position;
    return sightPoints(k).map((p) => segmentHits(cut, cam, p).map((h) => `${h.object.name}@${h.distance.toFixed(6)}`));
  }
  model.root.updateMatrixWorld(true);
  const reference = sightHits(0); // cylinder 1, the original section

  it('defaults to cylinder 1, exactly the original section (layout cutZ)', () => {
    expect(model.sectionCylinder).toBe(0);
    expect(L.sectionZ[0]).toBe(L.cutZ);
    expect(model.root.getObjectByName('block-frame')!.position.z).toBe(L.cutZ);
    for (let k = 0; k < n; k++) expect(L.sectionZ[k]).toBe(Math.max(...L.valves.filter((v) => v.cyl === k).map((v) => v.z)));
  });

  for (let k = 0; k < n; k++) {
    it(`cylinder ${k + 1}: the close-up looks into the opened chamber (no cut-housing hits) and cutRegion matches the rebuilt section`, () => {
      model.setSectionCylinder(k);
      model.root.updateMatrixWorld(true);
      expect(model.sectionCylinder).toBe(k);
      const frameGroup = model.root.getObjectByName('block-frame')!;
      expect(frameGroup.position.z).toBeCloseTo(L.sectionZ[k], 12);
      expect(model.root.getObjectsByProperty('name', 'block-frame')).toHaveLength(1); // the old set was swapped out

      // sight lines of the chamber close-up (the reviewer's measure: 0 / 6 / 4 / 6 hits before)
      const cut = model.root.getObjectByName('block-cut')!;
      const cam = model.cameraView('chamber', k).position;
      for (const p of openPoints(k)) {
        const hits = segmentHits(cut, cam, p);
        expect(hits.map((h) => h.object.name), `cylinder ${k + 1} → (${p.x.toFixed(3)}, ${p.y.toFixed(3)}, ${p.z.toFixed(3)})`).toEqual([]);
      }
      // ... and the whole view is cylinder 1's view moved here
      expect(sightHits(k)).toEqual(reference);
      expect(cam.z - model.cameraView('chamber', 0).position.z).toBeCloseTo(L.sectionZ[k] - L.sectionZ[0], 12);

      // the region handed to the gas: every cut-housing vertex near the opened chamber lies outside it
      // (in every cylinder frame), and the uncut housings reach into it
      const full = model.root.getObjectByName('block-full')!;
      const near = (v: THREE.Vector3): boolean => Math.abs(v.z - L.axisZ[k]) < 0.08 && v.y > L.deckY - 0.12 && v.y < L.deckY + 0.1;
      const vc = vertices(cut).filter(near);
      const vf = vertices(full).filter(near);
      expect(vc.length).toBeGreaterThan(1000);
      const p = new THREE.Vector3();
      for (let i = 0; i < n; i++) {
        const r = model.cutRegion(i)!;
        const toFrame = new THREE.Matrix4().copy(model.cylinderFrames[i].matrixWorld).invert();
        let inside = 0;
        for (const v of vc) if (inRegion(r, p.copy(v).applyMatrix4(toFrame), 1e-6)) inside++;
        expect(inside, `cylinder ${i + 1}: cut-housing vertices inside the region`).toBe(0);
        let removed = 0;
        for (const v of vf) if (inRegion(r, p.copy(v).applyMatrix4(toFrame), 1e-6)) removed++;
        expect(removed, `cylinder ${i + 1}: full-housing vertices inside the region`).toBeGreaterThan(100);
      }

      // the section cylinder's front valve region is open, the cylinders behind it stay closed, the ones
      // in front have their whole valve-side chamber open
      const front = L.valves.find((v) => v.cyl === k && v.z === L.sectionZ[k])!;
      const toK = new THREE.Matrix4().copy(model.cylinderFrames[k].matrixWorld).invert();
      expect(inRegion(model.cutRegion(k)!, new THREE.Vector3(front.x, L.deckY + 0.004, front.z + 0.01).applyMatrix4(toK))).toBe(true);
      for (let i = 0; i < n; i++) {
        if (i === k) continue;
        const pocket = new THREE.Vector3(L.cutSide * 0.02, -0.01, 0); // valve side of the bore column, at its axis
        expect(inRegion(model.cutRegion(i)!, pocket), `cylinder ${i + 1}`).toBe(i < k);
      }
    });
  }

  it('cutaway off: the section move is deferred (no rebuild) and applied when the cutaway comes back', () => {
    model.setSectionCylinder(0);
    model.setCutaway(false);
    const before = model.root.getObjectByName('block-frame')!;
    model.setSectionCylinder(2);
    expect(model.cutRegion(2)).toBeNull();
    expect(model.root.getObjectByName('block-frame')).toBe(before); // nothing rebuilt while hidden
    expect(model.root.getObjectByName('block-full')!.visible).toBe(true);
    model.setCutaway(true);
    const after = model.root.getObjectByName('block-frame')!;
    expect(after).not.toBe(before);
    expect(after.position.z).toBeCloseTo(L.sectionZ[2], 12);
    expect(model.root.getObjectByName('block-cut')!.visible).toBe(true);
    expect(model.root.getObjectByName('block-full')!.visible).toBe(false);
    expect(model.cutRegion(2)![1][3]).toBeCloseTo(L.sectionZ[2] - L.axisZ[2], 12);
  });

  it('clamps the index, ignores repeats, and disposes the geometry of a replaced section', () => {
    model.setSectionCylinder(0);
    const old = model.root.getObjectByName('block-frame')!;
    const geoms = new Set<THREE.BufferGeometry>();
    old.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.isMesh) geoms.add(m.geometry);
    });
    let disposed = 0;
    for (const g of geoms) g.addEventListener('dispose', () => disposed++);
    model.setSectionCylinder(0); // repeat: no rebuild
    expect(disposed).toBe(0);
    model.setSectionCylinder(99);
    expect(model.sectionCylinder).toBe(n - 1);
    expect(disposed).toBe(geoms.size);
    expect(old.parent).toBeNull();
    model.setSectionCylinder(-3);
    expect(model.sectionCylinder).toBe(0);
    model.setSectionCylinder(NaN);
    expect(model.sectionCylinder).toBe(0);
  });

  it('every cut variant keeps the shadow flags (section caps cast none)', () => {
    model.setSectionCylinder(1);
    model.root.getObjectByName('block-frame')!.traverse((o) => {
      if (!(o as THREE.Mesh).isMesh) return;
      const cap = o.name.includes('-section-');
      expect(o.castShadow, o.name).toBe(!cap);
      expect(o.receiveShadow, o.name).toBe(!cap);
    });
  });
});
