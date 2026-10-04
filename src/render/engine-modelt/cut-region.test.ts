/**
 * ModelTEngineModel.cutRegion: the housings' cut-away quadrant handed to CombustionVisuals.setCutRegion,
 * expressed in each cylinder frame (frames 2 and 4 are mirrored, scale.z = −1).
 *
 * Oracles (independent of the plane formula in model.ts):
 *  - the BUILT section geometry: every vertex of the cut block/head ('block-cut') around the chambers
 *    lies outside the region (accessories built with `cut: false` — generator, starter, brackets — stay
 *    whole but are far from every chamber), and the uncut block ('block-full') has vertices inside it;
 *  - the scene graph: random points in each cylinder frame are mapped to ROOT with the frame's own
 *    matrix (mirror included) and classified against the quadrant {cutSide·x > 0, z > z(block-frame)}.
 */
import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { MODEL_T, MODEL_T_CRUISE } from '../../physics/engines/model-t';
import type { CutPlanes } from '../combustion/index';
import { ModelTEngineModel } from './index';

function inRegion(planes: CutPlanes, p: THREE.Vector3, eps = 0): boolean {
  return planes.every(([nx, ny, nz, d]) => nx * p.x + ny * p.y + nz * p.z > d + eps);
}

/** World (ROOT) positions of every vertex of the meshes under `o`. */
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

describe('ModelTEngineModel.cutRegion', () => {
  const model = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
  model.root.updateMatrixWorld(true);
  const L = model.layout;
  const frames = model.cylinderFrames;
  const toFrame = frames.map((f) => new THREE.Matrix4().copy(f.matrixWorld).invert());

  it('gives one region of two half-spaces per cylinder while the cutaway is on, null when off', () => {
    expect(model.isCutaway).toBe(true);
    for (let i = 0; i < frames.length; i++) {
      const r = model.cutRegion(i);
      expect(r).not.toBeNull();
      expect(r!).toHaveLength(2);
      for (const pl of r!) expect(Math.hypot(pl[0], pl[1], pl[2])).toBeCloseTo(1, 12);
    }
    expect(model.cutRegion(frames.length)).toBeNull();
    model.setCutaway(false);
    for (let i = 0; i < frames.length; i++) expect(model.cutRegion(i)).toBeNull();
    model.setCutaway(true);
    expect(model.cutRegion(0)).not.toBeNull();
  });

  it('matches the built section: cut housings lie outside the region, the full housings reach into it', () => {
    const cut = model.root.getObjectByName('block-cut')!;
    const full = model.root.getObjectByName('block-full')!;
    // block and head around cylinder 1's chamber (8 cm of its axis; deck −12 cm … +10 cm)
    const near = (v: THREE.Vector3): boolean => Math.abs(v.z - L.axisZ[0]) < 0.08 && v.y > L.deckY - 0.12 && v.y < L.deckY + 0.1;
    const vc = vertices(cut).filter(near);
    const vf = vertices(full).filter(near);
    expect(vc.length).toBeGreaterThan(1000);
    for (let i = 0; i < frames.length; i++) {
      const r = model.cutRegion(i)!;
      const p = new THREE.Vector3();
      let inside = 0;
      for (const v of vc) if (inRegion(r, p.copy(v).applyMatrix4(toFrame[i]), 1e-6)) inside++;
      expect(inside, `cylinder ${i + 1}: cut-housing vertices inside the region`).toBe(0);
      let removed = 0;
      for (const v of vf) if (inRegion(r, p.copy(v).applyMatrix4(toFrame[i]), 1e-6)) removed++;
      expect(removed, `cylinder ${i + 1}: full-housing vertices inside the region`).toBeGreaterThan(100);
    }
  });

  it('classifies points like the ROOT quadrant through each (mirrored) frame matrix', () => {
    const zCut = model.root.getObjectByName('block-frame')!.position.z;
    expect(zCut).toBeCloseTo(L.cutZ, 12);
    expect(frames.filter((f) => f.matrixWorld.determinant() < 0)).toHaveLength(2);
    let seed = 12345;
    const rand = (): number => ((seed = (Math.imul(seed, 1103515245) + 12345) >>> 0) / 4294967296) - 0.5;
    const world = new THREE.Vector3();
    for (let i = 0; i < frames.length; i++) {
      const r = model.cutRegion(i)!;
      let inside = 0;
      for (let k = 0; k < 4000; k++) {
        // points around this cylinder's chamber and out past cylinder 1's cut (±0.5 m in z)
        const p = new THREE.Vector3(0.2 * rand(), 0.1 * rand(), rand());
        world.copy(p).applyMatrix4(frames[i].matrixWorld);
        const expected = L.cutSide * world.x > 0 && world.z > zCut;
        if (Math.abs(world.x) < 1e-9 || Math.abs(world.z - zCut) < 1e-9) continue;
        expect(inRegion(r, p), `cylinder ${i + 1} at (${p.x}, ${p.y}, ${p.z})`).toBe(expected);
        if (expected) inside++;
      }
      expect(inside).toBeGreaterThan(50);
    }
  });

  it('cuts through cylinder 1 (the chamber framed by the chamber view), not through the others', () => {
    // Cylinder 1's chamber axis (frame origin, below the roof) on the valve side is inside the region.
    const r1 = model.cutRegion(0)!;
    const valve = L.valves.find((v) => v.cyl === 0 && v.z >= L.cutZ - 1e-12)!;
    const pV = new THREE.Vector3(valve.x, L.deckY + 0.004, valve.z + 0.01).applyMatrix4(toFrame[0]);
    expect(inRegion(r1, pV)).toBe(true);
    // every other cylinder's chamber (bore column at its own axis) stays closed
    for (let i = 1; i < frames.length; i++) {
      const p = new THREE.Vector3(L.cutSide * 0.02, -0.01, 0);
      expect(inRegion(model.cutRegion(i)!, p)).toBe(false);
    }
  });
});
