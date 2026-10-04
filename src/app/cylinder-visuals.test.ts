/**
 * CylinderVisuals: one CombustionVisuals per cylinder frame of the REAL engine render models (Model T,
 * CFR), with the flow tracers and the chamber light on the featured cylinder only. Scene graph only (no
 * GPU): frames, PointLight count, tracer capacity, the shared gas-volume uniforms (render mode, cut
 * region) and the conductor's gas port.
 */
import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import type { CombustionVisuals, CutPlanes } from '../render/combustion/index';
import { createEngineModel, ModelTEngineModel } from '../render/index';
import { Conductor, GasFanOut, type SimPort } from './conductor';
import { CylinderVisuals, type CylinderFrameHost } from './cylinder-visuals';
import { createEmptySnapshot } from './sim-client';

function lights(root: THREE.Object3D): THREE.PointLight[] {
  const out: THREE.PointLight[] = [];
  root.traverse((o) => {
    if ((o as THREE.PointLight).isPointLight) out.push(o as THREE.PointLight);
  });
  return out;
}

/** The gas-volume uniforms of a visual (shared by all its proxy materials). */
function uniforms(v: CombustionVisuals): Record<string, THREE.IUniform> {
  return (v.root.getObjectByName('combustion-gas-volume') as THREE.Mesh<THREE.BufferGeometry, THREE.ShaderMaterial>).material.uniforms;
}

function cutOf(v: CombustionVisuals): CutPlanes | null {
  const u = uniforms(v);
  if (u.uCutOn.value !== 1) return null;
  const a = (u.uCut0.value as THREE.Vector4).toArray();
  const b = (u.uCut1.value as THREE.Vector4).toArray();
  return [a as [number, number, number, number], b as [number, number, number, number]];
}

const noSim: SimPort = {
  playbackTime: NaN,
  advance: () => null,
  stepDegrees: () => null,
  recentSnapshots: () => [],
  setOperatingPoint: () => undefined,
  reset: () => undefined,
  onCycle: () => undefined,
};

describe('CylinderVisuals with the Ford Model T model', () => {
  it('features one cylinder: tracers and the one chamber light there, every visual in its own frame', () => {
    const model = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const emitted: THREE.Object3D[] = [];
    const cv = new CylinderVisuals(MODEL_T, model, (r) => emitted.push(r), 1);
    const frames = model.cylinderFrames;
    expect(cv.visuals).toHaveLength(4);
    expect(cv.featured).toBe(1);
    cv.visuals.forEach((v, i) => {
      expect(v.cylinder).toBe(i);
      expect(v.root.parent).toBe(frames[i]);
      expect(v.tracers.capacity > 0).toBe(i === 1);
    });
    const ls = lights(model.root);
    expect(ls).toHaveLength(1);
    expect(cv.visuals[1].root.getObjectById(ls[0].id)).toBe(ls[0]);
    expect(emitted).toEqual(cv.visuals.map((v) => v.root));
    expect(cv.port).toBeInstanceOf(GasFanOut);
    expect((cv.port as GasFanOut).ports).toEqual(cv.visuals);
    cv.dispose();
    for (let i = 0; i < 4; i++) expect(frames[i].children.some((c) => c.name === 'combustion-visuals')).toBe(false);
    model.dispose();
  });

  it('moves tracers and light with the focus, rebuilding only the two affected cylinders and keeping mode and cut', () => {
    const model = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    const emitted: THREE.Object3D[] = [];
    const cv = new CylinderVisuals(MODEL_T, model, (r) => emitted.push(r), 0);
    // The conductor hands over the cut region at construction and on cutaway changes.
    const conductor = new Conductor(MODEL_T, noSim, model, cv.port);
    cv.visuals.forEach((v, i) => expect(cutOf(v)).toEqual(model.cutRegion(i)));
    conductor.handleView({ cutaway: true, mode: 'temperature' });
    cv.visuals.forEach((v) => expect(uniforms(v).uMode.value).toBe(1));

    const before = cv.visuals.slice();
    expect(cv.setFeatured(2, conductor.view.mode)).toBe(true);
    expect(cv.featured).toBe(2);
    const after = cv.visuals;
    expect(after[1]).toBe(before[1]);
    expect(after[3]).toBe(before[3]);
    expect(after[0]).not.toBe(before[0]);
    expect(after[2]).not.toBe(before[2]);
    expect(before[0].root.parent).toBeNull();
    expect(before[2].root.parent).toBeNull();
    expect(after[0].root.parent).toBe(model.cylinderFrames[0]);
    expect(after[2].root.parent).toBe(model.cylinderFrames[2]);
    after.forEach((v, i) => {
      expect(v.cylinder).toBe(i);
      expect(v.tracers.capacity > 0).toBe(i === 2);
      expect(uniforms(v).uMode.value).toBe(1);
      expect(cutOf(v)).toEqual(model.cutRegion(i));
    });
    const ls = lights(model.root);
    expect(ls).toHaveLength(1); // constant light count: no material recompiles
    expect(after[2].root.getObjectById(ls[0].id)).toBe(ls[0]);
    expect(emitted.slice(4)).toEqual([after[0].root, after[2].root]);
    expect((cv.port as GasFanOut).ports).toEqual(after);
    expect(cv.setFeatured(2, 'temperature')).toBe(false);

    // the fan-out drives the rebuilt visuals
    const s = createEmptySnapshot(4);
    s.t = 0.5;
    s.rpm = 1000;
    cv.port.update(s, 1 / 60, 1);
    after.forEach((v) => expect(v.state.t).toBe(0.5));

    // cutaway off: every visual (and a later rebuild) gets no cut region
    conductor.handleView({ cutaway: false, mode: 'physical' });
    cv.visuals.forEach((v) => {
      expect(cutOf(v)).toBeNull();
      expect(uniforms(v).uMode.value).toBe(0);
    });
    cv.setFeatured(3.7, conductor.view.mode); // clamped / floored → 3
    expect(cv.featured).toBe(3);
    cv.visuals.forEach((v) => {
      expect(cutOf(v)).toBeNull();
      expect(uniforms(v).uMode.value).toBe(0);
    });
    expect(lights(model.root)).toHaveLength(1);
    cv.dispose();
    model.dispose();
  });

  it('leaves everything as it was when a replacement cannot be built', () => {
    const model = new ModelTEngineModel(MODEL_T, MODEL_T_CRUISE);
    let failOn: THREE.Object3D | null = null;
    const host: CylinderFrameHost = { cylinderFrames: model.cylinderFrames, cutRegion: (i) => model.cutRegion(i) };
    const cv = new CylinderVisuals(MODEL_T, host, (r) => {
      if (r.parent === failOn) throw new Error('no emitters');
    }, 0);
    const before = cv.visuals.slice();
    failOn = model.cylinderFrames[1];
    expect(() => cv.setFeatured(1, 'physical')).toThrow('no emitters');
    expect(cv.featured).toBe(0);
    expect(cv.visuals).toEqual(before);
    expect(lights(model.root)).toHaveLength(1);
    // no orphan visuals left in either frame
    for (const f of model.cylinderFrames) expect(f.children.filter((c) => c.name === 'combustion-visuals')).toHaveLength(1);
    cv.dispose();
    model.dispose();
  });
});

describe('CylinderVisuals with the CFR F-1 (single cylinder, unchanged)', () => {
  it('builds one visual with tracers and light that is itself the gas port; never sets a cut region', () => {
    const model = createEngineModel(CFR_F1, CFR_RON_CONDITIONS);
    expect(model.cutRegion).toBeUndefined();
    const cv = new CylinderVisuals(CFR_F1, model, () => undefined, 3);
    expect(cv.visuals).toHaveLength(1);
    expect(cv.featured).toBe(0);
    const v = cv.visuals[0];
    expect(cv.port).toBe(v);
    expect(v.tracers.capacity).toBeGreaterThan(0);
    expect(lights(model.root)).toHaveLength(1);
    const conductor = new Conductor(CFR_F1, noSim, model, cv.port);
    conductor.handleView({ cutaway: true, mode: 'physical' });
    conductor.handleView({ cutaway: false, mode: 'physical' });
    expect(uniforms(v).uCutOn.value).toBe(0);
    expect(cv.setFeatured(2, 'physical')).toBe(false);
    expect(cv.visuals[0]).toBe(v);
    cv.dispose();
    model.dispose();
  });
});
