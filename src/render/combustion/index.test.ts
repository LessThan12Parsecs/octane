/**
 * Smoke test of the three.js wrapper in node (no WebGL: builds the scene graph,
 * runs updates through a synthetic cycle, switches modes, disposes).
 */
import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { cylinderAngleDeg } from '../../physics/core/engine-spec';
import { MODEL_T } from '../../physics/engines/model-t';
import { chamberShapeOf } from './chamber';
import { TRACER_CAPACITY } from './constants';
import { CombustionVisuals, CylinderSnapshotView, temperatureLegend } from './index';
import { engineSnap, flameSnap, modelTSnap, snap, TEST_R } from './test-utils';

function cycle(): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  let t = 0;
  const push = (s: EngineSnapshot) => { out.push(s); t += 2e-5; };
  for (let i = 0; i < 50; i++) push(snap({ t, phase: 'gas-exchange', intakeLift: 5e-3, intakeMassFlow: 0.02, pressure: 0.9e5 }));
  push(snap({ t, spark: { phase: 'charging' } }));
  push(snap({ t, spark: { phase: 'breakdown', breakdownVoltage: 12e3 } }));
  for (let i = 0; i < 10; i++) push(snap({ t, spark: { phase: 'arc', secondaryVoltage: 100, secondaryCurrent: 0.1, energyDelivered: 1e-3 * i } }));
  for (let i = 0; i < 30; i++) {
    push(flameSnap({
      t, flame: { radius: 0.002 + i * 1e-3, area: 1e-4 + i * 1e-4 },
      spark: { phase: 'glow', secondaryVoltage: 450, secondaryCurrent: 0.05 },
    }));
  }
  for (let i = 0; i < 20; i++) push(flameSnap({ t, heatReleaseRate: 2e6, knock: { autoignited: true, oscillation: 1e5 * Math.cos(i) } }));
  push(flameSnap({ t, flame: { stage: 'done' }, massFractionBurned: 1 }));
  return out;
}

describe('CombustionVisuals (scene graph, no GPU)', () => {
  it('builds under a cylinder frame, updates through a cycle, switches modes, disposes', () => {
    const cylinderFrame = new THREE.Object3D();
    const vis = new CombustionVisuals(CFR_F1);
    cylinderFrame.add(vis.root);
    expect(vis.root.parent).toBe(cylinderFrame);
    const names = new Set<string>();
    vis.root.traverse((o) => names.add(o.name));
    expect(names.has('combustion-gas-volume')).toBe(true);
    expect(names.has('combustion-flow-tracers')).toBe(true);
    expect(names.has('combustion-spark')).toBe(true);
    expect(names.has('combustion-chamber-light')).toBe(true);

    let sawFlash = false, sawLight = false, sawKnock = false;
    for (const s of cycle()) {
      vis.update(s, 1 / 60, 1e-3);
      sawFlash ||= vis.state.spark.flash[2] > 0;
      sawLight ||= vis.state.light.intensity > 0;
      sawKnock ||= vis.state.knock.active;
      const vol = vis.root.getObjectByName('combustion-gas-volume') as THREE.Mesh;
      // proxy spans the clearance height (minus the small gap)
      expect(vol.scale.y).toBeGreaterThan(0);
      expect(vol.scale.x).toBeLessThan(TEST_R * 1.01 + 1e-3);
    }
    expect(sawFlash && sawLight && sawKnock).toBe(true);
    expect(vis.tracers.count).toBeGreaterThan(0);

    vis.setMode('temperature');
    vis.update(flameSnap({ t: 1 }), 1 / 60, 1e-3);
    vis.setMode('physical');
    vis.setTemperatureRange(300, 2800);
    vis.setTracersVisible(false);
    vis.setChamberLightEnabled(false);
    vis.update(flameSnap({ t: 1.001 }), 1 / 60, 1e-3);
    const light = vis.root.getObjectByName('combustion-chamber-light') as THREE.PointLight;
    expect(light.intensity).toBe(0);
    vis.dispose();
    expect(vis.root.parent).toBeNull();
  });

  it('exports the legend helper', () => {
    expect(temperatureLegend().stops.length).toBeGreaterThan(8);
  });
});

describe('CombustionVisuals, Model T (four L-head cylinders)', () => {
  const T = chamberShapeOf(MODEL_T);
  const offsets = MODEL_T.layout.firingOffsetDeg;
  /** Engine snapshot at engine angle θ: cylinder i at its local angle, a flame only in `flameCyl`. */
  function engineAt(t: number, theta: number, flameCyl: number, sparks: number[] = [0, 0, 0, 0]): EngineSnapshot {
    const cyl = offsets.map((_, i) => {
      const th = cylinderAngleDeg(MODEL_T, i, theta);
      const x = 0.003 + 0.01 * i; // distinct piston positions
      const lit = i === flameCyl;
      return modelTSnap({
        t, thetaDeg: th, phase: lit ? 'combustion' : 'compression',
        temperatureBurned: lit ? 2300 : 0, temperatureUnburned: 650, massFractionBurned: lit ? 0.2 : 0, heatReleaseRate: lit ? 4e4 : 0,
        flame: lit ? { stage: 'turbulent', radius: 0.03, area: 2e-3, laminarSpeed: 0.4, turbulenceIntensity: 1.5 } : { stage: 'none', radius: 0 },
        spark: { phase: sparks[i] > 0 ? 'glow' : 'off', breakdownCount: sparks[i], breakdownVoltage: 8e3, secondaryVoltage: 300, secondaryCurrent: sparks[i] > 0 ? 0.03 : 0 },
      }, x);
    });
    return { ...engineSnap(cyl), t, thetaDeg: theta };
  }

  it('one instance per cylinder: three L-head proxies each, tracers and light only where asked', () => {
    const vis = offsets.map((_, i) => new CombustionVisuals(MODEL_T, { cylinder: i, tracers: i === 0, light: i === 0 }));
    for (const [i, v] of vis.entries()) {
      const names = new Set<string>();
      v.root.traverse((o) => names.add(o.name));
      expect(v.cylinder).toBe(i);
      expect(names.has('combustion-gas-volume')).toBe(true);
      expect(names.has('combustion-gas-volume-pocket')).toBe(true);
      expect(names.has('combustion-gas-volume-transfer')).toBe(true);
      expect(names.has('combustion-spark')).toBe(true);
      expect(names.has('combustion-flow-tracers')).toBe(i === 0);
      expect(names.has('combustion-chamber-light')).toBe(i === 0);
      expect(v.tracers.capacity).toBe(i === 0 ? TRACER_CAPACITY : 0);
    }
    for (const v of vis) v.dispose();
  });

  it('each instance reads its own cylinder record (local angle, crown, flame, spark), top level as fallback', () => {
    const vis = offsets.map((_, i) => new CombustionVisuals(MODEL_T, { cylinder: i, tracers: false }));
    const s0 = engineAt(0, -10, 2);
    for (const v of vis) v.update(s0, 1 / 60, 1e-3);
    for (const [i, v] of vis.entries()) {
      expect(v.state.flameVisible).toBe(i === 2);
      expect(v.state.h).toBeCloseTo(T.depthTDC + 0.003 + 0.01 * i, 12);
      const vol = v.root.getObjectByName('combustion-gas-volume') as THREE.Mesh;
      expect(vol.scale.y).toBeCloseTo(v.state.h - 4e-4, 12);
    }
    // cylinder 4's trembler buzzes: only instance 3 flashes
    const s1 = engineAt(1e-4, -9.4, 2, [0, 0, 0, 4]);
    for (const v of vis) v.update(s1, 1 / 60, 1e-3);
    expect(vis.map((v) => v.state.breakdownsSeen)).toEqual([0, 0, 0, 4]);
    expect(vis[3].state.spark.flash[2]).toBeGreaterThan(0);
    expect(vis[0].state.spark.flash[2]).toBe(0);
    // no per-cylinder records: every instance falls back to the top level (cylinder 1)
    const flat = { ...engineAt(2e-4, -9, 0), cylinders: undefined };
    for (const v of vis) v.update(flat, 1 / 60, 1e-3);
    for (const v of vis) expect(v.state.flameVisible).toBe(true);
    for (const v of vis) v.dispose();
  });

  it('runs a synthetic cycle with finite outputs; the per-cylinder view allocates nothing', () => {
    const v = new CombustionVisuals(MODEL_T, { cylinder: 1 });
    let t = 0;
    for (let k = 0; k < 200; k++, t += 2e-5) {
      const th = -360 + k * 3.6;
      const s = engineAt(t, th, 1, [0, k > 100 && k < 120 ? k - 100 : 0, 0, 0]);
      v.update(s, 1 / 60, 1e-3);
      const L = v.state.light;
      expect(L.position.every(Number.isFinite)).toBe(true);
      expect(Number.isFinite(L.intensity)).toBe(true);
    }
    expect(v.state.breakdownsSeen).toBeGreaterThan(0);
    const view = new CylinderSnapshotView();
    const s = engineAt(0, 0, 0);
    const a = view.select(s, 2), b = view.select(s, 3);
    expect(a).toBe(b); // the same reused object
    expect(b.thetaDeg).toBe(s.cylinders![3].thetaDeg);
    expect(b.flame).toBe(s.cylinders![3].flame); // by reference, not copied
    expect(view.select(s, 7)).toBe(s);
    v.setCutRegion([[-1, 0, 0, 0], [0, 0, 1, 0]]);
    v.setCutRegion(null);
    v.dispose();
  });
});
