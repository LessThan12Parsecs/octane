/**
 * Smoke test of the three.js wrapper in node (no WebGL: builds the scene graph,
 * runs updates through a synthetic cycle, switches modes, disposes).
 */
import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { CombustionVisuals, temperatureLegend } from './index';
import { flameSnap, snap, TEST_R } from './test-utils';

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
