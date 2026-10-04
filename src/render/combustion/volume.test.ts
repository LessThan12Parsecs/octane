/**
 * Gas-volume proxies (no GPU): the flat disc keeps its single back-face cylinder; the L-head gets the
 * bore cylinder, the closed pocket prism and the transfer arc, sharing one uniform set.
 */
import * as THREE from 'three';
import { describe, expect, it } from 'vitest';
import { MODEL_T } from '../../physics/engines/model-t';
import { boundingSize, chamberShapeOf } from './chamber';
import { CombustionVisualState } from './state';
import { modelTSnap, snap, testSpec } from './test-utils';
import { GasVolume, pocketPrismGeometry } from './volume';

const T = chamberShapeOf(MODEL_T);

describe('L-head pocket prism proxy', () => {
  it('is closed and wound outward: signed volume = inset rounded-rectangle area × height', () => {
    const g = 2e-4, y0 = T.deckY + g, y1 = T.roofY - g;
    const geo = pocketPrismGeometry(T, g, y0, y1);
    const p = geo.getAttribute('position');
    let V = 0;
    const a = new THREE.Vector3(), b = new THREE.Vector3(), c = new THREE.Vector3();
    for (let i = 0; i < p.count; i += 3) {
      a.fromBufferAttribute(p, i);
      b.fromBufferAttribute(p, i + 1);
      c.fromBufferAttribute(p, i + 2);
      V += a.dot(b.clone().cross(c)) / 6;
    }
    const W = T.ix1 - T.ix0 + 2 * T.rc - 2 * g, H = T.iz1 - T.iz0 + 2 * T.rc - 2 * g, r = T.rc - g;
    const area = W * H - (4 - Math.PI) * r * r;
    // polygonal corners (16 segments per quarter) cut the arcs slightly
    expect(V).toBeGreaterThan(0);
    expect(Math.abs(V - area * (y1 - y0)) / (area * (y1 - y0))).toBeLessThan(2e-3);
    geo.dispose();
  });
});

describe('GasVolume proxies and uniforms', () => {
  it('flat disc: one back-face cylinder, disc shader path', () => {
    const root = new THREE.Group();
    const gv = new GasVolume(testSpec(), root);
    expect(gv.meshes.length).toBe(1);
    expect(gv.meshes[0]).toBe(gv.mesh);
    expect(gv.mesh.material.side).toBe(THREE.BackSide);
    expect(gv.mesh.material.uniforms.uShape.value).toBe(0);
    expect((gv.mesh.material.uniforms.uValveSeat0.value as THREE.Vector2).toArray()).toEqual([0, -1]);
    gv.dispose();
  });

  it('L-head: bore (back faces), pocket (back faces), transfer arc (front faces) with one shared uniform set', () => {
    const root = new THREE.Group();
    const gv = new GasVolume(MODEL_T, root);
    expect(gv.meshes.map((m) => m.name)).toEqual([
      'combustion-gas-volume', 'combustion-gas-volume-pocket', 'combustion-gas-volume-transfer',
    ]);
    expect(gv.meshes.map((m) => m.material.side)).toEqual([THREE.BackSide, THREE.BackSide, THREE.FrontSide]);
    expect(gv.meshes.map((m) => m.material.uniforms.uPiece.value)).toEqual([0, 1, 2]);
    const u0 = gv.meshes[0].material.uniforms, u1 = gv.meshes[1].material.uniforms, u2 = gv.meshes[2].material.uniforms;
    expect(u1.uFlameR).toBe(u0.uFlameR);
    expect(u2.uKnockAmp).toBe(u0.uKnockAmp);
    expect(u1.uMeshToCyl).not.toBe(u0.uMeshToCyl);
    expect(u0.uShape.value).toBe(1);
    expect((u0.uPocketY.value as THREE.Vector2).toArray()).toEqual([T.deckY, T.roofY]);
    expect((u0.uValveSeat0.value as THREE.Vector2).toArray()).toEqual([MODEL_T.intakeValve.seatY, 1]);
    expect((u0.uValveSeat1.value as THREE.Vector2).toArray()).toEqual([MODEL_T.exhaustValve.seatY, 1]);

    // transfer arc: radius R + gap, spanning the bore-circle arc inside the pocket plan
    const strip = gv.meshes[2];
    strip.updateMatrix();
    const pos = strip.geometry.getAttribute('position');
    const v = new THREE.Vector3();
    let aMin = Infinity, aMax = -Infinity;
    for (let i = 0; i < pos.count; i++) {
      v.fromBufferAttribute(pos, i).applyMatrix4(strip.matrix);
      expect(Math.hypot(v.x, v.z)).toBeCloseTo(T.R + 2e-4, 8); // float32 vertices
      expect(v.y).toBeGreaterThanOrEqual(T.deckY);
      expect(v.y).toBeLessThanOrEqual(T.roofY);
      let a = Math.atan2(v.z, v.x);
      if (a < 0) a += 2 * Math.PI;
      aMin = Math.min(aMin, a);
      aMax = Math.max(aMax, a);
    }
    const [a0, a1] = T.transferArc;
    expect(aMin).toBeLessThan(a0);
    expect(aMax).toBeGreaterThan(a1);
    // the vertex shader's uMeshToCyl is the mesh's own transform
    expect((u2.uMeshToCyl.value as THREE.Matrix4).equals(strip.matrix)).toBe(true);

    // per-frame sync: bore proxy follows the crown, ortho back-off covers the chamber
    const st = new CombustionVisualState(MODEL_T);
    const s = modelTSnap({}, 0.004);
    st.update(s, 0.016, 1e-3);
    gv.sync(st, s);
    expect(gv.mesh.scale.y).toBeCloseTo(st.h - 4e-4, 12);
    expect((u0.uMeshToCyl.value as THREE.Matrix4).equals(gv.mesh.matrix)).toBe(true);
    expect(u0.uBackoff.value).toBeCloseTo(4 * boundingSize(T, st.h), 12);

    gv.setCutRegion([[-1, 0, 0, 0], [0, 0, 1, 0.01]]);
    expect(u0.uCutOn.value).toBe(1);
    expect((u0.uCut1.value as THREE.Vector4).toArray()).toEqual([0, 0, 1, 0.01]);
    gv.setCutRegion(null);
    expect(u0.uCutOn.value).toBe(0);
    gv.setMode('temperature');
    expect(gv.meshes.every((m) => m.material.toneMapped === false)).toBe(true);
    gv.dispose();
  });

  it('flat disc sync is unchanged: proxy height h − 2·gap at −h/2', () => {
    const gv = new GasVolume(testSpec(), new THREE.Group());
    const st = new CombustionVisualState(testSpec());
    const s = snap({ clearanceHeight: 0.0213 });
    st.update(s, 0.016, 1e-3);
    gv.sync(st, s);
    expect(gv.mesh.scale.y).toBeCloseTo(0.0213 - 4e-4, 15);
    expect(gv.mesh.position.y).toBeCloseTo(-0.0213 / 2, 15);
    gv.dispose();
  });
});
