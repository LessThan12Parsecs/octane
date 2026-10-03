import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { MeshBuilder, box, clipGeometryZ, lathe, planarFace, rect, ringXZ, signedVolume, totalArea } from './geometry';

function maxZ(g: THREE.BufferGeometry): number {
  const P = g.getAttribute('position');
  let m = -Infinity;
  for (let i = 0; i < P.count; i++) m = Math.max(m, P.getZ(i));
  return m;
}

/** Every triangle's geometric normal agrees with its stored vertex normals. */
function orientedConsistently(g: THREE.BufferGeometry): boolean {
  const P = g.getAttribute('position'), N = g.getAttribute('normal');
  const a = new THREE.Vector3(), b = new THREE.Vector3(), c = new THREE.Vector3(), n = new THREE.Vector3();
  for (let i = 0; i < P.count; i += 3) {
    a.fromBufferAttribute(P as THREE.BufferAttribute, i);
    b.fromBufferAttribute(P as THREE.BufferAttribute, i + 1);
    c.fromBufferAttribute(P as THREE.BufferAttribute, i + 2);
    const gn = b.sub(a).cross(c.sub(a));
    n.set(N.getX(i) + N.getX(i + 1) + N.getX(i + 2), N.getY(i) + N.getY(i + 1) + N.getY(i + 2), N.getZ(i) + N.getZ(i + 1) + N.getZ(i + 2));
    if (gn.dot(n) < 0) return false;
  }
  return true;
}

describe('MeshBuilder solids', () => {
  it('box is closed with outward normals', () => {
    const b = new MeshBuilder();
    box(b, [-1, -2, -3], [1, 2, 3]);
    const g = b.build();
    expect(signedVolume(g)).toBeCloseTo(2 * 4 * 6, 9);
    expect(orientedConsistently(g)).toBe(true);
  });

  it('lathe: CCW profile → outward solid, CW cavity → negative volume', () => {
    const b = new MeshBuilder();
    lathe(b, [{ r: 1, y: 0 }, { r: 2, y: 0 }, { r: 2, y: 3 }, { r: 1, y: 3 }], { segments: 256, closed: true });
    const g = b.build();
    expect(signedVolume(g) / (Math.PI * (4 - 1) * 3)).toBeCloseTo(1, 3);
    expect(orientedConsistently(g)).toBe(true);
    const c = new MeshBuilder();
    lathe(c, [{ r: 1, y: 0 }, { r: 1, y: 3 }, { r: 2, y: 3 }, { r: 2, y: 0 }], { segments: 256, closed: true });
    expect(signedVolume(c.build()) / (Math.PI * 3 * 3)).toBeCloseTo(-1, 3);
  });

  it('planar face with matching lathe hole stays watertight (volume of a holed plate)', () => {
    // plate 4×4×1 with a through-hole r=0.5 as a passage
    const b = new MeshBuilder();
    const seg = 64;
    const hole = ringXZ(0, 0, 0.5, seg);
    planarFace(b, rect(-2, -2, 2, 2), [hole], [0, 0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
    planarFace(b, rect(-2, -2, 2, 2), [hole], [0, 1, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
    planarFace(b, rect(0, -2, 1, 2), [], [2, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
    planarFace(b, rect(0, -2, 1, 2), [], [-2, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
    planarFace(b, rect(-2, 0, 2, 1), [], [0, 0, 2], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
    planarFace(b, rect(-2, 0, 2, 1), [], [0, 0, -2], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
    lathe(b, [{ r: 0.5, y: 1 }, { r: 0.5, y: 0 }], { segments: seg }); // passage, normals inward
    const g = b.build();
    const polyArea = 0.5 * seg * 0.25 * Math.sin((2 * Math.PI) / seg);
    expect(signedVolume(g)).toBeCloseTo(16 - polyArea, 6);
    expect(orientedConsistently(g)).toBe(true);
  });
});

describe('clipGeometryZ', () => {
  it('keeps exactly the z ≤ 0 half of a box surface', () => {
    const g = new THREE.BoxGeometry(2, 2, 2, 3, 3, 3);
    const c = clipGeometryZ(g, 0);
    expect(maxZ(c)).toBeLessThanOrEqual(1e-7);
    // back face (4) + half of the four side faces (4 × 2)
    expect(totalArea(c)).toBeCloseTo(12, 6);
    expect(orientedConsistently(c)).toBe(true);
  });
  it('clips at an arbitrary plane and interpolates normals', () => {
    const g = new THREE.SphereGeometry(1, 48, 24);
    const c = clipGeometryZ(g, 0.3);
    expect(maxZ(c)).toBeLessThanOrEqual(0.3 + 1e-6);
    const N = c.getAttribute('normal');
    for (let i = 0; i < N.count; i++) expect(Math.hypot(N.getX(i), N.getY(i), N.getZ(i))).toBeCloseTo(1, 5);
    // spherical cap area removed: 2πRh with h = 0.7
    expect(totalArea(c)).toBeCloseTo(4 * Math.PI - 2 * Math.PI * 0.7, 1);
  });
});
