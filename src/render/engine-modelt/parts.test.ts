import { describe, expect, it } from 'vitest';
import * as THREE from 'three';
import { MODEL_T } from '../../physics/engines/model-t';
import { merge, signedVolume } from '../engine/geometry';
import { computeModelTLayout } from './layout';
import { buildBlockShell } from './parts-block';
import { buildHeadShell, PLUG_SEG } from './parts-head';
import { polygonArea } from './outline';

/**
 * Closed, consistently oriented 2-manifold: after welding positions (1 µm grid) every directed edge
 * appears exactly once and its reverse exactly once. Returns the number of offending edges.
 */
function manifoldDefects(g: THREE.BufferGeometry): number {
  const P = g.getAttribute('position');
  const q = (i: number) => `${Math.round(P.getX(i) * 1e6)},${Math.round(P.getY(i) * 1e6)},${Math.round(P.getZ(i) * 1e6)}`;
  const edges = new Map<string, number>();
  for (let i = 0; i < P.count; i += 3) {
    const k = [q(i), q(i + 1), q(i + 2)];
    for (let e = 0; e < 3; e++) {
      const key = `${k[e]}|${k[(e + 1) % 3]}`;
      edges.set(key, (edges.get(key) ?? 0) + 1);
    }
  }
  let bad = 0;
  for (const [key, n] of edges) {
    const [a, b] = key.split('|');
    if (n !== 1 || edges.get(`${b}|${a}`) !== 1) bad++;
  }
  return bad;
}

/** Area of a regular n-gon of circumradius r (the discretised circles of the shells). */
const ngon = (r: number, n: number) => 0.5 * n * r * r * Math.sin((2 * Math.PI) / n);

describe('Model T static solids', () => {
  const L = computeModelTLayout(MODEL_T);

  it('head shell is a watertight 2-manifold with the analytic volume (slab − chambers − plug holes − jacket)', () => {
    const sh = buildHeadShell(L);
    const g = merge([sh.paint.build(), sh.machined.build(), sh.raw.build()]);
    expect(manifoldDefects(g)).toBe(0);
    const H = L.head;
    const y0 = L.headBottomY, y1 = L.headTopY;
    const Aslab = (H.xL - H.xR) * (H.zF - H.zB);
    let chambers = 0;
    for (const ch of L.chambers) {
      chambers += polygonArea(ch.union) * (L.pocketRoofY - y0) + polygonArea(ch.circle) * (L.cylOriginY - L.pocketRoofY);
    }
    const plugs = L.nCyl * ngon(H.plugPassageRadius, PLUG_SEG) * (y1 - L.pocketRoofY);
    const j = H.jacket;
    const jacket = ((j.x1 - j.x0) * (j.z1 - j.z0) - L.nCyl * ngon(H.plugTubeRadius, PLUG_SEG)) * (j.y1 - j.y0);
    const expected = Aslab * (y1 - y0) - chambers - plugs - jacket;
    expect(signedVolume(g) / expected).toBeCloseTo(1, 6); // float32 vertex positions
  });

  it('the head chamber reproduces the clearance volume of the physics chamber (Vc − crevice, 0.5 %)', () => {
    const lh = MODEL_T.geometry.lHead!;
    const ch = L.chambers[0];
    // at TDC: pocket region from the deck (valve seats) to its roof + the bore column above the crown
    const crownTop = L.deckY + lh.crownAboveDeckAtTDC;
    const vRender = polygonArea(ch.pocketRoof) * (L.pocketRoofY - L.deckY) + polygonArea(ch.circle) * (L.cylOriginY - crownTop);
    const A = (Math.PI * MODEL_T.geometry.bore ** 2) / 4;
    const Vc = (A * MODEL_T.geometry.stroke) / (MODEL_T.geometry.compressionRatio - 1);
    expect(Math.abs(vRender / (Vc - MODEL_T.geometry.creviceVolume) - 1)).toBeLessThan(0.005);
  });

  it('block shell (bores, crankcase, valve chamber, seats, ports, galleries, jacket) is a watertight 2-manifold', () => {
    const sh = buildBlockShell(L);
    const g = merge([sh.paint.build(), sh.machined.build(), sh.raw.build()]);
    expect(manifoldDefects(g)).toBe(0);
    const K = L.block;
    // outer prism minus the box/cylinder cavities; the cored ports/seats/galleries take the rest (< 0.6 L)
    const section = [[K.xR, K.yBot], [K.xLb, K.yBot], [K.xLb, K.yFlareBot], [K.xLu, K.yFlareTop], [K.xLu, K.yTop], [K.xR, K.yTop]] as [number, number][];
    const outer = polygonArea(section) * (K.zF - K.zB);
    const cc = K.cc;
    const ccV = (cc.x1 - cc.x0) * ((cc.zLo1 - cc.zLo0) * cc.yStep + (cc.zUp1 - cc.zUp0) * (cc.yTop - cc.yStep));
    const vc = K.vc;
    const vcV = (vc.x1 - K.xR) * (vc.y1 - vc.y0) * (vc.z1 - vc.z0);
    const j = K.jacket;
    const jV = (j.x1 - j.x0) * (j.y1 - j.y0) * (j.z1 - j.z0);
    const boreV = L.nCyl * ngon(L.boreRadius, 64) * (K.yTop - cc.yTop);
    const rest = outer - ccV - vcV - jV - boreV - signedVolume(g);
    expect(signedVolume(g)).toBeGreaterThan(0);
    expect(rest).toBeGreaterThan(0.1e-3);
    expect(rest).toBeLessThan(0.6e-3);
  });
});
