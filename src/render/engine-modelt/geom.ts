/**
 * Small geometry helpers for the Model T parts on top of render/engine/geometry.ts (MeshBuilder,
 * lathe, planarFace, sweepTube, extrude): z-axis lathes, plan/section extrusions from polygons, hulls of
 * circles, inverted (cavity) boxes, closed tube solids along a path, helical gears and flat belts.
 * Build-time only (allocates freely).
 */
import * as THREE from 'three';
import {
  MeshBuilder, type ProfilePoint, type TubeFrame, type V3, cross, extrude, lathe, norm, planarFace, rect, ringPoints, scale, sub, sweepTube,
} from '../engine/geometry';
import type { P2 } from './outline';

/** Lathe about an axis parallel to z through (x, y); profile (r, y = axial z). */
export function latheZ(profile: ProfilePoint[], opts: { segments?: number; closed?: boolean; x?: number; y?: number; phiStart?: number; phiLength?: number } = {}): THREE.BufferGeometry {
  const b = new MeshBuilder();
  lathe(b, profile, { segments: opts.segments ?? 48, closed: opts.closed, phiStart: opts.phiStart, phiLength: opts.phiLength });
  const g = b.build();
  g.rotateX(Math.PI / 2); // axial y → z
  g.translate(opts.x ?? 0, opts.y ?? 0, 0);
  return g;
}

/** Closed ring (tube) along z: rIn..rOut, z0..z1, axis through (x, y). */
export function ringSolidZ(z0: number, z1: number, rIn: number, rOut: number, x = 0, y = 0, seg = 48): THREE.BufferGeometry {
  if (rIn <= 1e-9) return latheZ([{ r: 0, y: z0 }, { r: rOut, y: z0 }, { r: rOut, y: z1 }, { r: 0, y: z1 }], { segments: seg, x, y });
  return latheZ([{ r: rIn, y: z0 }, { r: rOut, y: z0 }, { r: rOut, y: z1 }, { r: rIn, y: z1 }], { segments: seg, closed: true, x, y });
}

/** Closed ring (tube) along y around a vertical axis through (x, z). */
export function ringSolidY(y0: number, y1: number, rIn: number, rOut: number, x = 0, z = 0, seg = 48): THREE.BufferGeometry {
  const b = new MeshBuilder();
  if (rIn <= 1e-9) lathe(b, [{ r: 0, y: y0 }, { r: rOut, y: y0 }, { r: rOut, y: y1 }, { r: 0, y: y1 }], { segments: seg, center: [x, 0, z] });
  else lathe(b, [{ r: rIn, y: y0 }, { r: rOut, y: y0 }, { r: rOut, y: y1 }, { r: rIn, y: y1 }], { segments: seg, closed: true, center: [x, 0, z] });
  return b.build();
}

/** THREE.Shape from (u, v) polygons (holes are re-oriented as needed by three). */
export function polyShape(outline: readonly P2[], holes: readonly (readonly P2[])[] = []): THREE.Shape {
  const s = new THREE.Shape(outline.map(([u, v]) => new THREE.Vector2(u, v)));
  for (const h of holes) s.holes.push(new THREE.Path(h.map(([u, v]) => new THREE.Vector2(u, v))));
  return s;
}

/** Extrude a section in the (x, y) plane along z from z0 to z1. */
export function extrudeZ(shape: THREE.Shape, z0: number, z1: number, bevel = 0, curveSegments = 24): THREE.BufferGeometry {
  const g = extrude(shape, { depth: z1 - z0, bevel, bevelSegments: 1, curveSegments });
  g.translate(0, 0, z0);
  return g;
}

/** Extrude a plan polygon given in world (x, z) along +y from y0 to y1. */
export function extrudePlanY(outline: readonly P2[], holes: readonly (readonly P2[])[], y0: number, y1: number, bevel = 0): THREE.BufferGeometry {
  const flip = (p: readonly P2[]) => p.map(([x, z]) => [x, -z] as P2);
  const g = extrude(polyShape(flip(outline), holes.map(flip)), { depth: y1 - y0, bevel, bevelSegments: 1, curveSegments: 16 });
  g.rotateX(-Math.PI / 2); // (u, v, d) → (u, d, −v)
  g.translate(0, y0, 0);
  return g;
}

/** Convex hull (CCW) of circles sampled with `seg` points each. */
export function hullOfCircles(circles: { x: number; y: number; r: number }[], seg = 48): P2[] {
  const pts: P2[] = [];
  for (const c of circles) for (let k = 0; k < seg; k++) {
    const a = (k * 2 * Math.PI) / seg;
    pts.push([c.x + c.r * Math.cos(a), c.y + c.r * Math.sin(a)]);
  }
  pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cr = (o: P2, a: P2, b: P2) => (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo: P2[] = [];
  for (const q of pts) {
    while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop();
    lo.push(q);
  }
  const up: P2[] = [];
  for (let i = pts.length - 1; i >= 0; i--) {
    const q = pts[i];
    while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], q) <= 0) up.pop();
    up.push(q);
  }
  up.pop();
  lo.pop();
  return lo.concat(up);
}

/** Axis-aligned box as a CAVITY (normals pointing inward). */
export function boxCavity(b: MeshBuilder, min: V3, max: V3): void {
  const [x0, y0, z0] = min, [x1, y1, z1] = max;
  planarFace(b, rect(x0, z0, x1, z1), [], [0, y0, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(b, rect(x0, z0, x1, z1), [], [0, y1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(b, rect(y0, z0, y1, z1), [], [x0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(b, rect(y0, z0, y1, z1), [], [x1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(b, rect(x0, y0, x1, y1), [], [0, 0, z0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(b, rect(x0, y0, x1, y1), [], [0, 0, z1], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
}

/** Frames along a polyline with parallel-transported ring axes (first e1 from `up`). */
export function pathFrames(pts: V3[], radius: number | ((s: number) => number), up: V3 = [0, 1, 0]): TubeFrame[] {
  const n = pts.length;
  const len: number[] = [0];
  for (let i = 1; i < n; i++) len.push(len[i - 1] + Math.hypot(...sub(pts[i], pts[i - 1])));
  const total = len[n - 1] || 1;
  const frames: TubeFrame[] = [];
  let e1: V3 | null = null;
  for (let i = 0; i < n; i++) {
    const t = norm(sub(pts[Math.min(i + 1, n - 1)], pts[Math.max(i - 1, 0)]));
    if (!e1) {
      let p = cross(t, up);
      if (Math.hypot(...p) < 1e-6) p = cross(t, [1, 0, 0]);
      e1 = norm(cross(p, t));
    } else {
      const d = e1[0] * t[0] + e1[1] * t[1] + e1[2] * t[2];
      e1 = norm(sub(e1, scale(t, d)));
    }
    const r = typeof radius === 'number' ? radius : radius(len[i] / total);
    frames.push({ c: pts[i], e1, e2: cross(t, e1), r });
  }
  return frames;
}

/** Closed tube solid along frames: inner passage + outer skin + end annuli. */
export function tubeSolid(outer: MeshBuilder, inner: MeshBuilder, frames: TubeFrame[], wall: number, seg: number): void {
  sweepTube(inner, frames, seg, true);
  const of = frames.map((f) => ({ ...f, r: f.r + wall }));
  sweepTube(outer, of, seg, false);
  for (const [i, sgn] of [[0, -1], [frames.length - 1, 1]] as [number, number][]) {
    const a = ringPoints(frames[i], seg), b = ringPoints(of[i], seg);
    const f = frames[i];
    const n = scale(norm(cross(f.e1, f.e2)), sgn);
    for (let k = 0; k < seg; k++) {
      const k1 = (k + 1) % seg;
      outer.quad(a[k], a[k1], b[k1], b[k], n, n, n, n);
    }
  }
}

/** Solid rod (closed, capped) along frames. */
export function rodSolid(b: MeshBuilder, frames: TubeFrame[], seg: number): void {
  sweepTube(b, frames, seg, false);
  for (const [i, sgn] of [[0, -1], [frames.length - 1, 1]] as [number, number][]) {
    const f = frames[i];
    const ring = ringPoints(f, seg);
    const n = scale(norm(cross(f.e1, f.e2)), sgn);
    for (let k = 0; k < seg; k++) {
      const a = ring[k], d = ring[(k + 1) % seg];
      b.tri(f.c[0], f.c[1], f.c[2], a[0], a[1], a[2], d[0], d[1], d[2], n, n, n);
    }
  }
}

/** Gear outline (trapezoidal teeth) as a point list, tooth k centred at phase + 2πk/teeth. */
export function gearOutline(teeth: number, pitchR: number, module: number, phase = 0): P2[] {
  const ra = pitchR + module, rd = pitchR - 1.25 * module;
  const pa = (2 * Math.PI) / teeth;
  const out: P2[] = [];
  const pt = (r: number, a: number) => out.push([r * Math.cos(a), r * Math.sin(a)]);
  for (let k = 0; k < teeth; k++) {
    const c = phase + k * pa;
    pt(rd, c - 0.5 * pa);
    pt(rd, c - 0.3 * pa);
    pt(pitchR - 0.2 * module, c - 0.25 * pa);
    pt(ra, c - 0.11 * pa);
    pt(ra, c + 0.11 * pa);
    pt(pitchR - 0.2 * module, c + 0.25 * pa);
    pt(rd, c + 0.3 * pa);
  }
  return out;
}

/**
 * Closed helical ("spiral-cut") gear along z, centred at the origin: the tooth outline twists by
 * width·tan(β)/r_pitch over the face (sign = hand). Optional bore.
 */
export function helicalGear(teeth: number, pitchR: number, module: number, z0: number, z1: number, helix: number, hand: 1 | -1, phase = 0, boreR = 0, layers = 6): THREE.BufferGeometry {
  const b = new MeshBuilder();
  const twist = (hand * (z1 - z0) * Math.tan(helix)) / pitchR;
  const base = gearOutline(teeth, pitchR, module, phase);
  const n = base.length;
  const ring = (t: number): V3[] => {
    const a = twist * t, ca = Math.cos(a), sa = Math.sin(a), z = z0 + (z1 - z0) * t;
    return base.map(([x, y]) => [x * ca - y * sa, x * sa + y * ca, z] as V3);
  };
  let prev = ring(0);
  for (let L = 1; L <= layers; L++) {
    const cur = ring(L / layers);
    for (let k = 0; k < n; k++) {
      const k1 = (k + 1) % n;
      const a = prev[k], bb = prev[k1], c = cur[k1], d = cur[k];
      // outward normal ≈ edge direction × z
      const ex = bb[0] - a[0], ey = bb[1] - a[1];
      const nn = norm([ey, -ex, 0]);
      b.quad(a, bb, c, d, nn, nn, nn, nn);
    }
    prev = cur;
  }
  // end faces (with bore hole)
  for (const [t, nz] of [[0, -1], [1, 1]] as [number, number][]) {
    const pts = ring(t);
    const outline: P2[] = pts.map((p) => [p[0], p[1]]);
    const holes: P2[][] = [];
    if (boreR > 0) {
      // same vertex placement as the bore lathe below (x = r cos φ, y = −r sin φ after the rotation)
      const h: P2[] = [];
      for (let k = 0; k < 32; k++) {
        const aa = (k * 2 * Math.PI) / 32;
        h.push([boreR * Math.cos(aa), -boreR * Math.sin(aa)]);
      }
      holes.push(h);
    }
    planarFace(b, outline, holes, [0, 0, pts[0][2]], [1, 0, 0], [0, 1, 0], [0, 0, nz]);
  }
  if (boreR > 0) {
    const bb = new MeshBuilder();
    lathe(bb, [{ r: boreR, y: z1 }, { r: boreR, y: z0 }], { segments: 32 });
    const g = bb.build();
    g.rotateX(Math.PI / 2);
    b.append(g);
    g.dispose();
  }
  return b.build();
}

/** Flat belt (rectangular section) around a closed path in the x-y plane at z = zc. */
export function flatBelt(b: MeshBuilder, path: [number, number][], zc: number, width: number, thick: number): void {
  const n = path.length;
  const sec: V3[][] = [];
  const nor: V3[] = [];
  for (let i = 0; i < n; i++) {
    const p = path[i], q = path[(i + 1) % n], o = path[(i - 1 + n) % n];
    let tx = q[0] - o[0], ty = q[1] - o[1];
    const tl = Math.hypot(tx, ty) || 1;
    tx /= tl;
    ty /= tl;
    const nx = ty, ny = -tx; // right of travel = outward for a CCW loop
    sec.push([
      [p[0] + nx * thick, p[1] + ny * thick, zc - width / 2],
      [p[0] + nx * thick, p[1] + ny * thick, zc + width / 2],
      [p[0], p[1], zc + width / 2],
      [p[0], p[1], zc - width / 2],
    ]);
    nor.push([nx, ny, 0]);
  }
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    const faces: [number, V3, V3][] = [[0, nor[i], nor[j]], [1, [0, 0, 1], [0, 0, 1]], [2, scale(nor[i], -1), scale(nor[j], -1)], [3, [0, 0, -1], [0, 0, -1]]];
    for (const [f, ni, nj] of faces) {
      const f1 = (f + 1) % 4;
      b.quad(sec[i][f], sec[i][f1], sec[j][f1], sec[j][f], ni, ni, nj, nj);
    }
  }
}

/** Open-belt path around two circles (CCW), tangent lines included. */
export function beltPath(c1: P2, r1: number, c2: P2, r2: number, n1 = 64, n2 = 48): [number, number][] {
  const dx = c2[0] - c1[0], dy = c2[1] - c1[1];
  const d = Math.hypot(dx, dy);
  const phi = Math.atan2(dy, dx);
  const beta = Math.asin((r1 - r2) / d);
  const aP = phi + Math.PI / 2 - beta;
  const aM = phi - Math.PI / 2 + beta;
  const path: [number, number][] = [];
  let a0 = aP, a1 = aM;
  while (a1 < a0) a1 += 2 * Math.PI;
  for (let i = 0; i <= n1; i++) {
    const a = a0 + ((a1 - a0) * i) / n1;
    path.push([c1[0] + r1 * Math.cos(a), c1[1] + r1 * Math.sin(a)]);
  }
  a0 = aM;
  a1 = aP;
  while (a1 < a0) a1 += 2 * Math.PI;
  for (let i = 0; i <= n2; i++) {
    const a = a0 + ((a1 - a0) * i) / n2;
    path.push([c2[0] + r2 * Math.cos(a), c2[1] + r2 * Math.sin(a)]);
  }
  return path;
}
