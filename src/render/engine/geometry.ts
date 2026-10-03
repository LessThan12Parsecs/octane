/**
 * Procedural geometry helpers for the engine model.
 *
 * Cutaway strategy (no CSG, no stencil, no renderer clipping state needed):
 * static housings are modelled as CLOSED solids whose internal passages
 * (bore, ports, water jacket, bearing bores…) are surfaces whose normals point
 * into the passage. For the cutaway each solid is clipped by the plane z = 0
 * (`clipGeometryZ`) and drawn twice: front faces with the PBR material and
 * back faces with a flat "section" material (see materials.ts). Looking into a
 * cut solid, the first surface hit is an interior back face → it reads as a
 * solid section cap; looking into a passage, the passage wall's front face is
 * hit → the passage shows. So every triangle's winding must match its
 * intended normal: `MeshBuilder` orients triangles automatically.
 */
import * as THREE from 'three';
import { mergeGeometries, mergeVertices, toCreasedNormals } from 'three/examples/jsm/utils/BufferGeometryUtils.js';

export type V3 = [number, number, number];

/** Accumulates non-indexed triangles (position, normal, uv). */
export class MeshBuilder {
  readonly pos: number[] = [];
  readonly nor: number[] = [];
  readonly uv: number[] = [];

  /**
   * Add a triangle with per-vertex normals. The winding is flipped if the
   * geometric normal disagrees with the supplied normals.
   */
  tri(
    ax: number, ay: number, az: number,
    bx: number, by: number, bz: number,
    cx: number, cy: number, cz: number,
    na: V3, nb: V3, nc: V3,
    uva: [number, number] = [0, 0], uvb: [number, number] = [0, 0], uvc: [number, number] = [0, 0],
  ): void {
    const ux = bx - ax, uy = by - ay, uz = bz - az;
    const vx = cx - ax, vy = cy - ay, vz = cz - az;
    const gx = uy * vz - uz * vy, gy = uz * vx - ux * vz, gz = ux * vy - uy * vx;
    const area2 = gx * gx + gy * gy + gz * gz;
    if (area2 < 1e-24) return; // degenerate
    const sx = na[0] + nb[0] + nc[0], sy = na[1] + nb[1] + nc[1], sz = na[2] + nb[2] + nc[2];
    const flip = gx * sx + gy * sy + gz * sz < 0;
    if (!flip) {
      this.pos.push(ax, ay, az, bx, by, bz, cx, cy, cz);
      this.nor.push(...na, ...nb, ...nc);
      this.uv.push(...uva, ...uvb, ...uvc);
    } else {
      this.pos.push(ax, ay, az, cx, cy, cz, bx, by, bz);
      this.nor.push(...na, ...nc, ...nb);
      this.uv.push(...uva, ...uvc, ...uvb);
    }
  }

  /** Quad a-b-c-d (any winding) with per-vertex normals. */
  quad(a: V3, b: V3, c: V3, d: V3, na: V3, nb: V3, nc: V3, nd: V3): void {
    this.tri(a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2], na, nb, nc);
    this.tri(a[0], a[1], a[2], c[0], c[1], c[2], d[0], d[1], d[2], na, nc, nd);
  }

  /** Append all triangles of a geometry (orientation taken from its normals). */
  append(g: THREE.BufferGeometry): void {
    const src = g.index ? g.toNonIndexed() : g;
    const P = src.getAttribute('position'), N = src.getAttribute('normal');
    for (let i = 0; i < P.count; i += 3) {
      const n = (k: number): V3 => [N.getX(i + k), N.getY(i + k), N.getZ(i + k)];
      this.tri(
        P.getX(i), P.getY(i), P.getZ(i),
        P.getX(i + 1), P.getY(i + 1), P.getZ(i + 1),
        P.getX(i + 2), P.getY(i + 2), P.getZ(i + 2),
        n(0), n(1), n(2),
      );
    }
    if (src !== g) src.dispose();
  }

  get empty(): boolean {
    return this.pos.length === 0;
  }

  build(): THREE.BufferGeometry {
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.Float32BufferAttribute(this.pos, 3));
    g.setAttribute('normal', new THREE.Float32BufferAttribute(this.nor, 3));
    g.setAttribute('uv', new THREE.Float32BufferAttribute(this.uv, 2));
    g.computeBoundingSphere();
    return g;
  }
}

// ---------------------------------------------------------------------------
// Lathe (surface of revolution about a vertical axis)
// ---------------------------------------------------------------------------

export interface ProfilePoint {
  r: number;
  y: number;
  /** Smooth vertex: normals averaged with the neighbouring segment. Default: sharp. */
  smooth?: boolean;
}

export interface LatheOptions {
  segments?: number;
  center?: V3;
  /** Closed profile loop. Orientation decides normal direction (see below). */
  closed?: boolean;
  phiStart?: number;
  phiLength?: number;
}

/**
 * Revolve an (r, y) profile about the vertical axis through `center`.
 * Normal of segment p_i→p_{i+1} is (dy, −dr) in the (r, y) plane: traverse a
 * solid's cross-section COUNTER-clockwise (r to the right, y up) for outward
 * normals; traverse a cavity CLOCKWISE for normals pointing into the cavity.
 * Point at angle φ: (cx + r cos φ, cy + y, cz + r sin φ).
 */
export function lathe(b: MeshBuilder, profile: ProfilePoint[], opts: LatheOptions = {}): void {
  const seg = opts.segments ?? 64;
  const c = opts.center ?? [0, 0, 0];
  const closed = opts.closed ?? false;
  const phi0 = opts.phiStart ?? 0;
  const dphi = (opts.phiLength ?? Math.PI * 2) / seg;
  const n = profile.length;
  const nSeg = closed ? n : n - 1;
  const segN: [number, number][] = [];
  for (let i = 0; i < nSeg; i++) {
    const p = profile[i], q = profile[(i + 1) % n];
    const dr = q.r - p.r, dy = q.y - p.y;
    const l = Math.hypot(dr, dy) || 1;
    segN.push([dy / l, -dr / l]);
  }
  const vertexNormal = (i: number, s: number): [number, number] => {
    // normal at profile vertex i as used by segment s
    if (!profile[i].smooth) return segN[s];
    const adj: number[] = [];
    if (closed) adj.push((i - 1 + n) % n, i % n);
    else {
      if (i >= 1) adj.push(i - 1);
      if (i <= n - 2) adj.push(i);
    }
    let nr = 0, ny = 0;
    for (const t of adj) { nr += segN[t][0]; ny += segN[t][1]; }
    const l = Math.hypot(nr, ny) || 1;
    return [nr / l, ny / l];
  };
  for (let s = 0; s < nSeg; s++) {
    const i0 = s, i1 = (s + 1) % n;
    const p = profile[i0], q = profile[i1];
    const np = vertexNormal(i0, s), nq = vertexNormal(i1, s);
    for (let k = 0; k < seg; k++) {
      const a0 = phi0 + k * dphi, a1 = a0 + dphi;
      const c0 = Math.cos(a0), s0 = Math.sin(a0), c1 = Math.cos(a1), s1 = Math.sin(a1);
      const A: V3 = [c[0] + p.r * c0, c[1] + p.y, c[2] + p.r * s0];
      const Bv: V3 = [c[0] + p.r * c1, c[1] + p.y, c[2] + p.r * s1];
      const C: V3 = [c[0] + q.r * c1, c[1] + q.y, c[2] + q.r * s1];
      const D: V3 = [c[0] + q.r * c0, c[1] + q.y, c[2] + q.r * s0];
      const nA: V3 = [np[0] * c0, np[1], np[0] * s0];
      const nB: V3 = [np[0] * c1, np[1], np[0] * s1];
      const nC: V3 = [nq[0] * c1, nq[1], nq[0] * s1];
      const nD: V3 = [nq[0] * c0, nq[1], nq[0] * s0];
      b.quad(A, Bv, C, D, nA, nB, nC, nD);
    }
  }
}

/** Circle polygon matching `lathe` vertex placement (for holes in planar faces). */
export function ringXZ(cx: number, cz: number, r: number, seg: number, phi0 = 0): [number, number][] {
  const out: [number, number][] = [];
  for (let k = 0; k < seg; k++) {
    const a = phi0 + (k * 2 * Math.PI) / seg;
    out.push([cx + r * Math.cos(a), cz + r * Math.sin(a)]);
  }
  return out;
}

// ---------------------------------------------------------------------------
// Planar faces with holes
// ---------------------------------------------------------------------------

/**
 * Triangulate a 2-D outline with holes and place it in 3-D:
 * P(u, v) = origin + u·U + v·V. Triangles are oriented to `normal`.
 */
export function planarFace(
  b: MeshBuilder,
  outline: [number, number][],
  holes: [number, number][][],
  origin: V3, U: V3, V: V3, normal: V3,
): void {
  const contour = outline.map(([u, v]) => new THREE.Vector2(u, v));
  const hs = holes.map((h) => h.map(([u, v]) => new THREE.Vector2(u, v)));
  if (THREE.ShapeUtils.isClockWise(contour)) contour.reverse();
  for (const h of hs) if (!THREE.ShapeUtils.isClockWise(h)) h.reverse();
  const faces = THREE.ShapeUtils.triangulateShape(contour, hs);
  const all = contour.concat(...hs);
  const P = (i: number): V3 => {
    const p = all[i];
    return [origin[0] + p.x * U[0] + p.y * V[0], origin[1] + p.x * U[1] + p.y * V[1], origin[2] + p.x * U[2] + p.y * V[2]];
  };
  for (const f of faces) {
    const a = P(f[0]), c = P(f[1]), d = P(f[2]);
    b.tri(a[0], a[1], a[2], c[0], c[1], c[2], d[0], d[1], d[2], normal, normal, normal);
  }
}

/** Rectangle outline helper (CCW). */
export function rect(u0: number, v0: number, u1: number, v1: number): [number, number][] {
  return [[u0, v0], [u1, v0], [u1, v1], [u0, v1]];
}

// ---------------------------------------------------------------------------
// Swept tubes (passages, pipes)
// ---------------------------------------------------------------------------

export interface TubeFrame {
  c: V3;
  /** Unit vectors spanning the ring plane. */
  e1: V3;
  e2: V3;
  r: number;
}

/** Points of a ring: c + r (cos a e1 + sin a e2). */
export function ringPoints(f: TubeFrame, seg: number, phi0 = 0): V3[] {
  const out: V3[] = [];
  for (let k = 0; k < seg; k++) {
    const a = phi0 + (k * 2 * Math.PI) / seg;
    const ca = Math.cos(a), sa = Math.sin(a);
    out.push([
      f.c[0] + f.r * (ca * f.e1[0] + sa * f.e2[0]),
      f.c[1] + f.r * (ca * f.e1[1] + sa * f.e2[1]),
      f.c[2] + f.r * (ca * f.e1[2] + sa * f.e2[2]),
    ]);
  }
  return out;
}

/**
 * Sweep circles through a list of frames. `inward` = normals point toward
 * the tube axis (a passage); otherwise outward (a solid rod / pipe exterior).
 */
export function sweepTube(b: MeshBuilder, frames: TubeFrame[], seg: number, inward: boolean, phi0 = 0): void {
  const rings = frames.map((f) => ringPoints(f, seg, phi0));
  const sgn = inward ? -1 : 1;
  const normals = frames.map((f) => {
    const out: V3[] = [];
    for (let k = 0; k < seg; k++) {
      const a = phi0 + (k * 2 * Math.PI) / seg;
      const ca = Math.cos(a), sa = Math.sin(a);
      out.push([
        sgn * (ca * f.e1[0] + sa * f.e2[0]),
        sgn * (ca * f.e1[1] + sa * f.e2[1]),
        sgn * (ca * f.e1[2] + sa * f.e2[2]),
      ]);
    }
    return out;
  });
  for (let i = 0; i + 1 < frames.length; i++) {
    for (let k = 0; k < seg; k++) {
      const k1 = (k + 1) % seg;
      b.quad(rings[i][k], rings[i][k1], rings[i + 1][k1], rings[i + 1][k], normals[i][k], normals[i][k1], normals[i + 1][k1], normals[i + 1][k]);
    }
  }
}

/** Straight cylinder between two points (closed with caps unless `open`). */
export function cylinderBetween(
  b: MeshBuilder, p0: V3, p1: V3, r: number, seg: number,
  opts: { inward?: boolean; caps?: boolean } = {},
): void {
  const ax = sub(p1, p0);
  const t = norm(ax);
  const e1 = norm(perp(t));
  const e2 = cross(t, e1);
  const frames: TubeFrame[] = [
    { c: p0, e1, e2, r },
    { c: p1, e1, e2, r },
  ];
  sweepTube(b, frames, seg, !!opts.inward);
  if (opts.caps ?? !opts.inward) {
    for (const [c, s] of [[p0, -1], [p1, 1]] as [V3, number][]) {
      const ring = ringPoints({ c, e1, e2, r }, seg);
      const n: V3 = [t[0] * s, t[1] * s, t[2] * s];
      for (let k = 0; k < seg; k++) {
        const a = ring[k], d = ring[(k + 1) % seg];
        b.tri(c[0], c[1], c[2], a[0], a[1], a[2], d[0], d[1], d[2], n, n, n);
      }
    }
  }
}

/** Axis-aligned box as a closed solid (outward normals), optional holes per face are not supported here. */
export function box(b: MeshBuilder, min: V3, max: V3): void {
  const [x0, y0, z0] = min, [x1, y1, z1] = max;
  planarFace(b, rect(x0, z0, x1, z1), [], [0, y0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(b, rect(x0, z0, x1, z1), [], [0, y1, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(b, rect(y0, z0, y1, z1), [], [x0, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(b, rect(y0, z0, y1, z1), [], [x1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(b, rect(x0, y0, x1, y1), [], [0, 0, z0], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
  planarFace(b, rect(x0, y0, x1, y1), [], [0, 0, z1], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
}

// ---------------------------------------------------------------------------
// Extrusions
// ---------------------------------------------------------------------------

export interface ExtrudeOpts {
  depth: number;
  bevel?: number;
  bevelSegments?: number;
  curveSegments?: number;
  creaseAngle?: number;
}

/**
 * Extrude a shape along +z from z = 0 to z = depth (bevels stay inside that
 * range), with creased normals so curved outlines shade smoothly.
 */
export function extrude(shape: THREE.Shape | THREE.Shape[], o: ExtrudeOpts): THREE.BufferGeometry {
  const bev = o.bevel ?? 0;
  const g = new THREE.ExtrudeGeometry(shape, {
    depth: Math.max(o.depth - 2 * bev, 1e-5),
    bevelEnabled: bev > 0,
    bevelThickness: bev,
    bevelSize: bev,
    bevelOffset: -bev,
    bevelSegments: o.bevelSegments ?? 2,
    curveSegments: o.curveSegments ?? 24,
    steps: 1,
  });
  g.translate(0, 0, bev);
  return crease(g, o.creaseAngle ?? (35 * Math.PI) / 180);
}

/** Circle as a THREE.Path (for holes), CW so it works as a hole. */
export function circlePath(cx: number, cy: number, r: number, seg = 48): THREE.Path {
  const p = new THREE.Path();
  for (let k = 0; k <= seg; k++) {
    const a = -(k * 2 * Math.PI) / seg;
    const x = cx + r * Math.cos(a), y = cy + r * Math.sin(a);
    if (k === 0) p.moveTo(x, y);
    else p.lineTo(x, y);
  }
  return p;
}

export function circleShape(cx: number, cy: number, r: number, seg = 64): THREE.Shape {
  const s = new THREE.Shape();
  for (let k = 0; k <= seg; k++) {
    const a = (k * 2 * Math.PI) / seg;
    const x = cx + r * Math.cos(a), y = cy + r * Math.sin(a);
    if (k === 0) s.moveTo(x, y);
    else s.lineTo(x, y);
  }
  return s;
}

export function roundedRectShape(x0: number, y0: number, x1: number, y1: number, rad: number): THREE.Shape {
  const s = new THREE.Shape();
  const r = Math.min(rad, (x1 - x0) / 2, (y1 - y0) / 2);
  s.moveTo(x0 + r, y0);
  s.lineTo(x1 - r, y0);
  s.absarc(x1 - r, y0 + r, r, -Math.PI / 2, 0, false);
  s.lineTo(x1, y1 - r);
  s.absarc(x1 - r, y1 - r, r, 0, Math.PI / 2, false);
  s.lineTo(x0 + r, y1);
  s.absarc(x0 + r, y1 - r, r, Math.PI / 2, Math.PI, false);
  s.lineTo(x0, y0 + r);
  s.absarc(x0 + r, y0 + r, r, Math.PI, 1.5 * Math.PI, false);
  return s;
}

/**
 * Spur-gear outline (trapezoidal teeth with rounded flanks) centred at the
 * origin. Tooth k is centred at angle phase + 2πk/teeth.
 */
export function gearShape(teeth: number, pitchR: number, opts: { module?: number; holeR?: number; phase?: number } = {}): THREE.Shape {
  const m = opts.module ?? (2 * pitchR) / teeth;
  const ra = pitchR + m; // addendum
  const rd = pitchR - 1.25 * m; // dedendum
  const phase = opts.phase ?? 0;
  const pitchAng = (2 * Math.PI) / teeth;
  const s = new THREE.Shape();
  let first = true;
  const pt = (r: number, a: number) => {
    const x = r * Math.cos(a), y = r * Math.sin(a);
    if (first) { s.moveTo(x, y); first = false; } else s.lineTo(x, y);
  };
  for (let k = 0; k < teeth; k++) {
    const c = phase + k * pitchAng;
    // root → flank → tip → flank → root, as fractions of the pitch angle
    pt(rd, c - 0.5 * pitchAng);
    pt(rd, c - 0.3 * pitchAng);
    pt(pitchR - 0.2 * m, c - 0.25 * pitchAng);
    pt(pitchR, c - 0.2 * pitchAng);
    pt(ra, c - 0.1 * pitchAng);
    pt(ra, c + 0.1 * pitchAng);
    pt(pitchR, c + 0.2 * pitchAng);
    pt(pitchR - 0.2 * m, c + 0.25 * pitchAng);
    pt(rd, c + 0.3 * pitchAng);
  }
  s.closePath();
  if (opts.holeR && opts.holeR > 0) s.holes.push(circlePath(0, 0, opts.holeR));
  return s;
}

// ---------------------------------------------------------------------------
// Clipping by the section plane
// ---------------------------------------------------------------------------

/**
 * Clip a triangle mesh to the half-space n·p ≤ d (Sutherland–Hodgman per
 * triangle, attributes interpolated). Winding is preserved. Returns a new
 * non-indexed geometry with position/normal/uv.
 */
export function clipGeometryPlane(src: THREE.BufferGeometry, n: V3, d: number): THREE.BufferGeometry {
  const g = src.index ? src.toNonIndexed() : src;
  const P = g.getAttribute('position') as THREE.BufferAttribute;
  const N = g.getAttribute('normal') as THREE.BufferAttribute | undefined;
  const T = g.getAttribute('uv') as THREE.BufferAttribute | undefined;
  const pos: number[] = [], nor: number[] = [], uv: number[] = [];
  // vertex record: x y z nx ny nz u v
  const vin: number[][] = [[], [], []];
  const poly: number[][] = [];
  const dist = new Float64Array(3);
  const lerp = (a: number[], b: number[], t: number) => a.map((v, i) => v + (b[i] - v) * t);
  for (let i = 0; i < P.count; i += 3) {
    for (let k = 0; k < 3; k++) {
      const j = i + k;
      vin[k] = [
        P.getX(j), P.getY(j), P.getZ(j),
        N ? N.getX(j) : 0, N ? N.getY(j) : 0, N ? N.getZ(j) : 1,
        T ? T.getX(j) : 0, T ? T.getY(j) : 0,
      ];
      dist[k] = n[0] * vin[k][0] + n[1] * vin[k][1] + n[2] * vin[k][2] - d;
    }
    if (dist[0] <= 0 && dist[1] <= 0 && dist[2] <= 0) {
      for (const v of vin) { pos.push(v[0], v[1], v[2]); nor.push(v[3], v[4], v[5]); uv.push(v[6], v[7]); }
      continue;
    }
    if (dist[0] > 0 && dist[1] > 0 && dist[2] > 0) continue;
    poly.length = 0;
    for (let k = 0; k < 3; k++) {
      const a = vin[k], b = vin[(k + 1) % 3];
      const da = dist[k], db = dist[(k + 1) % 3];
      if (da <= 0) poly.push(a);
      if ((da <= 0) !== (db <= 0)) poly.push(lerp(a, b, da / (da - db)));
    }
    for (let k = 1; k + 1 < poly.length; k++) {
      for (const v of [poly[0], poly[k], poly[k + 1]]) {
        const l = Math.hypot(v[3], v[4], v[5]) || 1;
        pos.push(v[0], v[1], v[2]);
        nor.push(v[3] / l, v[4] / l, v[5] / l);
        uv.push(v[6], v[7]);
      }
    }
  }
  if (g !== src) g.dispose();
  const out = new THREE.BufferGeometry();
  out.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  out.setAttribute('normal', new THREE.Float32BufferAttribute(nor, 3));
  out.setAttribute('uv', new THREE.Float32BufferAttribute(uv, 2));
  out.computeBoundingSphere();
  return out;
}

/** Keep z ≤ zMax. */
export function clipGeometryZ(src: THREE.BufferGeometry, zMax = 0): THREE.BufferGeometry {
  return clipGeometryPlane(src, [0, 0, 1], zMax);
}

/**
 * Cutaway region of a static assembly (object-local coordinates):
 *  - 'half': remove z > 0;
 *  - 'quadrant': remove {side·x > 0 and z > 0} (a quarter section: two cut
 *    faces, z = 0 for side·x ≥ 0 and x = 0 for z ≥ 0).
 */
export type CutMode = { kind: 'half' } | { kind: 'quadrant'; side: 1 | -1 };

export function clipForCut(src: THREE.BufferGeometry, mode: CutMode): THREE.BufferGeometry {
  if (mode.kind === 'half') return clipGeometryZ(src, 0);
  const a = mode.side;
  const keepSide = clipGeometryPlane(src, [a, 0, 0], 0);
  const other = clipGeometryPlane(src, [-a, 0, 0], 0);
  const behind = clipGeometryZ(other, 0);
  other.dispose();
  const m = mergeGeometries([keepSide, behind], false);
  keepSide.dispose();
  behind.dispose();
  if (!m) throw new Error('mergeGeometries failed');
  m.computeBoundingSphere();
  return m;
}

/** Keep only position/normal/uv (non-indexed) so geometries can be merged. */
export function normalizeAttributes(g: THREE.BufferGeometry): THREE.BufferGeometry {
  const src = g.index ? g.toNonIndexed() : g;
  const out = new THREE.BufferGeometry();
  out.setAttribute('position', src.getAttribute('position'));
  if (!src.getAttribute('normal')) src.computeVertexNormals();
  out.setAttribute('normal', src.getAttribute('normal'));
  const uv = src.getAttribute('uv');
  out.setAttribute('uv', uv ?? new THREE.Float32BufferAttribute(new Float32Array(src.getAttribute('position').count * 2), 2));
  return out;
}

export function merge(geoms: THREE.BufferGeometry[]): THREE.BufferGeometry {
  const norm = geoms.map(normalizeAttributes);
  const m = mergeGeometries(norm, false);
  if (!m) throw new Error('mergeGeometries failed');
  m.computeBoundingSphere();
  return m;
}

/**
 * Smooth normals across shared vertices within a crease angle. three's
 * toCreasedNormals quantises positions to 1e-2 units, so work in mm.
 */
export function crease(g: THREE.BufferGeometry, angle = (35 * Math.PI) / 180): THREE.BufferGeometry {
  const src = g.index ? g.toNonIndexed() : g;
  if (src !== g) g.dispose();
  src.scale(1000, 1000, 1000);
  const out = toCreasedNormals(src, angle);
  out.scale(0.001, 0.001, 0.001);
  if (out !== src) src.dispose();
  return out;
}

export { mergeVertices };

/** Signed volume (divergence theorem). Positive for closed meshes with outward normals. */
export function signedVolume(g: THREE.BufferGeometry): number {
  const P = g.getAttribute('position');
  const idx = g.index;
  const n = idx ? idx.count : P.count;
  let v = 0;
  const get = (i: number) => (idx ? idx.getX(i) : i);
  for (let i = 0; i < n; i += 3) {
    const a = get(i), b = get(i + 1), c = get(i + 2);
    const ax = P.getX(a), ay = P.getY(a), az = P.getZ(a);
    const bx = P.getX(b), by = P.getY(b), bz = P.getZ(b);
    const cx = P.getX(c), cy = P.getY(c), cz = P.getZ(c);
    v += ax * (by * cz - bz * cy) - ay * (bx * cz - bz * cx) + az * (bx * cy - by * cx);
  }
  return v / 6;
}

export function totalArea(g: THREE.BufferGeometry): number {
  const P = g.getAttribute('position');
  let a = 0;
  for (let i = 0; i < P.count; i += 3) {
    const ux = P.getX(i + 1) - P.getX(i), uy = P.getY(i + 1) - P.getY(i), uz = P.getZ(i + 1) - P.getZ(i);
    const vx = P.getX(i + 2) - P.getX(i), vy = P.getY(i + 2) - P.getY(i), vz = P.getZ(i + 2) - P.getZ(i);
    a += 0.5 * Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
  }
  return a;
}

// ---------------------------------------------------------------------------
// small vector helpers (build-time only)
// ---------------------------------------------------------------------------

export function sub(a: V3, b: V3): V3 { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
export function add(a: V3, b: V3): V3 { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
export function scale(a: V3, s: number): V3 { return [a[0] * s, a[1] * s, a[2] * s]; }
export function dot(a: V3, b: V3): number { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
export function cross(a: V3, b: V3): V3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
export function norm(a: V3): V3 {
  const l = Math.hypot(a[0], a[1], a[2]) || 1;
  return [a[0] / l, a[1] / l, a[2] / l];
}
/** Some unit vector perpendicular to a. */
export function perp(a: V3): V3 {
  const ax = Math.abs(a[0]), ay = Math.abs(a[1]), az = Math.abs(a[2]);
  if (ax <= ay && ax <= az) return norm(cross(a, [1, 0, 0]));
  if (ay <= az) return norm(cross(a, [0, 1, 0]));
  return norm(cross(a, [0, 0, 1]));
}
