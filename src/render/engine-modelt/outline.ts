/**
 * Plan-view (x, z) polygons for the L-head chamber: the bore circle, the rounded-rectangle valve pocket
 * (engine-spec.ts LHeadChamberSpec) and their union / difference, discretised so that every derived
 * outline shares the SAME vertices (the 3D cavity built from them is watertight: walls, roof and the
 * gasket-face hole meet edge to edge). Pure math, no three.js.
 *
 * Orientation: polygons are counter-clockwise in the (x, z) coordinate plane (positive shoelace area
 * Σ x_i z_{i+1} − x_{i+1} z_i); the interior lies to the LEFT of each edge.
 */
export type P2 = [number, number];

/** Circle polygon with the same vertex placement as geometry.ts ringXZ / lathe: (cx + r cos a, cz + r sin a). */
export function circlePolygon(cx: number, cz: number, r: number, seg: number, phi0 = 0): P2[] {
  const out: P2[] = [];
  for (let k = 0; k < seg; k++) {
    const a = phi0 + (k * 2 * Math.PI) / seg;
    out.push([cx + r * Math.cos(a), cz + r * Math.sin(a)]);
  }
  return out;
}

/** Rounded rectangle (CCW), `cornerSeg` segments per quarter arc. */
export function roundedRectPolygon(xMin: number, xMax: number, zMin: number, zMax: number, rc: number, cornerSeg = 8): P2[] {
  const r = Math.max(0, Math.min(rc, (xMax - xMin) / 2, (zMax - zMin) / 2));
  const out: P2[] = [];
  const corner = (cx: number, cz: number, a0: number) => {
    for (let k = 0; k <= cornerSeg; k++) {
      const a = a0 + (k * Math.PI) / 2 / cornerSeg;
      out.push([cx + r * Math.cos(a), cz + r * Math.sin(a)]);
    }
  };
  if (r <= 0) return [[xMin, zMin], [xMax, zMin], [xMax, zMax], [xMin, zMax]];
  corner(xMax - r, zMin + r, -Math.PI / 2);
  corner(xMax - r, zMax - r, 0);
  corner(xMin + r, zMax - r, Math.PI / 2);
  corner(xMin + r, zMin + r, Math.PI);
  return out;
}

/** Shoelace area (positive for CCW). */
export function polygonArea(p: readonly P2[]): number {
  let s = 0;
  for (let i = 0; i < p.length; i++) {
    const a = p[i], b = p[(i + 1) % p.length];
    s += a[0] * b[1] - b[0] * a[1];
  }
  return s / 2;
}

/** Point strictly inside a convex CCW polygon. */
export function insideConvex(p: readonly P2[], x: number, z: number): boolean {
  for (let i = 0; i < p.length; i++) {
    const a = p[i], b = p[(i + 1) % p.length];
    if ((b[0] - a[0]) * (z - a[1]) - (b[1] - a[1]) * (x - a[0]) <= 0) return false;
  }
  return true;
}

/** Point inside a simple polygon (even–odd rule). */
export function insidePolygon(p: readonly P2[], x: number, z: number): boolean {
  let inside = false;
  for (let i = 0, j = p.length - 1; i < p.length; j = i++) {
    const a = p[i], b = p[j];
    if (a[1] > z !== b[1] > z && x < ((b[0] - a[0]) * (z - a[1])) / (b[1] - a[1]) + a[0]) inside = !inside;
  }
  return inside;
}

/** Minimum distance from a point to a polygon's boundary. */
export function distanceToPolygon(p: readonly P2[], x: number, z: number): number {
  let best = Infinity;
  for (let i = 0; i < p.length; i++) {
    const a = p[i], b = p[(i + 1) % p.length];
    const dx = b[0] - a[0], dz = b[1] - a[1];
    const t = Math.max(0, Math.min(1, ((x - a[0]) * dx + (z - a[1]) * dz) / (dx * dx + dz * dz || 1)));
    best = Math.min(best, Math.hypot(a[0] + t * dx - x, a[1] + t * dz - z));
  }
  return best;
}

export interface ChamberOutlines {
  /** Gasket-face opening: boundary of (bore disc ∪ pocket), CCW. */
  union: P2[];
  /** Pocket roof region: boundary of (pocket \ bore disc), CCW. */
  pocketRoof: P2[];
  /** The bore circle with the two crossing points inserted, CCW (walls above the pocket roof, bore roof). */
  circle: P2[];
}

/**
 * Union and difference outlines of a circle (centre cx, cz, radius r, `seg` vertices) and a convex
 * CCW polygon that crosses it exactly twice (the valve pocket overlapping one side of the bore).
 * The crossing points are computed on the TRUE circle against the polygon edges and inserted into both
 * boundaries, so all three outlines share their vertices.
 */
export function chamberOutlines(cx: number, cz: number, r: number, seg: number, pocket: readonly P2[]): ChamberOutlines {
  const n = pocket.length;
  const inC = (q: P2) => (q[0] - cx) ** 2 + (q[1] - cz) ** 2 < r * r;
  // pocket vertices inside the circle form one cyclic run
  const flags = pocket.map(inC);
  let start = -1;
  for (let i = 0; i < n; i++) if (flags[i] && !flags[(i - 1 + n) % n]) { start = i; break; }
  if (start < 0) throw new Error('chamberOutlines: pocket must overlap the bore circle on one side');
  let end = start;
  while (flags[(end + 1) % n]) end = (end + 1) % n;
  const cross = (a: P2, b: P2): P2 => {
    // segment a→b crosses the circle once (a and b on different sides)
    const dx = b[0] - a[0], dz = b[1] - a[1];
    const fx = a[0] - cx, fz = a[1] - cz;
    const A = dx * dx + dz * dz, B = 2 * (fx * dx + fz * dz), C = fx * fx + fz * fz - r * r;
    const disc = Math.sqrt(Math.max(B * B - 4 * A * C, 0));
    const t1 = (-B - disc) / (2 * A), t2 = (-B + disc) / (2 * A);
    const t = t1 >= -1e-12 && t1 <= 1 + 1e-12 ? t1 : t2;
    return [a[0] + t * dx, a[1] + t * dz];
  };
  const qIn = cross(pocket[(start - 1 + n) % n], pocket[start]); // pocket boundary enters the disc
  const qOut = cross(pocket[end], pocket[(end + 1) % n]); // and leaves it
  const ang = (q: P2) => Math.atan2(q[1] - cz, q[0] - cx);
  const twoPi = 2 * Math.PI;
  const norm = (a: number) => ((a % twoPi) + twoPi) % twoPi;
  const aIn = norm(ang(qIn)), aOut = norm(ang(qOut));
  // circle vertices: angles a_k = 2πk/seg; outside-pocket arc runs CCW from qIn to qOut
  const ccwArc = (from: number, to: number): P2[] => {
    const pts: { da: number; p: P2 }[] = [];
    const span = norm(to - from);
    for (let k = 0; k < seg; k++) {
      const a = (k * twoPi) / seg;
      const da = norm(a - from);
      if (da > 1e-9 && da < span - 1e-9) pts.push({ da, p: [cx + r * Math.cos(a), cz + r * Math.sin(a)] });
    }
    pts.sort((p, q) => p.da - q.da);
    return pts.map((e) => e.p);
  };
  const arcOut = ccwArc(aIn, aOut); // outside the pocket
  const arcIn = ccwArc(aOut, aIn); // inside the pocket (CCW from qOut to qIn)
  const pocketOutside: P2[] = [];
  for (let i = (end + 1) % n; i !== start; i = (i + 1) % n) pocketOutside.push([pocket[i][0], pocket[i][1]]);
  const union: P2[] = [qOut, ...pocketOutside, qIn, ...arcOut];
  const pocketRoof: P2[] = [qOut, ...pocketOutside, qIn, ...arcIn.slice().reverse()];
  const circle: P2[] = [qIn, ...arcOut, qOut, ...arcIn];
  return { union, pocketRoof, circle };
}
