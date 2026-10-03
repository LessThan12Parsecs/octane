/**
 * Rotating crank assembly (crank-local frame: crank axis = z, crank pin at
 * (0, r) — the group is rotated by −θ about z): main journals, two webs with
 * counterweights, crank pin, timing gear and the flywheel with belt grooves
 * and a degree scale (the CFR flywheel carries a timing scale).
 */
import * as THREE from 'three';
import type { EngineLayout } from './layout';
import type { EngineMaterials } from './materials';
import { MeshBuilder, extrude, gearShape, lathe, type ProfilePoint, crease, box } from './geometry';

function hull2(points: [number, number][]): [number, number][] {
  const p = points.slice().sort((a, b) => a[0] - b[0] || a[1] - b[1]);
  const cr = (o: [number, number], a: [number, number], b: [number, number]) =>
    (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
  const lo: [number, number][] = [];
  for (const q of p) {
    while (lo.length >= 2 && cr(lo[lo.length - 2], lo[lo.length - 1], q) <= 0) lo.pop();
    lo.push(q);
  }
  const up: [number, number][] = [];
  for (let i = p.length - 1; i >= 0; i--) {
    const q = p[i];
    while (up.length >= 2 && cr(up[up.length - 2], up[up.length - 1], q) <= 0) up.pop();
    up.push(q);
  }
  up.pop();
  lo.pop();
  return lo.concat(up);
}

/** Shaft segment along z (lathe about y, then rotated so y → z). */
export function shaftZ(z0: number, z1: number, radius: number, chamfer = 0.0008, seg = 48): THREE.BufferGeometry {
  const b = new MeshBuilder();
  const c = Math.min(chamfer, radius * 0.3, Math.abs(z1 - z0) * 0.3);
  const lo = Math.min(z0, z1), hi = Math.max(z0, z1);
  const prof: ProfilePoint[] = [
    { r: 0, y: lo },
    { r: radius - c, y: lo },
    { r: radius, y: lo + c },
    { r: radius, y: hi - c },
    { r: radius - c, y: hi },
    { r: 0, y: hi },
  ];
  lathe(b, prof, { segments: seg });
  const g = b.build();
  g.rotateX(Math.PI / 2); // y → z  (note: z' = y)
  return g;
}

export function buildCrank(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'crankshaft';
  const r = L.crankRadius;

  // --- webs (convex pear: pin boss + journal boss + counterweight sector) ---
  const pts: [number, number][] = [];
  const circ = (cx: number, cy: number, rad: number, a0 = 0, a1 = Math.PI * 2, n = 48) => {
    for (let i = 0; i <= n; i++) {
      const a = a0 + ((a1 - a0) * i) / n;
      pts.push([cx + rad * Math.cos(a), cy + rad * Math.sin(a)]);
    }
  };
  circ(0, r, L.webPinBossRadius);
  circ(0, 0, L.webJournalBossRadius);
  const h = L.counterweightHalfAngle;
  circ(0, 0, L.counterweightRadius, -Math.PI / 2 - h, -Math.PI / 2 + h, 64);
  // soften counterweight corners
  const cwr = L.counterweightRadius - 0.006;
  circ(0, 0, cwr, -Math.PI / 2 - h - 0.06, -Math.PI / 2 - h, 4);
  circ(0, 0, cwr, -Math.PI / 2 + h, -Math.PI / 2 + h + 0.06, 4);
  const hp = hull2(pts);
  const webShape = new THREE.Shape(hp.map(([x, y]) => new THREE.Vector2(x, y)));
  for (const s of [-1, 1]) {
    const w = extrude(webShape, { depth: L.webThickness, bevel: 0.0018, bevelSegments: 2 });
    w.translate(0, 0, s * L.webZ - L.webThickness / 2);
    const m = new THREE.Mesh(w, M.m.steelForged);
    m.name = 'crank-web';
    g.add(m);
  }

  // --- crank pin (runs into both webs) ---
  const pin = shaftZ(-L.webZ, L.webZ, L.crankPinRadius, 0.0005);
  pin.translate(0, r, 0);
  g.add(named(new THREE.Mesh(pin, M.m.steelPolished), 'crank-pin'));
  // fillet collars where the pin meets the webs
  for (const s of [-1, 1]) {
    const f = shaftZ(s * (L.crankPinLength / 2) - 0.0015, s * (L.crankPinLength / 2), L.crankPinRadius + 0.0025, 0.0012);
    f.translate(0, r, 0);
    g.add(new THREE.Mesh(f, M.m.steelPolished));
  }

  // --- main journals ---
  const zf0 = L.webZ + L.webThickness / 2 - 0.002;
  g.add(named(new THREE.Mesh(shaftZ(zf0, L.frontJournalEndZ, L.mainJournalRadius), M.m.steelPolished), 'main-journal-front'));
  g.add(named(new THREE.Mesh(shaftZ(-zf0, L.rearJournalEndZ, L.mainJournalRadius), M.m.steelPolished), 'main-journal-rear'));
  // front nose: keyed hub
  const nose = shaftZ(L.frontJournalEndZ - 0.002, L.frontJournalEndZ + 0.018, L.mainJournalRadius * 0.72);
  g.add(new THREE.Mesh(nose, M.m.steelForged));

  // --- timing gear (tooth at +γ toward the cam at θ = 0) ---
  const gamma = Math.atan2(L.camY, L.camX);
  const gearG = extrude(gearShape(L.crankGearTeeth, L.crankGearRadius, { phase: gamma, holeR: L.mainJournalRadius * 0.98 }), {
    depth: L.gearWidth, bevel: 0.0006, bevelSegments: 1, curveSegments: 4, creaseAngle: 0.5,
  });
  gearG.translate(0, 0, L.gearZ - L.gearWidth / 2);
  g.add(named(new THREE.Mesh(gearG, M.m.steelForged), 'crank-gear'));

  g.add(buildFlywheel(L, M));
  return g;
}

function named<T extends THREE.Object3D>(o: T, name: string): T {
  o.name = name;
  return o;
}

export function buildFlywheel(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'flywheel';
  const R = L.flywheelRadius;
  const W = L.flywheelWidth;
  const hubR = L.mainJournalRadius + 0.035;
  const hubW = W * 1.1;
  const webT = W * 0.32;
  const rimIn = R - 0.05;
  // Profile in (r, axial y) — CCW loop; axial y becomes z after rotation (y → z).
  const gr = 0.009; // groove depth
  const grooveHalfTop = 0.0075;
  const nG = 2;
  const outer: ProfilePoint[] = [];
  // rim outer surface with V grooves, from back (−W/2) to front (+W/2)
  outer.push({ r: R - 0.003, y: -W / 2 });
  outer.push({ r: R, y: -W / 2 + 0.003 });
  const pitch = W / (nG + 1);
  for (let i = 1; i <= nG; i++) {
    const c = -W / 2 + i * pitch;
    outer.push({ r: R, y: c - grooveHalfTop });
    outer.push({ r: R - gr, y: c - grooveHalfTop * 0.35 });
    outer.push({ r: R - gr, y: c + grooveHalfTop * 0.35 });
    outer.push({ r: R, y: c + grooveHalfTop });
  }
  outer.push({ r: R, y: W / 2 - 0.003 });
  outer.push({ r: R - 0.003, y: W / 2 });
  // Build the full CCW loop: start at hub bore bottom (back face) → …
  // Loop (r right, axial y up): back face goes outward (right), outer rim goes up (front), front face goes inward.
  const loop: ProfilePoint[] = [
    { r: 0, y: -hubW / 2 },
    { r: hubR, y: -hubW / 2 },
    { r: hubR, y: -webT / 2, smooth: false },
    { r: rimIn, y: -webT / 2 },
    { r: rimIn, y: -W / 2 },
    ...outer,
    { r: rimIn, y: W / 2 },
    { r: rimIn, y: webT / 2 },
    { r: hubR, y: webT / 2 },
    { r: hubR, y: hubW / 2 },
    { r: 0, y: hubW / 2 },
  ];
  // Separate materials: machined rim (outer), painted web+hub.
  const rimB = new MeshBuilder();
  const idxRim0 = 4, idxRim1 = 4 + outer.length + 1;
  lathe(rimB, loop.slice(idxRim0, idxRim1 + 1), { segments: 128 });
  const bodyB = new MeshBuilder();
  lathe(bodyB, loop.slice(0, idxRim0 + 1), { segments: 128 });
  lathe(bodyB, loop.slice(idxRim1), { segments: 128 });
  const rim = rimB.build();
  const body = bodyB.build();
  for (const geo of [rim, body]) {
    geo.rotateX(Math.PI / 2);
    geo.translate(0, 0, L.flywheelZ);
  }
  g.add(named(new THREE.Mesh(rim, M.m.ironMachined), 'flywheel-rim'));
  g.add(named(new THREE.Mesh(body, M.m.paint), 'flywheel-web'));

  // Degree scale on the front rim face: ticks every 5°, long ticks every 10°,
  // TDC (θ = 0, crank pin up) marked with a wide bright tick.
  const tickB = new MeshBuilder();
  const tdcB = new MeshBuilder();
  const zFace = L.flywheelZ + W / 2;
  for (let d = 0; d < 360; d += 5) {
    const long = d % 10 === 0;
    const len = d === 0 ? 0.03 : long ? 0.016 : 0.009;
    const wid = d === 0 ? 0.004 : 0.0013;
    // crank-local angle of the pin is +90°; θ increases clockwise → tick for θ=d at angle 90° − d
    const a = ((90 - d) * Math.PI) / 180;
    const bld = d === 0 ? tdcB : tickB;
    const tmp = new MeshBuilder();
    box(tmp, [-wid / 2, R - 0.004 - len, 0], [wid / 2, R - 0.004, 0.0006]);
    const geo = tmp.build();
    geo.rotateZ(a - Math.PI / 2);
    geo.translate(0, 0, zFace);
    bld.append(geo);
    geo.dispose();
  }
  g.add(named(new THREE.Mesh(tickB.build(), M.m.plated), 'flywheel-scale'));
  g.add(named(new THREE.Mesh(tdcB.build(), M.m.copper), 'flywheel-tdc'));

  // Hub bolts (6) on the front hub face.
  const boltB = new MeshBuilder();
  for (let k = 0; k < 6; k++) {
    const a = (k * Math.PI) / 3 + Math.PI / 6;
    const bx = (hubR - 0.012) * Math.cos(a), by = (hubR - 0.012) * Math.sin(a);
    lathe(boltB, [
      { r: 0, y: 0 }, { r: 0.0065, y: 0 }, { r: 0.0065, y: 0.005 }, { r: 0.005, y: 0.0065 }, { r: 0, y: 0.0065 },
    ], { segments: 6, center: [bx, 0, -by] });
  }
  const bolts = boltB.build();
  bolts.rotateX(Math.PI / 2);
  bolts.translate(0, 0, L.flywheelZ + hubW / 2);
  g.add(new THREE.Mesh(crease(bolts, 0.3), M.m.steelDark));
  return g;
}
