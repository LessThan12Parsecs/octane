/**
 * Static ROOT-frame housings: crankcase (walls, end plates with main-bearing
 * bosses, deck with the cylinder sleeve bore, intermediate main bearing,
 * tappet guide block, cam bearing bracket), bedplate, drive motor body and
 * V-belts, worm bearing blocks.
 */
import * as THREE from 'three';
import type { EngineLayout } from './layout';
import type { StaticSet } from './static-set';
import { MeshBuilder, type V3, box, circlePath, extrude, lathe } from './geometry';

/** Plan (x, z) shape for extrusion along +y: shape v = −z. */
function planShape(pts: [number, number][], holes: { x: number; z: number; r: number }[]): THREE.Shape {
  const s = new THREE.Shape(pts.map(([x, z]) => new THREE.Vector2(x, -z)));
  for (const h of holes) s.holes.push(circlePath(h.x, -h.z, h.r, 48));
  return s;
}

function extrudeY(shape: THREE.Shape, y0: number, y1: number, bevel = 0): THREE.BufferGeometry {
  const g = extrude(shape, { depth: y1 - y0, bevel, bevelSegments: 1, curveSegments: 32 });
  g.rotateX(-Math.PI / 2); // (u, v, d) → (u, d, −v)
  g.translate(0, y0, 0);
  return g;
}

function extrudeZ(shape: THREE.Shape, z0: number, z1: number, bevel = 0): THREE.BufferGeometry {
  const g = extrude(shape, { depth: z1 - z0, bevel, bevelSegments: 1, curveSegments: 32 });
  g.translate(0, 0, z0);
  return g;
}

/** Ring along z (bearing shell, boss): closed solid, outward normals. */
export function ringZ(z0: number, z1: number, rIn: number, rOut: number, cx = 0, cy = 0, seg = 64): THREE.BufferGeometry {
  const b = new MeshBuilder();
  lathe(b, [{ r: rIn, y: z0 }, { r: rOut, y: z0 }, { r: rOut, y: z1 }, { r: rIn, y: z1 }], { segments: seg, closed: true });
  const g = b.build();
  g.rotateX(Math.PI / 2); // lathe axial y → z
  g.translate(cx, cy, 0);
  return g;
}

export function buildCrankcase(L: EngineLayout, set: StaticSet): void {
  const c = L.crankcase;
  const xo = c.xHalfOuter, xi = c.xHalfInner;
  const fl = c.baseFlangeHalf;
  const flH = 0.018;
  const J = L.mainJournalRadius;
  const bore = J + 0.0025;

  // ---- U-profile side walls + floor (between the end plates) ----
  const u = new THREE.Shape();
  const rf = 0.03;
  u.moveTo(-xo, c.yDeckBottom);
  u.lineTo(-xo, c.yBase + flH + 0.02);
  u.quadraticCurveTo(-xo, c.yBase + flH, -xo - 0.02, c.yBase + flH);
  u.lineTo(-fl, c.yBase + flH);
  u.lineTo(-fl, c.yBase);
  u.lineTo(fl, c.yBase);
  u.lineTo(fl, c.yBase + flH);
  u.lineTo(xo + 0.02, c.yBase + flH);
  u.quadraticCurveTo(xo, c.yBase + flH, xo, c.yBase + flH + 0.02);
  u.lineTo(xo, c.yDeckBottom);
  u.lineTo(xi, c.yDeckBottom);
  u.lineTo(xi, c.yFloorInner + rf);
  u.quadraticCurveTo(xi, c.yFloorInner, xi - rf, c.yFloorInner);
  u.lineTo(-xi + rf, c.yFloorInner);
  u.quadraticCurveTo(-xi, c.yFloorInner, -xi, c.yFloorInner + rf);
  u.lineTo(-xi, c.yDeckBottom);
  u.closePath();
  set.add('paint', extrudeZ(u, c.zBackInner, c.zFrontInner));

  // ---- end plates (outer outline with bearing bores) ----
  const plate = (holes: { x: number; y: number; r: number }[]) => {
    const s = new THREE.Shape();
    s.moveTo(-xo, c.yDeckBottom);
    s.lineTo(-xo, c.yBase + flH + 0.02);
    s.quadraticCurveTo(-xo, c.yBase + flH, -xo - 0.02, c.yBase + flH);
    s.lineTo(-fl, c.yBase + flH);
    s.lineTo(-fl, c.yBase);
    s.lineTo(fl, c.yBase);
    s.lineTo(fl, c.yBase + flH);
    s.lineTo(xo + 0.02, c.yBase + flH);
    s.quadraticCurveTo(xo, c.yBase + flH, xo, c.yBase + flH + 0.02);
    s.lineTo(xo, c.yDeckBottom);
    s.closePath();
    for (const h of holes) s.holes.push(circlePath(h.x, h.y, h.r, 64));
    return s;
  };
  set.add('paint', extrudeZ(plate([{ x: 0, y: 0, r: bore }]), c.zFrontInner, c.zFrontOuter, 0.002));
  const camBore = L.camShaftRadius + 0.0008;
  set.add('paint', extrudeZ(plate([{ x: 0, y: 0, r: bore }, { x: L.camX, y: L.camY, r: camBore }]), c.zBackOuter, c.zBackInner, 0.002));

  // main-bearing bosses outside the plates + bearing shells through them
  const bossR = J + 0.024;
  set.add('paint', ringZ(c.zFrontOuter - 0.001, c.zFrontOuter + 0.012, bore, bossR));
  set.add('paint', ringZ(c.zBackOuter - 0.012, c.zBackOuter + 0.001, bore, bossR));
  set.add('bearing', ringZ(c.zFrontInner - 0.004, c.zFrontOuter + 0.012, J + 0.0002, bore), { section: 'light' });
  set.add('bearing', ringZ(c.zBackOuter - 0.012, c.zBackInner + 0.004, J + 0.0002, bore), { section: 'light' });
  // cam bearing boss on the back plate
  set.add('paint', ringZ(c.zBackOuter - 0.008, c.zBackOuter + 0.001, camBore, camBore + 0.012, L.camX, L.camY, 48));

  // ---- deck with the cylinder-sleeve bore and pushrod holes ----
  const pushHoles = L.valves.map((v) => {
    // pushrod x at deck height (rest), straight line tappet → rocker cup
    const yb = L.camY + L.camBaseRadius + L.tappetLength;
    const yt = L.headY(L.spec.geometry.compressionRatio) + v.pushrodTopRest[1];
    const t = ((c.yDeckBottom + c.yDeckTop) / 2 - yb) / (yt - yb);
    return { x: v.tappetX + (v.pushrodTopRest[0] - v.tappetX) * t, z: v.lobeZ, r: L.pushrodPassageRadius };
  });
  const deck = planShape(
    [[-xo, c.zBackOuter], [xo, c.zBackOuter], [xo, c.zFrontOuter], [-xo, c.zFrontOuter]],
    [{ x: 0, z: 0, r: L.spigotRadius + 0.003 }, ...pushHoles],
  );
  set.add('paint', extrudeY(deck, c.yDeckBottom, c.yDeckTop, 0.0015));

  // ---- intermediate main bearing pedestal ----
  {
    const w = J + 0.032;
    const s = new THREE.Shape();
    s.moveTo(-w, c.yFloorInner);
    s.lineTo(w, c.yFloorInner);
    s.lineTo(w, 0);
    s.absarc(0, 0, w, 0, Math.PI, false);
    s.lineTo(-w, c.yFloorInner);
    s.holes.push(circlePath(0, 0, bore, 64));
    const z0 = c.pedestalZ - c.pedestalThickness / 2, z1 = c.pedestalZ + c.pedestalThickness / 2;
    set.add('paint', extrudeZ(s, z0, z1, 0.0015));
    set.add('bearing', ringZ(z0 - 0.002, z1 + 0.002, J + 0.0002, bore), { section: 'light' });
  }

  // ---- tappet guide block (hangs from the deck) ----
  {
    const tg = L.tappetGuide;
    const zs = L.valves.map((v) => v.lobeZ);
    const x0 = L.camX - L.tappetRadius - 0.012, x1 = L.camX + L.tappetRadius + 0.012;
    const z0 = Math.min(...zs) - L.tappetRadius - 0.012, z1 = Math.max(...zs) + L.tappetRadius + 0.012;
    const s = planShape([[x0, z0], [x1, z0], [x1, z1], [x0, z1]], L.valves.map((v) => ({ x: v.tappetX, z: v.lobeZ, r: L.tappetRadius + 0.0004 })));
    set.add('paint', extrudeY(s, tg.yBottom, c.yDeckBottom + 0.0005, 0.0015));
  }

  // ---- cam bearing bracket (front end of the camshaft) ----
  {
    const w = L.camShaftRadius + 0.014;
    const s = new THREE.Shape();
    s.moveTo(L.camX - w, c.yDeckBottom + 0.0005);
    s.lineTo(L.camX - w, L.camY);
    s.absarc(L.camX, L.camY, w, Math.PI, 2 * Math.PI, false);
    s.lineTo(L.camX + w, c.yDeckBottom + 0.0005);
    s.closePath();
    s.holes.push(circlePath(L.camX, L.camY, L.camShaftRadius + 0.0008, 48));
    set.add('paint', extrudeZ(s, L.camZFront - 0.009, L.camZFront, 0.001));
  }

  // ---- bedplate ----
  {
    const b = new MeshBuilder();
    const xm = L.motor.x;
    const x0 = Math.min(xm - L.motor.radius - 0.08, -c.baseFlangeHalf - 0.12);
    const x1 = Math.max(xm + L.motor.radius + 0.08, c.baseFlangeHalf + 0.12);
    const z0 = L.flywheelZ - L.flywheelWidth / 2 - 0.08, z1 = c.zFrontOuter + 0.1;
    box(b, [x0, c.yBase - 0.045, z0], [x1, c.yBase, z1]);
    set.add('base', b.build(), { cut: false, section: false });
  }

  // ---- motor body (static), feet ----
  {
    const m = L.motor;
    const b = new MeshBuilder();
    const z0 = m.zCenter - m.length / 2, z1 = m.zCenter + m.length / 2;
    const R = m.radius;
    const prof: { r: number; y: number; smooth?: boolean }[] = [
      { r: 0.02, y: z0 - 0.02 }, { r: R * 0.55, y: z0 - 0.02 }, { r: R * 0.85, y: z0 - 0.004, smooth: true }, { r: R * 0.92, y: z0 },
    ];
    const nf = 9;
    for (let i = 0; i <= nf; i++) {
      const y = z0 + 0.01 + ((m.length - 0.02) * i) / nf;
      prof.push({ r: R * 0.93, y: y - 0.004 });
      prof.push({ r: R, y: y - 0.003 });
      prof.push({ r: R, y: y + 0.003 });
      prof.push({ r: R * 0.93, y: y + 0.004 });
    }
    prof.push({ r: R * 0.92, y: z1 }, { r: R * 0.85, y: z1 + 0.004, smooth: true }, { r: R * 0.55, y: z1 + 0.02 }, { r: 0.02, y: z1 + 0.02 });
    lathe(b, prof, { segments: 64 });
    const g = b.build();
    g.rotateX(Math.PI / 2);
    g.translate(m.x, m.y, 0);
    set.add('motor', g, { cut: false, section: false });
    const feet = new MeshBuilder();
    box(feet, [m.x - R * 0.8, c.yBase, m.zCenter - m.length * 0.4], [m.x + R * 0.8, c.yBase + 0.012, m.zCenter + m.length * 0.4]);
    box(feet, [m.x - R * 0.55, c.yBase, m.zCenter - m.length * 0.35], [m.x + R * 0.55, m.y - R * 0.6, m.zCenter + m.length * 0.35]);
    box(feet, [m.x - 0.035, m.y + R - 0.01, m.zCenter - 0.04], [m.x + 0.035, m.y + R + 0.045, m.zCenter + 0.04]); // terminal box
    set.add('motor', feet.build(), { cut: false, section: false });
    const shaft = new MeshBuilder();
    lathe(shaft, [{ r: 0, y: L.flywheelZ - 0.03 }, { r: 0.014, y: L.flywheelZ - 0.03 }, { r: 0.014, y: z0 - 0.02 }, { r: 0, y: z0 - 0.02 }], { segments: 32 });
    const sg = shaft.build();
    sg.rotateX(Math.PI / 2);
    sg.translate(m.x, m.y, 0);
    set.add('steelPolished', sg, { cut: false, section: false });
  }

  // ---- V-belts (open drive, flywheel ↔ motor pulley) ----
  {
    const b = new MeshBuilder();
    const m = L.motor;
    const r1 = L.flywheelGrooveRadius - 0.001, r2 = m.pulleyRadius - 0.001;
    const C1: [number, number] = [0, 0], C2: [number, number] = [m.x, m.y];
    const dx = C2[0] - C1[0], dy = C2[1] - C1[1];
    const d = Math.hypot(dx, dy);
    const phi = Math.atan2(dy, dx);
    const beta = Math.asin((r1 - r2) / d);
    // tangent points share the normal n with n·û = (r1 − r2)/d = sin β
    const aP = phi + Math.PI / 2 - beta;
    const aM = phi - Math.PI / 2 + beta;
    const path: [number, number][] = [];
    const n1 = 96, n2 = 48;
    // around C1 from aP to aM the long way (away from C2)
    let a0 = aP, a1 = aM;
    while (a1 < a0) a1 += 2 * Math.PI;
    for (let i = 0; i <= n1; i++) {
      const a = a0 + ((a1 - a0) * i) / n1;
      path.push([C1[0] + r1 * Math.cos(a), C1[1] + r1 * Math.sin(a)]);
    }
    // around C2 from aM to aP through φ (the side facing away from C1)
    a0 = aM; a1 = aP;
    while (a1 < a0) a1 += 2 * Math.PI;
    for (let i = 0; i <= n2; i++) {
      const a = a0 + ((a1 - a0) * i) / n2;
      path.push([C2[0] + r2 * Math.cos(a), C2[1] + r2 * Math.sin(a)]);
    }
    const pitch = L.flywheelWidth / 3;
    for (let k = 1; k <= 2; k++) {
      const zc = L.flywheelZ - L.flywheelWidth / 2 + k * pitch;
      sweepBelt(b, path, zc, 0.012, 0.009);
    }
    set.add('rubber', b.build(), { cut: false, section: false });
  }

  // ---- worm bearing blocks on the deck ----
  {
    const b = new MeshBuilder();
    const w = L.worm;
    for (const z of [w.zBack + 0.006, c.zFrontOuter - 0.012]) {
      box(b, [w.x - 0.014, c.yDeckTop, z - 0.008], [w.x + 0.014, w.y - w.radius * 0.55 - 0.0005, z + 0.008]);
      const g = new MeshBuilder();
      lathe(g, [{ r: w.radius * 0.55 + 0.0005, y: -0.008 }, { r: 0.014, y: -0.008 }, { r: 0.014, y: 0.008 }, { r: w.radius * 0.55 + 0.0005, y: 0.008 }], { segments: 32, closed: true });
      const gg = g.build();
      gg.rotateX(Math.PI / 2);
      gg.translate(w.x, w.y, z);
      b.append(gg);
      gg.dispose();
    }
    set.add('paint', b.build(), { cut: false, section: false });
  }
}

/** Sweep a (trapezoidal V-belt) cross-section along a closed path in the x-y plane at z = zc. */
function sweepBelt(b: MeshBuilder, path: [number, number][], zc: number, width: number, thick: number): void {
  const n = path.length;
  const sec: V3[][] = [];
  const nor: V3[][] = [];
  for (let i = 0; i < n; i++) {
    const p = path[i], q = path[(i + 1) % n], o = path[(i - 1 + n) % n];
    let tx = q[0] - o[0], ty = q[1] - o[1];
    const tl = Math.hypot(tx, ty) || 1;
    tx /= tl; ty /= tl;
    // outward normal of the loop (belt runs CCW-ish): pick left normal and check against centroid later
    const nx = ty, ny = -tx;
    const top = width / 2, bot = width * 0.32;
    const corners: V3[] = [
      [p[0] + nx * thick / 2, p[1] + ny * thick / 2, zc - top],
      [p[0] + nx * thick / 2, p[1] + ny * thick / 2, zc + top],
      [p[0] - nx * thick / 2, p[1] - ny * thick / 2, zc + bot],
      [p[0] - nx * thick / 2, p[1] - ny * thick / 2, zc - bot],
    ];
    sec.push(corners);
    nor.push([[nx, ny, 0], [0, 0, 1], [-nx, -ny, 0], [0, 0, -1]]);
  }
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    for (let f = 0; f < 4; f++) {
      const f1 = (f + 1) % 4;
      const nf: V3 = f === 0 ? nor[i][0] : f === 1 ? [0, 0, 1] : f === 2 ? nor[i][2] : [0, 0, -1];
      const nfj: V3 = f === 0 ? nor[j][0] : f === 1 ? [0, 0, 1] : f === 2 ? nor[j][2] : [0, 0, -1];
      b.quad(sec[i][f], sec[i][f1], sec[j][f1], sec[j][f], nf, nf, nfj, nfj);
    }
  }
}
