/**
 * En-bloc cylinder casting with the integral upper crankcase (one closed solid with every internal
 * cavity: four bores, the stepped crankcase cavity open to the pan, the valve (spring/tappet) chamber
 * open to the valve door, the eight seat countersinks with their cored ports — four exhaust ports
 * bending out to the manifold face and two siamesed intake galleries — and a water-jacket pocket in
 * front of cylinder 1), plus main-bearing caps and webs, cam bearings, the pressed-steel pan, the timing
 * gear cover, the valve door, the head gasket and the generator. World coordinates (ROOT).
 *
 * The block shell is a watertight 2-manifold (parts-block.test.ts): every hole in a planar face uses the
 * exact ring of the passage that meets it (geometry.ts lathe / sweepTube vertex placement).
 * Sizes not taken from layout.ts (fillets, chamfers, bosses, caps, fastener heads) are render-only
 * cosmetic choices — UNVERIFIED.
 */
import * as THREE from 'three';
import {
  MeshBuilder, type TubeFrame, type V3, box, cylinderBetween, lathe, planarFace, rect, ringPoints, ringXZ, sweepTube,
} from '../engine/geometry';
import { BORE_SEG, PORT_SEG, type ModelTLayout, type ModelTValveLayout } from './layout';
import type { FrameSet } from './frame-set';
import { boxCavity, extrudePlanY, extrudeZ, hullOfCircles, latheZ, polyShape, ringSolidZ } from './geom';
import type { P2 } from './outline';

export interface ShellBuilders {
  paint: MeshBuilder;
  machined: MeshBuilder;
  raw: MeshBuilder;
}

/** Height of the bottom of the seat countersink (top of the port throat), world y. */
export function coneBottomY(L: ModelTLayout): number {
  const V = L.valve;
  return L.block.yTop - (V.seatOuterRadius - V.throatRadius) * Math.tan(V.seatAngle);
}

/** Frame with ring axes e1 = ẑ × T, e2 = ẑ (ring placement shared with planar-face holes). */
function frameZ(c: V3, T: V3, r: number): TubeFrame {
  return { c, e1: [-T[1], T[0], 0], e2: [0, 0, 1], r };
}

/** Exhaust port: down from the seat throat, quarter bend toward the manifold face, straight out. */
export function exhaustPortFrames(L: ModelTLayout, v: ModelTValveLayout): TubeFrame[] {
  const K = L.block;
  const y0 = coneBottomY(L);
  const yc = K.portY;
  const Rb = y0 - 0.003 - yc;
  const pts: V3[] = [[v.x, y0, v.z], [v.x, yc + Rb, v.z]];
  const tans: V3[] = [[0, -1, 0], [0, -1, 0]];
  const nb = 12;
  for (let i = 1; i <= nb; i++) {
    const a = (i / nb) * (Math.PI / 2);
    pts.push([v.x - Rb + Rb * Math.cos(a), yc + Rb - Rb * Math.sin(a), v.z]);
    tans.push([-Math.sin(a), -Math.cos(a), 0]);
  }
  const x1 = v.x - Rb;
  const ns = 4;
  for (let i = 1; i <= ns; i++) {
    pts.push([x1 + ((K.xR - x1) * i) / ns, yc, v.z]);
    tans.push([-1, 0, 0]);
  }
  const s: number[] = [0];
  for (let i = 1; i < pts.length; i++) s.push(s[i - 1] + Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]));
  const total = s[s.length - 1];
  const r0 = L.valve.throatRadius, r1 = K.portRadius;
  return pts.map((p, i) => {
    const w = s[i] / total;
    const sm = w * w * (3 - 2 * w);
    return frameZ(p, tans[i], r0 + (r1 - r0) * sm);
  });
}

/** The closed block shell, split by material (exterior paint / machined deck + bores / as-cast cavities). */
export function buildBlockShell(L: ModelTLayout): ShellBuilders {
  const K = L.block;
  const paint = new MeshBuilder(), machined = new MeshBuilder(), raw = new MeshBuilder();
  const { xR, xLu, xLb, yFlareTop: yFT, yFlareBot: yFB, zF, zB, yTop, yBot } = K;
  const rb = L.boreRadius;
  const V = L.valve;
  const yCone = coneBottomY(L);

  // ---- deck: bores + seat countersinks ----
  const deckHoles = [
    ...L.axisZ.map((z) => ringXZ(0, z, rb, BORE_SEG)),
    ...L.valves.map((v) => ringXZ(v.x, v.z, V.seatOuterRadius, PORT_SEG)),
  ];
  planarFace(machined, rect(xR, zB, xLu, zF), deckHoles, [0, yTop, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  for (const z of L.axisZ) lathe(machined, [{ r: rb, y: yTop }, { r: rb, y: K.cc.yTop }], { segments: BORE_SEG, center: [0, 0, z] });

  // ---- seats, exhaust ports, siamesed intake galleries ----
  const faceHoles: [number, number][][] = [];
  const toYZ = (f: TubeFrame) => ringPoints(f, PORT_SEG).map((p) => [p[1], p[2]] as [number, number]);
  for (const v of L.valves) {
    lathe(machined, [{ r: V.seatOuterRadius, y: yTop }, { r: V.throatRadius, y: yCone }], { segments: PORT_SEG, center: [v.x, 0, v.z] });
    if (v.kind === 'exhaust') {
      const fr = exhaustPortFrames(L, v);
      sweepTube(raw, fr, PORT_SEG, true);
      faceHoles.push(toYZ(fr[fr.length - 1]));
    }
  }
  for (const gl of K.galleries) {
    const top = K.portRoofY, bot = K.portFloorY;
    const holes: P2[][] = [];
    for (const k of gl.valves) {
      const v = L.valves[k];
      const fr = [frameZ([v.x, yCone, v.z], [0, -1, 0], V.throatRadius), frameZ([v.x, top, v.z], [0, -1, 0], V.throatRadius)];
      sweepTube(raw, fr, PORT_SEG, true);
      holes.push(ringXZ(v.x, v.z, V.throatRadius, PORT_SEG));
    }
    planarFace(raw, rect(gl.x0, gl.z0, gl.x1, gl.z1), holes, [0, top, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
    planarFace(raw, rect(gl.x0, gl.z0, gl.x1, gl.z1), [], [0, bot, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
    planarFace(raw, rect(bot, gl.z0, top, gl.z1), [], [gl.x1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
    const out = [frameZ([gl.x0, K.portY, gl.zOutlet], [-1, 0, 0], K.portRadius), frameZ([xR, K.portY, gl.zOutlet], [-1, 0, 0], K.portRadius)];
    planarFace(raw, rect(bot, gl.z0, top, gl.z1), [toYZ(out[0])], [gl.x0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
    planarFace(raw, rect(gl.x0, bot, gl.x1, top), [], [0, 0, gl.z0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
    planarFace(raw, rect(gl.x0, bot, gl.x1, top), [], [0, 0, gl.z1], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
    sweepTube(raw, out, PORT_SEG, true);
    faceHoles.push(toYZ(out[1]));
  }

  // ---- exterior ----
  const vc = K.vc;
  planarFace(paint, rect(yBot, zB, yTop, zF), [...faceHoles, rect(vc.y0, vc.z0, vc.y1, vc.z1)], [xR, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(paint, rect(yFT, zB, yTop, zF), [], [xLu, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(paint, rect(yBot, zB, yFB, zF), [], [xLb, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  {
    const dx = xLb - xLu, dy = yFB - yFT;
    const l = Math.hypot(dx, dy);
    const n: V3 = [-dy / l, dx / l, 0];
    paint.quad([xLu, yFT, zB], [xLu, yFT, zF], [xLb, yFB, zF], [xLb, yFB, zB], n, n, n, n);
  }
  const section: P2[] = [[xR, yBot], [xLb, yBot], [xLb, yFB], [xLu, yFT], [xLu, yTop], [xR, yTop]];
  planarFace(paint, section, [], [0, 0, zF], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(paint, section, [], [0, 0, zB], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
  const cc = K.cc;
  planarFace(paint, rect(xR, zB, xLb, zF), [rect(cc.x0, cc.zLo0, cc.x1, cc.zLo1)], [0, yBot, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);

  // ---- crankcase cavity (open at the pan rail), stepped in z ----
  {
    const prof: P2[] = [[cc.zLo0, yBot], [cc.zLo1, yBot], [cc.zLo1, cc.yStep], [cc.zUp1, cc.yStep], [cc.zUp1, cc.yTop], [cc.zUp0, cc.yTop], [cc.zUp0, cc.yStep], [cc.zLo0, cc.yStep]];
    planarFace(raw, prof, [], [cc.x0, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0]);
    planarFace(raw, prof, [], [cc.x1, 0, 0], [0, 0, 1], [0, 1, 0], [-1, 0, 0]);
    const band = (z: number, y0: number, y1: number, nz: number) =>
      planarFace(raw, rect(cc.x0, y0, cc.x1, y1), [], [0, 0, z], [1, 0, 0], [0, 1, 0], [0, 0, nz]);
    band(cc.zLo1, yBot, cc.yStep, -1);
    band(cc.zUp1, cc.yStep, cc.yTop, -1);
    band(cc.zUp0, cc.yStep, cc.yTop, 1);
    band(cc.zLo0, yBot, cc.yStep, 1);
    planarFace(raw, rect(cc.x0, cc.zLo1, cc.x1, cc.zUp1), [], [0, cc.yStep, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
    planarFace(raw, rect(cc.x0, cc.zUp0, cc.x1, cc.zLo0), [], [0, cc.yStep, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
    planarFace(raw, rect(cc.x0, cc.zUp0, cc.x1, cc.zUp1), L.axisZ.map((z) => ringXZ(0, z, rb, BORE_SEG)), [0, cc.yTop, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  }

  // ---- valve chamber (open to the valve door on the right face) ----
  planarFace(raw, rect(xR, vc.z0, vc.x1, vc.z1), [], [0, vc.y1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(raw, rect(xR, vc.z0, vc.x1, vc.z1), [], [0, vc.y0, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(raw, rect(vc.y0, vc.z0, vc.y1, vc.z1), [], [vc.x1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(raw, rect(xR, vc.y0, vc.x1, vc.y1), [], [0, 0, vc.z0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(raw, rect(xR, vc.y0, vc.x1, vc.y1), [], [0, 0, vc.z1], [1, 0, 0], [0, 1, 0], [0, 0, -1]);

  // ---- water jacket in front of cylinder 1 (shows in the section) ----
  const j = K.jacket;
  boxCavity(raw, [j.x0, j.y0, j.z0], [j.x1, j.y1, j.z1]);
  return { paint, machined, raw };
}

/** Bearing block with a half-bore notch on its y = 0 edge (main caps below, webs above the split line). */
function notchedBlock(w: number, yFar: number, rB: number, z0: number, z1: number, seg = 32): THREE.BufferGeometry {
  const s = Math.sign(yFar);
  const pts: P2[] = [[-w, yFar], [w, yFar], [w, 0], [rB, 0]];
  for (let k = 1; k < seg; k++) {
    const a = (s * k * Math.PI) / seg; // notch on the block's side of the split line
    pts.push([rB * Math.cos(a), rB * Math.sin(a)]);
  }
  pts.push([-rB, 0], [-w, 0]);
  return extrudeZ(polyShape(pts), z0, z1, 0.0008);
}

/** U-shaped pan section (outer half-width xo, bottom yb, wall t, corner radius rc), as one polygon. */
function panSection(xo: number, yb: number, t: number, rc: number): { u: P2[]; plate: P2[] } {
  const arc = (cx: number, cy: number, r: number, a0: number, a1: number, n = 10): P2[] => {
    const o: P2[] = [];
    for (let k = 0; k <= n; k++) {
      const a = a0 + ((a1 - a0) * k) / n;
      o.push([cx + r * Math.cos(a), cy + r * Math.sin(a)]);
    }
    return o;
  };
  const outer: P2[] = [[-xo, 0], ...arc(-xo + rc, yb + rc, rc, Math.PI, 1.5 * Math.PI), ...arc(xo - rc, yb + rc, rc, 1.5 * Math.PI, 2 * Math.PI), [xo, 0]];
  const inner: P2[] = [[xo - t, 0], ...arc(xo - rc, yb + rc, rc - t, 2 * Math.PI, 1.5 * Math.PI), ...arc(-xo + rc, yb + rc, rc - t, 1.5 * Math.PI, Math.PI), [-xo + t, 0]];
  return { u: [...outer, ...inner], plate: outer };
}

/** Everything of the block assembly except the head, added to the block/head cut-frame set. */
export function buildBlock(L: ModelTLayout, set: FrameSet): void {
  const K = L.block;
  const sh = buildBlockShell(L);
  set.add('paint', sh.paint.build());
  set.add('ironMachined', sh.machined.build());
  set.add('castRaw', sh.raw.build());
  const C = L.crank;
  const J = C.journalRadius;
  const babbitt = L.rod.babbitt;
  const rB = J + babbitt;

  // ---- main bearings: babbitt shells, caps below the split line, centre web ----
  for (const m of C.mains) {
    const z0 = m.z - m.length / 2, z1 = m.z + m.length / 2;
    set.add('bearing', ringSolidZ(z0, z1, J + 0.0001, rB, 0, 0, 40), { section: 'light' });
    set.add('paint', notchedBlock(0.034, -(J + 0.6 * 0.0254), rB, z0 + 0.001, z1 - 0.001));
    const nuts = new MeshBuilder();
    for (const sx of [-1, 1]) lathe(nuts, [{ r: 0, y: 0 }, { r: 0.0075, y: 0 }, { r: 0.0075, y: 0.007 }, { r: 0, y: 0.007 }], { segments: 6, center: [sx * 0.025, -(J + 0.6 * 0.0254) - 0.007, m.z] });
    set.add('steelDark', nuts.build());
  }
  {
    const cm = C.mains[1];
    set.add('paint', notchedBlock(0.034, 0.03, rB, cm.z - cm.length / 2 + 0.001, cm.z + cm.length / 2 - 0.001));
    const web = new MeshBuilder();
    box(web, [K.cc.x0 + 0.0005, 0.024, -0.011], [K.cc.x1 - 0.0005, K.cc.yTop - 0.0005, 0.011]);
    set.add('paint', web.build());
  }
  // ---- cam bearings (split cast-iron front/centre, pressed rear bushing [FS par. 478]) ----
  for (const cb of L.cam.bearings) {
    const z0 = cb.z - cb.length / 2, z1 = cb.z + cb.length / 2;
    set.add('bearing', ringSolidZ(z0, z1, L.cam.journalRadius + 0.0001, L.cam.journalRadius + 0.0025, L.cam.x, L.cam.y, 32), { section: 'light' });
  }
  {
    const cb = L.cam.bearings[1];
    const rr = L.cam.journalRadius + 0.0025;
    set.add('paint', ringSolidZ(cb.z - cb.length / 2, cb.z + cb.length / 2, rr, rr + 0.007, L.cam.x, L.cam.y, 32));
    const arm = new MeshBuilder();
    box(arm, [K.cc.x0 + 0.0005, L.cam.y - 0.008, cb.z - cb.length / 2], [L.cam.x - rr - 0.004, L.cam.y + 0.008, cb.z + cb.length / 2]);
    set.add('paint', arm.build());
  }

  // ---- pressed-steel pan (lower crankcase) [Page18; AF15] ----
  {
    const xo = K.xLb + 0.002, t = 0.0025, yb = -0.105, rc = 0.035; // UNVERIFIED depth / radius
    const ps = panSection(xo, yb, t, rc);
    set.add('paint', extrudeZ(polyShape(ps.u), K.zB, K.zF - t));
    set.add('paint', extrudeZ(polyShape(ps.plate), K.zF - t, K.zF));
    const fl = xo + 0.012;
    set.add('paint', extrudePlanY(rect(-fl, K.zB, fl, K.zF), [rect(-xo + t, K.zB + 0.004, xo - t, K.zF - 0.004)], -0.003, 0));
  }

  // ---- timing-gear cover with the generator-pinion bulge, timer sleeve, crank-nose boss ----
  {
    const G = L.gears;
    const m = G.module;
    const circ = (dr: number) => [
      { x: 0, y: 0, r: G.crankPitchR + m + 0.006 + dr },
      { x: L.cam.x, y: L.cam.y, r: G.camPitchR + m + 0.006 + dr },
      { x: G.generator.x, y: G.generator.y, r: G.generator.pitchR + m + 0.006 + dr },
    ];
    const outer = hullOfCircles(circ(0));
    const inner = hullOfCircles(circ(-0.004));
    const zc1 = G.z1 + 0.006;
    set.add('paint', extrudeZ(polyShape(outer, [inner]), K.zF, zc1 - 0.004, 0.001));
    const circle = (cx: number, cy: number, r: number): P2[] => Array.from({ length: 40 }, (_, k) => [cx + r * Math.cos((k * 2 * Math.PI) / 40), cy + r * Math.sin((k * 2 * Math.PI) / 40)] as P2);
    set.add('paint', extrudeZ(polyShape(outer, [circle(0, 0, C.noseRadius + 0.0015), circle(L.cam.x, L.cam.y, 0.0058)]), zc1 - 0.004, zc1, 0.001));
    set.add('paint', ringSolidZ(zc1 - 0.001, L.timer.z0, 0.0058, 0.012, L.cam.x, L.cam.y, 32));
    set.add('paint', ringSolidZ(zc1 - 0.001, zc1 + 0.003, C.noseRadius + 0.0015, C.noseRadius + 0.011, 0, 0, 40));
  }

  // ---- valve door (one pressed-steel cover, 1921+ [McC]) ----
  {
    const vc = K.vc;
    const d = new MeshBuilder();
    box(d, [K.xR - 0.0025, vc.y0 - 0.008, vc.z0 - 0.01], [K.xR + 0.0002, vc.y1 + 0.008, vc.z1 + 0.01]);
    box(d, [K.xR - 0.0045, vc.y0 + 0.004, vc.z0 + 0.004], [K.xR - 0.0024, vc.y1 - 0.004, vc.z1 - 0.004]);
    set.add('paint', d.build());
    const nuts = new MeshBuilder();
    for (const z of [vc.z0 * 0.45, vc.z1 * 0.45]) {
      const g = latheZ([{ r: 0, y: 0 }, { r: 0.008, y: 0 }, { r: 0.008, y: 0.006 }, { r: 0, y: 0.006 }], { segments: 6 });
      g.rotateY(-Math.PI / 2);
      g.translate(K.xR - 0.0045, (vc.y0 + vc.y1) / 2, z);
      nuts.append(g);
      g.dispose();
    }
    set.add('steelDark', nuts.build());
  }

  // ---- copper-asbestos head gasket [Page18 pp. 258–259] ----
  {
    const H = L.head;
    set.add('copper', extrudePlanY(rect(H.xR, H.zB, H.xL, H.zF), L.chambers.map((c) => c.union), L.deckY, L.headBottomY), { section: 'light' });
  }

  // ---- generator (1919+ starter cars), driven off the large time gear [McC 1917-20] ----
  {
    const G = L.generator;
    const R = G.radius;
    const prof = [
      { r: 0, y: G.z0 - 0.012 }, { r: R * 0.6, y: G.z0 - 0.012 }, { r: R * 0.92, y: G.z0 - 0.002, smooth: true }, { r: R, y: G.z0 + 0.01 },
      { r: R, y: G.z1 - 0.03 }, { r: R + 0.004, y: G.z1 - 0.028 }, { r: R + 0.004, y: G.z1 - 0.004 }, { r: R * 0.7, y: G.z1 }, { r: 0, y: G.z1 },
    ];
    set.add('motor', latheZ(prof, { segments: 40, x: G.x, y: G.y }), { cut: false, section: false });
    const tb = new MeshBuilder();
    box(tb, [G.x - 0.016, G.y + R - 0.004, G.z0 + 0.05], [G.x + 0.016, G.y + R + 0.022, G.z0 + 0.1]);
    cylinderBetween(tb, [G.x, G.y + R + 0.022, G.z0 + 0.075], [G.x, G.y + R + 0.032, G.z0 + 0.075], 0.004, 12);
    set.add('motor', tb.build(), { cut: false, section: false });
  }
}
