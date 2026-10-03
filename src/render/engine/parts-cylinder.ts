/**
 * Cylinder assembly (cylinder frame; moves with the compression ratio):
 * barrel with water jacket and threaded spigot, flat head with seat inserts,
 * valve guides, intake/exhaust ports and pushrod passages, rocker pedestals
 * and cover, port pipes; plus the moving valve train parts carried by the
 * head (valves, springs, retainers, rockers) and the spark plug.
 */
import * as THREE from 'three';
import type { EngineLayout, ValveLayout } from './layout';
import type { EngineMaterials } from './materials';
import type { StaticSet } from './static-set';
import {
  MeshBuilder, type V3, type TubeFrame, lathe, planarFace, rect, ringXZ, ringPoints, sweepTube,
  cylinderBetween, box, extrude, crease, cross, norm, sub, add, scale, dot,
} from './geometry';

const SEG = 64; // must be even (port/pocket vertex matching)

export interface ValveParts {
  /** Valve + retainer; translate y by −lift. */
  moving: THREE.Group;
  /** Spring mesh; bottom at springSeatY, scale.y = length / installed length. */
  spring: THREE.Mesh;
}

export interface CylinderParts {
  valves: [ValveParts, ValveParts];
  /** Rocker groups positioned at the pivot; set rotation (α about local x after yaw). */
  rockers: [THREE.Group, THREE.Group];
  plug: THREE.Group;
}

/** Plan x of the pushrod where it passes through the head (rest pose). */
function pushrodXAt(L: EngineLayout, v: ValveLayout, yCyl: number, headY: number): number {
  const top: V3 = [v.pushrodTopRest[0], v.pushrodTopRest[1], v.pushrodTopRest[2]];
  const botY = L.camY + L.camBaseRadius + L.tappetLength - headY;
  const t = (yCyl - botY) / (top[1] - botY);
  return v.tappetX + (top[0] - v.tappetX) * t;
}

/** Frames along the port centreline: up from the seat, bend, out through the side wall. */
function portFrames(L: EngineLayout, v: ValveLayout): TubeFrame[] {
  const dir = v.portDir;
  const y0 = v.insertHeight;
  const Rb = Math.max(v.portY - y0, 0.004);
  const P0: V3 = [v.x, y0, v.z];
  const P3: V3 = [v.x + dir * Rb, v.portY, v.z];
  const k = 0.5523 * Rb; // quarter circle
  const P1: V3 = [v.x, y0 + k, v.z];
  const P2: V3 = [v.x + dir * (Rb - k), v.portY, v.z];
  const exitX = dir * L.headXHalf;
  const bn: V3 = [0, 0, 1];
  const frames: TubeFrame[] = [];
  const pts: V3[] = [];
  const tans: V3[] = [];
  const nb = 18;
  for (let i = 0; i <= nb; i++) {
    const t = i / nb, u = 1 - t;
    const p: V3 = [
      u * u * u * P0[0] + 3 * u * u * t * P1[0] + 3 * u * t * t * P2[0] + t * t * t * P3[0],
      u * u * u * P0[1] + 3 * u * u * t * P1[1] + 3 * u * t * t * P2[1] + t * t * t * P3[1],
      v.z,
    ];
    const d: V3 = [
      3 * u * u * (P1[0] - P0[0]) + 6 * u * t * (P2[0] - P1[0]) + 3 * t * t * (P3[0] - P2[0]),
      3 * u * u * (P1[1] - P0[1]) + 6 * u * t * (P2[1] - P1[1]) + 3 * t * t * (P3[1] - P2[1]),
      0,
    ];
    pts.push(p);
    tans.push(norm(d));
  }
  const nStraight = 5;
  for (let i = 1; i <= nStraight; i++) {
    const t = i / nStraight;
    pts.push([P3[0] + (exitX - P3[0]) * t, v.portY, v.z]);
    tans.push([dir, 0, 0]);
  }
  // exact end tangents
  tans[0] = [0, 1, 0];
  // cumulative length for the radius blend
  const s: number[] = [0];
  for (let i = 1; i < pts.length; i++) s.push(s[i - 1] + Math.hypot(pts[i][0] - pts[i - 1][0], pts[i][1] - pts[i - 1][1]));
  const total = s[s.length - 1];
  for (let i = 0; i < pts.length; i++) {
    const T = tans[i];
    const e1 = cross(bn, T);
    const w = s[i] / total;
    const sm = w * w * (3 - 2 * w);
    frames.push({ c: pts[i], e1, e2: bn, r: v.seatInnerRadius + (v.portExitRadius - v.seatInnerRadius) * sm });
  }
  return frames;
}

function hollowBoxOpenBottom(
  outer: MeshBuilder, inner: MeshBuilder,
  x0: number, x1: number, y0: number, y1: number, z0: number, z1: number, w: number,
): void {
  // outer skin
  planarFace(outer, rect(x0, z0, x1, z1), [], [0, y1, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(outer, rect(y0, z0, y1, z1), [], [x0, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(outer, rect(y0, z0, y1, z1), [], [x1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(outer, rect(x0, y0, x1, y1), [], [0, 0, z0], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
  planarFace(outer, rect(x0, y0, x1, y1), [], [0, 0, z1], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  // bottom rim
  planarFace(outer, rect(x0, z0, x1, z1), [rect(x0 + w, z0 + w, x1 - w, z1 - w)], [0, y0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  // inner skin (normals into the cavity)
  const X0 = x0 + w, X1 = x1 - w, Z0 = z0 + w, Z1 = z1 - w, Y1 = y1 - w;
  planarFace(inner, rect(X0, Z0, X1, Z1), [], [0, Y1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(inner, rect(y0, Z0, Y1, Z1), [], [X0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(inner, rect(y0, Z0, Y1, Z1), [], [X1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(inner, rect(X0, y0, X1, Y1), [], [0, 0, Z0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(inner, rect(X0, y0, X1, Y1), [], [0, 0, Z1], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
}

/** Tube solid (inner passage + outer skin + end annuli) along frames. */
function tubeSolid(outer: MeshBuilder, inner: MeshBuilder, frames: TubeFrame[], wall: number, seg: number): void {
  sweepTube(inner, frames, seg, true);
  const of = frames.map((f) => ({ ...f, r: f.r + wall }));
  sweepTube(outer, of, seg, false);
  for (const [i, sgn] of [[0, -1], [frames.length - 1, 1]] as [number, number][]) {
    const a = ringPoints(frames[i], seg), b = ringPoints(of[i], seg);
    const f = frames[i];
    const T = norm(cross(f.e1, f.e2));
    const n: V3 = scale(T, sgn);
    for (let k = 0; k < seg; k++) {
      const k1 = (k + 1) % seg;
      outer.quad(a[k], a[k1], b[k1], b[k], n, n, n, n);
    }
  }
}

export function buildCylinder(L: EngineLayout, M: EngineMaterials, set: StaticSet): CylinderParts {
  const T = L.headThickness;
  const X = L.headXHalf;
  const Z0 = L.headZMin, Z1 = L.headZMax;
  const rb = L.boreRadius;
  const Rj = L.jacketOuterRadius;
  const Lj = L.jacketDepth;
  const Lc = L.cylinderLength;
  const Rs = L.spigotRadius;

  // ---------------- barrel casting ----------------
  {
    const paint = new MeshBuilder();
    lathe(paint, [
      { r: rb, y: -Lc },
      { r: Rs - 0.001, y: -Lc },
      { r: Rs, y: -Lc + 0.001 },
      { r: Rs, y: -Lj - 0.003 },
      { r: Rs + 0.003, y: -Lj },
      { r: Rj - 0.003, y: -Lj },
      { r: Rj, y: -Lj + 0.003 },
      { r: Rj, y: -0.004 },
      { r: Rj - 0.004, y: 0 },
      { r: rb, y: 0 },
    ], { segments: 96 });
    set.add('paint', paint.build());
    const bore = new MeshBuilder();
    lathe(bore, [{ r: rb, y: 0 }, { r: rb, y: -Lc }], { segments: 96 });
    set.add('ironMachined', bore.build());
    // water jacket cavity (CW → normals into the cavity)
    const jc = L.jacketCavity;
    const jb = new MeshBuilder();
    lathe(jb, [
      { r: jc.rIn, y: jc.yBottom }, { r: jc.rIn, y: jc.yTop }, { r: jc.rOut, y: jc.yTop }, { r: jc.rOut, y: jc.yBottom },
    ], { segments: 96, closed: true });
    set.add('castRaw', jb.build());
    // external thread on the spigot (the worm wheel is threaded onto it)
    const th = new MeshBuilder();
    const yA = -Lc + 0.006, yB = -Lj - 0.006;
    const turns = (yB - yA) / L.threadPitch;
    const n = Math.max(8, Math.round(turns * 40));
    const frames: TubeFrame[] = [];
    const R = Rs + 0.0003;
    for (let i = 0; i <= n; i++) {
      const t = (i / n) * turns * 2 * Math.PI;
      const c: V3 = [R * Math.cos(t), yA + (yB - yA) * (i / n), R * Math.sin(t)];
      const tan = norm([-R * Math.sin(t), L.threadPitch / (2 * Math.PI), R * Math.cos(t)]);
      const e1: V3 = [Math.cos(t), 0, Math.sin(t)];
      frames.push({ c, e1, e2: norm(cross(tan, e1)), r: 0.0008 });
    }
    sweepTube(th, frames, 6, false);
    set.add('ironMachined', th.build(), { section: false });
  }

  // ---------------- head block ----------------
  const valves = L.valves;
  const passX = valves.map((v) => pushrodXAt(L, v, T / 2, L.headY(L.spec.geometry.compressionRatio)));
  const pr = L.pushrodPassageRadius;
  {
    const outline = rect(-X, Z0, X, Z1);
    const deck = new MeshBuilder();
    const deckHoles = [
      ...valves.map((v) => ringXZ(v.x, v.z, v.pocketRadius, SEG)),
      ...valves.map((v, i) => ringXZ(passX[i], v.lobeZ, pr, SEG)),
    ];
    planarFace(deck, outline, deckHoles, [0, 0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
    set.add('ironMachined', deck.build());

    const skin = new MeshBuilder();
    planarFace(skin, outline, valves.map((v, i) => ringXZ(passX[i], v.lobeZ, pr, SEG)), [0, T, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
    const pf = valves.map((v) => portFrames(L, v));
    for (const side of [1, -1] as const) {
      const holes: [number, number][][] = [];
      valves.forEach((v, i) => {
        if (v.portDir !== side) return;
        const last = pf[i][pf[i].length - 1];
        holes.push(ringPoints(last, SEG).map((p) => [p[1], p[2]] as [number, number]));
      });
      planarFace(skin, rect(0, Z0, T, Z1), holes, [side * X, 0, 0], [0, 1, 0], [0, 0, 1], [side, 0, 0]);
    }
    planarFace(skin, rect(-X, 0, X, T), [], [0, 0, Z0], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
    planarFace(skin, rect(-X, 0, X, T), [], [0, 0, Z1], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
    set.add('paint', skin.build());

    // pockets (seat-insert bores), ports, pushrod passages
    const pockets = new MeshBuilder();
    const ports = new MeshBuilder();
    valves.forEach((v, i) => {
      lathe(pockets, [
        { r: v.seatInnerRadius, y: v.insertHeight }, { r: v.pocketRadius, y: v.insertHeight }, { r: v.pocketRadius, y: 0 },
      ], { segments: SEG, center: [v.x, 0, v.z] });
      sweepTube(ports, pf[i], SEG, true);
      lathe(ports, [{ r: pr, y: T }, { r: pr, y: 0 }], { segments: SEG, center: [passX[i], 0, v.lobeZ] });
    });
    set.add('ironMachined', pockets.build());
    set.add('castRaw', ports.build());

    // seat inserts
    const ins = new MeshBuilder();
    for (const v of valves) {
      const lipR = v.headRadius + v.lipClearance;
      const ta = Math.tan(v.spec.seatAngle);
      lathe(ins, [
        { r: lipR, y: 0 },
        { r: v.pocketRadius, y: 0 },
        { r: v.pocketRadius, y: v.insertHeight },
        { r: v.seatInnerRadius, y: v.insertHeight },
        { r: v.seatInnerRadius, y: v.margin + v.faceHeight },
        { r: lipR, y: v.margin - v.lipClearance * ta },
      ], { segments: SEG, center: [v.x, 0, v.z], closed: true });
    }
    set.add('nickel', ins.build());

    // valve guides + spring seat washers
    const guides = new MeshBuilder();
    const washers = new MeshBuilder();
    for (const v of valves) {
      lathe(guides, [
        { r: v.stemRadius + 0.0003, y: v.guideBottomY },
        { r: v.guideOuterRadius, y: v.guideBottomY },
        { r: v.guideOuterRadius, y: v.guideTopY - 0.0006 },
        { r: v.guideOuterRadius - 0.0006, y: v.guideTopY },
        { r: v.stemRadius + 0.0003, y: v.guideTopY },
      ], { segments: 32, center: [v.x, 0, v.z], closed: true });
      lathe(washers, [
        { r: v.guideOuterRadius + 0.0002, y: T },
        { r: v.springMeanRadius + v.springWireRadius + 0.002, y: T },
        { r: v.springMeanRadius + v.springWireRadius + 0.002, y: v.springSeatY },
        { r: v.guideOuterRadius + 0.0002, y: v.springSeatY },
      ], { segments: 48, center: [v.x, 0, v.z], closed: true });
    }
    set.add('bearing', guides.build(), { section: 'light' });
    set.add('steelDark', washers.build());

    // head studs / nuts
    const nuts = new MeshBuilder();
    const nx = X - 0.013, nz0 = Z0 + 0.013, nz1 = Z1 - 0.013;
    for (const [x, z] of [[-nx, nz0], [nx, nz0], [-nx, nz1], [nx, nz1], [-nx, (nz0 + nz1) / 2], [nx, (nz0 + nz1) / 2]]) {
      lathe(nuts, [
        { r: 0, y: T }, { r: 0.0095, y: T }, { r: 0.0095, y: T + 0.009 }, { r: 0.0075, y: T + 0.0105 }, { r: 0, y: T + 0.0105 },
      ], { segments: 6, center: [x, 0, z], closed: false });
    }
    set.add('steelDark', crease(nuts.build(), 0.4));
  }

  // ---------------- spark-plug hole and boss ----------------
  const pl = L.plug;
  {
    const ax = pl.axis;
    const at = (d: number): V3 => sub(pl.gapCenter, scale(ax, d));
    const passage = new MeshBuilder();
    const bossStart = Math.max(pl.noseDistance, pl.seatDistance - 0.014);
    cylinderBetween(passage, at(pl.noseDistance), at(bossStart), pl.threadRadius + 0.0003, 32, { inward: true, caps: false });
    set.add('ironMachined', passage.build());
    const bo = new MeshBuilder(), bi = new MeshBuilder();
    const e1 = norm(cross(ax, Math.abs(ax[1]) < 0.9 ? [0, 1, 0] : [1, 0, 0]));
    const e2 = cross(e1, ax);
    // frames ordered from the chamber side outward; ring plane ⟂ axis
    const fr = [bossStart, pl.seatDistance].map((d): TubeFrame => ({ c: at(d), e1, e2, r: pl.threadRadius + 0.0003 }));
    tubeSolid(bo, bi, fr, pl.bossRadius - pl.threadRadius - 0.0003, 32);
    set.add('paint', bo.build());
    set.add('ironMachined', bi.build());
  }

  // ---------------- port pipes ----------------
  {
    const pipeOut = new MeshBuilder(), pipeIn = new MeshBuilder();
    const flOut = new MeshBuilder(), flIn = new MeshBuilder();
    for (const v of valves) {
      const dir = v.portDir;
      const exit = portFrames(L, v);
      const ef = exit[exit.length - 1];
      const x0 = dir * X;
      const straight = 0.03;
      const bendR = 0.035;
      const pts: V3[] = [];
      for (let i = 0; i <= 2; i++) pts.push([x0 + dir * straight * (i / 2), v.portY, v.z]);
      const cx = x0 + dir * straight, cz = v.z - bendR;
      for (let i = 1; i <= 12; i++) {
        const a = (i / 12) * (Math.PI / 2);
        pts.push([cx + dir * bendR * Math.sin(a), v.portY, cz + bendR * Math.cos(a)]);
      }
      const last = pts[pts.length - 1];
      for (let i = 1; i <= 3; i++) pts.push([last[0], v.portY, last[2] - 0.03 * i]);
      // parallel-transport the exit frame (matches the wall hole) along the path
      const fr: TubeFrame[] = [{ ...ef }];
      for (let i = 1; i < pts.length; i++) {
        const tan = norm(sub(pts[Math.min(i + 1, pts.length - 1)], pts[i - 1]));
        const prev = fr[i - 1];
        const e1 = norm(sub(prev.e1, scale(tan, dot(prev.e1, tan))));
        fr.push({ c: pts[i], e1, e2: cross(tan, e1), r: ef.r });
      }
      tubeSolid(pipeOut, pipeIn, fr, 0.0045, SEG);
      // flange at the wall
      const ffr: TubeFrame[] = [{ ...ef }, { ...ef, c: add(ef.c, [dir * 0.009, 0, 0]) }];
      tubeSolid(flOut, flIn, ffr, 0.016, SEG);
    }
    set.add('paint', pipeOut.build());
    set.add('castRaw', pipeIn.build());
    set.add('paint', flOut.build());
    set.add('castRaw', flIn.build());
  }

  // ---------------- rocker pedestals, shafts, cover ----------------
  const rockerW = 0.012;
  const hubR = 0.0095;
  const shaftR = 0.0045;
  {
    const ped = new MeshBuilder();
    const shafts = new MeshBuilder();
    for (const v of valves) {
      const [px, py, pz] = v.rockerPivot;
      const cheek = 0.008;
      const yaw = v.rockerYaw;
      const tmp = new MeshBuilder();
      for (const s of [-1, 1]) {
        const xa = s * (rockerW / 2 + 0.0008), xb = s * (rockerW / 2 + 0.0008 + cheek);
        box(tmp, [Math.min(xa, xb), T - py, -0.011], [Math.max(xa, xb), 0.008, 0.011]);
      }
      const g = tmp.build();
      g.rotateY(yaw);
      g.translate(px, py, pz);
      ped.append(g);
      g.dispose();
      const half = rockerW / 2 + 0.0008 + cheek + 0.003;
      const d: V3 = [Math.cos(yaw), 0, -Math.sin(yaw)]; // pivot axis = ŷ × u
      cylinderBetween(shafts, add([px, py, pz], scale(d, -half)), add([px, py, pz], scale(d, half)), shaftR, 24);
    }
    set.add('steelForged', ped.build());
    set.add('steelPolished', shafts.build());

    const co = new MeshBuilder(), ci = new MeshBuilder();
    const m = 0.004;
    hollowBoxOpenBottom(co, ci, -X + m, X - m, T, L.rockerCover.yTop, Z0 + m, Z1 - m, L.rockerCover.wall);
    set.add('paint', co.build());
    set.add('castRaw', ci.build());
  }

  // ---------------- moving: valves, springs, rockers ----------------
  const vparts = valves.map((v, i) => buildValve(v, i === 0 ? M.m.valve : M.m.valveExhaust, M)) as [ValveParts, ValveParts];
  const rockers = valves.map((v) => buildRocker(v, M, rockerW, hubR, shaftR)) as [THREE.Group, THREE.Group];
  const plug = buildSparkPlug(L, M);
  return { valves: vparts, rockers, plug };
}

function buildValve(v: ValveLayout, mat: THREE.Material, M: EngineMaterials): ValveParts {
  const moving = new THREE.Group();
  moving.name = `${v.key}-valve`;
  moving.position.set(v.x, 0, v.z);
  const yf = v.margin + v.faceHeight;
  const prof: { r: number; y: number; smooth?: boolean }[] = [
    { r: 0, y: 0 },
    { r: v.headRadius - 0.0004, y: 0 },
    { r: v.headRadius, y: 0.0004 },
    { r: v.headRadius, y: v.margin },
    { r: v.seatInnerRadius, y: yf },
  ];
  const nt = 10;
  for (let i = 1; i <= nt; i++) {
    const t = i / nt;
    prof.push({ r: v.stemRadius + (v.seatInnerRadius - v.stemRadius) * Math.pow(1 - t, 2.4), y: yf + v.tulipHeight * t, smooth: i < nt });
  }
  // keeper groove near the tip
  const yk = v.springSeatY + v.springInstalledLength + v.retainerThickness * 0.5;
  prof.push({ r: v.stemRadius, y: yk - 0.002 });
  prof.push({ r: v.stemRadius - 0.0007, y: yk - 0.0015 });
  prof.push({ r: v.stemRadius - 0.0007, y: yk + 0.0015 });
  prof.push({ r: v.stemRadius, y: yk + 0.002 });
  prof.push({ r: v.stemRadius, y: v.tipY - 0.0008 });
  prof.push({ r: v.stemRadius - 0.0008, y: v.tipY });
  prof.push({ r: 0, y: v.tipY });
  const b = new MeshBuilder();
  lathe(b, prof, { segments: 64 });
  moving.add(new THREE.Mesh(b.build(), mat));

  // shroud (masked arc on the port side of the valve head, standing in the throat)
  const arc = (v.spec.shroudArcDeg * Math.PI) / 180;
  if (arc > 0) {
    const sb = new MeshBuilder();
    const rOut = v.seatInnerRadius - 0.0004;
    const rIn = rOut - 0.0022;
    const h = v.spec.maxLift + 0.002;
    const phiC = v.spec.shroudDirection + Math.PI; // masked side is opposite the opening
    const phi0 = phiC - arc / 2;
    const segs = Math.max(8, Math.round((arc / (2 * Math.PI)) * 64));
    lathe(sb, [
      { r: rIn, y: yf - 0.0006 }, { r: rOut, y: yf - 0.0006 }, { r: rOut, y: yf + h }, { r: rIn, y: yf + h },
    ], { segments: segs, closed: true, phiStart: phi0, phiLength: arc });
    // end caps
    for (const a of [phi0, phi0 + arc]) {
      const c = Math.cos(a), s = Math.sin(a);
      const n: V3 = a === phi0 ? [s, 0, -c] : [-s, 0, c];
      const P = (r: number, y: number): V3 => [r * c, y, r * s];
      sb.quad(P(rIn, yf - 0.0006), P(rOut, yf - 0.0006), P(rOut, yf + h), P(rIn, yf + h), n, n, n, n);
    }
    moving.add(new THREE.Mesh(sb.build(), mat));
  }

  // retainer (moves with the valve)
  const ry = v.springSeatY + v.springInstalledLength;
  const rr = new MeshBuilder();
  lathe(rr, [
    { r: v.stemRadius + 0.0002, y: ry - 0.001 },
    { r: v.springMeanRadius - v.springWireRadius, y: ry },
    { r: v.springMeanRadius + v.springWireRadius + 0.0015, y: ry },
    { r: v.springMeanRadius + v.springWireRadius + 0.0015, y: ry + v.retainerThickness * 0.6 },
    { r: v.stemRadius + 0.0025, y: ry + v.retainerThickness },
    { r: v.stemRadius + 0.0002, y: ry + v.retainerThickness },
  ], { segments: 48, closed: true });
  moving.add(new THREE.Mesh(rr.build(), M.m.steelForged));

  // spring (static bottom, scaled in y)
  const sp = new MeshBuilder();
  const Ls = v.springInstalledLength;
  const turns = v.springTurns;
  const n = Math.round(turns * 28);
  const frames: TubeFrame[] = [];
  const wr = v.springWireRadius;
  const Rm = v.springMeanRadius;
  for (let i = 0; i <= n; i++) {
    const u = i / n;
    // closed ends: first/last 3/4 turn at wire-touching pitch
    const e = 0.75 / turns;
    const deadRise = 0.75 * 2.1 * wr;
    const span = Ls - 2 * wr;
    const rise = u < e ? (u / e) * deadRise : u > 1 - e ? span - ((1 - u) / e) * deadRise : deadRise + ((u - e) / (1 - 2 * e)) * (span - 2 * deadRise);
    const y = wr + rise;
    const t = u * turns * 2 * Math.PI;
    const c: V3 = [Rm * Math.cos(t), y, Rm * Math.sin(t)];
    const tan = norm([-Rm * Math.sin(t), (Ls / turns) / (2 * Math.PI), Rm * Math.cos(t)]);
    const e1: V3 = [Math.cos(t), 0, Math.sin(t)];
    frames.push({ c, e1, e2: norm(cross(tan, e1)), r: wr });
  }
  sweepTube(sp, frames, 10, false);
  const spring = new THREE.Mesh(sp.build(), M.m.steelDark);
  spring.name = `${v.key}-spring`;
  spring.position.set(v.x, v.springSeatY, v.z);
  return { moving, spring };
}

function buildRocker(v: ValveLayout, M: EngineMaterials, w: number, hubR: number, shaftR: number): THREE.Group {
  const g = new THREE.Group();
  g.name = `${v.key}-rocker`;
  g.position.set(v.rockerPivot[0], v.rockerPivot[1], v.rockerPivot[2]);
  g.rotation.order = 'YXZ';
  g.rotation.y = v.rockerYaw;
  const av = v.armValve, ap = v.armPushrod;
  // side profile in (u = local z, v = local y)
  const s = new THREE.Shape();
  const padR = 0.0065;
  s.moveTo(av + padR, padR);
  s.lineTo(av + 0.006, 0.011);
  s.quadraticCurveTo(av * 0.5, 0.017, 0, 0.019);
  s.quadraticCurveTo(-ap * 0.5, 0.017, -ap - 0.009, 0.013);
  s.lineTo(-ap - 0.009, 0.0035);
  s.lineTo(-ap + 0.009, 0.0035);
  s.quadraticCurveTo(-ap * 0.4, 0.002, -hubR * 0.8, -hubR * 0.6);
  s.absarc(0, 0, hubR, Math.PI + 0.7, -0.7, false);
  s.quadraticCurveTo(av * 0.4, 0.002, av - padR, padR);
  s.absarc(av, padR, padR, Math.PI, 2 * Math.PI, false);
  s.closePath();
  s.holes.push(new THREE.Path().absarc(0, 0, shaftR + 0.0002, 0, Math.PI * 2, true));
  const geo = extrude(s, { depth: w, bevel: 0.0012, bevelSegments: 2, curveSegments: 16 });
  geo.rotateY(-Math.PI / 2); // (u, v, d) → (x = −d, y = v, z = u)
  geo.translate(w / 2, 0, 0);
  g.add(new THREE.Mesh(geo, M.m.steelForged));
  // hub bushing
  const hb = new MeshBuilder();
  cylinderBetween(hb, [-w / 2 - 0.0005, 0, 0], [w / 2 + 0.0005, 0, 0], shaftR + 0.0012, 24);
  g.add(new THREE.Mesh(hb.build(), M.m.bearing));
  // adjusting screw + lock nut + ball at the pushrod end
  const sb = new MeshBuilder();
  lathe(sb, [{ r: 0, y: -0.001 }, { r: 0.0038, y: -0.001 }, { r: 0.0038, y: 0.024 }, { r: 0, y: 0.024 }], { segments: 16, center: [0, 0, -ap] });
  lathe(sb, [{ r: 0, y: 0.014 }, { r: 0.0068, y: 0.014 }, { r: 0.0068, y: 0.02 }, { r: 0.0055, y: 0.021 }, { r: 0, y: 0.021 }], { segments: 6, center: [0, 0, -ap] });
  g.add(new THREE.Mesh(crease(sb.build(), 0.4), M.m.steelDark));
  const ball = new THREE.Mesh(new THREE.SphereGeometry(0.0045, 20, 12), M.m.steelPolished);
  ball.position.set(0, 0.0025, -ap);
  g.add(ball);
  return g;
}

export function buildSparkPlug(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const pl = L.plug;
  const sp = L.spec.sparkPlug;
  const g = new THREE.Group();
  g.name = 'spark-plug';
  // basis: X = side (ground strut), Y = axis (into the chamber), Z = X × Y
  const X = new THREE.Vector3(...pl.side);
  const Y = new THREE.Vector3(...pl.axis);
  const Zb = new THREE.Vector3().crossVectors(X, Y);
  g.matrixAutoUpdate = true;
  g.quaternion.setFromRotationMatrix(new THREE.Matrix4().makeBasis(X, Y, Zb));
  g.position.set(...pl.gapCenter);

  const gap = sp.gap;
  const rce = sp.centerElectrodeDiameter / 2;
  const tr = pl.threadRadius;
  const rInShell = Math.max(tr - 0.0028, rce + 0.0025);
  const yNose = -pl.noseDistance;
  const ySeat = -pl.seatDistance;
  const ceProtrude = Math.min(0.0025, Math.max(0.0006, (pl.noseDistance - gap / 2) * 0.55));
  const yInsTip = -gap / 2 - ceProtrude;

  // centre electrode
  const ce = new MeshBuilder();
  lathe(ce, [{ r: 0, y: yInsTip - 0.006 }, { r: rce, y: yInsTip - 0.006 }, { r: rce, y: -gap / 2 }, { r: 0, y: -gap / 2 }], { segments: 24 });
  g.add(Object.assign(new THREE.Mesh(ce.build(), M.m.nickel), { name: 'center-electrode' }));

  // insulator: nose cone inside the shell + upper body with ribs
  const ins = new MeshBuilder();
  const yHexBot = ySeat - 0.0015;
  const hexH = Math.max(0.009, 0.55 * sp.threadDiameter);
  const yHexTop = yHexBot - hexH;
  const yCrimp = yHexTop - 0.004;
  const rUp = Math.max(0.0055, 0.38 * sp.threadDiameter);
  const upLen = Math.max(0.03, 2.1 * sp.threadDiameter);
  const prof: { r: number; y: number; smooth?: boolean }[] = [
    { r: 0, y: yInsTip },
    { r: rce + 0.0009, y: yInsTip },
    { r: rInShell - 0.0006, y: yNose - 0.006 },
    { r: rInShell - 0.0003, y: yCrimp + 0.002 },
    { r: rUp, y: yCrimp },
  ];
  const nr = 5;
  const ribTop = yCrimp - 0.004;
  const ribLen = upLen * 0.55;
  for (let i = 0; i < nr; i++) {
    const y0 = ribTop - (ribLen * i) / nr;
    const y1 = y0 - ribLen / nr;
    prof.push({ r: rUp, y: y0 - 0.0005, smooth: true });
    prof.push({ r: rUp - 0.0011, y: (y0 + y1) / 2, smooth: true });
    prof.push({ r: rUp, y: y1 + 0.0005, smooth: true });
  }
  prof.push({ r: rUp * 0.92, y: yCrimp - upLen });
  prof.push({ r: 0, y: yCrimp - upLen });
  lathe(ins, prof.reverse(), { segments: 48 });
  g.add(Object.assign(new THREE.Mesh(ins.build(), M.m.ceramic), { name: 'insulator' }));

  // shell: threaded body (nose → seat), hex, crimp
  const sh = new MeshBuilder();
  const shp: { r: number; y: number }[] = [{ r: rInShell, y: ySeat - 0.0015 }, { r: tr, y: ySeat - 0.0015 }, { r: tr, y: ySeat }];
  const pitch = 0.0015;
  for (let y = ySeat + pitch; y < yNose - pitch * 0.5; y += pitch) {
    shp.push({ r: tr - 0.0006, y: y - pitch / 2 });
    shp.push({ r: tr, y });
  }
  shp.push({ r: tr, y: yNose - 0.0003 });
  shp.push({ r: tr - 0.0003, y: yNose });
  shp.push({ r: rInShell, y: yNose });
  lathe(sh, shp, { segments: 48, closed: true });
  // hex
  const hexR = (0.6 * sp.threadDiameter) / Math.cos(Math.PI / 6);
  lathe(sh, [
    { r: rUp, y: yHexTop }, { r: hexR * 0.92, y: yHexTop }, { r: hexR, y: yHexTop + 0.0012 }, { r: hexR, y: yHexBot - 0.0012 },
    { r: hexR * 0.92, y: yHexBot }, { r: tr, y: yHexBot },
  ], { segments: 6, phiStart: Math.PI / 6 });
  lathe(sh, [{ r: rUp + 0.0004, y: yCrimp }, { r: hexR * 0.72, y: yHexTop }], { segments: 48 });
  g.add(Object.assign(new THREE.Mesh(crease(sh.build(), 0.6), M.m.plated), { name: 'plug-shell' }));

  // gasket washer
  const gw = new MeshBuilder();
  lathe(gw, [{ r: tr, y: yHexBot }, { r: tr + 0.0045, y: yHexBot }, { r: tr + 0.0045, y: ySeat }, { r: tr, y: ySeat }], { segments: 48, closed: true });
  g.add(new THREE.Mesh(gw.build(), M.m.copper));

  // terminal
  const tb = new MeshBuilder();
  const yT = yCrimp - upLen;
  lathe(tb, [{ r: 0, y: yT }, { r: 0.0022, y: yT }, { r: 0.0022, y: yT - 0.006 }, { r: 0.0035, y: yT - 0.007 }, { r: 0.0035, y: yT - 0.012 }, { r: 0, y: yT - 0.013 }].reverse(), { segments: 24 });
  g.add(new THREE.Mesh(tb.build(), M.m.steelPolished));

  // ground electrode: strut welded on the shell nose at +X, bent over the centre electrode
  const ge = new MeshBuilder();
  const tge = Math.max(0.0008, 0.5 * sp.groundElectrodeWidth);
  const wge = sp.groundElectrodeWidth;
  const rho = (rInShell + tr) / 2;
  box(ge, [rho - tge / 2, yNose - 0.0005, -wge / 2], [rho + tge / 2, gap / 2 + tge, wge / 2]);
  box(ge, [-(rce + 0.0009), gap / 2, -wge / 2], [rho + tge / 2, gap / 2 + tge, wge / 2]);
  g.add(Object.assign(new THREE.Mesh(ge.build(), M.m.nickel), { name: 'ground-electrode' }));
  return g;
}
