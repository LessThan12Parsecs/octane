/**
 * Side-valve train of the Model T: the two-piece valve standing head-up in the block (seat line at its
 * local origin, lift +y), the spring hanging from the valve-chamber roof onto a cup pinned through the
 * stem, the mushroom push rod (tappet) riding its lobe, the camshaft with eight lobes sharing one
 * three-arc outline per cam (each rotated by its phase), the 48-tooth spiral gear and the timer rotor on
 * the nose, and the timer (commutator) case turned by the spark lever.
 * Sizes not taken from layout.ts (fillets, chamfers, bosses, caps, fastener heads) are render-only
 * cosmetic choices — UNVERIFIED.
 */
import * as THREE from 'three';
import { MeshBuilder, type ProfilePoint, type TubeFrame, type V3, cross, cylinderBetween, lathe, norm, sweepTube, merge } from '../engine/geometry';
import { buildLobeProfile, type LobeProfile } from '../engine/cam-profile';
import { lobeGeometry } from '../engine/parts-valvetrain';
import type { ModelTLayout } from './layout';
import type { ModelTMaterials } from './materials';
import type { ValveCam } from './cam';
import { helicalGear, latheZ } from './geom';

/** Valve group: origin at the seat line on the valve axis (+y up); position.y = seat line + lift. */
export function buildValve(L: ModelTLayout, M: ModelTMaterials, mat: THREE.Material): THREE.Group {
  const V = L.valve;
  const g = new THREE.Group();
  const tanA = Math.tan(V.seatAngle);
  const top = V.topY;
  const yMarginBot = top - V.margin;
  const yFaceBot = yMarginBot - (V.headRadius - V.faceInnerRadius) * tanA;
  const yBoss = top - V.bossHeight;
  // CCW (r right, y up) from the stem end: outward normals for the head-up valve
  const prof: ProfilePoint[] = [
    { r: 0, y: V.stemEndY },
    { r: V.stemRadius - 0.0007, y: V.stemEndY },
    { r: V.stemRadius, y: V.stemEndY + 0.0007 },
    { r: V.stemRadius, y: yBoss },
  ];
  const nt = 8;
  for (let i = 1; i <= nt; i++) {
    const t = i / nt;
    prof.push({ r: V.stemRadius + (V.faceInnerRadius - V.stemRadius) * Math.pow(t, 2.2), y: yBoss + (yFaceBot - yBoss) * t, smooth: i < nt });
  }
  prof.push(
    { r: V.headRadius, y: yMarginBot },
    { r: V.headRadius, y: top - 0.0004 },
    { r: V.headRadius - 0.0004, y: top },
    { r: 0, y: top },
  );
  const b = new MeshBuilder();
  lathe(b, prof, { segments: 36 });
  const valve = new THREE.Mesh(b.build(), mat);
  valve.name = 'valve-body';
  g.add(valve);
  // spring cup and its pin through the stem (move with the valve)
  // the cup rests on the pin; its top is the spring seat (layout.spring.bottomRestY at zero lift)
  const pinR = L.spring.bottomRestY - L.cup.thickness - (L.seatLineY + V.pinY);
  const y0 = V.pinY + pinR, y1 = y0 + L.cup.thickness, y2 = y1 + 0.0015;
  const cup = new MeshBuilder();
  lathe(cup, [
    { r: V.stemRadius + 0.0002, y: y0 },
    { r: L.cup.radius, y: y0 },
    { r: L.cup.radius, y: y2 },
    { r: L.cup.radius - 0.0012, y: y2 },
    { r: L.cup.radius - 0.0012, y: y1 },
    { r: V.stemRadius + 0.0002, y: y1 },
  ], { segments: 32, closed: true });
  cylinderBetween(cup, [-(V.stemRadius + 0.003), V.pinY, 0], [V.stemRadius + 0.003, V.pinY, 0], pinR, 8);
  g.add(new THREE.Mesh(cup.build(), M.m.steelForged));
  return g;
}

/** Helical spring hanging from its fixed top (origin) down `installed`; scale.y = length / installed. */
export function buildSpring(L: ModelTLayout, M: ModelTMaterials): THREE.Mesh {
  const S = L.spring;
  const sp = new MeshBuilder();
  const turns = S.turns;
  const n = Math.round(turns * 18);
  const frames: TubeFrame[] = [];
  const wr = S.wireRadius, Rm = S.meanRadius, Ls = S.installed;
  const e = 1.5 / turns; // last 1-1/2 coils close-wound [T431]
  const dead = 1.5 * 2.05 * wr;
  const span = Ls - 2 * wr;
  for (let i = 0; i <= n; i++) {
    const u = i / n;
    const rise = u < e ? (u / e) * dead : u > 1 - e ? span - ((1 - u) / e) * dead : dead + ((u - e) / (1 - 2 * e)) * (span - 2 * dead);
    const y = -(wr + rise);
    const t = u * turns * 2 * Math.PI;
    const c: V3 = [Rm * Math.cos(t), y, Rm * Math.sin(t)];
    const tan = norm([-Rm * Math.sin(t), -(Ls / turns) / (2 * Math.PI), Rm * Math.cos(t)]);
    const e1: V3 = [Math.cos(t), 0, Math.sin(t)];
    frames.push({ c, e1, e2: norm(cross(tan, e1)), r: wr });
  }
  sweepTube(sp, frames, 7, false);
  const m = new THREE.Mesh(sp.build(), M.m.steelDark);
  m.name = 'valve-spring';
  return m;
}

/** Mushroom push rod: origin at the flat face, +y up. */
export function buildTappet(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const T = L.tappet;
  const g = new THREE.Group();
  const b = new MeshBuilder();
  lathe(b, [
    { r: 0, y: 0 }, { r: T.faceRadius - 0.0005, y: 0 }, { r: T.faceRadius, y: 0.0005 }, { r: T.faceRadius, y: T.footThickness - 0.0012 },
    { r: T.faceRadius - 0.002, y: T.footThickness, smooth: true }, { r: T.stemRadius + 0.001, y: T.footThickness + 0.0025, smooth: true },
    { r: T.stemRadius, y: T.footThickness + 0.005 }, { r: T.stemRadius, y: T.length - 0.0006 }, { r: T.stemRadius - 0.0006, y: T.length }, { r: 0, y: T.length },
  ], { segments: 32 });
  const m = new THREE.Mesh(b.build(), M.m.steelPolished);
  m.name = 'tappet-body';
  g.add(m);
  return g;
}

/** Lobe outline of a valve cam (built once per distinct cam; nose along +y at rotation 0). */
export function lobeProfileOf(cam: ValveCam): LobeProfile {
  const rn = cam.geometry ? cam.geometry.noseRadius : 0.003;
  return buildLobeProfile((th) => cam.tappetLiftFromCentre(th), cam.baseRadius, 360, 0.5 * rn);
}

export interface CamshaftParts {
  /** Rotates about +z (rotation.z = θ/2) at (cam.x, cam.y). */
  group: THREE.Group;
  /** One mesh per valve (layout order), rotated by its lobe phase inside the group. */
  lobes: THREE.Mesh[];
  /** Outline used by each lobe (shared between lobes of identical cams). */
  profiles: LobeProfile[];
}

export function buildCamshaft(L: ModelTLayout, M: ModelTMaterials): CamshaftParts {
  const C = L.cam;
  const g = new THREE.Group();
  g.name = 'camshaft';
  g.position.set(C.x, C.y, 0);
  const shaftG = latheZ([
    { r: 0, y: C.z0 }, { r: C.journalRadius - 0.0006, y: C.z0 }, { r: C.journalRadius, y: C.z0 + 0.0006 },
    { r: C.journalRadius, y: C.z1 - 0.012 }, { r: 0.0058, y: C.z1 - 0.010 }, { r: 0.0058, y: C.z1 }, { r: 0, y: C.z1 },
  ], { segments: 32 });
  g.add(Object.assign(new THREE.Mesh(shaftG, M.m.steelPolished), { name: 'cam-shaft' }));
  // lobes: geometry shared per distinct cam
  // identical geometric lobes (the Ford intake and exhaust cams are the same part) share one outline
  const geoms = new Map<unknown, { geo: THREE.BufferGeometry; prof: LobeProfile }>();
  const keyOf = (c: ValveCam): unknown => {
    const g = c.geometry;
    return g ? `${g.baseRadius}|${g.flankRadius}|${g.noseRadius}|${g.rise}` : c;
  };
  const lobes: THREE.Mesh[] = [];
  const profiles: LobeProfile[] = [];
  for (const v of L.valves) {
    const key = keyOf(v.cam);
    let e = geoms.get(key);
    if (!e) {
      const prof = lobeProfileOf(v.cam);
      e = { geo: lobeGeometry(prof, C.lobeWidth, 0), prof };
      geoms.set(key, e);
    }
    const m = new THREE.Mesh(e.geo, M.m.steelPolished);
    m.name = `lobe-${v.index}`;
    m.position.z = v.z;
    m.rotation.z = v.lobePhase;
    g.add(m);
    lobes.push(m);
    profiles.push(e.prof);
  }
  // 48-tooth spiral cam gear on its flange [Good22 p. 25]: a tooth gap faces the crank at θ = 0
  const G = L.gears;
  const gamma = Math.atan2(C.y, C.x);
  const gear = helicalGear(G.camTeeth, G.camPitchR, G.module, G.z0, G.z1, G.helixAngle, -1, gamma + Math.PI + Math.PI / G.camTeeth, C.journalRadius);
  g.add(Object.assign(new THREE.Mesh(gear, M.m.steelForged), { name: 'cam-gear' }));
  const hub = new MeshBuilder();
  lathe(hub, [{ r: 0, y: G.z0 - 0.006 }, { r: C.flangeRadius, y: G.z0 - 0.006 }, { r: C.flangeRadius, y: G.z0 + 0.0005 }, { r: 0, y: G.z0 + 0.0005 }], { segments: 32 });
  lathe(hub, [{ r: 0, y: G.z1 - 0.0005 }, { r: 0.0115, y: G.z1 - 0.0005 }, { r: 0.0115, y: G.z1 + 0.008 }, { r: 0, y: G.z1 + 0.008 }], { segments: 6 });
  for (let k = 0; k < 4; k++) {
    const a = gamma + Math.PI / 4 + (k * Math.PI) / 2;
    const rr = G.camPitchR * 0.55;
    lathe(hub, [{ r: 0, y: G.z1 - 0.001 }, { r: 0.006, y: G.z1 - 0.001 }, { r: 0.006, y: G.z1 + 0.004 }, { r: 0, y: G.z1 + 0.004 }], { segments: 6, center: [rr * Math.cos(a), 0, -rr * Math.sin(a)] });
  }
  const hg = hub.build();
  hg.rotateX(Math.PI / 2);
  g.add(new THREE.Mesh(hg, M.m.steelDark));
  // timer rotor on the nose: arm + roller, along cam-local +y (cylinder 1's contact at θ = −advance)
  const T = L.timer;
  const rot = new MeshBuilder();
  const zr = (T.z0 + T.z1) / 2;
  lathe(rot, [{ r: 0, y: T.z0 + 0.004 }, { r: 0.009, y: T.z0 + 0.004 }, { r: 0.009, y: T.z1 - 0.008 }, { r: 0, y: T.z1 - 0.008 }], { segments: 16 });
  const rg = rot.build();
  rg.rotateX(Math.PI / 2);
  const arm = new MeshBuilder();
  cylinderBetween(arm, [0, 0, zr], [0, T.rotorRadius - 0.004, zr], 0.0035, 10);
  cylinderBetween(arm, [0, T.rotorRadius - 0.003, zr - 0.006], [0, T.rotorRadius - 0.003, zr + 0.006], 0.0045, 16);
  const rotor = new THREE.Mesh(merge([rg, arm.build()]), M.m.steelPolished);
  rotor.name = 'timer-rotor';
  g.add(rotor);
  return { group: g, lobes, profiles };
}

/**
 * Timer (commutator) case at the cam axis: a shallow cup (shown without its dust cover) with a fibre
 * race and four contact segments, terminal posts and the spark-lever arm. rotation.z = −advance/2.
 * Segment i spans cam-local angles [90° + offset_i/2, + contactArc/2] in the case frame.
 */
export function buildTimerCase(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const T = L.timer;
  const g = new THREE.Group();
  g.name = 'timer-case';
  g.position.set(T.x, T.y, 0);
  const R = T.radius;
  const caseG = latheZ([
    { r: 0.0062, y: T.z0 }, { r: R, y: T.z0 }, { r: R, y: T.z1 }, { r: R - 0.0025, y: T.z1 }, { r: R - 0.0025, y: T.z0 + 0.003 }, { r: 0.0062, y: T.z0 + 0.003 },
  ], { segments: 48, closed: true });
  g.add(Object.assign(new THREE.Mesh(caseG, M.x.brass), { name: 'timer-shell' }));
  const rIn = T.rotorRadius + 0.0012;
  const race = latheZ([{ r: rIn, y: T.z0 + 0.008 }, { r: R - 0.0024, y: T.z0 + 0.008 }, { r: R - 0.0024, y: T.z1 - 0.006 }, { r: rIn, y: T.z1 - 0.006 }], { segments: 48, closed: true });
  g.add(new THREE.Mesh(race, M.x.fibre));
  // contact segments (latheZ maps lathe angle φ to world angle −φ)
  const segs: THREE.BufferGeometry[] = [];
  const arc = (T.contactArcDeg / 2) * (Math.PI / 180);
  const posts = new MeshBuilder();
  L.firingOffsetDeg.forEach((off) => {
    const lead = Math.PI / 2 + (off / 2) * (Math.PI / 180);
    segs.push(latheZ([{ r: rIn - 0.0003, y: T.z0 + 0.0095 }, { r: rIn + 0.0016, y: T.z0 + 0.0095 }, { r: rIn + 0.0016, y: T.z1 - 0.0075 }, { r: rIn - 0.0003, y: T.z1 - 0.0075 }], { segments: 10, closed: true, phiStart: -(lead + arc), phiLength: arc }));
    const mid = lead + arc / 2;
    const c = Math.cos(mid), s = Math.sin(mid);
    const zc = (T.z0 + T.z1) / 2;
    cylinderBetween(posts, [c * (R - 0.001), s * (R - 0.001), zc], [c * (R + 0.012), s * (R + 0.012), zc], 0.0028, 10);
    cylinderBetween(posts, [c * (R + 0.004), s * (R + 0.004), zc], [c * (R + 0.007), s * (R + 0.007), zc], 0.0055, 6);
  });
  g.add(new THREE.Mesh(merge(segs), M.x.brass));
  for (const s of segs) s.dispose();
  g.add(new THREE.Mesh(posts.build(), M.m.steelDark));
  // spark-lever arm (to the spark rod)
  const arm = new MeshBuilder();
  const a0 = (150 * Math.PI) / 180;
  cylinderBetween(arm, [Math.cos(a0) * R, Math.sin(a0) * R, T.z0 + 0.008], [Math.cos(a0) * (R + 0.045), Math.sin(a0) * (R + 0.045), T.z0 + 0.008], 0.0032, 10);
  cylinderBetween(arm, [Math.cos(a0) * (R + 0.045), Math.sin(a0) * (R + 0.045), T.z0 + 0.002], [Math.cos(a0) * (R + 0.045), Math.sin(a0) * (R + 0.045), T.z0 + 0.016], 0.0042, 10);
  g.add(new THREE.Mesh(arm.build(), M.m.steelForged));
  return g;
}
