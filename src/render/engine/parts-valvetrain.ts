/**
 * ROOT-frame moving valve-train and drive parts: camshaft (lobes rebuilt
 * from the learned lift profile), tappets, telescoping pushrods, the CR
 * worm wheel + worm + handwheel, and the motor pulley.
 */
import * as THREE from 'three';
import type { EngineLayout } from './layout';
import type { EngineMaterials } from './materials';
import { MeshBuilder, type V3, type TubeFrame, clipGeometryZ, cross, crease, cylinderBetween, extrude, gearShape, lathe, norm, sweepTube } from './geometry';
import type { LobeProfile } from './cam-profile';

export function lobeGeometry(p: LobeProfile, width: number, z: number): THREE.BufferGeometry {
  const pts: THREE.Vector2[] = [];
  for (let i = 0; i < p.outline.length; i += 2) pts.push(new THREE.Vector2(p.outline[i], p.outline[i + 1]));
  const shape = new THREE.Shape(pts);
  const g = extrude(shape, { depth: width, bevel: 0.0006, bevelSegments: 1, creaseAngle: (40 * Math.PI) / 180 });
  g.translate(0, 0, z - width / 2);
  return g;
}

export interface CamParts {
  /** Rotates about +z at (camX, camY). */
  group: THREE.Group;
  lobes: [THREE.Mesh, THREE.Mesh];
}

export function buildCamshaft(L: EngineLayout, M: EngineMaterials, profiles: [LobeProfile, LobeProfile]): CamParts {
  const g = new THREE.Group();
  g.name = 'camshaft';
  g.position.set(L.camX, L.camY, 0);
  const shaft = new MeshBuilder();
  const z0 = L.camZBack, z1 = L.camZFront + 0.002;
  lathe(shaft, [
    { r: 0, y: z0 }, { r: L.camShaftRadius - 0.0006, y: z0 }, { r: L.camShaftRadius, y: z0 + 0.0006 },
    { r: L.camShaftRadius, y: z1 - 0.0006 }, { r: L.camShaftRadius - 0.0006, y: z1 }, { r: 0, y: z1 },
  ], { segments: 40 });
  const sg = shaft.build();
  sg.rotateX(Math.PI / 2);
  g.add(new THREE.Mesh(sg, M.m.steelPolished));
  const mkLobe = (i: 0 | 1): THREE.Mesh => {
    const v = L.valves[i];
    const m: THREE.Mesh = new THREE.Mesh(lobeGeometry(profiles[i], L.lobeWidth, v.lobeZ), M.m.steelPolished);
    m.name = `${v.key}-lobe`;
    g.add(m);
    return m;
  };
  const lobes: [THREE.Mesh, THREE.Mesh] = [mkLobe(0), mkLobe(1)];
  // cam gear: a gap faces the crank at θ = 0 (half-pitch offset from the crank's tooth)
  const gamma = Math.atan2(L.camY, L.camX);
  const phase = gamma + Math.PI + Math.PI / L.camGearTeeth;
  const gear = extrude(gearShape(L.camGearTeeth, L.camGearRadius, { phase, holeR: 0 }), {
    depth: L.gearWidth, bevel: 0.0006, bevelSegments: 1, curveSegments: 4, creaseAngle: 0.5,
  });
  gear.translate(0, 0, L.gearZ - L.gearWidth / 2);
  g.add(Object.assign(new THREE.Mesh(gear, M.m.steelForged), { name: 'cam-gear' }));
  // gear web detail: hub + 4 lightening bosses so rotation reads
  const hub = new MeshBuilder();
  lathe(hub, [
    { r: 0, y: L.gearZ - L.gearWidth / 2 - 0.004 }, { r: L.camShaftRadius + 0.01, y: L.gearZ - L.gearWidth / 2 - 0.004 },
    { r: L.camShaftRadius + 0.01, y: L.gearZ + L.gearWidth / 2 + 0.004 }, { r: 0, y: L.gearZ + L.gearWidth / 2 + 0.004 },
  ], { segments: 32 });
  for (let k = 0; k < 4; k++) {
    const a = gamma + (k * Math.PI) / 2 + Math.PI / 4;
    const rr = L.camGearRadius * 0.55;
    lathe(hub, [
      { r: 0, y: L.gearZ + L.gearWidth / 2 }, { r: 0.006, y: L.gearZ + L.gearWidth / 2 }, { r: 0.006, y: L.gearZ + L.gearWidth / 2 + 0.005 }, { r: 0, y: L.gearZ + L.gearWidth / 2 + 0.005 },
    ], { segments: 6, center: [rr * Math.cos(a), 0, -rr * Math.sin(a)] });
  }
  const hg = hub.build();
  hg.rotateX(Math.PI / 2);
  g.add(new THREE.Mesh(crease(hg, 0.5), M.m.steelDark));
  return { group: g, lobes };
}

/** Flat-faced (mushroom) tappet; origin at the face, +y up. */
export function buildTappet(L: EngineLayout, M: EngineMaterials, faceRadius: number): THREE.Group {
  const g = new THREE.Group();
  g.name = 'tappet';
  const b = new MeshBuilder();
  const rb = L.tappetRadius, rf = faceRadius, h = L.tappetLength;
  lathe(b, [
    { r: 0, y: 0 }, { r: rf - 0.0005, y: 0 }, { r: rf, y: 0.0005 }, { r: rf, y: 0.004 }, { r: rb, y: 0.0065, smooth: true },
    { r: rb, y: h - 0.0006 }, { r: rb - 0.0006, y: h }, { r: 0.0052, y: h }, { r: 0.0048, y: h - 0.0035, smooth: true }, { r: 0, y: h - 0.0045 },
  ], { segments: 40 });
  g.add(new THREE.Mesh(b.build(), M.m.steelPolished));
  return g;
}

/** Straight rod along +y from the origin, length `len`, with ball ends. */
function rodY(len: number, r: number, balls: boolean, M: EngineMaterials, mat: THREE.Material): THREE.Group {
  const g = new THREE.Group();
  const b = new MeshBuilder();
  cylinderBetween(b, [0, 0, 0], [0, len, 0], r, 20);
  g.add(new THREE.Mesh(b.build(), mat));
  if (balls) {
    const s = new THREE.SphereGeometry(r * 1.35, 20, 12);
    const b0 = new THREE.Mesh(s, M.m.steelPolished);
    g.add(b0);
  }
  return g;
}

export interface PushrodParts {
  /** Lower rod: origin at the tappet ball, +y along the pushrod. */
  lower: THREE.Group;
  /** Upper tube: origin at the rocker cup, +y pointing DOWN the pushrod. */
  upper: THREE.Group;
}

export function buildPushrod(L: EngineLayout, M: EngineMaterials, lowerLen: number, upperLen: number): PushrodParts {
  const lower = rodY(lowerLen, L.pushrodRadius, true, M, M.m.steelPolished);
  lower.name = 'pushrod-lower';
  const upper = new THREE.Group();
  upper.name = 'pushrod-upper';
  const b = new MeshBuilder();
  const R = L.pushrodTubeRadius, ri = L.pushrodRadius + 0.0002;
  lathe(b, [
    { r: 0, y: 0 }, { r: R * 0.9, y: 0 }, { r: R, y: 0.003 }, { r: R, y: upperLen - 0.001 }, { r: R - 0.001, y: upperLen }, { r: ri, y: upperLen }, { r: ri, y: upperLen - 0.02 }, { r: 0, y: upperLen - 0.02 },
  ], { segments: 24 });
  upper.add(new THREE.Mesh(b.build(), M.m.steelForged));
  return { lower, upper };
}

export interface WormParts {
  /** Worm wheel for the full view (rotates about +y). */
  wheel: THREE.Mesh;
  /** Half wheel for the cutaway (rotated by angle mod tooth pitch). */
  wheelCut: THREE.Group;
  /** Worm + handwheel; rotates about +z. */
  worm: THREE.Group;
  toothPitch: number;
}

export function buildWormDrive(L: EngineLayout, M: EngineMaterials): WormParts {
  const w = L.wormWheel;
  const shape = gearShape(w.teeth, w.pitchRadius, { holeR: w.innerRadius });
  const geo = extrude(shape, { depth: w.thickness, bevel: 0.001, bevelSegments: 1, curveSegments: 2, creaseAngle: 0.5 });
  geo.rotateX(-Math.PI / 2); // extrusion along +y
  geo.translate(0, w.yBottom, 0);
  const wheel = new THREE.Mesh(geo, M.m.bronze);
  wheel.name = 'cr-worm-wheel';
  const wheelCut = new THREE.Group();
  wheelCut.name = 'cr-worm-wheel-cut';
  const cg = clipGeometryZ(geo, 0);
  wheelCut.add(new THREE.Mesh(cg, M.m.bronze));
  const sec = new THREE.Mesh(cg, M.sectionMaterial('light', { kind: 'half' }));
  sec.name = 'section';
  sec.renderOrder = 1;
  wheelCut.add(sec);

  const worm = new THREE.Group();
  worm.name = 'cr-worm';
  worm.position.set(L.worm.x, L.worm.y, 0);
  const wr = L.worm.radius;
  const z0 = L.worm.zBack, z1 = L.worm.zFront;
  const sb = new MeshBuilder();
  lathe(sb, [
    { r: 0, y: z0 - 0.01 }, { r: wr * 0.55, y: z0 - 0.01 }, { r: wr * 0.55, y: z1 }, { r: 0, y: z1 },
  ], { segments: 24 });
  lathe(sb, [
    { r: 0, y: -0.022 }, { r: wr * 0.62, y: -0.022 }, { r: wr * 0.62, y: 0.022 }, { r: 0, y: 0.022 },
  ], { segments: 24 });
  const sgeo = sb.build();
  sgeo.rotateX(Math.PI / 2);
  worm.add(new THREE.Mesh(sgeo, M.m.steelPolished));
  // single-start thread, pitch = wheel tooth pitch, around z = 0 (mesh point)
  const pitch = (2 * Math.PI * w.pitchRadius) / w.teeth;
  const thr = new MeshBuilder();
  const turns = 0.044 / pitch;
  const n = Math.round(turns * 32);
  const frames: TubeFrame[] = [];
  const R = wr * 0.72;
  for (let i = 0; i <= n; i++) {
    const t = (i / n - 0.5) * turns * 2 * Math.PI;
    // worm axis = z; helix in local x-y
    const c: V3 = [R * Math.cos(t), R * Math.sin(t), (t / (2 * Math.PI)) * pitch];
    const tan = norm([-R * Math.sin(t), R * Math.cos(t), pitch / (2 * Math.PI)]);
    const e1: V3 = [Math.cos(t), Math.sin(t), 0];
    frames.push({ c, e1, e2: norm(cross(tan, e1)), r: wr * 0.2 });
  }
  sweepTube(thr, frames, 8, false);
  worm.add(new THREE.Mesh(thr.build(), M.m.steelPolished));
  // handwheel at the front end
  const hw = new THREE.Group();
  hw.position.set(0, 0, z1);
  const rim = new THREE.Mesh(new THREE.TorusGeometry(0.05, 0.0055, 12, 48), M.m.steelForged);
  hw.add(rim);
  const spokes = new MeshBuilder();
  for (let k = 0; k < 3; k++) {
    const a = (k * 2 * Math.PI) / 3;
    cylinderBetween(spokes, [0, 0, 0], [0.05 * Math.cos(a), 0.05 * Math.sin(a), 0], 0.0032, 10);
  }
  lathe(spokes, [{ r: 0, y: -0.008 }, { r: 0.011, y: -0.008 }, { r: 0.011, y: 0.01 }, { r: 0, y: 0.01 }], { segments: 20 });
  const spg = spokes.build();
  hw.add(new THREE.Mesh(spg, M.m.steelForged));
  const knob = new MeshBuilder();
  cylinderBetween(knob, [0.05, 0, 0], [0.05, 0, 0.04], 0.006, 16);
  hw.add(new THREE.Mesh(knob.build(), M.m.rubber));
  worm.add(hw);
  return { wheel, wheelCut, worm, toothPitch: (2 * Math.PI) / w.teeth };
}

/** Motor pulley (rotates about +z at the motor axis). */
export function buildMotorPulley(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'motor-pulley';
  g.position.set(L.motor.x, L.motor.y, 0);
  const R = L.motor.pulleyRadius + 0.006;
  const W = L.flywheelWidth;
  const pitch = W / 3;
  const prof: { r: number; y: number }[] = [{ r: 0.014, y: L.flywheelZ - W / 2 }, { r: R, y: L.flywheelZ - W / 2 }];
  for (let k = 1; k <= 2; k++) {
    const c = L.flywheelZ - W / 2 + k * pitch;
    prof.push({ r: R, y: c - 0.0075 }, { r: R - 0.009, y: c - 0.0026 }, { r: R - 0.009, y: c + 0.0026 }, { r: R, y: c + 0.0075 });
  }
  prof.push({ r: R, y: L.flywheelZ + W / 2 }, { r: 0.014, y: L.flywheelZ + W / 2 });
  const b = new MeshBuilder();
  lathe(b, prof, { segments: 48 });
  const geo = b.build();
  geo.rotateX(Math.PI / 2);
  g.add(new THREE.Mesh(geo, M.m.ironMachined));
  // set screws / bolt heads so rotation reads
  const bb = new MeshBuilder();
  for (let k = 0; k < 3; k++) {
    const a = (k * 2 * Math.PI) / 3;
    const rr = R * 0.55;
    lathe(bb, [{ r: 0, y: 0 }, { r: 0.005, y: 0 }, { r: 0.005, y: 0.004 }, { r: 0, y: 0.004 }], { segments: 6, center: [rr * Math.cos(a), 0, -rr * Math.sin(a)] });
  }
  const bgeo = bb.build();
  bgeo.rotateX(Math.PI / 2);
  bgeo.translate(0, 0, L.flywheelZ + W / 2);
  g.add(new THREE.Mesh(crease(bgeo, 0.5), M.m.steelDark));
  return g;
}
