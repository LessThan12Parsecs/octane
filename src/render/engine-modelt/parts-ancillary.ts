/**
 * Ancillaries of the Model T engine: cast-iron Y intake manifold and exhaust log on the valve side,
 * updraft carburettor with its butterfly, flat-belt two-blade fan on a bracket at the front, starter on
 * the transmission cover, and the bell-shaped transmission cover with the stationary 16-coil magneto
 * ring inside. World coordinates; dimensions in layout.ts (positions UNVERIFIED, render only).
 */
import * as THREE from 'three';
import { MeshBuilder, type V3, box, cylinderBetween, lathe, merge } from '../engine/geometry';
import type { ModelTLayout } from './layout';
import type { ModelTMaterials } from './materials';
import type { FrameSet } from './frame-set';
import { beltPath, flatBelt, latheZ, pathFrames, ringSolidY, ringSolidZ, tubeSolid } from './geom';
import { spurGear } from './parts-crank';

/** Cubic Bézier polyline. */
function bezier(p0: V3, p1: V3, p2: V3, p3: V3, n: number): V3[] {
  const out: V3[] = [];
  for (let i = 0; i <= n; i++) {
    const t = i / n, u = 1 - t;
    const a = u * u * u, b = 3 * u * u * t, c = 3 * u * t * t, d = t * t * t;
    out.push([a * p0[0] + b * p1[0] + c * p2[0] + d * p3[0], a * p0[1] + b * p1[1] + c * p2[1] + d * p3[1], a * p0[2] + b * p1[2] + c * p2[2] + d * p3[2]]);
  }
  return out;
}

/** Closed ring around an axis parallel to x through (y, z), from x0 to x1. */
function ringSolidX(x0: number, x1: number, rIn: number, rOut: number, y: number, z: number, seg = 32): THREE.BufferGeometry {
  const g = ringSolidZ(x0, x1, rIn, rOut, 0, 0, seg);
  g.rotateY(Math.PI / 2); // z → x (z' = −x, x' = z)
  g.translate(0, y, z);
  return g;
}

export function buildManifolds(L: ModelTLayout, set: FrameSet): void {
  const K = L.block;
  const I = L.intake, E = L.exhaust;
  // ---- intake: carburettor flange → riser → T → two branches to the siamesed ports ----
  {
    const o = new MeshBuilder(), i = new MeshBuilder();
    tubeSolid(o, i, pathFrames([[I.x, I.flangeY + 0.004, 0], [I.x, I.teeY + 0.004, 0]], I.radius, [0, 0, 1]), I.wall, 28);
    for (const p of K.ports.filter((q) => q.kind === 'intake')) {
      const pts = bezier([I.x, I.teeY, 0], [I.x, I.teeY + 0.012, 0.55 * p.z], [K.xR - 0.04, K.portY, p.z], [K.xR - 0.0015, K.portY, p.z], 22);
      tubeSolid(o, i, pathFrames(pts, I.radius, [0, 1, 0]), I.wall, 28);
      o.append(ringSolidX(K.xR - 0.011, K.xR - 0.0005, I.radius, I.radius + I.wall + 0.006, K.portY, p.z, 28));
    }
    const jb = new MeshBuilder();
    const rj = I.radius + I.wall + 0.003;
    const prof = Array.from({ length: 13 }, (_, k) => {
      const a = -Math.PI / 2 + (k * Math.PI) / 12;
      return { r: rj * Math.cos(a), y: I.teeY + rj * Math.sin(a), smooth: k > 0 && k < 12 };
    });
    lathe(jb, prof, { segments: 28, center: [I.x, 0, 0] });
    o.append(jb.build());
    o.append(ringSolidY(I.flangeY - 0.002, I.flangeY + 0.006, I.radius, I.radius + I.wall + 0.012, I.x, 0, 28));
    set.add('paint', o.build());
    set.add('castRaw', i.build());
  }
  // ---- exhaust log above, stubs from the four exhaust ports, outlet at the rear turned down ----
  {
    const o = new MeshBuilder(), cav = new MeshBuilder();
    const Ro = E.radius + E.wall;
    lathe(o, [{ r: 0, y: E.z0 }, { r: Ro - 0.004, y: E.z0 }, { r: Ro, y: E.z0 + 0.004, smooth: true }, { r: Ro, y: E.z1 - 0.004, smooth: true }, { r: Ro - 0.004, y: E.z1 }, { r: 0, y: E.z1 }], { segments: 32 });
    lathe(cav, [{ r: 0, y: E.z0 + E.wall }, { r: 0, y: E.z1 - E.wall }, { r: E.radius, y: E.z1 - E.wall }, { r: E.radius, y: E.z0 + E.wall }], { segments: 32, closed: true });
    const log = o.build(), logCav = cav.build();
    for (const g of [log, logCav]) {
      g.rotateX(Math.PI / 2);
      g.translate(E.x, E.y, 0);
    }
    set.add('castRaw', log);
    set.add('castRaw', logCav);
    const so = new MeshBuilder(), si = new MeshBuilder();
    for (const p of K.ports.filter((q) => q.kind === 'exhaust')) {
      const pts = bezier([K.xR - 0.0015, K.portY, p.z], [K.xR - 0.03, K.portY, p.z], [E.x, E.y - 0.03, p.z], [E.x, E.y - 0.004, p.z], 18);
      tubeSolid(so, si, pathFrames(pts, K.portRadius, [0, 0, 1]), 0.0045, 28);
      so.append(ringSolidX(K.xR - 0.011, K.xR - 0.0005, K.portRadius, K.portRadius + 0.0105, K.portY, p.z, 28));
    }
    // outlet: down from the rear end of the log, then back toward the exhaust pipe
    const zo = E.z0 + 0.025;
    const pts: V3[] = [[E.x, E.y - 0.01, zo], ...bezier([E.x, E.y - 0.02, zo], [E.x, E.outletY - 0.01, zo], [E.x, E.outletY - 0.02, zo - 0.01], [E.x, E.outletY - 0.03, zo - 0.07], 14)];
    tubeSolid(so, si, pathFrames(pts, 0.017, [0, 0, 1]), 0.005, 28);
    set.add('castRaw', so.build());
    set.add('castRaw', si.build());
  }
}

export interface CarbParts {
  body: THREE.Group;
  /** Butterfly plate + throttle lever; rotation.z = plate angle ψ (0 = across the bore). */
  throttle: THREE.Group;
}

/** Updraft carburettor below the intake riser (Holley NH / Kingston L-4 class), brass body. */
export function buildCarburettor(L: ModelTLayout, M: ModelTMaterials): CarbParts {
  const C = L.carb;
  const body = new THREE.Group();
  body.name = 'carburettor';
  const y1 = C.topY, y0 = C.topY - C.bodyHeight;
  const parts: THREE.BufferGeometry[] = [];
  parts.push(ringSolidY(y0, y1, C.throatRadius, C.bodyRadius, C.x, C.z, 32));
  parts.push(ringSolidY(y1 - 0.007, y1, C.throatRadius, C.bodyRadius + 0.012, C.x, C.z, 32));
  // float bowl beside the throat
  const bx = C.x + C.bowlOffset[0], by = C.topY + C.bowlOffset[1], bz = C.z + C.bowlOffset[2];
  const bowl = new MeshBuilder();
  lathe(bowl, [
    { r: 0, y: by - 0.03 }, { r: C.bowlRadius - 0.004, y: by - 0.03 }, { r: C.bowlRadius, y: by - 0.026, smooth: true }, { r: C.bowlRadius, y: by + 0.022 },
    { r: C.bowlRadius - 0.006, y: by + 0.026 }, { r: 0, y: by + 0.026 },
  ], { segments: 32, center: [bx, 0, bz] });
  parts.push(bowl.build());
  const link = new MeshBuilder();
  box(link, [C.x - 0.012, by - 0.012, C.z], [C.x + 0.012, by + 0.012, bz - C.bowlRadius + 0.004]);
  parts.push(link.build());
  // air horn turning forward at the bottom (hot-air pipe from the exhaust stove fits here)
  const ho = new MeshBuilder(), hi = new MeshBuilder();
  const hp: V3[] = [[C.x, y0 + 0.002, C.z], ...Array.from({ length: 9 }, (_, k) => {
    const a = ((k + 1) / 9) * (Math.PI / 2);
    const R = 0.022;
    return [C.x, y0 - R * Math.sin(a), C.z - R + R * Math.cos(a)] as V3;
  })];
  hp.push([C.x, y0 - 0.022, C.z - 0.06]);
  tubeSolid(ho, hi, pathFrames(hp, C.throatRadius, [0, 0, 1]), 0.003, 28);
  parts.push(ho.build(), hi.build());
  body.add(new THREE.Mesh(merge(parts), M.x.brass));
  for (const p of parts) p.dispose();

  // butterfly: plate in the throat + lever outside, both on the shaft (axis along z)
  const throttle = new THREE.Group();
  throttle.name = 'throttle';
  throttle.position.set(C.x, y1 - 0.016, C.z);
  const plate = new MeshBuilder();
  lathe(plate, [{ r: 0, y: -0.0008 }, { r: C.throatRadius - 0.0004, y: -0.0008 }, { r: C.throatRadius - 0.0004, y: 0.0008 }, { r: 0, y: 0.0008 }], { segments: 28 });
  const plateG = plate.build();
  const lever = new MeshBuilder();
  cylinderBetween(lever, [0, 0, -C.bodyRadius - 0.012], [0, 0, C.bodyRadius + 0.006], 0.0025, 10);
  box(lever, [-0.004, -0.0025, C.bodyRadius + 0.006], [0.032, 0.0025, C.bodyRadius + 0.011]);
  cylinderBetween(lever, [0.029, 0, C.bodyRadius + 0.011], [0.029, 0, C.bodyRadius + 0.02], 0.0022, 8);
  throttle.add(new THREE.Mesh(plateG, M.x.brass), new THREE.Mesh(lever.build(), M.m.steelForged));
  return { body, throttle };
}

export interface FanParts {
  /** Rotates about +z at (fan.x, fan.y): rotation.z = −θ·ratio. */
  fan: THREE.Group;
  /** Static belt (the bracket is a housing part: buildFanBracket). */
  belt: THREE.Mesh;
}

export function buildFan(L: ModelTLayout, M: ModelTMaterials): FanParts {
  const F = L.fan;
  const fan = new THREE.Group();
  fan.name = 'fan';
  fan.position.set(F.x, F.y, 0);
  const zb = L.belt.z;
  const hub: THREE.BufferGeometry[] = [];
  hub.push(latheZ([{ r: 0, y: F.z - 0.012 }, { r: F.hubRadius, y: F.z - 0.012 }, { r: F.hubRadius, y: F.z + 0.008 }, { r: F.hubRadius * 0.5, y: F.z + 0.014, smooth: true }, { r: 0, y: F.z + 0.016 }], { segments: 32 }));
  hub.push(latheZ([{ r: 0.008, y: zb - L.belt.width / 2 - 0.003 }, { r: F.pulleyRadius, y: zb - L.belt.width / 2 - 0.003 }, { r: F.pulleyRadius - 0.001, y: zb }, { r: F.pulleyRadius, y: zb + L.belt.width / 2 + 0.003 }, { r: 0.008, y: zb + L.belt.width / 2 + 0.003 }], { segments: 40, closed: true }));
  hub.push(latheZ([{ r: 0, y: zb - 0.016 }, { r: 0.0075, y: zb - 0.016 }, { r: 0.0075, y: F.z }, { r: 0, y: F.z }], { segments: 16 }));
  fan.add(new THREE.Mesh(merge(hub), M.m.steelForged));
  for (const h of hub) h.dispose();
  // two pressed-steel blades with ≈ 28° pitch (UNVERIFIED)
  const blades = new MeshBuilder();
  for (let k = 0; k < F.blades; k++) {
    const tmp = new MeshBuilder();
    box(tmp, [-0.022, F.hubRadius * 0.6, -0.0012], [0.022, F.radius, 0.0012]);
    const g = tmp.build();
    g.rotateY((28 * Math.PI) / 180);
    g.rotateZ((k * 2 * Math.PI) / F.blades);
    g.translate(0, 0, F.z - 0.004);
    blades.append(g);
    g.dispose();
  }
  fan.add(Object.assign(new THREE.Mesh(blades.build(), M.m.paint), { name: 'fan-blades' }));

  // flat belt (crank pulley ↔ fan pulley)
  const bm = new MeshBuilder();
  flatBelt(bm, beltPath([0, 0], L.crank.pulleyRadius + 0.0005, [F.x, F.y], F.pulleyRadius + 0.0005), zb, L.belt.width, L.belt.thickness);
  const belt = new THREE.Mesh(bm.build(), M.x.leather);
  belt.name = 'fan-belt';
  return { fan, belt };
}

/** Fan bracket from the cylinder front to the fan spindle (static housing part, not cut). */
export function buildFanBracket(L: ModelTLayout, set: FrameSet): void {
  const F = L.fan;
  const zb = L.belt.z;
  const br = new MeshBuilder();
  cylinderBetween(br, [F.x, F.y, zb - 0.016], [F.x, F.y, zb - 0.032], 0.012, 20);
  box(br, [F.x - 0.009, L.deckY - 0.03, L.block.zF - 0.001], [F.x + 0.009, F.y + 0.006, zb - 0.026]);
  set.add('paint', br.build(), { cut: false, section: false });
}

/** Starter (1919+) on the left of the transmission cover, Bendix pinion at the ring gear. */
export function buildStarter(L: ModelTLayout, set: FrameSet): void {
  const S = L.starter;
  const R = S.radius;
  set.add('motor', latheZ([
    { r: 0, y: S.z0 }, { r: R * 0.85, y: S.z0 }, { r: R, y: S.z0 + 0.008, smooth: true }, { r: R, y: S.z1 - 0.012 }, { r: R * 0.92, y: S.z1 - 0.004, smooth: true },
    { r: R * 0.5, y: S.z1 }, { r: 0, y: S.z1 },
  ], { segments: 40, x: S.x, y: S.y }), { cut: false, section: false });
  set.add('motor', ringSolidZ(S.pinionZ - 0.01, S.z0 + 0.002, 0, R * 0.42, S.x, S.y, 28), { cut: false, section: false });
  set.add('steelForged', spurGear(S.pinionTeeth, (2 * S.pinionR) / S.pinionTeeth, S.pinionZ - 0.0095, S.pinionZ + 0.0095, S.x, S.y), { cut: false, section: false });
  const tb = new MeshBuilder();
  box(tb, [S.x - 0.014, S.y + R - 0.004, S.z1 - 0.07], [S.x + 0.014, S.y + R + 0.016, S.z1 - 0.035]);
  set.add('motor', tb.build(), { cut: false, section: false });
}

/** Transmission cover / flywheel housing (closed lathe shell) + the stationary magneto coil ring. */
export function buildBell(L: ModelTLayout, set: FrameSet): void {
  const B = L.bell;
  const inner = B.inner.map((p) => ({ r: p.r, y: p.z }));
  inner.push({ r: 0, y: inner[inner.length - 1].y });
  // outward offset of the inner profile
  const outer = inner.map((p, k) => {
    const a = inner[Math.max(k - 1, 0)], b = inner[Math.min(k + 1, inner.length - 1)];
    const dr = b.r - a.r, dz = b.y - a.y;
    const l = Math.hypot(dr, dz) || 1;
    return { r: p.r + (-dz / l) * B.wall, y: p.y + (dr / l) * B.wall };
  });
  outer[0] = { r: inner[0].r + B.wall, y: inner[0].y };
  outer[outer.length - 1] = { r: 0, y: inner[inner.length - 1].y - B.wall };
  const loop = [...outer.slice().reverse(), ...inner];
  const shell = latheZ(loop, { segments: 96, closed: true });
  set.add('paint', shell);
  // front wall (closes the bell against the block / pan)
  const zf = L.block.zB;
  set.add('paint', ringSolidZ(zf - B.wall, zf, B.frontInnerRadius, inner[0].r + B.wall, 0, 0, 96));
  // magneto coil ring [Dyke24 pp. 248–249]: plate + 16 spools (copper windings on iron cores)
  const C = L.coilRing;
  set.add('paint', ringSolidZ(C.plateZ0, C.plateZ1, C.rIn, C.rOut, 0, 0, 96));
  const wind: THREE.BufferGeometry[] = [], cores: THREE.BufferGeometry[] = [];
  for (let k = 0; k < C.count; k++) {
    const a = ((k + 0.5) * 2 * Math.PI) / C.count;
    const x = C.coilR * Math.cos(a), y = C.coilR * Math.sin(a);
    wind.push(ringSolidZ(C.coilZ0 + 0.002, C.plateZ0, C.coilRadius * 0.45, C.coilRadius, x, y, 20));
    cores.push(ringSolidZ(C.coilZ0, C.plateZ0, 0, C.coilRadius * 0.45, x, y, 14));
  }
  set.add('copper', merge(wind), { section: 'light' });
  set.add('steelDark', merge(cores));
  for (const g of [...wind, ...cores]) g.dispose();
}
