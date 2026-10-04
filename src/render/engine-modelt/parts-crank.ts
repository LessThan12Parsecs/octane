/**
 * Rotating and reciprocating parts of the Model T: the four-throw, three-main crankshaft without
 * counterweights (throws in one plane, crank-local frame: rotated by −θ about z), its 24-tooth spiral
 * timing gear and the flat-belt pulley; the flywheel with the 16 V magnets of the magneto, the 120-tooth
 * starter ring gear and the three triple-gear clusters of the planetary transmission; the I-beam rod with
 * its two-bolt babbitt big end and clamped small end; the long cast-iron piston with two rings above
 * and the oil ring below the pin.
 * Sizes not taken from layout.ts (fillets, chamfers, bosses, caps, fastener heads) are render-only
 * cosmetic choices — UNVERIFIED.
 */
import * as THREE from 'three';
import { MeshBuilder, type ProfilePoint, cylinderBetween, lathe, merge } from '../engine/geometry';
import type { ModelTLayout } from './layout';
import type { ModelTMaterials } from './materials';
import { extrudeZ, gearOutline, helicalGear, hullOfCircles, latheZ, polyShape } from './geom';
import type { P2 } from './outline';

function named<T extends THREE.Object3D>(o: T, name: string): T {
  o.name = name;
  return o;
}

/** Plain shaft along z (closed, chamfered). */
function shaft(z0: number, z1: number, r: number, x = 0, y = 0, chamfer = 0.0008, seg = 40): THREE.BufferGeometry {
  const lo = Math.min(z0, z1), hi = Math.max(z0, z1);
  const c = Math.min(chamfer, r * 0.3, (hi - lo) * 0.3);
  return latheZ([{ r: 0, y: lo }, { r: r - c, y: lo }, { r, y: lo + c }, { r, y: hi - c }, { r: r - c, y: hi }, { r: 0, y: hi }], { segments: seg, x, y });
}

/** Crank-gear tooth phase: a tooth points at the camshaft at θ = 0. */
export function crankGearPhase(L: ModelTLayout): number {
  return Math.atan2(L.cam.y, L.cam.x);
}

export function buildCrank(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'crankshaft';
  const C = L.crank;
  const r = L.crankRadius;
  const J = C.journalRadius;
  const half = C.pinLength / 2;
  const t = C.webThickness;
  const forged: THREE.BufferGeometry[] = [];
  const polished: THREE.BufferGeometry[] = [];
  for (let i = 0; i < L.nCyl; i++) {
    const a = L.throwAngle[i];
    const px = r * Math.sin(a), py = r * Math.cos(a);
    const z = L.axisZ[i];
    const web = polyShape(hullOfCircles([{ x: 0, y: 0, r: C.webBossRadius }, { x: px, y: py, r: C.webBossRadius }], 40));
    forged.push(extrudeZ(web, z + half, z + half + t, 0.0015));
    forged.push(extrudeZ(web, z - half - t, z - half, 0.0015));
    polished.push(shaft(z - half - 0.002, z + half + 0.002, C.pinRadius, px, py, 0.0005));
  }
  // journals / connecting shafts between the webs, front nose, rear flange
  const zs = L.axisZ;
  polished.push(shaft(zs[0] + half + t - 0.002, L.gears.z0, J));
  for (let i = 0; i + 1 < L.nCyl; i++) polished.push(shaft(zs[i + 1] + half + t - 0.002, zs[i] - half - t + 0.002, J));
  polished.push(shaft(C.flangeZ1 - 0.001, zs[L.nCyl - 1] - half - t + 0.002, J));
  forged.push(shaft(L.gears.z0 - 0.001, C.frontEndZ, C.noseRadius));
  forged.push(latheZ([{ r: 0, y: C.flangeZ0 }, { r: C.flangeRadius - 0.001, y: C.flangeZ0 }, { r: C.flangeRadius, y: C.flangeZ0 + 0.001 }, { r: C.flangeRadius, y: C.flangeZ1 - 0.001 }, { r: J + 0.004, y: C.flangeZ1 }, { r: 0, y: C.flangeZ1 }], { segments: 48 }));
  g.add(named(new THREE.Mesh(merge(forged), M.m.steelForged), 'crank-forging'));
  g.add(named(new THREE.Mesh(merge(polished), M.m.steelPolished), 'crank-journals'));
  for (const geo of [...forged, ...polished]) geo.dispose();

  // 24-tooth spiral timing gear [Good22 pp. 24–25]
  const G = L.gears;
  const gear = helicalGear(G.crankTeeth, G.crankPitchR, G.module, G.z0, G.z1, G.helixAngle, 1, crankGearPhase(L), C.noseRadius);
  g.add(named(new THREE.Mesh(gear, M.m.steelForged), 'crank-gear'));

  // flat-belt fan pulley with the starting-crank pin at the nose
  const pr = C.pulleyRadius;
  const w = C.pulleyZ1 - C.pulleyZ0;
  const pul = latheZ([
    { r: C.noseRadius, y: C.pulleyZ0 }, { r: pr - 0.002, y: C.pulleyZ0 }, { r: pr, y: C.pulleyZ0 + 0.002 }, { r: pr + 0.0008, y: C.pulleyZ0 + w / 2, smooth: true },
    { r: pr, y: C.pulleyZ1 - 0.002 }, { r: pr - 0.002, y: C.pulleyZ1 }, { r: pr * 0.45, y: C.pulleyZ1 }, { r: pr * 0.45, y: C.pulleyZ1 - w * 0.3 }, { r: C.noseRadius, y: C.pulleyZ1 - w * 0.3 },
  ], { segments: 56, closed: true });
  g.add(named(new THREE.Mesh(pul, M.m.ironMachined), 'crank-pulley'));
  const pin = new MeshBuilder();
  cylinderBetween(pin, [-0.016, 0, C.frontEndZ - 0.005], [0.016, 0, C.frontEndZ - 0.005], 0.0032, 12);
  g.add(new THREE.Mesh(pin.build(), M.m.steelDark));

  g.add(buildFlywheel(L, M));
  return g;
}

/** V-shaped magnet outline (apex toward the hub) in the face plane, centred at polar angle β. */
function magnetOutline(rApex: number, rPole: number, halfAngle: number, w: number, beta: number): P2[] {
  const A: P2 = [rApex, 0];
  const Pp: P2 = [rPole * Math.cos(halfAngle), rPole * Math.sin(halfAngle)];
  const Pm: P2 = [Pp[0], -Pp[1]];
  const dl = Math.hypot(Pp[0] - A[0], Pp[1] - A[1]);
  const d: P2 = [(Pp[0] - A[0]) / dl, (Pp[1] - A[1]) / dl];
  const np: P2 = [-d[1], d[0]]; // outer side of the + leg
  const nm: P2 = [-d[1], -d[0]]; // outer side of the − leg (mirror)
  const sinPhi = Math.abs(d[1]);
  const hw = w / 2;
  const pts: P2[] = [
    [rApex - hw / sinPhi, 0],
    [Pm[0] + nm[0] * hw, Pm[1] + nm[1] * hw],
    [Pm[0] - nm[0] * hw, Pm[1] - nm[1] * hw],
    [rApex + hw / sinPhi, 0],
    [Pp[0] - np[0] * hw, Pp[1] - np[1] * hw],
    [Pp[0] + np[0] * hw, Pp[1] + np[1] * hw],
  ];
  const c = Math.cos(beta), s = Math.sin(beta);
  return pts.map(([x, y]) => [x * c - y * s, x * s + y * c] as P2);
}

export function buildFlywheel(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'flywheel';
  const F = L.flywheel;
  const R = F.ringGear.pitchR - 1.25 * F.ringGear.module;
  // hub + web + disc as one closed section (r, z)
  const hubR = F.hubRadius;
  const prof: ProfilePoint[] = [
    { r: 0.0165, y: F.hubZ1 },
    { r: 0.0165, y: F.discZ0 },
    { r: R - 0.003, y: F.discZ0 },
    { r: R, y: F.discZ0 + 0.003 },
    { r: R, y: F.discZ1 - 0.003 },
    { r: R - 0.003, y: F.discZ1 },
    { r: hubR + 0.012, y: F.discZ1 },
    { r: hubR, y: F.discZ1 + 0.012, smooth: true },
    { r: hubR, y: F.hubZ1 - 0.002 },
    { r: hubR - 0.002, y: F.hubZ1 },
  ];
  g.add(named(new THREE.Mesh(latheZ(prof, { segments: 96, closed: true }), M.m.paint), 'flywheel-body'));
  // starter ring gear (120 T) on the rim
  const rg = F.ringGear;
  const ring = extrudeZ(polyShape(gearOutline(rg.teeth, rg.pitchR, rg.module), [hullOfCircles([{ x: 0, y: 0, r: R - 0.004 }], 96)]), rg.z0, rg.z1, 0.0006, 4);
  g.add(named(new THREE.Mesh(ring, M.m.steelForged), 'ring-gear'));
  // 16 V magnets on the front face [Dyke24 pp. 248–249]
  const mg = F.magnets;
  const half = (Math.PI / mg.count) * 0.68; // adjacent like poles nearly touch
  const mags: THREE.BufferGeometry[] = [];
  for (let k = 0; k < mg.count; k++) {
    const beta = (k * 2 * Math.PI) / mg.count;
    mags.push(extrudeZ(polyShape(magnetOutline(mg.rApex, mg.rPole, half, mg.legWidth, beta)), mg.z0, mg.z1, 0.0012, 2));
  }
  g.add(named(new THREE.Mesh(merge(mags), M.x.magnet), 'magnets'));
  for (const m of mags) m.dispose();
  // cap screws on the hub (4 [Good22])
  const sc = new MeshBuilder();
  for (let k = 0; k < 4; k++) {
    const a = Math.PI / 4 + (k * Math.PI) / 2;
    const rr = hubR - 0.012;
    sc.append(latheZ([{ r: 0, y: F.hubZ1 }, { r: 0.0075, y: F.hubZ1 }, { r: 0.0075, y: F.hubZ1 + 0.0065 }, { r: 0, y: F.hubZ1 + 0.0065 }], { segments: 6, x: rr * Math.cos(a), y: rr * Math.sin(a) }));
  }
  g.add(new THREE.Mesh(sc.build(), M.m.steelDark));
  // triple gears (planetary transmission) on the rear face
  const T = L.flywheel.tripleGears;
  const tg: THREE.BufferGeometry[] = [];
  const dz = (T.z1 - T.z0) / T.teeth.length;
  for (let k = 0; k < T.count; k++) {
    const a = Math.PI / 2 + (k * 2 * Math.PI) / T.count;
    const cx = T.r * Math.cos(a), cy = T.r * Math.sin(a);
    T.teeth.forEach((n, j) => {
      const geo = extrudeZ(polyShape(gearOutline(n, (n * T.module) / 2, T.module)), T.z1 - (j + 1) * dz + 0.001, T.z1 - j * dz, 0.0005, 2);
      geo.translate(cx, cy, 0);
      tg.push(geo);
    });
    tg.push(shaft(T.z0 - 0.004, F.discZ0 + 0.004, 0.0075, cx, cy));
  }
  g.add(named(new THREE.Mesh(merge(tg), M.m.steelForged), 'triple-gears'));
  for (const geo of tg) geo.dispose();
  return g;
}

/** Connecting rod: rod-local origin at the crank-pin centre, +y toward the wrist pin, z = crank axis. */
export function buildRod(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'connecting-rod';
  const R = L.rod;
  const len = R.length;
  const rIn = R.pinRadius + R.babbitt;
  const a = R.boltOffset + R.boltRadius * 1.6;
  const Rse = R.smallEndOuterRadius;
  const ySe = len - Rse * 0.75;
  const yS0 = R.bossHeight - 0.002;
  const circle = (cx: number, cy: number, r: number, n = 40): P2[] => Array.from({ length: n }, (_, k) => [cx + r * Math.cos((k * 2 * Math.PI) / n), cy + r * Math.sin((k * 2 * Math.PI) / n)] as P2);
  const forged: THREE.BufferGeometry[] = [];
  // big end (rod half + cap) with the bearing bore
  const be = hullOfCircles([
    { x: -a + 0.005, y: R.bossHeight - 0.005, r: 0.005 }, { x: a - 0.005, y: R.bossHeight - 0.005, r: 0.005 },
    { x: -a + 0.006, y: -R.capDepth + 0.006, r: 0.006 }, { x: a - 0.006, y: -R.capDepth + 0.006, r: 0.006 },
  ], 12);
  forged.push(extrudeZ(polyShape(be, [circle(0, 0, rIn, 48)]), -R.bigEndWidth / 2, R.bigEndWidth / 2, 0.001));
  // I-beam shank: web through the middle, flanges along both edges
  const shank: P2[] = [[-R.shankWidthBig / 2, yS0], [R.shankWidthBig / 2, yS0], [R.shankWidthSmall / 2, ySe], [-R.shankWidthSmall / 2, ySe]];
  forged.push(extrudeZ(polyShape(shank), -R.webThickness / 2, R.webThickness / 2, 0.0005));
  for (const s of [-1, 1]) {
    const f: P2[] = [[s * R.shankWidthBig / 2, yS0], [s * R.shankWidthSmall / 2, ySe], [s * (R.shankWidthSmall / 2 - R.flangeThickness), ySe], [s * (R.shankWidthBig / 2 - R.flangeThickness), yS0]];
    forged.push(extrudeZ(polyShape(f), -R.smallEndWidth / 2, R.smallEndWidth / 2, 0.0008));
  }
  // clamped small end with the pinch-bolt lug on the camshaft side [DB23]
  forged.push(extrudeZ(polyShape(circle(0, len, Rse, 40), [circle(0, len, R.wristPinRadius + 0.0001, 40)]), -R.smallEndWidth / 2, R.smallEndWidth / 2, 0.0008));
  const lug: P2[] = [[-Rse - 0.007, len - 0.009], [-Rse * 0.6, len - 0.012], [-Rse * 0.6, len + 0.012], [-Rse - 0.007, len + 0.009]];
  forged.push(extrudeZ(polyShape(lug), -R.smallEndWidth / 2 + 0.002, R.smallEndWidth / 2 - 0.002, 0.0006));
  g.add(named(new THREE.Mesh(merge(forged), M.m.steelForged), 'rod-body'));
  for (const geo of forged) geo.dispose();
  // babbitt lining
  g.add(new THREE.Mesh(extrudeZ(polyShape(circle(0, 0, rIn, 48), [circle(0, 0, R.pinRadius + 0.00005, 48)]), -R.bigEndWidth / 2 + 0.0005, R.bigEndWidth / 2 - 0.0005), M.m.bearing));
  // cap bolts + castellated nuts, pinch bolt
  const bolts = new MeshBuilder();
  for (const s of [-1, 1]) {
    const x = s * R.boltOffset, br = R.boltRadius;
    lathe(bolts, [{ r: 0, y: -R.capDepth - 0.009 }, { r: br, y: -R.capDepth - 0.009 }, { r: br, y: R.bossHeight }, { r: 0, y: R.bossHeight }], { segments: 12, center: [x, 0, 0] });
    lathe(bolts, [{ r: 0, y: R.bossHeight }, { r: br * 1.7, y: R.bossHeight }, { r: br * 1.7, y: R.bossHeight + 0.005 }, { r: 0, y: R.bossHeight + 0.005 }], { segments: 6, center: [x, 0, 0] });
    lathe(bolts, [{ r: 0, y: -R.capDepth - 0.008 }, { r: br * 1.7, y: -R.capDepth - 0.008 }, { r: br * 1.7, y: -R.capDepth }, { r: 0, y: -R.capDepth }], { segments: 6, center: [x, 0, 0] });
  }
  cylinderBetween(bolts, [-Rse - 0.0035, len - 0.016, 0], [-Rse - 0.0035, len + 0.016, 0], 0.0028, 10);
  g.add(new THREE.Mesh(bolts.build(), M.m.steelDark));
  return g;
}

/** Piston: piston-local origin on the wrist-pin axis, +y toward the crown. */
export function buildPiston(L: ModelTLayout, M: ModelTMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'piston';
  const P = L.piston;
  const R = P.radius;
  const CH = L.compressionHeight;
  const yb = P.skirtBottom;
  const rIn = R - P.grooves[0].depth - 0.0025;
  const outer: ProfilePoint[] = [{ r: R - 0.0006, y: yb }, { r: R, y: yb + 0.0006 }];
  for (const gr of P.grooves.slice().sort((p, q) => q.top - p.top)) {
    const y1 = CH - gr.top, y0 = y1 - gr.height;
    outer.push({ r: R, y: y0 }, { r: R - gr.depth, y: y0 }, { r: R - gr.depth, y: y1 }, { r: R, y: y1 });
  }
  outer.push({ r: R, y: CH - P.topChamfer }, { r: R - P.topChamfer, y: CH });
  const loop: ProfilePoint[] = [
    { r: rIn, y: yb },
    ...outer,
    { r: 0, y: CH },
    { r: 0, y: CH - P.crownThickness },
    { r: rIn - 0.004, y: CH - P.crownThickness, smooth: true },
    { r: rIn, y: CH - P.crownThickness - 0.004 },
  ];
  const b = new MeshBuilder();
  lathe(b, loop, { segments: 72, closed: true });
  g.add(named(new THREE.Mesh(b.build(), M.m.pistonIron), 'piston-body'));
  // pin bosses
  const zIn = L.rod.smallEndWidth / 2 + 0.001;
  const zOut = Math.sqrt(Math.max(rIn * rIn - P.bossRadius ** 2, 0)) + 0.001;
  const bosses = [
    latheZ([{ r: P.pinRadius, y: zIn }, { r: P.bossRadius, y: zIn }, { r: P.bossRadius, y: zOut }, { r: P.pinRadius, y: zOut }], { segments: 36, closed: true }),
    latheZ([{ r: P.pinRadius, y: -zOut }, { r: P.bossRadius, y: -zOut }, { r: P.bossRadius, y: -zIn }, { r: P.pinRadius, y: -zIn }], { segments: 36, closed: true }),
  ];
  g.add(new THREE.Mesh(merge(bosses), M.m.pistonIron));
  for (const geo of bosses) geo.dispose();
  // rings
  const rb = new MeshBuilder();
  for (const gr of P.grooves) {
    const y1 = CH - gr.top - 0.0001, y0 = CH - gr.top - gr.height + 0.0001;
    const ri = R - gr.depth + 0.0005, ro = R - 0.00003;
    lathe(rb, [{ r: ri, y: y0 }, { r: ro, y: y0 }, { r: ro, y: y1 }, { r: ri, y: y1 }], { segments: 72, closed: true });
  }
  g.add(named(new THREE.Mesh(rb.build(), M.m.steelDark), 'piston-rings'));
  // hollow wrist pin
  const hl = P.pinLength / 2;
  g.add(named(new THREE.Mesh(latheZ([
    { r: P.pinBore, y: -hl }, { r: P.pinRadius - 0.0005, y: -hl }, { r: P.pinRadius, y: -hl + 0.0005 },
    { r: P.pinRadius, y: hl - 0.0005 }, { r: P.pinRadius - 0.0005, y: hl }, { r: P.pinBore, y: hl },
  ], { segments: 36, closed: true }), M.m.steelPolished), 'wrist-pin'));
  return g;
}

/** Spur-gear disc (closed) along z centred at (x, y). */
export function spurGear(teeth: number, module: number, z0: number, z1: number, x = 0, y = 0, phase = 0, boreR = 0): THREE.BufferGeometry {
  const holes: P2[][] = boreR > 0 ? [hullOfCircles([{ x: 0, y: 0, r: boreR }], 32)] : [];
  const geo = extrudeZ(polyShape(gearOutline(teeth, (teeth * module) / 2, module, phase), holes), z0, z1, 0.0005, 2);
  geo.translate(x, y, 0);
  return geo;
}

