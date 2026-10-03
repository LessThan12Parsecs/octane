/**
 * Connecting rod (rod-local frame: origin at the crank-pin centre, +y toward
 * the wrist pin, z = crank axis) and piston (piston-local frame: origin on
 * the wrist-pin axis, +y toward the crown).
 */
import * as THREE from 'three';
import type { EngineLayout } from './layout';
import type { EngineMaterials } from './materials';
import { MeshBuilder, circlePath, extrude, lathe, type ProfilePoint, crease } from './geometry';

function centred(g: THREE.BufferGeometry, depth: number): THREE.BufferGeometry {
  g.translate(0, 0, -depth / 2);
  return g;
}

function annulusShape(cx: number, cy: number, rOut: number, rIn: number): THREE.Shape {
  const s = new THREE.Shape();
  const n = 64;
  for (let k = 0; k <= n; k++) {
    const a = (k * 2 * Math.PI) / n;
    const x = cx + rOut * Math.cos(a), y = cy + rOut * Math.sin(a);
    if (k === 0) s.moveTo(x, y);
    else s.lineTo(x, y);
  }
  s.holes.push(circlePath(cx, cy, rIn, n));
  return s;
}

export function buildRod(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'connecting-rod';
  const len = L.rodLength;
  const rIn = L.crankPinRadius + L.rodBearingThickness;
  const a = L.rodBoltOffset + L.rodBoltRadius * 1.7; // boss half width
  const hb = 0.022; // boss height above the split line
  const wb = L.rodShankWidthBig, ws = L.rodShankWidthSmall;
  const Rse = L.rodSmallEndOuterRadius;
  const yShank0 = hb + 0.018;
  const ySe = len - Math.sqrt(Math.max(Rse * Rse - (ws / 2) ** 2, 0));
  const rSeIn = L.wristPinRadius + L.rodBushingThickness;

  // Outline of rod body (upper big end + shank + small-end eye), CCW.
  const body = (halfShank: (y: number) => number) => {
    const s = new THREE.Shape();
    s.moveTo(a, 0);
    s.lineTo(a, hb);
    s.quadraticCurveTo(a, yShank0 - 0.004, halfShank(yShank0), yShank0 + 0.01);
    s.lineTo(halfShank(ySe), ySe);
    const a0 = Math.atan2(ySe - len, halfShank(ySe));
    s.absarc(0, len, Rse, a0, Math.PI - a0, false);
    s.lineTo(-halfShank(yShank0), yShank0 + 0.01);
    s.quadraticCurveTo(-a, yShank0 - 0.004, -a, hb);
    s.lineTo(-a, 0);
    s.lineTo(-rIn, 0);
    s.absarc(0, 0, rIn, Math.PI, 0, true);
    s.lineTo(a, 0);
    return s;
  };
  const half = (y: number) => {
    const t = Math.min(1, Math.max(0, (y - yShank0) / (ySe - yShank0)));
    return wb / 2 + (ws / 2 - wb / 2) * t;
  };
  // 1) thin web through the whole outline
  const webShape = body(half);
  webShape.holes.push(circlePath(0, len, rSeIn));
  g.add(new THREE.Mesh(centred(extrude(webShape, { depth: L.rodWebThickness, bevel: 0.0008, bevelSegments: 1 }), L.rodWebThickness), M.m.steelForged));

  // 2) flanges: strips along the shank edges, full small-end width
  const tf = L.rodFlangeThickness;
  for (const sgn of [-1, 1]) {
    const s = new THREE.Shape();
    const y0 = yShank0 + 0.004, y1 = ySe + 0.002;
    s.moveTo(sgn * half(y0), y0);
    s.lineTo(sgn * half(y1), y1);
    s.lineTo(sgn * (half(y1) - tf), y1);
    s.lineTo(sgn * (half(y0) - tf), y0);
    s.closePath();
    const f = centred(extrude(s, { depth: L.rodSmallEndWidth, bevel: 0.0012, bevelSegments: 2 }), L.rodSmallEndWidth);
    g.add(new THREE.Mesh(f, M.m.steelForged));
  }

  // 3) big-end upper boss (full width)
  const be = new THREE.Shape();
  be.moveTo(a, 0);
  be.lineTo(a, hb);
  be.quadraticCurveTo(a, yShank0 - 0.004, half(yShank0), yShank0 + 0.01);
  be.lineTo(-half(yShank0), yShank0 + 0.01);
  be.quadraticCurveTo(-a, yShank0 - 0.004, -a, hb);
  be.lineTo(-a, 0);
  be.lineTo(-rIn, 0);
  be.absarc(0, 0, rIn, Math.PI, 0, true);
  be.lineTo(a, 0);
  g.add(new THREE.Mesh(centred(extrude(be, { depth: L.rodBigEndWidth, bevel: 0.0014, bevelSegments: 2 }), L.rodBigEndWidth), M.m.steelForged));

  // 4) cap
  const cd = L.rodCapDepth;
  const cap2 = new THREE.Shape();
  cap2.moveTo(a, 0);
  cap2.lineTo(rIn, 0);
  cap2.absarc(0, 0, rIn, 0, -Math.PI, true);
  cap2.lineTo(-a, 0);
  cap2.lineTo(-a, -cd * 0.55);
  cap2.quadraticCurveTo(-a, -cd, -a * 0.45, -cd);
  cap2.lineTo(a * 0.45, -cd);
  cap2.quadraticCurveTo(a, -cd, a, -cd * 0.55);
  cap2.lineTo(a, 0);
  g.add(new THREE.Mesh(centred(extrude(cap2, { depth: L.rodBigEndWidth, bevel: 0.0014, bevelSegments: 2 }), L.rodBigEndWidth), M.m.steelForged));

  // 5) bearing shells + small-end bushing
  const shell = centred(extrude(annulusShape(0, 0, rIn, L.crankPinRadius + 0.00005), { depth: L.rodBigEndWidth - 0.001, curveSegments: 64 }), L.rodBigEndWidth - 0.001);
  g.add(new THREE.Mesh(shell, M.m.bearing));
  const bush = centred(extrude(annulusShape(0, len, rSeIn, L.wristPinRadius + 0.00005), { depth: L.rodSmallEndWidth - 0.0005 }), L.rodSmallEndWidth - 0.0005);
  g.add(new THREE.Mesh(bush, M.m.bearing));
  // small-end eye (full width)
  const eye = centred(extrude(annulusShape(0, len, Rse, rSeIn), { depth: L.rodSmallEndWidth, bevel: 0.0012, bevelSegments: 2 }), L.rodSmallEndWidth);
  g.add(new THREE.Mesh(eye, M.m.steelForged));

  // 6) cap bolts + nuts (along y)
  const bolts = new MeshBuilder();
  for (const sgn of [-1, 1]) {
    const x = sgn * L.rodBoltOffset;
    const br = L.rodBoltRadius;
    // shank + head above the boss, nut below the cap
    lathe(bolts, [
      { r: 0, y: -cd - 0.012 }, { r: br, y: -cd - 0.012 }, { r: br, y: hb }, { r: 0, y: hb },
    ], { segments: 16, center: [x, 0, 0] });
    lathe(bolts, [
      { r: 0, y: hb }, { r: br * 1.75, y: hb }, { r: br * 1.75, y: hb + 0.006 }, { r: br * 1.4, y: hb + 0.0075 }, { r: 0, y: hb + 0.0075 },
    ], { segments: 6, center: [x, 0, 0] });
    lathe(bolts, [
      { r: 0, y: -cd - 0.0085 }, { r: br * 1.4, y: -cd - 0.0085 }, { r: br * 1.75, y: -cd - 0.007 }, { r: br * 1.75, y: -cd }, { r: 0, y: -cd },
    ], { segments: 6, center: [x, 0, 0] });
  }
  g.add(new THREE.Mesh(crease(bolts.build(), 0.5), M.m.steelDark));
  return g;
}

export function buildPiston(L: EngineLayout, M: EngineMaterials): THREE.Group {
  const g = new THREE.Group();
  g.name = 'piston';
  const R = L.pistonRadius;
  const CH = L.compressionHeight;
  const yb = L.pistonSkirtBottom;
  const tw = L.pistonWallThickness;
  const tc = L.pistonCrownThickness;
  const grooves = L.ringGrooves;
  const beltBottom = CH - (grooves[grooves.length - 1].top + grooves[grooves.length - 1].height) - 0.002;
  const rBeltIn = R - grooves[0].depth - 0.0035;
  const rSkirtIn = R - tw;

  // Outer profile (bottom → top), with grooves, crown chamfer.
  const outer: ProfilePoint[] = [{ r: R - 0.0006, y: yb }, { r: R, y: yb + 0.0006 }];
  for (const gr of grooves.slice().reverse()) {
    const y1 = CH - gr.top; // groove upper edge
    const y0 = y1 - gr.height;
    outer.push({ r: R, y: y0 });
    outer.push({ r: R - gr.depth, y: y0 });
    outer.push({ r: R - gr.depth, y: y1 });
    outer.push({ r: R, y: y1 });
  }
  outer.push({ r: R, y: CH - 0.0008 });
  outer.push({ r: R - 0.0008, y: CH });
  const loop: ProfilePoint[] = [
    { r: rSkirtIn, y: yb },
    ...outer,
    { r: 0, y: CH },
    { r: 0, y: CH - tc },
    { r: rBeltIn - 0.004, y: CH - tc, smooth: true },
    { r: rBeltIn, y: CH - tc - 0.004 },
    { r: rBeltIn, y: beltBottom },
    { r: rSkirtIn, y: beltBottom - 0.004 },
  ];
  const b = new MeshBuilder();
  lathe(b, loop, { segments: 96, closed: true });
  g.add(Object.assign(new THREE.Mesh(b.build(), M.m.pistonIron), { name: 'piston-body' }));

  // Pin bosses (inside the skirt), along z.
  const zIn = L.rodSmallEndWidth / 2 + 0.0012;
  const zOut = Math.sqrt(Math.max(rSkirtIn * rSkirtIn - L.pinBossRadius ** 2, 0)) + 0.001;
  for (const s of [-1, 1]) {
    const bb = new MeshBuilder();
    lathe(bb, [
      { r: L.wristPinRadius, y: zIn }, { r: L.pinBossRadius, y: zIn }, { r: L.pinBossRadius, y: zOut }, { r: L.wristPinRadius, y: zOut },
    ], { segments: 48, closed: true });
    const geo = bb.build();
    geo.rotateX(Math.PI / 2);
    if (s < 0) geo.rotateY(Math.PI);
    g.add(new THREE.Mesh(geo, M.m.pistonIron));
  }

  // Rings in the grooves.
  const rb = new MeshBuilder();
  for (const gr of grooves) {
    const y1 = CH - gr.top - 0.0001, y0 = CH - gr.top - gr.height + 0.0001;
    const ri = R - gr.depth + 0.0004, ro = R - 0.00003;
    if (!gr.oil) {
      lathe(rb, [{ r: ri, y: y0 }, { r: ro, y: y0 }, { r: ro, y: y1 }, { r: ri, y: y1 }], { segments: 96, closed: true });
    } else {
      // oil control ring: two rails with a slotted expander between
      const h = y1 - y0;
      lathe(rb, [{ r: ri, y: y0 }, { r: ro, y: y0 }, { r: ro, y: y0 + h * 0.3 }, { r: ri, y: y0 + h * 0.3 }], { segments: 96, closed: true });
      lathe(rb, [{ r: ri, y: y1 - h * 0.3 }, { r: ro, y: y1 - h * 0.3 }, { r: ro, y: y1 }, { r: ri, y: y1 }], { segments: 96, closed: true });
      lathe(rb, [{ r: ri, y: y0 + h * 0.3 }, { r: ro - 0.0012, y: y0 + h * 0.3 }, { r: ro - 0.0012, y: y1 - h * 0.3 }, { r: ri, y: y1 - h * 0.3 }], { segments: 96, closed: true });
    }
  }
  g.add(Object.assign(new THREE.Mesh(rb.build(), M.m.steelDark), { name: 'piston-rings' }));

  // Wrist pin (hollow).
  const pb = new MeshBuilder();
  const hl = L.wristPinLength / 2;
  lathe(pb, [
    { r: L.wristPinBore, y: -hl }, { r: L.wristPinRadius - 0.0005, y: -hl }, { r: L.wristPinRadius, y: -hl + 0.0005 },
    { r: L.wristPinRadius, y: hl - 0.0005 }, { r: L.wristPinRadius - 0.0005, y: hl }, { r: L.wristPinBore, y: hl },
  ], { segments: 48, closed: true });
  const pin = pb.build();
  pin.rotateX(Math.PI / 2);
  g.add(Object.assign(new THREE.Mesh(pin, M.m.steelPolished), { name: 'wrist-pin' }));
  return g;
}
