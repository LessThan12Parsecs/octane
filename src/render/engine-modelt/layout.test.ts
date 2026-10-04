import { describe, expect, it } from 'vitest';
import { MODEL_T } from '../../physics/engines/model-t';
import { CFR_F1 } from '../../physics/engines/cfr';
import { computeModelTLayout, pistonDisplacementAt } from './layout';
import { distanceToPolygon, insideConvex, insidePolygon, polygonArea } from './outline';
import { plugXZ } from './parts-head';

const IN = 0.0254;

/** Mechanical clearances and period dimensions of the procedural Model T layout. */
describe('Model T layout', () => {
  const L = computeModelTLayout(MODEL_T);
  const lh = MODEL_T.geometry.lHead!;

  it('uses the spec geometry and the period vertical stack', () => {
    expect(L.bore).toBe(MODEL_T.geometry.bore);
    expect(L.crankRadius).toBe(MODEL_T.geometry.stroke / 2);
    expect(L.rodLength).toBe(MODEL_T.geometry.conRodLength);
    expect(L.axisZ).toEqual(MODEL_T.layout.axisZ);
    // block 10-5/8 in high from the crank centreline to the deck [DB23]; crown 5/16 in proud at TDC [FM19 A22]
    expect(L.deckY).toBeCloseTo(10.625 * IN, 9);
    expect(L.crownTdcY - L.deckY).toBeCloseTo((5 / 16) * IN, 12);
    // the cylinder frame origin is the head-cavity roof: crown at −h_TDC, deck at lHead.deckY
    expect(L.crownTdcY - L.cylOriginY).toBeCloseTo(lh.deckY + lh.crownAboveDeckAtTDC, 12);
    expect(L.deckY - L.cylOriginY).toBeCloseTo(lh.deckY, 12);
    expect(L.pocketRoofY - L.cylOriginY).toBeCloseTo(lh.pocket.roofY, 12);
    // high head: 2-11/16 in from the gasket face to the bolt-head seats
    expect(L.boltBossTopY - L.headBottomY).toBeCloseTo((2 + 11 / 16) * IN, 12);
    expect(L.head.bolts).toHaveLength(15);
  });

  it('rejects specs without an L-head chamber', () => {
    expect(() => computeModelTLayout(CFR_F1)).toThrow();
  });

  it('crank: throws 1&4 / 2&3 in one plane, pins clear of the journals, mains inside the crankcase', () => {
    const a = L.throwAngle;
    expect(Math.cos(a[0] - a[3])).toBeCloseTo(1, 12);
    expect(Math.cos(a[1] - a[2])).toBeCloseTo(1, 12);
    expect(Math.cos(a[0] - a[1])).toBeCloseTo(-1, 12);
    expect(L.crankRadius - L.crank.journalRadius - L.crank.pinRadius).toBeGreaterThan(0.015);
    // front and centre mains where the spec puts them; the rear one slid back to clear cylinder 4's web
    expect(L.crank.mains[0].z).toBeCloseTo(MODEL_T.layout.mainBearingZ[0], 12);
    expect(L.crank.mains[1].z).toBeCloseTo(MODEL_T.layout.mainBearingZ[1], 12);
    expect(Math.abs(L.crank.mains[2].z - MODEL_T.layout.mainBearingZ[2])).toBeLessThan(0.015);
    const webRear = L.axisZ[3] - L.crank.pinLength / 2 - L.crank.webThickness;
    expect(L.crank.mains[2].z + L.crank.mains[2].length / 2).toBeLessThan(webRear);
    // overall crank length 25-5/32 in [DB23]
    expect(L.crank.frontEndZ - L.crank.flangeZ0).toBeCloseTo((25 + 5 / 32) * IN, 9);
    // crank sweep inside the crankcase cavity
    const reach = L.crankRadius + L.rod.boltOffset + L.rod.boltRadius * 1.7;
    expect(L.block.cc.x1).toBeGreaterThan(reach);
    expect(-L.block.cc.x0).toBeGreaterThan(reach);
  });

  it('rods and pistons clear the camshaft and lobes over the whole revolution', () => {
    const C = L.cam;
    const lobeR = C.maxRadius;
    let minRodCam = Infinity, minSkirtLobe = Infinity;
    for (let th = 0; th < 360; th += 1) {
      const phi = (th * Math.PI) / 180;
      const px = L.crankRadius * Math.sin(phi), py = L.crankRadius * Math.cos(phi);
      const wy = L.wristPinTdcY - pistonDisplacementAt(L, th);
      const ux = -px, uy = wy - py; // rod axis
      const ul = Math.hypot(ux, uy);
      const ex = ux / ul, ey = uy / ul, nx = -ey, ny = ex;
      // big-end envelope: bolt heads / nuts and the cap corners (rod frame) — the rod is 1.5 in wide and sits
      // between the lobes only through the camshaft (journal) radius, so check against the shaft
      for (const [lx, ly] of [
        [L.rod.boltOffset + L.rod.boltRadius * 1.7, L.rod.bossHeight + 0.005], [-(L.rod.boltOffset + L.rod.boltRadius * 1.7), L.rod.bossHeight + 0.005],
        [L.rod.boltOffset + L.rod.boltRadius * 1.7, -L.rod.capDepth - 0.008], [-(L.rod.boltOffset + L.rod.boltRadius * 1.7), -L.rod.capDepth - 0.008],
      ] as [number, number][]) {
        const x = px + nx * lx + ex * ly, y = py + ny * lx + ey * ly;
        minRodCam = Math.min(minRodCam, Math.hypot(x - C.x, y - C.y) - C.journalRadius);
      }
      // skirt corner vs lobe cylinder (lobes sit beside the bore, ±(24 ± 11) mm from the cylinder axis)
      const skirtY = wy + L.piston.skirtBottom;
      const dx = Math.max(0, Math.abs(C.x) - L.piston.radius);
      const dy = Math.max(0, skirtY - C.y);
      minSkirtLobe = Math.min(minSkirtLobe, Math.hypot(dx, dy) - lobeR);
    }
    expect(minRodCam).toBeGreaterThan(0.0005);
    expect(minSkirtLobe).toBeGreaterThan(0.002);
  });

  it('valvetrain stack: seat → stem end → lash → push rod → cam, tappets straight over their lobes', () => {
    // stem end 4.974 in below the seat line [T424B], lash at the cam heel, push rod 2-11/32 in, base radius
    expect(L.seatLineY + L.valve.stemEndY - L.tappet.topRestY).toBeCloseTo(MODEL_T.intakeValve.lash!, 12);
    expect(L.tappet.topRestY - L.tappet.faceRestY).toBeCloseTo((2 + 11 / 32) * IN, 12);
    expect(L.tappet.faceRestY - L.cam.y).toBeCloseTo(0.406 * IN, 12);
    // cam-to-crank centre distance from the stack-up ≈ 3.9–4.0 in (valvetrain research: Dyke's gear diameters)
    const d = Math.hypot(L.cam.x, L.cam.y);
    expect(d).toBeGreaterThan(0.099);
    expect(d).toBeLessThan(0.102);
    expect(L.gears.crankPitchR + L.gears.camPitchR).toBeCloseTo(d, 12);
    expect(L.gears.camTeeth / L.gears.crankTeeth).toBe(2);
    // the 48-tooth gear lands near Dyke's 5-1/3 in (42 T straight-cut) diameter
    expect(2 * L.gears.camPitchR / IN).toBeGreaterThan(5.0);
    expect(2 * L.gears.camPitchR / IN).toBeLessThan(5.5);
    // spring installed at 2-1/8 in between the cup and the valve-chamber roof [FS par. 258]
    expect(L.spring.topY - L.spring.bottomRestY).toBeCloseTo((2 + 1 / 8) * IN, 12);
    expect(L.spring.topY).toBeCloseTo(L.block.vc.y1, 12);
    // valve chamber contains the cups / tappet tops over the full lift; ports above it
    const maxLift = MODEL_T.intakeValve.maxLift;
    expect(L.spring.bottomRestY + maxLift).toBeLessThan(L.block.vc.y1 - 0.02);
    expect(L.tappet.topRestY).toBeGreaterThan(L.block.vc.y0);
    expect(L.block.portFloorY).toBeGreaterThan(L.block.vc.y1 + 0.01);
    // tappet feet stay inside the crankcase cavity below the guide wall
    expect(L.tappet.faceRestY + L.tappet.footThickness + (L.cam.maxRadius - L.cam.baseRadius)).toBeLessThan(L.block.cc.yTop);
  });

  it('eight valves E-I-I-E-E-I-I-E on the right side (−x), mirrored where layout.mirrorZ, lobes phased by firing order', () => {
    const order = L.valves.slice().sort((a, b) => b.z - a.z).map((v) => v.kind[0].toUpperCase()).join('');
    expect(order).toBe('EIIEEIIE');
    for (const v of L.valves) {
      expect(v.x).toBeLessThan(0);
      // nose at the tappet when the cylinder's local angle is at the lobe centre
      const camAngleAtCentre = ((v.firingOffsetDeg + v.cam.centreDeg) * Math.PI) / 360;
      expect(Math.cos(camAngleAtCentre + v.lobePhase)).toBeCloseTo(1, 12);
    }
    // the cam bearings do not overlap the lobes; the shaft reaches the timer
    for (const b of L.cam.bearings) for (const v of L.valves) {
      expect(Math.abs(b.z - v.z)).toBeGreaterThan(b.length / 2 + L.cam.lobeWidth / 2);
    }
    expect(L.cam.z1).toBeGreaterThan(L.timer.z0);
    // six manifold ports: four exhaust + two siamesed intakes [Good22 p. 21]
    expect(L.block.ports.filter((p) => p.kind === 'exhaust')).toHaveLength(4);
    expect(L.block.ports.filter((p) => p.kind === 'intake')).toHaveLength(2);
  });

  it('valves sit beside the bores, inside the pocket, and clear the pocket roof at full lift', () => {
    const V = L.valve;
    for (const v of L.valves) {
      const zc = L.axisZ[v.cyl];
      const r = Math.hypot(v.x, v.z - zc);
      expect(r - V.seatOuterRadius).toBeGreaterThan(L.boreRadius + 0.003); // iron between seat and bore
      const pocket = L.pockets[v.cyl];
      expect(insideConvex(pocket, v.x, v.z)).toBe(true);
      expect(distanceToPolygon(pocket, v.x, v.z)).toBeGreaterThan(V.headRadius + 0.002);
      const top = L.seatLineY + V.topY + MODEL_T.intakeValve.maxLift;
      expect(L.pocketRoofY - top).toBeGreaterThan(0.005);
    }
  });

  it('chamber outlines: gasket opening = bore + pocket (18.3 in² measured), roof regions consistent', () => {
    const ch = L.chambers[0];
    const Abore = polygonArea(ch.circle);
    expect(polygonArea(ch.union)).toBeCloseTo(Abore + polygonArea(ch.pocketRoof), 12);
    expect(Math.abs(polygonArea(ch.pocketRoof) - 46.016e-4) / 46.016e-4).toBeLessThan(0.004);
    expect(Math.abs(polygonArea(ch.union) - 18.3 * IN * IN) / (18.3 * IN * IN)).toBeLessThan(0.01);
    // head-only chamber volume ≈ 284 cm³ (spec build-up; Tulsa 294 cm³ before the 9/1918 −1/16 in change)
    const Vhead = polygonArea(ch.union) * (L.pocketRoofY - L.headBottomY) + Abore * (L.cylOriginY - L.pocketRoofY);
    expect(Vhead * 1e6).toBeGreaterThan(275);
    expect(Vhead * 1e6).toBeLessThan(294);
  });

  it('spark plugs: gap in the pocket under the roof, between the valves, plug clear of the bore column', () => {
    const P = L.plug;
    for (let c = 0; c < L.nCyl; c++) {
      const [x, z] = plugXZ(L, c);
      expect(insidePolygon(L.chambers[c].pocketRoof, x, z)).toBe(true);
      expect(Math.hypot(x, z - L.axisZ[c]) - L.head.plugPassageRadius).toBeGreaterThan(L.boreRadius);
      for (const v of L.valves.filter((w) => w.cyl === c)) expect(Math.hypot(v.x - x, v.z - z)).toBeGreaterThan(L.valve.headRadius + 0.002);
    }
    expect(P.noseDistance).toBeCloseTo(L.pocketRoofY - (L.cylOriginY + P.gapCenter[1]), 12);
    expect(P.seatDistance).toBeGreaterThan(P.noseDistance + (L.headTopY - L.pocketRoofY));
  });

  it('flywheel magneto: magnets inside the bell, coil faces 1/32 in from the magnet poles, ring gear 120 T', () => {
    const F = L.flywheel;
    expect(F.ringGear.teeth).toBe(120);
    expect(F.magnets.count).toBe(16);
    expect(F.magnets.z1 + L.coilRing.gap).toBeCloseTo(L.coilRing.coilZ0, 12);
    expect(L.coilRing.gap).toBeCloseTo(IN / 32, 12);
    expect(F.radius).toBeLessThan(L.bell.inner[0].r - 0.005);
    expect(L.coilRing.rOut).toBeLessThan(L.bell.inner[0].r);
    // coil cores face the magnet poles
    expect(L.coilRing.coilR).toBeGreaterThan(F.magnets.rApex);
    expect(L.coilRing.coilR).toBeLessThan(F.magnets.rPole);
    // the flywheel stays behind the block and in front of the bell's taper
    expect(F.hubZ1).toBeLessThan(L.block.zB);
    expect(F.discZ0).toBeGreaterThan(L.bell.inner[1].z);
    expect(F.tripleGears.r + ((Math.max(...F.tripleGears.teeth) + 2) * F.tripleGears.module) / 2).toBeLessThan(L.bell.inner[2].r);
  });

  it('front: gears inside the timing cover, pulley/belt/fan/timer clear of one another', () => {
    const G = L.gears;
    expect(G.z0).toBeGreaterThan(L.block.zF);
    expect(L.crank.pulleyZ0).toBeGreaterThan(G.z1);
    expect(L.timer.z0).toBeGreaterThan(L.belt.z + L.belt.width / 2);
    expect(L.fan.z).toBeGreaterThan(L.timer.z1 - 0.01);
    // fan disc clear of the timer body
    expect(Math.hypot(L.timer.x - L.fan.x, L.timer.y - L.fan.y) - L.timer.radius).toBeGreaterThan(L.fan.radius);
    // generator body outside the block face, under the valve door
    expect(L.generator.x + L.generator.radius).toBeLessThan(L.block.xR);
    expect(L.generator.y + L.generator.radius).toBeLessThan(L.block.vc.y0 - 0.008);
  });

  it('cutaway frames: quarter section through cylinder 1’s front valve, bell frame behind the bell', () => {
    expect(L.cutSide).toBe(-1);
    const v1 = L.valves.filter((v) => v.cyl === 0).map((v) => v.z);
    expect(L.cutZ).toBe(Math.max(...v1));
    expect(L.bellCutZ).toBeLessThan(L.bell.inner[L.bell.inner.length - 1].z - L.bell.wall);
    for (const k of [0, 1, 2] as const) expect(L.bounds.max[k]).toBeGreaterThan(L.bounds.min[k]);
  });
});
