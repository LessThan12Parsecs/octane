/**
 * Detachable "high head" (1917–27): a cast-iron slab with one L-head chamber per cylinder (bore column
 * under a roof at the cylinder-frame origin + the valve pocket roof, engine-spec.ts LHeadChamberSpec), a
 * water jacket above the chambers with cast tubes for the four plugs, 15 bolt bosses, plug bosses and the
 * water-outlet gooseneck at the front. World coordinates.
 *
 * The head shell is a watertight 2-manifold built from the shared chamber outlines (outline.ts): the
 * gasket-face hole, the walls, the pocket roof and the bore roof meet edge to edge.
 * Sizes not taken from layout.ts (fillets, chamfers, bosses, caps, fastener heads) are render-only
 * cosmetic choices — UNVERIFIED.
 */
import { MeshBuilder, type V3, lathe, planarFace, rect, ringXZ } from '../engine/geometry';
import type { ShellBuilders } from './parts-block';
import type { FrameSet } from './frame-set';
import type { ModelTLayout } from './layout';
import type { P2 } from './outline';
import { pathFrames, ringSolidY, tubeSolid } from './geom';

export const PLUG_SEG = 32;

/** World (x, z) of cylinder c's plug axis. */
export function plugXZ(L: ModelTLayout, c: number): P2 {
  const g = L.plug.gapCenter;
  return [g[0], L.axisZ[c] + (L.mirror[c] ? -g[2] : g[2])];
}

/** Vertical walls along a CCW (x, z) outline from y0 to y1 with normals INTO the outline (cavity). */
function cavityWalls(b: MeshBuilder, outline: readonly P2[], y0: number, y1: number): void {
  for (let i = 0; i < outline.length; i++) {
    const p = outline[i], q = outline[(i + 1) % outline.length];
    const dx = q[0] - p[0], dz = q[1] - p[1];
    const l = Math.hypot(dx, dz) || 1;
    const n: V3 = [-dz / l, 0, dx / l];
    b.quad([p[0], y0, p[1]], [q[0], y0, q[1]], [q[0], y1, q[1]], [p[0], y1, p[1]], n, n, n, n);
  }
}

export function buildHeadShell(L: ModelTLayout): ShellBuilders {
  const H = L.head;
  const paint = new MeshBuilder(), machined = new MeshBuilder(), raw = new MeshBuilder();
  const y0 = L.headBottomY, y1 = L.headTopY;
  const outline = rect(H.xR, H.zB, H.xL, H.zF);
  const plugs = L.axisZ.map((_, c) => plugXZ(L, c));
  const rp = H.plugPassageRadius, rt = H.plugTubeRadius;

  planarFace(machined, outline, L.chambers.map((c) => c.union), [0, y0, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(paint, outline, plugs.map(([x, z]) => ringXZ(x, z, rp, PLUG_SEG)), [0, y1, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(paint, rect(y0, H.zB, y1, H.zF), [], [H.xR, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(paint, rect(y0, H.zB, y1, H.zF), [], [H.xL, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(paint, rect(H.xR, y0, H.xL, y1), [], [0, 0, H.zF], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(paint, rect(H.xR, y0, H.xL, y1), [], [0, 0, H.zB], [1, 0, 0], [0, 1, 0], [0, 0, -1]);

  // chambers: walls of the gasket opening up to the pocket roof, pocket roof (with the plug hole), bore column up to its roof
  L.chambers.forEach((ch, c) => {
    const [px, pz] = plugs[c];
    cavityWalls(raw, ch.union, y0, L.pocketRoofY);
    planarFace(raw, ch.pocketRoof, [ringXZ(px, pz, rp, PLUG_SEG)], [0, L.pocketRoofY, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
    cavityWalls(raw, ch.circle, L.pocketRoofY, L.cylOriginY);
    planarFace(raw, ch.circle, [], [0, L.cylOriginY, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
    lathe(machined, [{ r: rp, y: y1 }, { r: rp, y: L.pocketRoofY }], { segments: PLUG_SEG, center: [px, 0, pz] });
  });

  // water jacket over the chambers, crossed by the cast plug tubes
  const j = H.jacket;
  const tubes = plugs.map(([x, z]) => ringXZ(x, z, rt, PLUG_SEG));
  planarFace(raw, rect(j.x0, j.z0, j.x1, j.z1), tubes, [0, j.y0, 0], [1, 0, 0], [0, 0, 1], [0, 1, 0]);
  planarFace(raw, rect(j.x0, j.z0, j.x1, j.z1), tubes, [0, j.y1, 0], [1, 0, 0], [0, 0, 1], [0, -1, 0]);
  planarFace(raw, rect(j.y0, j.z0, j.y1, j.z1), [], [j.x0, 0, 0], [0, 1, 0], [0, 0, 1], [1, 0, 0]);
  planarFace(raw, rect(j.y0, j.z0, j.y1, j.z1), [], [j.x1, 0, 0], [0, 1, 0], [0, 0, 1], [-1, 0, 0]);
  planarFace(raw, rect(j.x0, j.y0, j.x1, j.y1), [], [0, 0, j.z0], [1, 0, 0], [0, 1, 0], [0, 0, 1]);
  planarFace(raw, rect(j.x0, j.y0, j.x1, j.y1), [], [0, 0, j.z1], [1, 0, 0], [0, 1, 0], [0, 0, -1]);
  for (const [x, z] of plugs) lathe(raw, [{ r: rt, y: j.y0 }, { r: rt, y: j.y1 }], { segments: PLUG_SEG, center: [x, 0, z] });
  return { paint, machined, raw };
}

export function buildHead(L: ModelTLayout, set: FrameSet): void {
  const H = L.head;
  const sh = buildHeadShell(L);
  set.add('paint', sh.paint.build());
  set.add('ironMachined', sh.machined.build());
  set.add('castRaw', sh.raw.build());

  // bolt bosses (2-11/16 in to the bolt-head underside) and hex heads
  const bosses = new MeshBuilder();
  const heads = new MeshBuilder();
  const hexR = H.boltHeadAcrossFlats / 2 / Math.cos(Math.PI / 6);
  for (const [x, z] of H.bolts) {
    bosses.append(ringSolidY(L.headTopY - 0.001, L.boltBossTopY, 0, H.boltBossRadius, x, z, 28));
    lathe(heads, [
      { r: 0, y: L.boltBossTopY }, { r: hexR * 1.02, y: L.boltBossTopY }, { r: hexR * 1.02, y: L.boltBossTopY + 0.0012 },
      { r: hexR, y: L.boltBossTopY + 0.0014 }, { r: hexR, y: L.boltBossTopY + H.boltHeadHeight - 0.001 },
      { r: hexR * 0.85, y: L.boltBossTopY + H.boltHeadHeight }, { r: 0, y: L.boltBossTopY + H.boltHeadHeight },
    ], { segments: 6, center: [x, 0, z], phiStart: Math.PI / 6 });
  }
  set.add('paint', bosses.build());
  set.add('steelDark', heads.build());

  // plug bosses around the plug holes
  const pb = new MeshBuilder();
  for (let c = 0; c < L.nCyl; c++) {
    const [x, z] = plugXZ(L, c);
    pb.append(ringSolidY(L.headTopY - 0.001, L.headTopY + 0.25 * 0.0254, H.plugPassageRadius, L.plug.bossRadius, x, z, PLUG_SEG));
  }
  set.add('paint', pb.build());

  // water-outlet gooseneck at the front top (to the radiator top hose)
  {
    const g = H.gooseneck;
    const fl = ringSolidY(L.headTopY - 0.001, L.headTopY + 0.006, g.boreRadius, g.flangeRadius, g.x, g.z, 40);
    set.add('paint', fl);
    const y0 = L.headTopY + 0.006;
    const pts: V3[] = [[g.x, y0, g.z]];
    const R = g.rise;
    const n = 14;
    pts.push([g.x, y0 + 0.004, g.z]);
    for (let i = 1; i <= n; i++) {
      const a = (i / n) * (Math.PI / 2);
      pts.push([g.x, y0 + 0.004 + R * Math.sin(a), g.z + R * (1 - Math.cos(a))]);
    }
    pts.push([g.x, y0 + 0.004 + R, g.z + g.reach]);
    const fr = pathFrames(pts, g.boreRadius, [0, 0, 1]);
    const o = new MeshBuilder(), i = new MeshBuilder();
    tubeSolid(o, i, fr, 0.0045, 32);
    set.add('paint', o.build());
    set.add('castRaw', i.build());
    // hose bead at the spigot
    const tip = pts[pts.length - 1];
    const bead = new MeshBuilder();
    lathe(bead, [{ r: g.boreRadius + 0.0045, y: -0.004 }, { r: g.boreRadius + 0.008, y: -0.001 }, { r: g.boreRadius + 0.008, y: 0.001 }, { r: g.boreRadius + 0.0045, y: 0.004 }], { segments: 32 });
    const bg = bead.build();
    bg.rotateX(Math.PI / 2);
    bg.translate(tip[0], tip[1], tip[2] - 0.012);
    set.add('paint', bg, { section: false });
  }
}
