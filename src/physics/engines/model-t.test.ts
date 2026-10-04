import { describe, expect, it } from 'vitest';
import lheadFx from '../../../test/fixtures/combustion_geometry_lhead_mc.json';
import { cylinderAngleDeg, firingOffsetsDeg, hasVariableCompressionRatio } from '../core/engine-spec';
import { chamberFixedVolume, createChamber, lHeadPlanMetrics, sliderCrankGeometry } from '../combustion/chamber';
import { SliderCrank } from '../mechanics/kinematics';
import { MODEL_T, MODEL_T_LAYOUT, MODEL_T_POCKET_PLAN_AREA } from './model-t';

const IN = 0.0254;
const g = MODEL_T.geometry;
const lh = g.lHead!;
const R = g.bore / 2;
const Ap = Math.PI * R * R;

describe('Model T spec: basic geometry', () => {
  it('3.750 × 4.000 in, 176.7 in³ (2.896 L) — Ford 1923 Data Book; McCalley', () => {
    expect(g.bore).toBeCloseTo(3.75 * IN, 12);
    expect(g.stroke).toBeCloseTo(4.0 * IN, 12);
    const Vd = Ap * g.stroke * MODEL_T.cylinders;
    expect(Vd * 1e3).toBeCloseTo(2.896, 3);
    expect(Vd / IN ** 3).toBeCloseTo(176.7, 1);
    expect(MODEL_T.cylinders).toBe(4);
    expect(g.conRodLength).toBeCloseTo(7 * IN, 12);
  });

  it('fixed compression ratio 3.98 (no variable-CR range)', () => {
    expect(g.compressionRatio).toBe(3.98);
    expect(g.compressionRatioRange).toEqual([3.98, 3.98]);
    expect(hasVariableCompressionRatio(MODEL_T)).toBe(false);
  });

  it("side-valve 'l-head' chamber: crown 5/16 in above the deck at TDC, ≈ 1 in to the head", () => {
    expect(g.chamber).toBe('l-head');
    expect(lh.crownAboveDeckAtTDC).toBeCloseTo((5 / 16) * IN, 12);
    expect(-lh.deckY - lh.crownAboveDeckAtTDC).toBeCloseTo(1.0 * IN, 12);
    expect(lh.deckY).toBeLessThan(lh.pocket.roofY);
    expect(lh.pocket.roofY).toBeLessThan(0);
    // pocket roof ≈ 12.9 mm above the deck (model-t.ts chamber note)
    expect((lh.pocket.roofY - lh.deckY) * 1e3).toBeGreaterThan(12);
    expect((lh.pocket.roofY - lh.deckY) * 1e3).toBeLessThan(14);
    expect(g.quenchCreviceVolume!).toBeLessThanOrEqual(g.creviceVolume);
  });
});

describe('Model T spec: chamber volume closure', () => {
  const m = lHeadPlanMetrics(g.bore, lh);

  it('MODEL_T_POCKET_PLAN_AREA is the exact pocket plan area (Green, chamber.ts) and the independent quadrature', () => {
    expect(MODEL_T_POCKET_PLAN_AREA).toBe(m.pocketArea);
    expect(Math.abs(lheadFx.modelT.pocketArea / m.pocketArea - 1)).toBeLessThan(1e-12);
  });

  it('the oracle fixture rebuilt the same chamber from the sourced inputs (deck, roof)', () => {
    expect(Math.abs(lheadFx.modelT.deckY - lh.deckY)).toBeLessThan(1e-15);
    expect(Math.abs(lheadFx.modelT.roofY - lh.pocket.roofY)).toBeLessThan(1e-13);
  });

  it('CR, crevice, pocket and h_TDC close the clearance volume (< 1e-9)', () => {
    const Vd = Ap * g.stroke;
    const Vc = Vd / (g.compressionRatio - 1);
    const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
    const closure = Ap * hTdc + chamberFixedVolume(g) + g.creviceVolume;
    expect(Math.abs(closure / Vc - 1)).toBeLessThan(1e-9);
    expect(Math.abs((Vd + closure) / closure / g.compressionRatio - 1)).toBeLessThan(1e-9);
  });

  it('SliderCrank with the chamber fixed volume: clearanceHeightTDC = −deckY − crown rise (bore-column depth)', () => {
    const kin = new SliderCrank(sliderCrankGeometry(g), g.compressionRatio);
    const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
    expect(Math.abs(kin.clearanceHeightTDC / hTdc - 1)).toBeLessThan(1e-9);
    expect(kin.fixedChamberVolume).toBe(m.fixedVolume);
    // volume = A_p h + V_crevice + V_fixed at any angle
    for (const th of [0, 0.7, 2, Math.PI]) {
      const V = kin.volume(th);
      expect(Math.abs((Ap * kin.clearanceHeight(th) + g.creviceVolume + m.fixedVolume) / V - 1)).toBeLessThan(1e-12);
    }
    // the trap: the bare spec geometry treats the pocket as part of a bore-sized disc (h_TDC ≈ 34 mm)
    const bare = new SliderCrank(g, g.compressionRatio);
    expect(bare.clearanceHeightTDC).toBeGreaterThan(0.033);
  });

  it('chamber: chamberVolume(h_TDC) = V_c − V_crevice; plan area ≈ 18.3 in² gasket opening [Tulsa]', () => {
    const ch = createChamber(MODEL_T);
    const Vc = (Ap * g.stroke) / (g.compressionRatio - 1);
    const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
    expect(Math.abs(ch.chamberVolume(hTdc) / (Vc - g.creviceVolume) - 1)).toBeLessThan(1e-9);
    expect(ch.fixedVolume).toBe(m.fixedVolume);
    expect(ch.planformArea / IN ** 2).toBeGreaterThan(18.3 * 0.98);
    expect(ch.planformArea / IN ** 2).toBeLessThan(18.3 * 1.01);
    expect(m.planformArea).toBeCloseTo(ch.planformArea, 15);
  });
});

describe('Model T spec: valves and plug in the chamber', () => {
  const valves = [MODEL_T.intakeValve, MODEL_T.exhaustValve];
  const pk = lh.pocket;
  const inRR = (x: number, z: number): boolean => {
    const c = pk.cornerRadius;
    const xc = Math.min(Math.max(x, pk.xMin + c), pk.xMax - c);
    const zc = Math.min(Math.max(z, pk.zMin + c), pk.zMax - c);
    return x >= pk.xMin && x <= pk.xMax && z >= pk.zMin && z <= pk.zMax && Math.hypot(x - xc, z - zc) <= c;
  };

  it('side valves seat in the deck and lift upward into the pocket', () => {
    for (const v of valves) {
      expect(v.seatY).toBe(lh.deckY);
      expect(v.liftDirection).toBe(1);
      expect(v.count).toBe(1);
      expect(v.headDiameter).toBeGreaterThan(v.seatInnerDiameter);
    }
  });

  it('valve heads lie inside the pocket outline and outside the bore', () => {
    for (const v of valves) {
      const [x, z] = v.position;
      const rv = v.headDiameter / 2;
      expect(Math.hypot(x, z) - rv).toBeGreaterThan(R);
      for (let k = 0; k < 360; k++) {
        const a = (k * Math.PI) / 180;
        expect(inRR(x + rv * Math.cos(a), z + rv * Math.sin(a))).toBe(true);
      }
    }
  });

  it('valve heads clear of each other', () => {
    const [a, b] = valves;
    const gap = Math.hypot(a.position[0] - b.position[0], a.position[1] - b.position[1]) - a.headDiameter / 2 - b.headDiameter / 2;
    expect(gap).toBeGreaterThan(0.005);
  });

  it('plug gap inside the chamber: over the pocket (outside the bore planform), between floor and roof, clear of the valve faces', () => {
    const [x, y, z] = MODEL_T.sparkPlug.gapCenter;
    expect(Math.hypot(x, z)).toBeGreaterThan(R); // d > R: the chamber model's outside-spark branch
    expect(inRR(x, z)).toBe(true);
    expect(y).toBeGreaterThan(lh.deckY);
    expect(y).toBeLessThan(lh.pocket.roofY);
    expect(y).toBeLessThanOrEqual(0);
    for (const v of valves) expect(Math.hypot(x - v.position[0], z - v.position[1])).toBeGreaterThan(v.headDiameter / 2);
    // the gap is no closer to the roof than half the electrode gap
    expect(lh.pocket.roofY - y).toBeGreaterThan(MODEL_T.sparkPlug.gap);
  });

  it('the chamber model accepts the spec (valves and plug validated) and covers the spark-to-far-corner distance', () => {
    const ch = createChamber(MODEL_T);
    const hTdc = -lh.deckY - lh.crownAboveDeckAtTDC;
    expect(ch.inscribedRadius(hTdc)).toBeCloseTo(lh.pocket.roofY - MODEL_T.sparkPlug.gapCenter[1], 12);
    // farthest point: the far side of the bore at the crown
    const [x, y] = MODEL_T.sparkPlug.gapCenter;
    expect(ch.maxRadius(hTdc)).toBeGreaterThanOrEqual(Math.hypot(R + Math.abs(x), Math.max(-y, Math.abs(-y - hTdc))) - 1e-12);
  });
});

describe('Model T spec: cylinder layout', () => {
  const L = MODEL_T_LAYOUT;

  it('arrays have one entry per cylinder (three mains)', () => {
    expect(MODEL_T.layout).toBe(L);
    expect(L.firingOrder).toHaveLength(4);
    expect(L.firingOffsetDeg).toHaveLength(4);
    expect(L.axisZ).toHaveLength(4);
    expect(L.mirrorZ).toHaveLength(4);
    expect(L.mainBearingZ).toHaveLength(3);
    expect([...L.firingOrder].sort()).toEqual([1, 2, 3, 4]);
  });

  it('firing offsets follow the firing order 1-2-4-3 at 180° intervals', () => {
    expect(L.firingOrder).toEqual([1, 2, 4, 3]);
    L.firingOrder.forEach((cyl, k) => expect(L.firingOffsetDeg[cyl - 1]).toBe(180 * k));
    expect(firingOffsetsDeg(MODEL_T)).toEqual(L.firingOffsetDeg);
    for (const off of L.firingOffsetDeg) expect(off >= 0 && off < 720).toBe(true);
    // each cylinder is at its own firing TDC when the engine angle equals its offset
    L.firingOffsetDeg.forEach((off, i) => expect(cylinderAngleDeg(MODEL_T, i, off)).toBe(0));
  });

  it('crank throws: 1 & 4 and 2 & 3 in phase (offsets equal mod 360°), the pairs 180° apart', () => {
    const o = L.firingOffsetDeg;
    expect(o[0] % 360).toBe(o[3] % 360);
    expect(o[1] % 360).toBe(o[2] % 360);
    expect(Math.abs((o[1] % 360) - (o[0] % 360))).toBe(180);
  });

  it('cylinder 1 at the front (+z); bore spacing 4-1/8, 5-1/4, 4-1/8 in, symmetric about z = 0; mains outboard', () => {
    const z = L.axisZ;
    expect(z[0]).toBeGreaterThan(z[1]);
    expect(z[1]).toBeGreaterThan(z[2]);
    expect(z[2]).toBeGreaterThan(z[3]);
    expect(z[0] - z[1]).toBeCloseTo(4.125 * IN, 12);
    expect(z[1] - z[2]).toBeCloseTo(5.25 * IN, 12);
    expect(z[2] - z[3]).toBeCloseTo(4.125 * IN, 12);
    expect(z[0] + z[3]).toBeCloseTo(0, 12);
    for (let i = 0; i < 3; i++) expect(z[i] - z[i + 1]).toBeGreaterThan(g.bore); // walls between bores
    expect(L.mainBearingZ[0]).toBeGreaterThan(z[0]);
    expect(L.mainBearingZ[2]).toBeLessThan(z[3]);
    expect(L.mainBearingZ[1]).toBeLessThan(z[1]);
    expect(L.mainBearingZ[1]).toBeGreaterThan(z[2]);
  });

  it('mirrored valve order E-I-I-E-E-I-I-E: cylinders 2 and 4 mirrored', () => {
    expect(L.mirrorZ).toEqual([false, true, false, true]);
    // cylinder 1 (unmirrored): exhaust in front (+z) of the intake
    expect(MODEL_T.exhaustValve.position[1]).toBeGreaterThan(MODEL_T.intakeValve.position[1]);
  });
});
