import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/gasex_lhead_flow.json';
import type { LHeadChamberSpec, ValveSpec } from '../core/engine-spec';
import { CFR_F1 } from '../engines/cfr';
import { MODEL_T } from '../engines/model-t';
import {
  createValveFlowModel,
  isSideValve,
  lHeadPocketTransfer,
  pocketBoreBoundaryLength,
  POCKET_TRANSFER_CD,
  seriesEffectiveArea,
  SLOT_CONTRACTION_COEFFICIENT,
  ValveFlowModel,
} from './valve-flow';

interface Case {
  name: string;
  headDiameter: number;
  seatInnerDiameter: number;
  seatAngle: number;
  stemDiameter: number;
  portDiameter: number;
  bore: number;
  pocket: LHeadChamberSpec['pocket'] extends infer P ? Omit<P, 'roofY'> : never;
  roofClearance: number;
  pocketHeight: number;
  boundaryLength: number;
  rows: { lift: number; area: number; stage: number; roofArea: number }[];
  transfer: { crownAboveDeck: number; transferArea: number; series: { cdaValve: number; cdT: number; cdaSeries: number }[] }[];
}
const cases = (fixture as unknown as { cases: Case[] }).cases;

const rel = (a: number, b: number) => Math.abs(a - b) / Math.abs(b);

const valveOf = (c: Case): ValveSpec => ({
  count: 1,
  headDiameter: c.headDiameter,
  seatInnerDiameter: c.seatInnerDiameter,
  seatAngle: c.seatAngle,
  stemDiameter: c.stemDiameter,
  maxLift: 0.9 * c.roofClearance,
  openDeg: -350,
  closeDeg: -130,
  timingLiftThreshold: 0,
  position: [-0.06, 0],
  shroudArcDeg: 0,
  shroudDirection: 0,
  liftDirection: 1,
});

describe('side valve in an L-head pocket: roof stage and pocket → bore transfer (oracle gasex_lhead_flow.py)', () => {
  for (const c of cases) {
    it(`${c.name}: minimum area with the roof stage, boundary arc, transfer area and series combination`, () => {
      const lHead: LHeadChamberSpec = { deckY: -0.03, crownAboveDeckAtTDC: 0.005, pocket: { ...c.pocket, roofY: -0.03 + c.pocketHeight } };
      const w = pocketBoreBoundaryLength(lHead, c.bore);
      expect(rel(w, c.boundaryLength)).toBeLessThan(1e-11);
      const m = new ValveFlowModel(valveOf(c), 'intake', {
        portDiameter: c.portDiameter,
        roofClearance: c.roofClearance,
        transfer: { width: w, height: c.pocketHeight },
      });
      for (const r of c.rows) {
        expect(rel(m.flowArea(r.lift), r.area)).toBeLessThan(1e-12);
        expect(m.stage(r.lift)).toBe(r.stage);
        expect(rel(m.roofArea(r.lift), r.roofArea)).toBeLessThan(1e-12);
      }
      for (const t of c.transfer) {
        const at = m.transferArea(t.crownAboveDeck);
        if (t.transferArea === 0) expect(at).toBe(0);
        else expect(rel(at, t.transferArea)).toBeLessThan(1e-11);
        for (const s of t.series) {
          const v = seriesEffectiveArea(s.cdaValve, s.cdT * at);
          if (s.cdaSeries === 0) expect(v).toBe(0);
          else expect(rel(v, s.cdaSeries)).toBeLessThan(1e-11);
        }
      }
    });
  }

  it('effective area: roof masks the curtain C_D like a shroud; transfer in series, direction-dependent C_D', () => {
    const c = cases[0];
    const spec = valveOf(c);
    const tr = { width: c.boundaryLength, height: c.pocketHeight };
    const bare = new ValveFlowModel(spec, 'intake', { portDiameter: c.portDiameter });
    const roof = new ValveFlowModel(spec, 'intake', { portDiameter: c.portDiameter, roofClearance: c.roofClearance });
    const full = new ValveFlowModel(spec, 'intake', { portDiameter: c.portDiameter, roofClearance: c.roofClearance, transfer: tr });
    const ex = new ValveFlowModel(spec, 'exhaust', { portDiameter: c.portDiameter, roofClearance: c.roofClearance, transfer: tr });
    for (const L of [5e-4, 2e-3, 5e-3, 8e-3, 0.0125]) {
      for (const rev of [false, true]) {
        const a0 = bare.effectiveArea(L, rev);
        const a1 = roof.effectiveArea(L, rev);
        expect(rel(a1, a0 * (roof.flowArea(L) / bare.flowArea(L)))).toBeLessThan(1e-13);
        // without crown height the transfer stage is off
        expect(full.effectiveArea(L, rev)).toBe(a1);
        const crown = 0.002;
        const cdT = rev ? POCKET_TRANSFER_CD.outflow : POCKET_TRANSFER_CD.inflow; // intake forward = into the bore
        expect(rel(full.effectiveArea(L, rev, crown), seriesEffectiveArea(a1, cdT * full.transferArea(crown)))).toBeLessThan(1e-14);
        expect(full.effectiveArea(L, rev, crown)).toBeLessThan(a1);
        const cdTx = rev ? POCKET_TRANSFER_CD.inflow : POCKET_TRANSFER_CD.outflow; // exhaust forward = out of the bore
        expect(rel(ex.effectiveArea(L, rev, crown), seriesEffectiveArea(ex.effectiveArea(L, rev), cdTx * ex.transferArea(crown)))).toBeLessThan(1e-14);
      }
    }
    expect(POCKET_TRANSFER_CD.outflow).toBeCloseTo(0.6110154703516573, 15); // π/(π+2)
    expect(SLOT_CONTRACTION_COEFFICIENT).toBe(POCKET_TRANSFER_CD.outflow);
    // the roof closes the flow continuously as the head reaches it
    expect(roof.flowArea(c.roofClearance)).toBe(0);
    expect(roof.effectiveArea(c.roofClearance * (1 - 1e-9), false)).toBeLessThan(1e-12);
    // crown above the pocket roof closes the transfer
    expect(full.effectiveArea(3e-3, false, c.pocketHeight + 1e-4)).toBe(0);
    expect(seriesEffectiveArea(2e-4, Infinity)).toBe(2e-4);
    expect(seriesEffectiveArea(0, 1e-4)).toBe(0);
  });

  it('createValveFlowModel: bit-identical for the CFR (overhead valves); side valves get roof + transfer', () => {
    for (const kind of ['intake', 'exhaust'] as const) {
      const v = kind === 'intake' ? CFR_F1.intakeValve : CFR_F1.exhaustValve;
      const port = kind === 'intake' ? CFR_F1.manifolds.intakePortDiameter : CFR_F1.manifolds.exhaustPortDiameter;
      const a = createValveFlowModel(v, kind, CFR_F1);
      const b = new ValveFlowModel(v, kind, { portDiameter: Math.max(port, v.stemDiameter * 1.01) });
      expect(isSideValve(v, CFR_F1)).toBe(false);
      expect(a.transfer).toBe(null);
      expect(a.roofClearance).toBe(Infinity);
      for (let L = 0; L < 8e-3; L += 1.37e-5) {
        for (const rev of [false, true]) {
          expect(a.effectiveArea(L, rev)).toBe(b.effectiveArea(L, rev));
          expect(a.effectiveArea(L, rev, 0.001)).toBe(b.effectiveArea(L, rev)); // no transfer for OHV
        }
        expect(a.flowArea(L)).toBe(b.flowArea(L));
        expect(a.stage(L)).toBe(b.stage(L));
      }
    }
    // Model T: side valves in the L-head pocket
    const lh = MODEL_T.geometry.lHead!;
    const m = createValveFlowModel(MODEL_T.intakeValve, 'intake', MODEL_T);
    expect(isSideValve(MODEL_T.intakeValve, MODEL_T)).toBe(true);
    expect(m.roofClearance).toBeCloseTo(lh.pocket.roofY - lh.deckY, 15);
    expect(m.roofClearance * 1e3).toBeCloseTo(12.9, 1);
    const tr = lHeadPocketTransfer(MODEL_T);
    expect(m.transfer!.width).toBe(tr.width);
    expect(tr.width * 1e3).toBeGreaterThan(70); // bore arc inside the pocket (≈ 80 mm)
    expect(tr.width * 1e3).toBeLessThan(90);
    // the roof does not bind at the Model T lift (seat-limited, stage 2, just below the 1-1/8 in port)
    expect(m.stage(MODEL_T.intakeValve.maxLift)).toBe(2);
    expect(m.roofArea(MODEL_T.intakeValve.maxLift)).toBeGreaterThan(m.portArea);
    // transfer near TDC (crown 5/16 in above the deck) vs mid-stroke
    const L = MODEL_T.intakeValve.maxLift;
    const mid = m.effectiveArea(L, false, lh.crownAboveDeckAtTDC - 0.05);
    const tdc = m.effectiveArea(L, false, lh.crownAboveDeckAtTDC);
    const none = m.effectiveArea(L, false);
    expect(mid).toBeLessThan(none);
    expect(mid / none).toBeGreaterThan(0.85);
    expect(tdc).toBeLessThan(mid);
    expect(m.transferArea(lh.crownAboveDeckAtTDC) * 1e6).toBeCloseTo(tr.width * (tr.height - lh.crownAboveDeckAtTDC) * 1e6, 9);
  });

  it('pocket boundary arc: whole circle, no intersection, symmetric plans', () => {
    const R = 0.04;
    const big: LHeadChamberSpec = { deckY: -0.03, crownAboveDeckAtTDC: 0, pocket: { xMin: -1, xMax: 1, zMin: -1, zMax: 1, cornerRadius: 0.01, roofY: -0.02 } };
    expect(rel(pocketBoreBoundaryLength(big, 2 * R), 2 * Math.PI * R)).toBeLessThan(1e-14);
    const far: LHeadChamberSpec = { ...big, pocket: { xMin: -0.2, xMax: -0.1, zMin: -0.02, zMax: 0.02, cornerRadius: 0.005, roofY: -0.02 } };
    expect(pocketBoreBoundaryLength(far, 2 * R)).toBe(0);
    // half plane x ≤ 0 (huge rectangle): half the circumference
    const half: LHeadChamberSpec = { ...big, pocket: { xMin: -1, xMax: 0, zMin: -1, zMax: 1, cornerRadius: 0, roofY: -0.02 } };
    expect(rel(pocketBoreBoundaryLength(half, 2 * R), Math.PI * R)).toBeLessThan(1e-12);
  });
});
