import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import { computeLayout } from './layout';
import { sliderCrank } from './kinematics';

/** Mechanical clearances of the procedural layout for the CFR spec. */
describe('CFR layout clearances', () => {
  const L = computeLayout(CFR_F1);
  const [crMin, crMax] = L.crRange;

  it('uses the spec geometry', () => {
    expect(L.bore).toBe(CFR_F1.geometry.bore);
    expect(L.crankRadius).toBe(CFR_F1.geometry.stroke / 2);
    expect(L.rodLength).toBe(CFR_F1.geometry.conRodLength);
    expect(L.boreRadius).toBe(CFR_F1.geometry.bore / 2);
    expect(L.pistonRadius).toBeLessThan(L.boreRadius);
    expect(L.boreRadius - L.pistonRadius).toBeLessThan(3e-4);
  });

  it('crown never reaches the head (positive clearance over the CR range)', () => {
    expect(L.headYMin - L.crownTdcY).toBeCloseTo(L.clearanceAtTdc(crMax), 12);
    expect(L.headYMin).toBeGreaterThan(L.crownTdcY);
  });

  it('piston skirt clears the counterweights at BDC', () => {
    const bdc = L.wristPinTdcY - sliderCrank(L, 180).pistonDisplacement;
    expect(bdc + L.pistonSkirtBottom).toBeGreaterThan(L.counterweightRadius + 0.005);
  });

  it('cylinder barrel bottom stays above the crank sweep and inside the deck over the full CR travel', () => {
    const barrelBottomLow = L.headYMin - L.cylinderLength;
    const barrelBottomHigh = L.headYMax - L.cylinderLength;
    const bigEndTop = L.crankRadius + L.rodBigEndOuterRadius;
    expect(barrelBottomLow).toBeGreaterThan(Math.max(L.counterweightRadius, bigEndTop) + 0.005);
    expect(barrelBottomHigh).toBeLessThan(L.crankcase.yDeckBottom + 1e-9);
    // jacket shoulder above the worm wheel at the lowest position
    expect(L.headYMin - L.jacketDepth).toBeGreaterThan(L.wormWheel.yBottom + L.wormWheel.thickness);
  });

  it('piston rings stay inside the bore at BDC for every CR', () => {
    for (const cr of [crMin, crMax]) {
      const h = L.clearanceAtTdc(cr) + L.stroke;
      const lowestRing = L.ringGrooves[L.ringGrooves.length - 1];
      expect(h + lowestRing.top + lowestRing.height).toBeLessThan(L.cylinderLength);
    }
  });

  it('crank sweep fits inside the crankcase', () => {
    expect(L.crankcase.xHalfInner).toBeGreaterThan(L.crankRadius + L.rodBoltOffset);
    expect(-L.crankcase.yFloorInner).toBeGreaterThan(L.crankRadius + L.rodCapDepth);
    expect(L.crankcase.zFrontInner).toBeGreaterThan(L.webZ + L.webThickness / 2);
  });

  it('cam lobes are behind the crank throw and clear the journal; gears mesh at 2:1', () => {
    for (const v of L.valves) expect(v.lobeZ + L.lobeWidth / 2).toBeLessThan(-(L.webZ + L.webThickness / 2));
    expect(L.camY - L.camBaseRadius * 1.4 - 0.006).toBeGreaterThan(L.mainJournalRadius);
    expect(L.crankGearRadius + L.camGearRadius).toBeCloseTo(Math.hypot(L.camX, L.camY), 12);
    expect(L.camGearRadius / L.crankGearRadius).toBeCloseTo(2, 12);
    expect(L.camGearTeeth).toBe(2 * L.crankGearTeeth);
  });

  it('pushrods clear the water jacket and worm wheel', () => {
    for (const v of L.valves) {
      const clear = Math.abs(v.lobeZ) - L.pushrodPassageRadius;
      expect(clear).toBeGreaterThan(L.jacketOuterRadius);
      expect(clear).toBeGreaterThan(L.wormWheel.outerRadius);
      expect(v.lobeZ).toBeGreaterThan(L.headZMin);
    }
  });

  it('rocker arms satisfy the rocker ratio and reach the valve axis', () => {
    for (const v of L.valves) {
      expect(v.armValve / v.armPushrod).toBeCloseTo(1.4, 12);
      const padZ = v.rockerPivot[2] + v.armValve * v.rockerDir[1];
      const padX = v.rockerPivot[0] + v.armValve * v.rockerDir[0];
      expect(padX).toBeCloseTo(v.x, 12);
      expect(padZ).toBeCloseTo(v.z, 12);
      expect(v.rockerPivot[1]).toBeCloseTo(v.tipY, 12);
      // spring never goes solid at max lift
      const coils = v.springTurns + 1;
      expect(v.springInstalledLength - v.spec.maxLift).toBeGreaterThan(coils * 2 * v.springWireRadius);
    }
  });

  it('ports stay inside the head', () => {
    for (const v of L.valves) {
      expect(v.portY + v.portExitRadius).toBeLessThan(L.headThickness);
      expect(v.portY - v.portExitRadius).toBeGreaterThan(0);
      expect(v.insertHeight).toBeLessThan(v.portY);
    }
  });

  it('spark plug: electrodes inside the chamber, seat outside the casting', () => {
    const p = L.plug;
    expect(p.noseDistance).toBeGreaterThan(CFR_F1.sparkPlug.gap / 2);
    expect(p.seatDistance).toBeGreaterThan(p.noseDistance);
    const dot = p.side[0] * p.axis[0] + p.side[1] * p.axis[1] + p.side[2] * p.axis[2];
    expect(Math.abs(dot)).toBeLessThan(1e-12);
    // jacket cavity kept clear of the plug boss
    expect(L.jacketCavity.yTop).toBeLessThan(p.gapCenter[1] - p.bossRadius);
  });
});
