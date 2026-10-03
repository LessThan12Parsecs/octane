import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../physics/engines/cfr';
import {
  crankAngleInRevolution,
  cycleDialAngleDeg,
  cycleDialPoint,
  cylinderVolumes,
  describeCrankAngle,
  isValveOpen,
  nearestTdc,
  strokeOf,
  strokeProgress,
  valveOpenInterval,
  wrapCycleDeg,
} from './engine-cycle';
import { dialArcPath } from './hud';

describe('crank-angle bookkeeping', () => {
  it('wraps into [-360, 360)', () => {
    expect(wrapCycleDeg(360)).toBe(-360);
    expect(wrapCycleDeg(-360)).toBe(-360);
    expect(wrapCycleDeg(375)).toBe(-345);
    expect(wrapCycleDeg(-361)).toBe(359);
    expect(wrapCycleDeg(1080)).toBe(-360);
    expect(wrapCycleDeg(12.5)).toBe(12.5);
  });

  it('names strokes with the firing-TDC convention', () => {
    expect(strokeOf(-360)).toBe('intake');
    expect(strokeOf(-180.01)).toBe('intake');
    expect(strokeOf(-180)).toBe('compression');
    expect(strokeOf(-0.01)).toBe('compression');
    expect(strokeOf(0)).toBe('power');
    expect(strokeOf(179.9)).toBe('power');
    expect(strokeOf(180)).toBe('exhaust');
    expect(strokeOf(359.9)).toBe('exhaust');
    expect(strokeOf(360)).toBe('intake');
  });

  it('stroke progress runs 0→1 over each stroke', () => {
    expect(strokeProgress(-180)).toBe(0);
    expect(strokeProgress(-90)).toBeCloseTo(0.5, 12);
    expect(strokeProgress(270)).toBeCloseTo(0.5, 12);
  });

  it('mechanical crank angle and nearest TDC', () => {
    expect(crankAngleInRevolution(-13)).toBeCloseTo(347, 12);
    expect(crankAngleInRevolution(-350)).toBeCloseTo(10, 12);
    expect(nearestTdc(-10)).toBe('firing');
    expect(nearestTdc(-300)).toBe('gas-exchange');
    expect(nearestTdc(330)).toBe('gas-exchange');
    expect(nearestTdc(180)).toBeNull();
  });

  it('describes angles relative to TDC/BDC', () => {
    expect(describeCrankAngle(-13)).toBe('13.0° BTDC');
    expect(describeCrankAngle(5.5)).toBe('5.5° ATDC');
    expect(describeCrankAngle(0)).toBe('TDC');
    expect(describeCrankAngle(-350)).toBe('10.0° ATDC');
    expect(describeCrankAngle(200)).toBe('20.0° ABDC');
    expect(describeCrankAngle(-146)).toBe('34.0° ABDC');
    expect(describeCrankAngle(140)).toBe('40.0° BBDC');
    expect(describeCrankAngle(180)).toBe('BDC');
  });

  it('cycle dial puts firing TDC at the top and strokes in quadrants', () => {
    expect(cycleDialAngleDeg(0)).toBe(0);
    expect(cycleDialAngleDeg(180)).toBe(90);
    expect(cycleDialAngleDeg(-180)).toBe(270);
    expect(cycleDialAngleDeg(-360)).toBe(180);
    const [x, y] = cycleDialPoint(0, 50, 50, 10);
    expect(x).toBeCloseTo(50, 12);
    expect(y).toBeCloseTo(40, 12);
    const [x2, y2] = cycleDialPoint(180, 50, 50, 10); // 3 o'clock
    expect(x2).toBeCloseTo(60, 12);
    expect(y2).toBeCloseTo(50, 12);
  });

  it('dial arcs use the large-arc flag beyond a half turn', () => {
    expect(dialArcPath(0, 180, 10)).toContain(' 0 0 1 ');
    expect(dialArcPath(-300, 200, 10)).toContain(' 0 1 1 ');
  });
});

describe('valve timing', () => {
  it('handles events straddling the cycle boundary', () => {
    const ev = { openDeg: 140, closeDeg: -345 }; // EVO 40° BBDC, EVC 15° ATDC (gas exchange)
    expect(valveOpenInterval(ev)).toEqual([140, 375]);
    expect(isValveOpen(ev, 200)).toBe(true);
    expect(isValveOpen(ev, -350)).toBe(true); // = 370
    expect(isValveOpen(ev, -340)).toBe(false);
    expect(isValveOpen(ev, 0)).toBe(false);
  });

  it('reads the CFR intake event from the spec', () => {
    const [o, c] = valveOpenInterval(CFR_F1.intakeValve);
    expect(c).toBeGreaterThan(o);
    const mid = (o + c) / 2;
    expect(isValveOpen(CFR_F1.intakeValve, mid)).toBe(true);
    expect(isValveOpen(CFR_F1.intakeValve, 20)).toBe(false); // power stroke
  });
});

describe('cylinder volumes', () => {
  it('Vmax / Vc equals the compression ratio', () => {
    for (const cr of [4, 7, 12.5, 18]) {
      const v = cylinderVolumes(CFR_F1, cr);
      expect(v.vmax / v.vc).toBeCloseTo(cr, 10);
    }
    const g = CFR_F1.geometry;
    expect(cylinderVolumes(CFR_F1, 7).vd).toBeCloseTo((Math.PI / 4) * g.bore * g.bore * g.stroke, 15);
  });
});
