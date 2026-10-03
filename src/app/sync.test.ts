import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../physics/core/engine-spec';
import { CFR_F1, clearanceHeightAtTDC, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { computeLayout } from '../render/engine/layout';
import { MockKinematics } from '../worker/mock-simulator';
import {
  clampOperatingPoint,
  compressionRatioFromSnapshot,
  FRAMING_REFERENCE_ASPECT,
  framingDistanceScale,
  FrameClock,
  parseUrlOptions,
  pistonTravel,
  RecentWindow,
} from './sync';

describe('compressionRatioFromSnapshot', () => {
  it('inverts the physics clearance height at TDC (crevice included) across the CR range', () => {
    for (const cr of [4, 5.5, 6.5, 7.55, 10, 18]) {
      const h = clearanceHeightAtTDC(cr, CFR_F1);
      for (const disp of [0, 0.01, 0.1143]) {
        expect(compressionRatioFromSnapshot(CFR_F1, { clearanceHeight: h + disp, pistonDisplacement: disp })).toBeCloseTo(cr, 10);
      }
    }
  });

  it("matches the engine model's head placement (headY differs by exactly the clearance)", () => {
    const L = computeLayout(CFR_F1);
    const cr = 7.3;
    const h = clearanceHeightAtTDC(cr, CFR_F1);
    const crBack = compressionRatioFromSnapshot(CFR_F1, { clearanceHeight: h, pistonDisplacement: 0 });
    expect(L.headY(crBack) - L.headY(18)).toBeCloseTo(h - clearanceHeightAtTDC(18, CFR_F1), 12);
  });

  it('agrees with the mock simulator kinematics at every crank angle', () => {
    const k = new MockKinematics(CFR_F1.geometry, 6.5);
    for (let th = -360; th < 360; th += 37) {
      const s = { clearanceHeight: k.clearanceHeight(th), pistonDisplacement: k.pistonDisplacement(th) };
      expect(compressionRatioFromSnapshot(CFR_F1, s)).toBeCloseTo(6.5, 9);
    }
  });

  it('handles a wrist-pin offset through the exact piston travel', () => {
    const spec: EngineSpec = { ...CFR_F1, geometry: { ...CFR_F1.geometry, pinOffset: 0.004 } };
    expect(pistonTravel(spec)).toBeGreaterThan(CFR_F1.geometry.stroke);
    const h = clearanceHeightAtTDC(8, spec);
    expect(compressionRatioFromSnapshot(spec, { clearanceHeight: h, pistonDisplacement: 0 })).toBeCloseTo(8, 10);
  });

  it('returns NaN for snapshots without usable geometry', () => {
    expect(compressionRatioFromSnapshot(CFR_F1, { clearanceHeight: 0, pistonDisplacement: 0 })).toBeNaN();
    expect(compressionRatioFromSnapshot(CFR_F1, { clearanceHeight: 0.01, pistonDisplacement: 0.02 })).toBeNaN();
    expect(compressionRatioFromSnapshot(CFR_F1, { clearanceHeight: NaN, pistonDisplacement: 0 })).toBeNaN();
  });
});

describe('FrameClock', () => {
  it('returns 0 on the first tick, then clamped wall seconds', () => {
    const c = new FrameClock(0.1);
    expect(c.tick(1000)).toBe(0);
    expect(c.tick(1016)).toBeCloseTo(0.016, 12);
    expect(c.tick(5000)).toBe(0.1); // background tab
    expect(c.tick(4990)).toBe(0); // non-monotonic timestamps
    c.reset();
    expect(c.tick(6000)).toBe(0);
  });
});

describe('RecentWindow', () => {
  it('uses the previous playback time while moving forward', () => {
    const w = new RecentWindow();
    expect(w.next(NaN)).toBe(-Infinity);
    expect(w.next(0)).toBe(-Infinity);
    expect(w.next(0.01)).toBe(0);
    expect(w.next(0.01)).toBe(0.01); // paused
    expect(w.next(0.02)).toBe(0.01);
  });

  it('re-opens the window after a backwards move or a reset', () => {
    const w = new RecentWindow();
    w.next(0.5);
    w.next(0.6);
    expect(w.next(0.55)).toBe(-Infinity); // step back
    expect(w.next(0.56)).toBe(0.55);
    expect(w.next(0.001)).toBe(-Infinity); // reset restarted the clock
    w.reset();
    expect(w.next(0.002)).toBe(-Infinity);
  });
});

describe('parseUrlOptions', () => {
  it('reads operating-point, playback and view knobs', () => {
    const o = parseUrlOptions('?cr=7.2&on=95&rpm=900&spark=20&phi=1.0&ts=0.01&paused=1&view=engine');
    expect(o.op).toEqual({
      compressionRatio: 7.2,
      fuel: { kind: 'PRF', octaneNumber: 95 },
      rpm: 900,
      sparkAdvanceDeg: 20,
      equivalenceRatio: 1,
    });
    expect(o.timeScale).toBe(0.01);
    expect(o.paused).toBe(true);
    expect(o.framing).toBe('engine');
  });

  it('ignores missing or invalid values', () => {
    const o = parseUrlOptions('?cr=abc&on=500&rpm=-3&ts=0&view=side&phi=');
    expect(o.op).toEqual({});
    expect(o.timeScale).toBeUndefined();
    expect(o.paused).toBeUndefined();
    expect(o.framing).toBeUndefined();
    expect(parseUrlOptions('').op).toEqual({});
    expect(parseUrlOptions('?paused=0').paused).toBe(false);
  });
});

describe('clampOperatingPoint', () => {
  it('clamps the CR to the spec range and leaves the rest alone', () => {
    const [lo, hi] = CFR_F1.geometry.compressionRatioRange;
    expect(clampOperatingPoint(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: 40 }).compressionRatio).toBe(hi);
    expect(clampOperatingPoint(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: 1.5 }).compressionRatio).toBe(lo);
    const same = clampOperatingPoint(CFR_F1, CFR_RON_CONDITIONS);
    expect(same).toEqual(CFR_RON_CONDITIONS);
    expect(same).not.toBe(CFR_RON_CONDITIONS);
  });
});

describe('framingDistanceScale', () => {
  it('pulls back only on viewports narrower than the reference aspect', () => {
    expect(framingDistanceScale(2.2)).toBe(1);
    expect(framingDistanceScale(FRAMING_REFERENCE_ASPECT)).toBe(1);
    expect(framingDistanceScale(FRAMING_REFERENCE_ASPECT / 2)).toBeCloseTo(2, 12);
    expect(framingDistanceScale(0.2)).toBe(2.5);
    expect(framingDistanceScale(NaN)).toBe(1);
    expect(framingDistanceScale(0)).toBe(1);
  });
});
