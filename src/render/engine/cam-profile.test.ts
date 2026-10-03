import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { ValveSpec } from '../../physics/core/engine-spec';
import {
  LiftProfileLearner, buildLobeProfile, followerPhi, lobeHeightAt, modelValveLift, wrapDeg720,
} from './cam-profile';

describe('wrapDeg720', () => {
  it('wraps into [−360, 360)', () => {
    expect(wrapDeg720(0)).toBe(0);
    expect(wrapDeg720(360)).toBe(-360);
    expect(wrapDeg720(375)).toBe(-345);
    expect(wrapDeg720(-361)).toBe(359);
    expect(wrapDeg720(-360)).toBe(-360);
    expect(wrapDeg720(1080 + 10)).toBe(-350);
  });
});

describe('modelValveLift', () => {
  const iv = CFR_F1.intakeValve;
  const ev = CFR_F1.exhaustValve;
  it('is zero outside the event and peaks at the middle', () => {
    expect(modelValveLift(iv, iv.openDeg - 1)).toBe(0);
    expect(modelValveLift(iv, iv.closeDeg + 1)).toBe(0);
    const mid = iv.openDeg + ((((iv.closeDeg - iv.openDeg) % 720) + 720) % 720) / 2;
    expect(modelValveLift(iv, mid)).toBeCloseTo(iv.maxLift, 12);
  });
  it('handles events that wrap through ±360 (exhaust closes after gas-exchange TDC)', () => {
    // exhaust: opens 140, closes −345 (≡ 375)
    expect(modelValveLift(ev, 300)).toBeGreaterThan(0);
    expect(modelValveLift(ev, -355)).toBeGreaterThan(0);
    expect(modelValveLift(ev, -340)).toBe(0);
    expect(modelValveLift(ev, 100)).toBe(0);
  });
  it('honours a non-zero timing lift threshold', () => {
    const v: ValveSpec = { ...iv, timingLiftThreshold: 0.0003 };
    expect(modelValveLift(v, v.openDeg)).toBeCloseTo(0.0003, 9);
    expect(modelValveLift(v, v.closeDeg)).toBeCloseTo(0.0003, 9);
  });
});

describe('LiftProfileLearner', () => {
  const f = (th: number) => modelValveLift(CFR_F1.exhaustValve, th);
  it('learns a lift curve from scattered samples', () => {
    const ln = new LiftProfileLearner();
    let seed = 1;
    const rnd = () => ((seed = (seed * 16807) % 2147483647) / 2147483647);
    expect(ln.stage()).toBe(0);
    for (let i = 0; i < 900; i++) {
      const th = -360 + 720 * rnd();
      ln.add(th, f(th));
    }
    expect(ln.stage()).toBeGreaterThanOrEqual(1);
    let maxErr = 0;
    for (let th = -360; th < 360; th += 0.37) maxErr = Math.max(maxErr, Math.abs(ln.sample(th) - f(th)));
    expect(maxErr).toBeLessThan(5e-5);
    // dense sweep fills every bin
    for (let th = -360; th < 360; th += 0.5) ln.add(th + 0.25, f(th + 0.25));
    expect(ln.stage()).toBe(2);
    maxErr = 0;
    for (let th = -360; th < 360; th += 0.37) maxErr = Math.max(maxErr, Math.abs(ln.sample(th) - f(th)));
    expect(maxErr).toBeLessThan(2e-6);
  });
  it('ignores isolated outliers once learned (e.g. a sample straddling the ±360° wrap)', () => {
    const ln = new LiftProfileLearner();
    for (let th = -360; th < 360; th += 0.5) ln.add(th + 0.25, f(th + 0.25));
    ln.add(0.5, 0.004); // bogus
    ln.add(-120.5, 0.003); // bogus
    expect(ln.sample(0.5)).toBeCloseTo(f(0.5), 9);
    expect(ln.sample(-120.5)).toBeCloseTo(f(-120.5), 9);
    expect(ln.resets).toBe(0);
  });

  it('restarts when the profile changes', () => {
    const ln = new LiftProfileLearner();
    for (let th = -360; th < 360; th += 0.5) ln.add(th, f(th));
    expect(ln.stage()).toBe(2);
    const g = (th: number) => 1.3 * f(th - 20);
    for (let k = 0; k < 3; k++) for (let th = -360; th < 360; th += 1.1) ln.add(th, g(th));
    expect(ln.resets).toBeGreaterThanOrEqual(1);
    for (let th = -360; th < 360; th += 0.5) ln.add(th, g(th));
    let maxErr = 0;
    for (let th = -360; th < 360; th += 0.7) maxErr = Math.max(maxErr, Math.abs(ln.sample(th) - g(th)));
    expect(maxErr).toBeLessThan(1e-5);
  });
});

describe('flat-follower cam lobe', () => {
  const ratio = 1.4;
  for (const v of [CFR_F1.intakeValve, CFR_F1.exhaustValve]) {
    it(`touches the tappet at the commanded lift for every crank angle (${v.openDeg}→${v.closeDeg})`, () => {
      const s = (th: number) => modelValveLift(v, th) / ratio;
      const p = buildLobeProfile(s, 0.024);
      expect(p.baseRadius).toBeGreaterThanOrEqual(0.024);
      let maxErr = 0;
      for (let th = -360; th < 360; th += 1.3) {
        const camAngle = (th * Math.PI) / 360; // cam turns +θ/2
        const h = lobeHeightAt(p.outline, camAngle);
        maxErr = Math.max(maxErr, Math.abs(h - (p.baseRadius + s(th))));
      }
      expect(maxErr).toBeLessThan(2e-5);
    });
  }

  it('grows the base circle when the lift curve is not realisable by a flat follower', () => {
    const v = CFR_F1.intakeValve;
    const s = (th: number) => modelValveLift(v, th) / ratio;
    const p = buildLobeProfile(s, 0.004, 720, 0.003);
    expect(p.baseRadius).toBeGreaterThan(0.004);
    // still exact everywhere after growing
    let maxErr = 0;
    for (let th = -360; th < 360; th += 2.1) {
      maxErr = Math.max(maxErr, Math.abs(lobeHeightAt(p.outline, (th * Math.PI) / 360) - (p.baseRadius + s(th))));
    }
    expect(maxErr).toBeLessThan(3e-5);
  });

  it('followerPhi points the lobe nose at the follower at peak lift', () => {
    const v = CFR_F1.intakeValve;
    const dur = (((v.closeDeg - v.openDeg) % 720) + 720) % 720;
    const peak = v.openDeg + dur / 2; // middle of the event (duration read from the spec)
    const phi = followerPhi(peak);
    const p = buildLobeProfile((th) => modelValveLift(v, th) / ratio, 0.024);
    // outline point farthest from the axis lies at polar angle ≈ phi
    let best = 0, bestAng = 0;
    for (let i = 0; i < p.outline.length; i += 2) {
      const r = Math.hypot(p.outline[i], p.outline[i + 1]);
      if (r > best) { best = r; bestAng = Math.atan2(p.outline[i + 1], p.outline[i]); }
    }
    const d = Math.atan2(Math.sin(bestAng - phi), Math.cos(bestAng - phi));
    expect(Math.abs(d)).toBeLessThan(0.02);
  });
});
