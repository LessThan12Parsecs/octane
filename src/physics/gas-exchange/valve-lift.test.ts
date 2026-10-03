import { describe, expect, it } from 'vitest';
import type { ValveSpec } from '../core/engine-spec';
import {
  CHOI2018_EXHAUST_LIFT_START_DEG,
  CHOI2018_EXHAUST_LIFT_ZERO_LASH,
  CHOI2018_INTAKE_LIFT_START_DEG,
  CHOI2018_INTAKE_LIFT_ZERO_LASH,
} from '../engines/cfr-valve-lift';
import {
  POLYDYNE_2_10_18_26,
  polydyneCam,
  polydyneCoefficients,
  ValveLiftProfile,
  valveLift,
  valveLiftProfile,
  valveLiftRate,
  wrapDeg720,
} from './valve-lift';

const IN = 0.0254;
/** CFR-like timing (numbers passed as inputs; the module never hard-codes them). */
const intake = { maxLift: 0.246 * IN, openDeg: -350, closeDeg: -146, timingLiftThreshold: 0 };
const exhaust = { maxLift: 0.246 * IN, openDeg: 140, closeDeg: -345, timingLiftThreshold: 0 };

describe('polydyne cam polynomial', () => {
  it('2-10-18-26 coefficients are the exact rationals (−195, 117, −65, 15)/128', () => {
    const c = polydyneCoefficients([2, 10, 18, 26]);
    const exact = [-195 / 128, 117 / 128, -65 / 128, 15 / 128];
    for (let i = 0; i < 4; i++) expect(c[i]).toBeCloseTo(exact[i], 14);
    expect(POLYDYNE_2_10_18_26.coeffs.length).toBe(4);
  });

  it('closes onto zero lift with zero 1st–(n−1)th derivatives (C^{n−1})', () => {
    for (const powers of [
      [2, 4],
      [2, 6, 10, 14],
      [2, 10, 18, 26],
      [2, 6, 10, 14, 18],
    ]) {
      const cam = polydyneCam(powers);
      const n = powers.length;
      for (let k = 0; k < n; k++) {
        let s = k === 0 ? 1 : 0;
        for (let i = 0; i < n; i++) {
          let f = 1;
          for (let j = 0; j < k; j++) f *= powers[i] - j;
          s += cam.coeffs[i] * f;
        }
        expect(Math.abs(s)).toBeLessThan(1e-11);
      }
    }
    expect(() => polydyneCoefficients([2, 3])).toThrow();
  });

  it('shape metrics of the default cam: monotone flanks, nose deceleration, fullness', () => {
    const p = new ValveLiftProfile({ maxLift: 1, openDeg: -90, closeDeg: 90, timingLiftThreshold: 0 });
    let prev = 1;
    let area = 0;
    const N = 20000;
    let maxAcc = -Infinity;
    for (let i = 1; i <= N; i++) {
      const x = i / N;
      const s = p.shape(x);
      expect(s).toBeLessThanOrEqual(prev + 1e-15);
      expect(s).toBeGreaterThanOrEqual(-1e-14); // rounding at x = 1 (lift() clamps at 0)
      prev = s;
      area += s / N;
      maxAcc = Math.max(maxAcc, p.shapeD2(x));
    }
    expect(p.shapeD2(0)).toBeCloseTo(-2 * (195 / 128), 12);
    expect(area).toBeCloseTo(0.553, 2);
    expect(maxAcc / -p.shapeD2(0)).toBeGreaterThan(3);
    expect(maxAcc / -p.shapeD2(0)).toBeLessThan(3.5);
  });
});

describe('ValveLiftProfile', () => {
  it('intake (no wrap): seated outside [open, close], max lift at the centre', () => {
    const p = new ValveLiftProfile(intake);
    expect(p.centerDeg).toBeCloseTo(-248, 12);
    expect(p.halfDurationDeg).toBeCloseTo(102, 12);
    expect(p.lift(-350)).toBe(0);
    expect(p.lift(-146)).toBe(0);
    expect(p.lift(-351)).toBe(0);
    expect(p.lift(0)).toBe(0);
    expect(p.lift(-349.9)).toBeGreaterThan(0);
    expect(p.lift(-248)).toBeCloseTo(intake.maxLift, 15);
    expect(p.liftRate(-248)).toBeCloseTo(0, 15);
    expect(p.isOpen(-300)).toBe(true);
    expect(p.isOpen(100)).toBe(false);
  });

  it('exhaust closing after the gas-exchange TDC wraps across ±360°', () => {
    const p = new ValveLiftProfile(exhaust);
    expect(p.centerDeg).toBeCloseTo(257.5, 12);
    expect(p.seatCloseDeg).toBeCloseTo(-345, 12);
    expect(p.lift(140)).toBe(0);
    expect(p.lift(-345)).toBe(0);
    expect(p.lift(-344)).toBe(0);
    expect(p.lift(-350)).toBeGreaterThan(0); // −350 ≡ 370 < 375
    expect(p.lift(359.99)).toBeGreaterThan(0);
    expect(p.lift(-360)).toBeGreaterThan(0);
    // periodicity
    for (const th of [-355, -300, 150, 257.5, 300, 359]) {
      expect(p.lift(th + 720)).toBeCloseTo(p.lift(th), 15);
      expect(p.lift(th - 720)).toBeCloseTo(p.lift(th), 15);
    }
    // overlap with the intake: both open between IVO (−350) and EVC (−345)
    const pi = new ValveLiftProfile(intake);
    expect(pi.lift(-347) > 0 && p.lift(-347) > 0).toBe(true);
    expect(wrapDeg720(375)).toBe(-345);
  });

  it('honours timingLiftThreshold: L(open) = L(close) = threshold', () => {
    const thr = 0.006 * IN; // SAE 0.006 in timing point (Heywood 1988 §6.3)
    const spec = { ...intake, timingLiftThreshold: thr };
    const p = new ValveLiftProfile(spec);
    expect(p.lift(spec.openDeg)).toBeCloseTo(thr, 14);
    expect(p.lift(spec.closeDeg)).toBeCloseTo(thr, 14);
    expect(p.halfDurationDeg).toBeGreaterThan(102);
    expect(p.lift(p.seatOpenDeg + 1e-9)).toBeLessThan(1e-12);
    expect(() => new ValveLiftProfile({ ...intake, timingLiftThreshold: intake.maxLift })).toThrow();
  });

  it('is C² (and C³ at the seat): analytic derivatives match finite differences', () => {
    for (const spec of [intake, exhaust, { ...intake, timingLiftThreshold: 1e-4 }]) {
      const p = new ValveLiftProfile(spec);
      const h = 1e-4;
      let maxErr1 = 0;
      let maxErr2 = 0;
      for (let th = -360; th < 360; th += 0.37) {
        const d1 = (p.lift(th + h) - p.lift(th - h)) / (2 * h);
        const d2 = (p.liftRate(th + h) - p.liftRate(th - h)) / (2 * h);
        maxErr1 = Math.max(maxErr1, Math.abs(d1 - p.liftRate(th)));
        maxErr2 = Math.max(maxErr2, Math.abs(d2 - p.liftAccel(th)));
      }
      // scales: rate ~ L/H ≈ 6e-5 m/deg, accel ~ L/H² ≈ 6e-6 m/deg²
      expect(maxErr1).toBeLessThan(1e-11);
      expect(maxErr2).toBeLessThan(1e-11);
      // seat: velocity, acceleration and jerk vanish (acceleration ∝ ε² just off the seat)
      const e = 1e-3;
      const a1 = Math.abs(p.liftAccel(p.seatOpenDeg + e));
      const a2 = Math.abs(p.liftAccel(p.seatOpenDeg + 2 * e));
      expect(a2 / a1).toBeCloseTo(4, 2);
      expect(Math.abs(p.liftRate(p.seatCloseDeg - e))).toBeLessThan(1e-15);
    }
  });

  it('valve lash: valve lift = cam rise − lash, finite seat velocity, timing met at the valve', () => {
    const lash = 0.008 * IN;
    const spec = { ...intake, maxLift: 0.23 * IN, lash };
    const p = new ValveLiftProfile(spec);
    expect(p.camMaxLift).toBeCloseTo(0.238 * IN, 15);
    expect(p.seatOpenDeg).toBeCloseTo(spec.openDeg, 9);
    expect(p.seatCloseDeg).toBeCloseTo(spec.closeDeg, 9);
    expect(p.lift(p.centerDeg)).toBeCloseTo(spec.maxLift, 15);
    expect(p.lift(spec.openDeg - 1e-6)).toBe(0);
    expect(p.isOpen(spec.openDeg - 1e-6)).toBe(false);
    expect(p.isOpen(spec.openDeg + 1e-6)).toBe(true);
    // the cam is already moving when the valve leaves the seat
    const v0 = p.liftRate(spec.openDeg + 1e-9);
    expect(v0).toBeGreaterThan(1e-6); // m/deg (CFR-like: ≈ 0.03 mm/deg)
    expect(p.lift(spec.openDeg + 1e-3) / 1e-3).toBeCloseTo(v0, 8);
    expect(p.liftRate(spec.closeDeg - 1e-9)).toBeLessThan(-1e-6);
    for (const th of [-340, -300, -248, -200, -160]) expect(p.camLift(th) - lash).toBeCloseTo(p.lift(th), 15);
    expect(p.camLift(spec.openDeg - 0.5)).toBeGreaterThan(0); // on the ramp, clearance not yet taken up
    // with a timing threshold on top of the lash
    const thr = 0.006 * IN;
    const q = new ValveLiftProfile({ ...spec, timingLiftThreshold: thr });
    expect(q.lift(spec.openDeg)).toBeCloseTo(thr, 14);
    expect(q.lift(spec.closeDeg)).toBeCloseTo(thr, 14);
    expect(q.lift(q.seatOpenDeg + 1e-9)).toBeLessThan(1e-12);
    expect(q.lift(q.seatOpenDeg - 1e-9)).toBe(0);
    // lash = 0 reduces to the original profile
    const r0 = new ValveLiftProfile({ ...intake, lash: 0 });
    const r1 = new ValveLiftProfile(intake);
    for (const th of [-349, -300, -200, -147]) expect(r0.lift(th)).toBe(r1.lift(th));
    // cache tracks lash edits
    const s2 = { ...spec };
    const a = valveLiftProfile(s2);
    s2.lash = 0;
    expect(valveLiftProfile(s2)).not.toBe(a);
    expect(valveLiftProfile(s2).lash).toBe(0);
  });

  it('CFR: with the 0.008 in running clearance the cam reproduces the measured lift (Choi et al. 2018)', () => {
    // Oracle: zero-lash lift of the ANL CFR engine, Choi et al. SAE 2018-01-0848 Fig. 3 (digitised
    // table, engines/cfr-valve-lift.ts); valve events at the 0.008 in hot clearance from their
    // Table 2 (IVO 376 ≡ −344, IVC −152, EVO 141, EVC 373 ≡ −347 CAD), lash per ASTM D2699 §10.3.2.1.
    const lash = 0.008 * IN;
    const measured = (table: readonly number[], start: number, th: number): number => {
      let u = (th - start) % 720;
      if (u < 0) u += 720;
      if (!(u < table.length - 1)) return 0;
      const i = Math.floor(u);
      const l0 = table[i] + (u - i) * (table[i + 1] - table[i]);
      return l0 > lash ? l0 - lash : 0;
    };
    const cases = [
      { kind: 'intake', table: CHOI2018_INTAKE_LIFT_ZERO_LASH, start: CHOI2018_INTAKE_LIFT_START_DEG, openDeg: -344, closeDeg: -152, rmsTol: 0.12e-3, intTol: 0.02 },
      { kind: 'exhaust', table: CHOI2018_EXHAUST_LIFT_ZERO_LASH, start: CHOI2018_EXHAUST_LIFT_START_DEG, openDeg: 141, closeDeg: -347, rmsTol: 0.35e-3, intTol: 0.08 },
    ] as const;
    for (const c of cases) {
      const maxLift = Math.max(...c.table) - lash;
      const base = { maxLift, openDeg: c.openDeg, closeDeg: c.closeDeg, timingLiftThreshold: 0 };
      const stats = (p: ValveLiftProfile) => {
        let e2 = 0;
        let n = 0;
        let iP = 0;
        let iM = 0;
        let early = 0; // lift integral over the first/last 20° of the event (IVO/IVC, EVO/EVC flow)
        let earlyM = 0;
        for (let th = -360; th < 360; th += 0.25) {
          const lp = p.lift(th);
          const lm = measured(c.table, c.start, th);
          iP += lp;
          iM += lm;
          if (lp > 0 || lm > 0) {
            e2 += (lp - lm) ** 2;
            n++;
          }
          const dOpen = wrapDeg720(th - c.openDeg);
          const dClose = wrapDeg720(c.closeDeg - th);
          if ((dOpen >= 0 && dOpen < 20) || (dClose >= 0 && dClose < 20)) {
            early += lp;
            earlyM += lm;
          }
        }
        return { rms: Math.sqrt(e2 / n), integral: iP / iM - 1, early: early / earlyM };
      };
      const withLash = stats(new ValveLiftProfile({ ...base, lash }));
      const noLash = stats(new ValveLiftProfile(base));
      expect(withLash.rms).toBeLessThan(c.rmsTol);
      expect(Math.abs(withLash.integral)).toBeLessThan(c.intTol);
      expect(withLash.early).toBeGreaterThan(0.7);
      expect(withLash.early).toBeLessThan(1.5);
      // the zero-velocity seat departure (lash ignored) starves the valve events of flow area
      expect(noLash.early).toBeLessThan(0.4);
      expect(noLash.rms).toBeGreaterThan(2 * withLash.rms);
    }
  });

  it('functional API caches per spec and tracks spec edits', () => {
    const spec: ValveSpec = {
      count: 1,
      headDiameter: 0.037,
      seatInnerDiameter: 0.034,
      seatAngle: Math.PI / 4,
      stemDiameter: 0.0087,
      maxLift: 0.00625,
      openDeg: -350,
      closeDeg: -146,
      timingLiftThreshold: 0,
      position: [-0.02, 0],
      shroudArcDeg: 180,
      shroudDirection: Math.PI / 2,
    };
    const a = valveLiftProfile(spec);
    expect(valveLiftProfile(spec)).toBe(a);
    expect(valveLift(spec, -248)).toBeCloseTo(0.00625, 15);
    expect(valveLiftRate(spec, -300)).toBeGreaterThan(0);
    expect(valveLiftRate(spec, -200)).toBeLessThan(0);
    spec.maxLift = 0.005;
    expect(valveLift(spec, -248)).toBeCloseTo(0.005, 15);
    expect(valveLiftProfile(spec)).not.toBe(a);
  });
});
