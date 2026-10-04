/** SparkCapture with trembler-coil showers: one event per timer contact. */
import { describe, expect, it } from 'vitest';
import type { EngineSnapshot, SparkPhase } from '../physics/core/snapshot';
import { SparkCapture } from './spark-capture';
import { makeSnapshot } from './test-helpers';

interface Seg {
  phase: SparkPhase;
  dur: number;
  timer: boolean;
}

const DT = 10e-6;
const RPM = 1000;

/**
 * A trembler stream sampled every DT: segments of [phase, duration, timer contact]. `count` reports
 * breakdownCount / firstSparkDeg like the physics does (or omits them, as a minimal simulator might).
 */
function stream(t0: number, segs: Seg[], opts: { count?: boolean; cycle?: number; theta0?: number } = {}): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  let t = t0;
  let n = 0;
  let prevDischarge = false;
  let first = NaN;
  const th0 = opts.theta0 ?? -30;
  for (const seg of segs) {
    const steps = Math.round(seg.dur / DT);
    for (let i = 0; i < steps; i++) {
      const theta = th0 + (t - t0) * 6 * RPM;
      const discharge = seg.phase === 'breakdown' || seg.phase === 'arc' || seg.phase === 'glow';
      if (discharge && !prevDischarge) {
        n++;
        if (!Number.isFinite(first)) first = theta;
      }
      prevDischarge = discharge;
      const s = makeSnapshot({ t, cycle: opts.cycle ?? 0, thetaDeg: theta, rpm: RPM });
      s.spark = {
        phase: seg.phase,
        primaryCurrent: seg.phase === 'charging' ? 4 : 0,
        secondaryVoltage: seg.phase === 'breakdown' ? 8e3 : seg.phase === 'arc' ? 100 : 0,
        secondaryCurrent: seg.phase === 'arc' ? 0.1 : 0,
        energyDelivered: 0.002 * n,
        breakdownVoltage: 7e3,
        timerClosed: seg.timer,
        pointsOpen: seg.phase !== 'charging',
      };
      if (opts.count) {
        s.spark.breakdownCount = n;
        s.spark.firstSparkDeg = first;
      }
      out.push(s);
      t += DT;
    }
  }
  return out;
}

/** Timer make, coil build-up, `sparks` buzz cycles, timer break, quiet. */
function contact(sparks: number, opts: { doneBetween?: boolean } = {}): Seg[] {
  const segs: Seg[] = [
    { phase: 'off', dur: 100e-6, timer: false },
    { phase: 'charging', dur: 3e-3, timer: true },
  ];
  for (let k = 0; k < sparks; k++) {
    segs.push({ phase: 'breakdown', dur: 10e-6, timer: true });
    segs.push({ phase: 'arc', dur: 60e-6, timer: true });
    segs.push({ phase: 'glow', dur: 300e-6, timer: true });
    if (opts.doneBetween) segs.push({ phase: 'done', dur: 200e-6, timer: true });
    if (k < sparks - 1) segs.push({ phase: 'charging', dur: 1.6e-3, timer: true });
  }
  segs.push({ phase: 'done', dur: 1.5e-3, timer: false });
  return segs;
}

describe('SparkCapture: trembler showers', () => {
  for (const count of [true, false]) {
    it(`groups a 6-spark shower into one event (${count ? 'with' : 'without'} breakdownCount/firstSparkDeg)`, () => {
      const cap = new SparkCapture();
      const run = stream(0, contact(6), { count });
      cap.ingest(run, Infinity);
      const ev = cap.display!;
      expect(ev.complete).toBe(true);
      expect(ev.shower).toBe(true);
      expect(ev.sparkCount).toBe(6);
      const firstBd = run.find((s) => s.spark.phase === 'breakdown')!;
      expect(ev.tBreakdown).toBe(firstBd.t);
      expect(ev.firstSparkDeg).toBeCloseTo(firstBd.thetaDeg, 9);
      expect(ev.tStart).toBeCloseTo(100e-6, 12); // timer make
      // ends when the contact opens (and the last glow has finished), then the tail
      const open = run.find((s) => s.t > ev.tBreakdown && s.spark.timerClosed === false)!;
      expect(ev.tEnd).toBe(open.t);
      expect(ev.energy).toBeCloseTo(0.012, 12);
    });
  }

  it("a 'done' between buzzes while the contact is closed does not split the shower", () => {
    const cap = new SparkCapture();
    cap.ingest(stream(0, contact(4, { doneBetween: true }), { count: true }), Infinity);
    expect(cap.display!.sparkCount).toBe(4);
    expect(cap.display!.complete).toBe(true);
  });

  it('two timer contacts are two events (next cycle of the cylinder)', () => {
    const cap = new SparkCapture();
    const a = stream(0, contact(3), { count: true, cycle: 0 });
    const b = stream(0.12, contact(5), { count: true, cycle: 1 });
    cap.ingest(a, Infinity);
    const first = cap.display!;
    expect(first.cycle).toBe(0);
    expect(first.sparkCount).toBe(3);
    cap.ingest(b, Infinity);
    expect(cap.display!.cycle).toBe(1);
    expect(cap.display!.sparkCount).toBe(5);
  });

  it('shows the live shower once it fires, the previous one before', () => {
    const cap = new SparkCapture();
    cap.ingest(stream(0, contact(3), { count: true, cycle: 0 }), Infinity);
    const b = stream(0.12, contact(5), { count: true, cycle: 1 });
    const iBd = b.findIndex((s) => s.spark.phase === 'breakdown');
    cap.ingest(b.slice(0, iBd - 1), Infinity);
    expect(cap.display!.cycle).toBe(0);
    expect(cap.live!.cycle).toBe(1);
    cap.ingest(b.slice(iBd - 1, iBd + 400), Infinity);
    expect(cap.display!.cycle).toBe(1);
    expect(cap.display!.complete).toBe(false);
    expect(cap.display!.sparkCount).toBeGreaterThanOrEqual(2);
  });

  it('a contact without breakdown is a failed event, closed when the contact opens', () => {
    const cap = new SparkCapture();
    cap.ingest(
      stream(0, [
        { phase: 'off', dur: 100e-6, timer: false },
        { phase: 'charging', dur: 2e-3, timer: true },
        { phase: 'off', dur: 1e-3, timer: false },
      ]),
      Infinity,
    );
    const ev = cap.display!;
    expect(ev.complete).toBe(true);
    expect(ev.sparkCount).toBe(0);
    expect(ev.tBreakdown).toBeNaN();
    expect(ev.firstSparkDeg).toBeNaN();
  });

  it('a contact that never reports opening is closed after maxShower', () => {
    const cap = new SparkCapture({ maxShower: 5e-3 });
    const segs = contact(8);
    segs.pop(); // no timer break
    cap.ingest(stream(0, segs, { count: true }), Infinity);
    expect(cap.display!.complete).toBe(true);
  });
});
