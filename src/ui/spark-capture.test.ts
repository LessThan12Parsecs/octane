import { describe, expect, it } from 'vitest';
import type { EngineSnapshot, SparkPhase } from '../physics/core/snapshot';
import { SparkCapture, sparkPhaseDurations, SPARK_PHASE_CODE } from './spark-capture';
import { makeSnapshot } from './test-helpers';

/** A spark event sampled every `dt`: [phase, duration] segments starting at t0. */
function sparkRun(t0: number, dt: number, segments: [SparkPhase, number][], cycle = 0): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  let t = t0;
  for (const [phase, dur] of segments) {
    const n = Math.round(dur / dt);
    for (let i = 0; i < n; i++) {
      const s = makeSnapshot({ t, cycle, thetaDeg: -13 + (t - t0) * 3600 });
      s.spark = {
        phase,
        primaryCurrent: phase === 'charging' ? 6 : 0,
        secondaryVoltage: phase === 'breakdown' ? 12e3 : phase === 'arc' ? 100 : phase === 'glow' ? 500 : 0,
        secondaryCurrent: phase === 'arc' ? 0.2 : phase === 'glow' ? 0.05 : 0,
        energyDelivered: phase === 'done' ? 0.03 : 0,
        breakdownVoltage: 11e3,
      };
      out.push(s);
      t += dt;
    }
  }
  return out;
}

const DT = 10e-6;
const EVENT: [SparkPhase, number][] = [
  ['off', 100e-6],
  ['charging', 3e-3],
  ['breakdown', 10e-6],
  ['arc', 50e-6],
  ['glow', 1.5e-3],
  ['done', 1e-3],
];

describe('SparkCapture', () => {
  it('records a full event and closes it after the tail', () => {
    const cap = new SparkCapture({ tail: 0.5e-3 });
    const run = sparkRun(0, DT, EVENT);
    expect(cap.display).toBeNull();
    cap.ingest(run, Infinity);
    const ev = cap.display!;
    expect(ev).not.toBeNull();
    expect(ev.complete).toBe(true);
    expect(ev.tStart).toBeCloseTo(100e-6, 12);
    expect(ev.tBreakdown).toBeCloseTo(3.1e-3, 9);
    expect(ev.tEnd).toBeCloseTo(3.1e-3 + 1.56e-3, 9);
    // Stops recording ~tail after the end.
    expect(ev.t.data[ev.length - 1] - ev.tEnd).toBeLessThan(0.5e-3 + 2 * DT);
    expect(ev.vSec.data[ev.t.length - 1]).toBe(0);
    expect(Math.max(...ev.vSec.view())).toBeCloseTo(12, 12); // kV
    expect(Math.max(...ev.iSec.view())).toBeCloseTo(200, 12); // mA
    expect(ev.energy).toBeCloseTo(0.03, 12);

    const d = sparkPhaseDurations(ev);
    expect(d.dwell).toBeCloseTo(3e-3, 9);
    expect(d.arc).toBeCloseTo(50e-6, 9);
    expect(d.glow).toBeCloseTo(1.5e-3, 9);
    expect(d.discharge).toBeCloseTo(1.56e-3, 9);
  });

  it('keeps showing the last event during the next dwell, then switches at breakdown', () => {
    const cap = new SparkCapture();
    const first = sparkRun(0, DT, EVENT, 0);
    cap.ingest(first, Infinity);
    const shown = cap.display!;
    const tNext = 0.2;
    const second = sparkRun(tNext, DT, EVENT, 1);
    const iBd = second.findIndex((s) => s.spark.phase === 'breakdown');
    cap.ingest(second.slice(0, iBd - 5), Infinity); // still charging
    expect(cap.display!.cycle).toBe(0);
    expect(cap.live!.cycle).toBe(1);
    expect(cap.display).toBe(shown);
    cap.ingest(second.slice(iBd - 5, iBd + 3), Infinity);
    expect(cap.display!.cycle).toBe(1);
    expect(cap.display!.complete).toBe(false);
  });

  it('only ingests up to the playhead', () => {
    const cap = new SparkCapture();
    const run = sparkRun(0, DT, EVENT);
    const iBd = run.findIndex((s) => s.spark.phase === 'breakdown');
    cap.ingest(run, run[iBd - 1].t);
    expect(cap.display).toBeNull();
    expect(cap.live!.length).toBe(iBd - 10);
  });

  it('keeps a failed (no breakdown) event', () => {
    const cap = new SparkCapture();
    cap.ingest(
      sparkRun(0, DT, [
        ['charging', 1e-3],
        ['done', 1e-3],
      ]),
      Infinity,
    );
    const ev = cap.display!;
    expect(ev.complete).toBe(true);
    expect(ev.tBreakdown).toBeNaN();
    expect(ev.tRef).toBe(ev.tStart);
  });

  it('closes an event that never reports done', () => {
    const cap = new SparkCapture({ maxAfterBreakdown: 1e-3 });
    cap.ingest(
      sparkRun(0, DT, [
        ['charging', 1e-3],
        ['arc', 3e-3],
      ]),
      Infinity,
    );
    expect(cap.display!.complete).toBe(true);
  });

  it('phase codes are distinct', () => {
    expect(new Set(Object.values(SPARK_PHASE_CODE)).size).toBe(6);
  });
});
