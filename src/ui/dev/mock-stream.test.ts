/** Sanity checks for the DEV-ONLY synthetic stream that drives the UI harness. */
import { describe, expect, it } from 'vitest';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { OperatingPoint } from '../../physics/core/operating-point';
import { CycleTraceStore } from '../trace-store';
import { SparkCapture } from '../spark-capture';
import { MockStream } from './mock-stream';

const OP: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: 101325,
  ambientTemperature: 298.15,
  relativeHumidity: 0.3,
  intakeMixtureTemperature: 325.15,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: 1.05,
  sparkAdvanceDeg: 13,
  dwellTime: 3e-3,
  compressionRatio: 7,
  egrFraction: 0,
  coolantTemperature: 373.15,
};

describe('MockStream (dev harness)', () => {
  it('produces a consistent, plausible stream the UI can digest', () => {
    const m = new MockStream(CFR_F1, OP);
    const snaps = [];
    for (let i = 0; i < 720 * 2 * 4; i++) snaps.push(m.step(0.25)); // two cycles
    let prevT = -1;
    for (const s of snaps) {
      expect(s.t).toBeGreaterThan(prevT);
      prevT = s.t;
      expect(s.thetaDeg).toBeGreaterThanOrEqual(-360);
      expect(s.thetaDeg).toBeLessThan(360);
      expect(s.pressure).toBeGreaterThan(0);
      expect(s.volume).toBeGreaterThan(0);
      expect(s.massFractionBurned).toBeGreaterThanOrEqual(0);
      expect(s.massFractionBurned).toBeLessThanOrEqual(1);
    }
    const peak = Math.max(...snaps.map((s) => s.pressure));
    expect(peak).toBeGreaterThan(15e5);
    expect(peak).toBeLessThan(120e5);
    expect(snaps[snaps.length - 1].cycle).toBe(2);
    const phases = new Set(snaps.map((s) => s.spark.phase));
    for (const p of ['charging', 'arc', 'glow', 'done'] as const) expect(phases.has(p)).toBe(true);
    const cyc = m.drainCycles();
    expect(cyc.length).toBe(2);
    expect(cyc[0].imepNet).toBeGreaterThan(0);

    const store = new CycleTraceStore(4);
    store.ingest(snaps, Infinity);
    expect(store.traces.length).toBe(3);
    const cap = new SparkCapture();
    cap.ingest(snaps, Infinity);
    expect(cap.display).not.toBeNull();
    expect(Number.isFinite(cap.display!.tBreakdown)).toBe(true);
  });
});
