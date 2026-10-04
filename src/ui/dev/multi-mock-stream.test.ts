/** Sanity checks for the DEV-ONLY multi-cylinder stream that drives the UI harness. */
import { describe, expect, it } from 'vitest';
import { cylinderAngleDeg } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { MODEL_T, MODEL_T_CRUISE } from '../../physics/engines/model-t';
import { CylinderViewPool } from '../cylinder-view';
import { wrapCycleDeg } from '../engine-cycle';
import { SparkCapture } from '../spark-capture';
import { MultiMockStream } from './multi-mock-stream';

describe('MultiMockStream (dev harness)', () => {
  it('emits 4-cylinder snapshots at the firing offsets, one shower per contact, and engine summaries', () => {
    const m = new MultiMockStream(MODEL_T, { ...MODEL_T_CRUISE, speedMode: 'fixed' });
    const snaps: EngineSnapshot[] = [];
    while (m.cycle < 2) snaps.push(m.step(m.inSparkWindow() ? 0.1 : 0.5));
    while (m.inSparkWindow()) snaps.push(m.step(0.1)); // end between contacts
    for (let k = 0; k < 10; k++) snaps.push(m.step(0.5));
    let prevT = -1;
    for (const s of snaps) {
      expect(s.t).toBeGreaterThan(prevT);
      prevT = s.t;
      expect(s.cylinders).toHaveLength(4);
      // (mod 720: the cores accumulate θ separately and may sit on either side of the wrap)
      for (let i = 0; i < 4; i++) expect(wrapCycleDeg(s.cylinders![i].thetaDeg - cylinderAngleDeg(MODEL_T, i, s.thetaDeg) + 1e-6)).toBeCloseTo(1e-6, 6);
      expect(s.cylinders![0].pressure).toBe(s.pressure);
    }
    // Every cylinder's timer contact makes one event of several sparks.
    const pool = new CylinderViewPool();
    for (let i = 0; i < 4; i++) {
      const cap = new SparkCapture();
      cap.ingest(pool.map(snaps, i), Infinity);
      const ev = cap.display!;
      expect(ev.shower).toBe(true);
      expect(ev.sparkCount).toBeGreaterThan(1);
      expect(ev.firstSparkDeg).toBeGreaterThan(-MODEL_T_CRUISE.sparkAdvanceDeg);
    }
    const cycles = m.drainCycles();
    expect(new Set(cycles.map((c) => c.cylinder))).toEqual(new Set([0, 1, 2, 3]));
    const withEngine = cycles.filter((c) => c.engine);
    expect(withEngine.length).toBeGreaterThan(0);
    expect(withEngine.every((c) => c.cylinder === 0)).toBe(true);
    expect(snaps.some((s) => s.firingCylinder !== undefined && s.firingCylinder >= 0)).toBe(true);
    expect(snaps.some((s) => (s.magnetoEmf ?? 0) > 1)).toBe(true);
  });
});
