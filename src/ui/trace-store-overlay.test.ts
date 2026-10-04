/** The all-cylinders pressure store (engine angle) and per-cylinder stores fed with cylinder views. */
import { describe, expect, it } from 'vitest';
import { cylinderAngleDeg } from '../physics/core/engine-spec';
import type { EngineSnapshot } from '../physics/core/snapshot';
import { MODEL_T } from '../physics/engines/model-t';
import { createEmptySnapshot } from '../app/sim-client';
import { CylinderViewPool } from './cylinder-view';
import { CH, CycleTraceStore, cylinderPressureSampler, overlayPressureChannel } from './trace-store';

const RPM = 1000;

/** 4-cylinder stream at 1° from engine θ = −360, cycle 0; each cylinder counts its own cycles. */
function run(n: number): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  const cyc = [0, 0, 0, 0];
  const prev = [0, 1, 2, 3].map((i) => cylinderAngleDeg(MODEL_T, i, -360));
  let cycle = 0;
  for (let k = 0; k < n; k++) {
    const raw = -360 + k;
    const th = cylinderAngleDeg(MODEL_T, 0, raw);
    if (k > 0 && th === -360) cycle++;
    const s = createEmptySnapshot(4);
    s.t = k / (6 * RPM);
    s.rpm = RPM;
    s.cycle = cycle;
    s.thetaDeg = th;
    s.cylinders!.forEach((c, i) => {
      const lt = cylinderAngleDeg(MODEL_T, i, raw);
      if (k > 0 && lt < prev[i]) cyc[i]++;
      prev[i] = lt;
      c.thetaDeg = lt;
      c.cycle = cyc[i];
      // a pressure peak at each cylinder's own firing TDC
      c.pressure = 1e5 + 20e5 * Math.exp(-((lt / 40) ** 2));
    });
    s.pressure = s.cylinders![0].pressure;
    out.push(s);
  }
  return out;
}

describe('cylinderPressureSampler store', () => {
  it('stores every cylinder’s pressure against the engine angle, split at the engine wrap', () => {
    const store = new CycleTraceStore(4, 256, cylinderPressureSampler(4));
    const snaps = run(1000);
    expect(store.ingest(snaps, Infinity)).toBe(1000);
    expect(store.traces.map((t) => t.cycle)).toEqual([0, 1]);
    const tr = store.traces[0];
    expect(tr.channels).toBe(6);
    expect(tr.n).toBe(720);
    // Cylinder 4 (offset 360) peaks at engine θ = 0 − 360 = −360 … i.e. at the start; cylinder 2 at 180.
    const th = tr.col(CH.theta);
    const peakAt = (i: number) => {
      const p = tr.col(overlayPressureChannel(i));
      let k = 0;
      for (let j = 1; j < p.length; j++) if (p[j] > p[k]) k = j;
      return th[k];
    };
    expect(peakAt(0)).toBe(0);
    expect(peakAt(1)).toBe(180);
    expect(peakAt(2)).toBe(-180); // cylinder 3: offset 540 ≡ −180
    expect(peakAt(3)).toBe(-360); // cylinder 4: offset 360
    expect(tr.col(overlayPressureChannel(1))[0]).toBeCloseTo(1 + 20 * Math.exp(-((180 / 40) ** 2)), 9);
  });

  it('single-cylinder snapshots fill cylinder 1 and leave the others NaN', () => {
    const store = new CycleTraceStore(2, 16, cylinderPressureSampler(2));
    const s = createEmptySnapshot();
    s.t = 0.001;
    s.pressure = 3e5;
    store.ingest([s], Infinity);
    expect(store.traces[0].cols[2][0]).toBe(3);
    expect(store.traces[0].cols[3][0]).toBeNaN();
  });
});

describe('per-cylinder stores fed with cylinder views', () => {
  it("split each cylinder's traces at that cylinder's own wrap, in its local angle", () => {
    const pool = new CylinderViewPool();
    const stores = [0, 1, 2, 3].map(() => new CycleTraceStore(4, 256));
    const snaps = run(1100);
    for (let i = 0; i < 4; i++) stores[i].ingest(pool.map(snaps, i), Infinity);
    // Cylinder 2 starts at local 180 (cycle 0), wraps after 180°: traces [0, 1, 2?]
    const lens = stores.map((st) => st.traces.map((t) => t.n));
    expect(lens[0]).toEqual([720, 380]);
    expect(lens[1]).toEqual([180, 720, 200]);
    expect(lens[2]).toEqual([540, 560]);
    expect(lens[3]).toEqual([360, 720, 20]);
    for (const st of stores) {
      for (const tr of st.traces) {
        const th = tr.col(CH.theta);
        for (let k = 1; k < th.length; k++) expect(th[k]).toBeGreaterThan(th[k - 1]);
        expect(th[th.length - 1]).toBeLessThan(360);
      }
    }
    // each cylinder's own peak lands at its local TDC
    const p = stores[2].traces[1].col(CH.p);
    const th = stores[2].traces[1].col(CH.theta);
    let k = 0;
    for (let j = 1; j < p.length; j++) if (p[j] > p[k]) k = j;
    expect(th[k]).toBe(0);
  });
});
