import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import { createEmptySnapshot } from '../app/sim-client';
import type { FromWorker, SimulatorOptions } from './protocol';
import { HOST_MAX_BATCH, type HostPort, SimWorkerHost } from './sim.worker';
import type { SimulatorLike } from './simulator-like';

const OP: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: 101325,
  ambientTemperature: 298,
  relativeHumidity: 0.5,
  intakeMixtureTemperature: 325,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: 1,
  sparkAdvanceDeg: 13,
  dwellTime: 3e-3,
  compressionRatio: 7,
  egrFraction: 0,
  coolantTemperature: 373,
};

/** Trivial simulator: one snapshot per ms, a cycle summary every 100 snapshots. */
class FakeSim implements SimulatorLike {
  n = -1;
  patches: Partial<OperatingPoint>[] = [];
  resets = 0;
  throwAt = Infinity;
  private pending: CycleSummary[] = [];
  get time(): number {
    return Math.max(0, this.n) * 1e-3;
  }
  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.patches.push(patch);
  }
  reset(): void {
    this.n = -1;
    this.resets++;
    this.pending = [];
  }
  advanceToNextSnapshot(): EngineSnapshot {
    this.n++;
    if (this.n >= this.throwAt) throw new Error('boom');
    const s = createEmptySnapshot();
    s.t = this.n * 1e-3;
    s.cycle = Math.floor(this.n / 100);
    if (this.n > 0 && this.n % 100 === 0) this.pending.push({ cycle: s.cycle - 1 } as CycleSummary);
    return s;
  }
  drainCycleSummaries(): CycleSummary[] {
    const p = this.pending;
    this.pending = [];
    return p;
  }
}

function harness(opts: Partial<SimulatorOptions> = {}, clock?: () => number) {
  const posted: FromWorker[] = [];
  const tasks: (() => void)[] = [];
  const port: HostPort = { post: (m) => posted.push(m), schedule: (cb) => tasks.push(cb), now: clock ?? (() => 0) };
  let sim!: FakeSim;
  const host = new SimWorkerHost(port, () => (sim = new FakeSim()));
  const options: SimulatorOptions = { snapshotEveryDeg: 0.5, bufferAheadSeconds: 0.05, ...opts };
  host.handle({ type: 'init', spec: {} as EngineSpec, operatingPoint: OP, options });
  const run = (max = 1000) => {
    let i = 0;
    while (tasks.length > 0 && i++ < max) tasks.shift()!();
  };
  const snapshots = () => posted.flatMap((m) => (m.type === 'snapshots' ? m.batch : []));
  return { host, posted, tasks, run, snapshots, sim: () => sim };
}

describe('SimWorkerHost', () => {
  it('posts ready on init and buffers up to demand + bufferAhead', () => {
    const h = harness();
    expect(h.posted[0]).toEqual({ type: 'ready' });
    h.run();
    const snaps = h.snapshots();
    expect(snaps[0].t).toBe(0);
    expect(snaps[snaps.length - 1].t).toBeCloseTo(0.05, 12);
    expect(h.tasks.length).toBe(0); // idle
    h.host.handle({ type: 'demand', t: 0.1 });
    h.run();
    expect(h.snapshots().at(-1)!.t).toBeGreaterThanOrEqual(0.15);
    expect(h.snapshots().at(-1)!.t).toBeLessThan(0.1515);
    // Snapshots are strictly increasing and contiguous across batches.
    const all = h.snapshots();
    for (let i = 1; i < all.length; i++) expect(all[i].t).toBeGreaterThan(all[i - 1].t);
  });

  it('accepts backwards demand (idles) and caps batch size', () => {
    const h = harness({ bufferAheadSeconds: 2 });
    h.run();
    const batches = h.posted.filter((m) => m.type === 'snapshots');
    expect(batches.length).toBeGreaterThan(1);
    for (const b of batches) if (b.type === 'snapshots') expect(b.batch.length).toBeLessThanOrEqual(HOST_MAX_BATCH);
    const n = h.snapshots().length;
    h.host.handle({ type: 'demand', t: -1.5 });
    h.run();
    expect(h.snapshots().length).toBe(n);
  });

  it('forwards cycle summaries after the snapshots that complete them', () => {
    const h = harness({ bufferAheadSeconds: 0.35 });
    h.run();
    const cycles = h.posted.filter((m) => m.type === 'cycle');
    expect(cycles.map((m) => (m.type === 'cycle' ? m.summary.cycle : -1))).toEqual([0, 1, 2]);
    const firstCycleIdx = h.posted.findIndex((m) => m.type === 'cycle');
    const before = h.posted.slice(0, firstCycleIdx).flatMap((m) => (m.type === 'snapshots' ? m.batch : []));
    expect(before.some((s) => s.cycle === 1)).toBe(true);
  });

  it('time-slices long runs so messages can interleave', () => {
    let now = 0;
    const h = harness({ bufferAheadSeconds: 0.2 }, () => (now += 1)); // every now() call = 1 ms
    h.tasks.shift()!();
    const first = h.posted.filter((m) => m.type === 'snapshots');
    expect(first.length).toBe(1);
    if (first[0].type === 'snapshots') expect(first[0].batch.length).toBeLessThan(200);
    expect(h.tasks.length).toBe(1); // re-scheduled
    h.run();
    expect(h.snapshots().at(-1)!.t).toBeCloseTo(0.2, 12);
  });

  it('forwards operating-point changes and restarts on reset with a fresh ready', () => {
    const h = harness();
    h.run();
    h.host.handle({ type: 'set-operating-point', patch: { throttle: 0.3 } });
    expect(h.sim().patches).toEqual([{ throttle: 0.3 }]);
    h.posted.length = 0;
    h.host.handle({ type: 'demand', t: 1 });
    h.host.handle({ type: 'reset' });
    expect(h.sim().resets).toBe(1);
    expect(h.posted[0]).toEqual({ type: 'ready' });
    h.run();
    const snaps = h.snapshots();
    expect(snaps[0].t).toBe(0);
    expect(snaps.at(-1)!.t).toBeCloseTo(0.05, 12); // demand reset to 0
  });

  it('reports simulator failures as log errors and stops', () => {
    const h = harness({ bufferAheadSeconds: 1 });
    h.sim().throwAt = 10;
    h.run();
    const err = h.posted.find((m) => m.type === 'log');
    expect(err && err.type === 'log' && err.level).toBe('error');
    expect(h.snapshots().length).toBe(10);
    expect(h.tasks.length).toBe(0);
  });

  it('runs the real (mock) simulator through the default factory', () => {
    const posted: FromWorker[] = [];
    const tasks: (() => void)[] = [];
    const host = new SimWorkerHost({ post: (m) => posted.push(m), schedule: (cb) => tasks.push(cb), now: () => 0 });
    host.handle({
      type: 'init',
      spec: CFR_F1,
      operatingPoint: OP,
      options: { snapshotEveryDeg: 1, bufferAheadSeconds: 0.25 },
    });
    while (tasks.length) tasks.shift()!();
    const snaps = posted.flatMap((m) => (m.type === 'snapshots' ? m.batch : []));
    expect(snaps.length).toBeGreaterThan(900);
    expect(snaps.at(-1)!.t).toBeGreaterThanOrEqual(0.25);
    expect(posted.some((m) => m.type === 'cycle')).toBe(true);
    // Messages survive structured cloning (plain data only).
    expect(structuredClone(snaps[500])).toEqual(snaps[500]);
  });
});
