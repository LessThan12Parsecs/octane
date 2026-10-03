import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import type { FromWorker, SimulatorOptions, ToWorker } from '../worker/protocol';
import { SimWorkerHost } from '../worker/sim.worker';
import {
  SimClient,
  SnapshotBuffer,
  type WorkerLike,
  copySnapshot,
  createEmptySnapshot,
  cumulativeCrankDeg,
  interpolateSnapshot,
  lerpCrankAngleDeg,
} from './sim-client';

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

function snap(t: number, cycle: number, thetaDeg: number, extra: (s: EngineSnapshot) => void = () => {}): EngineSnapshot {
  const s = createEmptySnapshot();
  s.t = t;
  s.cycle = cycle;
  s.thetaDeg = thetaDeg;
  s.rpm = 600;
  extra(s);
  return s;
}

describe('crank-angle interpolation', () => {
  it('lerps normally inside a cycle', () => {
    expect(lerpCrankAngleDeg(-10, 10, 0.25)).toBeCloseTo(-5, 12);
  });
  it('takes the short way across the ±360 wrap', () => {
    expect(lerpCrankAngleDeg(359.5, -360, 0.5)).toBeCloseTo(359.75, 12);
    expect(lerpCrankAngleDeg(359.5, -359.5, 0.5)).toBe(-360);
    expect(lerpCrankAngleDeg(359.5, -359.5, 0.75)).toBeCloseTo(-359.75, 12);
    expect(lerpCrankAngleDeg(359.5, -360, 1)).toBe(-360);
  });
  it('cumulative crank angle is continuous across cycles', () => {
    expect(cumulativeCrankDeg({ cycle: 0, thetaDeg: -360 })).toBe(0);
    expect(cumulativeCrankDeg({ cycle: 0, thetaDeg: 359.5 })).toBe(719.5);
    expect(cumulativeCrankDeg({ cycle: 1, thetaDeg: -360 })).toBe(720);
  });
});

describe('interpolateSnapshot', () => {
  const a = snap(1, 3, 10, (s) => {
    s.pressure = 10e5;
    s.phase = 'combustion';
    s.flame.stage = 'kernel';
    s.flame.center = [0, -0.002, 0.03];
    s.spark.phase = 'arc';
    s.knock.autoignited = false;
    s.burnedComposition.CO2 = 0.1;
  });
  const b = snap(2, 3, 12, (s) => {
    s.pressure = 20e5;
    s.phase = 'expansion';
    s.flame.stage = 'turbulent';
    s.flame.center = [0.002, -0.002, 0.03];
    s.spark.phase = 'glow';
    s.knock.autoignited = true;
    s.burnedComposition.CO2 = 0.2;
  });
  it('linear for scalars, discrete fields from the earlier sample', () => {
    const out = interpolateSnapshot(a, b, 0.25, createEmptySnapshot());
    expect(out.t).toBeCloseTo(1.25, 12);
    expect(out.thetaDeg).toBeCloseTo(10.5, 12);
    expect(out.pressure).toBeCloseTo(12.5e5, 6);
    expect(out.flame.center[0]).toBeCloseTo(0.0005, 12);
    expect(out.burnedComposition.CO2).toBeCloseTo(0.125, 12);
    expect(out.phase).toBe('combustion');
    expect(out.flame.stage).toBe('kernel');
    expect(out.spark.phase).toBe('arc');
    expect(out.knock.autoignited).toBe(false);
    expect(out.cycle).toBe(3);
  });
  it('does not alias the inputs', () => {
    const out = copySnapshot(a, createEmptySnapshot());
    expect(out).toEqual(a);
    expect(out.flame).not.toBe(a.flame);
    expect(out.flame.center).not.toBe(a.flame.center);
  });
  it('takes cycle and discrete fields from the later sample once past the wrap', () => {
    const e = snap(0, 4, 359.5, (s) => (s.flame.stage = 'done'));
    const f = snap(1, 5, -359.5, (s) => (s.flame.stage = 'none'));
    const before = interpolateSnapshot(e, f, 0.25, createEmptySnapshot());
    expect(before.thetaDeg).toBeCloseTo(359.75, 12);
    expect(before.cycle).toBe(4);
    expect(before.flame.stage).toBe('done');
    const after = interpolateSnapshot(e, f, 0.75, createEmptySnapshot());
    expect(after.thetaDeg).toBeCloseTo(-359.75, 12);
    expect(after.cycle).toBe(5);
    expect(after.flame.stage).toBe('none');
  });
});

describe('SnapshotBuffer', () => {
  const series = (n: number, t0 = 0) => Array.from({ length: n }, (_, i) => snap(t0 + i * 0.001, 0, -360 + i));

  it('keeps order, grows past its capacity, drops non-increasing samples', () => {
    const buf = new SnapshotBuffer(16);
    expect(buf.pushBatch(series(100))).toBe(100);
    expect(buf.capacity).toBeGreaterThanOrEqual(100);
    expect(buf.length).toBe(100);
    for (let i = 0; i < 100; i++) expect(buf.at(i).t).toBeCloseTo(i * 0.001, 12);
    expect(buf.push(snap(0.05, 0, 0))).toBe(false);
    expect(buf.push(snap(0.099, 0, 0))).toBe(false);
    expect(buf.length).toBe(100);
  });

  it('wraps around when trimmed and refilled', () => {
    const buf = new SnapshotBuffer(16);
    let t = 0;
    for (let round = 0; round < 50; round++) {
      for (let i = 0; i < 8; i++) buf.push(snap((t += 0.001), 0, 0));
      buf.trimBefore(t - 0.0055);
    }
    expect(buf.capacity).toBe(16);
    expect(buf.length).toBeLessThanOrEqual(7);
    for (let i = 1; i < buf.length; i++) expect(buf.at(i).t).toBeGreaterThan(buf.at(i - 1).t);
    expect(buf.latestTime).toBeCloseTo(t, 12);
  });

  it('finds, samples (clamped) and collects by time', () => {
    const buf = new SnapshotBuffer();
    buf.pushBatch(series(10));
    expect(buf.indexAtOrBefore(-1)).toBe(-1);
    expect(buf.indexAtOrBefore(0)).toBe(0);
    expect(buf.indexAtOrBefore(0.0035)).toBe(3);
    expect(buf.indexAtOrBefore(0.004)).toBe(4);
    expect(buf.indexAtOrBefore(1)).toBe(9);
    const out = createEmptySnapshot();
    expect(buf.sample(0.0035, out)!.thetaDeg).toBeCloseTo(-356.5, 9);
    expect(out.t).toBe(0.0035);
    expect(buf.sample(-5, out)!.t).toBe(0);
    expect(buf.sample(5, out)!.t).toBeCloseTo(0.009, 12);
    const got = buf.collect(0.002, 0.005, []);
    expect(got.map((s) => s.thetaDeg)).toEqual([-358, -357, -356, -355]);
    expect(buf.collect(0.0021, 0.0049, []).length).toBe(2);
    expect(new SnapshotBuffer().sample(0, out)).toBeNull();
  });

  it('trimBefore keeps the interpolation bracket', () => {
    const buf = new SnapshotBuffer();
    buf.pushBatch(series(10));
    expect(buf.trimBefore(0.0055)).toBe(5);
    expect(buf.earliestTime).toBeCloseTo(0.005, 12);
    expect(buf.trimBefore(0.0055)).toBe(0);
  });

  it('maps cumulative crank angle to time across a cycle wrap', () => {
    const buf = new SnapshotBuffer();
    buf.push(snap(0.1, 0, 358));
    buf.push(snap(0.2, 0, 359));
    buf.push(snap(0.3, 1, -360));
    buf.push(snap(0.4, 1, -359));
    const phi = cumulativeCrankDeg({ cycle: 0, thetaDeg: 359.5 });
    expect(buf.timeAtCumulativeAngle(phi)).toBeCloseTo(0.25, 12);
    expect(buf.timeAtCumulativeAngle(phi + 1)).toBeCloseTo(0.35, 12);
    expect(buf.timeAtCumulativeAngle(-1e9)).toBe(0.1);
    expect(buf.timeAtCumulativeAngle(1e9)).toBe(0.4);
  });
});

/** An in-process "worker": SimWorkerHost + MockSimulator driven synchronously, with structured cloning. */
class InProcessWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  readonly toHost: ToWorker[] = [];
  readonly toClient: FromWorker[] = [];
  readonly sent: ToWorker[] = [];
  private readonly tasks: (() => void)[] = [];
  terminated = false;
  readonly host = new SimWorkerHost({
    post: (m) => this.toClient.push(structuredClone(m)),
    schedule: (cb) => this.tasks.push(cb),
    now: () => 0,
  });
  postMessage(msg: ToWorker): void {
    const m = structuredClone(msg);
    this.sent.push(m);
    this.toHost.push(m);
  }
  terminate(): void {
    this.terminated = true;
  }
  /** Worker side: process inbound messages and run scheduled chunks (without delivering). */
  runWorker(): void {
    for (let i = 0; i < 10000 && (this.toHost.length || this.tasks.length); i++) {
      if (this.toHost.length) this.host.handle(this.toHost.shift()!);
      else this.tasks.shift()!();
    }
  }
  /** Main side: deliver everything the worker posted. */
  deliver(): void {
    while (this.toClient.length) this.onmessage?.({ data: this.toClient.shift()! } as MessageEvent<FromWorker>);
  }
  flush(): void {
    for (let i = 0; i < 20 && (this.toHost.length || this.tasks.length || this.toClient.length); i++) {
      this.runWorker();
      this.deliver();
    }
  }
}

function makeClient(options: Partial<SimulatorOptions> = {}) {
  const w = new InProcessWorker();
  const client = new SimClient(
    CFR_F1,
    OP,
    { snapshotEveryDeg: 0.5, bufferAheadSeconds: 0.5, ...options },
    { createWorker: () => w },
  );
  return { w, client };
}

describe('SimClient (worker-less)', () => {
  it('initialises, resolves ready and returns null until data arrives', async () => {
    const { w, client } = makeClient();
    expect(w.sent[0].type).toBe('init');
    expect(w.sent[1].type).toBe('demand');
    expect(client.advance(0.016, 1)).toBeNull();
    expect(Number.isNaN(client.playbackTime)).toBe(true);
    w.flush();
    await client.ready;
    const s = client.advance(0, 1)!;
    expect(s.t).toBe(0);
    expect(s.thetaDeg).toBe(-360);
    expect(client.playbackTime).toBe(0);
  });

  it('plays back at dtWall × timeScale with interpolation', () => {
    const { w, client } = makeClient();
    w.flush();
    client.advance(0, 0.01);
    const s = client.advance(0.1, 0.01)!; // 1 ms simulated = 3.6° at 600 rpm
    expect(client.playbackTime).toBeCloseTo(1e-3, 12);
    expect(s.t).toBeCloseTo(1e-3, 12);
    expect(s.thetaDeg).toBeCloseTo(-360 + 3.6, 6);
    const s2 = client.advance(0.05, 0.01)!;
    expect(s2.thetaDeg).toBeCloseTo(-360 + 5.4, 6);
    expect(s2).not.toBe(s); // alternating scratch objects
    expect(s.thetaDeg).toBeCloseTo(-360 + 3.6, 6); // previous result still intact
  });

  it('stalls at the newest buffered sample instead of running ahead', () => {
    const { w, client } = makeClient();
    w.flush();
    const until = client.bufferedUntil;
    expect(until).toBeGreaterThan(0.01);
    const s = client.advance(10, 1)!;
    expect(client.playbackTime).toBe(until);
    expect(s.t).toBe(until);
    // Worker is asked for more; after it responds playback can continue.
    w.flush();
    expect(client.bufferedUntil).toBeGreaterThan(until);
    client.advance(0.001, 1);
    expect(client.playbackTime).toBeGreaterThan(until);
  });

  it('keeps the worker look-ahead proportional to the time scale', () => {
    const { w, client } = makeClient({ bufferAheadSeconds: 1 });
    w.flush();
    for (let k = 0; k < 20; k++) {
      client.advance(1, 0.01); // 10 ms simulated per call
      w.flush();
    }
    // ahead = max(minAhead 0.02 s, 0.01 × 0.5 s) → ~0.02 s simulated (not the full 1 s)
    expect(client.playbackTime).toBeGreaterThan(0.1);
    expect(client.bufferedUntil - client.playbackTime).toBeLessThan(0.03);
    client.advance(0.001, 1);
    w.flush();
    expect(client.bufferedUntil - client.playbackTime).toBeGreaterThan(0.4);
  });

  it('steps by exact crank degrees while paused, forwards and backwards', () => {
    const { w, client } = makeClient();
    w.flush();
    const s0 = client.advance(0, 0)!;
    const phi0 = cumulativeCrankDeg(s0);
    const s1 = client.stepDegrees(1)!;
    expect(cumulativeCrankDeg(s1) - phi0).toBeCloseTo(1, 6);
    const s2 = client.stepDegrees(0.25)!;
    expect(cumulativeCrankDeg(s2) - phi0).toBeCloseTo(1.25, 6);
    const s3 = client.stepDegrees(-1)!;
    expect(cumulativeCrankDeg(s3) - phi0).toBeCloseTo(0.25, 6);
    // Stepping across the cycle wrap.
    for (let i = 0; i < 40; i++) {
      client.stepDegrees(18);
      w.flush();
    }
    const s4 = client.stepDegrees(0)!;
    expect(cumulativeCrankDeg(s4) - phi0).toBeCloseTo(720.25, 6);
    expect(s4.cycle).toBe(1);
    expect(s4.thetaDeg).toBeCloseTo(-359.75, 6);
  });

  it('delivers each cycle summary once, after playback passes the end of that cycle', () => {
    const { w, client } = makeClient();
    const got: CycleSummary[] = [];
    client.onCycle((c) => {
      got.push(c);
      // Playback is already in a later cycle when the summary arrives.
      const last = client.recentSnapshots(client.playbackTime - 0.01).at(-1)!;
      expect(last.cycle).toBeGreaterThan(c.cycle);
    });
    w.flush();
    // One cycle = 0.2 s at 600 rpm.
    for (let k = 0; k < 390; k++) {
      client.advance(0.01, 0.2); // 2 ms per call
      w.flush();
      if (k === 90) expect(got.length).toBe(0); // 0.18 s: cycle 0 not finished
    }
    expect(client.playbackTime).toBeCloseTo(0.78, 6);
    expect(got.map((c) => c.cycle)).toEqual([0, 1, 2]);
    expect(got[0].imepNet).toBeGreaterThan(3e5);
  });

  it('bounds memory: history is trimmed to ~historyCycles', () => {
    const { w, client } = makeClient();
    w.flush();
    for (let k = 0; k < 600; k++) {
      client.advance(0.01, 0.5);
      w.flush();
    }
    // 4 cycles of history at 0.5° + look-ahead + event refinements.
    expect(client.bufferedCount).toBeLessThan(4 * 1440 * 1.3 + 0.5 * 7200);
    const recent = client.recentSnapshots(client.playbackTime - 0.4);
    expect(recent.length).toBeGreaterThan(2 * 1440);
    expect(recent[0].t).toBeGreaterThanOrEqual(client.playbackTime - 0.4);
    expect(recent.at(-1)!.t).toBeLessThanOrEqual(client.playbackTime);
  });

  it('reset discards stale in-flight data and restarts at t = 0', () => {
    const { w, client } = makeClient();
    w.flush();
    for (let k = 0; k < 50; k++) {
      client.advance(0.01, 0.2);
      w.flush();
    }
    expect(client.playbackTime).toBeGreaterThan(0.05);
    // Worker produced more data that is still "in flight" when the client resets.
    client.advance(0.01, 0.2);
    w.runWorker();
    const inFlight = w.toClient.length;
    expect(inFlight).toBeGreaterThan(0);
    client.reset();
    expect(Number.isNaN(client.playbackTime)).toBe(true);
    expect(client.advance(0.01, 1)).toBeNull();
    w.deliver(); // stale messages arrive first — must be ignored
    expect(client.bufferedCount).toBe(0);
    w.flush(); // then reset → ready → fresh data
    const s = client.advance(0, 1)!;
    expect(s.t).toBe(0);
    expect(s.cycle).toBe(0);
  });

  it('forwards operating-point changes and terminates on dispose', () => {
    const { w, client } = makeClient();
    client.setOperatingPoint({ throttle: 0.4 });
    expect(w.sent.at(-1)).toEqual({ type: 'set-operating-point', patch: { throttle: 0.4 } });
    client.dispose();
    expect(w.terminated).toBe(true);
    expect(client.advance(0.1, 1)).toBeNull();
  });
});
