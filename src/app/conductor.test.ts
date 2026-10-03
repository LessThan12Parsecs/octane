import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1, CFR_RON_CONDITIONS, clearanceHeightAtTDC } from '../physics/engines/cfr';
import type { FromWorker, ToWorker } from '../worker/protocol';
import { SimWorkerHost } from '../worker/sim.worker';
import { Conductor, CR_EPSILON, type GasPort, type MechanismPort, type SimPort, type UiPort } from './conductor';
import { createEmptySnapshot, SimClient, type WorkerLike } from './sim-client';

// ---------------------------------------------------------------------------
// Fakes
// ---------------------------------------------------------------------------

function snap(t: number, cr: number, disp = 0): EngineSnapshot {
  const s = createEmptySnapshot();
  s.t = t;
  s.rpm = 600;
  s.thetaDeg = -360 + t * 3600;
  s.pistonDisplacement = disp;
  s.clearanceHeight = clearanceHeightAtTDC(cr, CFR_F1) + disp;
  return s;
}

class FakeSim implements SimPort {
  playbackTime = NaN;
  readonly calls: string[] = [];
  readonly advances: [number, number][] = [];
  readonly patches: Partial<OperatingPoint>[] = [];
  readonly recentSince: number[] = [];
  next: EngineSnapshot | null = null;
  private cycleCb: ((c: CycleSummary) => void) | null = null;
  advance(dtWall: number, timeScale: number): EngineSnapshot | null {
    this.advances.push([dtWall, timeScale]);
    if (this.next) this.playbackTime = this.next.t;
    return this.next;
  }
  stepDegrees(deg: number): EngineSnapshot | null {
    this.calls.push(`step ${deg}`);
    return this.next;
  }
  recentSnapshots(sinceT: number): EngineSnapshot[] {
    this.recentSince.push(sinceT);
    return this.next ? [this.next] : [];
  }
  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.patches.push(patch);
  }
  reset(): void {
    this.calls.push('reset');
    this.playbackTime = NaN;
  }
  onCycle(cb: (c: CycleSummary) => void): void {
    this.cycleCb = cb;
  }
  emitCycle(c: CycleSummary): void {
    this.cycleCb?.(c);
  }
}

class FakeEngine implements MechanismPort {
  compressionRatio: number;
  readonly crCalls: number[] = [];
  readonly updates: number[] = [];
  cutaway = true;
  constructor(cr: number) {
    this.compressionRatio = cr;
  }
  update(s: EngineSnapshot): void {
    this.updates.push(s.t);
  }
  setCompressionRatio(cr: number): void {
    this.crCalls.push(cr);
    this.compressionRatio = cr;
  }
  setCutaway(on: boolean): void {
    this.cutaway = on;
  }
}

class FakeGas implements GasPort {
  readonly updates: [number, number, number][] = [];
  mode = 'physical';
  update(s: EngineSnapshot, dtWall: number, timeScale: number): void {
    this.updates.push([s.t, dtWall, timeScale]);
  }
  setMode(mode: 'physical' | 'temperature'): void {
    this.mode = mode;
  }
}

class FakeUi implements UiPort {
  readonly updates: [number, number][] = [];
  readonly cycles: number[] = [];
  update(s: EngineSnapshot, recent: EngineSnapshot[]): void {
    this.updates.push([s.t, recent.length]);
  }
  pushCycle(c: CycleSummary): void {
    this.cycles.push(c.cycle);
  }
}

function setup(cr = 7) {
  const sim = new FakeSim();
  const engine = new FakeEngine(cr);
  const gas = new FakeGas();
  const ui = new FakeUi();
  const crEvents: [number, number][] = [];
  const views: string[] = [];
  let first = 0;
  const c = new Conductor(CFR_F1, sim, engine, gas, {
    onCompressionRatio: (a, b) => crEvents.push([a, b]),
    onView: (v) => views.push(`${v.cutaway}/${v.mode}`),
    onFirstSnapshot: () => first++,
  });
  c.attachUi(ui);
  return { c, sim, engine, gas, ui, crEvents, views, first: () => first };
}

// ---------------------------------------------------------------------------

describe('Conductor frame loop (fakes)', () => {
  it('does nothing downstream until data arrives', () => {
    const { c, engine, gas, ui, sim } = setup();
    expect(c.frame(0)).toBeNull();
    expect(c.frame(16)).toBeNull();
    expect(sim.advances).toEqual([
      [0, 0.02],
      [0.016, 0.02],
    ]);
    expect(engine.updates).toEqual([]);
    expect(gas.updates).toEqual([]);
    expect(ui.updates).toEqual([]);
  });

  it('drives engine, gas and UI with the same snapshot and the frame dt', () => {
    const { c, sim, engine, gas, ui, first } = setup();
    c.frame(0);
    sim.next = snap(0.001, 7);
    const s = c.frame(20)!;
    expect(s.t).toBe(0.001);
    expect(engine.updates).toEqual([0.001]);
    expect(gas.updates[0][0]).toBe(0.001);
    expect(gas.updates[0][1]).toBeCloseTo(0.02, 12);
    expect(gas.updates[0][2]).toBe(0.02);
    expect(ui.updates).toEqual([[0.001, 1]]);
    expect(first()).toBe(1);
    expect(c.lastFrameSeconds).toBeCloseTo(0.02, 12);
    c.frame(40);
    expect(first()).toBe(1);
    expect(c.frames).toBe(2);
  });

  it('clamps long frames and passes timeScale 0 while paused', () => {
    const { c, sim, gas } = setup();
    sim.next = snap(0.001, 7);
    c.frame(0);
    c.frame(10_000);
    expect(sim.advances[1][0]).toBe(0.1);
    c.handlePlayback({ timeScale: 0.5, paused: true });
    c.frame(10_016);
    expect(sim.advances[2][1]).toBe(0);
    expect(gas.updates[2][2]).toBe(0);
    c.handlePlayback({ timeScale: 0.5, paused: false });
    c.frame(10_032);
    expect(sim.advances[3][1]).toBe(0.5);
    expect(c.effectiveTimeScale).toBe(0.5);
  });

  it('asks for recent snapshots since the previous playback time, and from −∞ after a reset', () => {
    const { c, sim } = setup();
    sim.next = snap(0.001, 7);
    c.frame(0);
    sim.next = snap(0.002, 7);
    c.frame(16);
    expect(sim.recentSince).toEqual([-Infinity, 0.001]);
    c.handleReset();
    expect(sim.calls).toContain('reset');
    expect(c.snapshot).toBeNull();
    sim.next = snap(0.0005, 7);
    c.frame(32);
    expect(sim.recentSince[2]).toBe(-Infinity);
  });

  it('moves the mechanism to the CR found in the data, not the slider', () => {
    const { c, sim, engine, crEvents } = setup(7);
    c.handleOperatingPoint({ compressionRatio: 8 });
    expect(sim.patches).toEqual([{ compressionRatio: 8 }]);
    sim.next = snap(0.001, 7, 0.03);
    c.frame(0);
    expect(engine.crCalls).toEqual([]); // data still at CR 7
    sim.next = snap(0.2, 8, 0.05);
    c.frame(16);
    expect(engine.crCalls).toHaveLength(1);
    expect(engine.crCalls[0]).toBeCloseTo(8, 9);
    expect(crEvents[0][1]).toBe(7);
    // tiny round-off does not re-trigger
    sim.next = snap(0.21, 8 + CR_EPSILON / 10, 0.02);
    c.frame(32);
    expect(engine.crCalls).toHaveLength(1);
  });

  it('clamps a data CR outside the spec range', () => {
    const { c, sim, engine } = setup(7);
    sim.next = snap(0.001, 25);
    c.frame(0);
    expect(engine.compressionRatio).toBe(CFR_F1.geometry.compressionRatioRange[1]);
  });

  it('routes step, view and cycle summaries', () => {
    const { c, sim, engine, gas, ui, views } = setup();
    c.handleStep(-5);
    expect(sim.calls).toEqual(['step -5']);
    c.handleView({ cutaway: false, mode: 'temperature' });
    expect(engine.cutaway).toBe(false);
    expect(gas.mode).toBe('temperature');
    expect(views).toEqual(['false/temperature']);
    expect(c.view).toEqual({ cutaway: false, mode: 'temperature' });
    sim.emitCycle({ cycle: 3 } as CycleSummary);
    expect(ui.cycles).toEqual([3]);
  });
});

// ---------------------------------------------------------------------------
// Integration: real SimClient + worker host + MockSimulator (in process)
// ---------------------------------------------------------------------------

class InProcessWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  private readonly toHost: ToWorker[] = [];
  private readonly toClient: FromWorker[] = [];
  private readonly tasks: (() => void)[] = [];
  private readonly host = new SimWorkerHost({
    post: (m) => this.toClient.push(structuredClone(m)),
    schedule: (cb) => this.tasks.push(cb),
    now: () => 0,
  });
  postMessage(msg: ToWorker): void {
    this.toHost.push(structuredClone(msg));
  }
  terminate(): void {}
  flush(): void {
    for (let i = 0; i < 50 && (this.toHost.length || this.tasks.length || this.toClient.length); i++) {
      for (let k = 0; k < 10000 && (this.toHost.length || this.tasks.length); k++) {
        if (this.toHost.length) this.host.handle(this.toHost.shift()!);
        else this.tasks.shift()!();
      }
      while (this.toClient.length) this.onmessage?.({ data: this.toClient.shift()! } as MessageEvent<FromWorker>);
    }
  }
}

describe('Conductor + SimClient + MockSimulator', () => {
  it('plays cycles, forwards summaries and follows a CR change through the data', () => {
    const w = new InProcessWorker();
    const op: OperatingPoint = { ...CFR_RON_CONDITIONS, compressionRatio: 6.5 };
    const sim = new SimClient(CFR_F1, op, { snapshotEveryDeg: 1, bufferAheadSeconds: 0.25 }, { createWorker: () => w });
    const engine = new FakeEngine(op.compressionRatio);
    const gas = new FakeGas();
    const ui = new FakeUi();
    const c = new Conductor(CFR_F1, sim, engine, gas, {}, { timeScale: 1, paused: false });
    c.attachUi(ui);

    // ~0.5 s of simulated time (2.5 cycles at 600 rpm) in 25 ms frames at real time
    let now = 0;
    const run = (frames: number) => {
      for (let i = 0; i < frames; i++) {
        w.flush();
        c.frame(now);
        now += 25;
      }
    };
    run(22);
    expect(sim.playbackTime).toBeGreaterThan(0.4);
    expect(ui.updates.length).toBeGreaterThan(15);
    expect(ui.cycles.length).toBeGreaterThanOrEqual(1);
    // the data CR equals the operating point: no mechanism move
    expect(engine.crCalls).toEqual([]);
    const ingested = ui.updates.reduce((n, u) => n + u[1], 0);
    expect(ingested).toBeGreaterThan(300); // raw 1° samples delivered to the charts

    c.handleOperatingPoint({ compressionRatio: 8 });
    run(30);
    expect(engine.crCalls.length).toBeGreaterThanOrEqual(1);
    expect(engine.compressionRatio).toBeCloseTo(8, 6);
    sim.dispose();
  });
});
