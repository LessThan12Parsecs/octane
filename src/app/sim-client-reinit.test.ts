/**
 * SimClient.reset() after the simulator could not be constructed: the worker is asked to construct it
 * again with the CURRENT operating point (initial point + every patch since), so picking a supported
 * preset and pressing Reset recovers an engine whose start-up point the simulator rejects. A simulator
 * that was running restarts with 'reset' exactly as before.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import type { FromWorker, ToWorker } from '../worker/protocol';
import { SimWorkerHost } from '../worker/sim.worker';
import type { SimulatorLike } from '../worker/simulator-like';
import { createEmptySnapshot, SimClient, type WorkerLike } from './sim-client';

/** Rejects free-speed operation at construction (like a model lacking a load integration). */
class FixedSpeedOnly implements SimulatorLike {
  time = 0;
  constructor(readonly op: OperatingPoint) {
    if (op.speedMode === 'free') throw new Error('free-speed mode is not implemented for this engine');
  }
  setOperatingPoint(): void {}
  reset(): void {
    this.time = 0;
  }
  advanceToNextSnapshot(): EngineSnapshot {
    const s = createEmptySnapshot();
    s.t = this.time;
    s.rpm = this.op.rpm;
    this.time += 1e-3;
    return s;
  }
  drainCycleSummaries(): CycleSummary[] {
    return [];
  }
}

class InProcessWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  readonly received: ToWorker[] = [];
  readonly built: OperatingPoint[] = [];
  private readonly tasks: (() => void)[] = [];
  private readonly host: SimWorkerHost;
  constructor() {
    this.host = new SimWorkerHost(
      {
        post: (m) => this.onmessage?.({ data: structuredClone(m) } as MessageEvent<FromWorker>),
        schedule: (cb) => this.tasks.push(cb),
        now: () => 0,
      },
      (_spec: EngineSpec, op: OperatingPoint) => {
        this.built.push(op);
        return new FixedSpeedOnly(op);
      },
    );
  }
  postMessage(msg: ToWorker): void {
    this.received.push(msg);
    this.host.handle(structuredClone(msg));
  }
  terminate(): void {}
  run(): void {
    for (let i = 0; i < 100 && this.tasks.length; i++) this.tasks.shift()!();
  }
}

const OPTIONS = { snapshotEveryDeg: 1, bufferAheadSeconds: 0.05 };

describe('SimClient.reset after a simulator construction failure', () => {
  it('re-initialises with the operating point patched since, and data flows', () => {
    const w = new InProcessWorker();
    const client = new SimClient(MODEL_T, MODEL_T_CRUISE, OPTIONS, { createWorker: () => w });
    const errors: string[] = [];
    client.onError((m) => errors.push(m));
    expect(client.error).toContain('free-speed mode');
    w.run();
    expect(client.advance(0.1, 1)).toBeNull();

    // Pick a supported setting: rejected again until the simulator exists, but remembered.
    client.setOperatingPoint({ speedMode: 'fixed', rpm: 1600 });
    client.setOperatingPoint({ throttle: 1 });
    client.reset();
    const init = w.received.filter((m) => m.type === 'init');
    expect(init).toHaveLength(2);
    const op = (init[1] as Extract<ToWorker, { type: 'init' }>).operatingPoint;
    expect(op).toEqual({ ...MODEL_T_CRUISE, speedMode: 'fixed', rpm: 1600, throttle: 1 });
    expect(w.received.some((m) => m.type === 'reset')).toBe(false);
    w.run();
    const s = client.advance(0.01, 1);
    expect(s).not.toBeNull();
    expect(s!.rpm).toBe(1600);
    expect(errors).toEqual([]);

    // Running now: a later reset is a plain restart, not another construction.
    client.reset();
    expect(w.received.filter((m) => m.type === 'init')).toHaveLength(2);
    expect(w.received.filter((m) => m.type === 'reset')).toHaveLength(1);
    client.dispose();
  });

  it('keeps re-initialising while the construction still fails', () => {
    const w = new InProcessWorker();
    const client = new SimClient(MODEL_T, MODEL_T_CRUISE, OPTIONS, { createWorker: () => w });
    const errors: string[] = [];
    client.onError((m) => errors.push(m));
    client.reset(); // still free speed
    expect(errors).toHaveLength(1);
    client.setOperatingPoint({ speedMode: 'fixed' });
    client.reset();
    expect(w.built.map((o) => o.speedMode)).toEqual(['free', 'free', 'fixed']);
    w.run();
    expect(client.advance(0.01, 1)).not.toBeNull();
    client.dispose();
  });

  it('a simulator that started restarts with reset, as before (CFR)', () => {
    const w = new InProcessWorker();
    const client = new SimClient(CFR_F1, { ...CFR_RON_CONDITIONS, speedMode: 'fixed' }, OPTIONS, { createWorker: () => w });
    expect(client.error).toBeNull();
    client.setOperatingPoint({ sparkAdvanceDeg: 20 });
    client.reset();
    expect(w.received.map((m) => m.type)).toEqual(['init', 'demand', 'set-operating-point', 'reset', 'demand']);
    client.dispose();
  });
});
