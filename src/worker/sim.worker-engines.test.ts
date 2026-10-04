/** Worker host with registry engines: mock-simulator guard and error reporting to the client. */
import { describe, expect, it } from 'vitest';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import { SimClient, type WorkerLike } from '../app/sim-client';
import type { FromWorker, ToWorker } from './protocol';
import { createSimulator, mockUnsupportedMessage, SimWorkerHost } from './sim.worker';

class InProcessWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  readonly posted: FromWorker[] = [];
  private readonly tasks: (() => void)[] = [];
  private readonly host = new SimWorkerHost({
    post: (m) => {
      this.posted.push(m);
      this.onmessage?.({ data: structuredClone(m) } as MessageEvent<FromWorker>);
    },
    schedule: (cb) => this.tasks.push(cb),
    now: () => 0,
  });
  postMessage(msg: ToWorker): void {
    this.host.handle(structuredClone(msg));
  }
  terminate(): void {}
  run(): void {
    for (let i = 0; i < 100 && this.tasks.length; i++) this.tasks.shift()!();
  }
}

describe('createSimulator', () => {
  it('refuses a multi-cylinder engine on the single-cylinder mock with a clear message', () => {
    const opts = { snapshotEveryDeg: 1, bufferAheadSeconds: 0.1, simulator: 'mock' as const };
    expect(() => createSimulator(MODEL_T, MODEL_T_CRUISE, opts)).toThrow(mockUnsupportedMessage(MODEL_T));
    expect(mockUnsupportedMessage(MODEL_T)).toContain('4');
    expect(mockUnsupportedMessage(MODEL_T)).toContain('?sim=mock');
    // the CFR still runs on the mock
    const sim = createSimulator(CFR_F1, CFR_RON_CONDITIONS, opts);
    expect(sim.advanceToNextSnapshot().t).toBe(0);
  });

  it('the host turns the refusal into a log error and never announces ready', () => {
    const posted: FromWorker[] = [];
    const host = new SimWorkerHost({ post: (m) => posted.push(m), schedule: () => {}, now: () => 0 });
    host.handle({ type: 'init', spec: MODEL_T, operatingPoint: MODEL_T_CRUISE, options: { snapshotEveryDeg: 1, bufferAheadSeconds: 0.1, simulator: 'mock' } as never });
    expect(posted.some((m) => m.type === 'ready')).toBe(false);
    const err = posted.find((m) => m.type === 'log');
    expect(err && err.type === 'log' && err.level).toBe('error');
    expect(err && err.type === 'log' && err.message).toContain('mock simulator models a single cylinder');
  });

  it('SimClient surfaces the failure: onError fires and ready rejects', async () => {
    const w = new InProcessWorker();
    const errors: string[] = [];
    const client = new SimClient(
      MODEL_T,
      MODEL_T_CRUISE,
      { snapshotEveryDeg: 1, bufferAheadSeconds: 0.1, simulator: 'mock' } as never,
      { createWorker: () => w },
    );
    // The init already ran synchronously inside the constructor: attach afterwards and check the record.
    client.onError((m) => errors.push(m));
    expect(client.error).toContain('mock simulator models a single cylinder');
    await expect(client.ready).rejects.toThrow('single cylinder');
    w.run();
    expect(client.advance(0.1, 1)).toBeNull();
    expect(errors).toEqual([]); // nothing further after the first failure
    client.dispose();
  });
});
