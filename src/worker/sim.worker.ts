/**
 * Simulation Web Worker: hosts a SimulatorLike and streams snapshots to the
 * main thread following src/worker/protocol.ts.
 *
 * Flow control: after `init` (and after every `reset`) the worker posts
 * `{type:'ready'}` and then keeps the simulator running until its time reaches
 * `lastDemand.t + options.bufferAheadSeconds`, posting snapshot batches (plus
 * any completed cycle summaries) in time-sliced chunks so incoming messages
 * (`demand`, `set-operating-point`, `reset`) are handled promptly. A demand may
 * move backwards (the client shrinks its look-ahead in slow motion); the
 * worker then simply idles until playback catches up.
 *
 * Everything posted after a `ready` belongs to the run that `ready` announced
 * (postMessage is FIFO), which is how the client discards stale data after a reset.
 *
 * The logic lives in `SimWorkerHost` (testable without a Worker); the bottom of
 * the file wires it to the dedicated-worker global scope when loaded as a worker.
 */
import type { EngineSnapshot } from '../physics/core/snapshot';
import { EngineSimulator } from '../physics/cycle';
import type { CycleModelOptions } from '../physics/cycle';
import { MockSimulator, type MockSimulatorOptions } from './mock-simulator';
import type { FromWorker, SimulatorOptions, ToWorker } from './protocol';
import type { SimulatorFactory, SimulatorLike } from './simulator-like';

/**
 * Options the worker accepts in the `init` message: the protocol's SimulatorOptions plus
 *  - `simulator`: 'physics' (default, src/physics/cycle EngineSimulator) or 'mock'
 *    (MockSimulator — the lightweight stand-in, e.g. for UI work or slow machines);
 *  - any CycleModelOptions / MockSimulatorOptions overrides for the chosen simulator.
 * (Extra fields survive structured cloning, so the main thread passes them through SimClient
 * unchanged; the app sets `simulator: 'mock'` for the `?sim=mock` URL flag.)
 */
export type WorkerSimulatorOptions = SimulatorOptions &
  Partial<CycleModelOptions> &
  Partial<Omit<MockSimulatorOptions, keyof SimulatorOptions>> & {
    simulator?: 'physics' | 'mock';
  };

/**
 * THE simulator used by the worker: the physics EngineSimulator unless `simulator: 'mock'`. The spec
 * travels in the init message, and spec.id selects the physics model's per-engine option defaults
 * (cycle/options.ts); EngineSimulator handles multi-cylinder specs itself. The mock models a single
 * cylinder only: asking it for a multi-cylinder engine is an error (reported to the app as a worker log
 * error, which the app shows).
 */
export const createSimulator: SimulatorFactory = (spec, op, options) => {
  const o = options as WorkerSimulatorOptions;
  if (o.simulator === 'mock') {
    if (spec.cylinders > 1) throw new Error(mockUnsupportedMessage(spec));
    return new MockSimulator(spec, op, o);
  }
  return new EngineSimulator(spec, op, o);
};

/** Why the mock simulator cannot run a spec (multi-cylinder engines). */
export function mockUnsupportedMessage(spec: { name: string; cylinders: number }): string {
  return `The mock simulator models a single cylinder; ${spec.name} has ${spec.cylinders}. Use the physics simulator (drop ?sim=mock).`;
}

/** Environment the host runs in (the worker global scope, or a test harness). */
export interface HostPort {
  post(msg: FromWorker): void;
  /** Run `cb` as a new macrotask, so queued incoming messages are processed first. */
  schedule(cb: () => void): void;
  /** Monotonic wall clock, ms. */
  now(): number;
}

/** Max wall time spent producing one chunk before yielding, ms. */
export const HOST_SLICE_MS = 8;
/** Max snapshots per posted batch. */
export const HOST_MAX_BATCH = 512;

export class SimWorkerHost {
  private sim: SimulatorLike | null = null;
  private bufferAhead = 0;
  private demandT = 0;
  private emitted = 0;
  private scheduled = false;
  private failed = false;

  constructor(
    private readonly port: HostPort,
    private readonly factory: SimulatorFactory = createSimulator,
  ) {}

  /** Simulated time the host is buffering up to. */
  get targetTime(): number {
    return this.demandT + this.bufferAhead;
  }

  /** Simulated time of the last produced snapshot (NaN before init). */
  get simTime(): number {
    return this.sim ? this.sim.time : NaN;
  }

  handle(msg: ToWorker): void {
    try {
      switch (msg.type) {
        case 'init':
          this.sim = this.factory(msg.spec, msg.operatingPoint, msg.options);
          this.bufferAhead = Math.max(0, msg.options.bufferAheadSeconds);
          this.restart();
          break;
        case 'set-operating-point':
          if (!this.sim) {
            this.log('warn', 'set-operating-point before init ignored');
            return;
          }
          this.sim.setOperatingPoint(msg.patch);
          break;
        case 'demand':
          if (Number.isFinite(msg.t)) this.demandT = msg.t;
          this.kick();
          break;
        case 'reset':
          if (!this.sim) {
            this.log('warn', 'reset before init ignored');
            return;
          }
          this.sim.reset();
          this.restart();
          break;
      }
    } catch (err) {
      this.fail(err);
    }
  }

  /** Produce one chunk of snapshots (called from the scheduler). */
  pump(): void {
    this.scheduled = false;
    const sim = this.sim;
    if (!sim || this.failed) return;
    const target = this.targetTime;
    if (this.emitted > 0 && sim.time >= target) return;
    const start = this.port.now();
    const batch: EngineSnapshot[] = [];
    try {
      while (batch.length < HOST_MAX_BATCH) {
        batch.push(sim.advanceToNextSnapshot());
        if (sim.time >= target) break;
        if ((batch.length & 15) === 0 && this.port.now() - start > HOST_SLICE_MS) break;
      }
    } catch (err) {
      this.flush(batch);
      this.fail(err);
      return;
    }
    this.flush(batch);
    if (sim.time < target) this.kick();
  }

  private flush(batch: EngineSnapshot[]): void {
    if (batch.length > 0) {
      this.emitted += batch.length;
      this.port.post({ type: 'snapshots', batch });
    }
    const sim = this.sim;
    if (!sim) return;
    for (const summary of sim.drainCycleSummaries()) this.port.post({ type: 'cycle', summary });
  }

  private restart(): void {
    this.sim?.drainCycleSummaries();
    this.demandT = 0;
    this.emitted = 0;
    this.failed = false;
    this.port.post({ type: 'ready' });
    this.kick();
  }

  private kick(): void {
    if (this.scheduled || !this.sim || this.failed) return;
    this.scheduled = true;
    this.port.schedule(() => this.pump());
  }

  private fail(err: unknown): void {
    this.failed = true;
    const message = err instanceof Error ? `${err.message}\n${err.stack ?? ''}` : String(err);
    this.log('error', `simulator failed: ${message}`);
  }

  private log(level: 'info' | 'warn' | 'error', message: string): void {
    this.port.post({ type: 'log', level, message });
  }
}

// ---------------------------------------------------------------------------
// Worker bootstrap (only when actually running as a dedicated worker)
// ---------------------------------------------------------------------------

interface WorkerScopeLike {
  postMessage(msg: FromWorker): void;
  onmessage: ((ev: MessageEvent<ToWorker>) => void) | null;
}

function isWorkerScope(): boolean {
  const g = globalThis as { WorkerGlobalScope?: abstract new () => unknown };
  return typeof g.WorkerGlobalScope === 'function' && globalThis instanceof g.WorkerGlobalScope;
}

if (isWorkerScope()) {
  const scope = globalThis as unknown as WorkerScopeLike;
  // MessageChannel yields a macrotask without setTimeout's nested-timer clamping.
  const channel = new MessageChannel();
  let queue: (() => void)[] = [];
  channel.port1.onmessage = () => {
    const q = queue;
    queue = [];
    for (const cb of q) cb();
  };
  const host = new SimWorkerHost({
    post: (msg) => scope.postMessage(msg),
    schedule: (cb) => {
      queue.push(cb);
      if (queue.length === 1) channel.port2.postMessage(null);
    },
    now: () => performance.now(),
  });
  scope.onmessage = (ev) => host.handle(ev.data);
}
