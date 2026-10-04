/**
 * DEV ONLY — the REAL front end (App: Stage, the engine's render model, one CombustionVisuals per
 * cylinder, SimClient buffering + interpolation, Conductor, UI, engine picker) driven by a SYNTHETIC
 * snapshot stream instead of the physics, for visual checks of the 3D model, the in-cylinder visuals and
 * the multi-cylinder HUD/charts while the physics of an engine is still being integrated.
 *
 * The "worker" is the production SimWorkerHost (src/worker/sim.worker.ts) running on the main thread
 * behind a WorkerLike with structured-cloned, asynchronous messages, so the SimClient flow control is the
 * real one; its simulator is the UI harness's cartoon stream (src/ui/dev): MultiMockStream for
 * multi-cylinder specs (phase-shifted single-cylinder cartoons with a trembler-shower cartoon and
 * cumulative breakdownCount), MockStream for single-cylinder ones. NOT a physics model — never import it
 * from production code.
 *
 * Open http://localhost:<port>/src/app/dev/mock-app.html?engine=ford-model-t (all the app's URL knobs work:
 * &view=engine|chamber, &cyl=n, &ts=, &paused=1 …). Console: `octane` is the App.
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../../physics/core/snapshot';
import { MockStream } from '../../ui/dev/mock-stream';
import { MultiMockStream } from '../../ui/dev/multi-mock-stream';
import type { FromWorker, ToWorker } from '../../worker/protocol';
import { SimWorkerHost } from '../../worker/sim.worker';
import type { SimulatorLike } from '../../worker/simulator-like';
import { startApp, type App } from '../index';
import type { WorkerLike } from '../sim-client';

function fillFrontArea(f: EngineSnapshot['flame']): void {
  if (!(f.area > 0) && f.radius > 0) f.area = 2 * Math.PI * f.radius * f.radius;
}

/** SimulatorLike over the cartoon streams (cadence as in src/ui/dev/ui-harness.ts). */
class CartoonSimulator implements SimulatorLike {
  private readonly stream: MockStream | MultiMockStream;
  private rpm: number;
  private sparkDeg: number;

  constructor(spec: EngineSpec, op: OperatingPoint) {
    this.stream = spec.cylinders > 1 ? new MultiMockStream(spec, op) : new MockStream(spec, op);
    this.rpm = op.rpm;
    this.sparkDeg = -op.sparkAdvanceDeg;
  }

  get time(): number {
    return this.stream.t;
  }

  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.stream.setOperatingPoint(patch);
    if (patch.rpm !== undefined) this.rpm = patch.rpm;
    if (patch.sparkAdvanceDeg !== undefined) this.sparkDeg = -patch.sparkAdvanceDeg;
  }

  reset(): void {
    this.stream.reset();
  }

  advanceToNextSnapshot(): EngineSnapshot {
    const s = this.stream.step(this.cadenceDeg());
    // The cartoon reports no flame-front area; the renderer divides the heat release by it (areal
    // chemiluminescence), so give it a cartoon one: a hemisphere 2πr² (plug in the roof). Shape only.
    fillFrontArea(s.flame);
    if (s.cylinders) for (const c of s.cylinders) fillFrontArea(c.flame);
    return s;
  }

  drainCycleSummaries(): CycleSummary[] {
    return this.stream.drainCycles();
  }

  private cadenceDeg(): number {
    const s = this.stream;
    if (s instanceof MultiMockStream) return s.inSparkWindow() ? 20e-6 * 6 * this.rpm : 0.5; // 20 µs during timer contacts
    const th = s.theta;
    if (th >= this.sparkDeg - 1 && th < this.sparkDeg + 12) return 4e-6 * 6 * this.rpm; // 4 µs around the spark
    if (th >= -20 && th < 60) return 0.05; // knock ringing
    return 0.5;
  }
}

/** A Worker stand-in hosting SimWorkerHost on the main thread (asynchronous, structured-cloned messages). */
class InPageWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  private terminated = false;
  private readonly host = new SimWorkerHost(
    {
      post: (msg) => {
        const data = structuredClone(msg);
        setTimeout(() => {
          if (!this.terminated) this.onmessage?.({ data } as MessageEvent<FromWorker>);
        }, 0);
      },
      schedule: (cb) => {
        setTimeout(() => {
          if (!this.terminated) cb();
        }, 0);
      },
      now: () => performance.now(),
    },
    (spec, op) => new CartoonSimulator(spec, op),
  );

  postMessage(msg: ToWorker): void {
    const data = structuredClone(msg);
    setTimeout(() => {
      if (!this.terminated) this.host.handle(data);
    }, 0);
  }

  terminate(): void {
    this.terminated = true;
  }
}

const viewport = document.getElementById('viewport');
const uiContainer = document.getElementById('ui');
if (!viewport || !uiContainer) throw new Error('mock-app.html must provide #viewport and #ui');

const app: App | null = startApp({
  viewport,
  uiContainer,
  search: location.search,
  simClientConfig: { createWorker: () => new InPageWorker() },
});
(globalThis as { octane?: App }).octane = app ?? undefined;

if (import.meta.hot) import.meta.hot.dispose(() => app?.dispose());
