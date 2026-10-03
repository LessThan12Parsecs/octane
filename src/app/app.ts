/**
 * Octane front end: wires the simulation client, the 3D stage (engine
 * mechanism + in-cylinder combustion visuals) and the UI together and runs the
 * animation loop.
 *
 *   SimClient (worker) ──snapshots──▶ Conductor ──▶ EngineModel.update
 *                                               ├─▶ CombustionVisuals.update
 *                                               └─▶ UIController.update (+ pushCycle)
 *   UIController callbacks ──▶ Conductor (operating point, playback, step, view, reset)
 *
 * Starts on the CFR Research-method (RON) conditions in slow motion, framed on
 * the combustion chamber.
 */
import * as THREE from 'three';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1, CFR_RON_CONDITIONS } from '../physics/engines/cfr';
import { CombustionVisuals, createTemperatureLegendElement, EngineModel, recommendedCameraView, Stage } from '../render/index';
import { UIController } from '../ui/index';
import type { SimulatorOptions } from '../worker/protocol';
import { Conductor, type ViewState } from './conductor';
import { ViewportOverlay } from './overlay';
import { SimClient } from './sim-client';
import { clampOperatingPoint, framingDistanceScale, parseUrlOptions, type Framing } from './sync';

/** Default playback speed: 1/50 of real time (a 600 rpm burn of ~5 ms takes ~¼ s). */
export const DEFAULT_TIME_SCALE = 1 / 50;

/**
 * Snapshot cadence and look-ahead of the simulation worker. 0.5° gives smooth
 * kinematics and charts; the simulator adds its own samples around the spark
 * (µs) and during knock ringing. The look-ahead cap bounds how long an
 * operating-point change can wait behind data buffered before a slow-down
 * (0.15 s simulated → 7.5 s of wall time at 1/50×; SimClient keeps the
 * look-ahead far smaller than this while already in slow motion).
 */
export const DEFAULT_SIMULATOR_OPTIONS: SimulatorOptions = { snapshotEveryDeg: 0.5, bufferAheadSeconds: 0.15 };

/** Seconds without a first snapshot before the status card reports a stall. */
const STALL_WARNING_SECONDS = 10;

export interface AppOptions {
  /** Element with class `oct-viewport` that receives the canvas. */
  viewport: HTMLElement;
  /** Container for the UI overlay (controls, charts, HUD). */
  uiContainer: HTMLElement;
  spec?: EngineSpec;
  /** Initial operating point (default: CFR RON conditions). */
  operatingPoint?: OperatingPoint;
  simulatorOptions?: SimulatorOptions;
  /** `location.search` for the URL knobs (?cr=&on=&rpm=&spark=&phi=&ts=&paused=&view=). */
  search?: string;
}

export class App {
  readonly spec: EngineSpec;
  readonly stage: Stage;
  readonly engine: EngineModel;
  readonly combustion: CombustionVisuals;
  readonly sim: SimClient;
  readonly conductor: Conductor;
  readonly ui: UIController;
  private readonly overlay: ViewportOverlay;
  private framing: Framing;
  private raf = 0;
  private startedAt = NaN;
  private hasData = false;
  private failed = false;
  private disposed = false;

  constructor(opts: AppOptions) {
    const spec = (this.spec = opts.spec ?? CFR_F1);
    const url = parseUrlOptions(opts.search ?? '');
    const op = clampOperatingPoint(spec, { ...(opts.operatingPoint ?? CFR_RON_CONDITIONS), ...url.op });
    this.framing = url.framing ?? 'chamber';

    // If anything below throws (e.g. no WebGL 2), release what was already created.
    const created: { dispose(): void }[] = [];
    const track = <T extends { dispose(): void }>(x: T): T => {
      created.push(x);
      return x;
    };
    try {
      // ---- 3D ----
      this.stage = track(new Stage(opts.viewport));
      this.engine = track(new EngineModel(spec, op.compressionRatio));
      this.stage.scene.add(this.engine.root);
      this.combustion = track(new CombustionVisuals(spec));
      this.engine.cylinderFrame.add(this.combustion.root);
      this.stage.addEmitters(this.combustion.root); // flame, spark, tracers: the only bloom sources
      this.engine.root.updateMatrixWorld(true);
      this.stage.fitToBounds(new THREE.Box3().setFromObject(this.engine.root));
      this.stage.setFraming(this.framingView(), false);

      const legend = createTemperatureLegendElement();
      this.overlay = track(
        new ViewportOverlay(opts.viewport, {
          onFraming: (f) => this.setFraming(f),
          legend: legend.element,
        }),
      );
      this.overlay.setFraming(this.framing);
      this.overlay.setStatus('Starting the simulator…', 'Spawning the physics worker and running the first cycle.');

      // ---- simulation ----
      this.sim = track(new SimClient(spec, op, opts.simulatorOptions ?? DEFAULT_SIMULATOR_OPTIONS));
      this.conductor = new Conductor(
        spec,
        this.sim,
        this.engine,
        this.combustion,
        {
          onCompressionRatio: (cr, prev) => this.followCylinder(cr, prev),
          onView: (v) => this.applyView(v),
          onFirstSnapshot: () => this.onData(),
        },
        { timeScale: url.timeScale ?? DEFAULT_TIME_SCALE, paused: url.paused ?? false },
      );

      // ---- UI ----
      this.ui = track(
        new UIController(opts.uiContainer, {
          spec,
          initialOperatingPoint: op,
          onOperatingPointChange: this.conductor.handleOperatingPoint,
          onPlaybackChange: this.conductor.handlePlayback,
          onStep: this.conductor.handleStep,
          onViewChange: this.conductor.handleView,
          onReset: () => {
            this.conductor.handleReset();
            this.hasData = false;
            this.startedAt = NaN;
          },
        }),
      );
      this.conductor.attachUi(this.ui);
    } catch (err) {
      for (let i = created.length - 1; i >= 0; i--) {
        try {
          created[i].dispose();
        } catch {
          /* best effort */
        }
      }
      throw err;
    }
    this.sim.ready.catch((err: unknown) => this.fail('The simulation worker failed to start.', err));
    // The UI announces its initial playback state in a microtask; apply URL overrides after it.
    if (url.timeScale !== undefined || url.paused !== undefined) {
      queueMicrotask(() => {
        if (!this.disposed) this.ui.setPlayback({ timeScale: url.timeScale, paused: url.paused });
      });
    }

    window.addEventListener('keydown', this.onKey);
  }

  /** Start the animation loop. */
  start(): void {
    if (this.raf || this.disposed) return;
    this.raf = requestAnimationFrame(this.tick);
  }

  /** Stop the animation loop (the worker keeps its buffer). */
  stop(): void {
    if (this.raf) cancelAnimationFrame(this.raf);
    this.raf = 0;
  }

  get currentFraming(): Framing {
    return this.framing;
  }

  /** Frame the combustion chamber (close-up through the cutaway) or the whole engine. */
  setFraming(f: Framing, animate = true): void {
    this.framing = f;
    this.overlay.setFraming(f);
    this.stage.setFraming(this.framingView(), animate);
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.stop();
    window.removeEventListener('keydown', this.onKey);
    this.ui.dispose();
    this.sim.dispose();
    this.combustion.dispose();
    this.engine.dispose();
    this.overlay.dispose();
    this.stage.dispose();
  }

  // ------------------------------------------------------------------

  private readonly tick = (now: number): void => {
    if (this.disposed) return;
    this.raf = requestAnimationFrame(this.tick);
    if (Number.isNaN(this.startedAt)) this.startedAt = now;
    let s: EngineSnapshot | null = null;
    try {
      s = this.conductor.frame(now);
    } catch (err) {
      this.fail('The render loop stopped on an error.', err);
      this.stop();
      return;
    }
    if (!s && !this.hasData && !this.failed && (now - this.startedAt) / 1000 > STALL_WARNING_SECONDS) {
      this.overlay.setStatus(
        'Still waiting for the simulator…',
        'No snapshots have arrived yet. Check the browser console for worker errors.',
      );
    }
    this.stage.render(this.conductor.lastFrameSeconds);
  };

  private onData(): void {
    this.hasData = true;
    if (!this.failed) this.overlay.setStatus(null);
  }

  private fail(title: string, err: unknown): void {
    this.failed = true;
    const detail = err instanceof Error ? err.message : String(err);
    console.error(`[octane] ${title}`, err);
    this.overlay.setStatus(title, detail, true);
  }

  private applyView(v: ViewState): void {
    const temperature = v.mode === 'temperature';
    // False colour must reach the screen unchanged to match the legend.
    this.stage.setGrading(temperature ? 'exact' : 'filmic');
    this.overlay.setLegendVisible(temperature);
  }

  /** Keep the chamber in frame when the cylinder is raised/lowered for a new CR. */
  private followCylinder(cr: number, prev: number): void {
    if (this.framing !== 'chamber') return;
    const L = this.engine.layout;
    const dy = L.headY(cr) - L.headY(prev);
    if (Number.isFinite(dy) && dy !== 0) this.stage.shiftView(0, dy, 0);
  }

  /** Recommended view for the current framing, pulled back on narrow (portrait) viewports. */
  private framingView() {
    const v = recommendedCameraView(this.spec, this.engine.compressionRatio, this.framing);
    const k = framingDistanceScale(this.stage.camera.aspect);
    if (k !== 1) v.position.sub(v.target).multiplyScalar(k).add(v.target);
    return v;
  }

  private readonly onKey = (e: KeyboardEvent): void => {
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    const t = e.target as HTMLElement | null;
    if (t && (t.isContentEditable || /^(INPUT|SELECT|TEXTAREA)$/.test(t.tagName))) return;
    if (e.key === 'f' || e.key === 'F') {
      this.setFraming(this.framing === 'chamber' ? 'engine' : 'chamber');
      e.preventDefault();
    }
  };
}

/** Create and start the app; shows a readable message if WebGL is unavailable. */
export function startApp(opts: AppOptions): App | null {
  try {
    const app = new App(opts);
    app.start();
    return app;
  } catch (err) {
    console.error('[octane] failed to start', err);
    const box = document.createElement('div');
    box.style.cssText =
      'position:absolute;inset:0;display:grid;place-items:center;padding:16px;color:#e6e9ee;font:13px/1.5 system-ui,sans-serif;text-align:center';
    const msg = err instanceof Error ? err.message : String(err);
    box.textContent = `Octane could not start: ${msg}. A browser with WebGL 2 and module Web Workers is required.`;
    opts.viewport.appendChild(box);
    return null;
  }
}
