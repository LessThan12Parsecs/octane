/**
 * Octane front end: wires the simulation client, the 3D stage (engine
 * mechanism + in-cylinder combustion visuals) and the UI together and runs the
 * animation loop.
 *
 *   SimClient (worker) ──snapshots──▶ Conductor ──▶ EngineRenderModel.update
 *                                               ├─▶ CombustionVisuals.update (one per cylinder)
 *                                               └─▶ UIController.update (+ pushCycle)
 *   UIController callbacks ──▶ Conductor (operating point, playback, step, view, reset)
 *
 * The App is a SHELL that lives as long as the page — Stage (one WebGL context),
 * viewport overlay (engine picker, framing, status), rAF loop, keyboard and
 * camera framing — around a disposable EngineSession (session.ts) holding
 * everything engine-specific. setEngine(id) builds the new session first and
 * swaps it in atomically, so a failure (e.g. an engine whose 3D model is not
 * available) leaves the running engine untouched.
 *
 * Starts on the engine from ?engine=, the remembered choice or the default (the
 * CFR at its Research-method conditions) in slow motion, framed on the chamber. A
 * link with operating-point knobs but no ?engine= predates the engine picker and
 * opens the default engine (resolveStartEngine); after start-up the address bar
 * names the engine that started whenever it carries other parameters.
 *
 * Render-loop errors halt the animation loop (AnimationLoop, loop.ts) behind a
 * failure card; "Back to <engine>" (a successful session swap) and the UI's Reset
 * restart it, and an error that repeats right away halts it again.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineSnapshot } from '../physics/core/snapshot';
import { DEFAULT_ENGINE_ID, ENGINES, type EngineDefinition } from '../physics/engines/index';
import { createTemperatureLegendElement, Stage, type CameraView } from '../render/index';
import { definitionForSpec, engineTitle } from '../ui/index';
import { loadPref, savePref } from '../ui/prefs';
import type { SimulatorOptions } from '../worker/protocol';
import type { PlaybackState, ViewState } from './conductor';
import { engineChoices, initialOperatingPoint, isEngineId, resolveStartEngine } from './engines';
import { AnimationLoop } from './loop';
import { ViewportOverlay, type StatusAction } from './overlay';
import './register-engines';
import { EngineSession } from './session';
import type { SimClientConfig } from './sim-client';
import {
  clampOperatingPoint,
  errorSummary,
  framingDistanceScale,
  parseUrlOptions,
  searchWithEngine,
  searchWithStartEngine,
  type Framing,
} from './sync';

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
/** Consecutive starved / fed frames before the simulation-limited hint appears / goes. */
const LIMITED_FRAMES = 45;
/** localStorage key (prefs.ts) of the remembered engine. */
const ENGINE_PREF = 'engine';

export interface AppOptions {
  /** Element with class `oct-viewport` that receives the canvas. */
  viewport: HTMLElement;
  /** Container for the UI overlay (controls, charts, HUD). */
  uiContainer: HTMLElement;
  /** Engine registry id to start with (default: ?engine=, the remembered choice, else the CFR). */
  engineId?: string;
  /** Ad-hoc spec to run instead of a registry engine (tests, experiments). */
  spec?: EngineSpec;
  /** Initial operating point (default: the engine's default operating point). */
  operatingPoint?: OperatingPoint;
  simulatorOptions?: SimulatorOptions;
  /** `location.search` for the URL knobs (?engine=&cr=&on=&rpm=&spark=&phi=&ts=&paused=&view=&cyl=&sim=). */
  search?: string;
  /**
   * Write the engine back to the address bar (history.replaceState) on a switch, and after start-up when
   * the URL carries parameters but does not name the engine that started (default true).
   */
  syncUrl?: boolean;
  /** SimClient configuration of every session (e.g. a custom worker factory for dev pages and tests). */
  simClientConfig?: SimClientConfig;
}

export class App {
  readonly stage: Stage;
  private session: EngineSession;
  private readonly overlay: ViewportOverlay;
  private readonly uiContainer: HTMLElement;
  private readonly simulatorOptions: SimulatorOptions;
  private readonly simClientConfig: SimClientConfig | undefined;
  private readonly syncUrl: boolean;
  private framing: Framing;
  private focus = 0;
  private readonly loop: AnimationLoop;
  private startedAt = NaN;
  private hasData = false;
  private failed = false;
  private disposed = false;
  /** Last engine that produced data (target of "switch back" after a failure). */
  private lastGoodEngine: string | null = null;
  private starvedFrames = 0;
  private fedFrames = 0;
  private limitedShown = false;
  /** A dismissable notice (engine not available) is showing: data arriving must not hide it. */
  private notice = false;

  constructor(opts: AppOptions) {
    const search = opts.search ?? '';
    const url = parseUrlOptions(search);
    this.uiContainer = opts.uiContainer;
    this.syncUrl = opts.syncUrl ?? true;
    this.simClientConfig = opts.simClientConfig;
    this.framing = url.framing ?? 'chamber';
    const base = opts.simulatorOptions ?? DEFAULT_SIMULATOR_OPTIONS;
    // ?sim=mock: the lightweight single-cylinder stand-in (extra fields pass through to the worker).
    this.simulatorOptions = url.simulator ? ({ ...base, simulator: url.simulator } as SimulatorOptions) : base;

    const requested: EngineDefinition = opts.spec
      ? definitionForSpec(opts.spec)
      : resolveStartEngine(opts.engineId, search, loadPref<string | null>(ENGINE_PREF, null));
    this.focus = Math.max(0, Math.min(requested.spec.cylinders - 1, url.focusCylinder ?? 0));
    const op = initialOperatingPoint(requested, url.op, opts.operatingPoint);
    const playback: PlaybackState = { timeScale: url.timeScale ?? DEFAULT_TIME_SCALE, paused: url.paused ?? false };

    // ---- shell ----
    this.loop = new AnimationLoop(this.frame, (err, consecutive) =>
      this.fail(consecutive > 1 ? 'The render loop stopped again on an error.' : 'The render loop stopped on an error.', err),
    );
    this.stage = new Stage(opts.viewport);
    let session: EngineSession | null = null;
    let startError: { def: EngineDefinition; err: unknown } | null = null;
    let overlay: ViewportOverlay | null = null;
    try {
      const legend = createTemperatureLegendElement();
      this.overlay = overlay = new ViewportOverlay(opts.viewport, {
        onFraming: (f) => this.setFraming(f),
        legend: legend.element,
        engines: opts.spec ? [] : engineChoices(),
        engineId: requested.id,
        onEngine: (id) => this.setEngine(id),
      });
      this.overlay.setFraming(this.framing);

      // ---- first engine (fall back to the default one if it cannot be built) ----
      try {
        session = this.buildSession(requested, op, playback);
      } catch (err) {
        if (opts.spec || requested.id === DEFAULT_ENGINE_ID) throw err;
        startError = { def: requested, err };
        const def = ENGINES[DEFAULT_ENGINE_ID];
        this.focus = 0;
        session = this.buildSession(def, initialOperatingPoint(def), playback);
      }
    } catch (err) {
      overlay?.dispose();
      this.stage.dispose();
      throw err;
    }
    this.session = session;
    this.afterSwap(null, url.timeScale !== undefined || url.paused !== undefined ? playback : null, null);
    if (startError) this.reportEngineError(startError.def, startError.err);
    else if (!opts.spec) {
      // A copied address bar must reopen this engine, not the copier's remembered one.
      const next = searchWithStartEngine(typeof location !== 'undefined' ? location.search : search, session.def.id);
      if (next !== null) this.replaceSearch(next);
    }

    window.addEventListener('keydown', this.onKey);
  }

  // ------------------------------------------------------------------ public API

  /** The running engine's registry entry. */
  get engine(): EngineDefinition {
    return this.session.def;
  }
  get engineId(): string {
    return this.session.def.id;
  }
  get spec(): EngineSpec {
    return this.session.spec;
  }
  /** Current session parts (replaced by setEngine; do not keep references). */
  get sim() {
    return this.session.sim;
  }
  get conductor() {
    return this.session.conductor;
  }
  get ui() {
    return this.session.ui;
  }
  get mechanism() {
    return this.session.engine;
  }
  get combustion() {
    return this.session.gas;
  }

  /** Start the animation loop (idempotent). */
  start(): void {
    if (this.disposed) return;
    this.loop.start();
  }

  /** Stop the animation loop (the worker keeps its buffer). */
  stop(): void {
    this.loop.stop();
  }

  get currentFraming(): Framing {
    return this.framing;
  }

  get focusCylinder(): number {
    return this.focus;
  }

  /** Frame the combustion chamber (close-up through the cutaway) or the whole engine. */
  setFraming(f: Framing, animate = true): void {
    this.framing = f;
    this.overlay.setFraming(f);
    this.stage.setFraming(this.framingView(), animate);
  }

  /**
   * Focus cylinder (0-based): UI readouts and charts, the cutaway section (engines whose section can
   * move: the Model T's quarter section opens the focus cylinder's chamber), the 'chamber' close-up, and
   * the featured in-cylinder visuals (flow tracers and chamber light follow the focus).
   */
  setFocusCylinder(index: number, animate = true): void {
    const n = this.session.spec.cylinders;
    const i = Math.max(0, Math.min(n - 1, Math.floor(index)));
    if (i === this.focus) return;
    this.focus = i;
    this.session.ui.setFocusCylinder(i);
    try {
      this.session.setFocusCylinder(i);
    } catch (err) {
      // Section and/or visuals stay where they were; readouts and framing still follow the focus.
      console.error('[octane] could not move the cutaway section / flow tracers to the focus cylinder', err);
    }
    if (this.framing === 'chamber') this.stage.setFraming(this.framingView(), animate);
  }

  /**
   * Switch to another registry engine, keeping playback speed, view mode and camera framing. Returns
   * false (and keeps the running engine) if the id is unknown or the new engine cannot be built.
   */
  setEngine(id: string): boolean {
    if (this.disposed) return false;
    if (!isEngineId(id)) return false;
    const def = ENGINES[id as keyof typeof ENGINES];
    const old = this.session;
    if (def.id === old.def.id) return true;
    const playback = { ...old.conductor.playback };
    const view = { ...old.conductor.view };
    const prevFocus = this.focus;
    this.focus = Math.min(this.focus, def.spec.cylinders - 1);
    this.overlay.setEngineEnabled(false);
    let next: EngineSession;
    try {
      next = this.buildSession(def, initialOperatingPoint(def), playback);
    } catch (err) {
      this.focus = prevFocus;
      this.overlay.setEngineEnabled(true);
      this.overlay.setEngine(old.def.id);
      this.reportEngineError(def, err);
      return false;
    }
    this.session = next;
    old.dispose();
    this.afterSwap(old.def, playback, view);
    this.overlay.setEngineEnabled(true);
    return true;
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.loop.dispose();
    window.removeEventListener('keydown', this.onKey);
    this.session.dispose();
    this.overlay.dispose();
    this.stage.dispose();
  }

  // ------------------------------------------------------------------

  private buildSession(def: EngineDefinition, op: OperatingPoint, playback: PlaybackState): EngineSession {
    let session: EngineSession | null = null;
    session = new EngineSession({
      stage: this.stage,
      uiContainer: this.uiContainer,
      def,
      operatingPoint: clampOperatingPoint(def.spec, op),
      simulatorOptions: this.simulatorOptions,
      playback,
      focusCylinder: this.focus,
      hooks: {
        onCompressionRatio: (cr, prev) => this.followCylinder(cr, prev),
        onView: (v) => this.applyView(v),
        onFirstSnapshot: () => this.onData(session),
      },
      onReset: () => {
        if (session !== this.session) return;
        this.hasData = false;
        this.startedAt = NaN;
        if (this.failed) {
          // The simulator restarts (re-initialised if it could not be built): give it a fresh chance.
          this.failed = false;
          if (!this.notice) this.overlay.setStatus('Restarting the simulator…', 'Running the first cycle again with the current settings.');
        }
        this.loop.resume(); // a render-loop error halted it: try again (it halts again if the error repeats)
      },
      onFocusCylinder: (i) => this.setFocusCylinder(i),
      simClientConfig: this.simClientConfig,
    });
    const s = session;
    s.sim.onError((message) => {
      if (this.session === s) this.fail('The simulator stopped with an error.', new Error(message));
    });
    return s;
  }

  /** Stage, overlay, title and persistence after a session became current. */
  private afterSwap(previous: EngineDefinition | null, playback: PlaybackState | null, view: ViewState | null): void {
    const s = this.session;
    this.hasData = false;
    this.failed = false;
    this.notice = false;
    this.startedAt = NaN;
    this.starvedFrames = 0;
    this.fedFrames = 0;
    this.overlay.setLimited((this.limitedShown = false));
    this.stage.fitToBounds(s.bounds());
    this.stage.setFraming(this.framingView(), false);
    this.overlay.setEngine(s.def.id);
    this.overlay.setStatus(
      'Starting the simulator…',
      s.spec.cylinders > 1
        ? `Spawning the physics worker and running the first cycles of the ${s.def.label} (${s.spec.cylinders} cylinders).`
        : 'Spawning the physics worker and running the first cycle.',
    );
    document.title = engineTitle(s.def);
    // The UI announces its own initial playback / view state in a microtask; re-apply ours after it.
    if (playback || view) {
      queueMicrotask(() => {
        if (this.disposed || this.session !== s) return;
        if (playback) s.ui.setPlayback(playback);
        if (view) s.ui.setView(view);
      });
    }
    if (previous) {
      savePref(ENGINE_PREF, s.def.id);
      if (typeof location !== 'undefined') this.replaceSearch(searchWithEngine(location.search, s.def.id));
    }
    // After a render-loop error the card offered this swap as the way out: run the new session.
    this.loop.resume();
  }

  /** history.replaceState with a new search string (when syncUrl; path and hash kept). */
  private replaceSearch(search: string): void {
    if (!this.syncUrl || typeof history === 'undefined' || typeof location === 'undefined') return;
    try {
      history.replaceState(history.state, '', `${location.pathname}${search}${location.hash}`);
    } catch {
      /* sandboxed frames may refuse; the URL is a convenience */
    }
  }

  /** One animation frame (AnimationLoop; a throw halts the loop and shows the failure card). */
  private readonly frame = (now: number): void => {
    if (Number.isNaN(this.startedAt)) this.startedAt = now;
    const session = this.session;
    const s: EngineSnapshot | null = session.conductor.frame(now);
    if (!s && !this.hasData && !this.failed && (now - this.startedAt) / 1000 > STALL_WARNING_SECONDS) {
      this.overlay.setStatus(
        'Still waiting for the simulator…',
        'No snapshots have arrived yet. Check the browser console for worker errors.',
        false,
        this.switchBackActions(),
      );
    }
    this.updateLimited(s !== null && session.sim.starved);
    this.stage.render(session.conductor.lastFrameSeconds);
  };

  /** Show the simulation-limited hint after a sustained stall of playback behind the worker. */
  private updateLimited(starved: boolean): void {
    if (starved) {
      this.starvedFrames++;
      this.fedFrames = 0;
    } else {
      this.fedFrames++;
      this.starvedFrames = 0;
    }
    const show = this.limitedShown ? this.fedFrames < LIMITED_FRAMES : this.starvedFrames >= LIMITED_FRAMES;
    if (show !== this.limitedShown) this.overlay.setLimited((this.limitedShown = show));
  }

  private onData(session: EngineSession | null): void {
    if (!session || session !== this.session) return;
    this.hasData = true;
    this.lastGoodEngine = session.def.id;
    if (!this.failed && !this.notice) this.overlay.setStatus(null);
  }

  private fail(title: string, err: unknown): void {
    this.failed = true;
    this.notice = false;
    const detail = errorSummary(err instanceof Error ? err.message : String(err));
    console.error(`[octane] ${title}`, err);
    this.overlay.setStatus(title, detail, true, this.switchBackActions());
  }

  /** "Back to <engine>" for the last engine that worked (or the default engine). */
  private switchBackActions(): StatusAction[] {
    const current = this.session.def.id;
    const target = this.lastGoodEngine && this.lastGoodEngine !== current ? this.lastGoodEngine : current !== DEFAULT_ENGINE_ID ? DEFAULT_ENGINE_ID : null;
    if (!target || !isEngineId(target)) return [];
    const def = ENGINES[target as keyof typeof ENGINES];
    return [{ label: `Back to ${def.label}`, onClick: () => this.setEngine(def.id) }];
  }

  /** An engine could not be built: say so; the running engine continues. */
  private reportEngineError(def: EngineDefinition, err: unknown): void {
    const detail = errorSummary(err instanceof Error ? err.message : String(err));
    console.error(`[octane] could not start the ${def.label}`, err);
    this.notice = true;
    this.overlay.setStatus(`The ${def.label} is not available.`, `${detail}\nStill running the ${this.session.def.label}.`, true, [
      {
        label: 'Dismiss',
        onClick: () => {
          this.notice = false;
          this.overlay.setStatus(this.hasData || this.failed ? null : 'Starting the simulator…');
        },
      },
    ]);
  }

  private applyView(v: ViewState): void {
    const temperature = v.mode === 'temperature';
    // False colour must reach the screen unchanged to match the legend.
    this.stage.setGrading(temperature ? 'exact' : 'filmic');
    this.overlay.setLegendVisible(temperature);
  }

  /** Keep the chamber in frame when the cylinder is raised/lowered for a new CR (variable-CR engines). */
  private followCylinder(cr: number, prev: number): void {
    if (this.framing !== 'chamber') return;
    const dy = this.session.engine.chamberShift(cr, prev);
    if (Number.isFinite(dy) && dy !== 0) this.stage.shiftView(0, dy, 0);
  }

  /** Recommended view for the current framing, pulled back on narrow (portrait) viewports. */
  private framingView(): CameraView {
    const v = this.session.cameraView(this.framing, this.focus);
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
