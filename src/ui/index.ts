/**
 * UI layer: lil-gui controls (top-left), collapsible chart panel (right), HUD strip
 * (bottom). Consumes only EngineSnapshot / CycleSummary / EngineSpec / OperatingPoint.
 *
 * Layout contract (src/style.css): the UI root is a fixed full-window overlay that
 * lets pointer events through to the 3D canvas. Put the canvas (or its wrapper) in an
 * element with class `oct-viewport` and it will be sized to the area not covered by
 * the panel and HUD (the panel toggle dispatches a window 'resize' event afterwards).
 *
 * `update(s, recent)`: `recent` may be any window of raw (non-interpolated) snapshots
 * up to the playback time — e.g. `simClient.recentSnapshots(previousPlaybackTime)`.
 * Each snapshot is ingested once; those later than `s.t` are ignored until reached.
 *
 * Multi-cylinder engines: every cylinder keeps its own trace history and spark
 * capture (fed with cylinder views, cylinder-view.ts), so switching the FOCUS
 * cylinder (HUD strip, chart toolbar, Cycles tab, or setFocusCylinder) shows its
 * history at once; a further store holds every cylinder's pressure against the
 * engine angle for the all-cylinders chart.
 */
import '../style.css';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import type { EngineDefinition } from '../physics/engines/index';
import { TracePanel, type TraceFrame } from './charts/trace-panel';
import { Controls, type RenderMode } from './controls';
import { CyclesView } from './cycles-view';
import { clampCylinder, createCylinderViewTarget, cylinderCount, cylinderView, CylinderViewPool } from './cylinder-view';
import { definitionForSpec } from './engine-ui';
import { Hud } from './hud';
import { SidePanel } from './panel';
import { DEFAULT_TIME_SCALE, snapTimeScale, STEP_LARGE_DEG, STEP_SMALL_DEG, stepTimeScale } from './playback';
import { SparkCapture } from './spark-capture';
import { CycleTraceStore, cylinderPressureSampler, emptySelection, type TraceSelection } from './trace-store';

export type { RenderMode } from './controls';
export { definitionForSpec, engineTitle } from './engine-ui';

export interface UIControllerOptions {
  spec: EngineSpec;
  initialOperatingPoint: OperatingPoint;
  onOperatingPointChange: (patch: Partial<OperatingPoint>) => void;
  onPlaybackChange: (p: { timeScale: number; paused: boolean }) => void;
  onStep: (deg: number) => void;
  onViewChange: (v: { cutaway: boolean; mode: RenderMode }) => void;
  onReset: () => void;
  /** Engine registry entry (controls profile, presets, label); default: the spec's entry. */
  engine?: EngineDefinition;
  /** Initial focus cylinder, 0-based (multi-cylinder engines). */
  focusCylinder?: number;
  /** The user picked a focus cylinder (strip, chart toolbar or Cycles tab). */
  onFocusCylinderChange?: (index: number) => void;
}

/** Previous cycles kept for the charts. */
const HISTORY_CYCLES = 8;

export class UIController {
  private readonly root: HTMLElement;
  private readonly controls: Controls;
  private readonly hud: Hud;
  private readonly panel: SidePanel;
  private readonly traces: TracePanel;
  private readonly cycles: CyclesView;
  /** One trace store and spark capture per cylinder. */
  private readonly stores: CycleTraceStore[] = [];
  private readonly sparks: SparkCapture[] = [];
  /** Playhead view of each cylinder (reused). */
  private readonly views: EngineSnapshot[] = [];
  private readonly viewPool = new CylinderViewPool();
  /** All cylinders' pressure vs the engine angle (multi-cylinder engines only). */
  private readonly overlay: CycleTraceStore | null = null;
  private readonly overlaySel: TraceSelection = emptySelection();
  private readonly nCyl: number;
  private focus = 0;
  private readonly sel = emptySelection();
  private readonly frame: TraceFrame;
  private timeScale = DEFAULT_TIME_SCALE;
  private paused = false;
  private pending: Partial<OperatingPoint> | null = null;
  private flushHandle = 0;
  private readonly ro: ResizeObserver | null = null;
  private hudSizeHandle = 0;
  private disposed = false;
  private lastAdvance = NaN;
  private lastMarker = NaN;

  constructor(
    container: HTMLElement,
    private readonly opts: UIControllerOptions,
  ) {
    const spec = opts.spec;
    const def = opts.engine ?? definitionForSpec(spec);
    const n = (this.nCyl = cylinderCount(spec));
    this.focus = clampCylinder(opts.focusCylinder ?? 0, n);
    for (let i = 0; i < n; i++) {
      this.stores.push(new CycleTraceStore(HISTORY_CYCLES));
      this.sparks.push(new SparkCapture());
      this.views.push(createCylinderViewTarget());
    }
    if (n > 1) this.overlay = new CycleTraceStore(HISTORY_CYCLES, 2048, cylinderPressureSampler(n));
    const shower = spec.ignition.type === 'trembler-magneto';

    this.frame = {
      sel: this.sel,
      version: -1,
      playhead: null as unknown as EngineSnapshot, // set before first use in update()
      sparkDeg: -opts.initialOperatingPoint.sparkAdvanceDeg,
      compressionRatio: opts.initialOperatingPoint.compressionRatio,
      spark: null,
      sparkVersion: -1,
      lastCycle: null,
      overlay: this.overlay ? this.overlaySel : null,
      overlayVersion: -1,
      engineTheta: -360,
    };
    this.root = document.createElement('div');
    this.root.className = 'oct-ui';
    container.appendChild(this.root);

    const guiHost = document.createElement('div');
    guiHost.className = 'oct-gui-host';
    this.root.appendChild(guiHost);

    this.controls = new Controls(
      guiHost,
      spec,
      opts.initialOperatingPoint,
      {
        onOperatingPointChange: (patch) => this.queuePatch(patch),
        onViewChange: (v) => opts.onViewChange(v),
        onTogglePause: () => this.togglePause(),
        onTimeScale: (ts) => this.setTimeScale(ts, false),
        onStep: (deg) => this.step(deg),
        onReset: () => this.reset(),
      },
      def,
    );

    const pick = (i: number): void => {
      this.setFocusCylinder(i);
      opts.onFocusCylinderChange?.(this.focus);
    };
    this.panel = new SidePanel(this.root);
    this.traces = new TracePanel(this.panel.tracesHost, spec, this.panel.scroll, { cylinders: n, focus: this.focus, onFocus: pick, shower });
    this.cycles = new CyclesView(this.panel.cyclesHost, { cylinders: n, cylinder: this.focus, onCylinder: pick });
    this.panel.onChange(() => {
      if (this.panel.open && this.panel.tab === 'traces') this.traces.invalidate();
      if (this.panel.open && this.panel.tab === 'cycles') this.cycles.invalidate();
    });

    this.hud = new Hud(
      this.root,
      spec,
      {
        onTogglePause: () => this.togglePause(),
        onStep: (deg) => this.step(deg),
        onFaster: () => this.setTimeScale(stepTimeScale(this.timeScale, 1)),
        onSlower: () => this.setTimeScale(stepTimeScale(this.timeScale, -1)),
        onFocusCylinder: pick,
      },
      { focus: this.focus, brake: n > 1 || def.ui.loadModels.includes('vehicle'), vehicle: !!spec.vehicle },
    );
    this.hud.setLastCycle(null);
    this.syncPlaybackViews();
    this.syncOpViews();

    if (typeof ResizeObserver !== 'undefined') {
      this.ro = new ResizeObserver(() => {
        // On the next frame: the variable re-lays out the page, which inside the observer callback trips
        // the browser's "ResizeObserver loop completed with undelivered notifications" error.
        if (this.hudSizeHandle) return;
        this.hudSizeHandle = requestAnimationFrame(() => {
          this.hudSizeHandle = 0;
          if (this.disposed) return;
          document.documentElement.style.setProperty('--oct-hud-h', `${Math.ceil(this.hud.el.offsetHeight)}px`);
        });
      });
      this.ro.observe(this.hud.el);
    }
    window.addEventListener('keydown', this.onKey);

    // Tell the integrator the initial playback/view state once it has finished wiring.
    queueMicrotask(() => {
      if (this.disposed) return;
      opts.onPlaybackChange({ timeScale: this.timeScale, paused: this.paused });
      opts.onViewChange({ ...this.controls.view });
    });
  }

  // ------------------------------------------------------------------ public API

  /** Call every animation frame with the playback snapshot and recent raw snapshots. */
  update(s: EngineSnapshot, recent: EngineSnapshot[]): void {
    if (this.disposed) return;
    this.flush();
    const n = this.nCyl;
    for (let i = 0; i < n; i++) {
      const v = cylinderView(s, i, this.views[i]);
      const store = this.stores[i];
      if (store.checkDiscontinuity(v.t, v.cycle, v.rpm)) this.sparks[i].clear();
      // Single-cylinder streams (no cylinders[]) map to the raw snapshots themselves.
      const rv = s.cylinders ? this.viewPool.map(recent, i) : recent;
      store.ingest(rv, s.t);
      this.sparks[i].ingest(rv, s.t);
      if (!s.cylinders) break;
    }
    const ov = this.overlay;
    if (ov) {
      ov.checkDiscontinuity(s.t, s.cycle, s.rpm);
      ov.ingest(recent, s.t);
      ov.select(s.t, 0, this.overlaySel);
    }
    // A stream without cylinders[] only feeds cylinder 1's store.
    const f = s.cylinders ? this.focus : 0;
    const view = cylinderView(s, f, this.views[f]);
    this.stores[f].select(s.t, HISTORY_CYCLES - 1, this.sel);

    this.hud.update(view, s);
    const op = this.controls.operatingPoint;
    const spark = this.sparks[f];
    const fr = this.frame;
    fr.version = this.stores[f].version;
    fr.playhead = view;
    fr.sparkDeg = this.sparkMarkerDeg(f, op);
    fr.compressionRatio = op.compressionRatio;
    fr.spark = spark.display;
    fr.sparkVersion = spark.version;
    fr.lastCycle = this.cycles.latestFor(f);
    fr.overlayVersion = ov ? ov.version : -1;
    fr.engineTheta = s.thetaDeg;
    if (fr.sparkDeg !== this.lastMarker) {
      this.lastMarker = fr.sparkDeg;
      this.hud.setSparkAdvance(-fr.sparkDeg);
    }
    this.traces.update(fr, this.panel.open && this.panel.tab === 'traces');
    this.cycles.render(this.panel.open && this.panel.tab === 'cycles');
  }

  pushCycle(c: CycleSummary): void {
    if (this.disposed) return;
    this.cycles.push(c);
    if ((c.cylinder ?? 0) === this.focus) this.hud.setLastCycle(c);
    if (c.engine) this.hud.setEngineCycle(c.engine);
  }

  /** Current playback state (also delivered through onPlaybackChange). */
  get playback(): { timeScale: number; paused: boolean } {
    return { timeScale: this.timeScale, paused: this.paused };
  }

  /** Current view state (also delivered through onViewChange). */
  get view(): { cutaway: boolean; mode: RenderMode } {
    return { ...this.controls.view };
  }

  /** Current operating point as shown by the controls. */
  get operatingPoint(): OperatingPoint {
    return this.controls.operatingPoint;
  }

  /** Focus cylinder, 0-based. */
  get focusCylinder(): number {
    return this.focus;
  }

  /** Show cylinder `index` (0-based) in the HUD readouts, θ / p–V / spark charts and Cycles tab (no callback). */
  setFocusCylinder(index: number): void {
    const i = clampCylinder(index, this.nCyl);
    if (i === this.focus) return;
    this.focus = i;
    this.hud.setFocus(i);
    this.traces.setFocus(i);
    this.cycles.setCylinder(i);
    this.hud.setLastCycle(this.cycles.latestFor(i));
    this.lastMarker = NaN;
  }

  /** Forget chart history (call if the simulation is reset from outside the UI). */
  clearHistory(): void {
    for (const st of this.stores) st.clear();
    for (const sp of this.sparks) sp.clear();
    this.overlay?.clear();
    this.traces.invalidate();
  }

  /** Programmatic playback control (keeps controls and HUD in sync, fires onPlaybackChange). */
  setPlayback(p: Partial<{ timeScale: number; paused: boolean }>): void {
    if (p.timeScale !== undefined) this.timeScale = snapTimeScale(p.timeScale);
    if (p.paused !== undefined) this.paused = p.paused;
    this.emitPlayback();
  }

  /** Programmatic view control (keeps the controls in sync, fires onViewChange). */
  setView(v: { cutaway: boolean; mode: RenderMode }): void {
    this.controls.setView(v);
    this.opts.onViewChange({ ...this.controls.view });
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    window.removeEventListener('keydown', this.onKey);
    if (this.flushHandle) cancelAnimationFrame(this.flushHandle);
    if (this.hudSizeHandle) cancelAnimationFrame(this.hudSizeHandle);
    this.ro?.disconnect();
    this.traces.dispose();
    this.cycles.dispose();
    this.panel.dispose();
    this.hud.dispose();
    this.controls.dispose();
    this.root.remove();
    document.documentElement.style.removeProperty('--oct-hud-h');
  }

  // ------------------------------------------------------------------ internals

  /**
   * Spark marker (θ of the focus cylinder): the first gap breakdown of its latest ignition event as the
   * physics produced it — with a magneto and trembler coils the spark lags the lever (timer make) by the
   * coil's rpm-dependent firing time — falling back to the lever / commanded angle before any spark.
   */
  private sparkMarkerDeg(cyl: number, op: OperatingPoint): number {
    const ev = this.sparks[cyl].display;
    if (ev && Number.isFinite(ev.firstSparkDeg)) return ev.firstSparkDeg;
    return -op.sparkAdvanceDeg;
  }

  private queuePatch(patch: Partial<OperatingPoint>): void {
    this.pending = { ...(this.pending ?? {}), ...patch };
    this.syncOpViews();
    // Slider drags fire many changes per frame: coalesce to one message per frame.
    if (!this.flushHandle && typeof requestAnimationFrame !== 'undefined') {
      this.flushHandle = requestAnimationFrame(() => {
        this.flushHandle = 0;
        this.flush();
      });
    }
  }

  private flush(): void {
    if (!this.pending) return;
    const p = this.pending;
    this.pending = null;
    this.opts.onOperatingPointChange(p);
  }

  private syncOpViews(): void {
    const op = this.controls.operatingPoint;
    if (op.sparkAdvanceDeg !== this.lastAdvance) {
      this.lastAdvance = op.sparkAdvanceDeg;
      // Before the first spark of the focus cylinder the marker shows the lever.
      if (!this.sparks[this.focus].display) {
        this.lastMarker = -op.sparkAdvanceDeg;
        this.hud.setSparkAdvance(op.sparkAdvanceDeg);
      }
    }
  }

  private togglePause(): void {
    this.paused = !this.paused;
    this.emitPlayback();
  }

  private setTimeScale(ts: number, syncSlider = true): void {
    this.timeScale = snapTimeScale(ts);
    this.emitPlayback(syncSlider);
  }

  private step(deg: number): void {
    if (!this.paused) {
      this.paused = true;
      this.emitPlayback();
    }
    this.opts.onStep(deg);
  }

  private reset(): void {
    this.flush();
    this.clearHistory();
    this.cycles.clear();
    this.hud.setLastCycle(null);
    this.hud.setEngineCycle(null);
    this.opts.onReset();
  }

  private emitPlayback(syncSlider = true): void {
    this.syncPlaybackViews(syncSlider);
    this.opts.onPlaybackChange({ timeScale: this.timeScale, paused: this.paused });
  }

  private syncPlaybackViews(syncSlider = true): void {
    if (syncSlider) this.controls.setPlayback(this.timeScale, this.paused);
    else this.controls.setPlaybackButtons(this.paused);
    this.hud.setPlayback(this.timeScale, this.paused);
  }

  private readonly onKey = (e: KeyboardEvent): void => {
    if (e.metaKey || e.ctrlKey || e.altKey) return;
    const t = e.target as HTMLElement | null;
    if (t && (t.isContentEditable || /^(INPUT|SELECT|TEXTAREA)$/.test(t.tagName))) return;
    let handled = true;
    switch (e.key) {
      case ' ':
      case 'Spacebar':
        this.togglePause();
        break;
      case 'ArrowLeft':
        this.step(-(e.shiftKey ? STEP_LARGE_DEG : STEP_SMALL_DEG));
        break;
      case 'ArrowRight':
        this.step(e.shiftKey ? STEP_LARGE_DEG : STEP_SMALL_DEG);
        break;
      case '[':
        this.setTimeScale(stepTimeScale(this.timeScale, -1));
        break;
      case ']':
        this.setTimeScale(stepTimeScale(this.timeScale, 1));
        break;
      case 'c':
      case 'C':
        this.panel.setOpen(!this.panel.open);
        break;
      default:
        // 1–9: focus that cylinder (multi-cylinder engines).
        if (this.nCyl > 1 && /^[1-9]$/.test(e.key) && Number(e.key) <= this.nCyl) {
          this.setFocusCylinder(Number(e.key) - 1);
          this.opts.onFocusCylinderChange?.(this.focus);
        } else handled = false;
    }
    if (handled) {
      e.preventDefault();
      // A focused button would otherwise also "click" on the Space keyup.
      if (t instanceof HTMLButtonElement) t.blur();
    }
  };
}
