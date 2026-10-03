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
 */
import '../style.css';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { TracePanel, type TraceFrame } from './charts/trace-panel';
import { Controls, type RenderMode } from './controls';
import { CyclesView } from './cycles-view';
import { Hud } from './hud';
import { SidePanel } from './panel';
import { DEFAULT_TIME_SCALE, snapTimeScale, STEP_LARGE_DEG, STEP_SMALL_DEG, stepTimeScale } from './playback';
import { SparkCapture } from './spark-capture';
import { CycleTraceStore, emptySelection } from './trace-store';

export type { RenderMode } from './controls';

export interface UIControllerOptions {
  spec: EngineSpec;
  initialOperatingPoint: OperatingPoint;
  onOperatingPointChange: (patch: Partial<OperatingPoint>) => void;
  onPlaybackChange: (p: { timeScale: number; paused: boolean }) => void;
  onStep: (deg: number) => void;
  onViewChange: (v: { cutaway: boolean; mode: RenderMode }) => void;
  onReset: () => void;
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
  private readonly store = new CycleTraceStore(HISTORY_CYCLES);
  private readonly spark = new SparkCapture();
  private readonly sel = emptySelection();
  private readonly frame: TraceFrame;
  private timeScale = DEFAULT_TIME_SCALE;
  private paused = false;
  private pending: Partial<OperatingPoint> | null = null;
  private flushHandle = 0;
  private readonly ro: ResizeObserver | null = null;
  private disposed = false;
  private lastAdvance = NaN;

  constructor(
    container: HTMLElement,
    private readonly opts: UIControllerOptions,
  ) {
    this.frame = {
      sel: this.sel,
      version: -1,
      playhead: null as unknown as EngineSnapshot, // set before first use in update()
      sparkDeg: -opts.initialOperatingPoint.sparkAdvanceDeg,
      compressionRatio: opts.initialOperatingPoint.compressionRatio,
      spark: null,
      sparkVersion: -1,
      lastCycle: null,
    };
    this.root = document.createElement('div');
    this.root.className = 'oct-ui';
    container.appendChild(this.root);

    const guiHost = document.createElement('div');
    guiHost.className = 'oct-gui-host';
    this.root.appendChild(guiHost);

    this.controls = new Controls(guiHost, opts.spec, opts.initialOperatingPoint, {
      onOperatingPointChange: (patch) => this.queuePatch(patch),
      onViewChange: (v) => opts.onViewChange(v),
      onTogglePause: () => this.togglePause(),
      onTimeScale: (ts) => this.setTimeScale(ts, false),
      onStep: (deg) => this.step(deg),
      onReset: () => this.reset(),
    });

    this.panel = new SidePanel(this.root);
    this.traces = new TracePanel(this.panel.tracesHost, opts.spec, this.panel.scroll);
    this.cycles = new CyclesView(this.panel.cyclesHost);
    this.panel.onChange(() => {
      if (this.panel.open && this.panel.tab === 'traces') this.traces.invalidate();
      if (this.panel.open && this.panel.tab === 'cycles') this.cycles.invalidate();
    });

    this.hud = new Hud(this.root, opts.spec, {
      onTogglePause: () => this.togglePause(),
      onStep: (deg) => this.step(deg),
      onFaster: () => this.setTimeScale(stepTimeScale(this.timeScale, 1)),
      onSlower: () => this.setTimeScale(stepTimeScale(this.timeScale, -1)),
    });
    this.hud.setLastCycle(null);
    this.syncPlaybackViews();
    this.syncOpViews();

    if (typeof ResizeObserver !== 'undefined') {
      this.ro = new ResizeObserver(() => {
        document.documentElement.style.setProperty('--oct-hud-h', `${Math.ceil(this.hud.el.offsetHeight)}px`);
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
    if (this.store.checkDiscontinuity(s.t, s.cycle, s.rpm)) this.spark.clear();
    this.store.ingest(recent, s.t);
    this.spark.ingest(recent, s.t);
    this.store.select(s.t, HISTORY_CYCLES - 1, this.sel);

    this.hud.update(s);
    const op = this.controls.operatingPoint;
    const f = this.frame;
    f.version = this.store.version;
    f.playhead = s;
    f.sparkDeg = -op.sparkAdvanceDeg;
    f.compressionRatio = op.compressionRatio;
    f.spark = this.spark.display;
    f.sparkVersion = this.spark.version;
    f.lastCycle = this.cycles.latest;
    this.traces.update(f, this.panel.open && this.panel.tab === 'traces');
    this.cycles.render(this.panel.open && this.panel.tab === 'cycles');
  }

  pushCycle(c: CycleSummary): void {
    if (this.disposed) return;
    this.cycles.push(c);
    this.hud.setLastCycle(c);
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

  /** Forget chart history (call if the simulation is reset from outside the UI). */
  clearHistory(): void {
    this.store.clear();
    this.spark.clear();
    this.traces.invalidate();
  }

  /** Programmatic playback control (keeps controls and HUD in sync, fires onPlaybackChange). */
  setPlayback(p: Partial<{ timeScale: number; paused: boolean }>): void {
    if (p.timeScale !== undefined) this.timeScale = snapTimeScale(p.timeScale);
    if (p.paused !== undefined) this.paused = p.paused;
    this.emitPlayback();
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    window.removeEventListener('keydown', this.onKey);
    if (this.flushHandle) cancelAnimationFrame(this.flushHandle);
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
      this.hud.setSparkAdvance(op.sparkAdvanceDeg);
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
        handled = false;
    }
    if (handled) {
      e.preventDefault();
      // A focused button would otherwise also "click" on the Space keyup.
      if (t instanceof HTMLButtonElement) t.blur();
    }
  };
}
