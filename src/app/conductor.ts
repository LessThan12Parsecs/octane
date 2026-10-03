/**
 * Frame-by-frame orchestration of the front end, free of DOM and WebGL so it can
 * be unit-tested with fakes:
 *
 *   dtWall → sim.advance(dt, timeScale) → engine.update(s)
 *                                        → gas.update(s, dt, timeScale)
 *                                        → ui.update(s, sim.recentSnapshots(since))
 *
 * plus the UI callbacks (operating point, playback, stepping, view, reset) and
 * cycle summaries → ui.pushCycle.
 *
 * Compression ratio: the mechanism follows the CR *in the snapshot stream*
 * (h_TDC = clearanceHeight − pistonDisplacement), not the slider. The simulator
 * applies a CR change at the next cycle and playback runs behind the worker, so
 * following the data keeps the piston crown, the head and the gas volume in step
 * exactly — the view shows the CR the physics is actually running.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { compressionRatioFromSnapshot, FrameClock, RecentWindow } from './sync';

export type RenderMode = 'physical' | 'temperature';

export interface ViewState {
  cutaway: boolean;
  mode: RenderMode;
}

export interface PlaybackState {
  timeScale: number;
  paused: boolean;
}

/** What the conductor needs from SimClient. */
export interface SimPort {
  advance(dtWall: number, timeScale: number): EngineSnapshot | null;
  stepDegrees(deg: number): EngineSnapshot | null;
  recentSnapshots(sinceT: number): EngineSnapshot[];
  setOperatingPoint(patch: Partial<OperatingPoint>): void;
  reset(): void;
  onCycle(cb: (c: CycleSummary) => void): void;
  readonly playbackTime: number;
}

/** What the conductor needs from EngineModel. */
export interface MechanismPort {
  update(s: EngineSnapshot): void;
  setCompressionRatio(cr: number): void;
  setCutaway(on: boolean): void;
  readonly compressionRatio: number;
}

/** What the conductor needs from CombustionVisuals. */
export interface GasPort {
  update(s: EngineSnapshot, dtWall: number, timeScale: number): void;
  setMode(mode: RenderMode): void;
}

/** What the conductor needs from UIController. */
export interface UiPort {
  update(s: EngineSnapshot, recent: EngineSnapshot[]): void;
  pushCycle(c: CycleSummary): void;
}

export interface ConductorHooks {
  /** The mechanism moved to a new CR (from the data). */
  onCompressionRatio?(cr: number, previous: number): void;
  /** View state changed (cutaway / render mode). */
  onView?(v: ViewState): void;
  /** First snapshot after start or reset arrived. */
  onFirstSnapshot?(s: EngineSnapshot): void;
}

/** Minimum CR change the mechanism reacts to (filters interpolation round-off). */
export const CR_EPSILON = 1e-4;

export class Conductor {
  private ui: UiPort | null = null;
  private readonly clock: FrameClock;
  private readonly window = new RecentWindow();
  private _playback: PlaybackState;
  private _view: ViewState = { cutaway: true, mode: 'physical' };
  private _snapshot: EngineSnapshot | null = null;
  private waitingFirst = true;
  private _frames = 0;
  private _lastDt = 0;

  constructor(
    private readonly spec: EngineSpec,
    private readonly sim: SimPort,
    private readonly engine: MechanismPort,
    private readonly gas: GasPort,
    private readonly hooks: ConductorHooks = {},
    initialPlayback: PlaybackState = { timeScale: 0.02, paused: false },
    maxFrameSeconds = 0.1,
  ) {
    this._playback = { ...initialPlayback };
    this.clock = new FrameClock(maxFrameSeconds);
    sim.onCycle((c) => this.ui?.pushCycle(c));
  }

  /** Connect the UI (created after the conductor because its callbacks point here). */
  attachUi(ui: UiPort): void {
    this.ui = ui;
  }

  get playback(): Readonly<PlaybackState> {
    return this._playback;
  }
  get view(): Readonly<ViewState> {
    return this._view;
  }
  /** Snapshot shown in the latest frame (owned by the SimClient; do not keep). */
  get snapshot(): EngineSnapshot | null {
    return this._snapshot;
  }
  /** Frames rendered with data. */
  get frames(): number {
    return this._frames;
  }
  /** Wall-clock duration of the latest frame, s (clamped). */
  get lastFrameSeconds(): number {
    return this._lastDt;
  }
  /** Time scale actually applied to the playback clock (0 while paused). */
  get effectiveTimeScale(): number {
    return this._playback.paused ? 0 : this._playback.timeScale;
  }

  // ---------------------------------------------------------------- UI callbacks

  readonly handleOperatingPoint = (patch: Partial<OperatingPoint>): void => {
    this.sim.setOperatingPoint(patch);
  };

  readonly handlePlayback = (p: PlaybackState): void => {
    this._playback = { timeScale: Math.max(0, p.timeScale), paused: p.paused };
  };

  readonly handleStep = (deg: number): void => {
    // The UI pauses before stepping; the next frame (advance with timeScale 0) shows the result.
    this.sim.stepDegrees(deg);
  };

  readonly handleView = (v: ViewState): void => {
    this._view = { cutaway: v.cutaway, mode: v.mode };
    this.engine.setCutaway(v.cutaway);
    this.gas.setMode(v.mode);
    this.hooks.onView?.(this._view);
  };

  readonly handleReset = (): void => {
    this.sim.reset();
    this.window.reset();
    this.waitingFirst = true;
    this._snapshot = null;
  };

  // ---------------------------------------------------------------- frame

  /**
   * One animation frame. `nowMs` is the rAF timestamp (ms); returns the snapshot
   * shown (null while no data is available yet).
   */
  frame(nowMs: number): EngineSnapshot | null {
    const dt = this.clock.tick(nowMs);
    this._lastDt = dt;
    const ts = this.effectiveTimeScale;
    const s = this.sim.advance(dt, ts);
    this._snapshot = s;
    if (!s) return null;
    if (this.waitingFirst) {
      this.waitingFirst = false;
      this.hooks.onFirstSnapshot?.(s);
    }
    this.syncCompressionRatio(s);
    this.engine.update(s);
    this.gas.update(s, dt, ts);
    if (this.ui) {
      const since = this.window.next(this.sim.playbackTime);
      this.ui.update(s, this.sim.recentSnapshots(since));
    }
    this._frames++;
    return s;
  }

  private syncCompressionRatio(s: EngineSnapshot): void {
    const cr = compressionRatioFromSnapshot(this.spec, s);
    if (!Number.isFinite(cr)) return;
    const [lo, hi] = this.spec.geometry.compressionRatioRange;
    const target = Math.min(hi, Math.max(lo, cr));
    const prev = this.engine.compressionRatio;
    if (Math.abs(target - prev) <= CR_EPSILON) return;
    this.engine.setCompressionRatio(target);
    this.hooks.onCompressionRatio?.(this.engine.compressionRatio, prev);
  }
}
