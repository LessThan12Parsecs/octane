/**
 * In-cylinder gas, flame, spark and flow visuals.
 *
 * `CombustionVisuals.root` is parented to a cylinder frame of the engine model
 * (EngineRenderModel.cylinderFrames[i]), so all of this works in CYLINDER-frame
 * coordinates (m): origin at the centre of the head fire-deck face (L-head: on
 * the bore axis at the roof of the head cavity), +y toward the head; the gas is
 * the disc x²+z² ≤ R², −h ≤ y ≤ 0, or for an L-head that bore column plus the
 * valve pocket (chamber.ts).
 *
 * Multi-cylinder engines: one instance per cylinder (`{ cylinder: i }`) reads
 * `s.cylinders[i]` (its local crank angle, flame, spark, knock …) through a
 * zero-allocation view, falling back to the top-level fields; flow tracers run
 * only where `tracers` is not false (they dominate the CPU cost).
 *
 * Everything is driven by EngineSnapshot fields (+ EngineSpec geometry):
 *  - gas volume (one ray-marched pass): flame = sphere(flame.center,
 *    flame.radius) ∩ chamber with a wrinkled turbulent brush; physical mode
 *    shows chemiluminescence ∝ heat release per front area, faint gray-body
 *    burned gas, end-gas autoignition flash; temperature mode shows the zone
 *    temperatures in inferno false colour (legend: `temperatureLegend()`);
 *  - knock: (1,0) acoustic mode + expanding pressure front overlay;
 *  - spark: breakdown flash, arc vs glow column bent by the flow, glare;
 *  - flow tracers through the valves (rate ∝ mass flow) in simulated time;
 *  - a PointLight carrying the chamber's luminous intensity (lights the walls).
 *
 * Brightness: physical luminance × VISUAL_GAIN (constants.ts) — the one
 * labelled exposure constant.
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { TRACER_CAPACITY, VIS_TEMPERATURE_RANGE } from './constants';
import { CylinderSnapshotView } from './cylinder-view';
import { SparkView } from './spark-view';
import { CombustionVisualState, type CombustionMode } from './state';
import { TracerView } from './tracer-view';
import { TracerSystem } from './tracers';
import { GasVolume, type CutPlanes } from './volume';

export { temperatureLegend, temperatureColor, inferno, type TemperatureLegend } from './colour/colormap';
export { createTemperatureLegendElement, type TemperatureLegendElement } from './legend';
export { VISUAL_GAIN, VIS_TEMPERATURE_RANGE } from './constants';
export type { CombustionMode } from './state';
export type { CutPlanes } from './volume';
export { chamberShapeOf, chamberDepth, type ChamberShape } from './chamber';
export { CylinderSnapshotView } from './cylinder-view';

/** Per-instance options of CombustionVisuals. */
export interface CombustionVisualsOptions {
  /** 0-based cylinder this instance shows (reads `s.cylinders[cylinder]` when present). Default 0. */
  cylinder?: number;
  /** Run and draw the flow tracers (default true; the dominant CPU cost — enable on one cylinder). */
  tracers?: boolean;
  /** Add the chamber PointLight (default true; each light costs every lit material a loop iteration). */
  light?: boolean;
}

export class CombustionVisuals {
  readonly root: THREE.Group;
  /** Derived per-frame state (read-only for callers; handy for HUDs/debugging). */
  readonly state: CombustionVisualState;
  /** Flow tracers (capacity 0 and never stepped when the instance was built with `tracers: false`). */
  readonly tracers: TracerSystem;
  /** 0-based cylinder this instance shows. */
  readonly cylinder: number;

  private readonly volume: GasVolume;
  private readonly sparkView: SparkView;
  private readonly tracerView: TracerView;
  private readonly light: THREE.PointLight;
  private readonly cylView = new CylinderSnapshotView();
  private readonly tracersOn: boolean;
  private mode: CombustionMode = 'physical';
  private lightEnabled: boolean;

  /**
   * @param spec engine spec
   * @param opts.cylinder 0-based cylinder this instance shows (multi-cylinder engines read
   *   `s.cylinders[cylinder]`, falling back to the top-level fields; default 0)
   * @param opts.tracers  run the flow tracers (default true; multi-cylinder apps enable them on one
   *   featured cylinder only)
   * @param opts.light    add the chamber PointLight (default true)
   */
  constructor(spec: EngineSpec, readonly opts: CombustionVisualsOptions = {}) {
    this.root = new THREE.Group();
    this.root.name = 'combustion-visuals';
    this.cylinder = opts.cylinder ?? 0;
    this.tracersOn = opts.tracers !== false;
    this.lightEnabled = opts.light !== false;
    this.state = new CombustionVisualState(spec);
    this.tracers = new TracerSystem(spec, this.tracersOn ? TRACER_CAPACITY : 0);

    this.volume = new GasVolume(spec, this.root);
    this.sparkView = new SparkView(spec);
    this.tracerView = new TracerView(this.tracers);
    // Inverse-square light from the chamber's emitters. Windowed at ~1 bore so
    // it does not visibly leak onto the engine's exterior.
    this.light = new THREE.PointLight(0xffffff, 0, spec.geometry.bore, 2);
    this.light.name = 'combustion-chamber-light';
    this.light.castShadow = false;

    this.root.add(...this.volume.meshes);
    if (this.tracersOn) this.root.add(this.tracerView.points);
    this.root.add(this.sparkView.group);
    if (opts.light !== false) this.root.add(this.light);
  }

  /**
   * @param snapshot  snapshot for this frame (interpolated); multi-cylinder: the engine snapshot
   * @param dtWall    wall seconds since the previous frame
   * @param timeScale simulated seconds per wall second
   */
  update(snapshot: EngineSnapshot, dtWall: number, timeScale: number): void {
    const s = this.cylView.select(snapshot, this.cylinder);
    const st = this.state;
    st.update(s, dtWall, timeScale);
    if (this.tracersOn) {
      // A jump of more than one crank revolution (seek, long stall) cannot be
      // integrated meaningfully from two end points: restart the tracers.
      const rev = 60 / Math.max(Math.abs(s.rpm), 1);
      if (st.wasReset || st.dtSim > rev) this.tracers.reset();
      const tracerAllBurned = st.allBurned && s.phase !== 'gas-exchange';
      this.tracers.update(s, st.dtSim > rev ? 0 : st.dtSim, st.flow, st.flameCenter, st.flameRadius, tracerAllBurned);
    }

    this.volume.sync(st, s);
    this.sparkView.sync(st);
    if (this.tracersOn) this.tracerView.sync(st.Tb);

    const L = st.light;
    const on = this.lightEnabled && this.mode === 'physical' && L.intensity > 0;
    this.light.intensity = on ? L.intensity : 0;
    if (on) {
      this.light.color.setRGB(L.color[0], L.color[1], L.color[2], THREE.LinearSRGBColorSpace);
      this.light.position.set(L.position[0], L.position[1], L.position[2]);
    }
  }

  setMode(mode: 'physical' | 'temperature'): void {
    this.mode = mode;
    this.volume.setMode(mode);
  }

  /** False-colour range for 'temperature' mode, K (keep the UI legend in sync). */
  setTemperatureRange(min = VIS_TEMPERATURE_RANGE[0], max = VIS_TEMPERATURE_RANGE[1]): void {
    this.volume.setTemperatureRange(min, max);
  }

  setTracersVisible(on: boolean): void {
    this.tracerView.setVisible(on);
  }

  /** The chamber PointLight (flame/spark lighting the walls). On by default in 'physical' mode. */
  setChamberLightEnabled(on: boolean): void {
    this.lightEnabled = on && this.opts.light !== false;
  }

  /**
   * The engine's cut-away region in this cylinder's frame — {p : n·p > d} for every plane [nx, ny, nz, d]
   * (≤ 2) — or null (default). Only non-convex (L-head) chambers use it: gas behind a metal gap between
   * the bore column and the pocket is drawn only where that metal is cut away (never through the piston).
   */
  setCutRegion(planes: CutPlanes | null): void {
    this.volume.setCutRegion(planes);
  }

  /** Draw simple electrode tips (off by default; the engine model owns the plug). */
  setElectrodesVisible(on: boolean): void {
    this.sparkView.setElectrodesVisible(on);
  }

  dispose(): void {
    this.volume.dispose();
    this.sparkView.dispose();
    this.tracerView.dispose();
    this.light.dispose();
    this.root.removeFromParent();
  }
}
