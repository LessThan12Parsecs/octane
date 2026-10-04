/**
 * One engine's worth of front end: the 3D mechanism (EngineRenderModel from the render registry), one
 * in-cylinder CombustionVisuals per cylinder frame, the simulation client (its own worker), the
 * Conductor and the UI. The App keeps the Stage, overlay, animation loop and keyboard alive and swaps
 * whole sessions when the user picks another engine (App.setEngine).
 *
 * Featured cylinder (multi-cylinder engines): only its visuals run the flow tracers and carry the
 * chamber light; setFeaturedCylinder moves both (CylinderVisuals, cylinder-visuals.ts). The focus
 * cylinder (setFocusCylinder) also carries the cutaway section of engines that can move it (the Model T
 * opens the focus cylinder's chamber), so the featured visuals are never hidden inside closed metal.
 *
 * Construction is all-or-nothing: if any part throws (no render model registered for the engine yet,
 * no WebGL, …) everything already built is released and the stage is left as it was.
 */
import * as THREE from 'three';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineDefinition } from '../physics/engines/index';
import { createEngineModel, type CombustionVisuals, type EngineRenderModel, type Stage } from '../render/index';
import { UIController } from '../ui/index';
import type { SimulatorOptions } from '../worker/protocol';
import { Conductor, type ConductorHooks, type PlaybackState } from './conductor';
import { CylinderVisuals } from './cylinder-visuals';
import { SimClient, type SimClientConfig } from './sim-client';
import type { Framing } from './sync';

export interface EngineSessionOptions {
  stage: Stage;
  /** Container of the UI overlay (controls, charts, HUD). */
  uiContainer: HTMLElement;
  def: EngineDefinition;
  operatingPoint: OperatingPoint;
  simulatorOptions: SimulatorOptions;
  playback: PlaybackState;
  /** Focus cylinder (0-based): the cutaway section and the flow tracers start there; so does the UI. */
  focusCylinder: number;
  hooks: ConductorHooks;
  /** The UI's Reset button (after the conductor reset the simulation). */
  onReset: () => void;
  /** The user picked a focus cylinder in the UI. */
  onFocusCylinder: (index: number) => void;
  simClientConfig?: SimClientConfig;
}

export class EngineSession {
  readonly def: EngineDefinition;
  readonly engine: EngineRenderModel;
  readonly sim: SimClient;
  readonly conductor: Conductor;
  readonly ui: UIController;
  private readonly visuals: CylinderVisuals;
  private disposed = false;

  constructor(o: EngineSessionOptions) {
    const def = (this.def = o.def);
    const spec = def.spec;
    const op = o.operatingPoint;
    const created: { dispose(): void }[] = [];
    const track = <T extends { dispose(): void }>(x: T): T => {
      created.push(x);
      return x;
    };
    let root: THREE.Object3D | null = null;
    try {
      // ---- 3D ----
      const engine = (this.engine = track(createEngineModel(spec, op)));
      root = engine.root;
      o.stage.scene.add(engine.root);
      engine.setControls?.(op);
      engine.setSectionCylinder?.(o.focusCylinder); // before the visuals and the conductor read the cut regions
      const visuals = (this.visuals = track(new CylinderVisuals(spec, engine, (r) => o.stage.addEmitters(r), o.focusCylinder)));
      engine.root.updateMatrixWorld(true);

      // ---- simulation ----
      this.sim = track(new SimClient(spec, op, o.simulatorOptions, o.simClientConfig));
      this.conductor = new Conductor(spec, this.sim, engine, visuals.port, o.hooks, o.playback);

      // ---- UI ----
      this.ui = track(
        new UIController(o.uiContainer, {
          spec,
          engine: def,
          initialOperatingPoint: op,
          focusCylinder: o.focusCylinder,
          onOperatingPointChange: this.conductor.handleOperatingPoint,
          onPlaybackChange: this.conductor.handlePlayback,
          onStep: this.conductor.handleStep,
          onViewChange: this.conductor.handleView,
          onReset: () => {
            this.conductor.handleReset();
            o.onReset();
          },
          onFocusCylinderChange: o.onFocusCylinder,
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
      root?.removeFromParent();
      throw err;
    }
  }

  get spec() {
    return this.def.spec;
  }

  /** One in-cylinder visual per cylinder frame (index = cylinder − 1; replaced by setFeaturedCylinder). */
  get gas(): readonly CombustionVisuals[] {
    return this.visuals.visuals;
  }

  /** Cylinder (0-based) whose visuals run the flow tracers and carry the chamber light. */
  get featuredCylinder(): number {
    return this.visuals.featured;
  }

  /**
   * Feature another cylinder (clamped): its visuals take over the flow tracers and the chamber light;
   * only the previous and the new featured cylinder's visuals are rebuilt (render mode and cut region
   * carried over). Throws (leaving everything as it was) if a replacement cannot be built.
   */
  setFeaturedCylinder(index: number): void {
    if (this.disposed) return;
    if (this.visuals.setFeatured(index, this.conductor.view.mode)) this.engine.root.updateMatrixWorld(true);
  }

  /**
   * Focus cylinder `index` in 3D: move the cutaway section there (engines whose section can move; the
   * conductor re-pushes the cut regions), then feature it (setFeaturedCylinder). Throws if a step fails:
   * a failed section move changes nothing; a failed feature leaves the section moved and the previous
   * cylinder featured.
   */
  setFocusCylinder(index: number): void {
    if (this.disposed) return;
    this.conductor.setSectionCylinder(index);
    this.setFeaturedCylinder(index);
  }

  /** World bounds of the mechanism (for Stage.fitToBounds). */
  bounds(): THREE.Box3 {
    this.engine.root.updateMatrixWorld(true);
    return new THREE.Box3().setFromObject(this.engine.root);
  }

  /** Recommended camera for a framing ('chamber' frames cylinder `focus`). */
  cameraView(framing: Framing, focus: number) {
    return this.engine.cameraView(framing, focus);
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    this.ui.dispose();
    this.sim.dispose();
    this.visuals.dispose();
    this.engine.root.removeFromParent();
    this.engine.dispose();
  }
}
