/**
 * One engine's worth of front end: the 3D mechanism (EngineRenderModel from the render registry), one
 * in-cylinder CombustionVisuals per cylinder frame, the simulation client (its own worker), the
 * Conductor and the UI. The App keeps the Stage, overlay, animation loop and keyboard alive and swaps
 * whole sessions when the user picks another engine (App.setEngine).
 *
 * Construction is all-or-nothing: if any part throws (no render model registered for the engine yet,
 * no WebGL, …) everything already built is released and the stage is left as it was.
 */
import * as THREE from 'three';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineDefinition } from '../physics/engines/index';
import { CombustionVisuals, createEngineModel, type EngineRenderModel, type Stage } from '../render/index';
import { UIController } from '../ui/index';
import type { SimulatorOptions } from '../worker/protocol';
import { Conductor, GasFanOut, type ConductorHooks, type PlaybackState } from './conductor';
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
  /** Featured cylinder (0-based): its visuals run the flow tracers; the UI starts focused on it. */
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
  /** One per cylinder frame (index = cylinder − 1). */
  readonly gas: readonly CombustionVisuals[];
  readonly sim: SimClient;
  readonly conductor: Conductor;
  readonly ui: UIController;
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
      const frames = engine.cylinderFrames;
      const gas: CombustionVisuals[] = [];
      for (let i = 0; i < frames.length; i++) {
        const g = track(new CombustionVisuals(spec, { cylinder: i, tracers: i === o.focusCylinder || frames.length === 1 }));
        frames[i].add(g.root);
        o.stage.addEmitters(g.root); // flame, spark, tracers: the only bloom sources
        gas.push(g);
      }
      this.gas = gas;
      engine.root.updateMatrixWorld(true);

      // ---- simulation ----
      this.sim = track(new SimClient(spec, op, o.simulatorOptions, o.simClientConfig));
      this.conductor = new Conductor(spec, this.sim, engine, gas.length === 1 ? gas[0] : new GasFanOut(gas), o.hooks, o.playback);

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
    for (const g of this.gas) {
      g.root.removeFromParent();
      g.dispose();
    }
    this.engine.root.removeFromParent();
    this.engine.dispose();
  }
}
