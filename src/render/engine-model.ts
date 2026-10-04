/**
 * Engine render-model contract and registry. Every engine's 3D mechanism implements EngineRenderModel;
 * the app builds it with createEngineModel(spec, op) and never touches a concrete class.
 *
 * Frames: `root` is the world/ROOT frame (crank axis = world z through the origin, +y up, +z = front).
 * `cylinderFrames[i]` is cylinder i's frame (engine-spec.ts: the physics cylinder frame translated to
 * the cylinder's position, mirrored z → −z where spec.layout.mirrorZ[i]); in-cylinder visuals
 * (CombustionVisuals) are parented there. The model is animated ONLY from snapshot fields: per-cylinder
 * kinematics from `s.cylinders[i]` when present, else derived from θ and the spec (identical cylinders).
 */
import type * as THREE from 'three';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineSnapshot } from '../physics/core/snapshot';
import type { CutPlanes } from './combustion/index';
import { EngineModel, recommendedCameraView, type CameraView } from './engine/index';

export type { CameraView } from './engine/index';

export interface EngineRenderModel {
  readonly root: THREE.Group;
  /** One frame per cylinder (index = cylinder − 1). */
  readonly cylinderFrames: readonly THREE.Object3D[];
  /** Compression ratio the mechanism shows (variable-CR engines move the head). */
  readonly compressionRatio: number;
  readonly isCutaway: boolean;
  update(s: EngineSnapshot): void;
  /**
   * Operating-point inputs that move parts but are not physics state (spark lever → timer case,
   * hand throttle → butterfly). Optional; called whenever the operating point changes.
   */
  setControls?(op: Partial<OperatingPoint>): void;
  /** No-op for fixed-CR engines. */
  setCompressionRatio(cr: number): void;
  setCutaway(on: boolean): void;
  /**
   * The cut-away region of the housings in cylinder `cylinder`'s frame (0-based), for the in-cylinder
   * visuals (CombustionVisuals.setCutRegion): {p : n·p > d for every plane [nx, ny, nz, d]}, mirrored
   * frames included; null while the cutaway is off. Optional: engines whose chamber is convex (the CFR
   * disc) need not provide it, and the app then never sets a cut region.
   */
  cutRegion?(cylinder: number): CutPlanes | null;
  /**
   * Multi-cylinder cutaways: move the housings' section to cylinder `cylinder` (0-based, clamped) so its
   * chamber is open; cutRegion() then reports the moved region for every cylinder, and the 'chamber'
   * camera view of that cylinder looks through the section. Optional: single-cylinder engines (the CFR)
   * have one fixed section. The app calls it with the focus cylinder and re-pushes the cut regions.
   */
  setSectionCylinder?(cylinder: number): void;
  /** Camera framing: 'engine' = whole machine, 'chamber' = close-up of cylinder `cylinder` (0-based, default 0). */
  cameraView(framing: 'engine' | 'chamber', cylinder?: number): CameraView;
  /** Vertical shift of the chamber when the CR changes from `previous` to `cr`, m (0 for fixed-CR engines). */
  chamberShift(cr: number, previous: number): number;
  dispose(): void;
}

/** The CFR F-1 EngineModel behind the EngineRenderModel contract. */
export class CfrEngineRenderModel implements EngineRenderModel {
  readonly model: EngineModel;
  readonly cylinderFrames: readonly THREE.Object3D[];

  constructor(private readonly spec: EngineSpec, compressionRatio: number) {
    this.model = new EngineModel(spec, compressionRatio);
    this.cylinderFrames = [this.model.cylinderFrame];
  }

  get root(): THREE.Group {
    return this.model.root;
  }
  get compressionRatio(): number {
    return this.model.compressionRatio;
  }
  get isCutaway(): boolean {
    return this.model.isCutaway;
  }
  update(s: EngineSnapshot): void {
    this.model.update(s);
  }
  setCompressionRatio(cr: number): void {
    this.model.setCompressionRatio(cr);
  }
  setCutaway(on: boolean): void {
    this.model.setCutaway(on);
  }
  cameraView(framing: 'engine' | 'chamber'): CameraView {
    return recommendedCameraView(this.spec, this.model.compressionRatio, framing);
  }
  chamberShift(cr: number, previous: number): number {
    return this.model.layout.headY(cr) - this.model.layout.headY(previous);
  }
  dispose(): void {
    this.model.dispose();
  }
}

/**
 * Builders of non-CFR engines, keyed by EngineSpec.id (registered by their modules: the Ford Model T
 * ('ford-model-t') by render/engine-modelt/index.ts, which the render barrel render/index.ts imports for
 * that side effect — import createEngineModel from the barrel, or import the engine module yourself).
 */
const BUILDERS = new Map<string, (spec: EngineSpec, op: OperatingPoint) => EngineRenderModel>();

/** Register the render-model builder of an engine id (called by e.g. render/engine-modelt). */
export function registerEngineModel(id: string, build: (spec: EngineSpec, op: OperatingPoint) => EngineRenderModel): void {
  BUILDERS.set(id, build);
}

/**
 * Build the render model of a spec: the registered builder for spec.id, else the CFR model for a
 * single-cylinder flat-disc spec. Throws for any other spec rather than rendering a broken CFR.
 */
export function createEngineModel(spec: EngineSpec, op: OperatingPoint): EngineRenderModel {
  const b = spec.id ? BUILDERS.get(spec.id) : undefined;
  if (b) return b(spec, op);
  if (spec.cylinders === 1 && spec.geometry.chamber === 'flat-disc') return new CfrEngineRenderModel(spec, op.compressionRatio);
  throw new Error(`No render model registered for engine '${spec.id ?? spec.name}'`);
}
