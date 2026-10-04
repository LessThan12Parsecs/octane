/**
 * The in-cylinder visuals of one engine session: one CombustionVisuals per cylinder frame of the engine
 * render model (parented to that frame, so mirrored frames work), each registered as a bloom emitter.
 *
 * Featured cylinder: only its visuals run the flow tracers (their dominant CPU cost) and carry the
 * chamber PointLight, so the scene always holds exactly ONE combustion light — a change of the light
 * count would make three.js recompile every lit material. Tracer capacity and the light are fixed when a
 * CombustionVisuals is built (render/combustion/index.ts), so setFeatured rebuilds just the two instances
 * whose flags change (the previous and the new featured cylinder) and carries the view state the app
 * sets on them — render mode and cut region — over; the other cylinders are untouched. (The app never
 * changes the temperature range, tracer visibility, chamber-light switch or electrode tips from their
 * defaults, so a rebuilt instance needs nothing else.) A single-cylinder engine always features its one
 * cylinder, exactly as before the multi-cylinder work.
 */
import type * as THREE from 'three';
import type { EngineSpec } from '../physics/core/engine-spec';
import { CombustionVisuals, type CutPlanes } from '../render/combustion/index';
import { GasFanOut, type GasPort, type RenderMode } from './conductor';

/** What the visuals need from the engine render model (EngineRenderModel satisfies it). */
export interface CylinderFrameHost {
  readonly cylinderFrames: readonly THREE.Object3D[];
  cutRegion?(cylinder: number): CutPlanes | null;
}

export class CylinderVisuals {
  /** The conductor's gas port: the single visual, or a fan-out over all of them. */
  readonly port: GasPort;
  private readonly list: CombustionVisuals[] = [];
  private featuredIndex: number;
  private disposed = false;

  /**
   * @param addEmitters registers a built subtree as bloom emitters (Stage.addEmitters)
   * @param featured    0-based featured cylinder (clamped)
   */
  constructor(
    private readonly spec: EngineSpec,
    private readonly host: CylinderFrameHost,
    private readonly addEmitters: (root: THREE.Object3D) => void,
    featured = 0,
  ) {
    const n = host.cylinderFrames.length;
    this.featuredIndex = clampIndex(featured, n);
    try {
      for (let i = 0; i < n; i++) this.list.push(this.build(i, this.featuredIndex));
    } catch (err) {
      for (const g of this.list) g.dispose();
      throw err;
    }
    // A lone visual is the port itself (the CFR, exactly as before); the fan-out also routes cut regions.
    this.port = this.list.length === 1 && !host.cutRegion ? this.list[0] : new GasFanOut(this.list);
  }

  /** One visual per cylinder frame (index = cylinder − 1; entries are replaced by setFeatured). */
  get visuals(): readonly CombustionVisuals[] {
    return this.list;
  }

  /** 0-based cylinder whose visuals run the flow tracers and carry the chamber light. */
  get featured(): number {
    return this.featuredIndex;
  }

  /**
   * Feature cylinder `index` (clamped): rebuild the previous and the new featured cylinder's visuals with
   * swapped tracer/light flags, apply `mode` and the host's current cut region to them, and swap them into
   * the gas port. Both replacements are built before anything is swapped, so a failure leaves everything
   * as it was. Returns false when nothing changed.
   */
  setFeatured(index: number, mode: RenderMode): boolean {
    if (this.disposed) return false;
    const i = clampIndex(index, this.list.length);
    const prev = this.featuredIndex;
    if (i === prev) return false;
    const a = this.build(prev, i);
    let b: CombustionVisuals;
    try {
      b = this.build(i, i);
    } catch (err) {
      a.dispose();
      throw err;
    }
    this.featuredIndex = i;
    for (const g of [a, b]) {
      const k = g.cylinder;
      g.setMode(mode);
      if (this.host.cutRegion) g.setCutRegion(this.host.cutRegion(k));
      const old = this.list[k];
      this.list[k] = g;
      if (this.port instanceof GasFanOut) this.port.replace(k, g);
      old.dispose(); // also detaches it from the cylinder frame
    }
    return true;
  }

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    for (const g of this.list) g.dispose();
  }

  /** Visuals of cylinder `i` with the flags for featured cylinder `featured`, attached to its frame. */
  private build(i: number, featured: number): CombustionVisuals {
    const frames = this.host.cylinderFrames;
    const on = i === featured || frames.length === 1;
    const g = new CombustionVisuals(this.spec, { cylinder: i, tracers: on, light: on });
    try {
      frames[i].add(g.root);
      this.addEmitters(g.root); // flame, spark, tracers: the only bloom sources
    } catch (err) {
      g.dispose();
      throw err;
    }
    return g;
  }
}

function clampIndex(i: number, n: number): number {
  return Math.max(0, Math.min(n - 1, Math.floor(i) || 0));
}
