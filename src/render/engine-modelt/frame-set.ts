/**
 * A StaticSet whose local frame is ROOT translated by `origin` (the cut frame). The section shader
 * assumes the cut planes pass through the object's local origin with local x / z normals
 * (render/engine/materials.ts), so solids are authored in world coordinates and shifted into the frame
 * here; the built groups are parented to a Group placed at `origin`.
 */
import * as THREE from 'three';
import type { CutMode } from '../engine/geometry';
import type { EngineMaterials, MaterialKey } from '../engine/materials';
import { StaticSet, type StaticPartOptions } from '../engine/static-set';

export class FrameSet {
  private readonly set = new StaticSet();
  constructor(readonly origin: [number, number, number]) {}

  /** Add a solid given in WORLD coordinates. */
  add(key: MaterialKey, geom: THREE.BufferGeometry, o?: StaticPartOptions): void {
    geom.translate(-this.origin[0], -this.origin[1], -this.origin[2]);
    this.set.add(key, geom, o);
  }

  /** Build full/cut variants inside a group placed at `origin`. */
  build(mats: EngineMaterials, name: string, mode: CutMode): { group: THREE.Group; full: THREE.Group; cut: THREE.Group } {
    const group = new THREE.Group();
    group.name = `${name}-frame`;
    group.position.set(...this.origin);
    const { full, cut } = this.set.build(mats, name, mode);
    group.add(full, cut);
    return { group, full, cut };
  }
}
