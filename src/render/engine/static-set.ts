/**
 * Collects the static (non-moving relative to their parent frame) solids of
 * one frame, merges them per material and builds both display variants:
 *  - full: plain meshes;
 *  - cut:  geometry clipped to z ≤ 0 + back-face section meshes.
 */
import * as THREE from 'three';
import { type CutMode, clipForCut, merge } from './geometry';
import type { EngineMaterials, MaterialKey } from './materials';

export interface StaticPartOptions {
  /** Clip in cutaway mode (default true). */
  cut?: boolean;
  /**
   * Draw section caps for this solid in cutaway mode ('dark' = iron, 'light' = non-ferrous).
   * Only for CLOSED solids. Default 'dark'. `false` for open surfaces.
   */
  section?: 'dark' | 'light' | false;
}

interface Entry {
  key: MaterialKey;
  geom: THREE.BufferGeometry;
  cut: boolean;
  section: 'dark' | 'light' | false;
}

export class StaticSet {
  private entries: Entry[] = [];

  add(key: MaterialKey, geom: THREE.BufferGeometry, o: StaticPartOptions = {}): void {
    if (geom.getAttribute('position').count === 0) return;
    this.entries.push({ key, geom, cut: o.cut ?? true, section: o.section ?? 'dark' });
  }

  /** Build `full` and `cut` groups; disposes the source geometries. */
  build(mats: EngineMaterials, name: string, mode: CutMode = { kind: 'half' }): { full: THREE.Group; cut: THREE.Group } {
    const full = new THREE.Group();
    full.name = `${name}-full`;
    const cut = new THREE.Group();
    cut.name = `${name}-cut`;
    const byKey = new Map<MaterialKey, { full: THREE.BufferGeometry[]; cut: THREE.BufferGeometry[] }>();
    const sections: Record<'dark' | 'light', THREE.BufferGeometry[]> = { dark: [], light: [] };
    for (const e of this.entries) {
      let slot = byKey.get(e.key);
      if (!slot) byKey.set(e.key, (slot = { full: [], cut: [] }));
      slot.full.push(e.geom);
      if (e.cut) {
        const c = clipForCut(e.geom, mode);
        if (c.getAttribute('position').count > 0) {
          slot.cut.push(c);
          if (e.section) sections[e.section].push(c);
        }
      } else {
        slot.cut.push(e.geom);
      }
    }
    const disposeLater = new Set<THREE.BufferGeometry>();
    for (const [key, slot] of byKey) {
      const mat = mats.m[key];
      if (slot.full.length) {
        const mesh = new THREE.Mesh(merge(slot.full), mat);
        mesh.name = `${name}-${key}`;
        full.add(mesh);
      }
      if (slot.cut.length) {
        const mesh = new THREE.Mesh(merge(slot.cut), mat);
        mesh.name = `${name}-${key}-cut`;
        cut.add(mesh);
      }
      for (const g of slot.full) disposeLater.add(g);
      for (const g of slot.cut) disposeLater.add(g);
    }
    for (const tint of ['dark', 'light'] as const) {
      if (!sections[tint].length) continue;
      const mesh = new THREE.Mesh(merge(sections[tint]), mats.sectionMaterial(tint, mode));
      mesh.name = `${name}-section-${tint}`;
      mesh.renderOrder = 1;
      cut.add(mesh);
    }
    for (const g of disposeLater) g.dispose();
    this.entries = [];
    return { full, cut };
  }
}
