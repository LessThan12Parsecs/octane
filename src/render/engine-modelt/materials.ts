/**
 * Model T palette over the shared EngineMaterials (same keys, so StaticSet and the section caps work
 * unchanged), plus a few extra materials for parts outside the static sets. Each EngineMaterials
 * instance owns its materials, so re-tinting here never touches the CFR model.
 *
 * UNVERIFIED (render only): black enamel engine finish, natural cast-iron exhaust manifold, brass
 * carburettor and timer fittings, leather/fabric fan belt, cast-iron pistons.
 */
import * as THREE from 'three';
import { EngineMaterials } from '../engine/materials';

export class ModelTMaterials extends EngineMaterials {
  /** Materials not in the shared MaterialKey union (used by plain meshes only). */
  readonly x: {
    brass: THREE.MeshStandardMaterial;
    leather: THREE.MeshStandardMaterial;
    magnet: THREE.MeshStandardMaterial;
    coil: THREE.MeshStandardMaterial;
    fibre: THREE.MeshStandardMaterial;
    exhaustIron: THREE.MeshStandardMaterial;
  };

  constructor() {
    super();
    const m = this.m;
    m.paint.color.setHex(0x26282b);
    m.paint.metalness = 0.25;
    m.paint.roughness = 0.4;
    m.paint.envMapIntensity = 1.0;
    m.castRaw.color.setHex(0x55504a);
    m.motor.color.setHex(0x2b2d30);
    m.base.color.setHex(0x1f2124);
    const std = (color: number, metalness: number, roughness: number) => new THREE.MeshStandardMaterial({ color, metalness, roughness });
    this.x = {
      brass: std(0xc9a14a, 1.0, 0.3),
      leather: std(0x5a3b22, 0.0, 0.72),
      magnet: std(0x2e3135, 0.9, 0.38),
      coil: std(0xb06a34, 0.95, 0.4),
      fibre: std(0x3a2a20, 0.0, 0.6),
      exhaustIron: std(0x5e4a3f, 0.6, 0.78),
    };
  }

  override dispose(): void {
    super.dispose();
    for (const k of Object.keys(this.x) as (keyof ModelTMaterials['x'])[]) this.x[k].dispose();
  }
}
