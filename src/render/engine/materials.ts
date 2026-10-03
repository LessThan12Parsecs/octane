/**
 * PBR materials for the engine model (tuned for a RoomEnvironment env map
 * plus a key light). Colours are sRGB hex.
 */
import * as THREE from 'three';
import type { CutMode } from './geometry';

export type MaterialKey =
  | 'paint' // painted cast iron (crankcase, cylinder, head exterior)
  | 'castRaw' // as-cast iron interior (crankcase cavity, water jacket)
  | 'ironMachined' // bores, deck faces, bearing bores
  | 'steelPolished' // journals, pins, valve stems, pushrods
  | 'steelForged' // crank webs, connecting rod
  | 'steelDark' // piston rings, bolts, springs
  | 'aluminium' // piston
  | 'aluminiumCast' // rocker cover
  | 'bearing' // babbitt/bronze shells, bushings
  | 'bronze' // worm wheel
  | 'valve' // heat-resistant valve steel
  | 'valveExhaust'
  | 'ceramic' // spark-plug insulator
  | 'plated' // plug shell (zinc/nickel plated)
  | 'nickel' // electrodes
  | 'copper' // gaskets, washers
  | 'rubber' // belts, plug boot
  | 'motor' // motor housing paint
  | 'base' // bedplate
  | 'pistonIron'; // cast-iron piston

export class EngineMaterials {
  readonly m: Record<MaterialKey, THREE.MeshStandardMaterial>;
  private readonly sections = new Map<string, THREE.MeshStandardMaterial>();

  constructor() {
    const std = (color: number, metalness: number, roughness: number, extra: THREE.MeshStandardMaterialParameters = {}) =>
      new THREE.MeshStandardMaterial({ color, metalness, roughness, ...extra });
    this.m = {
      paint: std(0x3b4c58, 0.35, 0.52, { envMapIntensity: 0.9 }),
      castRaw: std(0x5b5752, 0.55, 0.82),
      ironMachined: std(0x9ea3a8, 0.9, 0.3),
      steelPolished: std(0xd3d7dc, 1.0, 0.16),
      steelForged: std(0x80868d, 0.95, 0.4),
      steelDark: std(0x3a3d41, 0.9, 0.34),
      aluminium: std(0xc4c8cc, 1.0, 0.36),
      aluminiumCast: std(0xa9adb0, 0.85, 0.55),
      bearing: std(0xc49a5a, 1.0, 0.3),
      bronze: std(0xa8703a, 1.0, 0.34),
      valve: std(0x9a968f, 1.0, 0.3),
      valveExhaust: std(0x7d7063, 0.95, 0.42),
      ceramic: new THREE.MeshPhysicalMaterial({ color: 0xf3f0e8, metalness: 0, roughness: 0.2, clearcoat: 0.7, clearcoatRoughness: 0.15 }),
      plated: std(0xd6d8d2, 1.0, 0.22),
      nickel: std(0xb8b2a6, 1.0, 0.28),
      copper: std(0xb8733a, 1.0, 0.28),
      rubber: std(0x1c1c1e, 0.0, 0.78),
      motor: std(0x2f4a3e, 0.3, 0.5),
      base: std(0x24292d, 0.4, 0.6),
      pistonIron: std(0x8f949a, 0.85, 0.34),
    };
  }

  /** Back-face "section cap" material for a cut mode ('dark' = ferrous, 'light' = non-ferrous tint). */
  sectionMaterial(tint: 'dark' | 'light', mode: CutMode): THREE.MeshStandardMaterial {
    const key = `${tint}:${mode.kind}:${mode.kind === 'quadrant' ? mode.side : 0}`;
    let m = this.sections.get(key);
    if (!m) {
      m = makeSectionMaterial(tint === 'dark' ? 0xb8472f : 0xc8763a, mode);
      this.sections.set(key, m);
    }
    return m;
  }

  dispose(): void {
    for (const k of Object.keys(this.m) as MaterialKey[]) this.m[k].dispose();
    for (const m of this.sections.values()) m.dispose();
    this.sections.clear();
  }
}

/**
 * Section material: renders BACK faces of a clipped closed solid as if they
 * were the planar cut face. Per fragment it finds which cut face the view ray
 * entered through (z = 0, or x = 0 for a quarter section), shades with that
 * face's normal and draws a faint engineering hatch in the cut plane.
 * Negative polygon offset so it wins depth ties against coincident front
 * faces of neighbouring solids (joints). Assumes the object's local frame is
 * the cut frame (the static assemblies are only translated).
 */
export function makeSectionMaterial(color: number, mode: CutMode): THREE.MeshStandardMaterial {
  const m = new THREE.MeshStandardMaterial({
    color,
    metalness: 0.15,
    roughness: 0.62,
    side: THREE.BackSide,
    polygonOffset: true,
    polygonOffsetFactor: -1,
    polygonOffsetUnits: -4,
  });
  const quad = mode.kind === 'quadrant';
  const side = mode.kind === 'quadrant' ? mode.side : 1;
  const varyings = 'varying vec3 vSecW;\nvarying vec3 vSecO;\nvarying vec3 vSecX;\nvarying vec3 vSecY;\nvarying vec3 vSecZ;';
  m.onBeforeCompile = (shader) => {
    shader.vertexShader = shader.vertexShader
      .replace('#include <common>', `#include <common>\n${varyings}`)
      .replace(
        '#include <begin_vertex>',
        [
          '#include <begin_vertex>',
          'vSecW = ( modelMatrix * vec4( transformed, 1.0 ) ).xyz;',
          'vSecO = ( modelMatrix * vec4( 0.0, 0.0, 0.0, 1.0 ) ).xyz;',
          'vSecX = normalize( mat3( modelMatrix ) * vec3( 1.0, 0.0, 0.0 ) );',
          'vSecY = normalize( mat3( modelMatrix ) * vec3( 0.0, 1.0, 0.0 ) );',
          'vSecZ = normalize( mat3( modelMatrix ) * vec3( 0.0, 0.0, 1.0 ) );',
        ].join('\n'),
      );
    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', `#include <common>\n${varyings}\nvec3 secNormalW;`)
      .replace(
        '#include <color_fragment>',
        [
          '#include <color_fragment>',
          '{',
          '  vec3 C = cameraPosition;',
          '  vec3 rd = vSecW - C;',
          '  // entry through the z = 0 face',
          '  float dz = dot( vSecZ, rd );',
          '  float tz = abs( dz ) > 1e-9 ? dot( vSecZ, vSecO - C ) / dz : 2.0;',
          '  vec3 hz = C + rd * tz;',
          `  bool okZ = tz > 0.0 && tz < 1.0 ${quad ? `&& ${side.toFixed(1)} * dot( hz - vSecO, vSecX ) >= 0.0` : ''};`,
          '  vec3 hit = hz;',
          '  secNormalW = vSecZ;',
          '  vec2 uv2 = vec2( dot( hz - vSecO, vSecX ), dot( hz - vSecO, vSecY ) );',
          quad
            ? [
                '  float dx = dot( vSecX, rd );',
                '  float tx = abs( dx ) > 1e-9 ? dot( vSecX, vSecO - C ) / dx : 2.0;',
                '  vec3 hx = C + rd * tx;',
                '  bool okX = tx > 0.0 && tx < 1.0 && dot( hx - vSecO, vSecZ ) >= 0.0;',
                '  if ( okX && ( !okZ || tx < tz ) ) {',
                `    secNormalW = ${side.toFixed(1)} * vSecX;`,
                '    uv2 = vec2( dot( hx - vSecO, vSecZ ), dot( hx - vSecO, vSecY ) );',
                '  }',
              ].join('\n')
            : '',
          '  float d = ( uv2.x + uv2.y ) / 0.0035;',
          '  float f = abs( fract( d ) - 0.5 );',
          '  float w = fwidth( d );',
          '  float line = 1.0 - smoothstep( 0.1 - w, 0.1 + w, f );',
          '  diffuseColor.rgb *= 1.0 - 0.18 * line;',
          '}',
        ].join('\n'),
      )
      .replace(
        '#include <normal_fragment_begin>',
        '#include <normal_fragment_begin>\nnormal = normalize( ( viewMatrix * vec4( secNormalW, 0.0 ) ).xyz );',
      );
  };
  m.customProgramCacheKey = () => `octane-section-v3-${mode.kind}-${side}`;
  return m;
}
