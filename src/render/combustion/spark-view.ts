/**
 * Spark discharge rendering: the column between the electrodes (bent by the
 * gas flow along a quadratic Bézier), a glare halo, and the breakdown flash.
 * Everything additive, driven by CombustionVisualState.spark.
 *
 * Column shading: an optically thin uniform cylinder seen from the side has a
 * line-of-sight path ∝ sqrt(1 − (b/a)²) = |n·v| on its front surface, so
 * the core brightness is uCore·|n·v|^k (k = 1 thin arc, k < 1 diffuse glow).
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { CombustionVisualState } from './state';

const CHANNEL_VERTEX = /* glsl */ `
uniform vec3 uA;
uniform vec3 uB;
uniform vec3 uC;     // Bézier control point
uniform float uRadius;
varying vec3 vNormalV;
varying vec3 vViewPos;
void main() {
  float t = position.y + 0.5;           // CylinderGeometry y ∈ [-0.5, 0.5]
  float s = 1.0 - t;
  vec3 p = s * s * uA + 2.0 * s * t * uC + t * t * uB;
  vec3 tan = normalize(2.0 * s * (uC - uA) + 2.0 * t * (uB - uC) + 1e-9);
  vec3 ref = abs(tan.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
  vec3 n1 = normalize(cross(tan, ref));
  vec3 n2 = cross(tan, n1);
  vec2 ring = normalize(position.xz + 1e-9);
  vec3 nrm = ring.x * n1 + ring.y * n2;
  // taper at the electrode roots (cathode spot / anode spot)
  float taper = 0.55 + 0.45 * sin(3.14159265 * clamp(t, 0.0, 1.0));
  vec3 world = p + nrm * uRadius * taper;
  vec4 mv = modelViewMatrix * vec4(world, 1.0);
  vViewPos = mv.xyz;
  vNormalV = normalize(normalMatrix * nrm);
  gl_Position = projectionMatrix * mv;
}
`;

const CHANNEL_FRAGMENT = /* glsl */ `
uniform vec3 uCore;
uniform float uSharp;
varying vec3 vNormalV;
varying vec3 vViewPos;
void main() {
  vec3 v = normalize(-vViewPos);
  float nv = abs(dot(normalize(vNormalV), v));
  if (isOrthographic) nv = abs(normalize(vNormalV).z);
  vec3 c = uCore * pow(nv, uSharp);
  // Hue-preserving exposure clip (like an over-exposed photo of a spark):
  // brightness saturates as 1 − e^{−m}; strongly over-exposed centres go white-hot.
  float m = max(max(c.r, c.g), c.b);
  vec3 hue = c / max(m, 1e-6);
  float white = clamp(log(max(m, 1.0)) / log(256.0), 0.0, 1.0) * pow(nv, 4.0);
  vec3 col = mix(hue, vec3(1.0), white) * (1.0 - exp(-m));
  gl_FragColor = vec4(col, 1.0);
  #include <colorspace_fragment>
}
`;

/**
 * Radial glare/flash falloff texture (no DOM). Returns the texture and the
 * mean texel value over the inscribed disc, so a sprite tinted with
 * `colour / discMean` has that mean radiance over the disc.
 */
function makeGlareTexture(size = 64): { tex: THREE.DataTexture; discMean: number } {
  const data = new Uint8Array(size * size * 4);
  const c = (size - 1) / 2;
  let sum = 0, cnt = 0;
  for (let j = 0; j < size; j++) {
    for (let i = 0; i < size; i++) {
      const r = Math.hypot(i - c, j - c) / c;
      const v = r >= 1 ? 0 : Math.pow(1 - r, 3) * 0.8 + Math.exp(-r * r * 40) * 0.2;
      const b = Math.min(255, Math.round(v * 255));
      const k = 4 * (j * size + i);
      data[k] = data[k + 1] = data[k + 2] = b;
      data[k + 3] = 255;
      if (r < 1) { sum += b / 255; cnt++; }
    }
  }
  const tex = new THREE.DataTexture(data, size, size, THREE.RGBAFormat, THREE.UnsignedByteType);
  tex.magFilter = THREE.LinearFilter;
  tex.minFilter = THREE.LinearFilter;
  tex.needsUpdate = true;
  return { tex, discMean: cnt > 0 ? sum / cnt : 1 };
}

export class SparkView {
  readonly group = new THREE.Group();
  private readonly core: THREE.Mesh<THREE.CylinderGeometry, THREE.ShaderMaterial>;
  private readonly glare: THREE.Sprite;
  private readonly flash: THREE.Sprite;
  private readonly glareTex: THREE.DataTexture;
  private readonly texScale: number;
  private readonly tmpV = new THREE.Vector3();
  private readonly camLocal = new THREE.Vector3();
  private readonly invRoot = new THREE.Matrix4();
  private readonly glareBase = new THREE.Vector3();
  private readonly flashBase = new THREE.Vector3();
  private readonly electrodes: THREE.Group;
  private readonly gapCenter: THREE.Vector3;

  constructor(spec: EngineSpec) {
    this.group.name = 'combustion-spark';
    const gc = spec.sparkPlug.gapCenter;
    this.gapCenter = new THREE.Vector3(gc[0], gc[1], gc[2]);

    const geo = new THREE.CylinderGeometry(1, 1, 1, 14, 24, true);
    const mat = new THREE.ShaderMaterial({
      name: 'SparkChannel',
      vertexShader: CHANNEL_VERTEX,
      fragmentShader: CHANNEL_FRAGMENT,
      uniforms: {
        uA: { value: new THREE.Vector3() },
        uB: { value: new THREE.Vector3() },
        uC: { value: new THREE.Vector3() },
        uRadius: { value: 5e-5 },
        uCore: { value: new THREE.Vector3() },
        uSharp: { value: 1 },
      },
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      side: THREE.FrontSide,
      // exposure clipping is done in the shader (hue-preserving)
      toneMapped: false,
    });
    this.core = new THREE.Mesh(geo, mat);
    this.core.frustumCulled = false;
    this.core.renderOrder = 12;
    this.core.visible = false;

    const gt = makeGlareTexture();
    this.glareTex = gt.tex;
    this.texScale = 1 / gt.discMean;
    const mkSprite = (): THREE.Sprite => {
      const m = new THREE.SpriteMaterial({
        map: this.glareTex,
        blending: THREE.AdditiveBlending,
        depthWrite: false,
        transparent: true,
        sizeAttenuation: true,
      });
      const sp = new THREE.Sprite(m);
      sp.renderOrder = 13;
      sp.visible = false;
      sp.position.copy(this.gapCenter);
      return sp;
    };
    this.glare = mkSprite();
    this.flash = mkSprite();
    this.flashBase.copy(this.gapCenter);
    this.glareBase.copy(this.gapCenter);
    // Glare is optical (in the camera), so pull the sprite toward the camera by
    // its radius: walls just behind the source must not clip the halo.
    this.glare.onBeforeRender = (_r, _s, cam) => this.pull(this.glare, this.glareBase, cam);
    this.flash.onBeforeRender = (_r, _s, cam) => this.pull(this.flash, this.flashBase, cam);

    this.electrodes = this.buildElectrodes(spec);
    this.electrodes.visible = false;
    this.group.add(this.core, this.glare, this.flash, this.electrodes);
  }

  private pull(sprite: THREE.Sprite, base: THREE.Vector3, camera: THREE.Camera): void {
    const parent = sprite.parent;
    if (!parent) return;
    this.invRoot.copy(parent.matrixWorld).invert();
    this.camLocal.setFromMatrixPosition(camera.matrixWorld).applyMatrix4(this.invRoot);
    const d = this.tmpV.copy(this.camLocal).sub(base);
    const L = d.length();
    const shift = Math.min(0.9 * sprite.scale.x * 0.5, 0.5 * L);
    sprite.position.copy(base).addScaledVector(d, L > 0 ? shift / L : 0);
    sprite.updateMatrixWorld();
  }

  /** Simple electrode tips (off by default: the engine model owns the plug). */
  setElectrodesVisible(on: boolean): void {
    this.electrodes.visible = on;
  }

  private buildElectrodes(spec: EngineSpec): THREE.Group {
    const g = new THREE.Group();
    g.name = 'combustion-electrodes';
    const sp = spec.sparkPlug;
    const axis = new THREE.Vector3(...sp.axis).normalize();
    const mat = new THREE.MeshStandardMaterial({ color: 0xb8b8b8, metalness: 0.9, roughness: 0.35 });
    // centre electrode: cylinder ending at gapCenter − axis·gap/2, extending back along −axis
    const len = 3 * sp.gap + 1e-3;
    const ce = new THREE.Mesh(new THREE.CylinderGeometry(sp.centerElectrodeDiameter / 2, sp.centerElectrodeDiameter / 2, len, 16), mat);
    const q = new THREE.Quaternion().setFromUnitVectors(new THREE.Vector3(0, 1, 0), axis);
    ce.quaternion.copy(q);
    ce.position.copy(this.gapCenter).addScaledVector(axis, -sp.gap / 2 - len / 2);
    // ground electrode: a bar across the gap end (J-gap)
    const w = sp.groundElectrodeWidth;
    const ge = new THREE.Mesh(new THREE.BoxGeometry(w, w * 0.5, w * 2.2), mat);
    ge.quaternion.copy(q);
    ge.position.copy(this.gapCenter).addScaledVector(axis, sp.gap / 2 + w * 0.25);
    g.add(ce, ge);
    return g;
  }

  sync(st: CombustionVisualState): void {
    const sp = st.spark;
    const u = this.core.material.uniforms;
    this.core.visible = sp.visible;
    if (sp.visible) {
      (u.uA.value as THREE.Vector3).fromArray(sp.start);
      (u.uB.value as THREE.Vector3).fromArray(sp.end);
      // Bézier control point so the curve's midpoint is displaced by `bow`
      (u.uC.value as THREE.Vector3)
        .set((sp.start[0] + sp.end[0]) / 2, (sp.start[1] + sp.end[1]) / 2, (sp.start[2] + sp.end[2]) / 2)
        .addScaledVector(this.tmpV.fromArray(sp.bow), 2);
      u.uRadius.value = sp.drawRadius;
      (u.uCore.value as THREE.Vector3).fromArray(sp.coreRadiance);
      // thin arc: sharp; glow: diffuse, broader profile
      u.uSharp.value = sp.kind === 'arc' ? 1.0 : 0.45;
    }

    const gl = sp.glare;
    const glareOn = sp.visible && gl[0] + gl[1] + gl[2] > 1e-4;
    this.glare.visible = glareOn;
    if (glareOn) {
      const m = this.glare.material;
      const k = this.texScale;
      m.color.setRGB(gl[0] * k, gl[1] * k, gl[2] * k, THREE.LinearSRGBColorSpace);
      this.glare.scale.setScalar(2 * sp.glareRadius);
      this.glareBase.set(
        (sp.start[0] + sp.end[0]) / 2 + sp.bow[0],
        (sp.start[1] + sp.end[1]) / 2 + sp.bow[1],
        (sp.start[2] + sp.end[2]) / 2 + sp.bow[2],
      );
    }
    const fl = sp.flash;
    const flashOn = fl[0] + fl[1] + fl[2] > 1e-4;
    this.flash.visible = flashOn;
    if (flashOn) {
      const k = this.texScale;
      this.flash.material.color.setRGB(fl[0] * k, fl[1] * k, fl[2] * k, THREE.LinearSRGBColorSpace);
      this.flash.scale.setScalar(2 * sp.flashRadius);
    }
  }

  dispose(): void {
    this.core.geometry.dispose();
    this.core.material.dispose();
    this.glare.material.dispose();
    this.flash.material.dispose();
    this.glareTex.dispose();
    this.electrodes.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.isMesh) {
        m.geometry.dispose();
        (m.material as THREE.Material).dispose();
      }
    });
  }
}
