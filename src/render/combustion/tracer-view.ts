/**
 * three.js Points for the flow tracers. Reads TracerSystem buffers in place
 * (positions are the simulation array itself) and computes per-particle
 * colours: fresh charge = cool cyan, hot residual/exhaust = warm (black-body
 * chromaticity of its temperature), flash as the flame front passes, then
 * burned-gas warm. Tracers are a flow visualisation, not light: normal
 * (premultiplied) alpha blending, no emission.
 */
import * as THREE from 'three';
import { blackbodyXYZFast, clipToGamut, xyzToLinearSrgb, type Vec3 } from './colour/cie';
import { TRACER_SIZE } from './constants';
import { FRONT_FLASH_S, TracerState, type TracerSystem } from './tracers';

const VERT = /* glsl */ `
attribute vec4 aColor;
uniform float uSize;
uniform float uScale;
varying vec4 vColor;
void main() {
  vColor = aColor;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mv;
  float s = isOrthographic ? uSize * uScale : uSize * uScale / max(-mv.z, 1e-4);
  gl_PointSize = aColor.a > 0.0 ? clamp(s, 1.5, 48.0) : 0.0;
}
`;

const FRAG = /* glsl */ `
varying vec4 vColor;
void main() {
  vec2 d = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(d, d);
  if (r2 > 1.0) discard;
  float a = vColor.a * (1.0 - smoothstep(0.35, 1.0, r2));
  gl_FragColor = vec4(vColor.rgb * a, a);
  #include <colorspace_fragment>
}
`;

const FRESH: Vec3 = [0.3, 0.72, 1.0];
const FRONT_FLASH: Vec3 = [0.75, 0.85, 1.0];
const COOL_EXHAUST: Vec3 = [0.55, 0.55, 0.58];

/** Warm display colour (max channel 1) of black-body light at T; grey below ~800 K. */
export function warmColor(T: number, out: Vec3, tmp: Vec3 = [0, 0, 0]): Vec3 {
  const Tc = Math.min(Math.max(T, 900), 4000);
  blackbodyXYZFast(Tc, tmp);
  clipToGamut(xyzToLinearSrgb(tmp, out));
  const m = Math.max(out[0], out[1], out[2]) || 1;
  const w = Math.min(Math.max((T - 700) / 700, 0), 1);
  for (let i = 0; i < 3; i++) out[i] = COOL_EXHAUST[i] + (out[i] / m - COOL_EXHAUST[i]) * w;
  return out;
}

export class TracerView {
  readonly points: THREE.Points<THREE.BufferGeometry, THREE.ShaderMaterial>;
  private readonly colors: Float32Array;
  private readonly posAttr: THREE.BufferAttribute;
  private readonly colAttr: THREE.BufferAttribute;
  private readonly size = new THREE.Vector2();
  private readonly burnedCol: Vec3 = [1, 0.6, 0.3];
  private readonly tmp: Vec3 = [0, 0, 0];
  private readonly tmp2: Vec3 = [0, 0, 0];

  constructor(private readonly sys: TracerSystem) {
    const geo = new THREE.BufferGeometry();
    this.posAttr = new THREE.BufferAttribute(sys.position, 3);
    this.posAttr.setUsage(THREE.DynamicDrawUsage);
    this.colors = new Float32Array(4 * sys.capacity);
    this.colAttr = new THREE.BufferAttribute(this.colors, 4);
    this.colAttr.setUsage(THREE.DynamicDrawUsage);
    geo.setAttribute('position', this.posAttr);
    geo.setAttribute('aColor', this.colAttr);
    geo.setDrawRange(0, 0);
    const mat = new THREE.ShaderMaterial({
      name: 'FlowTracers',
      vertexShader: VERT,
      fragmentShader: FRAG,
      uniforms: { uSize: { value: TRACER_SIZE }, uScale: { value: 500 } },
      transparent: true,
      depthWrite: false,
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneMinusSrcAlphaFactor,
      premultipliedAlpha: true,
    });
    this.points = new THREE.Points(geo, mat);
    this.points.name = 'combustion-flow-tracers';
    this.points.frustumCulled = false;
    this.points.renderOrder = 11;
    this.points.onBeforeRender = (renderer, _s, camera) => {
      renderer.getDrawingBufferSize(this.size);
      const proj = (camera as THREE.PerspectiveCamera).projectionMatrix;
      const worldScale = this.points.matrixWorld.getMaxScaleOnAxis();
      mat.uniforms.uScale.value = proj.elements[5] * this.size.y * 0.5 * worldScale;
    };
  }

  /** Recompute colours for the alive particles and flag buffers for upload. */
  sync(burnedT: number): void {
    const s = this.sys;
    const n = s.count;
    warmColor(burnedT, this.burnedCol, this.tmp);
    const c = this.colors;
    for (let i = 0; i < n; i++) {
      const k = 4 * i;
      const life = s.life[i] > 0 ? s.life[i] : 1;
      const x = s.age[i] / life;
      let a = 0.85 * Math.min(1, x / 0.03) * Math.min(1, (1 - x) / 0.2);
      let r: number, g: number, b: number;
      const bt = s.burnTime[i];
      if (!Number.isNaN(bt)) {
        const since = s.t - bt;
        const w = Math.min(Math.max(since / FRONT_FLASH_S, 0), 1);
        r = FRONT_FLASH[0] + (this.burnedCol[0] - FRONT_FLASH[0]) * w;
        g = FRONT_FLASH[1] + (this.burnedCol[1] - FRONT_FLASH[1]) * w;
        b = FRONT_FLASH[2] + (this.burnedCol[2] - FRONT_FLASH[2]) * w;
        a = Math.min(1, a * (1.6 - 0.6 * w));
      } else if (s.hot[i]) {
        warmColor(s.spawnT[i], this.tmp2, this.tmp);
        r = this.tmp2[0]; g = this.tmp2[1]; b = this.tmp2[2];
      } else {
        r = FRESH[0]; g = FRESH[1]; b = FRESH[2];
      }
      if (s.state[i] === TracerState.PortOut || s.state[i] === TracerState.PortIn) {
        // fade toward the far end of the port section (distance from the seat along the port axis)
        const v = s.valve[i];
        const L = s.portLen[v];
        const y = (s.position[3 * i + 1] - s.seatY[v]) * -s.openDir[v];
        a *= Math.min(1, Math.max(0, (1 - y / L) / 0.5));
      }
      c[k] = r; c[k + 1] = g; c[k + 2] = b; c[k + 3] = Math.max(a, 0);
    }
    this.points.geometry.setDrawRange(0, n);
    if (n > 0) {
      this.posAttr.clearUpdateRanges();
      this.posAttr.addUpdateRange(0, 3 * n);
      this.posAttr.needsUpdate = true;
      this.colAttr.clearUpdateRanges();
      this.colAttr.addUpdateRange(0, 4 * n);
      this.colAttr.needsUpdate = true;
    }
  }

  setVisible(on: boolean): void {
    this.points.visible = on;
  }

  dispose(): void {
    this.points.geometry.dispose();
    this.points.material.dispose();
  }
}
