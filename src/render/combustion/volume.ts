/**
 * three.js wrapper for the ray-marched chamber gas (see volume-shader.ts).
 * The proxy geometry is a cylinder slightly inside the chamber, drawn with
 * back faces, so (a) the camera may be inside it and (b) its fragments sit
 * just in front of the far walls and are depth-tested against opaque parts
 * (head, piston, liner when not cut away) that occlude the whole chamber.
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { infernoLUT } from './colour/colormap';
import {
  OVERLAY_KNOCK_KAPPA,
  OVERLAY_KNOCK_NEGATIVE,
  OVERLAY_KNOCK_POSITIVE,
  VIS_TEMPERATURE_KAPPA,
  VIS_TEMPERATURE_RANGE,
  VIS_UNBURNED_HAZE_COLOR,
  VIS_UNBURNED_HAZE_KAPPA,
} from './constants';
import type { CombustionMode, CombustionVisualState } from './state';
import { VOLUME_FRAGMENT, VOLUME_VERTEX } from './volume-shader';

/** Proxy shrink so its faces never coincide with the engine's liner/head/crown. */
const PROXY_RADIUS_FRACTION = 0.985;
const PROXY_GAP = 2e-4; // m

/** Valve-head thickness used for occlusion, as a fraction of the head diameter. Approximate. */
const VALVE_HEAD_THICKNESS_FRACTION = 0.12;

export class GasVolume {
  readonly mesh: THREE.Mesh<THREE.CylinderGeometry, THREE.ShaderMaterial>;
  private readonly lut: THREE.DataTexture;
  private readonly uniforms: Record<string, THREE.IUniform>;
  private readonly invRoot = new THREE.Matrix4();
  private readonly tmpV = new THREE.Vector3();
  private frame = 0;
  private mode: CombustionMode = 'physical';

  constructor(private readonly spec: EngineSpec, private readonly frameObject: THREE.Object3D) {
    const R = spec.geometry.bore / 2;
    const geo = new THREE.CylinderGeometry(1, 1, 1, 96, 1, false);
    this.lut = new THREE.DataTexture(infernoLUT(256), 256, 1, THREE.RGBAFormat, THREE.UnsignedByteType);
    this.lut.colorSpace = THREE.SRGBColorSpace;
    this.lut.magFilter = THREE.LinearFilter;
    this.lut.minFilter = THREE.LinearFilter;
    this.lut.wrapS = THREE.ClampToEdgeWrapping;
    this.lut.needsUpdate = true;

    const iv = spec.intakeValve, ev = spec.exhaustValve;
    this.uniforms = {
      uMeshToCyl: { value: new THREE.Matrix4() },
      uCamPos: { value: new THREE.Vector3(0, 0.2, 0.2) },
      uCamDir: { value: new THREE.Vector3(0, 0, -1) },
      uOrtho: { value: 0 },
      uR: { value: R },
      uH: { value: 0.01 },
      uValve0: { value: new THREE.Vector4(iv.position[0], iv.position[1], iv.headDiameter / 2, 0) },
      uValve1: { value: new THREE.Vector4(ev.position[0], ev.position[1], ev.headDiameter / 2, 0) },
      uValveThick0: { value: VALVE_HEAD_THICKNESS_FRACTION * iv.headDiameter },
      uValveThick1: { value: VALVE_HEAD_THICKNESS_FRACTION * ev.headDiameter },
      uMode: { value: 0 },
      uFlameOn: { value: 0 },
      uAllBurned: { value: 0 },
      uFlameC: { value: new THREE.Vector3() },
      uFlameR: { value: 0 },
      uBrush: { value: 1e-3 },
      uWrinkleAmp: { value: 0 },
      uWrinkleK: { value: 250 },
      uNoiseT: { value: 0 },
      uFrontJ: { value: new THREE.Vector3() },
      uBurnedJ: { value: new THREE.Vector3() },
      uEndGasJ: { value: new THREE.Vector3() },
      uHazeK: { value: VIS_UNBURNED_HAZE_KAPPA },
      uHazeC: { value: new THREE.Vector3(...VIS_UNBURNED_HAZE_COLOR) },
      uHazeU: { value: 1 },
      uHazeCore: { value: 0 },
      uTu: { value: 300 },
      uTb: { value: 300 },
      uTRange: { value: new THREE.Vector2(VIS_TEMPERATURE_RANGE[0], VIS_TEMPERATURE_RANGE[1]) },
      uTempK: { value: VIS_TEMPERATURE_KAPPA },
      uColormap: { value: this.lut },
      uKnockAmp: { value: 0 },
      uKnockTemporal: { value: 0 },
      uKnockRms: { value: 0 },
      uKnockAxis: { value: 0 },
      uKnockOrigin: { value: new THREE.Vector2() },
      uKnockRing: { value: 0 },
      uKnockRingW: { value: 0.005 },
      uKnockStanding: { value: 1 },
      uKnockK: { value: OVERLAY_KNOCK_KAPPA },
      uKnockPos: { value: new THREE.Vector3(...OVERLAY_KNOCK_POSITIVE) },
      uKnockNeg: { value: new THREE.Vector3(...OVERLAY_KNOCK_NEGATIVE) },
      uFrame: { value: 0 },
    };
    const mat = new THREE.ShaderMaterial({
      name: 'CombustionGasVolume',
      vertexShader: VOLUME_VERTEX,
      fragmentShader: VOLUME_FRAGMENT,
      uniforms: this.uniforms,
      side: THREE.BackSide,
      transparent: true,
      depthWrite: false,
      depthTest: true,
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneMinusSrcAlphaFactor,
      premultipliedAlpha: true,
    });
    this.mesh = new THREE.Mesh(geo, mat);
    this.mesh.name = 'combustion-gas-volume';
    this.mesh.renderOrder = 10;
    this.mesh.scale.set(PROXY_RADIUS_FRACTION * R, 0.01, PROXY_RADIUS_FRACTION * R);
    this.mesh.onBeforeRender = (_r, _s, camera) => this.beforeRender(camera);
  }

  setMode(mode: CombustionMode): void {
    if (mode === this.mode) return;
    this.mode = mode;
    this.uniforms.uMode.value = mode === 'temperature' ? 1 : 0;
    // False colour must match the legend → no tone mapping in 'temperature' mode.
    this.mesh.material.toneMapped = mode !== 'temperature';
    this.mesh.material.needsUpdate = true;
  }

  setTemperatureRange(min: number, max: number): void {
    (this.uniforms.uTRange.value as THREE.Vector2).set(min, Math.max(max, min + 1));
  }

  sync(st: CombustionVisualState, s: EngineSnapshot): void {
    const u = this.uniforms;
    const h = Math.max(st.h, 2 * PROXY_GAP + 1e-5);
    this.mesh.scale.y = h - 2 * PROXY_GAP;
    this.mesh.position.set(0, -h / 2, 0);
    this.mesh.updateMatrix();
    (u.uMeshToCyl.value as THREE.Matrix4).copy(this.mesh.matrix);
    u.uH.value = h;
    (u.uValve0.value as THREE.Vector4).w = Math.max(s.intakeLift, 0);
    (u.uValve1.value as THREE.Vector4).w = Math.max(s.exhaustLift, 0);

    u.uFlameOn.value = st.flameVisible ? 1 : 0;
    u.uAllBurned.value = st.allBurned ? 1 : 0;
    (u.uFlameC.value as THREE.Vector3).fromArray(st.flameCenter);
    u.uFlameR.value = st.flameRadius;
    u.uBrush.value = Math.max(st.brush.thickness, 1e-6);
    u.uWrinkleAmp.value = st.brush.wrinkleAmplitude;
    u.uWrinkleK.value = st.brush.wrinkleWavelength > 0 ? 1 / st.brush.wrinkleWavelength : 0;
    u.uNoiseT.value = st.noisePhase;
    (u.uFrontJ.value as THREE.Vector3).fromArray(st.frontJ);
    (u.uBurnedJ.value as THREE.Vector3).fromArray(st.burnedJ);
    (u.uEndGasJ.value as THREE.Vector3).fromArray(st.endGasJ);
    u.uHazeU.value = st.hazeUnburned;
    u.uHazeCore.value = st.hazeCore;
    u.uTu.value = st.Tu;
    u.uTb.value = st.Tb;

    const k = st.knock;
    u.uKnockAmp.value = k.amplitude > 1e-3 ? k.amplitude : 0;
    u.uKnockTemporal.value = k.temporal;
    u.uKnockRms.value = k.rms;
    u.uKnockAxis.value = k.axisAngle;
    (u.uKnockOrigin.value as THREE.Vector2).set(k.origin[0], k.origin[1]);
    u.uKnockRing.value = k.ringRadius;
    u.uKnockRingW.value = k.ringWidth;
    u.uKnockStanding.value = k.standing;
    this.frame = (this.frame + 1) % 64;
    u.uFrame.value = this.frame;
  }

  private beforeRender(camera: THREE.Camera): void {
    this.invRoot.copy(this.frameObject.matrixWorld).invert();
    const cp = this.uniforms.uCamPos.value as THREE.Vector3;
    cp.setFromMatrixPosition(camera.matrixWorld).applyMatrix4(this.invRoot);
    const cd = this.uniforms.uCamDir.value as THREE.Vector3;
    camera.getWorldDirection(this.tmpV);
    cd.copy(this.tmpV).transformDirection(this.invRoot);
    this.uniforms.uOrtho.value = (camera as THREE.OrthographicCamera).isOrthographicCamera ? 1 : 0;
  }

  dispose(): void {
    this.mesh.geometry.dispose();
    this.mesh.material.dispose();
    this.lut.dispose();
  }
}
