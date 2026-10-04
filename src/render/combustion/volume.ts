/**
 * three.js wrapper for the ray-marched chamber gas (see volume-shader.ts).
 * The proxy geometry is a cylinder slightly inside the chamber, drawn with
 * back faces, so (a) the camera may be inside it and (b) its fragments sit
 * just in front of the far walls and are depth-tested against opaque parts
 * (head, piston, liner when not cut away) that occlude the whole chamber.
 *
 * L-head chambers (non-convex: bore column ∪ valve pocket) use three proxies
 * sharing one set of uniforms: the bore cylinder (back faces), the pocket's
 * rounded-rectangle prism (back faces) and the transfer arc — the bore circle
 * where it opens into the pocket, just outside the piston — drawn with FRONT
 * faces (rays leaving the pocket into the piston's top land). Each is convex
 * (or crossed inward at most once), so each yields ≤ 1 fragment per pixel, and
 * the shader keeps only the fragment of the proxy owning the exit of the last
 * visible gas interval: one fragment per pixel, depth-tested at the gas exit.
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { boundingSize, chamberShapeOf, EXIT_BORE, EXIT_POCKET, EXIT_TRANSFER, type ChamberShape } from './chamber';
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

/** Angular overlap of the transfer-arc proxy beyond the ends of the arc, rad (covers tessellation). */
const TRANSFER_ARC_MARGIN = (2 * Math.PI) / 180;
/** Segments of the transfer arc per radian and of each pocket corner. */
const TRANSFER_ARC_SEGMENTS_PER_RAD = 24;
const POCKET_CORNER_SEGMENTS = 16;

/**
 * Cut-away region of the engine in the cylinder frame: the set {p : n·p > d} for every plane (≤ 2),
 * each plane given as [nx, ny, nz, d]. Metal there is not drawn, so gas behind it stays visible.
 */
export type CutPlanes = readonly (readonly [number, number, number, number])[];

/**
 * Closed rounded-rectangle prism (plan inset by `inset`), y ∈ [y0, y1], cylinder-frame coordinates,
 * every triangle wound counter-clockwise seen from outside.
 */
export function pocketPrismGeometry(s: ChamberShape, inset: number, y0: number, y1: number): THREE.BufferGeometry {
  const r = Math.max(s.rc - inset, 0);
  const shrink = Math.max(inset - s.rc, 0);
  const cx = [s.ix1 - shrink, s.ix0 + shrink, s.ix0 + shrink, s.ix1 - shrink];
  const cz = [s.iz1 - shrink, s.iz1 - shrink, s.iz0 + shrink, s.iz0 + shrink];
  const seg = r > 0 ? POCKET_CORNER_SEGMENTS : 1;
  const ox: number[] = [], oz: number[] = [];
  for (let k = 0; k < 4; k++) {
    for (let i = 0; i <= seg; i++) {
      if (r === 0 && i > 0) break;
      const a = ((k + i / seg) * Math.PI) / 2;
      ox.push(cx[k] + r * Math.cos(a));
      oz.push(cz[k] + r * Math.sin(a));
    }
  }
  const mx = 0.5 * (s.ix0 + s.ix1), mz = 0.5 * (s.iz0 + s.iz1);
  const pos: number[] = [];
  const tri = (a: number[], b: number[], c: number[], out: number[]): void => {
    const ux = b[0] - a[0], uy = b[1] - a[1], uz = b[2] - a[2];
    const vx = c[0] - a[0], vy = c[1] - a[1], vz = c[2] - a[2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    if (nx * out[0] + ny * out[1] + nz * out[2] >= 0) pos.push(...a, ...b, ...c);
    else pos.push(...a, ...c, ...b);
  };
  const n = ox.length;
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    if (Math.hypot(ox[j] - ox[i], oz[j] - oz[i]) < 1e-12) continue;
    tri([mx, y1, mz], [ox[i], y1, oz[i]], [ox[j], y1, oz[j]], [0, 1, 0]);
    tri([mx, y0, mz], [ox[i], y0, oz[i]], [ox[j], y0, oz[j]], [0, -1, 0]);
    const out = [0.5 * (ox[i] + ox[j]) - mx, 0, 0.5 * (oz[i] + oz[j]) - mz];
    tri([ox[i], y0, oz[i]], [ox[j], y0, oz[j]], [ox[j], y1, oz[j]], out);
    tri([ox[i], y0, oz[i]], [ox[j], y1, oz[j]], [ox[i], y1, oz[i]], out);
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(pos, 3));
  geo.computeBoundingSphere();
  return geo;
}

export class GasVolume {
  /** Bore-column proxy (the flat disc's only proxy). */
  readonly mesh: THREE.Mesh<THREE.CylinderGeometry, THREE.ShaderMaterial>;
  /** Every proxy mesh: the bore column, then (L-head) the pocket prism and the transfer arc. */
  readonly meshes: readonly THREE.Mesh<THREE.BufferGeometry, THREE.ShaderMaterial>[];
  readonly shape: ChamberShape;
  private readonly lut: THREE.DataTexture;
  /** Uniforms shared by every proxy material. */
  private readonly uniforms: Record<string, THREE.IUniform>;
  private readonly boreToCyl: THREE.Matrix4;
  private readonly invRoot = new THREE.Matrix4();
  private readonly tmpV = new THREE.Vector3();
  private frame = 0;
  private mode: CombustionMode = 'physical';

  constructor(private readonly spec: EngineSpec, private readonly frameObject: THREE.Object3D) {
    const shape = (this.shape = chamberShapeOf(spec));
    const lhead = shape.kind === 'l-head';
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
      uCamPos: { value: new THREE.Vector3(0, 0.2, 0.2) },
      uCamDir: { value: new THREE.Vector3(0, 0, -1) },
      uOrtho: { value: 0 },
      uR: { value: R },
      uH: { value: 0.01 },
      uValve0: { value: new THREE.Vector4(iv.position[0], iv.position[1], iv.headDiameter / 2, 0) },
      uValve1: { value: new THREE.Vector4(ev.position[0], ev.position[1], ev.headDiameter / 2, 0) },
      uValveSeat0: { value: new THREE.Vector2(iv.seatY ?? 0, iv.liftDirection ?? -1) },
      uValveSeat1: { value: new THREE.Vector2(ev.seatY ?? 0, ev.liftDirection ?? -1) },
      uValveThick0: { value: VALVE_HEAD_THICKNESS_FRACTION * iv.headDiameter },
      uValveThick1: { value: VALVE_HEAD_THICKNESS_FRACTION * ev.headDiameter },
      uMode: { value: 0 },
      uShape: { value: lhead ? 1 : 0 },
      uPocketRect: { value: new THREE.Vector4(shape.ix0, shape.ix1, shape.iz0, shape.iz1) },
      uPocketRad: { value: shape.rc },
      uPocketY: { value: new THREE.Vector2(shape.deckY, shape.roofY) },
      uBackoff: { value: 1 },
      uCutOn: { value: 0 },
      uCut0: { value: new THREE.Vector4(0, 0, 0, -1) },
      uCut1: { value: new THREE.Vector4(0, 0, 0, -1) },
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
      uKnockShape: { value: lhead ? 1 : 0 },
      uKnockLine: { value: new THREE.Vector3(shape.axis[0], shape.axis[1], shape.s0) },
      uKnockLen: { value: shape.s1 - shape.s0 },
      uKnockSign: { value: 1 },
      uKnockOrigin: { value: new THREE.Vector2() },
      uKnockRing: { value: 0 },
      uKnockRingW: { value: 0.005 },
      uKnockStanding: { value: 1 },
      uKnockK: { value: OVERLAY_KNOCK_KAPPA },
      uKnockPos: { value: new THREE.Vector3(...OVERLAY_KNOCK_POSITIVE) },
      uKnockNeg: { value: new THREE.Vector3(...OVERLAY_KNOCK_NEGATIVE) },
      uFrame: { value: 0 },
    };
    this.boreToCyl = new THREE.Matrix4();
    const makeMaterial = (piece: number, meshToCyl: THREE.Matrix4, side: THREE.Side): THREE.ShaderMaterial =>
      new THREE.ShaderMaterial({
        name: 'CombustionGasVolume',
        vertexShader: VOLUME_VERTEX,
        fragmentShader: VOLUME_FRAGMENT,
        uniforms: { ...this.uniforms, uMeshToCyl: { value: meshToCyl }, uPiece: { value: piece } },
        side,
        transparent: true,
        depthWrite: false,
        depthTest: true,
        blending: THREE.CustomBlending,
        blendEquation: THREE.AddEquation,
        blendSrc: THREE.OneFactor,
        blendDst: THREE.OneMinusSrcAlphaFactor,
        premultipliedAlpha: true,
      });
    const onBefore = (_r: THREE.WebGLRenderer, _s: THREE.Scene, camera: THREE.Camera): void => this.beforeRender(camera);

    this.mesh = new THREE.Mesh(geo, makeMaterial(EXIT_BORE, this.boreToCyl, THREE.BackSide));
    this.mesh.name = 'combustion-gas-volume';
    this.mesh.renderOrder = 10;
    this.mesh.scale.set(PROXY_RADIUS_FRACTION * R, 0.01, PROXY_RADIUS_FRACTION * R);
    this.mesh.onBeforeRender = onBefore;
    const meshes: THREE.Mesh<THREE.BufferGeometry, THREE.ShaderMaterial>[] = [this.mesh];

    if (lhead) {
      const y0 = shape.deckY + PROXY_GAP, y1 = shape.roofY - PROXY_GAP;
      const pocket = new THREE.Mesh(
        pocketPrismGeometry(shape, PROXY_GAP, y0, y1),
        makeMaterial(EXIT_POCKET, new THREE.Matrix4(), THREE.BackSide),
      );
      pocket.name = 'combustion-gas-volume-pocket';
      pocket.renderOrder = 10;
      pocket.onBeforeRender = onBefore;
      meshes.push(pocket);
      const [a0, a1] = shape.transferArc;
      if (a1 >= a0 && y1 > y0) {
        // CylinderGeometry vertex (r sin θ, y, r cos θ): plan angle atan2(z, x) = π/2 − θ
        const len = Math.min(a1 - a0 + 2 * TRANSFER_ARC_MARGIN, 2 * Math.PI);
        const segs = Math.max(8, Math.ceil(len * TRANSFER_ARC_SEGMENTS_PER_RAD));
        const strip = new THREE.CylinderGeometry(
          R + PROXY_GAP, R + PROXY_GAP, y1 - y0, segs, 1, true, Math.PI / 2 - a1 - TRANSFER_ARC_MARGIN, len,
        );
        const toCyl = new THREE.Matrix4().makeTranslation(0, 0.5 * (y0 + y1), 0);
        const transfer = new THREE.Mesh(strip, makeMaterial(EXIT_TRANSFER, toCyl, THREE.FrontSide));
        transfer.name = 'combustion-gas-volume-transfer';
        transfer.renderOrder = 10;
        transfer.position.set(0, 0.5 * (y0 + y1), 0);
        transfer.onBeforeRender = onBefore;
        meshes.push(transfer);
      }
    }
    this.meshes = meshes;
  }

  setMode(mode: CombustionMode): void {
    if (mode === this.mode) return;
    this.mode = mode;
    this.uniforms.uMode.value = mode === 'temperature' ? 1 : 0;
    // False colour must match the legend → no tone mapping in 'temperature' mode.
    for (const m of this.meshes) {
      m.material.toneMapped = mode !== 'temperature';
      m.material.needsUpdate = true;
    }
  }

  setTemperatureRange(min: number, max: number): void {
    (this.uniforms.uTRange.value as THREE.Vector2).set(min, Math.max(max, min + 1));
  }

  /** Cut-away region of the engine (cylinder frame), or null for none (metal gaps always stop the ray). */
  setCutRegion(planes: CutPlanes | null): void {
    const u = this.uniforms;
    const p0 = planes?.[0], p1 = planes?.[1];
    u.uCutOn.value = p0 ? 1 : 0;
    (u.uCut0.value as THREE.Vector4).set(p0?.[0] ?? 0, p0?.[1] ?? 0, p0?.[2] ?? 0, p0?.[3] ?? -1);
    (u.uCut1.value as THREE.Vector4).set(p1?.[0] ?? 0, p1?.[1] ?? 0, p1?.[2] ?? 0, p1?.[3] ?? -1);
  }

  sync(st: CombustionVisualState, s: EngineSnapshot): void {
    const u = this.uniforms;
    const h = Math.max(st.h, 2 * PROXY_GAP + 1e-5);
    this.mesh.scale.y = h - 2 * PROXY_GAP;
    this.mesh.position.set(0, -h / 2, 0);
    this.mesh.updateMatrix();
    this.boreToCyl.copy(this.mesh.matrix);
    u.uH.value = h;
    if (this.shape.kind === 'l-head') u.uBackoff.value = 4 * boundingSize(this.shape, h);
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
    u.uKnockSign.value = k.modeSign;
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
    for (const m of this.meshes) {
      m.geometry.dispose();
      m.material.dispose();
    }
    this.lut.dispose();
  }
}
