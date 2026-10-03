/**
 * The 3D stage: renderer, camera + orbit controls, image-based lighting, key/fill
 * lights with soft shadows, a fading floor, and the post-processing chain
 * (RenderPass → EmitterBloomPass (UnrealBloomPass on light sources only) → OutputPass).
 *
 * Colour pipeline: everything renders linear into a half-float, multisampled
 * target; OutputPass applies the tone mapping and the sRGB transfer once, at the end.
 *
 * Bloom is selective: it stands in for the glare of the combustion emitters
 * (flame front, spark, autoignition flash), not for the engine's specular
 * glints, which reach similar HDR values under the key light. Objects registered
 * with `addEmitters` are moved to EMITTER_LAYER; the bloom input is those objects
 * alone, depth-tested against a depth-only pre-pass of everything else, so a
 * flame hidden behind the (uncut) cylinder wall does not glow through it.
 *
 * Grading modes:
 *  - 'filmic' (default): ACES filmic tone mapping + bloom, for the physically
 *    exposed view (flame chemiluminescence, spark, gray-body glow).
 *  - 'exact': no tone mapping and no bloom, so false-colour output (the gas
 *    temperature colour map) reaches the screen unchanged and matches its legend.
 *    Scene lights are dimmed a little so the engine does not clip without the
 *    filmic shoulder.
 *
 * The stage knows nothing about engines or snapshots; the app places content in
 * `scene` and calls `fitToBounds` / `setFraming`.
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/examples/jsm/environments/RoomEnvironment.js';
import { EffectComposer } from 'three/examples/jsm/postprocessing/EffectComposer.js';
import { OutputPass } from 'three/examples/jsm/postprocessing/OutputPass.js';
import { FullScreenQuad, Pass } from 'three/examples/jsm/postprocessing/Pass.js';
import { RenderPass } from 'three/examples/jsm/postprocessing/RenderPass.js';
import { UnrealBloomPass } from 'three/examples/jsm/postprocessing/UnrealBloomPass.js';

/** Layer of light-emitting objects (bloom sources). The main camera sees layers 0–2. */
export const EMITTER_LAYER = 1;
/** Layer of decoration that must neither glow nor occlude emitters (the floor). */
export const DECOR_LAYER = 2;

/** A camera placement (world coordinates). Structurally compatible with EngineModel's CameraView. */
export interface CameraFraming {
  position: THREE.Vector3;
  target: THREE.Vector3;
  /** Vertical field of view, degrees. */
  fov: number;
  near: number;
  far: number;
}

export type Grading = 'filmic' | 'exact';

export interface StageOptions {
  /** Device-pixel-ratio cap (default 2). */
  maxPixelRatio?: number;
  /** MSAA samples of the HDR target (default 4; 0 disables). */
  msaaSamples?: number;
  /** Bloom settings (threshold = linear HDR luminance of the emitter image). */
  bloom?: { strength: number; radius: number; threshold: number };
  /** Bloom input resolution relative to the canvas (default 0.5). */
  bloomResolution?: number;
  /** Soft shadows from the key light (default true). */
  shadows?: boolean;
}

/**
 * Bloom on the emitter image. At the combustion visuals' VISUAL_GAIN the flame
 * front and spark are ≳ 8 display units, the burned-gas glow and the (non-
 * emissive) flow tracers stay below 1, so a threshold of 1 blooms only real light.
 */
const DEFAULT_BLOOM = { strength: 0.7, radius: 0.35, threshold: 1.0 };
/** Light scale applied in 'exact' grading (no filmic shoulder to absorb highlights). */
const EXACT_LIGHT_SCALE = 0.7;
const TRANSITION_SECONDS = 0.7;

const easeInOutCubic = (x: number): number => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2);

interface Transition {
  fromPos: THREE.Vector3;
  fromTarget: THREE.Vector3;
  fromFov: number;
  to: CameraFraming;
  elapsed: number;
  duration: number;
}

export class Stage {
  readonly renderer: THREE.WebGLRenderer;
  readonly scene: THREE.Scene;
  readonly camera: THREE.PerspectiveCamera;
  readonly controls: OrbitControls;
  readonly composer: EffectComposer;
  readonly bloomPass: UnrealBloomPass;
  readonly emitterBloom: EmitterBloomPass;
  readonly outputPass: OutputPass;

  private readonly container: HTMLElement;
  private readonly maxPixelRatio: number;
  private readonly keyLight: THREE.DirectionalLight;
  private readonly fillLight: THREE.DirectionalLight;
  private readonly rimLight: THREE.DirectionalLight;
  private readonly lightBase: { key: number; fill: number; rim: number; env: number };
  private readonly floor: THREE.Mesh<THREE.CircleGeometry, THREE.MeshStandardMaterial>;
  private readonly envTexture: THREE.Texture;
  private readonly backgroundTexture: THREE.CanvasTexture | null;
  private readonly floorAlpha: THREE.CanvasTexture | null;
  private readonly resizeObserver: ResizeObserver | null = null;
  private transition: Transition | null = null;
  private grading: Grading = 'filmic';
  private bloomWanted = true;
  private width = 0;
  private height = 0;
  private pixelRatio = 1;
  private disposed = false;

  constructor(container: HTMLElement, options: StageOptions = {}) {
    this.container = container;
    this.maxPixelRatio = options.maxPixelRatio ?? 2;
    const bloom = options.bloom ?? DEFAULT_BLOOM;
    const shadows = options.shadows ?? true;

    // ---- renderer ----
    const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'high-performance' });
    this.renderer = renderer;
    renderer.outputColorSpace = THREE.SRGBColorSpace;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = 1.0;
    renderer.localClippingEnabled = true;
    renderer.shadowMap.enabled = shadows;
    // PCFSoftShadowMap was folded into PCFShadowMap (soft via shadow.radius) in r18x.
    renderer.shadowMap.type = THREE.PCFShadowMap;
    this.pixelRatio = Math.min(window.devicePixelRatio || 1, this.maxPixelRatio);
    renderer.setPixelRatio(this.pixelRatio);
    renderer.domElement.classList.add('oct-canvas');
    renderer.domElement.setAttribute('aria-label', '3D view of the engine and its combustion chamber');
    container.appendChild(renderer.domElement);

    // ---- scene, environment ----
    const scene = (this.scene = new THREE.Scene());
    this.backgroundTexture = makeBackgroundTexture();
    scene.background = this.backgroundTexture ?? new THREE.Color(0x10141a);
    const pmrem = new THREE.PMREMGenerator(renderer);
    this.envTexture = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
    pmrem.dispose();
    scene.environment = this.envTexture;
    scene.environmentIntensity = 0.85;

    // ---- lights ----
    const key = (this.keyLight = new THREE.DirectionalLight(0xfff1e0, 2.2));
    key.name = 'key-light';
    key.position.set(1.3, 2.1, 1.7);
    key.castShadow = shadows;
    key.shadow.mapSize.set(2048, 2048);
    key.shadow.radius = 4;
    key.shadow.bias = -0.0004;
    key.shadow.normalBias = 0.0015;
    scene.add(key, key.target);

    const fill = (this.fillLight = new THREE.DirectionalLight(0xc8d8ff, 0.55));
    fill.name = 'fill-light';
    fill.position.set(-1.6, 0.9, 0.8);
    scene.add(fill, fill.target);

    const rim = (this.rimLight = new THREE.DirectionalLight(0xdfe8ff, 0.8));
    rim.name = 'rim-light';
    rim.position.set(-0.6, 1.4, -1.8);
    scene.add(rim, rim.target);
    this.lightBase = { key: key.intensity, fill: fill.intensity, rim: rim.intensity, env: scene.environmentIntensity };

    // ---- floor (fades into the background) ----
    this.floorAlpha = makeRadialAlphaTexture();
    const floorMat = new THREE.MeshStandardMaterial({
      color: 0x1b2129,
      roughness: 0.92,
      metalness: 0,
      transparent: this.floorAlpha !== null,
      alphaMap: this.floorAlpha,
      depthWrite: false,
      envMapIntensity: 0.4,
    });
    this.floor = new THREE.Mesh(new THREE.CircleGeometry(1, 96), floorMat);
    this.floor.name = 'floor';
    this.floor.rotation.x = -Math.PI / 2;
    this.floor.receiveShadow = shadows;
    this.floor.renderOrder = -1;
    this.floor.layers.set(DECOR_LAYER);
    scene.add(this.floor);

    // ---- camera + controls ----
    const camera = (this.camera = new THREE.PerspectiveCamera(34, 1, 0.005, 20));
    camera.layers.enable(EMITTER_LAYER);
    camera.layers.enable(DECOR_LAYER);
    camera.position.set(0.8, 0.8, 1.3);
    const controls = (this.controls = new OrbitControls(camera, renderer.domElement));
    controls.enableDamping = true;
    controls.dampingFactor = 0.09;
    controls.screenSpacePanning = true;
    controls.zoomToCursor = true;
    controls.minDistance = 0.02;
    controls.maxDistance = 6;
    controls.rotateSpeed = 0.7;
    controls.addEventListener('start', this.cancelTransition);

    // ---- post-processing ----
    const samples = Math.max(0, options.msaaSamples ?? 4);
    const target = new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType, samples });
    target.texture.name = 'Stage.hdr';
    this.composer = new EffectComposer(renderer, target);
    this.composer.addPass(new RenderPass(scene, camera));
    this.bloomPass = new UnrealBloomPass(new THREE.Vector2(1, 1), bloom.strength, bloom.radius, bloom.threshold);
    this.emitterBloom = new EmitterBloomPass(scene, camera, this.bloomPass, options.bloomResolution ?? 0.5);
    this.composer.addPass(this.emitterBloom);
    this.outputPass = new OutputPass();
    this.composer.addPass(this.outputPass);

    // ---- sizing ----
    window.addEventListener('resize', this.resize);
    if (typeof ResizeObserver !== 'undefined') {
      this.resizeObserver = new ResizeObserver(() => this.resize());
      this.resizeObserver.observe(container);
    }
    this.resize();
  }

  /**
   * Register light-emitting content (e.g. the combustion visuals' root): every
   * object in the subtree moves to EMITTER_LAYER and becomes a bloom source.
   * Call after the subtree is fully built.
   */
  addEmitters(root: THREE.Object3D): void {
    root.traverse((o) => o.layers.set(EMITTER_LAYER));
  }

  /**
   * Fit the shadow camera and the floor to the content's world bounds (call
   * after adding the engine; bounds need not be tight).
   */
  fitToBounds(box: THREE.Box3): void {
    if (box.isEmpty()) return;
    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const radius = 0.5 * size.length();

    const key = this.keyLight;
    const dir = key.position.clone().sub(key.target.position).normalize();
    key.target.position.copy(center);
    key.position.copy(center).addScaledVector(dir, radius * 3);
    const cam = key.shadow.camera;
    const r = radius * 1.05;
    cam.left = -r;
    cam.right = r;
    cam.top = r;
    cam.bottom = -r;
    cam.near = radius * 0.5;
    cam.far = radius * 5.5;
    cam.updateProjectionMatrix();

    this.fillLight.target.position.copy(center);
    this.rimLight.target.position.copy(center);

    this.floor.position.set(center.x, box.min.y - 0.0005, center.z);
    const fr = Math.max(size.x, size.z) * 1.6;
    this.floor.scale.set(fr, fr, 1);
  }

  /** Move the camera to a framing: animated (default) or immediately. */
  setFraming(view: CameraFraming, animate = true): void {
    if (!animate) {
      this.transition = null;
      this.applyFraming(view, 1, this.camera.position, this.controls.target, this.camera.fov);
      return;
    }
    this.transition = {
      fromPos: this.camera.position.clone(),
      fromTarget: this.controls.target.clone(),
      fromFov: this.camera.fov,
      to: {
        position: view.position.clone(),
        target: view.target.clone(),
        fov: view.fov,
        near: view.near,
        far: view.far,
      },
      elapsed: 0,
      duration: TRANSITION_SECONDS,
    };
    this.camera.near = Math.min(this.camera.near, view.near);
    this.camera.far = Math.max(this.camera.far, view.far);
    this.camera.updateProjectionMatrix();
  }

  /** Translate camera and orbit target together (e.g. to follow the cylinder when the CR changes). */
  shiftView(dx: number, dy: number, dz: number): void {
    this.camera.position.x += dx;
    this.camera.position.y += dy;
    this.camera.position.z += dz;
    this.controls.target.x += dx;
    this.controls.target.y += dy;
    this.controls.target.z += dz;
    const tr = this.transition;
    if (tr) {
      tr.to.position.x += dx;
      tr.to.position.y += dy;
      tr.to.position.z += dz;
      tr.to.target.x += dx;
      tr.to.target.y += dy;
      tr.to.target.z += dz;
    }
  }

  /** 'filmic' = ACES + bloom (physical view); 'exact' = no tone mapping, no bloom (false colour). */
  setGrading(g: Grading): void {
    if (g === this.grading) return;
    this.grading = g;
    const exact = g === 'exact';
    this.renderer.toneMapping = exact ? THREE.NoToneMapping : THREE.ACESFilmicToneMapping;
    const k = exact ? EXACT_LIGHT_SCALE : 1;
    this.keyLight.intensity = this.lightBase.key * k;
    this.fillLight.intensity = this.lightBase.fill * k;
    this.rimLight.intensity = this.lightBase.rim * k;
    this.scene.environmentIntensity = this.lightBase.env * k;
    this.emitterBloom.enabled = this.bloomWanted && !exact;
  }

  get currentGrading(): Grading {
    return this.grading;
  }

  setBloomEnabled(on: boolean): void {
    this.bloomWanted = on;
    this.emitterBloom.enabled = on && this.grading === 'filmic';
  }

  /** Advance camera transitions / damping and draw one frame. */
  render(dtWall: number): void {
    if (this.disposed) return;
    const tr = this.transition;
    if (tr) {
      tr.elapsed += Math.max(0, dtWall);
      const u = Math.min(1, tr.elapsed / tr.duration);
      const e = easeInOutCubic(u);
      this.applyFraming(tr.to, e, tr.fromPos, tr.fromTarget, tr.fromFov);
      if (u >= 1) this.transition = null;
    }
    this.controls.update(dtWall);
    this.composer.render(dtWall);
  }

  readonly resize = (): void => {
    if (this.disposed) return;
    const w = Math.max(1, Math.floor(this.container.clientWidth));
    const h = Math.max(1, Math.floor(this.container.clientHeight));
    const pr = Math.min(window.devicePixelRatio || 1, this.maxPixelRatio);
    if (w === this.width && h === this.height && pr === this.pixelRatio) return;
    this.width = w;
    this.height = h;
    if (pr !== this.pixelRatio) {
      this.pixelRatio = pr;
      this.renderer.setPixelRatio(pr);
      this.composer.setPixelRatio(pr);
    }
    this.renderer.setSize(w, h, false);
    this.composer.setSize(w, h);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  };

  dispose(): void {
    if (this.disposed) return;
    this.disposed = true;
    window.removeEventListener('resize', this.resize);
    this.resizeObserver?.disconnect();
    this.controls.removeEventListener('start', this.cancelTransition);
    this.controls.dispose();
    this.composer.dispose();
    this.emitterBloom.dispose();
    this.outputPass.dispose();
    this.floor.geometry.dispose();
    this.floor.material.dispose();
    this.floorAlpha?.dispose();
    this.backgroundTexture?.dispose();
    this.envTexture.dispose();
    this.keyLight.shadow.dispose();
    this.renderer.dispose();
    this.renderer.domElement.remove();
  }

  // ------------------------------------------------------------------

  private readonly cancelTransition = (): void => {
    const tr = this.transition;
    if (!tr) return;
    // Keep the near/far of the destination so a half-finished move still clips sensibly.
    this.camera.near = tr.to.near;
    this.camera.far = tr.to.far;
    this.camera.updateProjectionMatrix();
    this.transition = null;
  };

  private applyFraming(view: CameraFraming, e: number, fromPos: THREE.Vector3, fromTarget: THREE.Vector3, fromFov: number): void {
    const cam = this.camera;
    if (e >= 1) {
      cam.position.copy(view.position);
      this.controls.target.copy(view.target);
      cam.fov = view.fov;
      cam.near = view.near;
      cam.far = view.far;
    } else {
      cam.position.lerpVectors(fromPos, view.position, e);
      this.controls.target.lerpVectors(fromTarget, view.target, e);
      cam.fov = fromFov + (view.fov - fromFov) * e;
    }
    cam.updateProjectionMatrix();
  }
}

// ---------------------------------------------------------------------------
// Emitter-only bloom
// ---------------------------------------------------------------------------

/**
 * Bloom from the objects on EMITTER_LAYER only, added onto the HDR scene image.
 *   1. depth-only pre-pass of layer 0 (the engine; double-sided so section caps,
 *      which the engine draws from back faces, occlude too) into a small target;
 *   2. the emitters with their own materials on top (so they are occluded by the
 *      engine exactly as in the main image, minus MSAA);
 *   3. UnrealBloomPass on that image; its blurred mip composite (before its own
 *      blend) is added onto the scene with additive blending.
 */
export class EmitterBloomPass extends Pass {
  private readonly target: THREE.WebGLRenderTarget;
  private readonly depthOnly = new THREE.MeshBasicMaterial({ colorWrite: false, side: THREE.DoubleSide });
  private readonly addMaterial: THREE.MeshBasicMaterial;
  private readonly quad: FullScreenQuad;
  private readonly clearColor = new THREE.Color();

  constructor(
    private readonly scene: THREE.Scene,
    private readonly camera: THREE.Camera,
    private readonly bloom: UnrealBloomPass,
    private readonly resolutionScale = 0.5,
  ) {
    super();
    this.needsSwap = false;
    this.target = new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType });
    this.target.texture.name = 'EmitterBloom.input';
    // Pure additive (One, One): the bloom composite's alpha is not a coverage value.
    this.addMaterial = new THREE.MeshBasicMaterial({
      blending: THREE.CustomBlending,
      blendEquation: THREE.AddEquation,
      blendSrc: THREE.OneFactor,
      blendDst: THREE.OneFactor,
      depthTest: false,
      depthWrite: false,
      transparent: true,
      toneMapped: false,
    });
    this.quad = new FullScreenQuad(this.addMaterial);
  }

  override setSize(width: number, height: number): void {
    const w = Math.max(1, Math.round(width * this.resolutionScale));
    const h = Math.max(1, Math.round(height * this.resolutionScale));
    this.target.setSize(w, h);
    this.bloom.setSize(w, h);
  }

  override render(
    renderer: THREE.WebGLRenderer,
    _writeBuffer: THREE.WebGLRenderTarget,
    readBuffer: THREE.WebGLRenderTarget,
    deltaTime: number,
  ): void {
    const scene = this.scene;
    const cam = this.camera;
    const mask = cam.layers.mask;
    const background = scene.background;
    const override = scene.overrideMaterial;
    const autoClear = renderer.autoClear;
    const shadowAuto = renderer.shadowMap.autoUpdate;
    const clearAlpha = renderer.getClearAlpha();
    renderer.getClearColor(this.clearColor);

    renderer.autoClear = false;
    renderer.shadowMap.autoUpdate = false; // shadow maps were rendered by the main pass
    scene.background = null;
    renderer.setRenderTarget(this.target);
    renderer.setClearColor(0x000000, 0);
    renderer.clear(true, true, false);

    cam.layers.set(0);
    scene.overrideMaterial = this.depthOnly;
    renderer.render(scene, cam);
    cam.layers.set(EMITTER_LAYER);
    scene.overrideMaterial = null;
    renderer.render(scene, cam);

    cam.layers.mask = mask;
    scene.overrideMaterial = override;
    scene.background = background;
    renderer.shadowMap.autoUpdate = shadowAuto;
    renderer.setClearColor(this.clearColor, clearAlpha);

    // Bloom of the emitter image; the pure blurred composite is left in renderTargetsHorizontal[0].
    this.bloom.render(renderer, this.target, this.target, deltaTime, false);
    this.addMaterial.map = this.bloom.renderTargetsHorizontal[0].texture;
    renderer.setRenderTarget(this.renderToScreen ? null : readBuffer);
    this.quad.render(renderer);
    renderer.autoClear = autoClear;
  }

  override dispose(): void {
    this.target.dispose();
    this.depthOnly.dispose();
    this.addMaterial.dispose();
    this.quad.dispose();
    this.bloom.dispose();
  }
}

// ---------------------------------------------------------------------------
// Procedural textures (null outside a DOM, e.g. in node tests)
// ---------------------------------------------------------------------------

function makeCanvas(w: number, h: number): HTMLCanvasElement | null {
  if (typeof document === 'undefined') return null;
  const c = document.createElement('canvas');
  c.width = w;
  c.height = h;
  return c;
}

/** Soft vignette behind the engine (screen-space background plane). */
function makeBackgroundTexture(): THREE.CanvasTexture | null {
  const c = makeCanvas(256, 256);
  const g = c?.getContext('2d');
  if (!c || !g) return null;
  const grad = g.createRadialGradient(128, 110, 8, 128, 128, 190);
  grad.addColorStop(0, '#2a323d');
  grad.addColorStop(0.55, '#171c23');
  grad.addColorStop(1, '#0c0f13');
  g.fillStyle = grad;
  g.fillRect(0, 0, 256, 256);
  const t = new THREE.CanvasTexture(c);
  t.colorSpace = THREE.SRGBColorSpace;
  return t;
}

/** Radial alpha falloff for the floor disc (alphaMap samples the green channel). */
function makeRadialAlphaTexture(): THREE.CanvasTexture | null {
  const c = makeCanvas(128, 128);
  const g = c?.getContext('2d');
  if (!c || !g) return null;
  const grad = g.createRadialGradient(64, 64, 0, 64, 64, 64);
  grad.addColorStop(0, '#ffffff');
  grad.addColorStop(0.35, '#d0d0d0');
  grad.addColorStop(1, '#000000');
  g.fillStyle = grad;
  g.fillRect(0, 0, 128, 128);
  return new THREE.CanvasTexture(c);
}
