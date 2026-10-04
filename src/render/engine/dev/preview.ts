/**
 * Stand-alone visual check of the engine model (open
 * /src/render/engine/dev/preview.html on the Vite dev server). Drives the
 * model with an exact slider-crank + model valve-lift "snapshot", no physics.
 * Keys: space pause · c cutaway · [ ] compression ratio · 1/2/3 views · ,/. step 1°
 * URL params: ?theta=deg&view=engine|chamber|back|valves&cut=0|1&cr=8&paused=1&rpm=60
 * Model T: ?engine=modelt (keys a/z spark lever, t/g throttle, 1 engine / 2 chamber / 3 back views;
 * &cyl=k frames cylinder k in the chamber view; &percyl=0 drives it without s.cylinders).
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/examples/jsm/environments/RoomEnvironment.js';
import { CFR_F1 } from '../../../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../../../physics/engines/model-t';
import type { CylinderSnapshot, EngineSnapshot } from '../../../physics/core/snapshot';
import { cylinderAngleDeg } from '../../../physics/core/engine-spec';
import { EngineModel, recommendedCameraView } from '../index';
import { sliderCrank } from '../kinematics';
import { modelValveLift, wrapDeg720 } from '../cam-profile';
import { ModelTEngineModel } from '../../engine-modelt/model';
import { pistonDisplacementAt } from '../../engine-modelt/layout';

const q = new URLSearchParams(location.search);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.toneMappingExposure = 1.0;
renderer.shadowMap.enabled = true;
renderer.shadowMap.type = THREE.PCFShadowMap;
document.body.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x1b1e22);
const pmrem = new THREE.PMREMGenerator(renderer);
scene.environment = pmrem.fromScene(new RoomEnvironment(), 0.04).texture;
scene.environmentIntensity = 0.9;
const key = new THREE.DirectionalLight(0xfff4e6, 2.2);
key.position.set(1.2, 2.0, 1.6);
key.castShadow = true;
key.shadow.mapSize.set(2048, 2048);
key.shadow.camera.left = -1; key.shadow.camera.right = 1; key.shadow.camera.top = 1; key.shadow.camera.bottom = -1;
key.shadow.bias = -0.0004;
scene.add(key);

// ?size=WxH renders at a fixed resolution scaled to the window width (screenshots from narrow panes)
const fixed = (q.get('size') ?? '').split('x').map(Number);
const fixedSize = fixed.length === 2 && fixed[0] > 0 && fixed[1] > 0;
const camera = new THREE.PerspectiveCamera(34, innerWidth / innerHeight, 0.005, 30);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
function fit(): void {
  if (fixedSize) {
    renderer.setPixelRatio(1);
    renderer.setSize(fixed[0], fixed[1], false);
    renderer.domElement.style.width = '100%';
    renderer.domElement.style.height = 'auto';
    camera.aspect = fixed[0] / fixed[1];
  } else {
    camera.aspect = innerWidth / innerHeight;
    renderer.setSize(innerWidth, innerHeight);
  }
  camera.updateProjectionMatrix();
}
fit();
window.addEventListener('resize', fit);
const hud = document.getElementById('hud')!;
let theta = Number(q.get('theta') ?? -90);
let paused = q.get('paused') === '1';
const rpm = Number(q.get('rpm') ?? 30);

if (q.get('engine') === 'modelt') runModelT();
else runCfr();

function runModelT(): void {
  let spark = Number(q.get('spark') ?? MODEL_T_CRUISE.sparkAdvanceDeg);
  let throttle = Number(q.get('throttle') ?? MODEL_T_CRUISE.throttle);
  const model = new ModelTEngineModel(MODEL_T, { ...MODEL_T_CRUISE, sparkAdvanceDeg: spark, throttle });
  scene.add(model.root);
  const L = model.layout;
  const perCyl = q.get('percyl') !== '0';
  const cylIdx = Number(q.get('cyl') ?? 0);
  // spark-gap markers in every cylinder frame
  for (const f of model.cylinderFrames) {
    const mk = new THREE.Mesh(new THREE.SphereGeometry(MODEL_T.sparkPlug.gap * 0.6, 12, 8), new THREE.MeshBasicMaterial({ color: 0x66ccff }));
    mk.position.set(...MODEL_T.sparkPlug.gapCenter);
    f.add(mk);
  }
  const setView = (which: string) => {
    if (which === 'back') {
      camera.position.set(0.9, 0.55, -1.3);
      controls.target.set(0, 0.1, -0.2);
    } else {
      const v = model.cameraView(which === 'chamber' ? 'chamber' : 'engine', cylIdx);
      camera.position.copy(v.position);
      controls.target.copy(v.target);
      camera.fov = v.fov;
      camera.updateProjectionMatrix();
    }
    controls.update();
  };
  setView(q.get('view') ?? 'engine');
  if (q.get('cam') && q.get('tgt')) {
    const c = q.get('cam')!.split(',').map(Number), t = q.get('tgt')!.split(',').map(Number);
    camera.position.set(c[0], c[1], c[2]);
    controls.target.set(t[0], t[1], t[2]);
    if (q.get('fov')) { camera.fov = Number(q.get('fov')); camera.updateProjectionMatrix(); }
    controls.update();
  }
  model.setCutaway(q.get('cut') !== '0');
  const cyls: CylinderSnapshot[] = L.axisZ.map((_, i) => ({ index: i, thetaDeg: 0, cycle: 0, pistonDisplacement: 0, intakeLift: 0, exhaustLift: 0 }) as unknown as CylinderSnapshot);
  const s = { cylinders: perCyl ? cyls : undefined } as EngineSnapshot;
  const iv = L.valves.find((v) => v.kind === 'intake')!.cam, ev = L.valves.find((v) => v.kind === 'exhaust')!.cam;
  const makeSnap = (th: number): EngineSnapshot => {
    s.thetaDeg = wrapDeg720(th);
    for (let i = 0; i < L.nCyl; i++) {
      const ti = cylinderAngleDeg(MODEL_T, i, s.thetaDeg);
      const c = cyls[i];
      c.thetaDeg = ti;
      c.pistonDisplacement = pistonDisplacementAt(L, ti);
      c.intakeLift = iv.valveLift(ti);
      c.exhaustLift = ev.valveLift(ti);
    }
    s.pistonDisplacement = cyls[0].pistonDisplacement;
    s.intakeLift = cyls[0].intakeLift;
    s.exhaustLift = cyls[0].exhaustLift;
    return s;
  };
  window.addEventListener('keydown', (e: KeyboardEvent) => {
    if (e.key === ' ') paused = !paused;
    else if (e.key === 'c') model.setCutaway(!model.isCutaway);
    else if (e.key === '1') setView('engine');
    else if (e.key === '2') setView('chamber');
    else if (e.key === '3') setView('back');
    else if (e.key === ',') theta -= 1;
    else if (e.key === '.') theta += 1;
    else if (e.key === 'a') model.setControls({ sparkAdvanceDeg: (spark = Math.min(64.5, spark + 5)) });
    else if (e.key === 'z') model.setControls({ sparkAdvanceDeg: (spark = Math.max(-15.5, spark - 5)) });
    else if (e.key === 't') model.setControls({ throttle: (throttle = Math.min(1, throttle + 0.1)) });
    else if (e.key === 'g') model.setControls({ throttle: (throttle = Math.max(0, throttle - 0.1)) });
  });
  let last = performance.now();
  const frame = (now: number) => {
    const dt = Math.min(0.05, (now - last) / 1000);
    last = now;
    if (!paused) theta += rpm * 6 * dt;
    model.update(makeSnap(theta));
    controls.update();
    renderer.render(scene, camera);
    hud.textContent = `Model T  θ ${s.thetaDeg.toFixed(1)}°  spark ${spark.toFixed(1)}° BTDC  throttle ${(throttle * 100).toFixed(0)} %  ` +
      `IV1 ${(s.intakeLift * 1000).toFixed(2)} mm  EV1 ${(s.exhaustLift * 1000).toFixed(2)} mm  ` +
      `[space] pause  [c] cutaway  [1-3] views  [a z] spark  [t g] throttle  [, .] step`;
    requestAnimationFrame(frame);
  };
  requestAnimationFrame(frame);
  (window as unknown as { __engine: ModelTEngineModel }).__engine = model;
  (window as unknown as { __view: unknown }).__view = { camera, controls, setView };
}

function runCfr(): void {
  let cr = Number(q.get('cr') ?? 8);
  const model = new EngineModel(CFR_F1, cr);
  scene.add(model.root);
  const L = model.layout;

  // marker at the spark gap (what the combustion agent will draw there)
  const gapMarker = new THREE.Mesh(new THREE.SphereGeometry(CFR_F1.sparkPlug.gap * 0.4, 12, 8), new THREE.MeshBasicMaterial({ color: 0x66ccff }));
  gapMarker.position.set(...CFR_F1.sparkPlug.gapCenter);
  model.cylinderFrame.add(gapMarker);

  function setView(which: string): void {
    if (which === 'back') {
      camera.position.set(-0.9, 0.75, -1.2);
      controls.target.set(0, 0.25, -0.1);
    } else if (which === 'valves') {
      camera.position.set(-0.35, 0.75, 0.45);
      controls.target.set(0, L.headY(cr) + 0.08, -0.06);
    } else {
      const v = recommendedCameraView(CFR_F1, cr, which === 'chamber' ? 'chamber' : 'engine');
      camera.position.copy(v.position);
      controls.target.copy(v.target);
      camera.fov = v.fov;
      camera.updateProjectionMatrix();
    }
    controls.update();
  }
  setView(q.get('view') ?? 'engine');
  // explicit camera in CYLINDER-frame coordinates: ?cam=x,y,z&tgt=x,y,z
  if (q.get('cam') && q.get('tgt')) {
    const c = q.get('cam')!.split(',').map(Number), t = q.get('tgt')!.split(',').map(Number);
    const hy = L.headY(cr);
    camera.position.set(c[0], c[1] + hy, c[2]);
    controls.target.set(t[0], t[1] + hy, t[2]);
    if (q.get('fov')) { camera.fov = Number(q.get('fov')); camera.updateProjectionMatrix(); }
    controls.update();
  }
  model.setCutaway(q.get('cut') !== '0');

  const s = {} as EngineSnapshot;
  function makeSnap(th: number): EngineSnapshot {
    const sc = sliderCrank(L, th);
    s.thetaDeg = wrapDeg720(th);
    s.pistonDisplacement = sc.pistonDisplacement;
    s.rodAngle = sc.rodAngle;
    s.clearanceHeight = L.clearanceAtTdc(cr) + sc.pistonDisplacement;
    s.intakeLift = modelValveLift(CFR_F1.intakeValve, s.thetaDeg);
    s.exhaustLift = modelValveLift(CFR_F1.exhaustValve, s.thetaDeg);
    return s;
  }

  window.addEventListener('keydown', (e: KeyboardEvent) => {
    if (e.key === ' ') paused = !paused;
    else if (e.key === 'c') model.setCutaway(!model.isCutaway);
    else if (e.key === '[') { cr = Math.max(L.crRange[0], cr - 0.5); model.setCompressionRatio(cr); }
    else if (e.key === ']') { cr = Math.min(L.crRange[1], cr + 0.5); model.setCompressionRatio(cr); }
    else if (e.key === '1') setView('engine');
    else if (e.key === '2') setView('chamber');
    else if (e.key === '3') setView('back');
    else if (e.key === '4') setView('valves');
    else if (e.key === ',') theta -= 1;
    else if (e.key === '.') theta += 1;
  });

  let last = performance.now();
  function frame(now: number): void {
    const dt = Math.min(0.05, (now - last) / 1000);
    last = now;
    if (!paused) theta += rpm * 6 * dt;
    model.update(makeSnap(theta));
    controls.update();
    renderer.render(scene, camera);
    hud.textContent = `θ ${s.thetaDeg.toFixed(1)}°  CR ${cr.toFixed(1)}  h ${(s.clearanceHeight * 1000).toFixed(1)} mm  ` +
      `IV ${(s.intakeLift * 1000).toFixed(2)} mm  EV ${(s.exhaustLift * 1000).toFixed(2)} mm  ` +
      `[space] pause  [c] cutaway  [ [ ] ] CR  [1-4] views  [, .] step`;
    requestAnimationFrame(frame);
  }
  requestAnimationFrame(frame);
  (window as unknown as { __engine: EngineModel }).__engine = model;
}
