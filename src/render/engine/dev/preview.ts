/**
 * Stand-alone visual check of the engine model (open
 * /src/render/engine/dev/preview.html on the Vite dev server). Drives the
 * model with an exact slider-crank + model valve-lift "snapshot", no physics.
 * Keys: space pause · c cutaway · [ ] compression ratio · 1/2/3 views · ,/. step 1°
 * URL params: ?theta=deg&view=engine|chamber|back|valves&cut=0|1&cr=8&paused=1&rpm=60
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { RoomEnvironment } from 'three/examples/jsm/environments/RoomEnvironment.js';
import { CFR_F1 } from '../../../physics/engines/cfr';
import type { EngineSnapshot } from '../../../physics/core/snapshot';
import { EngineModel, recommendedCameraView } from '../index';
import { sliderCrank } from '../kinematics';
import { modelValveLift, wrapDeg720 } from '../cam-profile';

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

let cr = Number(q.get('cr') ?? 8);
const model = new EngineModel(CFR_F1, cr);
scene.add(model.root);
const L = model.layout;

// marker at the spark gap (what the combustion agent will draw there)
const gapMarker = new THREE.Mesh(new THREE.SphereGeometry(CFR_F1.sparkPlug.gap * 0.4, 12, 8), new THREE.MeshBasicMaterial({ color: 0x66ccff }));
gapMarker.position.set(...CFR_F1.sparkPlug.gapCenter);
model.cylinderFrame.add(gapMarker);

const camera = new THREE.PerspectiveCamera(34, innerWidth / innerHeight, 0.005, 30);
const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
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

let theta = Number(q.get('theta') ?? -90);
let paused = q.get('paused') === '1';
const rpm = Number(q.get('rpm') ?? 30);
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
window.addEventListener('resize', () => {
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
  renderer.setSize(innerWidth, innerHeight);
});

const hud = document.getElementById('hud')!;
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
