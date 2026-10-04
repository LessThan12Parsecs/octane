/**
 * DEV ONLY — L-head / multi-cylinder check page for the in-cylinder visuals: the Ford Model T spec
 * driven by SYNTHETIC snapshots (no physics — scripted, for looking at shapes): slider-crank kinematics,
 * sin² valve lifts, a scripted flame from the plug over the valves and a trembler spark shower
 * (cumulative breakdownCount) per cylinder at its firing offset, packed as cylinders[] like the worker.
 * Open http://localhost:<port>/src/render/combustion/dev/lhead-dev.html with `npm run dev`.
 * Query: ?mode=temperature  ?ts=<sim s per wall s>  ?theta=<deg>  ?paused  ?cut (cut region on).
 * Keys: m mode, space pause, . step.
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import { cylinderAngleDeg } from '../../../physics/core/engine-spec';
import type { CylinderSnapshot, EngineSnapshot } from '../../../physics/core/snapshot';
import { MODEL_T } from '../../../physics/engines/model-t';
import { chamberShapeOf } from '../chamber';
import { CombustionVisuals } from '../index';

const spec = MODEL_T;
const g = spec.geometry;
const R = g.bore / 2;
const shape = chamberShapeOf(spec);
const params = new URLSearchParams(location.search);
const rpm = 1000;
let timeScale = Number(params.get('ts') ?? 2e-3);
let paused = params.has('paused');
let mode: 'physical' | 'temperature' = (params.get('mode') as 'physical' | 'temperature') ?? 'physical';

// ------------------------------------------------------------ synthetic snapshots
const SPARK = -25; // first breakdown, deg
const BURN = 70; // burn duration, deg
const smooth = (u: number): number => (u <= 0 ? 0 : u >= 1 ? 1 : u * u * (3 - 2 * u));
function lift(th: number, open: number, close: number, max: number): number {
  let span = close - open;
  if (span <= 0) span += 720;
  let u = th - open;
  if (u < 0) u += 720;
  return u < span ? max * Math.sin((Math.PI * u) / span) ** 2 : 0;
}
function cylinderState(th: number, t: number, cycle: number, index: number): CylinderSnapshot {
  const r = g.stroke / 2, L = g.conRodLength, a = (th * Math.PI) / 180;
  const x = r * (1 - Math.cos(a)) + L - Math.sqrt(L * L - (r * Math.sin(a)) ** 2);
  const h = shape.depthTDC + x;
  const il = lift(th, spec.intakeValve.openDeg, spec.intakeValve.closeDeg, spec.intakeValve.maxLift);
  const el = lift(th, spec.exhaustValve.openDeg, spec.exhaustValve.closeDeg, spec.exhaustValve.maxLift);
  const u = (th - SPARK - 4) / BURN;
  const burning = u > 0 && u < 1;
  const xb = smooth(u);
  const rf = burning ? 0.002 + 0.14 * Math.pow(Math.max(u, 0), 0.8) : 0;
  const nbd = th >= SPARK && th < SPARK + 40 ? Math.floor((th - SPARK) / 2.5) + 1 : th >= SPARK + 40 ? 16 : 0;
  const sparkOn = th >= SPARK && th < SPARK + 40 && (th - SPARK) % 2.5 < 0.6;
  const gx = il > 0 || el > 0;
  return {
    index, thetaDeg: th, cycle, gasTorque: 0,
    pistonDisplacement: x, clearanceHeight: h, rodAngle: 0, intakeLift: il, exhaustLift: el,
    phase: gx ? 'gas-exchange' : th < SPARK ? 'compression' : burning ? 'combustion' : 'expansion',
    volume: Math.PI * R * R * h + shape.pocketVolume, pressure: 6e5 + 1.2e6 * xb * (burning ? 1 : 0),
    temperatureMean: 600 + 1200 * xb, temperatureUnburned: 650, temperatureBurned: xb > 0 ? 2300 : 0,
    massFractionBurned: gx ? 0 : th > SPARK ? xb : 0, mass: 4e-4, heatReleaseRate: burning ? 4e4 * Math.sin(Math.PI * u) : 0,
    heatLossRate: 0,
    flame: {
      stage: burning ? (u < 0.05 ? 'kernel' : u < 0.8 ? 'turbulent' : 'burnout') : th > SPARK + 4 && !gx ? 'done' : 'none',
      radius: rf, center: [...spec.sparkPlug.gapCenter] as [number, number, number], area: 4 * Math.PI * rf * rf * 0.4,
      laminarSpeed: 0.4, turbulentSpeed: 2, turbulenceIntensity: 1.5,
    },
    spark: {
      phase: sparkOn ? 'arc' : nbd > 0 && th < SPARK + 40 ? 'glow' : 'off', primaryCurrent: 3, secondaryVoltage: sparkOn ? 300 : 0,
      secondaryCurrent: sparkOn ? 0.05 : 0, energyDelivered: nbd * 2e-3, breakdownVoltage: 8e3, breakdownCount: nbd,
    },
    intakeMassFlow: il > 0 ? 0.02 * (il / spec.intakeValve.maxLift) : 0,
    exhaustMassFlow: el > 0 ? 0.03 * (el / spec.exhaustValve.maxLift) : 0,
    knock: { integral: 0, autoignited: false, oscillation: 0 },
    burnedComposition: { CO2: 0.12, H2O: 0.13, CO: 0.01, O2: 0.005, H2: 0.002, OH: 0.002, H: 0.0005, O: 0.0003, NO: 0.002, N2: 0.72 },
  };
}
const cylStates: CylinderSnapshot[] = [];
function snapshotAt(t: number): EngineSnapshot {
  const deg = (t * rpm * 6) % 720;
  const theta = deg - 360;
  const engCycle = Math.floor((t * rpm * 6) / 720);
  for (let i = 0; i < spec.cylinders; i++) {
    const th = cylinderAngleDeg(spec, i, theta);
    cylStates[i] = cylinderState(th, t, engCycle, i);
  }
  const c0 = cylStates[0];
  return {
    ...c0, t, cycle: engCycle, thetaDeg: theta, rpm, intakeManifoldPressure: 0.7e5, exhaustManifoldPressure: 1.05e5,
    gasTorque: 0, netTorque: 0, cylinders: cylStates.slice(),
  } as EngineSnapshot;
}

// ------------------------------------------------------------ scene
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.outputColorSpace = THREE.SRGBColorSpace;
document.body.prepend(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0d10);
const camera = new THREE.PerspectiveCamera(40, innerWidth / innerHeight, 0.005, 10);
const controls = new OrbitControls(camera, renderer.domElement);
camera.position.set(-0.32, 0.42, 0.42);
controls.target.set(-0.03, 0.27, 0.12);
controls.update();
scene.add(new THREE.HemisphereLight(0xbfd4ff, 0x202020, 0.4));
const sun = new THREE.DirectionalLight(0xffffff, 0.7);
sun.position.set(-0.3, 0.8, 0.5);
scene.add(sun);
const metal = new THREE.MeshStandardMaterial({ color: 0x8a8f96, metalness: 0.4, roughness: 0.6 });
const lineMat = new THREE.LineBasicMaterial({ color: 0x5c7080 });
const layout = spec.layout;
const vis: CombustionVisuals[] = [];
const pistons: THREE.Mesh[] = [];
for (let i = 0; i < spec.cylinders; i++) {
  const frame = new THREE.Object3D();
  frame.position.set(0, 0.3, layout.axisZ[i]);
  if (layout.mirrorZ[i]) frame.scale.z = -1;
  scene.add(frame);
  const v = new CombustionVisuals(spec, { cylinder: i, tracers: i === 0, light: i === 0 });
  v.setMode(mode);
  v.setElectrodesVisible(true);
  if (params.has('cut')) v.setCutRegion([[-1, 0, 0, 0], [0, 0, 1, 0]]); // quadrant x < 0, z > 0 removed
  frame.add(v.root);
  vis.push(v);
  // stand-ins: piston, outlines of the bore circle and pocket at the deck and roof
  const piston = new THREE.Mesh(new THREE.CylinderGeometry(R * 0.995, R * 0.995, 0.06, 64), metal);
  frame.add(piston);
  pistons.push(piston);
  for (const y of [shape.deckY, shape.roofY, 0]) {
    const pts: THREE.Vector3[] = [];
    for (let k = 0; k <= 96; k++) pts.push(new THREE.Vector3(R * Math.cos((k / 48) * Math.PI), y, R * Math.sin((k / 48) * Math.PI)));
    frame.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), lineMat));
    if (y === 0) continue;
    const rr: THREE.Vector3[] = [];
    const cx = [shape.ix1, shape.ix0, shape.ix0, shape.ix1], cz = [shape.iz1, shape.iz1, shape.iz0, shape.iz0];
    for (let c = 0; c < 4; c++) for (let k = 0; k <= 8; k++) {
      const a = ((c + k / 8) * Math.PI) / 2;
      rr.push(new THREE.Vector3(cx[c] + shape.rc * Math.cos(a), y, cz[c] + shape.rc * Math.sin(a)));
    }
    rr.push(rr[0].clone());
    frame.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(rr), lineMat));
  }
}

// ------------------------------------------------------------ loop
let simT = 0;
if (params.has('theta')) simT = ((Number(params.get('theta')) + 360) / (rpm * 6));
const hud = document.getElementById('hud')!;
window.addEventListener('keydown', (e: KeyboardEvent) => {
  if (e.key === 'm') { mode = mode === 'physical' ? 'temperature' : 'physical'; for (const v of vis) v.setMode(mode); }
  if (e.key === ' ') paused = !paused;
  if (e.key === '.') simT += 0.25 / (6 * rpm);
});
let last = performance.now();
let frame = 0;
function tick(now: number): void {
  const dtWall = Math.min((now - last) / 1000, 0.1);
  last = now;
  const ts = paused ? 0 : timeScale;
  simT += dtWall * ts;
  const s = snapshotAt(simT);
  for (let i = 0; i < vis.length; i++) {
    vis[i].update(s, dtWall, ts);
    pistons[i].position.y = -s.cylinders![i].clearanceHeight - 0.03;
  }
  renderer.render(scene, camera);
  if ((frame++ & 7) === 0) {
    hud.textContent = [
      `θ ${s.thetaDeg.toFixed(1)}°  mode ${mode}  (m: mode, space: pause, .: step)`,
      ...s.cylinders!.map((c, i) => `#${i + 1} θ ${c.thetaDeg.toFixed(1)}  ${c.phase.padEnd(12)} flame ${c.flame.stage} r ${(c.flame.radius * 1e3).toFixed(1)} mm  sparks ${c.spark.breakdownCount}  tracers ${vis[i].tracers.count}  h ${(vis[i].state.h * 1e3).toFixed(1)} mm`),
    ].join('\n');
  }
  requestAnimationFrame(tick);
}
addEventListener('resize', () => {
  renderer.setSize(innerWidth, innerHeight);
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
});
requestAnimationFrame(tick);
function look(pos: [number, number, number], target: [number, number, number]): void {
  camera.position.set(...pos);
  controls.target.set(...target);
  controls.update();
}
Object.assign(window as unknown as Record<string, unknown>, {
  vis, look, camera, controls, renderer, scene, THREE,
  setTheta: (th: number) => { simT = (th + 360) / (rpm * 6); },
  setTimeScale: (x: number) => { timeScale = x; },
  setPaused: (p: boolean) => { paused = p; },
});
