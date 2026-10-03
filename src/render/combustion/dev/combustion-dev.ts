/**
 * DEV ONLY — standalone page for the in-cylinder visuals driven by the
 * MockSimulator (main thread, simple linear-interpolating playback).
 * Open http://localhost:<port>/src/render/combustion/dev/combustion-dev.html with `npm run dev`.
 * The engine here is a crude stand-in (half liner, head, piston, valve discs);
 * the real mechanism is EngineModel's job.
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';
import GUI from 'lil-gui';
import { CFR_F1 } from '../../../physics/engines/cfr';
import type { OperatingPoint } from '../../../physics/core/operating-point';
import type { EngineSnapshot } from '../../../physics/core/snapshot';
import { MockSimulator } from '../../../worker/mock-simulator';
import { EngineModel, recommendedCameraView } from '../../engine/index';
import { CombustionVisuals, createTemperatureLegendElement } from '../index';

const spec = CFR_F1;
const R = spec.geometry.bore / 2;
const op: OperatingPoint = {
  speedMode: 'fixed', rpm: 600, loadTorque: 0, throttle: 1, ambientPressure: 101325,
  ambientTemperature: 298.15, relativeHumidity: 0.3, intakeMixtureTemperature: 325.15,
  fuel: { kind: 'PRF', octaneNumber: 90 }, equivalenceRatio: 1.05, sparkAdvanceDeg: 13,
  dwellTime: 3e-3, compressionRatio: 7, egrFraction: 0, coolantTemperature: 373.15,
};
const params = new URLSearchParams(location.search);
if (params.has('on')) op.fuel = { kind: 'PRF', octaneNumber: Number(params.get('on')) };
if (params.has('cr')) op.compressionRatio = Number(params.get('cr'));

// ---------------------------------------------------------------- playback
const sim = new MockSimulator(spec, op, { snapshotEveryDeg: 0.2, bufferAheadSeconds: 0.05, cyclicVariability: 1 });
const buf: EngineSnapshot[] = [sim.advanceToNextSnapshot()];
let playT = 0;
const ui = {
  timeScale: Number(params.get('ts') ?? 2e-3),
  paused: params.has('paused'),
  mode: (params.get('mode') ?? 'physical') as 'physical' | 'temperature',
  tracers: true,
  chamberLight: true,
  electrodes: true,
  octane: op.fuel.kind === 'PRF' ? op.fuel.octaneNumber : 90,
  compressionRatio: op.compressionRatio,
  sparkAdvance: op.sparkAdvanceDeg,
  rpm: op.rpm,
  jumpToSpark: () => jumpTo(-ui.sparkAdvance - 1.5),
  jumpToIntake: () => jumpTo(-330),
  jumpToTDC: () => jumpTo(8),
  step: () => { playT += (0.25 / (6 * ui.rpm)); },
};

function ensure(t: number): void {
  while (buf[buf.length - 1].t < t) buf.push(sim.advanceToNextSnapshot());
  while (buf.length > 2 && buf[1].t < t - 1e-3) buf.shift();
}

function jumpTo(theta: number): void {
  // play forward until the next time the crank passes theta
  let guard = 0;
  let last = sample(playT).thetaDeg;
  for (;;) {
    playT += 0.2 / (6 * ui.rpm);
    const th = sample(playT).thetaDeg;
    if ((last < theta && th >= theta) || guard++ > 20000) break;
    last = th;
  }
}

function lerp(a: number, b: number, w: number): number { return a + (b - a) * w; }
function interp(a: EngineSnapshot, b: EngineSnapshot, w: number): EngineSnapshot {
  const o = structuredClone(w < 0.5 ? a : b) as unknown as Record<string, unknown>;
  const A = a as unknown as Record<string, unknown>, B = b as unknown as Record<string, unknown>;
  for (const k of Object.keys(o)) {
    const x = A[k], y = B[k];
    if (typeof x === 'number' && typeof y === 'number' && k !== 'thetaDeg' && k !== 'cycle') o[k] = lerp(x, y, w);
    else if (x && y && typeof x === 'object' && !Array.isArray(x)) {
      const ox = o[k] as Record<string, unknown>;
      for (const kk of Object.keys(ox)) {
        const xx = (x as Record<string, unknown>)[kk], yy = (y as Record<string, unknown>)[kk];
        if (typeof xx === 'number' && typeof yy === 'number') ox[kk] = lerp(xx, yy, w);
        else if (Array.isArray(xx) && Array.isArray(yy)) ox[kk] = xx.map((v: number, i: number) => lerp(v, yy[i] as number, w));
      }
    }
  }
  let dth = b.thetaDeg - a.thetaDeg;
  if (dth < -360) dth += 720;
  let th = a.thetaDeg + dth * w;
  if (th >= 360) th -= 720;
  o.thetaDeg = th;
  return o as unknown as EngineSnapshot;
}
function sample(t: number): EngineSnapshot {
  ensure(t);
  let i = 0;
  while (i < buf.length - 2 && buf[i + 1].t <= t) i++;
  const a = buf[i], b = buf[i + 1] ?? a;
  const w = b.t > a.t ? Math.min(Math.max((t - a.t) / (b.t - a.t), 0), 1) : 0;
  return interp(a, b, w);
}

// ---------------------------------------------------------------- scene
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
renderer.setSize(innerWidth, innerHeight);
renderer.toneMapping = THREE.ACESFilmicToneMapping;
renderer.outputColorSpace = THREE.SRGBColorSpace;
document.body.prepend(renderer.domElement);
const scene = new THREE.Scene();
scene.background = new THREE.Color(0x0b0d10);
const camera = new THREE.PerspectiveCamera(40, innerWidth / innerHeight, 0.002, 5);
camera.position.set(0.02, 0.33, 0.17);
const controls = new OrbitControls(camera, renderer.domElement);
controls.target.set(0, 0.29, 0);
controls.update();
scene.add(new THREE.HemisphereLight(0xbfd4ff, 0x202020, 0.35));
const sun = new THREE.DirectionalLight(0xffffff, 0.6);
sun.position.set(0.3, 0.8, 0.5);
scene.add(sun);

// cylinder frame stand-in: head fire deck at world y = 0.3
const cylFrame = new THREE.Object3D();
cylFrame.position.set(0, 0.3, 0);
scene.add(cylFrame);
const metal = new THREE.MeshStandardMaterial({ color: 0x8a8f96, metalness: 0.5, roughness: 0.55, side: THREE.DoubleSide });
// back half of the liner (cutaway)
const liner = new THREE.Mesh(new THREE.CylinderGeometry(R * 1.0005, R * 1.0005, 0.14, 64, 1, true, Math.PI / 2, Math.PI), metal);
liner.position.y = -0.07;
cylFrame.add(liner);
const head = new THREE.Mesh(new THREE.CylinderGeometry(R * 1.3, R * 1.3, 0.03, 64), metal);
head.position.y = 0.015 + 1e-4;
cylFrame.add(head);
const piston = new THREE.Mesh(new THREE.CylinderGeometry(R * 0.999, R * 0.999, 0.06, 64), metal);
cylFrame.add(piston);
const mkValve = (v: typeof spec.intakeValve): THREE.Mesh => {
  const m = new THREE.Mesh(new THREE.CylinderGeometry(v.headDiameter / 2, v.headDiameter / 2 * 0.8, 0.12 * v.headDiameter, 48), metal);
  m.position.set(v.position[0], 0, v.position[1]);
  cylFrame.add(m);
  return m;
};
const iv = mkValve(spec.intakeValve);
const ev = mkValve(spec.exhaustValve);

// ?engine=1 → mount on the real EngineModel instead of the stand-in parts
const useEngine = params.has('engine');
let model: EngineModel | null = null;
const vis = new CombustionVisuals(spec);
if (useEngine) {
  model = new EngineModel(spec, op.compressionRatio);
  scene.add(model.root);
  model.cylinderFrame.add(vis.root);
  cylFrame.visible = false;
  const v = recommendedCameraView(spec, op.compressionRatio, 'chamber');
  camera.fov = v.fov; camera.near = v.near; camera.far = v.far; camera.updateProjectionMatrix();
  camera.position.copy(v.position);
  controls.target.copy(v.target);
  controls.update();
  ui.electrodes = false;
} else {
  cylFrame.add(vis.root);
  vis.setElectrodesVisible(true);
}
vis.setMode(ui.mode);

const legend = createTemperatureLegendElement();
document.getElementById('legend')!.append(legend.element);
legend.element.style.display = ui.mode === 'temperature' ? '' : 'none';

// ---------------------------------------------------------------- gui
const gui = new GUI({ title: 'combustion dev' });
gui.add(ui, 'timeScale', 1e-5, 1, 1e-5).name('time scale (sim s / wall s)');
gui.add(ui, 'paused');
gui.add(ui, 'step').name('step 0.25°');
gui.add(ui, 'mode', ['physical', 'temperature']).onChange((m: 'physical' | 'temperature') => {
  vis.setMode(m);
  legend.element.style.display = m === 'temperature' ? '' : 'none';
});
gui.add(ui, 'tracers').onChange((v: boolean) => vis.setTracersVisible(v));
gui.add(ui, 'chamberLight').onChange((v: boolean) => vis.setChamberLightEnabled(v));
gui.add(ui, 'electrodes').onChange((v: boolean) => vis.setElectrodesVisible(v));
gui.add(ui, 'octane', 0, 120, 1).onFinishChange((v: number) => sim.setOperatingPoint({ fuel: { kind: 'PRF', octaneNumber: v } }));
gui.add(ui, 'compressionRatio', 4, 12, 0.1).onFinishChange((v: number) => {
  sim.setOperatingPoint({ compressionRatio: v });
  model?.setCompressionRatio(v);
});
gui.add(ui, 'sparkAdvance', -10, 40, 0.5).onFinishChange((v: number) => sim.setOperatingPoint({ sparkAdvanceDeg: v }));
gui.add(ui, 'rpm', 300, 1800, 10).onFinishChange((v: number) => sim.setOperatingPoint({ rpm: v }));
gui.add(ui, 'jumpToIntake');
gui.add(ui, 'jumpToSpark');
gui.add(ui, 'jumpToTDC');

// ---------------------------------------------------------------- loop
const hud = document.getElementById('hud')!;
let last = performance.now();
let frame = 0;
function tick(now: number): void {
  const dtWall = Math.min((now - last) / 1000, 0.1);
  last = now;
  const ts = ui.paused ? 0 : ui.timeScale;
  playT += dtWall * ts;
  const s = sample(playT);
  if (model) model.update(s);
  piston.position.y = -s.clearanceHeight - 0.03;
  iv.position.y = -s.intakeLift + 0.06 * spec.intakeValve.headDiameter;
  ev.position.y = -s.exhaustLift + 0.06 * spec.exhaustValve.headDiameter;
  vis.update(s, dtWall, ts);
  renderer.render(scene, camera);
  if ((frame++ & 3) === 0) {
    const st = vis.state;
    hud.textContent = [
      `t ${(s.t * 1e3).toFixed(3)} ms  θ ${s.thetaDeg.toFixed(2)}°  cycle ${s.cycle}  ${s.phase}  p ${(s.pressure / 1e5).toFixed(2)} bar  h ${(s.clearanceHeight * 1e3).toFixed(1)} mm`,
      `flame ${s.flame.stage} r ${(s.flame.radius * 1e3).toFixed(2)} mm  A ${(s.flame.area * 1e4).toFixed(2)} cm²  SL ${s.flame.laminarSpeed.toFixed(2)}  u' ${s.flame.turbulenceIntensity.toFixed(2)}  xb ${s.massFractionBurned.toFixed(3)}  HRR ${(s.heatReleaseRate / 1e3).toFixed(1)} kW`,
      `Tu ${s.temperatureUnburned.toFixed(0)} K  Tb ${s.temperatureBurned.toFixed(0)} K  q'' ${(st.frontFlux / 1e6).toFixed(1)} MW/m²  δ ${(st.brush.thickness * 1e3).toFixed(2)} mm  a ${(st.brush.wrinkleAmplitude * 1e3).toFixed(2)} mm  φ≈${st.phi.toFixed(2)}`,
      `spark ${s.spark.phase}  V ${s.spark.secondaryVoltage.toFixed(0)}  I ${(s.spark.secondaryCurrent * 1e3).toFixed(1)} mA  E ${(s.spark.energyDelivered * 1e3).toFixed(2)} mJ  Vbd ${(s.spark.breakdownVoltage / 1e3).toFixed(1)} kV  P ${st.spark.power.toFixed(1)} W  ${st.spark.kind}`,
      `knock ${s.knock.autoignited ? 'AUTOIGNITED' : '-'}  LW ${s.knock.integral.toFixed(2)}  osc ${(s.knock.oscillation / 1e5).toFixed(2)} bar  env ${(st.knock.envelopePa / 1e5).toFixed(2)} bar  f ${(st.knock.frequency / 1e3).toFixed(2)} kHz  ${st.knock.rms ? 'rms' : 'resolved'}`,
      `ṁi ${(s.intakeMassFlow * 1e3).toFixed(2)} g/s  ṁe ${(s.exhaustMassFlow * 1e3).toFixed(2)} g/s  swirl ω ${st.flow.omega.toFixed(1)} rad/s  tracers ${vis.tracers.count}  light ${st.light.intensity.toExponential(2)}`,
    ].join('\n');
  }
  requestAnimationFrame(tick);
}
addEventListener('resize', () => {
  renderer.setSize(innerWidth, innerHeight);
  camera.aspect = innerWidth / innerHeight;
  camera.updateProjectionMatrix();
});
if (params.has('theta')) jumpTo(Number(params.get('theta')));
requestAnimationFrame(tick);
// debugging hooks
function look(pos: [number, number, number], target: [number, number, number]): void {
  camera.position.set(...pos);
  controls.target.set(...target);
  controls.update();
}
Object.assign(window as unknown as Record<string, unknown>, { vis, sim, ui, jumpTo, look, THREE, renderer, scene, camera, controls, model });
