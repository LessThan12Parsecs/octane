/**
 * DEV ONLY — standalone page for the UI layer driven by MockStream.
 * Open http://localhost:5173/src/ui/dev/ui-harness.html with `npm run dev`.
 * Mimics the SimClient contract loosely: raw snapshots at an adaptive cadence
 * (finer around the spark and knock), played back at timeScale × wall time.
 */
import '../../style.css';
import { CFR_F1 } from '../../physics/engines/cfr';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { UIController } from '../index';
import { MockStream } from './mock-stream';

const op: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: 101325,
  ambientTemperature: 298.15,
  relativeHumidity: 0.3,
  intakeMixtureTemperature: 325.15,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: 1.05,
  sparkAdvanceDeg: 13,
  dwellTime: 3e-3,
  compressionRatio: 7,
  egrFraction: 0,
  coolantTemperature: 373.15,
};

const mock = new MockStream(CFR_F1, op);
let timeScale = 1 / 50;
let paused = false;
let playT = 0;
/** Generated but not yet delivered snapshots. */
const pending: EngineSnapshot[] = [];
/** Delivered history for stepping back. */
const history: EngineSnapshot[] = [];
let current: EngineSnapshot = mock.step(0.5);
let sparkDeg = -op.sparkAdvanceDeg;

function cadenceDeg(theta: number, rpm: number): number {
  if (theta >= sparkDeg - 1 && theta < sparkDeg + 12) return 4e-6 * 6 * rpm; // 4 µs around the spark
  if (theta >= -20 && theta < 60) return 0.05; // resolve knock oscillation
  return 0.5;
}

function generateUntil(t: number): void {
  while (mock.t < t) {
    const s = mock.step(cadenceDeg(mock.theta, 600));
    pending.push(s);
  }
}

const viewport = document.getElementById('viewport')!;
viewport.innerHTML =
  '<div style="position:absolute;inset:0;display:grid;place-items:center;color:#6b7480;font:12px ui-monospace,Menlo,monospace">3D viewport (EngineModel renders here)</div>';

const ui = new UIController(document.getElementById('ui')!, {
  spec: CFR_F1,
  initialOperatingPoint: op,
  onOperatingPointChange: (patch) => {
    mock.setOperatingPoint(patch);
    if (patch.sparkAdvanceDeg !== undefined) sparkDeg = -patch.sparkAdvanceDeg;
    console.debug('[harness] op patch', patch);
  },
  onPlaybackChange: (p) => {
    timeScale = p.timeScale;
    paused = p.paused;
  },
  onStep: (deg) => {
    const dt = deg / (6 * current.rpm);
    playT = Math.max(0, playT + dt);
    advance(0);
  },
  onViewChange: (v) => console.debug('[harness] view', v),
  onReset: () => {
    mock.reset();
    pending.length = 0;
    history.length = 0;
    playT = 0;
  },
});

const recent: EngineSnapshot[] = [];

function advance(dtWall: number): void {
  if (!paused) playT += dtWall * timeScale;
  generateUntil(playT + 1e-3);
  recent.length = 0;
  while (pending.length && pending[0].t <= playT) {
    const s = pending.shift()!;
    recent.push(s);
    history.push(s);
  }
  if (history.length > 40000) history.splice(0, history.length - 40000);
  // Playback snapshot = last sample at or before playT (stepping back searches history).
  let s = recent.length ? recent[recent.length - 1] : current;
  if (s.t > playT) {
    for (let i = history.length - 1; i >= 0; i--) {
      if (history[i].t <= playT) {
        s = history[i];
        break;
      }
    }
  }
  current = s;
  for (const c of mock.drainCycles()) ui.pushCycle(c);
}

// URL knobs: ?warm=<sim s>&ts=<scale>&pause=1&tab=cycles&scroll=<px>&zoom=a,b&perf=1
const q = new URLSearchParams(location.search);
const warm = Number(q.get('warm') ?? 0);
if (warm > 0) {
  // Fast-forward: feed frames of 1/60 s wall at a large time scale.
  const saved = timeScale;
  timeScale = 0.5;
  while (playT < warm - 1e-12) {
    advance(Math.min(1 / 60, (warm - playT) / timeScale));
    ui.update(current, recent);
  }
  timeScale = saved;
}
if (q.get('ts')) ui.setPlayback({ timeScale: Number(q.get('ts')) });
if (q.get('pause') === '1') ui.setPlayback({ paused: true });
if (q.get('tab') === 'cycles') (document.querySelector('.oct-tabs button:nth-child(2)') as HTMLButtonElement | null)?.click();
if (q.get('scroll')) setTimeout(() => ((document.querySelector('.oct-panel-scroll') as HTMLElement).scrollTop = Number(q.get('scroll'))), 50);
if (q.get('zoom')) {
  const [a, b] = q.get('zoom')!.split(',').map(Number);
  (ui as unknown as { traces: { thetaZoom: { set(a: number, b: number): void } } }).traces.thetaZoom.set(a, b);
}

const perfAcc: number[] = [];
if (q.get('perf') === '1') {
  setTimeout(() => {
    const a = perfAcc.slice(5).sort((x, y) => x - y);
    const pre = document.createElement('pre');
    pre.style.cssText = 'position:fixed;left:320px;top:120px;color:#fff;background:#000;font:12px monospace;z-index:99;padding:6px';
    const avg = a.reduce((x, y) => x + y, 0) / a.length;
    pre.textContent = `ui.update over ${a.length} frames: mean ${avg.toFixed(2)} ms, p50 ${a[a.length >> 1]?.toFixed(2)} ms, p95 ${a[Math.floor(a.length * 0.95)]?.toFixed(2)} ms (sync part; uPlot redraw runs in a microtask)`;
    document.body.appendChild(pre);
  }, 2500);
}
let last = performance.now();
function frame(now: number): void {
  const dt = Math.min(0.1, (now - last) / 1000);
  last = now;
  advance(dt);
  const t0 = performance.now();
  ui.update(current, recent);
  if (perfAcc.length < 20000) perfAcc.push(performance.now() - t0);
  requestAnimationFrame(frame);
}
requestAnimationFrame(frame);

(window as unknown as { octaneUi: UIController }).octaneUi = ui;

