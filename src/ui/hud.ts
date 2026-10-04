/**
 * Bottom HUD strip: 720° cycle dial (strokes, valve events, spark timing, crank
 * needle), instantaneous state readouts, knock lamp and the playback transport.
 *
 * Multi-cylinder engines: the dial and the in-cylinder readouts follow the FOCUS
 * cylinder (its local angle), a firing-order strip shows every cylinder's stroke
 * (click a cylinder to focus it), and free-speed / vehicle engines get brake
 * torque-power and road-speed cells.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { CycleSummary, EngineCycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { firingStrip, type FiringStripEntry } from './cylinder-view';
import { STROKE_COLOR, SERIES, STATUS } from './charts/theme';
import {
  cycleDialAngleDeg,
  cycleDialPoint,
  describeCrankAngle,
  STROKE_LABEL,
  STROKES,
  strokeOf,
  strokeStartDeg,
  valveOpenInterval,
} from './engine-cycle';
import { GatedMean } from './cycle-mean';
import { formatTimeScale, STEP_LARGE_DEG, STEP_SMALL_DEG } from './playback';
import { fmt, fmtSimTime, fractionToPct, mpsToKmh, mpsToMph, paToBar, wattsToHp } from './units';

const SVG_NS = 'http://www.w3.org/2000/svg';
const C = 40; // dial centre
const R_STROKE = 34;
const R_IN = 28.5;
const R_EX = 24.5;

export interface HudCallbacks {
  onTogglePause: () => void;
  onStep: (deg: number) => void;
  onFaster: () => void;
  onSlower: () => void;
  /** A cylinder of the firing-order strip was clicked (0-based). */
  onFocusCylinder?: (index: number) => void;
}

export interface HudOptions {
  /** Initial focus cylinder, 0-based. */
  focus?: number;
  /** Show the brake torque / power cell (free-speed or multi-cylinder engines). */
  brake?: boolean;
  /** Show the road-speed cell (engines with a vehicle). */
  vehicle?: boolean;
}

const STROKE_INITIAL: Readonly<Record<string, string>> = { intake: 'I', compression: 'C', power: 'P', exhaust: 'E' };

/** Text node that only touches the DOM when its value changes. */
class Slot {
  private v = '';
  constructor(readonly el: HTMLElement) {}
  set(v: string): void {
    if (v !== this.v) {
      this.v = v;
      this.el.textContent = v;
    }
  }
}

function svg<K extends keyof SVGElementTagNameMap>(tag: K, attrs: Record<string, string | number>): SVGElementTagNameMap[K] {
  const el = document.createElementNS(SVG_NS, tag);
  for (const [k, v] of Object.entries(attrs)) el.setAttribute(k, String(v));
  return el;
}

/** SVG arc path on the cycle dial from θ0 to θ1 (θ1 > θ0, cycle degrees). */
export function dialArcPath(theta0: number, theta1: number, r: number, cx = C, cy = C): string {
  const span = Math.min(719.9, theta1 - theta0) / 2; // dial degrees
  const [x0, y0] = cycleDialPoint(theta0, cx, cy, r);
  const a1 = ((cycleDialAngleDeg(theta0) + span) * Math.PI) / 180;
  const x1 = cx + r * Math.sin(a1);
  const y1 = cy - r * Math.cos(a1);
  const large = span > 180 ? 1 : 0;
  return `M ${x0.toFixed(2)} ${y0.toFixed(2)} A ${r} ${r} 0 ${large} 1 ${x1.toFixed(2)} ${y1.toFixed(2)}`;
}

export class Hud {
  readonly el: HTMLElement;
  private readonly needle: SVGGElement;
  private readonly sparkTick: SVGLineElement;
  private readonly slots: Record<string, Slot> = {};
  private readonly xbBar: HTMLElement;
  private readonly lamp: HTMLElement;
  private readonly sparkDot: HTMLElement;
  private readonly playBtn: HTMLButtonElement;
  private readonly strokeEl: HTMLElement;
  private lastNeedle = NaN;
  private lastLamp = '';
  private lastSparkPhase = '';
  private lastStroke = '';
  private lastXb = -1;
  private readonly mapMean = new GatedMean();
  private readonly strip: { el: HTMLButtonElement; stroke: HTMLElement; key: string }[] = [];
  private readonly stripEntries: FiringStripEntry[] = [];
  private focus = 0;
  private engineCycle: EngineCycleSummary | null = null;

  constructor(
    host: HTMLElement,
    private readonly spec: EngineSpec,
    cb: HudCallbacks,
    opts: HudOptions = {},
  ) {
    this.el = document.createElement('footer');
    this.el.className = 'oct-hud';
    host.appendChild(this.el);
    this.focus = Math.max(0, opts.focus ?? 0);

    // ---- dial ----
    const dialWrap = document.createElement('div');
    dialWrap.className = 'oct-hud-dial';
    dialWrap.title = '720° cycle dial: firing TDC at the top. Outer ring: strokes; inner rings: intake (blue) and exhaust (magenta) valve-open periods; yellow tick: spark.';
    const s = svg('svg', { viewBox: '0 0 80 80', width: 72, height: 72, 'aria-hidden': 'true' });
    s.appendChild(svg('circle', { cx: C, cy: C, r: R_STROKE + 3.5, fill: 'rgba(255,255,255,0.03)', stroke: 'rgba(255,255,255,0.08)' }));
    for (const st of STROKES) {
      const a = strokeStartDeg(st);
      s.appendChild(svg('path', { d: dialArcPath(a + 1, a + 179, R_STROKE), stroke: STROKE_COLOR[st], 'stroke-width': 3, fill: 'none', 'stroke-linecap': 'butt' }));
    }
    const [io, ic] = valveOpenInterval(spec.intakeValve);
    const [eo, ec] = valveOpenInterval(spec.exhaustValve);
    s.appendChild(svg('path', { d: dialArcPath(io, ic, R_IN), stroke: STROKE_COLOR.intake, 'stroke-opacity': 0.55, 'stroke-width': 2.5, fill: 'none' }));
    s.appendChild(svg('path', { d: dialArcPath(eo, ec, R_EX), stroke: STROKE_COLOR.exhaust, 'stroke-opacity': 0.55, 'stroke-width': 2.5, fill: 'none' }));
    // TDC / BDC ticks
    for (const [th, label] of [
      [0, 'TDC'],
      [-180, ''],
      [180, ''],
      [-360, ''],
    ] as const) {
      const [x0, y0] = cycleDialPoint(th, C, C, R_STROKE + 5);
      const [x1, y1] = cycleDialPoint(th, C, C, R_STROKE - 5);
      s.appendChild(svg('line', { x1: x0, y1: y0, x2: x1, y2: y1, stroke: 'rgba(230,233,238,0.55)', 'stroke-width': th === 0 ? 1.5 : 1 }));
      if (label) {
        const t = svg('text', { x: C, y: 5.5, 'text-anchor': 'middle', 'font-size': 6, fill: 'rgba(230,233,238,0.7)', 'font-family': 'ui-monospace, Menlo, monospace' });
        t.textContent = label;
        s.appendChild(t);
      }
    }
    this.sparkTick = svg('line', { x1: C, y1: C, x2: C, y2: C, stroke: SERIES.yellow, 'stroke-width': 2, 'stroke-linecap': 'round' });
    s.appendChild(this.sparkTick);
    this.needle = svg('g', {});
    this.needle.appendChild(svg('line', { x1: C, y1: C + 5, x2: C, y2: C - R_STROKE + 2, stroke: '#e6e9ee', 'stroke-width': 1.6, 'stroke-linecap': 'round' }));
    this.needle.appendChild(svg('circle', { cx: C, cy: C - R_STROKE + 2, r: 2.2, fill: '#e6e9ee' }));
    s.appendChild(this.needle);
    s.appendChild(svg('circle', { cx: C, cy: C, r: 2.5, fill: '#12161c', stroke: '#e6e9ee', 'stroke-width': 1 }));
    dialWrap.appendChild(s);
    this.el.appendChild(dialWrap);

    // ---- crank angle ----
    const th = document.createElement('div');
    th.className = 'oct-hud-theta';
    th.innerHTML = `<div class="oct-hud-theta-val"><span data-k="theta">—</span><span class="oct-unit">°CA</span></div><div class="oct-hud-stroke"><span class="oct-stroke-name" data-k="stroke"></span><span class="oct-hud-rel" data-k="rel"></span></div>`;
    this.el.appendChild(th);
    this.strokeEl = th.querySelector('.oct-stroke-name') as HTMLElement;

    // ---- firing-order strip (multi-cylinder) ----
    if (spec.cylinders > 1) {
      const strip = document.createElement('div');
      strip.className = 'oct-hud-strip';
      strip.setAttribute('role', 'group');
      strip.setAttribute('aria-label', 'Cylinders in firing order');
      const cap = document.createElement('div');
      cap.className = 'oct-hud-label';
      cap.textContent = 'Firing order';
      const row = document.createElement('div');
      row.className = 'oct-hud-strip-row';
      for (const e of firingStrip(spec, makeIdleSnapshot(), this.stripEntries)) {
        const b = document.createElement('button');
        b.type = 'button';
        b.className = 'oct-cyl-chip';
        b.dataset.cyl = String(e.index);
        b.title = `Cylinder ${e.number}: click to focus the readouts, charts and close-up on it`;
        b.innerHTML = `<span class="oct-cyl-n"></span><span class="oct-cyl-s"></span>`;
        (b.firstElementChild as HTMLElement).textContent = String(e.number);
        const idx = e.index;
        b.addEventListener('click', () => {
          cb.onFocusCylinder?.(idx);
          b.blur();
        });
        row.appendChild(b);
        this.strip.push({ el: b, stroke: b.lastElementChild as HTMLElement, key: '' });
      }
      strip.append(cap, row);
      this.el.appendChild(strip);
    }

    // ---- readouts ----
    const cells = document.createElement('div');
    cells.className = 'oct-hud-cells';
    const cell = (k: string, label: string, unit: string, hint: string, extra = ''): void => {
      const d = document.createElement('div');
      d.className = `oct-hud-cell oct-hud-${k}`;
      d.title = hint;
      d.innerHTML = `<div class="oct-hud-label"></div><div class="oct-hud-value"><span data-k="${k}">—</span>${unit ? `<span class="oct-unit">${unit}</span>` : ''}</div>${extra}`;
      (d.firstElementChild as HTMLElement).textContent = label;
      cells.appendChild(d);
    };
    cell('rpm', 'Speed', 'rpm', 'Crankshaft speed');
    cell('map', 'Intake', 'bar', 'Intake manifold pressure averaged over the last intake-valve-open period — the pressure the cylinder breathes, set by the throttle');
    cell('last', 'Last cycle', '', 'Net IMEP and indicated efficiency of the last complete cycle');
    cell('p', 'Pressure', 'bar', 'Cylinder pressure (incl. knock oscillation)');
    cell('T', 'T unb / burned', 'K', 'Unburned- and burned-zone temperatures');
    cell('xb', 'Burned', '%', 'Burned mass fraction', '<div class="oct-hud-bar"><div class="oct-hud-bar-fill"></div></div>');
    cell('flame', 'Flame', '', 'Flame stage: kernel → turbulent → burn-out');
    cell('spark', 'Spark', '', 'Ignition phase: charging (dwell) → breakdown → arc → glow', '');
    cell('knock', 'Knock', '', 'End-gas Livengood–Wu integral (autoignition at 100 %) — lamp lights on autoignition');
    if (opts.brake) cell('brake', 'Brake', '', 'Brake torque and power at the crankshaft (engine output, excludes the car’s inertia): mean of the last engine cycle (instantaneous load torque until one is complete)');
    if (opts.vehicle) cell('road', 'Road speed', '', 'Vehicle speed (vehicle load model)');
    cell('cycle', 'Cycle', '', 'Completed cycles and simulated time');
    this.el.appendChild(cells);
    this.xbBar = cells.querySelector('.oct-hud-bar-fill') as HTMLElement;
    const knockVal = cells.querySelector('.oct-hud-knock .oct-hud-value') as HTMLElement;
    this.lamp = document.createElement('span');
    this.lamp.className = 'oct-lamp';
    knockVal.prepend(this.lamp);
    const sparkVal = cells.querySelector('.oct-hud-spark .oct-hud-value') as HTMLElement;
    this.sparkDot = document.createElement('span');
    this.sparkDot.className = 'oct-phase-dot';
    sparkVal.prepend(this.sparkDot);

    // ---- transport ----
    const tr = document.createElement('div');
    tr.className = 'oct-hud-transport';
    const btn = (label: string, title: string, fn: () => void, cls = ''): HTMLButtonElement => {
      const b = document.createElement('button');
      b.type = 'button';
      b.className = `oct-tbtn ${cls}`;
      b.innerHTML = label;
      b.title = title;
      b.addEventListener('click', fn);
      tr.appendChild(b);
      return b;
    };
    btn(`−${STEP_LARGE_DEG}°`, `Step back ${STEP_LARGE_DEG}° (Shift+←)`, () => cb.onStep(-STEP_LARGE_DEG));
    btn(`−${STEP_SMALL_DEG}°`, `Step back ${STEP_SMALL_DEG}° (←)`, () => cb.onStep(-STEP_SMALL_DEG));
    this.playBtn = btn('', 'Play / pause (Space)', cb.onTogglePause, 'oct-tbtn-play');
    btn(`+${STEP_SMALL_DEG}°`, `Step forward ${STEP_SMALL_DEG}° (→)`, () => cb.onStep(STEP_SMALL_DEG));
    btn(`+${STEP_LARGE_DEG}°`, `Step forward ${STEP_LARGE_DEG}° (Shift+→)`, () => cb.onStep(STEP_LARGE_DEG));
    const ts = document.createElement('div');
    ts.className = 'oct-hud-ts';
    const slower = document.createElement('button');
    slower.type = 'button';
    slower.className = 'oct-tbtn oct-tbtn-sm';
    slower.textContent = '−';
    slower.title = 'Slower ( [ )';
    slower.addEventListener('click', cb.onSlower);
    const tsVal = document.createElement('span');
    tsVal.className = 'oct-hud-ts-val';
    tsVal.dataset.k = 'ts';
    tsVal.title = 'Playback speed: simulated time per wall-clock time';
    const faster = document.createElement('button');
    faster.type = 'button';
    faster.className = 'oct-tbtn oct-tbtn-sm';
    faster.textContent = '+';
    faster.title = 'Faster ( ] )';
    faster.addEventListener('click', cb.onFaster);
    ts.append(slower, tsVal, faster);
    tr.appendChild(ts);
    this.el.appendChild(tr);

    this.el.querySelectorAll<HTMLElement>('[data-k]').forEach((e) => {
      this.slots[e.dataset.k!] = new Slot(e);
    });
    this.setPlayback(1, false);
    this.setFocus(this.focus);
  }

  /** Focus cylinder (0-based): the dial and readouts show it, its strip chip is outlined. */
  setFocus(index: number): void {
    this.focus = index;
    for (const c of this.strip) {
      const on = Number(c.el.dataset.cyl) === index;
      c.el.classList.toggle('is-focus', on);
      c.el.setAttribute('aria-pressed', String(on));
    }
    this.mapMean.reset();
    this.lastSparkPhase = '';
  }

  /** Latest engine-level cycle results (brake torque / power), or null. */
  setEngineCycle(e: EngineCycleSummary | null): void {
    this.engineCycle = e;
  }

  setSparkAdvance(advanceDeg: number): void {
    const [x0, y0] = cycleDialPoint(-advanceDeg, C, C, R_STROKE + 5);
    const [x1, y1] = cycleDialPoint(-advanceDeg, C, C, R_STROKE - 3);
    this.sparkTick.setAttribute('x1', x0.toFixed(2));
    this.sparkTick.setAttribute('y1', y0.toFixed(2));
    this.sparkTick.setAttribute('x2', x1.toFixed(2));
    this.sparkTick.setAttribute('y2', y1.toFixed(2));
  }

  setPlayback(timeScale: number, paused: boolean): void {
    this.slots.ts?.set(paused ? `${formatTimeScale(timeScale).replace(' (real time)', '')} · paused` : formatTimeScale(timeScale).replace(' (real time)', ' real time'));
    this.playBtn.innerHTML = paused
      ? '<svg viewBox="0 0 12 12" width="12" height="12" aria-hidden="true"><path d="M3 1.5 L10.5 6 L3 10.5 Z" fill="currentColor"/></svg>'
      : '<svg viewBox="0 0 12 12" width="12" height="12" aria-hidden="true"><rect x="2.5" y="1.5" width="2.6" height="9" fill="currentColor"/><rect x="6.9" y="1.5" width="2.6" height="9" fill="currentColor"/></svg>';
    this.playBtn.setAttribute('aria-label', paused ? 'Play' : 'Pause');
    this.el.classList.toggle('is-paused', paused);
  }

  setLastCycle(c: CycleSummary | null): void {
    this.slots.last?.set(c ? `${fmt(paToBar(c.imepNet), 2)} bar · η ${fmt(fractionToPct(c.indicatedEfficiency), 1)} %` : '—');
  }

  /**
   * @param s      the focus cylinder's view of the playback snapshot (cylinder-view.ts), or the snapshot
   * @param engine the playback snapshot itself (engine-level fields, every cylinder); default `s`
   */
  update(s: EngineSnapshot, engine: EngineSnapshot = s): void {
    const sl = this.slots;
    const a = cycleDialAngleDeg(s.thetaDeg);
    if (Math.abs(a - this.lastNeedle) > 0.05) {
      this.lastNeedle = a;
      this.needle.setAttribute('transform', `rotate(${a.toFixed(2)} ${C} ${C})`);
    }
    sl.theta.set(fmt(s.thetaDeg, 1));
    const st = strokeOf(s.thetaDeg);
    if (st !== this.lastStroke) {
      this.lastStroke = st;
      sl.stroke.set(STROKE_LABEL[st]);
      this.strokeEl.style.setProperty('--stroke', STROKE_COLOR[st]);
    }
    sl.rel.set(describeCrankAngle(s.thetaDeg));
    sl.rpm.set(fmt(s.rpm, 0));
    const map = this.mapMean.update(s.intakeManifoldPressure, s.t, s.intakeLift > 0);
    sl.map.set(Number.isFinite(map) ? fmt(paToBar(map), 2) : '—');
    sl.p.set(fmt(paToBar(s.pressure), 2));
    sl.T.set(`${s.temperatureUnburned > 0 ? fmt(s.temperatureUnburned, 0) : '—'} / ${s.temperatureBurned > 0 ? fmt(s.temperatureBurned, 0) : '—'}`);
    const xb = fractionToPct(s.massFractionBurned);
    sl.xb.set(fmt(xb, 1));
    const xbq = Math.round(xb * 2) / 2;
    if (xbq !== this.lastXb) {
      this.lastXb = xbq;
      this.xbBar.style.transform = `scaleX(${Math.min(1, Math.max(0, xb / 100)).toFixed(3)})`;
    }
    sl.flame.set(s.flame.stage);
    const sp = s.spark.phase;
    const count = s.spark.breakdownCount;
    const sparkKey = count === undefined ? sp : `${sp}|${count}`;
    if (sparkKey !== this.lastSparkPhase) {
      this.lastSparkPhase = sparkKey;
      // Trembler shower: phase plus the number of sparks so far in this timer contact.
      sl.spark.set(count !== undefined && count > 0 ? `${sp} · ${count}×` : sp);
      this.sparkDot.dataset.phase = sp;
    }
    const lw = s.knock.integral;
    const lampState = s.knock.autoignited ? 'knock' : lw >= 0.7 ? 'warn' : lw > 0 ? 'low' : 'idle';
    if (lampState !== this.lastLamp) {
      this.lastLamp = lampState;
      this.lamp.dataset.state = lampState;
      this.lamp.style.color = lampState === 'knock' ? STATUS.critical : lampState === 'warn' ? STATUS.warning : '';
    }
    sl.knock.set(s.knock.autoignited ? 'KNOCK' : `∫ ${fmt(Math.min(lw, 9.99) * 100, 0)} %`);
    sl.cycle.set(`#${engine.cycle} · ${fmtSimTime(s.t)}`);

    if (sl.brake) {
      const ec = this.engineCycle;
      if (ec) sl.brake.set(`${fmt(ec.brakeTorque, 1)} N·m · ${fmt(wattsToHp(ec.brakePower), 1)} hp`);
      else if (engine.loadTorque !== undefined) {
        const p = engine.loadTorque * ((2 * Math.PI * engine.rpm) / 60);
        sl.brake.set(`${fmt(engine.loadTorque, 1)} N·m · ${fmt(wattsToHp(p), 1)} hp`);
      } else sl.brake.set('—');
    }
    if (sl.road) {
      const v = engine.vehicleSpeed;
      sl.road.set(v === undefined ? '—' : `${fmt(mpsToMph(v), 1)} mph · ${fmt(mpsToKmh(v), 0)} km/h`);
    }
    if (this.strip.length) this.updateStrip(engine);
  }

  private updateStrip(engine: EngineSnapshot): void {
    const entries = firingStrip(this.spec, engine, this.stripEntries);
    for (let k = 0; k < entries.length && k < this.strip.length; k++) {
      const e = entries[k];
      const chip = this.strip[k];
      const key = `${e.stroke}|${e.sparking ? 1 : 0}|${e.timerClosed ? 1 : 0}`;
      if (key === chip.key) continue;
      chip.key = key;
      chip.el.style.setProperty('--stroke', STROKE_COLOR[e.stroke]);
      chip.el.dataset.stroke = e.stroke;
      chip.el.classList.toggle('is-spark', e.sparking);
      chip.el.classList.toggle('is-timer', e.timerClosed);
      chip.stroke.textContent = STROKE_INITIAL[e.stroke] ?? '';
      chip.el.setAttribute('aria-label', `Cylinder ${e.number}: ${STROKE_LABEL[e.stroke]}${e.sparking ? ', sparking' : ''}`);
    }
  }

  dispose(): void {
    this.el.remove();
  }
}

/** Minimal snapshot for laying out the strip before data arrives (engine angle −360°). */
function makeIdleSnapshot(): EngineSnapshot {
  return { thetaDeg: -360, spark: { phase: 'off' } } as EngineSnapshot;
}
