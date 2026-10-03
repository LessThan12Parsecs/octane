/**
 * The "Traces" tab: live crank-angle charts, the p–V diagram and the spark trace.
 * Only cards that are expanded, scrolled into view and on the active tab redraw,
 * and only when their inputs changed.
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { CycleSummary, EngineSnapshot } from '../../physics/core/snapshot';
import { cylinderVolumes, STROKE_LABEL, STROKES } from '../engine-cycle';
import { sparkPhaseDurations, columnMax, type SparkEvent } from '../spark-capture';
import { CH, polytropicIndex, type TraceSelection } from '../trace-store';
import { fmt, fmtMicros, fractionToPct, jToMJ, kgPerSToGPerS, m3ToCm3, mToMm, paToBar } from '../units';
import { Card } from './card';
import { ZoomGroup } from './plot-kit';
import { PvPlot } from './pv-plot';
import { SparkPlots } from './spark-plots';
import { SERIES, STATUS, STROKE_COLOR, withAlpha } from './theme';
import { ThetaPlot, type ThetaFrame } from './theta-plot';

export interface TraceFrame {
  sel: TraceSelection;
  /** Store version (changes when new samples arrive). */
  version: number;
  playhead: EngineSnapshot;
  sparkDeg: number;
  compressionRatio: number;
  spark: SparkEvent | null;
  sparkVersion: number;
  lastCycle: CycleSummary | null;
}

const THETA_VIEWS: { label: string; title: string; range: [number, number] }[] = [
  { label: '720°', title: 'Whole cycle (−360° … 360°)', range: [-360, 360] },
  { label: '±180°', title: 'Compression and expansion strokes', range: [-180, 180] },
  { label: 'Burn', title: 'Combustion window (−40° … 80°)', range: [-40, 80] },
  { label: 'TDC', title: 'Close-up around firing TDC (−10° … 40°) — knock oscillations', range: [-10, 40] },
];

export class TracePanel {
  readonly el: HTMLElement;
  readonly thetaZoom = new ZoomGroup(-360, 360);
  private readonly cards: Card[] = [];
  private readonly thetaPlots: { card: Card; plot: ThetaPlot }[] = [];
  private readonly pvCard: Card;
  private readonly pv: PvPlot;
  private readonly sparkCard: Card;
  private readonly spark: SparkPlots;
  private readonly cPressure: Card;
  private readonly cBurn: Card;
  private readonly cTemp: Card;
  private readonly cGas: Card;
  private readonly thetaCards: Card[];
  private readonly thetaFrame: ThetaFrame;
  private readonly io: IntersectionObserver | null = null;
  private readonly ro: ResizeObserver | null = null;
  private readonly viewButtons: HTMLButtonElement[] = [];
  private lastVersion = -1;
  private lastPlayheadT = NaN;
  private lastSparkVersion = -1;
  private lastSparkReadout = -1;
  private lastCR = NaN;
  private lastSparkDeg = NaN;
  private zoomDirty = true;
  private width = 380;
  private pvStats = { version: -2, nc: NaN, ne: NaN };

  constructor(
    host: HTMLElement,
    private readonly spec: EngineSpec,
    scrollRoot: HTMLElement,
  ) {
    this.el = document.createElement('div');
    this.el.className = 'oct-traces';
    host.appendChild(this.el);
    const hw = host.clientWidth; // card border: 2 px
    this.width = hw > 0 ? Math.max(200, Math.floor(hw - 2)) : 380;

    // --- θ-window toolbar ---
    const bar = document.createElement('div');
    bar.className = 'oct-toolbar';
    const lbl = document.createElement('span');
    lbl.className = 'oct-toolbar-label';
    lbl.textContent = 'Crank window';
    bar.appendChild(lbl);
    const seg = document.createElement('div');
    seg.className = 'oct-seg';
    for (const v of THETA_VIEWS) {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = v.label;
      b.title = v.title;
      b.addEventListener('click', () => this.thetaZoom.set(v.range[0], v.range[1]));
      seg.appendChild(b);
      this.viewButtons.push(b);
    }
    bar.appendChild(seg);
    const hint = document.createElement('span');
    hint.className = 'oct-toolbar-hint';
    hint.textContent = 'drag to zoom · double-click resets';
    bar.appendChild(hint);
    this.el.appendChild(bar);
    this.thetaZoom.onChange(() => {
      this.zoomDirty = true;
      this.syncViewButtons();
    });
    this.syncViewButtons();

    this.thetaFrame = { sel: { current: null, currentCount: 0, ghosts: [] }, playheadTheta: 0, sparkDeg: 0 };

    // --- cylinder pressure ---
    this.cPressure = this.addCard('pressure', 'Cylinder pressure', 'Cylinder pressure vs crank angle; faded traces are previous cycles (knocking cycles tinted).');
    this.cPressure.setSwatches([
      { label: 'previous cycles', color: withAlpha(SERIES.blue, 0.45) },
      { label: 'knocking cycle', color: STATUS.serious },
      { label: 'spark', color: SERIES.yellow, dashed: true },
    ]);
    this.addTheta(this.cPressure, {
      yLabel: 'p [bar]',
      height: 170,
      ghosts: 5,
      xLabels: true,
      sparkMarker: true,
      knockMarker: true,
      range: { includeMin: 0, minSpan: 2 },
      series: [{ label: 'p', ch: CH.p, color: SERIES.blue, unit: 'bar', decimals: 2 }],
    });

    // --- p–V ---
    this.pvCard = this.addCard('pv', 'log p – log V', 'Indicator diagram on log–log axes; the current cycle is coloured by stroke.');
    this.pvCard.setSwatches([
      ...STROKES.map((s) => ({ label: STROKE_LABEL[s], color: STROKE_COLOR[s] })),
      { label: 'previous cycles', color: 'rgba(163, 172, 185, 0.5)' },
    ]);
    this.pv = new PvPlot(this.pvCard.plotHost(), this.width, 200);
    this.cards.push(this.pvCard);

    // --- combustion ---
    this.cBurn = this.addCard('burn', 'Combustion', 'Burned mass fraction, end-gas Livengood–Wu knock integral (autoignition at 1) and heat-release rate.');
    this.addTheta(this.cBurn, {
      yLabel: 'x_b, knock int. [–]',
      height: 120,
      ghosts: 1,
      xLabels: false,
      sparkMarker: true,
      knockMarker: true,
      range: { includeMin: 0, includeMax: 1, pad: 0.02 },
      series: [
        { label: 'x_b burned', ch: CH.xb, color: SERIES.blue, unit: '', decimals: 3 },
        { label: 'knock integral', ch: CH.lw, color: SERIES.orange, unit: '', decimals: 3, dash: [4, 3] },
      ],
    });
    this.addTheta(this.cBurn, {
      yLabel: 'dQ/dθ [J/°]',
      height: 110,
      ghosts: 1,
      xLabels: true,
      sparkMarker: true,
      range: { includeMin: 0, minSpan: 1 },
      zeroLine: true,
      series: [{ label: 'heat release', ch: CH.hrr, color: SERIES.aqua, unit: 'J/°', decimals: 2, fill: true }],
    });

    // --- temperatures ---
    this.cTemp = this.addCard('temps', 'Gas temperatures', 'Unburned-zone, burned-zone and mass-averaged temperatures.');
    this.addTheta(this.cTemp, {
      yLabel: 'T [K]',
      height: 150,
      ghosts: 1,
      xLabels: true,
      sparkMarker: true,
      range: { includeMin: 250, minSpan: 100 },
      series: [
        { label: 'T unburned', ch: CH.Tu, color: SERIES.blue, unit: 'K', decimals: 0 },
        { label: 'T burned', ch: CH.Tb, color: SERIES.orange, unit: 'K', decimals: 0 },
        { label: 'T mean', ch: CH.Tm, color: SERIES.aqua, unit: 'K', decimals: 0, dash: [4, 3] },
      ],
    });

    // --- spark ---
    this.sparkCard = this.addCard('spark', 'Ignition', 'Coil secondary voltage vs the breakdown voltage the gap requires, gap current and primary current, around the spark event.');
    this.sparkCard.addAction('Fit', 'Show the whole captured event', () => this.spark.fit());
    this.sparkCard.addAction('Reset', 'Default window (−200 µs … 2.5 ms)', () => this.spark.zoom.reset());
    this.spark = new SparkPlots(this.sparkCard.body, this.width);
    this.spark.zoom.onChange(() => (this.sparkCard.stale = true));

    // --- gas exchange ---
    this.cGas = this.addCard('gas', 'Gas exchange', 'Valve lifts and valve mass flows (intake: + into cylinder; exhaust: + out of cylinder).');
    this.addTheta(this.cGas, {
      yLabel: 'lift [mm]',
      height: 100,
      ghosts: 0,
      xLabels: false,
      range: { includeMin: 0, minSpan: 1 },
      series: [
        { label: 'intake lift', ch: CH.liftIn, color: STROKE_COLOR.intake, unit: 'mm', decimals: 2 },
        { label: 'exhaust lift', ch: CH.liftEx, color: STROKE_COLOR.exhaust, unit: 'mm', decimals: 2 },
      ],
    });
    this.addTheta(this.cGas, {
      yLabel: 'ṁ [g/s]',
      height: 120,
      ghosts: 1,
      xLabels: true,
      zeroLine: true,
      range: { minSpan: 1 },
      series: [
        { label: 'intake flow', ch: CH.mIn, color: STROKE_COLOR.intake, unit: 'g/s', decimals: 1 },
        { label: 'exhaust flow', ch: CH.mEx, color: STROKE_COLOR.exhaust, unit: 'g/s', decimals: 1 },
      ],
    });

    this.thetaCards = [this.cPressure, this.cBurn, this.cTemp, this.cGas];

    if (typeof IntersectionObserver !== 'undefined') {
      this.io = new IntersectionObserver(
        (entries) => {
          for (const e of entries) {
            const card = this.cards.find((c) => c.el === e.target);
            if (!card) continue;
            const was = card.inView;
            card.inView = e.isIntersecting;
            if (!was && card.inView) card.stale = true;
          }
        },
        { root: scrollRoot, rootMargin: '80px 0px' },
      );
      for (const c of this.cards) this.io.observe(c.el);
    }
    if (typeof ResizeObserver !== 'undefined') {
      this.ro = new ResizeObserver(() => this.resize());
      this.ro.observe(this.el);
    }
  }

  private addCard(id: string, title: string, hint: string): Card {
    const c = new Card(this.el, id, title, hint);
    this.cards.push(c);
    return c;
  }

  private addTheta(card: Card, def: ConstructorParameters<typeof ThetaPlot>[1]): void {
    const plot = new ThetaPlot(card.plotHost(), def, this.thetaZoom, this.width);
    this.thetaPlots.push({ card, plot });
  }

  private syncViewButtons(): void {
    THETA_VIEWS.forEach((v, i) => {
      const on = v.range[0] === this.thetaZoom.min && v.range[1] === this.thetaZoom.max;
      this.viewButtons[i]?.classList.toggle('is-on', on);
    });
  }

  /** Plot width from the card body; 0 while hidden (display: none). */
  private measureWidth(): number {
    const body = this.cards[0]?.body ?? this.el;
    const w = body.clientWidth || this.el.clientWidth;
    return w > 0 ? Math.max(200, Math.floor(w)) : 0;
  }

  resize(): void {
    const w = this.measureWidth();
    if (w === 0 || w === this.width) return;
    this.width = w;
    for (const { plot } of this.thetaPlots) plot.resize(w);
    this.pv.resize(w);
    this.spark.resize(w);
    for (const c of this.cards) c.stale = true;
  }

  /** Redraw what is visible and changed. `visible` = the tab/panel is shown. */
  update(f: TraceFrame, visible: boolean): void {
    if (!visible) return;
    const dataChanged = f.version !== this.lastVersion || f.playhead.t !== this.lastPlayheadT;
    const markerChanged = f.sparkDeg !== this.lastSparkDeg;
    const zoomChanged = this.zoomDirty;
    this.lastVersion = f.version;
    this.lastPlayheadT = f.playhead.t;
    this.lastSparkDeg = f.sparkDeg;
    this.zoomDirty = false;

    const tf = this.thetaFrame;
    tf.sel = f.sel;
    tf.playheadTheta = f.playhead.thetaDeg;
    tf.sparkDeg = f.sparkDeg;

    const anyChange = dataChanged || markerChanged || zoomChanged;
    for (const c of this.thetaCards) c.redraw = c.active && (anyChange || c.stale);
    for (const { card, plot } of this.thetaPlots) if (card.redraw) plot.update(tf);
    for (const c of this.thetaCards) if (c.redraw) c.stale = false;

    if (this.pvCard.active && (dataChanged || this.pvCard.stale || f.compressionRatio !== this.lastCR)) {
      if (f.compressionRatio !== this.lastCR) {
        const v = cylinderVolumes(this.spec, f.compressionRatio);
        this.pv.setVolumeRange(m3ToCm3(v.vc), m3ToCm3(v.vmax));
        this.lastCR = f.compressionRatio;
      }
      this.pv.update({ sel: f.sel });
      this.pvCard.stale = false;
    }

    const sparkChanged = f.sparkVersion !== this.lastSparkVersion;
    if (this.sparkCard.active && (sparkChanged || this.sparkCard.stale)) {
      this.spark.update(f.spark);
      this.lastSparkVersion = f.sparkVersion;
      this.sparkCard.stale = false;
    }
    if (sparkChanged && f.sparkVersion !== this.lastSparkReadout) {
      this.lastSparkReadout = f.sparkVersion;
      this.updateSparkReadout(f.spark);
    }

    if (dataChanged) this.updateReadouts(f);
  }

  private updateReadouts(f: TraceFrame): void {
    const s = f.playhead;
    const c = f.lastCycle;
    const peak = c ? ` · p_max ${fmt(paToBar(c.peakPressure), 1)} @ ${fmt(c.peakPressureDeg, 1)}°` : '';
    this.cPressure.setReadout(`${fmt(paToBar(s.pressure), 2)} bar${peak}`);

    // Polytropic indices from the most recent complete cycle.
    const ref = f.sel.ghosts[0] ?? null;
    if (ref && this.pvStats.version !== ref.cycle) {
      this.pvStats.version = ref.cycle;
      this.pvStats.nc = polytropicIndex(ref, ref.n, -120, -40);
      this.pvStats.ne = polytropicIndex(ref, ref.n, 60, 120);
    }
    this.pvCard.setReadout(
      `${fmt(m3ToCm3(s.volume), 0)} cm³ · n_comp ${fmt(this.pvStats.nc, 3)} · n_exp ${fmt(this.pvStats.ne, 3)}`,
    );

    const ca50 = c ? ` · last CA50 ${fmt(c.ca50, 1)}°` : '';
    this.cBurn.setReadout(`x_b ${fmt(fractionToPct(s.massFractionBurned), 1)} %${ca50}`);
    this.cTemp.setReadout(
      `T_u ${s.temperatureUnburned > 0 ? fmt(s.temperatureUnburned, 0) : '—'} K · T_b ${s.temperatureBurned > 0 ? fmt(s.temperatureBurned, 0) : '—'} K`,
    );
    this.cGas.setReadout(
      `lift ${fmt(mToMm(s.intakeLift), 2)} / ${fmt(mToMm(s.exhaustLift), 2)} mm · ṁ ${fmt(kgPerSToGPerS(s.intakeMassFlow), 1)} / ${fmt(kgPerSToGPerS(s.exhaustMassFlow), 1)} g/s`,
    );
  }

  private updateSparkReadout(ev: SparkEvent | null): void {
    if (!ev) {
      this.sparkCard.setReadout('waiting for the first spark…');
      return;
    }
    const d = sparkPhaseDurations(ev);
    const vPeak = columnMax(ev.vSec);
    const broke = Number.isFinite(ev.tBreakdown);
    const parts = [
      broke ? `V_peak ${fmt(vPeak, 1)} kV` : `no breakdown (V_peak ${fmt(vPeak, 1)} kV)`,
      `E ${fmt(jToMJ(ev.energy), 1)} mJ`,
    ];
    if (d.arc > 0) parts.push(`arc ${fmtMicros(d.arc * 1e6)}`);
    if (d.glow > 0) parts.push(`glow ${fmtMicros(d.glow * 1e6)}`);
    if (!ev.complete) parts.push('live');
    this.sparkCard.setReadout(parts.join(' · '));
    this.sparkCard.readout.style.color = broke ? '' : STATUS.warning;
  }

  /** Mark all cards for redraw (e.g. when the panel becomes visible again). */
  invalidate(): void {
    for (const c of this.cards) c.stale = true;
    this.resize();
  }

  dispose(): void {
    this.io?.disconnect();
    this.ro?.disconnect();
    for (const { plot } of this.thetaPlots) plot.dispose();
    this.pv.dispose();
    this.spark.dispose();
    this.el.remove();
  }
}
