/**
 * Ignition-event time traces (µs–ms): secondary voltage with the required breakdown
 * voltage, secondary (gap) current and primary current, stacked with a shared,
 * zoomable time axis (t − t_breakdown). Discharge phases are shaded on the top plot.
 */
import { emptyExtent, extendExtent, niceRange, type Extent, type RangeOptions } from '../chart-math';
import { SPARK_PHASE_CODE, type SparkEvent } from '../spark-capture';
import { fillSeries } from '../trace-store';
import { fmt, fmtMicros, fmtMicrosTick } from '../units';
import { axis, drawVMarker, uPlot, withPlotClip, zoomCursor, ZoomGroup } from './plot-kit';
import { INK, SERIES, withAlpha } from './theme';

interface SubDef {
  yLabel: string;
  height: number;
  xLabels: boolean;
  series: { label: string; unit: string; decimals: number; color: string; dash?: number[]; key: 'vSec' | 'vBd' | 'iSec' | 'iPrim' }[];
  range: RangeOptions;
  phases?: boolean;
}

const PHASE_STYLE: Record<number, { color: string; label: string }> = {
  [SPARK_PHASE_CODE.charging]: { color: withAlpha(SERIES.yellow, 0.1), label: 'dwell' },
  [SPARK_PHASE_CODE.breakdown]: { color: 'rgba(230, 233, 238, 0.28)', label: 'breakdown' },
  [SPARK_PHASE_CODE.arc]: { color: withAlpha(SERIES.blue, 0.16), label: 'arc' },
  [SPARK_PHASE_CODE.glow]: { color: withAlpha(SERIES.violet, 0.14), label: 'glow' },
};

export const SPARK_DEFAULT_WINDOW_US: [number, number] = [-200, 2500];

const xTicks: uPlot.Axis.Values = (_u, splits) => splits.map((v) => (v == null ? '' : fmtMicrosTick(v)));

class SparkSubPlot {
  readonly u: uPlot;
  private readonly ys: (number | null)[][];
  private readonly data: uPlot.AlignedData;
  private readonly ext: Extent = [0, 1];
  private readonly range: Extent = [0, 1];
  ev: SparkEvent | null = null;
  n = 0;
  private hovering = false;

  constructor(
    host: HTMLElement,
    readonly def: SubDef,
    private readonly zoom: ZoomGroup,
    private readonly xs: number[],
    width: number,
  ) {
    this.ys = def.series.map(() => []);
    this.data = [xs, ...this.ys] as uPlot.AlignedData;
    const opts: uPlot.Options = {
      width,
      height: def.height,
      series: [
        { label: 't', value: (_u, v) => (v == null ? '—' : fmtMicros(v)) },
        ...def.series.map(
          (s): uPlot.Series => ({
            label: s.label,
            stroke: s.color,
            dash: s.dash,
            width: 1.5,
            spanGaps: false,
            points: { show: false },
            value: (_u, v) => (v == null ? '—' : `${fmt(v, s.decimals)} ${s.unit}`),
          }),
        ),
      ],
      legend: { show: true, live: true },
      cursor: zoomCursor('spark'),
      padding: [6, 22, def.xLabels ? 0 : 4, 0],
      scales: { x: { time: false, auto: false, range: () => [zoom.min, zoom.max] }, y: { range: () => [this.range[0], this.range[1]] } },
      axes: [
        axis({
          size: def.xLabels ? 34 : 8,
          label: def.xLabels ? 't − t_breakdown' : undefined,
          values: def.xLabels ? xTicks : () => [],
          ticks: { show: def.xLabels, stroke: INK.tick, width: 1, size: 3 },
        }),
        axis({ label: def.yLabel, size: 44 }),
      ],
      hooks: {
        ...zoom.hooks(),
        drawAxes: def.phases ? [(u) => this.drawPhases(u)] : [],
        draw: [(u) => drawVMarker(u, 0, 'rgba(230, 233, 238, 0.4)', [2, 3])],
      },
    };
    this.u = new uPlot(opts, this.data, host);
    zoom.add(this.u);
    this.u.over.addEventListener('mouseenter', () => (this.hovering = true));
    this.u.over.addEventListener('mouseleave', () => {
      this.hovering = false;
      this.syncLegend();
    });
  }

  update(ev: SparkEvent | null, n: number): void {
    this.ev = ev;
    this.n = n;
    const e = emptyExtent(this.ext);
    for (let k = 0; k < this.def.series.length; k++) {
      const s = this.def.series[k];
      if (ev) {
        const col = ev[s.key];
        fillSeries(col.data, n, this.ys[k]);
        extendExtent(this.xs, col.data, n, this.zoom.min, this.zoom.max, e);
      } else this.ys[k].length = 0;
    }
    niceRange(e, this.def.range, this.range);
    this.u.setData(this.data);
    this.syncLegend();
  }

  private syncLegend(): void {
    if (this.hovering || this.n === 0) return;
    // Show the latest sample while the event is live, else the peak-voltage region start.
    this.u.setLegend({ idx: this.n - 1 });
  }

  private drawPhases(u: uPlot): void {
    const ev = this.ev;
    if (!ev || this.n < 2) return;
    const dpr = uPlot.pxRatio;
    const ph = ev.phase.data;
    const xs = this.xs;
    withPlotClip(u, (ctx) => {
      const { top, height } = u.bbox;
      ctx.font = `${Math.round(9 * dpr)}px ui-monospace, Menlo, monospace`;
      ctx.textBaseline = 'bottom';
      let i0 = 0;
      for (let i = 1; i <= this.n; i++) {
        if (i < this.n && ph[i] === ph[i0]) continue;
        const style = PHASE_STYLE[ph[i0]];
        if (style) {
          const x0 = u.valToPos(xs[i0], 'x', true);
          const x1 = u.valToPos(xs[Math.min(i, this.n - 1)], 'x', true);
          const w = Math.max(1 * dpr, x1 - x0);
          ctx.fillStyle = style.color;
          ctx.fillRect(x0, top, w, height);
          if (w > 34 * dpr) {
            ctx.fillStyle = INK.secondary;
            ctx.fillText(style.label, x0 + 3 * dpr, top + 12 * dpr);
          }
        }
        i0 = i;
      }
    });
  }

  resize(width: number): void {
    if (Math.abs(this.u.width - width) < 1) return;
    this.u.setSize({ width, height: this.def.height });
  }

  dispose(): void {
    this.zoom.remove(this.u);
    this.u.destroy();
  }
}

/**
 * Default window of a trembler-coil shower, µs around the first breakdown: the coil needs ≈ 2–3.5 ms from
 * timer make to the first spark and then buzzes at ≈ 200 Hz for the rest of the contact [ECCT; Cool, via
 * engines/model-t.ts], 87° of crank ≈ 14.5 ms at 1000 rpm.
 */
export const SPARK_SHOWER_WINDOW_US: [number, number] = [-4000, 12000];

export class SparkPlots {
  readonly zoom: ZoomGroup;
  private readonly xs: number[] = [];
  private readonly plots: SparkSubPlot[];
  private lastEv: SparkEvent | null = null;

  constructor(host: HTMLElement, width: number, window: [number, number] = SPARK_DEFAULT_WINDOW_US) {
    this.zoom = new ZoomGroup(window[0], window[1]);
    const defs: SubDef[] = [
      {
        yLabel: 'V [kV]',
        height: 120,
        xLabels: false,
        phases: true,
        range: { includeMin: 0, minSpan: 1 },
        series: [
          { key: 'vSec', label: 'V secondary', unit: 'kV', decimals: 2, color: SERIES.blue },
          { key: 'vBd', label: 'V breakdown req.', unit: 'kV', decimals: 2, color: SERIES.orange, dash: [4, 3] },
        ],
      },
      {
        yLabel: 'I₂ [mA]',
        height: 80,
        xLabels: false,
        range: { includeMin: 0, minSpan: 1 },
        series: [{ key: 'iSec', label: 'I secondary', unit: 'mA', decimals: 1, color: SERIES.aqua }],
      },
      {
        yLabel: 'I₁ [A]',
        height: 96,
        xLabels: true,
        range: { includeMin: 0, minSpan: 0.5 },
        series: [{ key: 'iPrim', label: 'I primary', unit: 'A', decimals: 2, color: SERIES.yellow }],
      },
    ];
    this.plots = defs.map((d) => {
      const el = document.createElement('div');
      el.className = 'oct-plot';
      host.appendChild(el);
      return new SparkSubPlot(el, d, this.zoom, this.xs, width);
    });
  }

  update(ev: SparkEvent | null): void {
    const n = ev ? ev.length : 0;
    this.xs.length = n;
    if (ev) {
      const t = ev.t.data;
      const t0 = ev.tRef;
      for (let i = 0; i < n; i++) this.xs[i] = (t[i] - t0) * 1e6;
    }
    this.lastEv = ev;
    for (const p of this.plots) p.update(ev, n);
  }

  /** Zoom to show the whole captured event. */
  fit(): void {
    const n = this.xs.length;
    if (n < 2) return this.zoom.reset();
    this.zoom.set(this.xs[0], this.xs[n - 1]);
  }

  get event(): SparkEvent | null {
    return this.lastEv;
  }

  resize(width: number): void {
    for (const p of this.plots) p.resize(width);
  }

  dispose(): void {
    for (const p of this.plots) p.dispose();
  }
}
