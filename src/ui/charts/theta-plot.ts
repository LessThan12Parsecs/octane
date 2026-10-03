/**
 * A crank-angle (θ) plot: the current cycle as uPlot series (up to the playback
 * cursor), previous cycles as faded "ghost" traces drawn underneath, plus
 * playhead / spark / knock markers on top. All θ plots share one ZoomGroup.
 */
import { emptyExtent, extendExtent, niceRange, type Extent, type RangeOptions } from '../chart-math';
import { lowerBound, upperBound } from '../ring-buffer';
import { CH, fillSeries, type TraceSelection } from '../trace-store';
import { fmt } from '../units';
import { axis, drawDot, drawVMarker, strokeSeries, tickValues, uPlot, withPlotClip, zoomCursor, type ZoomGroup } from './plot-kit';
import { GHOST_KNOCK_RGB, INK, STATUS, SERIES, withAlpha } from './theme';

export interface ThetaSeriesDef {
  label: string;
  /** Trace channel (CH.*). */
  ch: number;
  color: string;
  unit: string;
  decimals: number;
  dash?: number[];
  width?: number;
  /** Draw this series for previous cycles too (default true). */
  ghost?: boolean;
  /** Light area fill down to zero. */
  fill?: boolean;
}

export interface ThetaPlotDef {
  /** Axis label, e.g. "p [bar]". */
  yLabel: string;
  height: number;
  series: ThetaSeriesDef[];
  range: RangeOptions;
  /** Number of previous cycles drawn as ghosts. */
  ghosts: number;
  /** Show x tick labels (only on the bottom plot of a stack). */
  xLabels: boolean;
  /** Draw the spark-timing marker. */
  sparkMarker?: boolean;
  /** Draw the knock-onset marker and colour knocking ghost cycles. */
  knockMarker?: boolean;
  /** Draw a y = 0 reference line. */
  zeroLine?: boolean;
}

export interface ThetaFrame {
  sel: TraceSelection;
  /** Crank angle at the playback cursor, deg. */
  playheadTheta: number;
  /** Spark timing (−advance), deg. */
  sparkDeg: number;
}

const X_AXIS_LABEL = 'θ [°CA]  (0 = firing TDC)';
const EMPTY = new Float64Array(0);

export class ThetaPlot {
  readonly u: uPlot;
  private readonly ys: (number | null)[][];
  private readonly data: uPlot.AlignedData;
  private frame: ThetaFrame | null = null;
  private readonly ext: Extent = [0, 1];
  private readonly range: Extent = [0, 1];
  private hovering = false;
  private count = 0;

  constructor(
    host: HTMLElement,
    private readonly def: ThetaPlotDef,
    private readonly zoom: ZoomGroup,
    width: number,
    syncKey = 'theta',
  ) {
    this.ys = def.series.map(() => []);
    this.data = [EMPTY, ...this.ys] as uPlot.AlignedData;

    const series: uPlot.Series[] = [
      {
        label: 'θ',
        value: (_u, v) => (v == null ? '—' : `${fmt(v, 1)}°`),
      },
      ...def.series.map(
        (s): uPlot.Series => ({
          label: s.label,
          stroke: s.color,
          width: s.width ?? 1.5,
          dash: s.dash,
          fill: s.fill ? withAlpha(s.color, 0.14) : undefined,
          spanGaps: false,
          points: { show: false },
          value: (_u, v) => (v == null ? '—' : `${fmt(v, s.decimals)} ${s.unit}`),
        }),
      ),
    ];

    const opts: uPlot.Options = {
      width,
      height: def.height,
      series,
      legend: { show: true, live: true },
      cursor: zoomCursor(syncKey),
      padding: [6, 12, def.xLabels ? 0 : 4, 0],
      scales: {
        // range() resolves null/pending x scales (e.g. setData right after construction) to the group window.
        x: { time: false, auto: false, range: () => [zoom.min, zoom.max] },
        y: { range: () => [this.range[0], this.range[1]] },
      },
      axes: [
        axis({
          show: true,
          size: def.xLabels ? 34 : 8,
          label: def.xLabels ? X_AXIS_LABEL : undefined,
          values: def.xLabels ? tickValues() : () => [],
          ticks: { show: def.xLabels, stroke: INK.tick, width: 1, size: 3 },
        }),
        axis({ label: def.yLabel, values: tickValues(), size: 44 }),
      ],
      hooks: {
        ...zoom.hooks(),
        drawAxes: [(u) => this.drawGhosts(u)],
        draw: [(u) => this.drawMarkers(u)],
      },
    };
    this.u = new uPlot(opts, this.data, host);
    zoom.add(this.u);
    this.u.over.addEventListener('mouseenter', this.onEnter);
    this.u.over.addEventListener('mouseleave', this.onLeave);
  }

  private readonly onEnter = (): void => {
    this.hovering = true;
  };

  private readonly onLeave = (): void => {
    this.hovering = false;
    this.syncLegend();
  };

  update(frame: ThetaFrame): void {
    this.frame = frame;
    const { sel } = frame;
    const cur = sel.current;
    const n = cur ? sel.currentCount : 0;
    this.count = n;
    const xmin = this.zoom.min;
    const xmax = this.zoom.max;

    const e = emptyExtent(this.ext);
    const xs = cur ? cur.cols[CH.theta] : null;
    for (let k = 0; k < this.def.series.length; k++) {
      const s = this.def.series[k];
      if (cur && xs) {
        fillSeries(cur.cols[s.ch], n, this.ys[k]);
        extendExtent(xs, cur.cols[s.ch], n, xmin, xmax, e);
      } else {
        this.ys[k].length = 0;
      }
      if (s.ghost !== false) {
        for (let g = 0; g < sel.ghosts.length && g < this.def.ghosts; g++) {
          const tr = sel.ghosts[g];
          extendExtent(tr.cols[CH.theta], tr.cols[s.ch], tr.n, xmin, xmax, e);
        }
      }
    }
    niceRange(e, this.def.range, this.range);
    this.data[0] = cur ? cur.col(CH.theta, n) : EMPTY;
    this.u.setData(this.data);
    this.syncLegend();
  }

  /** Legend shows values at the playback cursor unless the mouse is hovering. */
  private syncLegend(): void {
    if (this.hovering || this.count === 0) return;
    this.u.setLegend({ idx: this.count - 1 });
  }

  private drawGhosts(u: uPlot): void {
    const f = this.frame;
    if (!f || this.def.ghosts <= 0 || f.sel.ghosts.length === 0) return;
    const dpr = uPlot.pxRatio;
    const xmin = this.zoom.min;
    const xmax = this.zoom.max;
    withPlotClip(u, (ctx) => {
      ctx.lineJoin = 'round';
      const G = Math.min(this.def.ghosts, f.sel.ghosts.length);
      for (let g = G - 1; g >= 0; g--) {
        const tr = f.sel.ghosts[g];
        const alpha = 0.42 * Math.pow(0.62, g);
        const th = tr.cols[CH.theta];
        const i0 = Math.max(0, lowerBound(th, tr.n, xmin) - 1);
        const i1 = Math.min(tr.n, upperBound(th, tr.n, xmax) + 1);
        for (const s of this.def.series) {
          if (s.ghost === false) continue;
          const knockTint = this.def.knockMarker && tr.knocked;
          ctx.strokeStyle = knockTint ? `rgba(${GHOST_KNOCK_RGB}, ${alpha + 0.1})` : withAlpha(s.color, alpha);
          ctx.lineWidth = 1 * dpr;
          ctx.setLineDash(s.dash ? s.dash.map((d) => d * dpr) : []);
          strokeSeries(u, ctx, th, tr.cols[s.ch], i0, i1);
        }
      }
    });
  }

  private drawMarkers(u: uPlot): void {
    const f = this.frame;
    if (!f) return;
    if (this.def.zeroLine) {
      const sy = u.scales.y;
      if (sy.min != null && sy.max != null && sy.min < 0 && sy.max > 0) {
        const ctx = u.ctx;
        const py = Math.round(u.valToPos(0, 'y', true)) + 0.5;
        ctx.save();
        ctx.strokeStyle = INK.axis;
        ctx.lineWidth = uPlot.pxRatio;
        ctx.beginPath();
        ctx.moveTo(u.bbox.left, py);
        ctx.lineTo(u.bbox.left + u.bbox.width, py);
        ctx.stroke();
        ctx.restore();
      }
    }
    if (this.def.sparkMarker) drawVMarker(u, f.sparkDeg, withAlpha(SERIES.yellow, 0.85), [3, 3], 'spark');
    const cur = f.sel.current;
    if (this.def.knockMarker && cur && cur.knocked && Number.isFinite(cur.knockDeg) && cur.knockDeg <= f.playheadTheta) {
      drawVMarker(u, cur.knockDeg, STATUS.serious, [2, 2], 'knock', 1, 1);
    }
    drawVMarker(u, f.playheadTheta, 'rgba(230, 233, 238, 0.55)', null);
    if (cur && this.count > 0) {
      const i = this.count - 1;
      const x = cur.cols[CH.theta][i];
      for (const s of this.def.series) drawDot(u, x, cur.cols[s.ch][i], s.color, 3);
    }
  }

  resize(width: number): void {
    if (Math.abs(this.u.width - width) < 1) return;
    this.u.setSize({ width, height: this.def.height });
  }

  dispose(): void {
    this.zoom.remove(this.u);
    this.u.over.removeEventListener('mouseenter', this.onEnter);
    this.u.over.removeEventListener('mouseleave', this.onLeave);
    this.u.destroy();
  }
}
