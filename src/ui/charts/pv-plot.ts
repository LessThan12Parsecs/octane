/**
 * log p – log V indicator diagram. The loop is not a function of V, so uPlot only
 * provides the log axes/grid; the traces are drawn in hooks: previous cycles faded,
 * the current cycle coloured by stroke and traced live up to the playback cursor.
 */
import { emptyExtent, extendExtent, logTicks125, niceLogRange, type Extent } from '../chart-math';
import { strokeOf, type StrokeName } from '../engine-cycle';
import { CH, type TraceSelection } from '../trace-store';
import { axis, drawDot, strokeSeries, uPlot, withPlotClip } from './plot-kit';
import { STROKE_COLOR, withAlpha } from './theme';

export interface PvFrame {
  sel: TraceSelection;
}

const fmtTick = (v: number): string => {
  if (v >= 100) return v.toFixed(0);
  if (v >= 1) return Number(v.toPrecision(3)).toString();
  return Number(v.toPrecision(2)).toString();
};

export class PvPlot {
  readonly u: uPlot;
  private frame: PvFrame | null = null;
  private readonly ext: Extent = [Infinity, -Infinity];
  private readonly yRange: Extent = [0.5, 100];
  private xMin = 50;
  private xMax = 1000;
  private readonly data: uPlot.AlignedData;
  private readonly ghosts = 4;

  constructor(host: HTMLElement, width: number, private readonly height: number) {
    this.data = [
      [this.xMin, this.xMax],
      [null, null],
    ];
    const opts: uPlot.Options = {
      width,
      height,
      series: [{ label: 'V' }, { label: 'p', stroke: 'transparent', points: { show: false } }],
      legend: { show: false },
      cursor: { show: false, drag: { x: false, y: false } },
      padding: [6, 10, 0, 0],
      scales: {
        x: { time: false, distr: 3, log: 10, auto: false, range: () => [this.xMin, this.xMax] },
        y: { distr: 3, log: 10, range: () => [this.yRange[0], this.yRange[1]] },
      },
      axes: [
        axis({
          label: 'V [cm³]  (log)',
          size: 34,
          splits: (_u, _i, min, max) => logTicks125(min, max, []),
          values: (_u, splits) => splits.map((v) => (v == null ? '' : fmtTick(v))),
        }),
        axis({
          label: 'p [bar]  (log)',
          size: 44,
          splits: (_u, _i, min, max) => logTicks125(min, max, []),
          values: (_u, splits) => splits.map((v) => (v == null ? '' : fmtTick(v))),
        }),
      ],
      hooks: {
        drawAxes: [(u) => this.drawTraces(u)],
      },
    };
    this.u = new uPlot(opts, this.data, host);
    this.u.setScale('x', { min: this.xMin, max: this.xMax });
  }

  /** Set the V axis from the cylinder's clearance and maximum volumes (cm³). */
  setVolumeRange(vcCm3: number, vmaxCm3: number): void {
    const lo = Number((vcCm3 * 0.8).toPrecision(2));
    const hi = Number((vmaxCm3 * 1.2).toPrecision(2));
    if (lo === this.xMin && hi === this.xMax) return;
    this.xMin = lo;
    this.xMax = hi;
    this.data[0] = [lo, hi];
    this.u.setScale('x', { min: lo, max: hi });
  }

  update(frame: PvFrame): void {
    this.frame = frame;
    const { sel } = frame;
    const e = emptyExtent(this.ext);
    if (sel.current) extendExtent(null, sel.current.cols[CH.p], sel.currentCount, 0, 0, e);
    for (let g = 0; g < sel.ghosts.length && g < this.ghosts; g++) {
      const tr = sel.ghosts[g];
      extendExtent(null, tr.cols[CH.p], tr.n, 0, 0, e);
    }
    if (e[0] <= 0) e[0] = 0.1;
    niceLogRange(e, 0.05, this.yRange);
    this.u.setData(this.data);
  }

  private drawTraces(u: uPlot): void {
    const f = this.frame;
    if (!f) return;
    const dpr = uPlot.pxRatio;
    withPlotClip(u, (ctx) => {
      ctx.lineJoin = 'round';
      ctx.lineCap = 'round';
      const G = Math.min(this.ghosts, f.sel.ghosts.length);
      for (let g = G - 1; g >= 0; g--) {
        const tr = f.sel.ghosts[g];
        ctx.strokeStyle = withAlpha('#a3acb9', 0.3 * Math.pow(0.65, g));
        ctx.lineWidth = dpr;
        strokeSeries(u, ctx, tr.cols[CH.V], tr.cols[CH.p], 0, tr.n);
      }
      const cur = f.sel.current;
      const n = f.sel.currentCount;
      if (!cur || n < 1) return;
      const V = cur.cols[CH.V];
      const P = cur.cols[CH.p];
      const TH = cur.cols[CH.theta];
      ctx.lineWidth = 1.75 * dpr;
      let seg0 = 0;
      let stroke: StrokeName = strokeOf(TH[0]);
      for (let i = 1; i <= n; i++) {
        const s = i < n ? strokeOf(TH[i]) : null;
        if (s !== stroke) {
          ctx.strokeStyle = STROKE_COLOR[stroke];
          // Overlap one sample so consecutive stroke segments join.
          strokeSeries(u, ctx, V, P, seg0, Math.min(n, i + 1));
          seg0 = i;
          if (s) stroke = s;
        }
      }
      drawDot(u, V[n - 1], P[n - 1], STROKE_COLOR[strokeOf(TH[n - 1])], 4);
    });
  }

  resize(width: number): void {
    if (Math.abs(this.u.width - width) < 1) return;
    this.u.setSize({ width, height: this.height });
  }

  dispose(): void {
    this.u.destroy();
  }
}
