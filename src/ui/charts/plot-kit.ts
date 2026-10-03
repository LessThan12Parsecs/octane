/**
 * Shared uPlot plumbing: dark axis styling, x-zoom groups, canvas overlay helpers
 * (ghost traces, markers) and sizing.
 */
import uPlot from 'uplot';
import 'uplot/dist/uPlot.min.css';
import { FONT_MONO, FONT_UI, INK } from './theme';

export { uPlot };

export function axis(opts: Partial<uPlot.Axis> & { label?: string } = {}): uPlot.Axis {
  return {
    stroke: INK.secondary,
    font: FONT_MONO,
    labelFont: FONT_UI,
    labelSize: opts.label ? 16 : 0,
    labelGap: 0,
    size: 38,
    gap: 3,
    grid: { stroke: INK.grid, width: 1 },
    ticks: { stroke: INK.tick, width: 1, size: 3 },
    border: { show: true, stroke: INK.axis, width: 1 },
    ...opts,
  };
}

/** Tick formatter with a fixed number of decimals chosen from the tick increment. */
export function tickValues(suffix = ''): uPlot.Axis.Values {
  return (_u, splits, _ax, _space, incr) => {
    const d = incr >= 1 || incr === 0 ? 0 : Math.min(6, Math.ceil(-Math.log10(incr) - 1e-9));
    return splits.map((v) => (v == null ? '' : `${(Math.abs(v) < 1e-12 ? 0 : v).toFixed(d).replace('-', '−')}${suffix}`));
  };
}

/**
 * Keeps the x scale of several plots in lock-step (drag-zoom on any of them,
 * double-click or `reset()` to return to the default range).
 */
export class ZoomGroup {
  private readonly plots = new Set<uPlot>();
  private applying = false;
  min: number;
  max: number;
  private readonly listeners = new Set<() => void>();

  constructor(
    public defaultMin: number,
    public defaultMax: number,
  ) {
    this.min = defaultMin;
    this.max = defaultMax;
  }

  get zoomed(): boolean {
    return this.min !== this.defaultMin || this.max !== this.defaultMax;
  }

  /** Hooks for plot options: propagates user zooms to the group. */
  hooks(): uPlot.Hooks.Arrays {
    return {
      setScale: [
        (u, key) => {
          if (key !== 'x' || this.applying) return;
          const { min, max } = u.scales.x;
          if (min == null || max == null) return;
          if (min === this.min && max === this.max) return;
          this.set(min, max, u);
        },
      ],
    };
  }

  add(u: uPlot): void {
    this.plots.add(u);
    u.over.addEventListener('dblclick', this.onDblClick);
    this.applyTo(u);
  }

  remove(u: uPlot): void {
    this.plots.delete(u);
    u.over.removeEventListener('dblclick', this.onDblClick);
  }

  onChange(cb: () => void): () => void {
    this.listeners.add(cb);
    return () => this.listeners.delete(cb);
  }

  set(min: number, max: number, source?: uPlot): void {
    if (!(max > min)) return;
    this.min = min;
    this.max = max;
    this.applying = true;
    try {
      for (const u of this.plots) if (u !== source) this.applyTo(u);
    } finally {
      this.applying = false;
    }
    for (const cb of this.listeners) cb();
  }

  setDefault(min: number, max: number, resetView = true): void {
    const wasZoomed = this.zoomed;
    this.defaultMin = min;
    this.defaultMax = max;
    if (resetView || !wasZoomed) this.set(min, max);
  }

  reset(): void {
    this.set(this.defaultMin, this.defaultMax);
  }

  private applyTo(u: uPlot): void {
    const was = this.applying;
    this.applying = true;
    try {
      u.setScale('x', { min: this.min, max: this.max });
    } finally {
      this.applying = was;
    }
  }

  private readonly onDblClick = (): void => this.reset();
}

/** Cursor options: x drag-zoom (no default double-click autoscale), optional sync key. */
export function zoomCursor(syncKey?: string): uPlot.Cursor {
  return {
    drag: { x: true, y: false, setScale: true },
    bind: { dblclick: () => null },
    points: { size: 5, width: 1 },
    ...(syncKey ? { sync: { key: syncKey, setSeries: false } } : {}),
  };
}

/** Run `fn` with the context clipped to the plotting area. */
export function withPlotClip(u: uPlot, fn: (ctx: CanvasRenderingContext2D) => void): void {
  const ctx = u.ctx;
  const { left, top, width, height } = u.bbox;
  ctx.save();
  ctx.beginPath();
  ctx.rect(left, top, width, height);
  ctx.clip();
  try {
    fn(ctx);
  } finally {
    ctx.restore();
  }
}

/**
 * Stroke a polyline (xs[i], ys[i]·yScale), i ∈ [i0, i1), breaking at non-finite
 * values. Coordinates go through the plot's scales (canvas pixels).
 */
export function strokeSeries(
  u: uPlot,
  ctx: CanvasRenderingContext2D,
  xs: ArrayLike<number>,
  ys: ArrayLike<number>,
  i0: number,
  i1: number,
  yKey = 'y',
  xKey = 'x',
): void {
  ctx.beginPath();
  let pen = false;
  for (let i = i0; i < i1; i++) {
    const x = xs[i];
    const y = ys[i];
    if (!Number.isFinite(x) || !Number.isFinite(y)) {
      pen = false;
      continue;
    }
    const px = u.valToPos(x, xKey, true);
    const py = u.valToPos(y, yKey, true);
    if (pen) ctx.lineTo(px, py);
    else {
      ctx.moveTo(px, py);
      pen = true;
    }
  }
  ctx.stroke();
}

/** Vertical marker line at x (data units) with an optional small label at the top. */
export function drawVMarker(
  u: uPlot,
  x: number,
  color: string,
  dash: number[] | null,
  label?: string,
  lineWidth = 1,
  labelRow = 0,
): void {
  const sx = u.scales.x;
  if (!Number.isFinite(x) || sx.min == null || sx.max == null || x < sx.min || x > sx.max) return;
  const ctx = u.ctx;
  const px = Math.round(u.valToPos(x, 'x', true)) + 0.5;
  const { top, height } = u.bbox;
  const dpr = uPlot.pxRatio;
  ctx.save();
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth * dpr;
  if (dash) ctx.setLineDash(dash.map((d) => d * dpr));
  ctx.beginPath();
  ctx.moveTo(px, top);
  ctx.lineTo(px, top + height);
  ctx.stroke();
  if (label) {
    ctx.setLineDash([]);
    ctx.font = `${Math.round(9 * dpr)}px ui-monospace, Menlo, monospace`;
    ctx.fillStyle = color;
    ctx.textBaseline = 'top';
    ctx.textAlign = px > u.bbox.left + u.bbox.width - 40 * dpr ? 'right' : 'left';
    ctx.fillText(label, px + (ctx.textAlign === 'right' ? -3 : 3) * dpr, top + (2 + 11 * labelRow) * dpr);
  }
  ctx.restore();
}

/** Filled dot at (x, y) data coordinates. */
export function drawDot(u: uPlot, x: number, y: number, color: string, r = 3.5, yKey = 'y'): void {
  if (!Number.isFinite(x) || !Number.isFinite(y)) return;
  const ctx = u.ctx;
  const dpr = uPlot.pxRatio;
  const px = u.valToPos(x, 'x', true);
  const py = u.valToPos(y, yKey, true);
  ctx.save();
  ctx.fillStyle = color;
  ctx.strokeStyle = '#12161c';
  ctx.lineWidth = 2 * dpr;
  ctx.beginPath();
  ctx.arc(px, py, r * dpr, 0, 2 * Math.PI);
  ctx.stroke();
  ctx.fill();
  ctx.restore();
}

/** Width available for a plot inside `el` (content box). */
export function plotWidth(el: HTMLElement): number {
  return Math.max(160, Math.floor(el.clientWidth));
}
