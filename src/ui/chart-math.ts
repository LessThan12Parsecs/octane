/** Pure numeric helpers for chart ranging (no DOM, no uPlot). */

/** Running [min, max] accumulator. */
export type Extent = [number, number];

export function emptyExtent(out: Extent = [Infinity, -Infinity]): Extent {
  out[0] = Infinity;
  out[1] = -Infinity;
  return out;
}

export function isEmptyExtent(e: Extent): boolean {
  return !(e[0] <= e[1]);
}

/**
 * Widen `acc` by the finite values ys[i] (i < n) whose xs[i] lies in [xmin, xmax].
 * xs may be null when every sample is inside the window.
 */
export function extendExtent(
  xs: ArrayLike<number> | null,
  ys: ArrayLike<number>,
  n: number,
  xmin: number,
  xmax: number,
  acc: Extent,
  scale = 1,
): Extent {
  for (let i = 0; i < n; i++) {
    if (xs) {
      const x = xs[i];
      if (x < xmin || x > xmax) continue;
    }
    const y = ys[i] * scale;
    if (!Number.isFinite(y)) continue;
    if (y < acc[0]) acc[0] = y;
    if (y > acc[1]) acc[1] = y;
  }
  return acc;
}

/** 1-2-2.5-5 × 10ⁿ "nice" number ≥ |v| (sign preserved). */
export function niceCeil(v: number): number {
  if (!Number.isFinite(v) || v === 0) return v;
  if (v < 0) return -niceFloor(-v);
  const mag = Math.pow(10, Math.floor(Math.log10(v)));
  const f = v / mag;
  const steps = [1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10];
  for (const s of steps) if (f <= s + 1e-12) return s * mag;
  return 10 * mag;
}

/** Nice number ≤ v for v > 0 (sign preserved). */
export function niceFloor(v: number): number {
  if (!Number.isFinite(v) || v === 0) return v;
  if (v < 0) return -niceCeil(-v);
  const mag = Math.pow(10, Math.floor(Math.log10(v)));
  const f = v / mag;
  const steps = [10, 8, 6, 5, 4, 3, 2.5, 2, 1.5, 1.2, 1];
  for (const s of steps) if (f >= s - 1e-12) return s * mag;
  return mag;
}

export interface RangeOptions {
  /** Fixed lower / upper bound (overrides data). */
  min?: number;
  max?: number;
  /** Bounds that the data range always includes (e.g. 0 for a bar/pressure axis). */
  includeMin?: number;
  includeMax?: number;
  /** Fractional headroom added above/below the data span. */
  pad?: number;
  /** Minimum total span (prevents a flat line filling the whole plot). */
  minSpan?: number;
}

/**
 * Stable, readable y-range for live data: pads the data extent, snaps outward to
 * nice numbers (so the axis doesn't jitter every frame) and honours fixed bounds.
 */
export function niceRange(e: Extent, o: RangeOptions = {}, out: Extent = [0, 1]): Extent {
  let lo = e[0];
  let hi = e[1];
  if (!(lo <= hi)) {
    lo = o.min ?? o.includeMin ?? 0;
    hi = o.max ?? o.includeMax ?? lo + (o.minSpan ?? 1);
  }
  if (o.includeMin !== undefined) lo = Math.min(lo, o.includeMin);
  if (o.includeMax !== undefined) hi = Math.max(hi, o.includeMax);
  const minSpan = o.minSpan ?? 0;
  if (hi - lo < minSpan) {
    const mid = (hi + lo) / 2;
    lo = mid - minSpan / 2;
    hi = mid + minSpan / 2;
    if (o.includeMin !== undefined && lo < o.includeMin && e[0] >= o.includeMin) {
      lo = o.includeMin;
      hi = lo + minSpan;
    }
  }
  const pad = (o.pad ?? 0.05) * (hi - lo);
  if (o.min === undefined && lo !== o.includeMin) lo -= pad;
  if (o.max === undefined) hi += pad;
  lo = o.min ?? snapDown(lo, hi - lo);
  hi = o.max ?? snapUp(hi, hi - lo);
  if (!(hi > lo)) hi = lo + 1;
  out[0] = lo;
  out[1] = hi;
  return out;
}

/** Round v up to a multiple of a nice step of the given span. */
function snapUp(v: number, span: number): number {
  const step = niceFloor(Math.max(span, 1e-300) / 10);
  return Math.ceil(v / step - 1e-9) * step;
}

function snapDown(v: number, span: number): number {
  const step = niceFloor(Math.max(span, 1e-300) / 10);
  return Math.floor(v / step + 1e-9) * step;
}

/** Log-scale range: decade-ish nice bounds around a positive extent. */
export function niceLogRange(e: Extent, floor: number, out: Extent = [1, 10]): Extent {
  let lo = Math.max(e[0], floor);
  let hi = Math.max(e[1], lo * 1.0001);
  if (!(lo > 0) || !Number.isFinite(lo)) lo = floor;
  if (!Number.isFinite(hi) || hi <= lo) hi = lo * 10;
  // 1-2-5 grid in log space.
  out[0] = niceFloor(lo / 1.15);
  out[1] = niceCeil(hi * 1.15);
  return out;
}

/** Tick values on a log axis between min and max using a 1-2-5 sequence. */
export function logTicks125(min: number, max: number, out: number[] = []): number[] {
  out.length = 0;
  if (!(min > 0) || !(max > min)) return out;
  const e0 = Math.floor(Math.log10(min));
  const e1 = Math.ceil(Math.log10(max));
  const decades = Math.log10(max / min);
  const mult = decades > 3 ? [1] : decades > 1.3 ? [1, 2, 5] : [1, 1.5, 2, 3, 5, 7];
  for (let e = e0; e <= e1; e++) {
    const base = Math.pow(10, e);
    for (const m of mult) {
      const v = m * base;
      if (v >= min * (1 - 1e-9) && v <= max * (1 + 1e-9)) out.push(Number(v.toPrecision(6)));
    }
  }
  return out;
}
