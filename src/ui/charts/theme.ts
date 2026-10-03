/**
 * Chart palette and shared styling.
 *
 * Series colours are the dark-surface steps of a validated colour-blind-safe
 * categorical palette (blue, orange, aqua, yellow, magenta, violet). Checked with the
 * dataviz palette validator against the panel surface #12161c: every chart uses at
 * most three hues together in the validated order (all-pairs CVD ΔE ≥ 9), and the
 * four stroke colours pass the cyclic-adjacency check (intake↔compression↔power↔
 * exhaust↔intake). Identity is never colour-alone: every chart has a legend.
 */
import type { StrokeName } from '../engine-cycle';

export const SURFACE = '#12161c';

export const INK = {
  primary: '#e6e9ee',
  secondary: '#a3acb9',
  muted: '#6b7480',
  grid: 'rgba(255,255,255,0.06)',
  tick: 'rgba(255,255,255,0.16)',
  axis: 'rgba(255,255,255,0.22)',
} as const;

export const SERIES = {
  blue: '#3987e5',
  orange: '#d95926',
  aqua: '#199e70',
  yellow: '#c98500',
  magenta: '#d55181',
  violet: '#9085e9',
} as const;

/** Status colours — reserved for states (always paired with a label or icon). */
export const STATUS = {
  good: '#0ca30c',
  warning: '#fab219',
  serious: '#ec835a',
  critical: '#d03b3b',
} as const;

export const STROKE_COLOR: Readonly<Record<StrokeName, string>> = {
  intake: SERIES.blue,
  compression: SERIES.aqua,
  power: SERIES.yellow,
  exhaust: SERIES.magenta,
};

/** Colour of previous-cycle traces. */
export const GHOST_RGB = '163, 172, 185';
/** Colour of knocking previous-cycle traces (status: serious). */
export const GHOST_KNOCK_RGB = '236, 131, 90';

export const FONT_UI = '11px ui-sans-serif, -apple-system, "Segoe UI", Roboto, sans-serif';
export const FONT_MONO = '10px ui-monospace, "SF Mono", Menlo, Consolas, monospace';

/** Convert "#rrggbb" to "rgba(r, g, b, a)". */
export function withAlpha(hex: string, a: number): string {
  const m = /^#?([0-9a-f]{2})([0-9a-f]{2})([0-9a-f]{2})$/i.exec(hex);
  if (!m) return hex;
  return `rgba(${parseInt(m[1], 16)}, ${parseInt(m[2], 16)}, ${parseInt(m[3], 16)}, ${a})`;
}
