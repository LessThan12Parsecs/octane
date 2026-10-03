/**
 * Perceptually uniform 'inferno' colormap (van der Walt & Smith, matplotlib,
 * CC0) and the temperature-legend helpers exported for the UI.
 *
 * Table: every 5th entry of matplotlib 3.11 `_cm_listed._inferno_data`
 * (indices 0, 5, …, 255 → 52 knots at x = k/51), sRGB-encoded.
 * Linear interpolation between knots deviates < 0.004 from the full table.
 */
import { VIS_TEMPERATURE_RANGE } from '../constants';

const INFERNO_KNOTS: readonly (readonly [number, number, number])[] = [
  [0.001462, 0.000466, 0.013866], [0.007676, 0.006136, 0.046836], [0.019373, 0.015133, 0.088767],
  [0.037668, 0.025921, 0.132232], [0.06134, 0.03659, 0.177642], [0.087411, 0.044556, 0.224813],
  [0.116656, 0.047574, 0.272321], [0.149073, 0.045468, 0.317085], [0.183429, 0.040329, 0.354971],
  [0.217949, 0.036615, 0.383522], [0.25162, 0.037705, 0.403378], [0.284321, 0.043933, 0.416608],
  [0.316282, 0.05349, 0.425116], [0.347771, 0.064616, 0.430217], [0.379001, 0.076253, 0.432719],
  [0.410113, 0.087896, 0.433098], [0.441207, 0.099338, 0.431594], [0.472328, 0.110547, 0.428334],
  [0.503493, 0.121575, 0.423356], [0.534683, 0.132534, 0.416667], [0.565854, 0.143567, 0.408258],
  [0.59694, 0.154848, 0.398125], [0.627847, 0.166575, 0.386276], [0.658463, 0.178962, 0.372748],
  [0.688653, 0.192239, 0.357603], [0.718264, 0.206636, 0.340931], [0.747127, 0.222378, 0.322856],
  [0.775059, 0.239667, 0.303526], [0.801871, 0.258674, 0.283099], [0.827372, 0.279517, 0.26175],
  [0.851384, 0.30226, 0.239636], [0.873741, 0.326906, 0.216886], [0.894305, 0.353399, 0.193584],
  [0.912966, 0.381636, 0.169755], [0.929644, 0.411479, 0.145367], [0.944285, 0.442772, 0.120354],
  [0.956852, 0.475356, 0.094695], [0.967322, 0.509078, 0.068659], [0.975677, 0.543798, 0.043618],
  [0.981895, 0.579392, 0.02625], [0.985952, 0.61575, 0.025592], [0.987819, 0.652773, 0.045581],
  [0.987464, 0.690366, 0.07999], [0.984865, 0.728427, 0.120785], [0.980032, 0.766837, 0.166353],
  [0.973088, 0.805409, 0.216877], [0.964394, 0.843848, 0.273391], [0.954997, 0.881569, 0.337475],
  [0.947594, 0.917399, 0.410665], [0.947937, 0.949318, 0.491426], [0.961812, 0.975924, 0.571925],
  [0.988362, 0.998364, 0.644924],
];

/** Inferno at x ∈ [0, 1] (clamped), sRGB-encoded components in [0, 1]. */
export function inferno(x: number, out: [number, number, number] = [0, 0, 0]): [number, number, number] {
  const n = INFERNO_KNOTS.length - 1;
  const t = Math.min(Math.max(Number.isFinite(x) ? x : 0, 0), 1) * n;
  const i = Math.min(Math.floor(t), n - 1);
  const w = t - i;
  const a = INFERNO_KNOTS[i], b = INFERNO_KNOTS[i + 1];
  out[0] = a[0] + (b[0] - a[0]) * w;
  out[1] = a[1] + (b[1] - a[1]) * w;
  out[2] = a[2] + (b[2] - a[2]) * w;
  return out;
}

/** sRGB transfer (IEC 61966-2-1) decode: encoded → linear. */
export function srgbToLinear(c: number): number {
  return c <= 0.04045 ? c / 12.92 : Math.pow((c + 0.055) / 1.055, 2.4);
}

/** Map a temperature to [0, 1] on the legend's linear scale. */
export function temperatureToUnit(T: number, range: readonly [number, number] = VIS_TEMPERATURE_RANGE): number {
  const [lo, hi] = range;
  return Math.min(Math.max((T - lo) / (hi - lo), 0), 1);
}

/** False colour of a gas temperature (sRGB-encoded). */
export function temperatureColor(T: number, range: readonly [number, number] = VIS_TEMPERATURE_RANGE): [number, number, number] {
  return inferno(temperatureToUnit(T, range));
}

/**
 * 256×1 RGBA8 sRGB-encoded LUT (for a GPU texture; mark it SRGBColorSpace so
 * the sampler decodes to linear).
 */
export function infernoLUT(size = 256): Uint8Array {
  const data = new Uint8Array(size * 4);
  const c: [number, number, number] = [0, 0, 0];
  for (let i = 0; i < size; i++) {
    inferno(i / (size - 1), c);
    data[4 * i] = Math.round(c[0] * 255);
    data[4 * i + 1] = Math.round(c[1] * 255);
    data[4 * i + 2] = Math.round(c[2] * 255);
    data[4 * i + 3] = 255;
  }
  return data;
}

// ---------------------------------------------------------------------------
// Legend helper for the UI
// ---------------------------------------------------------------------------

export interface TemperatureLegend {
  /** Unit of the values. */
  unit: 'K';
  min: number;
  max: number;
  /** Colour stops along the bar, position 0..1 → CSS colour. */
  stops: { position: number; color: string }[];
  /** Ready-to-use CSS `linear-gradient(...)` (left → right = min → max). */
  cssGradient: string;
  /** Suggested tick values (K) and their positions 0..1. */
  ticks: { value: number; position: number; label: string }[];
  /** Short title for the legend. */
  title: string;
}

function hex(c: [number, number, number]): string {
  const h = (v: number) => Math.round(Math.min(Math.max(v, 0), 1) * 255).toString(16).padStart(2, '0');
  return `#${h(c[0])}${h(c[1])}${h(c[2])}`;
}

/**
 * Data for drawing the 'temperature' mode colour legend (matches the shader
 * exactly: linear in T over `range`, inferno).
 */
export function temperatureLegend(
  range: readonly [number, number] = VIS_TEMPERATURE_RANGE,
  nStops = 16,
): TemperatureLegend {
  const [min, max] = range;
  const stops: TemperatureLegend['stops'] = [];
  for (let i = 0; i <= nStops; i++) {
    const p = i / nStops;
    stops.push({ position: p, color: hex(inferno(p)) });
  }
  const cssGradient = `linear-gradient(to right, ${stops.map((s) => `${s.color} ${(s.position * 100).toFixed(1)}%`).join(', ')})`;
  const step = niceStep((max - min) / 5);
  const ticks: TemperatureLegend['ticks'] = [];
  for (let v = Math.ceil(min / step) * step; v <= max + 1e-9; v += step) {
    ticks.push({ value: v, position: (v - min) / (max - min), label: `${Math.round(v)}` });
  }
  return { unit: 'K', min, max, stops, cssGradient, ticks, title: 'Gas temperature (K)' };
}

function niceStep(raw: number): number {
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  const m = raw / p;
  return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p;
}
