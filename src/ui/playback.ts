/**
 * Playback time-scale helpers. The time scale is simulated seconds per wall-clock
 * second: 1 = real time, 1/5000 = the slowest (at 600 rpm one crank degree then
 * takes ~1.4 s of wall time, enough to watch the spark discharge and flame kernel).
 */

export const MAX_TIME_SCALE = 1;
export const MIN_TIME_SCALE = 1 / 5000;
export const STEP_SMALL_DEG = 0.5;
export const STEP_LARGE_DEG = 5;
export const DEFAULT_TIME_SCALE = 1 / 50;

/** 1-2-5 ladder of slow-down factors used by the [ / ] keys. */
export const TIME_SCALE_LADDER: readonly number[] = [
  1, 1 / 2, 1 / 5, 1 / 10, 1 / 20, 1 / 50, 1 / 100, 1 / 200, 1 / 500, 1 / 1000, 1 / 2000, 1 / 5000,
];

/** Round the slow-down factor 1/ts to two significant digits so the readout is exact. */
export function snapTimeScale(ts: number): number {
  if (!(ts > 0)) return MIN_TIME_SCALE;
  const clamped = Math.min(MAX_TIME_SCALE, Math.max(MIN_TIME_SCALE, ts));
  const n = 1 / clamped;
  const mag = Math.pow(10, Math.floor(Math.log10(n)) - 1);
  const snapped = Math.max(1, Math.round(n / mag) * mag);
  return Math.min(MAX_TIME_SCALE, Math.max(MIN_TIME_SCALE, 1 / Number(snapped.toPrecision(2))));
}

/**
 * Logarithmic slider mapping: position 0 → MIN_TIME_SCALE, 1 → real time.
 * The result is snapped (see snapTimeScale).
 */
export function sliderToTimeScale(u: number): number {
  const x = Math.min(1, Math.max(0, u));
  const lnMin = Math.log(MIN_TIME_SCALE);
  const lnMax = Math.log(MAX_TIME_SCALE);
  return snapTimeScale(Math.exp(lnMin + x * (lnMax - lnMin)));
}

export function timeScaleToSlider(ts: number): number {
  const lnMin = Math.log(MIN_TIME_SCALE);
  const lnMax = Math.log(MAX_TIME_SCALE);
  const c = Math.min(MAX_TIME_SCALE, Math.max(MIN_TIME_SCALE, ts));
  return (Math.log(c) - lnMin) / (lnMax - lnMin);
}

/** "1×" (real time) or "1/250×". */
export function formatTimeScale(ts: number): string {
  if (!(ts > 0)) return '—';
  if (ts >= 0.999) return '1× (real time)';
  const n = 1 / ts;
  const nStr = n >= 100 ? Math.round(n).toString() : Number(n.toPrecision(2)).toString();
  return `1/${nStr}×`;
}

/** Next ladder step: dir = +1 faster, −1 slower. */
export function stepTimeScale(ts: number, dir: 1 | -1): number {
  const L = TIME_SCALE_LADDER;
  if (dir > 0) {
    for (let i = L.length - 1; i >= 0; i--) if (L[i] > ts * (1 + 1e-9)) return L[i];
    return L[0];
  }
  for (let i = 0; i < L.length; i++) if (L[i] < ts * (1 - 1e-9)) return L[i];
  return L[L.length - 1];
}

/** Wall-clock seconds needed to play `deg` crank degrees at `rpm` and time scale `ts`. */
export function wallSecondsPerDegrees(deg: number, rpm: number, ts: number): number {
  if (!(rpm > 0) || !(ts > 0)) return NaN;
  return deg / (6 * rpm) / ts;
}
