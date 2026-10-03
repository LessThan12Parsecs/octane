/**
 * Unit conversions and number formatting for the UI. The physics side is SI
 * everywhere; everything the user sees goes through these helpers so the
 * conversions live in exactly one (tested) place.
 */

// ---------------------------------------------------------------- conversions

export const ZERO_CELSIUS_K = 273.15;
export const PA_PER_BAR = 1e5;
export const PA_PER_KPA = 1e3;
/** J per kWh. */
export const J_PER_KWH = 3.6e6;

export const kelvinToCelsius = (k: number): number => k - ZERO_CELSIUS_K;
export const celsiusToKelvin = (c: number): number => c + ZERO_CELSIUS_K;
export const paToBar = (pa: number): number => pa / PA_PER_BAR;
export const barToPa = (bar: number): number => bar * PA_PER_BAR;
export const paToKPa = (pa: number): number => pa / PA_PER_KPA;
export const kPaToPa = (kpa: number): number => kpa * PA_PER_KPA;
export const mToMm = (m: number): number => m * 1e3;
export const m3ToCm3 = (m3: number): number => m3 * 1e6;
export const kgToMg = (kg: number): number => kg * 1e6;
export const kgPerSToGPerS = (kgs: number): number => kgs * 1e3;
export const voltsToKV = (v: number): number => v * 1e-3;
export const ampsToMA = (a: number): number => a * 1e3;
export const jToMJ = (j: number): number => j * 1e3;
export const sToMs = (s: number): number => s * 1e3;
export const msToS = (ms: number): number => ms * 1e-3;
export const sToUs = (s: number): number => s * 1e6;
export const fractionToPct = (f: number): number => f * 100;
export const pctToFraction = (p: number): number => p / 100;

/** Indicated specific fuel consumption: kg/J → g/kWh. */
export const isfcToGPerKWh = (kgPerJ: number): number => kgPerJ * 1e3 * J_PER_KWH;

/** Crank degrees per second at `rpm` (360°/rev × rpm/60). */
export const degPerSecond = (rpm: number): number => 6 * rpm;

/** Duration of one crank degree, s. NaN for rpm ≤ 0. */
export const secondsPerDeg = (rpm: number): number => (rpm > 0 ? 1 / (6 * rpm) : NaN);

/** Duration of a 720° four-stroke cycle, s. NaN for rpm ≤ 0. */
export const cyclePeriod = (rpm: number): number => (rpm > 0 ? 120 / rpm : NaN);

/**
 * Heat-release rate W → J per crank degree (the conventional engine unit):
 * dQ/dθ = Q̇ / (dθ/dt). NaN when the engine is not turning.
 */
export const wattsToJPerDeg = (w: number, rpm: number): number => (rpm > 0 ? w / (6 * rpm) : NaN);

// -------------------------------------------------------------- psychrometry

/**
 * Saturation vapour pressure of water over liquid, Pa.
 * Buck (1981/1996 revision): e_s = 6.1121 exp((18.678 − T/234.5)(T/(257.14 + T))) hPa, T in °C.
 * A. L. Buck, "New equations for computing vapor pressure and enhancement factor",
 * J. Appl. Meteorol. 20 (1981) 1527–1532 (1996 coefficients).
 */
export function saturationVapourPressure(tK: number): number {
  const t = tK - ZERO_CELSIUS_K;
  return 611.21 * Math.exp((18.678 - t / 234.5) * (t / (257.14 + t)));
}

/** Molar-mass ratio M_H2O / M_air (18.01528 / 28.9647). */
export const EPSILON_WATER_AIR = 18.01528 / 28.9647;

/** Humidity ratio (kg water / kg dry air) for relative humidity `rh` (0..1) at T, p. */
export function humidityRatio(rh: number, tK: number, p: number): number {
  const pv = rh * saturationVapourPressure(tK);
  return (EPSILON_WATER_AIR * pv) / (p - pv);
}

/** Relative humidity (0..1) that gives humidity ratio `w` (kg/kg dry air) at T, p. */
export function relativeHumidityForHumidityRatio(w: number, tK: number, p: number): number {
  const pv = (w * p) / (EPSILON_WATER_AIR + w);
  return pv / saturationVapourPressure(tK);
}

// ----------------------------------------------------------------- formatting

export const MISSING = '—';

/** Fixed-decimal formatting with a dash for NaN/±Inf, and no "-0.0". */
export function fmt(v: number, decimals: number, fallback: string = MISSING): string {
  if (!Number.isFinite(v)) return fallback;
  const s = v.toFixed(decimals);
  // Avoid "-0.00" for tiny negative values.
  if (s.charCodeAt(0) === 45 /* - */ && Number(s) === 0) return s.slice(1);
  return s.replace('-', '−');
}

/** Like fmt but always shows a sign (+/−). */
export function fmtSigned(v: number, decimals: number, fallback: string = MISSING): string {
  if (!Number.isFinite(v)) return fallback;
  const s = fmt(v, decimals);
  return s.charCodeAt(0) === 0x2212 || Number(v.toFixed(decimals)) === 0 ? s : `+${s}`;
}

/** Significant-figure formatting (for values spanning decades), e.g. 0.00123 → "0.00123". */
export function fmtSig(v: number, sig: number, fallback: string = MISSING): string {
  if (!Number.isFinite(v)) return fallback;
  if (v === 0) return '0';
  const mag = Math.floor(Math.log10(Math.abs(v)));
  const decimals = Math.max(0, sig - 1 - mag);
  return fmt(v, Math.min(decimals, 12));
}

/** A duration given in microseconds, labelled µs below 1 ms and ms above. */
export function fmtMicros(us: number): string {
  if (!Number.isFinite(us)) return MISSING;
  const a = Math.abs(us);
  if (a === 0) return '0';
  if (a < 1000) return `${fmt(us, a < 10 ? 1 : 0)} µs`;
  const ms = us / 1000;
  return `${fmt(ms, Math.abs(ms) < 10 ? 2 : 1)} ms`;
}

/** Compact axis tick for a time in µs: "0", "500 µs", "1.5 ms" (thin space before the unit). */
export function fmtMicrosTick(us: number): string {
  if (!Number.isFinite(us)) return '';
  if (Math.abs(us) < 1e-9) return '0';
  const minus = (s: string): string => s.replace('-', '\u2212');
  if (Math.abs(us) < 1000) return `${minus(Math.round(us).toString())}\u2009µs`;
  return `${minus(Number((us / 1000).toPrecision(3)).toString())}\u2009ms`;
}

/** Simulated time in s → "12.345 s" / "85.2 ms". */
export function fmtSimTime(t: number): string {
  if (!Number.isFinite(t)) return MISSING;
  if (Math.abs(t) < 1) return `${fmt(t * 1e3, 2)} ms`;
  return `${fmt(t, 3)} s`;
}
