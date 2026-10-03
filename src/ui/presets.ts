/**
 * Standard CFR rating conditions (ASTM D2699 Research, D2700 Motor).
 *
 * If src/physics/engines/cfr.ts exports `CFR_RON_CONDITIONS` / `CFR_MON_CONDITIONS`
 * (Partial<OperatingPoint>), those win; otherwise the local defaults below are used.
 * The compression ratio is deliberately NOT part of a preset: in a rating it is the
 * variable adjusted to reach standard knock intensity for the fuel under test.
 */
import * as cfrModule from '../physics/engines/cfr';
import type { OperatingPoint } from '../physics/core/operating-point';
import { P_ATM } from '../physics/core/constants';
import { celsiusToKelvin, relativeHumidityForHumidityRatio } from './units';

export type PresetId = 'RON' | 'MON';

/** Ambient temperature assumed by the presets, K. */
const PRESET_AMBIENT_T = celsiusToKelvin(25);

/**
 * D2699/D2700 specify intake-air humidity as 3.56–7.12 g water / kg dry air; we use
 * the mid-point and express it as relative humidity at the preset ambient state.
 * UNVERIFIED: humidity band quoted from memory of ASTM D2699 Table (engine conditions).
 */
const PRESET_HUMIDITY_RATIO = 5.34e-3;

/**
 * Local RON (ASTM D2699) conditions: 600 rpm, 13° BTDC fixed spark, intake air 52 °C
 * (at standard barometric pressure), coolant 100 °C, WOT, fuel–air ratio for maximum knock.
 * UNVERIFIED: values from memory of ASTM D2699; the mixture temperature is not
 * specified by the method — we approximate the port mixture temperature with the
 * intake-air temperature (carburettor evaporative cooling is not separated here).
 */
export const LOCAL_RON_CONDITIONS: Partial<OperatingPoint> = {
  speedMode: 'fixed',
  rpm: 600,
  throttle: 1,
  sparkAdvanceDeg: 13,
  intakeMixtureTemperature: celsiusToKelvin(52),
  ambientPressure: P_ATM,
  ambientTemperature: PRESET_AMBIENT_T,
  relativeHumidity: relativeHumidityForHumidityRatio(PRESET_HUMIDITY_RATIO, PRESET_AMBIENT_T, P_ATM),
  // Fuel–air ratio is tuned for maximum knock in a rating; slightly rich for PRFs.
  equivalenceRatio: 1.05,
  egrFraction: 0,
  coolantTemperature: celsiusToKelvin(100),
};

/**
 * Local MON (ASTM D2700) conditions: 900 rpm, intake air 38 °C, mixture 149 °C,
 * coolant 100 °C, spark advance scheduled with compression ratio (see monSparkAdvance).
 * UNVERIFIED: values from memory of ASTM D2700.
 */
export const LOCAL_MON_CONDITIONS: Partial<OperatingPoint> = {
  speedMode: 'fixed',
  rpm: 900,
  throttle: 1,
  intakeMixtureTemperature: celsiusToKelvin(149),
  ambientPressure: P_ATM,
  ambientTemperature: PRESET_AMBIENT_T,
  relativeHumidity: relativeHumidityForHumidityRatio(PRESET_HUMIDITY_RATIO, PRESET_AMBIENT_T, P_ATM),
  equivalenceRatio: 1.05,
  egrFraction: 0,
  coolantTemperature: celsiusToKelvin(100),
};

/**
 * MON basic spark timing vs compression ratio: 26° BTDC at low CR falling to 14° BTDC
 * at high CR (D2700 ties timing to cylinder height). Linear between (CR 5, 26°) and
 * (CR 14, 14°), clamped.
 * UNVERIFIED: end points approximate the D2700 schedule; replace with the tabulated
 * relation when available.
 */
export function monSparkAdvance(compressionRatio: number): number {
  const cr0 = 5;
  const cr1 = 14;
  const a0 = 26;
  const a1 = 14;
  const f = Math.min(1, Math.max(0, (compressionRatio - cr0) / (cr1 - cr0)));
  return a0 + f * (a1 - a0);
}

function isOpPatch(v: unknown): v is Partial<OperatingPoint> {
  return typeof v === 'object' && v !== null && !Array.isArray(v);
}

/** Read an optional export from cfr.ts without a hard compile-time dependency on it. */
function fromCfr(name: string): Partial<OperatingPoint> | null {
  const v = (cfrModule as unknown as Record<string, unknown>)[name];
  return isOpPatch(v) ? v : null;
}

export interface PresetInfo {
  id: PresetId;
  label: string;
  /** Operating-point patch (compression ratio excluded unless the source sets it). */
  patch: Partial<OperatingPoint>;
  /** Whether spark advance follows the MON compression-ratio schedule. */
  sparkFollowsCR: boolean;
  source: 'cfr.ts' | 'local';
}

/**
 * Preset patch for `id`. For MON the spark advance is computed for the given CR
 * unless the cfr.ts export provides one.
 */
export function getPreset(id: PresetId, compressionRatio: number): PresetInfo {
  if (id === 'RON') {
    const ext = fromCfr('CFR_RON_CONDITIONS');
    return {
      id,
      label: 'CFR Research method (RON)',
      patch: { ...LOCAL_RON_CONDITIONS, ...(ext ?? {}) },
      sparkFollowsCR: false,
      source: ext ? 'cfr.ts' : 'local',
    };
  }
  const ext = fromCfr('CFR_MON_CONDITIONS');
  const patch: Partial<OperatingPoint> = { ...LOCAL_MON_CONDITIONS, ...(ext ?? {}) };
  const sparkFollowsCR = ext?.sparkAdvanceDeg === undefined;
  if (sparkFollowsCR) patch.sparkAdvanceDeg = monSparkAdvance(patch.compressionRatio ?? compressionRatio);
  return { id, label: 'CFR Motor method (MON)', patch, sparkFollowsCR, source: ext ? 'cfr.ts' : 'local' };
}
