/**
 * Engine registry: every engine Octane can simulate, with its default operating point, named presets
 * and the control ranges the UI offers. Pure data (no cycle-model imports): the cycle model's
 * per-engine option defaults (calibration, friction, knock sensor, …) live in cycle/options.ts, keyed by
 * EngineSpec.id, so the worker derives them from the spec it is given.
 */
import type { EngineSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS } from './cfr';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_FORD_DYNO } from './model-t';

export type EngineId = 'cfr-f1' | 'ford-model-t';

export interface EnginePreset {
  id: string;
  label: string;
  /** Fields applied over the current operating point. */
  op: Partial<OperatingPoint>;
  /** One-line explanation (tooltip). */
  note?: string;
}

/** Control ranges and labels of the operating-point UI for one engine. */
export interface EngineUiProfile {
  /** Speed slider range, rev/min. */
  rpmRange: [number, number];
  /** Spark slider range, deg BTDC (Model T: the spark-lever / timer-make range). */
  sparkRange: [number, number];
  /** Spark control label ('Spark advance' / 'Spark lever'). */
  sparkLabel: string;
  /** Load-torque slider range, N m. */
  loadRange: [number, number];
  /** Throttle control label. */
  throttleLabel: string;
  /** Speed-mode labels for 'fixed' and 'free'. */
  speedModeLabels: { fixed: string; free: string };
  /** Load models offered in 'free' mode. */
  loadModels: ('constant' | 'brake' | 'vehicle')[];
  /** Show the dwell-time control (inductive ignition). */
  showDwell: boolean;
  /** Show the compression-ratio control (variable-CR engines). */
  showCompressionRatio: boolean;
  /** Show the ignition-source switch (MAG / BAT). */
  showIgnitionSource: boolean;
  /** PRF octane-number slider label. */
  fuelLabel: string;
}

export interface EngineDefinition {
  id: EngineId;
  /** Short display name. */
  label: string;
  /** One-line description for the engine picker. */
  description: string;
  spec: EngineSpec;
  defaultOperatingPoint: OperatingPoint;
  presets: EnginePreset[];
  ui: EngineUiProfile;
}

export const ENGINES: Readonly<Record<EngineId, EngineDefinition>> = {
  'cfr-f1': {
    id: 'cfr-f1',
    label: 'CFR F-1',
    description: 'Waukesha CFR octane-rating engine (1931–): one cylinder, variable compression ratio',
    spec: CFR_F1,
    defaultOperatingPoint: CFR_RON_CONDITIONS,
    presets: [
      { id: 'RON', label: 'RON (ASTM D2699)', op: CFR_RON_CONDITIONS, note: '600 rpm, 52 °C air, 13° BTDC' },
      { id: 'MON', label: 'MON (ASTM D2700)', op: CFR_MON_CONDITIONS, note: '900 rpm, 149 °C mixture, spark follows CR' },
    ],
    ui: {
      rpmRange: [200, 3000],
      sparkRange: [-10, 60],
      sparkLabel: 'Spark advance',
      loadRange: [-20, 100],
      throttleLabel: 'Throttle',
      speedModeLabels: { fixed: 'Fixed (synchronous motor)', free: 'Free (load torque)' },
      loadModels: ['constant', 'brake'],
      showDwell: true,
      showCompressionRatio: true,
      showIgnitionSource: false,
      fuelLabel: 'PRF octane number',
    },
  },
  'ford-model-t': {
    id: 'ford-model-t',
    label: 'Ford Model T',
    description: 'Ford Model T (1924–25): 2.9 L side-valve four, magneto and trembler coils',
    spec: MODEL_T,
    defaultOperatingPoint: MODEL_T_CRUISE,
    presets: [
      {
        id: 'cruise',
        label: 'Cruise, high gear',
        op: MODEL_T_CRUISE,
        note: 'Level road, part throttle, spark lever part-advanced, on magneto',
      },
      {
        id: 'full-throttle',
        label: 'Full throttle, high gear',
        op: { ...MODEL_T_CRUISE, throttle: 1, sparkAdvanceDeg: 40 },
        note: 'Accelerate toward top speed (≈ 42–45 mph)',
      },
      {
        id: 'hill-low',
        label: 'Hill in low gear',
        op: { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'low', grade: 0.1 }, throttle: 1, sparkAdvanceDeg: 20, rpm: 1200 },
        note: '10 % grade in low (2.75:1)',
      },
      {
        id: 'idle',
        label: 'Idle (≈ 400 rpm)',
        op: { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'neutral', grade: 0 }, rpm: 400, throttle: 0.06, sparkAdvanceDeg: 0 },
        note: 'Declutched, throttle nearly shut, spark retarded',
      },
      {
        id: 'ford-dyno',
        label: 'Ford dyno, WOT 1600 rpm',
        op: MODEL_T_FORD_DYNO,
        note: "Speed held by a brake, wide open — Ford's 1918 rating condition (20 hp at 1500–1600 rpm)",
      },
      {
        id: 'battery',
        label: 'Cruise on battery (BAT)',
        op: { ...MODEL_T_CRUISE, ignitionSource: 'battery' },
        note: '6 V battery instead of the magneto: slower coil build-up retards the spark at speed',
      },
    ],
    ui: {
      rpmRange: [150, 2400],
      sparkRange: [-15.5, 64.5],
      sparkLabel: 'Spark lever (timer make)',
      loadRange: [0, 150],
      throttleLabel: 'Hand throttle',
      speedModeLabels: { fixed: 'Fixed (dynamometer)', free: 'Free (driving)' },
      loadModels: ['vehicle', 'brake', 'constant'],
      showDwell: false,
      showCompressionRatio: false,
      showIgnitionSource: true,
      fuelLabel: 'Fuel octane (PRF surrogate)',
    },
  },
};

export const DEFAULT_ENGINE_ID: EngineId = 'cfr-f1';

export const ENGINE_IDS = Object.keys(ENGINES) as EngineId[];

/** Registry lookup; unknown ids fall back to the default engine. */
export function getEngine(id: string | null | undefined): EngineDefinition {
  return (id && (ENGINES as Record<string, EngineDefinition>)[id]) || ENGINES[DEFAULT_ENGINE_ID];
}

/** Registry entry of a spec (by spec.id), or undefined for ad-hoc specs. */
export function engineOfSpec(spec: EngineSpec): EngineDefinition | undefined {
  return spec.id ? (ENGINES as Record<string, EngineDefinition>)[spec.id] : undefined;
}
