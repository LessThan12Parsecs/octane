/**
 * Engine-dependent UI configuration (pure, no DOM): which operating-point controls an engine shows,
 * their ranges and labels (EngineUiProfile from the engine registry), the load-model and gear choices,
 * and the definition to use for an ad-hoc spec.
 */
import { hasVariableCompressionRatio, type EngineSpec } from '../physics/core/engine-spec';
import { ENGINES, engineOfSpec, type EngineDefinition, type EngineUiProfile } from '../physics/engines/index';
import type { LoadModel, OpParams } from './op-binding';

/**
 * The registry definition of a spec, with the spec itself substituted (a registry spec modified by a
 * caller keeps its engine's profile). Specs without a registry id get the CFR profile adapted to what
 * the spec can do (CR control only with a CR range, dwell only for inductive ignition, MAG/BAT only for
 * a magneto).
 */
export function definitionForSpec(spec: EngineSpec): EngineDefinition {
  const def = engineOfSpec(spec);
  if (def) return def.spec === spec ? def : { ...def, spec };
  const base = ENGINES['cfr-f1'];
  return {
    ...base,
    label: spec.name,
    description: spec.name,
    spec,
    presets: [],
    ui: {
      ...base.ui,
      showCompressionRatio: hasVariableCompressionRatio(spec),
      showDwell: spec.ignition.type === 'inductive',
      showIgnitionSource: spec.ignition.type === 'trembler-magneto',
      loadModels: spec.vehicle ? ['constant', 'brake', 'vehicle'] : ['constant', 'brake'],
    },
  };
}

/** Window / controls title for an engine. */
export function engineTitle(def: Pick<EngineDefinition, 'label'>): string {
  return `Octane · ${def.label}`;
}

/** Whether the engine uses the CFR octane-rating presets (RON/MON with the MON spark-vs-CR schedule). */
export function usesCfrRatingPresets(def: Pick<EngineDefinition, 'id' | 'spec'>): boolean {
  return def.id === 'cfr-f1' && hasVariableCompressionRatio(def.spec);
}

export const LOAD_MODEL_LABELS: Readonly<Record<LoadModel, string>> = {
  constant: 'Constant torque',
  brake: 'Brake (T ∝ nᵏ)',
  vehicle: 'Vehicle (road load)',
};

/** Dropdown label → load model, for the models the profile offers. */
export function loadModelOptions(profile: Pick<EngineUiProfile, 'loadModels'>): Record<string, LoadModel> {
  const out: Record<string, LoadModel> = {};
  for (const m of profile.loadModels) out[LOAD_MODEL_LABELS[m]] = m;
  return out;
}

/** Dropdown label → gear name: the spec's vehicle gears (with their ratios) plus neutral (declutched). */
export function gearOptions(spec: EngineSpec): Record<string, string> {
  const out: Record<string, string> = {};
  const gears = spec.vehicle?.gears ?? {};
  for (const [name, ratio] of Object.entries(gears)) {
    const label = `${name.charAt(0).toUpperCase()}${name.slice(1)} (${Number(ratio.toFixed(3))}:1)`;
    out[label] = name;
  }
  out['Neutral (declutched)'] = 'neutral';
  return out;
}

/** Which operating-point controls are visible. */
export interface ControlVisibility {
  compressionRatio: boolean;
  dwell: boolean;
  ignitionSource: boolean;
  octaneNumber: boolean;
  loadModel: boolean;
  loadTorque: boolean;
  gear: boolean;
  grade: boolean;
  brake: boolean;
}

/** Visible controls for an engine and the current control values. */
export function controlVisibility(def: Pick<EngineDefinition, 'spec' | 'ui'>, p: Pick<OpParams, 'speedMode' | 'loadModel' | 'gear' | 'fuel'>): ControlVisibility {
  const ui = def.ui;
  const free = p.speedMode === 'free';
  const model: LoadModel = ui.loadModels.includes(p.loadModel) ? p.loadModel : 'constant';
  const vehicle = free && model === 'vehicle' && !!def.spec.vehicle;
  return {
    compressionRatio: ui.showCompressionRatio && hasVariableCompressionRatio(def.spec),
    dwell: ui.showDwell,
    ignitionSource: ui.showIgnitionSource,
    octaneNumber: p.fuel === 'PRF',
    loadModel: free && ui.loadModels.length > 1,
    loadTorque: free && !vehicle,
    gear: vehicle,
    grade: vehicle && p.gear !== 'neutral',
    brake: free && model === 'brake',
  };
}
