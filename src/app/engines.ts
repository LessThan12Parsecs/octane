/** Engine selection for the app (pure helpers; render-model registration is in register-engines.ts). */
import { ENGINE_IDS, ENGINES, getEngine, type EngineDefinition } from '../physics/engines/index';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineChoice } from './overlay';
import { clampOperatingPoint, hasOpKnobs, parseUrlOptions } from './sync';

/** Engines the picker offers, in registry order. */
export function engineChoices(): EngineChoice[] {
  return ENGINE_IDS.map((id) => ({ id, label: ENGINES[id].label, description: ENGINES[id].description }));
}

/** True for a registry id. */
export function isEngineId(id: string | null | undefined): boolean {
  return !!id && Object.prototype.hasOwnProperty.call(ENGINES, id);
}

/**
 * The engine to start with: the first registry id among the candidates (explicit option, ?engine=,
 * remembered preference), else the default engine.
 */
export function resolveEngine(...candidates: (string | null | undefined)[]): EngineDefinition {
  for (const c of candidates) if (isEngineId(c)) return getEngine(c);
  return getEngine(null);
}

/**
 * The engine the app starts with for the page's `search`: the explicit option, ?engine=, the remembered
 * choice, else the default engine — except that a URL carrying operating-point knobs (?cr= ?on= ?rpm=
 * ?spark= ?phi=) without a registry ?engine= ignores the remembered choice. Such links predate the engine
 * picker and were written for the default engine (the CFR F-1); applied to whichever engine this browser
 * ran last they would describe something else (e.g. the Model T with its CR clamped to 3.98).
 */
export function resolveStartEngine(explicit: string | null | undefined, search: string, remembered: string | null | undefined): EngineDefinition {
  const url = parseUrlOptions(search);
  const legacyKnobs = !isEngineId(url.engine) && hasOpKnobs(search);
  return resolveEngine(explicit, url.engine, legacyKnobs ? null : remembered);
}

/**
 * Initial operating point of an engine: its default with the URL overrides on top (only for the engine
 * the app starts with; a later switch starts from the new engine's defaults), CR clamped to the spec.
 */
export function initialOperatingPoint(def: EngineDefinition, overrides: Partial<OperatingPoint> = {}, base?: OperatingPoint): OperatingPoint {
  return clampOperatingPoint(def.spec, { ...(base ?? def.defaultOperatingPoint), ...overrides });
}
