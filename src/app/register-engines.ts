/**
 * Side-effect module: makes the render models of the non-CFR engines available.
 *
 * Each engine's 3D model registers itself with render/engine-model.ts (registerEngineModel) when its
 * module is evaluated. Evaluating every `src/render/engine-<name>/index.ts` here means an engine module
 * only has to exist to be picked up; without one, createEngineModel throws for that engine and the app
 * reports it and stays on (or returns to) a working engine. A module the render barrel already imports
 * is evaluated once either way.
 */
import.meta.glob('../render/engine-*/index.ts', { eager: true });

export {};
