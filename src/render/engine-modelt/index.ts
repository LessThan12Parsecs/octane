/**
 * Ford Model T mechanism (procedural three.js). Importing this module registers the builder for
 * EngineSpec.id 'ford-model-t' with createEngineModel (render/engine-model.ts); the render barrel
 * (render/index.ts) imports it for that side effect, so `createEngineModel(MODEL_T, op)` works for any
 * consumer of the barrel. registerModelTEngineModel() is idempotent for direct importers.
 */
import { registerEngineModel } from '../engine-model';
import { ModelTEngineModel } from './model';

export { ModelTEngineModel } from './model';
export { computeModelTLayout, pistonDisplacementAt, type ModelTLayout, type ModelTValveLayout } from './layout';
export { computeModelTPose, createModelTInput, createModelTPose, type ModelTKinematicInput, type ModelTPose } from './kinematics';
export { ValveCam, lobeCentreDeg, threeArcGeometry, threeArcLift, threeArcSupport, type ThreeArcGeometry } from './cam';

/** Engine id the Model T builder is registered for. */
export const MODEL_T_RENDER_ID = 'ford-model-t';

export function registerModelTEngineModel(): void {
  registerEngineModel(MODEL_T_RENDER_ID, (spec, op) => new ModelTEngineModel(spec, op));
}

registerModelTEngineModel();
