/**
 * Render layer barrel: the 3D stage, the engine mechanism and the in-cylinder
 * combustion visuals. Everything here consumes EngineSnapshot / EngineSpec only.
 */
export { Stage, type CameraFraming, type Grading, type StageOptions } from './scene';
export { EngineModel, recommendedCameraView, type CameraView, type EngineLayout, type EnginePose } from './engine/index';
export {
  chamberDepth,
  chamberShapeOf,
  CombustionVisuals,
  createTemperatureLegendElement,
  CylinderSnapshotView,
  temperatureLegend,
  VIS_TEMPERATURE_RANGE,
  type ChamberShape,
  type CombustionMode,
  type CombustionVisualsOptions,
  type CutPlanes,
  type TemperatureLegendElement,
} from './combustion/index';
export {
  CfrEngineRenderModel,
  createEngineModel,
  registerEngineModel,
  type EngineRenderModel,
} from './engine-model';
// Side effect: registers the Ford Model T builder with createEngineModel (spec id 'ford-model-t').
import './engine-modelt/index';
export { ModelTEngineModel, computeModelTLayout, type ModelTLayout } from './engine-modelt/index';
