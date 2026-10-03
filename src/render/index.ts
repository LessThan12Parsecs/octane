/**
 * Render layer barrel: the 3D stage, the engine mechanism and the in-cylinder
 * combustion visuals. Everything here consumes EngineSnapshot / EngineSpec only.
 */
export { Stage, type CameraFraming, type Grading, type StageOptions } from './scene';
export { EngineModel, recommendedCameraView, type CameraView, type EngineLayout, type EnginePose } from './engine/index';
export {
  CombustionVisuals,
  createTemperatureLegendElement,
  temperatureLegend,
  VIS_TEMPERATURE_RANGE,
  type CombustionMode,
  type TemperatureLegendElement,
} from './combustion/index';
