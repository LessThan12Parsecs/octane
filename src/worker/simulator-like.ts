/**
 * The minimal interface the simulation worker needs from a simulator.
 *
 * Implemented today by `MockSimulator` (./mock-simulator.ts) and later by the
 * real physics `EngineSimulator` (src/physics/cycle). The worker host
 * (./sim.worker.ts) only ever talks to this interface, so swapping the
 * implementation is a one-line change of the factory there.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import type { SimulatorOptions } from './protocol';

export interface SimulatorLike {
  /**
   * Change operating conditions. rpm / throttle / load apply immediately, the
   * rest at the start of the next cycle (θ = −360°).
   */
  setOperatingPoint(patch: Partial<OperatingPoint>): void;
  /** Restart from t = 0 with the current (latest) operating point. */
  reset(): void;
  /**
   * Advance the simulation to the next output instant and return a fresh
   * snapshot there. The first call after construction / reset returns the
   * initial state at t = 0. Returned objects are never mutated afterwards.
   */
  advanceToNextSnapshot(): EngineSnapshot;
  /** Cycle summaries completed since the previous drain (oldest first). */
  drainCycleSummaries(): CycleSummary[];
  /** Simulated time of the most recently returned snapshot, s (0 before the first). */
  readonly time: number;
}

/** Constructs a simulator; the worker host is parameterised by one of these. */
export type SimulatorFactory = (
  spec: EngineSpec,
  operatingPoint: OperatingPoint,
  options: SimulatorOptions,
) => SimulatorLike;
