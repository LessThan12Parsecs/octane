/**
 * Main thread <-> simulation worker message protocol.
 *
 * Flow control: the worker simulates ahead of the playback cursor. The main
 * thread plays snapshots back at `timeScale` × wall-clock (slow motion is the
 * normal way to watch combustion: at 600 rpm the whole burn lasts ~5 ms) and
 * periodically tells the worker how far playback has progressed with `demand`.
 * The worker keeps at most `bufferAheadSeconds` of simulated time buffered
 * beyond the last demanded time, then idles.
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';

export interface SimulatorOptions {
  /** Snapshot cadence in crank degrees (e.g. 0.5). */
  snapshotEveryDeg: number;
  /** Max simulated seconds to buffer ahead of playback. */
  bufferAheadSeconds: number;
}

export type ToWorker =
  | { type: 'init'; spec: EngineSpec; operatingPoint: OperatingPoint; options: SimulatorOptions }
  /** Change operating conditions; applied at the start of the next cycle (or immediately for rpm/throttle/load). */
  | { type: 'set-operating-point'; patch: Partial<OperatingPoint> }
  /** Playback has consumed snapshots up to simulated time `t`; keep buffering ahead of it. */
  | { type: 'demand'; t: number }
  | { type: 'reset' };

export type FromWorker =
  | { type: 'ready' }
  | { type: 'snapshots'; batch: EngineSnapshot[] }
  | { type: 'cycle'; summary: CycleSummary }
  | { type: 'log'; level: 'info' | 'warn' | 'error'; message: string };
