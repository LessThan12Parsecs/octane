/**
 * Small pure helpers used by the app loop (unit-tested in node).
 */
import type { EngineSpec } from '../physics/core/engine-spec';
import type { FuelSelection, OperatingPoint } from '../physics/core/operating-point';
import type { EngineSnapshot } from '../physics/core/snapshot';

/**
 * Exact piston travel TDC→BDC of the spec's slider crank, m (equals the stroke
 * without a wrist-pin offset). Same expression as the physics kinematics.
 */
export function pistonTravel(spec: EngineSpec): number {
  const g = spec.geometry;
  const a = g.stroke / 2;
  const l = g.conRodLength;
  const e = g.pinOffset;
  return Math.sqrt((l + a) ** 2 - e * e) - Math.sqrt((l - a) ** 2 - e * e);
}

/**
 * Geometric compression ratio implied by a snapshot of a flat-disc chamber:
 *   h_TDC = clearanceHeight − pistonDisplacement,
 *   CR = 1 + V_d / (A·h_TDC + V_crevice)
 * (the crevice is part of the clearance volume, as in the physics kinematics and
 * EngineModel's layout). Linear interpolation of the two fields preserves their
 * difference, so interpolated snapshots give the same CR. NaN if the snapshot
 * carries no usable geometry.
 */
export function compressionRatioFromSnapshot(
  spec: EngineSpec,
  s: Pick<EngineSnapshot, 'clearanceHeight' | 'pistonDisplacement'>,
): number {
  const hTdc = s.clearanceHeight - s.pistonDisplacement;
  if (!(hTdc > 0) || !Number.isFinite(hTdc)) return NaN;
  const B = spec.geometry.bore;
  const A = (Math.PI * B * B) / 4;
  const vd = A * pistonTravel(spec);
  return 1 + vd / (A * hTdc + Math.max(0, spec.geometry.creviceVolume));
}

/**
 * Wall-clock frame timer: converts rAF timestamps (ms) into frame durations (s),
 * clamped to [0, maxSeconds] so a background tab or a debugger pause does not
 * produce a huge simulated jump.
 */
export class FrameClock {
  private last = NaN;

  constructor(readonly maxSeconds = 0.1) {}

  tick(nowMs: number): number {
    const prev = this.last;
    this.last = nowMs;
    if (!Number.isFinite(prev) || !Number.isFinite(nowMs)) return 0;
    const dt = (nowMs - prev) / 1000;
    return dt > 0 ? Math.min(dt, this.maxSeconds) : 0;
  }

  reset(): void {
    this.last = NaN;
  }
}

/**
 * Chooses the `sinceT` for SimClient.recentSnapshots each frame: the previous
 * playback time while playback moves forward; −∞ after a reset or a backwards
 * step (the UI's trace store skips samples it has already ingested, so
 * re-delivering history is harmless, while a gap would lose data).
 */
export class RecentWindow {
  private prev = NaN;

  next(playbackTime: number): number {
    const prev = this.prev;
    this.prev = playbackTime;
    if (!Number.isFinite(playbackTime) || !Number.isFinite(prev) || playbackTime < prev) return -Infinity;
    return prev;
  }

  reset(): void {
    this.prev = NaN;
  }
}

// ---------------------------------------------------------------------------
// URL options
// ---------------------------------------------------------------------------

export type Framing = 'chamber' | 'engine';

export interface UrlOptions {
  /** Operating-point overrides (?cr=7.2&on=95&rpm=900&spark=20&phi=1.0). */
  op: Partial<OperatingPoint>;
  /** Initial time scale (?ts=0.01). */
  timeScale?: number;
  /** Start paused (?paused=1). */
  paused?: boolean;
  /** Initial camera framing (?view=engine). */
  framing?: Framing;
}

const num = (q: URLSearchParams, key: string): number | undefined => {
  const raw = q.get(key);
  if (raw === null || raw.trim() === '') return undefined;
  const v = Number(raw);
  return Number.isFinite(v) ? v : undefined;
};

/** Parse the developer/sharing URL knobs; invalid values are ignored. */
export function parseUrlOptions(search: string): UrlOptions {
  const q = new URLSearchParams(search);
  const op: Partial<OperatingPoint> = {};
  const cr = num(q, 'cr');
  if (cr !== undefined && cr > 1) op.compressionRatio = cr;
  const on = num(q, 'on');
  if (on !== undefined && on >= 0 && on <= 120) op.fuel = { kind: 'PRF', octaneNumber: on } satisfies FuelSelection;
  const rpm = num(q, 'rpm');
  if (rpm !== undefined && rpm > 0) op.rpm = rpm;
  const spark = num(q, 'spark');
  if (spark !== undefined) op.sparkAdvanceDeg = spark;
  const phi = num(q, 'phi');
  if (phi !== undefined && phi > 0) op.equivalenceRatio = phi;

  const out: UrlOptions = { op };
  const ts = num(q, 'ts');
  if (ts !== undefined && ts > 0) out.timeScale = ts;
  const paused = q.get('paused');
  if (paused !== null) out.paused = paused !== '0' && paused !== 'false';
  const view = q.get('view');
  if (view === 'engine' || view === 'chamber') out.framing = view;
  return out;
}

/** Aspect ratio the recommended camera views are composed for. */
export const FRAMING_REFERENCE_ASPECT = 1.6;

/**
 * Camera distance multiplier that keeps a view composed for
 * FRAMING_REFERENCE_ASPECT horizontally in frame on a narrower viewport (same
 * vertical FOV): k = ref / aspect, clamped to [1, maxScale].
 */
export function framingDistanceScale(aspect: number, maxScale = 2.5): number {
  if (!(aspect > 0)) return 1;
  return Math.min(maxScale, Math.max(1, FRAMING_REFERENCE_ASPECT / aspect));
}

/** Clamp an operating point's CR into the spec's range. */
export function clampOperatingPoint(spec: EngineSpec, op: OperatingPoint): OperatingPoint {
  const [lo, hi] = spec.geometry.compressionRatioRange;
  return { ...op, compressionRatio: Math.min(hi, Math.max(lo, op.compressionRatio)) };
}
