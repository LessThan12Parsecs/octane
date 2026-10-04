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
  /** Engine registry id (?engine=ford-model-t); validated by the app against the registry. */
  engine?: string;
  /** Focus cylinder, 0-based (?cyl=2 → cylinder 2 → 1). */
  focusCylinder?: number;
  /** Simulator implementation (?sim=mock: the lightweight single-cylinder stand-in). */
  simulator?: 'physics' | 'mock';
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
  const engine = q.get('engine')?.trim();
  if (engine && /^[a-z0-9][a-z0-9-]*$/i.test(engine)) out.engine = engine;
  const cyl = num(q, 'cyl');
  if (cyl !== undefined && Number.isInteger(cyl) && cyl >= 1) out.focusCylinder = cyl - 1;
  const sim = q.get('sim');
  if (sim === 'mock' || sim === 'physics') out.simulator = sim;
  return out;
}

/**
 * `search` with ?engine=<id> set (other parameters kept, in order). Operating-point knobs that belong to
 * the previous engine (cr, on, rpm, spark, phi) are dropped when `dropOpKnobs`, so a reload starts the new
 * engine from its own defaults. Returns the search string including its leading '?' ('' when empty).
 */
export function searchWithEngine(search: string, id: string, dropOpKnobs = true): string {
  const q = new URLSearchParams(search);
  if (dropOpKnobs) for (const k of OP_KNOBS) q.delete(k);
  q.delete('cyl');
  q.set('engine', id);
  const s = q.toString();
  return s ? `?${s}` : '';
}

/** URL keys that set operating-point fields (parseUrlOptions). */
export const OP_KNOBS: readonly string[] = ['cr', 'on', 'rpm', 'spark', 'phi'];

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

/**
 * A short, readable form of an error message for the status card: stack frames ("    at …") and
 * repeated lines dropped, at most `maxLines` lines and `maxChars` characters (the console keeps the rest).
 */
export function errorSummary(message: string, maxLines = 3, maxChars = 360): string {
  const out: string[] = [];
  for (const raw of message.split('\n')) {
    const line = raw.trim();
    if (!line || /^at\s/.test(line)) continue;
    const bare = line.replace(/^Error:\s*/, '');
    if (out.some((l) => l.includes(bare))) continue;
    out.push(line);
    if (out.length >= maxLines) break;
  }
  const s = out.join('\n');
  return s.length > maxChars ? `${s.slice(0, maxChars - 1)}…` : s;
}

/** Clamp an operating point's CR into the spec's range. */
export function clampOperatingPoint(spec: EngineSpec, op: OperatingPoint): OperatingPoint {
  const [lo, hi] = spec.geometry.compressionRatioRange;
  return { ...op, compressionRatio: Math.min(hi, Math.max(lo, op.compressionRatio)) };
}
