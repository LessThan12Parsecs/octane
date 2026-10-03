/**
 * Four-stroke cycle bookkeeping for the UI (pure, no DOM).
 *
 * Crank-angle convention (DESIGN.md): θ = 0 at firing TDC, one cycle spans
 * [-360°, 360°): intake −360→−180, compression −180→0, power 0→180,
 * exhaust 180→360.
 */
import type { EngineSpec, ValveSpec } from '../physics/core/engine-spec';

export type StrokeName = 'intake' | 'compression' | 'power' | 'exhaust';

export const STROKES: readonly StrokeName[] = ['intake', 'compression', 'power', 'exhaust'];

export const STROKE_LABEL: Readonly<Record<StrokeName, string>> = {
  intake: 'Intake',
  compression: 'Compression',
  power: 'Power',
  exhaust: 'Exhaust',
};

/** Wrap any crank angle into the cycle range [-360, 360). */
export function wrapCycleDeg(theta: number): number {
  if (!Number.isFinite(theta)) return theta;
  const w = (((theta + 360) % 720) + 720) % 720; // [0, 720)
  return w - 360;
}

/** Stroke containing crank angle θ (any value; wrapped first). */
export function strokeOf(theta: number): StrokeName {
  const w = wrapCycleDeg(theta);
  if (w < -180) return 'intake';
  if (w < 0) return 'compression';
  if (w < 180) return 'power';
  return 'exhaust';
}

/** Start angle (inclusive) of a stroke in the cycle range. */
export function strokeStartDeg(s: StrokeName): number {
  switch (s) {
    case 'intake':
      return -360;
    case 'compression':
      return -180;
    case 'power':
      return 0;
    case 'exhaust':
      return 180;
  }
}

/** Fraction (0..1) of the current stroke that has elapsed at θ. */
export function strokeProgress(theta: number): number {
  const w = wrapCycleDeg(theta);
  return (w - strokeStartDeg(strokeOf(w))) / 180;
}

/** Mechanical crank angle within one revolution, [0, 360), 0 = TDC. */
export function crankAngleInRevolution(theta: number): number {
  return ((theta % 360) + 360) % 360;
}

/**
 * Angle on a 720°-per-turn "cycle dial" in degrees clockwise from 12 o'clock:
 * firing TDC at the top, compression in the upper-left quadrant, power upper-right,
 * exhaust lower-right, intake lower-left.
 */
export function cycleDialAngleDeg(theta: number): number {
  const d = wrapCycleDeg(theta) / 2; // [-180, 180)
  return (d + 360) % 360;
}

/** Position on a dial of radius r centred at (cx, cy) for crank angle θ (SVG coords, y down). */
export function cycleDialPoint(theta: number, cx: number, cy: number, r: number): [number, number] {
  const a = (cycleDialAngleDeg(theta) * Math.PI) / 180;
  return [cx + r * Math.sin(a), cy - r * Math.cos(a)];
}

/**
 * Human description of θ relative to the nearest TDC, e.g. "13.0° BTDC", "5.5° ATDC",
 * "20.0° ABDC". Which TDC (firing / gas-exchange) is implied by the stroke.
 */
export function describeCrankAngle(theta: number, decimals = 1): string {
  const w = wrapCycleDeg(theta);
  const inRev = crankAngleInRevolution(w); // 0 = a TDC, 180 = a BDC
  const toTdc = inRev <= 180 ? inRev : inRev - 360; // (-180, 180]
  if (Math.abs(toTdc) <= 90) {
    const a = Math.abs(toTdc).toFixed(decimals);
    return Number(a) === 0 ? 'TDC' : toTdc < 0 ? `${a}° BTDC` : `${a}° ATDC`;
  }
  const toBdc = inRev - 180; // (-90, 90) here
  const a = Math.abs(toBdc).toFixed(decimals);
  return Number(a) === 0 ? 'BDC' : toBdc < 0 ? `${a}° BBDC` : `${a}° ABDC`;
}

/** Which TDC is nearest: 'firing', 'gas-exchange', or null (closer to a BDC). */
export function nearestTdc(theta: number): 'firing' | 'gas-exchange' | null {
  const w = wrapCycleDeg(theta);
  if (Math.abs(w) <= 90) return 'firing';
  if (w <= -270 || w >= 270) return 'gas-exchange';
  return null;
}

/**
 * Open interval of a valve in the cycle range as [open, close] with close > open
 * (close may exceed 360 when the event straddles the cycle boundary, e.g. an exhaust
 * valve opening at 140° and closing 15° after gas-exchange TDC → [140, 375]).
 */
export function valveOpenInterval(v: Pick<ValveSpec, 'openDeg' | 'closeDeg'>): [number, number] {
  const o = wrapCycleDeg(v.openDeg);
  let c = wrapCycleDeg(v.closeDeg);
  if (c <= o) c += 720;
  return [o, c];
}

/** Whether a valve with the given timing is open (per its timing threshold) at θ. */
export function isValveOpen(v: Pick<ValveSpec, 'openDeg' | 'closeDeg'>, theta: number): boolean {
  const [o, c] = valveOpenInterval(v);
  const w = wrapCycleDeg(theta);
  return (w >= o && w < c) || (w + 720 >= o && w + 720 < c);
}

export interface ValveEvents {
  /** Intake valve opens / closes, cycle degrees (wrapped). */
  ivo: number;
  ivc: number;
  /** Exhaust valve opens / closes, cycle degrees (wrapped). */
  evo: number;
  evc: number;
}

export function valveEvents(spec: EngineSpec): ValveEvents {
  return {
    ivo: wrapCycleDeg(spec.intakeValve.openDeg),
    ivc: wrapCycleDeg(spec.intakeValve.closeDeg),
    evo: wrapCycleDeg(spec.exhaustValve.openDeg),
    evc: wrapCycleDeg(spec.exhaustValve.closeDeg),
  };
}

/** Clearance volume and swept volume of the cylinder from geometry, m³. */
export function cylinderVolumes(spec: EngineSpec, compressionRatio: number): { vd: number; vc: number; vmax: number } {
  const b = spec.geometry.bore;
  const vd = (Math.PI / 4) * b * b * spec.geometry.stroke;
  const vc = vd / Math.max(compressionRatio - 1, 1e-6);
  return { vd, vc, vmax: vd + vc };
}
