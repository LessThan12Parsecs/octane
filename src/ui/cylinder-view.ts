/**
 * Multi-cylinder helpers for the UI (pure, no DOM).
 *
 * A multi-cylinder snapshot carries every cylinder in `s.cylinders[i]` at its LOCAL crank angle
 * (snapshot.ts: CylinderSnapshot); the top level is the engine state plus cylinder 1. The UI binds its
 * single-cylinder consumers (HUD readouts, trace store, spark capture) to one "focus" cylinder by
 * handing them a VIEW: an EngineSnapshot-shaped object whose per-cylinder fields come from
 * `s.cylinders[i]` (θ = local angle, cycle = that cylinder's cycle) and whose engine-level fields come
 * from `s`. Nested objects are shared by reference, so building a view allocates nothing.
 */
import { cylinderAngleDeg, type EngineSpec } from '../physics/core/engine-spec';
import type { CylinderSnapshot, EngineSnapshot } from '../physics/core/snapshot';
import { strokeOf, type StrokeName } from './engine-cycle';
import { isDischargePhase } from './spark-capture';

/** Number of cylinders whose data the UI shows (1 for a single-cylinder spec). */
export function cylinderCount(spec: EngineSpec): number {
  return Math.max(1, Math.floor(spec.cylinders));
}

/** Clamp a (possibly stale or invalid) cylinder index into [0, n). */
export function clampCylinder(i: number, n: number): number {
  if (!Number.isFinite(i)) return 0;
  return Math.min(Math.max(0, Math.floor(i)), Math.max(0, n - 1));
}

/** An empty view target (fields are overwritten by cylinderView). */
export function createCylinderViewTarget(): EngineSnapshot {
  return {} as EngineSnapshot;
}

/**
 * Cylinder `i` (0-based) of `s` as an EngineSnapshot written into `out` (no allocation). Returns `s`
 * itself when the snapshot has no `cylinders[]` (single-cylinder stream) or `i` is out of range.
 */
export function cylinderView(s: EngineSnapshot, i: number, out: EngineSnapshot): EngineSnapshot {
  const cyl = s.cylinders;
  if (!cyl || i < 0 || i >= cyl.length) return s;
  const c: CylinderSnapshot = cyl[i];
  // engine level
  out.t = s.t;
  out.rpm = s.rpm;
  out.intakeManifoldPressure = s.intakeManifoldPressure;
  out.exhaustManifoldPressure = s.exhaustManifoldPressure;
  out.netTorque = s.netTorque;
  out.frictionTorque = s.frictionTorque;
  out.loadTorque = s.loadTorque;
  out.vehicleSpeed = s.vehicleSpeed;
  out.firingCylinder = s.firingCylinder;
  out.magnetoEmf = s.magnetoEmf;
  out.cylinders = cyl;
  // this cylinder
  out.cycle = c.cycle;
  out.thetaDeg = c.thetaDeg;
  out.gasTorque = c.gasTorque;
  out.pistonDisplacement = c.pistonDisplacement;
  out.clearanceHeight = c.clearanceHeight;
  out.rodAngle = c.rodAngle;
  out.intakeLift = c.intakeLift;
  out.exhaustLift = c.exhaustLift;
  out.phase = c.phase;
  out.volume = c.volume;
  out.pressure = c.pressure;
  out.temperatureMean = c.temperatureMean;
  out.temperatureUnburned = c.temperatureUnburned;
  out.temperatureBurned = c.temperatureBurned;
  out.massFractionBurned = c.massFractionBurned;
  out.mass = c.mass;
  out.heatReleaseRate = c.heatReleaseRate;
  out.heatLossRate = c.heatLossRate;
  out.flame = c.flame;
  out.spark = c.spark;
  out.intakeMassFlow = c.intakeMassFlow;
  out.exhaustMassFlow = c.exhaustMassFlow;
  out.knock = c.knock;
  out.burnedComposition = c.burnedComposition;
  return out;
}

/**
 * Pool of view objects for mapping an array of raw snapshots to cylinder `i` views. The returned
 * array and its views are reused by the next map() call with the same pool.
 */
export class CylinderViewPool {
  private readonly views: EngineSnapshot[] = [];
  private readonly out: EngineSnapshot[] = [];

  map(snaps: readonly EngineSnapshot[], i: number): EngineSnapshot[] {
    const out = this.out;
    out.length = snaps.length;
    for (let k = 0; k < snaps.length; k++) {
      let v = this.views[k];
      if (!v) v = this.views[k] = createCylinderViewTarget();
      out[k] = cylinderView(snaps[k], i, v);
    }
    return out;
  }
}

/** One cylinder's entry of the firing-order strip. */
export interface FiringStripEntry {
  /** 0-based cylinder index. */
  index: number;
  /** Cylinder number (1-based, as stamped on the engine). */
  number: number;
  /** Local crank angle, deg. */
  thetaDeg: number;
  stroke: StrokeName;
  /** The gap is conducting (breakdown, arc or glow). */
  sparking: boolean;
  /** The ignition timer is grounding this cylinder's coil (trembler ignition). */
  timerClosed: boolean;
}

/** Firing order (1-based cylinder numbers) of a spec: layout.firingOrder, or [1] for one cylinder. */
export function firingOrder(spec: EngineSpec): number[] {
  const n = cylinderCount(spec);
  const fo = spec.layout?.firingOrder;
  if (fo && fo.length === n) return fo.slice();
  return Array.from({ length: n }, (_, i) => i + 1);
}

/**
 * The firing-order strip: every cylinder in firing order with its stroke and spark state. Local angles
 * come from `s.cylinders` when present, else from the engine angle and the spec's firing offsets.
 */
export function firingStrip(spec: EngineSpec, s: EngineSnapshot, out: FiringStripEntry[] = []): FiringStripEntry[] {
  const order = firingOrder(spec);
  out.length = order.length;
  for (let k = 0; k < order.length; k++) {
    const i = order[k] - 1;
    const c = s.cylinders?.[i];
    const th = c ? c.thetaDeg : cylinderAngleDeg(spec, i, s.thetaDeg);
    const sp = c ? c.spark : i === 0 ? s.spark : null;
    const e = out[k] ?? (out[k] = { index: i, number: i + 1, thetaDeg: 0, stroke: 'intake', sparking: false, timerClosed: false });
    e.index = i;
    e.number = i + 1;
    e.thetaDeg = th;
    e.stroke = strokeOf(th);
    e.sparking = sp ? isDischargePhase(sp.phase) : false;
    e.timerClosed = (sp?.timerClosed ?? false) || s.firingCylinder === i;
  }
  return out;
}
