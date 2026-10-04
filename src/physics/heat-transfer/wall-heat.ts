/**
 * Split of the gas-to-wall heat loss over the chamber surfaces (flat head with valve faces,
 * flat piston crown, exposed liner) and, when flame-geometry wetted areas are supplied, over
 * the burned and unburned zones.
 *
 * Per surface i with area A_i at temperature T_w,i (WallSpec):
 *   Q̇_i = h A_i (T_g − T_w,i) + c A_i (T_g⁴ − T_w,i⁴)       [W, positive gas → wall]
 * c is Annand's radiation constant (0 for Woschni/Hohenberg, which lump radiation into h).
 *
 * Two-zone split (DESIGN.md §Heat transfer): the same h (evaluated by the caller with the
 * mass-averaged temperature, as the correlations were fitted) is applied to each zone's own
 * temperature and wetted area — a common two-zone practice, e.g. Heywood (1988) §14.4
 * (UNVERIFIED section). The burned-wetted head area from the flame geometry covers the whole
 * head disc (valve faces included); its fraction f_h = wettedHead / (πB²/4) is applied
 * uniformly to the head deck and both valve faces (assumption: the valves lie anywhere under
 * the flame footprint with equal probability).
 *
 * All inputs SI; allocation-free given a reused result object.
 */
import type { ValveSpec, WallSpec, WallThermalResistance } from '../core/engine-spec';

/** Gas-exposed chamber areas, m². */
export interface ChamberAreas {
  /** Head fire deck EXCLUDING the valve faces. */
  head: number;
  /** Intake valve faces (all intake valves). */
  intakeValves: number;
  /** Exhaust valve faces (all exhaust valves). */
  exhaustValves: number;
  /** Piston crown. */
  piston: number;
  /** Exposed liner π B h(θ). */
  liner: number;
}

/** Allocate a zeroed {@link ChamberAreas}. */
export function newChamberAreas(): ChamberAreas {
  return { head: 0, intakeValves: 0, exhaustValves: 0, piston: 0, liner: 0 };
}

/**
 * Flat-disc ("pancake") chamber areas at clearance height h (m): head disc πB²/4 split into
 * the valve faces (count · πD_v²/4) and the remaining deck, piston crown πB²/4, liner πBh.
 * Writes into and returns `out`.
 */
export function flatChamberAreas(
  bore: number,
  clearanceHeight: number,
  intake: Pick<ValveSpec, 'count' | 'headDiameter'>,
  exhaust: Pick<ValveSpec, 'count' | 'headDiameter'>,
  out: ChamberAreas,
): ChamberAreas {
  const disc = 0.25 * Math.PI * bore * bore;
  out.intakeValves = intake.count * 0.25 * Math.PI * intake.headDiameter * intake.headDiameter;
  out.exhaustValves = exhaust.count * 0.25 * Math.PI * exhaust.headDiameter * exhaust.headDiameter;
  out.head = Math.max(0, disc - out.intakeValves - out.exhaustValves);
  out.piston = disc;
  out.liner = Math.PI * bore * (clearanceHeight > 0 ? clearanceHeight : 0);
  return out;
}

/** Burned-zone wetted areas from the flame geometry, m² (combustion/flame-geometry.ts). */
export interface BurnedWettedAreas {
  /** Head disc (incl. valve faces) inside the flame ball. */
  wettedHead: number;
  /** Piston crown inside the flame ball. */
  wettedPiston: number;
  /** Liner inside the flame ball. */
  wettedLiner: number;
}

/** Heat-loss rates, W (positive gas → wall). */
export interface WallHeatResult {
  total: number;
  head: number;
  intakeValves: number;
  exhaustValves: number;
  piston: number;
  liner: number;
  /** Part of `total` leaving the burned zone (0 in single-zone mode). */
  burned: number;
  /** Part of `total` leaving the unburned zone (= total in single-zone mode). */
  unburned: number;
}

/** Allocate a zeroed {@link WallHeatResult}. */
export function newWallHeatResult(): WallHeatResult {
  return { total: 0, head: 0, intakeValves: 0, exhaustValves: 0, piston: 0, liner: 0, burned: 0, unburned: 0 };
}

const q4 = (T: number, Tw: number): number => {
  const a = T * T;
  const b = Tw * Tw;
  return a * a - b * b;
};

/** Flux to one surface, W/m². */
function flux(h: number, c: number, T: number, Tw: number): number {
  return h * (T - Tw) + (c !== 0 ? c * q4(T, Tw) : 0);
}

/**
 * Single-zone wall heat loss, W (positive gas → wall), split over the surfaces.
 * @param h heat-transfer coefficient, W/(m² K)
 * @param T gas temperature, K
 * @param areas chamber areas, m²
 * @param walls surface temperatures, K
 * @param c Annand radiation constant, W/(m² K⁴) (default 0)
 */
export function wallHeatLoss(h: number, T: number, areas: ChamberAreas, walls: WallSpec, out: WallHeatResult, c = 0): number {
  out.head = areas.head * flux(h, c, T, walls.headTemperature);
  out.intakeValves = areas.intakeValves * flux(h, c, T, walls.intakeValveTemperature);
  out.exhaustValves = areas.exhaustValves * flux(h, c, T, walls.exhaustValveTemperature);
  out.piston = areas.piston * flux(h, c, T, walls.pistonTemperature);
  out.liner = areas.liner * flux(h, c, T, walls.linerTemperature);
  out.total = out.head + out.intakeValves + out.exhaustValves + out.piston + out.liner;
  out.burned = 0;
  out.unburned = out.total;
  return out.total;
}

/**
 * Two-zone wall heat loss, W (positive gas → wall): burned zone at T_b over the wetted areas,
 * unburned zone at T_u over the rest; `out.burned` / `out.unburned` give the zone split.
 * Wetted areas are clamped to the surface areas.
 */
export function wallHeatLossTwoZone(
  h: number,
  Tu: number,
  Tb: number,
  areas: ChamberAreas,
  wetted: BurnedWettedAreas,
  walls: WallSpec,
  out: WallHeatResult,
  c = 0,
): number {
  const headDisc = areas.head + areas.intakeValves + areas.exhaustValves;
  let fh = headDisc > 0 ? wetted.wettedHead / headDisc : 0;
  fh = fh > 0 ? (fh < 1 ? fh : 1) : 0;
  let fp = areas.piston > 0 ? wetted.wettedPiston / areas.piston : 0;
  fp = fp > 0 ? (fp < 1 ? fp : 1) : 0;
  let fl = areas.liner > 0 ? wetted.wettedLiner / areas.liner : 0;
  fl = fl > 0 ? (fl < 1 ? fl : 1) : 0;

  let qb = 0;
  let qu = 0;
  let b: number;
  let u: number;
  // head deck
  b = fh * areas.head * flux(h, c, Tb, walls.headTemperature);
  u = (1 - fh) * areas.head * flux(h, c, Tu, walls.headTemperature);
  out.head = b + u;
  qb += b;
  qu += u;
  b = fh * areas.intakeValves * flux(h, c, Tb, walls.intakeValveTemperature);
  u = (1 - fh) * areas.intakeValves * flux(h, c, Tu, walls.intakeValveTemperature);
  out.intakeValves = b + u;
  qb += b;
  qu += u;
  b = fh * areas.exhaustValves * flux(h, c, Tb, walls.exhaustValveTemperature);
  u = (1 - fh) * areas.exhaustValves * flux(h, c, Tu, walls.exhaustValveTemperature);
  out.exhaustValves = b + u;
  qb += b;
  qu += u;
  b = fp * areas.piston * flux(h, c, Tb, walls.pistonTemperature);
  u = (1 - fp) * areas.piston * flux(h, c, Tu, walls.pistonTemperature);
  out.piston = b + u;
  qb += b;
  qu += u;
  b = fl * areas.liner * flux(h, c, Tb, walls.linerTemperature);
  u = (1 - fl) * areas.liner * flux(h, c, Tu, walls.linerTemperature);
  out.liner = b + u;
  qb += b;
  qu += u;
  out.burned = qb;
  out.unburned = qu;
  out.total = qb + qu;
  return out.total;
}

// ---------------------------------------------------------------------------------------------
// Chamber-generic wall surfaces (combustion/chamber.ts: surfaceAreas / evaluate().burnedFraction)
//
// The same per-surface law over N_WALL_SURFACES = 6 surfaces held in Float64Arrays: the five of the
// flat-disc split plus the BLOCK (side-valve 'l-head' chambers: the water-jacketed block deck / valve-
// pocket floor between the bore and the valve seats, WallSpec.blockTemperature). Each surface carries
// its OWN burned fraction (an L-head valve face is wetted when the flame covers the pocket floor, not
// in proportion to the head roof). With the flat-disc areas (block 0) and the disc fractions of
// DiscChamber (valves at the head fraction) the sums run in the order of wallHeatLoss /
// wallHeatLossTwoZone (head, intake, exhaust, piston, liner, block), so the results are bit-identical
// to the 5-surface functions (wall-heat-surfaces.test.ts, combustion/chamber.test.ts).
// ---------------------------------------------------------------------------------------------

/** Surface index: head (flat head fire deck / L-head cavity roof and head-casting walls), valve faces excluded. */
export const WALL_HEAD = 0;
/** Surface index: piston (crown, plus the side above the deck when an L-head piston rises into the head). */
export const WALL_PISTON = 1;
/** Surface index: liner (bore wall between the crown and the head face / block deck). */
export const WALL_LINER = 2;
/** Surface index: intake valve faces (all intake valves). */
export const WALL_INTAKE_VALVE = 3;
/** Surface index: exhaust valve faces (all exhaust valves). */
export const WALL_EXHAUST_VALVE = 4;
/** Surface index: block deck / valve-pocket floor outside the valve faces ('l-head'; 0 area for a flat disc). */
export const WALL_BLOCK = 5;
/** Number of chamber wall surfaces of the generic path. */
export const N_WALL_SURFACES = 6;

/** Surface ids in index order (= the cycle model's per-surface heat ledgers, block appended last). */
export type WallSurfaceId = 'head' | 'piston' | 'liner' | 'intakeValve' | 'exhaustValve' | 'block';
export const WALL_SURFACE_IDS: readonly WallSurfaceId[] = Object.freeze([
  'head',
  'piston',
  'liner',
  'intakeValve',
  'exhaustValve',
  'block',
] as WallSurfaceId[]);

/** Summation order reproducing wallHeatLoss / wallHeatLossTwoZone bit for bit. */
const SUM_ORDER = [WALL_HEAD, WALL_INTAKE_VALVE, WALL_EXHAUST_VALVE, WALL_PISTON, WALL_LINER, WALL_BLOCK];

/** Allocate a zeroed per-surface array (areas, temperatures, fractions, …). */
export function newWallSurfaceArray(): Float64Array {
  return new Float64Array(N_WALL_SURFACES);
}

/** Surface temperatures, K, in WALL_* order; the block falls back to the liner temperature. Returns `out`. */
export function wallSurfaceTemperatures(walls: WallSpec, out: Float64Array): Float64Array {
  out[WALL_HEAD] = walls.headTemperature;
  out[WALL_PISTON] = walls.pistonTemperature;
  out[WALL_LINER] = walls.linerTemperature;
  out[WALL_INTAKE_VALVE] = walls.intakeValveTemperature;
  out[WALL_EXHAUST_VALVE] = walls.exhaustValveTemperature;
  out[WALL_BLOCK] = walls.blockTemperature ?? walls.linerTemperature;
  return out;
}

/**
 * Lumped-wall thermal resistances, K/W, in WALL_* order (WallSpec.thermalResistance); the block falls back
 * to the liner resistance (both water-jacketed block surfaces). Returns `out`.
 */
export function wallSurfaceResistances(R: WallThermalResistance, out: Float64Array): Float64Array {
  out[WALL_HEAD] = R.head;
  out[WALL_PISTON] = R.piston;
  out[WALL_LINER] = R.liner;
  out[WALL_INTAKE_VALVE] = R.intakeValve;
  out[WALL_EXHAUST_VALVE] = R.exhaustValve;
  out[WALL_BLOCK] = R.block ?? R.liner;
  return out;
}

/** {@link ChamberAreas} → per-surface areas (block 0). Returns `out`. */
export function chamberAreasToSurfaces(areas: ChamberAreas, out: Float64Array): Float64Array {
  out[WALL_HEAD] = areas.head;
  out[WALL_PISTON] = areas.piston;
  out[WALL_LINER] = areas.liner;
  out[WALL_INTAKE_VALVE] = areas.intakeValves;
  out[WALL_EXHAUST_VALVE] = areas.exhaustValves;
  out[WALL_BLOCK] = 0;
  return out;
}

/** Heat-loss rates of the generic path, W (positive gas → wall). */
export interface WallSurfaceHeatResult {
  total: number;
  /** Part of `total` leaving the burned zone (0 in single-zone mode). */
  burned: number;
  /** Part of `total` leaving the unburned zone (= total in single-zone mode). */
  unburned: number;
  /** Per-surface heat flow, W, WALL_* order. */
  surface: Float64Array;
}

/** Allocate a zeroed {@link WallSurfaceHeatResult}. */
export function newWallSurfaceHeatResult(): WallSurfaceHeatResult {
  return { total: 0, burned: 0, unburned: 0, surface: new Float64Array(N_WALL_SURFACES) };
}

/**
 * Single-zone wall heat loss over the generic surfaces, W (positive gas → wall).
 * @param h heat-transfer coefficient, W/(m² K)
 * @param T gas temperature, K
 * @param areas per-surface areas, m² (WALL_* order)
 * @param Tw per-surface temperatures, K ({@link wallSurfaceTemperatures})
 * @param c Annand radiation constant, W/(m² K⁴) (default 0)
 */
export function wallHeatLossSurfaces(
  h: number,
  T: number,
  areas: Float64Array,
  Tw: Float64Array,
  out: WallSurfaceHeatResult,
  c = 0,
): number {
  const q = out.surface;
  let total = 0;
  for (let k = 0; k < N_WALL_SURFACES; k++) {
    const i = SUM_ORDER[k];
    const qi = areas[i] * flux(h, c, T, Tw[i]);
    q[i] = qi;
    total = k === 0 ? qi : total + qi;
  }
  out.total = total;
  out.burned = 0;
  out.unburned = total;
  return total;
}

/**
 * Two-zone wall heat loss over the generic surfaces, W (positive gas → wall): on surface i the burned zone
 * (T_b) covers the fraction f_i of the area, the unburned zone (T_u) the rest. Fractions are clamped to
 * [0, 1] (NaN → 0). `out.burned` / `out.unburned` give the zone split.
 * @param burnedFraction per-surface burned fractions (chamber evaluate().burnedFraction)
 */
export function wallHeatLossTwoZoneSurfaces(
  h: number,
  Tu: number,
  Tb: number,
  areas: Float64Array,
  burnedFraction: Float64Array,
  Tw: Float64Array,
  out: WallSurfaceHeatResult,
  c = 0,
): number {
  const q = out.surface;
  let qb = 0;
  let qu = 0;
  for (let k = 0; k < N_WALL_SURFACES; k++) {
    const i = SUM_ORDER[k];
    let f = burnedFraction[i];
    f = f > 0 ? (f < 1 ? f : 1) : 0;
    const A = areas[i];
    const b = f * A * flux(h, c, Tb, Tw[i]);
    const u = (1 - f) * A * flux(h, c, Tu, Tw[i]);
    q[i] = b + u;
    qb += b;
    qu += u;
  }
  out.burned = qb;
  out.unburned = qu;
  out.total = qb + qu;
  return out.total;
}
