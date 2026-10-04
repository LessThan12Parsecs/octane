/**
 * Carburettor with a fixed venturi and a butterfly throttle in series (ManifoldSpec.venturiDiameter
 * present; the Model T Holley/Kingston class). Without venturiDiameter the cycle model keeps its
 * historical venturi-as-throttle orifice (the CFR: no throttle plate; ManifoldSpec.throttleDiameter
 * is the 9/16 in venturi), which this module does not touch.
 *
 * ── Restrictions ────────────────────────────────────────────────────────────────────────
 *  - Venturi: C_D,v · πD_v²/4 as a full-loss orifice (the same convention as the CFR venturi; C_D,v is
 *    the calibrated CycleModelOptions.venturiDischargeCoefficient or an engine calibration). An
 *    optional fixed series restriction (CycleModelOptions.intakeRestrictionArea: air horn, bends) is
 *    folded into it incompressibly, 1/(C_D A)² = 1/(C_D A)²_v + 1/(C_D A)²_r, as today.
 *  - Butterfly: geometric open area A_o(u) from throttle.ts throttleArea (Harrington–Bolt + leakage)
 *    with the plate data of ManifoldSpec.throttle. Its effective full-loss area follows from the
 *    vena contracta and the re-expansion into the bore downstream:
 *      A_c = C_c A_o,   Δp_loss = ½ρ v_c² (1 − A_c/A_bore)²   (Borda–Carnot sudden expansion)
 *      ⇒ (C_D A)_t = A_c / (1 − A_c/A_bore).
 *    C_c is the free-streamline contraction of the sharp-edged plate gaps, π/(π + 2) = 0.611 for a
 *    2-D slot (Kirchhoff 1869 / von Mises 1917; e.g. Batchelor 1967 §6.13 — UNVERIFIED section; it
 *    neglects the approach velocity, so it is exact only for small openings, where the throttle
 *    controls and chokes the flow). At small openings (C_D A)_t → 0.611 A_o; wide open (Model T plate:
 *    A_o/A_bore ≈ 0.75) the re-expansion recovers most of the head, (C_D A)_t ≈ 1.12 A_o.
 *    UNVERIFIED: no flow-bench C_D of a Model T butterfly was found.
 *  - The two are solved as compressible restrictions in series (series-orifice.ts: Newton on the
 *    pressure between venturi and plate), which stays exact when the throttle chokes at idle and
 *    returns ∂ṁ/∂p_manifold for the cycle model's stiffness bound.
 *
 * Hot paths allocate nothing; setOpening() is cheap (one throttleArea) and is meant to be called when
 * the throttle opening changes.
 */
import { DEG } from '../core/constants';
import type { ManifoldSpec, ThrottlePlateSpec } from '../core/engine-spec';
import { seriesOrificeFlow, type SeriesOrificeFlow } from './series-orifice';
import { throttleArea, type ThrottleOptions } from './throttle';
import { seriesEffectiveArea, SLOT_CONTRACTION_COEFFICIENT } from './valve-flow';

/**
 * Effective full-loss area of a butterfly throttle (Kirchhoff contraction + Borda–Carnot
 * re-expansion; file header), m².
 * @param diameter throttle bore, m
 * @param opening u ∈ [0, 1]
 * @param opts plate options (throttle.ts)
 * @param contraction jet contraction coefficient C_c (default π/(π + 2))
 */
export function throttleEffectiveArea(
  diameter: number,
  opening: number,
  opts: ThrottleOptions = {},
  contraction: number = SLOT_CONTRACTION_COEFFICIENT,
): number {
  const Ab = 0.25 * Math.PI * diameter * diameter;
  const Ac = contraction * throttleArea(diameter, opening, opts);
  if (!(Ac > 0)) return 0;
  const s = Ac / Ab;
  return s < 1 ? Ac / (1 - s) : Infinity;
}

/** ThrottleOptions of a spec plate (deg → rad). */
export function throttlePlateOptions(plate: ThrottlePlateSpec | undefined): ThrottleOptions {
  if (plate === undefined) return {};
  return { shaftRatio: plate.shaftRatio, closedAngle: plate.closedAngleDeg * DEG, leakageArea: plate.leakageArea };
}

/** True when the manifold spec describes a venturi with a separate butterfly (series model applies). */
export function hasSeparateThrottle(manifolds: ManifoldSpec): boolean {
  return manifolds.venturiDiameter !== undefined && manifolds.venturiDiameter > 0;
}

/** Options of {@link CarburettorFlowModel}. */
export interface CarburettorOptions {
  /** Venturi throat discharge coefficient (full-loss orifice convention). */
  venturiDischargeCoefficient: number;
  /** Fixed series restriction C_D·A folded into the venturi, m² (0/absent: none). */
  restrictionArea?: number;
  /** Throttle-gap jet contraction coefficient (default π/(π + 2)). */
  throttleContraction?: number;
}

/**
 * Venturi + butterfly carburettor between the ambient (side a) and the intake manifold (side b).
 * Allocation-free after construction.
 */
export class CarburettorFlowModel {
  /** Venturi effective area (incl. the folded series restriction), m². */
  readonly venturiCdA: number;
  /** Throttle bore, m. */
  readonly throttleDiameter: number;
  /** Throttle bore area, m². */
  readonly throttleBoreArea: number;
  /** Current throttle effective area, m². */
  throttleCdA = 0;
  /** Opening the throttle area was computed for (NaN = none). */
  opening = Number.NaN;
  private readonly plate: ThrottleOptions;
  private readonly cc: number;

  constructor(manifolds: ManifoldSpec, opts: CarburettorOptions) {
    const Dv = manifolds.venturiDiameter;
    if (Dv === undefined || !(Dv > 0)) throw new RangeError('CarburettorFlowModel: ManifoldSpec.venturiDiameter required');
    const Dt = manifolds.throttleDiameter;
    if (!(Dt > 0)) throw new RangeError('CarburettorFlowModel: throttleDiameter must be > 0');
    if (!(opts.venturiDischargeCoefficient > 0)) throw new RangeError('CarburettorFlowModel: venturi C_D must be > 0');
    const cdaV = opts.venturiDischargeCoefficient * 0.25 * Math.PI * Dv * Dv;
    const r = opts.restrictionArea ?? 0;
    this.venturiCdA = r > 0 ? seriesEffectiveArea(cdaV, r) : cdaV;
    this.throttleDiameter = Dt;
    this.throttleBoreArea = 0.25 * Math.PI * Dt * Dt;
    this.plate = throttlePlateOptions(manifolds.throttle);
    this.cc = opts.throttleContraction ?? SLOT_CONTRACTION_COEFFICIENT;
    if (!(this.cc > 0 && this.cc <= 1)) throw new RangeError('CarburettorFlowModel: throttleContraction must be in (0, 1]');
  }

  /** Set the throttle opening u ∈ [0, 1] (recomputes the throttle effective area when it changed). */
  setOpening(u: number): void {
    if (u === this.opening) return;
    this.throttleCdA = throttleEffectiveArea(this.throttleDiameter, u, this.plate, this.cc);
    this.opening = u;
  }

  /**
   * Mass flow ambient (a) → manifold (b) through venturi then throttle (signed; reverse flow passes
   * the throttle first). Writes `out` (drop-in OrificeFlow: mdot, dmdotdpa, dmdotdpb, …) and returns
   * out.mdot, kg/s.
   */
  flow(
    pa: number,
    Ta: number,
    Ra: number,
    gammaA: number,
    pb: number,
    Tb: number,
    Rb: number,
    gammaB: number,
    out: SeriesOrificeFlow,
  ): number {
    return seriesOrificeFlow(this.venturiCdA, this.throttleCdA, pa, Ta, Ra, gammaA, pb, Tb, Rb, gammaB, out);
  }
}
