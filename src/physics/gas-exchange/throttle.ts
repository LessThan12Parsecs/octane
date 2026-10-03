/**
 * Butterfly throttle open area and carburettor-venturi loss.
 *
 * ── Throttle geometry ───────────────────────────────────────────────────────────────────
 * Thin elliptical plate that exactly fills a bore of diameter D at the closed angle ψ₀ (angle
 * between the plate and the plane normal to the bore axis), rotating to ψ about a shaft of
 * diameter d = a·D across the bore. Projected on the bore cross-section the plate is an ellipse
 * of semi-axes D/2 (along the shaft) and (D/2)·k, k = cosψ/cosψ₀, and the shaft a band |y| ≤ aD/2.
 * Open area = bore − (plate ∪ shaft), which integrates in closed form to
 *
 *   A/(πD²/4) = 1 − cosψ/cosψ₀
 *               + (2/π)[ (a/cosψ)√(cos²ψ − a² cos²ψ₀) + (cosψ/cosψ₀) asin(a cosψ₀/cosψ)
 *                        − a√(1 − a²) − asin a ]                     for cosψ ≥ a cosψ₀,
 *   A/(πD²/4) = 1 − (2/π)[ a√(1 − a²) + asin a ]                     (plate hidden by the shaft).
 *
 * This is the Harrington & Bolt (1970) throttle-area formula used in Heywood (1988) §7 and in
 * Guzzella & Onder (2010); its structure — including the switch at ψ = acos(a cosψ₀) and the
 * reference "Harrington, D., Bolt, J., Analysis and digital simulation of carburetor metering,
 * SAE 700082" — was checked against IJCRT (2018) paper IJCRT1892376 (fetched; the task brief's
 * "SAE 700087" appears to be a typo). The expression itself is DERIVED here from the geometry
 * above and verified by direct numerical integration (tools/reference/gasex_throttle.py), so it
 * does not rely on the transcription. UNVERIFIED: the Heywood equation number (7.x).
 * Neglected: plate thickness and edge rounding, the shaft's own bosses.
 *
 * Opening u ∈ [0, 1] maps linearly to the plate angle ψ = ψ₀ + u (90° − ψ₀). A constant
 * closed-plate leakage area (clearance gap, shaft bushings) is added.
 * The throttle discharge coefficient (C_D ≈ 0.6–0.9, angle and pressure-ratio dependent) must be
 * applied by the caller from flow-bench data; none is assumed here.
 *
 * ── Venturi ─────────────────────────────────────────────────────────────────────────────
 * The CFR has no throttle plate (throttle = 1 always): the only intake restriction is the
 * carburettor venturi (ManifoldSpec.throttleDiameter holds the throat, 9/16 in for the CFR).
 * Throat depression Δp_t = ṁ²/(2ρ (C_D A_t)²); most of it is recovered in the diffuser. The
 * unrecovered stagnation-pressure loss is taken as a fraction of Δp_t (classical venturi tubes:
 * 5–20 %, UNVERIFIED — from memory of ISO 5167-4). CFR estimate at 600 rpm: mean ṁ ≈ 2.7 g/s,
 * A_t = 1.6 cm² → v_t ≈ 15 m/s, Δp_t ≈ 0.13 kPa (≈ 1 kPa at the pulsation peaks), unrecovered
 * ≲ 0.05 kPa (≈ 0.05 % of p_atm): negligible for the trapped mass. Model the venturi as a
 * fixed-area orifice (effectiveArea = C_D·A_t) when the intake plenum is simulated.
 */

/** Throttle options. */
export interface ThrottleOptions {
  /** Shaft diameter / bore diameter a (0 ≤ a < 1). */
  shaftRatio?: number;
  /** Closed-plate angle ψ₀, rad (0 ≤ ψ₀ < π/2). */
  closedAngle?: number;
  /** Constant leakage area added at all openings, m². */
  leakageArea?: number;
}

/**
 * Default throttle parameters.
 * UNVERIFIED: typical automotive values — shaft ≈ 12 % of the bore, closed angle 8°,
 * leakage = annular gap of 25 µm radial clearance on a 50 mm bore (≈ 3.9 mm²; scaled with D
 * by throttleArea when leakageArea is not given).
 */
export const THROTTLE_DEFAULTS = Object.freeze({
  shaftRatio: 0.12,
  closedAngle: (8 * Math.PI) / 180,
  /** Radial clearance used for the default leakage π D c, m. */
  leakageClearance: 25e-6,
});

/** Plate angle ψ (rad) for opening u ∈ [0, 1] (clamped), linear between ψ₀ and π/2. */
export function throttlePlateAngle(opening: number, closedAngle: number = THROTTLE_DEFAULTS.closedAngle): number {
  const u = opening > 0 ? (opening < 1 ? opening : 1) : 0;
  return closedAngle + u * (0.5 * Math.PI - closedAngle);
}

/**
 * Geometric (projected) open area of a butterfly throttle, m² (no leakage).
 * @param diameter bore D, m
 * @param plateAngle ψ, rad (ψ ≤ ψ₀ → 0; ψ ≥ π/2 → fully open)
 * @param shaftRatio a = d/D
 * @param closedAngle ψ₀, rad
 */
export function throttleGeometricArea(
  diameter: number,
  plateAngle: number,
  shaftRatio: number = THROTTLE_DEFAULTS.shaftRatio,
  closedAngle: number = THROTTLE_DEFAULTS.closedAngle,
): number {
  if (!(diameter > 0) || !(plateAngle > closedAngle)) return 0;
  const a = shaftRatio > 0 ? (shaftRatio < 1 ? shaftRatio : 1) : 0;
  const A0 = 0.25 * Math.PI * diameter * diameter;
  const c0 = Math.cos(closedAngle);
  const shaft = a * Math.sqrt(1 - a * a) + Math.asin(a);
  const psi = plateAngle < 0.5 * Math.PI ? plateAngle : 0.5 * Math.PI;
  const c = Math.cos(psi);
  if (!(c > a * c0)) return A0 * (1 - (2 / Math.PI) * shaft);
  const k = c / c0;
  const ak = (a * c0) / c; // a/k ≤ 1
  const inner = a * Math.sqrt(Math.max(0, 1 - ak * ak)) + k * Math.asin(ak);
  return A0 * (1 - k + (2 / Math.PI) * (inner - shaft));
}

/**
 * Throttle flow area at opening u ∈ [0, 1] (0 = closed, 1 = wide open), m²: geometric open
 * area + leakage. Multiply by the throttle C_D for the effective orifice area.
 */
export function throttleArea(diameter: number, opening: number, opts: ThrottleOptions = {}): number {
  const a = opts.shaftRatio ?? THROTTLE_DEFAULTS.shaftRatio;
  const psi0 = opts.closedAngle ?? THROTTLE_DEFAULTS.closedAngle;
  const leak = opts.leakageArea ?? Math.PI * diameter * THROTTLE_DEFAULTS.leakageClearance;
  return throttleGeometricArea(diameter, throttlePlateAngle(opening, psi0), a, psi0) + leak;
}

/**
 * Venturi throat depression p₀ − p_t ≈ ṁ²/(2ρ(C_D A_t)²), Pa (incompressible; valid while
 * Δp_t ≪ p₀). Signed like ṁ.
 * @param mdot mass flow, kg/s; @param rho upstream density, kg/m³
 * @param throatDiameter m; @param cd throat discharge coefficient (default 1)
 */
export function venturiThroatDepression(mdot: number, rho: number, throatDiameter: number, cd = 1): number {
  const A = cd * 0.25 * Math.PI * throatDiameter * throatDiameter;
  return (mdot * Math.abs(mdot)) / (2 * rho * A * A);
}

/**
 * Unrecovered stagnation-pressure loss across a venturi, Pa (signed like ṁ):
 * lossFraction × throat depression. lossFraction default 0.15 (UNVERIFIED, see header).
 */
export function venturiPressureLoss(mdot: number, rho: number, throatDiameter: number, lossFraction = 0.15, cd = 1): number {
  return lossFraction * venturiThroatDepression(mdot, rho, throatDiameter, cd);
}
