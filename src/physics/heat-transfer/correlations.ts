/**
 * Spatially averaged in-cylinder gas-to-wall heat-transfer correlations (quasi-steady,
 * Newton's law q = h (T − T_w)). All inputs/outputs SI (m, Pa, K, m/s, m³; h in W/(m² K));
 * the published unit systems are converted internally.
 *
 * ── Woschni (1967), SAE 670931 ──────────────────────────────────────────────────────────
 * As quoted by "Revising engine heat transfer", Annals of the Faculty of Engineering
 * Hunedoara VI(3), 2008, eq. (4) (fetched):
 *   h = 0.820 D^−0.2 p^0.8 w^0.8 T^−0.53   [kW/(m² K); D m, p MPa, T K, w m/s]
 *     = 820 · D^−0.2 · (p/1 MPa)^0.8 · w^0.8 · T^−0.53   W/(m² K)
 *   w = C₁ c_m + C₂ (V_s T₁/(p₁ V₁)) (p − p_mot)
 *   C₁ = 6.18 (gas exchange), 2.28 (compression, combustion, expansion);
 *   swirl variant: C₁ = 6.18 + 0.417 c_u/c_m and 2.28 + 0.308 c_u/c_m, with c_u = B ω_p/2 the
 *   tangential speed of the swirl-meter paddle wheel (ω_p its angular speed; ≈ ω_swirl·B/2 for a
 *   solid-body swirl) [ScienceDirect Topics "Average Heat Transfer Coefficient" (Heywood-based),
 *   search-result text, 2026-09-29: "C1 = 6.18 + 0.417 v_s/S_p … 2.28 + 0.308 v_s/S_p, v_s = Bω_p/2"];
 *   C₂ = 0 (gas exchange, compression), 3.24e-3 m/(s K) (combustion, expansion),
 *   6.22e-3 m/(s K) for divided-chamber (IDI) engines [xarin.com catoolRT docs, fetched];
 *   (p₁, T₁, V₁) a reference state (IVC or start of combustion), p_mot the motored pressure at
 *   the same crank angle.
 * The same review states Woschni's derivation assumed k ∝ T^0.75 and μ ∝ T^0.62.
 *
 * Heywood (1988) eq. 12.19 writes the correlation as h = 3.26 B^−0.2 p^0.8 T^−0.55 w^0.8 with p in
 * kPa (form quoted by the Chalmers thesis research.chalmers.se/publication/508534, eq. 12, fetched;
 * UNVERIFIED: the Heywood equation number). The constant is the same (820/1000^0.8 = 3.2645) but
 * the temperature exponent differs: at 1000–2500 K the Heywood form is 13–15 % lower. The default
 * here is the ORIGINAL Woschni exponent −0.53 (variant 'woschni1967'); 'heywood1988' is provided.
 *
 * ── Hohenberg (1979), SAE 790825 ────────────────────────────────────────────────────────
 *   h = 130 · V^−0.06 · (p/1 bar)^0.8 · T^−0.4 · (c_m + 1.4)^0.8   W/(m² K), V m³ (instantaneous)
 * UNVERIFIED: the leading constant — 130 from memory of SAE 790825; the xarin.com catoolRT
 * documentation (fetched again 2026-09-29) prints 129.8 (0.15 % lower — immaterial next to the
 * correlation's own scatter). Exponents and the 1.4 m/s offset agree with that page.
 *
 * Unit cross-check (review): Woschni's original 110 kcal/(m² h K) with p in at (kp/cm²) gives
 * 110 · 1.163 · (1/0.980665)^0.8 = 129.9 with p in bar = 0.8199 kW/(m² K) with p in MPa
 * = 3.264 W/(m² K) with p in kPa — the three published forms agree; UNVERIFIED: the original
 * kcal/at form (from memory), the bar/MPa/kPa forms are the fetched ones above.
 *
 * ── Annand (1963), Proc. IMechE 177(36) ────────────────────────────────────────────────
 *   q = a (k/B) Re^b (T − T_w) + c (T⁴ − T_w⁴),   Re = ρ c_m B/μ,   b = 0.7,
 *   a = 0.35 … 0.8 (increasing with charge-motion intensity),
 *   c = 0 (intake, compression), 0.075 σ (SI combustion/expansion), 0.576 σ (diesel)
 * [xarin.com catoolRT docs, fetched]. Gas k, μ, ρ at the bulk temperature (thermo/transport).
 * UNVERIFIED: the default a = 0.49 (a commonly used mid-range value; not checked in Annand 1963).
 *
 * Hot paths allocate nothing (callers should reuse the input objects).
 */
import { SIGMA_SB } from '../core/constants';

/** Woschni parameter sets. */
export type WoschniVariant = 'woschni1967' | 'heywood1988';

/**
 * Cycle phase for the Woschni velocity constants: 'gas-exchange' (valves open),
 * 'compression' (closed, before combustion), 'combustion' (combustion and expansion);
 * 'expansion' is accepted as a synonym of 'combustion' (Woschni's C₂ term applies to both), so
 * core/snapshot CylinderPhase values can be passed directly.
 */
export type WoschniPhase = 'gas-exchange' | 'compression' | 'combustion' | 'expansion';

/** Woschni C₁ (gas exchange / closed cycle) and swirl slopes, C₂ values (m/(s K)). */
export const WOSCHNI_CONSTANTS = Object.freeze({
  c1GasExchange: 6.18,
  c1Closed: 2.28,
  swirlGasExchange: 0.417,
  swirlClosed: 0.308,
  c2Combustion: 3.24e-3,
  c2DividedChamber: 6.22e-3,
  /** 0.820 kW/(m² K) with p in MPa → W/(m² K) with p in MPa. */
  cWoschni1967: 820,
  tExpWoschni1967: -0.53,
  /** 3.26 W/(m² K) with p in kPa (Heywood 1988). */
  cHeywood: 3.26,
  tExpHeywood: -0.55,
});

/** Inputs of {@link woschniVelocity} / {@link woschniCoefficient} (reuse one object; SI). */
export interface WoschniInputs {
  /** Cylinder bore D, m. */
  bore: number;
  /** Cylinder pressure p, Pa. */
  pressure: number;
  /** Mass-averaged gas temperature T, K. */
  temperature: number;
  /** Mean piston speed c_m = 2 S N, m/s. */
  meanPistonSpeed: number;
  /** Cycle phase (sets C₁, C₂). */
  phase: WoschniPhase;
  /** Motored pressure at the same crank angle, Pa (required for phase 'combustion'). */
  motoredPressure?: number;
  /** Displaced volume V_s, m³ (required for phase 'combustion'). */
  displacedVolume?: number;
  /** Reference-state pressure p₁ (e.g. IVC), Pa (required for phase 'combustion'). */
  refPressure?: number;
  /** Reference-state temperature T₁, K (required for phase 'combustion'). */
  refTemperature?: number;
  /** Reference-state volume V₁, m³ (required for phase 'combustion'). */
  refVolume?: number;
  /** Swirl ratio c_u/c_m (paddle-wheel tangential speed / mean piston speed); default 0. */
  swirlRatio?: number;
  /** C₂ override, m/(s K) (e.g. WOSCHNI_CONSTANTS.c2DividedChamber). */
  c2?: number;
  /** Correlation variant (default 'woschni1967'). */
  variant?: WoschniVariant;
}

/**
 * Woschni characteristic gas velocity w, m/s (≥ 0; the combustion term is clipped at p < p_mot).
 * Phase 'combustion' (or 'expansion') REQUIRES motoredPressure, displacedVolume and the reference state
 * (refPressure, refTemperature, refVolume, all > 0): a missing value would silently drop the
 * combustion term (≈ halving h after ignition), so it throws a RangeError instead.
 */
export function woschniVelocity(i: WoschniInputs): number {
  const s = i.swirlRatio ?? 0;
  const K = WOSCHNI_CONSTANTS;
  if (i.phase === 'gas-exchange') return (K.c1GasExchange + K.swirlGasExchange * s) * i.meanPistonSpeed;
  let w = (K.c1Closed + K.swirlClosed * s) * i.meanPistonSpeed;
  if (i.phase === 'combustion' || i.phase === 'expansion') {
    const pm = i.motoredPressure;
    const Vs = i.displacedVolume;
    const p1 = i.refPressure;
    const T1 = i.refTemperature;
    const V1 = i.refVolume;
    if (!(pm !== undefined && pm >= 0 && Vs !== undefined && Vs > 0 && p1 !== undefined && p1 > 0 && T1 !== undefined && T1 > 0 && V1 !== undefined && V1 > 0)) {
      throw new RangeError('woschniVelocity: phase "combustion" needs motoredPressure, displacedVolume, refPressure, refTemperature, refVolume');
    }
    const dp = i.pressure - pm;
    if (dp > 0) {
      const c2 = i.c2 ?? K.c2Combustion;
      w += (c2 * Vs * T1 * dp) / (p1 * V1);
    }
  }
  return w;
}

/** Woschni heat-transfer coefficient h, W/(m² K) (see file header for the variants). */
export function woschniCoefficient(i: WoschniInputs): number {
  const w = woschniVelocity(i);
  const K = WOSCHNI_CONSTANTS;
  const p = i.pressure > 0 ? i.pressure : 0;
  if (i.variant === 'heywood1988') {
    return K.cHeywood * Math.pow(i.bore, -0.2) * Math.pow(p * 1e-3, 0.8) * Math.pow(i.temperature, K.tExpHeywood) * Math.pow(w, 0.8);
  }
  return K.cWoschni1967 * Math.pow(i.bore, -0.2) * Math.pow(p * 1e-6, 0.8) * Math.pow(i.temperature, K.tExpWoschni1967) * Math.pow(w, 0.8);
}

/**
 * Polytropic motored-pressure estimate p_mot = p_ref (V_ref/V)^n, Pa, for Woschni's combustion
 * term when no motored trace is simulated. n is the polytropic exponent of the motored
 * compression/expansion (typically 1.30–1.35; caller's choice).
 */
export function woschniMotoredPressure(refPressure: number, refVolume: number, volume: number, n: number): number {
  return refPressure * Math.pow(refVolume / volume, n);
}

/**
 * Hohenberg (1979) coefficient h, W/(m² K).
 * @param volume instantaneous cylinder volume, m³
 * @param pressure Pa; @param temperature bulk gas temperature, K; @param meanPistonSpeed m/s
 */
export function hohenbergCoefficient(volume: number, pressure: number, temperature: number, meanPistonSpeed: number): number {
  const p = pressure > 0 ? pressure : 0;
  return 130 * Math.pow(volume, -0.06) * Math.pow(p * 1e-5, 0.8) * Math.pow(temperature, -0.4) * Math.pow(meanPistonSpeed + 1.4, 0.8);
}

/** Annand radiation constants c (W/(m² K⁴)): 0.075σ (SI), 0.576σ (diesel). */
export const ANNAND_RADIATION = Object.freeze({
  sparkIgnition: 0.075 * SIGMA_SB,
  diesel: 0.576 * SIGMA_SB,
});

/** Default Annand convective constants (a UNVERIFIED mid-range value, b = 0.7). */
export const ANNAND_DEFAULTS = Object.freeze({ a: 0.49, b: 0.7 });

/**
 * Annand (1963) convective coefficient h = a (k/B) Re^b, W/(m² K), Re = ρ c_m B/μ.
 * @param bore m; @param meanPistonSpeed m/s; @param density kg/m³; @param viscosity Pa s;
 * @param conductivity W/(m K); @param a, b correlation constants
 */
export function annandCoefficient(
  bore: number,
  meanPistonSpeed: number,
  density: number,
  viscosity: number,
  conductivity: number,
  a: number = ANNAND_DEFAULTS.a,
  b: number = ANNAND_DEFAULTS.b,
): number {
  const Re = (density * meanPistonSpeed * bore) / viscosity;
  return ((a * conductivity) / bore) * Math.pow(Re > 0 ? Re : 0, b);
}

/**
 * Annand wall heat flux q = h (T − T_w) + c (T⁴ − T_w⁴), W/m² (positive gas → wall).
 * @param h convective coefficient from {@link annandCoefficient}, W/(m² K)
 * @param c radiation constant, W/(m² K⁴) (0 during intake/compression)
 */
export function annandFlux(h: number, T: number, Tw: number, c = 0): number {
  const T2 = T * T;
  const W2 = Tw * Tw;
  return h * (T - Tw) + c * (T2 * T2 - W2 * W2);
}
