/**
 * Spark-gap breakdown voltage (the "required voltage" of the plug), SI units (m, Pa, K, V).
 *
 * At engine conditions the gap sits far on the RIGHT side of the Paschen curve
 * (p·d ≈ 10–60 bar·mm ≫ the Paschen minimum at ≈ 0.0075 bar·mm), where breakdown is a
 * streamer process and the breakdown voltage is almost linear in the reduced gas density
 * N·d (number density × gap). Three formulations are provided:
 *
 * 1. `breakdownVoltageUniformAir` (DEFAULT) — the semi-empirical uniform-field law for air
 *        V_b = E_c (p* d) + K √(p* d),     E_c = 24.36 kV/(cm bar),  K = 6.72 kV/(cm bar)^½,
 *    p* = p (T₀/T) the pressure reduced to T₀ = 293.15 K, in BAR, d in cm ("p in bar and d in
 *    cm, at 20 °C" — Kuffel 1961 / Dakin et al. 1974 as given in Kuffel, Zaengl & Kuffel,
 *    "High Voltage Engineering: Fundamentals", 2nd ed. 2000, §5). The reference is 1 bar, not
 *    1 atm: an earlier revision used the relative air density δ (1 atm), which lowered V_b by
 *    ≈ 1.1 % at engine densities (reviewer fix). The form follows from the streamer
 *    criterion ∫(α − η)dx = const with an effective ionisation coefficient
 *    (α − η)/N ∝ (E/N − (E/N)_c)² about the critical reduced field (E/N)_c ≈ 99 Td where
 *    ionisation balances attachment.
 *    UNVERIFIED against the textbook itself (not fetched): the coefficients and the "p in bar,
 *    d in cm, 20 °C" units are confirmed only by secondary summaries of Kuffel (1961). An
 *    independent, fetched set (Schumann form, relative density δ at 1 atm)
 *    U_b = 24.22 δd + 6.08 √(δd) kV (J. Li, Chongqing Univ., "Fundamentals of High Voltage
 *    Engineering", lecture 4-2) lies 4–7 % lower over δd = 0.05–1 cm, which bounds the
 *    uncertainty. Checks: 1 atm, 20 °C, 1 mm → 4.61 kV (the air Paschen curve gives
 *    ≈ 4.5 kV); 1 cm → 31.4 kV; valid for p*d ≳ 0.01 bar·cm (right of the Paschen minimum).
 *
 * 2. `breakdownVoltagePaschen` — classical Townsend/Paschen law with the Townsend first
 *    ionisation coefficient α/p = A exp(−B p/E):
 *        V_b = B (pd) / ln[A (pd) / ln(1 + 1/γ)],
 *    A = 14.6 /(cm Torr), B = 365 V/(cm Torr) for air (Cobine 1942, "Gaseous Conductors",
 *    as quoted by Shaffer, Zare & Askari, IMECE2021-73138, eqs. 1–2), γ = secondary-electron
 *    yield (0.01–0.02 for iron/steel). pd is evaluated at the reference density (p T₀/T) d
 *    because α depends on E/N. UNVERIFIED: the temperature to which Cobine's pressures are
 *    reduced (0 °C, common in the older literature, would lower V_b by ≈ 6 %); T₀ = 20 °C is
 *    used. At engine densities E/p ≈ 55 V/(cm Torr), below the 100–800 V/(cm Torr) fit
 *    range of A, B. Neglects attachment, so it overestimates V_b by ≈ 10–20 % at engine
 *    δd; kept as an independent first-principles cross-check.
 *
 * 3. `breakdownVoltagePashley` — engine/bomb correlation of Pashley, Stone & Roberts (2000),
 *    SAE 2000-01-0245, as quoted by the COBEM 2011 paper "A standard inductive ignition
 *    versus a high energy capacitive system" (ABCM, Proc. COBEM 2011, eq. 2,
 *    abcm.org.br/anais/cobem/2011/PDF/104101.pdf):
 *        V_b[kV] = 4.3 + 136 p/T + 324 p d/T      (p in bar, T in K, d in mm).
 *    The intercept is 4,3 (Brazilian decimal comma). Reviewer fix: an earlier revision read
 *    the garbled text layer as "3,4". Placing the glyphs by their PDF text-matrix
 *    x-positions gives the order "V = 4 , 3 + 136 p/T + 324 p d/T" (x = 77.6, 99.5, 105.4,
 *    108.2, 123.9, 143.8, 164.3, 184.0 pt). An independent quote of Pashley et al. has the
 *    same intercept, "V = 4.3 + (146·P/T) + (324·P/T·Dg)" (forum post citing SAE 2000-01-0245).
 *    UNVERIFIED: the p/T coefficient, 136 (COBEM) vs 146 (forum). The SAE paper itself was
 *    not accessible.
 *
 * Electrode temperature and polarity (not modelled quantitatively): ignition coils are wound
 * so that the hot centre electrode is the CATHODE (negative). Thermionic emission from the
 * hotter cathode aids the initial avalanche, so negative-centre polarity breaks down at a
 * somewhat lower voltage than positive polarity (standard practice, e.g. Bosch; the
 * magnitude — typically quoted as 10–20 % — is UNVERIFIED). Hot electrodes also heat the gap
 * gas, lowering N·d; the caller can pass the local gas temperature. Fast voltage rise causes
 * an impulse over-voltage (statistical + formative time lag): Michler et al. (KIT, "Influence
 * of the electrical parameters of the ignition system on the phases of spark ignition",
 * §4.1) measured 5 kV vs 3.8 kV static Paschen (factor 1.32) at 1 bar, 0.9 mm; this is
 * exposed as `impulseFactor` (default 1 = static breakdown).
 *
 * Composition: V_b depends on N·d and weakly on the gas (fuel vapour, CO₂/H₂O residuals).
 * The mole fractions `X` are accepted for API stability but no composition correction is
 * applied (air-equivalent dielectric strength) — see module report.
 */

import { K_BOLTZMANN } from '../core/constants';

/** Reference pressure of the relative air density δ, Pa (1 atm). */
export const RELATIVE_DENSITY_P0 = 101325;
/** Reference temperature of the relative air density δ, K (20 °C, HV-engineering convention). */
export const RELATIVE_DENSITY_T0 = 293.15;
/** Pressure unit of the Kuffel uniform-field air law, Pa (1 bar; see module doc). */
export const AIR_UNIFORM_P_UNIT = 1e5;

/**
 * Critical field of the uniform-field air law, V/m per bar (24.36 kV/(cm bar)) at 20 °C.
 * Kuffel (1961) via Kuffel, Zaengl & Kuffel (2000) §5 (see module doc).
 */
export const AIR_UNIFORM_EC = 24.36e5;
/**
 * √-term coefficient of the uniform-field air law, V/√m per √bar
 * (6.72 kV/(cm bar)^½ = 6.72e3 V/√(0.01 m) = 6.72e4 V/√m). Kuffel et al. (2000) §5.
 */
export const AIR_UNIFORM_K = 6.72e4;

/** Townsend A coefficient for air, 1/(m Pa)  (14.6 /(cm Torr); Cobine 1942 via Shaffer et al. 2021). */
export const TOWNSEND_A_AIR = 14.6 / (1e-2 * 133.322368);
/** Townsend B coefficient for air, V/(m Pa)  (365 V/(cm Torr); Cobine 1942 via Shaffer et al. 2021). */
export const TOWNSEND_B_AIR = 365 / (1e-2 * 133.322368);
/**
 * Secondary-electron emission yield γ for iron/steel cathodes in air (0.02, the value used by
 * Shaffer, Zare & Askari, IMECE2021-73138, with Cobine's A, B → normal cathode fall 295 V).
 */
export const GAMMA_SE_IRON = 0.02;

/** Available breakdown formulations. */
export type BreakdownModel = 'uniform-air' | 'paschen' | 'pashley';

/** Options of `breakdownVoltage`. */
export interface BreakdownOptions {
  /** Formulation (default 'uniform-air'). */
  model?: BreakdownModel;
  /** Impulse (dynamic) over-voltage factor ≥ 1 applied to the static value (default 1). */
  impulseFactor?: number;
  /** Secondary-electron yield for the 'paschen' model (default 0.02). */
  gammaSE?: number;
}

/** Relative air density δ = (p/p₀)(T₀/T) (dimensionless), p in Pa, T in K. */
export function relativeAirDensity(p: number, T: number): number {
  return (p / RELATIVE_DENSITY_P0) * (RELATIVE_DENSITY_T0 / T);
}

/** Gas number density N = p/(k_B T), 1/m³. */
export function numberDensity(p: number, T: number): number {
  return p / (K_BOLTZMANN * T);
}

/**
 * Pressure reduced to 20 °C in bar, p* = (p / 1 bar)(T₀/T) (dimensionless), p in Pa, T in K —
 * the density variable of the Kuffel uniform-field law (≡ 1.01325 δ).
 */
export function reducedPressureBar(p: number, T: number): number {
  return (p / AIR_UNIFORM_P_UNIT) * (RELATIVE_DENSITY_T0 / T);
}

/**
 * Uniform-field (streamer-regime) breakdown voltage of air, V:
 *   V_b = E_c p*d + K √(p*d),  p* in bar at 20 °C (Kuffel et al. 2000 §5; see module doc).
 * gap (m), p (Pa), T (K). Linear in N·d at large p*d; ≈ 3 kV/mm-scale at 1 atm.
 */
export function breakdownVoltageUniformAir(gap: number, p: number, T: number): number {
  const dd = reducedPressureBar(p, T) * gap; // p*·d, bar·m
  return AIR_UNIFORM_EC * dd + AIR_UNIFORM_K * Math.sqrt(dd);
}

/**
 * Classical Paschen/Townsend breakdown voltage, V:
 *   V_b = B pd / ln(A pd / ln(1 + 1/γ)),
 * pd evaluated at the reference density: (p T₀/T)·d (Pa m). Returns +Infinity left of the
 * Paschen asymptote (A pd ≤ ln(1 + 1/γ), no self-sustained Townsend breakdown).
 * Coefficients: Cobine (1942) via Shaffer et al. (2021) — see module doc.
 */
export function breakdownVoltagePaschen(
  gap: number,
  p: number,
  T: number,
  gammaSE: number = GAMMA_SE_IRON,
  A: number = TOWNSEND_A_AIR,
  B: number = TOWNSEND_B_AIR,
): number {
  const pd = p * (RELATIVE_DENSITY_T0 / T) * gap;
  const den = Math.log((A * pd) / Math.log(1 + 1 / gammaSE));
  return den > 0 ? (B * pd) / den : Infinity;
}

/** Intercept of the Pashley et al. (2000) correlation, kV (COBEM 2011 eq. 2: "4,3"). */
export const PASHLEY_INTERCEPT_KV = 4.3;
/** p/T coefficient of the Pashley correlation, kV K/bar. UNVERIFIED: 136 (COBEM 2011) vs 146 (forum quote). */
export const PASHLEY_P_COEFF = 136;
/** p·d/T coefficient of the Pashley correlation, kV K/(bar mm) (COBEM 2011 and forum quote agree). */
export const PASHLEY_PD_COEFF = 324;

/**
 * Pashley, Stone & Roberts (2000, SAE 2000-01-0245) linear engine correlation, V:
 *   V_b[kV] = 4.3 + 136 p/T + 324 p d/T (p bar, T K, d mm).
 * See module doc for the sources and the UNVERIFIED p/T coefficient.
 */
export function breakdownVoltagePashley(gap: number, p: number, T: number): number {
  const pBar = p * 1e-5;
  const dMm = gap * 1e3;
  return 1e3 * (PASHLEY_INTERCEPT_KV + (PASHLEY_P_COEFF * pBar) / T + (PASHLEY_PD_COEFF * pBar * dMm) / T);
}

/**
 * Normal-glow cathode fall from Townsend theory, V (von Engel; Cobine 1942 — quoted by
 * Shaffer, Zare & Askari, IMECE2021-73138, eq. 1):  V_n = 3 (B/A) ln(1 + 1/γ).
 * Air, γ = 0.02 → 295 V (Cobine's measured value for iron in air: 269 V).
 */
export function normalCathodeFall(
  gammaSE: number = GAMMA_SE_IRON,
  A: number = TOWNSEND_A_AIR,
  B: number = TOWNSEND_B_AIR,
): number {
  return 3 * (B / A) * Math.log(1 + 1 / gammaSE);
}

/**
 * Paschen-curve minimum of the Townsend law, V: V_min = e (B/A) ln(1 + 1/γ) (at
 * pd = e ln(1 + 1/γ)/A). Air on iron (Cobine A, B; γ = 0.02) → 267 V; the measured air
 * minimum is ≈ 330 V. No gas gap breaks down below it, whatever N·d.
 */
export function paschenMinimumVoltage(
  gammaSE: number = GAMMA_SE_IRON,
  A: number = TOWNSEND_A_AIR,
  B: number = TOWNSEND_B_AIR,
): number {
  return Math.E * (B / A) * Math.log(1 + 1 / gammaSE);
}

/** Paschen minimum of air on iron (γ = 0.02), V (≈ 267 V). */
export const PASCHEN_MINIMUM_AIR = paschenMinimumVoltage();

/**
 * Breakdown ("required") voltage of a spark gap, V.
 * The static value of every model is floored at the Paschen minimum (reviewer fix): the
 * uniform-field and Pashley laws are fits for the right branch of the Paschen curve and fall
 * below V_min at the very low N·d of a hot kernel at low pressure (e.g. 2500 K, 0.5 bar,
 * 0.2 mm → 0.24 kV), which previously triggered unbounded restrike cascades.
 * @param gap electrode gap, m
 * @param p gas pressure, Pa
 * @param T gas temperature in the gap, K (unburned-gas temperature for the first breakdown;
 *          the hot kernel temperature for a restrike)
 * @param _X mole fractions (accepted, not used: air-equivalent dielectric strength)
 * @param opts model selection and impulse factor
 */
export function breakdownVoltage(
  gap: number,
  p: number,
  T: number,
  _X?: Float64Array,
  opts?: BreakdownOptions,
): number {
  const model = opts?.model ?? 'uniform-air';
  const f = opts?.impulseFactor ?? 1;
  let v: number;
  if (model === 'paschen') v = breakdownVoltagePaschen(gap, p, T, opts?.gammaSE ?? GAMMA_SE_IRON);
  else if (model === 'pashley') v = breakdownVoltagePashley(gap, p, T);
  else v = breakdownVoltageUniformAir(gap, p, T);
  if (!(v >= PASCHEN_MINIMUM_AIR)) v = PASCHEN_MINIMUM_AIR;
  return f * v;
}
