/**
 * Turbulent flame: eddy-entrainment + laminar burn-up model (Blizard & Keck 1974; Tabaczynski,
 * Ferguson & Radhakrishnan 1977; in the form derived by Beretta, Rashidi & Keck from flame photography
 * and written by Keck 1982).
 *
 * Keck 1982, "Turbulent flame structure and speed in spark-ignition engines", 19th Symp. (Int.) on
 * Combustion, pp. 1451–1466 (fetched; equation numbers from that paper):
 *   dm_b/dt = ρ_u A_f s_ℓ + μ/τ_b                (4.1A)
 *   dμ/dt   = ρ_u A_f u_T − μ/τ_b                (4.1B)
 *   τ_b     = ℓ_T / s_ℓ                           (text below 4.2; ℓ_T = s_ℓ τ_b, eq. 3.7)
 * where μ is the mass entrained into the flame but not yet burned. With μ = m_e − m_b this is
 *   dm_e/dt = ρ_u A_f (u_T + s_ℓ),   dm_b/dt = ρ_u A_f s_ℓ + (m_e − m_b)/τ_b,
 * the formulation in DESIGN.md. Closures (Tabaczynski et al. 1977; Tabaczynski, Trinker & Shannon 1980,
 * Combust. Flame 39:111–121, abstract fetched: burned vortex tubes "propagate at a rate equal to
 * U′ + S_L"): the characteristic speed u_T is the turbulence intensity u′, and the burn-up length ℓ_T is
 * the Taylor microscale λ (the spacing of the Kolmogorov-scale vortex tubes), so τ_b = λ/S_L.
 * Keck's limiting cases, all exact consequences of the ODEs (tested):
 *   quiescent (u′ → 0):            s_b = s_ℓ                                    (4.6)
 *   quasi-steady (dμ/dt = 0):      s_b = u_T + s_ℓ,  μ = ρ_u A_f u_T τ_b        (4.7)
 *   initial (t ≪ τ_b):             s_b/s_ℓ = 1 + t/τ_T (constant A_f; t/(3τ_T) for spherical growth) (4.8)
 * with s_b = (dm_b/dt)/(ρ_u A_f) (eq. 3.14) and τ_T = ℓ_T/u_T. The steady flame-brush thickness is
 * (V_f − V_b)/A_f = μ/(ρ_u A_f) = u′τ_b (Keck eq. 4.3: r_f − r_b → u_T τ_b).
 *
 * Geometry coupling: the enflamed volume V_f = m_b/ρ_b + (m_e − m_b)/ρ_u gives the front radius through
 * FlameGeometry.radiusForVolume, and A_f = FlameGeometry.evaluate(r_f, h).frontArea (= dV_f/dr_f).
 * Rapid distortion of the unburned gas ahead of the front (u′ ∝ ρ_u^{1/3}, L ∝ ρ_u^{-1/3}; Tabaczynski
 * et al. 1977, Wong & Hoult 1979) is carried by the turbulence model: advance the K–k equations with
 * dlnRhoDt = d ln ρ_u/dt and take L from angularMomentumLengthScale (turbulence.ts).
 * Keck's own empirical closures (u_T ∝ ρ_u^{1/2}, ℓ_T ∝ ρ_u^{-3/4}; eq. 4.10, Fig. 15) are provided as
 * {@link keckCharacteristicSpeed} / {@link keckCharacteristicLength} for comparison and calibration.
 */

/** Inputs to {@link entrainmentRates} (reuse one object). SI. */
export interface EntrainmentInputs {
  /** Unburned-gas density ahead of the front, kg/m³. */
  rhoU: number;
  /** Flame-front area A_f, m². */
  frontArea: number;
  /** Turbulence intensity u′ of the unburned gas, m/s. */
  uPrime: number;
  /** Laminar burning velocity S_L of the unburned mixture, m/s. */
  SL: number;
  /** Entrained mass m_e, kg. */
  me: number;
  /** Burned mass m_b, kg. */
  mb: number;
  /**
   * Burn-up length: Taylor microscale λ, m (see {@link taylorMicroscale}). Infinity (u′ = 0) ⇒ no burn-up
   * term (laminar limit); λ ≤ 0 or NaN is treated the same way (invalid input, never instant burn-up).
   */
  lambda: number;
  /**
   * Total combustible mass (entrainment stops when m_e reaches it), kg. Optional (not part of the
   * DESIGN.md contract); omitted or Infinity = no cap. The front area → 0 at the far corner stops
   * entrainment continuously anyway; this is only a guard against integrator overshoot.
   */
  mTotal?: number;
}

/** Outputs of {@link entrainmentRates}. SI. */
export interface EntrainmentRates {
  /** dm_e/dt, kg/s. */
  dme: number;
  /** dm_b/dt, kg/s. */
  dmb: number;
  /** Burn-up time τ_b = λ/S_L, s (Infinity if S_L = 0). */
  tauB: number;
  /** Entrainment speed u′ + S_L, m/s. */
  entrainmentSpeed: number;
  /** Burning speed s_b = (dm_b/dt)/(ρ_u A_f), m/s (Keck 1982 eq. 3.14); 0 when A_f = 0. */
  burningSpeed: number;
}

/** Allocate an {@link EntrainmentInputs} (no flame). */
export const newEntrainmentInputs = (): EntrainmentInputs => ({
  rhoU: 1,
  frontArea: 0,
  uPrime: 0,
  SL: 0,
  me: 0,
  mb: 0,
  lambda: 1e-3,
  mTotal: Infinity,
});

/** Allocate an {@link EntrainmentRates}. */
export const newEntrainmentRates = (): EntrainmentRates => ({
  dme: 0,
  dmb: 0,
  tauB: Infinity,
  entrainmentSpeed: 0,
  burningSpeed: 0,
});

/**
 * Entrainment and burning rates (Keck 1982 eqs. 4.1A/B with u_T = u′, ℓ_T = λ). Pure, allocation-free;
 * writes `out` and returns it.
 */
export function entrainmentRates(inp: EntrainmentInputs, out: EntrainmentRates): EntrainmentRates {
  const Af = inp.frontArea > 0 ? inp.frontArea : 0;
  const rhoU = inp.rhoU;
  const SL = inp.SL > 0 ? inp.SL : 0;
  const uP = inp.uPrime > 0 ? inp.uPrime : 0;
  const uE = uP + SL;
  const flux = rhoU * Af; // kg/(m s)
  let dme = flux * uE;
  const mTot = inp.mTotal;
  if (mTot !== undefined && inp.me >= mTot) dme = 0; // everything entrained (front has reached the far walls)
  const tauB = SL > 0 && inp.lambda > 0 ? inp.lambda / SL : Infinity;
  const mu = inp.me > inp.mb ? inp.me - inp.mb : 0;
  const dmb = flux * SL + (tauB < Infinity ? mu / tauB : 0);
  out.dme = dme;
  out.dmb = dmb;
  out.tauB = tauB;
  out.entrainmentSpeed = uE;
  out.burningSpeed = flux > 0 ? dmb / flux : 0;
  return out;
}

/**
 * Taylor microscale λ, m. For isotropic turbulence ε = 15 ν u′²/λ² (Taylor 1935); with the model's
 * dissipation ε = A u′³/L this gives
 *   λ/L = (15/A)^{1/2} Re_L^{-1/2},  Re_L = u′L/ν
 * (Tennekes & Lumley 1972 §3.2 — UNVERIFIED section; Xu & Filipi 2020 eq. 30 — fetched). A must be
 * the constant of the
 * dissipation closure actually used: A = 1 with the K–k model's ε = u′³/L (turbulence.ts), which makes λ
 * exactly the microscale implied by the modelled ε.
 * @param L dissipation length scale, m
 * @param uPrime turbulence intensity, m/s
 * @param nu kinematic viscosity of the unburned gas, m²/s
 * @param A dissipation constant (default 1)
 * @returns λ (Infinity when u′ = 0)
 */
export function taylorMicroscale(L: number, uPrime: number, nu: number, A = 1): number {
  if (!(uPrime > 0)) return Infinity;
  return Math.sqrt((15 * nu * L) / (A * uPrime));
}

/**
 * Kolmogorov length η = (ν³/ε)^{1/4} with ε = A u′³/L, m (Infinity when u′ = 0).
 * @param L dissipation length scale, m
 * @param uPrime turbulence intensity, m/s
 * @param nu kinematic viscosity, m²/s
 * @param A dissipation constant (default 1)
 */
export function kolmogorovScale(L: number, uPrime: number, nu: number, A = 1): number {
  if (!(uPrime > 0)) return Infinity;
  const eps = (A * uPrime * uPrime * uPrime) / L;
  return Math.sqrt(Math.sqrt((nu * nu * nu) / eps));
}

/**
 * Effective turbulent burning velocity for display: s_b = (dm_b/dt)/(ρ_u A_f) (Keck 1982 eq. 3.14),
 * m/s. Equals S_L for a laminar flame and tends to u′ + S_L for a fully developed flame brush.
 */
export function effectiveTurbulentBurningVelocity(inp: EntrainmentInputs): number {
  const Af = inp.frontArea;
  if (!(Af > 0) || !(inp.rhoU > 0)) return 0;
  const SL = inp.SL > 0 ? inp.SL : 0;
  const tauB = SL > 0 && inp.lambda > 0 ? inp.lambda / SL : Infinity;
  const mu = inp.me > inp.mb ? inp.me - inp.mb : 0;
  return SL + (tauB < Infinity ? mu / (tauB * inp.rhoU * Af) : 0);
}

/**
 * Steady flame-brush thickness δ = u′ τ_b = u′λ/S_L, m (Keck 1982 eq. 4.3 asymptote, quasi-steady μ).
 */
export function steadyFlameBrushThickness(uPrime: number, lambda: number, SL: number): number {
  return SL > 0 ? (uPrime * lambda) / SL : Infinity;
}

/**
 * Enflamed volume V_f = m_b/ρ_b + (m_e − m_b)/ρ_u, m³ — burned gas plus entrained-but-unburned gas.
 * Invert with FlameGeometry.radiusForVolume to get the front radius.
 */
export function enflamedVolume(me: number, mb: number, rhoB: number, rhoU: number): number {
  const mu = me > mb ? me - mb : 0;
  return mb / rhoB + mu / rhoU;
}

/**
 * Keck (1982) eq. 4.10 empirical characteristic speed u_T = 0.08 ū_i (ρ_u/ρ_i)^{1/2}, m/s (±10 % over
 * the engines of his Fig. 14), with ū_i the mean inlet-gas speed ({@link keckMeanInletSpeed}).
 * @param meanInletSpeed ū_i, m/s
 * @param rhoU unburned density at the time of interest, kg/m³
 * @param rhoInlet inlet density ρ_i, kg/m³
 */
export function keckCharacteristicSpeed(meanInletSpeed: number, rhoU: number, rhoInlet: number): number {
  return 0.08 * meanInletSpeed * Math.sqrt(rhoU / rhoInlet);
}

/**
 * Keck (1982) Fig. 15 correlation ℓ_T = 0.8 L_IV (ρ_i/ρ_u)^{3/4}, m (L_IV = maximum inlet-valve lift).
 * @param maxIntakeLift L_IV, m
 * @param rhoU unburned density, kg/m³
 * @param rhoInlet inlet density ρ_i, kg/m³
 */
export function keckCharacteristicLength(maxIntakeLift: number, rhoU: number, rhoInlet: number): number {
  return 0.8 * maxIntakeLift * Math.pow(rhoInlet / rhoU, 0.75);
}

/**
 * Keck (1982) mean inlet-gas speed ū_i = ε_v (A_p/A_IV) 2 N S, m/s (text above eq. 4.10).
 * @param volumetricEfficiency ε_v
 * @param pistonArea A_p, m²
 * @param maxValveArea maximum open area of the inlet valve A_IV, m²
 * @param rpm engine speed N, rev/min
 * @param stroke S, m
 */
export function keckMeanInletSpeed(
  volumetricEfficiency: number,
  pistonArea: number,
  maxValveArea: number,
  rpm: number,
  stroke: number,
): number {
  return (volumetricEfficiency * (pistonArea / maxValveArea) * 2 * rpm * stroke) / 60;
}
