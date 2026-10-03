/**
 * 0-D in-cylinder turbulence: the K–k energy-cascade model of Poulos & Heywood (1983), with the
 * rapid-distortion compression source and an optional solid-body swirl equation.
 *
 * State (extensive, SI):
 *   K = ½ m U²      mean-flow kinetic energy, J (U = mean-flow speed)
 *   k = (3/2) m u′²  turbulent kinetic energy, J (u′ = rms turbulence intensity, isotropic)
 *   H = I ω_s        swirl angular momentum about the cylinder axis, kg m²/s (I = m B²/8, solid body)
 *
 * Equations (Poulos & Heywood 1983, SAE 830334, as reproduced by Xu & Filipi 2020, Front. Mech. Eng.
 * 6:46, eqs. 32–39 — re-verified by the reviewer from the JATS XML of that article; the original SAE text
 * could not be obtained):
 *   dK/dt = ½ ṁ_in v_in² − P − K ṁ_out/m                               (eq. 35)
 *   dk/dt = P − m ε − k ṁ_out/m  [+ (2/3) k d(ln ρ)/dt, see below]      (eq. 36)
 *   P = 0.3307 c_β (K/L) √(k/m),    ε = u′³/L = (2k/3m)^{3/2}/L        (eqs. 37, 38)
 *   L = 4V/(πB²) (instantaneous chamber height)                        (eq. 32; replaced here, see below)
 *   during combustion: L/L₀ = (ρ_u0/ρ_u)^{1/3}, u′/u′₀ = (ρ_u/ρ_u0)^{1/3}  (eq. 39)
 * The constant 0.3307 is exactly the eddy-viscosity production of a mean shear U/L:
 *   P = m ν_t (U/L)²,  ν_t = C_μ k_s²/ε = (9/4) C_μ u′L  ⇒  P = 4.5 C_μ √(2/3) (K/L)√(k/m) = 0.3307 (K/L)√(k/m)
 * with the standard C_μ = 0.09 (Launder & Spalding 1974), so c_β = 1 is the "no-tuning" value. (The task
 * brief's "0.3875 β" appears in no source found; 0.3307 is both the published value and the derivation.)
 * UNVERIFIED: C_μ = 0.09 (standard k–ε constant) not re-fetched from Launder & Spalding; it only explains
 * the published 0.3307, which is verified.
 *
 * Rapid distortion: for a compression much faster than the eddy turnover, Kelvin's circulation theorem
 * applied to an isotropically compressed eddy of fixed mass (length ∝ ρ^{-1/3}, circulation conserved)
 * gives u′ ∝ ρ^{1/3}, i.e. d(k/m)/dt = (2/3)(k/m) d ln ρ/dt — the isotropic rapid-distortion result
 * (Wong & Hoult 1979, SAE 790357: turbulence of the unburned gas "depends only on its initial value and
 * the degree of compression"; same (2/3)(k/ρ)dρ/dt source term in the 0-D k–ε model of Sjerić et al.
 * 2014, Thermal Science 18(1), eq. 1 — fetched). It is what makes u′ rise during compression.
 *
 * Length scale (CALIBRATED DEPARTURE from P&H): L here is the DISSIPATION length L ≡ u′³/ε of the model
 * (P&H's ε = u′³/L). Poulos & Heywood set it to the chamber height h = 4V/(πB²) (Xu & Filipi 2020 eq. 32,
 * a closed-cycle study initialised at IVC). Driven through the intake stroke with c_β = 1 that choice
 * leaves far too little dissipation: the synthetic motored CFR cycle of turbulence.test.ts gives
 * u′(TDC) = 1.87 S̄p (f_L = 1), 0.98 S̄p (f_L = 0.5), 0.51 S̄p (f_L = 0.25) at 600 rpm, CR 7 (reviewer re-run)
 * against the measured 0.45–0.6 S̄p without induction swirl and 0.4–0.5 S̄p with it (Bopp, Vafidis &
 * Whitelaw 1986, SAE 860023, abstract — fetched; volume-averaged, three components, motored engine).
 * Measured in-cylinder integral scales near TDC are "limited by the clearance height, scaling typically by
 * about 10–15 %" vertically and "5–12 %" of the bore horizontally (Aleiferis & Behringer 2017, Fuel 189,
 * 238–259, abstract — fetched), and ε = C_ε u′³/L_int with C_ε ≈ 0.4–0.5 for high-Re forced isotropic
 * turbulence (Sreenivasan 1998, Phys. Fluids 10:528; McComb et al. 2010, arXiv:1002.2131), so
 *   L = f_L · min(h, B),  f_L = L_int/(C_ε h) ≈ 0.12/0.5 ≈ 0.25.
 * f_L is therefore a calibration constant supported by two independent measurements, not a P&H value.
 * With f_L = 0.25 and the standard c_β = 1 the model gives u′(TDC) = 0.42–0.59 S̄p for CR 10–5,
 * independent of speed and of the valve discharge coefficient (see turbulence.test.ts). Known limitation:
 * the measured ratio grows with speed (0.45 → 0.6 over 300–2000 rpm, Bopp et al.); a speed-independent
 * ratio is intrinsic to this inviscid 0-D model.
 * During combustion the unburned-gas eddies conserve angular momentum: L = L_s (ρ_s/ρ)^{1/3}
 * ({@link angularMomentumLengthScale}; Tabaczynski, Ferguson & Radhakrishnan 1977, SAE 770647; used by
 * Poulos & Heywood during combustion per Xu & Filipi 2020); take the minimum with the geometric value.
 *
 * Validation targets: isotropic decay u′(t) = u′₀/(1 + u′₀t/(3L)), RDT u′ ∝ ρ^{1/3}, and u′(TDC) =
 * 0.45–0.6 × mean piston speed for a motored engine without induction swirl (Bopp, Vafidis & Whitelaw
 * 1986, SAE 860023, measured — fetched abstract).
 *
 * Zero-turbulence seed: P ∝ √k is not Lipschitz at k = 0, so k ≡ 0 is a spurious solution (turbulence
 * could never start from rest, or after an explicit integrator undershoots k below 0). The production
 * therefore uses k_eff = max(k, K_SEED_FRACTION·K), which selects the physical (maximal) solution and
 * changes nothing once u′ ≳ 1e-5 U.
 *
 * All functions are allocation-free; state/inputs/outputs are plain objects the caller reuses.
 */

/** Turbulence state (extensive, SI). */
export interface TurbulenceState {
  /** Mean-flow kinetic energy K = ½ m U², J. */
  K: number;
  /** Turbulent kinetic energy k = (3/2) m u′², J. */
  k: number;
  /** Swirl angular momentum about the cylinder axis, kg m²/s (0 when swirl is not modelled). */
  swirl: number;
}

/** Time derivatives written by {@link turbulenceDerivatives}. */
export interface TurbulenceRates {
  /** dK/dt, W. */
  dK: number;
  /** dk/dt, W. */
  dk: number;
  /** d(swirl)/dt, N m. */
  dSwirl: number;
  /** Production P (mean flow → turbulence), W (diagnostic). */
  production: number;
  /** Viscous dissipation m ε, W (diagnostic). */
  dissipation: number;
  /** Swirl wall-friction torque, N m (diagnostic, ≥ 0 opposes the swirl). */
  swirlFrictionTorque: number;
}

/** Instantaneous inputs to {@link turbulenceDerivatives} (reuse one object; all fields required). */
export interface TurbulenceInputs {
  /** Cylinder (or zone) gas mass, kg. */
  m: number;
  /** Total mass inflow rate (intake + backflow through any valve), kg/s, ≥ 0. */
  mDotIn: number;
  /** Total mass outflow rate, kg/s, ≥ 0. */
  mDotOut: number;
  /** Jet speed of the inflow, m/s (see {@link intakeJetVelocity}). */
  vIn: number;
  /** Dissipation length scale L ≡ u′³/ε, m (see {@link integralLengthScale}). */
  L: number;
  /**
   * Parcel compression rate (1/ρ) Dρ/Dt of the gas carrying the turbulence, 1/s. Closed homogeneous
   * charge: −(dV/dt)/V. Well-mixed cylinder with open valves: (1/ρ) dρ/dt = ṁ_net/m − (dV/dt)/V.
   * Unburned zone during combustion (isentropic): (1/(γ_u p)) dp/dt. 0 disables RDT.
   */
  dlnRhoDt: number;
  /** Angular-momentum flux of the inflow about the cylinder axis, N m (see {@link swirlTorqueFromJet}). */
  swirlTorqueIn: number;
  /** Gas density, kg/m³ (swirl wall friction). */
  rho: number;
  /** Gas dynamic viscosity, Pa s (swirl wall friction). */
  mu: number;
  /** Cylinder bore, m. */
  bore: number;
  /** Instantaneous clearance height, m (swirl wall friction on the liner). */
  h: number;
}

/** Model constants. */
export interface TurbulenceParams {
  /**
   * Production multiplier c_β of P = 0.3307 c_β (K/L)√(k/m). 1 = standard k–ε eddy viscosity
   * (C_μ = 0.09); Poulos & Heywood (1983) treat it as an adjustable constant.
   */
  cBeta: number;
  /** Rapid-distortion coefficient on d ln ρ/dt: 2/3 for isotropic compression; 0 disables. */
  cRdt: number;
  /** Apply the same rapid-distortion scaling to the mean-flow energy K (default false: P&H has none). */
  rdtMeanFlow: boolean;
  /** Fraction of the swirl wall-friction power fed to k (0 = all dissipated in the wall layer). */
  swirlFrictionToTurbulence: number;
}

/** Default constants (see {@link TurbulenceParams}). */
export const DEFAULT_TURBULENCE_PARAMS: Readonly<TurbulenceParams> = Object.freeze({
  cBeta: 1,
  cRdt: 2 / 3,
  rdtMeanFlow: false,
  swirlFrictionToTurbulence: 0,
});

/**
 * P = C_P c_β (K/L)√(k/m): C_P = 4.5 C_μ √(2/3) = 0.33068 with C_μ = 0.09 (Launder & Spalding 1974).
 * Verified against the published 0.3307 of Poulos & Heywood (1983) as reproduced in Xu & Filipi 2020,
 * Front. Mech. Eng. 6:46, eq. 37 (fetched JATS XML). The original SAE 830334 text was not accessible; the
 * task brief's "0.3875 β" is not supported by any source found.
 */
export const KK_PRODUCTION_CONSTANT = 4.5 * 0.09 * Math.sqrt(2 / 3);

/**
 * Mean skin-friction coefficient of a turbulent flat plate, c_f = 0.074 Re^{-1/5}
 * (Prandtl 1927 / Schlichting, "Boundary-Layer Theory", 5×10⁵ < Re < 10⁷; CFD-Online lists 0.0725 —
 * checked via search summaries only).
 * UNVERIFIED model choice: applied to the swirl wall friction with Re = ρωR²/μ, which is ~1e4–1e5 in a
 * CFR-size cylinder — below the correlation's validity range; the swirl sub-model is optional and
 * not validated against engine data.
 */
const CF_COEFF = 0.074;

/**
 * Seed of the production term: P uses k_eff = max(k, K_SEED_FRACTION · K) (see the module comment).
 * Numerical device, not a physical constant: it only matters while u′ < ~1e-5 × U.
 */
export const K_SEED_FRACTION = 1e-10;

/** Allocate a zeroed {@link TurbulenceState}. */
export const newTurbulenceState = (): TurbulenceState => ({ K: 0, k: 0, swirl: 0 });
/** Allocate a zeroed {@link TurbulenceRates}. */
export const newTurbulenceRates = (): TurbulenceRates => ({
  dK: 0,
  dk: 0,
  dSwirl: 0,
  production: 0,
  dissipation: 0,
  swirlFrictionTorque: 0,
});
/** Allocate a {@link TurbulenceInputs} with neutral values (no flow, no compression, no swirl). */
export const newTurbulenceInputs = (): TurbulenceInputs => ({
  m: 0,
  mDotIn: 0,
  mDotOut: 0,
  vIn: 0,
  L: 1,
  dlnRhoDt: 0,
  swirlTorqueIn: 0,
  rho: 1,
  mu: 1.8e-5,
  bore: 0.1,
  h: 0.1,
});

/** Turbulence intensity u′ = √(2k/(3m)), m/s, from the extensive TKE k (J) and mass m (kg). */
export function turbulenceIntensity(k: number, m: number): number {
  return k > 0 && m > 0 ? Math.sqrt((2 * k) / (3 * m)) : 0;
}

/** Mean-flow speed U = √(2K/m), m/s. */
export function meanFlowVelocity(K: number, m: number): number {
  return K > 0 && m > 0 ? Math.sqrt((2 * K) / m) : 0;
}

/** Moment of inertia of the charge in solid-body rotation about the cylinder axis, I = m B²/8, kg m². */
export function swirlInertia(m: number, bore: number): number {
  return (m * bore * bore) / 8;
}

/** Swirl angular velocity ω_s = H/I, rad/s. */
export function swirlAngularVelocity(swirl: number, m: number, bore: number): number {
  return m > 0 ? swirl / swirlInertia(m, bore) : 0;
}

/**
 * f_L in L = f_L · min(h, B): L_int/h ≈ 0.10–0.15 (Aleiferis & Behringer 2017, Fuel 189:238–259, abstract
 * fetched) divided by the dissipation coefficient C_ε = εL_int/u′³ ≈ 0.4–0.5.
 * UNVERIFIED: C_ε ≈ 0.4–0.5 comes from search-result summaries of Sreenivasan 1998 (Phys. Fluids 10:528)
 * and McComb et al. 2010; the primary texts could not be read (scanned PDF). C_ε is not universal
 * (≈ 0.5 shear flows, 0.7–1 grid turbulence), which is why f_L is exposed as a parameter.
 */
export const DEFAULT_LENGTH_SCALE_FRACTION = 0.25;

/**
 * Dissipation coefficient C_ε = ε L_int/u′³ that links the model's dissipation length L (ε = u′³/L) to
 * the INTEGRAL scale, L_int = C_ε L — the value 0.5 of the f_L derivation above (L_int/h ≈ 0.12,
 * f_L = 0.25). UNVERIFIED (see DEFAULT_LENGTH_SCALE_FRACTION): C_ε ≈ 0.4–0.5 from search summaries.
 */
export const DISSIPATION_COEFFICIENT = 0.5;

/**
 * Geometric dissipation length scale L = f_L · min(h, B), m (see the module comment).
 * @param h instantaneous clearance height, m
 * @param bore bore, m
 * @param fraction f_L (default {@link DEFAULT_LENGTH_SCALE_FRACTION})
 */
export function integralLengthScale(h: number, bore: number, fraction = DEFAULT_LENGTH_SCALE_FRACTION): number {
  return fraction * (h < bore ? h : bore);
}

/**
 * Length scale of eddies of fixed mass conserving angular momentum under isotropic compression from a
 * reference state (e.g. spark): L = L_ref (ρ_ref/ρ)^{1/3}, m — Poulos & Heywood 1983 during combustion,
 * as reproduced in Xu & Filipi 2020 eq. 39 (fetched). During combustion use min(this,
 * {@link integralLengthScale}). The exponent is kinematic (an eddy of fixed mass compressed isotropically
 * scales as ρ^{-1/3}).
 * UNVERIFIED: the further attribution of the idea to Tabaczynski, Ferguson & Radhakrishnan 1977
 * (SAE 770647) — paper not accessible.
 * @param LRef length scale at the reference state, m
 * @param rhoRef density at the reference state, kg/m³
 * @param rho current density (unburned gas), kg/m³
 */
export function angularMomentumLengthScale(LRef: number, rhoRef: number, rho: number): number {
  return LRef * Math.cbrt(rhoRef / rho);
}

/**
 * Rapid-distortion (angular-momentum-conservation) scaling of turbulence of fixed mass under isotropic
 * compression from density ρ₀ to ρ: u′ = u′₀ (ρ/ρ₀)^{1/3}, L = L₀ (ρ₀/ρ)^{1/3} (Wong & Hoult 1979;
 * Tabaczynski et al. 1977). Returns u′ (m/s); the matching L is {@link angularMomentumLengthScale}.
 */
export function rapidDistortionIntensity(uPrime0: number, rho0: number, rho: number): number {
  return uPrime0 * Math.cbrt(rho / rho0);
}

/**
 * Jet speed of a valve inflow, v = ṁ / (ρ A_eff), m/s, with A_eff = C_D × reference (curtain) area.
 * @param mDot mass flow rate, kg/s (sign ignored)
 * @param rho density of the flowing gas at the valve throat, kg/m³
 * @param effectiveArea C_D·A_ref, m²
 */
export function intakeJetVelocity(mDot: number, rho: number, effectiveArea: number): number {
  if (!(rho > 0) || !(effectiveArea > 0)) return 0;
  return Math.abs(mDot) / (rho * effectiveArea);
}

/**
 * Angular-momentum flux (N m) about the cylinder (+y) axis carried by a valve jet leaving at horizontal
 * position (x, z) with horizontal velocity v_h (cos ψ, sin ψ) in the x–z plane: the y-component of
 * ṁ (r × v) = ṁ (z v_x − x v_z).
 * @param mDot jet mass flow, kg/s (≥ 0)
 * @param vHorizontal horizontal jet speed, m/s
 * @param x jet-exit x (cylinder frame), m
 * @param z jet-exit z (cylinder frame), m
 * @param psi jet direction angle from +x toward +z, rad
 */
export function swirlTorqueFromJet(mDot: number, vHorizontal: number, x: number, z: number, psi: number): number {
  const vx = vHorizontal * Math.cos(psi);
  const vz = vHorizontal * Math.sin(psi);
  return mDot * (z * vx - x * vz);
}

/**
 * Right-hand side of the K–k(–swirl) model. Pure and allocation-free: reads `state`, `inputs`,
 * `params`, writes `out` (and returns it). Signature per DESIGN.md §Module contracts
 * (`turbulenceDerivatives(state, inputs, out)`), with the model constants as an optional 4th argument.
 * @param state K (J), k (J), swirl (kg m²/s)
 * @param inputs instantaneous gas/flow inputs (SI)
 * @param out rates to fill: dK, dk (W), dSwirl (N m) + diagnostics
 * @param params model constants (default {@link DEFAULT_TURBULENCE_PARAMS})
 */
export function turbulenceDerivatives(
  state: TurbulenceState,
  inputs: TurbulenceInputs,
  out: TurbulenceRates,
  params: Readonly<TurbulenceParams> = DEFAULT_TURBULENCE_PARAMS,
): TurbulenceRates {
  const m = inputs.m;
  const K = state.K > 0 ? state.K : 0;
  const k = state.k > 0 ? state.k : 0;
  const L = inputs.L;
  if (!(m > 0) || !(L > 0)) {
    out.dK = 0;
    out.dk = 0;
    out.dSwirl = 0;
    out.production = 0;
    out.dissipation = 0;
    out.swirlFrictionTorque = 0;
    return out;
  }
  const km = k / m; // specific TKE, J/kg
  const uP = Math.sqrt((2 / 3) * km); // u′
  // production with the k = 0 seed (see K_SEED_FRACTION): P ∝ √k is non-Lipschitz at 0
  const kSeed = K_SEED_FRACTION * K;
  const P = KK_PRODUCTION_CONSTANT * params.cBeta * (K / L) * Math.sqrt((k > kSeed ? k : kSeed) / m);
  const diss = (m * uP * uP * uP) / L; // m ε
  const outFrac = inputs.mDotOut > 0 ? inputs.mDotOut / m : 0;
  const mIn = inputs.mDotIn > 0 ? inputs.mDotIn : 0;
  let dK = 0.5 * mIn * inputs.vIn * inputs.vIn - P - K * outFrac;
  // linear terms act on the signed state so an explicit-integrator undershoot below 0 relaxes back
  let dk = P - diss - state.k * outFrac;
  if (params.cRdt !== 0) {
    dk += params.cRdt * k * inputs.dlnRhoDt;
    if (params.rdtMeanFlow) dK += params.cRdt * K * inputs.dlnRhoDt;
  }
  // ---- swirl: solid-body rotation, turbulent wall friction on head, piston and liner ----
  let dSwirl = inputs.swirlTorqueIn - state.swirl * outFrac;
  let Tf = 0;
  const H = state.swirl;
  if (H !== 0 && inputs.rho > 0 && inputs.mu > 0) {
    const R = 0.5 * inputs.bore;
    const omega = H / ((m * inputs.bore * inputs.bore) / 8);
    const aw = Math.abs(omega);
    const Re = (inputs.rho * aw * R * R) / inputs.mu; // tip speed × radius
    const cf = CF_COEFF * Math.pow(Re > 1 ? Re : 1, -0.2);
    // torque = ∫ ½ρ c_f (ωr)² r dA: two discs π ρ c_f ω² R⁵/5 each + liner π ρ c_f ω² R⁴ h
    Tf = Math.PI * inputs.rho * cf * omega * aw * R * R * R * R * (0.4 * R + inputs.h);
    dSwirl -= Tf;
    if (params.swirlFrictionToTurbulence > 0) dk += params.swirlFrictionToTurbulence * Tf * omega;
  }
  out.dK = dK;
  out.dk = dk;
  out.dSwirl = dSwirl;
  out.production = P;
  out.dissipation = diss;
  out.swirlFrictionTorque = Tf;
  return out;
}

/**
 * Initialise the turbulence state (e.g. at IVC for a closed-cycle-only run, or at the start of the
 * simulation), from an intensity u′, a mean-flow speed U and a swirl angular velocity ω_s.
 * When the gas-exchange strokes are simulated the K–k equations build their own IVC state and this
 * only seeds the first cycle (u′ = 0 is allowed: production is seeded, see {@link K_SEED_FRACTION}).
 * @param m trapped mass, kg
 * @param uPrime turbulence intensity, m/s
 * @param meanFlowSpeed mean-flow speed U, m/s
 * @param swirlOmega swirl angular velocity, rad/s
 * @param bore bore, m
 * @param out state to fill
 */
export function initTurbulenceState(
  m: number,
  uPrime: number,
  meanFlowSpeed: number,
  swirlOmega: number,
  bore: number,
  out: TurbulenceState,
): TurbulenceState {
  out.k = 1.5 * m * uPrime * uPrime;
  out.K = 0.5 * m * meanFlowSpeed * meanFlowSpeed;
  out.swirl = swirlInertia(m, bore) * swirlOmega;
  return out;
}
