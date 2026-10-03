/**
 * Spark-kernel growth from the breakdown channel to a self-sustaining flame, after Herweg &
 * Maly (1992), "A fundamental model for flame kernel formation in S.I. engines", SAE 922243.
 * SI units. The equations are as reproduced (fetched) in:
 *   [HM-R] "Choice of tuning parameters on 3D IC engine simulations using G-equation",
 *          SAE 18PFL-1037 (OSTI 1435273), eqs. 7–11 (DPIK/ANSYS Forte implementation of HM);
 *   [G19]  Gianetti, Sforza, Lucchini et al., AIP Conf. Proc. 2191:020087 (2019), eqs. 5–9, 14.
 *
 *  Kernel growth (sphere of radius r_k, density ρ_k, in unburned gas of density ρ_u):
 *     dr_k/dt = (ρ_u/ρ_k)(S_plasma + S_T,k)                                     [HM-R eq. 7]
 *     S_plasma = η Q̇_spk / (4π r_k² [ρ_u(u_k − h_u) + p ρ_u/ρ_k])               [HM-R eq. 8]
 *              = Q̇_net / (4π r_k² ρ_u (h_k − h_u))      (identity: u_k + p/ρ_k = h_k)
 *   with h_k − h_u the enthalpy needed to heat fresh mixture from T_u to the kernel
 *   temperature (same composition, thermo module). We integrate the kernel VOLUME:
 *     dV_k/dt = Q̇_net/(ρ_k Δh) + A_k (ρ_u/ρ_k) S_T,k,       Q̇_net = η Q̇_spk − Q̇_electrodes
 *   (the first term is independent of r_k, which keeps the tiny initial kernel non-stiff).
 *
 *  Kernel burning velocity                                                     [HM-R eq. 9]
 *     S_T,k/S_L = I₀ + √I₀ · √(u′/(u′+S_L)) · √(1 − e^{−r_k/l_I}) · √(1 − e^{−t(u′+S_L)/l_I}) · (u′/S_L)^{5/6}
 *   — only eddies smaller than the kernel wrinkle it (r_k/l_I factor) and the turbulence
 *   spectrum needs about one eddy turnover to act (time factor), t from the spark.
 *
 *  Stretch factor (HM-R eq. 10/11: strain + curvature, times the thermo-diffusive factor
 *  [1/Le + (Le−1)/Le · T_a/T_ad], which is Clavin's asymptotic Markstein number). We use the
 *  caller's Markstein length L_M (unburned-gas, w.r.t. total stretch) directly, i.e. the
 *  linear stretch law S_n = S_L − L_M κ with κ = κ_strain + κ_curv:
 *     I₀ = 1 − (L_M/l_F) [ (l_F/(15 l_I))^{1/2} (u′/S_L)^{3/2} + 2 (l_F/r_k)(ρ_u/ρ_k) ]
 *   where the first bracket term is the turbulent-strain Karlovitz number and the second the
 *   curvature Karlovitz number 2 σ S_L/r_k · l_F/S_L of the propagating spherical kernel.
 *   I₀ is clipped to [0, maxStretchFactor]; I₀ ≤ 0 means the kernel is below its critical
 *   radius r_c = 2σL_M/(1 − Ma·Ka_t) and cannot propagate without spark energy.
 *   Convention: flameThickness l_F = ν_u/S_L (Abdel-Gayed & Bradley), so the strain term is
 *   L_M·(u′/λ)/S_L with the Taylor microscale λ = √(15 ν l_I/u′). One ν is used for both
 *   strain measures (this term and the AGB R_L below): the caller's kinematicViscosity,
 *   else l_F·S_L, else thermo transport (reviewer change: the two could previously
 *   disagree when l_F ≠ ν/S_L, e.g. the combustion module's α_u/S_L thickness).
 *   l_I → 0 is clamped to MIN_INTEGRAL_SCALE (1 µm) — the quiescent limit u′ = 0 previously
 *   produced NaN (0/0 in the size/time factors).
 *
 *  Kernel state: the kernel gas is at the adiabatic flame state of the mixture
 *  (ρ_k = ρ_u/σ, σ = caller's expansion ratio; T_k = σ T_u M_b/M_u with M_b the molar mass
 *  of complete-combustion products) — the usual assumption of HM-type implementations
 *  (UNVERIFIED which kernel temperature HM used). For a non-combustible charge (S_L = 0 or
 *  σ ≤ 1) the kernel is an inert hot-gas sphere at `inertKernelTemperature` (default 3000 K,
 *  the peak glow-discharge gas temperature quoted by Heywood 1988 p. 429 — UNVERIFIED as a
 *  kernel-mean value). Note HM's S_plasma assumes the plasma-heated gas does not add its
 *  own chemical energy to the kernel (it is accounted for in S_T,k).
 *
 *  Initial kernel: the impulsive breakdown + capacitive-arc energy E₀ (to gas) creates a
 *  kernel of volume E₀/(ρ_k Δh) (the time integral of the S_plasma term over the ns–µs
 *  pulse). Shock-wave energy is regained within the kernel (Heywood 1988 p. 428: "most of
 *  this is regained since spherical blast waves transfer most of their energy to the gas
 *  within a small (~2 mm diameter) sphere").
 *
 *  Electrode heat loss (Shaffer, Luna, Wang, Egolfopoulos & Askari 2023, J. Phys. D 56:225501,
 *  eq. 8, after Keck 1981): energy in the thermal boundary layer on the wetted electrode area
 *     E_cond = γ/(γ−1) · p · δ · A_cond,   δ = (2/π)(1 − T_e/T_k)√(α t)
 *  so over a sub-step ΔE = γ/(γ−1) p (2/π)(1 − T_e/T_k)√α · A_cond(r_k) · (√(t+h) − √t).
 *  (Reviewer check: the 2/π prefactor is as printed in Shaffer et al. eq. 8 — re-extracted
 *  from the NSF-PAR full text; the constant-property semi-infinite conduction solution
 *  would give 2/√π, i.e. 1.77× more loss. Kept as cited.)
 *  A_cond is computed from the plug geometry (centre-electrode end face + flank, ground-strap
 *  face). Radiation is neglected (significant only above 6000 K; Shaffer et al. 2023 §3).
 *  The internal sub-step also limits the conduction loss to maxRelativeVolumeChange of the
 *  kernel volume (the loss rate ∝ A/√t is singular at t = 0).
 *
 *  Quench / misfire:
 *   - turbulent flame quench (Abdel-Gayed, Bradley & Lung 1989, Combust. Flame 76:213,
 *     fetched): K = 0.157 (u′/S_L)² R_L^{−1/2} (their eq. 11), complete quench when
 *     K·Le > 1.5 (R_L > 300) or K·R_L^{−1/2} > 0.079 (R_L < 300); R_L = u′ l_I/ν_u;
 *   - after the discharge has ended, a kernel that cannot propagate (I₀ ≤ 0 or quenched) and does
 *     not resume propagating within nonPropagatingSurvival × min(l_I/u′, r²/α_k) (fixer round 2;
 *     round 1: at once), or that shrinks below half its largest radius, is extinguished (misfire);
 *   - a kernel not handed off within `maxKernelTime` is declared a misfire.
 *  Hand-off to the turbulent entrainment model when r_k ≥ max(minHandoffRadius, C·l_I) and
 *  I₀ > 0. C = 1 by default: ANSYS Fluent's Herweg–Maly spark model grows the kernel "until
 *  the length scale is reached" — the local turbulent length scale, beyond which "the flame
 *  speed is affected by all the turbulent scales present" — and then switches the spark
 *  flame-speed model off (ANSYS Fluent Theory Guide 2025 R1, §10.1.2 "Spark Model Theory",
 *  fetched). ANSYS Forte instead uses C_m1 = 2.0 ([HM-R] eq. 13, "kernel-flame-to-G-equation
 *  switch constant, always equal to 2.0"). The hand-off time is located inside the internal
 *  sub-step (linear in V), so it does not depend on the caller's step.
 *  UNVERIFIED model choices: minHandoffRadius = 1 mm, maxStretchFactor = 3, maxKernelTime.
 */

import type { SparkPlugSpec } from '../core/engine-spec';
import { NS } from '../core/species';
import { completeCombustionProducts } from '../thermo/fuels';
import { mixGamma, mixHMass, mixMolarMass } from '../thermo/mixture';
import { mixKinematicViscosity, mixThermalDiffusivity } from '../thermo/transport';

/** Kernel life-cycle stage. */
export type KernelStage = 'none' | 'kernel' | 'handoff' | 'quenched';

/** Gas state seen by the kernel (unburned mixture around the gap). SI units. */
export interface KernelGasState {
  /** Pressure, Pa. */
  p: number;
  /** Unburned-gas temperature, K. */
  Tu: number;
  /** Unburned-gas density, kg/m³. */
  rhoU: number;
  /** Unburned mole fractions, Float64Array(NS). */
  X: Float64Array;
  /** Unstretched laminar burning velocity, m/s (0 = non-combustible). */
  SL: number;
  /** Markstein length (unburned side, total stretch), m. */
  marksteinLength: number;
  /** Laminar flame thickness l_F = ν_u/S_L, m. */
  flameThickness: number;
  /** Turbulence intensity u′, m/s. */
  uPrime: number;
  /** Integral length scale l_I, m (eddy size: HM size/time factors, hand-off radius, AGB R_L). */
  integralScale: number;
  /**
   * Dissipation length L_ε = u′³/ε, m, of the caller's turbulence model — the length in the Taylor
   * microscale λ = √(15 ν L_ε/u′) of the HM strain Karlovitz number. Default: integralScale (round-1
   * behaviour, which then had to receive L_ε and made the eddy-size uses ≈ 2× too large).
   */
  dissipationLength?: number;
  /** Expansion ratio σ = ρ_u/ρ_b at the adiabatic flame state. */
  expansionRatio: number;
  /** Lewis number of the deficient reactant (quench criterion; default 1). */
  lewisNumber?: number;
  /** Unburned kinematic viscosity ν_u, m²/s (default: from thermo transport). */
  kinematicViscosity?: number;
  /** Electrode surface temperature, K (default: options.electrodeTemperature). */
  electrodeTemperature?: number;
}

/** Kernel model constants (see module doc). */
export interface SparkKernelOptions {
  /** Hand-off radius as a multiple of l_I (default 1). */
  handoffIntegralScaleMultiple: number;
  /** Minimum hand-off radius, m (default 1 mm; UNVERIFIED). */
  minHandoffRadius: number;
  /** Kernel temperature for a non-combustible charge, K (default 3000; UNVERIFIED). */
  inertKernelTemperature: number;
  /** Kernel older than this without hand-off → misfire, s (default 20 ms; UNVERIFIED). */
  maxKernelTime: number;
  /** Upper clip of the stretch factor I₀ (negative Markstein lengths), default 3 (UNVERIFIED). */
  maxStretchFactor: number;
  /** Default electrode temperature, K (default 700; UNVERIFIED typical plug-tip value). */
  electrodeTemperature: number;
  /** Include electrode conduction losses (default true). */
  electrodeHeatLoss: boolean;
  /** Apply the Abdel-Gayed–Bradley–Lung turbulent quench criterion (default true). */
  turbulentQuench: boolean;
  /** Max relative volume change per internal sub-step (default 0.05). */
  maxRelativeVolumeChange: number;
  /**
   * Survival of a NON-PROPAGATING kernel after the discharge (I₀ ≤ 0 or AGB-quenched), in integral
   * eddy-turnover times l_I/u′ (see step): 0 = extinguished at once (round-1 behaviour).
   */
  nonPropagatingSurvival: number;
}

/**
 * Smallest integral length scale the kernel model uses, m (1 µm, below any laminar flame
 * thickness): l_I = 0 (quiescent input) is clamped to it so the Herweg–Maly size/time
 * factors and strain rate stay finite; for u′ = 0 the results do not depend on it.
 */
export const MIN_INTEGRAL_SCALE = 1e-6;

/** Default kernel options. */
export const DEFAULT_KERNEL_OPTIONS: Readonly<SparkKernelOptions> = Object.freeze({
  handoffIntegralScaleMultiple: 1,
  minHandoffRadius: 1e-3,
  inertKernelTemperature: 3000,
  maxKernelTime: 20e-3,
  maxStretchFactor: 3,
  electrodeTemperature: 700,
  electrodeHeatLoss: true,
  turbulentQuench: true,
  maxRelativeVolumeChange: 0.05,
  nonPropagatingSurvival: 1,
});

/**
 * Karlovitz number of Abdel-Gayed, Bradley & Lung (1989) eq. 11: K = 0.157 (u′/S_L)² R_L^−½,
 * R_L = u′ l_I/ν (u′ rms, m/s; S_L m/s; l_I m; ν m²/s). 0 for u′ ≤ 0 (quiescent limit;
 * K ∝ u′^{3/2} → 0), +∞ for S_L ≤ 0 or l_I ≤ 0 with u′ > 0.
 */
export function abdelGayedKarlovitz(uPrime: number, SL: number, integralScale: number, nu: number): number {
  if (!(uPrime > 0)) return 0;
  if (!(SL > 0) || !(integralScale > 0)) return Infinity;
  const RL = (uPrime * integralScale) / nu;
  const r = uPrime / SL;
  return (0.157 * r * r) / Math.sqrt(RL);
}

/**
 * Complete turbulent flame quench criterion of Abdel-Gayed, Bradley & Lung (1989):
 * K·Le > 1.5 for R_L > 300, K·R_L^−½ > 0.079 for R_L < 300.
 */
export function turbulentQuenchReached(uPrime: number, SL: number, integralScale: number, nu: number, Le: number): boolean {
  if (uPrime <= 0) return false;
  if (SL <= 0) return true;
  const RL = (uPrime * integralScale) / nu;
  const K = abdelGayedKarlovitz(uPrime, SL, integralScale, nu);
  return RL > 300 ? K * Le > 1.5 : K / Math.sqrt(RL) > 0.079;
}

/**
 * Herweg–Maly stretch factor I₀ (unclipped) with a Markstein length:
 *   I₀ = 1 − (L_M/l_F)[(l_F/(15 l_I))^½ (u′/S_L)^{3/2} + 2 (l_F/r)σ_k]
 *      = 1 − L_M [ (u′/S_L)^{3/2}/√(15 l_F l_I) + 2σ_k/r ]
 * (the second, division-free form is evaluated: finite for L_M = 0 or l_F → 0 without
 * strain). r, l_F, L_M, l_I in m; S_L, u′ in m/s. No turbulent strain for u′ ≤ 0; for
 * u′ > 0 with l_I ≤ 0 or l_F ≤ 0 the strain rate is infinite.
 */
export function herwegMalyStretchFactor(
  r: number,
  SL: number,
  lF: number,
  LM: number,
  uPrime: number,
  lI: number,
  sigmaK: number,
): number {
  if (LM === 0) return 1;
  const strain = uPrime > 0 ? Math.pow(uPrime / SL, 1.5) / Math.sqrt(15 * lF * lI) : 0;
  return 1 - LM * (strain + (2 * sigmaK) / r);
}

/**
 * Herweg–Maly kernel burning velocity S_T,k, m/s ([HM-R] eq. 9) for a clipped I₀ ≥ 0.
 * r, l_I in m; t = time since the spark, s. For l_I ≤ 0 the size and time factors take
 * their limit 1 (every eddy is smaller than the kernel).
 */
export function herwegMalyKernelSpeed(I0: number, r: number, t: number, SL: number, uPrime: number, lI: number): number {
  if (!(SL > 0) || !(I0 > 0)) return 0;
  const up = uPrime > 0 ? uPrime : 0;
  if (up === 0) return SL * I0;
  const fu = Math.sqrt(up / (up + SL));
  const fr = lI > 0 ? Math.sqrt(1 - Math.exp(-r / lI)) : 1;
  const ft = lI > 0 ? Math.sqrt(1 - Math.exp(-(t * (up + SL)) / lI)) : 1;
  return SL * (I0 + Math.sqrt(I0) * fu * fr * ft * Math.pow(up / SL, 5 / 6));
}

/** Area of a circle of radius rho inside the strip |x| ≤ halfW, m². */
function circleStripArea(rho: number, halfW: number): number {
  const x = Math.min(halfW, rho);
  return 2 * (x * Math.sqrt(Math.max(0, rho * rho - x * x)) + rho * rho * Math.asin(x / rho));
}

/**
 * Electrode surface area wetted by a kernel sphere of radius r centred in the gap, m²:
 * centre-electrode end face + flank (cylinder of diameter dc whose face is gap/2 from the
 * centre) and the ground-strap face (strip of width w, gap/2 on the other side).
 */
export function electrodeContactArea(r: number, gap: number, dc: number, w: number): number {
  const z0 = 0.5 * gap;
  if (r <= z0) return 0;
  const rho = Math.sqrt(r * r - z0 * z0);
  const a = 0.5 * dc;
  // centre electrode: end face + lateral flank
  const face = Math.PI * Math.min(rho, a) ** 2;
  const flank = rho > a ? 2 * Math.PI * a * (Math.sqrt(r * r - a * a) - z0) : 0;
  // ground strap face
  const strap = circleStripArea(rho, 0.5 * w);
  return face + flank + strap;
}

/**
 * Spherical spark kernel (Herweg & Maly 1992). Create once; `ignite` at breakdown, then
 * `step` every engine time step. Allocation-free in `step`.
 */
export class SparkKernel {
  readonly opts: SparkKernelOptions;
  readonly gap: number;
  readonly centreElectrodeDiameter: number;
  readonly groundElectrodeWidth: number;

  // ---- state / outputs ----
  stage: KernelStage = 'none';
  /** Kernel volume, m³. */
  volume = 0;
  /** Kernel radius, m. */
  radius = 0;
  /** Largest radius reached, m. */
  maxRadius = 0;
  /** Time since the first breakdown, s. */
  age = 0;
  /** Kernel gas temperature, K. */
  temperature = 0;
  /** Kernel gas density, kg/m³. */
  density = 0;
  /** Kernel (burned) mass ρ_k V_k, kg. */
  burnedMass = 0;
  /** Stretch factor I₀ (clipped). */
  stretchFactor = 1;
  /** Kernel burning velocity S_T,k, m/s. */
  turbulentSpeed = 0;
  /** Plasma velocity S_plasma (relative to unburned gas), m/s. */
  plasmaSpeed = 0;
  /** dr_k/dt at the end of the last step, m/s. */
  growthRate = 0;
  /** Strain / curvature Karlovitz numbers of the stretch factor. */
  karlovitzStrain = 0;
  karlovitzCurvature = 0;
  /** Abdel-Gayed–Bradley Karlovitz number K. */
  karlovitzAGB = 0;
  /** True if the AGB complete-quench criterion is met. */
  turbulentlyQuenched = false;
  /** Hand-off radius at the current turbulence, m. */
  handoffRadius = 0;
  /** Time of hand-off since breakdown, s (NaN if none). */
  handoffTime = NaN;
  /** Energy deposited in the kernel by the discharge, J. */
  energyDeposited = 0;
  /** Electrode conduction loss, J. */
  energyElectrodeLoss = 0;
  /** Why the kernel was quenched ('' if not). */
  quenchReason = '';
  /** Time since breakdown at which the kernel stopped propagating after the discharge (NaN: propagating). */
  blockedSince = NaN;

  /** True once the kernel was extinguished (misfire). */
  get misfire(): boolean {
    return this.stage === 'quenched';
  }

  // ---- per-step cached properties ----
  private sigmaK = 1;
  private dh = 1;
  private condCoeff = 0;
  /** Thermal diffusivity of the kernel gas, m²/s. */
  private alphaK = 0;
  /** Integral scale used by the kernel model (≥ MIN_INTEGRAL_SCALE), m. */
  private lI = MIN_INTEGRAL_SCALE;
  /** Dissipation length for the strain Karlovitz number (≥ MIN_INTEGRAL_SCALE), m. */
  private lEps = MIN_INTEGRAL_SCALE;
  /** Flame thickness l_F = ν_u/S_L consistent with the ν used by the AGB criterion, m. */
  private lF = 0;
  private readonly Xb = new Float64Array(NS);

  constructor(plug: SparkPlugSpec, opts?: Partial<SparkKernelOptions>) {
    this.opts = { ...DEFAULT_KERNEL_OPTIONS, ...opts };
    this.gap = plug.gap;
    this.centreElectrodeDiameter = plug.centerElectrodeDiameter;
    this.groundElectrodeWidth = plug.groundElectrodeWidth;
  }

  /** Forget the kernel (new cycle). */
  reset(): void {
    this.stage = 'none';
    this.volume = 0;
    this.radius = 0;
    this.maxRadius = 0;
    this.age = 0;
    this.temperature = 0;
    this.density = 0;
    this.burnedMass = 0;
    this.stretchFactor = 1;
    this.turbulentSpeed = 0;
    this.plasmaSpeed = 0;
    this.growthRate = 0;
    this.karlovitzStrain = 0;
    this.karlovitzCurvature = 0;
    this.karlovitzAGB = 0;
    this.turbulentlyQuenched = false;
    this.handoffRadius = 0;
    this.handoffTime = NaN;
    this.energyDeposited = 0;
    this.energyElectrodeLoss = 0;
    this.quenchReason = '';
    this.blockedSince = NaN;
  }

  /** Update the kernel-gas properties that depend only on the gas state (once per step). */
  private prepare(gas: KernelGasState): void {
    const o = this.opts;
    const combustible = gas.SL > 0 && gas.expansionRatio > 1 + 1e-9;
    let Tk: number;
    let Xk: Float64Array;
    if (combustible) {
      completeCombustionProducts(gas.X, this.Xb);
      Xk = this.Xb;
      this.sigmaK = gas.expansionRatio;
      Tk = (this.sigmaK * gas.Tu * mixMolarMass(Xk)) / mixMolarMass(gas.X);
    } else {
      Xk = gas.X;
      Tk = Math.max(o.inertKernelTemperature, gas.Tu);
      this.sigmaK = Tk / gas.Tu;
    }
    this.temperature = Tk;
    this.density = gas.rhoU / this.sigmaK;
    this.dh = Math.max(1, mixHMass(gas.X, Tk) - mixHMass(gas.X, gas.Tu));
    const alpha = mixThermalDiffusivity(Xk, Tk, gas.p);
    this.alphaK = alpha;
    if (o.electrodeHeatLoss) {
      const Te = gas.electrodeTemperature ?? o.electrodeTemperature;
      const g = mixGamma(Xk, Tk);
      this.condCoeff = Math.max(0, (g / (g - 1)) * gas.p * (2 / Math.PI) * (1 - Te / Tk) * Math.sqrt(alpha));
    } else {
      this.condCoeff = 0;
    }
    // l_I → 0 (quiescent bomb, "no turbulence" input) would give 0/0 in the size/time factors
    // and the strain rate; clamp to a scale below any flame thickness (results for u′ = 0
    // are independent of it).
    const lI = gas.integralScale > MIN_INTEGRAL_SCALE ? gas.integralScale : MIN_INTEGRAL_SCALE;
    this.lI = lI;
    const le = gas.dissipationLength !== undefined && gas.dissipationLength > MIN_INTEGRAL_SCALE ? gas.dissipationLength : lI;
    this.lEps = le;
    this.handoffRadius = Math.max(o.minHandoffRadius, o.handoffIntegralScaleMultiple * lI);
    if (combustible) {
      // One kinematic viscosity for BOTH strain measures (HM Karlovitz via l_F = ν/S_L and the
      // AGB R_L): the caller's ν, else ν = l_F·S_L from the caller's l_F, else thermo transport.
      const nu =
        gas.kinematicViscosity !== undefined && gas.kinematicViscosity > 0
          ? gas.kinematicViscosity
          : gas.flameThickness > 0
            ? gas.flameThickness * gas.SL
            : mixKinematicViscosity(gas.X, gas.Tu, gas.p);
      this.lF = nu / gas.SL;
      if (gas.uPrime > 0) {
        this.karlovitzAGB = abdelGayedKarlovitz(gas.uPrime, gas.SL, lI, nu);
        this.turbulentlyQuenched =
          o.turbulentQuench && turbulentQuenchReached(gas.uPrime, gas.SL, lI, nu, gas.lewisNumber ?? 1);
      } else {
        this.karlovitzAGB = 0;
        this.turbulentlyQuenched = false;
      }
    } else {
      this.lF = 0;
      this.karlovitzAGB = 0;
      this.turbulentlyQuenched = false;
    }
  }

  /** Clipped stretch factor and kernel speed at radius r, time t (sets diagnostics). */
  private flameSpeed(r: number, t: number, gas: KernelGasState): number {
    if (!(gas.SL > 0) || this.turbulentlyQuenched || gas.expansionRatio <= 1 + 1e-9) {
      this.stretchFactor = 0;
      return 0;
    }
    const lF = this.lF;
    const lI = this.lI;
    const up = gas.uPrime > 0 ? gas.uPrime : 0;
    this.karlovitzStrain = Math.sqrt(lF / (15 * this.lEps)) * Math.pow(up / gas.SL, 1.5);
    this.karlovitzCurvature = (2 * lF * this.sigmaK) / r;
    // I₀ = 1 − Ma (Ka_strain + Ka_curv), Ma = L_M/l_F (division-free: L_M/l_F · Ka = L_M κ/S_L); the
    // strain rate u′/λ with λ from the dissipation length (ε = u′³/L_ε)
    let I0 = herwegMalyStretchFactor(r, gas.SL, lF, gas.marksteinLength, up, this.lEps, this.sigmaK);
    if (!(I0 > 0)) I0 = 0;
    else if (I0 > this.opts.maxStretchFactor) I0 = this.opts.maxStretchFactor;
    this.stretchFactor = I0;
    return herwegMalyKernelSpeed(I0, r, t, gas.SL, up, lI);
  }

  /** dV/dt without electrode losses, m³/s. */
  private dVdt(V: number, t: number, power: number, gas: KernelGasState): number {
    const r = Math.cbrt((3 * V) / (4 * Math.PI));
    const A = 4 * Math.PI * r * r;
    const ST = this.flameSpeed(r, t, gas);
    return power / (this.density * this.dh) + A * this.sigmaK * ST;
  }

  private setVolume(V: number): void {
    this.volume = V;
    this.radius = Math.cbrt((3 * V) / (4 * Math.PI));
    this.burnedMass = this.density * V;
    if (this.radius > this.maxRadius) this.maxRadius = this.radius;
  }

  /**
   * Impulsive energy deposit (breakdown + capacitive arc), J to the gas. Creates the kernel
   * on the first call of a cycle; a restrike adds volume to an existing kernel.
   */
  deposit(energyToGas: number, gas: KernelGasState): void {
    if (this.stage === 'handoff' || this.stage === 'quenched' || energyToGas <= 0) return;
    this.prepare(gas);
    if (this.stage === 'none') {
      this.stage = 'kernel';
      this.age = 0;
      this.setVolume(energyToGas / (this.density * this.dh));
    } else {
      this.setVolume(this.volume + energyToGas / (this.density * this.dh));
    }
    this.energyDeposited += energyToGas;
  }

  /**
   * Advance the kernel by dt (s).
   * @param gasPower mean electrical power transferred to the gas during the step, W
   * @param sparkActive true while the discharge is still delivering energy (end of step)
   * @param gas unburned-gas state around the gap
   */
  step(dt: number, gasPower: number, sparkActive: boolean, gas: KernelGasState): void {
    if (this.stage !== 'kernel' || dt <= 0) return;
    // The hand-off radius max(r_min, l_I) follows the caller's gas state (l_I shrinks during
    // compression): take it LINEAR in time over the step, from its value at the previous step to
    // the new one, so the kernel/hand-off crossing is located inside the step whatever the caller's
    // step boundaries (a step-wise constant radius put the hand-off at a step end whenever l_I fell
    // below the kernel radius across a boundary — a crank-angle stair of CA50 / knock onset).
    const rhoPrev = this.handoffRadius;
    this.prepare(gas);
    const o = this.opts;
    const rhoDh = this.density * this.dh;
    const t0 = this.age;
    const rh1 = this.handoffRadius;
    const rh0 = rhoPrev > 0 ? rhoPrev : rh1;
    let rem = dt;
    let t = t0;
    let V = this.volume;
    let guard = 0;
    while (rem > 0 && guard++ < 100_000) {
      const k1 = this.dVdt(V, t, gasPower, gas);
      let h = rem;
      const hLim = (o.maxRelativeVolumeChange * V) / Math.max(1e-300, Math.abs(k1));
      if (hLim < h) h = Math.max(hLim, 1e-9);
      // electrode loss (∝ A/√t, singular at t = 0) limited to the same relative volume change:
      // coeff·A·(√(t+h) − √t) ≤ ε V ρΔh  ⇔  h ≤ (√t + q)² − t,  q = ε V ρΔh/(coeff·A)
      let Ac = 0;
      if (this.condCoeff > 0) {
        Ac = electrodeContactArea(this.radius, this.gap, this.centreElectrodeDiameter, this.groundElectrodeWidth);
        if (Ac > 0) {
          const sq = Math.sqrt(t) + (o.maxRelativeVolumeChange * V * rhoDh) / (this.condCoeff * Ac);
          const hLoss = sq * sq - t;
          if (hLoss < h) h = Math.max(hLoss, 1e-9);
        }
      }
      if (h > rem) h = rem;
      const Vp = Math.max(V + h * k1, 1e-30);
      const k2 = this.dVdt(Vp, t + h, gasPower, gas);
      let Vn = V + 0.5 * h * (k1 + k2);
      // electrode conduction loss over [t, t+h] (Keck boundary layer, exact in time)
      let dE = 0;
      if (this.condCoeff > 0) {
        const rm = Math.cbrt((3 * 0.5 * (V + Math.max(Vn, 0))) / (4 * Math.PI));
        Ac = electrodeContactArea(rm, this.gap, this.centreElectrodeDiameter, this.groundElectrodeWidth);
        dE = this.condCoeff * Ac * (Math.sqrt(t + h) - Math.sqrt(t));
        Vn -= dE / rhoDh;
      }
      if (Vn <= 1e-24) {
        this.energyElectrodeLoss += dE;
        t += h;
        this.energyDeposited += gasPower * (t - t0);
        this.age = t;
        this.setVolume(0);
        this.quench('kernel extinguished by electrode heat loss');
        return;
      }
      const vhA = this.handoffVolumeAt(t - t0, dt, rh0, rh1);
      const vhB = this.handoffVolumeAt(t + h - t0, dt, rh0, rh1);
      if (Vn >= vhB && V < vhA) {
        // Hand-off inside this sub-step: locate the crossing of V − V_ho(t) (linear) so the
        // hand-off time does not depend on the sub-step / caller step (reviewer fix).
        const f = (vhA - V) / (Vn - vhB + vhA - V);
        const rc = rh0 + ((rh1 - rh0) * (t + f * h - t0)) / dt;
        this.flameSpeed(rc, t + f * h, gas); // stretch factor at the hand-off radius
        if (this.stretchFactor > 0) {
          this.energyElectrodeLoss += f * dE;
          t += f * h;
          this.energyDeposited += gasPower * (t - t0);
          this.age = t;
          this.handoffRadius = rc;
          this.setVolume((4 / 3) * Math.PI * rc * rc * rc);
          this.stage = 'handoff';
          this.handoffTime = t;
          break;
        }
      }
      this.energyElectrodeLoss += dE;
      t += h;
      rem -= h;
      V = Vn;
      this.setVolume(V);
      const rhNow = rh0 + ((rh1 - rh0) * (t - t0)) / dt;
      if (this.radius >= rhNow && this.stretchFactor > 0) {
        // already beyond the hand-off radius (e.g. a large initial deposit)
        this.handoffRadius = rhNow;
        this.energyDeposited += gasPower * (t - t0);
        this.age = t;
        this.stage = 'handoff';
        this.handoffTime = t;
        break;
      }
    }
    if (this.stage === 'kernel') {
      this.energyDeposited += gasPower * (t - t0);
      this.age = t;
    }
    // growth diagnostics at the end of the step
    const A = 4 * Math.PI * this.radius * this.radius;
    this.turbulentSpeed = this.flameSpeed(this.radius, t, gas);
    this.plasmaSpeed = sparkActive ? gasPower / (A * gas.rhoU * this.dh) : 0;
    const lossRate = this.condCoeff > 0 && t > 0
      ? (this.condCoeff *
          electrodeContactArea(this.radius, this.gap, this.centreElectrodeDiameter, this.groundElectrodeWidth)) /
        (2 * Math.sqrt(t))
      : 0;
    const dV = (sparkActive ? gasPower : 0) / rhoDh + A * this.sigmaK * this.turbulentSpeed - lossRate / rhoDh;
    this.growthRate = dV / A;
    if (this.stage !== 'kernel') return;
    // ---- misfire checks (only once the discharge no longer feeds the kernel) ----
    // A kernel that cannot propagate (I₀ ≤ 0: stretch above critical, or AGB quench) is a hot-gas
    // pocket that no longer burns; it is extinguished when it has not resumed propagating within
    // nonPropagatingSurvival × min(l_I/u′, r²/α_k) — the integral eddy-turnover time (turbulent mixing
    // dilutes and cools it) or, in quiescent gas, its conductive cooling time (UNVERIFIED model
    // choice: a time scale, not a fitted constant) — or when it shrinks below half its largest radius. Its stretch changes during that time (compression raises S_L and lowers l_F, the
    // turbulence decays), so a marginal kernel may recover and hand off late. Round 1 extinguished it
    // at once: a pass/fail switch that a 14 % change of S_L (residual 3 % vs 9 %) flipped, producing
    // deterministic fire/misfire period-2 cycles at ordinary conditions (validation round 2).
    if (!sparkActive) {
      const blocked = this.turbulentlyQuenched || this.stretchFactor <= 0;
      if (blocked) {
        if (Number.isNaN(this.blockedSince)) this.blockedSince = t;
        const up = gas.uPrime > 0 ? gas.uPrime : 0;
        const tCond = this.alphaK > 0 ? (this.radius * this.radius) / this.alphaK : Infinity;
        const tLim = o.nonPropagatingSurvival * Math.min(up > 0 ? this.lI / up : Infinity, tCond);
        if (!(t - this.blockedSince < tLim)) {
          this.quench(
            this.turbulentlyQuenched
              ? 'turbulent flame quench (Abdel-Gayed–Bradley K·Le limit)'
              : 'kernel below critical radius (stretch factor I0 ≤ 0)',
          );
        }
      } else this.blockedSince = NaN;
      if (this.stage === 'kernel' && this.radius < 0.5 * this.maxRadius) this.quench('kernel shrinking (heat loss exceeds burning)');
    }
    if (this.stage === 'kernel' && this.age > o.maxKernelTime) this.quench('kernel did not reach hand-off in time');
  }

  /** Hand-off volume at time τ into a step of length dt, the radius linear from r0 to r1, m³. */
  private handoffVolumeAt(tau: number, dt: number, r0: number, r1: number): number {
    const r = r0 + ((r1 - r0) * tau) / dt;
    return (4 / 3) * Math.PI * r * r * r;
  }

  private quench(reason: string): void {
    this.stage = 'quenched';
    this.quenchReason = reason;
  }
}
