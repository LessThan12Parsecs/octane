/**
 * Engine friction: FMEP correlations and an instantaneous friction-torque model.
 *
 * 1. Patton, Nitschke & Heywood (1989), "Development and evaluation of a friction
 *    model for spark-ignition engines", SAE 890836 ("PNH") — component model.
 *    Equations and constants as reproduced in D. Sandoval, "An Improved Friction
 *    Model for Spark Ignition Engines", B.S. thesis, MIT 2002 (supervisor J. B.
 *    Heywood; basis of Sandoval & Heywood 2003, SAE 2003-01-0725):
 *      eq. 7.a crankshaft, 8.a reciprocating, 9.a ring gas loading, 10.a valvetrain,
 *      eq. 11 auxiliaries, Table 4.2 valvetrain constants, App. A.2 oil Vogel constants.
 *    Units INSIDE the correlation (as published): fmep kPa; B, S, D_b, L_b, L_v mm;
 *    N rev/min; mean piston speed S_p m/s. The public API converts to/from SI (Pa, m).
 *    Optional modifications from Sandoval (2002) §3.2–3.3 / Sandoval & Heywood (2003):
 *      - oil-viscosity scaling of EVERY hydrodynamic term by μ_scaling = √(μ/μ₀) =
 *        √(ν/ν₀) (Sandoval 2002 eq. 3 "the viscosity scaling that should be included in
 *        the hydrodynamic friction terms", eq. 4, eqs. 7.b–10.b and the total-fmep
 *        expression of App. A-2; checked on the rendered scan pages 11, 16 and 36):
 *        main and rod journal bearings, piston skirt, the 0.088·r_c gas-loading term,
 *        cam bearings and the oscillating-hydrodynamic term. ν₀ = 10.6 cSt, the 10W-30
 *        reference oil of Patton et al. at 90 °C (Sandoval 2002 §3.2; his App. A.3 sample
 *        input lists 10.3 cSt — we keep the value stated in the text);
 *      - boundary-friction speed function 1 + 1000/N → 1 + 500/N (rings, flat follower,
 *        oscillating mixed; Sandoval 2002 eqs. 8.b, 10.b).
 *    Not implemented: Sandoval's ring-tension (F_t/F_t0) / roughness (C_r) factors and the
 *    doubled K in the gas-loading exponent (engine-specific calibrations). Our
 *    `ringTensionFactor` multiplies only the ring-tension term 4.06e4/B² (Sandoval also
 *    applies F_t/F_t0 to the 0.182 gas-loading term).
 *    Pumping losses are NOT included: the cycle simulation computes pumping work
 *    from the gas-exchange pressures directly.
 *
 * 2. Chen & Flynn (1965), SAE 650733: FMEP = A + B·p_max + C·S̄p + D·S̄p², with the
 *    form as given by Pipitone (2009), SAE 2009-01-1984, eq. 1 (fetched), who writes
 *    it with engine speed n instead of S̄p (equivalent for a fixed stroke).
 *    UNVERIFIED: the Chen & Flynn paper itself (title "Development of a single
 *    cylinder compression ignition research engine") was not fetched.
 *    Coefficients are engine specific (no defaults).
 *
 * 3. Instantaneous friction torque (FrictionTorqueModel): the crank-angle-
 *    independent part (bearings, valvetrain, auxiliaries) is applied as a constant
 *    torque; ring friction as a Coulomb-type force of constant magnitude opposing
 *    piston motion (crank torque F_c|dx/dθ|); piston-skirt friction as a viscous
 *    force ∝ piston velocity. All are normalised so that the cycle-mean torque
 *    equals FMEP·V_d,total/(4π) (four-stroke; FMEP is the whole-engine value PNH returns,
 *    V_d,total = N·V_d) at the reference speed, up to the sign smoothing
 *    1 − ω/√(ω² + ε²) (3e-5 at 600 rpm with ε = 0.5 rad/s). For N cylinders the constant
 *    part carries the total displacement and the piston terms act on every cylinder at
 *    its own phase: T = −sgn(ω)[T_c + F_c Σ_i|x′(θ_i)|] − c_v Σ_i x′(θ_i)² ω with the
 *    per-cylinder F_c and c_v of the single-cylinder normalisation (⟨|x′|⟩ and ⟨x′²⟩ do
 *    not depend on the phase, so the cycle mean is exact for any firing offsets).
 *    The Coulomb term F_c|x′| has a slope discontinuity (not a jump) in θ at TDC/BDC,
 *    where the piston reverses — physical, and harmless for an explicit RK integrator.
 *    Modelling choice (documented, not from PNH): ring gas-loading friction is
 *    distributed like ring-tension friction (not weighted by the instantaneous
 *    cylinder pressure).
 *
 * Sign convention: friction torque is returned as a torque on the crank in the
 * rotation direction, i.e. NEGATIVE while ω > 0 (dynamics.ts convention).
 */

import type { SliderCrank } from './kinematics';

/** kPa → Pa. */
const KPA = 1e3;
/** m → mm. */
const MM = 1e3;

/**
 * Valvetrain layouts of PNH (Sandoval 2002, Table 4.2), plus 'L-head' (side valves: cam in the block
 * lifting each valve directly through a flat-footed "mushroom" tappet; no pushrod, no rocker).
 */
export type ValvetrainType = 'SOHC-finger' | 'SOHC-rocker' | 'SOHC-direct' | 'DOHC-finger' | 'DOHC-direct' | 'OHV' | 'L-head';

/**
 * PNH valvetrain constants [C_ff flat follower (kPa·mm), C_rf roller follower
 * (kPa·mm·min), C_oh oscillating hydrodynamic, C_om oscillating mixed (kPa)].
 * Source: Sandoval (2002) Table 4.2 (PNH 1989 values; 'DOHC-finger' was added by
 * Sandoval with C_om = 0.6 × SOHC-finger). All values re-read from the rendered scan.
 * UNVERIFIED: the OHV roller-follower entry is printed as "0.5" in Sandoval's table —
 * clearly a typo (the other C_rf are 0.005–0.0227; 0.5 would give an OHV roller follower
 * ≈ 30× the SOHC value); set to NaN here — use flat followers for OHV (the CFR has
 * flat tappets).
 * UNVERIFIED mapping 'L-head' = 'SOHC-direct' (PNH has no side-valve type). Reasoning: a side valve is
 * opened by a flat tappet sliding in a block bore and bearing directly on the valve stem — one
 * cam/flat-follower contact per valve at a motion ratio of 1, the follower side load taken by the
 * tappet bore, and the stem in a plain guide — which is the kinematics of a direct-acting (bucket)
 * follower, not of the OHV train (lifter + pushrod + rocker, motion ratio ≈ 1.5, extra pivots: C_ff 400,
 * C_om 32.1). PNH's constants embed 1980s spring loads; a light-spring engine (Ford Model T: 24–28 lb
 * installed, ≈ 32 lb open [Ford drawing T-431; Ford Service par. 258]) probably has less valvetrain
 * friction than this, a heavy-spring flathead more.
 */
export const PNH_VALVETRAIN_CONSTANTS: Readonly<Record<ValvetrainType, readonly [number, number, number, number]>> = {
  'SOHC-finger': [600, 0.0227, 0.2, 42.8],
  'SOHC-rocker': [400, 0.0151, 0.5, 21.4],
  'SOHC-direct': [200, 0.0076, 0.5, 10.7],
  'DOHC-finger': [600, 0.0227, 0.2, 25.8],
  'DOHC-direct': [133, 0.005, 0.5, 10.7],
  OHV: [400, Number.NaN, 0.5, 32.1],
  'L-head': [200, 0.0076, 0.5, 10.7],
};

/** Reference kinematic viscosity of the PNH calibration oil (10W-30 at 90 °C), cSt (Sandoval 2002 §3.2). */
export const PNH_REFERENCE_VISCOSITY_CST = 10.6;

/**
 * Vogel low-shear viscosity constants ν = k·exp(θ₁/(T + θ₂)) (cSt, T in °C) and the
 * high-shear/low-shear ratio μ_hs/μ_ls, from Sandoval (2002) App. A.2 (his ref. [8]).
 */
export const OIL_VOGEL = {
  '5W20': { k: 0.04576, theta1: 1224, theta2: 134.1, highShearRatio: 0.94 },
  '5W40': { k: 0.15, theta1: 1018.74, theta2: 125.91, highShearRatio: 0.8 },
  '10W30': { k: 0.1403, theta1: 869.72, theta2: 104.4, highShearRatio: 0.76 },
  '15W40': { k: 0.1223, theta1: 933.46, theta2: 103.89, highShearRatio: 0.9 },
  '20W50': { k: 0.0639, theta1: 1255.46, theta2: 117.7, highShearRatio: 0.84 },
  SAE10: { k: 0.0258, theta1: 1345.42, theta2: 144.58, highShearRatio: 1 },
  SAE30: { k: 0.0246, theta1: 1432.29, theta2: 132.94, highShearRatio: 1 },
  SAE50: { k: 0.0384, theta1: 1349.94, theta2: 115.16, highShearRatio: 1 },
} as const;
export type OilGrade = keyof typeof OIL_VOGEL;

/**
 * High-shear kinematic viscosity of an oil grade at temperature T (K), cSt
 * (Vogel equation × μ_hs/μ_ls; Sandoval 2002 eqs. 5–6).
 */
export function oilViscosityCst(grade: OilGrade, T: number): number {
  const o = OIL_VOGEL[grade];
  const tc = T - 273.15;
  return o.k * Math.exp(o.theta1 / (tc + o.theta2)) * o.highShearRatio;
}

/** Geometry/design inputs of the PNH model (SI units: m). */
export interface PnhFrictionInputs {
  bore: number;
  stroke: number;
  cylinders: number;
  /** Geometric compression ratio (ring gas-loading term). */
  compressionRatio: number;
  /** Main bearings: count, journal diameter (m), bearing length (m). */
  mainBearings: { count: number; diameter: number; length: number };
  /** Connecting-rod (big-end) bearings: count (= cylinders), diameter (m), length (m). */
  rodBearings: { count: number; diameter: number; length: number };
  /** Number of camshaft bearings. */
  camBearings: number;
  /** Total number of valves (all cylinders). */
  valves: number;
  /** Maximum valve lift L_v, m. */
  maxValveLift: number;
  valvetrain: ValvetrainType;
  follower: 'flat' | 'roller';
  /** Multiplier on the ring-tension term (e.g. ring count / 3); 1 = PNH. */
  ringTensionFactor?: number;
  /** Multiplier on the auxiliary (oil + water pump + alternator) term; 1 = PNH. */
  auxiliaryFactor?: number;
  /**
   * Oil viscosity ratio ν/ν₀ (operating oil vs the 10.6 cSt PNH reference); the hydrodynamic
   * terms are scaled by √(ν/ν₀) (Sandoval 2002 eq. 3, App. A-2). 1 = PNH.
   */
  viscosityRatio?: number;
  /** N₀ of the boundary-friction function 1 + N₀/N: 1000 (PNH 1989) or 500 (Sandoval & Heywood 2003). */
  boundaryRpm?: number;
}

/** PNH component breakdown, all Pa (reused output record). */
export interface PnhFmepBreakdown {
  crankSeals: number;
  mainBearings: number;
  turbulentDissipation: number;
  pistonSkirt: number;
  rings: number;
  rodBearings: number;
  ringGasLoading: number;
  camBearings: number;
  follower: number;
  oscillatingHydrodynamic: number;
  oscillatingMixed: number;
  auxiliaries: number;
  /** Crankshaft group (seals + main bearings + turbulent dissipation). */
  crankshaft: number;
  /** Reciprocating group (skirt + rings + rod bearings + gas loading). */
  reciprocating: number;
  /** Valvetrain group. */
  valvetrain: number;
  /** Crank-angle-independent part (crankshaft + rod bearings + valvetrain + auxiliaries). */
  constantPart: number;
  /** Piston-assembly part (skirt + rings + gas loading). */
  pistonPart: number;
  /** Total mechanical FMEP (no pumping). */
  total: number;
}

export const newPnhFmepBreakdown = (): PnhFmepBreakdown => ({
  crankSeals: 0,
  mainBearings: 0,
  turbulentDissipation: 0,
  pistonSkirt: 0,
  rings: 0,
  rodBearings: 0,
  ringGasLoading: 0,
  camBearings: 0,
  follower: 0,
  oscillatingHydrodynamic: 0,
  oscillatingMixed: 0,
  auxiliaries: 0,
  crankshaft: 0,
  reciprocating: 0,
  valvetrain: 0,
  constantPart: 0,
  pistonPart: 0,
  total: 0,
});

/**
 * Patton–Nitschke–Heywood (1989) mechanical FMEP, Pa, with component breakdown.
 * @param rpm engine speed, rev/min
 * @param intakePressure intake-manifold pressure p_i, Pa
 * @param ambientPressure p_a, Pa
 * @param out optional reused output record
 * @throws RangeError if rpm ≤ 0 or NaN (PNH boundary terms ∝ 1 + N₀/N; the correlation was fitted ≳ 1000 rpm)
 */
export function pnhFmep(
  inp: PnhFrictionInputs,
  rpm: number,
  intakePressure: number,
  ambientPressure: number,
  out: PnhFmepBreakdown = newPnhFmepBreakdown(),
): PnhFmepBreakdown {
  if (!(rpm > 0)) throw new RangeError('pnhFmep: rpm must be > 0 (the boundary terms ∝ 1 + N₀/N diverge at rest)');
  const B = inp.bore * MM;
  const S = inp.stroke * MM;
  const nc = inp.cylinders;
  const N = rpm;
  const Sp = (2 * inp.stroke * rpm) / 60; // m/s
  // μ_scaling = √(μ/μ₀) on every hydrodynamic term (Sandoval 2002 eq. 3 and App. A-2).
  const vr = Math.sqrt(inp.viscosityRatio ?? 1);
  const N0 = inp.boundaryRpm ?? 1000;
  const bnd = 1 + N0 / N;
  const B2S = B * B * S * nc;

  // Crankshaft (Sandoval 2002 eq. 7.a): seals (boundary), main-bearing hydrodynamic, turbulent dissipation.
  const Db = inp.mainBearings.diameter * MM;
  const Lb = inp.mainBearings.length * MM;
  const nb = inp.mainBearings.count;
  out.crankSeals = (1.22e5 * Db) / B2S;
  out.mainBearings = (3.03e-4 * vr * N * Db * Db * Db * Lb * nb) / B2S;
  out.turbulentDissipation = (1.35e-10 * Db * Db * N * N * nb) / nc;

  // Reciprocating (eq. 8.a): piston skirt (hydrodynamic), rings (mixed, 1 + N₀/N), rod bearings.
  const Dr = inp.rodBearings.diameter * MM;
  const Lr = inp.rodBearings.length * MM;
  out.pistonSkirt = (2.94e2 * vr * Sp) / B;
  out.rings = (4.06e4 * (inp.ringTensionFactor ?? 1) * bnd) / (B * B);
  out.rodBearings = (3.03e-4 * vr * N * Dr * Dr * Dr * Lr * inp.rodBearings.count) / B2S;
  // Ring gas-pressure loading (eq. 9.a), K = 2.38e-2 s/m; √(μ/μ₀) on the 0.088·r_c term (eq. 9.b, App. A-2).
  const rc = inp.compressionRatio;
  out.ringGasLoading = 6.89 * (intakePressure / ambientPressure) * (0.088 * vr * rc + 0.182 * Math.pow(rc, 1.33 - 2.38e-2 * Sp));

  // Valvetrain (eq. 10.a, Table 4.2): cam bearings (+4.12 kPa seal constant), follower, oscillating terms.
  const [Cff, Crf, Coh, Com] = PNH_VALVETRAIN_CONSTANTS[inp.valvetrain];
  if (inp.follower === 'roller' && !Number.isFinite(Crf)) {
    throw new RangeError(`pnhFmep: no roller-follower constant for valvetrain '${inp.valvetrain}'`);
  }
  const nv = inp.valves;
  const Lv = inp.maxValveLift * MM;
  out.camBearings = (244 * vr * N * inp.camBearings) / B2S + 4.12;
  out.follower = inp.follower === 'flat' ? (Cff * bnd * nv) / (S * nc) : (Crf * N * nv) / (S * nc);
  out.oscillatingHydrodynamic = (Coh * vr * Math.pow(Lv, 1.5) * Math.sqrt(N) * nv) / (B * S * nc);
  out.oscillatingMixed = (Com * bnd * Lv * nv) / (S * nc);

  // Auxiliaries (eq. 11): oil pump + water pump + non-charging alternator.
  out.auxiliaries = (inp.auxiliaryFactor ?? 1) * (6.23 + 5.22e-3 * N - 1.79e-7 * N * N);

  // kPa → Pa and groups
  out.crankSeals *= KPA;
  out.mainBearings *= KPA;
  out.turbulentDissipation *= KPA;
  out.pistonSkirt *= KPA;
  out.rings *= KPA;
  out.rodBearings *= KPA;
  out.ringGasLoading *= KPA;
  out.camBearings *= KPA;
  out.follower *= KPA;
  out.oscillatingHydrodynamic *= KPA;
  out.oscillatingMixed *= KPA;
  out.auxiliaries *= KPA;
  out.crankshaft = out.crankSeals + out.mainBearings + out.turbulentDissipation;
  out.reciprocating = out.pistonSkirt + out.rings + out.rodBearings + out.ringGasLoading;
  out.valvetrain = out.camBearings + out.follower + out.oscillatingHydrodynamic + out.oscillatingMixed;
  out.pistonPart = out.pistonSkirt + out.rings + out.ringGasLoading;
  out.constantPart = out.crankshaft + out.rodBearings + out.valvetrain + out.auxiliaries;
  out.total = out.crankshaft + out.reciprocating + out.valvetrain + out.auxiliaries;
  return out;
}

/** Chen–Flynn coefficients: FMEP = A + B·p_max + C·S̄p + D·S̄p² (SI: Pa, –, Pa·s/m, Pa·s²/m²). */
export interface ChenFlynnCoefficients {
  A: number;
  B: number;
  C: number;
  D: number;
}

/**
 * Chen & Flynn (1965, SAE 650733) FMEP, Pa.
 * @param pMax peak cylinder pressure of the cycle, Pa
 * @param meanPistonSpeed m/s
 */
export function chenFlynnFmep(c: ChenFlynnCoefficients, pMax: number, meanPistonSpeed: number): number {
  return c.A + c.B * pMax + c.C * meanPistonSpeed + c.D * meanPistonSpeed * meanPistonSpeed;
}

/** Cycle-mean friction torque of one four-stroke cylinder, FMEP·V_d/(4π), N m (magnitude). */
export function meanFrictionTorque(fmep: number, displacedVolume: number): number {
  return (fmep * displacedVolume) / (4 * Math.PI);
}

/**
 * Instantaneous friction torque on the crank, distributed over the cycle
 * (see file header, item 3). Allocation-free after construction.
 *
 *   T_f(θ, ω) = −sgn(ω)·[T_c + F_c Σ_i |x′_i|] − c_v Σ_i x′_i² ω
 *
 * T_c: constant part (bearings, valvetrain, auxiliaries) of the whole engine; F_c: Coulomb-type
 * piston friction force per cylinder (ring tension + ring gas loading, mixed/boundary lubrication);
 * c_v: viscous piston-skirt coefficient per cylinder (hydrodynamic, force ∝ piston velocity, so
 * its FMEP ∝ S̄p as in PNH). x′_i = dx/dθ of cylinder i at its own phase (single cylinder: the
 * kinematics' dx/dθ). Single-cylinder engines may use `torque(x′, ω)`; engines with N > 1 must use
 * `torqueCylinders` (one x′ would apply one piston's friction only).
 */
export class FrictionTorqueModel {
  /** Constant (crank-angle-independent) friction torque magnitude T_c of the whole engine, N m. */
  constantTorque = 0;
  /** Coulomb-type piston friction force magnitude F_c per cylinder, N. */
  coulombForce = 0;
  /** Viscous piston friction coefficient c_v per cylinder (force = c_v·dx/dt), N s/m. */
  viscousCoefficient = 0;
  /** Speed below which the Coulomb sign is smoothed, rad/s (avoids a discontinuity at ω = 0). */
  smoothingOmega: number;
  /** Number of cylinders N (identical, sharing the slider-crank). */
  readonly cylinders: number;
  /** Displaced volume per cylinder, m³. */
  readonly displacedVolume: number;
  /** Displaced volume of the engine N·V_d, m³ (the FMEP normalisation volume). */
  readonly totalDisplacedVolume: number;
  /** ⟨|x′|⟩ over a revolution = travel/π, m. */
  readonly meanAbsDxdTheta: number;
  /** ⟨x′²⟩ over a revolution, m². */
  readonly meanSqDxdTheta: number;

  /**
   * @param kin the engine's slider-crank (for the cycle normalisation of the piston terms)
   * @param smoothingOmega sign smoothing speed, rad/s
   * @param cylinders number of identical cylinders (EngineSpec.cylinders); FMEPs passed to setFmep /
   *   setFromPnh are whole-engine values normalised by N·V_d
   */
  constructor(kin: SliderCrank, smoothingOmega = 0.5, cylinders = 1) {
    if (!(Number.isInteger(cylinders) && cylinders >= 1)) throw new RangeError('FrictionTorqueModel: cylinders must be an integer ≥ 1');
    this.smoothingOmega = smoothingOmega;
    this.cylinders = cylinders;
    this.displacedVolume = kin.displacedVolume;
    this.totalDisplacedVolume = kin.displacedVolume * cylinders;
    this.meanAbsDxdTheta = kin.pistonTravel / Math.PI;
    // periodic integrand → the rectangle rule is spectrally accurate
    const n = 2048;
    let s2 = 0;
    for (let i = 0; i < n; i++) {
      const d = kin.dxdTheta((2 * Math.PI * i) / n);
      s2 += d * d;
    }
    this.meanSqDxdTheta = s2 / n;
  }

  /**
   * Set the friction level from whole-engine FMEP components (Pa) evaluated at crank speed
   * `omegaRef` (rad/s). Each component's cycle-mean torque equals fmep·N·V_d/(4π) at ω = omegaRef:
   *   T_c = fmep_c N V_d/(4π),  F_c = fmep_coul V_d/(4π⟨|x′|⟩),  c_v = fmep_visc V_d/(4π ω_ref ⟨x′²⟩)
   * (the per-cylinder piston terms act N times, once per cylinder).
   */
  setFmep(fmepConstant: number, fmepCoulomb: number, fmepViscous: number, omegaRef: number): void {
    const k = this.displacedVolume / (4 * Math.PI);
    this.constantTorque = fmepConstant * (this.totalDisplacedVolume / (4 * Math.PI));
    this.coulombForce = (fmepCoulomb * k) / this.meanAbsDxdTheta;
    this.viscousCoefficient = (fmepViscous * k) / (Math.abs(omegaRef) * this.meanSqDxdTheta);
  }

  /** Convenience: set from a PNH breakdown (skirt → viscous; rings + gas loading → Coulomb). */
  setFromPnh(b: PnhFmepBreakdown, omegaRef: number): void {
    this.setFmep(b.constantPart, b.rings + b.ringGasLoading, b.pistonSkirt, omegaRef);
  }

  /**
   * Friction torque on the crank in the rotation direction, N m (negative for ω > 0), of a
   * SINGLE-cylinder engine (cylinders = 1; use torqueCylinders otherwise).
   * @param dxdTheta kinematic dx/dθ at the current angle (SliderCrank.dxdTheta), m/rad
   * @param omega crank angular velocity, rad/s
   */
  torque(dxdTheta: number, omega: number): number {
    const e = this.smoothingOmega;
    const coul = this.constantTorque + this.coulombForce * Math.abs(dxdTheta);
    return (-coul * omega) / Math.sqrt(omega * omega + e * e) - this.viscousCoefficient * dxdTheta * dxdTheta * omega;
  }

  /**
   * Friction torque of all cylinders on the crank in the rotation direction, N m (negative for ω > 0):
   *   T = −sgn_ε(ω)[T_c + F_c Σ_i |x′_i|] − Σ_i c_v x′_i² ω
   * Bit-identical to `torque(x′, ω)` for one cylinder.
   * @param dxdThetas per-cylinder dx/dθ at each cylinder's own phase, m/rad (first `cylinders` entries
   *   used; MultiCylinderCrankTrain.update(θ) returns exactly this array)
   * @param omega crank angular velocity, rad/s
   */
  torqueCylinders(dxdThetas: ArrayLike<number>, omega: number): number {
    const e = this.smoothingOmega;
    const cv = this.viscousCoefficient;
    let sa = 0;
    let sv = 0;
    for (let i = 0; i < this.cylinders; i++) {
      const d = dxdThetas[i];
      sa += Math.abs(d);
      sv += cv * d * d;
    }
    const coul = this.constantTorque + this.coulombForce * sa;
    return (-coul * omega) / Math.sqrt(omega * omega + e * e) - sv * omega;
  }
}
