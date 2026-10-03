/**
 * Crank-train dynamics: gas torque, exact rigid-body inertia torque of the
 * piston + connecting rod, gravity, bearing/liner forces, and the crank angular
 * acceleration for 'free' speed mode.
 *
 * ── Model ───────────────────────────────────────────────────────────────────
 * One degree of freedom, the crank angle θ (kinematics.ts conventions). Bodies:
 *  - crank + flywheel + everything rigidly rotating with it: inertia J_rot about
 *    the crank axis (EngineSpec.masses.rotatingInertia; must NOT include the rod),
 *    assumed balanced about its axis (gravity-neutral);
 *  - piston assembly (piston + rings + wrist pin), mass m_p, translating on the
 *    wrist-pin line;
 *  - connecting rod as a full planar rigid body: mass m_r, CG at distance l_g
 *    from the big end on the big-end→small-end line, inertia I_g about its CG.
 *
 * Lagrange's equation with kinetic energy T = ½ (J_rot + J_m(θ)) ω², where the
 * reflected mechanism inertia is
 *     J_m(θ) = m_p y_p′² + m_r |r_g′|² + I_g β′²        (′ = d/dθ)
 * (y_p = wrist-pin height, r_g = rod-CG position, β = rod angle), gives
 *     (J_rot + J_m) α = T_gas + T_grav − ½ J_m′ ω² + T_ext.
 * Hence the mechanism inertia torque acting on the crank is
 *     T_inertia = −(J_m α + ½ J_m′ ω²)
 *               = −[m_p a_p·y_p′ + m_r a_g·r_g′ + I_g β̈ β′]     (virtual power form)
 * which is exactly 2π-periodic and integrates to zero over a revolution at
 * constant ω. (See e.g. Taylor 1985, "The Internal-Combustion Engine in Theory
 * and Practice" vol. 2 ch. 8, for the classical two-mass treatment.)
 *
 * The classical TWO-MASS approximation replaces the rod by m_r(1−λ) at the crank
 * pin and m_r λ at the wrist pin (λ = l_g/l); it is exact only if
 * I_g = m_r l_g (l − l_g). Its error is the rod "inertia couple"
 *     T_exact − T_2mass = −ΔI (β′² α + β′ β″ ω²),  ΔI = I_g − m_r l_g (l − l_g),
 * i.e. O(ΔI (a/l)² ω²). Both are provided; the test checks this identity.
 *
 * T_gas = (p_cyl − p_cc) A_p dx/dθ = (p_cyl − p_cc) dV/dθ (virtual work, exact
 * for any pin offset). Gravity acts along world −y (cylinder axis vertical).
 *
 * ── Sign conventions ────────────────────────────────────────────────────────
 *  - Torques are positive when they DRIVE the crank in its direction of
 *    rotation (θ increasing). ω > 0 in normal running.
 *  - gasForce > 0 pushes the piston away from the head.
 *  - Forces are world-frame (x, y) components (kinematics.ts frame).
 *    wristPinForce = force of the rod ON the piston; crankPinForce = force of the
 *    rod ON the crank pin; sideThrust = x-force of the liner ON the piston
 *    (> 0 means the piston bears on the −x, i.e. major-thrust, wall).
 *  - rodAxialForce > 0 = compression (at the small end).
 *  - Ring/skirt friction is not included in the force balance (friction.ts
 *    supplies it as a crank torque).
 *
 * The Newton–Euler force solution is an independent route to the crank torque;
 * `crankTorque` must equal gasTorque + inertiaTorque + gravityTorque (tested).
 *
 * Hot paths allocate nothing: `evaluate` writes into a caller-owned
 * CrankTrainState.
 */
import type { MassSpec } from '../core/engine-spec';
import { newKinematicState, type KinematicState, type SliderCrank } from './kinematics';

/** Standard gravity, m/s² (CGPM 1901 / ISO 80000-3 conventional value). */
export const G_STANDARD = 9.80665;

/** The subset of MassSpec used here (kg, m, kg m²). */
export type CrankTrainMasses = Pick<MassSpec, 'piston' | 'conRod' | 'conRodCgFromBigEnd' | 'conRodInertiaCg' | 'rotatingInertia'>;

/** Output record of CrankTrainDynamics.evaluate (SI; reused, never reallocated). */
export interface CrankTrainState {
  /** (p_cyl − p_crankcase)·A_p, N (> 0 pushes the piston away from the head). */
  gasForce: number;
  /** Gas-pressure torque on the crank, N m. */
  gasTorque: number;
  /** Exact rigid-body inertia torque of piston + rod on the crank, N m. */
  inertiaTorque: number;
  /** Two-mass-approximation inertia torque (same rotating/reciprocating split), N m. */
  inertiaTorqueTwoMass: number;
  /** Gravity torque of piston + rod on the crank, N m. */
  gravityTorque: number;
  /** Reflected mechanism inertia J_m(θ) (piston + rod, excluding J_rot), kg m². */
  mechanismInertia: number;
  /** dJ_m/dθ, kg m²/rad. */
  mechanismInertiaDerivative: number;
  /** Force of the rod on the piston (x, y), N. */
  wristPinForceX: number;
  wristPinForceY: number;
  /** Force of the rod on the crank pin (x, y), N. */
  crankPinForceX: number;
  crankPinForceY: number;
  /** Liner normal force on the piston (x), N. */
  sideThrust: number;
  /** Rod axial force at the small end, N (> 0 compression). */
  rodAxialForce: number;
  /** Crank torque from the Newton–Euler rod force (= gas + inertia + gravity), N m. */
  crankTorque: number;
}

/** Allocate a zeroed CrankTrainState (do this once, outside the hot loop). */
export const newCrankTrainState = (): CrankTrainState => ({
  gasForce: 0,
  gasTorque: 0,
  inertiaTorque: 0,
  inertiaTorqueTwoMass: 0,
  gravityTorque: 0,
  mechanismInertia: 0,
  mechanismInertiaDerivative: 0,
  wristPinForceX: 0,
  wristPinForceY: 0,
  crankPinForceX: 0,
  crankPinForceY: 0,
  sideThrust: 0,
  rodAxialForce: 0,
  crankTorque: 0,
});

export interface CrankTrainOptions {
  /** Gravitational acceleration along world −y, m/s² (0 disables gravity). Default G_STANDARD. */
  gravity?: number;
}

export class CrankTrainDynamics {
  readonly kinematics: SliderCrank;
  /** Piston assembly mass, kg. */
  readonly pistonMass: number;
  /** Connecting-rod mass, kg. */
  readonly rodMass: number;
  /** λ = l_g/l, rod CG position as a fraction of rod length from the big end. */
  readonly rodCgFraction: number;
  /** Rod moment of inertia about its CG, kg m². */
  readonly rodInertiaCg: number;
  /** Crank + flywheel (+ reflected driveline) inertia, kg m². */
  readonly rotatingInertia: number;
  /** Gravitational acceleration along −y, m/s². */
  readonly gravity: number;
  /** Two-mass split: rod mass lumped at the wrist pin m_r λ, kg. */
  readonly rodReciprocatingMass: number;
  /** Two-mass split: rod mass lumped at the crank pin m_r (1 − λ), kg. */
  readonly rodRotatingMass: number;
  /** ΔI = I_g − m_r l_g (l − l_g): rod inertia not represented by the two-mass model, kg m². */
  readonly rodInertiaDefect: number;
  /** Kinematics scratch (hot path, no allocation). */
  private readonly ks: KinematicState = newKinematicState();

  constructor(kinematics: SliderCrank, masses: CrankTrainMasses, options: CrankTrainOptions = {}) {
    const l = kinematics.rodLength;
    const lg = masses.conRodCgFromBigEnd;
    if (!(masses.piston >= 0 && masses.conRod >= 0 && masses.conRodInertiaCg >= 0 && masses.rotatingInertia > 0)) {
      throw new RangeError('CrankTrainDynamics: masses/inertias must be ≥ 0 (rotatingInertia > 0)');
    }
    if (!(lg >= 0 && lg <= l)) throw new RangeError('CrankTrainDynamics: rod CG must lie between the pins');
    this.kinematics = kinematics;
    this.pistonMass = masses.piston;
    this.rodMass = masses.conRod;
    this.rodCgFraction = lg / l;
    this.rodInertiaCg = masses.conRodInertiaCg;
    this.rotatingInertia = masses.rotatingInertia;
    this.gravity = options.gravity ?? G_STANDARD;
    this.rodReciprocatingMass = masses.conRod * this.rodCgFraction;
    this.rodRotatingMass = masses.conRod * (1 - this.rodCgFraction);
    this.rodInertiaDefect = masses.conRodInertiaCg - masses.conRod * lg * (l - lg);
  }

  /** Gas-pressure torque (p_cyl − p_crankcase)·dV/dθ, N m (pressures in Pa). */
  gasTorque(theta: number, pCyl: number, pCrankcase: number): number {
    return (pCyl - pCrankcase) * this.kinematics.dVdTheta(theta);
  }

  /** Reflected mechanism inertia J_m(θ) of piston + rod, kg m². */
  mechanismInertia(theta: number): number {
    const k = this.kinematics;
    const a = k.crankRadius;
    const phi = theta + k.tdcAxisAngle;
    const yp1 = -k.dxdTheta(theta);
    const b1 = k.dRodAngledTheta(theta);
    const lam = this.rodCgFraction;
    const gx1 = (1 - lam) * a * Math.cos(phi);
    const gy1 = (1 - lam) * -a * Math.sin(phi) + lam * yp1;
    return this.pistonMass * yp1 * yp1 + this.rodMass * (gx1 * gx1 + gy1 * gy1) + this.rodInertiaCg * b1 * b1;
  }

  /**
   * Full evaluation at crank angle θ (rad), speed ω (rad/s), acceleration α (rad/s²),
   * cylinder and crankcase pressures (Pa). Writes into and returns `out`.
   */
  evaluate(theta: number, omega: number, alpha: number, pCyl: number, pCrankcase: number, out: CrankTrainState): CrankTrainState {
    const k = this.kinematics;
    const ks = k.evaluate(theta, this.ks);
    const a = k.crankRadius;
    const l = k.rodLength;
    const sinp = ks.sinPhi;
    const cosp = ks.cosPhi;
    const x1 = ks.dxdTheta;
    const x2 = ks.d2xdTheta2;
    const b1 = ks.dBetadTheta;
    const b2 = ks.d2BetadTheta2;
    const sb = ks.sinBeta;
    const cb = ks.cosBeta;

    const mp = this.pistonMass;
    const mr = this.rodMass;
    const Ig = this.rodInertiaCg;
    const lam = this.rodCgFraction;
    const g = this.gravity;
    const w2 = omega * omega;

    // wrist pin: y_p' = −x', y_p'' = −x''
    const yp1 = -x1;
    const yp2 = -x2;
    // crank pin r_c = a(sinφ, cosφ)
    const cx = a * sinp;
    const cy = a * cosp;
    const cx1 = cy;
    const cy1 = -cx;
    const cx2 = -cx;
    const cy2 = -cy;
    // rod CG r_g = (1−λ) r_c + λ r_w
    const gx1 = (1 - lam) * cx1;
    const gy1 = (1 - lam) * cy1 + lam * yp1;
    const gx2 = (1 - lam) * cx2;
    const gy2 = (1 - lam) * cy2 + lam * yp2;

    // accelerations
    const ap = yp1 * alpha + yp2 * w2; // piston (along +y)
    const agx = gx1 * alpha + gx2 * w2;
    const agy = gy1 * alpha + gy2 * w2;
    const bdd = b1 * alpha + b2 * w2;

    // Lagrangian quantities
    const Jm = mp * yp1 * yp1 + mr * (gx1 * gx1 + gy1 * gy1) + Ig * b1 * b1;
    const dJm = 2 * (mp * yp1 * yp2 + mr * (gx1 * gx2 + gy1 * gy2) + Ig * b1 * b2);
    const Fg = (pCyl - pCrankcase) * k.boreArea;
    const Tgas = Fg * x1;
    const Tin = -(mp * ap * yp1 + mr * (agx * gx1 + agy * gy1) + Ig * bdd * b1);
    const Tgrav = -g * (mp * yp1 + mr * gy1);
    // two-mass: (m_p + m_r λ) at the wrist pin, m_r(1−λ) at the crank pin (|r_c'|² = a², const)
    const mrec = mp + this.rodReciprocatingMass;
    const Tin2 = -(mrec * ap * yp1 + this.rodRotatingMass * a * a * alpha);

    // Newton–Euler: piston (y), rod (translation + rotation about CG)
    const Fwy = mp * ap + Fg + mp * g; // force of rod on piston, y
    const dx = -l * sb; // d = r_w − r_c = l(−sinβ, cosβ)
    const dy = l * cb;
    const Fwx = (Ig * bdd + lam * mr * (dx * (agy + g) - dy * agx) + dx * Fwy) / dy;
    const Fcx = mr * agx + Fwx; // force of crank pin on rod
    const Fcy = mr * agy + Fwy + mr * g;

    out.gasForce = Fg;
    out.gasTorque = Tgas;
    out.inertiaTorque = Tin;
    out.inertiaTorqueTwoMass = Tin2;
    out.gravityTorque = Tgrav;
    out.mechanismInertia = Jm;
    out.mechanismInertiaDerivative = dJm;
    out.wristPinForceX = Fwx;
    out.wristPinForceY = Fwy;
    out.crankPinForceX = -Fcx; // rod on crank pin = −(crank pin on rod)
    out.crankPinForceY = -Fcy;
    out.sideThrust = -Fwx;
    out.rodAxialForce = Fwx * -sb + Fwy * cb;
    // drive torque (clockwise about +z is the rotation direction): T = r_c × F_cr,rod-on-crank reversed
    out.crankTorque = cx * Fcy - cy * Fcx;
    return out;
  }

  /**
   * Crank angular acceleration for 'free' speed mode, rad/s²:
   *   α = (T_gas + T_grav − ½ J_m′ ω² + T_ext) / (J_rot + J_m)
   * @param externalTorque sum of all other torques on the crank in the rotation
   *   direction, N m: e.g. −loadTorque (OperatingPoint.loadTorque opposes
   *   rotation) + friction torque (negative while ω > 0) + motor torque.
   */
  angularAcceleration(theta: number, omega: number, pCyl: number, pCrankcase: number, externalTorque: number): number {
    const k = this.kinematics;
    const ks = k.evaluate(theta, this.ks);
    const a = k.crankRadius;
    const sinp = ks.sinPhi;
    const cosp = ks.cosPhi;
    const x1 = ks.dxdTheta;
    const x2 = ks.d2xdTheta2;
    const b1 = ks.dBetadTheta;
    const b2 = ks.d2BetadTheta2;
    const lam = this.rodCgFraction;
    const yp1 = -x1;
    const yp2 = -x2;
    const gx1 = (1 - lam) * a * cosp;
    const gy1 = -(1 - lam) * a * sinp + lam * yp1;
    const gx2 = -(1 - lam) * a * sinp;
    const gy2 = -(1 - lam) * a * cosp + lam * yp2;
    const mp = this.pistonMass;
    const mr = this.rodMass;
    const Ig = this.rodInertiaCg;
    const Jm = mp * yp1 * yp1 + mr * (gx1 * gx1 + gy1 * gy1) + Ig * b1 * b1;
    const halfdJm = mp * yp1 * yp2 + mr * (gx1 * gx2 + gy1 * gy2) + Ig * b1 * b2;
    const Tgas = (pCyl - pCrankcase) * k.boreArea * x1;
    const Tgrav = -this.gravity * (mp * yp1 + mr * gy1);
    return (Tgas + Tgrav - halfdJm * omega * omega + externalTorque) / (this.rotatingInertia + Jm);
  }

  /**
   * Torque the speed-holding machine (CFR: belt-coupled synchronous motor/
   * dynamometer) must apply to the crank to impose acceleration α (0 for a
   * perfectly held speed), N m (> 0 = motoring, < 0 = absorbing):
   *   T_hold = (J_rot + J_m) α − (T_gas + T_grav − ½ J_m′ ω² + T_ext)
   */
  holdingTorque(theta: number, omega: number, alpha: number, pCyl: number, pCrankcase: number, externalTorque: number): number {
    const Jm = this.mechanismInertia(theta);
    const aFree = this.angularAcceleration(theta, omega, pCyl, pCrankcase, externalTorque);
    return (this.rotatingInertia + Jm) * (alpha - aFree);
  }

  /** Total kinetic energy ½ (J_rot + J_m(θ)) ω², J. */
  kineticEnergy(theta: number, omega: number): number {
    return 0.5 * (this.rotatingInertia + this.mechanismInertia(theta)) * omega * omega;
  }

  /** Gravitational potential energy of piston + rod relative to the crank axis, J. */
  potentialEnergy(theta: number): number {
    const k = this.kinematics;
    const phi = theta + k.tdcAxisAngle;
    const yp = k.pinHeight(theta);
    const yg = (1 - this.rodCgFraction) * k.crankRadius * Math.cos(phi) + this.rodCgFraction * yp;
    return this.gravity * (this.pistonMass * yp + this.rodMass * yg);
  }
}
