/**
 * External load on the crank in 'free' speed mode (OperatingPoint.load, core/operating-point.ts):
 * a load torque T_L(ω) and, for a vehicle, the inertia it reflects to the crank.
 *
 * ── Models ──────────────────────────────────────────────────────────────────
 *  - 'constant': T_L = loadTorque (any ω; the CFR free-mode behaviour, no inertia).
 *  - 'brake' (fan, water brake, Prony brake): T_L = loadTorque·(|n|/refRpm)^e·sgn(ω); e = 2 is the
 *    turbulent-drag law of a fan or hydraulic brake (torque ∝ ρ ω² D⁵). No inertia (the brake rotor is
 *    taken as part of J_rot).
 *  - 'vehicle': the engine drives EngineSpec.vehicle through gear ratio G (transmission) and G_f (final
 *    drive) with the clutch engaged and a rigid driveline:
 *      wheel speed ω_w = ω/(G G_f),  road speed v = s r ω/(G G_f)   (s = −1 in reverse, else +1)
 *      road load   F = m g (C_rr cos a·sgn_ε(v) + sin a) + ½ ρ C_d A v|v|,   a = atan(grade)
 *    (rolling resistance and drag oppose the motion, gravity acts down the slope; F is the force the
 *    tyres must push forward with, positive = resisting forward motion). The crank torque follows from
 *    the power balance with a constant driveline efficiency η:
 *      T_w = F s r/(G G_f)                        (lossless crank torque, = F v/ω)
 *      T_L = T_w/η   while the engine drives the car (T_w ω > 0)
 *      T_L = T_w·η   on overrun (the car drives the engine, T_w ω < 0)
 *    so the driveline always dissipates (P_loss = (1/η − 1)|P_w| driving, (1 − η)|P_w| overrun). The
 *    car's translating mass and the wheels' rotary inertia are reflected to the crank as the constant
 *      J_L = (m r² + J_wheels)/(G G_f)²
 *    (kinetic energy ½ m v² + ½ J_w ω_w² = ½ J_L ω²; η is not applied to inertia, so the coasting
 *    system conserves energy). Gear 'neutral' = declutched: T_L = 0, J_L = 0 on the crank, and the car
 *    COASTS on its own road load (the same rolling, grade and drag terms, no driveline):
 *      m_c dv/dt = −F(v),   m_c = m + J_wheels/r²   (coastSpeedAfter; the car's road speed is then a
 *    state of its own, tracked by the cycle model, which also applies the inelastic clutch engagement when a
 *    gear is selected — CycleModel.setOperatingPoint).
 *    Road-load formula: e.g. Gillespie, Fundamentals of Vehicle Dynamics (SAE 1992) ch. 4 —
 *    UNVERIFIED chapter reference; the terms are the standard rolling/aerodynamic/grade resistances.
 *
 * Sign convention: T_L > 0 OPPOSES rotation (OperatingPoint.loadTorque convention); the crank sees
 * −T_L. The free-speed RHS is
 *   α = crank.angularAcceleration(θ, ω, p, p_cc, −load.torque(ω) + T_friction, load.inertia()).
 *
 * Hot paths (torque, inertia, vehicleSpeed, coastSpeedAfter) allocate nothing; configure() validates and
 * precomputes.
 */
import type { VehicleSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import { R_UNIVERSAL } from '../core/constants';
import { G_STANDARD } from './dynamics';

/** Molar mass of dry air, kg/mol (U.S. Standard Atmosphere 1976, M₀ = 28.9644 kg/kmol). */
export const DRY_AIR_MOLAR_MASS = 28.9644e-3;

/** Resolved load model kinds ('vehicle' with gear 'neutral' resolves to 'neutral'). */
export type LoadModelKind = 'constant' | 'brake' | 'vehicle' | 'neutral';

/** The operating-point fields the load model reads. */
export type LoadOperatingPoint = Pick<OperatingPoint, 'loadTorque' | 'load' | 'ambientPressure' | 'ambientTemperature'>;

export interface LoadModelOptions {
  /** Air density for the aerodynamic drag, kg/m³. Default: dry air at the operating point's ambient p, T. */
  airDensity?: number;
  /**
   * Rotary inertia of the road wheels (+ tyres, brake drums, axle shafts) about their axes, all wheels
   * together, kg m² (reflected through the full ratio). Default 0 (not in VehicleSpec).
   */
  wheelInertia?: number;
  /**
   * Road speed below which the rolling-resistance sign is smoothed, m/s: sgn_ε(v) = v/√(v² + ε²)
   * (keeps the right-hand side continuous at v = 0). Default 0.05.
   */
  rollingSmoothingSpeed?: number;
  /** Gravitational acceleration, m/s². Default G_STANDARD. */
  gravity?: number;
}

/** Gear name meaning "declutched" (OperatingPoint LoadSpec 'vehicle'). */
export const NEUTRAL_GEAR = 'neutral';
/** Gear name driving the car backwards (its VehicleSpec ratio is given as a positive number). */
export const REVERSE_GEAR = 'reverse';

/**
 * Load torque and reflected inertia for one operating point. Construct once per engine, call
 * `configure(op)` whenever the operating point changes, then use torque/inertia/vehicleSpeed in the
 * integrator.
 */
export class LoadModel {
  /** The vehicle of the engine spec (required by the 'vehicle' load). */
  readonly vehicle: VehicleSpec | undefined;
  /** Rolling-resistance sign smoothing speed ε, m/s. */
  readonly rollingSmoothingSpeed: number;
  /** Gravitational acceleration, m/s². */
  readonly gravity: number;
  /** Wheel rotary inertia (all wheels), kg m². */
  readonly wheelInertia: number;
  private readonly airDensityOverride: number | undefined;

  /** Active model. */
  kind: LoadModelKind = 'constant';
  /** Engaged gear ('vehicle'/'neutral'), else ''. */
  gear = '';
  /** Road grade, rise/run ('vehicle'). */
  grade = 0;
  /** Air density used for drag, kg/m³ ('vehicle'). */
  airDensity = 0;
  /** Signed road travel per crank radian s·r/(G G_f), m/rad (0 unless a gear is engaged). */
  roadSpeedPerOmega = 0;
  /**
   * Translating mass of the vehicle with the wheels' rotary inertia, m_c = m + J_wheels/r², kg ('vehicle' and
   * 'neutral'; else 0). In gear J_L = m_c·roadSpeedPerOmega².
   */
  coastMass = 0;

  private t0 = 0;
  private omegaRef = 1;
  private exponent = 0;
  private jLoad = 0;
  private fGrade = 0;
  private fRoll = 0;
  private halfRhoCdA = 0;
  private eta = 1;
  private invEta = 1;

  constructor(vehicle?: VehicleSpec, options: LoadModelOptions = {}) {
    this.vehicle = vehicle;
    this.rollingSmoothingSpeed = options.rollingSmoothingSpeed ?? 0.05;
    this.gravity = options.gravity ?? G_STANDARD;
    this.wheelInertia = options.wheelInertia ?? 0;
    this.airDensityOverride = options.airDensity;
    if (!(this.rollingSmoothingSpeed > 0)) throw new RangeError('LoadModel: rollingSmoothingSpeed must be > 0');
    if (!(this.wheelInertia >= 0)) throw new RangeError('LoadModel: wheelInertia must be ≥ 0');
    if (this.airDensityOverride !== undefined && !(this.airDensityOverride >= 0)) throw new RangeError('LoadModel: airDensity must be ≥ 0');
  }

  /**
   * Select and precompute the load of an operating point (not a hot path). `op.load` absent =
   * { kind: 'constant' }.
   * @throws RangeError for an invalid brake law, a 'vehicle' load without EngineSpec.vehicle, an unknown
   *   gear or invalid vehicle data
   */
  configure(op: LoadOperatingPoint): this {
    const ls = op.load ?? { kind: 'constant' as const };
    this.t0 = op.loadTorque;
    this.gear = '';
    this.grade = 0;
    this.roadSpeedPerOmega = 0;
    this.jLoad = 0;
    this.airDensity = 0;
    this.coastMass = 0;
    if (ls.kind === 'constant') {
      this.kind = 'constant';
    } else if (ls.kind === 'brake') {
      if (!(ls.refRpm > 0) || !(ls.exponent >= 0) || !Number.isFinite(ls.exponent)) {
        throw new RangeError('LoadModel: brake load needs refRpm > 0 and a finite exponent ≥ 0');
      }
      this.kind = 'brake';
      this.omegaRef = (ls.refRpm * 2 * Math.PI) / 60;
      this.exponent = ls.exponent;
    } else {
      const veh = this.vehicle;
      if (!veh) throw new RangeError("LoadModel: 'vehicle' load but the engine spec has no vehicle");
      this.gear = ls.gear;
      this.grade = ls.grade;
      if (!Number.isFinite(ls.grade)) throw new RangeError('LoadModel: grade must be finite');
      if (ls.gear === NEUTRAL_GEAR) {
        // declutched: no torque or inertia on the crank; the coasting car's road load (coastSpeedAfter)
        if (!(veh.wheelRadius > 0 && veh.mass > 0)) throw new RangeError('LoadModel: wheel radius and mass must be > 0');
        if (!(veh.rollingResistance >= 0 && veh.dragArea >= 0)) throw new RangeError('LoadModel: rollingResistance and dragArea must be ≥ 0');
        this.setRoadLoad(veh, op, ls.grade);
        this.kind = 'neutral';
        return this;
      }
      // (own gears only: an inherited name — 'toString', 'constructor', '__proto__', … — is an unknown gear)
      const G = Object.hasOwn(veh.gears, ls.gear) ? veh.gears[ls.gear] : undefined;
      if (G === undefined) throw new RangeError(`LoadModel: unknown gear '${ls.gear}' (have ${Object.keys(veh.gears).join(', ')}, ${NEUTRAL_GEAR})`);
      if (!(G > 0 && veh.finalDrive > 0 && veh.wheelRadius > 0 && veh.mass > 0)) {
        throw new RangeError('LoadModel: gear ratio, final drive, wheel radius and mass must be > 0');
      }
      if (!(veh.drivelineEfficiency > 0 && veh.drivelineEfficiency <= 1)) throw new RangeError('LoadModel: drivelineEfficiency must be in (0, 1]');
      if (!(veh.rollingResistance >= 0 && veh.dragArea >= 0)) throw new RangeError('LoadModel: rollingResistance and dragArea must be ≥ 0');
      const ratio = G * veh.finalDrive;
      const s = ls.gear === REVERSE_GEAR ? -1 : 1;
      this.setRoadLoad(veh, op, ls.grade);
      this.kind = 'vehicle';
      this.roadSpeedPerOmega = (s * veh.wheelRadius) / ratio;
      this.jLoad = (veh.mass * veh.wheelRadius * veh.wheelRadius + this.wheelInertia) / (ratio * ratio);
      this.eta = veh.drivelineEfficiency;
      this.invEta = 1 / veh.drivelineEfficiency;
    }
    return this;
  }

  /** Road-load coefficients of the vehicle (in gear or declutched): air density, grade, rolling and drag terms, m_c. */
  private setRoadLoad(veh: VehicleSpec, op: LoadOperatingPoint, grade: number): void {
    const rho = this.airDensityOverride ?? (op.ambientPressure * DRY_AIR_MOLAR_MASS) / (R_UNIVERSAL * op.ambientTemperature);
    if (!(rho >= 0) || !Number.isFinite(rho)) throw new RangeError('LoadModel: invalid air density (ambient p, T)');
    const ang = Math.atan(grade);
    const mg = veh.mass * this.gravity;
    this.airDensity = rho;
    this.fGrade = mg * Math.sin(ang);
    this.fRoll = mg * veh.rollingResistance * Math.cos(ang);
    this.halfRhoCdA = 0.5 * rho * veh.dragArea;
    this.coastMass = veh.mass + this.wheelInertia / (veh.wheelRadius * veh.wheelRadius);
  }

  /** Inertia reflected to the crank, kg m² (vehicle in gear: (m r² + J_wheels)/(G G_f)²; else 0). */
  inertia(): number {
    return this.jLoad;
  }

  /** Road speed at crank speed ω (rad/s), m/s (negative in reverse; 0 unless a gear is engaged). */
  vehicleSpeed(omega: number): number {
    return this.roadSpeedPerOmega * omega;
  }

  /**
   * Road load at road speed v (m/s): the forward force the tyres must provide, N
   * (m g (C_rr cos a·sgn_ε(v) + sin a) + ½ ρ C_d A v|v|; 0 unless a gear is engaged).
   */
  roadLoadForce(v: number): number {
    if (this.kind !== 'vehicle') return 0;
    return this.resistance(v);
  }

  /**
   * Road speed of the DECLUTCHED vehicle (gear 'neutral') after coasting for h seconds from road speed v (m/s),
   * m/s: m_c dv/dt = −F(v) with the road load F of roadLoadForce (rolling, grade, drag; no driveline), classical
   * RK4 in ⌈h/h_max⌉ equal sub-steps, h_max = 0.5/λ with λ ≥ |∂(F/m_c)/∂v| (the smoothed rolling term near v = 0
   * and the drag term) — far inside RK4's stability limit (2.78/λ). Other loads (no coasting car): v unchanged.
   */
  coastSpeedAfter(v: number, h: number): number {
    if (this.kind !== 'neutral' || !(h > 0)) return v;
    const im = 1 / this.coastMass;
    const lam = (this.fRoll / this.rollingSmoothingSpeed + 2 * this.halfRhoCdA * Math.abs(v)) * im;
    const n = lam * h > 0.5 ? Math.ceil((lam * h) / 0.5) : 1;
    const hs = h / n;
    for (let i = 0; i < n; i++) {
      const k1 = -im * this.resistance(v);
      const k2 = -im * this.resistance(v + 0.5 * hs * k1);
      const k3 = -im * this.resistance(v + 0.5 * hs * k2);
      const k4 = -im * this.resistance(v + hs * k3);
      v += (hs / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
    }
    return v;
  }

  /** Road resistance of the vehicle at road speed v, N (m g (C_rr cos a·sgn_ε(v) + sin a) + ½ ρ C_d A v|v|). */
  private resistance(v: number): number {
    const e = this.rollingSmoothingSpeed;
    return this.fGrade + (this.fRoll * v) / Math.sqrt(v * v + e * e) + this.halfRhoCdA * v * Math.abs(v);
  }

  /** Load torque at crank speed ω (rad/s), N m (positive opposes rotation; the crank sees −torque). */
  torque(omega: number): number {
    switch (this.kind) {
      case 'constant':
        return this.t0;
      case 'brake': {
        const r = Math.abs(omega) / this.omegaRef;
        const f = this.exponent === 2 ? r * r : Math.pow(r, this.exponent);
        return omega >= 0 ? this.t0 * f : -this.t0 * f;
      }
      case 'vehicle': {
        const k = this.roadSpeedPerOmega;
        const tw = this.roadLoadForce(k * omega) * k;
        return tw * omega >= 0 ? tw * this.invEta : tw * this.eta;
      }
      default:
        return 0;
    }
  }

  /** Power absorbed by the load, T_L·ω, W. */
  power(omega: number): number {
    return this.torque(omega) * omega;
  }
}

/**
 * Steady free-running speed: the crank speed in [omegaLo, omegaHi] (rad/s) where the engine's mean
 * torque equals the load torque, T_e(ω) = T_L(ω), found by bisection on T_e − T_L (stable when
 * T_e − T_L falls through zero). Returns NaN when T_e − T_L does not change sign on the interval.
 * Not a hot path (calls `engineTorque` ~60 times).
 * @param engineTorque mean brake torque available at the crank (after friction), N m, as a function of ω
 */
export function equilibriumOmega(engineTorque: (omega: number) => number, load: LoadModel, omegaLo: number, omegaHi: number, tol = 1e-9): number {
  let lo = omegaLo;
  let hi = omegaHi;
  let flo = engineTorque(lo) - load.torque(lo);
  const fhi = engineTorque(hi) - load.torque(hi);
  if (!(flo * fhi <= 0)) return Number.NaN;
  for (let i = 0; i < 200 && hi - lo > tol * Math.max(1, Math.abs(hi)); i++) {
    const mid = 0.5 * (lo + hi);
    const fm = engineTorque(mid) - load.torque(mid);
    if (fm * flo > 0) {
      lo = mid;
      flo = fm;
    } else {
      hi = mid;
    }
  }
  return 0.5 * (lo + hi);
}
