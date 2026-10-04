/**
 * Trembler-coil vibrator: the spring-steel armature over the open end of the coil core that carries the
 * moving point (Ford/K-W coil, Williams patent US 1,092,417). SI units (m, m/s, A, s).
 *
 * 1-DOF model (core/engine-spec.ts TremblerVibratorSpec): armature travel x ≥ 0 toward the core,
 *     m ẍ = K I₁²/(g₀ − x)² − F₀ − k_s x − c ẋ,      0 ≤ x ≤ x_max,
 * with k_s = m ω_n², c = 2 ζ m ω_n, ω_n = 2π f_n, the spring preload F₀ = k_s x_p, and the reluctance
 * pull K I²/(g₀ − x)² (an inverse-square gap law; g₀ = `airGap` is an EFFECTIVE gap: it also stands for
 * the return path of the open core) with K set so that the pull equals the preload at x = 0 and
 * I₁ = I_p = `pullCurrent`. Dividing by m:
 *     ẍ = ω_n² [ x_p ((I/I_p)² (g₀/(g₀ − x))² − 1) − x ] − 2 ζ ω_n ẋ.
 * I_p is the STATIC PULL-IN current only if the magnetic stiffness at x = 0, 2F₀/g₀, is at least the
 * spring rate, i.e. x_p ≥ g₀/2; we take the marginal value x_p = g₀/2 (the contract has no preload
 * field), so that for |I| > I_p the armature leaves the rest stop with no static equilibrium short of
 * the core (snap action), and below I_p it stays put. Near x = 0 the equation is then
 * ẍ ≈ ω_n² ((I/I_p)² − 1)(x_p + x): the escape is set by ω_n and the over-drive, which makes the firing
 * current rise with dI/dt (fast 12 V ramps fire at ≈ 6.5 A, slow magneto pulses near I_p) as measured
 * (Kossor, "The Double Spark Doctrine Paradox" 2017; Cool386 coil-tester oscillograms).
 *
 * Contacts: the points stay closed while the armature takes up the cushion-spring clearance
 * (x < `breakTravel`: the upper point on its cushion spring follows the armature) and are open for
 * x > breakTravel; the coil integrator locates both crossings as events. Stops: the rest stop x = 0 and
 * the core x = `maxTravel` (point gap 1/32 in with the vibrator held down + the cushion clearance), with
 * a coefficient of restitution (UNVERIFIED default 0: inelastic — contact bounce is not modelled).
 * The armature does not feed back on the circuit (its motion's effect on L₁ is neglected), so the
 * circuit's energy ledger is unaffected.
 *
 * Integration: classical RK4 over a coil sub-step with I₁(t) linear between the sub-step ends
 * (ω_n h ≤ 0.02 for the sub-steps used); the stops clip the end state.
 */

import type { TremblerVibratorSpec } from '../core/engine-spec';

/** Allocation-free 1-DOF vibrator state + integrator. */
export class Vibrator {
  readonly spec: TremblerVibratorSpec;
  /** ω_n = 2π f_n, rad/s. */
  readonly omega: number;
  /** Preload deflection x_p = g₀/2, m. */
  readonly preload: number;
  /** Coefficient of restitution at the stops. */
  readonly restitution: number;
  /** Points-separation travel (cushion clearance), m. */
  readonly breakTravel: number;
  /** Travel to the core stop, m. */
  readonly maxTravel: number;

  /** Armature travel toward the core, m. */
  x = 0;
  /** Armature velocity, m/s. */
  v = 0;
  /** Trial end state of the last `trial` call (accepted by `commit`). */
  xNew = 0;
  vNew = 0;
  /** Largest travel since the last resetPeak, m (diagnostic). */
  peakTravel = 0;

  private readonly w2: number;
  private readonly cPull: number;
  private readonly cDamp: number;
  private readonly g0: number;
  private readonly ip: number;

  constructor(spec: TremblerVibratorSpec, restitution = 0) {
    if (!(spec.airGap > spec.maxTravel)) throw new Error('Vibrator: airGap must exceed maxTravel');
    if (!(spec.maxTravel > spec.breakTravel && spec.breakTravel > 0)) throw new Error('Vibrator: need 0 < breakTravel < maxTravel');
    this.spec = spec;
    this.omega = 2 * Math.PI * spec.naturalFrequency;
    this.preload = 0.5 * spec.airGap;
    this.restitution = restitution;
    this.breakTravel = spec.breakTravel;
    this.maxTravel = spec.maxTravel;
    this.g0 = spec.airGap;
    this.ip = spec.pullCurrent;
    this.w2 = this.omega * this.omega;
    // ω²·x_p·(I/I_p)²·(g₀/(g₀−x))² = cPull·I²/(g₀ − x)²
    this.cPull = (this.w2 * this.preload * spec.airGap * spec.airGap) / (spec.pullCurrent * spec.pullCurrent);
    this.cDamp = 2 * spec.dampingRatio * this.omega;
  }

  /** Back to rest. */
  reset(): void {
    this.x = 0;
    this.v = 0;
    this.xNew = 0;
    this.vNew = 0;
    this.peakTravel = 0;
  }

  /**
   * Armature acceleration at travel x, velocity v, primary current i, m/s². (The pull is evaluated at
   * min(x, maxTravel) so an RK stage overshooting the core stop stays finite.)
   */
  accel(x: number, v: number, i: number): number {
    const d = this.g0 - (x < this.maxTravel ? x : this.maxTravel);
    return (this.cPull * i * i) / (d * d) - this.w2 * (this.preload + x) - this.cDamp * v;
  }

  /**
   * Integrate from the current (x, v) over h with the primary current linear from i0 to i1 into
   * (xNew, vNew); the state itself is unchanged until `commit`.
   */
  trial(h: number, i0: number, i1: number): void {
    const x = this.x;
    const v = this.v;
    if (x === 0 && v === 0 && Math.abs(i0) <= this.ip && Math.abs(i1) <= this.ip) {
      // resting on the stop with the pull below the preload during the whole sub-step
      this.xNew = 0;
      this.vNew = 0;
      return;
    }
    const im = 0.5 * (i0 + i1);
    const hh = 0.5 * h;
    const a1 = this.accel(x, v, i0);
    const x2 = x + hh * v;
    const v2 = v + hh * a1;
    const a2 = this.accel(x2, v2, im);
    const x3 = x + hh * v2;
    const v3 = v + hh * a2;
    const a3 = this.accel(x3, v3, im);
    const x4 = x + h * v3;
    const v4 = v + h * a3;
    const a4 = this.accel(x4, v4, i1);
    let xn = x + (h / 6) * (v + 2 * v2 + 2 * v3 + v4);
    let vn = v + (h / 6) * (a1 + 2 * a2 + 2 * a3 + a4);
    if (xn <= 0) {
      xn = 0;
      if (vn < 0) vn = -this.restitution * vn;
    } else if (xn >= this.maxTravel) {
      xn = this.maxTravel;
      if (vn > 0) vn = -this.restitution * vn;
    }
    this.xNew = xn;
    this.vNew = vn;
  }

  /** Accept the trial state. */
  commit(): void {
    this.x = this.xNew;
    this.v = this.vNew;
    if (this.x > this.peakTravel) this.peakTravel = this.x;
  }
}
