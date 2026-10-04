/**
 * Primary supply of a trembler-coil ignition: a storage battery (DC EMF behind its internal
 * resistance) or the low-tension AC flywheel magneto, SI units (V, Ω, H, s, rad/s; angles in crank
 * degrees).
 *
 * Magneto (Ford Model T: 16 V-magnets on the flywheel sweeping 16 series coils on a stationary ring,
 * Ford Service Par. 281–283; 8 electrical cycles per crank revolution, Patterson & Coniff 2003) as an
 * EMF behind a series source impedance (core/engine-spec.ts MagnetoSpec):
 *     e(t) = k · ω(t) · sin(N · (θ(t) − φ) · π/180),       Z_s = R_s + jω_e L_s,
 * θ the ENGINE crank angle (cylinder 1; the EMF is periodic in 360°/N so θ mod 360° is implied),
 * ω the crank speed, N = cyclesPerRevolution, φ = phaseDeg (an EMF zero crossing, rising). The open-
 * circuit amplitude ∝ ω is Faraday's law for a fixed flux pattern; the source inductance L_s makes the
 * short-circuit current tend to k/(N L_s) at speed — the "constant current" Ford's test-stand pamphlet
 * describes (see engines/model-t.ts MODEL_T_IGNITION for the sourced constants).
 *
 * Within one caller step the crank angle and speed are taken LINEAR in time between the step-start and
 * step-end values the ignition system receives (`setStep`); the coil integrator samples e(t) at its own
 * sub-step ends (trapezoidal mean, which keeps its energy ledger exact for any e(t)).
 *
 * The source impedance is part of the connected coil's primary loop. A Ford timer grounds one coil at a
 * time (contact 87° < 180° between cylinders), so lumping R_s, L_s into the connected loop is exact for
 * a shared magneto; the selected supply (BAT/MAG) may only change while no current flows (`request` →
 * `applyRequested`, called by the coil at the timer make).
 */

import type { MagnetoSpec } from '../core/engine-spec';

/** Selected supply (OperatingPoint.ignitionSource). */
export type SupplyKind = 'magneto' | 'battery';

/** Battery data (TremblerMagnetoIgnitionSpec.battery). */
export interface BatterySpec {
  /** Open-circuit voltage, V. */
  voltage: number;
  /** Internal resistance, Ω. */
  internalResistance: number;
}

const DEG = Math.PI / 180;

/**
 * Open-circuit magneto EMF e = k ω sin(N (θ − φ) π/180), V, at ENGINE crank angle θ (deg) and crank
 * speed ω (rad/s).
 */
export function magnetoEmf(m: MagnetoSpec, thetaEngineDeg: number, omega: number): number {
  return m.emfConstant * omega * Math.sin(m.cyclesPerRevolution * (thetaEngineDeg - m.phaseDeg) * DEG);
}

/**
 * Switchable primary supply (battery and/or magneto) with the per-step angle/speed interpolation the
 * coil integrator samples. Allocation-free.
 */
export class PrimarySupply {
  readonly battery: BatterySpec | null;
  readonly magneto: MagnetoSpec | null;
  /** Supply in use (fixed while a coil is connected). */
  kind: SupplyKind;
  /** Supply selected by the driver; becomes `kind` at the next `applyRequested`. */
  requested: SupplyKind;

  // ---- step interpolation: θ(t) = th0 + dth·f, ω(t) = w0 + dw·f, f = (t − t0)/dt ----
  private t0 = 0;
  private invDt = 0;
  private th0 = 0;
  private dth = 0;
  private w0 = 0;
  private dw = 0;

  constructor(battery: BatterySpec | null, magneto: MagnetoSpec | null, kind: SupplyKind) {
    this.battery = battery;
    this.magneto = magneto;
    if (kind === 'magneto' && !magneto) throw new Error('PrimarySupply: no magneto');
    if (kind === 'battery' && !battery) throw new Error('PrimarySupply: no battery');
    this.kind = kind;
    this.requested = kind;
  }

  /** DC bench supply: a battery of EMF `voltage` behind `resistance` (Ω). */
  static dc(voltage: number, resistance = 0): PrimarySupply {
    return new PrimarySupply({ voltage, internalResistance: resistance }, null, 'battery');
  }

  /** Series source resistance of the supply in use, Ω. */
  get resistance(): number {
    return this.kind === 'magneto' ? this.magneto!.internalResistance : this.battery!.internalResistance;
  }

  /** Series source inductance of the supply in use, H (0 for the battery). */
  get inductance(): number {
    return this.kind === 'magneto' ? this.magneto!.internalInductance : 0;
  }

  /** Ask for another supply (takes effect at the next applyRequested, i.e. the next timer make). */
  request(kind: SupplyKind): void {
    if (kind === 'magneto' && !this.magneto) throw new Error('PrimarySupply: no magneto');
    if (kind === 'battery' && !this.battery) throw new Error('PrimarySupply: no battery');
    this.requested = kind;
  }

  /** Switch to the requested supply (call only while the primary carries no current). */
  applyRequested(): void {
    this.kind = this.requested;
  }

  /**
   * Engine angle and speed over the caller step [t0, t0 + dt] (linear in time). theta1 must be the
   * step-end angle UNWRAPPED forward from theta0 (deg); omega in rad/s.
   */
  setStep(t0: number, dt: number, theta0: number, theta1: number, omega0: number, omega1: number): void {
    this.t0 = t0;
    this.invDt = dt > 0 ? 1 / dt : 0;
    this.th0 = theta0;
    this.dth = theta1 - theta0;
    this.w0 = omega0;
    this.dw = omega1 - omega0;
  }

  /** Engine angle at time t from the current step interpolation, deg. */
  thetaAt(t: number): number {
    return this.th0 + this.dth * (t - this.t0) * this.invDt;
  }

  /** Crank speed at time t from the current step interpolation, rad/s. */
  omegaAt(t: number): number {
    return this.w0 + this.dw * (t - this.t0) * this.invDt;
  }

  /** EMF of the supply in use at (ignition-clock) time t, V. */
  emf(t: number): number {
    if (this.kind === 'battery') return this.battery!.voltage;
    const f = (t - this.t0) * this.invDt;
    const m = this.magneto!;
    return m.emfConstant * (this.w0 + this.dw * f) * Math.sin(m.cyclesPerRevolution * (this.th0 + this.dth * f - m.phaseDeg) * DEG);
  }

  /** Open-circuit magneto EMF at time t (0 without a magneto), V. */
  magnetoEmfAt(t: number): number {
    if (!this.magneto) return 0;
    return magnetoEmf(this.magneto, this.thetaAt(t), this.omegaAt(t));
  }
}
