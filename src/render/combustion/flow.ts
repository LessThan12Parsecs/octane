/**
 * Bulk in-cylinder gas motion for the tracer particles and the spark-channel
 * convection. The snapshot does not carry a velocity field, so this is a
 * small DERIVED 0D model driven only by snapshot fields + geometry:
 *
 *  - axial: uniform compression of a flat disc → v_y(y) = y · ḣ/h
 *    (0 at the head, piston speed at the crown), ḣ from successive snapshots;
 *  - swirl: solid-body rotation ω about the cylinder axis from an angular
 *    momentum balance  dH/dt = η_s Σ ṁ_in · v_h · arm − H(ṁ_out/m + 1/τ_s),
 *    with the jet speed from the curtain flow area and the moment arm from
 *    the shroud orientation (geometry.swirlArm); I = m R²/2;
 *  - turbulence intensity u' from the snapshot (or ½ S̄p fallback, Heywood 1988).
 *
 * Pure TypeScript (no three.js).
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import {
  CURTAIN_DISCHARGE_COEFF,
  INTEGRAL_SCALE_FRACTION,
  MAX_SWIRL_RATIO,
  SWIRL_DECAY_REVS,
  SWIRL_TRANSFER_EFFICIENCY,
  UPRIME_OVER_MEAN_PISTON_SPEED,
} from './constants';
import { openCurtainArea, swirlArm, type V3 } from './geometry';

/** Upper bound on estimated jet speeds, m/s (≈ sonic at intake conditions). */
export const MAX_JET_SPEED = 450;

/** Curtain jet speed estimate for a valve mass flow, m/s. */
export function curtainJetSpeed(massFlow: number, density: number, curtainArea: number): number {
  if (!(curtainArea > 1e-9) || !(density > 1e-6)) return 0;
  const v = Math.abs(massFlow) / (density * CURTAIN_DISCHARGE_COEFF * curtainArea);
  return Math.min(v, MAX_JET_SPEED);
}

export class GasFlowModel {
  readonly R: number;
  readonly stroke: number;
  private readonly armIn: number;
  private readonly armEx: number;
  private readonly cosSeatIn: number;
  private readonly cosSeatEx: number;

  /** Angular momentum of the charge about +y, kg m²/s. */
  H = 0;
  /** Solid-body swirl rate about +y, rad/s. */
  omega = 0;
  /** Clearance height h and its rate, m, m/s. */
  h = 0;
  hdot = 0;
  /** Turbulence intensity u', m/s. */
  uPrime = 0;
  /** Integral scale L_I, m. */
  integralScale = 1e-3;
  /** Cylinder gas density, kg/m³. */
  density = 1;
  /** Intake / exhaust curtain jet speeds (magnitude), m/s. */
  intakeJetSpeed = 0;
  exhaustJetSpeed = 0;
  private hasPrev = false;

  constructor(private readonly spec: EngineSpec) {
    this.R = spec.geometry.bore / 2;
    this.stroke = spec.geometry.stroke;
    this.armIn = swirlArm(spec.intakeValve);
    this.armEx = swirlArm(spec.exhaustValve);
    this.cosSeatIn = Math.cos(spec.intakeValve.seatAngle);
    this.cosSeatEx = Math.cos(spec.exhaustValve.seatAngle);
  }

  reset(): void {
    this.H = 0;
    this.omega = 0;
    this.hdot = 0;
    this.hasPrev = false;
  }

  /** Advance by dtSim simulated seconds to snapshot s. */
  update(s: EngineSnapshot, dtSim: number): void {
    const m = s.mass > 0 ? s.mass : 1e-4;
    this.density = s.volume > 0 ? m / s.volume : 1;
    if (this.hasPrev && dtSim > 0) this.hdot = (s.clearanceHeight - this.h) / dtSim;
    else if (!this.hasPrev) this.hdot = 0;
    this.h = s.clearanceHeight;
    this.hasPrev = true;

    const iv = this.spec.intakeValve, ev = this.spec.exhaustValve;
    this.intakeJetSpeed = curtainJetSpeed(s.intakeMassFlow, this.density, openCurtainArea(iv, s.intakeLift));
    this.exhaustJetSpeed = curtainJetSpeed(s.exhaustMassFlow, this.density, openCurtainArea(ev, s.exhaustLift));

    // Angular-momentum balance (inflow through either valve carries its jets' moment).
    if (dtSim > 0) {
      let source = 0;
      let outflow = 0;
      if (s.intakeMassFlow > 0) source += s.intakeMassFlow * this.intakeJetSpeed * this.cosSeatIn * this.armIn;
      else outflow += -s.intakeMassFlow;
      if (s.exhaustMassFlow < 0) source += -s.exhaustMassFlow * this.exhaustJetSpeed * this.cosSeatEx * this.armEx;
      else outflow += s.exhaustMassFlow;
      source *= SWIRL_TRANSFER_EFFICIENCY;
      const rpm = Math.max(Math.abs(s.rpm), 1);
      const tau = (SWIRL_DECAY_REVS * 60) / rpm;
      const k = outflow / m + 1 / tau;
      const e = Math.exp(-k * dtSim);
      this.H = this.H * e + (k > 1e-12 ? (source * (1 - e)) / k : source * dtSim);
    }
    const I = 0.5 * m * this.R * this.R;
    const wMax = (MAX_SWIRL_RATIO * 2 * Math.PI * Math.max(Math.abs(s.rpm), 1)) / 60;
    if (Math.abs(this.H) > wMax * I) this.H = Math.sign(this.H) * wMax * I;
    this.omega = this.H / I;

    this.integralScale = INTEGRAL_SCALE_FRACTION * Math.max(s.clearanceHeight, 1e-4);
    const up = s.flame.turbulenceIntensity;
    this.uPrime = up > 0
      ? up
      : UPRIME_OVER_MEAN_PISTON_SPEED * (2 * this.stroke * Math.abs(s.rpm)) / 60;
  }

  /** Bulk gas velocity (swirl + axial compression) at a cylinder-frame point, m/s. */
  velocityAt(x: number, y: number, z: number, out: V3): V3 {
    const w = this.omega;
    out[0] = w * z;
    out[2] = -w * x;
    out[1] = this.h > 1e-6 ? (y * this.hdot) / this.h : 0;
    return out;
  }
}
