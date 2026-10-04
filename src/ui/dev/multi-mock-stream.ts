/**
 * DEV ONLY — a multi-cylinder synthetic stream for the UI harness (src/ui/dev/ui-harness.html
 * ?engine=ford-model-t), so the 4-cylinder HUD, charts and Cycles tab can be exercised without the
 * physics. It is NOT a physics model: one MockStream cartoon per cylinder, started at that cylinder's
 * local angle (spec.layout.firingOffsetDeg) and stepped in lock-step, combined into multi-cylinder
 * snapshots (top level = cylinder 1 + engine fields, `cylinders[]` = every cylinder); for a
 * trembler-magneto spec each cylinder's spark is replaced by a cartoon shower on its timer contact, and
 * cylinder 1's summaries carry a cartoon EngineCycleSummary. Never import it from production code.
 */
import { firingOffsetsDeg, type EngineSpec } from '../../physics/core/engine-spec';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { CycleSummary, CylinderSnapshot, EngineCycleSummary, EngineSnapshot, SparkPhase } from '../../physics/core/snapshot';
import { wrapCycleDeg } from '../engine-cycle';
import { MockStream } from './mock-stream';

/** Cartoon trembler timing (shape only): first fire after the coil build-up, then a ≈ 500 Hz buzz. */
const FIRST_FIRE_S = 3.5e-3;
const BUZZ_PERIOD_S = 2e-3;

export class MultiMockStream {
  t = 0;
  private readonly cores: MockStream[];
  private readonly offsets: readonly number[];
  private op: OperatingPoint;
  private readonly cycles: CycleSummary[] = [];

  constructor(
    private readonly spec: EngineSpec,
    op: OperatingPoint,
  ) {
    this.op = { ...op };
    this.offsets = firingOffsetsDeg(spec);
    this.cores = this.offsets.map(() => new MockStream(spec, this.op));
    this.startCores();
  }

  /** Engine crank angle (cylinder 1's), deg. */
  get theta(): number {
    return this.cores[0].theta;
  }

  /** Engine cycle (cylinder 1's). */
  get cycle(): number {
    return this.cores[0].cycle;
  }

  setOperatingPoint(p: Partial<OperatingPoint>): void {
    this.op = { ...this.op, ...p };
    for (const c of this.cores) c.setOperatingPoint(p);
  }

  reset(): void {
    this.t = 0;
    for (const c of this.cores) c.reset();
    this.cycles.length = 0;
    this.startCores();
  }

  drainCycles(): CycleSummary[] {
    return this.cycles.splice(0);
  }

  /** True while any cylinder's ignition timer contact is closed (the harness samples densely then). */
  inSparkWindow(): boolean {
    for (let i = 0; i < this.cores.length; i++) if (this.timerPhase(this.cores[i].theta) >= 0) return true;
    return false;
  }

  /** Advance by dθ degrees and return the engine snapshot. */
  step(dDeg: number): EngineSnapshot {
    const snaps = this.cores.map((c) => c.step(dDeg));
    this.t += dDeg / (6 * this.op.rpm);
    const cyl: CylinderSnapshot[] = snaps.map((s, i) => {
      // Engine-level fields stay at the top level only.
      const { t: _t, rpm: _rpm, intakeManifoldPressure: _im, exhaustManifoldPressure: _em, netTorque: _nt, ...rest } = s;
      const c: CylinderSnapshot = { ...rest, index: i };
      if (this.spec.ignition.type === 'trembler-magneto') c.spark = this.trembler(c.thetaDeg, s.spark.energyDelivered);
      return c;
    });
    const top: EngineSnapshot = { ...snaps[0], ...cyl[0], t: this.t, cylinders: cyl };
    delete (top as Partial<CylinderSnapshot>).index;
    let firing = -1;
    cyl.forEach((c, i) => {
      if (c.spark.timerClosed) firing = i;
    });
    top.firingCylinder = firing;
    const ig = this.spec.ignition;
    if (ig.type === 'trembler-magneto') {
      const w = (2 * Math.PI * this.op.rpm) / 60;
      const crank = ((top.thetaDeg % 360) + 360) % 360;
      top.magnetoEmf = ig.magneto.emfConstant * w * Math.sin((ig.magneto.cyclesPerRevolution * (crank - ig.magneto.phaseDeg) * Math.PI) / 180);
    }
    const v = this.vehicleSpeed();
    if (v !== undefined) top.vehicleSpeed = v;
    top.loadTorque = this.roadTorque(v);
    top.frictionTorque = 12;
    for (let i = 0; i < this.cores.length; i++) {
      for (const c of this.cores[i].drainCycles()) {
        c.cylinder = i;
        if (i === 0) c.engine = this.engineSummary(c);
        this.cycles.push(c);
      }
    }
    return top;
  }

  // ------------------------------------------------------------------

  /** Pre-advance each core to its local angle at engine θ = −360. */
  private startCores(): void {
    this.cores.forEach((c, i) => {
      const th0 = wrapCycleDeg(-360 - this.offsets[i]);
      if (th0 > -360) c.step(th0 + 360);
    });
    this.t = 0;
  }

  /** Time since the timer contact closed, s (−1 outside the contact). */
  private timerPhase(theta: number): number {
    const ig = this.spec.ignition;
    if (ig.type !== 'trembler-magneto') return -1;
    const make = -this.op.sparkAdvanceDeg;
    let d = theta - make;
    if (d < 0) d += 720;
    if (d >= ig.timer.contactArcDeg) return -1;
    return d / (6 * this.op.rpm);
  }

  private trembler(theta: number, energy: number): EngineSnapshot['spark'] {
    const tau = this.timerPhase(theta);
    let phase: SparkPhase = 'off';
    let count = 0;
    let iPrim = 0;
    let vSec = 0;
    let iSec = 0;
    if (tau >= 0) {
      if (tau < FIRST_FIRE_S) {
        phase = 'charging';
        iPrim = 5 * (tau / FIRST_FIRE_S);
      } else {
        const k = Math.floor((tau - FIRST_FIRE_S) / BUZZ_PERIOD_S);
        const ph = tau - FIRST_FIRE_S - k * BUZZ_PERIOD_S;
        count = k + 1;
        if (ph < 10e-6) {
          phase = 'breakdown';
          vSec = 7e3;
          iSec = 50;
        } else if (ph < 60e-6) {
          phase = 'arc';
          vSec = 100;
          iSec = 0.15;
        } else if (ph < 400e-6) {
          phase = 'glow';
          vSec = 500;
          iSec = 0.04 * (1 - (ph - 60e-6) / 340e-6);
        } else {
          phase = 'charging';
          iPrim = 5 * ((ph - 400e-6) / (BUZZ_PERIOD_S - 400e-6));
        }
      }
    } else if (wrapCycleDeg(theta + this.op.sparkAdvanceDeg) > 0 && wrapCycleDeg(theta + this.op.sparkAdvanceDeg) < 180) {
      phase = 'done';
    }
    const make = -this.op.sparkAdvanceDeg;
    const first = tau >= FIRST_FIRE_S ? make + FIRST_FIRE_S * 6 * this.op.rpm : NaN;
    return {
      phase,
      primaryCurrent: iPrim,
      secondaryVoltage: vSec,
      secondaryCurrent: iSec,
      energyDelivered: count > 0 ? Math.max(energy, 2e-3 * count) : 0,
      breakdownVoltage: 6e3,
      breakdownCount: tau >= 0 ? count : 0,
      pointsOpen: phase !== 'charging',
      timerClosed: tau >= 0,
      firstSparkDeg: first,
      primaryVoltage: phase === 'breakdown' ? 300 : 0,
    };
  }

  private vehicleSpeed(): number | undefined {
    const load = this.op.load;
    const veh = this.spec.vehicle;
    if (!veh || this.op.speedMode !== 'free' || load?.kind !== 'vehicle' || load.gear === 'neutral') return undefined;
    const ratio = veh.gears[load.gear] ?? 1;
    return ((2 * Math.PI * this.op.rpm) / 60 / (ratio * veh.finalDrive)) * veh.wheelRadius;
  }

  /** Road load reflected to the crank (cartoon of the 'vehicle' model), N m. */
  private roadTorque(v: number | undefined): number {
    const veh = this.spec.vehicle;
    const load = this.op.load;
    if (!veh || v === undefined || load?.kind !== 'vehicle') return this.op.loadTorque;
    const g = 9.80665;
    const f = veh.mass * g * (veh.rollingResistance + load.grade) + 0.5 * 1.2 * veh.dragArea * v * v;
    const ratio = (veh.gears[load.gear] ?? 1) * veh.finalDrive;
    return (f * veh.wheelRadius) / (ratio * veh.drivelineEfficiency);
  }

  private engineSummary(c: CycleSummary): EngineCycleSummary {
    const g = this.spec.geometry;
    const vd = (Math.PI / 4) * g.bore * g.bore * g.stroke * this.spec.cylinders;
    const w = (2 * Math.PI * this.op.rpm) / 60;
    const indicated = (c.imepNet * vd) / (4 * Math.PI);
    const fmep = 1e5;
    const friction = (fmep * vd) / (4 * Math.PI);
    const brake = indicated - friction;
    return {
      rpmMean: this.op.rpm,
      indicatedTorque: indicated,
      frictionTorque: friction,
      brakeTorque: brake,
      brakePower: brake * w,
      bmep: (4 * Math.PI * brake) / vd,
      imepNet: c.imepNet,
      fmep,
      airMassFlow: 0,
      fuelMassFlow: (c.fuelMass * this.spec.cylinders * this.op.rpm) / 120,
      volumetricEfficiency: c.volumetricEfficiency,
      bsfc: (c.fuelMass * this.spec.cylinders * this.op.rpm) / 120 / Math.max(1, brake * w),
      brakeEfficiency: c.indicatedEfficiency * (brake / Math.max(1e-9, indicated)),
      loadTorque: brake,
      vehicleSpeed: this.vehicleSpeed(),
    };
  }
}
