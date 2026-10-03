/**
 * Flow-tracer particles for gas exchange (pure TypeScript, no three.js).
 *
 * Each particle stands for a fixed parcel of gas mass, so the spawn rate is
 * |ṁ_valve| / m_parcel — the particle count is the integrated mass flow.
 * Everything advances in SIMULATED time (slow motion shows true velocities).
 *
 * Paths (cylinder frame, valves open downward to y = −lift):
 *  - inflow  (intake ṁ > 0, exhaust ṁ < 0): down the port (radius D_v/2 above
 *    the seat) at Q/A_port, then out through the open arc of the curtain as a
 *    jet along the seat cone at v_jet = Q/(C_d A_curtain); the jet decays by
 *    quadratic drag into the bulk flow (swirl + axial compression) plus an
 *    Ornstein–Uhlenbeck turbulent velocity with u' and T_L = L_I/u'.
 *  - outflow (exhaust ṁ > 0, intake backflow ṁ < 0): hemispherical sink
 *    flow Q/(2π r²) toward the curtain, then up the port.
 * Particles crossing the flame sphere ∩ chamber are marked burned (and flash
 * as the front passes).
 */
import type { EngineSpec, ValveSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import {
  TRACER_CAPACITY,
  TRACER_INTAKE_LIFE_REVS,
  TRACER_MAX_DT,
  TRACER_MAX_SUBSTEPS,
  TRACER_PARCELS_PER_DISPLACED_CHARGE,
  TRACER_PORT_LENGTH_DIAMETERS,
} from './constants';
import { curtainJetSpeed, GasFlowModel, MAX_JET_SPEED } from './flow';
import { insideBurnedRegion, openCurtainArea, sampleOpenArc, type V3 } from './geometry';

export const TracerState = { PortIn: 1, InCylinder: 2, ToValve: 3, PortOut: 4 } as const;
export type TracerState = (typeof TracerState)[keyof typeof TracerState];

/** Reference charge density used to size tracer parcels, kg/m³. */
const RHO_REF = 1.15;
/** Jet decay length in valve seat diameters (quadratic drag). Approximate. */
const JET_DECAY_DIAMETERS = 2;
/** Duration of the "front passing" highlight, simulated s. */
export const FRONT_FLASH_S = 3e-4;

interface FlowSample {
  mi: number; me: number; li: number; le: number; h: number; rho: number; pim: number; pem: number; p: number; T: number;
}

function sampleFrom(s: EngineSnapshot, out: FlowSample): FlowSample {
  out.mi = s.intakeMassFlow;
  out.me = s.exhaustMassFlow;
  out.li = s.intakeLift;
  out.le = s.exhaustLift;
  out.h = s.clearanceHeight;
  out.rho = s.volume > 0 ? s.mass / s.volume : RHO_REF;
  out.pim = s.intakeManifoldPressure;
  out.pem = s.exhaustManifoldPressure;
  out.p = s.pressure;
  out.T = s.temperatureMean;
  return out;
}

function lerpSample(a: FlowSample, b: FlowSample, w: number, out: FlowSample): FlowSample {
  out.mi = a.mi + (b.mi - a.mi) * w;
  out.me = a.me + (b.me - a.me) * w;
  out.li = a.li + (b.li - a.li) * w;
  out.le = a.le + (b.le - a.le) * w;
  out.h = a.h + (b.h - a.h) * w;
  out.rho = a.rho + (b.rho - a.rho) * w;
  out.pim = a.pim + (b.pim - a.pim) * w;
  out.pem = a.pem + (b.pem - a.pem) * w;
  out.p = a.p + (b.p - a.p) * w;
  out.T = a.T + (b.T - a.T) * w;
  return out;
}

const newSample = (): FlowSample => ({ mi: 0, me: 0, li: 0, le: 0, h: 0, rho: RHO_REF, pim: 1e5, pem: 1e5, p: 1e5, T: 300 });

export class TracerSystem {
  readonly capacity: number;
  /** Alive particles occupy [0, count). */
  count = 0;
  /** Mass represented by one particle, kg. */
  readonly parcelMass: number;

  // ---- per-particle state (SoA) ----
  /** Positions, cylinder frame, m (also the render buffer). */
  readonly position: Float32Array;
  /** Jet (excess) velocity and turbulent fluctuation, m/s. */
  readonly jet: Float32Array;
  readonly turb: Float32Array;
  readonly age: Float32Array;
  readonly life: Float32Array;
  readonly state: Uint8Array;
  /** 0 = intake valve, 1 = exhaust valve. */
  readonly valve: Uint8Array;
  /** Curtain angle ψ (rad) assigned at spawn. */
  readonly psi: Float32Array;
  /** Simulated time of flame passage (NaN = unburned). */
  readonly burnTime: Float64Array;
  /** Gas temperature when spawned, K. */
  readonly spawnT: Float32Array;
  /** 1 if the parcel is hot residual/exhaust gas (vs fresh charge). */
  readonly hot: Uint8Array;

  /** Simulated time, s. */
  t = 0;

  private readonly R: number;
  private readonly valves: [ValveSpec, ValveSpec];
  /** Length of the port section above each valve seat (intake, exhaust), m. */
  readonly portLen: [number, number];
  private readonly accum = [0, 0];
  private prev: FlowSample = newSample();
  private cur: FlowSample = newSample();
  private readonly mid: FlowSample = newSample();
  private hasPrev = false;
  private rng: number;
  private readonly vtmp: V3 = [0, 0, 0];

  constructor(spec: EngineSpec, capacity = TRACER_CAPACITY, seed = 12345) {
    this.capacity = capacity;
    this.R = spec.geometry.bore / 2;
    this.valves = [spec.intakeValve, spec.exhaustValve];
    this.portLen = [
      TRACER_PORT_LENGTH_DIAMETERS * spec.intakeValve.seatInnerDiameter,
      TRACER_PORT_LENGTH_DIAMETERS * spec.exhaustValve.seatInnerDiameter,
    ];
    const Vd = (Math.PI / 4) * spec.geometry.bore ** 2 * spec.geometry.stroke;
    this.parcelMass = (RHO_REF * Vd) / TRACER_PARCELS_PER_DISPLACED_CHARGE;
    this.position = new Float32Array(3 * capacity);
    this.jet = new Float32Array(3 * capacity);
    this.turb = new Float32Array(3 * capacity);
    this.age = new Float32Array(capacity);
    this.life = new Float32Array(capacity);
    this.state = new Uint8Array(capacity);
    this.valve = new Uint8Array(capacity);
    this.psi = new Float32Array(capacity);
    this.burnTime = new Float64Array(capacity);
    this.spawnT = new Float32Array(capacity);
    this.hot = new Uint8Array(capacity);
    this.rng = seed >>> 0;
  }

  reset(): void {
    this.count = 0;
    this.accum[0] = this.accum[1] = 0;
    this.hasPrev = false;
  }

  private rand(): number {
    let t = (this.rng = (this.rng + 0x6d2b79f5) | 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  private gauss(): number {
    const u = Math.max(this.rand(), 1e-12), v = this.rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  /**
   * Advance to snapshot `s` (dtSim simulated seconds after the previous call).
   * `flow` must already be updated to `s`. `flame*` describe the burned sphere
   * (radius 0 = none); `allBurned` marks every in-cylinder parcel burned.
   */
  update(
    s: EngineSnapshot, dtSim: number, flow: GasFlowModel,
    flameCenter: ArrayLike<number>, flameRadius: number, allBurned: boolean,
  ): void {
    if (!this.hasPrev || dtSim <= 0) {
      sampleFrom(s, this.cur);
      if (!this.hasPrev) sampleFrom(s, this.prev);
      this.hasPrev = true;
      this.t = s.t;
      return;
    }
    // rotate samples
    const tmp = this.prev;
    this.prev = this.cur;
    this.cur = sampleFrom(s, tmp);

    const n = Math.min(TRACER_MAX_SUBSTEPS, Math.max(1, Math.ceil(dtSim / TRACER_MAX_DT)));
    const dt = dtSim / n;
    const rpm = Math.max(Math.abs(s.rpm), 1);
    const intakeLife = (TRACER_INTAKE_LIFE_REVS * 60) / rpm;
    const t0 = this.t;
    for (let k = 0; k < n; k++) {
      const smp = lerpSample(this.prev, this.cur, (k + 0.5) / n, this.mid);
      // walls are enforced at the END of the sub-step (piston position then)
      smp.h = this.prev.h + (this.cur.h - this.prev.h) * ((k + 1) / n);
      this.spawn(0, smp, dt, intakeLife);
      this.spawn(1, smp, dt, intakeLife);
      this.step(smp, dt, flow, flameCenter, flameRadius, allBurned, t0 + (k + 1) * dt);
    }
    this.t = s.t;
  }

  // -------------------------------------------------------------------------

  private spawn(vi: 0 | 1, smp: FlowSample, dt: number, intakeLife: number): void {
    const v = this.valves[vi];
    const lift = vi === 0 ? smp.li : smp.le;
    const mdot = vi === 0 ? smp.mi : smp.me;
    if (!(lift > 1e-6) || mdot === 0) return;
    // inflow into the cylinder: intake ṁ > 0 or exhaust ṁ < 0
    const inflow = vi === 0 ? mdot > 0 : mdot < 0;
    this.accum[vi] += (Math.abs(mdot) * dt) / this.parcelMass;
    let nNew = Math.floor(this.accum[vi]);
    this.accum[vi] -= nNew;
    const a = v.seatInnerDiameter / 2;
    const [vx, vz] = v.position;
    const portLen = this.portLen[vi];
    const P = this.position;
    while (nNew-- > 0 && this.count < this.capacity) {
      const i = this.count++;
      const i3 = 3 * i;
      this.age[i] = this.rand() * dt; // spread births across the sub-step
      this.valve[i] = vi;
      this.burnTime[i] = NaN;
      this.jet[i3] = this.jet[i3 + 1] = this.jet[i3 + 2] = 0;
      this.turb[i3] = this.turb[i3 + 1] = this.turb[i3 + 2] = 0;
      this.spawnT[i] = smp.T;
      if (inflow) {
        // enters from the top of the port above its curtain angle
        const psi = vi === 0 ? sampleOpenArc(v, this.rand()) : 2 * Math.PI * this.rand();
        this.psi[i] = psi;
        const rr = a * (0.55 + 0.4 * this.rand());
        P[i3] = vx + rr * Math.cos(psi);
        P[i3 + 1] = portLen * (0.9 + 0.1 * this.rand());
        P[i3 + 2] = vz + rr * Math.sin(psi);
        this.state[i] = TracerState.PortIn;
        this.life[i] = intakeLife;
        // intake charge is fresh; exhaust backflow is hot residual
        this.hot[i] = vi === 1 ? 1 : 0;
      } else {
        // outflow: start on a lower half-shell around the curtain centre
        const r0 = a * (1.3 + 1.4 * this.rand());
        const th = 2 * Math.PI * this.rand();
        const cz = this.rand(); // cos of polar angle from −y
        const sz = Math.sqrt(1 - cz * cz);
        let x = vx + r0 * sz * Math.cos(th);
        let z = vz + r0 * sz * Math.sin(th);
        let y = -lift / 2 - r0 * cz;
        const rxz = Math.hypot(x, z);
        const lim = 0.97 * this.R;
        if (rxz > lim) { x *= lim / rxz; z *= lim / rxz; }
        if (y < -smp.h) y = -smp.h * (0.2 + 0.7 * this.rand());
        if (y > 0) y = 0;
        P[i3] = x; P[i3 + 1] = y; P[i3 + 2] = z;
        this.psi[i] = Math.atan2(z - vz, x - vx);
        this.state[i] = TracerState.ToValve;
        this.life[i] = 0.5 * intakeLife;
        this.hot[i] = smp.T > 600 ? 1 : 0;
      }
    }
  }

  private kill(i: number): void {
    const j = --this.count;
    if (i === j) return;
    const i3 = 3 * i, j3 = 3 * j;
    for (let c = 0; c < 3; c++) {
      this.position[i3 + c] = this.position[j3 + c];
      this.jet[i3 + c] = this.jet[j3 + c];
      this.turb[i3 + c] = this.turb[j3 + c];
    }
    this.age[i] = this.age[j];
    this.life[i] = this.life[j];
    this.state[i] = this.state[j];
    this.valve[i] = this.valve[j];
    this.psi[i] = this.psi[j];
    this.burnTime[i] = this.burnTime[j];
    this.spawnT[i] = this.spawnT[j];
    this.hot[i] = this.hot[j];
  }

  private step(
    smp: FlowSample, dt: number, flow: GasFlowModel,
    fc: ArrayLike<number>, fr: number, allBurned: boolean, tNow: number,
  ): void {
    const P = this.position;
    const R = this.R;
    const h = smp.h;
    const up = flow.uPrime;
    const TL = up > 0 ? flow.integralScale / up : 1;
    const eT = Math.exp(-dt / TL);
    const sqT = up * Math.sqrt(Math.max(0, 1 - eT * eT));
    const v = this.vtmp;

    // per-valve flow quantities for this sub-step
    for (let i = 0; i < this.count; ) {
      const i3 = 3 * i;
      const vi = this.valve[i];
      const valve = this.valves[vi];
      const lift = vi === 0 ? smp.li : smp.le;
      const mdot = vi === 0 ? smp.mi : smp.me;
      const inflow = vi === 0 ? mdot > 0 : mdot < 0;
      const a = valve.seatInnerDiameter / 2;
      const [vx, vz] = valve.position;
      // port density: manifold pressure at cylinder temperature (approximate)
      const pm = vi === 0 ? smp.pim : smp.pem;
      const rhoPort = smp.p > 0 ? smp.rho * (pm / smp.p) : smp.rho;
      const Q = Math.abs(mdot) / Math.max(inflow ? rhoPort : smp.rho, 1e-3);
      const vPort = Math.min(Q / (Math.PI * a * a), MAX_JET_SPEED);
      let dead = false;

      this.age[i] += dt;
      switch (this.state[i]) {
        case TracerState.PortIn: {
          if (lift > 1e-6 && inflow) P[i3 + 1] -= vPort * dt;
          if (P[i3 + 1] <= 0) {
            // exit through the curtain as a jet along the seat cone
            const psi = this.psi[i];
            const c = Math.cos(psi), sn = Math.sin(psi);
            P[i3] = vx + a * c;
            P[i3 + 1] = -Math.max(lift, 1e-5) * this.rand();
            P[i3 + 2] = vz + a * sn;
            const vj = curtainJetSpeed(mdot, inflow ? rhoPort : smp.rho, openCurtainArea(valve, lift));
            const cb = Math.cos(valve.seatAngle), sb = Math.sin(valve.seatAngle);
            this.jet[i3] = vj * cb * c;
            this.jet[i3 + 1] = -vj * sb;
            this.jet[i3 + 2] = vj * cb * sn;
            this.state[i] = TracerState.InCylinder;
          }
          break;
        }
        case TracerState.ToValve: {
          if (!(lift > 1e-6) || inflow) {
            this.state[i] = TracerState.InCylinder;
            break;
          }
          const dx = vx - P[i3], dy = -lift / 2 - P[i3 + 1], dz = vz - P[i3 + 2];
          const r = Math.hypot(dx, dy, dz);
          if (r < 0.9 * a) {
            this.state[i] = TracerState.PortOut;
            P[i3 + 1] = 0;
            break;
          }
          const sp = Math.min(Q / (2 * Math.PI * r * r), MAX_JET_SPEED);
          const stepLen = Math.min(sp * dt, r);
          P[i3] += (dx / r) * stepLen;
          P[i3 + 1] += (dy / r) * stepLen;
          P[i3 + 2] += (dz / r) * stepLen;
          break;
        }
        case TracerState.PortOut: {
          P[i3 + 1] += vPort * dt;
          if (P[i3 + 1] > this.portLen[vi]) dead = true;
          break;
        }
        case TracerState.InCylinder: {
          const x = P[i3], y = P[i3 + 1], z = P[i3 + 2];
          flow.velocityAt(x, y, z, v);
          // jet decay: de/dt = −|e| e / (K D_v)
          const jx = this.jet[i3], jy = this.jet[i3 + 1], jz = this.jet[i3 + 2];
          const jm = Math.hypot(jx, jy, jz);
          if (jm > 0) {
            const f = 1 / (1 + (jm * dt) / (JET_DECAY_DIAMETERS * 2 * a));
            this.jet[i3] = jx * f; this.jet[i3 + 1] = jy * f; this.jet[i3 + 2] = jz * f;
          }
          // turbulence (Ornstein–Uhlenbeck)
          this.turb[i3] = this.turb[i3] * eT + sqT * this.gauss();
          this.turb[i3 + 1] = this.turb[i3 + 1] * eT + sqT * this.gauss();
          this.turb[i3 + 2] = this.turb[i3 + 2] * eT + sqT * this.gauss();
          let ux = v[0] + this.jet[i3] + this.turb[i3];
          let uy = v[1] + this.jet[i3 + 1] + this.turb[i3 + 1];
          let uz = v[2] + this.jet[i3 + 2] + this.turb[i3 + 2];
          // sinks: any valve with outflow draws the gas toward its curtain
          for (let w = 0; w < 2; w++) {
            const vw = this.valves[w];
            const lw = w === 0 ? smp.li : smp.le;
            const mw = w === 0 ? smp.mi : smp.me;
            const outw = w === 0 ? mw < 0 : mw > 0;
            if (!outw || !(lw > 1e-6)) continue;
            const aw = vw.seatInnerDiameter / 2;
            const dx = vw.position[0] - x, dy = -lw / 2 - y, dz = vw.position[1] - z;
            const r = Math.hypot(dx, dy, dz);
            if (r < 0.9 * aw) {
              this.state[i] = TracerState.PortOut;
              this.valve[i] = w;
              P[i3 + 1] = 0;
              break;
            }
            const sp = Math.min(Math.abs(mw) / Math.max(smp.rho, 1e-3) / (2 * Math.PI * r * r), MAX_JET_SPEED);
            ux += (dx / r) * sp; uy += (dy / r) * sp; uz += (dz / r) * sp;
          }
          if (this.state[i] !== TracerState.InCylinder) break;
          let nx = x + ux * dt, ny = y + uy * dt, nz = z + uz * dt;
          // walls: liner, head, piston (reflect the non-bulk velocities)
          const rxz = Math.hypot(nx, nz);
          const lim = 0.995 * R;
          if (rxz > lim) {
            const ex = nx / rxz, ez = nz / rxz;
            nx = ex * lim; nz = ez * lim;
            const jr = this.jet[i3] * ex + this.jet[i3 + 2] * ez;
            if (jr > 0) { this.jet[i3] -= 2 * jr * ex; this.jet[i3 + 2] -= 2 * jr * ez; }
            const tr = this.turb[i3] * ex + this.turb[i3 + 2] * ez;
            if (tr > 0) { this.turb[i3] -= 2 * tr * ex; this.turb[i3 + 2] -= 2 * tr * ez; }
          }
          if (ny > 0) {
            ny = 0;
            if (this.jet[i3 + 1] > 0) this.jet[i3 + 1] = -this.jet[i3 + 1];
            if (this.turb[i3 + 1] > 0) this.turb[i3 + 1] = -this.turb[i3 + 1];
          } else if (ny < -h) {
            ny = -h;
            if (this.jet[i3 + 1] < 0) this.jet[i3 + 1] = -this.jet[i3 + 1];
            if (this.turb[i3 + 1] < 0) this.turb[i3 + 1] = -this.turb[i3 + 1];
          }
          P[i3] = nx; P[i3 + 1] = ny; P[i3 + 2] = nz;
          break;
        }
      }
      // flame passage (any parcel in the chamber)
      const st = this.state[i];
      if ((st === TracerState.InCylinder || st === TracerState.ToValve) && Number.isNaN(this.burnTime[i])) {
        if (allBurned) this.burnTime[i] = tNow - 1;
        else if (fr > 0 && insideBurnedRegion(P[i3], P[i3 + 1], P[i3 + 2], R, h, fc, fr)) this.burnTime[i] = tNow;
      }
      if (dead || this.age[i] > this.life[i]) {
        this.kill(i);
        continue;
      }
      i++;
    }
  }
}
