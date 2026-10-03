/**
 * DEV ONLY — a crude synthetic snapshot source used by the UI harness
 * (src/ui/dev/ui-harness.html) to exercise the charts/HUD before the real physics
 * simulator exists. It is NOT a physics model: single-zone ideal gas with a Wiebe
 * burn, sinusoidal valve lifts, a cartoon ignition coil and a Douaud–Eyzat knock
 * integral with a damped 1st-circumferential-mode oscillation. Never import it
 * from production code.
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot, SparkPhase } from '../../physics/core/snapshot';
import { valveOpenInterval, wrapCycleDeg } from '../engine-cycle';

const R_AIR = 287;
const GAMMA_U = 1.33;
const GAMMA_B = 1.25;
const LHV = 44e6;

export class MockStream {
  t = 0;
  theta = -360;
  cycle = 0;
  private op: OperatingPoint;
  private readonly B: number;
  private readonly S: number;
  private readonly L: number;
  // per-cycle state
  private pIVC = 0.95e5;
  private TIVC = 340;
  private mass = 0;
  private burnStart = -15;
  private burnDur = 45;
  private lw = 0;
  private knockAt = NaN;
  private knockAmp = 0;
  private sparkT = NaN;
  private energy = 0;
  private summary: { peakP: number; peakDeg: number; work: number; knockDeg: number; mapo: number; ca: number[] } = { peakP: 0, peakDeg: 0, work: 0, knockDeg: NaN, mapo: 0, ca: [NaN, NaN, NaN] };
  private readonly cycles: CycleSummary[] = [];
  private prevV = NaN;
  private prevP = NaN;

  constructor(
    private readonly spec: EngineSpec,
    op: OperatingPoint,
  ) {
    this.op = { ...op };
    this.B = spec.geometry.bore;
    this.S = spec.geometry.stroke;
    this.L = spec.geometry.conRodLength;
    this.newCycle();
  }

  setOperatingPoint(p: Partial<OperatingPoint>): void {
    this.op = { ...this.op, ...p };
  }

  reset(): void {
    this.t = 0;
    this.theta = -360;
    this.cycle = 0;
    this.cycles.length = 0;
    this.newCycle();
  }

  drainCycles(): CycleSummary[] {
    return this.cycles.splice(0);
  }

  private get rpm(): number {
    return this.op.rpm;
  }

  private volume(theta: number): { V: number; x: number; rod: number } {
    const a = this.S / 2;
    const th = (theta * Math.PI) / 180;
    const vd = (Math.PI / 4) * this.B * this.B * this.S;
    const vc = vd / (this.op.compressionRatio - 1);
    const s = a * Math.cos(th) + Math.sqrt(this.L * this.L - a * a * Math.sin(th) ** 2);
    const x = this.L + a - s;
    return { V: vc + (Math.PI / 4) * this.B * this.B * x, x, rod: Math.asin((a * Math.sin(th)) / this.L) };
  }

  private newCycle(): void {
    const ivc = this.spec.intakeValve.closeDeg;
    this.pIVC = 0.92e5 * (0.3 + 0.7 * this.op.throttle);
    this.TIVC = this.op.intakeMixtureTemperature + 25;
    this.mass = (this.pIVC * this.volume(ivc).V) / (R_AIR * this.TIVC);
    // Cycle-to-cycle variation of the burn.
    const jitter = (Math.random() - 0.5) * 6;
    this.burnStart = -this.op.sparkAdvanceDeg + 6 + jitter * 0.5;
    this.burnDur = 42 + jitter;
    this.lw = 0;
    this.knockAt = NaN;
    this.knockAmp = 0;
    this.sparkT = NaN;
    this.energy = 0;
    this.summary = { peakP: 0, peakDeg: 0, work: 0, knockDeg: NaN, mapo: 0, ca: [NaN, NaN, NaN] };
  }

  private xb(theta: number): number {
    if (theta < this.burnStart) return 0;
    const z = (theta - this.burnStart) / this.burnDur;
    return 1 - Math.exp(-5 * Math.pow(Math.max(0, z), 3));
  }

  private lift(v: EngineSpec['intakeValve'], theta: number): number {
    const [o, c] = valveOpenInterval(v);
    let w = wrapCycleDeg(theta);
    if (w < o) w += 720;
    if (w > c) return 0;
    return v.maxLift * Math.sin((Math.PI * (w - o)) / (c - o)) ** 2;
  }

  /** Advance by dθ degrees and return the snapshot. */
  step(dDeg: number): EngineSnapshot {
    const rpm = this.rpm;
    const dt = dDeg / (6 * rpm);
    const prevTheta = this.theta;
    this.t += dt;
    this.theta += dDeg;
    if (this.theta >= 360) {
      this.finishCycle();
      this.theta -= 720;
      this.cycle++;
      this.newCycle();
    }
    const th = this.theta;
    const { V, x, rod } = this.volume(th);
    const ivc = this.spec.intakeValve.closeDeg;
    const evo = this.spec.exhaustValve.openDeg;
    const Vivc = this.volume(ivc).V;
    let p: number;
    let Tu: number;
    let Tb = 0;
    const xb = this.xb(th);
    const closed = th >= ivc && th < evo;
    const fuel = (this.mass * this.op.equivalenceRatio) / (14.7 + this.op.equivalenceRatio);
    if (closed) {
      const pc = this.pIVC * Math.pow(Vivc / V, GAMMA_U);
      // Heat release raises pressure: Δp ≈ (γ−1) Q x_b / V (constant-volume approximation per step).
      const Q = 0.82 * fuel * LHV * xb;
      const pComb = ((GAMMA_B - 1) * Q) / V;
      p = pc + pComb * 0.9;
      Tu = this.TIVC * Math.pow(Vivc / V, GAMMA_U - 1) * Math.pow(p / pc, (GAMMA_U - 1) / GAMMA_U);
      if (xb > 0.001) Tb = Math.min(2900, 2300 + 500 * Math.pow(p / 60e5, 0.2) - 300 * (1 - xb) + (Math.random() - 0.5) * 2);
    } else if (th >= evo && th < 360 - 180 + 190) {
      // blowdown then exhaust stroke
      const pEvo = this.pIVC * Math.pow(Vivc / this.volume(evo).V, GAMMA_U) * 3.2;
      p = 1.05e5 + (pEvo - 1.05e5) * Math.exp(-(th - evo) / 18);
      Tu = 0;
      Tb = 1100 * Math.exp(-(th - evo) / 400);
    } else {
      p = th < -180 || th > 300 ? 0.95e5 + 0.05e5 * Math.sin((th * Math.PI) / 90) : this.pIVC * Math.pow(Vivc / V, GAMMA_U);
      Tu = this.TIVC;
      Tb = th > 300 ? 900 : 0;
    }
    // --- knock (Douaud–Eyzat ignition delay on the end gas) ---
    let osc = 0;
    let autoign = false;
    if (closed && xb > 0 && xb < 0.99 && Number.isNaN(this.knockAt)) {
      const on = this.op.fuel.kind === 'PRF' ? this.op.fuel.octaneNumber : 95;
      const tau = 17.68e-3 * Math.pow(on / 100, 3.402) * Math.pow(p / 101325, -1.7) * Math.exp(3800 / Tu);
      this.lw += dt / tau;
      if (this.lw >= 1) {
        this.knockAt = th;
        this.knockAmp = 0.25e5 * (1 - xb) * (p / 40e5) * 6;
        this.summary.knockDeg = th;
      }
    }
    if (!Number.isNaN(this.knockAt) && th >= this.knockAt && closed) {
      autoign = true;
      const tk = (th - this.knockAt) / (6 * rpm);
      const c = Math.sqrt(1.3 * R_AIR * 2400);
      const f = (1.841 * c) / (Math.PI * this.B);
      osc = this.knockAmp * Math.exp(-tk / 0.8e-3) * Math.sin(2 * Math.PI * f * tk);
      this.summary.mapo = Math.max(this.summary.mapo, Math.abs(osc));
    }
    const pTot = p + osc;
    if (pTot > this.summary.peakP) {
      this.summary.peakP = pTot;
      this.summary.peakDeg = th;
    }
    if (Number.isFinite(this.prevV)) this.summary.work += 0.5 * (pTot + this.prevP) * (V - this.prevV);
    this.prevV = V;
    this.prevP = pTot;
    for (const [k, lvl] of [0.1, 0.5, 0.9].entries()) if (Number.isNaN(this.summary.ca[k]) && xb >= lvl) this.summary.ca[k] = th;

    // --- ignition coil cartoon ---
    const sparkDeg = -this.op.sparkAdvanceDeg;
    const dwellDeg = this.op.dwellTime * 6 * rpm;
    let phase: SparkPhase = 'off';
    let iPrim = 0;
    let vSec = 0;
    let iSec = 0;
    const vBd = 2.5e3 + 0.35e3 * (p / 1e5) * (300 / Math.max(Tu, 250)) * 3;
    if (th >= sparkDeg - dwellDeg && th < sparkDeg) {
      phase = 'charging';
      const tc = (th - (sparkDeg - dwellDeg)) / (6 * rpm);
      iPrim = 7 * (1 - Math.exp(-tc / 1.2e-3));
      vSec = 0.3e3 * Math.exp(-tc / 0.2e-3);
    } else if (th >= sparkDeg && th < sparkDeg + 0.02 * rpm) {
      if (Number.isNaN(this.sparkT)) this.sparkT = this.t;
      const ts = this.t - this.sparkT;
      const tRise = 25e-6;
      if (ts < tRise) {
        phase = 'charging';
        vSec = Math.min(vBd, 14e3 * Math.sin((Math.PI / 2) * (ts / tRise)));
        iPrim = 7 * (1 - ts / tRise);
      } else if (ts < tRise + 2e-6) {
        phase = 'breakdown';
        vSec = vBd;
        iSec = 150;
      } else if (ts < tRise + 60e-6) {
        phase = 'arc';
        vSec = 120;
        iSec = 0.2 * Math.exp(-(ts - tRise) / 40e-6) + 0.09;
      } else if (ts < tRise + 1.6e-3) {
        phase = 'glow';
        const g = (ts - tRise - 60e-6) / 1.54e-3;
        vSec = 600 + 200 * g;
        iSec = 0.09 * (1 - g);
      } else {
        phase = 'done';
      }
      this.energy = Math.min(40e-3, this.energy + vSec * iSec * dt);
    } else if (th >= sparkDeg) {
      phase = 'done';
    }

    const Vm = V;
    const Tm = closed ? (pTot * Vm) / (this.mass * R_AIR) : p > 0 ? Math.max(Tu, Tb * 0.8) : 0;
    const liftIn = this.lift(this.spec.intakeValve, th);
    const liftEx = this.lift(this.spec.exhaustValve, th);
    const hrr = closed ? 0.82 * fuel * LHV * ((this.xb(th + 0.05) - this.xb(th - 0.05)) / 0.1) * 6 * rpm : 0;
    void prevTheta;
    return {
      t: this.t,
      cycle: this.cycle,
      thetaDeg: th,
      rpm,
      pistonDisplacement: x,
      clearanceHeight: V / ((Math.PI / 4) * this.B * this.B),
      rodAngle: rod,
      intakeLift: liftIn,
      exhaustLift: liftEx,
      phase: closed ? (xb > 0 && xb < 0.99 ? 'combustion' : th < 0 ? 'compression' : 'expansion') : 'gas-exchange',
      volume: V,
      pressure: pTot,
      temperatureMean: Tm,
      temperatureUnburned: Tu,
      temperatureBurned: Tb,
      massFractionBurned: closed || th >= evo ? xb : 0,
      mass: this.mass,
      heatReleaseRate: hrr,
      heatLossRate: 0,
      flame: {
        stage: !closed ? 'none' : xb <= 0 ? (phase === 'arc' || phase === 'glow' ? 'kernel' : 'none') : xb < 0.02 ? 'kernel' : xb < 0.9 ? 'turbulent' : xb < 0.999 ? 'burnout' : 'done',
        radius: 0.04 * Math.cbrt(xb),
        center: this.spec.sparkPlug.gapCenter,
        area: 0,
        laminarSpeed: 0.4,
        turbulentSpeed: 2,
        turbulenceIntensity: 2,
      },
      spark: { phase, primaryCurrent: iPrim, secondaryVoltage: vSec, secondaryCurrent: iSec, energyDelivered: this.energy, breakdownVoltage: vBd },
      intakeMassFlow: liftIn > 0 ? 40 * (liftIn / this.spec.intakeValve.maxLift) * Math.sin((Math.PI * (th + 360)) / 200) * 1e-3 : 0,
      exhaustMassFlow: liftEx > 0 ? 60 * (liftEx / this.spec.exhaustValve.maxLift) * Math.exp(-Math.max(0, th - evo) / 60) * 1e-3 : 0,
      intakeManifoldPressure: 0.95e5,
      exhaustManifoldPressure: 1.05e5,
      knock: { integral: this.lw, autoignited: autoign, oscillation: osc },
      burnedComposition: { CO2: 0.12, H2O: 0.13, CO: 0.005, O2: 0.01, H2: 0.001, OH: 0.001, H: 0, O: 0, NO: 0.002, N2: 0.72 },
      gasTorque: 0,
      netTorque: 0,
    };
  }

  private finishCycle(): void {
    const vd = (Math.PI / 4) * this.B * this.B * this.S;
    const fuel = (this.mass * this.op.equivalenceRatio) / (14.7 + this.op.equivalenceRatio);
    const W = this.summary.work;
    this.cycles.push({
      cycle: this.cycle,
      imepGross: W / vd + 0.1e5,
      imepNet: W / vd,
      pmep: -0.1e5,
      peakPressure: this.summary.peakP,
      peakPressureDeg: this.summary.peakDeg,
      maxPressureRiseRate: 2.5e5,
      ca10: this.summary.ca[0],
      ca50: this.summary.ca[1],
      ca90: this.summary.ca[2],
      indicatedEfficiency: W / (fuel * LHV),
      isfc: fuel / Math.max(W, 1),
      trappedMass: this.mass,
      residualFraction: 0.06,
      volumetricEfficiency: 0.82,
      fuelMass: fuel,
      noPpm: 1800 + Math.random() * 200,
      coFraction: 0.004,
      knockOnsetDeg: this.summary.knockDeg,
      knockEndGasFraction: Number.isFinite(this.summary.knockDeg) ? 0.12 : NaN,
      mapo: this.summary.mapo,
      misfire: false,
      heatLoss: 90,
      indicatedWorkGross: W + 0.1e5 * vd,
    });
  }
}
