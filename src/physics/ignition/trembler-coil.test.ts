import { describe, expect, it } from 'vitest';
import type { MagnetoSpec, TremblerCoilSpec, TremblerVibratorSpec } from '../core/engine-spec';
import fx from '../../../test/fixtures/ignition_trembler.json';
import { SparkGap } from './discharge';
import { PrimarySupply } from './source';
import { TremblerCoil, type TremblerCoilOptions } from './trembler-coil';
import { Vibrator } from './vibrator';

const REPORT = !!process.env.IGNITION_REPORT;
const P = fx.params;

/** The oracle's coil (Ford/K-W 1916 data + fitted vibrator). */
function oracleCoil(vib: Partial<TremblerVibratorSpec> = {}): TremblerCoilSpec {
  const c = P.coil;
  const v = P.vibrator;
  return {
    primaryInductance: c.L1,
    primaryResistance: c.R1,
    secondaryInductance: c.L2,
    secondaryResistance: c.R2,
    secondaryCapacitance: c.C2,
    couplingCoefficient: c.k,
    condenserCapacitance: c.C1,
    vibrator: { pullCurrent: v.Ip, naturalFrequency: v.fn, dampingRatio: v.zeta, breakTravel: v.xb, maxTravel: v.xmax, airGap: v.g0, ...vib },
  };
}

const ORACLE_MAGNETO: MagnetoSpec = {
  cyclesPerRevolution: P.magneto.N,
  emfConstant: P.magneto.k,
  phaseDeg: P.magneto.phi,
  internalResistance: P.magneto.Rs,
  internalInductance: P.magneto.Ls,
};

/** Spark gap with the oracle's glow law at a fixed breakdown voltage (glow only, as the oracle). */
class LoggingGap extends SparkGap {
  readonly bdTimes: number[] = [];
  override startConduction(vAbs: number, iAbs: number, c2: number, t: number): number {
    const n = this.breakdownCount;
    const r = super.startConduction(vAbs, iAbs, c2, t);
    if (this.breakdownCount > n) this.bdTimes.push(t);
    return r;
  }
}
function oracleGap(g: { Vbd: number; p: number; l: number }, suppress = false): LoggingGap {
  const gap = new LoggingGap(g.l, { glowCathodeFall: P.glow.Vsheath, glowAnodeFall: 0, extinctionCurrent: P.glow.Iext, arcGlowCurrent: 1e3 });
  gap.pressure = g.p;
  gap.breakdownVoltage = g.Vbd;
  gap.suppressBreakdown = suppress;
  return gap;
}

interface Run {
  coil: TremblerCoil;
  gap: LoggingGap;
  trips: number[];
  tripI: number[];
  closes: number[];
  /** Ledger at the end of the caller step containing the first re-close. */
  atFirstClose: { E_sup: number; E_R: number; E_R2: number; E_pts: number; E_gap: number } | null;
}

/** Drive a coil with the timer closed on [0, tBreak) up to tEnd with caller step dt, logging events. */
function drive(coil: TremblerCoil, gap: LoggingGap, tEnd: number, dt: number, tBreak = Infinity): Run {
  const run: Run = { coil, gap, trips: [], tripI: [], closes: [], atFirstClose: null };
  let nt = 0;
  let lastClose = NaN;
  let t = 0;
  while (t < tEnd - 1e-12) {
    let h = Math.min(dt, tEnd - t);
    if (t < tBreak && t + h > tBreak) h = tBreak - t;
    coil.step(h, t < tBreak - 1e-15, gap);
    t += h;
    while (coil.tripCount > nt) {
      nt++;
      run.trips.push(coil.lastTripTime);
      run.tripI.push(coil.lastTripCurrent);
    }
    if (!Number.isNaN(coil.lastCloseTime) && coil.lastCloseTime !== lastClose) {
      lastClose = coil.lastCloseTime;
      run.closes.push(lastClose);
      if (!run.atFirstClose)
        run.atFirstClose = {
          E_sup: coil.energySupplied,
          E_R: coil.energyPrimaryResistance + coil.energyExternalResistance,
          E_R2: coil.energySecondaryResistance,
          E_pts: coil.energyPoints,
          E_gap: coil.energyGap,
        };
    }
  }
  return run;
}

function bench(V: number, tEnd: number, dt = 20e-6, opts?: Partial<TremblerCoilOptions>, suppress = false): Run {
  const coil = new TremblerCoil(oracleCoil(), PrimarySupply.dc(V, P.batteryResistance), 0, opts);
  return drive(coil, oracleGap(P.benchGap, suppress), tEnd, dt);
}

describe('Vibrator — analytic armature motion', () => {
  it('reaches the break travel with the energy-integral velocity and quadrature time (ζ = 0, constant current)', () => {
    const spec: TremblerVibratorSpec = { ...oracleCoil().vibrator, dampingRatio: 0 };
    const vib = new Vibrator(spec);
    const I = 5;
    const w2 = vib.omega ** 2;
    const xp = vib.preload;
    const g0 = spec.airGap;
    const ip = spec.pullCurrent;
    // ½v² = ∫₀ˣ ẍ ds = ω²[x_p (I/I_p)² g₀² (1/(g₀−x) − 1/g₀) − x_p x − x²/2]
    const vOf = (x: number): number => Math.sqrt(2 * w2 * (xp * (I / ip) ** 2 * g0 * g0 * (1 / (g0 - x) - 1 / g0) - xp * x - 0.5 * x * x));
    const xb = spec.breakTravel;
    // t(x_b) = ∫ dx/v with x = x_b u² (removes the 1/√x singularity at the start)
    let tq = 0;
    const n = 20000;
    for (let k = 0; k < n; k++) {
      const u = (k + 0.5) / n;
      tq += (2 * xb * u) / vOf(xb * u * u) / n;
    }
    const h = 20e-6;
    let t = 0;
    let tCross = NaN;
    let vCross = NaN;
    for (let k = 0; k < 5000 && Number.isNaN(tCross); k++) {
      vib.trial(h, I, I);
      if (vib.xNew >= xb) {
        const f = (xb - vib.x) / (vib.xNew - vib.x);
        tCross = t + f * h;
        vCross = vib.v + f * (vib.vNew - vib.v);
      }
      vib.commit();
      t += h;
    }
    expect(Math.abs(tCross / tq - 1)).toBeLessThan(2e-3);
    expect(Math.abs(vCross / vOf(xb) - 1)).toBeLessThan(5e-3);
    if (REPORT) console.log(`vibrator: t(x_b) ${tCross * 1e6} µs vs ${tq * 1e6} µs; v ${vCross} vs ${vOf(xb)} m/s`);
  });

  it('stays on its stop below the pull-in current and snaps through above it (no static equilibrium)', () => {
    const vib = new Vibrator(oracleCoil().vibrator);
    const ip = vib.spec.pullCurrent;
    for (let k = 0; k < 1000; k++) {
      vib.trial(20e-6, 0.999 * ip, 0.999 * ip);
      vib.commit();
    }
    expect(vib.x).toBe(0);
    for (let k = 0; k < 2000 && vib.x < vib.maxTravel; k++) {
      vib.trial(20e-6, 1.02 * ip, 1.02 * ip);
      vib.commit();
    }
    expect(vib.x).toBe(vib.maxTravel);
  });
});

describe('TremblerCoil — analytic limits', () => {
  it('DC ramp up to the first trip follows the RL law I = (V/R)(1 − e^{−Rt/L}) (Heywood eq. 9.53)', () => {
    // secondary uncoupled (k → 0): with k = 0.9 the lossless coil's make ring (L₂(1−k²)C₂ ≈ 12 kHz,
    // ±k√(L₂/L₁)·V ≈ ±0.4 kV on the plug) reflects ≈ ±0.1 A onto I₁ for milliseconds (the coupled ramp
    // is checked against the Radau oracle below)
    for (const V of [6, 12]) {
      const coil = new TremblerCoil({ ...oracleCoil(), couplingCoefficient: 1e-9 }, PrimarySupply.dc(V, 0.05), 0.1);
      const gap = oracleGap(P.benchGap);
      const R = coil.loopResistance;
      const L = coil.loopInductance;
      expect(R).toBeCloseTo(0.295 + 0.05 + 0.1, 12);
      let t = 0;
      let maxErr = 0;
      while (coil.tripCount === 0) {
        coil.step(20e-6, true, gap);
        t += 20e-6;
        if (coil.tripCount === 0) maxErr = Math.max(maxErr, Math.abs(coil.I1 / ((V / R) * (1 - Math.exp((-R * t) / L))) - 1));
      }
      expect(maxErr).toBeLessThan(1e-4);
      const tt = coil.firstTripTime;
      expect(coil.firstTripCurrent / ((V / R) * (1 - Math.exp((-R * tt) / L)))).toBeCloseTo(1, 4);
    }
  });

  it('sinusoid-driven series R–(L₁ + L_s) circuit matches the analytic response incl. the transient', () => {
    // magneto at 1000 rpm, make at a rising EMF zero; vibrator never trips, gap open
    const coil = new TremblerCoil(oracleCoil({ pullCurrent: 1e9 }), new PrimarySupply(null, ORACLE_MAGNETO, 'magneto'), 0.1);
    const gap = oracleGap(P.benchGap, true);
    const rpm = 1000;
    const w = (rpm * 2 * Math.PI) / 60;
    const we = ORACLE_MAGNETO.cyclesPerRevolution * w;
    const E = ORACLE_MAGNETO.emfConstant * w;
    const R = coil.loopResistance;
    const L = coil.loopInductance;
    expect(L).toBeCloseTo(3.3e-3 + 3.05e-3, 12);
    const Z = Math.hypot(R, we * L);
    const th = Math.atan2(we * L, R);
    const dt = 37e-6;
    let t = 0;
    let maxErr = 0;
    for (let k = 0; k < 400; k++) {
      coil.supply.setStep(t, dt, ORACLE_MAGNETO.phaseDeg + rpm * 6 * t, ORACLE_MAGNETO.phaseDeg + rpm * 6 * (t + dt), w, w);
      coil.invalidateEmf();
      coil.step(dt, true, gap);
      t += dt;
      const ia = (E / Z) * (Math.sin(we * t - th) + Math.sin(th) * Math.exp((-R * t) / L));
      maxErr = Math.max(maxErr, Math.abs(coil.I1 - ia));
    }
    expect(coil.tripCount).toBe(0);
    expect(maxErr / (E / Z)).toBeLessThan(2e-3);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });
});

describe('TremblerCoil — Radau oracle (tools/reference/ignition_trembler_oracle.py)', () => {
  const cases = [fx.bench6, fx.bench9, fx.bench12];
  for (const c of cases) {
    it(`DC bench ${c.V} V: first trip, firing current, spark train and re-close`, () => {
      const r = bench(c.V, c.t_end);
      expect(Math.abs(r.trips[0] - c.trips[0])).toBeLessThan(1e-6);
      // instantaneous firing current: the 20 µs ramp sub-step aliases the secondary make ring (±1.5 %)
      expect(Math.abs(r.tripI[0] / c.tripI[0] - 1)).toBeLessThan(0.02);
      // the spark train of the first trip: every breakdown before the re-close
      const nOracle = c.bds.filter((t) => t < c.closes[0]).length;
      const tsBds = r.gap.bdTimes.filter((t) => t < r.closes[0]);
      expect(tsBds.length).toBe(nOracle);
      for (let k = 0; k < nOracle; k++) expect(Math.abs(tsBds[k] - c.bds[k])).toBeLessThan(k === 0 ? 1e-6 : 8e-6);
      // re-close (depends on the force history through the whole train): ≤ 20 µs of ≈ 1.8–2.7 ms
      expect(Math.abs(r.closes[0] - c.closes[0])).toBeLessThan(20e-6);
      const a = r.atFirstClose!;
      const o = c.atFirstClose;
      expect(a.E_gap / o.E_gap).toBeCloseTo(1, 2);
      expect(a.E_R2 / o.E_R2).toBeCloseTo(1, 2);
      expect(Math.abs(a.E_sup / o.E_sup - 1)).toBeLessThan(5e-3);
      expect(Math.abs(a.E_R / o.E_R - 1)).toBeLessThan(5e-3);
      // the next buzz cycle (sensitive to the re-close state): within 2 % of its period
      expect(Math.abs(r.trips[1] - c.trips[1])).toBeLessThan(0.02 * (c.trips[1] - c.trips[0]));
      expect(Math.abs(r.coil.energyResidual())).toBeLessThan(1e-12);
      if (REPORT)
        console.log(
          `bench ${c.V} V: trip ${(r.trips[0] * 1e3).toFixed(4)} vs ${(c.trips[0] * 1e3).toFixed(4)} ms, close ${(r.closes[0] * 1e3).toFixed(4)} vs ${(c.closes[0] * 1e3).toFixed(4)} ms, ` +
            `bds ${tsBds.length}, E_gap ${(a.E_gap * 1e3).toFixed(3)} vs ${(o.E_gap * 1e3).toFixed(3)} mJ, trip 2 ${(r.trips[1] * 1e3).toFixed(3)} vs ${(c.trips[1] * 1e3).toFixed(3)} ms`,
        );
    });
  }

  for (const c of [fx.open6, fx.open12]) {
    it(`open circuit ${c.V} V: ring-up peaks after the first trip (${(c.peakV1).toFixed(0)} V on the points, ${(c.peakV2 / 1e3).toFixed(1)} kV)`, () => {
      const coil = new TremblerCoil(oracleCoil(), PrimarySupply.dc(c.V, P.batteryResistance), 0);
      const gap = oracleGap(P.benchGap, true);
      let pV1 = 0;
      let pV2 = 0;
      let tV2 = 0;
      for (let t = 0; t < 4.5e-3; t += 5e-6) {
        coil.step(5e-6, true, gap);
        if (coil.tripCount === 1) {
          pV1 = Math.max(pV1, Math.abs(coil.V1));
          if (Math.abs(coil.V2) > pV2) {
            pV2 = Math.abs(coil.V2);
            tV2 = coil.time - coil.firstTripTime;
          }
        }
      }
      expect(Math.abs(coil.firstTripTime - c.trip)).toBeLessThan(1e-6);
      // ∝ the firing current (±1.5 %, see the bench cases)
      expect(Math.abs(pV1 / c.peakV1 - 1)).toBeLessThan(0.02);
      expect(Math.abs(pV2 / c.peakV2 - 1)).toBeLessThan(0.02);
      expect(Math.abs(tV2 - c.tpeakV2)).toBeLessThan(6e-6);
      expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
    });
  }

  /** Magneto through the timer at constant speed, make at `make` deg, break after 87°. */
  function magnetoRun(c: { rpm: number; make: number; t_break: number; t_end: number }, dt = 20e-6): Run {
    const sup = new PrimarySupply(null, ORACLE_MAGNETO, 'magneto');
    const w = (c.rpm * 2 * Math.PI) / 60;
    sup.setStep(0, 1, c.make, c.make + c.rpm * 6, w, w);
    const coil = new TremblerCoil(oracleCoil(), sup, P.timerResistance);
    return drive(coil, oracleGap(P.engineGap), c.t_end, dt, c.t_break);
  }

  it('magneto 1000 rpm through the timer: trips, firing currents, first spark angle, timer-break energy', () => {
    const c = fx.magneto1000;
    const r = magnetoRun(c);
    expect(r.trips.length).toBe(c.trips.length);
    expect(Math.abs(r.trips[0] - c.trips[0])).toBeLessThan(1e-6);
    expect(Math.abs(r.tripI[0] / c.tripI[0] - 1)).toBeLessThan(0.02);
    for (let k = 1; k < c.trips.length; k++) expect(Math.abs(r.trips[k] - c.trips[k])).toBeLessThan(30e-6);
    const deg = c.make + c.rpm * 6 * r.gap.bdTimes[0];
    expect(Math.abs(deg - c.firstSparkDeg!)).toBeLessThan(0.03); // 5 µs at 1000 rpm
    expect(Math.abs(r.gap.bdTimes.length - c.bds.length)).toBeLessThanOrEqual(2);
    // the current broken by the timer depends on the buzz phase at the break (accumulated over 14.5 ms):
    // absolute tolerance, ≪ the ≈ 50 mJ of one spark (the break itself: analytic test below)
    expect(Math.abs(r.coil.energyTimerBreak - c.E_timer)).toBeLessThan(0.5e-3);
    expect(r.coil.energyGap / c.E_gap).toBeCloseTo(1, 1);
    expect(r.coil.energySupplied / c.E_sup).toBeCloseTo(1, 2);
    expect(Math.abs(r.coil.energyResidual())).toBeLessThan(1e-12);
    if (REPORT)
      console.log(`magneto 1000: trips ${r.trips.map((t) => (t * 1e3).toFixed(4))} vs ${c.trips.map((t) => (t * 1e3).toFixed(4))} ms, bds ${r.gap.bdTimes.length} vs ${c.bds.length}, E_timer ${r.coil.energyTimerBreak} vs ${c.E_timer}`);
  });

  it('slow magneto pulses (150 rpm): the firing current approaches the pull-in current', () => {
    const c = fx.slow150;
    const r = magnetoRun(c);
    expect(Math.abs(r.trips[0] - c.trips[0])).toBeLessThan(2e-6);
    expect(Math.abs(r.tripI[0] / c.tripI[0] - 1)).toBeLessThan(0.02);
  });

  it('magneto 600 rpm spark ladder: first-spark angle vs timer make (33 lever positions)', () => {
    const c = fx.ladder600;
    let maxErr = 0;
    for (let k = 0; k < c.makes.length; k++) {
      const make = c.makes[k];
      const r = magnetoRun({ rpm: c.rpm, make, t_break: 87 / (6 * c.rpm), t_end: 87 / (6 * c.rpm) });
      const deg = make + c.rpm * 6 * r.gap.bdTimes[0];
      maxErr = Math.max(maxErr, Math.abs(deg - c.firstSparkDeg[k]!));
    }
    expect(maxErr).toBeLessThan(0.02);
  });
});

describe('TremblerCoil — ledger, disconnected primary, timer break, caller step', () => {
  it('closes the energy ledger to ≤ 1e-12 J over 50 ms of buzzing (≥ 9 breaks, ≥ 80 breakdowns)', () => {
    const r = bench(6, 50e-3);
    expect(r.trips.length).toBeGreaterThanOrEqual(9);
    expect(r.gap.breakdownCount).toBeGreaterThanOrEqual(80);
    expect(Math.abs(r.coil.energyResidual())).toBeLessThan(1e-12);
    if (REPORT) console.log(`50 ms: ${r.trips.length} trips, ${r.gap.breakdownCount} breakdowns, residual ${r.coil.energyResidual()} J of ${r.coil.energySupplied} J`);
  });

  it('timer break: I₁ → 0 with the secondary flux linkage conserved, ½(L_p − M²/L₂)I₁² to the contact', () => {
    const coil = new TremblerCoil(oracleCoil(), PrimarySupply.dc(6, 0.05), 0.1);
    const gap = oracleGap(P.engineGap);
    for (let t = 0; t < 2e-3; t += 20e-6) coil.step(20e-6, true, gap);
    expect(coil.pointsOpen).toBe(false);
    const i1 = coil.I1;
    const i2 = coil.I2;
    const M = coil.mutualInductance;
    const L2 = coil.coil.secondaryInductance;
    const e0 = coil.energyTimerBreak;
    coil.step(0, false, gap);
    expect(coil.I1).toBe(0);
    expect(L2 * coil.I2).toBeCloseTo(M * i1 + L2 * i2, 12);
    expect((coil.energyTimerBreak - e0) / (0.5 * (coil.loopInductance - (M * M) / L2) * i1 * i1)).toBeCloseTo(1, 10);
    // the flux transferred to the secondary rings it up and fires the plug ("timer spark")
    for (let t = 0; t < 0.5e-3; t += 20e-6) coil.step(20e-6, false, gap);
    expect(gap.breakdownCount).toBeGreaterThan(0);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('holds I₁ ≡ 0 and the condenser charge while the timer is open, then settles (points close, coil skipped)', () => {
    const coil = new TremblerCoil(oracleCoil(), PrimarySupply.dc(6, 0.05), 0.1);
    const gap = oracleGap(P.benchGap);
    while (coil.tripCount === 0) coil.step(20e-6, true, gap);
    for (let k = 0; k < 10; k++) coil.step(20e-6, true, gap); // ring under way, points open
    expect(coil.pointsOpen).toBe(true);
    coil.step(0, false, gap);
    let held = true;
    let v1 = NaN;
    let ePts = NaN;
    while (coil.pointsOpen) {
      coil.step(20e-6, false, gap);
      held &&= coil.I1 === 0;
      if (Number.isNaN(v1) && coil.pointsOpen) v1 = coil.V1;
      if (coil.pointsOpen) held &&= coil.V1 === v1;
      ePts = coil.energyPoints;
    }
    expect(held).toBe(true);
    expect(coil.V1).toBe(0);
    expect(ePts).toBeGreaterThan(0);
    // the free L₂R₂C₂ ring decays with 2L₂/R₂ ≈ 13 ms and is snapped to rest below ½C₂(25 V)²
    for (let k = 0; k < 10000; k++) coil.step(20e-6, false, gap);
    expect(coil.I1).toBe(0);
    expect(coil.vibrator.x).toBe(0);
    expect(coil.I2).toBe(0);
    expect(coil.V2).toBe(0);
    expect(Math.abs(coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('is insensitive to the caller step (5, 50, 200 µs): trip, train, re-close', () => {
    const runs = [5e-6, 50e-6, 200e-6].map((dt) => bench(6, 6e-3, dt));
    const a = runs[0];
    for (const b of runs.slice(1)) {
      expect(Math.abs(b.trips[0] - a.trips[0])).toBeLessThan(0.5e-6);
      expect(b.gap.bdTimes.filter((t) => t < b.closes[0]).length).toBe(a.gap.bdTimes.filter((t) => t < a.closes[0]).length);
      expect(Math.abs(b.closes[0] - a.closes[0])).toBeLessThan(5e-6);
      expect(b.atFirstClose!.E_gap / a.atFirstClose!.E_gap).toBeCloseTo(1, 2);
    }
  });
});
