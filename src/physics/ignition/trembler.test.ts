/**
 * Trembler-magneto ignition (Ford Model T, MODEL_T_IGNITION) at system level: the fitted vibrator and
 * magneto against the sourced targets, spark trains, the timer/lever command, battery vs magneto, the
 * four-cylinder set-up and the CPU cost. Fit residuals are printed with IGNITION_REPORT=1.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSpec, TremblerCoilSpec, TremblerMagnetoIgnitionSpec } from '../core/engine-spec';
import { cylinderAngleDeg } from '../core/engine-spec';
import gasFixture from '../../../test/fixtures/ignition_kernel_gas.json';
import oracleFx from '../../../test/fixtures/ignition_trembler.json';
import { MODEL_T, MODEL_T_IGNITION } from '../engines/model-t';
import { DEFAULT_SPARK_GAP_OPTIONS, REIGNITION_WINDOW, SparkGap } from './discharge';
import { createIgnitionSystems, IgnitionSystem, magnetoEmf, PrimarySupply, TremblerCoil, tremblerCoilOf, tremblerTimerCommand, type IgnitionGasState } from './index';

const REPORT = !!process.env.IGNITION_REPORT;
const log = (s: string): void => {
  if (REPORT) console.log(s);
};
const IG = MODEL_T_IGNITION;

/** Model-T-like spark state of the Cantera fixture: CR 4.5, φ 0.6 (5.3 bar, 526 K). */
const C45 = gasFixture.cases.find((k) => k.cr === 4.5)!;
function gasT(SL = 0.3, Ma = 2, uPrime = 1, lI = 3e-3): IgnitionGasState {
  const lF = C45.nuU / SL;
  return {
    p: C45.p,
    Tu: C45.Tu,
    rhoU: C45.rhoU,
    X: Float64Array.from(C45.X),
    SL,
    marksteinLength: Ma * lF,
    flameThickness: lF,
    uPrime,
    integralScale: lI,
    expansionRatio: C45.sigma,
    kinematicViscosity: C45.nuU,
  };
}

/** Coil-tester bench: battery + 0.05 Ω straight onto the coil (no timer), 10 kV air gap. */
function benchRun(V: number, tEnd: number, coilSpec: TremblerCoilSpec = IG.coil) {
  const coil = new TremblerCoil(coilSpec, PrimarySupply.dc(V, IG.battery.internalResistance), 0);
  const gap = new SparkGap(3e-3);
  gap.pressure = 1e5;
  gap.breakdownVoltage = 10e3;
  const trips: number[] = [];
  const tripI: number[] = [];
  const tripE: number[] = [];
  const closes: number[] = [];
  let lastClose = NaN;
  for (let t = 0; t < tEnd - 1e-12; t += 20e-6) {
    coil.step(20e-6, true, gap);
    if (coil.tripCount > trips.length) {
      trips.push(coil.lastTripTime);
      tripI.push(coil.lastTripCurrent);
      tripE.push(coil.energyAtTrip);
    }
    if (!Number.isNaN(coil.lastCloseTime) && coil.lastCloseTime !== lastClose) closes.push((lastClose = coil.lastCloseTime));
  }
  return { coil, trips, tripI, tripE, closes };
}

/** Crank-driven first spark: constant speed, lever at `makeDeg` (timer make, local angle). */
function firstSpark(spec: TremblerMagnetoIgnitionSpec, rpm: number, makeDeg: number, source: 'magneto' | 'battery' = 'magneto', gas = gasT()) {
  const ign = new IgnitionSystem(spec, MODEL_T.sparkPlug, { ignitionSource: source });
  const cmd = tremblerTimerCommand(spec, -makeDeg);
  const w = (rpm * 2 * Math.PI) / 60;
  const dt = 20e-6;
  let th = makeDeg - 3;
  ign.step(dt, th, cmd, gas, w);
  for (let i = 0; i < 200000 && th < makeDeg + spec.timer.contactArcDeg; i++) {
    th += rpm * 6 * dt;
    const s = ign.step(dt, th, cmd, gas, w);
    if (s.breakdownCount > 0) return { deg: s.firstSparkDeg, delay: s.firstSparkDelay, trip: s.firstTripDelay, I: s.firstTripCurrent };
  }
  return { deg: NaN, delay: NaN, trip: NaN, I: NaN };
}

describe('MODEL_T_IGNITION — vibrator fit to the coil-tester data (ECCT; Cool386)', () => {
  it('uses the parameters of the Radau oracle / fit (tools/reference/ignition_trembler_oracle.py); oracle residuals ≤ 5 µs', () => {
    const p = oracleFx.params;
    const c = IG.coil;
    expect([c.primaryInductance, c.primaryResistance, c.secondaryInductance, c.secondaryResistance, c.secondaryCapacitance, c.couplingCoefficient, c.condenserCapacitance]).toEqual([
      p.coil.L1, p.coil.R1, p.coil.L2, p.coil.R2, p.coil.C2, p.coil.k, p.coil.C1,
    ]);
    const v = c.vibrator;
    expect([v.pullCurrent, v.naturalFrequency, v.dampingRatio, v.airGap]).toEqual([p.vibrator.Ip, p.vibrator.fn, p.vibrator.zeta, p.vibrator.g0]);
    expect(v.breakTravel).toBeCloseTo(p.vibrator.xb, 15);
    expect(v.maxTravel).toBeCloseTo(p.vibrator.xmax, 15);
    const m = IG.magneto;
    expect([m.cyclesPerRevolution, m.emfConstant, m.phaseDeg, m.internalResistance, m.internalInductance]).toEqual([p.magneto.N, p.magneto.k, p.magneto.phi, p.magneto.Rs, p.magneto.Ls]);
    expect(IG.timer.contactResistance).toBe(p.timerResistance);
    expect(IG.battery.internalResistance).toBe(p.batteryResistance);
    // the oracle's re-ignition law is the gap model's
    expect([DEFAULT_SPARK_GAP_OPTIONS.reignitionTime, REIGNITION_WINDOW, DEFAULT_SPARK_GAP_OPTIONS.transitionHysteresis]).toEqual([p.reignition.tau, p.reignition.window, p.reignition.hyst]);
    const f = oracleFx.fit;
    expect(Math.abs(f.achieved.t6 - f.targets.t6)).toBeLessThan(5e-6);
    expect(Math.abs(f.achieved.t9 - f.targets.t9)).toBeLessThan(5e-6);
    expect(Math.abs(f.achieved.t12 - f.targets.t12)).toBeLessThan(5e-6);
    // The vibrator in use was fitted before the re-ignition model (discharge.ts): each current zero then
    // re-broke the 10 kV bench gap and dumped ½C₂V_bd² = 2 mJ, which damped the ring; re-igniting at V_r(t)
    // the spark dissipates less and the points re-close ≈ 110 µs later. The refit vibrator (same targets, with
    // re-ignition: oracle VIB_REFIT, proposed for MODEL_T_IGNITION) meets all four again.
    expect(Math.abs(f.achieved.reclose6 - f.targets.reclose6)).toBeLessThan(0.12e-3);
    const rf = oracleFx.refit.achieved;
    expect(Math.abs(rf.t6 - f.targets.t6)).toBeLessThan(5e-6);
    expect(Math.abs(rf.t9 - f.targets.t9)).toBeLessThan(5e-6);
    expect(Math.abs(rf.t12 - f.targets.t12)).toBeLessThan(5e-6);
    expect(Math.abs(rf.reclose6 - f.targets.reclose6)).toBeLessThan(15e-6);
    const us = (a: { t6: number; t9: number; t12: number; reclose6: number }) =>
      `${((a.t6 - f.targets.t6) * 1e6).toFixed(1)} / ${((a.t9 - f.targets.t9) * 1e6).toFixed(1)} / ${((a.t12 - f.targets.t12) * 1e6).toFixed(1)} µs (6/9/12 V), re-close ${((a.reclose6 - f.targets.reclose6) * 1e6).toFixed(1)} µs`;
    log(`oracle fit residuals (vibrator in use): ${us(f.achieved)}; firing currents ${f.achieved.I6.toFixed(2)} / ${f.achieved.I9.toFixed(2)} / ${f.achieved.I12.toFixed(2)} A; buzz ${f.achieved.buzz6.toFixed(0)} Hz`);
    log(`oracle refit residuals (${JSON.stringify(oracleFx.refit.vibrator)}): ${us(rf)}; firing currents ${rf.I6.toFixed(2)} / ${rf.I9.toFixed(2)} / ${rf.I12.toFixed(2)} A; buzz ${rf.buzz6.toFixed(0)} Hz`);
  });

  const r6 = benchRun(6, 20e-3);
  const r9 = benchRun(9, 8e-3);
  const r12 = benchRun(12, 8e-3);
  /** The coil with the refit vibrator (oracle VIB_REFIT; contract request for MODEL_T_IGNITION). */
  const rv = oracleFx.refit.vibrator;
  const REFIT_COIL = { ...IG.coil, vibrator: { ...IG.coil.vibrator, pullCurrent: rv.Ip, naturalFrequency: rv.fn, dampingRatio: rv.zeta, airGap: rv.g0 } };
  const q6 = benchRun(6, 20e-3, REFIT_COIL);
  const q9 = benchRun(9, 8e-3, REFIT_COIL);
  const q12 = benchRun(12, 8e-3, REFIT_COIL);

  it('first fire after the timer make: 3.5 ms / 6 V, 2.5 ms / 9 V, 2.0 ms / 12 V at ≈ 5.0–5.4, 6.2, 6–7 A', () => {
    const rows = [
      [6, r6, 3.5e-3, [5.0, 5.45]],
      [9, r9, 2.5e-3, [5.8, 6.5]],
      [12, r12, 2.0e-3, [6.0, 7.0]],
    ] as const;
    for (const [V, r, tTarget, [iLo, iHi]] of rows) {
      expect(Math.abs(r.trips[0] - tTarget)).toBeLessThan(0.02e-3);
      expect(r.tripI[0]).toBeGreaterThan(iLo);
      expect(r.tripI[0]).toBeLessThan(iHi);
      log(`fit ${V} V: fire ${(r.trips[0] * 1e3).toFixed(3)} ms (target ${tTarget * 1e3}, residual ${((r.trips[0] - tTarget) * 1e6).toFixed(1)} µs) at ${r.tripI[0].toFixed(2)} A`);
    }
  });

  it('points re-close ≈ 1.8 ms after the 6 V fire; buzz ≈ 190–200 Hz; points stay open longer at higher current', () => {
    const check = (a: typeof r6, b: typeof r9, c: typeof r12, reclose: number, buzz: [number, number], label: string) => {
      const off6 = a.closes[0] - a.trips[0];
      expect(Math.abs(off6 - 1.8e-3)).toBeLessThan(reclose);
      const period = (a.trips[a.trips.length - 1] - a.trips[0]) / (a.trips.length - 1);
      expect(1 / period).toBeGreaterThan(buzz[0]);
      expect(1 / period).toBeLessThan(buzz[1]);
      const off9 = b.closes[0] - b.trips[0];
      const off12 = c.closes[0] - c.trips[0];
      expect(off9).toBeGreaterThan(off6);
      expect(off12).toBeGreaterThan(off9);
      log(`${label}: re-close ${(off6 * 1e3).toFixed(3)} ms after the 6 V fire (target 1.8, residual ${((off6 - 1.8e-3) * 1e6).toFixed(1)} µs); buzz ${(1 / period).toFixed(0)} Hz; open ${(off9 * 1e3).toFixed(2)} ms at 9 V, ${(off12 * 1e3).toFixed(2)} ms at 12 V`);
    };
    // refit vibrator: the fit targets
    check(q6, q9, q12, 0.05e-3, [180, 215], 'refit vibrator');
    // vibrator in use (fitted before the re-ignition model; see the oracle fit test): ≈ 0.14 ms late, buzz ≈ 175 Hz
    check(r6, r9, r12, 0.16e-3, [165, 215], 'vibrator in use');
  });

  it('stored energy at the fire ½L₁I² ≈ 47 mJ on 6 V (48–50 mJ at 5.4–5.5 A, ECCT)', () => {
    expect(r6.tripE[0]).toBeCloseTo(0.5 * IG.coil.primaryInductance * r6.tripI[0] ** 2, 12);
    expect(r6.tripE[0]).toBeGreaterThan(44e-3);
    expect(r6.tripE[0]).toBeLessThan(50e-3);
    // 5.4–5.5 A → 48.1–49.9 mJ with L₁ = 3.3 mH
    expect(0.5 * IG.coil.primaryInductance * 5.45 ** 2).toBeCloseTo(49e-3, 3);
  });

  it('open circuit: ≈ 400 V on the points, ≥ 20 kV on the plug (Cool386 ≈ 300 V; Ford 8–20 kV in service)', () => {
    const coil = new TremblerCoil(IG.coil, PrimarySupply.dc(6, IG.battery.internalResistance), 0);
    const gap = new SparkGap(3e-3);
    gap.suppressBreakdown = true;
    let v1 = 0;
    let v2 = 0;
    for (let t = 0; t < 4.5e-3; t += 10e-6) {
      coil.step(10e-6, true, gap);
      if (coil.tripCount === 1) {
        v1 = Math.max(v1, Math.abs(coil.V1));
        v2 = Math.max(v2, Math.abs(coil.V2));
      }
    }
    expect(v1).toBeGreaterThan(300);
    expect(v1).toBeLessThan(450);
    expect(v2).toBeGreaterThan(20e3);
    log(`open circuit 6 V: ${v1.toFixed(0)} V on the points, ${(v2 / 1e3).toFixed(1)} kV on the plug`);
  });

  it('slow magneto pulses (120–150 rpm, 4–5 V peak) fire at 3.0–4.4 A (Kossor, HCCT)', () => {
    for (const rpm of [120, 150]) {
      // make at a rising EMF zero crossing so the coil sees a whole slow pulse
      const r = firstSpark(IG, rpm, IG.magneto.phaseDeg);
      expect(r.I).toBeGreaterThan(3.0);
      expect(r.I).toBeLessThan(4.4);
      log(`slow pulses ${rpm} rpm: fire at ${r.I.toFixed(2)} A after ${(r.trip * 1e3).toFixed(2)} ms`);
    }
  });
});

describe('MODEL_T_IGNITION — magneto phase fit: Patterson & Coniff 600 rpm spark ladder', () => {
  it('full retard fires at 26.5° ATDC; the lever plateaus fire at 4° ATDC, 18.5° BTDC, 41° BTDC (22.5° steps)', () => {
    const full = firstSpark(IG, 600, -IG.timer.advanceRangeDeg[0]);
    expect(Math.abs(full.deg - 26.5)).toBeLessThan(0.3);
    // scan the lever: first-spark angle vs timer make; each plateau's stationary (minimum) angle
    const minima: number[] = [];
    let prev = NaN;
    let cur = Infinity;
    for (let make = 15.5; make >= -64.5; make -= 0.5) {
      const f = firstSpark(IG, 600, make).deg;
      if (!Number.isNaN(prev) && f < prev - 10) {
        minima.push(cur);
        cur = Infinity;
      }
      cur = Math.min(cur, f);
      prev = f;
    }
    const target = [26.5, 4, -18.5, -41];
    for (let k = 0; k < target.length; k++) expect(Math.abs(minima[k] - target[k])).toBeLessThan(0.3);
    log(`ladder 600 rpm: full retard ${full.deg.toFixed(2)}° (26.5); plateaus ${minima.slice(0, 4).map((x) => x.toFixed(2)).join(', ')} (26.5, 4, −18.5, −41); residuals ${minima.slice(0, 4).map((x, k) => (x - target[k]).toFixed(2)).join(', ')}`);
  });
});

describe('Trembler-magneto ignition — spark timing vs speed and supply', () => {
  it('on a 6 V battery the first-spark lag ≈ 6·rpm·t_fire grows with rpm (t_fire ≈ const ≈ 3.5–4 ms)', () => {
    const lag: number[] = [];
    const tf: number[] = [];
    for (const rpm of [400, 1000, 1800]) {
      const r = firstSpark(IG, rpm, -25, 'battery');
      lag.push(r.deg + 25);
      tf.push(r.delay);
      expect(r.deg + 25).toBeCloseTo(6 * rpm * r.delay, 1);
    }
    for (const t of tf) {
      expect(t).toBeGreaterThan(3.4e-3);
      expect(t).toBeLessThan(4.0e-3);
    }
    expect(lag[2] / lag[0]).toBeGreaterThan(4.2);
    expect(lag[2]).toBeGreaterThan(35);
    log(`battery 6 V lag: ${lag.map((x) => x.toFixed(1)).join(' / ')}° at 400 / 1000 / 1800 rpm (t_fire ${tf.map((x) => (x * 1e3).toFixed(2)).join(' / ')} ms)`);
  });

  it('on the magneto the firing time shortens with speed (stronger, faster pulses: "automatic advance")', () => {
    const tf: number[] = [];
    const deg: number[] = [];
    for (const rpm of [400, 1000, 1800]) {
      const r = firstSpark(IG, rpm, -25, 'magneto');
      tf.push(r.delay);
      deg.push(r.deg);
    }
    expect(tf[2]).toBeLessThan(tf[0]);
    const bat1800 = firstSpark(IG, 1800, -25, 'battery').delay;
    expect(tf[2]).toBeLessThan(bat1800);
    log(`magneto, lever 25°: first spark ${deg.map((x) => x.toFixed(1)).join(' / ')}° after ${tf.map((x) => (x * 1e3).toFixed(2)).join(' / ')} ms at 400 / 1000 / 1800 rpm (battery 1800 rpm: ${(bat1800 * 1e3).toFixed(2)} ms)`);
  });

  it('timer command: make at −advance (clamped to the lever range), break after the 87° contact', () => {
    const c = tremblerTimerCommand(IG, 25);
    expect(c.dwellStartDeg).toBe(-25);
    expect(c.sparkDeg).toBe(62);
    expect(tremblerTimerCommand(IG, 100).dwellStartDeg).toBe(-64.5);
    expect(tremblerTimerCommand(IG, -40).dwellStartDeg).toBe(15.5);
  });
});

describe('Trembler-magneto ignition — one event = one timer contact, spark train, kernel', () => {
  /** Crank-driven run over one engine cycle at constant speed; returns per-step observations. */
  function cycleRun(rpm: number, source: 'magneto' | 'battery', gas: IgnitionGasState, advance = 25, ign?: IgnitionSystem) {
    const sys = ign ?? new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: source });
    const cmd = tremblerTimerCommand(IG, advance);
    const w = (rpm * 2 * Math.PI) / 60;
    const dt = 20e-6;
    let th = -359;
    const phases: string[] = [];
    let misDuring = false;
    let maxTrips = 0;
    let maxBds = 0;
    let timerSeen = false;
    let ledger = 0;
    for (let i = 0; i * rpm * 6 * dt < 720; i++) {
      th += rpm * 6 * dt;
      if (th >= 360) th -= 720;
      const s = sys.step(dt, th, cmd, gas, w);
      if (phases[phases.length - 1] !== s.phase) phases.push(s.phase);
      if (s.timerClosed) {
        timerSeen = true;
        if (s.misfire) misDuring = true;
      }
      maxTrips = Math.max(maxTrips, s.pointsBreakCount);
      maxBds = Math.max(maxBds, s.breakdownCount);
      ledger = Math.max(ledger, Math.abs(sys.coil.energyResidual()));
    }
    return { ign: sys, phases, misDuring, maxTrips, maxBds, timerSeen, ledger };
  }

  it('400 rpm on the battery: a buzz of sparks during the 87° contact, counted as ONE event; ledger exact', () => {
    const r = cycleRun(400, 'battery', gasT(0.3, 2, 1));
    const s = r.ign.state;
    expect(r.timerSeen).toBe(true);
    // 87° at 400 rpm = 36 ms of contact ≈ 7 buzz cycles
    expect(r.maxTrips).toBeGreaterThanOrEqual(5);
    // ONE spark (breakdown) per points opening, + at most the timer-break spark and a fresh breakdown after
    // a recovery window; each spark oscillates with the condenser ring, re-igniting at its current zeros
    // (review: previously every current zero was booked as a new breakdown, ≈ 18 per trip)
    expect(r.maxBds).toBeGreaterThanOrEqual(r.maxTrips);
    expect(r.maxBds).toBeLessThanOrEqual(r.maxTrips + 2);
    expect(s.reignitionCount).toBeGreaterThan(5 * s.breakdownCount);
    expect(s.breakdownCount).toBe(r.maxBds); // cumulative over the train (no reset per buzz)
    expect(r.phases.slice(0, 3)).toEqual(['off', 'charging', 'breakdown']);
    expect(r.phases).toContain('glow');
    // the phase no longer toggles 'breakdown' / 'charging' at every current zero of a spark (ring-up +
    // breakdown step per trip, + the timer-break spark)
    expect(r.phases.filter((p) => p === 'breakdown').length).toBeLessThanOrEqual(2 * r.maxTrips + 2);
    expect(r.phases[r.phases.length - 1]).toBe('done');
    expect(r.misDuring).toBe(false);
    expect(s.kernel.stage).toBe('handoff');
    expect(s.misfire).toBe(false);
    expect(s.timerClosed).toBe(false);
    expect(s.firstSparkDeg).toBeGreaterThan(-25);
    expect(s.firstSparkDeg).toBeLessThan(-25 + 6 * 400 * 4.5e-3);
    expect(s.conductingTime).toBeGreaterThan(1e-3);
    expect(s.sparkDuration).toBe(s.conductingTime);
    expect(r.ledger).toBeLessThan(1e-12);
    // the next timer make starts a new event
    const r2 = cycleRun(400, 'battery', gasT(0.3, 2, 1), 25, r.ign);
    expect(r2.ign.state.breakdownCount).toBe(r2.maxBds);
    expect(r2.maxBds).toBeLessThan(2 * r.maxBds);
  });

  it('a lean charge that one spark cannot ignite is ignited by the train (kernel kept alive and re-fed)', () => {
    const gas = gasT(0.2, 2, 1);
    const run = (contact: number) => {
      const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: 'battery' });
      let misDuring = false;
      for (let t = 0; t < 40e-3; t += 20e-6) {
        const s = ign.stepTimed(20e-6, t < contact, gas);
        if (t < contact && s.misfire) misDuring = true;
      }
      return { s: ign.state, misDuring };
    };
    const single = run(3.8e-3); // the timer breaks right after the first trip
    const train = run(30e-3);
    expect(single.s.pointsBreakCount).toBe(1);
    expect(single.s.misfire).toBe(true);
    expect(single.s.kernel.stage).toBe('quenched');
    expect(train.s.pointsBreakCount).toBeGreaterThanOrEqual(4);
    expect(train.misDuring).toBe(false);
    expect(train.s.kernel.stage).toBe('handoff');
    expect(train.s.misfire).toBe(false);
    log(`lean charge: single spark → ${single.s.kernel.quenchReason}; train (${train.s.pointsBreakCount} trips, ${train.s.breakdownCount} breakdowns) → hand-off ${(train.s.kernel.handoffTime * 1e3).toFixed(2)} ms after the kernel's start`);
  });

  it('the timer-break spark re-seeds a kernel quenched during the train; re-ignitions never re-seed (review finding)', () => {
    const blocking = gasT(0.05, 4, 4);
    const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: 'battery' });
    let s = ign.stepTimed(20e-6, true, blocking);
    let t = 20e-6;
    for (; t < 40e-3 && s.kernel.stage !== 'quenched'; t += 20e-6) s = ign.stepTimed(20e-6, true, blocking);
    expect(s.kernel.stage).toBe('quenched');
    expect(s.kernel.quenchReason).toBe('kernel did not reach hand-off in time');
    expect(s.misfire).toBe(false); // no verdict while the timer is closed
    // contact held until the points are closed with the current back up (no new trip yet): then the timer
    // breaks — ideal break, the coil fires once more ("timer spark") — into gas that can now ignite
    const bd0 = s.breakdownCount;
    for (; t < 60e-3 && (s.pointsOpen || s.primaryCurrent < 2.5); t += 20e-6) s = ign.stepTimed(20e-6, true, blocking);
    expect(s.breakdownCount).toBe(bd0);
    expect(s.pointsOpen).toBe(false);
    expect(s.kernel.stage).toBe('quenched');
    const ignitable = gasT(0.3, 2, 1);
    s = ign.stepTimed(20e-6, false, ignitable);
    for (let k = 0; k < 25 && s.breakdownCount === bd0; k++) s = ign.stepTimed(20e-6, false, ignitable);
    expect(s.breakdownCount).toBe(bd0 + 1);
    expect(s.kernel.stage).toBe('kernel'); // re-seeded (was: kept quenched → misfire)
    for (let k = 0; k < 500; k++) s = ign.stepTimed(20e-6, false, ignitable);
    expect(s.kernel.stage).toBe('handoff');
    expect(s.misfire).toBe(false);
    log(`timer-break re-seed: quenched at ${(t * 1e3).toFixed(2)} ms, timer spark → ${s.kernel.stage}`);

    // a kernel lost DURING a spark (forced here) is not re-seeded by that spark's re-ignitions, only by the
    // next breakdown (the next trip)
    const ign2 = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: 'battery' });
    s = ign2.stepTimed(20e-6, true, ignitable);
    for (let k = 0; k < 2000 && s.reignitionCount < 2; k++) s = ign2.stepTimed(20e-6, true, ignitable);
    expect(s.breakdownCount).toBe(1);
    (ign2.kernel as { stage: string }).stage = 'quenched';
    let reigSteps = 0;
    for (let k = 0; k < 5000 && s.breakdownCount === 1; k++) {
      const rg = s.reignitionCount;
      s = ign2.stepTimed(20e-6, true, ignitable);
      if (s.breakdownCount === 1 && s.reignitionCount > rg) {
        reigSteps++;
        expect(s.kernel.stage).toBe('quenched');
      }
    }
    expect(reigSteps).toBeGreaterThan(3);
    expect(s.breakdownCount).toBe(2);
    expect(s.kernel.stage).toBe('kernel');
  });

  it('a non-combustible charge: no misfire verdict while the timer is closed, misfire after the break', () => {
    const r = cycleRun(1000, 'magneto', gasT(0.3, 2, 1, 3e-3), 25);
    void r;
    const gas = { ...gasT(), SL: 0, expansionRatio: 1 };
    const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: 'magneto' });
    const out = cycleRun(1000, 'magneto', gas, 25, ign);
    expect(out.misDuring).toBe(false);
    expect(out.maxBds).toBeGreaterThan(0);
    expect(ign.state.misfire).toBe(true);
  });

  it('reports the snapshot train fields: timer, points, magneto EMF, supply EMF, condenser voltage', () => {
    const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: 'magneto' });
    const cmd = tremblerTimerCommand(IG, 25);
    const rpm = 1000;
    const w = (rpm * 2 * Math.PI) / 60;
    let th = -40;
    let sawOpen = false;
    let maxV1 = 0;
    for (let i = 0; i < 2000; i++) {
      th += rpm * 6 * 20e-6;
      const s = ign.step(20e-6, th, cmd, gasT(), w);
      expect(s.magnetoEmf).toBeCloseTo(magnetoEmf(IG.magneto, th, w), 9);
      expect(s.sourceEmf).toBeCloseTo(s.magnetoEmf, 9);
      expect(s.timerClosed).toBe(th >= -25 && th < 62);
      sawOpen ||= s.pointsOpen;
      maxV1 = Math.max(maxV1, Math.abs(s.primaryVoltage));
      expect(s.switchState).toBe(s.pointsOpen ? 'open' : 'closed');
    }
    expect(sawOpen).toBe(true);
    expect(maxV1).toBeGreaterThan(100);
    // battery: the supply EMF is the battery voltage (selection applied at the next make)
    ign.ignitionSource = 'battery';
    expect(ign.ignitionSource).toBe('battery');
    for (let i = 0; i < 40000; i++) {
      th += rpm * 6 * 20e-6;
      if (th >= 360) th -= 720;
      ign.step(20e-6, th, cmd, gasT(), w);
    }
    expect(ign.state.sourceEmf).toBe(IG.battery.voltage);
  });
});

describe('Trembler-magneto ignition — four cylinders, shared magneto, per-cylinder coils', () => {
  /**
   * Run all cylinders over two engine cycles at 1000 rpm; returns each cylinder's first-spark local angle
   * of its LAST event (the run may start inside a cylinder's contact window: that first event is partial).
   */
  function engineCycle(spec: EngineSpec, advance = 25): number[] {
    const sys = createIgnitionSystems(spec, { ignitionSource: 'magneto' });
    const cmd = tremblerTimerCommand(spec.ignition as TremblerMagnetoIgnitionSpec, advance);
    const rpm = 1000;
    const w = (rpm * 2 * Math.PI) / 60;
    const dt = 20e-6;
    const first = sys.map(() => NaN);
    const prevCount = sys.map(() => 0);
    let th = -359;
    for (let i = 0; i * rpm * 6 * dt < 2 * 720; i++) {
      th += rpm * 6 * dt;
      if (th >= 360) th -= 720;
      for (let c = 0; c < sys.length; c++) {
        const s = sys[c].step(dt, cylinderAngleDeg(spec, c, th), cmd, gasT(), w, th);
        if (prevCount[c] === 0 && s.breakdownCount > 0) first[c] = s.firstSparkDeg;
        prevCount[c] = s.breakdownCount;
      }
    }
    return first;
  }

  it('createIgnitionSystems: one system per cylinder with its firing offset; equal coils fire at equal local angles', () => {
    const sys = createIgnitionSystems(MODEL_T);
    expect(sys.length).toBe(4);
    expect(sys.map((s) => s.firingOffsetDeg)).toEqual([0, 180, 540, 360]);
    expect(sys.map((s) => s.cylinder)).toEqual([0, 1, 2, 3]);
    // 8 magneto cycles per revolution: every 180° firing offset is a whole number of EMF periods
    const f = engineCycle(MODEL_T);
    for (const x of f) expect(Math.abs(x - f[0])).toBeLessThan(1e-6);
  });

  it('a mis-adjusted vibrator (coilOverrides) fires later on its cylinder only', () => {
    // (the contract's override type intersects TremblerCoilSpec & {vibrator: Partial<…>}, which demands a
    // full vibrator; tremblerCoilOf merges partial vibrators field by field at run time)
    const ig: TremblerMagnetoIgnitionSpec = { ...IG, coilOverrides: [{}, {}, { vibrator: { ...IG.coil.vibrator, pullCurrent: 4.2 } }] };
    expect(tremblerCoilOf(ig, 2).vibrator.pullCurrent).toBe(4.2);
    expect(tremblerCoilOf(ig, 2).vibrator.naturalFrequency).toBe(IG.coil.vibrator.naturalFrequency);
    expect(tremblerCoilOf(ig, 0)).toEqual(IG.coil);
    expect(tremblerCoilOf(IG, 0)).toBe(IG.coil);
    // a genuinely partial vibrator override is merged field by field
    const partial = { ...IG, coilOverrides: [{ condenserCapacitance: 0.3e-6, vibrator: { pullCurrent: 4.2 } }] } as unknown as TremblerMagnetoIgnitionSpec;
    expect(tremblerCoilOf(partial, 0)).toEqual({ ...IG.coil, condenserCapacitance: 0.3e-6, vibrator: { ...IG.coil.vibrator, pullCurrent: 4.2 } });
    const f = engineCycle({ ...MODEL_T, ignition: ig });
    expect(f[2]).toBeGreaterThan(f[0] + 0.5);
    expect(Math.abs(f[1] - f[0])).toBeLessThan(1e-6);
    expect(Math.abs(f[3] - f[0])).toBeLessThan(1e-6);
    log(`coil mismatch: cylinder 3 (I_p 4.2 A) fires at ${f[2].toFixed(2)}° vs ${f[0].toFixed(2)}° (Δ ${(f[2] - f[0]).toFixed(2)}° at 1000 rpm; HCCT-set coil sets spread ≈ 4° [Kossor])`);
  });
});

describe('Trembler-magneto ignition — robustness of the re-ignition events', () => {
  it('seeded random sweep (gas, flow, speed, lever, supply, step): finite, exact ledger, bounded events', () => {
    let seed = 20261003;
    const rnd = (): number => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    const lerp = (a: number, b: number): number => a + (b - a) * rnd();
    const logu = (a: number, b: number): number => Math.exp(lerp(Math.log(a), Math.log(b)));
    for (let trial = 0; trial < 16; trial++) {
      const p = logu(0.3e5, 30e5);
      const Tu = lerp(300, 900);
      const SL = rnd() < 0.2 ? 0 : logu(0.02, 2);
      const nu = C45.nuU * (Tu / C45.Tu) ** 1.7 * (C45.p / p);
      const lF = SL > 0 ? nu / SL : 1e-5;
      const gas: IgnitionGasState = {
        ...gasT(),
        p,
        Tu,
        rhoU: (p * 0.0298) / (8.314 * Tu),
        SL,
        marksteinLength: lerp(-2, 8) * lF,
        flameThickness: lF,
        uPrime: logu(0.05, 10),
        kinematicViscosity: nu,
        flowVelocity: rnd() < 0.4 ? 0 : logu(0.1, 30),
      };
      const rpm = logu(150, 2500);
      const source = rnd() < 0.5 ? 'magneto' : 'battery';
      const make = -lerp(-15, 64);
      const dt = logu(5e-6, 200e-6);
      const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: source });
      const cmd = tremblerTimerCommand(IG, -make);
      const w = (rpm * 2 * Math.PI) / 60;
      let th = make - 5;
      let finite = true;
      const t0 = performance.now();
      ign.step(dt, th, cmd, gas, w);
      while (th < make + IG.timer.contactArcDeg + 30) {
        th += rpm * 6 * dt;
        const s = ign.step(dt, th >= 360 ? th - 720 : th, cmd, gas, w);
        finite &&= Number.isFinite(s.secondaryVoltage) && Number.isFinite(s.energyToGas) && Number.isFinite(s.kernel.radius);
      }
      const ms = performance.now() - t0;
      const s = ign.state;
      const c = ign.trembler!;
      expect(finite).toBe(true);
      expect(Math.abs(c.energyResidual())).toBeLessThan(1e-11);
      expect(s.energyToGas + s.energyToElectrodes + s.energyRadiated).toBeCloseTo(s.energyDelivered, 12);
      expect(s.breakdownCount).toBeLessThan(1000);
      expect(s.reignitionCount).toBeLessThan(5000);
      expect(c.forcedSubsteps).toBeLessThan(50); // no zero-length event loops
      expect(ms).toBeLessThan(2000);
      log(`sweep ${trial}: ${source} ${rpm.toFixed(0)} rpm, ${(p / 1e5).toFixed(1)} bar ${Tu.toFixed(0)} K, flow ${(gas.flowVelocity ?? 0).toFixed(1)} m/s, dt ${(dt * 1e6).toFixed(0)} µs: ${s.pointsBreakCount} trips, ${s.breakdownCount} bds, ${s.reignitionCount} re-ignitions, forced ${c.forcedSubsteps}, ${ms.toFixed(0)} ms`);
    }
  });
});

describe('Trembler-magneto ignition — CPU per ignition event', () => {
  it('one cylinder over a full 720° cycle (20 µs engine steps) at 400 / 1000 / 1800 rpm', () => {
    const rows: string[] = [];
    let worst = 0;
    for (const source of ['magneto', 'battery'] as const) {
      for (const rpm of [400, 1000, 1800]) {
        const runOnce = () => {
          const ign = new IgnitionSystem(IG, MODEL_T.sparkPlug, { ignitionSource: source });
          const cmd = tremblerTimerCommand(IG, 25);
          const w = (rpm * 2 * Math.PI) / 60;
          let th = -359;
          for (let i = 0; i * rpm * 6 * 20e-6 < 720; i++) {
            th += rpm * 6 * 20e-6;
            if (th >= 360) th -= 720;
            ign.step(20e-6, th, cmd, gasT(), w);
          }
          return ign.state;
        };
        runOnce(); // warm-up
        const n = 5;
        const t0 = performance.now();
        let s = runOnce();
        for (let k = 1; k < n; k++) s = runOnce();
        const ms = (performance.now() - t0) / n;
        worst = Math.max(worst, ms);
        rows.push(`${source.padEnd(8)} ${String(rpm).padStart(4)} rpm: ${ms.toFixed(2)} ms/event (${s.pointsBreakCount} trips, ${s.breakdownCount} breakdowns, ${s.reignitionCount} re-ignitions, ${(s.energyToGas * 1e3).toFixed(1)} mJ to the gas, cycle ${(120 / rpm * 1e3).toFixed(0)} ms)`);
      }
    }
    log('\nCPU per trembler ignition event (incl. the 720° of engine steps):\n' + rows.join('\n'));
    expect(worst).toBeLessThan(60);
  });
});
