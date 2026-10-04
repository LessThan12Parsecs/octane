/**
 * Code-review regressions of the cycle model (Model T review, 'cycle' and 'gasx-mech-load' dimensions):
 *  - the load sanitiser and LoadModel accept only a vehicle's OWN gears (no Object.prototype names), and
 *    CycleModel.setOperatingPoint validates the load before committing anything (no poisoned state);
 *  - the dense snapshot times honour stepUntil's relative time tolerance: no stalled snapshot stream at a
 *    large model time (the reviewer's stall reproduced at t ≈ 3e4 s);
 *  - a compression-ratio change of a cylinder i ≥ 1 at its local wrap is booked in neither the cylinder nor
 *    the ENGINE summary (the CFR geometry × 4 with its variable CR).
 * Runtime ≈ 10–20 s (vitest).
 */
import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../core/engine-spec';
import type { EngineSnapshot } from '../core/snapshot';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_FORD_DYNO, MODEL_T_VEHICLE } from '../engines/model-t';
import { LoadModel } from '../mechanics/load';
import { CycleModel, Cylinder } from './cycle-model';
import { DENSE_SPARK_AFTER_BREAKDOWN, DENSE_SPARK_DT, DENSE_TRAIN_DT, DENSE_TRIP_WINDOW, EngineSimulator } from './engine-simulator';
import { sanitizeOperatingPoint } from './options';

const INHERITED = ['toString', 'constructor', '__proto__', 'hasOwnProperty', 'valueOf'];

describe('load sanitiser / LoadModel / setOperatingPoint: only the vehicle\'s own gears', () => {
  it('sanitizeOperatingPoint replaces an inherited gear name by the fallback load', () => {
    const fb = { ...MODEL_T_CRUISE, load: { kind: 'vehicle' as const, gear: 'low', grade: 0.02 } };
    for (const gear of INHERITED) {
      const s = sanitizeOperatingPoint(MODEL_T, { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear, grade: 0 } }, fb);
      expect(s.load, gear).toEqual({ kind: 'vehicle', gear: 'low', grade: 0.02 });
    }
    // own gears and neutral are kept; a gear with a non-finite or non-positive ratio is refused
    for (const gear of ['low', 'high', 'reverse', 'neutral']) {
      expect(sanitizeOperatingPoint(MODEL_T, { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear, grade: 0 } }).load).toEqual({ kind: 'vehicle', gear, grade: 0 });
    }
    const odd = { ...MODEL_T, vehicle: { ...MODEL_T_VEHICLE, gears: { low: 2.75, high: 1, broken: Number.NaN, zero: 0 } } };
    for (const gear of ['broken', 'zero']) {
      expect(sanitizeOperatingPoint(odd, { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear, grade: 0 } }, fb).load).toEqual(fb.load);
    }
  });

  it('LoadModel.configure rejects an inherited gear name as an unknown gear', () => {
    for (const gear of INHERITED) {
      expect(() => new LoadModel(MODEL_T_VEHICLE).configure({ ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear, grade: 0 } })).toThrow(/unknown gear/);
    }
  });

  it('setOperatingPoint: an inherited gear falls back to the pending load; a load that cannot be configured throws BEFORE anything is committed', () => {
    // (a) the sanitiser repairs the malformed gear: the model keeps its gear and keeps running
    const m = new CycleModel(MODEL_T, { ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'low', grade: 0 } }, { warmupCycles: 0 });
    m.stepUntil(Infinity, -300);
    for (const gear of INHERITED) {
      expect(() => m.setOperatingPoint({ load: { kind: 'vehicle', gear, grade: 0 } })).not.toThrow();
      expect(m.op.load).toEqual({ kind: 'vehicle', gear: 'low', grade: 0 });
      expect(m.load.gear).toBe('low');
    }
    // (b) a load that passes the sanitiser but cannot be configured (invalid vehicle data): the call throws, and
    // pendingOp, op, the load model and the state stay as they were — later patches and steps still work
    const bad: EngineSpec = { ...CFR_F1, vehicle: { ...MODEL_T_VEHICLE, drivelineEfficiency: 0 } };
    const c = new CycleModel(bad, CFR_RON_CONDITIONS, { warmupCycles: 0 });
    c.stepUntil(Infinity, -300);
    const pending = JSON.stringify(c.pendingOp);
    const op = JSON.stringify(c.op);
    const y = Array.from(c.y);
    expect(() => c.setOperatingPoint({ speedMode: 'free', load: { kind: 'vehicle', gear: 'high', grade: 0 } })).toThrow(RangeError);
    expect(JSON.stringify(c.pendingOp)).toBe(pending);
    expect(JSON.stringify(c.op)).toBe(op);
    expect(c.load.kind).toBe('constant');
    expect(Array.from(c.y)).toEqual(y);
    expect(() => c.setOperatingPoint({ throttle: 0.5 })).not.toThrow();
    expect(c.op.throttle).toBe(0.5);
    c.runCycles(1);
    expect(c.cycle).toBe(1);
  });
});

/**
 * Offset every model clock by T0 (the state after T0 seconds of running): the engine time and its cycle start, and
 * each cylinder's time stamps (tSparkCmd, tTimerBreak, tKnockOnset, …). The ignition systems keep their own clocks
 * (EngineSimulator converts with the offset model time − ignition time).
 */
function offsetClock(m: CycleModel, T0: number): void {
  m.t += T0;
  (m as unknown as { tEngineStart: number }).tEngineStart += T0;
  for (const c of m.cylinders) c.shiftTimes(T0);
}

/** n snapshots; returns their times. */
function times(sim: EngineSimulator, n: number): number[] {
  const out: number[] = [];
  for (let i = 0; i < n; i++) out.push(sim.advanceToNextSnapshot().t);
  return out;
}

const T0 = 3e4;

describe('dense snapshot times at a large model time (t = 3e4 s): no stalled stream', () => {
  // stepUntil returns without stepping when t ≥ t_target − 1e-14·t_target (3e-10 s at 3e4 s). The dense windows
  // used absolute bounds (1e-12 s window ends, 1e-6 of the grid interval = 5e-12 s for the 5 µs grid), so a step
  // ending in the band between them made nextDenseTime return a target stepUntil treats as reached: the same
  // snapshot came back on every call (reproduced by the reviewer after 48 Model T cycles at T0 = 3e4 s).
  it('CFR (inductive): t just below a 5 µs grid point and just below the end of the fine window', () => {
    const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 2, warmupCycles: 0 });
    const m = sim.model;
    const c = m.c0;
    let s: EngineSnapshot = sim.advanceToNextSnapshot();
    while (Number.isNaN(c.tSparkCmd)) s = sim.advanceToNextSnapshot();
    // inside the 5 µs window after the switch-off; wait for the breakdown (the window then ends 100 µs later)
    while (Number.isNaN(c.ign!.state.breakdownDelay)) s = sim.advanceToNextSnapshot();
    expect(s.t).toBeGreaterThan(c.tSparkCmd);
    offsetClock(m, T0);
    const ts = c.tSparkCmd;
    // (1) 1e-10 s below the next 5 µs grid point (inside the band (5e-12, 3e-10) s)
    const k = Math.floor((m.t - ts) / DENSE_SPARK_DT) + 1;
    const tGrid = ts + k * DENSE_SPARK_DT - 1e-10;
    m.t = tGrid;
    let tt = times(sim, 4);
    expect(tt[0]).toBeGreaterThan(tGrid);
    for (let i = 1; i < 4; i++) expect(tt[i]).toBeGreaterThan(tt[i - 1]);
    // (2) 1e-10 s below the end of the fine window (switch-off + breakdown delay + 100 µs)
    const tFineEnd = ts + c.ign!.state.breakdownDelay + DENSE_SPARK_AFTER_BREAKDOWN;
    expect(m.t).toBeLessThan(tFineEnd);
    m.t = tFineEnd - 1e-10;
    tt = times(sim, 4);
    expect(tt[0]).toBeGreaterThan(tFineEnd - 1e-10);
    for (let i = 1; i < 4; i++) expect(tt[i]).toBeGreaterThan(tt[i - 1]);
    // the stream then runs on through the rest of the cycle and the wrap at the large time
    let s2 = sim.advanceToNextSnapshot();
    let prev = s2.t;
    while (s2.cycle < 1) {
      s2 = sim.advanceToNextSnapshot();
      expect(s2.t).toBeGreaterThan(prev);
      prev = s2.t;
    }
  });

  it('Model T (trembler): t just below a 100 µs spark-train grid point and just below the end of the timer-break window', () => {
    const sim = new EngineSimulator(MODEL_T, MODEL_T_FORD_DYNO, { snapshotEveryDeg: 2, warmupCycles: 0 });
    const m = sim.model;
    // a cylinder in its spark train (timer contact made in this run: the train grid is anchored at tSparkCmd)
    let c: Cylinder | undefined;
    let n = 0;
    while (!c && n++ < 100_000) {
      sim.advanceToNextSnapshot();
      c = m.cylinders.find((x) => x.dwellSeen && !Number.isNaN(x.tSparkCmd) && Number.isNaN(x.tTimerBreak) && x.ign!.state.trainActive);
    }
    expect(c).toBeDefined();
    offsetClock(m, T0);
    // (1) 1e-10 s below the next train grid point (model-time anchored: tSparkCmd + k·100 µs)
    const ts = c!.tSparkCmd;
    const k = Math.floor((m.t - ts) / DENSE_TRAIN_DT) + 1;
    const tGrid = ts + k * DENSE_TRAIN_DT - 1e-10;
    m.t = tGrid;
    let tt = times(sim, 4);
    expect(tt[0]).toBeGreaterThan(tGrid);
    for (let i = 1; i < 4; i++) expect(tt[i]).toBeGreaterThan(tt[i - 1]);
    // (2) the timer break of that contact (a snapshot event), then 1e-10 s below the end of its 150 µs window
    n = 0;
    while (Number.isNaN(c!.tTimerBreak) && n++ < 100_000) sim.advanceToNextSnapshot();
    const tb = c!.tTimerBreak;
    expect(tb).toBeGreaterThan(T0);
    expect(m.t - tb).toBeLessThan(DENSE_TRIP_WINDOW);
    const tEnd = tb + DENSE_TRIP_WINDOW - 1e-10;
    m.t = tEnd;
    tt = times(sim, 4);
    expect(tt[0]).toBeGreaterThan(tEnd);
    for (let i = 1; i < 4; i++) expect(tt[i]).toBeGreaterThan(tt[i - 1]);
  });
});

/** CFR F-1 geometry with n cylinders and the network scaled by n (as multicylinder.test.ts). */
function cfrMulti(offsets: number[]): EngineSpec {
  const n = offsets.length;
  const mf = CFR_F1.manifolds;
  return {
    ...CFR_F1,
    id: undefined,
    name: `CFR F-1 geometry × ${n} (test)`,
    cylinders: n,
    layout: {
      firingOrder: offsets.map((_, i) => i + 1),
      firingOffsetDeg: offsets,
      axisZ: offsets.map((_, i) => 0.12 * (i - (n - 1) / 2)),
      mirrorZ: offsets.map(() => false),
      mainBearingZ: [],
    },
    manifolds: {
      ...mf,
      intakeVolume: mf.intakeVolume * n,
      exhaustVolume: mf.exhaustVolume * n,
      throttleDiameter: mf.throttleDiameter * Math.sqrt(n),
      exhaustOutletDiameter: mf.exhaustOutletDiameter * Math.sqrt(n),
    },
  };
}

describe('compression-ratio change of cylinders i ≥ 1: in neither summary (CFR geometry × 4, CR 6 → 8)', () => {
  it('the engine summary of the cycle with the local-wrap CR changes excludes their isentropic work', () => {
    const proto = Cylinder.prototype as unknown as { changeCompressionRatio(this: Cylinder, cr: number): void };
    const orig = proto.changeCompressionRatio;
    const booked: { cyl: number; engineCycle: number; dW: number }[] = [];
    proto.changeCompressionRatio = function (this: Cylinder, cr: number): void {
      const w0 = this.e.y[this.ix.W];
      orig.call(this, cr);
      booked.push({ cyl: this.index, engineCycle: this.e.cycle, dW: this.e.y[this.ix.W] - w0 });
    };
    try {
      const m = new CycleModel(cfrMulti([0, 180, 540, 360]), { ...CFR_RON_CONDITIONS, compressionRatio: 6 }, { warmupCycles: 1 });
      m.runCycles(1);
      expect(booked).toEqual([]);
      m.setOperatingPoint({ compressionRatio: 8 });
      // cycle K: the new CR is applied at its END (cylinder 0 at the engine wrap, after the engine summary)
      const LK0 = m.workLedger();
      const eK = m.runCycles(1).find((s) => s.engine)!.engine!;
      const LK1 = m.workLedger();
      expect(booked.map((b) => b.cyl)).toEqual([0]);
      const w0 = booked[0].dW;
      expect(Math.abs(w0)).toBeGreaterThan(1); // ≈ −4.4 J
      expect(eK.indicatedTorque * 4 * Math.PI).toBeCloseTo(LK1 - w0 - LK0, 6);
      // cycle K+1: cylinders 1–3 change at their local wraps INSIDE the engine cycle
      const eK1 = m.runCycles(1).find((s) => s.engine)!.engine!;
      const LK2 = m.workLedger();
      const inside = booked.filter((b) => b.cyl > 0);
      expect(inside.map((b) => b.cyl).sort()).toEqual([1, 2, 3]);
      const wInside = inside.reduce((a, b) => a + b.dW, 0);
      expect(Math.abs(wInside)).toBeGreaterThan(3 * 1); // ≈ −13 J
      // engine ∫p dV of cycle K+1 = ledger change − the CR-change work booked inside it
      expect(eK1.indicatedTorque * 4 * Math.PI).toBeCloseTo(LK2 - LK1 - wInside, 6);
      expect(m.cylinders.every((x) => x.kin.compressionRatio === 8)).toBe(true);
    } finally {
      proto.changeCompressionRatio = orig;
    }
  });
});
