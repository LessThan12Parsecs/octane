/**
 * Ford Model T against its period data (calibration phase, 2026-10): the calibrated simulator defaults
 * (cycle/calibration.ts MODEL_T_CALIBRATION, options.ts MODEL_T_CALIBRATED_FRICTION), public API only.
 *
 * [SLOW: ≈ 3–6 min. Runs only with OCTANE_VALIDATION=1, e.g.
 *    OCTANE_VALIDATION=1 npx vitest run test/validation/measured_modelt_data.test.ts ]
 *
 * Targets (test/fixtures/modelt_period_data.json, grades and caveats in modelt_validation_README.md):
 *  1. Ford's WOT brake-torque table (FSB Fig. 84, transmission output): the model's MAXIMUM brake torque over
 *     the spark-lever range at each speed (magneto; the lever→first-spark map is a magneto staircase, so the
 *     lever is scanned) at Ford's dyno condition MODEL_T_FORD_DYNO. 900–1800 rpm within ±10 % passes;
 *     500–700 rpm (+19…+29 %) is a documented disagreement (`it.fails`, README).
 *  2. Upton (J. SAE Aug. 1923) MBT a₀ = 0.108R/(1 + 0.001R) at his test condition (MODEL_T_UPTON_CONDITIONS:
 *     140 °F air and water, 29.35 inHg, an ideal timed spark = his Atwater Kent battery distributor, simulated
 *     with the CFR inductive coil on the Model T spec): model MBT within 6° at 400 / 800 / 1200 rpm; his
 *     Table 3 wide-open intake suction (MODEL_T_UPTON_TABLE3) at 800 / 1200 / 1400 rpm within ±15 %; his Fig. 4
 *     peak BMEP ≈ 70 psi.
 *  3. Motored compression 55–60 psig at cranking speed: a documented residual (+13 %, no blow-by; `it.fails`).
 *  4. Top speed 42–45 mph in high gear on a level road (fixed-speed WOT torque vs the 'vehicle' road load);
 *     the cruise preset (throttle 0.35, lever 50) settles at 25–35 mph without knock.
 *  5. Knock: none at the cruise settings; spark knock at WOT and low speed with the lever fully advanced
 *     (Ricardo: side-valve engines at ≈ 4:1 on 45–50 ON fuel detonated unless the spark was retarded).
 *  6. Idle preset ≈ 400 rpm declutched (free speed).
 * `it.fails` marks a target the model does not meet — flip to `it` when a model change makes it pass.
 */
import { beforeAll, describe, expect, it } from 'vitest';
import period from '../fixtures/modelt_period_data.json';
import { CycleModel } from '../../src/physics/cycle/index';
import type { EngineSpec } from '../../src/physics/core/engine-spec';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary, EngineCycleSummary } from '../../src/physics/core/snapshot';
import { CFR_F1 } from '../../src/physics/engines/cfr';
import { ENGINES } from '../../src/physics/engines';
import {
  MODEL_T,
  MODEL_T_CRUISE,
  MODEL_T_FORD_DYNO,
  MODEL_T_FORD_WOT_TABLE,
  MODEL_T_UPTON_CONDITIONS,
  MODEL_T_UPTON_TABLE3,
  modelTUptonMbtDeg,
} from '../../src/physics/engines/model-t';
import { LoadModel } from '../../src/physics/mechanics/load';

const RUN = !!process.env.OCTANE_VALIDATION;
const LBFT = 1.3558179483314004;
const PSI = 6894.757293168;
const INHG = 3386.389;
const MPH = 0.44704;
const VD_TOTAL = 4 * (Math.PI / 4) * MODEL_T.geometry.bore ** 2 * MODEL_T.geometry.stroke;

interface Run {
  e: EngineCycleSummary;
  s: CycleSummary[];
  /** cycle-mean intake-manifold pressure of the last cycle, Pa (NaN unless traced) */
  pIntMean: number;
}

/** Fixed-speed run: warm-up + 2 engine cycles, summaries of the last; optional trace of the last cycle. */
function run(spec: EngineSpec, op: OperatingPoint, trace = false, opts: Record<string, unknown> = {}): Run {
  const m = new CycleModel(spec, op, { warmupCycles: 2, ...opts });
  m.runCycles(1);
  if (trace) m.recordTrace();
  const s = m.runCycles(1);
  const e = s.find((x) => x.engine)!.engine!;
  let pIntMean = NaN;
  if (trace && m.trace) {
    const t = m.trace.t;
    const p = m.trace.pIntake;
    let a = 0;
    for (let k = 1; k < t.length; k++) a += 0.5 * (p[k] + p[k - 1]) * (t[k] - t[k - 1]);
    pIntMean = a / (t[t.length - 1] - t[0]);
  }
  return { e, s, pIntMean };
}

const ford = (rpm: number): number => MODEL_T_FORD_WOT_TABLE.find((r) => r[0] === rpm)![1];
const UPTON_SPEC: EngineSpec = { ...MODEL_T, ignition: CFR_F1.ignition } as EngineSpec;
const knocks = (s: CycleSummary[]): number => s.filter((x) => Number.isFinite(x.knockOnsetDeg)).length;
const mapo = (s: CycleSummary[]): number => Math.max(...s.map((x) => x.mapo)) / 1e5;

/** Least-squares parabola vertex of y(x) (x ascending, ≥ 3 points). */
function vertex(x: number[], y: number[]): number {
  const n = x.length;
  const xm = x.reduce((a, b) => a + b, 0) / n;
  let s0 = 0, s1 = 0, s2 = 0, s3 = 0, s4 = 0, t0 = 0, t1 = 0, t2 = 0;
  for (let k = 0; k < n; k++) {
    const u = x[k] - xm;
    s0 += 1; s1 += u; s2 += u * u; s3 += u ** 3; s4 += u ** 4;
    t0 += y[k]; t1 += u * y[k]; t2 += u * u * y[k];
  }
  // normal equations for y = a + b u + c u²
  const M = [[s0, s1, s2], [s1, s2, s3], [s2, s3, s4]];
  const v = [t0, t1, t2];
  const det = (A: number[][]) => A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  const D = det(M);
  const col = (j: number) => M.map((row, i) => row.map((c, k) => (k === j ? v[i] : c)));
  const b = det(col(1)) / D;
  const c = det(col(2)) / D;
  return xm - b / (2 * c);
}

describe.skipIf(!RUN)('Ford Model T: period data [SLOW]', () => {
  /** WOT lever scan per speed: the lever set covers the torque maximum of the magneto staircase (calibration scans). */
  const WOT_SPEEDS = [500, 700, 900, 1100, 1300, 1500, 1600, 1700, 1800];
  const leversFor = (rpm: number) => (rpm <= 900 ? [42, 47, 54.5, 62, 64.5] : [62, 64.5]);
  const wot = new Map<number, { lever: number; Tb: number; spark: number; knock: number; mapo: number }[]>();
  const best = new Map<number, { lever: number; Tb: number; spark: number; knock: number; mapo: number }>();
  const UPTON_SPEEDS = [400, 800, 1200];
  const mbt = new Map<number, { spark: number; Tb: number; bmepPsi: number }[]>();
  const suction = new Map<number, number>();
  let motored: Run;
  let top1900: Run;
  const cruise: { rpm: number; r: Run }[] = [];
  let idleRpm: number[] = [];
  let idleMisfires = 0;

  beforeAll(() => {
    for (const rpm of WOT_SPEEDS) {
      const rows = leversFor(rpm).map((lever) => {
        const r = run(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm, sparkAdvanceDeg: lever });
        const sp = r.s.reduce((a, x) => a + (x.sparkDeg ?? NaN), 0) / r.s.length;
        return { lever, Tb: r.e.brakeTorque / LBFT, spark: sp, knock: knocks(r.s), mapo: mapo(r.s) };
      });
      wot.set(rpm, rows);
      best.set(rpm, rows.reduce((a, b) => (b.Tb > a.Tb ? b : a)));
    }
    for (const rpm of UPTON_SPEEDS) {
      const a0 = modelTUptonMbtDeg(rpm);
      const rows = [-10, -5, 0, 5, 10].map((d) => {
        const r = run(UPTON_SPEC, { ...MODEL_T_UPTON_CONDITIONS, rpm, sparkAdvanceDeg: a0 + d }, d === 0);
        if (d === 0) suction.set(rpm, (MODEL_T_UPTON_CONDITIONS.ambientPressure - r.pIntMean) / INHG);
        return { spark: a0 + d, Tb: r.e.brakeTorque / LBFT, bmepPsi: r.e.bmep / PSI };
      });
      mbt.set(rpm, rows);
    }
    for (const rpm of [600, 1400]) {
      const r = run(UPTON_SPEC, { ...MODEL_T_UPTON_CONDITIONS, rpm, sparkAdvanceDeg: modelTUptonMbtDeg(rpm) }, true);
      suction.set(rpm, (MODEL_T_UPTON_CONDITIONS.ambientPressure - r.pIntMean) / INHG);
    }
    motored = run(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm: 150 }, false, { combustionModel: 'none', warmupCycles: 3 });
    top1900 = run(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm: 1900, sparkAdvanceDeg: 64.5 });
    for (const rpm of [1100, 1400]) cruise.push({ rpm, r: run(MODEL_T, { ...MODEL_T_CRUISE, speedMode: 'fixed', load: undefined, rpm }) });
    // idle preset, free speed in neutral
    const idle = ENGINES['ford-model-t'].presets.find((p) => p.id === 'idle')!;
    const im = new CycleModel(MODEL_T, { ...MODEL_T_CRUISE, ...idle.op } as OperatingPoint, { warmupCycles: 1 });
    for (let c = 0; c < 12; c++) {
      const s = im.runCycles(1);
      idleRpm.push(s.find((x) => x.engine)!.engine!.rpmMean);
      if (c >= 8) idleMisfires += s.filter((x) => x.misfire).length;
    }
    // log the comparison tables
    const lines = WOT_SPEEDS.map((rpm) => {
      const b = best.get(rpm)!;
      return `${rpm} rpm: max ${b.Tb.toFixed(1)} lb-ft at lever ${b.lever} (first spark ${(-b.spark).toFixed(1)}° BTDC, knock ${b.knock}/4) vs Ford ${ford(rpm)} (${((b.Tb / ford(rpm) - 1) * 100).toFixed(1)} %)`;
    });
    console.log(`[Model T] Ford WOT (max over lever):\n  ${lines.join('\n  ')}`);
  }, 1_800_000);

  describe('1. Ford WOT brake-torque table (FSB Fig. 84, grade B)', () => {
    it('the fixture and the spec carry the same table', () => {
      const ds = period.datasets.find((d) => d.key === 'ford_wot_1918')!;
      for (const [rpm, lbft] of MODEL_T_FORD_WOT_TABLE) expect((ds.data as { rpm: number; torqueLbFt: number }[]).find((r) => r.rpm === rpm)!.torqueLbFt).toBe(lbft);
    });
    it('900–1800 rpm: maximum brake torque over the lever range within ±10 % of Ford', () => {
      for (const rpm of WOT_SPEEDS.filter((r) => r >= 900)) {
        const e = best.get(rpm)!.Tb / ford(rpm) - 1;
        expect(Math.abs(e), `${rpm} rpm: ${(e * 100).toFixed(1)} %`).toBeLessThan(0.1);
      }
    });
    it.fails('500–700 rpm within ±10 % of Ford [model +19…+29 %: no mixture-preparation / blow-by physics, knock-limited spark in period practice — README]', () => {
      for (const rpm of [500, 700]) expect(Math.abs(best.get(rpm)!.Tb / ford(rpm) - 1)).toBeLessThan(0.1);
    });
    it('shape: peak torque at 600–1100 rpm, falling monotonically above 1100 rpm to 50–70 % of the peak at 1800 rpm (Ford 64 %)', () => {
      const peak = [...best.entries()].reduce((a, b) => (b[1].Tb > a[1].Tb ? b : a));
      expect(peak[0]).toBeGreaterThanOrEqual(600);
      expect(peak[0]).toBeLessThanOrEqual(1100);
      const hi = WOT_SPEEDS.filter((r) => r >= 1100).map((r) => best.get(r)!.Tb);
      for (let k = 1; k < hi.length; k++) expect(hi[k]).toBeLessThan(hi[k - 1]);
      const ratio = best.get(1800)!.Tb / peak[1].Tb;
      expect(ratio).toBeGreaterThan(0.5);
      expect(ratio).toBeLessThan(0.7);
    });
    it('peak brake power 18–22.5 hp at 1500–1700 rpm (Ford: 20 hp at 1500–1600, best engines 22.5 hp)', () => {
      const hp = [1500, 1600, 1700].map((rpm) => (best.get(rpm)!.Tb * rpm) / 5252.1);
      const p = Math.max(...hp);
      expect(p).toBeGreaterThan(18);
      expect(p).toBeLessThan(22.5);
    });
    it('the stock ignition cannot reach MBT above ≈ 1000 rpm: the torque maximum is at the lever stop from 1100 rpm', () => {
      for (const rpm of WOT_SPEEDS.filter((r) => r >= 1100 && r <= 1800)) expect(best.get(rpm)!.lever).toBe(64.5);
    });
  });

  describe("2. Upton 1923 (MBT relation grade C; Table 3 suction, Fig. 4 BMEP)", () => {
    it('best-torque spark advance within 6° of a₀ = 0.108R/(1 + 0.001R) at 400 / 800 / 1200 rpm', () => {
      for (const rpm of UPTON_SPEEDS) {
        const rows = mbt.get(rpm)!;
        const v = vertex(rows.map((r) => r.spark), rows.map((r) => r.Tb));
        console.log(`[Model T] Upton MBT ${rpm} rpm: model ${v.toFixed(1)}° vs ${modelTUptonMbtDeg(rpm).toFixed(1)}°`);
        expect(Math.abs(v - modelTUptonMbtDeg(rpm))).toBeLessThan(6);
      }
    });
    it('wide-open intake suction within ±15 % of Table 3 at 800 / 1200 / 1400 rpm', () => {
      for (const rpm of [800, 1200, 1400]) {
        const meas = MODEL_T_UPTON_TABLE3.filter((r) => r[0] === rpm).reduce((a, r) => Math.min(a, r[1]), Infinity);
        console.log(`[Model T] Upton suction ${rpm} rpm: model ${suction.get(rpm)!.toFixed(2)} vs ${meas} inHg`);
        expect(Math.abs(suction.get(rpm)! / meas - 1)).toBeLessThan(0.15);
      }
    });
    it.fails('wide-open intake suction within ±15 % of Table 3 at 600 rpm [model 1.61 vs 1.10 inHg]', () => {
      expect(Math.abs(suction.get(600)! / 1.1 - 1)).toBeLessThan(0.15);
    });
    it('peak BMEP at the MBT spark ≈ 70 psi (Fig. 4, read by eye) within ±10 % at 800 rpm', () => {
      const b = Math.max(...mbt.get(800)!.map((r) => r.bmepPsi));
      expect(Math.abs(b / 70 - 1)).toBeLessThan(0.1);
    });
  });

  describe('3. Motored compression pressure (grade C)', () => {
    it.fails('55–60 psig at 150 rpm, wide open [model 67.9 psig: no blow-by; CR 3.78 would give 62.8 — README]', () => {
      const g = (motored.s.reduce((a, x) => a + x.peakPressure, 0) / motored.s.length - MODEL_T_FORD_DYNO.ambientPressure) / PSI;
      expect(g).toBeGreaterThan(55);
      expect(g).toBeLessThan(60);
    });
    it('the residual is bounded: 55–70 psig (Ford 60, Motor Age 55, Dyke 64 at its stated CR 3.6)', () => {
      const g = (motored.s.reduce((a, x) => a + x.peakPressure, 0) / motored.s.length - MODEL_T_FORD_DYNO.ambientPressure) / PSI;
      expect(g).toBeGreaterThan(55);
      expect(g).toBeLessThan(70);
    });
  });

  describe('4. Vehicle (grade C)', () => {
    it('top speed 42–45 mph: WOT, lever 64.5, high gear, level road', () => {
      const load = new LoadModel(MODEL_T.vehicle).configure({ ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'high', grade: 0 } });
      const pts: [number, number][] = [1700, 1800].map((rpm) => [rpm, wot.get(rpm)!.find((r) => r.lever === 64.5)!.Tb * LBFT]);
      pts.push([1900, top1900.e.brakeTorque]);
      let eq = NaN;
      for (let k = 1; k < pts.length; k++) {
        const [r0, t0] = pts[k - 1];
        const [r1, t1] = pts[k];
        const f0 = t0 - load.torque((r0 * Math.PI) / 30);
        const f1 = t1 - load.torque((r1 * Math.PI) / 30);
        if (f0 >= 0 && f1 < 0) eq = r0 + ((r1 - r0) * f0) / (f0 - f1);
      }
      const mph = load.vehicleSpeed((eq * Math.PI) / 30) / MPH;
      console.log(`[Model T] top speed ${mph.toFixed(1)} mph at ${eq.toFixed(0)} rpm`);
      expect(mph).toBeGreaterThan(42);
      expect(mph).toBeLessThan(45);
    });
    it('cruise preset (throttle 0.35, lever 50) settles at 25–35 mph on a level road, without knock', () => {
      const load = new LoadModel(MODEL_T.vehicle).configure({ ...MODEL_T_CRUISE, load: { kind: 'vehicle', gear: 'high', grade: 0 } });
      const [a, b] = cruise;
      const fa = a.r.e.brakeTorque - load.torque((a.rpm * Math.PI) / 30);
      const fb = b.r.e.brakeTorque - load.torque((b.rpm * Math.PI) / 30);
      expect(fa).toBeGreaterThan(0);
      expect(fb).toBeLessThan(0);
      const rpm = a.rpm + ((b.rpm - a.rpm) * fa) / (fa - fb);
      const mph = load.vehicleSpeed((rpm * Math.PI) / 30) / MPH;
      console.log(`[Model T] cruise ${mph.toFixed(1)} mph at ${rpm.toFixed(0)} rpm`);
      expect(mph).toBeGreaterThan(25);
      expect(mph).toBeLessThan(35);
      for (const c of cruise) expect(knocks(c.r.s)).toBe(0);
    });
  });

  describe('5. Knock (PRF 45 surrogate of 1918–22 gasoline, ON ≈ 40–55)', () => {
    it('the period-fuel surrogate is inside the evidence range', () => {
      expect(MODEL_T_FORD_DYNO.fuel.kind).toBe('PRF');
      const on = (MODEL_T_FORD_DYNO.fuel as { octaneNumber: number }).octaneNumber;
      expect(on).toBeGreaterThanOrEqual(40);
      expect(on).toBeLessThanOrEqual(55);
    });
    it('spark knock at WOT and 500–900 rpm with the lever fully advanced; none at 1500–1800 rpm', () => {
      for (const rpm of [500, 700, 900]) expect(wot.get(rpm)!.find((r) => r.lever === 64.5)!.knock).toBeGreaterThan(0);
      for (const rpm of [1500, 1600, 1700, 1800]) expect(wot.get(rpm)!.find((r) => r.lever === 64.5)!.knock).toBe(0);
    });
  });

  describe('6. Idle preset', () => {
    it('settles at 350–460 rpm declutched (Ford: "about 400 R.P.M."), no misfire', () => {
      console.log(`[Model T] idle rpm: ${idleRpm.map((r) => r.toFixed(0)).join(' ')}`);
      const last = idleRpm[idleRpm.length - 1];
      expect(last).toBeGreaterThan(350);
      expect(last).toBeLessThan(460);
      expect(Math.abs(last - idleRpm[idleRpm.length - 3])).toBeLessThan(10);
      expect(idleMisfires).toBe(0);
    });
  });

  it('displacement used for the BMEP figures', () => {
    expect(VD_TOTAL).toBeCloseTo(2.8958e-3, 6);
  });
});
