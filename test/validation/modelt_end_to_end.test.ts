/**
 * Ford Model T end to end (Model T integration): the real MODEL_T spec — L-head chamber with its block-deck
 * surface, three-arc cams and side-valve flow, venturi + butterfly carburettor, four trembler coils on the
 * flywheel magneto, the rigid four-cylinder crank train with whole-engine friction — run through
 * CycleModel and EngineSimulator. Written for the uncalibrated model; the checks here are integration and
 * physics properties, not Ford's numbers (those are test/validation/measured_modelt_data.test.ts, against the
 * calibrated MODEL_T_CALIBRATION): finite, deterministic output; four identical cylinders giving equal cycles 180° apart in the
 * firing order 1-2-4-3; no misfire; the spark train of each timer contact; conservation of mass, species /
 * elements and energy of the whole four-cylinder network (the tolerances of the CFR and multi-cylinder
 * tests); the snapshot contract (cylinders[] at local angles, trembler spark fields, magneto EMF, firing
 * cylinder, torques) and the engine summary. Runtime ≈ 30–60 s (vitest).
 */
import { describe, expect, it } from 'vitest';
import { cylinderAngleDeg } from '../../src/physics/core/engine-spec';
import type { CycleSummary, EngineSnapshot } from '../../src/physics/core/snapshot';
import { NE, NS } from '../../src/physics/core/species';
import { MODEL_T, MODEL_T_FORD_DYNO } from '../../src/physics/engines/model-t';
import { elementMoles } from '../../src/physics/thermo/mixture';
import { MOLAR_MASS } from '../../src/physics/thermo/thermo';
import { CycleModel, I_HO, I_HV, I_LO, I_LV } from '../../src/physics/cycle/cycle-model';
import { EngineSimulator } from '../../src/physics/cycle/index';
import { DENSE_SPARK_DT } from '../../src/physics/cycle/engine-simulator';

/** System inventory + boundary ledgers of the whole network (all cylinders, both plenums). */
function ledgers(m: CycleModel) {
  const y = m.y;
  const LV = new Float64Array(NS);
  const LO = new Float64Array(NS);
  let mv = 0;
  let mo = 0;
  for (let k = 0; k < NS; k++) {
    LV[k] = y[I_LV + k];
    LO[k] = y[I_LO + k];
    mv += LV[k] * MOLAR_MASS[k];
    mo += LO[k] * MOLAR_MASS[k];
  }
  const sp = m.speciesInventory();
  const species = new Float64Array(NS);
  for (let k = 0; k < NS; k++) species[k] = sp[k] - LV[k] + LO[k];
  const b = m.elementInventory();
  const bV = elementMoles(LV);
  const bO = elementMoles(LO);
  const elements = new Float64Array(NE);
  for (let e = 0; e < NE; e++) elements[e] = b[e] - bV[e] + bO[e];
  return {
    mass: m.systemMass() - mv + mo,
    energy: m.systemEnergy() - y[I_HV] + y[I_HO] + m.heatLedger() + m.workLedger() - m.sparkEnergy,
    species,
    elements,
  };
}

const rel = (a: number, b: number): number => Math.abs(a / b - 1);
const LBFT = 1.3558179483314004;

/** Every number in a snapshot / summary tree is finite, except the documented NaN fields. */
function nonFinite(o: unknown, path = ''): string[] {
  const out: string[] = [];
  const walk = (v: unknown, p: string): void => {
    if (typeof v === 'number') {
      if (!Number.isFinite(v) && !/(firstSparkDeg|sparkDeg|ca10|ca50|ca90|knockOnsetDeg|isfc)$/.test(p)) out.push(`${p}=${v}`);
    } else if (Array.isArray(v)) v.forEach((x, i) => walk(x, `${p}[${i}]`));
    else if (v && typeof v === 'object') for (const [k, x] of Object.entries(v)) walk(x, `${p}.${k}`);
  };
  walk(o, path);
  return out;
}

describe("Model T at Ford's dyno condition (fixed 1600 rpm, wide open, lever fully advanced, PRF 45, φ 1.15)", () => {
  const m = new CycleModel(MODEL_T, MODEL_T_FORD_DYNO, { warmupCycles: 3 });
  const sums = m.runCycles(2);
  const last = sums.filter((s) => s.cycle === 1);
  const e = last.find((s) => s.cylinder === 0)!.engine!;

  it('runs: one finite summary per cylinder and cycle, no misfire, no failed closure solve', () => {
    expect(sums.length).toBe(8);
    for (let i = 0; i < 4; i++) expect(sums.filter((s) => s.cylinder === i).map((s) => s.cycle)).toEqual([0, 1]);
    for (const s of sums) {
      expect(nonFinite(s, `cyl${s.cylinder}`)).toEqual([]);
      expect(s.misfire).toBe(false);
      expect(Number.isFinite(s.ca50)).toBe(true);
    }
    expect(m.closureFailures).toBe(0);
  });

  it('four identical cylinders give equal cycles at their local angles (IMEP ≤ 2e-3, angles ≤ 0.5°)', () => {
    const r = last.find((s) => s.cylinder === 0)!;
    let worst = 0;
    for (const s of last) {
      worst = Math.max(worst, rel(s.imepNet, r.imepNet));
      expect(rel(s.imepNet, r.imepNet)).toBeLessThan(2e-3);
      expect(rel(s.peakPressure, r.peakPressure)).toBeLessThan(5e-3);
      expect(rel(s.trappedMass, r.trappedMass)).toBeLessThan(1e-3);
      expect(Math.abs(s.ca50 - r.ca50)).toBeLessThan(0.5);
      expect(Math.abs(s.sparkDeg! - r.sparkDeg!)).toBeLessThan(0.5);
    }
    console.log(`[Model T] 1600 rpm WOT: per-cylinder IMEP spread ${worst.toExponential(1)} after 3 warm-up cycles`);
  });

  it('the spark train of each timer contact: first breakdown after the make (lever 64.5°), several breakdowns', () => {
    for (const s of sums) {
      expect(s.sparkDeg!).toBeGreaterThan(-64.5);
      expect(s.sparkDeg!).toBeLessThan(0);
      expect(s.sparkCount!).toBeGreaterThan(1);
    }
  });

  it('engine summary: brake = indicated − friction at fixed speed, load = brake, plausible (calibrated) levels', () => {
    expect(e.rpmMean).toBeCloseTo(1600, 6);
    expect(e.brakeTorque).toBeCloseTo(e.indicatedTorque - e.frictionTorque, 9);
    expect(e.loadTorque).toBe(e.brakeTorque);
    const Vt = 4 * m.kin.displacedVolume;
    expect(rel(e.bmep, (4 * Math.PI * e.brakeTorque) / Vt)).toBeLessThan(1e-12);
    // plausibility bands (Ford: 65 lb-ft = 88 N m at 1600 rpm). Calibration phase: the friction is PNH with the
    // calibrated ring-pack and transmission-churning factors (options.ts MODEL_T_CALIBRATED_FRICTION, ≈ 1.13 bar
    // at 1600 rpm; pure PNH ≈ 0.66 bar) — the band was 0.5–0.8 bar for the uncalibrated model
    expect(e.brakeTorque).toBeGreaterThan(0.5 * 65 * LBFT);
    expect(e.brakeTorque).toBeLessThan(2 * 65 * LBFT);
    expect(e.fmep / 1e5).toBeGreaterThan(0.9);
    expect(e.fmep / 1e5).toBeLessThan(1.4);
    expect(e.volumetricEfficiency).toBeGreaterThan(0.4);
    expect(e.volumetricEfficiency).toBeLessThan(0.95);
    expect(e.brakeEfficiency).toBeGreaterThan(0.05);
    expect(e.brakeEfficiency).toBeLessThan(0.3);
    for (const s of last) {
      expect(s.peakPressure / 1e5).toBeGreaterThan(8);
      expect(s.peakPressure / 1e5).toBeLessThan(40);
      expect(s.residualFraction).toBeGreaterThan(0.02);
      expect(s.residualFraction).toBeLessThan(0.3);
    }
    console.log(
      `[Model T] 1600 rpm WOT: T_b ${e.brakeTorque.toFixed(1)} N m (${(e.brakeTorque / LBFT).toFixed(1)} lb-ft; Ford 65), ` +
        `IMEP ${(e.imepNet / 1e5).toFixed(2)} bar, FMEP ${(e.fmep / 1e5).toFixed(2)} bar, η_v ${e.volumetricEfficiency.toFixed(3)}, ` +
        `CA50 ${last[0].ca50.toFixed(1)}°, first spark ${last[0].sparkDeg!.toFixed(1)}° (${last[0].sparkCount} breakdowns)`,
    );
  });
});

describe('Model T: determinism', () => {
  it('two engines from the same operating point give bit-identical summaries and states', () => {
    const op = { ...MODEL_T_FORD_DYNO, rpm: 1200 };
    const a = new CycleModel(MODEL_T, op, { warmupCycles: 0 });
    const b = new CycleModel(MODEL_T, op, { warmupCycles: 0 });
    const sa = a.runCycles(1);
    const sb = b.runCycles(1);
    expect(sb).toEqual(sa);
    expect(Array.from(b.y)).toEqual(Array.from(a.y));
  });
});

describe('Model T: conservation of the four-cylinder network (L-head block ledger, carburettor, trembler)', () => {
  it('motored: mass, species ≤ 1e-12, energy drift ≤ 1e-7 J', () => {
    const m = new CycleModel(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm: 900 }, { combustionModel: 'none', warmupCycles: 0 });
    const l0 = ledgers(m);
    m.runCycles(2);
    const l1 = ledgers(m);
    expect(rel(l1.mass, l0.mass)).toBeLessThan(1e-12);
    let nTot = 0;
    for (let k = 0; k < NS; k++) nTot += Math.abs(l0.species[k]);
    for (let k = 0; k < NS; k++) expect(Math.abs(l1.species[k] - l0.species[k])).toBeLessThan(1e-12 * nTot);
    expect(Math.abs(l1.energy - l0.energy)).toBeLessThan(1e-7);
  });

  it('fired: mass and elements ≤ 1e-12, energy ≤ 1e-10 of the fuel energy', () => {
    const m = new CycleModel(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm: 1200 }, { warmupCycles: 1 });
    const l0 = ledgers(m);
    const s = m.runCycles(2);
    const l1 = ledgers(m);
    expect(rel(l1.mass, l0.mass)).toBeLessThan(1e-12);
    for (let e = 0; e < NE; e++) if (l0.elements[e] !== 0) expect(rel(l1.elements[e], l0.elements[e])).toBeLessThan(1e-12);
    const fuelEnergy = s.reduce((a, x) => a + x.fuelMass, 0) * 44e6;
    expect(fuelEnergy).toBeGreaterThan(0);
    expect(Math.abs(l1.energy - l0.energy) / fuelEnergy).toBeLessThan(1e-10);
    // the block-deck surface takes part of the wall heat (its own ledger, appended after the cylinder blocks)
    expect(m.ny).toBe(159 + 3 * 42 + 4);
  });
});

describe('Model T EngineSimulator: the snapshot contract', () => {
  const sim = new EngineSimulator(MODEL_T, MODEL_T_FORD_DYNO, { snapshotEveryDeg: 2, warmupCycles: 2 });
  const snaps: EngineSnapshot[] = [sim.advanceToNextSnapshot()];
  while (snaps[snaps.length - 1].cycle < 1) snaps.push(sim.advanceToNextSnapshot());
  const drained: CycleSummary[] = sim.drainCycleSummaries();
  const offsets = MODEL_T.layout.firingOffsetDeg;

  it('cylinders[] at their local angles, top level = cylinder 1, torques and the engine fields present and finite', () => {
    expect(snaps.length).toBeGreaterThan(360);
    for (const s of snaps) {
      expect(nonFinite(s, 'snap')).toEqual([]);
      const cs = s.cylinders!;
      expect(cs.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(cs[i].index).toBe(i);
        const exp = cylinderAngleDeg(MODEL_T, i, s.thetaDeg);
        expect(Math.abs(cs[i].thetaDeg - exp) < 1e-9 || Math.abs(Math.abs(cs[i].thetaDeg - exp) - 720) < 1e-9).toBe(true);
      }
      expect(cs[0].pressure).toBe(s.pressure);
      expect(cs[0].spark).toEqual(s.spark);
      expect(s.gasTorque).toBeCloseTo(cs.reduce((x, c) => x + c.gasTorque, 0), 9);
      expect(s.frictionTorque!).toBeGreaterThan(0);
      expect(s.loadTorque).toBe(s.netTorque); // fixed speed: the dynamometer absorbs the net torque
      expect(Number.isFinite(s.magnetoEmf!)).toBe(true);
      // the firing cylinder is the one whose timer contact is closed
      const closed = cs.filter((c) => c.spark.timerClosed).map((c) => c.index);
      expect(closed.length).toBeLessThanOrEqual(1);
      expect(s.firingCylinder).toBe(closed.length ? closed[0] : -1);
    }
  });

  it('trembler spark fields: one event per timer contact with a train of breakdowns', () => {
    for (let i = 0; i < 4; i++) {
      const sp = snaps.map((s) => s.cylinders![i].spark);
      for (const x of sp) {
        expect(typeof x.breakdownCount).toBe('number');
        expect(typeof x.pointsOpen).toBe('boolean');
        expect(typeof x.timerClosed).toBe('boolean');
        expect(typeof x.firstSparkDeg).toBe('number');
        expect(typeof x.primaryVoltage).toBe('number');
      }
      expect(Math.max(...sp.map((x) => x.breakdownCount!))).toBeGreaterThan(1);
      expect(sp.some((x) => x.pointsOpen)).toBe(true);
    }
  });

  it('the four cylinders fire 180° apart in the firing order 1-2-4-3 (engine angles of the first breakdowns)', () => {
    // engine angle of the first snapshot after each cylinder's first breakdown in engine cycle 0 (the trip
    // interrupt and the 5 µs train sampling put it within ≈ 0.05° of the breakdown)
    const first: number[] = [];
    const local: number[] = [];
    for (let i = 0; i < 4; i++) {
      // (the snapshot where the event's breakdown count leaves 0 — a train already running at t = 0 is skipped)
      const j = snaps.findIndex((x, k) => k > 0 && x.cycle === 0 && x.cylinders![i].spark.breakdownCount! > 0 && snaps[k - 1].cylinders![i].spark.breakdownCount === 0);
      expect(j).toBeGreaterThan(0);
      const s = snaps[j];
      first.push(s.thetaDeg);
      local.push(s.cylinders![i].spark.firstSparkDeg!);
      // the engine angle is the local angle shifted by the firing offset
      expect(Math.abs(cylinderAngleDeg(MODEL_T, i, s.thetaDeg) - s.cylinders![i].thetaDeg) % 720).toBeLessThan(1e-9);
    }
    const order = [0, 1, 3, 2];
    for (let k = 1; k < 4; k++) {
      const d = (((first[order[k]] - first[order[k - 1]]) % 720) + 720) % 720;
      expect(Math.abs(d - 180)).toBeLessThan(0.5);
      expect(Math.abs(local[order[k]] - local[order[0]])).toBeLessThan(0.5);
    }
    expect(offsets).toEqual([0, 180, 540, 360]);
  });

  it('dense 5 µs sampling after a vibrator trip; summaries carry the train and the engine summary', () => {
    const k = snaps.findIndex((s, j) => j > 0 && s.cylinders![0].spark.pointsOpen && !snaps[j - 1].cylinders![0].spark.pointsOpen);
    expect(k).toBeGreaterThan(0);
    for (let j = k; j < k + 6; j++) expect(snaps[j + 1].t - snaps[j].t).toBeLessThan(DENSE_SPARK_DT * 1.001);
    expect(drained.filter((s) => s.engine).length).toBeGreaterThanOrEqual(1);
    for (const s of drained) {
      expect(s.cylinder).toBeGreaterThanOrEqual(0);
      expect(Number.isFinite(s.sparkDeg!)).toBe(true);
      expect(s.sparkCount!).toBeGreaterThan(0);
    }
    console.log(`[Model T] EngineSimulator 1600 rpm, 2° grid: ${snaps.filter((s) => s.cycle === 0).length} snapshots in engine cycle 0`);
  });
});
