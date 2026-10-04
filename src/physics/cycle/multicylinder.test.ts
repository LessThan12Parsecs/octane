/**
 * Multi-cylinder cycle model (Model T integration). The L-head chamber, the trembler ignition and the
 * side-valve cams are separate work, so the engine here is the CFR F-1 geometry (flat disc, inductive
 * coil) replicated N times on one crank with shared intake and exhaust plenums; the intake/exhaust
 * network (plenum volumes, venturi and outlet areas) is scaled by N so that each cylinder sees the
 * single-cylinder network.
 *  - Oracle (analytic scaling limit): N cylinders firing IN PHASE on an N-times larger network are N
 *    copies of the single-cylinder engine — every per-cylinder summary must equal the CFR's.
 *  - Inline four, firing order 1-2-4-3 (offsets 0/180/540/360°): identical cylinders give identical
 *    per-cylinder cycles at the same LOCAL angles, i.e. shifted by their firing offsets.
 *  - Conservation of mass, species/elements and energy of the whole network, determinism, the snapshot
 *    contract (EngineSnapshot.cylinders) and the CPU cost per engine cycle.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSpec } from '../core/engine-spec';
import { cylinderAngleDeg } from '../core/engine-spec';
import type { CycleSummary, EngineSnapshot } from '../core/snapshot';
import { NE, NS } from '../core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { elementMoles } from '../thermo/mixture';
import { MOLAR_MASS } from '../thermo/thermo';
import { CycleModel, I_HO, I_HV, I_LO, I_LV } from './cycle-model';
import { EngineSimulator } from './index';

/** CFR F-1 geometry with n = offsets.length cylinders and the network scaled by n (an ad-hoc spec: CFR option defaults). */
function cfrMulti(offsets: number[]): EngineSpec {
  const n = offsets.length;
  const m = CFR_F1.manifolds;
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
      ...m,
      intakeVolume: m.intakeVolume * n,
      exhaustVolume: m.exhaustVolume * n,
      throttleDiameter: m.throttleDiameter * Math.sqrt(n),
      exhaustOutletDiameter: m.exhaustOutletDiameter * Math.sqrt(n),
    },
  };
}

const FOUR = cfrMulti([0, 180, 540, 360]);

/** System inventory + boundary ledgers of the whole network (all cylinders). */
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
const byCylinder = (s: CycleSummary[], i: number): CycleSummary[] => s.filter((x) => (x.cylinder ?? 0) === i);

describe('analytic scaling limit: N in-phase cylinders on an N× network = N single cylinders', () => {
  const ref = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, {}).runCycles(2);
  const m = new CycleModel(cfrMulti([0, 0, 0, 0]), CFR_RON_CONDITIONS, {});
  const sums = m.runCycles(2);

  it('reproduces every per-cylinder summary of the single-cylinder CFR', () => {
    expect(sums.length).toBe(8);
    let worst = 0;
    for (const s of sums) {
      const r = ref[s.cycle];
      expect(s.cylinder).toBeGreaterThanOrEqual(0);
      for (const k of ['imepNet', 'imepGross', 'peakPressure', 'trappedMass', 'volumetricEfficiency', 'noPpm', 'heatLoss'] as const) {
        expect(rel(s[k], r[k]), k).toBeLessThan(1e-9);
        worst = Math.max(worst, rel(s[k], r[k]));
      }
      for (const k of ['ca10', 'ca50', 'ca90', 'peakPressureDeg', 'knockOnsetDeg'] as const) expect(Math.abs(s[k] - r[k]), k).toBeLessThan(1e-6);
      expect(rel(s.mapo, r.mapo)).toBeLessThan(1e-7);
      expect(s.misfire).toBe(r.misfire);
    }
    console.log(`[multi-cylinder] in-phase ×4 vs single CFR: worst relative difference ${worst.toExponential(1)} (IMEP, p_max, m_trap, η_v, NO, Q)`);
  });

  it('the engine summary of N copies is N times the single cylinder', () => {
    const e = sums.find((x) => x.cycle === 1 && x.cylinder === 0)!.engine!;
    const r = ref[1];
    expect(rel(e.imepNet, r.imepNet)).toBeLessThan(1e-9);
    expect(rel(e.indicatedTorque, (4 * r.imepNet * m.kin.displacedVolume) / (4 * Math.PI))).toBeLessThan(1e-9);
    expect(rel(e.volumetricEfficiency, r.volumetricEfficiency)).toBeLessThan(1e-9);
  });
});

describe('inline four on the CFR geometry, firing order 1-2-4-3 (offsets 0/180/540/360°)', () => {
  // (the cylinders start the run at different points of their cycles; the inter-cylinder difference of
  // that cold start decays by ≈ 5–8× per cycle: 6 warm-up cycles bring it below 1e-7 by cycle 2)
  const m = new CycleModel(FOUR, CFR_RON_CONDITIONS, { warmupCycles: 6 });
  const sums = m.runCycles(3);

  it('emits one summary per cylinder and cycle, numbered by local cycle', () => {
    expect(sums.length).toBe(12);
    for (let i = 0; i < 4; i++) expect(byCylinder(sums, i).map((s) => s.cycle)).toEqual([0, 1, 2]);
    // cylinder 0's summaries carry the engine summary; the others do not
    for (const s of sums) expect(s.engine !== undefined).toBe(s.cylinder === 0);
  });

  it('identical cylinders give identical cycles at the same local angles (≤ 1e-6 in IMEP)', () => {
    const last = sums.filter((s) => s.cycle === 2);
    const r = last.find((s) => s.cylinder === 0)!;
    let worst = 0;
    for (const s of last) {
      worst = Math.max(worst, rel(s.imepNet, r.imepNet));
      expect(rel(s.imepNet, r.imepNet)).toBeLessThan(1e-6);
      expect(rel(s.peakPressure, r.peakPressure)).toBeLessThan(1e-5);
      expect(rel(s.trappedMass, r.trappedMass)).toBeLessThan(1e-6);
      expect(Math.abs(s.peakPressureDeg - r.peakPressureDeg)).toBeLessThan(0.01);
      expect(Math.abs(s.ca50 - r.ca50)).toBeLessThan(0.01);
      expect(Math.abs(s.knockOnsetDeg - r.knockOnsetDeg)).toBeLessThan(0.01);
      expect(s.misfire).toBe(false);
    }
    console.log(`[multi-cylinder] inline four, cycle 2: IMEP_n ${(r.imepNet / 1e5).toFixed(4)} bar, p_max ${(r.peakPressure / 1e5).toFixed(2)} bar @ ${r.peakPressureDeg.toFixed(2)}° (local), ` +
      `knock ${r.knockOnsetDeg.toFixed(2)}°, η_v ${r.volumetricEfficiency.toFixed(4)}; max per-cylinder IMEP spread ${worst.toExponential(1)}`);
    expect(m.closureFailures).toBe(0);
  });

  it('engine summary: brake = indicated − friction, fixed-speed load = brake, MEPs over the total displacement', () => {
    const s0 = sums.filter((s) => s.cylinder === 0)[2];
    const e = s0.engine!;
    const Vt = 4 * m.kin.displacedVolume;
    expect(e.rpmMean).toBeCloseTo(CFR_RON_CONDITIONS.rpm, 6);
    expect(e.brakeTorque).toBeCloseTo(e.indicatedTorque - e.frictionTorque, 9);
    expect(e.loadTorque).toBe(e.brakeTorque);
    expect(rel(e.bmep, (4 * Math.PI * e.brakeTorque) / Vt)).toBeLessThan(1e-12);
    expect(rel(e.brakePower, (e.brakeTorque * 2 * Math.PI * e.rpmMean) / 60)).toBeLessThan(1e-9);
    // the per-cylinder IMEPs (local cycles) and the engine's (engine cycle) agree at the periodic state
    const mean = sums.filter((s) => s.cycle === 2).reduce((a, s) => a + s.imepNet, 0) / 4;
    expect(rel(e.imepNet, mean)).toBeLessThan(1e-5);
    // cylinder 0's local cycle IS the engine cycle: same venturi window
    expect(e.volumetricEfficiency).toBe(s0.volumetricEfficiency);
    expect(e.fuelMassFlow).toBeGreaterThan(0);
    expect(e.brakeEfficiency).toBeGreaterThan(0.1);
    expect(e.brakeEfficiency).toBeLessThan(0.4);
    console.log(`[multi-cylinder] engine: T_i ${e.indicatedTorque.toFixed(2)} N m, T_f ${e.frictionTorque.toFixed(2)} N m, T_b ${e.brakeTorque.toFixed(2)} N m, ` +
      `P_b ${(e.brakePower / 1e3).toFixed(2)} kW, BMEP ${(e.bmep / 1e5).toFixed(3)} bar, η_b ${(100 * e.brakeEfficiency).toFixed(1)} %`);
  });

  it('cycle-mean friction torque = the cycle average of the instantaneous per-cylinder friction (quadrature)', () => {
    const om = (CFR_RON_CONDITIONS.rpm * 2 * Math.PI) / 60;
    const n = 7200;
    let s = 0;
    for (let j = 0; j < n; j++) {
      const th = -360 + (720 * (j + 0.5)) / n;
      for (let i = 0; i < 4; i++) s -= m.friction.torque(m.kin.dxdTheta(cylinderAngleDeg(FOUR, i, th) * (Math.PI / 180)), om);
    }
    expect(rel(m.meanFrictionTorque(om), s / n)).toBeLessThan(1e-6);
  });
});

describe('conservation of the four-cylinder network', () => {
  it('motored: mass, species and energy ledgers ≤ 1e-12', () => {
    const m = new CycleModel(FOUR, CFR_RON_CONDITIONS, { combustionModel: 'none', warmupCycles: 0, wallTemperatureModel: 'fixed' });
    const l0 = ledgers(m);
    m.runCycles(3);
    const l1 = ledgers(m);
    expect(rel(l1.mass, l0.mass)).toBeLessThan(1e-12);
    let nTot = 0;
    for (let k = 0; k < NS; k++) nTot += Math.abs(l0.species[k]);
    for (let k = 0; k < NS; k++) expect(Math.abs(l1.species[k] - l0.species[k])).toBeLessThan(1e-12 * nTot);
    // (energy relative to the network's absolute internal energy; the single-cylinder test allows 1e-7 J)
    expect(Math.abs(l1.energy - l0.energy)).toBeLessThan(1e-7);
  });

  it('fired: mass and elements ≤ 1e-12, energy ≤ 1e-10 of the fuel energy', () => {
    const m = new CycleModel(FOUR, CFR_RON_CONDITIONS, { warmupCycles: 1 });
    const l0 = ledgers(m);
    const s = m.runCycles(2);
    const l1 = ledgers(m);
    expect(rel(l1.mass, l0.mass)).toBeLessThan(1e-12);
    for (let e = 0; e < NE; e++) if (l0.elements[e] !== 0) expect(rel(l1.elements[e], l0.elements[e])).toBeLessThan(1e-12);
    const fuelEnergy = s.reduce((a, x) => a + x.fuelMass, 0) * 44e6;
    expect(fuelEnergy).toBeGreaterThan(0);
    expect(Math.abs(l1.energy - l0.energy) / fuelEnergy).toBeLessThan(1e-10);
  });
});

describe('four-cylinder EngineSimulator', () => {
  const opts = { snapshotEveryDeg: 2, bufferAheadSeconds: 0.1, warmupCycles: 2 };
  const take = (sim: EngineSimulator, n: number): EngineSnapshot[] => {
    const out: EngineSnapshot[] = [];
    for (let i = 0; i < n; i++) out.push(sim.advanceToNextSnapshot());
    return out;
  };
  // one whole engine cycle (every cylinder's spark, dense spark and knock sampling included)
  const a = new EngineSimulator(FOUR, CFR_RON_CONDITIONS, opts);
  const sa: EngineSnapshot[] = [a.advanceToNextSnapshot()];
  while (sa[sa.length - 1].cycle < 1) sa.push(a.advanceToNextSnapshot());

  it('is deterministic: two instances and a reset() reproduce the snapshot stream bit for bit', () => {
    const b = new EngineSimulator(FOUR, CFR_RON_CONDITIONS, opts);
    expect(take(b, sa.length)).toEqual(sa);
    a.setOperatingPoint({ sparkAdvanceDeg: 20 });
    take(a, 50);
    a.setOperatingPoint({ sparkAdvanceDeg: CFR_RON_CONDITIONS.sparkAdvanceDeg });
    a.reset();
    expect(a.time).toBe(0);
    expect(take(a, sa.length)).toEqual(sa);
  });

  it('snapshots carry every cylinder at its local angle; top level = cylinder 0, gas torque = sum', () => {
    for (const s of sa) {
      const cs = s.cylinders!;
      expect(cs.length).toBe(4);
      for (let i = 0; i < 4; i++) {
        expect(cs[i].index).toBe(i);
        const exp = cylinderAngleDeg(FOUR, i, s.thetaDeg);
        expect(Math.abs(cs[i].thetaDeg - exp) < 1e-9 || Math.abs(Math.abs(cs[i].thetaDeg - exp) - 720) < 1e-9).toBe(true);
      }
      expect(cs[0].thetaDeg).toBe(s.thetaDeg);
      expect(cs[0].pressure).toBe(s.pressure);
      expect(cs[0].volume).toBe(s.volume);
      expect(cs[0].flame).toEqual(s.flame);
      expect(cs[0].spark).toEqual(s.spark);
      expect(s.gasTorque).toBeCloseTo(cs.reduce((x, c) => x + c.gasTorque, 0), 9);
      expect(s.frictionTorque!).toBeGreaterThan(0);
      expect(s.loadTorque).toBe(s.netTorque); // fixed speed: the dynamometer absorbs the net torque
    }
    // the cylinders are phase-shifted copies: cylinder i at engine angle θ is cylinder 0 at its local
    // angle θ − φ_i (periodic state after the warm-up; grid snapshots of the engine cycle on the common 2° grid)
    const at = new Map<number, EngineSnapshot>();
    for (const s of sa) if (s.cycle === 0 && Number.isInteger(s.thetaDeg) && s.thetaDeg % 2 === 0) at.set(s.thetaDeg, s);
    let n = 0;
    let worst = 0;
    for (const s of at.values()) {
      for (let i = 1; i < 4; i++) {
        const c = s.cylinders![i];
        const r = at.get(c.thetaDeg)?.cylinders![0];
        if (!r) continue;
        expect(rel(c.volume, r.volume)).toBeLessThan(1e-12);
        expect(c.phase).toBe(r.phase);
        worst = Math.max(worst, rel(c.pressure, r.pressure));
        n++;
      }
    }
    expect(n).toBeGreaterThan(900);
    expect(worst).toBeLessThan(2e-3);
    console.log(`[multi-cylinder] ${sa.length} snapshots per engine cycle at 2°; phase-shifted pressures of cylinders 2–4 vs 1: max relative difference ${worst.toExponential(1)}`);
  });

  it('samples every cylinder densely around its spark (5 µs)', () => {
    for (let i = 0; i < 4; i++) {
      const k = sa.findIndex((s) => s.cylinders![i].spark.phase === 'breakdown');
      expect(k, `cylinder ${i}`).toBeGreaterThan(0);
      for (let j = k; j < k + 8; j++) expect(sa[j + 1].t - sa[j].t).toBeLessThan(5.01e-6);
    }
  });
});

describe('four-cylinder speed modes and cost', () => {
  it("free speed (temporary rigid-crank sum): accelerates under a light constant load; other loads are refused", () => {
    const m = new CycleModel(FOUR, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: 10 }, { warmupCycles: 1 });
    const rpm0 = m.rpm;
    const s = m.runCycles(1);
    expect(s.every((x) => Number.isFinite(x.imepNet))).toBe(true);
    expect(m.rpm).toBeGreaterThan(rpm0);
    expect(m.rpm).toBeLessThan(rpm0 * 1.5);
    expect(() => m.setOperatingPoint({ load: { kind: 'brake', refRpm: 600, exponent: 2 } })).toThrow(/multi-cylinder/);
  });

  it('free speed: the engine summary closes the crank energy balance (brake = load torque) for 1 and 4 cylinders', () => {
    // brake torque = indicated − friction − ΔE_kin/4π over the cycle must equal the constant load; the
    // residual is the friction-mean approximation at the cycle-mean speed of an accelerating crank
    // (also the regression case of a fine phase across θ = 360 in free mode)
    for (const [spec, load] of [[CFR_F1, 3], [FOUR, 40]] as const) {
      const m = new CycleModel(spec, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: load }, { warmupCycles: 1 });
      const es = m.runCycles(3).filter((s) => s.engine).map((s) => s.engine!);
      expect(es.length).toBe(3);
      for (const e of es) {
        expect(e.loadTorque).toBe(load);
        expect(Math.abs(e.brakeTorque - load)).toBeLessThan(2e-3 * e.frictionTorque);
      }
      expect(es[2].rpmMean).toBeGreaterThan(es[0].rpmMean); // light load: accelerating
    }
  });

  it('free speed: the four-cylinder snapshot stream runs through the engine wraps', () => {
    const sim = new EngineSimulator(FOUR, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: 40 }, { snapshotEveryDeg: 5, warmupCycles: 1 });
    let s = sim.advanceToNextSnapshot();
    let n = 0;
    while (s.cycle < 2 && n++ < 20_000) {
      const t = s.t;
      s = sim.advanceToNextSnapshot();
      expect(s.t).toBeGreaterThan(t);
    }
    expect(s.cycle).toBe(2);
    expect(s.loadTorque).toBe(40);
  });

  it('CPU per engine cycle at 600 and 1600 rpm (logged)', () => {
    const out: string[] = [];
    for (const rpm of [600, 1600]) {
      const m = new CycleModel(FOUR, { ...CFR_RON_CONDITIONS, rpm }, { warmupCycles: 2 });
      m.runCycles(1);
      const t0 = performance.now();
      const s = m.runCycles(2);
      const ms = (performance.now() - t0) / 2;
      out.push(`${rpm} rpm: ${ms.toFixed(0)} ms per engine cycle (real time ${(120000 / rpm).toFixed(0)} ms), IMEP_n ${(s[s.length - 1].imepNet / 1e5).toFixed(2)} bar`);
      expect(ms).toBeLessThan(5000); // generous for loaded CI machines (vitest transform ≈ 2× slower than node)
    }
    console.log(`[multi-cylinder] four cylinders: ${out.join('; ')}`);
  });
});
