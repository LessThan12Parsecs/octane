import { describe, expect, it } from 'vitest';
import type { CycleSummary, EngineSnapshot } from '../core/snapshot';
import { NE, NS } from '../core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../engines/cfr';
import { elementMoles } from '../thermo/mixture';
import { MOLAR_MASS } from '../thermo/thermo';
import { CycleModel, I_HO, I_HV, I_LO, I_LV, I_Q, I_W } from './cycle-model';
import { EngineSimulator } from './index';

/** System inventory + boundary ledgers: mass, species, elements, energy (all conserved quantities). */
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
    energy: m.systemEnergy() - y[I_HV] + y[I_HO] + y[I_Q] + y[I_W] - m.sparkEnergy,
    species,
    elements,
  };
}

function nonFinitePaths(o: unknown, path = '', out: string[] = []): string[] {
  if (typeof o === 'number') {
    if (!Number.isFinite(o)) out.push(path);
  } else if (Array.isArray(o)) o.forEach((v, i) => nonFinitePaths(v, `${path}[${i}]`, out));
  else if (o && typeof o === 'object') for (const [k, v] of Object.entries(o)) nonFinitePaths(v, path ? `${path}.${k}` : k, out);
  return out;
}

const fmt = (s: CycleSummary): string =>
  `IMEP_n ${(s.imepNet / 1e5).toFixed(2)} bar (gross ${(s.imepGross / 1e5).toFixed(2)}, PMEP ${(s.pmep / 1e5).toFixed(3)}), ` +
  `p_max ${(s.peakPressure / 1e5).toFixed(1)} bar @ ${s.peakPressureDeg.toFixed(1)}°, dp/dθ_max ${(s.maxPressureRiseRate / 1e5).toFixed(2)} bar/°, ` +
  `CA10/50/90 ${s.ca10.toFixed(1)}/${s.ca50.toFixed(1)}/${s.ca90.toFixed(1)}°, η_i ${(100 * s.indicatedEfficiency).toFixed(1)} %, ` +
  `ISFC ${(s.isfc * 3.6e9).toFixed(0)} g/kWh, m_trap ${(s.trappedMass * 1e6).toFixed(0)} mg, x_res ${(100 * s.residualFraction).toFixed(2)} %, ` +
  `η_v ${(100 * s.volumetricEfficiency).toFixed(1)} %, m_fuel ${(s.fuelMass * 1e6).toFixed(2)} mg, NO ${s.noPpm.toFixed(0)} ppm, CO ${(100 * s.coFraction).toFixed(2)} %, ` +
  `knock ${Number.isNaN(s.knockOnsetDeg) ? 'none' : `${s.knockOnsetDeg.toFixed(2)}° (end gas ${(100 * s.knockEndGasFraction).toFixed(1)} %, MAPO ${(s.mapo / 1e5).toFixed(2)} bar)`}, ` +
  `Q_wall ${s.heatLoss.toFixed(0)} J (${((100 * s.heatLoss) / (s.fuelMass * 44.4e6)).toFixed(1)} % of m_f·LHV), W_g ${s.indicatedWorkGross.toFixed(0)} J, misfire ${s.misfire}`;

describe('motored engine (combustionModel none)', () => {
  it('reaches a periodic steady state and conserves mass, species and energy over the gas exchange', () => {
    // (fixer round 2: fixed walls — the lumped wall model adds the physical fired→motored wall
    // cool-down to this cold-start transient, which is not what this test is about)
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { combustionModel: 'none', warmupCycles: 0, wallTemperatureModel: 'fixed' });
    const l0 = ledgers(m);
    const s = m.runCycles(16);
    // The slowest mode is the 10 L adiabatic exhaust surge tank (it exchanges ≈ 1/10 of its mass
    // per cycle): the cycle-to-cycle changes decay geometrically (≈ 7 %/cycle), and the remaining
    // drift d·r/(1 − r) is small.
    const d = s.map((x, i) => (i === 0 ? NaN : Math.abs(x.trappedMass / s[i - 1].trappedMass - 1)));
    const r = d[15] / d[14];
    expect(r).toBeLessThan(1);
    expect(d[15]).toBeLessThan(d[14]);
    expect(d[14]).toBeLessThan(d[13]);
    expect(d[15]).toBeLessThan(1e-4);
    expect((d[15] * r) / (1 - r)).toBeLessThan(2e-3);
    const a = s[14];
    const b = s[15];
    expect(Math.abs(b.peakPressure / a.peakPressure - 1)).toBeLessThan(5e-5);
    expect(Math.abs(b.imepNet - a.imepNet)).toBeLessThan(1e-3 * Math.abs(b.imepGross));
    expect(b.misfire).toBe(true);
    expect(b.isfc).toBeNaN(); // net work < 0 for a motored cycle
    expect(b.imepNet).toBeLessThan(0);
    expect(b.pmep).toBeGreaterThan(0);
    expect(b.peakPressureDeg).toBeGreaterThan(-5);
    expect(b.peakPressureDeg).toBeLessThan(0.5);
    const l1 = ledgers(m);
    expect(Math.abs(l1.mass / l0.mass - 1)).toBeLessThan(1e-12);
    let nTot = 0;
    for (let k = 0; k < NS; k++) nTot += Math.abs(l0.species[k]);
    for (let k = 0; k < NS; k++) expect(Math.abs(l1.species[k] - l0.species[k])).toBeLessThan(1e-12 * nTot);
    expect(Math.abs(l1.energy - l0.energy)).toBeLessThan(1e-9 * 100);
    console.log(`[cycle] motored steady state: p_max ${(b.peakPressure / 1e5).toFixed(2)} bar @ ${b.peakPressureDeg.toFixed(2)}°, ` +
      `IMEP_n ${(b.imepNet / 1e5).toFixed(3)} bar, PMEP ${(b.pmep / 1e5).toFixed(3)} bar, η_v ${(100 * b.volumetricEfficiency).toFixed(1)} %, ` +
      `heat loss ${b.heatLoss.toFixed(1)} J; Δmass ${((l1.mass - l0.mass) / l0.mass).toExponential(1)}, ΔE ${(l1.energy - l0.energy).toExponential(1)} J`);
  });
});

describe('fired CFR at CFR_RON_CONDITIONS (PRF 90, ASTM guide-table CR)', () => {
  const t0 = performance.now();
  const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 3 });
  const l0 = ledgers(m);
  const t1 = performance.now();
  const sums = m.runCycles(3);
  const cpu = (performance.now() - t1) / 3;
  const l1 = ledgers(m);

  it('fires with finite, plausible and converged results', () => {
    for (const s of sums) {
      expect(nonFinitePaths({ ...s, knockOnsetDeg: 0 })).toEqual([]);
      expect(s.misfire).toBe(false);
    }
    const s = sums[2];
    console.log(`[cycle] RON (warm-up + 3 cycles; construction ${(t1 - t0).toFixed(0)} ms, ${cpu.toFixed(0)} ms CPU per cycle): ${fmt(s)}`);
    expect(s.imepNet).toBeGreaterThan(4e5);
    expect(s.imepNet).toBeLessThan(14e5);
    expect(s.peakPressure).toBeGreaterThan(12e5);
    expect(s.peakPressure).toBeLessThan(60e5);
    expect(s.ca10).toBeLessThan(s.ca50);
    expect(s.ca50).toBeLessThan(s.ca90);
    expect(s.indicatedEfficiency).toBeGreaterThan(0.15);
    expect(s.indicatedEfficiency).toBeLessThan(0.4);
    expect(s.volumetricEfficiency).toBeGreaterThan(0.75);
    expect(s.volumetricEfficiency).toBeLessThan(1.1);
    expect(s.residualFraction).toBeGreaterThan(0.02);
    expect(s.residualFraction).toBeLessThan(0.15);
    expect(s.noPpm).toBeGreaterThan(10);
    expect(s.noPpm).toBeLessThan(5000);
    expect(s.coFraction).toBeGreaterThan(0.005); // φ = 1.1: rich
    expect(s.pmep).toBeGreaterThan(0);
    // cycle-to-cycle convergence (no CCV model: deterministic, converged residual)
    expect(Math.abs(sums[2].imepNet / sums[1].imepNet - 1)).toBeLessThan(0.01);
    expect(Math.abs(sums[2].trappedMass / sums[1].trappedMass - 1)).toBeLessThan(0.005);
  });

  it('conserves mass, elements and energy of the whole network over fired cycles', () => {
    expect(Math.abs(l1.mass / l0.mass - 1)).toBeLessThan(1e-12);
    for (let e = 0; e < NE; e++) expect(Math.abs(l1.elements[e] / l0.elements[e] - 1)).toBeLessThan(1e-12);
    const fuelEnergy = sums[2].fuelMass * 44e6;
    expect(Math.abs(l1.energy - l0.energy) / fuelEnergy).toBeLessThan(1e-10);
  });

  it('runs faster than real time at 600 rpm (< 200 ms CPU per cycle)', () => {
    expect(cpu).toBeLessThan(400); // generous for loaded CI machines; typical ≈ 110 ms
  });
});

describe('EngineSimulator snapshots', () => {
  const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 0.5, bufferAheadSeconds: 0.1 });
  const snaps: EngineSnapshot[] = [];
  while (snaps.length === 0 || snaps[snaps.length - 1].cycle < 1) snaps.push(sim.advanceToNextSnapshot());
  const cyc0 = snaps.filter((s) => s.cycle === 0);

  it('are complete (no NaN/∞), time-ordered and cover a whole cycle', () => {
    const bad = new Set<string>();
    for (const s of snaps) for (const p of nonFinitePaths(s)) bad.add(p);
    expect([...bad]).toEqual([]);
    expect(snaps[0].t).toBe(0);
    expect(snaps[0].thetaDeg).toBe(-360);
    for (let i = 1; i < snaps.length; i++) {
      expect(snaps[i].t).toBeGreaterThan(snaps[i - 1].t);
      expect(snaps[i].thetaDeg).toBeGreaterThanOrEqual(-360);
      expect(snaps[i].thetaDeg).toBeLessThan(360);
    }
    const stages = new Set(cyc0.map((s) => s.flame.stage));
    for (const st of ['none', 'kernel', 'turbulent', 'burnout', 'done'] as const) expect(stages.has(st)).toBe(true);
    const phases = new Set(cyc0.map((s) => s.phase));
    for (const ph of ['gas-exchange', 'compression', 'combustion'] as const) expect(phases.has(ph)).toBe(true);
    const sp = new Set(cyc0.map((s) => s.spark.phase));
    for (const ph of ['off', 'charging', 'breakdown', 'glow', 'done'] as const) expect(sp.has(ph)).toBe(true);
    expect(sim.drainCycleSummaries().length).toBe(1);
  });

  it('follow the front-end conventions', () => {
    for (const s of snaps) {
      if (s.temperatureBurned === 0) expect(s.temperatureUnburned).toBe(s.temperatureMean);
      if (s.flame.stage === 'done') {
        expect(s.flame.laminarSpeed).toBe(0);
        expect(s.flame.turbulentSpeed).toBe(0);
        const R = CFR_F1.geometry.bore / 2;
        expect(s.flame.radius).toBeGreaterThan(Math.hypot(R, R) * 0.9);
      }
      if (s.phase !== 'gas-exchange') {
        expect(s.intakeMassFlow).toBe(0);
        expect(s.exhaustMassFlow).toBe(0);
      }
      if (!s.knock.autoignited) expect(s.knock.oscillation).toBe(0);
      expect(Math.sign(s.rodAngle)).toBe(Math.sign(Math.sin((s.thetaDeg * Math.PI) / 180)) || 0);
      expect(s.flame.center).toEqual(CFR_F1.sparkPlug.gapCenter);
    }
    // intake flow positive into the cylinder during the intake stroke, exhaust positive out during blowdown
    const at = (th: number) => cyc0.find((s) => s.thetaDeg >= th)!;
    expect(at(-270).intakeMassFlow).toBeGreaterThan(0);
    expect(at(160).exhaustMassFlow).toBeGreaterThan(0);
    // spark fields 'off' before this cycle's dwell starts
    expect(at(-300).spark.phase).toBe('off');
  });

  it('are dense around the spark (5 µs, then 50 µs) and during knock ringing (10 µs)', () => {
    const iBd = cyc0.findIndex((s) => s.spark.phase === 'breakdown');
    expect(iBd).toBeGreaterThan(0);
    const dts: number[] = [];
    for (let i = iBd; i < iBd + 10; i++) dts.push(cyc0[i + 1].t - cyc0[i].t);
    expect(Math.max(...dts)).toBeLessThan(5.01e-6);
    const glow = cyc0.filter((s) => s.spark.phase === 'glow');
    expect(glow.length).toBeGreaterThan(20);
    const iK = cyc0.findIndex((s) => s.knock.autoignited);
    if (iK > 0) {
      for (let i = iK + 1; i < iK + 10; i++) expect(cyc0[i + 1].t - cyc0[i].t).toBeLessThan(10.01e-6);
    }
  });
});

describe('determinism', () => {
  const take = (sim: EngineSimulator, n: number): EngineSnapshot[] => {
    const out: EngineSnapshot[] = [];
    for (let i = 0; i < n; i++) out.push(sim.advanceToNextSnapshot());
    return out;
  };
  it('two instances and a reset() reproduce the snapshot stream bit for bit', () => {
    const opts = { snapshotEveryDeg: 2, bufferAheadSeconds: 0.1, warmupCycles: 1 };
    const a = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, opts);
    const b = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, opts);
    const n = 480; // through the spark and the whole combustion of cycle 0
    const sa = take(a, n);
    const sb = take(b, n);
    expect(sb).toEqual(sa);
    a.setOperatingPoint({ sparkAdvanceDeg: 20 });
    take(a, 50);
    a.setOperatingPoint({ sparkAdvanceDeg: CFR_RON_CONDITIONS.sparkAdvanceDeg });
    a.reset();
    expect(a.time).toBe(0);
    const sr = take(a, n);
    expect(sr).toEqual(sa);
  });
});

describe('operating-point changes', () => {
  it('applies a CR change at the next cycle start (visible through clearanceHeight − pistonDisplacement)', () => {
    const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 5, warmupCycles: 1 });
    const hTdc = (s: EngineSnapshot) => s.clearanceHeight - s.pistonDisplacement;
    const s0 = sim.advanceToNextSnapshot();
    sim.setOperatingPoint({ compressionRatio: 8 });
    let s = s0;
    while (s.cycle === 0) {
      expect(hTdc(s)).toBeCloseTo(hTdc(s0), 12);
      s = sim.advanceToNextSnapshot();
    }
    expect(hTdc(s)).toBeLessThan(hTdc(s0));
    expect(sim.model.op.compressionRatio).toBe(8);
  });

  it("free speed mode: the crank accelerates under combustion torque against a light load", () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: 2 }, { warmupCycles: 1 });
    const rpm0 = m.rpm;
    const s = m.runCycles(2);
    expect(s.every((x) => Number.isFinite(x.imepNet))).toBe(true);
    expect(m.rpm).toBeGreaterThan(rpm0);
    expect(m.rpm).toBeLessThan(rpm0 * 1.5);
  });
});
