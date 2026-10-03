/**
 * Validation round 1 (numerics): step-size convergence, event location, output-sampling and
 * warm-up independence of the cycle results. SLOW (≈ 30–60 s).
 *
 * Criteria (task statement): at default settings the change on refining the steps must be ≤ 0.2 %
 * for IMEP / peak pressure / NO / MAPO / max dp/dθ and ≤ 0.1° for CA10/50/90 and knock onset.
 *
 * Status after the round-1 fixer pass: all pass (Hermite peak / sliding-window dp/dθ / oscillator
 * peak tracking, in-step hand-off, warm-up plenum relaxation).
 * Known failures at the time of writing (validation round 1):
 *  - MAPO and maxPressureRiseRate are not converged / depend on where the steps (or the snapshot
 *    instants) fall: MAPO is the max of |p_osc| sampled at step ends (0.25° = 69 µs after the end-gas
 *    burn-out vs a 150 µs (1,0) period); maxPressureRiseRate is taken over ≥ 0.1° windows anchored at
 *    step ends.
 *  - the kernel → flame hand-off is applied at the end of the fine step in which it happened
 *    (cycle-model.ts ignitionSplit/onHandoff), so CA50 and the knock onset move in 0.05° stairs
 *    with the spark advance and MAPO jumps by ~2 %.
 *  - the default 3 warm-up cycles leave the adiabatic 10 L exhaust plenum ~260 K below its
 *    periodic state (800 K initial → 1273 K), so the first emitted RON cycle is off by 0.5 % IMEP,
 *    3.7 % peak pressure, 9 % NO, 0.55° knock onset, 29 % MAPO.
 */
import { describe, expect, it } from 'vitest';
import type { CycleSummary } from '../../src/physics/core/snapshot';
import { NS } from '../../src/physics/core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { completeCombustionProducts, freshCharge, fuelFromSelection, humidAir } from '../../src/physics/thermo/fuels';
import { mixMolarMass } from '../../src/physics/thermo/mixture';
import { CycleModel, EngineSimulator, MODE_OPEN, type ClosedCycleInit, type CycleModelOptions } from '../../src/physics/cycle/index';
import { I_W } from '../../src/physics/cycle/cycle-model';

function ronIvc(): ClosedCycleInit {
  const op = CFR_RON_CONDITIONS;
  const fuel = fuelFromSelection(op.fuel);
  const air = humidAir(op.ambientTemperature, op.ambientPressure, op.relativeHumidity);
  const fresh = freshCharge({ fuel, phi: op.equivalenceRatio, airX: air });
  const prod = completeCombustionProducts(fresh);
  const yRes = 0.06;
  const X = new Float64Array(NS);
  for (let k = 0; k < NS; k++) X[k] = ((1 - yRes) / mixMolarMass(fresh)) * fresh[k] + (yRes / mixMolarMass(prod)) * prod[k];
  return { startDeg: -152, endDeg: 141, T: 370, p: 1e5, X, residualMassFraction: yRes, uPrime: 5 };
}
const KNOCK_OP = { ...CFR_RON_CONDITIONS, compressionRatio: 9 };

interface Res {
  W: number;
  pmax: number;
  ca10: number;
  ca50: number;
  ca90: number;
  ko: number;
  no: number;
  mapo: number;
}

/** Closed-cycle run from the fixed RON IVC state; optional forced output sampling dtSample after knock onset. */
function closed(opts: Partial<CycleModelOptions>, sparkAdvanceDeg = 13, dtSample = 0): Res {
  const init = ronIvc();
  const m = new CycleModel(CFR_F1, { ...KNOCK_OP, sparkAdvanceDeg }, opts, init);
  let pmax = 0;
  while (!m.finished) {
    if (dtSample > 0 && m.knockOnset && m.t < m.tKnockOnset + 3e-3) m.stepUntil(m.t + dtSample, Infinity);
    else m.stepUntil(Infinity, init.endDeg + 1);
    if (m.p > pmax && m.mode !== MODE_OPEN) pmax = m.p;
  }
  return { W: m.y[I_W], pmax, ca10: m.ca10, ca50: m.ca50, ca90: m.ca90, ko: m.knockOnsetDeg, no: m.nNO, mapo: m.mapo };
}

const rel = (a: number, b: number): number => Math.abs(a / b - 1);

describe('step-size convergence at default settings', () => {
  it('closed knocking cycle (RON charge, CR 9): default vs 4× refined steps', () => {
    const d = closed({});
    const f = closed({ maxStepDeg: 0.0625, fineStepDeg: 0.0125, knockBurnStep: 0.5e-6 });
    const ff = closed({ maxStepDeg: 0.03125, fineStepDeg: 0.00625, knockBurnStep: 0.25e-6 });
    const line = (k: keyof Res, deg: boolean) =>
      deg ? `${k} ${(d[k] - ff[k]).toFixed(4)}° / ${(f[k] - ff[k]).toFixed(4)}°` : `${k} ${(100 * (d[k] / ff[k] - 1)).toFixed(3)} % / ${(100 * (f[k] / ff[k] - 1)).toFixed(3)} %`;
    console.log(
      `[validation] closed CR 9, error of default / 4×-refined vs 8×-refined: ` +
        [line('W', false), line('pmax', false), line('ca10', true), line('ca50', true), line('ca90', true), line('ko', true), line('no', false), line('mapo', false)].join(', '),
    );
    expect(rel(d.W, f.W)).toBeLessThan(2e-3);
    expect(rel(d.pmax, f.pmax)).toBeLessThan(2e-3);
    expect(rel(d.no, f.no)).toBeLessThan(2e-3);
    for (const k of ['ca10', 'ca50', 'ca90', 'ko'] as const) expect(Math.abs(d[k] - f[k])).toBeLessThan(0.1);
    expect(rel(d.mapo, f.mapo)).toBeLessThan(2e-3); // FAILS in round 1 (≈ 2 %)
  });

  it('fired RON cycle at the periodic state with burnRateMultiplier 2.6: default vs 4× refined steps', () => {
    const run = (o: Partial<CycleModelOptions>): CycleSummary => new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 12, burnRateMultiplier: 2.6, ...o }).runCycles(1)[0];
    const d = run({});
    const f = run({ maxStepDeg: 0.0625, fineStepDeg: 0.0125, knockBurnStep: 0.5e-6 });
    const pct = (k: keyof CycleSummary) => `${k} ${(100 * ((d[k] as number) / (f[k] as number) - 1)).toFixed(3)} %`;
    const dg = (k: keyof CycleSummary) => `${k} ${((d[k] as number) - (f[k] as number)).toFixed(4)}°`;
    console.log(`[validation] RON ×2.6 default vs 4× refined: ${[pct('imepNet'), pct('peakPressure'), pct('noPpm'), pct('mapo'), pct('maxPressureRiseRate'), dg('ca10'), dg('ca50'), dg('ca90'), dg('knockOnsetDeg')].join(', ')}`);
    expect(rel(d.imepNet, f.imepNet)).toBeLessThan(2e-3);
    expect(rel(d.peakPressure, f.peakPressure)).toBeLessThan(2e-3);
    expect(rel(d.noPpm, f.noPpm)).toBeLessThan(2e-3);
    for (const k of ['ca10', 'ca50', 'ca90', 'knockOnsetDeg'] as const) expect(Math.abs(d[k] - f[k])).toBeLessThan(0.1);
    expect(rel(d.mapo, f.mapo)).toBeLessThan(2e-3); // FAILS in round 1 (≈ 2 %)
    expect(rel(d.maxPressureRiseRate, f.maxPressureRiseRate)).toBeLessThan(2e-3); // FAILS in round 1 (≈ 8 %)
  });
});

describe('event location', () => {
  it('the kernel hand-off is located inside the step: knock onset and CA50 fall strictly with the spark advance', () => {
    const adv = [13.0, 13.02, 13.04, 13.06, 13.08, 13.1];
    const r = adv.map((a) => closed({}, a));
    console.log(`[validation] spark 13.00…13.10°: CA50 ${r.map((x) => x.ca50.toFixed(4)).join(' ')}; knock ${r.map((x) => x.ko.toFixed(4)).join(' ')}; MAPO ${r.map((x) => (x.mapo / 1e5).toFixed(3)).join(' ')} bar`);
    for (let i = 1; i < r.length; i++) {
      expect(r[i].ca50).toBeLessThan(r[i - 1].ca50 - 1e-3); // FAILS in round 1 (0.05° stairs)
      expect(r[i].ko).toBeLessThan(r[i - 1].ko - 1e-3);
    }
  });
});

describe('independence of output sampling and warm-up', () => {
  it('MAPO does not depend on the output sampling (default steps vs 1 µs forced sampling after onset)', () => {
    const d = closed({});
    const s = closed({}, 13, 1e-6);
    console.log(`[validation] MAPO default ${(d.mapo / 1e5).toFixed(4)} bar vs 1 µs sampling ${(s.mapo / 1e5).toFixed(4)} bar (${(100 * (d.mapo / s.mapo - 1)).toFixed(2)} %)`);
    expect(rel(d.mapo, s.mapo)).toBeLessThan(5e-3); // FAILS in round 1 (≈ −4 %)
  });

  it('CycleSummary does not depend on the snapshot cadence (runCycles vs EngineSimulator 0.5°)', () => {
    const opts = { warmupCycles: 3, burnRateMultiplier: 2.6 };
    const a = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 0.5, ...opts }).runCycles(1)[0];
    const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 0.5, ...opts });
    const b: CycleSummary[] = [];
    while (b.length < 1) {
      sim.advanceToNextSnapshot();
      b.push(...sim.drainCycleSummaries());
    }
    const s = b[0];
    console.log(`[validation] runCycles vs 0.5° snapshots: MAPO ${(a.mapo / 1e5).toFixed(4)} / ${(s.mapo / 1e5).toFixed(4)} bar, max dp/dθ ${(a.maxPressureRiseRate / 1e5).toFixed(4)} / ${(s.maxPressureRiseRate / 1e5).toFixed(4)} bar/°, NO ${a.noPpm.toFixed(2)} / ${s.noPpm.toFixed(2)} ppm, CA50 ${a.ca50.toFixed(4)} / ${s.ca50.toFixed(4)}°`);
    expect(rel(a.imepNet, s.imepNet)).toBeLessThan(1e-4);
    expect(Math.abs(a.ca50 - s.ca50)).toBeLessThan(0.01);
    expect(rel(a.noPpm, s.noPpm)).toBeLessThan(5e-3);
    expect(rel(a.mapo, s.mapo)).toBeLessThan(5e-3); // FAILS in round 1 (≈ 1.5 %)
    expect(rel(a.maxPressureRiseRate, s.maxPressureRiseRate)).toBeLessThan(5e-3); // FAILS in round 1 (≈ 6 %)
  });

  it('the first emitted cycle after the default warm-up is at the periodic steady state (RON)', () => {
    const first = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, {}).runCycles(1)[0];
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 0 });
    const s = m.runCycles(30);
    const pss = s[29];
    expect(Math.abs(s[29].imepNet / s[28].imepNet - 1)).toBeLessThan(2e-5); // (the 30th cycle is converged)
    console.log(
      `[validation] first emitted cycle (warmupCycles 3) vs periodic state: IMEP ${(100 * (first.imepNet / pss.imepNet - 1)).toFixed(2)} %, p_max ${(100 * (first.peakPressure / pss.peakPressure - 1)).toFixed(2)} %, ` +
        `NO ${(100 * (first.noPpm / pss.noPpm - 1)).toFixed(1)} %, knock onset ${(first.knockOnsetDeg - pss.knockOnsetDeg).toFixed(2)}°, MAPO ${(100 * (first.mapo / pss.mapo - 1)).toFixed(0)} %; exhaust plenum ${m.exhaust.state.T.toFixed(0)} K`,
    );
    expect(rel(first.imepNet, pss.imepNet)).toBeLessThan(1e-3); // FAILS in round 1 (−0.54 %)
    expect(rel(first.peakPressure, pss.peakPressure)).toBeLessThan(2e-3);
    expect(Math.abs(first.knockOnsetDeg - pss.knockOnsetDeg)).toBeLessThan(0.1);
    expect(rel(first.noPpm, pss.noPpm)).toBeLessThan(1e-2);
  });
});
