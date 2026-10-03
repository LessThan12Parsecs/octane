/**
 * Validation round 1 (numerics): robustness. Seeded random operating points over the whole
 * control envelope (φ 0.5–1.6, CR 4–18, 300–1500 rpm, spark 0–50° BTDC, PRF 0–100 and every
 * pure fuel, EGR 0–0.3, throttle 0.05–1, mixture 270–450 K, ambient 0.7–1.1 bar, free speed with
 * load) must run without exceptions and give finite summaries and finite snapshots (the render
 * contract); plus minimal reproductions of the NaN-snapshot defects found in round 1 and the
 * free-speed extremes. Runtime ≈ 30–60 s.
 *
 * Status after the round-1 fixer pass: all pass (fixes in cycle-model.ts afterStep / onAutoignition /
 * onIvc knockAvailable, engine-simulator.ts makeSnapshot).
 * Known failures at the time of writing (validation round 1):
 *  - 'snapshot pressure stays finite …': the knock oscillator is driven with q = NaN on the step
 *    where the burn-out merge happens during the end-gas burn (cycle-model.ts afterStep, qc =
 *    closure.heatOfReaction(0) = NaN because T_u = 0 in MODE_BURNED) → pressure / oscillation /
 *    gasTorque / netTorque NaN until EVO.
 *  - 'knock integral …non-PRF fuels': the PRF delay models return τ = NaN for CH4/C3H8/C2H5OH, the
 *    LW integral is NaN and EngineSnapshot.knock.integral is NaN.
 *  - 'knock-onset snapshot …': autoignition before any burned gas creates the burned zone in the
 *    split but the snapshot is taken before the closure is re-evaluated (T_b = 0, empty equilibrium
 *    result) → burnedComposition = 0/0.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSnapshot } from '../../src/physics/core/snapshot';
import type { FuelSelection, OperatingPoint } from '../../src/physics/core/operating-point';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { CycleModel, EngineSimulator, type CycleModelOptions } from '../../src/physics/cycle/index';

function nonFinite(o: unknown, path = '', out: string[] = []): string[] {
  if (typeof o === 'number') {
    if (!Number.isFinite(o)) out.push(path);
  } else if (Array.isArray(o)) o.forEach((v, i) => nonFinite(v, `${path}[${i}]`, out));
  else if (o && typeof o === 'object') for (const [k, v] of Object.entries(o)) nonFinite(v, path ? `${path}.${k}` : k, out);
  return out;
}

/** Snapshots of one full cycle (after `warmupCycles`); returns the set of non-finite field paths with the first angle. */
function snapshotScan(op: OperatingPoint, opts: Partial<CycleModelOptions> & { snapshotEveryDeg?: number }): Map<string, number> {
  const sim = new EngineSimulator(CFR_F1, op, { snapshotEveryDeg: 1, warmupCycles: 1, ...opts });
  const bad = new Map<string, number>();
  let s: EngineSnapshot = sim.advanceToNextSnapshot();
  let n = 0;
  while (s.cycle < 1 && n++ < 200_000) {
    for (const p of nonFinite(s)) if (!bad.has(p)) bad.set(p, s.thetaDeg);
    s = sim.advanceToNextSnapshot();
  }
  return bad;
}

const fmtBad = (b: Map<string, number>): string => [...b].map(([p, th]) => `${p}@${th.toFixed(2)}°`).join(', ');

describe('snapshot render contract (no NaN/∞) — minimal reproductions', () => {
  it('snapshot pressure stays finite when the end gas burns out after knock (RON, burnRateMultiplier 2.6)', () => {
    const bad = snapshotScan(CFR_RON_CONDITIONS, { burnRateMultiplier: 2.6 });
    if (bad.size) console.log(`[validation] RON ×2.6 non-finite snapshot fields: ${fmtBad(bad)}`);
    expect(fmtBad(bad)).toBe('');
  });

  it('snapshot pressure stays finite at the knock-limited CR 8 (uncalibrated burn rate)', () => {
    const bad = snapshotScan({ ...CFR_RON_CONDITIONS, compressionRatio: 8 }, {});
    if (bad.size) console.log(`[validation] RON CR 8 non-finite snapshot fields: ${fmtBad(bad)}`);
    expect(fmtBad(bad)).toBe('');
  });

  for (const species of ['CH4', 'C3H8', 'C2H5OH'] as const) {
    it(`knock integral is finite in snapshots for non-PRF fuels (${species})`, () => {
      const bad = snapshotScan({ ...CFR_RON_CONDITIONS, fuel: { kind: 'pure', species } }, { warmupCycles: 0 });
      if (bad.size) console.log(`[validation] ${species} non-finite snapshot fields: ${fmtBad(bad)}`);
      expect(fmtBad(bad)).toBe('');
    });
  }

  it('knock-onset snapshot has a finite burned composition when autoignition precedes the spark', () => {
    const op: OperatingPoint = {
      ...CFR_RON_CONDITIONS,
      rpm: 1338,
      throttle: 0.876,
      ambientPressure: 101673,
      intakeMixtureTemperature: 413.1,
      fuel: { kind: 'PRF', octaneNumber: 3 },
      equivalenceRatio: 0.863,
      sparkAdvanceDeg: 4.7,
      compressionRatio: 14.09,
      egrFraction: 0.075,
    };
    const bad = snapshotScan(op, { warmupCycles: 0 });
    if (bad.size) console.log(`[validation] PRF 3 autoignition-before-spark non-finite snapshot fields: ${fmtBad(bad)}`);
    expect(fmtBad(bad)).toBe('');
  });
});

describe('free speed extremes', () => {
  it('a load above the engine torque decelerates to the 60 rpm floor without NaN; a driving load accelerates', () => {
    const slow = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: 100 }, { warmupCycles: 0 });
    const s1 = slow.runCycles(6);
    expect(slow.rpm).toBeGreaterThan(59.5); // the stall guard acts per RHS evaluation (≈ 59.97 rpm observed)
    expect(slow.rpm).toBeLessThan(100);
    for (const s of s1) expect(nonFinite({ ...s, ca10: 0, ca50: 0, ca90: 0, isfc: 0, knockOnsetDeg: 0 })).toEqual([]);
    const fast = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: -30 }, { warmupCycles: 0 });
    fast.runCycles(4);
    expect(fast.rpm).toBeGreaterThan(800);
  });
});

function mulberry32(a: number): () => number {
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('seeded random sweep over the control envelope', () => {
  const PURE: FuelSelection[] = (['IC8H18', 'NC7H16', 'CH4', 'C3H8', 'C2H5OH'] as const).map((species) => ({ kind: 'pure', species }));
  const rnd = mulberry32(20260930);
  const U = (a: number, b: number): number => a + (b - a) * rnd();
  const ops: OperatingPoint[] = [];
  for (let i = 0; i < 24; i++) {
    const free = rnd() < 0.2;
    ops.push({
      ...CFR_RON_CONDITIONS,
      fuel: rnd() < 0.5 ? { kind: 'PRF', octaneNumber: Math.round(U(0, 100)) } : PURE[Math.floor(rnd() * PURE.length)],
      equivalenceRatio: U(0.5, 1.6),
      compressionRatio: U(4, 18),
      rpm: U(300, 1500),
      sparkAdvanceDeg: U(0, 50),
      egrFraction: U(0, 0.3),
      throttle: U(0.05, 1),
      intakeMixtureTemperature: U(270, 450),
      ambientPressure: U(0.7e5, 1.1e5),
      speedMode: free ? 'free' : 'fixed',
      loadTorque: free ? U(0, 15) : 0,
    });
  }
  it('no exceptions, finite summaries, finite snapshots (every 3rd point)', () => {
    const problems: string[] = [];
    let misfires = 0;
    ops.forEach((op, i) => {
      try {
        const m = new CycleModel(CFR_F1, op, { warmupCycles: 1 });
        for (const s of m.runCycles(2)) {
          // NaN is allowed by definition only for CA10/50/90 (not reached), isfc (W ≤ 0), knock onset (none)
          const nf = nonFinite({ ...s, ca10: 0, ca50: 0, ca90: 0, isfc: 0, knockOnsetDeg: 0 });
          if (nf.length) problems.push(`#${i} summary ${nf.join(',')}`);
          if (s.imepNet > 0 && !Number.isFinite(s.isfc)) problems.push(`#${i} isfc`);
          if (s.misfire) misfires++;
        }
        if (i % 3 === 0) {
          const bad = snapshotScan(op, { warmupCycles: 0 });
          if (bad.size) problems.push(`#${i} snapshot ${fmtBad(bad)}`);
        }
      } catch (e) {
        problems.push(`#${i} threw ${String(e).slice(0, 160)}`);
      }
    });
    console.log(`[validation] seeded sweep: ${ops.length} points, ${misfires}/${2 * ops.length} misfired cycles, ${problems.length} problems${problems.length ? ':\n  ' + problems.join('\n  ') : ''}`);
    expect(problems).toEqual([]);
  });
});
