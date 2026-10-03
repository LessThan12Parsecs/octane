/**
 * Validation round 2 (numerics, robustness, performance): regression tests for the defects found in
 * round 2 at the shipped defaults (CFR_CALIBRATION). Each test marked "FAILS in round 2" failed on the
 * code as delivered by the round-1 fixer pass; the numbers in the comments are the round-2 measurements.
 * Runtime ≈ 40–80 s (vitest).
 *
 *  1. Burned-zone NO is integrated with a FIRST-order operator split (cycle-model.ts afterStep §4:
 *     ZeldovichKinetics.advanceRateControlled over the whole step with the step-END equilibrium state).
 *     The in-cycle NO is 4.6 % off at 5° ATDC at the default 0.25° steps (halving each refinement), and
 *     the cycle-summary NO only converges by cancellation (formation-phase error +, decomposition-phase
 *     error −): refining only part of the cycle — exactly what the EngineSimulator's dense knock
 *     sampling (10 µs for 3 ms) does — moves the summary NO by +0.43…+0.49 % (> the 0.2 % criterion).
 *  2. KnockOscillator.setEndGasRegion projects the modes on the end-gas planform with a fixed 48 × 96
 *     midpoint quadrature: an end gas thinner than one cell row gives sourceShape = 0 and MAPO = 0
 *     EXACTLY although a knock onset is reported (RON PRF 90 at its guide-table CR 6.43: onset 13.06°,
 *     0.25 % end gas, MAPO 0 vs 0.116 bar with a 10× finer quadrature), and the coarse quadrature biases
 *     MAPO by 1–4 % elsewhere (CR 6.5: +3.5 %).
 *  3. The kinetic NO is swapped into the burned gas only at EVO (onEvo/swapNO) at constant U; for the
 *     rich CFR mixtures its oxygen comes from CO2 → CO, so T drops 11–15 K and p drops 0.8–1.2 % in one
 *     instant (energy is conserved, but the closed-phase thermodynamics never saw the NO).
 *  4. Free-speed mode: the step limit is set in TIME from ω at the step start (maxStepDeg/ω), the
 *     grid landing in ANGLE; while decelerating, the full step stops just short of the grid point and a
 *     sliver step (≪ 1e-3°) follows — 36 % of all steps (4,782 vs 3,079 per cycle), +50 % CPU.
 *  5. Warm-up (relaxPlenum) extrapolates plenum temperature/composition but not the plenum mass: at
 *     throttle ≤ 0.1 the intake plenum pumps down over 17–40 cycles and the first emitted cycle is
 *     26–30 % off in IMEP (the "converged first cycle" claim holds for throttle ≥ 0.2 only).
 *  6. Default knock model (Douaud–Eyzat, τ ∝ (ON/100)^3.4) gives τ = 0 for PRF 0 / n-heptane: the
 *     whole charge "autoignites" on the first step after IVC (−151.75°, 370 K, 1 bar), IMEP −6.5 bar,
 *     reported as misfire; PRF 10 at −74.6°. The LLNL tables put the onset at +8…+10° ATDC.
 *  7. sanitizeOperatingPoint passes NaN / undefined through clamp(): NaN rpm or φ throw inside the model
 *     (the worker then stops), an undefined spark advance silently disables the spark.
 *  8. Deterministic period-2 misfire (kernel quenched at r ≈ 1.5 mm after 42 mJ, then the low-residual
 *     cycle fires) at a moderate part-load point: 1200 rpm, throttle 0.5, φ 0.9 (S_L 0.78 m/s) — the
 *     ignition kernel's stretch/quench criterion (ignition/kernel.ts), not the cycle integrator.
 *  9. (passes) robustness regression over the UI control envelope and random mid-cycle operating-point
 *     changes / resets with every sub-model option.
 * The two-zone path itself is validated against an independent Cantera oracle in
 * numerics_two_zone_wiebe.test.ts (agreement 1e-9 … 1e-6, RK4 fourth order).
 */
import { describe, expect, it } from 'vitest';
import type { CycleSummary } from '../../src/physics/core/snapshot';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { CycleModel, EngineSimulator, MODE_OPEN, sanitizeOperatingPoint, type CycleModelOptions } from '../../src/physics/cycle/index';

const KNOCK_OP: OperatingPoint = { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 80 }, compressionRatio: 7 };
const rel = (a: number, b: number): number => Math.abs(a / b - 1);

describe('1. burned-zone NO integration (operator split)', () => {
  it('summary NO does not depend on the snapshot cadence (runCycles vs the worker path, knocking PRF 80 CR 7)', () => {
    const a = new EngineSimulator(CFR_F1, KNOCK_OP, { snapshotEveryDeg: 0.5 }).runCycles(1)[0];
    const sim = new EngineSimulator(CFR_F1, KNOCK_OP, { snapshotEveryDeg: 0.5 });
    let b: CycleSummary | undefined;
    while (!b) {
      sim.advanceToNextSnapshot();
      b = sim.drainCycleSummaries()[0];
    }
    console.log(`[validation] round 2: NO runCycles ${a.noPpm.toFixed(2)} ppm vs 0.5° snapshots ${b.noPpm.toFixed(2)} ppm (${(100 * (b.noPpm / a.noPpm - 1)).toFixed(3)} %)`);
    expect(rel(b.noPpm, a.noPpm)).toBeLessThan(2e-3); // FAILS in round 2 (+0.49 %)
  });

  it('in-cycle burned-zone NO at 5° ATDC is converged at the default steps (vs 16× refined)', () => {
    const nNOat = (o: Partial<CycleModelOptions>): number => {
      const m = new CycleModel(CFR_F1, KNOCK_OP, { knock: false, ...o });
      const c0 = m.cycle;
      m.stepUntil(Infinity, 360); // (wrap of the warm-up leftovers / spark interrupt)
      while (m.cycle === c0 || m.theta < 5 - 1e-9) {
        m.stepUntil(Infinity, m.cycle === c0 ? 360 : 5);
        if (m.cycle > c0 && Math.abs(m.theta - 5) < 1e-9) break;
      }
      return m.nNO;
    };
    const d = nNOat({});
    const f = nNOat({ maxStepDeg: 0.25 / 16, fineStepDeg: 0.05 / 16 });
    console.log(`[validation] round 2: burned-zone NO at 5° ATDC default vs 16× refined: ${(100 * (d / f - 1)).toFixed(2)} %`);
    expect(rel(d, f)).toBeLessThan(1e-2); // FAILS in round 2 (+4.6 %, first order)
  });
});

describe('2. MAPO: end-gas source-region quadrature', () => {
  /** MAPO with the shipped quadrature and with a 10× finer one (KnockOscillator private `quad`). */
  function mapoPair(cr: number): { d: CycleSummary; f: CycleSummary } {
    const d = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: cr }, {}).runCycles(1)[0];
    const mf = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: cr }, {});
    const osc = mf.knockOsc as unknown as { quad?: readonly [number, number] };
    if (osc.quad) osc.quad = [480, 960];
    const f = mf.runCycles(1)[0];
    return { d, f };
  }
  it('knock reported at the RON guide-table state (PRF 90, CR 6.43) has MAPO > 0 and converged', () => {
    const { d, f } = mapoPair(6.43);
    console.log(`[validation] round 2: RON CR 6.43 knock onset ${d.knockOnsetDeg.toFixed(2)}°, end gas ${(100 * d.knockEndGasFraction).toFixed(2)} %, MAPO ${(d.mapo / 1e5).toFixed(4)} bar (10× finer quadrature ${(f.mapo / 1e5).toFixed(4)} bar)`);
    expect(Number.isFinite(d.knockOnsetDeg)).toBe(true);
    expect(d.mapo).toBeGreaterThan(0); // FAILS in round 2 (exactly 0)
    expect(rel(d.mapo, f.mapo)).toBeLessThan(1e-2);
  });
  it('MAPO at CR 6.5 is within 1 % of the 10× finer source quadrature', () => {
    const { d, f } = mapoPair(6.5);
    console.log(`[validation] round 2: RON CR 6.5 MAPO ${(d.mapo / 1e5).toFixed(4)} vs ${(f.mapo / 1e5).toFixed(4)} bar (${(100 * (d.mapo / f.mapo - 1)).toFixed(2)} %)`);
    expect(rel(d.mapo, f.mapo)).toBeLessThan(1e-2); // FAILS in round 2 (+3.5 %)
  });
});

describe('3. EVO: kinetic-NO swap', () => {
  it('no pressure discontinuity at EVO (RON)', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, {});
    m.recordTrace();
    m.runCycles(1);
    const tr = m.trace!;
    const at: number[] = [];
    for (let i = 0; i < tr.theta.length; i++) if (Math.abs(tr.theta[i] - m.evoDeg) < 1e-9) at.push(i);
    const a = at.find((i) => tr.mode[i] !== MODE_OPEN)!;
    const b = at.find((i) => tr.mode[i] === MODE_OPEN)!;
    console.log(`[validation] round 2: EVO p ${(tr.pressure[a] / 1e5).toFixed(4)} → ${(tr.pressure[b] / 1e5).toFixed(4)} bar, T ${tr.Tmean[a].toFixed(1)} → ${tr.Tmean[b].toFixed(1)} K (U ${tr.U[a].toFixed(4)} → ${tr.U[b].toFixed(4)} J)`);
    expect(Math.abs(tr.U[b] - tr.U[a])).toBeLessThan(1e-9 * Math.abs(tr.U[a]));
    expect(rel(tr.pressure[b], tr.pressure[a])).toBeLessThan(1e-3); // FAILS in round 2 (−0.82 %)
  });
});

describe('4. free-speed stepping', () => {
  it('no sliver steps: < 5 % of the steps shorter than 1e-3° (free mode, load 5 N m)', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, speedMode: 'free', loadTorque: 5 }, {});
    m.recordTrace();
    m.runCycles(1);
    const th = m.trace!.theta;
    let n = 0;
    let sliver = 0;
    for (let i = 1; i < th.length; i++) {
      const d = th[i] - th[i - 1];
      if (d <= 0) continue;
      n++;
      if (d < 1e-3) sliver++;
    }
    console.log(`[validation] round 2: free mode ${n} steps per cycle, ${sliver} (${((100 * sliver) / n).toFixed(1)} %) shorter than 1e-3°`);
    expect(sliver / n).toBeLessThan(0.05); // FAILS in round 2 (36 %)
  });
});

describe('5. warm-up at part throttle', () => {
  it('the first emitted cycle at throttle 0.1 is within 1 % (IMEP) of the periodic state', () => {
    const op = { ...CFR_RON_CONDITIONS, throttle: 0.1 };
    const first = new CycleModel(CFR_F1, op, {}).runCycles(1)[0];
    const s = new CycleModel(CFR_F1, op, { warmupCycles: 0 }).runCycles(30);
    console.log(`[validation] round 2: throttle 0.1 first emitted IMEP ${(first.imepNet / 1e5).toFixed(3)} bar vs cycle 30 ${(s[29].imepNet / 1e5).toFixed(3)} bar`);
    expect(rel(s[29].imepNet, s[28].imepNet)).toBeLessThan(1e-3);
    expect(rel(first.imepNet, s[29].imepNet)).toBeLessThan(1e-2); // FAILS in round 2 (+26 %)
  });
});

describe('6. low-octane fuels with the default knock model', () => {
  for (const fuel of [{ kind: 'PRF', octaneNumber: 0 }, { kind: 'pure', species: 'NC7H16' }, { kind: 'PRF', octaneNumber: 10 }] as const) {
    it(`${JSON.stringify(fuel)} at RON conditions does not autoignite during early compression`, () => {
      const s = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel }, { warmupCycles: 1 }).runCycles(1)[0];
      console.log(`[validation] round 2: ${JSON.stringify(fuel)} knock onset ${s.knockOnsetDeg.toFixed(2)}°, IMEP ${(s.imepNet / 1e5).toFixed(2)} bar, misfire ${s.misfire}`);
      // (an end-gas autoignition before −30° at 52 °C intake and CR 6.43 is not physical; LLNL tables: +8…+10°)
      expect(!(s.knockOnsetDeg < -30)).toBe(true); // FAILS in round 2 (−151.75° / −151.75° / −74.6°)
    });
  }
});

describe('7. operating-point sanitiser', () => {
  it('non-finite / missing fields are repaired, not passed through', () => {
    const bad = { ...CFR_RON_CONDITIONS, rpm: NaN, equivalenceRatio: NaN, sparkAdvanceDeg: undefined as unknown as number };
    const s = sanitizeOperatingPoint(CFR_F1, bad);
    const vals = [s.rpm, s.equivalenceRatio, s.sparkAdvanceDeg];
    console.log(`[validation] round 2: sanitised rpm ${s.rpm}, φ ${s.equivalenceRatio}, spark ${s.sparkAdvanceDeg}`);
    for (const v of vals) expect(Number.isFinite(v)).toBe(true); // FAILS in round 2 (NaN, NaN, undefined)
  });
});

describe('8. misfire limit cycle at a moderate part-load point (ignition kernel)', () => {
  it('1200 rpm, throttle 0.5, φ 0.9: consecutive periodic cycles agree within 10 % IMEP', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, rpm: 1200, throttle: 0.5, equivalenceRatio: 0.9 }, {});
    const s = m.runCycles(4);
    console.log(`[validation] round 2: 1200 rpm thr 0.5 φ 0.9 IMEP ${s.map((x) => (x.imepNet / 1e5).toFixed(2)).join(' / ')} bar, misfire ${s.map((x) => (x.misfire ? 1 : 0)).join('')}`);
    expect(rel(s[3].imepNet, s[2].imepNet)).toBeLessThan(0.1); // FAILS in round 2 (−0.51 / 5.96 bar alternating)
  });
});

// ---------------------------------------------------------------------------------------------
// Robustness regression (passes in round 2): the UI control envelope (src/ui/controls.ts ranges)
// and random operating-point changes / resets mid-cycle with every sub-model option. No exception,
// no non-finite summary or snapshot field (NaN allowed by definition only for CA10/50/90, ISFC,
// knock onset). Round 2: 40 + 40 random points and 18 randomized runs (129k snapshots) were clean.
// ---------------------------------------------------------------------------------------------

function nonFinite(o: unknown, path = '', out: string[] = []): string[] {
  if (typeof o === 'number') {
    if (!Number.isFinite(o)) out.push(path);
  } else if (Array.isArray(o)) o.forEach((v, i) => nonFinite(v, `${path}[${i}]`, out));
  else if (o && typeof o === 'object') for (const [k, v] of Object.entries(o)) nonFinite(v, path ? `${path}.${k}` : k, out);
  return out;
}

function mulberry32(a: number): () => number {
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

describe('9. robustness over the UI envelope and mid-run changes (passes in round 2)', () => {
  const PURE = ['IC8H18', 'NC7H16', 'CH4', 'C3H8', 'C2H5OH'] as const;
  it('12 seeded points over the UI control ranges: no exception, finite summaries and snapshots', () => {
    const rnd = mulberry32(7);
    const U = (a: number, b: number): number => a + (b - a) * rnd();
    const problems: string[] = [];
    for (let i = 0; i < 12; i++) {
      const free = rnd() < 0.25;
      const op: OperatingPoint = {
        ...CFR_RON_CONDITIONS,
        fuel: rnd() < 0.5 ? { kind: 'PRF', octaneNumber: Math.round(U(0, 100)) } : { kind: 'pure', species: PURE[Math.floor(rnd() * 5)] },
        equivalenceRatio: U(0.4, 2.0),
        compressionRatio: U(4, 18),
        rpm: U(200, 3000),
        sparkAdvanceDeg: U(-10, 60),
        egrFraction: U(0, 0.4),
        throttle: U(0, 1),
        intakeMixtureTemperature: U(273, 473),
        ambientPressure: U(0.5e5, 1.1e5),
        ambientTemperature: U(243, 323),
        relativeHumidity: U(0, 1),
        coolantTemperature: U(293, 403),
        dwellTime: U(0.5e-3, 10e-3),
        speedMode: free ? 'free' : 'fixed',
        loadTorque: free ? U(-20, 100) : 0,
      };
      try {
        const sim = new EngineSimulator(CFR_F1, op, { snapshotEveryDeg: 2, warmupCycles: 1 });
        let n = 0;
        while (sim.model.cycle < 2 && n++ < 200_000) {
          const s = sim.advanceToNextSnapshot();
          const nf = nonFinite(s);
          if (nf.length) {
            problems.push(`#${i} snapshot ${nf.join(',')} @${s.thetaDeg.toFixed(1)}`);
            break;
          }
          for (const c of sim.drainCycleSummaries()) {
            const nf2 = nonFinite({ ...c, ca10: 0, ca50: 0, ca90: 0, isfc: 0, knockOnsetDeg: 0 });
            if (nf2.length) problems.push(`#${i} summary ${nf2.join(',')}`);
          }
        }
      } catch (e) {
        problems.push(`#${i} threw ${String(e).slice(0, 160)} ${JSON.stringify(op)}`);
      }
    }
    expect(problems).toEqual([]);
  });

  it('random operating-point changes and resets mid-cycle, every sub-model option', () => {
    const rnd = mulberry32(99);
    const U = (a: number, b: number): number => a + (b - a) * rnd();
    const OPTS: Partial<CycleModelOptions>[] = [
      {},
      { ignitionDelayModel: 'llnl-gasoline-2011' },
      { ignitionDelayModel: 'llnl-prf-v2', knockIntegral: 'two-stage' },
      { turbulentFlameClosure: 'keck1982' },
      { combustionModel: 'wiebe' },
      { combustionModel: 'instantaneous-at-tdc' },
    ];
    const patch = (): Partial<OperatingPoint> => {
      switch (Math.floor(rnd() * 8)) {
        case 0: return { rpm: U(200, 3000) };
        case 1: return { throttle: U(0, 1) };
        case 2: return { equivalenceRatio: U(0.4, 2) };
        case 3: return { sparkAdvanceDeg: U(-10, 60) };
        case 4: return { compressionRatio: U(4, 18) };
        case 5: return { fuel: rnd() < 0.5 ? { kind: 'PRF', octaneNumber: U(0, 100) } : { kind: 'pure', species: PURE[Math.floor(rnd() * 5)] } };
        case 6: return { speedMode: rnd() < 0.5 ? 'free' : 'fixed', loadTorque: U(-20, 100) };
        default: return { egrFraction: U(0, 0.4) };
      }
    };
    const problems: string[] = [];
    for (const o of OPTS) {
      const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, { snapshotEveryDeg: 1, warmupCycles: 1, ...o });
      try {
        let n = 0;
        while (sim.model.cycle < 4 && n++ < 100_000) {
          const s = sim.advanceToNextSnapshot();
          const nf = nonFinite(s);
          if (nf.length) {
            problems.push(`${JSON.stringify(o)} snapshot ${nf.join(',')} @${s.thetaDeg.toFixed(1)}`);
            break;
          }
          sim.drainCycleSummaries();
          if (rnd() < 0.005) sim.setOperatingPoint(patch());
          if (rnd() < 0.0004) sim.reset();
        }
      } catch (e) {
        problems.push(`${JSON.stringify(o)} threw ${String(e).slice(0, 160)}`);
      }
    }
    expect(problems).toEqual([]);
  });
});
