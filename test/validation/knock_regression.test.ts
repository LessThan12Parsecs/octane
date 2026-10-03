/**
 * Validation round 1 (knock & octane rating) — fast regression checks of the knock path of the
 * cycle model (src/physics/cycle). Each test states the defect it guards against; all run through
 * the public API only. Runtime ≈ 15–25 s.
 *
 *  1. The reported (knock-oscillation-including) pressure must stay finite after the burn-out
 *     merge of a knocking cycle. Found: cycle-model.ts step 5 (knock oscillator) evaluates
 *     `cl.heatOfReaction(this.Tu > 0 ? this.Tu : cl.Tu)` AFTER the burn-out merge of the same step
 *     has set Tu = cl.Tu = 0 (MODE_BURNED) → heatOfReaction(0) = NaN → every modal amplitude is NaN
 *     until EVO: EngineSnapshot.pressure / knock.oscillation are NaN from burn-out to EVO, and MAPO
 *     silently ignores everything after the merge (NaN comparisons are false).
 *  2. MAPO must stay within the physically observed range of CFR knock. CFR data: standard knock
 *     ≈ 0.6–0.9 bar (4–18 kHz, Hoth & Kolodziej 2021 Part 1 Tables 6–7, OSTI 1880351), 2.3 bar
 *     (AVL IndiCom, Pal et al. 2018 Fig. 1), heavy knock ≤ 3 bar mean / 8 bar single cycles
 *     (Rockstroh et al. 2018, SAE 2018-01-0210). Found: when the knock onset state sits at the NTC
 *     turning point of the ignition delay, |∂lnτ/∂T| → 0 and autoignitionBurnTime → τ_e ≈ 1 µs
 *     (knock.ts), the whole end gas is released as a constant-volume explosion, and MAPO jumps from
 *     0 to 60–115 bar between compression ratios 0.25 apart.
 *  3. CycleSummary.knockEndGasFraction ("end-gas mass fraction at onset") must be the unburned gas
 *     AHEAD of the flame, (m − m_e)/m. Found: it is m_u/m, which includes the entrained-but-unburned
 *     brush mass (m_e − m_b) behind the front — 4–9 percentage points too high at standard knock.
 *     (Fixer round 2: with the default knockBrushAutoignition the brush autoignites too and the reported
 *     fraction is m_u/m by design; this test pins knockBrushAutoignition false.)
 */
import { describe, expect, it } from 'vitest';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { CycleModel, EngineSimulator } from '../../src/physics/cycle/index';
import type { OperatingPoint } from '../../src/physics/core/operating-point';

const prf = (on: number, cr: number, patch: Partial<OperatingPoint> = {}): OperatingPoint => ({
  ...CFR_RON_CONDITIONS,
  fuel: { kind: 'PRF', octaneNumber: on },
  compressionRatio: cr,
  ...patch,
});

describe('knock path regressions', () => {
  it('reported pressure stays finite after the burn-out merge of a knocking cycle (snapshots + trace)', () => {
    // The author's own "burn-rate multiplier 2.6" example: PRF 90 at the ASTM guide-table CR knocks
    // at ≈ 23° and burns out at ≈ 32°.
    const sim = new EngineSimulator(CFR_F1, prf(90, 6.43), { snapshotEveryDeg: 0.5, burnRateMultiplier: 2.6 });
    const c0 = sim.model.cycle;
    let bad = 0;
    let firstBad = NaN;
    let n = 0;
    while (sim.model.cycle < c0 + 1) {
      const s = sim.advanceToNextSnapshot();
      n++;
      if (!Number.isFinite(s.pressure) || !Number.isFinite(s.knock.oscillation)) {
        bad++;
        if (Number.isNaN(firstBad)) firstBad = s.thetaDeg;
      }
    }
    const sum = sim.drainCycleSummaries()[0];
    expect(Number.isFinite(sum.knockOnsetDeg)).toBe(true); // the case does knock
    if (bad > 0) console.log(`[knock] ${bad} of ${n} snapshots with non-finite pressure, first at θ = ${firstBad.toFixed(2)}°`);
    expect(bad).toBe(0);
    // full-resolution trace of the next cycle
    const m = sim.model;
    m.recordTrace();
    m.runCycles(1);
    const tr = m.trace!;
    const nonFinite = tr.pressureReported.filter((v) => !Number.isFinite(v)).length;
    expect(nonFinite).toBe(0);
  });

  it('MAPO stays within the CFR range (< 20 bar) and varies continuously with CR near standard knock', () => {
    // Phasing-calibrated burn rate (CA50 ≈ 8–9° aTDC at RON, cf. Kolodziej & Wallner 2017 Fig. 9:
    // CA50 8.8°, knock point 12.0° for PRF98 at standard knock).
    const crs = [6.75, 7.0, 7.25, 7.5];
    const m = new CycleModel(CFR_F1, prf(90, crs[0]), { burnRateMultiplier: 3.9, warmupCycles: 3 });
    const mapo: number[] = [];
    for (const cr of crs) {
      m.setOperatingPoint({ compressionRatio: cr });
      m.runCycles(1);
      mapo.push(m.runCycles(1)[0].mapo);
    }
    console.log(`[knock] PRF 90 RON, burnRateMultiplier 3.9: MAPO(CR ${crs.join('/')}) = ${mapo.map((v) => (v / 1e5).toFixed(2)).join(' / ')} bar`);
    for (const v of mapo) expect(v).toBeLessThan(20e5);
    for (let i = 1; i < mapo.length; i++) {
      if (mapo[i - 1] > 0.5e5) expect(mapo[i] / mapo[i - 1]).toBeLessThan(4);
    }
  });

  it('knockEndGasFraction is the unburned mass ahead of the flame front, (m − m_e)/m at onset', () => {
    // fixer round 2: the default now autoignites ALL unburned gas (options.knockBrushAutoignition; the
    // reported fraction is the autoigniting gas, m_u/m — pinned in code-review.test.ts). The definition
    // guarded here — the gas ahead of the front when only that gas autoignites — is the option's
    // `false` setting (expectation unchanged; option pinned deliberately).
    const m = new CycleModel(CFR_F1, prf(90, 7.5), { burnRateMultiplier: 3.9, warmupCycles: 3, knockBrushAutoignition: false });
    m.recordTrace();
    const s = m.runCycles(1)[0];
    expect(Number.isFinite(s.knockOnsetDeg)).toBe(true);
    const tr = m.trace!;
    // the sample AT the onset (the step is redone to end exactly at the Livengood–Wu crossing and
    // recorded there). Round 1 took the last sample before it, one step earlier — with the calibrated
    // (faster) flame the front entrains ≈ 1 % of the charge per 0.25° step, which alone exceeded the
    // 0.01 tolerance (fixer pass: comparison point corrected, tolerance unchanged).
    let i = 0;
    while (i + 1 < tr.theta.length && tr.theta[i + 1] <= s.knockOnsetDeg + 1e-9) i++;
    const mTot = tr.mu[i] + tr.mb[i];
    const ahead = (mTot - tr.me[i]) / mTot;
    const withBrush = tr.mu[i] / mTot;
    console.log(`[knock] onset ${s.knockOnsetDeg.toFixed(2)}°: reported ${s.knockEndGasFraction.toFixed(4)}, ahead of front ${ahead.toFixed(4)}, m_u/m ${withBrush.toFixed(4)}`);
    expect(Math.abs(s.knockEndGasFraction - ahead)).toBeLessThan(0.01);
  });
});
