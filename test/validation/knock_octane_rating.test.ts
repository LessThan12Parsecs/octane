/**
 * Validation (knock & octane rating) — the headline check: does the simulated CFR reproduce the
 * ASTM D2699/D2700 knock-limited compression ratio vs PRF octane number?
 *
 * [SLOW: ≈ 60–120 s. Runs only with OCTANE_VALIDATION=1, e.g.
 *    OCTANE_VALIDATION=1 npx vitest run test/validation/knock_octane_rating.test.ts ]
 * Uses the simulator DEFAULTS (the one global calibration, src/physics/cycle/calibration.ts),
 * through the public API only.
 *
 * Knock-intensity criterion (revised in the round-1 fixer pass): standard knock = MAPO at the
 * D-1 pickup of STANDARD_KNOCK_MAPO = 0.67 bar, the measured MAPO at ASTM standard knock
 * (Rockstroh et al. 2018 standard point 0.67 bar; Hoth & Kolodziej 2021 Part 1 Tables 6–7:
 * 0.63–0.69 bar for PRF98 at standard RON, 4–18 kHz — test/fixtures/validation_knock_targets.json).
 * Round 1 used the unburned mass fraction AHEAD of the flame front at the Livengood–Wu onset,
 * self-calibrated at PRF 90, with a plausibility window 0.2–0.45 anchored on "knock point at
 * 60–70 % of heat release". That proxy is no longer usable: (i) with the calibrated burn model the
 * flame brush holds ≈ 10 % of the charge at the onset, so the fraction ahead of the front is not the
 * unburned fraction; (ii) the grade-A Choi 2018 PRF100 standard-knock trace analysed with the suite's
 * own apparent-heat-release routine gives 72 % burned at 10° and 83 % at the measured knock point
 * 11.0° (Pal 2018), i.e. 17–28 % unburned, not 30–40 %; (iii) the model's MAPO is now finite,
 * sampling-independent and validated (1.50 bar vs 1.43 bar measured at the Choi/Hoth PRF100
 * standard-knock state). MAPO rises by ≈ 5–7 bar per unit CR near the knock limit, so the
 * knock-limited CR hardly depends on the threshold (0.3 → 1.43 bar moves it by ≈ 0.2 CR).
 * Each (method, ON) gets its own model (a fuel change takes ≈ 5 cycles to reach the cylinder
 * through the 1 L intake plenum — the model now uses the TRAPPED mixture); CR changes are
 * followed by two transition cycles.
 *
 * Known model-form limitations (it.fails, see calibration.ts CFR_KNOCK_DELAY_MODEL): the model knocks
 * too early for PRF ≤ 70 in RON. Round 1 also had MON PRF 70–95 0.35–0.6 CR too high and a 4× too
 * weak spark response (both pass since the round-2 fixer pass). Round 2 (fixer): the whole scale sits
 * 0.14–0.45 CR low at ON 80–100 (RON 80 −0.45, RON 100 −0.37, MON 100 −0.38 fail ±0.35 and are left
 * failing as open items — the end-gas stratification ΔT is at its upper bound, calibration.ts).
 */
import { beforeAll, describe, expect, it } from 'vitest';
import targets from '../fixtures/validation_knock_targets.json';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS, cfrMonSparkAdvanceDeg } from '../../src/physics/engines/cfr';
import { CycleModel } from '../../src/physics/cycle/index';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary } from '../../src/physics/core/snapshot';

const RUN = !!process.env.OCTANE_VALIDATION;
type Method = 'RON' | 'MON';
const guideCR = (m: Method, on: number): number =>
  (targets.guideTable[m] as Record<string, { compressionRatio: number }>)[String(on)].compressionRatio;

/** Standard knock intensity as MAPO, Pa (Rockstroh 2018 standard point). */
const STANDARD_KNOCK_MAPO = targets.standardKnockState.Rockstroh2018_PRF90_standard.mapoPa;
/** ± tolerance on the knock-limited CR (≈ ±5 ON at 90–100 RON, ±3.5 ON at MON 90). */
const CR_TOL = 0.35;

/** Knock-limited CR (MAPO = STANDARD_KNOCK_MAPO) of one PRF, by bisection; one model per call. */
function knockLimitedCR(method: Method, on: number, patch: Partial<OperatingPoint> = {}): { cr: number; s: CycleSummary } {
  const base = method === 'RON' ? CFR_RON_CONDITIONS : CFR_MON_CONDITIONS;
  const opAt = (cr: number): OperatingPoint => {
    const op: OperatingPoint = { ...base, fuel: { kind: 'PRF', octaneNumber: on }, compressionRatio: cr, ...patch };
    if (method === 'MON' && patch.sparkAdvanceDeg === undefined) op.sparkAdvanceDeg = cfrMonSparkAdvanceDeg(cr);
    return op;
  };
  const g = guideCR(method, on);
  let lo = g - 1.5;
  let hi = g + 1.5;
  const m = new CycleModel(CFR_F1, opAt(g), { warmupCycles: 3 });
  let s: CycleSummary = m.runCycles(1)[0];
  const mapoAt = (cr: number): number => {
    m.setOperatingPoint(opAt(cr));
    m.runCycles(2); // transition
    s = m.runCycles(1)[0];
    return s.mapo;
  };
  for (let i = 0; i < 8; i++) {
    const mid = 0.5 * (lo + hi);
    if (mapoAt(mid) < STANDARD_KNOCK_MAPO) lo = mid;
    else hi = mid;
  }
  const cr = 0.5 * (lo + hi);
  mapoAt(cr);
  return { cr, s };
}

const fmt = (s: CycleSummary): string =>
  `onset ${Number.isFinite(s.knockOnsetDeg) ? s.knockOnsetDeg.toFixed(1) : '—'}° (end gas ${(100 * s.knockEndGasFraction).toFixed(1)} %), ` +
  `CA10/50/90 ${s.ca10.toFixed(1)}/${s.ca50.toFixed(1)}/${s.ca90.toFixed(1)}°, p_max ${(s.peakPressure / 1e5).toFixed(1)} bar @ ${s.peakPressureDeg.toFixed(1)}°, ` +
  `MAPO ${(s.mapo / 1e5).toFixed(2)} bar, dp/dθ_max ${(s.maxPressureRiseRate / 1e5).toFixed(2)} bar/°, gIMEP ${(s.imepGross / 1e5).toFixed(2)} bar`;

describe.skipIf(!RUN)('[slow] knock-limited CR vs PRF octane number (ASTM guide tables), simulator defaults', () => {
  const kl: Record<string, { cr: number; s: CycleSummary }> = {};
  const cases: [Method, number][] = [
    ['RON', 60],
    ['RON', 70],
    ['RON', 80],
    ['RON', 90],
    ['RON', 100],
    ['MON', 70],
    ['MON', 90],
    ['MON', 100],
  ];
  beforeAll(() => {
    for (const [m, on] of cases) {
      kl[`${m}${on}`] = knockLimitedCR(m, on);
      const r = kl[`${m}${on}`];
      console.log(`[knock] ${m} PRF ${on}: knock-limited CR ${r.cr.toFixed(2)} (guide ${guideCR(m, on).toFixed(2)}, ${(r.cr - guideCR(m, on)) >= 0 ? '+' : ''}${(r.cr - guideCR(m, on)).toFixed(2)}): ${fmt(r.s)}`);
    }
  }, 900_000);

  it('RON standard-knock state of PRF 90 at its knock-limited CR matches the measured CFR standard knock', () => {
    // Kolodziej & Wallner 2017 Fig. 9 (PRF98 at its guide CR): CA50 8.8°, knock point 12.0°;
    // Pal et al. 2018 Fig. 9: knock point 11.0 ± 1.0°.
    const s = kl.RON90.s;
    const k = targets.standardKnockState;
    expect(Number.isFinite(s.knockOnsetDeg)).toBe(true);
    expect(Math.abs(s.knockOnsetDeg - k.KW2017_PRF98.knockPointDeg)).toBeLessThan(4);
    expect(Math.abs(s.ca50 - k.KW2017_PRF98.ca50Deg)).toBeLessThan(4);
  });

  for (const [method, on, ok] of [
    ['RON', 80, true],
    ['RON', 90, true],
    ['RON', 100, true],
    ['MON', 100, true],
    ['RON', 60, false],
    ['RON', 70, false],
    ['MON', 70, true],
    ['MON', 90, true],
  ] as [Method, number, boolean][]) {
    const title = `${method} PRF ${on}: knock-limited CR within ±${CR_TOL} of the guide table`;
    const body = (): void => {
      expect(Math.abs(kl[`${method}${on}`].cr - guideCR(method, on))).toBeLessThan(CR_TOL);
    };
    // PRF ≤ 70 RON (−0.84 / −0.38 CR in round 1; −0.94 / −0.77 after the round-2 fixer pass):
    // model-form limitation (calibration.ts CFR_KNOCK_DELAY_MODEL). MON PRF 70 / 90 (+0.36 / +0.56 CR in
    // round 1) pass since the round-2 fixer pass (−0.16 / −0.14).
    if (ok) it(title, body);
    else it.fails(`${title} [known limitation]`, body);
  }

  it('octane scale: dCR/dON (PRF 80 → 100) within ±30 % of the guide table, RON and MON', () => {
    for (const m of ['RON', 'MON'] as Method[]) {
      const sm = (kl[`${m}100`].cr - kl[`${m}${m === 'RON' ? 80 : 70}`].cr) / (m === 'RON' ? 20 : 30);
      const sg = (guideCR(m, 100) - guideCR(m, m === 'RON' ? 80 : 70)) / (m === 'RON' ? 20 : 30);
      console.log(`[knock] ${m} slope ${sm.toFixed(4)} vs guide ${sg.toFixed(4)} CR/ON`);
      expect(Math.abs(sm / sg - 1)).toBeLessThan(0.3);
    }
  });

  // Round 1 (known limitation): the knock-limited CR moved only ≈ 0.015 CR per degree of spark (0.21
  // ON/°). Passes since the round-2 fixer pass (≈ 0.7 ON/°: all unburned gas autoignites, so the
  // autoigniting mass grows with advance). NB the target is KNOCKMETER-based (KU); with MAPO as the
  // metric Pal's own data give ≈ 0.26 bar/° against ≈ 10 bar per CR (≈ 0.37 ON/°).
  it('RON: one degree of spark advance is worth ≈ 1 ON of knock near standard knock (0.4–2.5 ON/°)', () => {
    // Pal et al. 2018 Fig. 1: +12.9 KU per degree advance around 13° bTDC; ASTM spread 12–15 KU/ON
    // → 0.86–1.07 ON/° (targets.sparkEquivalence). Evaluated on the knock-limited CR:
    // (∂CR_KL/∂θ_spark)/(∂CR_KL/∂ON).
    const a15 = knockLimitedCR('RON', 90, { sparkAdvanceDeg: 15 }).cr;
    const a11 = knockLimitedCR('RON', 90, { sparkAdvanceDeg: 11 }).cr;
    const o95 = knockLimitedCR('RON', 95).cr;
    const o85 = knockLimitedCR('RON', 85).cr;
    const perDeg = (a11 - a15) / 4;
    const perOn = (o95 - o85) / 10;
    const ratio = perDeg / perOn;
    console.log(`[knock] RON PRF 90: −dCR_KL/dθ_spark ${perDeg.toFixed(4)}/°, dCR_KL/dON ${perOn.toFixed(4)}/ON → ${ratio.toFixed(2)} ON per degree`);
    expect(ratio).toBeGreaterThan(0.4);
    expect(ratio).toBeLessThan(2.5);
  });
});
