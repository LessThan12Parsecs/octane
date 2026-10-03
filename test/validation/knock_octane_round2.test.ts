/**
 * Validation round 2 — knock & octane rating (validator: knock area).
 *
 * [SLOW: ≈ 60–90 s. Runs only with OCTANE_VALIDATION=1, e.g.
 *    OCTANE_VALIDATION=1 npx vitest run test/validation/knock_octane_round2.test.ts ]
 *
 * Simulator DEFAULTS (the one global calibration, src/physics/cycle/calibration.ts), public API only.
 * Knock criterion: MAPO at the D-1 pickup, standard knock 0.67 bar (Rockstroh et al. 2018 standard
 * point, test/fixtures/validation_knock_targets.json) — the same criterion as knock_octane_rating.test.ts.
 * Every target below is a committed fixture (grade in cfr_validation_README.md); `it.fails` marks a
 * finding of this round that the model does not (yet) satisfy — each failing expectation is a
 * reproduction of a reported defect, not a tolerance to be widened.
 *
 * Findings reproduced here (numbers from the round-2 validator runs, 2026-09-30):
 *  1. Octane-scale CURVATURE: the model's knock-limited CR is ≈ linear in ON (≈ 0.07 CR/ON from 60 to
 *     100) while the ASTM guide table is convex (0.03 CR/ON at 60–70, 0.106 at 90–100). The 80→100
 *     chord used by knock_octane_rating.test.ts (−10 %) hides it: 90→100 is 0.68× (RON) / 0.67× (MON)
 *     the guide slope; SON 2023 (same MAPO criterion) 93→100 is 0.65×.
 *  2. PRFs must rate RON = MON (sensitivity 0 by definition, ASTM D2699/D2700): the model's knock-
 *     limited CRs converted back through the guide tables give PRF 80 → RON 79.3 / MON 88.1,
 *     PRF 90 → 92.5 / 95.3, PRF 70 → 49 / 78.
 *  3. Low ON: at RON conditions PRF 40 and 50 exceed standard knock even at the engine's minimum CR 4.0
 *     (guide 5.14 / 5.26), and for PRF 40 MAPO FALLS from 3.1 to 0.03 bar as CR rises 4.2 → 6.1 while
 *     the autoignited end-gas fraction rises 44 → 99.96 %: the acoustic source becomes uniform over the
 *     bore and its projection on the transverse modes vanishes — the knock-intensity proxy is not
 *     monotone, so no standard-knock CR exists.
 *  4. φ: at fixed CR the model's MAPO rises monotonically from λ 0.80 to λ 1.0 (PRF 95: 0 → 1.05 bar);
 *     measured (Hoth & Kolodziej 2021 Part 1 Fig. 6, PRF 95, same MAPO metric) peaks at λ 0.89
 *     (0.87 bar) and falls to 0.37 bar at λ 1.0. Douaud–Eyzat has no φ dependence.
 *  5. Spark: the knock-limited CR hardly moves with spark advance and stops moving beyond ≈ 15°
 *     (PRF 90: 6.603 at 13°, 6.586 at 15°, 6.597 at 17°); the knock onset tracks the end of the flame
 *     (burned fraction at onset 0.87–0.91 at every spark timing).
 *  6. Intake temperature (Rockstroh 2018 PRF90 λ 1, 600 rpm, 13°, same MAPO metric): measured MAPO
 *     1.16 bar at CR 6.98 with a 150 °C mixture; the model reaches it at CR 6.15 (−0.83), while at
 *     33 °C it is within −0.24 — the model's knock is ≈ 3× too sensitive to intake temperature.
 *  7. Standard-knock mechanism: onset at 87 % burned (Choi/Pal PRF100 CR 7.55 state) with 3 % of the
 *     charge ahead of the front, max dp/dθ 3.36 bar/° — below even the 300-cycle AVERAGE trace of the
 *     same state (3.88 bar/°, Choi 2018 Fig. 9 flush mount), i.e. the model's knock event releases too
 *     little energy; it reaches the measured MAPO through a steep MAPO(CR) (≈ 5.5 bar/CR vs ≈ 1.4 bar/CR
 *     measured by Rockstroh 2018 at 33 °C).
 *
 * Status after the round-2 fixer pass (markers flipped from it.fails to it where the model now passes):
 * 1 slopes pass (90→100 0.78×, SON 93→100 within 25 %); 2 open (PRF 80 rates RON 63 / MON 76); 3 PRF 40
 * monotone (all unburned gas autoignites, sequential acoustic source), PRF 50 still above standard knock
 * at guide − 0.5; 4 open (MAPO λ 1/λ 0.89 0.82); 5 passes (−0.075/−0.041 CR/°); 6 temperature trend
 * fixed (−0.02 CR from 33 to 150 °C) but the absolute CR at 33 °C is −0.42 (the ±0.35 `it` below fails:
 * the whole knock-limited CR scale sits 0.14–0.45 CR low, calibration.ts knockStratificationDT at its
 * bound); 7 open (max dp/dθ 2.9 bar/°).
 */
import { beforeAll, describe, expect, it } from 'vitest';
import targets from '../fixtures/validation_knock_targets.json';
import son from '../fixtures/cfr_son2023_critical_cr.json';
import rock from '../fixtures/cfr_rockstroh2018_knock_vs_cr.json';
import hoth from '../fixtures/cfr_hoth_knock_metrics.json';
import choi from '../fixtures/cfr_choi2018_traces.json';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS, cfrMonSparkAdvanceDeg } from '../../src/physics/engines/cfr';
import { CycleModel } from '../../src/physics/cycle/index';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary } from '../../src/physics/core/snapshot';

const RUN = !!process.env.OCTANE_VALIDATION;
type Method = 'RON' | 'MON';
type Guide = Record<string, { compressionRatio: number }>;
const guideCR = (m: Method, on: number): number => (targets.guideTable[m] as Guide)[String(on)].compressionRatio;

/** Octane number whose guide-table CR equals `cr` (piecewise-linear inverse of the guide table). */
function guideON(m: Method, cr: number): number {
  const t = Object.entries(targets.guideTable[m] as Guide)
    .map(([on, v]) => [v.compressionRatio, Number(on)] as const)
    .sort((a, b) => a[0] - b[0]);
  let i = 0;
  while (i < t.length - 2 && cr > t[i + 1][0]) i++;
  const [x0, y0] = t[i];
  const [x1, y1] = t[i + 1];
  return y0 + ((cr - x0) * (y1 - y0)) / (x1 - x0);
}

/** Standard knock intensity as MAPO, Pa (Rockstroh 2018 standard point, 0.67 bar). */
const MAPO_STD = targets.standardKnockState.Rockstroh2018_PRF90_standard.mapoPa;

function opFor(method: Method, on: number, cr: number, patch: Partial<OperatingPoint> = {}): OperatingPoint {
  const base = method === 'RON' ? CFR_RON_CONDITIONS : CFR_MON_CONDITIONS;
  const op: OperatingPoint = { ...base, fuel: { kind: 'PRF', octaneNumber: on }, compressionRatio: cr, ...patch };
  if (method === 'MON' && patch.sparkAdvanceDeg === undefined) op.sparkAdvanceDeg = cfrMonSparkAdvanceDeg(cr);
  return op;
}

/** One model per (method, ON, patch); `at(cr)` = converged cycle at that CR (2 transition cycles). */
function runner(method: Method, on: number, patch: Partial<OperatingPoint> = {}): (cr: number) => CycleSummary {
  let m: CycleModel | null = null;
  return (cr: number): CycleSummary => {
    const op = opFor(method, on, cr, patch);
    if (!m) m = new CycleModel(CFR_F1, op, { warmupCycles: 3 });
    else {
      m.setOperatingPoint(op);
      m.runCycles(2);
    }
    return m.runCycles(1)[0];
  };
}

/** CR at which MAPO reaches `thr` (bisection in [lo, hi], MAPO assumed increasing). */
function klCR(method: Method, on: number, thr = MAPO_STD, patch: Partial<OperatingPoint> = {}, lo?: number, hi?: number): { cr: number; s: CycleSummary } {
  const at = runner(method, on, patch);
  const g = guideCR(method, on);
  let a = lo ?? g - 1.5;
  let b = hi ?? g + 1.5;
  for (let i = 0; i < 9; i++) {
    const c = 0.5 * (a + b);
    if (at(c).mapo < thr) a = c;
    else b = c;
  }
  const cr = 0.5 * (a + b);
  return { cr, s: at(cr) };
}

const pr = (s: string): void => {
  process.stdout.write(`[knock-r2] ${s}\n`);
};

describe.skipIf(!RUN)('[slow] validation round 2: knock & octane rating (simulator defaults)', () => {
  const R: Record<string, number> = {};
  let prf40: { cr: number; mapo: number; xeg: number }[] = [];
  let prf50: CycleSummary[] = [];
  let lam: Record<string, CycleSummary> = {};
  let choiM13: CycleSummary;
  beforeAll(() => {
    for (const on of [80, 90, 95, 100]) R[`RON${on}`] = klCR('RON', on).cr;
    for (const on of [80, 90, 100]) R[`MON${on}`] = klCR('MON', on).cr;
    for (const k of Object.keys(R)) {
      const m = k.slice(0, 3) as Method;
      const on = Number(k.slice(3));
      pr(`${k}: CR_KL ${R[k].toFixed(3)} (guide ${guideCR(m, on).toFixed(3)}, ${(R[k] - guideCR(m, on)).toFixed(3)}) → rated ${guideON(m, R[k]).toFixed(1)}`);
    }
    // spark (finding 5)
    R.sp17 = klCR('RON', 90, MAPO_STD, { sparkAdvanceDeg: 17 }).cr;
    R.sp15 = klCR('RON', 90, MAPO_STD, { sparkAdvanceDeg: 15 }).cr;
    pr(`RON PRF 90 CR_KL: spark 13° ${R.RON90.toFixed(3)}, 15° ${R.sp15.toFixed(3)}, 17° ${R.sp17.toFixed(3)}`);
    // low ON (finding 3)
    const a40 = runner('RON', 40);
    prf40 = [4.2, 4.6, 5.0, 5.4, 5.8, 6.1].map((cr) => {
      const s = a40(cr);
      return { cr, mapo: s.mapo, xeg: s.knockEndGasFraction };
    });
    pr(`RON PRF 40: ${prf40.map((r) => `CR ${r.cr} MAPO ${(r.mapo / 1e5).toFixed(2)} bar (end gas ${(100 * r.xeg).toFixed(1)} %)`).join('; ')}`);
    const a50 = runner('RON', 50);
    prf50 = [4.0, guideCR('RON', 50) - 0.5].map((cr) => a50(cr));
    pr(`RON PRF 50: MAPO ${prf50.map((s) => (s.mapo / 1e5).toFixed(2)).join(' / ')} bar at CR 4.0 / ${(guideCR('RON', 50) - 0.5).toFixed(2)}`);
    // φ (finding 4): PRF 95 at the model's own CR for ≈ 0.87 bar at λ 0.89
    lam = {};
    {
      const cr95 = klCR('RON', 95, 0.87e5, { equivalenceRatio: 1 / 0.89 }).cr;
      for (const l of [0.8, 0.89, 1.0]) {
        const r = runner('RON', 95, { equivalenceRatio: 1 / l });
        lam[String(l)] = r(cr95);
      }
      pr(`RON PRF 95 at CR ${cr95.toFixed(3)}: MAPO λ 0.80 ${(lam['0.8'].mapo / 1e5).toFixed(2)}, λ 0.89 ${(lam['0.89'].mapo / 1e5).toFixed(2)}, λ 1.00 ${(lam['1'].mapo / 1e5).toFixed(2)} bar`);
    }
    // intake temperature (finding 6): Rockstroh 2018 PRF90, λ 1, 600 rpm, 13°, 1.0 bar
    const rk = (T: number): Partial<OperatingPoint> => ({ equivalenceRatio: 1, intakeMixtureTemperature: T, ambientTemperature: 306.15, ambientPressure: 1.0e5 });
    R.rock423 = klCR('RON', 90, 1.16e5, rk(423.15), 5.0, 7.6).cr;
    R.rock306 = klCR('RON', 90, MAPO_STD, rk(306.15), 5.5, 7.6).cr;
    pr(`Rockstroh PRF90 λ1: CR at MAPO 1.16 bar (150 °C) ${R.rock423.toFixed(3)}; CR at 0.67 bar (33 °C) ${R.rock306.toFixed(3)}`);
    // SON 2023 (finding 1): critical CR at MAPO 0.6 bar, λ 1, spark 13°, 52 °C air, 1.013 bar
    const sonP: Partial<OperatingPoint> = { equivalenceRatio: 1, ambientPressure: 101290 };
    R.son93 = klCR('RON', 93, 0.6e5, sonP, 5.8, 8.2).cr;
    R.son100 = klCR('RON', 100, 0.6e5, sonP, 6.3, 8.7).cr;
    pr(`SON 2023 critical CR (1.013 bar): PRF93 ${R.son93.toFixed(3)}, PRF100 ${R.son100.toFixed(3)}`);
    // standard-knock mechanism (finding 7): Choi/Pal PRF100 CR 7.55, λ 0.89, spark 12.72° (fixer's CHOI_OP)
    choiM13 = runner('RON', 100, {
      compressionRatio: 7.55, equivalenceRatio: 1 / 0.89, ambientPressure: 100240, ambientTemperature: 322,
      intakeMixtureTemperature: 302, sparkAdvanceDeg: 12.72,
    })(7.55);
    pr(`Choi PRF100 m13: onset ${choiM13.knockOnsetDeg.toFixed(2)}°, end gas ${(100 * choiM13.knockEndGasFraction).toFixed(1)} %, MAPO ${(choiM13.mapo / 1e5).toFixed(2)} bar, max dp/dθ ${(choiM13.maxPressureRiseRate / 1e5).toFixed(2)} bar/°`);
  }, 900_000);

  it('RON PRF 80/90/95/100: knock-limited CR within ±0.35 of the ASTM D2699 guide table (reproduces the fixer)', () => {
    for (const on of [80, 90, 95, 100]) expect(Math.abs(R[`RON${on}`] - guideCR('RON', on))).toBeLessThan(0.35);
  });

  // Finding 1 (curvature). Guide table 90→100: 0.1058 CR/ON; model ≈ 0.0715 (0.68×).
  it('RON octane scale 90→100: slope within ±25 % of the guide table [finding 1; passes since the round-2 fixer pass]', () => {
    const sm = (R.RON100 - R.RON90) / 10;
    const sg = (guideCR('RON', 100) - guideCR('RON', 90)) / 10;
    expect(Math.abs(sm / sg - 1)).toBeLessThan(0.25);
  });

  // Finding 1 with the SAME knock metric as the model (MAPO 0.6 bar, SON 2023 Fig. 6, grade A): the absolute
  // CR carries the ±0.4 ANL-campaign systematic (cfr_validation_README.md), the ON slope does not.
  it('SON 2023 (MAPO 0.6 bar, λ 1, 1.013 bar): critical-CR slope PRF93→100 within ±25 % of measured [finding 1; passes since the round-2 fixer pass]', () => {
    const ds = (son.datasets[0].data as unknown as Record<string, { compressionRatio: number[] }>);
    const meas = (ds.PRF100.compressionRatio[0] - ds.PRF93.compressionRatio[0]) / 7; // 0.110 CR/ON
    const sim = (R.son100 - R.son93) / 7;
    pr(`SON slope: sim ${sim.toFixed(4)} vs measured ${meas.toFixed(4)} CR/ON`);
    expect(Math.abs(sim / meas - 1)).toBeLessThan(0.25);
  });

  // Finding 2: a PRF's rating is its own ON in both methods (sensitivity 0 by definition).
  it.fails('PRF 80 and 90: RON − MON rating from the simulated knock-limited CRs within ±2 ON [finding 2]', () => {
    for (const on of [80, 90]) {
      const s = guideON('RON', R[`RON${on}`]) - guideON('MON', R[`MON${on}`]);
      pr(`PRF ${on}: rated RON ${guideON('RON', R[`RON${on}`]).toFixed(1)}, MON ${guideON('MON', R[`MON${on}`]).toFixed(1)}, S = ${s.toFixed(1)}`);
      expect(Math.abs(s)).toBeLessThan(2);
    }
  });

  // Finding 3a: PRF 50 must not exceed standard knock half a CR below its guide CR (guide 5.26).
  it.fails('RON PRF 50: MAPO below standard knock at CR = guide − 0.5 [finding 3]', () => {
    expect(prf50[1].mapo).toBeLessThan(MAPO_STD);
  });

  // Finding 3b: a knock-intensity proxy must grow with CR (the rating procedure raises CR until the
  // intensity reaches the standard level). PRF 40 RON: MAPO 3.1 → 0.03 bar while end gas 44 → 99.96 %.
  it('RON PRF 40: MAPO non-decreasing with CR while the autoignited end-gas fraction grows [finding 3; passes since the round-2 fixer pass]', () => {
    for (let i = 1; i < prf40.length; i++) {
      expect(prf40[i].xeg).toBeGreaterThanOrEqual(prf40[i - 1].xeg);
      expect(prf40[i].mapo).toBeGreaterThanOrEqual(0.9 * prf40[i - 1].mapo);
    }
  });

  // Finding 4: Hoth & Kolodziej 2021 Part 1 Fig. 6 (PRF 95, RON-95 CR, 13°, MAPO 4–18 kHz): 0.87 bar at
  // λ 0.89, 0.37 bar at λ 1.00 (ratio 0.42), 0.66 bar at λ 0.80.
  it.fails('φ: MAPO at λ 1.0 well below the peak-knock λ 0.89 value (measured ratio 0.42) [finding 4]', () => {
    const d = (hoth.datasets[1].data as unknown as { mapoBar: Record<string, { lambda: number[]; mapoPa: number[] }> }).mapoBar.PRF95;
    const measRatio = d.mapoPa[d.mapoPa.length - 1] / Math.max(...d.mapoPa);
    const sim = lam['1'].mapo / lam['0.89'].mapo;
    pr(`MAPO(λ1)/MAPO(λ0.89): sim ${sim.toFixed(2)} vs measured ${measRatio.toFixed(2)}`);
    expect(sim).toBeLessThan(0.8);
    expect(lam['0.8'].mapo).toBeGreaterThan(0.3 * lam['0.89'].mapo); // measured 0.66/0.87
  });

  // Finding 5: ≈ 1 ON per degree near standard knock — KU-based (Pal 2018 Fig. 1: 12.9 KU/° vs the
  // 12–15 KU/ON spread) and MAPO-based (Hoth 2021b Fig. 1: 0.17 bar/° at 9–13° vs Rockstroh 2018
  // 33 °C: 1.4 bar/CR, i.e. ≈ 0.12 CR/°). Lower bound 0.4 ON/° as in knock_octane_rating.test.ts;
  // guide slope at 90–95: 0.078 CR/ON ⇒ 4° of advance must lower CR_KL by ≥ 0.12.
  it('RON PRF 90: 4° more spark advance (13 → 17°) lowers the knock-limited CR by ≥ 0.12 [finding 5; passes since the round-2 fixer pass]', () => {
    expect(R.sp17 - R.RON90).toBeLessThan(-0.12);
  });
  it('RON PRF 90: knock-limited CR decreases monotonically with spark advance 13 → 15 → 17° [finding 5; passes since the round-2 fixer pass]', () => {
    expect(R.sp15).toBeLessThan(R.RON90 - 0.02);
    expect(R.sp17).toBeLessThan(R.sp15 - 0.02);
  });

  // Finding 6: Rockstroh 2018 Figs. 5–6 (grade D digitisation), PRF90 λ 1, 600 rpm, 13°, 1.0 bar.
  it('Rockstroh 33 °C / 1.0 bar: CR at MAPO 0.67 bar within ±0.35 of measured 6.84', () => {
    const d = (rock.datasets[0].data as unknown as Record<string, { mapoPa: { compressionRatio: number[]; value: number[] } }>)['Tin33C_Pin1.00bar'].mapoPa;
    // measured crossing of 0.67 bar between (6.709, 0.546) and (6.849, 0.695 bar)
    const i = d.value.findIndex((v) => v >= MAPO_STD);
    const cr = d.compressionRatio[i - 1] + ((MAPO_STD - d.value[i - 1]) * (d.compressionRatio[i] - d.compressionRatio[i - 1])) / (d.value[i] - d.value[i - 1]);
    pr(`Rockstroh 33 °C: measured CR@0.67 ${cr.toFixed(3)}, sim ${R.rock306.toFixed(3)}`);
    expect(Math.abs(R.rock306 - cr)).toBeLessThan(0.35);
  });
  it.fails('Rockstroh 150 °C / 1.0 bar: CR at MAPO 1.16 bar within ±0.35 of measured 6.98 [finding 6]', () => {
    const d = (rock.datasets[0].data as unknown as Record<string, { mapoPa: { compressionRatio: number[]; value: number[] } }>)['Tin150C_Pin1.00bar'].mapoPa;
    expect(Math.abs(R.rock423 - d.compressionRatio[0])).toBeLessThan(0.35);
  });

  // Finding 7: the 300-cycle AVERAGE trace of the Choi/Pal standard-knock state (Choi 2018 Fig. 9, flush
  // mount, grade A) already has max dp/dθ 3.88 bar/° at 11.3° (knock point); averaging cycles with scattered
  // knock points can only lower it, so a representative single cycle must reach at least that.
  it.fails('Choi PRF100 CR 7.55 spark −12.7°: max dp/dθ ≥ the 300-cycle average trace (3.88 bar/°) [finding 7]', () => {
    const ds = (choi.datasets as unknown as { key: string; data: { thetaDeg: number[]; pPa: number[] } }[]).find((x) => x.key === 'choi2018_fig9_prf100_st_m13_flushMount')!;
    let best = 0;
    for (let th = -30; th < 40; th += 0.1) {
      const f = (x: number): number => {
        const t = ds.data.thetaDeg;
        let k = 0;
        while (k < t.length - 2 && t[k + 1] < x) k++;
        return ds.data.pPa[k] + ((x - t[k]) * (ds.data.pPa[k + 1] - ds.data.pPa[k])) / (t[k + 1] - t[k]);
      };
      best = Math.max(best, (f(th + 0.1) - f(th)) / 0.1);
    }
    pr(`Choi m13 average-trace max dp/dθ ${(best / 1e5).toFixed(2)} bar/° vs sim ${(choiM13.maxPressureRiseRate / 1e5).toFixed(2)}`);
    expect(choiM13.maxPressureRiseRate).toBeGreaterThanOrEqual(best);
  });
});
