/**
 * Validation round 2 — agreement with MEASURED CFR data (validator: measured-data area).
 *
 * [SLOW: ≈ 30–60 s. Runs only with OCTANE_VALIDATION=1, e.g.
 *    OCTANE_VALIDATION=1 npx vitest run test/validation/measured_cfr_round2.test.ts ]
 *
 * Simulator DEFAULTS (the one global calibration, src/physics/cycle/calibration.ts), public API only
 * (CycleModel). Complements measured_cfr_data.test.ts (round 1: burn angles, compression, fuel flow,
 * ASTM compression pressure) with the datasets and derived quantities that round did not use. Every
 * target is a committed fixture (grades in test/fixtures/cfr_validation_README.md). `it.fails` marks
 * a finding of this round that the model does not satisfy — flip to `it` when a fix makes it pass.
 *
 * "Apparent" heat release = single-zone net heat release with γ = 1.30 (as measured_cfr_data.test.ts),
 * evaluated with the SAME routine on the measured and the simulated pressure.
 *
 * Findings reproduced (numbers from the round-2 runs, 2026-09-30):
 *  A. Heat-loss partition: after the apparent-heat-release peak the model loses 2.1× (Choi ST +5) /
 *     1.7× (ST −13) the measured apparent energy up to 130° (Woschni ×1.3); Woschni ×0.8 reproduces
 *     the measured decline (63 vs 63 J, 152 vs 161 J) but then the apparent-heat-release PEAK is
 *     +11 % and work/fuel +12 %. At ×1.3 the peak is still +5…+7 % with the fuel flow matched
 *     (Choi PRF98 +1.6 %). The Woschni multiplier absorbs a phasing-independent sink the model lacks
 *     (crevice storage / unburned HC / blow-by, cf. no crevice model in DESIGN.md).
 *  B. MBT: the model's MBT at the Choi PRF100 state is 5° BTDC with CA50 = 16.8° aTDC; gIMEP at
 *     13° BTDC is 3.7 % BELOW that at 5° BTDC. Measured (Hoth & Kolodziej 2021b Fig. 3, PKL, RON-95
 *     height): gIMEP at −12.8° is 0 … +1.4 % ABOVE that at −4.7°, and Choi 2018 Fig. 9 gIMEP(ST −13) is
 *     +3.1 % above gIMEP(ST +5) (model +1.5 %). Follows from A (early phasing over-penalised).
 *  C. λ = 1 burn duration (Hoth 2025 Table 7, PRF 93–97, λ 1, CA50 = 12°): measured spark 13.2–13.5°
 *     BTDC; model 10.1° (burn ≈ 3° too short at λ 1), while at λ 0.89 it is ≈ 1–2° short (KW17 PRF98
 *     CA50 8.78°, Choi 2018 Fig. 2 apparent 8.48°; model 7.0° / 6.6°). With u_T = C_T·u′ and C_T = 5.25
 *     the entrainment speed is ≈ 88 % turbulence: the burn is too insensitive to S_L (λ, T_u).
 *  D. Blowdown: p(160°)/p(141°) = 0.761 (model) vs 0.814 measured (Choi ST +5; ST −13 0.776 vs 0.816):
 *     the exhaust valve's early effective area / the 0-D 10 L exhaust plenum (measured exhaust-port
 *     pulse 1.09 bar at 157°, model 1.01) empty the cylinder too fast.
 *  E. Intake port (Choi 2018 Fig. 2, PRF98, Kulite in the port): the 1 L plenum (UNVERIFIED volume)
 *     gives a dip to 0.882 bar at −250° vs 0.838 bar at −266° and recovers too slowly (0.958 vs 0.997
 *     bar at BDC): RMS 0.025 bar. A 0.25 L volume gives RMS 0.0145 bar (dip 0.854 at −269°, BDC 0.984)
 *     and keeps the ASTM venturi-size effect (7.8 % vs 8.0 %) and the fuel flow (+2.8 / +0.6 %).
 *  F. HCCI compression (Kalvakala et al., PRF90 λ 3, spark off): at the STATED CR 12.84 / 14.05 the
 *     model's compression pressure is +8 … +11 % high; it matches (±2 %) at CR 11.35 / 12.27, i.e. the
 *     stated CRs converted with a constant clearance-volume offset ΔV_c = 7.38 cm³ — exactly the offset
 *     implied by Hoth 2025 (RON-95 height = CR 7.26 on the ANL 2020s scale vs 6.82 by
 *     cfrCompressionRatioAtCounter). The ANL 2020s CRs (SON 2023, Hoth 2025, Kalvakala) must be
 *     converted before comparison; with it SON PRF100 at 1.013 bar is 7.30 vs model 7.29.
 *  G. Knock-limit pressure sensitivity (SON 2023, PRF100, λ 1, 13°, MAPO 0.6 bar): model −1.37 CR from
 *     1.013 to 1.5 bar vs measured −1.61 (−15 %; −1.40 after the ΔV_c conversion): passes.
 *  H. ASTM motored compression pressure (D2699 §10.3.17, A2.2.3, A2.3: warm up fired, SHUT DOWN, drain
 *     the carburettor, fit the gauge, restart MOTORED): the model applies the FIRED standard-knock wall
 *     temperatures of Pal 2018 Table 3 (head 525 K, piston 485 K) to this motored check and reads
 *     +4.5 / +4.7 / +4.6 % (RON 930 / 778 / 1061) and +5.4 % (MON 930). With motored walls near the
 *     100 °C jacket (385–400 K, UNVERIFIED estimate) it reads +2.2 … +2.6 %: half of the discrepancy the
 *     round-1 fix attributed to the counter→CR relation is the wall temperature. Walls do not depend
 *     on the operating point in the model (fired/motored, RON/MON, spark timing).
 */
import { beforeAll, describe, expect, it } from 'vitest';
import choi from '../fixtures/cfr_choi2018_traces.json';
import hcci from '../fixtures/cfr_hcci_autoignition.json';
import son from '../fixtures/cfr_son2023_critical_cr.json';
import { CycleModel } from '../../src/physics/cycle/index';
import type { CycleTrace } from '../../src/physics/cycle/index';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS, cfrCompressionRatioAtCounter } from '../../src/physics/engines/cfr';
import type { EngineSpec } from '../../src/physics/core/engine-spec';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary } from '../../src/physics/core/snapshot';

const RUN = !!process.env.OCTANE_VALIDATION;
const B = 0.08255;
const S = 0.1143;
const L_ROD = 0.254;
const VD = (Math.PI / 4) * B * B * S;
/** Clearance-volume offset between the ANL 2020s CR scale and cfrCompressionRatioAtCounter, m³:
 * Hoth 2025 (OSTI 2561394) RON-95 height = CR 7.26 vs 6.819 here (cfr_hoth_knock_metrics.json). */
const DVC_ANL2020S = VD / (cfrCompressionRatioAtCounter(805) - 1) - VD / (7.26 - 1);
const fromAnl2020s = (cr: number): number => 1 + VD / (VD / (cr - 1) + DVC_ANL2020S);

// ── numerics shared by measured and simulated traces (as measured_cfr_data.test.ts) ────────
function volume(thDeg: number, cr: number): number {
  const t = (thDeg * Math.PI) / 180;
  const a = S / 2;
  const x = a + L_ROD - (a * Math.cos(t) + Math.sqrt(L_ROD * L_ROD - (a * Math.sin(t)) ** 2));
  return VD / (cr - 1) + (Math.PI / 4) * B * B * x;
}
function interp(xs: ArrayLike<number>, ys: ArrayLike<number>, x: number): number {
  let lo = 0;
  let hi = xs.length - 1;
  if (x <= xs[0]) return ys[0];
  if (x >= xs[hi]) return ys[hi];
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid] <= x) lo = mid;
    else hi = mid;
  }
  return ys[lo] + ((x - xs[lo]) / (xs[hi] - xs[lo])) * (ys[hi] - ys[lo]);
}
function monotone(th: ArrayLike<number>, y: ArrayLike<number>): { th: number[]; y: number[] } {
  const o: number[] = [];
  const oy: number[] = [];
  for (let i = 0; i < th.length; i++) {
    if (o.length && th[i] <= o[o.length - 1]) {
      if (th[i] === o[o.length - 1]) oy[oy.length - 1] = y[i];
      continue;
    }
    o.push(th[i]);
    oy.push(y[i]);
  }
  return { th: o, y: oy };
}
/** Cumulative apparent net heat release (γ = 1.30, 1° box smoothing) from the spark, J, on a 0.2° grid. */
function apparentQ(th: number[], p: number[], cr: number, sparkDeg: number): { g: number[]; Q: number[] } {
  const d = 0.2;
  const g: number[] = [];
  for (let x = sparkDeg - 10; x < 140; x += d) g.push(x);
  const P0 = g.map((x) => interp(th, p, x));
  const P = P0.map((_, i) => {
    let s = 0;
    let n = 0;
    for (let k = i - 2; k <= i + 2; k++) if (k >= 0 && k < P0.length) (s += P0[k]), n++;
    return s / n;
  });
  const V = g.map((x) => volume(x, cr));
  const gam = 1.3;
  const dQ = g.map((_, i) => {
    const i0 = Math.max(0, i - 1);
    const i1 = Math.min(g.length - 1, i + 1);
    return ((gam / (gam - 1)) * P[i] * (V[i1] - V[i0]) + (1 / (gam - 1)) * V[i] * (P[i1] - P[i0])) / (g[i1] - g[i0]);
  });
  const Q: number[] = [0];
  for (let i = 1; i < g.length; i++) Q.push(Q[i - 1] + (g[i] >= sparkDeg ? 0.5 * (dQ[i] + dQ[i - 1]) * d : 0));
  return { g, Q };
}
function peakAndDecline(a: { g: number[]; Q: number[] }, end = 130): { peak: number; decline: number } {
  let iMax = 0;
  for (let i = 1; i < a.Q.length; i++) if (a.Q[i] > a.Q[iMax]) iMax = i;
  return { peak: a.Q[iMax], decline: a.Q[iMax] - interp(a.g, a.Q, end) };
}

// ── measured data ───────────────────────────────────────────────────────────────────────────
type Ds = { key: string; data: Record<string, unknown> };
const trace = (fixture: { datasets: unknown }, key: string, sub?: string): { th: number[]; p: number[] } => {
  const ds = (fixture.datasets as Ds[]).find((d) => d.key === key)!;
  const d = (sub ? ds.data[sub] : ds.data) as { thetaDeg: number[]; pPa: number[] };
  return { th: d.thetaDeg, p: d.pPa };
};
const CR_CHOI = 7.55;
const CHOI_OP: Partial<OperatingPoint> = {
  compressionRatio: CR_CHOI, equivalenceRatio: 1 / 0.89, fuel: { kind: 'PRF', octaneNumber: 100 },
  ambientPressure: 100240, ambientTemperature: 322, intakeMixtureTemperature: 302,
};
const CHOI98_OP: Partial<OperatingPoint> = {
  compressionRatio: cfrCompressionRatioAtCounter(873), equivalenceRatio: 1 / 0.9047, fuel: { kind: 'PRF', octaneNumber: 98 },
  ambientPressure: 99800, ambientTemperature: 322.9, intakeMixtureTemperature: 301.4, sparkAdvanceDeg: 13,
};
const KW98_OP: Partial<OperatingPoint> = {
  compressionRatio: cfrCompressionRatioAtCounter(904), equivalenceRatio: 1 / 0.8855, fuel: { kind: 'PRF', octaneNumber: 98 },
  ambientPressure: 96850, ambientTemperature: 314.5, intakeMixtureTemperature: 301.4, sparkAdvanceDeg: 13,
};

interface Run {
  s: CycleSummary;
  tr: CycleTrace;
  m: CycleModel;
  th: number[];
  p: number[];
}
function run(patch: Partial<OperatingPoint>, extra: Record<string, unknown> = {}, cycles = 4): Run {
  const op = { ...CFR_RON_CONDITIONS, fuel: { ...CFR_RON_CONDITIONS.fuel }, ...patch } as OperatingPoint;
  const m = new CycleModel(CFR_F1, op, { warmupCycles: 3, ...extra } as never);
  m.runCycles(cycles - 1);
  m.recordTrace();
  const s = m.runCycles(1)[0];
  const mm = monotone(m.trace!.theta, m.trace!.pressure);
  return { s, tr: m.trace!, m, th: mm.th, p: mm.y };
}

describe.skipIf(!RUN)('measured CFR data, round 2 [SLOW]', () => {
  const r: Record<string, Run> = {};
  beforeAll(() => {
    r.p5 = run({ ...CHOI_OP, sparkAdvanceDeg: -5 });
    r.m13 = run({ ...CHOI_OP, sparkAdvanceDeg: 12.72 }, { knock: false });
    r.choi98 = run(CHOI98_OP, {}, 6);
    r.kw98 = run(KW98_OP, {}, 5);
  }, 600_000);

  describe('A. heat-release / heat-loss partition (Choi 2018 Figs. 2, 9, grade A/B)', () => {
    const cases = [
      ['p5', 'choi2018_fig9_prf100_st_p5_flushMount', undefined, -5, CR_CHOI],
      ['m13', 'choi2018_fig9_prf100_st_m13_flushMount', undefined, 12.72, CR_CHOI],
      ['choi98', 'choi2018_fig2_prf98_ron', 'cylinder', 13, cfrCompressionRatioAtCounter(873)],
    ] as const;
    // fixed in the fixer round-2 pass (crevice zone + recalibrated heat transfer): ×1.33 / ×1.00
    it('A1: post-peak apparent energy decline (peak → 130°) within 0.7–1.4× measured [round 1: 2.1× / 1.7×]', () => {
      for (const [name, key, sub, spark, cr] of cases.slice(0, 2)) {
        const meas = trace(choi, key, sub);
        const a = peakAndDecline(apparentQ(meas.th, meas.p, cr, -spark));
        const b = peakAndDecline(apparentQ(r[name].th, r[name].p, cr, -spark));
        console.log(`[r2] ${name}: apparent decline peak→130° model ${b.decline.toFixed(0)} J vs measured ${a.decline.toFixed(0)} J (×${(b.decline / a.decline).toFixed(2)})`);
        expect(b.decline / a.decline).toBeGreaterThan(0.7);
        expect(b.decline / a.decline).toBeLessThan(1.4);
      }
    });
    it.fails('A2: apparent-heat-release peak within ±3 % of measured [model +4.6 … +7.4 % with the fuel flow matched]', () => {
      for (const [name, key, sub, spark, cr] of cases) {
        const meas = trace(choi, key, sub);
        const a = peakAndDecline(apparentQ(meas.th, meas.p, cr, -spark));
        const b = peakAndDecline(apparentQ(r[name].th, r[name].p, cr, -spark));
        console.log(`[r2] ${name}: apparent-HR peak model/measured ${(b.peak / a.peak).toFixed(3)}`);
        expect(Math.abs(b.peak / a.peak - 1)).toBeLessThan(0.03);
      }
    });
  });

  it.fails('B. MBT: gIMEP(13° BTDC) ≥ 0.995·gIMEP(5° BTDC) at the Choi PRF100 state (Hoth 2021b Fig. 3 PKL: +0…+1.4 %; Choi Fig. 9 ST −13 vs +5: +3.1 %) [model −3.7 %]', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, ...CHOI_OP, sparkAdvanceDeg: 13 } as OperatingPoint, { warmupCycles: 3, knock: false });
    const g = (sp: number): number => {
      m.setOperatingPoint({ sparkAdvanceDeg: sp });
      m.runCycles(2);
      return m.runCycles(1)[0].imepGross;
    };
    const g13 = g(13);
    const g5 = g(5);
    const gm5 = g(-5);
    console.log(`[r2] gIMEP 13°/5°/−5° BTDC ${(g13 / 1e5).toFixed(3)}/${(g5 / 1e5).toFixed(3)}/${(gm5 / 1e5).toFixed(3)} bar; 13 vs 5: ${(100 * (g13 / g5 - 1)).toFixed(1)} %`);
    expect(g13 / g5).toBeGreaterThan(0.995);
  });

  it.fails('C. burn phasing: PRF98 standard RON CA50 within 1.5° (KW17 Fig. 9 8.78°; Choi Fig. 2) and λ = 1 spark for CA50 = 12° within 1.5° of 13.2° BTDC (Hoth 2025 Table 7, PRF95) [model 7.0°; 10.1°]', () => {
    console.log(`[r2] KW17 PRF98 CA50 model ${r.kw98.s.ca50.toFixed(2)}° vs 8.78°; Choi PRF98 model ${r.choi98.s.ca50.toFixed(2)}°`);
    const op = { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 95 }, equivalenceRatio: 1.0, ambientPressure: 100000, compressionRatio: fromAnl2020s(7.34), sparkAdvanceDeg: 13 } as OperatingPoint;
    const m = new CycleModel(CFR_F1, op, { warmupCycles: 3, knock: false });
    let sp = 13;
    let s = m.runCycles(1)[0];
    for (let k = 0; k < 5 && Math.abs(s.ca50 - 12) > 0.1; k++) {
      sp += 0.9 * (s.ca50 - 12);
      m.setOperatingPoint({ sparkAdvanceDeg: sp });
      m.runCycles(2);
      s = m.runCycles(1)[0];
    }
    console.log(`[r2] Hoth 2025 λ 1 PRF95 CR ${op.compressionRatio.toFixed(2)}: spark for CA50 = 12° model ${sp.toFixed(2)}° vs 13.2° BTDC`);
    expect(Math.abs(r.kw98.s.ca50 - 8.78)).toBeLessThan(1.5);
    expect(Math.abs(sp - 13.2)).toBeLessThan(1.5);
  });

  it.fails('D. blowdown: p(160°)/p(141°) within 3 % of measured (Choi Fig. 9 flush: 0.814 / 0.816) [model 0.761 / 0.776]', () => {
    for (const [name, key] of [['p5', 'choi2018_fig9_prf100_st_p5_flushMount'], ['m13', 'choi2018_fig9_prf100_st_m13_flushMount']] as const) {
      const meas = trace(choi, key);
      const a = interp(meas.th, meas.p, 160) / interp(meas.th, meas.p, 141);
      const b = interp(r[name].th, r[name].p, 160) / interp(r[name].th, r[name].p, 141);
      console.log(`[r2] ${name}: p(160)/p(141) model ${b.toFixed(3)} vs measured ${a.toFixed(3)}`);
      expect(Math.abs(b / a - 1)).toBeLessThan(0.03);
    }
  });

  // fixed in the fixer round-2 pass (0.25 L plenum, cfr.ts): RMS 0.014 bar, dip 0.839 vs 0.839 bar
  it('E. intake-port pressure (Choi Fig. 2, PRF98): model intake-plenum pressure RMS < 0.02 bar and dip depth within 0.02 bar [round 1, 1 L plenum: 0.025 bar, 0.882 vs 0.838 bar]', () => {
    const meas = trace(choi, 'choi2018_fig2_prf98_ron', 'intakePort');
    const sim = monotone(r.choi98.tr.theta, r.choi98.tr.pIntake);
    let se = 0;
    let n = 0;
    let minS = Infinity;
    let minM = Infinity;
    for (let x = -360; x < 360; x += 1) {
      const ps = interp(sim.th, sim.y, x);
      const pm = interp(meas.th, meas.p, x);
      se += (ps - pm) ** 2;
      n++;
      minS = Math.min(minS, ps);
      minM = Math.min(minM, pm);
    }
    const rms = Math.sqrt(se / n);
    console.log(`[r2] intake port RMS ${(rms / 1e5).toFixed(4)} bar, min model ${(minS / 1e5).toFixed(3)} vs ${(minM / 1e5).toFixed(3)} bar`);
    expect(rms).toBeLessThan(0.02e5);
    expect(Math.abs(minS - minM)).toBeLessThan(0.02e5);
  });

  it('F. HCCI motored compression (Kalvakala BRON/BMON) matches within 3 % at the CR converted from the ANL 2020s scale (ΔV_c = Hoth 2025 offset); +8…+11 % at the stated CR', () => {
    for (const [key, cr, pin, T] of [['kalvakala_fig4_prf90_bron', 12.84, 130000, 325.15], ['kalvakala_fig5_prf90_bmon', 14.05, 100000, 422.15]] as const) {
      const meas = trace(hcci, key);
      const patch: Partial<OperatingPoint> = { equivalenceRatio: 1 / 3, fuel: { kind: 'PRF', octaneNumber: 90 }, ambientPressure: pin, ambientTemperature: Math.min(T, 325.15), intakeMixtureTemperature: T, relativeHumidity: 0.02 };
      const conv = run({ ...patch, compressionRatio: fromAnl2020s(cr) }, { combustionModel: 'none', knock: false }, 3);
      const stated = run({ ...patch, compressionRatio: cr }, { combustionModel: 'none', knock: false }, 3);
      const errs = [-40, -35, -30, -25].map((x) => interp(conv.th, conv.p, x) / interp(meas.th, meas.p, x) - 1);
      const errStated = interp(stated.th, stated.p, -30) / interp(meas.th, meas.p, -30) - 1;
      console.log(`[r2] ${key}: CR ${cr} → ${fromAnl2020s(cr).toFixed(2)}: p(−40…−25) errors ${errs.map((e) => (100 * e).toFixed(1)).join('/')} %; at the stated CR p(−30) ${(100 * errStated).toFixed(1)} %`);
      for (const e of errs) expect(Math.abs(e)).toBeLessThan(0.03);
      expect(errStated).toBeGreaterThan(0.05);
    }
  });

  it('G. knock-limit pressure sensitivity (SON 2023 Fig. 6, PRF100, λ 1, 13° BTDC, MAPO 0.6 bar): ΔCR(1.013 → 1.5 bar) within ±25 % of measured', () => {
    const ds = (son.datasets as unknown as { data: Record<string, { intakePressurePa: number[]; compressionRatio: number[] }> }[])[0].data.PRF100;
    const measLo = fromAnl2020s(ds.compressionRatio[0]);
    const measHi = fromAnl2020s(ds.compressionRatio[ds.compressionRatio.length - 1]);
    const klcr = (p: number): number => {
      let lo = 4.5;
      let hi = 9.5;
      const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 100 }, equivalenceRatio: 1, ambientPressure: p, sparkAdvanceDeg: 13, compressionRatio: 7 } as OperatingPoint, { warmupCycles: 3 });
      for (let i = 0; i < 9; i++) {
        const cr = 0.5 * (lo + hi);
        m.setOperatingPoint({ compressionRatio: cr });
        m.runCycles(2);
        if (m.runCycles(1)[0].mapo < 0.6e5) lo = cr;
        else hi = cr;
      }
      return 0.5 * (lo + hi);
    };
    const a = klcr(ds.intakePressurePa[0]);
    const b = klcr(ds.intakePressurePa[ds.intakePressurePa.length - 1]);
    console.log(`[r2] SON PRF100 critical CR 1.013 → 1.5 bar: model ${a.toFixed(2)} → ${b.toFixed(2)} (Δ ${(b - a).toFixed(2)}); measured (converted) ${measLo.toFixed(2)} → ${measHi.toFixed(2)} (Δ ${(measHi - measLo).toFixed(2)})`);
    expect(Math.abs((b - a) / (measHi - measLo) - 1)).toBeLessThan(0.25);
    expect(Math.abs(a - measLo)).toBeLessThan(0.35);
  });

  // Fixer round 2: the default model now computes the wall temperatures from the operating point
  // (lumped walls, options.wallTemperatureModel), so this motored check gets motored walls by default.
  // The round-1 comparison ("fired" = CFR_F1 defaults) is kept with the walls FIXED at the spec values
  // (Pal 2018 fired / the 385–400 K estimate), and the default is required to read like the motored walls.
  it('H. ASTM compression pressure within ±3 % with MOTORED wall temperatures (385–400 K, UNVERIFIED) and with the default (lumped) walls; fixed fired walls read +4.5…+5.4 %', () => {
    const cold: EngineSpec = { ...CFR_F1, walls: { ...CFR_F1.walls, headTemperature: 400, pistonTemperature: 400, linerTemperature: 385, intakeValveTemperature: 400, exhaustValveTemperature: 400 } };
    const psig = (spec: EngineSpec, meth: 'RON' | 'MON', counter: number, walls: 'fixed' | 'lumped' = 'fixed'): number => {
      const base = meth === 'RON' ? CFR_RON_CONDITIONS : CFR_MON_CONDITIONS;
      const patch: Partial<OperatingPoint> = meth === 'RON'
        ? { equivalenceRatio: 0, ambientPressure: 101325, ambientTemperature: 325.15, intakeMixtureTemperature: 325.15, compressionRatio: cfrCompressionRatioAtCounter(counter) }
        : { equivalenceRatio: 0, ambientPressure: 101325, ambientTemperature: 311.15, intakeMixtureTemperature: 422.15, compressionRatio: cfrCompressionRatioAtCounter(counter) };
      const m = new CycleModel(spec, { ...base, fuel: { ...base.fuel }, ...patch } as OperatingPoint, { warmupCycles: 8, combustionModel: 'none', wallTemperatureModel: walls });
      m.runCycles(2);
      return (m.runCycles(1)[0].peakPressure - 101325) / 6894.757;
    };
    for (const [meth, c, v] of [['RON', 930, 202.2], ['RON', 778, 169], ['RON', 1061, 241], ['MON', 930, 176.0]] as const) {
      const hot = psig(CFR_F1, meth, c);
      const mot = psig(cold, meth, c);
      const def = psig(CFR_F1, meth, c, 'lumped');
      console.log(`[r2] ASTM ${meth} ${c}: ${v} psig; model fixed fired walls ${hot.toFixed(1)} (${(100 * (hot / v - 1)).toFixed(1)} %), motored walls ${mot.toFixed(1)} (${(100 * (mot / v - 1)).toFixed(1)} %), default lumped walls ${def.toFixed(1)} (${(100 * (def / v - 1)).toFixed(1)} %)`);
      expect(Math.abs(mot / v - 1)).toBeLessThan(0.03);
      expect(hot / v - 1).toBeGreaterThan(mot / v - 1 + 0.015);
      expect(Math.abs(def / mot - 1)).toBeLessThan(0.01);
    }
  });
});
