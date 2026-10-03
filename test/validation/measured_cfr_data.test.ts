/**
 * Validation — agreement with MEASURED CFR data (test/fixtures/cfr_*.json, assembled from
 * Choi et al. 2018, Kolodziej & Wallner 2017, Pal et al. 2018, ASTM D2699/D2700; see
 * test/fixtures/cfr_validation_README.md).
 *
 * [SLOW: ≈ 60–90 s. Runs only with OCTANE_VALIDATION=1, e.g.
 *    OCTANE_VALIDATION=1 npx vitest run test/validation/measured_cfr_data.test.ts ]
 *
 * Public API only (CycleModel through src/physics/cycle/index). Two configurations:
 *  - DEFAULT: CFR_F1 with the shipped defaults, i.e. the ONE global calibration of
 *    src/physics/cycle/calibration.ts (round-1 fixer pass: burnRateMultiplier 5.25,
 *    taylorScaleMultiplier 2.0, Woschni ×1.3, intake-port heat transfer ×6 over the heated port
 *    (Dittus–Boelter, wall 387.8 K), venturi C_D 0.65, Douaud–Eyzat end-gas delay);
 *  - ROUND1: the uncalibrated round-1 state (all multipliers 1, adiabatic intake, ISO venturi C_D
 *    0.984, LLNL-2011 delay tables) — logged for before/after comparison only.
 * History: round 1 compared a validator-PROPOSED set that emulated port heat transfer and the
 * venturi restriction through spec/boundary overrides. Those sub-models now exist as options and
 * the calibration is the default, so the round-1 PROPOSED assertions (tolerances unchanged) now
 * apply to DEFAULT, and the round-1 `it.fails` that document DEFAULT disagreements that are fixed
 * were flipped to `it` (as this file instructed).
 * Tests that document a CURRENT DISAGREEMENT are `it.fails`: when a fix makes them pass, vitest
 * reports them as failing — flip them to `it` then.
 *
 * Measured burn metrics are "apparent" net heat release (single zone, γ = 1.30) evaluated with the
 * SAME routine on the measured and the simulated pressure, so definitional offsets cancel.
 */
import { beforeAll, describe, expect, it } from 'vitest';
import choi from '../fixtures/cfr_choi2018_traces.json';
import astmCp from '../fixtures/cfr_astm_compression_pressure.json';
import guide from '../fixtures/cfr_astm_guide_tables.json';
import pal from '../fixtures/cfr_pal2018_knock.json';
import { CycleModel } from '../../src/physics/cycle/index';
import type { CycleTrace } from '../../src/physics/cycle/index';
import {
  CFR_F1,
  CFR_MON_CONDITIONS,
  CFR_RON_CONDITIONS,
  cfrCompressionRatioAtCounter,
} from '../../src/physics/engines/cfr';
import type { EngineSpec } from '../../src/physics/core/engine-spec';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary } from '../../src/physics/core/snapshot';

const RUN = !!process.env.OCTANE_VALIDATION;
const PSI = 6894.757;
const B = 0.08255;
const S = 0.1143;
const L_ROD = 0.254;
const VD = (Math.PI / 4) * B * B * S;

// ── configurations ──────────────────────────────────────────────────────────────────────────
const ROUND1_OPTIONS = {
  burnRateMultiplier: 1,
  taylorScaleMultiplier: 1,
  woschniMultiplier: 1,
  intakePortHeatTransferMultiplier: 0,
  venturiDischargeCoefficient: 0.984,
  ignitionDelayModel: 'llnl-gasoline-2011',
};

interface Cfg {
  spec: EngineSpec;
  options: Record<string, unknown>;
}
const DEFAULT: Cfg = { spec: CFR_F1, options: {} };
const ROUND1: Cfg = { spec: CFR_F1, options: ROUND1_OPTIONS };

interface Run {
  s: CycleSummary;
  tr: CycleTrace;
  m: CycleModel;
}
function run(cfg: Cfg, base: OperatingPoint, patch: Partial<OperatingPoint>, extra: Record<string, unknown> = {}, cycles = 6, warmup = 3): Run {
  const op: OperatingPoint = { ...base, fuel: { ...base.fuel }, ...patch } as OperatingPoint;
  const m = new CycleModel(cfg.spec, op, { warmupCycles: warmup, ...cfg.options, ...extra } as never);
  if (cycles > 1) m.runCycles(cycles - 1);
  m.recordTrace();
  const s = m.runCycles(1)[0];
  return { s, tr: m.trace!, m };
}

// ── numerics shared by measured and simulated traces ───────────────────────────────────────
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
  const f = (x - xs[lo]) / (xs[hi] - xs[lo]);
  return ys[lo] + f * (ys[hi] - ys[lo]);
}
/** Monotonic (θ, p) of a recorded trace: post-event duplicates removed (the later sample kept). */
function tracePressure(tr: CycleTrace): { th: number[]; p: number[] } {
  const th: number[] = [];
  const p: number[] = [];
  for (let i = 0; i < tr.theta.length; i++) {
    if (th.length && tr.theta[i] <= th[th.length - 1]) {
      if (tr.theta[i] === th[th.length - 1]) p[p.length - 1] = tr.pressure[i];
      continue;
    }
    th.push(tr.theta[i]);
    p.push(tr.pressure[i]);
  }
  return { th, p };
}
/** Apparent net heat release (single zone, γ = 1.30, 1° box smoothing) → CA10/50/90 from the spark. */
function apparentBurn(th: number[], p: number[], cr: number, sparkDeg: number): { ca10: number; ca50: number; ca90: number } {
  const d = 0.2;
  const g: number[] = [];
  for (let x = sparkDeg - 10; x < 120; x += d) g.push(x);
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
    const dP = (P[i1] - P[i0]) / (g[i1] - g[i0]);
    const dV = (V[i1] - V[i0]) / (g[i1] - g[i0]);
    return (gam / (gam - 1)) * P[i] * dV + (1 / (gam - 1)) * V[i] * dP;
  });
  const i0 = g.findIndex((x) => x >= sparkDeg);
  const Q: number[] = [0];
  for (let i = i0 + 1; i < g.length; i++) Q.push(Q[Q.length - 1] + 0.5 * (dQ[i] + dQ[i - 1]) * d);
  let iMax = 0;
  for (let i = 1; i < Q.length; i++) if (Q[i] > Q[iMax]) iMax = i;
  const ca = (f: number): number => {
    const t = f * Q[iMax];
    for (let i = 1; i <= iMax; i++) if (Q[i] >= t) return g[i0 + i - 1] + ((t - Q[i - 1]) / (Q[i] - Q[i - 1])) * d;
    return NaN;
  };
  return { ca10: ca(0.1), ca50: ca(0.5), ca90: ca(0.9) };
}
function gimep(th: number[], p: number[], cr: number): number {
  let w = 0;
  for (let x = -180; x < 180 - 1e-9; x += 0.1) w += 0.5 * (interp(th, p, x) + interp(th, p, x + 0.1)) * (volume(x + 0.1, cr) - volume(x, cr));
  return w / VD;
}
function peak(th: number[], p: number[], lo = -30, hi = 90): { p: number; th: number } {
  let best = { p: -Infinity, th: NaN };
  for (let i = 0; i < th.length; i++) if (th[i] >= lo && th[i] <= hi && p[i] > best.p) best = { p: p[i], th: th[i] };
  return best;
}

// ── measured data ───────────────────────────────────────────────────────────────────────────
type ChoiDs = { key: string; data: { thetaDeg?: number[]; pPa?: number[] } };
const choiTrace = (key: string): { th: number[]; p: number[] } => {
  const ds = (choi.datasets as unknown as ChoiDs[]).find((d) => d.key === key)!;
  return { th: ds.data.thetaDeg!, p: ds.data.pPa! };
};
const CR_CHOI = 7.55; // stated (Choi 2018 p. 2, Pal 2018 p. 2)
/** Choi/Pal PRF100 campaign: CR 7.55 (≈ 29.6 inHg barometer by the D2699 compensation, UNVERIFIED inference), λ 0.89. */
const CHOI_OP: Partial<OperatingPoint> = {
  compressionRatio: CR_CHOI,
  equivalenceRatio: 1 / 0.89,
  fuel: { kind: 'PRF', octaneNumber: 100 },
  ambientPressure: 100240,
  ambientTemperature: 322,
  intakeMixtureTemperature: 302, // KW17 Fig. 3: 28.3 °C after the carburettor (PRF98); TPA port 29 °C (Pal Table 3 via Choi Fig. 11)
};

describe.skipIf(!RUN)('measured CFR data [SLOW]', () => {
  const r: Record<string, Run> = {};
  beforeAll(() => {
    for (const [name, cfg] of [['def', DEFAULT], ['r1', ROUND1]] as const) {
      r[`${name}:p5`] = run(cfg, CFR_RON_CONDITIONS, { ...CHOI_OP, sparkAdvanceDeg: -5 });
      r[`${name}:m13`] = run(cfg, CFR_RON_CONDITIONS, { ...CHOI_OP, sparkAdvanceDeg: 12.72 });
      // Choi 2018 Fig. 2 / Table 6 PRF98 standard RON (barometer ≈ 0.998 bar from the closed-valve intake-port
      // pressure → counter 867 + 6 = 873; UNVERIFIED inference) and KW17 PRF98 (IAT 41.3 °C → 28.6 inHg → 904)
      r[`${name}:choi98`] = run(cfg, CFR_RON_CONDITIONS, {
        compressionRatio: cfrCompressionRatioAtCounter(873), equivalenceRatio: 1 / 0.9047, fuel: { kind: 'PRF', octaneNumber: 98 },
        ambientPressure: 99800, ambientTemperature: 322.9, intakeMixtureTemperature: 301.4, sparkAdvanceDeg: 13,
      });
      r[`${name}:kw98`] = run(cfg, CFR_RON_CONDITIONS, {
        compressionRatio: cfrCompressionRatioAtCounter(904), equivalenceRatio: 1 / 0.8855, fuel: { kind: 'PRF', octaneNumber: 98 },
        ambientPressure: 96850, ambientTemperature: 314.5, intakeMixtureTemperature: 301.4, sparkAdvanceDeg: 13,
      });
    }
    r['def:m13:noknock'] = run(DEFAULT, CFR_RON_CONDITIONS, { ...CHOI_OP, sparkAdvanceDeg: 12.72 }, { knock: false });
  }, 900_000);

  describe('Choi 2018 Fig. 9, PRF100, CR 7.55 (flush-mount transducer, grade A)', () => {
    const m5 = choiTrace('choi2018_fig9_prf100_st_p5_flushMount');
    const m13 = choiTrace('choi2018_fig9_prf100_st_m13_flushMount');

    it('DEFAULT: ST +5 (non-knocking) burn, peak and gIMEP', () => {
      const a = apparentBurn(m5.th, m5.p, CR_CHOI, 5);
      const pm = peak(m5.th, m5.p);
      for (const name of ['r1', 'def']) {
        const sim = tracePressure(r[`${name}:p5`].tr);
        const b = apparentBurn(sim.th, sim.p, CR_CHOI, 5);
        const ps = peak(sim.th, sim.p);
        console.log(`[data] ${name} ST+5 apparent CA10/50/90 meas ${a.ca10.toFixed(1)}/${a.ca50.toFixed(1)}/${a.ca90.toFixed(1)} sim ${b.ca10.toFixed(1)}/${b.ca50.toFixed(1)}/${b.ca90.toFixed(1)}; ` +
          `p_max meas ${(pm.p / 1e5).toFixed(2)}@${pm.th.toFixed(1)} sim ${(ps.p / 1e5).toFixed(2)}@${ps.th.toFixed(1)}; gIMEP meas ${(gimep(m5.th, m5.p, CR_CHOI) / 1e5).toFixed(3)} sim ${(r[`${name}:p5`].s.imepGross / 1e5).toFixed(3)}`);
      }
      const sim = tracePressure(r['def:p5'].tr);
      const b = apparentBurn(sim.th, sim.p, CR_CHOI, 5);
      const ps = peak(sim.th, sim.p);
      expect(Math.abs(b.ca10 - a.ca10)).toBeLessThan(2.5);
      expect(Math.abs(b.ca50 - a.ca50)).toBeLessThan(2.5);
      expect(Math.abs(b.ca90 - a.ca90)).toBeLessThan(4.5);
      expect(Math.abs(ps.p / pm.p - 1)).toBeLessThan(0.1);
      expect(Math.abs(ps.th - pm.th)).toBeLessThan(3);
      expect(Math.abs(r['def:p5'].s.imepGross / gimep(m5.th, m5.p, CR_CHOI) - 1)).toBeLessThan(0.03);
    });

    it('DEFAULT: ST −13 (standard knock) early burn CA10/CA50 within 1.5°, trace RMS (−30…9°) < 3 % of peak', () => {
      const sim = tracePressure(r['def:m13:noknock'].tr);
      const a = apparentBurn(m13.th, m13.p, CR_CHOI, -12.72);
      const b = apparentBurn(sim.th, sim.p, CR_CHOI, -12.72);
      let se = 0;
      let n = 0;
      for (let x = -30; x <= 9; x += 0.25) (se += (interp(sim.th, sim.p, x) - interp(m13.th, m13.p, x)) ** 2), n++;
      const rms = Math.sqrt(se / n) / peak(m13.th, m13.p).p;
      console.log(`[data] ST−13 apparent CA10/50 meas ${a.ca10.toFixed(1)}/${a.ca50.toFixed(1)} sim ${b.ca10.toFixed(1)}/${b.ca50.toFixed(1)}; RMS(−30…9°) ${(100 * rms).toFixed(1)} % of peak`);
      expect(Math.abs(b.ca10 - a.ca10)).toBeLessThan(1.5);
      expect(Math.abs(b.ca50 - a.ca50)).toBeLessThan(1.5);
      expect(rms).toBeLessThan(0.03);
    });

    it('DEFAULT: compression pressure (motored part of ST +5) p(−60°), p(−30°), p(0°) within 4 %', () => {
      const sim = tracePressure(r['def:p5'].tr);
      for (const x of [-60, -30, 0]) expect(Math.abs(interp(sim.th, sim.p, x) / interp(m5.th, m5.p, x) - 1)).toBeLessThan(0.04);
    });

    it('DEFAULT: ST +5 apparent CA50 within 5° of the measured 29.1° (round 1: ≈ 32° late)', () => {
      const sim = tracePressure(r['def:p5'].tr);
      const a = apparentBurn(m5.th, m5.p, CR_CHOI, 5);
      const b = apparentBurn(sim.th, sim.p, CR_CHOI, 5);
      expect(Math.abs(b.ca50 - a.ca50)).toBeLessThan(5);
    });

    it('DEFAULT: knock onset at standard knock within 3° of the measured knock point 11.0 ± 1.0° (Pal 2018 Fig. 9)', () => {
      const meanKp = (pal.datasets[3] as unknown as { data: { statsOfSingleDots: { mean: number } } }).data.statsOfSingleDots.mean;
      const s = r['def:m13'].s;
      console.log(`[data] ST−13 knock onset sim ${s.knockOnsetDeg.toFixed(1)}° (end gas ${(100 * s.knockEndGasFraction).toFixed(1)} %; round 1 ${r['r1:m13'].s.knockOnsetDeg.toFixed(1)}°) vs measured knock point ${meanKp.toFixed(2)}°`);
      expect(Math.abs(s.knockOnsetDeg - meanKp)).toBeLessThan(3);
    });

    it('DEFAULT: MAPO at the standard-knock state within the measured 0.7–2.3 bar (Hoth 2021b 1.43 bar 4–18 kHz; Pal 2018 2.27 bar)', () => {
      // round 1: MAPO 0 (acoustic source vanished once the front radius covered the head plane)
      const s = r['def:m13'].s;
      console.log(`[data] ST−13 MAPO sim ${(s.mapo / 1e5).toFixed(2)} bar, max dp/dθ ${(s.maxPressureRiseRate / 1e5).toFixed(2)} bar/°`);
      expect(s.knockEndGasFraction).toBeGreaterThan(0.01);
      expect(s.mapo).toBeGreaterThan(0.7e5);
      expect(s.mapo).toBeLessThan(2.3e5);
    });
  });

  describe('breathing and efficiency: trapped fuel vs measured fuel flow (Choi 2018 Table 6, KW17 Fig. 2)', () => {
    // 600 rpm four-stroke: 5 cycles/s → kg/h = m_f·5·3600
    const kgh = (x: Run): number => x.m.fuelMassIvc * 5 * 3600;
    it('DEFAULT: fuel flow within 5 % (Choi 0.737 kg/h, KW17 0.747 kg/h; round 1: +11 % / +9 %)', () => {
      console.log(`[data] fuel kg/h — Choi98 r1 ${kgh(r['r1:choi98']).toFixed(3)} def ${kgh(r['def:choi98']).toFixed(3)} (0.737); KW98 r1 ${kgh(r['r1:kw98']).toFixed(3)} def ${kgh(r['def:kw98']).toFixed(3)} (0.747); ` +
        `T_IVC ST−13 r1 ${r['r1:m13'].m.TIvc.toFixed(0)} def ${r['def:m13'].m.TIvc.toFixed(0)} K (TPA 410); m_trap r1 ${(r['r1:m13'].s.trappedMass * 1e6).toFixed(0)} def ${(r['def:m13'].s.trappedMass * 1e6).toFixed(0)} mg (TPA 628)`);
      expect(Math.abs(kgh(r['def:choi98']) / 0.737 - 1)).toBeLessThan(0.05);
      expect(Math.abs(kgh(r['def:kw98']) / 0.747 - 1)).toBeLessThan(0.05);
    });
    // Pal 2018 Table 3 / Choi TPA trapped mass is a GT-Power MODEL value; with the measured fuel flow matched
    // the model traps ≈ +5 % more (T_IVC 383 K vs the TPA's 406–410 K; residual 6.5 % vs 6.0 %).
    it.fails('DEFAULT: trapped mass within 5 % of the TPA model value 628 mg [TPA is a model output]', () => {
      expect(Math.abs(r['def:m13'].s.trappedMass / 0.628e-3 - 1)).toBeLessThan(0.05);
    });
    it('DEFAULT: work per unit fuel within 4 % of Choi Table 6 (gIMEP 7.996 bar / 0.737 kg/h) and KW17 (round 1 proposed: +7 %)', () => {
      const x = r['def:choi98'];
      const ratio = (x.s.imepGross / kgh(x)) / (7.9957e5 / 0.737);
      const ratioKw = (r['def:kw98'].s.imepGross / kgh(r['def:kw98'])) / (7.892e5 / 0.747);
      console.log(`[data] gIMEP/fuel model/measured = ${ratio.toFixed(3)} (Choi98), KW98 ${ratioKw.toFixed(3)}`);
      expect(Math.abs(ratio - 1)).toBeLessThan(0.04);
      expect(Math.abs(ratioKw - 1)).toBeLessThan(0.04);
    });
  });

  describe('ASTM D2699/D2700 motored compression pressure (Table A2.2, Fig. 2; grade A)', () => {
    const pts: [('RON' | 'MON'), number, number][] = [
      ['RON', 930, 202.2], ['RON', 778, 169], ['RON', 1061, 241], ['MON', 930, 176.0], ['MON', 578, 120], ['MON', 1008, 194],
    ];
    void astmCp;
    const gauge = (cfg: Cfg, meth: 'RON' | 'MON', counter: number, venturiScale = 1): number => {
      const spec = venturiScale === 1 ? cfg.spec : { ...cfg.spec, manifolds: { ...cfg.spec.manifolds, throttleDiameter: cfg.spec.manifolds.throttleDiameter * venturiScale } };
      const base = meth === 'RON' ? CFR_RON_CONDITIONS : CFR_MON_CONDITIONS;
      // carburettor drained: air only (φ = 0 is accepted since the fixer pass; round 1 had to use φ 0.2)
      const patch: Partial<OperatingPoint> = meth === 'RON'
        ? { equivalenceRatio: 0, ambientPressure: 101325, ambientTemperature: 325.15, intakeMixtureTemperature: 325.15, compressionRatio: cfrCompressionRatioAtCounter(counter) }
        : { equivalenceRatio: 0, ambientPressure: 101325, ambientTemperature: 311.15, intakeMixtureTemperature: 422.15, compressionRatio: cfrCompressionRatioAtCounter(counter) };
      const x = run({ ...cfg, spec }, base, patch, { combustionModel: 'none' }, 3, 8);
      return (x.s.peakPressure - 101325) / PSI;
    };
    const res: Record<string, number> = {};
    beforeAll(() => {
      for (const [meth, c] of pts) {
        res[`def:${meth}:${c}`] = gauge(DEFAULT, meth, c);
        res[`r1:${meth}:${c}`] = gauge(ROUND1, meth, c);
      }
      // D2700 Fig. 2: 3/4 in venturi reads 190.1 psig vs 176.0 (9/16) at counter 930, 29.92 inHg
      res['def:MON:930:3/4'] = gauge(DEFAULT, 'MON', 930, (3 / 4) / (9 / 16));
      res['r1:MON:930:3/4'] = gauge(ROUND1, 'MON', 930, (3 / 4) / (9 / 16));
      console.log('[data] ASTM compression psig (ASTM / round 1 / default): ' +
        pts.map(([m, c, v]) => `${m} ${c}: ${v} / ${res[`r1:${m}:${c}`].toFixed(1)} / ${res[`def:${m}:${c}`].toFixed(1)}`).join('; ') +
        `; MON 930 3/4 in: 190.1 / ${res['r1:MON:930:3/4'].toFixed(1)} / ${res['def:MON:930:3/4'].toFixed(1)}`);
    }, 900_000);
    // The absolute level inherits the counter→CR relation (±0.4 CR between ANL campaigns ≈ ±5 % in
    // pressure, cfr_validation_README.md): the model reads +4…+6 % at every check point while it is
    // −1.7 % below the Choi 2018 fired-engine compression at the stated CR 7.55 (test above).
    it.fails('DEFAULT: all six check points within ±3 % [absolute level: counter→CR relation]', () => {
      for (const [meth, c, v] of pts) expect(Math.abs(res[`def:${meth}:${c}`] / v - 1)).toBeLessThan(0.03);
    });
    it('DEFAULT: all six check points within ±7 % (±0.4 CR calibration band), uniform to ±2 %', () => {
      const e = pts.map(([meth, c, v]) => res[`def:${meth}:${c}`] / v - 1);
      for (const x of e) expect(Math.abs(x)).toBeLessThan(0.07);
      expect(Math.max(...e) - Math.min(...e)).toBeLessThan(0.04);
    });
    it('DEFAULT: MON/RON compression-pressure ratio at counter 930 within 2 % (round 1: +5.8 %)', () => {
      const ratio = (res['def:MON:930'] / res['def:RON:930']) / (176.0 / 202.2);
      expect(Math.abs(ratio - 1)).toBeLessThan(0.02);
    });
    it('DEFAULT: MON venturi-size effect within 3.5 points (measured +8.0 %; round 1 model +2.7 %)', () => {
      const d = res['def:MON:930:3/4'] / res['def:MON:930'] - 1;
      expect(Math.abs(d - (190.1 / 176.0 - 1))).toBeLessThan(0.035);
    });
  });

  describe('octane scale with the LLNL-2011 detailed-chemistry delay (option; not the default)', () => {
    // Criterion: MAPO 0.67 bar (ASTM standard knock, Rockstroh 2018). The default (Douaud–Eyzat) is
    // validated in knock_octane_rating.test.ts.
    const guideCR = (on: number): number => {
      const t = (guide as unknown as { RON_9_16: { onToCounter: { octaneNumber: number[]; counter: number[] } } }).RON_9_16.onToCounter;
      return cfrCompressionRatioAtCounter(interp(t.octaneNumber, t.counter, on));
    };
    const klcr = (on: number): number => {
      let lo = 5;
      let hi = 11;
      const base: OperatingPoint = { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: on } };
      const m = new CycleModel(CFR_F1, { ...base, compressionRatio: 8 }, { warmupCycles: 3, ignitionDelayModel: 'llnl-gasoline-2011' });
      for (let it = 0; it < 8; it++) {
        const cr = 0.5 * (lo + hi);
        m.setOperatingPoint({ compressionRatio: cr });
        m.runCycles(2);
        if (m.runCycles(1)[0].mapo < 0.67e5) lo = cr;
        else hi = cr;
      }
      return 0.5 * (lo + hi);
    };
    it.fails('LLNL-2011 τ (single-stage Livengood–Wu): knock-limited CR of PRF 90 within ±0.35 and dCR/dON (80 → 100) within ±30 % of the guide table', () => {
      const c80 = klcr(80);
      const c90 = klcr(90);
      const c100 = klcr(100);
      const sm = (c100 - c80) / 20;
      const sg = (guideCR(100) - guideCR(80)) / 20;
      console.log(`[data] LLNL RON knock-limited CR (MAPO 0.67 bar): PRF80 ${c80.toFixed(2)} PRF90 ${c90.toFixed(2)} PRF100 ${c100.toFixed(2)}; slope ${sm.toFixed(4)} vs guide ${sg.toFixed(4)} /ON`);
      expect(Math.abs(c90 - guideCR(90))).toBeLessThan(0.35);
      expect(Math.abs(sm / sg - 1)).toBeLessThan(0.3);
    });
  });
});
