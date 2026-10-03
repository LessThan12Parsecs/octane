import { describe, expect, it } from 'vitest';
import { NS, SP } from '../core/species';
import type { FuelBlend } from '../thermo/fuels';
import { fuelFromSelection } from '../thermo/fuels';
import griFx from '../../../test/fixtures/flamespeed_cantera.json';
import prfFx from '../../../test/fixtures/flamespeed_prf_cantera.json';
import expFx from '../../../test/fixtures/flamespeed_experiments.json';
import fitFx from '../../../test/fixtures/flamespeed_fit.json';
import dilFx from '../../../test/fixtures/flamespeed_dilution.json';
import { EquilibriumSolver } from '../equilibrium/solver';
import {
  ACTIVATION_ENERGY_TABLE,
  activationEnergy,
  amirante2017,
  AMIRANTE_2017,
  binaryDiffusionCoefficient,
  burnedGasState,
  effectiveLewisNumber,
  flameThickness,
  flammabilityLimits,
  FUEL_LHV_MOLAR,
  gulder1984,
  interpolateTable,
  GULDER_1984,
  laminarFlameSpeed,
  laminarFlameSpeedModel,
  clearFlameSpeedCache,
  lewisNumbers,
  marksteinIntegrals,
  marksteinLength,
  marksteinLengths,
  metghalchiKeck1982,
  METGHALCHI_KECK_1982,
  mixtureDiffusionCoefficient,
  rhodesKeck1985Indolene,
  rhodesKeckDilutionFactor,
  TABLE_P,
  TABLE_PHI,
  TABLE_TU,
  unburnedComposition,
} from './laminar-flame-speed';

type Species = 'IC8H18' | 'NC7H16' | 'CH4' | 'C3H8' | 'C2H5OH';
const pure = (s: Species): FuelBlend => fuelFromSelection({ kind: 'pure', species: s });
const blend = (xIso: number): FuelBlend => {
  const X = new Float64Array(NS);
  X[SP.IC8H18] = xIso;
  X[SP.NC7H16] = 1 - xIso;
  return { label: `x_iso=${xIso}`, X };
};
const rel = (a: number, b: number): number => a / b - 1;
const stats = (e: number[]): { mean: number; max: number; bias: number } => ({
  mean: e.reduce((s, v) => s + Math.abs(v), 0) / e.length,
  max: e.reduce((s, v) => Math.max(s, Math.abs(v)), 0),
  bias: e.reduce((s, v) => s + v, 0) / e.length,
});
const pct = (v: number): string => `${(100 * v).toFixed(1)} %`;

interface Flame {
  fuel: string | Record<string, number>;
  set?: string;
  phi: number;
  Tu: number;
  p: number;
  ydil: number;
  SL: number;
  // absent for flames recovered from a run log (S_L only)
  Tb?: number;
  sigma?: number;
  Tad?: number;
  sigmaEq?: number;
  deltaT?: number;
  alphaU?: number;
  Ddef?: number;
}
const gri = griFx.flames as unknown as Flame[];
const prf = prfFx.flames as unknown as Flame[];
const single = (f: Flame): Species | null => {
  if (typeof f.fuel === 'string') return f.fuel as Species;
  const k = Object.keys(f.fuel);
  return k.length === 1 ? (k[0] as Species) : null;
};
const onGrid = (f: Flame): boolean =>
  TABLE_PHI.some((v) => Math.abs(v - f.phi) < 1e-9) &&
  TABLE_TU.includes(f.Tu) &&
  TABLE_P.some((v) => Math.abs(v - f.p) < 1);

describe('table interpolation (tensor-product cubic Hermite, C¹ extrapolation)', () => {
  // A function quadratic in each of φ, ln Tu, ln p (with cross terms) is reproduced exactly:
  // the 3-point node slopes are exact for quadratics and cubic Hermite then interpolates exactly.
  const q = (phi: number, T: number, p: number): number => {
    const a = phi - 1;
    const b = Math.log(T / 300);
    const c = Math.log(p / 1e5);
    return -1 + 0.3 * a - 2 * a * a + 1.7 * b + 0.4 * b * b - 0.3 * c - 0.01 * c * c + 0.08 * b * c + 0.2 * a * b * c;
  };
  const vals: number[] = [];
  for (const phi of TABLE_PHI) for (const T of TABLE_TU) for (const p of TABLE_P) vals.push(q(phi, T, p));
  const tab = Float64Array.from(vals);

  it('reproduces quadratics exactly inside the grid', () => {
    let worst = 0;
    for (let phi = 0.6; phi <= 1.6; phi += 0.037) {
      for (let T = 300; T <= 900; T += 47) {
        for (let lp = Math.log(1e5); lp <= Math.log(60e5); lp += 0.31) {
          worst = Math.max(worst, Math.abs(interpolateTable(tab, phi, T, Math.exp(lp)) - q(phi, T, Math.exp(lp))));
        }
      }
    }
    expect(worst).toBeLessThan(1e-11);
  });

  it('is C¹ across nodes and at the grid boundary (linear continuation outside)', () => {
    const h = 1e-6;
    for (const phi of [0.6, 0.7, 1.2, 1.4, 1.6]) {
      const l = (interpolateTable(tab, phi, 555, 7e5) - interpolateTable(tab, phi - h, 555, 7e5)) / h;
      const r = (interpolateTable(tab, phi + h, 555, 7e5) - interpolateTable(tab, phi, 555, 7e5)) / h;
      expect(Math.abs(l - r)).toBeLessThan(1e-4);
    }
    // beyond the grid: exactly linear in the extrapolated coordinate
    const f1 = interpolateTable(tab, 1.8, 600, 1e6);
    const f2 = interpolateTable(tab, 2.0, 600, 1e6);
    const f3 = interpolateTable(tab, 2.2, 600, 1e6);
    expect(Math.abs(f3 - 2 * f2 + f1)).toBeLessThan(1e-12);
  });
});

describe('tabulated detailed-chemistry burning velocities', () => {
  it('reproduce every undiluted grid flame of the oracles (interpolant passes through the nodes)', () => {
    const cases: [Species, Flame[]][] = [
      ['CH4', gri.filter((f) => f.fuel === 'CH4' && f.ydil === 0 && onGrid(f))],
      ['IC8H18', prf.filter((f) => f.set === 'grid' && single(f) === 'IC8H18' && f.ydil === 0)],
      ['NC7H16', prf.filter((f) => f.set === 'grid' && single(f) === 'NC7H16' && f.ydil === 0)],
    ];
    for (const [sp, rows] of cases) {
      expect(rows.length).toBeGreaterThan(100);
      let worst = 0;
      for (const f of rows) {
        const s = laminarFlameSpeedModel('default', sp, f.phi, f.Tu, f.p);
        worst = Math.max(worst, Math.abs(rel(s, f.SL)));
      }
      console.info(`[flamespeed] ${sp}: ${rows.length} grid flames reproduced to max ${(100 * worst).toFixed(3)} %`);
      expect(worst).toBeLessThan(5e-4); // tables are stored with 5 decimals of ln S_L
    }
  });

  it('interpolate off-grid flames (CH4 at 800 K / 40 bar; PRF blends and propane checks) to ≲ 3 %', () => {
    const off = gri.filter((f) => f.fuel === 'CH4' && f.ydil === 0 && !onGrid(f));
    expect(off.length).toBeGreaterThanOrEqual(3);
    const e = off.map((f) => rel(laminarFlameSpeedModel('default', 'CH4', f.phi, f.Tu, f.p), f.SL));
    const st = stats(e);
    console.info(`[flamespeed] CH4 off-grid (${off.length} flames): mean ${pct(st.mean)}, max ${pct(st.max)}`);
    expect(st.max).toBeLessThan(0.03);
    const jz = prf.filter((f) => f.set === 'jerzembeck');
    for (const sp of ['IC8H18', 'NC7H16'] as const) {
      const rows = jz.filter((f) => single(f) === sp);
      if (rows.length === 0) continue;
      const ej = rows.map((f) => rel(laminarFlameSpeedModel('default', sp, f.phi, f.Tu, f.p), f.SL));
      const sj = stats(ej);
      console.info(`[flamespeed] ${sp} off-grid at 373 K, 10–25 bar (${rows.length}): mean ${pct(sj.mean)}, max ${pct(sj.max)}`);
      expect(sj.max).toBeLessThan(0.03);
    }
  });
});

describe('comparison with published experiments (test/fixtures/flamespeed_experiments.json)', () => {
  const pts = expFx.points as { fuel: Species; phi: number; Tu: number; p: number; SL: number; dataset: string }[];
  const inRange = (r: { phi: number }): boolean => r.phi >= 0.7 && r.phi <= 1.4;
  // Jerzembeck et al. burned "air" with X_O2 = 0.205: convert model values with the ratio of
  // PRF_HT flames computed with that oxidiser and with O2 + 3.76 N2 (interpolated in φ, ln p).
  const o2 = (fuel: Species, phi: number, p: number): number => {
    const pair = (tag: string): Flame[] => prf.filter((f) => f.set === tag && single(f) === fuel);
    const a = pair('jerzembeck-o2');
    const b = pair('jerzembeck');
    if (a.length === 0) return 1;
    const ratioAt = (ph: number, pp: number): number => {
      const fa = a.find((f) => Math.abs(f.phi - ph) < 1e-9 && Math.abs(f.p - pp) < 1)!;
      const fb = b.find((f) => Math.abs(f.phi - ph) < 1e-9 && Math.abs(f.p - pp) < 1)!;
      return fa.SL / fb.SL;
    };
    const phis = [0.7, 0.8, 1.0, 1.2];
    let i = 0;
    while (i < phis.length - 2 && phi > phis[i + 1]) i++;
    const t = Math.min(Math.max((phi - phis[i]) / (phis[i + 1] - phis[i]), 0), 1);
    const s = Math.min(Math.max(Math.log(p / 10e5) / Math.log(2.5), 0), 1);
    const r = (pp: number): number => (1 - t) * ratioAt(phis[i], pp) + t * ratioAt(phis[i + 1], pp);
    return (1 - s) * r(10e5) + s * r(25e5);
  };

  it('default model vs measurements, per fuel and data set', () => {
    const limits: Record<string, number> = { CH4: 0.1, C3H8: 0.08, IC8H18: 0.1, NC7H16: 0.15 };
    for (const fuel of ['CH4', 'C3H8', 'IC8H18', 'NC7H16'] as const) {
      const rows = pts.filter((r) => r.fuel === fuel && inRange(r));
      const oxid = (r: (typeof rows)[number]): number =>
        r.dataset === 'Jerzembeck 09' && (fuel === 'IC8H18' || fuel === 'NC7H16') && r.phi <= 1.2 ? o2(fuel, r.phi, r.p) : 1;
      const e = rows.map((r) => rel(laminarFlameSpeed(pure(fuel), r.phi, r.Tu, r.p, 0) * oxid(r), r.SL));
      const st = stats(e);
      const byDs = new Map<string, number[]>();
      rows.forEach((r, i) => byDs.set(r.dataset, [...(byDs.get(r.dataset) ?? []), e[i]]));
      const detail = [...byDs].map(([k, v]) => `${k} ${pct(stats(v).bias)}`).join(', ');
      console.info(`[flamespeed] ${fuel} default vs ${rows.length} measurements (φ 0.7–1.4): mean |e| ${pct(st.mean)}, bias ${pct(st.bias)}, max ${pct(st.max)} — ${detail}`);
      expect(st.mean).toBeLessThan(limits[fuel]);
    }
  });

  it('high-pressure iso-octane / n-heptane data (Jerzembeck 2009, Kelley 2011): default vs published correlations', () => {
    for (const fuel of ['IC8H18', 'NC7H16'] as const) {
      const rows = pts.filter((r) => r.fuel === fuel && r.p > 1.5e5 && r.phi >= 0.7 && r.phi <= 1.2);
      expect(rows.length).toBeGreaterThan(10);
      const corr = (r: (typeof rows)[number]): number => (r.dataset === 'Jerzembeck 09' ? o2(fuel, r.phi, r.p) : 1);
      const models: [string, (r: (typeof rows)[number]) => number][] = [
        ['default (PRF_HT table)', (r) => laminarFlameSpeed(pure(fuel), r.phi, r.Tu, r.p, 0)],
        ['Metghalchi–Keck 1982 isooctane', (r) => metghalchiKeck1982(METGHALCHI_KECK_1982.isooctane, r.phi, r.Tu, r.p)],
        ['Gülder 1984 isooctane', (r) => gulder1984(GULDER_1984.IC8H18, r.phi, r.Tu, r.p)],
        ['Rhodes–Keck 1985 indolene', (r) => rhodesKeck1985Indolene(r.phi, r.Tu, r.p)],
      ];
      // all models describe standard air → the experiments' oxidiser correction applies to all
      const res = models.map(([name, fn]) => [name, stats(rows.map((r) => rel(fn(r) * corr(r), r.SL)))] as const);
      for (const [name, st] of res) {
        console.info(`[flamespeed] ${fuel} p ≥ 2 bar (${rows.length} pts) ${name}: mean |e| ${pct(st.mean)}, bias ${pct(st.bias)}, max ${pct(st.max)}`);
      }
      expect(res[0][1].mean).toBeLessThan(0.2);
    }
  });
});

describe('published correlations vs the Cantera oracles (report)', () => {
  it('CH4 and C3H8 against GRI-3.0 at engine-like conditions (Tu ≥ 500 K, p ≥ 5 bar, φ 0.8–1.2)', () => {
    for (const sp of ['CH4', 'C3H8'] as const) {
      const rows = gri.filter((f) => f.fuel === sp && f.ydil === 0 && f.Tu >= 500 && f.p >= 5e5 && f.phi >= 0.8 && f.phi <= 1.2);
      const models: [string, (f: Flame) => number][] = [
        ['Gülder 1984', (f) => gulder1984(GULDER_1984[sp], f.phi, f.Tu, f.p)],
        ['Amirante 2017', (f) => amirante2017(AMIRANTE_2017[sp], f.phi, f.Tu, f.p)],
      ];
      if (sp === 'C3H8') models.push(['Metghalchi–Keck 1982', (f) => metghalchiKeck1982(METGHALCHI_KECK_1982.propane, f.phi, f.Tu, f.p)]);
      models.push(['default', (f) => laminarFlameSpeed(pure(sp), f.phi, f.Tu, f.p, 0)]);
      for (const [name, fn] of models) {
        const st = stats(rows.map((f) => rel(fn(f), f.SL)));
        console.info(`[flamespeed] ${sp} vs GRI-3.0 (${rows.length}) ${name}: mean |e| ${pct(st.mean)}, bias ${pct(st.bias)}, max ${pct(st.max)}`);
      }
    }
    // GRI-3.0 is not a propane oracle: it over-predicts measured C3H8 S_L by ~20–27 %.
    const exp298 = (expFx.points as { fuel: string; phi: number; Tu: number; p: number; SL: number }[]).filter(
      (r) => r.fuel === 'C3H8' && Math.abs(r.phi - 1) < 0.02 && r.Tu === 298 && r.p === 1e5,
    );
    const g = gri.find((f) => f.fuel === 'C3H8' && f.phi === 1 && f.Tu === 300 && f.p === 1e5 && f.ydil === 0)!;
    const meas = exp298.reduce((s, r) => s + r.SL, 0) / exp298.length;
    console.info(`[flamespeed] C3H8 φ=1, 1 bar: GRI-3.0 ${g.SL.toFixed(3)} m/s (300 K) vs measured mean ${meas.toFixed(3)} m/s (298 K)`);
    expect(g.SL / meas).toBeGreaterThan(1.1);
  });

  it('Metghalchi–Keck (1982): Table 3A constants as printed; comparison with the per-φ fits of Table 2', () => {
    // At φ = φ_m, 298 K, 1 atm the correlation returns B_m exactly (Table 3A, verified on the scan).
    expect(metghalchiKeck1982(METGHALCHI_KECK_1982.isooctane, 1.13, 298, 101325)).toBeCloseTo(0.2632, 12);
    expect(metghalchiKeck1982(METGHALCHI_KECK_1982.propane, 1.08, 298, 101325)).toBeCloseTo(0.3422, 12);
    // Table 2 (individual power-law fits, S_u0 in cm/s at φ = 0.8 / 1.0 / 1.2): isooctane
    // 19.25 / 27.00 / 27.63. B_m, B_2 come from the paper's overall fit ("C_m are significantly
    // greater than the corresponding B_m obtained from the overall fit", p. 201), so Eq. (18)
    // sits 6–11 % below Table 2 — a property of the published correlation.
    const t2 = [[0.8, 0.1925], [1.0, 0.27], [1.2, 0.2763]];
    const e = t2.map(([phi, s]) => rel(metghalchiKeck1982(METGHALCHI_KECK_1982.isooctane, phi, 298, 101325), s));
    console.info(`[flamespeed] MK1982 Eq.(18) vs its Table 2 S_u0 (isooctane): ${e.map(pct).join(', ')}`);
    expect(Math.max(...e.map(Math.abs))).toBeLessThan(0.12);
  });
});

describe('dilution, blending and flammability limits', () => {
  it('thermal-theory dilution factor vs Cantera flames diluted with CO2/H2O/N2 products (Rhodes–Keck for reference)', () => {
    const groups = new Map<string, Flame[]>();
    const all = gri.concat(prf);
    for (const f of all) {
      const sp = single(f);
      if (!sp || f.phi !== 1 || (f.set !== undefined && f.set !== 'grid' && f.set !== 'dilution')) continue;
      if (!all.some((g) => single(g) === sp && g.ydil > 0 && g.Tu === f.Tu && g.p === f.p)) continue;
      const k = `${sp}@${f.Tu}K/${f.p / 1e5}bar`;
      groups.set(k, [...(groups.get(k) ?? []), f]);
    }
    const errs: number[] = [];
    const errsC3: number[] = [];
    const errsRK: number[] = [];
    for (const [k, rows] of groups) {
      const base = rows.find((f) => f.ydil === 0);
      if (!base) continue;
      const sp = single(base)!;
      const line: string[] = [];
      for (const f of rows.filter((r) => r.ydil > 0).sort((a, b) => a.ydil - b.ydil)) {
        const ratio = f.SL / base.SL;
        const model = laminarFlameSpeed(pure(sp), 1, f.Tu, f.p, f.ydil) / laminarFlameSpeed(pure(sp), 1, f.Tu, f.p, 0);
        // propane's dilution response is calibrated on PRF_HT (flamespeed_dilution.json); GRI-3.0 is
        // not a propane oracle (+25 % in S_L) and its diluted flames lose ≈ 0.035 less: cross-mechanism
        // check with its own tolerance
        (sp === 'C3H8' ? errsC3 : errs).push(model - ratio);
        errsRK.push(rhodesKeckDilutionFactor(f.ydil) - ratio);
        line.push(`x=${f.ydil}: oracle ${ratio.toFixed(3)} model ${model.toFixed(3)} RK ${rhodesKeckDilutionFactor(f.ydil).toFixed(3)}`);
      }
      console.info(`[flamespeed] dilution ${k}: ${line.join('; ')}`);
    }
    expect(errs.length).toBeGreaterThan(10);
    const st = stats(errs);
    const sr = stats(errsRK);
    console.info(`[flamespeed] dilution factor, absolute error in S_L(x)/S_L(0): thermal model max ${st.max.toFixed(3)} (C3H8 vs GRI-3.0 ${stats(errsC3).max.toFixed(3)}), Rhodes–Keck max ${sr.max.toFixed(3)}`);
    expect(st.max).toBeLessThan(0.04);
    expect(errsC3.length).toBeGreaterThan(3);
    expect(stats(errsC3).max).toBeLessThan(0.05);
  });

  // Paired diluted/undiluted flames with the model's diluent definition (complete-combustion
  // products of the SAME mixture) at φ = 0.6–1.4 (tools/reference/flamespeed_dilution.py).
  type DilFlame = Flame & { set: 'fit' | 'check'; fuel: Species };
  const dil = (dilFx as unknown as { flames: DilFlame[] }).flames;
  const dilRatios = (set: 'fit' | 'check'): { f: DilFlame; oracle: number; model: number }[] => {
    const out: { f: DilFlame; oracle: number; model: number }[] = [];
    for (const f of dil.filter((r) => r.set === set && r.ydil > 0)) {
      const b = dil.find((r) => r.set === set && r.fuel === f.fuel && r.phi === f.phi && r.Tu === f.Tu && r.p === f.p && r.ydil === 0);
      if (!b) continue;
      const fu = pure(f.fuel);
      // the propagation cut-off (flammabilityLimits) zeroes adiabatic flames of a few cm/s near the
      // limits (e.g. CH4 φ = 0.6, 300 K, 1 bar, x = 0.15: 4.6 cm/s): compare only where it is inactive
      const lim = flammabilityLimits(fu, f.Tu, f.ydil);
      if (f.phi < 1.06 * lim.lean || f.phi > 0.94 * lim.rich) continue;
      const model = laminarFlameSpeed(fu, f.phi, f.Tu, f.p, f.ydil) / laminarFlameSpeed(fu, f.phi, f.Tu, f.p, 0);
      out.push({ f, oracle: f.SL / b.SL, model });
    }
    return out;
  };

  it('dilution factor reproduces its fitting flames (φ 0.6–1.4, 1–60 bar) at the table nodes', () => {
    const r = dilRatios('fit');
    expect(r.length).toBeGreaterThan(100);
    const worst = Math.max(...r.map((v) => Math.abs(v.model / v.oracle - 1)));
    console.info(`[flamespeed] dilution fit set (${r.length} diluted flames): max rel. error ${pct(worst)}`);
    expect(worst).toBeLessThan(0.005);
  });

  it('dilution factor vs INDEPENDENT diluted flames off the fitting set (other φ, Tu, p, x)', () => {
    // A stoichiometric-only calibration missed these by up to 0.1 in S_L(x)/S_L(0) (review).
    const r = dilRatios('check');
    expect(r.length).toBeGreaterThanOrEqual(9);
    for (const v of r) {
      console.info(
        `[flamespeed] dilution check ${v.f.fuel} φ=${v.f.phi} ${v.f.Tu} K ${v.f.p / 1e5} bar x=${v.f.ydil}: ` +
          `oracle ${v.oracle.toFixed(3)} model ${v.model.toFixed(3)} Rhodes–Keck ${rhodesKeckDilutionFactor(v.f.ydil).toFixed(3)}`,
      );
    }
    const abs = r.map((v) => Math.abs(v.model - v.oracle));
    const relErr = r.map((v) => Math.abs(v.model / v.oracle - 1));
    const meanRel = relErr.reduce((a, b) => a + b, 0) / relErr.length;
    console.info(`[flamespeed] dilution check (${r.length}): max |ΔF| ${Math.max(...abs).toFixed(3)}, max rel. ${pct(Math.max(...relErr))}, mean rel. ${pct(meanRel)}`);
    // residual: E_a(φ, p) is calibrated along the engine-like (Tu, p) diagonal, so off-diagonal
    // states (e.g. CH4 φ = 1.1 at 600 K / 40 bar) keep up to ≈ 9 % error in S_L(x)/S_L(0)
    expect(Math.max(...abs)).toBeLessThan(0.04);
    expect(Math.max(...relErr)).toBeLessThan(0.1);
    expect(meanRel).toBeLessThan(0.035);
  });

  it('effective activation energy: lower off stoichiometry than at φ = 1 (the reason for E_a(φ))', () => {
    for (const f of [2, 0, 1]) {
      // CH4, IC8H18, NC7H16 at 20 bar, x = 0.15 node
      const e1 = activationEnergy(f, 1.0, 20e5, 0);
      expect(activationEnergy(f, 0.8, 20e5, 0)).toBeLessThan(e1);
    }
  });

  it('energy-fraction blending rule reproduces PRF_HT flames of iso-octane/n-heptane blends', () => {
    const bl = prf.filter((f) => f.set === 'blend');
    expect(bl.length).toBeGreaterThan(5);
    const e: number[] = [];
    for (const f of bl) {
      const x = (f.fuel as Record<string, number>).IC8H18;
      const pureIso = prf.find((g) => g.set === 'grid' && single(g) === 'IC8H18' && g.phi === f.phi && g.Tu === f.Tu && g.p === f.p)!;
      const pureHep = prf.find((g) => g.set === 'grid' && single(g) === 'NC7H16' && g.phi === f.phi && g.Tu === f.Tu && g.p === f.p)!;
      const q1 = x * FUEL_LHV_MOLAR[0];
      const q2 = (1 - x) * FUEL_LHV_MOLAR[1];
      const rule = (q1 * pureIso.SL + q2 * pureHep.SL) / (q1 + q2);
      e.push(rel(rule, f.SL));
      expect(Math.abs(rel(laminarFlameSpeed(blend(x), f.phi, f.Tu, f.p, 0), f.SL))).toBeLessThan(0.03);
    }
    const st = stats(e);
    console.info(`[flamespeed] energy-fraction rule vs PRF_HT blend flames (${bl.length}): mean ${pct(st.mean)}, max ${pct(st.max)}`);
    expect(st.max).toBeLessThan(0.03);
  });

  it('flammability limits: zero outside, positive inside, widen with temperature, narrow with dilution', () => {
    const f = pure('CH4');
    const lim = flammabilityLimits(f, 298.15);
    expect(lim.lean).toBeCloseTo((0.05 / 0.95) * 2 / 0.20939 * 0.99997, 2);
    expect(laminarFlameSpeed(f, lim.lean * 0.999, 298.15, 1e5, 0)).toBe(0);
    expect(laminarFlameSpeed(f, lim.rich * 1.001, 298.15, 1e5, 0)).toBe(0);
    expect(laminarFlameSpeed(f, 1, 300, 1e5, 0)).toBeGreaterThan(0.3);
    const hot = flammabilityLimits(f, 800);
    expect(hot.lean).toBeLessThan(lim.lean);
    expect(hot.rich).toBeGreaterThan(lim.rich);
    const dil = flammabilityLimits(f, 298.15, 0.2);
    expect(dil.lean).toBeGreaterThan(lim.lean);
    // PRF 90 lean limit between those of the pure components
    const p90 = flammabilityLimits(fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }), 298.15);
    const iso = flammabilityLimits(pure('IC8H18'), 298.15);
    const hep = flammabilityLimits(pure('NC7H16'), 298.15);
    expect(p90.lean).toBeGreaterThanOrEqual(Math.min(iso.lean, hep.lean) - 1e-12);
    expect(p90.lean).toBeLessThanOrEqual(Math.max(iso.lean, hep.lean) + 1e-12);
  });
});

describe('ethanol (Vancoillie et al. 2012 correlation, tabulated)', () => {
  const eth = pure('C2H5OH');
  const pts = (fitFx as unknown as { ethanolVancoillie2012: { points: { phi: number; Tu: number; p: number; SL: number }[] } })
    .ethanolVancoillie2012.points;

  it('table reproduces the correlation off the grid (400–800 K, 8–45 bar)', () => {
    const e = pts.map((r) => rel(laminarFlameSpeed(eth, r.phi, r.Tu, r.p, 0), r.SL));
    const st = stats(e);
    console.info(`[flamespeed] ethanol table vs Vancoillie 2012 (${pts.length} off-grid points): mean ${pct(st.mean)}, max ${pct(st.max)}`);
    // the correlation is cubic in p; four ln p nodes (1, 5, 20, 60 bar) reproduce it to ≲ 4 %
    // (its own fit to the mechanism: 7.4 % mean |residual|, Vancoillie et al. Table 6)
    expect(st.mean).toBeLessThan(0.03);
    expect(st.max).toBeLessThan(0.05);
  });

  it('pressure dependence is physical (the Gülder 1984 constants were not): S_L falls with p like other fuels', () => {
    // Gülder (β = −0.17, α = 1.75) gave S_L(373 K, 10 bar) ≈ S_L(300 K, 1 bar) and 1.45 × iso-octane.
    for (const phi of [0.8, 1.0, 1.2]) {
      const beta = Math.log(laminarFlameSpeed(eth, phi, 373, 10e5, 0) / laminarFlameSpeed(eth, phi, 373, 1e5, 0)) / Math.log(10);
      // measured ethanol β(φ=1) ≈ −0.26 (Hinton, Stone, Cracknell & Olm, Oxford bomb correlation, Table 3,
      // 0.7–17 bar); Bradley et al. (2009) steeper
      expect(beta).toBeLessThan(-0.15);
      expect(beta).toBeGreaterThan(-0.45);
      for (const [Tu, p] of [[373, 10e5], [700, 20e5], [800, 40e5]]) {
        const r = laminarFlameSpeed(eth, phi, Tu, p, 0) / laminarFlameSpeed(pure('IC8H18'), phi, Tu, p, 0);
        expect(r).toBeGreaterThan(0.85); // ethanol burns ~0–20 % faster than iso-octane
        expect(r).toBeLessThan(1.3);
      }
    }
    expect(gulder1984(GULDER_1984.C2H5OH, 1, 373, 10e5) / gulder1984(GULDER_1984.C2H5OH, 1, 300, 1e5)).toBeGreaterThan(0.95);
    // 1 bar (below the correlation's 5–85 bar range; power-law continued). Measured φ = 1, 1 bar
    // (Bradley et al., Konnov et al., Liao et al. as compiled in Vancoillie et al. 2012, Figs. 1–2,
    // read by eye): ≈ 0.36–0.42 m/s at 300 K, ≈ 0.52 m/s at 358 K, ≈ 0.85 m/s at 450 K.
    const at = (T: number): number => laminarFlameSpeed(eth, 1, T, 1e5, 0);
    console.info(`[flamespeed] ethanol φ = 1, 1 bar: ${at(300).toFixed(3)} / ${at(358).toFixed(3)} / ${at(450).toFixed(3)} m/s at 300 / 358 / 450 K`);
    expect(at(300)).toBeGreaterThan(0.31);
    expect(at(300)).toBeLessThan(0.46);
    expect(at(358)).toBeGreaterThan(0.44);
    expect(at(358)).toBeLessThan(0.6);
    expect(at(450)).toBeGreaterThan(0.72);
    expect(at(450)).toBeLessThan(0.98);
  });
});

describe('robustness: monotonicity, smoothness, no NaN', () => {
  const fuels: FuelBlend[] = [
    pure('IC8H18'), pure('NC7H16'), pure('CH4'), pure('C3H8'), pure('C2H5OH'),
    fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }),
  ];

  it('dense random sweep over (and beyond) the domain: finite, ≥ 0', () => {
    let seed = 12345;
    const rnd = (): number => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
    let n = 0;
    for (let i = 0; i < 200000; i++) {
      const f = fuels[i % fuels.length];
      const phi = 0.2 + 2.0 * rnd();
      const Tu = 250 + 750 * rnd();
      const p = Math.exp(Math.log(0.3e5) + (Math.log(150e5) - Math.log(0.3e5)) * rnd());
      const x = 0.4 * rnd();
      const s = laminarFlameSpeed(f, phi, Tu, p, x);
      if (!(Number.isFinite(s) && s >= 0 && s < 10)) n++;
    }
    expect(n).toBe(0);
    for (const bad of [NaN, -1, 0, Infinity]) {
      expect(laminarFlameSpeed(fuels[0], bad, 600, 1e6, 0)).toBe(0);
      expect(laminarFlameSpeed(fuels[0], 1, bad, 1e6, 0)).toBe(0);
      expect(laminarFlameSpeed(fuels[0], 1, 600, bad, 0)).toBe(0);
      expect(Number.isFinite(laminarFlameSpeed(fuels[0], 1, 600, 1e6, bad))).toBe(true);
    }
  });

  it('per-blend cache is keyed by composition: in-place mutation and fresh FuelBlend objects are safe', () => {
    // Review fix: the cache was keyed by the X array's identity, so mutating X returned the old
    // blend's S_L, and a fresh FuelBlend per call rebuilt the 480-node tables every call.
    const X = new Float64Array(NS);
    X[SP.IC8H18] = 0.9;
    X[SP.NC7H16] = 0.1;
    const b: FuelBlend = { label: 'mutable', X };
    const s90 = laminarFlameSpeed(b, 1, 600, 20e5, 0.1);
    X[SP.IC8H18] = 1;
    X[SP.NC7H16] = 0;
    const sIso = laminarFlameSpeed(b, 1, 600, 20e5, 0.1);
    expect(sIso).toBe(laminarFlameSpeed(pure('IC8H18'), 1, 600, 20e5, 0.1));
    expect(Math.abs(sIso / s90 - 1)).toBeGreaterThan(1e-3);
    clearFlameSpeedCache();
    const n = 20000;
    let sink = 0;
    const t0 = performance.now();
    for (let i = 0; i < n; i++) sink += laminarFlameSpeed(fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }), 0.95, 600 + (i % 7), 20e5, 0.05);
    const us = ((performance.now() - t0) * 1e3) / n;
    console.info(`[flamespeed] fresh FuelBlend per call: ${us.toFixed(2)} µs/call (incl. fuelFromSelection; sink ${sink.toFixed(0)})`);
    expect(us).toBeLessThan(10); // a table rebuild costs ≳ 50 µs
  });

  it('S_L rises with Tu, falls with dilution, falls with p (φ ≤ 1.1), and is continuous in φ, Tu, p', () => {
    for (const f of fuels) {
      for (const phi of [0.7, 0.9, 1.0, 1.1, 1.3]) {
        for (const p of [0.5e5, 1e5, 10e5, 40e5, 100e5]) {
          let prev = 0;
          for (let Tu = 300; Tu <= 950; Tu += 10) {
            const s = laminarFlameSpeed(f, phi, Tu, p, 0);
            expect(s).toBeGreaterThan(prev);
            prev = s;
          }
        }
        for (const Tu of [350, 600, 900]) {
          let prev = Infinity;
          for (let x = 0; x <= 0.5; x += 0.01) {
            const s = laminarFlameSpeed(f, phi, Tu, 20e5, x);
            if (prev > 0) expect(s).toBeLessThan(prev);
            else expect(s).toBe(0);
            prev = s;
          }
          if (phi <= 1.1) {
            prev = Infinity;
            for (let lp = Math.log(0.5e5); lp <= Math.log(100e5); lp += 0.05) {
              const s = laminarFlameSpeed(f, phi, Tu, Math.exp(lp), 0);
              expect(s).toBeLessThanOrEqual(prev * (1 + 1e-12));
              prev = s;
            }
          }
        }
      }
      // continuity in φ: no jumps at table nodes or at the limit cut-off (|dS/dφ| stays bounded)
      for (const [Tu, p] of [[300, 1e5], [650, 25e5], [900, 80e5]]) {
        const s: number[] = [];
        for (let phi = 0.4; phi <= 1.8; phi += 0.001) s.push(laminarFlameSpeed(f, phi, Tu, p, 0));
        const smax = Math.max(...s);
        let jump = 0;
        for (let i = 1; i < s.length; i++) jump = Math.max(jump, Math.abs(s[i] - s[i - 1]));
        expect(jump / smax).toBeLessThan(0.05); // per Δφ = 0.001 (steepest: the limit cut-off band)
      }
      // continuity in Tu and p across table nodes (relative change per small step)
      for (const phi of [0.8, 1.0, 1.2]) {
        let prev = laminarFlameSpeed(f, phi, 290, 3e5, 0);
        for (let Tu = 291; Tu <= 950; Tu += 1) {
          const s = laminarFlameSpeed(f, phi, Tu, 3e5, 0);
          expect(Math.abs(s / prev - 1)).toBeLessThan(0.01);
          prev = s;
        }
        prev = laminarFlameSpeed(f, phi, 600, 0.5e5, 0);
        for (let lp = Math.log(0.5e5) + 0.01; lp <= Math.log(100e5); lp += 0.01) {
          const s = laminarFlameSpeed(f, phi, 600, Math.exp(lp), 0);
          expect(Math.abs(s / prev - 1)).toBeLessThan(0.01);
          prev = s;
        }
      }
    }
  });
});

describe('flame thickness, burned gas, transport', () => {
  it('burnedGasState: tabulated equilibrium T_ad vs Cantera HP equilibrium; complete combustion is an upper bound', () => {
    const eT: number[] = [];
    const eS: number[] = [];
    const eC: number[] = [];
    for (const f of gri.concat(prf)) {
      const sp = single(f);
      if (!sp || f.Tad === undefined || f.sigmaEq === undefined) continue;
      if (f.set !== undefined && f.set !== 'grid' && f.set !== 'dilution') continue;
      const bg = burnedGasState(pure(sp), f.phi, f.Tu, f.p, f.ydil);
      eT.push(bg.Tb / f.Tad - 1);
      eS.push(bg.sigma / f.sigmaEq - 1);
      const cc = burnedGasState(pure(sp), f.phi, f.Tu, f.p, f.ydil, 'complete-combustion');
      eC.push(cc.Tb / f.Tad - 1);
      expect(cc.Tb).toBeGreaterThan(f.Tad - 15); // no dissociation → upper bound (thermo data differ slightly)
    }
    const st = stats(eT);
    const ss = stats(eS);
    const sc = stats(eC);
    console.info(`[flamespeed] T_b: table vs Cantera equilibrium mean ${pct(st.bias)}, max ${pct(st.max)}; σ: bias ${pct(ss.bias)}, max ${pct(ss.max)}; complete combustion: bias ${pct(sc.bias)}, max ${pct(sc.max)}`);
    expect(st.max).toBeLessThan(0.015); // air with/without Ar, GRI/PRF_HT vs NASA thermo, diluent definition
    expect(ss.max).toBeLessThan(0.04);
    expect(sc.max).toBeLessThan(0.15);
  });

  it('burnedGasState accepts the equilibrium module\'s EquilibriumSolver (HPEquilibriumSolver contract)', () => {
    const solver = new EquilibriumSolver();
    const e: number[] = [];
    for (const f of gri.concat(prf)) {
      const sp = single(f);
      if (!sp || f.Tad === undefined || (f.set !== undefined && f.set !== 'grid')) continue;
      if (f.Tu !== 500 && f.Tu !== 900) continue;
      const bg = burnedGasState(pure(sp), f.phi, f.Tu, f.p, f.ydil, solver);
      expect(bg.method).toBe('solver');
      e.push(bg.Tb / f.Tad - 1);
    }
    expect(e.length).toBeGreaterThan(100);
    const st = stats(e);
    console.info(`[flamespeed] T_b with EquilibriumSolver vs Cantera HP equilibrium (${e.length}): bias ${pct(st.bias)}, max ${pct(st.max)}`);
    expect(st.max).toBeLessThan(0.01); // dry air (Ar, CO2) vs O2 + 3.76 N2, NASA vs mechanism thermo
  });

  it('diffusive and Blint thickness vs the Cantera thermal thickness (T_b − T_u)/max|dT/dx|', () => {
    const eD: number[] = [];
    const eB: number[] = [];
    const eA: number[] = [];
    for (const f of gri.concat(prf)) {
      const sp = single(f);
      if (!sp || f.deltaT === undefined || f.alphaU === undefined || f.Tb === undefined) continue;
      const th = flameThickness(pure(sp), f.phi, f.Tu, f.p, f.SL, f.ydil, f.Tb);
      eA.push(rel(th.alphaU, f.alphaU));
      eD.push(rel(th.diffusive, f.deltaT));
      eB.push(rel(th.blint, f.deltaT));
    }
    const sa = stats(eA);
    const sb = stats(eB);
    console.info(`[flamespeed] α_u vs Cantera: max ${pct(sa.max)}; Blint vs thermal thickness: mean |e| ${pct(sb.mean)}, bias ${pct(sb.bias)}, max ${pct(sb.max)}; diffusive: bias ${pct(stats(eD).bias)}`);
    expect(sa.max).toBeLessThan(0.05); // different transport data (GRI/PRF_HT vs ours)
    expect(sb.mean).toBeLessThan(0.25);
    expect(stats(eD).bias).toBeLessThan(-0.7); // δ_D ≈ δ_T/5–10: the known under-estimate
  });

  it('Chapman–Enskog D_ij and mixture-averaged D_k,m match Cantera (same LJ data) to < 1.5 %', () => {
    let worst = 0;
    for (const r of fitFx.transport) {
      const f = pure(r.fuel as Species);
      const X = unburnedComposition(f, r.phi, 0, new Float64Array(NS));
      const k = SP[r.fuel as Species];
      worst = Math.max(worst, Math.abs(rel(binaryDiffusionCoefficient(k, SP.N2, r.Tu, r.p), r.DfuelN2)));
      worst = Math.max(worst, Math.abs(rel(binaryDiffusionCoefficient(SP.O2, SP.N2, r.Tu, r.p), r.DO2N2)));
      worst = Math.max(worst, Math.abs(rel(mixtureDiffusionCoefficient(X, k, r.Tu, r.p), r.Dfuel)));
      worst = Math.max(worst, Math.abs(rel(mixtureDiffusionCoefficient(X, SP.O2, r.Tu, r.p), r.DO2)));
      const le = lewisNumbers(f, r.phi, r.Tu, r.p);
      expect(Math.abs(rel(le.alphaU, r.alpha))).toBeLessThan(0.01);
    }
    console.info(`[flamespeed] diffusion coefficients vs Cantera: max ${pct(worst)}`);
    expect(worst).toBeLessThan(0.015);
  });

  it('Lewis numbers: CH4 ≈ 1, heavy fuels ≫ 1 (lean), O2 ≈ 1', () => {
    const ch4 = lewisNumbers(pure('CH4'), 0.8, 300, 1e5);
    const ic8 = lewisNumbers(pure('IC8H18'), 0.8, 373, 1e5);
    expect(ch4.fuel).toBeGreaterThan(0.9);
    expect(ch4.fuel).toBeLessThan(1.1);
    expect(ic8.fuel).toBeGreaterThan(2.5);
    expect(ic8.oxygen).toBeGreaterThan(0.95);
    expect(ic8.oxygen).toBeLessThan(1.2);
  });
});

describe('Markstein length (Bechtold–Matalon 2001)', () => {
  it('integrals: closed form J1, and J2 equals the alternative λ ≡ 1 form ∫₁^σ ln ξ/(ξ−1) dξ', () => {
    const quad = (fn: (x: number) => number, a: number, b: number, n = 200000): number => {
      let s = 0;
      const h = (b - a) / n;
      for (let i = 0; i < n; i++) s += fn(a + (i + 0.5) * h);
      return s * h;
    };
    for (const sigma of [1.5, 3, 7.5]) {
      const [j1, j2] = marksteinIntegrals(sigma, 0);
      expect(j1).toBeCloseTo(Math.log(sigma), 12);
      expect(Math.abs(j2 / quad((x) => Math.log(x) / (x - 1), 1, sigma) - 1)).toBeLessThan(1e-4);
      const [k1, k2] = marksteinIntegrals(sigma, 0.7);
      expect(Math.abs(k1 / quad((x) => Math.pow(x, -0.3), 1, sigma) - 1)).toBeLessThan(1e-6);
      expect(Math.abs(k2 / quad((x) => Math.log((sigma - 1) / (x - 1)) * Math.pow(x, -0.3), 1, sigma) - 1)).toBeLessThan(1e-3);
    }
    // σ → 1 limit: J1/(σ−1) → 1, J2/(σ−1) → 1  ⇒  Ma → 1 + β(Le−1)/2 (Matalon 2011)
    const [a1, a2] = marksteinIntegrals(1 + 1e-6, 0.7);
    expect(a1 / 1e-6).toBeCloseTo(1, 5);
    expect(a2 / 1e-6).toBeCloseTo(1, 4);
  });

  it('effective Lewis number limits', () => {
    expect(effectiveLewisNumber(2.8, 1.0, 1, 10)).toBeCloseTo(1.9, 12);
    expect(effectiveLewisNumber(2.8, 1.0, 0.5, 10)).toBeGreaterThan(2.6);
    expect(effectiveLewisNumber(2.8, 1.0, 1.6, 10)).toBeLessThan(1.25);
    expect(effectiveLewisNumber(2.8, 1.0, 3.0, 10)).toBeLessThan(1.1);
  });

  it('trends: heavy-fuel lean flames strongly stable, rich less; methane lean small; values in the measured range', () => {
    const iso = pure('IC8H18');
    const lean = marksteinLengths(iso, 0.8, 358, 1e5);
    const stoich = marksteinLengths(iso, 1.0, 358, 1e5);
    const rich = marksteinLengths(iso, 1.4, 358, 1e5);
    console.info(
      `[flamespeed] iso-octane 358 K, 1 bar: L_b = ${(1e3 * lean.burned).toFixed(2)} / ${(1e3 * stoich.burned).toFixed(2)} / ${(1e3 * rich.burned).toFixed(2)} mm at φ = 0.8 / 1.0 / 1.4 (β = ${stoich.beta.toFixed(1)}, Le_eff = ${stoich.leEff.toFixed(2)}, σ = ${stoich.sigma.toFixed(2)})`,
    );
    // Markstein NUMBER falls with φ (Le_eff: fuel-controlled lean, O2-controlled rich)
    expect(lean.numberBurned).toBeGreaterThan(stoich.numberBurned);
    expect(stoich.numberBurned).toBeGreaterThan(rich.numberBurned);
    expect(stoich.burned).toBeGreaterThan(0.3e-3);
    expect(stoich.burned).toBeLessThan(3e-3);
    // pressure reduces L_b (thinner flames)
    const hp = marksteinLengths(iso, 1.0, 358, 10e5);
    expect(hp.burned).toBeLessThan(stoich.burned);
    // Order-of-magnitude check against spherical-bomb data for stoichiometric iso-octane at
    // 10 bar (read approximately from the raster Fig. 4 of Jerzembeck et al. 2009): Bradley et al.
    // (1998, 358 K) S_b ≈ 1.81 → 1.74 m/s over κ = 0 → 360 1/s, i.e. L_b ≈ 0.19 mm; Jerzembeck
    // (373 K) ≈ 0. Asymptotic theory is expected at the upper end (Matalon 2011, slide 16).
    console.info(`[flamespeed] iso-octane φ = 1, 358 K, 10 bar: L_b = ${(1e3 * hp.burned).toFixed(3)} mm (Bradley et al. 1998 ≈ 0.19 mm)`);
    expect(hp.burned).toBeGreaterThan(0.19e-3 / 3);
    expect(hp.burned).toBeLessThan(0.19e-3 * 3);
    // methane: lean L_b < rich L_b (Le_F < 1 lean)
    const ch4 = pure('CH4');
    expect(marksteinLength(ch4, 0.7, 300, 1e5)).toBeLessThan(marksteinLength(ch4, 1.3, 300, 1e5));
    expect(ACTIVATION_ENERGY_TABLE.every((r) => r.every((e) => e > 1e5 && e < 4e5))).toBe(true);
  });
});

describe('performance', () => {
  it('laminarFlameSpeed ≲ 200 ns/call (PRF blend, φ fixed over the cycle, Tu/p/x varying)', () => {
    const f = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    const n = 300000;
    let sink = 0;
    const run = (phiVaries: boolean): number => {
      const t0 = performance.now();
      for (let i = 0; i < n; i++) {
        const phi = phiVaries ? 0.8 + (i % 50) * 0.01 : 0.95;
        sink += laminarFlameSpeed(f, phi, 500 + (i % 37) * 10, 1e6 + 37 * i, 0.05 + 1e-7 * (i % 100));
      }
      return ((performance.now() - t0) * 1e6) / n;
    };
    const ref = (): number => {
      const t0 = performance.now();
      for (let i = 0; i < n; i++) sink += Math.exp(1e-6 * i) + Math.log(1 + i);
      return ((performance.now() - t0) * 1e6) / n;
    };
    let fixed = Infinity;
    let varying = Infinity;
    let r = Infinity;
    for (let rep = 0; rep < 5; rep++) {
      fixed = Math.min(fixed, run(false));
      varying = Math.min(varying, run(true));
      r = Math.min(r, ref());
    }
    console.info(
      `[flamespeed] laminarFlameSpeed: ${fixed.toFixed(0)} ns/call at fixed φ (${(fixed / r).toFixed(1)} × exp+log), ` +
        `${varying.toFixed(0)} ns when φ changes every call; exp+log ${r.toFixed(1)} ns (sink ${sink.toFixed(0)})`,
    );
    expect(fixed / r).toBeLessThan(40); // ≈ 200 ns on an idle machine where exp + log ≈ 5–10 ns
    expect(fixed).toBeLessThan(2000);
  });
});
