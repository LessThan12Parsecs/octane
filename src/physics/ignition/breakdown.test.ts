import { describe, expect, it } from 'vitest';
import gasFixture from '../../../test/fixtures/ignition_kernel_gas.json';
import {
  AIR_UNIFORM_EC,
  AIR_UNIFORM_K,
  breakdownVoltage,
  breakdownVoltagePaschen,
  breakdownVoltagePashley,
  breakdownVoltageUniformAir,
  normalCathodeFall,
  PASCHEN_MINIMUM_AIR,
  paschenMinimumVoltage,
  PASHLEY_INTERCEPT_KV,
  reducedPressureBar,
  relativeAirDensity,
  TOWNSEND_A_AIR,
  TOWNSEND_B_AIR,
} from './breakdown';

const MM = 1e-3;
const ATM = 101325;
const TORR = 133.322368;

describe('uniform-field air breakdown (Schumann-type, Kuffel et al. 2000)', () => {
  it('reproduces the closed form at 1 bar, 20 °C, 1 mm (24.36·0.1 + 6.72·√0.1 kV)', () => {
    const v = breakdownVoltageUniformAir(1 * MM, 1e5, 293.15);
    expect(v).toBeCloseTo(1e3 * (24.36 * 0.1 + 6.72 * Math.sqrt(0.1)), 6);
    expect(v).toBeGreaterThan(4.4e3);
    expect(v).toBeLessThan(4.7e3);
  });

  it('uses the pressure in BAR reduced to 20 °C (Kuffel: "p in bar, d in cm"), not the 1-atm relative density', () => {
    // regression: an earlier revision evaluated the law with δ = (p/1 atm)(T₀/T)
    expect(reducedPressureBar(ATM, 293.15)).toBeCloseTo(1.01325, 12);
    expect(reducedPressureBar(ATM, 293.15) / relativeAirDensity(ATM, 293.15)).toBeCloseTo(1.01325, 12);
    const pd = 1.01325 * 0.1; // bar·cm at 1 atm, 1 mm
    expect(breakdownVoltageUniformAir(1 * MM, ATM, 293.15)).toBeCloseTo(1e3 * (24.36 * pd + 6.72 * Math.sqrt(pd)), 6);
  });

  it('is a few kV/mm at 1 atm, 300 K for engine-size gaps', () => {
    for (const d of [0.5, 0.7, 0.9, 1.1]) {
      const perMm = breakdownVoltageUniformAir(d * MM, ATM, 300) / d;
      expect(perMm).toBeGreaterThan(2.5e3);
      expect(perMm).toBeLessThan(6e3);
    }
  });

  it('lies 3–7 % above the independent Schumann-form set 24.22/6.08 kV (δ at 1 atm; J. Li, HV lecture 4-2)', () => {
    for (const dd of [0.05, 0.1, 0.2, 0.5, 1.0]) {
      // δd in cm → choose d = dd cm at δ = 1 (1 atm, 20 °C)
      const v = breakdownVoltageUniformAir(dd * 1e-2, ATM, 293.15);
      const li = 1e3 * (24.22 * dd + 6.08 * Math.sqrt(dd));
      expect(v / li - 1).toBeGreaterThan(0.03);
      expect(v / li - 1).toBeLessThan(dd >= 0.5 ? 0.05 : 0.07);
    }
  });

  it('depends on p and T only through the density (N·d scaling)', () => {
    const a = breakdownVoltageUniformAir(0.6 * MM, 12e5, 650);
    const b = breakdownVoltageUniformAir(0.6 * MM, 24e5, 1300);
    expect(a).toBeCloseTo(b, 8);
    const c = breakdownVoltageUniformAir(1.2 * MM, 6e5, 650); // same N·d
    expect(c).toBeCloseTo(a, 8);
  });

  it('is nearly linear in N·d at engine densities (slope → E_c)', () => {
    const d = 0.508 * MM;
    const T = 650;
    const v1 = breakdownVoltageUniformAir(d, 10e5, T);
    const v2 = breakdownVoltageUniformAir(d, 20e5, T);
    const v4 = breakdownVoltageUniformAir(d, 40e5, T);
    expect(v2 / v1).toBeGreaterThan(1.75);
    expect(v2 / v1).toBeLessThan(2.0);
    expect(v4 / v2).toBeGreaterThan(1.8);
    // local slope dV/d(δd) = E_c + K/(2√(δd)) → E_c (the critical reduced field) as N·d grows
    const dd = reducedPressureBar(40e5, T) * d;
    const slope = (breakdownVoltageUniformAir(d * 1.000001, 40e5, T) - v4) / (dd * 0.000001);
    expect(slope / (AIR_UNIFORM_EC + AIR_UNIFORM_K / (2 * Math.sqrt(dd)))).toBeCloseTo(1, 5);
    expect(slope / AIR_UNIFORM_EC).toBeLessThan(1.2);
  });
});

describe('Pashley, Stone & Roberts (2000) engine correlation', () => {
  it('V_b[kV] = 4.3 + 136 p/T + 324 p d/T (COBEM 2011 eq. 2; intercept "4,3", not 3.4)', () => {
    expect(PASHLEY_INTERCEPT_KV).toBe(4.3);
    // 1 bar, 300 K, 1 mm: 4.3 + 0.4533 + 1.08 = 5.833 kV
    expect(breakdownVoltagePashley(1 * MM, 1e5, 300)).toBeCloseTo(1e3 * (4.3 + 136 / 300 + 324 / 300), 6);
    // zero density → intercept only
    expect(breakdownVoltagePashley(1 * MM, 0, 300)).toBeCloseTo(4300, 9);
    // linear in p/T and in p·d/T
    const a = breakdownVoltagePashley(0.8 * MM, 10e5, 600);
    const b = breakdownVoltagePashley(0.8 * MM, 20e5, 1200);
    expect(a).toBeCloseTo(b, 9);
  });
});

describe('Paschen/Townsend law (Cobine A = 14.6 /(cm Torr), B = 365 V/(cm Torr))', () => {
  it('has the analytic minimum V_min = e(B/A)ln(1+1/γ) at pd = e ln(1+1/γ)/A', () => {
    const g = 0.02;
    const lnG = Math.log(1 + 1 / g);
    const vMin = Math.E * (TOWNSEND_B_AIR / TOWNSEND_A_AIR) * lnG;
    const pdMin = (Math.E * lnG) / TOWNSEND_A_AIR; // Pa m
    // scan pd around the minimum at the reference temperature
    let best = Infinity;
    let bestPd = 0;
    for (let f = 0.3; f < 3; f += 0.001) {
      const pd = pdMin * f;
      const v = breakdownVoltagePaschen(1e-3, pd / 1e-3, 293.15, g);
      if (v < best) {
        best = v;
        bestPd = pd;
      }
    }
    expect(best).toBeCloseTo(vMin, 1);
    expect(bestPd / pdMin).toBeCloseTo(1, 2);
    // ≈ 267 V at ≈ 0.73 Torr·cm for air on iron (measured air minimum ≈ 330 V)
    expect(vMin).toBeGreaterThan(250);
    expect(vMin).toBeLessThan(350);
    expect(pdMin / (TORR * 1e-2)).toBeCloseTo(0.73, 1);
  });

  it('returns +Infinity left of the Paschen asymptote', () => {
    expect(breakdownVoltagePaschen(1e-6, 1e3, 300)).toBe(Infinity);
  });

  it('normal cathode fall 3(B/A)ln(1+1/γ) = 295 V (air, iron, γ = 0.02)', () => {
    expect(normalCathodeFall()).toBeCloseTo(295, 0);
  });
});

describe('breakdown voltage at CFR-like spark conditions (Cantera-compressed charge)', () => {
  const cases = gasFixture.cases;
  it('is ≈ 9–25 kV for 0.5–0.9 mm gaps at spark timing, CR 7–10; Paschen within 30 %, Pashley within 40 %', () => {
    const rows: string[] = [];
    for (const c of cases.filter((k) => k.phi === 1)) {
      for (const d of [0.508, 0.9]) {
        const vu = breakdownVoltage(d * MM, c.p, c.Tu);
        const vp = breakdownVoltage(d * MM, c.p, c.Tu, undefined, { model: 'paschen' });
        const vs = breakdownVoltagePashley(d * MM, c.p, c.Tu);
        rows.push(`CR ${c.cr} d ${d} mm: uniform ${(vu / 1e3).toFixed(2)} kV, Paschen ${(vp / 1e3).toFixed(2)}, Pashley ${(vs / 1e3).toFixed(2)}`);
      }
    }
    if (process.env.IGNITION_REPORT) console.log(rows.join('\n'));
    for (const c of cases.filter((k) => k.phi === 1)) {
      for (const d of [0.508, 0.9]) {
        const vu = breakdownVoltage(d * MM, c.p, c.Tu);
        const vp = breakdownVoltage(d * MM, c.p, c.Tu, undefined, { model: 'paschen' });
        const vs = breakdownVoltagePashley(d * MM, c.p, c.Tu);
        expect(vu).toBeGreaterThan(8.5e3);
        expect(vu).toBeLessThan(25e3);
        expect(Math.abs(vp / vu - 1)).toBeLessThan(0.3);
        expect(Math.abs(vs / vu - 1)).toBeLessThan(0.4);
      }
    }
  });

  it('applies the impulse (dynamic over-voltage) factor multiplicatively', () => {
    const c = cases[0];
    const v = breakdownVoltage(0.5 * MM, c.p, c.Tu);
    expect(breakdownVoltage(0.5 * MM, c.p, c.Tu, undefined, { impulseFactor: 1.32 })).toBeCloseTo(1.32 * v, 6);
  });
});

describe('breakdown voltage floor (reviewer fix)', () => {
  it('never falls below the Paschen minimum e(B/A)ln(1+1/γ) ≈ 267 V, whatever the model', () => {
    expect(PASCHEN_MINIMUM_AIR).toBeCloseTo(paschenMinimumVoltage(0.02), 12);
    expect(PASCHEN_MINIMUM_AIR).toBeGreaterThan(260);
    expect(PASCHEN_MINIMUM_AIR).toBeLessThan(275);
    // hot kernel gas at low pressure: the uniform-field fit alone gives ≈ 0.24 kV
    const raw = breakdownVoltageUniformAir(0.2 * MM, 0.5e5, 2500);
    expect(raw).toBeLessThan(PASCHEN_MINIMUM_AIR);
    for (const model of ['uniform-air', 'paschen', 'pashley'] as const) {
      const v = breakdownVoltage(0.2 * MM, 0.5e5, 2500, undefined, { model });
      expect(v).toBeGreaterThanOrEqual(PASCHEN_MINIMUM_AIR);
      expect(Number.isFinite(v) || model === 'paschen').toBe(true);
    }
    expect(breakdownVoltage(0.2 * MM, 0, 300)).toBe(PASCHEN_MINIMUM_AIR);
    // engine densities are unaffected
    expect(breakdownVoltage(0.508 * MM, 10.85e5, 612)).toBe(breakdownVoltageUniformAir(0.508 * MM, 10.85e5, 612));
  });
});
