import { describe, expect, it } from 'vitest';
import { EL, NE, NS, SP } from '../core/species';
import fuelFx from '../../../test/fixtures/thermo_fuels.json';
import {
  completeCombustionMoles,
  completeCombustionProducts,
  dryAir,
  freshCharge,
  fuelFromSelection,
  fuelMolarMass,
  HEAT_OF_VAPORIZATION_298,
  humidAir,
  ISOOCTANE_DENSITY_60F,
  isooctaneLiquidDensity,
  lowerHeatingValue,
  lowerHeatingValueLiquid,
  lowerHeatingValueMolar,
  moistAirEnhancementFactor,
  NHEPTANE_DENSITY_60F,
  prfIsooctaneMoleFraction,
  stoichAirFuelRatio,
  stoichO2PerMolFuel,
  waterGasShiftK,
  waterSaturationPressure,
} from './fuels';
import { elementMoles, mixMolarMass, totalMass } from './mixture';
import { MOLAR_MASS } from './thermo';

const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);
const iso = fuelFromSelection({ kind: 'pure', species: 'IC8H18' });
const hep = fuelFromSelection({ kind: 'pure', species: 'NC7H16' });

function expectElementsConserved(Nin: Float64Array, Nout: Float64Array, tol = 1e-13): void {
  const a = elementMoles(Nin);
  const b = elementMoles(Nout);
  for (let e = 0; e < NE; e++) expect(Math.abs(a[e] - b[e])).toBeLessThanOrEqual(tol * Math.max(1, Math.abs(a[e])));
}

describe('fuel selection and PRF conversion', () => {
  it('NIST SRM 2214 iso-octane density equation reproduces its certified values', () => {
    expect(isooctaneLiquidDensity(15)).toBeCloseTo(695.969, 2);
    expect(isooctaneLiquidDensity(20)).toBeCloseTo(691.872, 3);
    expect(isooctaneLiquidDensity(25)).toBeCloseTo(687.753, 2);
    expect(ISOOCTANE_DENSITY_60F).toBeCloseTo(695.515, 3);
    expect(NHEPTANE_DENSITY_60F).toBe(687.597);
  });

  it('converts liquid-volume % to vapour mole fraction', () => {
    expect(prfIsooctaneMoleFraction(100)).toBe(1);
    expect(prfIsooctaneMoleFraction(0)).toBe(0);
    // Hand calculation for PRF 90: n_i = 0.9·695.515/0.114232 = 5479.74 mol/m³,
    // n_h = 0.1·687.597/0.100205 = 686.19 mol/m³  →  x_i = 0.88870
    expect(prfIsooctaneMoleFraction(90)).toBeCloseTo(0.8887, 4);
    // iso-octane is denser but heavier per molecule: x_i < volume fraction mid-range
    expect(prfIsooctaneMoleFraction(50)).toBeLessThan(0.5);
    const f = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    expect(f.label).toBe('PRF 90');
    expect(f.octaneNumber).toBe(90);
    expect(f.X[SP.IC8H18] + f.X[SP.NC7H16]).toBeCloseTo(1, 15);
    expect(() => fuelFromSelection({ kind: 'PRF', octaneNumber: 101 })).toThrow(RangeError);
    expect(iso.octaneNumber).toBe(100);
    expect(hep.octaneNumber).toBe(0);
    expect(fuelFromSelection({ kind: 'pure', species: 'CH4' }).octaneNumber).toBeUndefined();
    expect(fuelMolarMass(iso)).toBe(MOLAR_MASS[SP.IC8H18]);
  });
});

describe('air', () => {
  it('dry air = CIPM-2007 composition (renormalised), M ≈ 28.965 g/mol', () => {
    const a = dryAir();
    for (let k = 0; k < NS; k++) expect(a[k]).toBeCloseTo(fuelFx.dryAirX[k], 15);
    expect(rel(mixMolarMass(a), 28.96546e-3)).toBeLessThan(2e-5);
  });

  it('Wagner–Pruss saturation pressure reproduces the published check values', () => {
    // Wagner & Pruss (1993) Table 1
    expect(rel(waterSaturationPressure(273.16), 611.657)).toBeLessThan(1e-6);
    expect(rel(waterSaturationPressure(373.1243), 0.101325e6)).toBeLessThan(1e-6);
    expect(rel(waterSaturationPressure(647.096), 22.064e6)).toBeLessThan(1e-9);
    // Murphy & Koop (2005) supercooled branch joins to ~7e-8
    expect(rel(waterSaturationPressure(273.159999), waterSaturationPressure(273.16))).toBeLessThan(1e-6);
    // monotone across the branch switch and below
    let prev = 0;
    for (let T = 230; T < 400; T += 0.5) {
      const p = waterSaturationPressure(T);
      expect(p).toBeGreaterThan(prev);
      prev = p;
    }
  });

  it('humid air: x_H2O = RH f p_sat/p, rest dry air', () => {
    const T = 298.15;
    const p = 101325;
    const h = humidAir(T, p, 0.6);
    const f = moistAirEnhancementFactor(T, p);
    expect(f).toBeCloseTo(1.00062 + 3.14e-8 * p + 5.6e-7 * 25 * 25, 12);
    const xv = (0.6 * f * waterSaturationPressure(T)) / p;
    expect(h[SP.H2O]).toBeCloseTo(xv, 15);
    // steam tables: p_sat(25 °C) = 3.1699 kPa → x_v = 0.6 · 1.00415 · 3169.9 / 101325 = 0.018848
    expect(waterSaturationPressure(T)).toBeCloseTo(3169.9, 0);
    expect(xv).toBeCloseTo(0.018848, 5);
    expect(h.reduce((s, x) => s + x, 0)).toBeCloseTo(1, 15);
    expect(h[SP.O2] / h[SP.N2]).toBeCloseTo(dryAir()[SP.O2] / dryAir()[SP.N2], 14);
    expect(humidAir(T, p, 0)[SP.H2O]).toBe(0);
  });
});

describe('stoichiometry', () => {
  it('O2 demand and AFR', () => {
    expect(stoichO2PerMolFuel(iso)).toBe(12.5);
    expect(stoichO2PerMolFuel(hep)).toBe(11);
    expect(stoichO2PerMolFuel(fuelFromSelection({ kind: 'pure', species: 'C2H5OH' }))).toBe(3);
    const prf = fuelFromSelection({ kind: 'PRF', octaneNumber: 80 });
    const x = prf.X[SP.IC8H18];
    expect(stoichO2PerMolFuel(prf)).toBeCloseTo(12.5 * x + 11 * (1 - x), 13);
    const afr = stoichAirFuelRatio(iso, dryAir());
    expect(rel(afr, fuelFx.stoichAfrIsooctaneDryAir)).toBeLessThan(1e-12);
    expect(afr).toBeCloseTo(15.1, 1);
    console.info(`[fuels] stoichiometric AFR iso-octane/dry air = ${afr.toFixed(4)}`);
  });

  it('freshCharge: fuel/O2 ratio = φ/ν and EGR mass fraction honoured', () => {
    const air = humidAir(313.15, 101325, 0.4);
    const prf = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    const nu = stoichO2PerMolFuel(prf);
    for (const phi of [0.6, 1.0, 1.3]) {
      const X = freshCharge({ fuel: prf, phi, airX: air });
      const fuelMol = X[SP.IC8H18] + X[SP.NC7H16];
      expect(fuelMol / X[SP.O2]).toBeCloseTo(phi / nu, 14);
      expect(X.reduce((s, v) => s + v, 0)).toBeCloseTo(1, 15);
      expect(X[SP.H2O] / X[SP.N2]).toBeCloseTo(air[SP.H2O] / air[SP.N2], 14);
    }
    // EGR (default composition: complete-combustion products of the fuel–air mixture)
    const egr = 0.15;
    const Xe = freshCharge({ fuel: prf, phi: 1, airX: air, egrFraction: egr });
    // mass of fuel (only in the fresh part) / total mass = (1 − egr) × fuel mass fraction of fresh part
    const Xf = freshCharge({ fuel: prf, phi: 1, airX: air });
    const yFuel = (X: Float64Array): number =>
      (X[SP.IC8H18] * MOLAR_MASS[SP.IC8H18] + X[SP.NC7H16] * MOLAR_MASS[SP.NC7H16]) / totalMass(X);
    expect(yFuel(Xe)).toBeCloseTo((1 - egr) * yFuel(Xf), 14);
    expect(Xe[SP.CO2]).toBeGreaterThan(Xf[SP.CO2]);
    // explicit EGR composition
    const egrX = new Float64Array(NS);
    egrX[SP.N2] = 1;
    const Xn = freshCharge({ fuel: prf, phi: 1, airX: air, egrFraction: 0.2, egrX });
    expect(yFuel(Xn)).toBeCloseTo(0.8 * yFuel(Xf), 14);
    expect(() => freshCharge({ fuel: prf, phi: 1, airX: air, egrFraction: 1 })).toThrow(RangeError);
  });

  it('freshCharge extremes (φ 0…100, EGR up to 99.9 %) stay normalised and non-negative', () => {
    const prf = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    for (const phi of [0, 1e-9, 0.2, 3, 10, 100]) {
      for (const egrFraction of [0, 0.5, 0.999]) {
        const X = freshCharge({ fuel: prf, phi, airX: dryAir(), egrFraction });
        let s = 0;
        for (let k = 0; k < NS; k++) {
          expect(X[k]).toBeGreaterThanOrEqual(0);
          s += X[k];
        }
        expect(s).toBeCloseTo(1, 14);
      }
    }
    expect(freshCharge({ fuel: prf, phi: 0, airX: dryAir() })[SP.IC8H18]).toBe(0);
  });
});

describe('complete-combustion products', () => {
  const air = dryAir();

  it('lean and stoichiometric: CO2/H2O/O2/N2, elements conserved', () => {
    for (const phi of [0.3, 0.8, 1.0]) {
      const R = freshCharge({ fuel: iso, phi, airX: air });
      const P = completeCombustionMoles(R);
      expectElementsConserved(R, P);
      expect(P[SP.CO]).toBe(0);
      expect(P[SP.H2]).toBe(0);
      expect(P[SP.IC8H18]).toBe(0);
      const nf = R[SP.IC8H18];
      expect(P[SP.CO2]).toBeCloseTo(8 * nf + R[SP.CO2], 15);
      expect(P[SP.H2O]).toBeCloseTo(9 * nf, 15);
      expect(P[SP.O2]).toBeCloseTo(R[SP.O2] - 12.5 * nf, 14);
    }
  });

  it('rich: water-gas-shift equilibrium at 1740 K matches Cantera', () => {
    let worst = 0;
    for (const c of fuelFx.richProducts1740) {
      const R = Float64Array.from(c.reactantsX);
      const P = completeCombustionProducts(R);
      for (let k = 0; k < NS; k++) {
        worst = Math.max(worst, Math.abs(P[k] - c.productsX[k]));
        expect(Math.abs(P[k] - c.productsX[k]), `phi ${c.phi} k ${k}`).toBeLessThan(1e-9);
      }
      expectElementsConserved(R, completeCombustionMoles(R), 1e-13);
    }
    console.info(`[fuels] rich WGS products vs Cantera: max |ΔX| = ${worst.toExponential(2)}`);
  });

  it('beyond the CO limit burns only part of the fuel, conserving elements', () => {
    for (const phi of [3.5, 6]) {
      for (const fuel of [iso, fuelFromSelection({ kind: 'pure', species: 'C2H5OH' })]) {
        const R = freshCharge({ fuel, phi, airX: air });
        const P = completeCombustionMoles(R);
        expectElementsConserved(R, P, 1e-12);
        for (let k = 0; k < NS; k++) expect(P[k]).toBeGreaterThanOrEqual(0);
        expect(P[SP.O2]).toBe(0);
      }
    }
    // ethanol stays below its CO limit up to high φ; check a genuinely partial case
    const R = freshCharge({ fuel: iso, phi: 5, airX: air });
    const P = completeCombustionMoles(R);
    expect(P[SP.IC8H18]).toBeGreaterThan(0);
    expect(P[SP.IC8H18]).toBeLessThan(R[SP.IC8H18]);
  });

  it('re-allocates residual species (CO, H2, OH, NO …) in the reactants', () => {
    const R = new Float64Array(NS);
    R[SP.N2] = 0.7;
    R[SP.O2] = 0.2;
    R[SP.CO] = 0.01;
    R[SP.H2] = 0.005;
    R[SP.OH] = 0.001;
    R[SP.NO] = 0.002;
    R[SP.NC7H16] = 0.01;
    const P = completeCombustionMoles(R);
    expectElementsConserved(R, P);
    expect(P[SP.OH]).toBe(0);
    expect(P[SP.NO]).toBe(0);
  });

  it('seeded random reactant sweep: elements conserved, products finite and non-negative', () => {
    let seed = 777;
    const rnd = (): number => {
      seed = (seed * 1103515245 + 12345) & 0x7fffffff;
      return seed / 0x7fffffff;
    };
    const N = new Float64Array(NS);
    const P = new Float64Array(NS);
    for (let trial = 0; trial < 5000; trial++) {
      N.fill(0);
      const nsp = 1 + Math.floor(rnd() * NS);
      for (let j = 0; j < nsp; j++) N[Math.floor(rnd() * NS)] = Math.pow(10, -6 * rnd());
      completeCombustionMoles(N, P, 500 + 3000 * rnd());
      expectElementsConserved(N, P, 1e-11);
      for (let k = 0; k < NS; k++) expect(Number.isFinite(P[k]) && P[k] >= 0, `trial ${trial} k ${k}`).toBe(true);
    }
  });

  it('water-gas-shift constant vs Cantera', () => {
    for (const w of fuelFx.wgsK) expect(rel(waterGasShiftK(w.T), w.K)).toBeLessThan(1e-12);
    expect(waterGasShiftK(1740)).toBeCloseTo(3.5, 0); // Heywood (1988): K ≈ 3.5 near 1740 K
  });
});

describe('heating values', () => {
  it('gaseous LHV equals Cantera ΔH of complete combustion (298.15 K)', () => {
    const names = ['IC8H18', 'NC7H16', 'CH4', 'C3H8', 'C2H5OH'] as const;
    for (const n of names) {
      const f = fuelFromSelection({ kind: 'pure', species: n });
      const ref = fuelFx.lhvGas[n];
      expect(rel(lowerHeatingValueMolar(f), ref.molar), n).toBeLessThan(1e-9);
      expect(rel(lowerHeatingValue(f), ref.mass), n).toBeLessThan(1e-9);
    }
  });

  it('liquid-basis LHV of iso-octane ≈ 44.3 MJ/kg (Heywood 1988, App. D)', () => {
    const lhvL = lowerHeatingValueLiquid(iso);
    console.info(
      `[fuels] iso-octane LHV gas ${(lowerHeatingValue(iso) / 1e6).toFixed(3)} MJ/kg, liquid ${(lhvL / 1e6).toFixed(3)} MJ/kg; ` +
        `n-heptane liquid ${(lowerHeatingValueLiquid(hep) / 1e6).toFixed(3)} MJ/kg; ` +
        `ethanol liquid ${(lowerHeatingValueLiquid(fuelFromSelection({ kind: 'pure', species: 'C2H5OH' })) / 1e6).toFixed(3)} MJ/kg`,
    );
    expect(Math.abs(lhvL - 44.3e6)).toBeLessThan(0.1e6);
    expect(HEAT_OF_VAPORIZATION_298.IC8H18).toBeCloseTo(35.14e3, -1);
    expect(HEAT_OF_VAPORIZATION_298.NC7H16).toBeCloseTo(36.57e3, -1);
    // propane: liquid → IDEAL gas at 298.15 K; NIST WebBook EOS h_ig − h_liq,sat = 16.23 kJ/mol
    expect(Math.abs(HEAT_OF_VAPORIZATION_298.C3H8 - 16.23e3)).toBeLessThan(0.05e3);
    expect(Number.isNaN(lowerHeatingValueLiquid(fuelFromSelection({ kind: 'pure', species: 'CH4' })))).toBe(true);
    // PRF blend LHV is the mole-weighted mixture (ideal solution)
    const prf = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
    const x = prf.X[SP.IC8H18];
    expect(lowerHeatingValueMolar(prf)).toBeCloseTo(
      x * lowerHeatingValueMolar(iso) + (1 - x) * lowerHeatingValueMolar(hep),
      3,
    );
  });

  it('element vector helper sanity (C, H of iso-octane)', () => {
    const b = elementMoles(iso.X);
    expect(b[EL.C]).toBe(8);
    expect(b[EL.H]).toBe(18);
  });
});
