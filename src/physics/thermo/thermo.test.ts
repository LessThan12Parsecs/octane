import { describe, expect, it } from 'vitest';
import { ELEMENT_COUNTS, ELEMENTS, NS, SP, SPECIES } from '../core/species';
import speciesFx from '../../../test/fixtures/thermo_species.json';
import { mixSMolar } from './mixture';
import { CANTERA_NASA_GAS_P_REF, SPECIES_NASA7 } from './species-data';
import {
  evalAllCpH,
  evalAllG0RT,
  evalAllNondim,
  MOLAR_MASS,
  P_REF_THERMO,
  speciesCp,
  speciesCpR,
  speciesG0,
  speciesH,
  speciesHRT,
  speciesS0,
  speciesS0R,
  speciesTmax,
  speciesTmid,
  speciesTmin,
} from './thermo';
import { R_UNIVERSAL } from '../core/constants';

const relErr = (a: number, b: number, scale = Math.abs(b)): number => Math.abs(a - b) / Math.max(scale, 1e-300);

describe('species data', () => {
  it('matches SPECIES order and core ELEMENT_COUNTS', () => {
    expect(SPECIES_NASA7.map((r) => r.name)).toEqual([...SPECIES]);
    const sym: Record<string, string> = { C: 'C', H: 'H', O: 'O', N: 'N', AR: 'Ar' };
    for (let k = 0; k < NS; k++) {
      for (let e = 0; e < ELEMENTS.length; e++) {
        const n = SPECIES_NASA7[k].composition[sym[ELEMENTS[e]]] ?? 0;
        expect(n, `${SPECIES[k]} ${ELEMENTS[e]}`).toBe(ELEMENT_COUNTS[k][e]);
      }
    }
  });

  it('MOLAR_MASS equals Cantera molecular weights to 1e-12', () => {
    let worst = 0;
    for (let k = 0; k < NS; k++) {
      worst = Math.max(worst, relErr(MOLAR_MASS[k], speciesFx.molarMass[k]));
      expect(relErr(MOLAR_MASS[k], SPECIES_NASA7[k].canteraMolarMass)).toBeLessThan(1e-12);
    }
    expect(worst).toBeLessThan(1e-12);
  });

  it('P_REF_THERMO is the 1-bar standard state of the NASA TM-4513 data (and of the oracle phase)', () => {
    // McBride, Gordon & Reno (1993) TM-4513: ideal-gas standard state 1e5 Pa. Cantera's
    // nasa_gas.yaml labels the same fits 1 atm (ck2yaml default); the oracle rebuilds them at 1 bar.
    expect(P_REF_THERMO).toBe(1e5);
    expect(P_REF_THERMO).toBe(speciesFx.pRef);
    expect(CANTERA_NASA_GAS_P_REF).toBe(101325);
  });

  it('standard entropies at 298.15 K equal the CODATA 1-bar key values', () => {
    // CODATA Key Values for Thermodynamics (Cox, Wagman & Medvedev 1989), S°(298.15 K),
    // J/(mol K), standard-state pressure 1 bar (codata.info key1 table, fetched 2026-09-29).
    // A 1-atm basis would be R ln(1.01325) = 0.1094 J/(mol K) lower — this is what pins
    // P_REF_THERMO: s(T, p) = s° − R ln(p/P_REF_THERMO) must use the data's own p°.
    const codata: [number, number][] = [
      [SP.O2, 205.152], [SP.N2, 191.609], [SP.H2, 130.68], [SP.H2O, 188.835], [SP.CO, 197.66],
      [SP.CO2, 213.785], [SP.AR, 154.846], [SP.H, 114.717], [SP.O, 161.059], [SP.N, 153.301],
    ];
    let worst = 0;
    for (const [k, s] of codata) {
      const d = speciesS0(k, 298.15) - s;
      worst = Math.max(worst, Math.abs(d));
      expect(Math.abs(d), SPECIES[k]).toBeLessThan(0.01);
    }
    console.info(`[thermo] s°(298.15 K) vs CODATA (1 bar): max |Δ| = ${worst.toFixed(4)} J/(mol K)`);
    // Entropy of pure O2 at 1 atm, 298.15 K (ideal gas) = S°(1 bar) − R ln(101325/1e5)
    const X = new Float64Array(NS);
    X[SP.O2] = 1;
    expect(Math.abs(mixSMolar(X, 298.15, 101325) - (205.152 - R_UNIVERSAL * Math.log(1.01325)))).toBeLessThan(0.01);
  });
});

describe('species thermo vs Cantera (200–6000 K)', () => {
  it('cp/R, h/RT, s°/R to 1e-9 relative', () => {
    let worst = 0;
    const cpR = new Float64Array(NS);
    const hRT = new Float64Array(NS);
    const s0R = new Float64Array(NS);
    speciesFx.T.forEach((T, j) => {
      evalAllNondim(T, cpR, hRT, s0R);
      for (let k = 0; k < NS; k++) {
        const eCp = relErr(cpR[k], speciesFx.cpR[k][j]);
        // h/RT passes through ~0 near 298 K for elements: scale by max(|h/RT|, cp/R)
        const eH = relErr(hRT[k], speciesFx.hRT[k][j], Math.max(Math.abs(speciesFx.hRT[k][j]), speciesFx.cpR[k][j]));
        const eS = relErr(s0R[k], speciesFx.s0R[k][j]);
        worst = Math.max(worst, eCp, eH, eS);
        expect(eCp, `${SPECIES[k]} cp @${T}`).toBeLessThan(1e-9);
        expect(eH, `${SPECIES[k]} h @${T}`).toBeLessThan(1e-9);
        expect(eS, `${SPECIES[k]} s @${T}`).toBeLessThan(1e-9);
        // single-species kernels agree exactly with the batched hot path
        expect(relErr(speciesCpR(k, T), cpR[k])).toBeLessThan(1e-13);
        expect(relErr(speciesHRT(k, T), hRT[k], Math.abs(hRT[k]) + cpR[k])).toBeLessThan(1e-13);
        expect(relErr(speciesS0R(k, T), s0R[k])).toBeLessThan(1e-13);
      }
    });
    console.info(`[thermo] species max rel. error vs Cantera: ${worst.toExponential(2)}`);
  });

  it('dimensional wrappers and g° are consistent', () => {
    const g = new Float64Array(NS);
    for (const T of [250, 800, 1000, 1700, 4200]) {
      evalAllG0RT(T, g);
      for (let k = 0; k < NS; k++) {
        expect(relErr(speciesCp(k, T), R_UNIVERSAL * speciesCpR(k, T))).toBeLessThan(1e-15);
        const gk = speciesH(k, T) - T * speciesS0(k, T);
        expect(relErr(speciesG0(k, T), gk, Math.abs(gk) + R_UNIVERSAL * T)).toBeLessThan(1e-12);
        expect(relErr(g[k] * R_UNIVERSAL * T, gk, Math.abs(gk) + R_UNIVERSAL * T)).toBeLessThan(1e-12);
      }
    }
    const cpR = new Float64Array(NS);
    const hRT = new Float64Array(NS);
    const cpR2 = new Float64Array(NS);
    const hRT2 = new Float64Array(NS);
    const s = new Float64Array(NS);
    for (const T of [150, 700, 3000, 7000]) {
      evalAllCpH(T, cpR, hRT);
      evalAllNondim(T, cpR2, hRT2, s);
      for (let k = 0; k < NS; k++) {
        expect(cpR[k]).toBe(cpR2[k]);
        expect(hRT[k]).toBe(hRT2[k]);
      }
    }
  });
});

describe('extrapolation outside the fit range', () => {
  it('is continuous at Tmin/Tmax and thermodynamically consistent (dh/dT = cp, ds/dT = cp/T)', () => {
    for (let k = 0; k < NS; k++) {
      for (const Tb of [speciesTmin(k), speciesTmax(k)]) {
        for (const f of [speciesCpR, speciesHRT, speciesS0R]) {
          const a = f(k, Tb * (1 - 1e-12));
          const b = f(k, Tb * (1 + 1e-12));
          expect(Math.abs(a - b)).toBeLessThan(1e-8 * Math.max(1, Math.abs(a)));
        }
      }
      for (const T of [50, 120, 7000, 12000, 30000]) {
        const dT = 1e-3 * T;
        const dh = (speciesH(k, T + dT) - speciesH(k, T - dT)) / (2 * dT);
        const ds = (speciesS0(k, T + dT) - speciesS0(k, T - dT)) / (2 * dT);
        const cp = speciesCp(k, T);
        expect(relErr(dh, cp)).toBeLessThan(1e-7);
        expect(relErr(ds, cp / T)).toBeLessThan(1e-5);
        // constant cp beyond the range: no polynomial blow-up
        const cb = T < 200 ? speciesCp(k, speciesTmin(k)) : speciesCp(k, speciesTmax(k));
        expect(cp).toBe(cb);
      }
    }
  });
});

describe('Tmid discontinuity (documented property of NASA-7 fits)', () => {
  it('is small for every species', () => {
    let worstCp = 0;
    let worstH = 0;
    let worstS = 0;
    for (let k = 0; k < NS; k++) {
      const Tm = speciesTmid(k);
      if (Tm >= speciesTmax(k)) continue; // single-range species
      const lo = Tm;
      const hi = Tm * (1 + 1e-14);
      const dcp = Math.abs(speciesCpR(k, hi) - speciesCpR(k, lo)) / speciesCpR(k, lo);
      const dh = Math.abs(speciesH(k, hi) - speciesH(k, lo)); // J/mol
      const ds = Math.abs(speciesS0(k, hi) - speciesS0(k, lo)); // J/mol/K
      worstCp = Math.max(worstCp, dcp);
      worstH = Math.max(worstH, dh);
      worstS = Math.max(worstS, ds);
    }
    console.info(
      `[thermo] Tmid jumps: max |Δcp/cp| = ${worstCp.toExponential(2)}, max |Δh| = ${worstH.toFixed(3)} J/mol, ` +
        `max |Δs| = ${worstS.toExponential(2)} J/mol/K`,
    );
    expect(worstCp).toBeLessThan(1e-6);
    expect(worstH).toBeLessThan(1e-2);
    expect(worstS).toBeLessThan(1e-4);
  });
});
