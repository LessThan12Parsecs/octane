import { describe, expect, it } from 'vitest';
import { P_ATM, T_REF } from '../core/constants';
import { NS, SP } from '../core/species';
import { dryAir, freshCharge, fuelFromSelection } from '../thermo/fuels';
import stFx from '../../../test/fixtures/equilibrium_stress.json';
import { adiabaticFlameTemperature } from './flame';

/** Textbook "air": 1 O2 + 3.76 N2 (the convention of the textbook tables below). */
function textbookAir(): Float64Array {
  const X = new Float64Array(NS);
  X[SP.O2] = 1 / 4.76;
  X[SP.N2] = 3.76 / 4.76;
  return X;
}

const stoich = (species: 'CH4' | 'IC8H18' | 'C3H8', airX: Float64Array): Float64Array =>
  freshCharge({ fuel: fuelFromSelection({ kind: 'pure', species }), phi: 1, airX });

describe('adiabaticFlameTemperature: textbook values (298.15 K, 1 atm, stoichiometric)', () => {
  // UNVERIFIED: 2226 K (CH4) is the equilibrium adiabatic flame temperature tabulated by Turns,
  // "An Introduction to Combustion" (2nd ed.), Appendix B, Table B.1, for fuel–air (3.76 N2) at
  // 298 K and 1 atm; the table itself could not be fetched (value from memory).
  // Corroboration (fetched 2026-09-29): Cantera + GRI-Mech 3.0 (53 species) gives 2224.25 K for
  // stoichiometric CH4–air (arXiv:2503.11826, abstract); adding 15 minor NASA species (HO2, NO2,
  // N2O, …) to our 12 moves T_ad by < 0.02 K (equilibrium_stress.json, stress.test.ts).
  it('CH4–air ≈ 2226 K', () => {
    const f = adiabaticFlameTemperature(stoich('CH4', textbookAir()), T_REF, P_ATM, 'HP');
    console.log(`CH4-air (3.76 N2): T_ad = ${f.T.toFixed(2)} K; with CIPM dry air: ` +
      `${adiabaticFlameTemperature(stoich('CH4', dryAir()), T_REF, P_ATM, 'HP').T.toFixed(2)} K`);
    expect(f.converged).toBe(true);
    expect(Math.abs(f.T - 2226)).toBeLessThan(1);
  });

  // Review fix: this test used to compare iso-octane with Turns' 2275 K ± 5 K, but that Table B.1
  // row is n-octane (h_f = −208 447 kJ/kmol). NIST WebBook ΔfH°gas (fetched 2026-09-29, Prosen &
  // Rossini 1945): n-octane −208.4 ± 0.67 kJ/mol, 2,2,4-trimethylpentane −224.1 ± 1.3 kJ/mol, so
  // iso-octane burns ≈ 3.7 K cooler: n-octane via its NASA fit gives 2275.07 K (= Turns), iso-octane
  // 2271.41 K — both checked against the Cantera oracle in stress.test.ts. Here: iso-octane with our
  // own fuel thermo equals that oracle value.
  it('iso-octane–air = 2271.41 K (Turns 2275 K is n-octane)', () => {
    const f = adiabaticFlameTemperature(stoich('IC8H18', textbookAir()), T_REF, P_ATM, 'HP');
    console.log(`iC8H18-air (3.76 N2): T_ad = ${f.T.toFixed(2)} K; with CIPM dry air: ` +
      `${adiabaticFlameTemperature(stoich('IC8H18', dryAir()), T_REF, P_ATM, 'HP').T.toFixed(2)} K`);
    expect(f.converged).toBe(true);
    const ref = stFx.textbook.find((c) => c.label.startsWith('iso-C8H18'))!.T;
    expect(Math.abs(f.T - ref)).toBeLessThan(0.1);
  });

  it('constant-volume flame is hotter than constant-pressure, and p rises accordingly', () => {
    const X = stoich('C3H8', dryAir());
    const hp = adiabaticFlameTemperature(X, T_REF, P_ATM, 'HP');
    const uv = adiabaticFlameTemperature(X, T_REF, P_ATM, 'UV');
    expect(uv.converged).toBe(true);
    expect(uv.T).toBeGreaterThan(hp.T + 250);
    // p_final/p0 = (n_b T_b)/(n_r T0)
    expect(uv.p / P_ATM).toBeCloseTo((uv.molesPerMoleReactant * uv.T) / T_REF, 9);
    expect(hp.p).toBe(P_ATM);
  });

  it('peak flame temperature sits slightly rich of stoichiometric', () => {
    const fuel = fuelFromSelection({ kind: 'pure', species: 'IC8H18' });
    let best = 0;
    let bestPhi = 0;
    for (let phi = 0.9; phi <= 1.2 + 1e-9; phi += 0.01) {
      const f = adiabaticFlameTemperature(freshCharge({ fuel, phi, airX: dryAir() }), T_REF, P_ATM, 'HP');
      if (f.T > best) {
        best = f.T;
        bestPhi = phi;
      }
    }
    expect(bestPhi).toBeGreaterThan(1.0);
    expect(bestPhi).toBeLessThan(1.1);
  });
});
