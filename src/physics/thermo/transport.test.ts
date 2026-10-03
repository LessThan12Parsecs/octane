import { describe, expect, it } from 'vitest';
import { NS, SP, SPECIES } from '../core/species';
import trFx from '../../../test/fixtures/thermo_transport.json';
import { omega11LJ, omega22LJ } from './collision-integrals';
import {
  IAPWS_T_MAX,
  IAPWS_T_MIN,
  mixKinematicViscosity,
  mixThermalConductivity,
  mixThermalDiffusivity,
  mixTransport,
  mixViscosity,
  REDUCED_DIPOLE,
  speciesThermalConductivity,
  speciesThermalConductivityModifiedEucken,
  speciesThermalConductivityWarnatz,
  speciesViscosity,
  speciesViscosityChapmanEnskog,
  waterThermalConductivityIAPWS,
  waterViscosityIAPWS,
} from './transport';

const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);

describe('collision integrals', () => {
  it('Neufeld fits reproduce tabulated LJ values (HCB / Monchick–Mason δ*=0 column)', () => {
    // Monchick & Mason δ* = 0 column (= Hirschfelder–Curtiss–Bird LJ values)
    const tab22: [number, number][] = [
      [0.5, 2.2837], [0.6, 2.0838], [1, 1.5929], [2, 1.1757], [5, 0.92676], [10, 0.82435],
      [40, 0.67232], [50, 0.65099], [100, 0.5887],
    ];
    let worst = 0;
    for (const [t, v] of tab22) worst = Math.max(worst, rel(omega22LJ(t), v));
    console.info(`[transport] Neufeld Ω22* vs Monchick–Mason δ*=0 table: max ${(100 * worst).toFixed(2)} %`);
    expect(worst).toBeLessThan(5e-3);
    // Ω(1,1)* = Ω(2,2)*/A*, A*(T*=1) = 1.1063 (MM table): Ω11(1) = 1.5929/1.1063 = 1.4398
    expect(rel(omega11LJ(1), 1.5929 / 1.1063)).toBeLessThan(3e-3);
  });

  it('only H2O is polar in this species set; its δ* ≈ 1.22', () => {
    for (let k = 0; k < NS; k++) if (k !== SP.H2O) expect(REDUCED_DIPOLE[k]).toBe(0);
    expect(REDUCED_DIPOLE[SP.H2O]).toBeCloseTo(1.217, 2);
  });
});

describe('pure-species kinetic-theory model vs Cantera (same model, same parameters)', () => {
  it('Chapman–Enskog μ and Warnatz λ of every species (incl. the raw H2O model)', () => {
    let worstMu = 0;
    let worstLam = 0;
    let worstEucken = 0;
    let who = '';
    trFx.pureT.forEach((T, j) => {
      for (let k = 0; k < NS; k++) {
        const eMu = rel(speciesViscosityChapmanEnskog(k, T), trFx.pureViscosity[k][j]);
        const eLam = rel(speciesThermalConductivityWarnatz(k, T), trFx.pureConductivity[k][j]);
        const eE = rel(speciesThermalConductivityModifiedEucken(k, T), trFx.pureConductivity[k][j]);
        if (eLam > worstLam) who = `${SPECIES[k]}@${T}`;
        worstMu = Math.max(worstMu, eMu);
        worstLam = Math.max(worstLam, eLam);
        worstEucken = Math.max(worstEucken, eE);
        expect(eMu, `μ ${SPECIES[k]} ${T} K`).toBeLessThan(0.01);
        expect(eLam, `λ ${SPECIES[k]} ${T} K`).toBeLessThan(0.01);
        // every species except H2O uses exactly this model
        if (k !== SP.H2O) {
          expect(speciesViscosity(k, T)).toBe(speciesViscosityChapmanEnskog(k, T));
          expect(speciesThermalConductivity(k, T)).toBe(speciesThermalConductivityWarnatz(k, T));
        }
      }
    });
    console.info(
      `[transport] pure species 300–3000 K (Cantera fit error ≤ ${(100 * trFx.canteraFittingErrors['conductivity-max-rel-error']).toFixed(2)} %): max |Δμ/μ| = ${(100 * worstMu).toFixed(2)} %, ` +
        `max |Δλ/λ| = ${(100 * worstLam).toFixed(2)} % (${who}); modified Eucken would be ${(100 * worstEucken).toFixed(1)} %`,
    );
  });
});

describe('H2O: IAPWS dilute-gas correlations', () => {
  it('reproduce the IAPWS R15-11 check values and the oracle', () => {
    // IAPWS R15-11 (2011) Tables 4–5 (ρ = 0 / λ0 column), mW/(m K)
    expect(rel(waterThermalConductivityIAPWS(298.15), 18.4341883e-3)).toBeLessThan(1e-8);
    expect(rel(waterThermalConductivityIAPWS(647.35), 51.5764797e-3)).toBeLessThan(1e-8);
    expect(rel(waterThermalConductivityIAPWS(873.15), 79.1034659e-3)).toBeLessThan(1e-8);
    trFx.pureT.forEach((T, j) => {
      // exact inside the direct-use range; above it the scaled kinetic-theory value (the
      // oracle scales Cantera's fit of the same model, ≤ 0.5 % apart)
      const tol = T <= IAPWS_T_MAX ? 1e-13 : 0.01;
      expect(rel(speciesViscosity(SP.H2O, T), trFx.waterModelViscosity[j])).toBeLessThan(tol);
      expect(rel(speciesThermalConductivity(SP.H2O, T), trFx.waterModelConductivity[j])).toBeLessThan(tol);
    });
  });

  it('are continuous at the ends of their direct-use range and positive beyond', () => {
    for (const Tb of [IAPWS_T_MIN, IAPWS_T_MAX]) {
      for (const f of [speciesViscosity, speciesThermalConductivity]) {
        const a = f(SP.H2O, Tb * (1 - 1e-12));
        const b = f(SP.H2O, Tb * (1 + 1e-12));
        expect(rel(a, b)).toBeLessThan(1e-9);
      }
    }
    expect(speciesViscosity(SP.H2O, IAPWS_T_MAX)).toBe(waterViscosityIAPWS(IAPWS_T_MAX));
    let prevMu = 0;
    let prevLam = 0;
    for (let T = 100; T < 10000; T *= 1.02) {
      const mu = speciesViscosity(SP.H2O, T);
      const lam = speciesThermalConductivity(SP.H2O, T);
      expect(mu).toBeGreaterThan(prevMu);
      expect(lam).toBeGreaterThan(prevLam);
      prevMu = mu;
      prevLam = lam;
    }
  });
});

describe('pure-species transport vs reference data (physical accuracy)', () => {
  // NIST Chemistry WebBook fluid properties (REFPROP reference correlations: N2 — Lemmon &
  // Jacobsen 2004; CO2 — Laesecke & Muzny 2017 / Huber et al. 2016; H2O — IAPWS 2008/2011),
  // fetched 2026-09-29. Isobars: N2, CO2 at 0.1 MPa; H2O at 0.01 MPa (≈ dilute gas).
  // [T K, μ µPa s, λ W/(m K)]
  const N2: [number, number, number][] = [
    [300, 17.89, 0.025968], [600, 29.577, 0.04484], [900, 38.78, 0.06052],
    [1200, 46.767, 0.074679], [1500, 54.066, 0.088015], [1800, 60.95, 0.10088],
  ];
  const CO2: [number, number, number][] = [
    [300, 15.003, 0.016774], [500, 23.923, 0.032876], [700, 31.533, 0.048833],
    [900, 38.148, 0.063752], [1100, 44.073, 0.077531],
  ];
  const H2O: [number, number, number][] = [
    [400, 13.347, 0.02647], [600, 21.432, 0.046291], [800, 29.655, 0.06984],
    [1000, 37.611, 0.095812], [1200, 45.188, 0.12336],
  ];
  const check = (k: number, rows: [number, number, number][], tolMu: number, tolLam: number): string => {
    let wMu = 0;
    let wLam = 0;
    for (const [T, mu, lam] of rows) {
      const eMu = rel(speciesViscosity(k, T), mu * 1e-6);
      const eLam = rel(speciesThermalConductivity(k, T), lam);
      wMu = Math.max(wMu, eMu);
      wLam = Math.max(wLam, eLam);
      expect(eMu, `${SPECIES[k]} μ @${T}`).toBeLessThan(tolMu);
      expect(eLam, `${SPECIES[k]} λ @${T}`).toBeLessThan(tolLam);
    }
    return `${SPECIES[k]} μ ${(100 * wMu).toFixed(2)} %, λ ${(100 * wLam).toFixed(2)} %`;
  };

  it('N2, CO2 (kinetic theory) and H2O (IAPWS) against NIST reference correlations', () => {
    const lines = [
      // Kinetic-theory model limits (documented in transport.ts): the Eucken-type
      // internal-energy term over-predicts N2 λ by up to ~9 % at 1500–1800 K.
      check(SP.N2, N2, 0.015, 0.1),
      check(SP.CO2, CO2, 0.01, 0.05),
      // H2O: IAPWS dilute gas vs the full IAPWS formulation at 10 kPa (density effect ≤ 6e-4).
      // The Warnatz model it replaces is 20–40 % high in λ here.
      check(SP.H2O, H2O, 1e-3, 2e-3),
    ];
    console.info(`[transport] vs NIST reference data: ${lines.join('; ')}`);
    for (const [T, , lam] of H2O) expect(speciesThermalConductivityWarnatz(SP.H2O, T) / lam).toBeGreaterThan(1.15);
  });
});

describe('mixture transport (300–3000 K)', () => {
  it('matches the oracle: Wilke / Mathur–Saxena over Cantera pure-species values (H2O: IAPWS)', () => {
    const lines: string[] = [];
    let worst = 0;
    let worstCantera = 0;
    for (const c of trFx.mixtures) {
      const X = Float64Array.from(c.X);
      const hasWater = X[SP.H2O] > 0;
      let wMu = 0;
      let wLam = 0;
      let wAl = 0;
      let wC = 0;
      for (const st of c.states) {
        const mu = mixViscosity(X, st.T);
        const lam = mixThermalConductivity(X, st.T);
        const al = mixThermalDiffusivity(X, st.T, c.p);
        const eMu = rel(mu, st.modelViscosity);
        const eLam = rel(lam, st.modelConductivity);
        const eAl = rel(al, st.modelThermalDiffusivity);
        wMu = Math.max(wMu, eMu);
        wLam = Math.max(wLam, eLam);
        wAl = Math.max(wAl, eAl);
        expect(eMu, `${c.name} μ @${st.T}`).toBeLessThan(0.01);
        expect(eLam, `${c.name} λ @${st.T}`).toBeLessThan(0.01);
        expect(eAl, `${c.name} α @${st.T}`).toBeLessThan(0.01);
        // H2O-free mixtures: the oracle IS Cantera's mixture-averaged transport
        if (!hasWater) {
          const e = Math.max(rel(mu, st.viscosity), rel(lam, st.conductivity), rel(al, st.thermalDiffusivity));
          wC = Math.max(wC, e);
          expect(e, `${c.name} vs Cantera @${st.T}`).toBeLessThan(0.01);
        }
      }
      worst = Math.max(worst, wMu, wLam, wAl);
      worstCantera = Math.max(worstCantera, wC);
      lines.push(`  ${c.name}: μ ${(100 * wMu).toFixed(2)} %, λ ${(100 * wLam).toFixed(2)} %, α ${(100 * wAl).toFixed(2)} %`);
    }
    console.info(
      `[transport] max error vs oracle per mixture (H2O-free mixtures vs Cantera directly: ` +
        `${(100 * worstCantera).toFixed(2)} %):\n${lines.join('\n')}`,
    );
    expect(worst).toBeLessThan(0.01);
  });

  it('mixTransport agrees with the individual functions', () => {
    const c = trFx.mixtures[1];
    const X = Float64Array.from(c.X);
    const T = 750;
    const p = 2e6;
    const r = mixTransport(X, T, p);
    expect(r.viscosity).toBeCloseTo(mixViscosity(X, T), 15);
    expect(r.conductivity).toBeCloseTo(mixThermalConductivity(X, T), 15);
    expect(rel(r.thermalDiffusivity, mixThermalDiffusivity(X, T, p))).toBeLessThan(1e-14);
    expect(rel(r.kinematicViscosity, mixKinematicViscosity(X, T, p))).toBeLessThan(1e-14);
    expect(r.prandtl).toBeGreaterThan(0.6);
    expect(r.prandtl).toBeLessThan(0.8);
    // with water (burned gas)
    const Xb = Float64Array.from(trFx.mixtures[4].X);
    const rb = mixTransport(Xb, 2400, 4e6);
    expect(rel(rb.viscosity, mixViscosity(Xb, 2400))).toBeLessThan(1e-14);
    expect(rel(rb.conductivity, mixThermalConductivity(Xb, 2400))).toBeLessThan(1e-14);
  });

  it('pure-species limit of the mixture rules', () => {
    for (const k of [SP.N2, SP.H2O]) {
      const X = new Float64Array(NS);
      X[k] = 1;
      for (const T of [220, 900, 2800]) {
        expect(rel(mixViscosity(X, T), speciesViscosity(k, T))).toBeLessThan(1e-14);
        expect(rel(mixThermalConductivity(X, T), speciesThermalConductivity(k, T))).toBeLessThan(1e-14);
      }
    }
  });
});
