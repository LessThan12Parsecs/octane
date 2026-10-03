/**
 * Micro-benchmarks for the thermo hot paths. Run: npx vitest bench --run src/physics/thermo
 */
import { bench, describe } from 'vitest';
import { NS } from '../core/species';
import { evalAllNondim, evalAllCpH } from './thermo';
import { mixHMass, mixCpMass, temperatureFromU, mixUMass } from './mixture';
import { dryAir, freshCharge, fuelFromSelection, completeCombustionProducts } from './fuels';
import { mixTransport, mixViscosity } from './transport';

const cpR = new Float64Array(NS);
const hRT = new Float64Array(NS);
const s0R = new Float64Array(NS);
const fuel = fuelFromSelection({ kind: 'PRF', octaneNumber: 90 });
const charge = freshCharge({ fuel, phi: 1, airX: dryAir(), egrFraction: 0.05 });
const products = completeCombustionProducts(charge);
const out = { viscosity: 0, conductivity: 0, thermalDiffusivity: 0, kinematicViscosity: 0, prandtl: 0 };
let T = 300;
let sink = 0;
const next = (): number => (T = T > 2900 ? 300 : T + 13.7);

describe('thermo hot paths', () => {
  bench('evalAllNondim (17 species)', () => {
    evalAllNondim(next(), cpR, hRT, s0R);
    sink += cpR[3];
  });
  bench('evalAllCpH (17 species)', () => {
    evalAllCpH(next(), cpR, hRT);
    sink += hRT[3];
  });
  bench('mixHMass (fresh charge, 8 species)', () => {
    sink += mixHMass(charge, next());
  });
  bench('mixHMass (products, 6 species)', () => {
    sink += mixHMass(products, next());
  });
  bench('mixCpMass (fresh charge)', () => {
    sink += mixCpMass(charge, next());
  });
  const u = mixUMass(products, 2345);
  bench('temperatureFromU (products, warm start ±50 K)', () => {
    sink += temperatureFromU(products, u, 2300 + (next() % 100));
  });
  bench('mixViscosity (fresh charge)', () => {
    sink += mixViscosity(charge, next());
  });
  bench('mixTransport (fresh charge, μ λ α ν Pr)', () => {
    mixTransport(charge, next(), 1e6, out);
    sink += out.prandtl;
  });
});

export { sink };
