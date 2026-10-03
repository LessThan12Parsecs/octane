/**
 * Micro-benchmarks for the equilibrium hot paths. Run: npx vitest bench --run src/physics/equilibrium
 */
import { bench, describe } from 'vitest';
import { NE, SP } from '../core/species';
import { freshCharge, fuelFromSelection, humidAir } from '../thermo/fuels';
import { elementMoles } from '../thermo/mixture';
import { EquilibriumSolver, newEqProperties } from './solver';

const X = freshCharge({ fuel: fuelFromSelection({ kind: 'PRF', octaneNumber: 90 }), phi: 1, airX: humidAir(298, 1e5, 0.5), egrFraction: 0.08 });
const b = elementMoles(X, new Float64Array(NE));
for (let e = 0; e < NE; e++) b[e] *= 0.012;
const s = new EquilibriumSolver();
const props = newEqProperties();
const r0 = s.solveTP(b, 2600, 50e5);
const U0 = r0.U;
const V0 = r0.V;
const H0 = r0.H;
let k = 0;
let sink = 0;
// slowly varying states (≈ 0.25° crank steps of an expansion stroke)
const T = (i: number): number => 2600 - 2 * (i % 400);
const p = (i: number): number => 50e5 * Math.pow(T(i) / 2600, 4.5);

describe('equilibrium hot paths', () => {
  bench('solveTP, warm, slowly varying (T, p)', () => {
    sink += s.solveTP(b, T(k), p(k++)).X[SP.NO];
  });
  // burned zone of ≈ 0.35 g: ΔU = −0.2 J per step ≈ −0.4 K; Tguess = previous T (as the cycle does)
  let Tg = 2600;
  bench('solveUV, warm, slowly varying (U, V)', () => {
    const i = k++ % 400;
    Tg = s.solveUV(b, U0 - 0.2 * i, V0 * (1 + 0.004 * i), Tg).T;
    sink += Tg;
  });
  bench('solveHP, warm, slowly varying H', () => {
    Tg = s.solveHP(b, H0 - 0.2 * (k++ % 400), 50e5, Tg).T;
    sink += Tg;
  });
  bench('solveTP + properties()', () => {
    s.solveTP(b, T(k), p(k++));
    sink += s.properties(props).cp;
  });
  bench('solveTP, cold', () => {
    s.reset();
    sink += s.solveTP(b, T(k), p(k++)).T;
  });
});

export { sink };
