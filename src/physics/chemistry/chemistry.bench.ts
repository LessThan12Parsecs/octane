/**
 * Micro-benchmarks for the chemistry hot paths. Run: npx vitest bench --run src/physics/chemistry
 */
import { bench, describe } from 'vitest';
import fx from '../../../test/fixtures/chemistry_zeldovich.json';
import { N_EQ, NS, SP } from '../core/species';
import { CFR_F1 } from '../engines/cfr';
import { douaudEyzatTau } from './ignition-delay';
import { prfDetailedChemistry } from './ignition-delay-llnl';
import { KnockOscillator } from './knock';
import { LivengoodWuIntegrator } from './livengood-wu';
import { ZeldovichKinetics, zeldovichNORate } from './zeldovich';

const Xeq = new Float64Array(NS);
for (let k = 0; k < N_EQ; k++) Xeq[k] = fx.states[40].Xeq[k];
const out = new Float64Array(3);
const kin = new ZeldovichKinetics();
let T = 2000;
let sink = 0;
const nextT = (): number => (T = T > 2790 ? 1800 : T + 7.3);
let Tu = 700;
const nextTu = (): number => (Tu = Tu > 1090 ? 660 : Tu + 3.1);

describe('chemistry hot paths', () => {
  bench('zeldovichNORate (new T every call)', () => {
    sink += zeldovichNORate(nextT(), 4e6, Xeq, 1e-3);
  });
  bench('ZeldovichKinetics.advanceRateControlled (exact step)', () => {
    sink += kin.advanceRateControlled(nextT(), 4e6, Xeq, 0.4 * Xeq[SP.NO], 2e-5);
  });
  bench('ZeldovichKinetics.fullRates (incl. reverse via thermo)', () => {
    kin.fullRates(nextT(), 4e6, Xeq, out);
    sink += out[0];
  });
  bench('douaudEyzatTau', () => {
    sink += douaudEyzatTau(nextTu(), 3e6, 90);
  });
  bench('LLNL table lnTauAt (5-D multilinear)', () => {
    sink += prfDetailedChemistry.lnTauAt(nextTu(), 3e6, 1.05, 90, 0.06);
  });
  const lw = new LivengoodWuIntegrator();
  lw.reset();
  bench('LivengoodWuIntegrator.advance (log-mean)', () => {
    if (lw.integral > 1e6) lw.reset();
    sink += lw.advance(1e-5, 1e-3 * (1 + 0.1 * Math.sin(Tu)));
    nextTu();
  });
  const ko = new KnockOscillator(CFR_F1.geometry.bore);
  ko.setEndGasRegion(0, 0.04, 0.06);
  let q = 1e6;
  bench(`KnockOscillator.step (${ko.nModes} modes)`, () => {
    q = q > 5e6 ? 1e6 : q * 1.01;
    ko.step(2e-5, 950, 1.3, 1.1e-4, q);
    sink += ko.sensorPressure();
  });
  bench('KnockOscillator.setEndGasRegion (once per knock event)', () => {
    sink += ko.setEndGasRegion(0, 0.04, 0.06);
  });
});

export { sink };
