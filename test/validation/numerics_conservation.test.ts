/**
 * Validation round 1 (numerics): conservation of the whole network (cylinder + intake/exhaust
 * plenums + boundary ledgers) over 12 cycles from a cold start, checked at 7 points per cycle
 * (open and closed phases), for operating points the author's tests do not cover: burn-out after
 * knock (burnRateMultiplier 2.6), throttled + EGR + lean, free speed with load, Wiebe,
 * instantaneous combustion, motored, and operating-point changes mid-run (CR, fuel, φ, rpm, EGR).
 * Plus the closed-cycle heat and work ledgers against independent trapezoidal integrals.
 * Runtime ≈ 20–40 s.
 */
import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import { NE, NS } from '../../src/physics/core/species';
import { CFR_F1, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { completeCombustionProducts, freshCharge, fuelFromSelection, humidAir } from '../../src/physics/thermo/fuels';
import { elementMoles, mixMolarMass } from '../../src/physics/thermo/mixture';
import { MOLAR_MASS } from '../../src/physics/thermo/thermo';
import { CycleModel, I_HO, I_HV, I_LO, I_LV, I_Q, I_W } from '../../src/physics/cycle/cycle-model';
import { runClosedCycle, type CycleModelOptions } from '../../src/physics/cycle/index';

function ledgers(m: CycleModel) {
  const y = m.y;
  const LV = new Float64Array(NS);
  const LO = new Float64Array(NS);
  let mv = 0;
  let mo = 0;
  for (let k = 0; k < NS; k++) {
    LV[k] = y[I_LV + k];
    LO[k] = y[I_LO + k];
    mv += LV[k] * MOLAR_MASS[k];
    mo += LO[k] * MOLAR_MASS[k];
  }
  const b = m.elementInventory();
  const bV = elementMoles(LV);
  const bO = elementMoles(LO);
  const elements = new Float64Array(NE);
  for (let e = 0; e < NE; e++) elements[e] = b[e] - bV[e] + bO[e];
  return {
    mass: m.systemMass() - mv + mo,
    energy: m.systemEnergy() - y[I_HV] + y[I_HO] + y[I_Q] + y[I_W] - m.sparkEnergy,
    elements,
  };
}

const CASES: [string, Partial<OperatingPoint>, Partial<CycleModelOptions>, Partial<OperatingPoint>[]?][] = [
  ['RON, burnRateMultiplier 2.6 (knock + burn-out merge)', {}, { burnRateMultiplier: 2.6 }],
  ['throttle 0.1, φ 0.9, EGR 0.2', { throttle: 0.1, equivalenceRatio: 0.9, egrFraction: 0.2 }, {}],
  ['1500 rpm, CR 8, PRF 100', { rpm: 1500, compressionRatio: 8, fuel: { kind: 'PRF', octaneNumber: 100 } }, {}],
  ['free speed, load 5 N m', { speedMode: 'free', loadTorque: 5 }, {}],
  ['Wiebe', {}, { combustionModel: 'wiebe' }],
  ['instantaneous at TDC', {}, { combustionModel: 'instantaneous-at-tdc' }],
  ['motored', {}, { combustionModel: 'none' }],
  [
    'operating-point changes mid-run',
    {},
    {},
    [{ compressionRatio: 8 }, { fuel: { kind: 'pure', species: 'C2H5OH' }, equivalenceRatio: 0.9 }, { rpm: 900 }, { fuel: { kind: 'pure', species: 'CH4' }, egrFraction: 0.1 }],
  ],
];

describe('whole-network conservation over 12 cycles (cold start, plenums included)', () => {
  for (const [label, opP, opts, changes] of CASES) {
    it(label, () => {
      const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, ...opP }, { warmupCycles: 0, ...opts });
      const l0 = ledgers(m);
      let wM = 0;
      let wE = 0;
      let wEl = 0;
      for (let c = 0; c < 12; c++) {
        if (changes && c % 3 === 1) m.setOperatingPoint(changes[Math.floor(c / 3)] ?? {});
        for (const th of [-200, -100, 0, 60, 150, 300, 360]) {
          m.stepUntil(Infinity, th);
          const l = ledgers(m);
          wM = Math.max(wM, Math.abs(l.mass / l0.mass - 1));
          wE = Math.max(wE, Math.abs(l.energy - l0.energy));
          for (let e = 0; e < NE; e++) if (l0.elements[e] > 0) wEl = Math.max(wEl, Math.abs(l.elements[e] / l0.elements[e] - 1));
        }
      }
      console.log(`[validation] conservation, ${label}: mass ${wM.toExponential(1)}, elements ${wEl.toExponential(1)}, energy ${wE.toExponential(1)} J`);
      expect(wM).toBeLessThan(1e-12);
      expect(wEl).toBeLessThan(1e-12);
      expect(wE).toBeLessThan(1e-8); // J (fuel energy per cycle ≈ 2000 J)
    });
  }
});

describe('closed-cycle ledgers vs independent quadrature (heat transfer + knock, CR 9)', () => {
  it('∫Q̇ dt and ∫p dV agree with trapezoids of the traced rates', () => {
    const op = CFR_RON_CONDITIONS;
    const fuel = fuelFromSelection(op.fuel);
    const air = humidAir(op.ambientTemperature, op.ambientPressure, op.relativeHumidity);
    const fresh = freshCharge({ fuel, phi: op.equivalenceRatio, airX: air });
    const prod = completeCombustionProducts(fresh);
    const yRes = 0.06;
    const X = new Float64Array(NS);
    for (let k = 0; k < NS; k++) X[k] = ((1 - yRes) / mixMolarMass(fresh)) * fresh[k] + (yRes / mixMolarMass(prod)) * prod[k];
    const r = runClosedCycle(
      CFR_F1,
      { ...op, compressionRatio: 9 },
      { startDeg: -152, endDeg: 141, T: 370, p: 1e5, X, residualMassFraction: yRes, uPrime: 5 },
      { heatTransfer: true, knock: true },
    );
    const tr = r.trace;
    let Qt = 0;
    let Wt = 0;
    for (let i = 1; i < tr.t.length; i++) {
      const dt = tr.t[i] - tr.t[i - 1];
      if (dt > 0) Qt += 0.5 * (tr.heatLossRate[i] + tr.heatLossRate[i - 1]) * dt;
      Wt += 0.5 * (tr.pressure[i] + tr.pressure[i - 1]) * (tr.volume[i] - tr.volume[i - 1]);
    }
    const eQ = Qt / r.heatLoss - 1;
    const eW = Wt / r.work - 1;
    console.log(`[validation] closed-cycle ledgers: Q ${r.heatLoss.toFixed(3)} J (trapezoid ${eQ.toExponential(1)}), W ${r.work.toFixed(3)} J (trapezoid ${eW.toExponential(1)})`);
    expect(Math.abs(eQ)).toBeLessThan(2e-4);
    expect(Math.abs(eW)).toBeLessThan(2e-4);
    expect(Math.abs(r.U1 + r.work + r.heatLoss - r.sparkEnergy - r.U0) / (r.fuelMass * 44e6)).toBeLessThan(1e-12);
  });
});
