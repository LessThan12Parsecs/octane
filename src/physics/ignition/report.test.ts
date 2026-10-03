/**
 * Sensitivity tables of the ignition model (printed only with IGNITION_REPORT=1; the
 * assertions are coarse sanity checks). Conditions: CFR-like spark states of the Cantera
 * fixture (CR 7, 13° BTDC), CFR-like coil, 0.508 mm gap, injected flame inputs as in
 * index.test.ts.
 */
import { describe, expect, it } from 'vitest';
import type { IgnitionSystemSpec, SparkPlugSpec } from '../core/engine-spec';
import gasFixture from '../../../test/fixtures/ignition_kernel_gas.json';
import { breakdownVoltage, IgnitionSystem, type IgnitionGasState, type IgnitionSystemOptions } from './index';

const REPORT = !!process.env.IGNITION_REPORT;
const log = (s: string): void => {
  if (REPORT) console.log(s);
};

const COIL: IgnitionSystemSpec = {
  type: 'inductive',
  supplyVoltage: 13.5,
  primaryInductance: 4e-3,
  primaryResistance: 1.5,
  secondaryInductance: 40,
  secondaryResistance: 8e3,
  secondaryCapacitance: 60e-12,
  couplingCoefficient: 0.98,
  primaryCurrentLimit: 7,
  dwellTime: 3e-3,
};
const PLUG: SparkPlugSpec = {
  gapCenter: [0, -0.002, 0.04],
  gap: 0.508e-3,
  centerElectrodeDiameter: 2.5e-3,
  groundElectrodeWidth: 2.5e-3,
  axis: [0, 0, -1],
  threadDiameter: 18e-3,
};

function gasAt(cr: number, phi: number, over: Partial<IgnitionGasState> = {}): IgnitionGasState {
  const c = gasFixture.cases.find((k) => k.cr === cr && k.phi === phi)!;
  const SL = phi === 1 ? 0.6 : phi === 0.7 ? 0.28 : 0.15;
  const Ma = phi === 1 ? 3 : phi === 0.7 ? 3 : 6;
  const lF = c.nuU / SL;
  return {
    p: c.p,
    Tu: c.Tu,
    rhoU: c.rhoU,
    X: Float64Array.from(c.X),
    SL,
    marksteinLength: Ma * lF,
    flameThickness: lF,
    uPrime: 1.1,
    integralScale: 2e-3,
    expansionRatio: c.sigma,
    kinematicViscosity: c.nuU,
    ...over,
  };
}

function run(gas: IgnitionGasState, dwell: number, opts: IgnitionSystemOptions = {}, tRun = 8e-3) {
  const ign = new IgnitionSystem(COIL, PLUG, opts);
  const dt = 20e-6;
  for (let t = 0; t < dwell - 1e-12; t += dt) ign.stepTimed(dt, true, gas);
  for (let t = 0; t < tRun - 1e-12; t += dt) ign.stepTimed(dt, false, gas);
  return ign;
}

const f = (x: number, d = 2): string => (Number.isNaN(x) ? '   –  ' : x.toFixed(d).padStart(6));

describe('ignition sensitivity report', () => {
  it('spark energy (dwell) → discharge and kernel outcome, φ = 1 / 0.7 / 0.5', () => {
    log('\nφ    dwell  E_coil  E_gap  E_gas  E_bd   dur    V_bd   kernel    t_handoff  r_max');
    log('     ms     mJ      mJ     mJ     mJ     ms     kV               ms         mm');
    for (const phi of [1, 0.7, 0.5]) {
      for (const dwell of [0.4e-3, 0.7e-3, 1e-3, 1.5e-3, 2e-3, 3e-3]) {
        const ign = run(gasAt(7, phi), dwell);
        const s = ign.state;
        log(
          `${phi.toFixed(1)}  ${f(dwell * 1e3, 1)} ${f(s.energyAtSwitchOff * 1e3, 1)} ${f(s.energyDelivered * 1e3, 1)} ${f(
            s.energyToGas * 1e3,
            1,
          )} ${f(s.energyBreakdown * 1e3, 2)} ${f(s.sparkDuration * 1e3, 2)} ${f(s.lastBreakdownVoltage / 1e3, 2)}   ${s.kernel.stage.padEnd(8)} ${f(
            s.kernel.handoffTime * 1e3,
            3,
          )}    ${f(ign.kernel.maxRadius * 1e3, 2)}`,
        );
      }
    }
    expect(true).toBe(true);
  });

  it('turbulence intensity u′ and integral scale l_I → kernel hand-off time (φ = 1, 3 ms dwell)', () => {
    log('\nu′ m/s  l_I mm  r_handoff mm  t_handoff ms  I0     K_AGB');
    let prev = Infinity;
    for (const lI of [1e-3, 2e-3, 4e-3]) {
      for (const up of [0, 0.5, 1.1, 2, 4, 8]) {
        const ign = run(gasAt(7, 1, { uPrime: up, integralScale: lI }), 3e-3);
        const k = ign.state.kernel;
        log(`${f(up, 1)}  ${f(lI * 1e3, 1)}  ${f(k.handoffRadius * 1e3, 2)}        ${f(k.handoffTime * 1e3, 3)}       ${f(k.stretchFactor, 3)} ${f(k.karlovitz, 3)}  ${k.stage}`);
        if (lI === 4e-3 && up <= 2) {
          expect(k.handoffTime).toBeLessThan(prev + 1e-9);
          prev = k.handoffTime;
        }
      }
    }
  });

  it('Markstein number Ma = L/l_F → hand-off / misfire (φ = 1 and 0.7, 3 ms dwell)', () => {
    log('\nφ    Ma    t_handoff ms  stage');
    for (const phi of [1, 0.7]) {
      for (const Ma of [-1, 0, 1, 3, 5, 8, 12]) {
        const g = gasAt(7, phi);
        const ign = run({ ...g, marksteinLength: Ma * g.flameThickness }, 3e-3);
        log(`${phi.toFixed(1)}  ${f(Ma, 0)}  ${f(ign.state.kernel.handoffTime * 1e3, 3)}       ${ign.state.kernel.stage} ${ign.state.kernel.quenchReason}`);
      }
    }
  });

  it('compression ratio / gap → breakdown voltage, glow voltage, spark duration and gas share', () => {
    log('\nCR  gap mm  V_bd kV  V_glow(50 mA) V  dur ms  E_gap mJ  η_gas');
    for (const cr of [7, 10]) {
      for (const gap of [0.4e-3, 0.508e-3, 0.7e-3, 0.9e-3, 1.1e-3]) {
        const g = gasAt(cr, 1);
        const ign = new IgnitionSystem(COIL, { ...PLUG, gap });
        const dt = 20e-6;
        for (let t = 0; t < 3e-3 - 1e-12; t += dt) ign.stepTimed(dt, true, g);
        for (let t = 0; t < 6e-3 - 1e-12; t += dt) ign.stepTimed(dt, false, g);
        const s = ign.state;
        ign.gap.mode = 'glow';
        const vg = ign.gap.gapVoltage(0.05);
        log(`${cr}   ${f(gap * 1e3, 3)}  ${f(breakdownVoltage(gap, g.p, g.Tu) / 1e3, 2)}   ${f(vg, 0)}             ${f(s.sparkDuration * 1e3, 2)}  ${f(s.energyDelivered * 1e3, 1)}    ${f(s.energyToGas / s.energyDelivered, 3)}`);
      }
    }
    expect(true).toBe(true);
  });

  it('CPU cost of one full ignition event (3 ms dwell + 8 ms, 20 µs engine steps)', () => {
    const g = gasAt(7, 1);
    run(g, 3e-3); // warm-up
    const t0 = performance.now();
    const n = 20;
    for (let i = 0; i < n; i++) run(g, 3e-3);
    const ms = (performance.now() - t0) / n;
    log(`\nCPU per ignition event: ${ms.toFixed(2)} ms`);
    expect(ms).toBeLessThan(200);
  });
});
