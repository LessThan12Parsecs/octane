import { describe, expect, it } from 'vitest';
import type { IgnitionSystemSpec, SparkPlugSpec } from '../core/engine-spec';
import gasFixture from '../../../test/fixtures/ignition_kernel_gas.json';
import { breakdownVoltage, IgnitionSystem, type IgnitionGasState, type IgnitionSystemOptions } from './index';

/**
 * CFR-like inductive ignition (the numbers of CFR_F1.ignition at the time of writing — kept
 * local so these tests do not drift with engines/cfr.ts) and the ASTM D2699 plug gap.
 */
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

/**
 * Plausible flame inputs at the CFR-like spark states of the Cantera fixture (the laminar
 * flame-speed module is injected in the real simulator). Iso-octane at ≈ 11 bar, 610–640 K:
 * S_L(φ = 1) ≈ 0.6 m/s, 0.28 m/s at φ = 0.7, 0.15 m/s at φ = 0.5 (Metghalchi–Keck-type
 * scaling); Markstein numbers Ma = L/l_F ≈ 3 (φ = 1) rising for lean iso-octane (Le > 1);
 * l_F = ν_u/S_L; u′ ≈ ½ mean piston speed (600 rpm → 1.1 m/s); l_I = 2 mm.
 */
function gasAt(phi: number, over: Partial<IgnitionGasState> = {}): IgnitionGasState {
  const c = gasFixture.cases.find((k) => k.cr === 7 && k.phi === phi)!;
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

interface FireResult {
  ign: IgnitionSystem;
  handoff: number;
  misfire: boolean;
  energy: number;
  duration: number;
}

/** Charge for `dwell`, fire, and run `tRun` after switch-off with a fixed engine step. */
function fire(gas: IgnitionGasState, dwell = 3e-3, tRun = 6e-3, dt = 20e-6, opts: IgnitionSystemOptions = {}): FireResult {
  const ign = new IgnitionSystem(COIL, PLUG, opts);
  for (let t = 0; t < dwell - 1e-12; t += dt) ign.stepTimed(dt, true, gas);
  for (let t = 0; t < tRun - 1e-12; t += dt) ign.stepTimed(dt, false, gas);
  const s = ign.state;
  return {
    ign,
    handoff: s.kernel.handoffTime,
    misfire: s.misfire,
    energy: s.energyDelivered,
    duration: s.sparkDuration,
  };
}

describe('IgnitionSystem — CFR-like spark (CR 7, φ = 1, 13° BTDC, 3 ms dwell)', () => {
  const r = fire(gasAt(1));
  const s = r.ign.state;

  it('delivers 30–100 mJ to the gap over a 1–3 ms discharge', () => {
    expect(r.energy).toBeGreaterThan(30e-3);
    expect(r.energy).toBeLessThan(100e-3);
    expect(r.duration).toBeGreaterThan(1e-3);
    expect(r.duration).toBeLessThan(3.2e-3);
    expect(s.breakdownCount).toBe(1);
    expect(s.reignitionCount).toBe(0); // the inductive coil has no re-ignition model (discharge.ts)
    // breakdown ≈ 10 µs after switch-off at the density-dependent voltage
    expect(s.breakdownDelay).toBeGreaterThan(2e-6);
    expect(s.breakdownDelay).toBeLessThan(40e-6);
    // first breakdown at the unburned-gas density; afterwards the gap sits in hot kernel gas
    const g = gasAt(1);
    expect(s.lastBreakdownVoltage / breakdownVoltage(PLUG.gap, g.p, g.Tu)).toBeCloseTo(1, 2);
    expect(s.breakdownVoltage).toBeLessThan(0.5 * s.lastBreakdownVoltage);
    expect(s.phase).toBe('done');
  });

  it('splits the spark energy: breakdown ≈ 0.3–1 mJ, gas share ≈ 30–50 % (Maly & Vogel)', () => {
    expect(s.energyBreakdown).toBeGreaterThan(0.2e-3);
    expect(s.energyBreakdown).toBeLessThan(1.5e-3);
    const eta = s.energyToGas / s.energyDelivered;
    expect(eta).toBeGreaterThan(0.3);
    expect(eta).toBeLessThan(0.5);
    expect(s.energyToGas + s.energyToElectrodes + s.energyRadiated).toBeCloseTo(s.energyDelivered, 12);
  });

  it('closes the circuit energy balance to round-off over the whole event', () => {
    expect(Math.abs(r.ign.coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('grows the kernel to hand-off in ≈ 0.3–1 ms', () => {
    expect(s.kernel.stage).toBe('handoff');
    expect(r.handoff).toBeGreaterThan(0.3e-3);
    expect(r.handoff).toBeLessThan(1.0e-3);
    expect(s.kernel.radius).toBeGreaterThanOrEqual(s.kernel.handoffRadius);
    expect(s.misfire).toBe(false);
  });
});

describe('IgnitionSystem — misfire and turbulence sensitivity', () => {
  it('misfires for a φ = 0.5 charge with low spark energy (short dwell)', () => {
    const r = fire(gasAt(0.5), 0.5e-3);
    expect(r.energy).toBeLessThan(10e-3);
    expect(r.misfire).toBe(true);
    expect(r.ign.state.kernel.stage).toBe('quenched');
  });

  it('ignites φ = 0.7 with the full dwell but misfires with a weak spark', () => {
    expect(fire(gasAt(0.7), 3e-3).misfire).toBe(false);
    expect(fire(gasAt(0.7), 1e-3).misfire).toBe(true);
  });

  it('grows faster with higher u′ once the kernel is larger than the smallest eddies', () => {
    const lo = fire(gasAt(1, { uPrime: 0.5 }), 3e-3, 6e-3, 20e-6, { kernel: { minHandoffRadius: 4e-3 } });
    const hi = fire(gasAt(1, { uPrime: 3.0 }), 3e-3, 6e-3, 20e-6, { kernel: { minHandoffRadius: 4e-3 } });
    expect(hi.handoff).toBeLessThan(lo.handoff);
  });

  it('never breaks down with breakdown suppressed → misfire, open-circuit peak 20–45 kV', () => {
    const r = fire(gasAt(1), 3e-3, 3e-3, 20e-6, { suppressBreakdown: true });
    expect(r.ign.state.breakdownCount).toBe(0);
    expect(r.ign.state.peakSecondaryVoltage).toBeGreaterThan(20e3);
    expect(r.ign.state.peakSecondaryVoltage).toBeLessThan(45e3);
  });
});

describe('IgnitionSystem — crank-angle driven events', () => {
  it('closes the switch at dwellStartDeg, fires at sparkDeg and reports the phase sequence', () => {
    const ign = new IgnitionSystem(COIL, PLUG);
    const gas = gasAt(1);
    const rpm = 600;
    const degPerS = rpm * 6;
    const dt = 20e-6;
    const cmd = { dwellStartDeg: -13 - COIL.dwellTime * degPerS, sparkDeg: -13 };
    let theta = -60;
    ign.step(dt, theta, cmd, gas);
    const phases: string[] = [];
    let sparkAt = NaN;
    for (let i = 0; i < 2000; i++) {
      theta += degPerS * dt;
      if (theta >= 360) theta -= 720;
      const s = ign.step(dt, theta, cmd, gas);
      if (phases[phases.length - 1] !== s.phase) phases.push(s.phase);
      if (Number.isNaN(sparkAt) && s.breakdownCount > 0) sparkAt = theta;
    }
    expect(phases.slice(0, 4)).toEqual(['off', 'charging', 'breakdown', 'glow']);
    expect(phases).toContain('done');
    // breakdown within ~0.5° after the commanded spark angle
    expect(sparkAt).toBeGreaterThan(-13);
    expect(sparkAt).toBeLessThan(-12.5);
    // stored energy at switch-off ≈ ½ L₁ I₁²(3 ms)
    const i1 = (13.5 / 1.5) * (1 - Math.exp((-1.5 * 3e-3) / 4e-3));
    expect(ign.state.energyAtSwitchOff / (0.5 * 4e-3 * i1 * i1)).toBeCloseTo(1, 2);
  });
});

describe('IgnitionSystem — robustness (reviewer regressions)', () => {
  /** Crank-driven run over `cycles` engine cycles; returns the system and whether every output stayed finite. */
  function crankRun(
    coil: IgnitionSystemSpec,
    plug: SparkPlugSpec,
    gas: IgnitionGasState,
    rpm: number,
    dwell: number,
    dt: number,
    cycles = 1,
  ): { ign: IgnitionSystem; finite: boolean; firstBreakdownOvershoot: number } {
    const ign = new IgnitionSystem(coil, plug);
    const degPerS = rpm * 6;
    const cmd = { dwellStartDeg: -20 - dwell * degPerS, sparkDeg: -20 };
    let theta = -359;
    let finite = true;
    let firstBreakdownOvershoot = 0;
    const n = Math.ceil((cycles * 720) / (degPerS * dt));
    for (let i = 0; i < n; i++) {
      theta += degPerS * dt;
      if (theta >= 360) theta -= 720;
      const nb = ign.gap.breakdownCount;
      const s = ign.step(dt, theta, cmd, gas);
      // first breakdown of an event: voltage at breakdown vs the threshold in force (the
      // gap parameters are only updated at the start of the next step)
      if (nb === 0 && s.breakdownCount > 0 && ign.gap.breakdownCount === 1)
        firstBreakdownOvershoot = s.lastBreakdownVoltage / ign.gap.threshold(1) - 1;
      const k = s.kernel;
      finite &&= [s.primaryCurrent, s.secondaryVoltage, s.secondaryCurrent, s.energyDelivered, s.breakdownVoltage, s.gasPower, s.energyToGas, k.radius, k.volume, k.burnedMass, k.stretchFactor, k.turbulentSpeed, k.growthRate, k.karlovitz].every(Number.isFinite);
    }
    return { ign, finite, firstBreakdownOvershoot };
  }
  const air = Float64Array.from(gasFixture.cases[0].X);
  const lowDensityGas = (p: number, Tu: number, over: Partial<IgnitionGasState> = {}): IgnitionGasState => ({
    p,
    Tu,
    rhoU: (p * 0.0303) / (8.314 * Tu),
    X: air,
    SL: 0.4,
    marksteinLength: 2e-5,
    flameThickness: 5e-5,
    uPrime: 0.5,
    integralScale: 2e-3,
    expansionRatio: 6,
    kinematicViscosity: 2e-5,
    ...over,
  });

  it('finishes a make-spark in hot, low-density gas (regression: zero-time event loop that never returned)', () => {
    // 0.11 mm gap at 1.75 bar: the dwell's make voltage (≈ 1.3 kV) breaks the gap down, the
    // kernel heats it, and the secondary current (≈ 1 mA) is below the sustaining current.
    // (exact parameters of the random-sweep case that hung the previous revision)
    const coil = { ...COIL, secondaryCapacitance: 2.185726296641994e-11, couplingCoefficient: 0.9896605739020743, primaryCurrentLimit: 7.187444926472381, secondaryResistance: 5093.491414586463 };
    const plug = { ...PLUG, gap: 1.074793616460983e-4, centerElectrodeDiameter: 6.851047457728031e-4, groundElectrodeWidth: 2.185619159201956e-3 };
    const gas = lowDensityGas(174631.2759462537, 401.7235832102597, {
      X: Float64Array.from(gasFixture.cases[0].X),
      SL: 2.243711099879224,
      marksteinLength: 4.34894866267039e-5,
      flameThickness: 6.420390632334357e-6,
      uPrime: 0,
      integralScale: 5.304759899616837e-4,
      expansionRatio: 7.385446969419718,
      kinematicViscosity: 1.4405501727329186e-5,
      flowVelocity: 0,
      lewisNumber: 2.4995371089316905,
    });
    const t0 = performance.now();
    const { ign, finite } = crankRun(coil, plug, gas, 1885.2191650075838, 5.976221693633124e-3, 3.549288108885552e-5);
    expect(performance.now() - t0).toBeLessThan(5000);
    expect(finite).toBe(true);
    expect(ign.state.breakdownCount).toBeGreaterThan(0);
    expect(Math.abs(ign.coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('keeps restrikes bounded in hot, low-density gas with a slow gap flow (regression: 46 814 breakdowns)', () => {
    const coil = { ...COIL, secondaryCapacitance: 46e-12, couplingCoefficient: 0.92, primaryCurrentLimit: 10, secondaryResistance: 4.1e3 };
    const plug = { ...PLUG, gap: 0.18e-3 };
    const { ign, finite } = crankRun(coil, plug, lowDensityGas(0.55e5, 676, { flowVelocity: 0.4 }), 5750, 2.9e-3, 1.35e-6);
    expect(finite).toBe(true);
    expect(ign.state.breakdownCount).toBeGreaterThan(0);
    expect(ign.state.breakdownCount).toBeLessThan(1000);
    expect(Math.abs(ign.coil.energyResidual())).toBeLessThan(1e-12);
  });

  it('hand-off time is independent of the engine step (1–200 µs) to < 0.5 %', () => {
    const t: number[] = [];
    for (const dt of [1e-6, 20e-6, 50e-6, 200e-6]) t.push(fire(gasAt(1), 3e-3, 6e-3, dt).handoff);
    for (const x of t) expect(Math.abs(x / t[0] - 1)).toBeLessThan(5e-3);
  });

  it('seeded random sweep (p, T, gap, coil, flame inputs incl. zeros, dt): finite outputs, exact ledger', () => {
    let seed = 20260930;
    const rnd = (): number => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    const lerp = (a: number, b: number): number => a + (b - a) * rnd();
    const logu = (a: number, b: number): number => Math.exp(lerp(Math.log(a), Math.log(b)));
    for (let trial = 0; trial < 25; trial++) {
      const c = gasFixture.cases[Math.floor(rnd() * gasFixture.cases.length)];
      const p = logu(0.2e5, 100e5);
      const Tu = lerp(280, 1000);
      const SL = rnd() < 0.15 ? 0 : logu(0.01, 3);
      const nu = c.nuU * (Tu / c.Tu) ** 1.7 * (c.p / p);
      const lF = SL > 0 ? nu / SL : 1e-5;
      const gas: IgnitionGasState = {
        p,
        Tu,
        rhoU: (p * 0.0303) / (8.314 * Tu),
        X: Float64Array.from(c.X),
        SL,
        marksteinLength: lerp(-3, 12) * lF,
        flameThickness: lF,
        uPrime: rnd() < 0.15 ? 0 : logu(0.01, 20),
        integralScale: rnd() < 0.15 ? 0 : logu(1e-4, 2e-2),
        expansionRatio: lerp(1, 9),
        kinematicViscosity: nu,
        flowVelocity: rnd() < 0.5 ? 0 : logu(0.1, 30),
        lewisNumber: lerp(0.5, 3),
      };
      const coil = { ...COIL, secondaryCapacitance: logu(20e-12, 200e-12), couplingCoefficient: lerp(0.9, 0.995), primaryCurrentLimit: lerp(3, 12), secondaryResistance: logu(1e3, 2e4) };
      const plug = { ...PLUG, gap: logu(0.1e-3, 2e-3), centerElectrodeDiameter: logu(0.4e-3, 4e-3), groundElectrodeWidth: logu(0.5e-3, 4e-3) };
      const dt = logu(5e-6, 300e-6);
      const { ign, finite, firstBreakdownOvershoot } = crankRun(coil, plug, gas, lerp(300, 6000), lerp(0.3e-3, 6e-3), dt);
      expect(finite).toBe(true);
      expect(firstBreakdownOvershoot).toBeLessThan(0.02);
      expect(Math.abs(ign.coil.energyResidual())).toBeLessThan(1e-12);
      const s = ign.state;
      expect(s.energyToGas + s.energyToElectrodes + s.energyRadiated).toBeCloseTo(s.energyDelivered, 12);
      // re-ignition is a trembler-coil event only: the inductive path is unchanged (CFR bit-identical)
      expect(s.reignitionCount).toBe(0);
      expect(ign.gap.energyReignition).toBe(0);
    }
  });
});

describe('IgnitionSystem — high-current coil (arc phase; reviewer regression)', () => {
  it('breaks down at the threshold despite primary-clamp chatter, then arc → glow → extinction', () => {
    // turns ratio ≈ 32: I₂ ≈ 0.18 A > the 0.1 A arc–glow current. Regression: clamp events
    // at I₁ ≈ 0 shadowed the breakdown check and |V₂| overshot V_bd by 2.1 kV (21 %).
    const coil: IgnitionSystemSpec = { ...COIL, secondaryInductance: 4, secondaryResistance: 2e3 };
    const ign = new IgnitionSystem(coil, PLUG);
    const gas = gasAt(1);
    const dt = 20e-6;
    for (let t = 0; t < 3e-3 - 1e-12; t += dt) ign.stepTimed(dt, true, gas);
    const phases: string[] = [];
    for (let t = 0; t < 8e-3; t += dt) {
      const s = ign.stepTimed(dt, false, gas);
      if (phases[phases.length - 1] !== s.phase) phases.push(s.phase);
    }
    const s = ign.state;
    const vbd = breakdownVoltage(PLUG.gap, gas.p, gas.Tu);
    expect(s.breakdownCount).toBe(1);
    expect(s.peakSecondaryVoltage / vbd).toBeLessThan(1.01);
    expect(s.lastBreakdownVoltage / vbd).toBeLessThan(1.01);
    expect(phases).toEqual(['breakdown', 'arc', 'glow', 'done']);
    expect(ign.gap.energyArc).toBeGreaterThan(1e-3);
    expect(ign.gap.energyGlow).toBeGreaterThan(1e-3);
    expect(ign.gap.awaitingReentry).toBe(false);
    expect(Math.abs(ign.coil.energyResidual())).toBeLessThan(1e-12);
  });
});
