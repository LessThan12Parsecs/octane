import { describe, expect, it } from 'vitest';
import type { SparkPlugSpec } from '../core/engine-spec';
import gasFixture from '../../../test/fixtures/ignition_kernel_gas.json';
import { mixHMass } from '../thermo/mixture';
import {
  abdelGayedKarlovitz,
  electrodeContactArea,
  herwegMalyKernelSpeed,
  herwegMalyStretchFactor,
  MIN_INTEGRAL_SCALE,
  SparkKernel,
  turbulentQuenchReached,
  type KernelGasState,
} from './kernel';

const PLUG: SparkPlugSpec = {
  gapCenter: [0, -0.002, 0.04],
  gap: 0.508e-3,
  centerElectrodeDiameter: 2.5e-3,
  groundElectrodeWidth: 2.5e-3,
  axis: [0, 0, -1],
  threadDiameter: 18e-3,
};

function fixtureGas(i: number, over: Partial<KernelGasState> = {}): KernelGasState {
  const c = gasFixture.cases[i];
  const SL = over.SL ?? 0.6;
  return {
    p: c.p,
    Tu: c.Tu,
    rhoU: c.rhoU,
    X: Float64Array.from(c.X),
    SL,
    marksteinLength: 25e-6,
    flameThickness: c.nuU / SL,
    uPrime: 1.0,
    integralScale: 2e-3,
    expansionRatio: c.sigma,
    kinematicViscosity: c.nuU,
    ...over,
  };
}

describe('Herweg–Maly stretch factor and kernel speed', () => {
  it('I₀ = 1 − Ma(Ka_strain + Ka_curv) and vanishes at the critical radius', () => {
    const SL = 0.5;
    const lF = 10e-6;
    const LM = 30e-6;
    const up = 1.5;
    const lI = 2e-3;
    const sig = 4.5;
    const kaS = Math.sqrt(lF / (15 * lI)) * Math.pow(up / SL, 1.5);
    const rc = (2 * sig * LM) / (1 - (LM / lF) * kaS);
    expect(herwegMalyStretchFactor(rc, SL, lF, LM, up, lI, sig)).toBeCloseTo(0, 12);
    expect(herwegMalyStretchFactor(2 * rc, SL, lF, LM, up, lI, sig)).toBeGreaterThan(0);
    expect(herwegMalyStretchFactor(0.5 * rc, SL, lF, LM, up, lI, sig)).toBeLessThan(0);
    // no stretch sensitivity without a Markstein length
    expect(herwegMalyStretchFactor(1e-4, SL, lF, 0, up, lI, sig)).toBe(1);
  });

  it('S_T,k → S_L I₀ in quiescent gas and → the fully developed HM value for r, t ≫ l_I/u′', () => {
    expect(herwegMalyKernelSpeed(0.7, 1e-3, 1e-3, 0.4, 0, 2e-3)).toBeCloseTo(0.28, 12);
    const I0 = 0.8;
    const SL = 0.4;
    const up = 2;
    const full = SL * (I0 + Math.sqrt(I0) * Math.sqrt(up / (up + SL)) * Math.pow(up / SL, 5 / 6));
    expect(herwegMalyKernelSpeed(I0, 1, 10, SL, up, 2e-3)).toBeCloseTo(full, 9);
    // only eddies smaller than the kernel act: small kernels feel much less turbulence
    const small = herwegMalyKernelSpeed(I0, 0.2e-3, 10, SL, up, 2e-3);
    const large = herwegMalyKernelSpeed(I0, 4e-3, 10, SL, up, 2e-3);
    expect(large).toBeGreaterThan(small * 1.5);
  });

  it('Abdel-Gayed–Bradley–Lung (1989) Karlovitz number and quench limits', () => {
    const nu = 5e-6;
    // K = 0.157 (u'/S_L)² R_L^-1/2
    expect(abdelGayedKarlovitz(2, 0.5, 2e-3, nu)).toBeCloseTo((0.157 * 16) / Math.sqrt(800), 12);
    // R_L = 800 > 300: quench when K·Le > 1.5
    expect(turbulentQuenchReached(2, 0.5, 2e-3, nu, 1)).toBe(false); // K = 0.089
    expect(turbulentQuenchReached(8, 0.3, 2e-3, nu, 1)).toBe(true); // K ≈ 1.4·… > 1.5
    // R_L < 300: K R_L^-1/2 > 0.079
    expect(turbulentQuenchReached(0.5, 0.05, 1e-3, nu, 1)).toBe(true);
    expect(turbulentQuenchReached(0.5, 0.5, 1e-3, nu, 1)).toBe(false);
  });
});

describe('electrode contact area', () => {
  it('is zero until the kernel reaches the electrode faces, then two discs 2πρ²', () => {
    expect(electrodeContactArea(0.25e-3, 0.508e-3, 2.5e-3, 2.5e-3)).toBe(0);
    const r = 0.6e-3;
    const rho2 = r * r - 0.25e-3 * 0.25e-3;
    expect(electrodeContactArea(r, 0.5e-3, 2.5e-3, 2.5e-3)).toBeCloseTo(2 * Math.PI * rho2, 15);
  });

  it('adds the centre-electrode flank and the strap strip for large kernels', () => {
    const r = 3e-3;
    const z0 = 0.25e-3;
    const a = 1.25e-3;
    const rho = Math.sqrt(r * r - z0 * z0);
    const face = Math.PI * a * a;
    const flank = 2 * Math.PI * a * (Math.sqrt(r * r - a * a) - z0);
    const x = 1.25e-3;
    const strip = 2 * (x * Math.sqrt(rho * rho - x * x) + rho * rho * Math.asin(x / rho));
    expect(electrodeContactArea(r, 0.5e-3, 2.5e-3, 2.5e-3)).toBeCloseTo(face + flank + strip, 15);
  });
});

describe('SparkKernel', () => {
  it('uses the same h_u(T) as Cantera for the plasma-velocity denominator and T_k from σ', () => {
    for (const c of gasFixture.cases) {
      const X = Float64Array.from(c.X);
      const dh = mixHMass(X, c.Tad) - mixHMass(X, c.Tu);
      expect(Math.abs(dh / c.dhHeat - 1)).toBeLessThan(1e-6);
      const k = new SparkKernel(PLUG);
      const gas = fixtureGas(gasFixture.cases.indexOf(c));
      k.deposit(1e-3, gas);
      // kernel temperature from σ with complete-combustion M_b vs Cantera HP equilibrium T_ad
      expect(Math.abs(k.temperature / c.Tad - 1)).toBeLessThan(0.02);
      // initial kernel volume = E/(ρ_k Δh(T_k))
      const dhk = mixHMass(X, k.temperature) - mixHMass(X, c.Tu);
      expect(k.volume / (1e-3 / ((c.rhoU / c.sigma) * dhk))).toBeCloseTo(1, 12);
    }
  });

  it('inert kernel under constant power grows exactly as V₀ + P t/(ρ_k Δh) (Herweg–Maly plasma term)', () => {
    const c = gasFixture.cases[0];
    const gas = fixtureGas(0, { SL: 0, expansionRatio: 1, flameThickness: 1e-5 });
    const k = new SparkKernel(PLUG, { electrodeHeatLoss: false });
    k.deposit(1e-3, gas);
    const V0 = k.volume;
    const Tk = 3000;
    const rhoK = (c.rhoU * c.Tu) / Tk;
    const dh = mixHMass(gas.X, Tk) - mixHMass(gas.X, c.Tu);
    for (let i = 0; i < 50; i++) k.step(20e-6, 10, true, gas);
    expect(k.volume / (V0 + (10 * 1e-3) / (rhoK * dh))).toBeCloseTo(1, 10);
    expect(k.stage).toBe('kernel');
  });

  it('a laminar kernel grows at dr/dt = σ S_L I₀(r) without spark or losses', () => {
    const gas = fixtureGas(0, { uPrime: 0 });
    const k = new SparkKernel(PLUG, { electrodeHeatLoss: false, minHandoffRadius: 1 });
    k.deposit(0.3, gas); // large kernel (several mm)
    k.step(1e-6, 0, false, gas);
    const I0 = 1 - (2 * gas.marksteinLength * gas.expansionRatio) / k.radius;
    expect(k.stretchFactor).toBeCloseTo(I0, 4);
    expect(k.growthRate / (gas.expansionRatio * gas.SL * I0)).toBeCloseTo(1, 4);
  });

  it('is extinguished after the spark if it is below its critical radius (misfire)', () => {
    const gas = fixtureGas(0, { marksteinLength: 0.2e-3 });
    const k = new SparkKernel(PLUG);
    k.deposit(0.5e-3, gas);
    expect(k.radius).toBeLessThan(2 * gas.expansionRatio * gas.marksteinLength);
    // (fixer round 2: a non-propagating kernel survives up to min(l_I/u′ = 2 ms, its conductive cooling
    // time r²/α_k) before it is extinguished — round 1 extinguished it in the first step)
    k.step(20e-6, 0, false, gas);
    expect(k.stage).toBe('kernel');
    let n = 0;
    while (k.stage === 'kernel' && n++ < 1000) k.step(20e-6, 0, false, gas);
    expect(k.stage).toBe('quenched');
    expect(k.misfire).toBe(true);
    expect(k.quenchReason).toMatch(/critical radius/);
    expect(k.age).toBeGreaterThan(0.2e-3);
    expect(k.age).toBeLessThan(gas.integralScale / gas.uPrime + 40e-6);
    // round-1 behaviour on request
    const k1 = new SparkKernel(PLUG, { nonPropagatingSurvival: 0 });
    k1.deposit(0.5e-3, gas);
    k1.step(20e-6, 0, false, gas);
    expect(k1.stage).toBe('quenched');
  });

  it('hands off at r = max(1 mm, l_I) with a positive stretch factor', () => {
    const gas = fixtureGas(0, { integralScale: 1.5e-3 });
    const k = new SparkKernel(PLUG);
    k.deposit(2e-3, gas);
    let n = 0;
    while (k.stage === 'kernel' && n++ < 1000) k.step(10e-6, 5, true, gas);
    expect(k.stage).toBe('handoff');
    expect(k.radius).toBeGreaterThanOrEqual(1.5e-3);
    expect(k.radius).toBeLessThan(1.6e-3);
    expect(k.handoffTime).toBeGreaterThan(0);
    expect(k.burnedMass).toBeCloseTo((gas.rhoU / gas.expansionRatio) * k.volume, 15);
  });
});

describe('SparkKernel robustness and step independence (reviewer fixes)', () => {
  const finiteState = (k: SparkKernel): boolean =>
    [k.radius, k.volume, k.burnedMass, k.stretchFactor, k.turbulentSpeed, k.plasmaSpeed, k.growthRate, k.karlovitzAGB, k.energyElectrodeLoss].every(
      Number.isFinite,
    );

  it('pure helpers take their limits instead of 0/0 (u′ = 0, l_I = 0, L_M = 0, l_F = 0)', () => {
    expect(abdelGayedKarlovitz(0, 0.5, 0, 5e-6)).toBe(0);
    expect(abdelGayedKarlovitz(1, 0.5, 0, 5e-6)).toBe(Infinity);
    expect(herwegMalyKernelSpeed(0.8, 1e-3, 0, 0.5, 0, 0)).toBeCloseTo(0.4, 12);
    const full = 0.5 * (0.8 + Math.sqrt(0.8) * Math.sqrt(2 / 2.5) * Math.pow(2 / 0.5, 5 / 6));
    expect(herwegMalyKernelSpeed(0.8, 1e-3, 0, 0.5, 2, 0)).toBeCloseTo(full, 12);
    expect(herwegMalyStretchFactor(1e-4, 0.5, 0, 0, 2, 2e-3, 4)).toBe(1);
    // division-free form ≡ 1 − (L_M/l_F)(Ka_s + Ka_c)
    const lF = 8e-6, LM = 2e-5, up = 1.3, lI = 2e-3, SL = 0.5, sig = 4.4, r = 7e-4;
    const kaS = Math.sqrt(lF / (15 * lI)) * Math.pow(up / SL, 1.5);
    expect(herwegMalyStretchFactor(r, SL, lF, LM, up, lI, sig)).toBeCloseTo(1 - (LM / lF) * (kaS + (2 * lF * sig) / r), 12);
  });

  it('quiescent charge (u′ = 0) with l_I = 0 is finite and identical to any l_I > 0 at the same hand-off radius', () => {
    const run = (lI: number): SparkKernel => {
      const gas = fixtureGas(0, { uPrime: 0, integralScale: lI });
      const k = new SparkKernel(PLUG, { handoffIntegralScaleMultiple: 0 });
      k.deposit(1.5e-3, gas);
      for (let i = 0; i < 200 && k.stage === 'kernel'; i++) {
        k.step(20e-6, 5, true, gas);
        expect(finiteState(k)).toBe(true);
      }
      return k;
    };
    const a = run(0); // regression: NaN radius/volume before the fix
    const b = run(2e-3);
    expect(a.stage).toBe('handoff');
    expect(a.handoffTime).toBeCloseTo(b.handoffTime, 15);
    expect(MIN_INTEGRAL_SCALE).toBeLessThan(1e-5);
  });

  it('hand-off time does not depend on the caller step (located inside the sub-step)', () => {
    const times: number[] = [];
    for (const dt of [1e-6, 10e-6, 37e-6, 100e-6, 500e-6]) {
      const gas = fixtureGas(0, { integralScale: 1.5e-3 });
      const k = new SparkKernel(PLUG);
      k.deposit(1.5e-3, gas);
      while (k.stage === 'kernel') k.step(dt, 6, true, gas);
      expect(k.stage).toBe('handoff');
      expect(k.radius).toBeCloseTo(k.handoffRadius, 12);
      times.push(k.handoffTime);
    }
    for (const t of times) expect(Math.abs(t / times[0] - 1)).toBeLessThan(2e-3);
  });

  it('books the spark energy deposited up to the hand-off instant (E₀ + P·t_handoff)', () => {
    const gas = fixtureGas(0, { integralScale: 1.5e-3 });
    const k = new SparkKernel(PLUG);
    k.deposit(1.5e-3, gas);
    while (k.stage === 'kernel') k.step(100e-6, 6, true, gas);
    // regression: the partial step up to the hand-off was dropped from the diagnostic
    expect(k.energyDeposited).toBeCloseTo(1.5e-3 + 6 * k.handoffTime, 12);
  });

  it('uses one kinematic viscosity for both strain measures (the caller ν wins over l_F)', () => {
    const c = gasFixture.cases[0];
    const a = new SparkKernel(PLUG);
    const b = new SparkKernel(PLUG);
    const ga = fixtureGas(0);
    const gb = fixtureGas(0, { flameThickness: 1.4 * ga.flameThickness }); // e.g. α_u/S_L passed by mistake
    a.deposit(1.5e-3, ga);
    b.deposit(1.5e-3, gb);
    a.step(20e-6, 5, true, ga);
    b.step(20e-6, 5, true, gb);
    expect(b.karlovitzStrain).toBeCloseTo(a.karlovitzStrain, 15);
    expect(b.radius).toBeCloseTo(a.radius, 15);
    // without an explicit ν the caller's l_F defines it: ν = l_F S_L
    const gc = { ...fixtureGas(0), kinematicViscosity: undefined };
    const k = new SparkKernel(PLUG);
    k.deposit(1.5e-3, gc);
    k.step(20e-6, 5, true, gc);
    expect(k.karlovitzAGB).toBeCloseTo(abdelGayedKarlovitz(gc.uPrime, gc.SL, gc.integralScale, c.nuU), 12);
  });
});
