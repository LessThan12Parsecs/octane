import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import {
  MOCK_GAS_CONSTANT,
  MockKinematics,
  MockSimulator,
  type MockSimulatorOptions,
  MockValveLift,
  breakdownVoltage,
  burnedGasComposition,
  douaudEyzatDelay,
  mockFuel,
  orificeMassFlow,
  sphereCoverRadius,
  sphereDiscFrontArea,
  sphereDiscVolume,
  sphereRadiusForVolume,
  unburnedMixtureComposition,
} from './mock-simulator';

const OP: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: 101325,
  ambientTemperature: 298,
  relativeHumidity: 0.5,
  intakeMixtureTemperature: 325,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: 1.0,
  sparkAdvanceDeg: 13,
  dwellTime: 3e-3,
  compressionRatio: 7,
  egrFraction: 0,
  coolantTemperature: 373,
};
const OPTS: MockSimulatorOptions = { snapshotEveryDeg: 0.5, bufferAheadSeconds: 0.25 };

function simulate(
  patch: Partial<OperatingPoint> = {},
  cycles = 2,
  opts: Partial<MockSimulatorOptions> = {},
): { snaps: EngineSnapshot[]; summaries: CycleSummary[]; sim: MockSimulator } {
  const op = { ...OP, ...patch };
  const sim = new MockSimulator(CFR_F1, op, { ...OPTS, ...opts });
  const snaps: EngineSnapshot[] = [];
  const summaries: CycleSummary[] = [];
  for (;;) {
    const s = sim.advanceToNextSnapshot();
    snaps.push(s);
    summaries.push(...sim.drainCycleSummaries());
    if (s.cycle >= cycles) break;
  }
  return { snaps, summaries, sim };
}

const cycleSnaps = (snaps: EngineSnapshot[], c: number): EngineSnapshot[] => snaps.filter((s) => s.cycle === c);

describe('MockKinematics', () => {
  const g = CFR_F1.geometry;
  const k = new MockKinematics(g, 7);
  it('volume and displacement at the dead centres', () => {
    expect(k.volume(0)).toBeCloseTo(k.clearanceVolume, 12);
    expect(k.volume(-360)).toBeCloseTo(k.clearanceVolume, 12);
    expect(k.volume(180)).toBeCloseTo(k.clearanceVolume + k.displacedVolume, 12);
    expect(k.pistonDisplacement(180)).toBeCloseTo(g.stroke, 12);
    expect(k.displacedVolume / k.clearanceVolume + 1).toBeCloseTo(7, 10);
    expect(k.clearanceHeight(0)).toBeCloseTo((k.clearanceVolume - g.creviceVolume) / k.pistonArea, 12);
  });
  it('derivatives match finite differences', () => {
    const h = 1e-4; // deg
    for (const th of [-300, -120, -30, 10, 77, 200]) {
      const fd = (k.pistonDisplacement(th + h) - k.pistonDisplacement(th - h)) / (2 * h * (Math.PI / 180));
      expect(k.dxdTheta(th)).toBeCloseTo(fd, 7);
      const fd2 = (k.dxdTheta(th + h) - k.dxdTheta(th - h)) / (2 * h * (Math.PI / 180));
      expect(k.d2xdTheta2(th)).toBeCloseTo(fd2, 5);
    }
  });
  it('rod angle', () => {
    expect(k.rodAngle(90)).toBeCloseTo(Math.asin(g.stroke / 2 / g.conRodLength), 12);
    expect(k.rodAngle(0)).toBeCloseTo(0, 12);
  });
});

describe('MockValveLift', () => {
  it('intake: zero at the spec open/close angles, max mid-way', () => {
    const v = CFR_F1.intakeValve;
    const lift = new MockValveLift({ ...v, timingLiftThreshold: 0 });
    expect(lift.lift(v.openDeg)).toBeCloseTo(0, 12);
    expect(lift.lift(v.closeDeg)).toBeCloseTo(0, 12);
    expect(lift.lift(v.openDeg + 1)).toBeGreaterThan(0);
    const mid = v.openDeg + (v.closeDeg - v.openDeg) / 2;
    expect(lift.lift(mid)).toBeCloseTo(v.maxLift, 9);
    expect(lift.lift(0)).toBe(0);
  });
  it('exhaust handles the 720° wrap (EVC after gas-exchange TDC)', () => {
    const v = { ...CFR_F1.exhaustValve, openDeg: 140, closeDeg: -345, timingLiftThreshold: 0 };
    const lift = new MockValveLift(v);
    expect(lift.lift(139)).toBe(0);
    expect(lift.lift(150)).toBeGreaterThan(0);
    expect(lift.lift(359)).toBeGreaterThan(0);
    expect(lift.lift(-359)).toBeGreaterThan(0);
    expect(lift.lift(-346)).toBeGreaterThan(0);
    expect(lift.lift(-344)).toBe(0);
    expect(lift.lift(0)).toBe(0);
  });
  it('timing quoted at a lift threshold', () => {
    const v = { ...CFR_F1.intakeValve, timingLiftThreshold: 0.1e-3 };
    const lift = new MockValveLift(v);
    expect(lift.lift(v.openDeg)).toBeCloseTo(0.1e-3, 9);
    expect(lift.lift(v.closeDeg)).toBeCloseTo(0.1e-3, 9);
    expect(lift.lift(v.openDeg - 5)).toBeGreaterThan(0);
  });
});

describe('mock sub-models', () => {
  it('orifice flow: zero without Δp, choked flow independent of downstream', () => {
    expect(orificeMassFlow(1e-4, 1e5, 300, 1e5, 1.4)).toBe(0);
    expect(orificeMassFlow(1e-4, 1e5, 300, 2e5, 1.4)).toBe(0);
    const a = orificeMassFlow(1e-4, 5e5, 300, 1e5, 1.4);
    const b = orificeMassFlow(1e-4, 5e5, 300, 0.5e5, 1.4);
    expect(a).toBeCloseTo(b, 12);
    expect(orificeMassFlow(1e-4, 1e5, 300, 0.9e5, 1.4)).toBeLessThan(a);
    // Continuous at the linearisation boundary.
    const l = orificeMassFlow(1e-4, 1e5, 300, 0.99e5 - 1e-3, 1.4);
    const r = orificeMassFlow(1e-4, 1e5, 300, 0.99e5 + 1e-3, 1.4);
    expect(Math.abs(l - r) / l).toBeLessThan(1e-5);
  });

  it('sphere ∩ disc geometry', () => {
    const R = 0.04;
    const h = 0.02;
    const c: [number, number, number] = [0, -0.01, 0];
    const r = 0.004;
    expect(sphereDiscVolume(r, c, R, h) / ((4 / 3) * Math.PI * r ** 3)).toBeCloseTo(1, 6);
    expect(sphereDiscFrontArea(r, c, R, h) / (4 * Math.PI * r * r)).toBeCloseTo(1, 4);
    // Centred on the head plane: exactly half inside.
    const cHead: [number, number, number] = [0, 0, 0];
    expect(sphereDiscVolume(r, cHead, R, h) / ((2 / 3) * Math.PI * r ** 3)).toBeCloseTo(1, 6);
    expect(sphereDiscFrontArea(r, cHead, R, h) / (2 * Math.PI * r * r)).toBeCloseTo(1, 4);
    // Cover radius fills the disc.
    const cSide: [number, number, number] = [0.035, -0.002, 0];
    const rc = sphereCoverRadius(cSide, R, h);
    expect(sphereDiscVolume(rc, cSide, R, h) / (Math.PI * R * R * h)).toBeCloseTo(1, 3);
    // Inversion round-trip.
    for (const rr of [0.003, 0.012, 0.03, 0.06]) {
      const V = sphereDiscVolume(rr, cSide, R, h);
      expect(sphereRadiusForVolume(V, cSide, R, h)).toBeCloseTo(rr, 6);
    }
  });

  it('Douaud–Eyzat delay trends', () => {
    const tau = douaudEyzatDelay(90, 30e5, 800);
    expect(tau).toBeGreaterThan(1e-4);
    expect(tau).toBeLessThan(1e-1);
    expect(douaudEyzatDelay(100, 30e5, 800)).toBeGreaterThan(tau);
    expect(douaudEyzatDelay(90, 40e5, 800)).toBeLessThan(tau);
    expect(douaudEyzatDelay(90, 30e5, 850)).toBeLessThan(tau);
  });

  it('breakdown voltage: ~4.6 kV for 1 mm of air at 1 atm, rising with density and gap', () => {
    expect(breakdownVoltage(1e-3, 1.2041)).toBeGreaterThan(4.0e3);
    expect(breakdownVoltage(1e-3, 1.2041)).toBeLessThan(5.2e3);
    expect(breakdownVoltage(1e-3, 6)).toBeGreaterThan(breakdownVoltage(1e-3, 3));
    expect(breakdownVoltage(1.2e-3, 3)).toBeGreaterThan(breakdownVoltage(0.8e-3, 3));
  });

  it('fuels', () => {
    const iso = mockFuel({ kind: 'pure', species: 'IC8H18' });
    expect(iso.afrStoich).toBeGreaterThan(14.9);
    expect(iso.afrStoich).toBeLessThan(15.3);
    const prf = mockFuel({ kind: 'PRF', octaneNumber: 90 });
    expect(prf.octane).toBe(90);
    expect(prf.C).toBeGreaterThan(7);
    expect(prf.C).toBeLessThan(8);
    expect(mockFuel({ kind: 'pure', species: 'C2H5OH' }).afrStoich).toBeLessThan(9.5);
  });

  it('burned composition: normalised, lean → O2, rich → CO/H2, dissociation grows with T', () => {
    const f = mockFuel({ kind: 'PRF', octaneNumber: 90 });
    const z = () => ({ CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0 });
    const sum = (c: ReturnType<typeof z>) => Object.values(c).reduce((a, b) => a + b, 0);
    const lean = burnedGasComposition(f, 0.8, 2400, 30e5, NaN, z());
    const rich = burnedGasComposition(f, 1.3, 2400, 30e5, NaN, z());
    const hot = burnedGasComposition(f, 1.0, 2800, 30e5, NaN, z());
    const cool = burnedGasComposition(f, 1.0, 2000, 30e5, NaN, z());
    for (const c of [lean, rich, hot, cool]) {
      expect(sum(c)).toBeCloseTo(1, 12);
      for (const v of Object.values(c)) expect(v).toBeGreaterThanOrEqual(0);
    }
    expect(lean.O2).toBeGreaterThan(0.02);
    expect(rich.CO).toBeGreaterThan(0.02);
    expect(rich.H2).toBeGreaterThan(0.01);
    expect(hot.OH).toBeGreaterThan(cool.OH);
    expect(hot.CO).toBeGreaterThan(cool.CO);
    expect(hot.NO).toBeGreaterThan(1e-3);
    expect(hot.NO).toBeLessThan(2e-2);
    const unb = unburnedMixtureComposition(f, 1, 0, 0, z());
    expect(sum(unb)).toBeLessThan(1); // remainder = fuel vapour
    expect(sum(unb)).toBeGreaterThan(0.97);
    expect(unb.O2).toBeCloseTo(0.2095 * sum(unb), 2);
    expect(unb.CO2).toBe(0);
  });
});

describe('MockSimulator', () => {
  const { snaps, summaries } = simulate();
  const kin = new MockKinematics(CFR_F1.geometry, OP.compressionRatio);

  it('starts at t = 0, θ = −360, cycle 0 and advances monotonically', () => {
    expect(snaps[0].t).toBe(0);
    expect(snaps[0].thetaDeg).toBe(-360);
    expect(snaps[0].cycle).toBe(0);
    for (let i = 1; i < snaps.length; i++) {
      const a = snaps[i - 1];
      const b = snaps[i];
      expect(b.t).toBeGreaterThan(a.t);
      expect(b.thetaDeg).toBeGreaterThanOrEqual(-360);
      expect(b.thetaDeg).toBeLessThan(360);
      if (b.cycle === a.cycle) {
        expect(b.thetaDeg).toBeGreaterThan(a.thetaDeg);
        expect(b.thetaDeg - a.thetaDeg).toBeLessThanOrEqual(OPTS.snapshotEveryDeg + 1e-9);
      } else {
        expect(b.cycle).toBe(a.cycle + 1);
        expect(b.thetaDeg).toBe(-360);
      }
    }
  });

  it('emits every grid angle once per cycle (plus event refinements)', () => {
    const c1 = cycleSnaps(snaps, 1);
    const grid = new Set(c1.map((s) => Math.round((s.thetaDeg + 360) / 0.5 * 1e6) / 1e6).filter(Number.isInteger));
    expect(grid.size).toBe(1440);
    expect(c1.length).toBeGreaterThan(1440); // spark events
  });

  it('time is consistent with the fixed speed', () => {
    const c1 = cycleSnaps(snaps, 1);
    for (let i = 1; i < c1.length; i++) {
      const dt = c1[i].t - c1[i - 1].t;
      const dth = c1[i].thetaDeg - c1[i - 1].thetaDeg;
      expect(dt).toBeCloseTo(dth / (6 * OP.rpm), 9);
    }
  });

  it('snapshots are internally consistent', () => {
    for (const s of snaps) {
      expect(s.volume).toBeCloseTo(kin.clearanceVolume + kin.pistonArea * s.pistonDisplacement, 12);
      expect(s.clearanceHeight).toBeCloseTo(kin.clearanceHeightTdc + s.pistonDisplacement, 12);
      expect(s.pistonDisplacement).toBeCloseTo(kin.pistonDisplacement(s.thetaDeg), 12);
      const pThermo = s.pressure - s.knock.oscillation;
      expect((pThermo * s.volume) / (s.mass * MOCK_GAS_CONSTANT * s.temperatureMean)).toBeCloseTo(1, 9);
      const mix = (1 - s.massFractionBurned) * s.temperatureUnburned + s.massFractionBurned * s.temperatureBurned;
      expect(mix / s.temperatureMean).toBeCloseTo(1, 9);
      expect(s.massFractionBurned).toBeGreaterThanOrEqual(0);
      expect(s.massFractionBurned).toBeLessThanOrEqual(1);
      if (s.temperatureBurned === 0) expect(s.massFractionBurned).toBe(0);
      const valveOpen = s.intakeLift > 0 || s.exhaustLift > 0;
      // Gas exchange spans EVO → IVC; with a negative valve overlap (EVC before IVO,
      // as in the measured CFR timing) both valves are briefly shut inside it.
      if (valveOpen) expect(s.phase).toBe('gas-exchange');
      if (!valveOpen) {
        expect(s.intakeMassFlow).toBe(0);
        expect(s.exhaustMassFlow).toBe(0);
      }
      expect(s.flame.center).toEqual(CFR_F1.sparkPlug.gapCenter);
      for (const v of Object.values(s.burnedComposition)) {
        expect(v).toBeGreaterThanOrEqual(0);
        expect(v).toBeLessThanOrEqual(1);
      }
    }
  });

  it('burned fraction is monotone within the closed cycle; flame stages progress in order', () => {
    const order = ['none', 'kernel', 'turbulent', 'burnout', 'done'];
    const c1 = cycleSnaps(snaps, 1);
    let xb = 0;
    let stage = 0;
    for (const s of c1) {
      if (s.phase !== 'gas-exchange' || s.thetaDeg > 0) {
        expect(s.massFractionBurned).toBeGreaterThanOrEqual(xb - 1e-12);
        xb = s.massFractionBurned;
      }
      const k = order.indexOf(s.flame.stage);
      expect(k).toBeGreaterThanOrEqual(stage);
      stage = k;
    }
    expect(stage).toBe(4);
    expect(c1.some((s) => s.flame.stage === 'turbulent')).toBe(true);
    // Flame radius grows while burning.
    const burning = c1.filter((s) => s.flame.stage === 'turbulent');
    for (let i = 1; i < burning.length; i++) expect(burning[i].flame.radius).toBeGreaterThanOrEqual(burning[i - 1].flame.radius - 1e-9);
    expect(burning.every((s) => s.flame.area > 0)).toBe(true);
  });

  it('spark phases occur in order with plausible electrical values', () => {
    const order = ['off', 'charging', 'breakdown', 'arc', 'glow', 'done'];
    const c1 = cycleSnaps(snaps, 1);
    let k = 0;
    for (const s of c1) {
      const j = order.indexOf(s.spark.phase);
      expect(j).toBeGreaterThanOrEqual(k);
      k = j;
    }
    for (const ph of order) expect(c1.some((s) => s.spark.phase === ph)).toBe(true);
    const charging = c1.filter((s) => s.spark.phase === 'charging');
    expect(charging[0].thetaDeg).toBeCloseTo(-OP.sparkAdvanceDeg - OP.dwellTime * 6 * OP.rpm, 6);
    const lastCharge = charging[charging.length - 1];
    expect(lastCharge.spark.primaryCurrent).toBeGreaterThan(2);
    expect(lastCharge.spark.primaryCurrent).toBeLessThanOrEqual(CFR_F1.ignition.primaryCurrentLimit);
    const firstBd = c1.find((s) => s.spark.phase === 'breakdown')!;
    expect(firstBd.thetaDeg).toBeCloseTo(-OP.sparkAdvanceDeg, 9);
    const peakV = Math.max(...c1.filter((s) => s.spark.phase === 'breakdown').map((s) => s.spark.secondaryVoltage));
    expect(peakV).toBeGreaterThan(3e3);
    expect(peakV).toBeCloseTo(firstBd.spark.breakdownVoltage, -3);
    const glow = c1.filter((s) => s.spark.phase === 'glow');
    expect(glow[0].spark.secondaryCurrent).toBeGreaterThan(0.02);
    expect(glow[0].spark.secondaryCurrent).toBeLessThan(0.2);
    const glowDur = glow[glow.length - 1].t - glow[0].t;
    expect(glowDur).toBeGreaterThan(0.3e-3);
    expect(glowDur).toBeLessThan(6e-3);
    const done = c1.find((s) => s.spark.phase === 'done')!;
    expect(done.spark.energyDelivered).toBeGreaterThan(5e-3);
    expect(done.spark.energyDelivered).toBeLessThan(0.15);
    // The breakdown and arc are resolved in time (event snapshots within microseconds).
    const arc = c1.filter((s) => s.spark.phase === 'arc');
    expect(arc.length).toBeGreaterThanOrEqual(2);
  });

  it('produces a plausible fired cycle', () => {
    expect(summaries.map((c) => c.cycle)).toEqual([0, 1]);
    const c = summaries[1];
    expect(c.misfire).toBe(false);
    expect(c.imepNet).toBeGreaterThan(5e5);
    expect(c.imepNet).toBeLessThan(12e5);
    expect(c.pmep).toBeGreaterThan(0);
    expect(c.pmep).toBeLessThan(0.5e5);
    expect(c.peakPressure).toBeGreaterThan(25e5);
    expect(c.peakPressure).toBeLessThan(60e5);
    expect(c.ca10).toBeLessThan(c.ca50);
    expect(c.ca50).toBeLessThan(c.ca90);
    expect(c.ca50).toBeGreaterThan(0);
    expect(c.ca50).toBeLessThan(30);
    expect(c.indicatedEfficiency).toBeGreaterThan(0.2);
    expect(c.indicatedEfficiency).toBeLessThan(0.4);
    expect(c.volumetricEfficiency).toBeGreaterThan(0.6);
    expect(c.volumetricEfficiency).toBeLessThan(1.0);
    expect(c.residualFraction).toBeGreaterThan(0.02);
    expect(c.residualFraction).toBeLessThan(0.15);
    expect(c.noPpm).toBeGreaterThan(300);
    expect(c.noPpm).toBeLessThan(5000);
    expect(c.heatLoss).toBeGreaterThan(0.05 * c.fuelMass * 44e6);
    expect(c.heatLoss).toBeLessThan(0.4 * c.fuelMass * 44e6);
    expect(c.indicatedWorkGross).toBeCloseTo(c.imepGross * kin.displacedVolume, 6);
    const c1 = cycleSnaps(snaps, 1);
    const maxTb = Math.max(...c1.map((s) => s.temperatureBurned));
    expect(maxTb).toBeGreaterThan(2200);
    expect(maxTb).toBeLessThan(2900);
    const peak = c1.reduce((a, b) => (b.pressure > a.pressure ? b : a));
    expect(peak.pressure).toBeCloseTo(c.peakPressure, -4);
    // Ignition ~1 ms after the spark (kernel delay) before measurable burning.
    const firstBurn = c1.find((s) => s.massFractionBurned > 0)!;
    expect(firstBurn.t - c1.find((s) => s.spark.phase === 'arc')!.t).toBeGreaterThan(0.5e-3);
  });

  it('is deterministic and reset() restarts identically', () => {
    const a = new MockSimulator(CFR_F1, OP, OPTS);
    const b = new MockSimulator(CFR_F1, OP, OPTS);
    const sa: EngineSnapshot[] = [];
    for (let i = 0; i < 3000; i++) sa.push(a.advanceToNextSnapshot());
    for (let i = 0; i < 3000; i++) expect(b.advanceToNextSnapshot()).toEqual(sa[i]);
    a.reset();
    expect(a.time).toBe(0);
    for (let i = 0; i < 3000; i++) expect(a.advanceToNextSnapshot()).toEqual(sa[i]);
  });

  it('motored (lean misfire) cycle: no burn, peak pressure near TDC', () => {
    const { summaries: s, snaps: sn } = simulate({ equivalenceRatio: 0.4 }, 2);
    const c = s[1];
    expect(c.misfire).toBe(true);
    expect(Number.isNaN(c.ca50)).toBe(true);
    expect(Math.abs(c.peakPressureDeg)).toBeLessThan(5);
    expect(c.imepGross).toBeLessThan(0.3e5);
    expect(cycleSnaps(sn, 1).every((x) => x.massFractionBurned === 0 && x.flame.stage === 'none')).toBe(true);
  });

  it('knocks for a low-octane fuel at high CR, not for a high-octane fuel at low CR', () => {
    const { summaries: hi, snaps: sn } = simulate({ fuel: { kind: 'PRF', octaneNumber: 60 }, compressionRatio: 9 }, 2);
    const k = hi[1];
    expect(Number.isNaN(k.knockOnsetDeg)).toBe(false);
    expect(k.knockEndGasFraction).toBeGreaterThan(0);
    expect(k.mapo).toBeGreaterThan(0.2e5);
    const c1 = cycleSnaps(sn, 1);
    expect(c1.some((s) => s.knock.autoignited)).toBe(true);
    expect(Math.max(...c1.map((s) => Math.abs(s.knock.oscillation)))).toBeGreaterThan(0.1e5);
    // The ringing is resolved: ≥ 8 samples per period of the ~6–8 kHz mode.
    const ringing = c1.filter((s) => s.knock.oscillation !== 0);
    const dt = (ringing[ringing.length - 1].t - ringing[0].t) / (ringing.length - 1);
    expect(1 / dt).toBeGreaterThan(8 * 6000);
    const { summaries: lo } = simulate({ fuel: { kind: 'PRF', octaneNumber: 100 }, compressionRatio: 5 }, 2);
    expect(Number.isNaN(lo[1].knockOnsetDeg)).toBe(true);
    expect(lo[1].mapo).toBe(0);
  });

  it('knock-limited CR rises with octane number', () => {
    const endGas = (on: number, cr: number) =>
      simulate({ fuel: { kind: 'PRF', octaneNumber: on }, compressionRatio: cr }, 2, { cyclicVariability: 0 }).summaries[1]
        .knockEndGasFraction;
    expect(endGas(80, 6.5)).toBeGreaterThan(endGas(90, 6.5));
    expect(endGas(90, 6.5)).toBeGreaterThan(endGas(100, 6.5));
  });

  it('part throttle lowers manifold pressure, trapped mass and IMEP', () => {
    const { summaries: s, snaps: sn } = simulate({ throttle: 0.1 }, 2);
    const c1 = cycleSnaps(sn, 1);
    const pMan = c1.reduce((a, x) => a + x.intakeManifoldPressure, 0) / c1.length;
    expect(pMan).toBeLessThan(0.8 * OP.ambientPressure);
    expect(s[1].trappedMass).toBeLessThan(summaries[1].trappedMass);
    expect(s[1].imepNet).toBeLessThan(summaries[1].imepNet);
    expect(s[1].pmep).toBeGreaterThan(summaries[1].pmep);
  });

  it('applies rpm immediately and spark/CR changes at the next cycle', () => {
    const sim = new MockSimulator(CFR_F1, OP, OPTS);
    let s = sim.advanceToNextSnapshot();
    while (s.thetaDeg < -100) s = sim.advanceToNextSnapshot();
    sim.setOperatingPoint({ rpm: 900, sparkAdvanceDeg: 25, compressionRatio: 8 });
    s = sim.advanceToNextSnapshot();
    expect(s.rpm).toBe(900);
    const cycle0: EngineSnapshot[] = [];
    const cycle1: EngineSnapshot[] = [];
    while (s.cycle < 2) {
      (s.cycle === 0 ? cycle0 : cycle1).push(s);
      s = sim.advanceToNextSnapshot();
    }
    const bd0 = cycle0.find((x) => x.spark.phase === 'breakdown')!;
    const bd1 = cycle1.find((x) => x.spark.phase === 'breakdown')!;
    expect(bd0.thetaDeg).toBeCloseTo(-13, 6);
    expect(bd1.thetaDeg).toBeCloseTo(-25, 6);
    const k8 = new MockKinematics(CFR_F1.geometry, 8);
    const tdc1 = cycle1.find((x) => x.thetaDeg === 0)!;
    expect(tdc1.volume).toBeCloseTo(k8.clearanceVolume, 12);
    const tdc0 = cycle0.find((x) => x.thetaDeg === 0)!;
    expect(tdc0.volume).toBeCloseTo(kin.clearanceVolume, 12);
  });

  it('is fast enough for real-time playback', () => {
    const sim = new MockSimulator(CFR_F1, OP, OPTS);
    const t0 = performance.now();
    while (sim.time < 1.0) sim.advanceToNextSnapshot(); // 5 cycles at 600 rpm
    const wall = (performance.now() - t0) / 1000;
    expect(wall).toBeLessThan(1.0);
  });
});
