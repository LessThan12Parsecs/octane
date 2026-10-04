/**
 * Ford Model T presets and transients (Model T integration): every preset of the engine registry
 * (engines/index.ts — cruise and full throttle in high gear, hill in low, idle in neutral, Ford's fixed-speed
 * dyno point, cruise on the battery) runs in EngineSimulator without exceptions, stalls, misfires or
 * non-finite output; the engine keeps running through operating-point changes mid-run (gear, grade, spark
 * lever, hand throttle, MAG/BAT, fuel, speed mode) and a reset; in free speed with the vehicle load the
 * engine summary closes the crank energy balance (brake torque — the engine's output, ΔE_kin of the crank and
 * mechanisms only — = cycle-mean load torque + the car's kinetic-energy change / 4π, loadInertiaTorque); plus
 * the CPU cost per engine cycle of EngineSimulator at 0.5° snapshots
 * (logged). The model is uncalibrated: no Ford numbers are asserted here. The ≥ 10-cycle runs of every
 * preset are part of the integration report (scratch runs); here each preset runs one engine cycle after
 * a one-cycle warm-up. Runtime ≈ 60–120 s (vitest).
 */
import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../../src/physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../../src/physics/core/snapshot';
import { ENGINES } from '../../src/physics/engines/index';
import { MODEL_T, MODEL_T_CRUISE, MODEL_T_FORD_DYNO } from '../../src/physics/engines/model-t';
import { CycleModel } from '../../src/physics/cycle/cycle-model';
import { EngineSimulator } from '../../src/physics/cycle/index';

const def = ENGINES['ford-model-t'];
const MPS_PER_MPH = 0.44704;

/** Non-finite numbers of a snapshot / summary tree, except the documented NaN fields. */
function nonFinite(o: unknown, path = ''): string[] {
  const out: string[] = [];
  const walk = (v: unknown, p: string): void => {
    if (typeof v === 'number') {
      if (!Number.isFinite(v) && !/(firstSparkDeg|sparkDeg|ca10|ca50|ca90|knockOnsetDeg|isfc)$/.test(p)) out.push(`${p}=${v}`);
    } else if (Array.isArray(v)) v.forEach((x, i) => walk(x, `${p}[${i}]`));
    else if (v && typeof v === 'object') for (const [k, x] of Object.entries(v)) walk(x, `${p}.${k}`);
  };
  walk(o, path);
  return out;
}

/** Advance n snapshots, checking every snapshot's speed and every 4th snapshot's numbers. */
function advance(sim: EngineSimulator, n: number, bad: string[], label: string): EngineSnapshot {
  let s = sim.advanceToNextSnapshot();
  for (let i = 0; i < n; i++) {
    s = sim.advanceToNextSnapshot();
    if (!(s.rpm > 100)) bad.push(`${label}: rpm ${s.rpm} at t ${s.t}`);
    if (i % 4 === 0 && bad.length < 10) bad.push(...nonFinite(s, label));
  }
  return s;
}

describe('every Model T preset runs in EngineSimulator', () => {
  for (const p of def.presets) {
    it(`${p.id}: no exception, stall, misfire or non-finite output`, () => {
      const op: OperatingPoint = { ...def.defaultOperatingPoint, ...p.op };
      const sim = new EngineSimulator(def.spec, op, { snapshotEveryDeg: 2, warmupCycles: 1 });
      const bad: string[] = [];
      let s = sim.advanceToNextSnapshot();
      let n = 0;
      while (s.cycle < 1 && n < 200_000) {
        s = sim.advanceToNextSnapshot();
        if (!(s.rpm > 100)) bad.push(`${p.id}: rpm ${s.rpm}`);
        if (n % 4 === 0 && bad.length < 10) bad.push(...nonFinite(s, p.id));
        n++;
      }
      expect(s.cycle).toBe(1);
      const sums: CycleSummary[] = sim.drainCycleSummaries();
      expect(sums.length).toBeGreaterThanOrEqual(4);
      for (const x of sums) {
        bad.push(...nonFinite(x, `${p.id} cyl${x.cylinder}`));
        expect(x.misfire, `${p.id} cylinder ${x.cylinder}`).toBe(false);
        expect(x.sparkCount!).toBeGreaterThan(0);
      }
      expect(bad).toEqual([]);
      const e = sums.find((x) => x.engine)!.engine!;
      if (op.speedMode === 'free' && op.load?.kind === 'vehicle' && op.load.gear !== 'neutral') {
        // road speed from the gearing: v = r ω/(G G_f) (≈ 41 rpm per mph in high)
        const G = MODEL_T.vehicle.gears[op.load.gear] * MODEL_T.vehicle.finalDrive;
        const v = (MODEL_T.vehicle.wheelRadius * (s.rpm * 2 * Math.PI)) / 60 / G;
        expect(s.vehicleSpeed!).toBeCloseTo(v, 9);
        expect(e.vehicleSpeed!).toBeGreaterThan(0);
      }
      if (op.speedMode === 'free' && op.load?.kind === 'vehicle' && op.load.gear === 'neutral') expect(s.loadTorque).toBe(0);
      console.log(
        `[Model T preset] ${p.id}: ${s.rpm.toFixed(0)} rpm, T_b ${e.brakeTorque.toFixed(1)} N m, load ${e.loadTorque.toFixed(1)} N m, ` +
          `${((e.vehicleSpeed ?? 0) / MPS_PER_MPH).toFixed(1)} mph, η_v ${e.volumetricEfficiency.toFixed(3)}, first spark ${sums[0].sparkDeg!.toFixed(1)}° (${sums[0].sparkCount} breakdowns)`,
      );
    });
  }
});

describe('Model T: operating-point changes mid-run and reset', () => {
  it('gear, grade, lever, throttle, MAG/BAT, fuel and speed-mode changes keep the engine running', () => {
    const sim = new EngineSimulator(MODEL_T, MODEL_T_CRUISE, { snapshotEveryDeg: 4, warmupCycles: 1 });
    const bad: string[] = [];
    const steps: [string, Partial<OperatingPoint>][] = [
      ['low gear', { load: { kind: 'vehicle', gear: 'low', grade: 0 } }],
      ['5 % grade', { load: { kind: 'vehicle', gear: 'low', grade: 0.05 } }],
      ['lever 40', { sparkAdvanceDeg: 40 }],
      ['throttle 0.8', { throttle: 0.8 }],
      ['battery', { ignitionSource: 'battery' }],
      ['PRF 70, φ 1.0', { fuel: { kind: 'PRF', octaneNumber: 70 }, equivalenceRatio: 1.0 }],
      ['neutral', { load: { kind: 'vehicle', gear: 'neutral', grade: 0 }, throttle: 0.1 }],
      ['high gear, magneto', { load: { kind: 'vehicle', gear: 'high', grade: 0 }, ignitionSource: 'magneto', throttle: 0.5 }],
      ['fixed 1000 rpm', { speedMode: 'fixed', rpm: 1000 }],
      ['free again', { speedMode: 'free' }],
      ['brake load', { load: { kind: 'brake', refRpm: 1000, exponent: 2 }, loadTorque: 60 }],
    ];
    advance(sim, 120, bad, 'start');
    for (const [label, patch] of steps) {
      sim.setOperatingPoint(patch);
      advance(sim, 120, bad, label);
    }
    expect(sim.model.op.load).toEqual({ kind: 'brake', refRpm: 1000, exponent: 2 });
    expect(sim.model.load.kind).toBe('brake');
    // the dash switch applies at each coil's next timer make; the engine now runs on the magneto again
    for (const c of sim.model.cylinders) expect(c.ign!.ignitionSource).toBe('magneto');
    sim.reset();
    expect(sim.time).toBe(0);
    advance(sim, 200, bad, 'after reset');
    for (const x of sim.drainCycleSummaries()) bad.push(...nonFinite(x, `summary${x.cycle}`));
    expect(bad).toEqual([]);
  });
});

describe('Model T free speed with the vehicle load', () => {
  it('the engine summary closes the crank energy balance: brake torque = cycle-mean load torque + the car\'s ΔE_kin/4π', () => {
    const m = new CycleModel(MODEL_T, MODEL_T_CRUISE, { warmupCycles: 1 });
    const es = m.runCycles(2).filter((s) => s.engine).map((s) => s.engine!);
    expect(es.length).toBe(2);
    for (const e of es) {
      // brake = indicated − friction − ΔE_kin/4π with ΔE_kin of the ENGINE (crank + flywheel + mechanisms): its
      // output goes into the road load and into the car's kinetic energy (reflected through high gear,
      // loadInertiaTorque); residual: the friction mean at the cycle-mean speed of a changing speed. (Round 1 put
      // the car's inertia into ΔE_kin and asserted brake = load: brake torque was the road load, not the engine's.)
      expect(Math.abs(e.brakeTorque - (e.loadTorque + e.loadInertiaTorque!))).toBeLessThan(5e-3 * e.frictionTorque);
      expect(e.clutchLoss).toBe(0);
      expect(e.vehicleSpeed! / MPS_PER_MPH).toBeGreaterThan(5);
      expect(e.vehicleSpeed! / MPS_PER_MPH).toBeLessThan(60);
    }
    expect(m.load.inertia()).toBeGreaterThan(5 * MODEL_T.masses.rotatingInertia);
  });
});

describe('Model T EngineSimulator CPU cost', () => {
  it('ms per engine cycle at 0.5° snapshots, 600 and 1600 rpm (fixed speed, wide open; logged)', () => {
    const out: string[] = [];
    for (const rpm of [600, 1600]) {
      const sim = new EngineSimulator(MODEL_T, { ...MODEL_T_FORD_DYNO, rpm }, { snapshotEveryDeg: 0.5, warmupCycles: 1 });
      let s = sim.advanceToNextSnapshot();
      while (s.cycle < 1) s = sim.advanceToNextSnapshot();
      const t0 = performance.now();
      let n = 0;
      while (s.cycle < 2) {
        s = sim.advanceToNextSnapshot();
        n++;
      }
      const ms = performance.now() - t0;
      out.push(`${rpm} rpm: ${ms.toFixed(0)} ms per engine cycle (real time ${(120000 / rpm).toFixed(0)} ms), ${n} snapshots`);
      expect(ms).toBeLessThan(30_000); // generous for loaded machines (vitest ≈ 2× slower than node)
    }
    console.log(`[Model T] EngineSimulator, 0.5° snapshots: ${out.join('; ')}`);
  });
});
