/** Engine-level cycle metrics and period units. */
import { describe, expect, it } from 'vitest';
import type { EngineCycleSummary } from '../physics/core/snapshot';
import { MODEL_T_FORD_WOT_TABLE } from '../physics/engines/model-t';
import { CYCLE_METRICS, cyclesToCsv, ENGINE_METRICS, engineMetricMean } from './cycle-stats';
import { makeCycle } from './test-helpers';
import { MPS_PER_MPH, mpsToMph, NM_PER_LBFT, nmToLbft, W_PER_HP, wattsToHp } from './units';

function engine(over: Partial<EngineCycleSummary> = {}): EngineCycleSummary {
  return {
    rpmMean: 1000,
    indicatedTorque: 120,
    frictionTorque: 20,
    brakeTorque: 100,
    brakePower: (100 * 2 * Math.PI * 1000) / 60,
    bmep: 4.3e5,
    imepNet: 5.2e5,
    fmep: 0.9e5,
    airMassFlow: 0.02,
    fuelMassFlow: 0.0015,
    volumetricEfficiency: 0.7,
    bsfc: 0.0015 / ((100 * 2 * Math.PI * 1000) / 60),
    brakeEfficiency: 0.15,
    loadTorque: 100,
    vehicleSpeed: 15,
    ...over,
  };
}

const metric = (key: string) => ENGINE_METRICS.find((d) => d.key === key)!;

describe('period units (NIST SP 811 exact definitions)', () => {
  it('horsepower, pound-foot and mile per hour', () => {
    expect(W_PER_HP).toBeCloseTo(745.69987158227022, 9);
    expect(NM_PER_LBFT).toBeCloseTo(1.3558179483314004, 12);
    expect(MPS_PER_MPH).toBe(0.44704);
    expect(wattsToHp(745.69987158227022)).toBeCloseTo(1, 12);
    expect(nmToLbft(NM_PER_LBFT * 65)).toBeCloseTo(65, 12);
    expect(mpsToMph(0.44704 * 42)).toBeCloseTo(42, 12);
  });

  it("Ford's WOT table is self-consistent through these constants: hp = lb·ft × rpm / 5252", () => {
    // P = T ω: an oracle for the conversions. Ford rounded the torque and hp columns independently (low-speed
    // hp to 0.5 hp; up to 1.8 % apart at 1100 rpm); the 1400 rpm torque is illegible.
    for (const [rpm, lbft, hp] of MODEL_T_FORD_WOT_TABLE) {
      if (rpm === 1400) continue;
      const watts = lbft * NM_PER_LBFT * ((2 * Math.PI * rpm) / 60);
      expect(Math.abs(wattsToHp(watts) - hp), `${rpm} rpm`).toBeLessThanOrEqual(Math.max(0.25, 0.02 * hp));
    }
  });
});

describe('engine metrics', () => {
  it('convert an EngineCycleSummary to the table units', () => {
    const e = engine();
    expect(metric('rpm').value(e)).toBe(1000);
    expect(metric('brakeTorque').value(e)).toBe(100);
    expect(metric('brakeTorqueLbft').value(e)).toBeCloseTo(100 / NM_PER_LBFT, 12);
    expect(metric('brakePower').value(e)).toBeCloseTo(10.47197551, 6); // kW
    expect(metric('brakeHp').value(e)).toBeCloseTo(10471.97551 / W_PER_HP, 6);
    expect(metric('bmep').value(e)).toBeCloseTo(4.3, 12);
    expect(metric('fmep').value(e)).toBeCloseTo(0.9, 12);
    expect(metric('etaVEngine').value(e)).toBeCloseTo(70, 12);
    expect(metric('bsfc').value(e)).toBeCloseTo((0.0015 / 10471.97551) * 3.6e9, 3); // g/kWh
    expect(metric('etaB').value(e)).toBeCloseTo(15, 12);
    expect(metric('vehicleSpeed').value(e)).toBeCloseTo(15 / 0.44704, 12);
    expect(metric('vehicleSpeed').value(engine({ vehicleSpeed: undefined }))).toBeNaN();
    expect(new Set(ENGINE_METRICS.map((d) => d.key)).size).toBe(ENGINE_METRICS.length);
  });

  it('brake torque is described as the engine output, without the car’s inertia (EngineCycleSummary.brakeTorque)', () => {
    const hint = metric('brakeTorque').hint;
    expect(hint).toMatch(/indicated − friction/);
    expect(hint).toMatch(/engine output/);
    expect(hint).toMatch(/excludes the car’s inertia/);
  });

  it('means skip missing values', () => {
    const list = [engine({ brakeTorque: 90 }), engine({ brakeTorque: 110 }), engine({ vehicleSpeed: undefined })];
    expect(engineMetricMean(metric('brakeTorque'), list)).toBeCloseTo(100, 12);
    expect(engineMetricMean(metric('vehicleSpeed'), list)).toBeCloseTo(15 / 0.44704, 12);
    expect(engineMetricMean(metric('rpm'), [])).toBeNaN();
  });
});

describe('CSV export with cylinders and engine results', () => {
  it('adds a cylinder column and engine columns only when the data has them', () => {
    const csv = cyclesToCsv([
      makeCycle({ cycle: 4, cylinder: 0, engine: engine() }),
      makeCycle({ cycle: 4, cylinder: 1 }),
      makeCycle({ cycle: 4, cylinder: 3 }),
    ]);
    const [head, r1, r2, r3] = csv.trim().split('\n');
    const cols = head.split(',');
    expect(cols.slice(0, 2)).toEqual(['cycle', 'cylinder']);
    expect(cols.length).toBe(2 + CYCLE_METRICS.length + ENGINE_METRICS.length);
    expect(cols.at(-1)).toBe('engine Road speed [mph]');
    expect(r1.split(',')[1]).toBe('1');
    expect(r2.split(',')[1]).toBe('2');
    expect(r3.split(',')[1]).toBe('4');
    const bt = cols.indexOf('engine Brake torque [N·m]');
    expect(r1.split(',')[bt]).toBe('100');
    expect(r2.split(',')[bt]).toBe('');
  });

  it('single-cylinder data without engine results keeps the original columns', () => {
    const csv = cyclesToCsv([makeCycle({ cycle: 1 })]);
    expect(csv.split('\n')[0].split(',').length).toBe(CYCLE_METRICS.length + 1);
  });
});
