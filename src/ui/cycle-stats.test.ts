import { describe, expect, it } from 'vitest';
import { computeCycleStats, CYCLE_METRICS, cyclesToCsv, formatMetric, mean, metricMean, stdev } from './cycle-stats';
import { makeCycle } from './test-helpers';

const def = (key: string) => CYCLE_METRICS.find((d) => d.key === key)!;

describe('statistics', () => {
  it('mean / sample stdev ignore non-finite values', () => {
    expect(mean([1, 2, 3, NaN])).toBe(2);
    expect(stdev([2, 4, 4, 4, 5, 5, 7, 9])).toBeCloseTo(2.13809, 5);
    expect(stdev([1])).toBeNaN();
    expect(mean([])).toBeNaN();
  });

  it('COV and LNV of IMEP, knock share, misfires', () => {
    const imeps = [9, 10, 11, 10];
    const cycles = imeps.map((b, i) =>
      makeCycle({ cycle: i, imepNet: b * 1e5, knockOnsetDeg: i === 2 ? 12 : NaN, misfire: i === 3 }),
    );
    const st = computeCycleStats(cycles);
    expect(st.n).toBe(4);
    expect(st.imepMean).toBeCloseTo(10, 12);
    expect(st.imepStd).toBeCloseTo(Math.sqrt(2 / 3), 12);
    expect(st.imepCov).toBeCloseTo(100 * Math.sqrt(2 / 3) / 10, 10);
    expect(st.imepLnv).toBeCloseTo(90, 10);
    expect(st.knockFraction).toBe(0.25);
    expect(st.misfires).toBe(1);
    expect(st.etaMean).toBeCloseTo(32, 10);
  });

  it('handles an empty window', () => {
    const st = computeCycleStats([]);
    expect(st.n).toBe(0);
    expect(st.imepMean).toBeNaN();
    expect(st.knockFraction).toBeNaN();
  });
});

describe('metric definitions', () => {
  it('convert to engineering units', () => {
    const c = makeCycle({ imepNet: 9.5e5, isfc: 0.25 / 3.6e6, indicatedEfficiency: 0.335, coFraction: 0.0123, maxPressureRiseRate: 2.5e5 });
    expect(def('imepNet').value(c)).toBeCloseTo(9.5, 12);
    expect(def('isfc').value(c)).toBeCloseTo(250, 9);
    expect(def('etaI').value(c)).toBeCloseTo(33.5, 12);
    expect(def('co').value(c)).toBeCloseTo(1.23, 12);
    expect(def('dpdt').value(c)).toBeCloseTo(2.5, 12);
    expect(def('burn').value(c)).toBe(23);
    expect(formatMetric(def('imepNet'), c)).toBe('9.50');
    expect(formatMetric(def('knockDeg'), c)).toBe('—');
    expect(formatMetric(def('misfire'), makeCycle({ misfire: true }))).toBe('YES');
  });

  it('flags knocking and misfiring cycles', () => {
    expect(def('mapo').flag!(makeCycle({ knockOnsetDeg: 10 }))).toBe('knock');
    expect(def('mapo').flag!(makeCycle())).toBeNull();
    expect(def('misfire').flag!(makeCycle({ misfire: true }))).toBe('misfire');
  });

  it('covers every quantity the results table must show', () => {
    const keys = CYCLE_METRICS.map((d) => d.key);
    for (const k of ['imepNet', 'imepGross', 'pmep', 'etaI', 'isfc', 'pmax', 'pmaxDeg', 'ca10', 'ca50', 'ca90', 'burn', 'dpdt', 'knockDeg', 'mapo', 'no', 'co', 'xr', 'etaV', 'misfire'])
      expect(keys).toContain(k);
    expect(new Set(keys).size).toBe(keys.length);
  });

  it('metricMean skips NaN cycles', () => {
    const cs = [makeCycle({ knockOnsetDeg: 10 }), makeCycle({ knockOnsetDeg: NaN }), makeCycle({ knockOnsetDeg: 20 })];
    expect(metricMean(def('knockDeg'), cs)).toBe(15);
  });
});

describe('CSV export', () => {
  it('has a header with units and one row per cycle', () => {
    const csv = cyclesToCsv([makeCycle({ cycle: 7 }), makeCycle({ cycle: 8, ca50: NaN })]);
    const lines = csv.trim().split('\n');
    expect(lines).toHaveLength(3);
    expect(lines[0].startsWith('cycle,IMEP net [bar],')).toBe(true);
    expect(lines[0].split(',').length).toBe(CYCLE_METRICS.length + 1);
    expect(lines[1].startsWith('7,9.5,')).toBe(true);
    const ca50Col = 1 + CYCLE_METRICS.findIndex((d) => d.key === 'ca50');
    expect(lines[2].split(',')[ca50Col]).toBe('');
  });
});
