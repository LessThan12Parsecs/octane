/**
 * Cycle-results bookkeeping: per-cycle metric definitions (units, formatting,
 * flags), cycle-to-cycle statistics (COV of IMEP etc.) and CSV export. Pure.
 */
import type { CycleSummary } from '../physics/core/snapshot';
import { fmt, fractionToPct, isfcToGPerKWh, jToMJ, kgToMg, paToBar } from './units';

export type MetricFlag = 'knock' | 'misfire' | null;

export interface MetricDef {
  key: string;
  label: string;
  unit: string;
  /** Tooltip text. */
  hint: string;
  decimals: number;
  value: (c: CycleSummary) => number;
  /** Optional per-cell flag for highlighting. */
  flag?: (c: CycleSummary) => MetricFlag;
  /** Excluded from the "mean" column (e.g. booleans). */
  noMean?: boolean;
  /** Section heading shown above this row. */
  group?: string;
}

const knockFlag = (c: CycleSummary): MetricFlag => (Number.isFinite(c.knockOnsetDeg) ? 'knock' : null);

export const CYCLE_METRICS: readonly MetricDef[] = [
  { group: 'Work', key: 'imepNet', label: 'IMEP net', unit: 'bar', hint: 'Net indicated mean effective pressure (720°)', decimals: 2, value: (c) => paToBar(c.imepNet) },
  { key: 'imepGross', label: 'IMEP gross', unit: 'bar', hint: 'Gross IMEP (compression + expansion strokes)', decimals: 2, value: (c) => paToBar(c.imepGross) },
  { key: 'pmep', label: 'PMEP', unit: 'bar', hint: 'Pumping mean effective pressure', decimals: 3, value: (c) => paToBar(c.pmep) },
  { key: 'etaI', label: 'η indicated', unit: '%', hint: 'Net indicated thermal efficiency (fuel LHV basis)', decimals: 1, value: (c) => fractionToPct(c.indicatedEfficiency) },
  { key: 'isfc', label: 'ISFC', unit: 'g/kWh', hint: 'Net indicated specific fuel consumption', decimals: 0, value: (c) => isfcToGPerKWh(c.isfc) },
  { key: 'wGross', label: 'W gross', unit: 'J', hint: 'Gross indicated work per cycle', decimals: 1, value: (c) => c.indicatedWorkGross },
  { key: 'heatLoss', label: 'Q wall', unit: 'J', hint: 'Wall heat loss over the closed cycle', decimals: 1, value: (c) => c.heatLoss },
  { group: 'Pressure', key: 'pmax', label: 'p max', unit: 'bar', hint: 'Peak cylinder pressure', decimals: 1, value: (c) => paToBar(c.peakPressure) },
  { key: 'pmaxDeg', label: '∠ p max', unit: '°CA', hint: 'Crank angle of peak pressure (ATDC positive)', decimals: 1, value: (c) => c.peakPressureDeg },
  { key: 'dpdt', label: 'dp/dθ max', unit: 'bar/°', hint: 'Maximum rate of pressure rise', decimals: 2, value: (c) => paToBar(c.maxPressureRiseRate) },
  { group: 'Combustion', key: 'ca10', label: 'CA10', unit: '°CA', hint: '10 % mass fraction burned', decimals: 1, value: (c) => c.ca10 },
  { key: 'ca50', label: 'CA50', unit: '°CA', hint: '50 % mass fraction burned (combustion phasing)', decimals: 1, value: (c) => c.ca50 },
  { key: 'ca90', label: 'CA90', unit: '°CA', hint: '90 % mass fraction burned', decimals: 1, value: (c) => c.ca90 },
  { key: 'burn', label: 'Burn 10–90', unit: '°CA', hint: 'Main combustion duration CA90 − CA10', decimals: 1, value: (c) => c.ca90 - c.ca10 },
  { key: 'misfire', label: 'Misfire', unit: '', hint: 'Spark kernel failed to develop into a flame', decimals: 0, value: (c) => (c.misfire ? 1 : 0), flag: (c) => (c.misfire ? 'misfire' : null), noMean: true },
  { group: 'Knock', key: 'knockDeg', label: 'Knock onset', unit: '°CA', hint: 'End-gas autoignition angle (— = no knock)', decimals: 1, value: (c) => c.knockOnsetDeg, flag: knockFlag },
  { key: 'knockXu', label: 'End gas @ onset', unit: '%', hint: 'Unburned mass fraction at knock onset', decimals: 1, value: (c) => (Number.isFinite(c.knockOnsetDeg) ? fractionToPct(c.knockEndGasFraction) : NaN), flag: knockFlag },
  { key: 'mapo', label: 'MAPO', unit: 'bar', hint: 'Maximum amplitude of pressure oscillation', decimals: 2, value: (c) => paToBar(c.mapo), flag: knockFlag },
  { group: 'Emissions', key: 'no', label: 'NO', unit: 'ppm', hint: 'Exhaust NO at EVO (wet, mole basis)', decimals: 0, value: (c) => c.noPpm },
  { key: 'co', label: 'CO', unit: '%', hint: 'Exhaust CO at EVO (mole %)', decimals: 2, value: (c) => fractionToPct(c.coFraction) },
  { group: 'Breathing', key: 'xr', label: 'Residual', unit: '%', hint: 'Residual gas mass fraction at IVC', decimals: 1, value: (c) => fractionToPct(c.residualFraction) },
  { key: 'etaV', label: 'η volumetric', unit: '%', hint: 'Volumetric efficiency (ambient-referenced)', decimals: 1, value: (c) => fractionToPct(c.volumetricEfficiency) },
  { key: 'mTrap', label: 'Trapped mass', unit: 'mg', hint: 'Cylinder mass at IVC', decimals: 0, value: (c) => kgToMg(c.trappedMass) },
  { key: 'mFuel', label: 'Fuel mass', unit: 'mg', hint: 'Fuel mass burned this cycle', decimals: 2, value: (c) => kgToMg(c.fuelMass) },
];

/** Format one metric for one cycle ("—" for NaN). */
export function formatMetric(def: MetricDef, c: CycleSummary): string {
  if (def.key === 'misfire') return c.misfire ? 'YES' : 'no';
  return fmt(def.value(c), def.decimals);
}

export function mean(xs: ArrayLike<number>): number {
  let s = 0;
  let n = 0;
  for (let i = 0; i < xs.length; i++) {
    const v = xs[i];
    if (Number.isFinite(v)) {
      s += v;
      n++;
    }
  }
  return n ? s / n : NaN;
}

/** Sample standard deviation (n − 1) of the finite values. */
export function stdev(xs: ArrayLike<number>): number {
  const m = mean(xs);
  let s = 0;
  let n = 0;
  for (let i = 0; i < xs.length; i++) {
    const v = xs[i];
    if (Number.isFinite(v)) {
      s += (v - m) * (v - m);
      n++;
    }
  }
  return n > 1 ? Math.sqrt(s / (n - 1)) : NaN;
}

/** Mean of a metric over cycles (NaN values ignored). */
export function metricMean(def: MetricDef, cycles: readonly CycleSummary[]): number {
  let s = 0;
  let n = 0;
  for (const c of cycles) {
    const v = def.value(c);
    if (Number.isFinite(v)) {
      s += v;
      n++;
    }
  }
  return n ? s / n : NaN;
}

export interface CycleStats {
  n: number;
  /** bar */
  imepMean: number;
  imepStd: number;
  /** Coefficient of variation of net IMEP, % (the standard cyclic-variability metric). */
  imepCov: number;
  /** Lowest normalised value of IMEP, % (min / mean). */
  imepLnv: number;
  /** Fraction of cycles with knock, 0..1. */
  knockFraction: number;
  misfires: number;
  /** Mean CA50, °CA. */
  ca50Mean: number;
  /** Mean net indicated efficiency, %. */
  etaMean: number;
}

export function computeCycleStats(cycles: readonly CycleSummary[]): CycleStats {
  const n = cycles.length;
  const imep = new Float64Array(n);
  let knock = 0;
  let mis = 0;
  let min = Infinity;
  for (let i = 0; i < n; i++) {
    const c = cycles[i];
    imep[i] = paToBar(c.imepNet);
    if (Number.isFinite(imep[i]) && imep[i] < min) min = imep[i];
    if (Number.isFinite(c.knockOnsetDeg)) knock++;
    if (c.misfire) mis++;
  }
  const m = mean(imep);
  const sd = stdev(imep);
  const ca50 = cycles.map((c) => c.ca50);
  const eta = cycles.map((c) => fractionToPct(c.indicatedEfficiency));
  return {
    n,
    imepMean: m,
    imepStd: sd,
    imepCov: m !== 0 && Number.isFinite(sd) ? (100 * sd) / Math.abs(m) : NaN,
    imepLnv: m > 0 && Number.isFinite(min) ? (100 * min) / m : NaN,
    knockFraction: n ? knock / n : NaN,
    misfires: mis,
    ca50Mean: mean(ca50),
    etaMean: mean(eta),
  };
}

/** CSV of all metrics, one row per cycle (raw numbers in the displayed units). */
export function cyclesToCsv(cycles: readonly CycleSummary[]): string {
  const head = ['cycle', ...CYCLE_METRICS.map((d) => (d.unit ? `${d.label} [${d.unit}]` : d.label))];
  const esc = (s: string): string => (/[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s);
  const lines = [head.map(esc).join(',')];
  for (const c of cycles) {
    const row = [String(c.cycle)];
    for (const d of CYCLE_METRICS) {
      const v = d.value(c);
      row.push(Number.isFinite(v) ? String(Number(v.toPrecision(10))) : '');
    }
    lines.push(row.join(','));
  }
  return lines.join('\n') + '\n';
}

/** Energy in mJ with sensible decimals (spark energy readout). */
export const fmtMJ = (j: number): string => fmt(jToMJ(j), jToMJ(j) < 10 ? 1 : 0);
