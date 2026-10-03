/**
 * End-gas ignition-delay models τ(T, p, φ, fuel, x_res) for the Livengood–Wu knock integral.
 *
 *  - `douaudEyzat`: the classical PRF correlation of Douaud & Eyzat (1978), fitted to CFR
 *    engine knock data;
 *  - `TabulatedIgnitionDelay`: detailed chemistry — ln τ tabulated from constant-volume 0-D
 *    reactor simulations with real PRF kinetic mechanisms (see ignition-delay-llnl.ts for the
 *    LLNL tables; default prfDetailedChemistry), multilinear in (1000/T, ln p, φ, ON, x_res).
 *
 * Conventions: T unburned (end-gas) temperature, K; p pressure, Pa; φ fuel/air equivalence
 * ratio of the fresh charge; xResidual = MOLE fraction of residual (complete-combustion
 * products of the same fresh charge) in the unburned gas (see residualMoleFraction for the
 * conversion from a mass fraction); τ in s.
 */
import { P_ATM } from '../core/constants';
import { NS, SP } from '../core/species';
import { type FuelBlend, ISOOCTANE_DENSITY_60F, NHEPTANE_DENSITY_60F } from '../thermo/fuels';
import { MOLAR_MASS } from '../thermo/thermo';

/** A model of the (main, hot) autoignition delay of the unburned mixture. */
export interface IgnitionDelayModel {
  readonly id: string;
  /**
   * Ignition delay, s, at unburned-gas temperature T (K), pressure p (Pa), fresh-charge
   * equivalence ratio phi, fuel, and residual mole fraction xResidual (0..1).
   */
  tau(T: number, p: number, phi: number, fuel: FuelBlend, xResidual: number): number;
}

// ---------------------------------------------------------------------------------------
// PRF octane number of a fuel
// ---------------------------------------------------------------------------------------

/** Liquid molar densities (mol/m³) of iso-octane / n-heptane at 60 °F (thermo/fuels.ts data). */
const N_ISO = ISOOCTANE_DENSITY_60F / MOLAR_MASS[SP.IC8H18];
const N_HEP = NHEPTANE_DENSITY_60F / MOLAR_MASS[SP.NC7H16];

/**
 * PRF octane number (liquid-volume % iso-octane, 0..100) of a fuel blend: fuel.octaneNumber
 * when set, else inverted from its vapour IC8H18/NC7H16 mole fractions (inverse of
 * thermo/fuels.ts prfIsooctaneMoleFraction). NaN if the fuel contains other species.
 */
export function prfOctaneNumber(fuel: FuelBlend): number {
  if (fuel.octaneNumber !== undefined) return fuel.octaneNumber;
  const X = fuel.X;
  let other = 0;
  for (let k = 0; k < NS; k++) if (k !== SP.IC8H18 && k !== SP.NC7H16) other += X[k];
  const xi = X[SP.IC8H18];
  const xh = X[SP.NC7H16];
  if (other > 1e-12 * (xi + xh) || !(xi + xh > 0)) return NaN;
  // x_iso = v n_iso / (v n_iso + (1 − v) n_hep)  ⇒  v = x n_hep / (n_iso (1 − x) + x n_hep)
  const x = xi / (xi + xh);
  return (100 * x * N_HEP) / (N_ISO * (1 - x) + x * N_HEP);
}

/**
 * Convert a residual MASS fraction to the mole fraction used by the models:
 * x = (y/M_res) / (y/M_res + (1 − y)/M_fresh). Molar masses in kg/mol.
 */
export function residualMoleFraction(massFraction: number, molarMassResidual: number, molarMassFresh: number): number {
  const a = massFraction / molarMassResidual;
  return a / (a + (1 - massFraction) / molarMassFresh);
}

// ---------------------------------------------------------------------------------------
// Douaud & Eyzat (1978)
// ---------------------------------------------------------------------------------------

/**
 * Douaud & Eyzat (1978), "Four-octane-number method for predicting the anti-knock behavior
 * of fuels and engines", SAE 780080:  τ = 17.68 ms (ON/100)^3.402 (p/atm)^−1.7 exp(3800 K/T).
 * Form, constants and units (τ in ms, p in atm, T in K) as quoted by Heywood (1988) ch. 9 and
 * confirmed in fetched secondary sources: ANSYS Fluent 12 UDF Manual §2.3.10 example
 * (τ = 0.01768 (ON/100)^3.402 p^−1.7 exp(3800/T), "pressure in atm", τ in s), N. Sharma, KTH
 * thesis (2018) eq. (3), and the Leeds PhD thesis (etheses.whiterose.ac.uk/7554) eq. (1.6).
 * (Some papers print 17.69.) The original SAE paper itself was not accessible.
 * Independent of φ and residual. Fitted to CFR knock of PRFs around RON/MON conditions;
 * τ → 0 as ON → 0 (n-heptane), so it is not meaningful for low-octane PRFs.
 */
export function douaudEyzatTau(T: number, p: number, octaneNumber: number): number {
  return 17.68e-3 * Math.pow(octaneNumber / 100, 3.402) * Math.pow(p / P_ATM, -1.7) * Math.exp(3800 / T);
}

/** Douaud–Eyzat as an IgnitionDelayModel (ON from prfOctaneNumber; throws for non-PRF fuels). */
export const douaudEyzat: IgnitionDelayModel = Object.freeze({
  id: 'douaud-eyzat-1978',
  tau(T: number, p: number, _phi: number, fuel: FuelBlend, _xResidual: number): number {
    const on = prfOctaneNumber(fuel);
    if (Number.isNaN(on)) throw new RangeError(`Douaud–Eyzat needs a PRF octane number (fuel '${fuel.label}')`);
    return douaudEyzatTau(T, p, on);
  },
});

// ---------------------------------------------------------------------------------------
// Tabulated detailed chemistry
// ---------------------------------------------------------------------------------------

/** Uniform axis: node k at start + k·step, n nodes. */
export interface UniformAxis {
  start: number;
  step: number;
  n: number;
  unit?: string;
}

/** Axis given by its (strictly ascending) node values. */
export interface NodeAxis {
  nodes: number[];
  unit?: string;
}

/** Node values of a uniform or explicit axis. */
export function axisNodes(a: UniformAxis | NodeAxis): number[] {
  if ('nodes' in a) return a.nodes.slice();
  return Array.from({ length: a.n }, (_, k) => a.start + k * a.step);
}

/**
 * JSON layout written by tools/reference/chemistry_ignition_table.py. Flat index
 * ((((iT·nP + iP)·nPhi + iPhi)·nOn + iOn)·nRes + iRes); lnTau = ln(τ/s).
 */
export interface IgnitionDelayTable {
  description: string;
  inputHash?: string;
  reactor: string;
  tauMax: number;
  axes: {
    /** 1000/T, 1/K (ascending: uniform, or explicit nodes). */
    invT: UniformAxis | NodeAxis;
    /** ln(p/Pa) (ascending: uniform, or explicit nodes). */
    lnP: UniformAxis | NodeAxis;
    phi: number[];
    on: number[];
    xres: number[];
  };
  order: string[];
  lnTau: number[];
  lnTau1: number[];
  capped?: number[];
}

/**
 * Index/weight on an ascending axis WITHOUT clamping: outside the axis the edge interval is
 * returned with a weight < 0 or > 1 (for slope-limited linear extrapolation). NaN x → weight NaN.
 */
function locateExtrapolating(nodes: Float64Array, x: number, out: Float64Array, k: number, idx: Int32Array): void {
  const n = nodes.length;
  let i = 0;
  if (x >= nodes[n - 1]) i = n - 2;
  else while (i < n - 2 && x >= nodes[i + 1]) i++;
  idx[k] = i;
  out[k] = (x - nodes[i]) / (nodes[i + 1] - nodes[i]);
}

/** Index/weight on a non-uniform ascending axis with clamping to the end nodes. */
function locateClamped(nodes: Float64Array, x: number, out: Float64Array, k: number, idx: Int32Array): void {
  const n = nodes.length;
  if (!(x > nodes[0])) {
    idx[k] = 0;
    out[k] = 0;
    return;
  }
  if (!(x < nodes[n - 1])) {
    idx[k] = n - 2;
    out[k] = 1;
    return;
  }
  let i = 0;
  while (x >= nodes[i + 1]) i++;
  idx[k] = i;
  out[k] = (x - nodes[i]) / (nodes[i + 1] - nodes[i]);
}

/**
 * Multilinear interpolation of ln τ in (1000/T, ln p, φ, ON, x_res) over a detailed-chemistry
 * table (1000/T and ln p axes uniform or given as explicit ascending nodes). Outside the table:
 * φ and ON are clamped to the table range (x_res below 0 too); 1000/T and ln p are
 * extrapolated linearly from the edge interval with the edge slope limited to
 * ∂lnτ/∂(1000/T) ≥ 0 (τ grows as T falls below / shrinks as T rises above the table) and
 * ∂lnτ/∂ln p ≤ 0; x_res above its last node is extrapolated linearly (no limiter) up to twice
 * the last node and clamped beyond. Allocation-free per call.
 * Extrapolation is only a guard: the table must COVER the states the end gas actually visits.
 * Livengood–Wu integrates from IVC, so that includes the whole compression stroke
 * (≈ 350–650 K at 1–10 bar before the spark), where a 650 K / 10 bar-bounded table
 * extrapolated to τ 2–30× too short (see ignition-delay-llnl.ts).
 */
export class TabulatedIgnitionDelay implements IgnitionDelayModel {
  readonly id: string;
  readonly table: IgnitionDelayTable;
  private readonly lnTau: Float64Array;
  private readonly lnTau1: Float64Array;
  private readonly invT: Float64Array;
  private readonly lnP: Float64Array;
  private readonly phi: Float64Array;
  private readonly on: Float64Array;
  private readonly xres: Float64Array;
  private readonly sT: number;
  private readonly sP: number;
  private readonly sPhi: number;
  private readonly sOn: number;
  // scratch
  private readonly w = new Float64Array(5);
  private readonly ix = new Int32Array(5);

  constructor(id: string, table: IgnitionDelayTable) {
    this.id = id;
    this.table = table;
    const a = table.axes;
    this.invT = Float64Array.from(axisNodes(a.invT));
    this.lnP = Float64Array.from(axisNodes(a.lnP));
    for (const ax of [this.invT, this.lnP]) {
      for (let i = 1; i < ax.length; i++) if (!(ax[i] > ax[i - 1])) throw new RangeError('table axis not ascending');
      if (ax.length < 2) throw new RangeError('table axis needs ≥ 2 nodes');
    }
    const n = this.invT.length * this.lnP.length * a.phi.length * a.on.length * a.xres.length;
    if (table.lnTau.length !== n || table.lnTau1.length !== n) {
      throw new RangeError(`ignition table size ${table.lnTau.length} ≠ ${n}`);
    }
    if (table.order.join(',') !== 'invT,lnP,phi,on,xres') throw new RangeError('unexpected table order');
    this.lnTau = Float64Array.from(table.lnTau);
    this.lnTau1 = Float64Array.from(table.lnTau1);
    this.phi = Float64Array.from(a.phi);
    this.on = Float64Array.from(a.on);
    this.xres = Float64Array.from(a.xres);
    this.sPhi = a.on.length * a.xres.length;
    this.sOn = a.xres.length;
    this.sP = a.phi.length * this.sPhi;
    this.sT = this.lnP.length * this.sP;
  }

  /** ln(τ/s) of the main (hot) ignition. on = PRF octane number. */
  lnTauAt(T: number, p: number, phi: number, on: number, xRes: number): number {
    return this.interp(this.lnTau, T, p, phi, on, xRes);
  }

  /** ln(τ1/s) of the first (cool-flame) stage (= ln τ where ignition is single-stage). */
  lnTauFirstStageAt(T: number, p: number, phi: number, on: number, xRes: number): number {
    return this.interp(this.lnTau1, T, p, phi, on, xRes);
  }

  /** Main ignition delay, s, for PRF octane number `on`. */
  tauPRF(T: number, p: number, phi: number, on: number, xRes: number): number {
    return Math.exp(this.interp(this.lnTau, T, p, phi, on, xRes));
  }

  /** IgnitionDelayModel: main ignition delay, s (NaN for non-PRF fuels). */
  tau(T: number, p: number, phi: number, fuel: FuelBlend, xResidual: number): number {
    return Math.exp(this.interp(this.lnTau, T, p, phi, prfOctaneNumber(fuel), xResidual));
  }

  /** First-stage delay τ1, s (= τ where single-stage). */
  tauFirstStage(T: number, p: number, phi: number, fuel: FuelBlend, xResidual: number): number {
    return Math.exp(this.interp(this.lnTau1, T, p, phi, prfOctaneNumber(fuel), xResidual));
  }

  private interp(f: Float64Array, T: number, p: number, phi: number, on: number, xRes: number): number {
    if (Number.isNaN(on)) return NaN;
    const w = this.w;
    const ix = this.ix;
    // 1000/T and ln p: edge intervals outside the table, weights may leave [0, 1]
    locateExtrapolating(this.invT, 1000 / T, w, 0, ix);
    locateExtrapolating(this.lnP, Math.log(p), w, 1, ix);
    const iT = ix[0];
    const wT = w[0];
    const iP = ix[1];
    const wP = w[1];
    locateClamped(this.phi, phi, w, 2, ix);
    locateClamped(this.on, on, w, 3, ix);
    locateClamped(this.xres, xRes, w, 4, ix);
    // x_res above the table: linear extrapolation in ln τ up to twice the last node (the
    // chemical effect of residual need not be monotone — H2O/CO third-body and initiation
    // effects can shorten τ at low T — so no slope limiter), clamped beyond that.
    const xr = this.xres;
    const nx = xr.length;
    if (xRes > xr[nx - 1]) {
      const x = Math.min(xRes, 2 * xr[nx - 1]);
      w[4] = (x - xr[nx - 2]) / (xr[nx - 1] - xr[nx - 2]);
    }
    const base = ix[2] * this.sPhi + ix[3] * this.sOn + ix[4];
    const w2 = w[2];
    const w3 = w[3];
    const w4 = w[4];
    const sPhi = this.sPhi;
    const sOn = this.sOn;
    // value at the 4 (T, p) corners after trilinear interpolation in (φ, ON, x_res)
    const c00 = iT * this.sT + iP * this.sP + base;
    const v00 = tri(f, c00, sPhi, sOn, w2, w3, w4);
    const v01 = tri(f, c00 + this.sP, sPhi, sOn, w2, w3, w4);
    const v10 = tri(f, c00 + this.sT, sPhi, sOn, w2, w3, w4);
    const v11 = tri(f, c00 + this.sT + this.sP, sPhi, sOn, w2, w3, w4);
    // p direction (ln τ must not increase with p beyond the table)
    const a0 = edge(v00, v01, wP, -1);
    const a1 = edge(v10, v11, wP, -1);
    // T direction (ln τ must not decrease with 1000/T beyond the table)
    return edge(a0, a1, wT, +1);
  }
}

/** Trilinear interpolation of f at the (φ, ON, x_res) cell starting at c. */
function tri(f: Float64Array, c: number, sPhi: number, sOn: number, w2: number, w3: number, w4: number): number {
  const f000 = f[c];
  const f001 = f[c + 1];
  const f010 = f[c + sOn];
  const f011 = f[c + sOn + 1];
  const f100 = f[c + sPhi];
  const f101 = f[c + sPhi + 1];
  const f110 = f[c + sPhi + sOn];
  const f111 = f[c + sPhi + sOn + 1];
  // w4 > 1: linear extrapolation in x_res beyond the last node (see interp)
  const g00 = f000 + w4 * (f001 - f000);
  const g01 = f010 + w4 * (f011 - f010);
  const g10 = f100 + w4 * (f101 - f100);
  const g11 = f110 + w4 * (f111 - f110);
  const h0 = g00 + w3 * (g01 - g00);
  const h1 = g10 + w3 * (g11 - g10);
  return h0 + w2 * (h1 - h0);
}

/**
 * Linear interpolation v0 + w (v1 − v0) inside [0, 1]; outside, extrapolation with the edge
 * slope limited to sign `dir` (+1: slope ≥ 0, −1: slope ≤ 0).
 */
function edge(v0: number, v1: number, w: number, dir: number): number {
  const s = v1 - v0;
  if (w >= 0 && w <= 1) return v0 + w * s;
  const sl = dir > 0 ? Math.max(s, 0) : Math.min(s, 0);
  return w < 0 ? v0 + w * sl : v1 + (w - 1) * sl;
}

/**
 * Temperature sensitivity ∂lnτ/∂T (1/K) of any model by central differences (±dT, K).
 * Not a hot path (used once at knock onset for the end-gas burn time).
 */
export function ignitionDelayTemperatureSensitivity(
  model: IgnitionDelayModel,
  T: number,
  p: number,
  phi: number,
  fuel: FuelBlend,
  xResidual: number,
  dT = 1,
): number {
  return (
    (Math.log(model.tau(T + dT, p, phi, fuel, xResidual)) - Math.log(model.tau(T - dT, p, phi, fuel, xResidual))) /
    (2 * dT)
  );
}
