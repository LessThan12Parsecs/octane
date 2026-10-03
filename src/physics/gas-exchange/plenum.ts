/**
 * Zero-dimensional filling-and-emptying gas volume (intake / exhaust plenum, or any 0D node).
 *
 * Conserved state: species moles N (Float64Array(NS), mol) and internal energy U (J, absolute
 * NASA-convention energies incl. formation, so reacting mixtures need no special treatment).
 * The gas inside is at rest and perfectly mixed (static = stagnation state).
 *
 * Balance equations for a control volume V(t) with signed mass flows ṁ_i (> 0 INTO the volume)
 * that carry the upstream stagnation enthalpy h₀,i (J/kg) and mass fractions Y_i:
 *   dN_k/dt = Σ_i ṁ_i Y_i,k / M_k
 *   dU/dt   = Σ_i ṁ_i h₀,i + Q̇ − p dV/dt           (Q̇ > 0: heat INTO the gas)
 * For an outflow (ṁ_i < 0) the carried state is the volume's own: Y_i = Y, h₀,i = h(T)
 * (the filling-and-emptying formulation of Heywood 1988 ch. 14 — UNVERIFIED section number;
 * kinetic energy in the volume neglected).
 * Because every stream leaves one volume with exactly the state it enters the next with,
 * a network of plenums connected by restrictions conserves ΣN_k and ΣU to rounding.
 *
 * State recovery: T from U(N, T) = U (thermo.temperatureFromUMolar — linear in N, safeguarded
 * Newton warm-started from the previous T), then p = (ΣN) R_u T / V.
 *
 * Units SI. Methods marked hot allocate nothing.
 */
import { R_UNIVERSAL } from '../core/constants';
import { NS } from '../core/species';
import { MOLAR_MASS as MOLAR_MASS_IMPORT, mixStateMass, mixUMolar, temperatureFromUMolar, type MixState } from '../thermo';

// Module-local aliases: loops must not read imported bindings (live-binding getters under
// some module transforms, e.g. vitest's, cost a call per access).
const MOLAR_MASS = MOLAR_MASS_IMPORT;
const NSP = NS;

/** Thermodynamic state of a gas at rest (static = stagnation). All SI, mass-specific. */
export interface GasState {
  /** Temperature, K. */
  T: number;
  /** Pressure, Pa. */
  p: number;
  /** Density, kg/m³. */
  rho: number;
  /** Mole fractions (NS). */
  readonly X: Float64Array;
  /** Mass fractions (NS). */
  readonly Y: Float64Array;
  /** Mean molar mass, kg/mol. */
  molarMass: number;
  /** Specific gas constant, J/(kg K). */
  R: number;
  /** cp/cv (frozen). */
  gamma: number;
  /** cp, J/(kg K). */
  cp: number;
  /** Specific enthalpy (absolute), J/kg — the stagnation enthalpy of gas leaving at rest. */
  h: number;
  /** Specific internal energy (absolute), J/kg. */
  u: number;
}

/** Allocate a zeroed {@link GasState}. */
export function newGasState(): GasState {
  return {
    T: 0,
    p: 0,
    rho: 0,
    X: new Float64Array(NS),
    Y: new Float64Array(NS),
    molarMass: 0,
    R: 0,
    gamma: 0,
    cp: 0,
    h: 0,
    u: 0,
  };
}

const MIX: MixState = { cp: 0, cv: 0, h: 0, u: 0, gamma: 0, molarMass: 0 };

/** Fill the mass-specific properties of `out` from out.X, out.T, out.p (hot). */
function finishState(out: GasState): GasState {
  const X = out.X;
  const Y = out.Y;
  mixStateMass(X, out.T, MIX);
  let m = 0;
  for (let k = 0; k < NSP; k++) m += X[k] * MOLAR_MASS[k];
  const inv = 1 / m;
  for (let k = 0; k < NSP; k++) Y[k] = X[k] * MOLAR_MASS[k] * inv;
  out.molarMass = m;
  out.R = R_UNIVERSAL / m;
  out.gamma = MIX.gamma;
  out.cp = MIX.cp;
  out.h = MIX.h;
  out.u = MIX.u;
  out.rho = out.p / (out.R * out.T);
  return out;
}

/**
 * Gas state from temperature (K), pressure (Pa) and mole fractions X (normalised on copy).
 * Use for fixed boundaries (ambient, reservoir). Writes into and returns `out` (hot).
 */
export function gasStateFromTPX(T: number, p: number, X: Float64Array, out: GasState): GasState {
  let s = 0;
  for (let k = 0; k < NSP; k++) s += X[k];
  const inv = 1 / s;
  for (let k = 0; k < NSP; k++) out.X[k] = X[k] * inv;
  out.T = T;
  out.p = p;
  return finishState(out);
}

/**
 * Gas state of a volume V (m³) holding species moles N (mol) with internal energy U (J).
 * Tguess (K) warm-starts the temperature inversion. Writes into and returns `out` (hot).
 */
export function gasStateFromNU(N: Float64Array, U: number, V: number, Tguess: number, out: GasState): GasState {
  let n = 0;
  for (let k = 0; k < NSP; k++) n += N[k];
  const inv = 1 / n;
  for (let k = 0; k < NSP; k++) out.X[k] = N[k] * inv;
  // U(N, T) is linear in N: invert directly on the extensive quantities.
  const T = temperatureFromUMolar(N, U, Tguess > 0 ? Tguess : 1000);
  out.T = T;
  out.p = (n * R_UNIVERSAL * T) / V;
  return finishState(out);
}

/**
 * A 0D gas volume with its conserved state, derived state and derivative accumulators.
 *
 * Typical use inside an ODE right-hand side:
 *   plenum.N.set(y.subarray(...)); plenum.U = y[..]; plenum.updateState();
 *   plenum.beginRates();  plenum.addFlow(ṁ, upstreamState) …;  copy plenum.dNdt, plenum.dUdt.
 */
export class Plenum {
  /** Species moles, mol (conserved state). */
  readonly N = new Float64Array(NS);
  /** Internal energy (absolute), J (conserved state). */
  U = 0;
  /** Volume, m³. */
  volume: number;
  /** Derived state (valid after updateState / setTPX). */
  readonly state: GasState = newGasState();
  /** dN/dt accumulator, mol/s. */
  readonly dNdt = new Float64Array(NS);
  /** dU/dt accumulator, W. */
  dUdt = 0;

  /**
   * @param volume m³
   * @param T initial temperature, K; @param p initial pressure, Pa; @param X initial mole fractions
   */
  constructor(volume: number, T: number, p: number, X: Float64Array) {
    if (!(volume > 0)) throw new RangeError('Plenum: volume must be > 0');
    this.volume = volume;
    this.setTPX(T, p, X);
  }

  /** Reset the conserved state to (T, p, X) in the current volume. */
  setTPX(T: number, p: number, X: Float64Array): void {
    gasStateFromTPX(T, p, X, this.state);
    const n = (p * this.volume) / (R_UNIVERSAL * T);
    for (let k = 0; k < NSP; k++) this.N[k] = n * this.state.X[k];
    this.U = mixUMolar(this.N, T);
  }

  /** Recover T, p, composition and properties from (N, U, V) (hot; warm-started). */
  updateState(): GasState {
    return gasStateFromNU(this.N, this.U, this.volume, this.state.T, this.state);
  }

  /** Total mass, kg. */
  mass(): number {
    let m = 0;
    for (let k = 0; k < NSP; k++) m += this.N[k] * MOLAR_MASS[k];
    return m;
  }

  /**
   * Zero the accumulators and add the wall heat and boundary work.
   * @param heatRate Q̇ into the gas, W; @param dVdt dV/dt, m³/s (uses state.p)
   */
  beginRates(heatRate = 0, dVdt = 0): void {
    this.dNdt.fill(0);
    this.dUdt = heatRate - this.state.p * dVdt;
  }

  /**
   * Add a stream with an explicit carried state (hot).
   * @param mdot kg/s, > 0 into the volume
   * @param Y mass fractions of the carried gas (the upstream side's)
   * @param h0 stagnation enthalpy of the carried gas, J/kg (the upstream side's)
   */
  addStream(mdot: number, Y: Float64Array, h0: number): void {
    const d = this.dNdt;
    for (let k = 0; k < NSP; k++) d[k] += (mdot * Y[k]) / MOLAR_MASS[k];
    this.dUdt += mdot * h0;
  }

  /**
   * Add a bidirectional connection (hot): inflow (mdot > 0) carries `other`'s state,
   * outflow (mdot < 0) carries this volume's own state.
   * @param mdot kg/s, > 0 into this volume
   * @param other state of the volume/boundary on the other side of the restriction
   */
  addFlow(mdot: number, other: GasState): void {
    if (mdot > 0) this.addStream(mdot, other.Y, other.h);
    else if (mdot < 0) this.addStream(mdot, this.state.Y, this.state.h);
  }
}
