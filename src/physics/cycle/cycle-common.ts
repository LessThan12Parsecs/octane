/**
 * Shared definitions of the cycle model (cycle-model.ts = the engine, cylinder.ts = one cylinder):
 * the layout of the continuous state vector y, cylinder modes, event kinds, numerical constants,
 * the validation trace record and small allocation-free helpers. cycle-model.ts re-exports the
 * public ones (their historical home).
 *
 * ── State vector y (Float64Array) ─────────────────────────────────────────────────────────────
 * The single-cylinder layout of the CFR model is kept verbatim (indices I_*, length NY = 159): the
 * engine part (θ, ω, the two plenums, venturi/outlet ledgers, gross-inflow ledgers, port
 * conductance, intake EGR marker, venturi pressure derivative) interleaved with the block of
 * CYLINDER 0 (I_CN … I_CRB, the 42 entries listed in CYLINDER_BLOCK_OFFSETS). Cylinders i ≥ 1 of a
 * multi-cylinder engine append contiguous blocks of CYLINDER_BLOCK_SIZE entries after NY, in the
 * order of CylinderStateIndex (cylinderStateIndex).
 */
import { NE, NS, SP, EL, ELEMENT_COUNTS } from '../core/species';
import type { OperatingPoint } from '../core/operating-point';
import {
  douaudEyzat,
  douaudEyzatLLNL,
  prfLLNLGasoline2011,
  prfLLNLv2,
  type IgnitionDelayModel,
} from '../chemistry';
import { MOLAR_MASS } from '../thermo/thermo';
import { noSwapDeltas } from './closure';
import type { CycleModelOptions } from './options';

// =============================================================================================
// State layout
// =============================================================================================

export const I_TH = 0;
export const I_OM = 1;
export const I_CN = 2;
export const I_CU = I_CN + NS;
export const I_CBG = I_CU + 1;
export const I_UT = I_CBG + 1;
export const I_SU = I_UT + 1;
export const I_MU = I_SU + 1;
export const I_MB = I_MU + 1;
export const I_ME = I_MB + 1;
export const I_IN = I_ME + 1;
export const I_IU = I_IN + NS;
export const I_IBG = I_IU + 1;
export const I_EN = I_IBG + 1;
export const I_EU = I_EN + NS;
export const I_EBG = I_EU + 1;
export const I_TK = I_EBG + 1; // mean-flow kinetic energy K
export const I_TKE = I_TK + 1; // turbulent kinetic energy k
export const I_SW = I_TKE + 1;
export const I_W = I_SW + 1; // ∫ p dV
export const I_Q = I_W + 1; // ∫ net wall heat loss dt: cylinder walls − intake-port gain (W·s)
export const I_MIVI = I_Q + 1; // ∫ intake-valve flow into the cylinder
export const I_MIVO = I_MIVI + 1; // ∫ intake-valve backflow out of the cylinder
export const I_MEVO = I_MIVO + 1; // ∫ exhaust-valve flow out of the cylinder
export const I_MEVI = I_MEVO + 1; // ∫ exhaust-valve backflow into the cylinder
export const I_HV = I_MEVI + 1; // ∫ enthalpy flow through the venturi (+ into the intake plenum)
export const I_HO = I_HV + 1; // ∫ enthalpy flow through the outlet (+ out to ambient)
export const I_LV = I_HO + 1; // ∫ species moles through the venturi (+ into the plenum), NS
export const I_LO = I_LV + NS; // ∫ species moles through the outlet (+ out), NS
export const I_MK = I_LO + NS; // ∫ end-gas (autoignition) burn rate dt, kg
// gross INFLOW ledgers of the two plenums (warm-up acceleration, relaxPlenum): enthalpy (J),
// burned-gas scalar (kg), species moles (NS) of everything that entered the intake plenum (venturi
// inflow + intake-valve backflow) and the exhaust plenum (exhaust-valve outflow + outlet backflow)
export const I_GIH = I_MK + 1;
export const I_GIBG = I_GIH + 1;
export const I_GIN = I_GIBG + 1;
export const I_GEH = I_GIN + NS;
export const I_GEBG = I_GEH + 1;
export const I_GEN = I_GEBG + 1;
export const I_GIK = I_GEN + NS; // ∫ intake-port thermal conductance h·A dt (J/K), warm-up relaxation
// external-EGR mass scalars (kg) of the cylinder and the intake plenum: the part of the burned-gas
// scalar m_bg that entered as EGR through the carburettor (the products dilution of the TRAPPED
// charge is m_bg/m; the internal residual is (m_bg − m_egr)/m) — validation round 2
export const I_CEG = I_GIK + 1;
export const I_IEG = I_CEG + 1;
// ∫ −∂ṁ_venturi/∂p_plenum dt (kg/Pa): the intake plenum's pressure-mode relaxation rate (warm-up)
export const I_GIV = I_IEG + 1;
// ∫ gas-to-surface heat flow dt per surface (J): head, piston, liner, intake valves, exhaust valves
// (lumped wall-temperature model, updateWalls)
export const I_QS = I_GIV + 1;
// crevice zone (closed phase): unburned and burned (complete-combustion products) mass at the wall
// temperature and the cylinder pressure, kg
export const I_CRU = I_QS + 5;
export const I_CRB = I_CRU + 1;
/** Length of the single-cylinder state vector (the engine part + cylinder 0's block). */
export const NY = I_CRB + 1;

/** Absolute state indices of one cylinder's block. */
export interface CylinderStateIndex {
  /** Open-phase zone: species moles (NS), internal energy, burned-gas scalar. */
  CN: number;
  CU: number;
  CBG: number;
  /** Closed phase: U_tot, S_u, m_u, m_b, m_e. */
  UT: number;
  SU: number;
  MU: number;
  MB: number;
  ME: number;
  /** K–k turbulence: K, k, swirl. */
  TK: number;
  TKE: number;
  SW: number;
  /** Ledgers: ∫p dV, ∫(Q̇_wall − Q̇_port) dt, valve flows (in/out of both valves), end-gas burn. */
  W: number;
  Q: number;
  MIVI: number;
  MIVO: number;
  MEVO: number;
  MEVI: number;
  MK: number;
  /** External-EGR marker of the cylinder charge. */
  CEG: number;
  /** Per-surface heat ledgers (5). */
  QS: number;
  /** Crevice zone unburned / burned mass. */
  CRU: number;
  CRB: number;
}

/** Number of state entries per cylinder. */
export const CYLINDER_BLOCK_SIZE = NS + 25;

/** Cylinder 0: the legacy single-cylinder indices. */
export const CYLINDER0_STATE: Readonly<CylinderStateIndex> = Object.freeze({
  CN: I_CN, CU: I_CU, CBG: I_CBG, UT: I_UT, SU: I_SU, MU: I_MU, MB: I_MB, ME: I_ME, TK: I_TK, TKE: I_TKE, SW: I_SW,
  W: I_W, Q: I_Q, MIVI: I_MIVI, MIVO: I_MIVO, MEVO: I_MEVO, MEVI: I_MEVI, MK: I_MK, CEG: I_CEG, QS: I_QS, CRU: I_CRU, CRB: I_CRB,
});

/** State indices of cylinder `i` (0-based): the legacy indices for i = 0, an appended block otherwise. */
export function cylinderStateIndex(i: number): Readonly<CylinderStateIndex> {
  if (i === 0) return CYLINDER0_STATE;
  const b = NY + (i - 1) * CYLINDER_BLOCK_SIZE;
  const CN = b;
  const CU = CN + NS;
  return Object.freeze({
    CN, CU, CBG: CU + 1, UT: CU + 2, SU: CU + 3, MU: CU + 4, MB: CU + 5, ME: CU + 6, TK: CU + 7, TKE: CU + 8, SW: CU + 9,
    W: CU + 10, Q: CU + 11, MIVI: CU + 12, MIVO: CU + 13, MEVO: CU + 14, MEVI: CU + 15, MK: CU + 16, CEG: CU + 17,
    QS: CU + 18, CRU: CU + 23, CRB: CU + 24,
  });
}

/** Length of the state vector of an engine with n cylinders (without the block-surface heat ledgers). */
export function stateLength(n: number): number {
  return NY + (n - 1) * CYLINDER_BLOCK_SIZE;
}

/**
 * Index of cylinder i's 6th per-surface heat ledger (∫ gas-to-BLOCK heat flow dt, J: the block deck /
 * valve-pocket floor of an 'l-head' chamber, heat-transfer/wall-heat WALL_BLOCK) in an engine of n
 * cylinders. These ledgers exist only for chambers with a block surface and are APPENDED after every
 * cylinder block (one per cylinder), so a flat-disc engine keeps its layout (the CFR F-1: NY, five
 * surface ledgers) and the per-cylinder blocks keep CYLINDER_BLOCK_SIZE.
 */
export function blockLedgerIndex(n: number, i: number): number {
  return stateLength(n) + i;
}

/** Cylinder modes. */
export const MODE_OPEN = 0;
export const MODE_SINGLE = 1;
export const MODE_TWO = 2;
export const MODE_BURNED = 3;

// Event kinds (per cylinder, local crank angle)
export const EV_IVC = 1;
export const EV_EVO = 2;
export const EV_VALVE = 3; // IVO / EVC (step boundary only)
export const EV_BDC_START = 4; // −180
export const EV_BDC_END = 5; // +180
export const EV_TDC = 6; // 0 (instantaneous combustion)
export const EV_IGN = 7; // dwell start / spark (step boundary only)
export const EV_WIEBE = 8;
export const EV_WRAP = 9; // 360
export const EV_END = 10; // closed-cycle-only end
export const MAX_EVENTS = 16;

/** Exhaust outlet (pipe discharging to the ambient) discharge coefficient. UNVERIFIED: 1 (a plain
 * pipe exit has no vena contracta; its pressure drop at the CFR's ~3 g/s is ≈ 10 Pa anyway). */
export const OUTLET_DISCHARGE_COEFFICIENT = 1.0;
/** Burned-zone temperature below which the kinetic NO is frozen (no rate evaluation), K. The
 * extended-Zeldovich rates are negligible there (round-1 value kept). */
export const NO_FREEZE_T = 1000;
/** Relaxation time of the crevice-zone constraint keeper, s (numerical; ≫ the steps, ≪ the cycle). */
export const CREVICE_RELAX_TIME = 0.5e-3;
/** Warm-up: extra cycles (at most) while the intake-plenum mass changes by more than WARMUP_MASS_TOL per cycle. */
export const WARMUP_EXTRA_MAX = 20;
export const WARMUP_MASS_TOL = 1e-3;
/** Burned-mass fraction seeded at the Wiebe start (numerical, see Cylinder.handleEvents EV_WIEBE). */
export const WIEBE_SEED = 1e-5;
/** Lowest crank speed in 'free' mode, rev/min (numerical stall guard; the model needs ω > 0). */
export const FREE_MODE_MIN_RPM = 60;
export const FREE_MODE_MIN_OMEGA = (FREE_MODE_MIN_RPM * 2 * Math.PI) / 60;
/** Initial cylinder / exhaust temperatures of a cold start, K — initial conditions only (the
 * warm-up cycles erase them; not physical constants). */
export const INITIAL_CYLINDER_T = 900;
export const INITIAL_EXHAUST_T = 800;

/**
 * CycleSummary.maxPressureRiseRate = max over θ of [p(θ + Δ) − p(θ)]/Δ with Δ = DP_WINDOW_DEG
 * (the usual finite-difference definition for 0.1°-sampled indicator data), the window position
 * sampled every Δ/DP_SUBDIV so the result does not depend on where a steep (knock) rise falls
 * relative to a coarse grid.
 */
export const DP_WINDOW_DEG = 0.1;
export const DP_SUBDIV = 10;

// =============================================================================================
// Trace
// =============================================================================================

/** One sample per integration step of a recorded cycle (validation trace; cylinder 0). */
export interface CycleTrace {
  theta: number[];
  t: number[];
  volume: number[];
  /** Thermodynamic pressure, Pa. */
  pressure: number[];
  /** Pressure + knock oscillation at the pickup, Pa. */
  pressureReported: number[];
  Tu: number[];
  Tb: number[];
  Tmean: number[];
  mu: number[];
  mb: number[];
  xb: number[];
  me: number[];
  /** Total internal energy of the cylinder charge, J. */
  U: number[];
  /** Cumulative ∫p dV, J (whole run). */
  work: number[];
  /** Cumulative wall heat loss, J (whole run). */
  heatLoss: number[];
  /** Cumulative electrical energy to the gas minus kernel electrode loss, J (whole run). */
  sparkEnergy: number[];
  SL: number[];
  uPrime: number[];
  flameRadius: number[];
  frontArea: number[];
  /** Step-mean burn rate dm_b/dt, kg/s. */
  burnRate: number[];
  heatReleaseRate: number[];
  heatLossRate: number[];
  lwIntegral: number[];
  /** Burned-zone kinetic NO mole fraction. */
  xNO: number[];
  /**
   * Kinetic NO moles carried in the burned-zone thermodynamics by the closure solve of this sample
   * (closure.noKinetic; 0 when options.burnedNOThermo resolves to 'equilibrium').
   */
  nNOClosure: number[];
  mdotIntake: number[];
  mdotExhaust: number[];
  pIntake: number[];
  pExhaust: number[];
  mass: number[];
  mode: number[];
  /** Integral length scale and Taylor microscale, m. */
  L: number[];
  lambda: number[];
}

const TRACE_KEYS: readonly (keyof CycleTrace)[] = [
  'theta', 't', 'volume', 'pressure', 'pressureReported', 'Tu', 'Tb', 'Tmean', 'mu', 'mb', 'xb', 'me', 'U',
  'work', 'heatLoss', 'sparkEnergy', 'SL', 'uPrime', 'flameRadius', 'frontArea', 'burnRate', 'heatReleaseRate',
  'heatLossRate', 'lwIntegral', 'xNO', 'nNOClosure', 'mdotIntake', 'mdotExhaust', 'pIntake', 'pExhaust', 'mass', 'mode', 'L', 'lambda',
];

export function newCycleTrace(): CycleTrace {
  const t = {} as Record<keyof CycleTrace, number[]>;
  for (const k of TRACE_KEYS) t[k] = [];
  return t as CycleTrace;
}

/** Prescribed initial (IVC-like) state of a closed-cycle-only run. */
export interface ClosedCycleInit {
  /** Start angle (the "IVC"), crank deg. */
  startDeg: number;
  /** End angle (the "EVO"), crank deg. */
  endDeg: number;
  /** Charge temperature (K) and pressure (Pa) at startDeg. */
  T: number;
  p: number;
  /** Charge mole fractions (Float64Array(NS), or a species-name map); normalised on use. */
  X: Float64Array | Readonly<Record<string, number>>;
  /** Burned-gas (residual) mass fraction contained in X, for the S_L dilution (default 0). */
  residualMassFraction?: number;
  /**
   * Initial turbulence intensity u′ at startDeg, m/s. Default 0.5 × mean piston speed: the
   * motored-TDC value of Bopp, Vafidis & Whitelaw 1986 (0.45–0.6 S̄p, turbulence.ts); UNVERIFIED as
   * an IVC value (the full-cycle model computes it from the intake jet instead).
   */
  uPrime?: number;
}

/** Per-sub-model CPU time accumulators (ms) when options.profile is set. */
export interface CycleProfile {
  rhsOpen: number;
  rhsClosed: number;
  closure: number;
  flame: number;
  heat: number;
  ignition: number;
  knock: number;
  no: number;
  steps: number;
  rhsEvals: number;
  eqSolves: number;
}

export const now = (): number => performance.now();

// =============================================================================================
// Helpers
// =============================================================================================

/**
 * ∫₀¹ r(u) g(u) du for r exponential (r0 → r1, both ≥ 0) and g linear (g0 → g1): with b = ln(r1/r0),
 * r0 [g0 (e^b − 1)/b + (g1 − g0)(e^b (b − 1) + 1)/b²] (series for small |b|); trapezoid if r0 or r1 is 0.
 */
export function expLinearIntegral(r0: number, r1: number, g0: number, g1: number): number {
  if (!(r0 > 0 && r1 > 0)) return 0.5 * (r0 * g0 + r1 * g1);
  const b = Math.log(r1 / r0);
  const dg = g1 - g0;
  if (Math.abs(b) < 1e-4) return r0 * (g0 * (1 + b / 2 + (b * b) / 6) + dg * (0.5 + b / 3 + (b * b) / 8));
  const eb = Math.exp(b);
  return r0 * (g0 * ((eb - 1) / b) + (dg * (eb * (b - 1) + 1)) / (b * b));
}

/** Scratch: argument u of the last hermiteMax. */
export const HSCR = new Float64Array(1);

/** Cubic Hermite P(u) on [0, 1] with P(0) = p0, P′(0) = a0, P(1) = p1, P′(1) = a1. */
export function hermiteAt(p0: number, a0: number, p1: number, a1: number, u: number): number {
  const d = p1 - p0;
  const c2 = 3 * d - 2 * a0 - a1;
  const c3 = a0 + a1 - 2 * d;
  return p0 + u * (a0 + u * (c2 + u * c3));
}

/** Maximum of the cubic Hermite interpolant over u ∈ [0, 1] (argument in HSCR[0]). */
export function hermiteMax(p0: number, a0: number, p1: number, a1: number): number {
  let best = p0;
  HSCR[0] = 0;
  if (p1 > best) {
    best = p1;
    HSCR[0] = 1;
  }
  const d = p1 - p0;
  const c2 = 3 * d - 2 * a0 - a1;
  const c3 = a0 + a1 - 2 * d;
  // P′(u) = a0 + 2 c2 u + 3 c3 u²
  const A = 3 * c3;
  const B = 2 * c2;
  if (Math.abs(A) <= 1e-14 * (Math.abs(B) + Math.abs(a0))) {
    if (B !== 0) best = hermiteTry(p0, a0, c2, c3, -a0 / B, best);
  } else {
    const disc = B * B - 4 * A * a0;
    if (disc >= 0) {
      const sq = Math.sqrt(disc);
      const q = -0.5 * (B + (B >= 0 ? sq : -sq));
      if (q !== 0) {
        best = hermiteTry(p0, a0, c2, c3, q / A, best);
        best = hermiteTry(p0, a0, c2, c3, a0 / q, best);
      }
    }
  }
  return best;
}

function hermiteTry(p0: number, a0: number, c2: number, c3: number, u: number, best: number): number {
  if (!(u > 0 && u < 1)) return best;
  const v = p0 + u * (a0 + u * (c2 + u * c3));
  if (v > best) {
    HSCR[0] = u;
    return v;
  }
  return best;
}

/** Angle at which x crosses `lev` inside a step (linear in x), NaN if it does not. */
export function crossing(x0: number, x1: number, th0: number, th1: number, lev: number): number {
  return x1 >= lev && x0 < lev ? th0 + ((lev - x0) / (x1 - x0)) * (th1 - th0) : NaN;
}

export function cloneOp(op: OperatingPoint): OperatingPoint {
  const o = { ...op, fuel: { ...op.fuel } };
  if (op.load) o.load = { ...op.load };
  return o;
}

export function massOf(N: Float64Array): number {
  let m = 0;
  for (let k = 0; k < NS; k++) m += N[k] * MOLAR_MASS[k];
  return m;
}

export function elementMolesOf(N: Float64Array, b: Float64Array): Float64Array {
  b.fill(0);
  // ELEMENT_COUNTS via the thermo helper would allocate nothing either; inline for clarity
  for (let k = 0; k < NS; k++) {
    const n = N[k];
    if (n === 0) continue;
    b[EL.C] += n * EC[k * NE + EL.C];
    b[EL.H] += n * EC[k * NE + EL.H];
    b[EL.O] += n * EC[k * NE + EL.O];
    b[EL.N] += n * EC[k * NE + EL.N];
    b[EL.AR] += n * EC[k * NE + EL.AR];
  }
  return b;
}

const EC = new Float64Array(NS * NE);
for (let k = 0; k < NS; k++) for (let e = 0; e < NE; e++) EC[k * NE + e] = ELEMENT_COUNTS[k][e];

/**
 * Replace the NO moles of a burned-gas mole vector by nNO conserving elements (closure.ts
 * noSwapDeltas: N from N2; O from O2, then from CO2 → CO and H2O → H2 in proportion; a decrease
 * returns N2 + O2). The same swap is the kinetic-NO energy correction of the closure.
 */
export function swapNO(N: Float64Array, nNO: number): void {
  noSwapDeltas(N, nNO, SWAP_DN);
  N[SP.NO] += SWAP_DN[SP.NO];
  N[SP.N2] += SWAP_DN[SP.N2];
  N[SP.O2] += SWAP_DN[SP.O2];
  N[SP.CO2] += SWAP_DN[SP.CO2];
  N[SP.CO] += SWAP_DN[SP.CO];
  N[SP.H2O] += SWAP_DN[SP.H2O];
  N[SP.H2] += SWAP_DN[SP.H2];
}
const SWAP_DN = new Float64Array(NS);

export function resolveDelayModel(m: CycleModelOptions['ignitionDelayModel']): IgnitionDelayModel {
  if (typeof m === 'object') return m;
  if (m === 'douaud-eyzat') return douaudEyzat;
  if (m === 'douaud-eyzat-llnl') return douaudEyzatLLNL;
  if (m === 'llnl-prf-v2') return prfLLNLv2;
  return prfLLNLGasoline2011;
}
