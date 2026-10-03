/**
 * Chemical equilibrium of the burned gas by Gibbs-energy minimisation over the N_EQ product
 * species (ideal gas, no condensed phases).
 *
 * Formulation — NASA CEA (Gordon & McBride 1994, NASA RP-1311, ch. 2–3):
 *   unknowns ln n_j (every active species), ln n (total moles, a separate variable), and the
 *   modified Lagrange multipliers π_i = λ_i/(RT) of the active elements. With
 *     μ_j/RT = g°_j/RT + ln(n_j/n) + ln(p/p°)                        (ideal gas, RP-1311 §2.2)
 *   one Newton step at fixed (T, p) solves the reduced (l+1)×(l+1) system
 *     Σ_i Σ_j a_kj a_ij n_j π_i + Σ_j a_kj n_j Δln n = b°_k − b_k + Σ_j a_kj n_j μ_j/RT   (2.24)
 *     Σ_i Σ_j a_ij n_j π_i + (Σ_j n_j − n) Δln n = n − Σ_j n_j + Σ_j n_j μ_j/RT          (2.26)
 *   and back-substitutes  Δln n_j = −μ_j/RT + Σ_i a_ij π_i + Δln n       (RP-1311 §2.3)
 *   with the RP-1311 step control (eqs. 3.1–3.3): λ1 = 2/max(5|Δln n|, Δln n_j) over species
 *   with n_j/n > 1e-8 growing, λ2 keeps a growing trace species (n_j/n ≤ 1e-8) from passing
 *   1e-4 in one step, λ = min(1, λ1, λ2).
 *   At fixed (T, V) (Helmholtz minimisation, RP-1311 eqs. 2.45 ff.) μ_j/RT = g°_j/RT +
 *   ln(n_j R T/(V p°)), there is no ln n unknown and the system is l×l.
 *   The RP-1311 scan on NTRS has no text layer, so the equation numbers, constants and
 *   λ rules were checked against NASA's open-source CEA re-implementation, which cites them
 *   (github.com/nasa/cea, source/equilibrium.f90, fetched 2026-09-29: `size = 18.420681`,
 *   `FACTOR = −9.2103404`, EqSolver_compute_damped_update_factor "(Eq. 3.1)/(Eq. 3.2)/(Eq. 3.3)",
 *   matrix rows "(2.24/2.45)", "(2.26)", "(2.27)/(2.28)/(2.47)/(2.48)", partials "Table 2.3",
 *   "Table 2.4", "(Eq. 2.59)", γ_s "(Eq. 2.71/2.73)"), and re-checked (review, 2026-09-29) against
 *   the OCR text layer of the RP-1311 copy at shepherd.caltech.edu/EDL/PublicResources/sdt/refs/
 *   NASA-RP-1311-1.pdf: SIZE "−ln 10^−8 = 18.420681", eq. 3.1 λ1 = 2/max(5|Δln T|, 5|Δln n|,
 *   Δln n_j), eq. 3.2 λ2 with "9.2103404" (limits a trace species to X ≤ 1e-4 per step), eq. 3.3
 *   λ = min(1, λ1, λ2), eqs. 2.24/2.26 (TP reduced system), 2.50/2.51, 2.59, 2.70–2.73.
 *
 * Why this is robust at 300–4000 K with trace species far below 1e-30: the unknowns are
 * logarithms, and after every full step each species satisfies the stationarity condition
 * exactly for the current π, so trace species are exactly as converged as the element
 * potentials (no floor is needed; ln n_j is kept finite by a clamp at ln(n_j/n) ≥ −1400).
 *
 * Energy/entropy-constrained problems (HP, SP, UV) use an OUTER safeguarded Newton iteration
 * on T around the inner TP (or TV) solve; the slope is the analytic EQUILIBRIUM heat capacity
 * from the RP-1311 §2.5 derivative systems (Tables 2.3/2.4, eq. 2.59), and each new inner
 * solve is warm-started by the first-order predictor ln n_j += (∂ln n_j/∂ln T) Δln T.
 *
 * Internals are normalised to Σ_e b_e = 1 mol of atoms (everything is homogeneous of degree 1
 * in b), so tolerances are absolute on mole-fraction-like numbers.
 *
 * Species of absent elements are removed (e.g. no carbon → no CO/CO2); elements whose share of
 * the atoms is below EQ_B_REL_MIN are treated as absent. Carbon needs oxygen (only CO/CO2 carry
 * it): b_O ≤ b_C has no solution in this species set and returns converged = false.
 * Invalid inputs — a NaN/∞ element amount, a negative one beyond round-off (EQ_B_NEG_TOL),
 * T, p or V not finite and positive, a NaN/∞ H/U/S target — return converged = false with NaN
 * mole fractions instead of a silently "converged" answer.
 *
 * Hot path: after construction nothing is allocated. The returned EqResult is owned by the
 * solver and overwritten by the next call.
 */
import { R_UNIVERSAL } from '../core/constants';
import { EL, ELEMENT_COUNTS, N_EQ, NE, NS, SP } from '../core/species';
import { evalAllNondim, MOLAR_MASS, P_REF_THERMO } from '../thermo/thermo';
import { DENSE_STRIDE, luFactor, luSolve, scaleSymmetric } from './dense';

/** RP-1311 eq. (3.2) threshold −ln(1e-8): species with n_j/n below e^−SIZE are "trace". */
const SIZE = 18.420681;
/** RP-1311 eq. (3.2): ln(1e-4), the ceiling a growing trace species may reach in one step. */
const LN_1EM4 = -9.2103404;
/** ln(n_j/n) clamp that keeps ln n_j finite (e^−1400 underflows to 0 anyway). */
const LN_X_FLOOR = -1400;
/** Elements with b_e/Σb below this are treated as absent. */
export const EQ_B_REL_MIN = 1e-20;
/** Negative element amounts down to −EQ_B_NEG_TOL·Σb (round-off of an upstream difference) are
 *  treated as absent; anything more negative is non-physical and makes the problem infeasible. */
export const EQ_B_NEG_TOL = 1e-12;
/** Lower/upper temperature bracket of the HP/UV/SP outer iteration, K (NASA-7 fit range). */
export const EQ_T_MIN = 200;
export const EQ_T_MAX = 6000;
/** Newton convergence: every species' full step |Δln n_j| ≤ STEP_TOL ... */
const STEP_TOL = 5e-7;
/** ... or its mole-fraction change X_j |Δln n_j| ≤ ABS_TOL (species at the round-off level
 *  of the element balance, e.g. O2/CO/H2 of an exactly stoichiometric mixture at 300 K). */
const ABS_TOL = 1e-14;
/** Outer (temperature) iteration tolerance |ΔT|/T. The last Newton correction is applied by
 *  the first-order composition predictor, so the returned state misses its H/U/S target by
 *  O((T_TOL·T)²) only. */
const T_TOL = 1e-6;
/** Largest |Δln T| of one outer Newton step. */
const MAX_DLNT = 0.4;
/** Largest |Δln T| over which the first-order composition predictor is applied. */
const MAX_PREDICT_DLNT = 0.25;
/** Total moles of the CEA-style uniform cold start per mol of atoms. NASA CEA starts from
 *  n = 0.1 kmol/kg split evenly over the gas species (equilibrium.f90, EqSolution_init), ≈ 1.4 mol
 *  per mol of atoms for air; any O(1) value works (only a starting point). */
const COLD_N0 = 1;
const MAX_ITER_WARM = 60;
const MAX_ITER_COLD = 400;
const MAX_OUTER = 80;

const S = DENSE_STRIDE;

/** Result of an equilibrium solve. Owned by the solver: overwritten by the next call. */
export interface EqResult {
  /** Species moles, mol, Float64Array(NS) in SPECIES order; fuel entries (and species of
   *  absent elements) are 0. */
  readonly N: Float64Array;
  /** Mole fractions, Float64Array(NS). */
  readonly X: Float64Array;
  /** Temperature, K. */
  T: number;
  /** Pressure, Pa (for UV: n R T / V). */
  p: number;
  /** Volume n R T / p, m³. */
  V: number;
  /** Total moles Σ N_k, mol. */
  nTotal: number;
  /** Mass Σ N_k M_k, kg. */
  mass: number;
  /** Extensive enthalpy (absolute, NASA convention), J. */
  H: number;
  /** Extensive internal energy, J. */
  U: number;
  /** Extensive entropy, J/K. */
  S: number;
  /** True when the Newton iteration(s) met the tolerances. */
  converged: boolean;
  /** Total inner Newton iterations spent in this call. */
  iterations: number;
  /** Outer temperature iterations (HP/UV/SP; 0 for TP). */
  outerIterations: number;
  /** Element potentials π_e = −λ_e/(R T) (dimensionless; RP-1311's modified Lagrange multipliers,
   *  λ_e from eq. 2.8/2.10, so that μ_j/RT = Σ_e a_je π_e at equilibrium),
   *  Float64Array(NE) in ELEMENTS order; 0 for absent elements. */
  readonly pi: Float64Array;
}

/** Equilibrium derivatives and heat capacities of the current state (normalised internals). */
interface DerivScratch {
  /** (∂ln n/∂ln T)_p, (∂ln n/∂ln p)_T */
  dnT: number;
  dnP: number;
  /** Equilibrium cp/R and cv/R of the normalised system (per mol of atoms). */
  cpR: number;
  cvR: number;
}

/**
 * Gibbs-minimisation equilibrium solver for the burned-gas species [0, N_EQ). Holds its scratch
 * buffers and the last solution (used as the warm start of the next call when the set of
 * present elements is unchanged). Not re-entrant; use one instance per caller.
 */
export class EquilibriumSolver {
  /** Result of the most recent solve (overwritten by the next call). */
  readonly result: EqResult = {
    N: new Float64Array(NS),
    X: new Float64Array(NS),
    T: 0,
    p: 0,
    V: 0,
    nTotal: 0,
    mass: 0,
    H: 0,
    U: 0,
    S: 0,
    converged: false,
    iterations: 0,
    outerIterations: 0,
    pi: new Float64Array(NE),
  };

  // ---- active set (compact indices) -----------------------------------------------------
  private mask = -1;
  private nAct = 0;
  private nEl = 0;
  private readonly act = new Int32Array(N_EQ);
  private readonly elIdx = new Int32Array(NE);
  /** a[jj*NE + ii]: atoms of active element ii in active species jj. */
  private readonly a = new Float64Array(N_EQ * NE);
  /** Normalised element amounts (compact), Σ = 1 over all elements. */
  private readonly b0 = new Float64Array(NE);
  /** Normalisation: internal = physical × scale (scale = 1/Σb). */
  private scale = 1;
  private feasible = true;

  // ---- state (normalised) -----------------------------------------------------------------
  private readonly lnN = new Float64Array(N_EQ);
  private readonly nj = new Float64Array(N_EQ);
  private lnNt = 0;
  private hasState = false;
  private stateMask = -1;
  /** True if the state came from a TV solve (then lnNt must be refreshed before a TP solve). */
  private stateIsTV = false;

  // ---- species thermo at the cached temperature ----------------------------------------------
  private Tth = NaN;
  private readonly cpR = new Float64Array(NS);
  private readonly hRT = new Float64Array(NS);
  private readonly s0R = new Float64Array(NS);
  private readonly g = new Float64Array(NS);

  // ---- Newton scratch -----------------------------------------------------------------------
  private readonly M = new Float64Array(S * S);
  private readonly r = new Float64Array(S);
  private readonly r2 = new Float64Array(S);
  private readonly x = new Float64Array(S);
  private readonly x2 = new Float64Array(S);
  private readonly tmp = new Float64Array(S);
  private readonly d = new Float64Array(S);
  private readonly piv = new Int32Array(S);
  private readonly mu = new Float64Array(N_EQ);
  private readonly dln = new Float64Array(N_EQ);
  private readonly piC = new Float64Array(NE);

  // ---- derivatives of the current state ----------------------------------------------------
  /** (∂ln n_j/∂ln T)_p (TP state) or (∂ln n_j/∂ln T)_V (TV state). */
  private readonly dlnT = new Float64Array(N_EQ);
  /** (∂ln n_j/∂ln p)_T (TP state only). */
  private readonly dlnP = new Float64Array(N_EQ);
  private readonly der: DerivScratch = { dnT: 0, dnP: 0, cpR: 0, cvR: 0 };
  /** ∂π_i/∂ln T of the last derivTP/derivTV (compact element order). */
  private readonly dpiT = new Float64Array(NE);

  /** Forget the warm start (the next solve starts cold). */
  reset(): void {
    this.hasState = false;
  }

  // =========================================================================================
  // Public solves
  // =========================================================================================

  /**
   * Equilibrium at fixed temperature T (K) and pressure p (Pa) for element amounts b
   * (mol of atoms, Float64Array(NE) in ELEMENTS order). Warm-started from the previous solution.
   */
  solveTP(b: Float64Array, T: number, p: number): EqResult {
    const res = this.result;
    res.outerIterations = 0;
    res.iterations = 0;
    // invalid state (NaN, ≤ 0, ∞) or infeasible b: converged = false with NaN X (never stale)
    if (!(T > 0 && T < Infinity && p > 0 && p < Infinity) || !this.setup(b)) return this.infeasible(T, p);
    this.setThermo(T);
    const ok = this.innerTP(Math.log(p / P_REF_THERMO));
    this.writeResult(T, p, ok);
    return res;
  }

  /**
   * Equilibrium at fixed volume V (m³) and temperature T (K) for element amounts b (mol).
   * Returns the pressure n R T / V in result.p.
   */
  solveTV(b: Float64Array, T: number, V: number): EqResult {
    const res = this.result;
    res.outerIterations = 0;
    res.iterations = 0;
    if (!(T > 0 && T < Infinity && V > 0 && V < Infinity) || !this.setup(b)) return this.infeasible(T, NaN);
    this.setThermo(T);
    const ok = this.innerTV(Math.log(V * this.scale));
    const p = (this.sumN() * R_UNIVERSAL * T) / (V * this.scale);
    this.writeResult(T, p, ok);
    return res;
  }

  /**
   * Equilibrium at fixed total enthalpy H (J, absolute/NASA convention, for the moles b) and
   * pressure p (Pa). Tguess (K) seeds the outer Newton iteration on T.
   */
  solveHP(b: Float64Array, H: number, p: number, Tguess = 2000): EqResult {
    return this.solveOuter(0, b, H, p, Tguess);
  }

  /**
   * Equilibrium at fixed total internal energy U (J) and volume V (m³) for the moles b.
   * result.p = n R T / V.
   */
  solveUV(b: Float64Array, U: number, V: number, Tguess = 2000): EqResult {
    return this.solveOuter(1, b, U, V, Tguess);
  }

  /** Equilibrium at fixed total entropy S (J/K) and pressure p (Pa) for the moles b. */
  solveSP(b: Float64Array, Sx: number, p: number, Tguess = 2000): EqResult {
    return this.solveOuter(2, b, Sx, p, Tguess);
  }

  // =========================================================================================
  // Derivatives / properties of the last converged state
  // =========================================================================================

  /**
   * Equilibrium derivatives of the most recent solution (any solve kind; a UV/TV state is
   * evaluated at its own pressure). Writes into `out` and returns it:
   * mass-specific h, u (J/kg), s (J/(kg K)), cp, cv (J/(kg K), EQUILIBRIUM — composition
   * shifting), cpFrozen (J/(kg K)), M (kg/mol), rho (kg/m³), dlnV_dlnT = (∂lnV/∂lnT)_p,
   * dlnV_dlnp = (∂lnV/∂lnp)_T, gamma = cp/cv, gammaS = (∂ln p/∂ln ρ)_s, soundSpeed (m/s).
   * Derivative systems RP-1311 Tables 2.3/2.4 and cp,eq eq. (2.59) (numbers as cited in NASA
   * CEA equilibrium.f90); (∂lnV/∂lnT)_p = 1 + (∂ln n/∂ln T)_p, (∂lnV/∂lnp)_T = −1 + (∂ln n/∂ln p)_T
   * (RP-1311 eqs. 2.50/2.51, from V = nRT/p);
   * cv = cp + (pV/T)(∂lnV/∂lnT)_p²/(∂lnV/∂lnp)_T (RP-1311 eq. 2.70) and γ_s = −γ/(∂lnV/∂lnp)_T
   * (eqs. 2.71/2.73), as implemented in CEA EqPartials_compute_partials. Equation numbers checked
   * against the RP-1311 text layer of shepherd.caltech.edu/EDL/PublicResources/sdt/refs/
   * NASA-RP-1311-1.pdf (fetched 2026-09-29).
   */
  properties(out: EqProperties = newEqProperties()): EqProperties {
    const res = this.result;
    const T = res.T;
    const p = res.p;
    if (!this.feasible || !this.hasState) {
      out.valid = false;
      return out;
    }
    this.setThermo(T);
    const lnP = Math.log(p / P_REF_THERMO);
    if (!this.derivTP()) {
      out.valid = false;
      return out;
    }
    const nA = this.nAct;
    const act = this.act;
    let m = 0;
    let n = 0;
    let sh = 0;
    let ss = 0;
    let scp = 0;
    for (let jj = 0; jj < nA; jj++) {
      const k = act[jj];
      const w = this.nj[jj];
      n += w;
      m += w * MOLAR_MASS[k];
      sh += w * this.hRT[k];
      scp += w * this.cpR[k];
    }
    const lnn = Math.log(n);
    for (let jj = 0; jj < nA; jj++) {
      const w = this.nj[jj];
      if (w === 0) continue;
      const k = act[jj];
      ss += w * (this.s0R[k] - (this.lnN[jj] - lnn) - lnP);
    }
    const R = R_UNIVERSAL;
    const der = this.der;
    const invm = 1 / m;
    out.valid = res.converged;
    out.T = T;
    out.p = p;
    for (let k = 0; k < NS; k++) out.X[k] = res.X[k];
    out.M = m / n;
    out.rho = (p * out.M) / (R * T);
    out.h = R * T * sh * invm;
    out.u = R * T * (sh - n) * invm;
    out.s = R * ss * invm;
    out.cp = R * der.cpR * invm;
    out.cpFrozen = R * scp * invm;
    out.dlnV_dlnT = 1 + der.dnT;
    out.dlnV_dlnp = -1 + der.dnP;
    out.cv = out.cp + ((n * R * invm) * out.dlnV_dlnT * out.dlnV_dlnT) / out.dlnV_dlnp;
    out.gamma = out.cp / out.cv;
    out.gammaS = -out.gamma / out.dlnV_dlnp;
    out.soundSpeed = Math.sqrt((out.gammaS * n * R * T) * invm);
    return out;
  }

  /**
   * (∂ln N_k/∂ln T)_p and (∂ln N_k/∂ln p)_T of every species at the last solution (analytic,
   * RP-1311 Tables 2.3/2.4). Writes Float64Array(NS) outputs (0 for absent/fuel species).
   * Returns false if the derivative system is singular.
   */
  /**
   * (∂ln N_k/∂ln T)_p and (∂ln N_k/∂ln p)_T of ONE species k from the derivative system of the last
   * properties() / speciesDerivatives() call — no re-solve (hot-path accessor; 0 for inactive
   * species or before any derivative call). Writes [dlnT, dlnp] into out and returns it.
   */
  lastSpeciesLogDerivatives(k: number, out: Float64Array): Float64Array {
    out[0] = 0;
    out[1] = 0;
    for (let jj = 0; jj < this.nAct; jj++) {
      if (this.act[jj] === k) {
        out[0] = this.dlnT[jj];
        out[1] = this.dlnP[jj];
        break;
      }
    }
    return out;
  }

  speciesDerivatives(dlnNdlnT: Float64Array, dlnNdlnp: Float64Array): boolean {
    if (!this.feasible || !this.hasState) return false;
    this.setThermo(this.result.T);
    if (!this.derivTP()) return false;
    dlnNdlnT.fill(0);
    dlnNdlnp.fill(0);
    for (let jj = 0; jj < this.nAct; jj++) {
      dlnNdlnT[this.act[jj]] = this.dlnT[jj];
      dlnNdlnp[this.act[jj]] = this.dlnP[jj];
    }
    return true;
  }

  // =========================================================================================
  // Setup
  // =========================================================================================

  /**
   * Normalise b, determine the active element/species sets. Returns feasibility: false for a
   * non-finite entry (NaN/±Inf — never silently dropped), a negative entry below
   * −EQ_B_NEG_TOL·Σb (non-physical), no positive entry, or carbon without enough oxygen.
   * Round-off-level negatives (≥ −EQ_B_NEG_TOL·Σb) are treated as absent, like positive
   * amounts below EQ_B_REL_MIN·Σb.
   */
  private setup(b: Float64Array): boolean {
    let sb = 0;
    let neg = 0;
    for (let e = 0; e < NE; e++) {
      const v = b[e];
      if (v > 0) sb += v;
      else if (v < 0) neg -= v;
      else if (v !== 0) sb = NaN; // NaN entry
    }
    if (!(sb > 0) || !Number.isFinite(sb) || !(neg <= EQ_B_NEG_TOL * sb)) {
      this.feasible = false;
      return false;
    }
    const s = 1 / sb;
    this.scale = s;
    let mask = 0;
    for (let e = 0; e < NE; e++) if (b[e] * s > EQ_B_REL_MIN) mask |= 1 << e;
    // Carbon can only be carried by CO/CO2: it needs more O atoms than C atoms.
    if (mask & (1 << EL.C) && !(b[EL.O] > b[EL.C] * (1 + 1e-12))) {
      this.feasible = false;
      return false;
    }
    this.feasible = true;
    if (mask !== this.mask) this.buildActive(mask);
    for (let ii = 0; ii < this.nEl; ii++) this.b0[ii] = b[this.elIdx[ii]] * s;
    return true;
  }

  private buildActive(mask: number): void {
    this.mask = mask;
    let ne = 0;
    for (let e = 0; e < NE; e++) if (mask & (1 << e)) this.elIdx[ne++] = e;
    this.nEl = ne;
    let na = 0;
    for (let k = 0; k < N_EQ; k++) {
      const row = ELEMENT_COUNTS[k];
      let ok = true;
      for (let e = 0; e < NE; e++) if (row[e] !== 0 && !(mask & (1 << e))) ok = false;
      if (!ok) continue;
      this.act[na] = k;
      for (let ii = 0; ii < ne; ii++) this.a[na * NE + ii] = row[this.elIdx[ii]];
      na++;
    }
    this.nAct = na;
  }

  /** Species thermo (cp/R, h/RT, s°/R, g°/RT) at T, cached. */
  private setThermo(T: number): void {
    if (T === this.Tth) return;
    evalAllNondim(T, this.cpR, this.hRT, this.s0R);
    for (let k = 0; k < N_EQ; k++) this.g[k] = this.hRT[k] - this.s0R[k];
    this.Tth = T;
  }

  private sumN(): number {
    let s = 0;
    for (let jj = 0; jj < this.nAct; jj++) s += Math.exp(this.lnN[jj]);
    return s;
  }

  /** CEA-style uniform start (kind 0) or complete-combustion majors + 1e-6 minors (kind 1). */
  private coldStart(kind: number): void {
    const nA = this.nAct;
    if (kind === 0) {
      const lnj = Math.log(COLD_N0 / nA);
      for (let jj = 0; jj < nA; jj++) this.lnN[jj] = lnj;
      this.lnNt = Math.log(COLD_N0);
    } else {
      // Majors from the element pools (lean: CO2/H2O/O2; rich: CO first, then H2O, CO2, H2).
      const b = this.b0;
      let bC = 0;
      let bH = 0;
      let bO = 0;
      let bN = 0;
      let bAr = 0;
      for (let ii = 0; ii < this.nEl; ii++) {
        const e = this.elIdx[ii];
        if (e === EL.C) bC = b[ii];
        else if (e === EL.H) bH = b[ii];
        else if (e === EL.O) bO = b[ii];
        else if (e === EL.N) bN = b[ii];
        else if (e === EL.AR) bAr = b[ii];
      }
      let co = bC;
      let oLeft = bO - bC;
      const h2o = Math.min(bH / 2, Math.max(oLeft, 0));
      oLeft -= h2o;
      const co2 = Math.min(co, Math.max(oLeft, 0));
      co -= co2;
      oLeft -= co2;
      const h2 = bH / 2 - h2o;
      const o2 = Math.max(oLeft, 0) / 2;
      const tot = bN / 2 + bAr + co + co2 + h2o + h2 + o2;
      const floor = 1e-6 * tot;
      let nt = 0;
      for (let jj = 0; jj < nA; jj++) {
        const k = this.act[jj];
        let v = 0;
        if (k === SP.N2) v = bN / 2;
        else if (k === SP.O2) v = o2;
        else if (k === SP.AR) v = bAr;
        else if (k === SP.CO2) v = co2;
        else if (k === SP.H2O) v = h2o;
        else if (k === SP.CO) v = co;
        else if (k === SP.H2) v = h2;
        v = Math.max(v, floor);
        this.lnN[jj] = Math.log(v);
        nt += v;
      }
      this.lnNt = Math.log(nt);
    }
    this.hasState = true;
    this.stateMask = this.mask;
  }

  // =========================================================================================
  // Inner solves with fallbacks
  // =========================================================================================

  /** TP solve at the cached temperature with warm start + cold fallbacks. */
  private innerTP(lnP: number): boolean {
    const res = this.result;
    if (this.hasState && this.stateMask === this.mask) {
      if (this.stateIsTV) this.lnNt = Math.log(this.sumN());
      const it = this.newtonTP(lnP, MAX_ITER_WARM);
      res.iterations += Math.abs(it);
      if (it > 0) return this.finishTP();
    }
    for (let kind = 0; kind < 2; kind++) {
      this.coldStart(kind);
      const it = this.newtonTP(lnP, MAX_ITER_COLD);
      res.iterations += Math.abs(it);
      if (it > 0) return this.finishTP();
    }
    this.finishTP();
    return false;
  }

  /** TV solve (lnVh = ln of the normalised volume) at the cached temperature. */
  private innerTV(lnVh: number): boolean {
    const res = this.result;
    if (this.hasState && this.stateMask === this.mask) {
      const it = this.newtonTV(lnVh, MAX_ITER_WARM);
      res.iterations += Math.abs(it);
      if (it > 0) return this.finishTV();
    }
    for (let kind = 0; kind < 2; kind++) {
      this.coldStart(kind);
      const it = this.newtonTV(lnVh, MAX_ITER_COLD);
      res.iterations += Math.abs(it);
      if (it > 0) return this.finishTV();
    }
    this.finishTV();
    return false;
  }

  private finishTP(): boolean {
    this.stateIsTV = false;
    return true;
  }

  /** A TV state has no separate ln n variable: set it to ln Σn_j for later TP use. */
  private finishTV(): boolean {
    this.stateIsTV = true;
    this.lnNt = Math.log(this.sumN());
    return true;
  }

  // =========================================================================================
  // Newton iterations
  // =========================================================================================

  /**
   * Symmetric-scaled LU solve of the dim×dim system in M with right-hand side r → x.
   * `nRow` (≥ 0) marks the ln n row whose diagonal is ~0; it is scaled by 1/√(Σn_j).
   * Returns false if singular.
   */
  private factorScaled(dim: number, nRow: number, sumN: number): boolean {
    const M = this.M;
    const d = this.d;
    for (let i = 0; i < dim; i++) {
      const v = M[i * S + i];
      d[i] = i === nRow ? 1 / Math.sqrt(sumN) : v > 0 ? 1 / Math.sqrt(v) : 1;
    }
    scaleSymmetric(M, dim, d);
    return luFactor(M, dim, this.piv);
  }

  /** Solve with the factors from factorScaled: rhs (unscaled) → sol (unscaled). */
  private solveScaled(dim: number, rhs: Float64Array, sol: Float64Array): void {
    const d = this.d;
    const t = this.tmp;
    for (let i = 0; i < dim; i++) t[i] = rhs[i] * d[i];
    luSolve(this.M, dim, this.piv, t, sol);
    for (let i = 0; i < dim; i++) sol[i] *= d[i];
  }

  /**
   * Newton iterations at fixed (T, p) (T = cached thermo temperature). Returns +iterations on
   * convergence, −iterations on failure.
   */
  private newtonTP(lnP: number, maxIt: number): number {
    const nA = this.nAct;
    const nE = this.nEl;
    const dim = nE + 1;
    const act = this.act;
    const a = this.a;
    const g = this.g;
    const lnN = this.lnN;
    const nj = this.nj;
    const mu = this.mu;
    const dln = this.dln;
    const M = this.M;
    const r = this.r;
    const x = this.x2;
    const b0 = this.b0;
    for (let it = 1; it <= maxIt; it++) {
      const lnNt = this.lnNt;
      const nt = Math.exp(lnNt);
      let sumN = 0;
      for (let jj = 0; jj < nA; jj++) {
        const v = lnN[jj];
        const w = Math.exp(v);
        nj[jj] = w;
        sumN += w;
        mu[jj] = g[act[jj]] + v - lnNt + lnP;
      }
      // ---- assemble (2.24)/(2.26) ----
      for (let i = 0; i < dim; i++) {
        r[i] = 0;
        for (let j = 0; j < dim; j++) M[i * S + j] = 0;
      }
      let rn = 0;
      for (let jj = 0; jj < nA; jj++) {
        const w = nj[jj];
        if (w === 0) continue;
        const m = mu[jj];
        const base = jj * NE;
        for (let ii = 0; ii < nE; ii++) {
          const aw = a[base + ii] * w;
          if (aw === 0) continue;
          for (let kk = 0; kk <= ii; kk++) M[ii * S + kk] += aw * a[base + kk];
          M[ii * S + nE] += aw;
          r[ii] += aw * m;
        }
        rn += w * m;
      }
      let balBad = false;
      for (let ii = 0; ii < nE; ii++) {
        for (let kk = 0; kk < ii; kk++) M[kk * S + ii] = M[ii * S + kk];
        const bi = M[ii * S + nE];
        M[nE * S + ii] = bi;
        const res = b0[ii] - bi;
        r[ii] += res;
        if (Math.abs(res) > 1e-6 * b0[ii]) balBad = true;
      }
      M[nE * S + nE] = sumN - nt;
      r[nE] = rn + nt - sumN;
      if (!this.factorScaled(dim, nE, sumN)) return -it;
      this.solveScaled(dim, r, x);
      const dlnn = x[nE];
      if (!Number.isFinite(dlnn)) return -it;
      // ---- Δln n_j (back-substitution) and step control (RP-1311 eqs. 3.1–3.3) ----
      let l1 = 5 * Math.abs(dlnn);
      let lam2 = 1;
      for (let jj = 0; jj < nA; jj++) {
        const base = jj * NE;
        let s = dlnn - mu[jj];
        for (let ii = 0; ii < nE; ii++) s += a[base + ii] * x[ii];
        dln[jj] = s;
        if (s > 0) {
          const lx = lnN[jj] - lnNt;
          if (lx <= -SIZE) {
            const den = Math.abs(s - dlnn);
            if (den >= SIZE + LN_1EM4) {
              const l2 = Math.abs(LN_1EM4 - lx) / den;
              if (l2 < lam2) lam2 = l2;
            }
          } else if (s > l1) l1 = s;
        }
      }
      let lam = l1 > 2 ? 2 / l1 : 1;
      if (lam2 < lam) lam = lam2;
      // ---- convergence (full step small for every species) ----
      let conv = lam === 1 && !balBad && Math.abs(dlnn) <= STEP_TOL;
      if (conv) {
        const invSum = 1 / sumN;
        for (let jj = 0; jj < nA; jj++) {
          const dj = dln[jj];
          const ad = Math.abs(dj);
          // judge the absolute change on the larger of the old and new amounts
          if (ad > STEP_TOL && nj[jj] * invSum * ad * (dj > 0 ? Math.exp(dj) : 1) > ABS_TOL) {
            conv = false;
            break;
          }
        }
      }
      // ---- apply ----
      const lnNtNew = lnNt + lam * dlnn;
      for (let jj = 0; jj < nA; jj++) {
        let v = lnN[jj] + lam * dln[jj];
        if (v < lnNtNew + LN_X_FLOOR) v = lnNtNew + LN_X_FLOOR;
        lnN[jj] = v;
      }
      this.lnNt = lnNtNew;
      if (!Number.isFinite(lnNtNew)) return -it;
      if (conv) {
        for (let ii = 0; ii < nE; ii++) this.piC[ii] = x[ii];
        return it;
      }
    }
    return -maxIt;
  }

  /**
   * Newton iterations at fixed (T, V) (Helmholtz minimisation, RP-1311 eq. 2.45 ff.):
   * μ_j/RT = g°_j/RT + ln n_j + ln(R T/(V̂ p°)). Returns ±iterations as newtonTP.
   */
  private newtonTV(lnVh: number, maxIt: number): number {
    const nA = this.nAct;
    const nE = this.nEl;
    const dim = nE;
    const act = this.act;
    const a = this.a;
    const g = this.g;
    const lnN = this.lnN;
    const nj = this.nj;
    const mu = this.mu;
    const dln = this.dln;
    const M = this.M;
    const r = this.r;
    const x = this.x2;
    const b0 = this.b0;
    const lnC = Math.log((R_UNIVERSAL * this.Tth) / P_REF_THERMO) - lnVh;
    for (let it = 1; it <= maxIt; it++) {
      let sumN = 0;
      for (let jj = 0; jj < nA; jj++) {
        const v = lnN[jj];
        const w = Math.exp(v);
        nj[jj] = w;
        sumN += w;
        mu[jj] = g[act[jj]] + v + lnC;
      }
      const lnSum = Math.log(sumN);
      for (let i = 0; i < dim; i++) {
        r[i] = 0;
        for (let j = 0; j < dim; j++) M[i * S + j] = 0;
      }
      const bcur = this.r2;
      for (let ii = 0; ii < nE; ii++) bcur[ii] = 0;
      for (let jj = 0; jj < nA; jj++) {
        const w = nj[jj];
        if (w === 0) continue;
        const m = mu[jj];
        const base = jj * NE;
        for (let ii = 0; ii < nE; ii++) {
          const aw = a[base + ii] * w;
          if (aw === 0) continue;
          for (let kk = 0; kk <= ii; kk++) M[ii * S + kk] += aw * a[base + kk];
          bcur[ii] += aw;
          r[ii] += aw * m;
        }
      }
      let balBad = false;
      for (let ii = 0; ii < nE; ii++) {
        for (let kk = 0; kk < ii; kk++) M[kk * S + ii] = M[ii * S + kk];
        const res = b0[ii] - bcur[ii];
        r[ii] += res;
        if (Math.abs(res) > 1e-6 * b0[ii]) balBad = true;
      }
      if (!this.factorScaled(dim, -1, sumN)) return -it;
      this.solveScaled(dim, r, x);
      let l1 = 0;
      let lam2 = 1;
      for (let jj = 0; jj < nA; jj++) {
        const base = jj * NE;
        let s = -mu[jj];
        for (let ii = 0; ii < nE; ii++) s += a[base + ii] * x[ii];
        dln[jj] = s;
        if (!Number.isFinite(s)) return -it;
        if (s > 0) {
          const lx = lnN[jj] - lnSum;
          if (lx <= -SIZE) {
            if (s >= SIZE + LN_1EM4) {
              const l2 = Math.abs(LN_1EM4 - lx) / s;
              if (l2 < lam2) lam2 = l2;
            }
          } else if (s > l1) l1 = s;
        }
      }
      let lam = l1 > 2 ? 2 / l1 : 1;
      if (lam2 < lam) lam = lam2;
      let conv = lam === 1 && !balBad;
      if (conv) {
        const invSum = 1 / sumN;
        for (let jj = 0; jj < nA; jj++) {
          const dj = dln[jj];
          const ad = Math.abs(dj);
          // judge the absolute change on the larger of the old and new amounts
          if (ad > STEP_TOL && nj[jj] * invSum * ad * (dj > 0 ? Math.exp(dj) : 1) > ABS_TOL) {
            conv = false;
            break;
          }
        }
      }
      for (let jj = 0; jj < nA; jj++) {
        let v = lnN[jj] + lam * dln[jj];
        if (v < lnSum + LN_X_FLOOR) v = lnSum + LN_X_FLOOR;
        lnN[jj] = v;
      }
      if (conv) {
        for (let ii = 0; ii < nE; ii++) this.piC[ii] = x[ii];
        return it;
      }
    }
    return -maxIt;
  }

  // =========================================================================================
  // Derivatives (RP-1311 §2.5)
  // =========================================================================================

  /** Assemble Σ_j a_kj a_ij n_j (+ b column/row and a zero ln n diagonal if withN) at the
   *  current state; returns Σ n_j. Also refreshes nj. */
  private assembleG(withN: boolean): number {
    const nA = this.nAct;
    const nE = this.nEl;
    const dim = withN ? nE + 1 : nE;
    const a = this.a;
    const M = this.M;
    const nj = this.nj;
    for (let i = 0; i < dim; i++) for (let j = 0; j < dim; j++) M[i * S + j] = 0;
    let sumN = 0;
    for (let jj = 0; jj < nA; jj++) {
      const w = Math.exp(this.lnN[jj]);
      nj[jj] = w;
      sumN += w;
      if (w === 0) continue;
      const base = jj * NE;
      for (let ii = 0; ii < nE; ii++) {
        const aw = a[base + ii] * w;
        if (aw === 0) continue;
        for (let kk = 0; kk <= ii; kk++) M[ii * S + kk] += aw * a[base + kk];
        if (withN) M[ii * S + nE] += aw;
      }
    }
    for (let ii = 0; ii < nE; ii++) {
      for (let kk = 0; kk < ii; kk++) M[kk * S + ii] = M[ii * S + kk];
      if (withN) M[nE * S + ii] = M[ii * S + nE];
    }
    return sumN;
  }

  /**
   * (∂/∂ln T)_p and (∂/∂ln p)_T of the current state (RP-1311 Tables 2.3 and 2.4):
   *   Σ_i G_ki ∂π_i + b_k ∂ln n = −Σ_j a_kj n_j H_j/RT   |  b_k
   *   Σ_i b_i ∂π_i            = −Σ_j n_j H_j/RT          |  Σ_j n_j
   *   ∂ln n_j/∂ln T = H_j/RT + Σ_i a_ij ∂π_i + ∂ln n ;  ∂ln n_j/∂ln p = −1 + Σ_i a_ij ∂π_i + ∂ln n
   * and cp,eq/R = Σ n_j cp_j/R + Σ n_j (H_j/RT)(∂ln n_j/∂ln T)_p   (≡ eq. 2.59).
   * Uses the cached thermo (must be at the state's T). Returns false if singular.
   */
  private derivTP(): boolean {
    const nA = this.nAct;
    const nE = this.nEl;
    const dim = nE + 1;
    const a = this.a;
    const act = this.act;
    const hRT = this.hRT;
    const sumN = this.assembleG(true);
    const rT = this.r;
    const rP = this.r2;
    for (let ii = 0; ii <= nE; ii++) {
      rT[ii] = 0;
      rP[ii] = 0;
    }
    let sh = 0;
    for (let jj = 0; jj < nA; jj++) {
      const w = this.nj[jj];
      if (w === 0) continue;
      const wh = w * hRT[act[jj]];
      sh += wh;
      const base = jj * NE;
      for (let ii = 0; ii < nE; ii++) {
        const aij = a[base + ii];
        rT[ii] -= aij * wh;
        rP[ii] += aij * w;
      }
    }
    rT[nE] = -sh;
    rP[nE] = sumN;
    if (!this.factorScaled(dim, nE, sumN)) return false;
    const xT = this.x;
    const xP = this.x2;
    this.solveScaled(dim, rT, xT);
    this.solveScaled(dim, rP, xP);
    const dnT = xT[nE];
    const dnP = xP[nE];
    for (let ii = 0; ii < nE; ii++) this.dpiT[ii] = xT[ii];
    let cp = 0;
    for (let jj = 0; jj < nA; jj++) {
      const k = act[jj];
      const base = jj * NE;
      let sT = dnT + hRT[k];
      let sP = dnP - 1;
      for (let ii = 0; ii < nE; ii++) {
        sT += a[base + ii] * xT[ii];
        sP += a[base + ii] * xP[ii];
      }
      this.dlnT[jj] = sT;
      this.dlnP[jj] = sP;
      const w = this.nj[jj];
      cp += w * (this.cpR[k] + hRT[k] * sT);
    }
    const der = this.der;
    der.dnT = dnT;
    der.dnP = dnP;
    der.cpR = cp;
    // cv from the identity cv = cp + n (∂lnV/∂lnT)²/(∂lnV/∂lnp) (in units of R)
    const vT = 1 + dnT;
    der.cvR = cp + (sumN * vT * vT) / (dnP - 1);
    return Number.isFinite(cp);
  }

  /**
   * (∂/∂ln T)_V of the current (TV) state: Σ_i G_ki ∂π_i = −Σ_j a_kj n_j U_j/RT,
   * ∂ln n_j/∂ln T = U_j/RT + Σ_i a_ij ∂π_i, cv,eq/R = Σ n_j cv_j/R + Σ n_j (U_j/RT) ∂ln n_j/∂ln T.
   * (Constant-volume analogue of RP-1311 Table 2.3.) Returns false if singular.
   */
  private derivTV(): boolean {
    const nA = this.nAct;
    const nE = this.nEl;
    const a = this.a;
    const act = this.act;
    const hRT = this.hRT;
    const sumN = this.assembleG(false);
    const rT = this.r;
    for (let ii = 0; ii < nE; ii++) rT[ii] = 0;
    for (let jj = 0; jj < nA; jj++) {
      const w = this.nj[jj];
      if (w === 0) continue;
      const wu = w * (hRT[act[jj]] - 1);
      const base = jj * NE;
      for (let ii = 0; ii < nE; ii++) rT[ii] -= a[base + ii] * wu;
    }
    if (!this.factorScaled(nE, -1, sumN)) return false;
    const xT = this.x;
    this.solveScaled(nE, rT, xT);
    for (let ii = 0; ii < nE; ii++) this.dpiT[ii] = xT[ii];
    let cv = 0;
    for (let jj = 0; jj < nA; jj++) {
      const k = act[jj];
      const base = jj * NE;
      const uk = hRT[k] - 1;
      let sT = uk;
      for (let ii = 0; ii < nE; ii++) sT += a[base + ii] * xT[ii];
      this.dlnT[jj] = sT;
      cv += this.nj[jj] * (this.cpR[k] - 1 + uk * sT);
    }
    this.der.cvR = cv;
    return Number.isFinite(cv);
  }

  // =========================================================================================
  // Outer iteration on T (HP = 0, UV = 1, SP = 2)
  // =========================================================================================

  private solveOuter(kind: number, b: Float64Array, target: number, pv: number, Tguess: number): EqResult {
    const res = this.result;
    res.iterations = 0;
    res.outerIterations = 0;
    // A NaN/∞ target, pressure or volume (e.g. from a blown-up integrator step) must fail loudly;
    // without this check a NaN residual drove the bisection fallback to "converge" at 6000 K.
    if (!Number.isFinite(target) || !(pv > 0 && pv < Infinity) || !this.setup(b)) {
      return this.infeasible(Tguess, kind === 1 ? NaN : pv);
    }
    const s = this.scale;
    // normalised target in units of R (K for H, U; dimensionless for S)
    const tau = (target * s) / R_UNIVERSAL;
    const lnP = kind === 1 ? 0 : Math.log(pv / P_REF_THERMO);
    const lnVh = kind === 1 ? Math.log(pv * s) : 0;
    let T = Tguess > EQ_T_MIN && Tguess < EQ_T_MAX ? Tguess : 2000;
    let lo = EQ_T_MIN;
    let hi = EQ_T_MAX;
    let loSet = false;
    let hiSet = false;
    let ok = false;
    let innerOk = false;
    for (let outer = 1; outer <= MAX_OUTER; outer++) {
      res.outerIterations = outer;
      this.setThermo(T);
      innerOk = kind === 1 ? this.innerTV(lnVh) : this.innerTP(lnP);
      if (!innerOk) break;
      const dOk = kind === 1 ? this.derivTV() : this.derivTP();
      if (!dOk) break;
      // residual f(T) and slope f'(T), normalised, units of R
      let f = 0;
      let fp = 0;
      const nA = this.nAct;
      const act = this.act;
      if (kind === 0) {
        let sh = 0;
        for (let jj = 0; jj < nA; jj++) sh += this.nj[jj] * this.hRT[act[jj]];
        f = T * sh - tau;
        fp = this.der.cpR;
      } else if (kind === 1) {
        let su = 0;
        for (let jj = 0; jj < nA; jj++) su += this.nj[jj] * (this.hRT[act[jj]] - 1);
        f = T * su - tau;
        fp = this.der.cvR;
      } else {
        let ss = 0;
        let n = 0;
        for (let jj = 0; jj < nA; jj++) n += this.nj[jj];
        const lnn = Math.log(n);
        for (let jj = 0; jj < nA; jj++) {
          const w = this.nj[jj];
          if (w === 0) continue;
          ss += w * (this.s0R[act[jj]] - (this.lnN[jj] - lnn) - lnP);
        }
        f = ss - tau;
        fp = this.der.cpR / T;
      }
      if (!Number.isFinite(f)) break;
      if (f > 0) {
        hi = T;
        hiSet = true;
      } else {
        lo = T;
        loSet = true;
      }
      let Tn = T - f / fp;
      // Only a genuine Newton step may signal convergence; a fallback step (non-positive or
      // non-finite slope) bisects an evaluated bracket or else probes the unevaluated bound.
      const newton = fp > 0 && Number.isFinite(Tn);
      if (!newton) Tn = loSet && hiSet ? 0.5 * (lo + hi) : f > 0 ? lo : hi;
      // limit the step, keep it inside the bracket
      const up = T * Math.exp(MAX_DLNT);
      const dn = T * Math.exp(-MAX_DLNT);
      if (Tn > up) Tn = up;
      if (Tn < dn) Tn = dn;
      if (newton && Math.abs(Tn - T) <= T_TOL * T) {
        // Finalise at Tn with the first-order composition update: the remaining errors are
        // O(Δln T²) ≲ 1e-14 in ln n_j and in the energy/entropy target.
        this.predict(kind, Math.log(Tn / T));
        this.setThermo(Tn);
        ok = true;
        break;
      }
      if (!(Tn > lo && Tn < hi)) {
        if (!(loSet && hiSet)) {
          // target beyond a bound we have not evaluated yet: go to the bound
          Tn = f > 0 ? lo : hi;
          if (Tn === T) break;
        } else Tn = 0.5 * (lo + hi);
      }
      if (hi - lo <= T_TOL * lo && loSet && hiSet) {
        ok = true;
        break;
      }
      if (outer === MAX_OUTER) break; // keep the state consistent with the last solved T
      // first-order predictor for the next inner solve
      const dlt = Math.log(Tn / T);
      if (Math.abs(dlt) <= MAX_PREDICT_DLNT) this.predict(kind, dlt);
      T = Tn;
    }
    const p = kind === 1 ? (this.sumN() * R_UNIVERSAL * this.Tth) / (pv * s) : pv;
    this.writeResult(this.Tth, p, ok && innerOk);
    return res;
  }

  /**
   * First-order update of the state to ln T + dlt using the derivatives of the last derivTP
   * (kinds 0, 2: (∂/∂ln T)_p) or derivTV (kind 1: (∂/∂ln T)_V); also advances π. Element
   * balance is preserved to first order (Σ_j a_ij n_j ∂ln n_j/∂ln T = 0).
   */
  private predict(kind: number, dlt: number): void {
    const nA = this.nAct;
    if (kind === 1) {
      for (let jj = 0; jj < nA; jj++) this.lnN[jj] += this.dlnT[jj] * dlt;
    } else {
      const lnNt = this.lnNt + this.der.dnT * dlt;
      for (let jj = 0; jj < nA; jj++) {
        let v = this.lnN[jj] + this.dlnT[jj] * dlt;
        if (v > lnNt) v = lnNt;
        this.lnN[jj] = v;
      }
      this.lnNt = lnNt;
    }
    for (let ii = 0; ii < this.nEl; ii++) this.piC[ii] += this.dpiT[ii] * dlt;
  }

  // =========================================================================================
  // Output
  // =========================================================================================

  private writeResult(T: number, p: number, converged: boolean): void {
    const res = this.result;
    const N = res.N;
    const X = res.X;
    N.fill(0);
    X.fill(0);
    res.pi.fill(0);
    const inv = 1 / this.scale;
    let n = 0;
    let m = 0;
    let sh = 0;
    for (let jj = 0; jj < this.nAct; jj++) {
      const k = this.act[jj];
      const w = Math.exp(this.lnN[jj]);
      const v = w * inv;
      N[k] = v;
      n += v;
      m += v * MOLAR_MASS[k];
      sh += v * this.hRT[k];
    }
    const invn = 1 / n;
    const lnn = Math.log(n);
    const lnP = Math.log(p / P_REF_THERMO);
    const lnInv = Math.log(inv);
    let ss = 0;
    for (let jj = 0; jj < this.nAct; jj++) {
      const k = this.act[jj];
      X[k] = N[k] * invn;
      if (N[k] > 0) ss += N[k] * (this.s0R[k] - (this.lnN[jj] + lnInv - lnn) - lnP);
    }
    for (let ii = 0; ii < this.nEl; ii++) res.pi[this.elIdx[ii]] = this.piC[ii];
    res.T = T;
    res.p = p;
    res.nTotal = n;
    res.mass = m;
    res.V = (n * R_UNIVERSAL * T) / p;
    res.H = R_UNIVERSAL * T * sh;
    res.U = R_UNIVERSAL * T * (sh - n);
    res.S = R_UNIVERSAL * ss;
    res.converged = converged;
    if (!converged) this.hasState = false;
  }

  private infeasible(T: number, p: number): EqResult {
    const res = this.result;
    this.feasible = false; // properties()/speciesDerivatives() now report invalid (warm start kept)
    res.N.fill(0);
    res.X.fill(NaN);
    res.pi.fill(0);
    res.T = T;
    res.p = p;
    res.V = NaN;
    res.nTotal = NaN;
    res.mass = NaN;
    res.H = NaN;
    res.U = NaN;
    res.S = NaN;
    res.converged = false;
    return res;
  }
}

/** Equilibrium thermodynamic properties (output of equilibriumProperties / solver.properties). */
export interface EqProperties {
  /** False if the solve did not converge or the derivative system was singular. */
  valid: boolean;
  /** Temperature, K; pressure, Pa. */
  T: number;
  p: number;
  /** Equilibrium mole fractions, Float64Array(NS). */
  readonly X: Float64Array;
  /** Mean molar mass, kg/mol. */
  M: number;
  /** Density, kg/m³. */
  rho: number;
  /** Specific enthalpy / internal energy (absolute, NASA convention), J/kg. */
  h: number;
  u: number;
  /** Specific entropy, J/(kg K). */
  s: number;
  /** Equilibrium (composition-shifting) specific heats, J/(kg K). */
  cp: number;
  cv: number;
  /** Frozen-composition cp, J/(kg K). */
  cpFrozen: number;
  /** (∂ln V/∂ln T)_p at equilibrium (1 for a frozen ideal gas). */
  dlnV_dlnT: number;
  /** (∂ln V/∂ln p)_T at equilibrium (−1 for a frozen ideal gas). */
  dlnV_dlnp: number;
  /** cp/cv (equilibrium). */
  gamma: number;
  /** Isentropic exponent (∂ln p/∂ln ρ)_s at equilibrium. */
  gammaS: number;
  /** Equilibrium speed of sound √(γ_s p/ρ), m/s. */
  soundSpeed: number;
}

/** Allocate an EqProperties record. */
export function newEqProperties(): EqProperties {
  return {
    valid: false,
    T: 0,
    p: 0,
    X: new Float64Array(NS),
    M: 0,
    rho: 0,
    h: 0,
    u: 0,
    s: 0,
    cp: 0,
    cv: 0,
    cpFrozen: 0,
    dlnV_dlnT: 0,
    dlnV_dlnp: 0,
    gamma: 0,
    gammaS: 0,
    soundSpeed: 0,
  };
}
