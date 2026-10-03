/**
 * Algebraic closure of the closed-cycle cylinder model: recovers the common pressure p and the
 * zone temperatures from the conserved state.
 *
 * Conserved (integrated) closed-cycle state (see cycle-model.ts):
 *   U      total internal energy of the charge, J (absolute NASA energies incl. formation)
 *   S_u    extensive entropy of the unburned zone, J/K
 *   m_u    unburned mass, kg (frozen composition X_u: fresh charge + residual)
 *   m_b    burned mass, kg (element composition of the unburned charge, chemical equilibrium over
 *          the N_EQ product species at (T_b, p))
 *
 * Closure: given (U, S_u, m_u, m_b, V) solve for (p, T_b):
 *   T_u = T(s_u = S_u/m_u, p)                       frozen ideal-gas isentrope (thermo.temperatureFromS)
 *   F₁ = m_u v_u(T_u, p) + m_b v_b(T_b, p) − V = 0     volume constraint
 *   F₂ = m_u u_u(T_u)    + m_b u_b(T_b, p) − U = 0     energy constraint
 * by Newton's method with the analytic Jacobian (equilibrium derivatives from NASA RP-1311 §2.5
 * as implemented in equilibrium/solver.ts properties()):
 *   ∂V_u/∂p|_s = −V_u/(γ_u p)          ∂U_u/∂p|_s = V_u/γ_u         (du = T ds − p dv, frozen ideal gas)
 *   ∂V_b/∂T_b = V_b ε_T/T_b            ∂V_b/∂p = V_b ε_p/p          ε_T = (∂lnV/∂lnT)_p, ε_p = (∂lnV/∂lnp)_T
 *   ∂U_b/∂T_b = m_b (c_p,eq − p v_b ε_T/T_b)
 *   ∂U_b/∂p   = −m_b v_b (ε_T + ε_p)   (from (∂h/∂p)_T = v(1 − ε_T), DESIGN.md integration notes)
 * With S_u integrated as dS_u = s_u dm_u − Q̇_u/T_u dt the unburned zone is isentropic except for
 * its wall heat loss; with dU = −p dV − Q̇ + P_spark this is exactly the classical two-zone
 * energy equation set (e.g. Heywood 1988 §9.4 / §14.4 — UNVERIFIED section numbers) written in
 * conservative form, so energy is conserved to round-off by construction.
 *
 * The same Jacobian gives the time derivatives (ṗ, Ṫ_b) of the closure from the state rates
 * ({@link ZoneClosure.rates}) — used for the rapid-distortion term of the unburned-zone turbulence.
 *
 * Single-zone (no burned gas) and burned-only (no unburned gas) cases are separate, well-posed
 * 1-D / 2-D problems. Allocation-free after construction.
 */
import { R_UNIVERSAL } from '../core/constants';
import { NE, NS, SP } from '../core/species';
import { EquilibriumSolver, newEqProperties, type EqProperties } from '../equilibrium';
import { completeCombustionMoles } from '../thermo/fuels';
import {
  elementMoles,
  mixHMolar,
  mixSMass,
  mixStateMass,
  temperatureFromS,
  temperatureFromUMolar,
  type MixState,
} from '../thermo/mixture';
import { MOLAR_MASS, speciesCpR, speciesHRT } from '../thermo/thermo';

/** Species touched by the kinetic-NO swap ({@link noSwapDeltas}). */
const SWAP_SPECIES = [SP.NO, SP.N2, SP.O2, SP.CO2, SP.CO, SP.H2O, SP.H2] as const;

/**
 * Element-conserving change dn (mol, Float64Array(NS), only SWAP_SPECIES written; the caller zeroes
 * it once) that replaces the NO moles of a burned-gas mole vector N by nNO; returns d = nNO − N_NO.
 *  - d ≤ 0 (kinetic NO below equilibrium): the missing NO is N2 + O2 (½N2 + ½O2 per NO);
 *  - d > 0: N from N2; the O atoms from O2 first, the remainder from CO2 → CO and H2O → H2 IN
 *    PROPORTION to their amounts (the water–gas-shift partner re-equilibrates the oxygen deficit
 *    between the two; round 1 took it all from CO2, which added ≈ the NO amount of CO to a rich
 *    exhaust — validation round 2), each capped by what is available.
 */
export function noSwapDeltas(N: Float64Array, nNO: number, dn: Float64Array): number {
  for (const k of SWAP_SPECIES) dn[k] = 0;
  const d = nNO - N[SP.NO];
  dn[SP.NO] = d;
  dn[SP.N2] = -0.5 * d;
  if (d <= 0) {
    dn[SP.O2] = -0.5 * d;
    return d;
  }
  let o = d; // O atoms needed
  const fromO2 = Math.min(0.5 * o, N[SP.O2]);
  dn[SP.O2] = -fromO2;
  o -= 2 * fromO2;
  if (o > 0) {
    const c = N[SP.CO2];
    const w = N[SP.H2O];
    const tot = c + w;
    let fc = tot > 0 ? Math.min(o * (c / tot), c) : 0;
    let fw = tot > 0 ? Math.min(o - fc, w) : 0;
    if (fc + fw < o) fc = Math.min(c, o - fw); // (one of them exhausted)
    dn[SP.CO2] = -fc;
    dn[SP.CO] = fc;
    dn[SP.H2O] = -fw;
    dn[SP.H2] = fw;
  }
  return d;
}

/** Newton tolerance on the relative updates of p and T_b. */
const NEWTON_RTOL = 1e-10;
const NEWTON_MAX_ITER = 40;
/**
 * Relative residual regarded as round-off: volume vs V, energy vs Σ|U_i| (absolute NASA energies
 * of a few hundred J carry ≈ 1e-13 J of round-off, which a 1e-10 kg burned zone turns into
 * δT_b ≈ 1e-5 K — no iteration can go below that).
 */
const ROUNDOFF_RES = 256 * Number.EPSILON;
/**
 * Relative round-off floor of the ENERGY residual, referred to the zones' thermal energy scale
 * Σ m_i c_v,i T_i (+ |U|): u_u(T_u(s_u, p)) is a sum of species energies with formation terms that
 * cancel, and the isentrope inversion adds its own noise — measured ≈ 3e-13 … 8e-13 of that scale
 * (validation round 2: a 3e-10 kg kernel burned zone "failed" with F₂ at 3e-13 J while converged, and
 * T_b oscillated by 1e-6 K). 4096 ε ≈ 9e-13.
 */
const ENERGY_NOISE = 4096 * Number.EPSILON;

export class ZoneClosure {
  /** Unburned mole fractions (frozen), Float64Array(NS). */
  readonly Xu = new Float64Array(NS);
  /** Unburned moles per kg, mol/kg. */
  readonly nu = new Float64Array(NS);
  /** Element moles per kg of charge (unburned = burned), mol/kg. */
  readonly bu = new Float64Array(NE);
  /** Complete-combustion product moles per kg of charge (heat-release diagnostic), mol/kg. */
  readonly ncc = new Float64Array(NS);
  /** Unburned molar mass (kg/mol) and gas constant (J/(kg K)). */
  Mu = 0;
  Ru = 0;

  /** Burned-zone equilibrium solver (warm-started; owns its result). */
  readonly eq = new EquilibriumSolver();
  /** Equilibrium properties of the last burned-zone solve. */
  readonly props: EqProperties = newEqProperties();
  private readonly bb = new Float64Array(NE);
  private readonly mix: MixState = { cp: 0, cv: 0, h: 0, u: 0, gamma: 0, molarMass: 0 };

  // ---- results of the last solve (SI) ----
  p = 0;
  Tu = 0;
  Tb = 0;
  /** Unburned zone: specific entropy, internal energy, enthalpy, volume, cp, cv, γ (frozen). */
  su = 0;
  uu = 0;
  hu = 0;
  vu = 0;
  cpu = 0;
  cvu = 0;
  gu = 1.4;
  /** Burned zone (equilibrium): specific internal energy, enthalpy, volume, cp,eq, ε_T, ε_p, molar mass, γ_s, total moles. */
  ub = 0;
  hb = 0;
  vb = 0;
  cpb = 0;
  epsT = 1;
  epsP = -1;
  Mb = 0;
  gammaSb = 1.25;
  gammab = 1.25;
  /** Equilibrium c_v of the burned zone, J/(kg K). */
  cvb = 0;
  nb = 0;
  /** Jacobian of (F₁, F₂) w.r.t. (p, T_b) at the last iterate. */
  j11 = 0;
  j12 = 0;
  j21 = 0;
  j22 = 0;
  /** Time derivatives from {@link rates}. */
  dp = 0;
  dTb = 0;
  /**
   * Rate-controlled NO in the burned-zone thermodynamics (validation round 2): with noCoupled the
   * burned zone is the TP equilibrium with its NO replaced by the kinetic amount noKinetic (mol, whole
   * burned zone; noSwapDeltas), i.e. its internal energy and volume carry
   *   ΔU = Σ Δn_k u_k(T_b),  ΔV = Σ Δn_k R T_b/p
   * — the energy of the NO departure from equilibrium is stored in (formation) or released from (not
   * yet decomposed) the gas continuously, so the charge temperature/pressure do not jump when the
   * kinetic NO is written into the gas at EVO (round 1: −0.8…−1.2 % in p, −11…−15 K). Off: full
   * equilibrium (NO a passive tracer, the fuel–air-cycle and Cantera two-zone oracle definitions).
   */
  noCoupled = false;
  noKinetic = 0;
  /** Swap corrections of the last burned-zone evaluation and their (p, T_b) derivatives. */
  dUno = 0;
  dVno = 0;
  private dUnoT = 0;
  private dUnoP = 0;
  private dVnoT = 0;
  private dVnoP = 0;
  private readonly dnSwap = new Float64Array(NS);
  private readonly dlnNO = new Float64Array(2);

  /** Newton iterations and equilibrium solves of the last call; cumulative solve counter. */
  iterations = 0;
  converged = true;
  eqSolves = 0;

  /** Set the frozen unburned composition (mole fractions, normalised on copy). */
  setUnburned(X: Float64Array): void {
    let s = 0;
    for (let k = 0; k < NS; k++) s += X[k] > 0 ? X[k] : 0;
    let M = 0;
    for (let k = 0; k < NS; k++) {
      const x = X[k] > 0 ? X[k] / s : 0;
      this.Xu[k] = x;
      M += x * MOLAR_MASS[k];
    }
    this.Mu = M;
    this.Ru = R_UNIVERSAL / M;
    for (let k = 0; k < NS; k++) this.nu[k] = this.Xu[k] / M;
    elementMoles(this.nu, this.bu);
    completeCombustionMoles(this.nu, this.ncc);
    this.eq.reset();
  }

  /** Mass-specific unburned properties at temperature T and pressure p (sets uu, hu, vu, cp, cv, γ). */
  unburnedAtT(T: number, p: number): void {
    mixStateMass(this.Xu, T, this.mix);
    this.Tu = T;
    this.uu = this.mix.u;
    this.hu = this.mix.h;
    this.cpu = this.mix.cp;
    this.cvu = this.mix.cv;
    this.gu = this.mix.gamma;
    this.vu = (this.Ru * T) / p;
  }

  /** Unburned specific entropy at (T, p), J/(kg K). */
  entropyAt(T: number, p: number): number {
    return mixSMass(this.Xu, T, p);
  }

  /**
   * Chemical heat of reaction of the unburned mixture at T, J/kg: h_u(T) − h_cc(T) with the
   * complete-combustion products of thermo/fuels.completeCombustionMoles (lean: CO2/H2O/O2,
   * rich: water-gas shift frozen at 1740 K). Diagnostic (heat-release rate, knock source).
   */
  heatOfReaction(T: number): number {
    return mixHMolar(this.nu, T) - mixHMolar(this.ncc, T);
  }

  private burnedAt(mb: number, Tb: number, p: number): boolean {
    const bb = this.bb;
    for (let e = 0; e < NE; e++) bb[e] = this.bu[e] * mb;
    let r = this.eq.solveTP(bb, Tb, p);
    this.eqSolves++;
    if (!r.converged) {
      // cold restart once (a far warm start can stall the damped Newton)
      this.eq.reset();
      r = this.eq.solveTP(bb, Tb, p);
      this.eqSolves++;
    }
    const pr = this.eq.properties(this.props);
    this.Tb = Tb;
    this.ub = pr.u;
    this.hb = pr.h;
    this.vb = 1 / pr.rho;
    this.cpb = pr.cp;
    this.epsT = pr.dlnV_dlnT;
    this.epsP = pr.dlnV_dlnp;
    this.Mb = pr.M;
    this.gammaSb = pr.gammaS;
    this.gammab = pr.gamma;
    this.cvb = pr.cv;
    this.nb = r.nTotal;
    this.noCorrection(r.N, Tb, p);
    return r.converged && pr.valid;
  }

  /** ΔU, ΔV of the kinetic-NO swap at (T_b, p) and their derivatives (see noCoupled). */
  private noCorrection(N: Float64Array, Tb: number, p: number): void {
    if (!this.noCoupled) {
      this.dUno = 0;
      this.dVno = 0;
      this.dUnoT = 0;
      this.dUnoP = 0;
      this.dVnoT = 0;
      this.dVnoP = 0;
      return;
    }
    const dn = this.dnSwap;
    const d = noSwapDeltas(N, this.noKinetic, dn);
    const R = R_UNIVERSAL;
    let dU = 0;
    let dCv = 0;
    let dN = 0;
    for (const k of SWAP_SPECIES) {
      const x = dn[k];
      if (x === 0) continue;
      dU += x * R * Tb * (speciesHRT(k, Tb) - 1);
      dCv += x * R * (speciesCpR(k, Tb) - 1);
      dN += x;
    }
    this.dUno = dU;
    this.dVno = (dN * R * Tb) / p;
    // d = nNO − N_NO,eq(T_b, p): its derivatives from the equilibrium's (∂ln N_NO/∂ln T)_p, (∂ln N/∂ln p)_T
    this.eq.lastSpeciesLogDerivatives(SP.NO, this.dlnNO);
    const nEq = N[SP.NO];
    const ddT = (-nEq * this.dlnNO[0]) / Tb;
    const ddP = (-nEq * this.dlnNO[1]) / p;
    const perD = d !== 0 ? dU / d : 0; // Σ ν_k u_k of the swap stoichiometry
    const nPerD = d !== 0 ? dN / d : 0;
    this.dUnoT = dCv + ddT * perD;
    this.dUnoP = ddP * perD;
    this.dVnoT = (R / p) * (dN + Tb * ddT * nPerD);
    this.dVnoP = ((R * Tb) / p) * (ddP * nPerD - dN / p);
  }

  /**
   * Single zone (no burned gas): T from U/m (frozen composition), p = m R T/V.
   * @returns T (K)
   */
  solveSingle(U: number, m: number, V: number, Tguess: number): number {
    const T = temperatureFromUMolar(this.nu, U / m, Tguess > 0 ? Tguess : 1000);
    const p = (m * this.Ru * T) / V;
    this.p = p;
    this.unburnedAtT(T, p);
    this.Tb = 0;
    this.iterations = 1;
    this.converged = true;
    return T;
  }

  /**
   * Two zones (m_u > 0, m_b > 0): Newton on (p, T_b). Returns false if it did not converge
   * (the best iterate is kept).
   */
  solveTwoZone(U: number, S: number, mu: number, mb: number, V: number, pGuess: number, TbGuess: number, TuGuess: number): boolean {
    const su = S / mu;
    this.su = su;
    let p = pGuess;
    let Tb = TbGuess;
    let Tu = TuGuess;
    this.converged = false;
    let it = 0;
    // Pressure-only (volume-constraint) Newton steps first while the volume residual is large:
    // with a small burned mass the energy residual of a poor p guess would otherwise throw T_b
    // far off (U_u(p) is curved; the error lands entirely on m_b u_b).
    let coupled = false;
    let f2Prev = Infinity;
    for (; it < NEWTON_MAX_ITER; it++) {
      Tu = temperatureFromS(this.Xu, su, p, Tu);
      this.unburnedAtT(Tu, p);
      const okB = this.burnedAt(mb, Tb, p);
      const Vu = mu * this.vu;
      const Vb = mb * this.vb;
      const F1 = Vu + Vb + this.dVno - V;
      const F2 = mu * this.uu + mb * this.ub + this.dUno - U;
      const j11 = -Vu / (this.gu * p) + (Vb * this.epsP) / p + this.dVnoP;
      const j12 = (Vb * this.epsT) / Tb + this.dVnoT;
      const j21 = Vu / this.gu - Vb * (this.epsT + this.epsP) + this.dUnoP;
      const j22 = mb * (this.cpb - (p * this.vb * this.epsT) / Tb) + this.dUnoT;
      this.j11 = j11;
      this.j12 = j12;
      this.j21 = j21;
      this.j22 = j22;
      if (!coupled && Math.abs(F1) > 1e-7 * V) {
        let dp = -F1 / j11;
        if (!Number.isFinite(dp)) break;
        if (Math.abs(dp) > 0.5 * p) dp = 0.5 * p * Math.sign(dp);
        p += dp;
        continue;
      }
      coupled = true;
      // residuals at the round-off floor: the current iterate IS the solution (with a tiny burned
      // mass the T_b step tolerance below can sit under the energy round-off; round 1 counted such
      // solves as failures)
      const sumU = Math.abs(U) + Math.abs(mu * this.uu) + Math.abs(mb * this.ub);
      if (okB && Math.abs(F1) <= ROUNDOFF_RES * V && Math.abs(F2) <= ROUNDOFF_RES * sumU) {
        this.converged = true;
        it++;
        break;
      }
      // noise floor reached: the energy residual is within the round-off of the zones' thermal energy
      // scale and has stopped decreasing (validation round 2: a converged 3e-10 kg burned zone kept
      // oscillating at |F₂| ≈ 3e-13 J, above the Σ|U| floor, and was counted as a failed solve)
      const eScale = sumU + mu * this.cvu * Tu + mb * this.cvb * Tb;
      const f2 = Math.abs(F2);
      if (okB && Math.abs(F1) <= 16 * ROUNDOFF_RES * V && f2 <= ENERGY_NOISE * eScale && f2 >= 0.5 * f2Prev) {
        this.converged = true;
        it++;
        break;
      }
      f2Prev = f2;
      const det = j11 * j22 - j12 * j21;
      let dp = -(F1 * j22 - j12 * F2) / det;
      let dT = -(j11 * F2 - j21 * F1) / det;
      if (!(Number.isFinite(dp) && Number.isFinite(dT))) break;
      // damping: keep p > 0 and T_b in the equilibrium range
      const lim = Math.max(Math.abs(dp) / (0.5 * p), Math.abs(dT) / (0.3 * Tb), 1);
      dp /= lim;
      dT /= lim;
      p += dp;
      Tb += dT;
      if (Tb < 250) Tb = 250;
      if (Tb > 5900) Tb = 5900;
      // T_b tolerance: relative, or the round-off floor of the energy balance for a tiny burned
      // mass (δU ≈ ε·Σ|U_i| maps to δT_b = δU/∂U_b/∂T_b)
      const tolT = Math.max(NEWTON_RTOL * Tb, (ROUNDOFF_RES * sumU) / Math.abs(j22));
      if (okB && lim === 1 && Math.abs(dp) <= NEWTON_RTOL * p && Math.abs(dT) <= tolT) {
        this.converged = true;
        it++;
        break;
      }
    }
    this.iterations = it;
    this.p = p;
    this.Tb = Tb;
    return this.converged;
  }

  /** Burned zone only (m_u = 0): Newton on (p, T_b) for U = m_b u_b, V = m_b v_b (UV equilibrium). */
  solveBurnedOnly(U: number, mb: number, V: number, pGuess: number, TbGuess: number): boolean {
    let p = pGuess;
    let Tb = TbGuess;
    this.converged = false;
    let it = 0;
    for (; it < NEWTON_MAX_ITER; it++) {
      const okB = this.burnedAt(mb, Tb, p);
      const Vb = mb * this.vb;
      const F1 = Vb + this.dVno - V;
      const F2 = mb * this.ub + this.dUno - U;
      const j11 = (Vb * this.epsP) / p + this.dVnoP;
      const j12 = (Vb * this.epsT) / Tb + this.dVnoT;
      const j21 = -Vb * (this.epsT + this.epsP) + this.dUnoP;
      const j22 = mb * (this.cpb - (p * this.vb * this.epsT) / Tb) + this.dUnoT;
      this.j11 = j11;
      this.j12 = j12;
      this.j21 = j21;
      this.j22 = j22;
      if (okB && Math.abs(F1) <= ROUNDOFF_RES * V && Math.abs(F2) <= ENERGY_NOISE * (Math.abs(U) + Math.abs(mb * this.ub) + mb * this.cvb * Tb)) {
        this.converged = true;
        it++;
        break;
      }
      const det = j11 * j22 - j12 * j21;
      let dp = -(F1 * j22 - j12 * F2) / det;
      let dT = -(j11 * F2 - j21 * F1) / det;
      if (!(Number.isFinite(dp) && Number.isFinite(dT))) break;
      const lim = Math.max(Math.abs(dp) / (0.5 * p), Math.abs(dT) / (0.3 * Tb), 1);
      dp /= lim;
      dT /= lim;
      p += dp;
      Tb += dT;
      if (Tb < 250) Tb = 250;
      if (Tb > 5900) Tb = 5900;
      if (okB && lim === 1 && Math.abs(dp) <= NEWTON_RTOL * p && Math.abs(dT) <= NEWTON_RTOL * Tb) {
        this.converged = true;
        it++;
        break;
      }
    }
    this.iterations = it;
    this.p = p;
    this.Tb = Tb;
    this.Tu = 0;
    return this.converged;
  }

  /**
   * Time derivatives (ṗ, Ṫ_b) of the two-zone / burned-only closure from the state rates
   * (linearised constraints, same Jacobian): writes this.dp, this.dTb.
   * @param dU dU/dt (W); @param dSnet dS_u/dt − s_u dm_u/dt = −Q̇_u/T_u (W/K);
   * @param dmu dm_u/dt; @param dmb dm_b/dt (kg/s); @param dV dV/dt (m³/s)
   * @param mu unburned mass (0 in burned-only mode)
   */
  rates(dU: number, dSnet: number, dmu: number, dmb: number, dV: number, mu: number): void {
    let r1 = dV - this.vb * dmb;
    let r2 = dU - this.ub * dmb;
    if (mu > 0) {
      r1 -= this.vu * dmu + (this.vu / this.cpu) * dSnet;
      r2 -= this.uu * dmu + (this.Tu / this.gu) * dSnet;
    }
    const det = this.j11 * this.j22 - this.j12 * this.j21;
    this.dp = (r1 * this.j22 - this.j12 * r2) / det;
    this.dTb = (this.j11 * r2 - this.j21 * r1) / det;
  }

  /**
   * Adiabatic constant-pressure (HP) equilibrium flame temperature of the unburned charge at
   * (T_u, p) using `solver` (a separate instance so the burned-zone warm start is kept), K.
   */
  flameTemperatureHP(solver: EquilibriumSolver, Tu: number, p: number): number {
    const h = mixHMolar(this.nu, Tu); // J/kg (nu in mol/kg)
    const r = solver.solveHP(this.bu, h, p, 2300);
    return r.converged ? r.T : 2300;
  }

  /** Moles of NO per kg of unburned charge, mol/kg. */
  get noPerKg(): number {
    return this.nu[SP.NO];
  }
}
