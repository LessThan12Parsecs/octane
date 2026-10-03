/**
 * Dilute-gas transport properties (mixture-averaged), SI units.
 *
 * Pure species (Chapman–Enskog, first approximation; Hirschfelder, Curtiss & Bird 1954):
 *   μ_k = (5/16) √(π m_k k_B T) / (π σ_k² Ω(2,2)*(T*))                          [Pa s]
 *   ρD_kk/μ_k = (6/5) A*,  A* = Ω(2,2)* / Ω(1,1)*  (self-diffusion, exact identity)
 * with Lennard-Jones parameters from species-data.ts (GRI-Mech 3.0 / LLNL) and collision
 * integrals from collision-integrals.ts (Neufeld et al. 1972; Monchick–Mason 1961 polar
 * correction for H2O).
 *
 * Pure-species thermal conductivity: the Warnatz / Kee rotational-relaxation model
 * (Kee, Coltrin & Glarborg 2003, "Chemically Reacting Flow", ch. 12; the model of CHEMKIN
 * TRANFIT and Cantera — formula cross-checked line by line against Cantera 3.2.0
 * src/transport/GasTransport.cpp::fitProperties):
 *   λ = (μ/M) R (f_tr c_v,tr + f_rot c_v,rot + f_vib c_v,vib)        (c_v in units of R)
 *   f_tr = 5/2 (1 − (2/π)(c_v,rot/c_v,tr)(A/B)),  f_rot = (ρD/μ)(1 + (2/π) A/B),
 *   f_vib = ρD/μ,  A = 5/2 − ρD/μ,  B = Z_rot(T) + (2/π)(5/3 c_v,rot + ρD/μ),
 *   Z_rot(T) = Z_rot(298) F(298)/F(T),  Parker (1959) F(T) = 1 + (π^{3/2}/2)(ε/kT)^{1/2}
 *             + (π²/4 + 2)(ε/kT) + π^{3/2}(ε/kT)^{3/2}   (Kee et al. 2003 eq. 12.112)
 *   c_v,tr = 3/2, c_v,rot = 0 | 1 | 3/2 (atom | linear | nonlinear),
 *   c_v,vib = c_p/R − 5/2 − c_v,rot from the NASA polynomials.
 * UNVERIFIED: the Kee et al. equation numbers other than 12.112 (quoted by Cantera).
 * Design note: this model is used instead of the modified Eucken rule because it is
 * physically more complete (rotational relaxation, diffusion of internal energy) and tracks
 * the Cantera oracle to ≤ 0.6 % for every species, whereas modified Eucken deviates by up to
 * 10 % (OH, H2O). Modified Eucken is kept as speciesThermalConductivityModifiedEucken.
 *
 * WATER is the exception. For the strongly polar H2O molecule the model above (as in
 * Cantera/CHEMKIN with the GRI-Mech parameters) over-predicts the steam conductivity by
 * 20–40 % at 400–1200 K against the IAPWS 2011 standard (resonant exchange of rotational
 * energy between polar molecules hinders internal-energy diffusion; Mason & Monchick 1962),
 * and the viscosity is off by −4…+5 %. H2O therefore uses the IAPWS dilute-gas
 * correlations (waterViscosityIAPWS / waterThermalConductivityIAPWS below) — reference
 * data, not a model. The Chapman–Enskog/Warnatz values remain available as
 * speciesViscosityChapmanEnskog / speciesThermalConductivityWarnatz.
 * Other known model deviations vs reference correlations (NIST REFPROP via WebBook,
 * transport.test.ts): N2 μ ≤ 1.1 %, λ +4…+9 % at 900–1800 K; CO2 μ ≤ 0.6 %, λ −3…+4 %.
 *
 * Mixture rules: Wilke (1950), J. Chem. Phys. 18:517, for viscosity; Mathur, Tondon & Saxena
 * (1967), Mol. Phys. 12:569, for conductivity: λ = ½ (Σ X_k λ_k + 1/Σ X_k/λ_k).
 *
 * Hot paths allocate nothing (module-level scratch; not re-entrant — single worker thread).
 */
import { K_BOLTZMANN, N_AVOGADRO, R_UNIVERSAL } from '../core/constants';
import { NS, SP } from '../core/species';
import { buildPolarCorrection, interpTstar, omega11LJ, omega22LJ, type PolarCorrection } from './collision-integrals';
import { mixCpMass, mixMolarMass } from './mixture';
import { SPECIES_TRANSPORT } from './species-data';
import { MOLAR_MASS, speciesCpR } from './thermo';

/** Vacuum permittivity, F/m (CODATA 2018). */
const EPSILON_0 = 8.8541878128e-12;
/** 1 Debye in C m (1e-21 / c). */
const DEBYE = 1e-21 / 299792458;

const PI = Math.PI;
const PI15 = Math.pow(PI, 1.5);
const K_H2O = SP.H2O;

/** Parker (1959) temperature scaling of Z_rot, as in Kee et al. (2003) eq. 12.112. */
function parkerF(Ts: number): number {
  const r = 1 / Ts;
  return 1 + 0.5 * PI15 * Math.sqrt(r) + (0.25 * PI * PI + 2) * r + PI15 * r * Math.sqrt(r);
}

// ---- per-species constants ------------------------------------------------------------
const EPS_K = new Float64Array(NS); // ε/k_B, K
const VISC_PREF = new Float64Array(NS); // μ = VISC_PREF √T / Ω22
const CV_ROT = new Float64Array(NS); // c_v,rot / R
const ZROT_F298 = new Float64Array(NS); // Z_rot(298) F(298)
const POLAR: (PolarCorrection | null)[] = new Array(NS).fill(null);
/** Reduced dipole moment δ* of each species (0 for non-polar). */
export const REDUCED_DIPOLE = new Float64Array(NS);

for (let k = 0; k < NS; k++) {
  const t = SPECIES_TRANSPORT[k];
  const sigma = t.diameter * 1e-10;
  const m = MOLAR_MASS[k] / N_AVOGADRO;
  EPS_K[k] = t.wellDepth;
  VISC_PREF[k] = ((5 / 16) * Math.sqrt(PI * m * K_BOLTZMANN)) / (PI * sigma * sigma);
  CV_ROT[k] = t.geometry === 0 ? 0 : t.geometry === 1 ? 1 : 1.5;
  ZROT_F298[k] = t.rotRelax * parkerF(298 / t.wellDepth);
  if (t.dipole > 0) {
    const mu = t.dipole * DEBYE;
    const d = (0.5 * mu * mu) / (4 * PI * EPSILON_0 * t.wellDepth * K_BOLTZMANN * sigma * sigma * sigma);
    REDUCED_DIPOLE[k] = d;
    POLAR[k] = buildPolarCorrection(d);
  }
}

// Wilke pair constants: (M_j/M_k)^{1/4} and 1/√(8(1 + M_k/M_j)), row-major [k*NS + j].
const WILKE_M4 = new Float64Array(NS * NS);
const WILKE_DEN = new Float64Array(NS * NS);
for (let k = 0; k < NS; k++) {
  for (let j = 0; j < NS; j++) {
    WILKE_M4[k * NS + j] = Math.pow(MOLAR_MASS[j] / MOLAR_MASS[k], 0.25);
    WILKE_DEN[k * NS + j] = 1 / Math.sqrt(8 * (1 + MOLAR_MASS[k] / MOLAR_MASS[j]));
  }
}

// ---- Chapman–Enskog / Warnatz pure-species model -------------------------------------------

// Results of the last evalCollision call (module scratch; typed array — no boxed doubles):
// [0] Ω(2,2)*, [1] Ω(1,1)*
const COL = new Float64Array(2);

function evalCollision(k: number, T: number): void {
  const Ts = T / EPS_K[k];
  let o22 = omega22LJ(Ts);
  let o11 = omega11LJ(Ts);
  const pc = POLAR[k];
  if (pc !== null) {
    const r22 = interpTstar(pc.r22, Ts);
    const rA = interpTstar(pc.rA, Ts);
    o22 *= r22;
    o11 *= r22 / rA;
  }
  COL[0] = o22;
  COL[1] = o11;
}

/** Kee/Warnatz conductivity given μ and the current collision integrals. */
function keeConductivity(k: number, T: number, mu: number): number {
  const fInt = (1.2 * COL[0]) / COL[1]; // ρD_kk/μ_k
  const cvRot = CV_ROT[k];
  const cvVib = speciesCpR(k, T) - 2.5 - cvRot; // (c_v − c_v,tr − c_v,rot)/R
  let fTr = 2.5;
  let fRot = fInt;
  if (cvRot > 0) {
    const A = 2.5 - fInt;
    const B = ZROT_F298[k] / parkerF(T / EPS_K[k]) + (2 / PI) * ((5 / 3) * cvRot + fInt);
    const c1 = ((2 / PI) * A) / B;
    fTr = 2.5 * (1 - (c1 * cvRot) / 1.5);
    fRot = fInt * (1 + c1);
  }
  return ((mu / MOLAR_MASS[k]) * R_UNIVERSAL) * (fTr * 1.5 + fRot * cvRot + fInt * cvVib);
}

/** Chapman–Enskog (Lennard-Jones / Stockmayer) viscosity of species k at T (K), Pa s. */
export function speciesViscosityChapmanEnskog(k: number, T: number): number {
  evalCollision(k, T);
  return (VISC_PREF[k] * Math.sqrt(T)) / COL[0];
}

/**
 * Warnatz/Kee rotational-relaxation conductivity of species k at T (K) with the
 * Chapman–Enskog viscosity, W/(m K) — the Cantera/CHEMKIN pure-species model (see header).
 */
export function speciesThermalConductivityWarnatz(k: number, T: number): number {
  const mu = speciesViscosityChapmanEnskog(k, T);
  return keeConductivity(k, T, mu);
}

// ---- H2O: IAPWS dilute-gas correlations --------------------------------------------------

/**
 * Range in which the IAPWS dilute-gas formulas are used directly, K. IAPWS R12-08 (2008)
 * §2.4: the dilute-gas viscosity "behaves in a physically reasonable manner down to at
 * least 250 K" and its extrapolation "is physically reasonable up to at least 2500 K";
 * R15-11 (2011) §2.4: λ0 is reasonable down to 250 K and extrapolates reasonably above its
 * 1173.15 K validity limit (IAPWS warns that above ~1500 K the MEASURABLE conductivity
 * gains a dissociation (reactive) contribution, which λ0 excludes — the frozen conductivity
 * we need here). Outside [250 K, 2500 K] the Chapman–Enskog/Warnatz value is scaled by its
 * ratio to IAPWS at the nearer limit, so the property stays continuous and keeps the
 * kinetic-theory temperature dependence (the μ0 polynomial itself turns negative near 133 K).
 */
export const IAPWS_T_MIN = 250;
/** Upper limit of the direct use of the IAPWS dilute-gas formulas, K (see IAPWS_T_MIN). */
export const IAPWS_T_MAX = 2500;

/** IAPWS reference temperature T* = T_c, K (R12-08 eq. 1, R15-11 eq. 1). */
const T_STAR_W = 647.096;

/**
 * IAPWS 2008 dilute-gas viscosity of water vapour (IAPWS R12-08, eq. 11, Table 1):
 *   μ0/μ* = 100 √T̄ / Σ_{i=0..3} H_i / T̄^i,  T̄ = T/647.096 K,  μ* = 1e-6 Pa s,
 *   H = 1.67752, 2.20462, 0.6366564, −0.241605.
 * Check: agrees with the NIST WebBook (full IAPWS formulation at 0.01 MPa) to ≤ 6e-5 at
 * 600–1200 K. Pa s. Raw formula — use speciesViscosity(SP.H2O, T) outside [250, 2500] K.
 */
export function waterViscosityIAPWS(T: number): number {
  const t = T / T_STAR_W;
  const r = 1 / t;
  return (1e-4 * Math.sqrt(t)) / (1.67752 + r * (2.20462 + r * (0.6366564 - r * 0.241605)));
}

/**
 * IAPWS 2011 dilute-gas thermal conductivity of water vapour (IAPWS R15-11, eq. 16,
 * Table 1): λ0/λ* = √T̄ / Σ_{k=0..4} L_k / T̄^k,  λ* = 1e-3 W/(m K),
 *   L = 2.443221e-3, 1.323095e-2, 6.770357e-3, −3.454586e-3, 4.096266e-4.
 * Check values (R15-11 Tables 4–5, ρ = 0): 18.4341883 mW/(m K) at 298.15 K, 51.5764797 at
 * 647.35 K, 79.1034659 at 873.15 K. W/(m K). Raw formula (see waterViscosityIAPWS).
 */
export function waterThermalConductivityIAPWS(T: number): number {
  const t = T / T_STAR_W;
  const r = 1 / t;
  return (
    (1e-3 * Math.sqrt(t)) /
    (2.443221e-3 + r * (1.323095e-2 + r * (6.770357e-3 + r * (-3.454586e-3 + r * 4.096266e-4))))
  );
}

// Continuity factors at the IAPWS range limits (IAPWS / kinetic-theory model).
const H2O_MU_LO = waterViscosityIAPWS(IAPWS_T_MIN) / speciesViscosityChapmanEnskog(K_H2O, IAPWS_T_MIN);
const H2O_MU_HI = waterViscosityIAPWS(IAPWS_T_MAX) / speciesViscosityChapmanEnskog(K_H2O, IAPWS_T_MAX);
const H2O_LAM_LO = waterThermalConductivityIAPWS(IAPWS_T_MIN) / speciesThermalConductivityWarnatz(K_H2O, IAPWS_T_MIN);
const H2O_LAM_HI = waterThermalConductivityIAPWS(IAPWS_T_MAX) / speciesThermalConductivityWarnatz(K_H2O, IAPWS_T_MAX);

function waterViscosity(T: number): number {
  if (T < IAPWS_T_MIN) return H2O_MU_LO * speciesViscosityChapmanEnskog(K_H2O, T);
  if (T > IAPWS_T_MAX) return H2O_MU_HI * speciesViscosityChapmanEnskog(K_H2O, T);
  return waterViscosityIAPWS(T);
}

function waterConductivity(T: number): number {
  if (T < IAPWS_T_MIN) return H2O_LAM_LO * speciesThermalConductivityWarnatz(K_H2O, T);
  if (T > IAPWS_T_MAX) return H2O_LAM_HI * speciesThermalConductivityWarnatz(K_H2O, T);
  return waterThermalConductivityIAPWS(T);
}

// ---- pure species (the values the mixture rules use) ------------------------------------

/**
 * Pure-species dynamic viscosity of species k at T (K), Pa s: Chapman–Enskog, except H2O
 * (IAPWS 2008 dilute gas; see header).
 */
export function speciesViscosity(k: number, T: number): number {
  return k === K_H2O ? waterViscosity(T) : speciesViscosityChapmanEnskog(k, T);
}

/**
 * Pure-species thermal conductivity of species k at T (K), W/(m K): Warnatz/Kee model,
 * except H2O (IAPWS 2011 dilute gas; see header).
 */
export function speciesThermalConductivity(k: number, T: number): number {
  return k === K_H2O ? waterConductivity(T) : speciesThermalConductivityWarnatz(k, T);
}

/**
 * Pure-species thermal conductivity by the modified Eucken correlation, W/(m K):
 * λ M/μ = 1.32 c_v + 1.77 R (Svehla 1962, NASA TR R-132; Poling, Prausnitz & O'Connell
 * 2001, ch. 10), with the Chapman–Enskog μ. Reduces to the exact 15/4 (R/M) μ for
 * monatomic gases. For comparison only.
 * UNVERIFIED: equation numbering in Poling et al. (formula itself is standard).
 */
export function speciesThermalConductivityModifiedEucken(k: number, T: number): number {
  const mu = speciesViscosityChapmanEnskog(k, T);
  const cv = speciesCpR(k, T) - 1; // c_v/R
  return ((mu / MOLAR_MASS[k]) * R_UNIVERSAL) * (1.32 * cv + 1.77);
}

// ---- mixtures --------------------------------------------------------------------------

const scrMu = new Float64Array(NS);
const scrSqrtMu = new Float64Array(NS);
const scrLam = new Float64Array(NS);
const scrIdx = new Int32Array(NS);

/** Fill scratch with pure μ_k (and optionally λ_k) for the species with X_k > 0; returns count. */
function fillPure(X: Float64Array, T: number, withConductivity: boolean): number {
  let n = 0;
  const sqrtT = Math.sqrt(T);
  for (let k = 0; k < NS; k++) {
    if (!(X[k] > 0)) continue;
    let mu: number;
    if (k === K_H2O) {
      mu = waterViscosity(T);
      if (withConductivity) scrLam[k] = waterConductivity(T);
    } else {
      evalCollision(k, T);
      mu = (VISC_PREF[k] * sqrtT) / COL[0];
      if (withConductivity) scrLam[k] = keeConductivity(k, T, mu);
    }
    scrMu[k] = mu;
    scrSqrtMu[k] = Math.sqrt(mu);
    scrIdx[n++] = k;
  }
  return n;
}

/**
 * Mixture dynamic viscosity at T (K), Pa s — Wilke (1950):
 * μ = Σ_k X_k μ_k / Σ_j X_j Φ_kj,  Φ_kj = [1 + (μ_k/μ_j)^{1/2}(M_j/M_k)^{1/4}]² / [8(1 + M_k/M_j)]^{1/2}.
 * Species with X_k ≤ 0 are ignored; X need not be normalised.
 */
export function mixViscosity(X: Float64Array, T: number): number {
  const n = fillPure(X, T, false);
  return wilke(X, n);
}

function wilke(X: Float64Array, n: number): number {
  let mu = 0;
  for (let a = 0; a < n; a++) {
    const k = scrIdx[a];
    let den = 0;
    const row = k * NS;
    for (let b = 0; b < n; b++) {
      const j = scrIdx[b];
      const t = 1 + (scrSqrtMu[k] / scrSqrtMu[j]) * WILKE_M4[row + j];
      den += X[j] * t * t * WILKE_DEN[row + j];
    }
    mu += (X[k] * scrMu[k]) / den;
  }
  return mu;
}

/** Mixture thermal conductivity at T (K), W/(m K) — Mathur–Saxena average of the pure λ_k. */
export function mixThermalConductivity(X: Float64Array, T: number): number {
  const n = fillPure(X, T, true);
  return mathurSaxena(X, n);
}

function mathurSaxena(X: Float64Array, n: number): number {
  let s1 = 0;
  let s2 = 0;
  let sx = 0;
  for (let a = 0; a < n; a++) {
    const k = scrIdx[a];
    s1 += X[k] * scrLam[k];
    s2 += X[k] / scrLam[k];
    sx += X[k];
  }
  // Normalised form (identical to ½(Σ Xλ + 1/Σ X/λ) when Σ X = 1).
  return 0.5 * (s1 / sx + sx / s2);
}

/** Mixture thermal diffusivity λ/(ρ c_p) at T (K), p (Pa), m²/s (X normalised). */
export function mixThermalDiffusivity(X: Float64Array, T: number, p: number): number {
  const rho = (p * mixMolarMass(X)) / (R_UNIVERSAL * T);
  return mixThermalConductivity(X, T) / (rho * mixCpMass(X, T));
}

/** Mixture kinematic viscosity μ/ρ at T (K), p (Pa), m²/s (X normalised). */
export function mixKinematicViscosity(X: Float64Array, T: number, p: number): number {
  const rho = (p * mixMolarMass(X)) / (R_UNIVERSAL * T);
  return mixViscosity(X, T) / rho;
}

/** Output of mixTransport. */
export interface TransportProps {
  /** Dynamic viscosity, Pa s. */
  viscosity: number;
  /** Thermal conductivity, W/(m K). */
  conductivity: number;
  /** Thermal diffusivity λ/(ρ c_p), m²/s. */
  thermalDiffusivity: number;
  /** Kinematic viscosity μ/ρ, m²/s. */
  kinematicViscosity: number;
  /** Prandtl number μ c_p/λ. */
  prandtl: number;
}

/**
 * All mixture transport properties at T (K), p (Pa) in one pass (pure-species collision
 * integrals evaluated once; X normalised). Writes into `out` (allocation-free when supplied).
 */
export function mixTransport(
  X: Float64Array,
  T: number,
  p: number,
  out: TransportProps = { viscosity: 0, conductivity: 0, thermalDiffusivity: 0, kinematicViscosity: 0, prandtl: 0 },
): TransportProps {
  const n = fillPure(X, T, true);
  const mu = wilke(X, n);
  const lam = mathurSaxena(X, n);
  const rho = (p * mixMolarMass(X)) / (R_UNIVERSAL * T);
  const cp = mixCpMass(X, T);
  out.viscosity = mu;
  out.conductivity = lam;
  out.thermalDiffusivity = lam / (rho * cp);
  out.kinematicViscosity = mu / rho;
  out.prandtl = (mu * cp) / lam;
  return out;
}
