/**
 * Fuels, air and fresh-charge / product compositions. SI units (mol, kg, J, K, Pa).
 *
 * All composition vectors are mole fractions (Float64Array(NS), SPECIES order) unless the
 * name says N (moles). Fuels are handled as fully vaporised (gaseous) blends of the fuel
 * species in core/species.ts.
 */
import { R_UNIVERSAL, T_REF } from '../core/constants';
import type { FuelSelection } from '../core/operating-point';
import { EL, ELEMENT_COUNTS, FUEL_SPECIES, NE, NS, SP, type FuelSpeciesName } from '../core/species';
import { elementMoles, mixMolarMass } from './mixture';
import { MOLAR_MASS, speciesG0, speciesH } from './thermo';

/** A fully vaporised fuel blend. */
export interface FuelBlend {
  /** Display label, e.g. 'PRF 90' or 'iso-octane'. */
  label: string;
  /** Fuel-vapour mole fractions (Float64Array(NS), non-zero only at fuel species, sums to 1). */
  X: Float64Array;
  /** Octane number on the PRF scale where defined (PRF blends; iso-octane 100, n-heptane 0). */
  octaneNumber?: number;
}

// ---------------------------------------------------------------------------------------
// Liquid-fuel data
// ---------------------------------------------------------------------------------------

/**
 * Temperature at which PRF blends are proportioned by volume, K: 60 °F = 15.556 °C, the
 * petroleum-measurement reference temperature used by the ASTM D2699/D2700 rating methods.
 * The volume-% definition of the PRF scale is confirmed (ASTM D2699 scope/terminology via
 * secondary sources); UNVERIFIED: that the standard fixes the blending temperature at
 * 60 °F (the standard text was not accessible). The iso-octane/n-heptane density ratio
 * changes by ≲ 4e-4 between 15 °C and 25 °C (x_iso by ≲ 1e-4), so the choice is immaterial.
 */
export const PRF_BLEND_TEMPERATURE = 273.15 + (60 - 32) / 1.8;

/**
 * Liquid iso-octane density at t (°C), kg/m³ — NIST SRM 2214 certificate (issued
 * 20 Feb 2001), eq. (1):
 *   ρ(t) = 691.872 [1 − 1.18752e-3 (t − 20) − 6.2516e-7 (t − 20)²],  valid 15–25 °C
 * (certified 695.969 / 691.872 / 687.753 kg/m³ at 15 / 20 / 25 °C, ±0.035 kg/m³; certificate
 * re-fetched 2026-09-29 and the table values reproduced by the equation to < 1e-3 kg/m³).
 */
export function isooctaneLiquidDensity(tCelsius: number): number {
  const d = tCelsius - 20;
  return 691.872 * (1 - 1.18752e-3 * d - 6.2516e-7 * d * d);
}

/**
 * Liquid n-heptane density at 60 °F (288.706 K) and 0.101325 MPa, kg/m³ — NIST Chemistry
 * WebBook fluid properties (EOS of Tenji, Thol, Lemmon & Span 2018), fetched 2026-09-29:
 * 688.065 kg/m³ at 288.150 K, 687.597 kg/m³ at 288.706 K.
 */
export const NHEPTANE_DENSITY_60F = 687.597;

/** Liquid iso-octane density at 60 °F, kg/m³ (NIST SRM 2214 eq. (1); = 695.515). */
export const ISOOCTANE_DENSITY_60F = isooctaneLiquidDensity(PRF_BLEND_TEMPERATURE - 273.15);

/**
 * Enthalpy of vaporisation at 298.15 K, J/mol (liquid → ideal gas; used only for the
 * liquid-basis heating value). NaN where there is no liquid at 298.15 K.
 *  - IC8H18: 35.14 kJ/mol — Majer & Svoboda (1985) correlation A·exp(−βTr)(1−Tr)^β with
 *    A = 50.28 kJ/mol, β = 0.2668, Tc = 543.9 K (NIST WebBook); NIST average ΔvapH° = 35.1 ± 0.2.
 *  - NC7H16: 36.57 kJ/mol — Majer & Svoboda (1985), A = 53.66 kJ/mol, β = 0.2831, Tc = 540.2 K
 *    (NIST WebBook).
 *  - C2H5OH: 42.3 kJ/mol — NIST WebBook ΔvapH°, average of 12 values, ±0.4 kJ/mol.
 *  - C3H8: 16.25 kJ/mol — NIST WebBook ΔvapH° (Majer & Svoboda 1985): saturated liquid
 *    (p_sat = 0.952 MPa) → IDEAL gas at 298.15 K, which is the quantity needed here. Checked
 *    against the NIST WebBook EOS: h_ig(298.15 K) − h_liq,sat = 27.9216 − 11.6900 = 16.23 kJ/mol
 *    (the sat. liquid → sat. VAPOUR value, 14.80 kJ/mol, omits the vapour non-ideality).
 *    For the low-volatility liquids above, the difference is negligible (p_sat ≤ 8 kPa).
 *  - CH4: supercritical at 298.15 K (Tc = 190.6 K) → NaN.
 */
export const HEAT_OF_VAPORIZATION_298: Readonly<Record<FuelSpeciesName, number>> = {
  IC8H18: majerSvoboda(50.28e3, 0.2668, 543.9, T_REF),
  NC7H16: majerSvoboda(53.66e3, 0.2831, 540.2, T_REF),
  C2H5OH: 42.3e3,
  C3H8: 16.25e3,
  CH4: NaN,
};

/** Majer & Svoboda (1985) ΔvapH(T) = A exp(−β Tr) (1 − Tr)^β, J/mol (A in J/mol). */
function majerSvoboda(A: number, beta: number, Tc: number, T: number): number {
  const Tr = T / Tc;
  return A * Math.exp(-beta * Tr) * Math.pow(1 - Tr, beta);
}

const FUEL_LABEL: Readonly<Record<FuelSpeciesName, string>> = {
  IC8H18: 'iso-octane',
  NC7H16: 'n-heptane',
  CH4: 'methane',
  C3H8: 'propane',
  C2H5OH: 'ethanol',
};

// ---------------------------------------------------------------------------------------
// Fuel selection
// ---------------------------------------------------------------------------------------

/**
 * Mole fraction of iso-octane in the vapour of a PRF blend of octane number ON (liquid-volume
 * % iso-octane in n-heptane, 0..100). Complete vaporisation → vapour composition = liquid
 * mole fractions: n_i ∝ V_i ρ_i / M_i with the 60 °F liquid densities above (ideal
 * (additive) liquid volumes; the iso-octane/n-heptane excess volume is neglected).
 */
export function prfIsooctaneMoleFraction(octaneNumber: number): number {
  const v = octaneNumber / 100;
  const nIso = (v * ISOOCTANE_DENSITY_60F) / MOLAR_MASS[SP.IC8H18];
  const nHep = ((1 - v) * NHEPTANE_DENSITY_60F) / MOLAR_MASS[SP.NC7H16];
  return nIso / (nIso + nHep);
}

/** Build the vapour-phase FuelBlend for an operating-point fuel selection. */
export function fuelFromSelection(sel: FuelSelection): FuelBlend {
  const X = new Float64Array(NS);
  if (sel.kind === 'PRF') {
    const on = sel.octaneNumber;
    if (!(on >= 0 && on <= 100)) throw new RangeError(`PRF octane number ${on} outside 0..100`);
    const x = prfIsooctaneMoleFraction(on);
    X[SP.IC8H18] = x;
    X[SP.NC7H16] = 1 - x;
    return { label: `PRF ${on}`, X, octaneNumber: on };
  }
  X[SP[sel.species]] = 1;
  const blend: FuelBlend = { label: FUEL_LABEL[sel.species], X };
  if (sel.species === 'IC8H18') blend.octaneNumber = 100;
  else if (sel.species === 'NC7H16') blend.octaneNumber = 0;
  return blend;
}

/** Mean molar mass of the fuel vapour, kg/mol. */
export const fuelMolarMass = (fuel: FuelBlend): number => mixMolarMass(fuel.X);

// ---------------------------------------------------------------------------------------
// Air
// ---------------------------------------------------------------------------------------

/**
 * Dry-air reference composition, mole fractions — Picard, Davis, Gläser & Fujii (2008),
 * "Revised formula for the density of moist air (CIPM-2007)", Metrologia 45:149, Table 1:
 * N2 0.780848, O2 0.209390, Ar 0.009332, CO2 0.00040 (the remaining 30 µmol/mol — Ne, He,
 * CH4, Kr, H2, N2O, CO, Xe — is omitted and the four are renormalised; M = 28.9660 g/mol
 * with the Cantera atomic weights vs CIPM 28.96546 g/mol, 1.9e-5 relative).
 */
export const DRY_AIR_REFERENCE = { N2: 0.780848, O2: 0.20939, AR: 0.009332, CO2: 0.0004 } as const;

/** Dry air, mole fractions (new array; see DRY_AIR_REFERENCE). */
export function dryAir(out: Float64Array = new Float64Array(NS)): Float64Array {
  const a = DRY_AIR_REFERENCE;
  const s = a.N2 + a.O2 + a.AR + a.CO2;
  out.fill(0);
  out[SP.N2] = a.N2 / s;
  out[SP.O2] = a.O2 / s;
  out[SP.AR] = a.AR / s;
  out[SP.CO2] = a.CO2 / s;
  return out;
}

/** Triple point of water, K and Pa (Wagner & Pruss 1993, Table 1 check value). */
const T_TRIPLE = 273.16;
/** Critical point of water, K and Pa (Wagner & Pruss 1993, §2). */
const T_CRIT_W = 647.096;
const P_CRIT_W = 22.064e6;

/**
 * Saturation vapour pressure of water over LIQUID water at T (K), Pa.
 *  - 273.16 K ≤ T ≤ 647.096 K: Wagner & Pruss (1993), J. Phys. Chem. Ref. Data 22:783,
 *    eq. (1) (IAPWS saturation equation; check values 611.657 Pa at 273.16 K,
 *    101325 Pa at 373.1243 K reproduced in fuels.test.ts);
 *  - T < 273.16 K (supercooled liquid — WMO convention: relative humidity is always
 *    referred to liquid water): Murphy & Koop (2005), Q. J. R. Meteorol. Soc. 131:1539,
 *    eq. (10), valid 123–332 K (joins Wagner–Pruss to 7e-8 at the triple point;
 *    coefficients checked against the paper's eq. (10), PDF fetched 2026-09-29);
 *  - T > Tc: returns p_c.
 */
export function waterSaturationPressure(T: number): number {
  if (T >= T_CRIT_W) return P_CRIT_W;
  if (T >= T_TRIPLE) {
    const tau = 1 - T / T_CRIT_W;
    const t15 = tau * Math.sqrt(tau);
    const t3 = tau * tau * tau;
    const s =
      -7.85951783 * tau +
      1.84408259 * t15 +
      -11.7866497 * t3 +
      22.6807411 * t3 * Math.sqrt(tau) +
      -15.9618719 * t3 * tau +
      1.80122502 * t3 * t3 * t15;
    return P_CRIT_W * Math.exp((T_CRIT_W / T) * s);
  }
  const lnT = Math.log(T);
  return Math.exp(
    54.842763 - 6763.22 / T - 4.21 * lnT + 0.000367 * T +
      Math.tanh(0.0415 * (T - 218.8)) * (53.878 - 1331.22 / T - 9.44523 * lnT + 0.014025 * T),
  );
}

/**
 * Water-vapour enhancement factor of moist air f(p, t) = α + β p + γ t² — CIPM-2007
 * (Picard et al. 2008, eq. A1.2): α = 1.00062, β = 3.14e-8 Pa⁻¹, γ = 5.6e-7 K⁻², t in °C.
 * Stated validity 60–110 kPa, 15–27 °C; used as an extrapolation outside that range.
 */
export function moistAirEnhancementFactor(T: number, p: number): number {
  const t = T - 273.15;
  return 1.00062 + 3.14e-8 * p + 5.6e-7 * t * t;
}

/**
 * Humid air at temperature T (K), pressure p (Pa) and relative humidity RH (0..1, referred
 * to liquid water), mole fractions: x_H2O = RH f(p,T) p_sat(T) / p (CIPM-2007 eq. A1.3),
 * remainder dry air (DRY_AIR_REFERENCE). x_H2O is clamped to [0, 1].
 */
export function humidAir(T: number, p: number, RH: number, out: Float64Array = new Float64Array(NS)): Float64Array {
  let xv = (RH * moistAirEnhancementFactor(T, p) * waterSaturationPressure(T)) / p;
  xv = Math.min(1, Math.max(0, xv));
  dryAir(out);
  for (let k = 0; k < NS; k++) out[k] *= 1 - xv;
  out[SP.H2O] += xv;
  return out;
}

// ---------------------------------------------------------------------------------------
// Stoichiometry
// ---------------------------------------------------------------------------------------

const scratchB = new Float64Array(NE);

/**
 * O2 required for complete combustion (to CO2 + H2O) per mole of fuel vapour:
 * ν = C + H/4 − O/2 (atoms per mean fuel molecule). mol O2 / mol fuel.
 */
export function stoichO2PerMolFuel(fuel: FuelBlend): number {
  const b = elementMoles(fuel.X, scratchB);
  let n = 0;
  for (let k = 0; k < NS; k++) n += fuel.X[k];
  return (b[EL.C] + b[EL.H] / 4 - b[EL.O] / 2) / n;
}

/**
 * Stoichiometric air/fuel MASS ratio for the given air composition (mole fractions; include
 * H2O for humid air — the ratio is then per kg of humid air). kg air / kg fuel.
 */
export function stoichAirFuelRatio(fuel: FuelBlend, airX: Float64Array): number {
  const molAirPerMolFuel = stoichO2PerMolFuel(fuel) / airX[SP.O2];
  return (molAirPerMolFuel * mixMolarMass(airX)) / fuelMolarMass(fuel);
}

/** Options for freshCharge. */
export interface FreshChargeSpec {
  fuel: FuelBlend;
  /** Fuel/air equivalence ratio, defined against the O2 of the (fresh) air only. */
  phi: number;
  /** Air composition, mole fractions (dryAir() / humidAir()). */
  airX: Float64Array;
  /** External EGR as a MASS fraction of the total fresh charge (fuel + air + EGR), 0..1. */
  egrFraction?: number;
  /** EGR composition (mole fractions). Default: completeCombustionProducts of the fuel–air mixture. */
  egrX?: Float64Array;
}

/**
 * Fresh-charge mole fractions: air + fully vaporised fuel at φ (n_fuel/n_air =
 * φ x_O2,air / ν_O2), plus optional EGR given as a mass fraction of the total charge.
 */
export function freshCharge(spec: FreshChargeSpec, out: Float64Array = new Float64Array(NS)): Float64Array {
  const { fuel, phi, airX } = spec;
  const egr = spec.egrFraction ?? 0;
  if (!(phi >= 0)) throw new RangeError(`equivalence ratio ${phi} must be ≥ 0`);
  if (!(egr >= 0 && egr < 1)) throw new RangeError(`EGR mass fraction ${egr} outside [0, 1)`);
  const nFuel = (phi * airX[SP.O2]) / stoichO2PerMolFuel(fuel); // per mol air
  const fa = new Float64Array(NS);
  for (let k = 0; k < NS; k++) fa[k] = airX[k] + nFuel * fuel.X[k];
  // fa holds moles per mole of air; mass of that fuel–air mixture:
  let mFA = 0;
  let nFA = 0;
  for (let k = 0; k < NS; k++) {
    mFA += fa[k] * MOLAR_MASS[k];
    nFA += fa[k];
  }
  out.fill(0);
  for (let k = 0; k < NS; k++) out[k] = fa[k];
  let nTot = nFA;
  if (egr > 0) {
    let egrX = spec.egrX;
    if (!egrX) {
      const faX = new Float64Array(NS);
      for (let k = 0; k < NS; k++) faX[k] = fa[k] / nFA;
      egrX = completeCombustionProducts(faX);
    }
    let sx = 0;
    for (let k = 0; k < NS; k++) sx += egrX[k];
    const mEgr = (egr / (1 - egr)) * mFA;
    const nEgr = mEgr / (mixMolarMass(egrX) / sx); // egrX need not be normalised
    for (let k = 0; k < NS; k++) out[k] += (nEgr * egrX[k]) / sx;
    nTot += nEgr;
  }
  for (let k = 0; k < NS; k++) out[k] /= nTot;
  return out;
}

// ---------------------------------------------------------------------------------------
// Complete-combustion products
// ---------------------------------------------------------------------------------------

/**
 * Temperature at which the water-gas-shift equilibrium is frozen for rich products, K.
 * Heywood (1988), Internal Combustion Engine Fundamentals, ch. 4 (rich combustion-product
 * composition) uses K = 3.5, typical of ≈1740 K; our NASA data give K(1740 K) = 3.58.
 * UNVERIFIED: exact section/equation of the Heywood reference (the value 3.5 at 1740 K was
 * confirmed only via a secondary web source).
 */
export const WGS_FREEZE_TEMPERATURE = 1740;

/**
 * Equilibrium constant of CO2 + H2 ⇌ CO + H2O, K = x_CO x_H2O / (x_CO2 x_H2) = exp(−ΔG°/RT)
 * from the species data (Δn = 0, so pressure-independent). Dimensionless.
 */
export function waterGasShiftK(T: number): number {
  const dG = speciesG0(SP.CO, T) + speciesG0(SP.H2O, T) - speciesG0(SP.CO2, T) - speciesG0(SP.H2, T);
  return Math.exp(-dG / (R_UNIVERSAL * T));
}

const FUEL_IDX: readonly number[] = FUEL_SPECIES.map((s) => SP[s]);
const scratchBC = new Float64Array(NE);

/**
 * Complete-combustion product MOLES from reactant moles (element-conserving; extensive).
 * Element pools are taken from ALL reactant species (so residual/EGR CO, H2, OH, NO … are
 * re-allocated too). N → N2, Ar → Ar.
 *  - Lean/stoichiometric (O ≥ 2C + H/2): CO2 = C, H2O = H/2, O2 = (O − 2C − H/2)/2.
 *  - Rich: CO2, CO, H2O, H2 from the C, H, O balances plus the water-gas-shift equilibrium
 *    K(T_wgs) = n_CO n_H2O / (n_CO2 n_H2) (default T_wgs = WGS_FREEZE_TEMPERATURE), solved
 *    in closed form (quadratic).
 *  - Beyond the CO limit (O < C, i.e. φ ≳ 3 for alkanes): only the fraction of the fuel
 *    species that the oxygen can take to CO + H2 is burned; the rest stays as fuel.
 * Writes into outN if given (must not alias N). Returns outN.
 */
export function completeCombustionMoles(
  N: Float64Array,
  outN: Float64Array = new Float64Array(NS),
  wgsTemperature: number = WGS_FREEZE_TEMPERATURE,
): Float64Array {
  const b = elementMoles(N, scratchBC);
  outN.fill(0);
  let nC = b[EL.C];
  let nH = b[EL.H];
  let nO = b[EL.O];
  outN[SP.N2] = b[EL.N] / 2;
  outN[SP.AR] = b[EL.AR];

  if (nO < nC) {
    // Not even enough O to take all carbon to CO: burn only a fraction ζ of the fuel
    // species (all burned C → CO, burned H → H2); the rest is carried through as fuel.
    let cFuel = 0;
    let hFuel = 0;
    let oFuel = 0;
    for (const k of FUEL_IDX) {
      const n = N[k];
      if (n === 0) continue;
      const a = ELEMENT_COUNTS[k];
      cFuel += n * a[EL.C];
      hFuel += n * a[EL.H];
      oFuel += n * a[EL.O];
    }
    const cOther = nC - cFuel;
    const oOther = nO - oFuel;
    // O_burned = C_burned:  oOther + ζ oFuel = cOther + ζ cFuel
    const zeta = cFuel > oFuel ? Math.max(0, Math.min(1, (oOther - cOther) / (cFuel - oFuel))) : 0;
    for (const k of FUEL_IDX) outN[k] = N[k] * (1 - zeta);
    nC -= (1 - zeta) * cFuel;
    nH -= (1 - zeta) * hFuel;
    nO -= (1 - zeta) * oFuel;
  }

  const lean = nO >= 2 * nC + nH / 2;
  if (lean) {
    outN[SP.CO2] = nC;
    outN[SP.H2O] = nH / 2;
    outN[SP.O2] = (nO - 2 * nC - nH / 2) / 2;
  } else {
    // a = CO2, b = CO = C − a, c = H2O = O − C − a, d = H2 = H/2 − O + C + a
    // K a d = b c  →  (K − 1) a² + (K A + O) a − C D = 0,  A = H/2 − O + C, D = O − C
    const K = waterGasShiftK(wgsTemperature);
    const A = nH / 2 - nO + nC;
    const D = Math.max(0, nO - nC);
    const B = K * A + nC + D;
    const CD = nC * D;
    const disc = B * B + 4 * (K - 1) * CD;
    // Stable root of the quadratic (positive root, no cancellation):
    let a = CD > 0 ? (2 * CD) / (B + Math.sqrt(Math.max(0, disc))) : 0;
    a = Math.min(Math.max(a, 0), Math.min(nC, D));
    outN[SP.CO2] = a;
    outN[SP.CO] = nC - a;
    outN[SP.H2O] = Math.max(0, D - a);
    outN[SP.H2] = Math.max(0, A + a);
  }
  return outN;
}

/**
 * Complete-combustion products of a reactant mixture (mole fractions in → mole fractions
 * out). See completeCombustionMoles for the model. Intended for initial residual/EGR
 * guesses and fuel–air-cycle estimates, not for burned-gas equilibrium.
 */
export function completeCombustionProducts(
  reactantsX: Float64Array,
  out: Float64Array = new Float64Array(NS),
  wgsTemperature: number = WGS_FREEZE_TEMPERATURE,
): Float64Array {
  completeCombustionMoles(reactantsX, out, wgsTemperature);
  let s = 0;
  for (let k = 0; k < NS; k++) s += out[k];
  for (let k = 0; k < NS; k++) out[k] /= s;
  return out;
}

// ---------------------------------------------------------------------------------------
// Heating value
// ---------------------------------------------------------------------------------------

/**
 * Lower heating value of the GASEOUS fuel at 298.15 K per mole of fuel, J/mol:
 * −ΔH of  fuel(g) + ν O2 → C CO2 + H/2 H2O(g)  (+ ½N2 per N atom), from the species data.
 */
export function lowerHeatingValueMolar(fuel: FuelBlend): number {
  const T = T_REF;
  const b = elementMoles(fuel.X, scratchB);
  let n = 0;
  let hFuel = 0;
  for (let k = 0; k < NS; k++) {
    const x = fuel.X[k];
    if (x === 0) continue;
    n += x;
    hFuel += x * speciesH(k, T);
  }
  const nu = b[EL.C] + b[EL.H] / 4 - b[EL.O] / 2;
  const hReact = hFuel + nu * speciesH(SP.O2, T) + (b[EL.N] / 2) * speciesH(SP.N2, T);
  const hProd = b[EL.C] * speciesH(SP.CO2, T) + (b[EL.H] / 2) * speciesH(SP.H2O, T) +
    (b[EL.N] / 2) * speciesH(SP.N2, T);
  return (hReact - hProd) / n;
}

/** Lower heating value of the GASEOUS (vaporised) fuel at 298.15 K, J/kg fuel. */
export const lowerHeatingValue = (fuel: FuelBlend): number =>
  lowerHeatingValueMolar(fuel) / fuelMolarMass(fuel);

/**
 * Lower heating value of the LIQUID fuel at 298.15 K, J/kg fuel = gaseous LHV − Σ x_k
 * ΔvapH_k / M (ideal liquid solution). NaN if a component has no liquid state at 298 K.
 * For comparison with handbook values (e.g. iso-octane 44.3 MJ/kg, Heywood 1988 App. D).
 */
export function lowerHeatingValueLiquid(fuel: FuelBlend): number {
  let hv = 0;
  let n = 0;
  for (const name of FUEL_SPECIES) {
    const x = fuel.X[SP[name]];
    if (x === 0) continue;
    hv += x * HEAT_OF_VAPORIZATION_298[name];
    n += x;
  }
  return (lowerHeatingValueMolar(fuel) - hv / n) / fuelMolarMass(fuel);
}
