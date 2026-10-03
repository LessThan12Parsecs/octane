/**
 * Physical constants and unit helpers. Everything in src/physics is SI:
 * m, kg, s, K, Pa, J, mol (NOT kmol), W, V, A. Angles are radians internally.
 *
 * Atomic weights match Cantera 3.x so our results can be compared to the
 * Cantera oracle to ~1e-9 relative.
 */

/** Universal gas constant, J/(mol K) (Cantera value, CODATA 2018 derived). */
export const R_UNIVERSAL = 8.31446261815324;
/** Standard atmosphere, Pa. */
export const P_ATM = 101325;
/** One bar, Pa. */
export const P_BAR = 1e5;
/** Thermochemical reference temperature, K. */
export const T_REF = 298.15;
/** Boltzmann constant, J/K. */
export const K_BOLTZMANN = 1.380649e-23;
/** Avogadro constant, 1/mol. */
export const N_AVOGADRO = 6.02214076e23;
/** Elementary charge, C. */
export const Q_ELECTRON = 1.602176634e-19;
/** Stefan-Boltzmann constant, W/(m^2 K^4). */
export const SIGMA_SB = 5.670374419e-8;

export const DEG = Math.PI / 180;
export const RAD2DEG = 180 / Math.PI;

export const rpmToOmega = (rpm: number): number => (rpm * 2 * Math.PI) / 60;
export const omegaToRpm = (omega: number): number => (omega * 60) / (2 * Math.PI);

/** Atomic weights, kg/mol (Cantera 3.x values). */
export const ATOMIC_WEIGHT = {
  C: 12.011e-3,
  H: 1.008e-3,
  O: 15.999e-3,
  N: 14.007e-3,
  AR: 39.95e-3,
} as const;
