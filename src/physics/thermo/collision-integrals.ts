/**
 * Reduced collision integrals for Chapman–Enskog transport.
 *
 * Non-polar (Lennard-Jones 12-6) molecules: the empirical fits of
 *   Neufeld, Janzen & Aziz (1972), J. Chem. Phys. 57:1100 —
 *   Ω(2,2)* = 1.16145 T*^−0.14874 + 0.52487 e^(−0.77320 T*) + 2.16178 e^(−2.43787 T*)
 *             − 6.435e−4 T*^0.14874 sin(18.0323 T*^−0.76830 − 7.27371)
 *   Ω(1,1)* = 1.06036 T*^−0.15610 + 0.19300 e^(−0.47635 T*) + 1.03587 e^(−1.52996 T*)
 *             + 1.76474 e^(−3.89411 T*)
 * (Ω(2,2)* coefficients checked against CoolProp's implementation of NJA 1972, Ω(1,1)*
 * against a published reproduction of the NJA diffusion integral; stated validity
 * 0.3 ≤ T* ≤ 100, used slightly beyond — the fits are smooth power laws there; agreement
 * with the Monchick–Mason δ* = 0 column is ≤ 0.5 %, see transport.test.ts).
 *
 * Polar molecules (Stockmayer 12-6-3 potential, reduced dipole δ* = μ²/(2 ε σ³) in Gaussian
 * units = ½ μ²/(4π ε0 ε σ³) in SI): the Neufeld values are multiplied by the ratio
 * Ω(T*, δ*)/Ω(T*, 0) taken from the Monchick & Mason (1961) tables,
 *   L. Monchick & E.A. Mason, "Transport properties of polar gases",
 *   J. Chem. Phys. 35:1676 (1961),
 * transcribed programmatically from Cantera 3.2.0 src/transport/MMCollisionInt.cpp (the same
 * tables CHEMKIN's TRANFIT uses). Interpolation: cubic Lagrange in δ* over the four nearest
 * table columns (done once per species at load), quadratic in ln T* over three rows (runtime),
 * as in Cantera. The ratio form keeps δ* = 0 identical to Neufeld and makes the correction
 * vanish smoothly (→1) at high T*. (Brokaw's 1969 closed form Ω_D = Ω_D,LJ + 0.19 δ²/T* was
 * tried first; it is 3–4 % off the tables for H2O near 1000 K, so the tables are used.)
 */

/** δ* values of the table columns. */
const MM_DELTA: readonly number[] = [0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5];
const ND = 8;
const NT = 37;

/** Reduced temperatures T* of the table rows (37 values, 0.1 … 100). */
const MM_TSTAR: readonly number[] = [
  0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1,
  1.2, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5, 4, 5,
  6, 7, 8, 9, 10, 12, 14, 16, 18, 20,
  25, 30, 35, 40, 50, 75, 100,
];

/** Ω(2,2)*(T*, δ*) — rows: MM_TSTAR, columns: MM_DELTA. */
const MM_OMEGA22: readonly number[] = [
  4.1005, 4.266, 4.833, 5.742, 6.729, 8.624, 10.34, 11.89,
  3.2626, 3.305, 3.516, 3.914, 4.433, 5.57, 6.637, 7.618,
  2.8399, 2.836, 2.936, 3.168, 3.511, 4.329, 5.126, 5.874,
  2.531, 2.522, 2.586, 2.749, 3.004, 3.64, 4.282, 4.895,
  2.2837, 2.277, 2.329, 2.46, 2.665, 3.187, 3.727, 4.249,
  2.0838, 2.081, 2.13, 2.243, 2.417, 2.862, 3.329, 3.786,
  1.922, 1.924, 1.97, 2.072, 2.225, 2.614, 3.028, 3.435,
  1.7902, 1.795, 1.84, 1.934, 2.07, 2.417, 2.788, 3.156,
  1.6823, 1.689, 1.733, 1.82, 1.944, 2.258, 2.596, 2.933,
  1.5929, 1.601, 1.644, 1.725, 1.838, 2.124, 2.435, 2.746,
  1.4551, 1.465, 1.504, 1.574, 1.67, 1.913, 2.181, 2.451,
  1.3551, 1.365, 1.4, 1.461, 1.544, 1.754, 1.989, 2.228,
  1.28, 1.289, 1.321, 1.374, 1.447, 1.63, 1.838, 2.053,
  1.2219, 1.231, 1.259, 1.306, 1.37, 1.532, 1.718, 1.912,
  1.1757, 1.184, 1.209, 1.251, 1.307, 1.451, 1.618, 1.795,
  1.0933, 1.1, 1.119, 1.15, 1.193, 1.304, 1.435, 1.578,
  1.0388, 1.044, 1.059, 1.083, 1.117, 1.204, 1.31, 1.428,
  0.99963, 1.004, 1.016, 1.035, 1.062, 1.133, 1.22, 1.319,
  0.96988, 0.9732, 0.983, 0.9991, 1.021, 1.079, 1.153, 1.236,
  0.92676, 0.9291, 0.936, 0.9473, 0.9628, 1.005, 1.058, 1.121,
  0.89616, 0.8979, 0.903, 0.9114, 0.923, 0.9545, 0.9955, 1.044,
  0.87272, 0.8741, 0.878, 0.8845, 0.8935, 0.9181, 0.9505, 0.9893,
  0.85379, 0.8549, 0.858, 0.8632, 0.8703, 0.8901, 0.9164, 0.9482,
  0.83795, 0.8388, 0.8414, 0.8456, 0.8515, 0.8678, 0.8895, 0.916,
  0.82435, 0.8251, 0.8273, 0.8308, 0.8356, 0.8493, 0.8676, 0.8901,
  0.80184, 0.8024, 0.8039, 0.8065, 0.8101, 0.8201, 0.8337, 0.8504,
  0.78363, 0.784, 0.7852, 0.7872, 0.7899, 0.7976, 0.8081, 0.8212,
  0.76834, 0.7687, 0.7696, 0.7712, 0.7733, 0.7794, 0.7878, 0.7983,
  0.75518, 0.7554, 0.7562, 0.7575, 0.7592, 0.7642, 0.7711, 0.7797,
  0.74364, 0.7438, 0.7445, 0.7455, 0.747, 0.7512, 0.7569, 0.7642,
  0.71982, 0.72, 0.7204, 0.7211, 0.7221, 0.725, 0.7289, 0.7339,
  0.70097, 0.7011, 0.7014, 0.7019, 0.7026, 0.7047, 0.7076, 0.7112,
  0.68545, 0.6855, 0.6858, 0.6861, 0.6867, 0.6883, 0.6905, 0.6932,
  0.67232, 0.6724, 0.6726, 0.6728, 0.6733, 0.6743, 0.6762, 0.6784,
  0.65099, 0.651, 0.6512, 0.6513, 0.6516, 0.6524, 0.6534, 0.6546,
  0.61397, 0.6141, 0.6143, 0.6145, 0.6147, 0.6148, 0.6148, 0.6147,
  0.5887, 0.5889, 0.5894, 0.59, 0.5903, 0.5901, 0.5895, 0.5885,
];

/** A*(T*, δ*) = Ω(2,2)* / Ω(1,1)* — rows: MM_TSTAR, columns: MM_DELTA. */
const MM_ASTAR: readonly number[] = [
  1.0231, 1.066, 1.038, 1.04, 1.043, 1.05, 1.052, 1.051,
  1.0424, 1.045, 1.048, 1.052, 1.056, 1.065, 1.066, 1.064,
  1.0719, 1.067, 1.06, 1.055, 1.058, 1.068, 1.071, 1.071,
  1.0936, 1.087, 1.077, 1.069, 1.068, 1.075, 1.078, 1.078,
  1.1053, 1.098, 1.088, 1.08, 1.078, 1.082, 1.084, 1.084,
  1.1104, 1.104, 1.096, 1.089, 1.086, 1.089, 1.09, 1.09,
  1.1114, 1.107, 1.1, 1.095, 1.093, 1.095, 1.096, 1.095,
  1.1104, 1.107, 1.102, 1.099, 1.098, 1.1, 1.1, 1.099,
  1.1086, 1.106, 1.102, 1.101, 1.101, 1.105, 1.105, 1.104,
  1.1063, 1.104, 1.103, 1.103, 1.104, 1.108, 1.109, 1.108,
  1.102, 1.102, 1.103, 1.105, 1.107, 1.112, 1.115, 1.115,
  1.0985, 1.099, 1.101, 1.104, 1.108, 1.115, 1.119, 1.12,
  1.096, 1.096, 1.099, 1.103, 1.108, 1.116, 1.121, 1.124,
  1.0943, 1.095, 1.099, 1.102, 1.108, 1.117, 1.123, 1.126,
  1.0934, 1.094, 1.097, 1.102, 1.107, 1.116, 1.123, 1.128,
  1.0926, 1.094, 1.097, 1.099, 1.105, 1.115, 1.123, 1.13,
  1.0934, 1.095, 1.097, 1.099, 1.104, 1.113, 1.122, 1.129,
  1.0948, 1.096, 1.098, 1.1, 1.103, 1.112, 1.119, 1.127,
  1.0965, 1.097, 1.099, 1.101, 1.104, 1.11, 1.118, 1.126,
  1.0997, 1.1, 1.101, 1.102, 1.105, 1.11, 1.116, 1.123,
  1.1025, 1.103, 1.104, 1.105, 1.106, 1.11, 1.115, 1.121,
  1.105, 1.105, 1.106, 1.107, 1.108, 1.111, 1.115, 1.12,
  1.1072, 1.107, 1.108, 1.108, 1.109, 1.112, 1.115, 1.119,
  1.1091, 1.109, 1.109, 1.11, 1.111, 1.113, 1.115, 1.119,
  1.1107, 1.111, 1.111, 1.111, 1.112, 1.114, 1.116, 1.119,
  1.1133, 1.114, 1.113, 1.114, 1.114, 1.115, 1.117, 1.119,
  1.1154, 1.115, 1.116, 1.116, 1.116, 1.117, 1.118, 1.12,
  1.1172, 1.117, 1.117, 1.118, 1.118, 1.118, 1.119, 1.12,
  1.1186, 1.119, 1.119, 1.119, 1.119, 1.119, 1.12, 1.121,
  1.1199, 1.12, 1.12, 1.12, 1.12, 1.121, 1.121, 1.122,
  1.1223, 1.122, 1.122, 1.122, 1.122, 1.123, 1.123, 1.124,
  1.1243, 1.124, 1.124, 1.124, 1.124, 1.124, 1.125, 1.125,
  1.1259, 1.126, 1.126, 1.126, 1.126, 1.126, 1.126, 1.126,
  1.1273, 1.127, 1.127, 1.127, 1.127, 1.127, 1.127, 1.128,
  1.1297, 1.13, 1.13, 1.13, 1.13, 1.13, 1.13, 1.129,
  1.1339, 1.134, 1.134, 1.135, 1.135, 1.134, 1.134, 1.132,
  1.1364, 1.137, 1.137, 1.138, 1.139, 1.138, 1.137, 1.135,
];

const LN_TSTAR = new Float64Array(NT);
for (let i = 0; i < NT; i++) LN_TSTAR[i] = Math.log(MM_TSTAR[i]);

/** Neufeld–Janzen–Aziz (1972) Ω(2,2)* for the LJ 12-6 potential (viscosity integral). */
export function omega22LJ(Ts: number): number {
  const p = Math.pow(Ts, 0.14874);
  return (
    1.16145 / p +
    0.52487 * Math.exp(-0.7732 * Ts) +
    2.16178 * Math.exp(-2.43787 * Ts) -
    6.435e-4 * p * Math.sin(18.0323 * Math.pow(Ts, -0.7683) - 7.27371)
  );
}

/** Neufeld–Janzen–Aziz (1972) Ω(1,1)* for the LJ 12-6 potential (diffusion integral). */
export function omega11LJ(Ts: number): number {
  return (
    1.06036 / Math.pow(Ts, 0.1561) +
    0.193 * Math.exp(-0.47635 * Ts) +
    1.03587 * Math.exp(-1.52996 * Ts) +
    1.76474 * Math.exp(-3.89411 * Ts)
  );
}

/** Cubic Lagrange interpolation of table row `row` (columns = MM_DELTA) at δ*. */
function interpDelta(table: readonly number[], row: number, d: number): number {
  const dd = Math.min(Math.max(d, 0), MM_DELTA[ND - 1]);
  let j = 1;
  while (j < ND - 2 && MM_DELTA[j] < dd) j++;
  const j0 = Math.max(0, Math.min(ND - 4, j - 2));
  let s = 0;
  for (let a = j0; a < j0 + 4; a++) {
    let w = 1;
    for (let c = j0; c < j0 + 4; c++) if (c !== a) w *= (dd - MM_DELTA[c]) / (MM_DELTA[a] - MM_DELTA[c]);
    s += w * table[row * ND + a];
  }
  return s;
}

/**
 * Polar correction factors for a species of reduced dipole δ*, tabulated on the MM_TSTAR
 * grid: r22 = Ω22*(T*,δ*)/Ω22*(T*,0), rA = A*(T*,δ*)/A*(T*,0). Built once per species.
 */
export interface PolarCorrection {
  r22: Float64Array;
  rA: Float64Array;
}

/** Build the polar correction table for reduced dipole moment δ* (dimensionless). */
export function buildPolarCorrection(deltaStar: number): PolarCorrection {
  const r22 = new Float64Array(NT);
  const rA = new Float64Array(NT);
  for (let i = 0; i < NT; i++) {
    r22[i] = interpDelta(MM_OMEGA22, i, deltaStar) / MM_OMEGA22[i * ND];
    rA[i] = interpDelta(MM_ASTAR, i, deltaStar) / MM_ASTAR[i * ND];
  }
  return { r22, rA };
}

/**
 * Quadratic interpolation of a MM_TSTAR-gridded table in ln T* (three nearest rows, as in
 * Cantera's MMCollisionInt). Outside [0.1, 100] the end value is held.
 */
export function interpTstar(table: Float64Array, Ts: number): number {
  if (Ts <= MM_TSTAR[0]) return table[0];
  if (Ts >= MM_TSTAR[NT - 1]) return table[NT - 1];
  let i = 1;
  while (MM_TSTAR[i] < Ts) i++; // MM_TSTAR[i-1] < Ts <= MM_TSTAR[i]
  let i0 = i - 1;
  if (i0 > NT - 3) i0 = NT - 3;
  const x = Math.log(Ts);
  const x0 = LN_TSTAR[i0], x1 = LN_TSTAR[i0 + 1], x2 = LN_TSTAR[i0 + 2];
  const y0 = table[i0], y1 = table[i0 + 1], y2 = table[i0 + 2];
  return (
    (y0 * (x - x1) * (x - x2)) / ((x0 - x1) * (x0 - x2)) +
    (y1 * (x - x0) * (x - x2)) / ((x1 - x0) * (x1 - x2)) +
    (y2 * (x - x0) * (x - x1)) / ((x2 - x0) * (x2 - x1))
  );
}
