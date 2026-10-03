/**
 * Tiny dense linear algebra for the equilibrium Newton / derivative systems.
 *
 * The systems are at most (N_ELEMENTS + 1) × (N_ELEMENTS + 1) = 6 × 6, stored row-major in a
 * Float64Array with a fixed row stride of DENSE_STRIDE. Everything works in place on caller
 * scratch buffers — no allocation.
 *
 * The equilibrium matrices are symmetric but indefinite (the ln n row of the TP system has a
 * ~zero diagonal) and badly scaled (element amounts can differ by many decades), so the solve
 * is: symmetric Jacobi scaling D M D (D = diag(1/√|M_ii|), caller-supplied weights for rows
 * whose diagonal is not representative) → LU with partial pivoting → unscale.
 */

/** Row stride (and maximum dimension) of the dense matrices. */
export const DENSE_STRIDE = 8;

/**
 * In-place LU factorisation with partial (row) pivoting of the dim×dim row-major matrix A
 * (stride DENSE_STRIDE). On return A holds L (unit lower, below the diagonal) and U; piv[i] is
 * the original row now at position i. Returns false if a pivot is exactly zero or not finite
 * (numerically singular); the factors are then unusable.
 */
export function luFactor(A: Float64Array, dim: number, piv: Int32Array): boolean {
  const S = DENSE_STRIDE;
  for (let i = 0; i < dim; i++) piv[i] = i;
  for (let k = 0; k < dim; k++) {
    // pivot search in column k
    let p = k;
    let best = Math.abs(A[k * S + k]);
    for (let i = k + 1; i < dim; i++) {
      const v = Math.abs(A[i * S + k]);
      if (v > best) {
        best = v;
        p = i;
      }
    }
    if (!(best > 0) || !Number.isFinite(best)) return false;
    if (p !== k) {
      for (let j = 0; j < dim; j++) {
        const t = A[k * S + j];
        A[k * S + j] = A[p * S + j];
        A[p * S + j] = t;
      }
      const t = piv[k];
      piv[k] = piv[p];
      piv[p] = t;
    }
    const inv = 1 / A[k * S + k];
    for (let i = k + 1; i < dim; i++) {
      const f = A[i * S + k] * inv;
      if (f === 0) {
        A[i * S + k] = 0;
        continue;
      }
      A[i * S + k] = f;
      for (let j = k + 1; j < dim; j++) A[i * S + j] -= f * A[k * S + j];
    }
  }
  return true;
}

/**
 * Solve (LU) x = P b with factors from luFactor. `rhs` (length ≥ dim) is read; the solution is
 * written into `x` (must not alias rhs).
 */
export function luSolve(A: Float64Array, dim: number, piv: Int32Array, rhs: Float64Array, x: Float64Array): void {
  const S = DENSE_STRIDE;
  // forward substitution with the permutation
  for (let i = 0; i < dim; i++) {
    let s = rhs[piv[i]];
    for (let j = 0; j < i; j++) s -= A[i * S + j] * x[j];
    x[i] = s;
  }
  // back substitution
  for (let i = dim - 1; i >= 0; i--) {
    let s = x[i];
    for (let j = i + 1; j < dim; j++) s -= A[i * S + j] * x[j];
    x[i] = s / A[i * S + i];
  }
}

/**
 * Symmetric diagonal scaling in place: A ← D A D with D = diag(d). Use with luFactor, then
 * solve with scaled right-hand sides (D r) and unscale the solution (x = D y).
 */
export function scaleSymmetric(A: Float64Array, dim: number, d: Float64Array): void {
  const S = DENSE_STRIDE;
  for (let i = 0; i < dim; i++) {
    const di = d[i];
    for (let j = 0; j < dim; j++) A[i * S + j] *= di * d[j];
  }
}
