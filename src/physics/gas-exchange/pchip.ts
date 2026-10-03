/**
 * Monotone piecewise-cubic Hermite interpolation (PCHIP) of tabulated data:
 * Fritsch & Carlson (1980), SIAM J. Numer. Anal. 17(2):238–246, with the three-point
 * end-slope formula of Moler (2004), "Numerical Computing with MATLAB", §3.6 (pchip).
 * C¹, no overshoot for monotone data. Outside the table the end values are held (flat).
 * Construction allocates; evaluation does not.
 */
export class Pchip {
  private readonly x: Float64Array;
  private readonly y: Float64Array;
  private readonly d: Float64Array;
  private readonly n: number;

  /** @param xs strictly increasing abscissae (≥ 2), @param ys ordinates */
  constructor(xs: readonly number[], ys: readonly number[]) {
    const n = xs.length;
    if (n < 2 || ys.length !== n) throw new RangeError('Pchip: need ≥ 2 points and equal lengths');
    for (let i = 1; i < n; i++) if (!(xs[i] > xs[i - 1])) throw new RangeError('Pchip: x must increase');
    this.n = n;
    this.x = Float64Array.from(xs);
    this.y = Float64Array.from(ys);
    const h = new Float64Array(n - 1);
    const del = new Float64Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
      h[i] = xs[i + 1] - xs[i];
      del[i] = (ys[i + 1] - ys[i]) / h[i];
    }
    const d = new Float64Array(n);
    if (n === 2) {
      d[0] = d[1] = del[0];
    } else {
      for (let k = 1; k < n - 1; k++) {
        if (del[k - 1] * del[k] <= 0) d[k] = 0;
        else {
          const w1 = 2 * h[k] + h[k - 1];
          const w2 = h[k] + 2 * h[k - 1];
          d[k] = (w1 + w2) / (w1 / del[k - 1] + w2 / del[k]);
        }
      }
      d[0] = endSlope(h[0], h[1], del[0], del[1]);
      d[n - 1] = endSlope(h[n - 2], h[n - 3], del[n - 2], del[n - 3]);
    }
    this.d = d;
  }

  /** Interpolated value at t (flat extrapolation). */
  value(t: number): number {
    const x = this.x;
    const n = this.n;
    if (!(t > x[0])) return this.y[0];
    if (!(t < x[n - 1])) return this.y[n - 1];
    let i = 0;
    while (t >= x[i + 1]) i++;
    const h = x[i + 1] - x[i];
    const s = (t - x[i]) / h;
    const s2 = s * s;
    const s3 = s2 * s;
    return (
      (2 * s3 - 3 * s2 + 1) * this.y[i] +
      (s3 - 2 * s2 + s) * h * this.d[i] +
      (-2 * s3 + 3 * s2) * this.y[i + 1] +
      (s3 - s2) * h * this.d[i + 1]
    );
  }

  /** First abscissa. */
  get xMin(): number {
    return this.x[0];
  }
  /** Last abscissa. */
  get xMax(): number {
    return this.x[this.n - 1];
  }
  /** Last ordinate. */
  get yLast(): number {
    return this.y[this.n - 1];
  }
}

function endSlope(h0: number, h1: number, del0: number, del1: number): number {
  let d = ((2 * h0 + h1) * del0 - h0 * del1) / (h0 + h1);
  if (Math.sign(d) !== Math.sign(del0)) d = 0;
  else if (Math.sign(del0) !== Math.sign(del1) && Math.abs(d) > Math.abs(3 * del0)) d = 3 * del0;
  return d;
}
