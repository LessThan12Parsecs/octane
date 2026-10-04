/**
 * Two compressible restrictions in series between two gas volumes a and b, with a massless
 * intermediate node: a — [restriction A] — m — [restriction B] — b (carburettor venturi + butterfly
 * throttle; valve curtain + L-head pocket transfer).
 *
 * ── Model ───────────────────────────────────────────────────────────────────────────────
 * Each restriction is the quasi-steady isentropic nozzle of orifice.ts (subsonic or choked, low-Δp
 * regularisation). The jet of the first restriction is fully dissipated in the node (no pressure
 * recovery), so the node is a stagnation state at the intermediate pressure p_m; the flow is
 * adiabatic and frozen, so the node's stagnation temperature (and R, γ) are those of the upstream
 * volume. With the flow from u (upstream) to d (downstream) through restrictions 1 then 2,
 *   ṁ₁(p_u, p_m) = ṁ₂(p_m, p_d),   p_d ≤ p_m ≤ p_u,
 * f(p_m) = ṁ₁ − ṁ₂ is strictly decreasing (∂ṁ₁/∂p_m ≤ 0 < ∂ṁ₂/∂p_m), so the root is unique; it is
 * found by Newton's method safeguarded by bisection on the bracket [p_d, p_u], warm-started from the
 * previous call's node position (stored in the result). Every regime is covered without special
 * cases: restriction 2 choked (an idling throttle: ṁ₂ = C_D A₂ p_m Ψ* / √(RT), independent of p_d),
 * restriction 1 choked (ṁ₁ constant), both choked, and the regularised band near Δp = 0.
 * The flow derivatives follow from the implicit-function theorem. With a₁ = ∂ṁ₁/∂p_u,
 * b₁ = ∂ṁ₁/∂p_m, a₂ = ∂ṁ₂/∂p_m, b₂ = ∂ṁ₂/∂p_d:
 *   dṁ/dp_u = a₁ a₂ / (a₂ − b₁),   dṁ/dp_d = b₁ b₂ / (b₁ − a₂)
 * (series conductances; dṁ/dp_d = 0 as soon as either restriction chokes). They feed the cycle
 * model's stiffness bound exactly like orificeFlow's. The reported flow is that of the restriction
 * with the larger share of the drop (the other's small fractional drop carries the rounding of p_m).
 * Flow reverses with the pressure difference (b → a passes B, then A) and then carries side b's state.
 * Oracle: tools/reference/gasex_series_orifice.py (scipy brentq on the same two-nozzle problem,
 * derivatives by finite differences of the full solve).
 *
 * Hot paths allocate nothing (module-level scratch; not re-entrant — single worker thread).
 */
import { newOrificeFlow, ORIFICE_REG_DEFAULT, orificeFlow, type OrificeFlow } from './orifice';

/** Result of {@link seriesOrificeFlow}: an {@link OrificeFlow} (drop-in) plus the node state. */
export interface SeriesOrificeFlow extends OrificeFlow {
  /** Intermediate (node) pressure, Pa. */
  intermediatePressure: number;
  /**
   * Node position (p_m − p_d)/(p_u − p_d) ∈ [0, 1] of the last solve — the warm start of the next
   * call with this result object (NaN = none).
   */
  intermediateFraction: number;
  /** Restriction A (adjacent to side a) choked. */
  chokedA: boolean;
  /** Restriction B (adjacent to side b) choked. */
  chokedB: boolean;
  /** Isentropic jet speed in restriction A, m/s. */
  jetVelocityA: number;
  /** Isentropic jet speed in restriction B, m/s. */
  jetVelocityB: number;
  /** Orifice-pair evaluations of the last solve (1 = the warm start was already converged). */
  iterations: number;
}

/** Allocate a {@link SeriesOrificeFlow} result holder (reuse it: it carries the warm start). */
export function newSeriesOrificeFlow(): SeriesOrificeFlow {
  return {
    ...newOrificeFlow(),
    intermediatePressure: 0,
    intermediateFraction: Number.NaN,
    chokedA: false,
    chokedB: false,
    jetVelocityA: 0,
    jetVelocityB: 0,
    iterations: 0,
  };
}

const S1 = newOrificeFlow();
const S2 = newOrificeFlow();
/** Relative convergence tolerance on the flow balance |ṁ₁ − ṁ₂|/ṁ. */
const TOL = 1e-13;
/** Newton-step tolerance as a fraction of the total pressure difference. */
const STEP_TOL = 1e-11;
const MAX_IT = 60;

/**
 * Bidirectional flow through restriction A (adjacent to side a) and restriction B (adjacent to side
 * b) in series (see file header). Writes `out` and returns out.mdot (kg/s, + = a → b). The caller
 * carries the enthalpy/composition flux with the upstream side's state, as for orificeFlow.
 *
 * @param cdaA effective area of restriction A, m² (∞ = no restriction: single orifice B)
 * @param cdaB effective area of restriction B, m² (∞ = no restriction: single orifice A)
 * @param pa, Ta, Ra, gammaA side-a pressure (Pa), temperature (K), gas constant (J/(kg K)), γ
 * @param pb, Tb, Rb, gammaB side-b state
 * @param out result holder (also carries the warm start)
 * @param delta orifice regularisation half-width (fraction of each restriction's upstream pressure)
 */
export function seriesOrificeFlow(
  cdaA: number,
  cdaB: number,
  pa: number,
  Ta: number,
  Ra: number,
  gammaA: number,
  pb: number,
  Tb: number,
  Rb: number,
  gammaB: number,
  out: SeriesOrificeFlow,
  delta = ORIFICE_REG_DEFAULT,
): number {
  const aUp = pa >= pb;
  out.upstream = aUp ? 0 : 1;
  out.iterations = 0;
  if (!(cdaA > 0 && cdaB > 0) || !((aUp ? pa : pb) > 0)) return zero(out, aUp ? pa : pb);
  // single restriction when the other is absent (∞)
  if (cdaB === Infinity || cdaA === Infinity) {
    const single = cdaB === Infinity ? cdaA : cdaB;
    if (single === Infinity) throw new RangeError('seriesOrificeFlow: both restrictions infinite');
    const m = orificeFlow(single, pa, Ta, Ra, gammaA, pb, Tb, Rb, gammaB, out, delta);
    out.intermediatePressure = cdaB === Infinity ? pb : pa;
    out.intermediateFraction = Number.NaN;
    out.chokedA = cdaB === Infinity ? out.choked : false;
    out.chokedB = cdaB === Infinity ? false : out.choked;
    out.jetVelocityA = cdaB === Infinity ? out.jetVelocity : 0;
    out.jetVelocityB = cdaB === Infinity ? 0 : out.jetVelocity;
    return m;
  }
  const pu = aUp ? pa : pb;
  const pd = aUp ? pb : pa;
  const Tu = aUp ? Ta : Tb;
  const Ru = aUp ? Ra : Rb;
  const gu = aUp ? gammaA : gammaB;
  const c1 = aUp ? cdaA : cdaB; // first restriction along the flow
  const c2 = aUp ? cdaB : cdaA;
  const dp = pu - pd;
  // start: warm start, else the incompressible split Δp₁/Δp = c₂²/(c₁² + c₂²)
  const s0 = out.intermediateFraction;
  let p = s0 > 0 && s0 < 1 ? pd + s0 * dp : pd + (dp * (c1 * c1)) / (c1 * c1 + c2 * c2);
  // Newton steps below STEP_TOL·Δp end the solve (one more evaluation at the new point makes the
  // result exact to rounding: quadratic convergence); the flow balance itself is only resolvable to
  // ≈ ulp(p)/Δp₁ when one restriction takes a tiny share of the drop.
  const stepTol = STEP_TOL * dp > 1e-15 * pu ? STEP_TOL * dp : 1e-15 * pu;
  let lo = pd;
  let hi = pu;
  let m1 = 0;
  let m2 = 0;
  let evals = 0;
  for (let it = 0; it < MAX_IT; it++) {
    m1 = orificeFlow(c1, pu, Tu, Ru, gu, p, Tu, Ru, gu, S1, delta);
    m2 = orificeFlow(c2, p, Tu, Ru, gu, pd, Tu, Ru, gu, S2, delta);
    evals++;
    const f = m1 - m2;
    if (!(Math.abs(f) > TOL * (m1 + m2))) break;
    if (f > 0) lo = p;
    else hi = p;
    const fp = S1.dmdotdpb - S2.dmdotdpa; // b₁ − a₂ < 0
    const step = fp < 0 ? -f / fp : Number.NaN;
    if (Math.abs(step) <= stepTol) {
      // converged (the step may be below the resolution of p: test it before the bracket)
      const q = p + step;
      p = q < pd ? pd : q > pu ? pu : q;
      m1 = orificeFlow(c1, pu, Tu, Ru, gu, p, Tu, Ru, gu, S1, delta);
      m2 = orificeFlow(c2, p, Tu, Ru, gu, pd, Tu, Ru, gu, S2, delta);
      evals++;
      break;
    }
    let pn = p + step;
    if (!(pn > lo && pn < hi)) pn = 0.5 * (lo + hi); // Newton left the bracket (or fp ≥ 0): bisect
    const done = !(Math.abs(pn - p) > stepTol);
    p = pn;
    if (done) {
      m1 = orificeFlow(c1, pu, Tu, Ru, gu, p, Tu, Ru, gu, S1, delta);
      m2 = orificeFlow(c2, p, Tu, Ru, gu, pd, Tu, Ru, gu, S2, delta);
      evals++;
      break;
    }
  }
  out.iterations = evals;
  // the restriction taking the larger share of the drop resolves the flow best (the other one's
  // fractional pressure drop carries the rounding of p)
  const m = pu - p < p - pd ? m2 : m1;
  const a1 = S1.dmdotdpa;
  const b1 = S1.dmdotdpb;
  const a2 = S2.dmdotdpa;
  const b2 = S2.dmdotdpb;
  const den = a2 - b1; // > 0
  // (0 − x: no negative zeros in the derivatives)
  const dUp = den > 0 ? (a1 * a2) / den : 0;
  const dDown = den > 0 ? 0 - (b1 * b2) / den : 0;
  out.intermediatePressure = p;
  out.intermediateFraction = dp > 0 ? (p - pd) / dp : Number.NaN;
  out.pressureDrop = pd > 0 ? 1 - pd / pu : 1;
  out.choked = S1.choked || S2.choked;
  out.chokedA = aUp ? S1.choked : S2.choked;
  out.chokedB = aUp ? S2.choked : S1.choked;
  out.jetVelocityA = aUp ? S1.jetVelocity : S2.jetVelocity;
  out.jetVelocityB = aUp ? S2.jetVelocity : S1.jetVelocity;
  out.jetVelocity = S2.jetVelocity; // the jet entering the downstream volume
  if (aUp) {
    out.mdot = m;
    out.dmdotdpa = dUp;
    out.dmdotdpb = dDown;
  } else {
    out.mdot = -m;
    out.dmdotdpa = 0 - dDown;
    out.dmdotdpb = 0 - dUp;
  }
  return out.mdot;
}

function zero(out: SeriesOrificeFlow, p: number): number {
  out.mdot = 0;
  out.choked = false;
  out.pressureDrop = 0;
  out.jetVelocity = 0;
  out.dmdotdpa = 0;
  out.dmdotdpb = 0;
  out.intermediatePressure = p;
  out.chokedA = false;
  out.chokedB = false;
  out.jetVelocityA = 0;
  out.jetVelocityB = 0;
  return 0;
}
