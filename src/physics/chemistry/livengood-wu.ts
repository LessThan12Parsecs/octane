/**
 * Livengood–Wu autoignition integral (Livengood & Wu 1955, "Correlation of autoignition
 * phenomena in internal combustion engines and rapid compression machines", Proc. Combust.
 * Inst. 5:347):  I(t) = ∫₀ᵗ dt'/τ(T(t'), p(t'), …);  end-gas autoignition when I = 1.
 *
 * Quadrature: τ depends exponentially on the end-gas temperature, so over a step 1/τ is
 * interpolated EXPONENTIALLY (ln(1/τ) linear in t), whose exact integral is the logarithmic
 * mean  Δt (r1 − r0)/ln(r1/r0)  (r = 1/τ). It is exact for Arrhenius-like variation along a
 * smooth T(t), never negative, bounded by the trapezoid (which over-estimates convex growth),
 * and reduces to the trapezoid when r1 ≈ r0. A 'trapezoidal' mode is available for comparison.
 * The crossing time of I = 1 inside the step is found from the same interpolant, so the knock
 * onset is located to sub-step accuracy.
 *
 * LIMITATIONS (two-stage ignition): the integral assumes that the "fraction of the delay
 * consumed" adds linearly under changing conditions — true for a single-stage, induction-
 * controlled (chain-branching) process. PRF mixtures at 650–900 K ignite in two stages: the
 * first (cool-flame) stage releases heat and radicals that change the remaining delay, and
 * in the NTC range τ can even increase with T. A single integral of the total τ then
 * mis-weights the history (typically predicting knock too early after slow compressions
 * through the low-T region, and missing the pressure/temperature jump of the cool flame).
 * TwoStageLivengoodWu integrates the first-stage delay τ1 first and then the remaining
 * second-stage delay τ − τ1 (UNVERIFIED as a literature citation: this two-integral form is in
 * the spirit of two-stage knock models such as Hoepke et al. 2012, not checked). Neither
 * variant models the cool-flame heat release itself — that must come from the cycle model.
 */

/** Quadrature of 1/τ over a step. */
export type LivengoodWuQuadrature = 'logarithmic' | 'trapezoidal';

/**
 * Integral of r(t) over a step of length dt with end values r0, r1 ≥ 0 (1/s), and the
 * fraction f ∈ [0, 1] of the step at which the integral first reaches `need` (NaN if not).
 * Writes [integral, fraction] into out.
 */
function stepIntegral(r0: number, r1: number, dt: number, need: number, log: boolean, out: Float64Array): void {
  let I: number;
  const exp = log && r0 > 0 && r1 > 0 && Number.isFinite(r0) && Number.isFinite(r1);
  const b = exp ? Math.log(r1 / r0) : 0;
  if (exp && Math.abs(b) > 1e-8) I = (dt * (r1 - r0)) / b;
  else I = 0.5 * dt * (r0 + r1);
  out[0] = I;
  out[1] = NaN;
  if (!(I >= need)) return;
  if (!Number.isFinite(I)) {
    out[1] = 0;
    return;
  }
  let f: number;
  if (exp && Math.abs(b) > 1e-8) {
    // ∫₀^{f dt} r0 e^{b s/dt} ds = r0 dt (e^{b f} − 1)/b = need
    f = Math.log1p((need * b) / (r0 * dt)) / b;
  } else {
    // r0 f dt + (r1 − r0) f² dt / 2 = need
    const a = 0.5 * (r1 - r0) * dt;
    const bb = r0 * dt;
    if (Math.abs(a) < 1e-14 * Math.abs(bb)) f = need / bb;
    else f = (2 * need) / (bb + Math.sqrt(Math.max(0, bb * bb + 4 * a * need)));
  }
  out[1] = Math.min(1, Math.max(0, f));
}

/**
 * Livengood–Wu integrator. Call reset() at the start of each cycle (e.g. at IVC or spark),
 * then advance(dt, τ) once per step with τ (s) evaluated at the END of the step; the value at
 * the start of the step is remembered from the previous call (the first step after reset uses
 * τ at its end for both ends unless reset(τ0) is given). τ = Infinity is allowed (no
 * reactivity); τ ≤ 0 means immediate autoignition; τ = NaN yields a NaN integral. Allocation-free.
 */
export class LivengoodWuIntegrator {
  /** ∫dt/τ since reset (dimensionless). */
  integral = 0;
  /** Time since reset, s. */
  time = 0;
  /** Time since reset at which the integral reached 1, s (NaN before autoignition). */
  ignitionTime = NaN;
  readonly quadrature: LivengoodWuQuadrature;
  private rPrev = NaN;
  private readonly tmp = new Float64Array(2);

  constructor(quadrature: LivengoodWuQuadrature = 'logarithmic') {
    this.quadrature = quadrature;
  }

  /** Restart the integral; optionally give τ (s) at the start of the first step. */
  reset(tau0?: number): void {
    this.integral = 0;
    this.time = 0;
    this.ignitionTime = NaN;
    this.rPrev = tau0 === undefined ? NaN : rate(tau0);
  }

  /** True once the integral has reached 1. */
  get autoignited(): boolean {
    return this.integral >= 1;
  }

  /**
   * Advance by dt (s) with ignition delay tau (s) at the end of the step; returns the integral.
   * A NaN tau makes the integral NaN until the next reset() (never a spurious autoignition).
   */
  advance(dt: number, tau: number): number {
    const r1 = rate(tau);
    const r0 = Number.isNaN(this.rPrev) ? r1 : this.rPrev;
    this.rPrev = r1;
    if (dt > 0) {
      const need = 1 - this.integral;
      stepIntegral(r0, r1, dt, need, this.quadrature === 'logarithmic', this.tmp);
      if (need > 0 && this.tmp[0] >= need && Number.isNaN(this.ignitionTime)) {
        this.ignitionTime = this.time + this.tmp[1] * dt;
      }
      this.integral += this.tmp[0];
      this.time += dt;
    }
    return this.integral;
  }
}

/**
 * 1/τ (1/s): τ > 0 → 1/τ; τ = +∞ → 0 (no reactivity); τ ≤ 0 → +∞ (immediate autoignition);
 * τ = NaN (e.g. a PRF-only model asked about another fuel, or an unset state) → NaN, which
 * propagates into the integral so the failure is visible — it must NOT be read as "τ ≤ 0"
 * (that silently fired knock at the first step).
 */
const rate = (tau: number): number => (tau > 0 ? 1 / tau : tau <= 0 ? Infinity : NaN);

/**
 * Two-stage Livengood–Wu: I1 = ∫dt/τ1 until I1 = 1 (first-stage / cool-flame onset), then
 * I2 = ∫dt/(τ − τ1) until I2 = 1 (hot ignition). τ1 = τ for single-stage conditions, in which
 * case the second integral uses the remaining (τ − τ1 → 0 ⇒ immediate) — so for single-stage
 * data pass τ1 < τ only when a distinct first stage exists. See the file header (UNVERIFIED form).
 */
export class TwoStageLivengoodWu {
  readonly first: LivengoodWuIntegrator;
  readonly second: LivengoodWuIntegrator;
  constructor(quadrature: LivengoodWuQuadrature = 'logarithmic') {
    this.first = new LivengoodWuIntegrator(quadrature);
    this.second = new LivengoodWuIntegrator(quadrature);
  }
  reset(): void {
    this.first.reset();
    this.second.reset();
  }
  /** True once the second-stage integral has reached 1. */
  get autoignited(): boolean {
    return this.second.integral >= 1;
  }
  /** Progress measure: I1 before first-stage ignition, 1 + I2 after (reaches 2 at hot ignition). */
  get progress(): number {
    // written so that a NaN first integral reports NaN (not "first stage done")
    return this.first.integral >= 1 ? 1 + Math.min(this.second.integral, 1) : this.first.integral;
  }
  /**
   * Advance by dt (s) with total delay tau and first-stage delay tau1 (s) at the end of the
   * step. Returns progress.
   */
  advance(dt: number, tau: number, tau1: number): number {
    if (!(this.first.integral >= 1)) {
      this.first.advance(dt, tau1);
      if (this.first.integral >= 1) {
        // restart the second integral from the first-stage crossing time
        const tRest = this.first.time - this.first.ignitionTime;
        this.second.reset(Math.max(tau - tau1, 0));
        if (tRest > 0) this.second.advance(tRest, Math.max(tau - tau1, 0));
      }
    } else {
      this.second.advance(dt, Math.max(tau - tau1, 0));
    }
    return this.progress;
  }
}
