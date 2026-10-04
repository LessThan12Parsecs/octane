/**
 * Inductive (transistorised / Kettering) ignition circuit as coupled ODEs. SI units.
 *
 * Circuit (Heywood 1988 §9.5.2, Fig. 9-47/9-48; generic transformer model):
 *
 *   V_s ─ R₁ ─ L₁ (primary) ─ switch ∥ C₁ ∥ clamp ─ ground
 *                 ‖ M = k√(L₁L₂)
 *               L₂ (secondary) ─ R₂ ─ node V₂ ─ C₂ ∥ spark gap ─ ground
 *
 * State x = [I₁, I₂, V₁, V₂] (primary current, secondary current, primary-switch voltage,
 * secondary-capacitor voltage). Currents enter the dotted winding terminals:
 *   L₁ İ₁ + M İ₂ = V_s − R₁ I₁ − v_sw
 *   M İ₁ + L₂ İ₂ = V₂ − R₂ I₂
 *   C₁ V̇₁ = I₁ − I_clamp                          (switch open; V₁ ≡ 0 while closed)
 *   C₂ V̇₂ = −I₂ − I_gap
 * With these signs the stored energy is W = ½L₁I₁² + M I₁I₂ + ½L₂I₂² + ½C₁V₁² + ½C₂V₂² and
 *   Ẇ = V_s I₁ − R₁I₁² − R₂I₂² − v_drv I₁ − V_clamp I_clamp − V₂ I_gap.
 * The analytic dwell limit (R₂→∞ secondary) is Heywood eq. 9.53: I₁ = (V_s/R₁)(1 − e^{−R₁t/L₁}).
 *
 * Switch model:
 *  - dwell (closed): V₁ = 0; if the driver current limit I_lim is reached the switch
 *    drops the voltage v_drv needed to hold I₁ = I_lim (linear-mode IGBT; loss v_drv I₁);
 *  - open: primary current flows into the lumped primary-node capacitance C₁ (driver output
 *    + winding capacitance, or the condenser of a breaker system); when V₁ reaches the
 *    driver clamp voltage the clamp conducts (ignition IGBTs have an internal collector
 *    clamp, e.g. onsemi ISL9V3040: BV_CER = 400 V, "internal diodes provide voltage
 *    clamping"). The clamp limits the open-circuit secondary voltage to ≈ n·V_clamp. A
 *    reverse clamp (emitter–collector junction) limits negative collector swings.
 *
 * Gap (discharge.ts): open → C₂ dynamics with I_gap = 0 until |V₂| reaches the breakdown
 * voltage; conducting → V₂ is algebraic, V₂ = −sgn(I₂)·V_gap(|I₂|), because the discharge
 * time scale (ms) is ≫ the C₂–gap time scale (C₂ dV/dt ≈ 60 pF·500 V/ms = 3e-8 A ≪ I₂);
 * the gap current is I_gap = −I₂ − C₂V̇₂ and its energy is booked exactly as
 * −∫V₂I₂dt − ΔW_C₂. After an extinction the gap needs a fresh breakdown: the re-ignition of a
 * recovering channel (SparkGap.recovering / reignite, discharge.ts) is a trembler-coil event only —
 * this single inductive spark ends after a monotonic ms-long current decay.
 *
 * Integration: trapezoidal rule (A-stable, second order; for this linear circuit it
 * conserves the quadratic energy exactly: W¹ − W⁰ = h·P(x̄) with x̄ the step midpoint, so
 * the ledger closes to round-off). Conducting sub-steps linearise V_gap(I) about the
 * iterate (Newton, ≤ 3 iterations). Events (clamp on/off, breakdown, extinction) are
 * located by linear interpolation inside the sub-step and the sub-step is re-taken up to
 * the event. Sub-step sizes are mode dependent (0.1 µs during the ring-up, 1 µs while
 * conducting, 2 µs during the dwell, 5 µs in the post-discharge ringing tail); a quiescent
 * circuit is skipped entirely.
 *
 * Extensions beyond IgnitionSystemSpec (proposed contract additions): primary-node
 * capacitance C₁, forward and reverse clamp voltages (CoilOptions). Events are located as
 * the earliest of the sub-step. After two (near-)zero-length event sub-steps in a row (clamp
 * releasing and re-engaging at one instant) the clamp events are ignored for a sub-step,
 * after four all events (`forcedSubsteps` counts these; energy stays exact).
 */

import type { InductiveIgnitionSpec } from '../core/engine-spec';
import { GapTransition, type SparkGap } from './discharge';

/** Primary-switch state. */
export type PrimarySwitchState = 'closed' | 'open' | 'clamped';

/** Extra circuit constants not (yet) in IgnitionSystemSpec, plus integrator controls. */
export interface CoilOptions {
  /**
   * Lumped capacitance across the primary switch after it opens, F.
   * UNVERIFIED default 1 nF (IGBT output + winding capacitance; breaker-point systems use a
   * 0.2–0.3 µF condenser instead).
   */
  primaryCapacitance: number;
  /** Driver collector clamp voltage, V (default 400 V: ISL9V3040 BV_CER; Infinity = none). */
  primaryClampVoltage: number;
  /**
   * Reverse (negative-collector) clamp magnitude, V: the emitter–collector junction of an
   * ignition IGBT conducts when the collector swings negative. UNVERIFIED default 24 V
   * (typical BV_ECS rating of ignition IGBTs; Infinity = none).
   */
  primaryReverseClampVoltage: number;
  /** Max sub-step during the high-voltage ring-up (switch open, gap open), s. */
  ringUpSubstep: number;
  /** Max sub-step while the gap conducts, s. */
  conductingSubstep: number;
  /** Max sub-step during the dwell, s. */
  slowSubstep: number;
  /** Max sub-step of the post-discharge ringing tail (switch open, gap open, far from breakdown), s. */
  tailSubstep: number;
}

/** Default CoilOptions (see field docs). */
export const DEFAULT_COIL_OPTIONS: Readonly<CoilOptions> = Object.freeze({
  primaryCapacitance: 1e-9,
  primaryClampVoltage: 400,
  primaryReverseClampVoltage: 24,
  ringUpSubstep: 0.1e-6,
  conductingSubstep: 1e-6,
  slowSubstep: 2e-6,
  tailSubstep: 5e-6,
});

/**
 * Event sub-steps that advance less than this (s) count as "stalled" for the zero-length
 * event-loop guard (1 ps ≪ every circuit time scale; breakdown itself is instantaneous here).
 * Shared with the trembler coil (trembler-coil.ts).
 */
export const MIN_EVENT_PROGRESS = 1e-12;
/**
 * Residual |V₂| amplitude below which the open, ringing circuit is snapped to rest, V (see
 * settleIfQuiescent). Shared with the trembler coil (trembler-coil.ts).
 */
export const QUIESCENT_V2 = 25;

// primary / secondary modes (numeric for speed)
const P_CLOSED = 0;
const P_OPEN = 1;
const P_CLAMP = 2;
const S_OPEN = 0;
const S_COND = 1;

/** Solve the dense n×n system a·x = b in place (partial pivoting). a row-major n*n, result in b. */
export function solveDense(a: Float64Array, b: Float64Array, n: number): void {
  for (let c = 0; c < n; c++) {
    let piv = c;
    let best = Math.abs(a[c * n + c]);
    for (let r = c + 1; r < n; r++) {
      const v = Math.abs(a[r * n + c]);
      if (v > best) {
        best = v;
        piv = r;
      }
    }
    if (piv !== c) {
      for (let k = 0; k < n; k++) {
        const t = a[c * n + k];
        a[c * n + k] = a[piv * n + k];
        a[piv * n + k] = t;
      }
      const t = b[c];
      b[c] = b[piv];
      b[piv] = t;
    }
    const d = a[c * n + c];
    for (let r = c + 1; r < n; r++) {
      const f = a[r * n + c] / d;
      if (f === 0) continue;
      for (let k = c; k < n; k++) a[r * n + k] -= f * a[c * n + k];
      b[r] -= f * b[c];
    }
  }
  for (let r = n - 1; r >= 0; r--) {
    let s = b[r];
    for (let k = r + 1; k < n; k++) s -= a[r * n + k] * b[k];
    b[r] = s / a[r * n + r];
  }
}

/**
 * Coupled primary/secondary ignition-coil circuit with a spark-gap load.
 * All fields SI (A, V, J, s). Allocation-free `step`.
 */
export class IgnitionCoil {
  readonly spec: InductiveIgnitionSpec;
  readonly opts: CoilOptions;
  /** Mutual inductance M = k√(L₁L₂), H. */
  readonly mutualInductance: number;

  // ---- circuit state ----
  /** Primary current I₁, A. */
  I1 = 0;
  /** Secondary current I₂ (into the dotted terminal), A. */
  I2 = 0;
  /** Primary switch (collector) voltage V₁, V. */
  V1 = 0;
  /** Secondary-capacitor (plug) voltage V₂, V. */
  V2 = 0;
  /** Internal clock, s. */
  time = 0;
  /** Primary switch state. */
  switchState: PrimarySwitchState = 'open';
  /** True while the driver is current limiting. */
  currentLimited = false;

  // ---- energy ledger (since construction / resetLedger), J ----
  /** Energy drawn from the supply. */
  energySupplied = 0;
  /** Dissipated in R₁. */
  energyPrimaryResistance = 0;
  /** Dissipated in R₂. */
  energySecondaryResistance = 0;
  /** Dissipated in the driver (current limiting + C₁ discharge on switch closing). */
  energyDriver = 0;
  /** Dissipated in the collector clamp. */
  energyClamp = 0;
  /** Delivered to the spark gap. */
  energyGap = 0;
  /** Residual ringing energy removed when the quiescent circuit is snapped to rest (≤ ½C₂·(25 V)²). */
  energyTailTruncated = 0;
  /** Stored energy at the last resetLedger, J. */
  energyStoredRef = 0;
  /** Primary magnetic energy ½L₁I₁² at the last switch-off, J. */
  energyAtSwitchOff = 0;
  /** Diagnostic: sub-steps taken without event handling to break zero-length event loops. */
  forcedSubsteps = 0;
  /** Peak |V₂| since the last switch-off, V. */
  peakSecondaryVoltage = 0;
  /** Time of the peak |V₂|, s. */
  peakSecondaryVoltageTime = 0;

  private pmode = P_OPEN;
  /** Voltage at which the active clamp holds V₁ (+V_clamp or −V_reverse), V. */
  private clampLevel = 0;
  private smode = S_OPEN;
  /** sgn(I₂) frozen while conducting. */
  private sgn = 1;
  private readonly A = new Float64Array(16);
  private readonly b = new Float64Array(4);
  private readonly x0 = new Float64Array(4);
  private readonly x1 = new Float64Array(4);
  private readonly slope = new Float64Array(1);
  /** Driver voltage of the last limited sub-step, V. */
  private vDrv = 0;
  /** Time of the last primary switch-off, s. */
  private switchOffAt = -Infinity;
  /**
   * Time after switch-off during which the fast primary leakage mode (L₁(1−k²) with C₁,
   * ≈ 400 kHz for the CFR-like coil) must be resolved: 10 e-folding times 2L₁(1−k²)/R₁, s.
   */
  readonly fastModeDecayTime: number;
  /** Sub-step resolving the fast leakage mode, 0.5/ω_f = 0.5√(L₁(1−k²)C₁), s. */
  readonly fastModeSubstep: number;

  constructor(spec: InductiveIgnitionSpec, opts?: Partial<CoilOptions>) {
    this.spec = spec;
    this.opts = { ...DEFAULT_COIL_OPTIONS, ...opts };
    this.mutualInductance = spec.couplingCoefficient * Math.sqrt(spec.primaryInductance * spec.secondaryInductance);
    this.V1 = spec.supplyVoltage; // switch open at rest: collector sits at the supply voltage
    this.energyStoredRef = this.storedEnergy();
    const lLeak = spec.primaryInductance * (1 - spec.couplingCoefficient * spec.couplingCoefficient);
    this.fastModeDecayTime = spec.primaryResistance > 0 ? (20 * lLeak) / spec.primaryResistance : Infinity;
    this.fastModeSubstep = Math.max(this.opts.ringUpSubstep, 0.5 * Math.sqrt(lLeak * this.opts.primaryCapacitance));
  }

  /** Total stored energy W = ½L₁I₁² + M I₁I₂ + ½L₂I₂² + ½C₁V₁² + ½C₂V₂², J. */
  storedEnergy(): number {
    const s = this.spec;
    return (
      0.5 * s.primaryInductance * this.I1 * this.I1 +
      this.mutualInductance * this.I1 * this.I2 +
      0.5 * s.secondaryInductance * this.I2 * this.I2 +
      0.5 * this.opts.primaryCapacitance * this.V1 * this.V1 +
      0.5 * s.secondaryCapacitance * this.V2 * this.V2
    );
  }

  /** Energy-balance residual W − W_ref − (supplied − dissipated − gap), J (≈ 0). */
  energyResidual(): number {
    return (
      this.storedEnergy() -
      this.energyStoredRef -
      (this.energySupplied -
        this.energyPrimaryResistance -
        this.energySecondaryResistance -
        this.energyDriver -
        this.energyClamp -
        this.energyGap -
        this.energyTailTruncated)
    );
  }

  /** Zero the energy ledger (keeps the circuit state). */
  resetLedger(): void {
    this.energySupplied = 0;
    this.energyPrimaryResistance = 0;
    this.energySecondaryResistance = 0;
    this.energyDriver = 0;
    this.energyClamp = 0;
    this.energyGap = 0;
    this.energyTailTruncated = 0;
    this.energyStoredRef = this.storedEnergy();
  }

  /** True when the secondary is in the conducting (gap-clamped) mode. */
  get secondaryConducting(): boolean {
    return this.smode === S_COND;
  }

  /**
   * Advance the circuit by dt (s). `dwellOn` = driver switch closed (primary charging).
   * Switch edges take effect at the start of the step; the caller splits steps at events.
   */
  step(dt: number, dwellOn: boolean, gap: SparkGap): void {
    // ---- switch edges ----
    if (dwellOn && this.pmode !== P_CLOSED) {
      // Closing the switch discharges C₁ through it.
      this.energyDriver += 0.5 * this.opts.primaryCapacitance * this.V1 * this.V1;
      this.V1 = 0;
      this.pmode = P_CLOSED;
      this.switchState = 'closed';
    } else if (!dwellOn && this.pmode === P_CLOSED) {
      this.pmode = P_OPEN;
      this.switchState = 'open';
      this.currentLimited = false;
      this.energyAtSwitchOff = 0.5 * this.spec.primaryInductance * this.I1 * this.I1;
      this.switchOffAt = this.time;
      this.peakSecondaryVoltage = 0;
      this.peakSecondaryVoltageTime = this.time;
    }
    // ---- quiescent fast path (switch open, gap open, residual ringing negligible) ----
    if (this.pmode === P_OPEN && this.smode === S_OPEN && !gap.awaitingReentry && this.settleIfQuiescent()) {
      this.time += dt;
      return;
    }
    let remaining = dt;
    let stalls = 0;
    while (remaining > 1e-15) {
      const hMax = this.maxSubstep(gap);
      const h = remaining < hMax ? remaining : hMax;
      // Two (near-)zero-length event sub-steps in a row (e.g. a clamp releasing and
      // re-engaging at the same instant) → the next sub-steps ignore the primary clamp
      // events (gap events stay active); after four, all events. "Near-zero"
      // (< MIN_EVENT_PROGRESS) matters: an event chain advancing 1e-16 s per iteration is a
      // Zeno loop that never finishes the step (reviewer fix).
      const used = this.substep(h, gap, stalls < 2 ? 2 : stalls < 4 ? 1 : 0);
      if (stalls >= 2) this.forcedSubsteps++;
      if (used > MIN_EVENT_PROGRESS || used >= h) stalls = 0;
      else stalls++;
      remaining -= used;
    }
  }

  /**
   * Energy of the deviation from the rest state (I₁ = I₂ = V₂ = 0, V₁ = V_s), J — the
   * Lyapunov function of the open-switch circuit.
   */
  deviationEnergy(): number {
    const s = this.spec;
    const dv1 = this.V1 - s.supplyVoltage;
    return (
      0.5 * s.primaryInductance * this.I1 * this.I1 +
      this.mutualInductance * this.I1 * this.I2 +
      0.5 * s.secondaryInductance * this.I2 * this.I2 +
      0.5 * this.opts.primaryCapacitance * dv1 * dv1 +
      0.5 * s.secondaryCapacitance * this.V2 * this.V2
    );
  }

  /**
   * True if the open circuit is (or is snapped) at rest. The post-spark ringing of the
   * L₂C₂ tank decays only through R₁/R₂ here (no core loss: e-folding ≈ 2L₂/R₂ = 10 ms for
   * the CFR-like coil), so it would be integrated for ~0.2 s after every spark. Once its
   * deviation energy is below ½C₂(QUIESCENT_V2)² it can never again drive |V₂| above
   * QUIESCENT_V2 = 25 V — an order of magnitude below the Paschen minimum (≈ 270 V) that
   * any breakdown needs — so the state is snapped to rest and the (sub-10-nJ) remainder is
   * booked in `energyTailTruncated` (ledger stays exact). Reviewer change.
   */
  private settleIfQuiescent(): boolean {
    const s = this.spec;
    if (this.I1 === 0 && this.I2 === 0 && this.V2 === 0 && this.V1 === s.supplyVoltage) return true;
    if (this.deviationEnergy() > 0.5 * s.secondaryCapacitance * QUIESCENT_V2 * QUIESCENT_V2) return false;
    const w0 = this.storedEnergy();
    this.I1 = 0;
    this.I2 = 0;
    this.V2 = 0;
    this.V1 = s.supplyVoltage;
    this.energyTailTruncated += w0 - this.storedEnergy();
    return true;
  }

  private maxSubstep(gap: SparkGap): number {
    const o = this.opts;
    if (this.smode === S_COND) {
      // While the primary leakage mode still rings (first half of fastModeDecayTime after
      // switch-off) resolve it (ω_f h ≤ 0.5); trapezoidal steps with ω_f h ≫ 1 keep its energy
      // exact but under-damp it (numerical decay rate × 1/(1 + (ω_f h/2)²)).
      if (this.time - this.switchOffAt < 0.5 * this.fastModeDecayTime && this.fastModeSubstep < o.conductingSubstep)
        return this.fastModeSubstep;
      return o.conductingSubstep;
    }
    if (this.pmode === P_CLAMP) return o.ringUpSubstep;
    if (this.pmode === P_OPEN) {
      // Fast high-voltage ring-up right after switch-off (while the lightly damped primary
      // leakage mode is still ringing) …
      if (
        Math.abs(this.I1) > 1e-3 * this.spec.primaryCurrentLimit &&
        this.time - this.switchOffAt < this.fastModeDecayTime
      )
        return o.ringUpSubstep;
      // … or |V₂| could reach the conduction threshold within a few slow sub-steps
      // (open gap: |dV₂/dt| = |I₂|/C₂): resolve the slow L₂C₂ ring near its crest. Reviewer
      // change: the ring-up sub-step (0.1 µs) was used here for the whole post-spark /
      // no-breakdown ringing tail (tens of ms), costing up to ~10⁵ sub-steps per cycle.
      const thr = gap.threshold(Math.abs(this.I2));
      const reach = Math.abs(this.V2) + (4 * o.tailSubstep * Math.abs(this.I2)) / this.spec.secondaryCapacitance;
      if (reach > 0.9 * thr) return o.conductingSubstep < o.tailSubstep ? o.conductingSubstep : o.tailSubstep;
      return o.tailSubstep;
    }
    return o.slowSubstep;
  }

  /**
   * Trapezoidal solve of one sub-step h from x0 into x1 with the current modes.
   * @param limited solve with I₁¹ = I_lim fixed (unknown 0 becomes the driver voltage)
   */
  private solve(h: number, gap: SparkGap, limited: boolean): void {
    const s = this.spec;
    const L1 = s.primaryInductance;
    const L2 = s.secondaryInductance;
    const M = this.mutualInductance;
    const R1 = s.primaryResistance;
    const R2 = s.secondaryResistance;
    const C1 = this.opts.primaryCapacitance;
    const C2 = s.secondaryCapacitance;
    const x0 = this.x0;
    const x1 = this.x1;
    const A = this.A;
    const b = this.b;
    const hh = 0.5 * h;
    const I10 = x0[0];
    const I20 = x0[1];
    const V10 = x0[2];
    const V20 = x0[3];
    const Ilim = s.primaryCurrentLimit;
    // Newton iterations only needed while conducting (quadratic convergence: ≤ 3 suffice).
    const iters = this.smode === S_COND ? 3 : 1;
    let iStar = Math.abs(I20);
    for (let it = 0; it < iters; it++) {
      A.fill(0);
      // row 0: primary KVL
      if (limited) {
        // unknown 0 = mean driver voltage v̄_drv; I₁¹ = I_lim known
        A[0] = h;
        A[1] = M;
        b[0] = L1 * I10 + M * I20 + h * s.supplyVoltage - hh * R1 * (I10 + Ilim) - L1 * Ilim;
      } else {
        A[0] = L1 + hh * R1;
        A[1] = M;
        let r = L1 * I10 + M * I20 + h * s.supplyVoltage - hh * R1 * I10;
        if (this.pmode === P_OPEN) {
          A[2] = hh;
          r -= hh * V10;
        } else if (this.pmode === P_CLAMP) {
          r -= h * this.clampLevel;
        }
        b[0] = r;
      }
      // row 1: secondary KVL
      if (limited) {
        A[4] = 0;
        b[1] = M * I10 + L2 * I20 + hh * V20 - hh * R2 * I20 - M * Ilim;
      } else {
        A[4] = M;
        b[1] = M * I10 + L2 * I20 + hh * V20 - hh * R2 * I20;
      }
      A[5] = L2 + hh * R2;
      A[7] = -hh;
      // row 2: primary node
      if (this.pmode === P_OPEN) {
        A[8] = limited ? 0 : -hh;
        A[10] = C1;
        b[2] = C1 * V10 + hh * I10;
      } else {
        A[10] = 1;
        b[2] = this.pmode === P_CLAMP ? this.clampLevel : 0;
      }
      // row 3: secondary node
      if (this.smode === S_OPEN) {
        A[13] = hh;
        A[15] = C2;
        b[3] = C2 * V20 - hh * I20;
      } else {
        // V₂ = −sgn·V_g(|I₂|) linearised about i*: V₂ + g' I₂ = −sgn·(V_g(i*) − g' i*)
        const vg = gap.gapVoltageAndSlope(iStar, this.slope);
        const g = this.slope[0];
        A[13] = g;
        A[15] = 1;
        b[3] = -this.sgn * (vg - g * iStar);
      }
      solveDense(A, b, 4);
      if (limited) {
        this.vDrv = b[0];
        x1[0] = Ilim;
      } else {
        x1[0] = b[0];
      }
      x1[1] = b[1];
      x1[2] = b[2];
      x1[3] = b[3];
      if (this.smode === S_COND) {
        const iNew = Math.abs(x1[1]);
        if (Math.abs(iNew - iStar) <= 1e-12 + 1e-9 * iNew) break;
        iStar = iNew;
      }
    }
  }

  /** Solve a sub-step (with current limiting), returning into x1. */
  private solveStep(h: number, gap: SparkGap): boolean {
    this.solve(h, gap, false);
    if (this.pmode === P_CLOSED && this.x1[0] > this.spec.primaryCurrentLimit) {
      this.solve(h, gap, true);
      return true;
    }
    return false;
  }

  /** One sub-step of at most h with event handling; returns the time actually advanced. */
  private substep(h: number, gap: SparkGap, events: number): number {
    const x0 = this.x0;
    const x1 = this.x1;
    x0[0] = this.I1;
    x0[1] = this.I2;
    x0[2] = this.V1;
    x0[3] = this.V2;
    let limited = this.solveStep(h, gap);
    let hUsed = h;
    let event = 0; // 1 clamp-on, 2 clamp-off, 3 conduction start, 4 extinction
    // ---- event location: the EARLIEST event of the sub-step wins ----
    // (reviewer fix: the clamp events used to shadow the gap events in an if/else chain, so
    // clamp chatter at I₁ ≈ 0 plus the forced sub-steps could carry |V₂| kilovolts past the
    // breakdown threshold unchecked.) `events` = 2: all events; 1: gap events only (primary
    // clamp chatter suppressed); 0: none.
    if (events === 0) {
      this.commit(h, gap, limited);
      if (this.smode === S_COND) this.afterConducting(gap);
      return h;
    }
    let thBest = 2;
    let newClamp = 0;
    if (events >= 2) {
      const vFwd = this.opts.primaryClampVoltage;
      const vRev = -this.opts.primaryReverseClampVoltage;
      if (this.pmode === P_OPEN && (x1[2] > vFwd || x1[2] < vRev)) {
        newClamp = x1[2] > vFwd ? vFwd : vRev;
        thBest = (newClamp - x0[2]) / (x1[2] - x0[2]);
        event = 1;
      } else if (this.pmode === P_CLAMP && x1[0] * this.clampLevel <= 0) {
        // clamp current (∝ I₁ in the clamp's conduction direction) reached zero
        thBest = x0[0] / (x0[0] - x1[0]);
        event = 2;
      }
    }
    if (this.smode === S_OPEN) {
      const thr = gap.threshold(Math.abs(x1[1]));
      const v1a = Math.abs(x1[3]);
      if (v1a >= thr) {
        const v0a = Math.abs(x0[3]);
        const th = v1a > v0a ? (thr - v0a) / (v1a - v0a) : 1;
        if (th < thBest) {
          thBest = th;
          event = 3;
        }
      }
    } else {
      const iExt = gap.opts.extinctionCurrent;
      const i1 = x1[1] * this.sgn; // signed along the frozen direction
      if (i1 <= iExt) {
        const i0 = Math.abs(x0[1]);
        const th = i0 > i1 ? (i0 - iExt) / (i0 - i1) : 1;
        if (th < thBest) {
          thBest = th;
          event = 4;
        }
      }
    }
    if (event !== 0) hUsed = h * Math.min(1, Math.max(0, thBest));
    if (event !== 0 && hUsed < h) {
      if (hUsed <= 0) {
        // Event at the very start of the sub-step: x1 = x0.
        x1[0] = x0[0];
        x1[1] = x0[1];
        x1[2] = x0[2];
        x1[3] = x0[3];
        limited = false;
      } else {
        limited = this.solveStep(hUsed, gap);
      }
    }
    this.commit(hUsed, gap, limited);
    // ---- apply events ----
    if (event === 1) {
      // Snap V₁ onto the clamp level; the (tiny) C₁ energy difference goes to the clamp.
      const c1 = this.opts.primaryCapacitance;
      this.energyClamp += 0.5 * c1 * (this.V1 * this.V1 - newClamp * newClamp);
      this.V1 = newClamp;
      this.clampLevel = newClamp;
      this.pmode = P_CLAMP;
      this.switchState = 'clamped';
    } else if (event === 2) {
      this.pmode = P_OPEN;
      this.switchState = 'open';
    } else if (event === 3) {
      this.startConduction(gap);
    } else if (event === 4) {
      gap.extinguish(this.time);
      this.smode = S_OPEN;
    } else if (this.smode === S_COND) {
      this.afterConducting(gap);
    } else if (gap.awaitingReentry && Math.abs(this.I2) <= gap.opts.extinctionCurrent) {
      // arc→glow recharge: the inductor current decayed to the sustaining current before C₂
      // reached the glow voltage → the discharge is over (else `awaitingReentry` could stay
      // set for the rest of the cycle and the kernel would never see the spark end).
      gap.extinguish(this.time);
    }
    return hUsed;
  }

  /** Accept x1 as the new state and book the sub-step energies (exact trapezoidal ledger). */
  private commit(h: number, gap: SparkGap, limited: boolean): void {
    const s = this.spec;
    const x0 = this.x0;
    const x1 = this.x1;
    const i1m = 0.5 * (x0[0] + x1[0]);
    const i2m = 0.5 * (x0[1] + x1[1]);
    this.energySupplied += h * s.supplyVoltage * i1m;
    this.energyPrimaryResistance += h * s.primaryResistance * 0.25 * (x0[0] + x1[0]) * (x0[0] + x1[0]);
    this.energySecondaryResistance += h * s.secondaryResistance * i2m * i2m;
    if (this.pmode === P_CLOSED && limited) this.energyDriver += h * this.vDrv * i1m;
    if (this.pmode === P_CLAMP) this.energyClamp += h * this.clampLevel * i1m;
    if (this.smode === S_COND) {
      const c2 = s.secondaryCapacitance;
      const v2m = 0.5 * (x0[3] + x1[3]);
      const e = -h * v2m * i2m - 0.5 * c2 * (x1[3] * x1[3] - x0[3] * x0[3]);
      this.energyGap += e;
      gap.accumulate(e, Math.abs(x1[3]), Math.abs(x1[1]), h);
    }
    this.currentLimited = limited;
    this.I1 = x1[0];
    this.I2 = x1[1];
    this.V1 = this.pmode === P_CLOSED ? 0 : x1[2];
    this.V2 = x1[3];
    this.time += h;
    const va = Math.abs(this.V2);
    if (this.pmode !== P_CLOSED && va > this.peakSecondaryVoltage) {
      this.peakSecondaryVoltage = va;
      this.peakSecondaryVoltageTime = this.time;
    }
  }

  /** Open gap reached its threshold: start (or resume) conduction. */
  private startConduction(gap: SparkGap): void {
    const c2 = this.spec.secondaryCapacitance;
    const wasReentry = gap.awaitingReentry;
    const vAbs = Math.abs(this.V2);
    // Current direction: that of I₂, or (if ~0) the one that discharges C₂.
    this.sgn = this.I2 > 0 ? 1 : this.I2 < 0 ? -1 : this.V2 < 0 ? 1 : -1;
    const vAfter = gap.startConduction(vAbs, Math.abs(this.I2), c2, this.time);
    if (!gap.conducting) {
      // Capacitive spark only (circuit current below the sustaining current), or a re-entry
      // that found the current extinguished: the capacitor keeps its polarity and the gap
      // stays open (the dump, if any, was booked by the gap).
      const vNew = (this.V2 >= 0 ? 1 : -1) * vAfter;
      this.energyGap += 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
      this.V2 = vNew;
      this.smode = S_OPEN;
      return;
    }
    const vNew = -this.sgn * vAfter;
    // Capacitor energy released (breakdown: ½C₂(V_bd² − V_after²), booked by the gap itself;
    // re-entry: only the tiny interpolation mismatch, booked here to keep the ledger exact).
    const e = 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
    this.energyGap += e;
    if (wasReentry) gap.accumulate(e, vAfter, Math.abs(this.I2), 0);
    this.V2 = vNew;
    this.smode = S_COND;
  }

  /** Post-sub-step discharge transitions (mode changes, extinction, restrike). */
  private afterConducting(gap: SparkGap): void {
    const iAbs = this.I2 * this.sgn;
    const code = gap.afterConductingStep(iAbs > 0 ? iAbs : 0, this.time);
    if (code === GapTransition.None) return;
    if (code === GapTransition.VoltageDrop) {
      const c2 = this.spec.secondaryCapacitance;
      const vNew = -this.sgn * gap.gapVoltage(iAbs);
      const e = 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
      this.energyGap += e;
      gap.dump(e);
      this.V2 = vNew;
      return;
    }
    // Extinguish / Recharge / Restrike: back to the open-gap dynamics with V₂ continuous.
    this.smode = S_OPEN;
  }
}
