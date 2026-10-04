/**
 * Trembler (vibrator, "buzz") coil circuit, the Ford Model T ignition coil, as coupled ODEs. SI units.
 *
 * Circuit (one coil connected through its timer segment; Ford Service Par. 996–997, Boggess &
 * Patterson 1999):
 *
 *   e(t) ─ R_s, L_s (supply) ─ vibrator points ∥ C₁ (condenser) ─ R₁, L₁ (primary) ─ timer (R_t) ─ ground
 *                                                     ‖ M = k√(L₁L₂)
 *                                   L₂ (secondary) ─ R₂ ─ node V₂ ─ C₂ ∥ spark gap ─ ground
 *
 * State x = [I₁, I₂, V₁, V₂] (primary loop current, secondary current, condenser = points voltage,
 * secondary-capacitor voltage), as in coil.ts, with the loop self-inductance L_p = L₁ + L_s and loop
 * resistance R = R₁ + R_s + R_t (the supply's L_s couples to nothing: it enters the self-term only):
 *   L_p İ₁ + M İ₂ = e(t) − R I₁ − v_pts            (timer closed; v_pts = V₁ with the points open, 0 closed)
 *   I₁ ≡ 0                                         (timer open)
 *   M İ₁ + L₂ İ₂ = V₂ − R₂ I₂
 *   C₁ V̇₁ = I₁                                     (points open; V₁ ≡ 0 while closed)
 *   C₂ V̇₂ = −I₂ − I_gap
 * Stored energy W = ½L_p I₁² + M I₁I₂ + ½L₂I₂² + ½C₁V₁² + ½C₂V₂²;
 *   Ẇ = e I₁ − R I₁² − R₂I₂² − V₂ I_gap.
 * The supply EMF e(t) is DC (battery) or the AC magneto (source.ts). The trapezoidal rule with the
 * trapezoidal mean ē = (e⁰ + e¹)/2 of the EMF keeps the identity W¹ − W⁰ = h·P(x̄, ē) of coil.ts for a
 * time-varying source, so the ledger (supplied h·ē·Ī₁) closes to round-off.
 *
 * Interrupters:
 *  - vibrator points (vibrator.ts): open when the armature passes `breakTravel` toward the core, close
 *    when it returns; both crossings are EVENTS located inside the sub-step (earliest event wins, the
 *    sub-step is re-taken up to it, as for the clamp / gap events of coil.ts). Closing the points
 *    discharges the condenser through them (½C₁V₁² → `energyPoints`).
 *  - timer (the caller's `timerOn`): at the make the supply is connected (I₁ = 0, so the BAT/MAG
 *    selection and its impedance can change there); at the break the loop is interrupted IDEALLY:
 *    I₁ → 0 at once with the secondary flux linkage M I₁ + L₂ I₂ conserved (no capacitance across the
 *    roller contact — the condenser bridges the points, not the timer), so I₂ jumps by M I₁/L₂ and the
 *    primary leakage energy ½(L_p − M²/L₂) I₁² is dissipated in the timer contact (`energyTimerBreak`).
 *    A break with the points closed and the current up therefore fires the coil once more ("timer
 *    spark"); with the points open (the usual case between buzz sparks) little current is broken.
 *    While the timer is open the primary is DISCONNECTED: I₁ ≡ 0 and the condenser holds its charge.
 *  - the gap (discharge.ts) as in coil.ts: breakdown at the threshold, glow/arc with V₂ algebraic,
 *    extinction, restrike — every breakdown of the train goes through the same SparkGap.
 *  - re-ignition (discharge.ts module doc): with the gap conducting, the condenser ring drives the
 *    secondary current through zero every half period of the leakage–C₁ ring (≈ 50 µs on the battery,
 *    ≈ 125 µs with the magneto's L_s); the glow goes out at I_ext and the still-ionised channel re-ignites
 *    at the recovering voltage V_r(t) once |V₂| rises to it with a sustaining current (|I₂| > I_ext and
 *    V₂·I₂ < 0, i.e. |V₂| rising — EV_REIGNITION, located like the breakdown with V_r linear over the
 *    sub-step and the conditions' own crossings). One points opening therefore gives ONE spark
 *    (breakdown) that oscillates through re-ignitions, not a breakdown per current zero.
 *
 * Sub-steps (TremblerCoilOptions): 20 µs while the points are closed (primary L/R ≈ 7–10 ms;
 * trapezoidal is A-stable, the trip is an event and the make ring of the secondary is energy-exact),
 * 1 µs during the ring-up after a points opening and whenever |V₂| may reach the (re-ignition)
 * threshold, 4 µs while the gap conducts (primary ring with the gap conducting ≈ 10 kHz: ω h ≈ 0.25;
 * ≈ 90 % of the CPU is conduction), 2 µs in the open-points ring between sparks (it sets the armature's
 * return), 20 µs for the free decay of the secondary once the timer is open: then I₁ ≡ 0 and the
 * secondary is an isolated R₂–L₂–C₂ loop whose energy cannot grow, so |V₂| never exceeds
 * √(V₂² + (L₂/C₂)I₂²) and the 20 µs step is used whenever that bound is below 0.9 × the threshold (the
 * linear reach test overestimates the ring's peak ≈ 2.7×: review finding). With the timer open, the
 * secondary at rest and the armature at rest the coil is skipped (the ringing tail is snapped to rest
 * below ½C₂(25 V)², as in coil.ts).
 */

import type { TremblerCoilSpec } from '../core/engine-spec';
import { MIN_EVENT_PROGRESS, QUIESCENT_V2, solveDense } from './coil';
import { GapTransition, type SparkGap } from './discharge';
import type { PrimarySupply } from './source';
import { Vibrator } from './vibrator';

/** Integrator controls and model constants of the trembler coil that are not in the spec. */
export interface TremblerCoilOptions {
  /** Max sub-step while the timer and the points are closed (primary ramp), s. */
  closedSubstep: number;
  /** Max sub-step during the ring-up after a points opening / near the breakdown threshold, s. */
  ringUpSubstep: number;
  /** Duration after a points opening during which the ring-up sub-step is used, s. */
  ringUpWindow: number;
  /** Max sub-step while the gap conducts, s. */
  conductingSubstep: number;
  /** Max sub-step with the points open, the timer closed and the gap open (ring between sparks), s. */
  openSubstep: number;
  /** Max sub-step with the timer open and the gap open far from breakdown (free secondary decay), s. */
  idleSubstep: number;
  /** Coefficient of restitution of the armature at its stops (UNVERIFIED default 0: inelastic). */
  vibratorRestitution: number;
}

/**
 * Default TremblerCoilOptions (accuracy vs the Radau oracle: trembler-coil.test.ts). The 20 µs ramp
 * sub-step aliases the ≈ 12 kHz secondary "make" ring of the lossless coil (L₂(1−k²) with C₂, decaying
 * only through R₂): the trip TIME stays within 1 µs of the oracle, the instantaneous firing current within
 * ≈ 1.5 % (a real core damps that ring). The open-points ring between sparks needs 2 µs for the armature's
 * return (re-close within 5 µs; 5 µs sub-steps gave 15 µs). Conducting at 4 µs vs 2 µs changes the
 * re-close by < 3 µs and the gap energy of a train by < 0.2 %.
 */
export const DEFAULT_TREMBLER_COIL_OPTIONS: Readonly<TremblerCoilOptions> = Object.freeze({
  closedSubstep: 20e-6,
  ringUpSubstep: 1e-6,
  ringUpWindow: 0.15e-3,
  conductingSubstep: 4e-6,
  openSubstep: 2e-6,
  idleSubstep: 20e-6,
  vibratorRestitution: 0,
});

const S_OPEN = 0;
const S_COND = 1;

/** Event codes of a sub-step. */
const EV_NONE = 0;
const EV_BREAKDOWN = 3;
const EV_EXTINCTION = 4;
const EV_POINTS_OPEN = 5;
const EV_POINTS_CLOSE = 6;
const EV_REIGNITION = 7;

/**
 * One trembler coil + vibrator connected to a supply through its timer contact, with a spark-gap load.
 * All fields SI. Allocation-free `step`.
 */
export class TremblerCoil {
  readonly coil: TremblerCoilSpec;
  readonly opts: TremblerCoilOptions;
  readonly supply: PrimarySupply;
  readonly vibrator: Vibrator;
  /** Timer contact (roller + wiring) resistance, Ω. */
  readonly timerResistance: number;
  /** Mutual inductance M = k√(L₁L₂), H. */
  readonly mutualInductance: number;

  // ---- circuit state ----
  /** Primary loop current I₁, A. */
  I1 = 0;
  /** Secondary current I₂, A. */
  I2 = 0;
  /** Condenser (points) voltage V₁, V. */
  V1 = 0;
  /** Secondary-capacitor (plug) voltage V₂, V. */
  V2 = 0;
  /** Internal clock, s. */
  time = 0;
  /** Timer contact closed (supply connected). */
  timerClosed = false;
  /** Vibrator points open. */
  pointsOpen = false;

  // ---- energy ledger (since resetLedger), J ----
  /** Energy delivered by the supply EMF. */
  energySupplied = 0;
  /** Dissipated in the primary winding R₁. */
  energyPrimaryResistance = 0;
  /** Dissipated in the supply's internal resistance and the timer contact (R_s + R_t). */
  energyExternalResistance = 0;
  /** Dissipated in R₂. */
  energySecondaryResistance = 0;
  /** Condenser energy dissipated in the points when they close. */
  energyPoints = 0;
  /** Primary leakage energy dissipated in the timer contact at its break. */
  energyTimerBreak = 0;
  /** Delivered to the spark gap. */
  energyGap = 0;
  /** Residual ringing energy removed when the quiescent circuit is snapped to rest (≤ ½C₂·(25 V)²). */
  energyTailTruncated = 0;
  /** Stored energy at the last resetLedger, J. */
  energyStoredRef = 0;

  // ---- vibrator events (since resetEvent) ----
  /** Points openings. */
  tripCount = 0;
  /** Time of the first / last points opening, s (NaN: none). */
  firstTripTime = NaN;
  lastTripTime = NaN;
  /** |I₁| at the first / last points opening, A. */
  firstTripCurrent = NaN;
  lastTripCurrent = NaN;
  /** Time of the last points closing, s (NaN: none). */
  lastCloseTime = NaN;
  /** Primary magnetic energy ½L₁I₁² at the last points opening, J. */
  energyAtTrip = 0;
  /** Peak |V₁| and |V₂| since the last points opening, V. */
  peakPrimaryVoltage = 0;
  peakSecondaryVoltage = 0;
  /** Diagnostic: sub-steps taken without event handling to break zero-length event loops. */
  forcedSubsteps = 0;
  /** Diagnostic: integrator sub-steps taken since construction (CPU census). */
  substepCount = 0;

  private readonly L1: number;
  private readonly R1: number;
  private readonly L2: number;
  private readonly R2: number;
  private readonly C1: number;
  private readonly C2: number;
  /** Loop self-inductance L₁ + L_s of the connected supply, H. */
  private Lp: number;
  /** External loop resistance R_s + R_t (supply + timer contact), Ω. */
  private Rx: number;
  private smode = S_OPEN;
  private sgn = 1;
  private readonly A = new Float64Array(16);
  private readonly b = new Float64Array(4);
  private readonly x0 = new Float64Array(4);
  private readonly x1 = new Float64Array(4);
  private readonly slope = new Float64Array(1);
  /** Cached supply EMF at time `eAt`. */
  private eNow = 0;
  private eAt = NaN;

  constructor(coil: TremblerCoilSpec, supply: PrimarySupply, timerResistance: number, opts?: Partial<TremblerCoilOptions>) {
    this.coil = coil;
    this.supply = supply;
    this.timerResistance = timerResistance;
    this.opts = { ...DEFAULT_TREMBLER_COIL_OPTIONS, ...opts };
    this.vibrator = new Vibrator(coil.vibrator, this.opts.vibratorRestitution);
    this.L1 = coil.primaryInductance;
    this.R1 = coil.primaryResistance;
    this.L2 = coil.secondaryInductance;
    this.R2 = coil.secondaryResistance;
    this.C1 = coil.condenserCapacitance;
    this.C2 = coil.secondaryCapacitance;
    this.mutualInductance = coil.couplingCoefficient * Math.sqrt(coil.primaryInductance * coil.secondaryInductance);
    this.Lp = this.L1 + supply.inductance;
    this.Rx = supply.resistance + timerResistance;
    this.energyStoredRef = this.storedEnergy();
  }

  /** Loop self-inductance of the connected supply, L₁ + L_s, H. */
  get loopInductance(): number {
    return this.Lp;
  }

  /** Loop resistance of the connected supply, R₁ + R_s + R_t, Ω. */
  get loopResistance(): number {
    return this.R1 + this.Rx;
  }

  /** True when the secondary is in the conducting (gap-clamped) mode. */
  get secondaryConducting(): boolean {
    return this.smode === S_COND;
  }

  /** Total stored energy W = ½L_pI₁² + M I₁I₂ + ½L₂I₂² + ½C₁V₁² + ½C₂V₂², J. */
  storedEnergy(): number {
    return (
      0.5 * this.Lp * this.I1 * this.I1 +
      this.mutualInductance * this.I1 * this.I2 +
      0.5 * this.L2 * this.I2 * this.I2 +
      0.5 * this.C1 * this.V1 * this.V1 +
      0.5 * this.C2 * this.V2 * this.V2
    );
  }

  /** Energy-balance residual W − W_ref − (supplied − dissipated − gap), J (≈ 0). */
  energyResidual(): number {
    return (
      this.storedEnergy() -
      this.energyStoredRef -
      (this.energySupplied -
        this.energyPrimaryResistance -
        this.energyExternalResistance -
        this.energySecondaryResistance -
        this.energyPoints -
        this.energyTimerBreak -
        this.energyGap -
        this.energyTailTruncated)
    );
  }

  /** Zero the energy ledger (keeps the circuit state). */
  resetLedger(): void {
    this.energySupplied = 0;
    this.energyPrimaryResistance = 0;
    this.energyExternalResistance = 0;
    this.energySecondaryResistance = 0;
    this.energyPoints = 0;
    this.energyTimerBreak = 0;
    this.energyGap = 0;
    this.energyTailTruncated = 0;
    this.energyStoredRef = this.storedEnergy();
  }

  /** Zero the per-event vibrator counters (a new timer contact). */
  resetEvent(): void {
    this.tripCount = 0;
    this.firstTripTime = NaN;
    this.lastTripTime = NaN;
    this.firstTripCurrent = NaN;
    this.lastTripCurrent = NaN;
    this.lastCloseTime = NaN;
    this.energyAtTrip = 0;
    this.peakPrimaryVoltage = 0;
    this.peakSecondaryVoltage = 0;
  }

  /** Supply EMF at the current time (cached between sub-steps), V. */
  emfNow(): number {
    if (this.eAt !== this.time) {
      this.eNow = this.supply.emf(this.time);
      this.eAt = this.time;
    }
    return this.eNow;
  }

  /** Invalidate the cached EMF (the supply's step interpolation changed). */
  invalidateEmf(): void {
    this.eAt = NaN;
  }

  /**
   * Advance the circuit by dt (s) with the timer contact `timerOn`. Timer edges take effect at the
   * start of the step (the ignition system splits caller steps at the timer angles).
   */
  step(dt: number, timerOn: boolean, gap: SparkGap): void {
    if (timerOn && !this.timerClosed) {
      // make: connect the (possibly re-selected) supply; I₁ = 0 so the loop energy is unchanged
      this.supply.applyRequested();
      this.Lp = this.L1 + this.supply.inductance;
      this.Rx = this.supply.resistance + this.timerResistance;
      this.eAt = NaN;
      this.timerClosed = true;
    } else if (!timerOn && this.timerClosed) {
      this.breakTimer(gap);
    }
    // ---- quiescent fast path: disconnected, secondary and armature at rest ----
    if (!this.timerClosed && this.smode === S_OPEN && !gap.awaitingReentry && this.settleIfQuiescent()) {
      this.time += dt;
      return;
    }
    let remaining = dt;
    let stalls = 0;
    while (remaining > 1e-15) {
      const hMax = this.maxSubstep(gap);
      const h = remaining < hMax ? remaining : hMax;
      // zero-length event guard as in coil.ts: after two stalled sub-steps the points events are
      // ignored for a sub-step (gap events stay active), after four all events
      const used = this.substep(h, gap, stalls < 2 ? 2 : stalls < 4 ? 1 : 0);
      this.substepCount++;
      if (stalls >= 2) this.forcedSubsteps++;
      if (used > MIN_EVENT_PROGRESS || used >= h) stalls = 0;
      else stalls++;
      remaining -= used;
    }
  }

  /** Ideal timer break: I₁ → 0 with the secondary flux linkage conserved. */
  private breakTimer(gap: SparkGap): void {
    const w0 = this.storedEnergy();
    this.I2 += (this.mutualInductance * this.I1) / this.L2;
    this.I1 = 0;
    this.timerClosed = false;
    this.energyTimerBreak += w0 - this.storedEnergy();
    if (this.smode === S_COND && this.I2 * this.sgn <= gap.opts.extinctionCurrent) {
      gap.extinguish(this.time);
      this.smode = S_OPEN;
    }
  }

  /**
   * True if the disconnected circuit is (or is snapped) at rest: I₁ = 0 (timer open), armature at
   * rest with the points closed, condenser empty, and the secondary ringing below ½C₂(QUIESCENT_V2)²
   * (snapped to rest, the remainder booked in energyTailTruncated — see IgnitionCoil.settleIfQuiescent).
   */
  private settleIfQuiescent(): boolean {
    const vib = this.vibrator;
    if (this.pointsOpen || vib.x !== 0 || vib.v !== 0 || this.V1 !== 0) return false;
    if (this.I2 === 0 && this.V2 === 0) return true;
    const dev = 0.5 * this.L2 * this.I2 * this.I2 + 0.5 * this.C2 * this.V2 * this.V2;
    if (dev > 0.5 * this.C2 * QUIESCENT_V2 * QUIESCENT_V2) return false;
    this.I2 = 0;
    this.V2 = 0;
    this.energyTailTruncated += dev;
    return true;
  }

  private maxSubstep(gap: SparkGap): number {
    const o = this.opts;
    if (this.smode === S_COND) return o.conductingSubstep;
    // conduction threshold: the breakdown threshold, or the (lower, rising) re-ignition voltage
    let thr = gap.threshold(Math.abs(this.I2));
    if (gap.recovering(this.time)) thr = gap.reignitionVoltage(this.time, thr);
    if (!this.timerClosed && !gap.awaitingReentry) {
      // timer open: I₁ ≡ 0, the R₂–L₂–C₂ loop's energy cannot grow (trapezoidal: exactly), so |V₂| stays
      // below this amplitude; maxSubstep runs every sub-step, so a threshold that falls later is caught
      const amp = Math.sqrt(this.V2 * this.V2 + (this.L2 / this.C2) * this.I2 * this.I2);
      if (amp < 0.9 * thr) return o.idleSubstep;
    }
    const hc = !this.timerClosed ? o.idleSubstep : this.pointsOpen ? o.openSubstep : o.closedSubstep;
    // |V₂| could reach the conduction threshold within a few candidate sub-steps (|dV₂/dt| = |I₂|/C₂)
    const reach = Math.abs(this.V2) + (4 * hc * Math.abs(this.I2)) / this.C2;
    if (reach > 0.9 * thr) return o.ringUpSubstep;
    // ring-up right after a points opening (primary ring into the condenser + secondary ring)
    if (this.pointsOpen && this.time - this.lastTripTime < o.ringUpWindow) return o.ringUpSubstep;
    return hc;
  }

  /** Trapezoidal solve of one sub-step h from x0 into x1 with EMF e0 → e1 and the current modes. */
  private solve(h: number, gap: SparkGap, e0: number, e1: number): void {
    const L2 = this.L2;
    const M = this.mutualInductance;
    const R2 = this.R2;
    const C1 = this.C1;
    const C2 = this.C2;
    const Lp = this.Lp;
    const R = this.R1 + this.Rx;
    const x0 = this.x0;
    const x1 = this.x1;
    const A = this.A;
    const b = this.b;
    const hh = 0.5 * h;
    const I10 = x0[0];
    const I20 = x0[1];
    const V10 = x0[2];
    const V20 = x0[3];
    const iters = this.smode === S_COND ? 3 : 1;
    let iStar = Math.abs(I20);
    for (let it = 0; it < iters; it++) {
      A.fill(0);
      // row 0: primary loop KVL (timer closed) or I₁ = 0 (timer open)
      if (this.timerClosed) {
        A[0] = Lp + hh * R;
        A[1] = M;
        let r = Lp * I10 + M * I20 + hh * (e0 + e1) - hh * R * I10;
        if (this.pointsOpen) {
          A[2] = hh;
          r -= hh * V10;
        }
        b[0] = r;
      } else {
        A[0] = 1;
        b[0] = 0;
      }
      // row 1: secondary KVL
      A[4] = M;
      A[5] = L2 + hh * R2;
      A[7] = -hh;
      b[1] = M * I10 + L2 * I20 + hh * V20 - hh * R2 * I20;
      // row 2: condenser (points open) or V₁ = 0 (points closed)
      if (this.pointsOpen) {
        A[8] = -hh;
        A[10] = C1;
        b[2] = C1 * V10 + hh * I10;
      } else {
        A[10] = 1;
        b[2] = 0;
      }
      // row 3: secondary node (open) or the algebraic gap voltage (conducting), as in coil.ts
      if (this.smode === S_OPEN) {
        A[13] = hh;
        A[15] = C2;
        b[3] = C2 * V20 - hh * I20;
      } else {
        const vg = gap.gapVoltageAndSlope(iStar, this.slope);
        const g = this.slope[0];
        A[13] = g;
        A[15] = 1;
        b[3] = -this.sgn * (vg - g * iStar);
      }
      solveDense(A, b, 4);
      x1[0] = this.timerClosed ? b[0] : 0;
      x1[1] = b[1];
      x1[2] = this.pointsOpen ? b[2] : 0;
      x1[3] = b[3];
      if (this.smode === S_COND) {
        const iNew = Math.abs(x1[1]);
        if (Math.abs(iNew - iStar) <= 1e-12 + 1e-9 * iNew) break;
        iStar = iNew;
      }
    }
  }

  /** One sub-step of at most h with event handling; returns the time actually advanced. */
  private substep(h: number, gap: SparkGap, events: number): number {
    const x0 = this.x0;
    const x1 = this.x1;
    const vib = this.vibrator;
    x0[0] = this.I1;
    x0[1] = this.I2;
    x0[2] = this.V1;
    x0[3] = this.V2;
    const e0 = this.emfNow();
    let e1 = this.supply.emf(this.time + h);
    this.solve(h, gap, e0, e1);
    vib.trial(h, x0[0], x1[0]);
    let hUsed = h;
    let event = EV_NONE;
    if (events > 0) {
      let thBest = 2;
      // ---- vibrator points (primary events; suppressed first by the stall guard) ----
      if (events >= 2) {
        const xb = vib.breakTravel;
        if (!this.pointsOpen && vib.xNew > xb) {
          thBest = (xb - vib.x) / (vib.xNew - vib.x);
          event = EV_POINTS_OPEN;
        } else if (this.pointsOpen && vib.xNew < xb) {
          thBest = (vib.x - xb) / (vib.x - vib.xNew);
          event = EV_POINTS_CLOSE;
        }
      }
      // ---- gap ----
      if (this.smode === S_OPEN) {
        const thr = gap.threshold(Math.abs(x1[1]));
        const v1a = Math.abs(x1[3]);
        if (v1a >= thr) {
          const v0a = Math.abs(x0[3]);
          const th = v1a > v0a ? (thr - v0a) / (v1a - v0a) : 1;
          if (th < thBest) {
            thBest = th;
            event = EV_BREAKDOWN;
          }
        }
        // re-ignition of the recovering channel (at most at the breakdown time: V_r ≤ threshold; a tie —
        // V_r = threshold at the minimum-breakdown floor of hot low-density gas — is a re-ignition)
        const th = this.reignitionFraction(h, gap, thr, v1a);
        if (th < thBest || (th === thBest && event === EV_BREAKDOWN)) {
          thBest = th;
          event = EV_REIGNITION;
        }
      } else {
        const iExt = gap.opts.extinctionCurrent;
        const i1 = x1[1] * this.sgn;
        if (i1 <= iExt) {
          const i0 = Math.abs(x0[1]);
          const th = i0 > i1 ? (i0 - iExt) / (i0 - i1) : 1;
          if (th < thBest) {
            thBest = th;
            event = EV_EXTINCTION;
          }
        }
      }
      if (event !== EV_NONE) hUsed = h * Math.min(1, Math.max(0, thBest));
      if (event !== EV_NONE && hUsed < h) {
        if (hUsed <= 0) {
          x1[0] = x0[0];
          x1[1] = x0[1];
          x1[2] = x0[2];
          x1[3] = x0[3];
          vib.xNew = vib.x;
          vib.vNew = vib.v;
          e1 = e0;
        } else {
          e1 = this.supply.emf(this.time + hUsed);
          this.solve(hUsed, gap, e0, e1);
          vib.trial(hUsed, x0[0], x1[0]);
        }
      }
    }
    this.commit(hUsed, gap, e0, e1);
    vib.commit();
    // ---- apply events ----
    if (event === EV_POINTS_OPEN) {
      vib.x = vib.breakTravel;
      this.pointsOpen = true;
      if (this.timerClosed) {
        // a trip = an interruption of the primary current (counted only with the supply connected)
        const ia = Math.abs(this.I1);
        if (this.tripCount === 0) {
          this.firstTripTime = this.time;
          this.firstTripCurrent = ia;
        }
        this.tripCount++;
        this.lastTripTime = this.time;
        this.lastTripCurrent = ia;
        this.energyAtTrip = 0.5 * this.L1 * this.I1 * this.I1;
        this.peakPrimaryVoltage = 0;
        this.peakSecondaryVoltage = 0;
      }
    } else if (event === EV_POINTS_CLOSE) {
      vib.x = vib.breakTravel;
      this.pointsOpen = false;
      this.energyPoints += 0.5 * this.C1 * this.V1 * this.V1;
      this.V1 = 0;
      this.lastCloseTime = this.time;
    } else if (event === EV_BREAKDOWN) {
      this.startConduction(gap, false);
    } else if (event === EV_REIGNITION) {
      this.startConduction(gap, true);
    } else if (event === EV_EXTINCTION) {
      gap.extinguish(this.time);
      this.smode = S_OPEN;
    } else if (this.smode === S_COND) {
      this.afterConducting(gap);
    } else if (gap.awaitingReentry && Math.abs(this.I2) <= gap.opts.extinctionCurrent) {
      gap.extinguish(this.time);
    }
    return hUsed;
  }

  /** Accept x1 as the new state and book the sub-step energies (exact trapezoidal ledger). */
  private commit(h: number, gap: SparkGap, e0: number, e1: number): void {
    const x0 = this.x0;
    const x1 = this.x1;
    const i1m = 0.5 * (x0[0] + x1[0]);
    const i2m = 0.5 * (x0[1] + x1[1]);
    const i1s = h * i1m * i1m;
    this.energySupplied += h * 0.5 * (e0 + e1) * i1m;
    this.energyPrimaryResistance += this.R1 * i1s;
    this.energyExternalResistance += this.Rx * i1s;
    this.energySecondaryResistance += h * this.R2 * i2m * i2m;
    if (this.smode === S_COND) {
      const v2m = 0.5 * (x0[3] + x1[3]);
      const e = -h * v2m * i2m - 0.5 * this.C2 * (x1[3] * x1[3] - x0[3] * x0[3]);
      this.energyGap += e;
      gap.accumulate(e, Math.abs(x1[3]), Math.abs(x1[1]), h);
    }
    this.I1 = x1[0];
    this.I2 = x1[1];
    this.V1 = x1[2];
    this.V2 = x1[3];
    this.time += h;
    this.eNow = e1;
    this.eAt = this.time;
    const v1a = Math.abs(this.V1);
    if (v1a > this.peakPrimaryVoltage) this.peakPrimaryVoltage = v1a;
    const v2a = Math.abs(this.V2);
    if (v2a > this.peakSecondaryVoltage) this.peakSecondaryVoltage = v2a;
  }

  /**
   * Fraction θ ∈ [0, 1] of the trial sub-step h (x0 → x1) at which the recovering channel re-ignites, or 2
   * if it does not: |V₂| reaches V_r(t) with a sustaining current driving it (|I₂| > I_ext and V₂·I₂ < 0)
   * at the end of the trial; θ is the latest of the linear crossings of |V₂| − V_r, |I₂| − I_ext and V₂·I₂.
   * thr / v1a: breakdown threshold (current-independent while recovering) and |V₂| at the end of the trial.
   */
  private reignitionFraction(h: number, gap: SparkGap, thr: number, v1a: number): number {
    if (!gap.recovering(this.time)) return 2;
    const x0 = this.x0;
    const x1 = this.x1;
    const iExt = gap.opts.extinctionCurrent;
    const i1a = Math.abs(x1[1]);
    if (!(i1a > iExt) || !(x1[3] * x1[1] < 0)) return 2;
    const vr1 = gap.reignitionVoltage(this.time + h, thr);
    if (v1a < vr1) return 2;
    const v0a = Math.abs(x0[3]);
    const vr0 = gap.reignitionVoltage(this.time, thr);
    let th = v0a >= vr0 ? 0 : (vr0 - v0a) / (v1a - v0a - (vr1 - vr0));
    const i0a = Math.abs(x0[1]);
    if (!(i0a > iExt)) {
      const ti = (iExt - i0a) / (i1a - i0a);
      if (ti > th) th = ti;
    }
    if (!(x0[3] * x0[1] < 0)) {
      // V₂ or I₂ changed sign in the sub-step
      if (x0[3] * x1[3] <= 0 && x0[3] !== x1[3]) {
        const tv = x0[3] / (x0[3] - x1[3]);
        if (tv > th) th = tv;
      }
      if (x0[1] * x1[1] <= 0 && x0[1] !== x1[1]) {
        const ti = x0[1] / (x0[1] - x1[1]);
        if (ti > th) th = ti;
      }
    }
    return th < 0 ? 0 : th > 1 ? 1 : th;
  }

  /**
   * Open gap reached its threshold: start (or resume) conduction (as IgnitionCoil.startConduction); with
   * `reignition` the recovering channel re-ignites (SparkGap.reignite) instead of breaking down.
   */
  private startConduction(gap: SparkGap, reignition: boolean): void {
    const c2 = this.C2;
    const wasReentry = gap.awaitingReentry;
    const vAbs = Math.abs(this.V2);
    this.sgn = this.I2 > 0 ? 1 : this.I2 < 0 ? -1 : this.V2 < 0 ? 1 : -1;
    const vAfter = reignition ? gap.reignite(vAbs, Math.abs(this.I2), c2, this.time) : gap.startConduction(vAbs, Math.abs(this.I2), c2, this.time);
    if (!gap.conducting) {
      const vNew = (this.V2 >= 0 ? 1 : -1) * vAfter;
      this.energyGap += 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
      this.V2 = vNew;
      this.smode = S_OPEN;
      return;
    }
    const vNew = -this.sgn * vAfter;
    const e = 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
    this.energyGap += e;
    if (wasReentry) gap.accumulate(e, vAfter, Math.abs(this.I2), 0);
    this.V2 = vNew;
    this.smode = S_COND;
  }

  /** Post-sub-step discharge transitions (as IgnitionCoil.afterConducting). */
  private afterConducting(gap: SparkGap): void {
    const iAbs = this.I2 * this.sgn;
    const code = gap.afterConductingStep(iAbs > 0 ? iAbs : 0, this.time);
    if (code === GapTransition.None) return;
    if (code === GapTransition.VoltageDrop) {
      const c2 = this.C2;
      const vNew = -this.sgn * gap.gapVoltage(iAbs);
      const e = 0.5 * c2 * (this.V2 * this.V2 - vNew * vNew);
      this.energyGap += e;
      gap.dump(e);
      this.V2 = vNew;
      return;
    }
    this.smode = S_OPEN;
  }
}
