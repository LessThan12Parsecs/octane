/**
 * Spark-gap discharge physics: breakdown event, arc and glow phases, extinction and
 * restrike, and the split of the electrical energy between the gas (plasma) and the
 * electrodes. SI units (V, A, W, J, m, Pa, s).
 *
 * The gap is a nonlinear load on the secondary circuit (coil.ts). Its voltage–current
 * characteristic in the sustained phases is
 *     V_gap = V_sheath + V_column,
 *     V_column = C_col · l[mm] · p[bar]^0.51 · I[A]^n            (l = channel length)
 *   glow: C_col = 40.46 V, n = −0.32 — Kim & Anderson (1995), SAE 952459 ("spark
 *         anemometry"), as used in AKTIM (Duclos & Colin 2001) and LESI; coefficients as
 *         reported by Kazmouz, Scarcelli & Bresler (ANL, JENG-2023-0023, OSTI 2530141, eq. 1:
 *         R = A p^a1 l i^a2 with A = 40.46, a1 = 0.51, a2 = −1.32 ⇒ V = R i).
 *   arc:  C_col = 6.21 V, n = −0.75 — UNVERIFIED: values attributed to the AKTIM arc-phase
 *         column in a search-engine summary of an Energies (MDPI) 2024 paper; the paper itself
 *         could not be fetched. Only enters when the current exceeds the arc–glow threshold.
 *   V_sheath (electrode falls; power dissipated there goes mostly to the electrodes):
 *   glow: normal cathode fall from Townsend theory, 3(B/A)ln(1 + 1/γ) = 295 V for air on
 *         iron (Cobine 1942 via Shaffer et al. IMECE2021-73138 eq. 1; measured 269 V by Cobine,
 *         305–348 V zero-gap-extrapolated by Shaffer et al. 2023, J. Phys. D 56:225501,
 *         Table 4; Heywood 1988 p. 429 quotes 300–500 V) + anode fall (default 0,
 *         UNVERIFIED — AKTIM uses 18.75 V; Shaffer's zero-gap value includes it).
 *   arc:  15 V total sheath drop (Shaffer et al. 2023 §4.2: "15 ± 5 V", stainless steel,
 *         1–20 atm, independent of pressure and surface finish).
 *
 * Energy to the gas (thermal plasma) — Shaffer et al. (2023) eq. 3:
 *     Ė_th = (V_gap − V_sheath) · I   (sheath power → electrodes),
 * minus radiation (arc 5 %, glow < 1 % of the electrical power: Maly & Vogel 1979 via
 * Heywood 1988 Table 9.4). With the default constants this gives ≈ 30–45 % transfer in the
 * glow at engine pressures, consistent with Maly & Vogel's 30 % ("idealised, small
 * electrodes", Heywood Table 9.4).
 *
 * Breakdown event (ns time scale, treated as instantaneous): the secondary capacitance is
 * charged to V_bd when the gap breaks down. The energy ½C_plug V_bd² stored in the plug's own
 * capacitance (directly at the gap) is the BREAKDOWN energy, transferred to the plasma with
 * 94 % efficiency (Maly & Vogel 1979 via Heywood 1988 Table 9.4: plasma 94 %, electrodes 5 %,
 * radiation < 1 %; typical 0.3–1 mJ, ~10 ns, ~40 µm channel at up to 60 000 K). The rest of
 * the lumped secondary capacitance (coil + HT lead) discharges through the freshly formed
 * channel as the CAPACITIVE ARC (µs), ½(C₂ − C_plug)V_bd² − ½C₂V_after², transferred to the
 * plasma with the arc efficiency of 50 % (Heywood Table 9.4). No suppressor resistor is
 * assumed between the lead capacitance and the gap (CFR plug lead); with a suppressor most
 * of that energy would be dissipated in the resistor instead.
 *
 * Arc ↔ glow transition: glow above ≈ 100 mA turns into an arc (Heywood 1988 p. 431, citing
 * Maly); Michler et al. (KIT) observe a stable (inductive) arc down to ≈ 100 mA. Shaffer et
 * al. (2023) stress that the threshold depends on electrode temperature/finish (cooled
 * electrodes sustain glow up to 14 A at 13 atm, Cobine) — hence a parameter.
 * Extinction: the glow cannot be sustained below `extinctionCurrent`
 * (UNVERIFIED default 2 mA; the column voltage ∝ I^−0.32 diverges as I → 0, so the circuit
 * current collapses in finite time anyway). Restrike: if the channel (stretched by the flow,
 * l = d + 2 u t — geometric U-shaped channel convected at the gap-flow velocity u, the
 * principle of spark anemometry, Kim & Anderson 1995; factor 2 is geometric, UNVERIFIED)
 * needs more voltage than a fresh breakdown of the gap, it re-breaks across the gap.
 *
 * Consistency constraints (reviewer fixes, see minimumBreakdownVoltage/startConduction):
 *  - a (re)breakdown needs at least (1 + transitionHysteresis)·V_glow(d, I_ext): below it the
 *    discharge that follows cannot be sustained. Relevant only for hot, low-density gap gas
 *    (restrikes into a kernel at low pressure), where the breakdown laws fall to the
 *    cathode-fall scale; previously this produced breakdown/extinction cascades;
 *  - the capacitive discharge after breakdown ends at V_after = min(V_mode(max(I, I_ext)),
 *    V_bd) — the gap never returns energy to the circuit;
 *  - with a circuit current ≤ I_ext the breakdown is a capacitive spark only (dump booked,
 *    gap open again).
 *
 * Re-ignition of an extinguished channel (oscillatory discharges; recovering / reignitionVoltage /
 * reignite). A channel whose current went through zero microseconds ago has not recovered its
 * dielectric strength: the re-ignition ("restrike") voltage after a current zero rises from a low
 * value back toward the breakdown voltage of the gap — Slepian's "race" between the recovery voltage
 * of the circuit and the dielectric recovery of the arc space (Slepian 1928, Trans. AIEE 47:1398);
 * short a.c. arcs in air recover in two stages, a fast one while the space-charge layer deionises and a
 * slow thermal one as the gas density recovers (Browne 1936, "Dielectric Recovery of Short A-C Arcs",
 * PhD thesis, Caltech, abstract). Here the thermal stage is the hot-kernel breakdown voltage the owner
 * sets (IgnitionSystem.updateGap: breakdown law at the kernel temperature), and the fast stage is
 *     V_r(t) = V_ri + (V_bd − V_ri)·(1 − exp(−(t − t_ext)/τ_rec)),   t − t_ext < 5 τ_rec,
 * V_ri = minimumBreakdownVoltage() = (1 + h)·V_glow(d, I_ext): the voltage of a sustainable glow,
 * which includes a new normal cathode fall (the polarity reverses) — cf. the few-hundred-volt
 * near-cathode recovery right after current zero in short a.c. arcs (Slepian; figure UNVERIFIED);
 * τ_rec = reignitionTime (default 30 µs, UNVERIFIED: no measurement for a spark-plug gap was found;
 * order of the fast stage, 10–100 µs). Oscillatory capacitor-discharge ignitions deliver energy to the
 * plug on both half cycles "on a substantially continuous basis" (US 5,513,618) — the channel is
 * re-lit at each current zero, not broken down anew. A re-ignition is a CONTINUATION of the spark:
 * counted in reignitionCount (not breakdownCount), the C₂ dump through the existing channel is booked
 * as capacitive-arc energy (capacitiveArcEfficiency / arcRadiationFraction; energyReignition ⊂
 * energyCapacitiveArc), no breakdown phase, no stepBreakdownTime (the kernel is not re-seeded). It
 * needs a circuit current above extinctionCurrent driving |V| up (the caller's event logic; below it
 * the residual channel cannot carry a sustainable discharge). After the window the gap needs a fresh
 * breakdown (threshold()). Only the trembler-coil integrator (trembler-coil.ts) uses it: its condenser
 * ring drives the secondary current through zero every ≈ 50–125 µs from tens of mA. The inductive coil
 * (coil.ts) does not: its single spark ends after a monotonic ms-long current decay (the channel has
 * deionised at extinction; the coil's ringing after the spark line does not re-fire the plug in the
 * usual secondary oscilloscope patterns — UNVERIFIED for our coil), and the CFR results stay unchanged.
 */

import { normalCathodeFall } from './breakdown';

/** Discharge mode of the gap. 'open' = not conducting. */
export type DischargeMode = 'open' | 'arc' | 'glow';

/** Tunable constants of the gap model (all SI, see module doc for sources). */
export interface SparkGapOptions {
  /** Glow cathode fall, V (default: normal cathode fall 295 V, air on iron, γ = 0.02). */
  glowCathodeFall: number;
  /** Glow anode fall, V (default 0; UNVERIFIED). */
  glowAnodeFall: number;
  /** Total arc sheath (cathode + anode) drop, V (default 15; Shaffer et al. 2023). */
  arcSheathDrop: number;
  /** Glow column coefficient C_col, V (l in mm, p in bar, I in A). Default 40.46 (Kim & Anderson 1995). */
  glowColumnCoeff: number;
  /** Glow column current exponent n (default −0.32). */
  glowCurrentExponent: number;
  /** Arc column coefficient, V (default 6.21; UNVERIFIED). */
  arcColumnCoeff: number;
  /** Arc column current exponent (default −0.75; UNVERIFIED). */
  arcCurrentExponent: number;
  /** Column pressure exponent (default 0.51; Kim & Anderson 1995). */
  pressureExponent: number;
  /** Arc → glow transition current, A (default 0.1; Heywood 1988 p. 431). */
  arcGlowCurrent: number;
  /** Relative hysteresis of the glow → arc transition (glow→arc at arcGlowCurrent·(1+h)). */
  transitionHysteresis: number;
  /** Minimum current sustaining the discharge, A (default 2e-3; UNVERIFIED). */
  extinctionCurrent: number;
  /** Plug capacitance directly at the gap (part of C₂), F (default 10 pF; UNVERIFIED typical 5–15 pF). */
  plugCapacitance: number;
  /** Fraction of the breakdown energy transferred to the plasma (0.94; Maly & Vogel 1979). */
  breakdownEfficiency: number;
  /** Fraction of the capacitive-arc energy transferred to the plasma (0.5; Maly & Vogel 1979). */
  capacitiveArcEfficiency: number;
  /** Radiated fraction of arc electrical power (0.05; Maly & Vogel 1979). */
  arcRadiationFraction: number;
  /** Radiated fraction of glow electrical power (0; "< 1 %", Maly & Vogel 1979). */
  glowRadiationFraction: number;
  /** Maximum channel length as a multiple of the gap before a forced restrike (default 20). */
  maxChannelStretch: number;
  /**
   * Recovery time constant τ_rec of the re-ignition voltage after an extinction, s (default 30 µs;
   * UNVERIFIED, see the module doc). 0 disables re-ignition (full instantaneous dielectric recovery).
   */
  reignitionTime: number;
}

/** Default gap-model constants (sources in the module doc). */
export const DEFAULT_SPARK_GAP_OPTIONS: Readonly<SparkGapOptions> = Object.freeze({
  glowCathodeFall: normalCathodeFall(),
  glowAnodeFall: 0,
  arcSheathDrop: 15,
  glowColumnCoeff: 40.46,
  glowCurrentExponent: -0.32,
  arcColumnCoeff: 6.21,
  arcCurrentExponent: -0.75,
  pressureExponent: 0.51,
  arcGlowCurrent: 0.1,
  transitionHysteresis: 0.1,
  extinctionCurrent: 2e-3,
  plugCapacitance: 10e-12,
  breakdownEfficiency: 0.94,
  capacitiveArcEfficiency: 0.5,
  arcRadiationFraction: 0.05,
  glowRadiationFraction: 0,
  maxChannelStretch: 20,
  reignitionTime: 30e-6,
});

/** Re-ignition window after an extinction, in units of reignitionTime (V_r within 0.7 % of V_bd at its end). */
export const REIGNITION_WINDOW = 5;

/**
 * Result codes of SparkGap.afterConductingStep (consumed by the coil integrator).
 * Plain numeric constants (no enum) so the module stays erasable-syntax only.
 */
export const GapTransition = {
  /** Keep conducting in the same mode. */
  None: 0,
  /** Current fell below the sustaining current: channel extinguished. */
  Extinguish: 1,
  /** Arc → glow: the gap voltage must rise; C₂ recharges from the inductor until V_glow. */
  Recharge: 2,
  /** Glow → arc: the gap voltage drops; the excess C₂ energy is dumped into the gap. */
  VoltageDrop: 3,
  /** The stretched channel needs more than a fresh breakdown: re-break across the gap. */
  Restrike: 4,
} as const;
/** One of the GapTransition codes. */
export type GapTransitionCode = (typeof GapTransition)[keyof typeof GapTransition];

/** Smallest current used to evaluate the column law (regularisation), A. */
const I_REG = 1e-4;

/**
 * Stateful spark gap. The owner (IgnitionSystem) sets `gap`, `pressure`, `breakdownVoltage`
 * and `flowVelocity` before each circuit step; the coil integrator calls the event /
 * accumulation methods. Allocation-free.
 */
export class SparkGap {
  /** Model constants. */
  readonly opts: SparkGapOptions;
  /** Electrode gap, m. */
  gap: number;
  /** Gas pressure at the gap, Pa. */
  pressure = 101325;
  /** Required breakdown voltage for a (re)strike at the current gap-gas state, V. */
  breakdownVoltage = Infinity;
  /** Mean gas velocity across the gap (convects/stretches the channel), m/s. */
  flowVelocity = 0;
  /** When true the gap never breaks down (open-circuit test). */
  suppressBreakdown = false;

  // ---- state ----
  /** Current discharge mode. */
  mode: DischargeMode = 'open';
  /** True while an arc→glow recharge is pending (gap re-conducts at the glow voltage). */
  awaitingReentry = false;
  /** Channel length, m. */
  channelLength: number;
  /** Number of breakdowns (first + restrikes + fresh breakdowns after a recovery window) since reset. */
  breakdownCount = 0;
  /** Number of re-ignitions of a recovering channel (reignite) since reset; not breakdowns. */
  reignitionCount = 0;
  /** Simulated time of the first breakdown since reset, s (NaN if none). */
  firstBreakdownTime = NaN;
  /** Simulated time at which the discharge last extinguished, s (NaN if not yet). */
  extinctionTime = NaN;
  /** Voltage required at the last breakdown, V. */
  lastBreakdownVoltage = 0;
  /** Simulated time of the last breakdown since reset, s (NaN if none). Output only. */
  lastBreakdownTime = NaN;
  /** Total time the gap conducted since reset (sum over all discharges of a spark train), s. Output only. */
  conductingTime = 0;

  // ---- energy ledger (since reset), J ----
  /** Total electrical energy delivered to the gap (all phases). */
  energyTotal = 0;
  /** Breakdown-phase energy ½C_plug V_bd². */
  energyBreakdown = 0;
  /** Capacitive-arc energy (rest of C₂ dumped at breakdown; and the re-ignition dumps). */
  energyCapacitiveArc = 0;
  /** C₂ dumps at re-ignitions, J (included in energyCapacitiveArc). */
  energyReignition = 0;
  /** Energy of the sustained (inductive) arc phase. */
  energyArc = 0;
  /** Energy of the glow phase. */
  energyGlow = 0;
  /** Energy transferred to the gas (plasma), all phases. */
  energyToGas = 0;
  /** Energy lost to the electrodes (sheaths), all phases. */
  energyToElectrodes = 0;
  /** Energy radiated, all phases. */
  energyRadiated = 0;

  // ---- per-step outputs (reset by the owner with beginStep) ----
  /** Gas energy delivered impulsively (breakdown + capacitive arc) since beginStep, J. */
  stepImpulseGasEnergy = 0;
  /** Gas energy delivered by the sustained arc/glow since beginStep, J. */
  stepContinuousGasEnergy = 0;
  /** Electrical energy into the gap since beginStep, J. */
  stepElectricalEnergy = 0;
  /** Time of the first breakdown within the current step (s, absolute), NaN if none. */
  stepBreakdownTime = NaN;
  /** Last gap voltage magnitude, V, and current magnitude, A. */
  voltage = 0;
  current = 0;

  // cache of the l·p^0.51 factor
  private cacheP = NaN;
  private cacheL = NaN;
  private cacheLP = 0;
  // cache of minimumBreakdownVoltage()
  private cacheMinP = NaN;
  private cacheMinGap = NaN;
  private cacheMinV = 0;
  /** Start of the current re-ignition recovery window (the last extinction), s; NaN: none. */
  private recoveryStart = NaN;

  constructor(gap: number, opts?: Partial<SparkGapOptions>) {
    this.opts = { ...DEFAULT_SPARK_GAP_OPTIONS, ...opts };
    this.gap = gap;
    this.channelLength = gap;
  }

  /** Reset for a new ignition event (new cycle). */
  reset(): void {
    this.mode = 'open';
    this.awaitingReentry = false;
    this.channelLength = this.gap;
    this.breakdownCount = 0;
    this.reignitionCount = 0;
    this.firstBreakdownTime = NaN;
    this.extinctionTime = NaN;
    this.recoveryStart = NaN;
    this.lastBreakdownVoltage = 0;
    this.lastBreakdownTime = NaN;
    this.conductingTime = 0;
    this.energyTotal = 0;
    this.energyBreakdown = 0;
    this.energyCapacitiveArc = 0;
    this.energyReignition = 0;
    this.energyArc = 0;
    this.energyGlow = 0;
    this.energyToGas = 0;
    this.energyToElectrodes = 0;
    this.energyRadiated = 0;
    this.voltage = 0;
    this.current = 0;
    this.beginStep();
  }

  /** Clear the per-step outputs (call before each owner step). */
  beginStep(): void {
    this.stepImpulseGasEnergy = 0;
    this.stepContinuousGasEnergy = 0;
    this.stepElectricalEnergy = 0;
    this.stepBreakdownTime = NaN;
  }

  /** True while the gap carries a discharge. */
  get conducting(): boolean {
    return this.mode !== 'open' && !this.awaitingReentry;
  }

  /** Sheath (electrode-fall) voltage of a mode, V. */
  sheathVoltage(mode: DischargeMode): number {
    const o = this.opts;
    return mode === 'arc' ? o.arcSheathDrop : o.glowCathodeFall + o.glowAnodeFall;
  }

  /** l[mm]·p[bar]^0.51 factor of the column law (cached on pressure and channel length). */
  private lengthPressureFactor(): number {
    if (this.pressure !== this.cacheP || this.channelLength !== this.cacheL) {
      this.cacheP = this.pressure;
      this.cacheL = this.channelLength;
      this.cacheLP = this.channelLength * 1e3 * Math.pow(this.pressure * 1e-5, this.opts.pressureExponent);
    }
    return this.cacheLP;
  }

  /** Column voltage of a mode at current magnitude i (A), V. */
  columnVoltage(mode: DischargeMode, i: number): number {
    const o = this.opts;
    const ii = i > I_REG ? i : I_REG;
    const lp = this.lengthPressureFactor();
    return mode === 'arc'
      ? o.arcColumnCoeff * lp * Math.pow(ii, o.arcCurrentExponent)
      : o.glowColumnCoeff * lp * Math.pow(ii, o.glowCurrentExponent);
  }

  /** Gap voltage magnitude of a mode at current magnitude i (A), V. */
  modeVoltage(mode: DischargeMode, i: number): number {
    return this.sheathVoltage(mode) + this.columnVoltage(mode, i);
  }

  /** Gap voltage magnitude in the present mode, V. */
  gapVoltage(i: number): number {
    return this.modeVoltage(this.mode === 'arc' ? 'arc' : 'glow', i);
  }

  /** dV_gap/dI in the present mode at current magnitude i, V/A (negative: falling characteristic). */
  gapVoltageSlope(i: number): number {
    const o = this.opts;
    const ii = i > I_REG ? i : I_REG;
    const n = this.mode === 'arc' ? o.arcCurrentExponent : o.glowCurrentExponent;
    return (n * this.columnVoltage(this.mode === 'arc' ? 'arc' : 'glow', ii)) / ii;
  }

  /**
   * Gap voltage magnitude (return value, V) and its slope dV/dI (written to `slopeOut[0]`,
   * V/A) in the present mode with a single power evaluation (hot path of the coil solver).
   */
  gapVoltageAndSlope(i: number, slopeOut: Float64Array): number {
    const o = this.opts;
    const arc = this.mode === 'arc';
    const ii = i > I_REG ? i : I_REG;
    const n = arc ? o.arcCurrentExponent : o.glowCurrentExponent;
    const vc = (arc ? o.arcColumnCoeff : o.glowColumnCoeff) * this.lengthPressureFactor() * Math.pow(ii, n);
    slopeOut[0] = (n * vc) / ii;
    return (arc ? o.arcSheathDrop : o.glowCathodeFall + o.glowAnodeFall) + vc;
  }

  /** Mode the discharge takes at current magnitude i when it starts conducting. */
  private modeForCurrent(i: number): DischargeMode {
    return i >= this.opts.arcGlowCurrent ? 'arc' : 'glow';
  }

  /**
   * Smallest voltage at which the gap (re)breaks down, V: (1 + transitionHysteresis) × the
   * glow voltage of the unstretched gap at the extinction current, i.e. the highest voltage
   * a sustainable discharge across the gap needs. A breakdown below it would be followed by
   * a discharge that the circuit cannot sustain (V_after > V_bd, the "breakdown" would
   * return energy). This matters only for hot, low-density gap gas (e.g. a kernel at low
   * pressure) where the breakdown laws fall to the cathode-fall scale. Reviewer fix: without
   * it the model entered breakdown/extinction cascades (up to ~10⁶ restrikes, and a
   * zero-time event loop in the coil integrator).
   */
  minimumBreakdownVoltage(): number {
    if (this.pressure !== this.cacheMinP || this.gap !== this.cacheMinGap) {
      const o = this.opts;
      const col =
        o.glowColumnCoeff *
        this.gap *
        1e3 *
        Math.pow(this.pressure * 1e-5, o.pressureExponent) *
        Math.pow(o.extinctionCurrent > I_REG ? o.extinctionCurrent : I_REG, o.glowCurrentExponent);
      this.cacheMinP = this.pressure;
      this.cacheMinGap = this.gap;
      this.cacheMinV = (1 + o.transitionHysteresis) * (o.glowCathodeFall + o.glowAnodeFall + col);
    }
    return this.cacheMinV;
  }

  /**
   * Voltage magnitude at which an open (non-conducting) gap starts to conduct, V:
   * the glow sustaining voltage during an arc→glow recharge, else the breakdown voltage
   * (not below minimumBreakdownVoltage()).
   */
  threshold(iAbs: number): number {
    if (this.awaitingReentry) return this.modeVoltage('glow', iAbs);
    if (this.suppressBreakdown) return Infinity;
    const vMin = this.minimumBreakdownVoltage();
    return this.breakdownVoltage > vMin ? this.breakdownVoltage : vMin;
  }

  /**
   * True while the extinguished channel can re-ignite at reignitionVoltage (module doc): re-ignition
   * enabled, gap open (no glow re-entry pending, breakdown not suppressed) and t within REIGNITION_WINDOW
   * time constants of the last extinction. A Restrike (stretched channel) or a reset closes the window.
   * @param t absolute time, s
   */
  recovering(t: number): boolean {
    const tau = this.opts.reignitionTime;
    if (!(tau > 0) || this.mode !== 'open' || this.awaitingReentry || this.suppressBreakdown) return false;
    const s = t - this.recoveryStart;
    return s >= 0 && s < REIGNITION_WINDOW * tau;
  }

  /**
   * Re-ignition voltage of the recovering channel at time t, V:
   * V_r = V_ri + (vFull − V_ri)(1 − e^{−(t − t_ext)/τ_rec}) with V_ri = minimumBreakdownVoltage() and
   * vFull the breakdown threshold (threshold()). Only meaningful while recovering(t).
   */
  reignitionVoltage(t: number, vFull: number): number {
    const vMin = this.minimumBreakdownVoltage();
    if (!(vFull > vMin)) return vFull;
    return vMin + (vFull - vMin) * (1 - Math.exp(-(t - this.recoveryStart) / this.opts.reignitionTime));
  }

  /**
   * The recovering channel re-ignites at |V₂| = vAbs (≤ the breakdown threshold) with circuit current
   * magnitude iAbs (> extinctionCurrent; the caller's event condition). The discharge resumes in the mode
   * of the current and C₂ dumps down to its voltage through the existing channel: booked as
   * capacitive-arc energy (also energyReignition), counted in reignitionCount — not a breakdown (no
   * breakdown phase, breakdownCount / lastBreakdown* / stepBreakdownTime untouched). Returns the
   * gap-voltage magnitude after the event (≤ vAbs). With iAbs ≤ extinctionCurrent it is a capacitive
   * pulse only and the gap is open again (a new recovery window starts).
   * @param c2 lumped secondary capacitance, F; @param t absolute time, s
   */
  reignite(vAbs: number, iAbs: number, c2: number, t: number): number {
    const o = this.opts;
    this.channelLength = this.gap;
    this.recoveryStart = NaN;
    const iEff = iAbs > o.extinctionCurrent ? iAbs : o.extinctionCurrent;
    this.mode = this.modeForCurrent(iEff);
    let vAfter = this.modeVoltage(this.mode, iEff);
    if (vAfter > vAbs) vAfter = vAbs;
    const e = Math.max(0, 0.5 * c2 * (vAbs * vAbs - vAfter * vAfter));
    const gas = o.capacitiveArcEfficiency * e;
    const rad = o.arcRadiationFraction * e;
    this.energyCapacitiveArc += e;
    this.energyReignition += e;
    this.energyTotal += e;
    this.stepElectricalEnergy += e;
    this.energyToGas += gas;
    this.stepImpulseGasEnergy += gas;
    this.energyRadiated += rad;
    this.energyToElectrodes += e - gas - rad;
    this.reignitionCount++;
    this.voltage = vAfter;
    this.current = iAbs;
    if (iAbs <= o.extinctionCurrent) this.extinguish(t);
    return vAfter;
  }

  /**
   * The gap starts conducting at secondary-capacitor voltage |V₂| = vAbs with secondary
   * current magnitude iAbs. For a breakdown, dumps the capacitor energy above the new gap
   * voltage into the gap (breakdown + capacitive arc) and books it. Returns the gap-voltage
   * magnitude after the event (the new V₂ magnitude). If the circuit current cannot sustain
   * a discharge (iAbs ≤ extinctionCurrent) the breakdown is a capacitive spark only: the
   * dump is booked and the gap is open again afterwards (`conducting` false).
   * @param c2 lumped secondary capacitance, F; @param t absolute time, s
   */
  startConduction(vAbs: number, iAbs: number, c2: number, t: number): number {
    const o = this.opts;
    if (this.awaitingReentry) {
      // Recharge complete: resume in glow at exactly the glow voltage (no dump) — unless
      // the current has decayed below the sustaining current meanwhile.
      this.awaitingReentry = false;
      if (iAbs <= o.extinctionCurrent) {
        this.extinguish(t);
        return vAbs;
      }
      this.mode = 'glow';
      return this.modeVoltage('glow', iAbs);
    }
    this.channelLength = this.gap;
    // The capacitive discharge through the fresh channel ends at the sustaining voltage of
    // the discharge the circuit current can feed (at least the extinction current), and
    // never above the breakdown voltage (the gap cannot return energy).
    const iEff = iAbs > o.extinctionCurrent ? iAbs : o.extinctionCurrent;
    this.mode = this.modeForCurrent(iEff);
    let vAfter = this.modeVoltage(this.mode, iEff);
    if (vAfter > vAbs) vAfter = vAbs;
    const eDump = Math.max(0, 0.5 * c2 * (vAbs * vAbs - vAfter * vAfter));
    const cPlug = Math.min(o.plugCapacitance, c2);
    const eBd = Math.min(0.5 * cPlug * vAbs * vAbs, eDump);
    const eCap = eDump - eBd;
    const gasBd = o.breakdownEfficiency * eBd;
    const gasCap = o.capacitiveArcEfficiency * eCap;
    const radBd = 0.01 * eBd; // "< 1 %" (Heywood Table 9.4)
    const radCap = o.arcRadiationFraction * eCap;
    this.energyBreakdown += eBd;
    this.energyCapacitiveArc += eCap;
    this.energyTotal += eDump;
    this.stepElectricalEnergy += eDump;
    this.energyToGas += gasBd + gasCap;
    this.stepImpulseGasEnergy += gasBd + gasCap;
    this.energyRadiated += radBd + radCap;
    this.energyToElectrodes += eDump - gasBd - gasCap - radBd - radCap;
    this.breakdownCount++;
    this.lastBreakdownVoltage = vAbs;
    this.lastBreakdownTime = t;
    if (Number.isNaN(this.firstBreakdownTime)) this.firstBreakdownTime = t;
    if (Number.isNaN(this.stepBreakdownTime)) this.stepBreakdownTime = t;
    this.voltage = vAfter;
    this.current = iAbs;
    // capacitive spark only: nothing sustains a discharge after the dump
    if (iAbs <= o.extinctionCurrent) this.extinguish(t);
    return vAfter;
  }

  /**
   * Book a conducting sub-step: energy into the gap e = h·V̄·Ī − ΔW_C₂ (the gap current
   * includes the secondary-capacitor current), split into gas / sheath / radiation.
   * @param e electrical energy delivered to the gap in the sub-step, J
   * @param vAbs, iAbs end-of-step gap voltage and current magnitudes
   * @param h sub-step, s (channel stretching)
   */
  accumulate(e: number, vAbs: number, iAbs: number, h: number): void {
    const o = this.opts;
    const arc = this.mode === 'arc';
    const vSheath = this.sheathVoltage(this.mode);
    const v = vAbs > vSheath ? vAbs : vSheath;
    const radFrac = arc ? o.arcRadiationFraction : o.glowRadiationFraction;
    const fSheath = vSheath / v;
    const eSheath = e * fSheath;
    const eRad = e * radFrac;
    const eGas = e - eSheath - eRad;
    this.energyTotal += e;
    this.stepElectricalEnergy += e;
    this.conductingTime += h;
    if (arc) this.energyArc += e;
    else this.energyGlow += e;
    this.energyToGas += eGas;
    this.stepContinuousGasEnergy += eGas;
    this.energyToElectrodes += eSheath;
    this.energyRadiated += eRad;
    this.voltage = vAbs;
    this.current = iAbs;
    if (this.flowVelocity > 0) {
      const lMax = o.maxChannelStretch * this.gap;
      this.channelLength = Math.min(lMax, this.channelLength + 2 * this.flowVelocity * h);
    }
  }

  /** Book an instantaneous capacitor dump into the gap (glow→arc voltage drop), J. */
  dump(e: number): void {
    if (e <= 0) return;
    const o = this.opts;
    const gas = o.capacitiveArcEfficiency * e;
    const rad = o.arcRadiationFraction * e;
    this.energyTotal += e;
    this.stepElectricalEnergy += e;
    this.energyArc += e;
    this.energyToGas += gas;
    this.stepImpulseGasEnergy += gas;
    this.energyRadiated += rad;
    this.energyToElectrodes += e - gas - rad;
  }

  /**
   * Decide what happens after a conducting sub-step ending with current magnitude iAbs.
   * Mutates the mode for VoltageDrop / Recharge / Extinguish / Restrike.
   */
  afterConductingStep(iAbs: number, t: number): GapTransitionCode {
    const o = this.opts;
    if (iAbs <= o.extinctionCurrent) {
      this.extinguish(t);
      return GapTransition.Extinguish;
    }
    // Stretched channel needs more than re-breaking the gap → restrike. The threshold is
    // ≥ (1 + h)·V_glow(gap, I_ext) > V_glow(gap, I) for I > I_ext, so the fresh channel
    // never meets the criterion again before it has been stretched.
    if (this.gapVoltage(iAbs) > this.threshold(iAbs)) {
      this.mode = 'open';
      this.channelLength = this.gap;
      this.recoveryStart = NaN; // a fresh breakdown across the gap, not a re-ignition
      return GapTransition.Restrike;
    }
    if (this.mode === 'arc' && iAbs < o.arcGlowCurrent) {
      this.mode = 'glow';
      this.awaitingReentry = true;
      return GapTransition.Recharge;
    }
    if (this.mode === 'glow' && iAbs > o.arcGlowCurrent * (1 + o.transitionHysteresis)) {
      this.mode = 'arc';
      return GapTransition.VoltageDrop;
    }
    return GapTransition.None;
  }

  /** The discharge stops (current can no longer be sustained); the re-ignition recovery window opens. */
  extinguish(t: number): void {
    this.mode = 'open';
    this.awaitingReentry = false;
    this.channelLength = this.gap;
    this.extinctionTime = t;
    this.recoveryStart = t;
    this.current = 0;
  }
}
