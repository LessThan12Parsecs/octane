/**
 * Poppet-valve flow: geometric minimum flow area, discharge coefficients, effective area,
 * shrouded-valve masking and the swirl (angular-momentum) flux of the inflow jet.
 *
 * ── Geometry (Heywood 1988, §6.3.1, Fig. 6-12) ─────────────────────────────────────────
 * D_v = headDiameter (outer seat diameter), D = D_i = seatInnerDiameter (reference for L/D and
 * the curtain area), β = seatAngle, D_s = stemDiameter, D_p = port (throat) diameter
 * (default D_i, as drawn in Fig. 6-12). The seat width is derived as the RADIAL width
 *   w = (D_v − D_i)/2,   mean seat diameter D_m = D_v − w.
 * (w is radial: with it the stage-1 formula's D_v − 2w is exactly the inner seat diameter,
 * and the stage-1/2 frusta end on the seat edges — checked against the exact minimum-area
 * frustum, see below.) Minimum flow area per valve, three stages:
 *   1. L < w/(sinβ cosβ):   A_m = π L cosβ (D_v − 2w + (L/2) sin 2β)
 *      (frustum normal to the seat from the valve's inner seat edge)
 *   2. otherwise:           A_m = π D_m √((L − w tanβ)² + w²)
 *      (frustum from the valve's inner seat edge to the head seat's outer edge)
 *   3. port limited:        A_m = (π/4)(D_p² − D_s²)
 * Stages 1–2 are continuous at L = w/(sinβ cosβ) (slope jumps from πD_v cosβ to πD_m cosβ);
 * stage 3 is taken as min(stage 1|2, port area), which is exactly Heywood's stage-3 lift
 * criterion L > √(((D_p² − D_s²)/(4D_m))² − w²) + w tanβ and stays valid when the port is so
 * small that it limits already in stage 1.
 * Formulas verified against the transcription of Heywood §6.3.1 at
 * rgmracing.free.fr/luc/heywood1 (fetched 2026-09-29; numbered there (6-7)–(6-9)).
 * UNVERIFIED: the 1988 printed equation numbers (6.2)–(6.4) quoted in the task brief.
 * Oracle: tools/reference/gasex_valve_flow.py — independent formulas + the numerically exact
 * minimum-area conical frustum between the seat faces (Heywood exceeds it by ≤ 0.26 %).
 *
 * ── Shrouded (masked) valve ─────────────────────────────────────────────────────────────
 * A shroud of arc φ_s on the valve head (CFR intake: 180°, non-rotating) blocks that part of
 * the curtain up to the shroud height h_s (default ∞: the lip covers the whole lift). The
 * curtain-type area (stages 1–2) is multiplied by the open fraction f = 1 − φ_s/360 (above
 * h_s the masked arc adds the curtain area of the gap above the lip, (1 − f)·A_{1|2}(L − h_s));
 * the port limit is unchanged: A_m = min(f A_{1|2}(L) + …, A_port).
 *
 * ── Discharge coefficient ───────────────────────────────────────────────────────────────
 * Steady-flow data are given as the CURTAIN discharge coefficient C_Dc(L/D) (reference area
 * π D L, the Annand & Roe (1974) / Kastner et al. (1963) convention) per FLOW DIRECTION:
 * inflow into the cylinder (intake forward, exhaust reverse) and outflow from the cylinder
 * (exhaust forward, intake reverse). Modelling assumption (UNVERIFIED): the coefficient is
 * governed by the direction of the flow through the seat (jet separation at the seat edges),
 * not by which port the valve sits in. Above the last tabulated L/D the flow coefficient
 * C_F = C_Dc·4L/D is held constant (port-limited saturation). Tables are interpolated with
 * monotone PCHIP. The returned C_D refers to Heywood's A_m (unshrouded):
 *   C_D(L) = C_Dc(L/D) · π D L / A_m(L),   C_D(0) = C_Dc(0)/cosβ,
 * so C_D·A_m reproduces the tabulated effective area exactly; with a shroud the effective area
 * is C_D·A_m,shrouded (= f·C_Dc·πDL below the port limit).
 *
 * ── Swirl from the shrouded inflow jet ─────────────────────────────────────────────────
 * Inflow leaves the open curtain arc radially (from the valve axis), inclined by the seat
 * angle: horizontal speed v_h = v_j cosβ (UNVERIFIED: jet direction taken along the seat cone).
 * For the open arc of half-angle α = π(1 − φ_s/360) centred on the direction ψ (shroudDirection,
 * x–z plane from +x toward +z) and valve centre (x₀, z₀) (cylinder frame), integrating the
 * moment r × v of uniformly distributed mass flux over the arc gives the angular-momentum flux
 * about the cylinder axis (+y, right-handed):
 *   Ḣ_y = ṁ v_h (sin α/α) (z₀ cos ψ − x₀ sin ψ)      [kg m²/s², i.e. N m]
 * which vanishes for an unshrouded valve (α = π). The kinetic-energy flux ½ ṁ v_j² is returned
 * for the turbulence model (mean-flow K production).
 * Magnitude caution (review): this is the IDEAL jet impulse (uniform radial curtain jets, no jet
 * spreading or wall interaction). For the CFR at 600 rpm (180° shroud, 23 mm valve offset) a
 * crude quasi-steady filling estimate gives a frictionless swirl ratio ≈ 23 at IVC (swirl ratio
 * ≈ 0.19 s/m × the mass-averaged jet speed) — plausibly several times the real swirl; the
 * turbulence model should expect to scale swirlTorqueIn by a calibrated (< 1) swirl-momentum
 * efficiency. UNVERIFIED: no CFR swirl measurement was found to check against.
 *
 * Units SI (m, m², kg/s, m/s). Hot paths allocate nothing.
 */
import type { ValveSpec } from '../core/engine-spec';
import { Pchip } from './pchip';

/** Curtain discharge coefficient table C_Dc(L/D) (reference area π D_i L). */
export interface CurtainDischargeTable {
  /** L/D abscissae (D = seatInnerDiameter), strictly increasing, starting at 0. */
  readonly lOverD: readonly number[];
  /** C_Dc at each abscissa. */
  readonly cd: readonly number[];
  /** Provenance. */
  readonly source: string;
}

/**
 * Inflow (into the cylinder) curtain discharge coefficients: production small-block
 * Chevrolet 2-valve intake port on a SuperFlow SF-1020 flow bench (suction), as reported in
 * Kamili Zahidi et al. (2017), IOP Conf. Ser.: Mater. Sci. Eng. 257 012023, §4 text on Fig. 6:
 * C_D = 0.55 (L/D = 0), 0.54 (0.05), 0.61 (0.10), 0.57 (0.15), 0.42 (0.25); the curtain
 * reference is confirmed by their flow coefficient (throat area) also being 0.42 at L/D = 0.25
 * (fetched 2026-09-29). Not CFR data: the CFR's own flow-bench coefficients (Choi et al. 2018,
 * SAE 2018-01-0848, "personal communication") are unpublished.
 * UNVERIFIED (review, re-read of the fetched paper): their L/D uses the "valve diameter" (their
 * Table 2: 0.5 mm ↔ 0.025 for the 20 mm intake valve, i.e. the HEAD diameter) and their C_D uses
 * "the valve seat area A_k"; this module references both to seatInnerDiameter (contract). For a
 * valve with D_v/D_i = 1.05 (CFR-like) the two conventions differ by ≈ 5 % in curtain area and
 * L/D; the SB-Chevy valve dimensions needed to convert exactly are not given. The L/D = 0 value
 * (0.55) is their extrapolated/lowest-lift point, not a measurement at zero lift.
 */
export const INFLOW_CURTAIN_CD: CurtainDischargeTable = Object.freeze({
  lOverD: Object.freeze([0, 0.05, 0.1, 0.15, 0.25]),
  cd: Object.freeze([0.55, 0.54, 0.61, 0.57, 0.42]),
  source: 'Kamili Zahidi et al. 2017, IOP Conf. Ser. MSE 257 012023 (SB Chevy intake, Fig. 6)',
});

/**
 * Outflow (out of the cylinder) curtain discharge coefficients.
 * UNVERIFIED: representative exhaust-valve outflow curve from memory of the Annand & Roe (1974)
 * and Blair (1999) presentations — attached flow at low lift (C_Dc ≈ 0.7) decaying as the jet
 * separates and the port starts to limit (C_F → ≈ 0.48). Replace with measured data.
 */
export const OUTFLOW_CURTAIN_CD: CurtainDischargeTable = Object.freeze({
  lOverD: Object.freeze([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]),
  cd: Object.freeze([0.7, 0.68, 0.64, 0.58, 0.52, 0.46, 0.4]),
  source: 'UNVERIFIED representative outflow curve',
});

/** Which port the valve belongs to (sets the forward flow direction). */
export type ValveKind = 'intake' | 'exhaust';

/** Options of {@link ValveFlowModel}. */
export interface ValveFlowOptions {
  /** Port (throat) diameter for the stage-3 limit, m (default: seatInnerDiameter). */
  portDiameter?: number;
  /** Shroud lip height above the valve seat, m (default ∞: masks the whole lift). */
  shroudHeight?: number;
  /** Inflow curtain C_D table (default INFLOW_CURTAIN_CD). */
  inflow?: CurtainDischargeTable;
  /** Outflow curtain C_D table (default OUTFLOW_CURTAIN_CD). */
  outflow?: CurtainDischargeTable;
}

/** Output of {@link ValveFlowModel.inflowJet} (reuse one instance). */
export interface ValveJet {
  /** Angular-momentum flux about the cylinder +y axis carried INTO the cylinder, N m. */
  angularMomentumFlux: number;
  /** Kinetic-energy flux ½ ṁ v_j² into the cylinder, W. */
  kineticEnergyFlux: number;
  /** Horizontal (in-plane) jet speed v_j cosβ, m/s. */
  horizontalVelocity: number;
}

/** Allocate a {@link ValveJet}. */
export function newValveJet(): ValveJet {
  return { angularMomentumFlux: 0, kineticEnergyFlux: 0, horizontalVelocity: 0 };
}

class CurtainCd {
  private readonly p: Pchip;
  private readonly cfSat: number;
  private readonly xMax: number;
  constructor(t: CurtainDischargeTable) {
    this.p = new Pchip(t.lOverD, t.cd);
    this.xMax = this.p.xMax;
    this.cfSat = 4 * this.xMax * this.p.yLast;
  }
  /** C_Dc at L/D = x ≥ 0 (constant C_F beyond the table). */
  at(x: number): number {
    return x <= this.xMax ? this.p.value(x) : this.cfSat / (4 * x);
  }
}

/**
 * Precomputed flow model of one valve kind (all `count` valves together). Hot-path methods
 * are allocation-free.
 */
export class ValveFlowModel {
  /** Number of identical valves. */
  readonly count: number;
  /** Valve head (outer seat) diameter D_v, m. */
  readonly headDiameter: number;
  /** Inner seat diameter D_i (L/D and curtain reference), m. */
  readonly innerDiameter: number;
  /** Radial seat width w = (D_v − D_i)/2, m. */
  readonly seatWidth: number;
  /** Mean seat diameter D_m = D_v − w, m. */
  readonly meanSeatDiameter: number;
  /** Seat angle β, rad. */
  readonly seatAngle: number;
  /** Port flow area (π/4)(D_p² − D_s²) per valve, m². */
  readonly portArea: number;
  /** Lift of the stage 1 → 2 transition w/(sinβ cosβ), m. */
  readonly stage12Lift: number;
  /** Open (unshrouded) curtain fraction f = 1 − shroudArc/360. */
  readonly openFraction: number;
  /** Shroud lip height, m (∞ = whole lift). */
  readonly shroudHeight: number;
  /** Kind of valve (sets forward direction). */
  readonly kind: ValveKind;
  private readonly cb: number;
  private readonly sb: number;
  private readonly tb: number;
  private readonly cdIn: CurtainCd;
  private readonly cdOut: CurtainCd;
  private readonly swirlArm: number;

  /**
   * @param spec valve geometry and shroud (SI; angles rad except shroudArcDeg)
   * @param kind 'intake' (forward = into the cylinder) or 'exhaust' (forward = out)
   */
  constructor(spec: ValveSpec, kind: ValveKind, opts: ValveFlowOptions = {}) {
    const Dv = spec.headDiameter;
    const Di = spec.seatInnerDiameter;
    const b = spec.seatAngle;
    if (!(Di > 0 && Dv > Di)) throw new RangeError('ValveFlowModel: need headDiameter > seatInnerDiameter > 0');
    if (!(b > 0 && b < Math.PI / 2)) throw new RangeError('ValveFlowModel: seatAngle must be in (0, π/2)');
    if (!(spec.count >= 1)) throw new RangeError('ValveFlowModel: count must be ≥ 1');
    const Dp = opts.portDiameter ?? Di;
    const Ds = spec.stemDiameter;
    if (!(Dp > Ds && Ds >= 0)) throw new RangeError('ValveFlowModel: port diameter must exceed the stem');
    this.count = spec.count;
    this.kind = kind;
    this.headDiameter = Dv;
    this.innerDiameter = Di;
    this.seatWidth = (Dv - Di) / 2;
    this.meanSeatDiameter = Dv - this.seatWidth;
    this.seatAngle = b;
    this.cb = Math.cos(b);
    this.sb = Math.sin(b);
    this.tb = Math.tan(b);
    this.portArea = (Math.PI / 4) * (Dp * Dp - Ds * Ds);
    this.stage12Lift = this.seatWidth / (this.sb * this.cb);
    const arc = Math.min(Math.max(spec.shroudArcDeg, 0), 360);
    this.openFraction = 1 - arc / 360;
    this.shroudHeight = opts.shroudHeight ?? Infinity;
    this.cdIn = new CurtainCd(opts.inflow ?? INFLOW_CURTAIN_CD);
    this.cdOut = new CurtainCd(opts.outflow ?? OUTFLOW_CURTAIN_CD);
    // swirl lever: (sin α/α)(z₀ cosψ − x₀ sinψ), α = π f
    const alpha = Math.PI * this.openFraction;
    const sinc = alpha > 0 ? Math.sin(alpha) / alpha : 1;
    const [x0, z0] = spec.position;
    const psi = spec.shroudDirection;
    this.swirlArm = arc > 0 && arc < 360 ? sinc * (z0 * Math.cos(psi) - x0 * Math.sin(psi)) : 0;
  }

  /** Stage-1/2 (seat-limited, unmasked) area of ONE valve at lift L, m². */
  seatArea(L: number): number {
    if (!(L > 0)) return 0;
    if (L < this.stage12Lift) return Math.PI * L * this.cb * (this.innerDiameter + L * this.sb * this.cb);
    const a = L - this.seatWidth * this.tb;
    return Math.PI * this.meanSeatDiameter * Math.sqrt(a * a + this.seatWidth * this.seatWidth);
  }

  /** Heywood minimum-area stage at lift L: 0 closed, 1, 2 seat-limited, 3 port-limited (with shroud). */
  stage(L: number): 0 | 1 | 2 | 3 {
    if (!(L > 0)) return 0;
    if (this.maskedSeatArea(L) >= this.portArea) return 3;
    return L < this.stage12Lift ? 1 : 2;
  }

  private maskedSeatArea(L: number): number {
    const f = this.openFraction;
    let a = f * this.seatArea(L);
    if (f < 1 && L > this.shroudHeight) a += (1 - f) * this.seatArea(L - this.shroudHeight);
    return a;
  }

  /** Geometric minimum flow area of all `count` valves at lift L (shroud included), m². */
  flowArea(L: number): number {
    if (!(L > 0)) return 0;
    const a = this.maskedSeatArea(L);
    return this.count * (a < this.portArea ? a : this.portArea);
  }

  /** Curtain area count·π D_i L, m². */
  curtainArea(L: number): number {
    return L > 0 ? this.count * Math.PI * this.innerDiameter * L : 0;
  }

  /** True when the flow passes INTO the cylinder for the given direction flag. */
  private isInflow(reverse: boolean): boolean {
    return (this.kind === 'intake') !== reverse;
  }

  /** Curtain discharge coefficient C_Dc(L/D) for the flow direction. */
  curtainDischargeCoefficient(L: number, reverse: boolean): number {
    const x = L > 0 ? L / this.innerDiameter : 0;
    return this.isInflow(reverse) ? this.cdIn.at(x) : this.cdOut.at(x);
  }

  /**
   * Discharge coefficient referred to Heywood's (unmasked) minimum area A_m:
   * C_D = C_Dc · π D_i L / A_m(L); finite as L → 0 (C_Dc(0)/cosβ).
   * @param reverse false = the valve's normal direction (intake: into, exhaust: out of the cylinder)
   */
  dischargeCoefficient(L: number, reverse: boolean): number {
    const cdc = this.curtainDischargeCoefficient(L, reverse);
    const Di = this.innerDiameter;
    if (!(L > 0)) return cdc / this.cb;
    let ratio: number; // π D_i L / A_m
    const aSeat = this.seatArea(L);
    if (aSeat >= this.portArea) ratio = (Math.PI * Di * L) / this.portArea;
    else if (L < this.stage12Lift) ratio = Di / (this.cb * (Di + L * this.sb * this.cb));
    else ratio = (Math.PI * Di * L) / aSeat;
    return cdc * ratio;
  }

  /**
   * Effective flow area C_D·A_m of all valves (shroud included), m² — the CdA to pass to the
   * orifice model. Equals count·f·C_Dc·π D_i L below the port limit.
   */
  effectiveArea(L: number, reverse: boolean): number {
    if (!(L > 0)) return 0;
    const cdc = this.curtainDischargeCoefficient(L, reverse);
    const aSeat = this.seatArea(L);
    const aUnmasked = aSeat < this.portArea ? aSeat : this.portArea;
    if (!(aUnmasked > 0)) return 0; // denormal lift
    const aMasked = this.flowArea(L) / this.count;
    return this.count * cdc * Math.PI * this.innerDiameter * L * (aMasked / aUnmasked);
  }

  /**
   * Inflow jet source terms for the in-cylinder charge-motion model. Only flow INTO the
   * cylinder carries jet momentum (outflow removes the cylinder's own angular momentum and is
   * handled by the turbulence model). Writes into and returns `out`.
   * @param mdotIn mass flow into the cylinder through these valves, kg/s (≤ 0 → zeros)
   * @param jetVelocity isentropic throat jet speed, m/s (see orifice.isentropicJetVelocity)
   */
  inflowJet(mdotIn: number, jetVelocity: number, out: ValveJet): ValveJet {
    if (!(mdotIn > 0)) {
      out.angularMomentumFlux = 0;
      out.kineticEnergyFlux = 0;
      out.horizontalVelocity = 0;
      return out;
    }
    const vh = jetVelocity * this.cb;
    out.horizontalVelocity = vh;
    out.angularMomentumFlux = mdotIn * vh * this.swirlArm;
    out.kineticEnergyFlux = 0.5 * mdotIn * jetVelocity * jetVelocity;
    return out;
  }
}

// ---- functional API (contract) ---------------------------------------------------------

const N_KEY = 10;
interface Entry {
  key: Float64Array;
  model: ValveFlowModel;
}
const models = new WeakMap<ValveSpec, { intake?: Entry; exhaust?: Entry }>();

function writeKey(spec: ValveSpec, k: Float64Array): void {
  k[0] = spec.count;
  k[1] = spec.headDiameter;
  k[2] = spec.seatInnerDiameter;
  k[3] = spec.seatAngle;
  k[4] = spec.stemDiameter;
  k[5] = spec.shroudArcDeg;
  k[6] = spec.shroudDirection;
  k[7] = spec.position[0];
  k[8] = spec.position[1];
  k[9] = 0;
}
const KEY = new Float64Array(N_KEY);

/**
 * Cached {@link ValveFlowModel} (default options) for a spec object and kind; rebuilt when
 * the spec's geometry fields change. Allocation-free on cache hits.
 */
export function valveFlowModel(spec: ValveSpec, kind: ValveKind = 'intake'): ValveFlowModel {
  writeKey(spec, KEY);
  let slot = models.get(spec);
  if (slot === undefined) {
    slot = {};
    models.set(spec, slot);
  }
  const e = slot[kind];
  if (e !== undefined) {
    let same = true;
    for (let i = 0; i < N_KEY; i++) {
      if (e.key[i] !== KEY[i]) {
        same = false;
        break;
      }
    }
    if (same) return e.model;
  }
  const model = new ValveFlowModel(spec, kind);
  slot[kind] = { key: Float64Array.from(KEY), model };
  return model;
}

/**
 * Geometric minimum flow area (Heywood stages 1–3, shroud included) of all `spec.count`
 * valves at lift L (m), m². Port diameter = seatInnerDiameter.
 */
export function valveFlowArea(spec: ValveSpec, lift: number): number {
  return valveFlowModel(spec).flowArea(lift);
}

/**
 * Discharge coefficient referred to the unmasked Heywood minimum area, for lift L (m).
 * @param reverse flow opposite to the valve's normal direction
 * @param kind which port the valve belongs to (default 'intake')
 */
export function valveDischargeCoefficient(spec: ValveSpec, lift: number, reverse: boolean, kind: ValveKind = 'intake'): number {
  return valveFlowModel(spec, kind).dischargeCoefficient(lift, reverse);
}

/** Effective area C_D·A_m (m²) of all valves of a kind at lift L (m). */
export function valveEffectiveArea(spec: ValveSpec, lift: number, reverse: boolean, kind: ValveKind = 'intake'): number {
  return valveFlowModel(spec, kind).effectiveArea(lift, reverse);
}
