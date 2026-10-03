/**
 * MOCK engine simulator — a stand-in for the real physics `EngineSimulator`
 * while it is being built, so the renderer, charts and playback pipeline have
 * complete, internally consistent, physically plausible data to consume.
 *
 * THIS IS NOT THE PHYSICS MODEL. It is a compact quasi-dimensional sketch:
 *  - exact slider-crank kinematics from the EngineSpec (CR from the operating point);
 *  - sin² valve-lift profiles honouring the spec open/close angles and the 720° wrap;
 *  - 0D gas exchange: intake plenum filled through a butterfly throttle, exhaust
 *    plenum emptied through an outlet orifice, quasi-steady compressible valve flows
 *    (bidirectional, backflow into the port is re-inducted first);
 *  - constant-R ideal gas, T-dependent γ for unburned / burned gas, Woschni heat loss;
 *  - closed cycle: two-zone (unburned isentropic − heat loss, burned elements formed at
 *    constant-pressure adiabatic temperature), Wiebe burn rate starting ~1 ms after
 *    breakdown, pressure recovered exactly from the volume constraint;
 *  - inductive-ignition electrical sketch from the IgnitionSystemSpec (dwell charging,
 *    secondary ring-up to a Paschen/Schumann breakdown voltage, arc, glow);
 *  - flame radius = sphere at the spark gap clipped by the disc chamber, inverted from
 *    the burned-zone volume;
 *  - Douaud–Eyzat / Livengood–Wu knock with damped acoustic modes (Draper);
 *  - rough equilibrium burned-gas composition and NO relaxation.
 *
 * Deterministic (seeded cycle-to-cycle variability). Many constants below are
 * approximate and marked `UNVERIFIED` or `MOCK`; none of them should be copied
 * into src/physics without verification.
 */
import { DEG } from '../physics/core/constants';
import type { EngineGeometrySpec, EngineSpec, ValveSpec } from '../physics/core/engine-spec';
import type { FuelSelection, OperatingPoint } from '../physics/core/operating-point';
import type {
  CycleSummary,
  CylinderPhase,
  EngineSnapshot,
  SparkPhase,
} from '../physics/core/snapshot';
import type { SimulatorOptions } from './protocol';
import type { SimulatorLike } from './simulator-like';

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/** MOCK: single gas constant for every gas in the cylinder, J/(kg K). Keeps pV = mRT exact. */
export const MOCK_GAS_CONSTANT = 285;
const R_GAS = MOCK_GAS_CONSTANT;
/** Dry-air gas constant, J/(kg K) (for the ambient-referenced volumetric efficiency). */
const R_AIR = 287.05;
/** Moles of N2 (+Ar lumped) per mole O2 in air: 79.05/20.95. */
const N2_PER_O2 = 3.7733;
/** Molar masses, kg/mol. N2 includes the lumped argon. */
const M_O2 = 31.998e-3;
const M_N2L = 28.15e-3;
const M_AIR = 28.96e-3;
/** Discharge coefficients (MOCK). */
const CD_INTAKE = 0.6;
const CD_EXHAUST = 0.65;
const CD_THROTTLE = 0.75;
const CD_OUTLET = 0.8;
/** Closed butterfly plate angle (from the bore normal) and leakage area fraction (MOCK). */
const THROTTLE_CLOSED_ANGLE = 5 * DEG;
const THROTTLE_LEAK_FRACTION = 0.008;
/** Orifice flow is linearised for pressure ratios above this (keeps explicit steps stable). */
const PR_LIN = 0.99;
/** Wiebe efficiency parameter (99.9 % burned at the nominal duration) and form factor. */
const WIEBE_A = 6.908;
const WIEBE_M = 2;
/** Reference laminar speed for burn-duration scaling: iso-octane B_m (Metghalchi & Keck 1982). */
const SL_REF = 0.2632;
/** Cylinder acoustic mode coefficients α_mn (Draper 1938): (1,0) and (2,0). */
const ALPHA_10 = 1.841;
const ALPHA_20 = 3.054;
/** End-gas autoignition burn time constant, s (MOCK). */
const TAU_AUTOIGNITION = 1.0e-4;
/**
 * MOCK calibration multiplier on the Douaud–Eyzat delay so that onset with ~10 % end gas
 * (≈ "standard knock intensity") lands near the ASTM D2699 guide-table trend
 * (PRF 80 ≈ CR 5.9, PRF 90 ≈ CR 6.6, PRF 100 ≈ CR 7.7 — UNVERIFIED, from memory)
 * with this mock's burn rate and end-gas temperature history.
 */
export const KNOCK_DELAY_CALIBRATION = 1.4;
/** Knock ringing: amplitude coefficient and damping time constants of the two modes (MOCK). */
const KNOCK_AMPLITUDE_COEFF = 0.2;
const KNOCK_TAU_10 = 0.8e-3;
const KNOCK_TAU_20 = 0.4e-3;
const KNOCK_MODE2_REL = 0.4;
/** Spark: arc duration and voltage, glow base voltage, breakdown-current surge impedance (MOCK). */
const ARC_DURATION = 50e-6;
const ARC_VOLTAGE = 80;
const GLOW_VOLTAGE_BASE = 400;
const BREAKDOWN_WINDOW = 1e-6;
const BREAKDOWN_SURGE_IMPEDANCE = 200;
/** Kernel radius model (MOCK): hand-off radius and the kernel → turbulent stage radius. */
const KERNEL_HANDOFF_RADIUS = 2.0e-3;
const TURBULENT_STAGE_RADIUS = 5.0e-3;
/** Burned fraction at which the flame is considered to be in burn-out / done. */
const XB_BURNOUT = 0.9;
const XB_DONE = 0.9995;

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

const clamp = (x: number, lo: number, hi: number): number => (x < lo ? lo : x > hi ? hi : x);
/** Map any angle to [0, 720). */
export const mod720 = (x: number): number => ((x % 720) + 720) % 720;

/** MOCK γ(T) of the fresh charge (fuel vapour + air + a little residual). */
export function gammaUnburned(T: number): number {
  return clamp(1.355 - 6.0e-5 * (T - 300), 1.28, 1.36);
}
/** MOCK γ(T) of burned gas (dissociation-lowered at flame temperatures). */
export function gammaBurned(T: number): number {
  return clamp(1.3 - 3.5e-5 * (T - 1000), 1.2, 1.34);
}
const gammaMix = (T: number, yProd: number): number =>
  (1 - yProd) * gammaUnburned(T) + yProd * gammaBurned(T);
const cpOf = (g: number): number => (g * R_GAS) / (g - 1);

/** Deterministic PRNG (mulberry32). */
function mulberry32(seed: number): () => number {
  let a = seed >>> 0;
  return () => {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

// ---------------------------------------------------------------------------
// Kinematics
// ---------------------------------------------------------------------------

/** Exact slider-crank kinematics (with optional wrist-pin offset). θ in crank degrees. */
export class MockKinematics {
  readonly crankRadius: number;
  readonly rodLength: number;
  readonly pinOffset: number;
  readonly bore: number;
  readonly pistonArea: number;
  readonly displacedVolume: number;
  readonly clearanceVolume: number;
  /** Head-to-crown disc height at TDC (clearance volume minus crevices over the piston area), m. */
  readonly clearanceHeightTdc: number;
  private readonly sMax: number;

  constructor(g: EngineGeometrySpec, compressionRatio: number) {
    this.crankRadius = g.stroke / 2;
    this.rodLength = g.conRodLength;
    this.pinOffset = g.pinOffset;
    this.bore = g.bore;
    this.pistonArea = (Math.PI / 4) * g.bore * g.bore;
    this.displacedVolume = this.pistonArea * g.stroke;
    this.clearanceVolume = this.displacedVolume / (compressionRatio - 1);
    this.clearanceHeightTdc = Math.max(
      1e-5,
      (this.clearanceVolume - g.creviceVolume) / this.pistonArea,
    );
    const a = this.crankRadius;
    const l = this.rodLength;
    this.sMax = Math.sqrt((l + a) * (l + a) - this.pinOffset * this.pinOffset);
  }

  /** Piston crown displacement below its highest position, m. */
  pistonDisplacement(thetaDeg: number): number {
    const th = thetaDeg * DEG;
    const a = this.crankRadius;
    const w = a * Math.sin(th) - this.pinOffset;
    return this.sMax - (a * Math.cos(th) + Math.sqrt(this.rodLength * this.rodLength - w * w));
  }
  clearanceHeight(thetaDeg: number): number {
    return this.clearanceHeightTdc + this.pistonDisplacement(thetaDeg);
  }
  volume(thetaDeg: number): number {
    return this.clearanceVolume + this.pistonArea * this.pistonDisplacement(thetaDeg);
  }
  /** Rod angle from the cylinder axis, rad: asin((r sin θ − e)/l); same sign as sin θ for e = 0. */
  rodAngle(thetaDeg: number): number {
    const th = thetaDeg * DEG;
    return Math.asin((this.crankRadius * Math.sin(th) - this.pinOffset) / this.rodLength);
  }
  /** dx/dθ, m/rad (x = piston displacement). */
  dxdTheta(thetaDeg: number): number {
    const th = thetaDeg * DEG;
    const a = this.crankRadius;
    const s = Math.sin(th);
    const c = Math.cos(th);
    const w = a * s - this.pinOffset;
    const R = Math.sqrt(this.rodLength * this.rodLength - w * w);
    return a * s + (w * a * c) / R;
  }
  /** d²x/dθ², m/rad². */
  d2xdTheta2(thetaDeg: number): number {
    const th = thetaDeg * DEG;
    const a = this.crankRadius;
    const s = Math.sin(th);
    const c = Math.cos(th);
    const w = a * s - this.pinOffset;
    const R = Math.sqrt(this.rodLength * this.rodLength - w * w);
    return a * c + (a * a * c * c - a * w * s) / R + (w * w * a * a * c * c) / (R * R * R);
  }
}

// ---------------------------------------------------------------------------
// Valves, orifices, throttle
// ---------------------------------------------------------------------------

/**
 * sin² lift profile. The spec open/close angles are quoted at `timingLiftThreshold`;
 * the profile is stretched so that lift equals the threshold exactly there.
 * Handles the 720° wrap (e.g. EVC at −345° = 15° after gas-exchange TDC).
 */
export class MockValveLift {
  /** Angle where lift leaves zero (crank deg, may be < −360). */
  readonly openExt: number;
  /** Full duration from zero lift to zero lift, crank deg. */
  readonly fullDuration: number;
  readonly maxLift: number;

  constructor(v: ValveSpec) {
    this.maxLift = v.maxLift;
    let D = mod720(v.closeDeg - v.openDeg);
    if (D <= 0) D = 360;
    const thr = clamp(v.timingLiftThreshold, 0, 0.5 * v.maxLift);
    const u0 = v.maxLift > 0 ? Math.asin(Math.sqrt(thr / v.maxLift)) / Math.PI : 0;
    this.fullDuration = Math.min(719, D / (1 - 2 * u0));
    this.openExt = v.openDeg - u0 * this.fullDuration;
  }

  lift(thetaDeg: number): number {
    const psi = mod720(thetaDeg - this.openExt);
    if (psi >= this.fullDuration) return 0;
    const s = Math.sin((Math.PI * psi) / this.fullDuration);
    return this.maxLift * s * s;
  }
}

/** Effective flow area Cd·A of a poppet valve at `lift` (curtain- or port-limited), m². */
export function valveEffectiveArea(v: ValveSpec, lift: number, cd: number): number {
  if (lift <= 0) return 0;
  const Dv = v.seatInnerDiameter;
  const curtain = Math.PI * Dv * lift;
  const port = (Math.PI / 4) * (Dv * Dv - v.stemDiameter * v.stemDiameter);
  // MOCK: a shrouded valve loses about half the masked arc of curtain.
  const shroud = 1 - 0.5 * clamp(v.shroudArcDeg, 0, 360) / 360;
  return v.count * cd * Math.min(curtain, port) * shroud;
}

/** Butterfly throttle geometric open area incl. leakage, m² (Heywood 1988 §7.3 form, no shaft). */
export function throttleArea(diameter: number, opening: number): number {
  const Ab = (Math.PI / 4) * diameter * diameter;
  const psi = THROTTLE_CLOSED_ANGLE + clamp(opening, 0, 1) * (Math.PI / 2 - THROTTLE_CLOSED_ANGLE);
  return Ab * (1 - Math.cos(psi) / Math.cos(THROTTLE_CLOSED_ANGLE)) + THROTTLE_LEAK_FRACTION * Ab;
}

function psiFlow(pr: number, g: number): number {
  const crit = Math.pow(2 / (g + 1), g / (g - 1));
  if (pr <= crit) return Math.sqrt(g) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
  return Math.sqrt(((2 * g) / (g - 1)) * (Math.pow(pr, 2 / g) - Math.pow(pr, (g + 1) / g)));
}
const PSI_LIN_14 = psiFlow(PR_LIN, 1.4);
/** Slope factor of the linearised region: dṁ/dΔp = CdA · LIN_SLOPE / sqrt(R T). */
const LIN_SLOPE = PSI_LIN_14 / (1 - PR_LIN);

/**
 * Quasi-steady isentropic nozzle mass flow, kg/s ≥ 0, from (pUp, TUp) to pDown.
 * Linearised for pr > PR_LIN so that explicit integration stays stable near equilibrium.
 */
export function orificeMassFlow(
  cdA: number,
  pUp: number,
  TUp: number,
  pDown: number,
  gamma: number,
): number {
  if (cdA <= 0 || pUp <= pDown) return 0;
  const pr = pDown / pUp;
  const psi =
    pr > PR_LIN ? (psiFlow(PR_LIN, gamma) * (1 - pr)) / (1 - PR_LIN) : psiFlow(pr, gamma);
  return (cdA * pUp * psi) / Math.sqrt(R_GAS * TUp);
}

/** Signed flow, positive A → B. */
function signedFlow(cdA: number, pA: number, TA: number, pB: number, TB: number, g: number): number {
  return pA >= pB ? orificeMassFlow(cdA, pA, TA, pB, g) : -orificeMassFlow(cdA, pB, TB, pA, g);
}

// ---------------------------------------------------------------------------
// Fuels, flame speed, ignition delay, breakdown
// ---------------------------------------------------------------------------

export interface MockFuel {
  label: string;
  /** Atoms per molecule (blend average). */
  C: number;
  H: number;
  O: number;
  /** kg/mol */
  molarMass: number;
  /** Lower heating value, J/kg. UNVERIFIED (Heywood 1988 App. D, from memory). */
  lhv: number;
  /** Octane number used by the Douaud–Eyzat correlation. */
  octane: number;
  /** Stoichiometric air/fuel mass ratio. */
  afrStoich: number;
  /** Metghalchi–Keck S_L0 = Bm + Bφ (φ − φm)², m/s. */
  bm: number;
  bphi: number;
  phim: number;
}

interface PureFuelData {
  C: number;
  H: number;
  O: number;
  lhv: number;
  octane: number;
  bm: number;
  bphi: number;
  phim: number;
}
// Iso-octane / propane MK coefficients: Metghalchi & Keck (1982). Others: MOCK guesses.
const PURE_FUELS: Record<'IC8H18' | 'NC7H16' | 'CH4' | 'C3H8' | 'C2H5OH', PureFuelData> = {
  IC8H18: { C: 8, H: 18, O: 0, lhv: 44.3e6, octane: 100, bm: 0.2632, bphi: -0.8472, phim: 1.13 },
  NC7H16: { C: 7, H: 16, O: 0, lhv: 44.6e6, octane: 0, bm: 0.285, bphi: -0.9, phim: 1.1 },
  CH4: { C: 1, H: 4, O: 0, lhv: 50.0e6, octane: 120, bm: 0.36, bphi: -1.4, phim: 1.07 },
  C3H8: { C: 3, H: 8, O: 0, lhv: 46.4e6, octane: 112, bm: 0.3422, bphi: -1.3865, phim: 1.08 },
  C2H5OH: { C: 2, H: 6, O: 1, lhv: 26.9e6, octane: 109, bm: 0.39, bphi: -1.4, phim: 1.1 },
};
const fuelMolarMass = (C: number, H: number, O: number): number =>
  (12.011 * C + 1.008 * H + 15.999 * O) * 1e-3;

/** Fuel properties for the mock. PRF octane number = liquid-volume % iso-octane. */
export function mockFuel(sel: FuelSelection): MockFuel {
  let C: number;
  let H: number;
  let O: number;
  let lhv: number;
  let octane: number;
  let bm: number;
  let bphi: number;
  let phim: number;
  let label: string;
  if (sel.kind === 'PRF') {
    const on = clamp(sel.octaneNumber, 0, 100);
    const iso = PURE_FUELS.IC8H18;
    const hep = PURE_FUELS.NC7H16;
    // Liquid densities at 20 °C, kg/m³ (UNVERIFIED): iso-octane 691.9, n-heptane 683.8.
    const mIso = on * 691.9;
    const mHep = (100 - on) * 683.8;
    const nIso = mIso / fuelMolarMass(8, 18, 0);
    const nHep = mHep / fuelMolarMass(7, 16, 0);
    const xi = nIso / (nIso + nHep);
    const yi = mIso / (mIso + mHep);
    C = 8 * xi + 7 * (1 - xi);
    H = 18 * xi + 16 * (1 - xi);
    O = 0;
    lhv = yi * iso.lhv + (1 - yi) * hep.lhv;
    octane = on;
    bm = xi * iso.bm + (1 - xi) * hep.bm;
    bphi = xi * iso.bphi + (1 - xi) * hep.bphi;
    phim = xi * iso.phim + (1 - xi) * hep.phim;
    label = `PRF ${on.toFixed(1)}`;
  } else {
    const d = PURE_FUELS[sel.species];
    ({ C, H, O, lhv, octane, bm, bphi, phim } = d);
    label = sel.species;
  }
  const molarMass = fuelMolarMass(C, H, O);
  const o2 = C + H / 4 - O / 2;
  const afrStoich = (o2 * (1 + N2_PER_O2) * M_AIR) / molarMass;
  return { label, C, H, O, molarMass, lhv, octane, afrStoich, bm, bphi, phim };
}

/** Metghalchi–Keck (1982) laminar burning velocity with residual dilution, m/s. */
export function laminarFlameSpeed(
  fuel: MockFuel,
  phi: number,
  Tu: number,
  p: number,
  dilution: number,
): number {
  const sl0 = Math.max(0, fuel.bm + fuel.bphi * (phi - fuel.phim) * (phi - fuel.phim));
  const alpha = 2.18 - 0.8 * (phi - 1);
  const beta = -0.16 + 0.22 * (phi - 1);
  const dil = Math.max(0, 1 - 2.06 * Math.pow(clamp(dilution, 0, 1), 0.77));
  return sl0 * Math.pow(Math.max(Tu, 250) / 298, alpha) * Math.pow(Math.max(p, 1e3) / 101325, beta) * dil;
}
const sl0Of = (fuel: MockFuel, phi: number): number =>
  Math.max(0, fuel.bm + fuel.bphi * (phi - fuel.phim) * (phi - fuel.phim));

/**
 * Douaud & Eyzat (1978) end-gas ignition delay for PRF, s:
 * τ = 17.68 ms (ON/100)^3.402 (p/atm)^−1.7 exp(3800/T). Octane clamped ≥ 25 (outside fit).
 */
export function douaudEyzatDelay(octane: number, p: number, T: number): number {
  const on = Math.max(25, octane);
  return 17.68e-3 * Math.pow(on / 100, 3.402) * Math.pow(p / 101325, -1.7) * Math.exp(3800 / T);
}

/**
 * Uniform-field breakdown voltage of air (Schumann-type Paschen fit), V:
 * V[kV] = 24.36 (pd) + 6.72 sqrt(pd), pd in bar·cm at 20 °C, with pressure replaced by
 * relative density. UNVERIFIED (Kuffel, Zaengl & Kuffel, High Voltage Engineering).
 */
export function breakdownVoltage(gap: number, density: number): number {
  const pd = (Math.max(density, 1e-3) / 1.2041) * 1.01325 * gap * 100;
  return 1e3 * (24.36 * pd + 6.72 * Math.sqrt(pd));
}

// ---------------------------------------------------------------------------
// Flame geometry: sphere at the gap ∩ disc chamber
// ---------------------------------------------------------------------------

/** Area of the intersection of two circles (radii r1, r2, centre distance d). */
function circleIntersectionArea(r1: number, r2: number, d: number): number {
  if (r1 <= 0 || r2 <= 0) return 0;
  if (d >= r1 + r2) return 0;
  const rmin = Math.min(r1, r2);
  if (d <= Math.abs(r1 - r2)) return Math.PI * rmin * rmin;
  const c1 = clamp((d * d + r1 * r1 - r2 * r2) / (2 * d * r1), -1, 1);
  const c2 = clamp((d * d + r2 * r2 - r1 * r1) / (2 * d * r2), -1, 1);
  const k = (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2);
  return r1 * r1 * Math.acos(c1) + r2 * r2 * Math.acos(c2) - 0.5 * Math.sqrt(Math.max(0, k));
}
/** Angle (rad) of a circle (radius ρ, centre at distance d from the bore axis) lying inside the bore. */
function arcInsideBore(rho: number, R: number, d: number): number {
  if (rho <= 0) return d <= R ? 2 * Math.PI : 0; // pole of the sphere (Archimedes: dA = 2πr dy)
  if (d + rho <= R) return 2 * Math.PI;
  if (d <= 1e-12) return rho <= R ? 2 * Math.PI : 0;
  const c = (d * d + rho * rho - R * R) / (2 * d * rho);
  if (c >= 1) return 0;
  if (c <= -1) return 2 * Math.PI;
  return 2 * Math.acos(c);
}

const SIMPSON_N = 48;

/** Volume of the sphere (radius r, centre c in cylinder frame) inside x²+z² ≤ R², −h ≤ y ≤ 0. */
export function sphereDiscVolume(
  r: number,
  c: readonly [number, number, number],
  R: number,
  h: number,
): number {
  const ylo = Math.max(-h, c[1] - r);
  const yhi = Math.min(0, c[1] + r);
  if (yhi <= ylo || r <= 0) return 0;
  const d = Math.hypot(c[0], c[2]);
  const n = SIMPSON_N;
  const dy = (yhi - ylo) / n;
  let sum = 0;
  for (let i = 0; i <= n; i++) {
    const y = ylo + i * dy;
    const q = r * r - (y - c[1]) * (y - c[1]);
    const a = q > 0 ? circleIntersectionArea(Math.sqrt(q), R, d) : 0;
    sum += (i === 0 || i === n ? 1 : i % 2 === 1 ? 4 : 2) * a;
  }
  return (sum * dy) / 3;
}

/** Area of the sphere surface lying inside the chamber (flame front area), m². Uses dA = r dφ dy. */
export function sphereDiscFrontArea(
  r: number,
  c: readonly [number, number, number],
  R: number,
  h: number,
): number {
  const ylo = Math.max(-h, c[1] - r);
  const yhi = Math.min(0, c[1] + r);
  if (yhi <= ylo || r <= 0) return 0;
  const d = Math.hypot(c[0], c[2]);
  const n = SIMPSON_N;
  const dy = (yhi - ylo) / n;
  let sum = 0;
  for (let i = 0; i <= n; i++) {
    const y = ylo + i * dy;
    const q = r * r - (y - c[1]) * (y - c[1]);
    const phi = arcInsideBore(q > 0 ? Math.sqrt(q) : 0, R, d);
    sum += (i === 0 || i === n ? 1 : i % 2 === 1 ? 4 : 2) * phi;
  }
  return (r * sum * dy) / 3;
}

/** Radius beyond which the sphere covers the whole disc chamber. */
export function sphereCoverRadius(c: readonly [number, number, number], R: number, h: number): number {
  const d = Math.hypot(c[0], c[2]);
  const dyMax = Math.max(h + c[1], -c[1]);
  return Math.hypot(R + d, dyMax);
}

/** Invert sphereDiscVolume for r (bisection). */
export function sphereRadiusForVolume(
  V: number,
  c: readonly [number, number, number],
  R: number,
  h: number,
): number {
  if (V <= 0) return 0;
  const rMax = sphereCoverRadius(c, R, h);
  if (V >= Math.PI * R * R * h * (1 - 1e-6)) return rMax;
  let lo = 0;
  let hi = rMax;
  for (let i = 0; i < 32; i++) {
    const mid = 0.5 * (lo + hi);
    if (sphereDiscVolume(mid, c, R, h) < V) lo = mid;
    else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// ---------------------------------------------------------------------------
// Composition (rough equilibrium)
// ---------------------------------------------------------------------------

export type Composition = EngineSnapshot['burnedComposition'];

interface Majors {
  CO2: number;
  H2O: number;
  O2: number;
  CO: number;
  H2: number;
  N2: number;
}

/** Complete-combustion (lean) / water-gas-shift (rich, K = 3.5) products, mol per mol fuel. */
function majorProducts(f: MockFuel, phi: number, out: Majors): Majors {
  const aSt = f.C + f.H / 4 - f.O / 2;
  const a = aSt / phi;
  out.N2 = N2_PER_O2 * a;
  if (phi <= 1) {
    out.CO2 = f.C;
    out.H2O = f.H / 2;
    out.O2 = a - aSt;
    out.CO = 0;
    out.H2 = 0;
  } else {
    const K = 3.5; // UNVERIFIED: water-gas-shift constant at ~1740 K (Heywood 1988 §4.9.1)
    const alpha = 2 * a + f.O - 2 * f.C;
    const beta = f.H / 2 - alpha;
    const bq = K * alpha + f.C + beta;
    const d = (-bq + Math.sqrt(bq * bq + 4 * (K - 1) * f.C * beta)) / (2 * (K - 1));
    const co = clamp(d, 0, f.C);
    out.CO = co;
    out.CO2 = f.C - co;
    out.H2O = Math.max(0, alpha + co);
    out.H2 = Math.max(0, f.H / 2 - out.H2O);
    out.O2 = 0;
  }
  return out;
}

// Equilibrium constants (bar-based). CO2/H2O: log-linear fits to JANAF values (UNVERIFIED,
// from memory). OH, H, O, NO: ΔG ≈ ΔH°298 − T ΔS°298 (MOCK accuracy).
const kCO2 = (T: number): number => Math.pow(10, 4.232 - 14190 / T); // CO2 ⇌ CO + ½O2
const kH2O = (T: number): number => Math.pow(10, 2.674 - 12410 / T); // H2O ⇌ H2 + ½O2
const kOH = (T: number): number => Math.exp(1.9 - 4691 / T); // ½H2 + ½O2 ⇌ OH
const kH = (T: number): number => Math.exp(5.94 - 26221 / T); // ½H2 ⇌ H
const kO = (T: number): number => Math.exp(7.04 - 29973 / T); // ½O2 ⇌ O
const kNO = (T: number): number => Math.exp(1.49 - 10861 / T); // ½N2 + ½O2 ⇌ NO

const scratchMajors: Majors = { CO2: 0, H2O: 0, O2: 0, CO: 0, H2: 0, N2: 0 };

/**
 * Rough equilibrium burned-gas mole fractions at (T, p). If `xNO` is finite it replaces
 * the equilibrium NO (kinetically controlled). Returns the equilibrium NO in `.NO` otherwise.
 */
export function burnedGasComposition(
  f: MockFuel,
  phi: number,
  T: number,
  p: number,
  xNO: number,
  out: Composition,
): Composition {
  const mj = majorProducts(f, phi, scratchMajors);
  const nTot = mj.CO2 + mj.H2O + mj.O2 + mj.CO + mj.H2 + mj.N2;
  let XCO2 = mj.CO2 / nTot;
  let XH2O = mj.H2O / nTot;
  const XO2m = mj.O2 / nTot;
  const XCOm = mj.CO / nTot;
  const XH2m = mj.H2 / nTot;
  let XN2 = mj.N2 / nTot;
  const Tt = Math.max(T, 300);
  const P = Math.max(p, 1e3) / 1e5;
  const K1 = kCO2(Tt);
  const K2 = kH2O(Tt);
  const floor = Math.pow((0.5 * (K1 * XCO2 + K2 * XH2O)) / Math.sqrt(P), 2 / 3);
  let XO2: number;
  if (phi <= 1) {
    XO2 = 0.5 * (XO2m + Math.sqrt(XO2m * XO2m + 4 * floor * floor));
  } else {
    const rich = XCOm > 0 ? Math.pow((K1 * XCO2) / XCOm, 2) / P : floor;
    XO2 = Math.min(floor, rich);
  }
  XO2 = Math.max(XO2, 1e-30);
  const sq = Math.sqrt(XO2 * P);
  const XCO = Math.max(XCOm, (K1 * XCO2) / sq);
  const XH2 = Math.max(XH2m, (K2 * XH2O) / sq);
  const XOH = kOH(Tt) * Math.sqrt(XH2 * XO2);
  const XH = kH(Tt) * Math.sqrt(XH2 / P);
  const XO = kO(Tt) * Math.sqrt(XO2 / P);
  const XNOeq = kNO(Tt) * Math.sqrt(XN2 * XO2);
  const XNO = Number.isFinite(xNO) ? xNO : XNOeq;
  XCO2 = Math.max(0, XCO2 - (XCO - XCOm));
  XH2O = Math.max(0, XH2O - (XH2 - XH2m) - 0.5 * XOH - 0.5 * XH);
  XN2 = Math.max(0, XN2 - 0.5 * XNO);
  const s = XCO2 + XH2O + XCO + XO2 + XH2 + XOH + XH + XO + XNO + XN2;
  out.CO2 = XCO2 / s;
  out.H2O = XH2O / s;
  out.CO = XCO / s;
  out.O2 = XO2 / s;
  out.H2 = XH2 / s;
  out.OH = XOH / s;
  out.H = XH / s;
  out.O = XO / s;
  out.NO = XNO / s;
  out.N2 = XN2 / s;
  return out;
}

/**
 * Unburned mixture (fresh charge at φ + a mass fraction `yProd` of cold products),
 * mole fractions of the 10 reported species (the fuel vapour is the remainder).
 */
export function unburnedMixtureComposition(
  f: MockFuel,
  phi: number,
  yProd: number,
  xNOresidual: number,
  out: Composition,
): Composition {
  const a = (f.C + f.H / 4 - f.O / 2) / phi;
  const mR = f.molarMass + a * M_O2 + N2_PER_O2 * a * M_N2L; // kg per mol fuel
  const mj = majorProducts(f, phi, scratchMajors);
  const y = clamp(yProd, 0, 1);
  const wf = (1 - y) / mR;
  const wp = y / mR;
  const nProd = mj.CO2 + mj.H2O + mj.O2 + mj.CO + mj.H2 + mj.N2;
  const nTot = wf * (1 + a + N2_PER_O2 * a) + wp * nProd;
  const xNO = (wp * nProd * Math.max(0, xNOresidual)) / nTot;
  out.CO2 = (wp * mj.CO2) / nTot;
  out.H2O = (wp * mj.H2O) / nTot;
  out.CO = (wp * mj.CO) / nTot;
  out.O2 = (wf * a + wp * mj.O2) / nTot;
  out.H2 = (wp * mj.H2) / nTot;
  out.OH = 0;
  out.H = 0;
  out.O = 0;
  out.NO = xNO;
  out.N2 = (wf * N2_PER_O2 * a + wp * mj.N2) / nTot - xNO;
  return out;
}

// ---------------------------------------------------------------------------
// Simulator
// ---------------------------------------------------------------------------

export interface MockSimulatorOptions extends SimulatorOptions {
  /** Cycle-to-cycle variability scale: 0 = none, 1 = default (~5 % burn-duration σ). */
  cyclicVariability?: number;
  /** Silent cycles run after construction/reset so the first emitted cycle is converged. */
  warmupCycles?: number;
  /** Max integration step, crank degrees (default 0.1). */
  maxStepDeg?: number;
  /** Seed of the cycle-to-cycle variability PRNG. */
  seed?: number;
  /** Multiplier on the end-gas ignition delay (default KNOCK_DELAY_CALIBRATION). */
  knockDelayCalibration?: number;
}

type FlameStage = EngineSnapshot['flame']['stage'];

interface SparkTimeline {
  tSwitch: number;
  thetaSwitch: number;
  i1: number;
  i20: number;
  vMax: number;
  vBd: number;
  omega2: number;
  broke: boolean;
  tBd: number;
  i2AtBd: number;
  tArcStart: number;
  tArcEnd: number;
  iArcEnd: number;
  vGlow: number;
  tGlowEnd: number;
  eBreakdown: number;
  eArc: number;
  eGlow: number;
}

interface CyclePlan {
  mTrap: number;
  y0: number;
  yEgr0: number;
  mFuel: number;
  qTotal: number;
  qPerMass: number;
  pIvc: number;
  TIvc: number;
  VIvc: number;
  thetaIvc: number;
  ignitable: boolean;
  delay: number;
  duration: number;
  sl0: number;
  volEff: number;
}

/** Secondary-circuit inductive discharge: L2 dI/dt = −V − R2 I. */
function dischargeCurrent(i0: number, V: number, R2: number, L2: number, tau: number): number {
  if (R2 <= 0) return Math.max(0, i0 - (V * tau) / L2);
  const t2 = L2 / R2;
  return Math.max(0, (i0 + V / R2) * Math.exp(-tau / t2) - V / R2);
}
function dischargeCharge(i0: number, V: number, R2: number, L2: number, tau: number): number {
  if (R2 <= 0) return i0 * tau - (0.5 * V * tau * tau) / L2;
  const t2 = L2 / R2;
  return (i0 + V / R2) * t2 * (1 - Math.exp(-tau / t2)) - (V / R2) * tau;
}
function dischargeEndTime(i0: number, V: number, R2: number, L2: number): number {
  if (i0 <= 0) return 0;
  if (R2 <= 0) return (L2 * i0) / V;
  return (L2 / R2) * Math.log(1 + (i0 * R2) / V);
}

export class MockSimulator implements SimulatorLike {
  private readonly spec: EngineSpec;
  private readonly snapshotEveryDeg: number;
  private readonly ccvScale: number;
  private readonly warmupCycles: number;
  private readonly maxStepDeg: number;
  private readonly seed: number;
  private readonly knockCal: number;

  /** Operating point in effect for the current cycle, and the latest requested one. */
  private op: OperatingPoint;
  private pendingOp: OperatingPoint;
  private fuel: MockFuel;
  private kin: MockKinematics;
  private readonly iv: MockValveLift;
  private readonly ev: MockValveLift;
  private rng: () => number = mulberry32(1);

  // clock
  private t = 0;
  private theta = -360;
  private cycle = 0;
  private lastTime = 0;
  private initialPending = true;

  // cylinder (single zone: m, p; two-zone adds mb, Tu, Tb)
  private m = 0;
  private p = 0;
  private yProd = 1;
  private yEgr = 0;
  private closed = false;
  private twoZone = false;
  private mb = 0;
  private Tu = 300;
  private Tb = 0;

  // plenums
  private mi = 0;
  private me = 0;
  private Te = 700;
  private pocketM = 0;
  private pocketY = 0;
  private pocketYEgr = 0;
  private pocketT = 300;

  // last-step rates (for output)
  private mdotIv = 0;
  private mdotEv = 0;
  private qLossRate = 0;
  private hrr = 0;
  private mbRate = 0;

  // per-cycle: spark
  private thetaDwell = -Infinity;
  private thetaSpark = 0;
  private dwellStarted = false;
  private tDwell = 0;
  private sparkFired = false;
  private spark: SparkTimeline | null = null;
  private bdHandled = false;
  // per-cycle: combustion
  private plan: CyclePlan | null = null;
  private ignited = false;
  private tKernel = 0;
  private thetaWiebe = NaN;
  private xb = 0;
  private burnStarted = false;
  private burnDone = false;
  private xNO = 0;
  private xNOexhaust = 0;
  // per-cycle: knock
  private lw = 0;
  private knockOn = false;
  private ringing = false;
  private knockT0 = 0;
  private knockA0 = 0;
  private knockPh1 = 0;
  private knockPh2 = 0;
  private knockFreq = 7000;
  private knockOnsetDeg = NaN;
  private knockEndGas = 0;
  // per-cycle: accumulators
  private wNet = 0;
  private wGross = 0;
  private peakP = 0;
  private peakPDeg = 0;
  private maxDpDeg = 0;
  private ca10 = NaN;
  private ca50 = NaN;
  private ca90 = NaN;
  private heatLossClosed = 0;
  private noPpm = 0;
  private coEvo = 0;
  private mapo = 0;
  private pTotPrev = 0;

  private events: number[] = [];
  private summaries: CycleSummary[] = [];

  constructor(spec: EngineSpec, op: OperatingPoint, options: MockSimulatorOptions) {
    this.spec = spec;
    this.snapshotEveryDeg = clamp(options.snapshotEveryDeg > 0 ? options.snapshotEveryDeg : 0.5, 0.01, 30);
    this.ccvScale = Math.max(0, options.cyclicVariability ?? 1);
    this.warmupCycles = Math.max(0, Math.floor(options.warmupCycles ?? 3));
    this.maxStepDeg = clamp(options.maxStepDeg ?? 0.1, 0.005, 1);
    this.seed = options.seed ?? 0x0c7a4e;
    this.knockCal = Math.max(1e-3, options.knockDelayCalibration ?? KNOCK_DELAY_CALIBRATION);
    this.pendingOp = this.sanitize(cloneOp(op));
    this.op = cloneOp(this.pendingOp);
    this.fuel = mockFuel(this.op.fuel);
    this.kin = new MockKinematics(spec.geometry, this.op.compressionRatio);
    this.iv = new MockValveLift(spec.intakeValve);
    this.ev = new MockValveLift(spec.exhaustValve);
    this.initState();
  }

  get time(): number {
    return this.lastTime;
  }

  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    const next = { ...this.pendingOp, ...patch };
    if (patch.fuel) next.fuel = { ...patch.fuel };
    this.pendingOp = this.sanitize(next);
    // Immediate: speed, throttle, load.
    this.op.rpm = this.pendingOp.rpm;
    this.op.throttle = this.pendingOp.throttle;
    this.op.loadTorque = this.pendingOp.loadTorque;
    this.op.speedMode = this.pendingOp.speedMode;
  }

  reset(): void {
    this.op = cloneOp(this.pendingOp);
    this.fuel = mockFuel(this.op.fuel);
    this.kin = new MockKinematics(this.spec.geometry, this.op.compressionRatio);
    this.initState();
  }

  advanceToNextSnapshot(): EngineSnapshot {
    if (this.initialPending) {
      this.initialPending = false;
    } else {
      this.advanceInternal();
    }
    this.lastTime = this.t;
    return this.makeSnapshot();
  }

  drainCycleSummaries(): CycleSummary[] {
    const out = this.summaries;
    this.summaries = [];
    return out;
  }

  // -------------------------------------------------------------------------
  // Setup
  // -------------------------------------------------------------------------

  private sanitize(op: OperatingPoint): OperatingPoint {
    const [crMin, crMax] = this.spec.geometry.compressionRatioRange;
    op.rpm = clamp(op.rpm, 60, 6000);
    op.throttle = clamp(op.throttle, 0, 1);
    op.equivalenceRatio = clamp(op.equivalenceRatio, 0.2, 3);
    op.compressionRatio = clamp(op.compressionRatio, Math.max(1.5, crMin), crMax);
    op.sparkAdvanceDeg = clamp(op.sparkAdvanceDeg, -60, 120);
    op.dwellTime = clamp(op.dwellTime, 1e-4, 20e-3);
    op.egrFraction = clamp(op.egrFraction, 0, 0.6);
    op.ambientPressure = clamp(op.ambientPressure, 2e4, 3e5);
    op.ambientTemperature = clamp(op.ambientTemperature, 200, 400);
    op.intakeMixtureTemperature = clamp(op.intakeMixtureTemperature, 200, 500);
    return op;
  }

  private initState(): void {
    const op = this.op;
    this.t = 0;
    this.theta = -360;
    this.cycle = 0;
    this.lastTime = 0;
    this.initialPending = true;
    this.rng = mulberry32(this.seed);
    const V = this.kin.volume(-360);
    this.p = op.ambientPressure;
    this.m = (this.p * V) / (R_GAS * 800);
    this.yProd = 1;
    this.yEgr = 0;
    this.closed = false;
    this.twoZone = false;
    this.mb = 0;
    this.Tu = 800;
    this.Tb = 0;
    const pi = this.steadyManifoldPressure();
    this.mi = (pi * this.spec.manifolds.intakeVolume) / (R_GAS * op.intakeMixtureTemperature);
    this.Te = 700;
    this.me = (op.ambientPressure * this.spec.manifolds.exhaustVolume) / (R_GAS * this.Te);
    this.pocketM = 0;
    this.pocketY = 0;
    this.pocketYEgr = 0;
    this.pocketT = op.intakeMixtureTemperature;
    this.xNOexhaust = 0;
    this.summaries = [];
    this.startCycle();
    // Warm-up: run silently until the gas-exchange state is converged.
    const target = this.warmupCycles;
    while (this.cycle < target) this.advanceInternal();
    this.t = 0;
    this.lastTime = 0;
    this.cycle = 0;
    this.summaries = [];
    this.initialPending = true;
  }

  /** Quasi-steady manifold pressure: throttle flow = speed-density engine demand (ηv ≈ 0.8). */
  private steadyManifoldPressure(): number {
    const op = this.op;
    const Tm = op.intakeMixtureTemperature;
    const cdA = CD_THROTTLE * throttleArea(this.spec.manifolds.throttleDiameter, op.throttle);
    const demand = (pm: number): number =>
      (0.8 * pm * this.kin.displacedVolume * op.rpm) / (120 * R_GAS * Tm);
    let lo = 0.02 * op.ambientPressure;
    let hi = op.ambientPressure;
    for (let i = 0; i < 60; i++) {
      const mid = 0.5 * (lo + hi);
      const f = orificeMassFlow(cdA, op.ambientPressure, Tm, mid, 1.4) - demand(mid);
      if (f > 0) lo = mid;
      else hi = mid;
    }
    return 0.5 * (lo + hi);
  }

  /** Called at θ = −360: per-cycle state reset and spark schedule. */
  private startCycle(): void {
    const op = this.op;
    const degPerSec = 6 * op.rpm;
    this.thetaSpark = -op.sparkAdvanceDeg;
    this.thetaDwell = Math.max(-360, this.thetaSpark - op.dwellTime * degPerSec);
    this.dwellStarted = false;
    this.sparkFired = false;
    this.spark = null;
    this.bdHandled = false;
    this.plan = null;
    this.ignited = false;
    this.thetaWiebe = NaN;
    this.xb = 0;
    this.burnStarted = false;
    this.burnDone = false;
    this.xNO = 0;
    this.lw = 0;
    this.knockOn = false;
    this.ringing = false;
    this.knockA0 = 0;
    this.knockOnsetDeg = NaN;
    this.knockEndGas = 0;
    this.wNet = 0;
    this.wGross = 0;
    this.peakP = 0;
    this.peakPDeg = -360;
    this.maxDpDeg = 0;
    this.ca10 = NaN;
    this.ca50 = NaN;
    this.ca90 = NaN;
    this.heatLossClosed = 0;
    this.noPpm = 0;
    this.coEvo = 0;
    this.mapo = 0;
    this.pTotPrev = this.p;
    this.events = [-180, 180];
    if (this.thetaDwell > -360) this.events.push(this.thetaDwell);
    this.events.push(this.thetaSpark);
    this.events.sort((a, b) => a - b);
  }

  // -------------------------------------------------------------------------
  // Stepping
  // -------------------------------------------------------------------------

  private nextOutputAngle(): number {
    const th = this.theta;
    const d = this.snapshotEveryDeg;
    const eps = 1e-7;
    let target = -360 + (Math.floor((th + 360) / d + 1e-9) + 1) * d;
    if (target > 360) target = 360;
    for (let i = 0; i < this.events.length; i++) {
      const e = this.events[i];
      if (e > th + eps) {
        if (e < target) target = e;
        break;
      }
    }
    if (this.ringing) {
      // Resolve the knock ringing: ≥ 12 samples per period of the (1,0) mode.
      const dRing = ((6 * this.op.rpm) / this.knockFreq) / 12;
      if (th + dRing < target) target = th + dRing;
    }
    return Math.min(target, 360);
  }

  private advanceInternal(): void {
    const target = this.nextOutputAngle();
    while (this.theta < target - 1e-12) {
      let d = Math.min(target - this.theta, this.stepLimitDeg());
      if (target - this.theta - d < 1e-6) d = target - this.theta;
      const a = this.theta;
      const b = a + d >= target - 1e-12 ? target : a + d;
      this.step(a, b);
      this.theta = b;
      this.afterStep();
    }
    // Drop passed events.
    while (this.events.length > 0 && this.events[0] <= this.theta + 1e-7) this.events.shift();
    if (this.theta >= 360 - 1e-9) this.wrapCycle();
  }

  /** Max step (deg) for accuracy/stability at the current state. */
  private stepLimitDeg(): number {
    const degPerSec = 6 * this.op.rpm;
    let lim = this.maxStepDeg;
    // Stiffness of the orifice-coupled volumes (linearised conductances).
    const k = this.stiffness();
    if (k > 0) lim = Math.min(lim, (0.25 / k) * degPerSec);
    if (this.knockOn && this.xb < XB_DONE) lim = Math.min(lim, 0.2 * TAU_AUTOIGNITION * degPerSec);
    if (this.ringing) lim = Math.min(lim, degPerSec / this.knockFreq / 16);
    return Math.max(lim, 1e-6);
  }

  private stiffness(): number {
    const op = this.op;
    const s = this.spec;
    const Tm = op.intakeMixtureTemperature;
    const gThr =
      (CD_THROTTLE * throttleArea(s.manifolds.throttleDiameter, op.throttle) * LIN_SLOPE) /
      Math.sqrt(R_GAS * Tm);
    const aOut = (Math.PI / 4) * s.manifolds.exhaustOutletDiameter ** 2;
    const gOut = (CD_OUTLET * aOut * LIN_SLOPE) / Math.sqrt(R_GAS * this.Te);
    let gIv = 0;
    let gEv = 0;
    let kc = 0;
    if (!this.closed) {
      const th = this.theta;
      const V = this.kin.volume(th);
      const T = (this.p * V) / (this.m * R_GAS);
      const lin = Math.max(this.iv.lift(th), this.iv.lift(th + this.maxStepDeg));
      const lex = Math.max(this.ev.lift(th), this.ev.lift(th + this.maxStepDeg));
      gIv = (valveEffectiveArea(s.intakeValve, lin, CD_INTAKE) * LIN_SLOPE) / Math.sqrt(R_GAS * Math.min(T, Tm));
      gEv = (valveEffectiveArea(s.exhaustValve, lex, CD_EXHAUST) * LIN_SLOPE) / Math.sqrt(R_GAS * Math.min(T, this.Te));
      kc = (gammaMix(T, this.yProd) * R_GAS * T * (gIv + gEv)) / V;
    }
    const ki = (R_GAS * Tm * (gThr + gIv)) / s.manifolds.intakeVolume;
    const ke = (R_GAS * this.Te * (gEv + gOut)) / s.manifolds.exhaustVolume;
    return Math.max(kc, ki, ke);
  }

  private wiebe(theta: number): number {
    if (!(theta > this.thetaWiebe) || !this.plan) return 0;
    const u = (theta - this.thetaWiebe) / this.plan.duration;
    return 1 - Math.exp(-WIEBE_A * Math.pow(u, WIEBE_M + 1));
  }

  /** Woschni (1967) heat-transfer coefficient, W/(m² K). */
  private woschni(p: number, T: number, theta: number): number {
    const s = this.spec.geometry;
    const Sp = (2 * s.stroke * this.op.rpm) / 60;
    let w: number;
    if (!this.closed) {
      w = 6.18 * Sp;
    } else {
      w = 2.28 * Sp;
      const pl = this.plan;
      if (pl && this.burnStarted) {
        const V = this.kin.volume(theta);
        const pMot = pl.pIvc * Math.pow(pl.VIvc / V, 1.32);
        w += Math.max(0, ((3.24e-3 * this.kin.displacedVolume * pl.TIvc) / (pl.pIvc * pl.VIvc)) * (p - pMot));
      }
    }
    return 3.26 * Math.pow(s.bore, -0.2) * Math.pow(p / 1000, 0.8) * Math.pow(Math.max(T, 250), -0.55) * Math.pow(w, 0.8);
  }

  /** Wall area (m²) and Σ A·T_wall (m² K) at θ. */
  private wallAreas(theta: number, out: { A: number; AT: number }): void {
    const w = this.spec.walls;
    const dc = this.op.coolantTemperature - 373.15;
    const Ap = this.kin.pistonArea;
    const Al = Math.PI * this.kin.bore * this.kin.clearanceHeight(theta);
    out.A = 2 * Ap + Al;
    out.AT = Ap * (w.headTemperature + 0.5 * dc) + Ap * w.pistonTemperature + Al * (w.linerTemperature + dc);
  }
  private readonly areas = { A: 0, AT: 0 };

  private step(a: number, b: number): void {
    const op = this.op;
    const s = this.spec;
    const degPerSec = 6 * op.rpm;
    const dt = (b - a) / degPerSec;
    const V0 = this.kin.volume(a);
    const V1 = this.kin.volume(b);
    const dV = V1 - V0;
    const thMid = 0.5 * (a + b);
    const Tm = op.intakeMixtureTemperature;
    const pAmb = op.ambientPressure;
    const pInt = (this.mi * R_GAS * Tm) / s.manifolds.intakeVolume;
    const pExh = (this.me * R_GAS * this.Te) / s.manifolds.exhaustVolume;
    const pPrevThermo = this.p;

    // Plenum boundary flows.
    const cdAThr = CD_THROTTLE * throttleArea(s.manifolds.throttleDiameter, op.throttle);
    const mThr = signedFlow(cdAThr, pAmb, Tm, pInt, Tm, 1.4);
    const cdAOut = CD_OUTLET * (Math.PI / 4) * s.manifolds.exhaustOutletDiameter ** 2;
    const mOut = signedFlow(cdAOut, pExh, this.Te, pAmb, op.ambientTemperature, 1.33);

    this.wallAreas(thMid, this.areas);
    let mIv = 0;
    let mEv = 0;
    let dxb = 0;

    if (this.closed) {
      // Burn increments (Wiebe on the remaining charge + end-gas autoignition).
      if (this.plan && !this.burnDone) {
        const xw0 = this.wiebe(a);
        const xw1 = this.wiebe(b);
        if (xw1 > xw0 && xw0 < 1) dxb = ((1 - this.xb) * (xw1 - xw0)) / (1 - xw0);
        if (this.knockOn) dxb += (1 - this.xb - dxb) * (1 - Math.exp(-dt / TAU_AUTOIGNITION));
        dxb = clamp(dxb, 0, 1 - this.xb);
      }
      if (dxb > 0 && !this.twoZone) {
        this.twoZone = true;
        this.burnStarted = true;
        this.mb = 0;
        this.Tu = (this.p * V0) / (this.m * R_GAS);
        this.Tb = 0;
      }
    }

    if (this.twoZone) {
      this.twoZoneStep(dt, V0, V1, dxb, thMid);
    } else {
      const T = (this.p * V0) / (this.m * R_GAS);
      const g = gammaMix(T, this.yProd);
      const cp = cpOf(g);
      let Tin = Tm;
      if (!this.closed) {
        const cdAIv = valveEffectiveArea(s.intakeValve, this.iv.lift(thMid), CD_INTAKE);
        const cdAEv = valveEffectiveArea(s.exhaustValve, this.ev.lift(thMid), CD_EXHAUST);
        Tin = this.pocketM > 0 ? this.pocketT : Tm;
        mIv = signedFlow(cdAIv, pInt, Tin, this.p, T, g);
        mEv = signedFlow(cdAEv, this.p, T, pExh, this.Te, g);
      }
      const h = this.woschni(this.p, T, thMid);
      const qLoss = h * (this.areas.A * T - this.areas.AT);
      // Enthalpy flows (cylinder cp for all streams — MOCK simplification).
      let dmIn = 0;
      let inTSum = 0;
      let dProd = 0;
      let dEgr = 0;
      if (mIv > 0) {
        dmIn = mIv * dt;
        const fromPocket = Math.min(this.pocketM, dmIn);
        const fromPlenum = dmIn - fromPocket;
        inTSum = fromPocket * this.pocketT + fromPlenum * Tm;
        dProd += fromPocket * this.pocketY + fromPlenum * op.egrFraction;
        dEgr += fromPocket * this.pocketYEgr + fromPlenum * op.egrFraction;
        this.pocketM -= fromPocket;
      } else if (mIv < 0) {
        const dmOut = -mIv * dt;
        dProd -= dmOut * this.yProd;
        dEgr -= dmOut * this.yEgr;
        const pm = this.pocketM + dmOut;
        this.pocketY = (this.pocketM * this.pocketY + dmOut * this.yProd) / pm;
        this.pocketYEgr = (this.pocketM * this.pocketYEgr + dmOut * this.yEgr) / pm;
        this.pocketT = (this.pocketM * this.pocketT + dmOut * T) / pm;
        this.pocketM = pm;
      }
      let hFlow = 0;
      if (dmIn > 0) hFlow += cp * inTSum / dt;
      else hFlow += cp * mIv * T;
      if (mEv > 0) {
        hFlow -= cp * mEv * T;
        dProd -= mEv * dt * this.yProd;
        dEgr -= mEv * dt * this.yEgr;
      } else if (mEv < 0) {
        hFlow -= cp * mEv * this.Te;
        dProd += -mEv * dt;
      }
      const dp = ((g - 1) * (hFlow - qLoss) * dt - g * this.p * dV) / V0;
      const mNew = Math.max(1e-9, this.m + (mIv - mEv) * dt);
      const prod = clamp(this.yProd * this.m + dProd, 0, mNew);
      const egr = clamp(this.yEgr * this.m + dEgr, 0, prod);
      this.m = mNew;
      this.yProd = prod / mNew;
      this.yEgr = egr / mNew;
      this.p = Math.max(1e3, this.p + dp);
      if (this.closed && this.burnDone) this.relaxNO((this.p * V1) / (this.m * R_GAS), this.p, dt);
      this.qLossRate = qLoss;
      this.hrr = 0;
      this.mbRate = 0;
      if (this.closed) this.heatLossClosed += qLoss * dt;
    }

    // Plenums.
    this.mi = Math.max(1e-9, this.mi + (mThr - mIv) * dt);
    if (this.pocketM > this.mi) this.pocketM = this.mi;
    this.pocketT += (Tm - this.pocketT) * (1 - Math.exp(-dt / 0.01));
    const meNew = Math.max(1e-9, this.me + (mEv - mOut) * dt);
    if (mEv > 0) {
      const Tcyl = (this.p * V1) / (this.m * R_GAS);
      this.Te += ((mEv * dt) / meNew) * (Tcyl - this.Te);
    }
    this.Te += (600 - this.Te) * (1 - Math.exp(-dt / 0.05)); // MOCK manifold wall cooling
    this.me = meNew;
    this.mdotIv = mIv;
    this.mdotEv = mEv;

    this.t += dt;

    // Livengood–Wu end-gas integral (closed cycle, before autoignition).
    if (this.closed && !this.knockOn && this.plan && !this.burnDone) {
      const TuNow = this.twoZone ? this.Tu : (this.p * V1) / (this.m * R_GAS);
      const lw0 = this.lw;
      this.lw += dt / (this.knockCal * douaudEyzatDelay(this.fuel.octane, this.p, TuNow));
      if (this.lw >= 1 && 1 - this.xb > 0.002) this.onKnock(a + ((1 - lw0) / (this.lw - lw0)) * (b - a));
    }

    // Knock ringing.
    let osc = 0;
    if (this.ringing) {
      const Tb = this.twoZone ? this.Tb : (this.p * V1) / (this.m * R_GAS);
      const c = Math.sqrt(gammaBurned(Tb) * R_GAS * Tb);
      this.knockFreq = (ALPHA_10 * c) / (Math.PI * s.geometry.bore);
      this.knockPh1 += 2 * Math.PI * this.knockFreq * dt;
      this.knockPh2 += 2 * Math.PI * this.knockFreq * (ALPHA_20 / ALPHA_10) * dt;
      osc = this.knockOscillation();
      this.mapo = Math.max(this.mapo, Math.abs(osc));
      if (Math.exp(-(this.t - this.knockT0) / KNOCK_TAU_10) < 0.01) this.ringing = false;
    }

    // Cycle accumulators.
    const pTot = this.p + osc;
    const w = 0.5 * (this.pTotPrev + pTot) * dV;
    this.wNet += w;
    if (a >= -180 - 1e-9 && b <= 180 + 1e-9) this.wGross += w;
    if (pTot > this.peakP) {
      this.peakP = pTot;
      this.peakPDeg = b;
    }
    if (b - a > 1e-9) this.maxDpDeg = Math.max(this.maxDpDeg, (this.p - pPrevThermo) / (b - a));
    this.pTotPrev = pTot;
  }

  private twoZoneStep(dt: number, V0: number, V1: number, dxb: number, thMid: number): void {
    const pl = this.plan!;
    const m = this.m;
    const p = this.p;
    let mb = this.mb;
    const mu = m - mb;
    const dmb = Math.min(m * dxb, mu);
    const Tu = this.Tu;
    const Tb = mb > 0 ? this.Tb : 0;
    const gu = gammaUnburned(Tu);
    const cpu = cpOf(gu);
    let Tnew = (cpu * Tu + pl.qPerMass) / cpOf(gammaBurned(2400));
    Tnew = (cpu * Tu + pl.qPerMass) / cpOf(gammaBurned(Tnew));
    const gb = gammaBurned(mb > 0 ? Tb : Tnew);
    const cpb = cpOf(gb);
    const Vu = (mu * R_GAS * Tu) / p;
    const Vb = (mb * R_GAS * Tb) / p;
    const Tmean = (mu * Tu + mb * Tb) / m;
    const h = this.woschni(p, Tmean, thMid);
    const fb = Vb / Math.max(Vu + Vb, 1e-12);
    const Qu = h * (1 - fb) * (this.areas.A * Tu - this.areas.AT);
    const Qb = mb > 0 ? h * fb * (this.areas.A * Tb - this.areas.AT) : 0;
    const denom = Vu / gu + Vb / gb;
    const dp =
      (R_GAS * (Tnew - Tu) * dmb - p * (V1 - V0) - R_GAS * ((Qu * dt) / cpu + (Qb * dt) / cpb)) / denom;
    const pPred = Math.max(0.2 * p, p + dp);
    const ru = Math.pow(pPred / p, (gu - 1) / gu);
    const rb = Math.pow(pPred / p, (gb - 1) / gb);
    const muN = mu - dmb;
    let TuN = Tu * ru - (mu > 1e-12 ? (Qu * dt) / (mu * cpu) : 0);
    TuN = Math.max(200, TuN);
    const mbTb = (mb * Tb + Tnew * dmb) * rb - (Qb * dt) / cpb;
    mb += dmb;
    const TbN = mb > 0 ? Math.max(300, mbTb / mb) : 0;
    // NO: dilution by fresh burned gas, then relaxation toward equilibrium.
    if (mb > 0) {
      this.xNO *= (mb - dmb) / mb;
      this.relaxNO(TbN, pPred, dt);
    }
    this.mb = mb;
    this.Tu = TuN;
    this.Tb = TbN;
    this.p = (R_GAS * (muN * TuN + mb * TbN)) / V1;
    this.xb = mb / m;
    this.qLossRate = Qu + Qb;
    this.heatLossClosed += (Qu + Qb) * dt;
    this.hrr = dt > 0 ? (pl.qTotal * dxb) / dt : 0;
    this.mbRate = dt > 0 ? dmb / dt : 0;
  }

  /** Heywood (1988) eq. 11.? characteristic NO time τ = 8e−16 T exp(58300/T)/√p[atm] (UNVERIFIED). */
  private relaxNO(T: number, p: number, dt: number): void {
    if (T < 1500) return;
    const tau = (8e-16 * T * Math.exp(58300 / T)) / Math.sqrt(p / 101325);
    const eq = burnedGasComposition(this.fuel, this.op.equivalenceRatio, T, p, NaN, this.scratchComp).NO;
    this.xNO += (eq - this.xNO) * (1 - Math.exp(-dt / tau));
  }
  private readonly scratchComp: Composition = {
    CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0,
  };

  private knockOscillation(): number {
    if (this.knockA0 <= 0) return 0;
    const dt = this.t - this.knockT0;
    return (
      this.knockA0 *
      (Math.exp(-dt / KNOCK_TAU_10) * Math.sin(this.knockPh1) +
        KNOCK_MODE2_REL * Math.exp(-dt / KNOCK_TAU_20) * Math.sin(this.knockPh2))
    );
  }

  private onKnock(onsetDeg: number): void {
    this.knockOn = true;
    this.ringing = true;
    this.knockT0 = this.t;
    this.knockOnsetDeg = onsetDeg;
    this.knockEndGas = 1 - this.xb;
    this.knockA0 = Math.min(0.3 * this.p, KNOCK_AMPLITUDE_COEFF * (1 - this.xb) * this.p);
    this.knockPh1 = 0;
    this.knockPh2 = 0;
  }

  /** Discrete events at the end of a step. */
  private afterStep(): void {
    const th = this.theta;
    const eps = 1e-7;
    const op = this.op;
    const degPerSec = 6 * op.rpm;

    // Dwell start.
    if (!this.dwellStarted && th >= this.thetaDwell - eps) {
      this.dwellStarted = true;
      this.tDwell = this.t - (th - Math.max(this.thetaDwell, -360)) / degPerSec;
    }
    // Switch-open (spark command).
    if (!this.sparkFired && th >= this.thetaSpark - eps) {
      this.sparkFired = true;
      this.fireSpark();
    }
    // Breakdown.
    const sp = this.spark;
    if (sp && sp.broke && !this.bdHandled && this.t >= sp.tBd - 1e-12) {
      this.bdHandled = true;
      this.onBreakdown();
    }
    // Valve events.
    const closedNow = this.iv.lift(th) === 0 && this.ev.lift(th) === 0 && Math.abs(th) < 270;
    if (closedNow && !this.closed) this.onIvc();
    else if (!closedNow && this.closed) this.onEvo();
    // Burn completion.
    if (this.twoZone) {
      const wiebeEnd = this.plan && !Number.isNaN(this.thetaWiebe) ? this.thetaWiebe + this.plan.duration : Infinity;
      if (this.xb >= XB_DONE || (!this.knockOn && th >= wiebeEnd)) this.mergeZones();
    }
    // CA10/50/90 are tracked as xb crosses the levels (xb only changes in two-zone steps).
    if (this.burnStarted) {
      if (this.ca10 !== this.ca10 && this.xb >= 0.1) this.ca10 = th;
      if (this.ca50 !== this.ca50 && this.xb >= 0.5) this.ca50 = th;
      if (this.ca90 !== this.ca90 && this.xb >= 0.9) this.ca90 = th;
    }
  }

  private mergeZones(): void {
    // p is continuous; the single zone takes the mass-averaged temperature (same R).
    this.twoZone = false;
    this.burnDone = true;
    this.mb = this.m;
    this.yProd = 1;
    this.yEgr = 0;
  }

  private fireSpark(): void {
    const ig = this.spec.ignition;
    const op = this.op;
    const L1 = ig.primaryInductance;
    const R1 = Math.max(ig.primaryResistance, 1e-6);
    const L2 = ig.secondaryInductance;
    const C2 = ig.secondaryCapacitance;
    const R2 = ig.secondaryResistance;
    const k = ig.couplingCoefficient;
    const dwell = this.dwellStarted ? this.t - this.tDwell : 0;
    const i1 = Math.min(
      ig.primaryCurrentLimit,
      (ig.supplyVoltage / R1) * (1 - Math.exp((-dwell * R1) / L1)),
    );
    const V = this.kin.volume(this.theta);
    const rho = this.m / V;
    const vBd = breakdownVoltage(this.spec.sparkPlug.gap, rho);
    const vMax = k * i1 * Math.sqrt(L1 / C2);
    const i20 = k * i1 * Math.sqrt(L1 / L2);
    const omega2 = 1 / Math.sqrt(L2 * C2);
    const broke = vMax > vBd;
    const tS = this.t;
    const sp: SparkTimeline = {
      tSwitch: tS,
      thetaSwitch: this.theta,
      i1,
      i20,
      vMax,
      vBd,
      omega2,
      broke,
      tBd: Infinity,
      i2AtBd: 0,
      tArcStart: Infinity,
      tArcEnd: Infinity,
      iArcEnd: 0,
      vGlow: GLOW_VOLTAGE_BASE,
      tGlowEnd: Infinity,
      eBreakdown: 0,
      eArc: 0,
      eGlow: 0,
    };
    const degPerSec = 6 * op.rpm;
    const ev: number[] = [];
    if (broke) {
      const tRing = Math.asin(vBd / vMax) / omega2;
      sp.tBd = tS + tRing;
      sp.i2AtBd = i20 * Math.cos(omega2 * tRing);
      sp.eBreakdown = 0.5 * C2 * vBd * vBd;
      sp.tArcStart = sp.tBd + BREAKDOWN_WINDOW;
      const arcDur = Math.min(ARC_DURATION, dischargeEndTime(sp.i2AtBd, ARC_VOLTAGE, R2, L2));
      sp.tArcEnd = sp.tArcStart + arcDur;
      sp.iArcEnd = dischargeCurrent(sp.i2AtBd, ARC_VOLTAGE, R2, L2, arcDur);
      sp.eArc = ARC_VOLTAGE * dischargeCharge(sp.i2AtBd, ARC_VOLTAGE, R2, L2, arcDur);
      // MOCK: glow voltage rises with gas density.
      sp.vGlow = GLOW_VOLTAGE_BASE + 30 * (rho / 1.2041 - 1);
      const glowDur = dischargeEndTime(sp.iArcEnd, sp.vGlow, R2, L2);
      sp.tGlowEnd = sp.tArcEnd + glowDur;
      sp.eGlow = sp.vGlow * dischargeCharge(sp.iArcEnd, sp.vGlow, R2, L2, glowDur);
      const times = [
        0.5 * tRing,
        tRing,
        tRing + BREAKDOWN_WINDOW,
        tRing + 10e-6,
        tRing + 25e-6,
        sp.tArcEnd - tS,
        sp.tArcEnd - tS + 0.1e-3,
        sp.tArcEnd - tS + 0.3e-3,
        sp.tArcEnd - tS + 0.6e-3,
        sp.tGlowEnd - tS,
      ];
      for (const dt of times) ev.push(this.theta + dt * degPerSec);
    } else {
      const tHalf = Math.PI / omega2;
      ev.push(this.theta + 0.25 * tHalf * degPerSec, this.theta + tHalf * degPerSec);
    }
    this.spark = sp;
    this.addEvents(ev);
  }

  private addEvents(ev: number[]): void {
    for (const e of ev) if (e > this.theta && e < 360) this.events.push(e);
    this.events.sort((a, b) => a - b);
  }

  private onBreakdown(): void {
    const pl = this.plan;
    if (!this.closed || !pl || !pl.ignitable || this.burnDone) return;
    this.ignited = true;
    this.tKernel = this.spark!.tBd;
    this.thetaWiebe = this.theta + pl.delay * 6 * this.op.rpm;
    this.addEvents([this.thetaWiebe]);
  }

  private onIvc(): void {
    this.closed = true;
    const op = this.op;
    const f = this.fuel;
    const V = this.kin.volume(this.theta);
    const T = (this.p * V) / (this.m * R_GAS);
    const phi = op.equivalenceRatio;
    const far = phi / f.afrStoich;
    const y0 = this.yProd;
    const mFuel = (this.m * (1 - y0) * far) / (1 + far);
    const etaC = phi <= 1 ? 0.98 : 0.98 * (1 / phi + 0.35 * (1 - 1 / phi)); // MOCK combustion efficiency
    const qTotal = mFuel * f.lhv * etaC;
    // Cycle-to-cycle variability (always draw so the sequence is path-independent).
    const g1 = this.gauss();
    const g2 = this.gauss();
    const sl0 = sl0Of(f, phi);
    const slRatio = SL_REF / Math.max(sl0, 0.02);
    const ccv = this.ccvScale;
    const delay = 1.0e-3 * Math.pow(slRatio, 0.7) * (1 + 2 * y0) * Math.max(0.3, 1 + 0.08 * ccv * g2);
    const duration = clamp(
      48 * Math.sqrt(slRatio) * Math.pow(op.rpm / 600, 0.2) * (1 + 1.5 * y0) * Math.max(0.5, 1 + 0.05 * ccv * g1),
      20,
      150,
    );
    const mAir = (this.m * (1 - y0)) / (1 + far);
    const rhoAmb = op.ambientPressure / (R_AIR * op.ambientTemperature);
    this.plan = {
      mTrap: this.m,
      y0,
      yEgr0: this.yEgr,
      mFuel,
      qTotal,
      qPerMass: qTotal / this.m,
      pIvc: this.p,
      TIvc: T,
      VIvc: V,
      thetaIvc: this.theta,
      ignitable: sl0 > 0.05 && y0 < 0.45,
      delay,
      duration,
      sl0,
      volEff: mAir / (rhoAmb * this.kin.displacedVolume),
    };
    this.lw = 0;
    this.xNO = 0;
    // A spark that already happened while the valves were open cannot ignite this charge.
  }

  private onEvo(): void {
    if (this.twoZone) this.mergeZones();
    if (this.burnStarted) this.burnDone = true;
    const V = this.kin.volume(this.theta);
    const T = (this.p * V) / (this.m * R_GAS);
    if (this.burnStarted) {
      this.noPpm = this.xNO * this.xb * 1e6;
      this.coEvo = burnedGasComposition(this.fuel, this.op.equivalenceRatio, Math.max(T, 1800), this.p, this.xNO, this.scratchComp).CO * this.xb;
      this.xNOexhaust = this.xNO;
    } else {
      this.noPpm = 0;
      this.coEvo = 0;
      this.xNOexhaust = 0;
    }
    this.closed = false;
  }

  private gauss(): number {
    const u1 = Math.max(1e-12, this.rng());
    const u2 = this.rng();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  }

  private wrapCycle(): void {
    const pl = this.plan;
    const Vd = this.kin.displacedVolume;
    const fuelLhv = this.fuel.lhv;
    const mFuel = pl ? pl.mFuel : 0;
    const misfire = !this.ignited && this.xb < 0.5;
    const summary: CycleSummary = {
      cycle: this.cycle,
      imepGross: this.wGross / Vd,
      imepNet: this.wNet / Vd,
      pmep: (this.wGross - this.wNet) / Vd,
      peakPressure: this.peakP,
      peakPressureDeg: this.peakPDeg,
      maxPressureRiseRate: this.maxDpDeg,
      ca10: this.ca10,
      ca50: this.ca50,
      ca90: this.ca90,
      indicatedEfficiency: mFuel > 0 ? this.wNet / (mFuel * fuelLhv) : 0,
      isfc: this.wNet > 0 ? mFuel / this.wNet : NaN,
      trappedMass: pl ? pl.mTrap : 0,
      residualFraction: pl ? Math.max(0, pl.y0 - pl.yEgr0) : 0,
      volumetricEfficiency: pl ? pl.volEff : 0,
      fuelMass: this.burnStarted ? mFuel * Math.min(1, this.xb) : 0,
      noPpm: this.noPpm,
      coFraction: this.coEvo,
      knockOnsetDeg: this.knockOnsetDeg,
      knockEndGasFraction: this.knockEndGas,
      mapo: this.mapo,
      misfire,
      heatLoss: this.heatLossClosed,
      indicatedWorkGross: this.wGross,
    };
    this.summaries.push(summary);

    // New cycle: apply the pending operating point.
    const crOld = this.op.compressionRatio;
    this.op = cloneOp(this.pendingOp);
    this.fuel = mockFuel(this.op.fuel);
    if (this.op.compressionRatio !== crOld) {
      // Keep mass and temperature; p follows the new volume.
      const T = (this.p * this.kin.volume(360)) / (this.m * R_GAS);
      this.kin = new MockKinematics(this.spec.geometry, this.op.compressionRatio);
      this.p = (this.m * R_GAS * T) / this.kin.volume(-360);
    }
    this.theta = -360;
    this.cycle += 1;
    this.startCycle();
  }

  // -------------------------------------------------------------------------
  // Output
  // -------------------------------------------------------------------------

  private flameStage(radius: number): FlameStage {
    if (!this.ignited && !this.burnStarted) return 'none';
    if (this.burnDone) return 'done';
    if (this.xb >= XB_BURNOUT) return 'burnout';
    if (!this.burnStarted || radius < TURBULENT_STAGE_RADIUS) return 'kernel';
    return 'turbulent';
  }

  private kernelRadius(): number {
    const sp = this.spec.sparkPlug;
    if (!this.ignited || !this.plan) return 0;
    const r0 = 0.5 * sp.gap + 0.2e-3;
    const u = clamp((this.t - this.tKernel) / Math.max(this.plan.delay, 1e-6), 0, 1);
    return r0 + (KERNEL_HANDOFF_RADIUS - r0) * Math.pow(u, 0.6);
  }

  private sparkState(out: EngineSnapshot['spark'], rho: number): void {
    const ig = this.spec.ignition;
    const t = this.t;
    out.breakdownVoltage = breakdownVoltage(this.spec.sparkPlug.gap, rho);
    out.primaryCurrent = 0;
    out.secondaryVoltage = 0;
    out.secondaryCurrent = 0;
    out.energyDelivered = 0;
    let phase: SparkPhase = 'off';
    const sp = this.spark;
    if (!sp) {
      if (this.dwellStarted) {
        phase = 'charging';
        const R1 = Math.max(ig.primaryResistance, 1e-6);
        const tau1 = ig.primaryInductance / R1;
        const iRaw = (ig.supplyVoltage / R1) * (1 - Math.exp(-(t - this.tDwell) / tau1));
        out.primaryCurrent = Math.min(ig.primaryCurrentLimit, iRaw);
        const vL1 = iRaw >= ig.primaryCurrentLimit ? 0 : ig.supplyVoltage - R1 * iRaw;
        out.secondaryVoltage = ig.couplingCoefficient * Math.sqrt(ig.secondaryInductance / ig.primaryInductance) * vL1;
      }
      out.phase = phase;
      return;
    }
    const tau = t - sp.tSwitch;
    const L2 = ig.secondaryInductance;
    const R2 = ig.secondaryResistance;
    if (!sp.broke) {
      if (tau < Math.PI / sp.omega2) {
        phase = 'breakdown';
        out.secondaryVoltage = sp.vMax * Math.sin(sp.omega2 * tau);
        out.secondaryCurrent = sp.i20 * Math.cos(sp.omega2 * tau);
      } else {
        phase = 'done';
      }
    } else if (t < sp.tBd) {
      phase = 'breakdown';
      out.secondaryVoltage = sp.vMax * Math.sin(sp.omega2 * tau);
      out.secondaryCurrent = sp.i20 * Math.cos(sp.omega2 * tau);
    } else if (t < sp.tArcStart) {
      phase = 'breakdown';
      const u = t - sp.tBd;
      out.secondaryVoltage = ARC_VOLTAGE + (sp.vBd - ARC_VOLTAGE) * Math.exp(-u / 20e-9);
      out.secondaryCurrent = sp.i2AtBd + (sp.vBd / BREAKDOWN_SURGE_IMPEDANCE) * Math.exp(-u / 50e-9);
      out.energyDelivered = sp.eBreakdown;
    } else if (t < sp.tArcEnd) {
      phase = 'arc';
      const u = t - sp.tArcStart;
      out.secondaryVoltage = ARC_VOLTAGE;
      out.secondaryCurrent = dischargeCurrent(sp.i2AtBd, ARC_VOLTAGE, R2, L2, u);
      out.energyDelivered = sp.eBreakdown + ARC_VOLTAGE * dischargeCharge(sp.i2AtBd, ARC_VOLTAGE, R2, L2, u);
    } else if (t < sp.tGlowEnd) {
      phase = 'glow';
      const u = t - sp.tArcEnd;
      out.secondaryVoltage = sp.vGlow;
      out.secondaryCurrent = dischargeCurrent(sp.iArcEnd, sp.vGlow, R2, L2, u);
      out.energyDelivered = sp.eBreakdown + sp.eArc + sp.vGlow * dischargeCharge(sp.iArcEnd, sp.vGlow, R2, L2, u);
    } else {
      phase = 'done';
      out.energyDelivered = sp.eBreakdown + sp.eArc + sp.eGlow;
    }
    out.secondaryVoltage = Math.abs(out.secondaryVoltage);
    out.secondaryCurrent = Math.abs(out.secondaryCurrent);
    out.phase = phase;
  }

  private makeSnapshot(): EngineSnapshot {
    const spec = this.spec;
    const op = this.op;
    const th = this.theta;
    const kin = this.kin;
    const x = kin.pistonDisplacement(th);
    const h = kin.clearanceHeight(th);
    const V = kin.volume(th);
    const Tmean = (this.p * V) / (this.m * R_GAS);
    const osc = this.ringing ? this.knockOscillation() : 0;
    const pressure = this.p + osc;

    // Zone temperatures (T_mean = (1 − xb) Tu + xb Tb holds exactly in every mode).
    let Tu: number;
    let Tb: number;
    let xbOut: number;
    if (this.twoZone) {
      Tu = this.Tu;
      Tb = this.mb > 0 ? this.Tb : 0;
      xbOut = this.xb;
    } else if (this.burnDone) {
      Tu = Tmean;
      Tb = Tmean;
      xbOut = this.xb;
    } else {
      Tu = Tmean;
      Tb = 0;
      xbOut = 0;
    }

    // Flame.
    const R = 0.5 * spec.geometry.bore;
    const c = spec.sparkPlug.gapCenter;
    let radius = 0;
    let area = 0;
    if (this.twoZone && this.mb > 0) {
      const Vb = (this.mb * R_GAS * this.Tb) / this.p;
      radius = Math.max(this.kernelRadius(), sphereRadiusForVolume(Vb, c, R, h));
    } else if (this.burnDone) {
      radius = sphereCoverRadius(c, R, h);
    } else if (this.ignited) {
      radius = this.kernelRadius();
    }
    const stage = this.flameStage(radius);
    if (stage === 'kernel' || stage === 'turbulent' || stage === 'burnout') {
      area = sphereDiscFrontArea(radius, c, R, h);
    }
    const phi = op.equivalenceRatio;
    const dil = this.plan && this.closed ? this.plan.y0 : this.yProd;
    // No unburned gas left once the burn is done: no front, no burning velocities.
    const SL = this.burnDone ? 0 : laminarFlameSpeed(this.fuel, phi, Tu, this.p, dil);
    const Sp = (2 * spec.geometry.stroke * op.rpm) / 60;
    let uPrime: number;
    if (this.closed && this.plan) {
      const rhoIvc = this.plan.mTrap / this.plan.VIvc;
      uPrime = 0.27 * Sp * Math.cbrt(this.m / V / rhoIvc);
    } else {
      const aIv = valveEffectiveArea(spec.intakeValve, this.iv.lift(th), CD_INTAKE);
      const vJet = aIv > 0 ? Math.min(150, Math.abs(this.mdotIv) / ((this.m / V) * aIv)) : 0;
      uPrime = 0.5 * Sp + 0.08 * vJet;
    }
    let ST: number;
    if (this.burnDone) {
      ST = 0;
    } else if (stage === 'turbulent' || stage === 'burnout') {
      const rhoU = this.p / (R_GAS * Tu);
      ST = area > 1e-8 ? clamp(this.mbRate / (rhoU * area), SL, 30) : SL + uPrime;
    } else if (stage === 'kernel') {
      ST = SL * (1 + (0.5 * radius) / KERNEL_HANDOFF_RADIUS);
    } else {
      ST = SL + uPrime;
    }

    // Phase.
    let phase: CylinderPhase;
    const flameActive = stage === 'kernel' || stage === 'turbulent' || stage === 'burnout';
    if (!this.closed) phase = 'gas-exchange';
    else if (flameActive) phase = 'combustion';
    else if (this.burnDone || th > 0) phase = 'expansion';
    else phase = 'compression';

    // Composition.
    const comp: Composition = {
      CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0,
    };
    if (this.twoZone && this.mb > 0) {
      burnedGasComposition(this.fuel, phi, this.Tb, this.p, this.xNO, comp);
    } else if (this.burnDone) {
      burnedGasComposition(this.fuel, phi, Tmean, this.p, this.xNO, comp);
    } else {
      unburnedMixtureComposition(this.fuel, phi, this.yProd, this.xNOexhaust, comp);
    }

    // Spark.
    const spark: EngineSnapshot['spark'] = {
      phase: 'off',
      primaryCurrent: 0,
      secondaryVoltage: 0,
      secondaryCurrent: 0,
      energyDelivered: 0,
      breakdownVoltage: 0,
    };
    this.sparkState(spark, this.m / V);

    // Mechanics (fixed speed: α = 0).
    const omega = (op.rpm * 2 * Math.PI) / 60;
    const dxdth = kin.dxdTheta(th);
    const gasTorque = (pressure - op.ambientPressure) * kin.pistonArea * dxdth;
    const ms = spec.masses;
    const mRec = ms.piston + (ms.conRod * ms.conRodCgFromBigEnd) / spec.geometry.conRodLength;
    const inertiaTorque = -mRec * kin.d2xdTheta2(th) * omega * omega * dxdth;
    const fmep = 0.9e5 + 0.15e5 * (op.rpm / 1000); // MOCK friction MEP, Pa
    const frictionTorque = (fmep * kin.displacedVolume) / (4 * Math.PI);

    const pInt = (this.mi * R_GAS * op.intakeMixtureTemperature) / spec.manifolds.intakeVolume;
    const pExh = (this.me * R_GAS * this.Te) / spec.manifolds.exhaustVolume;

    return {
      t: this.t,
      cycle: this.cycle,
      thetaDeg: th,
      rpm: op.rpm,
      pistonDisplacement: x,
      clearanceHeight: h,
      rodAngle: kin.rodAngle(th),
      intakeLift: this.iv.lift(th),
      exhaustLift: this.ev.lift(th),
      phase,
      volume: V,
      pressure,
      temperatureMean: Tmean,
      temperatureUnburned: Tu,
      temperatureBurned: Tb,
      massFractionBurned: xbOut,
      mass: this.m,
      heatReleaseRate: this.hrr,
      heatLossRate: this.qLossRate,
      flame: {
        stage,
        radius,
        center: [c[0], c[1], c[2]],
        area,
        laminarSpeed: SL,
        turbulentSpeed: ST,
        turbulenceIntensity: uPrime,
      },
      spark,
      intakeMassFlow: this.closed || this.iv.lift(th) === 0 ? 0 : this.mdotIv,
      exhaustMassFlow: this.closed || this.ev.lift(th) === 0 ? 0 : this.mdotEv,
      intakeManifoldPressure: pInt,
      exhaustManifoldPressure: pExh,
      knock: {
        integral: this.lw,
        autoignited: this.knockOn,
        oscillation: osc,
      },
      burnedComposition: comp,
      gasTorque,
      netTorque: gasTorque + inertiaTorque - frictionTorque,
    };
  }
}

function cloneOp(op: OperatingPoint): OperatingPoint {
  return { ...op, fuel: { ...op.fuel } };
}
