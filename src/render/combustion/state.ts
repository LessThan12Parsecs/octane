/**
 * Per-frame derivation of everything the in-cylinder visuals draw, from the
 * EngineSnapshot only (+ EngineSpec geometry). Pure TypeScript, allocation-free
 * per frame, unit tested. The three.js views just copy these numbers into
 * uniforms / buffers.
 *
 * Radiometry: each light source's physical radiance (W m⁻² sr⁻¹) is turned
 * into luminance/colour through its emission spectrum (colour/emitters.ts)
 * and then scaled by the single exposure constant VISUAL_GAIN.
 */
import type { EngineSpec } from '../../physics/core/engine-spec';
import { ignitionSecondaryCapacitance } from '../../physics/core/engine-spec';
import type { EngineSnapshot, SparkPhase } from '../../physics/core/snapshot';
import { blackbodyXYZFast, xyzToLinearSrgb, clipToGamut, type Vec3 } from './colour/cie';
import {
  arcRGBPerWatt,
  autoignitionRGBPerWatt,
  breakdownRGBPerWatt,
  equivalenceRatioFromBurned,
  flameRGBPerWatt,
  glowRGBPerWatt,
  luminance,
} from './colour/emitters';
import {
  ARC_CURRENT_DENSITY,
  ARC_TEMPERATURE,
  AUTOIGNITION_CHEMILUMINESCENCE_EFFICIENCY,
  BREAKDOWN_RADIATIVE_EFFICIENCY,
  BREAKDOWN_TEMPERATURE,
  BURNED_GAS_VISIBLE_KAPPA,
  CHEMILUMINESCENCE_EFFICIENCY,
  GLOW_COLUMN_RADIUS,
  GLOW_REFERENCE_CURRENT,
  OVERLAY_KNOCK_FULL_SCALE_PA,
  SPARK_RADIATIVE_EFFICIENCY,
  VISUAL_GAIN,
  VIS_FLASH_PERSISTENCE_S,
  VIS_GLARE_FRACTION,
  VIS_GLARE_RADIUS,
  VIS_MIN_CHANNEL_RADIUS,
} from './constants';
import { GasFlowModel } from './flow';
import {
  brushParams,
  emptyBrush,
  knockModeFrequency,
  molarMassOf,
  soundSpeed,
  type BrushParams,
  type V3,
} from './geometry';

/** Minimum front area used when dividing heat release by area, m² (≈ a 0.3 mm kernel). */
export const MIN_FRONT_AREA = 1e-6;
/** Fastest e-folding time of the flame front's areal heat-release rate, simulated s (rate limiter). */
export const FRONT_FLUX_TIME_SCALE = 5e-4;
/** Minimum simulated exposure used for flash integration when paused/stepping, s. */
export const MIN_EXPOSURE_S = 1e-5;
/** Knock pattern: displayed as time-resolved only if a period spans ≥ ~7 frames. */
const KNOCK_MAX_CYCLES_PER_FRAME = 0.15;
/** Display clamp for emitters (keeps half-float targets finite). */
const DISPLAY_CLAMP = 5e3;
/** Typical sustaining voltages used only if the snapshot lacks a current (V). UNVERIFIED. */
const ARC_SUSTAIN_V = 100;
const GLOW_SUSTAIN_V = 450;

export type CombustionMode = 'physical' | 'temperature';

export interface KnockVisual {
  /** Autoignition has occurred this cycle. */
  active: boolean;
  /** Overlay strength 0..1 (envelope / OVERLAY_KNOCK_FULL_SCALE_PA). */
  amplitude: number;
  /** Envelope of |knock.oscillation|, Pa. */
  envelopePa: number;
  /** cos(2π f (t − t_onset)) when time-resolvable, else 0. */
  temporal: number;
  /** 1 when the oscillation is too fast to show at this playback speed (draw RMS pattern). */
  rms: number;
  /** Orientation of the (1,0) mode antinode, rad (toward the end-gas site). */
  axisAngle: number;
  /** End-gas autoignition site (x, z), m. */
  origin: [number, number];
  /** Radius of the expanding pressure front from the site, m. */
  ringRadius: number;
  ringWidth: number;
  /** 0 = pure expanding front, 1 = fully developed standing mode. */
  standing: number;
  /** (1,0) mode frequency, Hz; bulk sound speed, m/s. */
  frequency: number;
  soundSpeed: number;
  onsetTime: number;
}

export interface SparkVisual {
  visible: boolean;
  kind: 'none' | 'arc' | 'glow';
  /** Frame-averaged electrical power into the column, W. */
  power: number;
  /** Gap current, A. */
  current: number;
  /** Physical and drawn column radius, m. */
  radius: number;
  drawRadius: number;
  /** Electrode tips (plug side → ground electrode), cylinder frame, m. */
  start: V3;
  end: V3;
  /** Mid-point bow of the channel (convected by the flow), m. */
  bow: V3;
  /** Display radiance through the column centre (linear sRGB). */
  coreRadiance: V3;
  /** Glare halo display radiance (linear sRGB) and radius, m. */
  glare: V3;
  glareRadius: number;
  /** Breakdown flash glare, persisted (linear sRGB) and radius, m. */
  flash: V3;
  flashRadius: number;
}

export interface ChamberLight {
  /** three.js PointLight intensity (display units; = VISUAL_GAIN × luminous intensity in cd). */
  intensity: number;
  /** Normalised colour (max channel 1). */
  color: V3;
  position: V3;
}

const v3 = (): V3 => [0, 0, 0];

function setScaled(out: V3, rgb: ArrayLike<number>, k: number): V3 {
  out[0] = Math.min(rgb[0] * k, DISPLAY_CLAMP);
  out[1] = Math.min(rgb[1] * k, DISPLAY_CLAMP);
  out[2] = Math.min(rgb[2] * k, DISPLAY_CLAMP);
  return out;
}

export class CombustionVisualState {
  readonly flow: GasFlowModel;
  readonly R: number;
  readonly bore: number;

  // ---- time ----
  t = 0;
  /** Simulated seconds advanced by the last update. */
  dtSim = 0;
  /** True if the last update detected a discontinuity (time went backwards / reset). */
  wasReset = false;
  cycle = -1;

  // ---- chamber ----
  h = 0;
  volume = 0;

  // ---- flame / gas field ----
  flameVisible = false;
  /** Whole chamber is one zone (burned after combustion, or mixed during gas exchange). */
  allBurned = false;
  flameCenter: V3 = v3();
  flameRadius = 0;
  readonly brush: BrushParams = emptyBrush();
  /** Accumulated eddy-turnover phase ∫u'/L_I dt (drives the wrinkle noise). */
  noisePhase = 0;
  /** Front emission scale (display radiance for a face-on front), linear sRGB. */
  frontJ: V3 = v3();
  /** Burned-gas emission coefficient, display units per m. */
  burnedJ: V3 = v3();
  /** End-gas autoignition emission coefficient, display units per m. */
  endGasJ: V3 = v3();
  /** Haze visibility of the unburned region / of the "core" region (0..1). */
  hazeUnburned = 1;
  hazeCore = 0;
  /** Zone temperatures for the false-colour mode, K. */
  Tu = 300;
  Tb = 300;
  /** Equivalence ratio estimate from the burned composition (default 1). */
  phi = 1;
  /** Front heat release attributed to the flame / end gas, W. */
  hrrFront = 0;
  hrrEndGas = 0;
  /** Areal heat release rate of the front, W/m². */
  frontFlux = 0;

  readonly knock: KnockVisual = {
    active: false, amplitude: 0, envelopePa: 0, temporal: 0, rms: 0, axisAngle: Math.PI,
    origin: [0, 0], ringRadius: 0, ringWidth: 0, standing: 1, frequency: 0, soundSpeed: 0, onsetTime: 0,
  };

  readonly spark: SparkVisual;
  readonly light: ChamberLight = { intensity: 0, color: [1, 1, 1], position: v3() };

  // ---- internals ----
  private hasPrev = false;
  private lastPhase: SparkPhase = 'off';
  private lastEnergy = 0;
  private lastKind: 'none' | 'arc' | 'glow' = 'none';
  private pendingFlashJ = 0;
  private flashCd = 0;
  private qRef = 0;
  private readonly gapAxis: V3;
  private readonly gapCenter: V3;
  private readonly gap: number;
  private readonly secC: number;
  private readonly tmp: V3 = v3();
  private readonly tmp2: V3 = v3();
  private readonly turb: V3 = v3();
  private rngState = 0x9e3779b9;

  constructor(private readonly spec: EngineSpec) {
    this.flow = new GasFlowModel(spec);
    this.bore = spec.geometry.bore;
    this.R = this.bore / 2;
    const sp = spec.sparkPlug;
    const al = Math.hypot(sp.axis[0], sp.axis[1], sp.axis[2]) || 1;
    this.gapAxis = [sp.axis[0] / al, sp.axis[1] / al, sp.axis[2] / al];
    this.gapCenter = [sp.gapCenter[0], sp.gapCenter[1], sp.gapCenter[2]];
    this.gap = sp.gap;
    this.secC = ignitionSecondaryCapacitance(spec.ignition);
    const g2 = this.gap / 2;
    this.spark = {
      visible: false, kind: 'none', power: 0, current: 0, radius: 0, drawRadius: VIS_MIN_CHANNEL_RADIUS,
      start: [this.gapCenter[0] - this.gapAxis[0] * g2, this.gapCenter[1] - this.gapAxis[1] * g2, this.gapCenter[2] - this.gapAxis[2] * g2],
      end: [this.gapCenter[0] + this.gapAxis[0] * g2, this.gapCenter[1] + this.gapAxis[1] * g2, this.gapCenter[2] + this.gapAxis[2] * g2],
      bow: v3(), coreRadiance: v3(), glare: v3(), glareRadius: VIS_GLARE_RADIUS, flash: v3(), flashRadius: 1.6 * VIS_GLARE_RADIUS,
    };
  }

  /** Forget history (seek backwards / reset). */
  reset(): void {
    this.hasPrev = false;
    this.flow.reset();
    this.noisePhase = 0;
    this.lastPhase = 'off';
    this.lastEnergy = 0;
    this.lastKind = 'none';
    this.pendingFlashJ = 0;
    this.flashCd = 0;
    this.qRef = 0;
    this.knock.active = false;
    this.knock.envelopePa = 0;
    this.knock.amplitude = 0;
    this.spark.bow[0] = this.spark.bow[1] = this.spark.bow[2] = 0;
    this.spark.flash[0] = this.spark.flash[1] = this.spark.flash[2] = 0;
    this.cycle = -1;
  }

  /**
   * @param s        interpolated snapshot for this frame
   * @param dtWall   wall-clock seconds since the previous frame
   * @param timeScale simulated seconds per wall second
   */
  update(s: EngineSnapshot, dtWall: number, timeScale: number): void {
    let dtSim = this.hasPrev ? s.t - this.t : 0;
    this.wasReset = false;
    if (dtSim < 0 || (this.hasPrev && s.cycle < this.cycle)) {
      this.reset();
      this.wasReset = true;
      dtSim = 0;
    }
    this.dtSim = dtSim;
    this.t = s.t;
    this.cycle = s.cycle;
    this.hasPrev = true;
    this.h = s.clearanceHeight;
    this.volume = s.volume;

    this.flow.update(s, dtSim);
    this.updatePhi(s);
    this.updateKnock(s, dtWall, timeScale, dtSim);
    this.updateFlame(s, dtSim);
    this.updateSpark(s, dtWall, timeScale, dtSim);
    this.updateLight(s);
  }

  // -------------------------------------------------------------------------

  private updatePhi(s: EngineSnapshot): void {
    if (s.massFractionBurned > 0.01 && s.temperatureBurned > 0) {
      const phi = equivalenceRatioFromBurned(s.burnedComposition);
      if (Number.isFinite(phi) && phi > 0.2 && phi < 3) this.phi = phi;
    }
  }

  private updateFlame(s: EngineSnapshot, dtSim: number): void {
    const f = s.flame;
    const gasExchange = s.phase === 'gas-exchange';
    this.flameVisible = !gasExchange && (f.stage === 'kernel' || f.stage === 'turbulent' || f.stage === 'burnout') && f.radius > 0;
    this.allBurned = gasExchange || (!this.flameVisible && (f.stage === 'done' || s.massFractionBurned >= 0.999));
    this.flameCenter[0] = f.center[0];
    this.flameCenter[1] = f.center[1];
    this.flameCenter[2] = f.center[2];
    this.flameRadius = this.flameVisible ? f.radius : 0;

    const Tu = s.temperatureUnburned > 0 ? s.temperatureUnburned : s.temperatureMean;
    const TbRaw = s.temperatureBurned > 0 ? s.temperatureBurned : s.temperatureMean;
    if (gasExchange) {
      this.Tu = this.Tb = s.temperatureMean;
    } else {
      this.Tu = Tu;
      this.Tb = TbRaw;
    }

    brushParams({
      uPrime: f.turbulenceIntensity > 0 ? f.turbulenceIntensity : this.flow.uPrime,
      laminarSpeed: f.laminarSpeed,
      clearanceHeight: s.clearanceHeight,
      radius: this.flameRadius,
      pressure: s.pressure,
      unburnedTemperature: Tu,
    }, this.brush);
    if (this.brush.integralScale > 0 && dtSim > 0) this.noisePhase += (dtSim * this.flow.uPrime) / this.brush.integralScale;
    // keep the float small for the GPU
    if (this.noisePhase > 1e4) this.noisePhase -= 1e4;

    // ---- heat release partition: flame front vs end-gas autoignition ----
    // The front's areal heat-release rate q'' = ρ_u S_T Q evolves on the flame
    // time scale (≳ 0.5 ms); end-gas autoignition releases heat ~10× faster.
    // Track q'' with a rate limit in SIMULATED time and attribute any excess
    // to the end gas (robust to interpolated snapshots that blend a spike in
    // before the `autoignited` flag flips). Frozen once autoignition is flagged.
    const hrr = Math.max(0, s.heatReleaseRate);
    const A = Math.max(f.area, MIN_FRONT_AREA);
    if (!this.flameVisible) {
      if (!this.knock.active) this.qRef = 0;
    } else if (!this.knock.active && hrr > 0) {
      const q = hrr / A;
      if (!(this.qRef > 0) || dtSim <= 0) {
        if (!(this.qRef > 0)) this.qRef = q;
      } else {
        const g = 1 + dtSim / FRONT_FLUX_TIME_SCALE;
        this.qRef = Math.min(Math.max(q, this.qRef / g), this.qRef * g);
      }
    }
    this.hrrFront = this.flameVisible ? Math.min(hrr, this.qRef * A) : 0;
    this.hrrEndGas = hrr - this.hrrFront;
    if (!this.knock.active && !this.flameVisible) this.hrrEndGas = 0;
    this.frontFlux = this.flameVisible ? this.hrrFront / A : 0;

    // Front: radiance seen face-on = η q'' / (4π) (thin isotropic emitter).
    const flameRGB = flameRGBPerWatt(this.phi);
    setScaled(this.frontJ, flameRGB, VISUAL_GAIN * CHEMILUMINESCENCE_EFFICIENCY * this.frontFlux / (4 * Math.PI));

    // End gas: volumetric, spread over the unburned volume V_u (common p, ideal gas, equal M).
    const xb = Math.min(Math.max(s.massFractionBurned, 0), 1);
    const Vu = this.unburnedVolume(s, xb);
    const aiRGB = autoignitionRGBPerWatt();
    setScaled(this.endGasJ, aiRGB, Vu > 1e-9
      ? (VISUAL_GAIN * AUTOIGNITION_CHEMILUMINESCENCE_EFFICIENCY * this.hrrEndGas) / (4 * Math.PI * Vu)
      : 0);

    // Burned gas: gray body, j = κ B(T_b).
    const Tglow = gasExchange ? s.temperatureMean : (s.temperatureBurned > 0 ? s.temperatureBurned : (this.allBurned ? s.temperatureMean : 0));
    if (Tglow > 800) {
      blackbodyXYZFast(Tglow, this.tmp);
      clipToGamut(xyzToLinearSrgb(this.tmp, this.tmp2));
      setScaled(this.burnedJ, this.tmp2, VISUAL_GAIN * BURNED_GAS_VISIBLE_KAPPA);
    } else {
      this.burnedJ[0] = this.burnedJ[1] = this.burnedJ[2] = 0;
    }

    // Haze of the unburned charge (visual cue only).
    this.hazeUnburned = 1;
    this.hazeCore = gasExchange ? 1 - xb : 0;
  }

  /** V_u = V (1−x_b) T_u / ((1−x_b) T_u + x_b T_b). */
  unburnedVolume(s: EngineSnapshot, xb: number): number {
    const Tu = s.temperatureUnburned > 0 ? s.temperatureUnburned : s.temperatureMean;
    const Tb = s.temperatureBurned > 0 ? s.temperatureBurned : Tu;
    const den = (1 - xb) * Tu + xb * Tb;
    return den > 0 ? (s.volume * (1 - xb) * Tu) / den : 0;
  }

  private updateKnock(s: EngineSnapshot, dtWall: number, timeScale: number, dtSim: number): void {
    const k = this.knock;
    const c = soundSpeed(s.temperatureMean > 0 ? s.temperatureMean : 300, molarMassOf(s.burnedComposition));
    k.soundSpeed = c;
    k.frequency = knockModeFrequency(c, this.bore);

    if (s.knock.autoignited && !k.active) {
      k.active = true;
      k.onsetTime = s.t;
      // End-gas site: the liner point farthest from the flame centre.
      const cx = s.flame.center[0], cz = s.flame.center[2];
      const n = Math.hypot(cx, cz);
      const dx = n > 1e-6 ? -cx / n : -1, dz = n > 1e-6 ? -cz / n : 0;
      k.origin[0] = dx * this.R;
      k.origin[1] = dz * this.R;
      k.axisAngle = Math.atan2(dz, dx);
    } else if (!s.knock.autoignited && k.active) {
      k.active = false;
    }

    const osc = Math.abs(s.knock.oscillation);
    const tauEnv = k.frequency > 0 ? 3 / k.frequency : 5e-4;
    k.envelopePa = Math.max(osc, dtSim > 0 ? k.envelopePa * Math.exp(-dtSim / tauEnv) : k.envelopePa);
    if (k.envelopePa < 1) k.envelopePa = 0;
    k.amplitude = Math.min(1, k.envelopePa / OVERLAY_KNOCK_FULL_SCALE_PA);

    const since = k.active ? Math.max(s.t - k.onsetTime, 0) : Infinity;
    k.ringRadius = Number.isFinite(since) ? c * since : 4 * this.R;
    k.ringWidth = 0.12 * this.R;
    k.standing = Math.min(1, k.ringRadius / (2 * this.R));
    const cyclesPerFrame = k.frequency * timeScale * dtWall;
    if (cyclesPerFrame <= KNOCK_MAX_CYCLES_PER_FRAME) {
      k.rms = 0;
      k.temporal = Number.isFinite(since) ? Math.cos(2 * Math.PI * k.frequency * since) : 0;
    } else {
      k.rms = 1;
      k.temporal = 0;
    }
  }

  private rand(): number {
    // mulberry32
    let t = (this.rngState = (this.rngState + 0x6d2b79f5) | 0);
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  }

  private gauss(): number {
    const u = Math.max(this.rand(), 1e-12), v = this.rand();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
  }

  private updateSpark(s: EngineSnapshot, dtWall: number, timeScale: number, dtSim: number): void {
    const sp = this.spark;
    const phase = s.spark.phase;
    const discharging = phase === 'breakdown' || phase === 'arc' || phase === 'glow';
    const wasDischarging = this.lastPhase === 'breakdown' || this.lastPhase === 'arc' || this.lastPhase === 'glow';

    // ---- breakdown events (restrikes included) ----
    let breakdownElec = 0;
    if ((phase === 'breakdown' && this.lastPhase !== 'breakdown') || (discharging && !wasDischarging)) {
      const V = Math.max(s.spark.breakdownVoltage, 0);
      breakdownElec = 0.5 * this.secC * V * V;
      this.pendingFlashJ += BREAKDOWN_RADIATIVE_EFFICIENCY * breakdownElec;
      sp.bow[0] = sp.bow[1] = sp.bow[2] = 0;
      this.turb[0] = this.turb[1] = this.turb[2] = 0;
    }

    // ---- frame-averaged column power (what a camera integrates over one frame) ----
    // Uses the cumulative energy counter when the simulator fills it (true frame
    // average, so a 1 ms discharge in a 16 ms real-time frame shows ~1/16 as bright);
    // falls back to instantaneous V·I when paused or when the counter is absent.
    const E = s.spark.energyDelivered;
    let dE = E - this.lastEnergy;
    if (dE < 0) dE = Math.max(E, 0);
    const counterLive = E > 0 || this.lastEnergy > 0;
    this.lastEnergy = E;
    const Pinst = discharging ? Math.abs(s.spark.secondaryVoltage * s.spark.secondaryCurrent) : 0;
    let P: number;
    if (dtSim > 1e-12 && counterLive) {
      // The breakdown's capacitive energy is drawn as the flash, not as column glow.
      P = Math.max(0, dE - breakdownElec) / dtSim;
    } else {
      P = Pinst;
    }

    let kind: 'none' | 'arc' | 'glow' = phase === 'glow' ? 'glow' : phase === 'arc' || phase === 'breakdown' ? 'arc' : 'none';
    if (kind === 'none' && P > 0) kind = this.lastKind === 'none' ? 'glow' : this.lastKind;
    if (kind !== 'none') this.lastKind = kind;
    if (!discharging && P <= 0) this.lastKind = 'none';

    let I = Math.abs(s.spark.secondaryCurrent);
    if (!(I > 0) && P > 0) I = P / (kind === 'arc' ? ARC_SUSTAIN_V : GLOW_SUSTAIN_V);
    sp.current = I;
    sp.power = P;
    sp.kind = P > 0 ? kind : 'none';
    sp.visible = P > 0 && kind !== 'none';

    const a = kind === 'arc'
      ? Math.max(Math.sqrt(I / (Math.PI * ARC_CURRENT_DENSITY)), 2e-5)
      : GLOW_COLUMN_RADIUS * Math.min(Math.max(Math.sqrt(I / GLOW_REFERENCE_CURRENT), 0.5), 2);
    sp.radius = a;
    sp.drawRadius = Math.max(a, VIS_MIN_CHANNEL_RADIUS);

    // ---- channel convection by the gas flow (bow), with restrike ----
    if (sp.visible && dtSim > 0) {
      const u = this.flow.velocityAt(this.gapCenter[0], this.gapCenter[1], this.gapCenter[2], this.tmp);
      const up = this.flow.uPrime;
      const TL = up > 0 ? this.flow.integralScale / up : 1;
      const e = Math.exp(-dtSim / TL);
      const sq = up * Math.sqrt(Math.max(0, 1 - e * e));
      for (let i = 0; i < 3; i++) this.turb[i] = this.turb[i] * e + sq * this.gauss();
      const ux = u[0] + this.turb[0], uy = u[1] + this.turb[1], uz = u[2] + this.turb[2];
      const ax = this.gapAxis;
      const along = ux * ax[0] + uy * ax[1] + uz * ax[2];
      sp.bow[0] += (ux - along * ax[0]) * dtSim;
      sp.bow[1] += (uy - along * ax[1]) * dtSim;
      sp.bow[2] += (uz - along * ax[2]) * dtSim;
      // Channel stretched too far → restrike (short channel again). Approximate criterion.
      if (Math.hypot(sp.bow[0], sp.bow[1], sp.bow[2]) > 1.5 * this.gap) {
        sp.bow[0] = sp.bow[1] = sp.bow[2] = 0;
      }
    } else if (!sp.visible) {
      sp.bow[0] = sp.bow[1] = sp.bow[2] = 0;
    }
    const d = Math.hypot(sp.bow[0], sp.bow[1], sp.bow[2]);
    const len = this.gap * Math.sqrt(1 + (16 / 3) * (d / this.gap) ** 2);

    // ---- column radiance and glare ----
    const rgbPW = kind === 'arc' ? arcRGBPerWatt(ARC_TEMPERATURE) : glowRGBPerWatt();
    const lumPW = luminance(rgbPW);
    if (sp.visible) {
      const L = (SPARK_RADIATIVE_EFFICIENCY * P) / (2 * Math.PI * Math.PI * sp.drawRadius * len);
      setScaled(sp.coreRadiance, rgbPW, VISUAL_GAIN * L);
      const Icd = (lumPW * SPARK_RADIATIVE_EFFICIENCY * P) / (4 * Math.PI);
      const k = lumPW > 0 ? (VISUAL_GAIN * VIS_GLARE_FRACTION * Icd) / (Math.PI * sp.glareRadius ** 2) / lumPW : 0;
      setScaled(sp.glare, rgbPW, k);
    } else {
      sp.coreRadiance[0] = sp.coreRadiance[1] = sp.coreRadiance[2] = 0;
      sp.glare[0] = sp.glare[1] = sp.glare[2] = 0;
    }

    // ---- breakdown flash: integrate over the frame's simulated exposure, then persist ----
    const bdRGB = breakdownRGBPerWatt(BREAKDOWN_TEMPERATURE);
    const bdLum = luminance(bdRGB);
    let newCd = 0;
    if (this.pendingFlashJ > 0) {
      const exposure = Math.max(dtWall * timeScale, MIN_EXPOSURE_S);
      newCd = (bdLum * this.pendingFlashJ) / exposure / (4 * Math.PI);
      this.pendingFlashJ = 0;
    }
    if (dtWall > 0) {
      const e = Math.exp(-dtWall / VIS_FLASH_PERSISTENCE_S);
      this.flashCd = this.flashCd * e + newCd * (1 - e);
    } else {
      this.flashCd = Math.max(this.flashCd, newCd);
    }
    if (this.flashCd < 1e-9) this.flashCd = 0;
    const kf = bdLum > 0 ? (VISUAL_GAIN * VIS_GLARE_FRACTION * this.flashCd) / (Math.PI * sp.flashRadius ** 2) / bdLum : 0;
    setScaled(sp.flash, bdRGB, kf);

    this.lastPhase = phase;
  }

  /** c += I · rgb / luminance(rgb)  (unit-luminance chromaticity weighted by intensity). */
  private accumulateChroma(c: V3, rgb: ArrayLike<number>, I: number): void {
    const l = luminance(rgb);
    if (!(I > 0) || !(l > 0)) return;
    c[0] += (rgb[0] / l) * I;
    c[1] += (rgb[1] / l) * I;
    c[2] += (rgb[2] / l) * I;
  }

  private updateLight(s: EngineSnapshot): void {
    // Luminous intensities (cd) of the sources and their chromaticities.
    const flameRGB = flameRGBPerWatt(this.phi);
    const lf = luminance(flameRGB);
    const Iflame = (lf * CHEMILUMINESCENCE_EFFICIENCY * this.hrrFront) / (4 * Math.PI);
    const aiRGB = autoignitionRGBPerWatt();
    const la = luminance(aiRGB);
    const Iai = (la * AUTOIGNITION_CHEMILUMINESCENCE_EFFICIENCY * this.hrrEndGas) / (4 * Math.PI);
    const xb = Math.min(Math.max(s.massFractionBurned, 0), 1);
    const Vb = Math.max(s.volume - this.unburnedVolume(s, xb), 0);
    // burnedJ is display per metre = VISUAL_GAIN κ L_bb → I_b = κ L_bb V_b (cd)
    const lb = luminance(this.burnedJ) / VISUAL_GAIN;
    const Iburned = lb * (this.allBurned ? s.volume : Vb);
    const sp = this.spark;
    const sparkRGB = sp.kind === 'arc' ? arcRGBPerWatt(ARC_TEMPERATURE) : glowRGBPerWatt();
    const ls = luminance(sparkRGB);
    const Ispark = sp.visible ? (ls * SPARK_RADIATIVE_EFFICIENCY * sp.power) / (4 * Math.PI) : 0;
    const bdRGB = breakdownRGBPerWatt(BREAKDOWN_TEMPERATURE);
    const Iflash = this.flashCd;

    const c = this.tmp;
    c[0] = c[1] = c[2] = 0;
    this.accumulateChroma(c, flameRGB, Iflame);
    this.accumulateChroma(c, aiRGB, Iai);
    this.accumulateChroma(c, this.burnedJ, Iburned);
    this.accumulateChroma(c, sparkRGB, Ispark);
    this.accumulateChroma(c, bdRGB, Iflash);
    const m = Math.max(c[0], c[1], c[2]);
    const L = this.light;
    if (m > 0) {
      L.color[0] = c[0] / m;
      L.color[1] = c[1] / m;
      L.color[2] = c[2] / m;
      L.intensity = Math.min(VISUAL_GAIN * m, DISPLAY_CLAMP);
    } else {
      L.intensity = 0;
    }

    // Position: intensity-weighted between the burned-gas region and the spark gap.
    const Igas = Iflame + Iai + Iburned;
    const Igap = Ispark + Iflash;
    const px = this.allBurned ? 0 : this.flameCenter[0];
    const pz = this.allBurned ? 0 : this.flameCenter[2];
    const rr = Math.hypot(px, pz);
    const lim = 0.8 * this.R;
    const f = rr > lim ? lim / rr : 1;
    const gx = px * f, gz = pz * f, gy = -0.5 * this.h;
    const w = Igas + Igap > 0 ? Igap / (Igas + Igap) : 1;
    L.position[0] = gx + (this.gapCenter[0] - gx) * w;
    L.position[1] = gy + (this.gapCenter[1] - gy) * w;
    L.position[2] = gz + (this.gapCenter[2] - gz) * w;
  }
}

export type { Vec3 };
