import { describe, expect, it } from 'vitest';
import { flameRGBPerWatt, glowRGBPerWatt, luminance } from './colour/emitters';
import {
  CHEMILUMINESCENCE_EFFICIENCY,
  OVERLAY_KNOCK_FULL_SCALE_PA,
  VISUAL_GAIN,
  VIS_MIN_CHANNEL_RADIUS,
} from './constants';
import { MODEL_T } from '../../physics/engines/model-t';
import { chamberShapeOf, insideChamberShape } from './chamber';
import { axialModeFrequency, axialModeShape, knockModeFrequency } from './geometry';
import { CombustionVisualState } from './state';
import { flameSnap, modelTSnap, snap, testSpec, TEST_R } from './test-utils';

const fresh = () => new CombustionVisualState(testSpec());

describe('flame front emission', () => {
  it('no flame → no front emission', () => {
    const st = fresh();
    st.update(snap(), 0.016, 1e-3);
    expect(st.flameVisible).toBe(false);
    expect(st.frontJ).toEqual([0, 0, 0]);
    expect(st.light.intensity).toBe(0);
  });

  it('face-on front radiance = VISUAL_GAIN · rgb(φ) · η · (HRR/A) / 4π', () => {
    const st = fresh();
    const s = flameSnap({ heatReleaseRate: 1.2e5, flame: { area: 2e-3 } });
    st.update(s, 0.016, 1e-3);
    const phi = st.phi;
    const rgb = flameRGBPerWatt(phi);
    const k = (VISUAL_GAIN * CHEMILUMINESCENCE_EFFICIENCY * (1.2e5 / 2e-3)) / (4 * Math.PI);
    for (let i = 0; i < 3; i++) expect(st.frontJ[i]).toBeCloseTo(rgb[i] * k, 9);
    expect(st.frontFlux).toBeCloseTo(6e7, 3);
  });

  it('brightness ∝ heat release per front area', () => {
    const a = fresh(), b = fresh(), c = fresh();
    a.update(flameSnap({ heatReleaseRate: 1e5, flame: { area: 2e-3 } }), 0.016, 1e-3);
    b.update(flameSnap({ heatReleaseRate: 2e5, flame: { area: 2e-3 } }), 0.016, 1e-3);
    c.update(flameSnap({ heatReleaseRate: 2e5, flame: { area: 4e-3 } }), 0.016, 1e-3);
    expect(b.frontJ[2] / a.frontJ[2]).toBeCloseTo(2, 9);
    expect(c.frontJ[2] / a.frontJ[2]).toBeCloseTo(1, 9);
  });

  it('burned gas glows faintly (≪ the front) and warm (red > blue)', () => {
    const st = fresh();
    st.update(flameSnap(), 0.016, 1e-3);
    const pathFront = 1; // front: frontJ is already a radiance
    const burnedThroughChamber = st.burnedJ.map((v) => v * 0.05); // 5 cm of burned gas
    expect(luminance(burnedThroughChamber)).toBeLessThan(0.2 * luminance(st.frontJ) * pathFront);
    expect(st.burnedJ[0]).toBeGreaterThan(st.burnedJ[2]);
    expect(st.burnedJ[0]).toBeGreaterThan(0);
  });

  it('whole chamber is one zone after combustion and during gas exchange', () => {
    const st = fresh();
    st.update(flameSnap({ flame: { stage: 'done' }, massFractionBurned: 1 }), 0.016, 1e-3);
    expect(st.allBurned).toBe(true);
    expect(st.flameVisible).toBe(false);
    st.update(snap({ t: 1e-3, phase: 'gas-exchange', temperatureMean: 900 }), 0.016, 1e-3);
    expect(st.allBurned).toBe(true);
    expect(st.Tu).toBe(900);
    expect(st.Tb).toBe(900);
  });
});

describe('knock', () => {
  function runKnock(stepHRR: number[], flags: boolean[]) {
    const st = fresh();
    let t = 0;
    for (let i = 0; i < stepHRR.length; i++) {
      st.update(flameSnap({ t, heatReleaseRate: stepHRR[i], knock: { autoignited: flags[i], oscillation: flags[i] ? 1e5 : 0 } }), 0.016, 1e-4);
      t += 2e-6;
    }
    return st;
  }

  it('after autoignition the front keeps its areal rate; the excess is end-gas emission', () => {
    const st = runKnock([1e5, 1e5, 1e5, 2e6], [false, false, false, true]);
    expect(st.knock.active).toBe(true);
    expect(st.hrrFront).toBeCloseTo(1e5, 3);
    expect(st.hrrEndGas).toBeCloseTo(1.9e6, 0);
    expect(st.endGasJ[2]).toBeGreaterThan(0);
    // end-gas volumetric emission over a few cm outshines the front
    expect(luminance(st.endGasJ.map((v) => v * 0.02))).toBeGreaterThan(luminance(st.frontJ));
  });

  it('a spike that arrives before the autoignition flag (interpolated snapshot) is still end gas', () => {
    const st = runKnock([1e5, 1e5, 1e5, 2e6], [false, false, false, false]);
    expect(st.hrrFront).toBeLessThan(1.01e5);
    expect(st.hrrEndGas).toBeGreaterThan(1.8e6);
  });

  it('onset site is the liner point opposite the flame centre; envelope, frequency, temporal/RMS', () => {
    const st = fresh();
    st.update(flameSnap({ t: 0 }), 0.016, 1e-4);
    st.update(flameSnap({ t: 1e-5, knock: { autoignited: true, oscillation: -0.75 * OVERLAY_KNOCK_FULL_SCALE_PA } }), 0.016, 1e-4);
    const k = st.knock;
    expect(k.active).toBe(true);
    expect(k.origin[0]).toBeCloseTo(0, 12);
    expect(k.origin[1]).toBeCloseTo(-TEST_R, 12); // plug is at +z
    expect(k.axisAngle).toBeCloseTo(-Math.PI / 2, 12);
    expect(k.amplitude).toBeCloseTo(0.75, 9);
    expect(k.frequency).toBeGreaterThan(3000);
    expect(k.rms).toBe(0); // 1e-4 slow motion: ~0.01 periods per frame
    expect(k.ringRadius).toBe(0);
    // envelope decays after the oscillation stops (no new peaks)
    st.update(flameSnap({ t: 1e-5 + 5e-3, knock: { autoignited: true, oscillation: 0 } }), 0.016, 1e-4);
    expect(st.knock.envelopePa).toBeLessThan(0.01 * 0.75 * OVERLAY_KNOCK_FULL_SCALE_PA);
    expect(st.knock.standing).toBe(1);
    // real-time playback: ~100 periods per frame → RMS pattern
    const rt = fresh();
    rt.update(flameSnap({ t: 0, knock: { autoignited: true, oscillation: 1e5 } }), 0.016, 1);
    expect(rt.knock.rms).toBe(1);
  });
});

describe('spark', () => {
  it('breakdown flash triggers once per breakdown, integrates over the frame exposure and persists', () => {
    const run = (timeScale: number) => {
      const st = fresh();
      st.update(snap({ t: 0, spark: { phase: 'charging' } }), 0.016, timeScale);
      st.update(snap({ t: 1e-6, spark: { phase: 'breakdown', breakdownVoltage: 12e3 } }), 0.016, timeScale);
      return st;
    };
    const slow = run(1e-3), fast = run(1);
    expect(slow.spark.flash[2]).toBeGreaterThan(0);
    // same energy integrated over a 1000× longer exposure → 1000× dimmer
    expect(slow.spark.flash[2] / fast.spark.flash[2]).toBeCloseTo(1e3, 0);
    const f0 = slow.spark.flash[2];
    slow.update(snap({ t: 2e-6, spark: { phase: 'arc', secondaryVoltage: 100, secondaryCurrent: 0.1 } }), 0.016, 1e-3);
    expect(slow.spark.flash[2]).toBeLessThan(f0); // decays, no retrigger
    for (let i = 0; i < 60; i++) slow.update(snap({ t: 3e-6 + i * 1e-6, spark: { phase: 'done' } }), 0.016, 1e-3);
    expect(slow.spark.flash[2]).toBeLessThan(1e-3 * f0);
  });

  it('arc: thin column with radius ∝ √I; glow: diffuse violet column', () => {
    const st = fresh();
    st.update(snap({ t: 0, spark: { phase: 'charging' } }), 0.016, 1e-3);
    st.update(snap({ t: 0, spark: { phase: 'arc', secondaryVoltage: 100, secondaryCurrent: 0.1 } }), 0, 1e-3);
    expect(st.spark.kind).toBe('arc');
    const r1 = st.spark.radius;
    st.update(snap({ t: 0, spark: { phase: 'arc', secondaryVoltage: 100, secondaryCurrent: 0.4 } }), 0, 1e-3);
    expect(st.spark.radius / r1).toBeCloseTo(2, 6);
    expect(st.spark.drawRadius).toBeGreaterThanOrEqual(VIS_MIN_CHANNEL_RADIUS);
    const arcCore = [...st.spark.coreRadiance];
    st.update(snap({ t: 0, spark: { phase: 'glow', secondaryVoltage: 450, secondaryCurrent: 0.04 } }), 0, 1e-3);
    expect(st.spark.kind).toBe('glow');
    expect(st.spark.visible).toBe(true);
    const g = st.spark.coreRadiance;
    expect(g[2] / Math.max(g[1], 1e-12)).toBeGreaterThan(5); // violet
    expect(luminance(g)).toBeLessThan(luminance(arcCore));
    // colour follows the N2 band spectrum
    const gw = glowRGBPerWatt();
    expect(g[2] / g[0]).toBeCloseTo(gw[2] / gw[0], 6);
    // glow column is wider than the arc channel and thickens with current
    const rg = st.spark.radius;
    expect(rg).toBeGreaterThan(r1);
    st.update(snap({ t: 0, spark: { phase: 'glow', secondaryVoltage: 450, secondaryCurrent: 0.08 } }), 0, 1e-3);
    expect(st.spark.radius).toBeGreaterThan(rg);
  });

  it('frame-averaged power from the energy counter (camera exposure) vs instantaneous when paused', () => {
    const st = fresh();
    st.update(snap({ t: 0, spark: { phase: 'charging', energyDelivered: 0 } }), 0.016, 1);
    // real time: a 1 ms, 30 mJ discharge fully inside a 16 ms frame averages to ~1.9 W
    st.update(snap({ t: 0.016, spark: { phase: 'done', energyDelivered: 0.03 } }), 0.016, 1);
    expect(st.spark.power).toBeCloseTo(0.03 / 0.016, 6);
    st.update(snap({ t: 0.016, spark: { phase: 'glow', secondaryVoltage: 400, secondaryCurrent: 0.05, energyDelivered: 0.03 } }), 0, 1);
    expect(st.spark.power).toBeCloseTo(20, 9);
  });

  it('channel is convected by swirl and restrikes when over-stretched', () => {
    const st = fresh();
    st.flow.H = 0; // start quiescent; give it swirl via the flow model directly
    st.update(snap({ t: 0, spark: { phase: 'charging' } }), 0.016, 1e-3);
    st.update(snap({ t: 1e-6, spark: { phase: 'glow', secondaryVoltage: 400, secondaryCurrent: 0.05 } }), 0.016, 1e-3);
    // impose solid-body swirl of 300 rad/s (tangential ≈ 11 m/s at the gap)
    const I = 0.5 * 6e-4 * TEST_R * TEST_R;
    let maxBow = 0;
    for (let i = 0; i < 400; i++) {
      st.flow.H = 300 * I;
      st.update(snap({ t: 2e-6 + i * 1e-6, flame: { turbulenceIntensity: 1e-6 }, spark: { phase: 'glow', secondaryVoltage: 400, secondaryCurrent: 0.05 } }), 0.016, 1e-3);
      maxBow = Math.max(maxBow, Math.hypot(...st.spark.bow));
    }
    expect(maxBow).toBeGreaterThan(0.2e-3);
    expect(maxBow).toBeLessThanOrEqual(1.5 * 0.5e-3 + 1e-4);
  });
});

describe('time handling', () => {
  it('detects seeks backwards and resets transient state', () => {
    const st = fresh();
    st.update(flameSnap({ t: 0.1 }), 0.016, 1e-3);
    st.update(flameSnap({ t: 0.05 }), 0.016, 1e-3);
    expect(st.wasReset).toBe(true);
    expect(st.dtSim).toBe(0);
    st.update(flameSnap({ t: 0.0501 }), 0.016, 1e-3);
    expect(st.wasReset).toBe(false);
    expect(st.dtSim).toBeCloseTo(1e-4, 12);
  });

  it('chamber light: VISUAL_GAIN × luminous intensity of the flame', () => {
    const st = fresh();
    st.update(flameSnap({ temperatureBurned: 300 }), 0.016, 1e-3); // cold "burned" gas → only the front emits
    const I = (luminance(flameRGBPerWatt(st.phi)) * CHEMILUMINESCENCE_EFFICIENCY * st.hrrFront) / (4 * Math.PI);
    const c = st.light.color;
    expect(Math.max(...c)).toBeCloseTo(1, 12);
    // intensity × colour / luminance(colour) recovers I · gain
    expect((st.light.intensity * luminance(c))).toBeCloseTo(VISUAL_GAIN * I, 9);
  });
});

describe('L-head chamber (Model T)', () => {
  const T = chamberShapeOf(MODEL_T);
  const mt = () => new CombustionVisualState(MODEL_T);
  const burning = {
    phase: 'combustion' as const, temperatureUnburned: 650, temperatureBurned: 300, massFractionBurned: 0.05, heatReleaseRate: 3e4,
    flame: { stage: 'turbulent' as const, radius: 0.02, area: 1.5e-3, laminarSpeed: 0.4, turbulenceIntensity: 1 },
  };

  it('crown depth comes from pistonDisplacement + the spec (5/16 in above the deck at TDC), not clearanceHeight', () => {
    const st = mt();
    st.update(modelTSnap({ clearanceHeight: 0.5 }, 0), 0.016, 1e-3);
    expect(st.h).toBeCloseTo(0.0254, 12);
    expect(-st.h - MODEL_T.geometry.lHead!.deckY).toBeCloseTo((5 / 16) * 0.0254, 12);
    st.update(modelTSnap({ t: 1e-3, clearanceHeight: 0.5 }, 0.03), 0.016, 1e-3);
    expect(st.h).toBeCloseTo(0.0254 + 0.03, 12);
  });

  it('chamber light stays inside the chamber: over the pocket for the plug over the valves', () => {
    const st = mt();
    st.update(modelTSnap(burning, 0.002), 0.016, 1e-3);
    const p = st.light.position;
    expect(st.light.intensity).toBeGreaterThan(0);
    expect(insideChamberShape(p[0], p[1], p[2], T, st.h)).toBe(true);
    expect(Math.hypot(p[0], p[2])).toBeGreaterThan(T.R);
    // all burned: near the chamber's volume centroid, inside
    st.update(modelTSnap({ t: 1e-3, phase: 'expansion', temperatureBurned: 2200, massFractionBurned: 1, flame: { stage: 'done' } }, 0.03), 0.016, 1e-3);
    expect(st.allBurned).toBe(true);
    const q = st.light.position;
    expect(st.light.intensity).toBeGreaterThan(0);
    expect(insideChamberShape(q[0], q[1], q[2], T, st.h)).toBe(true);
  });

  it('knock: end gas at the far side of the bore; footprint axial mode f = c/(2L) (lower than the bore’s Draper mode)', () => {
    const st = mt();
    st.update(modelTSnap({ t: 0, ...burning }, 0.001), 0.016, 1e-4);
    st.update(modelTSnap({ t: 1e-5, ...burning, knock: { autoignited: true, oscillation: 5e4 } }, 0.001), 0.016, 1e-4);
    const k = st.knock;
    expect(k.active).toBe(true);
    expect(k.origin[0]).toBeCloseTo(T.R, 9);
    expect(k.origin[1]).toBeCloseTo(0, 9);
    const L = T.s1 - T.s0;
    expect(k.frequency).toBeCloseTo(k.soundSpeed / (2 * L), 9);
    expect(k.frequency).toBeLessThan(knockModeFrequency(k.soundSpeed, MODEL_T.geometry.bore));
    expect(k.modeSign).toBe(1);
    expect(axialModeShape(k.origin[0], k.origin[1], T.axis[0], T.axis[1], T.s0, L)).toBeCloseTo(1, 9);
    expect(k.ringRadius).toBe(0);
    expect(k.ringWidth).toBeCloseTo(0.12 * 0.5 * L, 12);
  });

  it('rigid-duct limit: axial mode c/(2L) with antinodes of opposite sign at the ends', () => {
    expect(axialModeFrequency(340, 0.5)).toBe(340);
    expect(axialModeShape(0.1, 0, 1, 0, 0.1, 0.2)).toBe(1);
    expect(axialModeShape(0.3, 5, 1, 0, 0.1, 0.2)).toBeCloseTo(-1, 12);
    expect(axialModeShape(0.2, -5, 1, 0, 0.1, 0.2)).toBeCloseTo(0, 12);
    // a uniform 1-D Neumann duct: ψ'' = −(π/L)² ψ (finite differences)
    const L = 0.2, d = 1e-4, s = 0.137;
    const f = (u: number) => axialModeShape(u, 0, 1, 0, 0, L);
    expect((f(s + d) - 2 * f(s) + f(s - d)) / (d * d)).toBeCloseTo(-((Math.PI / L) ** 2) * f(s), 3);
  });
});

describe('trembler spark showers (cumulative breakdownCount)', () => {
  const mt = (spec = MODEL_T) => new CombustionVisualState(spec);
  const trem = (t: number, n: number, phase: 'charging' | 'glow' | 'arc' = 'glow') =>
    modelTSnap({ t, spark: { phase, breakdownCount: n, breakdownVoltage: 8e3, secondaryVoltage: 300, secondaryCurrent: 0.03 } }, 0.001);

  it('every breakdown flashes exactly once at any playback speed; a new ignition event restarts the count', () => {
    const rt = mt();
    rt.update(trem(0, 0, 'charging'), 1 / 60, 1);
    rt.update(trem(1 / 60, 5), 1 / 60, 1); // five buzz sparks inside one real-time frame
    expect(rt.breakdownsSeen).toBe(5);
    const slow = mt();
    slow.update(trem(0, 0, 'charging'), 1 / 60, 1e-4);
    for (let i = 1; i <= 5; i++) slow.update(trem(i * 2e-6, i, i % 2 ? 'arc' : 'glow'), 1 / 60, 1e-4);
    expect(slow.breakdownsSeen).toBe(5);
    slow.update(trem(12e-6, 5), 1 / 60, 1e-4); // no new breakdown → no flash
    expect(slow.breakdownsSeen).toBe(5);
    rt.update(trem(2 / 60, 16), 1 / 60, 1);
    expect(rt.breakdownsSeen).toBe(16);
    rt.update(trem(3 / 60, 3), 1 / 60, 1); // next cylinder event: count restarted at 0
    expect(rt.breakdownsSeen).toBe(19);
  });

  it('flash energy = n · ½ C_sec V_bd² with the trembler coil’s secondary capacitance', () => {
    const run = (n: number, spec = MODEL_T) => {
      const st = mt(spec);
      st.update(trem(0, 0, 'charging'), 0.016, 1e-3);
      st.update(trem(1e-6, n), 0.016, 1e-3);
      return st.spark.flash[2];
    };
    const one = run(1);
    expect(one).toBeGreaterThan(0);
    expect(run(5) / one).toBeCloseTo(5, 9);
    const ig = MODEL_T.ignition;
    const doubled = { ...MODEL_T, ignition: { ...ig, coil: { ...ig.coil, secondaryCapacitance: 2 * ig.coil.secondaryCapacitance } } };
    expect(run(1, doubled) / one).toBeCloseTo(2, 9);
  });

  it('the first sample after a reset only latches the counter (a seek does not replay a shower)', () => {
    const st = mt();
    st.update(trem(0.5, 7), 0.016, 1e-3);
    expect(st.breakdownsSeen).toBe(0);
    expect(st.spark.flash[2]).toBe(0);
    st.update(trem(0.5 + 1e-6, 8), 0.016, 1e-3);
    expect(st.breakdownsSeen).toBe(1);
    st.update(trem(0.1, 9), 0.016, 1e-3); // seek backwards → reset → latch again
    expect(st.wasReset).toBe(true);
    expect(st.breakdownsSeen).toBe(1);
  });
});
