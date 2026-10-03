import { describe, expect, it } from 'vitest';
import { flameRGBPerWatt, glowRGBPerWatt, luminance } from './colour/emitters';
import {
  CHEMILUMINESCENCE_EFFICIENCY,
  OVERLAY_KNOCK_FULL_SCALE_PA,
  VISUAL_GAIN,
  VIS_MIN_CHANNEL_RADIUS,
} from './constants';
import { CombustionVisualState } from './state';
import { flameSnap, snap, testSpec, TEST_R } from './test-utils';

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
