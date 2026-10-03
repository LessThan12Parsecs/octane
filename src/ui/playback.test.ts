import { describe, expect, it } from 'vitest';
import {
  formatTimeScale,
  MAX_TIME_SCALE,
  MIN_TIME_SCALE,
  sliderToTimeScale,
  snapTimeScale,
  stepTimeScale,
  timeScaleToSlider,
  wallSecondsPerDegrees,
} from './playback';

describe('time scale', () => {
  it('slider ends map to 1/5000 and real time', () => {
    expect(sliderToTimeScale(0)).toBeCloseTo(MIN_TIME_SCALE, 15);
    expect(sliderToTimeScale(1)).toBe(MAX_TIME_SCALE);
    expect(sliderToTimeScale(-3)).toBeCloseTo(MIN_TIME_SCALE, 15);
    expect(sliderToTimeScale(7)).toBe(1);
  });

  it('slider mapping is logarithmic and round-trips (up to snapping)', () => {
    const mid = sliderToTimeScale(0.5);
    expect(1 / mid).toBeCloseTo(71, 0); // sqrt(5000) ≈ 70.7 → snapped to 71
    for (const u of [0.1, 0.33, 0.8]) {
      const ts = sliderToTimeScale(u);
      expect(timeScaleToSlider(ts)).toBeCloseTo(u, 2);
    }
  });

  it('snaps slow-down factors to two significant figures', () => {
    expect(1 / snapTimeScale(1 / 1234)).toBeCloseTo(1200, 9);
    expect(1 / snapTimeScale(1 / 250)).toBeCloseTo(250, 9);
    expect(1 / snapTimeScale(1 / 1.37)).toBeCloseTo(1.4, 9);
    expect(snapTimeScale(10)).toBe(1);
    expect(snapTimeScale(1e-9)).toBeCloseTo(MIN_TIME_SCALE, 15);
  });

  it('formats as a slow-down fraction', () => {
    expect(formatTimeScale(1)).toBe('1× (real time)');
    expect(formatTimeScale(1 / 250)).toBe('1/250×');
    expect(formatTimeScale(1 / 5000)).toBe('1/5000×');
    expect(formatTimeScale(1 / 1.4)).toBe('1/1.4×');
  });

  it('steps along the 1-2-5 ladder', () => {
    expect(stepTimeScale(1 / 50, -1)).toBeCloseTo(1 / 100, 15);
    expect(stepTimeScale(1 / 50, 1)).toBeCloseTo(1 / 20, 15);
    expect(stepTimeScale(1 / 71, 1)).toBeCloseTo(1 / 50, 15); // off-ladder value
    expect(stepTimeScale(1 / 71, -1)).toBeCloseTo(1 / 100, 15);
    expect(stepTimeScale(1, 1)).toBe(1);
    expect(stepTimeScale(1 / 5000, -1)).toBeCloseTo(1 / 5000, 15);
  });

  it('wall time to play crank degrees', () => {
    // 600 rpm: 1° = 1/3600 s; at 1/5000 → 1.39 s per degree
    expect(wallSecondsPerDegrees(1, 600, 1 / 5000)).toBeCloseTo(5000 / 3600, 12);
    expect(wallSecondsPerDegrees(1, 0, 1)).toBeNaN();
  });
});
