import { describe, expect, it } from 'vitest';
import { getPreset, LOCAL_MON_CONDITIONS, LOCAL_RON_CONDITIONS, monSparkAdvance } from './presets';
import { humidityRatio, kelvinToCelsius } from './units';

describe('CFR rating presets', () => {
  it('RON: 600 rpm, 13° BTDC, 52 °C intake, 100 °C coolant, WOT', () => {
    const p = LOCAL_RON_CONDITIONS;
    expect(p.rpm).toBe(600);
    expect(p.sparkAdvanceDeg).toBe(13);
    expect(kelvinToCelsius(p.intakeMixtureTemperature!)).toBeCloseTo(52, 9);
    expect(kelvinToCelsius(p.coolantTemperature!)).toBeCloseTo(100, 9);
    expect(p.throttle).toBe(1);
    expect(p.speedMode).toBe('fixed');
  });

  it('MON: 900 rpm, 149 °C mixture', () => {
    const p = LOCAL_MON_CONDITIONS;
    expect(p.rpm).toBe(900);
    expect(kelvinToCelsius(p.intakeMixtureTemperature!)).toBeCloseTo(149, 9);
  });

  it('humidity lands inside the ASTM 3.56–7.12 g/kg band', () => {
    const p = LOCAL_RON_CONDITIONS;
    const w = humidityRatio(p.relativeHumidity!, p.ambientTemperature!, p.ambientPressure!);
    expect(w).toBeGreaterThan(3.56e-3);
    expect(w).toBeLessThan(7.12e-3);
  });

  it('presets never set the compression ratio (it is the rating variable)', () => {
    expect(LOCAL_RON_CONDITIONS.compressionRatio).toBeUndefined();
    expect(LOCAL_MON_CONDITIONS.compressionRatio).toBeUndefined();
  });

  it('MON spark advance falls from 26° to 14° BTDC with CR and is clamped', () => {
    expect(monSparkAdvance(4)).toBe(26);
    expect(monSparkAdvance(5)).toBe(26);
    expect(monSparkAdvance(14)).toBe(14);
    expect(monSparkAdvance(18)).toBe(14);
    let prev = Infinity;
    for (let cr = 4; cr <= 18; cr += 0.5) {
      const a = monSparkAdvance(cr);
      expect(a).toBeLessThanOrEqual(prev);
      prev = a;
    }
  });

  it('getPreset: MON schedules spark with the current CR', () => {
    const mon = getPreset('MON', 8);
    expect(mon.patch.rpm).toBe(900);
    if (mon.sparkFollowsCR) expect(mon.patch.sparkAdvanceDeg).toBeCloseTo(monSparkAdvance(8), 12);
    const ron = getPreset('RON', 8);
    expect(ron.sparkFollowsCR).toBe(false);
    expect(ron.label).toContain('RON');
  });
});
