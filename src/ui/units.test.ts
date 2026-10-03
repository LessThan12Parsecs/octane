import { describe, expect, it } from 'vitest';
import {
  celsiusToKelvin,
  cyclePeriod,
  fmt,
  fmtMicros,
  fmtMicrosTick,
  fmtSig,
  fmtSigned,
  fmtSimTime,
  humidityRatio,
  isfcToGPerKWh,
  kelvinToCelsius,
  kPaToPa,
  paToBar,
  paToKPa,
  relativeHumidityForHumidityRatio,
  saturationVapourPressure,
  secondsPerDeg,
  wattsToJPerDeg,
} from './units';

describe('unit conversions', () => {
  it('temperature and pressure round-trip', () => {
    expect(kelvinToCelsius(celsiusToKelvin(52))).toBeCloseTo(52, 12);
    expect(celsiusToKelvin(149)).toBeCloseTo(422.15, 12);
    expect(paToBar(101325)).toBeCloseTo(1.01325, 12);
    expect(kPaToPa(paToKPa(98765))).toBeCloseTo(98765, 9);
  });

  it('ISFC kg/J → g/kWh', () => {
    // 250 g/kWh = 0.25 kg / 3.6e6 J
    expect(isfcToGPerKWh(0.25 / 3.6e6)).toBeCloseTo(250, 9);
  });

  it('heat-release rate W → J/°CA and crank timing', () => {
    // 600 rpm → 3600 °/s: 36 kW is 10 J/°
    expect(wattsToJPerDeg(36e3, 600)).toBeCloseTo(10, 12);
    expect(wattsToJPerDeg(1, 0)).toBeNaN();
    expect(secondsPerDeg(600)).toBeCloseTo(1 / 3600, 15);
    expect(cyclePeriod(600)).toBeCloseTo(0.2, 15);
    expect(cyclePeriod(-1)).toBeNaN();
  });

  it('saturation vapour pressure matches steam tables within 0.2 %', () => {
    // IAPWS reference values: 25 °C → 3169.9 Pa, 100 °C → 101418 Pa
    expect(saturationVapourPressure(298.15) / 3169.9 - 1).toBeLessThan(2e-3);
    expect(Math.abs(saturationVapourPressure(373.15) / 101418 - 1)).toBeLessThan(5e-3);
  });

  it('humidity ratio ↔ relative humidity round-trip', () => {
    const T = 298.15;
    const p = 101325;
    const w = humidityRatio(0.4, T, p);
    expect(w).toBeGreaterThan(0.007);
    expect(w).toBeLessThan(0.009);
    expect(relativeHumidityForHumidityRatio(w, T, p)).toBeCloseTo(0.4, 12);
  });
});

describe('formatting', () => {
  it('fmt handles NaN, -0 and uses a typographic minus', () => {
    expect(fmt(NaN, 2)).toBe('—');
    expect(fmt(Infinity, 1)).toBe('—');
    expect(fmt(-0.0001, 2)).toBe('0.00');
    expect(fmt(-1.5, 1)).toBe('−1.5');
    expect(fmt(3.14159, 3)).toBe('3.142');
  });

  it('fmtSigned', () => {
    expect(fmtSigned(2, 1)).toBe('+2.0');
    expect(fmtSigned(-2, 1)).toBe('−2.0');
    expect(fmtSigned(0, 1)).toBe('0.0');
  });

  it('fmtSig keeps significant figures', () => {
    expect(fmtSig(0.0012345, 3)).toBe('0.00123');
    expect(fmtSig(1234.5, 3)).toBe('1235');
    expect(fmtSig(0, 3)).toBe('0');
  });

  it('fmtMicros switches between µs and ms', () => {
    expect(fmtMicros(250)).toBe('250 µs');
    expect(fmtMicros(1500)).toBe('1.50 ms');
    expect(fmtMicros(-200)).toBe('−200 µs');
    expect(fmtMicros(0)).toBe('0');
    expect(fmtMicros(2.5)).toBe('2.5 µs');
  });

  it('fmtMicrosTick', () => {
    expect(fmtMicrosTick(0)).toBe('0');
    expect(fmtMicrosTick(500)).toBe('500\u2009µs');
    expect(fmtMicrosTick(2500)).toBe('2.5\u2009ms');
    expect(fmtMicrosTick(-200)).toBe('\u2212200\u2009µs');
  });

  it('fmtSimTime', () => {
    expect(fmtSimTime(0.0123)).toBe('12.30 ms');
    expect(fmtSimTime(12.3456)).toBe('12.346 s');
  });
});
