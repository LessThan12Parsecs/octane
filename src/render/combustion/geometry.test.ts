import { describe, expect, it } from 'vitest';
import { ALPHA_10, INTEGRAL_SCALE_FRACTION } from './constants';
import {
  besselJ1,
  brushParams,
  flameSignedDistance,
  insideBurnedRegion,
  insideChamber,
  insideFlameSphere,
  knockModeFrequency,
  knockModeShape,
  molarMassOf,
  openArc,
  openCurtainArea,
  raySphere,
  rayChamberInterval,
  sampleOpenArc,
  soundSpeed,
  swirlArm,
} from './geometry';
import { TEST_R, testSpec } from './test-utils';

const R = TEST_R;

describe('inside tests (used by the tracers)', () => {
  it('chamber disc: radius, head face and piston crown are inclusive boundaries', () => {
    const h = 0.02;
    expect(insideChamber(0, -0.01, 0, R, h)).toBe(true);
    expect(insideChamber(R, 0, 0, R, h)).toBe(true);
    expect(insideChamber(0, -h, 0, R, h)).toBe(true);
    expect(insideChamber(R * 1.0001, -0.01, 0, R, h)).toBe(false);
    expect(insideChamber(0, 1e-6, 0, R, h)).toBe(false);
    expect(insideChamber(0, -h - 1e-6, 0, R, h)).toBe(false);
    expect(insideChamber(R * 0.8, -0.01, R * 0.7, R, h)).toBe(false); // outside the circle, inside the square
  });

  it('burned region = flame sphere ∩ chamber (sphere at a side plug is clipped by head and liner)', () => {
    const c: [number, number, number] = [0, -0.002, R - 0.004];
    const r = 0.015, h = 0.02;
    expect(insideFlameSphere(0, -0.002, R - 0.004 - 0.014, c, r)).toBe(true);
    expect(insideBurnedRegion(0, -0.002, R - 0.004 - 0.014, R, h, c, r)).toBe(true);
    // inside the sphere but above the head face → not burned gas
    expect(insideFlameSphere(0, 0.005, R - 0.004, c, r)).toBe(true);
    expect(insideBurnedRegion(0, 0.005, R - 0.004, R, h, c, r)).toBe(false);
    // inside the sphere but outside the liner
    expect(insideBurnedRegion(0, -0.002, R + 0.005, R, h, c, r)).toBe(false);
    expect(insideFlameSphere(0, 0, 0, c, 0)).toBe(false);
    expect(flameSignedDistance(0, -0.002, R - 0.004 - 0.02, c, r)).toBeCloseTo(0.005, 12);
  });

  it('volume fraction of a centred sphere in the disc matches the analytic sphere-cap volume', () => {
    // Monte Carlo over the chamber box: sphere at head centre, r < R, cut by the head plane (hemisphere).
    const h = 0.02, r = 0.012;
    const c: [number, number, number] = [0, 0, 0];
    let seed = 7;
    const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    let inC = 0, inB = 0;
    for (let i = 0; i < 200000; i++) {
      const x = (2 * rnd() - 1) * R, z = (2 * rnd() - 1) * R, y = -rnd() * h;
      if (!insideChamber(x, y, z, R, h)) continue;
      inC++;
      if (insideBurnedRegion(x, y, z, R, h, c, r)) inB++;
    }
    const frac = inB / inC;
    const exact = ((2 / 3) * Math.PI * r ** 3) / (Math.PI * R * R * h);
    expect(Math.abs(frac - exact) / exact).toBeLessThan(0.03);
  });
});

describe('ray intersections (mirrored in the GLSL)', () => {
  it('chamber interval agrees with brute-force sampling along random rays', () => {
    const h = 0.018;
    let seed = 3;
    const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    for (let n = 0; n < 300; n++) {
      const ro = [(rnd() - 0.5) * 0.3, (rnd() - 0.5) * 0.3, (rnd() - 0.5) * 0.3];
      let rd = [rnd() - 0.5, rnd() - 0.5, rnd() - 0.5];
      const L = Math.hypot(rd[0], rd[1], rd[2]);
      rd = rd.map((v) => v / L);
      const iv = rayChamberInterval(ro, rd, R, h);
      let first = -1, last = -1;
      for (let t = 0; t < 0.6; t += 1e-4) {
        if (insideChamber(ro[0] + rd[0] * t, ro[1] + rd[1] * t, ro[2] + rd[2] * t, R, h)) {
          if (first < 0) first = t;
          last = t;
        }
      }
      if (first < 0) {
        expect(iv === null || iv[1] - iv[0] < 3e-4).toBe(true);
      } else {
        expect(iv).not.toBeNull();
        expect(Math.abs(iv![0] - first)).toBeLessThan(2e-4);
        expect(Math.abs(iv![1] - last)).toBeLessThan(2e-4);
      }
    }
  });

  it('camera inside the chamber starts at t = 0; vertical rays work', () => {
    const iv = rayChamberInterval([0, -0.01, 0], [0, -1, 0], R, 0.02);
    expect(iv![0]).toBe(0);
    expect(iv![1]).toBeCloseTo(0.01, 12);
    expect(rayChamberInterval([R * 1.1, 0.1, 0], [0, -1, 0], R, 0.02)).toBeNull();
  });

  it('sphere entry/exit', () => {
    const s = raySphere([0, 0, -1], [0, 0, 1], [0, 0, 0], 0.5)!;
    expect(s[0]).toBeCloseTo(0.5, 12);
    expect(s[1]).toBeCloseTo(1.5, 12);
    expect(raySphere([0, 1, -1], [0, 0, 1], [0, 0, 0], 0.5)).toBeNull();
  });
});

describe('turbulent flame brush', () => {
  const base = { uPrime: 2, laminarSpeed: 0.5, clearanceHeight: 0.02, radius: 0.03, pressure: 20e5, unburnedTemperature: 750 };

  it('grows with u′/S_L and is bounded by 2 L_I and by the radius', () => {
    let prev = 0;
    for (const up of [0, 0.5, 1, 2, 4, 8]) {
      const b = brushParams({ ...base, uPrime: up });
      expect(b.thickness).toBeGreaterThanOrEqual(prev);
      expect(b.thickness).toBeLessThanOrEqual(2 * INTEGRAL_SCALE_FRACTION * base.clearanceHeight + 1e-12);
      prev = b.thickness;
    }
    const lowSL = brushParams({ ...base, laminarSpeed: 0.25 });
    const hiSL = brushParams({ ...base, laminarSpeed: 1.0 });
    expect(lowSL.thickness).toBeGreaterThan(hiSL.thickness);
    const kernel = brushParams({ ...base, radius: 0.5e-3 });
    expect(kernel.thickness).toBeLessThanOrEqual(0.5e-3 + 1e-12);
    expect(kernel.wrinkleAmplitude).toBeLessThanOrEqual(0.25e-3 + 1e-12);
  });

  it('laminar limit: no turbulence → brush = laminar thickness, no wrinkles', () => {
    const b = brushParams({ ...base, uPrime: 0 });
    expect(b.thickness).toBeCloseTo(b.laminarThickness, 15);
    expect(b.wrinkleAmplitude).toBe(0);
    // δ_L ~ 10–100 µm at engine conditions and thinner at higher pressure
    expect(b.laminarThickness).toBeGreaterThan(5e-6);
    expect(b.laminarThickness).toBeLessThan(1e-4);
    expect(brushParams({ ...base, uPrime: 0, pressure: 40e5 }).laminarThickness).toBeLessThan(b.laminarThickness);
  });

  it('Blizard–Keck burn-up zone δ_b = λ(1 + u′/S_L) with λ = L sqrt(15/Re)', () => {
    const b = brushParams(base);
    const L = INTEGRAL_SCALE_FRACTION * base.clearanceHeight;
    const nu = 1.6e-5 * Math.pow(750 / 300, 1.7) * (1e5 / 20e5);
    const lambda = L * Math.sqrt(15 / ((2 * L) / nu));
    expect(b.taylorScale).toBeCloseTo(lambda, 12);
    const dev = 1 - Math.exp(-0.03 / L);
    expect(b.thickness).toBeCloseTo(b.laminarThickness + lambda * (1 + 2 / 0.5) * dev, 12);
    expect(b.wrinkleWavelength).toBeCloseTo(L, 15);
    // plausible engine numbers: λ ~ 0.1–1 mm, δ ~ 1–3 mm
    expect(lambda).toBeGreaterThan(1e-4);
    expect(lambda).toBeLessThan(1e-3);
  });

  it('kernel development factor rises from 0 to 1 with r/L_I', () => {
    expect(brushParams({ ...base, radius: 0 }).development).toBe(0);
    expect(brushParams({ ...base, radius: 1 }).development).toBeCloseTo(1, 10);
  });
});

describe('knock acoustics', () => {
  it('Bessel J1', () => {
    expect(besselJ1(1)).toBeCloseTo(0.4400505857, 8);
    expect(besselJ1(ALPHA_10)).toBeCloseTo(0.5818652, 6);
    expect(besselJ1(0)).toBe(0);
    // α_10 is a zero of J1′ (rigid-wall boundary condition)
    const d = (besselJ1(ALPHA_10 + 1e-5) - besselJ1(ALPHA_10 - 1e-5)) / 2e-5;
    expect(Math.abs(d)).toBeLessThan(1e-4);
  });

  it('Draper (1,0) frequency ≈ 6–7 kHz for the CFR bore with hot gas', () => {
    const c = soundSpeed(2200, 0.0285, 1.3);
    const f = knockModeFrequency(c, 2 * R);
    expect(f).toBeCloseTo((1.8412 * c) / (Math.PI * 2 * R), 6);
    expect(f).toBeGreaterThan(5500);
    expect(f).toBeLessThan(7500);
  });

  it('mode shape: unit antinode at the wall on its axis, node on the axis of the cylinder', () => {
    expect(knockModeShape(R, 0, R, 0)).toBeCloseTo(1, 6);
    expect(knockModeShape(-R, 0, R, 0)).toBeCloseTo(-1, 6);
    expect(knockModeShape(0, R, R, 0)).toBeCloseTo(0, 6);
    expect(knockModeShape(0, 0, R, 0)).toBe(0);
  });

  it('molar mass from composition', () => {
    const air = { CO2: 0, H2O: 0, CO: 0, O2: 0.21, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0.79 };
    expect(molarMassOf(air)).toBeCloseTo(0.21 * 31.9988e-3 + 0.79 * 28.0134e-3, 9);
  });
});

describe('valves and swirl geometry', () => {
  const spec = testSpec();

  it('open arc and curtain area of a 180° shrouded valve', () => {
    const a = openArc(spec.intakeValve);
    expect(a.halfWidth).toBeCloseTo(Math.PI / 2, 12);
    const L = 5e-3;
    expect(openCurtainArea(spec.intakeValve, L)).toBeCloseTo(0.5 * Math.PI * spec.intakeValve.seatInnerDiameter * L, 12);
    expect(openCurtainArea(spec.exhaustValve, L)).toBeCloseTo(Math.PI * spec.exhaustValve.seatInnerDiameter * L, 12);
    for (let i = 0; i <= 10; i++) {
      const psi = sampleOpenArc(spec.intakeValve, i / 10);
      expect(Math.abs(psi - spec.intakeValve.shroudDirection)).toBeLessThanOrEqual(Math.PI / 2 + 1e-12);
    }
  });

  it('swirl moment arm: sign and magnitude follow the shroud orientation', () => {
    // valve at x = −0.02 facing +z → counter-clockwise (+y) swirl, arm = 0.02·sinc(π/2) = 0.02·2/π
    expect(swirlArm(spec.intakeValve)).toBeCloseTo(0.02 * (2 / Math.PI), 12);
    expect(swirlArm({ ...spec.intakeValve, shroudDirection: -Math.PI / 2 })).toBeCloseTo(-0.02 * (2 / Math.PI), 12);
    // facing the cylinder axis (+x): jets point through the centre → no swirl
    expect(swirlArm({ ...spec.intakeValve, shroudDirection: 0 })).toBeCloseTo(0, 12);
    // unshrouded: symmetric → none
    expect(swirlArm(spec.exhaustValve)).toBe(0);
  });
});
