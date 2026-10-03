import { describe, expect, it } from 'vitest';
import {
  blackbodyXYZ,
  blackbodyXYZFast,
  chromaticity,
  cmfX,
  cmfY,
  cmfZ,
  planckRadiance,
  spectrumXYZPerWatt,
  xyzToLinearSrgb,
} from './cie';
import { inferno, infernoLUT, srgbToLinear, temperatureLegend, temperatureToUnit } from './colormap';
import {
  arcRGBPerWatt,
  equivalenceRatioFromBurned,
  flameRGBPerWatt,
  glowRGBPerWatt,
  luminance,
} from './emitters';

describe('CIE colour matching (Wyman et al. 2013 fit)', () => {
  it('matches tabulated CIE 1931 values at key wavelengths', () => {
    // CIE 1931 2° table: ȳ(555)=1.0002, x̄(600)=1.0622, z̄(450)=1.7471, x̄(442)≈0.3483
    expect(cmfY(555)).toBeCloseTo(1.0, 1);
    expect(Math.abs(cmfX(600) - 1.0622)).toBeLessThan(0.02);
    expect(Math.abs(cmfZ(450) - 1.7471)).toBeLessThan(0.05);
    expect(Math.abs(cmfX(440) - 0.3483)).toBeLessThan(0.03);
    expect(cmfY(300)).toBeLessThan(1e-3);
    expect(cmfY(800)).toBeLessThan(1e-3);
  });

  it('D65 white maps to sRGB (1,1,1)', () => {
    const rgb = xyzToLinearSrgb([0.95047, 1.0, 1.08883]);
    for (const c of rgb) expect(c).toBeCloseTo(1, 2);
  });
});

describe('black body', () => {
  it("Planck peak obeys Wien's law", () => {
    const T = 5000;
    let best = 0, lBest = 0;
    for (let l = 300; l <= 1200; l += 1) {
      const b = planckRadiance(l * 1e-9, T);
      if (b > best) { best = b; lBest = l; }
    }
    expect(Math.abs(lBest * 1e-9 * T - 2.8978e-3) / 2.8978e-3).toBeLessThan(0.005);
  });

  it('chromaticity of illuminant A (2856 K) and the 6504 K Planckian', () => {
    const [xa, ya] = chromaticity(blackbodyXYZ(2856));
    expect(Math.abs(xa - 0.4476)).toBeLessThan(0.004);
    expect(Math.abs(ya - 0.4074)).toBeLessThan(0.004);
    const [x6, y6] = chromaticity(blackbodyXYZ(6504));
    expect(Math.abs(x6 - 0.3135)).toBeLessThan(0.004);
    expect(Math.abs(y6 - 0.3237)).toBeLessThan(0.004);
  });

  it('absolute luminance: black body at the platinum point (2042 K) ≈ 60 cd/cm² (1948 candela definition)', () => {
    const Y = blackbodyXYZ(2042)[1];
    expect(Math.abs(Y - 6.0e5) / 6.0e5).toBeLessThan(0.05);
  });

  it('fast table agrees with direct integration', () => {
    const a: [number, number, number] = [0, 0, 0];
    for (const T of [800, 1234, 1800, 2400, 3100, 4400]) {
      const e = blackbodyXYZ(T);
      blackbodyXYZFast(T, a);
      for (let k = 0; k < 3; k++) expect(Math.abs(a[k] - e[k]) / e[k]).toBeLessThan(5e-3);
    }
  });
});

describe('emitter spectra', () => {
  it('spectrum XYZ per watt: a 555 nm line gives K_m·ȳ ≈ 683 cd/m² per W/m²/sr', () => {
    const xyz = spectrumXYZPerWatt([{ kind: 'line', nm: 555, fwhm: 1, weight: 1 }]);
    expect(Math.abs(xyz[1] - 683) / 683).toBeLessThan(0.01);
    // all-UV emission is invisible
    expect(spectrumXYZPerWatt([{ kind: 'line', nm: 309, fwhm: 4, weight: 1 }])[1]).toBeLessThan(0.5);
  });

  it('premixed flame front is blue and gets greener (C2*) when rich', () => {
    const s = flameRGBPerWatt(1);
    expect(s[2]).toBeGreaterThan(s[1]);
    expect(s[1]).toBeGreaterThanOrEqual(s[0]);
    const lean = flameRGBPerWatt(0.7), rich = flameRGBPerWatt(1.4);
    expect(rich[1] / rich[2]).toBeGreaterThan(lean[1] / lean[2]);
  });

  it('glow discharge is violet and much dimmer per watt than the arc', () => {
    const g = glowRGBPerWatt();
    expect(g[2]).toBeGreaterThan(g[0]);
    expect(g[0]).toBeGreaterThanOrEqual(g[1]);
    const a = arcRGBPerWatt(6000);
    expect(luminance(a)).toBeGreaterThan(10 * luminance(g));
    // arc is near-white
    expect(Math.min(a[0], a[1], a[2]) / Math.max(a[0], a[1], a[2])).toBeGreaterThan(0.5);
  });

  it('equivalence ratio from burned composition by element balance', () => {
    // C8H18 + a(O2 + 3.76 N2): stoich a = 12.5
    const mk = (a: number) => {
      const O2 = Math.max(a - 12.5, 0);
      const n = { CO2: 8, H2O: 9, CO: 0, O2, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 3.76 * a };
      const tot = Object.values(n).reduce((x, y) => x + y, 0);
      return Object.fromEntries(Object.entries(n).map(([k, v]) => [k, v / tot])) as typeof n;
    };
    expect(equivalenceRatioFromBurned(mk(12.5))).toBeCloseTo(1, 10);
    expect(equivalenceRatioFromBurned(mk(12.5 / 0.8))).toBeCloseTo(0.8, 10);
    // rich: C8H18 + 10 O2 → 5 CO2 + 3 CO + 7 H2O + 2 H2 (O atoms 20; stoich needs 25) → φ = 12.5/10
    const rich = { CO2: 5, H2O: 7, CO: 3, O2: 0, H2: 2, OH: 0, H: 0, O: 0, NO: 0 };
    expect(equivalenceRatioFromBurned(rich)).toBeCloseTo(1.25, 10);
    expect(equivalenceRatioFromBurned({ CO2: 0, H2O: 0, CO: 0, O2: 0.21, H2: 0, OH: 0, H: 0, O: 0, NO: 0 })).toBeNaN();
  });
});

describe('inferno colormap', () => {
  it('endpoints and midpoints match matplotlib', () => {
    expect(inferno(0)).toEqual([0.001462, 0.000466, 0.013866]);
    expect(inferno(1)).toEqual([0.988362, 0.998364, 0.644924]);
    // matplotlib inferno(x) = _inferno_data[int(x·256)]
    const ref: [number, [number, number, number]][] = [
      [0.25, [0.3415, 0.062325, 0.429425]],
      [0.5, [0.735683, 0.215906, 0.330245]],
      [0.75, [0.977092, 0.55085, 0.03905]],
    ];
    for (const [x, c] of ref) {
      const v = inferno(x);
      for (let k = 0; k < 3; k++) expect(Math.abs(v[k] - c[k])).toBeLessThan(0.012);
    }
  });

  it('is monotonic in luminance (perceptually ordered) and clamps', () => {
    let prev = -1;
    for (let i = 0; i <= 200; i++) {
      const c = inferno(i / 200);
      const L = 0.2126 * srgbToLinear(c[0]) + 0.7152 * srgbToLinear(c[1]) + 0.0722 * srgbToLinear(c[2]);
      expect(L).toBeGreaterThan(prev);
      prev = L;
    }
    expect(inferno(-3)).toEqual(inferno(0));
    expect(inferno(7)).toEqual(inferno(1));
    expect(inferno(Number.NaN)).toEqual(inferno(0));
  });

  it('LUT and temperature mapping', () => {
    const lut = infernoLUT(256);
    expect(lut.length).toBe(1024);
    expect(lut[4 * 255]).toBe(Math.round(0.988362 * 255));
    expect(temperatureToUnit(250, [250, 3000])).toBe(0);
    expect(temperatureToUnit(3000, [250, 3000])).toBe(1);
    expect(temperatureToUnit(1625, [250, 3000])).toBeCloseTo(0.5, 12);
  });

  it('legend helper', () => {
    const L = temperatureLegend([250, 3000]);
    expect(L.unit).toBe('K');
    expect(L.stops[0].color).toBe('#000004');
    expect(L.stops[L.stops.length - 1].position).toBe(1);
    expect(L.cssGradient.startsWith('linear-gradient(to right, #000004 0.0%')).toBe(true);
    expect(L.ticks.map((t) => t.value)).toEqual([500, 1000, 1500, 2000, 2500, 3000]);
    for (const t of L.ticks) expect(t.position).toBeCloseTo((t.value - 250) / 2750, 12);
  });
});
