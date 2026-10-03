import { describe, expect, it } from 'vitest';
import { DEG } from '../core/constants';
import type { WallSpec } from '../core/engine-spec';
import { ANNAND_RADIATION } from './correlations';
import {
  flatChamberAreas,
  newChamberAreas,
  newWallHeatResult,
  wallHeatLoss,
  wallHeatLossTwoZone,
  woschniCoefficient,
  type WoschniInputs,
} from './woschni';

// CFR-like inputs (passed in; the module never hard-codes engine numbers)
const IN = 0.0254;
const BORE = 3.25 * IN;
const STROKE = 4.5 * IN;
const ROD = 0.254;
const walls: WallSpec = {
  headTemperature: 525.2,
  pistonTemperature: 485.2,
  linerTemperature: 430.1,
  intakeValveTemperature: 464.9,
  exhaustValveTemperature: 519.6,
};
const valve = { count: 1, headDiameter: 1.471 * IN };

describe('wall heat-loss split', () => {
  const areas = flatChamberAreas(BORE, 0.02, valve, valve, newChamberAreas());
  const out = newWallHeatResult();

  it('flat chamber areas add up', () => {
    const disc = 0.25 * Math.PI * BORE * BORE;
    expect(areas.head + areas.intakeValves + areas.exhaustValves).toBeCloseTo(disc, 15);
    expect(areas.piston).toBeCloseTo(disc, 15);
    expect(areas.liner).toBeCloseTo(Math.PI * BORE * 0.02, 15);
    expect(flatChamberAreas(BORE, -1, valve, valve, newChamberAreas()).liner).toBe(0);
  });

  it('zero flux to a surface whose temperature equals the gas temperature (with radiation)', () => {
    const c = ANNAND_RADIATION.sparkIgnition;
    for (const [key, Tw] of [
      ['head', walls.headTemperature],
      ['piston', walls.pistonTemperature],
      ['liner', walls.linerTemperature],
      ['intakeValves', walls.intakeValveTemperature],
      ['exhaustValves', walls.exhaustValveTemperature],
    ] as const) {
      wallHeatLoss(800, Tw, areas, walls, out, c);
      expect(out[key]).toBe(0);
      wallHeatLossTwoZone(800, Tw, Tw, areas, { wettedHead: 1e-3, wettedPiston: 2e-3, wettedLiner: 1e-4 }, walls, out, c);
      expect(out[key]).toBe(0);
    }
    const uniform: WallSpec = { headTemperature: 450, pistonTemperature: 450, linerTemperature: 450, intakeValveTemperature: 450, exhaustValveTemperature: 450 };
    expect(wallHeatLoss(1000, 450, areas, uniform, out, c)).toBe(0);
    expect(wallHeatLoss(1000, 400, areas, uniform, out)).toBeLessThan(0); // wall heats the gas
  });

  it('single-zone sums; two-zone reduces to single-zone and splits by wetted area', () => {
    const h = 900;
    const T = 1500;
    const q1 = wallHeatLoss(h, T, areas, walls, out);
    const sum = out.head + out.piston + out.liner + out.intakeValves + out.exhaustValves;
    expect(q1).toBeCloseTo(sum, 9);
    expect(out.unburned).toBe(q1);
    const q2 = wallHeatLossTwoZone(h, T, T, areas, { wettedHead: 2e-3, wettedPiston: 1e-3, wettedLiner: 3e-4 }, walls, out);
    expect(q2).toBeCloseTo(q1, 9);
    // nothing wetted → all unburned; everything wetted → all burned
    wallHeatLossTwoZone(h, 700, 2400, areas, { wettedHead: 0, wettedPiston: 0, wettedLiner: 0 }, walls, out);
    expect(out.burned).toBe(0);
    const qu = wallHeatLoss(h, 700, areas, walls, newWallHeatResult());
    expect(out.total).toBeCloseTo(qu, 9);
    wallHeatLossTwoZone(h, 700, 2400, areas, { wettedHead: 1, wettedPiston: 1, wettedLiner: 1 }, walls, out);
    expect(out.unburned).toBe(0);
    expect(out.total).toBeCloseTo(wallHeatLoss(h, 2400, areas, walls, newWallHeatResult()), 9);
    // half the piston wetted: piston loss is the area-weighted mix
    wallHeatLossTwoZone(h, 700, 2400, areas, { wettedHead: 0, wettedPiston: 0.5 * areas.piston, wettedLiner: 0 }, walls, out);
    expect(out.piston).toBeCloseTo(0.5 * h * areas.piston * (2400 - walls.pistonTemperature) + 0.5 * h * areas.piston * (700 - walls.pistonTemperature), 9);
  });
});

/**
 * Magnitude sanity check on a synthetic CFR-like fired closed cycle (NOT a validation — the
 * system-level check of heat loss vs. fuel energy belongs to the cycle module): constant-γ
 * single zone, Wiebe burn, Woschni (1967). Expected from the literature: peak flux of order
 * 1–3 MW/m², 10–20 % (up to ~30 % at 600 rpm) of the fuel energy lost IVC → EVO.
 */
describe('CFR-like fired cycle: Woschni magnitude sanity', () => {
  it('peak flux 1–3 MW/m² and closed-cycle loss 8–35 % of fuel energy', () => {
    const CR = 7;
    const Ap = 0.25 * Math.PI * BORE * BORE;
    const a = STROKE / 2;
    const Vd = Ap * STROKE;
    const Vc = Vd / (CR - 1);
    const V = (th: number) => Vc + Ap * (ROD + a - a * Math.cos(th) - Math.sqrt(ROD * ROD - a * a * Math.sin(th) ** 2));
    const rpm = 600;
    const cm = (2 * STROKE * rpm) / 60;
    const R = 287;
    const g = 1.33;
    const cv = R / (g - 1);
    const thIVC = -146 * DEG;
    const thEVO = 140 * DEG;
    const pIVC = 0.95e5;
    const TIVC = 340;
    const VIVC = V(thIVC);
    const m = (pIVC * VIVC) / (R * TIVC);
    const far = 1.1 / 15.1; // iso-octane stoichiometric A/F ≈ 15.1
    const mf = (m * 0.95 * far) / (1 + far); // 5 % residual
    const LHV = 44.3e6;
    const Qfuel = mf * LHV;
    const Qrel = 0.93 * Qfuel; // rich mixture: incomplete oxidation (CO, H2)
    const thS = -13 * DEG;
    const dth = 55 * DEG;
    const xb = (th: number) => (th <= thS ? 0 : 1 - Math.exp(-5 * ((th - thS) / dth) ** 3));
    const areas = newChamberAreas();
    const res = newWallHeatResult();
    const wi: WoschniInputs = {
      bore: BORE,
      pressure: 0,
      temperature: 0,
      meanPistonSpeed: cm,
      phase: 'compression',
      displacedVolume: Vd,
      refPressure: pIVC,
      refTemperature: TIVC,
      refVolume: VIVC,
    };
    const omega = (rpm * 2 * Math.PI) / 60;
    // dT/dθ, and wall heat rate (W) and mean flux (W/m²) via closure
    let qw = 0;
    let flux = 0;
    const dTdth = (th: number, T: number): number => {
      const v = V(th);
      const p = (m * R * T) / v;
      wi.pressure = p;
      wi.temperature = T;
      wi.phase = th < thS ? 'compression' : 'combustion';
      wi.motoredPressure = pIVC * Math.pow(VIVC / v, g);
      const h = woschniCoefficient(wi);
      flatChamberAreas(BORE, v / Ap, valve, valve, areas); // disc chamber height V/A_p
      qw = wallHeatLoss(h, T, areas, walls, res);
      flux = qw / (areas.head + areas.intakeValves + areas.exhaustValves + areas.piston + areas.liner);
      const e = 1e-6;
      const dV = (V(th + e) - V(th - e)) / (2 * e);
      const dQ = (Qrel * (xb(th + e) - xb(th - e))) / (2 * e);
      return (dQ - p * dV - qw / omega) / (m * cv);
    };
    let T = TIVC;
    let Qw = 0;
    let peak = 0;
    let pmax = 0;
    const step = 0.1 * DEG;
    for (let th = thIVC; th < thEVO - 1e-12; th += step) {
      const k1 = dTdth(th, T);
      const q1 = qw;
      const k2 = dTdth(th + step / 2, T + (step / 2) * k1);
      const q2 = qw;
      const k3 = dTdth(th + step / 2, T + (step / 2) * k2);
      const q3 = qw;
      const k4 = dTdth(th + step, T + step * k3);
      const q4 = qw;
      T += (step / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
      Qw += ((step / omega) * (q1 + 2 * q2 + 2 * q3 + q4)) / 6;
      peak = Math.max(peak, flux);
      pmax = Math.max(pmax, (m * R * T) / V(th + step));
    }
    const lossFraction = Qw / Qfuel;
    console.info(
      `[heat] synthetic CFR cycle: p_max = ${(pmax / 1e5).toFixed(1)} bar, peak mean flux = ${(peak / 1e6).toFixed(2)} MW/m², ` +
        `closed-cycle loss = ${(100 * lossFraction).toFixed(1)} % of fuel energy`,
    );
    expect(pmax).toBeGreaterThan(20e5);
    expect(pmax).toBeLessThan(60e5);
    expect(peak).toBeGreaterThan(1e6);
    expect(peak).toBeLessThan(3e6);
    expect(lossFraction).toBeGreaterThan(0.08);
    expect(lossFraction).toBeLessThan(0.35);
  });
});
