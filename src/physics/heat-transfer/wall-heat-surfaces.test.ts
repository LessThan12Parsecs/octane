import { describe, expect, it } from 'vitest';
import type { WallSpec } from '../core/engine-spec';
import { ANNAND_RADIATION } from './correlations';
import {
  N_WALL_SURFACES,
  WALL_BLOCK,
  WALL_EXHAUST_VALVE,
  WALL_HEAD,
  WALL_INTAKE_VALVE,
  WALL_LINER,
  WALL_PISTON,
  WALL_SURFACE_IDS,
  chamberAreasToSurfaces,
  flatChamberAreas,
  newChamberAreas,
  newWallHeatResult,
  newWallSurfaceArray,
  newWallSurfaceHeatResult,
  wallHeatLoss,
  wallHeatLossSurfaces,
  wallHeatLossTwoZone,
  wallHeatLossTwoZoneSurfaces,
  wallSurfaceResistances,
  wallSurfaceTemperatures,
} from './wall-heat';

const IN = 0.0254;
const BORE = 3.25 * IN;
const walls: WallSpec = {
  headTemperature: 525.2,
  pistonTemperature: 485.2,
  linerTemperature: 430.1,
  intakeValveTemperature: 464.9,
  exhaustValveTemperature: 519.6,
};
const valve = { count: 1, headDiameter: 1.471 * IN };

describe('generic 6-surface wall heat (chamber-generic path)', () => {
  it('surface ids follow the cycle-model ledger order with the block last; block falls back to the liner', () => {
    expect(WALL_SURFACE_IDS).toEqual(['head', 'piston', 'liner', 'intakeValve', 'exhaustValve', 'block']);
    expect([WALL_HEAD, WALL_PISTON, WALL_LINER, WALL_INTAKE_VALVE, WALL_EXHAUST_VALVE, WALL_BLOCK]).toEqual([0, 1, 2, 3, 4, 5]);
    expect(N_WALL_SURFACES).toBe(6);
    const T = wallSurfaceTemperatures(walls, newWallSurfaceArray());
    expect(Array.from(T)).toEqual([525.2, 485.2, 430.1, 464.9, 519.6, 430.1]);
    expect(wallSurfaceTemperatures({ ...walls, blockTemperature: 450 }, newWallSurfaceArray())[WALL_BLOCK]).toBe(450);
    const Rr = wallSurfaceResistances({ head: 1, piston: 2, liner: 3, intakeValve: 4, exhaustValve: 5 }, newWallSurfaceArray());
    expect(Array.from(Rr)).toEqual([1, 2, 3, 4, 5, 3]);
    expect(wallSurfaceResistances({ head: 1, piston: 2, liner: 3, intakeValve: 4, exhaustValve: 5, block: 7 }, newWallSurfaceArray())[WALL_BLOCK]).toBe(7);
  });

  it('flat-disc areas and head-disc fractions reproduce wallHeatLoss / wallHeatLossTwoZone bit for bit', () => {
    const Tw = wallSurfaceTemperatures(walls, newWallSurfaceArray());
    const areas6 = newWallSurfaceArray();
    const f = newWallSurfaceArray();
    const q5 = newWallHeatResult();
    const q6 = newWallSurfaceHeatResult();
    let seed = 13;
    const rnd = () => (seed = (seed * 16807) % 2147483647) / 2147483647;
    for (let i = 0; i < 2000; i++) {
      const a = flatChamberAreas(BORE, 0.001 + 0.1 * rnd(), valve, valve, newChamberAreas());
      chamberAreasToSurfaces(a, areas6);
      const wet = { wettedHead: rnd() * 0.006, wettedPiston: rnd() * 0.006, wettedLiner: rnd() * 0.02 };
      const hc = 100 + 2000 * rnd();
      const Tu = 300 + 1000 * rnd();
      const Tb = 1500 + 1200 * rnd();
      const c = i % 2 ? ANNAND_RADIATION.sparkIgnition : 0;
      const headDisc = a.head + a.intakeValves + a.exhaustValves;
      const fh = Math.min(1, wet.wettedHead / headDisc);
      f[WALL_HEAD] = fh;
      f[WALL_INTAKE_VALVE] = fh;
      f[WALL_EXHAUST_VALVE] = fh;
      f[WALL_PISTON] = Math.min(1, wet.wettedPiston / a.piston);
      f[WALL_LINER] = Math.min(1, wet.wettedLiner / a.liner);
      f[WALL_BLOCK] = 0;
      wallHeatLossTwoZone(hc, Tu, Tb, a, wet, walls, q5, c);
      wallHeatLossTwoZoneSurfaces(hc, Tu, Tb, areas6, f, Tw, q6, c);
      expect(q6.total).toBe(q5.total);
      expect(q6.burned).toBe(q5.burned);
      expect(q6.unburned).toBe(q5.unburned);
      expect(q6.surface[WALL_HEAD]).toBe(q5.head);
      expect(q6.surface[WALL_INTAKE_VALVE]).toBe(q5.intakeValves);
      expect(q6.surface[WALL_EXHAUST_VALVE]).toBe(q5.exhaustValves);
      expect(q6.surface[WALL_PISTON]).toBe(q5.piston);
      expect(q6.surface[WALL_LINER]).toBe(q5.liner);
      expect(wallHeatLossSurfaces(hc, Tu, areas6, Tw, q6, c)).toBe(wallHeatLoss(hc, Tu, a, walls, q5, c));
    }
  });

  it('per-surface fractions: 0 → unburned only, 1 → burned only; the block is a surface of its own', () => {
    const Tw = wallSurfaceTemperatures({ ...walls, blockTemperature: 400 }, newWallSurfaceArray());
    const areas = Float64Array.from([0.02, 0.008, 0.003, 0.001, 0.001, 0.0025]);
    const out = newWallSurfaceHeatResult();
    const h = 800;
    const zero = newWallSurfaceArray();
    const one = newWallSurfaceArray().fill(1);
    wallHeatLossTwoZoneSurfaces(h, 700, 2400, areas, zero, Tw, out);
    expect(out.burned).toBe(0);
    expect(out.unburned).toBeCloseTo(wallHeatLossSurfaces(h, 700, areas, Tw, newWallSurfaceHeatResult()), 9);
    wallHeatLossTwoZoneSurfaces(h, 700, 2400, areas, one, Tw, out);
    expect(out.unburned).toBe(0);
    expect(out.burned).toBeCloseTo(wallHeatLossSurfaces(h, 2400, areas, Tw, newWallSurfaceHeatResult()), 9);
    // only the block burned: its share at T_b, everything else at T_u
    const fb = newWallSurfaceArray();
    fb[WALL_BLOCK] = 1;
    wallHeatLossTwoZoneSurfaces(h, 700, 2400, areas, fb, Tw, out);
    expect(out.burned).toBeCloseTo(h * areas[WALL_BLOCK] * (2400 - 400), 9);
    let sum = 0;
    for (let i = 0; i < N_WALL_SURFACES; i++) sum += out.surface[i];
    expect(sum).toBeCloseTo(out.total, 9);
    // fractions outside [0, 1] (or NaN) are clamped
    const bad = Float64Array.from([2, -1, Number.NaN, 0.5, 0.5, 0.5]);
    wallHeatLossTwoZoneSurfaces(h, 700, 2400, areas, bad, Tw, out);
    expect(Number.isFinite(out.total)).toBe(true);
    // zero flux to a surface at the gas temperature
    wallHeatLossTwoZoneSurfaces(h, 400, 400, areas, fb, Tw, out, ANNAND_RADIATION.sparkIgnition);
    expect(out.surface[WALL_BLOCK]).toBe(0);
  });
});
