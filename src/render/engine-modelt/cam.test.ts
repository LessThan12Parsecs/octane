import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/render_modelt_cam.json';
import { MODEL_T } from '../../physics/engines/model-t';
import { buildLobeProfile, lobeHeightAt, wrapDeg720 } from '../engine/cam-profile';
import { ValveCam, lobeCentreDeg, threeArcGeometry, threeArcLift } from './cam';
import { lobeProfileOf } from './parts-valvetrain';

const IN = 0.0254;
const DEG = Math.PI / 180;

/** Seat-to-seat duration (crank deg) of a symmetric lobe at a given lash, by bisection on the flank. */
function duration(g: ReturnType<typeof threeArcGeometry>, lash: number): number {
  let lo = 0, hi = g.f0;
  for (let k = 0; k < 80; k++) {
    const m = (lo + hi) / 2;
    if (threeArcLift(g, m) > lash) lo = m; else hi = m;
  }
  return (4 * lo) / DEG;
}

describe('three-arc flat-follower cam (Ford 1913–27 lobe)', () => {
  it('matches the brute-force support-function oracle (tools/reference/render_modelt_cam.py)', () => {
    for (const c of fixture.cases) {
      const g = threeArcGeometry(c);
      expect(g.f0 / DEG).toBeCloseTo(c.flankBaseAngleDeg, 6);
      c.crankDegFromCentre.forEach((th, k) => {
        expect(Math.abs(threeArcLift(g, (th / 2) * DEG) - c.grossLift[k])).toBeLessThan(2e-9);
      });
      for (const d of c.seatToSeat) expect(duration(g, d.lashIn * IN)).toBeCloseTo(d.durationDeg, 3);
    }
  });

  it('reproduces the MTFC Tulsa lift table and junction angles (design_stock.htm)', () => {
    // gross lift (in) at crank degrees from the lobe centreline, from the MTFC equations
    const table: [number, number][] = [
      [0, 0.2502], [10, 0.2478], [20, 0.2407], [30, 0.2289], [40, 0.2125], [50, 0.1916], [60, 0.1665],
      [70, 0.1372], [80, 0.104], [90, 0.0713], [100, 0.0445], [110, 0.0238], [120, 0.0095], [130, 0.0017], [137.2, 0],
    ];
    const g = threeArcGeometry({ baseRadius: 0.406 * IN, flankRadius: 1.2601 * IN, noseRadius: 0.0313 * IN, rise: 0.2502 * IN });
    for (const [th, L] of table) expect(Math.abs(threeArcLift(g, (th / 2) * DEG) / IN - L)).toBeLessThan(1.5e-4);
    expect(g.f1 / DEG).toBeCloseTo(40.3, 1);
    expect(g.f0 / DEG).toBeCloseTo(68.57, 1);
    // MTFC Table 2: seat-to-seat duration vs lash
    for (const [lash, dur] of [[0.01, 239.2], [0.02, 224.6], [0.0256, 218.0], [0.03, 213.4], [0.05, 195.5]]) {
      expect(Math.abs(duration(g, lash * IN) - dur)).toBeLessThan(0.15);
    }
  });

  it('the spec cam at its lash reproduces the Ford timing (≈ 218° seat to seat, lift = rise − lash)', () => {
    for (const v of [MODEL_T.intakeValve, MODEL_T.exhaustValve]) {
      const cam = new ValveCam(v, 0.01);
      expect(cam.geometry).not.toBeNull();
      const dur = ((v.closeDeg - v.openDeg) % 720 + 720) % 720;
      expect(Math.abs(duration(cam.geometry!, cam.lash) - dur)).toBeLessThan(0.5);
      expect(cam.valveLift(cam.centreDeg)).toBeCloseTo(v.maxLift, 9);
      // seated just outside the quoted events, open just inside them
      expect(cam.valveLift(v.openDeg - 0.6)).toBe(0);
      expect(cam.valveLift(v.openDeg + 0.6)).toBeGreaterThan(0);
      expect(cam.valveLift(v.closeDeg + 0.6)).toBe(0);
      expect(cam.valveLift(v.closeDeg - 0.6)).toBeGreaterThan(0);
    }
    expect(lobeCentreDeg(MODEL_T.intakeValve)).toBeCloseTo(-238.25, 9);
    expect(lobeCentreDeg(MODEL_T.exhaustValve)).toBeCloseTo(251.05, 9);
  });

  it('the rendered lobe outline is the exact lobe (support function within 2 µm)', () => {
    const cam = new ValveCam(MODEL_T.exhaustValve, 0.01);
    const prof = lobeProfileOf(cam);
    expect(prof.baseRadius).toBeCloseTo(cam.baseRadius, 12);
    for (let th = -360; th < 360; th += 7.3) {
      // follower along +y, cam turned by θ/2: height = R_b + s(θ − θ_c) for a lobe built around θ_c = 0
      const h = lobeHeightAt(prof.outline, (th * DEG) / 2);
      expect(Math.abs(h - (cam.baseRadius + cam.tappetLiftFromCentre(th)))).toBeLessThan(2e-6);
    }
    // flat-follower contact offset stays inside the 1 in push-rod foot (MTFC: ≥ 0.809 in needed)
    expect(2 * prof.maxContactOffset).toBeLessThan(1.0 * IN);
    expect(2 * prof.maxContactOffset).toBeGreaterThan(0.75 * IN);
  });

  it('falls back to the sin² model curve for other cam kinds', () => {
    const v = { ...MODEL_T.intakeValve, cam: undefined };
    const cam = new ValveCam(v, 0.0125);
    expect(cam.geometry).toBeNull();
    expect(cam.lash).toBe(0);
    expect(cam.baseRadius).toBe(0.0125);
    expect(cam.valveLift(cam.centreDeg)).toBeCloseTo(v.maxLift, 9);
    const prof = buildLobeProfile((th) => cam.tappetLiftFromCentre(th), cam.baseRadius);
    expect(Math.abs(lobeHeightAt(prof.outline, 0) - prof.baseRadius - v.maxLift)).toBeLessThan(5e-5);
    expect(wrapDeg720(cam.centreDeg)).toBe(cam.centreDeg);
  });
});
