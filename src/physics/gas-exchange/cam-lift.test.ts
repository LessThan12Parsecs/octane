import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/mechanics_cam_three_arc.json';
import type { CamSpec, ValveSpec } from '../core/engine-spec';
import {
  CHOI2018_EXHAUST_LIFT_START_DEG,
  CHOI2018_EXHAUST_LIFT_ZERO_LASH,
  CHOI2018_INTAKE_LIFT_START_DEG,
  CHOI2018_INTAKE_LIFT_ZERO_LASH,
} from '../engines/cfr-valve-lift';
import { CFR_F1, CFR_VALVE_LASH, cfrValveLiftChoi2018 } from '../engines/cfr';
import { MODEL_T, MODEL_T_VALVE_LASH } from '../engines/model-t';
import { TabulatedLiftProfile, ThreeArcFlatFollowerProfile, threeArcLobeGeometry } from './cam-lift';
import {
  createLiftProfile,
  polydyneCam,
  resolveValveLash,
  ValveLiftProfile,
  valveLift,
  valveLiftProfile,
  wrapDeg720,
} from './valve-lift';

const IN = 0.0254;
const DEGR = Math.PI / 180;

interface Row {
  camDeg: number;
  camLift: number;
  contactOffset: number;
}
interface Case {
  name: string;
  baseRadius: number;
  flankRadius: number;
  noseRadius: number;
  rise: number;
  baseAngleDeg: number;
  noseAngleDeg: number;
  maxContactOffset: number;
  rows: Row[];
  durations: { lash: number; durationCrankDeg: number }[];
}
const cases = (fixture as unknown as { cases: Case[] }).cases;

/** MTFC Tulsa stock lobe (design_stock.htm Table 1). */
const MTFC_LOBE = { baseRadius: 0.406 * IN, flankRadius: 1.2601 * IN, noseRadius: 0.0313 * IN, rise: 0.2502 * IN };
const timing = (lash: number) => ({ maxLift: 0.2246 * IN, openDeg: -347.3, closeDeg: -129.2, timingLiftThreshold: 0, lash });

describe('three-arc lobe on a flat follower: exact support-function lift', () => {
  for (const c of cases) {
    it(`${c.name}: arc layout, lift, contact offset and durations = numerical envelope of the polygonised outline`, () => {
      const lobe = { baseRadius: c.baseRadius, flankRadius: c.flankRadius, noseRadius: c.noseRadius, rise: c.rise };
      const g = threeArcLobeGeometry(lobe);
      expect(Math.abs(g.baseAngle / DEGR - c.baseAngleDeg)).toBeLessThan(1e-9);
      expect(Math.abs(g.noseAngle / DEGR - c.noseAngleDeg)).toBeLessThan(1e-9);
      expect(Math.abs(g.maxContactOffset - c.maxContactOffset) / c.maxContactOffset).toBeLessThan(1e-6);
      // lift at zero lash, centred at θ_c = 0 (cam angle φ ↔ crank 2φ)
      const p = new ThreeArcFlatFollowerProfile({ maxLift: c.rise, openDeg: -100, closeDeg: 100, timingLiftThreshold: 0 }, lobe);
      expect(p.centerDeg).toBe(0);
      for (const r of c.rows) {
        expect(Math.abs(p.camLift(2 * r.camDeg) - r.camLift)).toBeLessThan(2e-10);
        expect(Math.abs(p.camLift(-2 * r.camDeg) - r.camLift)).toBeLessThan(2e-10);
        // the touching vertex resolves the contact point to the polygon spacing (ρ/60000 on the flank)
        expect(Math.abs(p.contactOffset(2 * r.camDeg) - r.contactOffset)).toBeLessThan(c.flankRadius / 40000);
      }
      for (const d of c.durations) {
        const q = new ThreeArcFlatFollowerProfile({ maxLift: c.rise - d.lash, openDeg: -100, closeDeg: 100, timingLiftThreshold: 0, lash: d.lash }, lobe);
        expect(Math.abs(2 * q.halfDurationDeg - d.durationCrankDeg)).toBeLessThan(1e-6);
        expect(q.maxLift).toBeCloseTo(c.rise - d.lash, 15);
      }
    });
  }

  it('reproduces the MTFC Tulsa lift equations / table and duration-vs-lash Table 2 of the stock lobe', () => {
    // MTFC design_stock.htm: nose L = 0.2502 − 0.6250(1 − cos f), flank L = 0.8541[1 − cos(68.57° − |f|)],
    // f = cam degrees; gross lift (in) at crank degrees from the lobe centre (valvetrain research table).
    const table: [number, number][] = [
      [0, 0.2502], [10, 0.2478], [20, 0.2407], [30, 0.2289], [40, 0.2125], [50, 0.1916], [60, 0.1665],
      [70, 0.1372], [80, 0.104], [90, 0.0713], [100, 0.0445], [110, 0.0238], [120, 0.0095], [130, 0.0017],
    ];
    const p = new ThreeArcFlatFollowerProfile({ ...timing(0), openDeg: -137, closeDeg: 137 }, MTFC_LOBE);
    for (const [crank, L] of table) expect(Math.abs(p.camLift(crank) / IN - L)).toBeLessThan(1.5e-4);
    expect(p.geometry.noseAngle / DEGR).toBeCloseTo(40.3, 1);
    expect(p.geometry.baseAngle / DEGR).toBeCloseTo(68.57, 2);
    expect(p.camHalfDurationDeg).toBeCloseTo(137.1, 1); // "137.1° 0"
    // MTFC Table 2 (crank-degree seat-to-seat duration vs lash, in)
    const t2: [number, number][] = [
      [0.01, 239.2], [0.015, 231.2], [0.02, 224.6], [0.025, 218.7], [0.0256, 218.0], [0.03, 213.4], [0.035, 208.4],
      [0.04, 203.8], [0.05, 195.5], [0.065, 184.2],
    ];
    for (const [lash, dur] of t2) {
      const q = new ThreeArcFlatFollowerProfile(timing(lash * IN), MTFC_LOBE);
      expect(Math.abs(2 * q.halfDurationDeg - dur)).toBeLessThan(0.15);
    }
    // kinematics (Tilden 'Cam Design History' Model T row): max velocity 0.00706 in/cam-deg at the junction,
    // 0.0036 in/cam-deg at 0.025 in lift, flank acceleration +0.854 in/rad² at the base, nose −0.625 in/rad²,
    // minimum flat-follower diameter 0.809 in (vs the 1 in foot)
    const perCamDeg = (th: number) => Math.abs(p.liftRate(th)) * 2 / IN; // m/crank-deg → in/cam-deg
    const thJ = 2 * p.geometry.noseAngle / DEGR;
    expect(perCamDeg(thJ)).toBeCloseTo(0.00706, 4);
    const th025 = 2 * p.camAngleAtLift(0.025 * IN) / DEGR;
    expect(perCamDeg(th025)).toBeCloseTo(0.0036, 4);
    const perRad2 = (th: number) => p.liftAccel(th) / (Math.PI / 360) ** 2 / IN;
    expect(perRad2(2 * p.geometry.baseAngle / DEGR - 1e-6)).toBeCloseTo(0.8541, 4);
    expect(perRad2(0)).toBeCloseTo(-0.6249, 4);
    expect(2 * p.geometry.maxContactOffset / IN).toBeCloseTo(0.809, 2);
    expect(p.geometry.maxContactOffset).toBeLessThan(0.5 * IN);
  });

  it('valve lift = max(0, cam lift − lash); C¹ at the junctions, acceleration steps; analytic derivatives = FD', () => {
    const p = new ThreeArcFlatFollowerProfile(timing(MODEL_T_VALVE_LASH), MTFC_LOBE);
    const h = 1e-5;
    let maxE1 = 0;
    let maxE2 = 0;
    const thN = 2 * p.geometry.noseAngle / DEGR;
    for (let th = -360; th < 360; th += 0.173) {
      const l = p.lift(th);
      const lc = p.camLift(th);
      expect(l).toBeCloseTo(Math.max(0, lc - p.lash), 15);
      expect(p.isOpen(th)).toBe(l > 0);
      const x = Math.abs(wrapDeg720(th - p.centerDeg));
      if (Math.abs(x - p.halfDurationDeg) < 0.01 || Math.abs(x - thN) < 0.01) continue; // seat / junction kinks
      maxE1 = Math.max(maxE1, Math.abs((p.lift(th + h) - p.lift(th - h)) / (2 * h) - p.liftRate(th)));
      maxE2 = Math.max(maxE2, Math.abs((p.liftRate(th + h) - p.liftRate(th - h)) / (2 * h) - p.liftAccel(th)));
    }
    expect(maxE1).toBeLessThan(1e-11); // rate scale ≈ 1e-4 m/deg
    expect(maxE2).toBeLessThan(1e-9); // accel scale ≈ 4e-6 m/deg²
    // junction: lift and velocity continuous, acceleration jumps from −d cos φ_n to +(ρ − r_b)
    const thJ = p.centerDeg + thN;
    expect(Math.abs(p.lift(thJ + 1e-9) - p.lift(thJ - 1e-9))).toBeLessThan(1e-12);
    expect(Math.abs(p.liftRate(thJ + 1e-9) - p.liftRate(thJ - 1e-9))).toBeLessThan(1e-12);
    expect(p.liftAccel(thJ - 1e-6)).toBeLessThan(0);
    expect(p.liftAccel(thJ + 1e-6)).toBeGreaterThan(0);
    // seat: the valve leaves with the finite flank velocity
    expect(p.lift(p.seatOpenDeg - 1e-9)).toBe(0);
    expect(p.liftRate(p.seatOpenDeg + 1e-6)).toBeGreaterThan(1e-5);
    expect(p.liftRate(p.seatCloseDeg - 1e-6)).toBeLessThan(-1e-5);
    // symmetry about the lobe centre; 720° periodicity
    for (const d of [5, 40, 80, 105]) expect(p.lift(p.centerDeg + d)).toBeCloseTo(p.lift(p.centerDeg - d), 15);
    for (const th of [-300, -238, -150, 100]) expect(p.lift(th + 720)).toBeCloseTo(p.lift(th), 15);
    // a timing threshold above the seat: the quoted events are met at that lift (± timingErrorDeg)
    const thr = 0.006 * IN;
    const q = new ThreeArcFlatFollowerProfile({ ...timing(MODEL_T_VALVE_LASH), timingLiftThreshold: thr }, MTFC_LOBE);
    expect(q.lift(-347.3 - q.timingErrorDeg)).toBeCloseTo(thr, 13);
    expect(q.lift(-129.2 + q.timingErrorDeg)).toBeCloseTo(thr, 13);
    expect(q.centerDeg).toBeCloseTo(-238.25, 12);
  });

  it('rejects radii that admit no tangent three-arc lobe', () => {
    expect(() => threeArcLobeGeometry({ ...MTFC_LOBE, flankRadius: 0.3 * IN })).toThrow();
    expect(() => threeArcLobeGeometry({ ...MTFC_LOBE, rise: 0 })).toThrow();
    expect(() => new ThreeArcFlatFollowerProfile(timing(0.3 * IN), MTFC_LOBE)).toThrow();
  });
});

describe('Model T valve lift (MODEL_T cam + lash; valvetrain research)', () => {
  for (const kind of ['intake', 'exhaust'] as const) {
    it(`${kind}: ≈ 218° seat-to-seat and 0.225 in net lift at 0.0256 in, ≈ 212° at 1/32 in, events at the quoted angles`, () => {
      const v = kind === 'intake' ? MODEL_T.intakeValve : MODEL_T.exhaustValve;
      const p = createLiftProfile(v);
      expect(p).toBeInstanceOf(ThreeArcFlatFollowerProfile);
      const t = p as ThreeArcFlatFollowerProfile;
      expect(p.lash).toBe(MODEL_T_VALVE_LASH);
      expect(Math.abs(2 * t.halfDurationDeg - 218.0)).toBeLessThan(0.1);
      expect(Math.abs(p.maxLift / IN - 0.225)).toBeLessThan(0.001);
      expect(p.maxLift).toBeCloseTo(v.maxLift, 15); // the spec's quoted lift
      expect(p.lift(p.centerDeg)).toBeCloseTo(v.maxLift, 12);
      // quoted Ford events (218.1° / 217.9°) met within 0.1° each side
      expect(Math.abs(t.timingErrorDeg)).toBeLessThan(0.1);
      expect(Math.abs(wrapDeg720(p.seatOpenDeg - v.openDeg))).toBeLessThan(0.1);
      expect(Math.abs(wrapDeg720(p.seatCloseDeg - v.closeDeg))).toBeLessThan(0.1);
      // lobe centres (valvetrain research: intake −238.25°, exhaust +251.06°)
      expect(p.centerDeg).toBeCloseTo(kind === 'intake' ? -238.25 : 251.0, 1);
      // Ford's clearance limits: 1/32 in → 212.1°, 1/64 in → 230.4° (verified recomputation)
      const p32 = createLiftProfile(v, IN / 32) as ThreeArcFlatFollowerProfile;
      expect(Math.abs(2 * p32.halfDurationDeg - 212.1)).toBeLessThan(0.1);
      const p64 = createLiftProfile(v, IN / 64) as ThreeArcFlatFollowerProfile;
      expect(Math.abs(2 * p64.halfDurationDeg - 230.4)).toBeLessThan(0.1);
      expect(p32.maxLift).toBeLessThan(p.maxLift);
    });
  }

  it('the exhaust event wraps across the gas-exchange TDC; no overlap with the intake', () => {
    const iv = createLiftProfile(MODEL_T.intakeValve);
    const ev = createLiftProfile(MODEL_T.exhaustValve);
    expect(ev.lift(359)).toBeGreaterThan(0);
    expect(ev.lift(-359)).toBe(0);
    expect(ev.lift(-360 + 1e-6)).toBeGreaterThan(0); // closes 0.04° after TDC (geometric 217.9 vs quoted 217.9°)
    for (let th = -360; th < 360; th += 0.25) expect(iv.lift(th) > 0 && ev.lift(th) > 0).toBe(false);
  });
});

describe('tabulated lift (CamSpec "table")', () => {
  it('reproduces cfrValveLiftChoi2018 exactly (both valves, several clearances)', () => {
    const tabs = [
      { valve: 'intake', startDeg: CHOI2018_INTAKE_LIFT_START_DEG, zeroLashLift: CHOI2018_INTAKE_LIFT_ZERO_LASH },
      { valve: 'exhaust', startDeg: CHOI2018_EXHAUST_LIFT_START_DEG, zeroLashLift: CHOI2018_EXHAUST_LIFT_ZERO_LASH },
    ] as const;
    for (const tb of tabs) {
      for (const lash of [0, CFR_VALVE_LASH, 0.003 * IN]) {
        const p = new TabulatedLiftProfile({ startDeg: tb.startDeg, stepDeg: 1, zeroLashLift: tb.zeroLashLift }, lash);
        for (let th = -1080; th <= 1080; th += 0.0371) expect(p.lift(th)).toBe(cfrValveLiftChoi2018(tb.valve, th, lash));
        for (const th of [-360, -359.5, 0, 108, 359.999, 360, 467.5]) expect(p.lift(th)).toBe(cfrValveLiftChoi2018(tb.valve, th, lash));
        // seat events where the lift crosses zero; centre at the table maximum
        expect(p.lift(p.seatOpenDeg - 1e-7)).toBe(0);
        expect(p.lift(p.seatOpenDeg + 1e-3)).toBeGreaterThan(0);
        expect(p.lift(p.seatCloseDeg + 1e-7)).toBe(0);
        expect(p.lift(p.seatCloseDeg - 1e-3)).toBeGreaterThan(0);
        expect(p.maxLift).toBeCloseTo(Math.max(...tb.zeroLashLift) - lash, 15);
        expect(p.lift(p.centerDeg)).toBeCloseTo(p.maxLift, 15);
        // rate = segment slope
        for (const th of [-300.3, -200.7, 200.2, 300.6]) {
          const fd = (p.lift(th + 1e-4) - p.lift(th - 1e-4)) / 2e-4;
          if (p.lift(th) > 0) expect(Math.abs(fd - p.liftRate(th))).toBeLessThan(1e-9);
        }
      }
    }
    // with the CFR clearance the table's seat events are the measured Choi 2018 events (≈ −344/−152, 141/−347)
    const iv = new TabulatedLiftProfile({ startDeg: CHOI2018_INTAKE_LIFT_START_DEG, stepDeg: 1, zeroLashLift: CHOI2018_INTAKE_LIFT_ZERO_LASH }, CFR_VALVE_LASH);
    expect(Math.abs(iv.seatOpenDeg - CFR_F1.intakeValve.openDeg)).toBeLessThan(1);
    expect(Math.abs(iv.seatCloseDeg - CFR_F1.intakeValve.closeDeg)).toBeLessThan(1);
  });

  it('non-unit step and validation', () => {
    const t = { startDeg: -100, stepDeg: 2, zeroLashLift: [0, 1e-3, 3e-3, 4e-3, 3e-3, 1e-3, 0] };
    const p = new TabulatedLiftProfile(t, 5e-4);
    expect(p.lift(-97)).toBeCloseTo(2e-3 - 5e-4, 15);
    expect(p.camLift(-97)).toBeCloseTo(2e-3, 15);
    expect(p.seatOpenDeg).toBeCloseTo(-99, 12);
    expect(p.seatCloseDeg).toBeCloseTo(-89, 12);
    expect(p.centerDeg).toBe(-94);
    expect(p.liftRate(-97)).toBeCloseTo(1e-3, 15);
    expect(p.liftAccel(-97)).toBe(0);
    expect(() => new TabulatedLiftProfile({ ...t, stepDeg: 200 })).toThrow();
    expect(() => new TabulatedLiftProfile(t, 5e-3)).toThrow();
  });
});

describe('createLiftProfile / resolveValveLash / cam-aware functional API', () => {
  const cfrIv = CFR_F1.intakeValve;

  it('no cam: bit-identical to the historical cycle-model construction', () => {
    for (const lash of [0, CFR_VALVE_LASH]) {
      const a = createLiftProfile(cfrIv, lash);
      const b = new ValveLiftProfile({ ...cfrIv, lash });
      expect(a).toBeInstanceOf(ValveLiftProfile);
      for (let th = -360; th < 360; th += 0.0913) {
        expect(a.lift(th)).toBe(b.lift(th));
        expect(a.liftRate(th)).toBe(b.liftRate(th));
      }
      expect(a.seatOpenDeg).toBe(b.seatOpenDeg);
      expect(a.seatCloseDeg).toBe(b.seatCloseDeg);
    }
  });

  it('lash precedence: explicit override > ValveSpec.lash > 0; resolveValveLash: spec > option', () => {
    expect(resolveValveLash(cfrIv, CFR_VALVE_LASH)).toBe(CFR_VALVE_LASH); // CFR spec has no lash
    expect(resolveValveLash(MODEL_T.intakeValve, CFR_VALVE_LASH)).toBe(MODEL_T_VALVE_LASH);
    expect(resolveValveLash(cfrIv)).toBe(0);
    expect(createLiftProfile(cfrIv).lash).toBe(0);
    expect(createLiftProfile(MODEL_T.intakeValve).lash).toBe(MODEL_T_VALVE_LASH);
    expect(createLiftProfile(MODEL_T.intakeValve, 1e-4).lash).toBe(1e-4);
  });

  it('dispatches polydyne powers and tables', () => {
    const poly: ValveSpec = { ...cfrIv, cam: { kind: 'polydyne', powers: [2, 6] } };
    const a = createLiftProfile(poly, CFR_VALVE_LASH);
    const b = new ValveLiftProfile({ ...cfrIv, lash: CFR_VALVE_LASH }, polydyneCam([2, 6]));
    for (const th of [-340, -300, -248, -160]) expect(a.lift(th)).toBe(b.lift(th));
    const dflt = createLiftProfile({ ...cfrIv, cam: { kind: 'polydyne' } }, CFR_VALVE_LASH);
    expect(dflt.lift(-300)).toBe(new ValveLiftProfile({ ...cfrIv, lash: CFR_VALVE_LASH }).lift(-300));
    const cam: CamSpec = { kind: 'table', startDeg: CHOI2018_INTAKE_LIFT_START_DEG, stepDeg: 1, zeroLashLift: CHOI2018_INTAKE_LIFT_ZERO_LASH };
    const tab = createLiftProfile({ ...cfrIv, cam }, CFR_VALVE_LASH);
    expect(tab).toBeInstanceOf(TabulatedLiftProfile);
    expect(tab.lift(-250.3)).toBe(cfrValveLiftChoi2018('intake', -250.3));
  });

  it('functional API follows the spec cam and lash', () => {
    const v: ValveSpec = { ...MODEL_T.intakeValve };
    const ref = createLiftProfile(MODEL_T.intakeValve);
    expect(valveLift(v, -250)).toBe(ref.lift(-250));
    const a = valveLiftProfile(v);
    expect(valveLiftProfile(v)).toBe(a);
    v.cam = undefined;
    expect(valveLiftProfile(v)).not.toBe(a);
    expect(valveLiftProfile(v)).toBeInstanceOf(ValveLiftProfile);
  });
});
