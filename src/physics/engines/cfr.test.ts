import { describe, expect, it } from 'vitest';
import fixture from '../../../test/fixtures/mechanics_cfr_conditions.json';
import { SliderCrank } from '../mechanics/kinematics';
import { CrankTrainDynamics } from '../mechanics/dynamics';
import { pnhFmep } from '../mechanics/friction';
import liftFixture from '../../../test/fixtures/mechanics_cfr_valve_lift.json';
import {
  CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER,
  CFR_COUNTER_STEP,
  CFR_CRANKCASE_GAUGE_PRESSURE,
  CFR_F1,
  CFR_FRICTION,
  CFR_GUIDE_COUNTER_PRF90,
  CFR_KNOCK_PICKUP,
  CFR_MON_CONDITIONS,
  CFR_RON_CONDITIONS,
  CFR_VALVE_LASH,
  CFR_VALVE_TIMING_ASTM,
  CFR_VALVE_TIMING_CHOI2018,
  cfrCompressionRatioAtCounter,
  cfrCompressionRatioFromCounter,
  cfrCounterAtCompressionRatio,
  cfrCounterFromCompressionRatio,
  cfrMonSparkAdvanceDeg,
  cfrValveLiftChoi2018,
  clearanceHeightAtTDC,
  compressionRatioForClearanceHeight,
  cylinderRaise,
  displacedVolume,
  headFacePosition,
} from './cfr';

const fx = fixture as {
  humidity: { grainsPerLb: number; T: number; humidityRatio: number; relativeHumidity: number }[];
  ronMixture: { T_mixture: number };
  counter: {
    counts: number[];
    compressionRatio: number[];
    rigidRaise: {
      counterStep: number;
      basicCounter: number;
      effectiveClearanceHeightAtBasicCounter: number;
      guideTableCounters: Record<string, number>;
      guideTableCR: Record<string, number>;
      guideTableCRChoiPolynomial: Record<string, number>;
      monSparkTable: { counter: number; advanceDeg: number; compressionRatio: number }[];
    };
  };
};
interface LiftRec {
  startDeg: number;
  peakLiftIn: number;
  'cross_0.008': [number, number];
  'cross_0.054': [number, number];
}
const lf = liftFixture as unknown as { intake: LiftRec; exhaust: LiftRec };
const IN = 0.0254;

const g = CFR_F1.geometry;
const R = g.bore / 2;

/** Wrap a crank angle (deg) into [-360, 360). */
const wrap = (d: number) => ((((d + 360) % 720) + 720) % 720) - 360;

describe('CFR F-1 geometry', () => {
  it('bore, stroke, rod and displacement match the literature', () => {
    expect(g.bore).toBeCloseTo(0.08255, 12); // 82.55 mm
    expect(g.stroke).toBeCloseTo(0.1143, 12); // 114.3 mm
    expect(g.conRodLength).toBe(0.254);
    // 37.33 in³ [KW17] = 611.7 cm³
    expect(displacedVolume() / 1e-6).toBeCloseTo(611.7, 1);
    expect(displacedVolume() / (0.0254 ** 3)).toBeCloseTo(37.33, 2);
  });

  it('clearance height vs CR is consistent with SliderCrank and invertible', () => {
    for (const cr of [4, 5, 6.5, 7.55, 12, 18]) {
      const sc = new SliderCrank(g, cr);
      expect(clearanceHeightAtTDC(cr)).toBeCloseTo(sc.clearanceHeightTDC, 15);
      expect(compressionRatioForClearanceHeight(clearanceHeightAtTDC(cr))).toBeCloseTo(cr, 12);
    }
    // CR 18 → about 6.6 mm, CR 4 → about 38 mm
    expect(clearanceHeightAtTDC(18)).toBeGreaterThan(0.006);
    expect(clearanceHeightAtTDC(4)).toBeLessThan(0.04);
  });

  it('head-face position and cylinder raise: lowest at max CR, rises ≈ 31 mm from CR 18 to CR 4', () => {
    const [crMin, crMax] = g.compressionRatioRange;
    expect(cylinderRaise(crMax)).toBe(0);
    expect(cylinderRaise(crMin)).toBeCloseTo(clearanceHeightAtTDC(crMin) - clearanceHeightAtTDC(crMax), 15);
    expect(cylinderRaise(crMin)).toBeGreaterThan(0.028);
    expect(cylinderRaise(crMin)).toBeLessThan(0.034);
    let prev = Infinity;
    for (let cr = crMin; cr <= crMax; cr += 0.5) {
      const y = headFacePosition(cr);
      expect(y).toBeLessThan(prev);
      prev = y;
    }
    // head face sits above the crown at TDC by exactly h_TDC
    const sc = new SliderCrank(g, 7);
    expect(headFacePosition(7) - (sc.pinHeight(0) + g.compressionHeight)).toBeCloseTo(sc.clearanceHeightTDC, 15);
  });

  it('spark gap lies inside the chamber at every compression ratio, plug side-mounted', () => {
    const [x, y, z] = CFR_F1.sparkPlug.gapCenter;
    expect(Math.hypot(x, z)).toBeLessThan(R);
    expect(y).toBeLessThan(0);
    expect(-y + CFR_F1.sparkPlug.gap / 2).toBeLessThan(clearanceHeightAtTDC(g.compressionRatioRange[1]));
    const [ax, ay, az] = CFR_F1.sparkPlug.axis;
    expect(Math.hypot(ax, ay, az)).toBeCloseTo(1, 15);
    expect(ay).toBe(0); // horizontal plug
    // knock pickup on the opposite side of the chamber, inside the bore
    expect(Math.hypot(CFR_KNOCK_PICKUP.position[0], CFR_KNOCK_PICKUP.position[2])).toBeLessThan(R);
    expect(Math.sign(CFR_KNOCK_PICKUP.position[2])).toBe(-Math.sign(z));
  });

  it('ASTM D2699 §10.2.5: intake shroud toward the spark-plug side, swirl counter-clockwise seen from above', () => {
    const iv = CFR_F1.intakeValve;
    const [px, , pz] = CFR_F1.sparkPlug.gapCenter;
    // masked side (opposite the opening) points toward the plug
    const mx = Math.cos(iv.shroudDirection + Math.PI);
    const mz = Math.sin(iv.shroudDirection + Math.PI);
    expect(mx * (px - iv.position[0]) + mz * (pz - iv.position[1])).toBeGreaterThan(0);
    // charge leaves through the opening at the valve position: L_y = z v_x − x v_z > 0 ⇔ CCW about +y (viewed from +y)
    const vx = Math.cos(iv.shroudDirection);
    const vz = Math.sin(iv.shroudDirection);
    expect(iv.position[1] * vx - iv.position[0] * vz).toBeGreaterThan(0);
    // plug gap and knock pickup lie on the diameter perpendicular to the valves
    expect(Math.abs(px)).toBeLessThan(1e-12);
    expect(Math.abs(iv.position[1])).toBeLessThan(1e-12);
  });

  it('spark plug per ASTM D2699 §10.3.16: gap 0.51 mm', () => {
    expect(CFR_F1.sparkPlug.gap).toBeCloseTo(0.51e-3, 5);
  });

  it('valves fit in the bore without overlapping; intake 180° shroud', () => {
    const iv = CFR_F1.intakeValve;
    const ev = CFR_F1.exhaustValve;
    for (const v of [iv, ev]) {
      expect(Math.hypot(v.position[0], v.position[1]) + v.headDiameter / 2).toBeLessThan(R);
      expect(v.seatInnerDiameter).toBeLessThan(v.headDiameter);
    }
    const d = Math.hypot(iv.position[0] - ev.position[0], iv.position[1] - ev.position[1]);
    expect(d).toBeGreaterThan((iv.headDiameter + ev.headDiameter) / 2);
    expect(iv.shroudArcDeg).toBe(180);
    expect(ev.shroudArcDeg).toBe(0);
  });

  it('valve timing: ASTM nominal cam events (D2699 §10.2.3) and the measured seat-off events used by the spec', () => {
    const t = CFR_VALVE_TIMING_ASTM;
    expect(t.intakeOpenDeg).toBe(-360 + 10);
    expect(t.intakeCloseDeg).toBe(-180 + 34);
    expect(t.exhaustOpenDeg).toBe(180 - 40);
    expect(wrap(t.exhaustCloseDeg)).toBe(wrap(360 + 15));
    expect(wrap(t.exhaustCloseDeg - t.intakeOpenDeg)).toBe(5); // "camshaft overlap 5°" (D2699 Table 1)
    // the spec uses the measured valve events at the 0.008 in running clearance (Choi et al. 2018 Table 2)
    const iv = CFR_F1.intakeValve;
    const ev = CFR_F1.exhaustValve;
    expect(iv.openDeg).toBe(CFR_VALVE_TIMING_CHOI2018.intakeOpenDeg);
    expect(iv.closeDeg).toBe(CFR_VALVE_TIMING_CHOI2018.intakeCloseDeg);
    expect(ev.openDeg).toBe(CFR_VALVE_TIMING_CHOI2018.exhaustOpenDeg);
    expect(ev.closeDeg).toBe(CFR_VALVE_TIMING_CHOI2018.exhaustCloseDeg);
    expect([iv.timingLiftThreshold, ev.timingLiftThreshold]).toEqual([0, 0]);
    expect(wrap(ev.closeDeg - iv.openDeg)).toBe(-3); // no real overlap at running clearance
    // max lift = ASTM zero-lash valve lift 0.238 in (§10.2.4) − hot clearance 0.008 in (§10.3.2.1), not the 0.246 in lobe rise
    expect(iv.maxLift).toBeCloseTo(0.23 * IN, 12);
    expect(ev.maxLift).toBeCloseTo(0.23 * IN, 12);
    expect(CFR_FRICTION.maxValveLift).toBeCloseTo(0.23 * IN, 12);
  });

  it('masses build a valid crank train', () => {
    const dyn = new CrankTrainDynamics(new SliderCrank(g, 7), CFR_F1.masses);
    // realistic rod: I_g below the two-mass-equivalent value
    expect(dyn.rodInertiaDefect).toBeLessThan(0);
  });
});

describe('CFR digital counter (Choi et al. 2018, Fig. 4; ASTM D2699 A2.2)', () => {
  it('matches the published polynomial and inverts', () => {
    fx.counter.counts.forEach((c, i) => {
      expect(cfrCompressionRatioFromCounter(c)).toBeCloseTo(fx.counter.compressionRatio[i], 12);
      expect(cfrCounterFromCompressionRatio(fx.counter.compressionRatio[i])).toBeCloseTo(c, 8);
    });
  });

  it('rigid-raise model: 0.0007 in per digit, fitted intercept, consistent with the kinematics', () => {
    const rr = fx.counter.rigidRaise;
    expect(CFR_COUNTER_STEP).toBeCloseTo(rr.counterStep, 15);
    expect(CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER).toBeCloseTo(rr.effectiveClearanceHeightAtBasicCounter, 14);
    // dial indicator = 1.012 − counter/1410 in (D2699 Table A4 footnote) ⇔ 1/1410 in ≈ 0.0007 in per digit
    expect(1 / 1410 / 0.0007).toBeCloseTo(1, 1);
    for (const c of [264, 500, 726, 930, 1145, 1290]) {
      const cr = cfrCompressionRatioAtCounter(c);
      expect(cfrCounterAtCompressionRatio(cr)).toBeCloseTo(c, 9);
      // raising the cylinder by one digit changes the clearance volume by exactly A·0.0007 in
      const sc1 = new SliderCrank(g, cr);
      const sc2 = new SliderCrank(g, cfrCompressionRatioAtCounter(c + 1));
      expect((sc1.clearanceVolume - sc2.clearanceVolume) / sc1.boreArea).toBeCloseTo(CFR_COUNTER_STEP, 14);
    }
    // agrees with the oil-measured polynomial in the rating region (RON/MON PRF90 settings)
    for (const k of ['RON90', 'MON90']) {
      expect(cfrCompressionRatioAtCounter(rr.guideTableCounters[k])).toBeCloseTo(rr.guideTableCR[k], 12);
      expect(Math.abs(rr.guideTableCR[k] - rr.guideTableCRChoiPolynomial[k])).toBeLessThan(0.03);
    }
  });
});

describe('CFR measured valve lift (Choi et al. 2018 Fig. 3, digitised)', () => {
  it('peak lift ≈ ASTM 0.238 ± 0.002 in at zero lash', () => {
    for (const v of ['intake', 'exhaust'] as const) {
      let peak = 0;
      for (let th = -360; th < 360; th += 0.25) peak = Math.max(peak, cfrValveLiftChoi2018(v, th, 0));
      expect(peak / IN).toBeCloseTo(lf[v].peakLiftIn, 4);
      expect(Math.abs(peak / IN - 0.238)).toBeLessThan(0.003);
    }
  });

  it('crossings at the 0.008 in running clearance reproduce Choi Table 2 (≤ 1.5°) and the spec events', () => {
    const cross = (v: 'intake' | 'exhaust', thr: number, from: number): [number, number] => {
      let open = NaN;
      let close = NaN;
      let prev = cfrValveLiftChoi2018(v, from, 0) > thr;
      for (let k = 1; k <= 7200; k++) {
        const th = from + k * 0.1;
        const on = cfrValveLiftChoi2018(v, th, 0) > thr;
        if (on && !prev) open = th;
        if (!on && prev) close = th;
        prev = on;
      }
      return [wrap(open), wrap(close)];
    };
    const [io, ic] = cross('intake', 0.008 * IN, -300);
    const [eo, ec] = cross('exhaust', 0.008 * IN, 0);
    expect(Math.abs(io - CFR_F1.intakeValve.openDeg)).toBeLessThan(1.5);
    expect(Math.abs(ic - CFR_F1.intakeValve.closeDeg)).toBeLessThan(1.5);
    expect(Math.abs(eo - CFR_F1.exhaustValve.openDeg)).toBeLessThan(1.5);
    expect(Math.abs(ec - CFR_F1.exhaustValve.closeDeg)).toBeLessThan(1.5);
    // ASTM D2699 A2.1.2 timing check: 0.054 in rise 30 ± 2° after TDC (the ANL engine is in spec)
    const [io54] = cross('intake', 0.054 * IN, -300);
    expect(Math.abs(io54 - (-360 + 30))).toBeLessThan(2);
    expect(io).toBeCloseTo(lf.intake['cross_0.008'][0], 0);
  });

  it('lash subtracts from the zero-lash profile; periodic in 720°; zero outside the events', () => {
    const th = -250;
    expect(cfrValveLiftChoi2018('intake', th)).toBeCloseTo(cfrValveLiftChoi2018('intake', th, 0) - CFR_VALVE_LASH, 15);
    expect(cfrValveLiftChoi2018('intake', th + 720)).toBe(cfrValveLiftChoi2018('intake', th));
    expect(cfrValveLiftChoi2018('exhaust', -355)).toBeGreaterThan(0); // EVC after gas-exchange TDC (≡ 365°)
    expect(cfrValveLiftChoi2018('exhaust', 365)).toBe(cfrValveLiftChoi2018('exhaust', -355));
    expect(cfrValveLiftChoi2018('intake', 0)).toBe(0);
    expect(cfrValveLiftChoi2018('exhaust', 0)).toBe(0);
    expect(cfrValveLiftChoi2018('intake', Number.NaN)).toBe(0);
  });
});

describe('CFR standard operating conditions', () => {
  it('RON (ASTM D2699): 600 rpm, 13° BTDC, 52 °C air, PRF 90, fixed speed, WOT', () => {
    const op = CFR_RON_CONDITIONS;
    expect(op.rpm).toBe(600);
    expect(op.sparkAdvanceDeg).toBe(13);
    expect(op.ambientTemperature).toBeCloseTo(325.15, 12);
    expect(op.fuel).toEqual({ kind: 'PRF', octaneNumber: 90 });
    expect(op.speedMode).toBe('fixed');
    expect(op.throttle).toBe(1);
    expect(op.coolantTemperature).toBeCloseTo(373.15, 12);
    expect(op.equivalenceRatio).toBeGreaterThanOrEqual(1.05);
    expect(op.equivalenceRatio).toBeLessThanOrEqual(1.12);
    // guide-table standard-knock setting for 90 O.N. (D2699 Table A4.1: counter 726)
    expect(CFR_GUIDE_COUNTER_PRF90.RON).toBe(726);
    expect(op.compressionRatio).toBeCloseTo(fx.counter.rigidRaise.guideTableCR.RON90, 12);
    expect(op.compressionRatio).toBeGreaterThan(6.3);
    expect(op.compressionRatio).toBeLessThan(6.6);
    // humidity and mixture temperature from the Cantera script
    const h = fx.humidity.find((r) => r.grainsPerLb === 37.5 && Math.abs(r.T - 325.15) < 1e-9)!;
    expect(op.relativeHumidity).toBeCloseTo(h.relativeHumidity, 4);
    expect(op.intakeMixtureTemperature).toBeCloseTo(fx.ronMixture.T_mixture, 1);
  });

  it('MON (ASTM D2700): 900 rpm, 149 °C mixture, 38 °C air, CR-dependent spark 14–26° BTDC', () => {
    const op = CFR_MON_CONDITIONS;
    expect(op.rpm).toBe(900);
    expect(op.intakeMixtureTemperature).toBeCloseTo(422.15, 12);
    expect(op.ambientTemperature).toBeCloseTo(311.15, 12);
    expect(op.sparkAdvanceDeg).toBeCloseTo(cfrMonSparkAdvanceDeg(op.compressionRatio), 15);
    const h = fx.humidity.find((r) => r.grainsPerLb === 37.5 && Math.abs(r.T - 311.15) < 1e-9)!;
    expect(op.relativeHumidity).toBeCloseTo(h.relativeHumidity, 4);
    // guide-table standard-knock setting for 90 O.N. (D2700 Table A4.1, 9/16 in venturi: counter 749)
    expect(CFR_GUIDE_COUNTER_PRF90.MON).toBe(749);
    expect(op.compressionRatio).toBeCloseTo(fx.counter.rigidRaise.guideTableCR.MON90, 12);
    expect(op.sparkAdvanceDeg).toBeCloseTo(26 - (12 * (749 - 264)) / (1145 - 264), 9); // 19.4°
    // D2700 Table 3, every row (counter → CR by the rigid-raise model)
    for (const row of fx.counter.rigidRaise.monSparkTable) {
      expect(cfrMonSparkAdvanceDeg(row.compressionRatio)).toBeCloseTo(row.advanceDeg, 1);
    }
    expect(cfrMonSparkAdvanceDeg(4.5)).toBe(26);
    expect(cfrMonSparkAdvanceDeg(12)).toBe(14);
    let prev = Infinity;
    for (let cr = 4; cr <= 18; cr += 0.25) {
      const a = cfrMonSparkAdvanceDeg(cr);
      expect(a).toBeLessThanOrEqual(prev);
      prev = a;
    }
  });
});

describe('CFR friction (PNH with CFR inputs)', () => {
  it('FMEP at 600 / 900 rpm ≈ 1.22 / 1.26 bar, dominated by the piston group and valvetrain', () => {
    const ron = pnhFmep({ ...CFR_FRICTION, compressionRatio: CFR_RON_CONDITIONS.compressionRatio }, 600, 101325, 101325);
    const mon = pnhFmep({ ...CFR_FRICTION, compressionRatio: CFR_MON_CONDITIONS.compressionRatio }, 900, 101325, 101325);
    // model estimate only (no published CFR motoring-friction data found); regression band
    expect(ron.total).toBeGreaterThan(1.1e5);
    expect(ron.total).toBeLessThan(1.35e5);
    expect(mon.total).toBeGreaterThan(ron.total);
    expect(ron.reciprocating).toBeGreaterThan(ron.crankshaft);
    expect(CFR_FRICTION.viscosityRatio).toBeGreaterThan(4); // SAE 30 at 57 °C vs 10.6 cSt
    expect(CFR_FRICTION.viscosityRatio).toBeLessThan(5);
    // √(ν/ν₀) scaling (Sandoval 2002 eq. 3): skirt term = PNH × √4.37
    const pnh = pnhFmep({ ...CFR_FRICTION, viscosityRatio: 1 }, 600, 101325, 101325);
    expect(ron.pistonSkirt / pnh.pistonSkirt).toBeCloseTo(Math.sqrt(CFR_FRICTION.viscosityRatio as number), 12);
  });

  it('crankcase is at a slight vacuum (D2699 §10.3.10: 25–150 mm H2O)', () => {
    expect(CFR_CRANKCASE_GAUGE_PRESSURE).toBeLessThan(-25 * 9.80665);
    expect(CFR_CRANKCASE_GAUGE_PRESSURE).toBeGreaterThan(-150 * 9.80665);
  });
});
