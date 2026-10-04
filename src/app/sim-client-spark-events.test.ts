/**
 * interpolateSnapshot across ignition-event restarts: the cumulative per-event counters of the spark
 * record (spark.breakdownCount, spark.energyDelivered; snapshot.ts) restart at the next event, and a
 * blend across that restart would invent breakdowns / delivered energy that every consumer differencing
 * the counters (render/combustion/state.ts trembler flashes and column power, ui/spark-capture.ts) would
 * show. Oracle for the playback tests: the same differencing rule applied to the RAW stream.
 */
import { describe, expect, it } from 'vitest';
import type { EngineSnapshot, SparkPhase } from '../physics/core/snapshot';
import { MODEL_T } from '../physics/engines/model-t';
import { CombustionVisualState } from '../render/combustion/state';
import { createEmptyCylinderSnapshot, createEmptySnapshot, interpolateSnapshot, SnapshotBuffer, sparkEventRestarted } from './sim-client';

type Spark = EngineSnapshot['spark'];

function spark(p: Partial<Spark>): Spark {
  return {
    phase: 'off',
    primaryCurrent: 0,
    secondaryVoltage: 0,
    secondaryCurrent: 0,
    energyDelivered: 0,
    breakdownVoltage: 6e3,
    breakdownCount: undefined,
    pointsOpen: undefined,
    timerClosed: undefined,
    firstSparkDeg: undefined,
    primaryVoltage: undefined,
    ...p,
  };
}

function snap(t: number, sp: Partial<Spark>): EngineSnapshot {
  const s = createEmptySnapshot();
  s.t = t;
  s.thetaDeg = -40 + t; // no wrap
  s.spark = spark(sp);
  return s;
}

/** End of one event (timer open, 7 sparks, 30 mJ) → first sample of the next (timer make, counters reset). */
const endOfEvent = { phase: 'done' as SparkPhase, breakdownCount: 7, energyDelivered: 30e-3, timerClosed: false, pointsOpen: false, firstSparkDeg: -21.5 };
const newEvent = { phase: 'charging' as SparkPhase, breakdownCount: 0, energyDelivered: 0, timerClosed: true, pointsOpen: false, firstSparkDeg: NaN };

describe('sparkEventRestarted', () => {
  it('is true only when the cumulative breakdown counter went down', () => {
    expect(sparkEventRestarted(spark({ breakdownCount: 7 }), spark({ breakdownCount: 0 }))).toBe(true);
    expect(sparkEventRestarted(spark({ breakdownCount: 7 }), spark({ breakdownCount: 7 }))).toBe(false);
    expect(sparkEventRestarted(spark({ breakdownCount: 2 }), spark({ breakdownCount: 3 }))).toBe(false);
    // the CFR stream carries no counter: never a restart
    expect(sparkEventRestarted(spark({ energyDelivered: 1 }), spark({ energyDelivered: 0 }))).toBe(false);
    expect(sparkEventRestarted(spark({ breakdownCount: 3 }), spark({}))).toBe(false);
  });
});

describe('interpolateSnapshot: ignition event restarts', () => {
  it("takes the new event's spark record once past the earlier sample", () => {
    const a = snap(1, endOfEvent);
    const b = snap(2, newEvent);
    for (const u of [1e-9, 0.25, 0.5, 0.999, 1]) {
      const o = interpolateSnapshot(a, b, u, createEmptySnapshot()).spark;
      expect(o.breakdownCount, `u=${u}`).toBe(0);
      expect(o.energyDelivered, `u=${u}`).toBe(0);
      expect(o.timerClosed).toBe(true);
      expect(o.pointsOpen).toBe(false);
      expect(o.phase).toBe('charging');
      expect(o.firstSparkDeg).toBeNaN();
    }
    // exactly at the earlier sample: the earlier sample
    const o0 = interpolateSnapshot(a, b, 0, createEmptySnapshot()).spark;
    expect(o0.breakdownCount).toBe(7);
    expect(o0.energyDelivered).toBe(30e-3);
    expect(o0.timerClosed).toBe(false);
    expect(o0.phase).toBe('done');
    expect(o0.firstSparkDeg).toBe(-21.5);
    // continuous fields still blend
    a.spark.primaryCurrent = 0;
    b.spark.primaryCurrent = 4;
    expect(interpolateSnapshot(a, b, 0.25, createEmptySnapshot()).spark.primaryCurrent).toBeCloseTo(1, 15);
  });

  it('holds the counter (always an integer) and blends the energy within one event', () => {
    const a = snap(1, { phase: 'glow', breakdownCount: 2, energyDelivered: 4e-3, timerClosed: true, pointsOpen: true, firstSparkDeg: -21.5 });
    const b = snap(2, { phase: 'breakdown', breakdownCount: 3, energyDelivered: 6e-3, timerClosed: true, pointsOpen: true, firstSparkDeg: -21.5 });
    const o = interpolateSnapshot(a, b, 0.75, createEmptySnapshot()).spark;
    expect(o.breakdownCount).toBe(2);
    expect(o.phase).toBe('glow');
    expect(o.energyDelivered).toBeCloseTo(5.5e-3, 15);
    expect(o.firstSparkDeg).toBe(-21.5);
  });

  it('never blends energyDelivered across its reset on a counter-less (CFR) stream; the rest is unchanged', () => {
    const a = snap(1, { phase: 'done', energyDelivered: 40e-3 });
    const b = snap(2, { phase: 'charging', energyDelivered: 0 });
    const o = interpolateSnapshot(a, b, 0.4, createEmptySnapshot()).spark;
    expect(o.energyDelivered).toBe(0);
    expect(o.phase).toBe('done'); // discrete fields: earlier sample, as before
    expect(o.breakdownCount).toBeUndefined();
    expect(o.timerClosed).toBeUndefined();
    expect(interpolateSnapshot(a, b, 0, createEmptySnapshot()).spark.energyDelivered).toBe(40e-3);
    // growing energy still blends exactly as before
    const c = snap(3, { phase: 'arc', energyDelivered: 10e-3 });
    expect(interpolateSnapshot(b, c, 0.3, createEmptySnapshot()).spark.energyDelivered).toBe(0 + (10e-3 - 0) * 0.3);
  });

  it("applies per cylinder, each with its own event", () => {
    const a = createEmptySnapshot(4);
    const b = createEmptySnapshot(4);
    a.t = 1;
    b.t = 2;
    a.cylinders![2].spark = spark(endOfEvent);
    b.cylinders![2].spark = spark(newEvent);
    a.cylinders![1].spark = spark({ ...endOfEvent, breakdownCount: 4 });
    b.cylinders![1].spark = spark({ ...endOfEvent, breakdownCount: 5, energyDelivered: 32e-3 });
    const o = interpolateSnapshot(a, b, 0.5, createEmptySnapshot(4));
    expect(o.cylinders![2].spark.breakdownCount).toBe(0);
    expect(o.cylinders![2].spark.timerClosed).toBe(true);
    expect(o.cylinders![2].spark.energyDelivered).toBe(0);
    expect(o.cylinders![1].spark.breakdownCount).toBe(4);
    expect(o.cylinders![1].spark.energyDelivered).toBeCloseTo(31e-3, 15);
    expect(createEmptyCylinderSnapshot(0).spark.breakdownCount).toBeUndefined();
  });
});

// ---------------------------------------------------------------------------------------------------
// Playback of a trembler-like stream (shape only, not physics): one ignition event per 720° at
// 1000 rpm; timer make → 3.5 ms coil build-up → 10 breakdowns 2 ms apart, each followed by a 0.4 ms glow
// delivering 2 mJ; contact open after 25 ms; counters held until the next make resets them.
// Raw samples every 20 µs while the contact is closed, every 0.5° (83 µs) otherwise.
// ---------------------------------------------------------------------------------------------------

const RPM = 1000;
const EVENT_PERIOD = 120 / RPM; // s per 720°
const MAKE = 0.01; // s into each period
const BUILD = 3.5e-3;
const BUZZ = 2e-3;
const GLOW = 0.4e-3;
const SPARKS = 10;
const E_SPARK = 2e-3;
const CONTACT = 25e-3;

function rawStream(periods: number): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  const coarse = 0.5 / (6 * RPM);
  const fine = 20e-6;
  let t = 0;
  while (t < periods * EVENT_PERIOD) {
    const tp = t - Math.floor(t / EVENT_PERIOD) * EVENT_PERIOD;
    const tau = tp - MAKE; // since this period's make
    const ev = Math.floor(t / EVENT_PERIOD);
    let count = 0;
    let energy = 0;
    let phase: SparkPhase = ev > 0 || tau >= 0 ? 'done' : 'off';
    const closed = tau >= 0 && tau < CONTACT;
    if (tau < 0 && ev > 0) {
      // previous event's counters held until this period's make
      count = SPARKS;
      energy = E_SPARK * SPARKS;
    } else if (tau >= 0) {
      const since = tau - BUILD;
      const n = since >= 0 ? Math.min(SPARKS, Math.floor(since / BUZZ) + 1) : 0;
      count = n;
      const ph = Math.max(0, since - (n - 1) * BUZZ);
      energy = n > 0 ? E_SPARK * (n - 1 + Math.min(1, ph / GLOW)) : 0; // monotone within the event
      if (closed) phase = n === 0 || ph >= GLOW ? 'charging' : ph < fine ? 'breakdown' : 'glow';
    }
    const s = createEmptySnapshot();
    s.t = t;
    s.rpm = RPM;
    s.cycle = ev;
    s.thetaDeg = (t - ev * EVENT_PERIOD) * 6 * RPM - 360;
    s.spark = spark({
      phase,
      breakdownCount: count,
      energyDelivered: energy,
      timerClosed: closed,
      pointsOpen: closed && phase !== 'charging',
      firstSparkDeg: count > 0 ? -25 : NaN,
    });
    out.push(s);
    // land samples exactly on the breakdown instants and glow ends while the contact is closed
    t = closed || (tau < 0 && tau > -coarse) ? t + fine : t + coarse;
    t = Math.round(t / fine) * fine;
  }
  return out;
}

/** The consumers' differencing rule (state.ts): a decrease is a new event, counted from zero. */
function inc(prev: number, next: number): number {
  return next >= prev ? next - prev : Math.max(next, 0);
}

describe('playback across event restarts (trembler stream)', () => {
  const raw = rawStream(3);
  const buf = new SnapshotBuffer();
  buf.pushBatch(raw);

  for (const [label, frameDt] of [
    ['slow motion: many frames inside each restart bracket', 7.3e-6],
    ['about one frame per sample', 61e-6],
    ['real time: several breakdowns per frame', 1 / 60],
  ] as const) {
    it(`sees every breakdown once and no invented energy (${label})`, () => {
      const out = createEmptySnapshot();
      const vis = new CombustionVisualState(MODEL_T);
      const t0 = raw[0].t;
      const tEnd = raw[raw.length - 1].t;
      let first = true;
      let lastN = 0;
      let lastE = 0;
      let flashes = 0;
      let energy = 0;
      let latchedN = 0;
      let latchedE = 0;
      for (let t = t0 + 0.37 * frameDt; ; t += frameDt) {
        const tf = Math.min(t, tEnd);
        const s = buf.sample(tf, out)!;
        const sp = s.spark;
        expect(Number.isInteger(sp.breakdownCount)).toBe(true);
        // the discrete spark record is one raw sample's (a or b of the bracket), never a mixture
        const i = buf.indexAtOrBefore(tf);
        const cands = [buf.at(i), buf.at(Math.min(i + 1, buf.length - 1))].map((r) => r.spark);
        expect(cands.some((c) => c.breakdownCount === sp.breakdownCount && c.timerClosed === sp.timerClosed && c.phase === sp.phase)).toBe(true);
        vis.update(s, frameDt, 1);
        if (first) {
          latchedN = lastN = sp.breakdownCount!;
          latchedE = lastE = sp.energyDelivered;
          first = false;
        } else {
          flashes += inc(lastN, sp.breakdownCount!);
          energy += inc(lastE, sp.energyDelivered);
          lastN = sp.breakdownCount!;
          lastE = sp.energyDelivered;
        }
        if (tf >= tEnd) break;
      }
      // raw-stream oracle from the latched first frame on (same rule, every raw sample)
      let n0 = latchedN;
      let e0 = latchedE;
      let rawFlashes = 0;
      let rawEnergy = 0;
      const tFirst = t0 + 0.37 * frameDt;
      for (const r of raw) {
        if (r.t <= tFirst) continue;
        rawFlashes += inc(n0, r.spark.breakdownCount!);
        rawEnergy += inc(e0, r.spark.energyDelivered);
        n0 = r.spark.breakdownCount!;
        e0 = r.spark.energyDelivered;
      }
      // the stream really restarts its counters (the generator itself is monotone within an event)
      expect(raw.filter((r, j) => j > 0 && r.spark.energyDelivered < raw[j - 1].spark.energyDelivered)).toHaveLength(2);
      expect(rawFlashes).toBeGreaterThanOrEqual(2 * SPARKS);
      expect(flashes).toBe(rawFlashes);
      expect(vis.breakdownsSeen).toBe(rawFlashes);
      expect(energy).toBeCloseTo(rawEnergy, 12);
    });
  }
});
