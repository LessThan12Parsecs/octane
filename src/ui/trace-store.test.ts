import { describe, expect, it } from 'vitest';
import { makeRun, makeSnapshot } from './test-helpers';
import { CH, CycleTrace, CycleTraceStore, emptySelection, fillSeries, N_CH, polytropicIndex, sampleChannels } from './trace-store';

describe('sampleChannels', () => {
  it('converts to display units', () => {
    const s = makeSnapshot({
      pressure: 42e5,
      volume: 250e-6,
      temperatureBurned: 0,
      temperatureUnburned: 700,
      heatReleaseRate: 36e3,
      rpm: 600,
      intakeLift: 0.005,
      exhaustMassFlow: 0.012,
    });
    const out = sampleChannels(s, new Float64Array(N_CH));
    expect(out[CH.p]).toBeCloseTo(42, 12);
    expect(out[CH.V]).toBeCloseTo(250, 9);
    expect(out[CH.Tb]).toBeNaN(); // no burned zone
    expect(out[CH.Tu]).toBe(700);
    expect(out[CH.hrr]).toBeCloseTo(10, 12); // J/°CA
    expect(out[CH.liftIn]).toBeCloseTo(5, 12);
    expect(out[CH.mEx]).toBeCloseTo(12, 12);
  });
});

describe('CycleTraceStore.ingest', () => {
  it('ingests each snapshot once and never beyond the playhead', () => {
    const run = makeRun({ n: 100, dDeg: 1 });
    const st = new CycleTraceStore(4);
    const tPlay = run[49].t;
    expect(st.ingest(run, tPlay)).toBe(50);
    expect(st.ingest(run, tPlay)).toBe(0); // duplicates ignored
    expect(st.newest!.n).toBe(50);
    expect(st.ingest(run.slice(40, 70), run[69].t)).toBe(20);
    expect(st.newest!.n).toBe(70);
    expect(st.lastT).toBe(run[69].t);
  });

  it('starts a new trace on cycle change and on a θ wrap', () => {
    const st = new CycleTraceStore(4);
    const run = makeRun({ n: 800, dDeg: 1, theta0: -100 }); // wraps at +360 → cycle 1
    st.ingest(run, Infinity);
    expect(st.traces.length).toBe(2);
    expect(st.traces[0].cycle).toBe(0);
    expect(st.traces[0].lastTheta).toBe(359);
    expect(st.traces[1].cycle).toBe(1);
    expect(st.traces[1].cols[CH.theta][0]).toBe(-360);

    // θ wrap without a cycle-counter change still splits
    const st2 = new CycleTraceStore(4);
    const a = makeSnapshot({ t: 1, cycle: 5, thetaDeg: 350 });
    const b = makeSnapshot({ t: 2, cycle: 5, thetaDeg: -355 });
    st2.ingest([a, b], Infinity);
    expect(st2.traces.length).toBe(2);
  });

  it('keeps θ monotonic within a trace', () => {
    const st = new CycleTraceStore(4);
    const snaps = [10, 11, 10.5, 12].map((th, i) => makeSnapshot({ t: i + 1, thetaDeg: th }));
    st.ingest(snaps, Infinity);
    expect(Array.from(st.newest!.col(CH.theta))).toEqual([10, 11, 12]);
  });

  it('evicts the oldest cycle and reuses buffers', () => {
    const st = new CycleTraceStore(3, 16);
    st.ingest(makeRun({ n: 720 * 5, dDeg: 1 }), Infinity);
    expect(st.traces.map((t) => t.cycle)).toEqual([2, 3, 4]);
    expect(st.traces[0].n).toBe(720);
    expect(st.traces[0].capacity).toBeGreaterThanOrEqual(720);
  });

  it('records knock', () => {
    const st = new CycleTraceStore(3);
    const run = makeRun({
      n: 40,
      dDeg: 1,
      theta0: -20,
      fn: (s) => {
        s.knock.autoignited = s.thetaDeg >= 8;
      },
    });
    st.ingest(run, Infinity);
    expect(st.newest!.knocked).toBe(true);
    expect(st.newest!.knockDeg).toBe(8);
  });
});

describe('CycleTraceStore.select', () => {
  const st = new CycleTraceStore(6);
  const run = makeRun({ n: 720 * 3 + 100, dDeg: 1 });
  st.ingest(run, Infinity);

  it('returns the cycle under the playhead, truncated, plus older cycles newest-first', () => {
    const sel = st.select(run[720 * 3 + 50].t, 5);
    expect(sel.current!.cycle).toBe(3);
    expect(sel.currentCount).toBe(51);
    expect(sel.ghosts.map((g) => g.cycle)).toEqual([2, 1, 0]);
  });

  it('follows the playhead backwards (stepping back while paused)', () => {
    const out = emptySelection();
    st.select(run[720 * 2 + 10].t, 1, out);
    expect(out.current!.cycle).toBe(2);
    expect(out.currentCount).toBe(11);
    expect(out.ghosts.map((g) => g.cycle)).toEqual([1]);
  });

  it('is empty before any data', () => {
    expect(new CycleTraceStore().select(1, 3).current).toBeNull();
    expect(st.select(-1, 3).current).toBeNull();
  });
});

describe('CycleTraceStore.checkDiscontinuity', () => {
  it('clears after a reset (time jumps back before the history)', () => {
    const st = new CycleTraceStore(2);
    st.ingest(makeRun({ n: 720 * 4, dDeg: 1 }), Infinity);
    expect(st.checkDiscontinuity(st.lastT - 1e-3, 3, 600)).toBe(false); // small rewind
    expect(st.checkDiscontinuity(0, 0, 600)).toBe(true);
    expect(st.traces.length).toBe(0);
    expect(st.lastT).toBe(-Infinity);
  });

  it('clears on a long rewind even if data still covers it', () => {
    const st = new CycleTraceStore(10);
    const run = makeRun({ n: 720 * 5, dDeg: 1 });
    st.ingest(run, Infinity);
    expect(st.checkDiscontinuity(run[100].t, 0, 600)).toBe(true);
  });
});

describe('chart data helpers', () => {
  it('fillSeries maps non-finite values to null and reuses the array', () => {
    const out: (number | null)[] = [9, 9, 9, 9, 9];
    const src = new Float64Array([1, NaN, 3, Infinity]);
    const r = fillSeries(src, 4, out, 2);
    expect(r).toBe(out);
    expect(out).toEqual([2, null, 6, null]);
  });

  it('polytropicIndex recovers n from p V^n = const', () => {
    const tr = new CycleTrace(64);
    const s = new Float64Array(N_CH);
    const n = 1.32;
    for (let i = 0; i < 50; i++) {
      const V = 600 - i * 10;
      s[CH.theta] = -150 + i * 2;
      s[CH.V] = V;
      s[CH.p] = 1.0 * Math.pow(600 / V, n);
      tr.push(s);
    }
    expect(polytropicIndex(tr, tr.n, -120, -60)).toBeCloseTo(n, 10);
    expect(polytropicIndex(tr, tr.n, 100, 120)).toBeNaN();
  });
});
