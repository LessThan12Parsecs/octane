/**
 * SimClient / interpolateSnapshot with multi-cylinder streams (snapshot.cylinders[]) and the optional
 * snapshot fields (friction/load torque, vehicle speed, magneto EMF, trembler spark train).
 */
import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, CylinderSnapshot, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import { MODEL_T, MODEL_T_CRUISE } from '../physics/engines/model-t';
import type { FromWorker, ToWorker } from '../worker/protocol';
import {
  copySnapshot,
  createEmptyCylinderSnapshot,
  createEmptySnapshot,
  interpolateCylinder,
  interpolateSnapshot,
  SimClient,
  type WorkerLike,
} from './sim-client';

// ---------------------------------------------------------------------------
// A snapshot with EVERY field present (optional ones included). The types make the compiler reject
// this literal when the contract gains a field, so the leaf walk below always covers the contract.
// ---------------------------------------------------------------------------

type FullSpark = Required<EngineSnapshot['spark']>;
type FullCylinder = Omit<CylinderSnapshot, 'spark'> & { spark: FullSpark };
type FullSnapshot = Required<Omit<EngineSnapshot, 'spark' | 'cylinders'>> & { spark: FullSpark; cylinders: FullCylinder[] };

function fullSpark(): FullSpark {
  return {
    phase: 'arc',
    primaryCurrent: 0,
    secondaryVoltage: 0,
    secondaryCurrent: 0,
    energyDelivered: 0,
    breakdownVoltage: 0,
    breakdownCount: 0,
    reignitionCount: 0,
    pointsOpen: true,
    timerClosed: true,
    firstSparkDeg: 0,
    primaryVoltage: 0,
  };
}

function fullCylinder(index: number): FullCylinder {
  return {
    index,
    thetaDeg: 0,
    cycle: 0,
    gasTorque: 0,
    pistonDisplacement: 0,
    clearanceHeight: 0,
    rodAngle: 0,
    intakeLift: 0,
    exhaustLift: 0,
    phase: 'combustion',
    volume: 0,
    pressure: 0,
    temperatureMean: 0,
    temperatureUnburned: 0,
    temperatureBurned: 0,
    massFractionBurned: 0,
    mass: 0,
    heatReleaseRate: 0,
    heatLossRate: 0,
    flame: { stage: 'kernel', radius: 0, center: [0, 0, 0], area: 0, laminarSpeed: 0, turbulentSpeed: 0, turbulenceIntensity: 0 },
    spark: fullSpark(),
    intakeMassFlow: 0,
    exhaustMassFlow: 0,
    knock: { integral: 0, autoignited: false, oscillation: 0 },
    burnedComposition: { CO2: 0, H2O: 0, CO: 0, O2: 0, H2: 0, OH: 0, H: 0, O: 0, NO: 0, N2: 0 },
  };
}

function fullSnapshot(): FullSnapshot {
  const c = fullCylinder(0);
  const { index: _i, ...cyl } = c;
  void _i;
  return {
    ...cyl,
    t: 0,
    rpm: 0,
    intakeManifoldPressure: 0,
    exhaustManifoldPressure: 0,
    netTorque: 0,
    frictionTorque: 0,
    loadTorque: 0,
    vehicleSpeed: 0,
    firingCylinder: 0,
    magnetoEmf: 0,
    cylinders: [0, 1, 2, 3].map(fullCylinder),
  };
}

/** Leaves whose numeric value is discrete (taken from one sample, never blended). */
const DISCRETE_NUMERIC = /(^|\.)(cycle|index|breakdownCount|reignitionCount|firingCylinder|firstSparkDeg)$/;
/** Crank angles (wrap-aware lerp; the test keeps them away from the wrap). */
const ANGLE = /(^|\.)thetaDeg$/;

type Leaf = { path: string; get: (o: unknown) => unknown; set: (o: unknown, v: unknown) => void };

/** Every leaf (number / boolean / string) of a plain object tree, with accessors by path. */
function leaves(o: unknown, prefix = '', path: (string | number)[] = [], out: Leaf[] = []): Leaf[] {
  if (o !== null && typeof o === 'object') {
    for (const [k, v] of Object.entries(o as Record<string, unknown>)) {
      const key = Array.isArray(o) ? Number(k) : k;
      leaves(v, prefix ? `${prefix}.${k}` : k, [...path, key], out);
    }
    return out;
  }
  const p = path.slice();
  const walk = (x: unknown): Record<string | number, unknown> => {
    let cur = x as Record<string | number, unknown>;
    for (let i = 0; i < p.length - 1; i++) cur = cur[p[i]] as Record<string | number, unknown>;
    return cur;
  };
  out.push({ path: prefix, get: (x) => walk(x)[p[p.length - 1]], set: (x, v) => (walk(x)[p[p.length - 1]] = v) });
  return out;
}

describe('interpolateSnapshot: every numeric leaf', () => {
  // Distinct values for every numeric leaf of a and b (no wrap anywhere).
  const a = fullSnapshot();
  const b = fullSnapshot();
  const all = leaves(a);
  let k = 1;
  for (const leaf of all) {
    if (typeof leaf.get(a) !== 'number') continue;
    if (ANGLE.test(leaf.path)) {
      leaf.set(a, 10 + k);
      leaf.set(b, 12 + k);
    } else if (/(^|\.)index$/.test(leaf.path)) {
      leaf.set(b, leaf.get(a)); // identity of the cylinder
    } else {
      leaf.set(a, k);
      leaf.set(b, 1000 + 3 * k);
    }
    k++;
  }
  // Discrete leaves of b differ from a so we can tell which sample they came from.
  for (const leaf of leaves(b)) {
    const v = leaf.get(b);
    if (typeof v === 'boolean') leaf.set(b, !v);
    if (typeof v === 'string') leaf.set(b, leaf.path.endsWith('phase') && leaf.path.includes('spark') ? 'glow' : leaf.path.endsWith('stage') ? 'turbulent' : 'expansion');
  }
  const u = 0.25;

  for (const [label, target] of [
    ['preallocated target', () => createEmptySnapshot(4)],
    ['bare target (cylinders allocated on first use)', () => createEmptySnapshot()],
  ] as const) {
    it(`interpolates or carries every leaf (${label})`, () => {
      const out = interpolateSnapshot(a as EngineSnapshot, b as EngineSnapshot, u, target());
      const outLeaves = leaves(out);
      expect(outLeaves.map((l) => l.path).sort()).toEqual(all.map((l) => l.path).sort());
      for (const leaf of all) {
        const va = leaf.get(a);
        const vb = leaf.get(b);
        const vo = leaf.get(out);
        if (typeof va === 'number') {
          if (DISCRETE_NUMERIC.test(leaf.path)) expect(vo, leaf.path).toBe(va);
          else expect(vo as number, leaf.path).toBeCloseTo(va + ((vb as number) - va) * u, 9);
        } else {
          expect(vo, leaf.path).toBe(va); // discrete: earlier sample (no wrap)
        }
      }
    });
  }

  it('copySnapshot reproduces a full multi-cylinder snapshot without aliasing', () => {
    const out = copySnapshot(a as EngineSnapshot, createEmptySnapshot(4));
    expect(out).toEqual(a);
    expect(out.cylinders![2]).not.toBe(a.cylinders[2]);
    expect(out.cylinders![2].flame.center).not.toBe(a.cylinders[2].flame.center);
    expect(out.cylinders![2].spark).not.toBe(a.cylinders[2].spark);
  });
});

describe('interpolateSnapshot: multi-cylinder semantics', () => {
  it("chooses each cylinder's discrete fields by that cylinder's own wrap", () => {
    const a = createEmptySnapshot(4);
    const b = createEmptySnapshot(4);
    a.t = 1;
    b.t = 2;
    a.cycle = b.cycle = 7;
    a.thetaDeg = 179.5; // engine angle: no wrap
    b.thetaDeg = 180.5; // still engine cycle 7
    // Cylinder 2 (offset 180) is crossing its local wrap 359.5 → −359.5.
    const ca = a.cylinders![1];
    const cb = b.cylinders![1];
    ca.thetaDeg = 359.5;
    ca.cycle = 6;
    ca.phase = 'gas-exchange';
    ca.flame.stage = 'done';
    cb.thetaDeg = -359.5;
    cb.cycle = 7;
    cb.flame.stage = 'none';
    cb.spark.timerClosed = true;
    const before = interpolateSnapshot(a, b, 0.25, createEmptySnapshot(4));
    expect(before.cylinders![1].thetaDeg).toBeCloseTo(359.75, 12);
    expect(before.cylinders![1].cycle).toBe(6);
    expect(before.cylinders![1].flame.stage).toBe('done');
    const after = interpolateSnapshot(a, b, 0.75, createEmptySnapshot(4));
    expect(after.cylinders![1].thetaDeg).toBeCloseTo(-359.75, 12);
    expect(after.cylinders![1].cycle).toBe(7);
    expect(after.cylinders![1].flame.stage).toBe('none');
    expect(after.cylinders![1].spark.timerClosed).toBe(true);
    // The engine angle did not wrap: the engine cycle stays.
    expect(after.cycle).toBe(7);
    expect(after.thetaDeg).toBeCloseTo(180.25, 12);
  });

  it('reuses the preallocated cylinder objects (no allocation per call)', () => {
    const a = createEmptySnapshot(4);
    const b = createEmptySnapshot(4);
    b.t = 1;
    const out = createEmptySnapshot(4);
    const before = out.cylinders!.slice();
    const sparks = out.cylinders!.map((c) => c.spark);
    for (let i = 0; i < 3; i++) interpolateSnapshot(a, b, i / 3, out);
    expect(out.cylinders).toHaveLength(4);
    out.cylinders!.forEach((c, i) => {
      expect(c).toBe(before[i]);
      expect(c.spark).toBe(sparks[i]);
      expect(c.index).toBe(i);
    });
  });

  it('leaves optional fields undefined for a single-cylinder stream (CFR shape unchanged)', () => {
    const a = createEmptySnapshot();
    const b = createEmptySnapshot();
    b.t = 1;
    b.pressure = 2e5;
    const out = interpolateSnapshot(a, b, 0.5, createEmptySnapshot());
    expect(out.cylinders).toBeUndefined();
    expect(out.frictionTorque).toBeUndefined();
    expect(out.loadTorque).toBeUndefined();
    expect(out.vehicleSpeed).toBeUndefined();
    expect(out.magnetoEmf).toBeUndefined();
    expect(out.firingCylinder).toBeUndefined();
    expect(out.spark.breakdownCount).toBeUndefined();
    expect(out.spark.firstSparkDeg).toBeUndefined();
    expect(out.pressure).toBe(1e5);
  });

  it('interpolateCylinder clamps alpha like interpolateSnapshot', () => {
    const a = createEmptyCylinderSnapshot(2);
    const b = createEmptyCylinderSnapshot(2);
    a.pressure = 1;
    b.pressure = 3;
    expect(interpolateCylinder(a, b, 2, createEmptyCylinderSnapshot(0)).pressure).toBe(3);
    expect(interpolateCylinder(a, b, -1, createEmptyCylinderSnapshot(0)).index).toBe(2);
  });
});

// ---------------------------------------------------------------------------
// SimClient with a scripted worker
// ---------------------------------------------------------------------------

const OFFSETS = MODEL_T.layout.firingOffsetDeg;

const wrap = (t: number): number => ((((t + 360) % 720) + 720) % 720) - 360;

/**
 * A 4-cylinder stream at 1° steps from engine θ = −360, cycle 0. Each cylinder counts its own cycles
 * at its local wrap (starting at 0).
 */
function multiRun(n: number, rpm = 1000): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  const cyc = [0, 0, 0, 0];
  const prev = OFFSETS.map((o) => wrap(-360 - o));
  let cycle = 0;
  for (let i = 0; i < n; i++) {
    const thRaw = -360 + i;
    if (i > 0 && wrap(thRaw) === -360) cycle++;
    const s = createEmptySnapshot(4);
    s.t = i / (6 * rpm);
    s.rpm = rpm;
    s.cycle = cycle;
    s.thetaDeg = wrap(thRaw);
    s.cylinders!.forEach((c, k) => {
      const th = wrap(thRaw - OFFSETS[k]);
      if (i > 0 && th < prev[k]) cyc[k]++;
      prev[k] = th;
      c.thetaDeg = th;
      c.cycle = cyc[k];
      c.pressure = 1e5 + k;
    });
    out.push(s);
  }
  return out;
}

class ScriptedWorker implements WorkerLike {
  onmessage: ((ev: MessageEvent<FromWorker>) => void) | null = null;
  onerror: ((ev: ErrorEvent) => void) | null = null;
  readonly sent: ToWorker[] = [];
  postMessage(msg: ToWorker): void {
    this.sent.push(msg);
  }
  terminate(): void {}
  emit(msg: FromWorker): void {
    this.onmessage?.({ data: msg } as MessageEvent<FromWorker>);
  }
}

function modelTClient() {
  const w = new ScriptedWorker();
  const op: OperatingPoint = MODEL_T_CRUISE;
  const client = new SimClient(MODEL_T, op, { snapshotEveryDeg: 1, bufferAheadSeconds: 0.5 }, { createWorker: () => w });
  return { w, client };
}

describe('SimClient with a multi-cylinder stream', () => {
  it('interpolates cylinders[] in playback and keeps the scratch objects', () => {
    const { w, client } = modelTClient();
    w.emit({ type: 'ready' });
    w.emit({ type: 'snapshots', batch: multiRun(200) });
    const s0 = client.advance(0, 1)!;
    expect(s0.cylinders).toHaveLength(4);
    const s1 = client.advance(0.0005 / 6, 1)!; // half a degree at 1000 rpm, real time
    expect(s1.cylinders![3].thetaDeg).toBeCloseTo(wrap(-360 + 0.5 - OFFSETS[3]), 6);
    expect(s1.cylinders![2].pressure).toBe(1e5 + 2);
    const s2 = client.advance(0.0005 / 6, 1)!;
    expect(s2).toBe(s0); // alternating scratch objects carry their cylinder arrays
    expect(s2.cylinders![1].thetaDeg).toBeCloseTo(wrap(-359 - OFFSETS[1]), 6);
  });

  it("delivers every cylinder's summary once, after playback passes that cylinder's cycle end", () => {
    const { w, client } = modelTClient();
    const got: [number, number, number][] = [];
    client.onCycle((c) => got.push([c.cylinder ?? 0, c.cycle, client.playbackTime]));
    w.emit({ type: 'ready' });
    const run = multiRun(1500);
    w.emit({ type: 'snapshots', batch: run });
    // Local wraps (firing offsets 0/180/540/360): cylinder 2 at sample 180, 4 at 360, 3 at 540, 1 at 720.
    // The simulator runs ahead and reports each summary when that cylinder completes its cycle.
    for (const k of [1, 3, 2, 0]) w.emit({ type: 'cycle', summary: { cycle: 0, cylinder: k } as CycleSummary });
    const dt = 1 / (6 * 1000); // 1° per call at real time
    for (let i = 0; i < 1499; i++) {
      client.advance(dt, 1);
      // A repeated summary (e.g. resent) is ignored once delivered.
      if (i === 1000) w.emit({ type: 'cycle', summary: { cycle: 0, cylinder: 0 } as CycleSummary });
    }
    expect(got.map((g) => g[0])).toEqual([1, 3, 2, 0]);
    const wrapSample = [720, 180, 540, 360];
    for (const [k, , t] of got) {
      const i = Math.round(t * 6 * 1000);
      expect(i, `cylinder ${k + 1}`).toBeGreaterThanOrEqual(wrapSample[k]);
      expect(i, `cylinder ${k + 1}`).toBeLessThanOrEqual(wrapSample[k] + 1);
    }
  });

  it('reports worker errors through onError and rejects ready', async () => {
    const { w, client } = modelTClient();
    const errors: string[] = [];
    client.onError((m) => errors.push(m));
    w.emit({ type: 'log', level: 'error', message: 'simulator failed: no multi-cylinder model' });
    w.emit({ type: 'log', level: 'warn', message: 'just a warning' });
    expect(errors).toEqual(['simulator failed: no multi-cylinder model']);
    expect(client.error).toBe('simulator failed: no multi-cylinder model');
    await expect(client.ready).rejects.toThrow('no multi-cylinder model');
    client.dispose();
  });

  it('flags playback starvation when the simulator falls behind', () => {
    const { w, client } = modelTClient();
    w.emit({ type: 'ready' });
    w.emit({ type: 'snapshots', batch: multiRun(10) });
    client.advance(0, 1);
    expect(client.starved).toBe(false);
    client.advance(1 / 6000, 1);
    expect(client.starved).toBe(false);
    client.advance(1, 1); // far past the 9° buffered
    expect(client.starved).toBe(true);
    client.advance(0, 1); // paused / no movement requested
    expect(client.starved).toBe(false);
  });

  it('a single-cylinder client keeps cylinders undefined', () => {
    const w = new ScriptedWorker();
    const client = new SimClient(CFR_F1, MODEL_T_CRUISE, { snapshotEveryDeg: 1, bufferAheadSeconds: 0.5 }, { createWorker: () => w });
    w.emit({ type: 'ready' });
    const a = createEmptySnapshot();
    const b = createEmptySnapshot();
    b.t = 0.01;
    b.thetaDeg = -300;
    w.emit({ type: 'snapshots', batch: [a, b] });
    client.advance(0, 1);
    expect(client.advance(0.005, 1)!.cylinders).toBeUndefined();
  });
});
