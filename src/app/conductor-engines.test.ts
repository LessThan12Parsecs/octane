/** Conductor with a fixed-CR multi-cylinder engine (Model T), the gas fan-out and the controls hook. */
import { describe, expect, it } from 'vitest';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import { MODEL_T } from '../physics/engines/model-t';
import { Conductor, GasFanOut, type GasPort, type MechanismPort, type SimPort } from './conductor';
import { createEmptySnapshot } from './sim-client';

class Sim implements SimPort {
  playbackTime = NaN;
  next: EngineSnapshot | null = null;
  readonly patches: Partial<OperatingPoint>[] = [];
  advance(): EngineSnapshot | null {
    if (this.next) this.playbackTime = this.next.t;
    return this.next;
  }
  stepDegrees(): EngineSnapshot | null {
    return this.next;
  }
  recentSnapshots(): EngineSnapshot[] {
    return this.next ? [this.next] : [];
  }
  setOperatingPoint(patch: Partial<OperatingPoint>): void {
    this.patches.push(patch);
  }
  reset(): void {}
  onCycle(_cb: (c: CycleSummary) => void): void {}
}

class Mechanism implements MechanismPort {
  compressionRatio: number;
  readonly crCalls: number[] = [];
  readonly controls: Partial<OperatingPoint>[] = [];
  updates = 0;
  constructor(cr: number) {
    this.compressionRatio = cr;
  }
  update(): void {
    this.updates++;
  }
  setCompressionRatio(cr: number): void {
    this.crCalls.push(cr);
  }
  setCutaway(): void {}
  setControls(op: Partial<OperatingPoint>): void {
    this.controls.push(op);
  }
}

class Gas implements GasPort {
  readonly updates: number[] = [];
  mode = 'physical';
  update(s: EngineSnapshot): void {
    this.updates.push(s.t);
  }
  setMode(mode: 'physical' | 'temperature'): void {
    this.mode = mode;
  }
}

describe('Conductor with a fixed-CR engine', () => {
  it('never moves the mechanism for the CR, whatever the L-head clearance height implies', () => {
    const sim = new Sim();
    const engine = new Mechanism(MODEL_T.geometry.compressionRatio);
    const crEvents: number[] = [];
    const c = new Conductor(MODEL_T, sim, engine, new Gas(), { onCompressionRatio: (cr) => crEvents.push(cr) });
    for (let k = 0; k < 5; k++) {
      const s = createEmptySnapshot(4);
      s.t = 0.001 * (k + 1);
      s.rpm = 1000;
      // L-head: clearanceHeight is the bore-column depth, not Vc/A — the flat-disc inversion gives nonsense.
      s.clearanceHeight = 0.0254 + 0.01 * k;
      s.pistonDisplacement = 0.002 * k;
      sim.next = s;
      c.frame(16 * k);
    }
    expect(engine.updates).toBe(5);
    expect(engine.crCalls).toEqual([]);
    expect(crEvents).toEqual([]);
  });

  it('a variable-CR engine still follows the data (CFR)', () => {
    const sim = new Sim();
    const engine = new Mechanism(7);
    const c = new Conductor(CFR_F1, sim, engine, new Gas());
    const s = createEmptySnapshot();
    s.t = 0.001;
    s.rpm = 600;
    s.clearanceHeight = 0.02; // h_TDC = 2 cm → some CR ≠ 7
    sim.next = s;
    c.frame(0);
    expect(engine.crCalls).toHaveLength(1);
  });

  it('passes every operating-point change to the mechanism controls (spark lever, hand throttle)', () => {
    const sim = new Sim();
    const engine = new Mechanism(MODEL_T.geometry.compressionRatio);
    const c = new Conductor(MODEL_T, sim, engine, new Gas());
    c.handleOperatingPoint({ sparkAdvanceDeg: 30 });
    c.handleOperatingPoint({ throttle: 0.6, ignitionSource: 'battery' });
    expect(sim.patches).toEqual([{ sparkAdvanceDeg: 30 }, { throttle: 0.6, ignitionSource: 'battery' }]);
    expect(engine.controls).toEqual(sim.patches);
  });

  it('a mechanism without setControls is fine', () => {
    const sim = new Sim();
    const engine: MechanismPort = { compressionRatio: 7, update() {}, setCompressionRatio() {}, setCutaway() {} };
    const c = new Conductor(CFR_F1, sim, engine, new Gas());
    expect(() => c.handleOperatingPoint({ rpm: 900 })).not.toThrow();
    expect(sim.patches).toEqual([{ rpm: 900 }]);
  });
});

describe('GasFanOut', () => {
  it('drives every cylinder visual with the same snapshot and mode', () => {
    const gas = [new Gas(), new Gas(), new Gas(), new Gas()];
    const sim = new Sim();
    const c = new Conductor(MODEL_T, sim, new Mechanism(4), new GasFanOut(gas));
    const s = createEmptySnapshot(4);
    s.t = 0.002;
    sim.next = s;
    c.frame(0);
    c.handleView({ cutaway: true, mode: 'temperature' });
    for (const g of gas) {
      expect(g.updates).toEqual([0.002]);
      expect(g.mode).toBe('temperature');
    }
  });
});
