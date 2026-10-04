/** Conductor → in-cylinder visuals: the mechanism's cut-away region per cylinder frame (fakes, no three.js). */
import { describe, expect, it } from 'vitest';
import type { CycleSummary, EngineSnapshot } from '../physics/core/snapshot';
import { CFR_F1 } from '../physics/engines/cfr';
import { MODEL_T } from '../physics/engines/model-t';
import type { CutPlanes } from '../render/combustion/index';
import { Conductor, GasFanOut, type CylinderGasPort, type GasPort, type MechanismPort, type SimPort } from './conductor';

const sim: SimPort = {
  playbackTime: NaN,
  advance: () => null,
  stepDegrees: () => null,
  recentSnapshots: (): EngineSnapshot[] => [],
  setOperatingPoint: () => undefined,
  reset: () => undefined,
  onCycle: (_cb: (c: CycleSummary) => void) => undefined,
};

/** Cylinder i's region: one plane z > i (distinct per cylinder). */
const regionOf = (i: number): CutPlanes => [[0, 0, 1, i]];

class Mechanism implements MechanismPort {
  compressionRatio = 4;
  cut = true;
  setCompressionRatio(): void {}
  setCutaway(on: boolean): void {
    this.cut = on;
  }
  update(): void {}
}

class CutMechanism extends Mechanism {
  cutRegion(i: number): CutPlanes | null {
    return this.cut ? regionOf(i) : null;
  }
}

class CylGas implements CylinderGasPort {
  readonly regions: (CutPlanes | null)[] = [];
  constructor(readonly cylinder: number) {}
  update(): void {}
  setMode(): void {}
  setCutRegion(planes: CutPlanes | null): void {
    this.regions.push(planes);
  }
}

class PlainGas implements GasPort {
  update(): void {}
  setMode(): void {}
}

describe('Conductor cut-away region', () => {
  it('hands each cylinder its region at construction and on every cutaway change', () => {
    const gas = [0, 1, 2, 3].map((i) => new CylGas(i));
    const c = new Conductor(MODEL_T, sim, new CutMechanism(), new GasFanOut(gas));
    gas.forEach((g, i) => expect(g.regions).toEqual([regionOf(i)]));
    c.handleView({ cutaway: false, mode: 'physical' });
    gas.forEach((g, i) => expect(g.regions).toEqual([regionOf(i), null]));
    c.handleView({ cutaway: true, mode: 'temperature' });
    gas.forEach((g, i) => expect(g.regions).toEqual([regionOf(i), null, regionOf(i)]));
  });

  it('forwards through nested fan-outs and skips ports that cannot show a cut', () => {
    const inner = [new CylGas(2), new CylGas(3)];
    const a = new CylGas(0);
    new Conductor(MODEL_T, sim, new CutMechanism(), new GasFanOut([a, new PlainGas(), new GasFanOut(inner)]));
    expect(a.regions).toEqual([regionOf(0)]);
    expect(inner[0].regions).toEqual([regionOf(2)]);
    expect(inner[1].regions).toEqual([regionOf(3)]);
  });

  it('never sets a region when the mechanism reports none (the CFR)', () => {
    const g = new CylGas(0);
    const c = new Conductor(CFR_F1, sim, new Mechanism(), new GasFanOut([g]));
    c.handleView({ cutaway: false, mode: 'physical' });
    c.handleView({ cutaway: true, mode: 'physical' });
    c.setSectionCylinder(0); // no movable section: nothing happens
    expect(g.regions).toEqual([]);
  });

  it('moving the section to the focus cylinder re-pushes every cylinder’s moved region (review: section follows the focus)', () => {
    /** Section through cylinder k: cylinder i's region is one plane z > 10·k + i. */
    class MovableSection extends Mechanism {
      section = 0;
      readonly moves: number[] = [];
      setSectionCylinder(k: number): void {
        this.moves.push(k);
        this.section = k;
      }
      cutRegion(i: number): CutPlanes | null {
        return this.cut ? [[0, 0, 1, 10 * this.section + i]] : null;
      }
    }
    const mech = new MovableSection();
    const gas = [0, 1, 2, 3].map((i) => new CylGas(i));
    const c = new Conductor(MODEL_T, sim, mech, new GasFanOut(gas));
    c.setSectionCylinder(2);
    expect(mech.moves).toEqual([2]);
    gas.forEach((g, i) => expect(g.regions).toEqual([[[0, 0, 1, i]], [[0, 0, 1, 20 + i]]]));
    // cutaway off: the move still reaches the mechanism; the visuals get null (nothing cut)
    c.handleView({ cutaway: false, mode: 'physical' });
    c.setSectionCylinder(3);
    gas.forEach((g) => expect(g.regions.slice(2)).toEqual([null, null]));
    c.handleView({ cutaway: true, mode: 'physical' });
    gas.forEach((g, i) => expect(g.regions.at(-1)).toEqual([[0, 0, 1, 30 + i]]));
    // a mechanism that throws re-pushes nothing
    mech.setSectionCylinder = () => {
      throw new Error('no section');
    };
    const before = gas.map((g) => g.regions.length);
    expect(() => c.setSectionCylinder(1)).toThrow('no section');
    expect(gas.map((g) => g.regions.length)).toEqual(before);
  });

  it('GasFanOut.replace swaps one port and rejects other indices', () => {
    const gas = [new CylGas(0), new CylGas(1)];
    const fan = new GasFanOut(gas);
    const b = new CylGas(1);
    fan.replace(1, b);
    expect(fan.ports).toEqual([gas[0], b]);
    expect(gas).toEqual([gas[0], gas[1]]); // the caller's array is not aliased
    expect(() => fan.replace(2, b)).toThrow(RangeError);
    expect(() => fan.replace(-1, b)).toThrow(RangeError);
  });
});
