import { describe, expect, it } from 'vitest';
import { Pchip } from './pchip';

describe('Pchip (Fritsch–Carlson monotone cubic)', () => {
  it('interpolates the nodes, holds the ends flat, and reproduces a straight line exactly', () => {
    const xs = [0, 0.05, 0.1, 0.15, 0.25];
    const ys = [0.55, 0.54, 0.61, 0.57, 0.42];
    const p = new Pchip(xs, ys);
    xs.forEach((x, i) => expect(p.value(x)).toBeCloseTo(ys[i], 15));
    expect(p.value(-1)).toBe(0.55);
    expect(p.value(1)).toBe(0.42);
    const line = new Pchip([0, 1, 3, 4], [1, 3, 7, 9]);
    for (let t = 0; t <= 4; t += 0.1) expect(line.value(t)).toBeCloseTo(1 + 2 * t, 13);
  });

  it('no overshoot on monotone data; local extrema only at data extrema; C¹ at the nodes', () => {
    const p = new Pchip([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], [0.7, 0.68, 0.64, 0.58, 0.52, 0.46, 0.4]);
    let prev = Infinity;
    for (let t = 0; t <= 0.3; t += 1e-4) {
      const v = p.value(t);
      expect(v).toBeLessThanOrEqual(prev + 1e-15);
      prev = v;
    }
    const q = new Pchip([0, 0.05, 0.1, 0.15, 0.25], [0.55, 0.54, 0.61, 0.57, 0.42]);
    for (let t = 0; t <= 0.25; t += 1e-4) {
      const v = q.value(t);
      expect(v).toBeLessThanOrEqual(0.61 + 1e-15);
      expect(v).toBeGreaterThanOrEqual(0.42 - 1e-15);
    }
    for (const x of [0.05, 0.1, 0.15]) {
      const h = 1e-7;
      const left = (q.value(x) - q.value(x - h)) / h;
      const right = (q.value(x + h) - q.value(x)) / h;
      expect(Math.abs(left - right)).toBeLessThan(1e-4);
    }
    expect(() => new Pchip([0, 0], [1, 2])).toThrow();
  });
});
