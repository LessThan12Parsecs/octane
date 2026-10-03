import { describe, expect, it } from 'vitest';
import { emptyExtent, extendExtent, isEmptyExtent, logTicks125, niceCeil, niceFloor, niceLogRange, niceRange } from './chart-math';

describe('extents', () => {
  it('accumulates finite values inside the x window only', () => {
    const xs = [-10, 0, 10, 20];
    const ys = [5, NaN, -3, 100];
    const e = extendExtent(xs, ys, 4, -10, 10, emptyExtent());
    expect(e).toEqual([-3, 5]);
    expect(isEmptyExtent(emptyExtent())).toBe(true);
    expect(extendExtent(null, [2, 4], 2, 0, 0, emptyExtent(), 10)).toEqual([20, 40]);
  });
});

describe('nice numbers', () => {
  it('niceCeil / niceFloor bracket the value', () => {
    expect(niceCeil(37)).toBe(40);
    expect(niceCeil(0.0123)).toBeCloseTo(0.015, 15);
    expect(niceCeil(100)).toBe(100);
    expect(niceFloor(37)).toBe(30);
    expect(niceFloor(0.95)).toBeCloseTo(0.8, 15);
    expect(niceCeil(-37)).toBe(-30);
    expect(niceFloor(-37)).toBe(-40);
  });

  it('niceRange honours fixed and included bounds', () => {
    expect(niceRange([0.1, 0.9], { min: 0, max: 1 })).toEqual([0, 1]);
    const r = niceRange([1, 38.2], { includeMin: 0 });
    expect(r[0]).toBe(0);
    expect(r[1]).toBeGreaterThanOrEqual(38.2 * 1.05 - 1e-9);
    expect(r[1]).toBeLessThan(45);
    // Always contains the data.
    const r2 = niceRange([-12.3, 47.1], {});
    expect(r2[0]).toBeLessThanOrEqual(-12.3);
    expect(r2[1]).toBeGreaterThanOrEqual(47.1);
  });

  it('niceRange is stable for small data changes (no per-frame jitter)', () => {
    const a = niceRange([0, 38.2], { includeMin: 0 });
    const b = niceRange([0, 38.9], { includeMin: 0 });
    expect(a).toEqual(b);
  });

  it('niceRange enforces a minimum span and handles empty extents', () => {
    const r = niceRange([5, 5], { minSpan: 2 });
    expect(r[1] - r[0]).toBeGreaterThanOrEqual(2);
    expect(r[0]).toBeLessThanOrEqual(5);
    expect(r[1]).toBeGreaterThanOrEqual(5);
    const e = niceRange(emptyExtent(), { includeMin: 0, minSpan: 1 });
    expect(e[0]).toBe(0);
    expect(e[1]).toBeGreaterThan(0);
  });
});

describe('log axes', () => {
  it('niceLogRange brackets the data', () => {
    const r = niceLogRange([0.9, 42], 0.05);
    expect(r[0]).toBeLessThan(0.9);
    expect(r[1]).toBeGreaterThan(42);
    expect(niceLogRange([-1, 0], 0.1)[0]).toBeGreaterThan(0);
  });

  it('logTicks125 yields 1-2-5 ticks inside the range', () => {
    expect(logTicks125(0.5, 60)).toEqual([0.5, 1, 2, 5, 10, 20, 50]);
    expect(logTicks125(80, 900)).toEqual([100, 150, 200, 300, 500, 700]);
    expect(logTicks125(-1, 10)).toEqual([]);
  });
});
