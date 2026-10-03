import { describe, expect, it } from 'vitest';
import { Float64Column, lowerBound, RingBuffer, upperBound } from './ring-buffer';

describe('RingBuffer', () => {
  it('keeps the newest `capacity` items in order and reports evictions', () => {
    const rb = new RingBuffer<number>(3);
    expect(rb.push(1)).toBeUndefined();
    rb.push(2);
    rb.push(3);
    expect(rb.push(4)).toBe(1);
    expect(rb.length).toBe(3);
    expect(rb.toArray()).toEqual([2, 3, 4]);
    expect([...rb]).toEqual([2, 3, 4]);
    expect(rb.get(0)).toBe(2);
    expect(rb.latest()).toBe(4);
    expect(rb.latest(2)).toBe(2);
    expect(rb.latest(3)).toBeUndefined();
    expect(rb.get(-1)).toBeUndefined();
    rb.clear();
    expect(rb.length).toBe(0);
    expect(rb.latest()).toBeUndefined();
    rb.push(9);
    expect(rb.toArray()).toEqual([9]);
  });

  it('wraps many times', () => {
    const rb = new RingBuffer<number>(5);
    for (let i = 0; i < 103; i++) rb.push(i);
    expect(rb.toArray()).toEqual([98, 99, 100, 101, 102]);
  });

  it('rejects invalid capacities', () => {
    expect(() => new RingBuffer(0)).toThrow(RangeError);
    expect(() => new RingBuffer(1.5)).toThrow(RangeError);
  });
});

describe('Float64Column', () => {
  it('grows while preserving data', () => {
    const c = new Float64Column(2);
    for (let i = 0; i < 100; i++) c.push(i * 0.5);
    expect(c.length).toBe(100);
    expect(c.data.length).toBeGreaterThanOrEqual(100);
    expect(c.view()[99]).toBe(49.5);
    expect(c.view(10).length).toBe(10);
    c.clear();
    expect(c.length).toBe(0);
    expect(c.view().length).toBe(0);
  });
});

describe('binary search', () => {
  const a = new Float64Array([1, 2, 2, 3, 5]);
  it('upperBound: first index with value > x', () => {
    expect(upperBound(a, a.length, 0)).toBe(0);
    expect(upperBound(a, a.length, 2)).toBe(3);
    expect(upperBound(a, a.length, 4)).toBe(4);
    expect(upperBound(a, a.length, 9)).toBe(5);
    expect(upperBound(a, 3, 9)).toBe(3);
  });
  it('lowerBound: first index with value ≥ x', () => {
    expect(lowerBound(a, a.length, 2)).toBe(1);
    expect(lowerBound(a, a.length, 5)).toBe(4);
    expect(lowerBound(a, a.length, 6)).toBe(5);
  });
});
