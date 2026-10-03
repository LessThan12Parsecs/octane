import { describe, expect, it } from 'vitest';
import { GatedMean } from './cycle-mean';

/** Synthetic single-cylinder MAP: drawn down to 0.4 bar during a 0.08 s "valve-open" window every 0.2 s, ~1.0 bar otherwise. */
const period = 0.2;
const openWindow = 0.08;
const gate = (t: number): boolean => t % period < openWindow;
const map = (t: number): number => (gate(t) ? 0.4 + 0.1 * Math.sin((Math.PI * (t % period)) / openWindow) : 1.0);
const exactOpenMean = 0.4 + 0.1 * (2 / Math.PI);

describe('GatedMean', () => {
  it('reports the mean over the open window, independent of the sampling rate', () => {
    for (const dt of [2e-5, 1e-3]) {
      const f = new GatedMean();
      let m = NaN;
      for (let t = 0; t < 3 * period; t += dt) m = f.update(map(t), t, gate(t));
      expect(m).toBeCloseTo(exactOpenMean, 2);
    }
  });

  it('is NaN until a window has closed, then holds the value while the gate is shut', () => {
    const f = new GatedMean();
    expect(Number.isNaN(f.update(0.5, 0.01, true))).toBe(true);
    expect(Number.isNaN(f.update(0.5, 0.05, true))).toBe(true);
    const closed = f.update(1.0, 0.09, false);
    expect(Number.isFinite(closed)).toBe(true);
    expect(f.update(1.0, 0.15, false)).toBe(closed); // held while shut
  });

  it('abandons the window in progress when time steps backwards', () => {
    const f = new GatedMean();
    for (let t = 0; t < period; t += 1e-3) f.update(map(t), t, gate(t));
    const held = f.value;
    f.update(0.9, 0.01, true); // stepped back into an open window
    expect(f.value).toBe(held);
    f.reset();
    expect(Number.isNaN(f.value)).toBe(true);
  });
});
