/**
 * DOM helper for the 'temperature' mode legend (for the UI layer).
 * Data-only variant: `temperatureLegend()` in ./colour/colormap.
 */
import { temperatureLegend } from './colour/colormap';
import { VIS_TEMPERATURE_RANGE } from './constants';

export interface TemperatureLegendElement {
  element: HTMLElement;
  /** Update the range (e.g. after CombustionVisuals.setTemperatureRange). */
  setRange(range: readonly [number, number]): void;
}

/**
 * Build a compact horizontal colour bar with ticks, styled inline (no CSS
 * dependencies). Width follows the container; about 36 px tall.
 */
export function createTemperatureLegendElement(
  range: readonly [number, number] = VIS_TEMPERATURE_RANGE,
  doc: Document = document,
): TemperatureLegendElement {
  const root = doc.createElement('div');
  root.style.cssText = 'font: 11px/1.2 system-ui, sans-serif; color: inherit; min-width: 160px; user-select: none;';
  const title = doc.createElement('div');
  title.style.cssText = 'margin-bottom: 3px; opacity: 0.85;';
  const bar = doc.createElement('div');
  bar.style.cssText = 'height: 10px; border-radius: 2px; box-shadow: inset 0 0 0 1px rgba(127,127,127,0.35);';
  const ticks = doc.createElement('div');
  ticks.style.cssText = 'position: relative; height: 14px; margin-top: 2px;';
  root.append(title, bar, ticks);

  const setRange = (r: readonly [number, number]): void => {
    const L = temperatureLegend(r);
    title.textContent = L.title;
    bar.style.background = L.cssGradient;
    ticks.replaceChildren();
    for (const t of L.ticks) {
      const s = doc.createElement('span');
      s.textContent = t.label;
      const tx = t.position <= 0.02 ? '0' : t.position >= 0.98 ? '-100%' : '-50%';
      s.style.cssText = `position: absolute; left: ${(t.position * 100).toFixed(2)}%; transform: translateX(${tx}); opacity: 0.8;`;
      ticks.append(s);
    }
  };
  setRange(range);
  return { element: root, setRange };
}
