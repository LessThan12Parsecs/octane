/**
 * DOM overlay inside the 3D viewport (owned by the app, not the UI layer):
 * start-up / error status, camera-framing buttons and the temperature legend.
 */
import type { Framing } from './sync';
import './app.css';

export interface ViewportOverlayOptions {
  onFraming: (f: Framing) => void;
  /** Legend element for the 'temperature' render mode (from createTemperatureLegendElement). */
  legend?: HTMLElement;
}

export class ViewportOverlay {
  readonly el: HTMLElement;
  private readonly status: HTMLElement;
  private readonly statusTitle: HTMLElement;
  private readonly statusDetail: HTMLElement;
  private readonly legendBox: HTMLElement;
  private readonly buttons = new Map<Framing, HTMLButtonElement>();

  constructor(host: HTMLElement, opts: ViewportOverlayOptions) {
    this.el = document.createElement('div');
    this.el.className = 'oct-stage-overlay';

    this.status = document.createElement('div');
    this.status.className = 'oct-stage-status';
    this.status.setAttribute('role', 'status');
    this.statusTitle = document.createElement('span');
    this.statusTitle.className = 'oct-stage-status-title';
    this.statusDetail = document.createElement('span');
    this.statusDetail.className = 'oct-stage-status-detail';
    this.status.append(this.statusTitle, this.statusDetail);

    const dock = document.createElement('div');
    dock.className = 'oct-stage-dock';

    this.legendBox = document.createElement('div');
    this.legendBox.className = 'oct-stage-legend';
    this.legendBox.hidden = true;
    if (opts.legend) this.legendBox.append(opts.legend);

    const bar = document.createElement('div');
    bar.className = 'oct-framing';
    bar.setAttribute('role', 'group');
    bar.setAttribute('aria-label', 'Camera framing');
    const label = document.createElement('span');
    label.className = 'oct-framing-label';
    label.textContent = 'View';
    bar.append(label);
    for (const [id, text] of [
      ['chamber', 'Chamber'],
      ['engine', 'Engine'],
    ] as const) {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = text;
      b.title = `Frame the ${id === 'chamber' ? 'combustion chamber' : 'whole engine'} (F toggles)`;
      b.addEventListener('click', () => {
        opts.onFraming(id);
        b.blur();
      });
      this.buttons.set(id, b);
      bar.append(b);
    }
    const kbd = document.createElement('kbd');
    kbd.textContent = 'F';
    bar.append(kbd);

    dock.append(this.legendBox, bar);
    this.el.append(this.status, dock);
    host.appendChild(this.el);
  }

  /** Show a status card (null hides it). */
  setStatus(title: string | null, detail = '', error = false): void {
    if (title === null) {
      this.status.hidden = true;
      return;
    }
    this.status.hidden = false;
    this.status.classList.toggle('is-error', error);
    this.status.setAttribute('role', error ? 'alert' : 'status');
    this.statusTitle.textContent = title;
    this.statusDetail.textContent = detail;
  }

  setFraming(f: Framing): void {
    for (const [id, b] of this.buttons) {
      const on = id === f;
      b.classList.toggle('is-on', on);
      b.setAttribute('aria-pressed', String(on));
    }
  }

  setLegendVisible(on: boolean): void {
    this.legendBox.hidden = !on;
  }

  dispose(): void {
    this.el.remove();
  }
}
