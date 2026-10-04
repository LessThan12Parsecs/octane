/**
 * DOM overlay inside the 3D viewport (owned by the app, not the UI layer, so it
 * survives engine switches): start-up / error status, the engine picker,
 * camera-framing buttons, the temperature legend and the simulation-limited hint.
 */
import type { Framing } from './sync';
import './app.css';

/** An engine the picker offers. */
export interface EngineChoice {
  id: string;
  label: string;
  description?: string;
}

export interface ViewportOverlayOptions {
  onFraming: (f: Framing) => void;
  /** Legend element for the 'temperature' render mode (from createTemperatureLegendElement). */
  legend?: HTMLElement;
  /** Engines for the picker (omitted or fewer than two: no picker). */
  engines?: readonly EngineChoice[];
  /** Initially selected engine id. */
  engineId?: string;
  /** The user picked an engine. */
  onEngine?: (id: string) => void;
}

/** A button on the status card (e.g. "Back to CFR F-1"). */
export interface StatusAction {
  label: string;
  onClick: () => void;
}

export class ViewportOverlay {
  readonly el: HTMLElement;
  private readonly status: HTMLElement;
  private readonly statusTitle: HTMLElement;
  private readonly statusDetail: HTMLElement;
  private readonly statusActions: HTMLElement;
  private readonly legendBox: HTMLElement;
  private readonly limited: HTMLElement;
  private readonly buttons = new Map<Framing, HTMLButtonElement>();
  private readonly engineSelect: HTMLSelectElement | null = null;

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
    this.statusActions = document.createElement('span');
    this.statusActions.className = 'oct-stage-status-actions';
    this.statusActions.hidden = true;
    this.status.append(this.statusTitle, this.statusDetail, this.statusActions);

    this.limited = document.createElement('div');
    this.limited.className = 'oct-stage-limited';
    this.limited.hidden = true;
    this.limited.setAttribute('role', 'status');
    this.limited.textContent = 'Simulation-limited: playback waits for the physics (slow down with [ )';

    const dock = document.createElement('div');
    dock.className = 'oct-stage-dock';

    this.legendBox = document.createElement('div');
    this.legendBox.className = 'oct-stage-legend';
    this.legendBox.hidden = true;
    if (opts.legend) this.legendBox.append(opts.legend);

    const rows: HTMLElement[] = [this.legendBox];
    const engines = opts.engines ?? [];
    if (engines.length > 1) {
      const pick = document.createElement('label');
      pick.className = 'oct-framing oct-engine-pick';
      const label = document.createElement('span');
      label.className = 'oct-framing-label';
      label.textContent = 'Engine';
      const sel = document.createElement('select');
      sel.setAttribute('aria-label', 'Engine');
      for (const e of engines) {
        const o = document.createElement('option');
        o.value = e.id;
        o.textContent = e.label;
        if (e.description) o.title = e.description;
        sel.appendChild(o);
      }
      if (opts.engineId) sel.value = opts.engineId;
      sel.addEventListener('change', () => {
        opts.onEngine?.(sel.value);
        sel.blur();
      });
      pick.append(label, sel);
      this.engineSelect = sel;
      rows.push(pick);
    }

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
    rows.push(bar);

    dock.append(...rows);
    this.el.append(this.status, this.limited, dock);
    host.appendChild(this.el);
  }

  /** Show a status card (null hides it), optionally with action buttons. */
  setStatus(title: string | null, detail = '', error = false, actions: readonly StatusAction[] = []): void {
    if (title === null) {
      this.status.hidden = true;
      this.statusActions.replaceChildren();
      this.statusActions.hidden = true;
      return;
    }
    this.status.hidden = false;
    this.status.classList.toggle('is-error', error);
    this.status.setAttribute('role', error ? 'alert' : 'status');
    this.statusTitle.textContent = title;
    this.statusDetail.textContent = detail;
    this.statusActions.replaceChildren(
      ...actions.map((a) => {
        const b = document.createElement('button');
        b.type = 'button';
        b.textContent = a.label;
        b.addEventListener('click', () => a.onClick());
        return b;
      }),
    );
    this.statusActions.hidden = actions.length === 0;
  }

  /** Whether a status card is showing. */
  get statusVisible(): boolean {
    return !this.status.hidden;
  }

  setFraming(f: Framing): void {
    for (const [id, b] of this.buttons) {
      const on = id === f;
      b.classList.toggle('is-on', on);
      b.setAttribute('aria-pressed', String(on));
    }
  }

  /** Reflect the current engine in the picker (no callback). */
  setEngine(id: string): void {
    if (this.engineSelect && this.engineSelect.value !== id) this.engineSelect.value = id;
  }

  /** Enable / disable the engine picker (disabled while a switch is in progress). */
  setEngineEnabled(on: boolean): void {
    if (this.engineSelect) this.engineSelect.disabled = !on;
  }

  setLegendVisible(on: boolean): void {
    this.legendBox.hidden = !on;
  }

  /** Show the "simulation-limited" hint (playback is waiting for the simulator). */
  setLimited(on: boolean): void {
    if (this.limited.hidden === !on) return;
    this.limited.hidden = !on;
  }

  dispose(): void {
    this.el.remove();
  }
}
