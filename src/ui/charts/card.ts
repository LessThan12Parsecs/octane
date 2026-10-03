/** Collapsible chart card (title, live readout, actions, body). */
import { loadPref, savePref } from '../prefs';

export class Card {
  readonly el: HTMLElement;
  readonly body: HTMLElement;
  readonly readout: HTMLElement;
  readonly actions: HTMLElement;
  readonly legend: HTMLElement;
  private readonly toggleBtn: HTMLButtonElement;
  collapsed: boolean;
  /** In the scroll viewport (IntersectionObserver). */
  inView = true;
  /** Set when the card needs a redraw regardless of data changes. */
  stale = true;
  /** Scratch flag: redraw this frame. */
  redraw = false;
  private lastReadout = '';

  constructor(
    parent: HTMLElement,
    readonly id: string,
    title: string,
    hint: string,
  ) {
    this.el = document.createElement('section');
    this.el.className = 'oct-card';
    this.el.dataset.card = id;

    const head = document.createElement('header');
    head.className = 'oct-card-head';
    this.toggleBtn = document.createElement('button');
    this.toggleBtn.type = 'button';
    this.toggleBtn.className = 'oct-card-toggle';
    this.toggleBtn.title = hint;
    this.toggleBtn.innerHTML = `<span class="oct-chev" aria-hidden="true"></span><span class="oct-card-title"></span>`;
    (this.toggleBtn.querySelector('.oct-card-title') as HTMLElement).textContent = title;
    this.readout = document.createElement('span');
    this.readout.className = 'oct-card-readout';
    this.actions = document.createElement('span');
    this.actions.className = 'oct-card-actions';
    head.append(this.toggleBtn, this.readout, this.actions);

    this.legend = document.createElement('div');
    this.legend.className = 'oct-card-legend';

    this.body = document.createElement('div');
    this.body.className = 'oct-card-body';
    this.el.append(head, this.legend, this.body);
    parent.appendChild(this.el);

    const collapsedIds = loadPref<string[]>('collapsedCards', []);
    this.collapsed = collapsedIds.includes(id);
    this.applyCollapsed();
    this.toggleBtn.addEventListener('click', () => this.setCollapsed(!this.collapsed));
  }

  get active(): boolean {
    return !this.collapsed && this.inView;
  }

  setCollapsed(c: boolean): void {
    this.collapsed = c;
    this.applyCollapsed();
    const ids = new Set(loadPref<string[]>('collapsedCards', []));
    if (c) ids.add(this.id);
    else ids.delete(this.id);
    savePref('collapsedCards', [...ids]);
    this.stale = true;
  }

  private applyCollapsed(): void {
    this.el.classList.toggle('is-collapsed', this.collapsed);
    this.toggleBtn.setAttribute('aria-expanded', String(!this.collapsed));
  }

  setReadout(text: string): void {
    if (text === this.lastReadout) return;
    this.lastReadout = text;
    this.readout.textContent = text;
  }

  addAction(label: string, title: string, onClick: () => void): HTMLButtonElement {
    const b = document.createElement('button');
    b.type = 'button';
    b.className = 'oct-mini-btn';
    b.textContent = label;
    b.title = title;
    b.addEventListener('click', onClick);
    this.actions.appendChild(b);
    return b;
  }

  /** A legend of coloured swatches (for hook-drawn traces that uPlot's legend can't describe). */
  setSwatches(items: { label: string; color: string; dashed?: boolean }[]): void {
    this.legend.replaceChildren(
      ...items.map((it) => {
        const s = document.createElement('span');
        s.className = 'oct-swatch' + (it.dashed ? ' is-dashed' : '');
        s.style.setProperty('--sw', it.color);
        s.textContent = it.label;
        return s;
      }),
    );
  }

  plotHost(): HTMLElement {
    const el = document.createElement('div');
    el.className = 'oct-plot';
    this.body.appendChild(el);
    return el;
  }
}
