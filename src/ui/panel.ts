/**
 * Right-hand collapsible panel with "Traces" and "Cycles" tabs. Publishes its
 * footprint as CSS custom properties on <html> (--oct-panel-inset, --oct-hud-h)
 * so the 3D viewport (.oct-viewport) can avoid the panel and HUD.
 */
import { loadPref, savePref } from './prefs';

export type PanelTab = 'traces' | 'cycles';

export class SidePanel {
  readonly el: HTMLElement;
  readonly scroll: HTMLElement;
  readonly tracesHost: HTMLElement;
  readonly cyclesHost: HTMLElement;
  private readonly toggle: HTMLButtonElement;
  private readonly tabButtons = new Map<PanelTab, HTMLButtonElement>();
  open: boolean;
  tab: PanelTab;
  private readonly listeners = new Set<() => void>();
  private resizeTimer: ReturnType<typeof setTimeout> | null = null;

  constructor(host: HTMLElement) {
    const narrow = typeof window !== 'undefined' && window.innerWidth < 760;
    this.open = loadPref<boolean>('panelOpen', !narrow);
    if (narrow) this.open = false;
    this.tab = loadPref<PanelTab>('panelTab', 'traces') === 'cycles' ? 'cycles' : 'traces';

    this.el = document.createElement('aside');
    this.el.className = 'oct-panel';
    this.el.setAttribute('aria-label', 'Charts and cycle results');

    this.toggle = document.createElement('button');
    this.toggle.type = 'button';
    this.toggle.className = 'oct-panel-toggle';
    this.toggle.title = 'Show / hide the chart panel (C)';
    this.toggle.innerHTML = '<span class="oct-panel-toggle-icon" aria-hidden="true"></span><span class="oct-panel-toggle-label">Charts</span>';
    this.toggle.addEventListener('click', () => this.setOpen(!this.open));

    const head = document.createElement('header');
    head.className = 'oct-panel-head';
    const tabs = document.createElement('div');
    tabs.className = 'oct-tabs';
    tabs.setAttribute('role', 'tablist');
    for (const [id, label] of [
      ['traces', 'Traces'],
      ['cycles', 'Cycle results'],
    ] as const) {
      const b = document.createElement('button');
      b.type = 'button';
      b.setAttribute('role', 'tab');
      b.textContent = label;
      b.addEventListener('click', () => this.setTab(id));
      tabs.appendChild(b);
      this.tabButtons.set(id, b);
    }
    head.appendChild(tabs);

    this.scroll = document.createElement('div');
    this.scroll.className = 'oct-panel-scroll';
    this.tracesHost = document.createElement('div');
    this.tracesHost.className = 'oct-tab-page';
    this.tracesHost.setAttribute('role', 'tabpanel');
    this.cyclesHost = document.createElement('div');
    this.cyclesHost.className = 'oct-tab-page';
    this.cyclesHost.setAttribute('role', 'tabpanel');
    this.scroll.append(this.tracesHost, this.cyclesHost);

    this.el.append(this.toggle, head, this.scroll);
    host.appendChild(this.el);
    this.apply();
  }

  onChange(cb: () => void): () => void {
    this.listeners.add(cb);
    return () => this.listeners.delete(cb);
  }

  setOpen(open: boolean): void {
    if (open === this.open) return;
    this.open = open;
    savePref('panelOpen', open);
    this.apply();
    this.emit();
    // Let the 3D view re-fit once the slide transition has finished.
    if (this.resizeTimer) clearTimeout(this.resizeTimer);
    this.resizeTimer = setTimeout(() => {
      this.resizeTimer = null;
      window.dispatchEvent(new Event('resize'));
    }, 240);
  }

  setTab(tab: PanelTab): void {
    if (tab === this.tab) return;
    this.tab = tab;
    savePref('panelTab', tab);
    this.apply();
    this.emit();
  }

  private emit(): void {
    for (const cb of this.listeners) cb();
  }

  private apply(): void {
    this.el.classList.toggle('is-open', this.open);
    this.toggle.setAttribute('aria-expanded', String(this.open));
    document.documentElement.classList.toggle('oct-panel-open', this.open);
    for (const [id, b] of this.tabButtons) {
      const on = id === this.tab;
      b.classList.toggle('is-on', on);
      b.setAttribute('aria-selected', String(on));
    }
    this.tracesHost.hidden = this.tab !== 'traces';
    this.cyclesHost.hidden = this.tab !== 'cycles';
  }

  dispose(): void {
    if (this.resizeTimer) clearTimeout(this.resizeTimer);
    document.documentElement.classList.remove('oct-panel-open');
    this.el.remove();
  }
}
