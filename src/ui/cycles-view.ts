/**
 * "Cycles" tab: cycle-to-cycle statistics, an IMEP-vs-cycle sparkline and the
 * per-cycle results table (metrics as rows, most recent cycles as columns).
 *
 * Multi-cylinder engines: a cylinder selector picks whose cycles the tiles,
 * sparkline and table show; engine-level results (EngineCycleSummary: speed,
 * brake torque and power, BMEP, FMEP, volumetric efficiency, BSFC) appear as
 * extra tiles and an 'Engine' group of rows as soon as the simulator reports them.
 */
import type { CycleSummary, EngineCycleSummary } from '../physics/core/snapshot';
import { SERIES, STATUS, INK, withAlpha } from './charts/theme';
import {
  computeCycleStats,
  CYCLE_METRICS,
  cyclesToCsv,
  ENGINE_METRICS,
  engineMetricMean,
  formatMetric,
  metricMean,
  type CycleStats,
} from './cycle-stats';
import { RingBuffer } from './ring-buffer';
import { fmt, nmToLbft, paToBar, wattsToHp, wattsToKW } from './units';

export interface CyclesViewOptions {
  /** Cylinders of the engine (> 1 adds the cylinder selector). */
  cylinders?: number;
  /** Initially selected cylinder, 0-based. */
  cylinder?: number;
  /** The user picked a cylinder. */
  onCylinder?: (index: number) => void;
}

const TABLE_COLS = 8;
const HISTORY = 300;
const SPARK_CYCLES = 150;

export class CyclesView {
  readonly el: HTMLElement;
  private readonly perCylinder: RingBuffer<CycleSummary>[] = [];
  private selected = 0;
  /** Engine summaries by engine cycle number (attached to cylinder 1's summaries). */
  private readonly engineByCycle = new Map<number, EngineCycleSummary>();
  private readonly engineHistory = new RingBuffer<EngineCycleSummary>(HISTORY);
  private readonly engineEls: HTMLElement[] = [];
  private readonly engineCells: HTMLTableCellElement[][] = [];
  private readonly engineMeanCells: HTMLTableCellElement[] = [];
  private readonly cylButtons: HTMLButtonElement[] = [];
  private engineShown = false;
  private readonly statEls = new Map<string, HTMLElement>();
  private readonly canvas: HTMLCanvasElement;
  private readonly headCells: HTMLTableCellElement[] = [];
  private readonly meanCells: HTMLTableCellElement[] = [];
  private readonly cells: HTMLTableCellElement[][] = []; // [metric][col]
  private readonly empty: HTMLElement;
  private dirty = true;
  private readonly ro: ResizeObserver | null = null;

  constructor(host: HTMLElement, opts: CyclesViewOptions = {}) {
    const nCyl = Math.max(1, Math.floor(opts.cylinders ?? 1));
    for (let i = 0; i < nCyl; i++) this.perCylinder.push(new RingBuffer<CycleSummary>(HISTORY));
    this.selected = Math.min(nCyl - 1, Math.max(0, opts.cylinder ?? 0));
    this.el = document.createElement('div');
    this.el.className = 'oct-cycles';
    host.appendChild(this.el);

    // --- cylinder selector ---
    if (nCyl > 1) {
      const bar = document.createElement('div');
      bar.className = 'oct-toolbar';
      const lbl = document.createElement('span');
      lbl.className = 'oct-toolbar-label';
      lbl.textContent = 'Cylinder';
      const seg = document.createElement('div');
      seg.className = 'oct-seg oct-seg-cyl';
      seg.setAttribute('role', 'group');
      seg.setAttribute('aria-label', 'Cylinder whose cycles are shown');
      for (let i = 0; i < nCyl; i++) {
        const b = document.createElement('button');
        b.type = 'button';
        b.textContent = String(i + 1);
        b.title = `Cycle results of cylinder ${i + 1}`;
        b.addEventListener('click', () => {
          this.setCylinder(i);
          opts.onCylinder?.(i);
        });
        seg.appendChild(b);
        this.cylButtons.push(b);
      }
      bar.append(lbl, seg);
      this.el.appendChild(bar);
    }

    // --- stat tiles ---
    const stats = document.createElement('div');
    stats.className = 'oct-stats';
    const tiles: [string, string, string][] = [
      ['imep', 'IMEP net [bar]', 'Mean net IMEP ± standard deviation'],
      ['cov', 'COV of IMEP', 'Coefficient of variation of net IMEP (cyclic variability)'],
      ['lnv', 'LNV of IMEP', 'Lowest normalised value: min IMEP / mean IMEP'],
      ['eta', 'η indicated', 'Mean net indicated efficiency'],
      ['ca50', 'CA50 mean', 'Mean 50 % burn angle'],
      ['knock', 'Knocking cycles', 'Share of cycles with end-gas autoignition'],
      ['mis', 'Misfires', 'Cycles whose spark kernel failed'],
      ['n', 'Cycles in window', 'Cycles in the statistics window'],
    ];
    for (const [k, label, hint] of tiles) {
      const t = document.createElement('div');
      t.className = 'oct-stat';
      t.title = hint;
      t.innerHTML = `<div class="oct-stat-label"></div><div class="oct-stat-value">—</div>`;
      (t.firstElementChild as HTMLElement).textContent = label;
      this.statEls.set(k, t.lastElementChild as HTMLElement);
      stats.appendChild(t);
    }
    // Engine tiles (hidden until the simulator reports engine-level results).
    const engineTiles: [string, string, string][] = [
      ['bT', 'Brake torque [N·m]', 'Mean brake torque over the statistics window: engine output at the crankshaft, excludes the car’s inertia (lb·ft in brackets)'],
      ['bP', 'Brake power [kW]', 'Mean brake power over the statistics window (hp in brackets)'],
    ];
    for (const [k, label, hint] of engineTiles) {
      const t = document.createElement('div');
      t.className = 'oct-stat oct-stat-engine';
      t.title = hint;
      t.hidden = true;
      t.innerHTML = `<div class="oct-stat-label"></div><div class="oct-stat-value">—</div>`;
      (t.firstElementChild as HTMLElement).textContent = label;
      this.statEls.set(k, t.lastElementChild as HTMLElement);
      this.engineEls.push(t);
      stats.appendChild(t);
    }
    this.el.appendChild(stats);

    // --- sparkline ---
    const sw = document.createElement('div');
    sw.className = 'oct-sparkline';
    const cap = document.createElement('div');
    cap.className = 'oct-sparkline-cap';
    cap.innerHTML =
      `<span>IMEP net vs cycle</span>` +
      `<span class="oct-swatch" style="--sw:${SERIES.blue}">IMEP</span>` +
      `<span class="oct-swatch is-dashed" style="--sw:${INK.secondary}">mean ± σ</span>` +
      `<span class="oct-swatch is-dot" style="--sw:${STATUS.serious}">knock</span>` +
      `<span class="oct-swatch is-dot" style="--sw:${STATUS.critical}">misfire</span>`;
    this.canvas = document.createElement('canvas');
    this.canvas.className = 'oct-sparkline-canvas';
    this.canvas.setAttribute('role', 'img');
    this.canvas.setAttribute('aria-label', 'Net IMEP of recent cycles');
    sw.append(cap, this.canvas);
    this.el.appendChild(sw);

    // --- actions ---
    const actions = document.createElement('div');
    actions.className = 'oct-table-actions';
    const csv = document.createElement('button');
    csv.type = 'button';
    csv.className = 'oct-mini-btn';
    csv.textContent = 'Export CSV';
    csv.title = 'Download all recorded cycles as CSV';
    csv.addEventListener('click', () => this.exportCsv());
    const clear = document.createElement('button');
    clear.type = 'button';
    clear.className = 'oct-mini-btn';
    clear.textContent = 'Clear';
    clear.title = 'Forget recorded cycles';
    clear.addEventListener('click', () => this.clear());
    const note = document.createElement('span');
    note.className = 'oct-toolbar-hint';
    note.textContent = `newest first · last ${TABLE_COLS} cycles`;
    actions.append(note, csv, clear);
    this.el.appendChild(actions);

    // --- table ---
    const wrap = document.createElement('div');
    wrap.className = 'oct-table-wrap';
    const table = document.createElement('table');
    table.className = 'oct-table';
    const thead = table.createTHead();
    const hr = thead.insertRow();
    const th0 = document.createElement('th');
    th0.textContent = 'Metric';
    th0.className = 'oct-col-metric';
    hr.appendChild(th0);
    const thMean = document.createElement('th');
    thMean.textContent = 'mean';
    thMean.className = 'oct-col-mean';
    thMean.title = 'Mean over the statistics window';
    hr.appendChild(thMean);
    for (let j = 0; j < TABLE_COLS; j++) {
      const th = document.createElement('th');
      th.textContent = '';
      hr.appendChild(th);
      this.headCells.push(th);
    }
    const tbody = table.createTBody();
    for (const def of CYCLE_METRICS) {
      if (def.group) {
        const gr = tbody.insertRow();
        gr.className = 'oct-group-row';
        const gc = gr.insertCell();
        gc.colSpan = TABLE_COLS + 2;
        gc.textContent = def.group;
      }
      const row = tbody.insertRow();
      const name = row.insertCell();
      name.className = 'oct-col-metric';
      name.title = def.hint;
      name.innerHTML = `<span class="oct-metric-name"></span><span class="oct-metric-unit"></span>`;
      (name.firstElementChild as HTMLElement).textContent = def.label;
      (name.lastElementChild as HTMLElement).textContent = def.unit;
      const mc = row.insertCell();
      mc.className = 'oct-col-mean';
      this.meanCells.push(mc);
      const rowCells: HTMLTableCellElement[] = [];
      for (let j = 0; j < TABLE_COLS; j++) rowCells.push(row.insertCell());
      this.cells.push(rowCells);
    }
    // Engine rows: the engine summary of the same (engine) cycle as each column.
    const egr = tbody.insertRow();
    egr.className = 'oct-group-row';
    egr.hidden = true;
    const egc = egr.insertCell();
    egc.colSpan = TABLE_COLS + 2;
    egc.textContent = 'Engine';
    this.engineEls.push(egr);
    for (const def of ENGINE_METRICS) {
      const row = tbody.insertRow();
      row.hidden = true;
      this.engineEls.push(row);
      const name = row.insertCell();
      name.className = 'oct-col-metric';
      name.title = def.hint;
      name.innerHTML = `<span class="oct-metric-name"></span><span class="oct-metric-unit"></span>`;
      (name.firstElementChild as HTMLElement).textContent = def.label;
      (name.lastElementChild as HTMLElement).textContent = def.unit;
      const mc = row.insertCell();
      mc.className = 'oct-col-mean';
      this.engineMeanCells.push(mc);
      const rowCells: HTMLTableCellElement[] = [];
      for (let j = 0; j < TABLE_COLS; j++) rowCells.push(row.insertCell());
      this.engineCells.push(rowCells);
    }
    wrap.appendChild(table);
    this.el.appendChild(wrap);

    this.empty = document.createElement('div');
    this.empty.className = 'oct-empty';
    this.empty.textContent = 'No complete cycles yet — results appear after each 720° cycle.';
    this.el.appendChild(this.empty);

    if (typeof ResizeObserver !== 'undefined') {
      this.ro = new ResizeObserver(() => {
        this.dirty = true;
      });
      this.ro.observe(this.canvas);
    }
    if (nCyl > 1) this.setCylinder(this.selected);
  }

  /** Ring of the selected cylinder. */
  private get cycles(): RingBuffer<CycleSummary> {
    return this.perCylinder[this.selected];
  }

  /** Latest summary of the selected cylinder. */
  get latest(): CycleSummary | null {
    return this.cycles.latest() ?? null;
  }

  /** Latest summary of cylinder `index` (0-based). */
  latestFor(index: number): CycleSummary | null {
    return this.perCylinder[index]?.latest() ?? null;
  }

  /** Latest engine-level summary, or null. */
  get latestEngine(): EngineCycleSummary | null {
    return this.engineHistory.latest() ?? null;
  }

  /** Recorded summaries of the selected cylinder, oldest first. */
  get all(): CycleSummary[] {
    return this.cycles.toArray();
  }

  get cylinder(): number {
    return this.selected;
  }

  /** Show cylinder `index` (0-based). */
  setCylinder(index: number): void {
    const i = Math.min(this.perCylinder.length - 1, Math.max(0, Math.floor(index)));
    this.selected = i;
    this.cylButtons.forEach((b, k) => {
      b.classList.toggle('is-on', k === i);
      b.setAttribute('aria-pressed', String(k === i));
    });
    this.dirty = true;
  }

  push(c: CycleSummary): void {
    const k = Math.min(this.perCylinder.length - 1, Math.max(0, c.cylinder ?? 0));
    this.perCylinder[k].push(c);
    if (c.engine) {
      this.engineByCycle.set(c.cycle, c.engine);
      this.engineHistory.push(c.engine);
      if (this.engineByCycle.size > HISTORY) {
        const oldest = this.engineByCycle.keys().next().value;
        if (oldest !== undefined) this.engineByCycle.delete(oldest);
      }
    }
    this.dirty = true;
  }

  clear(): void {
    for (const r of this.perCylinder) r.clear();
    this.engineByCycle.clear();
    this.engineHistory.clear();
    this.dirty = true;
  }

  stats(): CycleStats {
    return computeCycleStats(this.window());
  }

  private window(): CycleSummary[] {
    const all = this.cycles.toArray();
    return all.slice(Math.max(0, all.length - SPARK_CYCLES));
  }

  /** Re-render if something changed and the tab is visible. */
  render(visible: boolean): void {
    if (!visible || !this.dirty) return;
    this.dirty = false;
    const win = this.window();
    const n = win.length;
    this.empty.style.display = n ? 'none' : '';

    const st = computeCycleStats(win);
    const set = (k: string, v: string): void => {
      const el = this.statEls.get(k);
      if (el && el.textContent !== v) el.textContent = v;
    };
    set('imep', n ? `${fmt(st.imepMean, 2)} ± ${fmt(st.imepStd, 2)}` : '—');
    set('cov', `${fmt(st.imepCov, 2)} %`);
    set('lnv', `${fmt(st.imepLnv, 1)} %`);
    set('eta', `${fmt(st.etaMean, 1)} %`);
    set('ca50', `${fmt(st.ca50Mean, 1)}°`);
    set('knock', n ? `${fmt(st.knockFraction * 100, 0)} %` : '—');
    set('mis', String(st.misfires));
    set('n', String(n));
    this.statEls.get('knock')?.classList.toggle('is-alert', st.knockFraction > 0);
    this.statEls.get('mis')?.classList.toggle('is-alert', st.misfires > 0);

    const eng = this.engineHistory.toArray().slice(-SPARK_CYCLES);
    const hasEngine = eng.length > 0;
    if (hasEngine !== this.engineShown) {
      this.engineShown = hasEngine;
      for (const e of this.engineEls) e.hidden = !hasEngine;
    }
    if (hasEngine) {
      const def = (key: string) => ENGINE_METRICS.find((d) => d.key === key)!;
      const tq = engineMetricMean(def('brakeTorque'), eng);
      const pw = engineMetricMean(def('brakePower'), eng);
      set('bT', `${fmt(tq, 1)} (${fmt(nmToLbft(tq), 1)})`);
      set('bP', `${fmt(pw, 2)} (${fmt(wattsToHp(pw * 1e3), 1)} hp)`);
    }

    // Table: newest first.
    for (let j = 0; j < TABLE_COLS; j++) {
      const c = this.cycles.latest(j);
      this.headCells[j].textContent = c ? `#${c.cycle}` : '';
      this.headCells[j].classList.toggle('is-newest', j === 0 && !!c);
      for (let r = 0; r < CYCLE_METRICS.length; r++) {
        const def = CYCLE_METRICS[r];
        const cell = this.cells[r][j];
        const txt = c ? formatMetric(def, c) : '';
        if (cell.textContent !== txt) cell.textContent = txt;
        const flag = c && def.flag ? def.flag(c) : null;
        cell.className = flag ? `is-${flag}` : '';
      }
    }
    for (let r = 0; r < CYCLE_METRICS.length; r++) {
      const def = CYCLE_METRICS[r];
      const txt = def.noMean || !n ? '' : fmt(metricMean(def, win), def.decimals);
      if (this.meanCells[r].textContent !== txt) this.meanCells[r].textContent = txt;
    }
    if (hasEngine) {
      for (let j = 0; j < TABLE_COLS; j++) {
        const c = this.cycles.latest(j);
        const e = c ? (c.engine ?? this.engineByCycle.get(c.cycle)) : undefined;
        for (let r = 0; r < ENGINE_METRICS.length; r++) {
          const txt = e ? fmt(ENGINE_METRICS[r].value(e), ENGINE_METRICS[r].decimals) : '';
          const cell = this.engineCells[r][j];
          if (cell.textContent !== txt) cell.textContent = txt;
        }
      }
      for (let r = 0; r < ENGINE_METRICS.length; r++) {
        const def = ENGINE_METRICS[r];
        const txt = fmt(engineMetricMean(def, eng), def.decimals);
        if (this.engineMeanCells[r].textContent !== txt) this.engineMeanCells[r].textContent = txt;
      }
    }
    this.drawSparkline(win, st);
  }

  private drawSparkline(win: CycleSummary[], st: CycleStats): void {
    const cv = this.canvas;
    const dpr = globalThis.devicePixelRatio || 1;
    const w = Math.max(100, cv.clientWidth);
    const h = Math.max(40, cv.clientHeight);
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr);
      cv.height = Math.round(h * dpr);
    }
    const ctx = cv.getContext('2d');
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    const n = win.length;
    const padL = 40;
    const padR = 6;
    const padT = 6;
    const padB = 14;
    const pw = w - padL - padR;
    const ph = h - padT - padB;
    ctx.font = '10px ui-monospace, Menlo, monospace';
    ctx.fillStyle = INK.muted;
    if (n === 0) {
      ctx.fillText('waiting for cycles…', padL, padT + ph / 2);
      return;
    }
    let lo = Infinity;
    let hi = -Infinity;
    const ys = win.map((c) => paToBar(c.imepNet));
    for (const y of ys) {
      if (!Number.isFinite(y)) continue;
      lo = Math.min(lo, y);
      hi = Math.max(hi, y);
    }
    if (!(hi >= lo)) {
      lo = 0;
      hi = 1;
    }
    if (hi - lo < 0.2) {
      const m = (hi + lo) / 2;
      lo = m - 0.1;
      hi = m + 0.1;
    }
    const span = hi - lo;
    lo -= span * 0.08;
    hi += span * 0.08;
    // Newest cycle at the right edge; at least 30 slots so early points don't stretch.
    const slots = Math.max(n, 30);
    const X = (i: number): number => padL + ((slots - n + i) / (slots - 1)) * pw;
    const Y = (v: number): number => padT + (1 - (v - lo) / (hi - lo)) * ph;

    // axes labels
    ctx.textAlign = 'right';
    ctx.textBaseline = 'middle';
    ctx.fillText(fmt(hi, 2), padL - 4, padT + 4);
    ctx.fillText(fmt(lo, 2), padL - 4, padT + ph - 4);
    ctx.textAlign = 'left';
    ctx.textBaseline = 'alphabetic';
    ctx.fillText(`#${win[0].cycle}`, padL, h - 2);
    ctx.textAlign = 'right';
    ctx.fillText(`#${win[n - 1].cycle}  (bar)`, w - padR, h - 2);

    // mean ± σ band
    if (Number.isFinite(st.imepMean)) {
      if (Number.isFinite(st.imepStd)) {
        ctx.fillStyle = withAlpha('#a3acb9', 0.08);
        const y0 = Y(st.imepMean + st.imepStd);
        const y1 = Y(st.imepMean - st.imepStd);
        ctx.fillRect(padL, y0, pw, y1 - y0);
      }
      ctx.strokeStyle = INK.secondary;
      ctx.setLineDash([4, 3]);
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(padL, Y(st.imepMean));
      ctx.lineTo(padL + pw, Y(st.imepMean));
      ctx.stroke();
      ctx.setLineDash([]);
    }

    // IMEP line
    ctx.strokeStyle = SERIES.blue;
    ctx.lineWidth = 1.5;
    ctx.lineJoin = 'round';
    ctx.beginPath();
    let pen = false;
    for (let i = 0; i < n; i++) {
      const y = ys[i];
      if (!Number.isFinite(y)) {
        pen = false;
        continue;
      }
      if (pen) ctx.lineTo(X(i), Y(y));
      else {
        ctx.moveTo(X(i), Y(y));
        pen = true;
      }
    }
    ctx.stroke();

    // knock / misfire markers
    for (let i = 0; i < n; i++) {
      const c = win[i];
      const y = Number.isFinite(ys[i]) ? ys[i] : lo;
      if (c.misfire) {
        ctx.strokeStyle = STATUS.critical;
        ctx.lineWidth = 1.5;
        const x = X(i);
        const yy = Y(y);
        ctx.beginPath();
        ctx.moveTo(x - 3, yy - 3);
        ctx.lineTo(x + 3, yy + 3);
        ctx.moveTo(x + 3, yy - 3);
        ctx.lineTo(x - 3, yy + 3);
        ctx.stroke();
      } else if (Number.isFinite(c.knockOnsetDeg)) {
        ctx.fillStyle = STATUS.serious;
        ctx.beginPath();
        ctx.arc(X(i), Y(y), 2.5, 0, 2 * Math.PI);
        ctx.fill();
      }
    }
  }

  private exportCsv(): void {
    // Every cylinder, in completion order (engine columns on the summaries that carry them).
    const all = this.perCylinder.length > 1 ? this.perCylinder.flatMap((r) => r.toArray()).sort((a, b) => a.cycle - b.cycle || (a.cylinder ?? 0) - (b.cylinder ?? 0)) : this.cycles.toArray();
    if (!all.length) return;
    const blob = new Blob([cyclesToCsv(all)], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `octane-cycles-${all[0].cycle}-${all[all.length - 1].cycle}.csv`;
    document.body.appendChild(a);
    a.click();
    a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  }

  invalidate(): void {
    this.dirty = true;
  }

  dispose(): void {
    this.ro?.disconnect();
    this.el.remove();
  }
}
