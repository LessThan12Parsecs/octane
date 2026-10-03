/**
 * lil-gui control panel: operating point (UI units), fuel, intake/ambient, CFR
 * rating presets, view options and playback.
 */
import GUI, { type Controller } from 'lil-gui';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import { FUEL_OPTIONS, mergeOp, paramsFromOp, patchForParam, type OpParamKey, type OpParams } from './op-binding';
import { formatTimeScale, sliderToTimeScale, STEP_LARGE_DEG, STEP_SMALL_DEG, timeScaleToSlider } from './playback';
import { getPreset, monSparkAdvance, type PresetId } from './presets';

export type RenderMode = 'physical' | 'temperature';

export interface ControlsCallbacks {
  onOperatingPointChange: (patch: Partial<OperatingPoint>) => void;
  onViewChange: (v: { cutaway: boolean; mode: RenderMode }) => void;
  onTogglePause: () => void;
  onTimeScale: (ts: number) => void;
  onStep: (deg: number) => void;
  onReset: () => void;
}

interface ViewParams {
  cutaway: boolean;
  mode: RenderMode;
}

interface PresetParams {
  monSchedule: boolean;
  source: string;
}

export class Controls {
  readonly gui: GUI;
  readonly params: OpParams;
  private op: OperatingPoint;
  readonly view: ViewParams = { cutaway: true, mode: 'physical' };
  private readonly presetParams: PresetParams = { monSchedule: false, source: '' };
  private readonly playback = { slider: 0, paused: false };
  private readonly ctrls = new Map<OpParamKey, Controller>();
  private readonly tsCtrl: Controller;
  private readonly pauseBtn: HTMLButtonElement;
  private readonly monCtrl: Controller;
  private readonly sourceCtrl: Controller;

  constructor(
    host: HTMLElement,
    spec: EngineSpec,
    initial: OperatingPoint,
    private readonly cb: ControlsCallbacks,
  ) {
    this.op = mergeOp(initial, {});
    this.params = paramsFromOp(this.op);
    const narrow = typeof window !== 'undefined' && window.innerWidth < 760;
    this.gui = new GUI({ container: host, title: 'Octane · CFR F-1', width: 284, closeFolders: false });
    this.gui.domElement.classList.add('oct-gui');
    if (narrow) this.gui.close();

    const p = this.params;
    const bind = (c: Controller, key: OpParamKey): Controller => {
      this.ctrls.set(key, c);
      c.onChange(() => this.changed(key));
      return c;
    };

    // --- engine & combustion ---
    const fEng = this.gui.addFolder('Engine');
    const [crMin, crMax] = spec.geometry.compressionRatioRange;
    bind(fEng.add(p, 'compressionRatio', crMin, crMax, 0.05).decimals(2).name('Compression ratio'), 'compressionRatio');
    bind(fEng.add(p, 'speedMode', { 'Fixed (synchronous motor)': 'fixed', 'Free (crank dynamics)': 'free' }).name('Speed control'), 'speedMode');
    bind(fEng.add(p, 'rpm', 200, 3000, 10).name('Speed [rpm]'), 'rpm');
    bind(fEng.add(p, 'loadTorque', -20, 100, 0.5).name('Load torque [N·m]'), 'loadTorque');
    const throttle = bind(fEng.add(p, 'throttlePct', 0, 100, 1).name('Throttle [%]'), 'throttlePct');
    throttle.domElement.title =
      'The real CFR runs wide open; here the throttle scales the carburettor venturi. ' +
      'With fixed speed (dynamometer) the rpm cannot change: watch the Intake (MAP) and Last-cycle IMEP readouts. ' +
      'Switch Speed control to Free to let the load change the speed.';
    bind(fEng.add(p, 'coolantC', 20, 130, 1).name('Coolant [°C]'), 'coolantC');

    const fIgn = this.gui.addFolder('Ignition');
    bind(fIgn.add(p, 'sparkAdvanceDeg', -10, 60, 0.5).name('Spark advance [° BTDC]'), 'sparkAdvanceDeg');
    bind(fIgn.add(p, 'dwellMs', 0.5, 10, 0.1).name('Dwell [ms]'), 'dwellMs');

    const fFuel = this.gui.addFolder('Fuel & mixture');
    bind(fFuel.add(p, 'fuel', FUEL_OPTIONS).name('Fuel'), 'fuel');
    bind(fFuel.add(p, 'octaneNumber', 0, 100, 0.5).name('PRF octane no.'), 'octaneNumber');
    bind(fFuel.add(p, 'equivalenceRatio', 0.4, 2.0, 0.01).name('Equivalence ratio φ'), 'equivalenceRatio');
    bind(fFuel.add(p, 'egrPct', 0, 40, 0.5).name('EGR [%]'), 'egrPct');

    const fAir = this.gui.addFolder('Intake & ambient');
    bind(fAir.add(p, 'intakeMixtureC', 0, 200, 1).name('Mixture temp. [°C]'), 'intakeMixtureC');
    bind(fAir.add(p, 'ambientKPa', 50, 110, 0.1).name('Ambient p [kPa]'), 'ambientKPa');
    bind(fAir.add(p, 'ambientC', -30, 50, 0.5).name('Ambient T [°C]'), 'ambientC');
    bind(fAir.add(p, 'humidityPct', 0, 100, 1).name('Rel. humidity [%]'), 'humidityPct');
    fAir.close();

    // --- presets ---
    const fPre = this.gui.addFolder('CFR rating presets');
    fPre.add({ ron: () => this.applyPreset('RON') }, 'ron').name('CFR Research method (RON)');
    fPre.add({ mon: () => this.applyPreset('MON') }, 'mon').name('CFR Motor method (MON)');
    this.monCtrl = fPre
      .add(this.presetParams, 'monSchedule')
      .name('MON spark tracks CR')
      .onChange(() => {
        if (this.presetParams.monSchedule) this.applyMonSpark();
      });
    this.sourceCtrl = fPre.add(this.presetParams, 'source').name('Conditions from').disable();
    this.sourceCtrl.hide();

    // --- view ---
    const fView = this.gui.addFolder('View');
    fView.add(this.view, 'cutaway').name('Cutaway').onChange(() => this.cb.onViewChange({ ...this.view }));
    fView
      .add(this.view, 'mode', { 'Physical (flame, spark)': 'physical', 'Gas temperature': 'temperature' })
      .name('Render mode')
      .onChange(() => this.cb.onViewChange({ ...this.view }));

    // --- playback ---
    const fPlay = this.gui.addFolder('Playback');
    this.tsCtrl = fPlay
      .add(this.playback, 'slider', 0, 1, 0.001)
      .name('Time scale')
      .onChange(() => {
        const ts = sliderToTimeScale(this.playback.slider);
        this.tsCtrl.name(`Time scale ${formatTimeScale(ts)}`);
        this.cb.onTimeScale(ts);
      });
    this.tsCtrl.domElement.classList.add('oct-log-slider');
    const row = document.createElement('div');
    row.className = 'oct-gui-row';
    const mk = (label: string, title: string, fn: () => void): HTMLButtonElement => {
      const b = document.createElement('button');
      b.type = 'button';
      b.textContent = label;
      b.title = title;
      b.addEventListener('click', fn);
      row.appendChild(b);
      return b;
    };
    mk(`−${STEP_LARGE_DEG}°`, `Step back ${STEP_LARGE_DEG}° (Shift+←)`, () => this.cb.onStep(-STEP_LARGE_DEG));
    mk(`−${STEP_SMALL_DEG}°`, `Step back ${STEP_SMALL_DEG}° (←)`, () => this.cb.onStep(-STEP_SMALL_DEG));
    this.pauseBtn = mk('Pause', 'Play / pause (Space)', () => this.cb.onTogglePause());
    this.pauseBtn.classList.add('is-primary');
    mk(`+${STEP_SMALL_DEG}°`, `Step forward ${STEP_SMALL_DEG}° (→)`, () => this.cb.onStep(STEP_SMALL_DEG));
    mk(`+${STEP_LARGE_DEG}°`, `Step forward ${STEP_LARGE_DEG}° (Shift+→)`, () => this.cb.onStep(STEP_LARGE_DEG));
    fPlay.$children.appendChild(row);
    fPlay.add({ reset: () => this.cb.onReset() }, 'reset').name('Reset simulation');

    this.syncVisibility();
  }

  get operatingPoint(): OperatingPoint {
    return this.op;
  }

  private changed(key: OpParamKey): void {
    const patch = patchForParam(this.params, key);
    if (key === 'compressionRatio' && this.presetParams.monSchedule) {
      const adv = monSparkAdvance(this.params.compressionRatio);
      this.params.sparkAdvanceDeg = Math.round(adv * 10) / 10;
      patch.sparkAdvanceDeg = this.params.sparkAdvanceDeg;
      this.ctrls.get('sparkAdvanceDeg')?.updateDisplay();
    }
    if (key === 'sparkAdvanceDeg' && this.presetParams.monSchedule) {
      // Manual spark override leaves the MON schedule.
      this.presetParams.monSchedule = false;
      this.monCtrl.updateDisplay();
    }
    this.op = mergeOp(this.op, patch);
    this.syncVisibility();
    this.cb.onOperatingPointChange(patch);
  }

  private applyMonSpark(): void {
    const adv = Math.round(monSparkAdvance(this.params.compressionRatio) * 10) / 10;
    this.params.sparkAdvanceDeg = adv;
    this.ctrls.get('sparkAdvanceDeg')?.updateDisplay();
    this.op = mergeOp(this.op, { sparkAdvanceDeg: adv });
    this.cb.onOperatingPointChange({ sparkAdvanceDeg: adv });
  }

  applyPreset(id: PresetId): void {
    const info = getPreset(id, this.params.compressionRatio);
    this.presetParams.monSchedule = info.sparkFollowsCR;
    this.presetParams.source = info.source === 'cfr.ts' ? 'engines/cfr.ts' : 'built-in (ASTM, approx.)';
    this.sourceCtrl.show();
    this.monCtrl.updateDisplay();
    this.sourceCtrl.updateDisplay();
    this.setOperatingPoint(mergeOp(this.op, info.patch));
    this.cb.onOperatingPointChange(info.patch);
  }

  /** Replace the displayed operating point (no callback). */
  setOperatingPoint(op: OperatingPoint): void {
    this.op = mergeOp(op, {});
    Object.assign(this.params, paramsFromOp(this.op, this.params));
    for (const c of this.ctrls.values()) c.updateDisplay();
    this.syncVisibility();
  }

  setPlayback(timeScale: number, paused: boolean): void {
    this.playback.slider = timeScaleToSlider(timeScale);
    this.playback.paused = paused;
    this.tsCtrl.name(`Time scale ${formatTimeScale(timeScale)}`);
    this.tsCtrl.updateDisplay();
    this.setPlaybackButtons(paused);
  }

  /** Update only the play/pause label (leaves the time-scale slider alone while dragging). */
  setPlaybackButtons(paused: boolean): void {
    this.playback.paused = paused;
    this.pauseBtn.textContent = paused ? 'Play' : 'Pause';
  }

  setView(v: ViewParams): void {
    this.view.cutaway = v.cutaway;
    this.view.mode = v.mode;
    this.gui.controllersRecursive().forEach((c) => {
      if (c.object === this.view) c.updateDisplay();
    });
  }

  private syncVisibility(): void {
    const p = this.params;
    this.ctrls.get('octaneNumber')?.show(p.fuel === 'PRF');
    this.ctrls.get('loadTorque')?.show(p.speedMode === 'free');
    this.ctrls.get('rpm')?.name(p.speedMode === 'free' ? 'Initial speed [rpm]' : 'Speed [rpm]');
  }

  dispose(): void {
    this.gui.destroy();
  }
}
