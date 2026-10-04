/**
 * lil-gui control panel: operating point (UI units), fuel, intake/ambient, the
 * engine's presets (CFR: the octane-rating presets), view options and playback.
 * Which controls appear, their ranges and labels come from the engine's
 * EngineUiProfile (engine registry; engine-ui.ts).
 */
import GUI, { type Controller } from 'lil-gui';
import type { EngineSpec } from '../physics/core/engine-spec';
import type { OperatingPoint } from '../physics/core/operating-point';
import type { EngineDefinition, EnginePreset } from '../physics/engines/index';
import { controlVisibility, definitionForSpec, engineTitle, gearOptions, loadModelOptions, usesCfrRatingPresets } from './engine-ui';
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

const CFR_THROTTLE_HINT =
  'The real CFR runs wide open; here the throttle scales the carburettor venturi. ' +
  'With fixed speed (dynamometer) the rpm cannot change: watch the Intake (MAP) and Last-cycle IMEP readouts. ' +
  'Switch Speed control to Free to let the load change the speed.';

const THROTTLE_HINT =
  'Butterfly throttle downstream of the carburettor venturi. With fixed speed the rpm is held: watch the Intake (MAP) and ' +
  'torque readouts. In Free mode the load decides the speed.';

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
  private readonly monCtrl: Controller | null = null;
  private readonly sourceCtrl: Controller | null = null;
  private readonly def: EngineDefinition;
  private readonly cfrPresets: boolean;

  constructor(
    host: HTMLElement,
    spec: EngineSpec,
    initial: OperatingPoint,
    private readonly cb: ControlsCallbacks,
    engine?: EngineDefinition,
  ) {
    const def = (this.def = engine ?? definitionForSpec(spec));
    const ui = def.ui;
    this.cfrPresets = usesCfrRatingPresets(def);
    this.op = mergeOp(initial, {});
    this.params = paramsFromOp(this.op);
    const narrow = typeof window !== 'undefined' && window.innerWidth < 760;
    this.gui = new GUI({ container: host, title: engineTitle(def), width: 284, closeFolders: false });
    this.gui.domElement.classList.add('oct-gui');
    if (narrow) this.gui.close();

    const p = this.params;
    const bind = (c: Controller, key: OpParamKey): Controller => {
      this.ctrls.set(key, c);
      c.onChange(() => this.changed(key));
      return c;
    };
    const vis = controlVisibility(def, p);

    // --- engine & combustion ---
    const fEng = this.gui.addFolder('Engine');
    if (vis.compressionRatio) {
      const [crMin, crMax] = spec.geometry.compressionRatioRange;
      bind(fEng.add(p, 'compressionRatio', crMin, crMax, 0.05).decimals(2).name('Compression ratio'), 'compressionRatio');
    }
    bind(fEng.add(p, 'speedMode', { [ui.speedModeLabels.fixed]: 'fixed', [ui.speedModeLabels.free]: 'free' }).name('Speed control'), 'speedMode');
    bind(fEng.add(p, 'rpm', ui.rpmRange[0], ui.rpmRange[1], 10).name('Speed [rpm]'), 'rpm');
    if (ui.loadModels.length > 1) bind(fEng.add(p, 'loadModel', loadModelOptions(ui)).name('Load'), 'loadModel');
    bind(fEng.add(p, 'loadTorque', ui.loadRange[0], ui.loadRange[1], 0.5).name('Load torque [N·m]'), 'loadTorque');
    if (ui.loadModels.includes('brake')) {
      bind(fEng.add(p, 'brakeRefRpm', ui.rpmRange[0], ui.rpmRange[1], 10).name('Brake ref. speed [rpm]'), 'brakeRefRpm').domElement.title =
        'The brake absorbs the load torque at this speed: T = T_load · (n / n_ref)^k';
      bind(fEng.add(p, 'brakeExponent', 0, 3, 0.1).name('Brake exponent k'), 'brakeExponent').domElement.title =
        'k = 0: constant torque; k = 2: fan or water brake';
    }
    if (ui.loadModels.includes('vehicle') && spec.vehicle) {
      bind(fEng.add(p, 'gear', gearOptions(spec)).name('Gear'), 'gear');
      bind(fEng.add(p, 'gradePct', -20, 20, 0.5).name('Road grade [%]'), 'gradePct');
    }
    const throttle = bind(fEng.add(p, 'throttlePct', 0, 100, 1).name(`${ui.throttleLabel} [%]`), 'throttlePct');
    throttle.domElement.title = this.cfrPresets ? CFR_THROTTLE_HINT : THROTTLE_HINT;
    bind(fEng.add(p, 'coolantC', 20, 130, 1).name('Coolant [°C]'), 'coolantC');

    const fIgn = this.gui.addFolder('Ignition');
    bind(fIgn.add(p, 'sparkAdvanceDeg', ui.sparkRange[0], ui.sparkRange[1], 0.5).name(`${ui.sparkLabel} [° BTDC]`), 'sparkAdvanceDeg');
    if (vis.dwell) bind(fIgn.add(p, 'dwellMs', 0.5, 10, 0.1).name('Dwell [ms]'), 'dwellMs');
    if (vis.ignitionSource) {
      bind(fIgn.add(p, 'ignitionSource', { 'Magneto (MAG)': 'magneto', 'Battery (BAT)': 'battery' }).name('Ignition switch'), 'ignitionSource').domElement.title =
        'Coil supply: the flywheel magneto (output rises with speed) or the 6 V battery';
    }

    const fFuel = this.gui.addFolder('Fuel & mixture');
    bind(fFuel.add(p, 'fuel', FUEL_OPTIONS).name('Fuel'), 'fuel');
    bind(fFuel.add(p, 'octaneNumber', 0, 100, 0.5).name(this.cfrPresets ? 'PRF octane no.' : ui.fuelLabel), 'octaneNumber');
    bind(fFuel.add(p, 'equivalenceRatio', 0.4, 2.0, 0.01).name('Equivalence ratio φ'), 'equivalenceRatio');
    bind(fFuel.add(p, 'egrPct', 0, 40, 0.5).name('EGR [%]'), 'egrPct');

    const fAir = this.gui.addFolder('Intake & ambient');
    bind(fAir.add(p, 'intakeMixtureC', 0, 200, 1).name('Mixture temp. [°C]'), 'intakeMixtureC');
    bind(fAir.add(p, 'ambientKPa', 50, 110, 0.1).name('Ambient p [kPa]'), 'ambientKPa');
    bind(fAir.add(p, 'ambientC', -30, 50, 0.5).name('Ambient T [°C]'), 'ambientC');
    bind(fAir.add(p, 'humidityPct', 0, 100, 1).name('Rel. humidity [%]'), 'humidityPct');
    fAir.close();

    // --- presets ---
    if (this.cfrPresets) {
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
    } else if (def.presets.length) {
      const fPre = this.gui.addFolder('Presets');
      for (const pr of def.presets) {
        const c = fPre.add({ apply: () => this.applyEnginePreset(pr) }, 'apply').name(pr.label);
        if (pr.note) c.domElement.title = pr.note;
      }
    }

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
      this.monCtrl?.updateDisplay();
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

  /** CFR rating preset (RON / MON). */
  applyPreset(id: PresetId): void {
    const info = getPreset(id, this.params.compressionRatio);
    this.presetParams.monSchedule = info.sparkFollowsCR;
    this.presetParams.source = info.source === 'cfr.ts' ? 'engines/cfr.ts' : 'built-in (ASTM, approx.)';
    this.sourceCtrl?.show();
    this.monCtrl?.updateDisplay();
    this.sourceCtrl?.updateDisplay();
    this.setOperatingPoint(mergeOp(this.op, info.patch));
    this.cb.onOperatingPointChange(info.patch);
  }

  /** A preset of the engine registry (fixed-CR engines never get a CR from it). */
  applyEnginePreset(pr: EnginePreset): void {
    const patch: Partial<OperatingPoint> = { ...pr.op };
    if (!controlVisibility(this.def, this.params).compressionRatio) delete patch.compressionRatio;
    this.setOperatingPoint(mergeOp(this.op, patch));
    this.cb.onOperatingPointChange(patch);
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
    const vis = controlVisibility(this.def, p);
    const show = (key: OpParamKey, on: boolean): void => {
      this.ctrls.get(key)?.show(on);
    };
    show('octaneNumber', vis.octaneNumber);
    show('loadModel', vis.loadModel);
    show('loadTorque', vis.loadTorque);
    show('brakeRefRpm', vis.brake);
    show('brakeExponent', vis.brake);
    show('gear', vis.gear);
    show('gradePct', vis.grade);
    this.ctrls.get('loadTorque')?.name(vis.brake ? 'Brake torque at n_ref [N·m]' : 'Load torque [N·m]');
    this.ctrls.get('rpm')?.name(p.speedMode === 'free' ? 'Initial speed [rpm]' : 'Speed [rpm]');
  }

  dispose(): void {
    this.gui.destroy();
  }
}
