/**
 * Validation round 1 — adversarial code review of src/physics/cycle/** (cycle model, EngineSimulator,
 * worker swap). Every test here demonstrates a defect found in review; they are EXPECTED TO FAIL until
 * the fixer addresses them (each `it` names the defect and the file/line of its cause).
 *
 * Status after the round-1 fixer pass: all pass. Tests whose scenario needs the round-1 (uncalibrated)
 * behaviour pin those options explicitly (the defaults are now the CFR calibration).
 *
 * Runtime: ≈ 20–40 s total on an idle machine (each case runs 1–3 cycles of the full model).
 * Scratch probes that produced the numbers quoted below live in tmp/cr/ (not part of the suite).
 */
import { describe, expect, it } from 'vitest';
import type { EngineSnapshot } from '../../src/physics/core/snapshot';
import { EL, FUEL_SPECIES, NS, SP } from '../../src/physics/core/species';
import { R_UNIVERSAL } from '../../src/physics/core/constants';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { CycleModel, EngineSimulator, MODE_OPEN } from '../../src/physics/cycle';
import { I_LV, I_MB, I_ME, I_MU } from '../../src/physics/cycle/cycle-model';
import { laminarFlameSpeed } from '../../src/physics/combustion';
import { mixMolarMass } from '../../src/physics/thermo';
import { MOLAR_MASS } from '../../src/physics/thermo/thermo';

/** Collect the snapshots of cycle 0 (after the warm-up) of an EngineSimulator. */
function cycleSnapshots(sim: EngineSimulator): EngineSnapshot[] {
  const out: EngineSnapshot[] = [];
  let s = sim.advanceToNextSnapshot();
  while (s.cycle < 1) {
    out.push(s);
    s = sim.advanceToNextSnapshot();
  }
  return out;
}

/** Round-2 model additions switched off (the round-1 physics these scenarios were written for). */
const ROUND2_OFF = {
  knockBrushAutoignition: false,
  knockSourceShells: 0,
  mapoBand: null,
  creviceModel: false,
  wallTemperatureModel: 'fixed',
  burnedNOThermo: 'equilibrium',
  kernelHandoffMultiple: 1,
  knockStratificationDT: 15,
} as const;
/** CFR F-1 with the round-1 1 L intake plenum. */
const SPEC_1L = { ...CFR_F1, manifolds: { ...CFR_F1.manifolds, intakeVolume: 1.0e-3 } };

describe('code review: knock oscillator NaN after the end-gas burn-out merge', () => {
  // cycle-model.ts afterStep(): the burn-out merge (step 3) switches to MODE_BURNED and refresh()
  // sets this.Tu = 0 and closure.Tu = 0; step 5 then evaluates
  //   qc = cl.heatOfReaction(this.Tu > 0 ? this.Tu : cl.Tu) = heatOfReaction(0) = NaN
  // with dmKnock > 0 (captured before the merge), so KnockOscillator η, η̇ become NaN until the EVO
  // reset. MAPO survives (NaN comparisons are false), but EngineSnapshot.pressure, knock.oscillation,
  // gasTorque and netTorque are NaN for the rest of the expansion (and CycleTrace.pressureReported).
  it('RON with the recommended burn-rate multiplier 2.6: every snapshot pressure is finite', () => {
    const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, {
      snapshotEveryDeg: 1,
      bufferAheadSeconds: 0.1,
      warmupCycles: 1,
      burnRateMultiplier: 2.6,
    });
    const snaps = cycleSnapshots(sim);
    const bad = snaps.filter((s) => !Number.isFinite(s.pressure) || !Number.isFinite(s.knock.oscillation));
    // measured: 262 of the 0.5°-cadence snapshots (θ = 32.5° … EVO) are NaN
    expect(bad.length, `first NaN at θ = ${bad[0]?.thetaDeg}`).toBe(0);
  });
});

describe('code review: crevice volume is outside the flame geometry — the charge never burns out', () => {
  // The crevice volume (1.5 cm³) is part of V (SliderCrank) and therefore of the unburned zone, but
  // the entrained volume V_e is mapped onto the disc of height h only (FlameGeometry.radiusForVolume
  // in evalClosed). V_e can never exceed A·h = V − V_crevice, so the last ≈ V_crev·ρ_u of unburned
  // gas is never entrained: A_f → 0.2–0.4 cm², the brush empties and m_u stalls at ≈ 0.8 % of the
  // charge. The burn-out merge (1e-5) never fires, burnDone stays false (x_b(EVO) = 0.9917 < 0.999).
  // DESIGN/known-issues say "crevice gas burns with the charge"; with creviceVolume = 0 the same cycle
  // reaches x_b = 1.000000 (measured), which isolates the cause.
  it('a non-knocking RON cycle burns the whole charge before EVO (x_b(EVO) ≥ 0.999)', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1, knock: false });
    const s = m.runCycles(1)[0];
    const xbEvo = s.fuelMass / m.fuelMassIvc;
    expect(xbEvo).toBeGreaterThan(0.999); // measured 0.99175
  });

  it('the flame of a non-knocking cycle reaches burnout/done before EVO and x_b is kept after EVO', () => {
    const sim = new EngineSimulator(CFR_F1, CFR_RON_CONDITIONS, {
      snapshotEveryDeg: 2,
      bufferAheadSeconds: 0.1,
      warmupCycles: 1,
      knock: false,
    });
    const snaps = cycleSnapshots(sim);
    const closedLate = snaps.filter((s) => s.thetaDeg > 100 && s.thetaDeg < 140);
    // measured: stage 'turbulent', phase 'combustion' up to 140°, flame radius 146 mm, A_f 0.3 cm²
    expect(closedLate.some((s) => s.flame.stage === 'burnout' || s.flame.stage === 'done')).toBe(true);
    const afterEvo = snaps.filter((s) => s.thetaDeg > 145 && s.thetaDeg < 300);
    // measured: massFractionBurned = 0 and flame stage 'none' for the whole exhaust stroke
    expect(afterEvo.every((s) => s.massFractionBurned > 0.9)).toBe(true);
  });

  it('MON: no "knock onset" is reported on the unburnable crevice remainder late in expansion', () => {
    // measured: knockOnsetDeg = 95.9°, knockEndGasFraction = 1.05 %, MAPO = 0 — the Livengood–Wu
    // integral keeps running on the ≈1 % leftover until it "autoignites" in the expansion stroke.
    const m = new CycleModel(CFR_F1, CFR_MON_CONDITIONS, { warmupCycles: 2 });
    const s = m.runCycles(1)[0];
    if (!Number.isNaN(s.knockOnsetDeg)) {
      expect(s.knockEndGasFraction).toBeGreaterThan(0.02);
      expect(s.mapo).toBeGreaterThan(0);
    }
  });
});

describe('code review: volumetric efficiency counts retained unburned charge as inducted air', () => {
  // cycle-model.ts onIvc(): fresh = m·(1 − y_res) with y_res the BURNED-gas scalar. When the gas
  // retained from the previous cycle is unburned (motored engine, misfire, partial burn) it is counted
  // as freshly inducted air → η_v > 1. Measured (RON conditions, 8 warm-up cycles):
  //   motored:        summary 1.0894 vs inducted (venturi ledger) 0.9324
  //   misfire φ 0.6:  summary 1.0931 vs inducted 0.9370
  //   fired:          summary 0.9607 vs inducted 0.9606 (agrees only because y_res is the residual)
  function inductedEv(m: CycleModel, run: () => void): number {
    const a = new Float64Array(NS);
    for (let k = 0; k < NS; k++) a[k] = m.y[I_LV + k];
    run();
    let mTot = 0;
    for (let k = 0; k < NS; k++) mTot += (m.y[I_LV + k] - a[k]) * MOLAR_MASS[k];
    let mNonAir = (m.y[I_LV + SP.H2O] - a[SP.H2O]) * MOLAR_MASS[SP.H2O];
    for (const f of FUEL_SPECIES) mNonAir += (m.y[I_LV + SP[f]] - a[SP[f]]) * MOLAR_MASS[SP[f]];
    const op = m.op;
    const xw = m.airX[SP.H2O];
    const rhoDry = ((mixMolarMass(m.airX) - xw * MOLAR_MASS[SP.H2O]) * op.ambientPressure) / (R_UNIVERSAL * op.ambientTemperature);
    return (mTot - mNonAir) / (rhoDry * m.kin.displacedVolume);
  }

  it('motored engine: summary η_v equals the dry air actually inducted per cycle', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { combustionModel: 'none', warmupCycles: 6 });
    const sums: { volumetricEfficiency: number }[] = [];
    const ev = inductedEv(m, () => {
      sums.push(...m.runCycles(1));
    });
    const s = sums[0];
    expect(s.volumetricEfficiency).toBeLessThan(1);
    expect(Math.abs(s.volumetricEfficiency / ev - 1)).toBeLessThan(0.01);
  });
});

describe('code review: combustion sub-models use the REQUESTED φ / fuel, not the trapped mixture', () => {
  // setOperatingPoint(φ, fuel) takes effect at the next cycle start, but the 1 L intake plenum still
  // holds ≈ 1.6 cycles of the old mixture, so the trapped charge changes over ≈ 5 cycles (modelled,
  // physical). evalClosed / ignitionSplit / tauNow / onIvc nevertheless evaluate S_L, Markstein and
  // Lewis numbers, the ignition delay and xResMole with this.op.equivalenceRatio and this.fuel.
  // Measured after {φ 1.1 → 0.8, PRF 90 → 100}: trapped φ (element balance) 1.039 / 0.938 / 0.879 /
  // 0.846 / 0.828 and iso-octane share 90.7 / 94.3 / 96.6 / 98.0 / 98.9 % in cycles 1–5, while S_L and
  // τ use φ = 0.8 and PRF 100 from cycle 1 (CA50 jumps from 36° to 52.8°).
  it('first IVC after a φ change: S_L is evaluated at the trapped equivalence ratio', () => {
    // (fixer round 2: the 1 L intake plenum this scenario relies on is pinned — CFR_F1 now has 0.25 L)
    const m = new CycleModel({ ...CFR_F1, manifolds: { ...CFR_F1.manifolds, intakeVolume: 1.0e-3 } }, CFR_RON_CONDITIONS, { warmupCycles: 2 });
    m.setOperatingPoint({ equivalenceRatio: 0.8 });
    m.runCycles(1); // finish the current cycle; the change applies at the wrap
    while (m.mode === MODE_OPEN || m.theta < -150) m.stepUntil(Infinity, -150);
    const b = m.closure.bu;
    // element-based φ of the trapped charge (humid-air H2O biases it by ≈ −0.2 %)
    const phiTrapped = (2 * b[EL.C] + 0.5 * b[EL.H]) / b[EL.O];
    expect(phiTrapped).toBeGreaterThan(0.95); // the charge is still near the old φ = 1.1 …
    const slTrapped = laminarFlameSpeed(m.fuel, phiTrapped, m.Tu, m.p, m.xDil);
    // … and the model's S_L must be the trapped mixture's (measured: S_L(0.8) ≈ 0.7 × S_L(1.04))
    expect(Math.abs(m.SL / slTrapped - 1)).toBeLessThan(0.03);
  });
});

describe('code review: NaN in snapshots for fuels without a delay model', () => {
  // The default knock model (LLNL PRF tables) returns τ = NaN for CH4 / C3H8 / C2H5OH (documented in
  // chemistry/index.ts); the cycle integrates it, so EngineSnapshot.knock.integral is NaN in every
  // closed-phase snapshot (the snapshot contract is "plain JSON-able", NaN is not). With
  // ignitionDelayModel 'douaud-eyzat' the same fuel throws RangeError out of tauNow() at the first IVC.
  it('methane: every snapshot is finite', () => {
    const sim = new EngineSimulator(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'pure', species: 'CH4' } }, {
      snapshotEveryDeg: 5,
      bufferAheadSeconds: 0.1,
      warmupCycles: 0,
    });
    const snaps = cycleSnapshots(sim);
    expect(snaps.filter((s) => !Number.isFinite(s.knock.integral)).length).toBe(0);
  });

  it("methane with ignitionDelayModel 'douaud-eyzat' does not throw", () => {
    expect(() => new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'pure', species: 'CH4' } }, { warmupCycles: 1, ignitionDelayModel: 'douaud-eyzat' })).not.toThrow();
  });
});

describe('code review: knock acoustic source vanishes when the flame circle covers the head plane', () => {
  // onAutoignition(): knockOsc.setEndGasRegion(gap.x, gap.z, flameRadius()) uses the FULL sphere
  // radius r_f as the planform circle (its head-plane section). Late in the burn r_f exceeds the plug-
  // to-far-wall distance while end gas remains near the piston: the planform end-gas area is 0 →
  // sourceShape = 0 → no excitation at all. Measured, RON conditions with PRF 95: onset 74.4° with
  // 4.9 % of the charge (1.4 % of V) autoigniting, planform end-gas area 0 % (22 % would lie outside
  // the flame at the piston plane), MAPO = 0.000 Pa. At PRF 90 the area used is 5.6 % vs an end-gas
  // volume fraction of 11.7 %.
  it('an autoignition of ≥ 1 % of the charge produces a non-zero MAPO', () => {
    // The scenario (late, expansion-stroke autoignition with the front radius beyond the head plane)
    // needs the slow round-1 burn; since the fixer pass the defaults are the calibrated CFR set, under
    // which PRF 95 does not knock at CR 6.43 at all — so the round-1 (uncalibrated) options are pinned.
    // (fixer round 2: the round-2 model additions are switched off too, with the 1 L plenum of round 1)
    const round1 = { burnRateMultiplier: 1, taylorScaleMultiplier: 1, woschniMultiplier: 1, intakePortHeatTransferMultiplier: 0, venturiDischargeCoefficient: 0.984, ignitionDelayModel: 'llnl-gasoline-2011' as const, ...ROUND2_OFF };
    const m = new CycleModel(SPEC_1L, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 95 } }, { warmupCycles: 2, ...round1 });
    const s = m.runCycles(1)[0];
    expect(Number.isNaN(s.knockOnsetDeg)).toBe(false);
    expect(s.knockEndGasFraction).toBeGreaterThan(0.01);
    expect(s.mapo).toBeGreaterThan(1000); // > 0.01 bar for a ≈ 100 J volumetric release
  });
});

describe('code review: CycleSummary.maxPressureRiseRate depends on the integration step', () => {
  // bookkeeping(): (p − p_ref)/(θ − θ_ref) over windows that close at the first step end ≥ 0.1° after
  // the previous one, i.e. windows of 0.1° … 0.1° + h. Same cycle, different maxStepDeg (measured,
  // cycle 2 at RON): 0.569 bar/° (0.25°), 0.662 bar/° (0.1°), 0.640 bar/° (0.05°); with burn-rate
  // multiplier 2.6: 5.08 / 4.78 / 5.34 bar/°. A summary metric should not move 15–20 % with the
  // step. (The analytic ṗ = closure.rates() is available at every evaluation.)
  it('RON: maxPressureRiseRate agrees within 5 % between maxStepDeg 0.25 and 0.1', () => {
    const a = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 2, maxStepDeg: 0.25 }).runCycles(1)[0];
    const b = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 2, maxStepDeg: 0.1, fineStepDeg: 0.05 }).runCycles(1)[0];
    // (IMEP, CA50 and the knock onset agree to < 0.1 % / 0.02° between these runs)
    expect(Math.abs(a.imepNet / b.imepNet - 1)).toBeLessThan(1e-3);
    expect(Math.abs(a.maxPressureRiseRate / b.maxPressureRiseRate - 1)).toBeLessThan(0.05);
  });
});

describe('code review (sanity): burned-gas scalar and zone bookkeeping', () => {
  // Round 1 recorded the then-documented deviation knockEndGasFraction = m_u/m (flame brush m_e − m_b
  // included, ≈ 2.2 % of the charge at the RON onset). Fixed in the round-1 fixer pass (knock
  // validator finding 8 / code-review finding 9): the reported end-gas fraction is the unburned gas
  // AHEAD of the front, (m − m_e)/m, and the autoignition burn-up consumes only that gas. This test
  // now pins the corrected definition (expectation updated deliberately).
  it('the knock end-gas fraction is the unburned mass ahead of the front, (m − m_e)/m at onset', () => {
    // (round-1 uncalibrated options pinned: this RON preset only knocks with the slow round-1 burn;
    // fixer round 2: the round-2 additions are switched off as well, with knockBrushAutoignition false —
    // the definition this test pins; the round-2 default is covered below)
    const round1 = { burnRateMultiplier: 1, taylorScaleMultiplier: 1, woschniMultiplier: 1, intakePortHeatTransferMultiplier: 0, venturiDischargeCoefficient: 0.984, ignitionDelayModel: 'llnl-gasoline-2011' as const, ...ROUND2_OFF };
    const m = new CycleModel(SPEC_1L, CFR_RON_CONDITIONS, { warmupCycles: 1, ...round1 });
    while (!m.knockOnset && m.theta < 140) m.stepUntil(Infinity, 140);
    expect(m.knockOnset).toBe(true);
    const mt = m.y[I_MU] + m.y[I_MB];
    expect(m.knockEndGasFraction).toBeCloseTo((mt - m.y[I_ME]) / mt, 12);
    expect(m.knockEndGasFraction).toBeLessThan(m.y[I_MU] / mt);
  });

  it('fixer round 2 default (knockBrushAutoignition): the autoigniting fraction is all unburned gas, m_u/m at onset', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: 6.6 }, { warmupCycles: 1 });
    while (!m.knockOnset && m.theta < 140) m.stepUntil(Infinity, 140);
    expect(m.knockOnset).toBe(true);
    const mt = m.y[I_MU] + m.y[I_MB];
    expect(m.knockEndGasFraction).toBeCloseTo(m.y[I_MU] / mt, 12);
    expect(m.knockEndGasFraction).toBeGreaterThan((mt - m.y[I_ME]) / mt);
  });
});
