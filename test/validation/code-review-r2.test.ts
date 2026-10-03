/**
 * Validation round 2 — adversarial code review of src/physics/cycle/** (cycle model, EngineSimulator,
 * worker swap) after the round-1 fixer pass. Every test here demonstrates a defect found in review
 * and is EXPECTED TO FAIL until the fixer addresses it (each `describe` names the cause, file/line).
 * The measured values quoted in the comments were taken on the shipped defaults (CFR_CALIBRATION).
 *
 * Runtime ≈ 7 s on an idle machine (each case runs 1–2 cycles of the full model after its warm-up).
 * Scratch probes that produced the numbers live in tmp/cr2/ (not part of the suite).
 */
import { describe, expect, it } from 'vitest';
import { DEG } from '../../src/physics/core/constants';
import type { EngineSnapshot } from '../../src/physics/core/snapshot';
import { FUEL_SPECIES, NS, SP } from '../../src/physics/core/species';
import { CFR_F1, CFR_MON_CONDITIONS, CFR_RON_CONDITIONS } from '../../src/physics/engines/cfr';
import { CycleModel, EngineSimulator, MODE_OPEN } from '../../src/physics/cycle';
import { I_CN, I_CU, I_IN, I_IU } from '../../src/physics/cycle/cycle-model';
import { newCrankTrainState } from '../../src/physics/mechanics';
import { KnockOscillator } from '../../src/physics/chemistry/knock';
import { completeCombustionProducts, freshCharge } from '../../src/physics/thermo/fuels';
import { mixSMolar, temperatureFromUMolar } from '../../src/physics/thermo';
import { MOLAR_MASS } from '../../src/physics/thermo/thermo';

/** Mass fraction of CO2 of a mole-fraction vector. */
function yCO2(X: Float64Array): number {
  let m = 0;
  for (let k = 0; k < NS; k++) m += X[k] * MOLAR_MASS[k];
  return (X[SP.CO2] * MOLAR_MASS[SP.CO2]) / m;
}

describe('knock: MAPO is exactly 0 for a detected knock with a small end gas (quadrature of the source region)', () => {
  // cycle-model.ts onAutoignition() (≈ line 2129) passes a planform circle whose complement has the
  // end-gas VOLUME fraction to KnockOscillator.setEndGasRegion (chemistry/knock.ts ≈ line 326), which
  // integrates on a fixed 48 × 96 midpoint grid. When the crescent is thinner than the outermost cell
  // midpoints (≈ 0.43 mm from the liner) NO cell is counted: aEg = 0 → sourceShape = 0 → MAPO = 0.
  // Measured at the RON guide-table point (PRF 90, CR 6.43): onset 13.07°, end gas 0.250 % of the charge,
  // planform fraction 0.0000 % (requested ≈ 0.083 %), MAPO 0.000 bar; with a 480 × 960 grid the same
  // cycle gives 0.116 bar. The planform fraction is quantised in 0.043 % steps (0.0859, 0.1718, 0.2577 %
  // at CR 6.46/6.50/6.55 vs 0.119/0.182/0.287 % resolved), so MAPO starts from exactly 0 with a jump.
  it('RON guide-table point: a knock onset with ≥ 0.2 % end gas produces a non-zero MAPO', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS);
    const s = m.runCycles(1)[0];
    expect(Number.isNaN(s.knockOnsetDeg)).toBe(false);
    expect(s.knockEndGasFraction).toBeGreaterThan(2e-3);
    expect(m.knockOsc.endGasAreaFraction, 'planform end-gas area fraction').toBeGreaterThan(0);
    expect(s.mapo).toBeGreaterThan(0.02e5); // resolved quadrature: 0.116 bar
  });
});

describe('knock: CycleSummary.mapo includes acoustic modes far above the measured MAPO band', () => {
  // KnockOscillator is built by the cycle (cycle-model.ts constructor ≈ line 663) with the default
  // maxAlpha 7.1: 16 modes up to ≈ 25 kHz at the end-gas sound speed (c ≈ 925 m/s), and
  // CycleSummary.mapo = max |Σ ψ_j η_j| over ALL of them (index.ts summary definitions). The MAPO the
  // model is calibrated/validated against is band-passed: 4–18 kHz (Hoth 2021, SON 2023, critical-CR
  // data, the 1.43 bar calibration target of knockStratificationDT), 6–20 kHz (Rockstroh 2018, the
  // 0.67 bar standard-knock definition) — test/fixtures/cfr_validation_README.md; the 0.1° indicated
  // data at 600 rpm cannot even resolve > 18 kHz. Measured (RON, PRF 90): MAPO with all modes vs modes
  // ≤ 17.5 kHz = 0.644 vs 0.420 bar (CR 6.6), 1.220 vs 0.832 (CR 6.7), 3.31 vs 2.48 (CR 7.0) — the
  // out-of-band modes add 30–55 % to the metric that defines standard knock. At the calibration point
  // (Choi/Pal PRF 100, CR 7.55, spark −12.72°) the shipped ΔT 15 K gives 1.50 bar with all modes (the
  // "match" to Hoth's 1.43 bar) but 1.05 bar in band; ΔT ≈ 10.5 K restores 1.43 bar in band.
  it('RON CR 6.7: MAPO equals the MAPO of the modes below 18 kHz within 5 %', () => {
    const op = { ...CFR_RON_CONDITIONS, compressionRatio: 6.7 };
    const full = new CycleModel(CFR_F1, op).runCycles(1)[0].mapo;
    const m = new CycleModel(CFR_F1, op);
    // α_max for 18 kHz at c ≈ 925 m/s: 18e3·π·B/c ≈ 5.05 → keep (1,0) (2,0) (0,1) (3,0) (both orientations)
    (m as unknown as { knockOsc: KnockOscillator }).knockOsc = new KnockOscillator(CFR_F1.geometry.bore, {
      decayTime: m.opts.knockDecayTime,
      sensor: m.opts.knockSensor,
      maxAlpha: 4.9,
    });
    const band = m.runCycles(1)[0].mapo;
    expect(band).toBeGreaterThan(0.5e5);
    expect(Math.abs(full / band - 1)).toBeLessThan(0.05); // measured 1.47
  });
});

describe('knock: knock intensity is non-monotonic in the octane number (uniform-source limit)', () => {
  // The modal source strength is ⟨ψ_j⟩ over the end-gas planform region (knock.ts), which → 0 as the
  // region → the whole bore (every Neumann mode has zero disc mean). onAutoignition() makes the region
  // the bore minus a circle about the plug with the end-gas volume fraction (MODE_SINGLE: the whole
  // bore). Heavier knock (more end gas) therefore rings LESS. Measured at RON conditions, CR 6.43:
  //   PRF 80 / 70 / 60 / 50 / 40: end gas 10 / 27 / 49 / 77 / 100 %, MAPO 4.18 / 7.55 / 7.49 / 3.84 / 0.00 bar
  // (PRF 40: the whole charge autoignites 2° after the spark, max dp/dθ 17.6 bar/° — yet MAPO 0).
  // A knockmeter on these fuels pegs; the model reports "no knock" for the most violent case.
  it('RON CR 6.43: MAPO(PRF 40) ≥ MAPO(PRF 70)', () => {
    const mapo = (on: number): number =>
      new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: on } }).runCycles(2)[1].mapo;
    const m70 = mapo(70);
    const m40 = mapo(40);
    expect(m70).toBeGreaterThan(1e5);
    expect(m40).toBeGreaterThanOrEqual(m70);
  });
});

describe('knock: low-octane PRFs "autoignite" long before the spark (Douaud–Eyzat τ → 0 as ON → 0)', () => {
  // chemistry/ignition-delay.ts douaudEyzatTau ∝ (ON/100)^3.402 → 0 for n-heptane-rich blends; the cycle
  // arms the Livengood–Wu integral from IVC (cycle-model.ts lwArmed) with the default model
  // (calibration.ts CFR_KNOCK_DELAY_MODEL), so the whole charge "autoignites" during compression at
  // T_u ≈ 370–450 K. Measured at RON conditions: PRF 0 onset −151.8° (at IVC), PRF 2 −150.0°, PRF 5
  // −121.2°, PRF 10 −74.6°, PRF 20 −41.5° — whole charge, MAPO 0, IMEP −6.45 / −6.44 / −5.60 / −1.99 /
  // +2.18 bar. No hydrocarbon autoignites in milliseconds at ≈ 400 K. The sanitiser accepts ON 0..100.
  it('PRF 5 at RON conditions: no end-gas autoignition before −60° (T_u < 450 K there)', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, fuel: { kind: 'PRF', octaneNumber: 5 } });
    const s = m.runCycles(2)[1];
    if (!Number.isNaN(s.knockOnsetDeg)) expect(s.knockOnsetDeg).toBeGreaterThan(-60);
    expect(s.imepNet).toBeGreaterThan(0);
  });
});

describe('dilution: x_dil uses the REQUESTED EGR fraction, not the trapped charge', () => {
  // cycle-model.ts onIvc() ≈ line 2659: xDil = yRes + (1 − yRes)·this.op.egrFraction. The fixer moved
  // φ and the fuel blend to the TRAPPED mixture (DESIGN.md: "Combustion sub-models … use the TRAPPED
  // mixture: a fuel or φ change reaches the cylinder over ≈ 5 cycles through the 1 L intake plenum")
  // but EGR still jumps at the first IVC. Measured after {egrFraction 0 → 0.2} (RON): first IVC
  // xDil 0.266 vs trapped products fraction 0.121 (from the trapped CO2), S_L 0.184 m/s, and that cycle
  // misfires; the trapped dilution only reaches the model's value after ≈ 7 cycles (0.121, 0.114,
  // 0.217, 0.181, 0.252, 0.202, 0.263 vs xDil 0.266, 0.212, 0.262, 0.211, 0.264, 0.211, 0.264).
  it('first IVC after an EGR step: xDil agrees with the trapped products fraction within 15 %', () => {
    // (fixer round 2: the premise "the plenum still delivers mostly the old charge" is the 1 L plenum of
    // round 1; CFR_F1 now has 0.25 L (Choi 2018 Fig. 2 intake-port pressure), so the 1 L plenum is pinned)
    const m = new CycleModel({ ...CFR_F1, manifolds: { ...CFR_F1.manifolds, intakeVolume: 1.0e-3 } }, CFR_RON_CONDITIONS, { warmupCycles: 2 });
    const fa = freshCharge({ fuel: m.fuel, phi: CFR_RON_CONDITIONS.equivalenceRatio, airX: m.airX });
    const yProd = yCO2(completeCombustionProducts(fa));
    const yAir = yCO2(fa);
    m.setOperatingPoint({ egrFraction: 0.2 });
    m.runCycles(1); // the change applies at this wrap
    while (m.mode === MODE_OPEN || m.theta < -150) m.stepUntil(Infinity, -150);
    const trapped = (yCO2(m.closure.Xu) - yAir) / (yProd - yAir);
    expect(trapped).toBeLessThan(0.15); // the plenum still delivers mostly the old charge
    expect(Math.abs(m.xDil / trapped - 1)).toBeLessThan(0.15);
  });
});

describe('gas exchange: the carburettor venturi sees the heated vaporised MIXTURE instead of the intake air', () => {
  // cycle-model.ts prepareCycleData() ≈ line 1008 builds `fresh` at op.intakeMixtureTemperature with the
  // fuel fully vaporised, and evaluate() ≈ line 1321 flows it through the 9/16 in venturi. On the CFR the
  // carburettor meters AIR at the intake-air temperature (RON: heated air 52 °C; MON: 38 °C air) and the
  // MON mixture heater (149 °C) is DOWNSTREAM of it (cfr.ts CFR_MON_CONDITIONS doc, D2700 §10.3.6–7).
  // At the same Δp the model's venturi passes 18 % less mass in MON (vaporised mixture at 422 K, M ≈ 30.5)
  // than air at 311 K + metered fuel, 0.8 % less in RON — so the ONE global venturi C_D (0.65,
  // calibration.ts, fitted to the MON venturi-size effect and the MON/RON compression-pressure ratio)
  // absorbs a MON-only placement error. Physically the heater setting cannot change the air flow the
  // carburettor passes at a given depression; in the model it changes it by 14 %.
  function airFlowAt(m: CycleModel, pPlenum: number): number {
    m.intake.setTPX(m.op.intakeMixtureTemperature, pPlenum, m.Xfresh);
    for (let k = 0; k < NS; k++) m.y[I_IN + k] = m.intake.N[k];
    m.y[I_IU] = m.intake.U;
    m.refresh();
    let mt = 0;
    let mf = 0;
    for (let k = 0; k < NS; k++) mt += m.Xfresh[k] * MOLAR_MASS[k];
    for (const s of FUEL_SPECIES) mf += m.Xfresh[SP[s]] * MOLAR_MASS[SP[s]];
    return m.mdotV * (1 - mf / mt); // air (with its humidity) through the carburettor, kg/s
  }
  it('MON: the carburettor air flow at a given depression does not depend on the mixture-heater setting', () => {
    const hot = new CycleModel(CFR_F1, CFR_MON_CONDITIONS, { warmupCycles: 0 });
    const cold = new CycleModel(CFR_F1, { ...CFR_MON_CONDITIONS, intakeMixtureTemperature: CFR_MON_CONDITIONS.ambientTemperature }, { warmupCycles: 0 });
    const p = CFR_MON_CONDITIONS.ambientPressure - 2000;
    const ratio = airFlowAt(hot, p) / airFlowAt(cold, p);
    // measured 0.858
    expect(Math.abs(ratio - 1)).toBeLessThan(0.02);
  });
});

describe('snapshot: gasTorque includes the synthesised knock oscillation at the PICKUP', () => {
  // engine-simulator.ts makeSnapshot(): m.dyn.evaluate(θ, ω, α, pressure, …) with pressure = p + p_osc(sensor).
  // The acoustic field is a sum of rigid-wall (Neumann) modes J_m(α r/R)·{cos, sin}(mθ) whose mean over
  // the piston face is exactly zero (m ≥ 1: angular integral; m = 0: ∫J₀(α₀ₙr/R) r dr ∝ J₁(α₀ₙ) = 0), so
  // it exerts no net force on the piston; the free-speed dynamics (evaluate(), this.p) correctly ignore
  // it, but the reported gas/net torque ring with ±A_p·p_osc·dx/dθ (measured up to 8.8 N m at MAPO 1.2 bar, on a mean gas torque of ≈ 40 N m).
  it('RON CR 6.7 (MAPO ≈ 1.2 bar): snapshot gasTorque equals the torque of the thermodynamic pressure', () => {
    const sim = new EngineSimulator(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: 6.7 }, { snapshotEveryDeg: 0.5, bufferAheadSeconds: 1 });
    const m = sim.model;
    const cts = newCrankTrainState();
    const pcc = m.op.ambientPressure + m.opts.crankcaseGaugePressure;
    let worst = 0;
    let ringing = 0;
    let s: EngineSnapshot = sim.advanceToNextSnapshot();
    while (s.cycle < 1) {
      if (s.knock.oscillation !== 0) {
        ringing++;
        const tq = m.dyn.evaluate(s.thetaDeg * DEG, (s.rpm * 2 * Math.PI) / 60, 0, s.pressure - s.knock.oscillation, pcc, cts).gasTorque;
        worst = Math.max(worst, Math.abs(s.gasTorque - tq));
      }
      s = sim.advanceToNextSnapshot();
    }
    expect(ringing).toBeGreaterThan(10);
    expect(worst).toBeLessThan(1e-6); // measured 8.76 N m
  });
});

describe('EVO: the kinetic-NO swap steps the charge temperature and pressure', () => {
  // cycle-model.ts onEvo() ≈ line 2760: the burned zone is carried as FULL equilibrium (equilibrium NO)
  // during the whole cycle, and at EVO swapNO() replaces its NO by the kinetic value at constant U. The
  // formation energy of the super-equilibrium NO (≈ 870 ppm at RON; +90 kJ/mol NO, and in the rich branch
  // CO2 → CO for the O atom, +283 kJ/mol) is debited all at once: in a fully burned RON cycle (x_b = 1)
  // T drops ≈ 10 K and p 0.82 % at the EVO instant (same θ, same U, same V), and the rich branch adds
  // ≈ 870 ppm CO to the exhaust. A state-preserving representation change should not move p.
  it('RON (fully burned): pressure is continuous across the EVO merge (< 0.1 %)', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS);
    m.recordTrace();
    m.runCycles(1);
    const tr = m.trace!;
    let k = -1;
    for (let i = 1; i < tr.theta.length; i++) if (tr.mode[i - 1] !== MODE_OPEN && tr.mode[i] === MODE_OPEN && Math.abs(tr.theta[i] - tr.theta[i - 1]) < 1e-9) k = i;
    expect(k).toBeGreaterThan(0);
    expect(tr.xb[k - 1]).toBeGreaterThan(0.999);
    const jump = tr.pressure[k] / tr.pressure[k - 1] - 1;
    expect(Math.abs(jump), `EVO pressure step ${(jump * 100).toFixed(3)} %`).toBeLessThan(1e-3);
  });
});

describe('summary: CA10/50/90 are NaN when the burn happens in an EVENT (instantaneous-at-tdc)', () => {
  // cycle-model.ts bookkeeping() ≈ line 2498 finds the x_b crossings inside a STEP from y0[I_MU]; the TDC
  // merge (handleEvents → EV_TDC → mergeToBurned) happens after bookkeeping at the step end, and the
  // next step starts at x_b = 1 — so no step ever contains the crossing and ca10/ca50/ca90 stay NaN
  // (CycleSummary: "NaN if not reached") although the whole charge burned at θ = 0.
  it("combustionModel 'instantaneous-at-tdc': CA10 = CA50 = CA90 = 0", () => {
    const s = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { combustionModel: 'instantaneous-at-tdc', warmupCycles: 1 }).runCycles(1)[0];
    expect(s.fuelMass).toBeGreaterThan(0);
    expect(s.ca10).toBeCloseTo(0, 6);
    expect(s.ca50).toBeCloseTo(0, 6);
    expect(s.ca90).toBeCloseTo(0, 6);
  });
});

describe("closure: the 'wiebe' start seed (1e-8 of the charge) makes the two-zone closure fail", () => {
  // cycle-model.ts handleEvents EV_WIEBE ≈ line 2597 seeds m_b = 1e-8·m_u (≈ 7e-12 kg); the first
  // closure solves on that zone fail (and the HP-flame-temperature retry too): 2–3 failures per cycle at
  // θ = −12.95° (RON, CR 9, knock on or off), counted in closureFailures, best iterate kept.
  it('wiebe, RON CR 9: no failed closure solve', () => {
    const m = new CycleModel(CFR_F1, { ...CFR_RON_CONDITIONS, compressionRatio: 9 }, { combustionModel: 'wiebe', knock: false, warmupCycles: 1 });
    const f0 = m.closureFailures;
    m.runCycles(1);
    expect(m.closureFailures - f0).toBe(0);
  });
});

describe('robustness: a NaN operating-point field corrupts the model before the sanitiser/friction throws', () => {
  // options.ts sanitizeOperatingPoint() clamps with `x < lo ? lo : x > hi ? hi : x`, which passes NaN.
  // CycleModel.setOperatingPoint({ rpm: NaN }) (≈ line 687) writes y[I_OM] = NaN, THEN updateFriction()
  // throws in pnhFmep — the state is left NaN (the worker host marks the simulator failed). φ = NaN is
  // accepted and throws at the next cycle start (freshCharge); spark advance NaN is accepted silently.
  it('setOperatingPoint({ rpm: NaN }) neither throws nor leaves a non-finite state', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 0 });
    let threw = false;
    try {
      m.setOperatingPoint({ rpm: Number.NaN });
    } catch {
      threw = true;
    }
    expect(threw).toBe(false);
    expect(Number.isFinite(m.rpm)).toBe(true);
  });
});

describe('compression-ratio change: the cylinder volume jumps at θ = −360 with no work (entropy decreases)', () => {
  // cycle-model.ts onWrap() ≈ line 2800: kin.setCompressionRatio() changes V_c instantaneously while the
  // cylinder state (N, U) is kept: CR 6.43 → 8 shrinks V by 22 % (112.65 → 87.39 cm³) at constant U and
  // T, p jumps 1.013 → 1.306 bar, and no p dV work is booked. For an adiabatic system with no work this
  // is an entropy DECREASE of R ln(1.289) per unit mole — a second-law violation (the real CFR moves the
  // cylinder slowly while running; either ramp V_c or book the compression isentropically).
  it('CR 6.43 → 8 at the wrap: the cylinder charge entropy does not decrease', () => {
    const m = new CycleModel(CFR_F1, CFR_RON_CONDITIONS, { warmupCycles: 1 });
    m.setOperatingPoint({ compressionRatio: 8 });
    const c0 = m.cycle;
    while (m.cycle === c0 && m.theta < 359.999) m.stepUntil(Infinity, 359.999);
    const entropy = (): number => {
      const N = m.y.slice(I_CN, I_CN + NS);
      let n = 0;
      for (let k = 0; k < NS; k++) n += N[k];
      const X = N.map((v) => v / n);
      const T = temperatureFromUMolar(N, m.y[I_CU], 500);
      return n * mixSMolar(X, T, m.p);
    };
    const s0 = entropy();
    const V0 = m.V;
    while (m.cycle === c0) m.stepUntil(Infinity, 360);
    m.refresh();
    const s1 = entropy();
    expect(m.V / V0).toBeLessThan(0.8); // the switch happened
    expect(s1).toBeGreaterThanOrEqual(s0 - 1e-9 * Math.abs(s0));
  });
});
