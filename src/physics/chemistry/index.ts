/**
 * Rate-controlled chemistry (public API of src/physics/chemistry): thermal NO (extended
 * Zeldovich), end-gas ignition delay (Douaud–Eyzat; LLNL detailed-chemistry tables — default
 * Mehl et al. 2011 PRF chemistry, alternative LLNL PRF v2),
 * Livengood–Wu knock integral, knock acoustics and end-gas burn-up.
 * See DESIGN.md §Module contracts → chemistry. SI units (mol, not kmol).
 *
 * Intended use by the cycle integrator (per time step dt, closed cycle):
 *  - NO: with the burned-zone equilibrium mole fractions Xeq at (T_b, p) from the equilibrium
 *    solver, advance the burned-zone NO mole fraction with
 *    `ZeldovichKinetics.advanceRateControlled(T_b, p, Xeq, xNO, dt)` (exact, stable at any dt),
 *    or use `zeldovichNORate` (mol/(m³ s), × V_b for dN_NO/dt) inside an explicit RK stage.
 *    Newly burned gas enters with the unburned NO (≈ 0).
 *  - Knock: `LivengoodWuIntegrator.reset()` at IVC; each step `advance(dt, τ)` with
 *    τ = model.tau(T_u, p, φ, fuel, x_res) at the END of the step (τ from `prfDetailedChemistry`
 *    or `douaudEyzat`). When `autoignited` first becomes true (onset time `ignitionTime`):
 *      τ_ab = autoignitionBurnTime(τ, ∂lnτ/∂T via ignitionDelayTemperatureSensitivity),
 *      burn the remaining unburned mass (in addition to the flame) with the exact, operator-split
 *      Δm = endGasBurnedMass(m_u, τ_ab, dt) per step — τ_ab ≈ 1–100 µs is usually shorter than
 *      the step, so the explicit rate endGasBurnRate(m_u, τ_ab) is unstable inside RK stages —
 *      `KnockOscillator.setEndGasRegion(flame centre x, z, flame radius)`, and from then on call
 *      `KnockOscillator.step(dt, c̄, γ, V, Q̇_endgas)` with Q̇_endgas = Δm·q_chem/dt (step-mean
 *      chemical heat-release rate of the end-gas burn, W) and c̄ = twoZoneSoundSpeed(...). Add
 *      `sensorPressure()` to the reported cylinder pressure (EngineSnapshot.knock.oscillation);
 *      MAPO = max |sensorPressure|. The step-mean Q̇ is a zero-order hold: a release faster than
 *      the step is smeared over dt, which filters mode j by |sinc(ω_j dt/2)| (0.94 for the (1,0)
 *      mode, 0.74 for (0,1) at dt = 28 µs, c = 1000 m/s). Sub-step the end-gas burn and
 *      KnockOscillator.step at dt ≲ 2 µs during the burn-up to keep this below 1 %.
 *    The ignition-delay model must cover the states the end gas visits from IVC on:
 *    prfDetailedChemistry spans 550–1100 K and 3.5–80 bar (compression stroke included);
 *    prfLLNLv2 only 650–1100 K / 10–80 bar and extrapolates below that (the same extrapolation
 *    of the 2011 table was 2–33× too reactive there, so v2 is not fit for full-cycle LW).
 *    A NaN τ (non-PRF fuel with a PRF model) makes the integral NaN rather than firing knock.
 */
export * from './zeldovich';
export * from './ignition-delay';
export * from './ignition-delay-llnl';
export * from './livengood-wu';
export * from './knock';
