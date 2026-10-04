# Octane — design

A physically detailed, browser-based combustion-engine simulator rendered with three.js.
The first engine is the **Waukesha CFR F-1** single-cylinder, variable-compression-ratio,
four-stroke spark-ignition engine — the engine used to *define* octane numbers. It is the
simplest real engine that still exercises all the combustion physics we care about (spark,
kernel, turbulent flame, knock), and it gives us a hard end-to-end validation target: the
knock-limited compression ratio vs. PRF octane number.

The second engine is the **Ford Model T** (1924–25 high-head engine): an inline four with side valves
(an L-head chamber offset over the valves), fixed CR ≈ 4, a flywheel magneto feeding four trembler
(vibrator) coils through a roller timer, a carburettor with a butterfly throttle, and a free-running
crank that can drive the car through its planetary gearbox. It exercises everything the CFR does not:
several cylinders on one crank sharing the manifolds, a non-cylindrical chamber, a multi-spark
electromechanical ignition and a real road load. The engine is chosen at run time (engine picker,
`?engine=ford-model-t`); see §Multi-cylinder engines and the Ford Model T.

## Goals

1. **Physics first.** Every subsystem is modelled from first principles as far as a real-time
   0D/quasi-dimensional model allows, with each empirical closure cited to the literature.
2. **Verified.** Every module is tested against an oracle: Cantera (thermo, equilibrium,
   kinetics, flame speed), published data, analytic limits, and conservation laws.
3. **Observable.** The 3D view is driven by the same state the physics produces
   (piston, valves, flame radius, spark phase…) — nothing is animated "for looks".

## Architecture

```
src/
  physics/            pure TypeScript, no DOM, no three.js; runs in a Web Worker and in vitest
    core/             shared contracts: constants, species set, EngineSpec, OperatingPoint, snapshots
    thermo/           NASA-polynomial species thermo, mixtures, fuels, transport
    equilibrium/      Gibbs-minimisation chemical equilibrium (burned gas)
    chemistry/        NO kinetics (extended Zeldovich), ignition delay, Livengood–Wu, knock acoustics
    combustion/       laminar flame speed, 0D turbulence, flame geometry, entrainment burning model
    ignition/         coil circuit, gap breakdown, arc/glow, spark-kernel growth
    gas-exchange/     valve lift, valve/throttle flow, plenums
    heat-transfer/    Woschni & co., wall heat flux split
    mechanics/        slider-crank kinematics, crank dynamics, friction
    engines/          EngineSpec instances (cfr.ts, model-t.ts) + the engine registry (index.ts)
    cycle/            the engine (crank, manifolds, integrator) + N cylinders that wire everything together
  worker/             Web Worker hosting the simulator + message protocol
  render/             three.js scene: engine-model.ts (render-model contract + registry), engine/ (CFR
                      mechanism), engine-modelt/ (Model T mechanism), combustion/ (in-cylinder visuals)
  app/                page shell + one disposable EngineSession per selected engine
  ui/                 controls, charts, playback
tools/reference/      Python (Cantera) scripts that generate oracle data → test/fixtures/*.json
test/fixtures/        committed oracle data
```

Rules:
- `src/physics/**` never imports from `render/`, `ui/`, `worker/`, or three.js.
- `render/` and `ui/` never call physics solvers; they consume `EngineSnapshot`/`CycleSummary`
  (src/physics/core/snapshot.ts) and read `EngineSpec` geometry.
- `src/physics/core/*` is the shared contract. Change it deliberately, not ad hoc.

## Conventions

- **Units:** SI everywhere in `src/physics` (m, kg, s, K, Pa, J, **mol** — not kmol). UI converts.
- **Angles:** radians internally. Fields ending in `Deg` are crank degrees.
- **Crank angle:** 0 = firing TDC; cycle spans [-360°, 360°):
  intake −360→−180, compression −180→0, expansion 0→180, exhaust 180→360.
  Spark advance is given in degrees *before* firing TDC (spark at θ = −advance).
- **Frames:** see `src/physics/core/engine-spec.ts`. Cylinder frame: origin at the centre of the
  head fire-deck face, +y toward the head, gas occupies −h(θ) ≤ y ≤ 0. World frame: +y up,
  crank axis = world z.
- **Composition vectors:** `Float64Array(NS)` indexed by `SP` in `core/species.ts`. Names tell the
  basis: `N` = moles (mol), `X` = mole fractions, `Y` = mass fractions. Element vectors are
  `Float64Array(NE)` in `ELEMENTS` order.
- **Thermo reference:** absolute enthalpies (formation included, NASA convention, 298.15 K
  reference). Standard-state pressure is **1 bar** (`P_REF_THERMO = 1e5`), the true basis of the
  NASA Glenn data. Cantera's `nasa_gas.yaml` mislabels it as 1 atm, so every Cantera oracle must
  build its species with `tools/reference/thermo_common.py: nasa_species()` / `build_phase()`.
- **Hot paths allocate nothing:** reuse scratch `Float64Array`s; no closures/objects per step.
- **Every empirical constant carries a citation comment** (author, year, eq./table). If a value
  comes from memory and could not be checked against a fetched source, mark it
  `// UNVERIFIED:` so reviewers target it.
- **Tests** are colocated `*.test.ts` (vitest). Oracle data lives in `test/fixtures/<module>_*.json`
  generated by `tools/reference/<module>_*.py` (run with `.venv/bin/python`; Cantera 3.2 is
  installed there). Commit both the script and its JSON output.

## Physics model (quasi-dimensional, two-zone)

The cylinder is a control volume V(θ) from slider-crank kinematics.

**Gas exchange (valves open):** single well-mixed zone with tracked species moles N (NS) and
internal energy U. Valve and throttle flows are quasi-steady compressible isentropic nozzle
flows with lift-dependent discharge coefficients, bidirectional (backflow carries cylinder gas
into the intake). Intake and exhaust plenums are 0D volumes with their own N, U. Fresh charge =
humid air + fully vaporised fuel at φ (+ optional EGR) at the intake mixture temperature
(the CFR is carburetted with a heated mixture). Residual gas is whatever stays in the cylinder.
The intake plenum (runner + port) exchanges heat with the port walls while gas flows through the
port (Dittus–Boelter, calibrated multiplier); the exhaust plenum is adiabatic.

**Closed cycle (IVC → EVO):** two zones at a common pressure p, plus a crevice zone:
- *unburned*: frozen composition (fresh charge + residual + external EGR), compressed with entropy
  changing only by wall heat loss (and the crevice gas it gets back);
- *burned*: chemical equilibrium over the 12 product species at (T_b, p), except NO, which is
  rate-controlled (extended Zeldovich, Strang-split) and carried in the burned-zone energy and
  volume (the equilibrium NO is swapped for the kinetic amount; `burnedNOThermo`);
- *crevice* (piston top land + ring grooves, `geometry.quenchCreviceVolume`; Namazian & Heywood
  1982): gas at the wall temperature and the cylinder pressure, filled from the zone at the crevice
  mouth as p rises and emptied back as p falls; its cooling is wall heat. No blow-by yet.
Total energy is conserved exactly: dU = −p dV − Q̇_wall dt + Σ h ṁ.

**Spark ignition:** inductive coil circuit (primary charge during dwell, secondary ring-up at
switch-off), gap breakdown when the secondary voltage exceeds the density-dependent breakdown
voltage (Paschen-type law), then arc/glow discharge phases delivering energy to the gas with
phase-dependent efficiency. A spherical kernel grows from the gap with spark-energy input +
stretch-affected laminar burning (Herweg–Maly type) until it is large enough to become the
turbulent flame; kernels that are quenched → misfire.

**Turbulent flame:** entrainment/burn-up model (Blizard–Keck / Tabaczynski):
ṁ_e = ρ_u A_f (u′ + S_L), ṁ_b = ρ_u A_f S_L + (m_e − m_b)/τ_b, τ_b = λ/S_L,
with S_L from fuel-specific correlations (Metghalchi–Keck, Gülder, …, with residual dilution),
u′ from a 0D K–k turbulence model (Poulos–Heywood) including rapid-distortion compression,
and A_f from the exact geometry of a sphere centred at the spark gap intersected with the
disc-shaped chamber at the current piston position.

**Knock:** Livengood–Wu integral of an end-gas ignition-delay model — by default Douaud–Eyzat
with the relative φ, residual and low-ON (< 80) sensitivities of the LLNL detailed-chemistry table
(`douaud-eyzat-llnl`; plain Douaud–Eyzat and the LLNL tables, single- or two-stage, are options) —
on the unburned gas. At autoignition ALL unburned gas (the end gas ahead of the front and the
flame-brush pockets, `knockBrushAutoignition`) burns on τ_ab = τ_e + the ignition-time spread of a
ΔT-stratified charge from the integral's history; the heat is released as a sequential front
sweeping the autoigniting region from the far periphery toward the flame (16 shells), which drives
the cylinder's acoustic modes (Draper: f_mn = α_mn c / (π B)) synthesised onto the pressure trace.
MAPO is taken over the modes in the measurement band (4–18 kHz, `mapoBand`).

**Heat transfer:** Woschni correlation (alternatives: Hohenberg, Annand) split over
head/piston/liner and burned/unburned wetted areas from the flame geometry. Wall temperatures follow
the operating point through a lumped wall model: each surface sits at T_coolant + R_i·Q̄_i, the
cycle-mean heat flow to it times a resistance fitted at the Pal 2018 reference state
(`walls.thermalResistance`), so motored, fired, RON and MON walls differ.

**Mechanics:** exact slider-crank kinematics; gas + reciprocating/rotating inertia torques;
friction (FMEP correlation); either fixed speed (CFR synchronous motor) or free crank dynamics.

**Integration:** time is the independent variable (so free-running speed changes work). The
state vector is integrated with an explicit adaptive or fixed-step Runge–Kutta at crank-angle
resolution (~0.1–0.5°), refined during the spark event; stiff sub-models (NO kinetics, coil
circuit) are sub-cycled / operator-split. Algebraic closures (T, p from conserved quantities,
equilibrium) are solved by Newton iteration with warm starts.

## Module contracts

Implementers may add exports but must provide at least these. All inputs/outputs SI.

### thermo (`src/physics/thermo/`)
- `species-data.ts` — GENERATED by `tools/reference/gen_species_data.py` from Cantera's
  `nasa_gas.yaml` (NASA Glenn database, NASA-7 fits) for every species in `SPECIES`.
- `thermo.ts`
  - `MOLAR_MASS: Float64Array` (kg/mol), `P_REF_THERMO: number` (Pa)
  - `speciesCp(k, T)`, `speciesH(k, T)`, `speciesS0(k, T)`, `speciesG0(k, T)` — molar, J/mol(/K)
  - `evalAllNondim(T, cpR, hRT, s0R)` — all species cp/R, h/(RT), s°/R at T in one pass (hot path)
  - Outside the fit range: constant-cp extrapolation (h, s consistent), never polynomial blow-up.
- `mixture.ts`
  - conversions: `molesToX(N, X?)`, `xToY(X, Y?)`, `yToX(Y, X?)`, `elementMoles(N, b?)`
  - `mixMolarMass(X)`, `mixGasConstant(X)` (J/kg/K)
  - molar: `mixCpMolar`, `mixHMolar`, `mixUMolar`, `mixSMolar(X, T, p)` (ideal mixing)
  - mass: `mixCpMass`, `mixCvMass`, `mixHMass`, `mixUMass`, `mixSMass(X, T, p)`, `mixGamma(X, T)`
  - inverses (Newton, robust): `temperatureFromH(X, h, Tguess)`, `temperatureFromU(X, u, Tguess)`,
    `temperatureFromS(X, s, p, Tguess)` (mass-specific inputs)
- `fuels.ts`
  - `interface FuelBlend { label: string; X: Float64Array /* fuel-vapour mole fractions */; octaneNumber?: number }`
  - `fuelFromSelection(sel: FuelSelection): FuelBlend` — PRF octane number is **liquid-volume %**
    iso-octane; convert with liquid densities and molar masses (cite).
  - `dryAir(): Float64Array`, `humidAir(T, p, RH): Float64Array`
  - `stoichO2PerMolFuel(fuel)`, `stoichAirFuelRatio(fuel, airX)` (mass)
  - `freshCharge({fuel, phi, airX, egrFraction?, egrX?}): Float64Array` (mole fractions)
  - `completeCombustionProducts(reactantsX): Float64Array` (for initial residual guesses)
  - `lowerHeatingValue(fuel)` (J/kg fuel, from the thermo data)
- `transport.ts` — `mixViscosity(X, T)`, `mixThermalConductivity(X, T)`, `mixThermalDiffusivity(X, T, p)`.

### equilibrium (`src/physics/equilibrium/`)
- `class EquilibriumSolver` (holds scratch buffers + warm start):
  - `solveTP(b, T, p): EqResult` — Gibbs minimisation over species [0, N_EQ) for element moles b
  - `solveHP(b, H, p, Tguess)`, `solveUV(b, U, V, Tguess)`, `solveSP(b, S, p, Tguess)`
    (extensive H/U/S/V for the given element moles)
  - `EqResult { N: Float64Array(NS) mol, X, T, p, converged, iterations }`
- `equilibriumProperties(solver, b, T, p)` → `{ X, M, h, u, s, cp, cv, dlnV_dlnT, dlnV_dlnp }`
  (mass-specific h/u/s/cp; derivatives at equilibrium)
- `adiabaticFlameTemperature(Xreactants, T0, p0, mode: 'HP'|'UV')`
- Robust for T 300–4000 K, p 0.1–300 bar, φ 0.2–3 (lean/rich extremes, trace species ≥ 1e-30).

### chemistry (`src/physics/chemistry/`)
- `zeldovich.ts` — `zeldovichNORate(T, p, Xeq, xNO): number` → d[NO]/dt, mol/(m³ s)
  (extended Zeldovich with equilibrium O/OH/H/N2/O2 and the forward/reverse form of
  Heywood eq. 11.8; rate constants cited).
- `ignition-delay.ts` — `interface IgnitionDelayModel { id: string; tau(T, p, phi, fuel: FuelBlend, xResidual): number }`
  with `douaudEyzat` (PRF, octane-number based) at minimum; tabulated detailed chemistry if
  a mechanism is obtainable.
- `livengood-wu.ts` — `class LivengoodWuIntegrator { reset(); advance(dt, tau): number /* integral */ }`
- `knock.ts` — `cylinderModeFrequency(m, n, c, bore)`, `class KnockOscillator` (damped acoustic modes
  excited by end-gas autoignition energy), end-gas rapid-burn rate.

### combustion (`src/physics/combustion/`)
- `laminar-flame-speed.ts` — `laminarFlameSpeed(fuel: FuelBlend, phi, Tu, p, xDiluent): number` (m/s)
- `turbulence.ts` — 0D K–k model: `turbulenceDerivatives(state, inputs, out)`, `turbulenceIntensity(k, m)`,
  `integralLengthScale(...)`
- `flame-geometry.ts` — `class FlameGeometry(bore, sparkCenter)`: `evaluate(r, h, out)` →
  `{ volume, frontArea, wettedHead, wettedPiston, wettedLiner }` for a sphere of radius r centred
  at the gap intersected with the disc of height h; `radiusForVolume(V, h)`; option
  `allowSparkOutsideBore` (gap outside the bore planform, d > R).
- `chamber.ts` — `createChamber(spec)` → `CombustionChamber` (the cycle model's ONLY geometry entry):
  `evaluate(r, h, out)` (volume, front area, per-surface wetted areas and burned fractions, crevice
  burned fraction), `radiusForVolume`, `maxRadius`, `inscribedRadius`, `chamberVolume(h)` (= V − V_crevice,
  consistent with SliderCrank), `meanDepth(h)` (turbulence height), `surfaceAreas(h, out)`. `DiscChamber`
  wraps FlameGeometry bit-identically (CFR); `LHeadChamber` = bore column ∪ (pocket rounded rectangle ∖
  bore disc) prism, dV/dr = A_f exact on the fast path, tables shared by identical cylinders.
- `entrainment.ts` — `entrainmentRates({rhoU, frontArea, uPrime, SL, me, mb, lambda}, out)`,
  `taylorMicroscale(L, uPrime, nu)`.

### ignition (`src/physics/ignition/`)
- `breakdown.ts` — `breakdownVoltage(gap, p, T, X?)` (V)
- `coil.ts` — `class IgnitionCoil(spec)`: `step(dt, dwellOn, gap)` — primary/secondary circuit with gap model
- `kernel.ts` — `class SparkKernel`: kernel radius/temperature/mass growth from spark energy + chemistry,
  quench/misfire, hand-off radius
- `index.ts` — `class IgnitionSystem(spec)` composing the above: `step(dt, thetaDeg, command, gas, omega?,
  thetaEngineDeg?)`; `createIgnitionSystems(spec)` (one per cylinder).
- Trembler-magneto (`source.ts`, `vibrator.ts`, `trembler-coil.ts`, wired in `index.ts`): primary supply = 6 V
  battery or the AC flywheel magneto E = kω sin(N(θ − φ)) behind R_s + L_s; a 1-DOF vibrator armature
  (pull ∝ I₁²/(g₀ − x)²) opens the points (events in the coil sub-step), the condenser rings, the gap
  fires — a spark TRAIN for as long as the timer contact lasts. One ignition event = one timer contact;
  the kernel is not declared misfired between sparks and can be re-seeded. Each points opening gives ONE
  spark (breakdown) whose condenser ring re-ignites the recovering channel at its current zeros at
  V_r(t) = V_ri + (V_bd − V_ri)(1 − e^{−t/τ_rec}) (discharge.ts, τ_rec 30 µs UNVERIFIED; counted in
  `reignitionCount`, booked as arc energy); only a breakdown — including the timer-break spark within
  0.5 ms — re-seeds a quenched kernel. `tremblerTimerCommand(spec,
  lever)`: the spark lever sets the timer MAKE; the first spark is an output (coil firing time).

### gas-exchange (`src/physics/gas-exchange/`)
- `valve-lift.ts` — `valveLift(spec, thetaDeg)`, `valveLiftRate(spec, thetaDeg)` (m/deg);
  `createLiftProfile(valve, lash)` → `LiftProfile` (polydyne default | exact three-arc lobe on a flat
  follower (`cam-lift.ts`) | measured table), `resolveValveLash(valve, optionLash)` (ValveSpec.lash wins).
- `valve-flow.ts` — `valveFlowArea(spec, lift)` (Heywood curtain/seat regimes), `valveDischargeCoefficient(spec, lift, reverse)`;
  `createValveFlowModel(valve, kind, spec)` adds the side-valve pocket-roof stage and the pocket→bore
  transfer restriction (crown height as third argument of `effectiveArea`).
- `carburettor.ts` — `CarburettorFlowModel`: venturi and butterfly as two compressible orifices in series
  (Newton on the intermediate pressure, correct when the throttle chokes) when `manifolds.venturiDiameter`
  is set; the CFR keeps its venturi-as-throttle.
- `orifice.ts` — `orificeMassFlow(CdA, p0, T0, R0, gamma0, pDown)` (kg/s ≥ 0, choked/unchoked)
- `throttle.ts` — `throttleArea(diameter, opening)` (butterfly geometry + leakage)
- `plenum.ts` — 0D filling/emptying volume with composition.

### heat-transfer (`src/physics/heat-transfer/`)
- `woschni.ts` — `woschniCoefficient(inputs)` (W/m²/K), plus `hohenbergCoefficient`, `annandFlux`;
  `wallHeatLoss(...)` split over head/piston/liner.
- `wall-heat.ts` — `wallHeatLossSurfaces` / `wallHeatLossTwoZoneSurfaces`: six surfaces (head, piston,
  liner, intake valve, exhaust valve, block deck) each with its own burned fraction; summed in the legacy
  order so the CFR's five-surface numbers are bit-identical.

### mechanics (`src/physics/mechanics/`) and engine spec (`src/physics/engines/cfr.ts`)
- `kinematics.ts` — `class SliderCrank(geometry, compressionRatio)`: `volume(θ)`, `dVdTheta(θ)`,
  `clearanceHeight(θ)`, `pistonDisplacement(θ)`, `rodAngle(θ)`, `pistonVelocity(θ, ω)`,
  `pistonAcceleration(θ, ω, α)`, `displacedVolume`, `clearanceVolume`, `linerArea(θ)`, …
- `dynamics.ts` — gas torque, reciprocating/rotating inertia torque, crank angular acceleration;
  `MultiCylinderCrankTrain.fromSpec(spec, kin)`: rigid inline crank, θ_i = θ − (offset_i mod 360°),
  per-cylinder pressures, load inertia.
- `friction.ts` — FMEP correlation; `FrictionTorqueModel(kin, ε, N)` + `torqueCylinders` (whole-engine FMEP
  over N·V_d — constructing it without N gives 1/N of the friction); PNH valvetrain type 'L-head'.
- `load.ts` — `LoadModel`: constant / brake (∝ n^k) / vehicle road load through the gearbox (reflected
  inertia, drive vs overrun efficiency) / neutral.
- `CFR_F1`, `MODEL_T: EngineSpec` — every number sourced or marked UNVERIFIED.

### cycle (`src/physics/cycle/`) — built after the modules
- `class EngineSimulator(spec, op)`: `advance(dt)`, `snapshot()`, emits `CycleSummary`.
- `cycle-model.ts` = the ENGINE (crank angle/speed, shared intake/exhaust plenums, carburettor/venturi and
  outlet, RK4 stepping with event landing, warm-up, engine wrap and `EngineCycleSummary`); `cylinder.ts` =
  one cylinder at its local angle θ_i = θ − layout.firingOffsetDeg[i] (state block, two-zone closure,
  ignition, flame, knock, NO, events, local wraps and its `CycleSummary`). Cylinder 0 keeps the
  single-cylinder state layout, so a one-cylinder spec runs the former model bit for bit.
- Options resolve by `EngineSpec.id` (options.ts `ENGINE_CYCLE_OPTION_DEFAULTS`): neutral physics defaults ←
  the engine's defaults (calibration set, friction, knock sensor/band, venturi C_D) ← the caller's options.

## Validation ladder

1. Species thermo = Cantera to 1e-9 rel.; mixture inverse functions round-trip.
2. Equilibrium = Cantera `equilibrate('TP'|'HP'|'UV')` (majors 1e-6 abs, T_ad ±0.5 K).
3. Laminar flame speed vs Cantera `FreeFlame` (CH4/GRI-3.0, H2) and published iso-octane data.
4. NO formation rate vs Cantera (thermal-NO sub-mechanism of GRI-3.0).
5. Energy conservation of the closed cycle (≤1e-6 relative drift without heat loss);
   motored cycle reversibility; fuel–air-cycle efficiency vs a Cantera fuel–air-cycle oracle.
6. CFR data: motored & fired pressure traces, IMEP, and knock-limited compression ratio vs PRF
   octane number (ASTM guide table).

## Integration notes (from the module reviews)

Facts the cycle model must respect. Each module's `index.ts` JSDoc has the details.

- **thermo:** use the molar functions on mole vectors N to get extensive totals
  (`temperatureFromUMolar(N, U)`). `xDiluent`/residual for flame speed = **mass fraction of
  complete-combustion products**; ignition delay takes the residual **mole** fraction
  (`residualMoleFraction` converter in chemistry). `FuelBlend` objects are immutable.
- **equilibrium:** `EquilibriumSolver` is warm-started and owns its result object (copy what you
  keep). Call `properties()` right after a solve instead of re-solving. For the burned zone,
  (∂h/∂p)_T = v(1 − ∂lnV/∂lnT). No condensed carbon: φ ≳ 3.1 is infeasible.
- **chemistry:** detailed-chemistry knock model = `prfDetailedChemistry` (LLNL 2011 gasoline-surrogate PRF
  chemistry, tabulated 550–1100 K, 3.5–80 bar) integrated with `LivengoodWuIntegrator` from IVC;
  `douaudEyzat` is the classic alternative; since validation round 2 the cycle model's default is
  `douaudEyzatLLNL` (Douaud–Eyzat with the LLNL table's relative φ/residual/low-ON sensitivities,
  calibration.ts `CFR_KNOCK_DELAY_MODEL`). After autoignition burn the end gas with the exact
  `endGasBurnedMass` and sub-step it and `KnockOscillator` at ≲ 2 µs. NO via
  `ZeldovichKinetics.advanceRateControlled` (exact step) on the burned-zone equilibrium state.
  Livengood–Wu is poor for two-stage (low-ON) fuels: up to ~16° early for PRF 0.
- **combustion:** S_L comes from tables of Cantera 1-D flames (GRI-3.0 CH4; PRF_HT for
  iso-octane/n-heptane; energy-fraction PRF blending; φ-dependent thermal dilution factor).
  The kernel needs `marksteinLengths(...).unburned` (NOT the burned-gas length) and
  l_F = ν_u/S_L. `FlameGeometry(bore, gap, {maxHeight: stroke + h_TDC(CR_min)})`; `dV/dr` equals
  the front area exactly on the fast path. Turbulence: K–k model with integral scale
  L = 0.25·min(h, B) (calibrated so u′(TDC) ≈ 0.4–0.6·S̄p); production is seeded so k can start
  from 0. Entrainment: Keck 1982 eqs 4.1A/B, u_T = u′, burn-up length = Taylor microscale —
  Keck's own empirical u_T at CFR conditions is ~2.6× our u′, so the burn rate is the prime
  calibration knob (calibrated: C_T = 5.3, C_λ = 3.0 — `cycle/calibration.ts`). S_L = 0 with u′ > 0 means extinction (entrained mass never burns).
- **ignition:** `IgnitionSystem.step(dt, thetaDegAtEndOfStep, {dwellStartDeg, sparkDeg}, gas)`,
  sub-steps the coil internally (~18 ms CPU per cycle). Hand-off radius max(1 mm, C·l_I) with the
  integral scale l_I = C_ε·L (C_ε 0.5) and C = 2 (calibration.ts); strain uses the dissipation length L.
  The real CFR ignition is capacitive discharge (`CFR_CDI_IGNITION`, values unverified); the
  inductive coil in `CFR_F1.ignition` is a stand-in.
- **gas exchange:** use `new ValveLiftProfile({...spec, lash: CFR_VALVE_LASH})` (lash-corrected; the
  bare spec gives the no-lash lift) or the measured `cfrValveLiftChoi2018`. Valve events at
  running clearance (Choi 2018): IVO −344, IVC −152, EVO 141, EVC −347 (3° with both valves
  shut). `orificeFlow` is smoothed near Δp = 0; the resulting equalisation time constant is
  ~80 µs for the CFR intake — keep explicit steps below ~2.5× that or treat flows implicitly.
  The ideal shrouded-valve jet gives swirl ratio ~23 at IVC (too strong): scale the swirl torque.
  The CFR has no throttle — `manifolds.throttleDiameter` is the 9/16 in carburettor venturi.
- **heat transfer:** `woschniCoefficient` defaults to the original 1967 form; it needs the IVC
  reference state and motored pressure in the combustion/expansion phases. Choi et al. 2018
  needed chamber-convection multipliers of 1–3 (and ×4 in the intake port) to match the CFR.
- **mechanics:** `SliderCrank.evaluate(θ, out)` is the hot-path call. Crevice volume is part of the
  clearance volume; the cycle model stores the gas of `geometry.quenchCreviceVolume` in a crevice
  zone. `CrankTrainDynamics` for free-speed mode; `FrictionTorqueModel` + `pnhFmep`
  with `CFR_FRICTION` (FMEP ≈ 1.2 bar at 600 rpm, model-only). Wall temperatures: lumped
  T_i = T_coolant + R_i·Q̄_i with R_i fitted to Pal 2018 Table 3 at the Choi/Pal standard-knock state.
- **Calibration policy:** only genuinely uncertain closures may be tuned (burn-rate multiplier /
  Taylor-scale constant, turbulence length-scale factor, Woschni multiplier, intake-port heat
  transfer, discharge coefficients, knock-intensity threshold / end-gas stratification), each within
  a literature-supported range, with ONE parameter set per ENGINE for all its operating points (CFR: RON,
  MON, all octane numbers), and every choice documented next to the parameter. The sets live in
  `src/physics/cycle/calibration.ts` (`CFR_CALIBRATION`, `MODEL_T_CALIBRATION`: value, range, source,
  evidence) and are selected by `EngineSpec.id`. Universal closures (Markstein ratio, kernel hand-off
  multiple, knock-delay model) are shared; an engine fitted to little data says so in its evidence.
- **cycle (validation round 1, fixer pass):**
  - Standard knock = MAPO ≈ 0.67 bar at the D-1 pickup (Rockstroh 2018; Hoth 2021 Table 6); the
    model's MAPO/peak pressure/max dp/dθ are sampling-independent (see cycle/index.ts definitions).
  - Combustion sub-models (S_L, Markstein/Lewis, τ, LHV) use the TRAPPED mixture: a fuel or φ change
    reaches the cylinder over a few cycles through the intake plenum (0.25 L since round 2).
  - The exhaust system is adiabatic (plenum ≈ 940 K vs 641 K measured thermocouple).
  - Known model-form limits of the calibrated model: the best single burn-rate multiplier differs
    between spark timings (u′ decays after TDC); Douaud–Eyzat knocks too early for PRF ≤ 70 (RON)
    and 0.35–0.6 CR too late in MON; the knock-limited CR responds ≈ 4× too weakly to spark advance;
    ASTM motored compression pressures read +4…+6 % (counter→CR relation, ±0.4 CR between ANL
    campaigns) while the Choi 2018 compression at the stated CR is −1.7 %.
  - Hot paths allocate no objects/arrays; V8 boxes doubles at non-inlined call boundaries
    (≈ 17 MB per 600 rpm cycle, GC ≈ 1 % of CPU).
- **cycle (validation round 2, fixer pass):**
  - Crevice zone (Namazian & Heywood 1982): top-land/ring gas at the mean piston/liner temperature
    and cylinder pressure, ṁ_cr = c·ṗ from the closure's rate linearity plus a 0.5 ms constraint
    keeper; unburned gas flows in during compression, burned gas later, everything returns at EVO.
    Active for the entrainment and motored models (the instantaneous test burn jumps p).
  - Walls: lumped (options.wallTemperatureModel 'lumped', default): motored walls sit a few K above
    the coolant, fired walls reproduce Pal 2018 Table 3 at the reference state; 'fixed' keeps the
    spec temperatures.
  - Rate-controlled NO is Strang-split and its energy is carried in the burned-zone closure
    (burnedNOThermo 'auto'); snapshot HRR is the instantaneous rate.
  - The venturi meters humid AIR at the intake-air temperature (the mixture heater is downstream);
    external EGR is a tracked intake scalar; the residual reported excludes EGR.
  - Knock: all unburned gas (end gas + flame-brush pockets) autoignites at the Livengood–Wu onset;
    the acoustic source is released sequentially from the periphery inward (16 shells); MAPO is the
    ideal 4–18 kHz band. The Herweg–Maly kernel survives non-propagating phases for
    min(l_I/u′, r²/α) before it is declared extinct; the Markstein length is scaled by 0.5 (Bradley 1998).
  - The ANL 2020s CR scale (SON 2023, Hoth 2025, Kalvakala) is converted with the clearance offset
    `CFR_ANL2020S_CLEARANCE_OFFSET` (7.38 cm³, `cfrCompressionRatioFromAnl2020s`).
  - Known model-form limits after round 2: the knock-limited CR scale sits 0.14–0.45 CR below the
    guide table at ON 80–100 (end-gas stratification at its bound; the model's MAPO differs 1.2× where
    the two ANL standard-knock states differ 2.2×), PRF ≤ 80 in RON knock 0.45–0.94 CR too early, the
    λ trend of MAPO is too weak; MBT ≈ 5° BTDC vs 9–10° measured; blowdown too fast (p(160)/p(141)
    −6.5 %); MON ASTM compression pressure +4.9 % (MON/RON +2.9 %: no mixture-heater flow loss).

## Multi-cylinder engines and the Ford Model T

**Contracts.** `EngineSpec.layout` gives the firing order, the firing-TDC offset of every cylinder
(Model T 1-2-4-3: [0, 180, 540, 360]), the cylinder positions along the crank and the render-only mirror
flags; absent = one cylinder at the origin. `geometry.chamber: 'l-head'` + `geometry.lHead` describe the
side-valve chamber (bore column from the crown up to the cavity roof y = 0 ∪ a valve-pocket prism from
the deck to the pocket roof; the Model T crown rises 5/16 in above the deck at TDC). `ValveSpec.cam`,
`lash`, `seatY`, `liftDirection` describe a side valve driven by a real cam. `ignition.type
'trembler-magneto'` is the Ford system. `EngineSpec.vehicle` + `OperatingPoint.load` ('constant' | 'brake'
| 'vehicle') give the free-speed load; `OperatingPoint.ignitionSource` is the dash MAG/BAT switch and
`sparkAdvanceDeg` is the spark LEVER (timer make) for a trembler engine.

**Snapshots.** A multi-cylinder snapshot carries `cylinders[i]` (every cylinder at its local angle, own
cycle counter and gas torque); the top-level per-cylinder fields are cylinder 1, and the engine-level
fields are `gasTorque` (sum), `netTorque`, `frictionTorque` (> 0 opposing), `loadTorque`, `vehicleSpeed`,
`firingCylinder`, `magnetoEmf`. Each cylinder emits its own `CycleSummary` (`summary.cylinder`);
cylinder 1's carries `summary.engine` (brake torque, power, BMEP, FMEP, η_v over N·V_d, bsfc). Brake
torque is the ENGINE output, indicated − friction − (J_rot + ΣJ_m)·ω̇: the car's inertia is downstream of
the clutch and appears separately (`loadInertiaTorque`; a gear change is an inelastic clutch engagement
conserving angular momentum through the gear train, its loss in `clutchLoss`). The CFR's
snapshot stream is unchanged (no `cylinders`).

**Engine registry and app.** `src/physics/engines/index.ts` maps ids to {spec, default operating point,
presets, UI profile}. The app shell (`src/app/app.ts`) keeps the WebGL stage alive and swaps a disposable
`EngineSession` (render model from `createEngineModel(spec, op)`, one CombustionVisuals per cylinder frame
with tracers and the chamber light on the featured cylinder, SimClient, Conductor, UI). The worker needs
no engine knowledge: the spec it receives selects the physics defaults. The Model T cutaway's quarter
section follows the focus cylinder (`EngineRenderModel.setSectionCylinder`). `?engine=<id>` selects the
engine; a link with operating-point knobs but no engine opens the default engine (the knobs were written
for it), otherwise the last engine used is remembered.

**Validation.** `test/fixtures/modelt_period_data.json` (`tools/reference/modelt_data_period.py`, README
alongside): Ford's WOT brake torque table, Upton's MBT relation, compression pressure, trembler-coil and
magneto behaviour, timer/lever range and vehicle performance — global targets only (no Model T pressure
trace exists in open sources). Known model-form limits: Euclidean flame sphere in the non-convex L-head
(the pocket burns slightly early), knock model weakest at the period fuel's low octane, PNH friction
extrapolated to a 1920s babbitt-bearing engine.
