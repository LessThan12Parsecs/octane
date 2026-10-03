/**
 * Detailed-chemistry PRF ignition delays, tabulated from Cantera 3.2 adiabatic constant-volume
 * 0-D reactor simulations from the unburned state (T, p) — see
 * tools/reference/chemistry_ignition_table.py for the exact mixture and ignition definitions
 * (τ = time of max dT/dt; τ1 = first-stage dT/dt maximum where distinct):
 *
 *  - prfLLNLGasoline2011 (DEFAULT): LLNL gasoline-surrogate detailed mechanism, PRF
 *    sub-mechanisms (Mehl, Pitz, Westbrook & Curran, Proc. Combust. Inst. 33 (2011) 193;
 *    LLNL-MI-536371; 1382 species after removing 7 pentene species without thermo data).
 *    Against the Fieweger et al. (1997) n-heptane shock-tube data (40 bar, φ = 1, 752–1143 K)
 *    it is 1.17–1.73× slow (geometric mean 1.33); it shows the iso-octane NTC at 40 bar.
 *  - prfLLNLv2 (comparison only): LLNL PRF v2 (Curran, Gaffuri, Pitz & Westbrook, Combust.
 *    Flame 129 (2002) 253; 1034 species): 1.07–2.44× slow vs Fieweger (geometric mean 1.57),
 *    ~2–3× slower than the 2011 model for PRF 80–100 below 750 K, and a spurious low-T
 *    sensitivity to trace CO (its CO + O2 → CO2 + O fit, 1.07e-15 T^7.13 exp(−13.3 kcal/RT),
 *    is ~700× the Baulch-type value at 650 K and acts as a radical initiator), which corrupts
 *    its rich-residual (x_res > 0, φ > 1) entries below ~750 K.
 *
 * Domains (φ 0.5–1.5, PRF 0–100, residual 0–15 % mole for both; outside see
 * TabulatedIgnitionDelay for the extrapolation rules):
 *  - prfLLNLGasoline2011: T 550–1100 K (1000/T nodes: 17 uniform 1100→650 K + 600, 550 K),
 *    p 3.54–80 bar (7 nodes uniform in ln p). It covers the compression stroke, where the
 *    Livengood–Wu integral starts. Accuracy of the interpolated ln τ vs direct Cantera:
 *    150 random points in 650–1100 K / 10–80 bar |Δln τ| median 0.026, p90 0.064, max 0.10;
 *    84 points in 550–700 K / 3–10 bar median 0.049, p90 0.116, max 0.255; along the 7
 *    engine-compression trajectories of test/fixtures/chemistry_lw_engine.json ≤ 0.041 for
 *    T ≥ 640 K and ≤ 0.16 at 578–624 K (coarse 550/600/650 K nodes; τ 0.2–1 s there).
 *    (The first version stopped at 650 K / 10 bar; its extrapolation was 2–33× too reactive
 *    below that — median |Δln τ| 0.32, max 3.5 — and put 0.05–0.4 of spurious integral on the
 *    end gas before it reached 650 K.) 14 lean, diluted entries at T ≤ 650 K, p ≤ 5.95 bar
 *    did not ignite within 100 s and hold τ = 100 s.
 *  - prfLLNLv2: T 650–1100 K, p 10–80 bar only (not extended; comparison use).
 * The JSON tables (0.17 / 0.11 MB) are imported statically so the worker bundle carries them.
 */
import gasoline2011 from './data/prf-ignition-llnl-gasoline-2011.json';
import prfV2 from './data/prf-ignition-llnl-prf-v2.json';
import { douaudEyzatTau, type IgnitionDelayModel, type IgnitionDelayTable, prfOctaneNumber, TabulatedIgnitionDelay } from './ignition-delay';
import type { FuelBlend } from '../thermo/fuels';

/** LLNL gasoline-surrogate (Mehl et al. 2011) PRF constant-volume ignition delays. */
export const prfLLNLGasoline2011 = new TabulatedIgnitionDelay(
  'llnl-gasoline-2011-cv',
  gasoline2011 as IgnitionDelayTable,
);

/** LLNL PRF v2 (Curran et al. 2002) constant-volume ignition delays — comparison only. */
export const prfLLNLv2 = new TabulatedIgnitionDelay('llnl-prf-v2-cv', prfV2 as IgnitionDelayTable);

/** Default detailed-chemistry PRF ignition-delay model. */
export const prfDetailedChemistry = prfLLNLGasoline2011;

// ---------------------------------------------------------------------------------------
// Douaud–Eyzat with detailed-chemistry corrections (fixer round 2)
// ---------------------------------------------------------------------------------------

/** Octane number down to which Douaud–Eyzat is used as fitted (see douaudEyzatLLNL). */
export const DE_LLNL_ON_MIN = 80;
/** Reference fresh-charge φ and residual mole fraction of the Douaud–Eyzat fit (see douaudEyzatLLNL). */
export const DE_LLNL_PHI_REF = 1.1;
export const DE_LLNL_XRES_REF = 0.06;
/** Temperature / pressure range in which the LLNL correction ratios are evaluated (clamped), K / Pa. */
const DE_LLNL_T_RANGE: readonly [number, number] = [600, 1100];
const DE_LLNL_P_RANGE: readonly [number, number] = [3.6e5, 80e5];

/**
 * Douaud–Eyzat (1978) end-gas delay with the RELATIVE sensitivities it lacks taken from the LLNL
 * detailed-chemistry table (prfLLNLGasoline2011):
 *   τ = τ_DE(T, p, ON*) · τ_L(T′, p′, φ, ON, x_res) / τ_L(T′, p′, φ_ref, ON*, x_ref),  ON* = max(ON, ON_min)
 *  - φ and residual: Douaud–Eyzat is a fit to CFR knock at the maximum-knock mixture of the rating
 *    procedure (φ_ref = 1.1, CFR_PEAK_KNOCK_PHI; x_ref = 0.06, the CFR's ≈ 6 % residual) and has no φ
 *    or dilution term (validation round 2: the model's knock rose monotonically toward λ = 1 where the
 *    CFR's falls 2.4×, Hoth & Kolodziej 2021a Fig. 6);
 *  - octane numbers below ON_min = 80: Douaud–Eyzat's (ON/100)^3.402 → 0 made PRF 0–10 "autoignite"
 *    at IVC (370 K, 1 bar) and PRF 60–70 knock too early; below ON_min the ON dependence is that of
 *    the detailed chemistry, anchored to Douaud–Eyzat at ON_min (80: lowest ON at which round 1
 *    verified the Douaud–Eyzat knock-limited CR against the ASTM guide table, −0.02 CR);
 *  - T′, p′ = T, p clamped to 600–1100 K / 3.6–80 bar (the table's accurate range; below 600 K the
 *    correction is frozen at its 600 K value — the table is coarse there and extrapolated below 550 K).
 * At (φ_ref, x_ref, ON ≥ ON_min) it IS Douaud–Eyzat. Model form (UNVERIFIED as a combination): the
 * absolute delay and its T, p and ON (≥ 80) dependence stay those of the CFR fit.
 */
export const douaudEyzatLLNL: IgnitionDelayModel = Object.freeze({
  id: 'douaud-eyzat-llnl-2011',
  tau(T: number, p: number, phi: number, fuel: FuelBlend, xResidual: number): number {
    const on = prfOctaneNumber(fuel);
    if (Number.isNaN(on)) throw new RangeError(`Douaud–Eyzat needs a PRF octane number (fuel '${fuel.label}')`);
    const onStar = on > DE_LLNL_ON_MIN ? on : DE_LLNL_ON_MIN;
    const Tc = T < DE_LLNL_T_RANGE[0] ? DE_LLNL_T_RANGE[0] : T > DE_LLNL_T_RANGE[1] ? DE_LLNL_T_RANGE[1] : T;
    const pc = p < DE_LLNL_P_RANGE[0] ? DE_LLNL_P_RANGE[0] : p > DE_LLNL_P_RANGE[1] ? DE_LLNL_P_RANGE[1] : p;
    const L = prfLLNLGasoline2011;
    const lnR = L.lnTauAt(Tc, pc, phi, on, xResidual) - L.lnTauAt(Tc, pc, DE_LLNL_PHI_REF, onStar, DE_LLNL_XRES_REF);
    return douaudEyzatTau(T, p, onStar) * Math.exp(lnR);
  },
});
