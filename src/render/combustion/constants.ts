/**
 * Constants for the in-cylinder visuals.
 *
 * Three kinds of numbers live here, and each is labelled:
 *
 *  1. PHYSICAL approximations (emission efficiencies, discharge channel sizes,
 *     turbulence rules of thumb). They have a citation or are marked
 *     `UNVERIFIED:` so a reviewer can target them. These set RELATIVE
 *     brightness between emitters (spark ≫ flame front ≫ burned gas).
 *  2. The single EXPOSURE constant `VISUAL_GAIN`. It maps physical luminance
 *     (cd/m²) to display units. It is the ONLY place light intensities are
 *     scaled for display; nothing else "boosts" a light source.
 *  3. VISUALISATION-ONLY settings that are not light at all (false-colour
 *     opacity, knock-pattern overlay scale, tracer-particle parcel size,
 *     display persistence). They are named `VIS_*` / `OVERLAY_*` / `TRACER_*`.
 */

// ---------------------------------------------------------------------------
// 2. Exposure (the only brightness gain)
// ---------------------------------------------------------------------------

/**
 * VISUAL GAIN — camera exposure: display units (linear sRGB, 1 = diffuse
 * white) per cd/m² of physical luminance.
 *
 * Chosen so a typical CFR turbulent flame front seen face-on
 * (q'' ≈ 6e7 W/m² → ≈ 8.6e4 cd/m² with the chemiluminescence model below)
 * has its (dominant) blue channel at ≈ 0.9 of display white. With this
 * exposure the spark channel (~1e8 cd/m²) saturates to white and the burned
 * gas (~1e3–1e4 cd/m²) is faint — what a camera exposed for the flame records.
 */
export const VISUAL_GAIN = 3e-6;

/** Luminous efficacy of radiation at 555 nm, lm/W (CIE definition). */
export const K_M = 683;

// ---------------------------------------------------------------------------
// 1. Physical approximations
// ---------------------------------------------------------------------------

/**
 * Fraction of the flame-front heat release radiated as chemiluminescence
 * (OH*, CH*, C2*, CO2*; 250–700 nm).
 * UNVERIFIED: order of magnitude 1e-4 for hydrocarbon premixed flames
 * (cf. Gaydon, "The Spectroscopy of Flames", 1974; Najm et al. 1998 report
 * CH* : OH* signals as tiny fractions of heat release). Sets absolute
 * flame brightness together with VISUAL_GAIN.
 */
export const CHEMILUMINESCENCE_EFFICIENCY = 2e-4;

/**
 * Same for end-gas autoignition (low/high-T chemistry: HCHO*, HCO*, CO2* continuum).
 * UNVERIFIED: assumed equal to the flame value.
 */
export const AUTOIGNITION_CHEMILUMINESCENCE_EFFICIENCY = 2e-4;

/**
 * Effective visible absorption (= emission, Kirchhoff) coefficient of hot,
 * soot-free burned gas, 1/m. Burned gas radiates mainly in the IR (H2O, CO2
 * bands); visible emission comes from band wings and the CO+O continuum.
 * UNVERIFIED: chosen so a few cm of 2400 K gas is ~1e-3 emissive in the visible.
 */
export const BURNED_GAS_VISIBLE_KAPPA = 0.05;

/**
 * Fraction of the gap's electrical power radiated in 300–800 nm by the arc /
 * glow column. UNVERIFIED: a few % (Maly & Vogel 1979 energy balance: most
 * energy goes to heat/conduction; radiation is a minor loss channel).
 */
export const SPARK_RADIATIVE_EFFICIENCY = 0.03;

/**
 * Fraction of the breakdown energy ½·C·V_bd² radiated in the visible.
 * UNVERIFIED: breakdown plasma (T ~ 3–6e4 K, Maly 1984) radiates strongly,
 * much of it in the UV; ~5 % assumed visible.
 */
export const BREAKDOWN_RADIATIVE_EFFICIENCY = 0.05;

/**
 * Arc column current density, A/m² (channel radius a = sqrt(I / (π j))).
 * UNVERIFIED: ~1e7 A/m² gives a ≈ 50 µm at 0.1 A, consistent with the
 * 40–100 µm arc channels imaged by Maly (1984).
 */
export const ARC_CURRENT_DENSITY = 1.3e7;

/**
 * Glow-discharge column radius at GLOW_REFERENCE_CURRENT, m; the column
 * cross-section scales with current (≈ constant current density, as in a
 * normal glow), so a ∝ √I. UNVERIFIED: 0.1–0.3 mm diffuse column at ~50 mA.
 */
export const GLOW_COLUMN_RADIUS = 1.5e-4;
export const GLOW_REFERENCE_CURRENT = 0.05;

/** Arc column temperature for its continuum colour, K. UNVERIFIED: 5000–7000 K (Maly 1984). */
export const ARC_TEMPERATURE = 6000;
/** Breakdown plasma temperature for its continuum colour, K. UNVERIFIED: ~3e4 K initially (Maly 1984). */
export const BREAKDOWN_TEMPERATURE = 25000;

/**
 * Integral length scale as a fraction of the instantaneous clearance height
 * for a flat disc chamber. UNVERIFIED: L_I ≈ 0.2 h near TDC is a common
 * quasi-dimensional closure (cf. Heywood 1988 §8.2; Tabaczynski et al. 1977).
 */
export const INTEGRAL_SCALE_FRACTION = 0.2;

/**
 * Fallback turbulence intensity when the snapshot has none, as a fraction of
 * mean piston speed. Heywood (1988) §8.2.2: u' ≈ 0.5 S̄p at TDC for open
 * chambers without swirl.
 */
export const UPRIME_OVER_MEAN_PISTON_SPEED = 0.5;

/** Kinematic viscosity of air at 300 K, 1 bar, m²/s (ν ∝ T^1.7 / p used for scaling). */
export const NU_AIR_REF = 1.6e-5;

/** First zero of J1' — (m, n) = (1, 0) circumferential acoustic mode, Draper (1938): f = α c / (π B). */
export const ALPHA_10 = 1.8412;

/** Ratio of specific heats used for the bulk sound speed during knock (hot mixed gas). UNVERIFIED: 1.28–1.32. */
export const KNOCK_GAMMA = 1.3;

/**
 * Intake/exhaust curtain discharge coefficient used only to estimate tracer
 * jet speed. UNVERIFIED: typical 0.6 at mid lift (Heywood 1988 Fig. 6-19).
 */
export const CURTAIN_DISCHARGE_COEFF = 0.6;

/**
 * Fraction of the curtain jets' ideal angular-momentum flux (ṁ·v_h·arm) that
 * ends up as solid-body swirl; the rest is lost in jet–wall/jet–jet
 * interaction during induction. UNVERIFIED: chosen so a 180° shrouded valve
 * gives a swirl ratio of order 2–4 at IVC, the range quoted for
 * shrouded-valve / directed-port engines (Heywood 1988 §8.3).
 */
export const SWIRL_TRANSFER_EFFICIENCY = 0.15;

/** Upper bound on the swirl ratio ω_s/ω_crank (guards the crude swirl model). */
export const MAX_SWIRL_RATIO = 6;

/** Swirl decay time in crank revolutions (wall shear). UNVERIFIED: swirl loses ~30 %/rev. */
export const SWIRL_DECAY_REVS = 2.5;

// ---------------------------------------------------------------------------
// 3. Visualisation-only settings (not light)
// ---------------------------------------------------------------------------

/** Default false-colour temperature range, K. */
export const VIS_TEMPERATURE_RANGE: readonly [number, number] = [250, 3000];

/** False-colour gas opacity per metre in 'temperature' mode (8 cm → ~90 %). */
export const VIS_TEMPERATURE_KAPPA = 30;

/**
 * Faint haze of the unburned charge in 'physical' mode (a visual cue — a
 * vapour/air mixture is transparent): extinction 1/m and in-scattered display colour.
 */
export const VIS_UNBURNED_HAZE_KAPPA = 0.8;
export const VIS_UNBURNED_HAZE_COLOR: readonly [number, number, number] = [0.08, 0.09, 0.1];

/** Knock pressure-pattern overlay: amplitude that maps to full overlay strength, Pa (≈ heavy knock MAPO in a CFR engine). */
export const OVERLAY_KNOCK_FULL_SCALE_PA = 1.5e5;
/** Knock overlay opacity per metre at full scale. */
export const OVERLAY_KNOCK_KAPPA = 120;
/** Knock overlay colours for compression (p' > 0) and rarefaction (p' < 0), display units. */
export const OVERLAY_KNOCK_POSITIVE: readonly [number, number, number] = [0.25, 1.0, 0.4];
export const OVERLAY_KNOCK_NEGATIVE: readonly [number, number, number] = [0.2, 0.45, 1.0];

/**
 * Fraction of a light source's luminous intensity scattered into a glare
 * halo (camera lens / eye). Optical, not an intensity boost: the halo carries
 * this fraction of the same physical luminous intensity.
 */
export const VIS_GLARE_FRACTION = 0.1;
/** Glare halo radius in the scene, m. */
export const VIS_GLARE_RADIUS = 4e-3;

/**
 * Display persistence of flashes (s of WALL time), like an eye/sensor
 * afterimage. The breakdown lasts ~10 ns; its energy is integrated over the
 * frame's simulated exposure (dtWall·timeScale) and then decays with this
 * time constant so it can be seen at all.
 */
export const VIS_FLASH_PERSISTENCE_S = 0.08;

/** Minimum on-screen radius of the discharge column, m (radiance scaled down to conserve power). */
export const VIS_MIN_CHANNEL_RADIUS = 6e-5;

/** Mass represented by one tracer particle (sets particle count ∝ mass flow). Computed per engine from this fraction of the displaced charge mass. */
export const TRACER_PARCELS_PER_DISPLACED_CHARGE = 1500;
/** Max tracer particles. */
export const TRACER_CAPACITY = 4096;
/** Intake tracer lifetime in crank revolutions after entering (lives through combustion). */
export const TRACER_INTAKE_LIFE_REVS = 1.3;
/** Height of the port section tracers traverse above the valve seat, in valve seat diameters. */
export const TRACER_PORT_LENGTH_DIAMETERS = 1.5;
/** Tracer point size, m (cylinder frame). */
export const TRACER_SIZE = 1.1e-3;
/** Max sub-steps per frame for tracer integration. */
export const TRACER_MAX_SUBSTEPS = 48;
/** Max tracer sub-step, s (simulated). */
export const TRACER_MAX_DT = 2e-4;
