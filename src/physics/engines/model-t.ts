/**
 * Ford Model T engine (1909–1927), represented as a c. 1924–25 engine: starter-type block, the post-1917
 * "high head" with the September 1918 chamber revision (CR ≈ 3.98), light rods and 1 lb 12 oz pistons,
 * 3/4 in magnets on a ring-gear flywheel, Ford/K-W trembler coils, Ford roller timer, Holley NH /
 * Kingston L-4 class carburettor, cast-iron manifolds, thermosyphon cooling.
 *
 * Inline four, side valves ("L-head" / flathead), 3.750 × 4.000 in, 176.7 in³ (2.896 L), firing order
 * 1-2-4-3, three main bearings, crank throws 1&4 / 2&3 in one plane, no counterweights.
 *
 * Every number carries a citation. Numbers that could not be checked against a fetched source are
 * marked `UNVERIFIED:` with the reasoning behind the estimate. Short keys (full references in
 * MODEL_T.sources):
 *   [DB23]   Ford Motor Co., 1923 Ford Dealers' Data Book ("Details of Ford Engine Construction"),
 *            transcribed in B. McCalley, MTFCA Model T Encyclopedia
 *   [Good22] A. A. Good, Ford Car, Truck and Tractor Repair, McGraw-Hill 1922 (archive.org FordRepair1922)
 *   [FM19]   Ford Motor Co., Ford Manual for Owners and Operators of Ford Cars and Trucks, 1919
 *            (Project Gutenberg #46206)
 *   [Dyke24] Dyke's Automobile and Gasoline Engine Encyclopedia, 13th ed., 1924
 *   [Page29] V. W. Pagé, Models T and A Ford Class, 1929 (Clymer reprint)
 *   [FSB]    Ford Service Bulletins (1919–), reprinted in D. R. Post (ed.), Model T Ford Service Bulletin
 *            Essentials, Post Motor Books 1966 (archive.org bwb_S0-DTZ-218; read via full-text search)
 *   [McC]    B. McCalley, MTFCA Model T Encyclopedia, 'E' pages (Ford drawing-change records)
 *   [Fahn]   M. Fahnestock, Know Your Model A Ford (1958, reprinting his c. 1928 Q&A) and Fast Ford Handbook
 *   [AF15]   H. L. Arnold & F. L. Faurote, Ford Methods and the Ford Shops, 1915
 *   [Tulsa]  MTFC Tulsa (F. Houston, L. Young et al.), technical pages at tildentechnologies.com/mtfctulsa
 *            ('Head Design', 'Piston Position vs Crankshaft Angle', dyno summaries) — club measurements
 *   [BP]     T. Boggess & R. Patterson, Model T ignition articles; R. Patterson & S. Coniff, "The Model T Ford
 *            Ignition System & Spark Timing" (fordmodelt.net) — restorer measurements
 *   [ECCT]   Electronic Coil Checking Tool (ECCT) V12 manual; M. Kossor, Model T coil articles
 *   [Cool]   cool386.com coil-tester pages (oscillograms of Ford coils on 6 V / 12 V)
 *   [Upton]  G. B. Upton (Cornell), spark-advance tests on a Ford engine, J. SAE Aug. 1923
 *   [BoM]    US Bureau of Mines Technical Paper 328 (1923) gasoline surveys
 *   [Page18] V. W. Pagé, The Model T Ford Car, Truck and Conversion Sets, Henley 1918
 *   [FM14]   Ford Motor Co., Ford Manual (1914): mixture "12 to 14 parts air"
 *
 * Research notes and the raw evidence (including the adversarial cross-checks) are summarised in
 * test/fixtures/modelt_validation_README.md.
 *
 * ── Frame / orientation ─────────────────────────────────────────────────────
 * +z = front of the engine (fan, timing gears, timer); cylinder 1 is at the front. The crank turns
 * clockwise seen from the front [FM19 hand-cranking practice; same convention as the CFR]. All eight
 * valves stand in one row on the RIGHT side of the engine (driver's right) [Page18 pp. 57–59; DB23]:
 * looking at the engine from the front, the right side of the car is on the viewer's left, i.e. −x in
 * the world/cylinder frame (+y up, +z front, right-handed). Valve order from the front
 * E-I-I-E-E-I-I-E (two siamesed intake ports for cylinders 1+2 and 3+4 between four exhaust ports
 * [valvetrain research; strongly supported, never stated explicitly — UNVERIFIED]): cylinders 1 and 3
 * have the exhaust valve in front (+z) of the intake, cylinders 2 and 4 are mirrored (layout.mirrorZ).
 *
 * ── Chamber ─────────────────────────────────────────────────────────────────
 * Cylinder frame origin: bore axis, at the roof of the head cavity over the bore. The piston crown
 * rises 5/16 in ABOVE the block deck at TDC [FM19 A22; Good22 pp. 46–47; Page29 p. 255], and the
 * distance from crown to head at TDC is "about 1 inch" [FSB, via full-text snippet — single source].
 * The rest of the clearance volume is the valve pocket over the valves (head cavity, "spade-shaped"
 * chambers [Wikipedia; EngineLabs 2017]); the gasket opening is 18.3 in² per cylinder [Tulsa 'Head
 * Design'] = bore 11.04 in² + pocket ≈ 7.3 in². With the pocket outline below that gives a pocket roof
 * ≈ 12.9 mm above the deck (head-only chamber ≈ 284 cm³ above a 1.1 mm gasket, between the 294 cm³ Tulsa
 * measured on a 1917–18 high head and the ≈ 275 cm³ after Ford's 1/16 in chamber reduction of 9-25-18
 * [McC] — consistent).
 */
import { DEG, P_ATM } from '../core/constants';
import type {
  CylinderLayoutSpec,
  EngineSpec,
  TremblerMagnetoIgnitionSpec,
  VehicleSpec,
} from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import type { PnhFrictionInputs } from '../mechanics/friction';
import { oilViscosityCst, PNH_REFERENCE_VISCOSITY_CST } from '../mechanics/friction';

const IN = 0.0254;
const LB = 0.45359237;
const OZ = LB / 16;

// ─────────────────────────────────────────────────────────────────────────────
// Basic geometry
// ─────────────────────────────────────────────────────────────────────────────

/** Bore, m: 3.750 in [DB23; Good22 pp. 19–20; Dyke24 Instr. 70; Page29 p. 402]. */
const BORE = 3.75 * IN;
/** Stroke, m: 4.000 in (2.000 in throw) [DB23; AF15 p. ~197; Page29 p. ~31]. */
const STROKE = 4.0 * IN;
const R_BORE = BORE / 2;
const A_PISTON = (Math.PI * BORE * BORE) / 4;
/** Displaced volume per cylinder, m³ (723.96 cm³; 176.7 in³ total [McC 'E'; Fahn]). */
const VD = A_PISTON * STROKE;

/**
 * Compression ratio 3.98:1 for the post-1916 high head [McC 'E' (Ford records): "4.5:1 approx., 4.1:1
 * after 1912, 3.98:1 after 1916"]; Gunnell, Standard Catalog of Ford (2007): 4.0:1; Allen, Ford Model T
 * (1987): 4.0:1 for the 20 bhp engine. Period spread: Dyke24 "3.6"; Tulsa measured 3.8 on a 1917–18 high
 * head (294 cm³). The chamber build-up for the post-9/1918 head (294 − 18.7 cm³ chamber reduction [McC]
 * − 45.9 cm³ crown protrusion + 13.1 cm³ gasket [Tulsa]) gives Vc = 242.5 cm³ → CR 3.985, consistent.
 * Sensitivity range for validation: 3.78–3.98 (Vc 243–261 cm³).
 */
const CR = 3.98;
/** Clearance volume per cylinder, m³. */
const VC = VD / (CR - 1);

/**
 * Crevice volumes, m³.
 *  - UNVERIFIED: top land: radial clearance 0.005–0.006 in (top 0.010–0.012 in under bore [Page29 p. 402;
 *    Dyke24]) × an assumed 1/4 in land height × πB ≈ 0.27 cm³; top-ring-groove back volume: groove
 *    1/4 × 13/64 in [Dyke24] with an eccentric ring 0.180 → 0.085 in thick [Dyke24 c. 1917–20], ≈ 1.9 mm
 *    mean back gap, only partly connected → ≈ 1.7 cm³ counted. Quench crevice ≈ 2.0 cm³.
 *  - UNVERIFIED: plug cavity of the long-body 1/2 in pipe-thread plug [Dyke 1917 Ford Supplement] ≈ 0.5 cm³.
 */
const QUENCH_CREVICE = 2.0e-6;
const CREVICE = QUENCH_CREVICE + 0.5e-6;

/** Crown rise above the block deck at TDC, m: 5/16 in [FM19 A22; Good22 pp. 46–47; Page29 p. 255; Tulsa]. */
const CROWN_ABOVE_DECK = (5 / 16) * IN;
/**
 * Bore-column depth at TDC (head roof over the bore → crown), m: "between the piston and the head … is
 * about 1 inch" [FSB, full-text snippet]. UNVERIFIED (single snippet; consistent with 18.3 in² open area
 * and Vc 243–261 cm³ only if the pocket is 10–14 mm deep, which it is).
 */
const H_TDC = 1.0 * IN;
/** Deck (gasket face of the block = valve-seat plane), cylinder frame y, m. */
const DECK_Y = -(H_TDC + CROWN_ABOVE_DECK);

/**
 * Valve pocket plan (cylinder frame, valves on −x): rounded rectangle covering both valve heads with
 * ≈ 5 mm margin, overlapping the bore (the part inside the bore disc belongs to the bore column).
 * UNVERIFIED: valve centres x = −68 mm, z = ±24 mm (no drawing found; bore radius 47.6 + ≈ 4 mm wall +
 * valve head radius 18.7 mm → |r| ≥ 70 mm, valve pitch ≈ 46–51 mm leaving ≈ 57 mm between the siamesed
 * intakes of cylinders 1 and 2 at 4-1/8 in bore spacing). The resulting gasket opening (bore + pocket)
 * is 117.3 cm² vs 118.1 cm² (18.3 in²) measured [Tulsa 'Head Design'].
 */
const POCKET = { xMin: -0.0915, xMax: -0.03, zMin: -0.0475, zMax: 0.0475, cornerRadius: 0.02 };
/**
 * Plan area of POCKET minus the bore disc, m² — numerical (4000² grid) integration of the outline above;
 * the chamber model's exact area must agree (engines/model-t.test.ts).
 */
export const MODEL_T_POCKET_PLAN_AREA = 46.03020157589403e-4;
/** Pocket height above the deck that closes the clearance volume (Vc = A_p·H_TDC + pocket + crevice), m. */
const POCKET_HEIGHT = (VC - CREVICE - A_PISTON * H_TDC) / MODEL_T_POCKET_PLAN_AREA;

// ─────────────────────────────────────────────────────────────────────────────
// Valvetrain
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Running valve clearance, m. Ford set the cold push-rod-to-stem gap at 1/64–1/32 in by stem length
 * [FM19; DB23], .022–.032 in in the Ford Service Course, .022–.028 in check gap [FSB; Dyke24]. 0.0256 in
 * reproduces the 218° piston-position timing with the three-arc lobe (valvetrain research, verified).
 */
export const MODEL_T_VALVE_LASH = 0.0256 * IN;

/**
 * Cam lobe: forged-integral three-arc profile, no ramp, acting on a 1 in flat-footed mushroom push rod
 * [DB23; valvetrain research]. Radii and rise at the 4-decimal precision of the MTFC Tulsa table of the
 * Ford-drawing lobe (design_stock.htm Table 1; via D. R. Post 1997 and M. C. Turkish 1946): base radius
 * 0.4060 in, flank radius 1.2601 in, nose radius 0.0313 in, rise 0.2502 in (Ford: heel 13/16 in,
 * heel-to-toe 1-1/16 in → 0.250 in). Same lobe for intake and exhaust. With these values the exact
 * flat-follower lift (gas-exchange/cam-lift.ts; envelope oracle tools/reference/mechanics_cam_three_arc.py)
 * reproduces MTFC's duration-vs-lash Table 2 within 0.11° and gives 218.04° seat-to-seat at the
 * 0.0256 in lash (Ford's 218.1°/217.9° piston-position timing); the earlier 3-decimal rounding
 * (0.406/1.260/0.031/0.250 in) gave 217.91°.
 */
const CAM = { kind: 'three-arc-flat-follower' as const, baseRadius: 0.406 * IN, flankRadius: 1.2601 * IN, noseRadius: 0.0313 * IN, rise: 0.2502 * IN };

/** Net valve lift at the running clearance, m: 0.2502 in rise − 0.0256 in lash = 0.2246 in (5.70 mm). */
const VALVE_LIFT = CAM.rise - MODEL_T_VALVE_LASH;

/**
 * Valve: two-piece, cast-iron head 1.47–1.48 in pressed on a steel stem, same part (T-424-B) for intake
 * and exhaust; 45° face seating directly in the block around a 1-5/16 in port; stem 0.311 in in plain
 * 0.3125 in holes [DB23; Good22; valvetrain research, verified against drawings T-424-B/T-431].
 */
const VALVE_HEAD = 1.47 * IN;
const VALVE_PORT = (5 / 16 + 1) * IN;
const VALVE_STEM = 0.311 * IN;
/** Valve centres (cylinder 1; −x = right side of the car). UNVERIFIED, see POCKET. */
const VALVE_X = -0.068;
const VALVE_Z = 0.024;

// ─────────────────────────────────────────────────────────────────────────────
// Ignition: flywheel magneto → roller timer → four Ford/K-W trembler coils
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Ignition. Coil values are Ford 1916 engineering data for the K-W coil (standard 1913–27) via [BP]:
 * primary 212 turns, 0.295 Ω, 3.3 mH (secondary open) / 0.6 mH (secondary shorted); secondary 16,600
 * turns, 3,300 Ω, 22.0 H; condenser ".40 – .45 MFD" (Ford spec). RL back-solves of independent
 * oscillograms give 3.1–3.9 mH [Cool; ECCT], supporting 3.3 mH.
 *  - couplingCoefficient: 1 − k² = 0.6/3.3 → k ≈ 0.90 (primary side) — derived; the secondary-side ratio
 *    gives ≈ 0.70 (UNVERIFIED which the coupled-circuit model should use). Kept at 0.90: with the condenser
 *    and C₂ the open-circuit ring after a 6 V fire peaks at 412 V on the points / 34 kV on the plug
 *    (12 V: 512 V / 42 kV) in ignition/trembler-coil.ts and its Radau oracle — the ≈ 300–400 V primary
 *    peak [Cool] including reflected C₂ (≈ 400 V lossless estimate, research verification) and ≥ the
 *    8–20 kV Ford quotes for the HT output; k = 0.85 / 0.95 change the peak by < 7 % (not refitted). The
 *    model has no core loss, so the open-circuit peaks are upper bounds.
 *  - secondaryCapacitance: UNVERIFIED typical 20–50 pF (winding + HT lead + plug).
 *  - vibrator: point gap 1/32 in with the armature held down, cushion-spring clearance 0.005 in [FM19;
 *    ECCT; BP] → breakTravel 0.127 mm, maxTravel ≈ 0.127 + 0.794 mm. pullCurrent, naturalFrequency,
 *    dampingRatio and airGap are FITTED (UNVERIFIED as individual numbers: no mechanical data exist) with
 *    the 1-DOF armature model of ignition/vibrator.ts (preload deflection g₀/2, i.e. pullCurrent is the
 *    static pull-in current) by least squares (tools/reference/ignition_trembler_oracle.py --fit, Radau
 *    model of the full circuit: battery + 0.05 Ω, no timer, 10 kV bench gap) to the DC coil-tester data
 *    [ECCT; Cool]: first fire after the make 3.5 ms on 6 V, 2.5 ms on 9 V, 2.0 ms on 12 V, points re-closing
 *    1.8 ms after the 6 V fire. Achieved (oracle): 3.503 / 2.496 / 1.998 ms at 5.34 / 5.96 / 6.52 A
 *    (targets ≈ 5.0–5.4 / 6.2 / 6–7 A — the currents are not fitted, they follow from the RL ramp),
 *    re-close 1.813 ms; emergent: buzz period 5.0 ms on 6 V (≈ 200 Hz; ≈ 190 Hz [Cool]), points open longer
 *    at higher current (2.3 ms at 9 V, 2.7 ms at 12 V; "higher firing current throws the vibrator open
 *    wider" [Kossor]), fire at 4.1–4.3 A on slow 120–150 rpm magneto pulses (HCCT fires at 3.0–4.4 A
 *    [Kossor]), ½L₁I² = 47 mJ at the 6 V fire (48–50 mJ at 5.4–5.5 A [ECCT]). The fit is degenerate along
 *    (f_n ↓, ζ ↑, g₀ ↑) — e.g. 3.47 A / 113 Hz / 0.56 / 1.56 mm fits equally well and gives the same
 *    observables to < 1 % — the solution with the lower damping and g₀ just above maxTravel (armature
 *    resting ≈ 1 mm above the core) is used. A blade-spring estimate (≈ 0.5 mm × 13 mm × 40 mm spring steel,
 *    ≈ 2 g armature) gives f_n ≈ 125 Hz (UNVERIFIED, own estimate), consistent.
 *  - magneto: 16 magnets / 16 coils, 8 cycles per crank revolution [Dyke24 pp. 248–249; FSB magneto table];
 *    open-circuit ≥ 7 V at 400 rpm (Ford service minimum), healthy ≈ 10 V at 407 rpm … 24.4 V at 1220 rpm
 *    (restorer, rms) → ≈ 0.0235 V_rms/rpm = 0.317 V_peak per rad/s; source ≈ 0.3 Ω DC [BP]. internalInductance
 *    FITTED (UNVERIFIED: assumes the table's amperes are short-circuit currents): least squares of
 *    |R_s + jω_e L_s| to V/I of the 1976 Gas Engine Magazine table (W. H. Payne; 400–1200 rpm) with R_s = 0.3 Ω
 *    → 3.05 mH (residuals −0.18…+0.17 Ω of 1.24–2.91 Ω; 3.15 mH in relative terms); short-circuit current at
 *    speed k/(N L_s) = 13 A peak = 9.2 A rms vs 7.9–9.0 A in the table. phaseDeg FITTED to Patterson's 600 rpm
 *    magneto spark ladder [BP] with the whole system (ignition/trembler.test.ts; timer 0.1 Ω, plug at 8 kV):
 *    first spark at full retard (make 15.5° ATDC) 26.54° ATDC (measured 26.5°); stationary first-spark angles
 *    of the lever plateaus 26.39° ATDC, 3.89° ATDC, 18.61° BTDC, 41.11° BTDC (measured 26.5, 4, −18.5, −41°)
 *    → φ = −8.5°: the EMF zero crossing 8.5° BTDC, i.e. the magneto "leads TDC" by 8.5° vs the 7° of Ford
 *    drawing T-701-C [BP] — consistent to 1.5°.
 *  - timer: each segment grounded for 87° crank (15.5° → 102.5° ATDC at full retard, measured) [BP]; spark
 *    lever travel 80° (28 notches), make at 15.5° ATDC full retard … 64.5° BTDC full advance with Ford's
 *    2-1/2 in gauge setting (1919–27) [BP; Ford Service Par. 126]. contactResistance UNVERIFIED (roller).
 *  - battery: 6 V three-cell storage battery on starter cars from 1919 [FM19 starting-system booklet];
 *    internal resistance UNVERIFIED (also the source resistance assumed for the coil-tester fit above).
 */
export const MODEL_T_IGNITION: TremblerMagnetoIgnitionSpec = {
  type: 'trembler-magneto',
  coil: {
    primaryInductance: 3.3e-3,
    primaryResistance: 0.295,
    secondaryInductance: 22,
    secondaryResistance: 3300,
    secondaryCapacitance: 40e-12, // UNVERIFIED
    couplingCoefficient: 0.9,
    condenserCapacitance: 0.43e-6,
    vibrator: {
      pullCurrent: 3.53, // UNVERIFIED: fitted to the DC firing times (see above)
      naturalFrequency: 133, // UNVERIFIED: fitted
      dampingRatio: 0.34, // UNVERIFIED: fitted
      breakTravel: 0.005 * IN,
      maxTravel: 0.005 * IN + (1 / 32) * IN,
      airGap: 1.16e-3, // UNVERIFIED: fitted (effective gap; armature rests ≈ 0.24 mm beyond its travel to the core)
    },
  },
  magneto: {
    cyclesPerRevolution: 8,
    emfConstant: 0.317,
    phaseDeg: -8.5, // UNVERIFIED: fitted to the [BP] 600 rpm spark ladder (EMF zero 8.5° BTDC; T-701-C: 7°)
    internalResistance: 0.3,
    internalInductance: 3.05e-3, // UNVERIFIED: fitted to the 1976 table's V/I (short-circuit assumption)
  },
  battery: { voltage: 6, internalResistance: 0.05 },
  timer: { contactArcDeg: 87, advanceRangeDeg: [-15.5, 64.5], contactResistance: 0.1 },
};

// ─────────────────────────────────────────────────────────────────────────────
// Cylinder arrangement
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Bore spacing 4-1/8, 5-1/4, 4-1/8 in (#1–#2, #2–#3 across the centre main, #3–#4) — UNVERIFIED (MTFCA forum
 * thread t=29195, search snippet only; Model A 4-1/4, 5-5/16, 4-1/4). Engine centred on z = 0, cylinder 1
 * at the front (+z). Firing order 1-2-4-3 [DB23 'Firing Order'; FM19]: firing TDCs 180° apart →
 * offsets [0, 180, 540, 360] for cylinders 1..4; throws 1&4 at 0°, 2&3 at 180°.
 * Main bearings: front, centre (between #2 and #3), rear; crank 25-5/32 in long, mains 2, 2-3/16 and
 * 3-1/8 in long [DB23; Page29 p. 403]. UNVERIFIED: front/rear main centres ≈ 66 mm outboard of cylinders
 * 1/4 (half rod journal 19 mm + web ≈ 22 mm + half main bearing).
 */
const Z1 = (4.125 + 5.25 / 2) * IN;
const Z2 = (5.25 / 2) * IN;
export const MODEL_T_LAYOUT: CylinderLayoutSpec = {
  firingOrder: [1, 2, 4, 3],
  firingOffsetDeg: [0, 180, 540, 360],
  axisZ: [Z1, Z2, -Z2, -Z1],
  mirrorZ: [false, true, false, true],
  mainBearingZ: [Z1 + 0.066, 0, -Z1 - 0.066],
};

// ─────────────────────────────────────────────────────────────────────────────
// Vehicle (1919–25 touring car) for the 'vehicle' load model
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Vehicle: 1919 touring with starter and demountable rims ≈ 1650–1750 lb [FSB shipping-weight table] + two
 * occupants → UNVERIFIED 900 kg laden (road-test cars 850–1090 kg). Planetary gearbox low 2.75:1, high
 * 1:1, reverse 4.0:1 (triple gears 27/33/24, sun gears 27/21/30); spiral-bevel axle 11/40 = 3.636:1;
 * 30 × 3-1/2 in rear tyres, rolling radius ≈ 0.379 m (Ford: 40 rpm per mph in high [FSB dyno table]).
 * UNVERIFIED: driveline efficiency 0.8 (chassis vs engine dyno: 15.1 hp at the wheels vs ≈ 20 hp
 * [Tulsa]); C_rr 0.012 (0.01–0.015 high-pressure tyres); C_d·A 1.86 m² (open car, C_d 0.8 × 2.32 m²).
 * Checks with the 'vehicle' load model (mechanics/load.test.ts): 40.96 rpm/mph in high; reflected inertia
 * 9.8 kg m² in high (15× J_rot), 1.3 kg m² in low; top speed on Ford's WOT torque curve 44.1 mph at 1806 rpm
 * (target 42–45 mph at ≈ 1720–1845 rpm); with Tulsa's 850 kg / C_rr 0.01 it reproduces their Ford-curve
 * road calculation (44.6 / 34.3 / 28.7 mph level / 5 % / 6.7 % vs 45 / 35 / 30; 8.1 % maximum grade in high).
 * Wheel rotary inertia (≈ 4 % of the car's mass, UNVERIFIED) is not in VehicleSpec (LoadModel option).
 */
export const MODEL_T_VEHICLE: VehicleSpec = {
  mass: 900,
  wheelRadius: 0.379,
  finalDrive: 40 / 11,
  gears: { low: 2.75, high: 1, reverse: 4.0 },
  drivelineEfficiency: 0.8,
  rollingResistance: 0.012,
  dragArea: 1.86,
};

// ─────────────────────────────────────────────────────────────────────────────
// Friction (PNH model inputs)
// ─────────────────────────────────────────────────────────────────────────────

/** Crankcase oil: "light/medium" engine oil shared with the planetary transmission. UNVERIFIED grade/temperature. */
export const MODEL_T_OIL_GRADE = 'SAE30' as const;
/** UNVERIFIED: sump/splash-oil temperature ≈ 70 °C (no oil cooler, oil in the flywheel pit). */
export const MODEL_T_OIL_TEMPERATURE = 343.15;

/**
 * PNH friction inputs (mechanics/friction.ts). 3 babbitt mains, all journals 1.248 in; main lengths 2,
 * 2-3/16, 3-1/8 in (mean used); rod bearings 1.495–1.505 in long [DB23; Page29 p. 403]; 3 cam bearings,
 * 8 valves, flat (mushroom) followers acting directly on the stems [DB23]; 3 rings (2 compression + 1 oil)
 * = PNH's basis [Dyke24]. No oil pump, no water pump; the belt fan and (1919+) generator are the only
 * auxiliaries → auxiliaryFactor UNVERIFIED 0.3. valvetrain 'L-head' (mechanics/friction.ts): PNH has no
 * side-valve type; it maps to the direct-acting flat-follower constants (no rocker or pushrod) — UNVERIFIED
 * mapping; the Model T's light valve springs (24–28 lb installed) suggest even less. PNH was fitted on 1980s
 * engines at ≳ 1000 rpm: the FMEP of a babbitt-bearing 1920s engine at 400–2000 rpm is an extrapolation
 * (model-only). Result (WOT, p_i = p_a): 0.60 / 0.58 / 0.66 / 0.73 bar at 400 / 1000 / 1600 / 2000 rpm
 * (piston group ≈ 60–70 %), η_m 0.85 / 0.89 / 0.85 / 0.76 against Ford's WOT brake table. Consistent with that
 * table only if the WOT net IMEP is 0.33 / 0.45 / 0.37 / 0.26 of the ideal fuel–air-cycle IMEP (12.0 bar for a
 * full cylinder at CR 3.98, φ 1.15, 330 K, 0.95 bar; Cantera, tools/reference/cycle_fuel_air_oracle.py), e.g.
 * volumetric efficiency ≈ 0.45–0.65 at 400–1600 rpm with η_i ≈ 0.75 η_fa; the period ALAM/SAE rating
 * assumption (η_m 0.75 at 1000 ft/min [Good22 pp. 37–38]) would mean ≈ 1.37 bar at 1500 rpm, about twice PNH.
 * Not in PNH and included in Ford's data (measured at the transmission output [Tulsa]): churning of the
 * magneto flywheel and the planetary gear in the oil bath, and the wide (1/4 in) cast-iron rings. Decide any
 * extra loss against the cycle model's own WOT IMEP, not here (friction.test.ts prints the comparison).
 */
export const MODEL_T_FRICTION: PnhFrictionInputs = {
  bore: BORE,
  stroke: STROKE,
  cylinders: 4,
  compressionRatio: CR,
  mainBearings: { count: 3, diameter: 1.248 * IN, length: ((2 + 2.1875 + 3.125) / 3) * IN },
  rodBearings: { count: 4, diameter: 1.248 * IN, length: 1.5 * IN },
  camBearings: 3,
  valves: 8,
  maxValveLift: VALVE_LIFT,
  valvetrain: 'L-head',
  follower: 'flat',
  ringTensionFactor: 1,
  auxiliaryFactor: 0.3,
  viscosityRatio: oilViscosityCst(MODEL_T_OIL_GRADE, MODEL_T_OIL_TEMPERATURE) / PNH_REFERENCE_VISCOSITY_CST,
  boundaryRpm: 1000,
};

// ─────────────────────────────────────────────────────────────────────────────
// The spec
// ─────────────────────────────────────────────────────────────────────────────

/** UNVERIFIED: thermosyphon coolant ≈ 87 °C (circulation starts ≈ 82 °C, boils at 100 °C at sea level [FM19; Manly 1917]). */
export const MODEL_T_COOLANT_TEMPERATURE = 360;

export const MODEL_T: EngineSpec & { ignition: TremblerMagnetoIgnitionSpec; layout: CylinderLayoutSpec; vehicle: VehicleSpec } = {
  id: 'ford-model-t',
  name: 'Ford Model T (1924–25, high head)',
  sources: [
    "[DB23] Ford Motor Co., 1923 Ford Dealers' Data Book — 'Details of Ford Engine Construction' (cylinder case, crankshaft, connecting rods, pistons, valves, firing order), transcribed in B. McCalley, MTFCA Model T Encyclopedia (mtfca.com/model_t_encyclopedia, Wayback 2020).",
    '[Good22] A. A. Good, Ford Car, Truck and Tractor Repair, McGraw-Hill 1922 (archive.org/details/FordRepair1922): pp. 19–29 (dimensions), pp. 36 (magneto gap), pp. 46–47 (piston rise 5/16 in, timing).',
    '[FM19] Ford Motor Co., Ford Manual for Owners and Operators of Ford Cars and Trucks, 1919 (Project Gutenberg #46206): Answers 19 (compression), 22 (valve timing, piston 5/16 in above the casting), 31/60 (plugs, 1/32 in gap), 54/68 (magneto).',
    "[Dyke24] Dyke's Automobile and Gasoline Engine Encyclopedia, 13th ed. 1924 (archive.org dykesautomobileg0013edunse_v1v1): Instr. 70 (piston measurements, rings), magneto pp. 248–249, altitude/compression table, carburettor floats Instr. 85.",
    '[Page29] V. W. Pagé, Models T and A Ford Class, 1929 (Clymer reprint, archive.org bwb_S0-ELB-139): pp. 255 (valve timing), 402–405 (pistons, crankshaft, bearings, magneto voltage).',
    '[FSB] Ford Service Bulletins (1919–), in D. R. Post (ed.), Model T Ford Service Bulletin Essentials, 1966 (archive.org bwb_S0-DTZ-218): WOT horsepower/torque table (Fig. 84), magneto table, light-design piston/rod weights, ≈1 in piston-to-head distance, shipping weights.',
    "[McC] B. McCalley, MTFCA Model T Encyclopedia, 'E' (engine specifications; compression ratio by year; head drawing T-401C change record).",
    '[Fahn] M. Fahnestock, Know Your Model A Ford (1958) and Fast Ford Handbook (archive.org bwb_S0-CBG-705, bwb_S0-ELG-708): rod length 7 in, crankshaft 19 lb, flywheel 32 lb, powerplant 420 lb.',
    '[AF15] H. L. Arnold & F. L. Faurote, Ford Methods and the Ford Shops, Engineering Magazine Co. 1915: crank throw, flywheel T-701 35 lb finished.',
    "[Tulsa] MTFC Tulsa technical pages (tildentechnologies.com/mtfctulsa): 'Head Design' (chamber volumes, 18.3 in² open area, measured CR), 'Piston Position vs Crankshaft Angle', dyno summaries.",
    "[BP] T. Boggess & R. Patterson; R. Patterson & S. Coniff, 'The Model T Ford Ignition System & Spark Timing' (fordmodelt.net/downloads/Model%20T%20Ignition.pdf): Ford 1916 coil data, timer dwell and lever range, magneto phase and spark ladder.",
    '[ECCT] ECCT V12 manual; M. Kossor coil articles; [Cool] cool386.com/tester (Ford coil oscillograms: 3.5 ms at 6 V, 2.0 ms at 12 V).',
    '[Upton] G. B. Upton, Cornell spark-advance tests on a Ford engine, J. SAE, Aug. 1923 (MBT ≈ 0.108R/(1 + 0.001R) deg).',
    '[BoM] US Bureau of Mines TP 328 (1923): 1920–21 gasoline survey averages.',
    '[Page18] V. W. Pagé, The Model T Ford Car, Truck and Conversion Sets, Henley 1918 (archive.org modeltfordcartr00paggoog): pp. 57–59 (valves side by side on one side), pp. 80–81 (magneto), pp. 258–259 (copper-asbestos gaskets).',
    "[FM14] Ford Manual (1914): carburettor mixture '12 to 14 parts air'.",
    'Valve order E-I-I-E-E-I-I-E, cylinder bore spacing and the spark-plug position: MTFCA forum material (L. Young, "High Compression Head with Plugs over Pistons?", 8 Jul 2016: the stock plug location is over the valves) — lower trust.',
  ],
  cycle: 'four-stroke-si',
  cylinders: 4,
  layout: MODEL_T_LAYOUT,
  geometry: {
    bore: BORE,
    stroke: STROKE,
    conRodLength: 7.0 * IN, // [DB23] "7 in between centers"; [Good22 p. 27]; [Fahn]; [Tulsa] piston-position table back-solves to 7.000 in
    compressionRatio: CR,
    compressionRatioRange: [CR, CR], // fixed head
    pinOffset: 0, // UNVERIFIED: no offset known
    // UNVERIFIED: 1.94 in = block height 10-5/8 in [DB23] (assumed measured from the crank centreline)
    // + 5/16 in rise − 2 in throw − 7 in rod.
    compressionHeight: (10.625 + 0.3125 - 2 - 7) * IN,
    chamber: 'l-head',
    lHead: {
      deckY: DECK_Y,
      crownAboveDeckAtTDC: CROWN_ABOVE_DECK,
      pocket: { ...POCKET, roofY: DECK_Y + POCKET_HEIGHT },
    },
    creviceVolume: CREVICE,
    quenchCreviceVolume: QUENCH_CREVICE,
  },
  sparkPlug: {
    // Over the valves (stock heads; MTFCA, L. Young 2016 — lower trust), in the pocket between the two valve
    // heads. UNVERIFIED: x = −60 mm (between the bore edge and the valve row), gap 4 mm below the pocket roof.
    gapCenter: [-0.06, DECK_Y + POCKET_HEIGHT - 0.004, 0],
    gap: (1 / 32) * IN, // [FM19] "1/32 in, about the thickness of a smooth dime"; [Dyke 1917 Ford Suppl.]
    centerElectrodeDiameter: 2.5e-3, // UNVERIFIED (Champion X)
    groundElectrodeWidth: 2.5e-3, // UNVERIFIED
    axis: [0, -1, 0], // vertical, screwed into the top of the head [FM19 A60]
    threadDiameter: 0.84 * IN, // 1/2 in pipe thread (0.840 in OD) [Dyke 1917 Ford Suppl.; Dyke24]
  },
  intakeValve: {
    count: 1,
    headDiameter: VALVE_HEAD,
    seatInnerDiameter: VALVE_PORT,
    seatAngle: 45 * DEG,
    stemDiameter: VALVE_STEM,
    maxLift: VALVE_LIFT,
    // Piston-position timing [FM19 A22; DB23; Page29 p. 255] converted with L = 7, r = 2 in (verified):
    // intake opens 12.7° ATDC (piston 1/16 in down), closes 50.8° ABDC (piston 3-1/8 in below the casting).
    openDeg: -347.3,
    closeDeg: -129.2,
    timingLiftThreshold: 0,
    position: [VALVE_X, -VALVE_Z],
    shroudArcDeg: 0,
    shroudDirection: 0,
    seatY: DECK_Y,
    liftDirection: 1,
    lash: MODEL_T_VALVE_LASH,
    cam: CAM,
  },
  exhaustValve: {
    count: 1,
    headDiameter: VALVE_HEAD,
    seatInnerDiameter: VALVE_PORT,
    seatAngle: 45 * DEG,
    stemDiameter: VALVE_STEM,
    maxLift: VALVE_LIFT,
    // Exhaust opens 37.9° BBDC (piston 3-3/8 in below the casting), closes exactly at TDC (piston 5/16 in
    // above the casting) [FM19 A22; DB23; Page29 p. 255]: no overlap.
    openDeg: 142.1,
    closeDeg: -360,
    timingLiftThreshold: 0,
    position: [VALVE_X, VALVE_Z],
    shroudArcDeg: 0,
    shroudDirection: 0,
    seatY: DECK_Y,
    liftDirection: 1,
    lash: MODEL_T_VALVE_LASH,
    cam: CAM,
  },
  manifolds: {
    // UNVERIFIED: cast-iron Y manifold (1-1/8 in passages, ≈ 0.5 m) + siamesed block ports ≈ 0.5 L.
    intakeVolume: 0.5e-3,
    // UNVERIFIED: exhaust log ≈ 0.6 L + 56 in × 1.5 in OD pipe ≈ 1.35 L + three-shell muffler (outer 5 × 12 in)
    // ≈ 3.9 L [valvetrain research] → ≈ 6 L lumped.
    exhaustVolume: 6e-3,
    // UNVERIFIED: Holley NH throttle bore ≈ 1 in (not found); venturi = Holley G "strangling tube" cut to
    // 23/32 in c. 1920–22 (Holley G proxy for the NH) [induction research].
    throttleDiameter: 1.0 * IN,
    venturiDiameter: (23 / 32) * IN,
    throttle: { shaftRatio: 0.2, closedAngleDeg: 10, leakageArea: 3e-6 }, // UNVERIFIED generic plate
    exhaustOutletDiameter: 1.25 * IN, // UNVERIFIED: muffler restriction as an equivalent outlet
    intakePortDiameter: 1.125 * IN, // block manifold ports 1-1/8 in [valvetrain research]
    exhaustPortDiameter: 1.125 * IN,
    intakePortLength: 0.1, // UNVERIFIED: siamesed port cast in the block next to the jacket
  },
  masses: {
    // Light-design piston "approximately 1 lb. 12 oz." [FSB] (whether pin and rings are included is unknown)
    // + UNVERIFIED pin ≈ 0.12 kg + 3 rings × 42.5 g [FSB: 12 rings 1 lb 2 oz].
    piston: 1.75 * LB + 0.12 + 3 * (18 * OZ) / 12,
    // Light rod: Ford rod 1 lb 9 oz [Fahn]; new design 6–7 oz lighter than the old [FSB] → 0.60 kg (0.51–0.71).
    conRod: 0.6,
    conRodCgFromBigEnd: (7.0 * IN) / 3, // UNVERIFIED: ≈ 2/3 rotating, 1/3 reciprocating (typical of the period)
    // UNVERIFIED: 0.8 × the two-mass-equivalent inertia m·l_g·(l − l_g).
    conRodInertiaCg: 0.8 * 0.6 * ((7.0 * IN) / 3) * ((2 * 7.0 * IN) / 3),
    // UNVERIFIED: flywheel casting 32–35 lb [Fahn; AF15] + 16 magnets 13.6 lb + ring gear 4 lb (forum
    // weighings) ≈ 0.55 kg m², crank (19 lb [Fahn]) ≈ 0.008, triple gears and clutch drum on the flywheel
    // ≈ 0.09 → 0.65 kg m² (transmission in high gear; the vehicle load adds the car).
    rotatingInertia: 0.65,
  },
  walls: {
    // UNVERIFIED estimates (no measurements found): thermosyphon jacket ≈ 87 °C, cast-iron head/piston, low
    // specific output. The 'fixed' wall model shifts them with the coolant temperature.
    headTemperature: 460,
    pistonTemperature: 530,
    linerTemperature: 410,
    intakeValveTemperature: 480,
    exhaustValveTemperature: 850,
    intakePortTemperature: 360,
    blockTemperature: 450,
    referenceCoolantTemperature: MODEL_T_COOLANT_TEMPERATURE,
  },
  ignition: MODEL_T_IGNITION,
  vehicle: MODEL_T_VEHICLE,
};

// ─────────────────────────────────────────────────────────────────────────────
// Operating points
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Period gasoline as a PRF surrogate: c. 1920–22 US straight-run gasoline, estimated ≈ 40–55 MON; Ricardo's
 * 1921 HUCRs of commercial petrols (4.3–6.0) interpolate to ≈ ON 15–65, median ≈ 40–47 [induction research,
 * estimate]. UNVERIFIED: PRF 45.
 */
export const MODEL_T_PERIOD_GASOLINE_ON = 45;

/**
 * Default: driving the car in high gear on a level road at part throttle, magneto ignition. φ 1.15: Ford's
 * adjustment (enrich until top speed without black smoke ≈ best power) and the period 12–14 parts air
 * [FM14; induction research] — UNVERIFIED. Intake mixture ≈ ambient + stove heating − fuel evaporation
 * (≈ 22 K, Ricardo) → UNVERIFIED 300 K at 15.6 °C ambient.
 */
export const MODEL_T_CRUISE: OperatingPoint = {
  speedMode: 'free',
  rpm: 1000,
  loadTorque: 0,
  load: { kind: 'vehicle', gear: 'high', grade: 0 },
  throttle: 0.35,
  ambientPressure: P_ATM,
  ambientTemperature: 288.7,
  relativeHumidity: 0.5,
  intakeMixtureTemperature: 300,
  fuel: { kind: 'PRF', octaneNumber: MODEL_T_PERIOD_GASOLINE_ON },
  equivalenceRatio: 1.15,
  sparkAdvanceDeg: 25,
  dwellTime: 3e-3, // unused (trembler ignition)
  ignitionSource: 'magneto',
  compressionRatio: CR,
  egrFraction: 0,
  coolantTemperature: MODEL_T_COOLANT_TEMPERATURE,
};

/** Ford's 1918 engine-dyno condition: wide-open throttle, speed held by the brake [FSB Fig. 84]. */
export const MODEL_T_FORD_DYNO: OperatingPoint = {
  ...MODEL_T_CRUISE,
  speedMode: 'fixed',
  rpm: 1600,
  load: undefined,
  throttle: 1,
  // Upton MBT ≈ 0.108R/(1 + 0.001R) = 66° at 1600 rpm; the lever's full advance is 64.5° timer make and the
  // coil firing time retards the spark further — the dyno engine ran with the lever fully advanced.
  sparkAdvanceDeg: 64.5,
};

/**
 * Ford's WOT power/torque table, "representative of the motors in general use" (ratings "as high as
 * 22 1/2 horsepower" were obtained), measured at the transmission output [FSB Fig. 84, via full-text
 * snippets; independently tabulated by Tulsa/Sigworth]. rpm, torque lb-ft, horsepower. The 1400 rpm torque
 * is illegible (≈ 74); Ford's torque and hp columns are rounded independently (≤ 1.5 % disagreement).
 */
export const MODEL_T_FORD_WOT_TABLE: readonly (readonly [number, number, number])[] = [
  [300, 35, 2], [400, 57, 4.5], [500, 69, 6.5], [600, 73, 8.5], [700, 78, 10.4], [800, 81, 12.33],
  [900, 83, 14.2], [1000, 82, 15.6], [1100, 81, 16.66], [1200, 79, 18.2], [1300, 77, 19],
  [1400, 74, 19.66], [1500, 70, 20], [1600, 65, 20], [1700, 60, 19.4], [1800, 53, 18.2], [1900, 47, 17],
];

/** Upton's MBT spark advance on a Ford engine at WOT, crank deg BTDC (data to ≈ 1500 rpm) [Upton]. */
export function modelTUptonMbtDeg(rpm: number): number {
  return (0.108 * rpm) / (1 + 0.001 * rpm);
}
