/**
 * Waukesha CFR F-1/F-2 single-cylinder, variable-compression-ratio octane-rating
 * engine (ASTM D2699 Research / D2700 Motor method).
 *
 * Every number carries a citation. Numbers that could not be checked against a
 * fetched source are marked `UNVERIFIED:` with the reasoning behind the estimate.
 * Short keys used below (full references in CFR_F1.sources):
 *   [D2699]   ASTM D2699-15a, Standard Test Method for Research Octane Number (public.resource.org copy,
 *             archive.org item gov.law.astm.D2699.15A): Table 1, §§7–10, Annex A2, Table A4.1
 *   [D2700]   ASTM D2700-14, Standard Test Method for Motor Octane Number (archive.org item
 *             gov.law.astm.D2700.14): §10, Table 2, Table 3, Table A4.1
 *   [Pal18]   Pal et al., SAE Int. J. Engines 11(6), 2018, SAE 2018-01-0187 (OSTI 1572720), Tables 1–4, Figs. 3, 11
 *   [Choi18]  Choi, Kolodziej, Wallner & Hoth, SAE 2018-01-0848 (OSTI 1501884), Tables 1–3, Figs. 3–4
 *   [KW17]    Kolodziej & Wallner, Combustion Engines 171(4):164–169, 2017 (OSTI 1394801), Tables 1–3
 *   [Kal22]   Kalvakala et al., ASME J. Energy Resour. Technol. (OSTI 1962089), Table 1
 *   [W850]    Waukesha Motor Co. Engine Bulletin 850 (first CFR brochure, 1929), reprinted in the
 *             ASME National Historic Landmark brochure "The Waukesha CFR Fuel Research Engine" (1980)
 *   [CFRE]    CFR Engines Inc., "F1/F2 Octane Rating" product page (cfrengines.com/f1-f2)
 *   [Sin]     Sinpar (octane-engine vendor) pages quoting ASTM D2699/D2700 operating conditions
 *   [OSTI1583125] "Characterization of low temperature reactions in the standard CFR engine" (ANL), Table 1
 *   [SON23]   "Development of a Supercharged Octane Number and a Supercharged Octane Index", SAE 2023-01-0251
 *             (ANL, OSTI 1969816), Tables 1–2
 *
 * ── Frame / orientation ─────────────────────────────────────────────────────
 * The crank turns clockwise seen from +z (mechanics/kinematics.ts); ASTM specifies
 * "clockwise rotation of the crankshaft when observed from the front of the engine"
 * [D2699 §10.3.1], so +z is the FRONT of the engine. The intake-valve shroud faces
 * "toward the spark plug side of the combustion chamber and the swirl is directed in a
 * counterclockwise direction if it could be observed from the top of the cylinder"
 * [D2699 §10.2.5]. With the plug on +z this requires the intake valve on +x with its
 * shroud opening toward −z (charge leaves tangentially at x > 0 moving −z → angular
 * momentum along +y = CCW seen from above). [Pal18] Fig. 11 (top view of the X-ray-scanned
 * ANL chamber) shows the same arrangement up to a rotation: plug and knockmeter cavity on
 * one diameter, intake and exhaust valves on the perpendicular one, shroud on the plug side.
 * Which of the two mirror-consistent layouts (intake left or right seen from the front) the
 * real engine has is UNVERIFIED; this choice keeps the plug on the front (+z) side.
 *
 * ── How the compression ratio is varied ─────────────────────────────────────
 * "Adjustable 4:1 to 18:1 by cranked worm shaft and worm wheel drive assembly in cylinder
 * clamping sleeve" [D2699 Table 1]: the cylinder, whose head is integral ("cast iron with
 * flat combustion surface and integral coolant jacket" [D2699 Table 1]), is raised or
 * lowered in the clamping sleeve; an "open rocker assembly with linkage for constant valve
 * clearance as C.R. changes" [D2699 Table 1] keeps the valve lash. The piston, rod and
 * crank are untouched: only the clearance height at TDC changes. Cylinder height is read on
 * a digital counter (one digit = 0.0007 in of cylinder movement, basic setting 930
 * [D2699 A2.2.1.1, A2.2.2]) or a dial indicator (reading = 1.012 in − counter/1410
 * [D2699 Table A4 footnote]); see cfrCompressionRatioAtCounter.
 */
import { DEG, P_ATM } from '../core/constants';
import type { EngineSpec, IgnitionSystemSpec } from '../core/engine-spec';
import type { OperatingPoint } from '../core/operating-point';
import type { PnhFrictionInputs } from '../mechanics/friction';
import { oilViscosityCst, PNH_REFERENCE_VISCOSITY_CST } from '../mechanics/friction';
import {
  CHOI2018_EXHAUST_LIFT_START_DEG,
  CHOI2018_EXHAUST_LIFT_ZERO_LASH,
  CHOI2018_INTAKE_LIFT_START_DEG,
  CHOI2018_INTAKE_LIFT_ZERO_LASH,
} from './cfr-valve-lift';

const IN = 0.0254;
const BORE = 3.25 * IN; // [D2699] Table 1 "3.250 (standard)" in; [Pal18] Table 2 82.55 mm; [W850] "3¼ × 4½"
const STROKE = 4.5 * IN; // [D2699] Table 1 "4.50" in; [Pal18] 114.3 mm
const R_BORE = BORE / 2;

/** Default (nominal) compression ratio used when no operating point overrides it. */
const CR_NOMINAL = 7.0;

/** Hot running valve clearance (lash), m: 0.008 ± 0.001 in, both valves [D2699 §10.3.2.1]. */
export const CFR_VALVE_LASH = 0.008 * IN;

/**
 * Valve lift at zero lash produced by the standard cam lobes, m: "The resulting valve lift
 * shall be 0.238 in. ± 0.002 in." [D2699 §10.2.4, A2.1] (lobe contour rise 0.246–0.250 in, of
 * which the first 0.008–0.010 in is the quieting ramp that takes up the clearance). The ANL
 * engine measures 0.2363 / 0.2406 in (intake / exhaust) [Choi18] Fig. 3 (cfr-valve-lift.ts).
 */
export const CFR_VALVE_LIFT_ZERO_LASH = 0.238 * IN;

/**
 * Inductive-coil parameters. The CFR does NOT use a plain inductive coil: the standard
 * ignition is "electronically triggered condenser discharge through coil to spark plug"
 * [D2699 Table 1] ([KW17], [Pal18], [Choi18]: "Capacitive discharge coil to spark"). The
 * original 1929 engine used "standard battery, interrupter and coil" (Kettering) ignition
 * with a magneto as an option [W850]. The contract (IgnitionSystemSpec) only models an
 * inductive system, so these are generic transistorised-coil values; see
 * CFR_CDI_IGNITION for the capacitive-discharge description.
 */
const INDUCTIVE_IGNITION: IgnitionSystemSpec = {
  type: 'inductive',
  // UNVERIFIED: 12 V lead-acid system with charging (typical 13.5–14 V running).
  supplyVoltage: 13.5,
  // UNVERIFIED: typical transistorised-coil primary 3–8 mH; 4 mH with 1.5 Ω gives τ = L/R = 2.7 ms.
  primaryInductance: 4e-3,
  // UNVERIFIED: winding + driver resistance of a typical TCI coil.
  primaryResistance: 1.5,
  // UNVERIFIED: turns ratio ≈ 100 → L2 = L1·n² = 40 H.
  secondaryInductance: 40,
  // UNVERIFIED: typical secondary winding 5–10 kΩ (+ no suppressor on the CFR plug lead).
  secondaryResistance: 8e3,
  // UNVERIFIED: coil + HT lead + plug capacitance, typically 40–100 pF.
  secondaryCapacitance: 60e-12,
  // UNVERIFIED: typical closed-core ignition coil coupling.
  couplingCoefficient: 0.98,
  // UNVERIFIED: typical driver current limit; stored energy ½L1I² ≈ 98 mJ at 7 A.
  primaryCurrentLimit: 7,
  // UNVERIFIED: typical dwell; primary current reaches ≈ 6 A in 3 ms.
  dwellTime: 3e-3,
};

/**
 * Capacitive-discharge ignition actually fitted to CFR F1/F2 engines ([D2699] Table 1
 * "electronically triggered condenser discharge through coil to spark plug"; [KW17],
 * [Pal18], [Choi18]). PROPOSED CONTRACT EXTENSION (IgnitionSystemSpec is inductive-only).
 * All circuit values are UNVERIFIED generic CDI numbers; only the system type is sourced.
 */
export interface CapacitiveDischargeIgnitionSpec {
  type: 'capacitive-discharge';
  /** Storage capacitor, F. */
  storageCapacitance: number;
  /** Capacitor charge voltage, V. */
  chargeVoltage: number;
  /** Coil primary inductance, H. */
  primaryInductance: number;
  /** Coil secondary inductance, H. */
  secondaryInductance: number;
  /** Secondary resistance, Ω. */
  secondaryResistance: number;
  /** Lumped secondary capacitance, F. */
  secondaryCapacitance: number;
  /** Magnetic coupling coefficient. */
  couplingCoefficient: number;
}

export const CFR_CDI_IGNITION: CapacitiveDischargeIgnitionSpec = {
  type: 'capacitive-discharge',
  storageCapacitance: 1.0e-6, // UNVERIFIED: typical CDI 0.5–2 µF
  chargeVoltage: 400, // UNVERIFIED: typical 300–450 V → ½CV² = 80 mJ stored
  primaryInductance: 0.3e-3, // UNVERIFIED: CDI coils have low primary inductance
  secondaryInductance: 3, // UNVERIFIED: turns ratio ≈ 100
  secondaryResistance: 5e3, // UNVERIFIED
  secondaryCapacitance: 60e-12, // UNVERIFIED
  couplingCoefficient: 0.98, // UNVERIFIED
};

/** 100 °C jacket, K (see CFR_COOLANT_TEMPERATURE; needed before CFR_F1). */
const CFR_COOLANT_TEMPERATURE_K = 100 + 273.15;

/**
 * Lumped wall model resistances, surface → coolant, K/W (fitted, see CFR_F1.walls). Fixer round 2
 * fit on the round-2 calibrated model (cycle/calibration.ts); UNVERIFIED as physical conductances.
 */
export const CFR_WALL_THERMAL_RESISTANCE = Object.freeze({
  head: 0.2742,
  piston: 0.1121,
  liner: 0.03526,
  intakeValve: 0.5162,
  exhaustValve: 0.8570,
});

export const CFR_F1: EngineSpec = {
  name: 'Waukesha CFR F-1/F-2 octane-rating engine',
  sources: [
    '[D2699] ASTM D2699-15a, "Standard Test Method for Research Octane Number of Spark-Ignition Engine Fuel" (incorporated by reference; public copy archive.org/details/gov.law.astm.D2699.15A) — Table 1 (cast-iron cylinder with flat combustion surface and integral jacket; CR 4:1–18:1 by worm shaft/worm wheel in the clamping sleeve; 3.250 × 4.50 in, 37.33 in³; stellite-faced intake valve with 180° shroud, plain exhaust valve; cast-iron flat-top piston; 1 chrome/ferrous + 3 ferrous compression rings + 1 oil ring; camshaft overlap 5°; 9/16 in venturi; condenser-discharge ignition; constant 13° btdc), §7.1 (V-belts to a power-absorption motor for constant speed; thermal-syphon jacket), §8.2 (SAE 30 oil, 9.3–12.5 cSt at 100 °C), §10.2–10.3 (600 ± 6 rpm; valve timing IVO 10.0 ± 2.5° atdc, IVC 34° abdc, EVO 40° bbdc, EVC 15.0 ± 2.5° atdc; valve lift 0.238 ± 0.002 in; shroud toward the spark plug, CCW swirl from above; clockwise rotation from the front; hot valve clearance 0.008 ± 0.001 in; oil 172–207 kPa, 57 ± 8 °C; jacket 100 ± 1.5 °C; IAT 52 ± 1 °C at 101.0 kPa; humidity 0.00356–0.00712 kg/kg dry air; crankcase 25–150 mm H2O vacuum; Champion D16 plug, gap 0.51 ± 0.13 mm), Annex A2 (cam lobe 0.248 in, quieting ramps at 0.008–0.010 in; timing check 0.054 in lifter rise at 30 ± 2°; 0.0007 in per counter digit; basic counter 930), Table A4.1 (guide table: RON 90 → counter 726).',
    '[D2700] ASTM D2700-14, "Standard Test Method for Motor Octane Number of Spark-Ignition Engine Fuel" (public copy archive.org/details/gov.law.astm.D2700.14) — §10.3.6–10.3.7 (IAT 38 ± 2.8 °C, mixture 149 ± 1 °C), §10.3.16 and Table 3 (spark 26° btdc at counter 264 … 14° btdc at counter 1145, linear), Table 2 (9/16 in venturi from sea level to 1600 ft), Table A4.1 (guide table, 9/16 in venturi: MON 90 → counter 749).',
    '[Pal18] Pal, P., Kolodziej, C., Choi, S., Som, S., Broatch, A., Gomez-Soriano, J., Wu, Y., Lu, T., See, Y.C., "Development of a Virtual CFR Engine Model for Knocking Combustion Analysis", SAE Int. J. Engines 11(6), 2018, SAE 2018-01-0187, OSTI 1572720 — Table 1 (MON spark 19–26° btdc), Table 2 (82.55 x 114.3 mm, con-rod 254 mm, 180° intake shroud non-rotating, rotating unshrouded exhaust valve, +5 CAD overlap, CD ignition), Table 3 (GT-Power FE wall temperatures), Table 4 (φ = 1.124, λ = 0.89 peak knock for iso-octane, CR 7.55), Fig. 3 (side-mounted plug and plug cavity, knockmeter cavity), Fig. 11 (top-view cut through the electrode: gap recessed ≈ 3.5 mm outside the bore in the plug passage, knockmeter cavity centre ≈ 3 mm inside the bore, valve centres ≈ ±23 mm; our pixel measurements, ±1 mm).',
    '[Choi18] Choi, S., Kolodziej, C., Wallner, T., Hoth, A., "Development and Validation of a Three Pressure Analysis (TPA) GT-Power Model of the CFR F1/F2 Engine for Estimating Cylinder Conditions", SAE 2018-01-0848, OSTI 1501884 — Table 2 (measured IVO/IVC 376/−152, EVO/EVC 141/373 CAD aTDC at the 0.008 in lash; valve reference diameters 1.346/1.356 in), Fig. 3 (zero-lash lift profiles, digitised in cfr-valve-lift.ts), Fig. 4 (CR vs digital counter), Table 3 (piston height 120.7 mm, cylinder length 140 mm, cast iron).',
    '[KW17] Kolodziej, C., Wallner, T., "Combustion characteristics of various fuels during research octane number testing on an instrumented CFR F1/F2 engine", Combustion Engines 171(4):164–169, 2017, doi:10.19206/CE-2017-427, OSTI 1394801 — Table 1 (RON/MON conditions), Table 2 (3.25 x 4.5 in, 37.33 in³, CR 4–18, 5 piston rings, CD ignition), Table 3 (worm gear CR adjustment, D-1 knock pickup + 501-C meter).',
    '[Kal22] Kalvakala, K., Pal, P., Wu, Y., Kukkadapu, G., et al., "Numerical analysis of fuel effects on advanced compression ignition using a cooperative fuel research engine CFD model", ASME J. Energy Resour. Technol., OSTI 1962089 — Table 1 (IVO 10° ATDC, IVC 34° ABDC, EVO 40° BBDC, EVC 15° ATDC).',
    '[W850] Waukesha Motor Co., Engine Bulletin 850 (1929) as reprinted in ASME, "The Waukesha CFR Fuel Research Engine — An International Historic Mechanical Engineering Landmark" (1980): worm-and-screw cylinder raising, cylinder–spark-advance linkage, cast-iron piston with five rings, sleeve bearings "double the normal dimensions", steam (evaporative) cooling ≈ 208 °F, battery/interrupter/coil ignition (magneto optional), belt drive to an induction motor with synchronous-motor constant-speed characteristic.',
    '[CFRE] CFR Engines Inc., F1/F2 Octane Rating product page (cfrengines.com/f1-f2): RON 600 rpm ±1 %, 13° BTDC; MON 900 rpm; jacket 100 ± 1.5 °C; oil 57 ± 8 °C, 172–207 kPa; synchronous motor 220/380/440 V 3-phase.',
    '[Sin] Sinpar pages quoting ASTM D2699/D2700 (valve timing, clearances, operating conditions) — superseded by [D2699]/[D2700].',
    '[OSTI1583125] Characterization of low temperature reactions in the standard CFR engine (ANL) — Table 1: MON spark timing 14–26° BTDC (varies with CR).',
    '[SON23] Development of a Supercharged Octane Number and a Supercharged Octane Index, SAE 2023-01-0251 (ANL), OSTI 1969816 — Table 1 (RON/MON conditions; MON spark timing variable with CR), Table 2 (displacement 0.612 L, one shrouded non-rotating intake valve, 5 CAD overlap, pancake chamber, side-mounted spark plug).',
    'Champion D16 spark plug: 18 mm thread, 12.7 mm reach (retailer catalogue data; the plug type itself is [D2699] §10.3.16).',
    'Patton, Nitschke & Heywood 1989 (SAE 890836) friction model as reproduced in Sandoval, D., "An Improved Friction Model for Spark Ignition Engines", B.S. thesis, MIT 2002 (Sandoval & Heywood 2003, SAE 2003-01-0725).',
    'NIST Chemistry WebBook: enthalpy of vaporization of 2,2,4-trimethylpentane (35.1 ± 0.2 kJ/mol) and n-heptane (36.4 kJ/mol at 303 K) — used for the RON mixture-temperature estimate.',
    'Bestel, D., Windom, B., et al., "3-D Modeling of the CFR Engine for the Investigation of Knock on Natural Gas", WSSCI Fall 2019 (OSTI 1597460): CFD geometry from X-ray scan incl. knockmeter cavity, J-gap spark plug, crevice.',
  ],
  cycle: 'four-stroke-si',
  cylinders: 1,
  geometry: {
    bore: BORE,
    stroke: STROKE,
    conRodLength: 0.254, // [Pal18] Table 2, [Choi18] Table 1, [Kal22] Table 1: 254 mm (= 10.0 in)
    compressionRatio: CR_NOMINAL, // nominal default only; the operating point sets the CR
    compressionRatioRange: [4, 18], // [D2699] Table 1 "Adjustable 4:1 to 18:1"; [KW17], [Pal18] Table 2
    // UNVERIFIED: no source mentions a wrist-pin offset; the CFR is a symmetric, non-offset design as far as known.
    pinOffset: 0,
    // UNVERIFIED (render only): piston is 120.7 mm tall [Choi18] Table 3 with five rings [D2699][KW17][W850]
    // above the pin; compression height estimated ≈ 58 % of the piston height.
    compressionHeight: 0.07,
    // "Cast iron with flat combustion surface" [D2699] Table 1; "Cast iron, flat 'pancake'" [KW17][Pal18];
    // piston "Cast iron, flat top" [D2699] Table 1.
    chamber: 'flat-disc',
    // UNVERIFIED: lumped dead volume ≈ 1.5 cm³ = top-land crevice (πB × ~8 mm × ~0.25 mm ≈ 0.52 cm³,
    // tight cast-iron-in-cast-iron fit) + ring-groove back volume (~0.2 cm³) + spark-plug and
    // knockmeter-port cavities connected to the chamber ([Pal18] Figs. 3, 11, Bestel 2019), ~0.8 cm³.
    creviceVolume: 1.5e-6,
    // UNVERIFIED: the narrow (flame-quenching) part of the above — top land + ring grooves ≈ 0.52 + 0.2
    // cm³ (same estimates); the plug and pickup cavities (≈ 0.8 cm³) are wide enough for the flame.
    quenchCreviceVolume: 0.72e-6,
  },
  sparkPlug: {
    // Side-mounted plug ([SON23] Table 2; [Pal18] Fig. 3 "spark plug & cavity") on the front (+z) side,
    // opposite the knockmeter port (see file header for the orientation argument).
    // [Pal18] Fig. 11 (cut through the electrode) shows the J-gap RECESSED ≈ 3.5 mm radially outside the
    // bore, inside the plug passage. The contract's flat-disc chamber (and FlameGeometry, which requires the
    // centre inside the bore) has no plug cavity, so the modelled gap centre is put just inside the bore
    // wall (1 mm): the kernel then touches the wall almost immediately, as the real recessed kernel does.
    // UNVERIFIED: depth below the fire deck (2 mm; must stay < h_TDC(18) = 6.4 mm).
    gapCenter: [0, -0.002, R_BORE - 0.001],
    // [D2699] §10.3.16.1, [D2700] §10.3.17.1: gap 0.51 ± 0.13 mm (0.020 ± 0.005 in).
    gap: 0.020 * IN,
    // UNVERIFIED: typical non-projected 18 mm plug (Champion D16-type) nickel centre electrode.
    centerElectrodeDiameter: 2.5e-3,
    // UNVERIFIED: typical J-gap ground strap width.
    groundElectrodeWidth: 2.5e-3,
    // Horizontal plug pointing radially inward from the +z side ([Pal18] Figs. 3, 11: axis horizontal, radial).
    axis: [0, 0, -1],
    // Champion D16 [D2699] §10.3.16: 18 mm thread per the Champion catalogue (retailer data).
    threadDiameter: 18e-3,
  },
  intakeValve: {
    count: 1,
    // UNVERIFIED: seat inner (reference) diameter + 2 × 1/32 in seat face; must fit the bore with the
    // valve centre 23 mm off-axis.
    headDiameter: 1.346 * IN + 0.0625 * IN,
    // [Choi18] Table 2 "valve reference diameter" (basis of the flow-bench flow coefficients): intake 1.346 in.
    seatInnerDiameter: 1.346 * IN,
    seatAngle: 45 * DEG, // UNVERIFIED: conventional 45° seat
    stemDiameter: (11 / 32) * IN, // UNVERIFIED: typical 11/32 in stem
    // Running valve lift = zero-lash lift 0.238 in [D2699 §10.2.4] − hot clearance 0.008 in [D2699 §10.3.2.1]
    // = 0.230 in (ANL engine: 0.2363 − 0.008 = 0.228 in, [Choi18] Fig. 3). NOT the 0.246–0.250 in cam-lobe rise.
    maxLift: CFR_VALVE_LIFT_ZERO_LASH - CFR_VALVE_LASH,
    // Actual seat-off / seat-on at the 0.008 in hot clearance, measured on the ASTM-timed ANL engine:
    // IVO 376 CAD ≡ −344°, IVC −152° ([Choi18] Table 2; our digitisation of its Fig. 3 gives −344.2 / −152.9 at
    // 0.008 in, and 0.054 in at 30.6° after TDC vs the ASTM check 30 ± 2° [D2699 A2.1.2]). The ASTM nominal
    // events (IVO 10° atdc, IVC 34° abdc; CFR_VALVE_TIMING_ASTM) are cam-lobe events on the quieting ramp
    // (0.008–0.010 in rise [D2699 A2.1]), ≈ 6° wider than the real valve opening at running clearance.
    openDeg: -344,
    closeDeg: -152,
    // Timing above is valve seat-off/seat-on with the running clearance (lift measured at the valve).
    timingLiftThreshold: 0,
    // Valve centres on the x diameter, ≈ 23 mm off-axis ([Pal18] Fig. 11, pixel measurement ±1 mm);
    // intake on +x — see file header for the side (UNVERIFIED mirror choice).
    position: [0.023, 0],
    // [D2699] §10.2.5: "180° shroud or protrusion just inside the valve face to direct the incoming fuel-air
    // charge and increase the turbulence"; the stem is pinned in a guide slot so the valve cannot rotate.
    shroudArcDeg: 180,
    // Masked half toward the spark plug (+z) [D2699 §10.2.5], so the opening faces −z → CCW swirl from above.
    shroudDirection: -Math.PI / 2,
  },
  exhaustValve: {
    count: 1,
    headDiameter: 1.356 * IN + 0.0625 * IN, // UNVERIFIED: as intake
    seatInnerDiameter: 1.356 * IN, // [Choi18] Table 2: exhaust reference diameter 1.356 in
    seatAngle: 45 * DEG, // UNVERIFIED
    stemDiameter: (11 / 32) * IN, // UNVERIFIED
    maxLift: CFR_VALVE_LIFT_ZERO_LASH - CFR_VALVE_LASH, // as intake ([Choi18] Fig. 3: 0.2406 − 0.008 = 0.233 in)
    // Seat-off/seat-on at the 0.008 in clearance: EVO 141°, EVC 373 CAD ≡ −347° ([Choi18] Table 2; digitised
    // Fig. 3: 142.2 / 371.9). ASTM nominal: EVO 40° bbdc (140°), EVC 15.0 ± 2.5° atdc (375°) [D2699 §10.2.3.2].
    // NB: at running clearance the measured valve events do not overlap (EVC 3° before IVO); the ASTM
    // "camshaft overlap 5°" [D2699 Table 1] is between the nominal cam events.
    openDeg: 141,
    closeDeg: -347,
    timingLiftThreshold: 0,
    position: [-0.023, 0], // [Pal18] Fig. 11 (≈ 23 mm), opposite the intake
    shroudArcDeg: 0, // "Stellite faced, plain type without shroud" [D2699 Table 1]; "No shroud, rotating" [Pal18][KW17]
    shroudDirection: 0,
  },
  manifolds: {
    // UNVERIFIED: carburettor-to-port volume incl. the MON mixture-heater housing. Round 1: 1 L (estimate);
    // fixer round 2: 0.25 L — the 0-D plenum then reproduces the measured intake-port pressure of [Choi18]
    // Fig. 2 (PRF98 RON: RMS 0.014 vs 0.025 bar at 1 L; dip 0.854 bar at −269° vs 0.838 at −266° measured;
    // measured-data validator r2 grid 0.25/0.5/1 L) while keeping the ASTM D2700 venturi-size effect.
    intakeVolume: 0.25e-3,
    // UNVERIFIED: volume. The ASTM exhaust line has a surge tank ([D2699] §10.3.11 "exhaust surge tank";
    // [Choi18]); estimated 10 L (≫ displacement, damps pulsations).
    exhaustVolume: 10e-3,
    // The CFR has NO throttle (always wide open); the only restriction is the carburettor venturi:
    // "A 9/16 in. (14.3 mm) venturi throat size shall be used regardless of ambient barometric pressure"
    // [D2699 §10.2.6] (MON: 9/16 in from sea level to 1600 ft, [D2700] Table 2).
    throttleDiameter: (9 / 16) * IN,
    exhaustOutletDiameter: 1.5 * IN, // UNVERIFIED
    intakePortDiameter: 1.25 * IN, // UNVERIFIED: slightly smaller than the 1.346 in valve reference diameter
    exhaustPortDiameter: 1.25 * IN, // UNVERIFIED
    // UNVERIFIED: heated length of the intake passage in the (water-cooled) head, estimated 0.1 m. Only
    // the product (port heat-transfer multiplier) × length matters (cycle/calibration.ts).
    intakePortLength: 0.1,
  },
  masses: {
    // UNVERIFIED: cast-iron piston "unusually long, and of heavy section" [W850], 120.7 mm tall [Choi18],
    // five rings [D2699]: shell 82.5 mm OD × 5 mm wall × 120.7 mm + crown + bosses ≈ 240 cm³ × 7.2 g/cm³
    // ≈ 1.7 kg, + wrist pin ≈ 0.15 kg + rings ≈ 0.1 kg.
    piston: 1.9,
    // UNVERIFIED: forged-steel rod, 254 mm centres, "long connecting rod" [W850]; typical 1.2–1.8 kg.
    conRod: 1.6,
    // UNVERIFIED: CG at ≈ 30 % of the rod length from the big end (typical).
    conRodCgFromBigEnd: 0.076,
    // UNVERIFIED: I_g ≈ 0.8·m_r·l_g·(l − l_g) (typical rods are slightly below the two-mass-equivalent inertia).
    conRodInertiaCg: 0.017,
    // UNVERIFIED: crank + flywheel ≈ 1.2 kg m² (≈ 0.4 m × 60 mm steel flywheel) + power-absorption motor
    // rotor reflected through the V-belts (J_motor·ratio², ≈ 0.03 × 3² at 600 rpm). In the octane tests
    // the belt-coupled motor holds the speed (600 ± 6 rpm, [D2699] §7.1, §10.2.1; [W850]; [CFRE]), so
    // this matters only in 'free' mode.
    rotatingInertia: 1.5,
  },
  walls: {
    // [Pal18] Table 3: GT-Power finite-element wall temperatures, RON conditions, PRF100 (iso-octane),
    // spark −13 CAD aTDC (standard knock). At spark 5 CAD aTDC: liner 426.8, head 493.2, piston 461.8,
    // intake valve 440.4, exhaust valve 501.9 K.
    headTemperature: 525.2,
    pistonTemperature: 485.2,
    linerTemperature: 430.1,
    intakeValveTemperature: 464.9,
    exhaustValveTemperature: 519.6,
    // [Pal18] Table 3 T_InPort (spark −13; 382.9 K at spark +5).
    intakePortTemperature: 387.8,
    // The temperatures above are the REFERENCE state (Pal/Choi PRF100, CR 7.55, λ 0.89, spark −12.72°,
    // 600 rpm, 100 °C jacket) of the lumped wall model (cycle-model.ts updateWalls): surface =
    // coolant + R·Q̄ with Q̄ the cycle-mean gas-to-surface heat flow. R fitted so that the calibrated
    // cycle model reproduces [Pal18] Table 3 at that state (tmp/fix2/fitwalls.ts, fixer round 2).
    // Check at spark +5 (Pal FE: head 493.2, piston 461.8, liner 426.8, IV 440.4, EV 501.9 K): see
    // cycle/calibration.ts — the model's heat load falls faster with retard than Pal's FE walls (no
    // exhaust-port or friction heat paths in the lumped model).
    referenceCoolantTemperature: CFR_COOLANT_TEMPERATURE_K,
    thermalResistance: CFR_WALL_THERMAL_RESISTANCE,
  },
  ignition: INDUCTIVE_IGNITION,
};

// ─────────────────────────────────────────────────────────────────────────────
// Geometry helpers (variable compression ratio)
// ─────────────────────────────────────────────────────────────────────────────

const areaOf = (spec: EngineSpec): number => (Math.PI * spec.geometry.bore ** 2) / 4;

/** Exact piston travel (TDC→BDC) of a spec's slider crank, m (= stroke without offset). */
function pistonTravel(spec: EngineSpec): number {
  const a = spec.geometry.stroke / 2;
  const l = spec.geometry.conRodLength;
  const e = spec.geometry.pinOffset;
  return Math.sqrt((l + a) ** 2 - e * e) - Math.sqrt((l - a) ** 2 - e * e);
}

/** Displaced volume of one cylinder, m³ (CFR: 611.7 cm³ = 37.33 in³, [D2699] Table 1, [KW17] Table 2). */
export function displacedVolume(spec: EngineSpec = CFR_F1): number {
  return areaOf(spec) * pistonTravel(spec);
}

/**
 * Flat-disc clearance height at TDC for compression ratio CR, m:
 *   h_TDC = (V_d/(CR − 1) − V_crevice)/A_p
 * (the crevice volume is part of the clearance volume, see mechanics/kinematics.ts).
 */
export function clearanceHeightAtTDC(compressionRatio: number, spec: EngineSpec = CFR_F1): number {
  return (displacedVolume(spec) / (compressionRatio - 1) - spec.geometry.creviceVolume) / areaOf(spec);
}

/**
 * World-frame height (along +y, crank axis at the origin) of the head fire-deck
 * face for compression ratio CR, m:
 *   y_head = y_pin,TDC + compressionHeight + h_TDC(CR),  y_pin,TDC = √((l + a)² − e²).
 * Use it to place the cylinder + head assembly in the renderer.
 */
export function headFacePosition(compressionRatio: number, spec: EngineSpec = CFR_F1): number {
  const g = spec.geometry;
  const a = g.stroke / 2;
  const yPin = Math.sqrt((g.conRodLength + a) ** 2 - g.pinOffset ** 2);
  return yPin + g.compressionHeight + clearanceHeightAtTDC(compressionRatio, spec);
}

/**
 * How far the cylinder/head assembly is raised at compression ratio CR relative
 * to its position at `referenceCR` (default: the maximum CR, i.e. the lowest
 * cylinder position), m (≥ 0 for CR ≤ referenceCR).
 */
export function cylinderRaise(compressionRatio: number, referenceCR?: number, spec: EngineSpec = CFR_F1): number {
  const ref = referenceCR ?? spec.geometry.compressionRatioRange[1];
  return clearanceHeightAtTDC(compressionRatio, spec) - clearanceHeightAtTDC(ref, spec);
}

/** Inverse of clearanceHeightAtTDC: compression ratio for a disc clearance height h (m). */
export function compressionRatioForClearanceHeight(h: number, spec: EngineSpec = CFR_F1): number {
  return 1 + displacedVolume(spec) / (areaOf(spec) * h + spec.geometry.creviceVolume);
}

// ─────────────────────────────────────────────────────────────────────────────
// Cylinder-height digital counter ↔ compression ratio
// ─────────────────────────────────────────────────────────────────────────────

/** Cylinder-height movement per digital-counter digit, m: 0.0007 in [D2699 A2.2.1.1]. */
export const CFR_COUNTER_STEP = 0.0007 * IN;

/** Basic (indexing) digital-counter reading [D2699 A2.2.2; D2700 §10.3.18.1]. */
export const CFR_BASIC_COUNTER = 930;

/**
 * Effective clearance height V_c/A_p at the basic counter reading 930, m. Least-squares fit
 * of the rigid-raise model V_c(c)/A_p = h₉₃₀ − 0.0007 in·(c − 930) to the oil-measured
 * clearance volumes of [Choi18] Fig. 4 (their polynomial, counter 400–1400):
 * tools/reference/mechanics_cfr_conditions.py (fixture counter.rigidRaise). Max height residual
 * 0.9 mm at the range ends; ≤ 0.02 CR from the polynomial at the RON/MON PRF90 settings.
 * Engine-specific (ANL CFR), but the counter is indexed per ASTM (basic compression pressure).
 */
export const CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER = 0.01742116238087605;

/**
 * Compression ratio at a (standard-barometer, uncompensated) digital-counter reading, rigid-raise
 * model: V_c(c) = A_p (h₉₃₀ − 0.0007 in·(c − 930)), CR = 1 + V_d/V_c. Exact kinematics of the worm
 * drive [D2699 A2.2.1.1] with one fitted constant (CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER); valid
 * for any counter reading (the ASTM tables span 264–1290), unlike the polynomial fit.
 */
export function cfrCompressionRatioAtCounter(counter: number): number {
  const h = CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER - CFR_COUNTER_STEP * (counter - CFR_BASIC_COUNTER);
  return 1 + STROKE / h;
}

/** Inverse of cfrCompressionRatioAtCounter (closed form). */
export function cfrCounterAtCompressionRatio(compressionRatio: number): number {
  return CFR_BASIC_COUNTER + (CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER - STROKE / (compressionRatio - 1)) / CFR_COUNTER_STEP;
}

/**
 * Oil-measured compression ratio vs cylinder-height digital counter of the ANL
 * CFR engine, [Choi18] Fig. 4 (R² = 1.000), valid for counter 400–1400:
 *   CR = 1.126e-8 c³ − 2.126e-5 c² + 1.694e-2 c + 1.024.
 * The published fit as is. Its local slope dV_c/dc varies by 2× over the range, which a rigid
 * cylinder raise (0.0007 in per digit, [D2699 A2.2.1.1]) cannot produce, so prefer
 * cfrCompressionRatioAtCounter (physically consistent, extrapolates) for modelling.
 */
export function cfrCompressionRatioFromCounter(counter: number): number {
  const c = counter;
  return ((1.126e-8 * c - 2.126e-5) * c + 1.694e-2) * c + 1.024;
}

/**
 * Clearance-volume offset between the compression ratios quoted by the 2020s ANL CFR campaigns
 * (Hoth & Kolodziej 2021/2025, SON 2023, Kalvakala et al.) and cfrCompressionRatioAtCounter, m³:
 * the 2025 paper (OSTI 2561394, test/fixtures/cfr_hoth_knock_metrics.json) gives CR 7.26 for the
 * standard RON-95 cylinder height, i.e. counter 805 → 6.819 here, so
 *   ΔV_c = V_d/(6.819 − 1) − V_d/(7.26 − 1) = 7.38 cm³
 * (validation round 2, measured-data validator: with it the Kalvakala motored HCCI compression
 * pressures match within ±2 % instead of +8…+11 %, and SON 2023 PRF100 at 1.013 bar gives 7.30).
 */
export const CFR_ANL2020S_CLEARANCE_OFFSET = (() => {
  const vd = (Math.PI / 4) * BORE * BORE * STROKE;
  return vd / (cfrCompressionRatioAtCounter(805) - 1) - vd / (7.26 - 1);
})();

/** Compression ratio on this model's scale for a CR quoted by the 2020s ANL campaigns (see CFR_ANL2020S_CLEARANCE_OFFSET). */
export function cfrCompressionRatioFromAnl2020s(crQuoted: number): number {
  const vd = (Math.PI / 4) * BORE * BORE * STROKE;
  return 1 + vd / (vd / (crQuoted - 1) + CFR_ANL2020S_CLEARANCE_OFFSET);
}

/** Inverse of cfrCompressionRatioFromCounter (Newton, monotonic on 400–1400). */
export function cfrCounterFromCompressionRatio(compressionRatio: number): number {
  let c = 900;
  for (let i = 0; i < 50; i++) {
    const f = cfrCompressionRatioFromCounter(c) - compressionRatio;
    const df = (3 * 1.126e-8 * c - 2 * 2.126e-5) * c + 1.694e-2;
    const dc = f / df;
    c -= dc;
    if (Math.abs(dc) < 1e-10) break;
  }
  return c;
}

// ─────────────────────────────────────────────────────────────────────────────
// Standard operating conditions
// ─────────────────────────────────────────────────────────────────────────────

/**
 * MON spark advance vs compression ratio, crank degrees BTDC.
 * [D2700] §10.3.16 and Table 3: the ignition timer is linked to the cylinder height; 26° btdc at
 * digital counter 264 (dial 0.825 in) falling linearly to 14° btdc at counter 1145 (dial 0.200 in),
 * 1° per ≈ 73.4 digits. The counter is converted with the rigid-raise model
 * (cfrCompressionRatioAtCounter): 26° at CR 4.91, 14° at CR 9.41. Outside the tabulated range the
 * advance is held at the end values (the table stops there; [OSTI1583125] "14–26").
 */
export function cfrMonSparkAdvanceDeg(compressionRatio: number): number {
  const c = cfrCounterAtCompressionRatio(compressionRatio);
  const adv = 26 + ((14 - 26) * (c - 264)) / (1145 - 264);
  return Math.min(26, Math.max(14, adv));
}

/** Crankcase oil temperature, K: 57 ± 8 °C [D2699 §10.3.4] ([CFRE]). */
export const CFR_OIL_TEMPERATURE = 57 + 273.15;

/** Coolant (jacket) temperature, K: 100 ± 1.5 °C, boiling thermal-syphon jacket [D2699 §7.1, §10.3.5] ([W850] "about 208 °F"). */
export const CFR_COOLANT_TEMPERATURE = 100 + 273.15;

/**
 * Crankcase gauge pressure, Pa: "less than zero (a vacuum) and is typically from 25 mm to
 * 150 mm (1 in. to 6 in.) of water less than atmospheric pressure" [D2699 §10.3.10]; mid-point
 * 87.5 mm H₂O (× 9.80665 Pa/mm H₂O). Use p_crankcase = p_ambient + this in the gas torque.
 */
export const CFR_CRANKCASE_GAUGE_PRESSURE = -87.5 * 9.80665;

/**
 * Intake-air humidity ratio, kg water / kg dry air: mid-point of the ASTM 0.00356–0.00712 kg/kg
 * (25–50 grains/lb dry air) band [D2699 §10.3.7] (37.5 gr/lb = 5.357e-3).
 */
export const CFR_INTAKE_HUMIDITY_RATIO = (37.5 * 64.79891e-6) / 0.45359237;

/**
 * Fuel–air equivalence ratio used for the standard presets. The rating is done at the fuel–air
 * ratio for maximum knock intensity (carburettor fuel level) [D2699 §10.3.19 / D2700 §10.3.19];
 * for iso-octane the ANL engine found λ = 0.89 (φ = 1.124) [Pal18] Table 4, [Choi18]. We use φ = 1.10
 * for PRF 90 (UNVERIFIED: the peak-knock φ of PRF 90 was not found; knock peaks are broad, ±0.05).
 * Measured peak-knock λ for PRF 93–100 (test/fixtures): 0.8855 PRF98 [KW17], 0.9047 PRF98 [Choi18
 * Table 6], 0.89 PRF100 [Pal18], 0.88 PRF93/95/97 (Hoth & Kolodziej 2021 Table 7) — φ 1.105–1.136;
 * 1.10 is kept (inside the ±0.05 plateau; validation round 1 knock finding 9, not changed).
 */
export const CFR_PEAK_KNOCK_PHI = 1.1;

/**
 * Guide-table digital counter readings for standard knock intensity of PRF 90 at 101.0 kPa (9/16 in
 * venturi, uncompensated): RON 726 [D2699 Table A4.1], MON 749 [D2700 Table A4.1].
 */
export const CFR_GUIDE_COUNTER_PRF90 = { RON: 726, MON: 749 } as const;

/** Research-method standard-knock CR for PRF 90: counter 726 → CR 6.43 (polynomial: 6.43). */
const RON_CR_PRF90 = cfrCompressionRatioAtCounter(CFR_GUIDE_COUNTER_PRF90.RON);
/** Motor-method standard-knock CR for PRF 90: counter 749 → CR 6.54 (polynomial: 6.52). */
const MON_CR_PRF90 = cfrCompressionRatioAtCounter(CFR_GUIDE_COUNTER_PRF90.MON);

/**
 * ASTM D2699 Research method, PRF 90.
 *  - 600 ± 6 rpm, 13° btdc "regardless of cylinder height", intake air 52 ± 1 °C at standard barometric
 *    pressure [D2699 Table 1, §10.2.1, §10.3.6, §10.3.15] ([KW17], [Pal18], [CFRE]).
 *  - relativeHumidity: 37.5 gr/lb dry air at 52 °C, 1 atm → RH = 0.0636 (ASHRAE psychrometrics with
 *    Cantera water psat, tools/reference/mechanics_cfr_conditions.py).
 *  - intakeMixtureTemperature: the RON mixture temperature is not controlled [KW17]; 302.3 K is the
 *    adiabatic mixing temperature of the 52 °C humid air with PRF90 liquid at 25 °C evaporating completely
 *    at φ = 1.1 (Cantera + NIST latent heats, same script) = 29 °C, inside the 15–45 °C range of intake-port
 *    temperatures measured during RON tests on the ANL engine ([Choi18] Fig. 11). UNVERIFIED as a model of
 *    the real port temperature (partial evaporation, manifold heat transfer).
 *  - compressionRatio 6.43: guide-table counter 726 for 90 O.N. [D2699 Table A4.1] converted with the
 *    ANL clearance-volume calibration (cfrCompressionRatioAtCounter). The knock-limited CR is an OUTPUT
 *    of validation step 6 (DESIGN.md); this is the value it must reproduce.
 *  - The belt-coupled power-absorption motor holds the speed: speedMode 'fixed' [D2699 §7.1].
 *  - No throttle: throttle = 1.
 */
export const CFR_RON_CONDITIONS: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 600,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: P_ATM,
  ambientTemperature: 52 + 273.15,
  relativeHumidity: 0.0636,
  intakeMixtureTemperature: 302.3,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: CFR_PEAK_KNOCK_PHI,
  sparkAdvanceDeg: 13,
  dwellTime: INDUCTIVE_IGNITION.dwellTime,
  compressionRatio: RON_CR_PRF90,
  egrFraction: 0,
  coolantTemperature: CFR_COOLANT_TEMPERATURE,
};

/**
 * ASTM D2700 Motor method, PRF 90.
 *  - 900 rpm, intake air 38 ± 2.8 °C, mixture 149 ± 1 °C (heater downstream of the carburettor)
 *    [D2700 §10.3.6–10.3.7] ([KW17] Table 1, [CFRE]).
 *  - relativeHumidity: 37.5 gr/lb at 38 °C, 1 atm → RH = 0.1308 (same script).
 *  - compressionRatio 6.54: guide-table counter 749 for 90 O.N. [D2700 Table A4.1] (NB: higher than the
 *    RON value — each method's guide table is set at its own standard knock intensity).
 *  - sparkAdvanceDeg = cfrMonSparkAdvanceDeg(CR) = 19.4° [D2700 Table 3]; it must follow CR changes.
 */
export const CFR_MON_CONDITIONS: OperatingPoint = {
  speedMode: 'fixed',
  rpm: 900,
  loadTorque: 0,
  throttle: 1,
  ambientPressure: P_ATM,
  ambientTemperature: 38 + 273.15,
  relativeHumidity: 0.1308,
  intakeMixtureTemperature: 149 + 273.15,
  fuel: { kind: 'PRF', octaneNumber: 90 },
  equivalenceRatio: CFR_PEAK_KNOCK_PHI,
  sparkAdvanceDeg: cfrMonSparkAdvanceDeg(MON_CR_PRF90),
  dwellTime: INDUCTIVE_IGNITION.dwellTime,
  compressionRatio: MON_CR_PRF90,
  egrFraction: 0,
  coolantTemperature: CFR_COOLANT_TEMPERATURE,
};

/**
 * Standard knock-meter pickup: CFR D-1 detonation pickup + 501-C detonation meter
 * ([KW17] Table 3, [Choi18]). The pickup screws into a port in the cylinder head
 * ([Choi18]: oil was poured "through the knockmeter pickup port in the head").
 * Location: on the diameter opposite the spark plug ([Pal18] Fig. 11: knockmeter cavity at the
 * chamber periphery opposite the plug); cavity centre ≈ 3 mm inside the bore wall, cavity
 * radius ≈ 8 mm (our pixel measurement of [Pal18] Fig. 11, ±1 mm). UNVERIFIED: pickup axis
 * vertical (into the head).
 */
export const CFR_KNOCK_PICKUP = {
  type: 'CFR D-1 detonation pickup (magnetostrictive) + 501-C detonation meter',
  /** Port centre on the head face, cylinder frame (x, y, z), m. */
  position: [0, 0, -(R_BORE - 0.003)] as [number, number, number],
  /** Pickup axis, cylinder frame (pointing away from the chamber). */
  axis: [0, 1, 0] as [number, number, number],
};

/**
 * Nominal ASTM valve events, crank degrees (firing-TDC convention): IVO 10.0 ± 2.5° atdc (−350),
 * IVC 34° abdc (−146), EVO 40° bbdc (140), EVC 15.0 ± 2.5° atdc (375 ≡ −345) [D2699 §10.2.3],
 * [Kal22] Table 1. These are CAM events (on the quieting ramps, [D2699 A2.1]); the valve itself
 * moves only once the 0.008 in clearance is taken up — see CFR_F1.intakeValve / exhaustValve.
 */
export const CFR_VALVE_TIMING_ASTM = {
  intakeOpenDeg: -350,
  intakeCloseDeg: -146,
  exhaustOpenDeg: 140,
  exhaustCloseDeg: -345,
} as const;

/**
 * Measured valve timing of the ANL CFR engine at the 0.008 in hot lash [Choi18]
 * Table 2, converted to the firing-TDC convention (crank deg): IVO 376 → −344,
 * IVC −152, EVO 141, EVC 373 → −347. (At this lash the valve events do not overlap:
 * EVC is 3° before IVO.) These are the values used in CFR_F1.
 */
export const CFR_VALVE_TIMING_CHOI2018 = {
  intakeOpenDeg: -344,
  intakeCloseDeg: -152,
  exhaustOpenDeg: 141,
  exhaustCloseDeg: -347,
  lash: CFR_VALVE_LASH,
} as const;

/**
 * Measured valve lift (m) of the ANL CFR engine at crank angle θ (deg, firing-TDC convention, any
 * value — wrapped into the 720° cycle) for a valve clearance `lash` (m, default the hot running
 * clearance 0.008 in): max(0, L₀(θ) − lash), with L₀ the zero-lash profile of [Choi18] Fig. 3
 * (cfr-valve-lift.ts, 1° table, linear interpolation). Allocation-free.
 */
export function cfrValveLiftChoi2018(valve: 'intake' | 'exhaust', thetaDeg: number, lash: number = CFR_VALVE_LASH): number {
  const table = valve === 'intake' ? CHOI2018_INTAKE_LIFT_ZERO_LASH : CHOI2018_EXHAUST_LIFT_ZERO_LASH;
  const start = valve === 'intake' ? CHOI2018_INTAKE_LIFT_START_DEG : CHOI2018_EXHAUST_LIFT_START_DEG;
  // wrap into [start, start + 720)
  let u = (thetaDeg - start) % 720;
  if (u < 0) u += 720;
  const n = table.length - 1;
  if (!(u < n)) return 0; // outside the event (tables end with an exact zero)
  const i = Math.floor(u);
  const f = u - i;
  const l0 = table[i] + f * (table[i + 1] - table[i]);
  return l0 > lash ? l0 - lash : 0;
}

/** Render-only dimensions of the cylinder assembly, m ([Choi18] Table 3; UNVERIFIED where noted). */
export const CFR_RENDER_GEOMETRY = {
  /** Piston overall height [Choi18] Table 3. */
  pistonHeight: 0.1207,
  /** Cylinder (liner) length [Choi18] Table 3. */
  cylinderLength: 0.14,
  /** Cylinder wall thickness [Choi18] Table 3. */
  cylinderWallThickness: 0.0066,
  /** Head deck thickness [Choi18] Table 3. */
  headDeckThickness: 0.0129,
  /** Number of piston rings: 1 chrome/ferrous + 3 ferrous compression + 1 slotted oil ring [D2699] Table 1, [KW17] Table 2. */
  pistonRings: 5,
};

// ─────────────────────────────────────────────────────────────────────────────
// Friction (PNH model inputs)
// ─────────────────────────────────────────────────────────────────────────────

/** Crankcase oil grade: "An SAE 30 viscosity grade oil … 9.3 to 12.5 cSt at 100 °C, viscosity index ≥ 85" [D2699 §8.2]. */
export const CFR_OIL_GRADE = 'SAE30' as const;

/**
 * Inputs of the Patton–Nitschke–Heywood friction model for the CFR
 * (mechanics/friction.ts). The CFR is known for high friction: cast-iron piston
 * with five rings [D2699][W850], sleeve bearings "of double the normal dimensions
 * for an engine of this size" [W850], monograde oil at a low sump temperature.
 *  - ringTensionFactor 5/3: PNH's ring term is calibrated on 3-ring pistons; the CFR has 5 [D2699 Table 1]
 *    (UNVERIFIED assumption: ring-tension friction ∝ number of rings).
 *  - viscosityRatio: SAE 30 [D2699 §8.2] (Vogel constants, Sandoval 2002 App. A.2) at the 57 °C sump
 *    temperature [D2699 §10.3.4] = 46.3 cSt vs the 10.6 cSt PNH reference → ν/ν₀ = 4.37; friction.ts applies
 *    √(ν/ν₀) = 2.09 to the hydrodynamic terms (Sandoval 2002 eq. 3).
 *  - auxiliaryFactor 0.5: UNVERIFIED — the CFR drives an oil pump (172–207 kPa [D2699 §10.3.3]) but no
 *    alternator (the belt-coupled motor starts it) and no water pump (thermal-syphon jacket [D2699 §7.1]).
 *  - Bearing dimensions, cam bearings, flat followers: UNVERIFIED estimates ("double the normal"
 *    sleeve bearings → L/D ≈ 0.8); OHV valvetrain (pushrods and an open rocker assembly [D2699 Table 1]).
 */
export const CFR_FRICTION: PnhFrictionInputs = {
  bore: BORE,
  stroke: STROKE,
  cylinders: 1,
  compressionRatio: CR_NOMINAL,
  mainBearings: { count: 2, diameter: 2.25 * IN, length: 1.75 * IN }, // UNVERIFIED
  rodBearings: { count: 1, diameter: 2.0 * IN, length: 1.5 * IN }, // UNVERIFIED
  camBearings: 2, // UNVERIFIED
  valves: 2,
  maxValveLift: CFR_VALVE_LIFT_ZERO_LASH - CFR_VALVE_LASH, // [D2699 §10.2.4, §10.3.2.1]
  valvetrain: 'OHV',
  follower: 'flat', // UNVERIFIED: flat tappets ("valve lifter" [D2699 A2.1.2])
  ringTensionFactor: 5 / 3,
  auxiliaryFactor: 0.5, // UNVERIFIED
  viscosityRatio: oilViscosityCst(CFR_OIL_GRADE, CFR_OIL_TEMPERATURE) / PNH_REFERENCE_VISCOSITY_CST,
  boundaryRpm: 1000, // PNH 1989
};
