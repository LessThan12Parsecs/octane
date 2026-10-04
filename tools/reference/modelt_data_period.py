"""Ford Model T engine validation data (period Ford service literature, period engineering tests, restorer
measurements), transcribed from the sources below into test/fixtures/modelt_period_data.json.

Nothing here is digitised from figures: every number is a transcription of a printed value (table, text or
service specification) or a closed-form relation printed by its author. The adversarial cross-checks behind
each value (independent second sources, corrections, disputes) are summarised per dataset in `notes`, and in
test/fixtures/modelt_validation_README.md.

Sources
  [FSB]   Ford Service Bulletins (1919-), reprinted in D. R. Post (ed.), Model T Ford Service Bulletin Essentials,
          Post Motor Books 1966 (archive.org bwb_S0-DTZ-218; access-restricted, read through Internet Archive
          full-text-search snippets): Fig. 84 WOT horsepower/torque table + text, magneto table, ~400 rpm idle.
  [Tulsa] MTFC Tulsa, DynoSummary.htm / power_and_torque.htm (tildentechnologies.com/mtfctulsa): the same Ford
          table (their 'Ford 1918' curve is Ford's values x 0.8 for comparison with chassis-dyno data).
  [Upton] G. B. Upton (Cornell), Ford-engine spark-advance tests, Journal of the SAE, Aug. 1923, p. 112.
  [FM19]  Ford Manual for Owners and Operators of Ford Cars and Trucks (1919), Project Gutenberg #46206.
  [DB23]  Ford Motor Co., 1923 Ford Dealers' Data Book (MTFCA Encyclopedia transcription).
  [MA22]  Motor Age 41(25), 22 Jun 1922, Q&A 'Ford Compression Ratios' (archive.org sim_motor-age_1922-06-22_41_25).
  [Dyke24] Dyke's Automobile and Gasoline Engine Encyclopedia, 13th ed. 1924 (altitude/compression table; magneto).
  [BP]    T. Boggess & R. Patterson; R. Patterson & S. Coniff, 'The Model T Ford Ignition System & Spark Timing'
          (fordmodelt.net/downloads/Model%20T%20Ignition.pdf).
  [ECCT]  ECCT V12 manual (Model T coil tester); [Cool] cool386.com/tester/tester.html (Ford-coil oscillograms);
          [Kossor] M. Kossor, Model T coil/HCCT oscillogram articles.
  [GEM76] W. H. Payne, letter in Gas Engine Magazine (1976): magneto volts/amperes/cycles vs rpm.

Writes test/fixtures/modelt_period_data.json.
"""
from __future__ import annotations

import json
import pathlib

IN = 0.0254
LBFT = 1.3558179483314004  # N m per lbf ft
HP = 745.69987158227022  # W per (mechanical) hp
PSI = 6894.757293168361
P_ATM = 101325.0
MPH = 0.44704

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "modelt_period_data.json"


def _round(obj, sig: int = 7):
    if isinstance(obj, float):
        return float(f"{obj:.{sig}g}")
    if isinstance(obj, dict):
        return {k: _round(v, sig) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_round(v, sig) for v in obj]
    return obj


# Ford WOT table [FSB Fig. 84]: rpm, car mph (high gear), torque lb-ft, hp. 1400 rpm torque illegible in the
# snippet (~74, from the hp column: 19.66 hp x 5252 / 1400 = 73.8 lb-ft).
FORD_WOT = [
    (300, 7.5, 35, 2.0), (400, 10, 57, 4.5), (500, 12.5, 69, 6.5), (600, 15, 73, 8.5), (700, 17.5, 78, 10.40),
    (800, 20, 81, 12.33), (900, 22.5, 83, 14.20), (1000, 25, 82, 15.60), (1100, 27.5, 81, 16.66),
    (1200, 30, 79, 18.20), (1300, 32.5, 77, 19.0), (1400, 35, 74, 19.66), (1500, 37.5, 70, 20.0),
    (1600, 40, 65, 20.0), (1700, 42.5, 60, 19.40), (1800, 45, 53, 18.20), (1900, 47.5, 47, 17.0),
]

DISPLACEMENT = 4 * 3.141592653589793 / 4 * (3.75 * IN) ** 2 * (4.0 * IN)


def ford_wot() -> dict:
    rows = []
    for rpm, mph, tq, hp in FORD_WOT:
        T = tq * LBFT
        rows.append({
            "rpm": rpm, "carSpeedMph": mph, "torqueLbFt": tq, "torqueNm": T, "horsepower": hp, "powerW": hp * HP,
            "bmepPa": 4 * 3.141592653589793 * T / DISPLACEMENT,
            "illegibleTorque": rpm == 1400,
        })
    return {
        "key": "ford_wot_1918",
        "kind": "brake torque and power vs speed, wide-open throttle",
        "source": "[FSB] Fig. 84 table and text (first Service Bulletins April 1919; test c. 1918); [Tulsa] tabulation",
        "grade": "B",
        "conditions": {
            "throttle": 1, "speedMode": "fixed (engine dynamometer)", "measuredAt": "transmission output (high gear)",
            "spark": "not stated (best attainable with the lever; the Ford ignition cannot reach MBT at speed)",
            "ignitionSource": "not stated (magneto assumed)", "fuel": "c. 1918 commercial gasoline",
        },
        "reported": {
            "peakTorque": {"rpm": 900, "torqueLbFt": 83, "torqueNm": 83 * LBFT},
            "peakPower": {"rpm": [1500, 1600], "horsepower": 20, "powerW": 20 * HP},
            "upperBoundRatingHp": 22.5,
        },
        "data": rows,
        "notes": [
            "Ford: 'While we have obtained ratings as high as 22 1/2 horsepower, we believe the figures given below "
            "are representative of the motors in general use' and 'These figures were obtained with a wide open "
            "throttle. They represent only the maximum power that can be developed at the given speeds.'",
            "Ford rounded the torque and hp columns independently (they disagree by up to 1.5 %).",
            "The 1400 rpm torque is illegible in the snippet; 74 lb-ft is back-computed from 19.66 hp.",
            "Car-speed column = Ford's gearing, 40 rpm per mph in high (30x3-1/2 tyres, 3.636 axle).",
            "Independent checks: Cornell 1923 Froude-brake tests ~70 psi peak BMEP at WOT, AFR 12-13 [Upton]; "
            "Gunnell (2007) 'Brake hp 20 at 1600 rpm, torque 83 at 900'.",
            "Tulsa's 'Ford 1918' curve (16.1 hp, 66.3 lb-ft) is this table x 0.8 (their driveline-loss "
            "comparison with chassis dynos); do not use it as engine data.",
        ],
    }


def upton_mbt() -> dict:
    rpm = [400, 600, 800, 1000, 1200, 1400, 1500]
    return {
        "key": "upton_mbt_1923",
        "kind": "MBT spark advance vs speed at WOT",
        "source": "[Upton] J. SAE Aug. 1923: a0 = 0.108 R / (1 + 0.001 R) crank degrees BTDC (R = rpm)",
        "grade": "C",
        "conditions": {"throttle": 1, "engine": "Ford Model T (Cornell test engine)", "dataRangeRpm": [0, 1500]},
        "relation": "sparkAdvanceDegBTDC = 0.108*rpm/(1 + 0.001*rpm)",
        "data": [{"rpm": r, "mbtDegBTDC": 0.108 * r / (1 + 0.001 * r)} for r in rpm],
        "derived": {
            "peakPressureAtMbtDegATDC": "a0/3 (Upton's theory: explosion time = ignition to PEAK pressure = (4/3) a0)",
            "partThrottleExtraAdvanceDeg": [10, 35],
            "torqueLossFor20DegErrorPercent": 10,
        },
        "notes": [
            "Upton's 'explosion time' is ignition to peak pressure (not a full burn duration); his 3/4 factor is "
            "a theoretical result. His data stop at about 1500 rpm; extrapolation above is the formula only.",
        ],
    }


def compression() -> dict:
    return {
        "key": "compression_pressure",
        "kind": "cranking / gauge compression pressure",
        "source": "[FM19] Answer 19; [DB23]; [MA22]; [Dyke24] altitude table",
        "grade": "C",
        "data": [
            {"source": "Ford Manual 1919 A19; 1923 Data Book", "gaugePsi": 60, "gaugePa": 60 * PSI,
             "note": "'Compression in cylinders 60 lbs.' (Data Book also: 'compressed to 40 to 60 lbs. pressure')"},
            {"source": "Motor Age 22 Jun 1922", "gaugePsi": 55, "gaugePa": 55 * PSI,
             "note": "'The gage compression on the new model Fords is approximately 55 lbs.'"},
            {"source": "Dyke's 1924 altitude table", "gaugePsi": 64, "gaugePa": 64 * PSI,
             "note": "'average Ford engine' at sea level, stated CR 3.6 (49 psi at 5000 ft, 30 psi at 14000 ft)"},
        ],
        "notes": [
            "Restorers typically read 35-45 psi when hand-cranking (low speed, leakage); Fahnestock's 30-35 psi "
            "is inconsistent with the geometry (probably a different test condition).",
            "Physics check (verification agent): CR 3.78-3.98 with intake closing 50.8 deg ABDC and a polytropic "
            "exponent 1.2-1.3 gives 49-62 psig.",
            "Test condition for the model: motored, wide-open throttle, cranking speed ~150-200 rpm (hand crank) "
            "or starter (~ 100-150 rpm) - not stated by Ford.",
        ],
    }


def coil() -> dict:
    return {
        "key": "trembler_coil",
        "kind": "Ford/K-W trembler coil electrical behaviour",
        "source": "[BP] Ford 1916 engineering data; [ECCT]; [Cool]; [Kossor]; Ford Service Par. 1002-1009",
        "grade": "B",
        "coil": {
            "primaryTurns": 212, "secondaryTurns": 16600, "primaryResistanceOhm": 0.295,
            "secondaryResistanceOhm": 3300, "primaryInductanceOpenH": 3.3e-3, "primaryInductanceShortedH": 0.6e-3,
            "secondaryInductanceOpenH": 22.0, "secondaryInductanceShortedH": 11.3,
            "condenserF": [0.40e-6, 0.45e-6], "pointGapM": (1 / 32) * IN, "cushionClearanceM": 0.005 * IN,
        },
        "firing": [
            {"supplyV": 6, "firingTimeS": 3.5e-3, "currentAtFireA": [5.0, 5.4]},
            {"supplyV": 9, "firingTimeS": 2.5e-3, "currentAtFireA": [6.2, 6.2]},
            {"supplyV": 12, "firingTimeS": 2.0e-3, "currentAtFireA": [6.0, 7.0]},
        ],
        "other": {
            "pointsRecloseAfterFireS_6V": 1.8e-3,
            "buzzFrequencyHz_6V": 190,
            "hcctFireCurrentA": [3.0, 4.4],
            "hcctAmmeterSettingA": 1.3,
            "primaryPeakV": [300, 400],
            "secondaryOpenCircuitV": [8000, 20000],
            "storedEnergyJ_at_5_53A": 49.5e-3,
            "coilToCoilSpreadS_12V": [0.0, 0.69e-3],
        },
        "notes": [
            "Firing time = timer make to the first spark, coil from rest (oscillograms, single coils).",
            "The opening current is not a fixed threshold: 3.0-4.4 A on slow magneto/HCCT pulses vs 5-7 A on DC "
            "(armature dynamics, dI/dt dependence) [Kossor].",
            "The HCCT 1.3 A is an ammeter reading on the hand-cranked tester, not the trip current.",
        ],
    }


def magneto() -> dict:
    return {
        "key": "magneto",
        "kind": "flywheel magneto output and spark ladder",
        "source": "[FSB] magneto table and service spec; [Dyke24]; [BP]; restorer measurements; [GEM76]",
        "grade": "C",
        "cyclesPerRevolution": 8,
        "serviceMinimum": {"rpm": 400, "openCircuitVrms": 7.0},
        "openCircuit_1911car_restorer": [
            {"rpm": 407, "Vrms": 10.0}, {"rpm": 610, "Vrms": 14.5}, {"rpm": 813, "Vrms": 18.0},
            {"rpm": 1016, "Vrms": 20.9}, {"rpm": 1220, "Vrms": 24.4}, {"rpm": 1423, "Vrms": 26.5},
        ],
        "table_GEM76": [
            {"rpm": 400, "V": 9.8, "A": 7.9, "Hz": 52.8}, {"rpm": 600, "V": 14.4, "A": 8.5, "Hz": 80.0},
            {"rpm": 800, "V": 18.8, "A": 8.8, "Hz": 106.4}, {"rpm": 1000, "V": 22.8, "A": 8.9, "Hz": 133.3},
            {"rpm": 1200, "V": 26.2, "A": 9.0, "Hz": 160.0},
        ],
        "sparkLadder600rpm": {
            "source": "[BP] (Patterson, Ford roller timer, HCCT-adjusted coils, cranking at full retard then advancing)",
            "firstFireDegATDC": [26.5, 4.0, -18.5, -41.0],
            "stepDeg": 22.5,
            "note": "negative = BTDC; firing positions are quantized to the magneto half-wave peaks",
        },
        "notes": [
            "Dyke's 1924: 'the frequency is 8 cycles per revolution'; Page 1929: 'primary current varies from 6 to 24 "
            "volts ... about 8 volts when engine is running throttled down'.",
            "The GEM76 1000 rpm frequency is printed as 146.4 (should be 133.3 for 8 cycles/rev); corrected here. "
            "The amperes column is assumed to be short-circuit current (|Z| ~ 1.2-2.9 ohm -> L_s ~ 3 mH, estimate).",
        ],
    }


def timer() -> dict:
    return {
        "key": "timer_and_lever",
        "kind": "ignition timer contact and spark-lever range",
        "source": "[BP] (measured on an original Ford roller timer set with Ford's 2-1/2 in gauge); Ford Service Par. 126",
        "grade": "C",
        "contactArcDeg": 87,
        "makeAtFullRetardDegATDC": 15.5,
        "makeAtFullAdvanceDegBTDC": 64.5,
        "leverTravelDeg": 80,
        "leverNotches": 28,
        "notes": ["Earlier (1917-18) 2-5/8 in setting: full-retard make about 7-8 deg ATDC (forum-derived, lower trust)."],
    }


def vehicle() -> dict:
    return {
        "key": "vehicle_performance",
        "kind": "road performance of the car (for the 'vehicle' load model)",
        "source": "[FSB] gearing and shipping weights; road tests and period recollections (see notes)",
        "grade": "C/D",
        "rpmPerMphHigh": 40,
        "topSpeedMph": [42, 45],
        "topSpeedMs": [42 * MPH, 45 * MPH],
        "gradeability": [
            {"gear": "high", "gradePercent": 5, "speedMph": 35},
            {"gear": "high", "gradePercent": 6.7, "speedMph": [30, 32]},
            {"gear": "high", "maxGradePercent": [8, 9], "vehicleMassLb": 1875},
        ],
        "fuelEconomyMpgUS": [13, 21],
        "notes": [
            "Wikipedia gives 42 mph; Ford's table puts 1700-1800 rpm at 42.5-45 mph. 'Anything up to 45 mph' is a "
            "personal recollection (Sellen, Motor Sport 1952), not a road test.",
            "Gradeability values were computed by the research agent from the Ford curve with the vehicle model - "
            "treat them as consistency checks, not measurements.",
        ],
    }


def main() -> None:
    payload = {
        "description": "Ford Model T (1917-1927 high-head engine) period validation data: Ford's WOT dyno table, "
                       "Upton's MBT spark relation, compression pressure, trembler-coil and magneto behaviour, timer "
                       "and lever range, vehicle performance.",
        "generator": "tools/reference/modelt_data_period.py",
        "conventions": {
            "units": "SI unless the key names the unit (Lb, Psi, Mph, Hp, Vrms)",
            "crankAngle": "deg; ATDC positive after TDC, BTDC values given as positive 'BTDC' numbers where named so",
            "grades": "A exact transcription of a primary table; B primary data read through snippets or minor gaps; "
                      "C secondary/single-source or derived relations; D recollection / trend only",
        },
        "datasets": [ford_wot(), upton_mbt(), compression(), coil(), magneto(), timer(), vehicle()],
    }
    OUT.write_text(json.dumps(_round(payload), indent=1, ensure_ascii=False) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
