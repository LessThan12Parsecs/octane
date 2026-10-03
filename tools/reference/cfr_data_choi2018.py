"""Choi, Kolodziej, Wallner & Hoth 2018 (SAE 2018-01-0848, OSTI 1501884): ANL CFR F1/F2 in-cylinder, intake-port
and exhaust-port pressure traces (cycle-averaged, full 720 deg cycle) and oil-measured CR vs digital counter.

Source: tools/reference/cfr_sources/osti_1501884_choi2018.pdf (accepted manuscript, https://www.osti.gov/servlets/purl/1501884).

Extraction: VECTOR-PATH DIGITISATION. The Excel charts are embedded as PDF drawing operators; every data point of a
series is a Bezier/line end point. Page coordinates are mapped to data with a least-squares linear fit through the
chart's own tick-mark segments (residuals reported). This is exact up to the 0.001 pt coordinate rounding of the PDF:
  Fig. 2 (p. 3): 0.004 deg, 0.0005 bar (cylinder) / 0.00001 bar (ports).
  Fig. 9 (p. 5): the four series are stored for the FULL cycle (-360..360 CAD) although the published axes clip
  them to -30..90 CAD (the clip path hides the rest); we extract the stored data. 0.0006 deg, 0.0004 bar.
  Fig. 4 (p. 3, CR vs counter) is NOT extracted: its ~800 marker paths are a dense rendering that follows the
  published cubic within 0.06-0.12 CR (not the 50-count oil measurements), and the cubic is already in cfr.ts.

Crank angle: the source uses CAD ATDC with 0 = firing TDC (same as the project); IVO is printed as 376 (i.e. -344).
Pressure: bar ABSOLUTE in the source (intake-stroke cylinder pressure ~0.9-1.05 bar, port pressures 0.9-1.3 bar;
the ANL Kulite port transducers are absolute, 2.0/3.5 bara, KW17 Table 4) -> Pa here, no offset applied.

Derived checks (computed here from the digitised traces, NOT published): gross IMEP over -180..180 and net IMEP
over the full cycle, using exact slider-crank dV (independent of the unknown clearance volume), bore 82.55 mm,
stroke 114.3 mm, rod 254 mm. The published gIMEP values (Fig. 9 insets) are listed alongside for validation.

Writes test/fixtures/cfr_choi2018_traces.json.
"""
from __future__ import annotations

import numpy as np

from cfr_data_common import (
    BAR, BORE, CONROD, STROKE, T0, Axis, conditions, drawing_points, page_drawings, page_words, prf, source_pdf,
    write_fixture, cr_rigid_raise, counter_rigid_raise,
)

PDF = source_pdf("osti_1501884_choi2018.pdf")
CITATION = ("Choi S., Kolodziej C., Wallner T., Hoth A., 'Development and validation of a three pressure analysis (TPA) "
            "GT-Power model of the CFR F1/F2 engine for estimating cylinder conditions', SAE Technical Paper 2018-01-0848 "
            "(2018), accepted manuscript OSTI 1501884")
URL = "https://www.osti.gov/servlets/purl/1501884"


def dvdtheta(theta_deg: np.ndarray) -> np.ndarray:
    """dV/dtheta (m^3/rad) of the CFR slider crank (no offset), independent of the clearance volume."""
    a = STROKE / 2.0
    l = CONROD
    th = np.radians(theta_deg)
    s = np.sin(th)
    return (np.pi * BORE ** 2 / 4.0) * a * s * (1.0 + a * np.cos(th) / np.sqrt(l * l - (a * s) ** 2))


def imep(theta: np.ndarray, p: np.ndarray, lo: float, hi: float) -> float:
    sel = (theta >= lo) & (theta <= hi)
    th, pp = theta[sel], p[sel]
    w = np.trapezoid(pp * dvdtheta(th), np.radians(th))
    return float(w / (np.pi * BORE ** 2 / 4.0 * STROKE))


def series_xy(dr: dict, xax: Axis, yax: Axis) -> tuple[np.ndarray, np.ndarray]:
    pts = drawing_points(dr)
    return xax(pts[:, 0]), yax(pts[:, 1])


def to_grid(x: np.ndarray, y: np.ndarray, step: float = 0.1) -> tuple[np.ndarray, np.ndarray]:
    """The stored series are already on a uniform ~0.1 deg grid; snap to exact multiples of `step` by linear
    interpolation (removes the <0.01 deg calibration jitter)."""
    order = np.argsort(x)
    x, y = x[order], y[order]
    lo = np.ceil(x[0] / step - 1e-6) * step
    hi = np.floor(x[-1] / step + 1e-6) * step
    g = np.round(np.arange(lo, hi + 0.5 * step, step), 4)
    return g, np.interp(g, x, y)


def peak(theta, p):
    i = int(np.argmax(p))
    return float(theta[i]), float(p[i])


def main() -> None:
    # ───────── Fig. 2: PRF98 standard RON test, cylinder / intake / exhaust pressure (p. 3) ─────────
    drs = page_drawings(PDF, 3)
    xt = [it[1].x for it in drs[71]["items"]]  # x ticks
    ylt = [it[1].y for it in drs[69]["items"]]  # left y ticks (cylinder, bar)
    yrt = [it[1].y for it in drs[67]["items"]]  # right y ticks (ports, bar)
    xax = Axis.fit(xt, [-360, -180, 0, 180, 360])
    yl = Axis.fit(ylt, [0, 10, 20, 30, 40, 50])
    yr = Axis.fit(yrt, [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
    # legend (drawings 76-78): solid = P-Cylinder, dotted [0 1.47] = P-Intake, dashed [1.47 1.47] = P-Exhaust
    assert drs[72]["dashes"] in ("[] 0", None) and "0 1.47" in drs[73]["dashes"] and "1.47" in drs[74]["dashes"]
    thc, pc = to_grid(*series_xy(drs[72], xax, yl))
    thi, pi_ = to_grid(*series_xy(drs[73], xax, yr), step=0.2)
    the, pe = to_grid(*series_xy(drs[74], xax, yr), step=0.2)
    prf98 = {
        "key": "choi2018_fig2_prf98_ron",
        "kind": "pressureTraceCycleAveraged",
        "figure": "Fig. 2", "page": 3,
        "caption": "Cycle-averaged instantaneous cylinder, intake, and exhaust pressure traces of RON test of PRF98 "
                   "(300-cycle ensemble average, used as TPA model input)",
        "conditions": conditions(
            engine="ANL CFR F1/F2", method="RON", fuel=prf(98), rpm=600, sparkAdvanceDegBTDC=13,
            compressionRatio=None,
            compressionRatioSource="not stated; guide-table counter for RON 98.0 is 867 at 29.92 inHg "
                                   f"(cfr_astm_guide_tables.json) -> CR {cr_rigid_raise(867):.3f} (rigid raise) before "
                                   "barometric compensation (the PRF100 CR 7.55 quoted for the same campaign implies ~+9 counts, i.e. a "
                                   "~29.6 inHg barometer -> counter ~876, CR ~7.22; UNVERIFIED inference)",
            digitalCounter=None, **{"lambda": 0.9047},
            intakeAirTemperatureK=None, intakeMixtureTemperatureK=None, coolantTemperatureK=373.0, oilTemperatureK=323.0,
            knockIntensityKU=50,
            notes="Standard RON rating conditions compensated for the day's barometer (Choi 2018 p. 2). lambda 0.9047 = "
                  "default measured value of the PRF98 standard RON DoE case (Table 5/6); gIMEP 7.9957 bar and fuel rate "
                  "0.7370 kg/h in the same case (Table 6). Coolant 373 K / oil 323 K are the TPA-model defaults (Table 5), "
                  "not measurements. Cylinder: AVL spark-plug transducer (GU13Z-24 per later ANL papers)."),
        "reported": {"gIMEPbar_Table6_default": 7.9957, "fuelRateKgPerH_Table6_default": 0.7370,
                     "lambda_Table6_default": 0.9047},
        "extraction": {
            "method": "vector (PDF drawing operators: cylinder = path #72, 7199 Bezier segments; intake = #73 and exhaust = "
                      "#74, 3599 segments each)",
            "axisCalibration": {"x": "5 tick segments -360..360", "yLeft": "6 ticks 0..50 bar", "yRight": "7 ticks 0.8..1.4 bar",
                                "maxResidualDeg": xax.resid, "maxResidualBarLeft": yl.resid, "maxResidualBarRight": yr.resid},
            "resolution": {"thetaDeg": 0.1, "pressureBarCylinder": 0.001 * yl.per_unit(), "pressureBarPorts": 0.001 * yr.per_unit()},
            "uncertainty": "digitisation < 0.01 deg / 0.001 bar; measurement: spark-plug transducer reads ~1-2% low in "
                           "combustion/expansion vs flush-mount (Choi Fig. 9), 300-cycle average of knocking cycles "
                           "smooths the ringing.",
            "resampled": "linear interpolation onto exact 0.1 deg (cylinder) / 0.2 deg (ports) grids",
        },
        "data": {
            "cylinder": {"thetaDeg": thc, "pPa": pc * BAR},
            "intakePort": {"thetaDeg": thi, "pPa": pi_ * BAR},
            "exhaustPort": {"thetaDeg": the, "pPa": pe * BAR},
        },
    }
    prf98["derivedFromTrace"] = {
        "peakPressurePa": peak(thc, pc * BAR)[1], "peakPressureThetaDeg": peak(thc, pc * BAR)[0],
        "gIMEPPa_minus180to180": imep(thc, pc * BAR, -180, 180),
        "nIMEPPa_fullCycle": imep(thc, pc * BAR, -360, 359.9),
        "pAtMinus180Pa": float(np.interp(-180, thc, pc * BAR)),
        "meanIntakePortPa": float(np.mean(pi_ * BAR)), "meanExhaustPortPa": float(np.mean(pe * BAR)),
        "note": "computed from the digitised trace with exact slider-crank dV/dtheta (B 82.55, S 114.3, L 254 mm)",
    }

    # ───────── Fig. 9: PRF100 spark plug vs flush-mount transducer, ST 5 and -13 CAD ATDC (p. 5) ─────────
    drs5 = page_drawings(PDF, 5)
    xt9 = [it[1].x for it in drs5[247]["items"]]
    yt9 = [it[1].y for it in drs5[245]["items"]]
    x9 = Axis.fit(xt9, [-30, 0, 30, 60, 90])
    y9 = Axis.fit(yt9, [0, 10, 20, 30, 40, 50])
    legend = {248: ("ST5a", "sparkPlug"), 249: ("ST5a", "flushMount"), 250: ("ST-13a", "sparkPlug"), 251: ("ST-13a", "flushMount")}
    # verify the legend styles (drawings 253-256: grey solid, grey dotted, black solid, black dotted)
    for tr, lg in zip((248, 249, 250, 251), (253, 254, 255, 256)):
        assert drs5[tr]["color"] == drs5[lg]["color"] and drs5[tr]["dashes"] == drs5[lg]["dashes"], (tr, lg)
    reported = {("ST-13a", "sparkPlug"): 7.938, ("ST-13a", "flushMount"): 8.137,
                ("ST5a", "sparkPlug"): 7.547, ("ST5a", "flushMount"): 7.892}
    fig9 = []
    for idx, (st, sensor) in legend.items():
        th, p = to_grid(*series_xy(drs5[idx], x9, y9))
        th = np.round(th, 4)
        sel = (th >= -360) & (th < 360)
        th, p = th[sel], p[sel]
        adv = 13 if st == "ST-13a" else -5
        entry = {
            "key": f"choi2018_fig9_prf100_{'st_m13' if adv == 13 else 'st_p5'}_{sensor}",
            "kind": "pressureTraceCycleAveraged",
            "figure": "Fig. 9", "page": 5,
            "caption": f"PRF100 (iso-octane), spark timing {'-13' if adv == 13 else '+5'} CAD ATDC, "
                       f"{'AVL spark-plug' if sensor == 'sparkPlug' else 'flush-mounted'} cylinder pressure transducer",
            "conditions": conditions(
                engine="ANL CFR F1/F2", method="RON-like spark sweep" if adv != 13 else "RON (standard knock)",
                fuel=prf(100), compressionRatio=7.55, compressionRatioSource="stated (Choi 2018 p. 2; Pal 2018 p. 2); oil calibration",
                rpm=600, sparkAdvanceDegBTDC=adv, **{"lambda": 0.89},
                knockIntensityKU=50 if adv == 13 else 0,
                notes="RON test conditions compensated for the actual barometer (IAT and CR adjusted per ASTM, values not "
                      "given); iso-octane at peak-knock lambda 0.89; ST -13 = standard knock (~50 KU), ST +5 = non-knocking. "
                      "Same data set as Pal et al. 2018 Figs. 1-2, 6 (cfr_pal2018_*.json). Knockmeter: D-1 pickup + 501-C."),
            "reported": {"gIMEPbar": reported[(st, sensor)],
                         "gIMEPsensorDifference": "+-1.24% at ST -13 (knocking), +-2.23% at ST 5 (Fig. 9 insets)"},
            "extraction": {
                "method": "vector (PDF path #%d, 7199 Bezier segments spanning the full cycle; the published axes clip "
                          "the view to -30..90 CAD)" % idx,
                "axisCalibration": {"maxResidualDeg": x9.resid, "maxResidualBar": y9.resid},
                "resolution": {"thetaDeg": 0.1, "pressureBar": 0.001 * y9.per_unit()},
                "uncertainty": "digitisation < 0.01 deg / 0.001 bar",
            },
            "data": {"thetaDeg": th, "pPa": p * BAR},
        }
        entry["derivedFromTrace"] = {
            "peakPressurePa": peak(th, p * BAR)[1], "peakPressureThetaDeg": peak(th, p * BAR)[0],
            "gIMEPPa_minus180to180": imep(th, p * BAR, -180, 180),
            "nIMEPPa_fullCycle": imep(th, p * BAR, -360, 359.9),
            "pAtMinus180Pa": float(np.interp(-180, th, p * BAR)),
            "pAtTDCPa": float(np.interp(0, th, p * BAR)),
            "gIMEPrelativeErrorVsReported": imep(th, p * BAR, -180, 180) / (reported[(st, sensor)] * BAR) - 1.0,
        }
        fig9.append(entry)

    # ───────── tables / text ─────────
    tables = {
        "key": "choi2018_tables",
        "kind": "reportedScalars",
        "table2_valveTiming_CADATDC_at_0p008in_lash": {"IVO": 376, "IVC": -152, "EVO": 141, "EVC": 373,
                                                        "note": "376/373 = -344/-347 in [-360,360)"},
        "table2_valveReferenceDiameterIn": {"intake": 1.346, "exhaust": 1.356},
        "table3_structureMm": {"headDeck": 12.9, "pistonTopDeck": 28.6, "pistonHeight": 120.7, "skirt": 5.0, "ring": 5.0,
                               "cylinderWall": 6.6, "cylinderLength": 140, "headToWaterJacketBottom": 121},
        "table4_coolingHTC": {"headCoolant": "5*rpm^0.8 W/m2K", "pistonOil": "1.5*rpm^0.8", "cylinderOil": "1.5*rpm^0.8",
                              "cylinderCoolant": "5*rpm^0.8"},
        "heatTransferCalibration": {"intakePortMultiplier": 4, "chamberConvectionMultiplierRange": [1, 3],
                                    "note": "TPA-model tuning to match VE (+-1%) and gIMEP; knock needs a local multiplier "
                                            "increase during ringing (Fig. 6)"},
        "validationRange": {"fuels": "PRF60-100, PRF-ethanol blends", "lambda": [0.885, 1.052],
                            "intakePortTemperatureC": [17, 40], "compressionRatio": [5.57, 7.34],
                            "sparkTimingCADATDC": [-13, 5]},
        "table6_PRF98_standardRON_default": {"lambda": 0.9047, "gIMEPbar": 7.9957, "fuelRateKgPerH": 0.7370,
                                             "inportMultiplier": 3.501, "combustionMultiplier": 1.921},
        "fig4Note": "CR vs counter: only the published cubic is usable (cfr.ts cfrCompressionRatioFromCounter); the plotted "
                    "markers are a dense curve rendering, not the individual oil measurements (checked, see generator docstring).",
        "tpaRegressions": {
            "TIVC_K": "0.5385 * T_port[degC] + 394.34 (Fig. 11, all fuels/lambdas, T_port 15-45 C)",
            "Tunburned_at_minus20_K": "0.6567 * T_port[degC] + 678.42 (Fig. 12)",
            "note": "TPA-model (GT-Power) outputs, not direct measurements; uncertainty +-5 C per +-1% lambda/fuel-rate error",
        },
    }

    out = {
        "description": "Choi et al. 2018 (ANL CFR): full-cycle cycle-averaged cylinder/intake/exhaust pressure (PRF98 RON), "
                       "PRF100 spark-timing traces with two transducers. SI units, Pa absolute, "
                       "crank angle deg with 0 = firing TDC.",
        "generator": "tools/reference/cfr_data_choi2018.py",
        "source": {"key": "Choi2018", "citation": CITATION, "url": URL, "osti": "1501884",
                   "localPdf": "tools/reference/cfr_sources/osti_1501884_choi2018.pdf"},
        "conventions": {"thetaDeg": "crank angle, deg, 0 = firing TDC, [-360, 360); source convention identical (CAD ATDC)",
                        "pPa": "absolute pressure, Pa (source bar absolute; no pegging change)"},
        "datasets": [prf98] + fig9 + [tables],
    }
    write_fixture("cfr_choi2018_traces.json", out)
    for e in [prf98] + fig9:
        d = e["derivedFromTrace"]
        print(e["key"], "pmax %.2f bar @ %.1f" % (d["peakPressurePa"] / BAR, d["peakPressureThetaDeg"]),
              "gIMEP %.3f bar" % (d["gIMEPPa_minus180to180"] / BAR), "nIMEP %.3f" % (d["nIMEPPa_fullCycle"] / BAR),
              "p(-180) %.3f" % (d["pAtMinus180Pa"] / BAR), "rep", e.get("reported"))


if __name__ == "__main__":
    main()
