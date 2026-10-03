"""CFR F-1 operating-condition and geometry reference numbers.

1. Intake-air humidity. ASTM D2699/D2700 specify 25–50 grains of water per lb of
   dry air (quoted e.g. by Sinpar, "CFR test engine unit – octane rating operation
   conditions"). Converted to humidity ratio W (kg/kg dry air) and to relative
   humidity at the intake-air temperature and 1 atm using
       p_w = W p / (0.621945 + W)           (ASHRAE Fundamentals 2017 ch. 1 eq. 20/22)
   with the saturation pressure of water from Cantera's PureFluid water model
   (ct.Water(), Reynolds 1979 equation of state).
2. Choi, Kolodziej, Wallner & Hoth (2018), SAE 2018-01-0848, Fig. 4: measured
   (oil-filling) compression ratio vs cylinder-height digital counter,
       CR = 1.126e-8 c³ − 2.126e-5 c² + 1.694e-2 c + 1.024     (400 ≤ c ≤ 1400).
   We tabulate it and the implied effective clearance height h_eff = S/(CR − 1)
   (= V_c/A_p), and fit h_eff linearly in c — the cylinder is raised by a worm
   gear, so the clearance height should be linear in the counter.
   Rigid-raise model: ASTM D2699-15a A2.2.1.1 states that one digital-counter digit is
   0.0007 in of cylinder-height movement (also Table A4 footnote: dial = 1.012 − c/1410 in),
   so h_eff(c) = h_930 − 0.0007 in·(c − 930) exactly (930 = basic counter setting,
   D2699 A2.2.2). The single unknown h_930 is fitted by least squares (in h, i.e. in
   clearance VOLUME, the quantity the oil method measures) to the Choi polynomial on its
   valid range 400–1400.
4. ASTM guide tables (standard knock intensity, 9/16 in venturi, 101.0 kPa, uncompensated
   digital counter): D2699-15a Table A4.1 RON 90.0 → 726, RON 100.0 → 919;
   D2700-14 Table A4.1 MON 90.0 → 749. D2700-14 Table 3: MON spark timing vs digital
   counter, 26° btdc at 264 … 14° btdc at 1145 (13 rows, linear). Converted to CR with
   the rigid-raise model (and, for comparison, the Choi polynomial).

3. RON mixture temperature estimate. In the RON test the intake AIR is heated to
   52 °C upstream of the carburettor; the mixture temperature is "not controlled"
   (Kolodziej & Wallner 2017, Table 1). Our model feeds a fully vaporised charge, so
   the port temperature is estimated by adiabatic mixing of humid air at 52 °C with
   liquid PRF at 25 °C that evaporates completely (φ = 1.1). Vapour enthalpies from
   Cantera's nasa_gas.yaml; latent heats at 298 K from the NIST Chemistry WebBook:
   2,2,4-trimethylpentane 35.1 ± 0.2 kJ/mol (average of 7 values), n-heptane
   36.4 kJ/mol at 303 K (Van Ness et al. 1967). PRF90 = 90 % iso-octane by liquid
   volume → mole fraction with 20 °C densities 0.6919 / 0.6837 g/cm³ (CRC Handbook;
   UNVERIFIED here — only a small effect on the estimate).

Writes test/fixtures/mechanics_cfr_conditions.json.
"""
from __future__ import annotations

import json
import pathlib

import cantera as ct
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
OUT = ROOT / "test" / "fixtures" / "mechanics_cfr_conditions.json"

GRAIN = 64.79891e-6  # kg (exact, international grain)
LB = 0.45359237  # kg (exact)
P_ATM = 101325.0
STROKE = 4.5 * 0.0254


def psat(T: float) -> float:
    w = ct.Water()
    w.TQ = T, 0.0
    return float(w.P)


def humidity(grains_per_lb: float, T: float, p: float = P_ATM) -> dict:
    W = grains_per_lb * GRAIN / LB
    pw = W * p / (0.621945 + W)
    ps = psat(T)
    return dict(grainsPerLb=grains_per_lb, T=T, p=p, humidityRatio=W, vapourPressure=pw, psat=ps, relativeHumidity=pw / ps)


def choi_cr(c: float) -> float:
    return 1.126e-8 * c**3 - 2.126e-5 * c**2 + 1.694e-2 * c + 1.024


COUNTER_STEP = 0.0007 * 0.0254  # m per digit, ASTM D2699-15a A2.2.1.1
BASIC_COUNTER = 930.0  # ASTM D2699-15a A2.2.2 / D2700-14 10.3.18.1
MON_SPARK_TABLE = [(264, 26), (337, 25), (410, 24), (484, 23), (556, 22), (630, 21), (704, 20),
                   (777, 19), (851, 18), (925, 17), (998, 16), (1072, 15), (1145, 14)]  # D2700-14 Table 3


def rigid_raise_fit() -> dict:
    c = np.arange(400.0, 1401.0, 10.0)
    h = STROKE / (np.array([choi_cr(x) for x in c]) - 1.0)
    h930 = float(np.mean(h + COUNTER_STEP * (c - BASIC_COUNTER)))  # LSQ intercept with the slope fixed
    hm = h930 - COUNTER_STEP * (c - BASIC_COUNTER)
    cr_model = 1.0 + STROKE / hm
    cr_choi = 1.0 + STROKE / h

    def cr_at(cnt: float) -> float:
        return 1.0 + STROKE / (h930 - COUNTER_STEP * (cnt - BASIC_COUNTER))

    guide = {"RON90": 726.0, "RON100": 919.0, "MON90": 749.0}
    return dict(
        counterStep=COUNTER_STEP,
        basicCounter=BASIC_COUNTER,
        effectiveClearanceHeightAtBasicCounter=h930,
        maxAbsHeightResidual=float(np.max(np.abs(hm - h))),
        maxAbsCrDifference=float(np.max(np.abs(cr_model - cr_choi))),
        guideTableCounters=guide,
        guideTableCR={k: cr_at(v) for k, v in guide.items()},
        guideTableCRChoiPolynomial={k: choi_cr(v) for k, v in guide.items()},
        monSparkTable=[dict(counter=float(cnt), advanceDeg=float(a), compressionRatio=cr_at(cnt)) for cnt, a in MON_SPARK_TABLE],
    )


def ron_mixture_temperature(on: float = 90.0, phi: float = 1.1, T_air: float = 325.15, T_fuel: float = 298.15,
                            grains: float = 37.5) -> dict:
    names = ["N2", "O2", "Ar", "CO2", "H2O", "C8H18,isooctane", "C7H16,n-heptane"]
    sp = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    gas = ct.Solution(thermo="ideal-gas", species=[sp[n] for n in names])
    # PRF liquid-volume % -> mole fraction
    n_iso = on * 0.6919 / 114.2285
    n_hep = (100 - on) * 0.6837 / 100.2019
    x_iso = n_iso / (n_iso + n_hep)
    # dry air (mole fractions; same as the project's dry air) + humidity ratio W
    air = {"N2": 0.780848, "O2": 0.209476, "Ar": 0.009365, "CO2": 0.000311}
    W = grains * GRAIN / LB
    M_air = sum(air[k] * gas.molecular_weights[gas.species_index(k)] for k in air) / 1000
    n_h2o = W / 0.01801528 * M_air  # mol water per mol dry air
    # stoichiometric O2 per mol fuel: iso-octane 12.5, n-heptane 11
    o2_st = x_iso * 12.5 + (1 - x_iso) * 11.0
    n_fuel = phi * air["O2"] / o2_st  # mol fuel per mol dry air
    comp = dict(air)
    comp["H2O"] = n_h2o
    comp["C8H18,isooctane"] = n_fuel * x_iso
    comp["C7H16,n-heptane"] = n_fuel * (1 - x_iso)
    tot = sum(comp.values())
    # enthalpy of the inlet streams (J), per 'tot' mol of mixture
    gas.TPX = T_air, P_ATM, {k: v for k, v in comp.items() if k in air or k == "H2O"}
    H_air = gas.enthalpy_mole / 1000 * (1 + n_h2o)  # J per mol dry air (Cantera: J/kmol)
    gas.TPX = T_fuel, P_ATM, {"C8H18,isooctane": x_iso, "C7H16,n-heptane": 1 - x_iso}
    h_fg = x_iso * 35.1e3 + (1 - x_iso) * 36.4e3
    H_fuel = n_fuel * (gas.enthalpy_mole / 1000 - h_fg)
    H = H_air + H_fuel
    gas.TPX = T_air, P_ATM, comp
    gas.HP = H / tot / gas.mean_molecular_weight * 1000, P_ATM
    return dict(octaneNumber=on, phi=phi, T_air=T_air, T_fuel=T_fuel, x_isooctane=x_iso, T_mixture=float(gas.T))


def main() -> None:
    hum = []
    for T in (52.0 + 273.15, 38.0 + 273.15):  # RON intake air, MON intake air (ASTM D2699 / D2700)
        for g in (25.0, 37.5, 50.0):
            hum.append(humidity(g, T))
    counts = np.arange(400.0, 1401.0, 50.0)
    crs = np.array([choi_cr(c) for c in counts])
    heff = STROKE / (crs - 1.0)
    A = np.vstack([np.ones_like(counts), counts]).T
    coef, *_ = np.linalg.lstsq(A, heff, rcond=None)
    resid = heff - A @ coef
    out = dict(
        description=(
            "CFR reference numbers: ASTM intake humidity (25-50 grains/lb dry air) as humidity ratio and RH at the "
            "intake-air temperature (Cantera water psat); Choi et al. 2018 (SAE 2018-01-0848) CR vs digital-counter "
            "polynomial, implied clearance height S/(CR-1) and its linear fit. SI units."
        ),
        generator="tools/reference/mechanics_cfr_conditions.py",
        cantera_version=ct.__version__,
        humidity=hum,
        ronMixture=ron_mixture_temperature(),
        counter=dict(
            counts=counts.tolist(),
            compressionRatio=crs.tolist(),
            effectiveClearanceHeight=heff.tolist(),
            linearFit=dict(intercept=float(coef[0]), slopePerCount=float(coef[1]), maxAbsResidual=float(np.max(np.abs(resid)))),
            rigidRaise=rigid_raise_fit(),
        ),
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, separators=(",", ":")) + "\n")
    print(f"wrote {OUT}")
    print(json.dumps(out["humidity"], indent=1))
    print(out["counter"]["linearFit"])
    print(json.dumps(out["counter"]["rigidRaise"], indent=1))
    print(out["ronMixture"])


if __name__ == "__main__":
    main()
