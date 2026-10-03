"""Motored CFR cylinder pressure at RON-relevant conditions (Neste Oyj CFR engine, data in Bhattacharya, Keskinen,
Larmi, Pal, Kroyan, Sarjovaara, Kaario, 'Fuel-air mixing in motored CFR engine at research octane number (RON)
relevant condition', Aalto University / ANL / Neste, accepted manuscript OSTI 2536666, Fig. 5).

The engine was motored immediately after fired RON-type operation on PRF87 (spark plug deactivated; 250 cycles
averaged, hot walls), 600 rpm, 1 atm intake, intake air 35 C (the 52 C nominal RON IAT 'adjusted to 35 C for
fit-for-use qualification'). The CR is NOT stated; the engine's CR relation is CR = 6345/(1850 - DCR) + 1 (Sec. 2).
If the rig sat at the PRF87 guide-table cylinder height (counter 685 at 29.92 inHg) this gives CR 6.45 (Neste
relation) / 6.23 (project rigid-raise relation) - an ASSUMPTION, flagged below.

Source: tools/reference/cfr_sources/osti_2536666_bhattacharya_motored.pdf, Fig. 5 (xref 72, 631x551 px RASTER, MATLAB).
The experiment is a blue DOTTED line: every dot is a sample of the curve. Digitisation: blue colour mask, connected
components (dots) -> centroids. Main panel: x from the tick marks (0, 180, 360, 540, 720 CAD; 0.611 px/deg ->
+-1 deg), y from ticks (0, 5, 10, 15 bar; 24.6 px/bar -> +-0.05 bar). Inset (330-400 CAD, 8-12 bar): x/y from the tick
label centres (2.34 px/deg, 37.5 px/bar -> +-0.5 deg, +-0.03 bar). Legend dots removed. Crank angle converted from
the source 0-720 CAD (firing TDC = 360) to the project convention: theta = CAD - 360.
Pressure: bar ABSOLUTE (the intake-stroke level is ~0.8-1.0 bar).

Writes test/fixtures/cfr_neste_motored.json.
"""
from __future__ import annotations

import numpy as np

from cfr_data_common import BAR, T0, Axis, blobs, conditions, cr_neste, cr_rigid_raise, load_image, prf, source_pdf, write_fixture

PDF = source_pdf("osti_2536666_bhattacharya_motored.pdf")


def main() -> None:
    img = load_image(PDF, 72).astype(int)
    R, G, B = img[:, :, 0], img[:, :, 1], img[:, :, 2]
    blue = (B > 150) & (B - R > 80) & (B - G > 60)
    xm = Axis.fit([139.0, 249.5, 359.5, 469.5, 579.5], [0, 180, 360, 540, 720])
    ym = Axis.fit([388.5, 265.5, 142.0, 19.0], [0, 5, 10, 15])
    xi = Axis.fit([200.5, 317.5], [350, 400])  # inset tick-label centres ('350', '400')
    yi = Axis.fit([222.0, 72.0], [8.0, 12.0])  # inset tick-label centres ('8', '12')
    inset = (110, 60, 330, 240)
    legend_y = 60  # legend dots at y ~ 45
    main_pts, inset_pts = [], []
    for b in blobs(blue, min_area=3):
        if b["cy"] < legend_y:
            continue
        if inset[0] <= b["cx"] <= inset[2] and inset[1] <= b["cy"] <= inset[3]:
            if b["area"] > 25:  # two dots merged -> skip
                continue
            inset_pts.append((float(xi(b["cx"])), float(yi(b["cy"]))))
        elif b["area"] <= 25:
            main_pts.append((float(xm(b["cx"])), float(ym(b["cy"]))))
    main_pts.sort()
    inset_pts.sort()
    th_m = np.array([p[0] for p in main_pts]) - 360.0
    p_m = np.array([p[1] for p in main_pts])
    th_i = np.array([p[0] for p in inset_pts]) - 360.0
    p_i = np.array([p[1] for p in inset_pts])
    counter_guess = 685  # ASTM D2699 guide table RON 87.0 (cfr_astm_guide_tables.json)
    ds = {
        "key": "neste_motored_ron87",
        "kind": "pressureTraceMotored",
        "figure": "Fig. 5", "page": 15,
        "caption": "Motored (spark off, hot) cylinder pressure, 250-cycle average, Neste CFR engine at RON-relevant conditions",
        "conditions": conditions(
            engine="Neste Oyj CFR F1 (AVL indicating spark plug)", method="motored (after fired RON operation)",
            fuel=prf(87), rpm=600, intakePressurePa=101325.0, intakeAirTemperatureK=35 + T0,
            compressionRatio=None,
            compressionRatioSource=f"NOT STATED. Assumption: PRF87 guide-table height (counter {counter_guess}) -> CR "
                                   f"{cr_neste(counter_guess):.2f} with the Neste relation CR = 6345/(1850 - DCR) + 1, "
                                   f"{cr_rigid_raise(counter_guess):.2f} with the project rigid-raise relation",
            notes="Carburettor still delivers PRF87 (the fuel-air mixing study): the motored charge contains fuel "
                  "(phi ~1.2 in the CFD after IVC); intake mixture temperature after the carburettor 292.5 K in the 1-D "
                  "model; wall temperatures in the CFD taken from Pal 2018. Treat CR as unknown (fit it from the "
                  "compression line) or as the assumption above."),
        "extraction": {"method": __doc__.split("Digitisation:")[1].split("Crank angle")[0].strip().replace("\n", " "),
                       "nMain": int(len(th_m)), "nInset": int(len(th_i)),
                       "uncertainty": {"mainThetaDeg": 1.0, "mainPressureBar": 0.05, "insetThetaDeg": 0.5, "insetPressureBar": 0.03}},
        "data": {"main": {"thetaDeg": th_m, "pPa": p_m * BAR}, "insetPeak": {"thetaDeg": th_i, "pPa": p_i * BAR}},
        "derived": {"peakPressurePa": float(max(p_i.max(), p_m.max()) * BAR),
                    "peakThetaDeg": float(th_i[np.argmax(p_i)]) if len(p_i) else None},
    }
    out = {
        "description": "Motored CFR cylinder pressure (Neste engine, RON-relevant conditions, PRF87 in the carburettor). "
                       "Pa absolute, deg with 0 = firing TDC.",
        "generator": "tools/reference/cfr_data_neste_motored.py",
        "source": {"key": "Bhattacharya2024", "citation": "Bhattacharya A., Keskinen K., Larmi M., Pal P., Kroyan Y., Sarjovaara T., "
                                                          "Kaario O., 'Fuel-air mixing in motored CFR engine at research octane "
                                                          "number (RON) relevant condition', accepted manuscript OSTI 2536666",
                   "url": "https://www.osti.gov/servlets/purl/2536666",
                   "localPdf": "tools/reference/cfr_sources/osti_2536666_bhattacharya_motored.pdf",
                   "experimentReference": "Keskinen K. et al., 'The impact of octane number boosters on knock characteristics "
                                          "in a CFR engine', SAE 2022-01-1082 (not open access)"},
        "datasets": [ds],
    }
    write_fixture("cfr_neste_motored.json", out)
    print(len(th_m), len(th_i), "peak", ds["derived"])
    print(np.round(th_m, 1).tolist())
    print(np.round(p_m, 2).tolist())


if __name__ == "__main__":
    main()
