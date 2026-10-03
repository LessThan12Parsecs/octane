"""Shared definitions for the thermo generator/oracle scripts (no side effects when run).

- Species list parsed from src/physics/core/species.ts (single source of truth).
- Mapping from our species names to Cantera nasa_gas.yaml names.
- The TRUE standard-state pressure of the NASA TM-4513 polynomials (1 bar) and a builder
  for Cantera Species objects that carry it (Cantera's nasa_gas.yaml mislabels it 1 atm —
  see P_REF_DATA below). Other oracle scripts should build their phases with
  `nasa_species()` / `build_phase()` so they agree with src/physics/thermo exactly.
- Lennard-Jones transport parameters: gri30.yaml where available, literature for the
  liquid-fuel species that GRI-Mech 3.0 does not contain.
- IAPWS dilute-gas viscosity / thermal conductivity of H2O (used by transport.ts for H2O).
- A Cantera phase builder that uses exactly the same thermo + transport data as the
  TypeScript implementation (for the oracle fixtures).
"""
from __future__ import annotations

import json
import math
import os
import re

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
SPECIES_TS = os.path.join(ROOT, "src", "physics", "core", "species.ts")
FIXTURES = os.path.join(ROOT, "test", "fixtures")

# ---------------------------------------------------------------------------------------
# Standard-state pressure of the NASA-7 data
# ---------------------------------------------------------------------------------------
# McBride, Gordon & Reno (1993), NASA TM-4513, §"Thermodynamic data": "For gases this is
# ideal gas at the standard pressure of 10^5 Pa (1 bar)", and the report adopts "the
# reference pressure for the ideal gases as 1 bar rather than 1 atmosphere" (fetched
# 2026-09-29 from https://archive.org/stream/nasa_techdoc_19940013151). Cantera's
# nasa_gas.yaml was converted by ck2yaml from the CHEMKIN-format file, which has no
# reference-pressure field, so Cantera assigns its default of 1 atm (101325 Pa) — a
# mislabel. The entropies prove it: s°(298.15 K) of the NASA fits equals the CODATA
# 1-bar key values to ≤ 0.007 J/(mol K), whereas a 1-atm basis would sit
# R ln(1.01325) = 0.1094 J/(mol K) lower (gen_species_data.py re-checks this).
P_REF_DATA = 1.0e5
CANTERA_NASA_GAS_P_REF = 101325.0

# CODATA Key Values for Thermodynamics (Cox, Wagman & Medvedev 1989), S°(298.15 K) in
# J/(mol K) at the standard-state pressure 100000 Pa (1 bar); fetched 2026-09-29 from
# https://www.codata.info/resources/databases/key1.html . Keys = our species names.
CODATA_S298 = {
    "O2": 205.152, "N2": 191.609, "H2": 130.680, "H2O": 188.835, "CO": 197.660,
    "CO2": 213.785, "AR": 154.846, "H": 114.717, "O": 161.059, "N": 153.301,
}


def read_species() -> list[str]:
    """Parse `export const SPECIES = [...] as const` from core/species.ts."""
    src = open(SPECIES_TS, encoding="utf-8").read()
    m = re.search(r"export const SPECIES = \[(.*?)\] as const", src, re.S)
    if not m:
        raise RuntimeError("could not find SPECIES in species.ts")
    names = re.findall(r"'([^']+)'", m.group(1))
    if not names:
        raise RuntimeError("empty SPECIES list")
    return names


# Our name -> name in Cantera's nasa_gas.yaml (NASA Glenn / McBride et al. 1993).
NASA_NAME = {
    "AR": "Ar",
    "IC8H18": "C8H18,isooctane",
    "NC7H16": "C7H16,n-heptane",
}


def nasa_name(sp: str) -> str:
    return NASA_NAME.get(sp, sp)


# Our name -> name in gri30.yaml (transport data source).
GRI_NAME = {"AR": "AR"}

# Transport parameters (CHEMKIN customary units) for species absent from GRI-Mech 3.0.
# Source: LLNL PRF mechanism transport file prf_tran_dat_v1b.txt (LLNL-MI-421252),
# Curran, Pitz & Westbrook, Combust. Flame 129 (2002) 253; identical entries in the LLNL
# n-heptane v3.1 transport file (LLNL-MI-536371, Mehl et al. 2011). Fetched 2026-09-29 from
# https://combustion.llnl.gov/sites/combustion/files/prf_tran_dat_v1b.txt :
#   IC8H18  2  458.5  6.414  0.0  0.0  1.0   ! WJP
#   NC7H16  2  459.6  6.253  0.0  0.0  1.0   ! TCPC
#   C2H5OH  2  470.6  4.410  0.0  0.0  1.5   ! NMM (Marinov, Int. J. Chem. Kinet. 31 (1999) 183)
LITERATURE_TRANSPORT = {
    "IC8H18": dict(geometry="nonlinear", well_depth=458.5, diameter=6.414, dipole=0.0,
                   polarizability=0.0, rotational_relaxation=1.0,
                   source="LLNL prf_tran_dat_v1b.txt (Curran et al. 2002), entry 'WJP'"),
    "NC7H16": dict(geometry="nonlinear", well_depth=459.6, diameter=6.253, dipole=0.0,
                   polarizability=0.0, rotational_relaxation=1.0,
                   source="LLNL prf_tran_dat_v1b.txt (Curran et al. 2002), entry 'TCPC'"),
    "C2H5OH": dict(geometry="nonlinear", well_depth=470.6, diameter=4.41, dipole=0.0,
                   polarizability=0.0, rotational_relaxation=1.5,
                   source="LLNL prf_tran_dat_v1b.txt, entry 'NMM' (Marinov 1999)"),
}

GEOMETRY_CODE = {"atom": 0, "linear": 1, "nonlinear": 2}


def _clean(x: float) -> float:
    """gri30 input_data round-trips through SI; recover the CK-file decimal (6 s.f.)."""
    return float(f"{x:.6g}")


def transport_params(sp: str, gri) -> dict:
    """CHEMKIN-unit transport parameters for species `sp` (gri = ct.Solution('gri30.yaml'))."""
    if sp in LITERATURE_TRANSPORT:
        return dict(LITERATURE_TRANSPORT[sp])
    gname = GRI_NAME.get(sp, sp)
    d = gri.species(gname).input_data["transport"]
    return dict(
        geometry=d["geometry"],
        well_depth=_clean(d["well-depth"]),
        diameter=_clean(d["diameter"]),
        dipole=_clean(d.get("dipole", 0.0)),
        polarizability=_clean(d.get("polarizability", 0.0)),
        rotational_relaxation=_clean(d.get("rotational-relaxation", 0.0)),
        source="gri30.yaml (GRI-Mech 3.0 transport data)",
    )


# ---------------------------------------------------------------------------------------
# IAPWS dilute-gas transport of H2O (same formulas as transport.ts)
# ---------------------------------------------------------------------------------------
# IAPWS R12-08 (2008), viscosity, eq. (11), Table 1; IAPWS R15-11 (2011), thermal
# conductivity, eq. (16), Table 1; T* = 647.096 K, mu* = 1e-6 Pa s, lambda* = 1e-3 W/(m K).
# Direct-use range in transport.ts (IAPWS: dilute-gas terms physically reasonable 250–2500 K);
# outside it the kinetic-theory value is scaled by its IAPWS ratio at the nearer limit.
IAPWS_T_MIN = 250.0
IAPWS_T_MAX = 2500.0
IAPWS_H = (1.67752, 2.20462, 0.6366564, -0.241605)
IAPWS_L = (2.443221e-3, 1.323095e-2, 6.770357e-3, -3.454586e-3, 4.096266e-4)


def water_viscosity_dilute(T: float) -> float:
    """IAPWS 2008 dilute-gas viscosity of H2O, Pa s."""
    t = T / 647.096
    return 1e-6 * 100.0 * math.sqrt(t) / sum(h / t ** i for i, h in enumerate(IAPWS_H))


def water_conductivity_dilute(T: float) -> float:
    """IAPWS 2011 dilute-gas thermal conductivity of H2O, W/(m K)."""
    t = T / 647.096
    return 1e-3 * math.sqrt(t) / sum(c / t ** k for k, c in enumerate(IAPWS_L))


# ---------------------------------------------------------------------------------------
# Cantera phase builders
# ---------------------------------------------------------------------------------------

def nasa_species(names: list[str] | None = None, trange: tuple[float, float] | None = None,
                 transport: bool = False):
    """Cantera Species objects (our names) with the nasa_gas.yaml NASA-7 coefficients and
    the CORRECT standard-state pressure P_REF_DATA (1 bar). Optionally narrow the declared
    validity range (same coefficients) and attach our LJ transport parameters."""
    import cantera as ct

    names = read_species() if names is None else names
    nasa = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    gri = ct.Solution("gri30.yaml") if transport else None
    out = []
    for sp in names:
        src = nasa[nasa_name(sp)]
        th = src.thermo
        s = ct.Species(sp, src.composition)
        tmin, tmax = th.min_temp, th.max_temp
        if trange is not None:
            tmin, tmax = max(trange[0], tmin), min(trange[1], tmax)
        s.thermo = ct.NasaPoly2(tmin, tmax, P_REF_DATA, th.coeffs)
        if transport:
            p = transport_params(sp, gri)
            tr = ct.GasTransportData()
            tr.set_customary_units(p["geometry"], p["diameter"], p["well_depth"],
                                   p["dipole"], p["polarizability"],
                                   p["rotational_relaxation"])
            s.transport = tr
        out.append(s)
    return out


def build_phase(transport: bool = False, trange: tuple[float, float] | None = None):
    """Cantera ideal-gas phase with our species, NASA-7 thermo from nasa_gas.yaml at the
    correct 1-bar standard state and (optionally) mixture-averaged transport with the same
    LJ parameters as the TS code.

    trange: optionally narrow the species' declared validity range (same coefficients, so
    values inside the range are unchanged). Cantera fits its pure-species transport
    properties with degree-4 polynomials in ln T over the phase's [Tmin, Tmax]; over the
    full 200–6000 K the conductivity fit error reaches ~6 % for the large fuel molecules at
    250 K, so the transport oracle uses a narrower range to keep the oracle itself accurate.
    """
    import cantera as ct

    kw = dict(thermo="ideal-gas", species=nasa_species(trange=trange, transport=transport))
    if transport:
        kw["transport_model"] = "mixture-averaged"
    return ct.Solution(**kw)


def write_fixture(name: str, data) -> None:
    os.makedirs(FIXTURES, exist_ok=True)
    path = os.path.join(FIXTURES, name)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"), allow_nan=False)
        f.write("\n")
    print("wrote", os.path.relpath(path, ROOT))


def rnd(x: float) -> float:
    """Full double precision (repr round-trips); NaN/inf rejected by json."""
    if not math.isfinite(x):
        raise ValueError(x)
    return float(x)


if __name__ == "__main__":
    # Library module; running it (e.g. from generate_all.sh) just prints the species list.
    print("thermo_common: species =", read_species())
