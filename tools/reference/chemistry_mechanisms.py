"""Obtain / convert the kinetic mechanisms used by the chemistry oracles (idempotent).

1. LLNL detailed PRF mechanism, version 2 (iso-octane + n-heptane):
     Curran, H.J., Gaffuri, P., Pitz, W.J., Westbrook, C.K., "A comprehensive modeling study
     of iso-octane oxidation", Combust. Flame 129 (2002) 253-280 (rate rules);
     Curran, Pitz, Westbrook, Callahan, Dryer, Proc. Combust. Inst. 27 (1998) 379-387 (PRF);
     file header: "Curran, H. J., Pitz, W. J., and Westbrook, C. K., 2002, UCRL-WEB-208393,
     Review and release date: December 3, 2004" (prf_2d, 9/13/04 WJP).
   Source: https://combustion.llnl.gov/archived-mechanisms/surrogates/prf-isooctane-n-heptane-mixture
     prf_2d_mech.txt, prf_2d_therm.txt  (fetched 2026-09-29).
   License: publicly released by LLNL (UCRL-WEB-208393); no explicit license text on the
   page. Cite the references above when using it.
   Converted with Cantera 3.2 ck2yaml (--permissive; transport not needed for 0-D).
   Post-processing: 'explicit-third-body-duplicates: mark-duplicate' is added to the phase
   so that the two LLNL reactions written with an explicit spectator H/OH
   (hocho + h => h2 + co2 + h, hocho + oh => h2o + co + oh) stay separate bimolecular
   reactions next to the genuine hocho(+M) decompositions (Cantera otherwise warns that they
   look like duplicate third-body reactions).

1b. LLNL gasoline-surrogate detailed mechanism (Mehl, Pitz, Westbrook & Curran, Proc. Combust.
   Inst. 33 (2011) 193; LLNL-MI-536371, ver. 1.0 2012-03-30; 1382 species after the patch
   below). It contains the UPDATED iso-octane and n-heptane sub-mechanisms (Mehl et al. 2009/
   2011 rate rules) and is the default for the tabulated PRF ignition delays because it
   reproduces the Fieweger et al. (1997) n-heptane shock-tube data better than PRF v2 (see
   chemistry_ignition_validation.py). Source:
   https://combustion.llnl.gov/mechanisms/surrogates/gasoline-surrogate
     ChemDetailed.inp.txt, gasoline_surrogate_therm.dat.txt (fetched 2026-09-29).
   Patch: 7 pentene low-T species have no thermo data in the published thermo file; they and
   the 32 reactions involving them are removed (GS_MISSING_THERMO); the Fortran-style exponent
   '0.0000+03' is read as 0.

2. Thermal-NO (extended Zeldovich) sub-mechanism of GRI-Mech 3.0 (Smith et al. 1999,
   http://combustion.berkeley.edu/gri-mech/), extracted from Cantera's bundled gri30.yaml:
     R178  N + NO <=> N2 + O      A = 2.7e13 cm3/mol/s, b = 0,    Ea = 355 cal/mol
     R179  N + O2 <=> NO + O      A = 9.0e9,            b = 1.0,  Ea = 6500
     R180  N + OH <=> NO + H      A = 3.36e13,          b = 0,    Ea = 385
   (gri30.yaml equation indices 177, 178, 179; values are checked by the script).
   Species: all 53 GRI species are kept (thermo only matters for N2 O2 O N NO OH H), so
   reverse rates use GRI-3.0 thermochemistry exactly as full GRI-3.0 does.
"""
from __future__ import annotations

import os
import subprocess
import sys
import urllib.request

from chemistry_common import LLNL_GS_DIR, LLNL_GS_YAML, LLNL_PRF_DIR, LLNL_PRF_YAML, ZELDOVICH_YAML

LLNL_BASE = "https://combustion.llnl.gov/sites/combustion/files/"
LLNL_FILES = ["prf_2d_mech.txt", "prf_2d_therm.txt", "prf_thermo_readme.txt"]


def fetch(url: str, dest: str) -> None:
    if os.path.exists(dest):
        return
    print("downloading", url)
    req = urllib.request.Request(url, headers={"User-Agent": "curl/8.7.1"})  # default UA gets 403
    with urllib.request.urlopen(req, timeout=60) as r, open(dest, "wb") as f:
        f.write(r.read())


def convert_llnl_prf() -> None:
    os.makedirs(LLNL_PRF_DIR, exist_ok=True)
    for name in LLNL_FILES:
        fetch(LLNL_BASE + name, os.path.join(LLNL_PRF_DIR, name))
    if os.path.exists(LLNL_PRF_YAML):
        return
    cmd = [sys.executable, "-m", "cantera.ck2yaml",
           "--input=" + os.path.join(LLNL_PRF_DIR, "prf_2d_mech.txt"),
           "--thermo=" + os.path.join(LLNL_PRF_DIR, "prf_2d_therm.txt"),
           "--output=" + LLNL_PRF_YAML, "--permissive", "--quiet", "--no-validate"]
    subprocess.run(cmd, check=True)
    txt = open(LLNL_PRF_YAML, encoding="utf-8").read()
    marker = "  kinetics: gas\n"
    if "explicit-third-body-duplicates" not in txt:
        txt = txt.replace(marker, marker + "  explicit-third-body-duplicates: mark-duplicate\n", 1)
    header = (
        "# LLNL detailed PRF mechanism v2 (prf_2d, Curran, Pitz & Westbrook, UCRL-WEB-208393,\n"
        "# released 2004-12-03; rate rules: Curran et al., Combust. Flame 129 (2002) 253).\n"
        "# Source: " + LLNL_BASE + "prf_2d_mech.txt / prf_2d_therm.txt\n"
        "# Converted by tools/reference/chemistry_mechanisms.py (Cantera ck2yaml, permissive).\n")
    with open(LLNL_PRF_YAML, "w", encoding="utf-8") as f:
        f.write(header + txt)
    print("wrote", LLNL_PRF_YAML)


GS_FILES = ["ChemDetailed.inp.txt", "gasoline_surrogate_therm.dat.txt"]
# Species of the 1-/2-pentene low-temperature sub-mechanism that appear in ChemDetailed.inp.txt
# but have no entry in the published gasoline_surrogate_therm.dat.txt (ck2yaml: "No thermo data
# found"); they and the 32 reactions that involve them are removed. They cannot affect PRF
# (iso-octane / n-heptane) chemistry except through negligible pentene side channels.
GS_MISSING_THERMO = ["C5H81OOH5-4", "C5H92O2-1", "C5H81OOH3-5O2", "C5H81OOH4-5O2",
                     "C5H81OOH5-4O2", "NC5D1KET34", "CY3C5H8O"]


def _species_tokens(eq: str) -> list[str]:
    import re

    eq = re.sub(r"\(\+[^)]*\)", "", eq.split("!")[0])
    toks = []
    for part in re.split(r"<=>|=>|=", eq):
        for t in part.split("+"):
            t = re.sub(r"^[0-9.]+", "", t.strip())
            if t:
                toks.append(t.upper())
    return toks


def patch_chemkin(src: str, dst: str, drop_species: list[str]) -> int:
    """Copy a Chemkin mechanism, fixing the Fortran-style '0.0000+03' exponent typo (= 0) and
    removing `drop_species` from SPECIES and every reaction (with its auxiliary lines) that
    involves them. Returns the number of reactions removed."""
    bad = {s.upper() for s in drop_species}
    out = []
    section = None
    drop = False
    n = 0
    for line in open(src, encoding="latin-1").read().split("\n"):
        line = line.replace("0.0000+03", "0.0000E+03")
        s = line.split("!")[0].strip()
        up = s.upper()
        if up.startswith("SPECIES"):
            section = "sp"
        elif up.startswith("REACTIONS"):
            section = "rx"
        elif up == "END":
            section = None
        elif section == "sp" and s:
            line = " ".join(t for t in s.split() if t.upper() not in bad)
        elif section == "rx":
            if "=" in s:
                fields = s.split()
                drop = any(t in bad for t in _species_tokens(" ".join(fields[:-3])))
                n += drop
            if drop:
                continue
        out.append(line)
    with open(dst, "w", encoding="utf-8") as f:
        f.write("\n".join(out))
    return n


def convert_llnl_gasoline() -> None:
    """LLNL gasoline-surrogate detailed mechanism (Mehl, Pitz, Westbrook & Curran, "Kinetic
    modeling of gasoline surrogate components and mixtures under engine conditions", Proc.
    Combust. Inst. 33 (2011) 193-200; LLNL-MI-536371, 'ver. 1.0 2012-03-30'), which contains the
    updated iso-octane / n-heptane (PRF) sub-mechanisms. Source:
    https://combustion.llnl.gov/mechanisms/surrogates/gasoline-surrogate (fetched 2026-09-29)."""
    os.makedirs(LLNL_GS_DIR, exist_ok=True)
    for name in GS_FILES:
        fetch(LLNL_BASE + name, os.path.join(LLNL_GS_DIR, name))
    if os.path.exists(LLNL_GS_YAML):
        return
    patched = os.path.join(LLNL_GS_DIR, "ChemDetailed_patched.inp")
    n = patch_chemkin(os.path.join(LLNL_GS_DIR, "ChemDetailed.inp.txt"), patched, GS_MISSING_THERMO)
    print(f"patched: removed {len(GS_MISSING_THERMO)} species without thermo and {n} reactions")
    cmd = [sys.executable, "-m", "cantera.ck2yaml", "--input=" + patched,
           "--thermo=" + os.path.join(LLNL_GS_DIR, "gasoline_surrogate_therm.dat.txt"),
           "--output=" + LLNL_GS_YAML, "--permissive", "--quiet", "--no-validate"]
    subprocess.run(cmd, check=True)
    os.remove(patched)
    txt = open(LLNL_GS_YAML, encoding="utf-8").read()
    header = (
        "# LLNL gasoline-surrogate detailed mechanism ver. 1.0 (2012-03-30), Mehl, Pitz, Westbrook &\n"
        "# Curran, Proc. Combust. Inst. 33 (2011) 193; LLNL-MI-536371.\n"
        "# Source: " + LLNL_BASE + "ChemDetailed.inp.txt / gasoline_surrogate_therm.dat.txt\n"
        "# Patched + converted by tools/reference/chemistry_mechanisms.py: removed pentene LT species\n"
        "# without thermo data (" + ", ".join(GS_MISSING_THERMO) + ") and their reactions.\n")
    with open(LLNL_GS_YAML, "w", encoding="utf-8") as f:
        f.write(header + txt)
    print("wrote", LLNL_GS_YAML)


def extract_zeldovich() -> None:
    import cantera as ct

    if os.path.exists(ZELDOVICH_YAML):
        return
    gri = ct.Solution("gri30.yaml")
    idx = [177, 178, 179]
    rxns = [gri.reaction(i) for i in idx]
    eqs = [r.equation for r in rxns]
    assert eqs == ["N + NO <=> N2 + O", "N + O2 <=> NO + O", "N + OH <=> H + NO"], eqs
    # GRI-Mech 3.0 values (cm3/mol/s, cal/mol) -> Cantera SI (m3/kmol/s, J/kmol)
    expect = [(2.7e13, 0.0, 355.0), (9.0e9, 1.0, 6500.0), (3.36e13, 0.0, 385.0)]
    for r, (A, b, Ea) in zip(rxns, expect):
        k = r.rate
        assert abs(k.pre_exponential_factor / (A * 1e-3) - 1) < 1e-12
        assert k.temperature_exponent == b
        assert abs(k.activation_energy / (Ea * 4184.0) - 1) < 1e-12
    sol = ct.Solution(thermo="ideal-gas", kinetics="gas", species=gri.species(), reactions=rxns,
                      name="zeldovich")
    sol.write_yaml(ZELDOVICH_YAML, header=True)
    txt = open(ZELDOVICH_YAML, encoding="utf-8").read()
    head = ("# Extended-Zeldovich (thermal NO) subset of GRI-Mech 3.0: gri30.yaml reactions\n"
            "# 178-180 (1-based) with all GRI-3.0 species thermo. Generated by\n"
            "# tools/reference/chemistry_mechanisms.py.\n")
    with open(ZELDOVICH_YAML, "w", encoding="utf-8") as f:
        f.write(head + txt)
    print("wrote", ZELDOVICH_YAML)


if __name__ == "__main__":
    convert_llnl_prf()
    convert_llnl_gasoline()
    extract_zeldovich()
