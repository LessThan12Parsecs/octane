"""Generate src/physics/thermo/species-data.ts from Cantera's bundled nasa_gas.yaml.

NASA-7 coefficients (NASA Glenn database; McBride, Gordon & Reno, NASA TM-4513, 1993) for
every species in core/species.ts SPECIES, in SPECIES order, plus the Lennard-Jones transport
parameters used by thermo/transport.ts (gri30.yaml where available, literature otherwise;
see thermo_common.py).

Checks (the script errors out rather than emitting inconsistent data):
  - every species uses the NASA7 model with 1 or 2 temperature regions (a single region is
    emitted with Tmid = Tmax and high = low, exactly as Cantera stores it);
  - all species share one reference pressure in Cantera;
  - the TRUE standard-state pressure of the data is determined from the data themselves:
    s°(298.15 K) of the species that have CODATA key values must match CODATA's 1-bar
    values (→ 1 bar) or sit R ln(1.01325) below them (→ 1 atm); anything else is an error.
    For nasa_gas.yaml this gives 1 bar (NASA TM-4513 states 1 bar), although Cantera labels
    the entries 1 atm (ck2yaml default). The data value is emitted as P_REF_NASA, Cantera's
    label as CANTERA_NASA_GAS_P_REF (see thermo_common.P_REF_DATA);
  - the Cantera composition is emitted so the TS test can check core ELEMENT_COUNTS.

Run: .venv/bin/python tools/reference/gen_species_data.py
"""
from __future__ import annotations

import os
import sys
from math import log as math_log

import cantera as ct

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from thermo_common import (  # noqa: E402
    CANTERA_NASA_GAS_P_REF, CODATA_S298, GEOMETRY_CODE, P_REF_DATA, ROOT, nasa_name, read_species,
    transport_params)

OUT = os.path.join(ROOT, "src", "physics", "thermo", "species-data.ts")


def fmt(x: float) -> str:
    r = repr(float(x))
    return r[:-2] if r.endswith(".0") and "e" not in r else r


def q(s: str) -> str:
    """Single-quoted TS string literal."""
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def fmt_list(xs) -> str:
    return "[" + ", ".join(fmt(x) for x in xs) + "]"


def main() -> None:
    species = read_species()
    nasa = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    gri = ct.Solution("gri30.yaml")

    p_ref = None
    s298: dict[str, float] = {}
    thermo_rows = []
    transport_rows = []
    for sp in species:
        cname = nasa_name(sp)
        if cname not in nasa:
            raise SystemExit(f"species {sp} ({cname}) not in nasa_gas.yaml")
        s = nasa[cname]
        inp = s.input_data["thermo"]
        if inp["model"] != "NASA7":
            raise SystemExit(f"{sp}: unsupported thermo model {inp['model']}")
        ranges = [float(t) for t in inp["temperature-ranges"]]
        data = [[float(c) for c in row] for row in inp["data"]]
        if len(ranges) == 2 and len(data) == 1:
            tmin, tmax = ranges
            tmid = tmax
            low = high = data[0]
        elif len(ranges) == 3 and len(data) == 2:
            tmin, tmid, tmax = ranges
            low, high = data
        else:
            raise SystemExit(
                f"{sp}: {len(ranges) - 1} temperature regions; thermo.ts supports 1 or 2 "
                "(extend the flat layout before adding this species)")
        # Cross-check against Cantera's own interpretation of the entry.
        th = s.thermo
        if not isinstance(th, ct.NasaPoly2):
            raise SystemExit(f"{sp}: Cantera did not build a NasaPoly2")
        c = th.coeffs  # [Tmid, high(7), low(7)]
        if abs(c[0] - tmid) > 0 or list(c[1:8]) != high or list(c[8:15]) != low:
            raise SystemExit(f"{sp}: coefficient cross-check against Cantera failed")
        if abs(th.min_temp - tmin) > 0 or abs(th.max_temp - tmax) > 0:
            raise SystemExit(f"{sp}: range cross-check failed")
        if p_ref is None:
            p_ref = th.reference_pressure
        elif th.reference_pressure != p_ref:
            raise SystemExit(f"{sp}: reference pressure {th.reference_pressure} != {p_ref}")
        if sp in CODATA_S298:
            s298[sp] = th.s(298.15) / 1000.0  # J/(mol K) at the data's own standard state
        comp = {k: int(round(v)) for k, v in s.composition.items()}
        thermo_rows.append(dict(name=sp, cname=cname, note=str(inp.get("note", "")).strip(),
                                tmin=tmin, tmid=tmid, tmax=tmax, low=low, high=high,
                                comp=comp, mw=s.molecular_weight / 1000.0))
        transport_rows.append((sp, transport_params(sp, gri)))

    # ---- the data's true standard-state pressure (CODATA 1-bar entropies) ----
    if p_ref != CANTERA_NASA_GAS_P_REF:
        raise SystemExit(f"Cantera now labels nasa_gas.yaml p_ref = {p_ref}; re-check P_REF_DATA")
    if len(s298) < 5:
        raise SystemExit("too few CODATA species to determine the standard-state pressure")
    r_u = ct.gas_constant / 1000.0
    shift = r_u * math_log(101325.0 / 1.0e5)  # s°(1 atm) = s°(1 bar) − shift
    dev_bar = max(abs(s298[k] - CODATA_S298[k]) for k in s298)
    dev_atm = max(abs(s298[k] - (CODATA_S298[k] - shift)) for k in s298)
    if dev_bar < 0.03 and dev_atm > 0.08:
        p_data = 1.0e5
    elif dev_atm < 0.03 and dev_bar > 0.08:
        p_data = 101325.0
    else:
        raise SystemExit(f"cannot determine standard-state pressure: dev(1 bar) {dev_bar}, "
                         f"dev(1 atm) {dev_atm} J/(mol K)")
    if p_data != P_REF_DATA:
        raise SystemExit(f"data standard state {p_data} Pa != thermo_common.P_REF_DATA")
    print(f"standard state: max |s°298 − CODATA(1 bar)| = {dev_bar:.4f} J/(mol K), "
          f"vs 1-atm hypothesis {dev_atm:.4f} → p_ref = {p_data:g} Pa (Cantera label {p_ref:g})")

    lines = []
    w = lines.append
    w("// GENERATED — do not edit.")
    w("// Generated by tools/reference/gen_species_data.py from Cantera "
      f"{ct.__version__} data files:")
    w("//   thermo:    nasa_gas.yaml — NASA Glenn database, NASA-7 fits (B.J. McBride, S. Gordon,")
    w("//              M.A. Reno, NASA TM-4513, 1993). `note` = source/date code of the fit.")
    w("//   transport: gri30.yaml (GRI-Mech 3.0) Lennard-Jones data; liquid-fuel species from the")
    w("//              LLNL PRF transport file (Curran, Pitz & Westbrook 2002) — see `source`.")
    w("// Regenerate: .venv/bin/python tools/reference/gen_species_data.py")
    w("")
    w("/**")
    w(" * Standard-state pressure of every polynomial below, Pa: 1 bar. NASA TM-4513 (McBride et al.")
    w(" * 1993) states the ideal-gas standard state is 10^5 Pa, and the generator verified it:")
    w(f" * s°(298.15 K) matches the CODATA 1-bar key values to {dev_bar:.4f} J/(mol K) (a 1-atm basis")
    w(f" * would be off by {shift:.4f}).")
    w(" */")
    w(f"export const P_REF_NASA = {fmt(p_data)};")
    w("")
    w("/**")
    w(" * Reference pressure that Cantera's nasa_gas.yaml ASSIGNS these same fits, Pa (ck2yaml's default")
    w(" * of 1 atm — a mislabel). Cantera oracles must rebuild the species with P_REF_NASA (see")
    w(" * tools/reference/thermo_common.py nasa_species()) to agree with this data set.")
    w(" */")
    w(f"export const CANTERA_NASA_GAS_P_REF = {fmt(p_ref)};")
    w("")
    w("/** One NASA-7 entry: cp/R = a0 + a1 T + a2 T² + a3 T³ + a4 T⁴,")
    w(" *  h/RT = a0 + a1 T/2 + a2 T²/3 + a3 T³/4 + a4 T⁴/5 + a5/T,")
    w(" *  s°/R = a0 ln T + a1 T + a2 T²/2 + a3 T³/3 + a4 T⁴/4 + a6.")
    w(" *  `low` applies for T ≤ tmid (Cantera convention), `high` for T > tmid. Temperatures in K. */")
    w("export interface Nasa7Record {")
    w("  /** Our species name (core/species.ts). */")
    w("  name: string;")
    w("  /** Name in nasa_gas.yaml. */")
    w("  canteraName: string;")
    w("  /** Source/date code of the fit (McBride et al. 1993). */")
    w("  note: string;")
    w("  tmin: number;")
    w("  tmid: number;")
    w("  tmax: number;")
    w("  low: readonly number[];")
    w("  high: readonly number[];")
    w("  /** Elemental composition as stored by Cantera (atoms per molecule). */")
    w("  composition: Readonly<Record<string, number>>;")
    w("  /** Cantera molecular weight, kg/mol (cross-check only; MOLAR_MASS is computed). */")
    w("  canteraMolarMass: number;")
    w("}")
    w("")
    w("/** NASA-7 data for every species, in SPECIES order. */")
    w("export const SPECIES_NASA7: readonly Nasa7Record[] = [")
    for r in thermo_rows:
        comp = ", ".join(f"{k}: {v}" for k, v in sorted(r["comp"].items()))
        w("  {")
        w(f"    name: {q(r['name'])}, canteraName: {q(r['cname'])}, note: {q(r['note'])},")
        w(f"    tmin: {fmt(r['tmin'])}, tmid: {fmt(r['tmid'])}, tmax: {fmt(r['tmax'])},")
        w(f"    low: {fmt_list(r['low'])},")
        w(f"    high: {fmt_list(r['high'])},")
        w(f"    composition: {{ {comp} }}, canteraMolarMass: {fmt(r['mw'])},")
        w("  },")
    w("];")
    w("")
    w("/** Lennard-Jones / Stockmayer transport parameters (CHEMKIN customary units). */")
    w("export interface TransportRecord {")
    w("  name: string;")
    w("  /** 0 = atom, 1 = linear, 2 = nonlinear. */")
    w("  geometry: 0 | 1 | 2;")
    w("  /** Potential well depth ε/k_B, K. */")
    w("  wellDepth: number;")
    w("  /** Collision diameter σ, Å (1e-10 m). */")
    w("  diameter: number;")
    w("  /** Permanent dipole moment, Debye. */")
    w("  dipole: number;")
    w("  /** Polarizability, Å³. */")
    w("  polarizability: number;")
    w("  /** Rotational relaxation collision number Z_rot at 298 K. */")
    w("  rotRelax: number;")
    w("  source: string;")
    w("}")
    w("")
    w("/** Transport parameters for every species, in SPECIES order. */")
    w("export const SPECIES_TRANSPORT: readonly TransportRecord[] = [")
    for sp, p in transport_rows:
        w(f"  {{ name: {q(sp)}, geometry: {GEOMETRY_CODE[p['geometry']]}, "
          f"wellDepth: {fmt(p['well_depth'])}, diameter: {fmt(p['diameter'])}, "
          f"dipole: {fmt(p['dipole'])}, polarizability: {fmt(p['polarizability'])}, "
          f"rotRelax: {fmt(p['rotational_relaxation'])}, source: {q(p['source'])} }},")
    w("];")
    w("")
    with open(OUT, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print("wrote", os.path.relpath(OUT, ROOT), f"({len(species)} species, p_ref = {p_data:g} Pa)")


if __name__ == "__main__":
    main()
