"""Convert the PRF_HT reduced iso-octane/n-heptane mechanism (Cantera .cti) to a JSON fixture
that Cantera 3.2 can load with ct.Solution(yaml=<text>) — WITHOUT executing the .cti file.

Mechanism: "PRF_HT" — the reduced high-temperature n-heptane/iso-octane mechanism (99 species,
601 reactions) of Jerzembeck, Peters, Pepiot-Desjardins & Pitsch, "Laminar burning velocities
at high pressure for primary reference fuels and gasoline: Experimental and numerical
investigation", Combust. Flame 156 (2009) 292-301, §4 (DRGEP reduction of the LLNL n-heptane
(Curran et al. 1998) and iso-octane (Curran et al. 2002) mechanisms; validated against laminar
burning velocities of C1-C4 alkanes and PRFs at 1 atm and the paper's own 10-25 bar / 373 K
spherical-bomb data). The .cti file is distributed by CERFACS at
https://www.cerfacs.fr/cantera/docs/mechanisms/iso-octane-air/Jerzembeck/PRF_HT.cti
(fetched 2026-09-29, 138 468 bytes, header "Mechanism contains 99 species and 601 reactions").

A .cti file is Python source; Cantera's cti2yaml runs it with exec(). To avoid executing a
downloaded file, this script parses it with `ast` and evaluates it with a restricted
interpreter that only allows literals, tuples/lists/dicts, unary minus, simple arithmetic on
numbers, the names OneAtm/OneBar, and calls to the cti2yaml DSL constructors. cti2yaml's
YAML writer is then reused unchanged.

Usage: .venv/bin/python tools/reference/flamespeed_prf_mechanism.py path/to/PRF_HT.cti
Writes test/fixtures/flamespeed_mech_prf_ht.json = {source, sha256, yaml}.
"""
from __future__ import annotations

import ast
import hashlib
import json
import operator
import os
import sys
import tempfile

import cantera as ct
from cantera import cti2yaml

ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
OUT = os.path.join(ROOT, "test", "fixtures", "flamespeed_mech_prf_ht.json")

ALLOWED_CALLS = {
    "units", "ideal_gas", "species", "NASA", "NASA9", "gas_transport", "state",
    "reaction", "three_body_reaction", "falloff_reaction", "Troe", "SRI", "Lindemann",
    "Arrhenius", "pdep_arrhenius", "chemically_activated_reaction", "standard_pressure",
}
ALLOWED_NAMES = {"OneAtm": cti2yaml.OneAtm, "OneBar": cti2yaml.OneBar}
BINOPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul,
          ast.Div: operator.truediv, ast.Pow: operator.pow}


def _eval(node: ast.AST):
    if isinstance(node, ast.Constant):
        if isinstance(node.value, (int, float, str, bool)) or node.value is None:
            return node.value
        raise ValueError(f"constant type {type(node.value)} not allowed")
    if isinstance(node, ast.Tuple):
        return tuple(_eval(e) for e in node.elts)
    if isinstance(node, ast.List):
        return [_eval(e) for e in node.elts]
    if isinstance(node, ast.Dict):
        return {_eval(k): _eval(v) for k, v in zip(node.keys, node.values)}
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.USub, ast.UAdd)):
        v = _eval(node.operand)
        if not isinstance(v, (int, float)):
            raise ValueError("unary op on non-number")
        return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.BinOp) and type(node.op) in BINOPS:
        a, b = _eval(node.left), _eval(node.right)
        if not (isinstance(a, (int, float)) and isinstance(b, (int, float))):
            raise ValueError("arithmetic on non-numbers")
        return BINOPS[type(node.op)](a, b)
    if isinstance(node, ast.Name):
        if node.id in ALLOWED_NAMES:
            return ALLOWED_NAMES[node.id]
        raise ValueError(f"name {node.id!r} not allowed")
    if isinstance(node, ast.Call):
        if not (isinstance(node.func, ast.Name) and node.func.id in ALLOWED_CALLS):
            raise ValueError(f"call {ast.dump(node.func)} not allowed")
        fn = getattr(cti2yaml, node.func.id)
        args = [_eval(a) for a in node.args]
        kwargs = {}
        for kw in node.keywords:
            if kw.arg is None:
                raise ValueError("**kwargs not allowed")
            kwargs[kw.arg] = _eval(kw.value)
        return fn(*args, **kwargs)
    raise ValueError(f"node {type(node).__name__} not allowed")


def safe_exec(text: str) -> None:
    tree = ast.parse(text, mode="exec")
    for stmt in tree.body:
        if not (isinstance(stmt, ast.Expr) and isinstance(stmt.value, ast.Call)):
            raise ValueError(f"line {stmt.lineno}: only top-level DSL calls are allowed")
        _eval(stmt.value)


# Linear molecules among the mechanism species (for species absent from GRI-Mech 3.0).
LINEAR = {"C2H", "C3H2", "HCCO"}


def fix_transport(yaml_text: str) -> str:
    """The distributed .cti keeps only the Lennard-Jones σ and ε/k of every species (the
    FlameMaster convention) and labels every species geom='atom', which Cantera 3 rejects for
    polyatomics. The mechanism's own σ and ε/k are kept. Species present in GRI-Mech 3.0
    (gri30.yaml) take geometry, dipole, polarizability and rotational relaxation from their
    gri30 records; the others get geometry from their structure (atom / linear / nonlinear), zero
    dipole and polarizability, and rotational relaxation Z_rot(298) = 1 (the value used for all
    C3+ species in the LLNL PRF transport file). Sensitivity of S_L to these secondary
    parameters (review check, iso-octane φ = 1, 300 K, 1 bar, 2-level refinement): Z_rot of the
    non-GRI species 1 → 5 or → 0: < 0.01 %; removing every dipole moment and polarizability
    (the gri30 values restored here): +0.66 %. Geometry follows from molecular structure."""
    from ruamel.yaml import YAML
    from io import StringIO
    y = YAML()
    doc = y.load(yaml_text)
    gri = {s.name: s for s in ct.Species.list_from_file("gri30.yaml")}
    for sp in doc["species"]:
        tr = sp["transport"]
        name = sp["name"]
        g = gri.get(name)
        if g is not None and g.transport is not None:
            gt = g.transport
            tr["geometry"] = gt.geometry
            if gt.dipole:
                tr["dipole"] = gt.dipole / 3.335640952e-30  # C m -> Debye
            if gt.polarizability:
                tr["polarizability"] = gt.polarizability * 1e30  # m^3 -> Å^3
            if gt.rotational_relaxation:
                tr["rotational-relaxation"] = gt.rotational_relaxation
        else:
            natoms = sum(sp["composition"].values())
            tr["geometry"] = "atom" if natoms == 1 else ("linear" if name in LINEAR or natoms == 2 else "nonlinear")
            if natoms > 1:
                tr["rotational-relaxation"] = 1.0
    buf = StringIO()
    y.dump(doc, buf)
    return buf.getvalue()


def main() -> None:
    src = sys.argv[1]
    raw = open(src, "rb").read()
    text = raw.decode("latin-1")
    # Route cti2yaml.convert's compile()/exec() through the restricted interpreter.
    cti2yaml.compile = lambda t, fn, mode: t  # type: ignore[attr-defined]
    cti2yaml.exec = safe_exec  # type: ignore[attr-defined]
    with tempfile.TemporaryDirectory() as tmp:
        out_yaml = os.path.join(tmp, "prf_ht.yaml")
        cti2yaml.convert(text=text, output_name=out_yaml)
        yaml_text = open(out_yaml, encoding="utf-8").read()
    yaml_text = fix_transport(yaml_text)
    gas = ct.Solution(yaml=yaml_text, name="gas")
    print(f"converted: {gas.n_species} species, {gas.n_reactions} reactions")
    data = dict(
        source="PRF_HT reduced mechanism of Jerzembeck, Peters, Pepiot-Desjardins & Pitsch, "
               "Combust. Flame 156 (2009) 292-301; .cti from https://www.cerfacs.fr/cantera/docs/"
               "mechanisms/iso-octane-air/Jerzembeck/PRF_HT.cti (fetched 2026-09-29), converted "
               "with tools/reference/flamespeed_prf_mechanism.py (restricted AST evaluation + "
               "Cantera %s cti2yaml writer)" % ct.__version__,
        sha256=hashlib.sha256(raw).hexdigest(),
        nSpecies=gas.n_species,
        nReactions=gas.n_reactions,
        yaml=yaml_text,
    )
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(data, fh, separators=(",", ":"))
        fh.write("\n")
    print(f"wrote {os.path.relpath(OUT, ROOT)}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Part of generate_all.sh: the mechanism fixture is committed; regenerate only on demand.
        print("flamespeed_prf_mechanism: pass the path to PRF_HT.cti to regenerate the fixture")
        sys.exit(0)
    main()
