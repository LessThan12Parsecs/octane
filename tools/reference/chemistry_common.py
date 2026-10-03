"""Shared helpers for the chemistry oracle scripts (no side effects on import).

- Mechanism paths (tools/reference/mechanisms/...).
- PRF fuel / fresh-charge / residual compositions built EXACTLY like src/physics/thermo/fuels.ts
  (dry air of Picard et al. 2008, PRF liquid-volume -> mole fraction with the 60 F liquid
  densities, complete-combustion residual with the water-gas shift frozen at 1740 K).
- A constant-volume adiabatic ignition-delay routine (Cantera IdealGasMoleReactor +
  AdaptivePreconditioner) returning the main (max dT/dt) and first-stage delays.
"""
from __future__ import annotations

import hashlib
import json
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
FIXTURES = os.path.join(ROOT, "test", "fixtures")
MECH_DIR = os.path.join(HERE, "mechanisms")
LLNL_PRF_DIR = os.path.join(MECH_DIR, "llnl_prf_v2")
LLNL_PRF_YAML = os.path.join(LLNL_PRF_DIR, "llnl_prf_v2.yaml")
LLNL_GS_DIR = os.path.join(MECH_DIR, "llnl_gasoline_surrogate")
LLNL_GS_YAML = os.path.join(LLNL_GS_DIR, "llnl_gasoline_surrogate_2011.yaml")
ZELDOVICH_YAML = os.path.join(MECH_DIR, "zeldovich_gri30.yaml")
CHEM_DATA_DIR = os.path.join(ROOT, "src", "physics", "chemistry", "data")

R_UNIVERSAL = 8.31446261815324
P_ATM = 101325.0

# Atomic weights, kg/mol (src/physics/core/constants.ts, Cantera 3.x values).
AW = {"C": 12.011e-3, "H": 1.008e-3, "O": 15.999e-3, "N": 14.007e-3, "AR": 39.95e-3}
M_IC8H18 = 8 * AW["C"] + 18 * AW["H"]
M_NC7H16 = 7 * AW["C"] + 16 * AW["H"]

# Liquid densities at 60 F (src/physics/thermo/fuels.ts): NIST SRM 2214 eq. (1) for
# iso-octane, NIST WebBook (Tenji et al. 2018 EOS) for n-heptane.
_T60F_C = (60 - 32) / 1.8
_d = _T60F_C - 20
RHO_ISO_60F = 691.872 * (1 - 1.18752e-3 * _d - 6.2516e-7 * _d * _d)
RHO_HEP_60F = 687.597

# Dry air (Picard et al. 2008, CIPM-2007 Table 1; as src/physics/thermo/fuels.ts dryAir()).
_AIR = {"N2": 0.780848, "O2": 0.20939, "AR": 0.009332, "CO2": 0.0004}
_s = sum(_AIR.values())
DRY_AIR = {k: v / _s for k, v in _AIR.items()}

WGS_FREEZE_TEMPERATURE = 1740.0  # K, as src/physics/thermo/fuels.ts


def prf_isooctane_mole_fraction(on: float) -> float:
    """Mole fraction of iso-octane in the vapour of PRF `on` (liquid-volume % iso-octane)."""
    v = on / 100.0
    n_iso = v * RHO_ISO_60F / M_IC8H18
    n_hep = (1 - v) * RHO_HEP_60F / M_NC7H16
    return n_iso / (n_iso + n_hep)


def _nasa_species():
    import cantera as ct

    return {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}


_WGS_K = None


def water_gas_shift_k(T: float = WGS_FREEZE_TEMPERATURE) -> float:
    """K = x_CO x_H2O / (x_CO2 x_H2) from nasa_gas.yaml (same data as the TS thermo)."""
    import cantera as ct

    sp = _nasa_species()
    g = ct.Solution(thermo="ideal-gas", species=[sp[n] for n in ("CO", "H2O", "CO2", "H2")])
    g.TP = T, P_ATM
    mu0 = g.standard_gibbs_RT  # g°/RT per species (J/kmol basis is irrelevant: dimensionless)
    return math.exp(-(mu0[0] + mu0[1] - mu0[2] - mu0[3]))


def fresh_charge(phi: float, on: float) -> dict[str, float]:
    """Air + PRF vapour at equivalence ratio phi (mole fractions, project species names).

    n_fuel / n_air = phi x_O2,air / nu,  nu = C + H/4 per mean fuel molecule (fuels.ts)."""
    x_iso = prf_isooctane_mole_fraction(on)
    nu = x_iso * (8 + 18 / 4) + (1 - x_iso) * (7 + 16 / 4)
    n_fuel = phi * DRY_AIR["O2"] / nu
    X = dict(DRY_AIR)
    if x_iso > 0:
        X["IC8H18"] = n_fuel * x_iso
    if x_iso < 1:
        X["NC7H16"] = n_fuel * (1 - x_iso)
    s = sum(X.values())
    return {k: v / s for k, v in X.items()}


def complete_combustion_products(X: dict[str, float]) -> dict[str, float]:
    """Complete-combustion products (mole fractions) exactly as fuels.ts
    completeCombustionMoles for C/H/O/N/Ar reactants with O >= C."""
    counts = {
        "N2": (0, 0, 0, 2, 0), "O2": (0, 0, 2, 0, 0), "AR": (0, 0, 0, 0, 1),
        "CO2": (1, 0, 2, 0, 0), "H2O": (0, 2, 1, 0, 0), "CO": (1, 0, 1, 0, 0),
        "H2": (0, 2, 0, 0, 0), "IC8H18": (8, 18, 0, 0, 0), "NC7H16": (7, 16, 0, 0, 0),
    }
    b = [0.0] * 5
    for k, x in X.items():
        for e in range(5):
            b[e] += x * counts[k][e]
    nC, nH, nO, nN, nAr = b
    out = {"N2": nN / 2, "AR": nAr}
    if nO >= 2 * nC + nH / 2:
        out.update(CO2=nC, H2O=nH / 2, O2=(nO - 2 * nC - nH / 2) / 2)
    else:
        K = water_gas_shift_k()
        A = nH / 2 - nO + nC
        D = max(0.0, nO - nC)
        B = K * A + nC + D
        CD = nC * D
        disc = B * B + 4 * (K - 1) * CD
        a = (2 * CD) / (B + math.sqrt(max(0.0, disc))) if CD > 0 else 0.0
        a = min(max(a, 0.0), min(nC, D))
        out.update(CO2=a, CO=nC - a, H2O=max(0.0, D - a), H2=max(0.0, A + a))
    out = {k: v for k, v in out.items() if v > 0}
    s = sum(out.values())
    return {k: v / s for k, v in out.items()}


def unburned_mixture(phi: float, on: float, x_res: float) -> dict[str, float]:
    """(1 - x_res) fresh charge + x_res complete-combustion residual (MOLE fractions)."""
    fc = fresh_charge(phi, on)
    res = complete_combustion_products(fc)
    X: dict[str, float] = {}
    for k, v in fc.items():
        X[k] = X.get(k, 0.0) + (1 - x_res) * v
    for k, v in res.items():
        X[k] = X.get(k, 0.0) + x_res * v
    return X


def to_mech_names(X: dict[str, float], lower: bool = True) -> dict[str, float]:
    """Project species names -> mechanism names (LLNL PRF v2 uses lower case)."""
    return {(k.lower() if lower else k): v for k, v in X.items() if v > 0}


def ignition_delay(gas, T0: float, p0: float, X: dict[str, float], t_max: float = 2.0,
                   rtol: float = 1e-6, atol: float = 1e-15) -> dict:
    """Adiabatic constant-volume ignition delay (Cantera IdealGasMoleReactor).

    Returns dict(tau, tau1, ignited, dT1):
      tau   time of max dT/dt (main / hot ignition), s (t_max if no ignition by t_max);
      tau1  time of the first-stage (cool-flame) dT/dt maximum, s, if a distinct first
            stage exists (a local dT/dt maximum with >= 5 K rise followed by a dip below
            30 % of it before the main event), else = tau;
      dT1   temperature rise at tau1, K (0 if single-stage).
    """
    import cantera as ct

    gas.TPX = T0, p0, X
    r = ct.IdealGasMoleReactor(gas, clone=False)  # Cantera 3.2: share the Solution
    net = ct.ReactorNet([r])
    net.preconditioner = ct.AdaptivePreconditioner()
    # Approximate (preconditioner-only) Jacobian: skipping third-body/falloff derivative terms
    # changes only the Krylov preconditioner, not the solution (tau changes ~4e-5 relative,
    # inside rtol), and saves ~30 % CPU (as in Cantera's preconditioned-integration example).
    net.derivative_settings = {"skip-third-bodies": True, "skip-falloff": True}
    net.rtol = rtol
    net.atol = atol
    ts = [0.0]
    Ts = [T0]
    rate_max = 0.0
    while True:
        t = net.step()
        T = r.T
        dt = t - ts[-1]
        if dt > 0:
            rate = (T - Ts[-1]) / dt
            rate_max = max(rate_max, rate)
        ts.append(t)
        Ts.append(T)
        if T > T0 + 400 and rate < 0.02 * rate_max:
            break
        if t >= t_max:
            break
    n = len(ts)
    rates = [(Ts[i + 1] - Ts[i]) / (ts[i + 1] - ts[i]) for i in range(n - 1)]
    im = max(range(n - 1), key=lambda i: rates[i])
    ignited = Ts[-1] > T0 + 400
    tau = 0.5 * (ts[im] + ts[im + 1]) if ignited else t_max
    # first stage
    tau1 = tau
    dT1 = 0.0
    run_max = 0.0
    i1 = -1
    for i in range(im):
        if rates[i] > run_max:
            run_max = rates[i]
            i1 = i
        elif i1 >= 0 and rates[i] < 0.3 * run_max and Ts[i1 + 1] - T0 > 5.0:
            tau1 = 0.5 * (ts[i1] + ts[i1 + 1])
            dT1 = Ts[i1 + 1] - T0
            break
    return dict(tau=tau, tau1=tau1, ignited=ignited, dT1=dT1)


def ignition_delay_robust(gas, T0: float, p0: float, X: dict[str, float], t_max: float = 2.0) -> dict:
    """ignition_delay with CVODES-failure retries at tighter tolerances (rtol 1e-6 -> 1e-8 ->
    1e-9). Near 550 K the hot-ignition event after seconds of slow chemistry occasionally makes
    CVODES fail its corrector test at rtol 1e-6 (observed at 550 K, 3 bar, PRF 0)."""
    import cantera as ct

    err = None
    for rtol, atol in ((1e-6, 1e-15), (1e-8, 1e-18), (1e-9, 1e-20)):
        try:
            return ignition_delay(gas, T0, p0, X, t_max=t_max, rtol=rtol, atol=atol)
        except ct.CanteraError as e:  # pragma: no cover - solver hiccup
            err = e
    raise err


def file_sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def write_json(path: str, data) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, separators=(",", ":"), allow_nan=False)
        f.write("\n")
    print("wrote", os.path.relpath(path, ROOT))
