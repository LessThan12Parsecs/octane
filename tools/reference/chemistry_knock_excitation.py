"""Excitation (heat-release) time of PRF end-gas autoignition from detailed chemistry.

Excitation time tau_e (Lutz, Kee, Miller, Dwyer & Oppenheim, "Dynamic effects of autoignition
centers for hydrogen and C1,2-hydrocarbon fuels", Proc. Combust. Inst. 22 (1988) 1683):
the time from 5 % of the peak heat-release rate to the peak, in an adiabatic constant-volume
0-D reactor. Used by src/physics/chemistry/knock.ts for the end-gas autoignition burn time.

Mechanism: LLNL gasoline surrogate 2011 (Mehl et al.; the default table mechanism, see
chemistry_mechanisms.py). Mixtures as chemistry_common.py.
Writes test/fixtures/chemistry_knock_excitation.json.
"""
from __future__ import annotations

import os

import cantera as ct
import numpy as np

from chemistry_common import FIXTURES, LLNL_GS_YAML, to_mech_names, unburned_mixture, write_json


def excitation(gas, T0, p0, X, t_max=0.5):
    gas.TPX = T0, p0, X
    r = ct.IdealGasMoleReactor(gas, clone=False)
    net = ct.ReactorNet([r])
    net.preconditioner = ct.AdaptivePreconditioner()
    net.rtol, net.atol = 1e-8, 1e-18
    ts, qs, Ts = [], [], []
    t = 0.0
    while t < t_max:
        t = net.step()
        g = r.phase
        q = -float(np.dot(g.net_production_rates, g.partial_molar_enthalpies))  # W/m3
        ts.append(t)
        qs.append(q)
        Ts.append(g.T)
        if g.T > T0 + 1000 and q < 1e-3 * max(qs):
            break
    ts, qs, Ts = np.array(ts), np.array(qs), np.array(Ts)
    im = int(np.argmax(qs))
    below = np.where(qs[:im] < 0.05 * qs[im])[0]
    i5 = int(below[-1]) if len(below) else 0
    # linear interpolation of the 5 % crossing
    t5 = ts[i5] + (0.05 * qs[im] - qs[i5]) * (ts[i5 + 1] - ts[i5]) / (qs[i5 + 1] - qs[i5])
    return dict(tau=float(ts[im]), tauE=float(ts[im] - t5), qMax=float(qs[im]),
                dTdtMax=float(np.max(np.diff(Ts) / np.diff(ts))))


def main():
    ct.suppress_thermo_warnings()
    gas = ct.Solution(LLNL_GS_YAML)
    out = []
    for T in (850.0, 950.0, 1050.0):
        for p in (30e5, 60e5):
            for on in (90.0, 100.0):
                X = to_mech_names(unburned_mixture(1.0, on, 0.05), lower=False)
                r = excitation(gas, T, p, X)
                r.update(T=T, p=p, on=on, phi=1.0, xres=0.05)
                out.append(r)
                print(f"T={T} p={p/1e5}bar ON={on}: tau={r['tau']*1e3:.3f} ms  tau_e={r['tauE']*1e6:.2f} us")
    te = np.array([o["tauE"] for o in out])
    write_json(os.path.join(FIXTURES, "chemistry_knock_excitation.json"),
               dict(definition="5% of peak heat-release rate to peak (Lutz et al. 1988), CV adiabatic, LLNL gasoline surrogate 2011 (Mehl et al.)",
                    cases=out, tauEMin=float(te.min()), tauEMax=float(te.max()),
                    tauEGeoMean=float(np.exp(np.mean(np.log(te))))))


if __name__ == "__main__":
    main()
