import { describe, expect, it } from 'vitest';
import { EL, N_EQ, NE, NS } from '../core/species';
import { dryAir, freshCharge, fuelFromSelection } from '../thermo/fuels';
import { elementMoles, mixCpMass, mixGamma } from '../thermo/mixture';
import propsFx from '../../../test/fixtures/equilibrium_props.json';
import { EquilibriumSolver, newEqProperties } from './solver';
import { equilibriumProperties } from './properties';

const rel = (a: number, b: number): number => Math.abs(a - b) / Math.abs(b);

function bOf(fuel: 'IC8H18' | 'CH4' | 'C2H5OH', phi: number, egr = 0): Float64Array {
  const X = freshCharge({ fuel: fuelFromSelection({ kind: 'pure', species: fuel }), phi, airX: dryAir(), egrFraction: egr });
  return elementMoles(X, new Float64Array(NE));
}

const STATES: [string, Float64Array, number, number][] = [];
for (const [lbl, b] of [
  ['iC8H18 φ=1', bOf('IC8H18', 1)],
  ['iC8H18 φ=0.7', bOf('IC8H18', 0.7)],
  ['iC8H18 φ=1.4 + 15% EGR', bOf('IC8H18', 1.4, 0.15)],
  ['CH4 φ=1', bOf('CH4', 1)],
  ['C2H5OH φ=2', bOf('C2H5OH', 2)],
] as const) {
  for (const T of [900, 1600, 2300, 2900, 3500]) for (const p of [0.5e5, 8e5, 70e5]) STATES.push([lbl, b, T, p]);
}

describe('equilibriumProperties: analytic derivatives vs finite differences of the solver', () => {
  it('cp, (∂lnV/∂lnT)_p, (∂lnV/∂lnp)_T, cv, γ_s over 75 states (≤ 2e-6 relative)', () => {
    const s = new EquilibriumSolver();
    const out = newEqProperties();
    let wCp = 0;
    let wVT = 0;
    let wVp = 0;
    let wCv = 0;
    let wGs = 0;
    for (const [, b, T, p] of STATES) {
      equilibriumProperties(s, b, T, p, out);
      expect(out.valid).toBe(true);
      const cp = out.cp;
      const vT = out.dlnV_dlnT;
      const vp = out.dlnV_dlnp;
      const cv = out.cv;
      const gs = out.gammaS;
      const m = s.result.mass;
      const S0 = s.result.S;
      const V0 = s.result.V;
      // (∂h/∂T)_p and (∂lnV/∂lnT)_p
      const dT = 1e-4 * T;
      const rp = s.solveTP(b, T + dT, p);
      const hP = rp.H / m;
      const lvP = Math.log(rp.V);
      const rm = s.solveTP(b, T - dT, p);
      const hM = rm.H / m;
      const lvM = Math.log(rm.V);
      wCp = Math.max(wCp, rel((hP - hM) / (2 * dT), cp));
      wVT = Math.max(wVT, rel((lvP - lvM) / Math.log((T + dT) / (T - dT)), vT));
      // (∂lnV/∂lnp)_T
      const ep = 1e-4;
      const lvp1 = Math.log(s.solveTP(b, T, p * (1 + ep)).V);
      const lvp2 = Math.log(s.solveTP(b, T, p * (1 - ep)).V);
      wVp = Math.max(wVp, rel((lvp1 - lvp2) / Math.log((1 + ep) / (1 - ep)), vp));
      // (∂u/∂T)_V via TV equilibria
      const uP = s.solveTV(b, T + dT, V0).U / m;
      const uM = s.solveTV(b, T - dT, V0).U / m;
      wCv = Math.max(wCv, rel((uP - uM) / (2 * dT), cv));
      // γ_s = (∂ln p/∂ln ρ)_s via SP equilibria
      const vs1 = s.solveSP(b, S0, p * (1 + ep), T).V;
      const vs2 = s.solveSP(b, S0, p * (1 - ep), T).V;
      wGs = Math.max(wGs, rel(Math.log((1 + ep) / (1 - ep)) / -Math.log(vs1 / vs2), gs));
    }
    console.log(
      `properties vs FD: cp ${wCp.toExponential(2)}, dlnV/dlnT ${wVT.toExponential(2)}, ` +
        `dlnV/dlnp ${wVp.toExponential(2)}, cv ${wCv.toExponential(2)}, γs ${wGs.toExponential(2)}`,
    );
    expect(wCp).toBeLessThan(2e-6);
    expect(wVT).toBeLessThan(2e-6);
    expect(wVp).toBeLessThan(2e-6);
    expect(wCv).toBeLessThan(2e-6);
    expect(wGs).toBeLessThan(2e-6);
  });

  it('speciesDerivatives (∂lnN_k/∂lnT)_p, (∂lnN_k/∂lnp)_T vs finite differences', () => {
    const s = new EquilibriumSolver();
    const dT = new Float64Array(NS);
    const dP = new Float64Array(NS);
    // Species ≥ 1e-6: h = 1e-5, strict. 1e-10–1e-6: near-stoichiometric O2/CO/H2/NO at ~1e-9 carry
    // ~1e-6 relative noise (their amounts are set by an ill-conditioned direction of the element
    // balance), so they are differenced over h = 2e-3.
    let worst = 0;
    let worstTrace = 0;
    const fd = (b: Float64Array, T: number, p: number, h: number, out: Float64Array[]): void => {
      out[0].set(s.solveTP(b, T * Math.exp(h), p).N);
      out[1].set(s.solveTP(b, T * Math.exp(-h), p).N);
      out[2].set(s.solveTP(b, T, p * Math.exp(h)).N);
      out[3].set(s.solveTP(b, T, p * Math.exp(-h)).N);
    };
    const A = [0, 1, 2, 3].map(() => new Float64Array(NS));
    const B = [0, 1, 2, 3].map(() => new Float64Array(NS));
    for (const [, b, T, p] of STATES) {
      s.solveTP(b, T, p);
      expect(s.speciesDerivatives(dT, dP)).toBe(true);
      const X0 = Float64Array.from(s.result.X);
      fd(b, T, p, 1e-5, A);
      fd(b, T, p, 2e-3, B);
      for (let k = 0; k < N_EQ; k++) {
        const x = X0[k];
        if (!(x > 1e-10)) continue;
        const [F, h] = x > 1e-6 ? [A, 1e-5] : [B, 2e-3];
        const fT = Math.log(F[0][k] / F[1][k]) / (2 * h);
        const fP = Math.log(F[2][k] / F[3][k]) / (2 * h);
        const e = Math.max(Math.abs(fT - dT[k]) / (1 + Math.abs(dT[k])), Math.abs(fP - dP[k]) / (1 + Math.abs(dP[k])));
        if (x > 1e-6) worst = Math.max(worst, e);
        else worstTrace = Math.max(worstTrace, e);
      }
    }
    console.log(`species derivatives vs FD: ${worst.toExponential(2)} (X > 1e-6), ${worstTrace.toExponential(2)} (1e-10 < X ≤ 1e-6)`);
    expect(worst).toBeLessThan(1e-5);
    expect(worstTrace).toBeLessThan(2e-3);
  });

  it('frozen limit at low T: cp → cp_frozen, dlnV/dlnT → 1, dlnV/dlnp → −1, γs → γ', () => {
    const s = new EquilibriumSolver();
    const out = newEqProperties();
    equilibriumProperties(s, bOf('IC8H18', 0.8), 500, 5e5, out);
    expect(out.valid).toBe(true);
    const cpF = mixCpMass(out.X, 500);
    expect(rel(out.cpFrozen, cpF)).toBeLessThan(1e-12);
    // not exactly frozen: NO/CO/H2 at ~1e-10 still shift with T (≈ 3e-8 of cp at 500 K)
    expect(rel(out.cp, cpF)).toBeLessThan(1e-6);
    expect(Math.abs(out.dlnV_dlnT - 1)).toBeLessThan(1e-12);
    expect(Math.abs(out.dlnV_dlnp + 1)).toBeLessThan(1e-12);
    expect(rel(out.gammaS, mixGamma(out.X, 500))).toBeLessThan(1e-6);
    expect(rel(out.soundSpeed ** 2, (out.gammaS * out.p) / out.rho)).toBeLessThan(1e-12);
  });

  it('Maxwell relations: (∂h/∂ln p)_T = p v (1 − (∂lnV/∂lnT)_p), (∂s/∂ln p)_T = −(p v/T)(∂lnV/∂lnT)_p', () => {
    // Exact thermodynamic identities (dh = T ds + v dp with (∂s/∂p)_T = −(∂v/∂T)_p). They tie the
    // analytic (∂lnV/∂lnT)_p to the pressure dependence of the equilibrium h and s, which the
    // two-zone model needs for the burned zone (h depends on p through dissociation).
    const s = new EquilibriumSolver();
    const out = newEqProperties();
    let wH = 0;
    let wS = 0;
    for (const [, b, T, p] of STATES) {
      if (T < 1600) continue; // below this the p-dependence is ~0 and the check is vacuous
      equilibriumProperties(s, b, T, p, out);
      const pv = out.p / out.rho;
      const m = s.result.mass;
      const e = 1e-4;
      const r1 = s.solveTP(b, T, p * Math.exp(e));
      const h1 = r1.H / m;
      const s1 = r1.S / m;
      const r2 = s.solveTP(b, T, p * Math.exp(-e));
      const dh = (h1 - r2.H / m) / (2 * e);
      const ds = (s1 - r2.S / m) / (2 * e);
      const dhRef = pv * (1 - out.dlnV_dlnT);
      const dsRef = (-pv / T) * out.dlnV_dlnT;
      wH = Math.max(wH, Math.abs(dh - dhRef) / pv);
      wS = Math.max(wS, Math.abs(ds - dsRef) / (pv / T));
    }
    console.log(`Maxwell: (∂h/∂lnp)_T ${wH.toExponential(2)} (×pv), (∂s/∂lnp)_T ${wS.toExponential(2)} (×pv/T)`);
    expect(wH).toBeLessThan(1e-6);
    expect(wS).toBeLessThan(1e-6);
  });

  it('equilibrium shifting raises cp far above frozen at high T (dissociation)', () => {
    const s = new EquilibriumSolver();
    const out = equilibriumProperties(s, bOf('IC8H18', 1), 3000, 1e5);
    expect(out.cp / out.cpFrozen).toBeGreaterThan(2);
    expect(out.dlnV_dlnT).toBeGreaterThan(1.1);
    expect(out.dlnV_dlnp).toBeLessThan(-1.005);
  });

  it('invalid (not thrown) for infeasible element sets', () => {
    const s = new EquilibriumSolver();
    const b = new Float64Array(NE);
    b[EL.C] = 1;
    b[EL.O] = 0.5;
    const out = equilibriumProperties(s, b, 2000, 1e5);
    expect(out.valid).toBe(false);
  });
});

describe('equilibriumProperties vs Cantera (equilibrium_props.json)', () => {
  it('h, u, s, M (1e-10) and FD cp, (∂lnV/∂lnT)_p, (∂lnV/∂lnp)_T (≤ 1e-6)', () => {
    const s = new EquilibriumSolver();
    const out = newEqProperties();
    let wH = 0;
    let wS = 0;
    let wM = 0;
    let wCp = 0;
    let wVT = 0;
    let wVp = 0;
    for (const c of propsFx.cases) {
      equilibriumProperties(s, Float64Array.from(c.b), c.T, c.p, out);
      expect(out.valid).toBe(true);
      const scaleE = out.cp * c.T; // energy scale for absolute-convention h, u
      wH = Math.max(wH, Math.abs(out.h - c.h) / scaleE, Math.abs(out.u - c.u) / scaleE);
      wS = Math.max(wS, Math.abs(out.s - c.s) / out.cp);
      wM = Math.max(wM, rel(out.M, c.M));
      wCp = Math.max(wCp, rel(out.cp, c.cp_fd));
      wVT = Math.max(wVT, rel(out.dlnV_dlnT, c.dlnV_dlnT_fd));
      wVp = Math.max(wVp, rel(out.dlnV_dlnp, c.dlnV_dlnp_fd));
    }
    console.log(
      `props vs Cantera (${propsFx.cases.length}): h,u ${wH.toExponential(2)} (×cpT), s ${wS.toExponential(2)} (×cp), ` +
        `M ${wM.toExponential(2)}, cp ${wCp.toExponential(2)}, dlnV/dlnT ${wVT.toExponential(2)}, dlnV/dlnp ${wVp.toExponential(2)}`,
    );
    expect(wH).toBeLessThan(1e-10);
    expect(wS).toBeLessThan(1e-10);
    expect(wM).toBeLessThan(1e-10);
    expect(wCp).toBeLessThan(1e-6);
    expect(wVT).toBeLessThan(1e-6);
    expect(wVp).toBeLessThan(1e-6);
  });

  it('solver.properties() after a UV solve equals equilibriumProperties at the resulting (T, p)', () => {
    const s = new EquilibriumSolver();
    const b = bOf('IC8H18', 1.1, 0.1);
    const r = s.solveTP(b, 2600, 50e5);
    const U = r.U;
    const V = r.V * 1.3;
    s.solveUV(b, U, V, 2500);
    const a = s.properties();
    const s2 = new EquilibriumSolver();
    const c = equilibriumProperties(s2, b, a.T, a.p);
    for (const key of ['h', 'u', 's', 'cp', 'cv', 'dlnV_dlnT', 'dlnV_dlnp', 'gammaS', 'M'] as const) {
      expect(rel(a[key], c[key]), key).toBeLessThan(1e-9);
    }
    expect(a.u * s.result.mass).toBeCloseTo(U, 3);
  });
});
