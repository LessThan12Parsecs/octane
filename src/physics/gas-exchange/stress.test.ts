/**
 * Seeded random sweeps over the gas-exchange hot paths: no NaN/Inf, physical signs and bounds,
 * exact pressure derivatives, and continuity (no jumps that would upset an ODE integrator) over
 * extreme inputs (1 kPa–300 bar, 200–3500 K, γ 1.1–1.7, areas 1e-9–1e-2 m², band widths
 * 1e-8–0.1, arbitrary valve geometry / shroud / lash / event timing).
 */
import { describe, expect, it } from 'vitest';
import { NS, SP } from '../core/species';
import { criticalPressureRatio, newOrificeFlow, orificeFlow, regularizedFluxFunction } from './orifice';
import { gasStateFromNU, gasStateFromTPX, newGasState } from './plenum';
import { throttleGeometricArea } from './throttle';
import { ValveFlowModel } from './valve-flow';
import { ValveLiftProfile } from './valve-lift';

/** Deterministic LCG in [0, 1). */
function rng(seed: number): () => number {
  let s = seed;
  return () => {
    s = (s * 1103515245 + 12345) % 2147483648;
    return s / 2147483648;
  };
}

describe('seeded stress sweeps (gas exchange)', () => {
  it('orificeFlow: finite, signed like Δp, monotone in p_up/p_down, exact derivatives', () => {
    const r = rng(1);
    const lu = (a: number, b: number) => Math.exp(Math.log(a) + r() * Math.log(b / a));
    const out = newOrificeFlow();
    const o2 = newOrificeFlow();
    let worstDeriv = 0;
    let bad = 0;
    for (let i = 0; i < 20000; i++) {
      const pa = lu(1e3, 3e7);
      const pb = r() < 0.3 ? pa * (1 + (r() - 0.5) * 2e-3) : lu(1e3, 3e7);
      const Ta = lu(200, 3500);
      const Tb = lu(200, 3500);
      const ga = 1.1 + 0.6 * r();
      const gb = 1.1 + 0.6 * r();
      const A = lu(1e-9, 1e-2);
      const d = r() < 0.5 ? 1e-3 : lu(1e-8, 0.1);
      const m = orificeFlow(A, pa, Ta, 287, ga, pb, Tb, 290, gb, out, d);
      for (const v of [out.mdot, out.jetVelocity, out.dmdotdpa, out.dmdotdpb, out.pressureDrop]) {
        if (!Number.isFinite(v)) bad++;
      }
      if (m !== 0 && Math.sign(m) !== Math.sign(pa - pb)) bad++;
      if (out.dmdotdpa < 0 || out.dmdotdpb > 0) bad++;
      // bounded by the choked flow of the upstream side
      const up = pa >= pb;
      const g = up ? ga : gb;
      const psiMax = Math.sqrt(g) * Math.pow(2 / (g + 1), (g + 1) / (2 * (g - 1)));
      if (Math.abs(m) > ((A * Math.max(pa, pb) * psiMax) / Math.sqrt(up ? 287 * Ta : 290 * Tb)) * (1 + 1e-12)) bad++;
      // ∂ṁ/∂p_a vs central differences (same side upstream on both sides of the stencil)
      const h = 1e-7 * Math.max(pa, pb);
      if (pa + h > pb === pa - h > pb) {
        const fd = (orificeFlow(A, pa + h, Ta, 287, ga, pb, Tb, 290, gb, o2, d) - orificeFlow(A, pa - h, Ta, 287, ga, pb, Tb, 290, gb, o2, d)) / (2 * h);
        if (Math.abs(fd) > 1e-25) worstDeriv = Math.max(worstDeriv, Math.abs(fd - out.dmdotdpa) / Math.abs(fd));
      }
    }
    expect(bad).toBe(0);
    expect(worstDeriv).toBeLessThan(1e-3); // FD truncation/round-off at the band edge and choke
  });

  it('regularised flux: continuous at the band edge and at the choke for any γ, δ', () => {
    const r = rng(2);
    for (let i = 0; i < 2000; i++) {
      const g = 1.05 + 0.65 * r();
      const d = Math.exp(Math.log(1e-8) + r() * Math.log(0.05 / 1e-8));
      const e = 1e-12;
      for (const x of [d, 1 - criticalPressureRatio(g)]) {
        const a = regularizedFluxFunction(x * (1 - e), g, d);
        const b = regularizedFluxFunction(x * (1 + e), g, d);
        expect(Math.abs(a - b) / b).toBeLessThan(1e-10);
      }
    }
  });

  it('valve effective area: finite, ≥ 0, continuous in lift for random geometry/shroud/port', () => {
    const r = rng(3);
    const lu = (a: number, b: number) => Math.exp(Math.log(a) + r() * Math.log(b / a));
    let bad = 0;
    let worstSlope = 0;
    for (let i = 0; i < 300; i++) {
      const Di = lu(0.005, 0.08);
      const Ds = Di * 0.5 * r();
      const spec = {
        count: 1 + Math.floor(3 * r()),
        headDiameter: Di * (1 + lu(1e-4, 0.3)),
        seatInnerDiameter: Di,
        seatAngle: ((5 + 80 * r()) * Math.PI) / 180,
        stemDiameter: Ds,
        maxLift: 0.3 * Di,
        openDeg: 0,
        closeDeg: 100,
        timingLiftThreshold: 0,
        position: [0.01, 0.002] as [number, number],
        shroudArcDeg: r() < 0.5 ? 0 : 360 * r(),
        shroudDirection: 6 * r(),
      };
      const m = new ValveFlowModel(spec, r() < 0.5 ? 'intake' : 'exhaust', {
        portDiameter: Math.max(Ds * 1.01, Di * (0.6 + 0.5 * r())),
        shroudHeight: r() < 0.5 ? undefined : 0.3 * Di * r(),
      });
      // Lipschitz bound: 4 × the curtain-area slope (the effective-area slope reaches ≈ 2.5× it in
      // the shrouded, port-limited regime; a jump would not shrink with the step)
      const slopeMax = 4 * Math.PI * Di * spec.count;
      for (const rev of [false, true]) {
        const N = 2000;
        const dL = (0.6 * Di) / N;
        let prev = 0;
        for (let k = 0; k <= N; k++) {
          const L = k * dL;
          const a = m.effectiveArea(L, rev);
          if (!(Number.isFinite(a) && a >= 0 && Number.isFinite(m.dischargeCoefficient(L, rev)))) bad++;
          if (k > 0) worstSlope = Math.max(worstSlope, Math.abs(a - prev) / (slopeMax * dL));
          prev = a;
        }
      }
    }
    expect(bad).toBe(0);
    expect(worstSlope).toBeLessThan(1);
  });

  it('valve lift: 0 ≤ L ≤ L_max, periodic, Lipschitz, timing honoured, with and without lash', () => {
    const r = rng(4);
    const lu = (a: number, b: number) => Math.exp(Math.log(a) + r() * Math.log(b / a));
    let built = 0;
    let bad = 0;
    let worstTiming = 0;
    for (let i = 0; i < 600; i++) {
      const maxLift = lu(1e-4, 0.02);
      const open = -360 + 720 * r();
      const dur = 20 + 400 * r();
      const thr = r() < 0.5 ? 0 : 0.2 * maxLift * r();
      const lash = r() < 0.5 ? 0 : 0.1 * maxLift * r();
      let p: ValveLiftProfile;
      try {
        p = new ValveLiftProfile({ maxLift, openDeg: open, closeDeg: open + dur, timingLiftThreshold: thr, lash });
      } catch {
        continue; // cam longer than 720° (threshold near the nose): rejected by design
      }
      built++;
      worstTiming = Math.max(worstTiming, Math.abs(p.lift(open) - thr) / maxLift, Math.abs(p.lift(open + dur) - thr) / maxLift);
      const vMax = (4 * p.camMaxLift) / p.halfDurationDeg; // |s'| ≤ 4 for the default cam
      let prev = p.lift(-800);
      for (let th = -800; th < 800; th += 0.37) {
        const l = p.lift(th);
        if (!(Number.isFinite(l) && l >= 0 && l <= maxLift * (1 + 1e-12))) bad++;
        if (Math.abs(l - prev) > vMax * 0.37 + 1e-15) bad++;
        if (Math.abs(p.lift(th + 720) - l) > 1e-12 * maxLift) bad++;
        prev = l;
      }
    }
    expect(bad).toBe(0);
    expect(worstTiming).toBeLessThan(1e-12);
    expect(built).toBeGreaterThan(400);
  });

  it('throttle area within [0, bore area] for any angle/shaft/closed angle', () => {
    const r = rng(5);
    const A0 = 0.25 * Math.PI * 0.05 * 0.05;
    let bad = 0;
    for (let i = 0; i < 50000; i++) {
      const A = throttleGeometricArea(0.05, 2 * r(), r() < 0.1 ? 0 : 0.999 * r(), 1.5 * r());
      if (!(Number.isFinite(A) && A >= 0 && A <= A0 * (1 + 1e-12))) bad++;
    }
    expect(bad).toBe(0);
  });

  it('plenum state recovery round-trips random compositions at 200–5000 K, 100 Pa–500 bar', () => {
    const r = rng(6);
    const s = newGasState();
    const s2 = newGasState();
    const X = new Float64Array(NS);
    const N = new Float64Array(NS);
    let worst = 0;
    for (let i = 0; i < 2000; i++) {
      X.fill(0);
      for (let k = 0; k < NS; k++) if (r() < 0.4) X[k] = r();
      X[SP.N2] += 0.1;
      const T = Math.exp(Math.log(200) + r() * Math.log(5000 / 200));
      const p = Math.exp(Math.log(100) + r() * Math.log(5e7 / 100));
      const V = 1e-3;
      gasStateFromTPX(T, p, X, s);
      const n = (p * V) / (8.31446261815324 * T);
      for (let k = 0; k < NS; k++) N[k] = n * s.X[k];
      const U = s.u * s.rho * V;
      gasStateFromNU(N, U, V, 100 + 6000 * r(), s2); // arbitrary (bad) warm start
      expect(Number.isFinite(s2.T) && Number.isFinite(s2.p)).toBe(true);
      worst = Math.max(worst, Math.abs(s2.T - T) / T, Math.abs(s2.p - p) / p);
    }
    expect(worst).toBeLessThan(1e-9);
  });
});
