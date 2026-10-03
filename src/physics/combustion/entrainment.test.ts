import { describe, expect, it } from 'vitest';
import {
  effectiveTurbulentBurningVelocity,
  enflamedVolume,
  entrainmentRates,
  keckCharacteristicLength,
  keckCharacteristicSpeed,
  keckMeanInletSpeed,
  kolmogorovScale,
  newEntrainmentInputs,
  newEntrainmentRates,
  steadyFlameBrushThickness,
  taylorMicroscale,
  type EntrainmentInputs,
} from './entrainment';
import { FlameGeometry } from './flame-geometry';
import { integralLengthScale, rapidDistortionIntensity } from './turbulence';

const PI = Math.PI;

/** RK4 on (m_e, m_b) with inputs refreshed by `set(me, mb, inp)` before each stage. */
function burn(
  inp: EntrainmentInputs,
  me0: number,
  mb0: number,
  dt: number,
  n: number,
  set: (me: number, mb: number, inp: EntrainmentInputs) => void,
  each?: (t: number, me: number, mb: number) => void,
): [number, number] {
  const r = newEntrainmentRates();
  let me = me0;
  let mb = mb0;
  const f = (e: number, b: number, o: number[]) => {
    set(e, b, inp);
    inp.me = e;
    inp.mb = b;
    entrainmentRates(inp, r);
    o[0] = r.dme;
    o[1] = r.dmb;
  };
  const k1 = [0, 0];
  const k2 = [0, 0];
  const k3 = [0, 0];
  const k4 = [0, 0];
  for (let i = 0; i < n; i++) {
    f(me, mb, k1);
    f(me + (dt / 2) * k1[0], mb + (dt / 2) * k1[1], k2);
    f(me + (dt / 2) * k2[0], mb + (dt / 2) * k2[1], k3);
    f(me + dt * k3[0], mb + dt * k3[1], k4);
    me += (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    mb += (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    each?.((i + 1) * dt, me, mb);
  }
  return [me, mb];
}

describe('turbulence length scales used by the burn-up model', () => {
  it('Taylor microscale: λ² = 15νu′²/ε with ε = A u′³/L, i.e. λ/L = (15/A)^½ Re_L^−½', () => {
    const L = 0.005;
    const u = 1.2;
    const nu = 3.7e-6;
    const lam = taylorMicroscale(L, u, nu);
    const eps = (u * u * u) / L;
    expect((lam * lam * eps) / (15 * nu * u * u)).toBeCloseTo(1, 13);
    expect(lam / L).toBeCloseTo(Math.sqrt(15) * Math.pow((u * L) / nu, -0.5), 13);
    expect(taylorMicroscale(L, u, nu, 0.5) / lam).toBeCloseTo(Math.SQRT2, 13);
    const eta = kolmogorovScale(L, u, nu);
    expect(eta).toBeCloseTo(Math.pow((nu * nu * nu) / eps, 0.25), 15);
    // engine conditions: η < λ < L (λ ≈ 0.5 mm, η ≈ 20 µm)
    expect(eta).toBeLessThan(lam);
    expect(lam).toBeLessThan(L);
    expect(taylorMicroscale(L, 0, nu)).toBe(Infinity);
  });
});

describe("model length/speed scales vs Keck's (1982) empirical correlations at CFR conditions", () => {
  it('λ from the K–k length scale is the same order as Keck ℓ_T; u′ vs Keck u_T (reported)', () => {
    const IN = 0.0254;
    const B = 3.25 * IN;
    const h = (4.5 * IN) / 6; // TDC, CR 7
    const rhoI = 101325 / (287 * 325);
    const rhoU = rhoI * Math.pow(7, 1 / 1.35); // isentropic compression to TDC
    const Tu = 325 * Math.pow(7, 0.35);
    // air-like viscosity μ ≈ 3.25e-5 Pa s (T/650 K)^0.7 — order of magnitude only (transport.ts owns μ)
    const nu = (3.25e-5 * Math.pow(Tu / 650, 0.7)) / rhoU;
    const uP = 1.15; // K–k model u′(TDC) at 600 rpm, CR 7 (turbulence.test.ts)
    const lam = taylorMicroscale(integralLengthScale(h, B), uP, nu);
    const lT = keckCharacteristicLength(0.246 * IN, rhoU, rhoI);
    const uI = keckMeanInletSpeed(0.9, (PI * B * B) / 4, PI * 1.2 * IN * 0.246 * IN, 600, 4.5 * IN);
    const uT = keckCharacteristicSpeed(uI, rhoU, rhoI);
    console.log(
      `[entrainment] CFR TDC (600 rpm, CR 7): λ = ${(lam * 1e3).toFixed(2)} mm vs Keck ℓ_T = ${(lT * 1e3).toFixed(2)} mm; ` +
        `u′ = ${uP} m/s vs Keck u_T = ${uT.toFixed(2)} m/s (ū_i = ${uI.toFixed(1)} m/s)`,
    );
    expect(lam / lT).toBeGreaterThan(1 / 4);
    expect(lam / lT).toBeLessThan(4);
  });
});

describe('Keck (1982) entrainment / burn-up equations — limiting cases', () => {
  it('laminar limit u′ → 0: dm_e/dt = dm_b/dt = ρ_u A_f S_L, s_b = S_L (Keck eq. 4.6)', () => {
    const inp = newEntrainmentInputs();
    Object.assign(inp, { rhoU: 8, frontArea: 1e-3, uPrime: 0, SL: 0.6, me: 1e-6, mb: 1e-6, lambda: Infinity });
    const r = entrainmentRates(inp, newEntrainmentRates());
    expect(r.dme).toBeCloseTo(8 * 1e-3 * 0.6, 15);
    expect(r.dmb).toBeCloseTo(8 * 1e-3 * 0.6, 15);
    expect(r.burningSpeed).toBeCloseTo(0.6, 15);
    expect(effectiveTurbulentBurningVelocity(inp)).toBeCloseTo(0.6, 15);
  });

  it('free laminar spherical flame expands at dr/dt = (ρ_u/ρ_b) S_L', () => {
    const rhoU = 8;
    const rhoB = 1.3;
    const SL = 0.5;
    const r0 = 1e-3;
    const inp = newEntrainmentInputs();
    Object.assign(inp, { rhoU, uPrime: 0, SL, lambda: Infinity });
    const m0 = (rhoB * 4 * PI * r0 ** 3) / 3;
    const radius = (me: number, mb: number) => Math.cbrt((3 * enflamedVolume(me, mb, rhoB, rhoU)) / (4 * PI));
    const T = 2e-3;
    const [me, mb] = burn(inp, m0, m0, T / 2000, 2000, (e, b, i) => {
      i.frontArea = 4 * PI * radius(e, b) ** 2;
    });
    expect(me).toBeCloseTo(mb, 15);
    expect(radius(me, mb) / (r0 + (rhoU / rhoB) * SL * T)).toBeCloseTo(1, 10);
  });

  it('constant area: μ(t) = ρ_u A u′τ_b(1 − e^{−t/τ_b}); s_b → u′ + S_L; brush thickness → u′λ/S_L (Keck 4.3, 4.7)', () => {
    const inp = newEntrainmentInputs();
    const rhoU = 6;
    const A = 2e-3;
    const u = 1.2;
    const SL = 0.7;
    const lam = 1e-3;
    Object.assign(inp, { rhoU, frontArea: A, uPrime: u, SL, lambda: lam });
    const tauB = lam / SL;
    const T = 8 * tauB;
    const [me, mb] = burn(inp, 0, 0, T / 4000, 4000, () => {});
    const muExact = rhoU * A * u * tauB * (1 - Math.exp(-T / tauB));
    expect((me - mb) / muExact).toBeCloseTo(1, 10);
    const r = entrainmentRates({ ...inp, me, mb }, newEntrainmentRates());
    expect(r.burningSpeed / (u + SL)).toBeCloseTo(1, 3); // (1 − e⁻⁸) u′/(u′ + S_L)
    const brush = (me - mb) / (rhoU * A); // (V_f − V_b)/A_f
    expect(brush / steadyFlameBrushThickness(u, lam, SL)).toBeCloseTo(1, 3);
    expect(steadyFlameBrushThickness(u, lam, SL)).toBeCloseTo((u * lam) / SL, 15);
  });

  it('initial burning: s_b/S_L = 1 + t/τ_T, τ_T = λ/u′ (Keck eq. 4.8, constant area)', () => {
    const inp = newEntrainmentInputs();
    Object.assign(inp, { rhoU: 6, frontArea: 2e-3, uPrime: 1.2, SL: 0.7, lambda: 1e-3 });
    const t = 1e-5; // ≪ τ_b = 1.43 ms
    const [me, mb] = burn(inp, 0, 0, t / 100, 100, () => {});
    const sb = effectiveTurbulentBurningVelocity({ ...inp, me, mb });
    expect((sb / 0.7 - 1) / (t / (1e-3 / 1.2))).toBeCloseTo(1, 2);
  });

  it('entrainment stops at the total mass; burn-up of the remaining pocket continues', () => {
    const inp = newEntrainmentInputs();
    Object.assign(inp, { rhoU: 6, frontArea: 1e-3, uPrime: 1.2, SL: 0.7, lambda: 1e-3, me: 1e-4, mb: 5e-5, mTotal: 1e-4 });
    const r = entrainmentRates(inp, newEntrainmentRates());
    expect(r.dme).toBe(0);
    expect(r.dmb).toBeCloseTo(6 * 1e-3 * 0.7 + 5e-5 / (1e-3 / 0.7), 15);
  });
});

describe('contract and degenerate inputs (review additions)', () => {
  it('DESIGN.md input shape {rhoU, frontArea, uPrime, SL, me, mb, lambda} works without mTotal', () => {
    const r = entrainmentRates(
      { rhoU: 6, frontArea: 1e-3, uPrime: 1.2, SL: 0.7, me: 2e-5, mb: 1e-5, lambda: 1e-3 },
      newEntrainmentRates(),
    );
    expect(r.dme).toBeCloseTo(6 * 1e-3 * 1.9, 15);
    expect(r.dmb).toBeCloseTo(6 * 1e-3 * 0.7 + 1e-5 / (1e-3 / 0.7), 15);
  });

  it('never returns NaN/Inf rates for zero or invalid scales', () => {
    const r = newEntrainmentRates();
    const cases: Array<Partial<EntrainmentInputs>> = [
      { SL: 0, uPrime: 0, lambda: Infinity },
      { SL: 0, uPrime: 1, lambda: 1e-3 },
      { SL: 0.5, uPrime: 1, lambda: 0 },
      { SL: 0.5, uPrime: 1, lambda: NaN },
      { frontArea: 0, me: 1e-5, mb: 0 },
      { rhoU: 0 },
      { me: 0, mb: 1e-5 }, // m_b > m_e (integrator round-off) ⇒ no negative burn-up
    ];
    for (const c of cases) {
      const inp = { ...newEntrainmentInputs(), rhoU: 6, frontArea: 1e-3, uPrime: 1, SL: 0.5, me: 2e-5, mb: 1e-5, ...c };
      entrainmentRates(inp, r);
      expect(Number.isFinite(r.dme) && Number.isFinite(r.dmb) && Number.isFinite(r.burningSpeed)).toBe(true);
      expect(r.dmb).toBeGreaterThanOrEqual(0);
      expect(r.dme).toBeGreaterThanOrEqual(0);
      expect(Number.isFinite(effectiveTurbulentBurningVelocity(inp))).toBe(true);
    }
  });
});

describe('closed-vessel burn in the CFR chamber at TDC (two-zone ideal gas + exact flame geometry)', () => {
  // CFR 3.25 in bore, TDC clearance height at CR 7 (4.5 in stroke): h = S/(CR − 1) = 19.05 mm.
  const IN = 0.0254;
  const B = 3.25 * IN;
  const R = B / 2;
  const h = (4.5 * IN) / 6;
  const V = PI * R * R * h;
  const gamma = 1.3;
  const p0 = 15e5;
  const Tu0 = 650;
  const Rg = 287;
  const q = 2.75e6; // J/kg of stoichiometric iso-octane–air mixture (LHV 44.3 MJ/kg ÷ (1 + 15.1))
  const rhoU0 = p0 / (Rg * Tu0);
  const m = rhoU0 * V;

  function run(center: [number, number, number], lambdaOf: (rhoU: number) => number, u0: number, SL: number) {
    const g = new FlameGeometry(B, center, { maxHeight: 0.03 });
    const inp = newEntrainmentInputs();
    inp.mTotal = m;
    const state = (me: number, mb: number) => {
      const p = p0 + ((gamma - 1) * mb * q) / V; // constant-volume energy balance
      const rhoU = rhoU0 * Math.pow(p / p0, 1 / gamma); // isentropic unburned gas
      const Vf = Math.min(V - (m - me) / rhoU, V); // burned + entrained-unburned volume
      return { p, rhoU, Vf };
    };
    const set = (me: number, mb: number, i: EntrainmentInputs) => {
      const s = state(me, mb);
      i.rhoU = s.rhoU;
      i.uPrime = rapidDistortionIntensity(u0, rhoU0, s.rhoU); // RDT of the unburned gas
      i.SL = SL;
      i.lambda = lambdaOf(s.rhoU);
      const rf = g.radiusForVolume(Math.max(s.Vf, 0), h);
      i.frontArea = g.frontArea(rf, h);
    };
    // 1 mm burned kernel
    const rhoB0 = (rhoU0 * Tu0) / (Tu0 + q / (Rg / (gamma - 1)));
    const mk = rhoB0 * (4 / 3) * PI * 1e-9;
    let t10 = NaN;
    let t90 = NaN;
    let ok = true;
    const dt = 2e-6;
    const [me, mb] = burn(inp, mk, mk, dt, 10000, set, (t, e, b) => {
      // (RK4 may overshoot m by ~1e-7 in the step where the front reaches the far corner)
      if (!(b <= e * (1 + 1e-12) && e <= m * (1 + 1e-6) && b >= 0)) ok = false;
      if (Number.isNaN(t10) && b >= 0.1 * m) t10 = t;
      if (Number.isNaN(t90) && b >= 0.9 * m) t90 = t;
    });
    return { me, mb, t10, t90, ok };
  }

  it('mass conservation m_b ≤ m_e ≤ m and burn-out; central ignition burns 10–90 % in a few ms', () => {
    const lam = () => 1e-3; // task: u′ ≈ 1.2 m/s, S_L ≈ 0.5–1 m/s, λ ≈ 1 mm
    const c = run([0, -h / 2, 0], lam, 1.2, 0.7);
    expect(c.ok).toBe(true);
    expect(c.mb / m).toBeGreaterThan(0.999);
    expect(c.t90 - c.t10).toBeGreaterThan(1e-3);
    expect(c.t90 - c.t10).toBeLessThan(10e-3);
    const s = run([0.9 * R, -0.002, 0], lam, 1.2, 0.7); // CFR side-mounted plug (draft position)
    expect(s.ok).toBe(true);
    expect(s.t90 - s.t10).toBeGreaterThan(c.t90 - c.t10); // longer flame travel from the side
    const sp = run([0, -0.002, R - 0.001], lam, 1.2, 0.7); // CFR_F1 spec gap centre (engines/cfr.ts), 1 mm from the liner
    expect(sp.ok).toBe(true);
    expect(sp.t90).toBeLessThan(20e-3); // slowest case: 1 mm from the wall the kernel is half-quenched by the liner
    expect(sp.t90 - sp.t10).toBeGreaterThan(s.t90 - s.t10);
    // laminar-only reference is much slower
    const l = run([0, -h / 2, 0], () => Infinity, 0, 0.7);
    expect(Number.isNaN(l.t90) || l.t90 - l.t10 > 2 * (c.t90 - c.t10)).toBe(true);
    console.log(
      `[entrainment] closed-vessel CFR TDC burn (u′₀ 1.2 m/s, S_L 0.7 m/s, λ 1 mm): 0–10 % ${(c.t10 * 1e3).toFixed(2)} ms, ` +
        `10–90 % ${((c.t90 - c.t10) * 1e3).toFixed(2)} ms central; ${((s.t90 - s.t10) * 1e3).toFixed(2)} ms side plug (0.9R), ` +
        `${((sp.t90 - sp.t10) * 1e3).toFixed(2)} ms CFR_F1 spec plug (R − 1 mm, 0–10 % ${(sp.t10 * 1e3).toFixed(2)} ms); ` +
        `laminar only: ${(l.mb / m * 100).toFixed(1)} % burned after 20 ms`,
    );
  });
});
