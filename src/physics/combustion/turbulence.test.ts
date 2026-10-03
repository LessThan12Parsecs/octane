import { describe, expect, it } from 'vitest';
import {
  DEFAULT_TURBULENCE_PARAMS,
  KK_PRODUCTION_CONSTANT,
  angularMomentumLengthScale,
  initTurbulenceState,
  intakeJetVelocity,
  integralLengthScale,
  meanFlowVelocity,
  newTurbulenceInputs,
  newTurbulenceRates,
  newTurbulenceState,
  rapidDistortionIntensity,
  swirlAngularVelocity,
  swirlTorqueFromJet,
  turbulenceDerivatives,
  turbulenceIntensity,
  type TurbulenceInputs,
  type TurbulenceParams,
  type TurbulenceState,
} from './turbulence';

/** Classical RK4 on the turbulence state with inputs refreshed by `setInputs(t)`. */
function integrate(
  st: TurbulenceState,
  inp: TurbulenceInputs,
  params: TurbulenceParams,
  t0: number,
  t1: number,
  n: number,
  setInputs: (t: number, inp: TurbulenceInputs) => void,
): void {
  const r = newTurbulenceRates();
  const s = newTurbulenceState();
  const dt = (t1 - t0) / n;
  const k1 = [0, 0, 0];
  const k2 = [0, 0, 0];
  const k3 = [0, 0, 0];
  const k4 = [0, 0, 0];
  const f = (t: number, K: number, k: number, H: number, o: number[]) => {
    setInputs(t, inp);
    s.K = K;
    s.k = k;
    s.swirl = H;
    turbulenceDerivatives(s, inp, r, params);
    o[0] = r.dK;
    o[1] = r.dk;
    o[2] = r.dSwirl;
  };
  let t = t0;
  for (let i = 0; i < n; i++) {
    f(t, st.K, st.k, st.swirl, k1);
    f(t + dt / 2, st.K + (dt / 2) * k1[0], st.k + (dt / 2) * k1[1], st.swirl + (dt / 2) * k1[2], k2);
    f(t + dt / 2, st.K + (dt / 2) * k2[0], st.k + (dt / 2) * k2[1], st.swirl + (dt / 2) * k2[2], k3);
    f(t + dt, st.K + dt * k3[0], st.k + dt * k3[1], st.swirl + dt * k3[2], k4);
    st.K += (dt / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
    st.k += (dt / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
    st.swirl += (dt / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
    t += dt;
  }
}

describe('K–k model building blocks', () => {
  it('production constant = 4.5 C_μ √(2/3) = 0.3307 (Xu & Filipi 2020 eq. 37)', () => {
    expect(KK_PRODUCTION_CONSTANT).toBeCloseTo(0.3307, 4);
  });

  it('u′ and U from the extensive energies', () => {
    expect(turbulenceIntensity(1.5 * 2 * 3 * 3, 2)).toBeCloseTo(3, 14);
    expect(meanFlowVelocity(0.5 * 2 * 7 * 7, 2)).toBeCloseTo(7, 14);
    const st = initTurbulenceState(0.8e-3, 1.2, 4, 50, 0.08255, newTurbulenceState());
    expect(turbulenceIntensity(st.k, 0.8e-3)).toBeCloseTo(1.2, 14);
    expect(meanFlowVelocity(st.K, 0.8e-3)).toBeCloseTo(4, 14);
    expect(swirlAngularVelocity(st.swirl, 0.8e-3, 0.08255)).toBeCloseTo(50, 12);
  });

  it('energy bookkeeping: dK + dk = ½ṁv² − mε − (K + k)ṁ_out/m + (2/3)k dlnρ/dt', () => {
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m: 7e-4, mDotIn: 0.02, mDotOut: 0.003, vIn: 40, L: 0.01, dlnRhoDt: 35 });
    const st: TurbulenceState = { K: 0.3, k: 0.05, swirl: 0 };
    const r = turbulenceDerivatives(st, inp, newTurbulenceRates()); // DESIGN.md 3-argument contract
    const uP = turbulenceIntensity(st.k, inp.m);
    const lhs = r.dK + r.dk;
    const rhs =
      0.5 * inp.mDotIn * inp.vIn ** 2 -
      (inp.m * uP ** 3) / inp.L -
      ((st.K + st.k) * inp.mDotOut) / inp.m +
      (2 / 3) * st.k * inp.dlnRhoDt;
    expect(lhs / rhs).toBeCloseTo(1, 13);
    expect(r.production).toBeCloseTo(KK_PRODUCTION_CONSTANT * (st.K / inp.L) * Math.sqrt(st.k / inp.m), 13);
  });

  it('length scales: f_L·min(h, B) and angular-momentum conservation', () => {
    expect(integralLengthScale(0.019, 0.08255)).toBeCloseTo(0.25 * 0.019, 15);
    expect(integralLengthScale(0.13, 0.08255)).toBeCloseTo(0.25 * 0.08255, 15);
    expect(angularMomentumLengthScale(0.004, 5, 40)).toBeCloseTo(0.002, 15);
    expect(rapidDistortionIntensity(1.5, 5, 40)).toBeCloseTo(3, 14);
    expect(intakeJetVelocity(0.03, 1.1, 3.6e-4)).toBeCloseTo(0.03 / (1.1 * 3.6e-4), 12);
  });

  it('jet angular-momentum flux = y-component of ṁ (r × v)', () => {
    // jet at z = −0.02 m moving along +x: (r × v)_y = z v_x − x v_z = −0.02 v
    expect(swirlTorqueFromJet(0.01, 30, 0, -0.02, 0)).toBeCloseTo(0.01 * 30 * -0.02, 15);
    expect(swirlTorqueFromJet(0.01, 30, 0.02, 0, Math.PI / 2)).toBeCloseTo(0.01 * 30 * -0.02, 15);
  });
});

describe('K–k model robustness (review additions)', () => {
  it('DESIGN.md contract: 3-argument call uses DEFAULT_TURBULENCE_PARAMS', () => {
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m: 7e-4, mDotIn: 0.02, vIn: 40, L: 0.01, dlnRhoDt: 35 });
    const st: TurbulenceState = { K: 0.3, k: 0.05, swirl: 0 };
    const a = turbulenceDerivatives(st, inp, newTurbulenceRates());
    const b = turbulenceDerivatives(st, inp, newTurbulenceRates(), DEFAULT_TURBULENCE_PARAMS);
    expect(a).toEqual(b);
  });

  it('k = 0 with mean flow present still produces turbulence; an undershoot k < 0 relaxes back', () => {
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m: 7e-4, L: 0.005, mDotOut: 0.01 });
    const r = newTurbulenceRates();
    turbulenceDerivatives({ K: 0.05, k: 0, swirl: 0 }, inp, r);
    expect(r.dk).toBeGreaterThan(0);
    expect(r.production).toBeGreaterThan(0);
    turbulenceDerivatives({ K: 0.05, k: -1e-6, swirl: 0 }, inp, r);
    expect(r.dk).toBeGreaterThan(0);
    // the seed is negligible once turbulence exists: P/P(no seed) = 1 for k ≫ 1e-10 K
    turbulenceDerivatives({ K: 0.05, k: 1e-6, swirl: 0 }, inp, r);
    expect(r.production / (KK_PRODUCTION_CONSTANT * (0.05 / 0.005) * Math.sqrt(1e-6 / 7e-4))).toBeCloseTo(1, 14);
  });

  it('degenerate inputs give finite zero-or-sane rates (no NaN/Inf)', () => {
    const inp = newTurbulenceInputs();
    const r = newTurbulenceRates();
    const cases: Array<[Partial<TurbulenceInputs>, TurbulenceState]> = [
      [{ m: 0, L: 0.01 }, { K: 1, k: 1, swirl: 1 }],
      [{ m: 1e-3, L: 0 }, { K: 1, k: 1, swirl: 1 }],
      [{ m: 1e-3, L: 0.01, rho: 0, mu: 0 }, { K: 0, k: 0, swirl: 0.01 }],
      [{ m: 1e-3, L: 0.01, mDotIn: 0, vIn: 0 }, { K: 0, k: 0, swirl: 0 }],
      [{ m: 1e-3, L: 1e-9, dlnRhoDt: -1e4 }, { K: 1e-3, k: 1e-3, swirl: -0.01 }],
    ];
    for (const [ci, st] of cases) {
      Object.assign(inp, newTurbulenceInputs(), ci);
      turbulenceDerivatives(st, inp, r);
      for (const v of [r.dK, r.dk, r.dSwirl, r.production, r.dissipation, r.swirlFrictionTorque]) {
        expect(Number.isFinite(v)).toBe(true);
      }
    }
  });
});

describe('K–k model dynamics', () => {
  it('isotropic decay: u′(t) = u′₀ / (1 + u′₀ t / (3L)) with P = 0', () => {
    const m = 1e-3;
    const L = 0.005;
    const u0 = 5;
    const st: TurbulenceState = { K: 0, k: 1.5 * m * u0 * u0, swirl: 0 };
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m, L });
    const T = 0.02;
    integrate(st, inp, DEFAULT_TURBULENCE_PARAMS, 0, T, 4000, () => {});
    const exact = u0 / (1 + (u0 * T) / (3 * L));
    expect(turbulenceIntensity(st.k, m) / exact).toBeCloseTo(1, 9);
  });

  it('rapid isotropic compression amplifies u′ as (ρ/ρ₀)^{1/3} (Wong & Hoult 1979)', () => {
    const m = 1e-3;
    const u0 = 1;
    const st: TurbulenceState = { K: 0, k: 1.5 * m * u0 * u0, swirl: 0 };
    const inp = newTurbulenceInputs();
    // compression ratio 8 in 1 ms with L = 1 m: eddy turnover L/u′ ≈ 1 s ≫ 1 ms
    const ratio = 8;
    const T = 1e-3;
    Object.assign(inp, { m, L: 1, dlnRhoDt: Math.log(ratio) / T });
    integrate(st, inp, DEFAULT_TURBULENCE_PARAMS, 0, T, 1000, () => {});
    const uT = turbulenceIntensity(st.k, m);
    expect(uT / (u0 * Math.cbrt(ratio))).toBeCloseTo(1, 3); // dissipation removes ~1e-3 over 1 ms
    // switching the source off leaves u′ unchanged apart from dissipation
    const st2: TurbulenceState = { K: 0, k: 1.5 * m * u0 * u0, swirl: 0 };
    integrate(st2, inp, { ...DEFAULT_TURBULENCE_PARAMS, cRdt: 0 }, 0, T, 1000, () => {});
    expect(turbulenceIntensity(st2.k, m)).toBeCloseTo(u0 / (1 + (u0 * T) / 3), 9);
  });

  it('mean-flow energy cascades into turbulence (closed system, no dissipation limit)', () => {
    const m = 1e-3;
    const st: TurbulenceState = { K: 0.05, k: 1e-6, swirl: 0 };
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m, L: 0.01 });
    const E0 = st.K + st.k;
    integrate(st, inp, DEFAULT_TURBULENCE_PARAMS, 0, 0.05, 20000, () => {});
    expect(st.K).toBeLessThan(0.05);
    expect(st.k).toBeGreaterThan(1e-6);
    expect(st.K + st.k).toBeLessThan(E0); // total only decreases (dissipation)
  });

  it('swirl decays by turbulent wall friction as ω(t) = (ω₀^−0.8 + 0.8 C t)^−1.25', () => {
    const m = 7e-4;
    const B = 0.08255;
    const R = B / 2;
    const h = 0.02;
    const rho = m / (Math.PI * R * R * h);
    const mu = 3e-5;
    const w0 = 400;
    const st = initTurbulenceState(m, 0, 0, w0, B, newTurbulenceState());
    const inp = newTurbulenceInputs();
    Object.assign(inp, { m, L: 0.005, rho, mu, bore: B, h });
    const T = 0.05;
    integrate(st, inp, DEFAULT_TURBULENCE_PARAMS, 0, T, 5000, () => {});
    // dω/dt = −(π ρ 0.074 (ρ R²/μ)^−0.2 R⁴ (0.4R + h) / I) ω^1.8,  I = m B²/8
    const C = (Math.PI * rho * 0.074 * Math.pow((rho * R * R) / mu, -0.2) * R ** 4 * (0.4 * R + h)) / ((m * B * B) / 8);
    const exact = Math.pow(Math.pow(w0, -0.8) + 0.8 * C * T, -1.25);
    const w = swirlAngularVelocity(st.swirl, m, B);
    expect(w / exact).toBeCloseTo(1, 8);
    expect(w).toBeLessThan(w0);
  });
});

describe('synthetic motored CFR cycle (self-contained slider-crank + quasi-steady intake)', () => {
  // CFR F-1 geometry (ASTM D2699 engine: 3.25 in bore × 4.5 in stroke; 10 in rod — draft value in
  // engines/cfr.ts) and the draft intake-valve data there. Only a driver for the turbulence model.
  const IN = 0.0254;
  const B = 3.25 * IN;
  const S = 4.5 * IN;
  const rod = 10 * IN;
  const Ap = (Math.PI * B * B) / 4;
  const DEG = Math.PI / 180;
  const IVO = -350;
  const IVC = -146;
  const Dv = 1.2 * IN;
  const Lmax = 0.246 * IN;

  /** @param u0 residual-gas u′ at IVO as a fraction of S̄p (0 = start from rest) */
  function motored(rpm: number, cr: number, cd: number, u0 = 0.2): { ratio: number; uTDC: number; Sp: number } {
    const a = S / 2;
    const Vc = (Ap * S) / (cr - 1);
    const V = (th: number) => Vc + Ap * (rod + a - (a * Math.cos(th) + Math.sqrt(rod * rod - (a * Math.sin(th)) ** 2)));
    const dVdth = (th: number) => (V(th + 1e-6) - V(th - 1e-6)) / 2e-6;
    const omega = (rpm * 2 * Math.PI) / 60;
    const Sp = (2 * S * rpm) / 60;
    const rhoIn = 101325 / (287 * 325); // carburetted mixture at 325 K, 1 atm
    const lift = (deg: number) => (deg <= IVO || deg >= IVC ? 0 : Lmax * Math.sin((Math.PI * (deg - IVO)) / (IVC - IVO)));
    const mIVC = rhoIn * V(IVC * DEG);
    const inp = newTurbulenceInputs();
    inp.bore = B;
    // time t ↔ crank angle θ = IVO + ωt
    const setInputs = (t: number, i: TurbulenceInputs) => {
      const th = IVO * DEG + omega * t;
      const deg = th / DEG;
      const Vt = V(th);
      const Vdot = dVdth(th) * omega;
      i.h = Vt / Ap;
      i.L = integralLengthScale(i.h, B);
      i.mDotIn = 0;
      i.mDotOut = 0;
      i.vIn = 0;
      if (deg < IVC) {
        // valve open: cylinder at intake density, flow follows the piston
        i.m = rhoIn * Vt;
        i.rho = rhoIn;
        const md = rhoIn * Vdot;
        if (md > 0) {
          i.mDotIn = md;
          i.vIn = intakeJetVelocity(md, rhoIn, cd * Math.PI * Dv * Math.max(lift(deg), 1e-5));
        } else i.mDotOut = -md;
        i.dlnRhoDt = 0;
      } else {
        i.m = mIVC;
        i.rho = mIVC / Vt;
        i.dlnRhoDt = -Vdot / Vt;
      }
    };
    const m0 = rhoIn * V(IVO * DEG);
    const st: TurbulenceState = { K: 0, k: 1.5 * m0 * (u0 * Sp) ** 2, swirl: 0 }; // residual turbulence
    const tEnd = (-IVO * DEG) / omega;
    integrate(st, inp, DEFAULT_TURBULENCE_PARAMS, 0, tEnd, 7000, setInputs);
    const uTDC = turbulenceIntensity(st.k, mIVC);
    return { ratio: uTDC / Sp, uTDC, Sp };
  }

  it("u′(TDC) = 0.4–0.6 × mean piston speed (Bopp, Vafidis & Whitelaw 1986: 0.45–0.6 without swirl)", () => {
    const rows: string[] = [];
    for (const cr of [5, 7, 10]) {
      for (const rpm of [600, 900]) {
        const r = motored(rpm, cr, 0.6);
        rows.push(`CR ${cr} ${rpm} rpm: u′(TDC) = ${r.uTDC.toFixed(3)} m/s = ${r.ratio.toFixed(3)} S̄p`);
        expect(r.ratio).toBeGreaterThan(0.4);
        expect(r.ratio).toBeLessThan(0.6);
      }
    }
    console.log('[turbulence] ' + rows.join('; '));
  });

  it('u′(TDC) scales with engine speed and forgets the intake details (valve C_D)', () => {
    const a = motored(600, 7, 0.6);
    const b = motored(900, 7, 0.6);
    expect(b.ratio / a.ratio).toBeCloseTo(1, 3);
    const c = motored(600, 7, 0.4);
    expect(Math.abs(c.ratio / a.ratio - 1)).toBeLessThan(0.05);
    // CFR RON conditions (600 rpm, CR ≈ 7): u′ ≈ 1.2 m/s
    expect(a.uTDC).toBeGreaterThan(1.0);
    expect(a.uTDC).toBeLessThan(1.35);
  });

  it('turbulence develops from rest (k = K = 0 at IVO) and forgets the residual u′ (review fix: P ∝ √k seed)', () => {
    // Before the fix, P = 0.3307 (K/L)√(k/m) kept k ≡ 0 for ever: the intake jet filled K but u′ stayed 0.
    const rest = motored(600, 7, 0.6, 0);
    const ref = motored(600, 7, 0.6, 0.2);
    const hot = motored(600, 7, 0.6, 1.0);
    expect(rest.ratio).toBeGreaterThan(0.4);
    expect(Math.abs(rest.ratio / ref.ratio - 1)).toBeLessThan(0.02);
    expect(Math.abs(hot.ratio / ref.ratio - 1)).toBeLessThan(0.02);
  });
});
