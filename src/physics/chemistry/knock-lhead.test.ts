import { describe, expect, it } from 'vitest';
import fx from '../../../test/fixtures/chemistry_knock_lhead_modes.json';
import { lensArea } from '../combustion/flame-geometry';
import type { EngineSpec } from '../core/engine-spec';
import { CFR_F1, CFR_KNOCK_PICKUP } from '../engines/cfr';
import { MODEL_T } from '../engines/model-t';
import {
  besselJPrimeZeros,
  type ChamberPlanform,
  createKnockOscillator,
  DiscEndGasProjector,
  DiscModeSet,
  discPlanform,
  GridEndGasProjector,
  GridModeSet,
  KnockOscillator,
  L_HEAD_MAPO_BAND,
  lHeadDepthAtTDC,
  lHeadKnockAcoustics,
  lHeadPlanform,
  symmetricEigen,
  virtualKnockSensor,
} from './knock';

const now = (): number => performance.now();
const B = MODEL_T.geometry.bore;
const lHead = MODEL_T.geometry.lHead!;
const gap = MODEL_T.sparkPlug.gapCenter;

describe('symmetric eigen-solver (Householder + implicit QL)', () => {
  it('diagonalises random symmetric matrices: A v = λ v, VᵀV = I, ascending', () => {
    let s = 12345;
    const rnd = (): number => {
      s = (Math.imul(s, 1103515245) + 12345) >>> 0;
      return s / 4294967296 - 0.5;
    };
    for (const n of [1, 2, 7, 40]) {
      const A = new Float64Array(n * n);
      for (let i = 0; i < n; i++) for (let j = 0; j <= i; j++) A[i * n + j] = A[j * n + i] = rnd();
      const V = Float64Array.from(A);
      const d = new Float64Array(n);
      symmetricEigen(V, n, d);
      for (let e = 0; e < n; e++) {
        if (e > 0) expect(d[e]).toBeGreaterThanOrEqual(d[e - 1]);
        for (let i = 0; i < n; i++) {
          let av = 0;
          for (let j = 0; j < n; j++) av += A[i * n + j] * V[j * n + e];
          expect(Math.abs(av - d[e] * V[i * n + e])).toBeLessThan(1e-12);
        }
        for (let f = 0; f < n; f++) {
          let dot = 0;
          for (let i = 0; i < n; i++) dot += V[i * n + e] * V[i * n + f];
          expect(Math.abs(dot - (e === f ? 1 : 0))).toBeLessThan(1e-12);
        }
      }
    }
  });
});

describe('GridModeSet: depth-weighted planform modes (analytic limits)', () => {
  it('flat disc: reproduces the Draper/Bessel α_mn (all 16 modes ≤ 7.1, degenerate pairs included) within 0.2 %', () => {
    const Bc = CFR_F1.geometry.bore;
    const g = new GridModeSet(Bc, discPlanform(Bc / 2, 0.01));
    const ref: number[] = [];
    for (let m = 0; m < 9; m++) for (const z of besselJPrimeZeros(m, 7.1)) for (let c = 0; c < (m === 0 ? 1 : 2); c++) ref.push(z);
    ref.sort((a, b) => a - b);
    expect(g.nModes).toBe(ref.length);
    expect(g.nModes).toBe(new DiscModeSet(Bc).nModes);
    for (let j = 0; j < g.nModes; j++) expect(Math.abs(g.alpha[j] / ref[j] - 1)).toBeLessThan(2e-3);
    // the same discretisation solved by scipy (ARPACK shift-invert) — the Krylov solver is exact
    fx.disc.alpha.forEach((a, j) => expect(Math.abs(g.alpha[j] / a - 1)).toBeLessThan(1e-9));
    fx.disc.besselReference.forEach((a, j) => expect(ref[j]).toBeCloseTo(a, 9));
    expect(Math.abs(g.stats.nullEigenvalue)).toBeLessThan(1e-9 * (2 / Bc) ** 2);
  });

  it('depth step: two-depth channel H1 tan(kL1) + H2 tan(kL2) = 0 within 0.2 %', () => {
    const c = fx.channel;
    const L = c.L1 + c.L2;
    const pf: ChamberPlanform = {
      xMin: 0,
      xMax: L,
      zMin: 0,
      zMax: c.width,
      yLo: [-c.H1, -c.H2],
      yHi: [0, 0],
      region: (x, z) => (x < 0 || x > L || z < 0 || z > c.width ? -1 : x < c.L1 ? 0 : 1),
      boundaryDistance: (x, z) => Math.min(Math.abs(x), Math.abs(x - c.L1), Math.abs(x - L), Math.abs(z), Math.abs(z - c.width)),
    };
    const g = new GridModeSet(0.1, pf, { cellSize: c.cellSize, maxAlpha: 8 });
    const k = Array.from(g.alpha, (a) => (2 * a) / 0.1);
    c.kExact.forEach((ke, j) => expect(Math.abs(k[j] / ke - 1)).toBeLessThan(2e-3));
    c.kSameGrid.slice(0, 5).forEach((ks, j) => expect(Math.abs(k[j] / ks - 1)).toBeLessThan(1e-9));
    // a uniform channel would have k = jπ/L: the step moves the fundamental by +11 %
    expect(k[0] / (Math.PI / L)).toBeGreaterThan(1.05);
  });

  it('modes are C-orthogonal and orthogonal to the uniform mode; modeShape interpolates the cell values', () => {
    const g = lHeadKnockAcoustics(MODEL_T).modes;
    const N = g.nModes;
    for (let i = 0; i < N; i++) {
      let s0 = 0;
      for (let c = 0; c < g.nCells; c++) s0 += g.capacity[c] * g.psi[c * N + i];
      expect(Math.abs(s0) / g.volume).toBeLessThan(1e-10);
      for (let j = i + 1; j < N; j++) {
        let s = 0;
        for (let c = 0; c < g.nCells; c++) s += g.capacity[c] * g.psi[c * N + i] * g.psi[c * N + j];
        expect(Math.abs(s) / g.volume / Math.sqrt(g.meanSquare[i] * g.meanSquare[j])).toBeLessThan(1e-8);
      }
    }
    // at a cell centre the interpolant is the cell value; outside the chamber 0
    const i = 40;
    const k = 30;
    const c = g.cellIndex[i + g.nx * k];
    expect(c).toBeGreaterThanOrEqual(0);
    const x = g.gx0 + (i + 0.5) * g.cellSize;
    const z = g.gz0 + (k + 0.5) * g.cellSize;
    for (let j = 0; j < N; j++) expect(g.modeShape(j, x, z)).toBeCloseTo(g.psi[c * N + j], 12);
    expect(g.modeShape(0, 0.06, 0.06)).toBe(0);
  });
});

describe('Model T L-head modes (oracle: scipy on the same grid and on finer grids)', () => {
  const ac = lHeadKnockAcoustics(MODEL_T);
  const g = ac.modes;

  it('fixture geometry is the current MODEL_T chamber', () => {
    expect(fx.geometry.bore).toBe(B);
    expect(fx.geometry.deckY).toBe(lHead.deckY);
    expect(fx.geometry.crownAboveDeckAtTDC).toBe(lHead.crownAboveDeckAtTDC);
    expect(fx.geometry.pocket).toEqual(lHead.pocket);
    expect(fx.geometry.gapCenter).toEqual(gap);
    expect(lHeadDepthAtTDC(lHead)).toBeCloseTo(fx.geometry.depthTDC, 15);
  });

  it('same discretisation as the scipy oracle: cells, volume, every α_j to 1e-9', () => {
    const s = fx.sameGrid;
    expect(g.cellSize).toBeCloseTo(s.cellSize, 15);
    expect(g.nCells).toBe(s.nCells);
    expect(g.columns.count).toBe(s.nColumns);
    expect(g.volume / s.volume - 1).toBeLessThan(1e-12);
    expect(g.nModes).toBe(s.alpha.length);
    s.alpha.forEach((a, j) => expect(Math.abs(g.alpha[j] / a - 1)).toBeLessThan(1e-9));
    s.meanSquare.forEach((m, j) => expect(Math.abs(g.meanSquare[j] / m - 1)).toBeLessThan(1e-7));
    expect(g.stats.maxResidual).toBeLessThan(1e-10);
  });

  it('chamber volume = clearance volume − crevice (exact geometry within 0.02 %)', () => {
    expect(Math.abs(g.volume / fx.geometry.chamberVolumeExact - 1)).toBeLessThan(2e-4);
    expect(Math.abs(g.volume / fx.endGasVolume.volume - 1)).toBeLessThan(2e-4);
  });

  it('grid-converged: α_j within 0.3 % of the B/192 solution; fundamental ≈ 3.8 kHz at 950 m/s (2–4 kHz expected)', () => {
    const fine = fx.converged[1];
    expect(fine.cellsPerBore).toBe(192);
    for (let j = 0; j < 12; j++) expect(Math.abs(g.alpha[j] / fine.alpha[j] - 1)).toBeLessThan(3e-3);
    const f1 = (g.alpha[0] * 950) / (Math.PI * B);
    console.log(
      `[knock] Model T L-head modes (B/64, ${g.nCells} cells, ${g.nModes} modes, build ${g.stats.buildMs.toFixed(0)} ms, Krylov ${g.stats.krylovDimension}): ` +
        `α = ${Array.from(g.alpha.slice(0, 6), (a) => a.toFixed(4)).join(', ')} …; f1 = ${f1.toFixed(0)} Hz at c = 950 m/s ` +
        `(bore (1,0): ${((1.8412 * 950) / (Math.PI * B)).toFixed(0)} Hz)`,
    );
    expect(f1).toBeGreaterThan(2000);
    expect(f1).toBeLessThan(4000);
    expect(f1).toBeCloseTo(3812, -1);
    // depth weighting matters: the uniform-depth planform's fundamental is 7.5 % lower
    const flat = new GridModeSet(B, { ...lHeadPlanform(B, lHead), yLo: [-0.01, -0.01], yHi: [0, 0] }, { maxAlpha: 2 });
    expect(flat.alpha[0] / g.alpha[0]).toBeLessThan(0.95);
  });

  it('bore-column depth at 15° ATDC moves the fundamental by +1 % (scipy B/96)', () => {
    const d15 = fx.depth15;
    const g15 = new GridModeSet(B, lHeadPlanform(B, lHead, d15.depth), { maxAlpha: 3.5 });
    for (let j = 0; j < 6; j++) expect(Math.abs(g15.alpha[j] / d15.alpha[j] - 1)).toBeLessThan(2e-3);
    expect(g15.alpha[0] / g.alpha[0] - 1).toBeGreaterThan(0.005);
    expect(g15.alpha[0] / g.alpha[0] - 1).toBeLessThan(0.02);
  });
});

describe('GridEndGasProjector (depth-weighted end gas outside the flame ball about the plug)', () => {
  const ac = lHeadKnockAcoustics(MODEL_T);
  const P = ac.projector;
  const g = ac.modes;

  it('end-gas volume outside r: same-grid scipy to 1e-12, exact chamber geometry within 0.2 % of V', () => {
    const s = fx.sameGrid;
    const e = fx.endGasVolume;
    s.radii.forEach((r, i) => {
      expect(Math.abs(P.outsideMeasure(r) - s.vOut[i]) / s.volume).toBeLessThan(1e-12);
      expect(Math.abs(P.outsideMeasure(r) / P.totalMeasure - e.vOut[i] / e.volume)).toBeLessThan(2e-3);
    });
    expect(P.outsideMeasure(0)).toBe(P.totalMeasure);
    expect(P.outsideMeasure(P.maxRadius)).toBe(0);
    // continuous and non-increasing in r
    let prev = Infinity;
    for (let r = 0; r <= P.maxRadius; r += P.maxRadius / 997) {
      const v = P.outsideMeasure(r);
      expect(v).toBeLessThanOrEqual(prev * (1 + 1e-14));
      prev = v;
    }
  });

  it('radiusForOutsideFraction inverts the end-gas volume (same 20 % radius as scipy)', () => {
    for (const f of [0.999, 0.5, 0.2, 0.01, 1e-4]) {
      const r = P.radiusForOutsideFraction(f);
      expect(Math.abs(P.outsideMeasure(r) / P.totalMeasure - f)).toBeLessThan(1e-12);
    }
    expect(Math.abs(P.radiusForOutsideFraction(0.2) - fx.sameGrid.r20)).toBeLessThan(1e-9);
    expect(P.radiusForOutsideFraction(1)).toBe(0);
    expect(P.radiusForOutsideFraction(0)).toBe(P.maxRadius);
  });

  it('sensor-weighted source shapes ψ_j(plug)·S_j: same-grid scipy to 1e-6, B/192 within 4 %', () => {
    const ko = createKnockOscillator(MODEL_T);
    ko.setEndGasOutsideRadius(P.radiusForOutsideFraction(0.2));
    const w = Array.from(ko.sourceShape, (s, j) => s * ko.psiSensor[j]);
    fx.sameGrid.weights.forEach((ref, j) => expect(Math.abs(w[j] - ref)).toBeLessThan(1e-6 * Math.max(1, Math.abs(ref))));
    const fine = fx.converged[1].weights;
    for (let j = 0; j < 12; j++) expect(Math.abs(w[j] - fine[j])).toBeLessThan(0.04 * Math.max(Math.abs(fine[j]), 0.1));
    // the plug sits on the planform's mirror plane z = 0: modes odd in z neither ring there nor are excited
    expect(Math.abs(ko.psiSensor[1])).toBeLessThan(1e-9);
    expect(Math.abs(ko.sourceShape[1])).toBeLessThan(1e-9);
  });

  it('a uniform release excites no mode; shells partition the end gas (volume-weighted mean = whole region)', () => {
    const ko = createKnockOscillator(MODEL_T);
    expect(ko.setEndGasOutsideRadius(0)).toBe(1);
    for (let j = 0; j < ko.nModes; j++) expect(Math.abs(ko.sourceShape[j])).toBeLessThan(1e-12);
    const rIn = P.radiusForOutsideFraction(0.3);
    ko.setEndGasOutsideRadius(rIn);
    const whole = Float64Array.from(ko.sourceShape);
    const K = 16;
    expect(ko.setSequentialEndGasOutsideRadius(rIn, K)).toBeCloseTo(0.3, 12);
    const mean = new Float64Array(ko.nModes);
    for (let k = 0; k < K; k++) {
      ko.setSequentialShell(k);
      for (let j = 0; j < ko.nModes; j++) mean[j] += ko.sourceShape[j] / K; // equal-volume shells
    }
    for (let j = 0; j < ko.nModes; j++) expect(Math.abs(mean[j] - whole[j])).toBeLessThan(1e-9 * Math.max(1, Math.abs(whole[j])));
    // shell 0 (farthest from the plug) differs from the last (next to the flame)
    ko.setSequentialShell(0);
    const s0 = ko.sourceShape[0];
    ko.setSequentialShell(K - 1);
    expect(Math.abs(s0 - ko.sourceShape[0])).toBeGreaterThan(0.1);
  });

  it('discrete Parseval: with the complete mode set Σ⟨ψ_j²⟩S_j² = Σc_i e_i²/(v²V) − 1; partial sums grow toward it', () => {
    // coarse grid so that every discrete mode is computed (Krylov space exhausted, restarted)
    const coarse = new GridModeSet(B, lHeadPlanform(B, lHead), { cellSize: B / 10, maxAlpha: 1e6 });
    expect(coarse.nModes).toBe(coarse.nCells - 1);
    const pc = new GridEndGasProjector(coarse, gap);
    const ko = new KnockOscillator(B, { modes: coarse, projector: pc, decayTime: 1e9 });
    const rEg = pc.radiusForOutsideFraction(0.25);
    const v = ko.setEndGasOutsideRadius(rEg);
    let sum = 0;
    for (let j = 0; j < ko.nModes; j++) sum += ko.meanSquare[j] * ko.sourceShape[j] ** 2;
    // end-gas volume of every cell from the sub-column chords (independent of the projector's sorting/suffix sums)
    const col = coarse.columns;
    const veg = new Float64Array(coarse.nCells);
    for (let q = 0; q < col.count; q++) {
      const w = Math.sqrt(Math.max(0, rEg * rEg - (col.x[q] - gap[0]) ** 2 - (col.z[q] - gap[2]) ** 2));
      const chord = Math.max(0, Math.min(col.yHi[q] - gap[1], w) - Math.max(col.yLo[q] - gap[1], -w));
      veg[col.cell[q]] += col.area[q] * (col.yHi[q] - col.yLo[q] - chord);
    }
    let ce2 = 0;
    for (let c = 0; c < coarse.nCells; c++) ce2 += veg[c] ** 2 / coarse.capacity[c];
    const exact = ce2 / (v * v * coarse.volume) - 1;
    expect(Math.abs(sum / exact - 1)).toBeLessThan(1e-9);
    // = (1 − v)/v when no cell is partly end gas; fractional coverage lowers it (−4.5 % at B/10)
    expect(exact).toBeLessThan((1 - v) / v);
    expect(exact).toBeGreaterThan(0.9 * ((1 - v) / v));
    // instantaneous release: modal energy → (γ−1)²(Q/V)²V(1 − v)/v/(2ρc²) (the default grid, more modes ⇒ closer)
    const ratios = [2.5, 5, 7.1].map((maxAlpha) => {
      const m = new GridModeSet(B, lHeadPlanform(B, lHead), { maxAlpha });
      const o = new KnockOscillator(B, { modes: m, projector: new GridEndGasProjector(m, gap), decayTime: 1e9 });
      const vv = o.setEndGasOutsideRadius(o.projector.radiusForOutsideFraction(0.15));
      const gm = 1.3;
      const V = 2.4e-4;
      const Q = 50;
      const dt = 1e-13;
      o.step(dt, 900, gm, V, Q / dt);
      o.step(dt, 900, gm, V, 0);
      const eExact = (V * (gm - 1) ** 2 * (Q / V) ** 2 * (1 - vv)) / vv / (2 * 20 * 900 * 900);
      return o.acousticEnergy(20, 900, V) / eExact;
    });
    expect(ratios[0]).toBeLessThan(ratios[1]);
    expect(ratios[1]).toBeLessThan(ratios[2]);
    expect(ratios[2]).toBeLessThan(1);
    expect(ratios[2]).toBeGreaterThan(0.5);
  });

  it('small end gas: finite, continuous source (no thin-crescent zero) and MAPO > 0', () => {
    const ko = createKnockOscillator(MODEL_T);
    const s: number[] = [];
    for (const f of [1e-2, 3e-3, 1e-3, 3e-4]) {
      ko.reset();
      expect(ko.setEndGasVolumeFraction(f, 16)).toBeCloseTo(f, 10);
      s.push(ko.sourceShape[0]);
      for (let i = 0; i < 400; i++) {
        ko.setSequentialPosition(Math.min(1, i / 20));
        ko.step(1e-6, 900, 1.3, 2.4e-4, i < 20 ? (f * 1500) / 20e-6 : 0);
      }
      expect(ko.peakSensorPressure).toBeGreaterThan(0);
    }
    for (let i = 1; i < s.length; i++) expect(Math.abs(s[i] / s[i - 1] - 1)).toBeLessThan(0.05);
  });
});

describe('Model T knock event and the engine factory', () => {
  it('createKnockOscillator: shared acoustics, virtual sensor at the plug, L-head band', () => {
    const ac = lHeadKnockAcoustics(MODEL_T);
    const cyl = [0, 1, 2, 3].map(() => createKnockOscillator(MODEL_T));
    for (const ko of cyl) {
      expect(ko.modes).toBe(ac.modes);
      expect(ko.projector).toBe(ac.projector);
    }
    expect(lHeadKnockAcoustics({ ...MODEL_T } as EngineSpec)).toBe(ac); // memoised by geometry
    expect(virtualKnockSensor(MODEL_T)).toEqual([gap[0], gap[2]]);
    const ko = cyl[0];
    expect(ko.psiSensor[0]).toBeCloseTo(ac.modes.modeShape(0, gap[0], gap[2]), 14);
    // band: every mode ≤ 18 kHz at the reference sound speed counts, including the fundamental
    ko.setBandReference(950);
    expect(ko.bandGain(0, 950)).toBe(1);
    expect(L_HEAD_MAPO_BAND[0]).toBeLessThan(ko.modeFrequency(0, 600));
    const unfiltered = createKnockOscillator(MODEL_T, { band: null });
    unfiltered.setBandReference(950);
    expect(unfiltered.bandGain(ko.nModes - 1, 950)).toBe(1);
    // planform-circle API is disc-only
    expect(() => ko.setEndGasRegion(0, 0, 0.01)).toThrow();
    expect(() => new KnockOscillator(B, { modes: ac.modes })).toThrow();
  });

  it('flat-disc chamber: identical to new KnockOscillator(bore, opts); setEndGasVolumeFraction = the cycle-model source', () => {
    const sensor: [number, number] = [CFR_KNOCK_PICKUP.position[0], CFR_KNOCK_PICKUP.position[2]];
    const a = createKnockOscillator(CFR_F1, { sensor, band: [4000, 18000] });
    const b = new KnockOscillator(CFR_F1.geometry.bore, { sensor, band: [4000, 18000] });
    expect(a.modes).toBeInstanceOf(DiscModeSet);
    expect(a.projector).toBeInstanceOf(DiscEndGasProjector);
    expect(Array.from(a.psiSensor)).toEqual(Array.from(b.psiSensor));
    expect(Array.from(a.meanSquare)).toEqual(Array.from(b.meanSquare));
    // cycle-model.ts onAutoignition: rIn = endGasCircleRadius(f) about the plug, then the sequential source
    const g = CFR_F1.sparkPlug.gapCenter;
    for (const f of [0.0025, 0.1, 0.6]) {
      const R = 0.5 * CFR_F1.geometry.bore;
      const d = Math.max(Math.hypot(g[0], g[2]), 1e-9);
      const target = (1 - f) * Math.PI * R * R;
      let lo = 0;
      let hi = R + d;
      for (let i = 0; i < 60; i++) {
        const mid = 0.5 * (lo + hi);
        if (lensArea(mid, R, d) < target) lo = mid;
        else hi = mid;
      }
      const rIn = 0.5 * (lo + hi);
      for (const K of [16, 0]) {
        const fa = a.setEndGasVolumeFraction(f, K);
        const fb = K > 0 ? b.setSequentialEndGasRegion(g[0], g[2], rIn, K) : b.setEndGasRegion(g[0], g[2], rIn);
        expect(fa).toBe(fb);
        expect(Array.from(a.sourceShape)).toEqual(Array.from(b.sourceShape));
        // an oscillator built without `spark` (code-review-r2.test.ts swaps one into the cycle) given the centre
        const c = new KnockOscillator(CFR_F1.geometry.bore, { sensor, band: [4000, 18000] });
        expect(c.setEndGasVolumeFraction(f, K, g)).toBe(fb);
        expect(Array.from(c.sourceShape)).toEqual(Array.from(b.sourceShape));
      }
    }
    const lt = createKnockOscillator(MODEL_T);
    expect(() => lt.setEndGasVolumeFraction(0.1, 16, gap)).not.toThrow();
    expect(() => lt.setEndGasVolumeFraction(0.1, 16, [0, 0, 0])).toThrow();
  });

  it('Model T: 10 % end gas, 150 J at the plug sensor — rings at the L-head fundamental; MAPO falls with τ_ab', () => {
    const ko = createKnockOscillator(MODEL_T, { band: null, decayTime: 1e-3 });
    const t0 = now();
    const n = 20;
    for (let i = 0; i < n; i++) ko.setEndGasVolumeFraction(0.1 + 1e-4 * i, 16);
    const perEvent = (now() - t0) / n;
    const t1 = now();
    for (let i = 0; i < n; i++) ko.setEndGasVolumeFraction(0.6 + 1e-4 * i, 16);
    const perEventLarge = (now() - t1) / n;
    console.log(`[knock] Model T knock source (16 shells): ${perEvent.toFixed(2)} ms per event at 10 % end gas, ${perEventLarge.toFixed(2)} ms at 60 % (CFR disc ≈ 1.5 ms)`);
    expect(perEvent).toBeLessThan(20);
    const V = 2.6e-4;
    const c = 900;
    const run = (tauAb: number): { mapo: number; f: number } => {
      ko.reset();
      ko.setEndGasVolumeFraction(0.1, 16);
      let m = 0;
      const dt = 1e-6;
      const sig: number[] = [];
      for (let i = 0; i < 3000; i++) {
        const t = (i + 0.5) * dt;
        ko.setSequentialPosition(1 - Math.exp(-t / tauAb));
        ko.step(dt, c, 1.3, V, (150 / tauAb) * Math.exp(-t / tauAb));
        const p = ko.sensorPressure();
        m = Math.max(m, Math.abs(p));
        sig.push(p);
      }
      // dominant spectral line of the sensor signal (DFT magnitude scan, 1–20 kHz in 5 Hz steps)
      let fPeak = 0;
      let aPeak = 0;
      for (let f = 1000; f <= 20000; f += 5) {
        let re = 0;
        let im = 0;
        for (let i = 0; i < sig.length; i++) {
          re += sig[i] * Math.cos(2 * Math.PI * f * i * dt);
          im += sig[i] * Math.sin(2 * Math.PI * f * i * dt);
        }
        if (re * re + im * im > aPeak) {
          aPeak = re * re + im * im;
          fPeak = f;
        }
      }
      return { mapo: m, f: fPeak };
    };
    const r = [1e-6, 20e-6, 50e-6, 100e-6].map(run);
    for (let i = 1; i < r.length; i++) expect(r[i].mapo).toBeLessThan(r[i - 1].mapo);
    const f1 = ko.modeFrequency(0, c);
    console.log(`[knock] Model T 150 J, 10 % end gas: MAPO ${r.map((x) => (x.mapo / 1e5).toFixed(1)).join(' / ')} bar (τ_ab 1/20/50/100 µs); ringing ${r[3].f.toFixed(0)} Hz vs f1 ${f1.toFixed(0)} Hz`);
    // the plug sensor sees mainly the fundamental (bore ↔ pocket sloshing)
    expect(Math.abs(r[3].f / f1 - 1)).toBeLessThan(0.03);
    expect(r[0].mapo / 1e5).toBeGreaterThan(5);
    expect(r[0].mapo / 1e5).toBeLessThan(100);
  });
});
