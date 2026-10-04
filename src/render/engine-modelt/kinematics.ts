/**
 * Snapshot → part poses of the Model T mechanism. Pure math, allocation-free (writes into a pose
 * created once by createModelTPose).
 *
 * Conventions (engine-spec.ts, render/engine/kinematics.ts):
 *  - θ = snapshot.thetaDeg = cylinder 1's angle (the engine angle); cylinder i runs at
 *    θ_i = cylinderAngleDeg(spec, i, θ). The crank turns clockwise seen from the front: crankshaft
 *    rotation about +z = −θ; cylinder i's pin sits at r·(sin θ_i, cos θ_i) (its throw is at the
 *    crank-local angle −firingOffset_i, so throws 1&4 / 2&3 come out 0° / 180° for 1-2-4-3).
 *  - Piston i: wrist pin at (0, y_pin,TDC − disp_i, axisZ[i]); disp_i from the snapshot
 *    (s.cylinders[i] / top level for cylinder 1) when given, else the exact slider-crank at θ_i.
 *    The rod is drawn crank pin → wrist pin.
 *  - Camshaft: gear-driven at half speed in the opposite sense, rotation +θ/2 (as the CFR); every lobe
 *    is the same outline rotated by its phase −(offset_i + θ_c)/2 (layout.ts), so its nose meets the
 *    tappet at θ_i = θ_c.
 *  - Side valves open UPWARD (+y) by the snapshot lift (else the cam model at θ_i). The tappet rides the
 *    lobe: face height = cam axis + support function of the lobe at the current cam angle, i.e. the
 *    gross cam lift at θ_i (the lash gap to the stem end closes while the valve is open).
 *  - Fan: belt-driven in the crank's sense at the pulley ratio, rotation −θ·ratio.
 */
import { cylinderAngleDeg } from '../../physics/core/engine-spec';
import { pistonDisplacementAt, type ModelTLayout } from './layout';

const DEG = Math.PI / 180;

export interface ModelTKinematicInput {
  thetaDeg: number;
  /** Per-cylinder piston displacement below TDC, m (NaN: derive from θ_i). */
  pistonDisplacement: Float64Array;
  /** Per-cylinder intake / exhaust valve lift, m (NaN: derive from the cam at θ_i). */
  intakeLift: Float64Array;
  exhaustLift: Float64Array;
}

export interface CylinderPoseMT {
  /** Local crank angle θ_i, deg. */
  thetaDeg: number;
  pistonDisplacement: number;
  /** Crank-pin centre (x, y), ROOT. */
  crankPin: [number, number];
  /** Wrist-pin centre height, ROOT. */
  wristPinY: number;
  /** Rod rotation about +z (rod-local +y from crank pin to wrist pin), rad. */
  rodRotation: number;
  /** |W − P| (= rod length for a consistent snapshot). */
  rodSpan: number;
  /** Crown height in the cylinder frame (= −h_i). */
  crownYCyl: number;
}

export interface ValvePoseMT {
  /** Valve lift, m (≥ 0). */
  lift: number;
  /** Gross tappet (cam) lift, m. */
  tappetLift: number;
  /** Tappet flat-face height, ROOT y. */
  tappetFaceY: number;
  /** Spring length / installed length. */
  springScale: number;
}

export interface ModelTPose {
  thetaDeg: number;
  crankAngle: number;
  camAngle: number;
  fanAngle: number;
  cylinders: CylinderPoseMT[];
  valves: ValvePoseMT[];
}

export function createModelTInput(L: ModelTLayout): ModelTKinematicInput {
  const f = () => new Float64Array(L.nCyl).fill(NaN);
  return { thetaDeg: 0, pistonDisplacement: f(), intakeLift: f(), exhaustLift: f() };
}

export function createModelTPose(L: ModelTLayout): ModelTPose {
  return {
    thetaDeg: 0,
    crankAngle: 0,
    camAngle: 0,
    fanAngle: 0,
    cylinders: L.axisZ.map(() => ({ thetaDeg: 0, pistonDisplacement: 0, crankPin: [0, 0] as [number, number], wristPinY: 0, rodRotation: 0, rodSpan: 0, crownYCyl: 0 })),
    valves: L.valves.map(() => ({ lift: 0, tappetLift: 0, tappetFaceY: 0, springScale: 1 })),
  };
}

/** Compute every moving-part pose. Allocation-free. */
export function computeModelTPose(L: ModelTLayout, input: ModelTKinematicInput, out: ModelTPose): ModelTPose {
  const th = Number.isFinite(input.thetaDeg) ? input.thetaDeg : 0;
  const r = L.crankRadius;
  out.thetaDeg = th;
  out.crankAngle = -th * DEG;
  out.camAngle = (th * DEG) / 2;
  out.fanAngle = -th * DEG * L.fan.ratio;
  for (let i = 0; i < L.nCyl; i++) {
    const c = out.cylinders[i];
    const ti = cylinderAngleDeg(L.spec, i, th);
    c.thetaDeg = ti;
    const phi = ti * DEG;
    const px = r * Math.sin(phi), py = r * Math.cos(phi);
    c.crankPin[0] = px;
    c.crankPin[1] = py;
    const d = input.pistonDisplacement[i];
    const disp = Number.isFinite(d) ? d : pistonDisplacementAt(L, ti);
    c.pistonDisplacement = disp;
    const wy = L.wristPinTdcY - disp;
    c.wristPinY = wy;
    const dx = -px, dy = wy - py;
    c.rodSpan = Math.hypot(dx, dy);
    c.rodRotation = Math.atan2(-dx, dy);
    c.crownYCyl = wy + L.compressionHeight - L.cylOriginY;
  }
  const sp = L.spring;
  for (let k = 0; k < L.valves.length; k++) {
    const v = L.valves[k];
    const vp = out.valves[k];
    const ti = out.cylinders[v.cyl].thetaDeg;
    const raw = v.kind === 'intake' ? input.intakeLift[v.cyl] : input.exhaustLift[v.cyl];
    const lift = Number.isFinite(raw) ? Math.max(0, raw) : v.cam.valveLift(ti);
    vp.lift = lift;
    const tl = v.cam.tappetLift(ti);
    vp.tappetLift = tl;
    vp.tappetFaceY = L.cam.y + v.cam.baseRadius + tl;
    vp.springScale = Math.max(0.2, (sp.installed - lift) / sp.installed);
  }
  return out;
}
