/**
 * Snapshot → part poses for the mechanism. Pure math (no three.js) so the
 * kinematic consistency can be unit-tested.
 *
 * Conventions (ROOT frame, see layout.ts; identical to physics/mechanics/kinematics.ts):
 *  - Crank angle θ (snapshot.thetaDeg, 0 = firing TDC). The crank turns
 *    CLOCKWISE seen from +z (the default viewing side): crank-pin axis angle
 *    φ = θ + φ_TDC (φ_TDC ≠ 0 only with a pin offset), P = r (sin φ, cos φ);
 *    crankshaft rotation about +z = −φ.
 *  - Wrist pin W = (e, y_TDC − pistonDisplacement), e = spec pin offset
 *    (the piston is placed from the snapshot's pistonDisplacement, so the
 *    crown is exactly where the physics puts it).
 *  - The rod is drawn from P toward W. For a consistent snapshot |W − P| = L and
 *    the rod rotation equals the physics rod angle β (sin β = (r sin φ − e)/L).
 *  - Camshaft: gear-driven (opposite sense) at half speed: rotation +θ/2.
 *  - Valves: lift from the snapshot. Rocker angle α = asin(lift / a_v); the
 *    pushrod end rises a_p sin α = lift / ROCKER_RATIO, which is the tappet lift.
 */
import type { EngineLayout, Vec3 } from './layout';

export interface KinematicInput {
  thetaDeg: number;
  pistonDisplacement: number;
  intakeLift: number;
  exhaustLift: number;
}

export interface ValvePose {
  /** Valve lift, m (≥ 0). */
  lift: number;
  /** Rocker rotation about its pivot axis, rad (positive = valve side down). */
  rockerAngle: number;
  /** Tappet (and pushrod) lift, m. */
  tappetLift: number;
  /** Tappet flat face height, ROOT y. */
  tappetFaceY: number;
  /** Pushrod ball ends, ROOT frame. */
  pushrodBottom: Vec3;
  pushrodTop: Vec3;
}

export interface EnginePose {
  thetaDeg: number;
  /** Crankshaft (and flywheel, crank gear) rotation about +z, rad. */
  crankAngle: number;
  /** Crank-pin centre (x, y), ROOT. */
  crankPin: [number, number];
  /** Wrist-pin centre (x, y), ROOT. */
  wristPin: [number, number];
  /** Rod rotation about +z so that rod-local +y points from crank pin to wrist pin. */
  rodRotation: number;
  /** |W − P| (equals the rod length for a consistent snapshot). */
  rodSpan: number;
  /** Piston crown height, ROOT y. */
  crownY: number;
  /** Cylinder-frame origin height (fire deck), ROOT y. */
  headY: number;
  /** Crown position in the cylinder frame (−h). */
  crownYCyl: number;
  /** Camshaft rotation about +z, rad. */
  camAngle: number;
  /** Motor pulley rotation about +z, rad. */
  motorAngle: number;
  valves: [ValvePose, ValvePose];
}

export function createPose(): EnginePose {
  const vp = (): ValvePose => ({
    lift: 0, rockerAngle: 0, tappetLift: 0, tappetFaceY: 0, pushrodBottom: [0, 0, 0], pushrodTop: [0, 0, 0],
  });
  return {
    thetaDeg: 0,
    crankAngle: 0,
    crankPin: [0, 0],
    wristPin: [0, 0],
    rodRotation: 0,
    rodSpan: 0,
    crownY: 0,
    headY: 0,
    crownYCyl: 0,
    camAngle: 0,
    motorAngle: 0,
    valves: [vp(), vp()],
  };
}

const DEG = Math.PI / 180;

/**
 * Compute every moving-part pose from the snapshot kinematic fields.
 * Allocation-free: writes into `out`.
 * @param headY world height of the fire deck (from the compression ratio)
 * @param baseRadius cam base-circle radius per valve [intake, exhaust]
 */
export function computePose(
  L: EngineLayout,
  headY: number,
  s: KinematicInput,
  baseRadius: readonly [number, number],
  out: EnginePose,
): EnginePose {
  const th = (Number.isFinite(s.thetaDeg) ? s.thetaDeg : 0) * DEG;
  const r = L.crankRadius;
  const phi = th + L.phiTdc;
  out.thetaDeg = s.thetaDeg;
  out.crankAngle = -phi;
  const px = r * Math.sin(phi);
  const py = r * Math.cos(phi);
  out.crankPin[0] = px;
  out.crankPin[1] = py;
  const disp = Number.isFinite(s.pistonDisplacement) ? s.pistonDisplacement : 0;
  const wx = L.pinOffset;
  const wy = L.wristPinTdcY - disp;
  out.wristPin[0] = wx;
  out.wristPin[1] = wy;
  const dx = wx - px, dy = wy - py;
  out.rodSpan = Math.hypot(dx, dy);
  out.rodRotation = Math.atan2(-dx, dy);
  out.crownY = wy + L.compressionHeight;
  out.headY = headY;
  out.crownYCyl = out.crownY - headY;
  out.camAngle = th / 2;
  out.motorAngle = -th * L.motor.ratio;

  for (let i = 0; i < 2; i++) {
    const v = L.valves[i];
    const vp = out.valves[i];
    const raw = i === 0 ? s.intakeLift : s.exhaustLift;
    const lift = Number.isFinite(raw) ? Math.max(0, raw) : 0;
    vp.lift = lift;
    const sa = Math.min(1, lift / v.armValve);
    const alpha = Math.asin(sa);
    vp.rockerAngle = alpha;
    const tl = v.armPushrod * sa;
    vp.tappetLift = tl;
    vp.tappetFaceY = L.camY + baseRadius[i] + tl;
    vp.pushrodBottom[0] = v.tappetX;
    vp.pushrodBottom[1] = vp.tappetFaceY + L.tappetLength;
    vp.pushrodBottom[2] = v.lobeZ;
    // Rocker cup: pivot + R_y(yaw) R_x(α) (0, 0, −a_p)
    const ca = Math.cos(alpha);
    const sy = Math.sin(v.rockerYaw), cy = Math.cos(v.rockerYaw);
    const lz = -v.armPushrod * ca;
    vp.pushrodTop[0] = v.rockerPivot[0] + lz * sy;
    vp.pushrodTop[1] = headY + v.rockerPivot[1] + v.armPushrod * sa;
    vp.pushrodTop[2] = v.rockerPivot[2] + lz * cy;
  }
  return out;
}

/**
 * Exact slider-crank for tests and previews (physics sign conventions):
 * piston displacement below TDC and rod angle β, sin β = (r sin φ − e)/L.
 */
export function sliderCrank(L: EngineLayout, thetaDeg: number): { pistonDisplacement: number; rodAngle: number } {
  const phi = thetaDeg * DEG + L.phiTdc;
  const r = L.crankRadius, l = L.rodLength, e = L.pinOffset;
  const px = r * Math.sin(phi), py = r * Math.cos(phi);
  const dx = px - e;
  const wy = py + Math.sqrt(Math.max(l * l - dx * dx, 0));
  return { pistonDisplacement: L.wristPinTdcY - wy, rodAngle: Math.asin(dx / l) };
}
