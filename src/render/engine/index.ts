/**
 * EngineModel — procedural, physically proportioned 3D model of the CFR F-1
 * mechanism, animated ONLY from snapshot fields:
 *   thetaDeg → crankshaft, flywheel, timing gears, camshaft (½ speed), motor pulley
 *   pistonDisplacement → piston (+ rod from crank pin to wrist pin)
 *   intakeLift / exhaustLift → valves, springs, rockers, pushrods, tappets
 *   setCompressionRatio → cylinder + head height (h_TDC = S/(CR−1)), worm gear
 *
 * Frames: `root` = world/ROOT frame (crank axis = z through the origin, +y up).
 * `cylinderFrame` is translated to the fire deck so its local coordinates are
 * the physics cylinder frame (engine-spec.ts); in-cylinder visuals go there.
 *
 * Cutaway: the crankcase is half-sectioned (z > 0 removed); the cylinder,
 * head and rocker cover get a quarter section (the quadrant {cutSide·x > 0,
 * z > 0} removed) so both the valve line (z = 0) and the spark-plug axis
 * (x = 0) are cut open. Section caps are drawn from back faces (geometry.ts),
 * so no renderer clipping state or stencil buffer is needed.
 *
 * The cam lobes are generated from the lift curve the physics produces:
 * a model curve from the ValveSpec first, replaced by the curve learned from
 * the snapshot stream (cam-profile.ts) so lobe and tappet stay in contact.
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { computeLayout, type EngineLayout } from './layout';
import { computePose, createPose, type EnginePose, type KinematicInput } from './kinematics';
import { EngineMaterials } from './materials';
import { StaticSet } from './static-set';
import { buildCrank } from './parts-crank';
import { buildPiston, buildRod } from './parts-rod-piston';
import { buildCylinder, type ValveParts } from './parts-cylinder';
import { buildCrankcase } from './parts-crankcase';
import {
  buildCamshaft, buildMotorPulley, buildPushrod, buildTappet, buildWormDrive, lobeGeometry,
  type CamParts, type PushrodParts, type WormParts,
} from './parts-valvetrain';
import { LiftProfileLearner, buildLobeProfile, modelValveLift, type LobeProfile } from './cam-profile';

export { computeLayout, type EngineLayout } from './layout';
export { computePose, createPose, sliderCrank, type EnginePose, type KinematicInput } from './kinematics';

const Y_UP = new THREE.Vector3(0, 1, 0);

export interface CameraView {
  position: THREE.Vector3;
  target: THREE.Vector3;
  /** Suggested vertical field of view, degrees. */
  fov: number;
  /** Suggested near/far planes, m. */
  near: number;
  far: number;
}

/**
 * Recommended camera framing (ROOT/world coordinates, assuming `root` is
 * placed at the scene origin). 'engine' frames the whole machine; 'chamber'
 * is a close-up of the cut-open combustion chamber and spark plug.
 */
export function recommendedCameraView(spec: EngineSpec, compressionRatio: number, framing: 'engine' | 'chamber' = 'engine'): CameraView {
  const L = computeLayout(spec);
  const hy = L.headY(compressionRatio);
  const a = L.cutSide; // look into the cut-away quadrant
  if (framing === 'chamber') {
    // look into the cut-away quadrant, slightly from above, at the chamber and the spark gap
    const gc = L.plug.gapCenter;
    return {
      target: new THREE.Vector3(gc[0] * 0.4, hy - 0.018, gc[2] * 0.3),
      position: new THREE.Vector3(a * 0.27, hy + 0.07, 0.23),
      fov: 32,
      near: 0.005,
      far: 20,
    };
  }
  return {
    target: new THREE.Vector3(a * 0.04, hy * 0.45, -0.03),
    position: new THREE.Vector3(a * 0.75, hy * 0.45 + 0.36, 1.3),
    fov: 34,
    near: 0.01,
    far: 30,
  };
}

export class EngineModel {
  readonly root: THREE.Group;
  readonly cylinderFrame: THREE.Object3D;
  readonly layout: EngineLayout;
  private readonly mats = new EngineMaterials();
  private readonly crank: THREE.Group;
  private readonly rod: THREE.Group;
  private readonly piston: THREE.Group;
  private readonly cam: CamParts;
  private readonly tappets: [THREE.Group, THREE.Group];
  private readonly pushrods: [PushrodParts, PushrodParts];
  private readonly valves: [ValveParts, ValveParts];
  private readonly rockers: [THREE.Group, THREE.Group];
  private readonly worm: WormParts;
  private readonly motorPulley: THREE.Group;
  private readonly staticRoot: { full: THREE.Group; cut: THREE.Group };
  private readonly staticCyl: { full: THREE.Group; cut: THREE.Group };
  private readonly _pose: EnginePose = createPose();
  private readonly learners: [LiftProfileLearner, LiftProfileLearner] = [new LiftProfileLearner(), new LiftProfileLearner()];
  private readonly builtStage: [number, number] = [0, 0];
  private readonly builtResets: [number, number] = [0, 0];
  private readonly baseRadius: [number, number];
  private readonly input: KinematicInput = { thetaDeg: 0, pistonDisplacement: 0, intakeLift: 0, exhaustLift: 0 };
  private readonly tmpDir = new THREE.Vector3();
  private readonly tmpDir2 = new THREE.Vector3();
  private cr: number;
  private headY: number;
  private cutaway = true;

  constructor(spec: EngineSpec, compressionRatio: number) {
    const L = (this.layout = computeLayout(spec));
    const M = this.mats;
    this.root = new THREE.Group();
    this.root.name = 'engine-model';
    const cyl = new THREE.Group();
    cyl.name = 'cylinder-frame';
    this.cylinderFrame = cyl;
    this.root.add(cyl);

    // ---- static housings ----
    const rootSet = new StaticSet();
    buildCrankcase(L, rootSet);
    this.staticRoot = rootSet.build(M, 'crankcase');
    this.root.add(this.staticRoot.full, this.staticRoot.cut);

    const cylSet = new StaticSet();
    const cp = buildCylinder(L, M, cylSet);
    this.staticCyl = cylSet.build(M, 'cylinder', { kind: 'quadrant', side: L.cutSide });
    cyl.add(this.staticCyl.full, this.staticCyl.cut);
    this.valves = cp.valves;
    this.rockers = cp.rockers;
    for (const v of cp.valves) cyl.add(v.moving, v.spring);
    for (const r of cp.rockers) cyl.add(r);
    cyl.add(cp.plug);

    // ---- rotating / reciprocating ----
    this.crank = buildCrank(L, M);
    this.rod = buildRod(L, M);
    this.piston = buildPiston(L, M);
    this.root.add(this.crank, this.rod, this.piston);

    // cam lobes from the model lift curve (replaced by the learned curve later)
    const profiles = L.valves.map((v) =>
      buildLobeProfile((th) => (modelValveLift(v.spec, th) * v.armPushrod) / v.armValve, L.camBaseRadius),
    ) as [LobeProfile, LobeProfile];
    this.baseRadius = [profiles[0].baseRadius, profiles[1].baseRadius];
    this.cam = buildCamshaft(L, M, profiles);
    this.root.add(this.cam.group);
    const faceR = Math.max(L.tappetFaceRadius, 1.3 * Math.max(profiles[0].maxContactOffset, profiles[1].maxContactOffset) + 0.002);
    this.tappets = [buildTappet(L, M, faceR), buildTappet(L, M, faceR)];
    this.root.add(...this.tappets);

    // telescoping pushrods sized for the full CR travel
    const spans = L.valves.map((v, i) => {
      const bot = new THREE.Vector3(v.tappetX, L.camY + this.baseRadius[i] + L.tappetLength, v.lobeZ);
      const top = (hy: number) => new THREE.Vector3(v.pushrodTopRest[0], hy + v.pushrodTopRest[1], v.pushrodTopRest[2]);
      return [top(L.headYMin).distanceTo(bot), top(L.headYMax).distanceTo(bot)] as const;
    });
    this.pushrods = spans.map(([smin, smax]) => {
      const lower = 0.55 * smin;
      return buildPushrod(L, M, lower, smax - lower + 0.04);
    }) as [PushrodParts, PushrodParts];
    for (const p of this.pushrods) this.root.add(p.lower, p.upper);

    this.worm = buildWormDrive(L, M);
    this.root.add(this.worm.wheel, this.worm.wheelCut, this.worm.worm);
    this.motorPulley = buildMotorPulley(L, M);
    this.root.add(this.motorPulley);

    this.root.traverse((o) => {
      if ((o as THREE.Mesh).isMesh) {
        // section caps are drawn at the far interior surface: no shadows (avoids self-shadow stripes)
        const section = o.name.includes('-section-') || o.name === 'section';
        o.castShadow = !section;
        o.receiveShadow = !section;
      }
    });

    this.cr = compressionRatio;
    this.headY = L.headY(compressionRatio);
    this.setCompressionRatio(compressionRatio);
    this.setCutaway(true);
    this.applyPose();
  }

  /** Latest computed part poses (read-only view). */
  get pose(): Readonly<EnginePose> {
    return this._pose;
  }

  get compressionRatio(): number {
    return this.cr;
  }

  setCompressionRatio(cr: number): void {
    const L = this.layout;
    const c = Math.min(Math.max(Number.isFinite(cr) ? cr : L.spec.geometry.compressionRatio, L.crRange[0]), L.crRange[1]);
    this.cr = c;
    this.headY = L.headY(c);
    this.cylinderFrame.position.set(0, this.headY, 0);
    this.cylinderFrame.updateMatrix();
    // worm wheel is threaded onto the cylinder spigot: rotation ∝ height change
    const omega = ((this.headY - L.headYMin) / L.threadPitch) * 2 * Math.PI;
    this.worm.wheel.rotation.y = omega;
    const p = this.worm.toothPitch;
    this.worm.wheelCut.rotation.y = ((omega % p) + p) % p;
    this.worm.worm.rotation.z = -omega * L.wormWheel.teeth;
    this.applyPose();
  }

  update(s: EngineSnapshot): void {
    const inp = this.input;
    inp.thetaDeg = s.thetaDeg;
    inp.pistonDisplacement = s.pistonDisplacement;
    inp.intakeLift = s.intakeLift;
    inp.exhaustLift = s.exhaustLift;
    this.learn(s.thetaDeg, s.intakeLift, s.exhaustLift);
    this.applyPose();
  }

  setCutaway(on: boolean): void {
    this.cutaway = on;
    this.staticRoot.full.visible = !on;
    this.staticRoot.cut.visible = on;
    this.staticCyl.full.visible = !on;
    this.staticCyl.cut.visible = on;
    this.worm.wheel.visible = !on;
    this.worm.wheelCut.visible = on;
  }

  get isCutaway(): boolean {
    return this.cutaway;
  }

  dispose(): void {
    const geoms = new Set<THREE.BufferGeometry>();
    this.root.traverse((o) => {
      const m = o as THREE.Mesh;
      if (m.isMesh && m.geometry) geoms.add(m.geometry);
    });
    for (const g of geoms) g.dispose();
    this.mats.dispose();
    this.root.removeFromParent();
  }

  // ------------------------------------------------------------------

  private learn(thetaDeg: number, li: number, le: number): void {
    const L = this.layout;
    for (let i = 0; i < 2; i++) {
      const ln = this.learners[i];
      ln.add(thetaDeg, i === 0 ? li : le);
      if (ln.resets !== this.builtResets[i]) {
        this.builtResets[i] = ln.resets;
        this.builtStage[i] = 0;
      }
      const st = ln.stage();
      if (st > this.builtStage[i]) {
        this.builtStage[i] = st;
        const v = L.valves[i];
        const prof = buildLobeProfile((th) => (ln.sample(th) * v.armPushrod) / v.armValve, L.camBaseRadius);
        this.baseRadius[i] = prof.baseRadius;
        const mesh = this.cam.lobes[i];
        mesh.geometry.dispose();
        mesh.geometry = lobeGeometry(prof, L.lobeWidth, v.lobeZ);
      }
    }
  }

  private applyPose(): void {
    const L = this.layout;
    const P = computePose(L, this.headY, this.input, this.baseRadius, this._pose);
    this.crank.rotation.z = P.crankAngle;
    this.rod.position.set(P.crankPin[0], P.crankPin[1], 0);
    this.rod.rotation.z = P.rodRotation;
    this.piston.position.set(P.wristPin[0], P.wristPin[1], 0);
    this.cam.group.rotation.z = P.camAngle;
    this.motorPulley.rotation.z = P.motorAngle;
    for (let i = 0; i < 2; i++) {
      const v = L.valves[i];
      const vp = P.valves[i];
      const parts = this.valves[i];
      parts.moving.position.y = -vp.lift;
      parts.spring.scale.y = Math.max(0.2, (v.springInstalledLength - vp.lift) / v.springInstalledLength);
      this.rockers[i].rotation.x = vp.rockerAngle;
      this.tappets[i].position.set(v.tappetX, vp.tappetFaceY, v.lobeZ);
      const pr = this.pushrods[i];
      const b = vp.pushrodBottom, t = vp.pushrodTop;
      const d = this.tmpDir.set(t[0] - b[0], t[1] - b[1], t[2] - b[2]).normalize();
      pr.lower.position.set(b[0], b[1], b[2]);
      pr.lower.quaternion.setFromUnitVectors(Y_UP, d);
      pr.upper.position.set(t[0], t[1], t[2]);
      pr.upper.quaternion.setFromUnitVectors(Y_UP, this.tmpDir2.copy(d).negate());
    }
  }
}
