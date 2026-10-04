/**
 * ModelTEngineModel — procedural 3D model of the Ford Model T engine behind the EngineRenderModel
 * contract (render/engine-model.ts), animated ONLY from snapshot fields:
 *   thetaDeg → crankshaft, flywheel + magnets + triple gears, timing gears, camshaft (½ speed) + timer
 *              rotor, fan (belt ratio)
 *   s.cylinders[i].pistonDisplacement / intakeLift / exhaustLift (cylinder 1: the top-level fields) →
 *              piston, rod, valve, spring of cylinder i; absent → exact slider-crank at θ_i and the lift
 *              curve learned from cylinder 1 (the cam model until learned)
 *   tappets ride their lobes (cam geometry at θ/2), so lobe and tappet always touch
 *   setControls(op): sparkAdvanceDeg → timer case, throttle → carburettor butterfly
 * cutRegion(i) hands the in-cylinder visuals the housings' cut-away quadrant in cylinder i's frame.
 * The compression ratio is fixed (setCompressionRatio is a no-op, chamberShift 0).
 *
 * Frames: `root` = ROOT/world. cylinderFrames[i] = the physics L-head cylinder frame of cylinder i
 * (origin on the bore axis at the head-cavity roof), mirrored z → −z where spec.layout.mirrorZ[i]; the
 * spark plugs live there. Cutaway: block, head, pan and front cover lose the quadrant
 * {valve side, z > z(cylinder 1's front valve)}; the transmission cover loses its valve-side half
 * (flywheel magnets and coil ring visible). Repeated parts share one geometry (Object3D.clone).
 */
import * as THREE from 'three';
import type { EngineSpec } from '../../physics/core/engine-spec';
import { cylinderAngleDeg } from '../../physics/core/engine-spec';
import type { OperatingPoint } from '../../physics/core/operating-point';
import type { EngineSnapshot } from '../../physics/core/snapshot';
import { throttlePlateAngle } from '../../physics/gas-exchange/throttle';
import type { CutPlanes } from '../combustion/index';
import type { EngineRenderModel } from '../engine-model';
import type { CameraView } from '../engine/index';
import { LiftProfileLearner } from '../engine/cam-profile';
import { buildSparkPlug } from '../engine/parts-cylinder';
import { computeModelTLayout, type ModelTLayout } from './layout';
import { computeModelTPose, createModelTInput, createModelTPose, type ModelTKinematicInput, type ModelTPose } from './kinematics';
import { ModelTMaterials } from './materials';
import { FrameSet } from './frame-set';
import { buildBlock } from './parts-block';
import { buildHead } from './parts-head';
import { buildCrank, buildPiston, buildRod } from './parts-crank';
import { buildCamshaft, buildSpring, buildTappet, buildTimerCase, buildValve } from './parts-valvetrain';
import { buildBell, buildCarburettor, buildFan, buildManifolds, buildStarter } from './parts-ancillary';

const DEG = Math.PI / 180;

export class ModelTEngineModel implements EngineRenderModel {
  readonly root: THREE.Group;
  readonly cylinderFrames: readonly THREE.Object3D[];
  readonly layout: ModelTLayout;
  readonly compressionRatio: number;
  private readonly mats = new ModelTMaterials();
  private readonly crank: THREE.Group;
  private readonly rods: THREE.Object3D[] = [];
  private readonly pistons: THREE.Object3D[] = [];
  private readonly cam: THREE.Group;
  private readonly valves: THREE.Object3D[] = [];
  private readonly springs: THREE.Object3D[] = [];
  private readonly tappets: THREE.Object3D[] = [];
  private readonly timerCase: THREE.Group;
  private readonly fan: THREE.Group;
  private readonly throttle: THREE.Group;
  private readonly blockSet: { full: THREE.Group; cut: THREE.Group };
  private readonly bellSet: { full: THREE.Group; cut: THREE.Group };
  private readonly input: ModelTKinematicInput;
  private readonly _pose: ModelTPose;
  private readonly learners: [LiftProfileLearner, LiftProfileLearner] = [new LiftProfileLearner(), new LiftProfileLearner()];
  /** Per cylinder: the housings' cut-away quadrant in that cylinder's frame (see cutRegion). */
  private readonly cutPlanes: readonly CutPlanes[];
  private cutaway = true;

  constructor(spec: EngineSpec, op?: Partial<OperatingPoint>) {
    const L = (this.layout = computeModelTLayout(spec));
    const M = this.mats;
    this.compressionRatio = spec.geometry.compressionRatio;
    this.input = createModelTInput(L);
    this._pose = createModelTPose(L);
    this.root = new THREE.Group();
    this.root.name = 'engine-model-t';

    // ---- cylinder frames (physics frames; mirrored where the valve order is) ----
    const frames: THREE.Object3D[] = [];
    for (let i = 0; i < L.nCyl; i++) {
      const f = new THREE.Group();
      f.name = `cylinder-frame-${i + 1}`;
      f.position.set(0, L.cylOriginY, L.axisZ[i]);
      if (L.mirror[i]) f.scale.set(1, 1, -1);
      this.root.add(f);
      frames.push(f);
    }
    this.cylinderFrames = frames;
    // The block/head/pan/cover quadrant {cutSide·x > 0, z > cutZ} (ROOT; FrameSet at z = cutZ) in each
    // cylinder frame p_ROOT = (x, cylOriginY + y, axisZ[i] + s_z·z), s_z = −1 where mirrored.
    this.cutPlanes = frames.map((_, i): CutPlanes => [
      [L.cutSide, 0, 0, 0],
      [0, 0, L.mirror[i] ? -1 : 1, L.cutZ - L.axisZ[i]],
    ]);

    // ---- static housings: block + head (+ pan, cover, manifolds) and the transmission cover ----
    const set = new FrameSet([0, 0, L.cutZ]);
    buildBlock(L, set);
    buildHead(L, set);
    buildManifolds(L, set);
    const fanParts = buildFan(L, M, set);
    buildStarter(L, set);
    const quad = { kind: 'quadrant' as const, side: L.cutSide };
    const blk = set.build(M, 'block', quad);
    this.blockSet = blk;
    this.root.add(blk.group);
    const bset = new FrameSet([0, 0, L.bellCutZ]);
    buildBell(L, bset);
    const bell = bset.build(M, 'bell', quad);
    this.bellSet = bell;
    this.root.add(bell.group);

    // ---- moving parts ----
    this.crank = buildCrank(L, M);
    this.root.add(this.crank);
    const rod0 = buildRod(L, M);
    const piston0 = buildPiston(L, M);
    for (let i = 0; i < L.nCyl; i++) {
      const rd = i === 0 ? rod0 : rod0.clone();
      rd.name = `cyl${i + 1}-rod`;
      const ps = i === 0 ? piston0 : piston0.clone();
      ps.name = `cyl${i + 1}-piston`;
      this.rods.push(rd);
      this.pistons.push(ps);
      this.root.add(rd, ps);
    }
    const cs = buildCamshaft(L, M);
    this.cam = cs.group;
    this.root.add(this.cam);
    const valveI = buildValve(L, M, M.m.valve);
    const valveE = buildValve(L, M, M.m.valveExhaust);
    const spring0 = buildSpring(L, M);
    const tappet0 = buildTappet(L, M);
    for (const v of L.valves) {
      const vg = (v.kind === 'intake' ? valveI : valveE).clone();
      vg.name = `valve-${v.index}`;
      vg.position.set(v.x, L.seatLineY, v.z);
      const sp = spring0.clone();
      sp.name = `spring-${v.index}`;
      sp.position.set(v.x, L.spring.topY, v.z);
      const tp = tappet0.clone();
      tp.name = `tappet-${v.index}`;
      tp.position.set(v.x, L.tappet.faceRestY, v.z);
      this.valves.push(vg);
      this.springs.push(sp);
      this.tappets.push(tp);
      this.root.add(vg, sp, tp);
    }
    this.timerCase = buildTimerCase(L, M);
    this.root.add(this.timerCase);
    this.fan = fanParts.fan;
    this.root.add(this.fan, fanParts.belt);
    const carb = buildCarburettor(L, M);
    this.throttle = carb.throttle;
    this.root.add(carb.body, carb.throttle);

    // spark plugs in the cylinder frames
    const plug0 = buildSparkPlug({ plug: L.plug, spec }, M);
    frames.forEach((f, i) => {
      const p = i === 0 ? plug0 : plug0.clone();
      p.name = `spark-plug-${i + 1}`;
      f.add(p);
    });

    this.root.traverse((o) => {
      if ((o as THREE.Mesh).isMesh) {
        const section = o.name.includes('-section-') || o.name === 'section';
        o.castShadow = !section;
        o.receiveShadow = !section;
      }
    });
    this.setControls({ sparkAdvanceDeg: op?.sparkAdvanceDeg ?? 0, throttle: op?.throttle ?? 0 });
    this.setCutaway(true);
    this.applyPose();
  }

  /** Latest computed part poses (read-only view). */
  get pose(): Readonly<ModelTPose> {
    return this._pose;
  }

  get isCutaway(): boolean {
    return this.cutaway;
  }

  update(s: EngineSnapshot): void {
    const L = this.layout;
    const inp = this.input;
    const th = s.thetaDeg;
    inp.thetaDeg = th;
    const cyl = s.cylinders;
    if (!cyl) this.learn(th, s.intakeLift, s.exhaustLift);
    const li = this.learners[0], le = this.learners[1];
    const readyI = li.isReady(), readyE = le.isReady();
    for (let i = 0; i < L.nCyl; i++) {
      const c = cyl ? cyl[i] : undefined;
      if (c) {
        inp.pistonDisplacement[i] = c.pistonDisplacement;
        inp.intakeLift[i] = c.intakeLift;
        inp.exhaustLift[i] = c.exhaustLift;
      } else if (i === 0) {
        inp.pistonDisplacement[i] = s.pistonDisplacement;
        inp.intakeLift[i] = s.intakeLift;
        inp.exhaustLift[i] = s.exhaustLift;
      } else {
        // identical cylinders: exact geometry at θ_i, lifts from the curve learned on cylinder 1
        const ti = cylinderAngleDeg(L.spec, i, th);
        inp.pistonDisplacement[i] = NaN;
        inp.intakeLift[i] = readyI ? li.sample(ti) : NaN;
        inp.exhaustLift[i] = readyE ? le.sample(ti) : NaN;
      }
    }
    this.applyPose();
  }

  setControls(op: Partial<OperatingPoint>): void {
    if (op.sparkAdvanceDeg !== undefined && Number.isFinite(op.sparkAdvanceDeg)) {
      // the case turns against the rotor by half the crank-angle advance (timer at cam speed)
      this.timerCase.rotation.z = (-op.sparkAdvanceDeg * DEG) / 2;
    }
    if (op.throttle !== undefined && Number.isFinite(op.throttle)) {
      this.throttle.rotation.z = throttlePlateAngle(op.throttle, this.layout.carb.plateClosedAngle);
    }
  }

  setCompressionRatio(_cr: number): void {
    // fixed head: nothing moves with the compression ratio
  }

  chamberShift(_cr: number, _previous: number): number {
    return 0;
  }

  setCutaway(on: boolean): void {
    this.cutaway = on;
    this.blockSet.full.visible = !on;
    this.blockSet.cut.visible = on;
    this.bellSet.full.visible = !on;
    this.bellSet.cut.visible = on;
  }

  /**
   * The cut-away region in cylinder `cylinder`'s frame for CombustionVisuals.setCutRegion: the removed
   * quadrant of block and head, {p : n·p > d for both planes}; null while the cutaway is off (and for an
   * index outside the engine). The transmission-cover half section is far from every chamber and is not
   * included. The returned planes are shared; do not modify them.
   */
  cutRegion(cylinder: number): CutPlanes | null {
    return this.cutaway ? (this.cutPlanes[cylinder] ?? null) : null;
  }

  cameraView(framing: 'engine' | 'chamber', cylinder = 0): CameraView {
    const L = this.layout;
    const { min, max } = L.bounds;
    if (framing === 'chamber') {
      const k = Math.min(Math.max(Math.round(Number.isFinite(cylinder) ? cylinder : 0), 0), L.nCyl - 1);
      const vz = L.valves.filter((v) => v.cyl === k).map((v) => v.z);
      const zc = k === 0 ? L.cutZ - 0.01 : Math.max(...vz);
      const target = new THREE.Vector3(L.cutSide * 0.045, L.deckY + 0.004, zc);
      const position = target.clone().add(new THREE.Vector3(L.cutSide * 0.22, 0.13, 0.26));
      return { position, target, fov: 30, near: 0.005, far: 20 };
    }
    const c = new THREE.Vector3((min[0] + max[0]) / 2, (min[1] + max[1]) / 2 + 0.02, (min[2] + max[2]) / 2 + 0.06);
    const radius = 0.5 * Math.hypot(max[0] - min[0], max[1] - min[1], max[2] - min[2]);
    const fov = 34;
    const dist = (radius / Math.sin((fov * DEG) / 2)) * 0.62;
    const dir = new THREE.Vector3(L.cutSide * 0.62, 0.42, 0.66).normalize();
    return { position: c.clone().addScaledVector(dir, dist), target: c, fov, near: 0.01, far: 30 };
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
    this.learners[0].add(thetaDeg, li);
    this.learners[1].add(thetaDeg, le);
  }

  private applyPose(): void {
    const L = this.layout;
    const P = computeModelTPose(L, this.input, this._pose);
    this.crank.rotation.z = P.crankAngle;
    this.cam.rotation.z = P.camAngle;
    this.fan.rotation.z = P.fanAngle;
    for (let i = 0; i < L.nCyl; i++) {
      const c = P.cylinders[i];
      const z = L.axisZ[i];
      this.rods[i].position.set(c.crankPin[0], c.crankPin[1], z);
      this.rods[i].rotation.z = c.rodRotation;
      this.pistons[i].position.set(0, c.wristPinY, z);
    }
    for (let k = 0; k < L.valves.length; k++) {
      const vp = P.valves[k];
      this.valves[k].position.y = L.seatLineY + vp.lift;
      this.springs[k].scale.y = vp.springScale;
      this.tappets[k].position.y = vp.tappetFaceY;
    }
  }
}
