/**
 * Procedural layout of the CFR-style engine mechanism: every dimension and
 * placement the 3D model uses, derived from the EngineSpec (bore, stroke,
 * rod length, compression height, valve/plug geometry, CR range). Pure data,
 * no three.js, so the kinematics and clearance logic can be unit-tested.
 *
 * Frames (see src/physics/core/engine-spec.ts):
 *  - ROOT/world: +y up, crank axis = z axis through the origin, cylinder axis = y axis.
 *    The section plane of the cutaway is z = 0; the removed half is z > 0
 *    (the default camera looks from +z).
 *  - CYL: origin at the centre of the head fire-deck face, +y toward the head.
 *    It is a pure translation of ROOT by (0, headY, 0), where headY depends on
 *    the compression ratio (the CFR raises the whole cylinder + head).
 *
 * Proportions that the spec does not define (journal sizes, wall thicknesses,
 * rocker ratio, cam base circle, flywheel size, …) are render-only choices,
 * scaled from the bore and stroke so the model stays proportioned if the spec
 * changes. They never feed back into the physics.
 */
import type { EngineSpec, ValveSpec } from '../../physics/core/engine-spec';

export type Vec3 = [number, number, number];

/** Rocker ratio (valve-side arm / pushrod-side arm). Render-only choice. UNVERIFIED for the CFR. */
export const ROCKER_RATIO = 1.4;

export interface ValveLayout {
  key: 'intake' | 'exhaust';
  spec: ValveSpec;
  /** Valve axis position on the fire deck, cylinder frame. */
  x: number;
  z: number;
  headRadius: number;
  seatInnerRadius: number;
  stemRadius: number;
  /** Cylindrical margin height of the valve head, m. */
  margin: number;
  /** Vertical extent of the conical seat face, m. */
  faceHeight: number;
  /** Height of the tulip (face top to stem), m. */
  tulipHeight: number;
  /** Radius of the pocket machined in the head for the seat insert. */
  pocketRadius: number;
  /** Radial clearance between valve margin and insert lip. */
  lipClearance: number;
  /** Seat insert height (pocket depth), m. The port starts at this height. */
  insertHeight: number;
  /** Stem tip height at zero lift (cylinder frame). */
  tipY: number;
  guideTopY: number;
  guideBottomY: number;
  guideOuterRadius: number;
  springSeatY: number;
  /** Spring length at zero lift. */
  springInstalledLength: number;
  springMeanRadius: number;
  springWireRadius: number;
  springTurns: number;
  retainerThickness: number;
  /** Port leaves toward +x (1) or −x (−1). */
  portDir: 1 | -1;
  portExitRadius: number;
  /** Height of the horizontal part of the port centreline. */
  portY: number;
  /** Rocker pivot (cylinder frame). At zero lift the valve pad, pivot and pushrod cup lie on y = pivot[1]. */
  rockerPivot: Vec3;
  /** Unit plan-view direction from the pushrod end toward the valve end. */
  rockerDir: [number, number];
  /** Yaw of the rocker (rotation about +y taking +z to rockerDir). */
  rockerYaw: number;
  /** Arm lengths pivot→valve pad and pivot→pushrod cup. a_v / a_p = ROCKER_RATIO. */
  armValve: number;
  armPushrod: number;
  /** Pushrod cup on the rocker at zero lift (cylinder frame). */
  pushrodTopRest: Vec3;
  /** Plan position (world x, z) of the tappet / cam lobe. */
  tappetX: number;
  lobeZ: number;
}

export interface EngineLayout {
  spec: EngineSpec;
  /**
   * Cutaway of the cylinder assembly removes the quadrant {cutSide·x > 0, z > 0}
   * (quarter section through the valve line and the spark-plug axis);
   * the crankcase is half-sectioned (z > 0 removed).
   */
  cutSide: 1 | -1;
  /** Crank-pin axis angle at TDC, rad (≠ 0 only with a pin offset). */
  phiTdc: number;
  /** Exact piston travel TDC→BDC, m (= stroke without pin offset). */
  pistonTravel: number;
  // ---- basic geometry ----
  bore: number;
  stroke: number;
  crankRadius: number;
  rodLength: number;
  pinOffset: number;
  compressionHeight: number;
  crRange: [number, number];

  // ---- crankshaft ----
  mainJournalRadius: number;
  crankPinRadius: number;
  crankPinLength: number;
  webThickness: number;
  /** z of the centre of each web (±). */
  webZ: number;
  webPinBossRadius: number;
  webJournalBossRadius: number;
  counterweightRadius: number;
  /** Half-angle of the counterweight sector, rad. */
  counterweightHalfAngle: number;
  frontJournalEndZ: number;
  rearJournalEndZ: number;

  // ---- connecting rod (rod frame: origin at crank pin, +y toward the wrist pin) ----
  rodBigEndWidth: number;
  rodBigEndOuterRadius: number;
  rodBearingThickness: number;
  rodBoltOffset: number;
  rodBoltRadius: number;
  rodCapDepth: number;
  rodSmallEndWidth: number;
  rodSmallEndOuterRadius: number;
  rodBushingThickness: number;
  rodShankWidthBig: number;
  rodShankWidthSmall: number;
  rodFlangeThickness: number;
  rodWebThickness: number;

  // ---- piston (piston frame: origin at the wrist-pin axis) ----
  pistonRadius: number;
  pistonCrownThickness: number;
  pistonSkirtBottom: number; // y below the pin axis (negative)
  pistonWallThickness: number;
  /** Ring grooves: top y (below crown, positive distance), height, radial depth. */
  ringGrooves: { top: number; height: number; depth: number; oil: boolean }[];
  wristPinRadius: number;
  wristPinBore: number;
  wristPinLength: number;
  pinBossRadius: number;

  // ---- vertical placement ----
  /** Wrist-pin height at TDC, world. */
  wristPinTdcY: number;
  crownTdcY: number;

  // ---- cylinder casting (cylinder frame) ----
  boreRadius: number;
  linerOuterRadius: number;
  jacketOuterRadius: number;
  jacketDepth: number;
  cylinderLength: number;
  spigotRadius: number;
  jacketCavity: { rIn: number; rOut: number; yTop: number; yBottom: number };
  threadPitch: number;

  // ---- head (cylinder frame) ----
  headThickness: number;
  headXHalf: number;
  headZMin: number;
  headZMax: number;
  valves: [ValveLayout, ValveLayout];
  pushrodRadius: number;
  pushrodTubeRadius: number;
  pushrodPassageRadius: number;
  rockerCover: { yBottom: number; yTop: number; wall: number };

  // ---- cam / tappets (world) ----
  camX: number;
  camY: number;
  camShaftRadius: number;
  lobeWidth: number;
  /** Nominal base-circle radius. The lobe builder may enlarge it for convexity. */
  camBaseRadius: number;
  tappetRadius: number;
  tappetFaceRadius: number;
  tappetLength: number;
  tappetGuide: { yBottom: number; yTop: number };
  camZFront: number;
  camZBack: number;
  gearZ: number;
  gearWidth: number;
  crankGearRadius: number;
  camGearRadius: number;
  crankGearTeeth: number;
  camGearTeeth: number;

  // ---- crankcase (world) ----
  crankcase: {
    xHalfOuter: number;
    xHalfInner: number;
    yBase: number;
    yFloorInner: number;
    yDeckBottom: number;
    yDeckTop: number;
    zFrontOuter: number;
    zFrontInner: number;
    zBackOuter: number;
    zBackInner: number;
    wall: number;
    baseFlangeHalf: number;
    pedestalZ: number;
    pedestalThickness: number;
  };

  // ---- flywheel + drive (world) ----
  flywheelRadius: number;
  flywheelWidth: number;
  flywheelZ: number;
  flywheelGrooveRadius: number;
  bedplateTop: number;
  motor: { x: number; y: number; radius: number; length: number; zCenter: number; pulleyRadius: number; ratio: number };

  // ---- CR drive (world) ----
  wormWheel: { innerRadius: number; outerRadius: number; pitchRadius: number; thickness: number; yBottom: number; teeth: number };
  worm: { x: number; y: number; radius: number; zFront: number; zBack: number };

  // ---- spark plug (cylinder frame) ----
  plug: {
    gapCenter: Vec3;
    axis: Vec3;
    /** Unit vector (⊥ axis) toward which the ground-electrode strut sits. */
    side: Vec3;
    /** Distance from the gap centre (along −axis) to where the plug enters the chamber wall (shell nose). */
    noseDistance: number;
    /** Distance from the gap centre (along −axis) to the seat face on the boss outside the casting. */
    seatDistance: number;
    threadRadius: number;
    bossRadius: number;
  };

  // ---- CR-dependent placement ----
  /**
   * Head-face height above the crown TDC position for a given CR:
   * h_TDC = travel/(CR−1) − V_crevice/A (same convention as the physics
   * kinematics: the crevice is part of the clearance volume).
   */
  clearanceAtTdc(cr: number): number;
  /** World y of the cylinder-frame origin (fire-deck face) for a given CR. */
  headY(cr: number): number;
  headYMin: number;
  headYMax: number;
}

/**
 * Disc clearance height at TDC: V_c = V_d/(CR−1) is the TOTAL clearance
 * volume (crevice included), h = (V_c − V_crevice)/A = travel/(CR−1) − V_crevice/A.
 */
export function clearanceAtTdc(travel: number, cr: number, creviceVolume = 0, bore = 1): number {
  const A = (Math.PI * bore * bore) / 4;
  return travel / (Math.max(cr, 1.0001) - 1) - Math.max(creviceVolume, 0) / A;
}

function clamp(v: number, lo: number, hi: number): number {
  return v < lo ? lo : v > hi ? hi : v;
}

/** Distance along −axis from `p` until the point leaves the casting envelope. */
function exitDistance(
  p: Vec3, axis: Vec3,
  inside: (x: number, y: number, z: number) => boolean,
  start: number, max: number,
): number {
  const step = 0.0002;
  let d = start;
  while (d < max) {
    const x = p[0] - axis[0] * d, y = p[1] - axis[1] * d, z = p[2] - axis[2] * d;
    if (!inside(x, y, z)) return d;
    d += step;
  }
  return max;
}

function normalize3(v: Vec3): Vec3 {
  const l = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / l, v[1] / l, v[2] / l];
}

function cross3(a: Vec3, b: Vec3): Vec3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

export function computeLayout(spec: EngineSpec): EngineLayout {
  const g = spec.geometry;
  const B = g.bore;
  const S = g.stroke;
  const r = S / 2;
  const L = g.conRodLength;
  const e = g.pinOffset;
  const CH = g.compressionHeight;
  const crRange: [number, number] = [
    Math.max(1.5, Math.min(g.compressionRatioRange[0], g.compressionRatioRange[1])),
    Math.max(g.compressionRatioRange[0], g.compressionRatioRange[1]),
  ];

  // ---- crankshaft ----
  const mainJournalRadius = 0.34 * B;
  // Pin must not overlap the journal (so the throw reads as a throw).
  const crankPinRadius = Math.max(0.2 * B, Math.min(0.3 * B, r - mainJournalRadius - 0.003));
  const crankPinLength = 0.46 * B;
  const webThickness = 0.26 * B;
  const webZ = crankPinLength / 2 + webThickness / 2;
  const webPinBossRadius = crankPinRadius + 0.15 * B;
  const webJournalBossRadius = mainJournalRadius + 0.14 * B;
  const counterweightRadius = r + 0.25 * B;

  // ---- rod ----
  const rodBearingThickness = 0.02 * B;
  const rodBigEndOuterRadius = crankPinRadius + rodBearingThickness + 0.13 * B;
  const rodBoltRadius = 0.06 * B;
  const rodBoltOffset = crankPinRadius + rodBearingThickness + rodBoltRadius + 0.035 * B;
  const rodCapDepth = crankPinRadius + rodBearingThickness + 0.16 * B;
  const rodBigEndWidth = crankPinLength - 0.002;
  const wristPinRadius = 0.14 * B;
  const rodBushingThickness = 0.02 * B;
  const rodSmallEndOuterRadius = wristPinRadius + rodBushingThickness + 0.1 * B;
  const rodSmallEndWidth = 0.36 * B;

  // ---- piston ----
  const pistonRadius = B / 2 - 0.00012;
  const ringH = 0.03 * B;
  const land = 0.05 * B;
  const topLand = 0.08 * B;
  const grooveDepth = 0.045 * B;
  // Five rings (CFR: [KW17] Table 2, [W850]) if the compression height allows, else four.
  const wristPinRadius0 = 0.14 * B;
  const nComp = topLand + 4 * (ringH + land) + 2 * ringH + wristPinRadius0 + 0.004 < CH ? 4 : 3;
  const ringGrooves: { top: number; height: number; depth: number; oil: boolean }[] = [];
  for (let i = 0; i < nComp; i++) ringGrooves.push({ top: topLand + i * (ringH + land), height: ringH, depth: grooveDepth, oil: false });
  ringGrooves.push({ top: topLand + nComp * (ringH + land), height: 2 * ringH, depth: grooveDepth, oil: true });

  // ---- vertical placement ----
  const wristPinTdcY = Math.sqrt(Math.max((L + r) ** 2 - e * e, 0));
  const crownTdcY = wristPinTdcY + CH;
  const phiTdc = Math.asin(Math.max(-1, Math.min(1, e / (L + r))));
  const pistonTravel = wristPinTdcY - Math.sqrt(Math.max((L - r) ** 2 - e * e, 0));
  const clr = (cr: number) => clearanceAtTdc(pistonTravel, cr, g.creviceVolume, B);
  const headY = (cr: number) => crownTdcY + clr(cr);
  const headYMin = headY(crRange[1]);
  const headYMax = headY(crRange[0]);
  const travel = headYMax - headYMin;

  // ---- cylinder casting ----
  const boreRadius = B / 2;
  const linerOuterRadius = boreRadius + 0.1 * B;
  const jacketOuterRadius = boreRadius + 0.36 * B;
  const jacketDepth = Math.max(1.05 * S, CH + 0.4 * S);
  const spigotRadius = boreRadius + 0.15 * B;
  const threadPitch = 0.004;

  // ---- worm wheel sits on the crankcase deck, around the spigot ----
  const wormWheelThickness = 0.28 * B;
  const wormWheelInner = spigotRadius + 0.002;
  const wormWheelOuter = spigotRadius + 0.32 * B;
  const deckGap = 0.004;
  const yDeckTop = headYMin - jacketDepth - wormWheelThickness - deckGap;
  const deckThickness = 0.24 * B;
  const yDeckBottom = yDeckTop - deckThickness;
  // The spigot must stay engaged in the deck over the full CR travel.
  const cylinderLength = jacketDepth + travel + wormWheelThickness + deckGap + deckThickness + 0.006;

  // ---- piston skirt: total piston height ≈ 1.46 B (CFR: 120.7 mm, [Choi18] Table 3) ----
  const pistonSkirtBottom = Math.min(-(0.4 * B), CH - 1.46 * B);

  // ---- head ----
  const headThickness = 0.9 * B;
  const headXHalf = jacketOuterRadius + 0.1 * B;

  // ---- valves ----
  const pushrodRadius = 0.06 * B;
  const pushrodTubeRadius = 0.085 * B;
  const pushrodPassageRadius = pushrodTubeRadius + 0.004;
  const clearOfCylinder = Math.max(jacketOuterRadius, wormWheelOuter) + pushrodPassageRadius + 0.012;

  const valveSpecs: [ValveSpec, ValveSpec] = [spec.intakeValve, spec.exhaustValve];
  const keys: ['intake', 'exhaust'] = ['intake', 'exhaust'];
  const camX = (valveSpecs[0].position[0] + valveSpecs[1].position[0]) / 2;

  const guideBoss = 0.07 * B;
  const retainerThickness = 0.06 * B;
  const tipProtrusion = 0.05 * B;

  // Pushrod plan positions: behind the cylinder (−z), one lane per valve.
  const lane = Math.max(0.36 * B, 2 * pushrodPassageRadius + 0.012);
  const zRocker0 = -clearOfCylinder;

  const valveGap = Math.hypot(valveSpecs[0].position[0] - valveSpecs[1].position[0], valveSpecs[0].position[1] - valveSpecs[1].position[1]);
  const valves = valveSpecs.map((vs, i): ValveLayout => {
    const headRadius = vs.headDiameter / 2;
    const seatInnerRadius = Math.min(vs.seatInnerDiameter / 2, headRadius - 0.0005);
    const stemRadius = vs.stemDiameter / 2;
    const margin = 0.04 * vs.headDiameter;
    const faceHeight = (headRadius - seatInnerRadius) * Math.tan(vs.seatAngle);
    const tulipHeight = 0.32 * vs.headDiameter;
    const lipClearance = 0.00025;
    // seat-insert pocket; kept clear of the neighbouring valve's pocket
    const pocketRadius = Math.max(headRadius + lipClearance + 0.0002, Math.min(headRadius + 0.004, valveGap / 2 - 0.0006));
    const insertHeight = margin + faceHeight + 0.006;
    const portDir: 1 | -1 = vs.position[0] > 1e-6 ? 1 : vs.position[0] < -1e-6 ? -1 : i === 0 ? -1 : 1;
    const portDiam = i === 0 ? spec.manifolds.intakePortDiameter : spec.manifolds.exhaustPortDiameter;
    const portExitRadius = clamp(portDiam / 2, 0.3 * seatInnerRadius, 0.36 * headThickness);
    const portY = Math.max(insertHeight + seatInnerRadius * 1.1, 0.46 * headThickness);
    const springMeanRadius = stemRadius + 0.11 * B;
    const springWireRadius = 0.022 * B;
    const springInstalledLength = Math.max(0.5 * vs.headDiameter, 0.03) + vs.maxLift * 1.2;
    const springSeatY = headThickness + guideBoss;
    const tipY = springSeatY + springInstalledLength + retainerThickness + tipProtrusion;
    const zPush = zRocker0 - i * lane;
    const d = Math.hypot(vs.position[0] - vs.position[0], vs.position[1] - zPush); // plan length pushrod cup → valve axis
    const armPushrod = d / (1 + ROCKER_RATIO);
    const armValve = d - armPushrod;
    const dirX = 0, dirZ = 1; // pushrod cup sits at the valve's x, behind it
    const pivot: Vec3 = [vs.position[0] - dirX * armValve, tipY, vs.position[1] - dirZ * armValve];
    return {
      key: keys[i],
      spec: vs,
      x: vs.position[0],
      z: vs.position[1],
      headRadius,
      seatInnerRadius,
      stemRadius,
      margin,
      faceHeight,
      tulipHeight,
      pocketRadius,
      lipClearance,
      insertHeight,
      tipY,
      guideTopY: springSeatY,
      guideBottomY: portY,
      guideOuterRadius: stemRadius + 0.055 * B,
      springSeatY,
      springInstalledLength,
      springMeanRadius,
      springWireRadius,
      springTurns: 6.5,
      retainerThickness,
      portDir,
      portExitRadius,
      portY,
      rockerPivot: pivot,
      rockerDir: [dirX, dirZ],
      rockerYaw: Math.atan2(dirX, dirZ),
      armValve,
      armPushrod,
      pushrodTopRest: [vs.position[0], tipY, zPush],
      tappetX: camX,
      lobeZ: zPush,
    };
  }) as [ValveLayout, ValveLayout];

  const headZMax = headXHalf;
  const headZMin = Math.min(-headXHalf, valves[1].lobeZ - pushrodPassageRadius - 0.012, valves[0].lobeZ - pushrodPassageRadius - 0.012);

  const maxTipY = Math.max(valves[0].tipY, valves[1].tipY);
  const rockerCover = { yBottom: headThickness, yTop: maxTipY + 0.5 * B, wall: 0.07 * B };

  // ---- cam, tappets, timing gears ----
  const lobeWidth = 0.19 * B;
  const maxTappetLift = Math.max(valves[0].spec.maxLift, valves[1].spec.maxLift) / ROCKER_RATIO;
  const camBaseRadius = Math.max(0.3 * B, 0.012);
  const camShaftRadius = 0.16 * B;
  const lobeMaxRadius = camBaseRadius * 1.35 + maxTappetLift;
  const camY = Math.max(mainJournalRadius + lobeMaxRadius + 0.012, 1.25 * B);
  const tappetRadius = 0.13 * B;
  const tappetFaceRadius = 0.2 * B;
  const tappetLength = 0.4 * B;
  const tappetGuide = {
    yBottom: camY + camBaseRadius * 1.35 + maxTappetLift + 0.012,
    yTop: Math.min(yDeckBottom - 0.004, camY + camBaseRadius * 1.35 + maxTappetLift + 0.012 + 0.34 * B),
  };
  const lobeZs = [valves[0].lobeZ, valves[1].lobeZ];
  const camZFront = Math.max(...lobeZs) + lobeWidth / 2 + 0.01;
  const gearWidth = 0.22 * B;
  const gearZ = Math.min(...lobeZs) - lobeWidth / 2 - gearWidth / 2 - 0.012;
  const camDist = Math.hypot(camX, camY);
  const crankGearRadius = camDist / 3;
  const camGearRadius = (2 * camDist) / 3;
  const module = 0.0032;
  const crankGearTeeth = Math.max(12, Math.round((2 * crankGearRadius) / module));
  const camGearTeeth = 2 * crankGearTeeth;

  // ---- crankcase ----
  const wall = 0.17 * B;
  const bigEndReach = Math.hypot(rodBoltOffset + rodBoltRadius * 1.6, rodCapDepth + 0.01);
  const sweep = Math.max(r + bigEndReach, counterweightRadius);
  const xHalfInner = sweep + 0.014;
  const xHalfOuter = xHalfInner + wall;
  const yFloorInner = -(sweep + 0.014);
  const zFrontInner = webZ + webThickness / 2 + 0.024;
  const zFrontOuter = zFrontInner + wall;
  const zBackInner = gearZ - gearWidth / 2 - 0.012;
  const zBackOuter = zBackInner - wall;
  const pedestalThickness = 0.3 * B;
  // Intermediate main bearing between the rear web and the front cam lobe.
  const pedestalZ = (-(webZ + webThickness / 2) + (Math.max(...lobeZs) + lobeWidth / 2)) / 2;

  // ---- flywheel ----
  const flywheelRadius = 2.8 * B;
  const flywheelWidth = 0.9 * B;
  const flywheelZ = zBackOuter - 0.03 - flywheelWidth / 2;
  const flywheelGrooveRadius = flywheelRadius - 0.012;
  const bedplateTop = -(flywheelRadius + 0.025);
  const yBase = bedplateTop;
  const frontJournalEndZ = zFrontOuter + 0.012;
  const rearJournalEndZ = flywheelZ - flywheelWidth / 2 - 0.006;

  const motorRatio = 3;
  const motorRadius = 1.3 * B;
  // cutSide is decided below from the valves; drive + CR worm go on the opposite (uncut) side
  const farSide = valveSpecs[0].position[0] > 1e-9 ? -1 : 1;
  const motor = {
    x: farSide * (flywheelRadius + motorRadius + 0.18),
    y: bedplateTop + motorRadius + 0.02,
    radius: motorRadius,
    length: 3.4 * B,
    zCenter: flywheelZ + 0.03 + 1.7 * B,
    pulleyRadius: flywheelGrooveRadius / motorRatio,
    ratio: motorRatio,
  };

  const wormWheelPitch = wormWheelOuter - 0.012;
  const wormWheelTeeth = Math.max(24, Math.round((2 * Math.PI * wormWheelPitch) / 0.008));
  const wormRadius = 0.2 * B;
  const worm = {
    x: farSide * (wormWheelPitch + wormRadius),
    y: yDeckTop + wormWheelThickness / 2,
    radius: wormRadius,
    zFront: zFrontOuter + 0.05,
    zBack: -0.5 * B,
  };

  // ---- spark plug ----
  const sp = spec.sparkPlug;
  const axis = normalize3(sp.axis);
  const gc: Vec3 = [sp.gapCenter[0], sp.gapCenter[1], sp.gapCenter[2]];
  const cutSide: 1 | -1 = valveSpecs[0].position[0] > 1e-9 ? 1 : -1;
  let side = cross3(axis, [0, 1, 0]);
  if (Math.hypot(...side) < 0.2) side = cross3(axis, [0, 0, 1]);
  side = normalize3(side);
  // ground strut away from the default viewing direction (from the cut-away quadrant)
  if (side[0] * cutSide + side[2] > 1e-9) side = [-side[0], -side[1], -side[2]];
  const insideChamber = (x: number, y: number, z: number) => x * x + z * z < boreRadius * boreRadius && y < 0;
  const insideCasting = (x: number, y: number, z: number) => {
    const inHead = y >= 0 && y <= headThickness && Math.abs(x) <= headXHalf && z >= headZMin && z <= headZMax;
    const rr = x * x + z * z;
    const inCyl = y < 0 && y >= -jacketDepth && rr <= jacketOuterRadius * jacketOuterRadius;
    return inHead || inCyl;
  };
  const noseDistance = exitDistance(gc, axis, insideChamber, 0, 0.2);
  const castingExit = exitDistance(gc, axis, (x, y, z) => insideCasting(x, y, z) || insideChamber(x, y, z), noseDistance, 0.4);
  const threadRadius = sp.threadDiameter / 2;
  const plug = {
    gapCenter: gc,
    axis,
    side,
    noseDistance,
    seatDistance: castingExit + 0.006,
    threadRadius,
    bossRadius: threadRadius + 0.007,
  };

  // Keep the water jacket clear of the plug boss when the plug goes through the barrel.
  let jacketTop = -0.006;
  if (gc[1] < 0 && Math.abs(axis[1]) < 0.7) {
    jacketTop = Math.min(jacketTop, gc[1] - plug.bossRadius - 0.004);
  }
  const jacketCavity = {
    rIn: linerOuterRadius,
    rOut: jacketOuterRadius - 0.1 * B,
    yTop: jacketTop,
    yBottom: -jacketDepth + 0.1 * B,
  };

  const crankcase = {
    xHalfOuter,
    xHalfInner,
    yBase,
    yFloorInner,
    yDeckBottom,
    yDeckTop,
    zFrontOuter,
    zFrontInner,
    zBackOuter,
    zBackInner,
    wall,
    baseFlangeHalf: xHalfOuter + 0.04,
    pedestalZ,
    pedestalThickness,
  };

  return {
    spec,
    cutSide,
    phiTdc,
    pistonTravel,
    bore: B,
    stroke: S,
    crankRadius: r,
    rodLength: L,
    pinOffset: e,
    compressionHeight: CH,
    crRange,
    mainJournalRadius,
    crankPinRadius,
    crankPinLength,
    webThickness,
    webZ,
    webPinBossRadius,
    webJournalBossRadius,
    counterweightRadius,
    counterweightHalfAngle: (68 * Math.PI) / 180,
    frontJournalEndZ,
    rearJournalEndZ,
    rodBigEndWidth,
    rodBigEndOuterRadius,
    rodBearingThickness,
    rodBoltOffset,
    rodBoltRadius,
    rodCapDepth,
    rodSmallEndWidth,
    rodSmallEndOuterRadius,
    rodBushingThickness,
    rodShankWidthBig: 0.4 * B,
    rodShankWidthSmall: 0.27 * B,
    rodFlangeThickness: 0.075 * B,
    rodWebThickness: 0.11 * B,
    pistonRadius,
    pistonCrownThickness: 0.12 * B,
    pistonSkirtBottom,
    pistonWallThickness: 0.05 * B,
    ringGrooves,
    wristPinRadius,
    wristPinBore: 0.55 * wristPinRadius,
    wristPinLength: B - 0.004,
    pinBossRadius: wristPinRadius + 0.07 * B,
    wristPinTdcY,
    crownTdcY,
    boreRadius,
    linerOuterRadius,
    jacketOuterRadius,
    jacketDepth,
    cylinderLength,
    spigotRadius,
    jacketCavity,
    threadPitch,
    headThickness,
    headXHalf,
    headZMin,
    headZMax,
    valves,
    pushrodRadius,
    pushrodTubeRadius,
    pushrodPassageRadius,
    rockerCover,
    camX,
    camY,
    camShaftRadius,
    lobeWidth,
    camBaseRadius,
    tappetRadius,
    tappetFaceRadius,
    tappetLength,
    tappetGuide,
    camZFront,
    camZBack: zBackOuter - 0.004,
    gearZ,
    gearWidth,
    crankGearRadius,
    camGearRadius,
    crankGearTeeth,
    camGearTeeth,
    crankcase,
    flywheelRadius,
    flywheelWidth,
    flywheelZ,
    flywheelGrooveRadius,
    bedplateTop,
    motor,
    wormWheel: {
      innerRadius: wormWheelInner,
      outerRadius: wormWheelOuter,
      pitchRadius: wormWheelPitch,
      thickness: wormWheelThickness,
      yBottom: yDeckTop + 0.0005,
      teeth: wormWheelTeeth,
    },
    worm,
    plug,
    clearanceAtTdc: clr,
    headY,
    headYMin,
    headYMax,
  };
}
