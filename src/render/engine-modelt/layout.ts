/**
 * Procedural layout of the Ford Model T mechanism: every dimension and placement the 3D model uses,
 * from the EngineSpec (bore, stroke, rod, compression height, L-head chamber, valves, plug, cylinder
 * layout) plus period dimensions that the spec does not carry. Pure data, no three.js, so clearances
 * and kinematics are unit-tested (layout.test.ts).
 *
 * Frames (engine-spec.ts): ROOT/world +y up, crank axis = z through the origin, +z = FRONT (fan, timing
 * gears, timer), cylinder 1 at the front, crank clockwise seen from the front. The valves stand in one row
 * on the RIGHT side of the car = −x. Cylinder i's frame = the physics L-head cylinder frame (origin on the
 * bore axis at the roof of the head cavity over the bore) translated to (0, cylOriginY, axisZ[i]) and
 * mirrored z → −z where layout.mirrorZ[i].
 *
 * Sources (short keys as in physics/engines/model-t.ts):
 *   [DB23]   1923 Ford Dealers' Data Book (MTFCA transcription)     [Good22] A. A. Good 1922
 *   [Page29] V. W. Pagé 1929                                       [Dyke24] Dyke's 1924
 *   [FM19]   Ford Manual 1919                                      [McC]    MTFCA Encyclopedia
 *   [P1917]  Ford Price List of Parts and Accessories, Model T, Aug 15 1917 (archive.org)
 *   [T424B]  Ford drawing T-424-B valve; [T431] drawing T-431 valve spring (The Henry Ford)
 *   [FS]     Ford Service (c. 1925–26) pars. 258–264;  [Tulsa] MTFC Tulsa technical pages
 * (gathered in the Model T research notes: geometry-masses, valvetrain). Everything not traceable to a
 * source is a render-only choice marked UNVERIFIED; none of it feeds back into the physics.
 */
import type { EngineSpec, ValveSpec } from '../../physics/core/engine-spec';
import { ValveCam } from './cam';
import { chamberOutlines, roundedRectPolygon, type ChamberOutlines, type P2 } from './outline';

export type Vec3 = [number, number, number];

const IN = 0.0254;
const DEG = Math.PI / 180;

/** Polygon resolution of the bore circles (block deck, bore walls, head cavities). */
export const BORE_SEG = 64;
/** Polygon resolution of valve seats / ports (even: lathe ↔ swept-tube ring matching). */
export const PORT_SEG = 32;

export interface ModelTValveLayout {
  /** 0..7, cylinder-major (cyl 1 intake, cyl 1 exhaust, cyl 2 intake, …). */
  index: number;
  /** 0-based cylinder. */
  cyl: number;
  kind: 'intake' | 'exhaust';
  spec: ValveSpec;
  cam: ValveCam;
  /** Valve axis, world x / z. */
  x: number;
  z: number;
  /** Cylinder-frame z of the spec (before mirroring). */
  zLocal: number;
  firingOffsetDeg: number;
  /** Rotation of this valve's lobe about the cam axis inside the camshaft group, rad: −(offset + θ_c)/2. */
  lobePhase: number;
}

export interface MainBearingLayout {
  /** Centre and length along z, m. */
  z: number;
  length: number;
}

export interface PlugLayoutMT {
  /** Gap centre, cylinder frame. */
  gapCenter: Vec3;
  axis: Vec3;
  side: Vec3;
  noseDistance: number;
  seatDistance: number;
  threadRadius: number;
  bossRadius: number;
}

export interface ModelTLayout {
  spec: EngineSpec;
  nCyl: number;
  bore: number;
  boreRadius: number;
  stroke: number;
  crankRadius: number;
  rodLength: number;
  compressionHeight: number;
  axisZ: number[];
  mirror: boolean[];
  firingOffsetDeg: number[];
  /** Crank-local angle of each throw (pin at r·(sin, cos) of it), rad = −firingOffset. */
  throwAngle: number[];

  // ---- vertical stack (world y) ----
  wristPinTdcY: number;
  crownTdcY: number;
  /** Block deck (gasket face of the block) = valve-seat plane. */
  deckY: number;
  crownAboveDeck: number;
  /** Bore-column depth at TDC (roof → crown), m. */
  hTdc: number;
  /** World y of the cylinder-frame origin (roof of the head cavity over the bore). */
  cylOriginY: number;
  gasketThickness: number;
  headBottomY: number;
  headTopY: number;
  boltBossTopY: number;
  /** Pocket roof, world y. */
  pocketRoofY: number;

  crank: {
    journalRadius: number;
    pinRadius: number;
    pinLength: number;
    webThickness: number;
    webBossRadius: number;
    mains: MainBearingLayout[];
    flangeRadius: number;
    flangeZ0: number;
    flangeZ1: number;
    noseRadius: number;
    frontEndZ: number;
    pulleyRadius: number;
    pulleyZ0: number;
    pulleyZ1: number;
  };

  rod: {
    length: number;
    bigEndWidth: number;
    pinRadius: number;
    babbitt: number;
    boltOffset: number;
    boltRadius: number;
    capDepth: number;
    bossHeight: number;
    smallEndWidth: number;
    smallEndOuterRadius: number;
    wristPinRadius: number;
    shankWidthBig: number;
    shankWidthSmall: number;
    flangeThickness: number;
    webThickness: number;
  };

  piston: {
    radius: number;
    length: number;
    /** Skirt bottom, piston frame (origin on the pin axis), m (< 0). */
    skirtBottom: number;
    topChamfer: number;
    crownThickness: number;
    wallThickness: number;
    grooves: { top: number; height: number; depth: number; oil: boolean }[];
    pinRadius: number;
    pinBore: number;
    pinLength: number;
    bossRadius: number;
  };

  valve: {
    headRadius: number;
    margin: number;
    faceInnerRadius: number;
    seatLineRadius: number;
    seatAngle: number;
    stemRadius: number;
    headThickness: number;
    bossHeight: number;
    bossRadius: number;
    /** Valve-local y (origin at the seat line) of the top face. */
    topY: number;
    /** Valve-local y of the stem end. */
    stemEndY: number;
    /** Valve-local y of the spring-seat pin. */
    pinY: number;
    seatOuterRadius: number;
    throatRadius: number;
  };
  valves: ModelTValveLayout[];
  /** World y of the valve seat line with the valve closed. */
  seatLineY: number;
  spring: { bottomRestY: number; topY: number; installed: number; meanRadius: number; wireRadius: number; turns: number };
  cup: { radius: number; thickness: number };
  lash: number;
  tappet: { faceRadius: number; footThickness: number; stemRadius: number; length: number; topRestY: number; faceRestY: number };

  cam: {
    x: number;
    y: number;
    baseRadius: number;
    maxRadius: number;
    journalRadius: number;
    lobeWidth: number;
    bearings: MainBearingLayout[];
    z0: number;
    z1: number;
    flangeRadius: number;
  };
  gears: {
    module: number;
    crankTeeth: number;
    camTeeth: number;
    crankPitchR: number;
    camPitchR: number;
    z0: number;
    z1: number;
    helixAngle: number;
    generator: { teeth: number; pitchR: number; x: number; y: number };
  };

  block: {
    xR: number;
    xLu: number;
    xLb: number;
    yFlareTop: number;
    yFlareBot: number;
    zF: number;
    zB: number;
    yTop: number;
    yBot: number;
    /** Crankcase cavity (open at the pan rail), stepped in z. */
    cc: { x0: number; x1: number; zLo0: number; zLo1: number; zUp0: number; zUp1: number; yStep: number; yTop: number };
    /** Valve (spring/tappet) chamber, open on the right face (valve door). */
    vc: { x1: number; y0: number; y1: number; z0: number; z1: number };
    /** Manifold-port centre height at the block face, port radius there. */
    portY: number;
    portRadius: number;
    /** Port roof (gallery top) and floor below the deck. */
    portRoofY: number;
    portFloorY: number;
    /** Siamesed intake galleries (one per intake pair). */
    galleries: { x0: number; x1: number; z0: number; z1: number; zOutlet: number; valves: number[] }[];
    /** Manifold ports: z and kind (4 exhaust + 2 siamesed intake). */
    ports: { z: number; kind: 'intake' | 'exhaust'; valve: number }[];
    jacket: { x0: number; x1: number; y0: number; y1: number; z0: number; z1: number };
  };

  head: {
    xR: number;
    xL: number;
    zF: number;
    zB: number;
    jacket: { x0: number; x1: number; y0: number; y1: number; z0: number; z1: number };
    plugPassageRadius: number;
    plugTubeRadius: number;
    bolts: P2[];
    boltBossRadius: number;
    boltHeadAcrossFlats: number;
    boltHeadHeight: number;
    gooseneck: { x: number; z: number; flangeRadius: number; boreRadius: number; rise: number; reach: number };
  };
  /** Chamber plan outlines per cylinder, world (x, z). */
  chambers: ChamberOutlines[];
  /** Pocket plan per cylinder, world (x, z) polygon. */
  pockets: P2[][];
  plug: PlugLayoutMT;

  flywheel: {
    hubRadius: number;
    hubZ0: number;
    hubZ1: number;
    discZ0: number;
    discZ1: number;
    radius: number;
    ringGear: { teeth: number; pitchR: number; module: number; z0: number; z1: number };
    magnets: { count: number; rApex: number; rPole: number; legWidth: number; z0: number; z1: number };
    tripleGears: { count: number; r: number; z0: number; z1: number; teeth: number[]; module: number };
  };
  coilRing: { rIn: number; rOut: number; plateZ0: number; plateZ1: number; coilR: number; coilRadius: number; coilZ0: number; coilZ1: number; gap: number; count: number };
  bell: { inner: { r: number; z: number }[]; wall: number; frontInnerRadius: number };

  fan: { x: number; y: number; z: number; radius: number; pulleyRadius: number; ratio: number; blades: number; hubRadius: number };
  belt: { z: number; width: number; thickness: number };
  timer: { x: number; y: number; z0: number; z1: number; radius: number; rotorRadius: number; contactArcDeg: number; advanceRangeDeg: [number, number] };
  intake: { x: number; teeY: number; flangeY: number; radius: number; wall: number };
  exhaust: { x: number; y: number; z0: number; z1: number; radius: number; wall: number; outletY: number };
  carb: { x: number; z: number; topY: number; throatRadius: number; bodyRadius: number; bodyHeight: number; bowlRadius: number; bowlOffset: Vec3; plateClosedAngle: number };
  generator: { x: number; y: number; radius: number; z0: number; z1: number };
  starter: { x: number; y: number; radius: number; z0: number; z1: number; pinionTeeth: number; pinionR: number; pinionZ: number };

  /**
   * Valve side (−1 = −x): the cutaway removes {cutSide·x > 0, z > sectionZ[k]} of block and head, the
   * quarter section through the front valve axis of the section cylinder k (the focus cylinder; cutZ =
   * sectionZ[0], cylinder 1's front valve).
   */
  cutSide: 1 | -1;
  cutZ: number;
  /** Per cylinder: z of its front (largest-z) valve axis, the section plane when it is the section cylinder. */
  sectionZ: number[];
  /** The bell (transmission cover) set removes its whole valve-side half: frame behind the bell. */
  bellCutZ: number;
  bounds: { min: Vec3; max: Vec3 };
}

function norm3(v: Vec3): Vec3 {
  const l = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / l, v[1] / l, v[2] / l];
}

/**
 * Layout of a Model T-type engine: inline, side valves ('l-head' chamber), 3 mains, cam low on the
 * valve side, timing gears + timer at the front, flywheel magneto inside a bell at the rear.
 * Throws if the spec has no L-head chamber or no cylinder layout.
 */
export function computeModelTLayout(spec: EngineSpec): ModelTLayout {
  const g = spec.geometry;
  const lh = g.lHead;
  const lay = spec.layout;
  if (g.chamber !== 'l-head' || !lh) throw new Error(`computeModelTLayout: '${spec.name}' has no L-head chamber`);
  if (!lay || lay.axisZ.length !== spec.cylinders) throw new Error(`computeModelTLayout: '${spec.name}' has no cylinder layout`);
  const n = spec.cylinders;
  const B = g.bore, S = g.stroke, r = S / 2, Lr = g.conRodLength, CH = g.compressionHeight;
  const rb = B / 2;
  const axisZ = lay.axisZ.slice();
  const mirror = axisZ.map((_, i) => !!lay.mirrorZ[i]);
  const offs = axisZ.map((_, i) => lay.firingOffsetDeg[i] ?? 0);

  // ---- vertical stack (no pin offset on the Model T: spec pinOffset 0) ----
  const wristPinTdcY = Math.sqrt(Math.max((Lr + r) ** 2 - g.pinOffset ** 2, 0));
  const crownTdcY = wristPinTdcY + CH;
  const crownAboveDeck = lh.crownAboveDeckAtTDC;
  const hTdc = -(lh.deckY + crownAboveDeck);
  const cylOriginY = crownTdcY + hTdc;
  const deckY = cylOriginY + lh.deckY;
  const pocketRoofY = cylOriginY + lh.pocket.roofY;
  // Copper-asbestos gasket ≈ 0.044 in compressed (Tulsa 0.8 in³ gasket volume / 18.3 in² opening, the
  // spec's chamber build-up) — UNVERIFIED (no primary thickness).
  const gasketThickness = 0.044 * IN;
  const headBottomY = deckY + gasketThickness;
  // High head: 2-11/16 in from the gasket face to the underside of the bolt heads (bolt bosses) [Lang's
  // 3001HP; McC T-483C1]. Body top 1/4 in below the boss tops — UNVERIFIED.
  const boltBossTopY = headBottomY + (2 + 11 / 16) * IN;
  const headTopY = boltBossTopY - 0.25 * IN;

  // ---- crankshaft ----
  // All journals 1.248 in [DB23; Page29 p. 403]; rod journals 1.495–1.505 in long [DB23]; mains 2, 2-3/16,
  // 3-1/8 in long [DB23]; overall 25-5/32 in [DB23; Page29]; no counterweights [Page29; MCCTA rules].
  const journalRadius = 1.248 * IN / 2;
  const pinRadius = journalRadius;
  const pinLength = 1.5 * IN;
  const webThickness = (5 / 8) * IN; // UNVERIFIED
  const webBossRadius = journalRadius + 0.375 * IN; // UNVERIFIED (slim webs, no counterweights)
  const mainLen = [2 * IN, (2 + 3 / 16) * IN, (3 + 1 / 8) * IN];
  const zWebFront = axisZ[0] + pinLength / 2 + webThickness; // front face of cylinder 1's front web
  const zWebRear = axisZ[n - 1] - pinLength / 2 - webThickness;
  const mbz = lay.mainBearingZ;
  const mains: MainBearingLayout[] = mbz.map((z, k) => ({ z, length: mainLen[Math.min(k, mainLen.length - 1)] }));
  // The spec's rear main (UNVERIFIED, ≈ 66 mm outboard) would overlap cylinder 4's web with the 3-1/8 in
  // rear bearing: keep the bearing length and slide it rearward until it clears the web by 1.5 mm.
  {
    const f = mains[0], b = mains[mains.length - 1];
    if (f.z - f.length / 2 < zWebFront + 0.0015) f.z = zWebFront + 0.0015 + f.length / 2;
    if (b.z + b.length / 2 > zWebRear - 0.0015) b.z = zWebRear - 0.0015 - b.length / 2;
  }
  const rearMain = mains[mains.length - 1];
  const flangeRadius = 2 * IN; // UNVERIFIED (flywheel on 4 cap screws + 2 dowels [Good22])
  const flangeZ0 = rearMain.z - rearMain.length / 2 - 0.3 * IN;
  const flangeZ1 = rearMain.z - rearMain.length / 2;
  const frontEndZ = flangeZ0 + (25 + 5 / 32) * IN;

  // ---- block outline [DB23]: 22-9/16 in long (outside of the cylinder chambers), 10-5/8 in high, 9-1/2 in
  // wide at the base. Deck = 10-5/8 in above the crank axis by the spec's compression-height derivation. ----
  const blockLen = (22 + 9 / 16) * IN;
  const zF = blockLen / 2, zB = -blockLen / 2;
  const xLb = (9.5 * IN) / 2;
  const xR = -xLb; // valve-side face (manifold face + valve door)
  const xLu = rb + 0.95 * IN; // UNVERIFIED: left wall over the cylinder bank (bore + jacket)
  const yFlareTop = 0.105, yFlareBot = 0.065; // UNVERIFIED

  // ---- valves ----
  // T-424-B [T424B; Good22 p. 21]: head 1.469–1.484 in, 3/16 in thick, head + boss 9/16 in, margin
  // .023–.031 in, 45° face from 1-17/64 in up to the head edge, seat (gauge) line 1.422 in, stem .311 in,
  // seat line → stem end 4.974 in, spring-pin hole 4-19/32 in below the seat line.
  // Block seat [Good22 pp. 20–21]: 1-5/16 in port throat, countersink 1-23/64 … 1-29/64 in.
  const vi = spec.intakeValve;
  const seatAngle = vi.seatAngle;
  const tanA = Math.tan(seatAngle);
  const valve = {
    headRadius: vi.headDiameter / 2,
    margin: 0.027 * IN,
    faceInnerRadius: (1 + 17 / 64) * IN / 2,
    seatLineRadius: 1.422 * IN / 2,
    seatAngle,
    stemRadius: vi.stemDiameter / 2,
    headThickness: (3 / 16) * IN,
    bossHeight: (9 / 16) * IN,
    bossRadius: (5 / 8) * IN / 2,
    topY: 0,
    stemEndY: -4.974 * IN,
    pinY: -(4 + 19 / 32) * IN,
    seatOuterRadius: (1 + 29 / 64) * IN / 2,
    throatRadius: vi.seatInnerDiameter / 2,
  };
  valve.topY = (valve.headRadius - valve.seatLineRadius) * tanA + valve.margin;
  // closed: the seat line sits on the 45° countersink below the deck
  const seatLineY = deckY - (valve.seatOuterRadius - valve.seatLineRadius) * tanA;
  // spring T-431 (1923–26): 12 coils of .105 in wire, 31/32 in OD [Dyke24], installed at 2-1/8 in [FS par. 258];
  // spring seat (cup) 1-1/8 in dia [P1917], pinned through the stem (7/64 in pin [P1917]); cup 1/16 in thick
  // UNVERIFIED. The spring's upper end bears on the block (valve-chamber roof) [FS pars. 208, 259].
  const wire = 0.105 * IN;
  const springOD = (31 / 32) * IN;
  const cup = { radius: (1 + 1 / 8) * IN / 2, thickness: (1 / 16) * IN };
  const pinYw = seatLineY + valve.pinY;
  const springBottomRest = pinYw + (7 / 64) * IN / 2 + cup.thickness;
  const springInstalled = (2 + 1 / 8) * IN;
  const springTop = springBottomRest + springInstalled;
  const spring = { bottomRestY: springBottomRest, topY: springTop, installed: springInstalled, meanRadius: (springOD - wire) / 2, wireRadius: wire / 2, turns: 12 };
  // Push rod (mushroom tappet) [Good22 pp. 21, 26; P1917]: 1 in flat foot, .436 in stem, 2-11/32 in long;
  // foot 5/32 in thick UNVERIFIED. Lash between push rod and stem end at the cam heel [FM19 A22].
  const lash = Math.max(0, vi.lash ?? 0);
  const stemEndRest = seatLineY + valve.stemEndY;
  const tappetTopRest = stemEndRest - lash;
  const tappetLength = (2 + 11 / 32) * IN;
  const tappetFaceRest = tappetTopRest - tappetLength;
  const tappet = { faceRadius: 0.5 * IN, footThickness: (5 / 32) * IN, stemRadius: 0.436 * IN / 2, length: tappetLength, topRestY: tappetTopRest, faceRestY: tappetFaceRest };

  // ---- camshaft: low on the valve side, tappets straight under the stems (ratio 1) ----
  const camI = new ValveCam(spec.intakeValve, 0.406 * IN);
  const camE = new ValveCam(spec.exhaustValve, 0.406 * IN);
  const camBase = camI.baseRadius;
  const camRise = camI.geometry ? camI.geometry.rise : vi.maxLift;
  const camX = vi.position[0];
  const camY = tappetFaceRest - camBase;
  const valves: ModelTValveLayout[] = [];
  for (let c = 0; c < n; c++) {
    for (const [kind, vs, vc] of [['intake', spec.intakeValve, camI], ['exhaust', spec.exhaustValve, camE]] as const) {
      const zl = vs.position[1];
      valves.push({
        index: valves.length,
        cyl: c,
        kind,
        spec: vs,
        cam: vc,
        x: vs.position[0],
        z: axisZ[c] + (mirror[c] ? -zl : zl),
        zLocal: zl,
        firingOffsetDeg: offs[c],
        lobePhase: -((offs[c] + vc.centreDeg) * DEG) / 2,
      });
    }
  }
  // Camshaft [Good22 pp. 21, 26; McC]: lobes 7/8 in wide, journals .748 in, bearings front 1.967, centre
  // 2-7/16, rear 1.750 in long, 22-23/32 in long, flange 1-3/4 in; three bearings (front, centre, rear).
  const lobeWidth = (7 / 8) * IN;
  const camJournalRadius = 0.748 * IN / 2;
  const vz = valves.map((v) => v.z);
  const lobeFront = Math.max(...vz) + lobeWidth / 2, lobeRear = Math.min(...vz) - lobeWidth / 2;
  const camBearings: MainBearingLayout[] = [
    { z: lobeFront + 0.003 + (1.967 * IN) / 2, length: 1.967 * IN },
    { z: 0, length: (2 + 7 / 16) * IN },
    { z: lobeRear - 0.003 - (1.75 * IN) / 2, length: 1.75 * IN },
  ];
  // Timing gears: 48 / 24 spiral-cut (c. 1919/20–27) [Good22 pp. 24–25]; width 1 in UNVERIFIED; helix 22°
  // UNVERIFIED. The module follows from the stack-up centre distance (≈ 3.9–4.0 in, valvetrain research).
  const camDist = Math.hypot(camX, camY);
  const crankTeeth = 24, camTeeth = 48;
  const module = (2 * camDist) / (crankTeeth + camTeeth);
  const gearZ0 = zF + 0.003, gearZ1 = gearZ0 + 1.0 * IN;
  // generator (1919+ starter cars) is driven by a pinion off the large time gear [McC 1917-20; FSB];
  // pinion teeth and angle UNVERIFIED.
  const genTeeth = 24;
  const genPitchR = (genTeeth * module) / 2;
  const genAng = 191 * DEG;
  const genDist = camTeeth * module / 2 + genPitchR;
  const generatorAxis = { x: camX + genDist * Math.cos(genAng), y: camY + genDist * Math.sin(genAng) };

  // ---- piston [DB23; Dyke24 Instr. 70]: cast iron, 3.81 in long, three 1/4 × 13/64 in rings, two above the
  // pin and the oil ring below it; top edge bevelled. Land heights UNVERIFIED (top land 1/4 in, 3/16 in land).
  const ringH = 0.25 * IN, grooveDepth = (13 / 64) * IN;
  const pistonPinRadius = 0.74 * IN / 2;
  const bossRadius = pistonPinRadius + 0.3 * IN; // UNVERIFIED
  const grooves = [
    { top: 0.25 * IN, height: ringH, depth: grooveDepth, oil: false },
    { top: 0.25 * IN + ringH + (3 / 16) * IN, height: ringH, depth: grooveDepth, oil: false },
    { top: CH + bossRadius + 0.06 * IN, height: ringH, depth: grooveDepth, oil: true },
  ];
  const pistonLength = 3.81 * IN;
  const piston = {
    radius: rb - 0.0035 * IN / 2, // skirt clearance .003–.004 in [DB23]
    length: pistonLength,
    skirtBottom: CH - pistonLength,
    topChamfer: (3 / 32) * IN, // UNVERIFIED bevel
    crownThickness: (5 / 16) * IN, // UNVERIFIED
    wallThickness: (3 / 16) * IN, // UNVERIFIED
    grooves,
    pinRadius: pistonPinRadius,
    pinBore: 0.45 * IN / 2, // UNVERIFIED (hollow seamless tube)
    pinLength: 3.5 * IN, // [DB23; Dyke24]
    bossRadius,
  };

  // ---- connecting rod [DB23; Good22 p. 27]: I-beam, babbitt big end 1.495 in wide (1/16 in babbitt from
  // 1915 [McC]), two-bolt cap, small end clamped on the pin (63/64–1-1/64 in wide). Section sizes UNVERIFIED.
  const babbitt = (1 / 16) * IN;
  const rodBoltRadius = (5 / 16) * IN / 2; // UNVERIFIED
  const rod = {
    length: Lr,
    bigEndWidth: 1.495 * IN,
    pinRadius,
    babbitt,
    boltOffset: pinRadius + babbitt + rodBoltRadius + 0.1 * IN,
    boltRadius: rodBoltRadius,
    capDepth: pinRadius + babbitt + 0.3 * IN,
    bossHeight: 0.3 * IN,
    smallEndWidth: 1.0 * IN,
    smallEndOuterRadius: pistonPinRadius + 0.24 * IN,
    wristPinRadius: pistonPinRadius,
    shankWidthBig: 0.85 * IN,
    shankWidthSmall: 0.62 * IN,
    flangeThickness: 0.16 * IN,
    webThickness: 0.16 * IN,
  };

  // ---- block internals ----
  const boreLen = 6.752 * IN; // [DB23]
  const ccTop = deckY - boreLen;
  const lobeR = camBase + camRise;
  const ccX0 = Math.min(camX - lobeR, camX - tappet.faceRadius) - 0.012;
  const ccX1 = r + rod.boltOffset + rod.boltRadius * 2 + 0.012;
  const zLo1 = Math.max(zWebFront, lobeFront) + 0.003;
  const zLo0 = Math.min(zWebRear, lobeRear) - 0.003;
  // piston skirts at BDC reach below the bore: the upper part of the cavity extends past the bores
  const skirtBdcY = wristPinTdcY - S + piston.skirtBottom;
  const zUp1 = axisZ[0] + rb + 0.007, zUp0 = axisZ[n - 1] - rb - 0.007;
  const yStep = Math.min(skirtBdcY - 0.006, ccTop - 0.006);
  // valve chamber: from above the tappet guides to the spring tops, open to the right face
  const vc = {
    x1: -(rb + 0.27 * IN), // UNVERIFIED inner wall (≈ 1/4 in of iron beside the bore)
    y0: ccTop + 0.45 * IN, // UNVERIFIED tappet-guide wall
    y1: springTop,
    z0: Math.min(...vz) - cup.radius - 0.012,
    z1: Math.max(...vz) + cup.radius + 0.012,
  };
  // Ports: six 1-1/8 in manifold ports (4 exhaust + 2 siamesed intake) [Good22 p. 21]; port centre height and
  // the cored shape UNVERIFIED (roof 1/2 in under the deck, floor ≈ 1.65 in).
  const portRadius = spec.manifolds.exhaustPortDiameter / 2;
  const portRoofY = deckY - 0.45 * IN;
  const portFloorY = deckY - 1.7 * IN;
  const portY = (portRoofY + portFloorY) / 2;
  const ports: { z: number; kind: 'intake' | 'exhaust'; valve: number }[] = [];
  const galleries: ModelTLayout['block']['galleries'] = [];
  for (const v of valves) if (v.kind === 'exhaust') ports.push({ z: v.z, kind: 'exhaust', valve: v.index });
  // siamese pairs: intake valves closer than 1.5 port diameters' worth of spacing share one port
  const intakes = valves.filter((v) => v.kind === 'intake').sort((a, b) => b.z - a.z);
  for (let k = 0; k < intakes.length; k++) {
    const a = intakes[k], b = intakes[k + 1];
    if (b && Math.abs(a.z - b.z) < 2.2 * B && !galleries.some((gl) => gl.valves.includes(a.index))) {
      const zo = (a.z + b.z) / 2;
      galleries.push({
        x0: camX - valve.throatRadius - 0.0015,
        x1: camX + valve.throatRadius + 0.0015,
        z0: Math.min(a.z, b.z) - valve.throatRadius - 0.004,
        z1: Math.max(a.z, b.z) + valve.throatRadius + 0.004,
        zOutlet: zo,
        valves: [a.index, b.index],
      });
      ports.push({ z: zo, kind: 'intake', valve: a.index });
      k++;
    } else if (!galleries.some((gl) => gl.valves.includes(a.index))) {
      ports.push({ z: a.z, kind: 'intake', valve: a.index });
    }
  }
  ports.sort((p, q) => q.z - p.z);
  const jacket = { x0: -0.04, x1: xLu - 0.008, y0: ccTop + 0.02, y1: deckY - 0.008, z0: axisZ[0] + rb + 0.005, z1: zF - 0.008 };

  // ---- head ----
  // 15 cap screws, 7/16-14 [FM19; Page29]. Pattern UNVERIFIED: rows beside the bores (+x), beyond the valve
  // pockets (−x) at every cylinder centre and both ends, plus three on the centreline (front, centre, rear).
  const pk = lh.pocket;
  const hxR = xR + 0.001, hxL = xLu - 0.001, hzF = zF - 0.001, hzB = zB + 0.001;
  const boltL = rb + 0.5 * IN;
  const boltR = pk.xMin - 0.5 * IN;
  const zEndF = zF - 0.55 * IN, zEndB = zB + 0.55 * IN;
  const bolts: P2[] = [];
  for (const x of [boltL, boltR]) {
    bolts.push([x, zEndF]);
    for (const z of axisZ) bolts.push([x, z]);
    bolts.push([x, zEndB]);
  }
  bolts.push([0, zEndF], [0, (axisZ[1] + axisZ[2]) / 2], [0, zEndB]);
  const plugThreadR = spec.sparkPlug.threadDiameter / 2;
  const plugPassageRadius = plugThreadR + 0.0003;
  const plugTubeRadius = plugPassageRadius + 0.16 * IN;
  const headJacket = { x0: hxR + 0.012, x1: hxL - 0.012, y0: cylOriginY + 0.006, y1: headTopY - 0.007, z0: hzB + 0.014, z1: hzF - 0.014 };

  // ---- chamber outlines (world plan) ----
  const pockets: P2[][] = [];
  const chambers: ChamberOutlines[] = [];
  for (let c = 0; c < n; c++) {
    const zc = axisZ[c];
    const z0 = mirror[c] ? zc - pk.zMax : zc + pk.zMin;
    const z1 = mirror[c] ? zc - pk.zMin : zc + pk.zMax;
    const poly = roundedRectPolygon(pk.xMin, pk.xMax, z0, z1, pk.cornerRadius, 10);
    pockets.push(poly);
    chambers.push(chamberOutlines(0, zc, rb, BORE_SEG, poly));
  }

  // ---- spark plug (cylinder frame): 1/2 in pipe thread long-body plug [Dyke 1917 Ford Suppl.] screwed
  // vertically through the head over the valve pocket; shell nose flush with the pocket roof. ----
  const sp = spec.sparkPlug;
  const axis = norm3(sp.axis);
  let side: Vec3 = [axis[1], -axis[0], 0];
  if (Math.hypot(...side) < 0.2) side = [0, -axis[2], axis[1]];
  side = norm3(side);
  side = Math.abs(axis[1]) > 0.9 ? [1, 0, 0] : side; // ground strut toward +x (away from the valve-side cut)
  const gapWorldY = cylOriginY + sp.gapCenter[1];
  const noseDistance = Math.max(0.0015, (pocketRoofY - gapWorldY) / Math.max(-axis[1], 1e-6));
  const plugBossH = 0.25 * IN; // UNVERIFIED boss on the head top
  const seatDistance = (headTopY + plugBossH - gapWorldY) / Math.max(-axis[1], 1e-6);
  const plug: PlugLayoutMT = {
    gapCenter: [sp.gapCenter[0], sp.gapCenter[1], sp.gapCenter[2]],
    axis,
    side,
    noseDistance,
    seatDistance,
    threadRadius: plugThreadR,
    bossRadius: plugThreadR + 0.3 * IN,
  };

  // ---- flywheel magneto [Dyke24 pp. 248–249; FM19 A54/68; McC]: 16 V magnets (3/4 in from Sep 1914) on the
  // flywheel's front face, 16 stationary coils on a ring 1/32 in in front of them; ring gear 120 teeth
  // (Bendix 10 T) [Dyke24]; 8 DP would give the 15-1/4 in OD of the forum snippet — UNVERIFIED. ----
  const ringTeeth = 120;
  const ringModule = IN / 8;
  const ringPitchR = (ringTeeth * ringModule) / 2;
  const coilGap = (1 / 32) * IN;
  const coilPlateZ1 = zB - 0.3 * IN; // UNVERIFIED: coil ring bolted in front of the flywheel
  const coilPlateZ0 = coilPlateZ1 - 0.25 * IN;
  const coilZ1 = coilPlateZ0, coilZ0 = coilZ1 - 0.75 * IN;
  const magZ1 = coilZ0 - coilGap, magZ0 = magZ1 - 0.75 * IN;
  const discZ1 = magZ0, discZ0 = discZ1 - 1.0 * IN; // UNVERIFIED disc thickness
  const flywheel = {
    hubRadius: flangeRadius + 0.004,
    hubZ0: flangeZ0 - 0.012,
    hubZ1: flangeZ0,
    discZ0,
    discZ1,
    radius: ringPitchR + ringModule,
    ringGear: { teeth: ringTeeth, pitchR: ringPitchR, module: ringModule, z0: discZ0, z1: discZ0 + 0.75 * IN },
    magnets: { count: 16, rApex: 4.4 * IN, rPole: 7.0 * IN, legWidth: 0.75 * IN, z0: magZ0, z1: magZ1 },
    // triple gears (27/33/24 T, planetary transmission) on three pins of the flywheel's rear face
    // [model-t.ts vehicle notes]; pin radius and module UNVERIFIED.
    tripleGears: { count: 3, r: 3.3 * IN, z0: discZ0 - 2.6 * IN, z1: discZ0 - 0.05 * IN, teeth: [27, 33, 24], module: 0.09 * IN },
  };
  const coilRing = {
    rIn: 5.2 * IN,
    rOut: 7.45 * IN,
    plateZ0: coilPlateZ0,
    plateZ1: coilPlateZ1,
    coilR: 6.3 * IN,
    coilRadius: 0.55 * IN,
    coilZ0,
    coilZ1,
    gap: coilGap,
    count: 16,
  };
  // bell-shaped transmission cover + flywheel housing (cast-iron cover over the pressed-steel pan) —
  // proportions UNVERIFIED.
  const bR = flywheel.radius + 0.4 * IN;
  const bell = {
    inner: [
      { r: bR, z: zB },
      { r: bR, z: discZ0 - 0.4 * IN },
      { r: 5.9 * IN, z: discZ0 - 3.0 * IN },
      { r: 5.9 * IN, z: discZ0 - 9.5 * IN },
      { r: 2.2 * IN, z: discZ0 - 11.0 * IN },
    ],
    wall: 0.22 * IN,
    frontInnerRadius: flangeRadius + 0.01,
  };

  // ---- front: crank pulley, fan, belt, timer ----
  // Flat fan belt from the crank pulley to the two-blade fan [Page18]; pulley sizes, fan diameter and
  // position UNVERIFIED (render only).
  const coverZ1 = gearZ1 + 0.006;
  const pulleyZ0 = coverZ1 + 0.002, pulleyZ1 = Math.min(frontEndZ, pulleyZ0 + 0.875 * IN);
  const pulleyRadius = 2.0 * IN;
  const fanPulleyR = 1.5 * IN;
  const crank = {
    journalRadius,
    pinRadius,
    pinLength,
    webThickness,
    webBossRadius,
    mains,
    flangeRadius,
    flangeZ0,
    flangeZ1,
    noseRadius: journalRadius * 0.92,
    frontEndZ,
    pulleyRadius,
    pulleyZ0,
    pulleyZ1,
  };
  const belt = { z: (pulleyZ0 + pulleyZ1) / 2, width: 0.75 * IN, thickness: 0.2 * IN };
  const fan = { x: 0, y: deckY + 0.6 * IN, z: pulleyZ1 + 1.4 * IN, radius: 6.0 * IN, pulleyRadius: fanPulleyR, ratio: pulleyRadius / fanPulleyR, blades: 2, hubRadius: 1.0 * IN };
  // Timer (commutator) on the front end of the camshaft [valvetrain research; FM19], case rotated by the spark
  // lever; four segments, each grounded for `contactArcDeg` crank degrees [BP]; diameter UNVERIFIED.
  const ig = spec.ignition;
  const timer = {
    x: camX,
    y: camY,
    z0: pulleyZ1 + 0.1 * IN,
    z1: pulleyZ1 + 1.45 * IN,
    radius: 1.4 * IN,
    rotorRadius: 1.05 * IN,
    contactArcDeg: ig.type === 'trembler-magneto' ? ig.timer.contactArcDeg : 87,
    advanceRangeDeg: (ig.type === 'trembler-magneto' ? ig.timer.advanceRangeDeg.slice() : [0, 60]) as [number, number],
  };
  // The camshaft is 22-23/32 in long [Good22]; the render shaft runs from the rear bearing to the timer on
  // its nose (axial stack-up UNVERIFIED: ≈ 5 % longer than the period figure).
  const camZ1 = timer.z0 + 0.012;
  const camZ0 = camBearings[2].z - camBearings[2].length / 2 - 0.001;

  // ---- manifolds + carburettor [valvetrain research; Page 1915 Fig. 8]: cast-iron Y intake from an updraft
  // carburettor at the bottom centre (between cylinders 2 and 3), exhaust log above it with the outlet at the
  // rear turned down. Centreline offsets UNVERIFIED. ----
  const manR = 1.125 * IN / 2;
  const intake = { x: xR - 1.8 * IN, teeY: portY - 2.2 * IN, flangeY: portY - 4.6 * IN, radius: manR, wall: 0.17 * IN };
  const exhaust = { x: xR - 2.0 * IN, y: portY + 1.8 * IN, z0: Math.min(...ports.map((p) => p.z)) - 0.9 * IN, z1: Math.max(...ports.map((p) => p.z)) + 0.5 * IN, radius: 0.8 * IN, wall: 0.2 * IN, outletY: portY - 1.2 * IN };
  const closedAngle = ((spec.manifolds.throttle?.closedAngleDeg ?? 8) * Math.PI) / 180;
  const carb = {
    x: intake.x,
    z: 0,
    topY: intake.flangeY,
    throatRadius: spec.manifolds.throttleDiameter / 2,
    bodyRadius: spec.manifolds.throttleDiameter / 2 + 0.22 * IN,
    bodyHeight: 2.6 * IN,
    bowlRadius: 1.25 * IN,
    bowlOffset: [0, -1.9 * IN, 1.9 * IN] as Vec3,
    plateClosedAngle: closedAngle,
  };
  const generator = { x: generatorAxis.x, y: generatorAxis.y, radius: 1.75 * IN, z0: zF - 6.5 * IN, z1: zF };
  // Starter (1919+) on the left of the transmission cover, Bendix 10-tooth pinion on the 120-tooth ring gear
  // [Dyke24]; body size and position UNVERIFIED.
  const pinionR = (10 * ringModule) / 2;
  const stAng = 32 * DEG;
  const stDist = ringPitchR + pinionR;
  const starter = { x: stDist * Math.cos(stAng), y: stDist * Math.sin(stAng), radius: 2.2 * IN, z0: zB - 0.1 * IN, z1: zB + 8.5 * IN, pinionTeeth: 10, pinionR, pinionZ: flywheel.ringGear.z0 + 0.375 * IN };

  // ---- cutaway: quarter section through the section cylinder's front valve axis on the valve side ----
  const sectionZ = axisZ.map((_, c) => Math.max(...valves.filter((v) => v.cyl === c).map((v) => v.z)));
  const cutSide: 1 | -1 = camX < 0 ? -1 : 1;
  const bellEnd = bell.inner[bell.inner.length - 1].z - bell.wall;

  const minX = Math.min(exhaust.x - exhaust.radius - exhaust.wall, generator.x - generator.radius, carb.x - carb.bowlRadius - 0.01, -bR - bell.wall);
  const maxX = Math.max(starter.x + starter.radius, bR + bell.wall, xLb);
  const minY = Math.min(-bR - bell.wall, -0.11);
  const maxY = Math.max(fan.y + fan.radius, headTopY + 0.06);
  const bounds = { min: [minX, minY, bellEnd] as Vec3, max: [maxX, maxY, Math.max(fan.z + 0.02, timer.z1)] as Vec3 };

  return {
    spec,
    nCyl: n,
    bore: B,
    boreRadius: rb,
    stroke: S,
    crankRadius: r,
    rodLength: Lr,
    compressionHeight: CH,
    axisZ,
    mirror,
    firingOffsetDeg: offs,
    throwAngle: offs.map((o) => -o * DEG),
    wristPinTdcY,
    crownTdcY,
    deckY,
    crownAboveDeck,
    hTdc,
    cylOriginY,
    gasketThickness,
    headBottomY,
    headTopY,
    boltBossTopY,
    pocketRoofY,
    crank,
    rod,
    piston,
    valve,
    valves,
    seatLineY,
    spring,
    cup,
    lash,
    tappet,
    cam: {
      x: camX,
      y: camY,
      baseRadius: camBase,
      maxRadius: camBase + camRise,
      journalRadius: camJournalRadius,
      lobeWidth,
      bearings: camBearings,
      z0: camZ0,
      z1: camZ1,
      flangeRadius: (1.75 * IN) / 2,
    },
    gears: {
      module,
      crankTeeth,
      camTeeth,
      crankPitchR: (crankTeeth * module) / 2,
      camPitchR: (camTeeth * module) / 2,
      z0: gearZ0,
      z1: gearZ1,
      helixAngle: 22 * DEG,
      generator: { teeth: genTeeth, pitchR: genPitchR, x: generatorAxis.x, y: generatorAxis.y },
    },
    block: {
      xR,
      xLu,
      xLb,
      yFlareTop,
      yFlareBot,
      zF,
      zB,
      yTop: deckY,
      yBot: 0,
      cc: { x0: ccX0, x1: ccX1, zLo0, zLo1, zUp0, zUp1, yStep, yTop: ccTop },
      vc,
      portY,
      portRadius,
      portRoofY,
      portFloorY,
      galleries,
      ports,
      jacket,
    },
    head: {
      xR: hxR,
      xL: hxL,
      zF: hzF,
      zB: hzB,
      jacket: headJacket,
      plugPassageRadius,
      plugTubeRadius,
      bolts,
      boltBossRadius: 0.4 * IN,
      boltHeadAcrossFlats: (5 / 8) * IN,
      boltHeadHeight: 0.3 * IN,
      gooseneck: { x: 0, z: zF - 1.6 * IN, flangeRadius: 1.25 * IN, boreRadius: 0.75 * IN, rise: 1.3 * IN, reach: 2.7 * IN },
    },
    chambers,
    pockets,
    plug,
    flywheel,
    coilRing,
    bell,
    fan,
    belt,
    timer,
    intake,
    exhaust,
    carb,
    generator,
    starter,
    cutSide,
    cutZ: sectionZ[0],
    sectionZ,
    bellCutZ: bellEnd - 0.01,
    bounds,
  };
}

/** Exact slider-crank for cylinder i: piston displacement below TDC at the local angle θ_i (deg), m. */
export function pistonDisplacementAt(L: ModelTLayout, thetaLocalDeg: number): number {
  const phi = thetaLocalDeg * DEG;
  const r = L.crankRadius, l = L.rodLength;
  const px = r * Math.sin(phi), py = r * Math.cos(phi);
  return L.wristPinTdcY - (py + Math.sqrt(Math.max(l * l - px * px, 0)));
}
