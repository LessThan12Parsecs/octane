/**
 * GLSL for the in-cylinder gas volume: a single ray-marched pass over the
 * chamber. The ray is clipped ANALYTICALLY to the chamber — the disc
 * (x²+z² ≤ R², −h ≤ y ≤ 0), or for an L-head the union of the bore column and
 * the valve pocket (chamber.ts rayChamberIntervals: ≤ 2 intervals) — and to the
 * open valve heads, then split at the flame-brush shells (sphere radii r + a
 * and r − δ − a) so the thin front is always sampled densely while the uniform
 * regions need 1 sample each (or a few when the knock overlay is on).
 *
 * L-head proxies (volume.ts): bore cylinder (back faces), pocket prism (back
 * faces) and the transfer arc where the pocket meets the bore circle (front
 * faces = inward crossings). Each gives ≤ 1 fragment per pixel; a fragment
 * survives only if its proxy (uPiece) owns the exit of the last visible gas
 * interval (EXIT_* in chamber.ts), so exactly one fragment integrates the ray
 * and it is depth-tested at the gas exit, like the disc's single proxy. A metal
 * gap between two intervals stops the ray unless it lies in the engine's
 * cut-away region (uCut*; the piston is never cut).
 *
 * Emission is integrated with exact per-segment extinction:
 *   L += T · E/κ · (1 − e^{−κ Δ}),   T *= e^{−κ Δ}
 * Front profile: progress b = clamp(−x/δ, 0, 1) behind the (wrinkled) leading
 * edge x = |p − c| − r + a·n(p); heat release ∝ 6 b(1−b)/δ (∫ = 1), so the
 * face-on radiance of the front equals uFrontJ exactly, whatever δ.
 *
 * Output is premultiplied (blend ONE, ONE_MINUS_SRC_ALPHA).
 */

export const VOLUME_VERTEX = /* glsl */ `
uniform mat4 uMeshToCyl;
varying vec3 vCyl;
void main() {
  vCyl = (uMeshToCyl * vec4(position, 1.0)).xyz;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}
`;

export const VOLUME_FRAGMENT = /* glsl */ `
precision highp float;

varying vec3 vCyl;

uniform vec3 uCamPos;
uniform vec3 uCamDir;
uniform float uOrtho;
uniform float uR;
uniform float uH;
uniform vec4 uValve0;   // x, z, head radius, lift
uniform vec4 uValve1;
uniform vec2 uValveSeat0; // seat plane y, lift direction (−1 opens down into the bore, +1 up into the pocket)
uniform vec2 uValveSeat1;
uniform float uValveThick0;
uniform float uValveThick1;
uniform int uMode;      // 0 physical, 1 temperature

// chamber shape
uniform float uShape;      // 0 flat disc, 1 L-head
uniform vec4 uPocketRect;  // pocket plan: inner rectangle (corner-arc centres) x0, x1, z0, z1
uniform float uPocketRad;  // pocket plan corner radius
uniform vec2 uPocketY;     // deck (valve-seat) plane y, pocket roof y
uniform float uPiece;      // proxy drawn by this material: 0 bore column, 1 pocket prism, 2 transfer arc
uniform float uBackoff;    // L-head orthographic ray back-off
uniform float uCutOn;      // 1: the region {n·p > d for both planes} is cut away (engine cutaway)
uniform vec4 uCut0;
uniform vec4 uCut1;

// flame
uniform float uFlameOn;
uniform float uAllBurned;
uniform vec3 uFlameC;
uniform float uFlameR;
uniform float uBrush;       // δ
uniform float uWrinkleAmp;  // a
uniform float uWrinkleK;    // 1 / wavelength
uniform float uNoiseT;      // eddy-turnover phase

// physical emission
uniform vec3 uFrontJ;
uniform vec3 uBurnedJ;
uniform vec3 uEndGasJ;
uniform float uHazeK;
uniform vec3 uHazeC;
uniform float uHazeU;
uniform float uHazeCore;

// temperature mode
uniform float uTu;
uniform float uTb;
uniform vec2 uTRange;
uniform float uTempK;
uniform sampler2D uColormap;

// knock overlay
uniform float uKnockAmp;
uniform float uKnockTemporal;
uniform float uKnockRms;
uniform float uKnockAxis;
uniform float uKnockShape;  // 0 (1,0) Bessel mode of the bore, 1 axial mode along the footprint
uniform vec3 uKnockLine;    // footprint axis (ex, ez) and its start s0
uniform float uKnockLen;    // footprint length along the axis
uniform float uKnockSign;   // sign of the axial mode at the end-gas site
uniform vec2 uKnockOrigin;
uniform float uKnockRing;
uniform float uKnockRingW;
uniform float uKnockStanding;
uniform float uKnockK;
uniform vec3 uKnockPos;
uniform vec3 uKnockNeg;

uniform float uFrame;

const float ALPHA10 = 1.8412;
const float J1_ALPHA10 = 0.58187;

float hash13(vec3 p3) {
  p3 = fract(p3 * 0.1031);
  p3 += dot(p3, p3.zyx + 31.32);
  return fract((p3.x + p3.y) * p3.z);
}

float vnoise(vec3 x) {
  vec3 i = floor(x);
  vec3 f = fract(x);
  vec3 u = f * f * (3.0 - 2.0 * f);
  float n000 = hash13(i);
  float n100 = hash13(i + vec3(1.0, 0.0, 0.0));
  float n010 = hash13(i + vec3(0.0, 1.0, 0.0));
  float n110 = hash13(i + vec3(1.0, 1.0, 0.0));
  float n001 = hash13(i + vec3(0.0, 0.0, 1.0));
  float n101 = hash13(i + vec3(1.0, 0.0, 1.0));
  float n011 = hash13(i + vec3(0.0, 1.0, 1.0));
  float n111 = hash13(i + vec3(1.0, 1.0, 1.0));
  float nx00 = mix(n000, n100, u.x);
  float nx10 = mix(n010, n110, u.x);
  float nx01 = mix(n001, n101, u.x);
  float nx11 = mix(n011, n111, u.x);
  return mix(mix(nx00, nx10, u.y), mix(nx01, nx11, u.y), u.z) * 2.0 - 1.0;
}

// Wrinkling of the front by the large eddies (wavelength ≈ L_I) and a finer octave.
// Evolves with the eddy-turnover phase so slow motion shows the physical rate.
float wrinkle(vec3 p) {
  vec3 q = p * uWrinkleK;
  float t = uNoiseT;
  float n = 0.62 * vnoise(q + vec3(0.31 * t, -0.23 * t, 0.17 * t))
          + 0.38 * vnoise(2.37 * q + vec3(-0.41 * t, 0.29 * t, 0.37 * t) + 17.0);
  return clamp(n * 1.6, -1.0, 1.0);
}

float besselJ1(float x) {
  float x2 = x * x;
  return x * (0.5 - x2 * (1.0 / 16.0 - x2 * (1.0 / 384.0 - x2 * (1.0 / 18432.0 - x2 * (1.0 / 1474560.0 - x2 / 176947200.0)))));
}

// Knock pressure field (normalised), signed. Standing mode + expanding front from the end-gas site.
// Disc: Draper (1,0) Bessel mode. L-head: lowest axial mode cos(π s/L) of the footprint (geometry.ts).
float knockField(vec3 p) {
  float shape;
  if (uKnockShape > 0.5) {
    float s = clamp((dot(p.xz, uKnockLine.xy) - uKnockLine.z) / uKnockLen, 0.0, 1.0);
    shape = uKnockSign * cos(3.14159265 * s);
  } else {
    float r = min(length(p.xz) / uR, 1.0);
    float phi = atan(p.z, p.x);
    shape = besselJ1(ALPHA10 * r) / J1_ALPHA10 * cos(phi - uKnockAxis);
  }
  float standing = mix(shape * uKnockTemporal, abs(shape) * 0.7071, uKnockRms);
  float d = length(p.xz - uKnockOrigin);
  float ring = exp(-pow((d - uKnockRing) / uKnockRingW, 2.0)) * (1.0 - uKnockStanding);
  return uKnockAmp * (uKnockStanding * standing + ring);
}

vec2 raySphere(vec3 ro, vec3 rd, vec3 c, float r) {
  vec3 oc = ro - c;
  float b = dot(oc, rd);
  float cc = dot(oc, oc) - r * r;
  float h = b * b - cc;
  if (h < 0.0) return vec2(1e9, -1e9);
  h = sqrt(h);
  return vec2(-b - h, -b + h);
}

// Entry parameter of the ray into a valve head (short vertical cylinder), or 1e9. The head occupies
// the thickness just below the plane y = seatY + dir·lift (overhead: hanging below the lift plane;
// side valve: its top face, flush with the deck when shut, risen by the lift).
float rayValve(vec3 ro, vec3 rd, vec4 v, vec2 seat, float thick) {
  if (v.w <= 1e-5) return 1e9;
  vec2 o = ro.xz - v.xy;
  float a = dot(rd.xz, rd.xz);
  float b = 2.0 * dot(o, rd.xz);
  float c = dot(o, o) - v.z * v.z;
  float s0 = -1e9, s1 = 1e9;
  if (a > 1e-12) {
    float disc = b * b - 4.0 * a * c;
    if (disc < 0.0) return 1e9;
    float sq = sqrt(disc);
    s0 = (-b - sq) / (2.0 * a);
    s1 = (-b + sq) / (2.0 * a);
  } else if (c > 0.0) {
    return 1e9;
  }
  float y1 = seat.x + seat.y * v.w;
  float y0 = y1 - thick;
  float t0, t1;
  if (abs(rd.y) > 1e-12) {
    t0 = (y0 - ro.y) / rd.y;
    t1 = (y1 - ro.y) / rd.y;
    if (t0 > t1) { float tmp = t0; t0 = t1; t1 = tmp; }
  } else {
    if (ro.y < y0 || ro.y > y1) return 1e9;
    t0 = -1e9; t1 = 1e9;
  }
  float en = max(s0, t0);
  float ex = min(s1, t1);
  return en < ex && ex > 0.0 ? en : 1e9;
}

// Accumulators
vec3 accL;
float accT;

void addSegment(vec3 E, float kappa, float dt) {
  if (dt <= 0.0) return;
  float tau = kappa * dt;
  if (tau > 1e-4) {
    float tr = exp(-tau);
    accL += accT * E / kappa * (1.0 - tr);
    accT *= tr;
  } else {
    accL += accT * E * dt;
    accT *= 1.0 - tau;
  }
}

// Evaluate emission / extinction at p. shell: evaluate the brush; else b = uniformB.
void sampleAt(vec3 p, bool shell, float uniformB, out vec3 E, out float kappa) {
  float b = uniformB;
  float prof = 0.0;
  if (shell) {
    float d = length(p - uFlameC);
    float n = uWrinkleAmp > 0.0 ? wrinkle(p) : 0.0;
    float x = d - uFlameR + n * uWrinkleAmp;
    b = clamp(-x / uBrush, 0.0, 1.0);
    prof = 6.0 * b * (1.0 - b) / uBrush;
  }
  if (uMode == 1) {
    float T = mix(uTu, uTb, b);
    float u = clamp((T - uTRange.x) / (uTRange.y - uTRange.x), 0.0, 1.0);
    vec3 col = texture2D(uColormap, vec2(u, 0.5)).rgb;
    E = col * uTempK;
    kappa = uTempK;
  } else {
    float haze = uHazeK * mix(uHazeU, uHazeCore, b);
    E = uFrontJ * prof + uBurnedJ * b + uEndGasJ * (1.0 - b) + uHazeC * haze;
    kappa = haze;
  }
  if (uKnockAmp > 0.0) {
    float k = knockField(p);
    // isobar bands every 1/6 of full scale make the field readable as a pressure map
    float band = 0.6 + 0.4 * smoothstep(0.3, 0.5, abs(fract(abs(k) * 6.0) - 0.5));
    float s = abs(k) * uKnockK;
    E += (k >= 0.0 ? uKnockPos : uKnockNeg) * s * band;
    kappa += 0.5 * s;
  }
}

void march(vec3 ro, vec3 rd, float t0, float t1, int n, bool shell, float uniformB, float jitter) {
  float len = t1 - t0;
  if (len <= 0.0) return;
  float dt = len / float(n);
  for (int i = 0; i < 32; i++) {
    if (i >= n) break;
    float t = t0 + (float(i) + jitter) * dt;
    vec3 E; float kappa;
    sampleAt(ro + rd * t, shell, uniformB, E, kappa);
    addSegment(E, kappa, dt);
    if (accT < 0.004) break;
  }
}

// Emission along [ta, tb], split at the flame-brush shells (front-to-back into accL/accT).
void integrateInterval(vec3 ro, vec3 rd, float ta, float tb, float jitter) {
  float o0 = tb, o1 = tb, i0 = tb, i1 = tb;
  if (uAllBurned > 0.5) {
    o0 = ta; i0 = ta; i1 = tb; o1 = tb;
  } else if (uFlameOn > 0.5) {
    float rOut = uFlameR + uWrinkleAmp;
    float rIn = max(uFlameR - uBrush - uWrinkleAmp, 0.0);
    vec2 so = raySphere(ro, rd, uFlameC, rOut);
    if (so.x < so.y) {
      vec2 si = rIn > 0.0 ? raySphere(ro, rd, uFlameC, rIn) : vec2(1e9, -1e9);
      if (si.x >= si.y) { float m = 0.5 * (so.x + so.y); si = vec2(m, m); }
      o0 = clamp(so.x, ta, tb);
      i0 = clamp(si.x, ta, tb);
      i1 = clamp(si.y, ta, tb);
      o1 = clamp(so.y, ta, tb);
    }
  }
  int nSmooth = uKnockAmp > 0.0 ? 8 : 1;
  int nShell = 24;
  float outsideB = uAllBurned > 0.5 ? 1.0 : 0.0;
  march(ro, rd, ta, o0, nSmooth, false, outsideB, jitter);
  march(ro, rd, o0, i0, nShell, true, 0.0, jitter);
  march(ro, rd, i0, i1, nSmooth, false, 1.0, jitter);
  march(ro, rd, i1, o1, nShell, true, 0.0, jitter);
  march(ro, rd, o1, tb, nSmooth, false, outsideB, jitter);
}

// ---- L-head chamber (TS twin: chamber.ts rayChamberIntervals) ----

// Ray ∩ y-slab [y0, y1]; empty: x > y.
vec2 raySlabY(float oy, float dy, float y0, float y1) {
  if (abs(dy) > 1e-12) {
    float a = (y0 - oy) / dy, b = (y1 - oy) / dy;
    return vec2(min(a, b), max(a, b));
  }
  return (oy < y0 || oy > y1) ? vec2(1e9, -1e9) : vec2(-1e9, 1e9);
}

// 2-D ray (xz) ∩ circle; empty: x > y.
vec2 rayCircle2(vec2 o, vec2 d, vec2 c, float r) {
  vec2 q = o - c;
  float a = dot(d, d);
  float cc = dot(q, q) - r * r;
  if (a < 1e-12) return cc <= 0.0 ? vec2(-1e9, 1e9) : vec2(1e9, -1e9);
  float b = 2.0 * dot(q, d);
  float disc = b * b - 4.0 * a * cc;
  if (disc < 0.0) return vec2(1e9, -1e9);
  float sq = sqrt(disc);
  return vec2((-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a));
}

// 2-D ray ∩ axis-aligned box [lo, hi]; empty: x > y.
vec2 rayBox2(vec2 o, vec2 d, vec2 lo, vec2 hi) {
  vec2 r = vec2(-1e9, 1e9);
  if (abs(d.x) > 1e-12) {
    float a = (lo.x - o.x) / d.x, b = (hi.x - o.x) / d.x;
    r = vec2(max(r.x, min(a, b)), min(r.y, max(a, b)));
  } else if (o.x < lo.x || o.x > hi.x) {
    return vec2(1e9, -1e9);
  }
  if (abs(d.y) > 1e-12) {
    float a = (lo.y - o.y) / d.y, b = (hi.y - o.y) / d.y;
    r = vec2(max(r.x, min(a, b)), min(r.y, max(a, b)));
  } else if (o.y < lo.y || o.y > hi.y) {
    return vec2(1e9, -1e9);
  }
  return r;
}

vec2 hullOf(vec2 acc, vec2 r) {
  return r.x <= r.y ? vec2(min(acc.x, r.x), max(acc.y, r.y)) : acc;
}

// 2-D ray ∩ rounded rectangle (convex: hull of two crossed boxes and four corner discs).
vec2 rayRoundRect(vec2 o, vec2 d, vec4 inner, float rad) {
  vec2 acc = vec2(1e9, -1e9);
  acc = hullOf(acc, rayBox2(o, d, vec2(inner.x - rad, inner.z), vec2(inner.y + rad, inner.w)));
  acc = hullOf(acc, rayBox2(o, d, vec2(inner.x, inner.z - rad), vec2(inner.y, inner.w + rad)));
  if (rad > 0.0) {
    acc = hullOf(acc, rayCircle2(o, d, inner.xz, rad));
    acc = hullOf(acc, rayCircle2(o, d, inner.yz, rad));
    acc = hullOf(acc, rayCircle2(o, d, inner.xw, rad));
    acc = hullOf(acc, rayCircle2(o, d, inner.yw, rad));
  }
  return acc;
}

// Gas intervals (≤ 2, sorted) and the proxy owning each exit (0 bore, 1 pocket, 2 transfer arc).
float gS0, gE0, gK0, gS1, gE1, gK1;
int gN;

void pushIv(float s, float e, float k) {
  if (e <= s) return;
  if (gN == 1 && s <= gE0) { if (e > gE0) { gE0 = e; gK0 = k; } return; }
  if (gN == 2 && s <= gE1) { if (e > gE1) { gE1 = e; gK1 = k; } return; }
  if (gN == 0) { gS0 = s; gE0 = e; gK0 = k; gN = 1; }
  else if (gN == 1) { gS1 = s; gE1 = e; gK1 = k; gN = 2; }
}

void lheadIntervals(vec3 ro, vec3 rd) {
  gN = 0;
  gS0 = 0.0; gE0 = -1.0; gK0 = 0.0; gS1 = 0.0; gE1 = -1.0; gK1 = 0.0;
  vec2 C = rayCircle2(ro.xz, rd.xz, vec2(0.0), uR);
  bool cHit = C.x <= C.y;
  vec2 A = vec2(1e9, -1e9);
  if (cHit) {
    vec2 sy = raySlabY(ro.y, rd.y, -uH, 0.0);
    A = vec2(max(C.x, sy.x), min(C.y, sy.y));
  }
  vec2 P = rayRoundRect(ro.xz, rd.xz, uPocketRect, uPocketRad);
  vec2 sp = raySlabY(ro.y, rd.y, uPocketY.x, uPocketY.y);
  P = vec2(max(P.x, sp.x), min(P.y, sp.y));
  bool pHit = P.y > P.x;
  // pieces in ray order: pocket before the bore, bore column, pocket after the bore
  if (pHit) {
    if (!cHit) pushIv(P.x, P.y, 1.0);
    else pushIv(P.x, min(P.y, C.x), C.x < P.y ? 2.0 : 1.0);
  }
  pushIv(A.x, A.y, 0.0);
  if (pHit && cHit) pushIv(max(P.x, C.y), P.y, 1.0);
  // keep t >= 0
  if (gN > 0 && gE0 <= 0.0) { gS0 = gS1; gE0 = gE1; gK0 = gK1; gN -= 1; }
  if (gN > 0 && gE0 <= 0.0) gN = 0;
  gS0 = max(gS0, 0.0);
  gS1 = max(gS1, 0.0);
}

bool inCut(vec3 p) {
  return dot(uCut0.xyz, p) > uCut0.w && dot(uCut1.xyz, p) > uCut1.w;
}

// Metal between two gas intervals hides the second one unless it is cut away (never the piston).
bool gapOpen(vec3 ro, vec3 rd, float t0, float t1) {
  if (uCutOn < 0.5) return false;
  vec3 m = ro + rd * (0.5 * (t0 + t1));
  if (dot(m.xz, m.xz) < uR * uR && m.y < -uH) return false;
  return inCut(ro + rd * t0) && inCut(ro + rd * t1);
}

void main() {
  vec3 ro, rd;
  if (uOrtho > 0.5) {
    rd = normalize(uCamDir);
    ro = vCyl - rd * (uShape > 0.5 ? uBackoff : 4.0 * (uR + uH));
  } else {
    ro = uCamPos;
    rd = normalize(vCyl - uCamPos);
  }

  if (uShape > 0.5) {
    lheadIntervals(ro, rd);
    if (gN == 0) discard;
    bool both = gN == 2 && gapOpen(ro, rd, gE0, gS1);
    // exactly one proxy fragment per pixel: the one at the exit of the last visible interval
    if (abs((both ? gK1 : gK0) - uPiece) > 0.5) discard;
    float tv = min(rayValve(ro, rd, uValve0, uValveSeat0, uValveThick0),
                   rayValve(ro, rd, uValve1, uValveSeat1, uValveThick1));
    accL = vec3(0.0);
    accT = 1.0;
    float jitter = hash13(vec3(gl_FragCoord.xy, uFrame));
    float ta = gS0;
    float tb = min(gE0, max(tv, ta));
    if (tb > ta) integrateInterval(ro, rd, ta, tb, jitter);
    if (both && tv >= gE0 && accT >= 0.004) {
      ta = gS1;
      tb = min(gE1, max(tv, ta));
      if (tb > ta) integrateInterval(ro, rd, ta, tb, jitter);
    }
    gl_FragColor = vec4(accL, 1.0 - accT);
    #include <tonemapping_fragment>
    #include <colorspace_fragment>
    return;
  }

  // ---- chamber disc ----
  float ta = -1e9, tb = 1e9;
  float a = dot(rd.xz, rd.xz);
  float b = 2.0 * dot(ro.xz, rd.xz);
  float c = dot(ro.xz, ro.xz) - uR * uR;
  if (a > 1e-12) {
    float disc = b * b - 4.0 * a * c;
    if (disc < 0.0) discard;
    float sq = sqrt(disc);
    ta = (-b - sq) / (2.0 * a);
    tb = (-b + sq) / (2.0 * a);
  } else if (c > 0.0) {
    discard;
  }
  if (abs(rd.y) > 1e-12) {
    float y0 = (-uH - ro.y) / rd.y;
    float y1 = (0.0 - ro.y) / rd.y;
    ta = max(ta, min(y0, y1));
    tb = min(tb, max(y0, y1));
  } else if (ro.y > 0.0 || ro.y < -uH) {
    discard;
  }
  ta = max(ta, 0.0);
  // ---- open valve heads occlude the gas behind them ----
  tb = min(tb, max(rayValve(ro, rd, uValve0, uValveSeat0, uValveThick0), ta));
  tb = min(tb, max(rayValve(ro, rd, uValve1, uValveSeat1, uValveThick1), ta));
  if (tb <= ta) discard;

  accL = vec3(0.0);
  accT = 1.0;
  float jitter = hash13(vec3(gl_FragCoord.xy, uFrame));
  integrateInterval(ro, rd, ta, tb, jitter);

  gl_FragColor = vec4(accL, 1.0 - accT);
  #include <tonemapping_fragment>
  #include <colorspace_fragment>
}
`;
