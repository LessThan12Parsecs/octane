/**
 * GLSL for the in-cylinder gas volume: a single ray-marched pass over the
 * chamber disc. The ray is clipped ANALYTICALLY to the chamber cylinder
 * (x²+z² ≤ R², −h ≤ y ≤ 0) and to the open valve heads, then split at the
 * flame-brush shells (sphere radii r + a and r − δ − a) so the thin front is
 * always sampled densely while the uniform regions need 1 sample each (or a
 * few when the knock overlay is on).
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
uniform float uValveThick0;
uniform float uValveThick1;
uniform int uMode;      // 0 physical, 1 temperature

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

// Knock pressure field (normalised), signed. Standing (1,0) mode + expanding front from the end-gas site.
float knockField(vec3 p) {
  float r = min(length(p.xz) / uR, 1.0);
  float phi = atan(p.z, p.x);
  float shape = besselJ1(ALPHA10 * r) / J1_ALPHA10 * cos(phi - uKnockAxis);
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

// Entry parameter of the ray into a valve head (short vertical cylinder), or 1e9.
float rayValve(vec3 ro, vec3 rd, vec4 v, float thick) {
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
  float y0 = -v.w - thick, y1 = -v.w;
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

void main() {
  vec3 ro, rd;
  if (uOrtho > 0.5) {
    rd = normalize(uCamDir);
    ro = vCyl - rd * 4.0 * (uR + uH);
  } else {
    ro = uCamPos;
    rd = normalize(vCyl - uCamPos);
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
  tb = min(tb, max(rayValve(ro, rd, uValve0, uValveThick0), ta));
  tb = min(tb, max(rayValve(ro, rd, uValve1, uValveThick1), ta));
  if (tb <= ta) discard;

  accL = vec3(0.0);
  accT = 1.0;
  float jitter = hash13(vec3(gl_FragCoord.xy, uFrame));

  // ---- split at the flame brush shells ----
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

  gl_FragColor = vec4(accL, 1.0 - accT);
  #include <tonemapping_fragment>
  #include <colorspace_fragment>
}
`;
