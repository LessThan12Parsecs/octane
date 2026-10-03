"""Shared helpers for the CFR F-1 experimental validation dataset (tools/reference/cfr_data_*.py).

Not a generator itself: running it only prints the counter -> CR relation at a few points.

Contents
--------
* Paths (repository root, source-PDF folder, fixture folder).
* Unit constants (exact or NIST values; see comments).
* Digital-counter -> compression-ratio relations, replicated from src/physics/engines/cfr.ts:
    - `cr_rigid_raise(c)`  == cfrCompressionRatioAtCounter (rigid cylinder raise, 0.0007 in/digit,
      ASTM D2699-15a A2.2.1.1; one constant fitted to Choi et al. 2018 Fig. 4, see
      tools/reference/mechanics_cfr_conditions.py). THE project relation.
    - `cr_choi_polynomial(c)` == cfrCompressionRatioFromCounter (Choi et al. 2018, SAE 2018-01-0848, Fig. 4).
    - `cr_neste(c)`: CR = 6345/(1850 - DCR) + 1, the relation used on the Neste Oyj CFR engine
      (Bhattacharya et al., "Fuel-air mixing in motored CFR engine at RON relevant condition", OSTI 2536666, Sec. 2).
  The TS originals are cross-checked by test in the README (vitest snippet) - the Python copies must stay identical.
* Vector-path digitisation helpers (PyMuPDF drawing operators -> data coordinates, exact up to the PDF's
  coordinate quantisation).
* Raster digitisation helpers (colour masks, column-wise curve tracing, marker blobs) with explicit pixel
  error estimates.
* JSON writer.
"""
from __future__ import annotations

import json
import math
import pathlib
from dataclasses import dataclass
from typing import Callable, Iterable, Sequence

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[2]
SOURCES = ROOT / "tools" / "reference" / "cfr_sources"
FIXTURES = ROOT / "test" / "fixtures"

# ───────────────────────────── units ─────────────────────────────
IN = 0.0254  # m, exact (international inch)
PSI = 6894.757293168361  # Pa per lbf/in², exact from lb = 0.45359237 kg, g_n = 9.80665 m/s², in = 0.0254 m
IN_HG = 3386.389  # Pa per inch of mercury at 0 °C (NIST SP 811, Appendix B.8: 3.386 389 E+03)
MM_HG = 133.322387415  # Pa per mmHg (conventional, NIST SP 811 B.8: 1.333 224 E+02)
BAR = 1.0e5  # Pa, exact
MM_H2O = 9.80665  # Pa per mm of water (conventional, g_n * 1000 kg/m³ * 1 mm)
T0 = 273.15  # K at 0 °C, exact

# ─────────────────────── CFR geometry / counter ───────────────────────
STROKE = 4.5 * IN  # m, ASTM D2699-15a Table 1 ("4.50 in"); same constant as src/physics/engines/cfr.ts
BORE = 3.25 * IN  # m, ASTM D2699-15a Table 1
CONROD = 0.254  # m, Pal et al. 2018 Table 2 / Choi et al. 2018 Table 1
COUNTER_STEP = 0.0007 * IN  # m per digit, ASTM D2699-15a A2.2.1.1 (= cfr.ts CFR_COUNTER_STEP)
BASIC_COUNTER = 930.0  # ASTM D2699-15a 10.3.17.1 / A2.2.2; D2700-14 10.3.18.1 (= cfr.ts CFR_BASIC_COUNTER)
# Effective clearance height at counter 930, m. Copied verbatim from src/physics/engines/cfr.ts
# (CFR_CLEARANCE_HEIGHT_AT_BASIC_COUNTER), itself fitted in tools/reference/mechanics_cfr_conditions.py to the
# oil-measured CR of Choi et al. 2018 Fig. 4.
H_930 = 0.01742116238087605


def cr_rigid_raise(counter: float) -> float:
    """cfr.ts cfrCompressionRatioAtCounter: CR = 1 + S / (h930 - 0.0007 in (c - 930))."""
    return 1.0 + STROKE / (H_930 - COUNTER_STEP * (counter - BASIC_COUNTER))


def counter_rigid_raise(cr: float) -> float:
    """cfr.ts cfrCounterAtCompressionRatio (closed-form inverse)."""
    return BASIC_COUNTER + (H_930 - STROKE / (cr - 1.0)) / COUNTER_STEP


def cr_choi_polynomial(counter: float) -> float:
    """cfr.ts cfrCompressionRatioFromCounter: Choi et al. 2018 Fig. 4 (valid 400-1400)."""
    c = counter
    return ((1.126e-8 * c - 2.126e-5) * c + 1.694e-2) * c + 1.024


def cr_neste(counter: float) -> float:
    """Neste Oyj CFR engine: CR = 6345/(1850 - DCR) + 1 (OSTI 2536666, Sec. 2)."""
    return 6345.0 / (1850.0 - counter) + 1.0


def counter_from_dial(dial_in: float) -> float:
    """ASTM D2699-15a Table A4.3 footnote A: equivalent digital counter = (1.012 - dial indicator) x 1410."""
    return (1.012 - dial_in) * 1410.0


def cr_all(counter: float | None) -> dict | None:
    if counter is None:
        return None
    return {
        "rigidRaise": cr_rigid_raise(counter),
        "choiPolynomial": cr_choi_polynomial(counter) if 400 <= counter <= 1400 else None,
    }


# ─────────────────────────── JSON output ───────────────────────────
def _round(obj, sig: int = 7):
    if isinstance(obj, float):
        if not math.isfinite(obj):
            return None
        if obj == 0.0:
            return 0.0
        return float(f"{obj:.{sig}g}")
    if isinstance(obj, (np.floating,)):
        return _round(float(obj), sig)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return [_round(v, sig) for v in obj.tolist()]
    if isinstance(obj, dict):
        return {k: _round(v, sig) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_round(v, sig) for v in obj]
    return obj


def write_fixture(name: str, payload: dict, sig: int = 7) -> pathlib.Path:
    """Write test/fixtures/<name> (compact JSON, floats rounded to `sig` significant digits)."""
    out = FIXTURES / name
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(_round(payload, sig), separators=(",", ":"), ensure_ascii=False) + "\n")
    print(f"wrote {out.relative_to(ROOT)} ({out.stat().st_size / 1024:.1f} kB)")
    return out


def source_pdf(name: str) -> pathlib.Path:
    p = SOURCES / name
    if not p.exists():
        raise FileNotFoundError(f"{p} missing - see test/fixtures/cfr_validation_README.md for the download URL")
    return p


# ─────────────────────────── conditions ───────────────────────────
def conditions(**kw) -> dict:
    """Operating conditions in SI with explicit nulls. Keys (all optional):
    engine, method ('RON'|'MON'|'RON-like'|'HCCI'|'motored'|...), fuel (dict), compressionRatio,
    compressionRatioSource, digitalCounter, rpm, sparkAdvanceDegBTDC, lambda, phi, intakeAirTemperatureK,
    intakeMixtureTemperatureK, intakePressurePa, barometricPressurePa, exhaustBackPressurePa,
    coolantTemperatureK, oilTemperatureK, humidityRatio, knockIntensityKU, notes."""
    keys = [
        "engine", "method", "fuel", "compressionRatio", "compressionRatioSource", "digitalCounter", "rpm",
        "sparkAdvanceDegBTDC", "lambda", "phi", "intakeAirTemperatureK", "intakeMixtureTemperatureK",
        "intakePressurePa", "barometricPressurePa", "exhaustBackPressurePa", "coolantTemperatureK",
        "oilTemperatureK", "humidityRatio", "knockIntensityKU", "notes",
    ]
    unknown = set(kw) - set(keys)
    if unknown:
        raise KeyError(f"unknown condition keys {unknown}")
    out = {k: kw.get(k) for k in keys}
    if out["lambda"] is not None and out["phi"] is None:
        out["phi"] = 1.0 / out["lambda"]
    if out["phi"] is not None and out["lambda"] is None:
        out["lambda"] = 1.0 / out["phi"]
    return out


def prf(on: float) -> dict:
    """Primary reference fuel: ON = liquid-volume % iso-octane in n-heptane (ASTM D2699-15a 3.1.19.1)."""
    return {"kind": "PRF", "octaneNumber": on, "label": f"PRF{on:g}"}


# ─────────────────────── vector digitisation ───────────────────────
@dataclass
class Axis:
    """Linear map page-coordinate -> data value, least-squares through (coord, value) pairs."""

    a: float
    b: float
    resid: float  # max |residual| of the calibration points, in data units
    log: bool = False

    @classmethod
    def fit(cls, coords: Sequence[float], values: Sequence[float], log: bool = False) -> "Axis":
        c = np.asarray(coords, float)
        v = np.log10(np.asarray(values, float)) if log else np.asarray(values, float)
        A = np.vstack([c, np.ones_like(c)]).T
        (a, b), *_ = np.linalg.lstsq(A, v, rcond=None)
        r = float(np.max(np.abs(A @ np.array([a, b]) - v))) if len(c) > 2 else 0.0
        return cls(float(a), float(b), r, log)

    def __call__(self, x):
        v = self.a * np.asarray(x, float) + self.b
        return 10.0 ** v if self.log else v

    def per_unit(self) -> float:
        """data units per page unit (|slope|)."""
        return abs(self.a)


def drawing_points(drawing: dict, control_points: bool = False) -> np.ndarray:
    """(N, 2) page coordinates of a PyMuPDF drawing path: segment end points in drawing order
    (Bezier control points only if requested). Consecutive duplicates removed."""
    pts: list[tuple[float, float]] = []
    for it in drawing["items"]:
        kind = it[0]
        if kind == "l":
            seq = [it[1], it[2]]
        elif kind == "c":
            seq = [it[1], it[2], it[3], it[4]] if control_points else [it[1], it[4]]
        elif kind == "re":
            r = it[1]
            seq = [pymupdf_point(r.x0, r.y0), pymupdf_point(r.x1, r.y1)]
        elif kind == "qu":
            q = it[1]
            seq = [q.ul, q.ur, q.lr, q.ll]
        else:
            continue
        for p in seq:
            xy = (float(p.x), float(p.y))
            if not pts or (abs(pts[-1][0] - xy[0]) > 1e-9 or abs(pts[-1][1] - xy[1]) > 1e-9):
                pts.append(xy)
    return np.asarray(pts, float)


def pymupdf_point(x, y):
    import pymupdf

    return pymupdf.Point(x, y)


def page_drawings(pdf: pathlib.Path, page_no: int) -> list[dict]:
    import pymupdf

    doc = pymupdf.open(pdf)
    return doc[page_no - 1].get_drawings()


def page_words(pdf: pathlib.Path, page_no: int, clip: tuple[float, float, float, float] | None = None) -> list[tuple]:
    import pymupdf

    doc = pymupdf.open(pdf)
    words = doc[page_no - 1].get_text("words")
    if clip is None:
        return words
    x0, y0, x1, y1 = clip
    return [w for w in words if w[0] >= x0 and w[1] >= y0 and w[2] <= x1 and w[3] <= y1]


def tick_positions(drawing: dict, axis: str) -> list[float]:
    """Positions of short tick segments in a drawing ('x' -> x coordinates of vertical ticks, 'y' -> y of horizontal)."""
    out = []
    for it in drawing["items"]:
        if it[0] != "l":
            continue
        p, q = it[1], it[2]
        out.append(float(p.x) if axis == "x" else float(p.y))
    return out


def marker_centres(drawings: Iterable[dict]) -> np.ndarray:
    """Centres of marker paths (bounding-box centres)."""
    c = []
    for d in drawings:
        r = d["rect"]
        c.append(((r.x0 + r.x1) / 2.0, (r.y0 + r.y1) / 2.0))
    return np.asarray(c, float)


def colour_eq(c1, c2, tol: float = 0.01) -> bool:
    if c1 is None or c2 is None:
        return c1 is c2
    return all(abs(a - b) <= tol for a, b in zip(c1, c2))


# ─────────────────────── raster digitisation ───────────────────────
def load_image(pdf: pathlib.Path, xref: int) -> np.ndarray:
    """RGB uint8 array (H, W, 3) of an embedded image (native resolution, no resampling)."""
    import pymupdf

    doc = pymupdf.open(pdf)
    pix = pymupdf.Pixmap(doc, xref)
    if pix.alpha or pix.n - pix.alpha != 3:
        pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)
    return arr[:, :, :3].copy()


def colour_mask(img: np.ndarray, rgb: Sequence[int], tol: float) -> np.ndarray:
    """Pixels within Euclidean RGB distance `tol` of `rgb`."""
    d = np.sqrt(np.sum((img.astype(float) - np.asarray(rgb, float)) ** 2, axis=2))
    return d <= tol


def trace_curve(mask: np.ndarray, x_range: tuple[int, int] | None = None, y_range: tuple[int, int] | None = None,
                mode: str = "centre", min_run: int = 1) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Column-wise curve trace of a boolean mask. For every pixel column with mask pixels (optionally in
    x/y windows): returns (x_px, y_px, half_width_px) where y is the centre (mean) of the LARGEST vertical
    run of mask pixels ('centre'), or its top/bottom edge. half_width is half the run length (the
    line-width contribution to the digitisation uncertainty)."""
    H, W = mask.shape
    x0, x1 = x_range if x_range else (0, W)
    y0, y1 = y_range if y_range else (0, H)
    xs, ys, hw = [], [], []
    for x in range(x0, x1):
        col = np.nonzero(mask[y0:y1, x])[0]
        if col.size == 0:
            continue
        # split into runs
        splits = np.nonzero(np.diff(col) > 1)[0]
        runs = np.split(col, splits + 1)
        run = max(runs, key=len)
        if len(run) < min_run:
            continue
        if mode == "centre":
            yc = 0.5 * (run[0] + run[-1])
        elif mode == "top":
            yc = run[0]
        else:
            yc = run[-1]
        xs.append(x)
        ys.append(yc + y0)
        hw.append(0.5 * (run[-1] - run[0] + 1))
    return np.asarray(xs, float), np.asarray(ys, float), np.asarray(hw, float)


def blobs(mask: np.ndarray, min_area: int = 4, max_area: int | None = None) -> list[dict]:
    """Connected components (8-connectivity) with centroid, bbox, area, fill ratio and centroid offset
    within the bbox (for marker-shape classification)."""
    from scipy import ndimage

    lab, n = ndimage.label(mask, structure=np.ones((3, 3), int))
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab), start=1):
        if sl is None:
            continue
        sub = lab[sl] == i
        area = int(sub.sum())
        if area < min_area or (max_area is not None and area > max_area):
            continue
        yy, xx = np.nonzero(sub)
        h, w = sub.shape
        cy, cx = yy.mean(), xx.mean()
        out.append(
            dict(
                cx=float(cx + sl[1].start), cy=float(cy + sl[0].start), area=area, w=w, h=h,
                x0=sl[1].start, y0=sl[0].start, fill=area / float(w * h),
                # centroid offset from bbox centre, normalised by bbox size (triangles are asymmetric)
                dx=float((cx - (w - 1) / 2.0) / w), dy=float((cy - (h - 1) / 2.0) / h),
            )
        )
    return out


def find_lines(mask: np.ndarray, axis: int, min_frac: float = 0.6) -> list[float]:
    """Positions (pixel centres) of long straight lines in a mask: axis=0 -> horizontal lines (returns y),
    axis=1 -> vertical lines (returns x). A row/column qualifies if >= min_frac of it is set; adjacent
    rows/columns are merged."""
    frac = mask.mean(axis=1 - axis) if axis == 0 else mask.mean(axis=0)
    idx = np.nonzero(frac >= min_frac)[0]
    if idx.size == 0:
        return []
    groups = np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)
    return [float(g.mean()) for g in groups]


def resample_uniform(x: np.ndarray, y: np.ndarray, step: float, x0: float | None = None, x1: float | None = None):
    """Linear resampling of a monotone-x curve onto a uniform grid (for compact fixtures)."""
    order = np.argsort(x)
    x, y = x[order], y[order]
    # merge duplicate x
    ux, inv = np.unique(x, return_inverse=True)
    uy = np.bincount(inv, weights=y) / np.bincount(inv)
    a = math.ceil((ux[0] if x0 is None else max(x0, ux[0])) / step - 1e-9) * step
    b = math.floor((ux[-1] if x1 is None else min(x1, ux[-1])) / step + 1e-9) * step
    grid = np.round(np.arange(a, b + 0.5 * step, step), 6)
    return grid, np.interp(grid, ux, uy)


if __name__ == "__main__":
    for c in (264, 726, 749, 930, 919, 1145):
        print(c, round(cr_rigid_raise(c), 4), round(cr_choi_polynomial(c), 4), round(cr_neste(c), 4))


# ─────────────────────── vector chart helpers ───────────────────────
_NUM_RE = __import__("re").compile(r"^[-−–]?\d+(\.\d+)?$")


def numeric_words(words: list[tuple], region: tuple[float, float, float, float]) -> list[tuple[float, float, float]]:
    """(x_centre, y_centre, value) of the numeric words whose centre lies inside region (x0, y0, x1, y1)."""
    out = []
    x0, y0, x1, y1 = region
    for w in words:
        t = w[4].replace("−", "-").replace("–", "-")
        if not _NUM_RE.match(w[4]):
            continue
        xc, yc = 0.5 * (w[0] + w[2]), 0.5 * (w[1] + w[3])
        if x0 <= xc <= x1 and y0 <= yc <= y1:
            out.append((xc, yc, float(t)))
    return out


def axis_from_labels(words: list[tuple], region: tuple[float, float, float, float], orient: str, log: bool = False) -> Axis:
    """Axis calibrated on the centres of the numeric tick labels found in `region` ('x' uses the label x centre,
    'y' the label y centre). Tick labels are centred on their ticks in Excel/matplotlib output, so the residual of
    the fit (Axis.resid) measures the label-centring error (typically < 0.3 pt)."""
    nw = numeric_words(words, region)
    assert len(nw) >= 2, (region, nw)
    c = [n[0] if orient == "x" else n[1] for n in nw]
    v = [n[2] for n in nw]
    return Axis.fit(c, v, log=log)


def short_segment_markers(drawings: list[dict], colour, max_len: float = 12.0,
                          exclude: list[tuple[float, float, float, float]] = (), gap: float = 0.8) -> np.ndarray:
    """Markers drawn as bundles of separate short line segments (Excel 'x', '+', '*' markers in some PDF exports):
    collect 1-segment (or 2-segment) stroked paths of `colour` shorter than max_len, cluster those whose bounding
    boxes touch (within `gap`), return cluster bounding-box centres. Clusters inside `exclude` boxes are dropped."""
    segs = []
    for d in drawings:
        if d["type"] not in ("s", "fs") or not colour_eq(d.get("color"), colour, 0.02):
            continue
        if not all(it[0] == "l" for it in d["items"]) or len(d["items"]) > 2:
            continue
        r = d["rect"]
        if max(r.width, r.height) > max_len:
            continue
        segs.append([r.x0, r.y0, r.x1, r.y1])
    boxes: list[list[float]] = []
    for s in segs:
        merged = False
        for b in boxes:
            if s[0] <= b[2] + gap and s[2] >= b[0] - gap and s[1] <= b[3] + gap and s[3] >= b[1] - gap:
                b[0], b[1], b[2], b[3] = min(b[0], s[0]), min(b[1], s[1]), max(b[2], s[2]), max(b[3], s[3])
                merged = True
                break
        if not merged:
            boxes.append(list(s))
    # second pass: merge boxes that now overlap
    changed = True
    while changed:
        changed = False
        for i in range(len(boxes)):
            for j in range(i + 1, len(boxes)):
                a, b = boxes[i], boxes[j]
                if a[0] <= b[2] + gap and a[2] >= b[0] - gap and a[1] <= b[3] + gap and a[3] >= b[1] - gap:
                    boxes[i] = [min(a[0], b[0]), min(a[1], b[1]), max(a[2], b[2]), max(a[3], b[3])]
                    boxes.pop(j)
                    changed = True
                    break
            if changed:
                break
    cen = []
    for b in boxes:
        cx, cy = 0.5 * (b[0] + b[2]), 0.5 * (b[1] + b[3])
        if any(e[0] <= cx <= e[2] and e[1] <= cy <= e[3] for e in exclude):
            continue
        cen.append((cx, cy))
    return np.asarray(sorted(cen), float)


def legend_names(drawings: list[dict], words: list[tuple], region: tuple[float, float, float, float],
                 pattern: str = r"^PRF\d+") -> dict[str, str]:
    """Map line style (dashes string, or colour) -> legend label for the 1-segment horizontal legend sample lines in
    `region`: the label is the matching word immediately to the right on the same text line."""
    import re as _re

    out: dict[str, str] = {}
    x0, y0, x1, y1 = region
    for d in drawings:
        if d["type"] != "s" or len(d["items"]) != 1 or d["items"][0][0] != "l":
            continue
        r = d["rect"]
        if not (x0 <= r.x0 and r.x1 <= x1 and y0 <= r.y0 and r.y1 <= y1) or r.height > 0.5 or r.width < 8:
            continue
        yc = 0.5 * (r.y0 + r.y1)
        cands = [w for w in words if _re.match(pattern, w[4]) and abs(0.5 * (w[1] + w[3]) - yc) < 3.5
                 and 0 <= w[0] - r.x1 < 12]
        if cands:
            out[style_key(d)] = min(cands, key=lambda w: w[0] - r.x1)[4]
    return out


def style_key(d: dict) -> str:
    col = tuple(round(c, 2) for c in (d.get("color") or ()))
    dash = (d.get("dashes") or "").strip()
    # dash lengths differ slightly between panels (scaled charts): keep the pattern structure, rounded to 0.5 pt
    import re as _re

    nums = [float(v) for v in _re.findall(r"[\d.]+", dash.split("]")[0])] if "[" in dash else []
    norm = tuple(round(v * 2) / 2 for v in nums)
    return f"{col}|{norm}"
