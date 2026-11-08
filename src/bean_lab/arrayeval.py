"""Evaluation of planar array layouts over a set of listener positions (foci).

Conventions: millimeters, panel centered at the origin, array in z = 0, normal +z. A listener position
is given from the panel center by distance R, azimuth alpha (toward +x) and elevation beta (toward +y):

    d = (sin(alpha) cos(beta), sin(beta), cos(alpha) cos(beta)),     F = R d

Power normalization. Every layout emits the same total power: each speaker has amplitude
A_N = A_REF sqrt(N_REF / N), so N A_N^2 / 2 is constant (A_REF = 0.05 at N_REF = 8). Levels are dB re one
speaker of amplitude A_REF at 1 m, so a coherent focus gains 10 log10(N / 8) over the "raw" equal-amplitude
reading (reported as `level_raw`).

Region of interest (ROI), common to every layout and focus: around each focus, directions in a local
frame (axis = the focus direction) within +/- half_az x +/- half_el degrees, on shells at fractions of
the focus distance; only directions within `max_normal_deg` of the panel normal are kept. Energies are
compensated for spherical spreading, E_c = E (r / R)^2, and weighted by the solid angle (cos of the local
elevation). Region statistics are therefore dimensionless shares of ROI energy: they are NOT total
radiated acoustic power.

Metrics per focus and band (energy sum over spectral lines, see broadband.py):
  * F5/F10/F15: share of ROI energy within 5/10/15 degrees (great circle) of the focus direction;
  * PSL_cap: focus level minus the highest compensated level outside 12 degrees (contrast to the
    secondary lobes); PSL_lobe (reference only): same, outside the dilated -6 dB main-lobe component;
  * widths at -3 dB along the local horizontal/vertical, areas of the -3/-6 dB components, depth along
    the focus direction, pedestal: all flagged `censored` when the region reaches the ROI limits
    (a lower bound, never an exact value);
  * line robustness: 10th percentile over the spectral lines of F10 and PSL_cap.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

from .arraygen import to_array_config
from .broadband import Spectrum, connected_component
from .campaign import Params
from .field import FieldGeometry
from .field3d import dilate, local_maxima
from .focus import SPEED_OF_SOUND, focus_delays

A_REF = 0.05
N_REF = 8
CAP_DEG = 10.0
EXCLUSION_DEG = 12.0
WINDOWS_DEG = (5.0, 10.0, 15.0)


def element_gain(n: int) -> float:
    """Per-speaker amplitude keeping the total emitted power equal to that of N_REF speakers of A_REF."""
    if n < 1:
        raise ValueError("n must be positive")
    return A_REF * math.sqrt(N_REF / n)


@dataclass(frozen=True)
class Focus:
    name: str
    role: str  # "train" | "val"
    distance_mm: float
    az_deg: float
    el_deg: float

    @property
    def direction(self) -> np.ndarray:
        a, b = math.radians(self.az_deg), math.radians(self.el_deg)
        return np.array([math.sin(a) * math.cos(b), math.sin(b), math.cos(a) * math.cos(b)])

    @property
    def position_mm(self) -> np.ndarray:
        return self.distance_mm * self.direction


def _f(name, role, r, az, el):
    return Focus(name, role, float(r), float(az), float(el))


# Four distances, five azimuths and three elevations are covered by each set; the validation set is
# never seen by the optimizer. Distances 600/800/1200/1500 mm, az 0/+-15/+-30, el 0/+-15 (mission grid).
TRAIN_FOCI = (_f("T1", "train", 800, 0, 0), _f("T2", "train", 600, 15, 0), _f("T3", "train", 1200, -15, 15),
              _f("T4", "train", 1500, 30, 0), _f("T5", "train", 800, -30, -15))
VALIDATION_FOCI = (_f("V1", "val", 600, -15, 15), _f("V2", "val", 800, 30, 15), _f("V3", "val", 1200, 15, 0),
                   _f("V4", "val", 1500, -30, -15), _f("V5", "val", 1200, 0, -15))
ALL_FOCI = TRAIN_FOCI + VALIDATION_FOCI

BANDS = {"b1535": (1500.0, 3500.0), "b24": (2000.0, 4000.0), "b18": (1000.0, 8000.0)}


@dataclass(frozen=True)
class RoiSpec:
    half_az: float = 60.0
    half_el: float = 40.0
    step_deg: float = 4.0  # angular step of the middle shell
    shells: tuple[float, ...] = (0.85, 1.0, 1.15)
    outer_stride: int = 2  # outer shells use every outer_stride-th direction (weights compensate)
    max_normal_deg: float = 85.0


FAST_ROI = RoiSpec()
REFERENCE_ROI = RoiSpec(step_deg=2.0, shells=(0.8, 0.9, 1.0, 1.1, 1.2))


class FocusRoi:
    """Sample points around one focus. Point 0 is the focus itself (weight 0)."""

    def __init__(self, focus: Focus, spec: RoiSpec, center_mm=(0.0, 0.0, 0.0)):
        self.focus, self.spec = focus, spec
        self.center = np.asarray(center_mm, dtype=float)
        d = focus.direction
        e_h = np.cross([0.0, 1.0, 0.0], d)
        e_h /= np.linalg.norm(e_h)
        e_v = np.cross(d, e_h)
        a = np.arange(-spec.half_az, spec.half_az + 1e-9, spec.step_deg)
        b = np.arange(-spec.half_el, spec.half_el + 1e-9, spec.step_deg)
        self.a_deg, self.b_deg = a, b
        ia0, ib0 = len(a) // 2, len(b) // 2
        pts, wts, gam, r2, shell, cell = [focus.position_mm], [0.0], [0.0], [1.0], [-1], [(-1, -1)]
        R = focus.distance_mm
        self.mid_grid = -np.ones((len(b), len(a)), dtype=int)
        for s_idx, s in enumerate(spec.shells):
            stride = 1 if math.isclose(s, 1.0) else spec.outer_stride
            for ib in range(len(b)):
                if (ib - ib0) % stride:
                    continue
                for ia in range(len(a)):
                    if (ia - ia0) % stride:
                        continue
                    A_, B_ = math.radians(a[ia]), math.radians(b[ib])
                    u = math.cos(B_) * (math.cos(A_) * d + math.sin(A_) * e_h) + math.sin(B_) * e_v
                    if u[2] < math.cos(math.radians(spec.max_normal_deg)):
                        continue
                    p = self.center + s * R * u
                    if stride == 1 and math.isclose(s, 1.0):
                        self.mid_grid[ib, ia] = len(pts)
                    pts.append(p)
                    wts.append(math.cos(B_) * stride * stride)
                    gam.append(math.degrees(math.acos(max(-1.0, min(1.0, math.cos(B_) * math.cos(A_))))))
                    r2.append((np.linalg.norm(p - self.center) / R) ** 2)
                    shell.append(s_idx)
                    cell.append((ib, ia))
        self.pts_mm = np.array(pts)
        self.pts_m = self.pts_mm / 1000.0
        self.w = np.array(wts)
        self.gamma = np.array(gam)
        self.r2 = np.array(r2)
        self.shell = np.array(shell)
        self.cell = np.array(cell)
        self.ia0, self.ib0 = ia0, ib0

    @property
    def n(self) -> int:
        return len(self.pts_mm)


def core_metrics(E_lines: np.ndarray, line_w: np.ndarray, roi: FocusRoi, excl_deg: float = EXCLUSION_DEG) -> dict:
    """Shares and contrast from per-line energies E_lines (K, n_points) at the ROI points (point 0 = focus)."""
    comp = roi.r2
    w = roi.w[1:]
    ec_lines = E_lines[:, 1:] * comp[1:]
    ec = line_w @ ec_lines
    tot = float(ec @ w)
    out = {}
    for win in WINDOWS_DEG:
        sel = roi.gamma[1:] <= win
        out[f"F{win:g}"] = float((ec * w)[sel].sum() / tot)
    sel10 = roi.gamma[1:] <= CAP_DEG
    totals_k = ec_lines @ w
    f10_k = (ec_lines * w)[:, sel10].sum(axis=1) / totals_k
    outside = roi.gamma[1:] > excl_deg
    e_focus_band = float(line_w @ E_lines[:, 0])
    psl_k = 10 * np.log10(E_lines[:, 0] / ec_lines[:, outside].max(axis=1))
    out.update(psl_cap=float(10 * math.log10(e_focus_band / ec[outside].max())),
               F10_lines=f10_k, psl_lines=psl_k, e_focus=e_focus_band, ec=ec, total_energy=tot)
    return out


@dataclass(frozen=True)
class Weights:
    out: float = 1.0  # mean energy outside the desired cap
    psl: float = 0.5  # mean normalized sidelobe penalty
    worst: float = 1.0  # worst focus of the outside energy
    spec: float = 0.5  # lower spectral percentile of the in-cap energy
    psl_target_db: float = 15.0
    line_percentile: float = 10.0

    def __post_init__(self) -> None:
        if min(self.out, self.psl, self.worst, self.spec) < 0 or self.psl_target_db <= 0:
            raise ValueError("weights must be nonnegative and psl_target_db positive")


def objective_components(per_focus: list[dict], weights: Weights) -> dict:
    """Normalized components (all in [0, 1], lower is better) and their weighted sum."""
    a = np.array([1.0 - m["F10"] for m in per_focus])
    b = np.clip((weights.psl_target_db - np.array([m["psl_cap"] for m in per_focus])) / weights.psl_target_db, 0, 1)
    c = np.array([1.0 - np.percentile(m["F10_lines"], weights.line_percentile) for m in per_focus])
    comp = dict(J_out=float(a.mean()), J_psl=float(b.mean()), J_worst=float(a.max()), J_spec=float(c.mean()))
    comp["J"] = (weights.out * comp["J_out"] + weights.psl * comp["J_psl"] + weights.worst * comp["J_worst"]
                 + weights.spec * comp["J_spec"])
    return comp


class FastEvaluator:
    """Incremental evaluator for the optimizer: caches the complex field sum of the whole layout.

    Moving element i only changes its own contribution. The focusing delays e_i = (d_max - d_i) / c enter
    the energy only through a phase that is common to all elements when d_max changes, so the cached
    contribution of element i at point x, line f and focus F is
        (gain / r) exp(-j 2 pi f (r - |F - p_i|) / c),
    depending on p_i alone. `reference_energies` below recomputes the same energies with the real delays
    and FieldGeometry; the two agree to rounding (tested).
    """

    def __init__(self, positions_mm: np.ndarray, foci, rois: list[FocusRoi], freqs: np.ndarray, line_w: np.ndarray,
                 weights: Weights = Weights(), speed_of_sound: float = SPEED_OF_SOUND):
        self.pos = np.column_stack([np.asarray(positions_mm, dtype=float), np.zeros(len(positions_mm))])
        self.n = len(self.pos)
        self.foci, self.rois = list(foci), rois
        self.freqs, self.line_w = np.asarray(freqs, dtype=float), np.asarray(line_w, dtype=float)
        self.weights, self.c = weights, speed_of_sound
        self.gain = element_gain(self.n)
        self.omega = 2 * np.pi * self.freqs
        self.evaluations = 0
        self.totals = [self._total(j) for j in range(len(self.rois))]

    def contribution(self, p_mm: np.ndarray, j: int) -> np.ndarray:
        roi, focus = self.rois[j], self.foci[j]
        r = np.linalg.norm(roi.pts_mm - p_mm, axis=1) / 1000.0
        d_f = np.linalg.norm(focus.position_mm - p_mm) / 1000.0
        phase = np.exp(-1j * self.omega[:, None] * (r[None, :] - d_f) / self.c)
        return (self.gain / r)[None, :] * phase

    def _total(self, j: int) -> np.ndarray:
        acc = np.zeros((len(self.freqs), self.rois[j].n), dtype=complex)
        for i in range(self.n):
            acc += self.contribution(self.pos[i], j)
        return acc

    def _metrics(self, totals) -> list[dict]:
        return [core_metrics(np.abs(t) ** 2, self.line_w, roi) for t, roi in zip(totals, self.rois)]

    def objective(self, totals=None) -> dict:
        self.evaluations += 1
        per_focus = self._metrics(self.totals if totals is None else totals)
        comp = objective_components(per_focus, self.weights)
        comp["per_focus"] = per_focus
        return comp

    def try_move(self, i: int, new_mm: np.ndarray):
        """(objective components, pending totals) if element i moved to new_mm (not committed)."""
        new3 = np.array([new_mm[0], new_mm[1], 0.0])
        totals = [t - self.contribution(self.pos[i], j) + self.contribution(new3, j) for j, t in enumerate(self.totals)]
        return self.objective(totals), totals, new3

    def commit(self, i: int, totals, new3) -> None:
        self.totals, self.pos[i] = totals, new3

    def refresh(self) -> None:
        """Recompute the cached totals from scratch (removes accumulated rounding)."""
        self.totals = [self._total(j) for j in range(len(self.rois))]


def fast_setup(positions_mm, foci, band: tuple[float, float], n_lines: int, spec: RoiSpec = FAST_ROI,
               weights: Weights = Weights(), distribution: str = "uniform") -> FastEvaluator:
    freqs, w = Spectrum(band[0], band[1], n=n_lines, distribution=distribution).lines()
    rois = [FocusRoi(f, spec) for f in foci]
    return FastEvaluator(positions_mm, foci, rois, freqs, w, weights)


# ---- reference evaluation ----------------------------------------------------------------------------------

def reference_energies(array, focus: Focus, roi: FocusRoi, freqs: np.ndarray, speed_of_sound: float = SPEED_OF_SOUND,
                       block: int = 20000):
    """(K, n_points) energies with the real focusing delays (FieldGeometry), amplitude element_gain(N).

    Points are processed in blocks of at most `block` so memory stays bounded for very fine ROIs.
    """
    n = len(array.elements)
    emission = {d.speaker: d.delay for d in focus_delays(array, tuple(float(v) for v in focus.position_mm / 1000),
                                                         speed_of_sound)}
    out = np.empty((len(freqs), roi.n))
    for lo in range(0, roi.n, block):
        geo = FieldGeometry(array, roi.pts_m[lo:lo + block], speed_of_sound=speed_of_sound, gain=element_gain(n))
        for k, f in enumerate(freqs):
            p = geo.field(float(f), emission)
            out[k, lo:lo + block] = p.real ** 2 + p.imag ** 2
    return out


def _db(x):
    return 10 * np.log10(x)


def _run_width(row_levels: np.ndarray, valid: np.ndarray, i0: int, thr: float, step: float):
    """Width (deg) of the contiguous run >= thr around i0, linear crossings; censored if it reaches the edge."""
    n = len(row_levels)

    def edge(direction):
        i = i0
        while 0 <= i + direction < n and valid[i + direction] and row_levels[i + direction] >= thr:
            i += direction
        j = i + direction
        if not (0 <= j < n) or not valid[j]:
            return i * step, True
        lo, hi = row_levels[i], row_levels[j]
        frac = (lo - thr) / (lo - hi) if lo != hi else 0.0
        return (i + direction * frac) * step, False

    left, cl = edge(-1)
    right, cr = edge(+1)
    return right - left, cl or cr


def map_metrics(roi: FocusRoi, ec_band: np.ndarray, R_mm: float) -> dict:
    """Shape of the focal region on the middle shell (reference resolution)."""
    nb, na = roi.mid_grid.shape
    valid = roi.mid_grid >= 0
    e_focus = ec_band_focus = None
    idx = np.where(valid, roi.mid_grid, 0)
    e_map = np.where(valid, ec_band[idx - 1], np.nan)  # ec_band excludes point 0
    e0 = ec_band_focus_val = e_map[roi.ib0, roi.ia0]
    lmap = _db(e_map / e0)
    step = roi.spec.step_deg
    out = dict(level_map_db=lmap)
    cosb = np.cos(np.radians(roi.b_deg))[:, None] * np.ones((1, na))
    for thr, tag in ((3.0, "3"), (6.0, "6")):
        above = valid & (np.nan_to_num(lmap, nan=-1e9) >= -thr)
        comp = connected_component(above, (roi.ib0, roi.ia0))
        padded = np.pad(valid & True, 1, constant_values=False)
        interior = padded[1:-1, 1:-1] & padded[:-2, 1:-1] & padded[2:, 1:-1] & padded[1:-1, :-2] & padded[1:-1, 2:]
        censored = bool(np.any(comp & ~interior))
        out[f"area_{tag}_deg2"] = float((cosb * step * step)[comp].sum())
        out[f"area_{tag}_censored"] = censored
        if tag == "6":
            comp6 = comp
    # -3 dB widths along the local horizontal and vertical through the focus (mm = angle in rad x R)
    wa, ca = _run_width(lmap[roi.ib0, :], valid[roi.ib0, :], roi.ia0, -3.0, step)
    wb, cb = _run_width(lmap[:, roi.ia0], valid[:, roi.ia0], roi.ib0, -3.0, step)
    out.update(width_h_mm=float(math.radians(wa) * R_mm), width_h_censored=bool(ca),
               width_v_mm=float(math.radians(wb) * R_mm), width_v_censored=bool(cb))
    outside_main = valid & ~dilate(comp6, 1)
    if outside_main.any():
        out["psl_lobe"] = float(-np.nanmax(np.where(outside_main, lmap, -np.inf)))
        mx = local_maxima(np.where(valid, lmap, -np.inf)) & outside_main
        order = np.argsort(-np.where(mx, lmap, -np.inf).ravel())[:3]
        lobes = []
        for o in order:
            ib, ia = divmod(int(o), na)
            if mx[ib, ia]:
                lobes.append(dict(az_offset_deg=float(roi.a_deg[ia]), el_offset_deg=float(roi.b_deg[ib]),
                                  level_db=float(lmap[ib, ia])))
        out["top_lobes"] = lobes
    else:
        out["psl_lobe"] = None
        out["top_lobes"] = []
    cap_out = valid & (np.hypot(*np.meshgrid(roi.a_deg, roi.b_deg)) > EXCLUSION_DEG)
    out["pedestal_median_db"] = float(np.nanmedian(lmap[cap_out])) if cap_out.any() else None
    return out


def depth_metric(roi: FocusRoi, ec_band: np.ndarray, spec: RoiSpec, R_mm: float) -> dict:
    """On-axis compensated level across shells; -3 dB depth in mm (censored at the shell limits)."""
    levels = {}
    for s_idx, s in enumerate(spec.shells):
        sel = (roi.shell[1:] == s_idx) & (roi.cell[1:, 0] == roi.ib0) & (roi.cell[1:, 1] == roi.ia0)
        if sel.any():
            levels[s] = float(ec_band[sel][0])
    e0 = levels[1.0]
    shells = sorted(levels)
    ok = [(_db(levels[s] / e0) >= -3.0) for s in shells]
    i = shells.index(1.0)
    lo = hi = i
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(shells) - 1 and ok[hi + 1]:
        hi += 1
    return dict(depth_mm=float((shells[hi] - shells[lo]) * R_mm), depth_censored=bool(lo == 0 or hi == len(shells) - 1),
                axis_levels_db={f"{s:g}": float(_db(levels[s] / e0)) for s in shells})


def evaluate_reference(array, foci=ALL_FOCI, bands=BANDS, spec: RoiSpec = REFERENCE_ROI, df: float = 50.0,
                       speed_of_sound: float = SPEED_OF_SOUND, keep_maps: bool = False) -> dict:
    """Metrics for every (focus, band) at reference resolution. Returns {"cells": {(focus, band): metrics}, ...}."""
    n = len(array.elements)
    t0 = time.perf_counter()
    cells = {}
    for focus in foci:
        roi = FocusRoi(focus, spec)
        for band, (lo, hi) in bands.items():
            freqs, w = Spectrum(lo, hi, df=df).lines()
            E = reference_energies(array, focus, roi, freqs, speed_of_sound)
            m = core_metrics(E, w, roi)
            ec_band = m["ec"]
            m.update(map_metrics(roi, ec_band, focus.distance_mm))
            m.update(depth_metric(roi, ec_band, spec, focus.distance_mm))
            coherent = element_gain(n) * sum(1.0 / (np.linalg.norm(focus.position_mm - np.append(np.array(e.position[:2]) * 1000, 0)) / 1000)
                                             for e in array.elements)
            m["level_focus_db"] = float(_db(m["e_focus"] / A_REF ** 2))
            m["level_coherent_db"] = float(20 * math.log10(coherent / A_REF))
            m["level_raw_db"] = m["level_focus_db"] - 10 * math.log10(n / N_REF)
            m["F10_p10_lines"] = float(np.percentile(m["F10_lines"], 10))
            m["psl_p10_lines"] = float(np.percentile(m["psl_lines"], 10))
            m["n_lines"] = len(freqs)
            if not keep_maps:
                m.pop("level_map_db", None)
            for k in ("ec", "F10_lines", "psl_lines"):
                m.pop(k, None)
            cells[(focus.name, band)] = m
    return dict(cells=cells, seconds=time.perf_counter() - t0, n=n)


def summarize(result: dict, foci_names: list[str], band: str, weights: Weights = Weights()) -> dict:
    """Mean / worst-case aggregates over a set of foci for one band."""
    ms = [result["cells"][(f, band)] for f in foci_names]
    F10 = np.array([m["F10"] for m in ms])
    psl = np.array([m["psl_cap"] for m in ms])
    p10 = np.array([m["F10_p10_lines"] for m in ms])
    comp = objective_components([dict(F10=m["F10"], psl_cap=m["psl_cap"], F10_lines=np.array([m["F10_p10_lines"]]))
                                 for m in ms], weights)
    return dict(F10_mean=float(F10.mean()), F10_worst=float(F10.min()), psl_mean=float(psl.mean()),
                psl_worst=float(psl.min()), F10_p10_lines_mean=float(p10.mean()), F10_p10_lines_worst=float(p10.min()),
                level_mean=float(np.mean([m["level_focus_db"] for m in ms])), J_ref=comp["J"])


# ---- Pareto -------------------------------------------------------------------------------------------------------

def pareto_mask(points: np.ndarray, maximize: list[bool]) -> np.ndarray:
    """True for non-dominated rows. `points` is (M, K); maximize[k] says whether larger is better."""
    pts = np.asarray(points, dtype=float)
    sign = np.where(np.array(maximize), 1.0, -1.0)
    s = pts * sign
    keep = np.ones(len(s), dtype=bool)
    for i in range(len(s)):
        dominated = np.all(s >= s[i], axis=1) & np.any(s > s[i], axis=1)
        keep[i] = not dominated.any()
    return keep
