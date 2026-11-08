"""Three-dimensional structure of the broadband energy field: orthogonal cuts through the focus.

Physics is not duplicated: every cut (and the optional volume) is evaluated with
`broadband.band_energy` on a `FieldGeometry` built for the cut's points, i.e. the energy
`sum_k w_k |P_k|^2` of the complex field with real focusing delays (see broadband.py, field.py).

Coordinates are those of the project: millimeters in the CLI, meters inside, array frame, no
rotation. Cuts: XY (z = focus z), XZ (y = focus y), YZ (x = focus x); in each cut the first named
axis is horizontal and the second vertical. All cuts pass through the focus.

Region of interest and indicators (all levels dB re one speaker at 1 m, energy):
  * Spreading compensation: Lc = L + 20 log10(r / R_f), r = distance from the aperture center,
    R_f = distance from the aperture center to the focus. The ROI is the shell 0.75 R_f <= r <= 1.25 R_f,
    which keeps the near field (dominated by 1/r) out of every statistic.
  * Focal region at -X dB: 4/6-connected region of the ROI containing the focus where Lc >= L_f - X.
    In a cut it is a *cut area* (cm2); in a volume it is a *volume* (liters). The two are never mixed.
  * Width along an axis: length of the run of ROI cells with Lc >= L_f - X on the axis line through
    the focus, flagged `limited` when the run ends at the window edge or the ROI shell instead of
    falling below the threshold (the true width is then larger).
  * Secondary maximum: highest local maximum of Lc in the ROI outside the (dilated) -3 dB focal
    region; contrast = L_f - that maximum.
  * Energy shares use E_c = E (r / R_f)^2, the energy compensated for spherical spreading, as a
    fraction of the ROI total: the share inside the focal region and the share outside it.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

from .array_config import ArrayConfig, load_array
from .broadband import (DISTRIBUTIONS, Spectrum, band_energy, connected_component, effective_level_db)
from .campaign import ARRAY_DIR, BAND, Params
from .field import FieldGeometry, difference_db
from .focus import SPEED_OF_SOUND, focus_delays

PLANES = {"XY": ("x", "y", "z"), "XZ": ("x", "z", "y"), "YZ": ("y", "z", "x")}  # (horizontal, vertical, normal)
AXIS = {"x": 0, "y": 1, "z": 2}
THRESHOLDS_DB = (3.0, 6.0)


# ---- windows and planes ---------------------------------------------------------------------------------------

def default_windows(focus_mm, front_axis: str = "y", half_width: float = 1000.0,
                    front_range=(100.0, 2000.0)) -> dict[str, tuple[float, float]]:
    """Axis ranges (mm): `front_range` along the beam axis, focus +/- half_width (snapped to 5 mm) elsewhere."""
    if front_axis not in ("y", "z"):
        raise ValueError("front axis must be 'y' or 'z'")
    out = {}
    for axis, i in AXIS.items():
        if axis == front_axis:
            out[axis] = (float(front_range[0]), float(front_range[1]))
        else:
            out[axis] = (round((focus_mm[i] - half_width) / 5) * 5.0, round((focus_mm[i] + half_width) / 5) * 5.0)
    return out


def axis_values(lo: float, hi: float, step: float, label: str = "axis") -> np.ndarray:
    if not (math.isfinite(lo) and math.isfinite(hi) and math.isfinite(step)) or step <= 0 or hi < lo:
        raise ValueError(f"{label}: invalid range or step")
    n = round((hi - lo) / step) + 1
    if abs((hi - lo) / step - (n - 1)) > 1e-9:
        raise ValueError(f"{label}: range {lo}..{hi} is not a multiple of the step {step}")
    return np.linspace(lo, hi, n) if n > 1 else np.array([float(lo)])


def plane_points(plane: str, focus_mm, windows, step_mm: float):
    """(u_mm, v_mm, points_m with shape (nv, nu, 3)) of the cut through the focus."""
    if plane not in PLANES:
        raise ValueError(f"plane must be one of {tuple(PLANES)}")
    ua, va, na = PLANES[plane]
    u = axis_values(*windows[ua], step_mm, ua)
    v = axis_values(*windows[va], step_mm, va)
    uu, vv = np.meshgrid(u, v)
    pts = np.empty(uu.shape + (3,))
    pts[..., AXIS[ua]] = uu
    pts[..., AXIS[va]] = vv
    pts[..., AXIS[na]] = focus_mm[AXIS[na]]
    return u, v, pts / 1000.0


def volume_points(focus_mm, windows, step_mm: float):
    xs = [axis_values(*windows[a], step_mm, a) for a in "xyz"]
    gx, gy, gz = np.meshgrid(*xs, indexing="ij")
    return xs, np.stack([gx, gy, gz], axis=-1) / 1000.0


def energies(array: ArrayConfig, focus_m, spectrum: Spectrum, points_m: np.ndarray, params: Params):
    """(simultaneous, focused) energy sum_k w_k |P_k|^2 at the points."""
    geo = FieldGeometry(array, points_m, speed_of_sound=params.speed_of_sound, gain=params.gain)
    freqs, weights = spectrum.lines()
    emission = {d.speaker: d.delay for d in focus_delays(array, focus_m, params.speed_of_sound)}
    return band_energy(geo, freqs, weights), band_energy(geo, freqs, weights, emission)


# ---- generic image helpers ------------------------------------------------------------------------------------------

def local_maxima(a: np.ndarray) -> np.ndarray:
    """Cells strictly above their 'previous' and at least equal to their 'next' 8-neighbors."""
    padded = np.pad(a, 1, constant_values=-np.inf)
    ny, nx = a.shape
    prev = [(-1, -1), (-1, 0), (-1, 1), (0, -1)]
    nxt = [(0, 1), (1, -1), (1, 0), (1, 1)]
    out = np.ones(a.shape, dtype=bool)
    for dy, dx in prev:
        out &= a > padded[1 + dy:1 + dy + ny, 1 + dx:1 + dx + nx]
    for dy, dx in nxt:
        out &= a >= padded[1 + dy:1 + dy + ny, 1 + dx:1 + dx + nx]
    return out


def dilate(mask: np.ndarray, cells: int) -> np.ndarray:
    out = mask.copy()
    for _ in range(cells):
        padded = np.pad(out, 1)
        nxt = out.copy()
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                nxt |= padded[1 + dy:1 + dy + out.shape[0], 1 + dx:1 + dx + out.shape[1]]
        out = nxt
    return out


def component_nd(mask: np.ndarray, seed: tuple) -> np.ndarray:
    """Face-connected component of an N-d boolean array containing `seed` (numpy frontier search)."""
    seen = np.zeros(mask.shape, dtype=bool)
    if not mask[seed]:
        return seen
    seen[seed] = True
    frontier = np.array(seed).reshape(mask.ndim, 1)
    shape = np.array(mask.shape).reshape(-1, 1)
    while frontier.shape[1]:
        parts = []
        for ax in range(mask.ndim):
            for step in (-1, 1):
                nb = frontier.copy()
                nb[ax] += step
                parts.append(nb)
        cand = np.concatenate(parts, axis=1)
        cand = cand[:, np.all((cand >= 0) & (cand < shape), axis=0)]
        keep = mask[tuple(cand)] & ~seen[tuple(cand)]
        cand = cand[:, keep]
        seen[tuple(cand)] = True
        frontier = np.unique(cand, axis=1)
    return seen


def _run(above: np.ndarray, allowed: np.ndarray, i0: int) -> tuple[int, int, bool, bool]:
    """Contiguous run of `above & allowed` around i0 -> (lo, hi, lo_limited, hi_limited)."""
    ok = above & allowed
    lo = hi = i0
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(ok) - 1 and ok[hi + 1]:
        hi += 1
    lo_lim = lo == 0 or (not allowed[lo - 1])
    hi_lim = hi == len(ok) - 1 or (not allowed[hi + 1])
    return lo, hi, bool(lo_lim), bool(hi_lim)


# ---- cut analysis ----------------------------------------------------------------------------------------------------------

def analyze_cut(plane: str, u_mm, v_mm, e_sim, e_foc, focus_mm, center_mm, l_f: float,
                params: Params, step_mm: float) -> dict:
    ua, va, na = PLANES[plane]
    uu, vv = np.meshgrid(u_mm, v_mm)
    coords = [None, None, None]
    coords[AXIS[ua]], coords[AXIS[va]] = uu, vv
    coords[AXIS[na]] = np.full_like(uu, focus_mm[AXIS[na]])
    r = np.sqrt(sum((coords[i] - center_mm[i]) ** 2 for i in range(3)))
    r_f = math.dist(focus_mm, center_mm)
    comp = 20 * np.log10(r / r_f)
    roi = (r >= (1 - BAND) * r_f) & (r <= (1 + BAND) * r_f)
    lev_f = effective_level_db(e_foc, params.gain, params.floor_db) + comp
    lev_s = effective_level_db(e_sim, params.gain, params.floor_db) + comp
    ec_f, ec_s = e_foc * (r / r_f) ** 2, e_sim * (r / r_f) ** 2
    out: dict = dict(normal_axis=na, normal_mm=float(focus_mm[AXIS[na]]), cells=int(roi.sum()),
                     step_mm=step_mm, roi_cells=int(roi.sum()))
    if not roi.any():
        return out
    fu, fv = focus_mm[AXIS[ua]], focus_mm[AXIS[va]]
    iu = int(np.argmin(np.abs(u_mm - fu)))
    iv = int(np.argmin(np.abs(v_mm - fv)))
    out["focus_on_grid"] = bool(abs(u_mm[iu] - fu) <= step_mm / 2 and abs(v_mm[iv] - fv) <= step_mm / 2)
    if not (out["focus_on_grid"] and roi[iv, iu]):
        return out
    # seed: best ROI cell in the 5 x 5 neighborhood of the focus (the focus need not be a grid node)
    sl_v = slice(max(iv - 2, 0), iv + 3)
    sl_u = slice(max(iu - 2, 0), iu + 3)
    nb = np.where(roi[sl_v, sl_u], lev_f[sl_v, sl_u], -np.inf)
    jv, ju = np.unravel_index(int(np.argmax(nb)), nb.shape)
    sv, su = sl_v.start + jv, sl_u.start + ju
    ny, nx = roi.shape
    comp3 = np.zeros_like(roi)
    for t in THRESHOLDS_DB:
        tag = f"{t:g}"
        above = lev_f >= l_f - t
        comp_mask = connected_component(roi & above, (sv, su))
        if comp_mask.any():
            rows, cols = np.where(comp_mask)
            padded = np.pad(roi, 1, constant_values=False)
            interior = padded[1:-1, 1:-1] & padded[:-2, 1:-1] & padded[2:, 1:-1] & padded[1:-1, :-2] & padded[1:-1, 2:]
            touches_window = bool(rows.min() == 0 or rows.max() == ny - 1 or cols.min() == 0 or cols.max() == nx - 1)
            touches_shell = bool(np.any(comp_mask & ~interior))
            out.update({
                f"area_{tag}_cm2": float(comp_mask.sum() * step_mm * step_mm / 100.0),
                f"u_extent_{tag}_mm": float((cols.max() - cols.min() + 1) * step_mm),
                f"v_extent_{tag}_mm": float((rows.max() - rows.min() + 1) * step_mm),
                f"truncated_{tag}": touches_window or touches_shell,
                f"energy_share_in_{tag}_foc": float(ec_f[comp_mask].sum() / ec_f[roi].sum()),
                f"energy_share_in_{tag}_sim": float(ec_s[comp_mask].sum() / ec_s[roi].sum()),
            })
            if t == 3.0:
                comp3 = comp_mask
        else:
            out.update({f"area_{tag}_cm2": None})
        # widths along the two axis lines through the focus
        lo, hi, lo_l, hi_l = _run(above[sv, :], roi[sv, :], su)
        out.update({f"width_u_{tag}_mm": float((hi - lo + 1) * step_mm), f"width_u_{tag}_limited": lo_l or hi_l})
        lo, hi, lo_l, hi_l = _run(above[:, su], roi[:, su], sv)
        out.update({f"width_v_{tag}_mm": float((hi - lo + 1) * step_mm), f"width_v_{tag}_limited": lo_l or hi_l})
        out[f"roi_pct_above_m{tag}_foc"] = float(100 * np.mean(lev_f[roi] >= l_f - t))
        out[f"roi_pct_above_m{tag}_sim"] = float(100 * np.mean(lev_s[roi] >= l_f - t))
    out["roi_pct_above_m10_foc"] = float(100 * np.mean(lev_f[roi] >= l_f - 10))
    out["roi_pct_above_m10_sim"] = float(100 * np.mean(lev_s[roi] >= l_f - 10))
    # secondary maxima: local maxima of the compensated level outside the dilated -3 dB region
    cand = local_maxima(np.where(np.isfinite(lev_f), lev_f, -np.inf)) & roi & ~dilate(comp3, 2)
    out.update(secondary_db=None, secondary_u_mm=None, secondary_v_mm=None, contrast_db=None)
    if cand.any():
        masked = np.where(cand, lev_f, -np.inf)
        jv, ju = np.unravel_index(int(np.argmax(masked)), masked.shape)
        out.update(secondary_db=float(lev_f[jv, ju]), secondary_u_mm=float(u_mm[ju]), secondary_v_mm=float(v_mm[jv]),
                   contrast_db=float(l_f - lev_f[jv, ju]))
    return out


# ---- scenario ---------------------------------------------------------------------------------------------------------------

@dataclass(frozen=True)
class CutScenario:
    series: str
    key: str
    label: str
    array: str
    focus_mm: tuple[float, float, float]
    spectrum: Spectrum
    front_axis: str = "y"


@dataclass
class CutPlane:
    name: str
    u_mm: np.ndarray
    v_mm: np.ndarray
    e_sim: np.ndarray
    e_foc: np.ndarray
    lev_sim: np.ndarray
    lev_foc: np.ndarray
    diff: np.ndarray
    normal_mm: float


@dataclass
class CutResult:
    scenario: CutScenario
    array: ArrayConfig
    planes: dict[str, CutPlane]
    indicators: dict
    windows: dict
    step_mm: float
    volume: dict | None = None


def aperture_center_mm(array: ArrayConfig) -> np.ndarray:
    return np.mean([e.position for e in array.elements], axis=0) * 1000.0


def run_cuts(sc: CutScenario, params: Params, *, step_mm: float = 5.0, windows=None, volume: bool = False,
             volume_step_mm: float = 25.0, array_dir: Path = ARRAY_DIR) -> CutResult:
    array = load_array(array_dir / sc.array)
    focus_m = tuple(v / 1000 for v in sc.focus_mm)
    windows = windows or default_windows(sc.focus_mm, sc.front_axis)
    center = aperture_center_mm(array)
    r_f = math.dist(sc.focus_mm, center)
    # exact levels at the focus (not grid dependent)
    e_s0, e_f0 = energies(array, focus_m, sc.spectrum, np.array([focus_m]), params)
    l_s = float(effective_level_db(e_s0, params.gain, params.floor_db)[0])
    l_f = float(effective_level_db(e_f0, params.gain, params.floor_db)[0])
    freqs, weights = sc.spectrum.lines()
    ind: dict = dict(scenario=f"{sc.series}/{sc.key}", array=sc.array, spectrum=sc.spectrum.label,
                     f_min_hz=sc.spectrum.f_min, f_max_hz=sc.spectrum.f_max, n_lines=len(freqs),
                     weights_sum=float(weights.sum()), focus_mm=list(sc.focus_mm), front_axis=sc.front_axis,
                     aperture_center_mm=[float(c) for c in center], focus_range_mm=r_f,
                     steering_angle_deg=math.degrees(math.acos(
                         max(-1.0, min(1.0, (sc.focus_mm[AXIS[sc.front_axis]] - center[AXIS[sc.front_axis]]) / r_f)))),
                     level_sim_focus_db=l_s, level_foc_focus_db=l_f, focusing_gain_db=l_f - l_s, windows=windows)
    planes: dict[str, CutPlane] = {}
    for name in PLANES:
        u, v, pts = plane_points(name, sc.focus_mm, windows, step_mm)
        e_s, e_f = energies(array, focus_m, sc.spectrum, pts, params)
        planes[name] = CutPlane(name, u, v, e_s, e_f,
                                effective_level_db(e_s, params.gain, params.floor_db),
                                effective_level_db(e_f, params.gain, params.floor_db),
                                difference_db(np.sqrt(e_f), np.sqrt(e_s), params.gain, params.floor_db),
                                float(sc.focus_mm[AXIS[PLANES[name][2]]]))
        ind[name] = analyze_cut(name, u, v, e_s, e_f, sc.focus_mm, center, l_f, params, step_mm)
    # widths per axis (the same axis line appears in two cuts; they agree by construction)
    ind["axis_width_mm"] = {t: {"x": ind["XY"].get(f"width_u_{t:g}_mm"), "y": ind["XY"].get(f"width_v_{t:g}_mm"),
                                "z": ind["XZ"].get(f"width_v_{t:g}_mm")} for t in THRESHOLDS_DB}
    ind["axis_width_limited"] = {t: {"x": ind["XY"].get(f"width_u_{t:g}_limited"), "y": ind["XY"].get(f"width_v_{t:g}_limited"),
                                     "z": ind["XZ"].get(f"width_v_{t:g}_limited")} for t in THRESHOLDS_DB}
    ind["axis_width_mm"] = {f"{t:g}": v for t, v in ind["axis_width_mm"].items()}
    ind["axis_width_limited"] = {f"{t:g}": v for t, v in ind["axis_width_limited"].items()}
    result = CutResult(sc, array, planes, ind, windows, step_mm)
    if volume:
        result.volume = run_volume(sc, array, params, windows, volume_step_mm, l_f, center)
        ind["volume"] = result.volume
    return result


def run_volume(sc: CutScenario, array: ArrayConfig, params: Params, windows, step_mm: float, l_f: float,
               center_mm) -> dict:
    """Focal volume (liters) and energy distribution on a coarse 3-D grid, same ROI shell as the cuts."""
    focus_m = tuple(v / 1000 for v in sc.focus_mm)
    axes, pts = volume_points(sc.focus_mm, windows, step_mm)
    e_s, e_f = energies(array, focus_m, sc.spectrum, pts, params)
    gx, gy, gz = np.meshgrid(*axes, indexing="ij")
    r = np.sqrt((gx - center_mm[0]) ** 2 + (gy - center_mm[1]) ** 2 + (gz - center_mm[2]) ** 2)
    r_f = math.dist(sc.focus_mm, center_mm)
    shell = (r >= (1 - BAND) * r_f) & (r <= (1 + BAND) * r_f)
    lev = effective_level_db(e_f, params.gain, params.floor_db) + 20 * np.log10(r / r_f)
    lev_s = effective_level_db(e_s, params.gain, params.floor_db) + 20 * np.log10(r / r_f)
    ec_f, ec_s = e_f * (r / r_f) ** 2, e_s * (r / r_f) ** 2
    seed = tuple(int(np.argmin(np.abs(a - c))) for a, c in zip(axes, sc.focus_mm))
    out: dict = dict(step_mm=step_mm, shape=list(shell.shape), shell_cells=int(shell.sum()),
                     shell_volume_liters=float(shell.sum() * step_mm ** 3 / 1e6))
    if not shell[seed]:
        return out
    for t in THRESHOLDS_DB:
        above = shell & (lev >= l_f - t)
        comp = component_nd(above, seed)
        idx = np.where(comp)
        tag = f"{t:g}"
        if comp.any():
            touches_box = any(i.min() == 0 or i.max() == s - 1 for i, s in zip(idx, comp.shape))
            ext = [float((i.max() - i.min() + 1) * step_mm) for i in idx]
            padded = np.pad(shell, 1, constant_values=False)
            interior = shell.copy()
            for ax in range(3):
                for s in (-1, 1):
                    interior &= np.roll(padded, s, axis=ax)[1:-1, 1:-1, 1:-1]
            out.update({f"volume_{tag}_liters": float(comp.sum() * step_mm ** 3 / 1e6),
                        f"extent_{tag}_mm": {"x": ext[0], "y": ext[1], "z": ext[2]},
                        f"truncated_{tag}": bool(touches_box or np.any(comp & ~interior)),
                        f"energy_share_in_{tag}_foc": float(ec_f[comp].sum() / ec_f[shell].sum()),
                        f"energy_share_in_{tag}_sim": float(ec_s[comp].sum() / ec_s[shell].sum())})
        out[f"shell_pct_above_m{tag}_foc"] = float(100 * np.mean(lev[shell] >= l_f - t))
        out[f"shell_pct_above_m{tag}_sim"] = float(100 * np.mean(lev_s[shell] >= l_f - t))
    return out


# ---- scenarios of the campaign ---------------------------------------------------------------------------------------------------------------

TONE = Spectrum(4000, 4000)
BAND_24 = Spectrum(2000, 4000)
BAND_18 = Spectrum(1000, 8000)
BAND_1535 = Spectrum(1500, 3500)
GEOMETRIES = {"L": ("line-8-60mm.json", "y"), "P": ("staggered-3-2-3.json", "z"),
              "W": ("staggered-3-2-3-wide.json", "z")}  # W = 3-2-3 scaled 3.5x: 420 x 420 mm, as wide as L
SPECTRA = {"tone": ("4 kHz tone", TONE), "b24": ("2-4 kHz", BAND_24), "b18": ("1-8 kHz", BAND_18)}
LATERAL_DEG = 20.0
LATERAL_RANGE_MM = 800.0


def broadside_focus(geometry: str) -> tuple[float, float, float]:
    """Focus 800 mm in front of the aperture center along the beam axis (y for L, z for P and W)."""
    array = load_array(ARRAY_DIR / GEOMETRIES[geometry][0])
    c = aperture_center_mm(array)
    return (float(c[0]), 800.0, 0.0) if GEOMETRIES[geometry][1] == "y" else (float(c[0]), float(c[1]), 800.0)


def lateral_focus(geometry: str) -> tuple[float, float, float]:
    """Same range from the aperture center, displaced LATERAL_DEG toward +x."""
    a = math.radians(LATERAL_DEG)
    array = load_array(ARRAY_DIR / GEOMETRIES[geometry][0])
    c = aperture_center_mm(array)
    dx, dn = LATERAL_RANGE_MM * math.sin(a), LATERAL_RANGE_MM * math.cos(a)
    return (float(c[0] + dx), dn, 0.0) if GEOMETRIES[geometry][1] == "y" else (float(c[0] + dx), float(c[1]), dn)


def cut_scenarios(series: str) -> list[CutScenario]:
    if series in ("A", "A0", "B", "W"):
        g = {"A": "L", "A0": "L", "B": "P", "W": "W"}[series]
        arr, front = GEOMETRIES[g]
        focus = (350.0, 800.0, 0.0) if series == "A" else broadside_focus(g)  # A: the focus used since Mission 03
        return [CutScenario(series, key, f"{g}: {label}", arr, focus, sp, front) for key, (label, sp) in SPECTRA.items()]
    if series == "C":
        out = []
        for g in ("L", "P", "W"):
            arr, front = GEOMETRIES[g]
            for key in ("tone", "b24"):
                label, sp = SPECTRA[key]
                out.append(CutScenario("C", f"{g}-{key}", f"{g}: {label}, {LATERAL_DEG:g} deg", arr, lateral_focus(g), sp, front))
        return out
    if series == "D":
        arr, front = GEOMETRIES["L"]
        focus = (800.0, 800.0, 0.0)
        return [CutScenario("D", key, f"L: {label}, focus x = 800", arr, focus, sp, front)
                for key, label, sp in (("b1535", "1.5-3.5 kHz", BAND_1535), ("b24", "2-4 kHz", BAND_24),
                                       ("b18", "1-8 kHz", BAND_18))]
    raise ValueError(f"unknown series {series!r}")


# ---- plots ---------------------------------------------------------------------------------------------------------------------------

def _speaker_marks(ax, array: ArrayConfig, plane: str, focus_mm):
    ua, va, na = PLANES[plane]
    inp, proj = [], []
    for e in array.elements:
        p = [c * 1000 for c in e.position]
        (inp if abs(p[AXIS[na]] - focus_mm[AXIS[na]]) < 1e-6 else proj).append((p[AXIS[ua]], p[AXIS[va]]))
    if proj:
        ax.plot(*zip(*proj), "v", markerfacecolor="none", markeredgecolor="dimgray", markersize=5)
    if inp:
        ax.plot(*zip(*inp), "s", color="white", markeredgecolor="black", markersize=5)
    return inp + proj


def decorate_cut(ax, array, plane, focus_mm, u_mm, v_mm, title=None, *, fontsize=8) -> None:
    ua, va, na = PLANES[plane]
    pts = _speaker_marks(ax, array, plane, focus_mm)
    ax.plot([focus_mm[AXIS[ua]]], [focus_mm[AXIS[va]]], "*", color="yellow", markeredgecolor="black", markersize=11)
    us = list(u_mm) + [p[0] for p in pts]
    vs = list(v_mm) + [p[1] for p in pts]
    pad = 0.03 * max(max(us) - min(us), max(vs) - min(vs))
    ax.set_xlim(min(us) - pad, max(us) + pad)
    ax.set_ylim(min(vs) - pad, max(vs) + pad)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{ua} (mm)", fontsize=fontsize)
    ax.set_ylabel(f"{va} (mm)", fontsize=fontsize)
    ax.tick_params(labelsize=6)
    if title:
        ax.set_title(title, fontsize=fontsize)


def _extent(u, v):
    du = (u[1] - u[0]) / 2 if len(u) > 1 else 0.5
    dv = (v[1] - v[0]) / 2 if len(v) > 1 else 0.5
    return u[0] - du, u[-1] + du, v[0] - dv, v[-1] + dv


NOTE = ("star = focus; white squares = speakers in the plane; gray triangles = speakers projected onto the plane. "
        "Levels outside the scale saturate (arrows on the colorbar); data are not clipped.")


def build_cut_plate(title: str, result: CutResult, params: Params):
    from matplotlib.figure import Figure

    from .fieldplot import LEVEL_CMAP, draw_difference, draw_level

    fig = Figure(figsize=(13.5, 12))
    gs = fig.add_gridspec(3, 4, width_ratios=[1, 1, 1, 0.05], wspace=0.25, hspace=0.18, top=0.92, bottom=0.05)
    vmax, vmin = params.vmax, params.vmax - params.range_db
    level_ims, diff_ims = [], []
    names = list(PLANES)
    rows = ("Simultaneous (energy)", "Focused (energy)", "Focused - simultaneous (dB)")
    for j, name in enumerate(names):
        cut = result.planes[name]
        ua, va, na = PLANES[name]
        for i in range(3):
            ax = fig.add_subplot(gs[i, j])
            if i == 0:
                im = draw_level(ax, cut.lev_sim, cut.u_mm, cut.v_mm, vmin, vmax)
                level_ims.append(im)
            elif i == 1:
                im = draw_level(ax, cut.lev_foc, cut.u_mm, cut.v_mm, vmin, vmax)
                level_ims.append(im)
            else:
                im = draw_difference(ax, cut.diff, cut.u_mm, cut.v_mm, params.diff_limit)
                diff_ims.append(im)
            decorate_cut(ax, result.array, name, result.scenario.focus_mm, cut.u_mm, cut.v_mm,
                         f"{name} cut, {na} = {cut.normal_mm:g} mm" if i == 0 else None)
            if j == 0:
                ax.annotate(rows[i], xy=(-0.28, 0.5), xycoords="axes fraction", rotation=90, va="center", fontsize=9)
    cb1 = fig.colorbar(level_ims[0], cax=fig.add_subplot(gs[0:2, 3]), extend="both")
    cb1.set_label(f"relative energy, dB re one speaker at 1 m (model, not SPL); scale {vmin:g}..{vmax:g} dB", fontsize=8)
    cb2 = fig.colorbar(diff_ims[0], cax=fig.add_subplot(gs[2, 3]), extend="both")
    cb2.set_label("focused - simultaneous (dB)", fontsize=8)
    fig.suptitle(f"{title}\n{NOTE}", fontsize=9, y=0.985)
    return fig, level_ims, diff_ims


def role_planes(front_axis: str) -> list[tuple[str, str]]:
    """(role, plane) in a fixed order: plane holding the beam axis and x, beam axis and the other axis, transverse."""
    if front_axis == "y":
        return [("beam axis + x", "XY"), ("beam axis + z", "YZ"), ("transverse", "XZ")]
    return [("beam axis + x", "XZ"), ("beam axis + y", "YZ"), ("transverse", "XY")]


def build_comparison_plate(title: str, results: list[CutResult], params: Params, labels: list[str] | None = None):
    """Focused-energy cuts of several scenarios side by side; rows by role of the plane."""
    from matplotlib.figure import Figure

    from .fieldplot import draw_level

    n = len(results)
    fig = Figure(figsize=(4.4 * n + 1.2, 12.5))
    gs = fig.add_gridspec(3, n + 1, width_ratios=[1] * n + [0.05], wspace=0.25, hspace=0.25, top=0.9, bottom=0.05)
    vmax, vmin = params.vmax, params.vmax - params.range_db
    ims = []
    for j, res in enumerate(results):
        for i, (role, name) in enumerate(role_planes(res.scenario.front_axis)):
            ax = fig.add_subplot(gs[i, j])
            cut = res.planes[name]
            ims.append(draw_level(ax, cut.lev_foc, cut.u_mm, cut.v_mm, vmin, vmax))
            head = (labels[j] if labels else res.scenario.label) + "\n" if i == 0 else ""
            decorate_cut(ax, res.array, name, res.scenario.focus_mm, cut.u_mm, cut.v_mm,
                         f"{head}{role}: {name} cut, {PLANES[name][2]} = {cut.normal_mm:g} mm", fontsize=8)
    cb = fig.colorbar(ims[0], cax=fig.add_subplot(gs[:, n]), extend="both")
    cb.set_label(f"focused relative energy, dB re one speaker at 1 m (model, not SPL); scale {vmin:g}..{vmax:g} dB", fontsize=8)
    fig.suptitle(f"{title}\n{NOTE}", fontsize=9, y=0.985)
    return fig, ims


# ---- output ---------------------------------------------------------------------------------------------------------------------------------

def _flatten(ind: dict, prefix: str = "") -> dict:
    out = {}
    for k, v in ind.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten(v, key + "."))
        elif isinstance(v, (list, tuple)):
            out[key] = ";".join(f"{x:.3f}" if isinstance(x, float) else str(x) for x in v)
        else:
            out[key] = v
    return out


def write_result(outdir: Path, r: CutResult, params: Params, individual: bool = True, flat: bool = False) -> Path:
    d = outdir if flat else outdir / r.scenario.series / r.scenario.key
    d.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for name, cut in r.planes.items():
        arrays.update({f"{name}_u_mm": cut.u_mm, f"{name}_v_mm": cut.v_mm, f"{name}_energy_simultaneous": cut.e_sim,
                       f"{name}_energy_focused": cut.e_foc, f"{name}_simultaneous_db": cut.lev_sim,
                       f"{name}_focused_db": cut.lev_foc, f"{name}_difference_db": cut.diff})
    freqs, weights = r.scenario.spectrum.lines()
    np.savez_compressed(d / "data.npz", frequencies_hz=freqs, weights=weights,
                        indicators=json.dumps(r.indicators), **arrays)
    if individual:
        fig, _, _ = build_cut_plate(f"{r.scenario.series}/{r.scenario.key}: {r.scenario.label}, focus "
                                    f"({r.scenario.focus_mm[0]:.1f}, {r.scenario.focus_mm[1]:.1f}, {r.scenario.focus_mm[2]:.1f}) mm",
                                    r, params)
        fig.savefig(d / "cuts.png", dpi=110, bbox_inches="tight")
    return d


def write_tables(outdir: Path, results: list[CutResult]) -> None:
    rows = [_flatten(r.indicators) for r in results]
    keys = list(dict.fromkeys(k for row in rows for k in row if not k.startswith("windows.")))
    with (outdir / "cut-indicators.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (outdir / "cut-indicators.json").write_text(json.dumps([r.indicators for r in results], indent=1) + "\n",
                                                encoding="utf-8")


def run_campaign(series: list[str], outdir: Path, params: Params, df: float, *, step_mm: float = 5.0,
                 volume: bool = False, volume_step_mm: float = 25.0, individual: bool = True, log=print,
                 array_dir: Path = ARRAY_DIR) -> list[CutResult]:
    outdir.mkdir(parents=True, exist_ok=True)
    results: dict[str, CutResult] = {}
    for s in series:
        for sc in cut_scenarios(s):
            sc = replace(sc, spectrum=replace(sc.spectrum, df=df))
            log(f"{sc.series}/{sc.key}: {sc.label}, focus {tuple(round(v, 1) for v in sc.focus_mm)}, "
                f"{len(sc.spectrum.lines()[0])} lines")
            r = run_cuts(sc, params, step_mm=step_mm, volume=volume, volume_step_mm=volume_step_mm, array_dir=array_dir)
            write_result(outdir, r, params, individual)
            results[f"{sc.series}/{sc.key}"] = r
    # comparison plates between geometries (same spectrum) and between bands (series D)
    groups = {
        "broadside-tone": [("A0", "tone"), ("B", "tone"), ("W", "tone")],
        "broadside-b24": [("A0", "b24"), ("B", "b24"), ("W", "b24")],
        "broadside-b18": [("A0", "b18"), ("B", "b18"), ("W", "b18")],
        "lateral-tone": [("C", "L-tone"), ("C", "P-tone"), ("C", "W-tone")],
        "lateral-b24": [("C", "L-b24"), ("C", "P-b24"), ("C", "W-b24")],
        "extreme-bands": [("D", "b1535"), ("D", "b24"), ("D", "b18")],
    }
    for name, members in groups.items():
        rs = [results[f"{s}/{k}"] for s, k in members if f"{s}/{k}" in results]
        if len(rs) == len(members):
            fig, _ = build_comparison_plate(name.replace("-", " "), rs, params)
            fig.savefig(outdir / f"compare-{name}.png", dpi=110, bbox_inches="tight")
    write_tables(outdir, list(results.values()))
    return list(results.values())


# ---- command line --------------------------------------------------------------------------------------------------------------------------------

def _fmt(v, f="{:.2f}"):
    return "n/a" if v is None else f.format(v)


def summary(results: list[CutResult]) -> None:
    print("scenario       L_foc  gain   width x/y/z at -3 dB (mm)        cut areas -3 dB XY/XZ/YZ (cm2)    contrast XY/XZ/YZ (dB)")
    for r in results:
        i = r.indicators
        w = i["axis_width_mm"]["3"]
        lim = i["axis_width_limited"]["3"]
        wt = "/".join(f"{_fmt(w[a], '{:.0f}')}{'+' if lim[a] else ''}" for a in "xyz")
        ar = "/".join(_fmt(i[p].get("area_3_cm2"), "{:.0f}") for p in PLANES)
        ct = "/".join(_fmt(i[p].get("contrast_db"), "{:.1f}") for p in PLANES)
        print(f"{i['scenario']:<13} {i['level_foc_focus_db']:6.2f} {i['focusing_gain_db']:+5.2f}   {wt:<30} {ar:<32} {ct}")


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--df", type=float, default=25.0, help="spectral line spacing, Hz")
    parser.add_argument("--step", type=float, default=5.0, help="cut grid step, mm")
    parser.add_argument("--volume", action="store_true", help="also compute a coarse 3-D volume (off by default)")
    parser.add_argument("--volume-step", type=float, default=25.0, help="volume grid step, mm")
    parser.add_argument("--vmax", type=float, default=20.0, help="top of the level color scale, dB (display only)")
    parser.add_argument("--range-db", type=float, default=40.0)
    parser.add_argument("--diff-limit", type=float, default=10.0)
    parser.add_argument("--speed-of-sound", type=float, default=SPEED_OF_SOUND)
    parser.add_argument("--gain", type=float, default=0.05)


def _params(args) -> Params:
    if not args.range_db > 0 or not args.diff_limit > 0 or not args.step > 0 or not args.volume_step > 0:
        raise ValueError("--range-db, --diff-limit, --step and --volume-step must be positive")
    return Params(speed_of_sound=args.speed_of_sound, gain=args.gain, vmax=args.vmax, range_db=args.range_db,
                  diff_limit=args.diff_limit)


def field_cuts_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab field-cuts",
                                     description="XY, XZ and YZ cuts of the broadband energy field through the focus")
    parser.add_argument("--array", type=Path, required=True)
    parser.add_argument("--target", "--focus", dest="target", type=float, nargs=3, metavar=("X", "Y", "Z"), required=True)
    parser.add_argument("--fmin", type=float, required=True)
    parser.add_argument("--fmax", type=float, required=True)
    parser.add_argument("--nfreq", type=int)
    parser.add_argument("--distribution", choices=DISTRIBUTIONS, default="uniform")
    parser.add_argument("--front-axis", choices=("y", "z"), help="beam axis (default: the one most aligned with the focus)")
    parser.add_argument("--half-width", type=float, default=1000.0, help="transverse half-width of the cuts, mm")
    parser.add_argument("--front-range", type=float, nargs=2, default=(100.0, 2000.0), metavar=("MIN", "MAX"))
    parser.add_argument("--outdir", type=Path, required=True)
    _common(parser)
    args = parser.parse_args(argv)
    try:
        params = _params(args)
        array = load_array(args.array)
        center = aperture_center_mm(array)
        rel = np.array(args.target) - center
        front = args.front_axis or ("y" if abs(rel[1]) >= abs(rel[2]) else "z")
        spectrum = Spectrum(args.fmin, args.fmax, df=args.df, n=args.nfreq, distribution=args.distribution)
        sc = CutScenario("X", "cuts", spectrum.label, args.array.name, tuple(args.target), spectrum, front)
        windows = default_windows(sc.focus_mm, front, args.half_width, args.front_range)
        result = run_cuts(sc, params, step_mm=args.step, windows=windows, volume=args.volume,
                          volume_step_mm=args.volume_step, array_dir=args.array.parent)
        write_result(args.outdir, result, params, flat=True)
    except (OSError, ValueError, KeyError, ImportError) as exc:
        parser.error(str(exc))
    print(f"Wrote cuts.png and data.npz to {args.outdir}")
    summary([result])
    return 0


def field3d_campaign_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab field3d-campaign",
                                     description="3-D campaign: A linear (350,800,0), A0 linear on-axis, B planar 3-2-3, W wide planar, C lateral, D lateral band")
    parser.add_argument("--series", nargs="+", choices=("A", "A0", "B", "W", "C", "D", "all"), default=["all"])
    parser.add_argument("--outdir", type=Path, default=Path("results/field3d"))
    parser.add_argument("--array-dir", type=Path, default=ARRAY_DIR)
    parser.add_argument("--no-individual", action="store_true")
    _common(parser)
    args = parser.parse_args(argv)
    series = ["A", "A0", "B", "W", "C", "D"] if "all" in args.series else list(dict.fromkeys(args.series))
    try:
        params = _params(args)
        results = run_campaign(series, args.outdir, params, args.df, step_mm=args.step, volume=args.volume,
                               volume_step_mm=args.volume_step, individual=not args.no_individual,
                               array_dir=args.array_dir)
    except (OSError, ValueError, KeyError, ImportError) as exc:
        parser.error(str(exc))
    print(f"Wrote cuts, comparison plates and indicators to {args.outdir}")
    summary(results)
    return 0
