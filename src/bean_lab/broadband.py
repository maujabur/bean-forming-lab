"""Broadband (multi-frequency) energy maps of a speaker array.

Physics. Every speaker emits the same stationary signal whose power spectrum is a weighted set of
sinusoidal lines (f_k, w_k), sum w_k = 1, delayed by its real focusing delay e_i (not a phase
shift for one frequency). Lines of different frequency are orthogonal over a long window, so
mean-square pressures add and phasors of different frequencies are never added:

    p_rms^2(x) = (1/2) sum_k w_k |P_k(x)|^2 ,   P_k = complex_field at f_k (field.py, amplitude A)

Each speaker emits the same total power A^2 / 2 whatever the spectrum (sum w_k = 1), so tone and
band scenarios are directly comparable. Levels are 10 log10 of p_rms^2 relative to one speaker
heard at 1 m (A^2 / 2); equivalently 20 log10 of the "effective amplitude"
a = sqrt(sum w_k |P_k|^2), which for a single line is |P| and reproduces the monochromatic maps.

This is the time-averaged energy map, not an instantaneous pressure map, and the average of an
ideal continuous spectrum: a finite realization of noise fluctuates around it (docs/broadband-study.md).

Quadrature. A band [f_min, f_max] is split in bins; each line sits at the middle of its bin and
carries the exact integral of the spectral density over the bin (uniform: bin width; pink, PSD
proportional to 1/f: ln of the edge ratio), normalized to sum 1. f_min == f_max is a tone.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .array_config import ArrayConfig, load_array
from .campaign import (ARC_LIMIT_DEG, ARC_STEP_DEG, ARRAY_DIR, BAND, BAND_ARCS, Params, Region,
                       Scenario, _arc_points, _flat, analyze_lobes, aperture_center, build_arcs,
                       build_plate, uniform_spacing)
from .field import FieldGeometry, difference_db, grid_axis, level_db
from .focus import SPEED_OF_SOUND, focus_delays

DISTRIBUTIONS = ("uniform", "pink")
THRESHOLDS_DB = (3.0, 6.0, 10.0)  # region above (focus level - X) is counted
SHARE_WINDOWS_DEG = (5.0, 10.0)


# ---- spectrum -------------------------------------------------------------------------------------

@dataclass(frozen=True)
class Spectrum:
    f_min: float
    f_max: float
    df: float = 25.0  # target line spacing in Hz (ignored when n is given)
    n: int | None = None  # number of lines; overrides df
    distribution: str = "uniform"

    def __post_init__(self) -> None:
        if not (math.isfinite(self.f_min) and math.isfinite(self.f_max)) or not 0 < self.f_min <= self.f_max:
            raise ValueError("frequencies must satisfy 0 < f_min <= f_max")
        if not (math.isfinite(self.df) and self.df > 0):
            raise ValueError("df must be finite and positive")
        if self.n is not None and (type(self.n) is not int or self.n < 1):
            raise ValueError("number of spectral lines must be a positive integer")
        if self.distribution not in DISTRIBUTIONS:
            raise ValueError(f"distribution must be one of {DISTRIBUTIONS}")

    @property
    def is_tone(self) -> bool:
        return self.f_min == self.f_max

    @property
    def center(self) -> float:
        return self.f_min if self.is_tone else math.sqrt(self.f_min * self.f_max)

    @property
    def label(self) -> str:
        if self.is_tone:
            return f"{self.f_min / 1000:g} kHz tone"
        tag = "" if self.distribution == "uniform" else f" ({self.distribution})"
        lo, hi = self.f_min, self.f_max
        fmt = (lambda f: f"{f / 1000:g} kHz") if lo >= 1000 else (lambda f: f"{f:g} Hz")
        return f"{fmt(lo)} to {fmt(hi)}{tag}"

    def lines(self) -> tuple[np.ndarray, np.ndarray]:
        """Line frequencies (Hz) and weights (sum exactly normalized to 1)."""
        if self.is_tone:
            return np.array([self.f_min]), np.array([1.0])
        n = self.n or max(1, math.ceil((self.f_max - self.f_min) / self.df - 1e-9))
        edges = np.linspace(self.f_min, self.f_max, n + 1)
        if self.distribution == "uniform":
            freqs = (edges[:-1] + edges[1:]) / 2
            weights = np.diff(edges)
        else:  # PSD ~ 1/f: power of a bin = ln(f_hi / f_lo); sample at the geometric mean
            freqs = np.sqrt(edges[:-1] * edges[1:])
            weights = np.diff(np.log(edges))
        return freqs, weights / weights.sum()


def band_energy(geometry: FieldGeometry, freqs: np.ndarray, weights: np.ndarray,
                emission=None) -> np.ndarray:
    """sum_k w_k |P_k|^2 at the geometry's points (amplitude units^2, not halved)."""
    total = np.zeros(geometry.shape)
    for f, w in zip(freqs, weights):
        p = geometry.field(float(f), emission)
        total += w * (p.real ** 2 + p.imag ** 2)
    return total


def effective_level_db(energy: np.ndarray, reference: float, floor_db: float) -> np.ndarray:
    """10 log10(energy / reference^2), clamped like level_db (floor only for display)."""
    return level_db(np.sqrt(energy), reference, floor_db)


# ---- scenarios -------------------------------------------------------------------------------------

@dataclass(frozen=True)
class BroadbandScenario:
    series: str
    key: str
    label: str
    array: str
    focus_mm: tuple[float, float, float]
    spectrum: Spectrum

    def as_scenario(self) -> Scenario:
        return Scenario(self.series, self.key, self.label, self.array, self.spectrum.center, self.focus_mm)


def broadband_scenarios(series: str, df: float = 25.0) -> list[BroadbandScenario]:
    base, focus = "line-8-60mm.json", (350.0, 800.0, 0.0)

    def sp(lo, hi):
        return Spectrum(lo, hi, df=df)

    if series == "B":
        return [BroadbandScenario("B", "B1", "B1 tone 4 kHz", base, focus, sp(4000, 4000)),
                BroadbandScenario("B", "B2", "B2 2-4 kHz", base, focus, sp(2000, 4000)),
                BroadbandScenario("B", "B3", "B3 1-8 kHz", base, focus, sp(1000, 8000)),
                BroadbandScenario("B", "B4", "B4 0.5-8 kHz", base, focus, sp(500, 8000))]
    if series == "G":
        out = []
        for band, (lo, hi) in (("B2", (2000, 4000)), ("B3", (1000, 8000))):
            for d in (30, 60, 120):
                out.append(BroadbandScenario("G", f"d{d:03d}-{band}", f"d = {d} mm, {lo / 1000:g}-{hi / 1000:g} kHz",
                                             f"line-8-{d}mm.json", focus, sp(lo, hi)))
        return out
    if series == "D":
        return [BroadbandScenario("D", f"x{x:g}", f"focus x = {x:g} mm, 1-8 kHz", base, (x, 800.0, 0.0), sp(1000, 8000))
                for x in (210.0, 350.0, 800.0)]
    raise ValueError(f"unknown series {series!r}")


# ---- helpers ------------------------------------------------------------------------------------------

def connected_component(mask: np.ndarray, seed: tuple[int, int]) -> np.ndarray:
    """4-connected component of `mask` containing `seed` (frontier BFS, numpy only)."""
    seen = np.zeros(mask.shape, dtype=bool)
    if not mask[seed]:
        return seen
    seen[seed] = True
    iy, ix = np.array([seed[0]]), np.array([seed[1]])
    ny, nx = mask.shape
    while iy.size:
        cy = np.concatenate([iy + 1, iy - 1, iy, iy])
        cx = np.concatenate([ix, ix, ix + 1, ix - 1])
        ok = (cy >= 0) & (cy < ny) & (cx >= 0) & (cx < nx)
        cy, cx = cy[ok], cx[ok]
        new = mask[cy, cx] & ~seen[cy, cx]
        cy, cx = cy[new], cx[new]
        seen[cy, cx] = True
        keys = np.unique(cy * nx + cx)
        iy, ix = keys // nx, keys % nx
    return seen


def _beamwidth_if_defined(lobes: dict) -> tuple[float | None, str]:
    """-3 dB width only for a closed lobe that stands out from the others by at least 3 dB."""
    if lobes["main_lobe_open"] or lobes["beamwidth_open"]:
        return None, "main lobe not bounded by minima inside the arc"
    if lobes["contrast_db"] is not None and lobes["contrast_db"] < 3.0:
        return None, "no single dominant lobe (highest sidelobe within 3 dB)"
    return lobes["beamwidth_deg"], "defined"


@dataclass
class BroadbandResult:
    scenario: Scenario  # carries series/key/label/focus for plate titles
    spec_scenario: BroadbandScenario
    array: ArrayConfig
    x_mm: np.ndarray
    y_mm: np.ndarray
    energy_sim: np.ndarray
    energy_foc: np.ndarray
    level_sim: np.ndarray
    level_foc: np.ndarray
    diff: np.ndarray
    indicators: dict
    arc_theta: np.ndarray
    arc_sim: np.ndarray
    arc_foc: np.ndarray
    freqs: np.ndarray
    weights: np.ndarray
    per_frequency_contrast: np.ndarray
    per_frequency_share10: np.ndarray
    grating: list = field(default_factory=list)


def map_energies(array, focus, spectrum: Spectrum, x, y, z, params: Params):
    """(energy_sim, energy_foc) on the grid; one frequency at a time, no field is stored."""
    kw = dict(speed_of_sound=params.speed_of_sound, gain=params.gain)
    emission = {d.speaker: d.delay for d in focus_delays(array, focus, params.speed_of_sound)}
    xx, yy = np.meshgrid(x, y)
    geo = FieldGeometry(array, np.stack([xx, yy, np.full_like(xx, z)], axis=-1), **kw)
    freqs, weights = spectrum.lines()
    return band_energy(geo, freqs, weights), band_energy(geo, freqs, weights, emission)


# ---- one scenario -----------------------------------------------------------------------------------------------

def run_broadband(bs: BroadbandScenario, region: Region, params: Params,
                  array_dir: Path = ARRAY_DIR) -> BroadbandResult:
    array = load_array(array_dir / bs.array)
    c, gain = params.speed_of_sound, params.gain
    kw = dict(speed_of_sound=c, gain=gain)
    focus = tuple(v / 1000 for v in bs.focus_mm)
    emission = {d.speaker: d.delay for d in focus_delays(array, focus, c)}
    z = region.z_mm / 1000
    x = grid_axis(region.x_mm[0] / 1000, region.x_mm[1] / 1000, region.nx, "x")
    y = grid_axis(region.y_mm[0] / 1000, region.y_mm[1] / 1000, region.ny, "y")
    freqs, weights = bs.spectrum.lines()

    center = aperture_center(array)
    cx, cy = float(center[0]), float(center[1])
    r_f = math.dist((focus[0], focus[1]), (cx, cy))
    theta_f = math.degrees(math.atan2(focus[0] - cx, focus[1] - cy))
    theta = np.arange(-ARC_LIMIT_DEG, ARC_LIMIT_DEG + ARC_STEP_DEG / 2, ARC_STEP_DEG)
    radii = r_f * np.linspace(1 - BAND, 1 + BAND, BAND_ARCS)

    xx, yy = np.meshgrid(x, y)
    geo_map = FieldGeometry(array, np.stack([xx, yy, np.full_like(xx, z)], axis=-1), **kw)
    geo_arc = FieldGeometry(array, _arc_points((cx, cy), r_f, theta, z), **kw)
    geo_band = FieldGeometry(array, np.concatenate([_arc_points((cx, cy), r, theta, z) for r in radii]), **kw)
    geo_focus = FieldGeometry(array, np.array([focus]), **kw)

    e_sim = np.zeros(geo_map.shape)
    e_foc = np.zeros(geo_map.shape)
    a_sim = np.zeros(len(theta))
    a_foc = np.zeros(len(theta))
    band_foc = np.zeros(geo_band.shape)
    focus_sim = focus_foc = 0.0
    contrast_k = np.full(len(freqs), np.nan)
    share10_k = np.zeros(len(freqs))
    win10 = np.abs(theta - theta_f) <= 10.0
    for k, (f, w) in enumerate(zip(freqs, weights)):
        f = float(f)
        ps, pf = geo_map.field(f), geo_map.field(f, emission)
        e_sim += w * (ps.real ** 2 + ps.imag ** 2)
        e_foc += w * (pf.real ** 2 + pf.imag ** 2)
        qs, qf = geo_arc.field(f), geo_arc.field(f, emission)
        ea_s, ea_f = qs.real ** 2 + qs.imag ** 2, qf.real ** 2 + qf.imag ** 2
        a_sim += w * ea_s
        a_foc += w * ea_f
        bf = geo_band.field(f, emission)
        band_foc += w * (bf.real ** 2 + bf.imag ** 2)
        s0, f0 = geo_focus.field(f)[0], geo_focus.field(f, emission)[0]
        focus_sim += w * abs(s0) ** 2
        focus_foc += w * abs(f0) ** 2
        lb = analyze_lobes(theta, 10 * np.log10(np.maximum(ea_f, 1e-6 * gain ** 2) / gain ** 2), reference_deg=theta_f)
        if lb["contrast_db"] is not None:
            contrast_k[k] = lb["contrast_db"]
        share10_k[k] = ea_f[win10].sum() / ea_f.sum()

    floor = params.floor_db
    lev_s = effective_level_db(e_sim, gain, floor)
    lev_f = effective_level_db(e_foc, gain, floor)
    diff = difference_db(np.sqrt(e_foc), np.sqrt(e_sim), gain, floor)
    arc_s = effective_level_db(a_sim, gain, floor)
    arc_f = effective_level_db(a_foc, gain, floor)
    l_s = float(effective_level_db(np.array(focus_sim), gain, floor))
    l_f = float(effective_level_db(np.array(focus_foc), gain, floor))

    lam_lo, lam_hi = c / bs.spectrum.f_max, c / bs.spectrum.f_min
    pitch = uniform_spacing(array)
    ind: dict = dict(
        scenario=f"{bs.series}/{bs.key}", f_min_hz=bs.spectrum.f_min, f_max_hz=bs.spectrum.f_max,
        distribution=bs.spectrum.distribution, n_lines=len(freqs),
        line_spacing_hz=0.0 if bs.spectrum.is_tone else float((freqs[1] - freqs[0]) if len(freqs) > 1 else 0.0),
        wavelength_min_mm=lam_lo * 1000, wavelength_max_mm=lam_hi * 1000,
        spacing_mm=None if pitch is None else pitch * 1000,
        spacing_over_wavelength_max=None if pitch is None else pitch / lam_lo,
        spacing_over_wavelength_min=None if pitch is None else pitch / lam_hi,
        aperture_mm=None if pitch is None else pitch * (len(array.elements) - 1) * 1000,
        focus_mm=list(bs.focus_mm), focus_range_mm=r_f * 1000, focus_angle_deg=theta_f,
        level_sim_focus_db=l_s, level_foc_focus_db=l_f, focusing_gain_db=l_f - l_s,
        weights_sum=float(weights.sum()))

    # arc at constant range through the focus, lobe containing the focus direction
    for tag, arc in (("foc", arc_f), ("sim", arc_s)):
        lb = analyze_lobes(theta, arc, reference_deg=theta_f)
        width, note = _beamwidth_if_defined(lb)
        ind.update({f"arc_{tag}_global_peak_db": float(arc.max()),
                    f"arc_{tag}_global_peak_deg": float(theta[int(np.argmax(arc))]),
                    f"arc_{tag}_lobe_peak_db": lb["peak_db"], f"arc_{tag}_lobe_peak_deg": lb["peak_deg"],
                    f"arc_{tag}_main_lobe_deg": [lb["null_left_deg"], lb["null_right_deg"]],
                    f"arc_{tag}_beamwidth3_deg": width, f"arc_{tag}_beamwidth_note": note,
                    f"arc_{tag}_sidelobe_db": lb["sidelobe_db"], f"arc_{tag}_sidelobe_deg": lb["sidelobe_deg"],
                    f"arc_{tag}_contrast_db": lb["contrast_db"]})
    for tag, e in (("foc", a_foc), ("sim", a_sim)):
        for w_deg in SHARE_WINDOWS_DEG:
            sel = np.abs(theta - theta_f) <= w_deg
            ind[f"arc_{tag}_share_{w_deg:g}deg"] = float(e[sel].sum() / e.sum())

    # sidelobes in the focal band: each circle excludes the lobe containing the focus direction
    band_e = band_foc.reshape(len(radii), len(theta))
    best = None
    for rad, e in zip(radii, band_e):
        lev = effective_level_db(e, gain, floor) + 20 * np.log10(rad / r_f)
        lb = analyze_lobes(theta, lev, reference_deg=theta_f)
        if lb["sidelobe_db"] is not None and (best is None or lb["sidelobe_db"] > best[0]):
            best = (lb["sidelobe_db"], float(rad), lb["sidelobe_deg"])
    ind.update(band_sidelobe_db=None, band_sidelobe_range_mm=None, band_sidelobe_deg=None,
               band_contrast_db=None)
    if best is not None:
        ind.update(band_sidelobe_db=best[0], band_sidelobe_range_mm=best[1] * 1000,
                   band_sidelobe_deg=best[2], band_contrast_db=l_f - best[0])

    # region statistics on the map: range-compensated levels in the focal band
    rr = np.hypot(xx - cx, yy - cy)
    roi = (rr >= (1 - BAND) * r_f) & (rr <= (1 + BAND) * r_f)
    comp = 20 * np.log10(np.maximum(rr, 1e-12) / r_f)
    ang = np.degrees(np.arctan2(xx - cx, yy - cy))
    dx = float(x[1] - x[0]) * 1000 if len(x) > 1 else 0.0
    dy = float(y[1] - y[0]) * 1000 if len(y) > 1 else 0.0
    ind["roi_cells"] = int(roi.sum())
    for tag, lev in (("foc", lev_f), ("sim", lev_s)):
        lc = lev + comp
        for t in THRESHOLDS_DB:
            ind[f"roi_pct_above_m{t:g}db_{tag}"] = float(100 * np.mean(lc[roi] >= l_f - t)) if roi.any() else None
    ind.update(focal_region_area_cm2=None, focal_region_radial_extent_mm=None,
               focal_region_angular_extent_deg=None)
    fx, fy = int(np.argmin(np.abs(x - focus[0]))), int(np.argmin(np.abs(y - focus[1])))
    near_focus = abs(x[fx] - focus[0]) <= dx / 2000 + 1e-12 and abs(y[fy] - focus[1]) <= dy / 2000 + 1e-12
    if roi.any() and near_focus:
        comp_mask = connected_component(roi & ((lev_f + comp) >= l_f - 3.0), (fy, fx))
        if comp_mask.any():
            ind.update(focal_region_area_cm2=float(comp_mask.sum() * dx * dy / 100),
                       focal_region_radial_extent_mm=float((rr[comp_mask].max() - rr[comp_mask].min()) * 1000),
                       focal_region_angular_extent_deg=float(ang[comp_mask].max() - ang[comp_mask].min()))
    iy, ix = np.unravel_index(int(np.argmax(lev_f)), lev_f.shape)
    ind.update(max_foc_db=float(lev_f[iy, ix]), max_foc_x_mm=float(x[ix] * 1000),
               max_foc_y_mm=float(y[iy] * 1000))

    # narrowband components of the same spectrum (arc analysis per line)
    defined = contrast_k[~np.isnan(contrast_k)]
    ind.update(
        narrowband_contrast_min_db=float(defined.min()) if defined.size else None,
        narrowband_contrast_median_db=float(np.median(defined)) if defined.size else None,
        narrowband_contrast_max_db=float(defined.max()) if defined.size else None,
        narrowband_fraction_below_3db=float(np.mean(defined < 3.0)) if defined.size else None,
        narrowband_share_10deg_min=float(share10_k.min()), narrowband_share_10deg_mean=float(share10_k.mean()),
        narrowband_share_10deg_max=float(share10_k.max()))
    return BroadbandResult(bs.as_scenario(), bs, array, x * 1000, y * 1000, e_sim, e_foc, lev_s, lev_f, diff,
                           ind, theta, arc_s, arc_f, freqs, weights, contrast_k, share10_k, [])


# ---- convergence ------------------------------------------------------------------------------------------------------

def convergence_study(bs: BroadbandScenario, region: Region, params: Params, dfs, reference_df: float,
                      array_dir: Path = ARRAY_DIR) -> list[dict]:
    """Max |difference| of the focused and difference maps (dB) against a finer spectral sampling."""
    array = load_array(array_dir / bs.array)
    x = grid_axis(region.x_mm[0] / 1000, region.x_mm[1] / 1000, region.nx, "x")
    y = grid_axis(region.y_mm[0] / 1000, region.y_mm[1] / 1000, region.ny, "y")
    focus = tuple(v / 1000 for v in bs.focus_mm)
    sp = bs.spectrum
    ref_s, ref_f = map_energies(array, focus, Spectrum(sp.f_min, sp.f_max, df=reference_df, distribution=sp.distribution),
                                x, y, region.z_mm / 1000, params)
    ref_f_db = effective_level_db(ref_f, params.gain, params.floor_db)
    ref_d = difference_db(np.sqrt(ref_f), np.sqrt(ref_s), params.gain, params.floor_db)
    rows = []
    for df in dfs:
        s2 = Spectrum(sp.f_min, sp.f_max, df=df, distribution=sp.distribution)
        e_s, e_f = map_energies(array, focus, s2, x, y, region.z_mm / 1000, params)
        f_db = effective_level_db(e_f, params.gain, params.floor_db)
        d = difference_db(np.sqrt(e_f), np.sqrt(e_s), params.gain, params.floor_db)
        rows.append(dict(scenario=f"{bs.series}/{bs.key}", df_hz=df, n_lines=len(s2.lines()[0]),
                         reference_df_hz=reference_df,
                         max_abs_focused_db=float(np.max(np.abs(f_db - ref_f_db))),
                         max_abs_difference_db=float(np.nanmax(np.abs(d - ref_d)))))
    return rows


# ---- output ---------------------------------------------------------------------------------------------------------------

def title_for(r: BroadbandResult) -> str:
    ind = r.indicators
    if r.scenario.series == "G":
        extra = f"d/lambda = {ind['spacing_over_wavelength_min']:.2f}-{ind['spacing_over_wavelength_max']:.2f}"
    elif r.scenario.series == "D":
        extra = f"angle {ind['focus_angle_deg']:+.1f} deg, range {ind['focus_range_mm']:.0f} mm"
    else:
        extra = f"d/lambda = {ind['spacing_over_wavelength_min']:.2f}-{ind['spacing_over_wavelength_max']:.2f}"
    return f"{r.scenario.label}\n{extra}"


def write_broadband(outdir: Path, r: BroadbandResult, params: Params, individual: bool = True,
                    flat: bool = False) -> None:
    d = outdir if flat else outdir / r.scenario.series / r.scenario.key
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "data.npz", x_mm=r.x_mm, y_mm=r.y_mm, energy_simultaneous=r.energy_sim,
                        energy_focused=r.energy_foc, simultaneous_db=r.level_sim, focused_db=r.level_foc,
                        difference_db=r.diff, arc_theta_deg=r.arc_theta, arc_simultaneous_db=r.arc_sim,
                        arc_focused_db=r.arc_foc, frequencies_hz=r.freqs, weights=r.weights,
                        per_frequency_contrast_db=r.per_frequency_contrast,
                        per_frequency_share10=r.per_frequency_share10, indicators=json.dumps(r.indicators))
    if not individual:
        return
    from .fieldplot import save_difference_map, save_level_map

    vmax, vmin = params.vmax, params.vmax - params.range_db
    sub = f"{r.scenario.label}"
    common = (r.array, r.scenario.focus_mm, r.x_mm, r.y_mm, 0.0)
    save_level_map(d / "field-simultaneous.png", r.level_sim, *common, f"Simultaneous (energy), {sub}", vmin, vmax)
    save_level_map(d / "field-focused.png", r.level_foc, *common, f"Focused (energy), {sub}", vmin, vmax)
    save_difference_map(d / "field-difference.png", r.diff, *common,
                        f"Focused - simultaneous (dB), {sub}", params.diff_limit)


def write_tables(outdir: Path, results: list[BroadbandResult]) -> None:
    rows = [_flat(r.indicators) for r in results]
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with (outdir / "broadband-indicators.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    (outdir / "broadband-indicators.json").write_text(
        json.dumps([r.indicators for r in results], indent=1) + "\n", encoding="utf-8")


PLATE_TITLES = {"B": "Series B - spectrum (8 speakers, 60 mm pitch, focus (350, 800) mm)",
                "G": "Series G - geometry (8 speakers; bands 2-4 kHz and 1-8 kHz; absolute focus (350, 800) mm)",
                "D": "Series D - focal direction (8 speakers, 60 mm pitch, 1-8 kHz, y = 800 mm)"}


def run_broadband_campaign(series: list[str], outdir: Path, region: Region, params: Params, df: float,
                           individual: bool = True, log=print, array_dir: Path = ARRAY_DIR):
    outdir.mkdir(parents=True, exist_ok=True)
    cache: dict = {}
    plates: dict[str, list[BroadbandResult]] = {}
    for s in series:
        results = []
        for bs in broadband_scenarios(s, df):
            ident = (bs.array, bs.focus_mm, bs.spectrum)
            if ident in cache:
                log(f"{bs.series}/{bs.key}: same as {cache[ident].scenario.series}/{cache[ident].scenario.key} (reused)")
                r = cache[ident]
            else:
                log(f"{bs.series}/{bs.key}: {bs.label}, {bs.spectrum.lines()[0].size} lines")
                r = run_broadband(bs, region, params, array_dir)
                cache[ident] = r
                write_broadband(outdir, r, params, individual)
            results.append(r)
        plates[s] = results
        plate = build_plate(PLATE_TITLES[s], results, params, titler=title_for)
        plate.fig.savefig(outdir / f"broadband-plate-{s}.png", dpi=140, bbox_inches="tight")
        build_arcs(PLATE_TITLES[s], results, titler=title_for).savefig(
            outdir / f"broadband-arcs-{s}.png", dpi=130, bbox_inches="tight")
    unique = list({id(r): r for rs in plates.values() for r in rs}.values())
    write_tables(outdir, unique)
    return unique


# ---- command line ------------------------------------------------------------------------------------------------------------------

def _print_summary(results: list[BroadbandResult]) -> None:
    print("scenario        L_foc  gain  arc_contrast  band_contrast  share10  roi>-6dB  area_cm2")
    for r in results:
        i = r.indicators
        f = lambda v, fmt="{:.2f}": "n/a" if v is None else fmt.format(v)
        print(f"{i['scenario']:<14} {i['level_foc_focus_db']:6.2f} {i['focusing_gain_db']:+5.2f}  "
              f"{f(i['arc_foc_contrast_db']):>11}  {f(i['band_contrast_db']):>12}  "
              f"{i['arc_foc_share_10deg']:7.3f}  {f(i['roi_pct_above_m6db_foc'], '{:.1f}'):>8}  "
              f"{f(i['focal_region_area_cm2'], '{:.1f}'):>8}")


def _spectrum_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--df", type=float, default=25.0, help="spectral line spacing in Hz")
    parser.add_argument("--vmax", type=float, default=20.0, help="top of the level color scale, dB (display only)")
    parser.add_argument("--range-db", type=float, default=40.0)
    parser.add_argument("--diff-limit", type=float, default=10.0)
    parser.add_argument("--speed-of-sound", type=float, default=SPEED_OF_SOUND)
    parser.add_argument("--gain", type=float, default=0.05, help="per-speaker amplitude at 1 m (sets the power)")


def broadband_map_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab broadband-map",
                                     description="Broadband energy maps (simultaneous, focused, difference)")
    parser.add_argument("--array", type=Path, required=True)
    parser.add_argument("--target", "--focus", dest="target", type=float, nargs=3, metavar=("X", "Y", "Z"),
                        required=True, help="focal point in mm")
    parser.add_argument("--fmin", type=float, required=True, help="lowest frequency, Hz")
    parser.add_argument("--fmax", type=float, required=True, help="highest frequency, Hz (= fmin for a tone)")
    parser.add_argument("--nfreq", type=int, help="number of spectral lines (overrides --df)")
    parser.add_argument("--distribution", choices=DISTRIBUTIONS, default="uniform",
                        help="power spectral density: uniform, or pink (1/f)")
    parser.add_argument("--x-range", type=float, nargs=2, default=(-700.0, 1300.0), metavar=("MIN", "MAX"))
    parser.add_argument("--y-range", type=float, nargs=2, default=(100.0, 2000.0), metavar=("MIN", "MAX"))
    parser.add_argument("--z", type=float, default=0.0, help="plane height in mm")
    parser.add_argument("--nx", type=int, default=401)
    parser.add_argument("--ny", type=int, default=381)
    parser.add_argument("--outdir", type=Path, required=True)
    _spectrum_args(parser)
    args = parser.parse_args(argv)
    try:
        if not args.range_db > 0 or not args.diff_limit > 0:
            raise ValueError("--range-db and --diff-limit must be positive")
        spectrum = Spectrum(args.fmin, args.fmax, df=args.df, n=args.nfreq, distribution=args.distribution)
        region = Region(x_mm=tuple(args.x_range), y_mm=tuple(args.y_range), nx=args.nx, ny=args.ny, z_mm=args.z)
        params = Params(speed_of_sound=args.speed_of_sound, gain=args.gain, vmax=args.vmax,
                        range_db=args.range_db, diff_limit=args.diff_limit)
        bs = BroadbandScenario("X", "map", f"{spectrum.label}", args.array.name, tuple(args.target), spectrum)
        result = run_broadband(bs, region, params, args.array.parent)
        write_broadband(args.outdir, result, params, flat=True)
    except (OSError, ValueError, KeyError, ImportError) as exc:
        parser.error(str(exc))
    print(f"Wrote broadband maps and data to {args.outdir}")
    i = result.indicators
    print(f"Spectrum {spectrum.label}, {i['n_lines']} lines, {args.distribution}; weights sum {i['weights_sum']:.12f}")
    _print_summary([result])
    return 0


def broadband_campaign_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab broadband-campaign",
                                     description="Broadband campaign: series B (spectrum), G (geometry), D (direction)")
    parser.add_argument("--series", nargs="+", choices=("B", "G", "D", "all"), default=["all"])
    parser.add_argument("--outdir", type=Path, default=Path("results/broadband"))
    parser.add_argument("--array-dir", type=Path, default=ARRAY_DIR)
    parser.add_argument("--nx", type=int, default=401)
    parser.add_argument("--ny", type=int, default=381)
    parser.add_argument("--no-individual", action="store_true")
    parser.add_argument("--convergence", action="store_true",
                        help="also write convergence.csv (spectral line spacing study on a map 4x coarser per axis)")
    _spectrum_args(parser)
    args = parser.parse_args(argv)
    series = ["B", "G", "D"] if "all" in args.series else list(dict.fromkeys(args.series))
    try:
        if not args.range_db > 0 or not args.diff_limit > 0:
            raise ValueError("--range-db and --diff-limit must be positive")
        region = Region(nx=args.nx, ny=args.ny)
        params = Params(speed_of_sound=args.speed_of_sound, gain=args.gain, vmax=args.vmax,
                        range_db=args.range_db, diff_limit=args.diff_limit)
        results = run_broadband_campaign(series, args.outdir, region, params, args.df,
                                         not args.no_individual, array_dir=args.array_dir)
        if args.convergence:
            coarse = Region(nx=max(11, (args.nx - 1) // 4 + 1), ny=max(11, (args.ny - 1) // 4 + 1))  # 101 x 96 by default
            rows = []
            for bs in broadband_scenarios("B", args.df)[1:]:
                rows += convergence_study(bs, coarse, params, (200.0, 100.0, 50.0, 25.0, 12.5), 6.25, args.array_dir)
            with (args.outdir / "convergence.csv").open("w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader()
                w.writerows(rows)
    except (OSError, ValueError, ImportError) as exc:
        parser.error(str(exc))
    print(f"Wrote plates, arcs, indicators and data to {args.outdir}")
    _print_summary(results)
    return 0
