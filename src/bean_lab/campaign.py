"""Parametric beamforming campaign: scenarios, indicators, comparison plates.

Reuses field.py (complex steady-state field) and fieldplot.py (axes-level drawing).
Conventions: millimeters in the array frame, top view of the plane z = region z; angles are
measured at the aperture center from the +y axis toward +x (0 deg = straight ahead).

Indicator definitions (all levels in dB re one speaker at 1 m):
  * focus level / focusing gain: field at the exact focus; gain = focused - simultaneous there.
  * global maximum: highest value over the whole region; usually in the near field.
  * arc: levels on the circle of radius R_f (aperture center -> focus) for |theta| <= 80 deg.
    A constant-range cut removes the geometric 1/r trend. Main lobe = lobe holding the arc peak,
    bounded by the first minima on each side; sidelobes = local maxima outside it;
    contrast = arc peak - highest sidelobe.
  * band sidelobe: the same arc analysis repeated on BAND_ARCS circles with radii from
    (1 - BAND) R_f to (1 + BAND) R_f. Each circle excludes its own main lobe (so the focused
    beam, an hourglass in the near field, is excluded as a whole) and levels are compensated for
    spherical spreading, L + 20 log10(r / R_f), so lobes at different ranges are comparable.
    The result is the highest compensated sidelobe in the band; contrast = focus level - it.
    These circles are not clipped to the image region. The band stays away from the near field.
Global maxima are searched on a grid `oversample` times finer than the images and refined
locally, so deep and narrow lobes at high frequency are not missed.
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
from .field import (DEFAULT_FLOOR_DB, complex_field, difference_db, field_maps, grid_axis,
                    level_db)
from .focus import SPEED_OF_SOUND, focus_delays

ARRAY_DIR = Path("configs/arrays")
BAND = 0.25  # sidelobe band: (1 - BAND) R_f .. (1 + BAND) R_f
BAND_ARCS = 21
ARC_LIMIT_DEG = 80.0
ARC_STEP_DEG = 0.02


@dataclass(frozen=True)
class Region:
    x_mm: tuple[float, float] = (-700.0, 1300.0)
    y_mm: tuple[float, float] = (100.0, 2000.0)
    nx: int = 401
    ny: int = 381
    z_mm: float = 0.0


@dataclass(frozen=True)
class Params:
    speed_of_sound: float = SPEED_OF_SOUND
    gain: float = 0.05
    vmax: float = 20.0  # display only
    range_db: float = 40.0  # display only
    diff_limit: float = 10.0  # display only
    floor_db: float = DEFAULT_FLOOR_DB
    oversample: int = 4  # indicator grid refinement factor


@dataclass(frozen=True)
class Scenario:
    series: str
    key: str
    label: str
    array: str  # file name under configs/arrays
    frequency: float
    focus_mm: tuple[float, float, float]


def scenarios(series: str) -> list[Scenario]:
    base = "line-8-60mm.json"
    if series == "F":
        return [Scenario("F", f"f{int(f):05d}", f"{f:g} Hz", base, float(f), (350.0, 800.0, 0.0))
                for f in (500, 1000, 2000, 4000, 6000, 8000, 12000, 16000)]
    if series == "S":
        return [Scenario("S", f"d{d:03d}", f"d = {d} mm (D = {7 * d} mm)", f"line-8-{d}mm.json",
                         4000.0, (350.0, 800.0, 0.0)) for d in (30, 60, 120)]
    if series == "D":
        return [Scenario("D", key, f"focus x = {x:g} mm", base, 4000.0, (x, 800.0, 0.0))
                for key, x in (("x+210", 210.0), ("x+350", 350.0), ("x+800", 800.0), ("x-380", -380.0))]
    raise ValueError(f"unknown series {series!r}")


# ---- geometry helpers ------------------------------------------------------------------

def aperture_center(array: ArrayConfig) -> np.ndarray:
    return np.mean([e.position for e in array.elements], axis=0)


def uniform_spacing(array: ArrayConfig) -> float | None:
    """Pitch in meters if the elements form a uniform line along x, else None."""
    pos = sorted((e.position for e in array.elements), key=lambda p: p[0])
    if len(pos) < 2 or any(abs(p[1] - pos[0][1]) > 1e-12 or abs(p[2] - pos[0][2]) > 1e-12 for p in pos):
        return None
    steps = np.diff([p[0] for p in pos])
    return float(steps[0]) if np.allclose(steps, steps[0], atol=1e-9) and steps[0] > 0 else None


def grating_angles(spacing: float, frequency: float, speed_of_sound: float, steer_deg: float,
                   max_order: int = 6) -> list[float]:
    """Far-field grating-lobe directions: sin(theta) = sin(steer) + m lambda / d, |sin| <= 1.

    An approximation here: the focus is in the near field of these apertures.
    """
    lam = speed_of_sound / frequency
    s0 = math.sin(math.radians(steer_deg))
    out = []
    for m in range(-max_order, max_order + 1):
        s = s0 + m * lam / spacing
        if m != 0 and abs(s) <= 1.0:
            out.append(math.degrees(math.asin(s)))
    return sorted(out)


# ---- lobe analysis ----------------------------------------------------------------------------

def analyze_lobes(theta: np.ndarray, level: np.ndarray, drop_db: float = 3.0,
                  reference_deg: float | None = None) -> dict:
    """Main lobe and sidelobes of a sampled angular pattern.

    The main lobe holds the global peak, or, if `reference_deg` is given, the lobe reached by
    climbing uphill from that direction (the intended beam direction).
    """
    n = len(level)
    if reference_deg is None:
        pk = int(np.argmax(level))
    else:
        pk = int(np.argmin(np.abs(theta - reference_deg)))
        while True:
            best = pk
            if pk > 0 and level[pk - 1] > level[best]:
                best = pk - 1
            if pk < n - 1 and level[pk + 1] > level[best]:
                best = pk + 1
            if best == pk:
                break
            pk = best
    lo = pk
    while lo > 0 and level[lo - 1] <= level[lo]:
        lo -= 1
    hi = pk
    while hi < n - 1 and level[hi + 1] <= level[hi]:
        hi += 1
    thr = level[pk] - drop_db
    a = pk
    while a > 0 and level[a - 1] >= thr:
        a -= 1
    b = pk
    while b < n - 1 and level[b + 1] >= thr:
        b += 1

    def cross(i_in: int, i_out: int) -> float:  # linear interpolation of the -drop crossing
        if i_in == i_out or level[i_in] == level[i_out]:
            return float(theta[i_in])
        t = (level[i_in] - thr) / (level[i_in] - level[i_out])
        return float(theta[i_in] + t * (theta[i_out] - theta[i_in]))

    left = cross(a, a - 1) if a > 0 else float(theta[0])
    right = cross(b, b + 1) if b < n - 1 else float(theta[-1])
    eps = 1e-9  # ignore floating-point noise on flat tails / the dB floor
    side = [i for i in range(1, n - 1) if (i < lo or i > hi)
            and level[i] > level[i - 1] + eps and level[i] >= level[i + 1] - eps]
    best = max(side, key=lambda i: level[i]) if side else None
    return dict(peak_deg=float(theta[pk]), peak_db=float(level[pk]), beamwidth_deg=right - left,
                beamwidth_open=bool(a == 0 or b == n - 1),
                null_left_deg=float(theta[lo]), null_right_deg=float(theta[hi]),
                main_lobe_open=bool(lo == 0 or hi == n - 1),
                sidelobe_deg=None if best is None else float(theta[best]),
                sidelobe_db=None if best is None else float(level[best]),
                contrast_db=None if best is None else float(level[pk] - level[best]),
                sidelobes=[(float(theta[i]), float(level[i])) for i in side])


def refine_peak(array: ArrayConfig, frequency: float, emission, kw: dict, x: float, y: float,
                z: float, half_width: float, bounds, *, n: int = 41, rounds: int = 4):
    """Local maximum of |P| near (x, y) by successive zooming (meters), kept inside `bounds`
    = (xmin, xmax, ymin, ymax) so a refined maximum never leaves the mapped region."""
    xmin, xmax, ymin, ymax = bounds
    for _ in range(rounds):
        gx = np.clip(np.linspace(x - half_width, x + half_width, n), xmin, xmax)
        gy = np.clip(np.linspace(y - half_width, y + half_width, n), ymin, ymax)
        xx, yy = np.meshgrid(gx, gy)
        pts = np.stack([xx, yy, np.full_like(xx, z)], axis=-1)
        a = np.abs(complex_field(array, pts, frequency, emission_delays=emission, **kw))
        iy, ix = np.unravel_index(int(np.argmax(a)), a.shape)
        x, y, best = float(gx[ix]), float(gy[iy]), float(a[iy, ix])
        half_width = 4 * (2 * half_width / (n - 1))
    return x, y, best


# ---- scenario analysis ---------------------------------------------------------------------------

@dataclass
class ScenarioResult:
    scenario: Scenario
    array: ArrayConfig
    x_mm: np.ndarray
    y_mm: np.ndarray
    simultaneous: np.ndarray  # complex (ny, nx)
    focused: np.ndarray
    level_sim: np.ndarray
    level_foc: np.ndarray
    diff: np.ndarray
    indicators: dict
    arc_theta: np.ndarray
    arc_sim: np.ndarray
    arc_foc: np.ndarray
    grating: list[float] = field(default_factory=list)


def _arc_points(center, radius, theta_deg, z):
    t = np.radians(theta_deg)
    return np.stack([center[0] + radius * np.sin(t), center[1] + radius * np.cos(t),
                     np.full_like(t, z)], axis=-1)


def run_scenario(sc: Scenario, region: Region, params: Params, array_dir: Path = ARRAY_DIR) -> ScenarioResult:
    array = load_array(array_dir / sc.array)
    f, c, gain = sc.frequency, params.speed_of_sound, params.gain
    kw = dict(speed_of_sound=c, gain=gain)
    focus = tuple(v / 1000 for v in sc.focus_mm)
    emission = {d.speaker: d.delay for d in focus_delays(array, focus, c)}
    z = region.z_mm / 1000
    x = grid_axis(region.x_mm[0] / 1000, region.x_mm[1] / 1000, region.nx, "x")
    y = grid_axis(region.y_mm[0] / 1000, region.y_mm[1] / 1000, region.ny, "y")
    maps = field_maps(array, focus, f, x, y, z, **kw)
    lev_s = level_db(maps.simultaneous, gain, params.floor_db)
    lev_f = level_db(maps.focused, gain, params.floor_db)
    diff = difference_db(maps.focused, maps.simultaneous, gain, params.floor_db)

    center = aperture_center(array)
    cx, cy = float(center[0]), float(center[1])
    r_f = math.dist((focus[0], focus[1]), (cx, cy))
    theta_f = math.degrees(math.atan2(focus[0] - cx, focus[1] - cy))
    lam = c / f
    pitch = uniform_spacing(array)
    ind: dict = dict(
        scenario=f"{sc.series}/{sc.key}", frequency_hz=f, wavelength_mm=lam * 1000,
        spacing_mm=None if pitch is None else pitch * 1000,
        spacing_over_wavelength=None if pitch is None else pitch / lam,
        aperture_mm=None if pitch is None else pitch * (len(array.elements) - 1) * 1000,
        focus_mm=list(sc.focus_mm), focus_range_mm=r_f * 1000, focus_angle_deg=theta_f)

    pf = complex_field(array, np.array(focus), f, emission_delays=emission, **kw)
    ps = complex_field(array, np.array(focus), f, **kw)
    l_f, l_s = float(level_db(pf, gain)), float(level_db(ps, gain))
    ind.update(level_sim_focus_db=l_s, level_foc_focus_db=l_f, focusing_gain_db=l_f - l_s)

    # arc at constant range
    theta = np.arange(-ARC_LIMIT_DEG, ARC_LIMIT_DEG + ARC_STEP_DEG / 2, ARC_STEP_DEG)
    arc_pts = _arc_points((cx, cy), r_f, theta, z)
    arc_f = level_db(complex_field(array, arc_pts, f, emission_delays=emission, **kw), gain, params.floor_db)
    arc_s = level_db(complex_field(array, arc_pts, f, **kw), gain, params.floor_db)
    lobes_f, lobes_s = analyze_lobes(theta, arc_f), analyze_lobes(theta, arc_s)
    for tag, lb in (("foc", lobes_f), ("sim", lobes_s)):
        ind.update({f"arc_{tag}_peak_db": lb["peak_db"], f"arc_{tag}_peak_deg": lb["peak_deg"],
                    f"arc_{tag}_beamwidth3_deg": lb["beamwidth_deg"],
                    f"arc_{tag}_beamwidth_open": lb["beamwidth_open"],
                    f"arc_{tag}_main_lobe_deg": [lb["null_left_deg"], lb["null_right_deg"]],
                    f"arc_{tag}_sidelobe_db": lb["sidelobe_db"],
                    f"arc_{tag}_sidelobe_deg": lb["sidelobe_deg"],
                    f"arc_{tag}_contrast_db": lb["contrast_db"]})
    grating = [] if pitch is None else grating_angles(pitch, f, c, theta_f)
    ind["grating_lobes_predicted_deg"] = grating

    # indicators on a finer grid
    os_ = params.oversample
    xf = grid_axis(x[0], x[-1], (region.nx - 1) * os_ + 1, "x")
    yf = grid_axis(y[0], y[-1], (region.ny - 1) * os_ + 1, "y")
    xx, yy = np.meshgrid(xf, yf)
    pts = np.stack([xx, yy, np.full_like(xx, z)], axis=-1)
    af = np.abs(complex_field(array, pts, f, emission_delays=emission, **kw))
    as_ = np.abs(complex_field(array, pts, f, **kw))
    del pts
    step = float(xf[1] - xf[0])
    bounds = (float(x[0]), float(x[-1]), float(y[0]), float(y[-1]))
    for tag, a, em in (("foc", af, emission), ("sim", as_, None)):
        iy, ix = np.unravel_index(int(np.argmax(a)), a.shape)
        gx, gy, val = refine_peak(array, f, em, kw, float(xf[ix]), float(yf[iy]), z, 1.5 * step, bounds)
        ind.update({f"max_{tag}_db": float(20 * np.log10(val / gain)),
                    f"max_{tag}_x_mm": gx * 1000, f"max_{tag}_y_mm": gy * 1000,
                    f"max_{tag}_range_mm": math.dist((gx, gy), (cx, cy)) * 1000,
                    f"max_{tag}_grid_db": float(20 * np.log10(a[iy, ix] / gain))})
    # sidelobes in the focal band, one circle per radius (each excludes its own main lobe)
    best = None
    for rad in r_f * np.linspace(1 - BAND, 1 + BAND, BAND_ARCS):
        pts_r = _arc_points((cx, cy), rad, theta, z)
        lev = level_db(complex_field(array, pts_r, f, emission_delays=emission, **kw), gain,
                       params.floor_db) + 20 * np.log10(rad / r_f)
        lb = analyze_lobes(theta, lev)
        if lb["sidelobe_db"] is not None and (best is None or lb["sidelobe_db"] > best[0]):
            best = (lb["sidelobe_db"], float(rad), lb["sidelobe_deg"])
    ind.update(band_sidelobe_db=None, band_sidelobe_range_mm=None, band_sidelobe_deg=None,
               band_sidelobe_x_mm=None, band_sidelobe_y_mm=None, band_contrast_db=None)
    if best is not None:
        db, rad, deg = best
        ind.update(band_sidelobe_db=db, band_sidelobe_range_mm=rad * 1000, band_sidelobe_deg=deg,
                   band_sidelobe_x_mm=(cx + rad * math.sin(math.radians(deg))) * 1000,
                   band_sidelobe_y_mm=(cy + rad * math.cos(math.radians(deg))) * 1000,
                   band_contrast_db=l_f - db)
    return ScenarioResult(sc, array, x * 1000, y * 1000, maps.simultaneous, maps.focused,
                          lev_s, lev_f, diff, ind, theta, arc_s, arc_f, grating)


# ---- plates ---------------------------------------------------------------------------------------------

def column_title(r: ScenarioResult) -> str:
    ind = r.indicators
    first = r.scenario.label
    if r.scenario.series == "F":
        extra = f"lambda = {ind['wavelength_mm']:.1f} mm, d/lambda = {ind['spacing_over_wavelength']:.2f}"
    elif r.scenario.series == "S":
        extra = f"d/lambda = {ind['spacing_over_wavelength']:.2f}"
    else:
        extra = f"angle {ind['focus_angle_deg']:+.1f} deg, range {ind['focus_range_mm']:.0f} mm"
    return f"{first}\n{extra}"


@dataclass
class Plate:
    fig: object
    level_images: list
    diff_images: list


def build_plate(title: str, results: list, params: Params, titler=None) -> Plate:
    from matplotlib.figure import Figure

    from .fieldplot import LEVEL_LABEL, decorate, draw_difference, draw_level

    titler = titler or column_title
    n = len(results)
    panel = 3.0
    fig = Figure(figsize=(panel * n + 1.4, panel * 3 + 1.3))
    gs = fig.add_gridspec(3, n + 1, width_ratios=[1] * n + [0.06], wspace=0.06, hspace=0.1,
                          top=0.90, bottom=0.06)
    vmax, vmin = params.vmax, params.vmax - params.range_db
    level_images, diff_images = [], []
    row_names = ("Simultaneous", "Focused", "Focused - simultaneous (dB)")
    for j, r in enumerate(results):
        for i in range(3):
            ax = fig.add_subplot(gs[i, j])
            if i == 0:
                im = draw_level(ax, r.level_sim, r.x_mm, r.y_mm, vmin, vmax)
                level_images.append(im)
            elif i == 1:
                im = draw_level(ax, r.level_foc, r.x_mm, r.y_mm, vmin, vmax)
                level_images.append(im)
            else:
                im = draw_difference(ax, r.diff, r.x_mm, r.y_mm, params.diff_limit)
                diff_images.append(im)
            decorate(ax, r.array, r.scenario.focus_mm, r.x_mm, r.y_mm,
                     titler(r) if i == 0 else None, legend=False, fontsize=8)
            ax.tick_params(labelsize=6)
            if j == 0:
                ax.set_ylabel(f"{row_names[i]}\ny (mm)", fontsize=8)
            else:
                ax.set_yticklabels([])
            if i == 2:
                ax.set_xlabel("x (mm)", fontsize=8)
            else:
                ax.set_xticklabels([])
    cax1 = fig.add_subplot(gs[0:2, n])
    cb1 = fig.colorbar(level_images[0], cax=cax1, extend="both")
    cb1.set_label(LEVEL_LABEL + f"; common scale {vmin:g}..{vmax:g} dB", fontsize=8)
    cb1.ax.tick_params(labelsize=7)
    cax2 = fig.add_subplot(gs[2, n])
    cb2 = fig.colorbar(diff_images[0], cax=cax2, extend="both")
    cb2.set_label("focused - simultaneous (dB)", fontsize=8)
    cb2.ax.tick_params(labelsize=7)
    fig.suptitle(f"{title}\nstar = focus, white squares = speakers (y = 0); free-field model, "
                 f"{params.gain:g} per speaker at 1 m; levels above/below the scale saturate",
                 fontsize=10, y=0.995)
    return Plate(fig, level_images, diff_images)


def build_arcs(title: str, results: list, titler=None):
    from matplotlib.figure import Figure

    titler = titler or column_title
    n = len(results)
    cols = min(n, 4)
    rows = math.ceil(n / cols)
    fig = Figure(figsize=(4.2 * cols, 3.2 * rows + 0.8))
    axes = fig.subplots(rows, cols, squeeze=False)
    top = max(float(max(r.arc_foc.max(), r.arc_sim.max())) for r in results) + 3
    for k, r in enumerate(results):
        ax = axes[k // cols][k % cols]
        ax.plot(r.arc_theta, r.arc_sim, "--", color="tab:gray", lw=1, label="simultaneous")
        ax.plot(r.arc_theta, r.arc_foc, "-", color="tab:red", lw=1.2, label="focused")
        ax.axvline(r.indicators["focus_angle_deg"], color="k", ls=":", lw=1, label="focus direction")
        for g in r.grating:
            ax.plot([g], [top - 50 + 1.5], "^", color="tab:blue", markersize=5)
        ax.set_ylim(top - 50, top)
        ax.set_xlim(-ARC_LIMIT_DEG, ARC_LIMIT_DEG)
        ax.set_title(titler(r).replace("\n", ", "), fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3)
        ax.set_xlabel("angle from +y axis (deg)", fontsize=8)
        ax.set_ylabel("dB re one speaker at 1 m", fontsize=8)
        if k == 0:
            ax.plot([], [], "^", color="tab:blue", label="far-field grating lobes (approx.)")
            ax.legend(fontsize=6, loc="lower left")
    for k in range(n, rows * cols):
        axes[k // cols][k % cols].axis("off")
    fig.suptitle(f"{title}: level on the circle of radius R_f (aperture center to focus)", fontsize=10)
    fig.tight_layout()
    return fig


# ---- output ---------------------------------------------------------------------------------------------------

def _flat(ind: dict) -> dict:
    out = {}
    for k, v in ind.items():
        if isinstance(v, list):
            out[k] = ";".join(f"{x:.2f}" if isinstance(x, float) else str(x) for x in v)
        else:
            out[k] = v
    return out


def write_scenario(outdir: Path, r: ScenarioResult, params: Params, individual: bool = True) -> None:
    d = outdir / r.scenario.series / r.scenario.key
    d.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(d / "data.npz", x_mm=r.x_mm, y_mm=r.y_mm, simultaneous=r.simultaneous,
                        focused=r.focused, simultaneous_db=r.level_sim, focused_db=r.level_foc,
                        difference_db=r.diff, arc_theta_deg=r.arc_theta, arc_simultaneous_db=r.arc_sim,
                        arc_focused_db=r.arc_foc, indicators=json.dumps(r.indicators))
    if not individual:
        return
    from .fieldplot import save_difference_map, save_level_map

    vmax, vmin = params.vmax, params.vmax - params.range_db
    sub = f"{r.scenario.series}: {r.scenario.label}, {r.scenario.frequency:g} Hz"
    common = (r.array, r.scenario.focus_mm, r.x_mm, r.y_mm, 0.0)
    save_level_map(d / "field-simultaneous.png", r.level_sim, *common, f"Simultaneous, {sub}", vmin, vmax)
    save_level_map(d / "field-focused.png", r.level_foc, *common, f"Focused, {sub}", vmin, vmax)
    save_difference_map(d / "field-difference.png", r.diff, *common,
                        f"Focused - simultaneous (dB), {sub}", params.diff_limit)


def write_tables(outdir: Path, results: list[ScenarioResult]) -> None:
    rows = [_flat(r.indicators) for r in results]
    with (outdir / "indicators.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (outdir / "indicators.json").write_text(
        json.dumps([r.indicators for r in results], indent=1) + "\n", encoding="utf-8")


def run_campaign(series: list[str], outdir: Path, region: Region, params: Params,
                 individual: bool = True, log=print, array_dir: Path = ARRAY_DIR) -> list[ScenarioResult]:
    outdir.mkdir(parents=True, exist_ok=True)
    all_results: list[ScenarioResult] = []
    for s in series:
        results = []
        for sc in scenarios(s):
            log(f"{sc.series}/{sc.key}: {sc.label}, focus {sc.focus_mm}")
            r = run_scenario(sc, region, params, array_dir)
            write_scenario(outdir, r, params, individual)
            results.append(r)
        titles = {"F": "Series F - frequency (8 speakers, 60 mm pitch, focus (350, 800) mm)",
                  "S": "Series S - spacing (8 speakers, 4 kHz, absolute focus (350, 800) mm)",
                  "D": "Series D - focal direction (8 speakers, 60 mm pitch, 4 kHz, y = 800 mm)"}
        plate = build_plate(titles[s], results, params)
        plate.fig.savefig(outdir / f"plate-{s}.png", dpi=140, bbox_inches="tight")
        build_arcs(titles[s], results).savefig(outdir / f"arcs-{s}.png", dpi=130, bbox_inches="tight")
        all_results += results
    write_tables(outdir, all_results)
    return all_results


def campaign_main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="bean-lab campaign",
                                     description="Parametric beamforming campaign (series F, S, D)")
    parser.add_argument("--series", nargs="+", choices=("F", "S", "D", "all"), default=["all"])
    parser.add_argument("--outdir", type=Path, default=Path("results/parametric"))
    parser.add_argument("--array-dir", type=Path, default=ARRAY_DIR)
    parser.add_argument("--nx", type=int, default=401)
    parser.add_argument("--ny", type=int, default=381)
    parser.add_argument("--oversample", type=int, default=4, help="indicator grid refinement factor")
    parser.add_argument("--vmax", type=float, default=20.0, help="top of the common level scale, dB")
    parser.add_argument("--range-db", type=float, default=40.0)
    parser.add_argument("--diff-limit", type=float, default=10.0)
    parser.add_argument("--no-individual", action="store_true", help="skip per-scenario PNGs")
    args = parser.parse_args(argv)
    series = ["F", "S", "D"] if "all" in args.series else list(dict.fromkeys(args.series))
    try:
        if args.oversample < 1 or not args.range_db > 0 or not args.diff_limit > 0:
            raise ValueError("--oversample >= 1, --range-db > 0 and --diff-limit > 0 are required")
        region = Region(nx=args.nx, ny=args.ny)
        params = Params(vmax=args.vmax, range_db=args.range_db, diff_limit=args.diff_limit,
                        oversample=args.oversample)
        results = run_campaign(series, args.outdir, region, params, not args.no_individual,
                               array_dir=args.array_dir)
    except (OSError, ValueError, ImportError) as exc:
        parser.error(str(exc))
    print(f"Wrote plates, arcs, indicators.csv/json and per-scenario data to {args.outdir}")
    print("scenario        L_sim   L_foc  gain   max_foc(pos mm)            arc: peak@deg  BW3  contrast")
    for r in results:
        i = r.indicators
        c = i["arc_foc_contrast_db"]
        print(f"{i['scenario']:<12} {i['level_sim_focus_db']:7.2f} {i['level_foc_focus_db']:7.2f} "
              f"{i['focusing_gain_db']:+6.2f}  {i['max_foc_db']:6.2f} ({i['max_foc_x_mm']:.0f},{i['max_foc_y_mm']:.0f})"
              f"   {i['arc_foc_peak_db']:6.2f}@{i['arc_foc_peak_deg']:+.1f}  {i['arc_foc_beamwidth3_deg']:5.1f}  "
              f"{'n/a' if c is None else f'{c:.2f}'}")
    return 0
