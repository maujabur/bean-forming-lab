"""PNG rendering of field maps (top view, XY plane). Requires matplotlib (extra: plot).

Only the color limits are adjustable; the data are never rescaled or normalized.
Axes-level drawing functions are shared by single-map PNGs and comparison plates.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .array_config import ArrayConfig

try:
    from matplotlib import colormaps
    from matplotlib.figure import Figure
except ImportError as exc:  # pragma: no cover
    raise ImportError("matplotlib is required for PNG output: pip install -e '.[plot]'") from exc

LEVEL_LABEL = "dB re one speaker at 1 m (model, not SPL)"
LEVEL_CMAP = "magma"
DIFF_CMAP = "RdBu_r"


def finite_max(*arrays: np.ndarray) -> float:
    return float(max(np.nanmax(a) for a in arrays))


def decorate(ax, array: ArrayConfig, focus_mm, x_mm, y_mm, title: str | None = None, *,
             mark_max: np.ndarray | None = None, legend: bool = True, note: str | None = None,
             fontsize: float = 10) -> None:
    """Speakers, focus, limits, aspect and labels (top view, x right, y up)."""
    sx = [e.position[0] * 1000 for e in array.elements]
    sy = [e.position[1] * 1000 for e in array.elements]
    ax.plot(sx, sy, "s", color="white", markeredgecolor="black", markersize=4, label="speakers")
    ax.plot([focus_mm[0]], [focus_mm[1]], "*", color="yellow", markeredgecolor="black",
            markersize=11, label="focus")
    if mark_max is not None:
        iy, ix = np.unravel_index(int(np.nanargmax(mark_max)), mark_max.shape)
        ax.plot([x_mm[ix]], [y_mm[iy]], "x", color="red", markersize=9, label="map maximum")
    pad = 0.04 * (y_mm[-1] - y_mm[0] if len(y_mm) > 1 else 1000.0)
    ax.set_ylim(min(y_mm[0], min(sy)) - pad, max(y_mm[-1], max(sy)) + pad)
    xs = list(x_mm) + sx
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title, fontsize=fontsize)
    if legend:
        ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
    if note:
        ax.text(0.01, 0.01, note, transform=ax.transAxes, fontsize=7, va="bottom",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))


def _extent(x_mm: np.ndarray, y_mm: np.ndarray) -> tuple[float, float, float, float]:
    dx = (x_mm[1] - x_mm[0]) / 2 if len(x_mm) > 1 else 0.5
    dy = (y_mm[1] - y_mm[0]) / 2 if len(y_mm) > 1 else 0.5
    return x_mm[0] - dx, x_mm[-1] + dx, y_mm[0] - dy, y_mm[-1] + dy


def draw_level(ax, levels_db: np.ndarray, x_mm, y_mm, vmin: float, vmax: float):
    if not vmin < vmax:
        raise ValueError("color range must have vmin < vmax")
    return ax.imshow(levels_db, origin="lower", extent=_extent(x_mm, y_mm), cmap=LEVEL_CMAP,
                     vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")


def draw_difference(ax, diff_db: np.ndarray, x_mm, y_mm, limit: float):
    if not limit > 0:
        raise ValueError("difference limit must be positive")
    # lightgray = undefined (both fields below the display floor)
    cmap = colormaps[DIFF_CMAP].with_extremes(bad="lightgray")
    return ax.imshow(np.ma.masked_invalid(diff_db), origin="lower", extent=_extent(x_mm, y_mm),
                     cmap=cmap, vmin=-limit, vmax=limit, interpolation="nearest", aspect="equal")


def _extend(levels_db: np.ndarray, vmin: float, vmax: float) -> str:
    """Colorbar arrows wherever data fall outside the display range (saturation)."""
    lo, hi = bool(np.nanmin(levels_db) < vmin), bool(np.nanmax(levels_db) > vmax)
    return "both" if lo and hi else "min" if lo else "max" if hi else "neither"


def save_level_map(path: Path, levels_db: np.ndarray, array: ArrayConfig, focus_mm, x_mm, y_mm,
                   z_mm: float, title: str, vmin: float, vmax: float, mark_max: bool = False) -> None:
    fig = Figure(figsize=(7.5, 8))
    ax = fig.subplots()
    im = draw_level(ax, levels_db, x_mm, y_mm, vmin, vmax)
    fig.colorbar(im, ax=ax, shrink=0.8, label=LEVEL_LABEL, extend=_extend(levels_db, vmin, vmax))
    decorate(ax, array, focus_mm, x_mm, y_mm, title, mark_max=levels_db if mark_max else None,
             note=f"plane z = {z_mm:g} mm; free-field point-source model")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    fig.savefig(path, dpi=130, bbox_inches="tight")


def save_difference_map(path: Path, diff_db: np.ndarray, array: ArrayConfig, focus_mm, x_mm, y_mm,
                        z_mm: float, title: str, limit: float) -> None:
    fig = Figure(figsize=(7.5, 8))
    ax = fig.subplots()
    im = draw_difference(ax, diff_db, x_mm, y_mm, limit)
    fig.colorbar(im, ax=ax, shrink=0.8, label="focused - simultaneous (dB)", extend="both")
    decorate(ax, array, focus_mm, x_mm, y_mm, title,
             note=f"plane z = {z_mm:g} mm; free-field point-source model")
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    fig.savefig(path, dpi=130, bbox_inches="tight")
