"""PNG rendering of field maps (top view, XY plane). Requires matplotlib (extra: plot).

Only the color limits are adjustable; the data are never rescaled or normalized.
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


def finite_max(*arrays: np.ndarray) -> float:
    return float(max(np.nanmax(a) for a in arrays))


def _decorate(ax, array: ArrayConfig, focus_mm, x_mm, y_mm, z_mm, title: str,
              mark_max: np.ndarray | None = None) -> None:
    sx = [e.position[0] * 1000 for e in array.elements]
    sy = [e.position[1] * 1000 for e in array.elements]
    ax.plot(sx, sy, "s", color="white", markeredgecolor="black", markersize=5, label="speakers")
    ax.plot([focus_mm[0]], [focus_mm[1]], "*", color="yellow", markeredgecolor="black",
            markersize=14, label="focus")
    if mark_max is not None:
        iy, ix = np.unravel_index(int(np.nanargmax(mark_max)), mark_max.shape)
        ax.plot([x_mm[ix]], [y_mm[iy]], "x", color="red", markersize=9, label="map maximum")
    pad = 0.04 * (y_mm[-1] - y_mm[0] if len(y_mm) > 1 else 1000.0)
    lo = min(y_mm[0], min(sy)) - pad
    hi = max(y_mm[-1], max(sy)) + pad
    ax.set_ylim(lo, hi)
    xs = list(x_mm) + sx
    ax.set_xlim(min(xs) - pad, max(xs) + pad)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=10)
    ax.legend(loc="upper right", fontsize=7, framealpha=0.85)
    ax.text(0.01, 0.01, f"plane z = {z_mm:g} mm; free-field point-source model",
            transform=ax.transAxes, fontsize=7, va="bottom",
            bbox=dict(facecolor="white", alpha=0.7, edgecolor="none"))


def _save(fig, path: Path) -> None:
    fig.savefig(path, dpi=130, bbox_inches="tight")


def _extent(x_mm: np.ndarray, y_mm: np.ndarray) -> tuple[float, float, float, float]:
    dx = (x_mm[1] - x_mm[0]) / 2 if len(x_mm) > 1 else 0.5
    dy = (y_mm[1] - y_mm[0]) / 2 if len(y_mm) > 1 else 0.5
    return x_mm[0] - dx, x_mm[-1] + dx, y_mm[0] - dy, y_mm[-1] + dy


def save_level_map(path: Path, levels_db: np.ndarray, array: ArrayConfig, focus_mm, x_mm, y_mm,
                   z_mm: float, title: str, vmin: float, vmax: float, mark_max: bool = False) -> None:
    if not vmin < vmax:
        raise ValueError("color range must have vmin < vmax")
    fig = Figure(figsize=(7.5, 8))
    ax = fig.subplots()
    im = ax.imshow(levels_db, origin="lower", extent=_extent(x_mm, y_mm), cmap="magma",
                   vmin=vmin, vmax=vmax, interpolation="nearest", aspect="equal")
    fig.colorbar(im, ax=ax, shrink=0.8, label=LEVEL_LABEL, extend="min")
    _decorate(ax, array, focus_mm, x_mm, y_mm, z_mm, title, levels_db if mark_max else None)
    _save(fig, path)


def save_difference_map(path: Path, diff_db: np.ndarray, array: ArrayConfig, focus_mm, x_mm, y_mm,
                        z_mm: float, title: str, limit: float) -> None:
    if not limit > 0:
        raise ValueError("difference limit must be positive")
    fig = Figure(figsize=(7.5, 8))
    ax = fig.subplots()
    # lightgray = undefined (both fields below the display floor)
    cmap = colormaps["RdBu_r"].with_extremes(bad="lightgray")
    im = ax.imshow(np.ma.masked_invalid(diff_db), origin="lower", extent=_extent(x_mm, y_mm),
                   cmap=cmap, vmin=-limit, vmax=limit, interpolation="nearest", aspect="equal")
    fig.colorbar(im, ax=ax, shrink=0.8, label="focused - simultaneous (dB)", extend="both")
    _decorate(ax, array, focus_mm, x_mm, y_mm, z_mm, title)
    _save(fig, path)
