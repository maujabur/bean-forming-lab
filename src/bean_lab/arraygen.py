"""Reproducible generators of planar loudspeaker layouts, with rigid mechanical constraints.

Coordinates: millimeters, panel centered at the origin, array in the plane z = 0, normal +z
(the project convention of Mission 06). A layout is a (N, 2) array of x, y centers.

Mechanical model (a study hypothesis, NOT a measured dimension of the hardware): every driver is a
disc of `diameter` mm (default 65), so
  * centers are at least `min_distance` apart (default = diameter), and
  * centers are at least `margin` from every panel edge (default = diameter / 2),
  i.e. centers lie in the inner rectangle [-w/2 + margin, w/2 - margin] x [-h/2 + margin, h/2 - margin].
Replace the diameter by the real driver size before any manufacturing decision.

A generator never moves a position to repair a violation. It either returns a valid layout or raises
InfeasibleGeometry. The only documented adjustments are explicit options: `adjust="scale"` for the
Vogel spiral (uniform radial scaling, factor returned in the metadata) and the repulsion relaxation
used only to seed the optimizer (`relaxed_random`). The final layout is always re-validated.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .array_config import ArrayConfig, Element

TOL = 1e-9
GOLDEN_ANGLE = math.pi * (3 - math.sqrt(5))


class InfeasibleGeometry(ValueError):
    """The requested layout cannot satisfy the mechanical constraints."""


@dataclass(frozen=True)
class Panel:
    width: float = 420.0
    height: float = 420.0
    margin: float = 32.5
    min_distance: float = 65.0
    diameter: float = 65.0  # reference outer diameter of a driver (hypothesis)

    def __post_init__(self) -> None:
        for name in ("width", "height", "margin", "min_distance", "diameter"):
            v = getattr(self, name)
            if not (isinstance(v, (int, float)) and math.isfinite(v)) or v < 0:
                raise ValueError(f"{name} must be a finite nonnegative number")
        if self.width <= 0 or self.height <= 0 or self.min_distance <= 0:
            raise ValueError("width, height and min_distance must be positive")
        if 2 * self.margin >= min(self.width, self.height):
            raise ValueError("margin leaves no room for any center")

    @property
    def inner(self) -> tuple[float, float, float, float]:
        return (-self.width / 2 + self.margin, self.width / 2 - self.margin,
                -self.height / 2 + self.margin, self.height / 2 - self.margin)

    @property
    def span(self) -> tuple[float, float]:
        x0, x1, y0, y1 = self.inner
        return x1 - x0, y1 - y0


def check_positions(pos: np.ndarray, panel: Panel) -> dict:
    """Constraint report: min center distance, min clearance to the panel edges, and violations."""
    pos = np.asarray(pos, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 2 or not np.all(np.isfinite(pos)) or len(pos) < 1:
        raise ValueError("positions must be a finite (N, 2) array")
    x0, x1, y0, y1 = panel.inner
    d = np.linalg.norm(pos[:, None, :] - pos[None, :, :], axis=-1) + np.diag(np.full(len(pos), np.inf))
    min_d = float(d.min()) if len(pos) > 1 else math.inf
    clear = float(min((pos[:, 0] - (-panel.width / 2)).min(), (panel.width / 2 - pos[:, 0]).min(),
                      (pos[:, 1] - (-panel.height / 2)).min(), (panel.height / 2 - pos[:, 1]).min()))
    problems = []
    if min_d < panel.min_distance - TOL:
        problems.append(f"minimum center distance {min_d:.3f} mm < {panel.min_distance:g} mm")
    if clear < panel.margin - TOL:
        problems.append(f"edge clearance {clear:.3f} mm < margin {panel.margin:g} mm")
    return dict(ok=not problems, n=len(pos), min_distance_mm=min_d, edge_clearance_mm=clear, problems=problems)


def validate_positions(pos: np.ndarray, panel: Panel) -> np.ndarray:
    report = check_positions(pos, panel)
    if not report["ok"]:
        raise InfeasibleGeometry("; ".join(report["problems"]))
    return np.asarray(pos, dtype=float)


def _grid_axis(count: int, lo: float, hi: float) -> np.ndarray:
    if count < 1:
        raise ValueError("grid counts must be positive")
    return np.array([(lo + hi) / 2]) if count == 1 else np.linspace(lo, hi, count)


def regular_grid(rows: int, cols: int, panel: Panel) -> np.ndarray:
    """rows x cols grid spread over the whole inner rectangle (the available envelope)."""
    x0, x1, y0, y1 = panel.inner
    xs, ys = _grid_axis(cols, x0, x1), _grid_axis(rows, y0, y1)
    gx, gy = np.meshgrid(xs, ys)
    return validate_positions(np.stack([gx.ravel(), gy.ravel()], axis=1), panel)


def jittered_grid(rows: int, cols: int, panel: Panel, max_shift: float, seed: int, attempts: int = 500) -> np.ndarray:
    """Regular grid with each element displaced by at most `max_shift` mm (uniform in the disc).

    Elements are moved one at a time in a random order and a displacement is accepted only if it keeps
    the element inside the panel limits and `min_distance` from every other element in its current
    state; after `attempts` rejections an element stays put only if its grid position is itself valid,
    otherwise the layout is rejected. The base grid must be valid: nothing is repaired.
    """
    if max_shift < 0 or not math.isfinite(max_shift):
        raise ValueError("max_shift must be finite and nonnegative")
    base = regular_grid(rows, cols, panel)
    rng = np.random.default_rng(seed)
    pos = base.copy()
    x0, x1, y0, y1 = panel.inner
    for i in rng.permutation(len(pos)):
        for _ in range(attempts):
            r, a = max_shift * math.sqrt(rng.random()), 2 * math.pi * rng.random()
            cand = base[i] + r * np.array([math.cos(a), math.sin(a)])
            if not (x0 - TOL <= cand[0] <= x1 + TOL and y0 - TOL <= cand[1] <= y1 + TOL):
                continue
            others = np.delete(pos, i, axis=0)
            if len(others) == 0 or np.linalg.norm(others - cand, axis=1).min() >= panel.min_distance - TOL:
                pos[i] = cand
                break
    return validate_positions(pos, panel)


def vogel(n: int, panel: Panel, radius: float | None = None, adjust: str = "none", phase_deg: float = 0.0):
    """Vogel (sunflower) spiral: r_i = R sqrt((i + 0.5) / N), theta_i = i * golden angle (+ phase).

    radius defaults to the disc inscribed in the inner rectangle. adjust="none": the spiral is validated
    as generated. adjust="scale": if the minimum distance is violated, all positions are scaled radially
    about the center by min_distance / min_distance_found (documented procedure) and re-validated.
    Returns (positions, metadata) with the scale factor actually applied.
    """
    if n < 1 or adjust not in ("none", "scale"):
        raise ValueError("n must be positive and adjust 'none' or 'scale'")
    sx, sy = panel.span
    R = min(sx, sy) / 2 if radius is None else float(radius)
    if not R > 0:
        raise ValueError("radius must be positive")
    i = np.arange(n)
    r = R * np.sqrt((i + 0.5) / n)
    th = i * GOLDEN_ANGLE + math.radians(phase_deg)
    pos = np.stack([r * np.cos(th), r * np.sin(th)], axis=1)
    scale = 1.0
    if adjust == "scale" and n > 1:
        d = check_positions(pos, Panel(panel.width, panel.height, panel.margin, 1e-9, panel.diameter))["min_distance_mm"]
        if d < panel.min_distance:
            scale = panel.min_distance / d * (1 + 1e-9)
            pos = pos * scale
    return validate_positions(pos, panel), dict(radius_mm=R, scale=scale, phase_deg=phase_deg)


def rings(counts: list[int], radii: list[float], panel: Panel, offsets_deg: list[float] | None = None) -> np.ndarray:
    """Concentric rings: ring k has counts[k] elements at radius radii[k] starting at offsets_deg[k]."""
    if len(counts) != len(radii) or not counts or any(c < 1 for c in counts):
        raise ValueError("counts and radii must be nonempty, equal length, counts >= 1")
    if offsets_deg is None:  # stagger rings by half an angular step of the previous one
        offsets_deg = [0.0] + [180.0 / counts[k] for k in range(1, len(counts))]
    pts = []
    for c, r, o in zip(counts, radii, offsets_deg):
        a = math.radians(o) + 2 * math.pi * np.arange(c) / c
        pts.append(np.stack([r * np.cos(a), r * np.sin(a)], axis=1))
    return validate_positions(np.concatenate(pts), panel)


def rings_auto(n: int, panel: Panel, n_rings: int | None = None, radius: float | None = None):
    """Rings with radii R k/m and counts proportional to the radius; tries m = 2, 3, 4 unless given.

    Returns (positions, metadata with the ring structure); raises InfeasibleGeometry if no m is valid.
    """
    sx, sy = panel.span
    R = min(sx, sy) / 2 if radius is None else float(radius)
    errors = []
    for m in ([n_rings] if n_rings else (2, 3, 4)):
        if m > n:
            continue
        weights = np.arange(1, m + 1, dtype=float)
        counts = np.maximum(1, np.floor(n * weights / weights.sum()).astype(int))
        for k in np.argsort(-(n * weights / weights.sum() - counts))[: n - counts.sum()]:
            counts[k] += 1
        radii = [R * (k + 1) / m for k in range(m)]
        try:
            return rings([int(c) for c in counts], radii, panel), dict(n_rings=m, counts=[int(c) for c in counts],
                                                                        radii_mm=radii)
        except InfeasibleGeometry as exc:
            errors.append(f"{m} rings: {exc}")
    raise InfeasibleGeometry("no ring structure is valid (" + "; ".join(errors) + ")")


def relaxed_random(n: int, panel: Panel, seed: int, iterations: int = 5000) -> np.ndarray:
    """Random positions pushed apart by pairwise repulsion until valid (seed for the optimizer only)."""
    rng = np.random.default_rng(seed)
    x0, x1, y0, y1 = panel.inner
    pos = np.stack([rng.uniform(x0, x1, n), rng.uniform(y0, y1, n)], axis=1)
    target = panel.min_distance * (1 + 1e-6)
    for _ in range(iterations):
        diff = pos[:, None, :] - pos[None, :, :]
        d = np.linalg.norm(diff, axis=-1) + np.eye(n) * 1e9
        close = d < target
        if not close.any():
            break
        push = np.where(close[..., None], diff / np.maximum(d, 1e-9)[..., None] * (target - d)[..., None] * 0.5, 0.0)
        pos = pos + push.sum(axis=1) + rng.normal(0, 1e-3, pos.shape)
        pos[:, 0], pos[:, 1] = np.clip(pos[:, 0], x0, x1), np.clip(pos[:, 1], y0, y1)
    return validate_positions(pos, panel)


@dataclass(frozen=True)
class GeometrySpec:
    family: str  # R, J, V, C, O (O only holds optimizer output)
    n: int
    panel: Panel
    seed: int = 0
    rows: int | None = None
    cols: int | None = None
    max_shift: float | None = None
    vogel_adjust: str = "scale"
    n_rings: int | None = None

    @property
    def key(self) -> str:
        return f"{self.family}{self.n}-{self.panel.width:g}x{self.panel.height:g}" + (
            f"-s{self.seed}" if self.family == "J" else "")


def default_grid(n: int) -> tuple[int, int]:
    """rows x cols of the regular reference grids: 2x4, 4x4, 4x8 (the project convention); otherwise the
    most square factorization with rows <= cols."""
    table = {8: (2, 4), 16: (4, 4), 32: (4, 8)}
    if n in table:
        return table[n]
    rows = max(r for r in range(1, int(math.isqrt(n)) + 1) if n % r == 0)
    return rows, n // rows


def build(spec: GeometrySpec):
    """(positions, metadata). Raises InfeasibleGeometry for an unsatisfiable request."""
    rows, cols = (spec.rows, spec.cols) if spec.rows and spec.cols else default_grid(spec.n)
    if spec.family in ("R", "J") and rows * cols != spec.n:
        raise ValueError("rows x cols must equal the number of elements")
    if spec.family == "R":
        return regular_grid(rows, cols, spec.panel), dict(rows=rows, cols=cols)
    if spec.family == "J":
        step = min(*(s / max(c - 1, 1) for s, c in zip(spec.panel.span, (cols, rows))))
        shift = spec.max_shift if spec.max_shift is not None else 0.25 * step
        return jittered_grid(rows, cols, spec.panel, shift, spec.seed), dict(rows=rows, cols=cols, max_shift_mm=shift,
                                                                              seed=spec.seed)
    if spec.family == "V":
        return vogel(spec.n, spec.panel, adjust=spec.vogel_adjust)
    if spec.family == "C":
        return rings_auto(spec.n, spec.panel, spec.n_rings)
    raise ValueError(f"family {spec.family!r} cannot be built directly (R, J, V, C)")


def to_array_config(name: str, pos: np.ndarray) -> ArrayConfig:
    return ArrayConfig(name, tuple(Element(i + 1, (float(x) / 1000, float(y) / 1000, 0.0)) for i, (x, y) in enumerate(pos)))


def positions_of(array: ArrayConfig) -> np.ndarray:
    """(N, 2) millimeters of a planar array (z must be 0)."""
    p = np.array([e.position for e in array.elements]) * 1000.0
    if np.any(np.abs(p[:, 2]) > 1e-9):
        raise ValueError("array is not in the z = 0 plane")
    return p[:, :2]


def save_array_json(path: str | Path, name: str, pos: np.ndarray, panel: Panel, meta: dict | None = None) -> None:
    """Array JSON compatible with load_array; extra keys (panel, meta) are ignored by the loader."""
    obj = {"name": name, "units": "mm",
           "panel": {"width": panel.width, "height": panel.height, "margin": panel.margin,
                     "min_distance": panel.min_distance, "diameter": panel.diameter,
                     "note": "driver diameter is a study hypothesis, not a measured dimension"},
           "meta": meta or {},
           "elements": [{"speaker": i + 1, "position": [float(x), float(y), 0.0]}
                        for i, (x, y) in enumerate(pos)]}
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=1) + "\n", encoding="utf-8")
