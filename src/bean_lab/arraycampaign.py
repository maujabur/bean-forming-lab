"""Campaigns comparing layout families, number of elements, aperture and optimized layouts.

Campaign A: families R, J, V, C with 8/16/32 elements on the 420 x 420 mm reference panel.
Campaign B: the same families at a common, larger 520 x 520 mm panel (the smallest where a 4 x 8 grid of
            65 mm drivers fits), 8/16/32 elements: effect of the number of elements at fixed envelope.
Campaign C: fixed element count (8, 16), families R and V, panels 260 to 840 mm: effect of the aperture.
Campaign D: optimized layouts (positions only) trained on a set of foci, validated on another.

Infeasible combinations are recorded with the reason and never enter the acoustic evaluation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .arrayeval import (A_REF, ALL_FOCI, BANDS, REFERENCE_ROI, TRAIN_FOCI, VALIDATION_FOCI, element_gain, evaluate_reference,
                        pareto_mask, summarize)
from .broadband import Spectrum
from .arraygen import (GeometrySpec, InfeasibleGeometry, Panel, build, check_positions, save_array_json,
                       to_array_config)
from .arrayopt import OptConfig, optimize_multi

ENVELOPES = {"E260": Panel(260, 260), "E420": Panel(420, 420), "E520": Panel(520, 520), "E630": Panel(630, 630),
             "E840": Panel(840, 840)}
FAMILY_COLORS = {"R": "tab:blue", "J": "tab:orange", "V": "tab:green", "C": "tab:red", "O": "black", "Ow": "tab:purple"}
BAND_LABEL = {"b1535": "1.5-3.5 kHz", "b24": "2-4 kHz", "b18": "1-8 kHz"}
TRAIN_NAMES = [f.name for f in TRAIN_FOCI]
VAL_NAMES = [f.name for f in VALIDATION_FOCI]


@dataclass
class Entry:
    id: str
    family: str  # R J V C O Ow
    n: int
    envelope: str
    campaigns: set
    seed: int = 0
    status: str = "ok"  # ok | infeasible
    reason: str = ""
    positions: np.ndarray | None = None
    meta: dict = field(default_factory=dict)
    ref: dict | None = None
    summary: dict | None = None  # {band: {"train": {...}, "val": {...}, "all": {...}}}
    seconds: float = 0.0

    @property
    def panel(self) -> Panel:
        return ENVELOPES[self.envelope]

    @property
    def area_mm2(self) -> float:
        return self.panel.width * self.panel.height


@dataclass(frozen=True)
class Profile:
    name: str
    ns: tuple
    j_seeds: tuple
    c_ns: tuple
    c_envelopes: tuple
    df: float
    opt_iterations: int
    restarts: tuple  # (init, seed)
    opt_ns: tuple
    wide_opt: bool
    cuts_df: float


PROFILES = {
    "quick": Profile("quick", (8, 16), (1,), (8, 16), ("E260", "E420", "E630"), 200.0, 800,
                     (("rand", 1), ("V", 0)), (8, 16), False, 200.0),
    "full": Profile("full", (8, 16, 32), (1, 2, 3), (8, 16), ("E260", "E420", "E630", "E840"), 50.0, 10000,
                    (("rand", 1), ("rand", 2), ("V", 0), ("C", 0), ("J", 1)), (8, 16, 32), True, 100.0),
}


def _make(family: str, n: int, env: str, seed: int, campaign: str) -> Entry:
    panel = ENVELOPES[env]
    e = Entry(f"{family}{n}-{env}" + (f"-s{seed}" if family == "J" else ""), family, n, env, {campaign}, seed)
    try:
        e.positions, meta = build(GeometrySpec(family, n, panel, seed=seed))
        e.meta = {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in meta.items()}
        e.meta["min_distance_mm"] = check_positions(e.positions, panel)["min_distance_mm"]
    except InfeasibleGeometry as exc:
        e.status, e.reason = "infeasible", str(exc)
    return e


def catalogue(profile: Profile) -> list[Entry]:
    entries: dict[str, Entry] = {}

    def add(e: Entry):
        if e.id in entries:
            entries[e.id].campaigns |= e.campaigns
        else:
            entries[e.id] = e

    for n in profile.ns:
        for fam in "RVC":
            add(_make(fam, n, "E420", 0, "A"))
            add(_make(fam, n, "E520", 0, "B"))
        for s in profile.j_seeds:
            add(_make("J", n, "E420", s, "A"))
            add(_make("J", n, "E520", s, "B"))
    for n in profile.c_ns:
        for env in profile.c_envelopes:
            for fam in "RV":
                add(_make(fam, n, env, 0, "C"))
    return list(entries.values())


def _hash(e: Entry, df: float) -> str:
    h = hashlib.sha256()
    h.update(np.round(e.positions, 9).tobytes())
    h.update(f"{df}|{REFERENCE_ROI}|{[f.name for f in ALL_FOCI]}".encode())
    return h.hexdigest()[:16]


def _cell_key(k) -> str:
    return f"{k[0]}|{k[1]}"


def evaluate_entry(e: Entry, df: float, cache_dir: Path | None, log=print) -> None:
    """Reference-resolution evaluation (all foci, all bands), cached on the layout content."""
    key = _hash(e, df)
    path = cache_dir / f"{e.id}.json" if cache_dir else None
    if path and path.exists():
        data = json.loads(path.read_text())
        if data.get("key") == key:
            e.ref = dict(cells={tuple(k.split("|")): v for k, v in data["cells"].items()}, seconds=data["seconds"], n=e.n)
            e.seconds = data["seconds"]
            log(f"  {e.id}: cached")
            e.summary = _summaries(e)
            return
    arr = to_array_config(e.id, e.positions)
    e.ref = evaluate_reference(arr, df=df, keep_maps=True)
    for cell in e.ref["cells"].values():
        if "level_map_db" in cell:
            cell["level_map_db"] = np.round(cell["level_map_db"], 2)
    e.seconds = e.ref["seconds"]
    e.summary = _summaries(e)
    log(f"  {e.id}: reference evaluation {e.seconds:.1f} s")
    if path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(dict(key=key, seconds=e.seconds, cells={_cell_key(k): v for k, v in e.ref["cells"].items()}),
                                   default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o)))


def _summaries(e: Entry) -> dict:
    out = {}
    for band in BANDS:
        out[band] = {name: summarize(e.ref, names, band) for name, names in
                     (("train", TRAIN_NAMES), ("val", VAL_NAMES), ("all", TRAIN_NAMES + VAL_NAMES))}
    return out


def objectives(e: Entry) -> dict:
    """Pareto coordinates: selectivity and robustness are means over the three bands on the validation foci."""
    val = [e.summary[b]["val"] for b in BANDS]
    return dict(selectivity=float(np.mean([v["F10_mean"] for v in val])),
                robust_focus=float(np.mean([v["F10_worst"] for v in val])),
                robust_lines=float(np.mean([v["F10_p10_lines_worst"] for v in val])),
                n=e.n, area_mm2=e.area_mm2)


def pareto_flags(entries: list[Entry], keys=("selectivity", "robust_focus", "n", "area_mm2"),
                 maximize=(True, True, False, False)) -> dict[str, bool]:
    ok = [e for e in entries if e.status == "ok" and e.summary]
    pts = np.array([[objectives(e)[k] for k in keys] for e in ok])
    mask = pareto_mask(pts, list(maximize))
    return {e.id: bool(m) for e, m in zip(ok, mask)}


# ---- optimization campaign (D) ------------------------------------------------------------------------------------

def run_optimization(profile: Profile, outdir: Path, log=print) -> tuple[list[Entry], list[dict]]:
    entries, restarts_log = [], []
    jobs = [(n, "E420", "b24", False) for n in profile.opt_ns] + [(32, "E520", "b24", False)]
    if profile.wide_opt:
        jobs += [(n, "E420", "b18", True) for n in profile.opt_ns if n <= 16]
    for n, env, band, wide in jobs:
        if n == 32 and env == "E520" and 32 not in profile.opt_ns:
            continue
        cfg = OptConfig(iterations=profile.opt_iterations, band=BANDS[band], n_lines=9 if band != "b18" else 28)
        panel = ENVELOPES[env]
        log(f"optimize N={n} {env} train band {BAND_LABEL[band]} ({len(profile.restarts)} restarts, {cfg.iterations} iterations)")
        results, skipped = optimize_multi(n, panel, list(profile.restarts), cfg, log)
        for r in results:
            restarts_log.append(dict(n=n, envelope=env, train_band=band, init=r.init, seed=r.seed, j_start=r.j_start,
                                     j_best=r.j_best, evaluations=r.evaluations, seconds=r.seconds, accepted=r.accepted,
                                     infeasible_proposals=r.infeasible_proposals, memory_bytes=r.memory_bytes,
                                     components=r.components_best))
        for init, seed, why in skipped:
            restarts_log.append(dict(n=n, envelope=env, train_band=band, init=init, seed=seed, skipped=why))
        if not results:
            continue
        best = min(results, key=lambda r: r.j_best)
        fam = "Ow" if wide else "O"
        e = Entry(f"{fam}{n}-{env}", fam, n, env, {"D"}, seed=best.seed, positions=best.positions,
                  meta=dict(init=best.init, seed=best.seed, train_band=band, j_best=best.j_best,
                            min_distance_mm=check_positions(best.positions, panel)["min_distance_mm"]))
        entries.append(e)
    return entries, restarts_log


# ---- plots --------------------------------------------------------------------------------------------------------------

def plot_layouts(entries: list[Entry], path: Path, title: str, cols: int = 4, subtitles: dict | None = None) -> None:
    from matplotlib.figure import Figure
    from matplotlib.patches import Circle, Rectangle

    n = len(entries)
    rows = math.ceil(n / cols)
    fig = Figure(figsize=(3.6 * cols, 3.7 * rows + 0.6))
    axes = fig.subplots(rows, cols, squeeze=False)
    for k in range(rows * cols):
        ax = axes[k // cols][k % cols]
        if k >= n:
            ax.axis("off")
            continue
        e = entries[k]
        pn = e.panel
        ax.add_patch(Rectangle((-pn.width / 2, -pn.height / 2), pn.width, pn.height, fill=False, ec="black", lw=1.2))
        x0, x1, y0, y1 = pn.inner
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, ec="gray", ls="--", lw=0.7))
        if e.status == "ok":
            for x, y in e.positions:
                ax.add_patch(Circle((x, y), pn.diameter / 2, fc=FAMILY_COLORS[e.family], alpha=0.35, ec="black", lw=0.6))
                ax.plot([x], [y], "k.", markersize=2)
            sub = f"min distance {e.meta.get('min_distance_mm', float('nan')):.1f} mm"
            if subtitles and e.id in subtitles:
                sub += "\n" + subtitles[e.id]
        else:
            ax.text(0, 0, "INFEASIBLE\n" + e.reason[:60], ha="center", va="center", fontsize=7, color="firebrick")
            sub = "not evaluated"
        m = 0.06 * max(pn.width, pn.height)
        ax.set_xlim(-pn.width / 2 - m, pn.width / 2 + m)
        ax.set_ylim(-pn.height / 2 - m, pn.height / 2 + m)
        ax.set_aspect("equal")
        ax.set_title(f"{e.id}\n{sub}", fontsize=8)
        ax.tick_params(labelsize=6)
    fig.suptitle(f"{title}\nto scale; driver disc = {Panel().diameter:g} mm (study hypothesis); solid = panel, dashed = center limits",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94), h_pad=2.5)
    fig.savefig(path, dpi=110, bbox_inches="tight")


def plot_lobe_maps(entries: list[Entry], path: Path, focus_name: str, band: str, title: str, cols: int = 4) -> None:
    from matplotlib.figure import Figure

    spec = REFERENCE_ROI
    roi_cache = {}
    n = len(entries)
    rows = math.ceil(n / cols)
    fig = Figure(figsize=(3.6 * cols + 1, 3.0 * rows + 0.8))
    axes = fig.subplots(rows, cols, squeeze=False)
    a = np.arange(-spec.half_az, spec.half_az + 1e-9, spec.step_deg)
    b = np.arange(-spec.half_el, spec.half_el + 1e-9, spec.step_deg)
    im = None
    for k in range(rows * cols):
        ax = axes[k // cols][k % cols]
        if k >= n:
            ax.axis("off")
            continue
        e = entries[k]
        lm = e.ref["cells"][(focus_name, band)].get("level_map_db")
        im = ax.imshow(np.array(lm, dtype=float), origin="lower", extent=(a[0] - 1, a[-1] + 1, b[0] - 1, b[-1] + 1),
                       cmap="magma", vmin=-30, vmax=0, aspect="equal")
        th = np.linspace(0, 2 * np.pi, 100)
        ax.plot(10 * np.cos(th), 10 * np.sin(th), "w--", lw=0.8)
        ax.plot(12 * np.cos(th), 12 * np.sin(th), "c:", lw=0.8)
        ax.set_title(f"{e.id}  ({e.n} el.)", fontsize=8)
        ax.tick_params(labelsize=6)
        ax.set_xlabel("azimuth offset (deg)", fontsize=7)
        ax.set_ylabel("elevation offset (deg)", fontsize=7)
    fig.colorbar(im, ax=axes, shrink=0.8, extend="min", label="compensated energy, dB re the focus")
    f = next(x for x in ALL_FOCI if x.name == focus_name)
    fig.suptitle(f"{title}\nfocus {focus_name}: {f.distance_mm:g} mm, az {f.az_deg:+g}, el {f.el_deg:+g}; {BAND_LABEL[band]}; "
                 "dashed = +/-10 deg cap, dotted = 12 deg exclusion (energy at ROI shell of radius = focus distance)", fontsize=8)
    fig.savefig(path, dpi=110, bbox_inches="tight")


def cut_levels(e: Entry, focus_mm, band: str, step_mm: float, df: float):
    """Energy levels (dB re A_REF) on the three planes through the focus; physics from field3d/broadband."""
    from .broadband import Spectrum, effective_level_db
    from .campaign import Params
    from .field3d import default_windows, energies, plane_points

    arr = to_array_config(e.id, e.positions)
    params = Params(gain=element_gain(e.n))
    lo, hi = BANDS[band]
    spec = Spectrum(lo, hi, df=df)
    windows = default_windows(focus_mm, "z")
    out = {}
    for name in ("XZ", "YZ", "XY"):
        u, v, pts = plane_points(name, focus_mm, windows, step_mm)
        _, e_f = energies(arr, tuple(c / 1000 for c in focus_mm), spec, pts, params)
        out[name] = (u, v, effective_level_db(e_f, A_REF, -60.0))
    return arr, out


def plot_cuts(entries: list[Entry], path: Path, band: str, df: float, step_mm: float, title: str) -> None:
    from matplotlib.figure import Figure

    from .field3d import PLANES, decorate_cut
    from .fieldplot import draw_level

    focus = tuple(float(v) for v in ALL_FOCI[0].position_mm)  # (0, 0, 800)
    n = len(entries)
    fig = Figure(figsize=(4.0 * n + 1.2, 12.3))
    gs = fig.add_gridspec(3, n + 1, width_ratios=[1] * n + [0.05], wspace=0.3, hspace=0.25, top=0.9, bottom=0.05)
    im = None
    roles = {"XZ": "beam axis + x (y = 0)", "YZ": "beam axis + y (x = 0)", "XY": "transverse (z = 800)"}
    for j, e in enumerate(entries):
        arr, cuts = cut_levels(e, focus, band, step_mm, df)
        for i, name in enumerate(("XZ", "YZ", "XY")):
            ax = fig.add_subplot(gs[i, j])
            u, v, lev = cuts[name]
            im = draw_level(ax, lev, u, v, -20.0, 20.0)
            decorate_cut(ax, arr, name, focus, u, v, (e.id + "\n" if i == 0 else "") + roles[name], fontsize=8)
    fig.colorbar(im, cax=fig.add_subplot(gs[:, n]), extend="both",
                 label="focused relative energy, dB re one A_REF speaker at 1 m (model); equal total power")
    fig.suptitle(f"{title}\nfocus (0, 0, 800) mm, {BAND_LABEL[band]}; white squares = speakers in the plane, gray triangles = projected",
                 fontsize=9, y=0.985)
    fig.savefig(path, dpi=100, bbox_inches="tight")


def plot_scaling(entries: list[Entry], path: Path) -> None:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(15, 4.6))
    ax1, ax2, ax3 = fig.subplots(1, 3)
    for env, ls in (("E520", "-"), ("E420", "--")):
        for fam in "RJVC":
            pts = {}
            for e in entries:
                if e.status == "ok" and e.envelope == env and e.family == fam and e.summary:
                    pts.setdefault(e.n, []).append(e.summary["b24"]["val"])
            ns = sorted(pts)
            if len(ns) < 2:
                continue
            for ax, key, lab in ((ax1, "F10_mean", "mean"), (ax2, "F10_worst", "worst"), (ax3, "psl_mean", "mean")):
                ax.plot(ns, [np.mean([p[key] for p in pts[n]]) for n in ns], ls, marker="o", color=FAMILY_COLORS[fam],
                        label=f"{fam} {env}")
    ax1.set_ylabel("share of ROI energy within +/-10 deg (validation foci, mean)")
    ax2.set_ylabel("same, worst validation focus")
    ax3.set_ylabel("contrast to secondary lobes, dB (validation foci, mean)")
    for ax in (ax1, ax2, ax3):
        ax.set_xlabel("number of elements")
        ax.set_xticks([8, 16, 32])
        ax.grid(alpha=0.3)
    ax1.legend(fontsize=7)
    fig.suptitle("Effect of the number of elements, 2-4 kHz, equal total power (solid: 520 mm panel, dashed: 420 mm panel)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=120, bbox_inches="tight")


def plot_regular_vs_irregular(entries: list[Entry], path: Path) -> None:
    from matplotlib.figure import Figure

    groups = [("E420", 8), ("E420", 16), ("E520", 8), ("E520", 16), ("E520", 32), ("E420", 32)]
    fig = Figure(figsize=(14, 4.8))
    ax = fig.subplots()
    width, x = 0.14, np.arange(len(groups))
    for k, fam in enumerate(("R", "J", "V", "C", "O")):
        vals = []
        for env, n in groups:
            sel = [e for e in entries if e.status == "ok" and e.summary and e.envelope == env and e.n == n and e.family == fam]
            vals.append(np.mean([e.summary["b24"]["val"]["F10_mean"] for e in sel]) if sel else np.nan)
        ax.bar(x + (k - 2) * width, vals, width, color=FAMILY_COLORS[fam], label=fam)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{n} el.\n{env[1:]} mm" for env, n in groups])
    ax.set_ylabel("share of ROI energy within +/-10 deg (validation foci, 2-4 kHz)")
    ax.legend(title="family (J: mean of seeds; missing = infeasible)", fontsize=8)
    ax.grid(axis="y", alpha=0.3)
    fig.suptitle("Regular (R) versus irregular (J, V, C) and optimized (O) layouts at matched element count and envelope")
    fig.savefig(path, dpi=120, bbox_inches="tight")


def plot_pareto(entries: list[Entry], flags: dict, path: Path) -> None:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(15, 5.2))
    ax1, ax2, ax3 = fig.subplots(1, 3)
    ok = [e for e in entries if e.status == "ok" and e.summary]
    for ax, xk, xl in ((ax1, "n", "number of elements"), (ax2, "area_mm2", "panel area (mm2)"),
                       (ax3, "robust_focus", "robustness: worst validation focus share")):
        for e in ok:
            o = objectives(e)
            nd = flags.get(e.id, False)
            ax.scatter(o[xk], o["selectivity"], s=40 + 40 * nd, c=FAMILY_COLORS[e.family], marker="o" if nd else "x",
                       edgecolors="black" if nd else None, linewidths=1.0, alpha=0.9)
            if nd:
                ax.annotate(e.id, (o[xk], o["selectivity"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(xl)
        ax.set_ylabel("selectivity: share within +/-10 deg (validation foci, mean of 3 bands)")
        ax.grid(alpha=0.3)
    from matplotlib.lines import Line2D
    ax1.legend(handles=[Line2D([], [], marker="o", ls="", color=c, label=f) for f, c in FAMILY_COLORS.items()]
               + [Line2D([], [], marker="x", ls="", color="gray", label="dominated")], fontsize=7)
    fig.suptitle("Trade-offs (4-objective non-dominated set circled: selectivity up, worst-focus share up, N down, area down)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(path, dpi=120, bbox_inches="tight")


# ---- tables ---------------------------------------------------------------------------------------------------------------

def write_tables(entries: list[Entry], outdir: Path, flags: dict) -> None:
    rows = []
    for e in entries:
        base = dict(id=e.id, family=e.family, n=e.n, envelope=e.envelope, panel_mm=f"{e.panel.width:g}x{e.panel.height:g}",
                    campaigns="".join(sorted(e.campaigns)), status=e.status, reason=e.reason, seed=e.seed,
                    min_distance_mm=e.meta.get("min_distance_mm"), reference_seconds=round(e.seconds, 2))
        if e.status != "ok":
            rows.append(base)
            continue
        o = objectives(e)
        base.update(selectivity=o["selectivity"], robust_focus=o["robust_focus"], robust_lines=o["robust_lines"],
                    pareto=flags.get(e.id, False))
        for band in BANDS:
            for name in ("train", "val"):
                s = e.summary[band][name]
                for k in ("F10_mean", "F10_worst", "psl_mean", "psl_worst", "F10_p10_lines_worst", "level_mean"):
                    base[f"{band}.{name}.{k}"] = s[k]
        rows.append(base)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with (outdir / "geometry-summary.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    (outdir / "geometry-summary.json").write_text(json.dumps(rows, indent=1, default=str) + "\n", encoding="utf-8")
    detail = [dict(id=e.id, band=b, focus=f, **{k: v for k, v in m.items() if not isinstance(v, (list, dict)) or k in ("top_lobes",)})
              for e in entries if e.status == "ok" and e.ref for (f, b), m in e.ref["cells"].items()]
    (outdir / "reference-cells.json").write_text(json.dumps(detail, default=str) + "\n", encoding="utf-8")


# ---- orchestration ---------------------------------------------------------------------------------------------------------

TERM_SECONDS = 3.1e-8  # measured cost of one (element x point x line x focus) term in the reference evaluation
ITERATION_SECONDS = 2.5e-3  # measured cost of one optimizer proposal (9 lines, 5 foci)


def estimate_seconds(profile: Profile) -> dict:
    """Rough planned cost, from the measured unit costs above (printed before a campaign starts)."""
    entries = [e for e in catalogue(profile) if e.status == "ok"]
    lines = sum(len(Spectrum(lo, hi, df=profile.df).lines()[0]) for lo, hi in BANDS.values())
    ref = sum(TERM_SECONDS * e.n * 5100 * len(ALL_FOCI) * lines for e in entries)
    jobs = len(profile.opt_ns) + (1 if 32 in profile.opt_ns else 0) + (2 if profile.wide_opt else 0)
    opt = jobs * len(profile.restarts) * profile.opt_iterations * ITERATION_SECONDS
    ref_opt = jobs * TERM_SECONDS * 20 * 5100 * len(ALL_FOCI) * lines
    return dict(geometries=len(entries), reference_seconds=ref + ref_opt, optimization_seconds=opt,
                total_minutes=(ref + ref_opt + opt + 120) / 60)


def export_finalists(entries: list[Entry], flags: dict, export_dir: Path) -> list[str]:
    ok = [e for e in entries if e.status == "ok" and e.summary]
    chosen = {e.id: e for e in ok if flags.get(e.id)}
    for e in ok:
        if e.family in ("O", "Ow"):
            chosen[e.id] = e
    for n in sorted({e.n for e in ok}):
        for env in sorted({e.envelope for e in ok}):
            pool = [e for e in ok if e.n == n and e.envelope == env and e.family in "RJVC"]
            if pool:
                best = max(pool, key=lambda e: objectives(e)["selectivity"])
                chosen[best.id] = best
    names = []
    for e in chosen.values():
        save_array_json(export_dir / f"{e.id}.json", e.id, e.positions, e.panel,
                        meta=dict(family=e.family, n=e.n, envelope=e.envelope, pareto=bool(flags.get(e.id)),
                                  selectivity=objectives(e)["selectivity"], note="reference diameter 65 mm is a hypothesis",
                                  **{k: v for k, v in e.meta.items() if isinstance(v, (int, float, str, bool))}))
        names.append(e.id)
    return sorted(names)


def run_campaign(profile: Profile, outdir: Path, export_dir: Path | None, log=print) -> list[Entry]:
    t_all = time.perf_counter()
    outdir.mkdir(parents=True, exist_ok=True)
    est = estimate_seconds(profile)
    log(f"profile {profile.name}: {est['geometries']} baseline geometries; estimated about {est['total_minutes']:.0f} minutes")
    cache = outdir / "reference"
    entries = catalogue(profile)
    for e in entries:
        if e.status == "ok":
            evaluate_entry(e, profile.df, cache, log)
        else:
            log(f"  {e.id}: infeasible ({e.reason})")
    opt_entries, restarts = run_optimization(profile, outdir, log)
    for e in opt_entries:
        evaluate_entry(e, profile.df, cache, log)
    entries += opt_entries
    flags = pareto_flags(entries)
    write_tables(entries, outdir, flags)
    log_data = dict(profile=profile.name, restarts=restarts, estimate=est,
                    total_seconds=time.perf_counter() - t_all,
                    reference_seconds={e.id: e.seconds for e in entries if e.status == "ok"})
    (outdir / "campaign-log.json").write_text(json.dumps(log_data, indent=1, default=str) + "\n", encoding="utf-8")
    finalists = export_finalists(entries, flags, export_dir) if export_dir else []
    plots = outdir / "plots"
    plots.mkdir(exist_ok=True)

    def sel(env, ns, fams="RJVCO"):
        return [e for e in entries if e.envelope == env and e.n in ns and e.family in tuple(fams)]

    order = {f: i for i, f in enumerate("RJVCO")}
    for tag, env in (("E420", "E420"), ("E520", "E520")):
        es = sorted(sel(env, profile.ns + (32,), "RJVCO"), key=lambda e: (e.n, order.get(e.family, 9), e.seed))
        plot_layouts(es, plots / f"layouts-{tag}.png", f"Layouts on the {env[1:]} x {env[1:]} mm panel")
    cs = sorted([e for e in entries if "C" in e.campaigns], key=lambda e: (e.n, e.family, e.panel.width))
    plot_layouts(cs, plots / "layouts-aperture.png", "Campaign C: fixed element count, growing aperture", cols=4)
    plot_scaling(entries, plots / "scaling-elements.png")
    plot_regular_vs_irregular(entries, plots / "regular-vs-irregular.png")
    plot_pareto(entries, flags, plots / "pareto.png")
    ok_by_id = {e.id: e for e in entries if e.status == "ok" and e.ref}
    for n, env in ((16, "E420"), (8, "E420"), (32, "E520")):
        pick = [i for i in (f"R{n}-{env}", f"J{n}-{env}-s1", f"V{n}-{env}", f"C{n}-{env}", f"O{n}-{env}", f"Ow{n}-{env}")
                if i in ok_by_id]
        if len(pick) >= 2:
            es = [ok_by_id[i] for i in pick]
            plot_lobe_maps(es, plots / f"lobes-{n}-{env}-T1.png", "T1", "b24", f"Angular lobe maps, {n} elements, {env}")
            plot_lobe_maps(es, plots / f"lobes-{n}-{env}-T4.png", "T4", "b24", f"Angular lobe maps, {n} elements, {env}")
            plot_cuts(es, plots / f"cuts-{n}-{env}.png", "b24", profile.cuts_df, 5.0, f"Orthogonal cuts, {n} elements, {env}")
    log(f"done in {time.perf_counter() - t_all:.0f} s; finalists exported: {len(finalists)}")
    return entries


# ---- figures for the optimized layouts ---------------------------------------------------------------------------

_ID = re.compile(r"^(R|J|V|C)(\d+)-(E\d+)(?:-s(\d+))?$")


def entry_from_id(entry_id: str, rows: dict, export_dir: Path) -> Entry:
    """Rebuild a baseline layout from its generator, or load an optimized one from its exported JSON."""
    from .array_config import load_array
    from .arraygen import positions_of

    r = rows[entry_id]
    if r["family"] in ("O", "Ow"):
        e = Entry(entry_id, r["family"], int(r["n"]), r["envelope"], {"D"})
        e.positions = positions_of(load_array(export_dir / f"{entry_id}.json"))
        e.meta = dict(min_distance_mm=check_positions(e.positions, e.panel)["min_distance_mm"])
        return e
    m = _ID.match(entry_id)
    return _make(m.group(1), int(m.group(2)), m.group(3), int(m.group(4) or 0), "")


def best_non_optimized(rows: dict, n: int, env: str, exclude=("R",)) -> str | None:
    pool = [i for i, r in rows.items() if r["status"] == "ok" and r["n"] == n and r["envelope"] == env
            and r["family"] in ("R", "J", "V", "C") and r["family"] not in exclude]
    return max(pool, key=lambda i: rows[i]["b24.val.F10_mean"]) if pool else None


def plot_generalization(rows: dict, path: Path) -> None:
    """Validation selectivity by band for the regular grid, the best non-optimized layout and the optimized ones."""
    from matplotlib.figure import Figure

    groups = [(8, "E420"), (16, "E420"), (32, "E520")]
    fig = Figure(figsize=(15, 4.8))
    axes = fig.subplots(1, 3)
    bands = list(BANDS)
    for ax, (n, env) in zip(axes, groups):
        ids = [f"R{n}-{env}"] if f"R{n}-{env}" in rows else []
        b = best_non_optimized(rows, n, env)
        if b:
            ids.append(b)
        for i in (f"O{n}-{env}", f"Ow{n}-{env}", f"O{n}-E420"):
            if i in rows and i not in ids:
                ids.append(i)
        if not ids:
            ax.axis("off")
            continue
        width = 0.8 / len(ids)
        for k, i in enumerate(ids):
            fam = rows[i]["family"]
            vals = [rows[i][f"{band}.val.F10_mean"] for band in bands]
            x = np.arange(len(bands)) + (k - (len(ids) - 1) / 2) * width
            ax.bar(x, vals, width, color=FAMILY_COLORS[fam], alpha=0.85 if fam != "R" else 0.6, label=i,
                   hatch="//" if i.endswith("E420") and env == "E520" else None, edgecolor="black", linewidth=0.5)
            if fam in ("O", "Ow"):  # training-foci value in the training band: no overfitting to the foci
                tb = "b18" if fam == "Ow" else "b24"
                ax.plot([x[bands.index(tb)]], [rows[i][f"{tb}.train.F10_mean"]], "kD", markersize=5)
        ax.set_xticks(range(len(bands)))
        ax.set_xticklabels([BAND_LABEL[b] + ("\n(training band of O)" if b == "b24" else "") for b in bands], fontsize=8)
        ax.set_title(f"{n} elements, {env[1:]} mm panel" + (" (+ O32 on 420 mm, hatched)" if n == 32 else ""), fontsize=9)
        ax.set_ylabel("share of ROI energy within +/-10 deg, validation foci")
        ax.grid(axis="y", alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("Optimized layouts versus the regular grid and the best non-optimized layout (diamond = training foci in the training band)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=120, bbox_inches="tight")


def plot_convergence(histories: dict, path: Path) -> None:
    from matplotlib.figure import Figure

    fig = Figure(figsize=(15, 4.4))
    axes = fig.subplots(1, len(histories))
    for ax, (title, hist) in zip(np.atleast_1d(axes), histories.items()):
        for label, h in hist.items():
            its = [a for a, _, _ in h]
            ax.plot(its, [c for _, _, c in h], marker=".", label=f"{label}: {h[-1][2]:.4f}")
        ax.set_title(title, fontsize=9)
        ax.set_xlabel("iteration")
        ax.set_ylabel("best objective J so far (training foci, lower is better)")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, title="restart: final J")
    fig.suptitle("Convergence of the restarts: very different starting layouts end at nearly the same objective", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(path, dpi=120, bbox_inches="tight")


def optimization_figures(outdir: Path, export_dir: Path, profile: Profile = PROFILES["full"], log=print) -> None:
    """Figures of the optimized layouts from a finished campaign (no re-evaluation; restarts are re-run for the curves)."""
    from .arrayopt import OptConfig, optimize

    rows = {r["id"]: r for r in json.loads((outdir / "geometry-summary.json").read_text())}
    plots = outdir / "plots"
    plots.mkdir(exist_ok=True)
    seq = []
    subtitles = {}
    for n, env in ((8, "E420"), (16, "E420"), (32, "E520")):
        ids = [f"R{n}-{env}"]
        b = best_non_optimized(rows, n, env)
        if b:
            ids.append(b)
        for i in (f"O{n}-{env}", f"Ow{n}-{env}", f"O{n}-E420"):
            if i in rows and i not in ids:
                ids.append(i)
        for i in ids[:4]:
            seq.append(entry_from_id(i, rows, export_dir))
            subtitles[i] = f"selectivity (validation, 2-4 kHz) {rows[i]['b24.val.F10_mean']:.3f}"
    plot_layouts(seq, plots / "layouts-optimized.png", "Optimized layouts next to the regular grid and the best non-optimized layout", cols=4,
                 subtitles=subtitles)
    plot_generalization(rows, plots / "generalization.png")
    logged = {(r["n"], r["envelope"], r["train_band"], r["init"], r["seed"]): r
              for r in json.loads((outdir / "campaign-log.json").read_text())["restarts"] if "j_best" in r}
    hist = {}
    for n, env in ((8, "E420"), (16, "E420"), (32, "E520")):
        cfg = OptConfig(iterations=profile.opt_iterations, band=BANDS["b24"], n_lines=9)
        h = {}
        for init, seed in profile.restarts:
            if (n, env, "b24", init, seed) not in logged:
                continue
            res = optimize(n, ENVELOPES[env], init, seed, cfg)
            ref = logged[(n, env, "b24", init, seed)]["j_best"]
            log(f"  {n} {env} {init}{seed}: J {res.j_best:.6f} (logged {ref:.6f}, difference {abs(res.j_best - ref):.1e})")
            h[f"{init}{seed}"] = [(0, res.j_start, res.j_start)] + res.history
        hist[f"{n} elements, {env[1:]} mm, 2-4 kHz"] = h
    plot_convergence(hist, plots / "convergence.png")


# ---- command line --------------------------------------------------------------------------------------------------------------

def _panel_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--panel", type=float, nargs=2, default=(420.0, 420.0), metavar=("W", "H"), help="panel size, mm")
    p.add_argument("--margin", type=float, default=32.5, help="minimum distance from a center to the panel edge, mm")
    p.add_argument("--min-distance", type=float, default=65.0, help="minimum distance between centers, mm")
    p.add_argument("--diameter", type=float, default=65.0, help="driver outer diameter: a study hypothesis, mm")


def _panel(a) -> Panel:
    return Panel(a.panel[0], a.panel[1], a.margin, a.min_distance, a.diameter)


def array_generate_main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="bean-lab array-generate", description="Generate a planar layout (R, J, V, C)")
    p.add_argument("--family", choices=("R", "J", "V", "C"), required=True)
    p.add_argument("--n", type=int, required=True)
    _panel_args(p)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--rows", type=int)
    p.add_argument("--cols", type=int)
    p.add_argument("--max-shift", type=float, help="J: maximum displacement from the grid, mm (default 0.25 pitch)")
    p.add_argument("--vogel-adjust", choices=("none", "scale"), default="scale")
    p.add_argument("--n-rings", type=int)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args(argv)
    try:
        panel = _panel(a)
        pos, meta = build(GeometrySpec(a.family, a.n, panel, a.seed, a.rows, a.cols, a.max_shift, a.vogel_adjust, a.n_rings))
        report = check_positions(pos, panel)
        save_array_json(a.out, a.out.stem, pos, panel, meta=dict(family=a.family, seed=a.seed, **{
            k: v for k, v in meta.items() if isinstance(v, (int, float, str, list))}))
    except InfeasibleGeometry as exc:
        p.error(f"infeasible: {exc}")
    except (OSError, ValueError) as exc:
        p.error(str(exc))
    print(f"Wrote {a.out}: {a.family}{a.n}, minimum center distance {report['min_distance_mm']:.2f} mm, "
          f"edge clearance {report['edge_clearance_mm']:.2f} mm; extent {np.ptp(pos, axis=0).round(1).tolist()} mm")
    return 0


def array_evaluate_main(argv: list[str]) -> int:
    from .array_config import load_array
    from .arraygen import positions_of

    p = argparse.ArgumentParser(prog="bean-lab array-evaluate", description="Reference-resolution evaluation of a layout file")
    p.add_argument("--array", type=Path, required=True)
    p.add_argument("--foci", choices=("train", "val", "all"), default="all")
    p.add_argument("--bands", nargs="+", choices=tuple(BANDS), default=list(BANDS))
    p.add_argument("--df", type=float, default=50.0, help="spectral line spacing, Hz")
    p.add_argument("--json", type=Path, help="write the per-focus metrics here")
    _panel_args(p)
    a = p.parse_args(argv)
    try:
        arr = load_array(a.array)
        pos = positions_of(arr)
        rep = check_positions(pos, _panel(a))
        if not rep["ok"]:
            raise InfeasibleGeometry("; ".join(rep["problems"]) + " (use --panel/--min-distance/--margin to match the layout)")
        foci = {"train": TRAIN_FOCI, "val": VALIDATION_FOCI, "all": ALL_FOCI}[a.foci]
        res = evaluate_reference(arr, foci, {b: BANDS[b] for b in a.bands}, df=a.df)
    except InfeasibleGeometry as exc:
        p.error(f"invalid layout: {exc}")
    except (OSError, ValueError) as exc:
        p.error(str(exc))
    names = [f.name for f in foci]
    print(f"{a.array.name}: {len(pos)} elements, min distance {rep['min_distance_mm']:.1f} mm; "
          f"gain per speaker {element_gain(len(pos)):.4f} (equal total power); {res['seconds']:.1f} s")
    print("band          level(dB)  F10 mean/worst   contrast mean/worst (dB)   line-p10 F10 worst")
    for b in a.bands:
        s = summarize(res, names, b)
        print(f"{BAND_LABEL[b]:<12} {s['level_mean']:8.2f}   {s['F10_mean']:.3f}/{s['F10_worst']:.3f}       "
              f"{s['psl_mean']:6.2f}/{s['psl_worst']:6.2f}            {s['F10_p10_lines_worst']:.3f}")
    if a.json:
        a.json.write_text(json.dumps({_cell_key(k): {kk: vv for kk, vv in v.items() if kk != "level_map_db"}
                                      for k, v in res["cells"].items()}, indent=1, default=str) + "\n", encoding="utf-8")
    return 0


def array_optimize_main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="bean-lab array-optimize", description="Optimize speaker positions (positions only)")
    p.add_argument("--n", type=int, required=True)
    _panel_args(p)
    p.add_argument("--init", nargs="+", default=["rand", "V"], choices=("rand", "R", "J", "V", "C"))
    p.add_argument("--seeds", nargs="+", type=int, default=[1])
    p.add_argument("--iterations", type=int, default=3000)
    p.add_argument("--band", choices=tuple(BANDS), default="b24", help="training band")
    p.add_argument("--lines", type=int, default=9, help="spectral lines used during the search")
    p.add_argument("--out", type=Path, required=True, help="output array JSON (best restart)")
    a = p.parse_args(argv)
    try:
        panel = _panel(a)
        cfg = OptConfig(iterations=a.iterations, band=BANDS[a.band], n_lines=a.lines)
        results, skipped = optimize_multi(a.n, panel, [(i, s) for i in a.init for s in a.seeds], cfg)
        if not results:
            raise InfeasibleGeometry("no initial layout could be built: " + "; ".join(x[2] for x in skipped))
        best = min(results, key=lambda r: r.j_best)
        save_array_json(a.out, a.out.stem, best.positions, panel, meta=dict(
            family="O", init=best.init, seed=best.seed, j_best=best.j_best, train_band=a.band, iterations=a.iterations,
            components=best.components_best, evaluations=best.evaluations, note="optimized on the training foci only"))
    except InfeasibleGeometry as exc:
        p.error(f"infeasible: {exc}")
    except (OSError, ValueError) as exc:
        p.error(str(exc))
    print(f"Wrote {a.out}: best of {len(results)} restarts, J {best.j_best:.4f} (init {best.init}, seed {best.seed}); "
          f"components {{{', '.join(f'{k}={v:.3f}' for k, v in best.components_best.items())}}}")
    return 0


def array_campaign_main(argv: list[str]) -> int:
    p = argparse.ArgumentParser(prog="bean-lab array-campaign", description="Layout campaigns A-D and Pareto analysis")
    p.add_argument("--profile", choices=tuple(PROFILES), default="quick")
    p.add_argument("--outdir", type=Path, default=Path("results/arrayopt"))
    p.add_argument("--export-dir", type=Path, help="write finalist layouts (array JSON) here")
    p.add_argument("--estimate", action="store_true", help="print the estimated cost and exit")
    p.add_argument("--optimization-figures", action="store_true",
                   help="only redraw the optimized-layout figures from a finished campaign in --outdir (uses --export-dir)")
    a = p.parse_args(argv)
    profile = PROFILES[a.profile]
    est = estimate_seconds(profile)
    print(f"profile {profile.name}: {est['geometries']} baseline geometries, reference spacing {profile.df:g} Hz, "
          f"{profile.opt_iterations} iterations x {len(profile.restarts)} restarts; estimated about {est['total_minutes']:.0f} minutes")
    if a.estimate:
        return 0
    if a.optimization_figures:
        try:
            optimization_figures(a.outdir, a.export_dir or Path("configs/arrays/optimized"), profile)
        except (OSError, ValueError, KeyError) as exc:
            p.error(str(exc))
        return 0
    try:
        entries = run_campaign(profile, a.outdir, a.export_dir)
    except (OSError, ValueError, ImportError) as exc:
        p.error(str(exc))
    ok = [e for e in entries if e.status == "ok"]
    print(f"{len(ok)} geometries evaluated, {len(entries) - len(ok)} infeasible; results in {a.outdir}")
    return 0
