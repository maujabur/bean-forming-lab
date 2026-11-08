"""Position optimizer: stochastic search with single-element moves (simulated annealing + greedy polish).

The only variables are the speaker positions (2N coordinates in the panel plane). Amplitudes and filters
are never touched; total emitted power is fixed by arrayeval.element_gain.

Why this method. The feasible set (inside the panel, centers >= min_distance apart) is non-convex and,
for 32 drivers on a 420 mm panel, quite dense; a move of one element is accepted only if it keeps every
constraint, so each state visited is a valid layout (no repair, no penalty term). Moving one element
changes one contribution of the cached complex field, so a proposal costs a few milliseconds whatever
N is (arrayeval.FastEvaluator). The search is stochastic and restarted from different initial layouts
(random relaxed, Vogel, rings, jittered grid): it finds good layouts, never a certified global optimum.
"""

from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np

from .arraygen import (InfeasibleGeometry, GeometrySpec, Panel, build, check_positions, relaxed_random,
                       validate_positions)
from .arrayeval import FAST_ROI, TRAIN_FOCI, RoiSpec, Weights, fast_setup

SLACK_MM = 1e-6  # proposals keep this much extra clearance so that file round trips stay valid


@dataclass(frozen=True)
class OptConfig:
    iterations: int = 3000
    step_mm: float = 25.0  # initial std of a local move
    step_min_mm: float = 2.0
    t0: float = 0.02  # initial temperature, in units of the objective
    t1: float = 2e-4
    jump_prob: float = 0.03  # probability of a global jump to a random position
    polish_fraction: float = 0.15  # final fraction of iterations with T = 0 and small steps
    refresh_every: int = 500
    band: tuple[float, float] = (2000.0, 4000.0)
    n_lines: int = 9
    weights: Weights = Weights()
    spec: RoiSpec = FAST_ROI
    foci: tuple = TRAIN_FOCI

    def __post_init__(self) -> None:
        if self.iterations < 1 or self.step_mm <= 0 or self.step_min_mm <= 0 or not 0 <= self.polish_fraction < 1:
            raise ValueError("invalid optimizer configuration")


@dataclass
class OptResult:
    positions: np.ndarray
    init: str
    seed: int
    j_start: float
    j_best: float
    components_start: dict
    components_best: dict
    history: list = field(default_factory=list)  # (iteration, J current, J best)
    evaluations: int = 0
    accepted: int = 0
    infeasible_proposals: int = 0
    seconds: float = 0.0
    memory_bytes: int = 0


def initial_positions(init: str, n: int, panel: Panel, seed: int) -> np.ndarray:
    """Valid starting layout; raises InfeasibleGeometry when the family cannot be built (no silent repair)."""
    if init == "rand":
        return relaxed_random(n, panel, seed)
    if init in ("R", "J", "V", "C"):
        return build(GeometrySpec(init, n, panel, seed=seed))[0]
    raise ValueError(f"unknown initializer {init!r}")


def _feasible(pos: np.ndarray, i: int, new: np.ndarray, panel: Panel) -> bool:
    x0, x1, y0, y1 = panel.inner
    if not (x0 <= new[0] <= x1 and y0 <= new[1] <= y1):
        return False
    d = np.linalg.norm(pos - new, axis=1)
    d[i] = np.inf
    return bool(d.min() >= panel.min_distance + SLACK_MM)


def optimize(n: int, panel: Panel, init: str, seed: int, cfg: OptConfig = OptConfig()) -> OptResult:
    rng = np.random.default_rng(seed)
    t_start = time.perf_counter()
    pos = initial_positions(init, n, panel, seed)
    ev = fast_setup(pos, cfg.foci, cfg.band, cfg.n_lines, cfg.spec, cfg.weights)
    cur = ev.objective()
    j_cur = cur["J"]
    best_pos, best = pos.copy(), cur
    res = OptResult(pos.copy(), init, seed, j_cur, j_cur, {k: v for k, v in cur.items() if k.startswith("J")},
                    {k: v for k, v in cur.items() if k.startswith("J")})
    x0, x1, y0, y1 = panel.inner
    polish_start = int(cfg.iterations * (1 - cfg.polish_fraction))
    hist_every = max(cfg.iterations // 20, 1)
    for it in range(cfg.iterations):
        if it and it % hist_every == 0:  # recorded before the proposal, so a rejected proposal never skips a checkpoint
            res.history.append((it, float(j_cur), float(best["J"])))
        frac = it / max(cfg.iterations - 1, 1)
        in_polish = it >= polish_start
        temp = 0.0 if in_polish else cfg.t0 * (cfg.t1 / cfg.t0) ** (it / max(polish_start, 1))
        sigma = cfg.step_mm * (cfg.step_min_mm / cfg.step_mm) ** frac
        i = int(rng.integers(n))
        if not in_polish and rng.random() < cfg.jump_prob:
            new = np.array([rng.uniform(x0, x1), rng.uniform(y0, y1)])
        else:
            new = ev.pos[i, :2] + rng.normal(0.0, sigma, 2)
        if not _feasible(ev.pos[:, :2], i, new, panel):
            res.infeasible_proposals += 1
            continue
        cand, totals, new3 = ev.try_move(i, new)
        dj = cand["J"] - j_cur
        if dj <= 0 or (temp > 0 and rng.random() < math.exp(-dj / temp)):
            ev.commit(i, totals, new3)
            j_cur = cand["J"]
            res.accepted += 1
            if j_cur < best["J"] - 1e-12:
                best, best_pos = cand, ev.pos[:, :2].copy()
        if (it + 1) % cfg.refresh_every == 0:
            ev.refresh()
            j_cur = ev.objective()["J"]
    res.history.append((cfg.iterations, float(j_cur), float(best["J"])))
    # the best state is re-evaluated from scratch (no accumulated rounding) before it is reported
    final = fast_setup(best_pos, cfg.foci, cfg.band, cfg.n_lines, cfg.spec, cfg.weights)
    comp = final.objective()
    validate_positions(best_pos, panel)
    res.positions = best_pos
    res.j_best = float(comp["J"])
    res.components_best = {k: v for k, v in comp.items() if k.startswith("J")}
    res.evaluations = ev.evaluations + final.evaluations
    res.seconds = time.perf_counter() - t_start
    res.memory_bytes = int(sum(t.nbytes for t in ev.totals) * 2 + sum(r.pts_mm.nbytes * 3 for r in ev.rois))
    return res


def optimize_multi(n: int, panel: Panel, inits: list[tuple[str, int]], cfg: OptConfig = OptConfig(), log=print):
    """Run several restarts; layouts that cannot be built are reported and skipped. Returns (results, skipped)."""
    results, skipped = [], []
    for init, seed in inits:
        try:
            r = optimize(n, panel, init, seed, cfg)
        except InfeasibleGeometry as exc:
            skipped.append((init, seed, str(exc)))
            log(f"  init {init} seed {seed}: infeasible ({exc})")
            continue
        log(f"  init {init} seed {seed}: J {r.j_start:.4f} -> {r.j_best:.4f} "
            f"({r.evaluations} evals, {r.seconds:.1f} s, accepted {r.accepted})")
        results.append(r)
    return results, skipped
