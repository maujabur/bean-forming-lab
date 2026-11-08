import json
import math
from pathlib import Path

import numpy as np
import pytest

import bean_lab.arrayopt as arrayopt
from bean_lab.array_config import load_array
from bean_lab.arraycampaign import (ENVELOPES, PROFILES, Entry, catalogue, estimate_seconds, evaluate_entry, export_finalists,
                                    objectives, pareto_flags)
from bean_lab.arrayeval import (A_REF, ALL_FOCI, BANDS, FAST_ROI, N_REF, REFERENCE_ROI, TRAIN_FOCI, VALIDATION_FOCI, Focus,
                                FocusRoi, Weights, core_metrics, element_gain, evaluate_reference, fast_setup, map_metrics,
                                objective_components, pareto_mask, reference_energies, summarize, _run_width)
from bean_lab.arraygen import (GOLDEN_ANGLE, GeometrySpec, InfeasibleGeometry, Panel, build, check_positions, jittered_grid,
                               positions_of, regular_grid, relaxed_random, rings, rings_auto, save_array_json,
                               to_array_config, validate_positions, vogel)
from bean_lab.arrayopt import OptConfig, _feasible, initial_positions, optimize, optimize_multi
from bean_lab.broadband import Spectrum
from bean_lab.cli import main

P420, P520 = ENVELOPES["E420"], ENVELOPES["E520"]


# ---- panel and constraints ------------------------------------------------------------------------------------------

def test_panel_geometry_and_validation():
    assert P420.inner == (-177.5, 177.5, -177.5, 177.5) and P420.span == (355.0, 355.0)
    for bad in (dict(width=0), dict(margin=-1), dict(min_distance=0), dict(width=float("nan")), dict(width=60, height=60, margin=30)):
        with pytest.raises(ValueError):
            Panel(**bad)


def test_constraint_report_detects_each_violation():
    ok = np.array([[0.0, 0.0], [65.0, 0.0]])
    assert check_positions(ok, P420)["ok"] and check_positions(ok, P420)["min_distance_mm"] == pytest.approx(65.0)
    close = check_positions(np.array([[0.0, 0.0], [64.9, 0.0]]), P420)
    assert not close["ok"] and "minimum center distance" in close["problems"][0]
    edge = check_positions(np.array([[178.0, 0.0], [0.0, 0.0]]), P420)  # 32 mm from the edge < 32.5 mm margin
    assert not edge["ok"] and "edge clearance" in edge["problems"][0]
    with pytest.raises(InfeasibleGeometry):
        validate_positions(np.array([[0.0, 0.0], [10.0, 0.0]]), P420)
    for bad in (np.zeros((3,)), np.array([[np.nan, 0.0]])):
        with pytest.raises(ValueError):
            check_positions(bad, P420)


# ---- generators -----------------------------------------------------------------------------------------------------------

def test_regular_grid_counts_positions_and_feasibility():
    g = regular_grid(2, 4, P420)
    assert g.shape == (8, 2)
    xs = np.linspace(-177.5, 177.5, 4)
    assert np.allclose(sorted(set(np.round(g[:, 0], 9))), xs) and set(np.round(g[:, 1], 9)) == {-177.5, 177.5}
    assert regular_grid(4, 4, P420).shape == (16, 2)
    with pytest.raises(InfeasibleGeometry, match="minimum center distance 50.7"):
        regular_grid(4, 8, P420)  # 355 / 7 = 50.7 mm < 65 mm: reported, never shrunk
    g32 = regular_grid(4, 8, P520)  # 455 mm span / 7 = 65 mm exactly
    assert g32.shape == (32, 2) and check_positions(g32, P520)["min_distance_mm"] == pytest.approx(65.0)


def test_every_generator_is_reproducible_and_valid():
    for fam in "RJVC":
        a, ma = build(GeometrySpec(fam, 16, P420, seed=3))
        b, mb = build(GeometrySpec(fam, 16, P420, seed=3))
        assert np.array_equal(a, b) and len(a) == 16 and check_positions(a, P420)["ok"]
    j1, _ = build(GeometrySpec("J", 8, P420, seed=1))
    j2, _ = build(GeometrySpec("J", 8, P420, seed=2))
    assert not np.array_equal(j1, j2)
    r1, r2 = relaxed_random(16, P420, 5), relaxed_random(16, P420, 5)
    assert np.array_equal(r1, r2) and check_positions(r1, P420)["ok"] and not np.array_equal(r1, relaxed_random(16, P420, 6))


def test_jittered_grid_limits():
    base = regular_grid(4, 4, P420)
    assert np.array_equal(jittered_grid(4, 4, P420, 0.0, 1), base)  # zero shift = the regular grid
    shift = 25.0
    j = jittered_grid(4, 4, P420, shift, 7)
    assert np.all(np.linalg.norm(j - base, axis=1) <= shift + 1e-9) and check_positions(j, P420)["ok"]
    with pytest.raises(InfeasibleGeometry):
        jittered_grid(4, 8, P420, 10.0, 1)  # the base grid itself is invalid: no repair
    with pytest.raises(ValueError):
        jittered_grid(4, 4, P420, -1.0, 1)


def test_vogel_formula_and_documented_adjustment():
    pos, meta = vogel(16, P420, adjust="none")
    i = np.arange(16)
    r = 177.5 * np.sqrt((i + 0.5) / 16)
    assert np.allclose(np.hypot(pos[:, 0], pos[:, 1]), r)
    assert np.allclose(np.arctan2(pos[:, 1], pos[:, 0]), np.angle(np.exp(1j * i * GOLDEN_ANGLE)))
    assert meta["scale"] == 1.0
    # 32 elements violate the minimum distance as generated: rejected, or scaled by the documented factor
    with pytest.raises(InfeasibleGeometry):
        vogel(32, P420, adjust="none")
    big, meta = vogel(32, P520, adjust="scale")
    assert meta["scale"] > 1 and check_positions(big, P520)["min_distance_mm"] == pytest.approx(65.0, abs=1e-6)
    with pytest.raises(InfeasibleGeometry):
        vogel(32, P420, adjust="scale")  # the scaled spiral no longer fits the panel
    with pytest.raises(ValueError):
        vogel(8, P420, adjust="repair")


def test_rings_counts_radii_and_infeasibility():
    pos = rings([3, 5], [88.75, 177.5], P420)
    assert len(pos) == 8 and np.allclose(sorted(np.hypot(pos[:, 0], pos[:, 1]))[:3], 88.75)
    auto, meta = rings_auto(16, P420)
    assert sum(meta["counts"]) == 16 and meta["n_rings"] in (2, 3, 4) and len(auto) == 16
    with pytest.raises(InfeasibleGeometry, match="no ring structure"):
        rings_auto(32, P420)
    with pytest.raises(ValueError):
        rings([3], [10.0, 20.0], P420)


def test_infeasible_spec_is_reported_for_32_elements_on_420():
    for fam in "RJVC":
        with pytest.raises(InfeasibleGeometry):
            build(GeometrySpec(fam, 32, P420))


def test_json_round_trip_is_exact_and_loadable(tmp_path):
    pos, _ = build(GeometrySpec("V", 16, P420))
    path = tmp_path / "v16.json"
    save_array_json(path, "v16", pos, P420, meta={"family": "V"})
    arr = load_array(path)
    back = positions_of(arr)
    assert np.allclose(back, pos, atol=1e-9) and check_positions(back, P420)["ok"]
    obj = json.loads(path.read_text())
    assert obj["panel"]["diameter"] == 65.0 and "hypothesis" in obj["panel"]["note"]
    from bean_lab.array_config import ArrayConfig, Element
    with pytest.raises(ValueError):  # a layout out of the z = 0 plane is not a panel layout
        positions_of(ArrayConfig("tilted", (Element(1, (0.0, 0.0, 0.1)),)))


# ---- foci, power normalization, ROI ----------------------------------------------------------------------------------------------

def test_foci_sets_are_disjoint_and_cover_the_mission_grid():
    assert not ({(f.distance_mm, f.az_deg, f.el_deg) for f in TRAIN_FOCI} & {(f.distance_mm, f.az_deg, f.el_deg) for f in VALIDATION_FOCI})
    assert not ({f.name for f in TRAIN_FOCI} & {f.name for f in VALIDATION_FOCI})
    for foci in (TRAIN_FOCI, VALIDATION_FOCI):
        assert {f.distance_mm for f in foci} == {600.0, 800.0, 1200.0, 1500.0}
    assert {abs(f.az_deg) for f in ALL_FOCI} == {0, 15, 30} and {f.el_deg for f in ALL_FOCI} == {-15, 0, 15}
    for f in ALL_FOCI:
        assert np.linalg.norm(f.position_mm) == pytest.approx(f.distance_mm) and np.linalg.norm(f.direction) == pytest.approx(1)
    f = Focus("x", "val", 800, 30, 0)
    assert f.position_mm == pytest.approx([400.0, 0.0, 800 * math.cos(math.radians(30))])
    assert Focus("x", "val", 1000, 0, 15).position_mm[1] == pytest.approx(1000 * math.sin(math.radians(15)))


def test_total_emitted_power_is_constant_across_element_counts():
    for n in (1, 8, 16, 32, 64):
        assert n * element_gain(n) ** 2 / 2 == pytest.approx(N_REF * A_REF ** 2 / 2)
    assert element_gain(N_REF) == A_REF
    with pytest.raises(ValueError):
        element_gain(0)


def test_coherent_focus_level_follows_the_power_normalization():
    res = {}
    for n in (8, 16):
        pos, _ = build(GeometrySpec("V", n, P420))
        res[n] = evaluate_reference(to_array_config("t", pos), foci=[TRAIN_FOCI[0]], bands={"b24": BANDS["b24"]}, df=500)
    for n, r in res.items():
        c = r["cells"][("T1", "b24")]
        assert c["level_focus_db"] == pytest.approx(c["level_coherent_db"], abs=1e-9)  # delays align every frequency
        assert c["level_raw_db"] == pytest.approx(c["level_focus_db"] - 10 * math.log10(n / 8))
    # at equal total power the coherent focus of 16 speakers is about 10 log10(2) = 3 dB above 8 speakers
    diff = res[16]["cells"][("T1", "b24")]["level_focus_db"] - res[8]["cells"][("T1", "b24")]["level_focus_db"]
    assert diff == pytest.approx(3.01, abs=0.15)


def test_roi_structure():
    f = Focus("t", "train", 800, 15, 0)
    roi = FocusRoi(f, FAST_ROI)
    assert np.allclose(roi.pts_mm[0], f.position_mm) and roi.w[0] == 0 and roi.gamma[0] == 0
    r = np.linalg.norm(roi.pts_mm[1:], axis=1) / 800
    assert set(np.round(r, 9)) == {0.85, 1.0, 1.15}
    assert np.all(roi.pts_mm[1:, 2] / np.linalg.norm(roi.pts_mm[1:], axis=1) >= math.cos(math.radians(85)) - 1e-12)
    assert np.all(roi.w[1:] > 0) and roi.gamma[1:].max() <= 75  # great-circle angle within the +/-60 x +/-40 window
    i0 = roi.mid_grid[roi.ib0, roi.ia0]
    assert np.allclose(roi.pts_mm[i0], f.position_mm) and roi.gamma[i0] == pytest.approx(0)
    # independent geometry: great-circle angle to the focus direction and solid-angle weight of the middle shell
    d = f.direction
    e_h = np.cross([0.0, 1.0, 0.0], d)
    e_h /= np.linalg.norm(e_h)
    e_v = np.cross(d, e_h)
    mid = roi.shell == 1
    u = roi.pts_mm[mid] / np.linalg.norm(roi.pts_mm[mid], axis=1)[:, None]
    assert np.allclose(roi.gamma[mid], np.degrees(np.arccos(np.clip(u @ d, -1, 1))), atol=1e-9)
    assert np.allclose(roi.w[mid], np.cos(np.arcsin(np.clip(u @ e_v, -1, 1))), atol=1e-12)
    outer = roi.shell == 0
    uo = roi.pts_mm[outer] / np.linalg.norm(roi.pts_mm[outer], axis=1)[:, None]
    assert np.allclose(roi.w[outer], 4 * np.cos(np.arcsin(np.clip(uo @ e_v, -1, 1))), atol=1e-12)  # stride 2: 4x area
    fine = FocusRoi(f, REFERENCE_ROI)
    assert fine.n > 3 * roi.n and len(fine.shell[1:]) > 0 and set(fine.shell[1:]) == {0, 1, 2, 3, 4}


def test_isotropic_source_has_analytic_shares_and_censored_regions():
    arr = to_array_config("one", np.zeros((1, 2)))  # one speaker at the panel center
    f = TRAIN_FOCI[0]
    roi = FocusRoi(f, REFERENCE_ROI)
    freqs, w = Spectrum(2000, 4000, n=5).lines()
    m = core_metrics(reference_energies(arr, f, roi, freqs), w, roi)
    wt = roi.w[1:]
    for win in (5.0, 10.0, 15.0):  # compensated energy is uniform: the share is the solid-angle weight ratio
        assert m[f"F{win:g}"] == pytest.approx(wt[roi.gamma[1:] <= win].sum() / wt.sum(), rel=1e-9)
    assert m["psl_cap"] == pytest.approx(0.0, abs=1e-9)  # nothing stands out
    mm = map_metrics(roi, m["ec"], f.distance_mm)
    assert mm["area_3_censored"] and mm["area_6_censored"] and mm["width_h_censored"] and mm["width_v_censored"]


def test_run_width_interpolation_and_censoring():
    lv = np.array([-10.0, -6.0, -2.0, 0.0, -2.0, -6.0, -10.0])
    valid = np.ones(7, dtype=bool)
    w, cens = _run_width(lv, valid, 3, -3.0, 1.0)
    assert w == pytest.approx(2.5) and not cens  # linear crossings of -3 dB at +-1.25 samples
    w2, cens2 = _run_width(lv, valid, 3, -12.0, 1.0)
    assert cens2 and w2 == pytest.approx(6.0)
    v2 = valid.copy()
    v2[5:] = False
    assert _run_width(lv, v2, 3, -12.0, 1.0)[1]  # ends at the domain boundary, not at a drop


# ---- evaluator coherence --------------------------------------------------------------------------------------------------------------

def test_fast_path_equals_reference_path_incremental_and_full():
    pos, _ = build(GeometrySpec("V", 16, P420))
    arr = to_array_config("t", pos)
    ev = fast_setup(pos, TRAIN_FOCI, (2000, 4000), 9)
    freqs, w = Spectrum(2000, 4000, n=9).lines()
    ref = objective_components([core_metrics(reference_energies(arr, f, roi, freqs), w, roi)
                                for f, roi in zip(TRAIN_FOCI, ev.rois)], Weights())
    got = ev.objective()
    for k in ("J", "J_out", "J_psl", "J_worst", "J_spec"):
        assert got[k] == pytest.approx(ref[k], abs=1e-12), k  # real delays vs the phase-only update rule
    new = pos[3] + np.array([12.0, -7.0])
    cand, totals, new3 = ev.try_move(3, new)
    pos2 = pos.copy()
    pos2[3] = new
    assert cand["J"] == pytest.approx(fast_setup(pos2, TRAIN_FOCI, (2000, 4000), 9).objective()["J"], abs=1e-12)
    ev.commit(3, totals, new3)
    ev.refresh()
    assert ev.objective()["J"] == pytest.approx(cand["J"], abs=1e-12)


def test_blocked_evaluation_equals_unblocked():
    pos, _ = build(GeometrySpec("R", 8, P420))
    arr = to_array_config("t", pos)
    f = ALL_FOCI[2]
    roi = FocusRoi(f, FAST_ROI)
    freqs = np.array([2500.0, 3500.0])
    assert np.array_equal(reference_energies(arr, f, roi, freqs, block=100), reference_energies(arr, f, roi, freqs, block=10 ** 6))


def test_fast_and_reference_resolution_agree_within_sampling_error():
    pos, _ = build(GeometrySpec("V", 16, P420))
    arr = to_array_config("t", pos)
    freqs, w = Spectrum(2000, 4000, n=9).lines()
    for f in TRAIN_FOCI[:3]:
        a = core_metrics(reference_energies(arr, f, FocusRoi(f, FAST_ROI), freqs), w, FocusRoi(f, FAST_ROI))
        b = core_metrics(reference_energies(arr, f, FocusRoi(f, REFERENCE_ROI), freqs), w, FocusRoi(f, REFERENCE_ROI))
        assert abs(a["F10"] - b["F10"]) < 0.04  # 4 deg versus 2 deg sampling of the same quantity
        assert abs(a["psl_cap"] - b["psl_cap"]) < 1.5


def test_mirror_invariance_of_the_metrics():
    pos, _ = build(GeometrySpec("J", 16, P420, seed=4))
    mirrored = pos * np.array([-1.0, 1.0])  # x -> -x
    f, fm = Focus("a", "train", 800, 15, 10), Focus("b", "train", 800, -15, 10)
    freqs, w = Spectrum(2000, 4000, n=5).lines()
    out = []
    for p, focus in ((pos, f), (mirrored, fm)):
        roi = FocusRoi(focus, FAST_ROI)
        out.append(core_metrics(reference_energies(to_array_config("t", p), focus, roi, freqs), w, roi))
    for k in ("F5", "F10", "F15", "psl_cap", "e_focus"):
        assert out[0][k] == pytest.approx(out[1][k], rel=1e-9, abs=1e-12)


def test_objective_components_are_normalized_and_use_the_worst_case():
    good = dict(F10=0.9, psl_cap=20.0, F10_lines=np.array([0.85, 0.9, 0.95]))
    bad = dict(F10=0.2, psl_cap=-3.0, F10_lines=np.array([0.1, 0.2, 0.3]))
    c = objective_components([good, bad], Weights())
    assert all(0 <= c[k] <= 1 for k in ("J_out", "J_psl", "J_worst", "J_spec"))
    assert c["J_worst"] == pytest.approx(0.8) and c["J_out"] == pytest.approx(0.45)  # worst case >= mean
    assert c["J_psl"] == pytest.approx((0 + 1) / 2)
    only_worst = objective_components([good, bad], Weights(out=0, psl=0, worst=1, spec=0))
    assert only_worst["J"] == pytest.approx(0.8)
    with pytest.raises(ValueError):
        Weights(out=-1)


def test_summarize_and_pareto_mask():
    cells = {("A", "b24"): dict(F10=0.5, psl_cap=10.0, F10_p10_lines=0.4, level_focus_db=20.0),
             ("B", "b24"): dict(F10=0.3, psl_cap=6.0, F10_p10_lines=0.2, level_focus_db=18.0)}
    s = summarize(dict(cells=cells), ["A", "B"], "b24")
    assert s["F10_mean"] == pytest.approx(0.4) and s["F10_worst"] == pytest.approx(0.3) and s["psl_worst"] == pytest.approx(6.0)
    pts = np.array([[0.5, 8, 100.0], [0.4, 8, 100.0], [0.6, 16, 100.0], [0.5, 8, 90.0], [0.5, 8, 100.0]])
    mask = pareto_mask(pts, [True, False, False])  # max selectivity, min N, min area
    assert mask.tolist() == [False, False, True, True, False]
    assert pareto_mask(np.array([[1.0, 1.0], [1.0, 1.0]]), [True, True]).tolist() == [True, True]  # ties do not dominate
    assert pareto_mask(np.array([[1.0, 2.0]]), [True, False]).tolist() == [True]


# ---- optimizer -----------------------------------------------------------------------------------------------------------------------------

def test_feasible_move_checks_box_and_distance():
    pos = np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
    assert _feasible(pos, 0, np.array([10.0, 10.0]), P420)
    assert not _feasible(pos, 0, np.array([40.0, 0.0]), P420)  # 60 mm from element 1
    assert not _feasible(pos, 0, np.array([200.0, 0.0]), P420)  # outside the center limits
    assert _feasible(pos, 0, np.array([0.0, 0.0]), P420)  # an element does not conflict with itself


def test_optimizer_is_reproducible_valid_and_never_worse():
    cfg = OptConfig(iterations=200)
    a = optimize(8, P420, "V", 0, cfg)
    b = optimize(8, P420, "V", 0, cfg)
    assert np.array_equal(a.positions, b.positions) and a.j_best == b.j_best
    assert check_positions(a.positions, P420)["ok"] and a.j_best <= a.j_start + 1e-12 and a.evaluations > 100
    c = optimize(8, P420, "rand", 1, cfg)
    d = optimize(8, P420, "rand", 2, cfg)
    assert not np.array_equal(c.positions, d.positions)
    assert a.memory_bytes > 0 and a.seconds > 0 and len(a.history) > 5


def test_optimizer_uses_training_foci_only(monkeypatch):
    seen = []
    real = arrayopt.fast_setup

    def spy(positions, foci, *args, **kwargs):
        seen.append([f.name for f in foci])
        return real(positions, foci, *args, **kwargs)

    monkeypatch.setattr(arrayopt, "fast_setup", spy)
    optimize(8, P420, "V", 0, OptConfig(iterations=20))
    assert seen and all(set(names) <= {f.name for f in TRAIN_FOCI} for names in seen)
    assert OptConfig().foci == TRAIN_FOCI


def test_optimization_improves_the_training_objective_and_keeps_every_constraint():
    cfg = OptConfig(iterations=600)
    r = optimize(8, P420, "rand", 3, cfg)
    assert r.j_best < r.j_start - 0.02
    assert check_positions(r.positions, P420)["min_distance_mm"] >= 65.0


def test_initializers_and_infeasible_restarts_are_reported():
    assert len(initial_positions("rand", 32, P420, 1)) == 32
    with pytest.raises(InfeasibleGeometry):
        initial_positions("V", 32, P420, 0)
    with pytest.raises(ValueError):
        initial_positions("hex", 8, P420, 0)
    logs = []
    results, skipped = optimize_multi(32, P420, [("V", 0), ("rand", 1)], OptConfig(iterations=30), logs.append)
    assert len(results) == 1 and results[0].init == "rand" and skipped[0][0] == "V" and "infeasible" in logs[0]
    with pytest.raises(ValueError):
        OptConfig(iterations=0)


# ---- campaign bookkeeping ---------------------------------------------------------------------------------------------------------------------

def test_catalogue_records_infeasible_combinations_and_never_evaluates_them():
    entries = {e.id: e for e in catalogue(PROFILES["full"])}
    for fam in "RV":
        assert entries[f"{fam}32-E420"].status == "infeasible" and entries[f"{fam}32-E420"].positions is None
        assert entries[f"{fam}32-E420"].reason
    assert entries["R32-E520"].status == "ok" and entries["C32-E520"].status == "ok"
    assert entries["R16-E420"].campaigns >= {"A", "C"} and "B" in entries["R16-E520"].campaigns
    assert {k for k in entries if k.startswith("J8-E420")} == {"J8-E420-s1", "J8-E420-s2", "J8-E420-s3"}
    for e in entries.values():
        if e.status == "ok":
            assert len(e.positions) == e.n and check_positions(e.positions, e.panel)["ok"]
    assert len(entries) == len(set(entries))  # unique ids
    est = estimate_seconds(PROFILES["quick"])
    assert est["total_minutes"] > 0 and est["geometries"] > 10


def test_pareto_flags_on_constructed_entries_and_export(tmp_path):
    def mk(id_, n, env, sel, rob, fam="V"):
        e = Entry(id_, fam, n, env, {"A"})
        e.positions = regular_grid(2, 4, ENVELOPES[env]) if n == 8 else regular_grid(4, 4, ENVELOPES[env])
        e.summary = {b: {"val": dict(F10_mean=sel, F10_worst=rob, F10_p10_lines_worst=rob), "train": {}, "all": {}}
                     for b in BANDS}
        return e

    a, b, c, d = mk("a", 8, "E420", 0.4, 0.3), mk("b", 8, "E420", 0.3, 0.2), mk("c", 16, "E420", 0.5, 0.4), mk("d", 8, "E260", 0.3, 0.2)
    flags = pareto_flags([a, b, c, d])
    assert flags == {"a": True, "b": False, "c": True, "d": True}  # b is dominated by a; d wins on area
    assert objectives(a)["area_mm2"] == 420 * 420
    names = export_finalists([a, b, c, d], flags, tmp_path)
    assert {"a", "c", "d"} <= set(names)
    for n in names:
        arr = load_array(tmp_path / f"{n}.json")
        assert check_positions(positions_of(arr), ENVELOPES["E420" if n != "d" else "E260"])["ok"]


def test_evaluate_entry_caches_on_content(tmp_path):
    e = Entry("R8-E420", "R", 8, "E420", {"A"}, positions=regular_grid(2, 4, P420))
    logs = []
    evaluate_entry(e, 2000.0, tmp_path, logs.append)
    first = {k: v["F10"] for k, v in e.ref["cells"].items()}
    assert e.summary["b24"]["val"]["F10_mean"] > 0 and len(first) == 30 and "reference evaluation" in logs[0]
    e2 = Entry("R8-E420", "R", 8, "E420", {"A"}, positions=regular_grid(2, 4, P420))
    evaluate_entry(e2, 2000.0, tmp_path, logs.append)
    assert "cached" in logs[1] and {k: v["F10"] for k, v in e2.ref["cells"].items()} == pytest.approx(first)
    e3 = Entry("R8-E420", "R", 8, "E420", {"A"}, positions=regular_grid(2, 4, P420) + 1.0)  # different layout, same id
    evaluate_entry(e3, 2000.0, tmp_path, logs.append)
    assert "reference evaluation" in logs[2]


# ---- command line ----------------------------------------------------------------------------------------------------------------------------------

def test_cli_generate_evaluate_optimize(tmp_path, capsys):
    out = tmp_path / "v16.json"
    assert main(["array-generate", "--family", "V", "--n", "16", "--out", str(out)]) == 0
    assert "minimum center distance" in capsys.readouterr().out
    assert check_positions(positions_of(load_array(out)), P420)["ok"]
    assert main(["array-evaluate", "--array", str(out), "--foci", "val", "--df", "1000"]) == 0
    text = capsys.readouterr().out
    assert "16 elements" in text and "2-4 kHz" in text and "equal total power" in text
    opt = tmp_path / "o8.json"
    assert main(["array-optimize", "--n", "8", "--init", "V", "--seeds", "0", "--iterations", "60", "--out", str(opt)]) == 0
    assert "J" in capsys.readouterr().out
    assert json.loads(opt.read_text())["meta"]["family"] == "O" and check_positions(positions_of(load_array(opt)), P420)["ok"]
    assert main(["array-campaign", "--profile", "quick", "--estimate"]) == 0
    assert "estimated about" in capsys.readouterr().out


@pytest.mark.parametrize("args", [["--family", "R", "--n", "32"], ["--family", "V", "--n", "32"],
                                  ["--family", "R", "--n", "8", "--rows", "3", "--cols", "3"],
                                  ["--family", "J", "--n", "8", "--max-shift", "-1"],
                                  ["--family", "R", "--n", "8", "--min-distance", "0"]])
def test_cli_generate_rejects_infeasible_or_invalid_requests(tmp_path, args, capsys):
    with pytest.raises(SystemExit) as exc:
        main(["array-generate", *args, "--out", str(tmp_path / "x.json")])
    assert exc.value.code == 2 and not (tmp_path / "x.json").exists()


def test_cli_evaluate_refuses_layouts_that_break_the_constraints(tmp_path):
    bad = tmp_path / "bad.json"
    save_array_json(bad, "bad", np.array([[0.0, 0.0], [10.0, 0.0]]), Panel(420, 420, 32.5, 1.0))
    with pytest.raises(SystemExit) as exc:
        main(["array-evaluate", "--array", str(bad), "--df", "1000"])
    assert exc.value.code == 2


def test_committed_finalists_are_valid_and_loadable_by_the_existing_engine():
    folder = Path(__file__).resolve().parents[1] / "configs/arrays/optimized"
    files = sorted(folder.glob("*.json"))
    assert len(files) >= 10
    for path in files:
        obj = json.loads(path.read_text())
        pn = obj["panel"]
        panel = Panel(pn["width"], pn["height"], pn["margin"], pn["min_distance"], pn["diameter"])
        pos = positions_of(load_array(path))
        assert len(pos) == obj["meta"]["n"] and check_positions(pos, panel)["ok"], path.name
        assert "hypothesis" in pn["note"]
    # the engine of the earlier missions reads them directly
    from bean_lab.field import complex_field
    arr = load_array(folder / "O16-E420.json")
    assert abs(complex_field(arr, np.array([0.0, 0.0, 0.8]), 3000.0, gain=0.05)) > 0
