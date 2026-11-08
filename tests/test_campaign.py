import json
import math
from pathlib import Path

import numpy as np
import pytest

from bean_lab.array_config import load_array
from bean_lab.campaign import (ARC_LIMIT_DEG, BAND, Params, Region, Scenario, analyze_lobes,
                               aperture_center, build_arcs, build_plate, grating_angles,
                               refine_peak, run_campaign, run_scenario, scenarios, uniform_spacing)
from bean_lab.cli import main
from bean_lab.field import complex_field
from bean_lab.focus import SPEED_OF_SOUND

ROOT = Path(__file__).resolve().parents[1]
ARRAYS = ROOT / "configs/arrays"
C = SPEED_OF_SOUND
SMALL = Region(nx=41, ny=39)  # 50 mm steps: the campaign's foci (x = 210, 350, 800, -380; y = 800) are on-grid
FAST = Params(oversample=1)


def far_field_arc(spacing_mm, freq, steer_deg=0.0, radius=500.0):
    """Arc at a far-field radius; returns (theta, level_dB). Independent of campaign internals."""
    from bean_lab.array_config import ArrayConfig, Element
    from bean_lab.focus import focus_delays
    arr = ArrayConfig("t", tuple(Element(i + 1, (spacing_mm * i / 1000, 0.0, 0.0)) for i in range(8)))
    cx = arr.elements[-1].position[0] / 2
    theta = np.arange(-80, 80.0001, 0.02)
    t = np.radians(theta)
    pts = np.stack([cx + radius * np.sin(t), radius * np.cos(t), np.zeros_like(t)], axis=-1)
    em = None
    if steer_deg:
        focus = (cx + radius * math.sin(math.radians(steer_deg)), radius * math.cos(math.radians(steer_deg)), 0.0)
        em = {d.speaker: d.delay for d in focus_delays(arr, focus)}
    p = complex_field(arr, pts, freq, emission_delays=em, gain=1.0)
    return theta, 20 * np.log10(np.abs(p) * radius)  # remove the common 1/R factor


# ---- scenario definitions -----------------------------------------------------------------

def test_scenarios_match_the_campaign_specification():
    f = scenarios("F")
    assert [s.frequency for s in f] == [500, 1000, 2000, 4000, 6000, 8000, 12000, 16000]
    assert {s.array for s in f} == {"line-8-60mm.json"} and {s.focus_mm for s in f} == {(350.0, 800.0, 0.0)}
    s = scenarios("S")
    assert [load_array(ARRAYS / x.array).elements[-1].position[0] for x in s] == pytest.approx([0.21, 0.42, 0.84])
    assert {x.frequency for x in s} == {4000.0} and {x.focus_mm for x in s} == {(350.0, 800.0, 0.0)}
    d = scenarios("D")
    assert [x.focus_mm[0] for x in d] == [210, 350, 800, -380] and {x.frequency for x in d} == {4000.0}
    with pytest.raises(ValueError):
        scenarios("X")


def test_spacing_apertures_and_center():
    for name, pitch in (("line-8-30mm", 0.03), ("line-8-60mm", 0.06), ("line-8-120mm", 0.12)):
        arr = load_array(ARRAYS / f"{name}.json")
        assert uniform_spacing(arr) == pytest.approx(pitch)
        assert aperture_center(arr)[0] == pytest.approx(3.5 * pitch)
    assert uniform_spacing(load_array(ARRAYS / "staggered-3-2-3.json")) is None


def test_800_and_minus_380_are_exact_mirrors_about_the_aperture_center():
    cx = aperture_center(load_array(ARRAYS / "line-8-60mm.json"))[0] * 1000
    assert cx == pytest.approx(210)
    assert 2 * cx - 800 == pytest.approx(-380)  # the lateral-opposite focus IS the mirror of x = 800


# ---- grating lobes and lobe analysis, against closed-form expectations ----------------------------------

def test_grating_angles_closed_form():
    lam = C / 8000
    assert grating_angles(0.06, 8000, C, 0.0) == pytest.approx(
        [-math.degrees(math.asin(lam / 0.06)), math.degrees(math.asin(lam / 0.06))])
    assert grating_angles(0.03, 4000, C, 17.0) == []  # lambda/d = 2.86: never visible
    assert grating_angles(0.06, 4000, C, 0.0) == []  # d < lambda: none at broadside
    g = grating_angles(0.06, 4000, C, math.degrees(math.atan2(0.59, 0.8)))
    assert g == pytest.approx([-56.7], abs=0.1)  # focus x = 800: sin(36.4) - 1.4286


def test_far_field_first_null_and_sidelobe_of_uniform_line():
    theta, level = far_field_arc(60, 2000)  # d/lambda = 0.35: no grating lobes
    lb = analyze_lobes(theta, level)
    lam, n, d = C / 2000, 8, 0.06
    first_null = math.degrees(math.asin(lam / (n * d)))
    assert lb["null_right_deg"] == pytest.approx(first_null, abs=0.05)
    assert lb["null_left_deg"] == pytest.approx(-first_null, abs=0.05)
    # independent array factor |sin(N x) / (N sin x)|, x = pi d sin(theta) / lambda
    x = np.pi * d * np.sin(np.radians(theta)) / lam
    with np.errstate(all="ignore"):
        af = np.abs(np.sin(n * x) / (n * np.sin(x)))
    af[np.isnan(af)] = 1.0
    outside = np.abs(theta) > first_null + 0.05
    expected_sll = -20 * np.log10(af[outside].max())
    assert lb["contrast_db"] == pytest.approx(expected_sll, abs=0.1)
    assert 12.5 < lb["contrast_db"] < 13.3  # classic ~ -12.8 dB for N = 8


def test_beamwidth_matches_the_array_factor_half_power_point():
    n, freq = 8, 2000.0
    lam = C / freq
    for spacing in (30, 60, 120):
        d = spacing / 1000
        lo, hi = 1e-9, math.pi / n  # bisection for |sin(N x) / (N sin x)| = 1/sqrt(2)
        for _ in range(100):
            mid = (lo + hi) / 2
            lo, hi = (mid, hi) if abs(math.sin(n * mid) / (n * math.sin(mid))) > 2 ** -0.5 else (lo, mid)
        expected = 2 * math.degrees(math.asin(lo * lam / (math.pi * d)))
        theta, level = far_field_arc(spacing, freq)
        assert analyze_lobes(theta, level)["beamwidth_deg"] == pytest.approx(expected, abs=0.1), spacing


def test_grating_lobe_appears_where_predicted_with_main_lobe_level():
    theta, level = far_field_arc(120, 4000)  # d/lambda = 1.4, broadside
    lb = analyze_lobes(theta, level)
    predicted = grating_angles(0.12, 4000, C, 0.0)
    assert predicted == pytest.approx([-45.6, 45.6], abs=0.1)
    for g in predicted:
        i = int(np.argmin(np.abs(theta - g)))
        window = slice(max(i - 100, 0), i + 100)
        assert level[window].max() == pytest.approx(lb["peak_db"], abs=0.05)  # equal to the main lobe
        assert theta[window][np.argmax(level[window])] == pytest.approx(g, abs=0.1)
    assert lb["contrast_db"] < 0.1


def test_steered_far_field_peak_follows_the_steering_angle():
    theta, level = far_field_arc(60, 2000, steer_deg=25.0)
    assert analyze_lobes(theta, level)["peak_deg"] == pytest.approx(25.0, abs=0.1)


def test_analyze_lobes_on_a_synthetic_pattern():
    theta = np.linspace(-80, 80, 1601)
    level = -40.0 + 30 * np.exp(-(theta / 10) ** 2) + 12 * np.exp(-((theta - 40) / 5) ** 2) \
        + 6 * np.exp(-((theta + 50) / 5) ** 2)
    lb = analyze_lobes(theta, level)
    assert lb["peak_deg"] == pytest.approx(0, abs=0.1)
    assert lb["sidelobe_deg"] == pytest.approx(40, abs=0.2)
    assert lb["contrast_db"] == pytest.approx(30 - 12, abs=0.2)
    assert len(lb["sidelobes"]) == 2
    assert lb["beamwidth_deg"] == pytest.approx(2 * 10 * math.sqrt(-math.log(0.9)), rel=0.01)  # 30 (1 - exp(-u^2)) = 3 dB
    mono = analyze_lobes(theta, -np.abs(theta))  # single monotone lobe: no sidelobes
    assert mono["contrast_db"] is None


def test_refine_peak_finds_the_known_maximum_and_respects_bounds():
    lam = C / 4000
    arr = load_array(ARRAYS / "line-8-60mm.json")
    # two in-phase sources, mic line at y = 1 m, symmetric about their midpoint: maximum at x = midpoint
    from bean_lab.array_config import ArrayConfig, Element
    pair = ArrayConfig("p", (Element(1, (-lam / 2, 0.0, 0.0)), Element(2, (lam / 2, 0.0, 0.0))))
    x, y, v = refine_peak(pair, 4000.0, None, dict(speed_of_sound=C, gain=1.0), 0.004, 1.0, 0.0, 0.01,
                          (-1.0, 1.0, 0.5, 1.5))
    assert x == pytest.approx(0.0, abs=2e-6)
    assert v == pytest.approx(2 / math.hypot(lam / 2, y), rel=1e-6) or v >= 2 / math.hypot(lam / 2, 1.5)
    # bounds: a single source's |P| grows toward it; the search stops at the boundary
    single = ArrayConfig("s", (Element(1, (0.0, 0.0, 0.0)),))
    _, y, _ = refine_peak(single, 4000.0, None, dict(speed_of_sound=C, gain=1.0), 0.0, 1.0, 0.0, 0.5,
                          (-1.0, 1.0, 0.9, 2.0))
    assert y >= 0.9 - 1e-12


# ---- scenario indicators ----------------------------------------------------------------------------------

def run_small(sc, region=SMALL, params=FAST):
    return run_scenario(sc, region, params, ARRAYS)


def test_indicators_are_consistent_with_maps_and_closed_forms():
    sc = scenarios("D")[1]  # focus (350, 800)
    r = run_small(sc)
    ind = r.indicators
    arr = load_array(ARRAYS / sc.array)
    focus = (0.35, 0.8, 0.0)
    coherent = 20 * math.log10(sum(1 / math.dist(focus, e.position) for e in arr.elements))
    assert ind["level_foc_focus_db"] == pytest.approx(coherent, abs=1e-9)  # delays align all phases
    ix, iy = int(np.argmin(np.abs(r.x_mm - 350))), int(np.argmin(np.abs(r.y_mm - 800)))
    assert (r.x_mm[ix], r.y_mm[iy]) == pytest.approx((350, 800))  # focus is on the 50 mm grid
    assert r.level_foc[iy, ix] == pytest.approx(ind["level_foc_focus_db"], abs=1e-9)
    assert r.level_sim[iy, ix] == pytest.approx(ind["level_sim_focus_db"], abs=1e-9)
    assert r.diff[iy, ix] == pytest.approx(ind["focusing_gain_db"], abs=1e-9)
    assert ind["focusing_gain_db"] == pytest.approx(ind["level_foc_focus_db"] - ind["level_sim_focus_db"])
    assert ind["wavelength_mm"] == pytest.approx(85.75) and ind["spacing_over_wavelength"] == pytest.approx(60 / 85.75)
    assert ind["focus_range_mm"] == pytest.approx(math.hypot(140, 800)) and ind["aperture_mm"] == pytest.approx(420)
    assert ind["focus_angle_deg"] == pytest.approx(math.degrees(math.atan2(140, 800)))
    # arc peak is the focus direction for a (finite-range) focused beam
    assert ind["arc_foc_peak_deg"] == pytest.approx(ind["focus_angle_deg"], abs=0.3)
    assert ind["arc_foc_peak_db"] == pytest.approx(ind["level_foc_focus_db"], abs=0.15)


def test_global_maximum_is_at_least_the_grid_maximum_and_is_not_the_focus():
    r = run_small(scenarios("D")[1])
    ind = r.indicators
    assert ind["max_foc_db"] >= r.level_foc.max() - 1e-9
    assert ind["max_sim_db"] >= r.level_sim.max() - 1e-9
    assert ind["max_foc_db"] >= ind["max_foc_grid_db"] - 1e-9
    assert ind["max_foc_db"] > ind["level_foc_focus_db"] + 3  # near-field maximum dominates
    assert math.dist((ind["max_foc_x_mm"], ind["max_foc_y_mm"]), (350, 800)) > 300
    assert 100 - 1e-6 <= ind["max_foc_y_mm"] <= 2000 + 1e-6 and -700 - 1e-6 <= ind["max_foc_x_mm"] <= 1300 + 1e-6


def test_band_sidelobe_is_outside_the_beam_and_consistent_with_the_arc():
    r = run_small(scenarios("D")[1], Region(nx=21, ny=19), Params(oversample=1))
    ind = r.indicators
    assert ind["band_sidelobe_db"] is not None
    cx = 210.0
    rng = math.hypot(ind["band_sidelobe_x_mm"] - cx, ind["band_sidelobe_y_mm"])
    assert rng == pytest.approx(ind["band_sidelobe_range_mm"], rel=1e-9)
    r_f = ind["focus_range_mm"]
    assert (1 - BAND) * r_f - 1e-6 <= rng <= (1 + BAND) * r_f + 1e-6
    assert ind["band_contrast_db"] == pytest.approx(ind["level_foc_focus_db"] - ind["band_sidelobe_db"])
    # the circle through the focus is one of the band circles, so the band can only be worse
    # than the arc analysis, by at most the residual spreading-compensation error
    assert ind["band_contrast_db"] <= ind["arc_foc_contrast_db"] + 0.3
    assert ind["band_contrast_db"] > ind["arc_foc_contrast_db"] - 6
    # the hourglass main lobe is excluded: the contrast stays near the ~13 dB uniform-array value
    assert 8 < ind["band_contrast_db"] < 14


def test_mirrored_foci_give_mirrored_fields_and_equal_indicators():
    # region symmetric about the aperture center x = 210 mm
    region = Region(x_mm=(-700.0, 1120.0), nx=37, ny=19)
    arr = "line-8-60mm.json"
    a = run_small(Scenario("D", "a", "a", arr, 4000.0, (800.0, 800.0, 0.0)), region)
    b = run_small(Scenario("D", "b", "b", arr, 4000.0, (-380.0, 800.0, 0.0)), region)
    assert np.allclose(a.level_foc, b.level_foc[:, ::-1], atol=1e-9)
    assert np.allclose(a.level_sim, b.level_sim[:, ::-1], atol=1e-9)
    assert np.allclose(a.arc_foc, b.arc_foc[::-1], atol=1e-9)
    for key in ("level_foc_focus_db", "level_sim_focus_db", "focusing_gain_db", "arc_foc_peak_db",
                "arc_foc_beamwidth3_deg", "arc_foc_contrast_db", "max_foc_db", "band_sidelobe_db",
                "band_contrast_db"):
        assert a.indicators[key] == pytest.approx(b.indicators[key], abs=1e-6), key
    assert a.indicators["max_foc_x_mm"] + b.indicators["max_foc_x_mm"] == pytest.approx(420, abs=1e-3)
    assert a.indicators["arc_foc_peak_deg"] == pytest.approx(-b.indicators["arc_foc_peak_deg"], abs=1e-9)
    # a mirror of the x = 350 focus (x = 70) behaves the same way, as the mission suggested
    c = run_small(Scenario("D", "c", "c", arr, 4000.0, (350.0, 800.0, 0.0)), region)
    d = run_small(Scenario("D", "d", "d", arr, 4000.0, (70.0, 800.0, 0.0)), region)
    assert np.allclose(c.level_foc, d.level_foc[:, ::-1], atol=1e-9)


def test_run_is_reproducible():
    sc = scenarios("S")[2]
    r1, r2 = run_small(sc), run_small(sc)
    for key in ("simultaneous", "focused", "level_sim", "level_foc", "arc_foc", "arc_sim"):
        assert np.array_equal(getattr(r1, key), getattr(r2, key))
    assert np.array_equal(r1.diff, r2.diff, equal_nan=True)
    assert r1.indicators == r2.indicators


def test_oversampling_never_lowers_the_maximum():
    sc = scenarios("F")[-1]  # 16 kHz: 5 mm images undersample the pattern
    coarse = run_small(sc, Region(nx=81, ny=77), Params(oversample=1))
    fine = run_small(sc, Region(nx=81, ny=77), Params(oversample=4))
    assert fine.indicators["max_foc_db"] >= coarse.indicators["max_foc_db"] - 1e-6
    assert np.array_equal(coarse.focused, fine.focused)  # the image data do not depend on it


# ---- behavior with frequency, spacing and direction ----------------------------------------------------------

@pytest.fixture(scope="module")
def series_results():
    tiny = Region(nx=5, ny=5)
    return {s: [run_small(sc, tiny) for sc in scenarios(s)] for s in "FSD"}


def test_focus_level_is_frequency_independent_but_gain_grows_with_it(series_results):
    f = series_results["F"]
    levels = [r.indicators["level_foc_focus_db"] for r in f]
    assert max(levels) - min(levels) < 1e-9  # coherent sum of 1/r depends only on geometry
    gains = [r.indicators["focusing_gain_db"] for r in f]
    assert all(g > 0 for g in gains) and gains[0] < 0.5 and gains[-1] > 8  # delays matter more at high f


def test_beamwidth_narrows_with_frequency_and_aperture(series_results):
    bw = [r.indicators["arc_foc_beamwidth3_deg"] for r in series_results["F"]]
    assert all(a > b for a, b in zip(bw, bw[1:]))
    assert bw[0] > 60 and bw[3] == pytest.approx(9.5, abs=0.5)
    bws = [r.indicators["arc_foc_beamwidth3_deg"] for r in series_results["S"]]
    assert bws[0] > bws[1] > bws[2]


def test_sidelobe_contrast_collapses_exactly_when_grating_lobes_are_predicted(series_results):
    for r in series_results["F"] + series_results["S"] + series_results["D"]:
        ind = r.indicators
        contrast = ind["arc_foc_contrast_db"]
        if ind["grating_lobes_predicted_deg"]:
            assert contrast is not None and contrast < 3.0, ind["scenario"]
        elif contrast is not None:
            assert contrast > 10.0, ind["scenario"]
    assert series_results["D"][2].indicators["arc_foc_contrast_db"] < 0.5  # focus x = 800: grating lobe ~ main lobe


def test_extreme_focal_direction_has_larger_main_lobe_and_grating_lobe(series_results):
    d = {r.scenario.key: r.indicators for r in series_results["D"]}
    assert d["x+210"]["arc_foc_beamwidth3_deg"] < d["x+350"]["arc_foc_beamwidth3_deg"] < d["x+800"]["arc_foc_beamwidth3_deg"]
    assert d["x+800"]["arc_foc_sidelobe_deg"] == pytest.approx(-56.7, abs=2.5)
    assert d["x-380"]["arc_foc_sidelobe_deg"] == pytest.approx(56.7, abs=2.5)


# ---- images, scales and outputs ---------------------------------------------------------------------------------

def test_plate_uses_one_absolute_scale_and_plots_the_data():
    params = Params(oversample=1, vmax=20.0, range_db=40.0, diff_limit=10.0)
    results = [run_small(sc, Region(nx=41, ny=39), params) for sc in scenarios("S")]
    plate = build_plate("test", results, params)
    assert len(plate.level_images) == 6 and len(plate.diff_images) == 3
    for im in plate.level_images:
        assert im.get_clim() == (-20.0, 20.0)
    for im in plate.diff_images:
        assert im.get_clim() == (-10.0, 10.0)
    for k, r in enumerate(results):
        assert np.array_equal(np.asarray(plate.level_images[2 * k].get_array()), r.level_sim)
        assert np.array_equal(np.asarray(plate.level_images[2 * k + 1].get_array()), r.level_foc)
        assert np.array_equal(np.ma.filled(plate.diff_images[k].get_array(), np.nan), r.diff, equal_nan=True)
    build_arcs("t", results)


def test_changing_display_range_changes_no_data():
    sc = scenarios("D")[0]
    a = run_small(sc, params=Params(oversample=1, vmax=20))
    b = run_small(sc, params=Params(oversample=1, vmax=5, range_db=10, diff_limit=3))
    assert np.array_equal(a.level_foc, b.level_foc) and a.indicators == b.indicators


def test_campaign_outputs_are_complete_and_consistent(tmp_path):
    logs = []
    results = run_campaign(["S"], tmp_path, SMALL, FAST, log=logs.append, array_dir=ARRAYS)
    assert len(results) == 3 and len(logs) == 3
    for name in ("plate-S.png", "arcs-S.png"):
        assert (tmp_path / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    for r in results:
        d = tmp_path / "S" / r.scenario.key
        for png in ("field-simultaneous.png", "field-focused.png", "field-difference.png"):
            assert (d / png).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
        z = np.load(d / "data.npz")
        assert np.array_equal(z["focused"], r.focused)
        assert json.loads(str(z["indicators"])) == json.loads(json.dumps(r.indicators))
    rows = json.loads((tmp_path / "indicators.json").read_text())
    assert [x["scenario"] for x in rows] == ["S/d030", "S/d060", "S/d120"]
    csv_lines = (tmp_path / "indicators.csv").read_text().splitlines()
    assert len(csv_lines) == 4 and csv_lines[0].startswith("scenario,frequency_hz")


def test_cli_campaign_runs_and_prints_table(tmp_path, capsys):
    rc = main(["campaign", "--series", "S", "--nx", "21", "--ny", "19", "--oversample", "1",
               "--outdir", str(tmp_path / "c"), "--array-dir", str(ARRAYS), "--no-individual"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "S/d030" in out and "S/d120" in out and "contrast" in out
    assert (tmp_path / "c/plate-S.png").exists() and not (tmp_path / "c/S/d030/field-focused.png").exists()
    assert (tmp_path / "c/S/d030/data.npz").exists()


@pytest.mark.parametrize("extra", [["--oversample", "0"], ["--range-db", "0"], ["--diff-limit", "-1"],
                                   ["--nx", "1"], ["--array-dir", "/nonexistent"]])
def test_cli_campaign_rejects_invalid_arguments(tmp_path, extra):
    with pytest.raises(SystemExit) as exc:
        main(["campaign", "--series", "S", "--nx", "21", "--ny", "19", "--outdir", str(tmp_path / "x"),
              "--array-dir", str(ARRAYS), *extra])
    assert exc.value.code == 2
