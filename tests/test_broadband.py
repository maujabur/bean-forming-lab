import json
import math
from pathlib import Path

import numpy as np
import pytest

from bean_lab.array_config import ArrayConfig, Element, load_array
from bean_lab.broadband import (BroadbandScenario, Spectrum, _beamwidth_if_defined, band_energy,
                                broadband_scenarios, connected_component, convergence_study,
                                effective_level_db, map_energies, run_broadband, run_broadband_campaign)
from bean_lab.campaign import BAND, Params, Region, run_scenario, scenarios
from bean_lab.cli import main
from bean_lab.field import FieldGeometry, complex_field, field_maps, grid_axis, level_db
from bean_lab.focus import SPEED_OF_SOUND, focus_delays
from bean_lab.microphone import receive

ROOT = Path(__file__).resolve().parents[1]
ARRAYS = ROOT / "configs/arrays"
LINE8 = ARRAYS / "line-8-60mm.json"
C = SPEED_OF_SOUND
FS = 48000
A = 0.05
FOCUS = (0.35, 0.8, 0.0)
SMALL = Region(nx=41, ny=39)  # 50 mm steps; the campaign foci are on-grid
PARAMS = Params(oversample=1)


def emission_for(array, focus=FOCUS):
    return {d.speaker: d.delay for d in focus_delays(array, focus)}


def make_array(*positions):
    return ArrayConfig("t", tuple(Element(i + 1, tuple(map(float, p))) for i, p in enumerate(positions)))


# ---- spectrum and weights -----------------------------------------------------------------------

@pytest.mark.parametrize("spec", [Spectrum(1000, 8000), Spectrum(2000, 4000, df=7.3), Spectrum(500, 8000, n=3),
                                  Spectrum(1000, 8000, distribution="pink"), Spectrum(4000, 4000), Spectrum(100, 101, n=1)])
def test_weights_are_normalized_and_positive(spec):
    f, w = spec.lines()
    assert w.sum() == pytest.approx(1.0, abs=1e-15) and np.all(w > 0)
    assert len(f) == len(w) and np.all(f >= spec.f_min) and np.all(f <= spec.f_max)


def test_uniform_band_has_equal_weights_and_bin_centers():
    f, w = Spectrum(1000, 2000, df=250).lines()
    assert f.tolist() == [1125, 1375, 1625, 1875] and w.tolist() == [0.25] * 4
    f, w = Spectrum(500, 8000, df=25).lines()
    assert len(f) == 300 and f[0] == 512.5 and f[-1] == 7987.5


def test_pink_weights_follow_one_over_f_power():
    spec = Spectrum(1000, 8000, n=3, distribution="pink")  # linear bins: edges 1000, 3333.3, 5666.7, 8000
    f, w = spec.lines()
    edges = np.linspace(1000, 8000, 4)
    assert w == pytest.approx(np.log(edges[1:] / edges[:-1]) / math.log(8))  # power of a bin = ln(f_hi / f_lo)
    f, w = Spectrum(1000, 8000, n=300, distribution="pink").lines()
    inside = w[(f > 1000) & (f < 2000)].sum()
    assert inside == pytest.approx(math.log(2) / math.log(8), abs=2e-3)
    assert w[0] > w[-1]  # more weight per Hz-bin at low frequency


def test_degenerate_band_is_a_tone_whatever_the_sampling():
    for spec in (Spectrum(4000, 4000), Spectrum(4000, 4000, n=9), Spectrum(4000, 4000, df=1)):
        f, w = spec.lines()
        assert f.tolist() == [4000.0] and w.tolist() == [1.0]
    assert Spectrum(4000, 4000).is_tone and not Spectrum(4000, 4001).is_tone


@pytest.mark.parametrize("kwargs", [dict(f_min=0, f_max=10), dict(f_min=-1, f_max=10), dict(f_min=10, f_max=5),
                                    dict(f_min=1, f_max=math.nan), dict(f_min=1, f_max=math.inf),
                                    dict(f_min=1, f_max=10, df=0), dict(f_min=1, f_max=10, df=math.nan),
                                    dict(f_min=1, f_max=10, n=0), dict(f_min=1, f_max=10, n=2.5),
                                    dict(f_min=1, f_max=10, distribution="white")])
def test_invalid_spectra(kwargs):
    with pytest.raises(ValueError):
        Spectrum(**kwargs)


# ---- energy combination ----------------------------------------------------------------------------

def test_single_line_reproduces_the_monochromatic_field():
    arr = load_array(LINE8)
    pts = np.array([[0.1, 0.5, 0.0], [0.35, 0.8, 0.0], [-0.3, 1.2, 0.0]])
    em = emission_for(arr)
    geo = FieldGeometry(arr, pts, gain=A)
    for emission in (None, em):
        p = complex_field(arr, pts, 4000.0, emission_delays=emission, gain=A)
        e = band_energy(geo, np.array([4000.0]), np.array([1.0]), emission)
        assert np.allclose(e, np.abs(p) ** 2, rtol=1e-14, atol=0)
        assert np.allclose(effective_level_db(e, A, -60.0), level_db(p, A), atol=1e-12)


def test_distinct_frequencies_add_in_power_not_in_phasor():
    arr = load_array(LINE8)
    pt = np.array([[0.2, 0.9, 0.0]])
    geo = FieldGeometry(arr, pt, gain=A)
    f1, f2, w1, w2 = 3000.0, 5300.0, 0.3, 0.7
    p1, p2 = geo.field(f1)[0], geo.field(f2)[0]
    e = band_energy(geo, np.array([f1, f2]), np.array([w1, w2]))[0]
    assert e == pytest.approx(w1 * abs(p1) ** 2 + w2 * abs(p2) ** 2, rel=1e-14)
    assert abs(e - abs(math.sqrt(w1) * p1 + math.sqrt(w2) * p2) ** 2) > 1e-6 * e  # not a phasor sum


def test_equal_total_power_per_speaker_means_equal_energy_for_a_single_source():
    arr = make_array((0, 0, 0))
    pt = np.array([[0.0, 2.0, 0.0]])
    geo = FieldGeometry(arr, pt, gain=A)
    for spec in (Spectrum(4000, 4000), Spectrum(2000, 4000), Spectrum(500, 8000), Spectrum(1000, 8000, distribution="pink")):
        f, w = spec.lines()
        # |P|^2 = (A / r)^2 at every frequency, so the mean square is spectrum independent
        assert band_energy(geo, f, w)[0] == pytest.approx((A / 2.0) ** 2, rel=1e-12)


def test_splitting_spectral_power_into_more_lines_does_not_change_the_result():
    arr = load_array(LINE8)
    pts = np.stack(np.meshgrid(np.linspace(-0.3, 0.8, 5), np.linspace(0.2, 1.5, 4)), -1)
    pts = np.concatenate([pts, np.zeros(pts.shape[:2] + (1,))], axis=-1)
    geo = FieldGeometry(arr, pts, gain=A)
    em = emission_for(arr)
    f, w = np.array([2500.0, 6100.0]), np.array([0.4, 0.6])
    base = band_energy(geo, f, w, em)
    k = 7
    f_split, w_split = np.repeat(f, k), np.repeat(w / k, k)
    assert np.allclose(band_energy(geo, f_split, w_split, em), base, rtol=1e-13)


def test_integration_converges_to_the_closed_form_of_a_two_source_interference():
    # two sources: E(f) = g1^2 + g2^2 + 2 g1 g2 cos(2 pi f dt); band mean has a closed form
    arr = make_array((0, 0, 0), (0.3, 0, 0))
    pt = np.array([[0.1, 0.9, 0.0]])
    geo = FieldGeometry(arr, pt, gain=1.0)
    r1, r2 = math.dist((0.1, 0.9, 0), (0, 0, 0)), math.dist((0.1, 0.9, 0), (0.3, 0, 0))
    g1, g2, dt = 1 / r1, 1 / r2, (r1 - r2) / C
    lo, hi = 1000.0, 8000.0
    exact = g1 ** 2 + g2 ** 2 + 2 * g1 * g2 * (math.sin(2 * math.pi * hi * dt) - math.sin(2 * math.pi * lo * dt)) / (
        2 * math.pi * dt * (hi - lo))
    errors = []
    for n in (20, 40, 80, 160, 320):
        f, w = Spectrum(lo, hi, n=n).lines()
        errors.append(abs(band_energy(geo, f, w)[0] - exact))
    assert errors[-1] < 1e-6 * exact
    assert all(a > 3 * b for a, b in zip(errors, errors[1:]))  # midpoint rule: error ~ 1/n^2


# ---- real delays: coherence at the focus for every frequency ------------------------------------------

def test_focus_is_coherent_at_every_frequency_of_the_band():
    arr = load_array(LINE8)
    em = emission_for(arr)
    f, _ = Spectrum(500, 8000, df=250).lines()
    for freq in f:
        phases = []
        for e in arr.elements:
            r = math.dist(FOCUS, e.position)
            z = np.exp(-2j * math.pi * freq * (em[e.speaker] + r / C))  # independent of the field code
            phases.append(z)
        assert max(abs(z - phases[0]) for z in phases) < 1e-9


def test_focus_level_is_the_coherent_sum_for_any_spectrum():
    arr = load_array(LINE8)
    coherent = 20 * math.log10(sum(1 / math.dist(FOCUS, e.position) for e in arr.elements))
    for bs in broadband_scenarios("B", df=500.0) + broadband_scenarios("D", df=500.0):
        if bs.focus_mm != (350.0, 800.0, 0.0):
            continue
        r = run_broadband(bs, Region(nx=11, ny=11), PARAMS, ARRAYS)
        assert r.indicators["level_foc_focus_db"] == pytest.approx(coherent, abs=1e-9), bs.key
        assert r.indicators["weights_sum"] == pytest.approx(1.0, abs=1e-14)


# ---- scenario-level behavior --------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def b_results():
    return {bs.key: run_broadband(bs, SMALL, PARAMS, ARRAYS) for bs in broadband_scenarios("B", df=100.0)}


def test_b1_reproduces_missions_03_and_04(b_results):
    b1 = b_results["B1"]
    m = run_scenario(scenarios("F")[3], SMALL, PARAMS, ARRAYS)  # 4 kHz tone, same geometry and focus
    assert np.allclose(b1.level_foc, m.level_foc, atol=1e-10) and np.allclose(b1.level_sim, m.level_sim, atol=1e-10)
    assert np.allclose(b1.diff, m.diff, atol=1e-10, equal_nan=True)
    for key in ("level_foc_focus_db", "level_sim_focus_db", "focusing_gain_db"):
        assert b1.indicators[key] == pytest.approx(m.indicators[key], abs=1e-10)
    assert b1.indicators["arc_foc_contrast_db"] == pytest.approx(m.indicators["arc_foc_contrast_db"], abs=1e-9)
    assert b1.indicators["arc_foc_beamwidth3_deg"] == pytest.approx(m.indicators["arc_foc_beamwidth3_deg"], abs=1e-9)
    assert b1.indicators["band_contrast_db"] == pytest.approx(m.indicators["band_contrast_db"], abs=1e-9)
    # Mission 03 field-map numbers
    arr = load_array(LINE8)
    fm = field_maps(arr, FOCUS, 4000.0, grid_axis(-0.7, 1.3, 41, "x"), grid_axis(0.1, 2.0, 39, "y"), 0.0, gain=A)
    assert np.allclose(b1.energy_foc, np.abs(fm.focused) ** 2, rtol=1e-12)


def test_band_map_is_the_weighted_average_of_tone_maps():
    arr = load_array(LINE8)
    spec = Spectrum(2000, 4000, n=5)
    f, w = spec.lines()
    x, y = grid_axis(-0.7, 1.3, 21, "x"), grid_axis(0.1, 2.0, 19, "y")
    e_sim, e_foc = map_energies(arr, FOCUS, spec, x, y, 0.0, PARAMS)
    ref_s = sum(wk * np.abs(field_maps(arr, FOCUS, fk, x, y, 0.0, gain=A).simultaneous) ** 2 for fk, wk in zip(f, w))
    ref_f = sum(wk * np.abs(field_maps(arr, FOCUS, fk, x, y, 0.0, gain=A).focused) ** 2 for fk, wk in zip(f, w))
    assert np.allclose(e_sim, ref_s, rtol=1e-12) and np.allclose(e_foc, ref_f, rtol=1e-12)


def test_percent_thresholds_are_nested_and_shares_are_convex_combinations(b_results):
    for r in b_results.values():
        i = r.indicators
        for tag in ("foc", "sim"):
            a, b, c = (i[f"roi_pct_above_m{t}db_{tag}"] for t in (3, 6, 10))
            assert 0 <= a <= b <= c <= 100
        assert i["narrowband_share_10deg_min"] - 1e-12 <= i["arc_foc_share_10deg"] <= i["narrowband_share_10deg_max"] + 1e-12
        assert i["arc_foc_share_5deg"] <= i["arc_foc_share_10deg"] <= 1
        assert i["focal_region_area_cm2"] is None or i["focal_region_area_cm2"] > 0


def test_broadband_smooths_peaks_but_not_the_energy_budget(b_results):
    b1, b3 = b_results["B1"].indicators, b_results["B3"].indicators
    # contrast of the strongest competing lobe improves, as the narrowband worst cases are averaged
    assert b3["arc_foc_contrast_db"] > b3["narrowband_contrast_min_db"] + 3
    assert b3["narrowband_contrast_min_db"] < 3 and b3["narrowband_fraction_below_3db"] > 0.2
    # energy share in the focal window cannot beat the best narrowband line (convex combination)
    assert b3["arc_foc_share_10deg"] <= b3["narrowband_share_10deg_max"]
    assert b1["narrowband_fraction_below_3db"] == 0.0


def test_lateral_mirror_symmetry_in_wide_band():
    region = Region(x_mm=(-700.0, 1120.0), nx=37, ny=19)  # symmetric about the aperture center x = 210 mm
    spec = Spectrum(1000, 8000, n=14)
    results = [run_broadband(BroadbandScenario("D", f"x{x:g}", "m", "line-8-60mm.json", (x, 800.0, 0.0), spec),
                             region, PARAMS, ARRAYS) for x in (800.0, -380.0)]
    a, b = results
    assert np.allclose(a.level_foc, b.level_foc[:, ::-1], atol=1e-9)
    assert np.allclose(a.level_sim, b.level_sim[:, ::-1], atol=1e-9)
    assert np.allclose(a.arc_foc, b.arc_foc[::-1], atol=1e-9)
    for key in ("level_foc_focus_db", "focusing_gain_db", "arc_foc_contrast_db", "band_contrast_db",
                "arc_foc_share_10deg", "roi_pct_above_m6db_foc", "focal_region_area_cm2"):
        assert a.indicators[key] == pytest.approx(b.indicators[key], abs=1e-6), key
    assert a.indicators["arc_foc_lobe_peak_deg"] == pytest.approx(-b.indicators["arc_foc_lobe_peak_deg"], abs=1e-9)


def test_reproducible_and_independent_of_display_scales():
    bs = broadband_scenarios("B", df=500.0)[2]
    r1 = run_broadband(bs, Region(nx=21, ny=19), PARAMS, ARRAYS)
    r2 = run_broadband(bs, Region(nx=21, ny=19), Params(oversample=1, vmax=5.0, range_db=10.0, diff_limit=3.0), ARRAYS)
    for key in ("energy_sim", "energy_foc", "level_sim", "level_foc", "arc_sim", "arc_foc"):
        assert np.array_equal(getattr(r1, key), getattr(r2, key)), key
    assert np.array_equal(r1.diff, r2.diff, equal_nan=True) and r1.indicators == r2.indicators


def test_spectral_convergence_study():
    bs = broadband_scenarios("B")[2]  # 1-8 kHz
    rows = convergence_study(bs, Region(nx=41, ny=39), PARAMS, (400.0, 200.0, 100.0, 50.0, 25.0), 12.5, ARRAYS)
    errs = [r["max_abs_focused_db"] for r in rows]
    assert errs[0] > errs[-1] and all(a > b for a, b in zip(errs[:-1], errs[1:]))
    assert errs[-1] < 0.1 and rows[-1]["max_abs_difference_db"] < 0.1  # default spacing is converged to ~0.1 dB
    assert [r["n_lines"] for r in rows] == [18, 35, 70, 140, 280]


# ---- lobe helpers ---------------------------------------------------------------------------------------------------

def test_connected_component():
    m = np.zeros((6, 6), dtype=bool)
    m[1:4, 1:4] = True
    m[5, 5] = True
    comp = connected_component(m, (2, 2))
    assert comp.sum() == 9 and not comp[5, 5]
    assert connected_component(m, (0, 0)).sum() == 0
    assert connected_component(np.ones((3, 4), dtype=bool), (0, 0)).sum() == 12


def test_beamwidth_is_withheld_without_a_dominant_closed_lobe():
    ok = dict(main_lobe_open=False, beamwidth_open=False, contrast_db=10.0, beamwidth_deg=7.0)
    assert _beamwidth_if_defined(ok) == (7.0, "defined")
    assert _beamwidth_if_defined(dict(ok, main_lobe_open=True))[0] is None
    assert _beamwidth_if_defined(dict(ok, beamwidth_open=True))[0] is None
    assert _beamwidth_if_defined(dict(ok, contrast_db=1.0))[0] is None
    assert _beamwidth_if_defined(dict(ok, contrast_db=None))[0] == 7.0  # no sidelobes at all: width is defined


def test_wide_band_withholds_width_where_grating_lobes_compete(b_results):
    assert b_results["B1"].indicators["arc_foc_beamwidth3_deg"] == pytest.approx(9.5, abs=0.1)
    assert b_results["B1"].indicators["arc_foc_beamwidth_note"] == "defined"
    b4 = b_results["B4"].indicators
    assert (b4["arc_foc_beamwidth3_deg"] is None) == (b4["arc_foc_beamwidth_note"] != "defined")


# ---- comparison with the virtual microphone --------------------------------------------------------------------------

@pytest.mark.parametrize("mic", [(0.35, 0.8, 0.0), (0.0, 0.8, 0.0), (0.7, 1.4, 0.0)])
@pytest.mark.parametrize("focused", [False, True])
def test_orthogonal_multisine_through_the_virtual_microphone_matches_the_energy_model(mic, focused):
    arr = load_array(LINE8)
    em = emission_for(arr) if focused else None
    lines = np.arange(1000.0, 8001.0, 300.0)  # multiples of 20 Hz: orthogonal over a 0.05 s window
    w = np.full(len(lines), 1 / len(lines))
    rng = np.random.default_rng(3)
    phases = rng.uniform(0, 2 * np.pi, len(lines))  # arbitrary: cross terms vanish for orthogonal lines
    t = np.arange(round(0.25 * FS)) / FS
    src = sum(A * math.sqrt(wk) * np.sin(2 * np.pi * fk * t + ph) for fk, wk, ph in zip(lines, w, phases))
    assert np.mean(src[2400:-2400] ** 2) == pytest.approx(A ** 2 / 2, rel=1e-3)  # unit emitted power
    rec = receive({e.speaker: src for e in arr.elements}, arr, mic, FS, emission_delays=em).signal
    seg = rec[4800:4800 + 4800]  # 0.1 s window = 2 x 0.05 s, after the transient, before the tail
    measured = np.mean(seg ** 2)
    geo = FieldGeometry(arr, np.array([mic]), gain=A)
    predicted = 0.5 * band_energy(geo, lines, w, em)[0]
    assert measured == pytest.approx(predicted, rel=2e-3)


def test_finite_noise_realization_fluctuates_around_the_ideal_spectrum():
    arr = load_array(LINE8)
    em = emission_for(arr)
    mic = (0.1, 1.2, 0.0)
    f, w = Spectrum(1000, 8000, df=25).lines()
    predicted = 0.5 * band_energy(FieldGeometry(arr, np.array([mic]), gain=A), f, w, em)[0]
    stats = {}
    for duration in (0.25, 1.0):
        n = round(duration * FS)
        ratios = []
        for seed in range(12):
            x = np.random.default_rng(seed).standard_normal(n)
            spec = np.fft.rfft(x)
            fr = np.fft.rfftfreq(n, 1 / FS)
            spec[(fr < 1000) | (fr > 8000)] = 0
            y = np.fft.irfft(spec, n) * (A / math.sqrt(2)) / math.sqrt(7000 / (FS / 2))
            rec = receive({e.speaker: y for e in arr.elements}, arr, mic, FS, emission_delays=em).signal
            ratios.append(np.mean(rec[2400:n - 2400] ** 2) / predicted)
        stats[duration] = (np.mean(ratios), np.std(ratios))
    for mean, std in stats.values():
        assert mean == pytest.approx(1.0, abs=0.04)  # unbiased on average...
        assert std > 0.005  # ...but a single realization is not the ideal spectrum
    assert stats[0.25][1] > 1.5 * stats[1.0][1]  # scatter shrinks roughly as 1 / sqrt(duration)


# ---- campaign and CLI -------------------------------------------------------------------------------------------------

def test_campaign_reuses_identical_scenarios_and_writes_outputs(tmp_path):
    logs = []
    results = run_broadband_campaign(["B", "D"], tmp_path, Region(nx=21, ny=19), PARAMS, 1000.0, log=logs.append,
                                     array_dir=ARRAYS)
    assert any("reused" in line and "D/x350" in line for line in logs)  # same as B3
    assert len(results) == 6  # B1-B4 + D x210, x800
    for name in ("broadband-plate-B.png", "broadband-arcs-B.png", "broadband-plate-D.png", "broadband-arcs-D.png"):
        assert (tmp_path / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    rows = json.loads((tmp_path / "broadband-indicators.json").read_text())
    assert [r["scenario"] for r in rows] == ["B/B1", "B/B2", "B/B3", "B/B4", "D/x210", "D/x800"]
    assert (tmp_path / "broadband-indicators.csv").read_text().splitlines()[0].startswith("scenario,")
    z = np.load(tmp_path / "B/B3/data.npz")
    assert np.isclose(z["weights"].sum(), 1.0) and z["energy_focused"].shape == (19, 21)
    assert (tmp_path / "B/B3/field-focused.png").exists()


def test_plates_use_one_absolute_scale(tmp_path):
    from bean_lab.broadband import PLATE_TITLES, title_for
    from bean_lab.campaign import build_plate
    results = [run_broadband(bs, Region(nx=21, ny=19), PARAMS, ARRAYS) for bs in broadband_scenarios("D", 2000.0)]
    plate = build_plate(PLATE_TITLES["D"], results, PARAMS, titler=title_for)
    assert all(im.get_clim() == (-20.0, 20.0) for im in plate.level_images)
    assert all(im.get_clim() == (-10.0, 10.0) for im in plate.diff_images)
    for k, r in enumerate(results):
        assert np.array_equal(np.asarray(plate.level_images[2 * k].get_array()), r.level_sim)
        assert np.array_equal(np.asarray(plate.level_images[2 * k + 1].get_array()), r.level_foc)


def test_cli_broadband_map(tmp_path, capsys):
    rc = main(["broadband-map", "--array", str(LINE8), "--target", "350", "800", "0", "--fmin", "2000",
               "--fmax", "4000", "--nfreq", "9", "--nx", "21", "--ny", "19", "--outdir", str(tmp_path / "m")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "9 lines" in out and "weights sum 1.000000000000" in out and "19.76" in out
    z = np.load(tmp_path / "m/data.npz")
    assert len(z["frequencies_hz"]) == 9
    for name in ("field-simultaneous.png", "field-focused.png", "field-difference.png"):
        assert (tmp_path / "m" / name).exists()
    # a tone through the same command matches the monochromatic map
    main(["broadband-map", "--array", str(LINE8), "--target", "350", "800", "0", "--fmin", "4000", "--fmax", "4000",
          "--nx", "21", "--ny", "19", "--outdir", str(tmp_path / "t")])
    zt = np.load(tmp_path / "t/data.npz")
    arr = load_array(LINE8)
    fm = field_maps(arr, FOCUS, 4000.0, grid_axis(-0.7, 1.3, 21, "x"), grid_axis(0.1, 2.0, 19, "y"), 0.0, gain=A)
    assert np.allclose(zt["energy_focused"], np.abs(fm.focused) ** 2, rtol=1e-12)


@pytest.mark.parametrize("extra", [["--fmin", "5000", "--fmax", "1000"], ["--fmin", "0"], ["--nfreq", "0"],
                                   ["--df", "-1"], ["--range-db", "0"], ["--gain", "0"], ["--nx", "1"],
                                   ["--speed-of-sound", "0"]])
def test_cli_broadband_map_rejects_invalid_arguments(tmp_path, extra):
    args = ["broadband-map", "--array", str(LINE8), "--target", "350", "800", "0", "--fmin", "1000", "--fmax", "2000",
            "--nx", "21", "--ny", "19", "--outdir", str(tmp_path / "x")]
    with pytest.raises(SystemExit) as exc:
        main(args + extra)
    assert exc.value.code == 2


def test_cli_broadband_campaign_with_convergence(tmp_path, capsys):
    rc = main(["broadband-campaign", "--series", "B", "--nx", "21", "--ny", "19", "--df", "1000", "--no-individual",
               "--convergence", "--array-dir", str(ARRAYS), "--outdir", str(tmp_path / "c")])
    assert rc == 0
    assert "B/B3" in capsys.readouterr().out
    lines = (tmp_path / "c/convergence.csv").read_text().splitlines()
    assert lines[0].startswith("scenario,df_hz") and len(lines) == 1 + 3 * 5
    assert not (tmp_path / "c/B/B1/field-focused.png").exists() and (tmp_path / "c/B/B1/data.npz").exists()


def test_cli_broadband_campaign_rejects_bad_options(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["broadband-campaign", "--series", "B", "--array-dir", "/nonexistent", "--outdir", str(tmp_path / "x"),
              "--nx", "21", "--ny", "19", "--df", "1000"])
    assert exc.value.code == 2


def test_earlier_commands_unchanged(tmp_path):
    assert main(["field-map", "--array", str(LINE8), "--target", "350", "800", "0", "--nx", "11", "--ny", "11",
                 "--outdir", str(tmp_path / "f")]) == 0


def test_range_compensation_removes_pure_spreading(tmp_path):
    # one omnidirectional source at the aperture center: after compensating 1/r the whole focal
    # band is exactly as strong as the focus, so every cell counts and the -3 dB region is the band
    (tmp_path / "one.json").write_text(json.dumps({"name": "one", "units": "mm", "elements": [
        {"speaker": 1, "position": [0, 0, 0]}]}))
    bs = BroadbandScenario("X", "one", "one", "one.json", (0.0, 800.0, 0.0), Spectrum(1000, 8000, n=5))
    r = run_broadband(bs, Region(x_mm=(-1000.0, 1000.0), y_mm=(100.0, 2000.0), nx=81, ny=77), PARAMS, tmp_path)
    i = r.indicators
    assert i["roi_pct_above_m3db_foc"] == pytest.approx(100.0)
    assert i["focal_region_area_cm2"] == pytest.approx(i["roi_cells"] * 25 * 25 / 100)
    assert i["focusing_gain_db"] == pytest.approx(0.0, abs=1e-9)  # nothing to focus with one source
    assert i["arc_foc_contrast_db"] is None and i["band_contrast_db"] is None


def test_reference_series_are_tones_with_the_same_geometry_and_foci():
    refs = broadband_scenarios("R")
    assert [r.key for r in refs] == ["d030", "d120", "x210", "x800"]
    assert all(r.spectrum.is_tone and r.spectrum.f_min == 4000.0 for r in refs)
    pairs = {(r.array, r.focus_mm) for r in refs}
    gd = {(b.array, b.focus_mm) for b in broadband_scenarios("G") + broadband_scenarios("D")}
    assert pairs <= gd  # every reference tone has a broadband counterpart
