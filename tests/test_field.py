import cmath
import json
import math
from pathlib import Path

import numpy as np
import pytest

from bean_lab.array_config import ArrayConfig, Element, load_array
from bean_lab.cli import main
from bean_lab.focus import SPEED_OF_SOUND, focus_delays
from bean_lab.field import (complex_field, difference_db, field_maps, grid_axis, level_db)
from bean_lab.fieldplot import save_level_map
from bean_lab.microphone import receive

ROOT = Path(__file__).resolve().parents[1]
LINE8 = ROOT / "configs/arrays/line-8-60mm.json"
C = SPEED_OF_SOUND
FS = 48000


def make_array(*positions):
    return ArrayConfig("t", tuple(Element(i + 1, tuple(map(float, p))) for i, p in enumerate(positions)))


def P(array, point, f=4000.0, **kw):
    return complex_field(array, np.array(point, dtype=float), f, **kw)


def reference_field(positions, delays, point, f, gain):
    """Independent scalar re-implementation with cmath (no shared code with field.py)."""
    total = 0j
    for pos, e in zip(positions, delays):
        r = math.sqrt(sum((a - b) ** 2 for a, b in zip(point, pos)))
        total += gain / r * cmath.exp(-2j * math.pi * f * (e + r / C))
    return total


# ---- single source and pairs ------------------------------------------------------

@pytest.mark.parametrize("r", [0.1, 0.5, 1.0, 2.0, 10.0])
def test_single_source_amplitude_and_phase(r):
    arr = make_array((0, 0, 0))
    p = P(arr, (r, 0, 0), 1000.0, gain=0.05)
    assert abs(p) == pytest.approx(0.05 / r, rel=1e-12)  # 1/r with 1 m reference
    expected_phase = -2 * math.pi * 1000.0 * r / C
    assert cmath.phase(p / cmath.exp(1j * expected_phase)) == pytest.approx(0, abs=1e-9)


def test_reference_distance_is_one_meter():
    arr = make_array((0, 0, 0))
    assert abs(P(arr, (0, 1.0, 0), gain=0.05)) == pytest.approx(0.05)  # 0 dB point
    assert level_db(P(arr, (0, 1.0, 0), gain=0.05), 0.05) == pytest.approx(0.0, abs=1e-12)
    assert level_db(P(arr, (0, 2.0, 0), gain=0.05), 0.05) == pytest.approx(-6.0206, abs=1e-3)


def test_two_sources_in_phase_add_and_in_antiphase_cancel_on_bisector():
    arr = make_array((-0.1, 0, 0), (0.1, 0, 0))
    pt = (0.0, 1.5, 0.0)
    r = math.hypot(0.1, 1.5)
    inphase = P(arr, pt, 1000.0)
    assert abs(inphase) == pytest.approx(2 / r, rel=1e-12)
    anti = P(arr, pt, 1000.0, emission_delays={2: 0.5 / 1000.0})  # half a period
    assert abs(anti) < 1e-12 * 2 / r


def test_half_wavelength_path_difference_residual_is_analytic():
    f = 1000.0
    lam = C / f
    arr = make_array((0, 0, 0), (0, 0, lam / 2))
    p = P(arr, (0, 0, -1.0), f)  # r1 = 1, r2 = 1 + lam/2
    assert abs(p) == pytest.approx(abs(1 - 1 / (1 + lam / 2)), rel=1e-9)


# ---- focusing ------------------------------------------------------------------------

def test_focus_delays_make_all_phasors_coincide_at_the_focus():
    arr = load_array(LINE8)
    focus = (0.35, 0.8, 0.0)
    delays = {d.speaker: d.delay for d in focus_delays(arr, focus)}
    omega = 2 * math.pi * 4000.0
    phases, amps = [], []
    for e in arr.elements:
        r = math.dist(focus, e.position)
        z = 0.05 / r * cmath.exp(-1j * omega * (delays[e.speaker] + r / C))
        phases.append(z / abs(z))
        amps.append(abs(z))
    assert max(abs(p - phases[0]) for p in phases) < 1e-9  # same phase for every speaker
    p = P(arr, focus, 4000.0, emission_delays=delays, gain=0.05)
    assert abs(p) == pytest.approx(sum(amps), rel=1e-9)  # fully coherent sum
    assert 20 * math.log10(abs(p) / 0.05) == pytest.approx(
        20 * math.log10(sum(1 / math.dist(focus, e.position) for e in arr.elements)), abs=1e-9)


def test_focused_field_never_exceeds_coherent_sum_and_beats_simultaneous_at_focus():
    arr = load_array(LINE8)
    focus = (0.35, 0.8, 0.0)
    x, y = grid_axis(-0.5, 1.0, 31, "x"), grid_axis(0.1, 2.0, 39, "y")
    m = field_maps(arr, focus, 4000.0, x, y, 0.0, gain=0.05)
    xx, yy = np.meshgrid(x, y)
    bound = sum(0.05 / np.hypot(xx - e.position[0], yy - e.position[1]) for e in arr.elements)
    assert np.all(np.abs(m.focused) <= bound * (1 + 1e-12))
    assert np.all(np.abs(m.simultaneous) <= bound * (1 + 1e-12))
    pf = P(arr, focus, 4000.0, gain=0.05, emission_delays={
        d.speaker: d.delay for d in focus_delays(arr, focus)})
    ps = P(arr, focus, 4000.0, gain=0.05)
    assert abs(pf) > abs(ps)


def test_symmetry_of_symmetric_array():
    arr = load_array(LINE8)  # symmetric about x = 0.21
    x, y = grid_axis(-0.5, 0.92, 29, "x"), grid_axis(0.1, 2.0, 20, "y")  # symmetric about 0.21
    for focus in [(0.21, 0.8, 0.0)]:
        m = field_maps(arr, focus, 4000.0, x, y, 0.0, gain=0.05)
        assert np.allclose(m.simultaneous, m.simultaneous[:, ::-1], atol=1e-12)
        assert np.allclose(m.focused, m.focused[:, ::-1], atol=1e-12)
    off = field_maps(arr, (0.35, 0.8, 0.0), 4000.0, x, y, 0.0, gain=0.05)
    assert np.allclose(off.simultaneous, off.simultaneous[:, ::-1], atol=1e-12)  # unfocused: still symmetric
    assert not np.allclose(off.focused, off.focused[:, ::-1], atol=1e-6)  # off-center focus breaks it


# ---- numerical values and orientation ------------------------------------------------------

def test_map_matches_independent_implementation_and_orientation():
    arr = load_array(LINE8)
    positions = [e.position for e in arr.elements]
    focus = (0.35, 0.8, 0.0)
    delays = [d.delay for d in focus_delays(arr, focus)]
    x, y = grid_axis(-0.5, 1.0, 7, "x"), grid_axis(0.1, 2.0, 5, "y")
    m = field_maps(arr, focus, 4000.0, x, y, 0.0, gain=0.05)
    assert m.simultaneous.shape == (5, 7)  # (ny, nx)
    for iy in range(5):
        for ix in range(7):
            pt = (x[ix], y[iy], 0.0)
            assert m.simultaneous[iy, ix] == pytest.approx(
                reference_field(positions, [0] * 8, pt, 4000.0, 0.05), rel=1e-10, abs=1e-14)
            assert m.focused[iy, ix] == pytest.approx(
                reference_field(positions, delays, pt, 4000.0, 0.05), rel=1e-10, abs=1e-14)


def test_experiment_levels_at_focus_are_analytic():
    arr = load_array(LINE8)
    focus = (0.35, 0.8, 0.0)
    rs = [math.dist(focus, e.position) for e in arr.elements]
    focused_db = 20 * math.log10(sum(1 / r for r in rs))
    m = field_maps(arr, focus, 4000.0, np.array([0.35]), np.array([0.8]), 0.0, gain=0.05)
    assert level_db(m.focused, 0.05)[0, 0] == pytest.approx(focused_db, abs=1e-9)
    assert focused_db == pytest.approx(19.76, abs=0.01)
    assert level_db(m.simultaneous, 0.05)[0, 0] == pytest.approx(12.49, abs=0.01)  # documented value


# ---- dB handling of near-zero values ------------------------------------------------------

def test_exact_cancellation_is_finite_in_db_and_difference_is_nan_when_both_vanish():
    ref = 0.05
    zero = np.array([0j, 1e-30 + 0j])
    db = level_db(zero, ref, -60.0)
    assert np.all(np.isfinite(db)) and np.all(db == -60.0)
    one = np.array([ref * 1j, ref + 0j])
    d = difference_db(np.array([0j, ref]), np.array([0j, ref]), ref)
    assert np.isnan(d[0]) and d[1] == 0.0
    d = difference_db(np.array([ref]), np.array([0j]), ref)  # one side below floor: clamped, finite
    assert d[0] == pytest.approx(60.0)
    assert np.all(np.isfinite(level_db(one, ref)))


def test_near_complete_cancellation_in_a_real_field():
    arr = make_array((-0.1, 0, 0), (0.1, 0, 0))
    p = P(arr, np.array([[0.0, 1.0, 0.0]]), 1000.0, emission_delays={2: 0.5e-3}, gain=0.05)
    assert abs(p[0]) < 1e-15
    assert level_db(p, 0.05)[0] == -60.0  # clamped for display only; |P| itself is kept


# ---- consistency with the virtual microphone -----------------------------------------------------

@pytest.mark.parametrize("point", [(0.35, 0.8, 0.0), (-0.2, 0.5, 0.0), (0.9, 1.5, 0.0), (0.21, 0.12, 0.0)])
@pytest.mark.parametrize("focused", [False, True])
def test_matches_virtual_microphone_amplitude_and_phase(point, focused):
    arr = load_array(LINE8)
    f, gain = 4000.0, 0.05
    em = {d.speaker: d.delay for d in focus_delays(arr, (0.35, 0.8, 0.0))} if focused else None
    t = np.arange(round(0.5 * FS)) / FS
    src = gain * np.sin(2 * np.pi * f * t)
    rec = receive({e.speaker: src for e in arr.elements}, arr, point, FS, emission_delays=em).signal
    i0, i1 = 4800, 4800 + 12 * 1000  # whole periods (12 samples), past all onsets
    n = np.arange(i0, i1)
    w = 2 * np.pi * f * n / FS
    s, c = np.mean(rec[i0:i1] * np.sin(w)), np.mean(rec[i0:i1] * np.cos(w))
    amp, phase = 2 * math.hypot(s, c), math.atan2(c, s)  # y = amp * sin(w + phase)
    p = P(arr, point, f, emission_delays=em, gain=gain)
    assert amp == pytest.approx(abs(p), rel=2e-4)
    assert abs(cmath.exp(1j * phase) - cmath.exp(1j * cmath.phase(p))) < 2e-4
    assert amp / math.sqrt(2) == pytest.approx(abs(p) / math.sqrt(2), rel=2e-4)  # RMS convention


# ---- invalid input ----------------------------------------------------------------------------------

def test_point_on_a_speaker_is_rejected():
    arr = load_array(LINE8)
    with pytest.raises(ValueError, match="coincides with speaker 3"):
        P(arr, (0.12, 0.0, 0.0))
    with pytest.raises(ValueError, match="coincides"):
        complex_field(arr, np.zeros((2, 3, 3)) + [0.06, 0.0, 0.0], 4000.0)
    assert abs(P(arr, (0.12, 0.001, 0.0))) > 0  # 1 mm away is valid


@pytest.mark.parametrize("kwargs", [{"frequency": 0}, {"frequency": -1}, {"frequency": math.nan},
                                    {"speed_of_sound": 0}, {"gain": 0}, {"gain": math.inf},
                                    {"reference_distance": -1}, {"emission_delays": {99: 0.0}},
                                    {"emission_delays": {1: -1e-6}}, {"emission_delays": {1: math.nan}}])
def test_invalid_parameters(kwargs):
    arr = load_array(LINE8)
    args = dict(frequency=4000.0)
    args.update(kwargs)
    f = args.pop("frequency")
    with pytest.raises(ValueError):
        complex_field(arr, np.array([0.3, 1.0, 0.0]), f, **args)


@pytest.mark.parametrize("pts", [np.array([0.0, 1.0]), np.array([[math.nan, 1.0, 0.0]]),
                                 np.array([[math.inf, 1.0, 0.0]])])
def test_invalid_points(pts):
    with pytest.raises(ValueError):
        complex_field(load_array(LINE8), pts, 4000.0)


@pytest.mark.parametrize("args", [(0, 1, 0), (0, 1, 1.5), (1, 0, 5), (math.nan, 1, 5), (0, 1, 1),
                                  (1, 1, 5), (0, math.inf, 5)])
def test_invalid_grid_axes(args):
    with pytest.raises(ValueError):
        grid_axis(*args, "x")


def test_grid_axis_single_point_and_orientation():
    assert grid_axis(0.5, 0.5, 1, "x").tolist() == [0.5]
    assert grid_axis(0.0, 1.0, 5, "x").tolist() == [0.0, 0.25, 0.5, 0.75, 1.0]


# ---- plotting and CLI ---------------------------------------------------------------------------------

def run_cli(tmp_path, *extra, outdir="out"):
    return main(["field-map", "--array", str(LINE8), "--target", "350", "800", "0",
                 "--frequency", "4000", "--gain", "0.05", "--x-range", "-500", "1000",
                 "--y-range", "100", "2000", "--nx", "31", "--ny", "39",
                 "--outdir", str(tmp_path / outdir), *extra])


def test_cli_writes_three_pngs_and_data(tmp_path, capsys):
    assert run_cli(tmp_path, "--probe", "0", "800", "0") == 0
    out = tmp_path / "out"
    for name in ("field-simultaneous.png", "field-focused.png", "field-difference.png"):
        data = (out / name).read_bytes()
        assert data[:8] == b"\x89PNG\r\n\x1a\n" and len(data) > 5000
    z = np.load(out / "field-data.npz")
    assert z["focused"].shape == (39, 31) and z["focused"].dtype == np.complex128
    assert z["x_mm"][0] == -500 and z["y_mm"][-1] == 2000
    params = json.loads(str(z["parameters"]))
    assert params["frequency_hz"] == 4000.0 and params["focus_mm"] == [350, 800, 0]
    arr = load_array(LINE8)
    ix = int(np.argmin(np.abs(z["x_mm"] - 0)))
    iy = int(np.argmin(np.abs(z["y_mm"] - 100)))
    assert z["simultaneous"][iy, ix] == pytest.approx(P(arr, (0.0, 0.1, 0.0), 4000.0, gain=0.05), rel=1e-10)
    text = capsys.readouterr().out
    assert "focus" in text and "19.76" in text and "Map maximum" in text and "probe 1" in text


def test_color_range_options_do_not_change_data(tmp_path, capsys):
    run_cli(tmp_path, outdir="a")
    run_cli(tmp_path, "--vmax", "5", "--range-db", "10", "--diff-limit", "3", outdir="b")
    a, b = np.load(tmp_path / "a/field-data.npz"), np.load(tmp_path / "b/field-data.npz")
    for key in ("simultaneous", "focused", "simultaneous_db", "focused_db"):
        assert np.array_equal(a[key], b[key])
    assert np.array_equal(a["difference_db"], b["difference_db"], equal_nan=True)
    assert "Level color scale: -5.00 to 5.00 dB; difference scale +/-3 dB" in capsys.readouterr().out


def test_first_two_maps_share_one_absolute_color_scale(tmp_path, capsys):
    run_cli(tmp_path)
    z = np.load(tmp_path / "out/field-data.npz")
    vmax = max(z["simultaneous_db"].max(), z["focused_db"].max())
    assert f"Level color scale: {vmax - 40:.2f} to {vmax:.2f} dB" in capsys.readouterr().out


def test_cli_grid_includes_requested_experiment_defaults(tmp_path):
    assert main(["field-map", "--array", str(LINE8), "--target", "350", "800", "0",
                 "--outdir", str(tmp_path / "d")]) == 0
    z = np.load(tmp_path / "d/field-data.npz")
    assert z["focused"].shape == (191, 151)  # 151 x 191 points (x by y)
    assert (z["x_mm"][0], z["x_mm"][-1], z["y_mm"][0], z["y_mm"][-1]) == (-500, 1000, 100, 2000)


@pytest.mark.parametrize("extra", [["--nx", "0"], ["--ny", "-3"], ["--x-range", "5", "1"],
                                   ["--frequency", "0"], ["--gain", "0"], ["--range-db", "0"],
                                   ["--diff-limit", "-1"], ["--speed-of-sound", "0"],
                                   ["--x-range", "0", "1000", "--y-range", "0", "0", "--ny", "5"],
                                   ["--y-range", "-10", "10", "--x-range", "0", "0", "--nx", "1", "--ny", "3"]])
def test_cli_rejects_invalid_arguments(tmp_path, extra):
    with pytest.raises(SystemExit) as exc:
        run_cli(tmp_path, *extra)
    assert exc.value.code == 2


def test_cli_rejects_grid_point_on_speaker(tmp_path):
    with pytest.raises(SystemExit):
        main(["field-map", "--array", str(LINE8), "--target", "350", "800", "0",
              "--x-range", "0", "100", "--y-range", "-100", "100", "--nx", "3", "--ny", "3",
              "--outdir", str(tmp_path / "o")])


def test_level_map_rejects_empty_color_range(tmp_path):
    arr = load_array(LINE8)
    with pytest.raises(ValueError):
        save_level_map(tmp_path / "x.png", np.zeros((2, 2)), arr, (0, 0, 0), np.array([0.0, 1.0]),
                       np.array([0.0, 1.0]), 0.0, "t", 1.0, 1.0)


def test_other_commands_unchanged(tmp_path):
    assert main(["listen", "--array", str(LINE8), "--mode", "simultaneous", "--mic", "100", "500", "0",
                 "--duration", "0.2", "--gain", "0.05", "--wav", str(tmp_path / "m.wav")]) == 0
