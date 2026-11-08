import json
import math
import wave
from pathlib import Path

import numpy as np
import pytest

import bean_lab.microphone as microphone_module
from bean_lab.analysis import check_full_scale, db_ratio, peak_abs, rms, steady_window
from bean_lab.array_config import ArrayConfig, Element, load_array
from bean_lab.cli import main
from bean_lab.focus import SPEED_OF_SOUND, focus_delays, synthesize_speakers
from bean_lab.fractional_delay import HALF_TAPS, apply_delays
from bean_lab.microphone import propagation_paths, receive
from bean_lab.signal import faded_tone

ROOT = Path(__file__).resolve().parents[1]
STAGGERED = ROOT / "configs/arrays/staggered-3-2-3.json"
FS = 48000
C = SPEED_OF_SOUND


def make_array(*positions, speakers=None):
    speakers = speakers or range(1, len(positions) + 1)
    return ArrayConfig("t", tuple(Element(s, tuple(map(float, p))) for s, p in zip(speakers, positions)))


def sine(freq, seconds=1.0, amp=1.0):
    return amp * np.sin(2 * np.pi * freq * np.arange(round(seconds * FS)) / FS)


def phase_delay_samples(y, freq, start, stop):
    """Delay of sin(2*pi*f*(m-D)/FS) in samples, unambiguous for D < one period.

    Quadrature projection over [start, stop), which must hold whole periods.
    """
    w = 2 * np.pi * freq * np.arange(start, stop) / FS
    s, c = np.sum(y[start:stop] * np.sin(w)), np.sum(y[start:stop] * np.cos(w))
    return math.atan2(-c, s) / (2 * np.pi * freq) * FS


# ---- single source: delay and 1/r (analytic expectations) ---------------------

def test_single_impulse_delay_and_gain():
    arr = make_array((0, 0, 0))
    x = np.zeros(300)
    x[100] = 1.0
    rec = receive({1: x}, arr, (3.0, 4.0, 0.0), FS)  # 3-4-5 triangle: r = 5 m
    assert rec.paths[0].distance == pytest.approx(5.0)
    expected = 100 + 5.0 / C * FS  # 100 + 699.708... samples
    idx = np.arange(len(rec.signal))
    # Windowed-sinc impulse: sum is the DC gain (1/r = 0.2), centroid gives the arrival time.
    assert np.sum(rec.signal) == pytest.approx(0.2, rel=1e-3)
    assert np.sum(idx * rec.signal) / np.sum(rec.signal) == pytest.approx(expected, abs=5e-3)


def test_propagation_delay_is_not_removed_and_tail_is_kept():
    arr = make_array((0, 0, 0))
    x = np.zeros(50)
    x[0] = 1.0
    rec = receive({1: x}, arr, (0.0, 0.0, 1.0), FS)
    d = 1.0 / C * FS  # 139.9 samples
    assert len(rec.signal) == 50 + math.ceil(d) + HALF_TAPS
    first = np.flatnonzero(rec.signal)[0]
    assert first >= int(d) - HALF_TAPS  # nothing before the earliest sinc tap
    assert np.all(rec.signal[:first] == 0)
    assert np.abs(rec.signal[-3:]).max() < 1e-3  # ring-down fully inside the output


@pytest.mark.parametrize("r", [0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
def test_one_over_r_attenuation(r):
    arr = make_array((0, 0, 0))
    rec = receive({1: sine(1000)}, arr, (0.0, r, 0.0), FS)
    # Unit sine, 1 m reference: steady-state RMS is (1/r)/sqrt(2).
    assert rms(rec.signal[2000:40000]) == pytest.approx(1.0 / r / math.sqrt(2), rel=2e-4)


def test_reference_distance_scales_gain():
    arr = make_array((0, 0, 0))
    a = receive({1: sine(1000)}, arr, (0, 2.0, 0), FS, reference_distance=1.0)
    b = receive({1: sine(1000)}, arr, (0, 2.0, 0), FS, reference_distance=0.5)
    assert rms(b.signal[2000:40000]) == pytest.approx(0.5 * rms(a.signal[2000:40000]), rel=1e-6)


@pytest.mark.parametrize("r", [0.0343, 0.1, 0.1234])
def test_phase_matches_travel_time(r):
    arr = make_array((0, 0, 0))
    rec = receive({1: sine(1000)}, arr, (r, 0.0, 0.0), FS)
    est = phase_delay_samples(rec.signal, 1000, 480, 480 + 48 * 40)
    # Delay is under one period (48 samples) so the phase is unambiguous; 1e-3 sample = 7 um of air.
    assert est == pytest.approx(r / C * FS, abs=1e-3)


def test_emission_and_propagation_delays_add():
    arr = make_array((0, 0, 0))
    rec = receive({1: sine(1000)}, arr, (0.1, 0, 0), FS, emission_delays={1: 0.3 / FS})
    est = phase_delay_samples(rec.signal, 1000, 480, 480 + 48 * 40)
    assert est == pytest.approx(0.1 / C * FS + 0.3, abs=1e-3)


# ---- sums of sources -----------------------------------------------------------

def test_coherent_sum_of_symmetric_sources():
    arr = make_array((-0.05, 0, 0), (0.05, 0, 0))
    rec = receive({1: sine(1000), 2: sine(1000)}, arr, (0, 0, 1.0), FS)
    r = math.hypot(0.05, 1.0)
    assert rms(rec.signal[2000:40000]) == pytest.approx(2 / r / math.sqrt(2), rel=2e-4)


def test_half_wavelength_difference_cancels_to_analytic_residual():
    lam = C / 1000
    arr = make_array((0, 0, 0), (0, 0, lam / 2))
    rec = receive({1: sine(1000), 2: sine(1000)}, arr, (0, 0, -1.0), FS)  # r1 = 1, r2 = 1 + lam/2
    expected = abs(1 / 1.0 - 1 / (1.0 + lam / 2)) / math.sqrt(2)
    assert rms(rec.signal[2000:40000]) == pytest.approx(expected, rel=1e-3)


def test_receive_is_linear_for_arbitrary_signals():
    arr = load_array(STAGGERED)
    rng = np.random.default_rng(1)
    a, b = rng.standard_normal(500) * 0.01, rng.standard_normal(500) * 0.01
    mic = (0.07, 0.01, 0.8)
    both = receive({1: a, 5: b}, arr, mic, FS).signal
    only_a = receive({1: a}, arr, mic, FS).signal
    only_b = receive({5: b}, arr, mic, FS).signal
    n = max(len(only_a), len(only_b), len(both))
    pad = lambda s: np.pad(s, (0, n - len(s)))
    assert np.allclose(pad(both), pad(only_a) + pad(only_b), atol=1e-12)


# ---- focusing at and away from the target -------------------------------------------

def two_speaker_case():
    lam = C / 1000
    arr = make_array((0, 0, 0), (lam / 2, 0, 0))  # x-axis, half a wavelength apart
    target = (-1.0, 0.0, 0.0)  # r1 = 1, r2 = 1 + lam/2
    delays = {d.speaker: d.delay for d in focus_delays(arr, target)}
    return arr, target, delays, lam


def test_focusing_adds_coherently_at_target_where_simultaneous_cancels():
    arr, target, delays, lam = two_speaker_case()
    sources = {1: sine(1000), 2: sine(1000)}
    simult = receive(sources, arr, target, FS)
    focused = receive(sources, arr, target, FS, emission_delays=delays)
    g_sum = 1.0 + 1.0 / (1.0 + lam / 2)
    g_diff = 1.0 - 1.0 / (1.0 + lam / 2)
    assert rms(focused.signal[3000:40000]) == pytest.approx(g_sum / math.sqrt(2), rel=1e-3)
    assert rms(simult.signal[3000:40000]) == pytest.approx(g_diff / math.sqrt(2), rel=1e-3)
    assert rms(focused.signal[3000:40000]) > 10 * rms(simult.signal[3000:40000])


def test_off_target_level_never_exceeds_coherent_sum_and_differs_from_target():
    arr, target, delays, lam = two_speaker_case()
    sources = {1: sine(1000), 2: sine(1000)}
    mic = (0.0, 1.0, 0.0)  # perpendicular to the array axis
    off = receive(sources, arr, mic, FS, emission_delays=delays)
    on = receive(sources, arr, target, FS, emission_delays=delays)
    bound = sum(p.gain for p in off.paths) / math.sqrt(2)  # triangle inequality
    level = rms(off.signal[3000:40000])
    assert level <= bound * (1 + 1e-6)
    assert level < 0.8 * rms(on.signal[3000:40000])  # focusing is NOT isotropic


def test_staggered_array_8khz_focus_raises_level_at_target():
    arr = load_array(STAGGERED)
    x = faded_tone(sample_rate=FS, duration=1.0, frequency=8000, gain=0.05)
    tgt = (0.06, 0.06, 1.0)
    em = {d.speaker: d.delay for d in focus_delays(arr, tgt)}
    src = {e.speaker: x for e in arr.elements}
    a = receive(src, arr, tgt, FS)
    b = receive(src, arr, tgt, FS, emission_delays=em)
    i0, i1 = steady_window([[p.arrival_time for p in r.paths] for r in (a, b)], 1.0, FS, 8000)
    assert rms(b.signal[i0:i1]) > rms(a.signal[i0:i1])


def test_single_pass_delay_is_at_least_as_accurate_as_two_stage():
    arr = load_array(STAGGERED)
    f = 16000.0
    x = sine(f, 1.0, 0.05)
    delays = focus_delays(arr, (0.06, 0.06, 1.0))
    em = {d.speaker: d.delay for d in delays}
    mic = (0.26, 0.06, 1.0)
    one = receive({d.speaker: x for d in delays}, arr, mic, FS, emission_delays=em)
    two = receive(synthesize_speakers(x, delays, FS), arr, mic, FS)
    m = np.arange(len(one.signal))
    ref = sum(p.gain * 0.05 * np.sin(2 * np.pi * f * (m - p.arrival_time * FS) / FS) for p in one.paths)
    s = slice(500, 40000)
    err_one = np.max(np.abs(one.signal[s] - ref[s]))
    err_two = np.max(np.abs(two.signal[s] - ref[s]))
    assert err_one < 1e-5 and err_two < 1e-5  # measured ~2e-6 / ~3e-6 (amplitude 0.05 per speaker)
    assert err_one <= err_two


# ---- ordering, independence from PCM ---------------------------------------------------

def test_independent_of_json_and_source_order(tmp_path):
    src = json.loads(STAGGERED.read_text())
    p = tmp_path / "rev.json"
    p.write_text(json.dumps(dict(src, elements=list(reversed(src["elements"])))))
    a, b = load_array(STAGGERED), load_array(p)
    x = faded_tone(sample_rate=FS, duration=0.1, frequency=1000, gain=0.05)
    mic = (0.1, 0.0, 0.7)
    s1 = {e.speaker: x for e in a.elements}
    s2 = dict(reversed(list(s1.items())))
    em = {d.speaker: d.delay for d in focus_delays(a, (0.06, 0.06, 1.0))}
    r1 = receive(s1, a, mic, FS, emission_delays=em)
    r2 = receive(s2, b, mic, FS, emission_delays=dict(reversed(list(em.items()))))
    assert np.array_equal(r1.signal, r2.signal)
    assert r1.paths == r2.paths


def test_microphone_has_no_hardware_dependency(tmp_path):
    assert "hardware_config" not in vars(microphone_module)
    assert not hasattr(microphone_module, "HardwareConfig")
    # the CLI needs no hardware file and refuses a --hardware option
    wav = tmp_path / "m.wav"
    assert main(["listen", "--array", str(STAGGERED), "--mode", "simultaneous",
                 "--mic", "60", "60", "1000", "--duration", "0.2", "--wav", str(wav)]) == 0
    with pytest.raises(SystemExit):
        main(["listen", "--array", str(STAGGERED), "--hardware", "x.json", "--mic", "0", "0", "1"])


# ---- analysis helpers ---------------------------------------------------------------------

def test_steady_window_is_inside_all_sources_and_whole_periods():
    arrivals = [[0.0030, 0.0032], [0.0031, 0.0034]]
    i0, i1 = steady_window(arrivals, 1.0, FS, 1000.0)
    margin = 0.005 + HALF_TAPS / FS
    assert i0 / FS >= 0.0034 + margin
    assert i1 / FS <= 0.0030 + 1.0 - margin + 1e-9
    assert (i1 - i0) % 48 == 0  # whole 1 kHz periods (48 samples)


def test_steady_window_rejects_too_short_duration():
    with pytest.raises(ValueError):
        steady_window([[0.003]], 0.005, FS, 1000.0)


def test_rms_of_sine_over_whole_periods_and_db():
    assert rms(sine(1000, 1.0, 0.5)) == pytest.approx(0.5 / math.sqrt(2), rel=1e-9)
    assert peak_abs(np.array([0.1, -0.7])) == 0.7
    assert db_ratio(2.0, 1.0) == pytest.approx(6.0206, abs=1e-3)
    with pytest.raises(ValueError):
        db_ratio(0.0, 1.0)
    with pytest.raises(ValueError):
        rms(np.array([]))


def test_check_full_scale():
    assert check_full_scale(np.array([0.5, -1.0])) == 1.0
    with pytest.raises(ValueError, match="full scale"):
        check_full_scale(np.array([0.5, 1.0001]))


# ---- invalid inputs -------------------------------------------------------------------------

def test_rejects_microphone_on_a_speaker():
    arr = make_array((0.1, 0.2, 0.3))
    with pytest.raises(ValueError, match="coincides"):
        receive({1: np.ones(4)}, arr, (0.1, 0.2, 0.3), FS)
    # but a point 1 mm away is a valid (very loud) position, there is no artificial minimum
    assert propagation_paths(arr, (0.1, 0.2, 0.301))[0].distance == pytest.approx(0.001)


@pytest.mark.parametrize("mic", [(0, 0), (0, 0, 0, 0), (0, math.nan, 1), (math.inf, 0, 1), ("a", 0, 1)])
def test_rejects_bad_mic(mic):
    with pytest.raises(ValueError):
        receive({1: np.ones(4)}, make_array((0, 0, 0)), mic, FS)


@pytest.mark.parametrize("kwargs", [{"speed_of_sound": 0}, {"speed_of_sound": math.nan},
                                    {"reference_distance": 0}, {"reference_distance": -1},
                                    {"emission_delays": {1: -1e-3}},
                                    {"emission_delays": {1: math.nan}},
                                    {"emission_delays": {9: 0.0}}])
def test_rejects_bad_parameters(kwargs):
    with pytest.raises(ValueError):
        receive({1: np.ones(4)}, make_array((0, 0, 0)), (0, 0, 1), FS, **kwargs)


def test_rejects_bad_sources():
    arr = make_array((0, 0, 0))
    for sources in ({}, {1: np.array([])}, {1: np.ones((2, 2))}, {1: np.array([1.0, math.nan])},
                    {7: np.ones(4)}):
        with pytest.raises(ValueError):
            receive(sources, arr, (0, 0, 1), FS)
    with pytest.raises(ValueError):
        receive({1: np.ones(4)}, arr, (0, 0, 1), 0)


# ---- CLI: WAV integrity, clipping, experiments ----------------------------------------------------

def read_mono(path):
    with wave.open(str(path)) as w:
        assert (w.getnchannels(), w.getsampwidth(), w.getframerate()) == (1, 2, 48000)
        return np.frombuffer(w.readframes(w.getnframes()), dtype="<i2")


def listen(tmp_path, *extra, wav="m.wav"):
    out = tmp_path / wav
    rc = main(["listen", "--array", str(STAGGERED), "--target", "60", "60", "1000",
               "--frequency", "1000", "--duration", "1", "--gain", "0.05", *extra, "--wav", str(out)])
    return rc, out


def test_cli_listen_wav_matches_model(tmp_path, capsys):
    rc, out = listen(tmp_path, "--mic", "60", "60", "1000")
    assert rc == 0
    pcm = read_mono(out)
    arr = load_array(STAGGERED)
    x = faded_tone(sample_rate=FS, duration=1.0, frequency=1000, gain=0.05)
    em = {d.speaker: d.delay for d in focus_delays(arr, (0.06, 0.06, 1.0))}
    ref = receive({e.speaker: x for e in arr.elements}, arr, (0.06, 0.06, 1.0), FS, emission_delays=em)
    longest = max(p.arrival_time for p in ref.paths) * FS
    assert len(pcm) == len(ref.signal) == 48000 + math.ceil(longest) + HALF_TAPS
    assert np.max(np.abs(pcm / 32767 - ref.signal)) <= 0.5 / 32767 + 1e-12  # int16 quantization only
    text = capsys.readouterr().out
    assert "distance_mm" in text and "dB re one speaker at 1 m" in text and "1000.450" in text


def test_cli_listen_modes_differ_and_focused_requires_target(tmp_path):
    _, a = listen(tmp_path, "--mode", "simultaneous", "--mic", "60", "60", "1000", wav="a.wav")
    _, b = listen(tmp_path, "--mode", "focused", "--mic", "60", "60", "1000", wav="b.wav")
    assert not np.array_equal(read_mono(a), read_mono(b))
    with pytest.raises(SystemExit):
        main(["listen", "--array", str(STAGGERED), "--mode", "focused", "--mic", "0", "0", "1000"])


def test_cli_listen_refuses_to_clip_and_writes_nothing(tmp_path):
    out = tmp_path / "clip.wav"
    with pytest.raises(SystemExit) as exc:  # 8 speakers at ~1 m with gain 0.2 sum above full scale
        main(["listen", "--array", str(STAGGERED), "--mode", "simultaneous", "--mic", "60", "60", "300",
              "--gain", "0.5", "--wav", str(out)])
    assert exc.value.code == 2
    assert not out.exists()


@pytest.mark.parametrize("extra", [["--gain", "0"], ["--gain", "2"], ["--frequency", "30000"],
                                   ["--duration", "0"], ["--speed-of-sound", "-1"],
                                   ["--window", "5", "6"], ["--mic", "0", "0", "0"]])
def test_cli_listen_rejects_invalid_arguments(tmp_path, extra):
    args = ["listen", "--array", str(STAGGERED), "--mode", "simultaneous", "--mic", "60", "60", "1000",
            "--duration", "0.5"]
    # later options override earlier ones in argparse
    with pytest.raises(SystemExit) as exc:
        main(args + extra)
    assert exc.value.code == 2


def compare(tmp_path, *extra):
    return main(["compare", "--array", str(STAGGERED), "--target", "60", "60", "1000",
                 "--off-target", "260", "60", "1000", "--frequency", "1000", "--duration", "2",
                 "--gain", "0.05", "--outdir", str(tmp_path / "exp"), *extra])


def test_cli_compare_writes_four_wavs_without_normalization(tmp_path, capsys):
    assert compare(tmp_path) == 0
    levels = {}
    for k in "ABCD":
        pcm = read_mono(tmp_path / "exp" / f"exp-{k}.wav")
        assert np.max(np.abs(pcm)) < 32767
        levels[k] = pcm
    text = capsys.readouterr().out
    for k in "ABCD":
        assert f"exp-{k}.wav" in text
    assert "B - A" in text and "D - C" in text
    # Same absolute scale: peaks reproduce the model, not a per-file normalization.
    arr = load_array(STAGGERED)
    x = faded_tone(sample_rate=FS, duration=2.0, frequency=1000, gain=0.05)
    src = {e.speaker: x for e in arr.elements}
    ref = receive(src, arr, (0.06, 0.06, 1.0), FS).signal
    assert np.max(np.abs(levels["A"] / 32767 - ref)) <= 0.5 / 32767 + 1e-12


def test_cli_compare_listen_reports_same_rms_with_common_window(tmp_path, capsys):
    compare(tmp_path)
    table = capsys.readouterr().out
    window = next(l for l in table.splitlines() if l.startswith("Common window"))
    t0, t1 = window.split()[2], window.split()[4]
    row_d = next(l for l in table.splitlines() if l.strip().startswith("D "))
    rc = main(["listen", "--array", str(STAGGERED), "--mode", "focused", "--target", "60", "60", "1000",
               "--mic", "260", "60", "1000", "--frequency", "1000", "--duration", "2", "--gain", "0.05",
               "--window", t0, t1])
    assert rc == 0
    assert row_d.split()[3] in capsys.readouterr().out  # identical 6-decimal RMS


def test_cli_compare_refuses_to_clip_and_writes_no_files(tmp_path):
    with pytest.raises(SystemExit) as exc:
        main(["compare", "--array", str(STAGGERED), "--target", "60", "60", "300",
              "--off-target", "260", "60", "1000", "--gain", "0.5", "--outdir", str(tmp_path / "exp")])
    assert exc.value.code == 2
    assert not (tmp_path / "exp").exists()


def test_existing_commands_unchanged(tmp_path):
    assert main(["tone", "--channel", "1", "--seconds", "0.1", "--output", str(tmp_path / "t.wav")]) == 0
    assert main(["focus", "--array", str(STAGGERED), "--hardware",
                 str(ROOT / "configs/hardware/cmedia-0d8c-0102.json"), "--target", "60", "60", "1000",
                 "--duration", "0.1", "--wav", str(tmp_path / "f.wav")]) == 0
