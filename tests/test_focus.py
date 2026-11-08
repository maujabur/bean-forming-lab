import json
import math
import wave
from pathlib import Path

import numpy as np
import pytest

from bean_lab.array_config import ArrayConfig, Element, load_array
from bean_lab.cli import main
from bean_lab.focus import SPEED_OF_SOUND, focus_delays, route_to_pcm, synthesize_speakers
from bean_lab.fractional_delay import apply_delays, delay_signal, output_length
from bean_lab.hardware_config import load_hardware

ROOT = Path(__file__).resolve().parents[1]
FS = 48000
C = SPEED_OF_SOUND


def make_array(*positions, speakers=None):
    speakers = speakers or range(1, len(positions) + 1)
    return ArrayConfig("t", tuple(Element(s, tuple(p)) for s, p in zip(speakers, positions)))


def estimate_delay_samples(y, freq, start, stop):
    """Delay (samples, within one period) of a sin(2*pi*freq*(m-D)/FS) via quadrature projection.

    The window [start, stop) should hold an integer number of periods.
    """
    m = np.arange(start, stop)
    w = 2 * np.pi * freq * m / FS
    s, c = np.sum(y[start:stop] * np.sin(w)), np.sum(y[start:stop] * np.cos(w))
    # y = sin(w - phi) = sin w cos phi - cos w sin phi
    phi = math.atan2(-c, s)
    return phi / (2 * np.pi * freq) * FS


# ---- geometry -------------------------------------------------------------

def test_known_distances():
    arr = make_array((0, 0, 0), (3, 0, 0))
    d = focus_delays(arr, (0, 4, 0))
    assert d[0].distance == pytest.approx(4.0)
    assert d[1].distance == pytest.approx(5.0)  # 3-4-5 triangle


def test_symmetric_array_central_focus_has_zero_delays():
    arr = make_array((-0.1, 0, 0), (0.1, 0, 0), (0, -0.1, 0), (0, 0.1, 0))
    d = focus_delays(arr, (0, 0, 1.0))
    assert all(x.delay == pytest.approx(0, abs=1e-15) for x in d)


def test_lateral_focus_delays_follow_formula():
    arr = make_array((0, 0, 0), (0.1, 0, 0))
    d = focus_delays(arr, (0.1, 0, 1.0))  # above speaker 2, which is the nearer one
    far, near = d[0], d[1]
    assert near.distance == pytest.approx(1.0)
    assert far.distance == pytest.approx(math.hypot(0.1, 1.0))
    assert far.delay == 0.0
    assert near.delay == pytest.approx((far.distance - 1.0) / C)
    assert near.delay > 0  # nearer speaker emits later


def test_delays_nonnegative_and_farthest_is_zero():
    arr = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    for focus in [(0.06, 0.06, 1.0), (0.3, -0.2, 0.5), (0.0, 0.0, 0.0), (-1.0, 0.7, 2.0)]:
        d = focus_delays(arr, focus)
        assert all(x.delay >= 0 for x in d)
        assert min(x.delay for x in d) == 0.0
        assert max(d, key=lambda x: x.distance).delay == 0.0


def test_arrival_times_aligned():
    arr = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    d = focus_delays(arr, (0.05, 0.02, 0.4))
    arrivals = [x.delay + x.distance / C for x in d]
    assert max(arrivals) - min(arrivals) < 1e-12


def test_independent_of_json_order(tmp_path):
    src = json.loads((ROOT / "configs/arrays/staggered-3-2-3.json").read_text())
    shuffled = dict(src, elements=list(reversed(src["elements"])))
    p = tmp_path / "rev.json"
    p.write_text(json.dumps(shuffled))
    focus = (0.06, 0.06, 1.0)
    assert focus_delays(load_array(p), focus) == focus_delays(
        load_array(ROOT / "configs/arrays/staggered-3-2-3.json"), focus)


def test_speed_of_sound_scales_delays():
    arr = make_array((0, 0, 0), (0.1, 0, 0))
    a = focus_delays(arr, (0.1, 0, 1.0), 343.0)[1].delay
    b = focus_delays(arr, (0.1, 0, 1.0), 686.0)[1].delay
    assert a == pytest.approx(2 * b)


@pytest.mark.parametrize("focus", [(0, 0), (0, 0, math.nan), (0, math.inf, 0), ("a", 0, 0)])
def test_invalid_focus(focus):
    with pytest.raises(ValueError):
        focus_delays(make_array((0, 0, 0)), focus)


@pytest.mark.parametrize("c", [0, -343, math.nan, math.inf])
def test_invalid_speed_of_sound(c):
    with pytest.raises(ValueError):
        focus_delays(make_array((0, 0, 0)), (0, 0, 1), c)


# ---- fractional delay -------------------------------------------------------

@pytest.mark.parametrize("freq", [500.0, 1000.0, 8000.0, 16000.0])
@pytest.mark.parametrize("delay", [0.0, 0.1, 0.25, 0.5, 0.77, 3.3, 40.5])
def test_fractional_delay_waveform_matches_analytic(freq, delay):
    n = 4800
    x = np.sin(2 * np.pi * freq * np.arange(n) / FS)
    y = delay_signal(x, delay, output_length(n, delay))
    m = np.arange(len(y))
    ref = np.sin(2 * np.pi * freq * (m - delay) / FS)
    steady = slice(int(delay) + 100, n - 100)  # away from onset/offset transients
    # Measured worst case ~1.3e-5 (Kaiser beta 8.6, 64 taps); 1e-4 leaves margin.
    assert np.max(np.abs(y[steady] - ref[steady])) < 1e-4


@pytest.mark.parametrize("delay", [0.05, 0.3, 0.5, 0.9, 7.25])
def test_fractional_delay_phase_precision(delay):
    freq, n = 1000.0, 4800  # period = 48 samples, so delays < 48 are unambiguous
    x = np.sin(2 * np.pi * freq * np.arange(n) / FS)
    y = delay_signal(x, delay, output_length(n, delay))
    est = estimate_delay_samples(y, freq, 480, 480 + 48 * 40)
    # 1e-3 sample = 21 ns = 7 um of air, ~1000x finer than one sample; measured error is far lower.
    assert est == pytest.approx(delay, abs=1e-3)


def test_integer_delay_is_exact_shift():
    x = np.random.default_rng(0).standard_normal(200)
    y = delay_signal(x, 5, output_length(200, 5))
    assert np.array_equal(y[5:205], x)
    assert np.all(y[:5] == 0) and np.all(y[205:] == 0)


def test_impulse_delay_centroid():
    x = np.zeros(200)
    x[50] = 1.0
    d = 10.37
    y = delay_signal(x, d, output_length(200, d))
    idx = np.arange(len(y))
    # Interpolated impulse: centroid sits at onset + delay. Sinc tails are odd-symmetric so
    # the first moment is accurate to ~1e-3 samples.
    assert np.sum(idx * y) / np.sum(y) == pytest.approx(50 + d, abs=5e-3)


def test_output_length_and_head_tail_policy():
    x = np.ones(1000)
    out = apply_delays(x, [0.0, 1e-3], FS)  # 48 samples
    assert out.shape == (output_length(1000, 48), 2)
    assert out.shape[0] == 1000 + 48 + 32
    assert np.all(out[:16, 1] == 0)  # silence before delayed onset
    assert np.allclose(out[48 + 100:1000 - 100, 1], 1.0, atol=1e-4)
    assert np.abs(out[-5:, :]).max() < 1e-6  # tail fully decayed, not truncated


def test_channels_keep_relative_timing():
    freq, n = 1000.0, 4800
    x = np.sin(2 * np.pi * freq * np.arange(n) / FS)
    delays = [0.0, 0.31 / FS, 2.8 / FS]
    out = apply_delays(x, delays, FS)
    ests = [estimate_delay_samples(out[:, i], freq, 480, 480 + 48 * 40) for i in range(3)]
    assert np.diff(ests) == pytest.approx(np.diff([d * FS for d in delays]), abs=1e-3)


@pytest.mark.parametrize("bad", [-1e-6, math.nan, math.inf])
def test_apply_delays_rejects_bad_delay(bad):
    with pytest.raises(ValueError):
        apply_delays(np.ones(10), [0.0, bad], FS)


def test_apply_delays_rejects_bad_signal():
    with pytest.raises(ValueError):
        apply_delays(np.ones((4, 2)), [0.0], FS)
    with pytest.raises(ValueError):
        apply_delays(np.array([]), [0.0], FS)
    with pytest.raises(ValueError):
        apply_delays(np.array([1.0, math.nan]), [0.0], FS)


# ---- synthesis and routing ---------------------------------------------------

def test_focused_signals_arrive_in_phase_at_focus():
    arr = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    freq, n = 1000.0, 4800
    focus = (0.06, 0.06, 1.0)
    delays = focus_delays(arr, focus)
    x = np.sin(2 * np.pi * freq * np.arange(n) / FS)
    sig = synthesize_speakers(x, delays, FS)
    # total delay in samples = emission delay + propagation time; must be equal for all speakers
    # (modulo one period, since the phase estimate is ambiguous beyond it).
    period = FS / freq
    totals = []
    for d in delays:
        emit = estimate_delay_samples(sig[d.speaker], freq, 480, 480 + 48 * 40) % period
        totals.append((emit + d.distance / C * FS) % period)
    spread = np.ptp(np.unwrap(np.array(totals) * 2 * np.pi / period) * period / (2 * np.pi))
    assert spread < 1e-3


def test_route_places_signals_on_mapped_pcm_channels():
    hw = load_hardware(ROOT / "configs/hardware/cmedia-0d8c-0102.json")
    sig = {s: np.full(10, 0.01 * s) for s in range(1, 9)}
    out = route_to_pcm(sig, hw)
    assert out.shape == (10, 8)
    for speaker, pcm in hw.speaker_to_pcm.items():
        assert np.all(out[:, pcm - 1] == pytest.approx(0.01 * speaker))


def test_route_leaves_unmapped_channels_silent():
    hw = load_hardware(ROOT / "configs/hardware/cmedia-0d8c-0102.json")
    out = route_to_pcm({1: np.full(4, 0.1)}, hw)  # speaker 1 -> PCM 2
    assert np.all(out[:, 1] == 0.1)
    assert np.count_nonzero(out) == 4


def test_route_rejects_unmapped_speaker_and_mismatched_lengths():
    hw = load_hardware(ROOT / "configs/hardware/cmedia-0d8c-0102.json")
    with pytest.raises(KeyError):
        route_to_pcm({99: np.zeros(4)}, hw)
    with pytest.raises(ValueError):
        route_to_pcm({1: np.zeros(4), 2: np.zeros(5)}, hw)
    with pytest.raises(ValueError):
        route_to_pcm({}, hw)


def test_route_refuses_to_clip_silently():
    hw = load_hardware(ROOT / "configs/hardware/cmedia-0d8c-0102.json")
    with pytest.raises(ValueError, match="full scale"):
        route_to_pcm({1: np.array([0.0, 1.2])}, hw)


def test_geometry_unaffected_by_hardware_mapping():
    # Same focus, same geometry => same delays regardless of which hardware map is used.
    arr = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    assert not hasattr(arr, "speaker_to_pcm")
    assert all(not hasattr(e, "pcm") for e in arr.elements)


# ---- CLI ---------------------------------------------------------------------

def run_focus(tmp_path, *extra):
    wav = tmp_path / "f.wav"
    args = ["focus", "--array", str(ROOT / "configs/arrays/staggered-3-2-3.json"),
            "--hardware", str(ROOT / "configs/hardware/cmedia-0d8c-0102.json"),
            "--target", "60", "60", "1000", "--frequency", "1000", "--duration", "0.5",
            "--gain", "0.05", "--wav", str(wav), *extra]
    return main(args), wav


def test_cli_focus_writes_valid_multichannel_wav(tmp_path, capsys):
    rc, wav = run_focus(tmp_path)
    assert rc == 0
    with wave.open(str(wav)) as w:
        assert (w.getnchannels(), w.getsampwidth(), w.getframerate()) == (8, 2, 48000)
        frames = w.getnframes()
        data = np.frombuffer(w.readframes(frames), dtype="<i2").reshape(frames, 8)
    assert frames == 24000 + 1 + 32  # tone + ceil(max delay 0.44 samples) + HALF_TAPS
    assert np.all(np.abs(data).max(axis=0) > 0)  # all eight channels carry signal
    assert np.abs(data).max() <= round(0.05 * 32767) + 5  # nominal gain: no clipping, no rescale
    out = capsys.readouterr().out
    assert "distance_mm" in out and "1003.594" in out


def test_cli_focus_channel_mapping(tmp_path):
    # Use a custom hardware map with only speaker 4 mapped from the array subset.
    arr = tmp_path / "a.json"
    arr.write_text(json.dumps({"name": "x", "units": "mm", "elements": [
        {"speaker": 4, "position": [0, 0, 0]}]}))
    hw = tmp_path / "h.json"
    hw.write_text(json.dumps({"name": "h", "device": "d", "channels": 8,
                              "speaker_to_pcm": {"4": 7}}))
    wav = tmp_path / "o.wav"
    assert main(["focus", "--array", str(arr), "--hardware", str(hw), "--target", "0", "0", "500",
                 "--duration", "0.1", "--wav", str(wav)]) == 0
    with wave.open(str(wav)) as w:
        n = w.getnframes()
        data = np.frombuffer(w.readframes(n), dtype="<i2").reshape(n, 8)
    assert np.count_nonzero(data[:, 6]) > 0
    assert np.count_nonzero(np.delete(data, 6, axis=1)) == 0


@pytest.mark.parametrize("extra", [["--gain", "1.5"], ["--gain", "0"], ["--frequency", "30000"],
                                   ["--duration", "0"], ["--speed-of-sound", "0"]])
def test_cli_focus_rejects_invalid_arguments(tmp_path, extra):
    with pytest.raises(SystemExit) as exc:
        run_focus(tmp_path, *extra)
    assert exc.value.code == 2


def test_cli_focus_rejects_missing_files_and_bad_target(tmp_path):
    with pytest.raises(SystemExit):
        main(["focus", "--array", str(tmp_path / "nope.json"), "--target", "0", "0", "1"])
    with pytest.raises(SystemExit):
        main(["focus", "--array", str(ROOT / "configs/arrays/linear-4.json"), "--target", "0", "0"])
    with pytest.raises(SystemExit):
        main(["focus", "--array", str(ROOT / "configs/arrays/linear-4.json"),
              "--target", "0", "nan", "1000"])


def test_existing_commands_still_work(tmp_path):
    wav = tmp_path / "t.wav"
    assert main(["tone", "--speaker", "1", "--hardware",
                 str(ROOT / "configs/hardware/cmedia-0d8c-0102.json"),
                 "--seconds", "0.1", "--output", str(wav)]) == 0
    assert main(["pulse", "--channel", "3", "--seconds", "0.5", "--output", str(wav)]) == 0
