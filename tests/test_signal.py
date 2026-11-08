import numpy as np
import pytest

from bean_lab.signal import pulse, to_s16le, tone


def test_single_tone_channel():
    data = tone(channels=8, sample_rate=48000, duration=0.1, active=(7,), frequency=1000)
    assert data.shape == (4800, 8)
    assert np.max(np.abs(data[:, :7])) == 0
    assert np.max(np.abs(data[:, 7])) > 0


def test_fractional_tone_delay():
    data = tone(channels=2, sample_rate=48000, duration=0.01, active=(0, 1),
                frequency=1000, delays_seconds=(0, 0.00025))
    assert np.all(data[:12, 1] == 0)
    assert np.max(np.abs(data[:, 1])) > 0
    assert not np.allclose(data[:, 0], data[:, 1])


def test_pulse():
    data = pulse(channels=8, active=(0, 7))
    assert np.count_nonzero(data) == 2
    assert data.shape[1] == 8


def test_int16_interleaving():
    pcm = np.array([[0.0, 1.0], [-1.0, 0.0]])
    assert np.frombuffer(to_s16le(pcm), dtype="<i2").tolist() == [0, 32767, -32767, 0]


def test_nyquist_guard():
    with pytest.raises(ValueError):
        tone(frequency=24000)
