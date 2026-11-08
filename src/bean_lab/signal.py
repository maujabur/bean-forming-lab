"""Deterministic multichannel signal generation.

Channels here are PCM frame indices, NOT speaker-test -s channel numbers.
All output is interleaved signed 16-bit PCM, at 48 kHz by default.
"""

from __future__ import annotations

import math
import wave
from pathlib import Path

import numpy as np


def tone(
    *,
    channels: int = 8,
    sample_rate: int = 48000,
    duration: float = 2.0,
    frequency: float = 1000.0,
    active: tuple[int, ...] = (0,),
    gain: float = 0.15,
    delays_seconds: tuple[float, ...] | None = None,
) -> np.ndarray:
    """Generate an (N, channels) float array.

    delays_seconds are independent per-channel positive delays in seconds.
    Fractional delays are computed analytically for this single sinusoid.
    Channel numbers in active are zero-based PCM indices.
    """
    if channels < 1 or sample_rate < 1 or duration <= 0:
        raise ValueError("channels, sample_rate and duration must be positive")
    if not 0 < frequency < sample_rate / 2:
        raise ValueError("frequency must be below Nyquist")
    if not 0 < gain <= 1:
        raise ValueError("gain must be in (0, 1]")
    if any(c < 0 or c >= channels for c in active):
        raise ValueError("active channel outside requested channel count")
    if delays_seconds is None:
        delays_seconds = (0.0,) * channels
    if len(delays_seconds) != channels or any(not math.isfinite(x) or x < 0 for x in delays_seconds):
        raise ValueError("delays_seconds must contain one nonnegative finite value per channel")

    frames = round(duration * sample_rate)
    out = np.zeros((frames, channels), dtype=np.float64)
    t = np.arange(frames, dtype=np.float64) / sample_rate
    # Gate the beginning of each delayed tone, to avoid pre-echo.
    for c in active:
        shifted = t - delays_seconds[c]
        out[:, c] = np.where(shifted >= 0, gain * np.sin(2 * np.pi * frequency * shifted), 0)
    return out


def pulse(
    *,
    channels: int = 8,
    sample_rate: int = 48000,
    duration: float = 1.0,
    active: tuple[int, ...] = (0,),
    amplitude: float = 0.15,
    position: float = 0.1,
) -> np.ndarray:
    """Single sample impulse on selected PCM channel indices."""
    if channels < 1 or sample_rate < 1 or duration <= 0:
        raise ValueError("channels, sample_rate and duration must be positive")
    if not 0 < amplitude <= 1 or not 0 <= position < duration:
        raise ValueError("invalid amplitude or pulse position")
    if any(c < 0 or c >= channels for c in active):
        raise ValueError("active channel outside requested channel count")
    out = np.zeros((round(duration * sample_rate), channels), dtype=np.float64)
    out[round(position * sample_rate), list(active)] = amplitude
    return out


def to_s16le(samples: np.ndarray) -> bytes:
    """Encode floating-point PCM with clipping protection."""
    if samples.ndim != 2 or not np.all(np.isfinite(samples)):
        raise ValueError("expected finite (frames, channels) samples")
    return np.rint(np.clip(samples, -1.0, 1.0) * 32767).astype("<i2").tobytes()


def write_wav(path: str | Path, samples: np.ndarray, sample_rate: int = 48000) -> None:
    with wave.open(str(path), "wb") as wav:
        wav.setnchannels(samples.shape[1])
        wav.setsampwidth(2)
        wav.setframerate(sample_rate)
        wav.writeframes(to_s16le(samples))
