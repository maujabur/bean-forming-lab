"""Numerical measurements on received signals."""

from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np

from .fractional_delay import HALF_TAPS
from .signal import FADE_SECONDS


def check_full_scale(samples: np.ndarray) -> float:
    """Return the absolute peak; raise if it exceeds full scale (no silent clipping)."""
    peak = peak_abs(samples)
    if peak > 1.0:
        raise ValueError(f"peak {peak:.4f} exceeds full scale; lower the gain")
    return peak


def peak_abs(x: np.ndarray) -> float:
    return float(np.max(np.abs(x)))


def rms(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0:
        raise ValueError("empty window")
    return float(np.sqrt(np.mean(x * x)))


def db_ratio(value: float, reference: float) -> float:
    """20*log10(value/reference); both must be positive (otherwise no valid reference)."""
    if not (value > 0 and reference > 0):
        raise ValueError("dB ratio needs two positive amplitudes")
    return 20.0 * math.log10(value / reference)


def steady_window(arrival_times: Sequence[Sequence[float]], duration: float, sample_rate: int,
                  frequency: float | None = None,
                  fade_seconds: float = FADE_SECONDS) -> tuple[int, int]:
    """Common sample window [i0, i1) in which every source of every case is fully on.

    arrival_times: for each case, the arrival time (emission + propagation) of each speaker.
    A source of duration `duration` with fades is fully on from arrival + margin to
    arrival + duration - margin, margin = fade + interpolation half-length. The window is
    the intersection over all sources of all cases. If `frequency` is given the length is
    trimmed to a whole number of periods, so RMS of a sine is not biased by a partial cycle.
    """
    margin = fade_seconds + HALF_TAPS / sample_rate
    flat = [a for case in arrival_times for a in case]
    if not flat:
        raise ValueError("no arrival times")
    t0 = max(flat) + margin
    t1 = min(flat) + duration - margin
    i0 = math.ceil(t0 * sample_rate)
    i1 = math.floor(t1 * sample_rate)
    if i1 - i0 < 1:
        raise ValueError("no common steady-state window; increase --duration")
    if frequency is not None:
        cycles = math.floor((i1 - i0) / sample_rate * frequency)
        if cycles < 1:
            raise ValueError("steady-state window shorter than one period")
        i1 = i0 + round(cycles * sample_rate / frequency)
    return i0, i1
