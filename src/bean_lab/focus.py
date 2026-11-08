"""Static acoustic focusing: geometry -> delays -> per-speaker signals -> PCM routing.

The geometry stage knows speaker identities only, never PCM channels.
Coordinates are XYZ meters (as normalized by array_config).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .array_config import ArrayConfig
from .fractional_delay import apply_delays
from .hardware_config import HardwareConfig

SPEED_OF_SOUND = 343.0  # m/s


@dataclass(frozen=True)
class SpeakerDelay:
    speaker: int
    distance: float  # meters, to the focal point
    delay: float  # seconds, nonnegative; the farthest speaker has exactly 0


def focus_delays(array: ArrayConfig, focus: tuple[float, float, float],
                 speed_of_sound: float = SPEED_OF_SOUND) -> tuple[SpeakerDelay, ...]:
    """Delays making all wavefronts arrive at `focus` simultaneously.

    d_i = |focus - p_i|; delay_i = (d_max - d_i) / c. Result is sorted by speaker
    number, hence independent of element order in the JSON.
    """
    if len(focus) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) for v in focus):
        raise ValueError("focus must contain three finite numbers")
    if not math.isfinite(speed_of_sound) or speed_of_sound <= 0:
        raise ValueError("speed of sound must be finite and positive")
    if not array.elements:
        raise ValueError("array has no elements")
    dist = {e.speaker: math.dist(focus, e.position) for e in array.elements}
    d_max = max(dist.values())
    return tuple(SpeakerDelay(s, dist[s], (d_max - dist[s]) / speed_of_sound)
                 for s in sorted(dist))


def synthesize_speakers(mono: np.ndarray, delays: tuple[SpeakerDelay, ...],
                        sample_rate: int) -> dict[int, np.ndarray]:
    """Delay `mono` for each speaker (fractional delays, no gain change).

    All signals have equal length: len(mono) + ceil(max delay) + HALF_TAPS samples;
    silence before the delayed onset, nothing truncated at the tail.
    """
    if not delays:
        raise ValueError("no speaker delays given")
    out = apply_delays(mono, [d.delay for d in delays], sample_rate)
    return {d.speaker: out[:, i] for i, d in enumerate(delays)}


def route_to_pcm(signals: dict[int, np.ndarray], hardware: HardwareConfig) -> np.ndarray:
    """Place each speaker signal on its PCM channel; unused channels stay silent.

    Raises ValueError if any sample magnitude exceeds 1.0 (no silent clipping or
    normalization: relative channel levels are never altered).
    """
    if not signals:
        raise ValueError("no signals to route")
    frames = {len(s) for s in signals.values()}
    if len(frames) != 1:
        raise ValueError("all speaker signals must have the same length")
    out = np.zeros((frames.pop(), hardware.channels), dtype=np.float64)
    for speaker, signal in signals.items():
        out[:, hardware.pcm_for(speaker) - 1] = signal  # KeyError if unmapped
    peak = float(np.max(np.abs(out)))
    if peak > 1.0:
        raise ValueError(f"peak {peak:.4f} exceeds full scale; lower the gain")
    return out
