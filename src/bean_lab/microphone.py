"""Virtual omnidirectional microphone in free field.

p(t) = sum_i (r_ref / r_i) * s_i(t - delay_i - r_i / c)

Sources are identified by speaker number (never PCM channels). s_i is the pressure
waveform the speaker would produce at the reference distance r_ref (1 m by default).
Emission delay and propagation delay are summed and applied in ONE fractional-delay
operation per speaker, so no interpolation filter is cascaded.

Time reference: t = 0 is sample 0 of the emitted signals. The propagation delay is
never removed; the recording starts at t = 0 and its silence before the first arrival
is part of the result. Free field, point sources, no directivity, no reflections.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .array_config import ArrayConfig
from .focus import SPEED_OF_SOUND
from .fractional_delay import HALF_TAPS, delay_signal, output_length

REFERENCE_DISTANCE = 1.0  # meters
MIN_DISTANCE = 1e-9  # meters; below this the mic is considered coincident with a speaker


@dataclass(frozen=True)
class SourcePath:
    speaker: int
    distance: float  # m, speaker -> microphone
    travel_time: float  # s, distance / c
    emission_delay: float  # s, applied before emission
    gain: float  # r_ref / r

    @property
    def arrival_time(self) -> float:
        return self.emission_delay + self.travel_time


@dataclass(frozen=True)
class Reception:
    signal: np.ndarray  # mono, 1-D, t = 0 at first sample
    paths: tuple[SourcePath, ...]  # sorted by speaker


def propagation_paths(array: ArrayConfig, mic: tuple[float, float, float],
                      emission_delays: Mapping[int, float] | None = None,
                      speed_of_sound: float = SPEED_OF_SOUND,
                      reference_distance: float = REFERENCE_DISTANCE,
                      speakers: list[int] | None = None) -> tuple[SourcePath, ...]:
    """Geometry-only part: distances, travel times, 1/r gains (sorted by speaker)."""
    if len(mic) != 3 or any(type(v) not in (int, float) or not math.isfinite(v) for v in mic):
        raise ValueError("microphone position must contain three finite numbers")
    if not math.isfinite(speed_of_sound) or speed_of_sound <= 0:
        raise ValueError("speed of sound must be finite and positive")
    if not math.isfinite(reference_distance) or reference_distance <= 0:
        raise ValueError("reference distance must be finite and positive")
    emission_delays = emission_delays or {}
    wanted = sorted(e.speaker for e in array.elements) if speakers is None else sorted(speakers)
    paths = []
    for speaker in wanted:
        try:
            pos = array.position_for(speaker)
        except KeyError:
            raise ValueError(f"speaker {speaker} not in array") from None
        r = math.dist(mic, pos)
        if r < MIN_DISTANCE:
            raise ValueError(f"microphone coincides with speaker {speaker}")
        delay = float(emission_delays.get(speaker, 0.0))
        if not math.isfinite(delay) or delay < 0:
            raise ValueError(f"emission delay of speaker {speaker} must be finite and nonnegative")
        paths.append(SourcePath(speaker, r, r / speed_of_sound, delay, reference_distance / r))
    return tuple(paths)


def receive(sources: Mapping[int, np.ndarray], array: ArrayConfig,
            mic: tuple[float, float, float], sample_rate: int, *,
            emission_delays: Mapping[int, float] | None = None,
            speed_of_sound: float = SPEED_OF_SOUND,
            reference_distance: float = REFERENCE_DISTANCE) -> Reception:
    """Mono signal heard by a virtual microphone at `mic` (meters).

    sources: speaker -> emitted signal (arbitrary, 1-D, finite); speakers of the array
    that are absent are silent. emission_delays: optional extra per-speaker delay (s) to
    apply to the source, e.g. focusing delays; combined with propagation in one pass.
    Output length = max source length + ceil(largest total delay) + HALF_TAPS (tail kept).
    No normalization or clipping is applied.
    """
    if sample_rate < 1:
        raise ValueError("sample_rate must be positive")
    if not sources:
        raise ValueError("no sources given")
    emission_delays = emission_delays or {}
    unknown = set(emission_delays) - set(sources)
    if unknown:
        raise ValueError(f"emission delays for speakers without a source: {sorted(unknown)}")
    clean = {}
    for speaker, s in sources.items():
        s = np.asarray(s, dtype=np.float64)
        if s.ndim != 1 or len(s) == 0 or not np.all(np.isfinite(s)):
            raise ValueError(f"source of speaker {speaker} must be a nonempty finite 1-D signal")
        clean[speaker] = s
    paths = propagation_paths(array, mic, emission_delays, speed_of_sound,
                              reference_distance, speakers=list(clean))
    delays = {p.speaker: p.arrival_time * sample_rate for p in paths}
    frames = output_length(max(len(s) for s in clean.values()), max(delays.values()), HALF_TAPS)
    out = np.zeros(frames, dtype=np.float64)
    for p in paths:  # fixed (sorted) order => deterministic floating-point sum
        out += p.gain * delay_signal(clean[p.speaker], delays[p.speaker], frames)
    return Reception(out, paths)
