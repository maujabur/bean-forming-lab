"""Steady-state complex pressure field of a speaker array (free field, point sources).

Each speaker emits A * sin(w (t - e_i)) with A = `gain` at the reference distance, so at
a point at distance r_i the contribution is A * (r_ref / r_i) * sin(w (t - e_i - r_i / c)).
The phasor P is defined by p(t) = |P| sin(w t + arg P), i.e.

    P(x) = sum_i  A * (r_ref / r_i) * exp(-j w (e_i + r_i / c))

This is exactly the steady state of microphone.receive for a sinusoidal source: |P| is the
amplitude and arg P the phase of the virtual-microphone signal. |P| / sqrt(2) is the RMS.

Levels are relative: 0 dB is one speaker heard at the reference distance (amplitude A).
Nothing here is a calibrated sound pressure level.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass

import numpy as np

from .array_config import ArrayConfig
from .focus import SPEED_OF_SOUND, focus_delays
from .microphone import MIN_DISTANCE, REFERENCE_DISTANCE

DEFAULT_FLOOR_DB = -60.0  # display floor, dB re one speaker at the reference distance


def complex_field(array: ArrayConfig, points: np.ndarray, frequency: float, *,
                  emission_delays: Mapping[int, float] | None = None,
                  speed_of_sound: float = SPEED_OF_SOUND, gain: float = 1.0,
                  reference_distance: float = REFERENCE_DISTANCE) -> np.ndarray:
    """Complex amplitude P at `points` (meters, shape (..., 3)); returns shape (...)."""
    points = np.asarray(points, dtype=np.float64)
    if points.ndim < 1 or points.shape[-1] != 3 or not np.all(np.isfinite(points)):
        raise ValueError("points must be finite XYZ coordinates, shape (..., 3)")
    if not math.isfinite(frequency) or frequency <= 0:
        raise ValueError("frequency must be finite and positive")
    if not math.isfinite(speed_of_sound) or speed_of_sound <= 0:
        raise ValueError("speed of sound must be finite and positive")
    if not math.isfinite(gain) or gain <= 0:
        raise ValueError("gain must be finite and positive")
    if not math.isfinite(reference_distance) or reference_distance <= 0:
        raise ValueError("reference distance must be finite and positive")
    if not array.elements:
        raise ValueError("array has no elements")
    emission_delays = emission_delays or {}
    known = {e.speaker for e in array.elements}
    if set(emission_delays) - known:
        raise ValueError("emission delays for speakers not in the array")
    omega = 2 * np.pi * frequency
    total = np.zeros(points.shape[:-1], dtype=np.complex128)
    for element in sorted(array.elements, key=lambda e: e.speaker):
        r = np.linalg.norm(points - np.asarray(element.position), axis=-1)
        if np.any(r < MIN_DISTANCE):
            bad = points[r < MIN_DISTANCE][0]
            raise ValueError(f"point {tuple(bad)} coincides with speaker {element.speaker}")
        delay = float(emission_delays.get(element.speaker, 0.0))
        if not math.isfinite(delay) or delay < 0:
            raise ValueError(f"emission delay of speaker {element.speaker} must be finite and >= 0")
        total += gain * (reference_distance / r) * np.exp(-1j * omega * (delay + r / speed_of_sound))
    return total


def grid_axis(lo: float, hi: float, n: int, label: str) -> np.ndarray:
    if type(n) is not int or n < 1:
        raise ValueError(f"{label}: point count must be a positive integer")
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo > hi:
        raise ValueError(f"{label}: range must be finite with min <= max")
    if n == 1:
        if lo != hi:
            raise ValueError(f"{label}: a single point needs min == max")
        return np.array([float(lo)])
    if lo == hi:
        raise ValueError(f"{label}: min must be below max for more than one point")
    return np.linspace(lo, hi, n)


@dataclass(frozen=True)
class FieldMaps:
    x: np.ndarray  # meters, (nx,)
    y: np.ndarray  # meters, (ny,)
    z: float  # meters
    simultaneous: np.ndarray  # complex, (ny, nx)
    focused: np.ndarray  # complex, (ny, nx)
    reference: float  # amplitude of 0 dB (= gain)


def field_maps(array: ArrayConfig, focus: tuple[float, float, float], frequency: float,
               x: np.ndarray, y: np.ndarray, z: float = 0.0, *,
               speed_of_sound: float = SPEED_OF_SOUND, gain: float = 1.0,
               reference_distance: float = REFERENCE_DISTANCE) -> FieldMaps:
    """Simultaneous and focused complex fields on the plane at height z (meters)."""
    emission = {d.speaker: d.delay for d in focus_delays(array, focus, speed_of_sound)}
    xx, yy = np.meshgrid(x, y)  # (ny, nx)
    points = np.stack([xx, yy, np.full_like(xx, z)], axis=-1)
    kw = dict(speed_of_sound=speed_of_sound, gain=gain, reference_distance=reference_distance)
    return FieldMaps(np.asarray(x), np.asarray(y), float(z),
                     complex_field(array, points, frequency, **kw),
                     complex_field(array, points, frequency, emission_delays=emission, **kw),
                     gain)


def level_db(p: np.ndarray, reference: float, floor_db: float = DEFAULT_FLOOR_DB) -> np.ndarray:
    """20 log10(|P| / reference), with |P| clamped to the floor so exact zeros stay finite.

    Only the dB array is clamped; keep the complex data for exact values.
    """
    if not reference > 0 or not math.isfinite(floor_db):
        raise ValueError("reference must be positive and floor finite")
    floor = reference * 10 ** (floor_db / 20)
    return 20 * np.log10(np.maximum(np.abs(p), floor) / reference)


def difference_db(focused: np.ndarray, simultaneous: np.ndarray, reference: float,
                  floor_db: float = DEFAULT_FLOOR_DB) -> np.ndarray:
    """Level(focused) - level(simultaneous), NaN where both are below the floor."""
    out = level_db(focused, reference, floor_db) - level_db(simultaneous, reference, floor_db)
    floor = reference * 10 ** (floor_db / 20)
    out[(np.abs(focused) < floor) & (np.abs(simultaneous) < floor)] = np.nan
    return out
