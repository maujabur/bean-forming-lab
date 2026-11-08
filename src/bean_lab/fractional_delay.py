"""Offline fractional delay by Kaiser-windowed sinc interpolation.

A delay of D samples (D = n + f, n integer, 0 <= f < 1) is applied as an integer
shift n plus a short FIR with taps h[r] = sinc(r - f) * kaiser(r - f), r = -H+1..H.
For integer D the kernel degenerates exactly to a unit impulse. All channels use
the same kernel length, so the group delay of every channel is exactly its
requested delay and inter-channel timing is preserved.

Output policy (explicit, see apply_delay):
  * head: zeros (no wrap-around, nothing before the delayed signal starts);
  * tail: never truncated; length = len(x) + ceil(max_delay) + HALF_TAPS samples.
"""

from __future__ import annotations

import math

import numpy as np

HALF_TAPS = 32  # H: kernel has 2*H taps
KAISER_BETA = 8.6  # ~ -90 dB stopband for the window


def kernel(frac: float, half_taps: int = HALF_TAPS, beta: float = KAISER_BETA) -> np.ndarray:
    """Taps h[r], r = -H+1..H, approximating a delay of `frac` samples (0 <= frac < 1)."""
    r = np.arange(-half_taps + 1, half_taps + 1, dtype=np.float64)
    if frac == 0:  # exact unit impulse (np.sinc leaves ~1e-17 residue at integers)
        return (r == 0).astype(np.float64)
    x = r - frac
    arg = np.clip(1.0 - (x / half_taps) ** 2, 0.0, None)
    window = np.i0(beta * np.sqrt(arg)) / np.i0(beta)
    return np.sinc(x) * window


def output_length(n_frames: int, max_delay_samples: float, half_taps: int = HALF_TAPS) -> int:
    return n_frames + math.ceil(max_delay_samples) + half_taps


def delay_signal(x: np.ndarray, delay_samples: float, out_frames: int,
                 half_taps: int = HALF_TAPS) -> np.ndarray:
    """Delay mono `x` by `delay_samples` (>= 0, may be fractional) into out_frames samples."""
    if x.ndim != 1:
        raise ValueError("expected a mono (1-D) signal")
    if not math.isfinite(delay_samples) or delay_samples < 0:
        raise ValueError("delay must be finite and nonnegative")
    n = int(math.floor(delay_samples))
    h = kernel(delay_samples - n, half_taps)
    full = np.convolve(x, h)  # full[k] = sum_r' h[r'] x[k - r'], r' = r + H - 1
    # y[m] = full[m - n + H - 1]
    out = np.zeros(out_frames, dtype=np.float64)
    offset = n - half_taps + 1
    m_start = max(0, offset)
    m_end = min(out_frames, len(full) + offset)
    if m_end > m_start:
        out[m_start:m_end] = full[m_start - offset:m_end - offset]
    return out


def apply_delays(x: np.ndarray, delays_seconds: list[float] | tuple[float, ...],
                 sample_rate: int) -> np.ndarray:
    """Return an (frames, len(delays)) array; column i is x delayed by delays_seconds[i].

    No gain change is applied. Frames = output_length(len(x), max delay in samples).
    """
    if sample_rate < 1:
        raise ValueError("sample_rate must be positive")
    x = np.asarray(x, dtype=np.float64)
    if x.ndim != 1 or len(x) == 0 or not np.all(np.isfinite(x)):
        raise ValueError("expected a nonempty finite mono signal")
    if len(delays_seconds) == 0:
        raise ValueError("at least one delay is required")
    if any(not math.isfinite(d) or d < 0 for d in delays_seconds):
        raise ValueError("delays must be finite and nonnegative")
    delays = [float(d) * sample_rate for d in delays_seconds]
    frames = output_length(len(x), max(delays))
    return np.stack([delay_signal(x, d, frames) for d in delays], axis=1)
