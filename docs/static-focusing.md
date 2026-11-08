# Static acoustic focusing

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

`bean-lab focus` computes per-speaker delays so that wavefronts from all speakers arrive at one point in space at the same instant, and renders a multichannel WAV (or plays it). The focal point is fixed for the whole signal. Moving focus is not implemented.

## Pipeline (each stage is a separate function)

1. `focus.focus_delays(array, focus, c)`: geometry only. Knows speaker IDs, not PCM channels.
2. `fractional_delay.apply_delays(mono, delays, rate)`: applies the delays to a mono signal.
3. `focus.route_to_pcm(signals, hardware)`: speaker ID -> PCM channel via `hardware_config` (unchanged mapping).
4. `signal.write_wav` / `aplay`: export or playback (same paths as `tone` and `pulse`).

## Coordinates

Array JSON files are in millimeters and are normalized to XYZ meters by `array_config`. The `--target X Y Z` option is in **millimeters, in the same frame as the array** (grid cell `[0,0]` is the origin, grid arrays lie at z = 0). It is an arbitrary point in space, not a grid cell. Internally everything is meters.

## Delay computation

```
d_i     = || focus - p_i ||            (m)
t_i     = d_i / c                      (s), c = 343 m/s by default (--speed-of-sound)
delay_i = (d_max - d_i) / c    >= 0
```

The farthest speaker has delay 0; nearer speakers emit later, so that `delay_i + d_i / c` is the same for every speaker. Results are sorted by speaker number and do not depend on element order in the JSON.

## Fractional delay technique

Kaiser-windowed sinc interpolation (`fractional_delay.py`): a delay of `D = n + f` samples is an integer shift `n` plus a 64-tap FIR (`HALF_TAPS = 32`, Kaiser beta 8.6) whose taps are `sinc(r - f) * kaiser(r - f)`. For `f = 0` the kernel is an exact unit impulse. Every channel uses the same kernel length, so relative timing between channels equals the requested relative delay.

At 48 kHz one sample is 20.83 us (7.15 mm of air).

Measured error, steady-state, unit-amplitude sine, fractional parts 0.01 to 0.99 (the `test_fractional_delay_*` tests assert looser bounds: 1e-4 waveform, 1e-3 sample phase):

| Frequency | Max waveform error |
|-----------|--------------------|
| 1 to 16 kHz | < 1.3e-5 (about -98 dB) |
| 19 kHz | 3.8e-5 |
| 21 kHz | 6.2e-5 |
| 22 kHz | 3.5e-4 |

Limitations:

- Accuracy degrades near Nyquist (24 kHz); keep content below about 20 kHz.
- Offline only (whole signal in memory). Not a streaming processor, and delays cannot change during playback.
- The kernel is symmetric (non-causal): there is a small amount of pre-ringing of up to 32 samples (0.67 ms) before a sharp onset. The `focus` command applies a 5 ms raised-cosine fade to the tone to avoid clicks and to keep these transients negligible.
- The FIR passband is not perfectly flat near Nyquist, so amplitudes of different channels can differ slightly when fractional parts differ (see table).

## Output length, padding and level policy

- Head: silence (zeros) before each delayed onset; no wrap-around.
- Tail: not truncated. Frames = `len(input) + ceil(max_delay_samples) + 32`. The extra frames contain the delayed tail and the filter ring-down.
- All channels have equal length; unused PCM channels are silent.
- No normalization and no per-channel gain change: relative levels are exactly the input level on every active channel. If any sample exceeds full scale, `route_to_pcm` raises an error instead of clipping. Choose `--gain` accordingly (`0.05` by default for headphone safety; `0 < gain <= 1`).

## Generate a focused WAV (no hardware needed)

```bash
.venv/bin/bean-lab focus \
  --array configs/arrays/staggered-3-2-3.json \
  --hardware configs/hardware/cmedia-0d8c-0102.json \
  --target 60 60 1000 \
  --frequency 1000 \
  --duration 2 \
  --gain 0.05 \
  --wav /tmp/bean-focus.wav
```

Without `--wav`, the signal is played through `aplay` on the hardware device (same mechanism as `tone`; `--device` overrides it). Paths for `--hardware` default to `configs/hardware/cmedia-0d8c-0102.json`, relative to the current directory. Unlike `tone`/`pulse`, `focus` is a subcommand with its own options (`--target`, `--wav` instead of `--output`, `--duration` instead of `--seconds`).

## Reading the output

```
speaker  pcm  distance_mm  delay_us  delay_samples
      1    2     1003.594      0.00          0.000
      4    7     1000.450      9.17          0.440
```

- `distance_mm`: speaker to focal point.
- `delay_us` / `delay_samples`: emission delay. Zero means the speaker is (one of) the farthest and starts first.
- `pcm`: 1-based PCM channel carrying that speaker.

The staggered 3-2-3 array has a 30 mm pitch, so for a target 1 m away the delay differences are only a few microseconds (a fraction of one sample). That is why fractional delays are required.

## Physical limitations

**Matching arrival times does not guarantee spatial isolation of the audio.** Sound is still emitted by every speaker in every direction. Delays only make the signals add in phase at the focal point; elsewhere the result depends on frequency, geometry, spacing, speaker directivity and room acoustics. With a 30 mm-pitch array and a 1 kHz tone (wavelength 343 mm), the aperture is small compared with the wavelength and any focusing is weak; the focal "spot" is much larger than the array. Do not expect a localized spot; the actual field has not been measured.

Not modeled or compensated: speaker frequency/phase response, per-channel gain mismatches, DAC/USB timing offsets between channels, air absorption, 1/r amplitude differences, reflections, temperature (c varies about 0.6 m/s per degree C). Channel timing of the interface is unmeasured (see README).
