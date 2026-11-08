# Virtual microphone

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

`bean-lab listen` simulates what an ideal omnidirectional microphone at any XYZ point would receive from the array, and writes it as a **mono** WAV. `bean-lab compare` runs four controlled experiments (A to D) and prints their levels. No audio hardware or hardware JSON is needed.

**Everything here is the behavior predicted by a simple model, not an acoustic measurement.** Audible or numerical differences between the WAVs show what the model predicts, nothing more.

## Physical model

```
p(t) = sum_i  (r_ref / r_i) * s_i(t - e_i - r_i / c)
```

- `s_i`: signal emitted by speaker `i`, defined as the pressure waveform it produces **at the reference distance** `r_ref = 1 m`. This makes `r_ref / r_i` dimensionless: a unit sine heard 1 m from one speaker has unit amplitude, 2 m gives 0.5.
- `r_i`: distance speaker to microphone; `c = 343 m/s` (`--speed-of-sound`).
- `e_i`: extra emission delay, e.g. the focusing delay from `docs/static-focusing.md`; zero in simultaneous mode.

Assumptions: point, omnidirectional sources; free field; no reflections, obstacles, absorption, directivity or transducer response; every speaker has the same amplitude factor `A_i = 1`. Microphone positions closer than 1 nm to a speaker are rejected (a singularity; below float resolution for mm input). There is no other distance limit: a microphone 10 mm away is valid but very loud, and the clipping check will say so.

## Coordinates

As in the focusing document: the CLI takes millimeters in the array frame (grid cell `[0,0]` is the origin, grid arrays at z = 0); internally meters. `--target` is where the array is focused, `--mic` is where it is heard. They are independent.

## Time reference and length

- `t = 0` is sample 0 of the emitted tone, at every speaker. The recording starts at `t = 0`.
- **Propagation delay is not removed.** A source 1 m away is silent for the first 2.9 ms (about 140 samples) of the recording. Emission and propagation delays are preserved and add.
- The tail is not truncated. Frames = `len(source) + ceil(largest total delay in samples) + 32`.
- Emission delay and propagation delay are applied as **one** fractional delay per speaker (windowed-sinc, see the focusing document), never as two successive filters. Measured on the staggered array, 16 kHz sine, focused, off-axis mic: single pass max error 2e-6 versus 3e-6 for two cascaded stages (amplitude 0.05 per speaker); both are negligible, but a single pass was never worse and uses one filter instead of two.

## Safety and levels

- No normalization of any kind: every WAV is written at the absolute scale predicted by the model, so files can be compared by level.
- If the received signal exceeds full scale (|x| > 1) the command stops with an error and **writes nothing**; lower `--gain`. `compare` checks all four cases before writing any file.
- Because of `1/r`, a microphone nearer than 1 m is louder than a speaker's nominal level, and eight coherent speakers add up to about eight times the nominal amplitude.

## Measurements printed

For each speaker: distance (mm), travel time (ms), emission delay (us), arrival time (emission + travel, ms), `1/r` gain.
For the signal: RMS, peak, level in dBFS (relative to amplitude 1.0), and RMS in dB relative to **one speaker's sine heard at 1 m** (`gain / sqrt(2)`; this reference is valid for the tone sources used by the CLI; eight speakers in perfect coherence at 1 m give +18.06 dB).

### Analysis window

RMS is never taken over the silence before arrival or over the fades. The window is the interval in which **every** speaker of **every** case is fully on, after the 5 ms fade and the 32-sample interpolation margin:

```
start = max(all arrival times) + 5 ms + 32/fs
end   = min(all arrival times) + duration - 5 ms - 32/fs
```

trimmed to a whole number of periods of `--frequency`, so a partial cycle does not bias the RMS. `compare` uses the intersection over A to D, so all four levels are measured over the same samples. `listen` uses its own window unless you pass `--window T0 T1` (seconds); `compare` prints the window to reuse.

## The four experiments

All with the same array, frequency, duration and gain. Here the target is (60, 60, 1000) mm and the off-target microphone is (260, 60, 1000) mm, 200 mm to the side.

| Exp | Emission | Microphone |
|---|---|---|
| A | simultaneous | target |
| B | focused | target |
| C | simultaneous | off target |
| D | focused | off target |

```bash
.venv/bin/bean-lab compare --array configs/arrays/staggered-3-2-3.json \
  --target 60 60 1000 --off-target 260 60 1000 \
  --frequency 1000 --duration 2 --gain 0.05 --outdir /tmp/exp-1000
```

The same four cases one by one (window from the `compare` output):

```bash
for m in "simultaneous 60 60 1000 A" "focused 60 60 1000 B" \
         "simultaneous 260 60 1000 C" "focused 260 60 1000 D"; do
  set -- $m
  .venv/bin/bean-lab listen --array configs/arrays/staggered-3-2-3.json \
    --mode $1 --target 60 60 1000 --mic $2 $3 $4 \
    --frequency 1000 --duration 2 --gain 0.05 \
    --window 0.008687500 1.996687500 --wav /tmp/bean-mic-$5.wav
done
```

### Results (RMS over the common window)

1 kHz, window 0.0086875 to 1.9966875 s:

| Exp | RMS | dBFS | Peak |
|---|---|---|---|
| A simultaneous, target | 0.282096 | -10.99 | 0.398692 |
| B focused, target | 0.282179 | -10.99 | 0.398386 |
| C simultaneous, off | 0.273107 | -11.27 | 0.386231 |
| D focused, off | 0.273175 | -11.27 | 0.386204 |

Differences: B - A = +0.003 dB, D - C = +0.002 dB, A - C = +0.281 dB, B - D = +0.282 dB.

8 kHz (`--frequency 8000`, `--outdir /tmp/exp-8000`), window 0.0086875 to 1.9971875 s:

| Exp | RMS | dBFS | Peak |
|---|---|---|---|
| A | 0.276893 | -11.15 | 0.380475 |
| B | 0.282181 | -10.99 | 0.398392 |
| C | 0.095359 | -20.41 | 0.128158 |
| D | 0.093980 | -20.54 | 0.126223 |

Differences: B - A = +0.164 dB, D - C = -0.126 dB, A - C = +9.259 dB, B - D = +9.550 dB.

## Interpreting the numbers

- The staggered array is only 120 mm wide and the target is 1 m away, so the emission delays are under 10 us. At 1 kHz that is about 3 degrees of phase (9 us of a 1000 us period); at 8 kHz about 26 degrees. Focusing changes the level at the target by only 0.003 dB at 1 kHz and about 0.16 dB at 8 kHz. Simultaneous emission is already almost in phase on the array axis, so the focusing gain is small. Focusing cannot exceed the coherent sum (about 18.06 dB re one speaker at 1 m, `RMS = 8 * 0.05 / sqrt(2) = 0.2828`); B reaches 0.2822.
- The off-target microphone is quieter at 8 kHz mainly because of the array's own interference pattern, which is the same with and without focusing (A - C and B - D are similar). D - C is the effect of the delays alone at that position: it can be positive or negative. Focusing does not guarantee less sound elsewhere, nor more at every instant or frequency.
- A level difference of a fraction of a dB is below what the model's own simplifications (no reflections, ideal sources) justify interpreting physically.
- Also compare peaks: a coherent sum has a peak close to `sqrt(2)` times its RMS for tones; large departures suggest interference with a different phase relationship.

## Comparing WAVs

1. Generate with the same `--array`, `--frequency`, `--duration`, `--gain`, `--rate`; vary only `--mode`/`--mic`.
2. Use `compare`, or `listen --window` with identical `T0 T1` for every file.
3. Do not normalize in an editor; relative levels are the result. Play at low volume (peaks of 0.4 of full scale are already loud).
4. Because the tone is a single frequency, a mono WAV cannot reveal spatial behavior by listening alone; use the printed levels.

## Limitations

Free field only; no directivity, no frequency or phase response of real speakers or amplifier channels (the eight outputs have unmeasured gain and timing); no air absorption or temperature dependence (c is a parameter); sources are point-like although the real speakers are not; the microphone is a point, omnidirectional and ideal. Offline only, static focus, fixed `c`. Not binaural: this is a single point pressure signal.

## Conclusion

The virtual microphone gives a consistent free-field prediction (`1/r` with a 1 m reference, summed propagation and emission delays in one fractional delay), and its levels can be compared without normalizing the files. In the four experiments the effect of focusing on the 120 mm staggered array is small: +0.003 dB at the target at 1 kHz and +0.164 dB at 8 kHz, because the emission delays are under 10 microseconds. The level differences between target and off-target positions (0.281 dB at 1 kHz, 9.259 dB at 8 kHz) come mostly from the array's own interference pattern, not from the delays. These numbers are predictions of an ideal model, not an acoustic measurement; the next step is to measure a few points (the focus and one null) with a real microphone.
