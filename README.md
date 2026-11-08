# Bean Forming Lab

Open experimental laboratory for multichannel audio, loudspeaker arrays, beam steering and near-field focusing.

## Current hardware

- C-Media CM106-like USB device: `0d8c:0102`, identified by ALSA as `ICUSBAUDIO7D`.
- USB PCM modes reported: **8/6/4/2 channels, S16_LE at 44.1/48 kHz** (plus stereo 96 kHz).
- Four stereo analog connectors: Front L/R, Side L/R, Rear L/R, Center/LFE. A fifth 'line out' appears to mirror Front L/R.
- All eight outputs were heard through headphones, including high-frequency content from LFE. **Frequency flatness and phase coherence have not yet been measured.**

## Installation (Debian / Linux)

```sh
sudo apt install python3-venv alsa-utils
python3 -m venv .venv
. .venv/bin/activate
pip install -e '.[test]'
pytest
```

## First tests (headphones, low volume)

The CLI addresses **raw ALSA PCM channel indices**, unlike `speaker-test -s`, whose channel selection order is different. **Do not assume the physical jack mapping until verified.**

```sh
bean-lab tone --channel 1 --gain 0.05

for n in 1 2 3 4 5 6 7 8; do
    echo "PCM channel $n"
    bean-lab tone --channel "$n" --seconds 1 --gain 0.05
done

bean-lab tone --channel 8 --frequency 3000 --output test-ch8.wav
bean-lab pulse --channel 4 --output impulse-ch4.wav
```

Default device: `hw:CARD=ICUSBAUDIO7D,DEV=0` at **48 kHz, S16_LE, 8 channels**, without software mixing or format conversion. The direct hardware path may be unavailable while PipeWire or another application holds the device.

**Safety:** keep gain low during headphone testing, especially with tones and impulses. Do not connect loudspeaker power-amplifier outputs to sound-card input sockets. Analog loopback measurements require verified safe input levels.

## Known caveats

1. USB descriptors do **not** establish analog gain matching, phase equality, or clock behavior.
2. `speaker-test -s` uses logical speaker selection, **not** interleaved PCM indexing.
3. LFE was audible at high frequencies, but no calibrated response has yet been measured.
4. `speaker-test -f 14000` printed `Sine wave rate is 8000.0000Hz` in the observed test; our generator checks frequency against Nyquist.

## Roadmap

- Verify physical PCM index-to-jack mapping.
- Measure output amplitudes, polarity and channel timing.
- Add per-channel gain and calibrated fractional delays for broadband signals.
- Model steering/focusing geometries and compare with microphone measurements.
## Static focusing
`bean-lab focus` delays one tone per speaker so wavefronts coincide at a chosen point, and writes a multichannel WAV with no hardware attached. See [docs/static-focusing.md](docs/static-focusing.md).
## Virtual microphone
`bean-lab listen` and `bean-lab compare` simulate what an ideal omnidirectional microphone would receive from the array (free field, 1/r), as mono WAV, with no hardware. See [docs/virtual-microphone.md](docs/virtual-microphone.md).
## Field maps
`bean-lab field-map` renders 2-D steady-state pressure maps (simultaneous, focused, difference in dB) of a tone from the array, free-field model. See [docs/field-map.md](docs/field-map.md). PNG output needs matplotlib (`pip install -e '.[plot]'`).
## Parametric study
`bean-lab campaign` runs the frequency, spacing and focal-direction series and writes comparison plates and indicators. Results and discussion: [docs/parametric-study.md](docs/parametric-study.md).
## Broadband study
`bean-lab broadband-map` and `bean-lab broadband-campaign` compute time-averaged energy maps for band-limited signals (real focusing delays, mean-square combination of lines). Results and discussion: [docs/broadband-study.md](docs/broadband-study.md).
## Three-dimensional field
`bean-lab field-cuts` and `bean-lab field3d-campaign` compute XY/XZ/YZ cuts (and an optional coarse volume) of the broadband energy field through the focus, comparing the line array with planar layouts. Results and discussion: [docs/three-dimensional-field.md](docs/three-dimensional-field.md).
