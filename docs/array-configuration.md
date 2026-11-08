# JSON speaker array configuration

Array configurations live in `configs/arrays/`. Positions use meters and Cartesian `[x, y, z]` coordinates.

```json
{
  "name": "linear-4",
  "units": "meters",
  "elements": [
    {"id": "SPK1", "pcm_channel": 2, "position": [-0.09, 0, 0]},
    {"id": "SPK2", "pcm_channel": 1, "position": [-0.03, 0, 0]},
    {"id": "SPK3", "pcm_channel": 8, "position": [0.03, 0, 0]},
    {"id": "SPK4", "pcm_channel": 7, "position": [0.09, 0, 0]}
  ]
}
```

`pcm_channel` is **1-based raw PCM channel indexing**, measured using `bean-lab tone --channel N`, not `speaker-test -s`. The sample configuration uses the observed Front and Side output mapping of one C-Media interface. It is not a universal mapping.

## Test by array element ID

```sh
bean-lab tone --array configs/arrays/linear-4.json --element SPK1 --gain 0.05
bean-lab pulse --array configs/arrays/linear-4.json --element SPK4 --output spk4.wav
```

The loader validates unique element IDs, unique PCM assignments, channel range, and finite 3D positions. Other channels remain silent. Position data is **not yet used to calculate delays**; that is a separate upcoming feature. No calibration parameters are assumed.

A physical output is not the same as an array element. Reassigning speakers should require editing the JSON, not the DSP source code.
