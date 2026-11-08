# Speaker identity, geometry, and hardware mapping

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

The lab uses **stable speaker numbers** (1–8) throughout experiments and DSP. There are two independent JSON files:

- `configs/arrays/linear-4.json`: physical speaker positions in meters. No PCM channels.
- `configs/hardware/cmedia-0d8c-0102.json`: interface ALSA device and **speaker → PCM** routing. No positions.

The C-Media mapping uses an **example assignment** of speakers 1–8 to output contacts in the following order: Front L, Front R, Side L, Side R, Rear L, Rear R, Center/Sub L, Center/Sub R. The observed raw PCM channel indices for these contacts are `2, 1, 8, 7, 5, 6, 4, 3` respectively. This is **not** a claim that the user's physical speaker labels have already been assigned to those connectors.

Example hardware JSON:

```json
{
  "name": "cmedia-0d8c-0102",
  "device": "hw:CARD=ICUSBAUDIO7D,DEV=0",
  "channels": 8,
  "speaker_to_pcm": {
    "1": 2, "2": 1, "3": 8, "4": 7,
    "5": 5, "6": 6, "7": 4, "8": 3
  }
}
```

The example `linear-4.json` has four physical positions identified by `speaker` 1–4. Positions are independent of PCM channels.

## Test

```sh
bean-lab tone --speaker 1 --gain 0.05
bean-lab tone --speaker 3 --hardware configs/hardware/cmedia-0d8c-0102.json --gain 0.05
bean-lab tone --speaker 2 --array configs/arrays/linear-4.json --gain 0.05
```

For raw diagnostics bypassing speaker mapping:

```sh
bean-lab tone --channel 8 --gain 0.05
```

Changing the sound interface or cabling changes **only the hardware JSON**, provided speaker IDs and physical geometry stay the same. The array geometry is not yet used to compute steering/focusing delays.

The loader rejects duplicate speaker IDs, duplicate PCM assignments, out-of-range channels and malformed geometry. Audio output is always sent as an interleaved PCM frame in **PCM order**, not speaker-number order.
