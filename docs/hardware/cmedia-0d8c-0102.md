# Verified analog channel mapping — C-Media USB 7.1

Hardware: USB `0d8c:0102`, ALSA card `ICUSBAUDIO7D`.

**Method:** direct 8-channel S16_LE / 48 kHz ALSA PCM playback with `bean-lab tone --channel N`, observed with headphones at each physical output. Indices are **1-based PCM channel numbers**, not `speaker-test -s` selectors.

| Physical output | Left contact | Right contact |
| --- | ---: | ---: |
| Plug 1 — Front | PCM 2 | PCM 1 |
| Plug 2 — Side | PCM 8 | PCM 7 |
| Plug 3 — Rear | PCM 5 | PCM 6 |
| Plug 4 — Center/Sub | PCM 4 | PCM 3 |

An additional `Line Out` appears to mirror the Front L/R signal; this is **not** an independently verified ninth/tenth output.

This mapping was observed on one physical interface. Do not assume it is universal across other devices with the same USB ID.

## Naming policy

Keep these distinct:

1. **PCM index**: 1–8, addressing the interleaved audio stream.
2. **Physical output**: plug number and left/right contact.
3. **Array element**: arbitrary loudspeaker position in an experimental geometry.

Array element assignment must be explicit and configurable. Do not silently equate `Front`, `Side`, `LFE`, etc. with beamforming array positions.

## Validation command

```sh
for n in {1..8}; do
  echo "Testing PCM channel $n"
  bean-lab tone --channel "$n" --seconds 1 --gain 0.05
done
```

Use a low headphone level. Analog frequency response, channel-to-channel phase, gain matching, and clock stability remain unmeasured.
