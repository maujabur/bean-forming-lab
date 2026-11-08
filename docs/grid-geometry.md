# Speaker geometry JSON

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

Two input formats; one internal representation: XYZ coordinates in **meters**. All JSON coordinates and pitches use **millimeters** (`"units": "mm"`).

## Integer grid

`configs/arrays/staggered-3-2-3.json` models:

```text
1 . 2 . 3
. . . . .
. 4 . 5 .
. . . . .
6 . 7 . 8
```

Each speaker has `"cell": [column, row]`. Cell indices are **nonnegative integers**. Cell `[0, 0]` is the **top-left origin**, not the acoustic center. X increases right and Y increases down, viewed from the front; Z=0 for the flat grid.

The `grid` object accepts:
- `pitch`: optional default spacing on both axes, in mm.
- `pitch_x`: optional horizontal override.
- `pitch_y`: optional vertical override.

Both effective pitches must be defined, finite, and strictly positive. Thus `{"pitch":30}`, `{"pitch":30,"pitch_y":40}`, and `{"pitch_x":25,"pitch_y":40}` are all valid. A grid has a fixed origin; there is no `origin_mm` field.

The pitch measures **distance between integer grid coordinates**, not enclosure width. With pitch 30 mm, columns 0 and 2 are 60 mm apart.

## Free XYZ positions

Omit `grid` and provide `"position": [x, y, z]` in millimeters for each speaker. Example: `{"speaker":3,"position":[32.5,58,0]}`. Free coordinates may be negative; the coordinate origin is simply a reference, not an acoustic center.

Do not mix `cell` and `position` in one array. Both formats are normalized by `load_array()` to XYZ **meters**, exposed through `position_for(speaker)`.

Speaker IDs are stable and independent of hardware PCM channel mapping. The loader does not enforce enclosure clearance or mechanical mounting constraints.
