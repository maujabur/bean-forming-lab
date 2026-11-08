# Integer grid geometry

For interchangeable MDF/3D-printed speaker mounts, speaker geometry can be expressed as integer grid cells instead of direct XYZ coordinates.

`configs/arrays/staggered-3-2-3.json` models eight speakers in this front-view pattern:

```text
1 . 2 . 3
. . . . .
. 4 . 5 .
. . . . .
6 . 7 . 8
```

`grid.pitch_mm` is the distance between adjacent **integer coordinate positions**, not necessarily the width of a speaker enclosure. For example, `pitch_mm: 30` yields 60 mm center spacing along each occupied row and 30 mm row offset. Choose the pitch only after measuring the enclosure and mounting constraints.

`cell: [column, row]` identifies the speaker center. X increases right, Y increases down in front view, and Z is constant for this flat grid. Optional `grid.origin_mm: [x, y, z]` offsets the array origin. The loader converts millimeters to meters and exposes the same `position_for(speaker)` API as the original `units: meters` XYZ format.

**Hardware routing is unaffected.** Speaker IDs remain 1–8 and are resolved separately through the hardware JSON. The software does not yet enforce enclosure clearance, mounting-hole geometry, or batten placement; these are mechanical design constraints.

The pattern is a prototype layout, not a validated acoustic optimum.
