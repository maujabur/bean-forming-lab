# Three-dimensional field of the arrays

*[maujabur/bean-forming-lab](https://github.com/maujabur/bean-forming-lab) · made with Kiro*

Orthogonal cuts (XY, XZ, YZ) of the broadband energy field through the focus, for the 8-speaker line array (L), the staggered 3-2-3 array (P) and a wide planar reference (W). Free-field ideal point-source model: **predictions in dB relative to one speaker at 1 m, not measurements and not dB SPL.** The energy model, the real focusing delays and the spectral quadrature are those of `docs/broadband-study.md` (a tone is a one-line spectrum).

## Geometries

| Tag | File | Elements | Extent (mm) | Aperture center (mm) | Nearest neighbors (mm) | Plane of the array | Beam axis used |
|---|---|---|---|---|---|---|---|
| L | `configs/arrays/line-8-60mm.json` | 8 on a line along x | 420 x 0 x 0 | (210, 0, 0) | 60 | line (no plane) | +y (by convention since Mission 01; any direction perpendicular to x is equivalent) |
| P | `configs/arrays/staggered-3-2-3.json` | 3-2-3 staggered | 120 x 120 x 0 | (60, 60, 0) | 60 | **z = 0 (XY plane)**, normal +z | +z |
| W | `configs/arrays/staggered-3-2-3-wide.json` (new) | the same file with grid pitch 105 mm instead of 30 mm | 420 x 420 x 0 | (210, 210, 0) | 210 | z = 0, normal +z | +z |

Findings from reading the files (nothing in the existing files was modified):

- The 3-2-3 file is on a 30 mm grid (cells 0, 2, 4 in rows 0, 2, 4) and lies in the **XY plane at z = 0**. Its normal is therefore +z, not +y: the focus "800 mm in front of the center along the normal" is **(60, 60, 800) mm**. It has a 120 x 120 mm aperture, 3.5 times smaller than the linear array's 420 mm, and a mirror symmetry in x and in y (not a 90-degree one).
- The linear array has **axial symmetry** about its own axis (the x axis): its field depends only on x and on the distance to that axis. This single fact decides most of the results below.
- **W** is a controlled comparison, not a design proposal: the 3-2-3 topology scaled by exactly 3.5 so that its largest dimension (420 mm) equals the linear array's. Criteria: same element count, same topology as P (isolates the aperture: P to W), same largest dimension as L (isolates the dimensionality: L to W). The price of a fixed element count: the nearest-neighbor spacing becomes 210 mm, so W is grating-free only below c/d = 1.6 kHz. Differences between L and W therefore mix dimensionality with the 3.5-fold larger spacing, and differences between P and W mix aperture with spacing; the report states which where it matters.

## Scenarios and reproduction

Planes always pass through the focus. In each cut the first named axis is horizontal. Windows: the beam axis from 100 to 2000 mm, the other axes focus +/- 1000 mm; 5 mm grid (401 x 381 or 401 points per axis). Optional 3-D volume on the same box with a 25 mm grid. c = 343 m/s, 0.05 per speaker at 1 m, line spacing 25 Hz.

| Series | Array | Focus (mm) | Spectra |
|---|---|---|---|
| A | L | (350, 800, 0) (the focus used since Mission 03) | tone 4 kHz, 2-4 kHz, 1-8 kHz |
| A0 | L | (210, 800, 0), on the axis of symmetry of the aperture | same |
| B | P | (60, 60, 800) | same |
| W | W | (210, 210, 800) | same |
| C | L, P, W | 20 degrees toward +x at 800 mm from the aperture center: L (483.6, 751.8, 0), P (333.6, 60, 751.8), W (483.6, 210, 751.8) | tone 4 kHz, 2-4 kHz |
| D | L | (800, 800, 0) (36.4 degrees, 994 mm) | 1.5-3.5 kHz, 2-4 kHz, 1-8 kHz |

Series A0 (not in the mission text) gives the linear array a focus on its own axis, so that L, P and W are compared with the same 0-degree steering.

```bash
.venv/bin/bean-lab field3d-campaign --volume --outdir results/field3d     # everything, about 15 to 20 minutes
.venv/bin/bean-lab field3d-campaign --series B W --outdir results/field3d  # some series (add --volume for volumes)
.venv/bin/bean-lab field-cuts --array configs/arrays/staggered-3-2-3.json --target 60 60 800 \
    --fmin 2000 --fmax 4000 --outdir /tmp/cuts-P                           # one scenario (cuts.png, data.npz)
```

`results/` is git-ignored (306 MB: energy and dB cuts, per-scenario `cuts.png`). Committed: the comparison plates and a few single-scenario plates (`docs/images/field3d/`), indicators (`docs/data/field3d-indicators.csv`, `.json`). Display options (`--vmax`, `--range-db`, `--diff-limit`) change no data (tested). Common color scale -20..20 dB for energy, +/-10 dB for differences, arrows on the colorbar where the display saturates.

## Indicators (definitions)

Levels are dB re one speaker at 1 m (energy). To keep the near field (which dominates any raw level through 1/r) out of the statistics, every region indicator is evaluated in a **region of interest (ROI)**: the shell 0.75 R_f <= r <= 1.25 R_f around the aperture center (R_f = center-to-focus distance), within 80 degrees of the beam axis, on the **spreading-compensated** level L + 20 log10(r / R_f).

- *Focal region at -X dB* (X = 3, 6): the connected region of the ROI containing the focus where the compensated level is at least L_focus - X.
  - In a cut it is a **cut area** (cm2) with its extents; in a volume it is a **volume** (liters). A `+` marks a region that touches the cut window, the ROI shell, or the volume box: the number is then a lower bound.
  - The ROI shell is 0.5 R_f thick (400 mm at 800 mm): any axial extent above about 400 mm is shell-limited by definition. It says "at least the shell", not a measured depth.
- *Width along an axis at -X dB*: length of the run of ROI cells above L_focus - X on the axis line through the focus (the same line appears in two cuts and agrees, tested); `+` when the run ended at the shell or window rather than at a drop below the threshold.
- *Interior secondary maximum*: the highest interior local maximum of the compensated level outside the (dilated) -3 dB region; contrast = L_focus - that maximum. A lobe whose compensated level keeps changing along the range has no such maximum (`n/a`).
- *In-plane arc*: for planes that contain the aperture center and the beam axis, the Mission 04/05 constant-range arc (radius R_f, within 80 degrees) inside that plane: contrast between the lobe at the focus direction and the highest other lobe, and the share of the arc energy within +/-10 degrees of the focus.
- *Energy share*: compensated energy inside the focal region as a fraction of the ROI total; given for the focused field and for simultaneous emission in the same region.
- *Cells above focus -X dB*: share of the ROI cells (or shell voxels) at or above L_focus - X (X = 3, 6, 10), focused and simultaneous.

## Plates

Linear versus planar, same steering, one plate per spectrum (columns L, P, W; rows: the plane that holds the beam axis and x, the plane that holds the beam axis and the other axis, and the transverse plane):

![Broadside, 4 kHz tone](images/field3d/compare-broadside-tone.png)

![Broadside, 2-4 kHz](images/field3d/compare-broadside-b24.png)

![Broadside, 1-8 kHz](images/field3d/compare-broadside-b18.png)

Lateral focus (20 degrees):

![Lateral, 4 kHz tone](images/field3d/compare-lateral-tone.png)

![Lateral, 2-4 kHz](images/field3d/compare-lateral-b24.png)

Extreme lateral focus of the linear array, three bands (1.5-3.5, 2-4, 1-8 kHz):

![Extreme focus, bands](images/field3d/compare-extreme-bands.png)

Complete plates (simultaneous, focused, difference for the three cuts) of selected scenarios:

![A, tone: the ring](images/field3d/cuts-A-tone.png)

![A, 2-4 kHz](images/field3d/cuts-A-b24.png)

![B, 2-4 kHz](images/field3d/cuts-B-b24.png)

![W, 2-4 kHz](images/field3d/cuts-W-b24.png)

![D, 1.5-3.5 kHz](images/field3d/cuts-D-b1535.png)

## Results

Widths are at -3 dB and -6 dB along x / y / z, in mm; `+` = limited by the shell or window (lower bound).

### Series A and A0 (linear array)

| Scenario | Spectrum | Focus (mm) | Range (mm) / angle from beam axis (deg) | Focus level (dB) | Gain (dB) | Width x / y / z at -3 dB (mm) | at -6 dB (mm) |
|---|---|---|---|---|---|---|---|
| A0/tone | 4 kHz tone | (210.0, 800.0, 0.0) | 800 / 0.0 | 19.88 | +2.52 | 135 / 405+ / 1205+ | 175 / 405+ / 1205+ |
| A0/b24 | 2 kHz to 4 kHz | (210.0, 800.0, 0.0) | 800 / 0.0 | 19.88 | +1.41 | 175 / 405+ / 1205+ | 245 / 405+ / 1205+ |
| A0/b18 | 1 kHz to 8 kHz | (210.0, 800.0, 0.0) | 800 / 0.0 | 19.88 | +3.07 | 115 / 405+ / 1205+ | 205 / 405+ / 1205+ |
| A/tone | 4 kHz tone | (350.0, 800.0, 0.0) | 812 / 9.9 | 19.76 | +7.26 | 135 / 415+ / 1215+ | 190 / 415+ / 1215+ |
| A/b24 | 2 kHz to 4 kHz | (350.0, 800.0, 0.0) | 812 / 9.9 | 19.76 | +6.31 | 180 / 415+ / 1215+ | 255 / 415+ / 1215+ |
| A/b18 | 1 kHz to 8 kHz | (350.0, 800.0, 0.0) | 812 / 9.9 | 19.76 | +5.81 | 125 / 415+ / 1215+ | 210 / 415+ / 1215+ |

| Scenario | Cut | Area -3 dB (cm2) | Area -6 dB (cm2) | Energy share in -6 dB region, focused | simultaneous | ROI cells above focus -6 dB (%) | Interior secondary contrast (dB) | In-plane arc contrast (dB) | In-plane arc share +/-10 deg |
|---|---|---|---|---|---|---|---|---|---|
| A0/tone | XY | 522+ | 714+ | 0.723 | 0.470 | 8.0 | 12.35 | 12.36 | 0.785 |
| A0/tone | XZ | 1705+ | 2312+ | 0.852 | 0.580 | 20.4 | 12.32 | - | - |
| A0/tone | YZ | 8935+ | 8935+ | 1.000 | 1.000 | 100.0 | n/a | n/a | 0.125 |
| A0/b24 | XY | 701+ | 984+ | 0.774 | 0.606 | 11.0 | n/a | 14.87 | 0.810 |
| A0/b24 | XZ | 2271+ | 3167+ | 0.850 | 0.704 | 28.0 | n/a | - | - |
| A0/b24 | YZ | 8935+ | 8935+ | 1.000 | 1.000 | 100.0 | n/a | n/a | 0.125 |
| A0/b18 | XY | 472+ | 814+ | 0.420 | 0.298 | 9.1 | n/a | 9.17 | 0.475 |
| A0/b18 | XZ | 1530+ | 2625+ | 0.733 | 0.562 | 23.2 | n/a | - | - |
| A0/b18 | YZ | 8935+ | 8935+ | 1.000 | 1.000 | 100.0 | n/a | n/a | 0.125 |
| A/tone | XY | 544+ | 746+ | 0.713 | 0.243 | 8.4 | 12.32 | 12.33 | 0.748 |
| A/tone | XZ | 1777+ | 2415+ | 0.854 | 0.263 | 19.7 | 12.24 | - | - |
| A/tone | YZ | 9192+ | 9192+ | 1.000 | 1.000 | 100.0 | n/a | - | - |
| A/b24 | XY | 732+ | 1027+ | 0.775 | 0.359 | 11.5 | 19.54 | 14.83 | 0.803 |
| A/b24 | XZ | 2366+ | 3299+ | 0.853 | 0.394 | 26.9 | 14.91 | - | - |
| A/b24 | YZ | 9192+ | 9192+ | 1.000 | 1.000 | 100.0 | n/a | - | - |
| A/b18 | XY | 492+ | 851+ | 0.447 | 0.194 | 9.6 | n/a | 9.05 | 0.489 |
| A/b18 | XZ | 1598+ | 2734+ | 0.718 | 0.329 | 22.3 | n/a | - | - |
| A/b18 | YZ | 9192+ | 9192+ | 1.000 | 1.000 | 100.0 | n/a | - | - |

| Scenario | Focal volume -3 dB (liters) | -6 dB (liters) | ROI shell volume in the box (liters) | Energy share in -6 dB volume, focused | simultaneous | Shell above focus -6 dB, focused (%) | simultaneous (%) | Bounding extent of the -6 dB volume x / y / z (mm) |
|---|---|---|---|---|---|---|---|---|
| A0/tone | 121+ | 164+ | 1357 | 0.805 | 0.520 | 12.1 | 12.3 | 225 / 900 / 1975 |
| A0/b24 | 160+ | 222+ | 1357 | 0.813 | 0.646 | 16.3 | 17.8 | 325 / 900 / 1975 |
| A0/b18 | 105+ | 185+ | 1357 | 0.564 | 0.411 | 13.6 | 19.4 | 275 / 900 / 1975 |
| A/tone | 125+ | 169+ | 1406 | 0.797 | 0.255 | 12.0 | 12.2 | 250 / 900 / 1975 |
| A/b24 | 166+ | 229+ | 1406 | 0.810 | 0.361 | 16.3 | 17.9 | 325 / 900 / 1975 |
| A/b18 | 111+ | 191+ | 1406 | 0.565 | 0.250 | 13.6 | 19.2 | 275 / 900 / 1975 |

### Series B and W (planar arrays)

| Scenario | Spectrum | Focus (mm) | Range (mm) / angle from beam axis (deg) | Focus level (dB) | Gain (dB) | Width x / y / z at -3 dB (mm) | at -6 dB (mm) |
|---|---|---|---|---|---|---|---|
| B/tone | 4 kHz tone | (60.0, 60.0, 800.0) | 800 / 0.0 | 19.97 | +0.06 | 405 / 345 / 405+ | 555 / 465 / 405+ |
| B/b24 | 2 kHz to 4 kHz | (60.0, 60.0, 800.0) | 800 / 0.0 | 19.97 | +0.04 | 555 / 465 / 405+ | 815 / 665 / 405+ |
| B/b18 | 1 kHz to 8 kHz | (60.0, 60.0, 800.0) | 800 / 0.0 | 19.97 | +0.10 | 355 / 305 / 405+ | 645 / 655 / 405+ |
| W/tone | 4 kHz tone | (210.0, 210.0, 800.0) | 800 / 0.0 | 19.63 | +12.85 | 115 / 105 / 405+ | 165 / 135 / 405+ |
| W/b24 | 2 kHz to 4 kHz | (210.0, 210.0, 800.0) | 800 / 0.0 | 19.63 | +5.68 | 155 / 135 / 405+ | 225 / 185 / 405+ |
| W/b18 | 1 kHz to 8 kHz | (210.0, 210.0, 800.0) | 800 / 0.0 | 19.63 | +5.71 | 105 / 95 / 405+ | 185 / 175 / 405+ |

| Scenario | Cut | Area -3 dB (cm2) | Area -6 dB (cm2) | Energy share in -6 dB region, focused | simultaneous | ROI cells above focus -6 dB (%) | Interior secondary contrast (dB) | In-plane arc contrast (dB) | In-plane arc share +/-10 deg |
|---|---|---|---|---|---|---|---|---|---|
| B/tone | XY | 1114 | 2208 | 0.722 | 0.716 | 19.5 | n/a | - | - |
| B/tone | XZ | 1573+ | 2154+ | 0.821 | 0.815 | 24.1 | n/a | 11.37 | 0.551 |
| B/tone | YZ | 1348+ | 1824+ | 0.610 | 0.608 | 23.7 | n/a | 5.97 | 0.463 |
| B/b24 | XY | 2090 | 4633 | 0.729 | 0.727 | 41.0 | n/a | - | - |
| B/b24 | XZ | 2122+ | 3012+ | 0.843 | 0.839 | 33.7 | n/a | 13.95 | 0.432 |
| B/b24 | YZ | 1810+ | 2528+ | 0.675 | 0.672 | 28.3 | n/a | 7.98 | 0.400 |
| B/b18 | XY | 904 | 2935 | 0.495 | 0.492 | 26.0 | n/a | - | - |
| B/b18 | XZ | 1411+ | 2447+ | 0.637 | 0.626 | 27.4 | n/a | 8.97 | 0.411 |
| B/b18 | YZ | 1201+ | 2473+ | 0.435 | 0.434 | 72.4 | n/a | 4.57 | 0.291 |
| W/tone | XY | 96 | 185 | 0.067 | 0.004 | 16.6 | 0.12 | - | - |
| W/tone | XZ | 459+ | 637+ | 0.193 | 0.021 | 32.6 | 0.70 | 1.03 | 0.215 |
| W/tone | YZ | 390+ | 537+ | 0.127 | 0.017 | 47.0 | 0.01 | 0.01 | 0.144 |
| W/b24 | XY | 174 | 346 | 0.136 | 0.043 | 17.5 | 3.19 | - | - |
| W/b24 | XZ | 625+ | 880+ | 0.307 | 0.103 | 37.1 | 5.01 | 5.10 | 0.333 |
| W/b24 | YZ | 532+ | 741+ | 0.157 | 0.071 | 78.9 | 3.07 | 3.07 | 0.173 |
| W/b18 | XY | 79 | 234 | 0.076 | 0.030 | 10.2 | 4.46 | - | - |
| W/b18 | XZ | 415+ | 747+ | 0.219 | 0.087 | 8.4 | 6.37 | 6.23 | 0.260 |
| W/b18 | YZ | 348+ | 807+ | 0.135 | 0.095 | 92.4 | 3.88 | 4.08 | 0.154 |

| Scenario | Focal volume -3 dB (liters) | -6 dB (liters) | ROI shell volume in the box (liters) | Energy share in -6 dB volume, focused | simultaneous | Shell above focus -6 dB, focused (%) | simultaneous (%) | Bounding extent of the -6 dB volume x / y / z (mm) |
|---|---|---|---|---|---|---|---|---|
| B/tone | 44+ | 83+ | 1357 | 0.465 | 0.465 | 6.2 | 6.2 | 625 / 575 / 450 |
| B/b24 | 79+ | 161+ | 1357 | 0.558 | 0.557 | 11.9 | 11.9 | 875 / 775 / 475 |
| B/b18 | 36+ | 108+ | 1357 | 0.232 | 0.230 | 17.5 | 17.2 | 725 / 775 / 450 |
| W/tone | 4+ | 7+ | 1357 | 0.021 | 0.002 | 17.0 | 16.8 | 175 / 175 / 425 |
| W/b24 | 7+ | 14+ | 1357 | 0.039 | 0.013 | 19.3 | 3.5 | 275 / 225 / 425 |
| W/b18 | 3+ | 10+ | 1357 | 0.023 | 0.010 | 8.2 | 0.5 | 225 / 275 / 425 |

### Series C (lateral focus, 20 degrees)

| Scenario | Spectrum | Focus (mm) | Range (mm) / angle from beam axis (deg) | Focus level (dB) | Gain (dB) | Width x / y / z at -3 dB (mm) | at -6 dB (mm) |
|---|---|---|---|---|---|---|---|
| C/L-tone | 4 kHz tone | (483.6, 751.8, 0.0) | 800 / 20.0 | 19.92 | +14.68 | 145 / 380+ / 1195+ | 200 / 430+ / 1195+ |
| C/L-b24 | 2 kHz to 4 kHz | (483.6, 751.8, 0.0) | 800 / 20.0 | 19.92 | +10.55 | 195 / 430+ / 1195+ | 275 / 430+ / 1195+ |
| C/P-tone | 4 kHz tone | (333.6, 60.0, 751.8) | 800 / 20.0 | 19.97 | +6.53 | 470 / 345 / 430+ | 650+ / 465 / 430+ |
| C/P-b24 | 2 kHz to 4 kHz | (333.6, 60.0, 751.8) | 800 / 20.0 | 19.97 | +3.34 | 645+ / 465 / 430+ | 740+ / 665 / 430+ |
| C/W-tone | 4 kHz tone | (483.6, 210.0, 751.8) | 800 / 20.0 | 19.68 | +8.46 | 130 / 95 / 370+ | 180 / 135 / 425+ |
| C/W-b24 | 2 kHz to 4 kHz | (483.6, 210.0, 751.8) | 800 / 20.0 | 19.68 | +9.63 | 175 / 135 / 415+ | 245 / 185 / 430+ |

| Scenario | Cut | Area -3 dB (cm2) | Area -6 dB (cm2) | Energy share in -6 dB region, focused | simultaneous | ROI cells above focus -6 dB (%) | Interior secondary contrast (dB) | In-plane arc contrast (dB) | In-plane arc share +/-10 deg |
|---|---|---|---|---|---|---|---|---|---|
| C/L-tone | XY | 551+ | 753+ | 0.662 | 0.085 | 11.5 | 12.26 | 12.25 | 0.580 |
| C/L-tone | XZ | 1785+ | 2423+ | 0.852 | 0.081 | 17.7 | 12.36 | - | - |
| C/L-tone | YZ | 8209+ | 8864+ | 1.000 | 1.000 | 100.0 | n/a | - | - |
| C/L-b24 | XY | 739+ | 1037+ | 0.776 | 0.117 | 12.9 | 19.69 | 14.72 | 0.766 |
| C/L-b24 | XZ | 2378+ | 3310+ | 0.858 | 0.111 | 24.2 | 14.77 | - | - |
| C/L-b24 | YZ | 8864+ | 8864+ | 1.000 | 1.000 | 100.0 | n/a | - | - |
| C/P-tone | XY | 1308 | 2662+ | 0.762 | 0.327 | 19.5 | 11.47 | - | - |
| C/P-tone | XZ | 1683+ | 2318+ | 0.842 | 0.434 | 28.8 | 11.46 | 11.47 | 0.501 |
| C/P-tone | YZ | 1412+ | 1934+ | 0.597 | 0.783 | 22.0 | 5.99 | - | - |
| C/P-b24 | XY | 2531+ | 4368+ | 0.798 | 0.488 | 32.0 | 14.08 | - | - |
| C/P-b24 | XZ | 2282+ | 3286+ | 0.854 | 0.624 | 40.8 | n/a | 14.07 | 0.401 |
| C/P-b24 | YZ | 1914+ | 2694+ | 0.706 | 0.784 | 30.4 | 8.01 | - | - |
| C/W-tone | XY | 104 | 202 | 0.061 | 0.026 | 16.7 | 0.15 | - | - |
| C/W-tone | XZ | 485+ | 670+ | 0.269 | 0.104 | 25.0 | 0.10 | 0.10 | 0.281 |
| C/W-tone | YZ | 322+ | 514+ | 0.152 | 0.084 | 32.1 | 0.50 | - | - |
| C/W-b24 | XY | 192 | 384 | 0.118 | 0.027 | 17.7 | 3.09 | - | - |
| C/W-b24 | XZ | 657+ | 926+ | 0.372 | 0.095 | 33.4 | 4.97 | 4.97 | 0.349 |
| C/W-b24 | YZ | 496+ | 750+ | 0.188 | 0.086 | 55.3 | 3.50 | - | - |

| Scenario | Focal volume -3 dB (liters) | -6 dB (liters) | ROI shell volume in the box (liters) | Energy share in -6 dB volume, focused | simultaneous | Shell above focus -6 dB, focused (%) | simultaneous (%) | Bounding extent of the -6 dB volume x / y / z (mm) |
|---|---|---|---|---|---|---|---|---|
| C/L-tone | 117+ | 160+ | 1295 | 0.777 | 0.084 | 12.8 | 12.8 | 300 / 850 / 1925 |
| C/L-b24 | 156+ | 219+ | 1295 | 0.816 | 0.113 | 16.9 | 18.6 | 375 / 875 / 1925 |
| C/P-tone | 47+ | 89+ | 1295 | 0.448 | 0.229 | 7.2 | 6.5 | 650 / 575 / 550 |
| C/P-b24 | 85+ | 175+ | 1295 | 0.576 | 0.400 | 13.5 | 12.4 | 900 / 775 / 650 |
| C/W-tone | 4+ | 8+ | 1295 | 0.024 | 0.009 | 15.9 | 15.8 | 275 / 175 / 425 |
| C/W-b24 | 7+ | 15+ | 1295 | 0.044 | 0.010 | 16.6 | 3.4 | 325 / 225 / 425 |

### Series D (extreme focus x = 800 mm, linear array)

| Scenario | Spectrum | Focus (mm) | Range (mm) / angle from beam axis (deg) | Focus level (dB) | Gain (dB) | Width x / y / z at -3 dB (mm) | at -6 dB (mm) |
|---|---|---|---|---|---|---|---|
| D/b1535 | 1.5 kHz to 3.5 kHz | (800.0, 800.0, 0.0) | 994 / 36.4 | 18.12 | +15.41 | 405 / 530+ / 1495+ | 590+ / 615+ / 1495+ |
| D/b24 | 2 kHz to 4 kHz | (800.0, 800.0, 0.0) | 994 / 36.4 | 18.12 | +18.23 | 330 / 465 / 1405 | 480 / 565+ / 1495+ |
| D/b18 | 1 kHz to 8 kHz | (800.0, 800.0, 0.0) | 994 / 36.4 | 18.12 | +16.45 | 215 / 310 / 1095 | 380 / 520+ / 1495+ |

| Scenario | Cut | Area -3 dB (cm2) | Area -6 dB (cm2) | Energy share in -6 dB region, focused | simultaneous | ROI cells above focus -6 dB (%) | Interior secondary contrast (dB) | In-plane arc contrast (dB) | In-plane arc share +/-10 deg |
|---|---|---|---|---|---|---|---|---|---|
| D/b1535 | XY | 1582+ | 2260+ | 0.830 | 0.043 | 25.1 | n/a | 16.15 | 0.691 |
| D/b1535 | XZ | 4902+ | 6488+ | 0.898 | 0.036 | 29.8 | 15.90 | - | - |
| D/b1535 | YZ | 11768+ | 12911+ | 0.995 | 0.989 | 98.1 | n/a | - | - |
| D/b24 | XY | 1310+ | 1843+ | 0.827 | 0.029 | 20.4 | 19.97 | 6.60 | 0.615 |
| D/b24 | XZ | 4124+ | 5614+ | 0.890 | 0.025 | 25.8 | 19.62 | - | - |
| D/b24 | YZ | 10442+ | 12239+ | 0.984 | 0.952 | 93.0 | n/a | - | - |
| D/b18 | XY | 878+ | 1508+ | 0.622 | 0.033 | 16.7 | n/a | 8.73 | 0.486 |
| D/b18 | XZ | 2807+ | 4700+ | 0.686 | 0.034 | 21.6 | 17.57 | - | - |
| D/b18 | YZ | 6814+ | 11614+ | 0.962 | 0.691 | 88.3 | n/a | - | - |

| Scenario | Focal volume -3 dB (liters) | -6 dB (liters) | ROI shell volume in the box (liters) | Energy share in -6 dB volume, focused | simultaneous | Shell above focus -6 dB, focused (%) | simultaneous (%) | Bounding extent of the -6 dB volume x / y / z (mm) |
|---|---|---|---|---|---|---|---|---|
| D/b1535 | 344+ | 483+ | 1782 | 0.854 | 0.038 | 27.1 | 26.8 | 650 / 1000 / 2025 |
| D/b24 | 286+ | 399+ | 1782 | 0.851 | 0.026 | 22.4 | 22.0 | 600 / 975 / 2025 |
| D/b18 | 194+ | 328+ | 1782 | 0.637 | 0.033 | 18.4 | 22.9 | 525 / 950 / 2025 |

## Answers

### 1. Does the linear array produce a well-localized three-dimensional focus?

**No. It produces a ring (a torus), not a point.** Every point of the circle of radius 800 mm around the x axis through x = 350 mm has exactly the same field (verified: the spread of the energy over 13 points of that circle is below 1e-12 relative). The XY cut, the only plane used in Missions 03 to 05, crosses that circle at a single point and hides this. In the other two cuts the structure is plain:

- In the XZ cut at the focus range the focal band is 115 to 180 mm wide in x but at least 1.2 m long in z (limited only by the shell); in the images it is a curved band, an arc of the ring.
- In the YZ cut through the focus (x = 350 mm in series A, 210 mm in A0) the field depends only on the distance to the array axis: concentric rings, and in series A and A0 the whole ROI of that cut lies within 3 dB of the focus level (100 % of the cells).
- The -3 and -6 dB focal volumes are 105 to 229 liters and are truncated by the 2 m box in z: they are lower bounds of a structure that closes on itself, whose circumference is 2 pi x 800 mm = 5 m.

Localization exists in one direction only: along x (the aperture direction), with a width of 115 to 180 mm at -3 dB, of the order of the diffraction width lambda R / D (about 165 mm at 4 kHz). The focusing gain over simultaneous emission (+1.4 to +7.3 dB) is a gain on that ring, not on a point.

### 2. How does the focal region extend perpendicular to the line?

Around the line, as a ring whose radius is the distance of the focus from the array axis: 800 mm in series A, A0 and D, 751.8 mm in series C. Along the ring the field is constant. Its radial thickness is set by the depth of focus: for the foci on the y side (A, A0, C) the -3 dB extent along the beam is at least the 400 mm shell. The extreme lateral focus (series D) has the same structure: the z-width of the -3 dB band is 1.1 to 1.5 m and the YZ cut shows the same concentric rings.

### 3. Does the planar array improve the spatial confinement?

It removes the ring and confines the energy in both transverse directions, but with the 120 mm aperture it does not focus at 800 mm (its spot is about three times wider than the line's x-width).

- **P (3-2-3, 120 mm):** the transverse -3 dB spot is bounded (not truncated): 405 x 345 mm for the tone and 555 x 465 mm for 2-4 kHz, about 3 times the x-width of the linear array. The focusing gain over simultaneous emission is +0.06 dB (tone) and +0.04 dB (2-4 kHz): the focus is so close to the axis of a small aperture that every speaker is practically equidistant from it, so the delays change nothing. The "focus" is a broad beam. Its -6 dB volume (83 to 161 liters, shell-limited in depth) is smaller than the line's ring, but it is still elongated along the beam (extent >= 450 mm).
- **W (3-2-3 scaled to 420 mm):** the transverse spot is 105 to 155 mm in x and 95 to 135 mm in y (all three spectra) and the -6 dB volume is only 7 to 14 liters. It is the only geometry here with a compact three-dimensional peak. But 8 elements spread over 420 mm are 210 mm apart: grating lobes are as strong as the main lobe for the tone (contrast 0.0 to 1.0 dB), 3 to 5 dB for 2-4 kHz and 4 to 6 dB for 1-8 kHz, and the focal volume holds only 2 to 4 % of the ROI energy. The peak is localized; the energy is not.

### 4. Which differences come from the geometry and which from the aperture?

- **Aperture (P to W, same topology, aperture x 3.5):** transverse spot 405 x 345 mm to 115 x 105 mm (tone), consistent with the 3.5 factor; focusing gain +0.06 to +12.9 dB; but the spacing also grows 3.5-fold, which creates the grating lobes. With 8 elements a larger aperture cannot be had without a larger spacing: that trade-off is the aperture effect.
- **Dimensionality (L to W, same largest dimension of 420 mm):** the same x-width (135 mm against 115 mm, tone), but the ring becomes a spot: the second transverse direction goes from at least 1.2 m (L, z) to 95-135 mm (W, y). The grating lobes of W come from its spacing, which L (60 mm) does not have in 2 to 4 kHz.
- **Spacing at equal topology and equal aperture is not tested:** it would need more elements (see below).
- L against P (same 60 mm neighbor distance, aperture 420 against 120 mm) is an aperture comparison, plus the ring.

### 5. Does the 1.5-3.5 kHz band improve the extreme lateral focus?

**Yes for selectivity, with a wider focal region.** At 36.4 degrees of steering the grating-free limit for the 60 mm pitch is 3.59 kHz. The 2-4 kHz and 1-8 kHz bands contain lines above it. In the XY plane of series D (the plane of the arc analysis of Missions 04 and 05):

| Band | In-plane arc contrast (dB) | Share of arc energy within +/-10 deg | Highest other lobe (dB @ deg) | x-width at -3 dB (mm) | Focal volume -6 dB (liters) | Gain (dB) |
|---|---|---|---|---|---|---|
| 1.5-3.5 kHz | 16.1 | 0.691 | 1.97 @ +3.1 | 405 | 483+ | +15.4 |
| 2-4 kHz | 6.6 | 0.615 | 11.5 @ -75.2 | 330 | 399+ | +18.2 |
| 1-8 kHz | 8.7 | 0.486 | 9.4 @ -64.7 | 215 | 329+ | +16.5 |

Mission 05's 1-8 kHz result (contrast 8.73 dB, share 0.486) is reproduced; the 4 kHz tone had 0.11 dB. The lower band removes the grating lobe (the highest competitor is now a weak lobe near broadside, 2 dB at +3 degrees) and raises the energy share from 0.49 to 0.69, but the main lobe is wider (405 mm against 215 mm in x). Nothing changes in the 3-D structure: it is still a ring (z-width 1.1 to 1.5 m, volumes truncated by the box), now larger.

### 6. Which limitations remain even with a planar geometry?

- **Depth of focus.** All apertures here are small compared with the 0.8 to 1 m range: the -3 dB axial extent is at least the 400 mm shell for P and W. A 420 mm aperture at 800 mm does not form a point focus along the beam axis at these wavelengths, it forms a beam that narrows around the focus.
- **Aperture against spacing for a fixed element count.** The compact peak of W comes with grating lobes. A small aperture with fine spacing (P) has no grating lobes and no focus. Eight elements cannot give both.
- **Frequency range.** For W the grating-free range is below 1.6 kHz, where the 420 mm aperture is too small for a compact spot.
- **Model.** Free field, ideal sources, no driver directivity, baffle or room, a finite set of 5 mm cuts and a 25 mm volume, regions truncated by the shell and the box (marked `+`), and indicators that depend on the definitions above.

## Preparing for larger arrays (16, 32 elements)

Nothing was implemented; the following was checked in the code and with tests/benchmarks.

**The simulator has no eight-element limit.** `array_config` accepts any number of elements (positive integer IDs, free XYZ or grid), `focus_delays`, `FieldGeometry`, `band_energy`, `microphone.receive` and the campaign/plot modules loop over `array.elements`. Tested: lines of 16 and 32, a 4 x 4 grid and 32 random positions give the exact coherent level at the focus and run through the virtual microphone (`test_simulator_supports_arbitrary_element_counts`). Eight appears only as a *default* (`--channels 8` and `channels=8` in the raw tone/pulse generator) and in the scenario definitions (`campaign.scenarios`, `field3d.cut_scenarios`), which reference 8-element files. Cost is linear in the number of elements (the 401 x 381 grid, one spectral line):

| Elements | Geometry setup (s) | Time per line (s) | Memory of the distances (MB) |
|---|---|---|---|
| 8 | 0.04 | 0.026 | 20 |
| 16 | 0.08 | 0.049 | 39 |
| 32 | 0.14 | 0.096 | 78 |
| 64 | 0.29 | 0.194 | 156 |

Memory also scales with the number of grid points: a 505 000-point volume with 32 elements needs about 260 MB, with 64 about 0.5 GB. Beyond that, `FieldGeometry` should evaluate the points in blocks. The arc and band analyses scale in the same way. The window-limited indicators and the 400 mm ROI shell are defined with respect to the focus range, not the number of elements.

**The physical limits are different and are not simulator limits.** `hardware_config` and `focus.route_to_pcm` map speakers to the PCM channels of one device (`channels` from the JSON). The C-Media interface has 8 channels, so 16 or 32 speakers need either a larger multichannel interface or several interfaces. Separate USB interfaces have independent clocks: their relative drift and start offset are not controlled and would act as uncontrolled delays of a few microseconds or more (the focusing delays here are of that order). That is a hardware and synchronization design question, not a change of the simulator.

**Recommendations for a future campaign** (to be run with the existing machinery and these indicators):

1. **Regular arrays:** 16 on a line (60 mm: 900 mm aperture) and a 4 x 4 grid (60 mm: 180 x 180 mm, and 120 mm: 360 x 360 mm), plus 32 as a 4 x 8 grid: shows the ring against spot trade-off with the same spacing, which the 8-element campaign could not isolate.
2. **Non-uniform layouts:** a Vogel (sunflower) spiral, concentric rings and pseudo-random layouts with a minimum spacing, to break the periodicity that creates grating lobes. The W result (compact peak, strong grating lobes) is exactly what a non-periodic layout with the same aperture should improve.
3. **Metrics:** keep the ROI-based shares and contrasts. Add a lobe-level criterion for ring-shaped regions (azimuthal energy distribution around the beam axis) because bounding boxes and widths mislead for rings.
4. **Practical:** compute driver sizes before choosing a pitch; start from simulated designs that keep the nearest-neighbor distance above the driver diameter.

## Limitations

- Free-field ideal point sources of equal amplitude; no driver directivity, baffle, diffraction, room, transducer response or mutual coupling.
- Time-averaged energy of an ideal spectrum (see `docs/broadband-study.md` for the difference with a noise realization); no instantaneous pressure.
- The ROI is a modelling choice (0.75 to 1.25 of the focus range, within 80 degrees of the beam axis) and truncates every axial extent to at most about 400 mm; volumes are limited to a 2 m box with a 25 mm grid. All truncated numbers are lower bounds and marked `+`.
- The interior secondary-maximum criterion only sees lobes with a local maximum inside the ROI; the in-plane arc complements it for the planes that hold the aperture center.
- Cuts are sampled on a 5 mm grid (8.6 points per wavelength at 8 kHz); the volumes use 25 mm (adequate for the broadband energy, coarse for a 4 kHz tone: the tone volumes are indicative).
- Only three geometries and one element count; no optimization of the layout.

## Implications for the physical array

- A **line** (the first prototype layout) cannot localize in three dimensions: it focuses on a ring. That is acceptable for experiments in a plane containing the line, but the result must not be read as a point focus. A measurement microphone moved off the XY plane along the ring will see the same level.
- The **3-2-3 file as it is (120 mm)** is a beam former, not a focuser, at 800 mm. Closer foci would be expected to make its delays matter (not simulated here).
- A **wide planar** layout is the only one here that gives a compact peak, at the cost of grating lobes with 8 elements. The next experiment should compare it with non-uniform layouts, not with more regular scalings.
- For the extreme lateral steering of the line (36 degrees), the **1.5-3.5 kHz** band is the best of the three tested.
- None of this is a recommendation of a final geometry: mechanical size, driver directivity and the unmeasured channel mismatch are missing.

## Conclusion

The XY plane used until Mission 05 hid the main fact: the **line array focuses on a ring** around its own axis, not on a point. Along the ring the field is constant (checked numerically to 1e-12) and the focal band is 115 to 180 mm wide in the aperture direction but at least 1.2 m long around it. The 3-2-3 planar array removes the ring but its 120 mm aperture is too small to focus at 800 mm (+0.04 dB over simultaneous emission). Scaling it to 420 mm gives the only compact peak (a spot of about 100 to 150 mm), at the price of strong grating lobes with eight elements, so the peak is localized but the energy is not. For the extreme lateral focus of the line, the 1.5-3.5 kHz band is the best of the three tested (arc contrast 16.1 dB against 6.6 and 8.7 dB). No geometry here confines the energy along the beam within the sampled range. The simulator has no eight-element limit; the limits for larger arrays are those of the hardware.
