"""Speaker geometry, independent of audio interface wiring.

Supports explicit XYZ meters and integer-grid positions with a pitch in mm.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Element:
    speaker: int
    position: tuple[float, float, float]


@dataclass(frozen=True)
class ArrayConfig:
    name: str
    elements: tuple[Element, ...]

    def position_for(self, speaker: int) -> tuple[float, float, float]:
        for element in self.elements:
            if element.speaker == speaker:
                return element.position
        raise KeyError(f"Speaker {speaker} not in array")


def load_array(path: str | Path) -> ArrayConfig:
    with Path(path).open(encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError("array must be a JSON object")
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("array name must be a nonempty string")

    grid = obj.get("grid")
    if grid is not None:
        if obj.get("units") != "millimeters" or not isinstance(grid, dict):
            raise ValueError("grid arrays require units='millimeters' and a grid object")
        pitch = grid.get("pitch_mm")
        if type(pitch) not in (int, float) or not math.isfinite(pitch) or pitch <= 0:
            raise ValueError("grid.pitch_mm must be positive and finite")
        if "origin_mm" in grid:
            origin = grid["origin_mm"]
            if not isinstance(origin, list) or len(origin) != 3 or any(
                type(v) not in (int, float) or not math.isfinite(v) for v in origin
            ):
                raise ValueError("grid.origin_mm must contain three finite numbers")
        else:
            origin = [0, 0, 0]
    elif obj.get("units") != "meters":
        raise ValueError("array must use units='meters' or a millimeter grid")

    raw_elements = obj.get("elements")
    if not isinstance(raw_elements, list) or not raw_elements:
        raise ValueError("elements must be a nonempty list")
    elements = []
    speakers = set()
    occupied = set()
    for idx, raw in enumerate(raw_elements):
        if not isinstance(raw, dict):
            raise ValueError(f"element {idx} must be an object")
        speaker = raw.get("speaker")
        if type(speaker) is not int or speaker < 1:
            raise ValueError(f"element {idx}: speaker must be a positive integer")
        if speaker in speakers:
            raise ValueError(f"duplicate speaker {speaker}")
        if grid is not None:
            cell = raw.get("cell")
            if not isinstance(cell, list) or len(cell) != 2 or any(type(v) is not int for v in cell):
                raise ValueError(f"speaker {speaker}: cell must contain two integers")
            if tuple(cell) in occupied:
                raise ValueError(f"duplicate grid cell {cell}")
            occupied.add(tuple(cell))
            position = ((origin[0] + cell[0] * pitch) / 1000,
                        (origin[1] + cell[1] * pitch) / 1000,
                        origin[2] / 1000)
        else:
            position = raw.get("position")
            if not isinstance(position, list) or len(position) != 3 or any(
                type(x) not in (int, float) or not math.isfinite(x) for x in position
            ):
                raise ValueError(f"speaker {speaker}: position must contain three finite numbers")
        speakers.add(speaker)
        elements.append(Element(speaker, tuple(float(x) for x in position)))
    return ArrayConfig(name, tuple(elements))
