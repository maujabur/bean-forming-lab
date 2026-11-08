"""Speaker geometry, independent of audio interface wiring.

Both free XYZ millimeters and integer-grid layouts normalize to XYZ meters.
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


def _positive_pitch(value: object, label: str) -> float:
    if type(value) not in (int, float) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{label} must be a finite positive number")
    return float(value)


def load_array(path: str | Path) -> ArrayConfig:
    with Path(path).open(encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict) or obj.get("units") != "mm":
        raise ValueError("array must be a JSON object with units='mm'")
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("array name must be a nonempty string")

    grid = obj.get("grid")
    if grid is not None:
        if not isinstance(grid, dict):
            raise ValueError("grid must be an object")
        if "origin_mm" in grid:
            raise ValueError("grid origin is fixed at [0, 0]; omit origin_mm")
        default = grid.get("pitch")
        pitch_x = _positive_pitch(grid.get("pitch_x", default), "pitch_x")
        pitch_y = _positive_pitch(grid.get("pitch_y", default), "pitch_y")

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
            if "position" in raw:
                raise ValueError(f"speaker {speaker}: cannot mix cell and position")
            cell = raw.get("cell")
            if not isinstance(cell, list) or len(cell) != 2 or any(
                type(v) is not int or v < 0 for v in cell
            ):
                raise ValueError(f"speaker {speaker}: cell must contain two nonnegative integers")
            if tuple(cell) in occupied:
                raise ValueError(f"duplicate grid cell {cell}")
            occupied.add(tuple(cell))
            position = (cell[0] * pitch_x / 1000, cell[1] * pitch_y / 1000, 0.0)
        else:
            if "cell" in raw:
                raise ValueError(f"speaker {speaker}: cell requires grid")
            position = raw.get("position")
            if not isinstance(position, list) or len(position) != 3 or any(
                type(x) not in (int, float) or not math.isfinite(x) for x in position
            ):
                raise ValueError(f"speaker {speaker}: position must contain three finite numbers")
            position = tuple(float(v) / 1000 for v in position)
        speakers.add(speaker)
        elements.append(Element(speaker, tuple(float(x) for x in position)))
    return ArrayConfig(name, tuple(elements))
