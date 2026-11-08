"""Speaker geometry, independent of audio interface wiring."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Element:
    speaker: int  # stable 1-based physical speaker identity
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
    if not isinstance(obj, dict) or obj.get("units") != "meters":
        raise ValueError("array must be a JSON object with units='meters'")
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("array name must be a nonempty string")
    raw_elements = obj.get("elements")
    if not isinstance(raw_elements, list) or not raw_elements:
        raise ValueError("elements must be a nonempty list")
    elements = []
    speakers = set()
    for idx, raw in enumerate(raw_elements):
        if not isinstance(raw, dict):
            raise ValueError(f"element {idx} must be an object")
        speaker = raw.get("speaker")
        position = raw.get("position")
        if type(speaker) is not int or speaker < 1:
            raise ValueError(f"element {idx}: speaker must be a positive integer")
        if speaker in speakers:
            raise ValueError(f"duplicate speaker {speaker}")
        if not isinstance(position, list) or len(position) != 3 or any(
            type(x) not in (int, float) or not math.isfinite(x) for x in position
        ):
            raise ValueError(f"speaker {speaker}: position must contain three finite numbers")
        speakers.add(speaker)
        elements.append(Element(speaker, tuple(float(x) for x in position)))
    return ArrayConfig(name, tuple(elements))
