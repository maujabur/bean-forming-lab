"""Load and validate speaker array-to-PCM channel assignments."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Element:
    id: str
    pcm_channel: int  # 1-based, as shown in the CLI
    position: tuple[float, float, float]


@dataclass(frozen=True)
class ArrayConfig:
    name: str
    elements: tuple[Element, ...]

    def pcm_for(self, element_id: str) -> int:
        for element in self.elements:
            if element.id == element_id:
                return element.pcm_channel
        raise KeyError(f"Unknown element: {element_id}")

    def element_for(self, pcm_channel: int) -> Element:
        for element in self.elements:
            if element.pcm_channel == pcm_channel:
                return element
        raise KeyError(f"No element assigned to PCM channel {pcm_channel}")


def load_array(path: str | Path, *, max_channels: int = 8) -> ArrayConfig:
    """Load JSON array configuration. No implicit channel sharing or remapping."""
    with Path(path).open(encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError("array configuration must be a JSON object")
    if obj.get("units") != "meters":
        raise ValueError("units must be 'meters'")
    name = obj.get("name")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("name must be a nonempty string")
    raw_elements = obj.get("elements")
    if not isinstance(raw_elements, list) or not raw_elements:
        raise ValueError("elements must be a nonempty list")
    if type(max_channels) is not int or max_channels < 1:
        raise ValueError("max_channels must be a positive integer")

    elements = []
    ids = set()
    channels = set()
    for idx, raw in enumerate(raw_elements):
        if not isinstance(raw, dict):
            raise ValueError(f"element {idx} must be an object")
        element_id = raw.get("id")
        pcm = raw.get("pcm_channel")
        position = raw.get("position")
        if not isinstance(element_id, str) or not element_id.strip():
            raise ValueError(f"element {idx}: id must be a nonempty string")
        if type(pcm) is not int or not 1 <= pcm <= max_channels:
            raise ValueError(f"{element_id}: pcm_channel must be 1..{max_channels}")
        if not isinstance(position, list) or len(position) != 3 or any(
            type(x) not in (int, float) or not math.isfinite(x) for x in position
        ):
            raise ValueError(f"{element_id}: position must contain three finite numbers")
        if element_id in ids:
            raise ValueError(f"duplicate element id: {element_id}")
        if pcm in channels:
            raise ValueError(f"PCM channel {pcm} assigned more than once")
        ids.add(element_id)
        channels.add(pcm)
        elements.append(Element(element_id, pcm, tuple(float(x) for x in position)))
    return ArrayConfig(name, tuple(elements))
