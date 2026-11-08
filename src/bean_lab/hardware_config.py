"""Hardware-specific speaker-to-PCM routing, separate from array geometry."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class HardwareConfig:
    name: str
    device: str
    channels: int
    speaker_to_pcm: dict[int, int]

    def pcm_for(self, speaker: int) -> int:
        try:
            return self.speaker_to_pcm[speaker]
        except KeyError:
            raise KeyError(f"Speaker {speaker} has no hardware mapping") from None


def load_hardware(path: str | Path) -> HardwareConfig:
    with Path(path).open(encoding="utf-8") as f:
        obj = json.load(f)
    if not isinstance(obj, dict):
        raise ValueError("hardware configuration must be an object")
    name, device, channels = obj.get("name"), obj.get("device"), obj.get("channels")
    if not isinstance(name, str) or not name.strip():
        raise ValueError("hardware name must be a nonempty string")
    if not isinstance(device, str) or not device.strip():
        raise ValueError("device must be a nonempty ALSA device name")
    if type(channels) is not int or channels < 1:
        raise ValueError("channels must be a positive integer")
    raw = obj.get("speaker_to_pcm")
    if not isinstance(raw, dict) or not raw:
        raise ValueError("speaker_to_pcm must be a nonempty object")
    mapping: dict[int, int] = {}
    used_pcm = set()
    for key, pcm in raw.items():
        if not key.isdecimal() or int(key) < 1:
            raise ValueError(f"invalid speaker ID: {key}")
        speaker = int(key)
        if speaker in mapping:
            raise ValueError(f"duplicate speaker ID: {speaker}")
        if type(pcm) is not int or not 1 <= pcm <= channels:
            raise ValueError(f"speaker {speaker}: PCM channel must be 1..{channels}")
        if pcm in used_pcm:
            raise ValueError(f"PCM channel {pcm} assigned to multiple speakers")
        mapping[speaker] = pcm
        used_pcm.add(pcm)
    return HardwareConfig(name, device, channels, mapping)
