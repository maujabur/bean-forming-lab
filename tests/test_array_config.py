import json
from pathlib import Path

import pytest

from bean_lab.array_config import load_array
from bean_lab.hardware_config import load_hardware


ROOT = Path(__file__).resolve().parents[1]


def write(tmp_path, name, payload):
    path = tmp_path / name
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_array_and_hardware_are_independent():
    array = load_array(ROOT / "configs/arrays/linear-4.json")
    hw = load_hardware(ROOT / "configs/hardware/cmedia-0d8c-0102.json")
    assert array.position_for(1) == (0.0, 0.0, 0.0)
    assert hw.pcm_for(1) == 2
    assert hw.pcm_for(8) == 3
    assert len(array.elements) == 4


def test_rejects_duplicate_array_speaker(tmp_path):
    data = {"name": "bad", "units": "mm", "elements": [
        {"speaker": 1, "position": [0, 0, 0]},
        {"speaker": 1, "position": [1, 0, 0]},
    ]}
    with pytest.raises(ValueError, match="duplicate"):
        load_array(write(tmp_path, "array.json", data))


@pytest.mark.parametrize("mapping", [
    {"1": 2, "2": 2},
    {"1": 0},
    {"1": 9},
    {"0": 1},
])
def test_rejects_bad_hardware_mapping(tmp_path, mapping):
    data = {"name": "test", "device": "hw:1,0", "channels": 8, "speaker_to_pcm": mapping}
    with pytest.raises(ValueError):
        load_hardware(write(tmp_path, "hardware.json", data))
