import json
from pathlib import Path

import pytest

from bean_lab.array_config import load_array


def write(tmp_path, payload):
    path = tmp_path / "array.json"
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def valid():
    return {"name": "test", "units": "meters", "elements": [
        {"id": "A", "pcm_channel": 2, "position": [0, 0, 0]},
        {"id": "B", "pcm_channel": 7, "position": [0.06, 0, 0]},
    ]}


def test_example():
    path = Path(__file__).resolve().parents[1] / "configs/arrays/linear-4.json"
    config = load_array(path)
    assert len(config.elements) == 4
    assert config.pcm_for("SPK1") == 2
    assert config.element_for(7).id == "SPK4"


@pytest.mark.parametrize("change", [
    lambda d: d["elements"][1].update(pcm_channel=2),
    lambda d: d["elements"][1].update(id="A"),
    lambda d: d["elements"][0].update(pcm_channel=0),
    lambda d: d["elements"][0].update(position=[0, 0]),
    lambda d: d.update(units="centimeters"),
    lambda d: d["elements"][0].update(position=[float("nan"), 0, 0]),
])
def test_rejects_invalid(tmp_path, change):
    data = valid()
    change(data)
    with pytest.raises(ValueError):
        load_array(write(tmp_path, data))
