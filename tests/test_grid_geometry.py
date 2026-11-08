import json
from pathlib import Path

import pytest

from bean_lab.array_config import load_array


ROOT = Path(__file__).resolve().parents[1]


def test_staggered_geometry():
    array = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    assert len(array.elements) == 8
    assert array.position_for(1) == (0.0, 0.0, 0.0)
    assert array.position_for(4) == (0.03, 0.06, 0.0)
    assert array.position_for(8) == (0.12, 0.12, 0.0)


def test_legacy_meter_geometry():
    array = load_array(ROOT / "configs/arrays/linear-4.json")
    assert array.position_for(1) == (-0.09, 0.0, 0.0)


@pytest.mark.parametrize("change", [
    lambda d: d["elements"][1].update(cell=[0, 0]),
    lambda d: d["elements"][1].update(cell=[1.5, 0]),
    lambda d: d["grid"].update(pitch_mm=0),
    lambda d: d["grid"].update(pitch_mm=float("nan")),
    lambda d: d["elements"][1].update(speaker=1),
])
def test_grid_rejects_invalid(tmp_path, change):
    data = {"name": "test", "units": "millimeters",
            "grid": {"pitch_mm": 30},
            "elements": [{"speaker": 1, "cell": [0, 0]},
                         {"speaker": 2, "cell": [2, 0]}]}
    change(data)
    path = tmp_path / "test.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    with pytest.raises(ValueError):
        load_array(path)
