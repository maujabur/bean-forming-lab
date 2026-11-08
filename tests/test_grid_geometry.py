import json
from pathlib import Path

import pytest

from bean_lab.array_config import load_array

ROOT = Path(__file__).resolve().parents[1]


def load_payload(tmp_path, data):
    path = tmp_path / "array.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return load_array(path)


def base():
    return {"name": "test", "units": "mm", "grid": {"pitch": 30},
            "elements": [{"speaker": 1, "cell": [0, 0]},
                         {"speaker": 2, "cell": [2, 1]}]}


def test_staggered_geometry():
    array = load_array(ROOT / "configs/arrays/staggered-3-2-3.json")
    assert len(array.elements) == 8
    assert array.position_for(1) == (0.0, 0.0, 0.0)
    assert array.position_for(4) == (0.03, 0.06, 0.0)
    assert array.position_for(8) == (0.12, 0.12, 0.0)


def test_free_mm_geometry():
    array = load_array(ROOT / "configs/arrays/linear-4.json")
    assert array.position_for(1) == (0.0, 0.0, 0.0)
    assert array.position_for(4) == (0.18, 0.0, 0.0)


@pytest.mark.parametrize(("grid", "expected"), [
    ({"pitch": 30}, (0.06, 0.03, 0.0)),
    ({"pitch": 30, "pitch_y": 40}, (0.06, 0.04, 0.0)),
    ({"pitch": 30, "pitch_x": 25, "pitch_y": 40}, (0.05, 0.04, 0.0)),
    ({"pitch_x": 25, "pitch_y": 40}, (0.05, 0.04, 0.0)),
])
def test_pitch_resolution(tmp_path, grid, expected):
    data = base()
    data["grid"] = grid
    assert load_payload(tmp_path, data).position_for(2) == expected


@pytest.mark.parametrize("change", [
    lambda d: d["elements"][1].update(cell=[0, 0]),
    lambda d: d["elements"][1].update(cell=[1.5, 0]),
    lambda d: d["elements"][1].update(cell=[-1, 0]),
    lambda d: d["elements"][1].update(cell=[True, 0]),
    lambda d: d["grid"].update(pitch=0),
    lambda d: d["grid"].update(pitch=float("nan")),
    lambda d: d.update(grid={"pitch_x": 25}),
    lambda d: d["elements"][1].update(speaker=1),
    lambda d: d["elements"][1].update(position=[0, 0, 0]),
    lambda d: d["grid"].update(origin_mm=[0, 0, 0]),
])
def test_invalid_grid(tmp_path, change):
    data = base()
    change(data)
    with pytest.raises(ValueError):
        load_payload(tmp_path, data)


def test_rejects_cell_without_grid(tmp_path):
    data = base()
    del data["grid"]
    with pytest.raises(ValueError):
        load_payload(tmp_path, data)


def test_rejects_old_units(tmp_path):
    data = base()
    data["units"] = "meters"
    with pytest.raises(ValueError):
        load_payload(tmp_path, data)
