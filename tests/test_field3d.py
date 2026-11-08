import json
import math
from pathlib import Path

import numpy as np
import pytest

from bean_lab.array_config import ArrayConfig, Element, load_array
from bean_lab.broadband import Spectrum, band_energy, effective_level_db, map_energies
from bean_lab.campaign import BAND, Params, Region, run_scenario, scenarios
from bean_lab.cli import main
from bean_lab.field import FieldGeometry, complex_field, field_maps, grid_axis
from bean_lab.field3d import (analyze_cut, AXIS, PLANES, CutScenario, _run, axis_values, component_nd, cut_scenarios, default_windows,
                              dilate, energies, lateral_focus, local_maxima, plane_points, run_cuts, volume_points,
                              broadside_focus, build_comparison_plate, build_cut_plate)
from bean_lab.focus import SPEED_OF_SOUND, focus_delays
from bean_lab.microphone import receive

ROOT = Path(__file__).resolve().parents[1]
ARRAYS = ROOT / "configs/arrays"
LINE8 = ARRAYS / "line-8-60mm.json"
PLANAR = ARRAYS / "staggered-3-2-3.json"
WIDE = ARRAYS / "staggered-3-2-3-wide.json"
C = SPEED_OF_SOUND
A = 0.05
PARAMS = Params()
STEP = 50.0  # coarse cuts keep the tests fast; every campaign focus is on a 50 mm node


def scenario(array, focus, spectrum, front="y"):
    return CutScenario("T", "t", "t", array.name, focus, spectrum, front)


def pts_of(array, plane, focus, step=STEP, front=None):
    front = front or ("y" if "line" in array.name else "z")
    return plane_points(plane, focus, default_windows(focus, front), step)


# ---- geometry facts of the arrays -----------------------------------------------------------------

def test_array_facts_the_campaign_depends_on():
    lin, pla, wide = (load_array(p) for p in (LINE8, PLANAR, WIDE))
    for arr, center, extent in ((lin, (210, 0, 0), (420, 0, 0)), (pla, (60, 60, 0), (120, 120, 0)),
                                (wide, (210, 210, 0), (420, 420, 0))):
        p = np.array([e.position for e in arr.elements]) * 1000
        assert p.mean(0) == pytest.approx(center) and np.ptp(p, axis=0) == pytest.approx(extent)
    assert all(e.position[2] == 0 for e in pla.elements)  # the 3-2-3 file lies in the z = 0 plane: its normal is +z
    p_pla = np.array([e.position for e in pla.elements])
    p_wide = np.array([e.position for e in wide.elements])
    assert np.allclose(p_wide, 3.5 * p_pla)  # W is exactly the 3-2-3 file scaled by 3.5
    assert broadside_focus("L") == (210.0, 800.0, 0.0)
    assert broadside_focus("P") == pytest.approx((60.0, 60.0, 800.0))
    assert broadside_focus("W") == pytest.approx((210.0, 210.0, 800.0))


def test_lateral_foci_keep_the_range_and_the_angle():
    for g in ("L", "P", "W"):
        arr = load_array(ARRAYS / {"L": "line-8-60mm.json", "P": "staggered-3-2-3.json", "W": "staggered-3-2-3-wide.json"}[g])
        c = np.mean([e.position for e in arr.elements], axis=0) * 1000
        f = np.array(lateral_focus(g))
        assert np.linalg.norm(f - c) == pytest.approx(800.0)
        front = 1 if g == "L" else 2
        assert math.degrees(math.acos((f - c)[front] / 800)) == pytest.approx(20.0)
        assert (f - c)[0] == pytest.approx(800 * math.sin(math.radians(20)))


def test_series_definitions():
    assert [s.key for s in cut_scenarios("A")] == ["tone", "b24", "b18"]
    assert cut_scenarios("A")[0].focus_mm == (350.0, 800.0, 0.0)
    assert cut_scenarios("A0")[0].focus_mm == (210.0, 800.0, 0.0)
    assert cut_scenarios("B")[0].focus_mm == pytest.approx((60.0, 60.0, 800.0)) and cut_scenarios("B")[0].front_axis == "z"
    d = cut_scenarios("D")
    assert [(s.spectrum.f_min, s.spectrum.f_max) for s in d] == [(1500, 3500), (2000, 4000), (1000, 8000)]
    assert {s.focus_mm for s in d} == {(800.0, 800.0, 0.0)}
    assert len(cut_scenarios("C")) == 6
    with pytest.raises(ValueError):
        cut_scenarios("Q")


# ---- planes ------------------------------------------------------------------------------------------------

def test_planes_pass_through_the_focus_and_use_project_coordinates():
    focus = (350.0, 800.0, 0.0)
    w = default_windows(focus, "y")
    assert w == {"x": (-650.0, 1350.0), "y": (100.0, 2000.0), "z": (-1000.0, 1000.0)}
    for name, (ua, va, na) in PLANES.items():
        u, v, pts = plane_points(name, focus, w, 50.0)
        assert pts.shape == (len(v), len(u), 3)
        assert np.all(pts[..., AXIS[na]] == focus[AXIS[na]] / 1000)  # constant normal coordinate: the focus value
        assert np.allclose(pts[0, :, AXIS[ua]] * 1000, u) and np.allclose(pts[:, 0, AXIS[va]] * 1000, v)
        assert focus[AXIS[ua]] in u and focus[AXIS[va]] in v  # the focus is a node of the cut


def test_windows_and_steps_are_validated():
    with pytest.raises(ValueError):
        axis_values(0, 100, 30)
    with pytest.raises(ValueError):
        axis_values(0, 100, 0)
    with pytest.raises(ValueError):
        axis_values(100, 0, 10)
    with pytest.raises(ValueError):
        default_windows((0, 0, 0), "x")
    with pytest.raises(ValueError):
        plane_points("AB", (0, 0, 0), default_windows((0, 800, 0)), 50)


def test_point_on_a_speaker_is_rejected():
    arr = load_array(LINE8)
    windows = {"x": (-100.0, 500.0), "y": (-100.0, 100.0), "z": (-100.0, 100.0)}
    u, v, pts = plane_points("XY", (0.0, 0.0, 0.0), windows, 20.0)  # the grid hits (0, 0, 0) = speaker 1
    with pytest.raises(ValueError, match="coincides"):
        energies(arr, (0.0, 0.5, 0.0), Spectrum(4000, 4000), pts, PARAMS)


# ---- equivalence with earlier missions ----------------------------------------------------------------------------

def test_xy_cut_equals_the_mission_03_map_and_the_mission_05_energy_map():
    arr = load_array(LINE8)
    focus = (350.0, 800.0, 0.0)
    windows = {"x": (-700.0, 1300.0), "y": (100.0, 2000.0), "z": (-1000.0, 1000.0)}
    u, v, pts = plane_points("XY", focus, windows, 50.0)
    e_s, e_f = energies(arr, tuple(c / 1000 for c in focus), Spectrum(4000, 4000), pts, PARAMS)
    fm = field_maps(arr, (0.35, 0.8, 0.0), 4000.0, grid_axis(-0.7, 1.3, 41, "x"), grid_axis(0.1, 2.0, 39, "y"), 0.0, gain=A)
    assert np.allclose(e_s, np.abs(fm.simultaneous) ** 2, rtol=1e-12) and np.allclose(e_f, np.abs(fm.focused) ** 2, rtol=1e-12)
    band = Spectrum(2000, 4000, df=100)
    b_s, b_f = energies(arr, (0.35, 0.8, 0.0), band, pts, PARAMS)
    m_s, m_f = map_energies(arr, (0.35, 0.8, 0.0), band, grid_axis(-0.7, 1.3, 41, "x"), grid_axis(0.1, 2.0, 39, "y"), 0.0, PARAMS)
    assert np.allclose(b_s, m_s, rtol=1e-12) and np.allclose(b_f, m_f, rtol=1e-12)


def test_focus_levels_reproduce_missions_03_04_05():
    r = run_cuts(scenario(LINE8, (350.0, 800.0, 0.0), Spectrum(4000, 4000)), PARAMS, step_mm=10.0, array_dir=ARRAYS)
    m = run_scenario(scenarios("F")[3], Region(nx=41, ny=39), Params(oversample=1), ARRAYS)
    for k in ("level_sim_focus_db", "level_foc_focus_db", "focusing_gain_db"):
        assert r.indicators[k] == pytest.approx(m.indicators[k], abs=1e-9)
    # secondary-lobe contrast in the cut agrees with the arc analysis of Mission 04 (two independent definitions)
    assert r.indicators["XY"]["contrast_db"] == pytest.approx(m.indicators["arc_foc_contrast_db"], abs=0.6)
    assert r.indicators["XY"]["secondary_db"] < r.indicators["level_foc_focus_db"] - 10
    # cut plane values equal the Mission 04 maps where the grids meet
    ix, iy = int(np.argmin(np.abs(r.planes["XY"].u_mm - 350))), int(np.argmin(np.abs(r.planes["XY"].v_mm - 800)))
    assert r.planes["XY"].lev_foc[iy, ix] == pytest.approx(m.indicators["level_foc_focus_db"], abs=1e-9)


def test_cut_values_equal_direct_point_evaluation_in_xyz():
    arr = load_array(PLANAR)
    focus = (60.0, 60.0, 800.0)
    spec = Spectrum(2000, 4000, df=200)
    f, w = spec.lines()
    em = {d.speaker: d.delay for d in focus_delays(arr, tuple(c / 1000 for c in focus))}
    rng = np.random.default_rng(5)
    for name in PLANES:
        u, v, pts = pts_of(arr, name, focus, 100.0)
        e_s, e_f = energies(arr, tuple(c / 1000 for c in focus), spec, pts, PARAMS)
        for _ in range(6):
            iv, iu = int(rng.integers(len(v))), int(rng.integers(len(u)))
            xyz = pts[iv, iu]  # an explicit XYZ point
            direct_s = sum(wk * abs(complex_field(arr, xyz, fk, gain=A)) ** 2 for fk, wk in zip(f, w))
            direct_f = sum(wk * abs(complex_field(arr, xyz, fk, emission_delays=em, gain=A)) ** 2 for fk, wk in zip(f, w))
            assert e_s[iv, iu] == pytest.approx(direct_s, rel=1e-12)
            assert e_f[iv, iu] == pytest.approx(direct_f, rel=1e-12)


def test_cuts_are_slices_of_the_3d_grid():
    arr = load_array(LINE8)
    focus = (350.0, 800.0, 0.0)
    windows = default_windows(focus, "y")
    xs, vpts = volume_points(focus, windows, 100.0)
    spec = Spectrum(3000, 3000)
    ev_s, ev_f = energies(arr, tuple(c / 1000 for c in focus), spec, vpts, PARAMS)
    iz = int(np.argmin(np.abs(xs[2] - 0.0)))
    iy = int(np.argmin(np.abs(xs[1] - 800.0)))
    ix = int(np.argmin(np.abs(xs[0] - 350.0)))
    for name, sl in (("XY", ev_f[:, :, iz].T), ("XZ", ev_f[:, iy, :].T), ("YZ", ev_f[ix, :, :].T)):
        u, v, pts = plane_points(name, focus, windows, 100.0)
        _, e_f = energies(arr, tuple(c / 1000 for c in focus), spec, pts, PARAMS)
        assert np.allclose(sl, e_f, rtol=1e-12), name


# ---- symmetry ------------------------------------------------------------------------------------------------------------

def test_linear_array_field_is_axially_symmetric_so_the_focus_is_a_ring():
    arr = load_array(LINE8)
    focus_m = (0.35, 0.8, 0.0)
    spec = Spectrum(1000, 8000, df=500)
    rho, x = 0.8, 0.35
    phi = np.linspace(0, 2 * np.pi, 13)
    ring = np.stack([np.full_like(phi, x), rho * np.cos(phi), rho * np.sin(phi)], axis=-1)  # circle around the x axis
    e_s, e_f = energies(arr, focus_m, spec, ring, PARAMS)
    assert np.ptp(e_f) < 1e-12 * e_f.mean() and np.ptp(e_s) < 1e-12 * e_s.mean()
    # in the YZ cut the energy depends only on the distance to the x axis
    u, v, pts = plane_points("YZ", (350.0, 800.0, 0.0), default_windows((350.0, 800.0, 0.0)), 100.0)
    _, ef = energies(arr, focus_m, spec, pts, PARAMS)
    iz0 = int(np.argmin(np.abs(v)))
    for iz in (3, 7, 12):
        for iy in (2, 9, 14):
            rho_cell = math.hypot(u[iy], v[iz])
            if u[0] <= rho_cell <= u[-1]:
                jy = int(np.argmin(np.abs(u - rho_cell)))
                if abs(u[jy] - rho_cell) < 1e-9:
                    assert ef[iz, iy] == pytest.approx(ef[iz0, jy], rel=1e-10)
    r = run_cuts(scenario(LINE8, (350.0, 800.0, 0.0), spec), PARAMS, step_mm=STEP, array_dir=ARRAYS)
    z_width = r.indicators["axis_width_mm"]["3"]["z"]
    x_width = r.indicators["axis_width_mm"]["3"]["x"]
    assert r.indicators["axis_width_limited"]["3"]["z"] and z_width > 4 * x_width  # extended along the ring, not a point


def test_planar_array_has_x_and_y_mirror_symmetry():
    arr = load_array(PLANAR)
    focus = (60.0, 60.0, 800.0)  # on the symmetry axes of the 3-2-3 file
    spec = Spectrum(2000, 4000, df=250)
    focus_m = tuple(c / 1000 for c in focus)
    u, v, pts = pts_of(arr, "XY", focus, 50.0)
    e_s, e_f = energies(arr, focus_m, spec, pts, PARAMS)
    assert np.allclose(e_f, e_f[:, ::-1], rtol=1e-9) and np.allclose(e_f, e_f[::-1, :], rtol=1e-9)
    assert np.allclose(e_s, e_s[:, ::-1], rtol=1e-9) and np.allclose(e_s, e_s[::-1, :], rtol=1e-9)
    # the two axial cuts are NOT equivalent: the layout has no 90-degree symmetry
    u, v, pxz = pts_of(arr, "XZ", focus, 50.0)
    u2, v2, pyz = pts_of(arr, "YZ", focus, 50.0)
    _, f_xz = energies(arr, focus_m, spec, pxz, PARAMS)
    _, f_yz = energies(arr, focus_m, spec, pyz, PARAMS)
    assert f_xz.shape == f_yz.shape and not np.allclose(f_xz, f_yz, rtol=1e-3)


def test_wide_array_mirror_symmetry_and_mirrored_lateral_foci():
    arr = load_array(WIDE)
    spec = Spectrum(2000, 2000)
    f1, f2 = (210.0 + 300.0, 210.0, 800.0), (210.0 - 300.0, 210.0, 800.0)  # mirror images through x = 210
    w = {"x": (-790.0, 1210.0), "y": (-790.0, 1210.0), "z": (100.0, 2000.0)}  # symmetric about x = y = 210
    for plane in ("XY",):
        u, v, p1 = plane_points(plane, f1, w, 50.0)
        _, _, p2 = plane_points(plane, f2, w, 50.0)
        a1 = energies(arr, tuple(c / 1000 for c in f1), spec, p1, PARAMS)[1]
        a2 = energies(arr, tuple(c / 1000 for c in f2), spec, p2, PARAMS)[1]
        assert np.allclose(a1, a2[:, ::-1], rtol=1e-9)


# ---- delays, normalization ------------------------------------------------------------------------------------------------------

@pytest.mark.parametrize("path,focus", [(LINE8, (350.0, 800.0, 0.0)), (PLANAR, (60.0, 60.0, 800.0)), (WIDE, (210.0, 210.0, 800.0)),
                                        (PLANAR, (333.7, 60.0, 751.8))])
def test_focus_is_coherent_for_every_spectrum(path, focus):
    arr = load_array(path)
    coherent = 20 * math.log10(sum(1 / math.dist([c / 1000 for c in focus], e.position) for e in arr.elements))
    for spec in (Spectrum(4000, 4000), Spectrum(2000, 4000, df=200), Spectrum(1000, 8000, df=500), Spectrum(1500, 3500, df=250)):
        front = "y" if "line" in path.name else "z"
        r = run_cuts(CutScenario("T", "t", "t", path.name, focus, spec, front), PARAMS, step_mm=100.0, array_dir=ARRAYS)
        assert r.indicators["level_foc_focus_db"] == pytest.approx(coherent, abs=1e-9)
        assert r.indicators["weights_sum"] == pytest.approx(1.0, abs=1e-14)


def test_spectral_normalization_on_cut_points():
    arr = ArrayConfig("one", (Element(1, (0.0, 0.0, 0.0)),))
    pts = np.array([[[0.0, 2.0, 0.0], [0.5, 1.0, 0.3]]])
    for spec in (Spectrum(4000, 4000), Spectrum(2000, 4000, df=100), Spectrum(1000, 8000, df=500)):
        e_s, _ = energies(arr, (0.0, 1.0, 0.0), spec, pts, PARAMS)
        r = np.linalg.norm(pts, axis=-1)
        assert np.allclose(e_s, (A / r) ** 2, rtol=1e-12)  # one speaker: same energy for every spectrum


# ---- indicators -------------------------------------------------------------------------------------------------------------------

def test_single_speaker_makes_the_whole_shell_the_focal_region():
    arr = ArrayConfig("one", (Element(1, (0.0, 0.0, 0.0)),))
    focus = (0.0, 800.0, 0.0)
    sc = CutScenario("T", "one", "one", "one.json", focus, Spectrum(1000, 3000, df=500), "y")
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "one.json").write_text(json.dumps({"name": "one", "units": "mm", "elements": [{"speaker": 1, "position": [0, 0, 0]}]}))
        r = run_cuts(sc, PARAMS, step_mm=20.0, array_dir=Path(d), volume=True, volume_step_mm=50.0)
    i = r.indicators["XY"]
    u, v = r.planes["XY"].u_mm, r.planes["XY"].v_mm
    uu, vv = np.meshgrid(u, v)
    rr = np.hypot(uu, vv)
    shell = (rr >= 600) & (rr <= 1000)
    # no directivity: after compensating 1/r every ROI cell is within 0 dB of the focus level
    assert i["roi_cells"] == int(shell.sum()) and i["area_3_cm2"] == pytest.approx(shell.sum() * 400 / 100)
    assert i["roi_pct_above_m3_foc"] == pytest.approx(100.0) and i["contrast_db"] is None
    assert r.indicators["focusing_gain_db"] == pytest.approx(0.0, abs=1e-9)
    vol = r.volume
    assert vol["volume_3_liters"] == pytest.approx(vol["shell_volume_liters"])  # whole shell; area never reported as volume
    assert vol["truncated_3"] and vol["energy_share_in_3_foc"] == pytest.approx(1.0)


def test_indicator_structure_and_consistency():
    r = run_cuts(scenario(PLANAR, (60.0, 60.0, 800.0), Spectrum(2000, 4000, df=500), "z"), PARAMS, step_mm=STEP, array_dir=ARRAYS)
    i = r.indicators
    for name in PLANES:
        c = i[name]
        assert c["area_6_cm2"] >= c["area_3_cm2"] > 0
        assert c["energy_share_in_6_foc"] >= c["energy_share_in_3_foc"]
        assert 0 <= c["energy_share_in_3_sim"] <= 1 and 0 <= c["energy_share_in_3_foc"] <= 1
        assert c["roi_pct_above_m3_foc"] <= c["roi_pct_above_m6_foc"] <= c["roi_pct_above_m10_foc"] <= 100
        assert c["u_extent_3_mm"] <= c["u_extent_6_mm"]
        assert c["width_u_3_mm"] <= c["width_u_6_mm"] and c["width_v_3_mm"] <= c["width_v_6_mm"]
        if c["contrast_db"] is not None:
            assert c["secondary_db"] == pytest.approx(i["level_foc_focus_db"] - c["contrast_db"])
    # the same axis line appears in two cuts: widths must agree
    assert i["XY"]["width_u_3_mm"] == i["XZ"]["width_u_3_mm"]  # x axis line
    assert i["XY"]["width_v_3_mm"] == i["YZ"]["width_u_3_mm"]  # y axis line
    assert i["XZ"]["width_v_3_mm"] == i["YZ"]["width_v_3_mm"]  # z axis line
    assert set(i["axis_width_mm"]["3"]) == {"x", "y", "z"}


def test_focal_widths_scale_with_aperture_and_wavelength():
    w = {}
    for key, path, focus in (("P", PLANAR, (60.0, 60.0, 800.0)), ("W", WIDE, (210.0, 210.0, 800.0))):
        r = run_cuts(CutScenario("T", key, key, path.name, focus, Spectrum(2000, 2000), "z"), PARAMS, step_mm=10.0, array_dir=ARRAYS)
        w[key] = r.indicators["axis_width_mm"]["3"]["x"]
    assert w["P"] > 2.5 * w["W"]  # 3.5x wider aperture: much narrower focal spot (diffraction: ~ lambda R / D)


def test_reproducible_and_independent_of_display_scales():
    sc = scenario(LINE8, (350.0, 800.0, 0.0), Spectrum(2000, 4000, df=500))
    a = run_cuts(sc, PARAMS, step_mm=STEP, array_dir=ARRAYS)
    b = run_cuts(sc, Params(vmax=3.0, range_db=7.0, diff_limit=2.0), step_mm=STEP, array_dir=ARRAYS)
    for name in PLANES:
        for attr in ("e_sim", "e_foc", "lev_sim", "lev_foc"):
            assert np.array_equal(getattr(a.planes[name], attr), getattr(b.planes[name], attr))
        assert np.array_equal(a.planes[name].diff, b.planes[name].diff, equal_nan=True)
    assert a.indicators == b.indicators


# ---- helpers ---------------------------------------------------------------------------------------------------------------------------

def test_local_maxima_dilate_component_run():
    a = np.array([[0, 1, 0, 0], [1, 5, 1, 0], [0, 1, 0, 3], [0, 0, 0, 3]], dtype=float)
    mx = local_maxima(a)
    assert mx[1, 1] and mx[2, 3] ^ mx[3, 3] and mx.sum() == 2  # plateau counted once
    m = np.zeros((5, 5), dtype=bool)
    m[2, 2] = True
    assert dilate(m, 1).sum() == 9 and dilate(m, 2).sum() == 25
    vol = np.zeros((4, 4, 4), dtype=bool)
    vol[1:3, 1:3, 1:3] = True
    vol[0, 0, 0] = True
    comp = component_nd(vol, (1, 1, 1))
    assert comp.sum() == 8 and not comp[0, 0, 0]
    assert component_nd(vol, (3, 3, 3)).sum() == 0
    above = np.array([0, 1, 1, 1, 0, 1], dtype=bool)
    ok = np.array([1, 1, 1, 1, 1, 1], dtype=bool)
    assert _run(above, ok, 2) == (1, 3, False, False)  # ended by falling below the threshold
    assert _run(above, np.array([1, 1, 1, 0, 1, 1], dtype=bool), 2) == (1, 2, False, True)  # ended by the ROI shell
    assert _run(np.ones(4, dtype=bool), np.ones(4, dtype=bool), 1) == (0, 3, True, True)  # window edges


# ---- plates, CLI, scaling -----------------------------------------------------------------------------------------------------------------

def test_cut_plate_and_comparison_plate_share_scales():
    params = Params()
    rs = [run_cuts(sc, params, step_mm=100.0, array_dir=ARRAYS) for sc in (
        CutScenario("T", "l", "L tone", LINE8.name, (210.0, 800.0, 0.0), Spectrum(4000, 4000), "y"),
        CutScenario("T", "p", "P tone", PLANAR.name, (60.0, 60.0, 800.0), Spectrum(4000, 4000), "z"))]
    fig, lims, dims = build_cut_plate("t", rs[0], params)
    assert all(im.get_clim() == (-20.0, 20.0) for im in lims) and all(im.get_clim() == (-10.0, 10.0) for im in dims)
    assert np.array_equal(np.asarray(lims[1].get_array()), rs[0].planes["XY"].lev_foc)
    fig2, ims = build_comparison_plate("cmp", rs, params)
    assert len(ims) == 6 and all(im.get_clim() == (-20.0, 20.0) for im in ims)
    assert np.array_equal(np.asarray(ims[3].get_array()), rs[1].planes["XZ"].lev_foc)  # role: beam axis + x for a z beam


def test_cli_field_cuts_and_campaign(tmp_path, capsys):
    rc = main(["field-cuts", "--array", str(PLANAR), "--target", "60", "60", "800", "--fmin", "2000", "--fmax", "4000",
               "--df", "500", "--step", "100", "--volume", "--volume-step", "100", "--outdir", str(tmp_path / "c")])
    assert rc == 0
    out = capsys.readouterr().out
    assert "X/cuts" in out and "width x/y/z" in out
    z = np.load(tmp_path / "c/data.npz")
    assert z["XY_focused_db"].shape == z["XY_energy_focused"].shape and (tmp_path / "c/cuts.png").exists()
    ind = json.loads(str(z["indicators"]))
    assert ind["front_axis"] == "z" and "volume" in ind and ind["volume"]["step_mm"] == 100.0
    rc = main(["field3d-campaign", "--series", "D", "--step", "100", "--df", "1000", "--no-individual",
               "--array-dir", str(ARRAYS), "--outdir", str(tmp_path / "k")])
    assert rc == 0
    assert (tmp_path / "k/compare-extreme-bands.png").exists() and (tmp_path / "k/cut-indicators.csv").exists()
    rows = json.loads((tmp_path / "k/cut-indicators.json").read_text())
    assert [r["scenario"] for r in rows] == ["D/b1535", "D/b24", "D/b18"]


@pytest.mark.parametrize("extra", [["--step", "30"], ["--step", "0"], ["--fmin", "5000", "--fmax", "1000"], ["--gain", "0"],
                                   ["--range-db", "0"], ["--volume-step", "-1"], ["--half-width", "-5"]])
def test_cli_field_cuts_rejects_invalid_arguments(tmp_path, extra):
    args = ["field-cuts", "--array", str(LINE8), "--target", "350", "800", "0", "--fmin", "2000", "--fmax", "4000",
            "--df", "1000", "--step", "100", "--outdir", str(tmp_path / "x")]
    with pytest.raises(SystemExit) as exc:
        main(args + extra)
    assert exc.value.code == 2


def test_cli_field_cuts_rejects_focus_on_a_speaker(tmp_path):
    with pytest.raises(SystemExit):
        main(["field-cuts", "--array", str(LINE8), "--target", "0", "0", "0", "--fmin", "2000", "--fmax", "2000",
              "--step", "100", "--front-range", "0", "2000", "--outdir", str(tmp_path / "x")])


# ---- Mission 06 section 10: no artificial limit of eight elements in the simulator ------------------------------------------------------------

@pytest.mark.parametrize("n,layout", [(16, "line"), (32, "line"), (16, "grid"), (32, "random")])
def test_simulator_supports_arbitrary_element_counts(n, layout):
    rng = np.random.default_rng(n)
    if layout == "line":
        pos = [(0.03 * i, 0.0, 0.0) for i in range(n)]
    elif layout == "grid":
        pos = [(0.05 * (i % 4), 0.05 * (i // 4), 0.0) for i in range(n)]
    else:
        pos = [tuple(rng.uniform(0, 0.4, 2)) + (0.0,) for _ in range(n)]
    arr = ArrayConfig("big", tuple(Element(i + 1, p) for i, p in enumerate(pos)))
    focus = (0.2, 0.2, 0.9)
    coherent = 20 * math.log10(sum(1 / math.dist(focus, e.position) for e in arr.elements))
    geo = FieldGeometry(arr, np.array([focus]), gain=A)
    em = {d.speaker: d.delay for d in focus_delays(arr, focus)}
    e = band_energy(geo, np.array([1000.0, 3000.0, 6000.0]), np.array([0.2, 0.5, 0.3]), em)
    assert effective_level_db(e, A, -60.0)[0] == pytest.approx(coherent, abs=1e-9)
    rec = receive({e_.speaker: np.sin(2 * np.pi * 3000 * np.arange(2000) / 48000) * 0.01 for e_ in arr.elements},
                  arr, (0.1, 0.3, 1.2), 48000, emission_delays=em)
    assert len(rec.paths) == n and np.all(np.isfinite(rec.signal))


def test_energy_share_uses_spreading_compensated_energy():
    # synthetic field: spherical spreading from the origin, 20 dB lower outside a +/-10 degree sector
    u = np.arange(-1000.0, 1000.1, 10.0)
    v = np.arange(100.0, 2000.1, 10.0)
    uu, vv = np.meshgrid(u, v)
    r = np.hypot(uu, vv)
    ang = np.degrees(np.arctan2(uu, vv))
    gain = np.where(np.abs(ang) <= 10, 1.0, 0.01)
    e = 1e6 * (A / r) ** 2 * gain  # +60 dB keeps every level above the display floor
    focus = (0.0, 800.0, 0.0)
    l_f = float(effective_level_db(np.array(1e6 * (A / 800.0) ** 2), A, -60.0))
    out = analyze_cut("XY", u, v, e, e, focus, np.zeros(3), l_f, PARAMS, 10.0)
    roi = (r >= 600) & (r <= 1000)
    sector = roi & (np.abs(ang) <= 10)
    expected = sector.sum() / (sector.sum() + 0.01 * (roi & ~sector).sum())
    assert out["energy_share_in_3_foc"] == pytest.approx(expected, rel=1e-9)  # E_c is uniform inside the sector
    assert out["area_3_cm2"] == pytest.approx(sector.sum() * 100 / 100)
    assert out["contrast_db"] is None or out["contrast_db"] > 15  # the 20 dB step is not counted as a lobe of the focus
