"""Tests for print-quality G-code helpers and emission."""

import math

from pyslicer.gcode.cooling import layer_speed_scales, path_time_seconds
from pyslicer.gcode.overhang import (
    apply_overhang_caps,
    overhang_feed_cap_mm_min,
    overlap_percent_at_point,
)
from pyslicer.gcode.pressure import pressure_advance_gcode
from pyslicer.gcode.seam import (
    prepare_contour_segments,
    rear_vertex_index,
    sharpest_corner_index,
)
from pyslicer.gcode.writer import write_gcode
from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.mesh.layer import Layer
from pyslicer.mesh.model import Model


def _square(z=0.2, size=10.0, origin=(0.0, 0.0)):
    ox, oy = origin
    c = Contour()
    c.zlevel = z
    pts = [
        (ox, oy),
        (ox + size, oy),
        (ox + size, oy + size),
        (ox, oy + size),
        (ox, oy),
    ]
    for i in range(len(pts) - 1):
        c.segments.append(
            Line.withVerticies(
                Vertex(pts[i][0], pts[i][1], z),
                Vertex(pts[i + 1][0], pts[i + 1][1], z),
            )
        )
    return c


def test_pressure_advance_dialects():
    assert pressure_advance_gcode("none", 0.05) is None
    assert pressure_advance_gcode("klipper", None) is None
    assert (
        pressure_advance_gcode("klipper", 0.05)
        == "SET_PRESSURE_ADVANCE ADVANCE=0.050000"
    )
    assert pressure_advance_gcode("marlin", 0.18) == "M900 K0.180000"
    assert pressure_advance_gcode("rrf", 0.05) == "M572 D0 S0.050000"


def test_pa_emitted_in_header(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.firmware = "klipper"
    m.pressure_advance = 0.042
    m.seam_position = "none"
    m.outer_perimeter_speed = 3000
    m.inner_perimeter_speed = 4800
    layer = Layer()
    layer.z = 0.2
    layer.perimeters = [[_square()]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = []
    out = tmp_path / "pa.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert "SET_PRESSURE_ADVANCE ADVANCE=0.042000" in text


def test_pa_not_emitted_when_firmware_none(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.firmware = "none"
    m.pressure_advance = 0.05
    m.seam_position = "none"
    layer = Layer()
    layer.z = 0.2
    layer.perimeters = [[_square()]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = []
    out = tmp_path / "nopa.gcode"
    write_gcode(m, str(out))
    assert "SET_PRESSURE_ADVANCE" not in out.read_text()
    assert "M900" not in out.read_text()
    assert "M572" not in out.read_text()


def test_seam_rear_prefers_min_y():
    c = _square()
    idx = rear_vertex_index(c.segments)
    v = c.segments[idx].verticies[0]
    assert math.isclose(v.y, 0.0)


def test_seam_aligned_picks_corner():
    c = _square()
    idx = sharpest_corner_index(c.segments)
    assert 0 <= idx < len(c.segments)


def test_seam_rotation_changes_start(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.retract_amount = 0
    m.seam_position = "rear"
    m.max_corner_speed = 100000
    m.max_accel = 1e9
    m.outer_perimeter_speed = 3000
    layer = Layer()
    layer.z = 0.2
    # Start naturally at (0,0); rear should also be a min-y vertex.
    layer.perimeters = [[_square()]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = []
    segs, _ = prepare_contour_segments(_square(), "rear")
    assert math.isclose(segs[0].verticies[0].y, 0.0)


def test_overlap_inside_is_100():
    poly = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    assert overlap_percent_at_point(5, 5, [poly], 0.5) == 100.0


def test_overlap_far_outside_is_0():
    poly = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    assert overlap_percent_at_point(20, 5, [poly], 0.5) == 0.0


def test_overhang_cap_zero_percent_uses_floor_not_cruise():
    """F1: overhang_speed_0=0 must clamp to floor, not restore full cruise."""
    m = Model()
    m.overhang_speed_0 = 0
    m.slow_down_min_speed = 600
    cap = overhang_feed_cap_mm_min(0.0, 3000.0, m)
    assert math.isclose(cap, 600.0)
    assert cap < 3000.0


def test_seam_gap_consumes_multiple_segments():
    """F3: gap larger than the last edge must walk backward."""
    from pyslicer.gcode.seam import apply_seam_gap

    c = _square(size=10.0)
    # Perimeter length 40; gap 15 should remove one full 10mm edge + 5mm of next.
    out = apply_seam_gap(c.segments, 15.0)
    total = sum(s.magnitude() for s in out)
    assert math.isclose(total, 25.0, abs_tol=1e-6)
    assert len(out) == 3


def test_aligned_seam_sticky_uses_world_xy(tmp_path):
    """F2: aligned seam sticks to world-space point across different island sizes."""
    m = Model()
    m.print_temperature = 200
    m.retract_amount = 0
    m.seam_position = "aligned"
    m.max_corner_speed = 100000
    m.max_accel = 1e9
    m.outer_perimeter_speed = 3000
    # Layer 0: square starting naturally at (0,0); sharpest is first 90° corner.
    layer0 = Layer()
    layer0.z = 0.2
    layer0.perimeters = [[_square(size=10)]]
    # Layer 1: larger square with more vertices — sticky XY should stay near (0,0).
    layer1 = Layer()
    layer1.z = 0.4
    big = _square(z=0.4, size=20)
    layer1.perimeters = [[big]]
    m.layers = [layer0, layer1]
    m.infillAtLayer[0.2] = []
    m.infillAtLayer[0.4] = []
    out = tmp_path / "sticky.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    # First extrude travel of each outer perimeter should start near sticky corner.
    parts = text.split(";; Perimeter 0")
    assert len(parts) >= 3
    for part in parts[1:]:
        travel = [ln for ln in part.splitlines() if ln.startswith("G1 X") and "E" not in ln]
        assert travel, part
        # First XY travel after perimeter header
        line = travel[0]
        assert "X0.000000" in line or "X0.0" in line
        assert "Y0.000000" in line or "Y0.0" in line


def test_wipe_emitted_in_gcode(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.retract_amount = 0
    m.seam_position = "none"
    m.wipe_distance = 2.0
    m.max_corner_speed = 100000
    m.max_accel = 1e9
    layer = Layer()
    layer.z = 0.2
    layer.perimeters = [[_square()]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = []
    out = tmp_path / "wipe.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert ";; wipe" in text


def test_firmware_without_pa_emits_comment(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.firmware = "klipper"
    m.pressure_advance = None
    m.seam_position = "none"
    layer = Layer()
    layer.z = 0.2
    layer.perimeters = [[_square()]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = []
    out = tmp_path / "pa_unset.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert "pressure_advance unset" in text
    assert "SET_PRESSURE_ADVANCE" not in text


def test_print_speed_floor_shared():
    from pyslicer.gcode.cooling import print_speed_floor_mm_min

    m = Model()
    m.slow_down_min_speed = 600
    assert print_speed_floor_mm_min(m) == 600.0
    assert overhang_feed_cap_mm_min(0.0, 3000.0, m) >= 600.0



def test_overhang_slows_unsupported_segment():
    m = Model()
    m.enable_dynamic_overhang_speeds = True
    m.overhang_speed_0 = 20
    m.overhang_speed_25 = 40
    m.overhang_speed_50 = 60
    m.overhang_speed_75 = 80
    m.slow_down_min_speed = 600
    m.nozzle_diameter = 0.4
    m.layerHeight = 0.2
    prev = [(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)]
    # Segment hanging off to the right of previous square
    seg = Line.withVerticies(Vertex(12, 5, 0.2), Vertex(14, 5, 0.2))
    feeds = apply_overhang_caps([seg], [3000.0], [prev], m, 3000.0)
    assert feeds[0] < 3000.0
    assert feeds[0] >= 600.0


def test_cooling_spares_outer_wall():
    m = Model()
    m.min_layer_time = 10.0
    m.dont_slow_down_outer_wall = True
    m.slow_down_min_speed = 600
    m.outer_perimeter_speed = 3000
    m.inner_perimeter_speed = 4800
    m.infill_speed = 4200
    # outer 2s, inner 1s, infill 1s → total 4s < 10s
    o, i, f = layer_speed_scales(m, 2.0, 1.0, 1.0)
    assert math.isclose(o, 1.0)
    assert i < 1.0
    assert f < 1.0


def test_cooling_emits_slower_infill(tmp_path):
    m = Model()
    m.print_temperature = 200
    m.retract_amount = 0
    m.seam_position = "none"
    m.min_layer_time = 3.0
    m.dont_slow_down_outer_wall = True
    m.slow_down_min_speed = 600
    m.outer_perimeter_speed = 3000
    m.inner_perimeter_speed = 4800
    m.infill_speed = 4800
    m.max_corner_speed = 100000
    m.max_accel = 1e9
    layer = Layer()
    layer.z = 0.2
    outer = _square(size=10)
    inner = _square(size=8, origin=(1, 1))
    layer.perimeters = [[outer], [inner]]
    m.layers = [layer]
    m.infillAtLayer[0.2] = [
        Line.withVerticies(Vertex(2, 2, 0.2), Vertex(8, 2, 0.2)),
    ]
    out = tmp_path / "cool.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    peri0 = text.split(";; Perimeter 0")[1].split(";; Infill")[0]
    infill = text.split(";; Infill")[1]
    assert "F3000.000000" in peri0
    # Infill cruise should be scaled below 4800
    assert "F4800.000000" not in infill.split("M106")[0]


def test_path_time_seconds():
    segs = [Line.withVerticies(Vertex(0, 0, 0), Vertex(60, 0, 0))]
    # 60 mm at 60 mm/min → 60 s
    assert math.isclose(path_time_seconds(segs, 60.0), 60.0)
