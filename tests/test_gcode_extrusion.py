"""GCode extrusion math — known-answer ΔE and segment sequences."""

import io
import math
import re

from pyslicer.gcode.parse import parse_gcode
from pyslicer.gcode.writer import extrude, retract, write_segment_gcode
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.mesh.model import Model


def _model(**kwargs):
    m = Model()
    m.retract_amount = kwargs.get("retract_amount", 3.0)
    m.retract_speed = 6200
    m.unretract_speed = 3400
    m.default_print_speed = 4200
    m.layerHeight = kwargs.get("layer_height", 0.2)
    m.nozzle_diameter = kwargs.get("nozzle", 0.5)
    m.filament_diameter = kwargs.get("filament", 1.75)
    return m


def _expected_delta_e(length, layer_h, nozzle, filament):
    volume = length * layer_h * nozzle
    area = math.pi * (filament / 2.0) ** 2
    return volume / area


def test_pi_r_squared():
    m = _model(filament=2.0)
    assert math.isclose(m.pi_r_squared(), math.pi)


def test_volume_extruded_unit_length():
    m = _model(layer_height=0.2, nozzle=0.5)
    seg = Line.withVerticies(Vertex(0, 0, 0), Vertex(1, 0, 0))
    assert math.isclose(m.volume_extruded(seg), 0.1)


def test_absolute_e_matches_hand_calculation():
    """ΔE = L · h · d / (π r²) for a 10 mm segment."""
    m = _model(layer_height=0.2, nozzle=0.4, filament=1.75, retract_amount=0.0)
    length = 10.0
    seg = Line.withVerticies(Vertex(0, 0, 0), Vertex(length, 0, 0))
    expected = _expected_delta_e(length, 0.2, 0.4, 1.75)
    assert math.isclose(m.volume_extruded(seg) / m.pi_r_squared(), expected, rel_tol=1e-12)

    buf = io.StringIO()
    e = write_segment_gcode(m, buf, 1, seg, 0.0, expected, feed=4200)[0]
    assert math.isclose(e, expected, abs_tol=1e-9)
    cmds = parse_gcode(buf.getvalue())
    assert len(cmds) == 1
    assert cmds[0][0] == "G1"
    assert math.isclose(cmds[0][1]["E"], expected, abs_tol=1e-6)
    assert math.isclose(cmds[0][1]["X"], length, abs_tol=1e-6)


def test_retract_decreases_e_by_exact_amount():
    m = _model(retract_amount=3.0)
    buf = io.StringIO()
    e = retract(m, 10.0, buf)
    assert math.isclose(e, 7.0)
    cmd, params = parse_gcode(buf.getvalue())[0]
    assert cmd == "G1"
    assert math.isclose(params["E"], 7.0)
    assert "X" not in params


def test_extrude_increases_e_by_exact_amount():
    m = _model(retract_amount=3.0)
    buf = io.StringIO()
    e = extrude(m, 7.0, buf)
    assert math.isclose(e, 10.0)
    assert math.isclose(parse_gcode(buf.getvalue())[0][1]["E"], 10.0)


def test_retract_noop_when_zero():
    m = _model(retract_amount=0.0)
    buf = io.StringIO()
    assert retract(m, 5.0, buf) == 5.0
    assert buf.getvalue() == ""


def test_write_segment_idx0_exact_sequence():
    """idx==0: retract → travel XY → unretract → extrude move with F and E."""
    m = _model(retract_amount=1.0, layer_height=0.2, nozzle=0.5, filament=1.75)
    seg = Line.withVerticies(Vertex(1, 2, 0), Vertex(4, 2, 0))  # length 3
    delta = _expected_delta_e(3.0, 0.2, 0.5, 1.75)
    buf = io.StringIO()
    e0 = 5.0
    e, last_f = write_segment_gcode(m, buf, 0, seg, e0, delta, feed=4200)
    text = buf.getvalue()
    lines = [ln for ln in text.splitlines() if ln.strip()]
    assert lines[0] == "G1 F6200.000000 E4.000000"  # retract 5→4
    assert lines[1] == ";; travel move"
    assert lines[2] == "G1 X1.000000 Y2.000000"
    assert lines[3] == "G1 F3400.000000 E5.000000"  # unretract 4→5
    assert lines[4].startswith("G1 F4200.000000 X4.000000 Y2.000000 E")
    final_e = float(re.search(r"E([0-9.]+)", lines[4]).group(1))
    assert math.isclose(final_e, 5.0 + delta, abs_tol=1e-6)
    assert math.isclose(e, 5.0 + delta, abs_tol=1e-9)
    assert math.isclose(last_f, 4200.0)


def test_write_segment_idx_gt0_extrude_only_no_travel():
    m = _model(retract_amount=1.0)
    buf = io.StringIO()
    seg = Line.withVerticies(Vertex(1, 2, 0), Vertex(4, 2, 0))
    e, _ = write_segment_gcode(m, buf, 1, seg, 2.0, 0.5, feed=4200, last_feed=4200)
    text = buf.getvalue()
    assert "travel" not in text
    assert "retract" not in text.lower() or ";;" not in text
    assert math.isclose(e, 2.5)
    cmds = parse_gcode(text)
    assert len(cmds) == 1
    assert math.isclose(cmds[0][1]["E"], 2.5)
    assert "F" not in cmds[0][1]


def test_write_segment_emits_f_when_feed_changes():
    m = _model(retract_amount=0.0)
    buf = io.StringIO()
    seg = Line.withVerticies(Vertex(0, 0, 0), Vertex(1, 0, 0))
    e, last = write_segment_gcode(m, buf, 1, seg, 0.0, 0.1, feed=1500, last_feed=4200)
    assert math.isclose(last, 1500)
    assert "F1500" in buf.getvalue()
