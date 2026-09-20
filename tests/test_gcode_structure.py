"""Full GCode structure / correctness tests."""

import math
from pathlib import Path

import pytest

from pyslicer.gcode.parse import parse_gcode, parse_line
from pyslicer.gcode.writer import write_gcode
from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.mesh.layer import Layer
from pyslicer.mesh.model import Model


def _square_contour(z=0.2, size=10.0):
    c = Contour()
    c.zlevel = z
    pts = [(0, 0), (size, 0), (size, size), (0, size), (0, 0)]
    for i in range(len(pts) - 1):
        c.segments.append(
            Line.withVerticies(
                Vertex(pts[i][0], pts[i][1], z),
                Vertex(pts[i + 1][0], pts[i + 1][1], z),
            )
        )
    c.closed = True
    return c


def _tiny_model(with_infill=True, retract_amount=1.0, min_retract=3.0):
    m = Model()
    m.print_temperature = 200.0
    m.retract_amount = retract_amount
    m.layerHeight = 0.2
    m.nozzle_diameter = 0.4
    m.filament_diameter = 1.75
    m.minimum_retract_travel = min_retract

    layer = Layer()
    layer.z = 0.2
    contour = _square_contour(0.2)
    layer.perimeters = [[contour]]
    layer.contours = []
    layer.infill = []
    m.layers = [layer]

    if with_infill:
        s1 = Line.withVerticies(Vertex(1, 1, 0.2), Vertex(2, 1, 0.2))
        s2 = Line.withVerticies(Vertex(2, 1, 0.2), Vertex(3, 1, 0.2))
        s3 = Line.withVerticies(Vertex(8, 5, 0.2), Vertex(9, 5, 0.2))
        m.infillAtLayer[0.2] = [s1, s2, s3]
        layer.infill = m.infillAtLayer[0.2]
    else:
        m.infillAtLayer[0.2] = []
    return m


def test_parse_line_g1():
    cmd, params = parse_line("G1 X1.5 Y2.0 E3.25")
    assert cmd == "G1"
    assert params["X"] == 1.5
    assert params["E"] == 3.25


def test_parse_ignores_comments():
    assert parse_line(";; New Layer") is None
    assert parse_line("G1 X1 ; comment")[0] == "G1"


def test_preamble_and_epilogue(tmp_path):
    m = _tiny_model()
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert "M109 S200.0" in text
    assert "G90" in text
    assert "G21" in text
    assert "M82" in text
    assert "G1 F200 E8" in text
    assert "M104 S0 ; Heat off 200.0C" in text
    assert "G28 X0 Y0" in text
    assert "M84" in text


def test_layer_structure(tmp_path):
    m = _tiny_model()
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert ";; New Layer Z: 0.2" in text
    assert ";; Perimeter" in text
    assert "Z0.200000" in text
    assert ";; Infill" in text


def test_e_monotonic_except_retracts(tmp_path):
    m = _tiny_model()
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    body = out.read_text().split("M107")[0]
    last_e = None
    for cmd, params in parse_gcode(body):
        if "E" not in params:
            continue
        e = params["E"]
        if last_e is not None and e < last_e:
            assert cmd in ("G1", "G92")
            if cmd == "G1":
                assert "X" not in params and "Y" not in params
        last_e = e


def test_infill_retract_when_travel_at_or_above_threshold(tmp_path):
    """Gap from (3,1) to (8,5) ≈ 6.4 mm ≥ 3.0 → retract comment + E dip."""
    m = _tiny_model(min_retract=3.0, retract_amount=1.0)
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    assert ";; retract" in text
    # Travel distance comment should be ≥ threshold
    assert "infill travel move. distance:" in text


def test_infill_no_retract_when_travel_below_threshold(tmp_path):
    """Raise threshold above the gap so retract must not fire."""
    m = _tiny_model(min_retract=100.0, retract_amount=1.0)
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    # First travel (idx 0) has distance 0 — no retract; gap travel also below 100
    assert ";; retract" not in text


def test_missing_infill_key_raises(tmp_path):
    m = _tiny_model()
    del m.infillAtLayer[0.2]
    with pytest.raises(KeyError, match="no infillAtLayer"):
        write_gcode(m, str(tmp_path / "out.gcode"))


def test_golden_snapshot(tmp_path):
    m = _tiny_model(retract_amount=0.0)
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    text = out.read_text()
    fixture = Path(__file__).parent / "fixtures" / "tiny_square.gcode"
    assert fixture.exists(), "golden fixture must be committed; do not auto-create"
    got = parse_gcode(text)
    expected = parse_gcode(fixture.read_text())
    assert [c[0] for c in got] == [c[0] for c in expected]
    for (gc, gp), (ec, ep) in zip(got, expected):
        assert gc == ec
        for k in ep:
            assert math.isclose(gp[k], ep[k], rel_tol=1e-5, abs_tol=1e-5)


def test_perimeter_e_accumulates_per_segment_length(tmp_path):
    """Four 10 mm sides → total ΔE = 4 * (10*h*d)/(πr²) after unretract."""
    m = _tiny_model(with_infill=True, retract_amount=0.0)
    m.infillAtLayer[0.2] = []
    out = tmp_path / "out.gcode"
    write_gcode(m, str(out))
    body = out.read_text().split(";; Infill")[0]
    cmds = parse_gcode(body)
    e_vals = [p["E"] for c, p in cmds if c == "G1" and "E" in p and "X" in p]
    assert len(e_vals) == 4
    length = 10.0
    delta = (length * m.layerHeight * m.nozzle_diameter) / m.pi_r_squared()
    for i, e in enumerate(e_vals):
        assert math.isclose(e, (i + 1) * delta, abs_tol=1e-5)
