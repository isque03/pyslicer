"""Tests for G-code HTML preview generation."""

from pathlib import Path

from pyslicer.preview import (
    layers_to_toolpaths_3d,
    parse_gcode_layers,
    render_gcode_html,
    write_gcode_preview,
)


SAMPLE = """M109 S180
G90
;; New Layer Z: 1.0
G1 F2400 Z1.0
G1 X0 Y0
G1 X10 Y0 E1.0
G1 X10 Y10 E2.0
G1 X0 Y10
;; New Layer Z: 2.0
G1 F2400 Z2.0
G1 X1 Y1
G1 X2 Y1 E2.5
"""


def test_parse_layers_separates_extrude_and_travel():
    layers = parse_gcode_layers(SAMPLE)
    assert len(layers) == 2
    assert layers[0]["z"] == 1.0
    assert len(layers[0]["extrude"]) >= 1
    assert len(layers[0]["travel"]) >= 1


def test_render_html_is_flat_no_card_chrome():
    doc = render_gcode_html(SAMPLE, title="Test cube")
    assert "Test cube" in doc
    assert "border-radius" not in doc
    assert "box-shadow" not in doc
    assert "path class=\"extrude\"" in doc or 'class="extrude"' in doc


def test_render_html_includes_threejs_orbit_viewer():
    doc = render_gcode_html(SAMPLE, title="3D")
    assert 'id="viewer3d"' in doc
    assert "OrbitControls" in doc
    assert "three@0.170.0" in doc
    assert "window.PYSLICER_TOOLPATHS" in doc


def test_layers_to_toolpaths_3d_uses_layer_z():
    layers = parse_gcode_layers(SAMPLE)
    paths = layers_to_toolpaths_3d(layers)
    assert paths["extrude"]
    # First extrude segment should sit on Z=1.0
    seg = paths["extrude"][0]
    assert seg[2] == 1.0 and seg[5] == 1.0


def test_write_gcode_preview(tmp_path):
    gcode = tmp_path / "out.gcode"
    html = tmp_path / "out.html"
    gcode.write_text(SAMPLE)
    write_gcode_preview(gcode, html, title="Saved")
    text = html.read_text()
    assert text.startswith("<!DOCTYPE html>")
    assert "Saved" in text
    assert Path(html).exists()
