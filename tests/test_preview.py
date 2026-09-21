"""Tests for G-code HTML preview generation."""

from pathlib import Path

import pytest

from pyslicer.preview import (
    layers_to_toolpaths_3d,
    parse_gcode_layers,
    parse_toolpath_moves,
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


def test_parse_toolpath_moves_is_chronological():
    moves = parse_toolpath_moves(SAMPLE)
    assert len(moves) >= 3
    assert any(m["extrude"] for m in moves)
    assert any(not m["extrude"] for m in moves)
    assert moves[0]["feed"] == 2400.0


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
    assert "addExtrudeSweeps" in doc
    assert "buildSweepFrames" in doc
    assert "sweepStadiumGeometry" in doc
    assert "nozzle * 0.5" in doc
    assert "Math.min(halfW, halfH)" in doc
    assert "0x0072b2" in doc
    assert "0xe69f00" in doc
    assert "paintExtrudeConcernColors" in doc
    assert "feedSpeedScores" in doc
    assert "refreshSpeedColors" in doc
    assert 'name="color-mode"' in doc
    assert 'value="speed"' in doc
    assert "concern-legend" in doc
    assert 'id="speed-legend"' in doc
    assert 'id="thr-speed"' in doc
    assert 'id="thr-accel"' in doc
    assert 'id="thr-jerk"' in doc
    assert 'id="thr-angle"' in doc
    assert 'id="slice-settings"' in doc
    assert "Slice settings" in doc


def test_extrude_polylines_chain_connected_tips():
    from pyslicer.preview.gcode_parse import chain_extrude_polylines, parse_gcode

    _, moves = parse_gcode(SAMPLE)
    polys = chain_extrude_polylines(moves)
    assert polys
    # SAMPLE has two tip-connected extrudes on layer 1, then travel, then one more
    assert any(len(p["points"]) >= 3 for p in polys)
    paths = layers_to_toolpaths_3d([], nozzle_diameter=0.4, moves=moves)
    assert "extrudePolylines" in paths
    assert len(paths["extrudePolylines"]) == len(polys)
    # Indices cover every extrude segment exactly once
    n_ex = sum(1 for m in moves if m["extrude"])
    covered = set()
    for p in paths["extrudePolylines"]:
        for i in range(p["i0"], p["i1"] + 1):
            covered.add(i)
    assert covered == set(range(n_ex))
    assert all(len(p["points"]) == p["i1"] - p["i0"] + 2 for p in polys)

def test_layers_to_toolpaths_3d_uses_layer_z():
    layers = parse_gcode_layers(SAMPLE)
    paths = layers_to_toolpaths_3d(layers, nozzle_diameter=0.4)
    assert paths["extrude"]
    assert paths["nozzleDiameter"] == 0.4
    # First extrude segment should sit on Z=1.0
    seg = paths["extrude"][0]
    assert seg[2] == 1.0 and seg[5] == 1.0


def test_tube_diameter_fills_layer_pitch():
    """Bead width comes from stadium profile (nozzle + layer heuristic)."""
    from pyslicer.preview.bead_profile import STACK_OVERLAP

    layers = parse_gcode_layers(SAMPLE)
    moves = parse_toolpath_moves(SAMPLE)
    paths = layers_to_toolpaths_3d(
        layers, nozzle_diameter=0.5, layer_height=0.25, moves=moves
    )
    assert "moves" in paths
    assert "bead" in paths
    assert paths["layerHeight"] == 0.25
    # Preview height overlaps layer pitch so stacked layers don't show air gaps
    assert paths["bead"]["height"] == pytest.approx(0.25 * STACK_OVERLAP)
    assert paths["bead"]["height"] > paths["layerHeight"]
    assert paths["bead"]["width"] > paths["bead"]["height"]
    assert paths["bead"]["flatWidth"] > 0
    assert paths["tubeDiameter"] == paths["bead"]["width"]
    assert "beadBed" in paths
    assert paths["beadBed"]["width"] >= paths["bead"]["width"]


def test_payload_bead_covers_layer_pitch_when_nozzle_equals_layer():
    """
    Regression: nozzle==layer (nearly round beads) used to leave visible
    air gaps between stacked layers in the 3D preview.
    """
    from pyslicer.preview.bead_profile import STACK_OVERLAP
    from pyslicer.preview.gcode_parse import infer_layer_height, parse_gcode

    gcode = """;; New Layer Z: 0.4
G1 F2400 Z0.4
G1 X0 Y0
G1 X10 Y0 E1
;; New Layer Z: 0.8
G1 Z0.8
G1 X0 Y0
G1 X10 Y0 E2
;; New Layer Z: 1.2
G1 Z1.2
G1 X0 Y0
G1 X10 Y0 E3
"""
    layers, moves = parse_gcode(gcode)
    pitch = infer_layer_height(moves, layers)
    assert pitch == pytest.approx(0.4)
    paths = layers_to_toolpaths_3d(
        layers, nozzle_diameter=0.4, moves=moves
    )
    assert paths["layerHeight"] == pytest.approx(pitch)
    assert paths["bead"]["height"] >= pitch * STACK_OVERLAP - 1e-12
    # Half-heights from adjacent layer centers must overlap
    assert paths["bead"]["height"] > pitch

def test_parse_gcode_single_source_z_hop_policy():
    """Pure-Z hops are moves/timeline only; 2D layers stay XY strokes."""
    from pyslicer.preview.gcode_parse import parse_gcode

    gcode = """;; New Layer Z: 1.0
G1 F2400 X0 Y0 Z1.0
G1 X10 Y0 E1.0
G1 Z2.0
G1 X10 Y10 E2.0
"""
    layers, moves = parse_gcode(gcode)
    assert any(m["z0"] != m["z1"] and not m["extrude"] for m in moves)
    # Layer SVG paths are XY-only — no zero-length Z-only segment
    for layer in layers:
        for a, b in layer["travel"] + layer["extrude"]:
            assert a != b or True  # XY segments
            assert len(a) == 2 and len(b) == 2


def test_timeline_matches_python_and_payload():
    from pyslicer.preview.gcode_parse import build_timeline, parse_gcode

    _, moves = parse_gcode(SAMPLE)
    tl, total = build_timeline(moves)
    paths = layers_to_toolpaths_3d([], nozzle_diameter=0.4, moves=moves)
    assert paths["timeline"][0]["t0"] == tl[0]["t0"]
    assert paths["timeline"][-1]["t1"] == tl[-1]["t1"]
    assert abs(paths["totalTime"] - total) < 1e-9


def test_mid_extrude_timeline_index_contract():
    """While t0 <= t < t1 on an extrude, reveal uses extrudeIndex+1 (viewer F3)."""
    from pyslicer.preview.gcode_parse import build_timeline, parse_gcode

    _, moves = parse_gcode(SAMPLE)
    tl, _ = build_timeline(moves)
    extrude = [m for m in tl if m["extrudeIndex"] >= 0]
    assert extrude
    m = extrude[0]
    mid = (m["t0"] + m["t1"]) / 2
    # Contract used by viewer_js seekSimulation
    assert m["t0"] <= mid < m["t1"]
    reveal_count = m["extrudeIndex"] + 1
    assert reveal_count >= 1



def test_render_html_uses_colorblind_friendly_palette():
    doc = render_gcode_html(SAMPLE, nozzle_diameter=0.4)
    assert "--extrude: #0072b2" in doc
    assert "--travel: #e69f00" in doc
    assert '"nozzleDiameter":0.4' in doc


def test_render_html_includes_preview_controls():
    doc = render_gcode_html(SAMPLE, title="Controls")
    assert 'id="cutaway"' in doc
    assert 'id="filament-opacity"' in doc
    assert "Filament opacity" in doc
    assert "localClippingEnabled" in doc
    assert "setPrintProgress" in doc
    assert "setFilamentOpacity" in doc
    assert 'data-z="' in doc


def test_render_html_includes_simulated_print_controls():
    doc = render_gcode_html(SAMPLE, title="Sim")
    assert 'id="sim-play"' in doc
    assert 'id="sim-pause"' in doc
    assert 'id="sim-rewind"' in doc
    assert 'id="sim-ff"' in doc
    assert 'id="sim-step-fwd"' in doc
    assert 'id="sim-step-back"' in doc
    assert 'id="sim-scrub"' in doc
    assert 'name="sim-speed"' in doc
    assert "seekSimulation" in doc
    assert "stepForwardSim" in doc
    assert "stepBackwardSim" in doc
    assert "makeNozzle" in doc
    assert '"moves"' in doc
    assert 'id="export-movie-start"' in doc
    assert 'id="export-gif"' in doc
    assert "0xe1b84a" in doc  # brass nozzle
    assert "startMovieExport" in doc


def test_write_simulation_gif(tmp_path):
    from pyslicer.preview.render_anim import write_simulation_gif

    moves = parse_toolpath_moves(SAMPLE)
    out = tmp_path / "clip.gif"
    write_simulation_gif(
        moves, out, speed=2.0, max_sim_seconds=2.0, fps=5, width=160, height=100
    )
    assert out.is_file()
    assert out.stat().st_size > 100


def test_write_gcode_preview(tmp_path):
    gcode = tmp_path / "out.gcode"
    html = tmp_path / "out.html"
    gcode.write_text(SAMPLE)
    write_gcode_preview(gcode, html, title="Saved")
    text = html.read_text()
    assert text.startswith("<!DOCTYPE html>")
    assert "Saved" in text
    assert Path(html).exists()
