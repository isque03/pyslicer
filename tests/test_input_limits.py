"""Regression coverage for slicer and preview input limits."""

import os
from pathlib import Path
from struct import pack
from types import SimpleNamespace

import pytest

from pyslicer import demo
from pyslicer.cli import main
from pyslicer.config import ConfigError, apply_config_to_model, load_layered_config
from pyslicer.geometry.contour import Contour
from pyslicer.infill import simple_linear_infill
from pyslicer.limits import InputError
from pyslicer.mesh import Model, read_file
from pyslicer.mesh import stl_io
from pyslicer.preview import gcode_parse, parse_gcode, write_gcode_preview
from pyslicer.preview import render_anim
from pyslicer.slicing import process

FIXTURE = Path(__file__).parent / "fixtures" / "cube_10mm.stl"


@pytest.mark.parametrize("option,value", [
    ("-l", "0"), ("-l", "-0.1"), ("-l", "nan"),
    ("--outer-speed", "1e308"), ("--outer-speed", "inf"), ("-n", "65"),
])
def test_bad_cli_overrides_fail_before_creating_output(tmp_path, capsys, option, value):
    output = tmp_path / "output.gcode"
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURE), str(output), option, value])
    assert exc.value.code == 2
    assert "invalid config" in capsys.readouterr().err.lower()
    assert not output.exists()


@pytest.mark.parametrize("config", [
    {"layer_height": float("nan")}, {"processing_threads": 33},
    {"processing_threads": 2.0}, {"min_layer_time": float("inf")},
])
def test_invalid_yaml_values(config):
    with pytest.raises(ConfigError):
        apply_config_to_model(Model(), config)


def test_large_config_file_rejected(tmp_path, monkeypatch):
    monkeypatch.setattr("pyslicer.config.MAX_CONFIG_BYTES", 16)
    path = tmp_path / "big.yaml"
    path.write_text("#" * 17)
    with pytest.raises(ConfigError, match="byte limit"):
        load_layered_config([path])


def test_config_file_cannot_be_used_as_output(tmp_path, capsys):
    path = tmp_path / "profile.yaml"
    path.write_text("layer_height: 0.2\n")
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURE), str(path), "--config", f"{path},{path}"])
    assert exc.value.code == 2
    assert "paths must be distinct" in capsys.readouterr().err
    assert path.read_text() == "layer_height: 0.2\n"


def test_stl_size_and_facet_limits_precede_mesh_allocation(monkeypatch):
    model = Model()
    monkeypatch.setattr(stl_io, "MAX_FACETS", 11)
    with pytest.raises(InputError, match="facet limit"):
        read_file(FIXTURE, model)
    assert not model.facets

    monkeypatch.setattr(stl_io, "MAX_STL_BYTES", FIXTURE.stat().st_size - 1)
    with pytest.raises(InputError, match="byte limit"):
        read_file(FIXTURE, model)


def test_stl_binary_and_ascii_coordinates_checked(tmp_path):
    binary = tmp_path / "binary.stl"
    record = pack("<12fH", 0, 0, 1, float("nan"), 0, 0, 1, 0, 0, 0, 1, 0, 0)
    binary.write_bytes(bytes(80) + pack("<I", 1) + record)
    ascii_file = tmp_path / "ascii.stl"
    ascii_file.write_text("facet normal 0 0 1\nvertex 0 200000 0\n")
    for path in (binary, ascii_file):
        with pytest.raises(InputError, match="coordinates must be finite"):
            read_file(path, Model())


def test_ascii_stl_duplicate_close_does_not_accumulate_edges():
    lines = ["facet normal 0 0 1", "vertex 0 0 0", "vertex 1 0 0", "vertex 0 1 0",
             "endfacet", "endfacet"]
    model = Model()
    with pytest.raises(InputError, match="already closed"):
        stl_io.process_ascii_stl(lines, model)
    assert len(model.facets[0].edges) == 3


def test_direct_layer_and_work_limits(monkeypatch):
    for step in (0, -1, float("inf"), float("nan")):
        with pytest.raises(InputError, match="positive and finite"):
            list(process.layer_range(0, 10, step))
    with pytest.raises(InputError, match="layer limit"):
        list(process.layer_range(0, 10, 0.00001))
    assert list(process.layer_range(0, 10, 3)) == [0, 3, 6, 9, 10]
    model = read_file(FIXTURE, Model())
    monkeypatch.setattr(process, "MAX_FACET_LAYERS", 1)
    with pytest.raises(InputError, match="facet-layer work limit"):
        process.slice_model(model)


def test_infill_checks_spacing_and_work_before_expanding_rows():
    square = Contour.from_path([(0, 0), (1, 0), (1, 1), (0, 1)], 0)
    with pytest.raises(InputError, match="spacing"):
        simple_linear_infill([square], 0, spacing=0)
    with pytest.raises(InputError, match="work limit"):
        simple_linear_infill([square], 0, spacing=0.00001)


@pytest.mark.parametrize("layer_z", ["1", ""])
def test_preview_move_limit(monkeypatch, layer_z):
    monkeypatch.setattr(gcode_parse, "MAX_PREVIEW_MOVES", 1)
    with pytest.raises(InputError, match="move limit"):
        parse_gcode(f";; New Layer Z: {layer_z}\nG1 X0 Y0\nG1 X1 Y0\nG1 X2 Y0\n")


@pytest.mark.parametrize("gcode", [";; New Layer Z: 200000\n", ";; New Layer Z: 1\nG1 X200000\n"])
def test_preview_coordinate_limit(gcode):
    with pytest.raises(InputError, match="coordinates must be finite"):
        parse_gcode(gcode)


def test_preview_feed_limit_prevents_derived_overflow():
    with pytest.raises(InputError, match="feed must be within"):
        parse_gcode(";; New Layer Z: 1\nG1 F" + "9" * 172 + "\n")


def test_preview_cannot_overwrite_input_through_hardlink(tmp_path):
    gcode = tmp_path / "path.gcode"
    output = tmp_path / "output.html"
    gcode.write_text(";; New Layer Z: 1\n")
    os.link(gcode, output)
    with pytest.raises(InputError, match="distinct"):
        write_gcode_preview(gcode, output)
    assert gcode.read_text() == ";; New Layer Z: 1\n"


def test_gif_requires_supported_python_before_downloading(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(render_anim, "sys", SimpleNamespace(version_info=(3, 10)))
    assert demo.main(["--samples-dir", str(tmp_path), "--write-readme-gif"]) == 1
    assert "requires Python 3.11 or newer" in capsys.readouterr().err
    assert not list(tmp_path.iterdir())


@pytest.mark.parametrize("installed", ["12.3.0", "13.0.0rc1"])
def test_gif_rejects_outdated_pillow(monkeypatch, installed):
    monkeypatch.setattr(render_anim, "sys", SimpleNamespace(version_info=(3, 11)))
    monkeypatch.setattr(render_anim, "version", lambda _: installed)
    with pytest.raises(ImportError, match="Pillow 13.0.0 or newer"):
        render_anim.require_gif_support()
