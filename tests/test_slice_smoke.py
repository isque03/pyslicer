"""End-to-end slice smoke test on a tiny cube STL."""

from pathlib import Path

from pyslicer.gcode.parse import parse_gcode
from pyslicer.mesh import Model, read_file
from pyslicer.slicing import slice_model
from pyslicer.infill import simple_linear_infill


FIXTURE = Path(__file__).parent / "fixtures" / "cube_10mm.stl"


def test_slice_cube_produces_layers_and_gcode(tmp_path):
    model = Model()
    model.layerHeight = 2.0
    model.number_perimeters = 1
    model.filament_diameter = 1.75
    read_file(str(FIXTURE), model)

    assert model.facet_vertices is not None
    assert len(model.facets) == 12
    assert model.zmax > model.zmin

    layers = slice_model(model)
    assert len(layers) >= 2

    # At least one layer should have perimeters
    layers_with = [ly for ly in layers if ly.perimeters and ly.perimeters[0]]
    assert len(layers_with) >= 1

    for ly in layers:
        model.infillAtLayer[ly.z] = []
        if ly.perimeters and ly.perimeters[-1]:
            ly.infill = simple_linear_infill(
                ly.perimeters[-1], ly.z, spacing=2.0, min_extrude=0.5, angle=45.0
            )
        else:
            ly.infill = []
        model.layers.append(ly)
        model.infillAtLayer[ly.z] = ly.infill

    out = tmp_path / "cube.gcode"
    model.writeGCode(str(out))
    text = out.read_text()
    commands = parse_gcode(text)
    assert any(c == "G1" and "E" in p for c, p in commands)
    assert ";; New Layer Z:" in text
