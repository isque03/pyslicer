"""STL I/O and Z-range correctness."""

import math
from pathlib import Path

from pyslicer.mesh.model import Model
from pyslicer.mesh.stl_io import read_file, update_model_z_range, process_ascii_stl


FIXTURE = Path(__file__).parent / "fixtures" / "cube_10mm.stl"


def test_update_z_range_descending_first_vertices():
    """F3 regression: first Z values descending must still capture zmax."""
    m = Model()
    for z in [10.0, 5.0, 0.0, 0.0, 0.0]:
        update_model_z_range(z, m)
    assert math.isclose(m.zmin, 0.0)
    assert math.isclose(m.zmax, 10.0)


def test_update_z_range_ascending():
    m = Model()
    for z in [0.0, 5.0, 10.0]:
        update_model_z_range(z, m)
    assert math.isclose(m.zmin, 0.0)
    assert math.isclose(m.zmax, 10.0)


def test_binary_cube_stl(tmp_path=None):
    m = Model()
    read_file(str(FIXTURE), m)
    assert len(m.facets) == 12
    assert m.facet_vertices is not None
    assert m.facet_vertices.shape == (12, 3, 3)
    assert math.isclose(m.zmin, 0.0)
    assert math.isclose(m.zmax, 10.0)


def test_ascii_stl_parsing():
    ascii_stl = """solid cube
  facet normal 0 0 -1
    outer loop
      vertex 0 0 0
      vertex 1 0 0
      vertex 1 1 0
    endloop
  endfacet
  facet normal 0 0 1
    outer loop
      vertex 0 0 2
      vertex 1 1 2
      vertex 1 0 2
    endloop
  endfacet
endsolid cube
"""
    m = Model()
    process_ascii_stl(ascii_stl.splitlines(keepends=True), m)
    m.rebuild_facet_arrays()
    assert len(m.facets) == 2
    assert m.triangles == 2
    assert m.endLoop == 2
    assert math.isclose(m.zmin, 0.0)
    assert math.isclose(m.zmax, 2.0)
