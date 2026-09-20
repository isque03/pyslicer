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


def test_binary_stl_with_solid_header(tmp_path):
    """Many binary STLs start with b'solid…' and must not be parsed as ASCII."""
    from struct import pack

    path = tmp_path / "solid_header.stl"
    # One triangle, binary layout, header begins with "solid"
    with open(path, "wb") as f:
        f.write(b"solid binary that is not ascii" + b"\0" * 50)
        f.write(pack("<I", 1))
        f.write(pack("<fff", 0, 0, 1))  # normal
        f.write(pack("<fff", 0, 0, 0))
        f.write(pack("<fff", 1, 0, 0))
        f.write(pack("<fff", 0, 1, 1))
        f.write(pack("<H", 0))
    m = Model()
    read_file(str(path), m)
    assert len(m.facets) == 1
    assert math.isclose(m.zmax, 1.0)


def test_padded_binary_stl_loads_facets(tmp_path):
    """Trailing bytes after a valid binary STL must not force empty ASCII parse."""
    from struct import pack

    path = tmp_path / "padded.stl"
    with open(path, "wb") as f:
        f.write(b"solid padded" + b"\0" * 68)
        f.write(pack("<I", 1))
        f.write(pack("<fff", 0, 0, 1))
        f.write(pack("<fff", 0, 0, 0))
        f.write(pack("<fff", 1, 0, 0))
        f.write(pack("<fff", 0, 1, 2))
        f.write(pack("<H", 0))
        f.write(b"\0" * 17)  # trailing slack
    m = Model()
    read_file(str(path), m)
    assert len(m.facets) == 1
    assert math.isclose(m.zmax, 2.0)


def test_empty_stl_raises(tmp_path):
    path = tmp_path / "empty.stl"
    path.write_text("solid empty\nendsolid empty\n", encoding="utf-8")
    m = Model()
    try:
        read_file(str(path), m)
        assert False, "expected ValueError"
    except ValueError as exc:
        assert "no triangles" in str(exc).lower()


def test_ascii_stl_via_process_ascii():
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
