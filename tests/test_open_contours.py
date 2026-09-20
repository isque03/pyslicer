"""process_layer open-contour hard-fail behavior (F12)."""

import pytest

from pyslicer.geometry.facet import Facet
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.mesh.model import Model
from pyslicer.slicing.process import process_layer


def test_open_contours_raise_instead_of_force_close(monkeypatch):
    """Refuse to invent closing segments for open contours."""
    model = Model()
    facet = Facet.withVerticies(
        Vertex(0, 0, 0), Vertex(1, 0, 2), Vertex(0, 1, 2)
    )
    model.facets = [facet]
    model.rebuild_facet_arrays()

    def fake_join(segments):
        from pyslicer.geometry.contour import Contour

        c = Contour()
        c.segments = [
            Line.withVerticies(Vertex(0, 0, 1), Vertex(1, 0, 1)),
            Line.withVerticies(Vertex(1, 0, 1), Vertex(1, 1, 1)),
        ]
        c.closed = False
        segments.clear()
        return c

    monkeypatch.setattr(
        "pyslicer.slicing.process.find_intersecting_lines",
        lambda facet, z: [Line.withVerticies(Vertex(0, 0, 0), Vertex(1, 0, 2))],
    )
    monkeypatch.setattr(
        "pyslicer.slicing.process.get_intersecting_points",
        lambda lines, z: [Vertex(0.5, 0, z), Vertex(0.6, 0.1, z)],
    )
    monkeypatch.setattr("pyslicer.slicing.process.join_segments", fake_join)

    with pytest.raises(ValueError, match="Refusing to invent closing segments"):
        process_layer(1.0, model)
