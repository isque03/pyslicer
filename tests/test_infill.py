"""Infill path generation — known segment hits."""

import math

from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.infill.linear import simple_linear_infill


def _square(z=0.0, size=10.0):
    c = Contour.from_path(
        [(0, 0), (size, 0), (size, size), (0, size), (0, 0)], z
    )
    c.closed = True
    return c


def test_infill_produces_segments_inside_square():
    paths = simple_linear_infill(
        [_square()], 0.0, spacing=2.0, min_extrude=0.5, angle=0.0
    )
    assert len(paths) >= 1
    for p in paths:
        assert p.magnitude() > 0.5
        # Endpoints should lie near the square interior/boundary
        for v in p.verticies:
            assert -0.5 <= v.x <= 10.5
            assert -0.5 <= v.y <= 10.5


def test_infill_empty_contours():
    assert simple_linear_infill([], 0.0) == []
