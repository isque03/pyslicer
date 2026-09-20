"""Known-answer geometry and math tests."""

import math

import numpy as np
import pytest

from pyslicer.geometry import (
    Contour,
    Line,
    Math,
    Vertex,
    dot_product,
    multiply_by_scalar_2d,
    subtract_vectors,
)
from pyslicer.geometry.exceptions import (
    CoincidentLines,
    NonIntersectingLines,
    ParallelLines,
)
from pyslicer.infill.linear import _ray_segment_u, bounds
from pyslicer.slicing.plane import get_intersecting_points, slice_at
from pyslicer.mesh.model import Model
from pyslicer.geometry.facet import Facet


def _square(z=0.0, size=10.0, cw=False):
    """Axis-aligned square; cw=False ⇒ counter-clockwise."""
    c = Contour()
    c.zlevel = z
    if cw:
        pts = [(0, 0), (0, size), (size, size), (size, 0), (0, 0)]
    else:
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


def test_empty_contour_not_closed():
    assert not Contour().is_closed()


def test_contour_closedness():
    c = Contour()
    c.segments.append(Line.withVerticies(Vertex(0, 0), Vertex(2, 2)))
    c.segments.append(Line.withVerticies(Vertex(2, 0), Vertex(0, 2)))
    assert not c.is_closed()
    c.segments.append(Line.withVerticies(Vertex(0, 2), Vertex(0, 0)))
    assert c.is_closed()


def test_winding_area_ccw_square():
    c = _square(cw=False)
    assert math.isclose(c.winding_area(), 100.0, abs_tol=1e-9)
    assert c.clockwise() is False


def test_winding_area_cw_square():
    c = _square(cw=True)
    assert math.isclose(c.winding_area(), -100.0, abs_tol=1e-9)
    assert c.clockwise() is True


def test_normal_perpendicular():
    l1 = Line.withVerticies(Vertex(1.0, 0.0, 0.0), Vertex(32.543, 10.4999, 0.0))
    n = l1.normal()
    v1 = subtract_vectors(l1.verticies[1], l1.verticies[0])
    v2 = subtract_vectors(n.verticies[1], n.verticies[0])
    assert Math.float_eq(dot_product(v1, v2), 0.0)


def test_unit_vector_magnitude():
    l1 = Line.withVerticies(Vertex(1.0, 0.0, 0.0), Vertex(32.543, 10.4999, 0.0))
    assert Math.float_eq(l1.unitVector().magnitude(), 1.0)


def test_line_intersection_known_point():
    lint = Line.withVerticies(Vertex(0, 0), Vertex(2, 2))
    lint2 = Line.withVerticies(Vertex(2, 0), Vertex(0, 2))
    vint = lint.intersect2D(lint2)
    assert vint == Vertex(1.0, 1.0, 0.0)


def test_parallel_lines_raise():
    a = Line.withVerticies(Vertex(0, 0), Vertex(1, 0))
    b = Line.withVerticies(Vertex(0, 1), Vertex(1, 1))
    with pytest.raises(ParallelLines):
        a.intersect2DU(b)


def test_non_intersecting_segments_raise():
    a = Line.withVerticies(Vertex(0, 0), Vertex(1, 0))
    b = Line.withVerticies(Vertex(0, 1), Vertex(1, 2))
    with pytest.raises(NonIntersectingLines):
        a.intersect2DU(b)


def test_coincident_lines_raise():
    a = Line.withVerticies(Vertex(0, 0), Vertex(2, 0))
    b = Line.withVerticies(Vertex(0, 0), Vertex(2, 0))
    with pytest.raises(CoincidentLines):
        a.intersect2DUExtend(b)


def test_distance_to_midpoint_is_exact():
    line = Line.withVerticies(Vertex(0, 0, 0), Vertex(2, 0, 0))
    # Point 3 units above midpoint (1,0)
    d = line.distance(Vertex(1.0, 3.0, 0.0))
    assert math.isclose(d, 3.0, abs_tol=1e-9)


def test_distance_outside_segment_returns_minus_one():
    line = Line.withVerticies(Vertex(0, 0, 0), Vertex(2, 0, 0))
    assert line.distance(Vertex(5.0, 1.0, 0.0)) == -1


def test_magnitude_of_3_4_5_triangle():
    line = Line.withVerticies(Vertex(0, 0, 0), Vertex(3, 4, 0))
    assert math.isclose(line.magnitude(), 5.0)
    assert math.isclose(line.magnitudeSquared(), 25.0)


def test_point_on_line_midpoint():
    line = Line.withVerticies(Vertex(0, 0, 0), Vertex(10, 0, 0))
    p = line.pointOnLine(0.5)
    assert math.isclose(p.x, 5.0) and math.isclose(p.y, 0.0)


def test_multiply_by_scalar_2d_preserves_z():
    v = multiply_by_scalar_2d(Vertex(2, 4, 7), 0.5)
    assert math.isclose(v.x, 1.0) and math.isclose(v.y, 2.0) and math.isclose(v.z, 7.0)


def test_rotate_90_about_origin():
    rot = Line.withVerticies(Vertex(0, 0), Vertex(1, 0))
    rot.rotate(90.0, Vertex(0, 0))
    assert Math.float_eq(rot.verticies[1].x, 0.0)
    assert Math.float_eq(rot.verticies[1].y, 1.0)


def test_rotate_90_about_midpoint():
    rot = Line.withVerticies(Vertex(0, 0), Vertex(1, 0))
    rot.rotate(90.0, rot.pointOnLine(0.5))
    assert Math.float_eq(rot.verticies[0].x, 0.5)
    assert Math.float_eq(rot.verticies[1].y, 0.5)


def test_remove_duplicate_intersections():
    class I:
        def __init__(self, u):
            self.uparam = u

    items = [I(0.0), I(0.5), I(0.5 + 1e-9), I(1.0)]
    Math.remove_duplicate_intersections(items)
    assert len(items) == 3


def test_bounding_box_of_square():
    box = _square().bounding_box()
    assert math.isclose(box.verticies[0].x, 0.0)
    assert math.isclose(box.verticies[0].y, 0.0)
    assert math.isclose(box.verticies[1].x, 10.0)
    assert math.isclose(box.verticies[1].y, 10.0)


def test_maybe_intersect_disjoint():
    a = _square()
    b = Contour.from_path([(20, 20), (30, 20), (30, 30), (20, 30)], 0.0)
    assert a.maybe_intersect(b) is False


def test_maybe_intersect_overlapping():
    a = _square()
    b = Contour.from_path([(5, 5), (15, 5), (15, 15), (5, 15)], 0.0)
    assert a.maybe_intersect(b) is True


def test_intersect_brute_force_crossing():
    a = Contour()
    a.segments = [Line.withVerticies(Vertex(0, 1), Vertex(2, 1))]
    b = Contour()
    b.segments = [Line.withVerticies(Vertex(1, 0), Vertex(1, 2))]
    pts = a.intersect_brute_force(b)
    assert len(pts) == 1
    assert pts[0] == Vertex(1.0, 1.0, 0.0)


def test_plane_cut_horizontal_edge():
    """Edge from z=0 to z=2 cut at z=1 → midpoint (1,0,1)."""
    line = Line.withVerticies(Vertex(0, 0, 0), Vertex(2, 0, 2))
    pts = get_intersecting_points([line], 1.0)
    assert len(pts) == 1
    assert math.isclose(pts[0].x, 1.0)
    assert math.isclose(pts[0].y, 0.0)
    assert math.isclose(pts[0].z, 1.0)


def test_ray_segment_u_known_hits():
    # Horizontal ray y=1 from x=-1..3; square bottom/top edges at y=0 and y=2 miss;
    # vertical sides at x=0 and x=2 hit at u corresponding to x=0 and x=2 on ray.
    ray0 = np.array([-1.0, 1.0])
    ray1 = np.array([3.0, 1.0])
    segs = np.array(
        [
            [[0.0, 0.0], [0.0, 2.0]],  # left side
            [[2.0, 0.0], [2.0, 2.0]],  # right side
            [[0.0, 0.0], [2.0, 0.0]],  # bottom (parallel, miss)
        ],
        dtype=np.float64,
    )
    uas = _ray_segment_u(ray0, ray1, segs)
    # Ray length in x is 4; x=0 → u=1/4=0.25; x=2 → u=3/4=0.75
    assert math.isclose(uas[0], 0.25, abs_tol=1e-9)
    assert math.isclose(uas[1], 0.75, abs_tol=1e-9)
    assert np.isnan(uas[2])


def test_slice_at_filters_by_z():
    model = Model()
    low = Facet.withVerticies(Vertex(0, 0, 0), Vertex(1, 0, 0), Vertex(0, 1, 0.5))
    high = Facet.withVerticies(Vertex(0, 0, 5), Vertex(1, 0, 5), Vertex(0, 1, 6))
    model.facets = [low, high]
    model.rebuild_facet_arrays()
    hit = slice_at(model, 0.25)
    assert hit == [low]
    assert slice_at(model, 5.5) == [high]
    assert slice_at(model, 3.0) == []


def test_bounds_of_contours():
    b = bounds([_square(z=1.0)], 1.0)
    assert math.isclose(b.verticies[0].x, 0.0)
    assert math.isclose(b.verticies[1].x, 10.0)
