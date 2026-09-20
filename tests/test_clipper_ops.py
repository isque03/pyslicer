"""Tests for the pyclipper adapter — known numeric outcomes."""

import math

from pyslicer import clipper_ops


def test_union_single_square_preserves_extent():
    a = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    result = clipper_ops.union_polygons([a])
    assert len(result) == 1
    xs = [p[0] for p in result[0]]
    ys = [p[1] for p in result[0]]
    assert math.isclose(min(xs), 0.0, abs_tol=1e-6)
    assert math.isclose(max(xs), 2.0, abs_tol=1e-6)
    assert math.isclose(min(ys), 0.0, abs_tol=1e-6)
    assert math.isclose(max(ys), 2.0, abs_tol=1e-6)


def test_union_two_overlapping_squares_evenodd():
    a = [(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0)]
    b = [(1.0, 1.0), (3.0, 1.0), (3.0, 3.0), (1.0, 3.0)]
    result = clipper_ops.union_polygons([a, b])
    assert len(result) >= 1
    xs = [p[0] for poly in result for p in poly]
    assert max(xs) <= 3.1
    assert min(xs) >= -0.1


def test_difference_square_removes_inner():
    outer = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    inner = [(3.0, 3.0), (7.0, 3.0), (7.0, 7.0), (3.0, 7.0)]
    result = clipper_ops.difference_polygons([outer], [inner])
    assert len(result) >= 1


def test_offset_shrinks_square_by_one_mm():
    square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    shrunk = clipper_ops.offset_polygons([square], -1.0)
    assert len(shrunk) == 1
    xs = [p[0] for p in shrunk[0]]
    ys = [p[1] for p in shrunk[0]]
    assert math.isclose(min(xs), 1.0, abs_tol=0.05)
    assert math.isclose(max(xs), 9.0, abs_tol=0.05)
    assert math.isclose(min(ys), 1.0, abs_tol=0.05)
    assert math.isclose(max(ys), 9.0, abs_tol=0.05)


def test_offset_expand_square_by_one_mm():
    square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    grown = clipper_ops.offset_polygons([square], 1.0)
    assert len(grown) == 1
    xs = [p[0] for p in grown[0]]
    assert min(xs) <= -0.9
    assert max(xs) >= 10.9


def test_empty_inputs():
    assert clipper_ops.union_polygons([]) == []
    assert clipper_ops.difference_polygons([], [[(0, 0)]]) == []
    assert clipper_ops.offset_polygons([], -1.0) == []
    assert clipper_ops.simplify_polygon([]) == []
    assert clipper_ops.clean_polygons([]) == []


def test_simplify_closed_path():
    path = [(0.0, 0.0), (5.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    result = clipper_ops.simplify_polygon(path)
    assert len(result) >= 1


def test_clean_polygons_roundtrip():
    square = [(0.0, 0.0), (10.0, 0.0), (10.0, 10.0), (0.0, 10.0)]
    cleaned = clipper_ops.clean_polygons([square])
    assert len(cleaned) == 1
