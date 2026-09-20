"""Direct unit tests for join_segments / find_closest_segment."""

from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.slicing.segments import find_closest_segment, join_segments


def _seg(x0, y0, x1, y1, z=0.0):
    return Line.withVerticies(Vertex(x0, y0, z), Vertex(x1, y1, z))


def test_join_closed_square_consumes_all_and_closes():
    segments = [
        _seg(0, 0, 1, 0),
        _seg(1, 0, 1, 1),
        _seg(1, 1, 0, 1),
        _seg(0, 1, 0, 0),
    ]
    contour = join_segments(segments)
    assert contour is not None
    assert contour.is_closed()
    assert len(contour.segments) == 4
    assert segments == []


def test_join_reverses_when_closest_is_tail():
    # Start with first segment; next only matches if reversed
    segments = [
        _seg(0, 0, 1, 0),
        _seg(2, 0, 1, 0),  # tail at (1,0) — must reverse to continue
        _seg(2, 0, 2, 1),
        _seg(2, 1, 0, 1),
        _seg(0, 1, 0, 0),
    ]
    contour = join_segments(segments)
    assert contour is not None
    assert len(contour.segments) >= 2
    # After join, chain should be tip-connected
    for i in range(len(contour.segments) - 1):
        a = contour.segments[i].verticies[1]
        b = contour.segments[i + 1].verticies[0]
        assert abs(a.x - b.x) < 1e-6
        assert abs(a.y - b.y) < 1e-6


def test_join_open_gap_leaves_remaining_and_not_closed():
    segments = [
        _seg(0, 0, 1, 0),
        _seg(1, 0, 2, 0),
        _seg(10, 10, 11, 10),  # far away — not stitched into first contour
    ]
    contour = join_segments(segments)
    assert contour is not None
    assert not contour.is_closed()
    assert len(contour.segments) == 2
    assert len(segments) == 1


def test_find_closest_prefers_head_within_eps():
    tip = _seg(0, 0, 1, 0)
    candidates = [_seg(1.01, 0, 2, 0), _seg(5, 5, 6, 6)]
    closest = find_closest_segment(tip, candidates)
    assert closest is candidates[0]
