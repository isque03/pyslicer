"""Unit tests for corner / accel feed planning."""

import math

from pyslicer.gcode.motion import junction_speed, plan_contour_feeds, turn_angle
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex


def _seg(x0, y0, x1, y1, z=0.0):
    return Line.withVerticies(Vertex(x0, y0, z), Vertex(x1, y1, z))


def test_turn_angle_straight():
    a = _seg(0, 0, 1, 0)
    b = _seg(1, 0, 2, 0)
    assert turn_angle(a, b) < 1e-6


def test_turn_angle_right_angle():
    a = _seg(0, 0, 1, 0)
    b = _seg(1, 0, 1, 1)
    assert math.isclose(turn_angle(a, b), math.pi / 2, rel_tol=1e-9)


def test_turn_angle_reverse():
    a = _seg(0, 0, 1, 0)
    b = _seg(1, 0, 0, 0)
    assert math.isclose(turn_angle(a, b), math.pi, rel_tol=1e-9)


def test_junction_speed_ninety_degrees():
    v = junction_speed(300.0, math.pi / 2)
    assert math.isclose(v, 300.0, rel_tol=1e-6)


def test_junction_speed_straight_unlimited():
    assert junction_speed(300.0, 0.0) == float("inf")


def test_junction_speed_reverse_zero():
    assert junction_speed(300.0, math.pi) == 0.0


def test_plan_l_shape_slows_at_corner():
    """L-corner: both legs capped by junction limits (conservative single-F)."""
    segs = [
        _seg(0, 0, 20, 0),
        _seg(20, 0, 20, 0.5),  # short leg after 90°
    ]
    feeds = plan_contour_feeds(
        segs, cruise_f=4200, max_corner_speed=300, max_accel=1000, closed=False
    )
    assert len(feeds) == 2
    assert feeds[0] <= 300 + 1e-6
    assert feeds[1] <= 300 + 1e-6
    assert feeds[0] < 4200
    assert feeds[1] < 4200


def test_plan_colinear_keeps_cruise():
    segs = [
        _seg(0, 0, 10, 0),
        _seg(10, 0, 20, 0),
        _seg(20, 0, 30, 0),
    ]
    feeds = plan_contour_feeds(
        segs, cruise_f=4200, max_corner_speed=300, max_accel=1000, closed=False
    )
    # After accelerating from rest, later segments should reach cruise
    assert feeds[-1] >= 4190


def test_plan_closed_square_corners_limit_short_sides():
    """Tiny square: sides too short to reach cruise given corner limit."""
    s = 0.4
    segs = [
        _seg(0, 0, s, 0),
        _seg(s, 0, s, s),
        _seg(s, s, 0, s),
        _seg(0, s, 0, 0),
    ]
    feeds = plan_contour_feeds(
        segs, cruise_f=4200, max_corner_speed=300, max_accel=1000, closed=True
    )
    assert all(f < 2000 for f in feeds)


def test_planned_feeds_keep_concern_within_thresholds():
    """Emitted F at a sharp corner should not exceed accel/jerk budgets."""
    from pyslicer.preview.gcode_parse import (
        corner_concern_metrics,
        scores_from_concern_metrics,
    )

    segs = [
        _seg(0, 0, 15, 0),
        _seg(15, 0, 15, 15),
        _seg(15, 15, 0, 15),
        _seg(0, 15, 0, 0),
    ]
    max_accel = 1000.0
    max_jerk = 20.0
    max_corner = 300.0
    cruise = 4200.0
    feeds = plan_contour_feeds(
        segs,
        cruise_f=cruise,
        max_corner_speed=max_corner,
        max_accel=max_accel,
        max_jerk=max_jerk,
        min_corner_angle_deg=20.0,
        closed=True,
    )
    moves = []
    z = 0.2
    pts = [(0, 0), (15, 0), (15, 15), (0, 15), (0, 0)]
    for i, feed in enumerate(feeds):
        x0, y0 = pts[i]
        x1, y1 = pts[i + 1]
        moves.append(
            {
                "x0": x0,
                "y0": y0,
                "z0": z,
                "x1": x1,
                "y1": y1,
                "z1": z,
                "extrude": True,
                "feed": feed,
            }
        )
    scores = scores_from_concern_metrics(
        corner_concern_metrics(moves),
        max_speed=cruise / 60.0,
        max_accel=max_accel,
        max_jerk=max_jerk,
        min_angle_deg=20.0,
    )
    assert max(scores) < 0.05
    assert max(feeds) <= cruise + 1e-6
    assert max(feeds) <= max_corner + 1e-6
