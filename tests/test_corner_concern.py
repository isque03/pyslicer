"""Tests for preview corner-concern scoring."""

from pyslicer.preview.gcode_parse import (
    corner_concern_metrics,
    corner_concern_scores,
    layers_to_toolpaths_3d,
    scores_from_concern_metrics,
)


def _move(x0, y0, x1, y1, *, feed=4200, extrude=True, z=0.2):
    return {
        "x0": x0,
        "y0": y0,
        "z0": z,
        "x1": x1,
        "y1": y1,
        "z1": z,
        "extrude": extrude,
        "feed": feed,
    }


def test_concern_straight_low():
    moves = [
        _move(0, 0, 10, 0, feed=6000),
        _move(10, 0, 20, 0, feed=6000),
    ]
    scores = corner_concern_scores(moves, max_accel=1.0, min_angle_deg=5)
    assert scores[0] < 0.05
    assert scores[1] < 0.05


def test_concern_right_angle_high_vs_straight():
    corner = [
        _move(0, 0, 10, 0, feed=6000),
        _move(10, 0, 10, 10, feed=6000),
    ]
    straight = [
        _move(0, 0, 10, 0, feed=6000),
        _move(10, 0, 20, 0, feed=6000),
    ]
    # Low accel budget → 90° at 100 mm/s exceeds; straight does not
    c = corner_concern_scores(
        corner, max_speed=200, max_accel=500.0, max_jerk=200.0, min_angle_deg=5
    )
    s = corner_concern_scores(
        straight, max_speed=200, max_accel=500.0, max_jerk=200.0, min_angle_deg=5
    )
    assert max(c) > max(s)
    assert max(c) > 0.5


def test_concern_at_limit_is_cool():
    """Feeds that exactly meet accel/jerk budgets should not light up."""
    # 90°: a = v² → v = sqrt(1000) ≈ 31.62 mm/s → F = 1897
    # jerk = v*√2 ≈ 44.7; set max_jerk high so accel binds
    v = 1000.0**0.5
    feed = v * 60.0
    moves = [
        _move(0, 0, 10, 0, feed=feed),
        _move(10, 0, 10, 10, feed=feed),
    ]
    scores = scores_from_concern_metrics(
        corner_concern_metrics(moves),
        max_speed=200,
        max_accel=1000,
        max_jerk=100,
        min_angle_deg=5,
    )
    assert max(scores) < 0.05


def test_concern_thresholds_gate_by_angle():
    moves = [
        _move(0, 0, 10, 0, feed=6000),
        _move(10, 0, 10, 10, feed=6000),
    ]
    metrics = corner_concern_metrics(moves)
    hot = scores_from_concern_metrics(
        metrics, max_speed=70, max_accel=500, max_jerk=5, min_angle_deg=10
    )
    cold = scores_from_concern_metrics(
        metrics, max_speed=70, max_accel=500, max_jerk=5, min_angle_deg=100
    )
    assert max(hot) > 0.5
    assert max(cold) == 0.0


def test_concern_in_toolpath_payload():
    moves = [
        _move(0, 0, 1, 0, feed=6000),
        _move(1, 0, 1, 1, feed=6000),
    ]
    layers = [{"z": 0.2, "extrude": [((0, 0), (1, 0))], "travel": []}]
    payload = layers_to_toolpaths_3d(
        layers,
        moves=moves,
        max_accel=100.0,
        planning_limits={
            "maxSpeed": 70,
            "outerSpeed": 30,
            "innerSpeed": 60,
            "infillSpeed": 80,
            "maxAccel": 100,
            "maxJerk": 5,
            "minAngleDeg": 5,
            "maxCornerSpeed": 5,
        },
    )
    assert "concern" in payload
    assert "concernMetrics" in payload
    assert "planningLimits" in payload
    assert payload["planningLimits"]["outerSpeed"] == 30
    assert len(payload["concern"]) == len(payload["extrude"])
    assert max(payload["concern"]) > 0


def test_parse_planning_limits_from_gcode_comment():
    from pyslicer.preview.gcode_parse import parse_planning_limits_from_gcode

    gcode = (
        "; pyslicer planning: outer=30.0mm/s inner=60.0mm/s infill=80.0mm/s "
        "accel=1000mm/s^2 jerk=20.0mm/s corner=5.0mm/s min_angle=20deg\n"
        "G1 X0 Y0\n"
    )
    lim = parse_planning_limits_from_gcode(gcode)
    assert lim is not None
    assert lim["outerSpeed"] == 30.0
    assert lim["innerSpeed"] == 60.0
    assert lim["infillSpeed"] == 80.0
    assert lim["maxAccel"] == 1000.0
    assert lim["maxJerk"] == 20.0
    assert lim["maxCornerSpeed"] == 5.0
    assert lim["minAngleDeg"] == 20.0
