"""Corner / accel-aware feedrate planning for extrusion paths."""

from __future__ import annotations

import math


def _direction_xy(segment) -> tuple[float, float] | None:
    dx = segment.verticies[1].x - segment.verticies[0].x
    dy = segment.verticies[1].y - segment.verticies[0].y
    mag = math.hypot(dx, dy)
    if mag < 1e-12:
        return None
    return (dx / mag, dy / mag)


def turn_angle(seg_a, seg_b) -> float:
    """Exterior turn angle in radians between two consecutive segments.

    0 = colinear continuation, π/2 = square corner, π = reverse.
    """
    u = _direction_xy(seg_a)
    v = _direction_xy(seg_b)
    if u is None or v is None:
        return 0.0
    dot = max(-1.0, min(1.0, u[0] * v[0] + u[1] * v[1]))
    return math.acos(dot)


def junction_speed(v_square_corner: float, angle_rad: float) -> float:
    """Max speed through a corner (same units as ``v_square_corner``).

    Uses a square-corner-velocity / junction-deviation style scale: a 90°
    corner is limited to ``v_square_corner``; gentler corners allow higher
    speeds; a full reversal forces near-zero.
    """
    if angle_rad < 1e-8:
        return float("inf")
    if angle_rad > math.pi - 1e-8:
        return 0.0
    # GRBL/Klipper convention: junction_cos = -dot(u, v) = -cos(angle)
    junction_cos = -math.cos(angle_rad)
    sin_d2 = math.sqrt(max(0.5 * (1.0 - junction_cos), 0.0))
    if sin_d2 >= 1.0 - 1e-12:
        return float("inf")
    if sin_d2 <= 1e-12:
        return 0.0
    factor = (math.sqrt(2.0) - 1.0) * sin_d2 / (1.0 - sin_d2)
    return v_square_corner * math.sqrt(max(factor, 0.0))


def junction_max_speed_mm_s(
    angle_rad: float,
    *,
    max_corner_speed_mm_s: float,
    max_accel: float,
    max_jerk: float,
    min_angle_rad: float,
) -> float:
    """Strictest junction speed (mm/s) from SCV + accel/jerk thresholds.

    Accel/jerk caps use the same proxies as the concern view:
    ``a = v²(1 - cos θ)``, ``Δv = 2 v sin(θ/2)``. Angles below
    ``min_angle_rad`` skip those caps (SCV still applies).
    """
    if angle_rad < 1e-8:
        return float("inf")
    if angle_rad > math.pi - 1e-8:
        return 0.0

    limits: list[float] = []
    scv = junction_speed(max_corner_speed_mm_s, angle_rad)
    if math.isfinite(scv):
        limits.append(scv)

    if angle_rad >= min_angle_rad:
        one_m_cos = 1.0 - math.cos(angle_rad)
        if one_m_cos > 1e-12 and max_accel > 0:
            limits.append(math.sqrt(max_accel / one_m_cos))
        sin_half = math.sin(angle_rad / 2.0)
        if sin_half > 1e-12 and max_jerk > 0:
            limits.append(max_jerk / (2.0 * sin_half))

    if not limits:
        return float("inf")
    return min(limits)


def plan_contour_feeds(
    segments,
    cruise_f: float,
    max_corner_speed: float,
    max_accel: float,
    *,
    max_jerk: float = 20.0,
    min_corner_angle_deg: float = 20.0,
    closed: bool = False,
) -> list[float]:
    """Return per-segment cruise feeds (mm/min) under corner and accel limits.

    Junction speeds are capped so concern proxies stay within ``max_accel``
    (mm/s²) and ``max_jerk`` (mm/s corner Δv) for turns ≥ ``min_corner_angle_deg``.
    Emitted F is also capped by those junction limits (conservative single-F
    moves) so the G-code itself stays inside the same thresholds the preview
    visualizes.
    """
    n = len(segments)
    if n == 0:
        return []
    if n == 1:
        return [float(cruise_f)]

    cruise_s = max(float(cruise_f), 0.0) / 60.0
    corner_s = max(float(max_corner_speed), 0.0) / 60.0
    a = max(float(max_accel), 1e-6)
    jerk = max(float(max_jerk), 1e-9)
    min_angle = math.radians(max(float(min_corner_angle_deg), 0.0))
    lengths = [max(seg.magnitude(), 1e-9) for seg in segments]

    def _geom_cap(angle: float) -> float:
        return junction_max_speed_mm_s(
            angle,
            max_corner_speed_mm_s=corner_s,
            max_accel=a,
            max_jerk=jerk,
            min_angle_rad=min_angle,
        )

    # Geometric junction caps at the start of each segment (mm/s)
    geom_cap = [cruise_s] * n
    for i in range(1, n):
        jv = _geom_cap(turn_angle(segments[i - 1], segments[i]))
        geom_cap[i] = cruise_s if not math.isfinite(jv) else min(cruise_s, jv)

    if closed and n >= 2:
        jv = _geom_cap(turn_angle(segments[-1], segments[0]))
        geom_cap[0] = cruise_s if not math.isfinite(jv) else min(cruise_s, jv)
    else:
        geom_cap[0] = cruise_s  # open path: no prior corner at start

    # Reachability speeds (may be lower than geom caps after accel pass)
    vertex_v = list(geom_cap)
    if not closed:
        vertex_v[0] = 0.0  # start from rest after travel

    def _backward_pass() -> None:
        for i in range(n - 1, 0, -1):
            max_prev = math.sqrt(vertex_v[i] ** 2 + 2.0 * a * lengths[i - 1])
            if vertex_v[i - 1] > max_prev:
                vertex_v[i - 1] = max_prev
        if closed and n >= 2:
            max_prev = math.sqrt(vertex_v[0] ** 2 + 2.0 * a * lengths[n - 1])
            if vertex_v[n - 1] > max_prev:
                vertex_v[n - 1] = max_prev

    def _forward_pass() -> None:
        for i in range(n - 1):
            max_next = math.sqrt(vertex_v[i] ** 2 + 2.0 * a * lengths[i])
            if vertex_v[i + 1] > max_next:
                vertex_v[i + 1] = max_next
        if closed and n >= 2:
            max_next = math.sqrt(vertex_v[n - 1] ** 2 + 2.0 * a * lengths[n - 1])
            if vertex_v[0] > max_next:
                vertex_v[0] = max_next

    _backward_pass()
    _forward_pass()
    if closed:
        _backward_pass()
        _forward_pass()

    feeds: list[float] = []
    for i in range(n):
        v_start = vertex_v[i]
        if closed:
            v_end = vertex_v[(i + 1) % n]
            g_end = geom_cap[(i + 1) % n]
        elif i + 1 < n:
            v_end = vertex_v[i + 1]
            g_end = geom_cap[i + 1]
        else:
            v_end = 0.0
            g_end = cruise_s

        peak_sq = (2.0 * a * lengths[i] + v_start**2 + v_end**2) / 2.0
        peak = math.sqrt(max(peak_sq, 0.0))
        peak = max(peak, v_start, v_end)

        # Cap emitted F by geometric junction limits so concern(F, θ) stays in budget
        g_start = geom_cap[i] if (closed or i > 0) else cruise_s
        seg_v = min(cruise_s, peak, g_start, g_end)
        feeds.append(seg_v * 60.0)
    return feeds
