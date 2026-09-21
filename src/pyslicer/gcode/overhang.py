"""Dynamic overhang speed from previous-layer support overlap."""

from __future__ import annotations

import math

from pyslicer.gcode.cooling import print_speed_floor_mm_min


def _point_in_ring(x: float, y: float, ring) -> bool:
    """Ray-cast even-odd for a single closed ring of (x, y) points."""
    n = len(ring)
    if n < 3:
        return False
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = float(ring[i][0]), float(ring[i][1])
        xj, yj = float(ring[j][0]), float(ring[j][1])
        if ((yi > y) != (yj > y)) and (
            x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi
        ):
            inside = not inside
        j = i
    return inside


def point_in_polygons(x: float, y: float, polygons) -> bool:
    for poly in polygons:
        if _point_in_ring(x, y, poly):
            return True
    return False


def _dist_point_to_segment(px, py, ax, ay, bx, by) -> float:
    dx, dy = bx - ax, by - ay
    len_sq = dx * dx + dy * dy
    if len_sq < 1e-24:
        return math.hypot(px - ax, py - ay)
    t = max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / len_sq))
    qx, qy = ax + t * dx, ay + t * dy
    return math.hypot(px - qx, py - qy)


def min_distance_to_polygons(x: float, y: float, polygons) -> float:
    best = float("inf")
    for poly in polygons:
        n = len(poly)
        if n < 2:
            continue
        for i in range(n):
            a = poly[i]
            b = poly[(i + 1) % n]
            d = _dist_point_to_segment(
                x, y, float(a[0]), float(a[1]), float(b[0]), float(b[1])
            )
            if d < best:
                best = d
    return 0.0 if not math.isfinite(best) else best


def overlap_percent_at_point(
    x: float, y: float, prev_polygons, extrusion_width: float
) -> float:
    """100% = fully supported by previous layer; 0% = full overhang/bridge."""
    width = max(float(extrusion_width), 1e-6)
    if not prev_polygons:
        return 0.0
    if point_in_polygons(x, y, prev_polygons):
        return 100.0
    d = min_distance_to_polygons(x, y, prev_polygons)
    return max(0.0, min(100.0, 100.0 * (1.0 - d / width)))


def segment_overlap_percent(segment, prev_polygons, extrusion_width: float) -> float:
    a, b = segment.verticies[0], segment.verticies[1]
    mx, my = 0.5 * (a.x + b.x), 0.5 * (a.y + b.y)
    return overlap_percent_at_point(mx, my, prev_polygons, extrusion_width)


def overhang_speed_curve(model) -> list[tuple[float, float]]:
    """Control points (overlap_percent, speed_fraction 0–1)."""
    return [
        (0.0, float(model.overhang_speed_0) / 100.0),
        (25.0, float(model.overhang_speed_25) / 100.0),
        (50.0, float(model.overhang_speed_50) / 100.0),
        (75.0, float(model.overhang_speed_75) / 100.0),
        (100.0, 1.0),
    ]


def interpolate_speed_fraction(overlap: float, curve: list[tuple[float, float]]) -> float:
    overlap = max(0.0, min(100.0, float(overlap)))
    for i in range(len(curve) - 1):
        o0, s0 = curve[i]
        o1, s1 = curve[i + 1]
        if overlap <= o1 or i == len(curve) - 2:
            if abs(o1 - o0) < 1e-12:
                return s1
            t = (overlap - o0) / (o1 - o0)
            return s0 + t * (s1 - s0)
    return curve[-1][1]


def overhang_feed_cap_mm_min(
    overlap: float,
    outer_cruise_mm_min: float,
    model,
) -> float:
    """Max feed (mm/min) for a given overlap; never below slow_down_min_speed."""
    frac = interpolate_speed_fraction(overlap, overhang_speed_curve(model))
    capped = float(outer_cruise_mm_min) * max(frac, 0.0)
    floor = print_speed_floor_mm_min(model)
    return max(capped, floor)


def apply_overhang_caps(
    segments,
    feeds_mm_min: list[float],
    prev_polygons,
    model,
    outer_cruise_mm_min: float,
) -> list[float]:
    if not model.enable_dynamic_overhang_speeds:
        return feeds_mm_min
    width = model.offset()
    out = []
    for seg, feed in zip(segments, feeds_mm_min):
        ov = segment_overlap_percent(seg, prev_polygons, width)
        cap = overhang_feed_cap_mm_min(ov, outer_cruise_mm_min, model)
        out.append(min(float(feed), cap))
    return out


def previous_layer_outer_paths(layers, layer_index: int):
    """Outer shell (index 0) paths from the previous layer, as float rings."""
    if layer_index <= 0 or layer_index >= len(layers):
        return []
    prev = layers[layer_index - 1]
    if not prev.perimeters or not prev.perimeters[0]:
        return []
    paths = []
    for contour in prev.perimeters[0]:
        if not contour.is_closed() or len(contour.segments) < 3:
            continue
        paths.append(contour.to_path())
    return paths
