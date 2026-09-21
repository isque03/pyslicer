"""Seam placement, seam gap, and wipe helpers for closed perimeters."""

from __future__ import annotations

import math

from pyslicer.gcode.motion import turn_angle
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex


def _vertex_at(segments, index: int) -> Vertex:
    return segments[index].verticies[0]


def sharpest_corner_index(segments) -> int:
    """Index of the vertex with the largest exterior turn (closed loop)."""
    n = len(segments)
    if n < 3:
        return 0
    best_i = 0
    best_angle = -1.0
    for i in range(n):
        prev = segments[(i - 1) % n]
        cur = segments[i]
        ang = turn_angle(prev, cur)
        if ang > best_angle:
            best_angle = ang
            best_i = i
    return best_i


def nearest_vertex_index(segments, x: float, y: float) -> int:
    best_i = 0
    best_d = float("inf")
    for i, seg in enumerate(segments):
        v = seg.verticies[0]
        d = (v.x - x) ** 2 + (v.y - y) ** 2
        if d < best_d:
            best_d = d
            best_i = i
    return best_i


def rear_vertex_index(segments) -> int:
    """Prefer lowest Y, then lowest X (back of bed in typical orientation)."""
    best_i = 0
    best = (float("inf"), float("inf"))
    for i, seg in enumerate(segments):
        v = seg.verticies[0]
        key = (v.y, v.x)
        if key < best:
            best = key
            best_i = i
    return best_i


def rotate_segments(segments, start_index: int):
    """Rotate segment list so ``start_index`` becomes the first segment."""
    if not segments or start_index <= 0:
        return list(segments)
    n = len(segments)
    start_index %= n
    return list(segments[start_index:]) + list(segments[:start_index])


def apply_seam_gap(segments, gap_mm: float):
    """Shorten the end of a closed loop by ``gap_mm``, walking backward.

    Consumes whole trailing segments when needed. Leaves at least one segment.
    If the requested gap exceeds the path length minus one segment, the gap is
    applied as far as possible (explicit clamp).
    """
    if gap_mm <= 0 or len(segments) < 2:
        return list(segments)
    out = list(segments)
    remaining = float(gap_mm)
    while remaining > 1e-9 and len(out) > 1:
        last = out[-1]
        length = last.magnitude()
        if length <= remaining + 1e-9:
            remaining -= length
            out.pop()
            continue
        t = (length - remaining) / length
        a, b = last.verticies[0], last.verticies[1]
        nx = a.x + (b.x - a.x) * t
        ny = a.y + (b.y - a.y) * t
        nz = a.z
        out[-1] = Line.withVerticies(a, Vertex(nx, ny, nz))
        remaining = 0.0
    return out


def choose_seam_index(
    segments,
    mode: str,
    *,
    last_xy: tuple[float, float] | None = None,
    sticky_xy: tuple[float, float] | None = None,
) -> int:
    mode = (mode or "none").lower()
    if mode == "none" or len(segments) < 2:
        return 0
    if mode == "aligned":
        # Sticky world-space point survives topology changes across islands/layers.
        if sticky_xy is not None:
            return nearest_vertex_index(segments, sticky_xy[0], sticky_xy[1])
        return sharpest_corner_index(segments)
    if mode == "nearest":
        if last_xy is not None:
            return nearest_vertex_index(segments, last_xy[0], last_xy[1])
        return 0
    if mode == "rear":
        return rear_vertex_index(segments)
    return 0


def prepare_contour_segments(
    contour,
    mode: str,
    *,
    last_xy: tuple[float, float] | None = None,
    sticky_xy: tuple[float, float] | None = None,
    seam_gap: float = 0.0,
):
    """Return (segments, seam_start_xy) ready for emission."""
    segs = list(contour.segments)
    if not contour.is_closed() or (mode or "none").lower() == "none":
        start = segs[0].verticies[0] if segs else None
        xy = (start.x, start.y) if start is not None else None
        return segs, xy
    idx = choose_seam_index(
        segs, mode, last_xy=last_xy, sticky_xy=sticky_xy
    )
    segs = rotate_segments(segs, idx)
    segs = apply_seam_gap(segs, seam_gap)
    start = segs[0].verticies[0]
    return segs, (start.x, start.y)


def contour_centroid(segments) -> tuple[float, float]:
    if not segments:
        return (0.0, 0.0)
    sx = sy = 0.0
    for seg in segments:
        sx += seg.verticies[0].x
        sy += seg.verticies[0].y
    n = len(segments)
    return (sx / n, sy / n)


def wipe_targets(
    segments,
    *,
    wipe_distance: float,
    wipe_on_loops: bool,
) -> list[tuple[float, float]]:
    """Points to travel (no extrusion) after finishing a closed loop.

    Order: optional inward tuck, then wipe back along the printed path.
    """
    if not segments:
        return []
    points: list[tuple[float, float]] = []
    end = segments[-1].verticies[1]
    if wipe_on_loops:
        cx, cy = contour_centroid(segments)
        dx, dy = cx - end.x, cy - end.y
        mag = math.hypot(dx, dy)
        tuck = min(0.2, wipe_distance if wipe_distance > 0 else 0.2)
        if mag > 1e-9 and tuck > 0:
            points.append((end.x + dx / mag * tuck, end.y + dy / mag * tuck))
    if wipe_distance > 0:
        remaining = wipe_distance
        # Walk backward along the path from the end.
        for seg in reversed(segments):
            if remaining <= 1e-9:
                break
            a, b = seg.verticies[0], seg.verticies[1]
            length = seg.magnitude()
            if length < 1e-12:
                continue
            use = min(remaining, length)
            t = use / length
            # From b toward a
            wx = b.x + (a.x - b.x) * t
            wy = b.y + (a.y - b.y) * t
            points.append((wx, wy))
            remaining -= use
    return points
