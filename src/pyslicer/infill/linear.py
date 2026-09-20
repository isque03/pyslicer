"""Linear infill generation (numpy-accelerated ray–segment hits)."""

import logging

import numpy as np

from pyslicer.geometry.intersection import Intersection
from pyslicer.geometry.line import Line
from pyslicer.geometry.math_utils import Math
from pyslicer.geometry.vertex import Vertex

logger = logging.getLogger(__name__)


def bounds(contours, z):
    current_min_x = sys_float_max = float("inf")
    current_max_x = -sys_float_max
    current_min_y = current_min_x
    current_max_y = current_max_x
    for contour in contours:
        for line in contour.segments:
            for vertex in line.verticies:
                current_min_x = min(current_min_x, vertex.x)
                current_min_y = min(current_min_y, vertex.y)
                current_max_x = max(current_max_x, vertex.x)
                current_max_y = max(current_max_y, vertex.y)
    return Line.withVerticies(
        Vertex(current_min_x, current_min_y, z),
        Vertex(current_max_x, current_max_y, z),
    )


def _segment_arrays(contours):
    """Pack all contour segments into (M,2,2) arrays with contour/index maps."""
    segs = []
    meta = []  # (contour, segment_index)
    for contour in contours:
        for idx, segment in enumerate(contour.segments):
            segs.append(
                [
                    [segment.verticies[0].x, segment.verticies[0].y],
                    [segment.verticies[1].x, segment.verticies[1].y],
                ]
            )
            meta.append((contour, idx))
    if not segs:
        return np.empty((0, 2, 2)), []
    return np.asarray(segs, dtype=np.float64), meta


def _ray_segment_u(ray_p0, ray_p1, segs):
    """
    Vectorized Bourke intersection: return ua for each segment, nan if no hit.
    ray_p0/ray_p1: (2,), segs: (M,2,2)
    """
    if len(segs) == 0:
        return np.empty(0)
    x1, y1 = ray_p0
    x2, y2 = ray_p1
    x3 = segs[:, 0, 0]
    y3 = segs[:, 0, 1]
    x4 = segs[:, 1, 0]
    y4 = segs[:, 1, 1]
    denom = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
    numa = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)
    numb = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)
    ua = np.full(len(segs), np.nan)
    ub = np.full(len(segs), np.nan)
    valid = np.abs(denom) > 1e-12
    ua[valid] = numa[valid] / denom[valid]
    ub[valid] = numb[valid] / denom[valid]
    hit = valid & (ua >= 0) & (ua <= 1) & (ub >= 0) & (ub <= 1)
    out = np.full(len(segs), np.nan)
    out[hit] = ua[hit]
    return out


def next_intersection(intersection, intersection_row):
    candidates = []
    for int_in_row in intersection_row:
        if int_in_row.contour is intersection.contour and not int_in_row.used:
            if intersection.intersectionIndex % 2 == intersection_row.index(int_in_row) % 2:
                continue
            candidates.append(int_in_row)
    if len(candidates) == 0:
        return None
    nseg = len(intersection.contour.segments)
    answer = min(
        candidates,
        key=lambda x: min(
            abs(x.segmentIndex - intersection.segmentIndex),
            nseg - abs(x.segmentIndex - intersection.segmentIndex),
        ),
    )
    answer.used = True
    return answer


def simple_linear_infill(
    contours, z, spacing=2.0, min_extrude=1.0, angle=45.0, name="unnamed"
):
    result = []
    if not contours:
        return result
    bounding = bounds(contours, z)
    p0 = bounding.pointOnLine(-1.0)
    p1 = bounding.pointOnLine(2.0)
    x1, y1 = p0.x, p0.y
    x2, y2 = p1.x, p1.y
    segs, meta = _segment_arrays(contours)
    intersection_rows = []
    row_count = 0
    cur_y = float(y1)

    while cur_y <= y2:
        v1 = Vertex(x1, cur_y, z)
        v2 = Vertex(x2, cur_y, z)
        l = Line.withVerticies(v1, v2)
        l.rotate(angle, l.pointOnLine(0.5))
        ray0 = np.array([l.verticies[0].x, l.verticies[0].y])
        ray1 = np.array([l.verticies[1].x, l.verticies[1].y])
        uas = _ray_segment_u(ray0, ray1, segs)
        intersections = []
        for i, ua in enumerate(uas):
            if np.isnan(ua):
                continue
            contour, seg_idx = meta[i]
            intersection = Intersection()
            intersection.uparam = float(ua)
            intersection.row = row_count
            intersection.ray = l.copy()
            intersection.contour = contour
            intersection.segmentIndex = seg_idx
            intersection.z = z
            intersections.append(intersection)

        row_count += 1
        if len(intersections) % 2 == 1:
            logger.error(
                "Intersection at z: %s y:%s resulted in odd number (%s) of intersections. ",
                z,
                cur_y,
                len(intersections),
            )
            Math.remove_duplicate_intersections(intersections)

        sort_reverse = row_count % 2 == 0
        intersections.sort(key=lambda x: x.uparam, reverse=sort_reverse)

        if len(intersections) == 1:
            logger.warning("Only one intersection, skipping it")
            intersection_rows.append([])
        else:
            intersection_rows.append(intersections)
        cur_y += spacing

    keep_going = True
    i = 0
    rows_processed = 0
    while keep_going:
        try:
            start_place = next(obj for obj in intersection_rows[i] if not obj.used)
            start_place_idx = intersection_rows[i].index(start_place)
            direction = -1 if start_place_idx % 2 else 1
            stop = intersection_rows[i][start_place_idx + direction]
            stop.used = True
            start = start_place.point()
            start_place.used = True
            end = stop.point()
            path = Line.withVerticies(start, end)
            if path.magnitude() > min_extrude:
                result.append(path)

            while stop is not None:
                if stop.row + 1 > len(intersection_rows) - 1:
                    next_stop = None
                else:
                    stop.intersectionIndex = intersection_rows[stop.row].index(stop)
                    next_stop = next_intersection(
                        stop, intersection_rows[stop.row + 1]
                    )
                if next_stop is not None:
                    try:
                        next_stop_point = next_stop.point()
                        next_idx = intersection_rows[stop.row + 1].index(next_stop)
                        if next_idx % 2:
                            direction = -1
                        else:
                            direction = 1
                        if next_idx + direction - 1 > len(
                            intersection_rows[stop.row + 1]
                        ):
                            direction = -1
                        if next_idx == 0:
                            direction = 1
                        next_stop_h = intersection_rows[stop.row + 1][
                            next_idx + direction
                        ]
                        next_stop_h.used = True
                        path = Line.withVerticies(next_stop_point, next_stop_h.point())
                        next_stop = next_stop_h
                        if path.magnitude() > min_extrude:
                            result.append(path)
                    except IndexError as exc:
                        logger.error(
                            "nextIntersection: index error %s %s",
                            exc,
                            intersection_rows[stop.row + 1],
                        )
                stop = next_stop
        except StopIteration:
            rows_processed += 1
        except IndexError:
            pass
        i += 1
        if i > len(intersection_rows) - 1:
            i = 0
        if rows_processed == len(intersection_rows):
            keep_going = False
    return result


# Back-compat
simpleLinearInfill = simple_linear_infill
nextIntersection = next_intersection
