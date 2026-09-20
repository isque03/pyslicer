"""Segment joining / contour stitching (spatial-hash closest-endpoint)."""

import logging
from collections import defaultdict
from math import sqrt

import numpy as np

from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line
from pyslicer.geometry.math_utils import Math

logger = logging.getLogger(__name__)

_CLOSE_EPS = 0.1  # max squared-distance for a head/tail stitch
_CLOSE_DIST = sqrt(_CLOSE_EPS)
# Cell size so ±2 Moore neighborhood covers the stitch ball (2D; Z is planar).
_GRID = _CLOSE_DIST / 2.0
_INV_GRID = 1.0 / _GRID
_NEIGHBOR_OFFS = tuple(
    (dx, dy) for dx in (-2, -1, 0, 1, 2) for dy in (-2, -1, 0, 1, 2)
)


def _qkey2(x: float, y: float) -> tuple[int, int]:
    return (round(x * _INV_GRID), round(y * _INV_GRID))


def find_closest_segment(segment, segments):
    """Return the open segment whose endpoint is closest to segment's tail.

    Kept for callers/tests; prefer join_segments which batches endpoint arrays.
    """
    if not segments:
        return None
    tip = np.array(
        [segment.verticies[1].x, segment.verticies[1].y, segment.verticies[1].z],
        dtype=np.float64,
    )
    heads = np.empty((len(segments), 3), dtype=np.float64)
    tails = np.empty((len(segments), 3), dtype=np.float64)
    for i, s in enumerate(segments):
        v0, v1 = s.verticies[0], s.verticies[1]
        heads[i] = (v0.x, v0.y, v0.z)
        tails[i] = (v1.x, v1.y, v1.z)
    d_head = np.sum((heads - tip) ** 2, axis=1)
    d_tail = np.sum((tails - tip) ** 2, axis=1)
    i_head = int(np.argmin(d_head))
    i_tail = int(np.argmin(d_tail))
    if d_head[i_head] <= d_tail[i_tail]:
        min_d, closest, reverse = float(d_head[i_head]), i_head, False
    else:
        min_d, closest, reverse = float(d_tail[i_tail]), i_tail, True
    if min_d <= _CLOSE_EPS:
        if reverse:
            segments[closest].reverse()
        return segments[closest]
    return None


def find_connected_segment(segment, segments):
    for x in range(len(segments)):
        if segments[x].verticies[0] == segment.verticies[0]:
            segments[x].reverse()
            return segments[x]
        if segments[x].verticies[1] == segment.verticies[0]:
            return segments[x]
    return None


def _endpoints_close(a, b) -> bool:
    return (
        Math.float_eq(a.x, b.x)
        and Math.float_eq(a.y, b.y)
        and Math.float_eq(a.z, b.z)
    )


def join_segments(segments):
    """Join an unordered set of segments into a closed contour.

    Mutates ``segments`` in place, removing pieces consumed by this contour.
    Uses a 2D spatial hash over endpoints (slice plane is constant-Z) so each
    stitch is O(1) expected rather than a full O(n) nearest-neighbor scan.
    """
    if segments is None or len(segments) == 0:
        return None

    segs = list(segments)
    n = len(segs)
    hx = [0.0] * n
    hy = [0.0] * n
    hz = [0.0] * n
    tx = [0.0] * n
    ty = [0.0] * n
    tz = [0.0] * n
    head_buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
    tail_buckets: dict[tuple[int, int], list[int]] = defaultdict(list)
    for i, s in enumerate(segs):
        v0, v1 = s.verticies[0], s.verticies[1]
        hx[i], hy[i], hz[i] = v0.x, v0.y, v0.z
        tx[i], ty[i], tz[i] = v1.x, v1.y, v1.z
        head_buckets[_qkey2(v0.x, v0.y)].append(i)
        tail_buckets[_qkey2(v1.x, v1.y)].append(i)
    alive = [True] * n

    contour = Contour()
    contour.segments = [segs[0]]
    alive[0] = False
    tipx, tipy, tipz = tx[0], ty[0], tz[0]
    used_tail_keys = {(round(tipx, 7), round(tipy, 7), round(tipz, 7))}

    while True:
        kx, ky = _qkey2(tipx, tipy)
        best_d = _CLOSE_EPS + 1.0
        best_i = -1
        best_reverse = False
        for ox, oy in _NEIGHBOR_OFFS:
            nk = (kx + ox, ky + oy)
            for i in head_buckets.get(nk, ()):
                if not alive[i]:
                    continue
                dx = hx[i] - tipx
                dy = hy[i] - tipy
                dz = hz[i] - tipz
                d = dx * dx + dy * dy + dz * dz
                if d < best_d:
                    best_d = d
                    best_i = i
                    best_reverse = False
            for i in tail_buckets.get(nk, ()):
                if not alive[i]:
                    continue
                dx = tx[i] - tipx
                dy = ty[i] - tipy
                dz = tz[i] - tipz
                d = dx * dx + dy * dy + dz * dz
                if d < best_d:
                    best_d = d
                    best_i = i
                    best_reverse = True
        if best_i < 0 or best_d > _CLOSE_EPS:
            break

        chosen = best_i
        if best_reverse:
            segs[chosen].reverse()
            hx[chosen], tx[chosen] = tx[chosen], hx[chosen]
            hy[chosen], ty[chosen] = ty[chosen], hy[chosen]
            hz[chosen], tz[chosen] = tz[chosen], hz[chosen]

        end_key = (round(tx[chosen], 7), round(ty[chosen], 7), round(tz[chosen], 7))
        if end_key in used_tail_keys:
            logger.warning(
                "already exists a segment(s) with same start point. skipping"
            )
        else:
            contour.segments.append(segs[chosen])
            used_tail_keys.add(end_key)
        alive[chosen] = False
        tipx, tipy, tipz = tx[chosen], ty[chosen], tz[chosen]

    used_ids = {id(segs[i]) for i in range(n) if not alive[i]}
    segments[:] = [s for s in segments if id(s) not in used_ids]

    if _endpoints_close(
        contour.segments[-1].verticies[1], contour.segments[0].verticies[0]
    ):
        contour.closed = True
    if not contour.is_closed():
        gap = Line.withVerticies(
            contour.segments[-1].verticies[1], contour.segments[0].verticies[0]
        )
        logger.error("Created a contour that is not closed. Gap=%s", gap.magnitude())
    return contour


# Back-compat
findClosestSegment = find_closest_segment
findConnectedSegment = find_connected_segment
joinSegments = join_segments
