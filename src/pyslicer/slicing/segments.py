"""Segment joining / contour stitching (numpy closest-endpoint)."""

import logging

import numpy as np

from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line
from pyslicer.geometry.math_utils import Math

logger = logging.getLogger(__name__)

_CLOSE_EPS = 0.1  # max squared-distance for a head/tail stitch


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
    Uses one-shot endpoint arrays + an alive mask so we never rebuild arrays or
    pay for list.remove / Vertex.__eq__ on every stitch.
    """
    if segments is None or len(segments) == 0:
        return None

    segs = list(segments)
    n = len(segs)
    heads = np.empty((n, 3), dtype=np.float64)
    tails = np.empty((n, 3), dtype=np.float64)
    for i, s in enumerate(segs):
        v0, v1 = s.verticies[0], s.verticies[1]
        heads[i, 0], heads[i, 1], heads[i, 2] = v0.x, v0.y, v0.z
        tails[i, 0], tails[i, 1], tails[i, 2] = v1.x, v1.y, v1.z
    alive = np.ones(n, dtype=bool)

    contour = Contour()
    contour.segments = [segs[0]]
    alive[0] = False
    tip = tails[0].copy()

    while True:
        idx = np.flatnonzero(alive)
        if idx.size == 0:
            break
        d_head = np.sum((heads[idx] - tip) ** 2, axis=1)
        d_tail = np.sum((tails[idx] - tip) ** 2, axis=1)
        jh = int(np.argmin(d_head))
        jt = int(np.argmin(d_tail))
        if d_head[jh] <= d_tail[jt]:
            min_d = float(d_head[jh])
            chosen = int(idx[jh])
            reverse = False
        else:
            min_d = float(d_tail[jt])
            chosen = int(idx[jt])
            reverse = True
        if min_d > _CLOSE_EPS:
            break
        if reverse:
            segs[chosen].reverse()
            h = heads[chosen].copy()
            t = tails[chosen].copy()
            heads[chosen], tails[chosen] = t, h
        # Skip if this would duplicate a segment ending at the same tip
        end = tails[chosen]
        dup = False
        for prev in contour.segments:
            pv = prev.verticies[1]
            if (
                Math.float_eq(pv.x, end[0])
                and Math.float_eq(pv.y, end[1])
                and Math.float_eq(pv.z, end[2])
            ):
                dup = True
                break
        if dup:
            logger.warning(
                "already exists a segment(s) with same start point. skipping"
            )
        else:
            contour.segments.append(segs[chosen])
        alive[chosen] = False
        tip = tails[chosen].copy()

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
