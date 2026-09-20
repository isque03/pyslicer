"""Segment joining / contour stitching (numpy closest-endpoint)."""

import logging
import sys

import numpy as np

from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line

logger = logging.getLogger(__name__)


def find_closest_segment(segment, segments):
    """Return the open segment whose endpoint is closest to segment's tail."""
    if not segments:
        return None

    tip = np.array(
        [segment.verticies[1].x, segment.verticies[1].y, segment.verticies[1].z],
        dtype=np.float64,
    )
    heads = np.array(
        [[s.verticies[0].x, s.verticies[0].y, s.verticies[0].z] for s in segments],
        dtype=np.float64,
    )
    tails = np.array(
        [[s.verticies[1].x, s.verticies[1].y, s.verticies[1].z] for s in segments],
        dtype=np.float64,
    )
    d_head = np.sum((heads - tip) ** 2, axis=1)
    d_tail = np.sum((tails - tip) ** 2, axis=1)
    i_head = int(np.argmin(d_head))
    i_tail = int(np.argmin(d_tail))
    if d_head[i_head] <= d_tail[i_tail]:
        min_d = float(d_head[i_head])
        closest = i_head
        reverse = False
    else:
        min_d = float(d_tail[i_tail])
        closest = i_tail
        reverse = True

    if closest > -1 and min_d <= 0.1:
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


def join_segments(segments):
    """Join an unordered set of segments into a closed contour."""
    if segments is None or len(segments) == 0:
        return None
    contour = Contour()
    contour.segments = []
    contour.segments.append(segments[0])
    last_segment = segments[0]
    segments.remove(last_segment)
    while True:
        connected = find_closest_segment(last_segment, segments)
        if connected is None:
            if last_segment.verticies[1] == contour.segments[0].verticies[0]:
                contour.closed = True
            break
        v1 = connected.verticies[1]
        duplicate_starts = [line for line in contour.segments if line.verticies[1] == v1]
        if len(duplicate_starts) > 0:
            logger.warning(
                "already exists a segment(s) with same start point. skipping"
            )
        else:
            contour.segments.append(connected)
        segments.remove(connected)
        last_segment = connected
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
