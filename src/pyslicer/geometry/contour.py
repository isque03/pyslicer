"""Closed or open contour (polygon ring)."""

import logging
import sys

from pyslicer.constants import CLIPPER_SCALE_FACTOR
from pyslicer.geometry.exceptions import (
    CoincidentLines,
    NonIntersectingLines,
    ParallelLines,
)
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex

logger = logging.getLogger(__name__)


class Contour:
    def __init__(self):
        self.segments = []
        self.zlevel = 0.0
        self.closed = False

    @classmethod
    def from_path(cls, path, zlevel):
        """Build a contour from a list of (x, y) float coordinates."""
        self = cls()
        self.zlevel = zlevel
        if not path:
            return self
        n = len(path)
        for i in range(n):
            p = path[i]
            p1 = path[(i + 1) % n]
            x0, y0 = float(p[0]), float(p[1])
            x1, y1 = float(p1[0]), float(p1[1])
            self.segments.append(
                Line.withVerticies(Vertex(x0, y0, zlevel), Vertex(x1, y1, zlevel))
            )
        return self

    @classmethod
    def withClipperPolygon(cls, clipper_polygon, zlevel):
        """Build from clipper-scaled integer (x,y) path. Prefer Contour.from_path for mm floats."""
        if not clipper_polygon:
            return cls.from_path([], zlevel)
        sample = clipper_polygon[0]
        if hasattr(sample, "x"):
            path = [
                (float(p.x) / CLIPPER_SCALE_FACTOR, float(p.y) / CLIPPER_SCALE_FACTOR)
                for p in clipper_polygon
            ]
        else:
            path = [
                (float(p[0]) / CLIPPER_SCALE_FACTOR, float(p[1]) / CLIPPER_SCALE_FACTOR)
                for p in clipper_polygon
            ]
        return cls.from_path(path, zlevel)

    def __str__(self):
        return f"Contour: {id(self)}"

    def is_closed(self):
        if len(self.segments) == 0:
            return False
        return self.segments[0].verticies[0] == self.segments[-1].verticies[1]

    def has_disconnected_segments(self):
        if not self.is_closed():
            print("## Segment not closed, therefore it's disconnected.")
        for x in range(len(self.segments)):
            if x < len(self.segments) - 1:
                cur = self.segments[x]
                nxt = self.segments[x + 1]
                if not cur.verticies[1] == nxt.verticies[0]:
                    print("## Contour has vertices that are not within tolerance.")

    def winding_area(self):
        """Signed shoelace area. Positive ⇒ counter-clockwise; negative ⇒ clockwise."""
        total = 0.0
        for segment in self.segments:
            x0, y0 = segment.verticies[0].x, segment.verticies[0].y
            x1, y1 = segment.verticies[1].x, segment.verticies[1].y
            total += x0 * y1 - x1 * y0
        return total / 2.0

    def clockwise(self):
        """True when polygon winds clockwise (negative shoelace area)."""
        return self.winding_area() < 0.0

    def intersect_brute_force(self, other):
        points = []
        for sega in self.segments:
            for segb in other.segments:
                try:
                    points.append(sega.intersect2D(segb))
                except NonIntersectingLines:
                    pass
        return points

    def maybe_intersect(self, other):
        other_box = other.bounding_box()
        box = self.bounding_box()
        if box.verticies[1].x < other_box.verticies[0].x:
            return False
        if other_box.verticies[1].x < box.verticies[0].x:
            return False
        if box.verticies[1].y < other_box.verticies[0].y:
            return False
        if other_box.verticies[1].y < box.verticies[0].y:
            return False
        return True

    def bounding_box(self):
        current_min_x = sys.float_info.max
        current_max_x = -sys.float_info.max
        current_min_y = current_min_x
        current_max_y = current_max_x
        z = self.segments[0].verticies[0].z
        for line in self.segments:
            for vertex in line.verticies:
                current_min_x = min(current_min_x, vertex.x)
                current_min_y = min(current_min_y, vertex.y)
                current_max_x = max(current_max_x, vertex.x)
                current_max_y = max(current_max_y, vertex.y)
        return Line.withVerticies(
            Vertex(current_min_x, current_min_y, z),
            Vertex(current_max_x, current_max_y, z),
        )

    def to_path(self):
        """Return float (x, y) path for clipper ops (one point per segment start + close)."""
        points = []
        for line in self.segments:
            points.append((line.verticies[0].x, line.verticies[0].y))
        if self.segments:
            line = self.segments[-1]
            points.append((line.verticies[1].x, line.verticies[1].y))
        return points

    def clipper_points(self):
        """Legacy: scaled integer (x, y) tuples for clipper."""
        return [
            (v.clipper_x(), v.clipper_y())
            for v in (
                [s.verticies[0] for s in self.segments]
                + ([self.segments[-1].verticies[1]] if self.segments else [])
            )
        ]
