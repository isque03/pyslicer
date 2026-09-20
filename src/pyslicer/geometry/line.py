"""Line segment geometry."""

import math

from pyslicer.geometry.exceptions import (
    CoincidentLines,
    NonIntersectingLines,
    ParallelLines,
)
from pyslicer.geometry.math_utils import Math
from pyslicer.geometry.vector_ops import (
    add_vectors,
    multiply_by_scalar,
    multiply_by_scalar_2d,
    subtract_vectors,
)
from pyslicer.geometry.vertex import Vertex


class Line:
    """A line segment of two points."""

    def __init__(self):
        self.verticies = []
        self.facets = []

    @classmethod
    def withVerticies(cls, va, vb):
        self = cls()
        self.verticies.append(va)
        self.verticies.append(vb)
        return self

    def key(self):
        format_one = (
            "[%.6f,%.6f,%.6f] [%.6f,%.6f,%.6f]"
            % (
                self.verticies[0].x,
                self.verticies[0].y,
                self.verticies[0].z,
                self.verticies[1].x,
                self.verticies[1].y,
                self.verticies[1].z,
            )
        )
        format_two = (
            "[%.6f,%.6f,%.6f] [%.6f,%.6f,%.6f]"
            % (
                self.verticies[1].x,
                self.verticies[1].y,
                self.verticies[1].z,
                self.verticies[0].x,
                self.verticies[0].y,
                self.verticies[0].z,
            )
        )
        if self.verticies[0].x > self.verticies[1].x:
            return format_one
        if self.verticies[0].x < self.verticies[1].x:
            return format_two
        if self.verticies[0].y > self.verticies[1].y:
            return format_one
        return format_two

    def __str__(self):
        return (
            "%s [%.15f,%.15f,%.15f] [%.15f,%.15f,%.15f]"
            % (
                id(self),
                self.verticies[0].x,
                self.verticies[0].y,
                self.verticies[0].z,
                self.verticies[1].x,
                self.verticies[1].y,
                self.verticies[1].z,
            )
        )

    def __eq__(self, other):
        return (
            self.verticies[0] == other.verticies[0]
            and self.verticies[1] == other.verticies[1]
        ) or (
            self.verticies[0] == other.verticies[1]
            and self.verticies[1] == other.verticies[0]
        )

    def __ne__(self, other):
        return not self.__eq__(other)

    def reverse(self):
        if self.verticies is None:
            return
        self.verticies[0], self.verticies[1] = self.verticies[1], self.verticies[0]

    def magnitude(self):
        if len(self.verticies) != 2:
            raise Exception(
                "Cannot calculate magnitude of segment not having exactly two vertices"
            )
        return math.sqrt(self.magnitudeSquared())

    def magnitudeSquared(self):
        return (
            (self.verticies[0].x - self.verticies[1].x) ** 2
            + (self.verticies[0].y - self.verticies[1].y) ** 2
            + (self.verticies[0].z - self.verticies[1].z) ** 2
        )

    def distance(self, vertex):
        magsqr = self.magnitudeSquared()
        if magsqr <= 0.0:
            return 0.0
        u = (
            (vertex.x - self.verticies[0].x)
            * (self.verticies[1].x - self.verticies[0].x)
            + (vertex.y - self.verticies[0].y)
            * (self.verticies[1].y - self.verticies[0].y)
        ) / magsqr
        if u < 0.0 or u > 1.0:
            return -1
        x = self.verticies[0].x + u * (self.verticies[1].x - self.verticies[0].x)
        y = self.verticies[0].y + u * (self.verticies[1].y - self.verticies[0].y)
        point = Vertex(x, y, vertex.z)
        tmp = Line.withVerticies(vertex, point)
        return tmp.magnitude()

    def pointOnLine(self, u):
        return add_vectors(
            self.verticies[0],
            multiply_by_scalar(
                subtract_vectors(self.verticies[1], self.verticies[0]), u
            ),
        )

    def intersect2D(self, line):
        ua = self.intersect2DU(line)
        x1, x2 = self.verticies[0].x, self.verticies[1].x
        y1, y2 = self.verticies[0].y, self.verticies[1].y
        return Vertex(x1 + ua * (x2 - x1), y1 + ua * (y2 - y1), self.verticies[0].z)

    def intersect2DUExtend(self, line):
        x4, x3 = line.verticies[1].x, line.verticies[0].x
        x2, x1 = self.verticies[1].x, self.verticies[0].x
        y4, y3 = line.verticies[1].y, line.verticies[0].y
        y2, y1 = self.verticies[1].y, self.verticies[0].y

        denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
        numeratora = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)
        numeratorb = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)

        if (
            Math.float_eq(denominator, 0.0)
            and Math.float_eq(numeratora, 0.0)
            and Math.float_eq(numeratorb, 0.0)
        ):
            raise CoincidentLines("Lines are coincient")
        if Math.float_eq(denominator, 0.0):
            raise ParallelLines("Lines are parallel and do not intersect")

        return numeratora / denominator, numeratorb / denominator

    def intersect2DU(self, line):
        ua, ub = self.intersect2DUExtend(line)
        if ua < 0 or ua > 1 or ub < 0 or ub > 1:
            raise NonIntersectingLines("Segments do not intersect with 0<u<1")
        return ua

    def normal(self):
        dx = self.verticies[1].x - self.verticies[0].x
        dy = self.verticies[1].y - self.verticies[0].y
        return Line.withVerticies(
            Vertex(-dy, dx, self.verticies[0].z),
            Vertex(dy, -dx, self.verticies[0].z),
        )

    def unitVector(self):
        mag = self.magnitude()
        v1 = Vertex(
            self.verticies[0].x / mag,
            self.verticies[0].y / mag,
            self.verticies[0].z / mag,
        )
        v2 = Vertex(
            self.verticies[1].x / mag,
            self.verticies[1].y / mag,
            self.verticies[1].z / mag,
        )
        return Line.withVerticies(v1, v2)

    def scale2D(self, scale_factor):
        self.verticies[0] = multiply_by_scalar_2d(self.verticies[0], scale_factor)
        self.verticies[1] = multiply_by_scalar_2d(self.verticies[1], scale_factor)

    def copy(self):
        v1 = Vertex(self.verticies[0].x, self.verticies[0].y, self.verticies[0].z)
        v2 = Vertex(self.verticies[1].x, self.verticies[1].y, self.verticies[1].z)
        return Line.withVerticies(v1, v2)

    def translate(self, vector):
        self.verticies[0].x += vector.x
        self.verticies[0].y += vector.y
        self.verticies[0].z += vector.z
        self.verticies[1].x += vector.x
        self.verticies[1].y += vector.y
        self.verticies[1].z += vector.z
        return self

    def offset(self, distance):
        n = self.normal()
        p = n.pointOnLine(0.5)
        start_point = 1 if distance < 0 else 0
        t = Line.withVerticies(p, n.verticies[start_point]).unitVector()
        t.scale2D(distance)
        return self.translate(Vertex(t.verticies[1].x, t.verticies[1].y))

    def rotate(self, angle, point):
        angle_rad = math.radians(angle)
        c, s = math.cos(angle_rad), math.sin(angle_rad)
        ox, oy = point.x, point.y
        for i in (0, 1):
            dx = self.verticies[i].x - ox
            dy = self.verticies[i].y - oy
            self.verticies[i].x = c * dx - s * dy + ox
            self.verticies[i].y = s * dx + c * dy + oy
