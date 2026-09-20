"""3D / 2D vertex (point)."""

from pyslicer.constants import CLIPPER_SCALE_FACTOR
from pyslicer.geometry.math_utils import Math


class Vertex:
    """A point in 3D space."""

    def __init__(self, x, y, z=0.0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __str__(self):
        return f"{id(self)} [{self.x:.15f},{self.y:.15f},{self.z:.15f}] "

    def isequal2(self, a, b):
        return Math.float_eq(a, b)

    def clipper_x(self):
        return int(self.x * CLIPPER_SCALE_FACTOR)

    def clipper_y(self):
        return int(self.y * CLIPPER_SCALE_FACTOR)

    # Back-compat
    clipperX = clipper_x
    clipperY = clipper_y

    def __eq__(self, other):
        if other is None:
            return False
        return (
            self.isequal2(self.x, other.x)
            and self.isequal2(self.y, other.y)
            and self.isequal2(self.z, other.z)
        )

    def __key(self):
        return (self.x, self.y, self.z)

    def __hash__(self):
        return hash(self.__key())
