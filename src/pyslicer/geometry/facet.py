"""Triangle facet."""

from pyslicer.geometry.math_utils import Math


class Facet:
    def __init__(self):
        self.verticies = []
        self.edges = []

    @classmethod
    def withVerticies(cls, v1, v2, v3):
        self = cls()
        self.verticies.append(v1)
        self.verticies.append(v2)
        self.verticies.append(v3)
        return self

    def isCoplanar(self, z):
        for vertex in self.verticies:
            if not Math.float_eq(vertex.z, z):
                return False
        return True

    def inZPlane(self):
        return self.isCoplanar(self.verticies[0].z)
