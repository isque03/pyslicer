"""Infill intersection bookkeeping."""


class Intersection:
    def __init__(self):
        self.uparam = -1.0
        self.contour = None
        self.segmentIndex = -1
        self.row = -1
        self.ray = None
        self.used = False
        self.z = 0.0
        self.intersectionIndex = -1

    def point(self):
        return self.ray.pointOnLine(self.uparam)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}: u: {self.uparam} "
            f"segmentIndex: {self.segmentIndex} row: {self.row} "
            f"contour: {self.contour}"
        )
