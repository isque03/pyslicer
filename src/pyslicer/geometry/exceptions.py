"""Geometry exceptions."""


class ParallelLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NonIntersectingLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class CoincidentLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
