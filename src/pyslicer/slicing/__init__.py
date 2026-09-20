"""Slicing package exports."""

from pyslicer.slicing.perimeters import make_perimeters, simplify_contours
from pyslicer.slicing.plane import slice_at, sliceAt
from pyslicer.slicing.process import slice_model
from pyslicer.slicing.segments import join_segments

# Legacy name used by CLI
slice = slice_model

__all__ = [
    "make_perimeters",
    "simplify_contours",
    "slice_at",
    "sliceAt",
    "slice_model",
    "slice",
    "join_segments",
]
