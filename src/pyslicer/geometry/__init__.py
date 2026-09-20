"""Public geometry exports."""

from pyslicer.geometry.contour import Contour
from pyslicer.geometry.exceptions import (
    CoincidentLines,
    NonIntersectingLines,
    ParallelLines,
)
from pyslicer.geometry.facet import Facet
from pyslicer.geometry.intersection import Intersection
from pyslicer.geometry.line import Line
from pyslicer.geometry.math_utils import Math
from pyslicer.geometry.vector_ops import (
    add_vectors,
    addVectors,
    dot_product,
    dotProduct,
    multiply_by_scalar,
    multiply_by_scalar_2d,
    multiplyByScalar,
    multiplyByScalar2D,
    subtract_vectors,
    subtractVectors,
)
from pyslicer.geometry.vertex import Vertex

__all__ = [
    "Contour",
    "CoincidentLines",
    "NonIntersectingLines",
    "ParallelLines",
    "Facet",
    "Intersection",
    "Line",
    "Math",
    "Vertex",
    "add_vectors",
    "addVectors",
    "dot_product",
    "dotProduct",
    "multiply_by_scalar",
    "multiplyByScalar",
    "multiply_by_scalar_2d",
    "multiplyByScalar2D",
    "subtract_vectors",
    "subtractVectors",
]
