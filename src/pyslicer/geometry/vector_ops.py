"""Vector arithmetic helpers."""

from pyslicer.geometry.vertex import Vertex


def dot_product(vector_a, vector_b):
    return (
        vector_a.x * vector_b.x
        + vector_a.y * vector_b.y
        + vector_a.z * vector_b.z
    )


def subtract_vectors(vector_a, vector_b):
    return Vertex(
        vector_a.x - vector_b.x,
        vector_a.y - vector_b.y,
        vector_a.z - vector_b.z,
    )


def multiply_by_scalar(vector, scalar):
    return Vertex(vector.x * scalar, vector.y * scalar, vector.z * scalar)


def multiply_by_scalar_2d(vector, scalar):
    return Vertex(vector.x * scalar, vector.y * scalar, vector.z)


def add_vectors(vector_a, vector_b):
    return Vertex(
        vector_a.x + vector_b.x,
        vector_a.y + vector_b.y,
        vector_a.z + vector_b.z,
    )


# Back-compat aliases
dotProduct = dot_product
subtractVectors = subtract_vectors
multiplyByScalar = multiply_by_scalar
multiplyByScalar2D = multiply_by_scalar_2d
addVectors = add_vectors
