"""Plane–mesh intersection and Z filtering (numpy-accelerated)."""

import logging

import numpy as np

from pyslicer.geometry.line import Line
from pyslicer.geometry.math_utils import Math
from pyslicer.geometry.vertex import Vertex

logger = logging.getLogger(__name__)


def intersects(facet, z):
    if len(facet.verticies) != 3:
        raise Exception(
            f"Found facet not having 3 verticies. OMG {len(facet.verticies)} vertices found "
        )
    zs = [v.z for v in facet.verticies]
    return min(zs) <= z <= max(zs)


def find_intersecting_lines(facet, z):
    intersections = []
    v = facet.verticies
    if v[0].z > z and Math.float_eq(v[1].z, z) and Math.float_eq(v[2].z, z):
        return intersections
    if v[1].z > z and Math.float_eq(v[0].z, z) and Math.float_eq(v[2].z, z):
        return intersections
    if v[2].z > z and Math.float_eq(v[0].z, z) and Math.float_eq(v[1].z, z):
        return intersections

    for a, b in ((0, 1), (0, 2), (1, 2)):
        line = get_intersecting_line(v[a], v[b], z)
        if line is not None:
            intersections.append(line)
    return intersections


def get_intersecting_line(vertexa, vertexb, z):
    if vertexa.z == z and vertexb.z == z:
        return None
    if max(vertexa.z, vertexb.z) >= z and min(vertexb.z, vertexa.z) <= z:
        return Line.withVerticies(vertexa, vertexb)
    return None


def slice_at(model, zslice):
    """Find all facets that intersect the plane at zslice (numpy Z filter)."""
    logger.info("================= slicing at %f =================", zslice)
    if model.facet_vertices is None:
        model.rebuild_facet_arrays()

    verts = model.facet_vertices
    if verts is None or len(verts) == 0:
        return []

    zvals = verts[:, :, 2]
    zmin = zvals.min(axis=1)
    zmax = zvals.max(axis=1)
    mask = (zmin <= zslice) & (zmax >= zslice)
    indices = np.nonzero(mask)[0]
    return [model.facets[int(i)] for i in indices]


def get_intersecting_points(lines, z):
    """Intersect line segments with horizontal plane z (vectorized when possible)."""
    if not lines:
        return []

    # Batch endpoints
    p0 = np.array([[ln.verticies[0].x, ln.verticies[0].y, ln.verticies[0].z] for ln in lines])
    p1 = np.array([[ln.verticies[1].x, ln.verticies[1].y, ln.verticies[1].z] for ln in lines])
    dz = p1[:, 2] - p0[:, 2]
    points = []
    for i in range(len(lines)):
        if dz[i] == 0.0:
            points.extend(lines[i].verticies)
            continue
        s = (z - p0[i, 2]) / dz[i]
        if s < 0.0 or s > 1.0:
            raise Exception(
                "Found a non intersecting line (s not in range [0-1]) "
                "in list that should only have intersecting lines."
            )
        ps = Vertex(
            p0[i, 0] + s * (p1[i, 0] - p0[i, 0]),
            p0[i, 1] + s * (p1[i, 1] - p0[i, 1]),
            p0[i, 2] + s * (p1[i, 2] - p0[i, 2]),
        )
        if ps not in points:
            points.append(ps)
    return points


# Back-compat
sliceAt = slice_at
findIntersectingLines = find_intersecting_lines
getIntersectingLine = get_intersecting_line
getIntersectingPoints = get_intersecting_points
