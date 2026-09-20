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
    logger.debug("================= slicing at %f =================", zslice)
    if model.facet_vertices is None:
        model.rebuild_facet_arrays()

    verts = model.facet_vertices
    if verts is None or len(verts) == 0:
        return []

    zmin = model.facet_zmin
    zmax = model.facet_zmax
    if zmin is None or zmax is None:
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


_EPS = 1e-7
_EDGE_PAIRS = ((0, 1), (0, 2), (1, 2))


def segments_at_plane(model, zslice):
    """Cut the mesh at ``zslice`` and return segment ``Line``s.

    Uses cached ``facet_vertices`` / Z bounds and scalar edge tests so we avoid
    building per-edge ``Line`` objects (same semantics as the old facet loop).
    """
    if model.facet_vertices is None:
        model.rebuild_facet_arrays()

    verts = model.facet_vertices
    if verts is None or len(verts) == 0:
        return []

    z = float(zslice)
    zmin = model.facet_zmin
    zmax = model.facet_zmax
    cand = np.nonzero((zmin <= z) & (zmax >= z))[0]
    if cand.size == 0:
        return []

    zvals = verts[cand, :, 2]
    coplanar = (
        (np.abs(zvals[:, 0] - z) <= _EPS)
        & (np.abs(zvals[:, 1] - z) <= _EPS)
        & (np.abs(zvals[:, 2] - z) <= _EPS)
    )
    on = np.abs(zvals - z) <= _EPS
    skip_base = (
        ((zvals[:, 0] > z) & on[:, 1] & on[:, 2])
        | ((zvals[:, 1] > z) & on[:, 0] & on[:, 2])
        | ((zvals[:, 2] > z) & on[:, 0] & on[:, 1])
    )
    idx = cand[~(coplanar | skip_base)]
    if idx.size == 0:
        return []

    segments = []
    append = segments.append
    with_verts = Line.withVerticies
    vertex = Vertex
    for i in idx:
        tri = verts[i]
        pts = []
        for a, b in _EDGE_PAIRS:
            pa = tri[a]
            pb = tri[b]
            za = float(pa[2])
            zb = float(pb[2])
            if za == z and zb == z:
                continue
            if not (max(za, zb) >= z and min(za, zb) <= z):
                continue
            dz = zb - za
            if dz == 0.0:
                px, py, pz = float(pa[0]), float(pa[1]), float(pa[2])
            else:
                s = (z - za) / dz
                if s < 0.0 or s > 1.0:
                    raise Exception(
                        "Found a non intersecting line (s not in range [0-1]) "
                        "in list that should only have intersecting lines."
                    )
                px = float(pa[0] + s * (pb[0] - pa[0]))
                py = float(pa[1] + s * (pb[1] - pa[1]))
                pz = float(pa[2] + s * (pb[2] - pa[2]))
            if any(
                abs(q[0] - px) <= _EPS
                and abs(q[1] - py) <= _EPS
                and abs(q[2] - pz) <= _EPS
                for q in pts
            ):
                continue
            pts.append((px, py, pz))
        if len(pts) == 2:
            append(
                with_verts(
                    vertex(pts[0][0], pts[0][1], pts[0][2]),
                    vertex(pts[1][0], pts[1][1], pts[1][2]),
                )
            )
    return segments


# Back-compat
sliceAt = slice_at
findIntersectingLines = find_intersecting_lines
getIntersectingLine = get_intersecting_line
getIntersectingPoints = get_intersecting_points
segmentsAtPlane = segments_at_plane
