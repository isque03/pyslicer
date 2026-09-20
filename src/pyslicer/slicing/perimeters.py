"""Perimeter generation via pyclipper offsets."""

import logging

from pyslicer import clipper_ops
from pyslicer.geometry.contour import Contour
from pyslicer.geometry.line import Line

logger = logging.getLogger(__name__)


def simplify_contour(contour, tolerance=0.001, slicer=False):
    if slicer:
        result = clipper_ops.simplify_polygon(contour.to_path())
        if len(result) == 0:
            return None
        return Contour.from_path(result[0], contour.zlevel)

    new_contour = Contour()
    new_contour.segments = []
    new_contour.zlevel = contour.zlevel
    while True:
        i = 0
        while i < len(contour.segments) - 1:
            l = Line.withVerticies(
                contour.segments[i].verticies[0],
                contour.segments[i + 1].verticies[1],
            )
            distance = l.distance(contour.segments[i].verticies[1])
            if distance != -1 and distance <= tolerance:
                new_contour.segments.append(l)
                i += 2
            else:
                new_contour.segments.append(contour.segments[i])
                i += 1
        if i == len(contour.segments) - 1:
            new_contour.segments.append(contour.segments[-1])
        if len(new_contour.segments) == len(contour.segments):
            break
        contour = new_contour
        new_contour = Contour()
        new_contour.zlevel = contour.zlevel
        new_contour.segments = []
    return new_contour


def simplify_contours(contours, tolerance=0.001):
    return [simplify_contour(c, tolerance=tolerance) for c in contours]


def make_perimeters(pre_contours, model, zcur):
    paths = [c.to_path() for c in pre_contours if c is not None]
    if not paths:
        return []
    unioned = clipper_ops.union_polygons(paths)
    perimeters = []
    offset_amount = model.offset() / 2.0
    for _ in range(model.number_perimeters + 1):
        perimeter = []
        simple_offset = clipper_ops.offset_polygons(unioned, -offset_amount)
        for poly in simple_offset:
            perimeter.append(Contour.from_path(poly, zcur))
        perimeters.append(perimeter)
        offset_amount += model.offset() * model.perimeter_overlap_percent
    return perimeters


# Back-compat
simplifyContour = simplify_contour
simplifyContours = simplify_contours
