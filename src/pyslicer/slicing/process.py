"""Layer processing and multithreaded slice orchestration."""

import logging
from functools import partial
from multiprocessing.dummy import Pool

from pyslicer.geometry.line import Line
from pyslicer.mesh.layer import Layer
from pyslicer.slicing.perimeters import make_perimeters, simplify_contours
from pyslicer.slicing.plane import (
    find_intersecting_lines,
    get_intersecting_points,
    slice_at,
)
from pyslicer.slicing.segments import join_segments

logger = logging.getLogger(__name__)


def layer_range(start, end, step):
    while start < end:
        yield start
        start += step
    if start > end:
        start = end
    yield start


def process_layer(zcur, model, mintolerance=0.001):
    facets = slice_at(model, zcur)
    logger.debug("processLayer: found %d facets at %f ", len(facets), zcur)
    segments_at_z = []

    for facet in facets:
        if facet.isCoplanar(zcur):
            continue
        lines = find_intersecting_lines(facet, zcur)
        points = get_intersecting_points(lines, zcur)
        if len(points) == 2:
            segment = Line()
            segment.verticies.extend(points)
            segments_at_z.append(segment)

    contours_from_stl = []
    segments_at_z.reverse()
    open_contours = []
    while len(segments_at_z) > 0:
        contour = join_segments(segments_at_z)
        if contour is None:
            break
        contour.zlevel = zcur
        if contour.closed:
            contours_from_stl.append(contour)
        else:
            open_contours.append(contour)

    if open_contours:
        gaps = []
        for oc in open_contours:
            gap = Line.withVerticies(
                oc.segments[-1].verticies[1], oc.segments[0].verticies[0]
            ).magnitude()
            gaps.append(gap)
        raise ValueError(
            f"Open contours at z={zcur}: {len(open_contours)} contour(s) "
            f"with gap magnitudes {gaps}. Refusing to invent closing segments."
        )

    pre_contours = simplify_contours(
        contours_from_stl,
        tolerance=model.nozzle_diameter / model.simplification_factor,
    )
    perimeters = make_perimeters(pre_contours, model, zcur)

    layer = Layer()
    layer.contours = []
    layer.perimeters = perimeters
    layer.infill = []
    layer.z = zcur
    return layer


def process_layer_print_exceptions(zcur, model, mintolerance=0.001):
    return process_layer(zcur, model, mintolerance)


def slice_model(model):
    logger.info(
        "slicing %s %f %f %f",
        model.name,
        model.zmin,
        model.zmax,
        model.layerHeight,
    )
    pool = Pool(model.processing_threads)
    partial_process = partial(process_layer_print_exceptions, model=model)
    results = pool.map(
        partial_process, layer_range(model.zmin, model.zmax, model.layerHeight)
    )
    pool.close()
    pool.join()
    return results


# Back-compat
slice = slice_model
layerRange = layer_range
