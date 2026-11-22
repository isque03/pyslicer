"""Layer processing and multithreaded slice orchestration."""

import logging
import math
from functools import partial
from multiprocessing.dummy import Pool

from pyslicer.config import validate_model_settings
from pyslicer.geometry.line import Line
from pyslicer.limits import InputError, MAX_FACET_LAYERS, MAX_FACETS, MAX_LAYERS
from pyslicer.mesh.layer import Layer
from pyslicer.slicing.perimeters import make_perimeters, simplify_contours
from pyslicer.slicing.plane import segments_at_plane
from pyslicer.slicing.segments import join_segments

logger = logging.getLogger(__name__)


def layer_range(start, end, step):
    if not all(math.isfinite(value) for value in (start, end, step)) or step <= 0:
        raise InputError("Layer height must be positive and finite, with finite Z bounds")
    if start > end:
        raise InputError("Layer Z bounds must be ordered")
    if (end - start) / step > MAX_LAYERS - 1:
        raise InputError(f"Slicing exceeds the {MAX_LAYERS:,}-layer limit")

    layers = 0
    while start < end:
        if layers == MAX_LAYERS - 1:
            raise InputError(f"Slicing exceeds the {MAX_LAYERS:,}-layer limit")
        yield start
        next_z = start + step
        if next_z <= start:
            raise InputError("Layer height is too small to advance Z")
        start = next_z
        layers += 1
    yield end


def process_layer(zcur, model, mintolerance=0.001):
    segments_at_z = segments_at_plane(model, zcur)
    logger.debug("processLayer: found %d segments at %f ", len(segments_at_z), zcur)

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
    validate_model_settings(model)
    if len(model.facets) > MAX_FACETS:
        raise InputError(f"STL exceeds the {MAX_FACETS:,}-facet limit")
    levels = list(layer_range(model.zmin, model.zmax, model.layerHeight))
    if len(model.facets) * len(levels) > MAX_FACET_LAYERS:
        raise InputError(f"Slicing exceeds the {MAX_FACET_LAYERS:,} facet-layer work limit")

    logger.info(
        "slicing %s %f %f %f",
        model.name,
        model.zmin,
        model.zmax,
        model.layerHeight,
    )
    partial_process = partial(process_layer_print_exceptions, model=model)
    with Pool(model.processing_threads) as pool:
        return pool.map(partial_process, levels)


# Back-compat
slice = slice_model
layerRange = layer_range
