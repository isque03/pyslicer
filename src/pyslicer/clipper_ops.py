"""Thin pyclipper adapter for polygon boolean ops and offsets.

All public functions take model-unit (mm) float paths. Scaled clipper
integers are never accepted — call Contour.clipper_points only for
legacy introspection, not as clipper_ops input.
"""

import pyclipper

from pyslicer.constants import CLIPPER_SCALE_FACTOR

SCALE = CLIPPER_SCALE_FACTOR


def _to_paths(polygons):
    """Convert mm float paths to clipper integers. Empty → []."""
    if not polygons:
        return []
    return [pyclipper.scale_to_clipper(path, SCALE) for path in polygons]


def _from_paths(paths):
    """Return float (x, y) paths in mm."""
    if not paths:
        return []
    return pyclipper.scale_from_clipper(paths, SCALE)


def union_polygons(paths):
    """Union subject paths. paths: list of float (x,y) rings in mm."""
    if not paths:
        return []
    pc = pyclipper.Pyclipper()
    pc.AddPaths(_to_paths(paths), pyclipper.PT_SUBJECT, True)
    result = pc.Execute(
        pyclipper.CT_UNION, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD
    )
    return _from_paths(result)


def difference_polygons(subject_paths, clip_paths):
    if not subject_paths:
        return []
    pc = pyclipper.Pyclipper()
    pc.AddPaths(_to_paths(subject_paths), pyclipper.PT_SUBJECT, True)
    if clip_paths:
        pc.AddPaths(_to_paths(clip_paths), pyclipper.PT_CLIP, True)
    result = pc.Execute(
        pyclipper.CT_DIFFERENCE, pyclipper.PFT_EVENODD, pyclipper.PFT_EVENODD
    )
    return _from_paths(result)


def offset_polygons(paths, delta):
    """Offset by delta in mm. Negative shrinks."""
    if not paths:
        return []
    pco = pyclipper.PyclipperOffset()
    for path in _to_paths(paths):
        pco.AddPath(path, pyclipper.JT_MITER, pyclipper.ET_CLOSEDPOLYGON)
    result = pco.Execute(int(delta * SCALE))
    return _from_paths(result)


def simplify_polygon(path):
    """Simplify a single mm float path; return list of float paths."""
    if not path:
        return []
    scaled = pyclipper.scale_to_clipper(path, SCALE)
    result = pyclipper.SimplifyPolygon(scaled, pyclipper.PFT_EVENODD)
    return _from_paths(result)


def clean_polygons(paths, distance=1.415):
    if not paths:
        return []
    cleaned = pyclipper.CleanPolygons(_to_paths(paths), distance)
    return _from_paths(cleaned)
