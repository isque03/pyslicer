"""STL file I/O (ASCII and binary)."""

import logging
import math
from struct import unpack

from pyslicer.geometry.facet import Facet
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex
from pyslicer.limits import InputError, MAX_COORDINATE_MM, MAX_FACETS, MAX_STL_BYTES, open_input

logger = logging.getLogger(__name__)


def update_model_z_range(zcur, model):
    if zcur < model.zmin:
        model.zmin = zcur
    if zcur > model.zmax:
        model.zmax = zcur


def _checked_vertex(coords):
    vertex = Vertex(*coords)
    if any(
        not math.isfinite(value) or abs(value) > MAX_COORDINATE_MM
        for value in (vertex.x, vertex.y, vertex.z)
    ):
        raise InputError(f"STL coordinates must be finite and within ±{MAX_COORDINATE_MM:,} mm")
    return vertex


def _check_facet_capacity(model, count):
    if len(model.facets) + count > MAX_FACETS:
        raise InputError(f"STL exceeds the {MAX_FACETS:,}-facet limit")


def _finish_ascii_facet(model):
    if not model.facets or len(model.facets[-1].verticies) != 3:
        raise InputError("ASCII STL facet must contain three vertices")
    facet = model.facets[-1]
    if facet.edges:
        raise InputError("ASCII STL facet was already closed")
    if facet.inZPlane():
        for i in range(3):
            edge = Line.withVerticies(facet.verticies[i], facet.verticies[(i + 1) % 3])
            edge = model.edges.setdefault(edge.key(), edge)
            edge.facets.append(facet)
            facet.edges.append(edge)


def process_line(line, model):
    if isinstance(line, bytes):
        line = line.decode("utf-8", errors="replace")
    stripped = line.strip()
    if stripped == "outer loop":
        model.triangles += 1
    elif stripped == "endloop":
        model.endLoop += 1
    elif stripped.startswith("facet"):
        _check_facet_capacity(model, 1)
        model.facets.append(Facet())
    elif stripped.startswith("endfacet"):
        _finish_ascii_facet(model)
    elif stripped.startswith("vertex"):
        points = stripped.split(maxsplit=4)[1:]
        if not model.facets or len(model.facets[-1].verticies) >= 3 or len(points) != 3:
            raise InputError("Invalid vertex in ASCII STL facet")
        vertex = _checked_vertex(map(float, points))
        model.facets[-1].verticies.append(vertex)
        update_model_z_range(vertex.z, model)


def process_ascii_stl(file, model):
    size = 0
    for line in file:
        size += len(line)
        if size > MAX_STL_BYTES:
            raise InputError(f"STL file exceeds the {MAX_STL_BYTES:,}-byte limit")
        process_line(line, model)
    if any(len(facet.verticies) != 3 for facet in model.facets):
        raise InputError("ASCII STL facet must contain three vertices")


def read_triangle(f, model):
    record = f.read(50)
    if len(record) != 50:
        raise InputError("Truncated binary STL facet")
    values = unpack("<12fH", record)
    vertices = [_checked_vertex(values[offset:offset + 3]) for offset in (3, 6, 9)]
    model.facets.append(Facet.withVerticies(*vertices))
    for vertex in vertices:
        update_model_z_range(vertex.z, model)


def process_binary_stl(fp, model):
    header = fp.read(84)
    if len(header) != 84:
        raise InputError("Truncated binary STL header")
    num_triangles = unpack("<I", header[80:])[0]
    _check_facet_capacity(model, num_triangles)
    logger.info("Num Triangles: %s", num_triangles)
    for _ in range(num_triangles):
        read_triangle(fp, model)


def _looks_like_binary_stl(fp) -> bool:
    """Recognize binary STL by size, including padded files and solid headers."""
    fp.seek(0, 2)
    size = fp.tell()
    if size < 84:
        return False
    fp.seek(80)
    count_bytes = fp.read(4)
    if len(count_bytes) != 4:
        return False
    num_triangles = unpack("<I", count_bytes)[0]
    return size >= 84 + num_triangles * 50


def read_file(name, model):
    with open_input(name, MAX_STL_BYTES, "STL file") as infile:
        model.name = str(name)
        is_binary = _looks_like_binary_stl(infile)
        infile.seek(0)
        if is_binary:
            process_binary_stl(infile, model)
        else:
            process_ascii_stl(infile, model)
    if not model.facets:
        raise InputError(f"STL contained no triangles: {name}")
    model.rebuild_facet_arrays()
    return model


# Back-compat aliases
readFile = read_file
processASCIISTL = process_ascii_stl
processBinarySTL = process_binary_stl
