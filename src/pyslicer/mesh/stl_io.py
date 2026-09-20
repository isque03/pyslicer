"""STL file I/O (ASCII and binary)."""

import logging
from struct import unpack

from pyslicer.geometry.facet import Facet
from pyslicer.geometry.line import Line
from pyslicer.geometry.vertex import Vertex

logger = logging.getLogger(__name__)


def update_model_z_range(zcur, model):
    if zcur < model.zmin:
        model.zmin = zcur
    if zcur > model.zmax:
        model.zmax = zcur


def process_line(line, model):
    if isinstance(line, bytes):
        line = line.decode("utf-8", errors="replace")
    stripped = line.strip()
    if stripped == "outer loop":
        model.triangles += 1
    elif stripped == "endloop":
        model.endLoop += 1
    elif stripped.startswith("facet"):
        facet = Facet()
        model.facets.append(facet)
    elif stripped.startswith("endfacet"):
        if model.facets[-1].isCoplanar(model.facets[-1].verticies[0].z):
            for i in range(3):
                nxt = i + 1
                if nxt > 2:
                    nxt = 0
                edge = Line.withVerticies(
                    model.facets[-1].verticies[i],
                    model.facets[-1].verticies[nxt],
                )
                try:
                    edge = model.edges[edge.key()]
                except KeyError:
                    model.edges[edge.key()] = edge
                if model.facets[-1].inZPlane():
                    edge.facets.append(model.facets[-1])
                model.facets[-1].edges.append(edge)
    elif stripped.startswith("vertex"):
        coords = stripped.split()
        v = Vertex(float(coords[1]), float(coords[2]), float(coords[3]))
        if model.facets[-1].verticies is None:
            model.facets[-1].verticies = []
        model.facets[-1].verticies.append(v)
        update_model_z_range(float(coords[3]), model)


def process_ascii_stl(file, model):
    for line in file:
        process_line(line, model)


def read_point(f):
    vertex = unpack("<fff", f.read(12))
    return Vertex(vertex[0], vertex[1], vertex[2])


def read_triangle(f, model):
    f.seek(12, 1)
    p1 = read_point(f)
    p2 = read_point(f)
    p3 = read_point(f)
    f.seek(2, 1)
    facet = Facet()
    facet.verticies.append(p1)
    facet.verticies.append(p2)
    facet.verticies.append(p3)
    model.facets.append(facet)
    update_model_z_range(p1.z, model)
    update_model_z_range(p2.z, model)
    update_model_z_range(p3.z, model)


def process_binary_stl(fp, model):
    header = fp.read(80)
    num_triangles = unpack("<i", fp.read(4))[0]
    logger.info(header)
    logger.info("Num Triangles: %s", num_triangles)
    for _ in range(num_triangles):
        read_triangle(fp, model)


def _looks_like_binary_stl(fp) -> bool:
    """True when the file is a binary STL (even if the header says 'solid').

    Accepts trailing padding after the triangle records. Rejects sizes that
    cannot hold the declared triangle count.
    """
    fp.seek(0, 2)
    size = fp.tell()
    if size < 84:
        return False
    fp.seek(80)
    num_triangles = unpack("<I", fp.read(4))[0]
    # Guard absurd counts that would imply a multi-GB triangle block
    if num_triangles > 100_000_000:
        return False
    expected = 84 + num_triangles * 50
    return size >= expected


def read_file(name, model):
    with open(name, "rb") as infile:
        model.name = infile.name
        if _looks_like_binary_stl(infile):
            infile.seek(0)
            process_binary_stl(infile, model)
        else:
            infile.seek(0)
            process_ascii_stl(infile, model)
    if not model.facets:
        raise ValueError(f"STL contained no triangles: {name}")
    model.rebuild_facet_arrays()
    return model


# Back-compat aliases
readFile = read_file
processASCIISTL = process_ascii_stl
processBinarySTL = process_binary_stl
