"""CLI entry point for pyslicer."""

import argparse
import logging
import sys
from datetime import datetime
from operator import attrgetter

from pyslicer import clipper_ops
from pyslicer.geometry.contour import Contour
from pyslicer.infill import simple_linear_infill
from pyslicer.mesh import Model, read_file
from pyslicer.slicing import slice_model
from pyslicer.timer import Timer


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="An STL slicer for producing 3D printer gcode."
    )
    parser.add_argument("stl", help="STL file to be sliced.")
    parser.add_argument("output", help="Output path for generated gcode.")
    parser.add_argument(
        "-p",
        "--perimeters_only",
        help="Generate only perimeters, no infill.",
        action="store_true",
    )
    parser.add_argument(
        "-a",
        "--append_perimeters",
        help="Include original stl perimeters with no offsetting in the output.",
        action="store_true",
    )
    parser.add_argument("-v", "--verbose", help="Be verbose.", action="store_true")
    parser.add_argument(
        "-o",
        "--perimeter_overlap_percent",
        help="Percent overlap between perimeters. Smaller results in more overlap.",
        type=float,
        default=1.0,
    )
    parser.add_argument(
        "-n",
        "--num-perimeters",
        help="Minimum number of perimeters.",
        type=int,
        default=3,
    )
    parser.add_argument(
        "-d",
        "--filament-diameter",
        help="Measured diameter of filament to at least two decimal places.",
        type=float,
        default=1.75,
    )
    parser.add_argument(
        "-l",
        "--layer-height",
        help="Slicing layer height",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--html-preview",
        metavar="PATH",
        help="Also write an HTML toolpath preview to this path.",
    )
    return parser


def _compute_roofs_and_overhangs(layers):
    for layer_num in range(1, len(layers)):
        zcur = layers[layer_num].z
        if not layers[layer_num].perimeters:
            continue
        current = [c.to_path() for c in layers[layer_num].perimeters[-1]]
        previous = []
        if layers[layer_num - 1].perimeters:
            previous = [
                c.to_path()
                for c in layers[layer_num - 1].perimeters[-1]
                if c.is_closed()
            ]
        nxt = []
        if layer_num < len(layers) - 2 and layers[layer_num + 1].perimeters:
            nxt = [c.to_path() for c in layers[layer_num + 1].perimeters[-1]]

        roofs = clipper_ops.difference_polygons(previous, current)
        for poly in roofs:
            layers[layer_num - 1].roofs.append(Contour.from_path(poly, zcur))

        previous_minus_roof = clipper_ops.difference_polygons(previous, roofs)
        for poly in previous_minus_roof:
            layers[layer_num - 1].normalInfill.append(Contour.from_path(poly, zcur))

        overhang = clipper_ops.difference_polygons(nxt, current)
        for poly in overhang:
            if layer_num + 1 < len(layers):
                layers[layer_num + 1].overhang.append(Contour.from_path(poly, zcur))


def run(args):
    logger = logging.getLogger()
    ch = logging.StreamHandler(sys.stdout)
    logger.addHandler(ch)
    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    sys.setrecursionlimit(6000)
    logger.info("Reading...")
    with Timer() as total_time:
        with Timer() as read_time:
            model = Model()
            model.perimeters_only = args.perimeters_only
            model.append_perimeters = args.append_perimeters
            model.perimeter_overlap_percent = args.perimeter_overlap_percent
            model.number_perimeters = args.num_perimeters
            model.filament_diameter = args.filament_diameter
            model.layerHeight = args.layer_height
            read_file(args.stl, model)
        logger.info("File read took %s seconds", read_time.secs)

        logger.info("Slicing...")
        with Timer() as slice_time:
            layers = slice_model(model)
        logger.warning(
            "Slicing took %s seconds. %s layers", slice_time.secs, len(layers)
        )

        layers.sort(key=attrgetter("z"))
        _compute_roofs_and_overhangs(layers)

        logger.info("Starting infill...")
        with Timer() as infill_time:
            for layer_num in range(len(layers)):
                z = layers[layer_num].z
                angle = (
                    model.infill_angle if layer_num % 2 else -model.infill_angle
                )
                perimeters = layers[layer_num].perimeters
                if not perimeters or not perimeters[-1]:
                    layers[layer_num].infill = []
                    continue
                if layer_num < 4 or layer_num > len(layers) - 4:
                    layers[layer_num].infill = simple_linear_infill(
                        perimeters[-1],
                        z,
                        spacing=model.offset(),
                        min_extrude=1.0,
                        angle=angle,
                        name="normal",
                    )
                else:
                    infill = simple_linear_infill(
                        layers[layer_num].normalInfill,
                        z,
                        spacing=2.0,
                        min_extrude=1.0,
                        angle=angle,
                        name="normal",
                    )
                    layers[layer_num].infill = infill
                    roof_infill = simple_linear_infill(
                        layers[layer_num].roofs,
                        z,
                        spacing=model.offset(),
                        min_extrude=0.5,
                        angle=angle,
                        name="roof",
                    )
                    layers[layer_num].infill.extend(roof_infill)
        logger.info("Infill time %s", infill_time.secs)

        with Timer() as model_time:
            for result in layers:
                zcur = result.z
                model.layers.append(result)
                for contour in result.contours:
                    contour.zlevel = zcur
                    model.contours.append(contour)
                model.infillAtLayer[zcur] = result.infill
        logger.warning("Model setup took %s seconds", model_time.secs)

        with Timer() as gcode_time:
            model.writeGCode(args.output)
        logger.warning("GCode generation took %s seconds", gcode_time.secs)

        if args.html_preview:
            from pyslicer.preview import write_gcode_preview

            preview_path = write_gcode_preview(
                args.output,
                args.html_preview,
                subtitle=f"From {args.stl}",
                nozzle_diameter=model.nozzle_diameter,
                layer_height=model.layerHeight,
            )
            logger.info("HTML preview wrote %s", preview_path)

    logger.info("TOTAL TIME %s seconds", total_time.secs)
    logger.info(
        "facets: %d  facet end: %d  zmin: %f zmax: %f ",
        len(model.facets),
        model.endLoop,
        model.zmin,
        model.zmax,
    )


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    main()
