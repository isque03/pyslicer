"""CLI entry point for pyslicer."""

import argparse
import logging
import sys
from operator import attrgetter

from pyslicer import clipper_ops
from pyslicer.config import (
    apply_settings_to_model,
    canonicalize_config,
    load_layered_config,
    parse_config_paths,
)
from pyslicer.geometry.contour import Contour
from pyslicer.infill import simple_linear_infill
from pyslicer.mesh import Model, read_file
from pyslicer.slicing import slice_model
from pyslicer.timer import Timer

# Seeded for every CLI run so --config and no-config share the same baseline.
_CLI_HISTORICAL_DEFAULTS = {
    "perimeters_only": False,
    "append_perimeters": False,
    "perimeter_overlap_percent": 1.0,
    "number_perimeters": 3,
    "filament_diameter": 1.75,
    "layerHeight": 0.1,
}

_CLI_SETTING_DESTS = frozenset(
    {
        "perimeters_only",
        "append_perimeters",
        "perimeter_overlap_percent",
        "num_perimeters",
        "filament_diameter",
        "layer_height",
        "outer_speed",
        "inner_speed",
        "infill_speed",
        "max_corner_speed",
        "max_accel",
        "max_jerk",
        "min_corner_angle",
    }
)


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="An STL slicer for producing 3D printer gcode."
    )
    parser.add_argument("stl", help="STL file to be sliced.")
    parser.add_argument("output", help="Output path for generated gcode.")
    parser.add_argument(
        "--config",
        metavar="PATHS",
        help=(
            "Comma-separated YAML config files, layered left-to-right "
            "(later files override earlier). CLI flags always win."
        ),
    )
    parser.add_argument(
        "-p",
        "--perimeters-only",
        "--perimeters_only",
        dest="perimeters_only",
        help="Generate only perimeters, no infill.",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-a",
        "--append-perimeters",
        "--append_perimeters",
        dest="append_perimeters",
        help="Include original stl perimeters with no offsetting in the output.",
        action=argparse.BooleanOptionalAction,
        default=argparse.SUPPRESS,
    )
    parser.add_argument("-v", "--verbose", help="Be verbose.", action="store_true")
    parser.add_argument(
        "-o",
        "--perimeter_overlap_percent",
        help="Percent overlap between perimeters. Smaller results in more overlap.",
        type=float,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-n",
        "--num-perimeters",
        help="Minimum number of perimeters.",
        type=int,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-d",
        "--filament-diameter",
        help="Measured diameter of filament to at least two decimal places.",
        type=float,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "-l",
        "--layer-height",
        help="Slicing layer height",
        type=float,
        default=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--outer-speed",
        type=float,
        default=argparse.SUPPRESS,
        help="Outer perimeter print speed (mm/s).",
    )
    parser.add_argument(
        "--inner-speed",
        type=float,
        default=argparse.SUPPRESS,
        help="Inner perimeter print speed (mm/s).",
    )
    parser.add_argument(
        "--infill-speed",
        type=float,
        default=argparse.SUPPRESS,
        help="Infill print speed (mm/s).",
    )
    parser.add_argument(
        "--max-corner-speed",
        type=float,
        default=argparse.SUPPRESS,
        help="Max speed at a 90° corner (mm/s).",
    )
    parser.add_argument(
        "--max-accel",
        type=float,
        default=argparse.SUPPRESS,
        help="Max print acceleration for feed planning (mm/s²).",
    )
    parser.add_argument(
        "--max-jerk",
        type=float,
        default=argparse.SUPPRESS,
        help="Max corner speed change Δv for feed planning (mm/s).",
    )
    parser.add_argument(
        "--min-corner-angle",
        type=float,
        default=argparse.SUPPRESS,
        help="Min turn angle (deg) before accel/jerk caps apply.",
    )
    parser.add_argument(
        "--html-preview",
        metavar="PATH",
        help="Also write an HTML toolpath preview to this path.",
    )
    return parser


def _cli_overrides_to_model_settings(args) -> dict:
    """Map explicitly passed CLI flags to Model-native settings."""
    raw = {k: v for k, v in vars(args).items() if k in _CLI_SETTING_DESTS}
    if not raw:
        return {}
    # Dest names match YAML keys (num_perimeters, layer_height, outer_speed, …).
    return canonicalize_config(raw)


def configure_model(args) -> Model:
    """Build a Model with CLI historical defaults → YAML layers → CLI overrides."""
    model = Model()
    apply_settings_to_model(model, _CLI_HISTORICAL_DEFAULTS)
    config_arg = getattr(args, "config", None)
    if config_arg:
        paths = parse_config_paths(config_arg)
        apply_settings_to_model(model, load_layered_config(paths))
    apply_settings_to_model(model, _cli_overrides_to_model_settings(args))
    return model


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
            model = configure_model(args)
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
                if model.perimeters_only:
                    layers[layer_num].infill = []
                    continue
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
                max_accel=model.max_accel,
                planning_limits={
                    # User-facing speeds are mm/s (G-code F is mm/min)
                    "maxSpeed": model.outer_perimeter_speed / 60.0,
                    "outerSpeed": model.outer_perimeter_speed / 60.0,
                    "innerSpeed": model.inner_perimeter_speed / 60.0,
                    "infillSpeed": model.infill_speed / 60.0,
                    "maxAccel": model.max_accel,
                    "maxJerk": model.max_jerk,
                    "minAngleDeg": model.min_corner_angle,
                    "maxCornerSpeed": model.max_corner_speed / 60.0,
                },
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
    from pyslicer.config import ConfigError

    parser = build_arg_parser()
    args = parser.parse_args(argv)
    try:
        run(args)
    except ConfigError as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    main()
