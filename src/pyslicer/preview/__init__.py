"""Browser preview helpers for G-code output."""

from pyslicer.preview.bead_profile import bead_cross_section, extrusion_width
from pyslicer.preview.gcode_html import (
    layers_to_toolpaths_3d,
    parse_gcode_layers,
    parse_toolpath_moves,
    render_gcode_html,
    write_gcode_preview,
)
from pyslicer.preview.gcode_parse import build_timeline, parse_gcode
from pyslicer.preview.render_anim import write_simulation_gif

__all__ = [
    "bead_cross_section",
    "build_timeline",
    "extrusion_width",
    "layers_to_toolpaths_3d",
    "parse_gcode",
    "parse_gcode_layers",
    "parse_toolpath_moves",
    "render_gcode_html",
    "write_gcode_preview",
    "write_simulation_gif",
]
