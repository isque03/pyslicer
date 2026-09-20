"""Deposited FDM bead cross-section (stadium / capsule approximation)."""

from __future__ import annotations

import math

# Same extrusion-width heuristic as ``Model.offset()``.
_WIDTH_FACTOR = 1.0 - (math.pi / 4.0)


def extrusion_width(nozzle_diameter: float, layer_height: float) -> float:
    """Expected bead width (mm) from nozzle diameter and layer height."""
    if nozzle_diameter <= 0 or layer_height <= 0:
        raise ValueError(
            f"nozzle_diameter and layer_height must be > 0 "
            f"(got nozzle={nozzle_diameter}, layer={layer_height})"
        )
    width = float(nozzle_diameter) + float(layer_height) * _WIDTH_FACTOR
    # Stadium needs width >= height (else it becomes a circle of diameter=height).
    # Also keep width at least the nozzle so side bulges stay visible.
    return max(width, float(layer_height), float(nozzle_diameter))


def bead_cross_section(
    nozzle_diameter: float,
    layer_height: float,
    *,
    on_bed: bool = False,
) -> dict:
    """
    Approximate how filament looks after extrusion.

    Real beads are not round tubes when layer height < nozzle diameter. The
    nozzle leaves a **flat top**, the bed/prior layer leaves a **flat bottom**,
    and plastic bulges into **rounded sides** (a stadium / capsule section).

    Returns mm dimensions for a 2D profile in the plane perpendicular to travel:

    - ``height`` — vertical extent (≈ layer height)
    - ``width`` — full horizontal extent including side bulges
    - ``sideRadius`` — radius of each side semicircle (``height / 2``)
    - ``flatWidth`` — length of the flat top/bottom between the side arcs
    - ``onBed`` — if True, width is increased slightly (first-layer squash)

    Raises ValueError for non-positive inputs.
    """
    height = float(layer_height)
    width = extrusion_width(nozzle_diameter, layer_height)
    if on_bed:
        # First layer is typically pressed flatter/wider onto the build plate
        width *= 1.06
    side_radius = height / 2.0
    flat_width = max(width - height, 0.0)
    return {
        "height": height,
        "width": width,
        "sideRadius": side_radius,
        "flatWidth": flat_width,
        "onBed": bool(on_bed),
    }
