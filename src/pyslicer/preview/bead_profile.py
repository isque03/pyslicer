"""Deposited FDM bead cross-section (stadium / capsule approximation)."""

from __future__ import annotations

import math

# Same extrusion-width heuristic as ``Model.offset()``.
_WIDTH_FACTOR = 1.0 - (math.pi / 4.0)

# Scale the whole stadium slightly so stacked layers overlap in the preview.
# Applied uniformly so nozzle>layer flattening (flat top/bottom) is preserved.
STACK_OVERLAP = 1.25


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

    - ``height`` — vertical extent (layer height × stack overlap)
    - ``width`` — full horizontal extent including side bulges
    - ``sideRadius`` — radius of each side semicircle (``height / 2``)
    - ``flatWidth`` — length of the flat top/bottom between the side arcs
    - ``onBed`` — if True, width is increased slightly (first-layer squash)

    Raises ValueError for non-positive inputs.
    """
    layer_h = float(layer_height)
    # Physical stadium from true layer height (this is the flattening).
    height = layer_h
    width = extrusion_width(nozzle_diameter, layer_h)
    if on_bed:
        # First layer is typically pressed flatter/wider onto the build plate
        width *= 1.06
    width = max(width, height)

    # Uniform preview scale: seal stacks without erasing width/height flattening.
    height *= STACK_OVERLAP
    width *= STACK_OVERLAP

    side_radius = height / 2.0
    flat_width = max(width - height, 0.0)
    return {
        "height": height,
        "width": width,
        "sideRadius": side_radius,
        "flatWidth": flat_width,
        "onBed": bool(on_bed),
    }
