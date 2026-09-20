"""Tests for FDM bead cross-section approximation."""

import math

import pytest

from pyslicer.mesh.model import Model
from pyslicer.preview.bead_profile import (
    STACK_OVERLAP,
    bead_cross_section,
    extrusion_width,
)


def test_extrusion_width_matches_model_offset():
    m = Model()
    m.nozzle_diameter = 0.5
    m.layerHeight = 0.25
    assert math.isclose(extrusion_width(0.5, 0.25), m.offset())


def test_bead_is_stadium_not_circle_when_nozzle_gt_layer():
    bead = bead_cross_section(0.5, 0.25)
    assert math.isclose(bead["height"], 0.25 * STACK_OVERLAP)
    assert bead["width"] > bead["height"]
    assert math.isclose(bead["sideRadius"], bead["height"] / 2)
    assert math.isclose(bead["flatWidth"], bead["width"] - bead["height"])
    assert bead["onBed"] is False


def test_bead_height_overlaps_layer_pitch():
    """Stacked layers must overlap — height strictly exceeds layer pitch."""
    for nozzle, layer in [(0.4, 0.4), (0.5, 0.2), (0.4, 0.2), (0.6, 0.3)]:
        bead = bead_cross_section(nozzle, layer)
        assert bead["height"] > layer
        assert bead["height"] >= layer * STACK_OVERLAP - 1e-12


def test_near_circular_bead_still_stacks_without_air_gap():
    """
    When nozzle ≈ layer height the stadium is nearly round. Circles that only
    touch look gappy; overlap must keep a positive neck between layer centers.
    """
    layer = 0.4
    bead = bead_cross_section(0.4, layer)
    assert bead["height"] > layer
    r = bead["height"] / 2.0
    # Chord width at the mid-plane between adjacent layer centers:
    neck = 2.0 * math.sqrt(r * r - (layer / 2.0) ** 2)
    assert neck > layer * 0.5  # visible seal, not a vanishing point
    assert bead["flatWidth"] >= 0.0


def test_bead_degenerates_to_circle_when_width_clamped_to_height():
    # Large layer vs small nozzle: width clamps up to (inflated) height.
    bead = bead_cross_section(0.2, 0.4)
    assert math.isclose(bead["width"], bead["height"])
    assert math.isclose(bead["flatWidth"], 0.0)


def test_bead_on_bed_is_slightly_wider():
    normal = bead_cross_section(0.5, 0.2, on_bed=False)
    bed = bead_cross_section(0.5, 0.2, on_bed=True)
    assert bed["width"] > normal["width"]
    assert bed["onBed"] is True


def test_bead_rejects_non_positive():
    with pytest.raises(ValueError):
        bead_cross_section(0.5, 0.0)
    with pytest.raises(ValueError):
        extrusion_width(-1, 0.2)
