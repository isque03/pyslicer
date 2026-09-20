"""Tests for FDM bead cross-section approximation."""

import math

import pytest

from pyslicer.mesh.model import Model
from pyslicer.preview.bead_profile import bead_cross_section, extrusion_width


def test_extrusion_width_matches_model_offset():
    m = Model()
    m.nozzle_diameter = 0.5
    m.layerHeight = 0.25
    assert math.isclose(extrusion_width(0.5, 0.25), m.offset())


def test_bead_is_stadium_not_circle_when_nozzle_gt_layer():
    bead = bead_cross_section(0.5, 0.25)
    assert math.isclose(bead["height"], 0.25)
    assert bead["width"] > bead["height"]
    assert math.isclose(bead["sideRadius"], 0.125)
    assert math.isclose(bead["flatWidth"], bead["width"] - bead["height"])
    assert bead["onBed"] is False


def test_bead_degenerates_to_circle_when_width_equals_height():
    # Force equal via very large layer relative to nozzle still clamps width>=height
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
