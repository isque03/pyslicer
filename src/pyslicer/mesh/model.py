"""Slicer model state and print settings."""

import math
import sys

import numpy as np


class Model:
    def __init__(self):
        self.facets = []
        self.layers = []
        self.facetsAtLayer = {}
        self.name = "new model"
        self.triangles = 0
        self.endLoop = 0
        self.zmin = sys.float_info.max
        self.zmax = -sys.float_info.max
        self.nozzle_diameter = 0.5
        self.layerHeight = 1.0 * self.nozzle_diameter
        self.retract_amount = 3.0
        self.retract_speed = 6200
        self.unretract_speed = 3400
        self.default_print_speed = 4200  # mm/min (G-code F); 70 mm/s
        self.outer_perimeter_speed = 4200
        self.inner_perimeter_speed = 4200
        self.infill_speed = 4200
        self.default_travel_speed = 8000
        self.default_z_speed = 2400
        # Corner planning (stored as G-code F units mm/min; CLI/UI use mm/s)
        self.max_corner_speed = 300  # mm/min ≈ 5 mm/s at 90° (SCV)
        self.max_accel = 1000.0  # mm/s²
        self.max_jerk = 20.0  # mm/s corner Δv proxy
        self.min_corner_angle = 20.0  # degrees; gentler turns uncapped by accel/jerk
        self.filament_diameter = 1.75
        self.number_perimeters = 2
        self.processing_threads = 4
        self.infill_density = 0.25
        self.infill_angle = 45.0
        self.width_over_height = 1.9
        self.perimeters_only = True
        self.append_perimeters = True
        self.perimeter_overlap_percent = 1.0
        self.minimum_retract_travel = 3.0
        self.print_temperature = 180.0
        self.simplification_factor = 4.0
        self.min_contour_area = 50000.00
        self.min_extrude = 1.25 * self.nozzle_diameter
        self.contours = []
        self.vertices = []
        self.edges = {}
        self.infillAtLayer = {}
        self.edgeList = []
        # Numpy mesh cache: (N, 3, 3) float64 vertex positions
        self.facet_vertices = None
        self.facet_zmin = None
        self.facet_zmax = None

    def rebuild_facet_arrays(self):
        """Cache facet vertices as a numpy array for fast Z filtering."""
        n = len(self.facets)
        if n == 0:
            self.facet_vertices = np.empty((0, 3, 3), dtype=np.float64)
            self.facet_zmin = np.empty((0,), dtype=np.float64)
            self.facet_zmax = np.empty((0,), dtype=np.float64)
            return
        arr = np.empty((n, 3, 3), dtype=np.float64)
        for i, facet in enumerate(self.facets):
            for j, v in enumerate(facet.verticies):
                arr[i, j, 0] = v.x
                arr[i, j, 1] = v.y
                arr[i, j, 2] = v.z
        self.facet_vertices = arr
        zvals = arr[:, :, 2]
        self.facet_zmin = zvals.min(axis=1)
        self.facet_zmax = zvals.max(axis=1)

    def offset(self):
        """Extrusion width from nozzle diameter and layer height."""
        return self.nozzle_diameter + (self.layerHeight * (1 - (math.pi / 4.0)))

    def volume_extruded(self, segment):
        length = segment.magnitude()
        return length * self.layerHeight * self.nozzle_diameter

    def pi_r_squared(self):
        return math.pi * math.pow((self.filament_diameter / 2.0), 2.0)

    def writeGCode(self, filename):
        from pyslicer.gcode.writer import write_gcode

        write_gcode(self, filename)
