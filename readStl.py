#!/usr/bin/env pypy
import sys

from datetime import datetime

from clipper import *
from timer import Timer

# import resource
from multiprocessing.dummy import Pool
from functools import partial
import math
import traceback
import threading
import logging
from operator import itemgetter, attrgetter
from copy import copy
from struct import unpack
import argparse


def dumpstacks(signal, frame):
    id2name = dict([(th.ident, th.name) for th in threading.enumerate()])
    code = []
    for threadId, stack in sys._current_frames().items():
        code.append("\n# Thread: %s(%d)" %
                    (id2name.get(threadId, ""), threadId))
        for filename, lineno, name, line in traceback.extract_stack(stack):
            code.append('File: "%s", line %d, in %s' %
                        (filename, lineno, name))
            if line:
                code.append("  %s" % (line.strip()))
    print "\n".join(code)


# signal.signal(signal.SIGQUIT, dumpstacks)

CLIPPER_SCALE_FACTOR = 10000000


class Intersection:
    def __init__(self):
        self.uparam = -1.0
        self.contour = None
        self.segmentIndex = -1
        self.row = -1
        self.ray = None
        self.used = False

    def next(self):
        return None

    def prev(self):
        return None

    def point(self):
        return self.ray.pointOnLine(self.uparam)

    def __repr__(self):
        return '{}: u: {} segmentIndex: {} row: {} contour: {}'.format(self.__class__.__name__,
                                                                       self.uparam,
                                                                       self.segmentIndex,
                                                                       self.row,
                                                                       self.contour)


'''
A Layer contains the results of slicing a model at a given layer. 
This includes contours, zlevel, roof, overhang and bridge information.
'''


class Layer:
    def __init__(self):
        self.contours = None
        self.z = None
        self.infill = None
        self.roofs = []
        self.overhang = []
        self.normalInfill = []

    def to_svg(self):

        with open('layer-{}.svg'.format(self.z), "w") as f:
            f.write('<?xml version="1.0" standalone="no"?>\n')
            f.write(
                '<svg width="1000mm" height="1000mm" viewBox="-50 -50 100 100" version="1.1" xmlns="http://www.w3.org/2000/svg">\n')
            # f.write('<path d="')
            # for i in range(len(polys)):
            # f.write(' M{} {}'.format(polys[i][0].x, polys[i][0].y))
            # f.write(' L{} {}'.format(polys[i][1].x, polys[i][1].y))
            #f.write('" stroke="red" fill="transparent" stroke-width=".01"/>\n')
            for x in range(len(self.perimeters)):
                import random

                color = lambda: random.randint(0, 255)
                r = color()
                g = color()
                b = color()
                for contour in self.perimeters[x]:
                    f.write('<path d="')
                    for path in contour.segments:
                        f.write(' M{} {}'.format(path.verticies[0].x, path.verticies[0].y))
                        f.write(' L{} {}'.format(path.verticies[1].x, path.verticies[1].y))

                    f.write(
                        '" stroke="#{:02x}{:02X}{:02X}" fill="transparent" stroke-width=".1"/>\n'.format(r, g, b))
                    f.write(
                        '<text x="{0}" y="{1}" font-size="0.5">P:{2} {0}x,{1}y</text>'.format(
                            contour.segments[0].verticies[0].x,
                            contour.segments[0].verticies[0].y,
                            x))
            f.write('</svg>')


class Model:
    def __init__(self):
        self.facets = []
        self.layers = []
        self.facetsAtLayer = {}
        self.name = 'new model'
        self.triangles = 0
        self.endLoop = 0
        self.zmin = sys.maxint
        self.zmax = -sys.maxint - 1
        self.nozzle_diameter = 0.5
        self.layerHeight = 1.0 * self.nozzle_diameter
        self.retract_amount = 3.0
        self.retract_speed = 6200
        self.unretract_speed = 3400
        self.default_print_speed = 4200
        self.default_travel_speed = 8000
        self.default_z_speed = 2400
        self.filament_diameter = 1.75
        self.number_perimeters = 2
        self.processing_threads = 4
        self.infill_density = 0.25
        self.infill_angle=45.0
        self.width_over_height = 1.9
        self.perimeters_only = True
        self.append_perimeters = True
        self.perimeter_overlap_percent = 1.0
        self.minimum_retract_travel=3.0
        # Print Temp in degrees C
        self.print_temperature = 180.0
        # nozzle diameter divided by this equals the tolerance for curve simplification
        self.simplification_factor = 4.0
        # remove contours with a winding area smaller than this
        self.min_contour_area = 50000.00
        # shortest distance to extrude when infilling
        self.min_extrude = 1.25 * self.nozzle_diameter
        self.contours = []
        # experiement
        self.vertices = []
        self.edges = {}
        self.infillAtLayer = {}
        self.edgeList = []


    def offset(self):
        """
        Actual width of extruded material
        :rtype : the calculated extrusion width based on the current layer height and feed rate
        """
        width = self.nozzle_diameter + (self.layerHeight * (1 - (math.pi / 4.0)))
        # print 'Calculated width {}'.format(width)
        return width

    def countOpenFacets(self):
        # 0,1 1,2 2,0
        openFaces = 0
        for facet in self.facets:
            f0 = self.facetForPoints(facet.verticies[0], facet.verticies[1], facet)
            if f0 is None: openFaces += 1
            f1 = self.facetForPoints(facet.verticies[1], facet.verticies[2], facet)
            if f1 is None: openFaces += 1
            f2 = self.facetForPoints(facet.verticies[2], facet.verticies[0], facet)
            if f2 is None: openFaces += 1
        return openFaces

    def facetForPoints(self, p1, p2, otherFacet):
        result = None
        for facet in model.facets:
            if facet == otherFacet:
                continue
            if facet.verticies[0] == p1 and facet.verticies[1] == p2:
                return facet
            elif facet.verticies[1] == p1 and facet.verticies[0] == p2:
                return facet
            elif facet.verticies[1] == p1 and facet.verticies[2] == p2:
                return facet
            elif facet.verticies[2] == p1 and facet.verticies[1] == p2:
                return facet
            elif facet.verticies[2] == p1 and facet.verticies[0] == p2:
                return facet
            elif facet.verticies[0] == p1 and facet.verticies[2] == p2:
                return facet

        return result

    def writeSegmentGCode(self, f, idx, segment, e, extrude_amount):
        if idx == 0:
            e = self.retract(e, f)
            f.write(';; travel move \n')
            f.write('G1 X{0:.6f} Y{1:.6f}\n'.format(
                segment.verticies[0].x, segment.verticies[0].y))
            e = self.extrude(e, f)
            # f.write(';; idx: {} \n'.format(idx))
            e += extrude_amount
            f.write('G1 F{3:.6f} X{0:.6f} Y{1:.6f} E{2:.6f}\n'.format(
                segment.verticies[1].x, segment.verticies[1].y, e, self.default_print_speed))
        else:
            # f.write(';; idx: {}\n '.format(idx))
            e += extrude_amount
            f.write('G1 X{0:.6f} Y{1:.6f} E{2:.6f}\n'.format(
                segment.verticies[1].x, segment.verticies[1].y, e))
        return e

    def retract(self, e, f):
        if model.retract_amount > 0.0:
            e -= model.retract_amount
            f.write('G1 F{0:.6f} E{1:.6f}\n'.format(model.retract_speed, e))
        return e

    def extrude(self, e, f):
        if model.retract_amount > 0.0:
            e += model.retract_amount
            f.write('G1 F{0:.6f} E{1:.6f}\n'.format(model.unretract_speed, e))
        return e

    def volume_extruded(self, segment):
        length = segment.magnitude()
        # calculate volume of material needed for segment
        volume = length * self.layerHeight * self.nozzle_diameter
        return volume

    def pi_r_squared(self):
        pirsquared = math.pi * math.pow((self.filament_diameter / 2.0), 2.0)
        return pirsquared

    def writeGCode(self, filename):



        logger.info('Writing {0}'.format(filename))
        f = open(filename, 'w')
        f.write('M109 S{0} ; Heat up to {0}C\n'.format(self.print_temperature))
        f.write('G90       ; Use absolute coordinates\n')
        f.write('G21       ; Set units to millimeters\n')
        f.write('M106 S0   ; Fan Off\n')
        f.write('G28       ; Home all axes\n')
        f.write('G92 E0    ; Zero extruder\n')
        f.write('M82       ; Use absolute distances for extrusion\n')
        f.write('G92 E0    ; Zero extruder\n')
        f.write('G1 F200 E8  ; prime extruder\n')
        f.write('G92 E0    ; Zero extruder\n')
        f.write('M117 Printing.\n')
        f.write('\n')
        e = 0.0000
        lastZ = None
        logger.debug(
            'size of infill at layer is {0} model is {1}'.format(len(self.infillAtLayer), len(self.contours)))
        # we need this later for extruder calculations. do it once here.
        pirsquared = self.pi_r_squared()
        for layer in self.layers:
            f.write(';; New Layer Z: {0} \n'.format(layer.z))
            contourIndex = 0
            while True:
                contoursWritten = 0
                for x in xrange(len(layer.perimeters) -1, -1, -1):
                    #for contour in layer.perimeters[x]:
                    if len(layer.perimeters[x]) < contourIndex+1:
                       continue
                    f.write(';; Perimeter {}\n'.format(x))
                    contoursWritten += 1
                    contour = layer.perimeters[x][contourIndex]
                    if (len(contour.segments) < 1):
                        continue
                    if lastZ is None or (lastZ != contour.zlevel):
                        e = self.retract(e, f)
                        f.write('G1 F{1:.6f} Z{0:.6f}\n'.format(contour.zlevel, self.default_z_speed))
                        e = self.extrude(e, f)
                        # print '===== {0} ====='.format(contour.zlevel)
                    lastZ = contour.zlevel
                    f.write(';; Contour {} Area: {}\n'.format(contourIndex, abs(contour.winding_area())))
                    f.write('G1 F{0:.6f}\n'.format(self.default_print_speed))
                    for idx, segment in enumerate(contour.segments):
                        volume = self.volume_extruded(segment)
                        # get the amount of filament needed to extrude that volume of material
                        extrude_amount = (volume / pirsquared)
                        # f.write(';; length {} volume {} pirsquared {}\n'.format(length, volume, pirsquared))
                        e = self.writeSegmentGCode(f, idx, segment, e, extrude_amount)
                contourIndex += 1
                # check if we are done
                if contoursWritten == 0:
                    break

            # write infill at this z level
            try:
                infill = self.infillAtLayer[lastZ]
                f.write(';; Infill\n')
                prevSegment = None
                for idx, segment in enumerate(infill):
                    # Don't need this travel move step all the time
                    travelMove = False
                    if idx == 0: travelMove = True
                    if (prevSegment is not None and (prevSegment.verticies[1] != segment.verticies[0])):
                        travelMove = True
                    if travelMove:
                        # Calc travel distance
                        travel_distance = 0.0
                        if prevSegment:
                            l = Line.withVerticies(prevSegment.verticies[1], segment.verticies[0])
                            travel_distance = l.magnitude()
                        if travel_distance >= self.minimum_retract_travel:
                            # or exceed some distance
                            f.write(';; retract \n')
                            e = self.retract(e, f)
                        f.write(';;infill travel move. distance: {0:.6f} \n'.format(travel_distance))
                        f.write('G1 F{2:.6f} X{0:.6f} Y{1:.6f}\n'.format(
                            segment.verticies[0].x, segment.verticies[0].y, self.default_travel_speed))
                        if travel_distance >= self.minimum_retract_travel:
                            e = self.extrude(e, f)
                        # Back to slower print speed
                        f.write('G1 F{0:.6f}\n'.format(self.default_print_speed))

                    # get the amount of filament needed to extrude that volume of material
                    volume = self.volume_extruded(segment)
                    e += (volume / pirsquared)
                    f.write('G1 X{0:.6f} Y{1:.6f} E{2:.6f}\n'.format(
                        segment.verticies[1].x, segment.verticies[1].y, e))
                    prevSegment = segment

            except KeyError:
                logger.error('GCODE output: no infill at layer {0}'.format(lastZ))



                # We can't just close the contour since it might not really connect back to the start
                # typically if we get a bad (non-manifold) stl
                # close contour
                # print 'closing contour.'
                # f.write(';; closing contour')
                # f.write('G1 X{0:.6f} Y{1:.6f} E{2:.6f}\n'.format(
                # contour.segments[0].verticies[1].x,
                # contour.segments[0].verticies[1].y, e))
            if False:
                f.write(';; Roofs\n')
                for contour in layer.roofs:
                    for idx, segment in enumerate(contour.segments):
                        extrude_amount = (volume / pirsquared)
                        self.writeSegmentGCode(f, idx, segment, e, extrude_amount)
                        e += 1.0
                f.write(';; Overhang\n')
                for contour in layer.overhang:
                    for idx, segment in enumerate(contour.segments):
                        extrude_amount = (volume / pirsquared)
                        self.writeSegmentGCode(f, idx, segment, e, extrude_amount)
                        e += 1.0

        # END
        f.write('M107    ; Fan off\n')
        f.write('M104 S0 ; Heat off {0}C\n')
        f.write('M140 S0 ; Bed heat off\n')
        f.write('G91     ; Relative\n')
        f.write('G1 E-1 F400 ; Retract\n')
        f.write('G1 Z+1.0 E-5 X-20 Y-20 F9000\n')
        f.write('G28 X0 Y0\n')
        f.write('M84     ; Motors off\n')
        f.write('G90     ; Absolute\n')
        f.close()


class Facet:
    def __init__(self):
        self.verticies = []
        self.edges = []

    @classmethod
    def withVerticies(clazz, v1, v2, v3):
        self = clazz()
        self.verticies.append(v1)
        self.verticies.append(v2)
        self.verticies.append(v3)
        return self

    def isCoplanar(self, z):
        coplanar = True
        for vertex in self.verticies:
            if not vertex.isequal2(vertex.z, z):
                coplanar = False
            if not coplanar:
                break
        return coplanar

    def inZPlane(self):
        return self.isCoplanar(self.verticies[0].z)


class Vertex:
    '''represents a vertex or a point'''

    def __init__(self, x, y, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return "%s [%.15f,%.15f,%.15f] " % (id(self), self.x, self.y, self.z)

    def isequal2(self, a, b):
        return Math.floatEq(a, b)

    def clipperX(self):
        return int(self.x * CLIPPER_SCALE_FACTOR)

    def clipperY(self):
        return int(self.y * CLIPPER_SCALE_FACTOR)

    def __eq__(self, other):
        return (self.isequal2(self.x, other.x)
                and self.isequal2(self.y, other.y)
                and self.isequal2(self.z, other.z))

    def __key(self):
        return (self.x, self.y, self.z)

    def __hash__(self):
        # return hash(self.x) ^ hash(self.y) ^ hash(self.z)
        return hash(self.__key())

        # precision = Decimal('1.00000')
        # return (isinstance(other, self.__class__)
        # and self.x.quantize(precision)
        # .compare(other.x.quantize(precision)) == 0
        # and self.y.quantize(precision)
        #        .compare(other.y.quantize(precision)) == 0
        #        and self.z.quantize(precision)
        #        .compare(other.z.quantize(precision)) == 0)


class Math:
    @classmethod
    def floatEq(self, a, b):
        eps1 = sys.float_info.epsilon
        eps2 = sys.float_info.epsilon
        eps1 = 0.0000001
        eps2 = eps1
        # if abs(a - b) <= eps1 * (abs(a) + abs(b)) + eps2:
        if abs(a - b) <= eps1:
            return True
        else:
            return False

    '''
    removes duplicate intersections
    '''''

    @classmethod
    def remove_duplicate_intersections(self, mylist):
        if mylist:
            mylist.sort(key=lambda x: x.uparam)
            last = mylist[-1]
            for i in range(len(mylist) - 2, -1, -1):
                if Math.floatEq(last.uparam, mylist[i].uparam):
                    del mylist[i]
                else:
                    last = mylist[i]


class Line:
    '''represents a line segment of two points'''

    def __init__(self):
        self.verticies = []
        self.facets = []  # referencing facets

    @classmethod
    def withVerticies(clz, va, vb):
        self = clz()
        self.verticies.append(va)
        self.verticies.append(vb)
        return self

    def key(self):
        formatOne = (
            "[%.6f,%.6f,%.6f] [%.6f,%.6f,%.6f]" % (self.verticies[0].x,
                                                   self.verticies[
                                                       0].y,
                                                   self.verticies[
                                                       0].z,
                                                   self.verticies[
                                                       1].x,
                                                   self.verticies[
                                                       1].y,
                                                   self.verticies[
                                                       1].z
            ))
        formatTwo = (
            "[%.6f,%.6f,%.6f] [%.6f,%.6f,%.6f]" % (self.verticies[1].x,
                                                   self.verticies[
                                                       1].y,
                                                   self.verticies[
                                                       1].z,
                                                   self.verticies[
                                                       0].x,
                                                   self.verticies[
                                                       0].y,
                                                   self.verticies[
                                                       0].z
            ))
        if self.verticies[0].x > self.verticies[1].x:
            return formatOne
        elif self.verticies[0].x < self.verticies[1].x:
            return formatTwo
        elif self.verticies[0].y > self.verticies[1].y:
            return formatOne
        else:
            return formatTwo

    def __str__(self):
        return (
            "%s [%.15f,%.15f,%.15f] [%.15f,%.15f,%.15f]" %
            (id(self), self.verticies[0].x,
             self.verticies[
                 0].y,
             self.verticies[
                 0].z,
             self.verticies[
                 1].x,
             self.verticies[
                 1].y,
             self.verticies[
                 1].z
            ))

    def __eq__(self, other):
        return ((self.verticies[0] ==
                 other.verticies[0]
                 and self.verticies[1] ==
                 other.verticies[1]) or
                (self.verticies[0] ==
                 other.verticies[1]
                 and self.verticies[1] ==
                 other.verticies[0]))

    def __ne__(self, other):
        return not self.__eq__(other)

    def reverse(self):
        if (self.verticies is None):
            return
        first = self.verticies[0]
        self.verticies[0] = self.verticies[1]
        self.verticies[1] = first

    '''
    Calculatea the magnitude of this line (vector)
    '''

    def magnitude(self):
        if (len(self.verticies) != 2):
            raise Exception(
                'Cannot calculate magnitude of segment not having exactly two vertices')
        return math.sqrt(self.magnitudeSquared())

    '''
    Calculates the magnitude squared of this line.
    '''

    def magnitudeSquared(self):
        return ((self.verticies[0].x - self.verticies[1].x) ** 2 +
                (self.verticies[0].y - self.verticies[1].y) ** 2 +
                (self.verticies[0].z - self.verticies[1].z) ** 2)

    '''
    Calculate the perpendicular distance from a given point to this line.
    '''

    def distance(self, vertex):
        magsqr = self.magnitudeSquared()
        if magsqr <= 0.0:
            return 0.0
        u = (((vertex.x - self.verticies[0].x) * (self.verticies[1].x - self.verticies[0].x)
              + (vertex.y - self.verticies[0].y) * (self.verticies[1].y - self.verticies[0].y))
             / magsqr)
        if (u < 0.0 or u > 1.0):
            return -1
            # outside of the line segment
        x = self.verticies[0].x + u * \
                                  (self.verticies[1].x - self.verticies[0].x)
        y = self.verticies[0].y + u * \
                                  (self.verticies[1].y - self.verticies[0].y)
        point = Vertex(x, y, vertex.z)
        tmpLine = Line()
        tmpLine.verticies.append(vertex)
        tmpLine.verticies.append(point)
        return tmpLine.magnitude()

    '''
    Returns a point at a given u parameter on this line where 0 <= u <= 1
    '''

    def pointOnLine(self, u):
        # intersection point is P(s) = P0 + s(P1-P0)
        ps = addVectors(self.verticies[0], multiplyByScalar(
            subtractVectors(self.verticies[1], self.verticies[0]), u))
        return ps

    '''
    Returns the 2D intersection point between the two given segment
    or raise exceptions if the segments do not intersect.
    Based off the explaination by Paul Bourke.
    '''

    def intersect2D(self, line):
        ua = self.intersect2DU(line)
        x1 = self.verticies[0].x
        x2 = self.verticies[1].x
        y1 = self.verticies[0].y
        y2 = self.verticies[1].y

        x = x1 + ua * (x2 - x1)
        y = y1 + ua * (y2 - y1)
        return Vertex(x, y, self.verticies[0].z)

    '''
    Returns the u parameters of the 2D intersections of the given segment
    and this line for any u value. A u value not between 0 and 1
    does not fall on the segment, but falls on the ray created by extending
    the segment.
    '''

    def intersect2DUExtend(self, line):
        '''

        x = x1 + ua (x2 - x1)
        y = y1 + ua (y2 - y1)

        The denominators for the equations for ua and ub are the same.
        If the denominator for the equations for ua and ub is 0 then the two
        lines are parallel.
        If the denominator and numerator for the equations for ua and ub are
        0 then the two lines are coincident.
        The equations apply to lines, if the intersection of line segments 
        is required then it is only necessary to test if ua and ub lie 
        between 0 and 1. Whichever one lies within that range then the 
        corresponding line segment contains the intersection point. If both
        lie within the range of 0 to 1 then the intersection point is within
        both line segments.

        '''

        x4 = line.verticies[1].x
        x3 = line.verticies[0].x
        x2 = self.verticies[1].x
        x1 = self.verticies[0].x

        y4 = line.verticies[1].y
        y3 = line.verticies[0].y
        y2 = self.verticies[1].y
        y1 = self.verticies[0].y

        denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1)
        numeratora = (x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)
        numeratorb = (x2 - x1) * (y1 - y3) - (y2 - y1) * (x1 - x3)

        # Test for coincident lines
        if (Math.floatEq(denominator, 0.0)
            and Math.floatEq(numeratora, 0.0)
            and Math.floatEq(numeratorb, 0.0)):
            raise CoincidentLines('Lines are coincient')
            # return Vertex((x1+x2)/2,(y1+y2)/2,self.verticies[0].z)
            # return 0.5

        # Test for parallel lines
        if (Math.floatEq(denominator, 0.0)):
            raise ParallelLines('Lines are parallel and do not intersect')

        # Test that intersection point is within both segments (0<u<1)
        ua = numeratora / denominator
        ub = numeratorb / denominator
        return ua, ub


    '''
    Returns the u parameter of the 2D intersection of the given segment
    and this line for 0 < u < 1
    '''

    def intersect2DU(self, line):
        ua, ub = self.intersect2DUExtend(line)
        if (ua < 0 or ua > 1 or ub < 0 or ub > 1):
            raise NonIntersectingLines('Segments do not intersect with 0<u<1')
        return ua

    def normal(self):
        dx = self.verticies[1].x - self.verticies[0].x
        dy = self.verticies[1].y - self.verticies[0].y
        # (-dy, dx) and (dy, -dx).
        return Line.withVerticies(Vertex(-dy, dx, self.verticies[0].z),
                                  Vertex(dy, -dx, self.verticies[0].z))

    def unitVector(self):
        '''
        Return a new line which is the unit vector of this line
        '''
        mag = self.magnitude()
        v1 = Vertex(self.verticies[0].x / mag, self.verticies[
            0].y / mag, self.verticies[0].z / mag)
        v2 = Vertex(self.verticies[1].x / mag, self.verticies[
            1].y / mag, self.verticies[1].z / mag)
        return Line.withVerticies(v1, v2)

    def scale2D(self, scaleFactor):
        self.verticies[0] = multiplyByScalar2D(self.verticies[0], scaleFactor)
        self.verticies[1] = multiplyByScalar2D(self.verticies[1], scaleFactor)

    def copy(self):
        v1 = Vertex(
            self.verticies[0].x, self.verticies[0].y, self.verticies[0].z)
        v2 = Vertex(
            self.verticies[1].x, self.verticies[1].y, self.verticies[1].z)
        return Line.withVerticies(v1, v2)

    def translate(self, vector):
        '''
        Translate this line by the vector represented by
        the given vertex.
        '''
        self.verticies[0].x += vector.x
        self.verticies[0].y += vector.y
        self.verticies[0].z += vector.z
        self.verticies[1].x += vector.x
        self.verticies[1].y += vector.y
        self.verticies[1].z += vector.z
        return self

    '''
    Offset this line by a given distance.
    A distance with a nagative value will offset in the opposite direction
    '''

    def offset(self, distance):
        # get normal
        n = self.normal()
        # pick the inside
        p = n.pointOnLine(0.5)
        # get unit vector in the direction we want
        startPoint = 0
        if distance < 0:
            startPoint = 1
        t = Line.withVerticies(p, n.verticies[startPoint]).unitVector()
        # scale it to the inset distance
        t.scale2D(distance)
        translationVector = Vertex(t.verticies[1].x, t.verticies[1].y)
        return self.translate(translationVector)

    '''
    Rotate this line by a given angle in degrees about a given point
    '''

    def rotate(self, angle, point):
        angle_rad = math.radians(angle)
        rotationMatrix = [[math.cos(angle_rad), -math.sin(angle_rad)],
                          [math.sin(angle_rad), math.cos(angle_rad)]]
        rotated = []
        # oR.x = oP.x + (o.x - oP.x) * cos(theta) - (o.y - oP.y) * sin(theta)
        class expando(object): pass

        oP = expando()

        oP.x = point.x
        oP.y = point.y

        rotated.append((rotationMatrix[0][0] * (self.verticies[0].x - oP.x)) + (
            rotationMatrix[0][1] * (self.verticies[0].y - oP.y)) + oP.x)
        rotated.append((rotationMatrix[1][0] * (self.verticies[0].x - oP.x)) + (
            rotationMatrix[1][1] * (self.verticies[0].y - oP.y)) + oP.y)
        self.verticies[0].x = rotated[0]
        self.verticies[0].y = rotated[1]

        rotated = []
        rotated.append((rotationMatrix[0][0] * (self.verticies[1].x - oP.x)) + (
            rotationMatrix[0][1] * (self.verticies[1].y - oP.y)) + oP.x)
        rotated.append((rotationMatrix[1][0] * (self.verticies[1].x - oP.x)) + (
            rotationMatrix[1][1] * (self.verticies[1].y - oP.y)) + oP.y)
        self.verticies[1].x = rotated[0]
        self.verticies[1].y = rotated[1]


class ParallelLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NonIntersectingLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class CoincidentLines(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class Contour:
    def __init__(self):
        self.segments = []
        self.zlevel = 0.0
        self.closed = False

    @classmethod
    def withClipperPolygon(clazz, clipper_polygon, zlevel):
        self = clazz()
        self.zlevel = zlevel
        if clipper_polygon[0] != clipper_polygon[-1]:
            pass
            # logger.error('z: {} Given clipper polygon not closed.'.format(zlevel))
        for i in range(len(clipper_polygon)):
            # TODO: need to get p[next]
            p = clipper_polygon[i]
            nextidx = i + 1
            # should I really connect back to start again?
            if nextidx > len(clipper_polygon) - 1:
                nextidx = 0
                # continue

            p1 = clipper_polygon[nextidx]
            l = Line.withVerticies(Vertex(float(p.x) / CLIPPER_SCALE_FACTOR, float(p.y) / CLIPPER_SCALE_FACTOR, zlevel),
                                   Vertex(float(p1.x) / CLIPPER_SCALE_FACTOR, float(p1.y) / CLIPPER_SCALE_FACTOR,
                                          zlevel))
            # logger.info('Adding segment from clipper {0}'.format(l))
            self.segments.append(l)
        return self

    def __str__(self):
        return 'Contour: {0}'.format(id(self))

    '''
    Tests if this contour is closed. Meaning it's last segment's tail matches the 
    first segments head.
    '''

    def is_closed(self):
        if len(self.segments) == 0: return False
        return self.segments[0].verticies[0] == self.segments[-1].verticies[1]

    '''
    checks if this contour has any disconnected segments (which is generally a bad thing)
    '''

    def has_disconnected_segments(self):
        # first check that it's closed
        if not self.is_closed():
            print('## Segment not closed, therefore it\'s disconnected.')
        for x in range(len(self.segments)):
            if x < len(self.segments) - 1:
                cur = self.segments[x]
                next = self.segments[x + 1]
                # make sure that the tail of cur
                # matches the head of next
                if not cur.verticies[1] == next.verticies[0]:
                    print('## Contour has vertices that are not within tolerance.')


    '''
    Calculate the signed area of this polygon
    '''

    def winding_area(self):
        total = 0.0
        for segment in self.segments:
            total += segment.verticies[0].x * segment.verticies[1].y
            total += segment.verticies[1].x * segment.verticies[0].y
        return total / 2.0

    '''
    Return True if this contour's direction is clockwise.
    '''

    def clockwise(self):
        if self.winding_area() > 0.0:
            return True
        else:
            return False

    def intersect_brute_force(self, other):
        '''
        Returns the intersection points, if any, between this
        contour and another given contour.
        Uses brute force O(n2) approach.
        '''
        points = []
        for i in range(len(self.segments)):
            sega = self.segments[i]
            for j in range(len(other.segments)):
                segb = other.segments[j]
                try:
                    point = sega.intersect2D(segb)
                    points.append(point)
                except NonIntersectingLines:
                    pass
        return points

    def maybe_intersect(self, other):
        other_box = other.bounding_box()
        box = self.bounding_box()
        disjoint = False
        if (box.verticies[1].x < other_box.verticies[0].x):
            disjoint = True
        elif (other_box.verticies[1].x < box.verticies[0].x):
            disjoint = True
        elif (box.verticies[1].y < other_box.verticies[0].y):
            disjoint = True
        elif (box.verticies[1].y < other_box.verticies[0].y):
            disjoint = True
        return not disjoint

    def bounding_box(self):
        '''
        Return a bounding box representing the extents of this contour
        The box itself is a diagonal of the box
        '''

        currentMinX = sys.float_info.max
        currentMaxX = -sys.float_info.max
        currentMinY = currentMinX
        currentMaxY = currentMaxX
        z = self.segments[0].verticies[0].z
        for line in self.segments:
            for vertex in line.verticies:
                if vertex.x < currentMinX:
                    currentMinX = vertex.x
                if vertex.y < currentMinY:
                    currentMinY = vertex.y
                if vertex.x > currentMaxX:
                    currentMaxX = vertex.x
                if vertex.y > currentMaxY:
                    currentMaxY = vertex.y
        boundingDiagonal = Line.withVerticies(
            Vertex(currentMinX, currentMinY, z),
            Vertex(currentMaxX, currentMaxY, z)
        )
        return boundingDiagonal

    def inset(self, distance):
        '''
        Brute force creates an inset contour from this contour
        by the given distance. This does not handle sharp edges
        gracefully.
        '''
        inset = Contour()
        if inset is inset.clockwise():
            logger.warn(
                'Contour is not counter clockwise at z {0}'.format(
                    inset.zlevel)
            )
        for segment in self.segments:
            copy = segment.copy()
            copy.offset
            if self.zlevel == 0.0:
                logger.warn(
                    '\nOrigin Line: {0}\n Offset Line {1}\n Distance {2}'.format(
                        segment, copy, distance))
            # TODO: Find intersection points and trim
            last = None
            try:
                last = inset.segments[-1]
            except IndexError:
                pass
            if last is not None:
                try:
                    ua, ub = last.intersect2DUExtend(copy)
                    intersection_point = last.pointOnLine(ua)
                    if self.zlevel == 0.0:
                        logger.warn('Before trimming at {0}\nLast: {1} \nCopy {2}'.format(
                            last.pointOnLine(ua), last, copy))
                    last.verticies[0] = intersection_point
                    copy.verticies[1] = intersection_point
                    if self.zlevel == 0.0:
                        logger.warn(
                            'After trimming ua {0}\nLast: {1} \nCopy {2}'.format(ua, last, copy))
                except (ParallelLines):
                    logger.warn(
                        'Lines do not intersect. ParallelLines\n Last: {0} \n Copy {1}'.format(last, copy))
            inset.segments.append(copy)
        # need to trim last to first
        try:
            ua, ub = inset.segments[0].intersect2DUExtend(inset.segments[-1])
            intersection_point = inset.segments[0].pointOnLine(ua)
            inset.segments[0].verticies[1] = intersection_point
            inset.segments[-1].verticies[0] = intersection_point
        except CoincidentLines:
            logger.warn('\n{0} and \n{1} are coincident. no intersection.'.format(
                inset.segments[0], inset.segments[-1]))
        return inset

    def clipper_points(self):
        points = []
        # if not self.is_closed():
        # logger.error('Request for clipper_points from polygon that is not closed.')
        for line in self.segments:
            points.append(Point(line.verticies[0].clipperX(),
                                line.verticies[0].clipperY()))
        # Now the last point
        if len(self.segments) == 0:
            return points
        line = self.segments[-1]
        points.append(Point(line.verticies[1].clipperX(),
                            line.verticies[1].clipperY()))
        return points


def sortFacets():
    pass


def intersects(facet, z):
    '''Tests if the facet intersects the infinite plane at z'''

    intersects = False
    zcoords = []
    if (len(facet.verticies) != 3):
        raise Exception(
            'Found facet not having 3 verticies. OMG %d vertices found ' %
            len(facet.verticies))
    for vertex in facet.verticies:
        zcoords.append(vertex.z)
    zmin = min(zcoords)
    zmax = max(zcoords)
    if (z <= zmax and z >= zmin):
        # print 'found intersection %.15f <= %.15f <= %.15f ' % (zmin,z,zmax)
        intersects = True
    return intersects


def findIntersectingLines(facet, z):
    '''returns a list of lines in the given facet that intersect an\
     infinite plane at the given z coordinate.

    Keyword arguments:
    facet - The facet in which to search for intersecting lines
    z - The z coordinate of the infinite plane intersecting the facet
    '''

    # test za,zb
    intersections = []
    # skip this facet if it has one point above the slicing plane and two on
    if facet.verticies[0].z > z and Math.floatEq(facet.verticies[1].z, z) and Math.floatEq(facet.verticies[2].z, z):
        return intersections
    if facet.verticies[1].z > z and Math.floatEq(facet.verticies[0].z, z) and Math.floatEq(facet.verticies[2].z, z):
        return intersections
    if facet.verticies[2].z > z and Math.floatEq(facet.verticies[0].z, z) and Math.floatEq(facet.verticies[1].z, z):
        return intersections

    line = getIntersectingLine(facet.verticies[0], facet.verticies[1], z)
    if line is not None:
        intersections.append(line)
    line = None
    line = getIntersectingLine(facet.verticies[0], facet.verticies[2], z)
    if line is not None:
        intersections.append(line)
    line = None
    line = getIntersectingLine(facet.verticies[1], facet.verticies[2], z)
    if line is not None:
        intersections.append(line)
    return intersections


def getIntersectingLine(vertexa, vertexb, z):
    '''Returns a line from vertexa to vertexb if that line would \
    intersect an infinite plane at z'''

    line = None
    # Ignore facet edges that lie completely in the slicing plane

    if vertexa.z == z and vertexb.z == z:
        return line
    if (max(vertexa.z, vertexb.z) >= z and min(vertexb.z, vertexa.z) <= z):
        line = Line()
        line.verticies = []
        line.verticies.append(vertexa)
        line.verticies.append(vertexb)
    return line


def sliceAt(model, zslice):
    """find all facets that intersect this slice

    Keyword arguments:
    model - the model to be sliced
    zslice - the z coordinate at which to find intersecting facets

    Returns a list of facets that intersects the infinite plane at the
    given z coordinates.
    """
    logger = logging.getLogger()
    logger.info('================= slicing at %f =================' % (zslice))
    slicedFacets = []

    try:
        for facet in model.facets:
            if (intersects(facet, zslice) is True):
                slicedFacets.append(facet)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))

    return slicedFacets


def dotProduct(vectorA, vectorB):
    '''Calculate the dot product of two vectors'''
    return (vectorA.x * vectorB.x +
            vectorA.y * vectorB.y +
            vectorA.z * vectorB.z)


def subtractVectors(vectorA, vectorB):
    '''Subtract two vectors'''

    x = vectorA.x - vectorB.x
    y = vectorA.y - vectorB.y
    z = vectorA.z - vectorB.z
    result = Vertex(x, y, z)
    return result


def multiplyByScalar(vector, scalar):
    '''Multiply a vector by a scalar value'''

    x = vector.x * scalar
    y = vector.y * scalar
    z = vector.z * scalar
    result = Vertex(x, y, z)
    return result


def multiplyByScalar2D(vector, scalar):
    '''Multiply a vector by a scalar value only in 2D space'''

    x = vector.x * scalar
    y = vector.y * scalar
    z = vector.z
    result = Vertex(x, y, z)
    return result


def addVectors(vectorA, vectorB):
    '''Add two vectors'''

    x = vectorA.x + vectorB.x
    y = vectorA.y + vectorB.y
    z = vectorA.z + vectorB.z
    result = Vertex(x, y, z)
    return result


def getIntersectingPoints(lines, z):
    '''Returns a list of points at the intersection of the given \
    list of lines and in infinite plane at z. Paul Bourke: http://paulbourke.net/geometry/pointlineplane/'''

    # x-x1 / x2-x1 = y-y1 / y2-y1 = z - z1 / z2 - z1
    # P = A + t(B-A) eq of a line
    # N dot (P-C) = 0 equation of a plane
    # N - normal, P,C are points on the plane
    # Ax + By + Cz = D eq of a plane
    # points defining the slicing plane: 0,0,z 1,0,z 0,1,z
    # normal to the plane: 0,0,z to 0,0,z+1
    # A=0,B=0,C=(z+1)
    # 0x + 0y + (z+1)z=D
    # D=-(z+1)
    # z=Az+(Bz-Az)*t (z of equation of line)
    #
    points = []
    shouldLog = False
    lastS = -1
    for line in lines:
        u = subtractVectors(line.verticies[1], line.verticies[0])
        n = Vertex(0.0, 0.0, 1.0)
        nDotu = dotProduct(n, u)
        if shouldLog:
            logger.debug('processing line {0}'.format(line))
        if nDotu == 0:
            pass
            # no single point, all points on the line are in the cutting plane
            logger.debug('Line %s lies in cutting plane %f ' % (line, z))
            # print 'appending {0}'.format(line)
            # add both verticies
            points.extend(line.verticies)
        else:
            # v0: a point on the cutting plane

            v0 = Vertex(0.0, 0.0, z)
            s = (dotProduct(n, subtractVectors(v0, line.verticies[0])) /
                 dotProduct(n, subtractVectors(line.verticies[1],
                                               line.verticies[0])))
            if shouldLog:
                logger.debug('s %f ' % s)

            if (s < 0.0 or s > 1.0):
                raise Exception(
                    'Found a non intersecting line (s not in range [0-1]) \
                    in list that should only have intersecting lines.')

            # intersection point is P(s) = P0 + s(P1-P0)
            ps = addVectors(line.verticies[0], multiplyByScalar(
                subtractVectors(line.verticies[1], line.verticies[0]), s))
            '''
            if lastS == -1:
                lastS = s
            else:
                if lastS == 0.0 or lastS == 1.0:
                    if s != 0.0 and s != 1.0:
                        logger.info("Found case where intersection was midpoint and endpoint {} {}".format(points[0], ps))
                    if lastS != 0.0 and lastS != 1.0:
                        if s == 0.0 and s == 1.0:
                            logger.info("Found case where intersection was midpoint and endpoint {} {}".format(points[0], ps))
            '''

            if ps not in points:
                points.append(ps)
                if shouldLog:
                    logger.debug('appending p(%f)=%s' % (s, ps))
    return points


def layerRange(start, end, step):
    """range generator for layers"""
    while start < end:
        yield start
        start += step
    if start > end:
        start = end
    yield start


def circularRange(start, end, max, step=1):
    '''Range that increases until reaching end or max. If it reaches max
    it will reset at 0 and continue until end. Both start and end must be
    less than max.'''
    if (end < start): step = -step
    if start > max or end > max:
        raise Exception('The start and end parameters must both be less than max.')
    while True:
        yield start
        # we are done
        if start == end: break
        start += step
        # roll-over
        if start > max: start = 0


def frange(start, end, step):
    while start <= end:
        yield start
        start += step


def joinSegments(segments):
    logger = logging.getLogger()
    '''joins an unordered set of segments into a closed contour.'''

    if (segments is None or len(segments) == 0):
        return
    contour = Contour()

    contour.segments = []
    contour.segments.append(segments[0])
    lastSegment = segments[0]
    shouldLog = False
    segments.remove(lastSegment)
    while True:
        # connected = findConnectedSegment(lastSegment, segments)
        connected = findClosestSegment(lastSegment, segments)
        if (connected is None):
            if shouldLog:
                logger.debug(
                    'no remaining segments connect to {0}'.format(lastSegment))
            # try to close contour
            if lastSegment.verticies[1] == contour.segments[0].verticies[0]:
                # ?closingSegment = Line.withVerticies(lastSegment.verticies[0], contour.segments[0].verticies[1])
                # contour.segments.append(closingSegment)
                contour.closed = True
            break
        else:
            v1 = connected.verticies[1]
            duplicateStarts = filter(
                lambda line: line.verticies[1] == v1, contour.segments)
            if len(duplicateStarts) > 0:
                logger.warn('already exists a segment(s) with same start point. skipping'.
                            format(duplicateStarts[0], connected, len(connected.verticies)))
            else:
                contour.segments.append(connected)

            segments.remove(connected)
            lastSegment = connected
    if not contour.is_closed():
        l = Line.withVerticies(contour.segments[-1].verticies[1], contour.segments[0].verticies[0])
        logger.error('Created a contour that is not closed. Gap={}'.format(l.magnitude()))
    return contour


'''
Does a search of the segments in
'segments' looking for those where the tail matches
the head of the given segment. Will reverse a segments
if matching heads are found in order to create a head to tail
match.
'''


def findConnectedSegment(segment, segments):
    # print 'length of segments %d ' % len(segments)
    shouldLog = False
    for x in xrange(len(segments)):
        # print('comparing %s to %s [%d] ' %
        # (segment.verticies[1], segments[x].verticies[0], x))
        allowReverse = True
        if (allowReverse and segments[x].verticies[0] == segment.verticies[0]):
            segments[x].reverse()
            if shouldLog:
                logger.debug(
                    'found reversed match: segments[%s].vertices[1] == %s \
                     and segment.vertices[0] %s segments[x] %s'
                    % (x, segments[x].verticies[1], segment.verticies[0], segments[x]))
            # found a reversed line, turn it around and add it to the list
            return segments[x]

        if (segments[x].verticies[1] == segment.verticies[0]):
            if shouldLog:
                logger.debug('matching head to tail: {0} and {1}'.format(
                    segments[x].verticies[1], segment.verticies[0]))
            return segments[x]

    return None


'''
  Does a head to tail search returning the segment which is closest to the given segmentIndex
'''


def findClosestSegment(segment, segments):
    minDistance = sys.float_info.max
    closestSegment = -1
    for x in xrange(len(segments)):
        dsquared = Line.withVerticies(segment.verticies[1], segments[x].verticies[0]).magnitudeSquared()
        if dsquared < minDistance:
            minDistance = dsquared
            closestSegment = x
        dsquared = Line.withVerticies(segment.verticies[1], segments[x].verticies[1]).magnitudeSquared()
        if dsquared < minDistance:
            minDistance = dsquared
            segments[x].reverse()
            closestSegment = x

    if closestSegment > -1 and minDistance <= 0.1:
        return segments[closestSegment]
    else:
        return None


'''
Returns a list of unorder segments which represent
the exterior edges of the given set of facets
'''


def findExteriorEdges(facets):
    exteriorEdges = []
    for facet in facets:
        for edge in facet.edges:
            if len(edge.facets) == 1:
                exteriorEdges.append(edge)
    return exteriorEdges
    # ugh brute force
    for facet in facets:
        facetidx = facets.index(facet)
        for vertex in facet.verticies:
            # does any other facet share this segment?

            sharedEdge = False
            segment = Line()
            idx = facet.verticies.index(vertex)
            nextidx = idx + 1
            if nextidx > len(facet.verticies) - 1:
                nextidx = 0
            segment.verticies.append(vertex)
            segment.verticies.append(facet.verticies[nextidx])

            for innerfacet in facets:
                if (facets.index(innerfacet) == facetidx):
                    continue
                for innervertex in innerfacet.verticies:
                    innersegment = Line()
                    innersegment.verticies.append(innervertex)
                    innernext = innerfacet.verticies.index(innervertex) + 1
                    if (innernext > len(innerfacet.verticies) - 1):
                        innernext = 0
                    innersegment.verticies.append(
                        innerfacet.verticies[innernext])
                    # print ' checking: {0} to
                    # {1}'.format(segment,innersegment)
                    if segment == innersegment:
                        sharedEdge = True
                    if sharedEdge:
                        break
            if not sharedEdge:
                exteriorEdges.append(segment)
    return exteriorEdges


'''
Reduces the number of points in a contour to only those
required to keep a given tolerance.
'''


def simplifyContour(contour, tolerance=0.001, slicer=False):
    if slicer:
        clip = []
        # clip.append()
        fillType = PolyFillType.EvenOdd
        result = SimplifyPolygon(contour.clipper_points(), fillType)
        if len(result) == 0:
            return None
        return Contour.withClipperPolygon(result[0], contour.zlevel)
    newContour = Contour()
    newContour.segments = []
    newContour.zlevel = contour.zlevel
    # if not contour.is_closed():
    # logger.error('simplifyContour given an open contour.')
    while True:
        try:
            i = 0
            # print ' len(contour.segments) {0}'.format(len(contour.segments))
            while (i < len(contour.segments) - 1):
                # print 'segments {4} i {2} {0} i+1 {3} {1}'.format(contour.segments[i],
                # contour.segments[i+1],i,i+1,len(contour.segments))
                l = Line()
                l.verticies.append(contour.segments[i].verticies[0])
                l.verticies.append(contour.segments[i + 1].verticies[1])
                distance = l.distance(contour.segments[i].verticies[1])
                # print ' distance {0} to {1} is {2} l.magnitude() {3}'.format(l,
                # contour.segments[i].verticies[0],distance,l.magnitude())
                if (distance != -1 and distance <= tolerance):
                    newContour.segments.append(l)
                    i += 2
                else:
                    newContour.segments.append(contour.segments[i])
                    i += 1
                    # if we didn't simplify out the last segment
                    # then just add it in so our contour won't have a gap
                    # if i >= (len(contour.segments) -1):
                    #    logger.info(' appending last segment to newContour')
                    #    newContour.segments.append(contour.segments[-1])
            # print 'append (last) {0}'.format(contour.segments[-1])
            # test if there are multiple segments with the same start point
            # v1 = contour.segments[-1].verticies[0]
            # duplicateStarts = filter(
            # lambda line: line.verticies[0] == v1, newContour.segments)
            #if len(duplicateStarts) > 0:
            #    logger.debug('(simplifyContour) already exists a segment(s) {0} with same start point as {1} i= {2}'.
            #        format(duplicateStarts[0], contour.segments[-1], i))
            #else:
            # PROBLEM
            if i == len(contour.segments) - 1:
                newContour.segments.append(contour.segments[-1])
            #d = Line.withVerticies(contour.segments[-1].verticies[1], contour.segments[0].verticies[0]).magnitude()
            #if  d > 2.0:
            #   logger.warn("Contour not closed before simplification. magnitude {}".format(d))

            # print '{0} len: newContour {1} contour {2} '.format(
            #    current_process(),
            #    len(newContour.segments),len(contour.segments))
            # means we couldn't find any more segments within tolerance
            if (len(newContour.segments) == len(contour.segments)):
                break
            else:
                contour = newContour
                newContour = Contour()
                newContour.zlevel = contour.zlevel
                newContour.segments = []
        except:
            raise Exception(
                "".join(traceback.format_exception(*sys.exc_info())))
    if not newContour.is_closed():
        pass
        # logger.error('simplifyContour returning an open contour.')
    return newContour


'''
Simplifies the given list of contours using simplifyContour function.
'''


def simplifyContours(contours, tolerance=0.001):
    partialSimplifyContour = partial(simplifyContour, tolerance=tolerance)
    return map(partialSimplifyContour, contours)


def processLayerPrintExceptions(zcur, model, mintolerance=0.001):
    try:
        return process_layer(zcur, model, mintolerance)
    except:
        # Put all exception text into an exception and raise that
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


'''
Slices the given model at the given z coordinates. Returns the
contours found as a list.
'''


def offset_contours(contours, distance, z):
    clipper_union = []
    clip = []
    clipper = Clipper()
    for contour in contours:
        clip.append(contours.clipper_points())
    clipper.AddPolygons(clip, PolyType.Subject)
    pft = PolyFillType.EvenOdd
    result = clipper.Execute(ClipType.Union, clipper_union, pft, pft)
    offset = OffsetPolygons(clipper_union, int(-distance) * CLIPPER_SCALE_FACTOR)
    offset = CleanPolygons(offset)
    offset_contours = []
    for poly in offset:
        offset_contour = Contour.withClipperPolygon(poly, z)
        offset_contours.append(offset_contour)
    return offset_contours


def make_perimeters(pre_contours, model, zcur):
    clip = []
    clipper_union = []
    clipper = Clipper()
    for contour in pre_contours:
        clip.append(contour.clipper_points())
    clipper.AddPolygons(clip, PolyType.Subject)
    pft = PolyFillType.EvenOdd
    result = clipper.Execute(ClipType.Union, clipper_union, pft, pft)
    perimeters = []
    # plus one adds an offset for infill
    offsetAmount = (model.offset() / 2.0)
    for i in range(model.number_perimeters + 1):
        perimeter = []
        simpleOffset = OffsetPolygons(clipper_union, int(-offsetAmount * CLIPPER_SCALE_FACTOR))
        # simpleOffset = CleanPolygons(simpleOffset)
        for poly in simpleOffset:
            offsetContour = Contour.withClipperPolygon(poly, zcur)
            perimeter.append(offsetContour)
        perimeters.append(perimeter)
        offsetAmount += (model.offset() * model.perimeter_overlap_percent)
    return perimeters


def do_offset(contours, contours_from_stl, model, pre_contours, zcur):
    clip = []
    clipper = Clipper()
    if model.append_perimeters:
        contours.extend(contours_from_stl)
    for contour in pre_contours:
        # inset = contour.inset(0.5)
        # inset.zlevel = contour.zlevel
        # No. We want to offset first by nozzle diameter
        # contours.append(contour)
        # contours.append(inset)
        if contour is not None:
            clip.append(contour.clipper_points())
    clipper_union = []  # exterior
    # if len(pre_contours) > 1:
    clipper.AddPolygons(clip, PolyType.Subject)
    pft = PolyFillType.EvenOdd
    result = clipper.Execute(ClipType.Union, clipper_union, pft, pft)
    # else:
    # clipper_union = clip
    # create perimiters
    multiplier = 0
    # try something new here
    # assume that the width with increase in proportion to the decrease in layer height
    # offsetAmount = (model.nozzle_diameter / model.layerHeight) * model.nozzle_diameter
    overlapPercent = 0.5
    offsetAmount = model.offset() * overlapPercent
    for i in range(model.number_perimeters):
        if i > 0 and model.perimeter_overlap_percent:
            offsetAmount += offsetAmount + (model.offset() * model.perimeter_overlap_percent)
        multiplier = 1 + i
        simpleOffset = OffsetPolygons(clipper_union, int(-offsetAmount * CLIPPER_SCALE_FACTOR))
        simpleOffset = CleanPolygons(simpleOffset)
        if len(clipper_union) > 0 and len(simpleOffset) < 1:
            logger.warn("{} {} {} no polygons returned after offset operation".format(zcur, len(clip), int(
                -offsetAmount * CLIPPER_SCALE_FACTOR)))
        for poly in simpleOffset:
            offsetContourB = Contour.withClipperPolygon(poly, zcur)
            contours.append(offsetContourB)

    # offset for infill
    offsetAmount += (model.offset() * model.perimeter_overlap_percent)
    infillOffset = OffsetPolygons(clipper_union, int(-offsetAmount * CLIPPER_SCALE_FACTOR))
    toInfill = []
    for poly in infillOffset:
        infillOffsetContour = Contour.withClipperPolygon(poly, zcur)
        toInfill.append(infillOffsetContour)
    return toInfill


def process_layer(zcur, model, mintolerance=0.001):
    # global model
    facets = sliceAt(model, zcur)
    logger = logging.getLogger()
    logger.debug('processLayer: found %d facets at %f ' % (len(facets), zcur))
    segmentsAtZ = []
    coplanarFacetsAtZ = []

    for facet in facets:
        if (facet.isCoplanar(zcur)):
            # coplanarFacetsAtZ.append(facet)
            continue
        lines = findIntersectingLines(facet, zcur)
        # find the points on the intersecting lines that make up our slice
        # contour(s)
        points = getIntersectingPoints(lines, zcur)
        # print "z: %f facet: %s num intersecting lines %d num intersecting points %d" % (zcur,facet,len(lines),len(points))
        # if (len(lines) != len(points)):
        # raise Exception("intersection point count (%d) does not equal the number of intersecting lines (%d)" % (len(points),len(lines)))
        # print "%f # of intersecting points: %d " % (zcur,len(points))
        segment = Line()
        # segment.verticies = []
        shouldLog = False
        if shouldLog:
            logger.debug('{0} len(points) {1}'.format(zcur, len(points)))

        if len(points) == 2:
            segment.verticies.extend(points)
            # just add everything even if there are duplicates
            segmentsAtZ.append(segment)

            # TODO: Hash table is much faster than not in test
            '''if segment not in segmentsAtZ:
                if shouldLog:
                    logger.debug(
                        'append segment {0} len {1}'.format(segment, segment.magnitude()))
                segmentsAtZ.append(segment)
            else:
                if shouldLog:
                    logger.debug(
                        'segment already exists {0} len {1}'.format(segment, segment.magnitude()))
            '''

            # keep searching till we have sorted
            # all segments. there could be multiple
            # contours in a single slice (holes,etc)

    # Deal with facets at this layer. These
    # are exterior faces of the model in the z plane
    # alternative methods are instead to do boolean csg operations
    # on this slice and the one above and below to detect finished
    # infill areas and areas requiring support

    # clipper seems to obviate doing this
    # exteriorEdges = findExteriorEdges(coplanarFacetsAtZ)
    exteriorContours = []
    exteriorEdges = []
    # while len(exteriorEdges) > 0:
    # exteriorContour = joinSegments(exteriorEdges)
    # exteriorContour.zlevel = zcur
    #    exteriorContours.append(exteriorContour)

    # Handle perimiters

    # free some memory since we are not using these anymore
    #del model.facets[:]

    contours_from_stl = []
    # why? 
    segmentsAtZ.reverse()
    open_contours = []
    while len(segmentsAtZ) > 0:
        contour = joinSegments(segmentsAtZ)
        contour.zlevel = zcur
        if contour.closed:
            contours_from_stl.append(contour)
        else:
            open_contours.append(contour)
    # link up open contours
    if len(open_contours) > 0:
        logger.error("THERE ARE OPEN CONTOURS z:{}".format(zcur))
    for ocontour in open_contours:
        lastP = ocontour.segments[-1].verticies[1]
        startP = ocontour.segments[0].verticies[0]
        distanceSquared = Line.withVerticies(lastP, startP).magnitudeSquared()
        foundCloser = False
        # find another contour with a closer start point
        for inner_contour in open_contours:
            if inner_contour == ocontour:
                continue
            istartP = inner_contour.segments[0].verticies[0]
            d2 = Line.withVerticies(lastP, istartP).magnitudeSquared()
            if (d2 < distanceSquared):
                logger.error("found another contour that is closer. need to handle this case.")
                #foundCloser = True
        if not foundCloser:
            ocontour.segments.append(Line.withVerticies(lastP, startP))
            ocontour.closed = True
            # TEMP TO SEE IF THESE ARE THE ISSUE contours_from_stl.append(ocontour)

    pre_contours = simplifyContours(contours_from_stl, tolerance=model.nozzle_diameter / model.simplification_factor)
    # sort contours by area

    #offset = offset_contours(pre_contours, model.nozzle_diameter / 2.0, zcur)
    contours = []
    perimeters = make_perimeters(pre_contours, model, zcur)
    #for perimeter in perimeters:
    #    perimeter.sort(key=lambda contour: abs(contour.winding_area()))
    infill = []

    layer = Layer()
    layer.contours = contours
    layer.perimeters = perimeters
    layer.infill = infill
    layer.z = zcur
    #dump perimeters
    #layer.to_svg()
    return layer


'''
Returns a line from minx,miny to maxx,maxy of the given set of contours
'''


def bounds(contours, z):
    currentMinX = sys.float_info.max
    currentMaxX = -sys.float_info.max
    currentMinY = currentMinX
    currentMaxY = currentMaxX

    for contour in contours:
        for line in contour.segments:
            for vertex in line.verticies:
                if vertex.x < currentMinX:
                    currentMinX = vertex.x
                if vertex.y < currentMinY:
                    currentMinY = vertex.y
                if vertex.x > currentMaxX:
                    currentMaxX = vertex.x
                if vertex.y > currentMaxY:
                    currentMaxY = vertex.y

    boundingDiagonal = Line.withVerticies(Vertex(currentMinX, currentMinY, z),
                                          Vertex(currentMaxX, currentMaxY, z))
    return boundingDiagonal


'''
Returns a list of tuples containing the contours that potentially intersect,
but need more detailed testing
'''


def maybe_intersect(contours):
    results = []
    for i in range(len(contours)):
        for j in range(len(contours)):
            disjoint = False
            if i == j:
                pass
            else:
                disjoint = contours[i].maybe_intersect(contours[j])
            if not disjoint:
                results.append((contours[i], contours[j]))
    return results


def simpleLinearInfill(contours, z, spacing=2.0, minExtrude=1.0, angle=45.0, name='unnamed'):
    result = []
    # logger.info("z: {} # contours {}".format(z, len(contours)))
    boundingDiagonal = bounds(contours, z)
    p0 = boundingDiagonal.pointOnLine(-1.0)
    p1 = boundingDiagonal.pointOnLine(2.0)
    x1 = p0.x
    y1 = p0.y
    x2 = p1.x
    y2 = p1.y
    intersectionRows = []
    rowCount = 0
    curY = copy(y1)
    contourIntersections = {}
    polys = []

    while curY <= y2:
        v1 = Vertex(x1, curY, z)
        v2 = Vertex(x2, curY, z)
        intersections = []
        l = Line.withVerticies(v1, v2)

        # l.scale2D(4.0)
        # rotate about the line's midpoint
        l.rotate(angle, l.pointOnLine(0.5))
        poly = [Point(l.verticies[0].x, l.verticies[0].y), Point(l.verticies[1].x, l.verticies[1].y)]
        polys.append(poly)
        int_point_u_params = []
        for contour in contours:
            intersectionList = []
            try:
                intersectionList = contourIntersections[contour]
            except KeyError:
                # first time seeing this contour, add a list
                contourIntersections[contour] = intersectionList
            for segment in contour.segments:
                try:
                    u_param = l.intersect2DU(segment)
                    if (u_param < 0.0 or u_param > 1.0):
                        logger.info('u_param {} curY {}'.format(u_param, curY))
                    intersection = Intersection()
                    intersection.uparam = u_param
                    intersection.row = rowCount
                    intersection.ray = copy(l)
                    intersection.contour = contour
                    intersection.segmentIndex = contour.segments.index(segment)
                    intersection.z = z
                    intersections.append(intersection)
                    intersectionList.append(intersection.segmentIndex)
                # except ParallelLines:
                # logger.info('parallel lines z: {}'.format(z))
                except (ParallelLines, NonIntersectingLines, CoincidentLines):
                    continue

        rowCount += 1

        if len(intersections) % 2 == 1:
            logger.error("Intersection at z: {} y:{} resulted in odd number ({}) of intersections. ".format(z, curY,
                                                                                                            len(
                                                                                                                intersections)))
            # good chance we picked up a duplicate at some point
            Math.remove_duplicate_intersections(intersections)
            logger.error("After set conversion, len is {}".format(len(intersections)))

        sortReverse = False
        if (rowCount % 2 == 0): sortReverse = True
        intersections.sort(key=lambda x: x.uparam, reverse=sortReverse)
        # logger.info('sorted intersections {}'.format(intersections))

        if len(intersections) == 1:
            logger.warn('Only one intersection, skipping it')
            # empty list, pretent we never saw the intersection
            intersectionRows.append([])
        else:
            intersectionRows.append(intersections)
        curY += spacing

    # path following
    keepGoing = True
    i = 0;
    rowsProcessed = 0
    direction = 1
    while keepGoing:
        try:
            # logger.info ('z: {0} i: {1} len intersectionRows[{1}]: {2}'.format(z, i,len(intersectionRows[i])))
            startPlace = next(obj for obj in intersectionRows[i] if obj.used == False)
            # logger.info ('z: {} startPlace: {}'.format(z, startPlace))
            # logger.info ('z: {} intersectionRows[{}]: {}'.format(z, i, intersectionRows[i]))
            # process starting here
            startPlaceIdx = intersectionRows[i].index(startPlace)
            if startPlaceIdx % 2:
                direction = -1
            else:
                direction = 1
            stop = intersectionRows[i][startPlaceIdx + direction]
            # logger.info ('z: {} stop: {}'.format(z, stop))
            stop.used = True
            start = startPlace.point()
            startPlace.used = True
            end = stop.point()
            path = Line.withVerticies(start, end)
            if path.magnitude() > minExtrude: result.append(path)
            counter = 0

            while stop is not None:
                counter += 1
                if stop.row + 1 > len(intersectionRows) - 1:
                    nextStop = None
                else:
                    stop.intersectionIndex = intersectionRows[stop.row].index(stop)
                    nextStop = nextIntersection(stop, intersectionRows[stop.row + 1])

                #logger.info ('z: {} next intesection:{} from: {}'.format(z, nextStop, stop))   
                if nextStop is not None:
                    try:
                        nextStopPoint = nextStop.point()
                        path = Line.withVerticies(stop.point(), nextStopPoint)
                        #logger.info('z: {} nextIntersection: adding path {} to {}'.format(z, stop, nextStop))
                        #result.append(path)
                        # now add the "horizontal" path
                        curIdx = intersectionRows[stop.row].index(stop)
                        nextIdx = intersectionRows[stop.row + 1].index(nextStop)
                        #logger.info('nextIntersection:  stop idx on row is {} next: {}'.format(curIdx, nextIdx))

                        # if we are at an odd index, go backwards, else go forwards. never want to connect 1-2, 3-4 etc...
                        if nextIdx % 2:
                            direction = -1
                        else:
                            direction = 1

                        if nextIdx + direction - 1 > len(intersectionRows[stop.row + 1]):
                            direction = -1
                        if nextIdx == 0:
                            direction = 1

                        #logger.info('z: {} direction: {} len(intersectionRows): {} nextIdx: {}'.format(z, direction, len(intersectionRows[stop.row+1]), nextIdx))
                        nextStopH = intersectionRows[stop.row + 1][nextIdx + direction]
                        nextStopH.used = True
                        path = Line.withVerticies(nextStopPoint, nextStopH.point())
                        #logger.info('z: {} nextIntersection: adding h-path {} to {}'.format(z, nextStop, nextStopH))
                        nextStop = nextStopH
                        if path.magnitude() > minExtrude: result.append(path)
                    except IndexError, e:
                        logger.error('nextIntersection: index error {} {}'.format(e, intersectionRows[stop.row + 1]))
                stop = nextStop
                #if counter > 100: raise Exception('Hit max count')
        except StopIteration:
            # this row is empty
            rowsProcessed += 1
        except IndexError:
            pass
        # keepGoing = False
        # for row in intersectionRows:
        # if any(x.used == False for x in row):
        # keepGoing = True
        i += 1
        if i > len(intersectionRows) - 1: i = 0
        # Keep going if there are any unused intersections
        if rowsProcessed == len(intersectionRows):
            keepGoing = False
            i = 0
    if False:
        with open('infill-{}-{}.svg'.format(name, z), "w") as f:
            f.write('<?xml version="1.0" standalone="no"?>\n')
            f.write(
                '<svg width="200mm" height="200mm" viewBox="-100 -100 200 200" version="1.1" xmlns="http://www.w3.org/2000/svg">\n')
            f.write('<path d="')
            for i in range(len(polys)):
                f.write(' M{} {}'.format(polys[i][0].x, polys[i][0].y))
                f.write(' L{} {}'.format(polys[i][1].x, polys[i][1].y))
            f.write('" stroke="red" fill="transparent" stroke-width=".01"/>\n')

            f.write('<path d="')
            for path in result:
                f.write(' M{} {}'.format(path.verticies[0].x, path.verticies[0].y))
                f.write(' L{} {}'.format(path.verticies[1].x, path.verticies[1].y))
            f.write('" stroke="blue" fill="transparent" stroke-width=".1"/>\n')

            f.write('</svg>')
    return result


def lineIntersectsContour(line, contour):
    result = False
    for segment in contour.segments:
        try:
            point = line.intersect2DU(segment)
            # if point > 0.0 and point < 1.0:
            result = True
            break
        except (ParallelLines, NonIntersectingLines, CoincidentLines):
            continue

    return result


'''
 Gets the next intersection point on a contour
'''


def nextIntersection(intersection, intersectionRow):
    intersectedSegment = intersection.contour.segments[intersection.segmentIndex]
    # TODO: find intersection with same contour and closest intersection index
    candidates = []
    for intInRow in intersectionRow:
        if intInRow.contour is intersection.contour and intInRow.used == False:
            # possible candidate next intersections
            # make sure we only connect odd index to odd and even index to even
            if intersection.intersectionIndex % 2 == intersectionRow.index(intInRow) % 2:
                continue
            else:
                candidates.append(intInRow)
    # logger.info('nextIntersection() candidates {}'.format(candidates))
    # find candidate with closest intersection index to the given intersection
    if len(candidates) == 0:
        return None
    # need to test "circular distance" since two indecies may be closer depending on the direction traveled along the contour
    # logger.info('z: {} len segments {}'.format(intersection.z, len(intersection.contour.segments)))
    # cases when this doesn't work, like if there are lots of little segments    
    answer = min(candidates, key=lambda x: min(abs(x.segmentIndex - intersection.segmentIndex),
                                               len(intersection.contour.segments) - abs(
                                                   x.segmentIndex - intersection.segmentIndex)))
    # this is just linear distance...
    # answer = min(candidates, key=lambda x: Line.withVerticies(intersection.point(),x.point()).magnitude() )
    # answer = min(candidates, key=lambda x:min(abs(x)))
    answer.used = True
    #logger.info('nextIntersection() Closest intersection to {} is {}'.format(intersection, answer))    
    return answer


'''
Generate infill for the given set of contours as part of the given model
'''


def infillContours(contours, z, spacing=0.5):
    llogger = logger
    boundingDiagonal = bounds(contours, z)
    x1 = boundingDiagonal.verticies[0].x
    y1 = boundingDiagonal.verticies[0].y
    x2 = boundingDiagonal.verticies[1].x
    y2 = boundingDiagonal.verticies[0].y
    v1 = Vertex(x1, y1, z)
    v2 = Vertex(x2, y2, z)
    llogger.info('Num Contours at z {0} is {1}'.format(z, len(contours)))
    infillPaths = []
    pointToSegment = {}
    intersectionsAtY = {}
    contourIntersections = {}
    # keep track of the max number of intersections at any Y
    # so we know how many passes to take
    maxInfilPathsAtY = 0
    reverseRay = False
    while v1.y <= boundingDiagonal.verticies[1].y:
        llogger.debug('v1 {0} boundingDiagonal.verticies[1].y {1}'.format(
            v1, boundingDiagonal.verticies[1].y))
        pathsAtY = []
        intersectionsAtY[v1.y] = pathsAtY
        l = Line.withVerticies(v1, v2)
        if reverseRay: l.reverse()
        reverseRay = not reverseRay
        points = []
        uParamToSegment = {}
        contourHitCounts = {}
        for contour in contours:
            for segment in contour.segments:
                try:
                    point = l.intersect2DU(segment)
                    # try:
                    # curCount = contourHitCounts[contour]
                    # curCount = curCount + 1
                    # llogger.info('hit count {}'.format(curCount))
                    # except KeyError:
                    #    contourHitCounts[contour] = 1
                    #    llogger.info('First hit on contour {}'.format(contour))

                    if point not in points:
                        points.append(point)
                        uParamToSegment[point] = (segment, contour)
                    llogger.debug('contour {5} intersection u param {0} y {1} z {2} with line {3} of length {4}'.format(
                        point, v1.y, z, segment, segment.magnitude(), contour))
                except (ParallelLines, NonIntersectingLines, CoincidentLines):
                    continue

        llogger.debug(
            'found {0} intersections at y {1} z {2}'.format(len(points), v1.y, z))
        points.sort()

        lastP = None

        for i in range(len(points)):
            if (i % 2 != 0):
                lastP = l.pointOnLine(points[i - 1])
                p = l.pointOnLine(points[i])
                llogger.debug('Y:{0} Scanline from u {1} to u {2}'.format(v1.y, points[i - 1], points[i]))
                # track the segment which caused the intersection
                # so we can trace back over it later
                (testSegment, testContour) = uParamToSegment[points[i]]
                infillPath = Line.withVerticies(lastP, p)
                pointToSegment[p] = (testSegment, testContour, infillPath)
                (prevtestSegment, prevtestContour) = uParamToSegment[points[i - 1]]
                pointToSegment[lastP] = (prevtestSegment, prevtestContour, infillPath)
                llogger.debug('adding path {0} i: {1}'.format(infillPath, i))
                intersectionList = []
                try:
                    intersectionList = contourIntersections[testContour]
                except KeyError:
                    # first time seeing this contour, add a list 
                    contourIntersections[contour] = intersectionList
                # Tuple of index of the segment and the intersection point
                # (
                llogger.debug('adding intersection with contour {0} at index {1} for point {2}'.format(
                    testContour, testContour.segments.index(testSegment), p))
                intersectionList.append((testContour.segments.index(testSegment), p))
                llogger.debug('adding (previous) intersection with contour {0} at index {1} for point {2}'.format(
                    prevtestContour, prevtestContour.segments.index(prevtestSegment), lastP))
                intersectionList.append((prevtestContour.segments.index(prevtestSegment), lastP))
                infillPaths.append(infillPath)
                pathsAtY.append(infillPath)
                pathsLen = len(pathsAtY)
                if (pathsLen > maxInfilPathsAtY): maxInfilPathsAtY = pathsLen

        v1.y += spacing
        v2.y += spacing


    # Now make another pass to connect all the infill segments
    # The should be connected by using the points of the contour
    # the intersect between each 'scan line' rather than just straght lines
    llogger.debug("Count of pointToSegment {0}".format(len(pointToSegment)))
    hatchedInfill = []
    reversePath = True
    reverseConnectionTravel = True
    nextInfillPath = None
    previousPathY = None
    usedConnections = []
    usedPoints = []
    while True:
        # this may have to alternate between 1 and 0
        # don't need this anymore because of call to reverse()
        if reversePath:
            vertexIndex = 0
        else:
            vertexIndex = 1

        try:
            if not nextInfillPath:
                path = infillPaths.pop()
                reverseConnectionTravel = not reverseConnectionTravel
                llogger.info('removed through pop() {0}'.format(path))
            else:
                path = nextInfillPath
                llogger.info('attempting to remove {0}'.format(nextInfillPath))
                try:
                    infillPaths.remove(nextInfillPath)
                except ValueError:
                    llogger.error('Already removed {0}'.format(nextInfillPath))
                    path = infillPaths.pop()
                    reverseConnectionTravel = not reverseConnectionTravel
                    llogger.info('Caught ValueError: removed through pop() {0}'.format(path))
            # don't allow connections back to this start point
            usedPoints.append(path.verticies[0])
            llogger.info('hatching path {0} vertex index {1} reverse: {2}'.format(path, vertexIndex, reversePath))
        except IndexError:
            # list is now empty
            llogger.info('No more infill paths left to pop')
            break

        # actually this doesn't work when there are multiple
        # objects. unless we can process them separately(?)
        # instead we always reverse
        # if previousPathY != path.verticies[1].y:
        # TRYING TO REVERSE RAYS WHEN THEY ARE CREATED INSTEAD
        # reversePath = not reversePath
        # if (reversePath): path.reverse()

        hatchedInfill.append(path)

        lastPoint = path.verticies[1]
        if lastPoint in usedPoints:
            llogger.info('ALREADY USED POINT {}'.format(lastPoint))
        else:
            usedPoints.append(lastPoint)
        llogger.info('intersection point (lastPoint) {0}'.format(lastPoint))
        connectingSegment = None
        try:
            connectingSegment, connectingContour, infillPath = pointToSegment[lastPoint]
            llogger.info('connectingContour {0} connectingSegment {1}'.format(connectingContour, connectingSegment))
            intersectionList = contourIntersections[connectingContour]
            # sort the intersection list, but not so many times!
            intersectionList.sort(key=itemgetter(0))
            # llogger.info('INTERSECTION POINTS')
            # for inter in intersectionList:
            #    llogger.info('sort on contour index: {} intersection point {}'.format(inter[0],inter[1]))
            #for interSectionTuple in intersectionList:
            #    llogger.info('intersection idx {0} point {1}'.format(interSectionTuple[0],interSectionTuple[1]))

            # this triggers an error when the list contains a tuple
            #intersectionPoint = intersectionList.index(lastPoint)

            # find the index in the list where our current point is
            # the next item in the list is the next intersection index
            # we should take all points between the two
            startIndex = -1
            connectingSegmentIdx = -1
            for idx, tup in enumerate(intersectionList):
                if tup[1] == lastPoint:
                    startIndex = idx
                    connectingSegmentIdx = tup[0]
                    break

            llogger.info('Start index {0}'.format(startIndex))
            #startIndex = connectingContour.segments.index(connectingSegment)
            # now append the section of the connectingContour from connectingSegment1
            # to 'next' connecting segment
            #if reversePath:
            #    nextidx = idx+1
            #else:
            #    nextidx = idx -1
            # TODO: Somehow need to know if we are hitting a hole wall
            # and do idx-2 instead
            direction = -1
            if reverseConnectionTravel: direction = 1
            nextidx = startIndex + direction
            if (nextidx > len(intersectionList) - 1): nextidx = 0
            llogger.info('next index is {} reverseConnectionTravel {}'.format(nextidx, reverseConnectionTravel))
            if nextidx < 0:
                llogger.info('nextidx is negative. done with current pass.')
                # Is this correct?
                nextInfillPath = None
            else:
                try:
                    next_tup = intersectionList[nextidx]
                    #prev_tup = intersectionList[startIndex+1]
                    llogger.info('Found next intersections at idx {0}'.format(next_tup[0]))
                    #llogger.info('previous  intersections at idx {0}. connectingSegmentIdx {1} '.format(prev_tup[0], connectingSegmentIdx))
                except IndexError:
                    llogger.info('No more intersections at idx {0}'.format(nextidx))
                    break

                llogger.info('{0} intersectionlist index {2} CONNECTING SEGMENT IDX {1}  NEXT SEG IDX {3}'.format(i,
                                                                                                                  connectingSegmentIdx,
                                                                                                                  startIndex,
                                                                                                                  next_tup[
                                                                                                                      0]))
                # check if we are re-useding a connection. that's bad.
                if connectingSegmentIdx in usedConnections and connectingSegmentIdx != next_tup[0]:
                    llogger.info('START INDEX {} already used'.format(startIndex))
                    #nextInfillPath = None
                    #break
                elif next_tup[0] in usedConnections and next_tup[0] != connectingSegmentIdx:
                    llogger.info('NEXT INDEX {} already used'.format(next_tup[0]))
                    #nextInfillPath = None
                    #break
                else:
                    usedConnections.append(startIndex)
                    usedConnections.append(next_tup[0])
                # Slic3r and Cura basically does this and doesn't try to follow
                # the curves (add 'or True' to the statement below to behave like Slic3r, Cura)
                if next_tup[0] == connectingSegmentIdx:
                    # don't reverse this case
                    if reverseConnectionTravel:
                        nextidx = startIndex - 1
                        next_tup = intersectionList[nextidx]
                    connecting_path = infillPath = Line.withVerticies(lastPoint, next_tup[1])
                    llogger.info('Next intersection on same segment. Straight shot. {0}'.format(connecting_path))
                    if next_tup[1] in usedPoints:
                        llogger.info('ALREADY USED POINT {}'.format(next_tup[1]))
                    else:
                        usedPoints.append(next_tup[1])
                        # Don't add the path if it contains an already
                        # used point... But actually we probably have a bug
                        hatchedInfill.append(connecting_path)
                else:
                    # this is a multi-segment path
                    llogger.info('multi-segment path from {0} to {1} clockwise={2}'.format(connectingSegmentIdx,
                                                                                           next_tup[0],
                                                                                           connectingContour.clockwise()))

                    # deal with continuing look back around to a starting point in consistent direction
                    ## Changed step to -1
                    # need to get the number of segments in the connecting contour as rangeMax
                    rangeMax = len(connectingContour.segments)
                    # when going the other direction we need prev_tup, not next?
                    for i in circularRange(connectingSegmentIdx, next_tup[0], rangeMax):
                        llogger.info('circular range start {0} end {1} max {2} result {3} '.format(connectingSegmentIdx,
                                                                                                   next_tup[0],
                                                                                                   rangeMax, i))
                        if i == connectingSegmentIdx:
                            # wrong hatchedInfill.append(Line.withVerticies(next_tup[1],connectingContour.segments[i].verticies[1]))
                            seg = Line.withVerticies(lastPoint, connectingContour.segments[i].verticies[0])
                            llogger.info('{0} (start) adding segment: {1}'.format(i, seg))
                            hatchedInfill.append(seg)
                        elif i == next_tup[0]:
                            seg = Line.withVerticies(next_tup[1], connectingContour.segments[i].verticies[1])
                            llogger.info('{0} (end) adding segment: {1}'.format(i, seg))
                            hatchedInfill.append(seg)
                        else:
                            try:
                                logger.info('{0} adding segment: {1}'.format(i, connectingContour.segments[i]))
                                hatchedInfill.append(connectingContour.segments[i])
                            except IndexError:
                                logger.info('No segment at index {0} in connecting contour.'.format(i))
                                pass


                # now get the next ray (infillPath) that caused this intersection
                connectingSegment, connectingContour, nextInfillPath = pointToSegment[next_tup[1]]
                if nextInfillPath == path:
                    llogger.warn('Next infill path is same as current path.. Skipping.')
                    nextInfillPath = None
                llogger.info('Next infill path {} connectingSegment {} connectingContour {}'.format(nextInfillPath,
                                                                                                    connectingSegment,
                                                                                                    connectingContour))


            # Note: This won't work with rotated rays we'll need to know if the lines
            # are the same original ray
            previousPathY = path.verticies[1].y
        except KeyError:
            llogger.info('No pointToSegment for point {0}'.format(lastPoint))
            pass

    return hatchedInfill


def slice(model):
    logger.info('slicing %s %f %f %f' %
                (model.name, model.zmin, model.zmax, model.layerHeight))
    pool = Pool(model.processing_threads)
    partialProcessLayer = partial(processLayerPrintExceptions, model=model)
    results = pool.map(partialProcessLayer, layerRange(
        model.zmin, model.zmax, model.layerHeight))
    pool.close()
    logger.debug('waiting for threads to complete...')
    pool.join()
    logger.debug('threads done')
    return results
    # for zcur in layerRange(model.zmin, model.zmax, model.layerHeight):
    # processLayer(zcur)


def processLine(line, model):
    stripped = line.strip()
    if (stripped == 'outer loop'):
        model.triangles += 1
    elif (stripped == 'endloop'):
        model.endLoop += 1
    elif (stripped.startswith('facet')):
        facet = Facet()
        model.facets.append(facet)
    elif (stripped.startswith('endfacet')):
        if model.facets[-1].isCoplanar(model.facets[-1].verticies[0].z):
            for i in range(3):
                next = i + 1
                if next > 2:
                    next = 0
                edge = Line.withVerticies(
                    model.facets[-1].verticies[i], model.facets[-1].verticies[next])
                # try:
                # idx = model.edgeList.index(edge)
                # print 'edge {0} found at {1} {2}'.format(edge,idx,model.edgeList[idx])
                # except ValueError:
                # model.edgeList.append(edge)

                try:
                    edge = model.edges[edge.key()]
                    # print 'found edge in edges dictionary {0}'.format(edge)
                except KeyError:
                    # print 'could not find edge {0}'.format(edge.key())
                    model.edges[edge.key()] = edge

                if (model.facets[-1].inZPlane()):
                    edge.facets.append(model.facets[-1])
                # print 'edge facet count {0}'.format(len(edge.facets))
                model.facets[-1].edges.append(edge)

    elif (stripped.startswith('vertex')):
        coords = stripped.split()
        v = Vertex(float(coords[1]), float(coords[2]), float(coords[3]))
        if (model.facets[-1].verticies is None):
            model.facets[-1].verticies = []
        model.facets[-1].verticies.append(v)
        zcur = float(coords[3])
        updateModelZRange(zcur, model)


def updateModelZRange(zcur, model):
    if (zcur < model.zmin):
        model.zmin = zcur
    elif (zcur > model.zmax):
        model.zmax = zcur


def readFile(name, model):
    with open(name, "rb") as infile:
        model.name = infile.name
        firstLine = infile.readline()
        if not firstLine.startswith('solid'):
            infile.seek(0)
            processBinarySTL(infile, model)
        else:
            infile.seek(0)
            processASCIISTL(infile, model)

    return model


def processASCIISTL(file, model):
    for line in file:
        processLine(line, model)


def processBinarySTL(fp, model):
    header = fp.read(80)
    num_triangles = unpack('<i', fp.read(4))[0]
    logger.info(header)
    logger.info('Num Triangles: {}'.format(num_triangles))
    for i in range(num_triangles):
        read_triangle(fp, model)


def read_triangle(f, model):
    # skip the normals
    f.seek(12, 1)
    # read three triangle verticies
    p1 = read_point(f)
    p2 = read_point(f)
    p3 = read_point(f)
    # skip Attribute byte count
    f.seek(2, 1)
    facet = Facet()
    facet.verticies.append(p1)
    facet.verticies.append(p2)
    facet.verticies.append(p3)
    model.facets.append(facet)
    zcur = p1.z
    updateModelZRange(zcur, model)


def read_point(f):
    vertex = unpack("<fff", f.read(12))
    v = Vertex(vertex[0], vertex[1], vertex[2])
    return v


'''
Test functions
'''


def test():
    # test contour closed
    c = Contour()
    if c.is_closed():
        raise Exception('An empty contour should not be considered closed.')
        # test line intersection
    seg1 = Line.withVerticies(Vertex(0.0, 0.0, 0.0),
                              Vertex(2.0, 2.0, 0.0))
    seg2 = Line.withVerticies(Vertex(2.0, 0.0, 0.0),
                              Vertex(0.0, 2.0, 0.0))
    c.segments.append(seg1)
    c.segments.append(seg2)
    if c.is_closed():
        raise Exception('This contour is not closed, but is_closed returned true.')
    seg1 = Line.withVerticies(Vertex(0.0, 2.0, 0.0),
                              Vertex(0.0, 0.0, 0.0))
    c.segments.append(seg1)
    if not c.is_closed():
        raise Exception('This contour is closed, but is_closed returned false.')

    l1 = Line.withVerticies(Vertex(1.0, 0.0, 0.0),
                            Vertex(32.543, 10.4999, 0.0))
    # test normal
    l1normal = l1.normal()
    v1 = subtractVectors(l1.verticies[1], l1.verticies[0])
    v2 = subtractVectors(l1normal.verticies[1], l1normal.verticies[0])
    if dotProduct(v1, v2) != 0.0:
        raise Exception('Dot product of vector and its normal must be zero')

    # test unit vector
    unitVector = l1.unitVector()
    if not Math.floatEq(unitVector.magnitude(), 1.0):
        raise Exception(
            'Magnitude of unit vector must be 1.0 was {0:.45f}'.format(unitVector.magnitude()))

    # get a point at u=.5 should be 0.5,0.0,0.0
    result = l1.pointOnLine(1)
    print 'point at u=.5 {0}'.format(result)

    # test line intersection
    lint = Line.withVerticies(Vertex(0.0, 0.0, 0.0),
                              Vertex(2.0, 2.0, 0.0))
    lint2 = Line.withVerticies(Vertex(2.0, 0.0, 0.0),
                               Vertex(0.0, 2.0, 0.0))
    vint = lint.intersect2D(lint2)
    if not vint == Vertex(1.0, 1.0, 0.0):
        raise Exception('intersection should be {0} was {1}'.format(
            Vertex(1.0, 1.0, 0.0), vint))
    print 'intersection between {0} and {1} is {2}'.format(
        lint, lint2, vint)

    line = Line()
    v1 = Vertex(1.0, 1.0, 1.0)
    v2 = Vertex(2.0, 1.0, 1.0)
    line.verticies.append(v1)
    line.verticies.append(v2)

    l2 = Line()
    l2.verticies.append(v1)
    l2.verticies.append(v2)
    if (line != l2):
        raise Exception('line != l2')
    length = line.magnitude()

    point = Vertex(1.5, 3.6173063927, 1.0)

    distance = line.distance(point)
    if (distance < 1.0):
        raise Exception('distance should be 2 but was {0}'.format(distance))
    print 'distance {0} length {1}'.format(distance, length)

    # find edges test
    va = Vertex(0.0, 0.0, 0.0)
    vb = Vertex(1.0, 0.0, 0.0)
    vc = Vertex(1.0, 1.0, 0.0)
    vd = Vertex(0.0, 1.0, 0.0)
    ve = Vertex(1.0, 2.0, 0.0)
    vf = Vertex(0.0, -1.0, 0.0)
    verticies = [va, vb, vc, vd, ve, vf]
    # test if equal can find an item by index
    vdd = Vertex(0.0, 1.0, 0.0000000000000)
    vddidx = verticies.index(vdd)
    if vddidx != 3:
        raise Exception('vddidx should be 3')
    f1 = Facet.withVerticies(va, vb, vc)
    f2 = Facet.withVerticies(va, vc, vd)
    f3 = Facet.withVerticies(vc, ve, vd)
    f4 = Facet.withVerticies(va, vb, vf)

    # Rotation Test
    rotLine = Line.withVerticies(Vertex(0.0, 0.0, 0.0), Vertex(1.0, 0.0, 0.0))
    rotLine.rotate(90.0, Vertex(0.0, 0.0, 0.0))
    print 'rotLine {}'.format(rotLine)
    if not Math.floatEq(rotLine.verticies[1].x, 0.0) or not Math.floatEq(rotLine.verticies[1].y, 1.0):
        raise Exception('Rotate test failed {}'.format(rotLine))

    # Rotation Test
    rotLine = Line.withVerticies(Vertex(0.0, 0.0, 0.0), Vertex(1.0, 0.0, 0.0))
    rotLine.rotate(90.0, rotLine.pointOnLine(0.5))
    print 'rotLine2 {}'.format(rotLine)
    if not Math.floatEq(rotLine.verticies[0].x, 0.5) or not Math.floatEq(rotLine.verticies[1].y, 0.5):
        raise Exception('Rotate test failed {}'.format(rotLine))


        # this test no long applies since the
        # function requires setup done by reading the stl file
        # testEdgesList = [f1, f2, f3, f4]
        # exteriorEdges = findExteriorEdges(testEdgesList)
        # if (len(exteriorEdges) != len(testEdgesList) + 2):
        # raise Exception(
        # 'Exterior edges count should be {1}. Instead was {0}'.format(len(exteriorEdges),
        # len(testEdgesList) + 2))
        # raise Exception('tests complete.')


if __name__ == "__main__":
    # logging.basicConfig(level=logging.INFO)
    # logger = logging.getLogger()
    parser = argparse.ArgumentParser(description='An STL slicer for producing 3D printer gcode.')
    parser.add_argument('stl', help="STL file to be sliced.")
    parser.add_argument('output', help="Output path for generated gcode.")
    parser.add_argument('-p', '--perimeters_only', help="Generate only perimeters, no infill.", action="store_true")
    parser.add_argument('-a', '--append_perimeters',
                        help="Include original stl perimeters with no offsetting in the output.", action="store_true")
    parser.add_argument('-v', '--verbose', help="Be verbose.", action="store_true")
    parser.add_argument('-o', '--perimeter_overlap_percent',
                        help="Percent overlap between perimeters. Smaller results in more overlap.", type=float)
    parser.add_argument('-n', '--num-perimeters', help="Minimum number of perimeters.", type=int)
    parser.add_argument('-d', '--filament-diameter',
                        help="Measured diameter of filament to at least two decimal places.", type=float, default=1.75)
    parser.add_argument('-l', '--layer-height', help="Slicing layer height", type=float, default=0.1)
    parser.set_defaults(perimeter_overlap_percent=1.0)
    parser.set_defaults(num_perimeters=3)
    args = parser.parse_args()
    logger = logging.getLogger()

    ch = logging.StreamHandler(sys.stdout)
    logger.addHandler(ch)
    logger.setLevel(logging.INFO)

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    sys.setrecursionlimit(6000)
    # global model
    # test()
    dt = datetime.now()
    t0 = dt.microsecond
    # print('Memory usage: %s (MB)'
    # % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024))
    logger.info('Reading...')
    with Timer() as totalTime:
        with Timer() as readTime:
            model = Model()
            model.perimeters_only = args.perimeters_only
            model.append_perimeters = args.append_perimeters
            model.perimeter_overlap_percent = args.perimeter_overlap_percent
            model.number_perimeters = args.num_perimeters
            model.filament_diameter = args.filament_diameter
            model.layerHeight = args.layer_height
            readFile(args.stl, model)
        logger.info('File read took {0} seconds'.format(readTime.secs))
        # check for open facets in model
        #num_open = model.countOpenFacets()
        #logger.info('Facets:  open: {}  total: {}'.format(num_open, len(model.facets)))

        logger.info('Slicing...')
        # list of tuples containing contours,infill
        with Timer() as sliceTime:
            layers = slice(model)
        logger.warn('Slicing took {0} seconds. {1} layers'.format(sliceTime.secs, len(layers)))


        # sort results list by zlevel
        layers.sort(key=attrgetter('z'))

        roofs = []

        for layerNum in xrange(1, len(layers)):
            logger.info('Z: {}'.format(layers[layerNum].z))
            zcur = layers[layerNum].z
            previousLayerPolys = []
            currentLayerPolys = []
            nextLayerPolys = []

            for contour in layers[layerNum].perimeters[-1]:
                currentLayerPolys.append(contour.clipper_points())

            for contour in layers[layerNum - 1].perimeters[-1]:
                if contour.is_closed():
                    previousLayerPolys.append(contour.clipper_points())

            if layerNum < len(layers) - 2:
                for contour in layers[layerNum + 1].perimeters[-1]:
                    nextLayerPolys.append(contour.clipper_points())

            cp = Clipper()
            cp.AddPolygons(currentLayerPolys, PolyType.Clip)
            cp.AddPolygons(previousLayerPolys, PolyType.Subject)
            sln = []
            pft = PolyFillType.EvenOdd
            clipperResult = cp.Execute(ClipType.Difference, sln, pft, pft)
            clean_sln = CleanPolygons(sln)
            if clipperResult:
                for poly in sln:
                    solutionContour = Contour.withClipperPolygon(poly, zcur)
                    layers[layerNum - 1].roofs.append(solutionContour)
                    # can't do this since it screws with results resultsList[layerNum][0].append(solutionContour)
            #logger.info(
            #    'clipperResult {} solution len {} len(clip) {} len (subject) {}'.format(clipperResult, len(clean_sln),
            #                                                                            len(currentLayerPolys),
            #                                                                            len(previousLayerPolys)))

            # Subtract the roof from the current layer contours to give us the area for "normal" infill
            cp = Clipper()
            cp.AddPolygons(previousLayerPolys, PolyType.Subject)
            cp.AddPolygons(clean_sln, PolyType.Clip)
            previousMinusRoof = []
            result = cp.Execute(ClipType.Difference, previousMinusRoof, pft, pft)
            if result:
                for poly in previousMinusRoof:
                    contour = Contour.withClipperPolygon(poly, zcur)
                    layers[layerNum - 1].normalInfill.append(contour)

            # now do bridges / overhang
            cp = Clipper()
            cp.AddPolygons(currentLayerPolys, PolyType.Clip)
            cp.AddPolygons(nextLayerPolys, PolyType.Subject)
            sln = []
            result = cp.Execute(ClipType.Difference, sln, pft, pft)
            if result:
                for poly in sln:
                    contour = Contour.withClipperPolygon(poly, zcur)
                    layers[layerNum + 1].overhang.append(contour)

        logger.info('Starting infill...')
        with Timer() as infillTime:
            for layerNum in xrange(len(layers)):
                z = layers[layerNum].z
                if layerNum % 2:
                    angle = model.infill_angle
                else:
                    angle = -model.infill_angle
                # TODO: Calculate roof thickness rather than hard coding a number of layers
                if layerNum < 4 or layerNum > len(layers) - 4:
                    logger.info('z: {} perimeter count {} contour count {}'.format(z, len(layers[layerNum].perimeters),
                                                                                   len(layers[layerNum].perimeters[
                                                                                       -1])))
                    layers[layerNum].infill = simpleLinearInfill(layers[layerNum].perimeters[-1], z,
                                                                 spacing=model.offset(), minExtrude=1.0,
                                                                 angle=angle, name='normal')
                else:
                    infill = simpleLinearInfill(layers[layerNum].normalInfill, z, spacing=2.0, minExtrude=1.0,
                                                angle=angle, name='normal')
                    layers[layerNum].infill = infill
                    roofInfill = simpleLinearInfill(layers[layerNum].roofs, z, spacing=model.offset(),
                                                    minExtrude=0.5, angle=angle, name='roof')
                    layers[layerNum].infill.extend(roofInfill)
        # print layers for debugging
        printLayers = [1, 2, 3, 6, 9]
        scale = (1.0 / float(CLIPPER_SCALE_FACTOR)) * 10.0
        if False:
            for x in printLayers:
                layer = layers[x]
                svg = layer.to_svg()
                filename = 'layer-{}.svg'.format(x)
                svg.SaveToFile(filename, invScale=scale)

        logger.info('Infill time {}'.format(infillTime.secs))

        with Timer() as modelTime:
            counter = 0
            for result in layers:
                counter += 1
                zcur = result.z
                #layer = Layer()
                #layer.zlevel = zcur
                #layer.contours = result[0]
                #logger.info('{} layer.contours len {}'.format(layer.zlevel, len(layer.contours)))
                #layer.infill = result[1]
                model.layers.append(result)
                for contour in result.contours:
                    contour.zlevel = zcur
                    model.contours.append(contour)
                    if zcur != contour.zlevel:
                        raise Exception('zcur != zlevel')
                        #logger.info("z: {} zcur: {} counter: {} appending contour".format(contour.zlevel, zcur, counter))
                model.infillAtLayer[zcur] = result.infill

        logger.warn('Model setup took {0} seconds'.format(modelTime.secs))

        with Timer() as gcodeTime:
            model.writeGCode(args.output)

        logger.warn('GCode generation took {0} seconds'.format(gcodeTime.secs))

    logger.info('TOTAL TIME {} seconds'.format(totalTime.secs))
    dt1 = datetime.now()
    t1 = dt1.microsecond
    logger.info('time: %f  facets: %d  facet end: %d  zmin: %f zmax: %f ' %
                ((t1 - t0), len(model.facets), model.endLoop, model.zmin, model.zmax))
    # print('Memory usage: %s (MB)' %
    #    (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 / 1024))
