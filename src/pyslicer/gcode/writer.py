"""GCode emission for a sliced model."""

import logging

from pyslicer.geometry.line import Line

logger = logging.getLogger(__name__)


def retract(model, e, f):
    if model.retract_amount > 0.0:
        e -= model.retract_amount
        f.write(f"G1 F{model.retract_speed:.6f} E{e:.6f}\n")
    return e


def extrude(model, e, f):
    if model.retract_amount > 0.0:
        e += model.retract_amount
        f.write(f"G1 F{model.unretract_speed:.6f} E{e:.6f}\n")
    return e


def write_segment_gcode(model, f, idx, segment, e, extrude_amount):
    if idx == 0:
        e = retract(model, e, f)
        f.write(";; travel move\n")
        f.write(
            f"G1 X{segment.verticies[0].x:.6f} Y{segment.verticies[0].y:.6f}\n"
        )
        e = extrude(model, e, f)
        e += extrude_amount
        f.write(
            f"G1 F{model.default_print_speed:.6f} "
            f"X{segment.verticies[1].x:.6f} Y{segment.verticies[1].y:.6f} "
            f"E{e:.6f}\n"
        )
    else:
        e += extrude_amount
        f.write(
            f"G1 X{segment.verticies[1].x:.6f} Y{segment.verticies[1].y:.6f} "
            f"E{e:.6f}\n"
        )
    return e


def write_gcode(model, filename):
    logger.info("Writing %s", filename)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"M109 S{model.print_temperature} ; Heat up to {model.print_temperature}C\n")
        f.write("G90       ; Use absolute coordinates\n")
        f.write("G21       ; Set units to millimeters\n")
        f.write("M106 S0   ; Fan Off\n")
        f.write("G28       ; Home all axes\n")
        f.write("G92 E0    ; Zero extruder\n")
        f.write("M82       ; Use absolute distances for extrusion\n")
        f.write("G92 E0    ; Zero extruder\n")
        f.write("G1 F200 E8  ; prime extruder\n")
        f.write("G92 E0    ; Zero extruder\n")
        f.write("M117 Printing.\n")
        f.write("\n")
        e = 0.0000
        last_z = None
        pirsquared = model.pi_r_squared()

        for layer in model.layers:
            f.write(f";; New Layer Z: {layer.z} \n")
            contour_index = 0
            while True:
                contours_written = 0
                for x in range(len(layer.perimeters) - 1, -1, -1):
                    if len(layer.perimeters[x]) < contour_index + 1:
                        continue
                    f.write(f";; Perimeter {x}\n")
                    contours_written += 1
                    contour = layer.perimeters[x][contour_index]
                    if len(contour.segments) < 1:
                        continue
                    if last_z is None or last_z != contour.zlevel:
                        e = retract(model, e, f)
                        f.write(
                            f"G1 F{model.default_z_speed:.6f} Z{contour.zlevel:.6f}\n"
                        )
                        e = extrude(model, e, f)
                    last_z = contour.zlevel
                    f.write(
                        f";; Contour {contour_index} Area: {abs(contour.winding_area())}\n"
                    )
                    f.write(f"G1 F{model.default_print_speed:.6f}\n")
                    for idx, segment in enumerate(contour.segments):
                        volume = model.volume_extruded(segment)
                        extrude_amount = volume / pirsquared
                        e = write_segment_gcode(
                            model, f, idx, segment, e, extrude_amount
                        )
                contour_index += 1
                if contours_written == 0:
                    break

            if last_z is None:
                continue
            try:
                infill = model.infillAtLayer[last_z]
            except KeyError as exc:
                raise KeyError(
                    f"GCODE output: no infillAtLayer entry for z={last_z}"
                ) from exc
            f.write(";; Infill\n")
            prev_segment = None
            for idx, segment in enumerate(infill):
                travel_move = idx == 0
                if prev_segment is not None and (
                    prev_segment.verticies[1] != segment.verticies[0]
                ):
                    travel_move = True
                if travel_move:
                    travel_distance = 0.0
                    if prev_segment:
                        travel_distance = Line.withVerticies(
                            prev_segment.verticies[1], segment.verticies[0]
                        ).magnitude()
                    if travel_distance >= model.minimum_retract_travel:
                        f.write(";; retract \n")
                        e = retract(model, e, f)
                    f.write(
                        f";;infill travel move. distance: {travel_distance:.6f} \n"
                    )
                    f.write(
                        f"G1 F{model.default_travel_speed:.6f} "
                        f"X{segment.verticies[0].x:.6f} "
                        f"Y{segment.verticies[0].y:.6f}\n"
                    )
                    if travel_distance >= model.minimum_retract_travel:
                        e = extrude(model, e, f)
                    f.write(f"G1 F{model.default_print_speed:.6f}\n")

                volume = model.volume_extruded(segment)
                e += volume / pirsquared
                f.write(
                    f"G1 X{segment.verticies[1].x:.6f} "
                    f"Y{segment.verticies[1].y:.6f} E{e:.6f}\n"
                )
                prev_segment = segment

        f.write("M107    ; Fan off\n")
        f.write(f"M104 S0 ; Heat off {model.print_temperature}C\n")
        f.write("M140 S0 ; Bed heat off\n")
        f.write("G91     ; Relative\n")
        f.write("G1 E-1 F400 ; Retract\n")
        f.write("G1 Z+1.0 E-5 X-20 Y-20 F9000\n")
        f.write("G28 X0 Y0\n")
        f.write("M84     ; Motors off\n")
        f.write("G90     ; Absolute\n")
