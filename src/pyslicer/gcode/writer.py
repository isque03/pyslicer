"""GCode emission for a sliced model."""

import logging

from pyslicer.gcode.motion import plan_contour_feeds
from pyslicer.geometry.line import Line

logger = logging.getLogger(__name__)

# Emit a new F when planned feed differs by more than this (mm/min)
_FEED_EPS = 0.5


def print_speed_for_shell(model, shell_index: int) -> float:
    """Outer shell (index 0) vs inner shells."""
    if shell_index == 0:
        return float(model.outer_perimeter_speed)
    return float(model.inner_perimeter_speed)


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


def write_segment_gcode(model, f, idx, segment, e, extrude_amount, feed, last_feed=None):
    """Write one extrusion segment. ``feed`` is mm/min. Returns (e, last_feed)."""
    feed = float(feed)
    if idx == 0:
        e = retract(model, e, f)
        f.write(";; travel move\n")
        f.write(
            f"G1 X{segment.verticies[0].x:.6f} Y{segment.verticies[0].y:.6f}\n"
        )
        e = extrude(model, e, f)
        e += extrude_amount
        f.write(
            f"G1 F{feed:.6f} "
            f"X{segment.verticies[1].x:.6f} Y{segment.verticies[1].y:.6f} "
            f"E{e:.6f}\n"
        )
        return e, feed

    e += extrude_amount
    if last_feed is None or abs(feed - last_feed) > _FEED_EPS:
        f.write(
            f"G1 F{feed:.6f} "
            f"X{segment.verticies[1].x:.6f} Y{segment.verticies[1].y:.6f} "
            f"E{e:.6f}\n"
        )
        return e, feed
    f.write(
        f"G1 X{segment.verticies[1].x:.6f} Y{segment.verticies[1].y:.6f} "
        f"E{e:.6f}\n"
    )
    return e, last_feed if last_feed is not None else feed


def _plan_feeds(model, segments, cruise_f, closed=False):
    return plan_contour_feeds(
        segments,
        cruise_f,
        model.max_corner_speed,
        model.max_accel,
        max_jerk=getattr(model, "max_jerk", 20.0),
        min_corner_angle_deg=getattr(model, "min_corner_angle", 20.0),
        closed=closed,
    )


def write_gcode(model, filename):
    logger.info("Writing %s", filename)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(f"M109 S{model.print_temperature} ; Heat up to {model.print_temperature}C\n")
        f.write("G90       ; Use absolute coordinates\n")
        f.write("G21       ; Set units to millimeters\n")
        f.write(
            f"; pyslicer planning: "
            f"outer={model.outer_perimeter_speed / 60.0:.1f}mm/s "
            f"inner={model.inner_perimeter_speed / 60.0:.1f}mm/s "
            f"infill={model.infill_speed / 60.0:.1f}mm/s "
            f"accel={model.max_accel:.0f}mm/s^2 "
            f"jerk={model.max_jerk:.1f}mm/s "
            f"corner={model.max_corner_speed / 60.0:.1f}mm/s "
            f"min_angle={model.min_corner_angle:.0f}deg\n"
        )
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
                    cruise = print_speed_for_shell(model, x)
                    feeds = _plan_feeds(
                        model,
                        contour.segments,
                        cruise,
                        closed=contour.is_closed(),
                    )
                    last_feed = None
                    f.write(f"G1 F{cruise:.6f}\n")
                    last_feed = cruise
                    for idx, segment in enumerate(contour.segments):
                        volume = model.volume_extruded(segment)
                        extrude_amount = volume / pirsquared
                        seg_feed = feeds[idx] if idx < len(feeds) else cruise
                        e, last_feed = write_segment_gcode(
                            model,
                            f,
                            idx,
                            segment,
                            e,
                            extrude_amount,
                            seg_feed,
                            last_feed=last_feed,
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
            run_start = 0
            cruise = float(model.infill_speed)
            last_feed = None

            def _flush_infill_run(start, end, e_val, last_f):
                """Plan and write contiguous infill segments [start, end)."""
                if start >= end:
                    return e_val, last_f
                run = infill[start:end]
                feeds = _plan_feeds(model, run, cruise, closed=False)
                lf = last_f
                for local_i, segment in enumerate(run):
                    volume = model.volume_extruded(segment)
                    extrude_amount = volume / pirsquared
                    seg_feed = feeds[local_i] if local_i < len(feeds) else cruise
                    # After a travel we already set F; treat as continuation (idx>0)
                    # unless this is the first extrude after travel with no prior F.
                    e_val += extrude_amount
                    if lf is None or abs(seg_feed - lf) > _FEED_EPS:
                        f.write(
                            f"G1 F{seg_feed:.6f} "
                            f"X{segment.verticies[1].x:.6f} "
                            f"Y{segment.verticies[1].y:.6f} E{e_val:.6f}\n"
                        )
                        lf = seg_feed
                    else:
                        f.write(
                            f"G1 X{segment.verticies[1].x:.6f} "
                            f"Y{segment.verticies[1].y:.6f} E{e_val:.6f}\n"
                        )
                return e_val, lf

            for idx, segment in enumerate(infill):
                travel_move = idx == 0
                if prev_segment is not None and (
                    prev_segment.verticies[1] != segment.verticies[0]
                ):
                    travel_move = True
                if travel_move:
                    # Finish previous run before traveling
                    e, last_feed = _flush_infill_run(run_start, idx, e, last_feed)
                    run_start = idx
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
                    f.write(f"G1 F{cruise:.6f}\n")
                    last_feed = cruise
                prev_segment = segment

            e, last_feed = _flush_infill_run(run_start, len(infill), e, last_feed)

        f.write("M107    ; Fan off\n")
        f.write(f"M104 S0 ; Heat off {model.print_temperature}C\n")
        f.write("M140 S0 ; Bed heat off\n")
        f.write("G91     ; Relative\n")
        f.write("G1 E-1 F400 ; Retract\n")
        f.write("G1 Z+1.0 E-5 X-20 Y-20 F9000\n")
        f.write("G28 X0 Y0\n")
        f.write("M84     ; Motors off\n")
        f.write("G90     ; Absolute\n")
