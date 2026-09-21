"""Minimum layer-time cooling slowdown (prefer sparing outer walls)."""

from __future__ import annotations


def print_speed_floor_mm_min(model) -> float:
    """Never emit print F below ``Model.slow_down_min_speed`` (mm/min).

    Design rule: Model always defines this field (default 600 = 10 mm/s).
    Callers must not invent a different floor.
    """
    floor = float(model.slow_down_min_speed)
    if floor <= 0:
        raise ValueError("slow_down_min_speed must be > 0 (mm/min)")
    return floor


def path_time_seconds(segments, feed_mm_min: float) -> float:
    feed = max(float(feed_mm_min), 1e-9)
    total_len = sum(max(seg.magnitude(), 0.0) for seg in segments)
    return (total_len / feed) * 60.0


def layer_speed_scales(model, outer_time_s: float, inner_time_s: float, infill_time_s: float):
    """Return (outer_scale, inner_scale, infill_scale) in (0, 1].

    Scales multiply cruise feeds. Outer is left at 1.0 when
    ``dont_slow_down_outer_wall`` is set, until inner+infill hit the floor.
    """
    min_t = float(model.min_layer_time)
    if min_t <= 0:
        return 1.0, 1.0, 1.0

    total = outer_time_s + inner_time_s + infill_time_s
    if total >= min_t - 1e-9 or total <= 1e-12:
        return 1.0, 1.0, 1.0

    floor_f = print_speed_floor_mm_min(model)
    spare_outer = bool(model.dont_slow_down_outer_wall)

    def _floor_scale(cruise: float) -> float:
        return min(1.0, floor_f / max(float(cruise), 1e-9))

    if spare_outer:
        adjustable_t = inner_time_s + infill_time_s
        needed = min_t - outer_time_s
        if adjustable_t <= 1e-12:
            scale = min(1.0, total / min_t)
            return max(scale, _floor_scale(model.outer_perimeter_speed)), 1.0, 1.0

        # Slow adjustable so their time grows to ``needed`` (scale < 1).
        scale = adjustable_t / max(needed, 1e-12)
        scale = min(1.0, scale)
        inner_floor = _floor_scale(model.inner_perimeter_speed)
        infill_floor = _floor_scale(model.infill_speed)
        adj_floor = max(inner_floor, infill_floor)
        if scale >= adj_floor - 1e-12:
            return 1.0, max(scale, inner_floor), max(scale, infill_floor)

        # Hit the floor on internals; slow outer for any remaining deficit.
        inner_scale = inner_floor
        infill_scale = infill_floor
        new_adj = inner_time_s / max(inner_scale, 1e-12) + infill_time_s / max(
            infill_scale, 1e-12
        )
        still_need = min_t - outer_time_s - new_adj
        if still_need <= 0:
            return 1.0, inner_scale, infill_scale
        outer_scale = outer_time_s / max(outer_time_s + still_need, 1e-12)
        outer_scale = min(1.0, outer_scale)
        outer_scale = max(outer_scale, _floor_scale(model.outer_perimeter_speed))
        return outer_scale, inner_scale, infill_scale

    # Uniform cooling: scale everything
    scale = total / min_t
    scale = min(1.0, scale)
    for cruise in (
        model.outer_perimeter_speed,
        model.inner_perimeter_speed,
        model.infill_speed,
    ):
        scale = max(scale, _floor_scale(cruise))
    return scale, scale, scale


def estimate_layer_feature_times(layer, model, infill_segments) -> tuple[float, float, float]:
    """Rough time estimate using feature cruise (no corner planner)."""
    outer_t = inner_t = 0.0
    for shell_i, shells in enumerate(layer.perimeters):
        cruise = (
            float(model.outer_perimeter_speed)
            if shell_i == 0
            else float(model.inner_perimeter_speed)
        )
        for contour in shells:
            t = path_time_seconds(contour.segments, cruise)
            if shell_i == 0:
                outer_t += t
            else:
                inner_t += t
    infill_t = path_time_seconds(infill_segments or [], float(model.infill_speed))
    return outer_t, inner_t, infill_t
