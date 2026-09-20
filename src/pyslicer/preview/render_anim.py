"""Render a print-simulation GIF from G-code moves (Pillow, no browser)."""

from __future__ import annotations

import math
from pathlib import Path

from pyslicer.preview.gcode_parse import (
    DEFAULT_FEED_MM_MIN,
    DEFAULT_NOZZLE_DIAMETER_MM,
    build_timeline,
)

# Wong palette + brass nozzle (matches the HTML viewer)
_COLOR_BG = (250, 250, 250)
_COLOR_EXTRUDE = (0, 114, 178)
_COLOR_TRAVEL = (230, 159, 0)
_COLOR_NOZZLE = (230, 190, 80)
_COLOR_NOZZLE_EDGE = (180, 130, 40)
_COLOR_GRID = (220, 220, 220)

# GIF export truncates long prints so README/assets stay small (design rule).
DEFAULT_GIF_MAX_SIM_SECONDS = 24.0
DEFAULT_GIF_MAX_FRAMES = 120


def _bounds(moves: list[dict]) -> tuple[float, float, float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for m in moves:
        xs.extend((m["x0"], m["x1"]))
        ys.extend((m["y0"], m["y1"]))
        zs.extend((m["z0"], m["z1"]))
    if not xs:
        raise ValueError("Cannot render GIF: no moves in toolpath")
    return min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)


def _project(x: float, y: float, z: float) -> tuple[float, float]:
    """Isometric-ish projection (G-code Z up)."""
    ang = math.radians(35.0)
    xr = x * math.cos(ang) - y * math.sin(ang)
    yr = x * math.sin(ang) + y * math.cos(ang)
    return xr, z - 0.55 * yr


def write_simulation_gif(
    moves: list[dict],
    out_path: str | Path,
    *,
    width: int = 480,
    height: int = 300,
    fps: int = 10,
    speed: float = 2.0,
    max_sim_seconds: float | None = DEFAULT_GIF_MAX_SIM_SECONDS,
    nozzle_diameter: float = DEFAULT_NOZZLE_DIAMETER_MM,
) -> Path:
    """
    Write an animated GIF of the print simulation.

    Playback uses the shared ``build_timeline`` durations at ``speed`` (e.g. 2.0 = 2×).
    If ``max_sim_seconds`` is set, only the first that many seconds of *print*
    time are shown (truncated; status should mention this to callers).
    """
    try:
        from PIL import Image, ImageDraw
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "Pillow is required to render GIFs. Install with: pip install pillow"
        ) from exc

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    timeline, total_time = build_timeline(moves)
    end_t = total_time
    truncated = False
    if max_sim_seconds is not None:
        if total_time > float(max_sim_seconds):
            truncated = True
        end_t = min(total_time, float(max_sim_seconds))

    wall_duration = end_t / max(speed, 1e-6)
    n_frames = max(2, int(math.ceil(wall_duration * fps)))
    n_frames = min(n_frames, DEFAULT_GIF_MAX_FRAMES)

    xmin, xmax, ymin, ymax, zmin, zmax = _bounds(moves)
    corners = [
        _project(x, y, z)
        for x in (xmin, xmax)
        for y in (ymin, ymax)
        for z in (zmin, zmax)
    ]
    px = [c[0] for c in corners]
    py = [c[1] for c in corners]
    pad = 8.0
    min_px, max_px = min(px) - pad, max(px) + pad
    min_py, max_py = min(py) - pad, max(py) + pad
    span_x = max(max_px - min_px, 1e-3)
    span_y = max(max_py - min_py, 1e-3)
    scale = 0.90 * min(width / span_x, height / span_y)
    ox = width / 2 - scale * (min_px + max_px) / 2
    oy = height / 2 + scale * (min_py + max_py) / 2

    def to_screen(x: float, y: float, z: float) -> tuple[int, int]:
        u, v = _project(x, y, z)
        return int(ox + scale * u), int(oy - scale * v)

    base = Image.new("RGB", (width, height), _COLOR_BG)
    draw = ImageDraw.Draw(base)
    for g in range(0, 11):
        t = g / 10.0
        x0 = xmin + (xmax - xmin) * t
        y0 = ymin + (ymax - ymin) * t
        draw.line(
            [to_screen(x0, ymin, zmin), to_screen(x0, ymax, zmin)],
            fill=_COLOR_GRID,
            width=1,
        )
        draw.line(
            [to_screen(xmin, y0, zmin), to_screen(xmax, y0, zmin)],
            fill=_COLOR_GRID,
            width=1,
        )

    filament = base.copy()
    frames: list[Image.Image] = []
    move_i = 0
    nozzle_r = max(2, int(scale * nozzle_diameter * 0.9))

    for fi in range(n_frames):
        t = end_t * (fi / max(n_frames - 1, 1))
        fdraw = ImageDraw.Draw(filament)

        while move_i < len(timeline) and timeline[move_i]["t1"] <= t:
            m = timeline[move_i]
            a = to_screen(m["x0"], m["y0"], m["z0"])
            b = to_screen(m["x1"], m["y1"], m["z1"])
            if m["extrude"]:
                fdraw.line([a, b], fill=_COLOR_EXTRUDE, width=max(2, nozzle_r // 2))
            elif m["travelIndex"] >= 0:
                fdraw.line([a, b], fill=_COLOR_TRAVEL, width=1)
            move_i += 1

        frame = filament.copy()
        fframe = ImageDraw.Draw(frame)
        nx = ny = nz = None
        if move_i < len(timeline):
            m = timeline[move_i]
            if t < m["t0"]:
                if move_i > 0:
                    prev = timeline[move_i - 1]
                    nx, ny, nz = prev["x1"], prev["y1"], prev["z1"]
            else:
                u = (t - m["t0"]) / (m["t1"] - m["t0"])
                x = m["x0"] + (m["x1"] - m["x0"]) * u
                y = m["y0"] + (m["y1"] - m["y0"]) * u
                z = m["z0"] + (m["z1"] - m["z0"]) * u
                a = to_screen(m["x0"], m["y0"], m["z0"])
                b = to_screen(x, y, z)
                # Mid-move: show in-progress filament (matches viewer F3)
                if m["extrude"]:
                    fframe.line([a, b], fill=_COLOR_EXTRUDE, width=max(2, nozzle_r // 2))
                elif m["travelIndex"] >= 0:
                    fframe.line([a, b], fill=_COLOR_TRAVEL, width=1)
                nx, ny, nz = x, y, z
        elif timeline:
            last = timeline[-1]
            nx, ny, nz = last["x1"], last["y1"], last["z1"]

        if nx is not None:
            cx, cy = to_screen(nx, ny, nz)
            fframe.ellipse(
                [cx - nozzle_r, cy - nozzle_r, cx + nozzle_r, cy + nozzle_r],
                fill=_COLOR_NOZZLE,
                outline=_COLOR_NOZZLE_EDGE,
            )
            hi = max(1, nozzle_r // 3)
            fframe.ellipse(
                [cx - hi, cy - nozzle_r + 1, cx + 1, cy - nozzle_r + 1 + hi + 1],
                fill=(255, 240, 180),
            )
        frames.append(frame)

    duration_ms = int(1000 / max(fps, 1))
    frames[0].save(
        out_path,
        save_all=True,
        append_images=frames[1:],
        duration=duration_ms,
        loop=0,
        optimize=True,
    )
    # Stash truncation note for callers (attribute on Path is unusual; print instead)
    if truncated:
        print(
            f"GIF truncated to first {max_sim_seconds:g}s of print time "
            f"(full sim {total_time:.1f}s); feed default {DEFAULT_FEED_MM_MIN:g} mm/min"
        )
    return out_path
