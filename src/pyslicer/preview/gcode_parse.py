"""G-code parsing for HTML preview / simulation (single source of truth)."""

from __future__ import annotations

import math
import re

_G1_RE = re.compile(r"([XYZEF])(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_LAYER_RE = re.compile(r"Z:\s*([-\d.]+)")

# Design rule: when G-code omits F, simulation uses this feed (mm/min).
DEFAULT_FEED_MM_MIN = 2400.0

# Design rule: default nozzle diameter matches Model.nozzle_diameter.
DEFAULT_NOZZLE_DIAMETER_MM = 0.5


def parse_gcode(gcode: str) -> tuple[list[dict], list[dict]]:
    """
    One pass over G-code → (layers, moves).

    Layers (for 2D SVG): XY segments only, bucketed by ``;; New Layer Z:``.
    Pure-Z hops are omitted from layer paths (they have no XY stroke).

    Moves (for 3D / simulation): every G1 that changes X, Y, or Z after the
    first layer marker, including pure-Z travel.
    """
    layers: list[dict] = []
    moves: list[dict] = []
    current: dict | None = None
    x = y = z = None
    e = 0.0
    last_e = 0.0
    feed = DEFAULT_FEED_MM_MIN

    for raw in gcode.splitlines():
        if raw.startswith(";; New Layer Z:"):
            m = _LAYER_RE.search(raw)
            zval = float(m.group(1)) if m else None
            current = {"z": zval, "extrude": [], "travel": []}
            layers.append(current)
            if z is None and zval is not None:
                z = zval
            continue

        line = raw.split(";")[0].strip()
        if not line or current is None:
            continue
        if not line.upper().startswith("G1"):
            continue

        params = {k.upper(): float(v) for k, v in _G1_RE.findall(line)}
        if "F" in params:
            feed = params["F"]
        nx = params.get("X", x)
        ny = params.get("Y", y)
        nz = params.get("Z", z)
        ne = params.get("E", e)

        if x is None and nx is not None:
            x = nx
        if y is None and ny is not None:
            y = ny
        if z is None and nz is not None:
            z = nz

        xy_moved = (
            x is not None
            and y is not None
            and nx is not None
            and ny is not None
            and ("X" in params or "Y" in params)
            and (nx != x or ny != y)
        )
        any_moved = (
            x is not None
            and y is not None
            and z is not None
            and nx is not None
            and ny is not None
            and nz is not None
            and ("X" in params or "Y" in params or "Z" in params)
            and (nx != x or ny != y or nz != z)
        )
        extruding = "E" in params and ne > last_e

        if any_moved:
            moves.append(
                {
                    "x0": x,
                    "y0": y,
                    "z0": z,
                    "x1": nx,
                    "y1": ny,
                    "z1": nz,
                    "extrude": extruding,
                    "feed": feed,
                }
            )

        if xy_moved:
            seg = ((x, y), (nx, ny))
            if extruding:
                current["extrude"].append(seg)
            else:
                current["travel"].append(seg)

        x = nx if nx is not None else x
        y = ny if ny is not None else y
        z = nz if nz is not None else z
        e = ne if ne is not None else e
        if "E" in params:
            last_e = e

    return layers, moves


def parse_gcode_layers(gcode: str) -> list[dict]:
    """Back-compat: layers only."""
    layers, _ = parse_gcode(gcode)
    return layers


def parse_toolpath_moves(gcode: str) -> list[dict]:
    """Back-compat: chronological moves only."""
    _, moves = parse_gcode(gcode)
    return moves


def build_timeline(moves: list[dict]) -> tuple[list[dict], float]:
    """
    Attach t0/t1 (seconds at 1×) and extrude/travel indices for simulation.

    Duration = path_length_mm / (feed_mm_min / 60).
    """
    timeline: list[dict] = []
    acc = 0.0
    extrude_idx = 0
    travel_idx = 0
    for m in moves:
        dx = m["x1"] - m["x0"]
        dy = m["y1"] - m["y0"]
        dz = m["z1"] - m["z0"]
        length = math.sqrt(dx * dx + dy * dy + dz * dz)
        feed = max(float(m.get("feed") or DEFAULT_FEED_MM_MIN), 1.0)
        dur = max(length / (feed / 60.0), 1e-4)
        xy_travel = (not m["extrude"]) and (dx != 0.0 or dy != 0.0)
        entry = {
            "x0": m["x0"],
            "y0": m["y0"],
            "z0": m["z0"],
            "x1": m["x1"],
            "y1": m["y1"],
            "z1": m["z1"],
            "extrude": bool(m["extrude"]),
            "feed": feed,
            "t0": acc,
            "t1": acc + dur,
            "extrudeIndex": extrude_idx if m["extrude"] else -1,
            "travelIndex": travel_idx if xy_travel else -1,
        }
        if entry["extrudeIndex"] >= 0:
            extrude_idx += 1
        if entry["travelIndex"] >= 0:
            travel_idx += 1
        timeline.append(entry)
        acc += dur
    return timeline, max(acc, 1e-3)


def infer_layer_height(moves: list[dict], layers: list[dict] | None = None) -> float | None:
    """Most common positive Z step between successive layer planes (mm)."""
    zs: list[float] = []
    if layers:
        for layer in layers:
            if layer.get("z") is not None:
                zs.append(float(layer["z"]))
    if not zs and moves:
        zs = sorted({float(m["z0"]) for m in moves} | {float(m["z1"]) for m in moves})
    else:
        zs = sorted(set(zs))
    deltas: list[float] = []
    for i in range(len(zs) - 1):
        d = zs[i + 1] - zs[i]
        if d > 1e-6:
            deltas.append(round(d, 6))
    if not deltas:
        return None
    # Mode (most frequent delta)
    counts: dict[float, int] = {}
    for d in deltas:
        counts[d] = counts.get(d, 0) + 1
    return max(counts.items(), key=lambda kv: kv[1])[0]


def layers_to_toolpaths_3d(
    layers: list[dict],
    *,
    nozzle_diameter: float = DEFAULT_NOZZLE_DIAMETER_MM,
    layer_height: float | None = None,
    moves: list[dict] | None = None,
    timeline: list[dict] | None = None,
    total_time: float | None = None,
) -> dict:
    """Build the JSON payload embedded in the HTML viewer."""
    if moves is not None:
        extrude = [
            [m["x0"], m["y0"], m["z0"], m["x1"], m["y1"], m["z1"]]
            for m in moves
            if m["extrude"]
        ]
        # Policy: pure-Z hops stay in moves/timeline but not in travel[] lines
        travel = [
            [m["x0"], m["y0"], m["z0"], m["x1"], m["y1"], m["z1"]]
            for m in moves
            if not m["extrude"] and (m["x0"] != m["x1"] or m["y0"] != m["y1"])
        ]
        if timeline is None:
            timeline, total_time = build_timeline(moves)
    else:
        extrude = []
        travel = []
        for i, layer in enumerate(layers):
            z = float(layer["z"]) if layer["z"] is not None else float(i)
            for a, b in layer["extrude"]:
                extrude.append([a[0], a[1], z, b[0], b[1], z])
            for a, b in layer["travel"]:
                travel.append([a[0], a[1], z, b[0], b[1], z])

    inferred = layer_height
    if inferred is None:
        inferred = infer_layer_height(moves or [], layers)

    from pyslicer.preview.bead_profile import bead_cross_section

    layer_h = float(inferred) if inferred is not None else float(nozzle_diameter)
    bead = bead_cross_section(float(nozzle_diameter), layer_h, on_bed=False)
    bead_bed = bead_cross_section(float(nozzle_diameter), layer_h, on_bed=True)

    payload: dict = {
        "extrude": extrude,
        "travel": travel,
        "nozzleDiameter": float(nozzle_diameter),
        "layerHeight": layer_h,
        "tubeDiameter": bead["width"],
        "bead": bead,
        "beadBed": bead_bed,
    }
    if moves is not None:
        payload["moves"] = moves
    if timeline is not None:
        payload["timeline"] = timeline
        payload["totalTime"] = float(total_time if total_time is not None else 0.0)
    return payload
