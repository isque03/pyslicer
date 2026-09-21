"""G-code parsing for HTML preview / simulation (single source of truth)."""

from __future__ import annotations

import math
import re

_G1_RE = re.compile(r"([XYZEF])(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_LAYER_RE = re.compile(r"Z:\s*([-\d.]+)")
_PLANNING_RE = re.compile(
    r";\s*pyslicer planning:\s*(.+)",
    re.IGNORECASE,
)
_PLANNING_KV_RE = re.compile(
    r"(outer|inner|infill|accel|jerk|corner|min_angle)\s*=\s*([-\d.]+)",
    re.IGNORECASE,
)

# Design rule: when G-code omits F, simulation uses this feed (mm/min).
DEFAULT_FEED_MM_MIN = 2400.0

# Design rule: default nozzle diameter matches Model.nozzle_diameter.
DEFAULT_NOZZLE_DIAMETER_MM = 0.5


def parse_planning_limits_from_gcode(gcode: str) -> dict | None:
    """Parse ``; pyslicer planning: ...`` comment into user-facing mm/s limits."""
    for raw in gcode.splitlines():
        m = _PLANNING_RE.search(raw)
        if not m:
            continue
        vals = {k.lower(): float(v) for k, v in _PLANNING_KV_RE.findall(m.group(1))}
        if not vals:
            return None
        outer = vals.get("outer")
        return {
            "maxSpeed": outer if outer is not None else 70.0,
            "outerSpeed": outer if outer is not None else 70.0,
            "innerSpeed": vals.get("inner", outer if outer is not None else 70.0),
            "infillSpeed": vals.get("infill", outer if outer is not None else 70.0),
            "maxAccel": vals.get("accel", 1000.0),
            "maxJerk": vals.get("jerk", 20.0),
            "maxCornerSpeed": vals.get("corner", 5.0),
            "minAngleDeg": vals.get("min_angle", 20.0),
        }
    return None


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


def _pts_close(a: list[float], b: list[float], eps: float = 1e-6) -> bool:
    return (
        abs(a[0] - b[0]) <= eps
        and abs(a[1] - b[1]) <= eps
        and abs(a[2] - b[2]) <= eps
    )


def chain_extrude_polylines(
    moves: list[dict] | None = None,
    *,
    segments: list[list[float]] | None = None,
) -> list[dict]:
    """
    Chain tip-connected extrude moves into polylines for CAD-style sweeps.

    Each polyline is ``{"points": [[x,y,z], ...], "i0": int, "i1": int}`` where
    ``i0``/``i1`` are inclusive indices into the chronological extrude list
    (same order as ``extrude`` in the viewer payload). Travel / disconnected
    tips break the chain.
    """
    if moves is not None:
        segs: list[list[float]] = [
            [m["x0"], m["y0"], m["z0"], m["x1"], m["y1"], m["z1"]]
            for m in moves
            if m["extrude"]
        ]
    elif segments is not None:
        segs = segments
    else:
        return []

    polylines: list[dict] = []
    pts: list[list[float]] | None = None
    i0 = 0
    for i, s in enumerate(segs):
        p0 = [float(s[0]), float(s[1]), float(s[2])]
        p1 = [float(s[3]), float(s[4]), float(s[5])]
        if pts is None:
            pts = [p0, p1]
            i0 = i
        elif _pts_close(pts[-1], p0):
            pts.append(p1)
        else:
            polylines.append({"points": pts, "i0": i0, "i1": i - 1})
            pts = [p0, p1]
            i0 = i
    if pts is not None:
        polylines.append({"points": pts, "i0": i0, "i1": len(segs) - 1})
    return polylines


def corner_concern_metrics(moves: list[dict]) -> dict:
    """Per-extrude-segment junction metrics for interactive concern coloring.

    For each extrude segment, records the worst connected junction:
    - ``feed``: mm/min modal G-code F for that extrude segment (from ``G1 F`` /
      carry-forward), never from planning-comment outer/inner/infill speeds.
    - ``turnRad``: max turn angle at either endpoint (0 = straight)
    - ``accelRaw``: ``v² · (1 - cos θ)`` with ``v`` in mm/s (centripetal proxy)
    - ``jerkRaw``: ``v · 2 · sin(θ/2)`` in mm/s (instantaneous Δv proxy)
    """
    extrude_moves = [m for m in moves if m.get("extrude")]
    n = len(extrude_moves)
    feeds = [
        max(float(m.get("feed") or DEFAULT_FEED_MM_MIN), 1.0) for m in extrude_moves
    ]
    turn_rad = [0.0] * n
    accel_raw = [0.0] * n
    jerk_raw = [0.0] * n
    if n == 0:
        return {
            "feed": feeds,
            "turnRad": turn_rad,
            "accelRaw": accel_raw,
            "jerkRaw": jerk_raw,
        }

    def _feed_mm_s(m: dict) -> float:
        return max(float(m.get("feed") or DEFAULT_FEED_MM_MIN), 1.0) / 60.0

    def _dir(m: dict) -> tuple[float, float, float] | None:
        dx = m["x1"] - m["x0"]
        dy = m["y1"] - m["y0"]
        dz = m["z1"] - m["z0"]
        mag = math.sqrt(dx * dx + dy * dy + dz * dz)
        if mag < 1e-12:
            return None
        return (dx / mag, dy / mag, dz / mag)

    def _share_vertex(a: dict, b: dict) -> bool:
        return (
            abs(a["x1"] - b["x0"]) <= 1e-6
            and abs(a["y1"] - b["y0"]) <= 1e-6
            and abs(a["z1"] - b["z0"]) <= 1e-6
        )

    def _apply_junction(i_a: int, i_b: int, theta: float, v: float) -> None:
        cos = math.cos(theta)
        accel = (v * v) * (1.0 - cos)
        jerk = v * 2.0 * math.sin(theta / 2.0)
        for i in (i_a, i_b):
            if theta > turn_rad[i]:
                turn_rad[i] = theta
            if accel > accel_raw[i]:
                accel_raw[i] = accel
            if jerk > jerk_raw[i]:
                jerk_raw[i] = jerk

    for i in range(n - 1):
        a = extrude_moves[i]
        b = extrude_moves[i + 1]
        if not _share_vertex(a, b):
            continue
        ua = _dir(a)
        ub = _dir(b)
        if ua is None or ub is None:
            continue
        dot = max(-1.0, min(1.0, ua[0] * ub[0] + ua[1] * ub[1] + ua[2] * ub[2]))
        theta = math.acos(dot)
        v = max(_feed_mm_s(a), _feed_mm_s(b))
        _apply_junction(i, i + 1, theta, v)

    return {
        "feed": feeds,
        "turnRad": turn_rad,
        "accelRaw": accel_raw,
        "jerkRaw": jerk_raw,
    }


def scores_from_concern_metrics(
    metrics: dict,
    *,
    max_speed: float = 70.0,
    max_accel: float = 1000.0,
    max_jerk: float = 20.0,
    min_angle_deg: float = 20.0,
) -> list[float]:
    """Map raw metrics to [0, 1] using user thresholds (sharp corners only).

    ``max_speed`` is mm/s (print-speed convention). Segment feeds in ``metrics``
    are mm/min (G-code F).
    """
    feeds = metrics.get("feed") or []
    n = len(feeds)
    if n == 0:
        return []
    turn = metrics["turnRad"]
    accel = metrics["accelRaw"]
    jerk = metrics["jerkRaw"]
    min_angle = math.radians(max(float(min_angle_deg), 0.0))
    max_speed = max(float(max_speed), 1e-9)
    max_accel = max(float(max_accel), 1e-9)
    max_jerk = max(float(max_jerk), 1e-9)
    out: list[float] = []
    for i in range(n):
        if turn[i] < min_angle:
            out.append(0.0)
            continue
        speed_part = (feeds[i] / 60.0) / max_speed
        accel_part = accel[i] / max_accel
        jerk_part = jerk[i] / max_jerk
        # 0 while within thresholds; ramps only when limits are exceeded
        excess = max(speed_part, accel_part, jerk_part)
        out.append(min(1.0, max(0.0, excess - 1.0)))
    return out


def corner_concern_scores(
    moves: list[dict],
    *,
    max_accel: float | None = None,
    max_speed: float = 70.0,
    max_jerk: float = 20.0,
    min_angle_deg: float = 20.0,
) -> list[float]:
    """Per-extrude-segment concern in [0, 1] from feed and turn severity.

    ``max_speed`` is mm/s.
    """
    metrics = corner_concern_metrics(moves)
    if not metrics["feed"]:
        return []
    if max_accel is None:
        # Fall back to 95th-percentile accel so old call sites stay useful
        positive = sorted(s for s in metrics["accelRaw"] if s > 0)
        if not positive:
            return [0.0] * len(metrics["feed"])
        idx = min(len(positive) - 1, int(0.95 * (len(positive) - 1)))
        max_accel = max(positive[idx], 1e-9)
    return scores_from_concern_metrics(
        metrics,
        max_speed=max_speed,
        max_accel=max_accel,
        max_jerk=max_jerk,
        min_angle_deg=min_angle_deg,
    )


def layers_to_toolpaths_3d(
    layers: list[dict],
    *,
    nozzle_diameter: float = DEFAULT_NOZZLE_DIAMETER_MM,
    layer_height: float | None = None,
    moves: list[dict] | None = None,
    timeline: list[dict] | None = None,
    total_time: float | None = None,
    max_accel: float | None = None,
    planning_limits: dict | None = None,
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
        polylines = chain_extrude_polylines(moves)
        metrics = corner_concern_metrics(moves)
        concern = corner_concern_scores(moves, max_accel=max_accel)
    else:
        extrude = []
        travel = []
        for i, layer in enumerate(layers):
            z = float(layer["z"]) if layer["z"] is not None else float(i)
            for a, b in layer["extrude"]:
                extrude.append([a[0], a[1], z, b[0], b[1], z])
            for a, b in layer["travel"]:
                travel.append([a[0], a[1], z, b[0], b[1], z])
        polylines = chain_extrude_polylines(segments=extrude)
        metrics = {
            "feed": [DEFAULT_FEED_MM_MIN] * len(extrude),
            "turnRad": [0.0] * len(extrude),
            "accelRaw": [0.0] * len(extrude),
            "jerkRaw": [0.0] * len(extrude),
        }
        concern = [0.0] * len(extrude)

    inferred = layer_height
    if inferred is None:
        inferred = infer_layer_height(moves or [], layers)

    from pyslicer.preview.bead_profile import bead_cross_section

    layer_h = float(inferred) if inferred is not None else float(nozzle_diameter)
    bead = bead_cross_section(float(nozzle_diameter), layer_h, on_bed=False)
    bead_bed = bead_cross_section(float(nozzle_diameter), layer_h, on_bed=True)

    payload: dict = {
        "extrude": extrude,
        "extrudePolylines": polylines,
        "travel": travel,
        "concern": concern,
        "concernMetrics": metrics,
        "nozzleDiameter": float(nozzle_diameter),
        "layerHeight": layer_h,
        "tubeDiameter": bead["width"],
        "bead": bead,
        "beadBed": bead_bed,
    }
    if planning_limits:
        payload["planningLimits"] = {
            "maxSpeed": float(planning_limits.get("maxSpeed", 70)),  # mm/s
            "outerSpeed": float(
                planning_limits.get("outerSpeed", planning_limits.get("maxSpeed", 70))
            ),
            "innerSpeed": float(
                planning_limits.get("innerSpeed", planning_limits.get("maxSpeed", 70))
            ),
            "infillSpeed": float(
                planning_limits.get("infillSpeed", planning_limits.get("maxSpeed", 70))
            ),
            "maxAccel": float(planning_limits.get("maxAccel", 1000)),
            "maxJerk": float(planning_limits.get("maxJerk", 20)),
            "minAngleDeg": float(planning_limits.get("minAngleDeg", 20)),
            "maxCornerSpeed": float(planning_limits.get("maxCornerSpeed", 5)),
        }
    if moves is not None:
        payload["moves"] = moves
    if timeline is not None:
        payload["timeline"] = timeline
        payload["totalTime"] = float(total_time if total_time is not None else 0.0)
    return payload
