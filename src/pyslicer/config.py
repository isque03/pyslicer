"""Layered YAML configuration for pyslicer print settings."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any, Mapping

import yaml
from jsonschema import Draft202012Validator

# CLI-facing mm/s keys → Model attribute (stored as mm/min).
# Values are converted (*60) during canonicalize; merged dict is Model-native.
_SPEED_MM_S_TO_ATTR = {
    "outer_speed": "outer_perimeter_speed",
    "inner_speed": "inner_perimeter_speed",
    "infill_speed": "infill_speed",
    "max_corner_speed": "max_corner_speed",
    "slow_down_min_speed": "slow_down_min_speed",
}

# Model-native mm/min speed keys that share an attribute with an mm/s alias.
_SPEED_MM_MIN_ATTRS = frozenset(
    {
        "outer_perimeter_speed",
        "inner_perimeter_speed",
    }
)

# Non-speed aliases → Model attribute (same units).
_KEY_ALIASES = {
    "layer_height": "layerHeight",
    "num_perimeters": "number_perimeters",
}

# Groups of YAML keys that target the same Model attribute (mutually exclusive).
_MUTUAL_EXCLUSION_GROUPS = (
    frozenset({"outer_speed", "outer_perimeter_speed"}),
    frozenset({"inner_speed", "inner_perimeter_speed"}),
    frozenset({"layer_height", "layerHeight"}),
    frozenset({"num_perimeters", "number_perimeters"}),
)

# Print settings that may appear in YAML (Model attribute names).
_MODEL_SETTING_KEYS = frozenset(
    {
        "nozzle_diameter",
        "layerHeight",
        "retract_amount",
        "retract_speed",
        "unretract_speed",
        "default_print_speed",
        "outer_perimeter_speed",
        "inner_perimeter_speed",
        "infill_speed",
        "default_travel_speed",
        "default_z_speed",
        "max_corner_speed",
        "max_accel",
        "max_jerk",
        "min_corner_angle",
        "filament_diameter",
        "number_perimeters",
        "processing_threads",
        "infill_density",
        "infill_angle",
        "width_over_height",
        "perimeters_only",
        "append_perimeters",
        "perimeter_overlap_percent",
        "minimum_retract_travel",
        "print_temperature",
        "simplification_factor",
        "min_contour_area",
        "min_extrude",
        "firmware",
        "pressure_advance",
        "enable_dynamic_overhang_speeds",
        "overhang_speed_0",
        "overhang_speed_25",
        "overhang_speed_50",
        "overhang_speed_75",
        "min_layer_time",
        "slow_down_min_speed",
        "dont_slow_down_outer_wall",
        "seam_position",
        "wipe_distance",
        "wipe_on_loops",
        "seam_gap",
    }
)

_ALLOWED_YAML_KEYS = (
    _MODEL_SETTING_KEYS | frozenset(_KEY_ALIASES) | frozenset(_SPEED_MM_S_TO_ATTR)
)


class ConfigError(ValueError):
    """Invalid YAML config content or unknown keys."""


def _shallow_merge(base: dict[str, Any], overlay: Mapping[str, Any]) -> dict[str, Any]:
    merged = dict(base)
    merged.update(overlay)
    return merged


@lru_cache(maxsize=1)
def load_config_schema() -> dict[str, Any]:
    """Load the packaged JSON Schema for pyslicer YAML configs."""
    text = resources.files("pyslicer").joinpath("config.schema.json").read_text(
        encoding="utf-8"
    )
    return json.loads(text)


def validate_config(data: Mapping[str, Any], *, source: str | Path | None = None) -> None:
    """Validate config data against the JSON Schema. Raises ConfigError on failure."""
    validator = Draft202012Validator(load_config_schema())
    errors = sorted(validator.iter_errors(data), key=lambda e: list(e.path))
    if not errors:
        return
    where = f" in {source}" if source is not None else ""
    messages = []
    for err in errors:
        path = ".".join(str(p) for p in err.path) or "(root)"
        messages.append(f"{path}: {err.message}")
    raise ConfigError(f"Invalid config{where}: " + "; ".join(messages))


def canonicalize_config(
    data: Mapping[str, Any], *, source: str | Path | None = None
) -> dict[str, Any]:
    """Map YAML keys to Model attributes with Model-native units.

    CLI-style speed keys (outer_speed, inner_speed, infill_speed, max_corner_speed,
    slow_down_min_speed) are mm/s and converted to mm/min.
    outer_perimeter_speed / inner_perimeter_speed are already mm/min.
    Dual keys that target the same attribute raise ConfigError.

    Note: YAML has no separate mm/min twin for infill_speed / max_corner_speed /
    slow_down_min_speed; those keys always mean mm/s (matching the CLI flags).
    """
    where = f" in {source}" if source is not None else ""
    for group in _MUTUAL_EXCLUSION_GROUPS:
        present = group & data.keys()
        if len(present) > 1:
            keys = ", ".join(sorted(present))
            raise ConfigError(
                f"Conflicting config keys{where}: {keys} set the same setting"
            )

    canonical: dict[str, Any] = {}
    for key, value in data.items():
        if key in _SPEED_MM_S_TO_ATTR:
            attr = _SPEED_MM_S_TO_ATTR[key]
            canonical[attr] = float(value) * 60.0
        elif key in _SPEED_MM_MIN_ATTRS:
            canonical[key] = float(value)
        elif key in _KEY_ALIASES:
            canonical[_KEY_ALIASES[key]] = value
        else:
            canonical[key] = value
    return canonical


def load_config_file(path: str | Path) -> dict[str, Any]:
    """Load one YAML file and return Model-native canonical settings."""
    path = Path(path)
    try:
        with path.open(encoding="utf-8") as f:
            data = yaml.safe_load(f)
    except OSError as exc:
        raise ConfigError(f"Cannot read config {path}: {exc}") from exc
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"Config root must be a mapping: {path}")
    validate_config(data, source=path)
    return canonicalize_config(data, source=path)


def load_layered_config(paths: list[str | Path]) -> dict[str, Any]:
    """Load YAML files in order; later layers replace overlapping Model settings."""
    if not paths:
        raise ConfigError("No config paths provided")
    merged: dict[str, Any] = {}
    for path in paths:
        layer = load_config_file(path)
        merged = _shallow_merge(merged, layer)
    return merged


def parse_config_paths(value: str) -> list[str]:
    """Split a comma-separated --config value into path strings."""
    paths = [p.strip() for p in value.split(",") if p.strip()]
    if not paths:
        raise ConfigError("No config paths provided")
    return paths


def apply_settings_to_model(model: Any, settings: Mapping[str, Any]) -> None:
    """Apply Model-native settings (already canonicalized) onto a Model."""
    for key, value in settings.items():
        setattr(model, key, value)


def apply_config_to_model(model: Any, config: Mapping[str, Any]) -> None:
    """Validate raw YAML keys, canonicalize to Model-native units, and apply.

    For already-canonical settings from load_layered_config, use
    apply_settings_to_model instead (avoids double-converting speeds).
    """
    validate_config(config)
    apply_settings_to_model(model, canonicalize_config(config))
