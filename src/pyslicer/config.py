"""Layered YAML configuration for pyslicer print settings."""

from __future__ import annotations

import json
from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any, Mapping

import yaml
from jsonschema import Draft202012Validator

# CLI-facing aliases in mm/s → Model attribute (stored as mm/min)
_SPEED_MM_S_ALIASES = {
    "outer_speed": "outer_perimeter_speed",
    "inner_speed": "inner_perimeter_speed",
    "infill_speed": "infill_speed",
    "max_corner_speed": "max_corner_speed",
}

# Key aliases that map 1:1 without unit conversion
_KEY_ALIASES = {
    "layer_height": "layerHeight",
    "num_perimeters": "number_perimeters",
}

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
    }
)

_ALLOWED_YAML_KEYS = (
    _MODEL_SETTING_KEYS | frozenset(_KEY_ALIASES) | frozenset(_SPEED_MM_S_ALIASES)
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


def load_config_file(path: str | Path) -> dict[str, Any]:
    """Load a single YAML config file. Root must be a mapping."""
    path = Path(path)
    with path.open(encoding="utf-8") as f:
        data = yaml.safe_load(f)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"Config root must be a mapping: {path}")
    validate_config(data, source=path)
    return data


def load_layered_config(paths: list[str | Path]) -> dict[str, Any]:
    """Load YAML files in order; later layers replace overlapping keys."""
    merged: dict[str, Any] = {}
    for path in paths:
        layer = load_config_file(path)
        merged = _shallow_merge(merged, layer)
    return merged


def parse_config_paths(value: str) -> list[str]:
    """Split a comma-separated --config value into path strings."""
    return [p.strip() for p in value.split(",") if p.strip()]


def apply_config_to_model(model: Any, config: Mapping[str, Any]) -> None:
    """Apply a merged config dict onto a Model instance.

    CLI-style speed aliases (outer_speed, etc.) are mm/s and converted to
    Model mm/min. Model attribute names for speeds use Model's native mm/min.
    """
    validate_config(config)
    for key, value in config.items():
        if key in _SPEED_MM_S_ALIASES:
            attr = _SPEED_MM_S_ALIASES[key]
            setattr(model, attr, float(value) * 60.0)
            continue
        attr = _KEY_ALIASES.get(key, key)
        setattr(model, attr, value)
