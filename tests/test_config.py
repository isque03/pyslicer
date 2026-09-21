"""Layered YAML config loading and apply-to-model."""

from pathlib import Path

import pytest

from pyslicer import config as config_mod
from pyslicer.config import (
    ConfigError,
    apply_config_to_model,
    load_config_schema,
    load_layered_config,
    parse_config_paths,
    validate_config,
)
from pyslicer.mesh import Model


def test_parse_config_paths():
    assert parse_config_paths("a.yaml,b.yaml") == ["a.yaml", "b.yaml"]
    assert parse_config_paths(" a.yaml , b.yaml ") == ["a.yaml", "b.yaml"]


def test_layer_merge(tmp_path):
    a = tmp_path / "a.yaml"
    b = tmp_path / "b.yaml"
    a.write_text("layer_height: 0.2\nnum_perimeters: 2\n")
    b.write_text("layer_height: 0.3\nnum_perimeters: 4\n")
    merged = load_layered_config([a, b])
    assert merged["layer_height"] == 0.3
    assert merged["num_perimeters"] == 4


def test_unknown_key(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("not_a_real_setting: 1\n")
    with pytest.raises(ConfigError, match="Invalid config"):
        load_layered_config([bad])


def test_schema_type_error(tmp_path):
    bad = tmp_path / "types.yaml"
    bad.write_text("layer_height: not_a_number\n")
    with pytest.raises(ConfigError, match="layer_height"):
        load_layered_config([bad])


def test_schema_exclusive_minimum(tmp_path):
    bad = tmp_path / "zero.yaml"
    bad.write_text("layer_height: 0\n")
    with pytest.raises(ConfigError, match="layer_height"):
        load_layered_config([bad])


def test_non_mapping_root(tmp_path):
    bad = tmp_path / "list.yaml"
    bad.write_text("- 1\n- 2\n")
    with pytest.raises(ConfigError, match="mapping"):
        load_layered_config([bad])


def test_apply_speed_mm_s_alias():
    model = Model()
    apply_config_to_model(model, {"outer_speed": 70})
    assert model.outer_perimeter_speed == 4200.0


def test_apply_model_only_key():
    model = Model()
    apply_config_to_model(model, {"print_temperature": 210})
    assert model.print_temperature == 210


def test_apply_layer_height_alias():
    model = Model()
    apply_config_to_model(model, {"layer_height": 0.2})
    assert model.layerHeight == 0.2


def test_empty_yaml(tmp_path):
    empty = tmp_path / "empty.yaml"
    empty.write_text("")
    assert load_layered_config([empty]) == {}


def test_schema_matches_allowed_keys():
    schema = load_config_schema()
    assert set(schema["properties"]) == set(config_mod._ALLOWED_YAML_KEYS)
    assert schema.get("additionalProperties") is False


def test_validate_config_ok():
    validate_config({"layer_height": 0.2, "perimeters_only": True})
