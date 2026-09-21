"""Layered YAML config loading and apply-to-model."""

from pathlib import Path

import pytest

from pyslicer import config as config_mod
from pyslicer.config import (
    ConfigError,
    apply_config_to_model,
    apply_settings_to_model,
    canonicalize_config,
    load_config_schema,
    load_layered_config,
    parse_config_paths,
    validate_config,
)
from pyslicer.mesh import Model


def test_parse_config_paths():
    assert parse_config_paths("a.yaml,b.yaml") == ["a.yaml", "b.yaml"]
    assert parse_config_paths(" a.yaml , b.yaml ") == ["a.yaml", "b.yaml"]


def test_parse_config_paths_empty_raises():
    with pytest.raises(ConfigError, match="No config paths"):
        parse_config_paths(",")
    with pytest.raises(ConfigError, match="No config paths"):
        parse_config_paths("  ")


def test_layer_merge(tmp_path):
    a = tmp_path / "a.yaml"
    b = tmp_path / "b.yaml"
    a.write_text("layer_height: 0.2\nnum_perimeters: 2\n")
    b.write_text("layer_height: 0.3\nnum_perimeters: 4\n")
    merged = load_layered_config([a, b])
    assert merged["layerHeight"] == 0.3
    assert merged["number_perimeters"] == 4


def test_cross_alias_layer_override(tmp_path):
    a = tmp_path / "a.yaml"
    b = tmp_path / "b.yaml"
    a.write_text("outer_speed: 50\n")
    b.write_text("outer_perimeter_speed: 1800\n")
    merged = load_layered_config([a, b])
    assert merged == {"outer_perimeter_speed": 1800.0}
    model = Model()
    apply_settings_to_model(model, merged)
    assert model.outer_perimeter_speed == 1800.0

    c = tmp_path / "c.yaml"
    d = tmp_path / "d.yaml"
    c.write_text("outer_perimeter_speed: 1800\n")
    d.write_text("outer_speed: 50\n")
    merged2 = load_layered_config([c, d])
    assert merged2 == {"outer_perimeter_speed": 3000.0}


def test_same_file_dual_keys_rejected(tmp_path):
    bad = tmp_path / "dual.yaml"
    bad.write_text("layer_height: 0.2\nlayerHeight: 0.3\n")
    with pytest.raises(ConfigError, match="Conflicting"):
        load_layered_config([bad])


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


def test_missing_config_file(tmp_path):
    missing = tmp_path / "nope.yaml"
    with pytest.raises(ConfigError, match="Cannot read config"):
        load_layered_config([missing])


def test_apply_speed_mm_s_alias():
    model = Model()
    apply_config_to_model(model, {"outer_speed": 70})
    assert model.outer_perimeter_speed == 4200.0


def test_infill_speed_is_mm_s():
    model = Model()
    apply_config_to_model(model, {"infill_speed": 70})
    assert model.infill_speed == 4200.0


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


def test_load_layered_config_empty_list():
    with pytest.raises(ConfigError, match="No config paths"):
        load_layered_config([])


def test_schema_matches_allowed_keys():
    schema = load_config_schema()
    assert set(schema["properties"]) == set(config_mod._ALLOWED_YAML_KEYS)
    assert schema.get("additionalProperties") is False


def test_validate_config_ok():
    validate_config({"layer_height": 0.2, "perimeters_only": True})


def test_canonicalize_outer_speed():
    assert canonicalize_config({"outer_speed": 50}) == {
        "outer_perimeter_speed": 3000.0
    }


def test_quality_config_round_trip(tmp_path):
    cfg = tmp_path / "quality.yaml"
    cfg.write_text(
        "\n".join(
            [
                "outer_speed: 50",
                "inner_speed: 80",
                "firmware: klipper",
                "pressure_advance: 0.05",
                "seam_position: aligned",
                "wipe_distance: 0.4",
                "wipe_on_loops: true",
                "seam_gap: 0.1",
                "enable_dynamic_overhang_speeds: true",
                "overhang_speed_0: 15",
                "overhang_speed_25: 25",
                "overhang_speed_50: 50",
                "overhang_speed_75: 75",
                "min_layer_time: 8",
                "slow_down_min_speed: 10",
                "dont_slow_down_outer_wall: true",
                "retract_amount: 0.8",
            ]
        )
        + "\n"
    )
    merged = load_layered_config([cfg])
    model = Model()
    apply_settings_to_model(model, merged)
    assert model.outer_perimeter_speed == 3000.0
    assert model.inner_perimeter_speed == 4800.0
    assert model.firmware == "klipper"
    assert model.pressure_advance == 0.05
    assert model.seam_position == "aligned"
    assert model.wipe_distance == 0.4
    assert model.wipe_on_loops is True
    assert model.enable_dynamic_overhang_speeds is True
    assert model.min_layer_time == 8
    assert model.slow_down_min_speed == 600.0  # 10 mm/s → mm/min
    assert model.dont_slow_down_outer_wall is True


def test_firmware_enum_rejected(tmp_path):
    bad = tmp_path / "bad_fw.yaml"
    bad.write_text("firmware: smoothie\n")
    with pytest.raises(ConfigError, match="firmware"):
        load_layered_config([bad])


def test_model_default_outer_slower_than_inner():
    m = Model()
    assert m.outer_perimeter_speed < m.inner_perimeter_speed
