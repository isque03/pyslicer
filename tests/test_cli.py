"""CLI orchestration coverage."""

from pathlib import Path

import pytest

from pyslicer.cli import build_arg_parser, configure_model, main, run
from pyslicer.gcode.parse import parse_gcode
from pyslicer.timer import Timer


FIXTURE = Path(__file__).parent / "fixtures" / "cube_10mm.stl"


def test_arg_parser_defaults():
    p = build_arg_parser()
    args = p.parse_args([str(FIXTURE), "out.gcode"])
    assert not hasattr(args, "layer_height")
    assert not hasattr(args, "num_perimeters")
    assert not hasattr(args, "filament_diameter")
    assert args.verbose is False
    assert args.html_preview is None
    assert args.config is None


def test_cli_run_end_to_end(tmp_path):
    out = tmp_path / "cube.gcode"
    args = build_arg_parser().parse_args(
        [
            str(FIXTURE),
            str(out),
            "-l",
            "2.0",
            "-n",
            "1",
            "-p",
        ]
    )
    run(args)
    assert out.exists()
    text = out.read_text()
    cmds = parse_gcode(text)
    assert any(c == "G1" and "E" in p for c, p in cmds)
    assert ";; New Layer Z:" in text


def test_main_entrypoint(tmp_path):
    out = tmp_path / "cube2.gcode"
    main([str(FIXTURE), str(out), "-l", "5.0", "-n", "1"])
    assert out.exists()


def test_timer_context():
    with Timer(verbose=False) as t:
        pass
    assert t.secs >= 0.0
    assert t.msecs >= 0.0


def test_partial_config_keeps_historical_defaults(tmp_path):
    cfg = tmp_path / "temp.yaml"
    cfg.write_text("print_temperature: 210\n")
    args = build_arg_parser().parse_args(
        [str(FIXTURE), str(tmp_path / "out.gcode"), "--config", str(cfg)]
    )
    model = configure_model(args)
    assert model.print_temperature == 210
    assert model.layerHeight == 0.1
    assert model.number_perimeters == 3
    assert model.perimeters_only is False
    assert model.append_perimeters is False


def test_config_cli_wins(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("layer_height: 0.2\n")
    args = build_arg_parser().parse_args(
        [str(FIXTURE), str(tmp_path / "out.gcode"), "--config", str(cfg), "-l", "0.4"]
    )
    model = configure_model(args)
    assert model.layerHeight == 0.4


def test_no_perimeters_only_overrides_config(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("perimeters_only: true\n")
    args = build_arg_parser().parse_args(
        [
            str(FIXTURE),
            str(tmp_path / "out.gcode"),
            "--config",
            str(cfg),
            "--no-perimeters-only",
        ]
    )
    model = configure_model(args)
    assert model.perimeters_only is False


def test_config_model_only_via_cli(tmp_path):
    cfg = tmp_path / "temp.yaml"
    cfg.write_text("print_temperature: 210\nlayer_height: 2.0\nnum_perimeters: 1\n")
    out = tmp_path / "cube.gcode"
    args = build_arg_parser().parse_args(
        [str(FIXTURE), str(out), "--config", str(cfg), "-p"]
    )
    model = configure_model(args)
    assert model.print_temperature == 210
    assert model.layerHeight == 2.0
    assert model.perimeters_only is True


def test_run_with_config(tmp_path):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text("layer_height: 2.0\nnum_perimeters: 1\nperimeters_only: true\n")
    out = tmp_path / "cube.gcode"
    args = build_arg_parser().parse_args(
        [str(FIXTURE), str(out), "--config", str(cfg)]
    )
    run(args)
    assert out.exists()


def test_main_invalid_config_exits(tmp_path):
    cfg = tmp_path / "bad.yaml"
    cfg.write_text("layer_height: not_a_number\n")
    out = tmp_path / "out.gcode"
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURE), str(out), "--config", str(cfg)])
    assert exc.value.code == 2


def test_main_missing_config_exits(tmp_path):
    out = tmp_path / "out.gcode"
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURE), str(out), "--config", str(tmp_path / "missing.yaml")])
    assert exc.value.code == 2


def test_main_empty_config_list_exits(tmp_path):
    out = tmp_path / "out.gcode"
    with pytest.raises(SystemExit) as exc:
        main([str(FIXTURE), str(out), "--config", ","])
    assert exc.value.code == 2


def test_no_config_historical_defaults(tmp_path):
    args = build_arg_parser().parse_args(
        [str(FIXTURE), str(tmp_path / "out.gcode")]
    )
    model = configure_model(args)
    assert model.layerHeight == 0.1
    assert model.number_perimeters == 3
    assert model.filament_diameter == 1.75
    assert model.perimeters_only is False
    assert model.append_perimeters is False
