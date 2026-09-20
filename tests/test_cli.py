"""CLI orchestration coverage."""

from pathlib import Path

from pyslicer.cli import build_arg_parser, main, run
from pyslicer.gcode.parse import parse_gcode
from pyslicer.timer import Timer


FIXTURE = Path(__file__).parent / "fixtures" / "cube_10mm.stl"


def test_arg_parser_defaults():
    p = build_arg_parser()
    args = p.parse_args([str(FIXTURE), "out.gcode"])
    assert args.layer_height == 0.1
    assert args.num_perimeters == 3
    assert args.filament_diameter == 1.75


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
