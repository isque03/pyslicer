"""Tests for the Benchy demo helper (no network in the default path)."""

from pathlib import Path

from pyslicer.demo import BENCHY_URL, build_arg_parser, default_samples_dir, open_in_browser


def test_demo_arg_parser_defaults():
    args = build_arg_parser().parse_args([])
    assert args.layer_height == 0.4
    assert args.num_perimeters == 2
    assert args.no_open is False
    assert BENCHY_URL.startswith("https://")


def test_default_samples_dir_points_at_repo_samples():
    d = default_samples_dir()
    assert d.name == "samples"
    assert (d.parent / "pyproject.toml").is_file()


def test_open_in_browser_uses_file_uri(monkeypatch, tmp_path):
    opened = {}

    def fake_open(uri):
        opened["uri"] = uri
        return True

    monkeypatch.setattr("webbrowser.open", fake_open)
    html = tmp_path / "preview.html"
    html.write_text("<html></html>", encoding="utf-8")
    open_in_browser(html)
    assert opened["uri"].startswith("file:")
    assert "preview.html" in opened["uri"]


def test_open_in_browser_raises_when_launch_fails(monkeypatch, tmp_path):
    monkeypatch.setattr("webbrowser.open", lambda uri: False)
    monkeypatch.setattr("sys.platform", "linux")

    def boom(*args, **kwargs):
        class R:
            returncode = 1

        return R()

    monkeypatch.setattr("subprocess.run", boom)
    html = tmp_path / "preview.html"
    html.write_text("<html></html>", encoding="utf-8")
    try:
        open_in_browser(html)
        assert False, "expected RuntimeError"
    except RuntimeError as exc:
        assert "Failed to open browser" in str(exc)
