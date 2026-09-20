"""Download #3DBenchy, slice it, and open the HTML preview in a browser.

Works on macOS, Windows, and Linux (stdlib webbrowser + urllib).
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import urllib.request
from pathlib import Path

BENCHY_URL = (
    "https://raw.githubusercontent.com/CreativeTools/3DBenchy/"
    "master/Single-part/3DBenchy.stl"
)
BENCHY_LICENSE = """#3DBenchy by Creative Tools — https://www.3dbenchy.com/
License: CC0 1.0 Universal (public domain dedication)
Source: https://github.com/CreativeTools/3DBenchy
"""


def default_samples_dir() -> Path:
    """Prefer the repo ``samples/`` folder when running from a source checkout."""
    here = Path(__file__).resolve()
    # src/pyslicer/demo.py → repo root is parents[2]
    if len(here.parents) >= 3:
        repo = here.parents[2]
        if (repo / "pyproject.toml").is_file():
            return repo / "samples"
    return Path.cwd() / "samples"


def download_benchy(dest: Path, *, force: bool = False) -> Path:
    """Download the Benchy STL (and a short license note) if needed."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    license_path = dest.with_name("3DBenchy.LICENSE.txt")
    if not license_path.exists() or force:
        license_path.write_text(BENCHY_LICENSE, encoding="utf-8")

    if dest.exists() and dest.stat().st_size > 0 and not force:
        print(f"Using existing {dest}")
        return dest

    print(f"Downloading #3DBenchy from\n  {BENCHY_URL}")
    print(f"  → {dest}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        urllib.request.urlretrieve(BENCHY_URL, tmp)
        tmp.replace(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
        raise
    print(f"Downloaded {dest.stat().st_size:,} bytes")
    return dest


def slice_benchy(
    stl: Path,
    gcode: Path,
    html: Path,
    *,
    layer_height: float,
    num_perimeters: int,
) -> None:
    """Run the slicer CLI and write G-code + HTML preview."""
    gcode.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "-m",
        "pyslicer",
        str(stl),
        str(gcode),
        "-l",
        str(layer_height),
        "-n",
        str(num_perimeters),
        "--html-preview",
        str(html),
    ]
    print("Slicing:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def open_in_browser(path: Path) -> None:
    """Open a local file in the default browser (macOS / Windows / Linux).

    Raises RuntimeError if no browser could be launched.
    """
    import webbrowser

    uri = path.resolve().as_uri()
    print(f"Opening {uri}")
    if webbrowser.open(uri):
        return
    system = sys.platform
    if system == "darwin":
        result = subprocess.run(["open", str(path.resolve())], check=False)
        if result.returncode == 0:
            return
    elif system.startswith("win"):
        import os

        try:
            os.startfile(str(path.resolve()))  # type: ignore[attr-defined]
            return
        except OSError as exc:
            raise RuntimeError(f"Failed to open browser for {path}") from exc
    else:
        result = subprocess.run(["xdg-open", str(path.resolve())], check=False)
        if result.returncode == 0:
            return
    raise RuntimeError(f"Failed to open browser for {path}")


def write_readme_gif(
    gcode_path: Path,
    gif_path: Path,
    *,
    speed: float = 2.0,
    max_sim_seconds: float = 24.0,
    nozzle_diameter: float = 0.5,
) -> Path:
    """Render a 2× simulation GIF from G-code for the project README."""
    from pyslicer.preview.gcode_html import parse_toolpath_moves
    from pyslicer.preview.render_anim import write_simulation_gif

    text = gcode_path.read_text(encoding="utf-8")
    moves = parse_toolpath_moves(text)
    print(f"Rendering README GIF ({len(moves)} moves, {speed}×) → {gif_path}")
    return write_simulation_gif(
        moves,
        gif_path,
        speed=speed,
        max_sim_seconds=max_sim_seconds,
        nozzle_diameter=nozzle_diameter,
        width=480,
        height=300,
        fps=10,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Download #3DBenchy (CC0), slice it with pyslicer, "
            "and open the HTML toolpath viewer."
        )
    )
    p.add_argument(
        "--samples-dir",
        type=Path,
        default=None,
        help="Directory for STL / G-code / preview (default: repo samples/ or ./samples).",
    )
    p.add_argument(
        "-l",
        "--layer-height",
        type=float,
        default=0.4,
        help="Layer height in mm (default: 0.4).",
    )
    p.add_argument(
        "-n",
        "--num-perimeters",
        type=int,
        default=2,
        help="Number of perimeters (default: 2).",
    )
    p.add_argument(
        "--force-download",
        action="store_true",
        help="Re-download the STL even if it already exists.",
    )
    p.add_argument(
        "--no-open",
        action="store_true",
        help="Do not open a browser after writing the preview.",
    )
    p.add_argument(
        "--skip-slice",
        action="store_true",
        help="Skip slicing; only download (and open existing preview if present).",
    )
    p.add_argument(
        "--write-readme-gif",
        type=Path,
        nargs="?",
        const=Path("docs/benchy_demo.gif"),
        default=None,
        help=(
            "Also write a 2× simulation GIF (default path: docs/benchy_demo.gif). "
            "Requires Pillow."
        ),
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    samples = args.samples_dir or default_samples_dir()
    samples = samples.expanduser().resolve()

    stl = samples / "3DBenchy.stl"
    gcode = samples / "3DBenchy.gcode"
    html = samples / "3DBenchy_preview.html"

    try:
        download_benchy(stl, force=args.force_download)
    except Exception as exc:
        print(f"Download failed: {exc}", file=sys.stderr)
        return 1

    if not args.skip_slice:
        try:
            slice_benchy(
                stl,
                gcode,
                html,
                layer_height=args.layer_height,
                num_perimeters=args.num_perimeters,
            )
        except subprocess.CalledProcessError as exc:
            print(f"Slice failed with exit code {exc.returncode}", file=sys.stderr)
            return exc.returncode or 1

    if not html.is_file():
        print(f"Preview not found: {html}", file=sys.stderr)
        return 1

    if args.write_readme_gif is not None:
        gif_path = args.write_readme_gif
        if not gif_path.is_absolute():
            # Resolve relative to repo root when possible
            root = samples.parent if samples.name == "samples" else Path.cwd()
            gif_path = (root / gif_path).resolve()
        if not gcode.is_file():
            print(f"G-code not found for GIF: {gcode}", file=sys.stderr)
            return 1
        try:
            write_readme_gif(gcode, gif_path, speed=2.0, max_sim_seconds=24.0)
            print(f"Wrote {gif_path}")
        except ImportError as exc:
            print(exc, file=sys.stderr)
            return 1

    if not args.no_open:
        try:
            open_in_browser(html)
        except RuntimeError as exc:
            print(exc, file=sys.stderr)
            return 1
    else:
        print(f"Preview written to {html}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
