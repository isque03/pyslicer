"""CLI: python -m pyslicer.preview input.gcode output.html"""

import argparse
from pathlib import Path

from pyslicer.preview.gcode_html import write_gcode_preview


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description="Build a flat HTML preview of G-code toolpaths."
    )
    parser.add_argument("gcode", help="Input G-code file.")
    parser.add_argument("html", help="Output HTML path.")
    parser.add_argument(
        "-t",
        "--title",
        default=None,
        help="Page title (default: derived from the G-code filename).",
    )
    parser.add_argument(
        "-s",
        "--subtitle",
        default=None,
        help="Optional subtitle under the title.",
    )
    return parser


def main(argv=None):
    args = build_arg_parser().parse_args(argv)
    path = write_gcode_preview(
        args.gcode,
        args.html,
        title=args.title,
        subtitle=args.subtitle,
    )
    print(f"Wrote {Path(path).resolve()}")


if __name__ == "__main__":
    main()
