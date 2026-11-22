"""Input and planning ceilings for local slicing and preview."""

import os
import stat
from pathlib import Path

MAX_STL_BYTES = 64 * 1024 * 1024
MAX_FACETS = 500_000
MAX_COORDINATE_MM = 100_000
MAX_FEED_MM_MIN = 60_000_000
MAX_LAYERS = 10_000
MAX_FACET_LAYERS = 250_000_000
MAX_CONFIG_BYTES = 64 * 1024
MAX_CONFIG_FILES = 16
MAX_PREVIEW_BYTES = 32 * 1024 * 1024
MAX_PREVIEW_MOVES = 250_000
MAX_INFILL_ROWS = 50_000
MAX_INFILL_CHECKS = 10_000_000


class InputError(ValueError):
    """Invalid or oversized input to the slicer or preview."""


def open_input(path, max_bytes, label):
    """Open a size-checked regular input file in binary mode."""
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_NONBLOCK", 0) | getattr(os, "O_BINARY", 0))
    stream = os.fdopen(fd, "rb")
    try:
        info = os.fstat(stream.fileno())
        if not stat.S_ISREG(info.st_mode):
            raise InputError(f"{label} must be a regular file")
        if info.st_size > max_bytes:
            raise InputError(f"{label} exceeds the {max_bytes:,}-byte limit")
    except (OSError, InputError):
        stream.close()
        raise
    return stream


def read_text_input(path, max_bytes, label):
    with open_input(path, max_bytes, label) as stream:
        data = stream.read(max_bytes + 1)
    if len(data) > max_bytes:
        raise InputError(f"{label} exceeds the {max_bytes:,}-byte limit")
    return data.decode("utf-8")


def require_distinct_paths(*paths):
    """Reject collisions between inputs and outputs, including links."""
    resolved_paths = set()
    file_ids = set()
    for path in paths:
        resolved = Path(path).resolve()
        if resolved in resolved_paths:
            raise InputError("Input and output paths must be distinct")
        resolved_paths.add(resolved)
        try:
            info = resolved.stat()
        except FileNotFoundError:
            continue
        identity = (info.st_dev, info.st_ino)
        if info.st_ino and identity in file_ids:
            raise InputError("Input and output paths must be distinct")
        file_ids.add(identity)
