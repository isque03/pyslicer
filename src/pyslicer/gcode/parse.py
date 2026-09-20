"""Lightweight GCode line tokenizer for tests and validation."""

import re

_TOKEN_RE = re.compile(
    r"^([GMT]\d+)(?:\s+(.*))?$", re.IGNORECASE
)
_PARAM_RE = re.compile(r"([A-Za-z])(-?\d+(?:\.\d+)?)")


def parse_line(line):
    """Parse one GCode line into (command, params_dict) or None for comments/blank."""
    raw = line.strip()
    if not raw or raw.startswith(";") or raw.startswith("("):
        return None
    # Strip trailing comment
    if ";" in raw:
        raw = raw.split(";", 1)[0].strip()
    if not raw:
        return None
    m = _TOKEN_RE.match(raw)
    if not m:
        return None
    cmd = m.group(1).upper()
    params = {}
    rest = m.group(2) or ""
    for pm in _PARAM_RE.finditer(rest):
        params[pm.group(1).upper()] = float(pm.group(2))
    return cmd, params


def parse_gcode(text):
    """Return list of (command, params) for all motion/command lines."""
    result = []
    for line in text.splitlines():
        parsed = parse_line(line)
        if parsed is not None:
            result.append(parsed)
    return result
