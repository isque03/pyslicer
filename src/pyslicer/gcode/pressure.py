"""Firmware pressure / linear advance start G-code."""

from __future__ import annotations


def pressure_advance_gcode(firmware: str, value: float | None) -> str | None:
    """Return a single PA/LA command line, or None if nothing should be emitted."""
    if value is None:
        return None
    fw = (firmware or "none").lower()
    if fw == "none":
        return None
    v = float(value)
    if fw == "klipper":
        return f"SET_PRESSURE_ADVANCE ADVANCE={v:.6f}"
    if fw == "marlin":
        return f"M900 K{v:.6f}"
    if fw == "rrf":
        return f"M572 D0 S{v:.6f}"
    return None
