"""Create a 10x10x10 mm cube as binary STL for smoke tests."""

from pathlib import Path
from struct import pack


def write_cube_stl(path: Path, size: float = 10.0) -> None:
    """Axis-aligned cube from (0,0,0) to (size,size,size)."""
    s = size
    # 12 triangles, 2 per face
    faces = [
        # bottom z=0
        ((0, 0, 0), (s, 0, 0), (s, s, 0)),
        ((0, 0, 0), (s, s, 0), (0, s, 0)),
        # top z=s
        ((0, 0, s), (s, s, s), (s, 0, s)),
        ((0, 0, s), (0, s, s), (s, s, s)),
        # front y=0
        ((0, 0, 0), (s, 0, s), (s, 0, 0)),
        ((0, 0, 0), (0, 0, s), (s, 0, s)),
        # back y=s
        ((0, s, 0), (s, s, 0), (s, s, s)),
        ((0, s, 0), (s, s, s), (0, s, s)),
        # left x=0
        ((0, 0, 0), (0, s, 0), (0, s, s)),
        ((0, 0, 0), (0, s, s), (0, 0, s)),
        # right x=s
        ((s, 0, 0), (s, 0, s), (s, s, s)),
        ((s, 0, 0), (s, s, s), (s, s, 0)),
    ]
    with open(path, "wb") as f:
        f.write(b"\0" * 80)
        f.write(pack("<I", len(faces)))
        for a, b, c in faces:
            f.write(pack("<fff", 0.0, 0.0, 0.0))  # normal
            for p in (a, b, c):
                f.write(pack("<fff", *p))
            f.write(pack("<H", 0))


if __name__ == "__main__":
    out = Path(__file__).with_name("cube_10mm.stl")
    write_cube_stl(out)
    print(f"wrote {out}")
