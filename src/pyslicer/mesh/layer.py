"""Layer result of slicing at a Z height."""

import random


class Layer:
    def __init__(self):
        self.contours = None
        self.z = None
        self.infill = None
        self.perimeters = []
        self.roofs = []
        self.overhang = []
        self.normalInfill = []

    def to_svg(self):  # pragma: no cover — debug dump only
        with open(f"layer-{self.z}.svg", "w", encoding="utf-8") as f:
            f.write('<?xml version="1.0" standalone="no"?>\n')
            f.write(
                '<svg width="1000mm" height="1000mm" viewBox="-50 -50 100 100" '
                'version="1.1" xmlns="http://www.w3.org/2000/svg">\n'
            )
            for x in range(len(self.perimeters)):
                r = random.randint(0, 255)
                g = random.randint(0, 255)
                b = random.randint(0, 255)
                for contour in self.perimeters[x]:
                    f.write('<path d="')
                    for path in contour.segments:
                        f.write(
                            f" M{path.verticies[0].x} {path.verticies[0].y}"
                        )
                        f.write(
                            f" L{path.verticies[1].x} {path.verticies[1].y}"
                        )
                    f.write(
                        f'" stroke="#{r:02x}{g:02X}{b:02X}" '
                        f'fill="transparent" stroke-width=".1"/>\n'
                    )
            f.write("</svg>")
