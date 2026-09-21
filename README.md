# pyslicer

pyslicer turns a 3D model file (STL) into printer instructions (G-code) for a filament 3D printer.

It reads the model, cuts it into thin layers, draws the outer walls and fill pattern, then writes a G-code file your printer can run.

![Simulated #3DBenchy print at 2×](docs/benchy_demo.gif)

## Requirements

- Python 3.10 or newer
- A working pip install

## Install

From the project folder:

```bash
pip install -e ".[dev]"
```

That installs pyslicer and the tools used to run tests.

## Use it

```bash
python -m pyslicer model.stl output.gcode
```

Or, after install:

```bash
pyslicer model.stl output.gcode
```

Replace `model.stl` with your model path and `output.gcode` with where you want the G-code saved.

### Common options

| Option | What it does | Default |
|--------|----------------|---------|
| `-l` / `--layer-height` | Thickness of each layer, in mm | `0.1` |
| `-n` / `--num-perimeters` | How many outer wall loops to print | `3` |
| `-d` / `--filament-diameter` | Filament width, in mm | `1.75` |
| `-o` / `--perimeter_overlap_percent` | How much wall loops overlap (smaller means more overlap) | `1.0` |
| `-p` / `--perimeters_only` | Walls only; skip fill | off |
| `-a` / `--append_perimeters` | Also keep the original outline from the STL | off |
| `-v` / `--verbose` | Print more detail while running | off |

Example with a thicker layer and two walls:

```bash
pyslicer model.stl output.gcode -l 0.2 -n 2
```

### YAML config files

You can put print settings in one or more YAML files and pass them with `--config`.
Files are layered left-to-right: later files override overlapping settings and add new ones.
Explicit CLI flags always win over config (including `--no-perimeters-only` /
`--no-append-perimeters` to force booleans off).

```yaml
# base.yaml
layer_height: 0.2
num_perimeters: 3
nozzle_diameter: 0.4
outer_speed: 50

# fine.yaml
outer_speed: 40
print_temperature: 205
```

```bash
pyslicer model.stl out.gcode --config base.yaml,fine.yaml -l 0.15
```

Precedence: CLI historical defaults → config layers → CLI flags.
Speeds that match CLI names (`outer_speed`, `inner_speed`, `infill_speed`,
`max_corner_speed`) are in mm/s. You may instead set `outer_perimeter_speed` /
`inner_perimeter_speed` in mm/min (not both forms for the same setting).
`infill_speed` and `max_corner_speed` have no mm/min YAML twin—they always mean mm/s.
Other Model print settings use Model attribute names and native units (for example
`print_temperature` in °C).

### Print quality settings

Defaults bias **outer walls slower than inner** (50 vs 80 mm/s) for surface finish.
Recover print time on inner walls and infill; keep the outer/inner gap modest if
pressure advance is not tuned (large flow drops can blob at the seam).

| Key | Role |
|-----|------|
| `firmware` | `none` / `klipper` / `marlin` / `rrf` — dialect for optional PA emit |
| `pressure_advance` | Firmware-native PA/LA value; emitted only when `firmware` ≠ `none` |
| `seam_position` | `aligned` / `nearest` / `rear` / `none` |
| `wipe_distance` / `wipe_on_loops` / `seam_gap` | Seam cleanup (portable `G1` moves) |
| `enable_dynamic_overhang_speeds` | Overlap % vs previous layer → lower `F` |
| `overhang_speed_0`…`_75` | Percent of outer cruise at 0/25/50/75% overlap |
| `min_layer_time` | Seconds; `0` disables cooling slowdown |
| `dont_slow_down_outer_wall` | Prefer slowing inner/infill (default true) |
| `slow_down_min_speed` | Floor in mm/s for cooling/overhang slowdowns |

Pressure advance is machine-side and **not** portable across firmwares:

- Klipper: `SET_PRESSURE_ADVANCE ADVANCE=…`
- Marlin: `M900 K…` (requires `LIN_ADVANCE` in the firmware build)
- RRF: `M572 D0 S…`

Values are not interchangeable. Default `firmware: none` emits nothing.
Do not use slicer coasting / “advanced pressure” hacks with Klipper.

Available keys and types are defined in
[`src/pyslicer/config.schema.json`](src/pyslicer/config.schema.json)
(validated automatically when you pass `--config`).
In VS Code / yaml-language-server, point a file at the schema like this:

```yaml
# yaml-language-server: $schema=./src/pyslicer/config.schema.json
layer_height: 0.2
outer_speed: 50
inner_speed: 80
```

## Preview G-code in a browser

After you have a `.gcode` file:

```bash
python -m pyslicer.preview output.gcode preview.html
```

Or while slicing, write both at once:

```bash
pyslicer model.stl output.gcode --html-preview preview.html
```

Open the HTML file in any browser. You get a rotatable 3D view of the toolpaths
(drag to spin, scroll to zoom), per-layer 2D drawings, and the full G-code listing.

## Demo: slice #3DBenchy and open the viewer

Downloads the CC0 [#3DBenchy](https://www.3dbenchy.com/) model (if needed),
slices it, writes an HTML preview, and opens it in your default browser
(macOS, Windows, and Linux):

```bash
python -m pyslicer.demo
```

Or, after install:

```bash
pyslicer-demo
```

Useful flags: `--no-open`, `--force-download`, `-l 0.4`, `-n 2`,
`--write-readme-gif` (writes `docs/benchy_demo.gif` at 2×; needs Pillow).

In the HTML viewer you can **Record movie** (WebM) or **Export GIF** of the simulation.

## Tests

```bash
pytest
```

To also check that most of the code is covered by tests:

```bash
pytest --cov=pyslicer
```


## Project layout

- `src/pyslicer/` — the slicer code (geometry, slicing, fill, G-code)
- `tests/` — unit tests, including checks that G-code numbers and commands are correct
- `pyproject.toml` — package and dependency settings

## License

MIT. See [LICENSE](LICENSE).
