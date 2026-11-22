# Hosted slicing and previews

PySlicer rejects non-finite model coordinates and settings, invalid layer heights, and inputs above its processing ceilings. The limits defend the local CLI. Consult the printer profile for physical limits. The ceilings admit the bundled Benchy at 0.4 mm and 0.2 mm layer heights.

| Input | Processing ceilings |
| --- | --- |
| STL | 64 MiB, 500,000 facets, and coordinates within ±100,000 mm |
| Slicing | 10,000 layers, 250 million facet-layer checks, 64 configured perimeters, and 32 threads |
| Infill | 50,000 rows and 10 million potential intersections per infill pass |
| Config | 16 YAML files of 64 KiB each |
| G-code preview | 32 MiB, 10,000 layers, 250,000 moves, and coordinates within ±100,000 mm |

## Run a worker per job

- Pin the audited wheel and dependency versions in the worker image. Launch PySlicer as an unprivileged process with no outbound network access, a read-only application filesystem, and an isolated scratch directory. Give the process only the current user's uploaded file and use worker-generated names and trusted configuration. Pass input and output as an argument list.
- Set a wall-clock timeout, CPU and memory limits, and a disk quota. Apply upload and concurrency limits before launching a process. Terminate the job and discard partial outputs on failures or limit violations. Treat the CLI ceilings as defense in depth and tune service quotas from representative workloads.
- Generate preview data separately through `pyslicer.preview.parse_gcode` or `pyslicer.preview.layers_to_toolpaths_3d`. Validate and size-limit its JSON in the worker and render it with the host application's existing viewer. PySlicer's standalone HTML preview loads browser modules from a CDN. Serve downloaded HTML as an attachment or on a separate origin from applications that hold credentials.

## Validate for the target machine

The preview parser models a subset of motion commands for visualization. A machine-ready export needs a validated profile for the actual printer and firmware. Use a separate profile validator to review the entire G-code output, including startup and shutdown, homing and coordinate modes, extrusion, travel bounds, temperature and speed limits. Hand off validated files for operator review before printing.
