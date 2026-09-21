"""G-code HTML preview: assemble 2D layers, 3D orbit view, and listing."""

from __future__ import annotations

import html
import json
from pathlib import Path

from pyslicer.preview.gcode_parse import (
    DEFAULT_NOZZLE_DIAMETER_MM,
    layers_to_toolpaths_3d,
    parse_gcode,
    parse_gcode_layers,
    parse_toolpath_moves,
)
from pyslicer.preview.viewer_js import viewer_script

# Three.js r170 ES modules (CDN)
_THREE_IMPORTMAP = """{
  "imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/"
  }
}"""


def _bounds(layers: list[dict]) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for layer in layers:
        for kind in ("extrude", "travel"):
            for a, b in layer[kind]:
                xs.extend([a[0], b[0]])
                ys.extend([a[1], b[1]])
    if not xs:
        raise ValueError("No toolpath geometry to bound for preview")
    return min(xs), max(xs), min(ys), max(ys)


def _path_d(segs: list[tuple]) -> str:
    parts = []
    for a, b in segs:
        parts.append(f"M{a[0]:.4f},{-a[1]:.4f} L{b[0]:.4f},{-b[1]:.4f}")
    return " ".join(parts)


def render_gcode_html(
    gcode: str,
    *,
    title: str = "pyslicer preview",
    subtitle: str | None = None,
    nozzle_diameter: float = DEFAULT_NOZZLE_DIAMETER_MM,
    layer_height: float | None = None,
    max_accel: float | None = None,
    planning_limits: dict | None = None,
) -> str:
    """Return a self-contained HTML document with 2D layers and a 3D orbit view."""
    from pyslicer.preview.gcode_parse import parse_planning_limits_from_gcode

    layers, moves = parse_gcode(gcode)
    limits = planning_limits or parse_planning_limits_from_gcode(gcode)
    paths_3d = layers_to_toolpaths_3d(
        layers,
        nozzle_diameter=nozzle_diameter,
        layer_height=layer_height,
        moves=moves,
        max_accel=max_accel if max_accel is not None else (
            limits.get("maxAccel") if limits else None
        ),
        planning_limits=limits,
    )
    paths_json = json.dumps(paths_3d, separators=(",", ":"))

    minx, maxx, miny, maxy = _bounds(layers)
    pad = 1.0
    vb_w = (maxx - minx) + 2 * pad
    vb_h = (maxy - miny) + 2 * pad
    vb = f"{minx - pad} {-(maxy + pad)} {vb_w} {vb_h}"

    layer_blocks = []
    drawn = 0
    for i, layer in enumerate(layers):
        if not layer["extrude"] and not layer["travel"]:
            continue
        drawn += 1
        ez = layer["z"] if layer["z"] is not None else i
        n_e = len(layer["extrude"])
        n_t = len(layer["travel"])
        layer_blocks.append(
            f"""
      <article class="layer" data-z="{ez}">
        <h2>Z {ez}</h2>
        <p class="meta">{n_e} extrude · {n_t} travel</p>
        <svg viewBox="{vb}" xmlns="http://www.w3.org/2000/svg" aria-label="Layer at Z {ez}">
          <path class="travel" d="{_path_d(layer['travel'])}" fill="none"/>
          <path class="extrude" d="{_path_d(layer['extrude'])}" fill="none"/>
        </svg>
      </article>"""
        )

    if not subtitle:
        subtitle = (
            f"{drawn} layers with toolpaths · {len(gcode.splitlines())} lines of G-code"
        )

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{html.escape(title)}</title>
<script type="importmap">
{_THREE_IMPORTMAP}
</script>
<style>
  :root {{
    --text: #1a1a1a;
    --muted: #5c5c5c;
    --extrude: #0072b2;
    --travel: #e69f00;
    --bg: #fafafa;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0;
    background: var(--bg);
    color: var(--text);
    font-family: "Source Sans 3", "Helvetica Neue", Helvetica, Arial, sans-serif;
    line-height: 1.45;
  }}
  .wrap {{
    max-width: min(100%, 1680px);
    margin: 0 auto;
    padding: 2.5rem clamp(1.25rem, 3vw, 3rem) 4rem;
  }}
  header {{ margin-bottom: 2.5rem; }}
  header h1 {{
    margin: 0 0 0.4rem;
    font-size: 1.65rem;
    font-weight: 600;
    letter-spacing: -0.02em;
  }}
  header .sub {{
    margin: 0;
    color: var(--muted);
    max-width: 40rem;
  }}
  .legend {{
    margin-top: 1rem;
    color: var(--muted);
    font-size: 0.9rem;
  }}
  .legend span + span {{ margin-left: 1.25rem; }}
  .swatch {{
    display: inline-block;
    margin-right: 0.35rem;
    vertical-align: middle;
  }}
  .swatch.extrude {{
    width: 0.85rem;
    height: 0.45rem;
    background: var(--extrude);
  }}
  .swatch.travel {{
    width: 1.1rem;
    height: 2px;
    background: var(--travel);
  }}
  .concern-legend {{
    display: none;
    margin-top: 0.75rem;
    align-items: center;
    gap: 0.5rem;
    font-size: 0.85rem;
    color: var(--muted);
  }}
  .concern-legend.visible {{
    display: flex;
  }}
  .concern-legend .bar {{
    flex: 0 0 8rem;
    height: 0.55rem;
    background: linear-gradient(90deg, #0072b2, #f0e442, #d55e00);
  }}
  .color-mode {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.55rem 0.75rem;
    font-size: 0.8rem;
    color: var(--text);
  }}
  .color-mode label {{
    display: inline-flex;
    align-items: center;
    gap: 0.25rem;
    cursor: pointer;
  }}
  .color-mode input {{
    margin: 0;
    accent-color: var(--extrude);
  }}
  #concern-thresholds[hidden] {{
    display: none !important;
  }}
  .thr-row {{
    display: flex;
    flex-direction: column;
    gap: 0.25rem;
  }}
  .thr-row .control-head {{
    font-size: 0.8rem;
  }}
  .slice-settings {{
    margin: 0;
    padding: 0;
    list-style: none;
    font-size: 0.8rem;
    color: var(--text);
    font-variant-numeric: tabular-nums;
  }}
  .slice-settings li {{
    display: flex;
    justify-content: space-between;
    gap: 0.75rem;
    padding: 0.15rem 0;
    border-bottom: 1px solid #ececec;
  }}
  .slice-settings li:last-child {{
    border-bottom: none;
  }}
  .slice-settings .k {{
    color: var(--muted);
  }}
  .slice-settings .v {{
    text-align: right;
    white-space: nowrap;
  }}
  .slice-settings-missing {{
    margin: 0;
    font-size: 0.8rem;
    color: var(--muted);
  }}
  .viewer-block {{
    margin-bottom: 3rem;
  }}
  .viewer-row {{
    display: grid;
    grid-template-columns: minmax(0, 1fr) minmax(11rem, 14rem);
    gap: clamp(1.25rem, 3vw, 2.5rem);
    align-items: start;
  }}
  @media (max-width: 720px) {{
    .viewer-row {{ grid-template-columns: 1fr; }}
  }}
  .viewer-block h2,
  .layers-heading,
  .gcode-block h2 {{
    margin: 0 0 0.75rem;
    font-size: 0.8rem;
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: var(--muted);
  }}
  .viewer-hint {{
    margin: 0 0 0.85rem;
    font-size: 0.85rem;
    color: var(--muted);
  }}
  #viewer3d {{
    width: 100%;
    min-height: 320px;
  }}
  #viewer3d canvas {{
    display: block;
    width: 100%;
    height: auto;
  }}
  .viewer-controls {{
    display: flex;
    flex-direction: column;
    gap: 1.5rem;
    padding-top: 0.15rem;
  }}
  .viewer-controls h2 {{
    margin: 0 0 0.35rem;
  }}
  .control {{
    display: flex;
    flex-direction: column;
    gap: 0.4rem;
  }}
  .control .control-head {{
    display: flex;
    justify-content: space-between;
    align-items: baseline;
    gap: 0.75rem;
    font-size: 0.85rem;
    color: var(--text);
  }}
  .control .control-meta {{
    color: var(--muted);
    font-size: 0.8rem;
    font-variant-numeric: tabular-nums;
    white-space: nowrap;
  }}
  .control .control-help {{
    margin: 0;
    font-size: 0.75rem;
    color: var(--muted);
  }}
  .control input[type="range"] {{
    width: 100%;
    margin: 0;
    accent-color: var(--extrude);
  }}
  .control.travel-accent input[type="range"] {{
    accent-color: var(--travel);
  }}
  .sim-transport {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.4rem;
  }}
  .sim-transport button {{
    font: inherit;
    font-size: 0.8rem;
    padding: 0.35rem 0.55rem;
    color: var(--text);
    background: transparent;
    border: 1px solid #c8c8c8;
    cursor: pointer;
  }}
  .sim-transport button:hover:not(:disabled) {{
    border-color: var(--text);
  }}
  .sim-transport button:disabled {{
    opacity: 0.4;
    cursor: default;
  }}
  .sim-speeds {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.55rem 0.75rem;
    font-size: 0.8rem;
    color: var(--text);
  }}
  .sim-speeds label {{
    display: inline-flex;
    align-items: center;
    gap: 0.25rem;
    cursor: pointer;
  }}
  .sim-speeds input {{
    margin: 0;
    accent-color: var(--extrude);
  }}
  .layout {{
    display: grid;
    grid-template-columns: 1fr minmax(18rem, 0.55fr);
    gap: clamp(2rem, 4vw, 4rem);
    align-items: start;
  }}
  @media (max-width: 820px) {{
    .layout {{ grid-template-columns: 1fr; gap: 2.5rem; }}
  }}
  @media (min-width: 1400px) {{
    .layers {{
      grid-template-columns: repeat(auto-fill, minmax(13rem, 1fr));
    }}
  }}
  .layers {{
    display: grid;
    grid-template-columns: repeat(auto-fill, minmax(11rem, 1fr));
    column-gap: 1.75rem;
    row-gap: 2rem;
  }}
  .layer h2 {{
    margin: 0;
    font-size: 0.8rem;
    font-weight: 600;
    letter-spacing: 0.04em;
    text-transform: uppercase;
    color: var(--muted);
  }}
  .layer .meta {{
    margin: 0.15rem 0 0.55rem;
    font-size: 0.75rem;
    color: var(--muted);
  }}
  .layer svg {{
    display: block;
    width: 100%;
    height: auto;
  }}
  .layer path.extrude {{
    stroke: var(--extrude);
    stroke-width: 0.16;
    stroke-linecap: round;
  }}
  .layer path.travel {{
    stroke: var(--travel);
    stroke-width: 0.08;
  }}
  .gcode-block pre {{
    margin: 0;
    padding: 0;
    max-height: 70vh;
    overflow: auto;
    font-family: "Source Code Pro", ui-monospace, Menlo, monospace;
    font-size: 0.72rem;
    line-height: 1.5;
    color: var(--text);
    white-space: pre;
    background: transparent;
    border: none;
  }}
</style>
</head>
<body>
  <div class="wrap">
    <header>
      <h1>{html.escape(title)}</h1>
      <p class="sub">{html.escape(subtitle)}</p>
      <p class="legend">
        <span><i class="swatch extrude"></i>Extrude</span>
        <span><i class="swatch travel"></i>Travel</span>
      </p>
      <p class="concern-legend" id="concern-legend" aria-hidden="true">
        <span>Low</span><span class="bar" role="img" aria-label="Concern gradient"></span><span>High</span>
      </p>
    </header>

    <section class="viewer-block" aria-label="3D toolpath view">
      <h2>3D view</h2>
      <p class="viewer-hint">Drag to rotate · scroll to zoom · right-drag to pan · filament uses stadium bead cross-section</p>
      <div class="viewer-row">
        <div id="viewer3d"></div>
        <aside class="viewer-controls" aria-label="Preview controls">
          <div class="control">
            <div class="control-head">
              <span>Slice settings</span>
            </div>
            <ul class="slice-settings" id="slice-settings" aria-label="Limits used to generate this G-code"></ul>
            <p class="slice-settings-missing" id="slice-settings-missing" hidden>
              No planning comment found in this G-code.
            </p>
            <p class="control-help">Values used when this file was sliced (feeds planned to stay within these).</p>
          </div>
          <div class="control">
            <div class="control-head">
              <span>Color mode</span>
            </div>
            <div class="color-mode" role="radiogroup" aria-label="Color mode">
              <label><input type="radio" name="color-mode" value="path" checked/> Path type</label>
              <label><input type="radio" name="color-mode" value="concern"/> Concern</label>
            </div>
            <p class="control-help">Concern paints the solid bead mesh where speed, accel, or corner Δv exceed your thresholds at sharp turns. Sliders recolor only — re-slice to change G-code feeds.</p>
          </div>
          <div class="control" id="concern-thresholds" hidden>
            <div class="control-head">
              <span>Concern thresholds</span>
            </div>
            <div class="thr-row">
              <div class="control-head">
                <label for="thr-speed">Max speed</label>
                <span class="control-meta"><span id="thr-speed-val">70</span> mm/s</span>
              </div>
              <input id="thr-speed" type="range" min="10" max="200" step="1" value="70"
                aria-describedby="thr-help"/>
            </div>
            <div class="thr-row">
              <div class="control-head">
                <label for="thr-accel">Max accel</label>
                <span class="control-meta"><span id="thr-accel-val">1000</span> mm/s²</span>
              </div>
              <input id="thr-accel" type="range" min="100" max="10000" step="50" value="1000"/>
            </div>
            <div class="thr-row">
              <div class="control-head">
                <label for="thr-jerk">Max corner Δv</label>
                <span class="control-meta"><span id="thr-jerk-val">20</span> mm/s</span>
              </div>
              <input id="thr-jerk" type="range" min="1" max="80" step="1" value="20"/>
            </div>
            <div class="thr-row">
              <div class="control-head">
                <label for="thr-angle">Min turn angle</label>
                <span class="control-meta"><span id="thr-angle-val">20</span>°</span>
              </div>
              <input id="thr-angle" type="range" min="0" max="90" step="1" value="20"/>
            </div>
            <p id="thr-help" class="control-help">Defaults match the limits used when this G-code was planned. At those values the mesh should stay cool if planning succeeded.</p>
          </div>
          <div class="control">
            <div class="control-head">
              <span>Simulated print</span>
              <span class="control-meta"><span id="sim-pct">100%</span> · <span id="sim-z">Z —</span></span>
            </div>
            <div class="sim-transport" role="group" aria-label="Playback">
              <button type="button" id="sim-rewind" title="Rewind to start">Rewind</button>
              <button type="button" id="sim-step-back" title="Previous G-code move">Step back</button>
              <button type="button" id="sim-play" title="Play">Play</button>
              <button type="button" id="sim-pause" title="Pause">Pause</button>
              <button type="button" id="sim-step-fwd" title="Next G-code move">Step forward</button>
              <button type="button" id="sim-ff" title="Jump ahead ~5%">Fast forward</button>
            </div>
            <p class="control-meta" id="sim-move" aria-live="polite">Move —</p>
            <input id="sim-scrub" type="range" min="0" max="100" step="0.1" value="100"
              aria-label="Playback position"/>
            <div class="sim-speeds" role="radiogroup" aria-label="Playback speed">
              <label><input type="radio" name="sim-speed" value="0.25"/> 0.25×</label>
              <label><input type="radio" name="sim-speed" value="0.5"/> 0.5×</label>
              <label><input type="radio" name="sim-speed" value="1" checked/> 1×</label>
              <label><input type="radio" name="sim-speed" value="2"/> 2×</label>
              <label><input type="radio" name="sim-speed" value="4"/> 4×</label>
            </div>
            <p class="control-help">Step forward/back one G-code move to verify the nozzle path (including corners).</p>
            <div class="sim-transport" role="group" aria-label="Export">
              <button type="button" id="export-movie-start" title="Record WebM movie">Record movie</button>
              <button type="button" id="export-movie-stop" title="Stop and download movie" disabled>Stop movie</button>
              <button type="button" id="export-gif" title="Export GIF clip at current speed">Export GIF</button>
            </div>
            <p id="export-status" class="control-help" aria-live="polite"></p>
          </div>
          <div class="control">
            <div class="control-head">
              <label for="cutaway">Print progress</label>
              <span class="control-meta"><span id="cutaway-pct">100%</span> · <span id="cutaway-z">Z —</span></span>
            </div>
            <input id="cutaway" type="range" min="0" max="100" step="1" value="100"
              aria-describedby="cutaway-help"/>
            <p id="cutaway-help" class="control-help">Cutaway to how far the print would be at this percent.</p>
          </div>
          <div class="control travel-accent">
            <div class="control-head">
              <label for="filament-opacity">Filament opacity</label>
              <span class="control-meta" id="opacity-pct">100%</span>
            </div>
            <input id="filament-opacity" type="range" min="5" max="100" step="1" value="100"
              aria-describedby="opacity-help"/>
            <p id="opacity-help" class="control-help">Opacity of extruded filament.</p>
          </div>
        </aside>
      </div>
    </section>

    <div class="layout">
      <section aria-label="Layer toolpaths">
        <h2 class="layers-heading">Layers</h2>
        <div class="layers">
          {''.join(layer_blocks) if layer_blocks else '<p class="sub">No toolpaths found.</p>'}
        </div>
      </section>
      <section class="gcode-block" aria-label="G-code listing">
        <h2>G-code</h2>
        <pre>{html.escape(gcode)}</pre>
      </section>
    </div>
  </div>
  <script>
    window.PYSLICER_TOOLPATHS = {paths_json};
  </script>
  <script type="module">
{viewer_script()}
  </script>
</body>
</html>
"""


def write_gcode_preview(
    gcode_path: str | Path,
    html_path: str | Path,
    *,
    title: str | None = None,
    subtitle: str | None = None,
    nozzle_diameter: float = DEFAULT_NOZZLE_DIAMETER_MM,
    layer_height: float | None = None,
    max_accel: float | None = None,
    planning_limits: dict | None = None,
) -> Path:
    """Read a G-code file and write an HTML preview. Returns the HTML path."""
    gcode_path = Path(gcode_path)
    html_path = Path(html_path)
    text = gcode_path.read_text(encoding="utf-8")
    doc = render_gcode_html(
        text,
        title=title or f"pyslicer — {gcode_path.name}",
        subtitle=subtitle,
        nozzle_diameter=nozzle_diameter,
        layer_height=layer_height,
        max_accel=max_accel,
        planning_limits=planning_limits,
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(doc, encoding="utf-8")
    return html_path


# Re-exports for callers that imported parsers from gcode_html
__all__ = [
    "parse_gcode",
    "parse_gcode_layers",
    "parse_toolpath_moves",
    "layers_to_toolpaths_3d",
    "render_gcode_html",
    "write_gcode_preview",
]
