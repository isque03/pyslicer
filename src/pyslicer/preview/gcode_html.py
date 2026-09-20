"""G-code HTML preview: 2D layers, 3D orbit view, and G-code listing."""

from __future__ import annotations

import html
import json
import re
from pathlib import Path

_G1_RE = re.compile(r"([XYZEF])(-?\d+(?:\.\d+)?)", re.IGNORECASE)
_LAYER_RE = re.compile(r"Z:\s*([-\d.]+)")

# Three.js r170 ES modules (CDN)
_THREE_IMPORTMAP = """{
  "imports": {
    "three": "https://cdn.jsdelivr.net/npm/three@0.170.0/build/three.module.js",
    "three/addons/": "https://cdn.jsdelivr.net/npm/three@0.170.0/examples/jsm/"
  }
}"""


def parse_gcode_layers(gcode: str) -> list[dict]:
    """
    Split G-code into layers with extrude and travel segments.

    Each layer: {z, extrude: [((x0,y0),(x1,y1)), ...], travel: [...]}
    """
    layers: list[dict] = []
    current = None
    x = y = z = None
    e = 0.0
    last_e = 0.0

    for raw in gcode.splitlines():
        if raw.startswith(";; New Layer Z:"):
            m = _LAYER_RE.search(raw)
            zval = float(m.group(1)) if m else None
            current = {"z": zval, "extrude": [], "travel": []}
            layers.append(current)
            continue

        line = raw.split(";")[0].strip()
        if not line or current is None:
            continue
        if not line.upper().startswith("G1"):
            continue

        params = {k.upper(): float(v) for k, v in _G1_RE.findall(line)}
        nx = params.get("X", x)
        ny = params.get("Y", y)
        nz = params.get("Z", z)
        ne = params.get("E", e)

        if (
            nx is not None
            and ny is not None
            and x is not None
            and y is not None
            and ("X" in params or "Y" in params)
        ):
            seg = ((x, y), (nx, ny))
            if "E" in params and ne > last_e:
                current["extrude"].append(seg)
            else:
                current["travel"].append(seg)

        x = nx if nx is not None else x
        y = ny if ny is not None else y
        z = nz if nz is not None else z
        e = ne if ne is not None else e
        if "E" in params:
            last_e = e

    return layers


def layers_to_toolpaths_3d(layers: list[dict]) -> dict:
    """Flatten layers into 3D polylines for the orbit viewer (mm coords)."""
    extrude: list[list[float]] = []
    travel: list[list[float]] = []
    for i, layer in enumerate(layers):
        z = float(layer["z"]) if layer["z"] is not None else float(i)
        for a, b in layer["extrude"]:
            extrude.append([a[0], a[1], z, b[0], b[1], z])
        for a, b in layer["travel"]:
            travel.append([a[0], a[1], z, b[0], b[1], z])
    return {"extrude": extrude, "travel": travel}


def _bounds(layers: list[dict]) -> tuple[float, float, float, float]:
    xs: list[float] = []
    ys: list[float] = []
    for layer in layers:
        for kind in ("extrude", "travel"):
            for a, b in layer[kind]:
                xs.extend([a[0], b[0]])
                ys.extend([a[1], b[1]])
    if not xs:
        return 0.0, 10.0, 0.0, 10.0
    return min(xs), max(xs), min(ys), max(ys)


def _path_d(segs: list[tuple]) -> str:
    parts = []
    for a, b in segs:
        parts.append(f"M{a[0]:.4f},{-a[1]:.4f} L{b[0]:.4f},{-b[1]:.4f}")
    return " ".join(parts)


def _viewer_script() -> str:
    return r"""
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const data = window.PYSLICER_TOOLPATHS;
const mount = document.getElementById('viewer3d');
if (!mount || !data) throw new Error('Missing 3D mount or toolpath data');

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xfafafa);

const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 5000);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
mount.appendChild(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.rotateSpeed = 0.85;

scene.add(new THREE.AmbientLight(0xffffff, 0.85));
const key = new THREE.DirectionalLight(0xffffff, 0.45);
key.position.set(2, 4, 3);
scene.add(key);

function addSegments(segments, color, linewidth) {
  if (!segments.length) return null;
  const positions = new Float32Array(segments.length * 6);
  for (let i = 0; i < segments.length; i++) {
    const s = segments[i];
    // G-code Z-up → Three.js Y-up
    const o = i * 6;
    positions[o] = s[0]; positions[o+1] = s[2]; positions[o+2] = -s[1];
    positions[o+3] = s[3]; positions[o+4] = s[5]; positions[o+5] = -s[4];
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  const mat = new THREE.LineBasicMaterial({ color, linewidth });
  const lines = new THREE.LineSegments(geo, mat);
  scene.add(lines);
  return lines;
}

addSegments(data.travel || [], 0x9a9a9a, 1);
addSegments(data.extrude || [], 0x0a7a4b, 2);

const box = new THREE.Box3().setFromObject(scene);
const size = box.getSize(new THREE.Vector3());
const center = box.getCenter(new THREE.Vector3());
const span = Math.max(size.x, size.y, size.z, 1);
controls.target.copy(center);
camera.position.set(center.x + span * 1.4, center.y + span * 1.1, center.z + span * 1.4);
camera.near = span / 200;
camera.far = span * 40;
camera.updateProjectionMatrix();
controls.update();

const grid = new THREE.GridHelper(span * 2.2, 16, 0xd0d0d0, 0xe8e8e8);
grid.position.set(center.x, box.min.y - 0.01, center.z);
scene.add(grid);

function resize() {
  const w = mount.clientWidth;
  const h = Math.max(320, Math.round(w * 0.62));
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
resize();
window.addEventListener('resize', resize);

(function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
})();
"""


def render_gcode_html(
    gcode: str,
    *,
    title: str = "pyslicer preview",
    subtitle: str | None = None,
) -> str:
    """Return a self-contained HTML document with 2D layers and a 3D orbit view."""
    layers = parse_gcode_layers(gcode)
    paths_3d = layers_to_toolpaths_3d(layers)
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
      <article class="layer">
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
    --extrude: #0a7a4b;
    --travel: #9a9a9a;
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
    width: 1.1rem;
    height: 2px;
    margin-right: 0.35rem;
    vertical-align: middle;
  }}
  .swatch.extrude {{ background: var(--extrude); }}
  .swatch.travel {{ background: var(--travel); }}
  .viewer-block {{
    margin-bottom: 3rem;
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
    </header>

    <section class="viewer-block" aria-label="3D toolpath view">
      <h2>3D view</h2>
      <p class="viewer-hint">Drag to rotate · scroll to zoom · right-drag to pan</p>
      <div id="viewer3d"></div>
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
{_viewer_script()}
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
) -> Path:
    """Read a G-code file and write an HTML preview. Returns the HTML path."""
    gcode_path = Path(gcode_path)
    html_path = Path(html_path)
    text = gcode_path.read_text(encoding="utf-8")
    doc = render_gcode_html(
        text,
        title=title or f"pyslicer — {gcode_path.name}",
        subtitle=subtitle,
    )
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(doc, encoding="utf-8")
    return html_path
