"""Three.js viewer script embedded in the HTML preview."""

def viewer_script() -> str:
    return r"""
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const data = window.PYSLICER_TOOLPATHS;
const mount = document.getElementById('viewer3d');
if (!mount || !data) throw new Error('Missing 3D mount or toolpath data');

// Wong palette: blue extrude / yellow travel (color-blind friendly)
const COLOR_EXTRUDE = 0x0072b2;
const COLOR_TRAVEL = 0xe69f00;
// Bright polished brass nozzle
const COLOR_NOZZLE = 0xe1b84a;
const COLOR_NOZZLE_TIP = 0xffd978;

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xfafafa);

const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 5000);
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
renderer.localClippingEnabled = true;
mount.appendChild(renderer.domElement);

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.rotateSpeed = 0.85;

scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const key = new THREE.DirectionalLight(0xfff2d9, 0.95);
key.position.set(2, 4, 3);
scene.add(key);
const fill = new THREE.DirectionalLight(0xddeeff, 0.35);
fill.position.set(-2, 1, -1);
scene.add(fill);
const rim = new THREE.DirectionalLight(0xffffff, 0.4);
rim.position.set(0, 3, -4);
scene.add(rim);

// Clip away geometry above print-progress height (Three Y = G-code Z)
const clipPlane = new THREE.Plane(new THREE.Vector3(0, -1, 0), 1e6);

function toThree(x, y, z) {
  // G-code Z-up → Three.js Y-up
  return new THREE.Vector3(x, z, -y);
}

function zRange(segmentsList) {
  let zMin = Infinity;
  let zMax = -Infinity;
  for (const segments of segmentsList) {
    for (const s of segments) {
      zMin = Math.min(zMin, s[2], s[5]);
      zMax = Math.max(zMax, s[2], s[5]);
    }
  }
  if (!Number.isFinite(zMin)) {
    zMin = 0;
    zMax = 1;
  }
  if (zMax <= zMin) zMax = zMin + 1;
  return { zMin, zMax };
}

function addTravelLines(segments, color) {
  if (!segments.length) return null;
  const positions = new Float32Array(segments.length * 6);
  for (let i = 0; i < segments.length; i++) {
    const s = segments[i];
    const o = i * 6;
    const a = toThree(s[0], s[1], s[2]);
    const b = toThree(s[3], s[4], s[5]);
    positions[o] = a.x; positions[o+1] = a.y; positions[o+2] = a.z;
    positions[o+3] = b.x; positions[o+4] = b.y; positions[o+5] = b.z;
  }
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geo.setDrawRange(0, segments.length * 2);
  const mat = new THREE.LineBasicMaterial({
    color,
    clippingPlanes: [clipPlane],
  });
  const lines = new THREE.LineSegments(geo, mat);
  scene.add(lines);
  return { lines, geometry: geo, total: segments.length };
}

function makeStadiumBeadGeometry(bead) {
  // Stadium / capsule cross-section (real FDM bead):
  //   flat top  ← nozzle face
  //   flat bottom ← bed / prior layer
  //   semicircle left & right ← plastic bulge
  //
  // Shape coords: x = sideways (width), y = up (layer height).
  // Extrude along Z (= path length), then rotate so length is local +X.
  const h = Math.max(Number(bead.height) || 0.4, 0.05);
  const w = Math.max(Number(bead.width) || h, h);
  const r = h / 2;
  const halfFlat = Math.max(w - h, 0) / 2;
  const arcSegs = 20;

  const pts = [];
  // Bottom flat (left → right)
  pts.push(new THREE.Vector2(-halfFlat, -r));
  pts.push(new THREE.Vector2(halfFlat, -r));
  // Right semicircle (bottom → top)
  for (let i = 1; i <= arcSegs; i++) {
    const a = -Math.PI / 2 + (Math.PI * i) / arcSegs;
    pts.push(new THREE.Vector2(halfFlat + Math.cos(a) * r, Math.sin(a) * r));
  }
  // Top flat (right → left)
  pts.push(new THREE.Vector2(-halfFlat, r));
  // Left semicircle (top → bottom)
  for (let i = 1; i <= arcSegs; i++) {
    const a = Math.PI / 2 + (Math.PI * i) / arcSegs;
    pts.push(new THREE.Vector2(-halfFlat + Math.cos(a) * r, Math.sin(a) * r));
  }

  const shape = new THREE.Shape(pts);
  const geo = new THREE.ExtrudeGeometry(shape, {
    depth: 1,
    bevelEnabled: false,
    steps: 1,
    curveSegments: arcSegs,
  });
  geo.translate(0, 0, -0.5);
  geo.rotateY(-Math.PI / 2);
  geo.computeVertexNormals();
  return geo;
}

function addExtrudeBeads(segments, color, bead) {
  if (!segments.length) return null;
  const profile = bead || { height: 0.4, width: 0.5 };
  const geo = makeStadiumBeadGeometry(profile);
  const mat = new THREE.MeshStandardMaterial({
    color,
    roughness: 0.55,
    metalness: 0.02,
    flatShading: false,
    transparent: true,
    opacity: 1,
    depthWrite: true,
    clippingPlanes: [clipPlane],
  });
  const mesh = new THREE.InstancedMesh(geo, mat, segments.length);
  mesh.count = segments.length;
  const dummy = new THREE.Object3D();
  const xAxis = new THREE.Vector3(1, 0, 0);
  // Slight vertical overlap so stacked flats don't show light leaks
  const yScale = 1.02;

  for (let i = 0; i < segments.length; i++) {
    const s = segments[i];
    const a = toThree(s[0], s[1], s[2]);
    const b = toThree(s[3], s[4], s[5]);
    const mid = a.clone().add(b).multiplyScalar(0.5);
    const dir = new THREE.Vector3(b.x - a.x, 0, b.z - a.z);
    const len = dir.length();
    const zSpan = Math.abs(b.y - a.y);
    if (len < 1e-8 && zSpan < 1e-8) {
      dummy.scale.set(0, 0, 0);
      dummy.quaternion.identity();
    } else if (len < 1e-8) {
      // Rare pure-Z extrude: stand the bead on end
      dummy.position.copy(mid);
      dummy.scale.set(Math.max(zSpan, 0.05), yScale, 1);
      dummy.quaternion.setFromAxisAngle(new THREE.Vector3(0, 0, 1), Math.PI / 2);
    } else {
      dummy.position.copy(mid);
      // Scale only along path (local X). Profile Y/Z already in mm.
      dummy.scale.set(len, yScale, 1);
      dummy.quaternion.setFromUnitVectors(xAxis, dir.normalize());
    }
    dummy.updateMatrix();
    mesh.setMatrixAt(i, dummy.matrix);
  }
  mesh.instanceMatrix.needsUpdate = true;
  scene.add(mesh);
  return { mesh, material: mat, total: segments.length };
}

function makeNozzle(nozzleDiameter) {
  const d = Math.max(Number(nozzleDiameter) || 0.5, 0.2);
  const group = new THREE.Group();
  // Bright brass — avoid metalness≈1 without an env map (renders black)
  const bodyMat = new THREE.MeshStandardMaterial({
    color: COLOR_NOZZLE,
    roughness: 0.35,
    metalness: 0.55,
    emissive: 0x3a2a08,
    emissiveIntensity: 0.25,
  });
  const tipMat = new THREE.MeshStandardMaterial({
    color: COLOR_NOZZLE_TIP,
    roughness: 0.28,
    metalness: 0.6,
    emissive: 0x4a3500,
    emissiveIntensity: 0.3,
  });
  const body = new THREE.Mesh(
    new THREE.CylinderGeometry(d * 1.1, d * 0.55, d * 3.2, 16),
    bodyMat
  );
  body.position.y = d * 2.4;
  const tip = new THREE.Mesh(new THREE.ConeGeometry(d * 0.55, d * 1.4, 16), tipMat);
  tip.rotation.x = Math.PI;
  tip.position.y = d * 0.7;
  group.add(body);
  group.add(tip);
  group.visible = false;
  scene.add(group);
  return group;
}

const travelObj = addTravelLines(data.travel || [], COLOR_TRAVEL);
const extrudeObj = addExtrudeBeads(
  data.extrude || [],
  COLOR_EXTRUDE,
  data.bead || {
    height: data.layerHeight || 0.4,
    width: data.tubeDiameter || data.nozzleDiameter || 0.5,
  }
);
const nozzle = makeNozzle(data.nozzleDiameter);
const { zMin, zMax } = zRange([data.extrude || [], data.travel || []]);

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

// --- Print simulation timeline (precomputed in Python; do not recompute here) ---
const timeline = data.timeline || [];
const totalTime = Math.max(Number(data.totalTime) || 0, 1e-3);

const sim = {
  playing: false,
  speed: 1,
  t: totalTime,
  mode: 'cutaway', // 'cutaway' | 'simulate'
};

function setFilamentReveal(extrudeCount, travelCount) {
  if (extrudeObj) {
    extrudeObj.mesh.count = Math.max(0, Math.min(extrudeObj.total, extrudeCount));
  }
  if (travelObj) {
    const n = Math.max(0, Math.min(travelObj.total, travelCount));
    travelObj.geometry.setDrawRange(0, n * 2);
  }
}

function showAllFilament() {
  if (extrudeObj) extrudeObj.mesh.count = extrudeObj.total;
  if (travelObj) travelObj.geometry.setDrawRange(0, travelObj.total * 2);
}

function setZCutaway(pct) {
  const t = Math.max(0, Math.min(100, Number(pct))) / 100;
  let cutoff = zMin + (zMax - zMin) * t;
  if (t <= 0) cutoff = zMin - 0.01;
  if (t >= 1) cutoff = zMax + span;
  clipPlane.constant = cutoff;

  const zLabel = document.getElementById('cutaway-z');
  const pctLabel = document.getElementById('cutaway-pct');
  if (pctLabel) pctLabel.textContent = `${Math.round(t * 100)}%`;
  if (zLabel) {
    const shownZ = t <= 0 ? zMin : Math.min(zMax, Math.max(zMin, zMin + (zMax - zMin) * t));
    zLabel.textContent = `Z ${shownZ.toFixed(2)}`;
  }

  document.querySelectorAll('.layer[data-z]').forEach((el) => {
    const z = Number(el.dataset.z);
    el.hidden = Number.isFinite(z) && z > cutoff + 1e-6;
  });
}

function setPrintProgress(pct) {
  sim.mode = 'cutaway';
  sim.playing = false;
  updatePlayButtons();
  nozzle.visible = false;
  showAllFilament();
  setZCutaway(pct);
}

function setFilamentOpacity(pct) {
  const opacity = Math.max(0.05, Math.min(100, Number(pct))) / 100;
  if (extrudeObj && extrudeObj.material) {
    extrudeObj.material.opacity = opacity;
    extrudeObj.material.depthWrite = opacity > 0.95;
    extrudeObj.material.needsUpdate = true;
  }
  const label = document.getElementById('opacity-pct');
  if (label) label.textContent = `${Math.round(opacity * 100)}%`;
}

function seekSimulation(timeSec) {
  sim.mode = 'simulate';
  sim.t = Math.max(0, Math.min(totalTime, timeSec));
  clipPlane.constant = zMax + span;

  let eCount = 0;
  let tCount = 0;
  let pos = null;
  let curZ = zMin;

  if (!timeline.length) {
    nozzle.visible = false;
    setFilamentReveal(0, 0);
    updateSimLabels();
    return;
  }

  pos = toThree(timeline[0].x0, timeline[0].y0, timeline[0].z0);
  for (const m of timeline) {
    if (sim.t >= m.t1) {
      if (m.extrudeIndex >= 0) eCount = m.extrudeIndex + 1;
      if (m.travelIndex >= 0) tCount = m.travelIndex + 1;
      pos = toThree(m.x1, m.y1, m.z1);
      curZ = m.z1;
    } else if (sim.t >= m.t0) {
      const u = (sim.t - m.t0) / (m.t1 - m.t0);
      const x = m.x0 + (m.x1 - m.x0) * u;
      const y = m.y0 + (m.y1 - m.y0) * u;
      const z = m.z0 + (m.z1 - m.z0) * u;
      pos = toThree(x, y, z);
      curZ = z;
      // Reveal the in-progress segment so filament tracks the nozzle (not one behind)
      if (m.extrudeIndex >= 0) eCount = m.extrudeIndex + 1;
      if (m.travelIndex >= 0) tCount = m.travelIndex + 1;
      break;
    } else {
      break;
    }
  }

  setFilamentReveal(eCount, tCount);
  if (pos) {
    nozzle.position.copy(pos);
    nozzle.visible = true;
  }

  const pathPct = (sim.t / totalTime) * 100;
  const zPct = ((curZ - zMin) / (zMax - zMin)) * 100;
  const cutawayInput = document.getElementById('cutaway');
  if (cutawayInput) cutawayInput.value = String(Math.round(Math.max(0, Math.min(100, zPct))));
  if (simScrub) simScrub.value = String(pathPct);
  // Path-based filament reveal; hide 2D layers above the nozzle
  clipPlane.constant = zMax + span;
  document.querySelectorAll('.layer[data-z]').forEach((el) => {
    const z = Number(el.dataset.z);
    el.hidden = Number.isFinite(z) && z > curZ + 1e-6;
  });
  const zLabel = document.getElementById('cutaway-z');
  const pctLabel = document.getElementById('cutaway-pct');
  if (pctLabel) pctLabel.textContent = `${Math.round(Math.max(0, Math.min(100, zPct)))}%`;
  if (zLabel) zLabel.textContent = `Z ${curZ.toFixed(2)}`;

  updateSimLabels(pathPct, curZ);
}

function updateSimLabels(pathPct, curZ) {
  const pctEl = document.getElementById('sim-pct');
  const zEl = document.getElementById('sim-z');
  if (pctEl) {
    const p = pathPct != null ? pathPct : (sim.t / totalTime) * 100;
    pctEl.textContent = `${Math.round(p)}%`;
  }
  if (zEl && curZ != null) zEl.textContent = `Z ${curZ.toFixed(2)}`;
}

function updatePlayButtons() {
  const playBtn = document.getElementById('sim-play');
  const pauseBtn = document.getElementById('sim-pause');
  if (playBtn) playBtn.disabled = sim.playing;
  if (pauseBtn) pauseBtn.disabled = !sim.playing;
}

function playSim() {
  if (sim.t >= totalTime - 1e-6) sim.t = 0;
  sim.playing = true;
  sim.mode = 'simulate';
  updatePlayButtons();
  seekSimulation(sim.t);
}

function pauseSim() {
  sim.playing = false;
  updatePlayButtons();
}

function rewindSim() {
  sim.playing = false;
  updatePlayButtons();
  seekSimulation(0);
}

function fastForwardSim() {
  sim.mode = 'simulate';
  seekSimulation(Math.min(totalTime, sim.t + Math.max(totalTime * 0.05, 2)));
}

const cutawayInput = document.getElementById('cutaway');
const opacityInput = document.getElementById('filament-opacity');
if (cutawayInput) {
  cutawayInput.addEventListener('input', () => setPrintProgress(cutawayInput.value));
}
if (opacityInput) {
  opacityInput.addEventListener('input', () => setFilamentOpacity(opacityInput.value));
  setFilamentOpacity(opacityInput.value);
}

document.getElementById('sim-play')?.addEventListener('click', playSim);
document.getElementById('sim-pause')?.addEventListener('click', pauseSim);
document.getElementById('sim-rewind')?.addEventListener('click', rewindSim);
document.getElementById('sim-ff')?.addEventListener('click', fastForwardSim);

document.querySelectorAll('input[name="sim-speed"]').forEach((el) => {
  el.addEventListener('change', () => {
    if (el.checked) sim.speed = Number(el.value) || 1;
  });
});
const checkedSpeed = document.querySelector('input[name="sim-speed"]:checked');
if (checkedSpeed) sim.speed = Number(checkedSpeed.value) || 1;

const simScrub = document.getElementById('sim-scrub');
if (simScrub) {
  simScrub.addEventListener('input', () => {
    sim.playing = false;
    updatePlayButtons();
    seekSimulation((Number(simScrub.value) / 100) * totalTime);
  });
}

// Start in full cutaway view (complete print)
setPrintProgress(cutawayInput ? cutawayInput.value : 100);
updatePlayButtons();
updateSimLabels(100, zMax);

function downloadBlob(blob, filename) {
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 2000);
}

function setExportStatus(text) {
  const el = document.getElementById('export-status');
  if (el) el.textContent = text || '';
}

const exportState = {
  recording: false,
  recorder: null,
  chunks: [],
};

function pickRecorderMime() {
  const candidates = [
    'video/webm;codecs=vp9',
    'video/webm;codecs=vp8',
    'video/webm',
  ];
  for (const m of candidates) {
    if (window.MediaRecorder && MediaRecorder.isTypeSupported(m)) return m;
  }
  return '';
}

function startMovieExport() {
  if (exportState.recording) return;
  const mime = pickRecorderMime();
  if (!mime) {
    setExportStatus('Movie export needs MediaRecorder (try Chrome/Firefox/Edge).');
    return;
  }
  const canvas = renderer.domElement;
  const stream = canvas.captureStream(30);
  exportState.chunks = [];
  const rec = new MediaRecorder(stream, { mimeType: mime, videoBitsPerSecond: 4_000_000 });
  rec.ondataavailable = (e) => {
    if (e.data && e.data.size) exportState.chunks.push(e.data);
  };
  rec.onstop = () => {
    const blob = new Blob(exportState.chunks, { type: mime.split(';')[0] });
    downloadBlob(blob, 'pyslicer-print.webm');
    setExportStatus(`Saved movie (${Math.round(blob.size / 1024)} KB)`);
    exportState.recording = false;
    updateExportButtons();
  };
  exportState.recorder = rec;
  exportState.recording = true;
  updateExportButtons();
  setExportStatus('Recording… press Stop movie when done');
  rewindSim();
  playSim();
  rec.start(200);
}

function stopMovieExport() {
  if (!exportState.recording || !exportState.recorder) return;
  pauseSim();
  if (exportState.recorder.state !== 'inactive') exportState.recorder.stop();
}

function updateExportButtons() {
  const startBtn = document.getElementById('export-movie-start');
  const stopBtn = document.getElementById('export-movie-stop');
  if (startBtn) startBtn.disabled = exportState.recording;
  if (stopBtn) stopBtn.disabled = !exportState.recording;
}

async function exportGifClip() {
  const btn = document.getElementById('export-gif');
  if (btn) btn.disabled = true;
  pauseSim();
  try {
    const { GIFEncoder, quantize, applyPalette } = await import(
      'https://cdn.jsdelivr.net/npm/gifenc@1.0.3/dist/gifenc.esm.js'
    );
    const speed = sim.speed || 2;
    const maxSimSec = Math.min(totalTime, 24);
    setExportStatus(
      `Rendering GIF… (truncated to first ${maxSimSec.toFixed(1)}s of print @ ${speed}×)`
    );
    const fps = 10;
    const wallSec = maxSimSec / speed;
    const nFrames = Math.min(120, Math.max(2, Math.ceil(wallSec * fps)));
    const w = renderer.domElement.width;
    const h = renderer.domElement.height;
    const gif = GIFEncoder();
    const delay = Math.round(1000 / fps);

    for (let i = 0; i < nFrames; i++) {
      const t = maxSimSec * (i / Math.max(nFrames - 1, 1));
      seekSimulation(t);
      controls.update();
      renderer.render(scene, camera);
      const tmp = document.createElement('canvas');
      tmp.width = w;
      tmp.height = h;
      const tctx = tmp.getContext('2d');
      tctx.drawImage(renderer.domElement, 0, 0);
      const imageData = tctx.getImageData(0, 0, w, h);
      const palette = quantize(imageData.data, 256);
      const index = applyPalette(imageData.data, palette);
      gif.writeFrame(index, w, h, { palette, delay });
      setExportStatus(`Rendering GIF… ${i + 1}/${nFrames} (max ${maxSimSec.toFixed(1)}s @ ${speed}×)`);
      await new Promise((r) => setTimeout(r, 0));
    }
    gif.finish();
    const bytes = gif.bytes();
    downloadBlob(new Blob([bytes], { type: 'image/gif' }), 'pyslicer-print.gif');
    setExportStatus(`Saved GIF (${Math.round(bytes.length / 1024)} KB); truncated to ${maxSimSec.toFixed(1)}s @ ${speed}×`);
  } catch (err) {
    console.error(err);
    setExportStatus('GIF export failed (network/CDN or browser limits).');
  } finally {
    if (btn) btn.disabled = false;
    updateExportButtons();
  }
}

document.getElementById('export-movie-start')?.addEventListener('click', startMovieExport);
document.getElementById('export-movie-stop')?.addEventListener('click', stopMovieExport);
document.getElementById('export-gif')?.addEventListener('click', exportGifClip);
updateExportButtons();

window.PYSLICER_VIEWER = {
  seekSimulation,
  playSim,
  pauseSim,
  rewindSim,
  totalTime: () => totalTime,
  setSpeed: (s) => { sim.speed = s; },
  startMovieExport,
  stopMovieExport,
  exportGifClip,
};

function resize() {
  const w = mount.clientWidth;
  const h = Math.max(320, Math.round(w * 0.62));
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
resize();
window.addEventListener('resize', resize);

let lastTs = performance.now();
(function animate(now) {
  requestAnimationFrame(animate);
  const dt = Math.min(0.1, (now - lastTs) / 1000);
  lastTs = now;
  if (sim.playing && sim.mode === 'simulate') {
    const next = sim.t + dt * sim.speed;
    if (next >= totalTime) {
      seekSimulation(totalTime);
      pauseSim();
    } else {
      seekSimulation(next);
      if (simScrub) simScrub.value = String((sim.t / totalTime) * 100);
    }
  }
  controls.update();
  renderer.render(scene, camera);
})(lastTs);
"""
