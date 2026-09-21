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

function concernColor(t) {
  // Blue → yellow → red (Wong-ish)
  const u = Math.max(0, Math.min(1, Number(t) || 0));
  const c = new THREE.Color();
  if (u < 0.5) {
    const a = new THREE.Color(0x0072b2);
    const b = new THREE.Color(0xf0e442);
    c.copy(a).lerp(b, u * 2);
  } else {
    const a = new THREE.Color(0xf0e442);
    const b = new THREE.Color(0xd55e00);
    c.copy(a).lerp(b, (u - 0.5) * 2);
  }
  return c;
}

function scoresFromConcernMetrics(metrics, thr) {
  const feeds = metrics.feed || [];
  const n = feeds.length;
  const turn = metrics.turnRad || [];
  const accel = metrics.accelRaw || [];
  const jerk = metrics.jerkRaw || [];
  const minA = (Math.max(Number(thr.minAngleDeg) || 0, 0) * Math.PI) / 180;
  const maxSpeed = Math.max(Number(thr.maxSpeed) || 70, 1e-9);
  const maxAccel = Math.max(Number(thr.maxAccel) || 1000, 1e-9);
  const maxJerk = Math.max(Number(thr.maxJerk) || 20, 1e-9);
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    if ((turn[i] || 0) < minA) {
      out[i] = 0;
      continue;
    }
    // feeds are G-code F (mm/min); maxSpeed is mm/s
    const speedPart = ((feeds[i] || 0) / 60) / maxSpeed;
    const accelPart = (accel[i] || 0) / maxAccel;
    const jerkPart = (jerk[i] || 0) / maxJerk;
    // 0 while within thresholds; color only when limits are exceeded
    const excess = Math.max(speedPart, accelPart, jerkPart);
    out[i] = Math.min(1, Math.max(0, excess - 1));
  }
  return out;
}

function stadiumProfile2d(bead) {
  // Stadium / capsule: flat top/bottom, semicircle sides (x=sideways, y=up).
  const h = Math.max(Number(bead.height) || 0.4, 0.05);
  const w = Math.max(Number(bead.width) || h, h);
  const r = h / 2;
  const halfFlat = Math.max(w - h, 0) / 2;
  const arcSegs = 12;
  const pts = [];
  pts.push(new THREE.Vector2(-halfFlat, -r));
  pts.push(new THREE.Vector2(halfFlat, -r));
  for (let i = 1; i <= arcSegs; i++) {
    const a = -Math.PI / 2 + (Math.PI * i) / arcSegs;
    pts.push(new THREE.Vector2(halfFlat + Math.cos(a) * r, Math.sin(a) * r));
  }
  pts.push(new THREE.Vector2(-halfFlat, r));
  for (let i = 1; i <= arcSegs; i++) {
    const a = Math.PI / 2 + (Math.PI * i) / arcSegs;
    pts.push(new THREE.Vector2(-halfFlat + Math.cos(a) * r, Math.sin(a) * r));
  }
  return pts;
}

function buildSweepFrames(points, joinRadius) {
  // CAD round joins: at each corner, rotate the sweep frame at the vertex so
  // the outer bead edge is a circular arc of radius ≈ bead half-width centered
  // on the G-code tip. joinRadius sets how many arc subdivisions (via angle).
  const frames = [];
  if (points.length < 2) return frames;
  const push = (p, T) => {
    const t = T.clone();
    if (t.lengthSq() < 1e-14) t.set(1, 0, 0);
    else t.normalize();
    frames.push({ p: p.clone(), T: t });
  };

  const T0 = new THREE.Vector3().subVectors(points[1], points[0]);
  push(points[0], T0);

  for (let i = 1; i < points.length - 1; i++) {
    const A = points[i - 1];
    const B = points[i];
    const C = points[i + 1];
    const d1 = new THREE.Vector3().subVectors(B, A);
    const d2 = new THREE.Vector3().subVectors(C, B);
    const len1 = d1.length();
    const len2 = d2.length();
    if (len1 < 1e-9 || len2 < 1e-9) continue;
    d1.multiplyScalar(1 / len1);
    d2.multiplyScalar(1 / len2);
    const cos = Math.max(-1, Math.min(1, d1.dot(d2)));
    const phi = Math.acos(cos);
    if (phi < 0.02) {
      push(B, d2);
      continue;
    }
    const axis = new THREE.Vector3().crossVectors(d1, d2);
    if (axis.lengthSq() < 1e-14) {
      push(B, d2);
      continue;
    }
    axis.normalize();
    // Denser arc for larger turns / larger nozzles
    const nSegs = Math.max(
      4,
      Math.ceil((phi / (Math.PI / 16)) * Math.max(1, joinRadius / 0.2))
    );
    push(B, d1);
    for (let s = 1; s < nSegs; s++) {
      const q = new THREE.Quaternion().setFromAxisAngle(axis, (phi * s) / nSegs);
      push(B, d1.clone().applyQuaternion(q));
    }
    push(B, d2);
  }

  const Tend = new THREE.Vector3().subVectors(
    points[points.length - 1],
    points[points.length - 2]
  );
  push(points[points.length - 1], Tend);
  return frames;
}

function sweepStadiumGeometry(frames, profile, up) {
  // Solid stadium sweep: side wall + filled cross-sections (no hollow tube).
  const nRings = frames.length;
  const nProf = profile.length;
  if (nRings < 2 || nProf < 3) return null;

  // Layout: [ring0 profile | ring0 center | ring1 profile | ring1 center | ...]
  const stride = nProf + 1;
  const positions = new Float32Array(nRings * stride * 3);
  const indices = [];
  const B = new THREE.Vector3();
  const N = new THREE.Vector3();
  const altUp = new THREE.Vector3(1, 0, 0);

  for (let i = 0; i < nRings; i++) {
    const T = frames[i].T;
    B.crossVectors(T, up);
    if (B.lengthSq() < 1e-10) B.crossVectors(T, altUp);
    B.normalize();
    N.crossVectors(B, T).normalize();
    if (N.dot(up) < 0) {
      N.negate();
      B.negate();
    }
    const P = frames[i].p;
    const base = i * stride;
    for (let j = 0; j < nProf; j++) {
      const pr = profile[j];
      const o = (base + j) * 3;
      positions[o] = P.x + B.x * pr.x + N.x * pr.y;
      positions[o + 1] = P.y + B.y * pr.x + N.y * pr.y;
      positions[o + 2] = P.z + B.z * pr.x + N.z * pr.y;
    }
    const co = (base + nProf) * 3;
    positions[co] = P.x;
    positions[co + 1] = P.y;
    positions[co + 2] = P.z;
  }

  // Side walls between rings
  for (let i = 0; i < nRings - 1; i++) {
    for (let j = 0; j < nProf; j++) {
      const j2 = (j + 1) % nProf;
      const a = i * stride + j;
      const b = i * stride + j2;
      const c = (i + 1) * stride + j;
      const d = (i + 1) * stride + j2;
      indices.push(a, c, b, b, c, d);
    }
  }
  // Solid caps: fan from ring center to profile (both ends + every ring for join solidity)
  for (let i = 0; i < nRings; i++) {
    const center = i * stride + nProf;
    for (let j = 0; j < nProf; j++) {
      const j2 = (j + 1) % nProf;
      const a = i * stride + j;
      const b = i * stride + j2;
      // Alternate winding so both ends face outward-ish; DoubleSide covers the rest
      if (i === 0) indices.push(center, b, a);
      else indices.push(center, a, b);
    }
  }

  const indicesPerSeg = nProf * 6; // side wall only; caps handled separately for reveal
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geo.setIndex(indices);
  geo.computeVertexNormals();
  geo.userData.indicesPerSeg = indicesPerSeg;
  geo.userData.nSegs = nRings - 1;
  geo.userData.sideIndexCount = (nRings - 1) * indicesPerSeg;
  return geo;
}

function polylineFromSegments(segments) {
  // Fallback when payload lacks extrudePolylines.
  const polys = [];
  let pts = null;
  let i0 = 0;
  for (let i = 0; i < segments.length; i++) {
    const s = segments[i];
    const a = [s[0], s[1], s[2]];
    const b = [s[3], s[4], s[5]];
    if (!pts) {
      pts = [a, b];
      i0 = i;
    } else {
      const t = pts[pts.length - 1];
      if (
        Math.abs(t[0] - a[0]) < 1e-6 &&
        Math.abs(t[1] - a[1]) < 1e-6 &&
        Math.abs(t[2] - a[2]) < 1e-6
      ) {
        pts.push(b);
      } else {
        polys.push({ points: pts, i0, i1: i - 1 });
        pts = [a, b];
        i0 = i;
      }
    }
  }
  if (pts) polys.push({ points: pts, i0, i1: segments.length - 1 });
  return polys;
}

function addExtrudeSweeps(polylines, color, bead, nozzleDiameter) {
  if (!polylines.length) return null;
  const profile = stadiumProfile2d(bead || { height: 0.4, width: 0.5 });
  // Round-join / end-cap ball: match bead half-extent (not √2 miter tip).
  const nozzle = Math.max(Number(nozzleDiameter) || Number(bead && bead.width) || 0.5, 0.1);
  const filletR = nozzle * 0.5;
  const halfW = Math.max(filletR, (Number(bead && bead.width) || nozzle) * 0.5);
  const halfH = Math.max(Number(bead && bead.height) || nozzle, 0.1) * 0.5;
  const joinR = Math.min(halfW, halfH);
  const up = new THREE.Vector3(0, 1, 0);
  const pathColor = new THREE.Color(color);
  const mat = new THREE.MeshStandardMaterial({
    color: pathColor,
    roughness: 0.55,
    metalness: 0.02,
    flatShading: false,
    transparent: true,
    opacity: 1,
    depthWrite: true,
    side: THREE.DoubleSide,
    clippingPlanes: [clipPlane],
    vertexColors: false,
  });

  // Build per-polyline sweeps, then merge into one mesh (one draw call).
  const items = [];
  const positions = [];
  const indices = [];
  let vertBase = 0;
  let indexBase = 0;
  let totalSegs = 0;
  const joinSphere = new THREE.SphereGeometry(joinR, 12, 10);
  const joinPos = joinSphere.getAttribute('position');
  const joinIdx = joinSphere.index;
  const nProf = profile.length;
  const stride = nProf + 1;

  const sorted = polylines.slice().sort((a, b) => (a.i0 | 0) - (b.i0 | 0));
  for (const poly of sorted) {
    const raw = (poly.points || []).map((p) => toThree(p[0], p[1], p[2]));
    if (raw.length < 2) continue;
    const frames = buildSweepFrames(raw, filletR);
    const geo = sweepStadiumGeometry(frames, profile, up);
    if (!geo) continue;

    const vertStart = vertBase;
    const pos = geo.getAttribute('position');
    for (let i = 0; i < pos.count; i++) {
      positions.push(pos.getX(i), pos.getY(i), pos.getZ(i));
    }
    const idx = geo.index;
    let indexCount = idx.count;
    for (let i = 0; i < indexCount; i++) {
      indices.push(idx.getX(i) + vertBase);
    }
    vertBase += pos.count;

    const nRings = frames.length;
    const joins = [];
    // Sphere joins only at turning corners (not open tips — those looked bulbous).
    for (let vi = 1; vi < raw.length - 1; vi++) {
      const d1 = new THREE.Vector3().subVectors(raw[vi], raw[vi - 1]);
      const d2 = new THREE.Vector3().subVectors(raw[vi + 1], raw[vi]);
      if (d1.lengthSq() < 1e-12 || d2.lengthSq() < 1e-12) continue;
      d1.normalize();
      d2.normalize();
      if (d1.dot(d2) > 0.998) continue; // ~colinear
      const c = raw[vi];
      const v0 = vertBase;
      for (let i = 0; i < joinPos.count; i++) {
        positions.push(
          joinPos.getX(i) + c.x,
          joinPos.getY(i) + c.y,
          joinPos.getZ(i) + c.z
        );
      }
      for (let i = 0; i < joinIdx.count; i++) {
        indices.push(joinIdx.getX(i) + v0);
      }
      const i0 = poly.i0 | 0;
      joins.push({
        vertStart: v0,
        vertCount: joinPos.count,
        segA: i0 + vi - 1,
        segB: i0 + vi,
      });
      vertBase += joinPos.count;
      indexCount += joinIdx.count;
    }

    const i0 = poly.i0 | 0;
    const i1 = poly.i1 | 0;
    const nEdges = Math.max(1, i1 - i0 + 1);
    items.push({
      i0,
      i1,
      nEdges,
      vertStart,
      nRings,
      stride,
      nProf,
      joins,
      spineSegs: geo.userData.nSegs,
      indicesPerSeg: geo.userData.indicesPerSeg,
      indexStart: indexBase,
      indexCount,
    });
    indexBase += indexCount;
    totalSegs += nEdges;
    geo.dispose();
  }
  joinSphere.dispose();

  if (!items.length) return null;

  const merged = new THREE.BufferGeometry();
  merged.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  const vertCount = positions.length / 3;
  const colors = new Float32Array(vertCount * 3);
  for (let i = 0; i < vertCount; i++) {
    colors[i * 3] = pathColor.r;
    colors[i * 3 + 1] = pathColor.g;
    colors[i * 3 + 2] = pathColor.b;
  }
  merged.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
  merged.setIndex(indices);
  merged.computeVertexNormals();
  const mesh = new THREE.Mesh(merged, mat);
  scene.add(mesh);
  return {
    mesh,
    geometry: merged,
    items,
    material: mat,
    total: totalSegs,
    pathColor,
  };
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
const extrudePolylines =
  data.extrudePolylines && data.extrudePolylines.length
    ? data.extrudePolylines
    : polylineFromSegments(data.extrude || []);
const extrudeObj = addExtrudeSweeps(
  extrudePolylines,
  COLOR_EXTRUDE,
  data.bead || {
    height: data.layerHeight || 0.4,
    width: data.tubeDiameter || data.nozzleDiameter || 0.5,
  },
  data.nozzleDiameter
);
const concernMetrics = data.concernMetrics || {
  feed: [],
  turnRad: [],
  accelRaw: [],
  jerkRaw: [],
};
const planningLimits = data.planningLimits || {};
const concernThresholds = {
  maxSpeed: Number(planningLimits.maxSpeed) || 70,
  maxAccel: Number(planningLimits.maxAccel) || 1000,
  maxJerk: Number(planningLimits.maxJerk) || 20,
  minAngleDeg: Number(planningLimits.minAngleDeg) || 20,
};
let colorMode = 'path';

function fmtNum(v, digits) {
  const n = Number(v);
  if (!Number.isFinite(n)) return '—';
  return digits == null ? String(n) : n.toFixed(digits);
}

function populateSliceSettings() {
  const list = document.getElementById('slice-settings');
  const missing = document.getElementById('slice-settings-missing');
  if (!list) return;
  const pl = data.planningLimits;
  if (!pl || Object.keys(pl).length === 0) {
    list.hidden = true;
    if (missing) missing.hidden = false;
    return;
  }
  if (missing) missing.hidden = true;
  list.hidden = false;
  const rows = [
    ['Outer wall', `${fmtNum(pl.outerSpeed ?? pl.maxSpeed, 1)} mm/s`],
    ['Inner wall', `${fmtNum(pl.innerSpeed ?? pl.maxSpeed, 1)} mm/s`],
    ['Infill', `${fmtNum(pl.infillSpeed ?? pl.maxSpeed, 1)} mm/s`],
    ['Max accel', `${fmtNum(pl.maxAccel, 0)} mm/s²`],
    ['Max corner Δv', `${fmtNum(pl.maxJerk, 1)} mm/s`],
    ['90° corner speed', `${fmtNum(pl.maxCornerSpeed, 1)} mm/s`],
    ['Min turn angle', `${fmtNum(pl.minAngleDeg, 0)}°`],
  ];
  list.innerHTML = rows
    .map(([k, v]) => `<li><span class="k">${k}</span><span class="v">${v}</span></li>`)
    .join('');
}
populateSliceSettings();
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
    let indexCount = 0;
    for (const item of extrudeObj.items) {
      const edges = Math.max(0, Math.min(item.nEdges, extrudeCount - item.i0));
      if (edges <= 0) break;
      if (edges >= item.nEdges) {
        indexCount = item.indexStart + item.indexCount;
      } else {
        const spineSegs = Math.max(
          1,
          Math.round((edges / item.nEdges) * item.spineSegs)
        );
        indexCount = item.indexStart + spineSegs * item.indicesPerSeg;
      }
    }
    extrudeObj.geometry.setDrawRange(0, indexCount);
  }
  if (travelObj) {
    const n = Math.max(0, Math.min(travelObj.total, travelCount));
    travelObj.geometry.setDrawRange(0, n * 2);
  }
}

function showAllFilament() {
  if (extrudeObj) {
    let n = 0;
    for (const item of extrudeObj.items) n += item.indexCount;
    extrudeObj.geometry.setDrawRange(0, n);
  }
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

function paintVertexRange(colors, vertStart, vertCount, c) {
  for (let i = 0; i < vertCount; i++) {
    const o = (vertStart + i) * 3;
    colors[o] = c.r;
    colors[o + 1] = c.g;
    colors[o + 2] = c.b;
  }
}

function paintExtrudeConcernColors(scores) {
  if (!extrudeObj) return;
  const colors = extrudeObj.geometry.getAttribute('color');
  if (!colors) return;
  const arr = colors.array;
  for (const item of extrudeObj.items) {
    const nRings = item.nRings | 0;
    const stride = item.stride | 0;
    for (let r = 0; r < nRings; r++) {
      const t = nRings <= 1 ? 0 : r / (nRings - 1);
      const edge =
        item.i0 + Math.min(item.nEdges - 1, Math.floor(t * item.nEdges + 1e-9));
      const c = concernColor(scores[edge] || 0);
      paintVertexRange(arr, item.vertStart + r * stride, stride, c);
    }
    for (const j of item.joins || []) {
      const s = Math.max(scores[j.segA] || 0, scores[j.segB] || 0);
      paintVertexRange(arr, j.vertStart, j.vertCount, concernColor(s));
    }
  }
  colors.needsUpdate = true;
}

function readConcernThresholds() {
  const speedEl = document.getElementById('thr-speed');
  const accelEl = document.getElementById('thr-accel');
  const jerkEl = document.getElementById('thr-jerk');
  const angleEl = document.getElementById('thr-angle');
  if (speedEl) concernThresholds.maxSpeed = Number(speedEl.value) || 70;
  if (accelEl) concernThresholds.maxAccel = Number(accelEl.value) || 1000;
  if (jerkEl) concernThresholds.maxJerk = Number(jerkEl.value) || 20;
  if (angleEl) concernThresholds.minAngleDeg = Number(angleEl.value) || 0;
  const speedVal = document.getElementById('thr-speed-val');
  const accelVal = document.getElementById('thr-accel-val');
  const jerkVal = document.getElementById('thr-jerk-val');
  const angleVal = document.getElementById('thr-angle-val');
  if (speedVal) speedVal.textContent = String(Number(concernThresholds.maxSpeed));
  if (accelVal) accelVal.textContent = String(Math.round(concernThresholds.maxAccel));
  if (jerkVal) jerkVal.textContent = String(Number(concernThresholds.maxJerk));
  if (angleVal) angleVal.textContent = String(Math.round(concernThresholds.minAngleDeg));
}

function applyPlanningLimitsToSliders() {
  const speedEl = document.getElementById('thr-speed');
  const accelEl = document.getElementById('thr-accel');
  const jerkEl = document.getElementById('thr-jerk');
  const angleEl = document.getElementById('thr-angle');
  if (speedEl) speedEl.value = String(concernThresholds.maxSpeed);
  if (accelEl) accelEl.value = String(concernThresholds.maxAccel);
  if (jerkEl) jerkEl.value = String(concernThresholds.maxJerk);
  if (angleEl) angleEl.value = String(concernThresholds.minAngleDeg);
  readConcernThresholds();
}

function refreshConcernColors() {
  readConcernThresholds();
  const scores = scoresFromConcernMetrics(concernMetrics, concernThresholds);
  paintExtrudeConcernColors(scores);
}

function setColorMode(mode) {
  colorMode = mode === 'concern' ? 'concern' : 'path';
  const concern = colorMode === 'concern';
  const thrPanel = document.getElementById('concern-thresholds');
  if (thrPanel) thrPanel.hidden = !concern;
  if (travelObj && travelObj.lines) {
    travelObj.lines.visible = !concern;
  }
  if (extrudeObj && extrudeObj.mesh && extrudeObj.material) {
    extrudeObj.mesh.visible = true;
    const opacity = opacityInput
      ? Math.max(0.05, Number(opacityInput.value) / 100)
      : 1;
    extrudeObj.material.opacity = opacity;
    extrudeObj.material.transparent = true;
    extrudeObj.material.depthWrite = opacity > 0.95;
    if (concern) {
      extrudeObj.material.vertexColors = true;
      extrudeObj.material.color.set(0xffffff);
      refreshConcernColors();
    } else {
      extrudeObj.material.vertexColors = false;
      extrudeObj.material.color.copy(extrudeObj.pathColor);
    }
    extrudeObj.material.needsUpdate = true;
  }
  const legend = document.getElementById('concern-legend');
  if (legend) {
    legend.classList.toggle('visible', concern);
    legend.setAttribute('aria-hidden', concern ? 'false' : 'true');
  }
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

  updateSimLabels(pathPct, curZ, pos);
}

function moveIndexAtTime(t) {
  if (!timeline.length) return -1;
  for (let i = 0; i < timeline.length; i++) {
    if (t <= timeline[i].t1 + 1e-9) return i;
  }
  return timeline.length - 1;
}

function updateSimLabels(pathPct, curZ, pos) {
  const pctEl = document.getElementById('sim-pct');
  const zEl = document.getElementById('sim-z');
  const moveEl = document.getElementById('sim-move');
  if (pctEl) {
    const p = pathPct != null ? pathPct : (sim.t / totalTime) * 100;
    pctEl.textContent = `${Math.round(p)}%`;
  }
  if (zEl && curZ != null) zEl.textContent = `Z ${curZ.toFixed(2)}`;
  if (moveEl && timeline.length) {
    const i = moveIndexAtTime(sim.t);
    const m = timeline[i];
    const kind = m.extrude ? 'extrude' : 'travel';
    const xy = pos
      ? `X${pos.x.toFixed(2)} Y${(-pos.z).toFixed(2)}`
      : `X${m.x1.toFixed(2)} Y${m.y1.toFixed(2)}`;
    moveEl.textContent = `Move ${i + 1}/${timeline.length} · ${kind} · ${xy}`;
  } else if (moveEl) {
    moveEl.textContent = 'Move —';
  }
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
  pauseSim();
  seekSimulation(Math.min(totalTime, sim.t + Math.max(totalTime * 0.05, 2)));
}

function stepForwardSim() {
  if (!timeline.length) return;
  pauseSim();
  sim.mode = 'simulate';
  const i = moveIndexAtTime(sim.t);
  const m = timeline[i];
  // Finish current move, or advance to the end of the next move
  if (sim.t < m.t1 - 1e-9) {
    seekSimulation(m.t1);
  } else if (i + 1 < timeline.length) {
    seekSimulation(timeline[i + 1].t1);
  }
}

function stepBackwardSim() {
  if (!timeline.length) return;
  pauseSim();
  sim.mode = 'simulate';
  const i = moveIndexAtTime(sim.t);
  const m = timeline[i];
  // Snap to start of this move, or previous move start
  if (sim.t > m.t0 + 1e-9) {
    seekSimulation(m.t0);
  } else if (i > 0) {
    seekSimulation(timeline[i - 1].t0);
  } else {
    seekSimulation(0);
  }
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

document.querySelectorAll('input[name="color-mode"]').forEach((el) => {
  el.addEventListener('change', () => {
    if (el.checked) setColorMode(el.value);
  });
});
['thr-speed', 'thr-accel', 'thr-jerk', 'thr-angle'].forEach((id) => {
  const el = document.getElementById(id);
  if (!el) return;
  el.addEventListener('input', () => {
    if (colorMode === 'concern') refreshConcernColors();
    else readConcernThresholds();
  });
});
applyPlanningLimitsToSliders();
const checkedColor = document.querySelector('input[name="color-mode"]:checked');
setColorMode(checkedColor ? checkedColor.value : 'path');

document.getElementById('sim-play')?.addEventListener('click', playSim);
document.getElementById('sim-pause')?.addEventListener('click', pauseSim);
document.getElementById('sim-rewind')?.addEventListener('click', rewindSim);
document.getElementById('sim-ff')?.addEventListener('click', fastForwardSim);
document.getElementById('sim-step-fwd')?.addEventListener('click', stepForwardSim);
document.getElementById('sim-step-back')?.addEventListener('click', stepBackwardSim);

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
  stepForwardSim,
  stepBackwardSim,
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
