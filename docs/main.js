import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

// ---------------------------------------------------------------------------
// Physics: U(x,y;lambda,eps) = 1 + lambda/r1 + (1+eps)/r2 + (1-eps)/r3, N = 1/U.
// Source layout and lambda-at-P1 convention match src/mp_model.jl exactly;
// eps=0 reduces to the symmetric (1,1,lambda) family.
// ---------------------------------------------------------------------------

const SOURCES = [
  { x: 1, y: 0, label: "λ" },
  { x: -0.5, y: Math.sqrt(3) / 2, label: "1+ε" },
  { x: -0.5, y: -Math.sqrt(3) / 2, label: "1−ε" },
];

const LAMBDA_MINUS = 0.6738774744707380113762233701706991;
const LAMBDA_PLUS = 1.1362522106646809015183512524900;

function masses(lambda, eps) {
  return [lambda, 1 + eps, 1 - eps];
}

function U(x, y, lambda, eps) {
  const m = masses(lambda, eps);
  let u = 1;
  for (let i = 0; i < 3; i++) {
    const dx = x - SOURCES[i].x, dy = y - SOURCES[i].y;
    const r = Math.max(Math.hypot(dx, dy), 1e-6);
    u += m[i] / r;
  }
  return u;
}

function gradU(x, y, lambda, eps) {
  const m = masses(lambda, eps);
  let gx = 0, gy = 0;
  for (let i = 0; i < 3; i++) {
    const dx = x - SOURCES[i].x, dy = y - SOURCES[i].y;
    const r2 = Math.max(dx * dx + dy * dy, 1e-12);
    const r3 = r2 * Math.sqrt(r2);
    gx += -m[i] * dx / r3;
    gy += -m[i] * dy / r3;
  }
  return [gx, gy];
}

function hessU(x, y, lambda, eps) {
  const m = masses(lambda, eps);
  let uxx = 0, uxy = 0, uyy = 0;
  for (let i = 0; i < 3; i++) {
    const dx = x - SOURCES[i].x, dy = y - SOURCES[i].y;
    const r2 = Math.max(dx * dx + dy * dy, 1e-12);
    const r = Math.sqrt(r2), r3 = r2 * r, r5 = r3 * r2;
    uxx += m[i] * (3 * dx * dx / r5 - 1 / r3);
    uyy += m[i] * (3 * dy * dy / r5 - 1 / r3);
    uxy += m[i] * (3 * dx * dy / r5);
  }
  return [uxx, uxy, uyy];
}

function N(x, y, lambda, eps) {
  return 1 / U(x, y, lambda, eps);
}

// 2x2 symmetric eigen-decomposition of [[a,b],[b,d]], closed form.
// Returns {mu1 (larger), v1, mu2 (smaller), v2}, unit eigenvectors.
function eigSym2(a, b, d) {
  const tr = a + d, diff = (a - d) / 2;
  const disc = Math.sqrt(diff * diff + b * b);
  const mu1 = tr / 2 + disc;
  const mu2 = tr / 2 - disc;
  function vecFor(mu) {
    let vx, vy;
    if (Math.abs(b) > 1e-10) { vx = b; vy = mu - a; }
    else if (a >= d) { vx = 1; vy = 0; }
    else { vx = 0; vy = 1; }
    const nrm = Math.hypot(vx, vy) || 1;
    return [vx / nrm, vy / nrm];
  }
  return { mu1, v1: vecFor(mu1), mu2, v2: vecFor(mu2) };
}

// ---------------------------------------------------------------------------
// Critical-point finder: coarse grid seeding + 2D Newton refinement.
// This is a plain, uncertified solver for a smooth real-time picture; the
// repository's proof/ scripts certify the same structure rigorously.
// ---------------------------------------------------------------------------

const DOMAIN = { xmin: -1.35, xmax: 1.65, ymin: -1.45, ymax: 1.45 };
const EXCLUDE_R = 0.12;

// Newton refinement of a single seed toward a critical point of U(.,.,lambda,eps).
// Returns {x,y} on success, null if it diverges or the Hessian goes singular.
function newtonRefine(x0, y0, lambda, eps) {
  let x = x0, y = y0;
  for (let iter = 0; iter < 25; iter++) {
    const [gx, gy] = gradU(x, y, lambda, eps);
    const [uxx, uxy, uyy] = hessU(x, y, lambda, eps);
    const det = uxx * uyy - uxy * uxy;
    if (Math.abs(det) < 1e-10) return null;
    const dx = (uyy * gx - uxy * gy) / det;
    const dy = (-uxy * gx + uxx * gy) / det;
    x -= dx; y -= dy;
    if (!isFinite(x) || !isFinite(y) || Math.abs(x) > 3 || Math.abs(y) > 3) return null;
    if (Math.hypot(dx, dy) < 1e-11) return { x, y };
  }
  return null;
}

// Previous frame's converged points, reused as extra Newton seeds so a pair
// that is close together (near a bifurcation) keeps resolving correctly even
// though the coarse grid below can no longer distinguish them as two minima.
let previousPoints = [];

function findCriticalPoints(lambda, eps) {
  const NX = 90, NY = 90;
  const grid = new Float64Array((NX + 1) * (NY + 1));
  const idx = (i, j) => i * (NY + 1) + j;

  for (let i = 0; i <= NX; i++) {
    const x = DOMAIN.xmin + ((DOMAIN.xmax - DOMAIN.xmin) * i) / NX;
    for (let j = 0; j <= NY; j++) {
      const y = DOMAIN.ymin + ((DOMAIN.ymax - DOMAIN.ymin) * j) / NY;
      let excluded = false;
      for (const s of SOURCES) {
        if ((x - s.x) ** 2 + (y - s.y) ** 2 < EXCLUDE_R * EXCLUDE_R) { excluded = true; break; }
      }
      if (excluded) { grid[idx(i, j)] = Infinity; continue; }
      const [gx, gy] = gradU(x, y, lambda, eps);
      grid[idx(i, j)] = gx * gx + gy * gy;
    }
  }

  const seeds = [];
  for (let i = 1; i < NX; i++) {
    for (let j = 1; j < NY; j++) {
      const v = grid[idx(i, j)];
      if (!isFinite(v)) continue;
      let isMin = true;
      for (let di = -1; di <= 1 && isMin; di++) {
        for (let dj = -1; dj <= 1; dj++) {
          if (di === 0 && dj === 0) continue;
          if (grid[idx(i + di, j + dj)] < v) { isMin = false; break; }
        }
      }
      if (isMin) {
        seeds.push([
          DOMAIN.xmin + ((DOMAIN.xmax - DOMAIN.xmin) * i) / NX,
          DOMAIN.ymin + ((DOMAIN.ymax - DOMAIN.ymin) * j) / NY,
        ]);
      }
    }
  }
  // Extra seeds around each previously-tracked point: not just the point
  // itself, but a small ring around it. A fold/pitchfork births two critical
  // points arbitrarily close together, and while they're still closer than
  // one grid cell apart, the grid scan above sees only a single local
  // minimum of |grad U|^2 and can only ever find one of them; a lone seed at
  // the old point converges back to the same single point every time. The
  // ring gives Newton a chance to fall into the *other* one's basin as soon
  // as they're separated by more than ~PERTURB_R, well before the grid can
  // resolve them on its own.
  const PERTURB_R = 0.03;
  for (const p of previousPoints) {
    seeds.push([p.x, p.y]);
    for (let k = 0; k < 8; k++) {
      const theta = (k / 8) * 2 * Math.PI;
      seeds.push([p.x + PERTURB_R * Math.cos(theta), p.y + PERTURB_R * Math.sin(theta)]);
    }
  }

  const results = [];
  for (const [x0, y0] of seeds) {
    const hit = newtonRefine(x0, y0, lambda, eps);
    if (!hit) continue;
    const { x, y } = hit;
    if (SOURCES.some((s) => Math.hypot(x - s.x, y - s.y) < 0.05)) continue;
    if (results.some((r) => Math.hypot(r.x - x, r.y - y) < 0.01)) continue;

    const [uxx, uxy, uyy] = hessU(x, y, lambda, eps);
    const det = uxx * uyy - uxy * uxy;
    const trace = uxx + uyy;
    const type = det < 0 ? "S1" : trace > 0 ? "S2" : "S0";
    results.push({ x, y, type });
  }
  previousPoints = results;
  return results;
}

// ---------------------------------------------------------------------------
// Gradient-flow integration: separatrices (stable/unstable manifolds) and
// basin-of-attraction classification. Both follow the SAME flow (steepest
// ascent of U == steepest descent of N, since N=1/U is monotonic in U, so
// dN/dt<0 wherever dU/dt>0); "unstable" / "stable" below always refers to
// the descending-lapse-flow convention the paper's Figure 1 uses.
// ---------------------------------------------------------------------------

function integrateFlow(x0, y0, lambda, eps, dir, otherPoints, opts = {}) {
  const step = opts.step ?? 0.014;
  const maxSteps = opts.maxSteps ?? 240;
  const pts = [[x0, y0]];
  let x = x0, y = y0;
  let end = { type: "maxsteps" };
  for (let i = 0; i < maxSteps; i++) {
    let hitSource = -1;
    for (let s = 0; s < SOURCES.length; s++) {
      if (Math.hypot(x - SOURCES[s].x, y - SOURCES[s].y) < 0.1) { hitSource = s; break; }
    }
    if (hitSource >= 0) { end = { type: "source", index: hitSource }; break; }
    if (otherPoints) {
      const hit = otherPoints.find(
        (p) => Math.hypot(p.x - x, p.y - y) < 0.06 && Math.hypot(p.x - x0, p.y - y0) > 0.08
      );
      if (hit) { end = { type: "critpoint", pointType: hit.type }; break; }
    }
    if (x < DOMAIN.xmin + 0.03 || x > DOMAIN.xmax - 0.03 || y < DOMAIN.ymin + 0.03 || y > DOMAIN.ymax - 0.03) {
      end = { type: "domain" };
      break;
    }
    const [gx, gy] = gradU(x, y, lambda, eps);
    const gn = Math.hypot(gx, gy) || 1;
    x += dir * step * (gx / gn);
    y += dir * step * (gy / gn);
    pts.push([x, y]);
  }
  return { pts, end };
}

// From an S1 saddle: 2 unstable branches (solid, run to a source — these ARE
// the basin-boundary separatrices) + 2 stable branches (dashed, traced via
// the reversed flow from the saddle, typically running up to an S2 hilltop
// or off toward the domain boundary).
function computeManifoldBranches(saddleX, saddleY, lambda, eps, otherPoints) {
  const [uxx, uxy, uyy] = hessU(saddleX, saddleY, lambda, eps);
  const { mu1, v1, mu2, v2 } = eigSym2(uxx, uxy, uyy);
  const posV = mu1 > 0 ? v1 : v2; // U-increasing direction -> N-decreasing -> unstable for descending-N flow
  const negV = mu1 > 0 ? v2 : v1; // U-decreasing direction -> N-increasing -> stable for descending-N flow
  const seedEps = 0.02;

  const unstable = [+1, -1].map((sgn) =>
    integrateFlow(
      saddleX + sgn * seedEps * posV[0],
      saddleY + sgn * seedEps * posV[1],
      lambda, eps, +1, otherPoints
    )
  );
  const stable = [+1, -1].map((sgn) =>
    integrateFlow(
      saddleX + sgn * seedEps * negV[0],
      saddleY + sgn * seedEps * negV[1],
      lambda, eps, -1, otherPoints
    )
  );
  return { unstable, stable };
}

// Coarse basin-of-attraction grid: which source each sample point's
// ascending-U (= descending-N) flow ends up closest to. Cheap, approximate,
// intended only for a color overlay, not a certified partition.
const BASIN_NX = 42, BASIN_NY = 42;
function computeBasinGrid(lambda, eps) {
  const labels = new Int8Array((BASIN_NX + 1) * (BASIN_NY + 1));
  let k = 0;
  for (let i = 0; i <= BASIN_NX; i++) {
    const x0 = DOMAIN.xmin + ((DOMAIN.xmax - DOMAIN.xmin) * i) / BASIN_NX;
    for (let j = 0; j <= BASIN_NY; j++, k++) {
      const y0 = DOMAIN.ymin + ((DOMAIN.ymax - DOMAIN.ymin) * j) / BASIN_NY;
      let x = x0, y = y0, label = -1;
      for (let s = 0; s < 34; s++) {
        let hitSource = -1;
        for (let si = 0; si < SOURCES.length; si++) {
          if (Math.hypot(x - SOURCES[si].x, y - SOURCES[si].y) < 0.12) { hitSource = si; break; }
        }
        if (hitSource >= 0) { label = hitSource; break; }
        const [gx, gy] = gradU(x, y, lambda, eps);
        const gn = Math.hypot(gx, gy) || 1;
        x += 0.05 * (gx / gn);
        y += 0.05 * (gy / gn);
      }
      if (label < 0) {
        // didn't converge in the step budget: fall back to nearest source
        let bestD = Infinity;
        for (let si = 0; si < SOURCES.length; si++) {
          const d = Math.hypot(x - SOURCES[si].x, y - SOURCES[si].y);
          if (d < bestD) { bestD = d; label = si; }
        }
      }
      labels[k] = label;
    }
  }
  return labels;
}
function basinLabelAt(labels, x, y) {
  let i = Math.round(((x - DOMAIN.xmin) / (DOMAIN.xmax - DOMAIN.xmin)) * BASIN_NX);
  let j = Math.round(((y - DOMAIN.ymin) / (DOMAIN.ymax - DOMAIN.ymin)) * BASIN_NY);
  i = Math.min(Math.max(i, 0), BASIN_NX);
  j = Math.min(Math.max(j, 0), BASIN_NY);
  return labels[i * (BASIN_NY + 1) + j];
}

// ---------------------------------------------------------------------------
// Three.js scene
// ---------------------------------------------------------------------------

const HEIGHT_SCALE = 3.4;
const SEG = 100;

const canvas = document.getElementById("scene");
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));

const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera(42, 1, 0.1, 100);
const DEFAULT_CAM_POS = new THREE.Vector3(3.4, 3.6, 0.0);
const DEFAULT_CAM_TARGET = new THREE.Vector3(0.0, 0.6, 0.0);
camera.position.copy(DEFAULT_CAM_POS);

const controls = new OrbitControls(camera, renderer.domElement);
controls.target.copy(DEFAULT_CAM_TARGET);
controls.enableDamping = true;
controls.dampingFactor = 0.08;
controls.minDistance = 1.5;
controls.maxDistance = 12;

scene.add(new THREE.AmbientLight(0xffffff, 0.55));
const key = new THREE.DirectionalLight(0xffffff, 1.1);
key.position.set(3, 6, 4);
scene.add(key);
const rim = new THREE.DirectionalLight(0x6ea8fe, 0.35);
rim.position.set(-4, 2, -3);
scene.add(rim);

// --- surface geometry (built once, positions/colors updated per frame) ---

const geometry = new THREE.PlaneGeometry(
  DOMAIN.xmax - DOMAIN.xmin,
  DOMAIN.ymax - DOMAIN.ymin,
  SEG,
  SEG
);
geometry.rotateX(-Math.PI / 2);
const posAttr = geometry.attributes.position;
const colorArr = new Float32Array(posAttr.count * 3);
geometry.setAttribute("color", new THREE.BufferAttribute(colorArr, 3));

const cxOffset = (DOMAIN.xmin + DOMAIN.xmax) / 2;
const cyOffset = (DOMAIN.ymin + DOMAIN.ymax) / 2;

const lowColor = new THREE.Color(0x1b3a6b);
const midColor = new THREE.Color(0x3f8fd6);
const highColor = new THREE.Color(0xf4e6c1);

function heightColor(n, out) {
  const t = Math.min(Math.max(n, 0), 1);
  if (t < 0.5) out.copy(lowColor).lerp(midColor, t / 0.5);
  else out.copy(midColor).lerp(highColor, (t - 0.5) / 0.5);
  return out;
}

const BASIN_COLORS = [new THREE.Color(0xdc6b3a), new THREE.Color(0x3aa3dc), new THREE.Color(0x6bc46a)];
const tmpBasinColor = new THREE.Color();

const material = new THREE.MeshStandardMaterial({
  vertexColors: true,
  roughness: 0.65,
  metalness: 0.05,
  side: THREE.DoubleSide,
  flatShading: false,
});
const surface = new THREE.Mesh(geometry, material);
scene.add(surface);

const tmpColor = new THREE.Color();

function updateSurface(lambda, eps, colorMode, basinLabels) {
  for (let i = 0; i < posAttr.count; i++) {
    const x = posAttr.getX(i) + cxOffset;
    const z = posAttr.getZ(i) + cyOffset;
    const n = N(x, z, lambda, eps);
    posAttr.setY(i, n * HEIGHT_SCALE);

    if (colorMode === "basin" && basinLabels) {
      const lbl = basinLabelAt(basinLabels, x, z);
      const base = BASIN_COLORS[lbl] ?? lowColor;
      const shade = 0.55 + 0.45 * Math.min(Math.max(n, 0), 1);
      tmpBasinColor.copy(base).multiplyScalar(shade);
      colorArr[i * 3] = tmpBasinColor.r;
      colorArr[i * 3 + 1] = tmpBasinColor.g;
      colorArr[i * 3 + 2] = tmpBasinColor.b;
    } else {
      heightColor(n, tmpColor);
      colorArr[i * 3] = tmpColor.r;
      colorArr[i * 3 + 1] = tmpColor.g;
      colorArr[i * 3 + 2] = tmpColor.b;
    }
  }
  posAttr.needsUpdate = true;
  geometry.attributes.color.needsUpdate = true;
  geometry.computeVertexNormals();
}

// --- source markers + labels ---

function makeLabelSprite(text) {
  const size = 160;
  const c = document.createElement("canvas");
  c.width = c.height = size;
  const ctx = c.getContext("2d");
  ctx.font = "bold 46px sans-serif";
  ctx.fillStyle = "#0b0e14";
  ctx.textAlign = "center";
  ctx.textBaseline = "middle";
  ctx.fillText(text, size / 2, size / 2 + 4);
  const tex = new THREE.CanvasTexture(c);
  const mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthTest: false, depthWrite: false });
  const sprite = new THREE.Sprite(mat);
  sprite.scale.set(0.26, 0.26, 1);
  sprite.renderOrder = 999;
  return sprite;
}

const sourceGroup = new THREE.Group();
scene.add(sourceGroup);
const sourceLabels = [];
for (const s of SOURCES) {
  const g = new THREE.Group();
  const sphere = new THREE.Mesh(
    new THREE.SphereGeometry(0.055, 24, 24),
    new THREE.MeshStandardMaterial({ color: 0xffffff, emissive: 0x333333, roughness: 0.3 })
  );
  g.add(sphere);
  const label = makeLabelSprite(s.label);
  label.position.set(0, 0.24, 0);
  g.add(label);
  g.position.set(s.x, 0.02, s.y);
  sourceGroup.add(g);
  sourceLabels.push(label);
}

function updateSourceLabels(eps) {
  const texts = [SOURCES[0].label, Math.abs(eps) < 1e-6 ? "1" : "1+ε", Math.abs(eps) < 1e-6 ? "1" : "1−ε"];
  sourceLabels.forEach((sprite, i) => {
    if (sprite.userData.lastText === texts[i]) return;
    sprite.userData.lastText = texts[i];
    const c = document.createElement("canvas");
    c.width = c.height = 160;
    const ctx = c.getContext("2d");
    ctx.font = "bold 44px sans-serif";
    ctx.fillStyle = "#0b0e14";
    ctx.textAlign = "center";
    ctx.textBaseline = "middle";
    ctx.fillText(texts[i], 80, 84);
    sprite.material.map.dispose();
    sprite.material.map = new THREE.CanvasTexture(c);
    sprite.material.needsUpdate = true;
  });
}

// --- critical point markers + manifold branches (persistent identity) ---
//
// Points are re-detected from scratch every frame, so naively assigning
// "point i -> marker i" makes a marker jump to a different physical point
// whenever the detected order changes — most visibly right at a
// bifurcation, where points are born, merge, or swap position in the array.
// Instead each marker keeps track of the point it last represented and,
// every frame, claims whichever surviving point is nearest to where it was;
// only a marker with no nearby point left goes idle. This keeps a
// continuously-existing critical point locked to the same mesh instance
// (smooth motion, no flicker) and lets its manifold Line objects live on
// the same persistent slot, so the separatrices move smoothly too.

const critGroup = new THREE.Group();
scene.add(critGroup);
const critColors = { S1: 0xf0a04b, S2: 0x5bc8f5, S0: 0x7be08e };
const MAX_CRIT_MARKERS = 8;
const MAX_MATCH_JUMP = 0.4; // world units; beyond this, treat as a different point

function makeManifoldLine(dashed) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute("position", new THREE.BufferAttribute(new Float32Array(3 * 3), 3));
  const mat = dashed
    ? new THREE.LineDashedMaterial({ color: 0xf0a04b, dashSize: 0.045, gapSize: 0.03, transparent: true, opacity: 0.85 })
    : new THREE.LineBasicMaterial({ color: 0xf0a04b, transparent: true, opacity: 0.95 });
  const line = new THREE.Line(geo, mat);
  line.visible = false;
  line.frustumCulled = false;
  critGroup.add(line);
  return line;
}

const markerState = [];
for (let i = 0; i < MAX_CRIT_MARKERS; i++) {
  const mesh = new THREE.Mesh(
    new THREE.OctahedronGeometry(0.055, 0),
    new THREE.MeshStandardMaterial({ color: 0xffffff, emissive: 0x222222, roughness: 0.35 })
  );
  mesh.visible = false;
  critGroup.add(mesh);
  markerState.push({
    mesh,
    x: 0,
    y: 0,
    type: null,
    active: false,
    branches: [makeManifoldLine(false), makeManifoldLine(false), makeManifoldLine(true), makeManifoldLine(true)],
  });
}

function setLineFromPath(line, pathPts, lambda, eps, color) {
  const n = pathPts.length;
  const arr = new Float32Array(n * 3);
  for (let i = 0; i < n; i++) {
    const [x, y] = pathPts[i];
    const h = N(x, y, lambda, eps) * HEIGHT_SCALE;
    arr[i * 3] = x;
    arr[i * 3 + 1] = h + 0.015;
    arr[i * 3 + 2] = y;
  }
  line.geometry.dispose();
  line.geometry = new THREE.BufferGeometry();
  line.geometry.setAttribute("position", new THREE.BufferAttribute(arr, 3));
  if (line.material.isLineDashedMaterial) line.computeLineDistances();
  line.material.color.setHex(color);
  line.visible = true;
}

function updateCriticalPoints(points, lambda, eps, showManifolds) {
  const claimed = new Array(points.length).fill(false);

  for (const st of markerState) {
    if (!st.active) continue;
    let best = -1, bestD = MAX_MATCH_JUMP;
    for (let i = 0; i < points.length; i++) {
      if (claimed[i]) continue;
      const d = Math.hypot(points[i].x - st.x, points[i].y - st.y);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best >= 0) {
      st.x = points[best].x; st.y = points[best].y; st.type = points[best].type;
      claimed[best] = true;
    } else {
      st.active = false;
    }
  }

  for (let i = 0; i < points.length; i++) {
    if (claimed[i]) continue;
    const free = markerState.find((st) => !st.active);
    if (!free) continue;
    free.active = true;
    free.x = points[i].x; free.y = points[i].y; free.type = points[i].type;
  }

  for (const st of markerState) {
    if (!st.active) {
      st.mesh.visible = false;
      for (const b of st.branches) b.visible = false;
      continue;
    }
    const h = N(st.x, st.y, lambda, eps) * HEIGHT_SCALE;
    st.mesh.position.set(st.x, h + 0.03, st.y);
    st.mesh.material.color.setHex(critColors[st.type] ?? 0xffffff);
    st.mesh.visible = true;

    if (!showManifolds || st.type !== "S1") {
      for (const b of st.branches) b.visible = false;
      continue;
    }
    const { unstable, stable } = computeManifoldBranches(st.x, st.y, lambda, eps, points);
    setLineFromPath(st.branches[0], unstable[0].pts, lambda, eps, 0xf0a04b);
    setLineFromPath(st.branches[1], unstable[1].pts, lambda, eps, 0xf0a04b);
    setLineFromPath(st.branches[2], stable[0].pts, lambda, eps, 0x9db4d1);
    setLineFromPath(st.branches[3], stable[1].pts, lambda, eps, 0x9db4d1);
  }
}

// ---------------------------------------------------------------------------
// UI wiring
// ---------------------------------------------------------------------------

const slider = document.getElementById("lambda-slider");
const lambdaValueEl = document.getElementById("lambda-value");
const epsSlider = document.getElementById("eps-slider");
const epsValueEl = document.getElementById("eps-value");
const phaseLabelEl = document.getElementById("phase-label");
const phaseDetailEl = document.getElementById("phase-detail");
const countsEl = document.getElementById("counts-readout");
const playBtn = document.getElementById("play-btn");
const resetBtn = document.getElementById("reset-btn");
const manifoldToggle = document.getElementById("manifold-toggle");
const colorModeButtons = document.querySelectorAll("#color-mode .seg-btn");

let colorMode = "height";
let basinLabelsCache = null;
let basinCacheKey = null;

colorModeButtons.forEach((btn) => {
  btn.addEventListener("click", () => {
    colorMode = btn.dataset.mode;
    colorModeButtons.forEach((b) => b.classList.toggle("active", b === btn));
    refresh(parseFloat(slider.value), parseFloat(epsSlider.value));
  });
});
manifoldToggle.addEventListener("change", () => refresh(parseFloat(slider.value), parseFloat(epsSlider.value)));

function phaseFor(lambda, eps, pts) {
  const n1 = pts.filter((p) => p.type === "S1").length;
  const n2 = pts.filter((p) => p.type === "S2").length;
  let label = n2 > 0 ? `${n1}S₁+${n2}S₂` : `${n1}S₁`;
  if (Math.abs(eps) < 1e-6) {
    if (lambda < LAMBDA_MINUS) return { label, detail: "analytic prediction (λ < λ₋)" };
    if (lambda > LAMBDA_PLUS) return { label, detail: "analytic prediction (λ > λ₊)" };
    return { label, detail: "analytic prediction (λ₋ < λ < λ₊)" };
  }
  return { label, detail: `live-detected (ε=${eps.toFixed(3)}, cusp-unfolded)` };
}

function refresh(lambda, eps) {
  if (colorMode === "basin") {
    const key = `${lambda.toFixed(3)}|${eps.toFixed(3)}`;
    if (key !== basinCacheKey) {
      basinLabelsCache = computeBasinGrid(lambda, eps);
      basinCacheKey = key;
    }
  }
  updateSurface(lambda, eps, colorMode, basinLabelsCache);
  updateSourceLabels(eps);
  const pts = findCriticalPoints(lambda, eps);
  updateCriticalPoints(pts, lambda, eps, manifoldToggle.checked);

  lambdaValueEl.textContent = lambda.toFixed(3);
  epsValueEl.textContent = eps.toFixed(3);
  const phase = phaseFor(lambda, eps, pts);
  phaseLabelEl.textContent = phase.label;
  phaseDetailEl.textContent = phase.detail;

  const n1 = pts.filter((p) => p.type === "S1").length;
  const n2 = pts.filter((p) => p.type === "S2").length;
  const n0 = pts.filter((p) => p.type === "S0").length;
  let line = `${pts.length} critical points — ${n1}×S1 + ${n2}×S2`;
  if (n0) line += ` + ${n0}×S0`;
  countsEl.textContent = line;
}

slider.addEventListener("input", () => refresh(parseFloat(slider.value), parseFloat(epsSlider.value)));
epsSlider.addEventListener("input", () => {
  basinCacheKey = null;
  refresh(parseFloat(slider.value), parseFloat(epsSlider.value));
});

let playing = false;
let playDir = 1;
let playLambda = parseFloat(slider.value);
// Slow the sweep down near the two bifurcation thresholds. A fold/pitchfork
// is instantaneous in lambda: the critical-point count changes in a single
// step no matter how fine that step is, so this can't make the transition
// itself gradual — but a coarse fixed step (0.0025) can jump clean over the
// narrow window where two points are visibly close together right before
// annihilating, making it look like they vanish with no warning. Slowing
// down there lets that close-approach actually render.
const BIFURCATION_SLOWDOWN_WINDOW = 0.02;
function speedFactorNear(lambda) {
  const d = Math.min(Math.abs(lambda - LAMBDA_MINUS), Math.abs(lambda - LAMBDA_PLUS));
  if (d >= BIFURCATION_SLOWDOWN_WINDOW) return 1;
  return 0.12 + 0.88 * (d / BIFURCATION_SLOWDOWN_WINDOW);
}

function playStep() {
  if (!playing) return;
  // Advance a dedicated JS number, not the <input> element's own .value:
  // reading .value back every frame round-trips the number through the
  // browser's own decimal serialization (which drops precision the fixed
  // toFixed(4) string doesn't always round-trip exactly), and near a
  // bifurcation the per-frame step legitimately shrinks to a few
  // ten-thousandths — small enough that a lossy round trip can make two
  // consecutive frames serialize to the identical displayed value and the
  // sweep never numerically progresses past it.
  const step = 0.0025 * speedFactorNear(playLambda);
  playLambda += playDir * step;
  const max = parseFloat(slider.max), min = parseFloat(slider.min);
  if (playLambda >= max) { playLambda = max; playDir = -1; }
  if (playLambda <= min) { playLambda = min; playDir = 1; }
  slider.value = playLambda.toFixed(4);
  refresh(playLambda, parseFloat(epsSlider.value));
  requestAnimationFrame(playStep);
}
playBtn.addEventListener("click", () => {
  playing = !playing;
  playBtn.classList.toggle("active", playing);
  playBtn.textContent = playing ? "⏸ Pause" : "▶ Sweep λ";
  if (playing) {
    playLambda = parseFloat(slider.value);
    requestAnimationFrame(playStep);
  }
});

resetBtn.addEventListener("click", () => {
  camera.position.copy(DEFAULT_CAM_POS);
  controls.target.copy(DEFAULT_CAM_TARGET);
  controls.update();
});

// ---------------------------------------------------------------------------
// Resize + render loop
// ---------------------------------------------------------------------------

function resize() {
  const wrap = document.getElementById("canvas-wrap");
  const w = wrap.clientWidth, h = wrap.clientHeight;
  renderer.setSize(w, h, false);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
window.addEventListener("resize", resize);
resize();

function animate() {
  requestAnimationFrame(animate);
  controls.update();
  renderer.render(scene, camera);
}

refresh(parseFloat(slider.value), parseFloat(epsSlider.value));
animate();
