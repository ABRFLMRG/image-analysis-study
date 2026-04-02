const GROUND_TRUTH_PATH = '../point_registration/ground_truth_coords.csv';
const PARTICIPANT_PATH = '../point_registration/coordinate_data_deidentified.csv';
const DATASET_ORDER = ['nuclei1', 'nuclei2', 'nuclei3', 'nuclei4', 'fish1', 'fish2', 'fish3', 'fish4'];

const DATASET_TO_GT_PATH = {
  nuclei1: '../../ground truth/nuclei/out_c00_dr90_label.tif',
  nuclei2: '../../ground truth/nuclei/out_c90_dr90_label.tif',
  nuclei3: '../../ground truth/nuclei/out_c00_dr10_label.tif',
  nuclei4: '../../ground truth/nuclei/out_c90_dr10_label.tif',
  fish1: '../../ground truth/FISH/celegans_dyn-90_ceff-0_label.ics.ome.tiff',
  fish2: '../../ground truth/FISH/celegans_dyn-90_ceff-90_label.ics.ome.tiff',
  fish3: '../../ground truth/FISH/celegans_dyn-10_ceff-0_label.ics.ome.tiff',
  fish4: '../../ground truth/FISH/celegans_dyn-10_ceff-90_label.ics.ome.tiff',
};
const GT_PATH_TO_DATASET = Object.fromEntries(
  Object.entries(DATASET_TO_GT_PATH).map(([dataset, path]) => [path, dataset])
);

const participantLabel = document.getElementById('participantLabel');
const datasetLabel = document.getElementById('datasetLabel');
const statusEl = document.getElementById('status');
const participantListEl = document.getElementById('participantList');
const datasetListEl = document.getElementById('datasetList');
const flipXYInput = document.getElementById('flipXY');
const xyScaleInput = document.getElementById('xyScale');
const zScaleInput = document.getElementById('zScale');
const matchThresholdInput = document.getElementById('matchThreshold');
const resetZoomBtn = document.getElementById('resetZoomBtn');
const computeMetricsBtn = document.getElementById('computeMetricsBtn');
const downloadSettingsBtn = document.getElementById('downloadSettingsBtn');
const uploadSettingsBtn = document.getElementById('uploadSettingsBtn');
const settingsFileInput = document.getElementById('settingsFileInput');
const metricsOutput = document.getElementById('metricsOutput');
const canvases = {
  xy: document.getElementById('xyPlot'),
  xz: document.getElementById('xzPlot'),
  yz: document.getElementById('yzPlot'),
};
const PLANE_CONFIG = {
  xy: { axisA: 'x', axisB: 'y' },
  xz: { axisA: 'x', axisB: 'z' },
  yz: { axisA: 'y', axisB: 'z' },
};

const groundTruthByDataset = new Map();
const participantData = new Map();
let participants = [];
let participantIndex = 0;
let datasetIndex = 0;
let flipXY = false;
let xyScale = 1;
let zScale = 1;
let matchThreshold = 0.2;
let zoomState = { xy: null, xz: null, yz: null };
let activeParticipantId = null;
let activeDatasetName = null;
const transformSettings = {};

setupControls();
setupParticipantList();
setupDatasetList();
setupCanvasInteractions();
setupActionButtons();
renderDatasetList();
init();

async function init() {
  try {
    const [gtText, participantText] = await Promise.all([
      fetch(GROUND_TRUTH_PATH).then((res) => res.text()),
      fetch(PARTICIPANT_PATH).then((res) => res.text()),
    ]);

    parseGroundTruth(gtText);

    parseParticipantData(participantText);

    participants = Array.from(participantData.keys()).sort();
    if (participants.length === 0) {
      statusEl.textContent = 'No participant data found.';
      return;
    }
    renderParticipantList();
    updateView();
    window.addEventListener('keydown', handleKeyPress);
  } catch (error) {
    console.error(error);
    statusEl.textContent = 'Failed to load data. Check console for details.';
  }
}

function parseGroundTruth(text) {
  const rows = parseCSV(text);
  for (const row of rows) {
    const datasetName = GT_PATH_TO_DATASET[row.path || row['path']];
    if (!datasetName) continue;
    const point = {
      x: Number(row.x ?? row['x']),
      y: Number(row.y ?? row['y']),
      z: Number(row.z ?? row['z']),
    };
    if (!Number.isFinite(point.x) || !Number.isFinite(point.y) || !Number.isFinite(point.z)) {
      continue;
    }
    if (!groundTruthByDataset.has(datasetName)) {
      groundTruthByDataset.set(datasetName, []);
    }
    groundTruthByDataset.get(datasetName).push(point);
  }
}

function parseParticipantData(text) {
  const rows = parseCSV(text);
  for (const row of rows) {
    const csvPath = row.csv_path || row['csv_path'];
    if (!csvPath) continue;
    const baseName = csvPath.split('/').pop()?.replace('.csv', '');
    if (!baseName) continue;

    const datasetName = DATASET_ORDER.find((name) => baseName.endsWith(`_${name}`));
    if (!datasetName) continue;

    const participant = baseName.slice(0, -(datasetName.length + 1));
    if (!participant) continue;

    const point = {
      x: Number(row.x ?? row['x']),
      y: Number(row.y ?? row['y']),
      z: Number(row.z ?? row['z']),
    };

    if (!Number.isFinite(point.x) || !Number.isFinite(point.y) || !Number.isFinite(point.z)) {
      continue;
    }

    if (!participantData.has(participant)) {
      participantData.set(participant, new Map());
    }

    const datasetMap = participantData.get(participant);
    if (!datasetMap.has(datasetName)) {
      datasetMap.set(datasetName, []);
    }

    datasetMap.get(datasetName).push(point);
  }
}

function parseCSV(text) {
  const lines = text.split(/\r?\n/).filter((line) => line.trim().length > 0);
  if (lines.length === 0) return [];
  const headers = lines[0].split(',').map((h) => h.trim());
  return lines.slice(1).map((line) => {
    const cells = line.split(',');
    const row = {};
    headers.forEach((header, index) => {
      row[header || String(index)] = cells[index]?.trim();
    });
    return row;
  });
}

function setupControls() {
  flipXYInput.addEventListener('change', () => {
    flipXY = flipXYInput.checked;
    saveCurrentSettings();
    updateView();
  });

  const handleScaleChange = () => {
    const nextValue = parsePositiveNumber(xyScaleInput.value, xyScale);
    if (nextValue !== xyScale) {
      xyScale = nextValue;
    }
    xyScaleInput.value = xyScale;
    saveCurrentSettings();
    updateView();
  };
  xyScaleInput.addEventListener('change', handleScaleChange);

  const handleZScaleChange = () => {
    const nextValue = parsePositiveNumber(zScaleInput.value, zScale);
    if (nextValue !== zScale) {
      zScale = nextValue;
    }
    zScaleInput.value = zScale;
    saveCurrentSettings();
    updateView();
  };
  zScaleInput.addEventListener('change', handleZScaleChange);

  matchThresholdInput.addEventListener('change', () => {
    const nextValue = parsePositiveNumber(matchThresholdInput.value, matchThreshold);
    if (nextValue !== matchThreshold) {
      matchThreshold = nextValue;
    }
    matchThresholdInput.value = matchThreshold;
  });
}

function setupParticipantList() {
  participantListEl.addEventListener('change', () => {
    const selected = participantListEl.value;
    const index = participants.indexOf(selected);
    if (index !== -1) {
      setParticipantIndex(index);
    }
  });
}

function renderParticipantList() {
  participantListEl.innerHTML = participants
    .map((id) => `<option value='${id}'>${id}</option>`)
    .join('');
  participantListEl.size = Math.min(12, Math.max(participants.length, 4));
  participantListEl.value = participants[participantIndex];
}

function setupDatasetList() {
  datasetListEl.addEventListener('change', () => {
    const selected = datasetListEl.value;
    const index = DATASET_ORDER.indexOf(selected);
    if (index !== -1) {
      setDatasetIndex(index);
    }
  });
}

function renderDatasetList() {
  datasetListEl.innerHTML = DATASET_ORDER
    .map((name, idx) => `<option value='${name}'>${idx + 1}. ${name}</option>`)
    .join('');
  datasetListEl.size = DATASET_ORDER.length;
  datasetListEl.value = DATASET_ORDER[datasetIndex];
}

function setParticipantIndex(index) {
  if (participants.length === 0) return;
  const normalized = ((index % participants.length) + participants.length) % participants.length;
  if (normalized === participantIndex) return;
  saveCurrentSettings();
  participantIndex = normalized;
  updateView();
}

function moveParticipant(delta) {
  setParticipantIndex(participantIndex + delta);
}

function setDatasetIndex(index) {
  const normalized = ((index % DATASET_ORDER.length) + DATASET_ORDER.length) % DATASET_ORDER.length;
  if (normalized === datasetIndex) return;
  saveCurrentSettings();
  datasetIndex = normalized;
  updateView();
}

function getSettingsBucket(participantId) {
  if (!transformSettings[participantId]) {
    transformSettings[participantId] = {};
  }
  return transformSettings[participantId];
}

function getSettingsForSelection(participantId, datasetName) {
  const bucket = getSettingsBucket(participantId);
  if (!bucket[datasetName]) {
    bucket[datasetName] = { flip: false, xyScale: 1, zScale: 1 };
  }
  return bucket[datasetName];
}

function applySettingsForSelection(participantId, datasetName) {
  const settings = getSettingsForSelection(participantId, datasetName);
  flipXY = settings.flip;
  xyScale = settings.xyScale;
  zScale = settings.zScale;
  flipXYInput.checked = flipXY;
  xyScaleInput.value = xyScale;
  zScaleInput.value = zScale;
}

function saveCurrentSettings() {
  if (!activeParticipantId || !activeDatasetName) return;
  const settings = getSettingsForSelection(activeParticipantId, activeDatasetName);
  settings.flip = flipXY;
  settings.xyScale = xyScale;
  settings.zScale = zScale;
}

function setupCanvasInteractions() {
  Object.entries(canvases).forEach(([key, canvas]) => {
    canvas.addEventListener('wheel', (event) => handleCanvasWheel(event, key));
  });
}

function setupActionButtons() {
  resetZoomBtn.addEventListener('click', () => {
    resetZoomState();
    updateView();
  });

  computeMetricsBtn.addEventListener('click', handleComputeMetrics);

  downloadSettingsBtn.addEventListener('click', downloadSettings);

  uploadSettingsBtn.addEventListener('click', () => {
    settingsFileInput.value = '';
    settingsFileInput.click();
  });

  settingsFileInput.addEventListener('change', handleSettingsUpload);
}

function parsePositiveNumber(value, fallback) {
  const parsed = Number(value);
  if (!Number.isFinite(parsed) || parsed <= 0) {
    return fallback;
  }
  return parsed;
}

function handleKeyPress(event) {
  if (participants.length === 0 || isEditingInput(event.target)) return;
  if (event.shiftKey && /^[1-8]$/.test(event.key)) {
    const numericIndex = Number(event.key) - 1;
    if (numericIndex >= 0 && numericIndex < DATASET_ORDER.length) {
      event.preventDefault();
      setDatasetIndex(numericIndex);
    }
    return;
  }

  switch (event.key) {
    case 'ArrowUp':
      event.preventDefault();
      moveParticipant(-1);
      break;
    case 'ArrowDown':
      event.preventDefault();
      moveParticipant(1);
      break;
    default:
      break;
  }
}

function updateView(options = {}) {
  const { forceApplySettings = false } = options;
  if (participants.length === 0) return;
  const participantId = participants[participantIndex];
  const datasetName = DATASET_ORDER[datasetIndex];
  const participantChanged = participantId !== activeParticipantId;
  const datasetChanged = datasetName !== activeDatasetName;

  if (participantChanged) {
    resetZoomState();
  }

  if (forceApplySettings || participantChanged || datasetChanged || activeParticipantId === null) {
    applySettingsForSelection(participantId, datasetName);
  }

  activeParticipantId = participantId;
  activeDatasetName = datasetName;

  participantLabel.textContent = participantId;
  datasetLabel.textContent = datasetName;
  participantListEl.value = participantId;
  datasetListEl.value = datasetName;

  const participantDatasets = participantData.get(participantId);
  const selectedPoints = participantDatasets?.get(datasetName) ?? [];
  const groundTruthPoints = groundTruthByDataset.get(datasetName) ?? [];

  if (!participantDatasets?.has(datasetName)) {
    statusEl.textContent = `No ${datasetName} data for ${participantId}.`;
  } else {
    statusEl.textContent = `${selectedPoints.length} points shown for ${participantId} / ${datasetName} (ground truth: ${groundTruthPoints.length}).`;
  }

  drawPlots(selectedPoints, groundTruthPoints);
}

function drawPlots(datasetPoints, groundTruthPoints) {
  const scaledDataset = prepareDatasetPoints(datasetPoints);
  drawPlane('xy', scaledDataset, groundTruthPoints);
  drawPlane('xz', scaledDataset, groundTruthPoints);
  drawPlane('yz', scaledDataset, groundTruthPoints);
}

function scalePoint(point) {
  return {
    x: point.x * xyScale,
    y: point.y * xyScale,
    z: point.z * zScale,
  };
}

function prepareDatasetPoints(points) {
  const transformedDataset = flipXY ? swapAxes(points, 'x', 'y') : points;
  return transformedDataset.map(scalePoint);
}

function drawPlane(planeKey, datasetPoints, groundTruthPoints) {
  const { axisA, axisB } = PLANE_CONFIG[planeKey];
  const canvas = canvases[planeKey];
  const ctx = canvas.getContext('2d');
  ctx.clearRect(0, 0, canvas.width, canvas.height);

  const gtPoints = groundTruthPoints;
  const dsPoints = datasetPoints;

  const points = [...gtPoints, ...dsPoints];
  if (points.length === 0) {
    ctx.fillStyle = '#94a3b8';
    ctx.fillText('No data', 10, 20);
    return;
  }

  const margins = 30;
  const computedBounds = computeBounds(points, axisA, axisB);
  const viewport = ensureViewport(planeKey, computedBounds);
  const { bounds } = viewport;
  const minA = bounds.minA;
  const maxA = bounds.maxA;
  const minB = bounds.minB;
  const maxB = bounds.maxB;

  const rangeA = maxA - minA || 1;
  const rangeB = maxB - minB || 1;

  const scale = (value, min, range, size) => ((value - min) / range) * (size - margins * 2) + margins;

  // axis grid
  ctx.strokeStyle = '#e2e8f0';
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.rect(margins, margins, canvas.width - margins * 2, canvas.height - margins * 2);
  ctx.stroke();

  drawHorizontalTicks(
    ctx,
    generateTicks(minA, maxA),
    (value) => scale(value, minA, rangeA, canvas.width),
    canvas.height,
    margins,
    rangeA
  );
  drawVerticalTicks(
    ctx,
    generateTicks(minB, maxB),
    (value) => canvas.height - scale(value, minB, rangeB, canvas.height),
    margins,
    rangeB
  );

  const drawPoints = (pointsToDraw, color, size) => {
    ctx.fillStyle = color;
    for (const point of pointsToDraw) {
      const x = scale(point[axisA], minA, rangeA, canvas.width);
      const y = canvas.height - scale(point[axisB], minB, rangeB, canvas.height);
      ctx.beginPath();
      ctx.arc(x, y, size, 0, Math.PI * 2);
      ctx.fill();
    }
  };

  drawPoints(gtPoints, 'rgba(59, 130, 246, 0.7)', 3);
  drawPoints(dsPoints, 'rgba(239, 68, 68, 0.9)', 4);
}

function swapAxes(points, axisA, axisB) {
  return points.map((point) => ({
    ...point,
    [axisA]: point[axisB],
    [axisB]: point[axisA],
  }));
}

function computeBounds(points, axisA, axisB) {
  const valuesA = points.map((p) => p[axisA]);
  const valuesB = points.map((p) => p[axisB]);
  const minA = Math.min(...valuesA);
  const maxA = Math.max(...valuesA);
  const minB = Math.min(...valuesB);
  const maxB = Math.max(...valuesB);
  return { minA, maxA, minB, maxB };
}

function ensureViewport(planeKey, bounds) {
  const existing = zoomState[planeKey];
  if (!existing) {
    zoomState[planeKey] = {
      bounds: { ...bounds },
      baseBounds: { ...bounds },
    };
    return zoomState[planeKey];
  }

  const ratios = deriveViewportRatios(existing.bounds, existing.baseBounds);
  existing.baseBounds = { ...bounds };
  existing.bounds = applyViewportRatios(ratios, existing.baseBounds);
  return existing;
}

function resetZoomState() {
  zoomState = { xy: null, xz: null, yz: null };
}

function handleCanvasWheel(event, planeKey) {
  const viewport = zoomState[planeKey];
  if (!viewport) return;
  event.preventDefault();
  const rect = event.currentTarget.getBoundingClientRect();
  const relX = clamp((event.clientX - rect.left) / rect.width, 0, 1);
  const relY = clamp((event.clientY - rect.top) / rect.height, 0, 1);
  const zoomFactor = Math.exp(-event.deltaY * 0.001);
  zoomAxis(viewport, 'A', relX, zoomFactor);
  // invert Y ratio since canvas origin top-left
  zoomAxis(viewport, 'B', 1 - relY, zoomFactor);
  updateView();
}

function zoomAxis(viewport, axisKey, ratio, zoomFactor) {
  const minKey = axisKey === 'A' ? 'minA' : 'minB';
  const maxKey = axisKey === 'A' ? 'maxA' : 'maxB';
  const bounds = viewport.bounds;
  const base = viewport.baseBounds;
  const baseRange = base[maxKey] - base[minKey];
  if (baseRange === 0) {
    bounds[minKey] = base[minKey];
    bounds[maxKey] = base[maxKey];
    return;
  }
  const currentRange = bounds[maxKey] - bounds[minKey] || 1;
  const minRange = baseRange * 0.01;
  const desiredRange = clamp(currentRange / zoomFactor, minRange, baseRange);
  const focus = bounds[minKey] + ratio * currentRange;
  let newMin = focus - ratio * desiredRange;
  let newMax = newMin + desiredRange;

  if (newMin < base[minKey]) {
    const diff = base[minKey] - newMin;
    newMin += diff;
    newMax += diff;
  }
  if (newMax > base[maxKey]) {
    const diff = newMax - base[maxKey];
    newMin -= diff;
    newMax -= diff;
  }

  newMin = clamp(newMin, base[minKey], base[maxKey] - minRange);
  newMax = clamp(newMax, newMin + minRange, base[maxKey]);

  bounds[minKey] = newMin;
  bounds[maxKey] = newMax;
}

function clamp(value, min, max) {
  return Math.min(Math.max(value, min), max);
}

function isEditingInput(element) {
  if (!element) return false;
  if (element.isContentEditable) return true;
  const tag = element.tagName?.toLowerCase();
  return tag === 'input' || tag === 'textarea' || tag === 'select';
}

function generateTicks(min, max, maxTicks = 5) {
  if (!Number.isFinite(min) || !Number.isFinite(max)) return [];
  if (max - min === 0) return [min];
  const rawStep = Math.abs(max - min) / Math.max(1, maxTicks);
  const magnitude = Math.pow(10, Math.floor(Math.log10(rawStep)));
  const normalized = rawStep / magnitude;
  let step;
  if (normalized < 1.5) step = 1;
  else if (normalized < 3) step = 2;
  else if (normalized < 7) step = 5;
  else step = 10;
  step *= magnitude;
  const ticks = [];
  let tick = Math.ceil(min / step) * step;
  while (tick <= max + step / 2) {
    ticks.push(Number.parseFloat(tick.toFixed(6)));
    tick += step;
  }
  return ticks;
}

function drawHorizontalTicks(ctx, ticks, toX, canvasHeight, margins, range) {
  ctx.save();
  ctx.strokeStyle = '#94a3b8';
  ctx.fillStyle = '#475467';
  ctx.font = '12px Inter, sans-serif';
  ctx.textAlign = 'center';
  ctx.textBaseline = 'top';
  const axisY = canvasHeight - margins;
  ticks.forEach((value) => {
    const x = toX(value);
    ctx.beginPath();
    ctx.moveTo(x, axisY);
    ctx.lineTo(x, axisY + 6);
    ctx.stroke();
    ctx.fillText(`${formatTickLabel(value, range)} um`, x, axisY + 8);
  });
  ctx.restore();
}

function drawVerticalTicks(ctx, ticks, toY, margins, range) {
  ctx.save();
  ctx.strokeStyle = '#94a3b8';
  ctx.fillStyle = '#475467';
  ctx.font = '12px Inter, sans-serif';
  ctx.textAlign = 'right';
  ctx.textBaseline = 'middle';
  const axisX = margins;
  ticks.forEach((value) => {
    const y = toY(value);
    ctx.beginPath();
    ctx.moveTo(axisX - 6, y);
    ctx.lineTo(axisX, y);
    ctx.stroke();
    ctx.fillText(`${formatTickLabel(value, range)} um`, axisX - 8, y);
  });
  ctx.restore();
}

function formatTickLabel(value, range) {
  const magnitude = Math.abs(range);
  let decimals = 0;
  if (magnitude < 10) decimals = 2;
  else if (magnitude < 100) decimals = 1;
  return value.toFixed(decimals);
}

function deriveViewportRatios(bounds, base) {
  return {
    axisA: deriveAxisRatio(bounds.minA, bounds.maxA, base.minA, base.maxA),
    axisB: deriveAxisRatio(bounds.minB, bounds.maxB, base.minB, base.maxB),
  };
}

function deriveAxisRatio(min, max, baseMin, baseMax) {
  const baseRange = baseMax - baseMin;
  if (!Number.isFinite(baseRange) || baseRange === 0) {
    return { center: 0.5, size: 1 };
  }
  const clampedMin = clamp(min, baseMin, baseMax);
  const clampedMax = clamp(max, baseMin, baseMax);
  const range = Math.max(clampedMax - clampedMin, baseRange * 0.01);
  const center = ((clampedMin + clampedMax) / 2 - baseMin) / baseRange;
  const size = range / baseRange;
  return {
    center: clamp(center, 0, 1),
    size: clamp(size, 0.01, 1),
  };
}

function applyViewportRatios(ratios, base) {
  const axisA = applyAxisRatio(ratios.axisA, base.minA, base.maxA);
  const axisB = applyAxisRatio(ratios.axisB, base.minB, base.maxB);
  return {
    minA: axisA.min,
    maxA: axisA.max,
    minB: axisB.min,
    maxB: axisB.max,
  };
}

function applyAxisRatio(ratio, baseMin, baseMax) {
  const baseRange = baseMax - baseMin;
  if (!Number.isFinite(baseRange) || baseRange === 0) {
    return { min: baseMin, max: baseMax };
  }
  const size = clamp(ratio.size, 0.01, 1) * baseRange;
  const half = size / 2;
  let center = baseMin + clamp(ratio.center, 0, 1) * baseRange;
  let min = center - half;
  let max = center + half;
  if (min < baseMin) {
    const diff = baseMin - min;
    min += diff;
    max += diff;
  }
  if (max > baseMax) {
    const diff = max - baseMax;
    min -= diff;
    max -= diff;
  }
  min = clamp(min, baseMin, baseMax - baseRange * 0.01);
  max = clamp(max, min + baseRange * 0.01, baseMax);
  return { min, max };
}

async function handleComputeMetrics() {
  if (!activeParticipantId || !activeDatasetName) return;
  const participantDatasets = participantData.get(activeParticipantId);
  const datasetPoints = participantDatasets?.get(activeDatasetName) ?? [];
  const groundTruthPoints = groundTruthByDataset.get(activeDatasetName) ?? [];
  if (datasetPoints.length === 0 || groundTruthPoints.length === 0) {
    metricsOutput.textContent = 'Need both participant and ground-truth points.';
    return;
  }
  metricsOutput.textContent = 'Computing assignment...';
  await new Promise((resolve) => requestAnimationFrame(resolve));
  try {
    const transformed = prepareDatasetPoints(datasetPoints);
    const { mse, jaccard, matchedCount } = computeAssignmentMetrics(transformed, groundTruthPoints);
    if (!matchedCount) {
      metricsOutput.textContent = 'No matches found.';
      return;
    }
    const mseLabel = Number.isFinite(mse) ? `${mse.toFixed(3)} um^2` : 'n/a';
    const jaccardLabel = `${(jaccard * 100).toFixed(2)}%`;
    metricsOutput.textContent = `Matches: ${matchedCount} · MSE: ${mseLabel} · Jaccard: ${jaccardLabel}`;
  } catch (error) {
    console.error(error);
    metricsOutput.textContent = 'Assignment failed. See console for details.';
  }
}

function computeAssignmentMetrics(datasetPoints, groundTruthPoints) {
  const data = datasetPoints.map((p) => [p.x, p.y, p.z]);
  const gt = groundTruthPoints.map((p) => [p.x, p.y, p.z]);
  const size = Math.max(data.length, gt.length);
  const penalty = 1e12;
  const thresholdSq = (matchThreshold || 0.2) ** 2;
  const costMatrix = Array.from({ length: size }, (_, row) => {
    return Array.from({ length: size }, (_, col) => {
      if (row < data.length && col < gt.length) {
        const dx = data[row][0] - gt[col][0];
        const dy = data[row][1] - gt[col][1];
        const dz = data[row][2] - gt[col][2];
        const distSq = dx * dx + dy * dy + dz * dz;
        return distSq <= thresholdSq ? distSq : penalty;
      }
      return penalty;
    });
  });

  const assignment = hungarian(costMatrix);
  let matchedCount = 0;
  let totalCost = 0;
  assignment.forEach((col, row) => {
    if (row < data.length && col < gt.length) {
      const cost = costMatrix[row][col];
      if (cost < penalty / 2) {
        matchedCount += 1;
        totalCost += cost;
      }
    }
  });

  const mse = matchedCount ? totalCost / matchedCount : Number.NaN;
  const jaccard = matchedCount ? matchedCount / (data.length + gt.length - matchedCount) : 0;
  return { mse, jaccard, matchedCount };
}

function hungarian(costMatrix) {
  const n = costMatrix.length;
  const m = costMatrix[0]?.length || 0;
  const size = Math.max(n, m);
  const cost = Array.from({ length: size }, (_, i) => (
    Array.from({ length: size }, (_, j) => (i < n && j < m ? costMatrix[i][j] : 1e12))
  ));
  const u = Array(size + 1).fill(0);
  const v = Array(size + 1).fill(0);
  const p = Array(size + 1).fill(0);
  const way = Array(size + 1).fill(0);
  for (let i = 1; i <= size; i += 1) {
    p[0] = i;
    let j0 = 0;
    const minv = Array(size + 1).fill(Infinity);
    const used = Array(size + 1).fill(false);
    do {
      used[j0] = true;
      const i0 = p[j0];
      let delta = Infinity;
      let j1 = 0;
      for (let j = 1; j <= size; j += 1) {
        if (used[j]) continue;
        const cur = cost[i0 - 1][j - 1] - u[i0] - v[j];
        if (cur < minv[j]) {
          minv[j] = cur;
          way[j] = j0;
        }
        if (minv[j] < delta) {
          delta = minv[j];
          j1 = j;
        }
      }
      for (let j = 0; j <= size; j += 1) {
        if (used[j]) {
          u[p[j]] += delta;
          v[j] -= delta;
        } else {
          minv[j] -= delta;
        }
      }
      j0 = j1;
    } while (p[j0] !== 0);
    do {
      const j1 = way[j0];
      p[j0] = p[j1];
      j0 = j1;
    } while (j0 !== 0);
  }
  const assignment = Array(size).fill(-1);
  for (let j = 1; j <= size; j += 1) {
    if (p[j] > 0 && p[j] <= size) {
      assignment[p[j] - 1] = j - 1;
    }
  }
  return assignment;
}

function downloadSettings() {
  const rows = ['participant,dataset,flip_xy,xy_scale,z_scale'];
  Object.entries(transformSettings).forEach(([participant, datasets]) => {
    Object.entries(datasets).forEach(([datasetName, settings]) => {
      rows.push(
        `${participant},${datasetName},${settings.flip ? 1 : 0},${settings.xyScale},${settings.zScale}`
      );
    });
  });
  const blob = new Blob([rows.join('\n')], { type: 'text/csv' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  link.href = url;
  link.download = 'dataset-settings.csv';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
  URL.revokeObjectURL(url);
}

function handleSettingsUpload(event) {
  const file = event.target.files?.[0];
  if (!file) return;
  const reader = new FileReader();
  reader.onload = () => {
    importSettingsFromCSV(reader.result);
  };
  reader.readAsText(file);
}

function importSettingsFromCSV(text = '') {
  const lines = text.split(/\r?\n/).filter((line) => line.trim().length > 0);
  if (lines.length === 0) return;
  const header = lines[0].split(',').map((h) => h.trim().toLowerCase());
  const colIndex = {
    participant: header.indexOf('participant'),
    dataset: header.indexOf('dataset'),
    flip: header.findIndex((h) => h === 'flip' || h === 'flip_xy'),
    xy: header.findIndex((h) => h.includes('xy')),
    z: header.findIndex((h) => h.includes('z')),
  };
  if (colIndex.participant === -1) colIndex.participant = 0;
  if (colIndex.dataset === -1) colIndex.dataset = 1;
  if (colIndex.flip === -1) colIndex.flip = 2;
  if (colIndex.xy === -1) colIndex.xy = 3;
  if (colIndex.z === -1) colIndex.z = 4;
  for (let i = 1; i < lines.length; i += 1) {
    const cells = lines[i].split(',');
    const participantId = cells[colIndex.participant]?.trim();
    const datasetName = cells[colIndex.dataset]?.trim();
    if (!participantId || !datasetName) continue;
    if (!DATASET_ORDER.includes(datasetName)) continue;
    const settings = getSettingsForSelection(participantId, datasetName);
    const flipValue = cells[colIndex.flip]?.trim().toLowerCase();
    settings.flip = flipValue === '1' || flipValue === 'true';
    const xyValue = Number(cells[colIndex.xy]);
    if (Number.isFinite(xyValue) && xyValue > 0) {
      settings.xyScale = xyValue;
    }
    const zValue = Number(cells[colIndex.z]);
    if (Number.isFinite(zValue) && zValue > 0) {
      settings.zScale = zValue;
    }
  }
  if (activeParticipantId && activeDatasetName) {
    applySettingsForSelection(activeParticipantId, activeDatasetName);
    updateView({ forceApplySettings: true });
    saveCurrentSettings();
  }
}
