import { SafeAudition } from "./audio.mjs";
import { CrashEngine } from "./engine.mjs";
import { FitControls } from "./fit_controls.mjs";
import { matchedModelLevelDb } from "./level_match.mjs";
import { readReferences } from "./references.mjs";
import { ReferenceBrowser } from "./reference_browser.mjs";
import { SpectrogramView } from "./spectrogram.mjs";
import { downloadFit, readFit, snapshotState } from "./state.mjs";
import { Tooltips } from "./tooltips.mjs";

const byId = id => document.getElementById(id);
new Tooltips();
const state = {
  reference: null, references: [], synthesis: null, snapshots: [],
  macros: [], activeSnapshotId: null,
  event: {
    strength: 0.8, location: 0.8, hardness: 0.65, implement: 1,
    contactSpread: 0.2, seed: 17,
  },
  eventDefaults: {
    strength: 0.8, location: 0.8, hardness: 0.65, implement: 1,
    contactSpread: 0.2,
  },
  modelLevelTouched: false,
  analysis: {
    size: 2048, hop: 512, window: "hann", floorDb: -160,
    dynamicRangeDb: 90,
  },
};
const view = new SpectrogramView(byId("spectrogram"));
view.setSettings({
  mode: byId("view-mode").value,
  dynamicRangeDb: state.analysis.dynamicRangeDb,
});
const audition = new SafeAudition(setStatus);
const worker = new Worker("analysis_worker.mjs", { type: "module" });
const renderWorker = new Worker("render_worker.mjs", { type: "module" });
const generations = { reference: 0, synthesis: 0 };
const pendingAnalysis = new Set();
const referenceSpectrumCache = new Map();
const referenceSpectrumCacheLimit = 8;
let engine;
let renderTimer;
let renderGeneration = 0;
let renderInFlight = false;
let queuedRender = null;
let pendingLevelMatch = null;
let fitControls;
let referenceBrowser;
let restoringReference = false;

const implementFamilies = [
  {
    value: 0, label: "Bristle stiffness",
    endpoints: ["soft / fine", "stiff / coarse"],
  },
  {
    value: .5, label: "Mallet firmness",
    endpoints: ["soft", "hard"],
  },
  {
    value: 1, label: "Tip hardness",
    endpoints: ["soft / round", "hard / sharp"],
  },
];

function selectedImplement(value) {
  return implementFamilies.reduce((best, item) =>
    Math.abs(item.value - value) < Math.abs(best.value - value) ? item : best,
  );
}

function paintPerformanceControls() {
  const family = selectedImplement(state.event.implement);
  state.event.implement = family.value;
  document.querySelectorAll('input[name="implement"]').forEach(input => {
    input.checked = Number(input.value) === family.value;
  });
  byId("character-label").textContent = family.label;
  const endpoints = byId("character-endpoints").querySelectorAll("i");
  endpoints[0].textContent = family.endpoints[0];
  endpoints[1].textContent = family.endpoints[1];
  byId("hardness").value = state.event.hardness;
  byId("hardness").nextElementSibling.textContent =
    state.event.hardness.toFixed(2);
  byId("contact-spread").value = state.event.contactSpread;
  byId("contact-spread").nextElementSibling.textContent =
    state.event.contactSpread.toFixed(2);
}

function setStatus(message) { byId("status").textContent = message; }
function updateColourCeiling() {
  const peak = state.referenceSpectrum?.peakDb;
  byId("colour-ceiling").textContent = Number.isFinite(peak)
    ? `Reference ceiling ${peak.toFixed(1)} dBFS`
    : "Reference ceiling …";
}
function refreshModelLevelControl() {
  fitControls?.refresh("model_level_db");
}

function requestLevelMatch(resetToFactory = false) {
  if (!engine || !state.reference) return;
  if (resetToFactory) state.macros[0] = engine.macros[0].defaultValue;
  state.modelLevelTouched = false;
  pendingLevelMatch = { referenceId: state.reference.id };
  refreshModelLevelControl();
  audition.setMacros(state.macros);
  setStatus("Matching model level…");
  renderSynthesis();
}

function buildControls() {
  state.macros = engine.macros.map(item => item.defaultValue);
  fitControls = new FitControls({
    descriptors: engine.macros, state,
    onChange: key => {
      if (key === "model_level_db") {
        state.modelLevelTouched = true;
        pendingLevelMatch = null;
      }
      scheduleRender();
    },
    onLevelReset: () => requestLevelMatch(),
  });
  fitControls.build();
  applyBodyUiMode(fitControls.bodyMode(), false);
}

function applyBodyUiMode(mode, notify = true) {
  const selected = mode === "legacy" ? "legacy" : "unified";
  document.body.dataset.bodyMode = selected;
  byId("body-ui-mode").value = selected;
  fitControls?.setBodyMode(selected, notify);
}

function analyze(kind, samples, sampleRate, cacheKey) {
  const generation = ++generations[kind];
  pendingAnalysis.add(kind);
  setStatus(`Analyzing ${[...pendingAnalysis].join(" + ")}…`);
  worker.postMessage({
    generation, kind, samples: samples.slice(), sampleRate,
    settings: state.analysis, cacheKey,
  });
}

function referenceSpectrumKey(reference) {
  const source = reference.sha256 ?? reference.id;
  const { size, hop, window, floorDb } = state.analysis;
  return `${source}|${size}|${hop}|${window}|${floorDb}`;
}

function cacheReferenceSpectrum(key, spectrum) {
  referenceSpectrumCache.delete(key);
  referenceSpectrumCache.set(key, spectrum);
  while (referenceSpectrumCache.size > referenceSpectrumCacheLimit) {
    referenceSpectrumCache.delete(referenceSpectrumCache.keys().next().value);
  }
}

function analyzeReference() {
  if (!state.reference) return;
  const key = referenceSpectrumKey(state.reference);
  const cached = referenceSpectrumCache.get(key);
  if (!cached) {
    analyze(
      "reference", state.reference.samples, state.reference.sampleRate, key,
    );
    return;
  }
  ++generations.reference;
  pendingAnalysis.delete("reference");
  cacheReferenceSpectrum(key, cached);
  state.referenceSpectrum = cached;
  view.setData("reference", cached);
  updateColourCeiling();
  if (!pendingAnalysis.size) setStatus("Ready");
}

worker.onmessage = ({ data }) => {
  if (data.generation !== generations[data.kind]) return;
  pendingAnalysis.delete(data.kind);
  if (data.error) { setStatus(data.error); return; }
  view.setData(data.kind, data.result);
  state[`${data.kind}Spectrum`] = data.result;
  if (data.kind === "reference" && data.cacheKey) {
    cacheReferenceSpectrum(data.cacheKey, data.result);
    updateColourCeiling();
  }
  if (!pendingAnalysis.size) setStatus("Ready");
};

function scheduleRender() {
  clearTimeout(renderTimer);
  audition.setMacros(state.macros);
  setStatus("Rendering…");
  renderTimer = setTimeout(renderSynthesis, 220);
}

function renderSynthesis() {
  clearTimeout(renderTimer);
  const sampleRate = state.reference?.sampleRate ?? 48000;
  const duration = byId("render-seconds").value;
  const seconds = duration === "reference"
    ? (state.reference?.duration ?? 6) : Number(duration);
  const request = {
    generation: ++renderGeneration, sampleRate, seconds,
    macros: [...state.macros], event: { ...state.event },
  };
  if (renderInFlight) queuedRender = request;
  else dispatchRender(request);
}

function dispatchRender(request) {
  renderInFlight = true;
  renderWorker.postMessage(request);
}

renderWorker.onmessage = ({ data }) => {
  renderInFlight = false;
  if (data.generation === renderGeneration) {
    if (data.error) {
      pendingLevelMatch = null;
      setStatus(data.error);
    }
    else {
      if (pendingLevelMatch &&
          pendingLevelMatch.referenceId === state.reference?.id) {
        const descriptor = engine.macros[0];
        const matched = matchedModelLevelDb({
          currentDb: state.macros[0],
          reference: state.reference.samples,
          referenceSampleRate: state.reference.sampleRate,
          synthesis: data.samples,
          synthesisSampleRate: data.sampleRate,
          minimumDb: descriptor.minimum,
          maximumDb: descriptor.maximum,
          referenceStartSeconds: state.reference.cell?.onset_seconds ?? 0,
        });
        pendingLevelMatch = null;
        if (Math.abs(matched - state.macros[0]) > 0.01) {
          state.macros[0] = matched;
          refreshModelLevelControl();
          audition.setMacros(state.macros);
          queuedRender = null;
          setStatus("Applying matched model level…");
          renderSynthesis();
          return;
        }
      }
      state.synthesis = data.samples;
      analyze("synthesis", state.synthesis, data.sampleRate);
      drawWaveform();
      byId("render-time").textContent = `${data.elapsedMs.toFixed(0)} ms DSP`;
    }
  }
  if (queuedRender) {
    const request = queuedRender;
    queuedRender = null;
    dispatchRender(request);
  }
};

renderWorker.onerror = event => {
  renderInFlight = false;
  pendingLevelMatch = null;
  setStatus(event.message);
};

function decimate(samples, sampleRate, limit = 5000) {
  const step = Math.max(1, Math.floor(samples.length / limit));
  const x = []; const y = [];
  for (let index = 0; index < samples.length; index += step) {
    x.push(index / sampleRate); y.push(samples[index]);
  }
  return { x, y };
}

function drawWaveform() {
  const traces = [];
  if (state.reference) traces.push({
    ...decimate(state.reference.samples, state.reference.sampleRate),
    name: "Reference", mode: "lines", line: { color: "#e8b45a", width: 1 },
  });
  if (state.synthesis) traces.push({
    ...decimate(state.synthesis, state.reference?.sampleRate ?? 48000),
    name: "Synthesis", mode: "lines", line: { color: "#68a7d8", width: 1 },
  });
  Plotly.react("waveform", traces, {
    margin: { l: 44, r: 12, t: 8, b: 28 }, paper_bgcolor: "#0d1015",
    plot_bgcolor: "#0d1015", font: { color: "#94a0af", size: 10 },
    xaxis: {
      title: "seconds", gridcolor: "#252d38", zeroline: false,
      fixedrange: true,
    },
    yaxis: {
      title: "amplitude", gridcolor: "#252d38",
      zerolinecolor: "#46515e", fixedrange: true,
    },
    dragmode: false,
    legend: { orientation: "h", x: 0.75, y: 1 }, uirevision: "waveform-v1",
  }, {
    responsive: true, displaylogo: false, displayModeBar: false,
    scrollZoom: false,
  });
}

function setReference(reference) {
  const firstReference = !state.reference;
  state.reference = reference;
  fitControls?.decayEditor?.refresh();
  audition.setTrim(reference.corpus?.auditionTrimDb ?? 0);
  byId("calibration-trim").textContent =
    `Corpus trim ${audition.trimDb >= 0 ? "+" : ""}${audition.trimDb.toFixed(1)} dB`;
  if (reference.cell) {
    Object.assign(state.eventDefaults, {
      strength: reference.cell.strength,
      location: reference.cell.location,
      hardness: reference.cell.hardness,
      implement: reference.cell.implement ?? 1,
      contactSpread: reference.cell.contactSpread ?? .2,
    });
    Object.assign(state.event, {
      ...state.eventDefaults,
      seed: reference.cell.seed,
    });
    paintPerformanceControls();
  }
  if (firstReference) {
    byId("view-mode").value = "mirror";
    view.setSettings({ mode: "mirror" });
  }
  analyzeReference();
  drawWaveform(); view.reset();
  if (!restoringReference) {
    if (state.modelLevelTouched) scheduleRender();
    else requestLevelMatch(true);
  }
}

function renderSnapshotList() {
  const parent = byId("snapshots"); parent.replaceChildren();
  for (const item of state.snapshots) {
    const button = document.createElement("button");
    button.className = `snapshot-chip${item.fit.id === state.activeSnapshotId ? " active" : ""}`;
    button.textContent = item.fit.name;
    button.onclick = () => restore(item).catch(error => setStatus(String(error)));
    parent.append(button);
  }
}

async function restore(item) {
  const fit = item.fit;
  const sameReference = !fit.reference ||
    (fit.reference.sha256 && fit.reference.sha256 === state.reference?.sha256) ||
    fit.reference.id === state.reference?.id;
  if (!sameReference) {
    restoringReference = true;
    try {
      if (!await referenceBrowser?.selectSavedReference(fit.reference)) {
        throw new Error(
          "This fit's reference is not in the local corpus; load that WAV first",
        );
      }
    } finally {
      restoringReference = false;
    }
  }
  state.macros.splice(0, state.macros.length, ...fit.controls.macros);
  state.modelLevelTouched = true;
  pendingLevelMatch = null;
  Object.assign(state.event, fit.controls.event);
  Object.assign(state.analysis, fit.controls.analysis);
  state.activeSnapshotId = fit.id;
  if (item.audio) {
    state.synthesis = item.audio.slice();
    analyze(
      "synthesis", state.synthesis, fit.reference?.sampleRate ??
        state.reference?.sampleRate ?? 48000,
    );
  }
  buildPageValues();
  audition.setMacros(state.macros);
  drawWaveform();
  renderSnapshotList();
  if (!item.audio) scheduleRender();
}

function buildPageValues() {
  fitControls.build();
  applyBodyUiMode(fitControls.bodyMode(), false);
  paintPerformanceControls();
  byId("colour-range").value = state.analysis.dynamicRangeDb;
  byId("colour-range").nextElementSibling.textContent =
    `${state.analysis.dynamicRangeDb} dB`;
  view.setSettings({ dynamicRangeDb: state.analysis.dynamicRangeDb });
}

async function initialize() {
  engine = await CrashEngine.create();
  buildControls();
  audition.setMacros(state.macros);
  audition.initialize().catch(error => setStatus(String(error)));
  referenceBrowser = new ReferenceBrowser({
    corpus: byId("reference-corpus"),
    articulation: byId("reference-articulation"),
    velocity: byId("reference-velocity"),
    repeat: byId("reference-repeat"),
  }, setReference, setStatus);
  await referenceBrowser.initialize();
  byId("body-ui-mode").onchange = event =>
    applyBodyUiMode(event.currentTarget.value);
  byId("reference-files").onchange = async event => {
    try {
      const loaded = await readReferences(event.target.files);
      state.references.push(...loaded);
      if (loaded.length) setReference(loaded[0]);
    } catch (error) { setStatus(String(error)); }
  };
  byId("play-reference").onclick = () => state.reference &&
    audition.play(state.reference.samples, state.reference.sampleRate)
      .catch(error => setStatus(String(error)));
  byId("play-synthesis").onclick = () =>
    audition.trigger({ ...state.event }).catch(error => setStatus(String(error)));
  byId("stop").onclick = () => audition.stop();
  byId("master").oninput = event => {
    audition.setMaster(event.target.value);
    event.target.nextElementSibling.textContent = `${event.target.value} dB`;
  };
  byId("master").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = -12;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
  byId("size-meta").oninput = event => {
    const value = Number(event.target.value);
    event.target.nextElementSibling.textContent = value.toFixed(2);
    fitControls.applySizeMeta(value);
  };
  byId("size-meta").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = .5;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
  byId("hardness").oninput = event => {
    state.event.hardness = Number(event.target.value);
    event.target.nextElementSibling.textContent = state.event.hardness.toFixed(2);
    scheduleRender();
  };
  byId("hardness").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = state.eventDefaults.hardness;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
  byId("contact-spread").oninput = event => {
    state.event.contactSpread = Number(event.currentTarget.value);
    event.currentTarget.nextElementSibling.textContent =
      state.event.contactSpread.toFixed(2);
    scheduleRender();
  };
  byId("contact-spread").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = state.eventDefaults.contactSpread;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
  document.querySelectorAll('input[name="implement"]').forEach(input => {
    input.onchange = event => {
      state.event.implement = Number(event.currentTarget.value);
      paintPerformanceControls();
      scheduleRender();
    };
  });
  byId("strike-pad").onpointerdown = event => {
    const bounds = event.currentTarget.getBoundingClientRect();
    state.event.location = Math.max(0, Math.min(1,
      (event.clientX - bounds.left) / bounds.width));
    state.event.strength = Math.max(0.02, Math.min(1,
      1 - (event.clientY - bounds.top) / bounds.height));
    state.event.seed += 1;
    audition.trigger({ ...state.event }).catch(error => setStatus(String(error)));
  };
  bindAnalysisControls(); bindAnalysisDivider(); bindSnapshotControls();
  renderSynthesis();
  setInterval(() => {
    const latency = audition.latencyMs ? ` · ${audition.latencyMs.toFixed(0)} ms` : "";
    const underflows = audition.underflows ? ` · xruns ${audition.underflows}` : "";
    const output = Number.isFinite(audition.outputDb)
      ? ` · out ${audition.outputDb.toFixed(0)} dBFS` : "";
    const live = audition.state === "off" ? "" :
      ` · ${(audition.sampleRate / 1000).toFixed(1)} kHz` +
      ` · ${audition.state} · hits ${audition.triggerCount}`;
    byId("limiter").textContent =
      `Limiter ${audition.reduction.toFixed(1)} dB${latency}${output}${live}${underflows}`;
    byId("live-commit").textContent = audition.macroCommitPending
      ? `Preparing live DSP… ${audition.macroCommitElapsedMs.toFixed(0)} ms`
      : audition.macroCommitMs > 0
        ? `Live DSP ready · ${audition.macroCommitMs.toFixed(0)} ms`
        : "Live DSP idle";
  }, 100);
}

function bindAnalysisDivider() {
  const divider = byId("analysis-divider");
  const analysis = document.querySelector(".analysis");
  const storageKey = "tf-spectrogram-height-v2";
  const apply = height => {
    const maximum = Math.min(620, Math.max(220, .65 * window.innerHeight));
    const value = Math.round(Math.max(110, Math.min(maximum, height)));
    analysis.style.setProperty("--spectrogram-height", `${value}px`);
    divider.setAttribute("aria-valuenow", value);
    return value;
  };
  let height = apply(Number(localStorage.getItem(storageKey)) || 300);
  let drag = null;
  divider.onpointerdown = event => {
    if (event.button !== 0) return;
    event.preventDefault();
    drag = { pointerId: event.pointerId, y: event.clientY, height };
    divider.setPointerCapture(event.pointerId);
  };
  divider.onpointermove = event => {
    if (!drag || event.pointerId !== drag.pointerId) return;
    height = apply(drag.height + event.clientY - drag.y);
  };
  const finish = event => {
    if (!drag || event.pointerId !== drag.pointerId) return;
    drag = null;
    localStorage.setItem(storageKey, height);
    if (divider.hasPointerCapture(event.pointerId)) {
      divider.releasePointerCapture(event.pointerId);
    }
  };
  divider.onpointerup = finish;
  divider.onpointercancel = finish;
  divider.ondblclick = event => {
    event.preventDefault();
    height = apply(300);
    localStorage.setItem(storageKey, height);
  };
  divider.onkeydown = event => {
    if (!["ArrowUp", "ArrowDown", "Home"].includes(event.key)) return;
    event.preventDefault();
    height = apply(event.key === "Home" ? 300 :
      height + (event.key === "ArrowDown" ? 12 : -12));
    localStorage.setItem(storageKey, height);
  };
}

function bindAnalysisControls() {
  byId("view-mode").onchange = event => view.setSettings({ mode: event.target.value });
  const update = () => {
    state.analysis.size = Number(byId("fft-size").value);
    state.analysis.hop = Math.round(state.analysis.size *
      (1 - Number(byId("overlap").value)));
    state.analysis.window = byId("window").value;
    analyzeReference();
    if (state.synthesis) analyze("synthesis", state.synthesis,
      state.reference?.sampleRate ?? 48000);
  };
  byId("fft-size").onchange = update; byId("overlap").onchange = update;
  byId("window").onchange = update;
  byId("render-seconds").onchange = scheduleRender;
  byId("colour-range").oninput = event => {
    state.analysis.dynamicRangeDb = Number(event.target.value);
    event.target.nextElementSibling.textContent = `${event.target.value} dB`;
    view.setSettings({ dynamicRangeDb: state.analysis.dynamicRangeDb });
  };
  byId("reset-view").onclick = () => view.reset();
  const resetSelect = (id, value, eventName = "change") => {
    byId(id).ondblclick = event => {
      event.preventDefault();
      event.currentTarget.value = value;
      event.currentTarget.dispatchEvent(new Event(eventName));
    };
  };
  resetSelect("view-mode", "mirror");
  resetSelect("window", "hann");
  resetSelect("fft-size", "2048");
  resetSelect("overlap", "0.75");
  resetSelect("render-seconds", "reference");
  byId("colour-range").ondblclick = event => {
    event.preventDefault();
    event.currentTarget.value = 90;
    event.currentTarget.dispatchEvent(new Event("input"));
  };
}

function bindSnapshotControls() {
  byId("snapshot").onclick = () => {
    const fit = snapshotState(state, byId("snapshot-name").value || "Candidate");
    state.snapshots.push({ fit, audio: state.synthesis?.slice() });
    state.activeSnapshotId = fit.id; renderSnapshotList();
  };
  byId("save-fit").onclick = () => {
    const active = state.snapshots.find(item => item.fit.id === state.activeSnapshotId);
    downloadFit(active?.fit ?? snapshotState(state, byId("snapshot-name").value));
  };
  byId("load-fit").onchange = async event => {
    try {
      const fit = await readFit(event.target.files[0], engine.macros);
      const item = { fit }; state.snapshots.push(item); await restore(item);
    } catch (error) { setStatus(String(error)); }
  };
}

initialize().catch(error => setStatus(String(error)));
