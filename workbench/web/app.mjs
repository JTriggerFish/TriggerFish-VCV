import { SafeAudition } from "./audio.mjs";
import { CrashEngine } from "./engine.mjs";
import { matchedModelLevelDb } from "./level_match.mjs";
import { readReferences } from "./references.mjs";
import { ReferenceBrowser } from "./reference_browser.mjs";
import { SpectrogramView } from "./spectrogram.mjs";
import { downloadFit, readFit, snapshotState } from "./state.mjs";

const byId = id => document.getElementById(id);
const state = {
  reference: null, references: [], synthesis: null, snapshots: [],
  macros: [], activeSnapshotId: null,
  event: { strength: 0.8, location: 0.8, hardness: 0.65, seed: 17 },
  eventDefaults: { strength: 0.8, location: 0.8, hardness: 0.65 },
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

function setStatus(message) { byId("status").textContent = message; }
function valueText(descriptor, value) {
  return `${Number(value).toFixed(descriptor.unit === "dB" ? 1 : 3)}${descriptor.unit ? ` ${descriptor.unit}` : ""}`;
}

function refreshModelLevelControl() {
  const descriptor = engine.macros[0];
  const parent = byId("model-level");
  const input = parent.querySelector("input");
  const output = parent.querySelector("output");
  if (input) input.value = state.macros[0];
  if (output) output.textContent = valueText(descriptor, state.macros[0]);
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

function slider(descriptor, parent, onInput) {
  const row = document.createElement("label");
  row.className = "slider-row";
  row.append(descriptor.name);
  const input = document.createElement("input");
  input.type = "range"; input.min = descriptor.minimum; input.max = descriptor.maximum;
  input.step = descriptor.unit === "dB"
    ? 0.1 : (descriptor.maximum - descriptor.minimum) / 400;
  input.value = state.macros[descriptor.index];
  const output = document.createElement("output");
  output.textContent = valueText(descriptor, input.value);
  input.oninput = () => {
    state.macros[descriptor.index] = Number(input.value);
    if (descriptor.index === 0) {
      state.modelLevelTouched = true;
      pendingLevelMatch = null;
    }
    output.textContent = valueText(descriptor, input.value);
    onInput?.();
  };
  input.ondblclick = event => {
    event.preventDefault();
    if (descriptor.index === 0) {
      requestLevelMatch();
      return;
    }
    state.macros[descriptor.index] = descriptor.defaultValue;
    input.value = descriptor.defaultValue;
    output.textContent = valueText(descriptor, descriptor.defaultValue);
    onInput?.();
  };
  row.append(input, output); parent.append(row);
}

function pad(first, second, parent) {
  const row = document.createElement("div"); row.className = "xy";
  const pad = document.createElement("div"); pad.className = "xy-pad";
  const labels = document.createElement("div"); labels.className = "xy-labels";
  labels.innerHTML = `<b>${first.name}</b><span>${second.name}</span>`;
  const update = event => {
    const bounds = pad.getBoundingClientRect();
    const x = Math.max(0, Math.min(1, (event.clientX - bounds.left) / bounds.width));
    const y = Math.max(0, Math.min(1, (event.clientY - bounds.top) / bounds.height));
    state.macros[first.index] = first.minimum + x * (first.maximum - first.minimum);
    state.macros[second.index] = second.minimum + (1 - y) *
      (second.maximum - second.minimum);
    paint(); scheduleRender();
  };
  const paint = () => {
    const x = (state.macros[first.index] - first.minimum) /
      (first.maximum - first.minimum);
    const y = 1 - (state.macros[second.index] - second.minimum) /
      (second.maximum - second.minimum);
    pad.style.setProperty("--x", `${100 * x}%`);
    pad.style.setProperty("--y", `${100 * y}%`);
  };
  pad.onpointerdown = event => {
    pad.setPointerCapture(event.pointerId); update(event);
    pad.onpointermove = update;
  };
  pad.onpointerup = () => { pad.onpointermove = null; };
  pad.ondblclick = event => {
    event.preventDefault();
    state.macros[first.index] = first.defaultValue;
    state.macros[second.index] = second.defaultValue;
    paint(); scheduleRender();
  };
  paint(); row.append(pad, labels); parent.append(row);
}

function buildControls() {
  state.macros = engine.macros.map(item => item.defaultValue);
  slider(engine.macros[0], byId("model-level"), scheduleRender);
  pad(engine.macros[1], engine.macros[2], byId("macro-pads"));
  pad(engine.macros[3], engine.macros[4], byId("macro-pads"));
  pad(engine.macros[5], engine.macros[6], byId("macro-pads"));
  for (const index of [7, 8, 9]) {
    slider(engine.macros[index], byId("decay-controls"), scheduleRender);
  }
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
    xaxis: { title: "seconds", gridcolor: "#252d38", zeroline: false },
    yaxis: { title: "amplitude", gridcolor: "#252d38", zerolinecolor: "#46515e" },
    legend: { orientation: "h", x: 0.75, y: 1 }, uirevision: "waveform-v1",
  }, { responsive: true, displaylogo: false, scrollZoom: true });
}

function setReference(reference) {
  const firstReference = !state.reference;
  state.reference = reference;
  audition.setTrim(reference.corpus?.auditionTrimDb ?? 0);
  byId("calibration-trim").textContent =
    `Corpus trim ${audition.trimDb >= 0 ? "+" : ""}${audition.trimDb.toFixed(1)} dB`;
  if (reference.cell) {
    Object.assign(state.eventDefaults, {
      strength: reference.cell.strength,
      location: reference.cell.location,
      hardness: reference.cell.hardness,
    });
    Object.assign(state.event, {
      ...state.eventDefaults,
      seed: reference.cell.seed,
    });
    byId("hardness").value = state.event.hardness;
    byId("hardness").nextElementSibling.textContent = state.event.hardness.toFixed(2);
  }
  if (firstReference) {
    byId("view-mode").value = "mirror";
    view.setSettings({ mode: "mirror" });
  }
  analyzeReference();
  drawWaveform(); view.reset();
  if (state.modelLevelTouched) scheduleRender();
  else requestLevelMatch(true);
}

function renderSnapshotList() {
  const parent = byId("snapshots"); parent.replaceChildren();
  for (const item of state.snapshots) {
    const button = document.createElement("button");
    button.className = `snapshot-chip${item.fit.id === state.activeSnapshotId ? " active" : ""}`;
    button.textContent = item.fit.name;
    button.onclick = () => restore(item);
    parent.append(button);
  }
}

function restore(item) {
  const fit = item.fit;
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
  document.querySelectorAll("#model-level,#macro-pads,#decay-controls")
    .forEach(element => element.replaceChildren());
  buildControlsFromState();
  byId("hardness").value = state.event.hardness;
  byId("hardness").nextElementSibling.textContent = state.event.hardness.toFixed(2);
  byId("colour-range").value = state.analysis.dynamicRangeDb;
  byId("colour-range").nextElementSibling.textContent =
    `${state.analysis.dynamicRangeDb} dB`;
  view.setSettings({ dynamicRangeDb: state.analysis.dynamicRangeDb });
}

function buildControlsFromState() {
  slider(engine.macros[0], byId("model-level"), scheduleRender);
  pad(engine.macros[1], engine.macros[2], byId("macro-pads"));
  pad(engine.macros[3], engine.macros[4], byId("macro-pads"));
  pad(engine.macros[5], engine.macros[6], byId("macro-pads"));
  for (const index of [7, 8, 9]) slider(
    engine.macros[index], byId("decay-controls"), scheduleRender,
  );
}

async function initialize() {
  engine = await CrashEngine.create();
  buildControls();
  audition.setMacros(state.macros);
  audition.initialize().catch(error => setStatus(String(error)));
  const referenceBrowser = new ReferenceBrowser({
    corpus: byId("reference-corpus"),
    articulation: byId("reference-articulation"),
    velocity: byId("reference-velocity"),
    repeat: byId("reference-repeat"),
  }, setReference, setStatus);
  await referenceBrowser.initialize();
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
  byId("strike-pad").onpointerdown = event => {
    const bounds = event.currentTarget.getBoundingClientRect();
    state.event.location = Math.max(0, Math.min(1,
      (event.clientX - bounds.left) / bounds.width));
    state.event.strength = Math.max(0.02, Math.min(1,
      1 - (event.clientY - bounds.top) / bounds.height));
    state.event.seed += 1;
    audition.trigger({ ...state.event }).catch(error => setStatus(String(error)));
  };
  bindAnalysisControls(); bindSnapshotControls();
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
  }, 100);
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
      const fit = await readFit(event.target.files[0]);
      const item = { fit }; state.snapshots.push(item); restore(item);
    } catch (error) { setStatus(String(error)); }
  };
}

initialize().catch(error => setStatus(String(error)));
