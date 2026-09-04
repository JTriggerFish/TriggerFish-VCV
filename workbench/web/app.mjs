import { SafeAudition } from "./audio.mjs";
import { bindAnalysisControls } from "./analysis_controls.mjs";
import { PercussionEngine } from "./engine.mjs";
import { FitControls } from "./fit_controls.mjs";
import {
  calibrationParameterValues, calibrationPatch,
} from "./instrument_calibrations.mjs";
import { KickControls } from "./kick_controls.mjs";
import { MembraneControls } from "./membrane_controls.mjs";
import {
  createAcousticKickPatch, createTomPatch, membranePresetValues,
} from "./membrane_patch.mjs";
import { SnareControls } from "./snare_controls.mjs";
import { modelLevelMatch } from "./level_match.mjs";
import { PerformanceControls } from "./performance_controls.mjs";
import { recipeAdapter } from "./recipe_adapter.mjs";
import { RecipeController } from "./recipe_controller.mjs";
import { readReferences } from "./references.mjs";
import { ReferenceBrowser } from "./reference_browser.mjs";
import { RoutingController } from "./routing_controller.mjs";
import { SettingsController } from "./settings.mjs";
import {
  alignedReferenceWindow, SpectrogramView,
} from "./spectrogram.mjs";
import {
  downloadFit, fitMacroValues, readFit, snapshotState,
} from "./state.mjs";
import { Tooltips } from "./tooltips.mjs";
import { drawWaveform } from "./waveform_view.mjs";

const byId = id => document.getElementById(id);
new Tooltips();
const state = {
  reference: null, references: [], synthesis: null, snapshots: [],
  macros: [], patch: null, recipeIndex: 0, recipeKey: "metal.cymbal.v1",
  activeSnapshotId: null,
  event: {
    strength: 0.8, location: 0.8, hardness: 0.65, implement: 1,
    contactSpread: 0.2, constraint: 0, seed: 17,
  },
  eventDefaults: {
    strength: 0.8, location: 0.8, hardness: 0.65, implement: 1,
    contactSpread: 0.2, constraint: 0,
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
let midiTriggerSeed = 0;
const settings = new SettingsController({
  audition,
  onStatus: setStatus,
  onMidiNote: ({ velocity }) => {
    const event = {
      ...state.event,
      strength: Math.max(.01, velocity),
      seed: (state.event.seed + ++midiTriggerSeed) >>> 0,
    };
    audition.trigger(event).catch(error => setStatus(String(error)));
  },
});
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
let performanceControls;
let recipeController;
let routingController;
let restoringReference = false;
let levelMatchWarning = "";

function setStatus(message) { byId("status").textContent = message; }
function setReadyIfIdle() {
  if (!pendingAnalysis.size && !renderInFlight && !renderTimer &&
      !pendingLevelMatch) {
    setStatus(levelMatchWarning ? `Ready · ${levelMatchWarning}` : "Ready");
  }
}
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
  levelMatchWarning = "";
  pendingLevelMatch = { referenceId: state.reference.id };
  refreshModelLevelControl();
  audition.setMacros(state.macros);
  setStatus("Matching model level…");
  renderSynthesis();
}

function buildControls(resetValues = true) {
  if (resetValues)
    state.macros = engine.parameters.map(item => item.defaultValue);
  const ControlType = state.recipeKey === "metal.cymbal.v1"
    ? FitControls : state.recipeKey === "drum.membrane.v1"
      ? MembraneControls : state.recipeKey === "drum.snare.v1"
        ? SnareControls : KickControls;
  fitControls = new ControlType({
    descriptors: engine.parameters, state,
    onChange: key => {
      levelMatchWarning = "";
      if (key === "model_level_db") {
        state.modelLevelTouched = true;
        pendingLevelMatch = null;
      }
      scheduleRender();
    },
    onLevelReset: () => requestLevelMatch(),
    onPreset: key => applyMembranePreset(key),
  });
  fitControls.build();
}

function applyMembranePreset(key) {
  const values = membranePresetValues(key, engine.parameters);
  state.macros.splice(0, state.macros.length, ...values);
  state.patch = key === "acousticKick"
    ? createAcousticKickPatch(engine.parameters)
    : createTomPatch(engine.parameters, values);
  state.modelLevelTouched = true;
  pendingLevelMatch = null;
  routingController?.setPatch(state.patch);
  routingController?.refreshPresentation();
  buildPageValues();
  audition.setRecipe(
    state.recipeIndex, state.macros,
    recipeAdapter(state.recipeKey).routing(state.patch),
  );
  scheduleRender(false);
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
  setReadyIfIdle();
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
  setReadyIfIdle();
};

function scheduleRender(updateLive = true) {
  clearTimeout(renderTimer);
  if (updateLive) audition.setMacros(state.macros);
  setStatus("Rendering…");
  renderTimer = setTimeout(() => {
    renderTimer = undefined;
    renderSynthesis();
  }, 220);
}

function renderSynthesis() {
  clearTimeout(renderTimer);
  renderTimer = undefined;
  const sampleRate = state.reference?.sampleRate ?? 48000;
  const duration = byId("render-seconds").value;
  const seconds = duration === "reference"
    ? (state.reference?.duration ?? 6) : Number(duration);
  const request = {
    generation: ++renderGeneration, recipeIndex: state.recipeIndex,
    sampleRate, seconds, parameters: [...state.macros],
    routing: recipeAdapter(state.recipeKey).routing(state.patch),
    event: { ...state.event },
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
        const match = modelLevelMatch({
          currentDb: state.macros[0],
          reference: state.reference.samples,
          referenceSampleRate: state.reference.sampleRate,
          synthesis: data.samples,
          synthesisSampleRate: data.sampleRate,
          minimumDb: descriptor.minimum,
          maximumDb: descriptor.maximum,
          referenceStartSeconds: state.reference.cell?.onset_seconds ?? 0,
        });
        const matched = match.appliedDb;
        pendingLevelMatch = null;
        if (match.clipped) {
          const shortfall = Math.max(0, match.requestedDb - match.appliedDb);
          levelMatchWarning = shortfall > .05
            ? `Model is ${shortfall.toFixed(1)} dB too quiet at its 0 dB ceiling — preset is not calibrated`
            : "Reference level is outside the attenuation-only model range";
        }
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
      drawWaveform(state);
      byId("render-time").textContent = `${data.elapsedMs.toFixed(0)} ms DSP`;
    }
  }
  if (queuedRender) {
    const request = queuedRender;
    queuedRender = null;
    dispatchRender(request);
  } else setReadyIfIdle();
};

renderWorker.onerror = event => {
  renderInFlight = false;
  pendingLevelMatch = null;
  setStatus(event.message);
};

function setReference(reference) {
  const firstReference = !state.reference;
  state.reference = reference;
  fitControls?.decayEditor?.refresh();
  audition.setTrim(reference.corpus?.auditionTrimDb ?? 0);
  byId("calibration-trim").textContent =
    `Corpus monitor gain ${audition.trimDb >= 0 ? "+" : ""}${audition.trimDb.toFixed(1)} dB`;
  if (reference.cell) {
    Object.assign(state.eventDefaults, {
      strength: reference.cell.strength,
      location: reference.cell.location,
      hardness: reference.cell.hardness,
      implement: reference.cell.implement ?? 1,
      contactSpread: reference.cell.contactSpread ?? .2,
      constraint: reference.cell.constraint ?? 0,
    });
    Object.assign(state.event, {
      ...state.eventDefaults,
      seed: reference.cell.seed,
    });
    performanceControls.paint();
  }
  if (firstReference) {
    byId("view-mode").value = "mirror";
    view.setSettings({ mode: "mirror" });
  }
  analyzeReference();
  const referenceWindow = alignedReferenceWindow(
    reference.duration, reference.cell?.onset_seconds ?? 0,
  );
  drawWaveform(state);
  view.reset(referenceWindow.duration, referenceWindow.offset);
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
  const recipe = engine.recipes.find(entry =>
    entry.key === fit.instrument.recipe);
  if (!recipe) throw new Error(`recipe is unavailable: ${fit.instrument.recipe}`);
  if (recipe.index !== state.recipeIndex) {
    recipeController.remember();
    recipeController.activate(recipe.index);
  }
  state.patch = structuredClone(fit.instrument);
  state.macros.splice(0, state.macros.length,
    ...fitMacroValues(fit, engine.parameters));
  state.modelLevelTouched = true;
  pendingLevelMatch = null;
  Object.assign(state.event, fit.controls.event);
  Object.assign(state.analysis, fit.controls.analysis);
  state.activeSnapshotId = fit.id;
  routingController.setPatch(state.patch);
  if (item.audio) {
    state.synthesis = item.audio.slice();
    analyze(
      "synthesis", state.synthesis, fit.reference?.sampleRate ??
        state.reference?.sampleRate ?? 48000,
    );
  }
  buildPageValues();
  audition.setRecipe(
    state.recipeIndex, state.macros,
    recipeAdapter(state.recipeKey).routing(state.patch),
  );
  drawWaveform(state);
  renderSnapshotList();
  if (!item.audio) scheduleRender();
}

function buildPageValues() {
  fitControls.build();
  performanceControls.paint();
  byId("colour-range").value = state.analysis.dynamicRangeDb;
  byId("colour-range").nextElementSibling.textContent =
    `${state.analysis.dynamicRangeDb} dB`;
  view.setSettings({ dynamicRangeDb: state.analysis.dynamicRangeDb });
}

async function initialize() {
  engine = await PercussionEngine.create();
  recipeController = new RecipeController({
    engine, state, audition,
    getRoutingController: () => routingController,
    buildControls, buildPageValues,
    onChanged: () => {
      // Reference targets carry a collection-level monitoring calibration. Keep it
      // fixed when neighbouring velocity cells are selected; per-cell matching
      // would erase the source velocity curve.
      state.modelLevelTouched = true;
      pendingLevelMatch = null;
      renderSynthesis();
    },
  });
  recipeController.populate();
  buildControls();
  state.patch = recipeAdapter(state.recipeKey).create(
    engine.parameters, state.macros,
  );
  routingController = new RoutingController({
    state, engine, audition, scheduleRender, setStatus,
  });
  routingController.bind();
  performanceControls = new PerformanceControls({
    state, audition, scheduleRender, setStatus,
  });
  performanceControls.bind();
  audition.setRecipe(
    state.recipeIndex, state.macros,
    recipeAdapter(state.recipeKey).routing(state.patch),
  );
  settings.bind();
  audition.initialize().catch(error => setStatus(String(error)));
  referenceBrowser = new ReferenceBrowser({
    corpus: byId("reference-corpus"),
    articulation: byId("reference-articulation"),
    velocity: byId("reference-velocity"),
    repeat: byId("reference-repeat"),
  }, setReference, setStatus);
  await referenceBrowser.initialize();
  bindCalibrationPresets();
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
  bindAnalysisControls({
    state, view, analyzeReference,
    analyzeSynthesis: () => analyze(
      "synthesis", state.synthesis, state.reference?.sampleRate ?? 48000),
    scheduleRender,
  });
  bindSnapshotControls();
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
        ? `Live DSP ready · ${audition.macroCommitMs.toFixed(0)} ms` +
          ` · install ${audition.installationMs.toFixed(2)} ms`
        : "Live DSP idle";
    settings.paintAudioStatus();
  }, 100);
}

function bindCalibrationPresets() {
  const select = byId("instrument-calibration");
  const calibrations = referenceBrowser.calibrationPresets();
  select.replaceChildren(
    new Option("Choose a reference target…", ""),
    ...calibrations.map(item => new Option(item.name, item.id)),
  );
  select.onchange = async () => {
    const calibration = calibrations.find(item => item.id === select.value);
    if (!calibration) return;
    try {
      const recipe = engine.recipes.find(
        item => item.key === calibration.recipe);
      if (!recipe) throw new Error(`recipe is unavailable: ${calibration.recipe}`);
      recipeController.remember();
      recipeController.activate(recipe.index);
      state.macros.splice(0, state.macros.length,
        ...calibrationParameterValues(
          calibration, engine.parameters, { strict: true },
        ));
      state.patch = calibrationPatch(
        calibration, engine.parameters, state.macros, state.patch,
      );
      state.patch = recipeAdapter(state.recipeKey).withValues(
        state.patch, engine.parameters, state.macros,
      );
      // A named reference target owns one fixed collection-level trim. Loading
      // its reference must not silently replace that value with a per-cell
      // match.
      state.modelLevelTouched = true;
      routingController.setPatch(state.patch);
      routingController.refreshPresentation();
      buildPageValues();
      audition.setRecipe(
        state.recipeIndex, state.macros,
        recipeAdapter(state.recipeKey).routing(state.patch),
      );
      const loaded = await referenceBrowser.selectSavedReference({
        corpus: { id: calibration.corpusId }, cell: calibration,
      });
      if (!loaded) throw new Error(`reference is unavailable: ${calibration.name}`);
    } catch (error) {
      setStatus(String(error));
    }
  };
}

function bindSnapshotControls() {
  byId("snapshot").onclick = () => {
    const fit = snapshotState(
      state, byId("snapshot-name").value || "Candidate", engine.macros,
    );
    state.snapshots.push({ fit, audio: state.synthesis?.slice() });
    state.activeSnapshotId = fit.id; renderSnapshotList();
  };
  byId("save-fit").onclick = () => {
    const active = state.snapshots.find(item => item.fit.id === state.activeSnapshotId);
    downloadFit(active?.fit ?? snapshotState(
      state, byId("snapshot-name").value, engine.macros,
    ));
  };
  byId("load-fit").onchange = async event => {
    try {
      const fit = await readFit(event.target.files[0], recipeKey => {
        const recipe = engine.recipes.find(item => item.key === recipeKey);
        if (!recipe) throw new Error(`recipe is unavailable: ${recipeKey}`);
        return engine.descriptorsForRecipe(recipe.index);
      });
      const item = { fit }; state.snapshots.push(item); await restore(item);
    } catch (error) { setStatus(String(error)); }
  };
}

initialize().catch(error => setStatus(String(error)));
