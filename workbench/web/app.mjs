import { SafeAudition } from "./audio.mjs";
import { paintLimiterMeter } from "./limiter_meter.mjs";
import { reportError } from "./error_banner.mjs";
import { mountDiagnosticAudition } from "./diagnostic_audition.mjs";
import { bindAnalysisControls } from "./analysis_controls.mjs";
import { PercussionEngine } from "./engine.mjs";
import { FitControls } from "./fit_controls.mjs";
import {
  calibrationParameterValues, calibrationPatch, calibrationEvent, calibrationDisplayName,
} from "./instrument_calibrations.mjs";
import { KickControls } from "./kick_controls.mjs";
import { MembraneControls } from "./membrane_controls.mjs";
import {
  createTomPatch, membranePresetValues,
} from "./membrane_patch.mjs";
import { SnareControls } from "./snare_controls.mjs";
import { PerformanceControls } from "./performance_controls.mjs";
import { recipeAdapter } from "./recipe_adapter.mjs";
import { RecipeController } from "./recipe_controller.mjs";
import { readReferences, setReferenceGain } from "./references.mjs";
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
import { mountTextureTrials } from "./texture_trials.mjs";

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
mountDiagnosticAudition(audition, setStatus);
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
let fitControls;
let referenceBrowser;
let performanceControls;
let recipeController;
let routingController;
let restoringReference = false;

function setStatus(message) {
  byId("status").textContent = String(message);
  if (message instanceof Error || /^(?:\w*Error:)/.test(String(message))) reportError(message);
}
function setReadyIfIdle() {
  if (!pendingAnalysis.size && !renderInFlight && !renderTimer) {
    setStatus("Ready");
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

function buildControls(resetValues = true) {
  fitControls?.destroy?.();
  if (resetValues)
    state.macros = engine.parameters.map(item => item.defaultValue);
  const ControlType = state.recipeKey === "metal.cymbal.v1"
    ? FitControls : state.recipeKey === "drum.membrane.v1"
      ? MembraneControls : state.recipeKey === "drum.snare.v1"
        ? SnareControls : KickControls;
  fitControls = new ControlType({
    descriptors: engine.parameters, state,
    decayHold: {
      read: () => ({parameters:state.macros.slice(), event:{...state.event},
        sampleRate:state.reference?.sampleRate ?? 48000, recipeIndex:state.recipeIndex,
        reference:state.reference?.id,
        routing:recipeAdapter(state.recipeKey).routing(state.patch)}),
      onError: reportError,
    },
    onChange: () => scheduleRender(),
    onLevelReset: () => {
      state.macros[0] = engine.macros[0].defaultValue;
      refreshModelLevelControl();
      scheduleRender();
    },
    onPreset: key => applyMembranePreset(key),
  });
  fitControls.build();
}

function applyMembranePreset(key) {
  const values = membranePresetValues(key, engine.parameters);
  state.macros.splice(0, state.macros.length, ...values);
  state.patch = createTomPatch(engine.parameters, values);
  routingController?.setPatch(state.patch);
  routingController?.refreshPresentation();
  buildPageValues();
  audition.setRecipe(
    state.recipeIndex, state.macros,
    recipeAdapter(state.recipeKey).routing(state.patch),
  );
  scheduleRender(false);
}

function analyze(kind, samples, sampleRate, cacheKey, preview = false, renderId) {
  const generation = ++generations[kind];
  pendingAnalysis.add(kind);
  setStatus(`Analyzing ${[...pendingAnalysis].join(" + ")}…`);
  const copy = samples.slice();
  worker.postMessage({
    generation, kind, samples: copy, sampleRate,
    settings: state.analysis, cacheKey, preview, renderId,
  }, [copy.buffer]);
}

function referenceSpectrumKey(reference) {
  const source = reference.sha256 ?? reference.id;
  const { size, hop, window, floorDb } = state.analysis;
  return `${source}|${reference.referenceGainDb ?? 0}|${size}|${hop}|${window}|${floorDb}`;
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
  fitControls?.refreshRadiation?.();
  view.setData("reference", cached);
  updateColourCeiling();
  setReadyIfIdle();
}

worker.onmessage = ({ data }) => {
  if (data.generation !== generations[data.kind]) return;
  pendingAnalysis.delete(data.kind);
  if (data.error) { setStatus(new Error(data.error)); return; }
  view.setData(data.kind, data.result);
  state[`${data.kind}Spectrum`] = data.result;
  fitControls?.refreshRadiation?.();
  if (data.kind === "reference" && data.cacheKey) {
    cacheReferenceSpectrum(data.cacheKey, data.result);
    updateColourCeiling();
  }
  setReadyIfIdle();
};

function invalidateRender() {
  clearTimeout(renderTimer);
  renderTimer = undefined;
  ++renderGeneration;
  ++generations.synthesis;
  pendingAnalysis.delete("synthesis");
  renderWorker.postMessage({cancel: true});
  renderInFlight = false;
  state.synthesisCurrent = false;
}

function scheduleRender(updateLive = true) {
  invalidateRender();
  if (updateLive) audition.setMacros(state.macros);
  setStatus("Rendering…");
  renderTimer = setTimeout(() => {
    renderTimer = undefined;
    renderSynthesis();
  }, 60);
}

function renderSynthesis() {
  clearTimeout(renderTimer);
  renderTimer = undefined;
  ++generations.synthesis;
  pendingAnalysis.delete("synthesis");
  state.synthesisCurrent = false;
  const sampleRate = state.reference?.sampleRate ?? 48000;
  const duration = byId("render-seconds").value;
  const seconds = duration === "reference"
    ? (state.reference?.duration ?? 6) : Number(duration);
  const request = {
    generation: ++renderGeneration, recipeIndex: state.recipeIndex,
    sampleRate, seconds, parameters: [...state.macros],
    previewFrames: Math.max(Math.ceil(.125 * sampleRate), state.analysis.size / 2 + state.analysis.hop),
    routing: recipeAdapter(state.recipeKey).routing(state.patch),
    event: { ...state.event },
  };
  dispatchRender(request);
}

function dispatchRender(request) {
  renderInFlight = true;
  renderWorker.postMessage(request);
}

renderWorker.onmessage = ({ data }) => {
  if (data.generation !== renderGeneration) return;
  renderInFlight = Boolean(data.preview);
  if (data.error) {
    setStatus(new Error(data.error));
    return;
  }
  analyze("synthesis", data.samples, data.sampleRate, undefined, data.preview, data.generation);
  if (!data.preview) {
    state.synthesis = data.samples;
    state.synthesisCurrent = true;
    drawWaveform(state);
  }
  byId("render-time").textContent = data.preview
    ? `Rendering ${(data.samples.length / data.sampleRate).toFixed(1)} / ${data.seconds.toFixed(1)} s…`
    : `${data.elapsedMs.toFixed(0)} ms DSP · ${(1000 * data.seconds / data.elapsedMs).toFixed(1)}× real time`;
  setReadyIfIdle();
};

function renderWorkerFailed(message) {
  renderInFlight = false;
  setStatus(new Error(message));
}
function spectrumWorkerFailed(message) {
  pendingAnalysis.clear();
  setStatus(new Error(message));
}
renderWorker.onerror = event => renderWorkerFailed(event.message);
worker.onerror = event => spectrumWorkerFailed(`Spectrogram worker: ${event.message}`);
worker.addEventListener("messageerror", () =>
  spectrumWorkerFailed("Could not decode the spectrogram worker response"));
renderWorker.addEventListener("messageerror", () =>
  renderWorkerFailed("Could not decode the render worker response"));

function setReference(reference) {
  const firstReference = !state.reference;
  state.reference = reference;
  fitControls?.decayEditor?.refresh();
  paintReferenceGain();
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
    // Selecting another layer must never alter the saved instrument gain.
    scheduleRender();
  }
}

function paintReferenceGain() {
  const input = byId("reference-gain");
  input.disabled = !state.reference;
  input.value = state.reference?.referenceGainDb ?? 0;
  input.nextElementSibling.textContent = `${Number(input.value).toFixed(1)} dB`;
}

function editReferenceGain(gainDb) {
  if (!state.reference) return;
  state.reference = setReferenceGain(state.reference, gainDb);
  referenceBrowser?.rememberGain(state.reference);
  paintReferenceGain();
  analyzeReference();
  drawWaveform(state);
}

byId("reference-gain").oninput = event => editReferenceGain(Number(event.target.value));
byId("reference-gain").ondblclick = event => {
  event.preventDefault();
  editReferenceGain(state.reference?.corpus?.referenceGainDb ?? 0);
};

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
  // A cached snapshot must also cancel any older in-flight preview/render.
  invalidateRender();
  const fit = item.fit;
  const sameSource = !fit.reference ||
    (fit.reference.sha256 && fit.reference.sha256 === state.reference?.sha256) ||
    fit.reference.id === state.reference?.id;
  if (sameSource && fit.reference && state.reference) {
    state.reference = setReferenceGain(state.reference, fit.reference.referenceGainDb ?? 0);
    referenceBrowser?.rememberGain(state.reference);
    paintReferenceGain();
  }
  const sameReference = sameSource;
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
  Object.assign(state.event, fit.controls.event);
  Object.assign(state.analysis, fit.controls.analysis);
  analyzeReference();
  state.activeSnapshotId = fit.id;
  routingController.setPatch(state.patch);
  if (item.audio) {
    state.synthesis = item.audio.slice();
    state.synthesisCurrent = true;
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
  mountTextureTrials(byId("texture-trial"), async item => {
    await restore(item);
    byId("instrument-calibration").value = "gong-standard";
  }, setStatus);
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
    if(document.hidden || !fitControls?.eqEditors?.length)return;
    const live=audition.readOutputSpectrum();
    const changed=Boolean(live)||Boolean(state.liveEqHistogram);
    state.liveEqHistogram=live;
    state.liveEqSampleRate=live ? audition.sampleRate : null;
    if(changed)fitControls.eqEditors.forEach(editor=>editor.background());
  }, 50);
  byId("limiter-reset").onclick = () => {
    audition.limiterMeter.reset();
    paintLimiterMeter(byId("limiter-meter"), audition.limiterMeter.read());
  };
  setInterval(() => {
    paintLimiterMeter(byId("limiter-meter"), audition.limiterMeter.read());
    const latency = audition.latencyMs ? ` · ${audition.latencyMs.toFixed(0)} ms` : "";
    const underflows = audition.underflows ? ` · xruns ${audition.underflows}` : "";
    const output = Number.isFinite(audition.outputDb)
      ? ` · out ${audition.outputDb.toFixed(0)} dBFS` : "";
    const live = audition.state === "off" ? "" :
      ` · ${(audition.sampleRate / 1000).toFixed(1)} kHz` +
      ` · ${audition.state} · hits ${audition.triggerCount}`;
    const input = Number.isFinite(audition.inputPeakDb)
      ? ` · limiter in ${audition.inputPeakDb.toFixed(1)} dBFS` : "";
    byId("limiter").textContent =
      `Audio${input}${output}${latency}${live}${underflows}`;
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
  let selectionGeneration = 0;
  const select = byId("instrument-calibration");
  const calibrations = referenceBrowser.calibrationPresets();
  select.replaceChildren(
    new Option("Choose a reference target…", ""),
    ...calibrations.map(item => new Option(calibrationDisplayName(item), item.id)),
  );
  select.onchange = async () => {
    const generation = ++selectionGeneration;
    const calibration = calibrations.find(item => item.id === select.value);
    if (!calibration) return;
    byId("texture-trial").value = "";
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
      // Loading a target never computes a reference/model level correction.
      routingController.setPatch(state.patch);
      routingController.refreshPresentation();
      buildPageValues();
      audition.setRecipe(
        state.recipeIndex, state.macros,
        recipeAdapter(state.recipeKey).routing(state.patch),
      );
      const loading = referenceBrowser.selectSavedReference({
        corpus: { id: calibration.corpusId }, cell: calibration,
      });
      const referenceGeneration = referenceBrowser.generation;
      const loaded = await loading;
      if (generation !== selectionGeneration ||
          referenceGeneration !== referenceBrowser.generation ||
          state.recipeKey !== calibration.recipe) return;
      if (!loaded) throw new Error(`reference is unavailable: ${calibration.name}`);
      const fittedEvent = calibrationEvent(calibration);
      if (fittedEvent) {
        Object.assign(state.event, fittedEvent);
        performanceControls.paint();
        scheduleRender();
      }
      byId("snapshot-name").value = calibrationDisplayName(calibration);
    } catch (error) {
      if (generation === selectionGeneration) setStatus(error);
    }
  };
}

function bindSnapshotControls() {
  byId("snapshot").onclick = () => {
    const fit = snapshotState(
      state, byId("snapshot-name").value || "Candidate", engine.macros,
    );
    state.snapshots.push({ fit, audio: state.synthesisCurrent ? state.synthesis?.slice() : undefined });
    state.activeSnapshotId = fit.id; renderSnapshotList();
  };
  byId("save-fit").onclick = () => {
    // Save the visible controls, including edits made after selecting a snapshot.
    downloadFit(snapshotState(
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
