// Persistent developer-only renderer. JSON lines in; float32 PCM out as base64.
// Uses the very same parameter adapter, source decoding and Wasm as the UI.
import { createInterface } from "node:readline";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { webcrypto, createHash } from "node:crypto";
import { readFile } from "node:fs/promises";

globalThis.crypto ??= webcrypto;
const site = resolve(process.argv[2] ?? "build/workbench-wasm/site");
const base = process.argv[3] ?? "http://127.0.0.1:8765";
const load = name => import(pathToFileURL(resolve(site, name)));
const { PercussionEngine } = await load("engine.mjs");
const { readRemoteReference } = await load("references.mjs");
const { calibrationParameterValues, calibrationPatch } = await load("instrument_calibrations.mjs");
const { recipeAdapter } = await load("recipe_adapter.mjs");
const { snapshotState, fitMacroValues, validateFit } = await load("state.mjs");
const pcm = samples => Buffer.from(samples.buffer, samples.byteOffset, samples.byteLength).toString("base64");
let engine, reference, event, defaults, patch;

async function initialize(id) {
  engine?.destroy();
  const catalog = await (await fetch(`${base}/api/reference-corpora`)).json();
  const corpus = catalog.corpora.find(item => item.calibration?.id === id);
  if (!corpus) throw new Error(`Unknown reference start: ${id}`);
  const calibration = corpus.calibration;
  const cell = corpus.cells.find(item => ["articulation", "velocity", "repeat"]
    .every(key => item[key] === calibration[key]));
  reference = await readRemoteReference(corpus, { ...cell, url: new URL(cell.url, base).href });
  engine = await PercussionEngine.create(reference.sampleRate);
  engine.setRecipe(engine.recipes.find(item => item.key === calibration.recipe).index);
  defaults = calibrationParameterValues(calibration, engine.parameters);
  const adapter = recipeAdapter(calibration.recipe);
  patch = calibrationPatch(calibration, engine.parameters, defaults,
    adapter.create(engine.parameters, defaults));
  event = Object.fromEntries(["strength", "location", "hardness", "implement",
    "contactSpread", "seed"].map(key => [key, cell[key]]));
  event.constraint = cell.constraint ?? 0;
  return {
    recipeKey: calibration.recipe,
    descriptors: engine.parameters, values: defaults, event,
    rendererSha256: createHash("sha256").update(await readFile(resolve(site, "triggerfish-percussion.wasm"))).digest("hex"),
    reference: { ...reference, samples: undefined, rawSamples: undefined },
    pcm: pcm(reference.samples),
  };
}

function parameters(overrides = {}) {
  const result = defaults.slice();
  for (const [key, value] of Object.entries(overrides)) {
    const descriptor = engine.parameters.find(item => item.key === key);
    if (!descriptor || !Number.isFinite(value) || value < descriptor.minimum || value > descriptor.maximum)
      throw new Error(`Invalid workbench parameter: ${key}=${value}`);
    result[descriptor.index] = value;
  }
  return result;
}

function render(request) {
  const values = parameters(request.parameters);
  const adapter = recipeAdapter(patch.recipe);
  engine.reset();
  const started = performance.now();
  const samples = engine.render({ seconds: request.seconds ?? reference.duration,
    ...event, seed: request.seed ?? event.seed, parameters: values,
    routing: adapter.routing(patch) });
  return { pcm: pcm(samples), elapsedMs: performance.now() - started };
}

function snapshot(request) {
  const macros = parameters(request.parameters);
  const result = snapshotState({
    reference, event, macros, patch,
    analysis: { size: 4096, hop: 512, window: "hann", floorDb: -160, dynamicRangeDb: 80 },
  }, request.name ?? "Fitting candidate — not listening-approved", engine.parameters);
  return { fit: result };
}

function renderSequence(request) {
  const seconds = request.seconds ?? 3;
  const hits = request.hits ?? [];
  if (!(seconds > 0 && seconds <= 30) || hits.length > 256)
    throw new Error("Invalid sequence duration or hit count");
  engine.reset();
  engine.setConfiguration(parameters(request.parameters), recipeAdapter(patch.recipe).routing(patch));
  const samples = new Float32Array(Math.round(seconds * engine.sampleRate));
  let cursor = 0;
  for (const hit of hits) {
    const frame = Math.round(hit.time * engine.sampleRate);
    if (!Number.isFinite(frame) || frame < cursor || frame >= samples.length)
      throw new Error("Sequence hits must be ordered and within the render");
    if (frame > cursor) engine.processTo(samples, cursor, frame - cursor);
    engine.trigger({ ...event, ...hit });
    cursor = frame;
  }
  if (cursor < samples.length) engine.processTo(samples, cursor, samples.length - cursor);
  return { pcm: pcm(samples) };
}

function renderSnapshot(request) {
  const fit = validateFit(request.fit, engine.parameters);
  engine.reset();
  const samples = engine.render({
    seconds: request.seconds ?? reference.duration, ...fit.controls.event,
    parameters: fitMacroValues(fit, engine.parameters),
    routing: recipeAdapter(fit.instrument.recipe).routing(fit.instrument),
  });
  return { pcm: pcm(samples) };
}

for await (const line of createInterface({ input: process.stdin })) {
  try {
    const request = JSON.parse(line);
    const result = request.command === "initialize" ? await initialize(request.id)
      : request.command === "snapshot" ? snapshot(request)
        : request.command === "renderSnapshot" ? renderSnapshot(request)
          : request.command === "renderSequence" ? renderSequence(request) : render(request);
    process.stdout.write(JSON.stringify(result) + "\n");
  } catch (error) { process.stdout.write(JSON.stringify({ error: String(error) }) + "\n"); }
}
engine?.destroy();
