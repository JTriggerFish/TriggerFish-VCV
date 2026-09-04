// Audit the exact reference starts exposed by the browser workbench.
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { webcrypto } from "node:crypto";

globalThis.crypto ??= webcrypto;

const site = resolve(process.argv[2] ?? "build/workbench-wasm/site");
const baseUrl = new URL(process.argv[3] ?? "http://127.0.0.1:8765");
const importFromSite = name => import(pathToFileURL(resolve(site, name)));
const { PercussionEngine } = await importFromSite("engine.mjs");
const { decodeWav } = await importFromSite("references.mjs");
const { centeredErrorStatistics, stft } = await importFromSite("analysis.mjs");
const { calibrationParameterValues, calibrationPatch } =
  await importFromSite("instrument_calibrations.mjs");
const { recipeAdapter } = await importFromSite("recipe_adapter.mjs");

const selectedId = process.env.TF_AUDIT_ID ?? "";
const parameterOverrides = JSON.parse(
  process.env.TF_AUDIT_OVERRIDES ?? "{}",
);
const routingOverrides = JSON.parse(
  process.env.TF_AUDIT_ROUTES ?? "{}",
);

const Regions = [
  ["contact", 0, .015], ["presence", .015, .12],
  ["body", .12, .6], ["tail", .6, 1.5],
];
const BandEdges = [20, 150, 400, 1000, 3000, 8000, 16000, 22000];

function energy(samples) {
  let result = 0;
  for (const sample of samples) result += sample * sample;
  return result;
}

function peak(samples) {
  let result = 0;
  for (const sample of samples) result = Math.max(result, Math.abs(sample));
  return result;
}

function ratioDb(numerator, denominator) {
  return 10 * Math.log10(Math.max(numerator, 1e-30) /
    Math.max(denominator, 1e-30));
}

function bandProfile(samples, sampleRate) {
  const spectrum = stft(samples, sampleRate, {
    size: 2048, hop: 256, window: "hann", floorDb: -180,
  });
  const result = [];
  for (const [region, start, end] of Regions) {
    const firstFrame = Math.max(0, Math.ceil(start * sampleRate / spectrum.hop));
    const lastFrame = Math.min(
      spectrum.frames, Math.ceil(end * sampleRate / spectrum.hop),
    );
    for (let band = 0; band + 1 < BandEdges.length; ++band) {
      const low = BandEdges[band];
      const high = Math.min(BandEdges[band + 1], sampleRate / 2);
      if (high <= low || lastFrame <= firstFrame) continue;
      const firstBin = Math.max(0, Math.ceil(low * spectrum.size / sampleRate));
      const lastBin = Math.min(
        spectrum.bins, Math.ceil(high * spectrum.size / sampleRate),
      );
      let power = 0;
      let count = 0;
      for (let frame = firstFrame; frame < lastFrame; ++frame) {
        for (let bin = firstBin; bin < lastBin; ++bin) {
          power += 10 ** (spectrum.values[frame * spectrum.bins + bin] / 10);
          ++count;
        }
      }
      result.push({
        region, band: `${low}-${high}`,
        db: 10 * Math.log10(Math.max(power / Math.max(1, count), 1e-30)),
      });
    }
  }
  return result;
}

function compareProfiles(reference, synthesis) {
  const referencePeakDb = Math.max(...reference.map(item => item.db));
  const errors = reference.map((item, index) => ({
    region: item.region, band: item.band,
    referenceDb: item.db, synthesisDb: synthesis[index].db,
    errorDb: synthesis[index].db - item.db,
  }));
  // Do not let reference-noise-floor cells dominate the score. They remain in
  // `errors` for inspection but do not participate in summary statistics.
  const active = errors.filter(item => item.referenceDb >= referencePeakDb - 80);
  const summary = centeredErrorStatistics(active.map(item => item.errorDb));
  const regional = Object.fromEntries(Regions.map(([name]) => {
    const selected = active.filter(item => item.region === name);
    if (!selected.length) return [name, null];
    const statistics = centeredErrorStatistics(
      selected.map(item => item.errorDb),
    );
    return [name, {
      biasDb: statistics.mean,
      shapeRmseDb: statistics.rmse,
    }];
  }));
  return {
    activeSpectralCells: active.length,
    spectralBiasDb: summary.mean, shapeRmseDb: summary.rmse, regional, errors,
  };
}

async function loadReference(corpus, cell) {
  const response = await fetch(new URL(cell.url, baseUrl));
  if (!response.ok) throw new Error(`could not load ${cell.label}`);
  return decodeWav(await response.arrayBuffer(), cell.label);
}

async function audit(corpus) {
  const calibration = corpus.calibration;
  const cell = corpus.cells.find(item =>
    item.articulation === calibration.articulation &&
    item.velocity === calibration.velocity &&
    item.repeat === calibration.repeat);
  if (!cell) throw new Error(`missing reference cell for ${calibration.id}`);
  const reference = await loadReference(corpus, cell);
  const engine = await PercussionEngine.create(reference.sampleRate);
  try {
    const recipe = engine.recipes.find(item => item.key === calibration.recipe);
    if (!recipe) throw new Error(`missing recipe ${calibration.recipe}`);
    engine.setRecipe(recipe.index);
    const parameters = calibrationParameterValues(
      calibration, engine.parameters, { strict: true },
    );
    for (const [key, value] of Object.entries(parameterOverrides)) {
      const descriptor = engine.parameters.find(item => item.key === key);
      if (!descriptor) throw new Error(`unknown parameter override: ${key}`);
      parameters[descriptor.index] = Number(value);
    }
    let patch = recipeAdapter(recipe.key).create(engine.parameters, parameters);
    patch = calibrationPatch(
      calibration, engine.parameters, parameters, patch,
    );
    patch = recipeAdapter(recipe.key).withValues(
      patch, engine.parameters, parameters,
    );
    for (const [id, enabled] of Object.entries(routingOverrides)) {
      const route = patch.connections.find(item => item.id === id);
      if (!route) throw new Error(`unknown routing override: ${id}`);
      if (typeof enabled !== "boolean")
        throw new Error(`routing override must be boolean: ${id}`);
      route.enabled = enabled;
    }
    const onset = Math.max(0, Math.round(
      (cell.onset_seconds ?? 0) * reference.sampleRate,
    ));
    const alignedReference = reference.samples.slice(onset);
    const seconds = alignedReference.length / reference.sampleRate;
    const synthesis = engine.render({
      seconds, parameters,
      routing: recipeAdapter(recipe.key).routing(patch),
      strength: cell.strength, location: cell.location,
      hardness: cell.hardness, implement: cell.implement,
      contactSpread: cell.contactSpread, constraint: cell.constraint ?? 0,
      seed: cell.seed,
    });
    const levelErrorDb = ratioDb(energy(synthesis), energy(alignedReference));
    return {
      id: calibration.id, name: calibration.name,
      recipe: recipe.key, reference: cell.label,
      durationSeconds: seconds,
      modelLevelDb: parameters[0],
      levelErrorDb,
      peakErrorDb: 20 * Math.log10(
        Math.max(peak(synthesis), 1e-30) /
        Math.max(peak(alignedReference), 1e-30),
      ),
      ...compareProfiles(
        bandProfile(alignedReference, reference.sampleRate),
        bandProfile(synthesis, reference.sampleRate),
      ),
      acceptedByListening: false,
    };
  } finally {
    engine.destroy();
  }
}

const response = await fetch(new URL("/api/reference-corpora", baseUrl));
if (!response.ok) throw new Error("could not load reference catalog");
const corpora = (await response.json()).corpora.filter(item =>
  item.calibration && (!selectedId || item.calibration.id === selectedId));
if (selectedId && corpora.length !== 1)
  throw new Error(`unknown calibration id: ${selectedId}`);
const starts = [];
for (const corpus of corpora) starts.push(await audit(corpus));
console.log(JSON.stringify({
  note: "Diagnostics only. A numerical result never promotes a listening start.",
  bandEdgesHz: BandEdges, regions: Regions, starts,
}, null, 2));
