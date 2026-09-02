import { readFile } from "node:fs/promises";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { webcrypto } from "node:crypto";

globalThis.crypto ??= webcrypto;

const site = resolve(process.argv[2] ?? "build/workbench-wasm/site");
const referencePath = resolve(process.argv[3]);
const { CrashEngine } = await import(pathToFileURL(resolve(site, "engine.mjs")));
const { decodeWav } = await import(pathToFileURL(resolve(site, "references.mjs")));
const { stft } = await import(pathToFileURL(resolve(site, "analysis.mjs")));

const file = await readFile(referencePath);
const buffer = file.buffer.slice(file.byteOffset, file.byteOffset + file.byteLength);
const reference = await decodeWav(buffer, referencePath);
const engine = await CrashEngine.create(reference.sampleRate);
const descriptors = new Map(engine.macros.map(item => [item.key, item]));
const defaults = engine.macros.map(item => item.defaultValue);
const event = {
  strength: 96 / 127, location: 1, hardness: .6, implement: 1,
  contactSpread: .2, seed: 1396978476,
};

const regions = [
  ["contact", 0, .015], ["presence", .015, .08],
  ["body", .08, .3], ["decay", .3, 1],
];
const edges = [100, 300, 700, 1500, 3000, 6000, 12000, 20000];

function set(macros, key, value) {
  macros[descriptors.get(key).index] = value;
}

function profile(samples, onsetSeconds) {
  const onset = Math.round(onsetSeconds * reference.sampleRate);
  const count = Math.round(1.1 * reference.sampleRate);
  const aligned = samples.slice(onset, onset + count);
  const spectrum = stft(aligned, reference.sampleRate, {
    size: 1024, hop: 128, window: "hann", floorDb: -180,
  });
  const rows = [];
  for (const [region, start, end] of regions) {
    const firstFrame = Math.max(0, Math.ceil(start * reference.sampleRate / spectrum.hop));
    const lastFrame = Math.min(spectrum.frames, Math.ceil(end * reference.sampleRate / spectrum.hop));
    for (let band = 0; band + 1 < edges.length; ++band) {
      const firstBin = Math.ceil(edges[band] * spectrum.size / reference.sampleRate);
      const lastBin = Math.min(
        spectrum.bins, Math.ceil(edges[band + 1] * spectrum.size / reference.sampleRate),
      );
      let power = 0;
      let count = 0;
      for (let frame = firstFrame; frame < lastFrame; ++frame) {
        for (let bin = firstBin; bin < lastBin; ++bin) {
          const db = spectrum.values[frame * spectrum.bins + bin];
          power += 10 ** (db / 10);
          ++count;
        }
      }
      rows.push({
        region, band: `${edges[band]}-${edges[band + 1]}`,
        db: 10 * Math.log10(Math.max(power / Math.max(1, count), 1e-18)),
      });
    }
  }
  return rows;
}

const referenceProfile = profile(reference.samples, .05);

function compare(candidate) {
  const offsets = candidate.map((row, index) => referenceProfile[index].db - row.db);
  const gainDb = offsets.reduce((sum, value) => sum + value, 0) / offsets.length;
  const errors = candidate.map((row, index) => row.db + gainDb - referenceProfile[index].db);
  const regionRmse = Object.fromEntries(regions.map(([name]) => {
    const selected = errors.filter((_, index) => candidate[index].region === name);
    return [name, Math.sqrt(selected.reduce((sum, value) => sum + value * value, 0) / selected.length)];
  }));
  const regionBias = Object.fromEntries(regions.map(([name]) => {
    const selected = errors.filter((_, index) => candidate[index].region === name);
    return [name, selected.reduce((sum, value) => sum + value, 0) / selected.length];
  }));
  const rmse = Math.sqrt(errors.reduce((sum, value) => sum + value * value, 0) / errors.length);
  const signedBandErrorDb = Object.fromEntries(regions.map(([name]) => [
    name,
    errors.filter((_, index) => candidate[index].region === name),
  ]));
  return { rmse, gainDb, regionRmse, regionBias, signedBandErrorDb };
}

function render(values) {
  const macros = [...defaults];
  for (const [key, value] of Object.entries(values)) set(macros, key, value);
  const samples = engine.render({ seconds: 1.1, macros, ...event });
  return compare(profile(samples, 0));
}

const cleanHighAnchors = Object.fromEntries(
  [...descriptors.values()]
    .filter(item => item.key.startsWith("resolved_frequency_") &&
      item.defaultValue >= 7000)
    .map(item => [item.key.replace("frequency", "turbulence"), 0]),
);
const trials = [
  ["coherent anchors", { field_turbulence: 0 }],
  ["low turbulence", { field_turbulence: .35 }],
  ["maximum turbulence", { field_turbulence: 1 }],
  ["no phase diffusion", { field_phase_bandwidth: 0 }],
  ["no neighbour exchange", { field_exchange: 0 }],
  ["no bloom all-pass diffusion", { bloom_diffusion: 0 }],
  ["clean high anchors", cleanHighAnchors],
];
const candidates = trials.map(([name, values]) => ({
  name, values, ...render(values),
}));
candidates.sort((a, b) => a.rmse - b.rmse);
console.log(JSON.stringify({
  reference: referencePath,
  event,
  factory: render({}),
  candidates,
}, null, 2));
