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

const paint = (levels) => Object.fromEntries(
  levels.map((value, index) => [`dense_wash_level_${index}`, value]),
);
const paintFrequencies = Object.fromEntries(
  [120, 500, 1500, 3000, 6000, 10000, 16000, 21000]
    .map((value, index) => [`dense_wash_frequency_${index}`, value]),
);
const decay = (seconds) => Object.fromEntries(
  seconds.map((value, index) => [`body_decay_seconds_${index}`, value]),
);
const decayFrequencies = Object.fromEntries(
  [150, 500, 1500, 6000, 16000]
    .map((value, index) => [`body_decay_frequency_${index}`, value]),
);
const base = {
  body_tone_wash: 1,
  dense_mode_density: 2,
};
const contact = {
  direct_gain: 1.25, direct_colour_frequency: 6500,
  direct_colour_gain: 6, direct_colour_q: .65,
  impact_tone_noise: .92,
};
const spectral = {
  ...paintFrequencies,
  ...paint([-6, -5, -1, 2, 4, 3, -3, -6]),
};
const shapedDecay = {
  ...decayFrequencies, ...decay([3.5, 1.8, 3.2, 5, .7]),
};
const turbulence = {
  turbulence_amount: .2, turbulence_persistence: 1.35,
  turbulence_frequency_0: 300, turbulence_frequency_1: 3500,
  turbulence_frequency_2: 10000,
  turbulence_level_0: -5, turbulence_level_1: 2,
  turbulence_level_2: 3,
};
const refinedSpectral = {
  ...paintFrequencies,
  ...paint([-7, -4.5, -2.5, 0, 2.5, 2, -4, -7]),
};
const refinedDecay = {
  ...decayFrequencies, ...decay([3.5, 2.8, 2.5, 4.5, .35]),
};
const gentleContact = {
  direct_gain: .95, direct_colour_frequency: 7000,
  direct_colour_gain: 4, direct_colour_q: .7,
  impact_tone_noise: .9,
};
const trials = [
  ["wash only", { ...base }],
  ["contact", { ...base, ...contact }],
  ["spectral body", { ...base, ...spectral }],
  ["decay", { ...base, ...shapedDecay }],
  ["turbulence", { ...base, ...turbulence }],
  ["contact + spectrum", { ...base, ...contact, ...spectral }],
  ["spectrum + decay", { ...base, ...spectral, ...shapedDecay }],
  ["contact + spectrum + decay", {
    ...base, ...contact, ...spectral, ...shapedDecay,
  }],
  ["refined spectrum + decay", {
    ...base, ...refinedSpectral, ...refinedDecay,
  }],
  ["refined + gentle contact", {
    ...base, ...gentleContact, ...refinedSpectral, ...refinedDecay,
  }],
  ["spectrum + turbulence", { ...base, ...spectral, ...turbulence }],
  ["combined", {
    ...base, ...contact, ...spectral, ...shapedDecay, ...turbulence,
  }],
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
