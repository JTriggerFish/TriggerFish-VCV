import { pathToFileURL } from "node:url";
import { resolve } from "node:path";
import { performance } from "node:perf_hooks";

import { stft } from "../web/analysis.mjs";

const enginePath = resolve(
  process.argv[2] ?? "build/workbench-wasm/site/engine.mjs",
);
const { PercussionEngine } = await import(pathToFileURL(enginePath));
const sampleRate = 48000;
const seconds = 10;
const event = {
  strength: .8, location: .8, hardness: .65, implement: .75,
  contactSpread: .2, seed: 17,
};

const renders = [];
let crashSamples;
for (const recipeIndex of [0, 1]) {
  const engine = await PercussionEngine.create(sampleRate, recipeIndex);
  const parameters = engine.parameters.map(
    descriptor => descriptor.defaultValue,
  );
  engine.render({ seconds: 1, parameters, ...event });
  const times = [];
  let samples;
  for (let repeat = 0; repeat < 5; ++repeat) {
    const started = performance.now();
    samples = engine.render({
      seconds, parameters, ...event, seed: event.seed + repeat,
    });
    times.push(performance.now() - started);
  }
  if (recipeIndex === 0) crashSamples = samples;
  const renderMedian = median(times);
  renders.push({
    recipe: engine.recipes[recipeIndex].key,
    medianMs: renderMedian,
    minimumMs: Math.min(...times),
    realtimeFraction: renderMedian / (1000 * seconds),
    nanosecondsPerSample: renderMedian * 1e6 / samples.length,
  });
}

function median(values) {
  const sorted = [...values].sort((first, second) => first - second);
  return sorted[Math.floor(sorted.length / 2)];
}

const analyses = [];
for (const size of [1024, 2048, 4096, 8192]) {
  const times = [];
  let result;
  for (let repeat = 0; repeat < 3; ++repeat) {
    const started = performance.now();
    result = stft(crashSamples, sampleRate, {
      size, hop: size / 4, window: "hann", floorDb: -160,
    });
    times.push(performance.now() - started);
  }
  analyses.push({
    size, frames: result.frames, bins: result.bins,
    medianMs: median(times), minimumMs: Math.min(...times),
    outputMiB: result.values.byteLength / (1024 * 1024),
  });
}

console.log(JSON.stringify({
  sampleRate, seconds,
  renders,
  stft: analyses,
}, null, 2));
