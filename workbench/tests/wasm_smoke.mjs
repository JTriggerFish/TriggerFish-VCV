import { pathToFileURL } from "node:url";

const modulePath = process.argv[2];
if (!modulePath) throw new Error("expected the generated module path");

const createModule = (await import(pathToFileURL(modulePath).href)).default;
const wasm = await createModule();
if (wasm._tf_crash_api_version() !== 4) throw new Error("unexpected API version");
if (wasm._tf_crash_macro_count() !== 82) throw new Error("unexpected macro count");
if (wasm.UTF8ToString(wasm._tf_crash_macro_name(0)) !== "Model level") {
  throw new Error("macro metadata is unavailable");
}

const frameCount = 8192;
const outputPointer = wasm._malloc(frameCount * Float32Array.BYTES_PER_ELEMENT);
const handle = wasm._tf_crash_create(48000);
if (!handle || !outputPointer) throw new Error("renderer allocation failed");

function render(seed, location, blockSize) {
  if (!wasm._tf_crash_reset(handle)) throw new Error("reset failed");
  if (!wasm._tf_crash_trigger(handle, 0.8, location, 0.65, 0.75, 0.2, seed)) {
    throw new Error("trigger failed");
  }
  const result = new Float32Array(frameCount);
  for (let first = 0; first < frameCount; first += blockSize) {
    const count = Math.min(blockSize, frameCount - first);
    if (!wasm._tf_crash_process(handle, outputPointer, count)) {
      throw new Error("render failed");
    }
    const start = outputPointer >>> 2;
    result.set(wasm.HEAPF32.subarray(start, start + count), first);
  }
  return result;
}

function equal(first, second) {
  return first.every((sample, index) => Object.is(sample, second[index]));
}

function difference(first, second) {
  let energy = 0;
  for (let index = 0; index < first.length; ++index) {
    const delta = first[index] - second[index];
    energy += delta * delta;
  }
  return energy;
}

const first = render(17, 0.8, 256);
const repeated = render(17, 0.8, 256);
const whole = render(17, 0.8, frameCount);
const variation = render(18, 0.8, 256);

if (!first.every(Number.isFinite)) throw new Error("non-finite output");
if (!equal(first, repeated)) throw new Error("render is not deterministic");
if (!equal(first, whole)) throw new Error("render depends on block size");
if (!(difference(first, variation) > 1e-8)) throw new Error("seed has no effect");

const peak = first.reduce((value, sample) => Math.max(value, Math.abs(sample)), 0);
const energy = first.reduce((value, sample) => value + sample * sample, 0);
const absoluteSum = first.reduce((value, sample) => value + Math.abs(sample), 0);
const earlyEnergy = first.subarray(0, 1024)
  .reduce((value, sample) => value + sample * sample, 0);
if (!(peak > 0 && energy > 0)) throw new Error("render is silent");
const quieterLevel = wasm._tf_crash_macro_default(0) - 6;
if (!wasm._tf_crash_macro_set(handle, 0, quieterLevel) ||
    !wasm._tf_crash_macro_commit(handle)) {
  throw new Error("macro update failed");
}
const quiet = render(17, 0.8, 256);
const quietEnergy = quiet.reduce((value, sample) => value + sample * sample, 0);
if (!(quietEnergy / energy > 0.24 && quietEnergy / energy < 0.26)) {
  throw new Error("model-level macro has the wrong gain law");
}

wasm._tf_crash_destroy(handle);
wasm._free(outputPointer);
console.log(JSON.stringify({
  api: 4, frames: frameCount, peak, energy, absoluteSum, earlyEnergy,
}));
