import { pathToFileURL } from "node:url";

const modulePath = process.argv[2];
if (!modulePath) throw new Error("expected the generated module path");

const createModule = (await import(pathToFileURL(modulePath).href)).default;
const wasm = await createModule();
if (wasm._tf_crash_api_version() !== 1 ||
    wasm._tf_percussion_api_version() !== 1) {
  throw new Error("unexpected API version");
}
if (wasm._tf_percussion_recipe_count() !== 4) {
  throw new Error("unexpected recipe count");
}
if (wasm.UTF8ToString(wasm._tf_percussion_recipe_key(1)) !==
    "drum.kick-fm.v1") throw new Error("kick recipe is unavailable");
if (wasm._tf_crash_route_count() !== 5) throw new Error("unexpected route count");
if (wasm._tf_crash_macro_count() !== 167) throw new Error("unexpected macro count");
if (wasm.UTF8ToString(wasm._tf_crash_macro_name(0)) !== "Model level") {
  throw new Error("macro metadata is unavailable");
}

const frameCount = 8192;
const outputPointer = wasm._malloc(frameCount * Float32Array.BYTES_PER_ELEMENT);
const handle = wasm._tf_crash_create(48000);
if (!handle || !outputPointer) throw new Error("renderer allocation failed");
if (wasm._tf_percussion_parameter_count(handle) !== 126) {
  throw new Error("metallic recipe still exposes legacy no-op controls");
}

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

function renderShared(target, seed, blockSize) {
  if (!wasm._tf_percussion_reset(target) ||
      !wasm._tf_percussion_trigger(target, .8, .8, .65, .75, .2, seed)) {
    throw new Error("shared renderer trigger failed");
  }
  const result = new Float32Array(frameCount);
  for (let first = 0; first < frameCount; first += blockSize) {
    const count = Math.min(blockSize, frameCount - first);
    if (!wasm._tf_percussion_process(target, outputPointer, count)) {
      throw new Error("shared renderer failed");
    }
    const start = outputPointer >>> 2;
    result.set(wasm.HEAPF32.subarray(start, start + count), first);
  }
  return result;
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
const preparedSize = wasm._tf_percussion_prepared_size(handle);
const preparedPointer = wasm._malloc(preparedSize);
if (!preparedPointer || !wasm._tf_percussion_export_prepared(
  handle, preparedPointer, preparedSize,
)) throw new Error("prepared recipe export failed");
const restored = wasm._tf_percussion_create_unprepared(0, 48000);
if (!restored || !wasm._tf_percussion_apply_prepared(
  restored, preparedPointer, preparedSize,
)) throw new Error("prepared recipe install failed");
if (!equal(renderShared(handle, 63, 128), renderShared(restored, 63, 128))) {
  throw new Error("prepared Wasm recipe differs from the configured recipe");
}
wasm.HEAPU8[preparedPointer] ^= 0xff;
if (wasm._tf_percussion_apply_prepared(
  restored, preparedPointer, preparedSize,
)) throw new Error("corrupt prepared Wasm recipe was accepted");
wasm._tf_percussion_destroy(restored);
wasm._free(preparedPointer);
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

const kick = wasm._tf_percussion_create(1, 48000);
if (!kick || wasm._tf_percussion_parameter_count(kick) !== 15 ||
    wasm._tf_percussion_route_count(kick) !== 3) {
  throw new Error("kick recipe allocation failed");
}
if (!wasm._tf_percussion_trigger(kick, .8, .5, .6, 1, .2, 91) ||
    !wasm._tf_percussion_process(kick, outputPointer, frameCount)) {
  throw new Error("kick recipe render failed");
}
const kickStart = outputPointer >>> 2;
const kickAudio = wasm.HEAPF32.slice(kickStart, kickStart + frameCount);
if (!kickAudio.every(Number.isFinite) ||
    !kickAudio.some(sample => Math.abs(sample) > 1e-6)) {
  throw new Error("kick recipe output is invalid");
}
wasm._tf_percussion_destroy(kick);

const membrane = wasm._tf_percussion_create(2, 48000);
if (!membrane || wasm._tf_percussion_parameter_count(membrane) !== 34 ||
    wasm._tf_percussion_route_count(membrane) !== 5) {
  throw new Error("membrane recipe allocation failed");
}
if (!wasm._tf_percussion_trigger(membrane, .8, .35, .6, 1, .2, 92) ||
    !wasm._tf_percussion_process(membrane, outputPointer, frameCount)) {
  throw new Error("membrane recipe render failed");
}
const membraneAudio = wasm.HEAPF32.slice(kickStart, kickStart + frameCount);
if (!membraneAudio.every(Number.isFinite) ||
    !membraneAudio.some(sample => Math.abs(sample) > 1e-6)) {
  throw new Error("membrane recipe output is invalid");
}
wasm._tf_percussion_destroy(membrane);

const snare = wasm._tf_percussion_create(3, 48000);
if (!snare || wasm._tf_percussion_parameter_count(snare) !== 52 ||
    wasm._tf_percussion_route_count(snare) !== 7) {
  throw new Error("snare recipe allocation failed");
}
if (!wasm._tf_percussion_trigger(snare, .8, .4, .65, 1, .2, 93) ||
    !wasm._tf_percussion_process(snare, outputPointer, frameCount)) {
  throw new Error("snare recipe render failed");
}
const snareAudio = wasm.HEAPF32.slice(kickStart, kickStart + frameCount);
if (!snareAudio.every(Number.isFinite) ||
    !snareAudio.some(sample => Math.abs(sample) > 1e-6)) {
  throw new Error("snare recipe output is invalid");
}
wasm._tf_percussion_destroy(snare);
wasm._free(outputPointer);
console.log(JSON.stringify({
  api: 1, frames: frameCount, peak, energy, absoluteSum, earlyEnergy,
}));
