// Render any compiled workbench recipe to a mono 32-bit float WAV.
import { writeFile } from "node:fs/promises";
import { pathToFileURL } from "node:url";

const [modulePath, recipeKey, outputPath, ...arguments_] = process.argv.slice(2);
if (!modulePath || !recipeKey || !outputPath) {
  throw new Error(
    "usage: node render_percussion_recipe.mjs MODULE RECIPE OUTPUT.wav " +
    "[seconds=N] [strength=N] [location=N] [hardness=N] [implement=N] " +
    "[contactSpread=N] [seed=N] [parameter=value ...]",
  );
}

const options = new Map(arguments_.map(argument => {
  const split = argument.indexOf("=");
  if (split < 1) throw new Error(`invalid argument: ${argument}`);
  return [argument.slice(0, split), Number(argument.slice(split + 1))];
}));
const option = (key, fallback) => options.has(key) ? options.get(key) : fallback;
const renderOptions = {
  sampleRate: option("sampleRate", 48000),
  seconds: option("seconds", 2),
  strength: option("strength", .8),
  location: option("location", .5),
  hardness: option("hardness", .65),
  implement: option("implement", 1),
  contactSpread: option("contactSpread", .2),
  seed: option("seed", 17),
};
const createModule = (await import(pathToFileURL(modulePath).href)).default;
const wasm = await createModule();

let recipe = -1;
for (let index = 0; index < wasm._tf_percussion_recipe_count(); ++index) {
  if (wasm.UTF8ToString(wasm._tf_percussion_recipe_key(index)) === recipeKey)
    recipe = index;
}
if (recipe < 0) throw new Error(`unknown recipe: ${recipeKey}`);
const sampleRate = renderOptions.sampleRate;
const handle = wasm._tf_percussion_create(recipe, sampleRate);
if (!handle) throw new Error("could not create percussion renderer");

for (let index = 0; index < wasm._tf_percussion_parameter_count(handle); ++index) {
  const key = wasm.UTF8ToString(wasm._tf_percussion_parameter_key(handle, index));
  if (options.has(key) &&
      !wasm._tf_percussion_parameter_set(handle, index, options.get(key))) {
    throw new Error(`could not set ${key}`);
  }
  options.delete(key);
}
for (const key of [
  "seconds", "sampleRate", "strength", "location", "hardness", "implement",
  "contactSpread", "seed",
]) options.delete(key);
if (options.size) throw new Error(`unknown options: ${[...options.keys()].join(", ")}`);
if (!wasm._tf_percussion_commit(handle)) throw new Error("commit failed");
if (!wasm._tf_percussion_trigger(
  handle, renderOptions.strength, renderOptions.location,
  renderOptions.hardness, renderOptions.implement,
  renderOptions.contactSpread, renderOptions.seed,
)) throw new Error("trigger failed");

const frameCount = Math.round(renderOptions.seconds * sampleRate);
const samples = new Float32Array(frameCount);
const blockSize = 512;
const pointer = wasm._malloc(blockSize * Float32Array.BYTES_PER_ELEMENT);
for (let first = 0; first < frameCount; first += blockSize) {
  const count = Math.min(blockSize, frameCount - first);
  if (!wasm._tf_percussion_process(handle, pointer, count))
    throw new Error("render failed");
  const offset = pointer >>> 2;
  samples.set(wasm.HEAPF32.subarray(offset, offset + count), first);
}
wasm._free(pointer);
wasm._tf_percussion_destroy(handle);

const header = Buffer.alloc(44);
header.write("RIFF", 0);
header.writeUInt32LE(36 + samples.byteLength, 4);
header.write("WAVEfmt ", 8);
header.writeUInt32LE(16, 16);
header.writeUInt16LE(3, 20);
header.writeUInt16LE(1, 22);
header.writeUInt32LE(sampleRate, 24);
header.writeUInt32LE(sampleRate * 4, 28);
header.writeUInt16LE(4, 32);
header.writeUInt16LE(32, 34);
header.write("data", 36);
header.writeUInt32LE(samples.byteLength, 40);
await writeFile(outputPath, Buffer.concat([
  header, Buffer.from(samples.buffer, samples.byteOffset, samples.byteLength),
]));
console.log(`${recipeKey}: ${frameCount} frames at ${sampleRate} Hz -> ${outputPath}`);
