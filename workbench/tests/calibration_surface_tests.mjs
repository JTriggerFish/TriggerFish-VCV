// Test shipped presets against the real compiled descriptors, never a surface
// inferred from the presets themselves (which cannot detect missing controls).
import assert from "node:assert/strict";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { WasmPercussionEngine } from "../web/wasm_engine_core.mjs";
import { calibrationParameterValues, calibrationPatch } from "../web/instrument_calibrations.mjs";

const { default: createModule } = await import(pathToFileURL(resolve(process.argv[2])));
const module = await createModule();
const recipes = Array.from({length: module._tf_percussion_recipe_count()}, (_, index) => ({
  index, key: module.UTF8ToString(module._tf_percussion_recipe_key(index)),
}));
for (const [id, recipe, parameter_preset] of [
  ["crash-standard", "metal.cymbal.v1"],
  ["ride-standard", "metal.cymbal.v1"],
  ["gong-standard", "metal.cymbal.v1"],
  ["hihat-standard", "metal.cymbal.v1"],
  ["snare-standard", "drum.snare.v1"],
  ["kick-standard", "drum.kick.v1", "kick"],
]) {
  const engine = new WasmPercussionEngine(module, 48000, recipes.find(r => r.key === recipe).index);
  try {
    const target = {id, recipe, parameter_preset};
    const values = calibrationParameterValues(target, engine.parameters);
    assert.equal(values.length, engine.parameters.length, id);
    const patch = calibrationPatch(target, engine.parameters, values, null);
    assert.equal(patch.recipe, recipe, id);
    engine.setParameters(values);
    console.log(`${id}: ${values.length} compiled controls validated`);
  } finally { engine.destroy(); }
}
