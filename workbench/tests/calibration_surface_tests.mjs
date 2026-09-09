// Test shipped presets against the real compiled descriptors, never a surface
// inferred from the presets themselves (which cannot detect missing controls).
import assert from "node:assert/strict";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import { WasmPercussionEngine } from "../web/wasm_engine_core.mjs";
import { calibrationParameterValues, calibrationPatch } from "../web/instrument_calibrations.mjs";
import { bloomTimingValues, BloomTimingKeys } from "../web/bloom_timing_meta.mjs";
import { recipeAdapter } from "../web/recipe_adapter.mjs";
import { helpFor } from "../web/fit_control_help.mjs";
import { readFile } from "node:fs/promises";
import { validateFit } from "../web/state.mjs";

for (const key of ["field_packet_spread", "field_satellite_density",
  "field_phase_bandwidth", "bloom_rate", "bloom_energy_acceleration",
  "resolved_level_0", "resolved_turbulence_0", "resolved_allocation_0",
  "body_decay_seconds_0"]) {
  assert.ok(helpFor(key).length > 60, `Missing sound-shaping help: ${key}`);
  assert.ok(!/squared norm|spectral-density gradient|eigenmode/.test(helpFor(key)));
}
assert.match(helpFor("field_packet_spread"), /spacing, not how many/);
assert.match(helpFor("field_phase_bandwidth"), /too much can sound like hiss/);

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
    if (recipe === "metal.cymbal.v1") {
      const fit = JSON.parse(await readFile(new URL('../web/crash_calibration.fit.json', import.meta.url)));
      const oldFit = structuredClone(fit);
      const p = oldFit.instrument.nodes.find(n=>n.id==='observation').parameters;
      for(const suffix of ['eq_enabled','low_cut','colour_frequency','colour_gain','high_cut']) {
        const oldSuffix=suffix==='eq_enabled'?'radiation_enabled':suffix;
        p['body_'+oldSuffix]=p['output_'+suffix];
        p['direct_'+oldSuffix]=p['output_'+suffix];
        delete p['output_'+suffix];
      }
      assert.deepEqual(validateFit(oldFit,engine.parameters),fit,
        'old saved fits import into the exact current patch ownership and surface');
      for (const [key, expected] of [["field_distribution",3], ["field_doublet_split",1.25],
        ["field_beat_depth",.3], ["field_beat_rate_tilt",.25]]) {
        const descriptor = engine.parameters.find(d=>d.key===key);
        assert.ok(Math.abs(descriptor.defaultValue-expected)<1e-6, `Gentle reset default: ${key}`);
        // Calibration values are fitted, not obliged to equal reset defaults.
      }
      const original = Object.fromEntries(engine.parameters.map(d=>[d.key,values[d.index]]));
      const {values: changes} = bloomTimingValues(original,.5,engine.parameters);
      const expanded = engine.parameters.map(d=>changes[d.key] ?? original[d.key]);
      const adapter = recipeAdapter(recipe);
      const stored = JSON.parse(JSON.stringify(adapter.withValues(patch,engine.parameters,expanded)));
      const restored = adapter.values(stored,engine.parameters);
      assert.deepEqual(restored,expanded,"saved patch needs no hidden timing state");
      for(const d of engine.parameters)
        if(!BloomTimingKeys.includes(d.key)) assert.equal(restored[d.index],values[d.index]);
      engine.setParameters(restored);
      const edited=restored.slice();
      for(const [key,value] of [['field_distribution',3],['field_doublet_split',1.25],['field_beat_depth',.15],['field_beat_rate_tilt',.25],['resolved_allocation_0',0]])
        edited[engine.parameters.find(d=>d.key===key).index]=value;
      const textured=JSON.parse(JSON.stringify(adapter.withValues(patch,engine.parameters,edited)));
      assert.deepEqual(adapter.values(textured,engine.parameters),edited,'packet controls round-trip without hidden state');
      engine.setParameters(edited);
    }
    console.log(`${id}: ${values.length} compiled controls validated`);
  } finally { engine.destroy(); }
}
