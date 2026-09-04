import createModule from "./triggerfish-percussion.mjs";
import { WasmPercussionEngine } from "./wasm_engine_core.mjs";

export class PercussionEngine extends WasmPercussionEngine {
  static async create(sampleRate = 48000, recipeIndex = 0) {
    const module = await createModule();
    const engine = new PercussionEngine(module, sampleRate, recipeIndex);
    engine.recipes = Array.from(
      { length: module._tf_percussion_recipe_count() }, (_, index) => ({
        index,
        key: module.UTF8ToString(module._tf_percussion_recipe_key(index)),
        name: module.UTF8ToString(module._tf_percussion_recipe_name(index)),
      }),
    );
    return engine;
  }

  render({
    seconds, strength, location, hardness, implement, contactSpread, seed,
    constraint, parameters, macros, routing,
  }) {
    this.setConfiguration(parameters ?? macros, routing);
    this.trigger({
      strength, location, hardness, implement, contactSpread, seed,
      constraint,
    });
    const frames = Math.max(1, Math.round(seconds * this.sampleRate));
    const result = new Float32Array(frames);
    this.processTo(result, 0, frames);
    return result;
  }
}
