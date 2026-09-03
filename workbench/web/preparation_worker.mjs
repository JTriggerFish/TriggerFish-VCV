import { PercussionEngine } from "./engine.mjs";

let engine;

self.onmessage = async ({ data }) => {
  try {
    const started = performance.now();
    if (!engine) {
      engine = await PercussionEngine.create(data.sampleRate, data.recipeIndex);
    } else {
      engine.setRecipe(data.recipeIndex);
      engine.setSampleRate(data.sampleRate);
    }
    engine.stageConfiguration(data.values, data.routing);
    const prepared = engine.exportPrepared();
    self.postMessage({
      generation: data.generation,
      elapsedMs: performance.now() - started,
      prepared: prepared.buffer,
    }, [prepared.buffer]);
  } catch (error) {
    self.postMessage({
      generation: data.generation, error: String(error),
    });
  }
};
