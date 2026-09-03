import { PercussionEngine } from "./engine.mjs";

let engine;

self.onmessage = async ({ data }) => {
  const { generation, recipeIndex, sampleRate, seconds, parameters, macros,
    routing, event } = data;
  try {
    const started = performance.now();
    if (!engine) engine = await PercussionEngine.create(sampleRate, recipeIndex);
    else {
      engine.setRecipe(recipeIndex);
      engine.setSampleRate(sampleRate);
    }
    const samples = engine.render({
      seconds, parameters: parameters ?? macros, routing, ...event,
    });
    self.postMessage({
      generation, samples, sampleRate, elapsedMs: performance.now() - started,
    }, [samples.buffer]);
  } catch (error) {
    self.postMessage({ generation, error: String(error) });
  }
};
