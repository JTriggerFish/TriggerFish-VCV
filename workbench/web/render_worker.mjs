import { CrashEngine } from "./engine.mjs";

let engine;

self.onmessage = async ({ data }) => {
  const { generation, sampleRate, seconds, macros, event } = data;
  try {
    const started = performance.now();
    if (!engine) engine = await CrashEngine.create(sampleRate);
    else engine.setSampleRate(sampleRate);
    const samples = engine.render({ seconds, macros, ...event });
    self.postMessage({
      generation, samples, sampleRate, elapsedMs: performance.now() - started,
    }, [samples.buffer]);
  } catch (error) {
    self.postMessage({ generation, error: String(error) });
  }
};
