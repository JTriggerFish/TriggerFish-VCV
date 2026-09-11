import { PercussionEngine } from "./engine.mjs";
import { renderOffline } from "./offline_render.mjs";

let engine;
let pending = null;
let running = false;
let revision = 0;
// MessageChannel yields to edits without the repeated-timer minimum delay.
const channel = new MessageChannel();
let resume;
channel.port1.onmessage = () => { const done = resume; resume = null; done(); };
const yieldTask = () => new Promise(resolve => {
  resume = resolve; channel.port2.postMessage(null);
});

self.onmessage = ({ data }) => {
  ++revision;
  pending = data.cancel ? null : data;
  if (!running) drain();
};

async function drain() {
  running = true;
  while (pending) {
    const request = pending, current = revision;
    pending = null;
    await render(request, () => current !== revision);
  }
  running = false;
}

async function render(request, cancelled) {
  const { generation, recipeIndex, sampleRate, seconds } = request;
  try {
    const started = performance.now();
    if (!engine) engine = await PercussionEngine.create(sampleRate, recipeIndex);
    else {
      engine.setRecipe(recipeIndex);
      engine.setSampleRate(sampleRate);
    }
    if (cancelled()) return;
    const samples = await renderOffline(engine, request, {
      cancelled, yieldTask,
      progress: samples => self.postMessage({
        generation, samples, sampleRate, seconds, preview: true,
        elapsedMs: performance.now() - started,
      }, [samples.buffer]),
    });
    if (!samples) return;
    self.postMessage({
      generation, samples, sampleRate, seconds, elapsedMs: performance.now() - started,
    }, [samples.buffer]);
  } catch (error) {
    self.postMessage({ generation, error: String(error) });
  }
}
