// Exact DSP, chunked only to service newer edits and publish display previews.
export async function renderOffline(engine, request, {
  cancelled = () => false, progress = () => {},
  yieldTask = () => new Promise(resolve => setTimeout(resolve, 0)),
  now = () => performance.now(),
} = {}) {
  engine.setConfiguration(request.parameters ?? request.macros, request.routing);
  engine.trigger(request.event);
  const frames = Math.max(1, Math.round(request.seconds * request.sampleRate));
  const samples = new Float32Array(frames);
  const block = 2048;
  let lastYield = now(), lastPreview = lastYield, firstPreview = true;
  for (let offset = 0; offset < frames;) {
    if (cancelled()) return null;
    const count = Math.min(block, frames - offset);
    engine.processTo(samples, offset, count);
    offset += count;
    const time = now();
    if (offset < frames && (firstPreview
      ? offset >= (request.previewFrames ?? request.sampleRate * .125)
      : time - lastPreview >= 120)) {
      progress(samples.slice(0, offset));
      firstPreview = false; lastPreview = time;
    }
    if (time - lastYield >= 12) {
      await yieldTask(); lastYield = now();
    }
  }
  return cancelled() ? null : samples;
}
