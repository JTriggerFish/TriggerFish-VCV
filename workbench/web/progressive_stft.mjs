import {stft} from "./analysis.mjs";

// Only complete windows are reusable; a new render ID invalidates the cache.
export class ProgressiveStft {
  cache = null;

  analyze(samples, sampleRate, settings, renderId, preview) {
    const size = settings.size ?? 2048, hop = settings.hop ?? size / 4;
    const key = JSON.stringify([renderId, sampleRate, size, hop,
      settings.window ?? "hann", settings.floorDb ?? -140]);
    const frames = preview
      ? Math.max(0, 1 + Math.floor((samples.length - size / 2) / hop))
      : Math.max(1, 1 + Math.ceil((samples.length - 1) / hop));
    const previous = renderId !== undefined && this.cache?.key === key &&
      this.cache.result.frames <= frames ? this.cache.result : null;
    const firstFrame = previous?.frames ?? 0;
    const part = stft(samples, sampleRate, {...settings, firstFrame,
      frameCount: frames - firstFrame});
    const values = new Float32Array(frames * part.bins);
    if (previous) values.set(previous.values);
    values.set(part.values, firstFrame * part.bins);
    const result = {...part, values, frames, incomplete: Boolean(preview),
      peakDb: Math.max(part.peakDb, previous?.peakDb ?? part.floorDb)};
    this.cache = preview && renderId !== undefined ? {key, result} : null;
    return result;
  }
}
