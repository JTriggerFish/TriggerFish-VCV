import { stft } from "./analysis.mjs";
import { ProgressiveStft } from "./progressive_stft.mjs";

const synthesis = new ProgressiveStft();

self.onmessage = ({ data }) => {
  const { generation, kind, samples, sampleRate, settings, cacheKey, preview, renderId } = data;
  try {
    const result = kind === "synthesis"
      ? synthesis.analyze(samples, sampleRate, settings, renderId, preview)
      : stft(samples, sampleRate, settings);
    // Keep cached frames in the worker; transfer a separate UI copy.
    const message = preview ? {...result, values: result.values.slice()} : result;
    self.postMessage(
      { generation, kind, result: message, cacheKey }, [message.values.buffer],
    );
  } catch (error) {
    self.postMessage({ generation, kind, error: String(error), cacheKey });
  }
};
