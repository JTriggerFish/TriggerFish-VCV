import { stft } from "./analysis.mjs";

self.onmessage = ({ data }) => {
  const { generation, kind, samples, sampleRate, settings, cacheKey } = data;
  try {
    const result = stft(samples, sampleRate, settings);
    self.postMessage(
      { generation, kind, result, cacheKey }, [result.values.buffer],
    );
  } catch (error) {
    self.postMessage({ generation, kind, error: String(error), cacheKey });
  }
};
