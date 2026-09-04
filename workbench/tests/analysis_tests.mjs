import assert from "node:assert/strict";
import {
  centeredErrorStatistics, fft, stft, windowSamples,
} from "../web/analysis.mjs";
import {
  matchedModelLevelDb, modelLevelMatch,
} from "../web/level_match.mjs";
import { waveformEnvelope } from "../web/waveform_view.mjs";
import {
  alignedReferenceWindow, normalizationCeilingDb, wheelPanSeconds,
} from "../web/spectrogram.mjs";

const impulseReal = new Float64Array(8);
impulseReal[0] = 1;
const impulseImaginary = new Float64Array(8);
fft(impulseReal, impulseImaginary);
for (const value of impulseReal) assert.ok(Math.abs(value - 1) < 1e-12);
for (const value of impulseImaginary) assert.ok(Math.abs(value) < 1e-12);

assert.deepEqual(centeredErrorStatistics([8, 8, 8]), { mean: 8, rmse: 0 });
assert.deepEqual(centeredErrorStatistics([]), { mean: 0, rmse: 0 });

const hann = windowSamples("hann", 1024);
assert.ok(Math.abs(hann[0]) < 1e-15 && Math.abs(hann.at(-1)) < 1e-15);

const sampleRate = 48000;
const tone = Float32Array.from(
  { length: 4096 }, (_, index) => 0.5 * Math.sin(2 * Math.PI * 3000 * index / sampleRate),
);
const spectrum = stft(tone, sampleRate, { size: 1024, hop: 256, window: "hann" });
const middleFrame = Math.floor(spectrum.frames / 2);
let peakBin = 0;
for (let bin = 1; bin < spectrum.bins; ++bin) {
  if (spectrum.values[middleFrame * spectrum.bins + bin] >
      spectrum.values[middleFrame * spectrum.bins + peakBin]) peakBin = bin;
}
assert.equal(peakBin, 64);
const peakDb = spectrum.values[middleFrame * spectrum.bins + peakBin];
assert.ok(Math.abs(peakDb - 20 * Math.log10(0.5)) < 0.05);
assert.ok(Math.abs(spectrum.peakDb - 20 * Math.log10(0.5)) < 0.05);

const reference = Float32Array.from(tone, value => 0.5 * value);
const synthesis = Float32Array.from(tone, value => 0.25 * value);
const matched = matchedModelLevelDb({
  currentDb: -36, reference, referenceSampleRate: sampleRate,
  synthesis, synthesisSampleRate: sampleRate,
  minimumDb: -60, maximumDb: 0,
});
assert.ok(Math.abs(matched - (-36 + 20 * Math.log10(2))) < 1e-6);

const delayedReference = new Float32Array(tone.length + 480);
delayedReference.set(reference, 480);
const onsetMatched = matchedModelLevelDb({
  currentDb: -36, reference: delayedReference,
  referenceSampleRate: sampleRate, referenceStartSeconds: .01,
  synthesis, synthesisSampleRate: sampleRate,
  minimumDb: -60, maximumDb: 0,
});
assert.ok(Math.abs(onsetMatched - (-36 + 20 * Math.log10(2))) < 1e-6);

const clipped = modelLevelMatch({
  currentDb: -3, reference: tone, referenceSampleRate: sampleRate,
  synthesis: Float32Array.from(tone, value => value * .01),
  synthesisSampleRate: sampleRate, minimumDb: -60, maximumDb: 0,
});
assert.equal(clipped.appliedDb, 0);
assert.equal(clipped.clipped, true);
assert.ok(clipped.requestedDb > 30);

const envelope = waveformEnvelope(
  Float32Array.from([0, -1, .5, 0, .25, -.75, 0, 0]), 4, .5, 2,
);
assert.deepEqual(envelope.y, [-1, .5, .25, -.75]);
assert.deepEqual(envelope.x, [-.25, 0, .5, .75]);

assert.deepEqual(alignedReferenceWindow(2, .125), {
  duration: 1.875, offset: .125,
});
assert.deepEqual(alignedReferenceWindow(2, -1), {
  duration: 2, offset: 0,
});

const forward = wheelPanSeconds(0, 120, 0, 1000, 500, 8);
const backward = wheelPanSeconds(0, -120, 0, 1000, 500, 8);
assert.equal(forward + backward, 0);
assert.equal(normalizationCeilingDb({ peakDb: -18 }, { peakDb: 6 }), -18);
assert.equal(normalizationCeilingDb(null, { peakDb: -9 }), -9);
console.log("workbench analysis tests passed");
