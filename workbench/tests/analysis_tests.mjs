import assert from "node:assert/strict";
import {
  centeredErrorStatistics, fft, stft, windowSamples,
} from "../web/analysis.mjs";
import { waveformEnvelope } from "../web/waveform_view.mjs";
import { calibrateReference, setReferenceGain } from "../web/references.mjs";
import {
  alignedReferenceWindow, normalizationCeilingDb, wheelPanSeconds,
} from "../web/spectrogram.mjs";

const impulseReal = new Float64Array(8);
const rawReference = { sha256: "source", samples: new Float32Array([.001, -.002]) };
const calibratedReference = calibrateReference(rawReference, 40);
assert.equal(calibratedReference.sha256, rawReference.sha256);
assert.ok(Math.abs(calibratedReference.samples[0] - .1) < 1e-7);
assert.equal(calibratedReference.samples[1] / calibratedReference.samples[0], -2);
assert.ok(Math.abs(rawReference.samples[0] - .001) < 1e-9);
assert.equal(calibratedReference.referenceGainDb, 40);
assert.throws(() => calibrateReference(calibratedReference, 40), /already calibrated/);
const quieterReference = setReferenceGain(calibratedReference, 20);
assert.ok(Math.abs(quieterReference.samples[0] - .01) < 1e-8);
assert.deepEqual(setReferenceGain(quieterReference, 40).samples, calibratedReference.samples);
assert.deepEqual(setReferenceGain(quieterReference, 0).samples, rawReference.samples);
assert.throws(() => setReferenceGain(rawReference, NaN), /invalid reference gain/);
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
