import assert from 'node:assert/strict';
import {overlaySpectrum, writeEdgeFraction} from '../web/spectrogram_history.mjs';
const spectrum = (frames, value, incomplete = false) => ({
  frames, bins: 3, size: 4, hop: 1, sampleRate: 4, floorDb: -160,
  values: new Float32Array(frames * 3).fill(value), incomplete,
});
const old = spectrum(8, -20), fresh = spectrum(3, -40, true);
const merged = overlaySpectrum(old, fresh);
assert.equal(merged.frames, 8);
assert.equal(merged.writeFrames, 3);
assert.deepEqual([...merged.values], [...fresh.values, ...old.values.slice(9)]);
assert.ok(old.values.every(x => x === -20), 'Display history must not mutate analysis data');
const restarted = overlaySpectrum(merged, spectrum(1, -60, true));
assert.equal(restarted.values[0], -60);
assert.equal(restarted.values[3], -40, 'Keep the preceding preview after a new edit');
assert.equal(restarted.values[9], -20);
const complete = spectrum(2, -30);
assert.equal(overlaySpectrum(restarted, complete), complete, 'Completed short render removes old tail');
assert.equal(overlaySpectrum(null, fresh), fresh);
assert.equal(writeEdgeFraction(merged, {start: 0, end: 2}), .25);
assert.equal(writeEdgeFraction(merged, {start: 0, end: 2}, .25), .125);
assert.equal(writeEdgeFraction(merged, {start: 1, end: 2}), null);
assert.equal(writeEdgeFraction(complete, {start: 0, end: 2}), null);
const changedGrid = {...spectrum(1, -60, true), hop: 2};
const resampled = overlaySpectrum(old, changedGrid);
assert.equal(resampled.frames, 4);
assert.equal(resampled.values[3], -20);
console.log('Spectrogram history: retained tail, interrupted updates, resolution and write edge pass');
